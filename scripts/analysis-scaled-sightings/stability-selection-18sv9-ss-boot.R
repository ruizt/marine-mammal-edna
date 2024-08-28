library(tidyverse)
library(fs)
library(modelr)
library(spls)
library(pls)
library(collapse)
data_dir <- 'data/processed/'
out_dir <- 'rslt/stability-selection/18sv9-ss/'
dir_create(out_dir)

## DATA INPUTS -----------------------------------------------------------------

# read in and combine 18sv9 data and scaled sightings
paste(data_dir, 'ncog18sv9.RData', sep = '') |> load()
paste(data_dir, 'mm-sightings.RData', sep = '') |> load() 
whales <- inner_join(sightings, edna, by = 'cruise')

# dimensions
p <- whales |> dplyr::select(starts_with('asv')) |> ncol()
n.obs <- nrow(whales)

## HYPERPARAMETER SPECIFICATION ------------------------------------------------

# hyperparameter grids
eta_grid_res <- 50
eta_grid <- rev(1 - exp(seq(log(0.075), log(0.6), length = eta_grid_res)))
ncomp_grid <- 4:12

## SELECTION ON BOOTSTRAP SAMPLES ----------------------------------------------

# generate boostrap samples
nboot <- 100
set.seed(82624)
bsamples <- bootstrap(whales, nboot)

# models to fit
model_grid <- expand_grid(species = c('bm', 'bp', 'mn'),
                          ncomp = ncomp_grid,
                          eta = eta_grid)

# function to fit spls model and extract selected asvs
fit_fn <- function(.species, .eta, .ncomp){
  
  out <- lapply(1:nboot, function(j){
    
    # retreive training bootstrap sample
    .train <- bsamples$strap[j][[1]] |> as_tibble()
    
    # separate predictors and response
    x.train <- .train |> as.data.frame() |> select(starts_with('asv'))
    y.train <- .train |> as.data.frame() |> pull(.species)
    
    # fit spls model with specified parameters
    fit <- spls(x.train, y.train, 
                K = .ncomp, eta = .eta, 
                scale.x = F, scale.y = F)
    
    # selected variables
    sel.asv <- fit$projection |> rownames()
    
    # outputs
    out <- tibble(species = .species,
                  eta = .eta,
                  ncomp = .ncomp,
                  .id = bsamples$.id[j],
                  sel.asv = list(sel.asv))
    paste('species = ', .species, 
          ', eta = ', round(.eta, 4), 
          ', ncomp = ', .ncomp, 
          ', boostrap sample ', bsamples$.id[j],
          sep = '') |>
      print()
    return(out)
  }) %>%
    Reduce(bind_rows, .)
  
  return(out)
}

paste(out_dir, 'rslt-boot/', sep = '_') |> dir_create()
for(i in 1:nrow(model_grid)){
  fit_fn(.species = model_grid$species[i],
         .eta = model_grid$eta[i],
         .ncomp = model_grid$ncomp[i]) |>
    write_rds(file = paste(out_dir, 'rslt-boot/', 
                           model_grid$species[i], 
                           model_grid$ncomp[i], 
                           round(model_grid$eta[i], 4), 
                           '.rds', sep = '_'))
}

## CONSTRUCT CANDIDATE STABLE SETS ---------------------------------------------

# retreive selected asvs from procedure above
sel_asvs <- paste(out_dir, '_rslt-boot/', sep = '') |>
  dir_ls() |>
  lapply(read_rds) %>%
  Reduce(bind_rows, .)

# compute selection frequencies from bootstrap samples for each model configuration
sel_freq <- sel_asvs |>
  unnest(sel.asv) |>
  group_by(species, eta, ncomp, sel.asv) |>
  count() |>
  ungroup()

# stability threshold (minimax selection prob.)
pi.max <- 0.7

# upper bound on expected no. false positives
EV.max <- 4

# limits for average number of selected asvs
q.max <- sqrt((2*pi.max - 1)*p*EV.max)
q.min <- 10

# generate intervals spanning eta grid (sparsity hyperparameter)
eta_intervals <- expand_grid(eta.min.ix = seq(1, length(eta_grid), by = 2),
            eta.max.ix = rev(eta.min.ix)) |>
  filter(eta.max.ix - eta.min.ix > 4) |>
  mutate(eta.max = eta_grid[eta.max.ix],
         eta.min = eta_grid[eta.min.ix])

# identify hyperparameter ranges maintaining specified EVmax for each ncomp
candidate_ranges <- sel_asvs |>
  mutate(n.asv = map(sel.asv, length)) |>
  unnest(n.asv) |>
  select(species, ncomp, eta, .id, n.asv) |>
  fgroup_by(species, ncomp, eta) |>
  fselect(n.asv) |>
  fmean() |>
  nest(model.rslt = -c(eta, ncomp)) |>
  mutate(interval = list(eta_intervals)) |>
  unnest(cols = interval) |>
  filter(eta >= eta.min, eta <= eta.max) |>
  unnest(cols = model.rslt) |>
  fgroup_by(species, ncomp, eta.min, eta.max) |>
  fsummarize(q = fmean(n.asv),
             sd = fsd(n.asv)) |>
  filter(q < q.max, q > q.min, sd/q < 1)  |>
  select(species, ncomp, eta.min, eta.max)

# compute stable sets from each candidate parameter configuration
# (then identify unique stable sets and match to *narrowest* eta range)
candidate_sets <- lapply(1:nrow(candidate_ranges), function(j){
  .parms <- slice(candidate_ranges, j)
  out <- sel_freq |> 
    fsubset(eta >= .parms$eta.min &
              eta <= .parms$eta.max &
              ncomp == .parms$ncomp &
              n >= nboot*pi.max &
              species == .parms$species) |>
    pull(sel.asv) |>
    unique()
  print(j)
  return(out)
}) %>%
  mutate(candidate_ranges, ss = .) |>
  group_by(species, ncomp, ss) |>
  slice_min(eta.max - eta.min) |>
  ungroup() |>
  filter(map(ss, length) > ncomp)

# inspect
candidate_sets$ncomp

## LEAVE ONE OUT PROCEDURE TO PICK PREDICTION-OPTIMAL CANDIDATE ----------------

# read in leave one out partitions
loo_partitions <- paste(data_dir, 
                        '_combined-partitions/partitions-18sv9-ss.rds', 
                        sep = '') |> 
  read_rds()

# generate logratio predictions on held out observations for each candidate set
## LONG RUNTIME: ~5min
loo_preds <- lapply(1:nrow(candidate_sets), function(i){
  
  # extract candidate stable set
  .candidate <- slice(candidate_sets, i)
  
  # fit pls model to every partition using specified candidate
  out <- lapply(1:nrow(loo_partitions), function(j){
    
    # extract data partition
    .partition <- slice(loo_partitions, j)
    
    # extract training data
    x.train <- .partition$train[[1]] |> select(any_of(.candidate$ss[[1]]))
    y.train <- .partition$train[[1]] |> pull(.candidate$species)
    train <- bind_cols(y = y.train, x.train)
    
    # extract test observation
    x.test <- .partition$test[[1]] |> select(any_of(.candidate$ss[[1]]))
    y.test <- .partition$test[[1]] |> pull(.candidate$species)
    
    # fit pls model
    .fit <- plsr(y ~ ., data = train, ncomp = .candidate$ncomp, 
                 scale = F, center = T)
    
    # compute prediction
    .pred <- predict(.fit, x.test)[, , paste(.candidate$ncomp, 'comps')]
    
    # outputs
    out <- .candidate |>
      select(species, ncomp, eta.min, eta.max) |>
      bind_cols(obs.id = .partition$test.cruise,
                obs.season = .partition$test.season,
                obs.lr = y.test,
                pred.lr = .pred)
    
    c(i, j) |> print()
    return(out)
  }) %>%
    Reduce(bind_rows, .)
}) %>%
  Reduce(bind_rows, .)

# rearrange raw scaled sighting data for join with fitted values
sightings_raw_long <- sightings_raw |>
  pivot_longer(c(bp, bm, mn), names_to = 'species', values_to = 'obs.ss')

# seasonal means from training data in leave one out partitions
paste(data_dir, '_cv/mm-sightings-partitions.RData', sep = '') |> load()
ss_means_long <- loo_sightings |> 
  select(test.cruise, test.season, seasonal.means) |>
  unnest(seasonal.means) |>
  filter(test.season == season) |>
  pivot_longer(starts_with('log'), 
               names_to = 'species',
               values_to = 'seasonal.mean') |>
  mutate(species = str_remove(species, 'log.') |> str_remove('.mean')) |>
  select(-season) 

# back-transform logratio predictions to original scale
loo_pred_df <- loo_preds |>
  left_join(sightings_raw_long, 
            join_by(species, obs.id == cruise, obs.season == season)) |>
  left_join(ss_means_long, 
            join_by(species,
                    obs.season == test.season, 
                    obs.id == test.cruise)) |>
  mutate(obs.ss.imp = exp(obs.lr + seasonal.mean),
         pred.ss = exp(pred.lr + seasonal.mean)) 

# choose prediction-optimal candidate set for each whale species
best_settings <- loo_pred_df |>
  group_by(species, ncomp, eta.min, eta.max) |>
  summarize(rmspe.lr = sqrt(mean((obs.lr - pred.lr)^2)),
            rmspe.ss = sqrt(mean((obs.ss - pred.ss)^2)),
            cor.lr = cor(obs.lr, pred.lr),
            cor.ss = cor(obs.ss, pred.ss)) |>
  ungroup() |>
  group_by(species) |>
  slice_min(rmspe.ss)

# inspect
best_settings

# plot (x-y)
loo_pred_df |>
  inner_join(best_settings, join_by(species, ncomp, eta.min, eta.max)) |>
  select(species, obs.id, obs.ss, pred.ss) |>
  ggplot(aes(x = obs.ss, y = pred.ss)) +
  facet_wrap(~species) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  scale_x_log10() +
  scale_y_log10()

# plot (series)
loo_pred_df |>
  inner_join(best_settings, join_by(species, ncomp, eta.min, eta.max)) |>
  select(species, obs.id, obs.ss, pred.ss) |>
  mutate(cruise.ym = ym(obs.id)) |>
  pivot_longer(c(obs.ss, pred.ss)) |>
  arrange(species, cruise.ym) |>
  ggplot(aes(x = cruise.ym, y = value, linetype = name)) +
  facet_wrap(~species, ncol = 1) +
  geom_path()
   
# recover stable sets
best_ss <- candidate_sets |>
  inner_join(best_settings, join_by(species, ncomp, eta.min, eta.max)) 

# inspect
best_ss |> select(ss)

## ASSESS STABLE SET CONSISTENCY -----------------------------------------------

# generate sets of bootstrap samples from each leave-one-out training partition
nboot_val <- 50
set.seed(82724)
val_bsamples <- loo_partitions |>
  select(test.cruise, train) |>
  mutate(bsamples = map(train, ~bootstrap(.x, n = nboot_val))) |>
  select(test.cruise, bsamples) |>
  unnest(bsamples)

# recover hyperparameter settings used to construct optimal stable sets
best_ss_model_grid <- best_ss |>
  select(species, ncomp, eta.min, eta.max) |>
  expand_grid(eta = eta_grid) |>
  filter(eta <= eta.max, eta >= eta.min) |>
  select(species, ncomp, eta)

# function to fit spls model and extract selected asvs
val_fit_fn <- function(.species, .eta, .ncomp){
  
  out <- lapply(1:nrow(val_bsamples), function(j){
    
    # retreive training bootstrap sample
    .train <- val_bsamples$strap[j][[1]] |> as_tibble()
    
    # separate predictors and response
    x.train <- .train |> as.data.frame() |> select(starts_with('asv'))
    y.train <- .train |> as.data.frame() |> pull(.species)
    
    # fit spls model with specified parameters
    fit <- spls(x.train, y.train, 
                K = .ncomp, eta = .eta, 
                scale.x = F, scale.y = F)
    
    # selected variables
    sel.asv <- fit$projection |> rownames()
    
    # outputs
    out <- tibble(species = .species,
                  eta = .eta,
                  ncomp = .ncomp,
                  boot.id = val_bsamples$.id[j],
                  obs.id = val_bsamples$test.cruise[j],
                  sel.asv = list(sel.asv))
    paste('species = ', .species, 
          ', eta = ', round(.eta, 4), 
          ', ncomp = ', .ncomp, 
          ', boostrap sample ', j,
          sep = '') |>
      print()
    return(out)
  }) %>%
    Reduce(bind_rows, .)
  
  return(out)
}

# fit spls models to each bootstrap sample
paste(out_dir, 'validation-boot/', sep = '_') |> dir_create()
for(i in 1:nrow(best_ss_model_grid)){
  val_fit_fn(.species = best_ss_model_grid$species[i],
         .eta = best_ss_model_grid$eta[i],
         .ncomp = best_ss_model_grid$ncomp[i]) |>
    write_rds(file = paste(out_dir, 'validation-boot/', 
                           best_ss_model_grid$species[i], 
                           best_ss_model_grid$ncomp[i], 
                           round(best_ss_model_grid$eta[i], 4), 
                           '.rds', sep = '_'))
}

# read in results
val_rslt <- paste(out_dir, 'validation-boot/', sep = '_') |> 
  dir_ls() |>
  lapply(read_rds) %>%
  Reduce(bind_rows, .)

# function to compute "soft intersection" of sets in <set_list>
intersect_fn <- function(set_list, thresh){
  out <- tibble(asv = Reduce(c, set_list)) |> 
    group_by(asv) |> 
    count() |>
    filter(n >= thresh*length(set_list)) |>
    pull(asv)
  return(out)
}

# function to compute no. of elements in union of sets in <set_list>
union_fn <- function(set_list){
  out <- Reduce(c, set_list) |> unique()
  return(out)
}

# compute stable sets
val_rslt |>
  unnest(sel.asv) |>
  group_by(species, obs.id, eta, sel.asv) |> 
  count() |>
  ungroup() |>
  fsubset(n > pi.max*nboot_val) |>
  select(species, obs.id, sel.asv) |>
  group_by(species, obs.id) |>
  nest(sel.asv = sel.asv) |>
  mutate(sel.asv = map(sel.asv, ~pull(.x, sel.asv))) |>
  group_by(species) |>
  summarize(int = intersect_fn(sel.asv, thresh = 0.7) |> list(),
            un = union_fn(sel.asv) |> list()) |>
  mutate(j.index = map2(int, un, ~length(.x)/length(.y))) |>
  unnest(j.index)


