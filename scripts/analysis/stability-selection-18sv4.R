library(tidyverse)
library(magrittr)
library(fs)
library(modelr)
library(spls)
library(pls)
library(collapse)
data_dir <- 'data/processed/'
out_dir <- 'rslt/stability-selection/18sv4/'
dir_create(out_dir)

## DATA INPUTS -----------------------------------------------------------------

# read in and combine 18sv4 data and density estimates
paste(data_dir, 'ncog18sv4.RData', sep = '') |> load()
paste(data_dir, 'density-estimates.RData', sep = '') |> load() 
whales <- inner_join(dens, edna, by = 'cruise')

# dimensions
p <- whales |> dplyr::select(starts_with('asv')) |> ncol()
n.obs <- nrow(whales)

# read in leave one out partitions
loo_partitions <- paste(data_dir, 
                        '_combined-partitions/partitions-18sv4-dens.rds', 
                        sep = '') |> 
  read_rds()

## HYPERPARAMETER SPECIFICATION ------------------------------------------------

# hyperparameter grids
eta_grid_res <- 55
eta_grid <- rev(1 - exp(seq(log(0.075), log(0.6), length = eta_grid_res)))
ncomp_grid <- 4:12

# models to fit
model_grid <- expand_grid(species = c('bm', 'bp', 'mn'),
                          ncomp = ncomp_grid,
                          eta = eta_grid)

# number of models to fit
nrow(model_grid)*nrow(loo_partitions)

# export
write_rds(model_grid, file = paste(out_dir, 'model-grid.rds', sep = ''))

## FIT MODELS TO TRAINING PARTITIONS -------------------------------------------
# 
# # for testing
# j <- 1
# .species <- model_grid$species[1]
# .eta <- model_grid$eta[1]
# .ncomp <- model_grid$ncomp[1]

# function to fit spls model and extract selected asvs
fit_fn <- function(.species, .eta, .ncomp){
  
  out <- lapply(1:nrow(loo_partitions), function(j){
    
    # retrieve training partition
    .train <- loo_partitions$train[j][[1]]
    
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
                  obs.id = loo_partitions$test.id[j],
                  sel.asv = list(sel.asv))
    paste('species = ', .species, 
          ', eta = ', round(.eta, 4), 
          ', ncomp = ', .ncomp, 
          ', partition ', loo_partitions$test.id[j],
          sep = '') |>
      print()
    return(out)
  }) %>%
    Reduce(bind_rows, .)
  
  return(out)
}

# fit every model specified by model grid to every data partition
## LONG RUNTIME ~1hr
paste(out_dir, 'spls-fits/', sep = '_') |> dir_create()
for(i in 1:nrow(model_grid)){
  fit_fn(.species = model_grid$species[i],
         .eta = model_grid$eta[i],
         .ncomp = model_grid$ncomp[i]) |>
    write_rds(file = paste(out_dir, 'spls-fits/', 
                           model_grid$species[i], 
                           model_grid$ncomp[i], 
                           round(model_grid$eta[i], 4), 
                           '.rds', sep = '_'))
}

## CONSTRUCT STABLE SETS -------------------------------------------------------

# retreive selected asvs from procedure above
sel_asvs <- paste(out_dir, '_spls-fits/', sep = '') |>
  dir_ls() |>
  lapply(read_rds) %>%
  Reduce(bind_rows, .)

# compute selection frequencies across partitions for each model configuration
sel_freq <- sel_asvs |>
  unnest(sel.asv) |>
  group_by(species, eta, ncomp, sel.asv) |>
  count() |>
  ungroup()

# compute average number of asvs selected across partitions for each configuration
avg_n_asv <- sel_asvs |>
  mutate(n.asv = map(sel.asv, length)) |>
  unnest(n.asv) |>
  select(species, ncomp, eta, obs.id, n.asv) |>
  fgroup_by(species, ncomp, eta) |>
  fselect(n.asv) |>
  fmean()

# stability threshold (minimax selection prob.)
pi.max <- 0.8

# upper bound on expected no. false positives
EV.max <- 1

# limits for average number of selected asvs
q.max <- sqrt((2*pi.max - 1)*p*EV.max)
q.min <- 15

# generate intervals spanning eta grid (sparsity hyperparameter)
eta_intervals <- expand_grid(eta.min.ix = seq(1, length(eta_grid), by = 2),
                             eta.max.ix = rev(eta.min.ix)) |>
  filter(eta.max.ix - eta.min.ix > 4) |>
  mutate(eta.max = eta_grid[eta.max.ix],
         eta.min = eta_grid[eta.min.ix])

# identify hyperparameter ranges maintaining specified EVmax for all partitions
candidate_ranges <- avg_n_asv |>
  nest(model.rslt = -eta) |>
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
candidate_sets_list <- lapply(1:nrow(candidate_ranges), function(j){
  .parms <- slice(candidate_ranges, j)
  out <- sel_freq |> 
    fsubset(eta >= .parms$eta.min &
              eta <= .parms$eta.max &
              ncomp == .parms$ncomp &
              n >= n.obs*pi.max &
              species == .parms$species) |>
    pull(sel.asv) |>
    unique()
  print(j)
  return(out)
}) 

# identify unique candidate sets at least twice as large as the number of components
candidate_sets <- candidate_sets_list %>%
  mutate(candidate_ranges, ss = .) |>
  group_by(across(-c(eta.min, eta.max))) |>
  slice_min(eta.max - eta.min) |>
  ungroup() |>
  mutate(n.asv = map(ss, length)) |>
  filter(n.asv >= 2*ncomp) |>
  select(species, ncomp, eta.min, eta.max, ss)


## LEAVE ONE OUT PROCEDURE TO PICK PREDICTION-OPTIMAL CANDIDATE ----------------

# generate logratio predictions on held out observations for each stable set
loo_preds <- lapply(1:nrow(candidate_sets), function(i){
  
  # extract stable set
  .ss <- slice(candidate_sets, i)
  
  # fit pls model to every partition using specified stable set
  out <- lapply(1:nrow(loo_partitions), function(j){
    
    # extract data partition
    .partition <- slice(loo_partitions, j)
    
    # extract training data
    x.train <- .partition$train[[1]] |> select(any_of(.ss$ss[[1]]))
    y.train <- .partition$train[[1]] |> pull(.ss$species)
    train <- bind_cols(y = y.train, x.train)
    
    # extract test observation
    x.test <- .partition$test[[1]] |> select(any_of(.ss$ss[[1]]))
    y.test <- .partition$test[[1]] |> pull(.ss$species)
    
    # fit pls model
    .fit <- plsr(y ~ ., data = train, ncomp = .ss$ncomp, 
                 scale = F, center = T)
    
    # compute prediction
    .pred <- predict(.fit, x.test)[, , paste(.ss$ncomp, 'comps')]
    
    # outputs
    out <- .ss |>
      bind_cols(obs.id = .partition$test.id,
                obs.season = .partition$test.season,
                obs.lr = y.test,
                pred.lr = .pred)
    
    c(i, j) |> print()
    return(out)
  }) %>%
    Reduce(bind_rows, .)
}) %>%
  Reduce(bind_rows, .)

# rearrange density data for join with fitted values
dens_long <- dens_raw |>
  select(cruise, season, species, estimate) |>
  rename(obs.dens = estimate)

# seasonal means from training data in leave one out partitions
paste(data_dir, '_partitions/density-partitions.RData', sep = '') |> load()
dens_means_long <- loo_dens |> 
  select(test.id, test.season, seasonal.means) |>
  unnest(seasonal.means) |>
  filter(test.season == season) |>
  pivot_longer(starts_with('log'), 
               names_to = 'species',
               values_to = 'seasonal.mean') |>
  mutate(species = str_remove(species, 'log.') |> str_remove('.mean')) |>
  select(-season) 

# back-transform logratio predictions to original scale
loo_pred_df <- loo_preds |>
  left_join(dens_long, 
            join_by(species, obs.id == cruise, obs.season == season)) |>
  left_join(dens_means_long, 
            join_by(species,
                    obs.season == test.season, 
                    obs.id == test.id)) |>
  mutate(obs.dens.imp = exp(obs.lr + seasonal.mean),
         pred.dens = exp(pred.lr + seasonal.mean)) 

# choose prediction-optimal candidate for each species
stable_sets <- loo_pred_df |>
  group_by(species, ncomp, eta.min, eta.max, ss) |>
  summarize(rmspe.dens = (obs.dens - pred.dens)^2 |> mean() |> sqrt(),
            rmspe.lr = (obs.lr - pred.lr)^2 |> mean() |> sqrt(),
            cor.dens = cor(obs.dens, pred.dens),
            cor.lr = cor(obs.lr, pred.lr)) |>
  group_by(species) |>
  slice_min(rmspe.dens)

# extract predictions from optimal candidates
best_loo_preds <- inner_join(stable_sets, loo_pred_df)

# # plot (x-y)
# best_loo_preds |>
#   select(species, obs.id, obs.dens, pred.dens, obs.lr, pred.lr) |>
#   ggplot(aes(x = obs.dens, y = pred.dens)) +
#   facet_wrap(~species) +
#   geom_point() +
#   geom_abline(intercept = 0, slope = 1) +
#   scale_x_log10() +
#   scale_y_log10()
# 
# # plot (series)
# best_loo_preds |>
#   select(species, obs.id, obs.dens, pred.dens, obs.lr, pred.lr) |>
#   mutate(cruise.ym = ym(obs.id)) |>
#   pivot_longer(c(obs.dens, pred.dens)) |>
#   arrange(species, cruise.ym) |>
#   ggplot(aes(x = cruise.ym, y = value, linetype = name)) +
#   facet_wrap(~species, ncol = 1) +
#   geom_path()

## BOOTSTRAP PREDICTION INTERVALS ----------------------------------------------

# function to fit models to .reps boostrap samples of training partition
boot_fn <- function(.species, .ncomp, .asv, .reps){
  
  rslt <- lapply(1:nrow(loo_partitions), function(j){
    
    # select data partition
    .partition <- slice(loo_partitions, j)
    
    # extract training and test data
    x.train <- dplyr::select(.partition$train[[1]], any_of(.asv))
    y.train <- pull(.partition$train[[1]], {{.species}})
    x.test <- dplyr::select(.partition$test[[1]], any_of(.asv))
    y.test <- pull(.partition$test[[1]], {{.species}})
    train <- bind_cols(y = y.train, x.train)
    test <- bind_cols(y = y.test, x.test)
    
    # generate predictions from bootstrap samples
    out <- bootstrap(train, n = .reps) |>
      mutate(pred.lr = map(strap, function(.data){
        fit <- plsr(y ~ ., ncomp = .ncomp, 
                    data = as.data.frame(.data), 
                    scale = F, center = T)
        pred <- predict(fit, test)[, , paste(.ncomp, 'comps')]
      })) |>
      unnest(pred.lr) |>
      select(.id, pred.lr) |>
      mutate(test.id = .partition$test.id,
             species = .species)
    
    print(j)
    return(out)
  }) 
  
  out <- Reduce(bind_rows, rslt)
  
  return(out)
}

# fit stable set models to 1k bootstrap samples for each partition
## LONG RUNTIME: ~5min
paste(out_dir, '_bootstrap-preds/', sep = '') |> dir_create()
for(j in 1:3){
  boot_fn(.species = stable_sets$species[j],
          .ncomp = stable_sets$ncomp[j],
          .asv = stable_sets$ss[j][[1]],
          .reps = 1000) |>
    write_rds(paste(out_dir, 
                    '_bootstrap-preds/model-', 
                    stable_sets$species[j], 
                    '.rds', 
                    sep = ''))
}

# read in results
boot_preds <- paste(out_dir, '_bootstrap-preds/', sep = '') |> 
  dir_ls() |>
  lapply(read_rds) %>%
  Reduce(bind_rows, .)

# compute quantiles of predictions
pred_quantiles <- loo_dens |> 
  select(test.id, test.season, seasonal.means) |>
  unnest(seasonal.means) |>
  filter(test.season == season) |>
  pivot_longer(starts_with('log'), 
               names_to = 'test.species',
               values_to = 'train.seasonal.mean') |>
  mutate(test.species = str_remove(test.species, 'log.') |> str_remove('.mean')) |>
  select(-season) |>
  right_join(boot_preds, 
             join_by(test.id, 
                     test.species == species)) |>
  mutate(pred.ss = exp(pred.lr + train.seasonal.mean)) |>
  group_by(test.id, test.season, test.species) |>
  summarize(across(starts_with('pred.'),
                   .fns = list(qlo = ~quantile(.x, 0.05),
                               qhi = ~quantile(.x, 0.95)),
                   .names = '{.col}.{.fn}'))

# join with leave one out predictions
pred_df <- loo_pred_df |>
  inner_join(stable_sets) |>
  left_join(pred_quantiles, 
            join_by(obs.id == test.id, 
                    obs.season == test.season, 
                    species == test.species)) |>
  select(species, starts_with('obs'), starts_with('pred'), -ends_with('.imp')) 

## EXPORT STABLE SETS AND LEAVE ONE OUT PREDICTIONS ----------------------------

write_rds(pred_df, file = paste(out_dir, 'loo-preds.rds', sep = ''))
write_rds(stable_sets, file = paste(out_dir, 'stable-sets.rds', sep = ''))
