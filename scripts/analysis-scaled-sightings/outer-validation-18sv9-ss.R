library(tidyverse)
library(fs)
library(modelr)
library(spls)
library(pls)
library(collapse)
library(magrittr)
data_dir <- 'data/processed/'
out_dir <- 'rslt/outer-validation/18sv9-ss/'
dir_create(out_dir)
paste(out_dir, 'nested-loo-selection/', sep = '_') |> dir_create()

## DATA INPUTS -----------------------------------------------------------------

# read in and combine 18sv9 data and scaled sightings
paste(data_dir, 'ncog18sv9.RData', sep = '') |> load()
paste(data_dir, 'mm-sightings.RData', sep = '') |> load() 
whales <- inner_join(sightings, edna, by = 'cruise')

# dimensions
p <- whales |> dplyr::select(starts_with('asv')) |> ncol()
n.obs <- nrow(whales)

# read in nested leave one out partitions
partitions_nested <- paste(data_dir, 
                           '_combined-partitions/nested-partitions-18sv9-ss.rds', 
                           sep = '') |> 
  read_rds()

## HYPERPARAMETER SPECIFICATION ------------------------------------------------

# hyperparameter grids
eta_grid_res <- 25
eta_grid <- rev(1 - exp(seq(log(0.075), log(0.55), length = eta_grid_res)))
ncomp_grid <- seq(4, 10, by = 1)

## SELECTION ON NESTED LEAVE ONE OUT PARTITIONS --------------------------------

# models to fit
model_grid <- expand_grid(species = c('bm', 'bp', 'mn'),
                          ncomp = ncomp_grid,
                          eta = eta_grid)

# check no. models to fit
nrow(partitions_nested)*nrow(model_grid)

# function to fit spls model and extract selected asvs
fit_fn <- function(.species, .eta, .ncomp){
  
  out <- lapply(1:nrow(partitions_nested), function(j){
    
    # retrieve training data
    .train <- partitions_nested$train[j][[1]]
    
    # separate predictors and response
    x.train <- .train |> as.data.frame() |> select(starts_with('asv'))
    y.train <- .train |> as.data.frame() |> pull(.species)
    
    # fit spls model with specified parameters
    .fit <- spls(x.train, y.train, 
                K = .ncomp, eta = .eta, 
                scale.x = F, scale.y = F)
    
    # selected variables
    .asv <- .fit$projection |> rownames()
    
    # outputs
    out <- tibble(species = .species,
                  eta = .eta,
                  ncomp = .ncomp,
                  outer.id = partitions_nested$outer.id[j],
                  inner.id = partitions_nested$inner.id[j],
                  asv = list(.asv))
    paste('species = ', .species, 
          ', eta = ', round(.eta, 4), 
          ', ncomp = ', .ncomp, 
          ', partition ', partitions_nested$outer.id[j],
          '-', partitions_nested$inner.id[j],
          sep = '') |>
      print()
    return(out)
  }) %>%
    Reduce(bind_rows, .)
  
  return(out)
}

# fit each model on every nested partition
## LONG RUNTIME: 8-10hr
for(i in 1:nrow(model_grid)){
  fit_fn(.species = model_grid$species[i],
         .eta = model_grid$eta[i],
         .ncomp = model_grid$ncomp[i]) |>
    write_rds(file = paste(out_dir, 'nested-loo-selection/', 
                           model_grid$species[i], 
                           model_grid$ncomp[i], 
                           round(model_grid$eta[i], 4), 
                           '.rds', sep = '_'))
}

## CONSTRUCT CANDIDATE STABLE SETS FOR EACH OUTER PARTITION --------------------

# retrieve selected asvs from procedure above
sel_asvs <- paste(out_dir, '_nested-loo-selection/', sep = '') |>
  dir_ls() |>
  lapply(read_rds) %>%
  Reduce(bind_rows, .)

# compute selection frequencies from bootstrap samples for each model configuration
sel_freq <- sel_asvs |>
  unnest(asv) |>
  group_by(outer.id, species, eta, ncomp, asv) |>
  count() |>
  ungroup()

# stability threshold (minimax selection prob.)
pi.max <- 0.8

# upper bound on expected no. false positives
EV.max <- 0.5

# limits for average number of selected asvs
q.max <- sqrt((2*pi.max - 1)*p*EV.max)
q.min <- 15

# generate intervals spanning eta grid (sparsity hyperparameter)
eta_intervals <- expand_grid(eta.min.ix = seq(1, length(eta_grid), by = 2),
                             eta.max.ix = rev(eta.min.ix)) |>
  filter(eta.max.ix - eta.min.ix > 2) |>
  mutate(eta.max = eta_grid[eta.max.ix],
         eta.min = eta_grid[eta.min.ix])

# identify hyperparameter ranges maintaining specified EVmax for each ncomp
candidate_ranges <- sel_asvs |>
  mutate(n.asv = map(asv, length)) |>
  unnest(n.asv) |>
  select(species, ncomp, eta, outer.id, n.asv) |>
  fgroup_by(outer.id, species, ncomp, eta) |>
  fselect(n.asv) |>
  fmean() |>
  nest(model.rslt = -eta) |>
  mutate(interval = list(eta_intervals)) |>
  unnest(cols = interval) |>
  filter(eta >= eta.min, eta <= eta.max) |>
  unnest(cols = model.rslt) |>
  fgroup_by(outer.id, species, ncomp, eta.min, eta.max) |>
  fsummarize(q = fmean(n.asv),
             sd = fsd(n.asv)) |>
  filter(q < q.max, q > q.min, sd/q < 1)  |>
  select(outer.id, species, ncomp, eta.min, eta.max, q)

# # filter to ranges that are common across all outer partitions
# candidate_ranges <- candidate_ranges_all %>%
#   group_by(species, ncomp, eta.min, eta.max) |>
#   count() |>
#   filter(n == 25) |>
#   expand_grid(outer.id = partitions_nested$outer.id |> unique()) |>
#   select(outer.id, species, ncomp, eta.min, eta.max)

# compute stable sets from each candidate parameter configuration
## LONG RUNTIME ~10min
candidate_sets_list <- lapply(1:nrow(candidate_ranges), function(j){
  .parms <- slice(candidate_ranges, j)
  out <- sel_freq |> 
    fsubset(outer.id == .parms$outer.id &
              eta >= .parms$eta.min &
              eta <= .parms$eta.max &
              ncomp == .parms$ncomp &
              n >= 24*pi.max &
              species == .parms$species) |>
    pull(asv) |>
    unique()
  print(j)
  return(out)
}) 

# identify unique sets for each outer partition and match to narrowest eta range
candidate_sets <- candidate_sets_list %>%
  mutate(candidate_ranges, ss = .) |>
  group_by(across(-c(eta.min, eta.max))) |>
  slice_min(eta.max - eta.min) |>
  ungroup() |>
  mutate(ss.size = map(ss, length)) |>
  filter(ss.size > q.min) |>
  select(-ss.size)

# inspect
candidate_sets

# export result
write_rds(candidate_sets, 
          file = paste(out_dir, 'candidate-sets.rds', sep = ''))

## LEAVE ONE OUT PROCEDURE TO PICK PREDICTION-OPTIMAL CANDIDATE ----------------

candidate_sets <- paste(out_dir, 'candidate-sets.rds', sep = '') |> read_rds()

# i <- j <- 1
# logratio predictions on inner held out observations for each candidate
## LONG RUNTIME: ~10min
nested_loo_preds <- lapply(1:nrow(candidate_sets), function(i){
  
  # extract candidate stable set
  .candidate <- slice(candidate_sets, i)
  
  # extract corresponding collection of data partitions
  .partitions <- filter(partitions_nested, outer.id == .candidate$outer.id)
  
  out <- lapply(1:nrow(.partitions), function(j){
  
    # extract partition
    .partition <- slice(.partitions, j)
      
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
    .pred <- predict(.fit, x.test, type = 'response', ncomp = .candidate$ncomp)[, , 1]
    
    # outputs
    out <- .candidate |>
      select(outer.id, species, ncomp, eta.min, eta.max) |>
      bind_cols(inner.id = .partition$inner.id, 
                obs.lr = y.test,
                pred.lr = .pred)
    
    c(i, j) |> print()
    return(out)
  }) %>%
    Reduce(bind_rows, .)

  return(out)
}) %>%
  Reduce(bind_rows, .)
  
# rearrange raw scaled sighting data for join with fitted values
sightings_raw_long <- sightings_raw |>
  pivot_longer(c(bp, bm, mn), names_to = 'species', values_to = 'obs.ss')

# seasonal means from training data in leave one out partitions
ss_means_long <- partitions_nested |>
  select(outer.id, inner.id, seasonal.means) |>
  left_join(select(sightings_raw, cruise, season), 
            join_by(inner.id == cruise)) |>
  rename(inner.season = season) |>
  unnest(seasonal.means) |>
  filter(inner.season == season) |>
  pivot_longer(starts_with('log'), 
               names_to = 'species',
               values_to = 'seasonal.mean') |>
  mutate(species = str_remove(species, 'log.') |> str_remove('.mean')) |>
  select(-season) 
  
# back-transform logratio predictions to original scale
loo_pred_df <- nested_loo_preds |>
  left_join(sightings_raw_long, 
            join_by(species, inner.id == cruise)) |>
  left_join(ss_means_long, 
            join_by(outer.id, inner.id, species, season == inner.season)) |>
  mutate(obs.ss.imp = exp(obs.lr + seasonal.mean),
         pred.ss = exp(pred.lr + seasonal.mean)) 

# choose prediction-optimal candidate set for each whale species
best_settings <- loo_pred_df |>
  group_by(outer.id, species, ncomp, eta.min, eta.max) |>
  summarize(inner.rmspe = sqrt(mean((obs.ss - pred.ss)^2)),
            inner.cor = cor(obs.ss, pred.ss)) |>
  ungroup() |>
  group_by(outer.id, species) |>
  slice_min(inner.rmspe) |>
  ungroup() |>
  arrange(species, outer.id) |>
  print(n = 100)

# retrieve stable sets
best_sets <- candidate_sets |>
  inner_join(best_settings)

# export
write_rds(best_sets, file = paste(out_dir, 'best-sets.rds', sep = ''))

## ASSESS STABLE SET CONSISTENCY -----------------------------------------------

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

# general consistency across all candidate sets considered
best_sets |> 
  select(species, ss) |>
  group_by(species) |>
  summarize(int = intersect_fn(ss, 0.5) |> list(),
            un = union_fn(ss) |> list()) |>
  mutate(j.index = map2(int, un, ~length(.x)/length(.y))) |>
  unnest(j.index) 

## "OUTER" PREDICTIONS ---------------------------------------------------------

# read in nonnested leave one out partitions
partitions <- paste(data_dir, 
                           '_combined-partitions/partitions-18sv9-ss.rds', 
                           sep = '') |> 
  read_rds()

# fit models on outer training partitions and compute predictions
preds_outer <- lapply(1:nrow(best_sets), function(i){

  # select stable set
  .ss <- slice(best_sets, i)
  
  # extract corresponding outer partition
  .partition <- partitions |> filter(test.id == .ss$outer.id)
  
  # extract training data and column filter to stable set asvs
  x.train <- .partition$train[[1]] |> select(any_of(.ss$ss[[1]]))
  y.train <- .partition$train[[1]] |> pull(.ss$species)
  .train <- bind_cols(y = y.train, x.train)
  
  # extract test data and column filter to stable set asvs
  x.test <- .partition$test[[1]] |> select(any_of(.ss$ss[[1]]))
  y.test <- .partition$test[[1]] |> pull(.ss$species)
  
  # fit pls model on training data
  .fit <- plsr(y ~ ., data = .train, ncomp = .ss$ncomp, scale = F, center = T)
  
  # compute prediction on outer test partition
  .pred <- predict(.fit, x.test, ncomp = .ss$ncomp)[, , 1]
  
  # outputs
  out <- .ss |>
    bind_cols(obs.lr = y.test,
              pred.lr = .pred)
  
  print(i)
  return(out)
}) %>%
  Reduce(bind_rows, .)

# seasonal means from outer partitions
ss_means_long_outer <- partitions |>
  select(test.id, seasonal.means) |>
  left_join(select(sightings_raw, cruise, season), 
            join_by(test.id == cruise)) |>
  rename(outer.season = season) |>
  unnest(seasonal.means) |>
  filter(outer.season == season) |>
  pivot_longer(starts_with('log'), 
               names_to = 'species',
               values_to = 'seasonal.mean') |>
  mutate(species = str_remove(species, 'log.') |> str_remove('.mean')) |>
  select(-season) 

# back-transform outer predictions to original scale
preds_outer_df <- preds_outer |>
  left_join(ss_means_long_outer, 
            join_by(species, outer.id == test.id)) |>
  left_join(sightings_raw_long,
            join_by(species, outer.id == cruise)) |>
  mutate(obs.ss.imp = exp(obs.lr + seasonal.mean),
         pred.ss = exp(pred.lr + seasonal.mean)) 

# export
write_rds(preds_outer_df, file = paste(out_dir, 'loo-preds-outer.rds', sep = ''))

# summarize
preds_outer_df |>
  group_by(species) |>
  summarize(outer.cor = cor(obs.ss, pred.ss),
            outer.rmspe = sqrt(mean((obs.ss - pred.ss)^2)),
            avg.inner.cor = mean(inner.cor),
            avg.inner.rmspe = mean(inner.rmspe))

