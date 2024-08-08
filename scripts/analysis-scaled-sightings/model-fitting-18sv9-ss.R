library(tidyverse)
library(spls)
library(pls)
library(modelr)
library(magrittr)
library(fs)
library(collapse)
data_dir <- 'data/processed/'
loo_dir <- 'rslt/loocv/18sv9-ss/'
out_dir <- 'rslt/models/scaled-sightings/'
dir_create(out_dir)

## DATA INPUTS -----------------------------------------------------------------

# read in and combine 18sv9 data and scaled sightings
paste(data_dir, 'ncog18sv9.RData', sep = '') |> load()
paste(data_dir, 'mm-sightings-2024-07-27.RData', sep = '') |> load() 
whales <- inner_join(sightings, edna, by = 'cruise') 

# dimensions
p <- whales |> dplyr::select(starts_with('asv')) |> ncol()
n.obs <- nrow(whales)

## CONSTRUCT CANDIDATE STABLE SETS ---------------------------------------------

# read in selection frequencies and fit metrics from LOOCV
sel_freq <- paste(loo_dir, 'selection-frequencies.rds', sep = '') |> read_rds()
metrics <- paste(loo_dir, 'metrics.rds', sep = '') |> read_rds()
loo_partitions <- paste(loo_dir, 'partitions.rds', sep = '') |> read_rds()

# chosen stability threshold (minimax selection prob.)
pi.max <- 0.8

# criterion: range of ncomp/eta to search should have an average number
# of selected asvs not exceeding 100; controls EV at 1
EV.max <- 0.5
q.max <- sqrt((2*pi.max - 1)*p*EV.max)
q.max

# criterion: range to search should not include regions having models with too
# few asv's selected
q.min <- 10

# recover eta and ncomp unique values from loo metrics
eta_grid <- metrics |> pull(eta) |> unique() |> sort()
ncomp_grid <- metrics |> pull(ncomp) |> unique() |> sort()

# generate sets of candidate eta intervals
eta_intervals <- expand_grid(eta.min = seq_range(eta_grid, by = 0.025),
                             eta.max = rev(eta.min)) |>
  filter(eta.max - eta.min > 0.1) 

# function to compute stable sets from a parameter (eta/ncomp) range
stableset_fn <- function(parms, .species){
  out <- sel_freq |> 
    filter(eta >= parms$eta.min,
           eta <= parms$eta.max,
           ncomp == parms$ncomp,
           n >= n.obs*pi.max,
           species == .species) |>
    distinct(asv) |>
    pull()
  
  return(out)
}

# identify largest stable sets maintaining EV control for each ncomp
candidate_sets <- metrics |>
  fgroup_by(species, ncomp, eta) |>
  fselect(n.asv) |>
  fmean() |>
  nest(metrics = -c(eta, ncomp)) |>
  mutate(interval = list(eta_intervals)) |>
  unnest(cols = interval) |>
  filter(eta >= eta.min, eta <= eta.max) |>
  unnest(cols = metrics) |>
  fgroup_by(species, ncomp, eta.min, eta.max) |>
  fsummarize(q = fmean(n.asv),
             sd = fsd(n.asv)) |>
  filter(q < q.max, q > q.min, sd/q < 1) |>
  mutate(setting = row_number()) |>
  nest(parms = c(ncomp, eta.min, eta.max)) |>
  mutate(ss = map2(parms, species, stableset_fn)) |>
  unnest(parms) |>
  # group_by(species, ncomp) |>
  # slice_max(ss.size) |>
  group_by(species, ncomp, ss) |>
  slice_min(eta.max - eta.min) |>
  ungroup()

# inspect
candidate_sets


## DETERMINE NUMBER OF LATENT COMPONENTS ---------------------------------------

# compute leave one out predictions for each candidate set
## LONG RUNTIME: ~10-15min
loo_preds <- loo_partitions |>
  expand_grid(candidate_sets) |>
  filter(ncomp < map(ss, length)) |>
  mutate(train = map2(train, species, ~select(.x, {.y}, starts_with('asv')) |>
                        rename(y = {.y})),
         test = map2(test, species, ~select(.x, {.y}, starts_with('asv')) |>
                       rename(y = {.y}))) |>
  mutate(train = map2(train, ss, ~select(.x, y, any_of(.y))),
         test = map2(test, ss, ~select(.x, y, any_of(.y)))) |>
  mutate(fit = map2(train, ncomp, ~plsr(y ~ ., data = .x, ncomp = .y, 
                                        scale = F, center = T)),
         pred.out = map2(fit, test, ~predict(.x, .y)),
         pred = map2(ncomp, pred.out, ~.y[, , paste(.x, 'comps')]),
         y = map(test, ~pull(.x, y))) |>
  select(-pred.out) |>
  unnest(c(pred, y)) 

# identify prediction-optimal settings
best_settings <- loo_preds |> 
  group_by(species, setting) |>
  summarize(mspe = mean((y - pred)^2),
            cor = cor(y, pred)) |>
  slice_min(mspe)

# retreive prediction-optimal stable sets and hyperparameter settings
best_ss <- candidate_sets |>
  left_join(best_settings, join_by(species), suffix = c('', '.best')) |>
  filter(setting == setting.best) |>
  select(species, ncomp, eta.min, eta.max, ss, mspe, cor)
best_ss

## FIT MODELS ON STABLE SETS ---------------------------------------------------
  
# function to retrieve data given a stable set
ssdata_fn <- function(stableset, species){
  whales |>
    rename(y = {{species}}) |>
    select(y, all_of(stableset))
}

# function to fit model given data and number of components
fit_fn <- function(.data, .ncomp){
  out <- plsr(y ~ ., data = .data, ncomp = .ncomp,
              scale = F, center = T)
  return(out)
}


# fit models to stable sets
fitted_models <- best_ss |>
  select(species, ncomp, ss) |>
  mutate(data = map2(ss, species, ssdata_fn),
         fit = map2(data, ncomp, fit_fn)) 

# inspect
fitted_models

## MODEL SUMMARIES -------------------------------------------------------------

# rearrange seasonal means for join with fitted values
ss_means_long <- ss_means |>
  pivot_longer(-season, names_to = 'species', values_to = 'seasonal.mean') |>
  mutate(species = str_remove(species, 'log.') |> str_remove('.imp.mean'))

# rearrange raw scaled sighting data for join with fitted values
sightings_raw_long <- sightings_raw |>
  pivot_longer(c(bp, bm, mn), names_to = 'species', values_to = 'ss.obs')

# compute fitted values and back-transform to fitted scaled sightings
fit_df <- fitted_models |>
  mutate(y = map(data, ~pull(.x, y)),
         fitted = map2(fit, ncomp, ~fitted(.x)[, , paste(.y, 'comps')]),
         cruise = pull(whales, cruise) |> list()) |>
  select(cruise, species, ncomp, y, fitted) |>
  unnest(everything()) |>
  left_join(sightings_raw_long, join_by(species, cruise)) |>
  left_join(ss_means_long, join_by(season, species)) |>
  mutate(ss.fit = exp(fitted + seasonal.mean),
         lr.resid = y - fitted,
         ss.resid = ss.obs - ss.fit) 

# compute adjusted r squared for both logratio scale and original scale
## (a bit awkward on original scale but familiar statistic)
fit_metrics <- fit_df |>
  group_by(species) |>
  summarize(adj.rsq.lr = (1 - ((n.obs - 1)/(n.obs - unique(ncomp) - 1))*var(lr.resid)/var(y)),
            adj.rsq.ss = (1 - ((n.obs - 1)/(n.obs - unique(ncomp) - 1))*var(ss.resid)/var(ss.obs))) |>
  left_join(fitted_models, join_by(species)) |>
  mutate(n.asv = map(ss, length)) |>
  unnest(n.asv) |>
  select(species, n.asv, adj.rsq.lr, adj.rsq.ss)

# inspect
fit_metrics

## PREDICTIVE PERFORMANCE ------------------------------------------------------

# retrieve leave-one-out predictions for each model from prior computations
loo_preds_best <- loo_preds |>
  left_join(best_settings, join_by(species), suffix = c('', '.best')) |>
  filter(setting == setting.best) |>
  select(test.cruise, species, y, pred)

# back-transform predictions to original scale
loo_pred_df <- loo_sightings |> 
  select(test.cruise, test.season, seasonal.means) |>
  unnest(seasonal.means) |>
  filter(test.season == season) |>
  pivot_longer(starts_with('log'), 
               names_to = 'test.species',
               values_to = 'train.seasonal.mean') |>
  mutate(test.species = str_remove(test.species, 'log.') |> str_remove('.mean')) |>
  select(-season) |>
  inner_join(loo_preds_best, join_by(test.cruise, test.species == species)) |>
  mutate(ss.pred = exp(pred + train.seasonal.mean),
         ss.imp = exp(y + train.seasonal.mean)) |>
  inner_join(sightings_raw_long,
             join_by(test.cruise == cruise, test.species == species)) 

# compute prediction metrics  
pred_metrics <- loo_pred_df |>
  group_by(test.species) |>
  summarize(rmspe.lr = mean((y - pred)^2) |> sqrt(),
            rmspe.lr.naive = mean((y - train.seasonal.mean)^2) |> sqrt(),
            cor.lr = cor(y, pred),
            rmspe.ss = mean((ss.obs - ss.pred)^2) |> sqrt(),
            cor.ss = cor(ss.obs, ss.pred),
            rmspe.ss.naive = mean((ss.obs - exp(train.seasonal.mean))^2) |> sqrt()) |>
  pivot_longer(-test.species) |>
  pivot_wider(names_from = test.species, values_from = value) |>
  separate(name, into = c('metric', 'scale', 'model')) |>
  mutate(model = replace_na(model, 'pls')) |>
  arrange(metric, scale, model)

# inspect
pred_metrics

# export
save(list = c('fitted_models', 
              'loo_preds',
              'loo_preds_best', 
              'loo_pred_df', 
              'fit_df', 
              'fit_metrics', 
              'pred_metrics'),
     file = paste(out_dir, 'fitted-models-18sv9-ss.RData', sep = ''))


## COMPARE WITH PREDICTION-OPTIMAL SPLS ----------------------------------------

# find prediction optimal spls models from loocv metrics
loo_best <- metrics |>
  group_by(ncomp, eta, species) |>
  summarize(mspe = mean(sq.pred.err),
            cor = cor(pred, pred + pred.err),
            df = mean(n.asv)) |>
  group_by(species) |>
  slice_min(mspe)

# retrieve prediction optimal hyperparameter configurations
parms <- loo_best |>
  dplyr::select(species, ncomp, eta) |>
  nest(parms = c(ncomp, eta))

# function to fit spls models
fit_spls_fn <- function(.data, .parms){
  x <- dplyr::select(.data, -y)
  y <- pull(.data, y)
  out <- spls(x, y, K = .parms$ncomp, eta = .parms$eta,
              scale.x = F, scale.y = F)
  return(out)
}

# fit models
fit_pos <- whales |>
  dplyr::select(-cruise) |>
  pivot_longer(c(bm, bp, mn), names_to = 'species', values_to = 'y') |>
  group_by(species) |>
  nest(data = -species) |>
  left_join(parms) |>
  mutate(fit = map2(data, parms, fit_spls_fn),
         df = map(fit, ~nrow(.x$projection)),
         fitted = map(fit, ~predict(.x, type = 'fit')[, 1]),
         y = map(data, ~pull(.x, y)),
         adj.rsq = map2(y, fitted, ~(1 - (24/16)*var(.x - .y)/var(.y))))

# compute fit metrics
loo_best |>
  mutate(rmse = sqrt(mspe),
         pred.corr = cor) |>
  dplyr::select(species, rmse, pred.corr) |>
  left_join(unnest(fit_pos, c(df, adj.rsq))) |>
  dplyr::select(species, rmse, df, adj.rsq, pred.corr)

# compare with stability selection
pred_metrics
