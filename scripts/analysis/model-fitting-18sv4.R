library(tidyverse)
library(spls)
library(pls)
library(modelr)
library(fs)
library(collapse)
data_dir <- 'data/processed/'
stbl_dir <- 'rslt/stability-selection/18sv4/'
out_dir <- 'rslt/models/'
fs::dir_create(out_dir)

## DATA INPUTS -----------------------------------------------------------------

# read in and combine 18sv4 data and density estimates
paste(data_dir, 'ncog18sv4.RData', sep = '') |> load()
paste(data_dir, 'density-estimates.RData', sep = '') |> load() 
whales <- inner_join(dens, edna, by = 'cruise') 

# dimensions
p <- whales |> dplyr::select(starts_with('asv')) |> ncol()
n.obs <- nrow(whales)

## FIT MODELS ON STABLE SETS ---------------------------------------------------

# stable sets
stable_sets <- paste(stbl_dir, 'stable-sets.rds', sep = '') |> read_rds()

# function to retrieve data given a stable set
subset_fn <- function(stableset, species){
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
fitted_models <- stable_sets |>
  mutate(data = map2(ss, species, subset_fn),
         fit = map2(data, ncomp, fit_fn),
         coef = map2(fit, ncomp, ~coef(.x)[, , paste(.y, 'comps')])) |>
  select(species, ss, ncomp, data, fit, coef)

## MODEL METRICS ---------------------------------------------------------------

# rearrange seasonal means for join with fitted values
dens_means_long <- dens_means |>
  pivot_longer(-season, names_to = 'species', values_to = 'seasonal.mean') |>
  mutate(species = str_remove(species, 'log.') |> str_remove('.imp.mean'))

# rearrange density data for join with fitted values
dens_long <- dens_raw |>
  select(cruise, season, species, estimate) |>
  rename(dens.obs = estimate)

# compute fitted values and back-transform to fitted scaled sightings
fit_df <- fitted_models |>
  mutate(lr.obs = map(data, ~pull(.x, y)),
         lr.fit = map2(fit, ncomp, ~fitted(.x)[, , paste(.y, 'comps')]),
         cruise = pull(whales, cruise) |> list()) |>
  select(species, cruise, ncomp, lr.obs, lr.fit) |>
  unnest(everything()) |>
  left_join(dens_long, join_by(species, cruise)) |>
  left_join(dens_means_long, join_by(season, species)) |>
  mutate(dens.fit = exp(lr.fit + seasonal.mean),
         lr.resid = lr.obs - lr.fit,
         dens.resid = dens.obs - dens.fit) |>
  select(species, cruise, ncomp, starts_with('lr'), starts_with('dens'))

# compute adjusted r squared for both logratio scale and original scale
## (a bit awkward on original scale but familiar statistic)
model_metrics <- fit_df |>
  group_by(species) |>
  summarize(adj.rsq.lr = (1 - ((n.obs - 1)/(n.obs - unique(ncomp) - 1))*var(lr.resid)/var(lr.obs)),
            adj.rsq.dens = (1 - ((n.obs - 1)/(n.obs - unique(ncomp) - 1))*var(dens.resid)/var(dens.obs))) |>
  left_join(stable_sets, join_by(species)) |>
  mutate(n.asv = map(ss, length)) |>
  unnest(n.asv) |>
  select(species, n.asv, ends_with('.lr'), ends_with('.dens'))

# inspect
model_metrics

## EXPORT ----------------------------------------------------------------------

save(list = c('whales',
              'fitted_models',
              'fit_df',
              'model_metrics'),
     file = paste(out_dir, 'fitted-models-18sv4.RData', sep = ''))
