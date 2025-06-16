library(tidyverse)

# directories
data_dir <- 'data/processed/'
out_dir <- 'rslt/models/density/'

# read in scaled sighting data
paste(data_dir, 'density-estimates.RData', sep = '') |> load()
paste(data_dir, '_partitions/density-partitions.RData', sep = '') |> load()

## NAIVE PREDICTION 1: SEASONAL LAG METHOD -------------------------------------

# carry forward last observed seasonal value
naive_pred_lag <- dens_raw |>
  select(cruise, season, species, estimate) |>
  rename(obs.dens = estimate) |>
  arrange(species, season, cruise) |>
  group_by(season) |>
  mutate(pred.dens = lag(obs.dens, n = 1)) |>
  drop_na() |>
  arrange(species, cruise) |>
  group_by(species) |>
  summarize(rmspe.lag = mean((obs.dens - pred.dens)^2) |> sqrt())

## NAIVE PREDICTION 2: SEASONAL MEAN METHOD ------------------------------------

mean_preds <- loo_dens |>
  select(test.id, test.season, seasonal.means) |>
  unnest(seasonal.means) |>
  filter(season == test.season) |>
  mutate(across(ends_with('mean'), exp)) |>
  rename(bp = log.bp.mean,
         bm = log.bm.mean,
         mn = log.mn.mean) |>
  pivot_longer(c(bp, bm, mn), 
               names_to = 'species', 
               values_to = 'pred.dens') |>
  select(test.id, species, pred.dens) |>
  rename(cruise = test.id) 

naive_pred_mean <- dens_raw |>
  select(cruise, season, species, estimate) |>
  rename(obs.dens = estimate) |>
  left_join(mean_preds) |>
  arrange(species, cruise) |>
  group_by(species) |>
  summarize(rmspe.mean = mean((obs.dens - pred.dens)^2) |> sqrt())

## COMBINE AND EXPORT ----------------------------------------------------------

inner_join(naive_pred_lag, naive_pred_mean) |>
  write_rds(file = paste(out_dir, 'naive-preds.rds', sep = ''))
