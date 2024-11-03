library(tidyverse)

# directories
data_dir <- 'data/processed/'
out_dir <- 'rslt/models/scaled-sightings/'

# read in scaled sighting data
paste(data_dir, 'mm-sightings.RData', sep = '') |> load()
paste(data_dir, '_partitions/mm-sightings-partitions.RData', sep = '') |> load()

## NAIVE PREDICTION 1: SEASONAL LAG METHOD -------------------------------------

# carry forward last observed seasonal value
naive_pred_lag <- sightings_raw |>
  pivot_longer(-c(cruise, season), 
               names_to = 'species', 
               values_to = 'obs.ss') |>
  arrange(species, season, cruise) |>
  group_by(season) |>
  mutate(pred.ss = lag(obs.ss, n = 1)) |>
  drop_na() |>
  arrange(species, cruise) |>
  group_by(species) |>
  summarize(rmspe.lag = mean((obs.ss - pred.ss)^2) |> sqrt())

## NAIVE PREDICTION 2: SEASONAL MEAN METHOD ------------------------------------

mean_preds <- loo_sightings |>
  select(test.id, test.season, seasonal.means) |>
  unnest(seasonal.means) |>
  filter(season == test.season) |>
  mutate(across(ends_with('mean'), exp)) |>
  rename(bp = log.bp.mean,
         bm = log.bm.mean,
         mn = log.mn.mean) |>
  pivot_longer(c(bp, bm, mn), names_to = 'species', values_to = 'pred.ss') |>
  select(test.id, species, pred.ss) |>
  rename(cruise = test.id) 

naive_pred_mean <- sightings_raw |>
  pivot_longer(-c(cruise, season), 
               names_to = 'species', 
               values_to = 'obs.ss') |>
  left_join(mean_preds) |>
  arrange(species, cruise) |>
  group_by(species) |>
  summarize(rmspe.mean = mean((obs.ss - pred.ss)^2) |> sqrt())

## COMBINE AND EXPORT ----------------------------------------------------------

inner_join(naive_pred_lag, naive_pred_mean) |>
  write_rds(file = paste(out_dir, 'naive-preds.rds', sep = ''))
