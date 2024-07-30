library(tidyverse)
library(modelr)
out_dir <- 'data/processed/'
in_dir <- 'data/_raw/'

# cruises present in edna data
cruises_of_interest <- paste(in_dir, "NCOG_sample_log_DNA_stvx_meta_2014-2020.csv", sep = '') |> 
  read_csv() |>
  rename_with(tolower) |>
  pull(cruise) |>
  unique() |>
  as.character()

## SCALED SIGHTING DATA --------------------------------------------------------

# import whale sighting data ( y_{i} )
sightings_raw <- paste(in_dir, "CC_on_effort_scaled_sightings_Ceta_SCB.csv", sep = '') |>
  read_csv() |>
  rename_with(tolower) |>
  rename_with(~gsub('_.*', '', .x)) |>
  mutate(cruise = str_replace(cruiseid, "CC", "20")) |>
  dplyr::select(cruise, season, bp, bm, mn) |>
  filter(cruise %in% cruises_of_interest)

# impute zeroes with uniform random numbers up to seasonal minima
set.seed(30924)
ss_imputed <- sightings_raw |>
  group_by(season) |>
  summarize(across(where(is.numeric), 
                   .fns = list(min = ~min(.x[.x > 0], na.rm = T)),
                   .names = '{.col}.{.fn}')) |>
  right_join(sightings_raw, by = 'season') |>
  filter(cruise %in% cruises_of_interest) |>
  mutate(bp.imp = if_else(bp == 0, 
                          runif(n = length(bp == 0), min = 0, max = bp.min),
                          bp),
         bm.imp = if_else(bm == 0,
                          runif(n = length(bm == 0), min = 0, max = bm.min),
                          bm),
         mn.imp = if_else(mn == 0,
                          runif(n = length(mn == 0), min = 0, max = mn.min),
                          mn))

# seasonal means ( log[g_{y}(i)] )
ss_means <- ss_imputed |>
  dplyr::select(cruise, season, ends_with('imp')) |>
  group_by(season) |>
  summarize(across(ends_with('imp'), 
                   .fns = list(mean = ~mean(log(.x), na.rm = T)),
                   .names = 'log.{.col}.{.fn}'))

# de-trend imputed sightings ( log[y_{i}/g_{y}(i)] )
sightings <- ss_means |>
  right_join(ss_imputed, join_by(season)) |>
  mutate(cruise = cruise,
         bp = log(bp.imp) - log.bp.imp.mean,
         bm = log(bm.imp) - log.bm.imp.mean,
         mn = log(mn.imp) - log.mn.imp.mean,
         .keep = "none")

# seasonal means from leave-one-out procedure
loo_ss_means <- ss_imputed |>
  select(cruise, season, ends_with('imp')) |>
  rename_with(~str_remove(.x, '.imp')) |>
  crossv_loo() |>
  mutate(train = map(train, as.data.frame),
         test = map(test, as.data.frame),
         test.cruise = map(test, ~pull(.x, cruise)),
         test.season = map(test, ~pull(.x, season))) |>
  unnest(c(train, test.cruise, test.season)) |>
  group_by(.id, test.cruise, test.season, season) |>
  summarize(across(c(bp, bm, mn), 
                   .fns = list(mean = ~mean(log(.x), na.rm = T)),
                   .names = 'log.{.col}.{.fn}'),
            .groups = 'drop') |>
  nest(seasonal.means = c(season, starts_with('log')))

# process data for leave-one-out partitions
detrend_fn <- function(.data, .means){
left_join(.data, .means, join_by(season)) |>
  mutate(cruise = cruise,
         bp = log(bp) - log.bp.mean,
         bm = log(bm) - log.bm.mean,
         mn = log(mn) - log.mn.mean,
         .keep = "none")
}

loo_sightings <- ss_imputed |>
  select(cruise, season, ends_with('imp')) |>
  rename_with(~str_remove(.x, '.imp')) |>
  crossv_loo() |>
  mutate(train = map(train, as.data.frame),
         test = map(test, as.data.frame),
         test.cruise = map(test, ~pull(.x, cruise)),
         test.season = map(test, ~pull(.x, season))) |>
  unnest(c(test.cruise, test.season)) |>
  left_join(loo_ss_means, join_by(.id, test.cruise, test.season)) |>
  mutate(train = map2(train, seasonal.means, detrend_fn),
         test = map2(test, seasonal.means, detrend_fn)) |>
  select(.id, test.cruise, test.season, seasonal.means, train, test)

# export
save(list = c('sightings_raw', 'sightings', 'ss_means', 'loo_sightings'), 
     file = paste(out_dir, 'mm-sightings-', today(), '.RData', sep = ''))
