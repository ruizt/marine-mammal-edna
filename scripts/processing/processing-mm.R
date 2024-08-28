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

# export processed data
save(list = c('sightings_raw', 'sightings', 'ss_means'), 
     file = paste(out_dir, 'mm-sightings.RData', sep = ''))

## LEAVE ONE OUT PARTITIONS ----------------------------------------------------

# generate leave one out partitions from imputed data
loo_raw <- ss_imputed |>
  select(cruise, season, ends_with('imp')) |>
  rename_with(~str_remove(.x, '.imp')) |>
  crossv_loo() |>
  mutate(train = map(train, as.data.frame),
         test = map(test, as.data.frame),
         test.id = map(test, ~pull(.x, cruise)),
         test.season = map(test, ~pull(.x, season))) |>
  unnest(c(test.id, test.season)) |>
  select(test.id, test.season, train, test)

# compute seasonal means from training partitions
loo_ss_means <- loo_raw |>
  unnest(train) |>
  group_by(test.id, test.season, season) |>
  summarize(across(c(bp, bm, mn), 
                   .fns = list(mean = ~mean(log(.x), na.rm = T)),
                   .names = 'log.{.col}.{.fn}'),
            .groups = 'drop') |>
  nest(seasonal.means = c(season, starts_with('log')))

# adjust for seasonal averages
loo_sightings <- loo_raw |>
  left_join(loo_ss_means, join_by(test.id, test.season)) |>
  mutate(train = map2(train, seasonal.means, detrend_fn),
         test = map2(test, seasonal.means, detrend_fn)) |>
  select(test.id, test.season, seasonal.means, train, test)

# export validation partitions
cv_dir <- paste(out_dir, '_cv/', sep = '')
fs::dir_create(cv_dir)
save(list = 'loo_sightings', 
     file = paste(cv_dir, 'mm-sightings-partitions.RData', sep = ''))


## NESTED LEAVE ONE OUT PARTITIONS ---------------------------------------------

# generate nested leave one out partitions
loo_nested_raw <- ss_imputed |>
  select(cruise, season, ends_with('imp')) |>
  rename_with(~str_remove(.x, '.imp')) |>
  crossv_loo(id = 'id.outer') |>
  rename(test.outer = test) |>
  mutate(cv.inner = map(train, ~crossv_loo(as.data.frame(.x), id = 'id.inner'))) |>
  select(-train) |>
  unnest(cv.inner) |>
  mutate(train.raw = map(train, as.data.frame),
         outer.id = map(test.outer, ~as_tibble(.x) |> pull(cruise)),
         inner.id = map(test, ~as_tibble(.x) |> pull(cruise))) |>
  unnest(ends_with('.id')) |>
  select(outer.id, inner.id, train.raw)

# compute seasonal means from training partitions
loo_means_nested <- loo_nested_raw |>
  unnest(train.raw) |>
  group_by(outer.id, inner.id, season) |>
  summarize(across(c(bp, bm, mn), 
                   .fns = list(mean = ~mean(log(.x), na.rm = T)),
                   .names = 'log.{.col}.{.fn}'),
            .groups = 'drop') |>
  nest(seasonal.means = c(season, starts_with('log')))

# adjust for seasonal averages
loo_sightings_nested <- loo_nested_raw |>
  left_join(loo_means_nested, 
            join_by(inner.id, outer.id)) |>
  mutate(train = map2(train.raw, seasonal.means, detrend_fn)) |>
  select(-train.raw, -seasonal.means)

# export validation partitions
cv_dir <- paste(out_dir, '_cv/', sep = '')
fs::dir_create(cv_dir)
save(list = 'loo_sightings_nested', 
     file = paste(cv_dir, 'mm-sightings-nested-partitions.RData', sep = ''))
