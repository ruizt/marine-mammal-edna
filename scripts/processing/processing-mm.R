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

# function to de-trend data for leave-one-out partitions
detrend_fn <- function(.data, .means){
left_join(.data, .means, join_by(season)) |>
  mutate(cruise = cruise,
         bp = log(bp) - log.bp.mean,
         bm = log(bm) - log.bm.mean,
         mn = log(mn) - log.mn.mean,
         .keep = "none")
}

# seasonal means for "outer" leave-one-out procedure
loo_ss_means_outer <- ss_imputed |>
  select(cruise, season, ends_with('imp')) |>
  rename_with(~str_remove(.x, '.imp')) |>
  crossv_loo() |>
  mutate(train = map(train, as.data.frame),
         test = map(test, as.data.frame),
         test.cruise = map(test, ~pull(.x, cruise)),
         test.season = map(test, ~pull(.x, season))) |>
  unnest(c(train, test.cruise, test.season)) |>
  group_by(test.cruise, test.season, season) |>
  summarize(across(c(bp, bm, mn), 
                   .fns = list(mean = ~mean(log(.x), na.rm = T)),
                   .names = 'log.{.col}.{.fn}'),
            .groups = 'drop') |>
  nest(seasonal.means = c(season, starts_with('log')))

# data partitions for "outer" leave one out validation
loo_sightings_outer <- ss_imputed |>
  select(cruise, season, ends_with('imp')) |>
  rename_with(~str_remove(.x, '.imp')) |>
  crossv_loo() |>
  mutate(train = map(train, as.data.frame),
         test = map(test, as.data.frame),
         test.cruise = map(test, ~pull(.x, cruise)),
         test.season = map(test, ~pull(.x, season))) |>
  unnest(c(test.cruise, test.season)) |>
  left_join(loo_ss_means_outer, join_by(test.cruise, test.season)) |>
  mutate(train = map2(train, seasonal.means, detrend_fn),
         test = map2(test, seasonal.means, detrend_fn)) |>
  select(test.cruise, test.season, seasonal.means, train, test)

# seasonal means for "inner" leave-one-out partitions
loo_ss_means_inner <- ss_imputed |>
  select(cruise, season, ends_with('imp')) |>
  rename_with(~str_remove(.x, '.imp')) |>
  crossv_loo(id = 'id.outer') |>
  rename(test.outer = test) |>
  mutate(cv.inner = map(train, ~crossv_loo(as.data.frame(.x), id = 'id.inner'))) |>
  select(-train) |>
  unnest(cv.inner) |>
  mutate(train = map(train, as.data.frame),
         test = map(test, as.data.frame),
         test.outer = map(test.outer, as.data.frame),
         test.inner.cruise = map(test, ~pull(.x, cruise)),
         test.inner.season = map(test, ~pull(.x, season)),
         test.outer.cruise = map(test.outer, ~pull(.x, cruise))) |>
  unnest(c(train, test.outer.cruise, 
           test.inner.cruise, test.inner.season)) |>
  group_by(test.inner.cruise, test.inner.season,
           test.outer.cruise, season) |>
  summarize(across(c(bp, bm, mn), 
                   .fns = list(mean = ~mean(log(.x), na.rm = T)),
                   .names = 'log.{.col}.{.fn}'),
            .groups = 'drop') |>
  nest(seasonal.means = c(season, starts_with('log'))) |>
  select(test.outer.cruise, test.inner.cruise, test.inner.season, seasonal.means)

# data partitions for "inner" leave one out validation
loo_sightings_inner <- ss_imputed |>
  select(cruise, season, ends_with('imp')) |>
  rename_with(~str_remove(.x, '.imp')) |>
  crossv_loo(id = 'id.outer') |>
  rename(test.outer = test) |>
  mutate(cv.inner = map(train, ~crossv_loo(as.data.frame(.x), id = 'id.inner'))) |>
  select(-train) |>
  unnest(cv.inner) |>
  mutate(train = map(train, as.data.frame),
         test = map(test, as.data.frame),
         test.outer = map(test.outer, as.data.frame),
         test.inner.cruise = map(test, ~pull(.x, cruise)),
         test.inner.season = map(test, ~pull(.x, season)),
         test.outer.cruise = map(test.outer, ~pull(.x, cruise))) |>
  unnest(c(test.outer.cruise, test.inner.cruise, test.inner.season)) |>
  left_join(loo_ss_means_inner, 
            join_by(test.outer.cruise, test.inner.cruise, test.inner.season)) |>
  mutate(train.inner = map2(train, seasonal.means, detrend_fn),
         test.inner = map2(test, seasonal.means, detrend_fn)) |>
  select(test.outer.cruise, test.inner.cruise, test.inner.season, 
         seasonal.means, train.inner, test.inner)


# export processed data
save(list = c('sightings_raw', 'sightings', 'ss_means'), 
     file = paste(out_dir, 'mm-sightings.RData', sep = ''))

# export validation partitions
val_dir <- paste(out_dir, '_cv-partitions/', sep = '')
fs::dir_create(val_dir)
save(list = c('loo_sightings_outer'), 
     file = paste(val_dir, 'mm-sightings-cv-outer.RData', sep = ''))
save(list = c('loo_sightings_inner'), 
     file = paste(val_dir, 'mm-sightings-cv-inner.RData', sep = ''))
