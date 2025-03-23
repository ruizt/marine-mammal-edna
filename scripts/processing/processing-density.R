library(tidyverse)
library(lubridate)
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

## DENSITY ESTIMATES -----------------------------------------------------------

# import whale density estimates ( y_{i} )
dens_in <- lapply(c('BM', 'BP', 'MN'),
       function(.sp){
         paste(in_dir, .sp, '.csv', sep = '') |>
           read_csv() |>
           mutate(species = .sp) |>
           filter(Label != 'Total')
       }) |>
  bind_rows() |>
  mutate(cruise = str_remove(Label, 'CC'),
         species = tolower(species),
         yr = ym(cruise) |> year(),
         qtr = ym(cruise) |> quarter()) |>
  rename_with(tolower) |>
  select(cruise, yr, qtr, species, estimate, se, cv, lcl, ucl, df) 

inspect <- dens_in |>
  count(yr, qtr) |>
  filter(n != 3)

qtr.replace <- inner_join(dens_in, inspect, join_by(yr, qtr)) |>
  distinct(cruise, yr, qtr) |>
  arrange(cruise) |>
  group_by(yr, qtr) |>
  slice_max(cruise) |>
  mutate(qtr.rep = qtr + 1) |>
  ungroup()

dens_raw <- left_join(dens_in, qtr.replace, join_by(cruise, yr, qtr)) |>
  mutate(season = if_else(is.na(qtr.rep), qtr, qtr.rep) |> 
           factor(labels = c('winter', 'spring', 'summer', 'fall')),
         cruise = str_remove(cruise, '-')) |>
  select(-qtr.rep) |>
  filter(cruise %in% cruises_of_interest)

# impute zeroes with uniform random numbers up to seasonal minima
set.seed(30924)
dens_wide <- dens_raw |>
  select(cruise, yr, qtr, season, species, estimate) |>
  pivot_wider(names_from = species, values_from = estimate) 

dens_imputed <- dens_wide |>
  group_by(season) |>
  summarize(across(c(bm, bp, mn), 
                   .fns = list(min = ~min(.x[.x > 0], na.rm = T)),
                   .names = '{.col}.{.fn}')) |>
  right_join(dens_wide, join_by(season)) |>
  mutate(bp.imp = if_else(bp == 0, 
                          runif(n = length(bp == 0), min = 0, max = bp.min),
                          bp),
         bm.imp = if_else(bm == 0,
                          runif(n = length(bm == 0), min = 0, max = bm.min),
                          bm),
         mn.imp = if_else(mn == 0,
                          runif(n = length(mn == 0), min = 0, max = mn.min),
                          mn)) |>
  select(cruise, yr, qtr, season, bm, bp, mn, bm.imp, bp.imp, mn.imp)

# seasonal means ( log[g_{y}(i)] )
dens_means <- dens_imputed |>
  dplyr::select(cruise, season, ends_with('imp')) |>
  group_by(season) |>
  summarize(across(ends_with('imp'), 
                   .fns = list(mean = ~mean(log(.x), na.rm = T)),
                   .names = 'log.{.col}.{.fn}'))

# de-trend imputed densities ( log[y_{i}/g_{y}(i)] )
dens <- dens_means |>
  right_join(dens_imputed, join_by(season)) |>
  mutate(cruise = cruise,
         bp = log(bp.imp) - log.bp.imp.mean,
         bm = log(bm.imp) - log.bm.imp.mean,
         mn = log(mn.imp) - log.mn.imp.mean,
         .keep = "none")

# export processed data
save(list = c('dens_raw', 'dens', 'dens_means'), 
     file = paste(out_dir, 'density-estimates.RData', sep = ''))

## LEAVE ONE OUT PARTITIONS ----------------------------------------------------

# generate leave one out partitions from imputed data
loo_raw <- dens_imputed |>
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
loo_dens_means <- loo_raw |>
  unnest(train) |>
  group_by(test.id, test.season, season) |>
  summarize(across(c(bp, bm, mn), 
                   .fns = list(mean = ~mean(log(.x), na.rm = T)),
                   .names = 'log.{.col}.{.fn}'),
            .groups = 'drop') |>
  nest(seasonal.means = c(season, starts_with('log')))

# # variation in seasonal means across partitions
# loo_dens_means |>
#   unnest(seasonal.means) |>
#   group_by(season) |>
#   summarize(across(ends_with('.mean'), sd)) |>
#   mutate(across(ends_with('.mean'), exp))

# function to de-trend data for leave-one-out partitions
detrend_fn <- function(.data, .means){
  left_join(.data, .means, join_by(season)) |>
    mutate(cruise = cruise,
           bp = log(bp) - log.bp.mean,
           bm = log(bm) - log.bm.mean,
           mn = log(mn) - log.mn.mean,
           .keep = "none")
}

# adjust for seasonal averages
loo_dens <- loo_raw |>
  left_join(loo_dens_means, join_by(test.id, test.season)) |>
  mutate(train = map2(train, seasonal.means, detrend_fn),
         test = map2(test, seasonal.means, detrend_fn)) |>
  select(test.id, test.season, seasonal.means, train, test)

# export validation partitions
cv_dir <- paste(out_dir, '_partitions/', sep = '')
fs::dir_create(cv_dir)
save(list = 'loo_dens', 
     file = paste(cv_dir, 'density-partitions.RData', sep = ''))


## NESTED LEAVE ONE OUT PARTITIONS ---------------------------------------------

# generate nested leave one out partitions
loo_nested_raw <- dens_imputed |>
  select(cruise, season, ends_with('imp')) |>
  rename_with(~str_remove(.x, '.imp')) |>
  crossv_loo(id = 'id.outer') |>
  rename(test.outer = test) |>
  mutate(cv.inner = map(train, ~crossv_loo(as.data.frame(.x), id = 'id.inner'))) |>
  select(-train) |>
  unnest(cv.inner) |>
  mutate(train.raw = map(train, as.data.frame),
         test.raw = map(test, as.data.frame),
         outer.id = map(test.outer, ~as_tibble(.x) |> pull(cruise)),
         inner.id = map(test, ~as_tibble(.x) |> pull(cruise))) |>
  unnest(ends_with('.id')) |>
  select(outer.id, inner.id, train.raw, test.raw) 

# compute seasonal means from training partitions
loo_means_nested <- loo_nested_raw |>
  unnest(train.raw) |>
  group_by(outer.id, inner.id, season) |>
  summarize(across(c(bp, bm, mn), 
                   .fns = list(mean = ~mean(log(.x), na.rm = T)),
                   .names = 'log.{.col}.{.fn}'),
            .groups = 'drop') |>
  nest(seasonal.means = c(season, starts_with('log')))

# function to de-trend data for leave-one-out partitions
detrend_fn <- function(.data, .means){
  left_join(.data, .means, join_by(season)) |>
    mutate(cruise = cruise,
           bp = log(bp) - log.bp.mean,
           bm = log(bm) - log.bm.mean,
           mn = log(mn) - log.mn.mean,
           .keep = "none") 
}

# adjust for seasonal averages
loo_dens_nested <- loo_nested_raw |>
  left_join(loo_means_nested, 
            join_by(inner.id, outer.id)) |>
  mutate(train = map2(train.raw, seasonal.means, detrend_fn), 
         test = map2(test.raw, seasonal.means, detrend_fn)) |>
  select(-train.raw, -test.raw)

# export validation partitions
cv_dir <- paste(out_dir, '_partitions/', sep = '')
fs::dir_create(cv_dir)
save(list = 'loo_dens_nested', 
     file = paste(cv_dir, 'density-nested-partitions.RData', sep = ''))
