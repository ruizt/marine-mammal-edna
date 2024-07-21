library(tidyverse)
out_dir <- 'data/processed/'
in_dir <- 'data/_raw/'

## SCALED SIGHTING DATA --------------------------------------------------------

# import whale sighting data ( y_{i} )
scaled_sightings <- paste(in_dir, "CC_on_effort_scaled_sightings_Ceta_SCB.csv", sep = '') |>
  read_csv() |>
  rename_with(tolower) |>
  rename_with(~gsub('_.*', '', .x)) |>
  mutate(cruise = str_replace(cruiseid, "CC", "20")) |>
  dplyr::select(cruise, season, bp, bm, mn)

# impute zeroes with uniform random numbers up to seasonal minima (no notation yet)
set.seed(30924)
ss_imputed <- scaled_sightings |>
  group_by(season) |>
  summarize(across(where(is.numeric), 
                   .fns = list(min = ~min(.x[.x > 0], na.rm = T)),
                   .names = '{.col}.{.fn}')) |>
  right_join(scaled_sightings, by = 'season') |>
  mutate(bp.imp = if_else(bp == 0, 
                          runif(n = length(bp == 0), min = 0, max = bp.min),
                          bp),
         bm.imp = if_else(bm == 0,
                          runif(n = length(bm == 0), min = 0, max = bm.min),
                          bm),
         mn.imp = if_else(mn == 0,
                          runif(n = length(mn == 0), min = 0, max = mn.min),
                          mn))

# seasonal de-trending ( log[y_{i}/g_{y}(i)] )
sightings <- ss_imputed |>
  dplyr::select(cruise, season, ends_with('imp')) |>
  group_by(season) |>
  summarize(across(ends_with('imp'), 
                   .fns = list(mean = ~mean(log(.x), na.rm = T)),
                   .names = 'log.{.col}.{.fn}')) |>
  right_join(ss_imputed) |>
  mutate(cruise = cruise,
         bp = log(bp.imp) - log.bp.imp.mean,
         bm = log(bm.imp) - log.bm.imp.mean,
         mn = log(mn.imp) - log.mn.imp.mean,
         .keep = "none")

# export
save(sightings, 
     file = paste(out_dir, 'mm-sightings-', today(), '.RData', sep = ''))
