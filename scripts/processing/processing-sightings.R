library(tidyverse)
library(lubridate)
in_dir <- 'data/_raw/'
out_dir <- 'data/processed/'
fs::dir_create(out_dir)

sightings <- paste(in_dir, 'CalCOFI_2004-2021_CombinedSightings.csv', sep = '') |>
  read_csv() |>
  rename_with(~str_remove_all(.x, '[:punct:]') |> 
                str_squish() |> 
                str_replace_all(' ', '.') |> 
                tolower()) |>
  mutate(species = tolower(species.1)) |>
  filter(adjusted.both.on.effort.and.on.transect == 'ON',
         species %in% c('bm', 'bp', 'mn')) |>
  mutate(datetime = mdy_hm(datetime.local),
         year = year(datetime)) |>
  select(cruise, year, datetime, species, declat, declong, best) |>
  filter(year >= 2014) 

save(sightings, file = paste(out_dir, 'sightings.RData', sep = ''))
