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
         year = year(datetime),
         season = factor(season, levels = c('winter', 'spring', 'summer', 'fall'))) |>
  select(datetime, cruise, year, season, species, declat, declong, best) |>
  filter(year >= 2014) 


effort_dates <- paste(in_dir, 'CalCOFI_2004-2025_Effort_OnTransectOnEffortONLY.xlsx', sep = '') |>
  readxl::read_xlsx() |>
  rename_with(~str_remove_all(.x, '[:punct:]') |> 
                str_squish() |> 
                str_replace_all(' ', '.') |> 
                tolower()) |>
  mutate(start.date = date(startdatetimelocal),
         end.date = date(enddatetimelocal),
         cruise = str_remove_all(cruise, '-') |> 
           str_trunc(4, side = 'left', ellipsis = '') %>%
           str_c('CC', .)) |>
  filter(cruise %in% unique(sightings$cruise)) |>
  distinct(year, season, cruise, start.date, end.date) |>
  group_by(cruise) |>
  summarize(start.date = min(start.date, na.rm = T),
            end.date = max(end.date, na.rm = T),
            .groups = 'drop')

save(list = c('sightings', 'effort_dates'), file = paste(out_dir, 'sightings.RData', sep = ''))
