library(tidyverse)
load('data/edna_16s_processed.RData')
mammal_raw <- read_csv('data/CalCOFI_2004-2021_CombinedSightings.csv')

# name repair
names <- colnames(mammal_raw)
newnames <- names %>%
  str_remove_all('[:punct:]') %>%
  str_squish() %>%
  tolower() %>%
  str_replace_all(' ', '.')
colnames(mammal_raw) <- newnames

# counts of individuals/sightings by cruise
event_counts <- mammal_raw %>% 
  select(cruise, season, datetime.local, declat, declong, unique.id, 
         adjusted.both.on.effort.and.on.transect, species.1, best) %>%
  filter(adjusted.both.on.effort.and.on.transect == 'ON') %>%
  mutate(across(starts_with('species'), tolower)) %>%
  filter(species.1 %in% c('mn', 'bp', 'bm')) %>%
  mutate(date = lubridate::mdy_hm(datetime.local),
         year = year(date)) %>%
  group_by(season, year, species.1) %>%
  summarize(n.individuals = sum(best),
            n.sightings = n(),
            .groups = 'drop') 

# individuals
event_counts %>%
  filter(species.1 == 'bm') %>%
  select(-n.sightings, -species.1) %>%
  pivot_wider(names_from = season, 
              values_from = starts_with('n'))

# sightings
event_counts %>%
  filter(species.1 == 'bm') %>%
  select(-n.individuals, -species.1) %>%
  pivot_wider(names_from = season, 
              values_from = starts_with('n'))

# cruises
mammal_cruise_ids <- mammal_raw %>%
  filter(adjusted.both.on.effort.and.on.transect == 'ON') %>%
  mutate(across(starts_with('species'), tolower)) %>%
  filter(species.1 %in% c('mn', 'bp', 'bm')) %>%
  pull(cruise) %>%
  unique() %>%
  sort() %>%
  str_remove('CC')

edna_cruise_ids <- edna_data %>%
  pull(cruise) %>%
  str_trunc(width = 4, side = 'left', ellipsis = '')

# 24 common cruises
intersect(edna_cruise_ids, mammal_cruise_ids) %>% length()

setdiff(edna_cruise_ids, mammal_cruise_ids)
