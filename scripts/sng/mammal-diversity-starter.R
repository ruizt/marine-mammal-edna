library(tidyverse)
library(vegan)
load('data/ncog-16s.RData')

## DATASETS

# read in mammal data
mammal <- read_csv('data/CC_on_effort_scaled_sightings_Ceta_SCB.csv') %>%
  mutate(cruise = str_replace(cruiseID, 'CC', 'X20')) %>%
  select(cruise, Bp_scaled, Bm_scaled, Mn_scaled) %>%
  rename(bp = Bp_scaled,
         bm = Bm_scaled,
         mn = Mn_scaled)

# same preprocessing as before for 16s edna
# common transects sampled across cruises (filter to these for common area) 
sampling_counts <- metadata %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') %>%
  count(Cruise, line) %>%
  mutate(line = as.numeric(line)) %>%
  filter((76 < line) & (line < 93.4)) 

# transect IDs
transects <- sampling_counts %>% pull(line) %>% unique()

# filter to common transects 
edna16s_reads <- edna_samples %>%
  filter(as.numeric(line) %in% transects)

# split into sample info and asvs
asv <- edna16s_reads %>%
  select(starts_with('asv'))
sampinfo <- edna16s_reads %>%
  select(-starts_with('asv'))

# alpha diversity measures
alpha_div <- sampinfo %>%
  bind_cols(alpha.div.sh = diversity(asv, index = 'shannon'))

## EXPLORATORY

# join mammal and diversity data
div_mammal <- alpha_div %>%
  filter(depthm < 150) %>%
  mutate(depth = cut_number(as.numeric(depthm), n = 2)) %>%
  group_by(cruise, depth) %>%
  summarize(alpha.mean = mean(alpha.div.sh),
            alpha.sd = sd(alpha.div.sh)) %>%
  left_join(mammal, by = 'cruise') %>%
  ungroup()

# diversity against bp by depth
div_mammal %>%
ggplot(aes(x = alpha.mean, y = bp)) +
  geom_point() +
  facet_wrap(~depth)
