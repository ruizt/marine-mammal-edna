# RUN ONCE, THEN COMMENT OUT
# source('scripts/sng/edna-preprocessing.R')
# rm(list = ls())

# script starts here
library(tidyverse)
library(vegan)
load('data/ncog-16s.RData')

## metadata summaries

# number of samples per cruise (27 cruises)
metadata %>% count(Cruise, name = 'n.samples')

# number of samples by transect and cruise (15 transects)
metadata %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') %>%
  count(Cruise, line) %>%
  arrange(line, Cruise) %>%
  pivot_wider(names_from = line, values_from = n) %>% 
  print(n = 27)

# number of samples by station, transect, and cruise (91 stations)
metadata %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') %>%
  count(Cruise, line, sta)

# number of samples by depth and cruise
metadata %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') %>%
  mutate(depth.rng = cut_interval(Depthm, n = 2)) %>%
  count(Cruise, depth.rng) %>%
  pivot_wider(names_from = depth.rng, values_from = n)

# common transects sampled across cruises (filter to these for common area) 
sampling_counts <- metadata %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') %>%
  count(Cruise, line) %>%
  mutate(line = as.numeric(line)) %>%
  filter((76 < line) & (line < 93.4)) 

# transect IDs
transects <- sampling_counts %>% pull(line) %>% unique()

# number of observations per transect/cruise for common transects
sampling_counts %>%
  arrange(line, Cruise) %>%
  pivot_wider(names_from = line, values_from = n) %>% 
  print(n = 27)

## 16s edna summaries

cruises <- edna_samples %>% pull(cruise) %>% unique()

# filter to common transects 
edna16s_reads <- edna_samples %>%
  filter(as.numeric(line) %in% transects)

# split into sample info and asvs
asv <- edna16s_reads %>%
  select(starts_with('asv'))
sampinfo <- edna16s_reads %>%
  select(-starts_with('asv'))


# asv counts by sample
sampinfo %>%
  bind_cols(num.asvs = specnumber(asv))

# alpha diversity measures
alpha_div <- sampinfo %>%
  bind_cols(alpha.div.sh = diversity(asv, index = 'shannon'),
            alpha.div.si = diversity(asv, index = 'simpson'),
            alpha.div.fi = fisher.alpha(asv))

# diversity by depth
alpha_div %>%
  mutate(depth = as.numeric(depthm),
         depth.fac = cut_number(depth, 3)) %>%
  ggplot(aes(x = alpha.div.sh, y = depth.fac)) + 
  geom_boxplot()

# diversity distribution by cruise
alpha_div %>%
  ggplot(aes(y = alpha.div.sh, x = cruise)) +
  geom_boxplot(coef = 2, outlier.shape = NA) +
  scale_y_continuous(limits = c(3, 6)) +
  theme(axis.text.x = element_text(angle = 90))

# seasonal fluctuations?
alpha_div %>%
  group_by(cruise) %>%
  summarize(alpha = mean(alpha.div.sh)) %>%
  ggplot(aes(x = cruise, y = alpha, group = 1)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90))

# further explorations??



# pairwise beta diversity for one cruise
edna16s_reads %>%
  filter(cruise == cruises[1], as.numeric(depthm) < 30) %>%
  betadiver() %>%
  plot()

# read docs
?betadiver

# indices supported
betadiver(help = T)

# further explorations?