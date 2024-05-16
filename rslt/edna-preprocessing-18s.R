library(tidyverse)

# read in metadata
metadata <- read_csv("data/NCOG_sample_log_DNA_meta_2014-2020.csv")

# read in 18s reads
edna_in <- read_tsv('data/NCOG_18sV9_asv_count_tax_S.tsv') %>%
  mutate(short.id = paste('asv', row_number(), sep = '.'))

# retain taxon names
taxa <- edna_in %>% 
  dplyr::select(where(is.character), silva_Confidence, pr2_Confidence)

# transpose and clean up names
edna_raw <- edna_in %>%
  dplyr::select(-starts_with('silva'), -starts_with('pr2'), -Feature.ID) %>%
  pivot_longer(-short.id, names_to = 'sample.id', values_to = 'read.count') %>%
  pivot_wider(names_from = short.id, values_from = read.count) %>% 
  separate(sample.id, 
           into = c('cruise', 'line', 'sta', 'depthm'), 
           sep = '_',
           remove = F)

# inspect samples without line or station assignment (verify)
edna_raw %>%
  filter(if_any(c(cruise, line, sta, depthm), is.na)) %>%
  pull(sample.id)

# extract genuine samples
edna_samples <- edna_raw %>%
  filter(if_all(c(cruise, line, sta, depthm), ~!is.na(.x))) 

# output
save(list = c('metadata', 'taxa', 'edna_samples'), file = 'data/ncog-18s.RData')
