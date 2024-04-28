library(tidyverse)
library(compositions)
load('data/ncog-18s.RData')

# dplyr::select columns of interest based on nonzero frequency across samples
cols_of_interest <- edna_samples %>%
  dplyr::select(starts_with('asv')) %>%
  as.matrix() %>%
  apply(2, function(.x){mean(.x > 0)}) %>%
  t() %>%
  as_tibble() %>%
  gather(col, prop.nz) %>%
  filter(prop.nz > 0.5) %>% 
  pull(col)

edna_samples %>% 
  dplyr::select(all_of(cols_of_interest)) %>%
  mutate(across(everything(), ~(.x > 0))) %>%
  rowSums() %>%
  bind_cols(dplyr::select(edna_samples, 1:5)) %>%
  mutate(prop.asv.nz = `...1`/length(cols_of_interest)) %>%
  arrange(prop.asv.nz) %>%
  mutate(q = row_number()/nrow(edna_samples)) %>%
  ggplot(aes(x = prop.asv.nz, y = q)) +
  geom_step() +
  labs(y = 'proportion of samples', x = 'proportion of nonzero reads')

# zero imputation (bayesian multiplicative, martin-fernandez 2015, luz calle 2019)
imputation_out <- edna_samples %>% # use for relative abundances (maybe)
  dplyr::select(all_of(cols_of_interest)) %>%
  zCompositions::cmultRepl(label = 0, 
                           method = 'GBM', 
                           output = 'prop',
                           z.warning = 1)

edna_imputed <- edna_samples %>%
  dplyr::select(1:5) %>%
  bind_cols(imputation_out)

# aggregate (efficient)
edna_agg <- edna_imputed %>%
  group_by(cruise, line, sta) %>%
  summarize(across(starts_with('asv'), ~mean(log(.x))), .groups = 'drop') %>%
  group_by(cruise) %>%
  summarize(across(starts_with('asv'), ~mean(.x)))

# interpretation: average proportion of asv.XXX across cruise is ZZZ
edna_data <- edna_agg %>% mutate(exp(pick(-cruise)))

# compute closure using isometric logratio transform
cruises <- edna_data %>% select(cruise)

edna_data_clo <- edna_data %>% select(-cruise) %>% 
  ilr() %>% 
  ilrInv() %>% as.data.frame() 

edna_data_clo <- cbind(cruises, edna_data_clo)

# save 
save(list = 'edna_data', file = 'data/edna_18s_processed.RData')
