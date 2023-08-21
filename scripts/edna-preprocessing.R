library(tidyverse)

# read in metadata
metadata_raw <- read_csv("data/NCOG_sample_log_DNA_meta_2014-2020.csv")

# 27 cruises total
metadata_raw %>%
  count(Cruise)

# station counts by transect and cruise; 15 total transects
metadata_raw %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') %>%
  count(Cruise, line) %>%
  arrange(line, Cruise) %>%
  pivot_wider(names_from = line, values_from = n) %>% 
  print(n = 27)

# common transects sampled across (most) cruises: 
sampling_counts <- metadata_raw %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') %>%
  count(Cruise, line) %>%
  mutate(line = as.numeric(line)) %>%
  filter((76 < line) & (line < 93.4)) 

transects <- sampling_counts %>% pull(line) %>% unique()

sampling_counts %>%
  arrange(line, Cruise) %>%
  pivot_wider(names_from = line, values_from = n) %>% 
  print(n = 27)


# read in 16s reads
edna_in <- read_tsv('data/NCOG_16S_asv_count_tax_S.tsv') %>%
  mutate(short.id = paste('asv', row_number(), sep = '.'))
taxa <- edna_in %>% 
  dplyr::select(where(is.character), silva_Confidence)

# transpose and clean up names
edna_raw <- edna_in %>%
  dplyr::select(-starts_with('silva'), -Feature.ID) %>%
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
imputation_out <- edna_samples %>% 
  dplyr::select(all_of(cols_of_interest)) %>%
  zCompositions::cmultRepl(label = 0, 
            method = 'GBM', 
            output = 'prop',
            z.warning = 1)

edna_imputed <- edna_samples %>%
  dplyr::select(1:5) %>%
  bind_cols(imputation_out)

# aggregate (verbose)
counts_depth <- edna_imputed %>%
  group_by(cruise, line, sta) %>%
  count(name = "n_depth")

counts_cruise <- counts_depth %>%
  group_by(cruise) %>%
  count(name = "n_sta")

edna_agg <- edna_imputed %>%
  left_join(counts_depth) %>%
  left_join(counts_cruise) %>%
  mutate(obs.weight = (1/n_depth)*(1/n_sta)) %>%
  # dplyr::select(cruise, line, sta, depthm, obs.weight)
  group_by(cruise) %>%
  summarize(across(starts_with('asv'), ~exp(weighted.mean(log(.x), obs.weight))))

# # aggregate (efficient)
# edna_agg <- edna_imputed %>%
#   group_by(cruise, line, sta) %>%
#   summarize(across(starts_with('asv'), ~mean(log(.x))), .groups = 'drop') %>%
#   group_by(cruise) %>%
#   summarize(across(starts_with('asv'), ~exp(mean(.x))))

# interpretation: average proportion of asv.XXX across survey is ZZZ
# problem: average proportions don't follow sum constraint exactly... rescale???
edna_agg %>% dplyr::select(-cruise) %>% rowSums()

# probably yes, note effect of imputation on aggregated compositions
edna_agg %>% 
  dplyr::select(-cruise) %>% 
  summarize(across(everything(), max)) %>%
  gather() %>%
  arrange(desc(value))

# *much* smaller maxima
edna_imputed %>% 
  dplyr::select(starts_with('asv')) %>% 
  summarize(across(everything(), max)) %>%
  gather() %>%
  arrange(desc(value))

# proposal (think/discuss/analyze)
edna_data <- select(edna_agg, cruise) %>%
  bind_cols(select(edna_agg, -cruise)/rowSums(select(edna_agg, -cruise)))

edna_data %>% select(-cruise) %>% rowSums()
edna_data %>% select(-cruise) %>% 
  summarize(across(everything(), max)) %>%
  gather() %>%
  arrange(desc(value))
