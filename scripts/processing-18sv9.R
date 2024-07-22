library(tidyverse)
library(compositions)
library(zCompositions)
library(fs)
out_dir <- 'data/processed/'
in_dir <- 'data/_raw/'
  
## IMPORT READ COUNTS ----------------------------------------------------------

# read in metadata
metadata <- paste(in_dir, "NCOG_sample_log_DNA_stvx_meta_2014-2020.csv", sep = '') |> 
  read_csv() |>
  rename_with(tolower) |>
  rename_with(~str_replace_all(.x, '_', '.'))

# read in 18s reads
edna_in <- paste(in_dir, 'NCOG_18sV9_asv_count_tax_S.tsv', sep = '') |>
  read_tsv() |>
  mutate(short.id = paste('asv', row_number(), sep = '.'))

# retain taxon names
taxa <- edna_in |> 
  dplyr::select(where(is.character), silva_Confidence, pr2_Confidence) |>
  rename_with(tolower) |>
  rename_with(~str_replace_all(.x, '_', '.')) |>
  separate(silva.taxon, 
           into = c('d', 'p', 'c', 'o', 'f', 'g'), sep = ';') |>
  dplyr::select(feature.id, short.id, silva.confidence, d, p, c, o, f, g) |>
  mutate(across(everything(), ~str_replace(.x, '._', '') |> str_remove_all('_') |> str_trim())) 

# retrieve sample ids of samples to exclude
exclude_samples <- paste(in_dir, 'exclude_samples.txt', sep = '') |>
  read_delim(delim = '/n', skip = 3, col_names = 'sample.id') |>
  pull(sample.id)

# rearrange read counts and cross reference sample ids with metadata 
## (data dimensions: 1144 x 50413;  1144 obs, 50408 ASVs)
edna_raw <- edna_in |>
  dplyr::select(-starts_with('silva'), -starts_with('pr2'), -Feature.ID) |>
  rename_with(~str_remove_all(.x, 'X')) |>
  pivot_longer(-short.id, names_to = 'sample.id', values_to = 'read.count') |>
  pivot_wider(names_from = short.id, values_from = read.count) |> 
  filter(sample.id %in% metadata$sample.name,
         !str_ends(sample.id, '_S'), 
         !str_ends(sample.id, '_T'),
         !(sample.id %in% exclude_samples))


## FILTERING -------------------------------------------------------------------

# select asv's that appear in at least 1% of samples and at most 99% of samples
## (data dimensions: 1 x 6832; 6832 ASVs)
cols_of_interest <- edna_raw |>
  dplyr::select(starts_with('asv')) |>
  as.matrix() |>
  apply(2, function(.x){mean(.x > 0)}) |>
  t() |>
  as_tibble() |>
  gather(col, prop.nz) |>
  filter(prop.nz >= 0.01, prop.nz <= 0.99) |> 
  pull(col)


# identify surface and dmc samples
## (data dimensions: 956 x 2; 956 samples)
surface_dmc_samples <- metadata |> 
  dplyr::select(sample.name, chlora) |>
  separate(sample.name, 
           into = c('cruise', 'line', 'sta', 'depth', 'misc'),
           sep = '_') |>
  filter(is.na(misc)) |>
  dplyr::select(-misc) |>
  group_by(cruise, line, sta) |>
  mutate(max.chla = chlora == max(chlora),
         min.depth = as.numeric(depth) == min(as.numeric(depth))) |>
  filter(max.chla + min.depth == 1) |>
  mutate(depth.fac = factor(max.chla, labels = c('max.chla', 'surface'))) |>
  dplyr::select(-where(is.logical), -chlora) |>
  unite(sample.id, c(cruise, line, sta, depth), sep = '_')

# identify samples in which at least 1% of asv's of interest are present
## (data dimensions: 1 x 948;  948 obs)
samples_of_interest <- edna_raw |> 
  dplyr::select(sample.id, all_of(cols_of_interest)) |>
  filter(sample.id %in% pull(surface_dmc_samples, sample.id)) |>
  mutate(across(where(is.numeric), ~ .x > 0)) |>
  mutate(total = rowMeans(pick(starts_with('asv')))) |>
  dplyr::select(sample.id, total) |>
  filter(total >= 0.01) |>
  pull(sample.id)


# implement filtering criteria above
## (data dimensions: 948 x 6832)
edna_filtered <- edna_raw |>
  dplyr::select(sample.id, all_of(cols_of_interest)) |>
  filter(sample.id %in% samples_of_interest)

# # check total number of zero reads 
# ## (91.5% of the data will be imputed)
# edna_filtered |>
#   gather() |>
#   summarize(nvals = n(),
#             num.nz = sum(value != 0),
#             prop.nz = mean(value != 0))
# 
# # check asv presence/absence frequencies across samples after filtering
# ## (all asvs are present in at least 0.1% of samples)
# edna_filtered |>
#   as.matrix() |>
#   apply(2, function(.x){mean(.x > 0)}) |>
#   t() |>
#   as_tibble() |>
#   gather(col, prop.nz) |>
#   arrange(prop.nz) |>
#   summarize(n = sum(prop.nz > 0.001),
#             prop = n/n())

## IMPUTATION ------------------------------------------------------------------

# impute zeros ( x_{ijkl} )
imputation_out <- edna_filtered |> 
  dplyr::select(-sample.id) |>
  zCompositions::cmultRepl(label = 0, 
                           method = 'GBM', 
                           output = 'prop',
                           z.warning = 0.999)

# bind imputed values to sample info
## (data dimensions: 948 x 6834;  948 obs, 6832 ASVs)
edna_imputed <- edna_filtered |>
  dplyr::select(sample.id) |>
  bind_cols(imputation_out) |>
  left_join(surface_dmc_samples, by = 'sample.id') |>
  dplyr::select(sample.id, depth.fac, starts_with('asv'))

# save imputed data
save(edna_imputed, 
     file = paste(out_dir, '_intermediates/ncog18sv9-imputed-', today(), '.RData', sep = ''))

## AGGREGATION -----------------------------------------------------------------

load(paste(out_dir, '_intermediates/ncog18sv9-imputed-2024-07-20.RData', sep = ''))

# grid search to identify depth averaging weights optimizing alpha diversity



# average first over depth, then over station, then over transect ( x_{ij} )
# interpretation: average proportion of asv.XX across cruise is ZZ
edna_aggregated <- edna_imputed |>
  mutate(depth.wt = ifelse(depth.fac == "surface", 0.09, 0.91)) |>
  separate(sample.id, into = c('cruise', 'line', 'sta', 'depth'), sep = '_') |>
  group_by(cruise, line, sta) |>
  summarize(across(starts_with('asv'), 
                   ~weighted.mean(log(.x), depth.wt)), 
            .groups = 'drop') |>
  group_by(cruise, line) |>
  summarize(across(starts_with('asv'), ~mean(.x))) |>
  group_by(cruise) |>
  summarize(across(starts_with('asv'), ~mean(.x))) |> 
  mutate(exp(pick(-cruise)))

## CLR AND SEASONAL DE-TRENDING ------------------------------------------------

# clr transformation ( log[z_{ij}] = log[x_{ij}/g_{i}] )
edna_clr <- edna_aggregated |> 
  dplyr::select(-cruise) |> 
  clr() |>
  as_tibble()

# center wrt seasonal (geometric) mean ( log[z_{ij}/g_{z_{j}}(i)] )
edna <- edna_aggregated |>
  dplyr::select(cruise) |>
  mutate(qtr = ym(cruise) |> quarter()) |>
  mutate(qtr = if_else(cruise == '201806', 3, qtr)) |>
  bind_cols(edna_clr) |>
  group_by(qtr) |>
  mutate(across(starts_with('asv'), ~.x - mean(.x))) |>
  ungroup() |>
  dplyr::select(cruise, starts_with('asv'))

## EXPORT ----------------------------------------------------------------------

# filter metadata to samples of interest
sample_metadata <- metadata %>% filter(sample.name %in% samples_of_interest) |>
  rename(sample.id = sample.name)
attr(sample_metadata, 'spec') <- NULL
attr(sample_metadata, 'problems') <- NULL

# filter asv taxonomy to asvs of interest
asv_taxa <- taxa %>% filter(short.id %in% colnames(edna))

save(list = c('sample_metadata', 
              'asv_taxa',
              'edna'), 
     file = paste(out_dir, 'ncog18sv9-', today(), '.RData', sep = ''))
