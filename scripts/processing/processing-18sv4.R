library(tidyverse)
library(modelr)
library(magrittr)
library(fs)
library(vegan)
library(collapse)
out_dir <- 'data/processed/'
in_dir <- 'data/_raw/'

## IMPORT READ COUNTS ----------------------------------------------------------

# read in metadata
metadata <- paste(in_dir, "NCOG_sample_log_DNA_stvx_meta_2014-2020.csv", sep = '') |> 
  read_csv() |>
  rename_with(tolower) |>
  rename_with(~str_replace_all(.x, '_', '.'))

# read in 18sV4 reads   # Rows: 42032 Columns: 1552
edna_in <- paste(in_dir, 'NCOG_18sV4_asv_count_tax.tsv', sep = '') |>
  read_tsv() |>
  mutate(short.id = paste('asv', row_number(), sep = '.'))

# retain taxon names   # Warning: pieces discarded and filled with NA
taxa <- edna_in |> 
  select(where(is.character), silva_Confidence, pr2_Confidence) |>
  rename_with(tolower) |>
  rename_with(~str_replace_all(.x, '_', '.')) |>
  separate(silva.taxon, 
           into = c('d', 'p', 'c', 'o', 'f', 'g'), sep = ';') |>
  select(feature.id, short.id, silva.confidence, d, p, c, o, f, g) |>
  mutate(across(everything(), 
                ~str_replace(.x, '._', '') |> str_remove_all('_') |> str_trim())) 

# retrieve sample ids of samples to exclude
## (drop 16 samples)
exclude_samples <- paste(in_dir, 'exclude_samples.txt', sep = '') |>
  read_delim(delim = '/n', skip = 3, col_names = 'sample.id') |>
  pull(sample.id)

# rearrange read counts and cross reference sample ids with metadata 
## (data dimensions: 1156 x 42032;  1156 obs, 42022 ASVs)
edna_raw <- edna_in |>
  select(-starts_with('silva'), -starts_with('pr2'), -Feature.ID) |>
  rename_with(~str_remove_all(.x, 'X')) |>
  pivot_longer(-short.id, names_to = 'sample.id', values_to = 'read.count') |>
  pivot_wider(names_from = short.id, values_from = read.count) |> 
  filter(sample.id %in% metadata$sample.name,
         !(sample.id %in% exclude_samples)) |>
  mutate(sterivex = str_ends(sample.id, 'S'),
         sample.id = str_remove_all(sample.id, '_S')) |>
  group_by(sample.id) |>
  slice_max(sterivex) |>
  select(-sterivex) |>
  ungroup()

## FILTERING -------------------------------------------------------------------

# identify surface and dmc samples
## (data dimensions: 956 x 2; 956 samples)
## Warning message: Expected 5 pieces. Missing pieces filled with `NA` in 1145 rows
surface_dmc_samples <- metadata |> 
  select(sample.name, chlora) |>
  separate(sample.name, 
           into = c('cruise', 'line', 'sta', 'depth', 'misc'),
           sep = '_') |>
  filter(is.na(misc)) |>
  select(-misc) |>
  group_by(cruise, line, sta) |>
  mutate(max.chla = chlora == max(chlora),
         min.depth = as.numeric(depth) == min(as.numeric(depth))) |>
  filter(max.chla + min.depth == 1) |>
  mutate(depth.fac = factor(max.chla, labels = c('max.chla', 'surface'))) |>
  select(-where(is.logical), -chlora) |>
  unite(sample.id, c(cruise, line, sta, depth), sep = '_')

# select asv's that appear in at least 1% of samples and at most 99% of samples
## (data dimensions: 1 x 6457; 6457 ASVs)
cols_of_interest <- edna_raw |>
  filter(sample.id %in% pull(surface_dmc_samples, sample.id)) |>
  select(starts_with('asv')) |>
  as.matrix() |>
  apply(2, function(.x){mean(.x > 0)}) |>
  t() |>
  as_tibble() |>
  gather(col, prop.nz) |>
  filter(prop.nz >= 0.01, prop.nz <= 0.99) |> 
  pull(col)

# identify samples in which at least 1% of asv's of interest are present
## (data dimensions: 1 x 925;  925 obs)
samples_of_interest <- edna_raw |> 
  select(sample.id, all_of(cols_of_interest)) |>
  filter(sample.id %in% pull(surface_dmc_samples, sample.id)) |>
  mutate(across(where(is.numeric), ~ .x > 0)) |>
  mutate(total = rowMeans(pick(starts_with('asv')))) |>
  select(sample.id, total) |>
  filter(total >= 0.01) |>
  pull(sample.id)


# implement filtering criteria above
## (data dimensions: 925 x 6458; 925 obs, 6457 ASVs)
edna_filtered <- edna_raw |>
  select(sample.id, all_of(cols_of_interest)) |>
  filter(sample.id %in% samples_of_interest)

# identify stations with only one observation (i.e. sample at one depth only)
## (removing 15 samples)
samples_to_drop <- edna_filtered |>
  separate(sample.id, 
           into = c('cruise', 'line', 'sta', 'depth'),
           sep = '_') |>
  group_by(cruise, line, sta) |>
  count() |>
  filter(n != 2) |>
  unite(sample.id, c(cruise, line, sta), sep = '_') |>
  pull(sample.id)

# remove single-depth samples and asvs present in under 1% of remaining samples
## (data dimensions: 910 x 6458; 910 obs, 6457 ASVs)
edna_filtered %<>% 
  filter(!(str_trunc(sample.id, 18, ellipsis = '') %in% samples_to_drop)) |>
  select(where(~mean(.x > 0) > 0.01))

# determine limit for z.warning: check maximum column sparsity (prop. zeroes)
z.max <- edna_filtered |> 
  pivot_longer(starts_with('asv')) |>
  group_by(name) |>
  summarize(prop.z = mean(value == 0)) |>
  pull(prop.z) |>
  max()

## IMPUTATION ------------------------------------------------------------------

# impute zeros ( x_{ijkl} )
## LONG RUNTIME: ~5-10min, Number of adjusted imputations:  664269 
imputation_out <- edna_filtered |> 
  select(-sample.id) |>
  zCompositions::cmultRepl(label = 0, 
                           method = 'GBM', 
                           output = 'prop',
                           z.warning = 1.001*z.max)

# bind imputed values to sample info
## (data dimensions: 910 x 6459; 910 obs, 6457 ASVs + 1 ID + 1 depth indicator)
edna_imputed <- edna_filtered |>
  select(sample.id) |>
  bind_cols(imputation_out) |>
  left_join(surface_dmc_samples, by = 'sample.id') |>
  select(sample.id, depth.fac, starts_with('asv'))



# save imputed data
## REMOVE DIRECTORY FROM VERSION CONTROL BEFORE FINALIZING
fs::dir_create(paste(out_dir, '_imputed/', sep = ''))
save(edna_imputed, 
     file = paste(out_dir, '_imputed/ncog18sv4-imputed.RData', sep = ''))

## AGGREGATION -----------------------------------------------------------------

# read in imputed data
load(paste(out_dir, '_imputed/ncog18sv4-imputed.RData', sep = ''))

# weight function for depth aggregation
weight_fn <- function(depth.range, w){
  return(ifelse(depth.range == "surface", w, 1 - w))
}

# candidate weights
w.vec <- seq(0.05,0.99, by = 0.01)

# grid search to identify depth averaging weights optimizing alpha diversity
gs.rslt <- sapply(w.vec, function(w){
  
  # weighted aggregation to cruise level (first depth, then station, then line)
  # then compute average alpha diversity across cruises
  out <- edna_imputed |>
    separate(sample.id, 
             into = c('cruise', 'line', 'sta', 'depth'),
             sep = '_') |>
    group_by(cruise) |>
    fmutate(across(.cols = is.numeric, log)) |>
    mutate(weight = weight_fn(depth.fac, w)) |>
    fgroup_by(cruise, line, sta) |>
    fselect(-depth.fac, -depth) |>
    fmean(weight, keep.w = F) |>
    fgroup_by(cruise, line) |>
    fselect(-sta) |>
    fmean() |>
    fgroup_by(cruise) |>
    fselect(-line) |>
    fmean() |>
    fmutate(across(.cols = is.numeric, exp)) |>
    fselect(-cruise) |>
    diversity(index = 'shannon') |>
    mean()
  
  return(out)
})

# optimal weight
# max average alpha diversity at surface weight 0.48 (max Chlorophyll weight ?)
w.opt <- w.vec[which.max(gs.rslt)]

# average first over depth, then over station, then over transect ( x_{ij} )
# interpretation: average proportion of asv.XX across cruise is ZZ
edna_aggregated <- edna_imputed |>
  separate(sample.id, 
           into = c('cruise', 'line', 'sta', 'depth'),
           sep = '_') |>
  group_by(cruise) |>
  fmutate(across(.cols = is.numeric, log)) |>
  mutate(weight = weight_fn(depth.fac, w.opt)) |>
  fgroup_by(cruise, line, sta) |>
  fselect(-depth.fac, -depth) |>
  fmean(weight, keep.w = F) |>
  fgroup_by(cruise, line) |>
  fselect(-sta) |>
  fmean() |>
  fgroup_by(cruise) |>
  fselect(-line) |>
  fmean() |>
  fmutate(across(.cols = is.numeric, exp)) 

## CLR AND SEASONAL DE-TRENDING ------------------------------------------------

# clr transformation ( log[z_{ij}] = log[x_{ij}/g_{i}] )
edna_clr <- edna_aggregated |> 
  select(-cruise) |> 
  compositions::clr() |>
  as_tibble() |>
  bind_cols(select(edna_aggregated, cruise))

# center wrt seasonal (geometric) mean ( log[z_{ij}/g_{z_{j}}(i)] )
edna <- edna_clr |>
  mutate(qtr = ym(cruise) |> quarter()) |>
  mutate(qtr = if_else(cruise == '201806', 3, qtr)) |>
  group_by(qtr) |>
  mutate(across(starts_with('asv'), ~.x - mean(.x))) |>
  ungroup() |>
  select(cruise, starts_with('asv'))

# re-process seasonal adjustments for leave one out runs
# .train <- edna_clr |> slice(-1)
# .test <- edna_clr |> slice(1)
adj_fn <- function(.train, .test = NULL){
  if(is.null(.test)){
    out <- .train |>
      as.data.frame() |>
      mutate(qtr = ym(cruise) |> quarter()) |>
      mutate(qtr = if_else(cruise == '201806', 3, qtr)) |>
      group_by(qtr) |>
      mutate(across(starts_with('asv'), ~.x - mean(.x))) |>
      ungroup() |>
      select(cruise, starts_with('asv'))
  }else{
    avgs <- .train |>
      as.data.frame() |>
      mutate(qtr = ym(cruise) |> quarter()) |>
      mutate(qtr = if_else(cruise == '201806', 3, qtr)) |>
      group_by(qtr) |>
      summarize(across(starts_with('asv'), mean)) |>
      pivot_longer(-qtr, names_to = 'asv', values_to = 'seasonal.mean')
    
    out <- .test |>
      as.data.frame() |>
      mutate(qtr = ym(cruise) |> quarter()) |>
      mutate(qtr = if_else(cruise == '201806', 3, qtr)) |>
      pivot_longer(-c(qtr, cruise), names_to = 'asv', values_to = 'value') |>
      left_join(avgs, by = c('qtr', 'asv')) |>
      mutate(adj.value = value - seasonal.mean) |>
      pivot_wider(id_cols = cruise, values_from = adj.value, names_from = asv)
  }
  
  return(out)
}

# generate data partitions for leave one out cross validation
## LONG RUNTIME: ~5min
loo_edna <- crossv_loo(edna_clr) |>
  rename(train.raw = train, 
         test.raw = test) |>
  mutate(train = map(train.raw, adj_fn),
         test = map2(train.raw, test.raw, adj_fn),
         test.cruise = map(test, ~pull(.x, cruise))) |>
  select(.id, test.cruise, train, test) |>
  unnest(c(test.cruise))

# # generate data partitions for nested leave one out cross validation
# loo_edna_nested <- crossv_loo(edna_clr, id = 'id.outer') |>
#   rename(test.outer.raw = test) |>
#   mutate(cv.inner = map(train, ~crossv_loo(as.data.frame(.x), id = 'id.inner'))) |>
#   select(-train) |>
#   unnest(cv.inner) |>
#   rename(test.inner.raw = test,
#          train.raw = train) |>
#   mutate(train.inner = map(train.raw, adj_fn),
#          test.inner = map2(train.raw, test.inner.raw, adj_fn),
#          test.inner.cruise = map(test.inner, ~pull(.x, cruise)),
#          test.outer.cruise = map(test.outer.raw, ~as.data.frame(.x) |> pull(cruise))) |>
#   select(test.outer.cruise, test.inner.cruise, test.inner, train.inner) |>
#   unnest(ends_with('cruise'))


## EXPORT ----------------------------------------------------------------------

# filter metadata to samples of interest
sample_metadata <- metadata %>% filter(sample.name %in% samples_of_interest) |>
  rename(sample.id = sample.name)
attr(sample_metadata, 'spec') <- NULL
attr(sample_metadata, 'problems') <- NULL

# filter asv taxonomy to asvs of interest
asv_taxa <- taxa %>% filter(short.id %in% colnames(edna))

# export processed data
save(list = c('sample_metadata', 
              'asv_taxa',
              'edna'), 
     file = paste(out_dir, 'ncog18sv4.RData', sep = ''))

# export validation partitions
cv_dir <- paste(out_dir, '_cv/', sep = '')
fs::dir_create(cv_dir)
save(list = c('loo_edna'),
     file = paste(cv_dir, 'ncog18sv4-partitions.RData', sep = ''))
