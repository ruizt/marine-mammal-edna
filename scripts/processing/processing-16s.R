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

# read in 16s reads   # Rows: 37903 Columns: 1544 
edna_in <- paste(in_dir, 'NCOG_16S_asv_count_tax_S.tsv', sep = '') |>
  read_tsv() |>
  mutate(short.id = paste('asv', row_number(), sep = '.'))

# retain taxon names      
taxa <- edna_in |> 
  select(where(is.character), silva_Confidence) |>
  rename_with(tolower) |>
  rename_with(~str_replace_all(.x, '_', '.')) |>
  separate(silva.taxon, 
           into = c('d', 'p', 'c', 'o', 'f', 'g'), sep = ';') |>
  select(feature.id, short.id, silva.confidence, d, p, c, o, f, g) |>
  mutate(across(everything(), 
                ~str_replace(.x, '._', '') |> str_remove_all('_') |> str_trim())) 

# retrieve sample ids of samples to exclude
## (drop up to 16 samples)
exclude_samples <- paste(in_dir, 'exclude_samples.txt', sep = '') |>
  read_delim(delim = '/n', skip = 3, col_names = 'sample.id') |>
  pull(sample.id)

# rearrange read counts and cross reference sample ids with metadata 
## (data dimensions: 1155 x 37904;  1155 obs, 37903 ASVs)
edna_raw <- edna_in |>
  select(-starts_with('silva'), -Feature.ID) |>
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
## (data dimensions: 918 x 2; 918 samples)
## Warning message: Expected 5 pieces. Missing pieces filled with `NA` in 1145 rows
surface_dmc_samples <- metadata |> 
  select(sample.name, chlora) |>
  separate(sample.name, 
           into = c('cruise', 'line', 'sta', 'depth', 'misc'),
           sep = '_') |>
  filter(is.na(misc)) |>
  select(-misc) |>
  filter(as.numeric(line) > 74) |>
  group_by(cruise, line, sta) |>
  mutate(max.chla = chlora == max(chlora),
         min.depth = as.numeric(depth) == min(as.numeric(depth))) |>
  filter(max.chla + min.depth == 1) |>
  mutate(depth.fac = factor(max.chla, labels = c('max.chla', 'surface'))) |>
  select(-where(is.logical), -chlora) |>
  unite(sample.id, c(cruise, line, sta, depth), sep = '_')

# select asv's that appear in at least 1% of samples and at most 99% of samples
## (data dimensions: 1 x 5951; 5951 ASVs)
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
## (data dimensions: 1 x 913;  913 obs)
samples_of_interest <- edna_raw |> 
  select(sample.id, all_of(cols_of_interest)) |>
  filter(sample.id %in% pull(surface_dmc_samples, sample.id)) |>
  mutate(across(where(is.numeric), ~ .x > 0)) |>
  mutate(total = rowMeans(pick(starts_with('asv')))) |>
  select(sample.id, total) |>
  filter(total >= 0.01) |>
  pull(sample.id)


# implement filtering criteria above
## (data dimensions: 913 x 5952; 949 obs, 5951 ASVs)
edna_filtered <- edna_raw |>
  select(sample.id, all_of(cols_of_interest)) |>
  filter(sample.id %in% samples_of_interest)

# identify stations with only one observation (i.e. sample at one depth only)
## (removing 3 samples)
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
## (data dimensions: 910 x 5952; 910 obs, 5951 ASVs)
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
## LONG RUNTIME: ~5-10min
imputation_out <- edna_filtered |> 
  select(-sample.id) |>
  zCompositions::cmultRepl(label = 0, 
                           method = 'GBM', 
                           output = 'prop',
                           z.warning = 1.001*z.max)

# bind imputed values to sample info
edna_imputed <- edna_filtered |>
  select(sample.id) |>
  bind_cols(imputation_out) |>
  left_join(surface_dmc_samples, by = 'sample.id') |>
  select(sample.id, depth.fac, starts_with('asv'))

# save imputed data
fs::dir_create(paste(out_dir, '_imputed/', sep = ''))
save(edna_imputed, 
     file = paste(out_dir, '_imputed/ncog16s-imputed.RData', sep = ''))

## AGGREGATION -----------------------------------------------------------------

# read in imputed data
load(paste(out_dir, '_imputed/ncog16s-imputed.RData', sep = ''))

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
# max average alpha diversity at surface weight 0.79max (max Chlorophyll weight ?)
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

# filter metadata to samples of interest
sample_metadata <- metadata %>% 
  filter(sample.name %in% samples_of_interest) |>
  rename(sample.id = sample.name)
attr(sample_metadata, 'spec') <- NULL
attr(sample_metadata, 'problems') <- NULL

# filter asv taxonomy to asvs of interest
asv_taxa <- taxa %>% filter(short.id %in% colnames(edna))

# export processed data
save(list = c('sample_metadata', 
              'asv_taxa',
              'edna'), 
     file = paste(out_dir, 'ncog16s.RData', sep = ''))

## LEAVE ONE OUT PARTITIONS ----------------------------------------------------

# generate leave one out partitions
loo_raw <- crossv_loo(edna_clr) |>
  rename(train.raw = train, 
         test.raw = test) 

# re-process seasonal adjustments for leave one out partitions
adj_fn <- function(.train, .test = NULL){
  
  if(is.null(.test)){ # for processing training data
    
    # subtract seasonal averages columnwise (on log scale)
    out <- .train |>
      as.data.frame() |>
      fmutate(qtr = ym(cruise) |> quarter()) |>
      fmutate(qtr = if_else(cruise == '201806', 3, qtr)) |>
      fgroup_by(qtr) |>
      fmutate(across(-c(cruise, qtr), function(.x){.x - mean(.x)})) |>
      ungroup() |>
      fselect(-qtr)
    
  }else{ # for processing test data
    
    # get seasonal averages from training data
    avgs <- .train |>
      as.data.frame() |>
      fmutate(qtr = ym(cruise) |> quarter()) |>
      fmutate(qtr = if_else(cruise == '201806', 3, qtr)) |>
      fgroup_by(qtr) |>
      fsummarize(across(-cruise, fmean)) |>
      pivot_longer(-qtr, names_to = 'asv', values_to = 'seasonal.mean')
    
    # subtract training seasonal averages columnwise from test data (on log scale)
    out <- .test |>
      as.data.frame() |>
      fmutate(qtr = ym(cruise) |> quarter()) |>
      fmutate(qtr = if_else(cruise == '201806', 3, qtr)) |>
      pivot_longer(-c(qtr, cruise), names_to = 'asv', values_to = 'value') |>
      left_join(avgs, by = c('qtr', 'asv')) |>
      fmutate(adj.value = value - seasonal.mean) |>
      pivot_wider(id_cols = cruise, values_from = adj.value, names_from = asv)
  }
  
  return(out)
}

# adjust for seasonal averages by partition
loo_edna <- lapply(1:nrow(loo_raw), function(j){
  .partition <- slice(loo_raw, j)
  .train <- adj_fn(.partition$train.raw[[1]])
  .test <- adj_fn(.partition$train.raw[[1]], .partition$test.raw[[1]])
  
  out <- tibble(test.id = .test$cruise,
                test = .test |> list(), 
                train = .train |> list())
  print(j)
  return(out)
}) %>%
  Reduce(bind_rows, .)

# export validation partitions
cv_dir <- paste(out_dir, '_partitions/', sep = '')
fs::dir_create(cv_dir)
save(list = c('loo_edna'),
     file = paste(cv_dir, 'ncog16s-partitions.RData', sep = ''))

## NESTED LEAVE ONE OUT PARTITIONS ---------------------------------------------

# generate data partitions for nested leave one out cross validation
loo_nested_raw <- crossv_loo(edna_clr, id = 'id.outer') |>
  rename(test.outer.raw = test) |>
  mutate(cv.inner = map(train, ~crossv_loo(as.data.frame(.x), id = 'id.inner'))) |>
  select(-train) |>
  unnest(cv.inner) |>
  rename(test.inner.raw = test,
         train.raw = train) |>
  select(test.outer.raw, test.inner.raw, train.raw)

# adjust for seasonal averages
loo_edna_nested <- lapply(1:nrow(loo_nested_raw), function(j){
  .partition <- slice(loo_nested_raw, j)
  .train <- adj_fn(.partition$train.raw[[1]])
  .test <- adj_fn(.partition$train.raw[[1]], .partition$test.inner.raw[[1]])
  .outer <- .partition$test.outer.raw[[1]] |> as_tibble() |> pull(cruise)
  .inner <- .partition$test.inner.raw[[1]] |> as_tibble() |> pull(cruise)
  
  out <- tibble(outer.id = .outer,
                inner.id = .inner,
                test = .test |> list(),
                train = .train |> list())
  print(j)
  return(out)
}) %>%
  Reduce(bind_rows, .)

# export processed data
save(list = c('loo_edna_nested'),
     file = paste(cv_dir, 'ncog16s-nested-partitions.RData', sep = ''))

