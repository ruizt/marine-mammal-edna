library(tidyverse)
library(modelr)
library(magrittr)
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
  select(where(is.character), silva_Confidence, pr2_Confidence) |>
  rename_with(tolower) |>
  rename_with(~str_replace_all(.x, '_', '.')) |>
  separate(silva.taxon, 
           into = c('d', 'p', 'c', 'o', 'f', 'g'), sep = ';') |>
  select(feature.id, short.id, silva.confidence, d, p, c, o, f, g) |>
  mutate(across(everything(), ~str_replace(.x, '._', '') |> str_remove_all('_') |> str_trim())) 

# retrieve sample ids of samples to exclude
exclude_samples <- paste(in_dir, 'exclude_samples.txt', sep = '') |>
  read_delim(delim = '/n', skip = 3, col_names = 'sample.id') |>
  pull(sample.id)

# rearrange read counts and cross reference sample ids with metadata 
## (data dimensions: 1156 x 50409;  1156 obs, 50408 ASVs)
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
## (data dimensions: 1 x 8875; 8875 ASVs)
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
## (data dimensions: 1 x 944;  944 obs)
samples_of_interest <- edna_raw |> 
  select(sample.id, all_of(cols_of_interest)) |>
  filter(sample.id %in% pull(surface_dmc_samples, sample.id)) |>
  mutate(across(where(is.numeric), ~ .x > 0)) |>
  mutate(total = rowMeans(pick(starts_with('asv')))) |>
  select(sample.id, total) |>
  filter(total >= 0.01) |>
  pull(sample.id)


# implement filtering criteria above
## (data dimensions: 944 x 8876; 944 obs, 8875 ASVs)
edna_filtered <- edna_raw |>
  select(sample.id, all_of(cols_of_interest)) |>
  filter(sample.id %in% samples_of_interest)

# remove asvs present in under 1% of samples after filtering
edna_filtered %<>% select(where(~mean(.x > 0) > 0.01))

## IMPUTATION ------------------------------------------------------------------

# impute zeros ( x_{ijkl} )
imputation_out <- edna_filtered |> 
  select(-sample.id) |>
  zCompositions::cmultRepl(label = 0, 
                           method = 'GBM', 
                           output = 'prop',
                           z.warning = 0.99)

# bind imputed values to sample info
## (data dimensions: 944 x 8876; 944 obs, 8875 ASVs)
edna_imputed <- edna_filtered |>
  select(sample.id) |>
  bind_cols(imputation_out) |>
  left_join(surface_dmc_samples, by = 'sample.id') |>
  select(sample.id, depth.fac, starts_with('asv'))

# save imputed data
fs::dir_create(paste(out_dir, '_intermediates/', sep = ''))
save(edna_imputed, 
     file = paste(out_dir, '_intermediates/ncog18sv9-imputed-', 
                  today(), '.RData', sep = ''))

## AGGREGATION -----------------------------------------------------------------

load(paste(out_dir, '_intermediates/ncog18sv9-imputed-2024-07-25.RData', sep = ''))

# grid search to identify depth averaging weights optimizing alpha diversity
library(vegan)
# weight function
binned_weight_fn <- function(depth.range, weight1){
  return(ifelse(depth.range == "surface", w1, 1 - w1))
}

# grid search
gs_results <- data.frame(0,0,0,0,0,0,0,0,0)
names(gs_results) = c("w1", "w2", "avgDiv", "SD", "Min", "25", "50","75", "Max")

w1.vec = seq(0.05,0.99, by = 0.01)

for (w1 in w1.vec) {
  
  w2 = (1 - w1)
  
  edna_agg <- edna_imputed |>
    separate(sample.id, 
             into = c('cruise', 'line', 'sta', 'depth'),
             sep = '_') |> 
    drop_na() |> 
    mutate(weight = binned_weight_fn(depth.fac, w1)) |>
    group_by(cruise, line, sta) |>
    summarize(across(starts_with('asv'), 
                     ~weighted.mean(log(.x), weight)), ## automatically normalizes weights
              .groups = 'drop') |>
    group_by(cruise, line) |>
    summarize(across(starts_with('asv'), ~mean(.x))) |>
    group_by(cruise) |>
    summarize(across(starts_with('asv'), ~mean(.x)))
  
  edna_data <- edna_agg |> mutate(exp(pick(-cruise))) 
  
  # split into sample info and asvs
  asv <- edna_data |> 
    dplyr::select(starts_with('asv'))
  
  sampinfo <- edna_data |> 
    dplyr::select(-starts_with('asv'))
  
  # alpha diversity measures
  alpha_divs <- sampinfo %>%
    bind_cols(alpha.div.sh = diversity(asv, index = 'shannon')) 
  
  alpha_div <- alpha_divs |> 
    slice(-14) |> 
    summarise(avgDiv = mean(alpha.div.sh),
              sd = sd(alpha.div.sh),
              min = min(alpha.div.sh),,
              twentyfifth = quantile(alpha.div.sh, .25),
              fifty = quantile(alpha.div.sh, .50),
              seventyfifth = quantile(alpha.div.sh, .75),
              max = max(alpha.div.sh))
  
  newRow = c(w1,w2, alpha_div$avgDiv, alpha_div$sd, alpha_div$min, alpha_div$twentyfifth, alpha_div$fifty, alpha_div$seventyfifth, alpha_div$max)
  
  print(newRow)
  
  gs_results <- rbind(gs_results, setNames(newRow,names(gs_results)))
  
}


gs_results <- gs_results |> 
  slice(-1)

gs_results |> 
  filter(avgDiv == max(gs_results$avgDiv))

# GRID SEARCH RESULTS:
# Maximum Average Alpha Diversity at surface weight of 0.45 and max Chlorophyll weight of 0.55
# different from before

gs_results |> 
  ggplot(aes(y = avgDiv, x = w1)) +
  geom_line() + 
  labs(x= "Weight of Surface Measurements", y = "Average Alpha Diversity") + 
  geom_abline(intercept = 6.235241, slope = 0)

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
loo_edna <- crossv_loo(edna_clr) |>
  rename(train.raw = train, 
         test.raw = test) |>
  mutate(train = map(train.raw, adj_fn),
         test = map2(train.raw, test.raw, adj_fn),
         test.cruise = map(test, ~pull(.x, cruise))) |>
  select(.id, test.cruise, train, test) |>
  unnest(c(test.cruise))

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
              'edna',
              'loo_edna'), 
     file = paste(out_dir, 'ncog18sv9-', today(), '.RData', sep = ''))
