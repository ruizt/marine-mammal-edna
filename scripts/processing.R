library(tidyverse)
library(compositions)
library(zCompositions)

## 18SV9 EDNA DATA -------------------------------------------------------------

## IMPORT READ COUNTS

# read in metadata
metadata <- read_csv("data/raw/NCOG_sample_log_DNA_meta_2014-2020.csv")

# read in 18s reads
edna_in <- read_tsv('data/raw/NCOG_18sV9_asv_count_tax_S.tsv') |>
  mutate(short.id = paste('asv', row_number(), sep = '.'))

# retain taxon names
taxa <- edna_in |> 
  dplyr::select(where(is.character), silva_Confidence, pr2_Confidence)

# transpose and clean up names 
## (data dimensions: 1536 x 50413;  1536 obs, 50408 ASVs)
edna_raw <- edna_in |>
  dplyr::select(-starts_with('silva'), -starts_with('pr2'), -Feature.ID) |>
  pivot_longer(-short.id, names_to = 'sample.id', values_to = 'read.count') |>
  pivot_wider(names_from = short.id, values_from = read.count) |> 
  separate(sample.id, 
           into = c('cruise', 'line', 'sta', 'depthm'), 
           sep = '_',
           remove = F)

# # inspect samples without line or station assignment
# edna_raw |>
#   filter(if_any(c(cruise, line, sta, depthm), is.na)) |>
#   pull(sample.id)

# extract genuine samples
## (data dimensions: 1486 x 50413;  1486 obs, 50408 ASVs)
edna_samples <- edna_raw |>
  filter(if_all(c(cruise, line, sta, depthm), ~!is.na(.x))) 

## FILTERING AND IMPUTATION

# select asv's that appear in at least 5% of samples and less than 90% of samples
## (data dimensions: 1 x 3248;  3248 ASVs)
cols_of_interest <- edna_samples |>
  dplyr::select(starts_with('asv')) |>
  as.matrix() |>
  apply(2, function(.x){mean(.x > 0)}) |>
  t() |>
  as_tibble() |>
  gather(col, prop.nz) |>
  filter(prop.nz > 0.05, prop.nz < 0.9) |> 
  pull(col)

# filter samples in which fewer than 10% of asv's of interest are present
## (data dimensions: 1 x 1185;  1185 obs)
rows_of_interest <- edna_samples |> 
  dplyr::select(all_of(cols_of_interest)) |>
  mutate(across(where(is.numeric), ~ .x > 0)) |>
  mutate(total = rowMeans(pick(starts_with('asv'))),
         rownum = row_number()) |>
  dplyr::select(total, rownum) |>
  filter(total > 0.1) |>
  pull(rownum)

# # check total number of zero reads (how much of the data will be imputed?)
# edna_samples |>
#   dplyr::select(all_of(cols_of_interest)) |>
#   slice(rows_of_interest) |>
#   gather() |>
#   summarize(nvals = n(),
#             num.nz = sum(value != 0),
#             prop.nz = mean(value != 0))

# # check asv presence/absence frequencies across samples after filtering
# edna_samples |> 
#   dplyr::select(all_of(cols_of_interest)) |>
#   slice(rows_of_interest) |>
#   as.matrix() |>
#   apply(2, function(.x){mean(.x > 0)}) |>
#   t() |>
#   as_tibble() |>
#   gather(col, prop.nz) |>
#   arrange(prop.nz) |>
#   summarize(n = sum(prop.nz > 0.01),
#             prop = n/n())

# impute zeros ( x_{ijkl} )
## (data dimensions: 1185 x 3248;  1185 obs, 3248 ASVs)
imputation_out <- edna_samples |> 
  dplyr::select(all_of(cols_of_interest)) |>
  slice(rows_of_interest) |>
  zCompositions::cmultRepl(label = 0, 
                           method = 'GBM', 
                           output = 'prop',
                           z.warning = 0.99)

# bind imputed values to sample info
## (data dimensions: 1185 x 3253;  1185 obs, 3248 ASVs)
edna_imputed <- edna_samples |> 
  slice(rows_of_interest) |>
  dplyr::select(1:5) |>
  bind_cols(imputation_out)

# binning surface vs deep measurements
min_depths <- edna_imputed |> 
  group_by(cruise, line, sta) |> 
  summarise(min.depth = min(as.numeric(depthm)), 
            md.binary = 1)

edna_imputed <- edna_imputed |> 
  mutate(depthm = as.numeric(depthm)) |> 
  left_join(min_depths, by = join_by(cruise,line,sta,depthm == min.depth))


edna_imputed <- edna_imputed |> 
  mutate(md.binary = replace_na(md.binary, 0))

edna_imputed <- edna_imputed |> 
  mutate(depth.range = case_when(md.binary == 1 ~ "Surface", 
                                 md.binary == 0 ~ "Deep"))

edna_imputed$depth.range <- factor(edna_imputed$depth.range, levels = c("Surface", "Deep"))

edna_imputed <- edna_imputed |> 
  mutate(depth.range = if_else(depth.range == "Surface" & depthm > 30, "Deep", depth.range))



# save imputed data, for grid search (remove later)
save(edna_imputed, file = paste('data/edna-imputed-', today(), '.RData', sep = ''))

## AGGREGATION

# function for weighting by depth
depth_weight_fn <- function(depth.range){
  return(ifelse(depth.range == "Surface", 0.09, 0.91))
}

# average first over depth, then over station, then over transect ( x_{ij} )
# interpretation: average proportion of asv.XX across cruise is ZZ
edna_aggregated <- edna_imputed |>
  mutate(depthm = as.numeric(depthm),
         weight = depth_weight_fn(depth.range)) |>
  group_by(cruise, line, sta) |>
  summarize(across(starts_with('asv'), 
                   ~weighted.mean(log(.x), weight)), 
            .groups = 'drop') |>
  group_by(cruise, line) |>
  summarize(across(starts_with('asv'), ~mean(.x))) |>
  group_by(cruise) |>
  summarize(across(starts_with('asv'), ~mean(.x))) |> 
  mutate(exp(pick(-cruise)))

## SEASONAL DE-TRENDING

# clr transformation ( log[z_{ij}] = log[x_{ij}/g_{i}] )
edna_clr <- edna_aggregated |> 
  dplyr::select(-cruise) |> 
  clr() |>
  as_tibble()

# center wrt seasonal (geometric) mean ( log[z_{ij}/g_{z_{j}}(i)] )
edna_clr_adj <- edna_aggregated |>
  dplyr::select(cruise) |>
  mutate(qtr = str_remove(cruise, 'X') |> ym() |> quarter()) |>
  mutate(qtr = if_else(cruise == 'X201806', 3, qtr)) |>
  bind_cols(edna_clr) |>
  group_by(qtr) |>
  mutate(across(starts_with('asv'), ~.x - mean(.x))) |>
  ungroup() |>
  dplyr::select(cruise, starts_with('asv'))

# export
save(list = c('metadata', 
              'taxa', 
              'edna_samples', 
              'edna_imputed', 
              'edna_aggregated', 
              'edna_clr_adj'), 
     file = paste('data/ncog-18s-processed-', today(), '.RData', sep = ''))

## SCALED SIGHTING DATA --------------------------------------------------------

# import whale sighting data ( y_{i} )
density_estimates <- read_csv("data/raw/CC_on_effort_scaled_sightings_Ceta_SCB.csv") |>
  rename_with(tolower) |>
  rename_with(~gsub('_.*', '', .x)) |>
  mutate(cruise = str_replace(cruiseid, "CC", "X20")) |>
  dplyr::select(cruise, season, bp, bm, mn)

# impute zeroes with uniform random numbers up to seasonal minima (no notation yet)
set.seed(30924)
density_estimates_imputed <- density_estimates |>
  group_by(season) |>
  summarize(across(where(is.numeric), 
                   .fns = list(min = ~min(.x[.x > 0], na.rm = T)),
                   .names = '{.col}.{.fn}')) |>
  right_join(density_estimates, by = 'season') |>
  mutate(bp.imp = if_else(bp == 0, 
                          runif(n = length(bp == 0), min = 0, max = bp.min),
                          bp),
         bm.imp = if_else(bm == 0,
                          runif(n = length(bm == 0), min = 0, max = bm.min),
                          bm),
         mn.imp = if_else(mn == 0,
                          runif(n = length(mn == 0), min = 0, max = mn.min),
                          mn))

# seasonal de-trending ( log[y_{i}/g_{y}(i)] )
log_density_estimates_adj <- density_estimates_imputed |>
  dplyr::select(cruise, season, ends_with('imp')) |>
  group_by(season) |>
  summarize(across(ends_with('imp'), 
                   .fns = list(mean = ~mean(log(.x), na.rm = T)),
                   .names = 'log.{.col}.{.fn}')) |>
  right_join(density_estimates_imputed) |>
  mutate(cruise = cruise,
         bp = log(bp.imp) - log.bp.imp.mean,
         bm = log(bm.imp) - log.bm.imp.mean,
         mn = log(mn.imp) - log.mn.imp.mean,
         .keep = "none")

# export
save(list = c('density_estimates', 'density_estimates_imputed', 'log_density_estimates_adj'), 
     file = paste('data/ceta-density-processed-', today(), '.RData', sep = ''))

