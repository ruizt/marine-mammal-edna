library(tidyverse)
library(compositions)
library(zCompositions)
load('data/ncog-18s.RData')

## FILTERING AND IMPUTATION

# dplyr::select columns of interest based on nonzero frequency across samples
cols_of_interest <- edna_samples |>
  dplyr::select(starts_with('asv')) |>
  as.matrix() |>
  apply(2, function(.x){mean(.x > 0)}) |>
  t() |>
  as_tibble() |>
  gather(col, prop.nz) |>
  filter(prop.nz > 0.3) |> # adjustable?
  pull(col)

# filter to samples with at least 25% nonzero reads
rows_of_interest <- edna_samples |> 
  dplyr::select(all_of(cols_of_interest)) |>
  mutate(across(where(is.numeric), ~ .x > 0)) |>
  mutate(total = rowMeans(pick(starts_with('asv'))),
         rownum = row_number()) |>
  dplyr::select(total, rownum) |>
  filter(total > 0.25) |>
  pull(rownum)

# zero imputation (bayesian multiplicative, martin-fernandez 2015)
imputation_out <- edna_samples |> 
  dplyr::select(all_of(cols_of_interest)) |>
  slice(rows_of_interest) |>
  zCompositions::cmultRepl(label = 0, 
                           method = 'GBM', 
                           output = 'prop',
                           z.warning = 0.75)

# bind imputed values to sample info
edna_imputed <- edna_samples |> 
  slice(rows_of_interest) |>
  dplyr::select(1:5) |>
  bind_cols(imputation_out)

save(edna_imputed, file = 'data/edna-samples-imputed-18s.RData')
load('data/edna-samples-imputed-18s.RData')

## AGGREGATION

edna_imputed |>
  mutate(depth = as.numeric(depthm)) |> 
  ggplot(aes(x = depth)) +
  geom_histogram(bins = 20) +
  facet_wrap(~cruise)

edna_imputed |>
  mutate(depth.fac = cut(as.numeric(depthm), 
                         breaks = c(0, 20, 50, 100, 200),
                         right = F)) |>
  group_by(cruise, line, sta, depth.fac) |>
  count() |>
  arrange(desc(n)) |>
  group_by(cruise, depth.fac) |>
  count() |>
  ggplot(aes(x = depth.fac, y = n)) +
  geom_col() +
  facet_wrap(~cruise)

edna_imputed |>
  group_by(cruise, line, sta) |>
  count() |>
  pull(n) |>
  unique()

# aggregate to cruise level
depth_weight_fn <- function(depth){
  dgamma(depth, shape = 2, scale = 5)
}

edna_agg <- edna_imputed |>
  mutate(depthm = as.numeric(depthm),
         weight = depth_weight_fn(depthm)) |>
  group_by(cruise, line, sta) |>
  summarize(across(starts_with('asv'), 
                   ~weighted.mean(log(.x), weight)), 
            .groups = 'drop') |>
  group_by(cruise, line) |>
  summarize(across(starts_with('asv'), ~mean(.x))) |>
  group_by(cruise) |>
  summarize(across(starts_with('asv'), ~mean(.x)))

# interpretation: average proportion of asv.XXX across cruise is ZZZ
edna_data <- edna_agg |> mutate(exp(pick(-cruise)))

# save
save(list = 'edna_data', file = 'data/edna_18s_processed.RData')
