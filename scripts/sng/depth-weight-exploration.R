library(tidyverse)
library(compositions)
library(zCompositions)
library(vegan)
load('data/ncog-18s.RData')

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

# aggregate to cruise level
depth_weight_fn <- function(depth, alpha, beta){
  dgamma(depth, shape = alpha, scale = beta) ##CHANGE????
}

##Grid search over shape and scale to maximize alpha diversity
#full list of summary statistics (Mean, median, sd, quartiles, min, max, etc)
alpha.vec <- seq(2,11,1)
beta.vec <- seq(1,15,1)

avgDivsDF <- data.frame(0,0,0)
names(avgDivsDF) = c("alpha", "beta", "avgDiv")


for (a in alpha.vec){
  for (b in beta.vec){
edna_agg <- edna_imputed |>
  mutate(depthm = as.numeric(depthm),
         weight = depth_weight_fn(depthm, a, b)) |>
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
  summarise(avgDiv = mean(alpha.div.sh))

newRow = c(a,b, alpha_div$avgDiv)

print(newRow)

avgDivsDF <- rbind(avgDivsDF, setNames(newRow,names(avgDivsDF)))
  }
}


avgDivsDF <- avgDivsDF |> 
  slice(-1)

avgDivsDF |> 
  filter((alpha-1)* beta >= 10)

library(rgl)

with(avgDivsDF, plot3d(x = alpha, y = beta, z = avgDiv))
