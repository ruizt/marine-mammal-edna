library(tidyverse)
library(modelr)

# make sure to source scripts/edna-aggregation-18s.R first, then:
load('data/edna-samples-imputed-18s.RData')

depth_weight_fn <- function(depth, alph, bet){
  dgamma(depth, shape = alph, scale = bet)
}

## STRATEGY 1

alpha_grid <- seq(1, 11, length = 3)
beta_grid <- seq(1, 15, length = 2)
alpha_divs <- matrix(rep(NA, 6), nrow = length(alpha_grid))

for(i in 1:3){
  for(j in 1:2){
    # aggregate
    edna_agg <- edna_imputed |>
      mutate(depthm = as.numeric(depthm),
             weight = depth_weight_fn(depthm, 
                                      alpha_grid[i],
                                      beta_grid[j])) |>
      group_by(cruise, line, sta) |>
      summarize(across(starts_with('asv'), 
                       ~weighted.mean(log(.x), weight)), 
                .groups = 'drop') |>
      group_by(cruise, line) |>
      summarize(across(starts_with('asv'), ~mean(.x))) |>
      group_by(cruise) |>
      summarize(across(starts_with('asv'), ~mean(.x)))
   
    # diversity
    alpha_divs[i, j] <- i + j
  }
}

## STRATEGY 2

par_grid <- expand_grid(alpha = seq(1, 11, length = 3),
                        beta = seq(1, 15, length = 2))

edna_agg_fn <- function(alph, bet){
  out <- edna_imputed |>
  mutate(depthm = as.numeric(depthm),
         weight = depth_weight_fn(depthm, 
                                  alph, bet)) |>
  group_by(cruise, line, sta) |>
  summarize(across(starts_with('asv'), 
                   ~weighted.mean(log(.x), weight)), 
            .groups = 'drop') |>
  group_by(cruise, line) |>
  summarize(across(starts_with('asv'), ~mean(.x))) |>
  group_by(cruise) |>
  summarize(across(starts_with('asv'), ~mean(.x)))
}

par_grid |>
  mutate(data = map2(alpha, beta, edna_agg_fn))
