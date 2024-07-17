library(tidyverse)
library(magrittr)
library(modelr)
library(spls)

load('data/ncog-18s-processed-2024-07-05.RData')
load('data/ceta-density-processed-2024-07-05.RData')

# combine seasonally adjusted density estimates and seasonally adjusted edna data
whales <- inner_join(log_density_estimates_adj, edna_clr_adj, by = 'cruise')

## LEAVE ONE OUT CROSS VALIDATION ----------------------------------------------

# function to fit spls model and compute evaluation metrics
loocv_fn <- function(.train, .test, .var, .parms){
  # fit model
  x.train <- dplyr::select(.train, starts_with('asv'))
  y.train <- pull(.train, {{.var}})
  fit <- spls(x.train, y.train, 
              K = .parms['ncomp'], eta = .parms['eta'], 
              scale.x = F, scale.y = F)
  
  # summaries
  y.train.hat <- predict(fit, type = 'fit')
  train.resid <- y.train - y.train.hat
  n <- length(y.train)
  p <- ncol(fit$projection)
  df <- nrow(fit$projection)
  rsq <- as.numeric(1 - (n - p)*var(train.resid)/((n - 1)*var(y.train)))
  sel.asv <- fit$projection |> rownames()
  
  # squared prediction error
  x.test <- dplyr::select(.test, starts_with('asv'))
  y.test <- pull(.test, {{.var}})
  y.test.hat <- predict(fit, newx = x.test, type = 'fit')
  test.resid <- y.test - y.test.hat
  test.spe <- as.numeric(test.resid^2)
  
  # outputs
  out <- list(parms = .parms,
              species = deparse(substitute(.var)),
              model = fit,
              sel.asv = sel.asv,
              metrics = c(n.asv = df,
                          adj.rsq = rsq,
                          pred = y.test.hat[,1],
                          pred.err = test.resid[,1],
                          sq.pred.err = test.spe))
  return(out)
}

# loocv_fn(whales, whales[1, ], 'bm', c(eta = 0.5, ncomp = 3))

# hyperparameter grids
eta_grid_res <- 50
eta_grid <- rev(1 - seq(0.01, 0.95, length = eta_grid_res)^2)
ncomp_grid <- 1:13
obs_grid <- 1:25
species_grid <- c("bm", "bp", "mn")

# data partitions
data_partitions <- crossv_loo(whales)

# combine eta and observation grids for iteration
settings <- expand_grid(eta = eta_grid,
                        obs = obs_grid) |>
  mutate(setting = row_number()) |>
  select(setting, eta, obs)

dir <- paste('rslt/loocv', today(), sep = '/')
fs::dir_create(dir)

for(.species in species_grid){
  path <- paste(dir, '/_', .species, sep = '')
  fs::dir_create(path)
  lapply(ncomp_grid, function(.ncomp){
    rslt <- lapply(settings$setting, function(i){
      setting <- slice(settings, i)
      obs <- setting$obs
      eta <- setting$eta
      train <- data_partitions$train[obs][[1]] %>% as.data.frame() 
      test <- data_partitions$test[obs][[1]] %>% as.data.frame() 
      
      out <- loocv_fn(train, test, .var = .species, .parms = c(ncomp = .ncomp, eta = eta))
      paste('eta = ', round(eta, 4), ', obs = ', obs, sep = '') |> print()
      return(out)
    })
    
    sapply(settings$setting, function(i){c(obs = i, rslt[[i]]$parms, rslt[[i]]$metrics)}) |>
      t() |> as_tibble() |>
      write_rds(file = paste(path, '/', .ncomp, 'comp-metrics.rds', sep = ''))
    
    lapply(settings$setting, function(i){rslt[[i]]$sel.asv}) |>
      write_rds(file = paste(path, '/', .ncomp, 'comp-selected.rds', sep = ''))
    
    lapply(settings$setting, function(i){rslt[[i]]$model}) |>
      write_rds(file = paste(path, '/', .ncomp, 'comp-models.rds', sep = ''))
  })
}

# read in and combine metrics
loo_metrics <- lapply(species_grid, function(.species){
  path <- paste(dir, '/_', .species, sep = '')
  out <- lapply(ncomp_grid, function(.ncomp){
    read_rds(file = paste(path, '/', .ncomp, 'comp-metrics.rds', sep = ''))
  }) %>% 
    Reduce(bind_rows, .) %>%
    mutate(species = .species) %>%
    left_join(settings, by = c('eta', 'obs'))
  return(out)
}) %>%
  Reduce(bind_rows, .)

write_rds(loo_metrics, file = paste(dir, 'metrics.rds', sep = '/'))

# read in and combine selected asvs
loo_sel_asv <- lapply(species_grid, function(.species){
  path <- paste(dir, '/_', .species, sep = '')
  asv_list <- lapply(ncomp_grid, function(.ncomp){
    read_rds(file = paste(path, '/', .ncomp, 'comp-selected.rds', sep = ''))
  }) %>%
    Reduce(c, .)
  
  out <- expand_grid(ncomp = ncomp_grid, 
                     eta = eta_grid, 
                     obs = obs_grid) |>
    mutate(asv = asv_list,
           species = .species)
  
  return(out)
}) %>%
  Reduce(bind_rows, .) 

# tabulate selection frequencies
loo_sel_freq <- loo_sel_asv |>
  mutate(df = map(asv, length)) |>
  unnest(df) |>
  filter(df > 1) |>
  unnest(asv) |>
  # mutate(asv = str_remove(asv, '\\.')) |>
  group_by(species, ncomp, eta, asv) |>
  count() |>
  ungroup()

write_rds(loo_sel_freq, file = paste(dir, 'selection-frequencies.rds', sep = '/'))
  