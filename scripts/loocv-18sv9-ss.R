library(tidyverse)
library(magrittr)
library(modelr)
library(spls)

## DATA INPUTS -----------------------------------------------------------------

# load edna and sighting data
load('data/processed/ncog18sv9.RData')
load('data/processed/mm-sightings-2024-07-27.RData')

# combine seasonally adjusted scaled sightings and seasonally adjusted edna data
loo_partitions <- left_join(loo_sightings, 
          loo_edna, 
          join_by(test.cruise), 
          suffix = c('.mm', '.edna')) |>
  select(-starts_with('.id')) |>
  mutate(train = map2(train.mm, train.edna, ~inner_join(.x, .y, join_by('cruise'))),
         test = map2(test.mm, test.edna, ~inner_join(.x, .y, join_by('cruise')))) |>
  select(test.cruise, test.season, train, test)

loo_partitions
  
## HYPERPARAMETER SPECIFICATION ------------------------------------------------

# hyperparameter grids
eta_grid_res <- 50
eta_grid <- rev(1 - exp(seq(log(0.1), log(0.55), length = eta_grid_res)))
ncomp_grid <- 4:8
obs_grid <- 1:25
species_grid <- c("bm", "bp", "mn")

# combine eta and observation grids for iteration
settings <- expand_grid(eta = eta_grid,
                        obs = obs_grid) |>
  mutate(setting = row_number()) |>
  select(setting, eta, obs)

## LEAVE ONE OUT CROSS VALIDATION ----------------------------------------------

# directory to store files
dir <- 'rslt/loocv/18sv9-ss'
fs::dir_create(dir)

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
              # model = fit,
              sel.asv = sel.asv,
              metrics = c(n.asv = df,
                          adj.rsq = rsq,
                          pred = y.test.hat[,1],
                          pred.err = test.resid[,1],
                          sq.pred.err = test.spe))
  return(out)
}

# execute
for(.species in species_grid){
  path <- paste(dir, '/_', .species, sep = '')
  fs::dir_create(path)
  lapply(ncomp_grid, function(.ncomp){
    rslt <- lapply(settings$setting, function(i){
      setting <- slice(settings, i)
      obs <- setting$obs
      eta <- setting$eta
      train <- loo_partitions$train[obs][[1]] %>% as.data.frame() 
      test <- loo_partitions$test[obs][[1]] %>% as.data.frame() 
      
      out <- loocv_fn(train, 
                      test, 
                      .var = .species, 
                      .parms = c(ncomp = .ncomp, eta = eta))
      paste('eta = ', round(eta, 4), ', obs = ', obs, sep = '') |> print()
      return(out)
    })
    
    sapply(settings$setting, function(i){
      c(obs = i, rslt[[i]]$parms, rslt[[i]]$metrics)
    }) |>
      t() |> 
      as_tibble() |>
      write_rds(file = paste(path, '/', .ncomp, 'comp-metrics.rds', sep = ''))
    
    lapply(settings$setting, function(i){rslt[[i]]$sel.asv}) |>
      write_rds(file = paste(path, '/', .ncomp, 'comp-selected.rds', sep = ''))
    
    # lapply(settings$setting, function(i){rslt[[i]]$model}) |>
    #   write_rds(file = paste(path, '/', .ncomp, 'comp-models.rds', sep = ''))
  })
}

## LOOCV SUMMARIES -------------------------------------------------------------

# read in and combine loocv run metrics
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

write_rds(loo_sel_freq, 
          file = paste(dir, 'selection-frequencies.rds', sep = '/'))

# store partitions
write_rds(loo_partitions, 
          file = paste(dir, 'partitions.rds', sep = '/'))
