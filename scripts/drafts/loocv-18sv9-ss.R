library(tidyverse)
library(fs)
library(modelr)
library(spls)
partition_dir <- 'data/processed/_combined-partitions'
out_dir <- 'rslt/loocv/18sv9-ss/'
dir_create(out_dir)

## HYPERPARAMETER SPECIFICATION ------------------------------------------------

# hyperparameter grids
eta_grid_res <- 50
eta_grid <- rev(1 - exp(seq(log(0.075), log(0.6), length = eta_grid_res)))
ncomp_grid <- 3:12

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
              sel.asv = sel.asv,
              metrics = c(n.asv = df,
                          adj.rsq = rsq,
                          pred = y.test.hat[,1],
                          pred.err = test.resid[,1],
                          sq.pred.err = test.spe))
  return(out)
}

## LEAVE ONE OUT CROSS VALIDATION ----------------------------------------------

# read in validation partitions
loo_partitions <- paste(partition_dir, 'partitions-18sv9-ss.rds') |> read_rds()

# iteration grids
partition_ix <- 1:25
species_grid <- c("bm", "bp", "mn")

# index eta and partition combos for innermost iteration
settings <- expand_grid(eta = eta_grid,
                        part_ix = partition_ix) |>
  mutate(setting = row_number()) |>
  select(setting, eta, part_ix)

# iterate over species
for(.species in species_grid){
  
  # define species-specific path
  path <- paste(out_dir, '/_', .species, sep = '')
  dir_create(path)
  
  # iterate over numbers of components
  lapply(ncomp_grid, function(.ncomp){
    
    # iterate over 'settings' (data partitions/sparsity parameter combo)
    rslt <- lapply(settings$setting, function(i){
      
      # identify setting
      setting <- slice(settings, i)
      partition <- setting$part_ix
      eta <- setting$eta
      
      # retrieve corresponding data partition
      train <- loo_partitions$train[partition][[1]] %>% as.data.frame() 
      test <- loo_partitions$test[partition][[1]] %>% as.data.frame() 
      
      # fit model and return fit metrics, prediction, and selected asvs
      out <- loocv_fn(train, 
                      test, 
                      .var = .species, 
                      .parms = c(ncomp = .ncomp, eta = eta))
      
      # print setting info in console
      paste('species = ', .species,
            ', ncomp = ', .ncomp,
            ', eta = ', round(eta, 4), 
            ', partition = ', partition, sep = '') |> print()
      return(out)
    })
    
    # write setting info and fit/prediction metrics to file
    sapply(settings$setting, function(i){
      c(setting = i, 
        rslt[[i]]$parms, 
        rslt[[i]]$metrics)
    }) |>
      t() |> 
      as_tibble() |>
      write_rds(file = paste(path, '/', .ncomp, 'comp-metrics.rds', sep = ''))
    
    # write selected ASVs to file
    lapply(settings$setting, function(i){rslt[[i]]$sel.asv}) |>
      write_rds(file = paste(path, '/', .ncomp, 'comp-selected.rds', sep = ''))
  })
}

## GENERATE SUMMARY OUTPUTS ----------------------------------------------------

# read in and combine loocv run metrics across species
loo_metrics <- lapply(species_grid, function(.species){
  
  # identify species-specific path
  path <- paste(out_dir, '/_', .species, sep = '')
  
  # read and combine fit/prediction metrics across numbers of components
  out <- lapply(ncomp_grid, function(.ncomp){
    read_rds(file = paste(path, '/', .ncomp, 'comp-metrics.rds', sep = ''))
  }) %>% 
    Reduce(bind_rows, .) %>%
    mutate(species = .species) %>%
    left_join(settings, join_by(setting, eta))
  return(out)
}) %>%
  Reduce(bind_rows, .)

# export to file
write_rds(loo_metrics, file = paste(out_dir, 'metrics.rds', sep = '/'))

# read in and combine selected asvs across species
loo_sel_asv <- lapply(species_grid, function(.species){
  
  # identify species-specific path
  path <- paste(out_dir, '/_', .species, sep = '')
  
  # read and combine selected asv lists across numbers of components
  asv_list <- lapply(ncomp_grid, function(.ncomp){
    read_rds(file = paste(path, '/', .ncomp, 'comp-selected.rds', sep = '')) 
  }) %>%
    Reduce(c, .)
  
  # format as tibble
  out <- expand_grid(ncomp = ncomp_grid, 
                     eta = eta_grid, 
                     part_ix = partition_ix) |>
    mutate(asv = asv_list,
           species = .species)
  
  return(out)
}) %>%
  Reduce(bind_rows, .) 

# export to file
write_rds(loo_sel_asv, file = paste(out_dir, 'sel-asv.rds', sep = '/'))

# tabulate selection frequencies
loo_sel_freq <- loo_sel_asv |>
  mutate(df = map(asv, length)) |>
  unnest(df) |>
  filter(df > 1) |>
  unnest(asv) |>
  group_by(species, ncomp, eta, asv) |>
  count() |>
  ungroup()

# write selection frequencies
write_rds(loo_sel_freq, 
          file = paste(out_dir, 'selection-frequencies.rds', sep = '/'))
