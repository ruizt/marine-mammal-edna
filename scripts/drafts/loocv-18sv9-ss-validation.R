library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(fs)
library(readr)
library(modelr)
library(spls)
library(parallel)

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



## NESTED LEAVE ONE OUT PROCEDURE FOR MODEL VALIDATION -------------------------

# clear outer loocv objects from environment
rm(list = setdiff(ls(), c('eta_grid', 'ncomp_grid', 'loocv_fn')))

# directory to store files
dir <- 'rslt/loocv/18sv9-ss/inner'
dir_create(dir)

# # read in and combine validation partitions
# load('data/processed/_cv-partitions/ncog18sv9-cv-inner.RData')
# load('data/processed/_cv-partitions/mm-sightings-cv-inner.RData')
# loo_partitions <- inner_join(loo_sightings_inner,
#                              loo_edna_inner, 
#                              join_by(test.outer.cruise, test.inner.cruise), 
#                              suffix = c('.mm', '.edna')) |>
#   mutate(train.inner = map2(train.inner.mm, train.inner.edna, 
#                             ~inner_join(.x, .y, join_by('cruise'))),
#          test.inner = map2(test.inner.mm, test.inner.edna, 
#                            ~inner_join(.x, .y, join_by('cruise')))) |>
#   select(test.outer.cruise, test.inner.cruise, test.inner.season,
#          train.inner, test.inner)
# 
# # store partitions
# write_rds(loo_partitions, 
#           file = paste(dir, 'partitions.rds', sep = '/'))

# read in combined validation partitions
loo_partitions <- read_rds(paste(dir, 'partitions.rds', sep = '/'))

# iteration grids
species_grid <- c("bm", "bp", "mn")
outer_partition_grid <- loo_partitions |>  pull(test.outer.cruise) |> unique()
inner_partition_ix <- 1:24

# index inner partition and eta grid combos for innermost iterations
settings <- expand_grid(eta = eta_grid,
                        part_ix = inner_partition_ix) |>
  mutate(setting = row_number()) |>
  select(setting, eta, part_ix)


# uncomment for testing
# .outer <- outer_partition_grid[1]
# .species <- species_grid[1]
# .ncomp <- ncomp_grid[1]
# i <- 1

# iterate over outermost data partitions
for(.outer in outer_partition_grid){
  
  # create outer-partition-specific filepath
  path <- paste(dir, '/_', .outer, sep = '')
  dir_create(path)
  
  # filter to subset of loocv partitions
  subpartitions <- filter(loo_partitions, test.outer.cruise == .outer)
  
  # iterate over species
  for(.species in species_grid){
    
    # create species-specific filepath as subpath
    subpath <- paste(path, .species, sep = '/')
    dir_create(subpath)
    
    # iterate over numbers of components
    lapply(ncomp_grid, function(.ncomp){
      
      # iterate over sparsity parameter and inner data partitions
      rslt <- lapply(settings$setting, function(i){
        
        # identify data partition and sparsity parameter
        setting <- slice(settings, i)
        partition <- setting$part_ix
        eta <- setting$eta
        
        # retrieve corresponding data partition
        train <- subpartitions$train.inner[partition][[1]] %>% as.data.frame() 
        test <- subpartitions$test.inner[partition][[1]] %>% as.data.frame() 
        
        # fit model and return fit, prediction, and selection metrics
        out <- loocv_fn(train, 
                        test, 
                        .var = .species, 
                        .parms = c(ncomp = .ncomp, eta = eta))
        
        # print status update to console
        paste('species = ', .species,
              ', outer = ', .outer,
              ', ncomp = ', .ncomp,
              ', eta = ', round(eta, 4), 
              ', inner = ', partition, sep = '') |> print()
        return(out)
      })
      
      # combine and write fit and prediction metrics to file
      sapply(settings$setting, function(i){
        c(setting = i, 
          rslt[[i]]$parms, 
          rslt[[i]]$metrics)
      }) |>
        t() |> 
        as_tibble() |>
        write_rds(file = paste(subpath, '/', .ncomp, 'comp-metrics.rds', sep = ''))
      
      # combine and write selected asv lists to file
      lapply(settings$setting, function(i){rslt[[i]]$sel.asv}) |>
        write_rds(file = paste(subpath, '/', .ncomp, 'comp-selected.rds', sep = ''))
    })
  }
}

# read in and combine selected asvs across species
subdirs <- list.files(dir) |> str_remove_all('_')
outer_partition_grid <- subdirs[-(26:28)]
# .outer <- outer_partition_grid[1]
# .species <- species_grid[1]

loo_sel_asv <- lapply(outer_partition_grid, function(.outer){
  lapply(species_grid, function(.species){
    
    path <- paste(dir, '/_', .outer, '/', .species, sep = '')
    inner_part_id <- loo_partitions |> 
      filter(test.outer.cruise == .outer) |>
      pull(test.inner.cruise)
    
    # read and combine selected asv lists across numbers of components
    asv_list <- lapply(ncomp_grid, function(.ncomp){
      read_rds(file = paste(path, '/', .ncomp, 'comp-selected.rds', sep = ''))
    }) %>%
      Reduce(c, .)
    
    # format as tibble
    out <- expand_grid(ncomp = ncomp_grid, 
                       eta = eta_grid, 
                       inner.part.id = inner_part_id) |>
      mutate(asv = asv_list,
             species = .species,
             outer.part.id = .outer)
    
    return(out)
  }) %>%
    Reduce(bind_rows, .) 
}) %>%
  Reduce(bind_rows, .) 

# export to file
write_rds(loo_sel_asv, file = paste(dir, 'sel-asv.rds', sep = '/'))

# tabulate selection frequencies
loo_sel_freq <- loo_sel_asv |>
  mutate(df = map(asv, length)) |>
  unnest(df) |>
  filter(df > 1) |>
  unnest(asv) |>
  group_by(outer.part.id, species, ncomp, eta, asv) |>
  count() |>
  ungroup()

# write selection frequencies
write_rds(loo_sel_freq, 
          file = paste(dir, 'selection-frequencies.rds', sep = '/'))

# read in and combine loocv run metrics across species
loo_metrics <- lapply(outer_partition_grid, function(.outer){
  lapply(species_grid, function(.species){
    
    path <- paste(dir, '/_', .outer, '/', .species, sep = '')
    inner_part_id <- loo_partitions |> 
      filter(test.outer.cruise == .outer) |>
      pull(test.inner.cruise)
    
    # read and combine fit/prediction metrics across numbers of components
    out <- lapply(ncomp_grid, function(.ncomp){
      read_rds(file = paste(path, '/', .ncomp, 'comp-metrics.rds', sep = '')) 
    }) %>% 
      Reduce(bind_rows, .) |>
      mutate(species = .species,
             outer.part.id = .outer,
             inner.part.id = rep(inner_part_id, length(ncomp_grid)*length(eta_grid))) 
    return(out)
  }) %>%
    Reduce(bind_rows, .)
}) %>%
  Reduce(bind_rows, .)

# export to file
write_rds(loo_metrics, file = paste(dir, 'metrics.rds', sep = '/'))

