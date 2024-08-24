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
library(rslurm)
out_dir <- "/Users/truiz01/Documents/projects/2020-2025/marine-mammal-edna/rslt/loocv/18sv9-ss/_validation"
setwd(out_dir)

## HYPERPARAMETER SPECIFICATION ------------------------------------------------

# hyperparameter grids
eta_grid_res <- 50
eta_grid <- rev(1 - exp(seq(log(0.075), log(0.6), length = eta_grid_res)))
ncomp_grid <- 3:12

# read in combined validation partitions
loo_partitions <- read_rds('partitions.rds')
# partition_ix <- 1:100

# index partitions consecutively
partition_ix <- 1:nrow(loo_partitions)

# form iteration grid (combinations of species and hyperparameter settings)
iter_grid <- expand_grid(eta = eta_grid, 
              ncomp = ncomp_grid, 
              species = c('bm', 'bp', 'mn')) %>%
  rename_with(~paste('.', .x, sep = ''))

# function to iterate over grid
loocv_fn <- function(.eta, .ncomp, .species){
  
  rslt <- mclapply(partition_ix, function(j){
    
    # select data partition
    .partition <- slice(loo_partitions, j)
    
    # extract training and test data
    x.train <- dplyr::select(.partition$train.inner[[1]], starts_with('asv'))
    y.train <- pull(.partition$train.inner[[1]], {{.species}})
    
    # fit model
    fit <- spls(x.train, y.train, 
                K = .ncomp, eta = .eta, 
                scale.x = F, scale.y = F)
    
    # summaries
    y.train.hat <- predict(fit, type = 'fit')
    train.resid <- y.train - y.train.hat
    n <- length(y.train)
    p <- ncol(fit$projection)
    df <- nrow(fit$projection)
    rsq <- as.numeric(1 - (n - p)*var(train.resid)/((n - 1)*var(y.train)))
    sel.asv <- fit$projection %>% rownames()
    
    # squared prediction error
    x.test <- dplyr::select(.partition$test.inner[[1]], starts_with('asv'))
    y.test <- pull(.partition$test.inner[[1]], .species)
    y.test.hat <- predict(fit, newx = x.test, type = 'fit')
    test.resid <- y.test - y.test.hat
    test.spe <- as.numeric(test.resid^2)
    
    # outputs
    out <- tibble(outer.id = .partition$test.outer.cruise,
                  inner.id = .partition$test.inner.cruise,
                  species = .species,
                  eta = .eta,
                  ncomp = .ncomp,
                  sel.asv = list(sel.asv),
                  n.asv = df,
                  adj.rsq = rsq,
                  pred = y.test.hat[,1],
                  pred.err = test.resid[,1],
                  sq.pred.err = test.spe)
    return(out)
  }) %>%
    Reduce(bind_rows, .)
  
  return(rslt)
}

# start <- Sys.time()
# loocv_fn(.eta, .ncomp, .species)
# end <- Sys.time()
# end - start

# for testing
partition_ix <- 1:20
iter_grid_sub <- iter_grid %>% slice(1:100)
sjob <- slurm_apply(f = loocv_fn, params = iter_grid_sub,
                    jobname = 'vtest',
                    global_objects = c('loo_partitions', 'partition_ix'),
                    nodes = 2, 
                    cpus_per_node = 1,
                    submit = F)

# rslt <- get_slurm_out(sjob, outtype = 'raw', wait = F)
