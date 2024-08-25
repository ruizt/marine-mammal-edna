library(tibble)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(fs)
library(readr)
library(modelr)
library(spls)
library(doParallel)
library(foreach)
# library(rslurm)
part_dir <- "rslt/loocv/18sv9-ss/_validation"
out_dir <- "rslt/loocv/18sv9-ss/_validation/results"

## HYPERPARAMETER SPECIFICATION ------------------------------------------------

# hyperparameter grids
eta_grid_res <- 50
eta_grid <- rev(1 - exp(seq(log(0.075), log(0.6), length = eta_grid_res)))
ncomp_grid <- 3:12

# read in combined validation partitions
loo_partitions <- paste(part_dir, '/partitions.rds', sep = '') |> read_rds()
# partition_ix <- 1:100

# index partitions consecutively
partition_ix <- 1:nrow(loo_partitions)

# form iteration grid (combinations of species and hyperparameter settings)
iter_grid <- expand_grid(eta = eta_grid, 
                         ncomp = ncomp_grid, 
                         species = c('bm', 'bp', 'mn'))

# function to iterate over grid
loocv_fn <- function(.eta, .ncomp, .species){
  
  rslt <- lapply(partition_ix, function(j){
    
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
  }) 
  
  out <- rslt %>%
    Reduce(bind_rows, .)
  
  return(out)
}

# chunk_ix <- seq(1, nrow(iter_grid), length = 16) |> round()
# partition_ix <- 1:5
chunk_ix <- seq(1, 200, length = detectCores() + 1) |> round()
registerDoParallel(cores=detectCores() + 1)
foreach(k = 1:(length(chunk_ix) - 1)) %dopar% {
  rslt <- lapply(chunk_ix[k]:(chunk_ix[k + 1] - 1),
                  function(i){
                    loocv_fn(iter_grid$eta[i],
                             iter_grid$ncomp[i],
                             iter_grid$species[i])
                  })
  out <- Reduce(bind_rows, rslt)
  # out <- rnorm(100)
  write_rds(out, paste(out_dir, '/rslt_', k, '.rds', sep = ''))
}
