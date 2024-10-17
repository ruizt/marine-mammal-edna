library(tidyverse)
library(fs)
library(modelr)
library(spls)
library(pls)
library(collapse)
library(magrittr)
data_dir <- 'data/processed/'
stbl_dir <- 'rslt/stability-selection/18sv4-ss/'
out_dir <- 'rslt/nested-validation/18sv4-ss/'
dir_create(out_dir)
paste(out_dir, 'inner-spls-fits/', sep = '_') |> dir_create()

## DATA INPUTS -----------------------------------------------------------------

# read in and combine 18sv4 data and scaled sightings
paste(data_dir, 'ncog18sv4.RData', sep = '') |> load()
paste(data_dir, 'mm-sightings.RData', sep = '') |> load() 
whales <- inner_join(sightings, edna, by = 'cruise')

# dimensions
p <- whales |> dplyr::select(starts_with('asv')) |> ncol()
n.obs <- nrow(whales)

# read in nested leave one out partitions
partitions_nested <- paste(data_dir, 
                           '_combined-partitions/nested-partitions-18sv4-ss.rds', 
                           sep = '') |> 
  read_rds()

# import optimal settings from stability selection
stable_sets <- paste(stbl_dir, 'stable-sets.rds', sep = '') |> read_rds()

# retrieve relevant portion of model grid
model_grid <- paste(stbl_dir, 'model-grid.rds', sep = '') |> 
  read_rds() |>
  left_join(stable_sets, join_by(species), suffix = c('', '.opt')) |>
  filter(ncomp == ncomp.opt,
         eta <= eta.max,
         eta >= eta.min) |>
  select(species, ncomp, eta)

## FIT MODELS TO INNER PARTITIONS ----------------------------------------------

# check no. models to fit
nrow(partitions_nested)*nrow(model_grid)

# function to fit spls model on inner partitions and extract selected asvs
fit_fn <- function(.species, .eta, .ncomp){
  
  out <- lapply(1:nrow(partitions_nested), function(j){
    
    # retrieve training data
    .train <- partitions_nested$train[j][[1]]
    
    # separate predictors and response
    x.train <- .train |> as.data.frame() |> select(starts_with('asv'))
    y.train <- .train |> as.data.frame() |> pull(.species)
    
    # fit spls model with specified parameters
    .fit <- spls(x.train, y.train, 
                 K = .ncomp, eta = .eta, 
                 scale.x = F, scale.y = F)
    
    # selected variables
    .asv <- .fit$projection |> rownames()
    
    # outputs
    out <- tibble(species = .species,
                  eta = .eta,
                  ncomp = .ncomp,
                  outer.id = partitions_nested$outer.id[j],
                  inner.id = partitions_nested$inner.id[j],
                  asv = list(.asv))
    paste('species = ', .species, 
          ', eta = ', round(.eta, 4), 
          ', ncomp = ', .ncomp, 
          ', partition ', partitions_nested$outer.id[j],
          '-', partitions_nested$inner.id[j],
          sep = '') |>
      print()
    return(out)
  }) %>%
    Reduce(bind_rows, .)
  
  return(out)
}

# fit each model on every nested partition
## LONG RUNTIME
for(i in 1:nrow(model_grid)){
  fit_fn(.species = model_grid$species[i],
         .eta = model_grid$eta[i],
         .ncomp = model_grid$ncomp[i]) |>
    write_rds(file = paste(out_dir, 'inner-spls-fits/', 
                           model_grid$species[i], 
                           model_grid$ncomp[i], 
                           round(model_grid$eta[i], 4), 
                           '.rds', sep = '_'))
}

## STABILITY SELECTION ON INNER PARTITIONS -------------------------------------

# retrieve selected asvs from procedure above
sel_asvs <- paste(out_dir, '_inner-spls-fits/', sep = '') |>
  dir_ls() |>
  lapply(read_rds) %>%
  Reduce(bind_rows, .)

# compute selection frequencies across inner partitions for each model configuration
sel_freq <- sel_asvs |>
  unnest(asv) |>
  group_by(outer.id, species, eta, ncomp, asv) |>
  count() |>
  ungroup()

# stability threshold (minimax selection prob.)
pi.max <- 0.8

# compute stable sets for each outer partition
validation_stable_sets <- sel_freq |>
  filter(n >= pi.max*24) |>
  select(outer.id, species, asv) |>
  nest(asv = asv) |>
  mutate(asv = map(asv, ~pull(.x, asv))) |>
  arrange(species, outer.id) 

# export
write_rds(validation_stable_sets, 
          file = paste(out_dir, 'validation-stable-sets.rda', sep = ''))

# ## ASSESS STABLE SET CONSISTENCY -----------------------------------------------
# 
# # function to compute "soft intersection" of sets in <set_list>
# intersect_fn <- function(set_list, thresh){
#   out <- tibble(asv = Reduce(c, set_list)) |> 
#     group_by(asv) |> 
#     count() |>
#     filter(n >= thresh*length(set_list)) |>
#     pull(asv)
#   return(out)
# }
# 
# # function to compute no. of elements in union of sets in <set_list>
# union_fn <- function(set_list){
#   out <- Reduce(c, set_list) |> unique()
#   return(out)
# }
# 
# # general consistency across all candidate sets considered
# inner_join(candidate_sets, best_settings) |>
#   select(species, ss)  |>
#   group_by(species) |>
#   summarize(int = intersect_fn(ss, 0.5) |> list(),
#             un = union_fn(ss) |> list()) |>
#   mutate(j.index = map2(int, un, ~length(.x)/length(.y))) |>
#   unnest(j.index) 


