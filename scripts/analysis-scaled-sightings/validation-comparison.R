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
data_dir <- 'data/processed/'
paste(data_dir, 'ncog18sv9.RData', sep = '') |> load()
paste(data_dir, 'mm-sightings.RData', sep = '') |> load() 
whales <- inner_join(sightings, edna, by = 'cruise') 
p <- whales |> dplyr::select(starts_with('asv')) |> ncol()
n.obs <- nrow(whales)

## VALIDATION FOR STABILITY SELECTION ------------------------------------------

# read in results of leave one out procedure
loo_dir <- 'rslt/loocv/18sv9-ss/inner/'
loo_metrics <- read_rds(paste(loo_dir, 'metrics.rds', sep = ''))
loo_sel_asv <- read_rds(paste(loo_dir, 'sel-asv.rds', sep = ''))
loo_sel_freq <- read_rds(paste(loo_dir, 'selection-frequencies.rds', sep = ''))
loo_partitions <- read_rds(paste(loo_dir, 'partitions.rds', sep = ''))

# retrieve hyperparameter and partition grids
eta_grid <- loo_metrics |> pull(eta) |> unique()
outer_partition_grid <- loo_metrics |> pull(outer.part.id) |> unique()

# fix stability threshold (minimax selection prob.)
pi.max <- 0.8

# fix max average number of selected variables (controls false positive rate)
EV.max <- 0.5
q.max <- sqrt((2*pi.max - 1)*p*EV.max)
q.max

# fix min average number of selected variables (models are not too sparse)
q.min <- 10

# generate sets of candidate eta intervals
eta_intervals <- expand_grid(eta.min = seq_range(eta_grid, by = 0.025),
                             eta.max = rev(eta.min)) |>
  filter(eta.max - eta.min > 0.1) 

# function to compute stable sets from a parameter (eta/ncomp) range
stableset_fn <- function(parms){
  out <- loo_sel_freq |> 
    filter(eta >= parms$eta.min,
           eta <= parms$eta.max,
           ncomp == parms$ncomp,
           n >= n.obs*pi.max,
           species == parms$species,
           outer.part.id == parms$outer.part.id) |>
    pull(asv) |>
    unique()
  
  return(out)
}

# # testing
# stableset_fn(data.frame(eta.min = 0.55, eta.max = 0.65, ncomp = 6, species = "bm", outer.part.id = outer_partition_grid[1]))

# numbers of asvs selected for each data partition and hyperparameter
df <- loo_metrics |> 
  select(outer.part.id, inner.part.id, species, ncomp, eta, n.asv) 

# identify intervals for the sparsity parameter that satisfy qmin, qmax criteria
parm_ranges <- df |>
  group_by(outer.part.id, species, ncomp, eta) |>
  summarize(q = mean(n.asv)) |>
  mutate(interval = list(eta_intervals)) |>
  unnest(cols = interval) |>
  filter(eta >= eta.min, eta <= eta.max) |>
  group_by(outer.part.id, species, ncomp, eta.min, eta.max) |>
  summarize(qmax = max(q),
            qmean = mean(q),
            qsd = sd(q)) |>
  filter(qmean < q.max, qmean > q.min, qsd/qmean < 1) |>
  ungroup() |>
  group_by(species, ncomp, eta.min, eta.max) |>
  count() |>
  filter(n == 25) |> # only pick eta ranges identified for *every* data partition
  select(-n) |>
  mutate(outer.part.id = list(outer_partition_grid)) |>
  unnest(outer.part.id) |>
  ungroup()

# for each parameter range, compute the stable set (in chunks)
fs::dir_create(paste(dir, '/stable-sets', sep = ''))
ix_partitions <- seq(1, nrow(parm_ranges), length = 50) |> round()
for(j in 1:(length(ix_partitions) - 1)){
  stable_sets <- mclapply(ix_partitions[j]:ix_partitions[j + 1], function(i){
    slice(parm_ranges, i) |> stableset_fn()
  }, mc.cores = detectCores())
  write_rds(stable_sets, file = paste(dir, '/stable-sets/ss-chunk-', j, '.rds', sep = ''))
  print(j)
}

# read in and combine results from above
stable_sets <- lapply(1:49, function(j){
  ss <- read_rds(paste(loo_dir, '/stable-sets/ss-chunk-', j, '.rds', sep = ''))
  if(j < 49){
    out <- ss[-length(ss)]
  }else{
    out <- ss
  }
}) %>%
  Reduce(c, .)

# function to compute no. of elements included in at least <thresh>% of sets in <set list>
intersect_fn <- function(set_list, thresh){
  fint <- tibble(asv = Reduce(c, set_list)) |> 
    group_by(asv) |> 
    count() |>
    filter(n >= thresh*length(set_list)) |>
    pull(asv)
  out <- length(fint)
  return(out)
}

# function to compute no. of elements in union of sets in <set_list>
union_fn <- function(set_list){
  int <- Reduce(c, set_list) |> unique()
  out <- length(int)
  return(out)
}

# consistency of stable sets across outer partitions
ss_consistency <- parm_ranges |>
  mutate(ss = stable_sets) |>
  group_by(species, ncomp, eta.min, eta.max) |>
  summarize(int.length = intersect_fn(ss, 1),
            fint.length = intersect_fn(ss, 0.8),
            un.length = union_fn(ss),
            avg.length = mean(sapply(ss, length)),
            n = n()) |>
  mutate(int.to.avg = int.length/avg.length,
         fint.to.avg = fint.length/avg.length,
         fj.index = fint.length/un.length) 

# table of fuzzy jaccard index by ncomp, species
fjix_tbl_ss <- ss_consistency |>
  group_by(species, ncomp) |>
  summarize(mean.fj.index = mean(fj.index)) |>
  mutate(metric = 'J index',
         method = 'SS') |>
  pivot_wider(names_from = 'ncomp', values_from = 'mean.fj.index')

# table of 80%-intersection to average size of stable set ratio
ita_tbl_ss <- ss_consistency |>
  group_by(species, ncomp) |>
  summarize(mean.fint = mean(fint.to.avg)) |>
  mutate(metric = 'IA ratio',
         method = 'SS') |>
  pivot_wider(names_from = 'ncomp', values_from = 'mean.fint')

# average sizes of core sets
core_tbl_ss <- ss_consistency |>
  group_by(species, ncomp) |>
  summarize(mean.fint = mean(fint.length)) |>
  mutate(metric = 'Avg. core set size',
         method = 'SS') |>
  pivot_wider(names_from = 'ncomp', values_from = 'mean.fint')


## SPLS VALIDATION -------------------------------------------------------------

# read in results of leave one out procedure
loo_dir <- 'rslt/loocv/18sv9-ss/outer/'
loo_metrics <- read_rds(paste(loo_dir, 'metrics.rds', sep = ''))
loo_sel_asv <- read_rds(paste(loo_dir, 'sel-asv.rds', sep = ''))
loo_sel_freq <- read_rds(paste(loo_dir, 'selection-frequencies.rds', sep = ''))
loo_partitions <- read_rds(paste(loo_dir, 'partitions.rds', sep = ''))

# compute consistency of selected asv's across data partitions
spls_consistency <- loo_sel_asv |>
  group_by(species, ncomp, eta) |>
  summarize(int.length = intersect_fn(asv, 1),
            fint.length = intersect_fn(asv, 0.8),
            un.length = union_fn(asv),
            avg.length = mean(sapply(asv, length)),
            n = n()) |>
  mutate(int.to.avg = int.length/avg.length,
         fint.to.avg = fint.length/avg.length,
         fj.index = fint.length/un.length)

# table of fuzzy jaccard index by ncomp, species
fjix_tbl_spls <- spls_consistency |>
  group_by(species, ncomp) |>
  summarize(mean.fj.index = mean(fj.index)) |>
  mutate(metric = 'J index',
         method = 'SPLS') |>
  pivot_wider(names_from = 'ncomp', values_from = 'mean.fj.index')

# table of 80%-intersection to average size of stable set ratio
ita_tbl_spls <- spls_consistency |>
  group_by(species, ncomp) |>
  summarize(mean.fint = mean(fint.to.avg)) |>
  mutate(metric = 'IA ratio',
         method = 'SPLS') |>
  pivot_wider(names_from = 'ncomp', values_from = 'mean.fint')

# average sizes of core sets
core_tbl_spls <- spls_consistency |>
  group_by(species, ncomp) |>
  summarize(mean.fint = mean(un.length)) |>
  mutate(metric = 'Avg. core set size',
         method = 'SPLS') |>
  pivot_wider(names_from = 'ncomp', values_from = 'mean.fint')

## RESULTS TABLE ---------------------------------------------------------------

bind_rows(fjix_tbl_ss, fjix_tbl_spls, ita_tbl_ss, ita_tbl_spls) |>
  write_csv('rslt/loocv/18sv9-ss/validation-comparison.csv')
