## CREATE EXAMPLE DATA FOR PACKAGE VIGNETTE -----------------------------------
#
# Backend script run once by developers; outputs committed to the repository so
# end users need no additional data downloads or auxiliary script execution.
#
# Produces: data/processed/example-data-16s.RData  (~141 KB, xz-compressed)
#
# Contents
# --------
# Input data (for the short vignette's data-loading step):
#   example_16s       27 x 501 tibble  cruise + 500 ASV columns
#   example_density   25 x 4 tibble    cruise + log-ratio density (BM/BP/MN)
#   example_asvs      character[500]   ASV IDs in example_16s
#   all_stable        character[93]    stable ASV IDs from full 16S analysis
#                     (included so vignette can verify feature recovery)
#   dens_means        4 x 4 tibble     seasonal log-ratio means (back-transform)
#   dens_raw          75 x 11 tibble   original-scale density estimates
#
# Pre-cached stability selection output for BM (blue whale), full eta grid:
#   sel_freq_bm       tibble  selection frequencies per (eta, ncomp, ASV)
#   avg_n_bm          tibble  mean ASVs selected per (eta, ncomp) config
#   candidate_sets_bm tibble  stable-set candidates with ASV lists
#
# Feature-set subsetting strategy
# --------------------------------
# The short vignette uses one marker (16S) and one species (BM). 500 total
# ASVs were chosen to keep vignette runtimes short. The subset guarantees all
# 93 ASVs identified as stable in the full analysis are present, then pads with
# prevalence-stratified non-stable ASVs. This ensures stability selection on
# the example data recovers a subset of the true stable features with zero
# false positives.
#
# Note: the stability selection false-positive bound q.max scales as
# sqrt(0.6 * p), so with p = 500 (vs. 6234 in the full analysis) the algorithm
# is more conservative and identifies smaller stable sets (~15-25 vs. 35-48
# ASVs). This partial recovery is a theoretical consequence of the reduced
# feature space and is consistent with the sensitivity to data perturbations
# observed in the full analysis. The vignette notes this explicitly.
#
# Pre-cached stability selection
# --------------------------------
# Stability selection for BM using the full hyperparameter grid (55 eta x
# 9 ncomp x 25 LOO folds = 12,375 fits) takes ~3 minutes on the 500-ASV
# example data. The aggregated output (sel_freq_bm, avg_n_bm, candidate_sets_bm)
# is cached here so the vignette can load it directly and demonstrate the
# downstream steps interactively without triggering a multi-minute computation
# at render time. The full computation code is shown in the vignette with a
# note indicating it is not executed.

library(tidyverse)
library(spls)
library(pls)
library(collapse)
data_dir <- 'data/processed/'
stbl_dir <- 'rslt/stability-selection/16s/'
out_file <- paste0(data_dir, 'example-data-16s.RData')

## LOAD SOURCE DATA ------------------------------------------------------------

load(paste0(data_dir, 'ncog16s.RData'))
load(paste0(data_dir, 'density-estimates.RData'))
stable_sets <- read_rds(paste0(stbl_dir, 'stable-sets.rds'))

## FEATURE SUBSET --------------------------------------------------------------

all_stable <- stable_sets$ss |> unlist() |> unique()
cat("Stable ASVs from full analysis:", length(all_stable), "\n")

asv_cols   <- edna |> select(starts_with('asv'))
prevalence <- colMeans(asv_cols > 0)
non_stable <- setdiff(names(prevalence), all_stable)
n_total    <- 500
n_sample   <- n_total - length(all_stable)

set.seed(8243)
sampled_non_stable <- tibble(asv  = non_stable,
                             prev = prevalence[non_stable]) |>
  mutate(quintile = ntile(prev, 5)) |>
  group_by(quintile) |>
  slice_sample(n = ceiling(n_sample / 5)) |>
  ungroup() |>
  slice_head(n = n_sample) |>
  pull(asv)

example_asvs <- c(all_stable, sampled_non_stable)
cat("Example feature set: ", length(example_asvs),
    "(", length(all_stable), "stable +", length(sampled_non_stable), "non-stable )\n")

## INPUT DATA OBJECTS ----------------------------------------------------------

example_16s     <- edna |> select(cruise, all_of(example_asvs))
example_density <- dens
# dens_means and dens_raw loaded from density-estimates.RData above

## PRE-CACHED STABILITY SELECTION (BM, full eta grid) -------------------------

whales_ex <- inner_join(example_density, example_16s, by = 'cruise')
n_ex      <- nrow(whales_ex)

loo_ex <- tibble(
  test.id = whales_ex$cruise,
  train   = map(whales_ex$cruise, ~filter(whales_ex, cruise != .x)),
  test    = map(whales_ex$cruise, ~filter(whales_ex, cruise == .x))
)

eta_grid  <- rev(1 - exp(seq(log(0.075), log(0.6), length = 55)))
ncomp_grid <- 4:12
model_grid_bm <- expand_grid(species = 'bm', ncomp = ncomp_grid, eta = eta_grid)
cat("Fitting", nrow(model_grid_bm) * nrow(loo_ex), "sPLS models for BM...\n")

t0 <- proc.time()
sel_asvs_bm_raw <- pmap(model_grid_bm, function(species, ncomp, eta) {
  lapply(seq_len(nrow(loo_ex)), function(j) {
    train <- loo_ex$train[[j]]
    x_tr  <- train |> select(starts_with('asv')) |> as.data.frame()
    y_tr  <- train |> pull(bm)
    fit   <- spls(x_tr, y_tr, K = ncomp, eta = eta, scale.x = FALSE, scale.y = FALSE)
    tibble(species = 'bm', eta = eta, ncomp = ncomp,
           obs.id  = loo_ex$test.id[[j]],
           sel.asv = list(rownames(fit$projection)))
  }) |> bind_rows()
}) |> bind_rows()
cat("Runtime:", round((proc.time() - t0)[['elapsed']] / 60, 1), "min\n")

# Selection frequencies (aggregated; sufficient for all downstream vignette steps)
sel_freq_bm <- sel_asvs_bm_raw |>
  unnest(sel.asv) |>
  group_by(species, eta, ncomp, sel.asv) |>
  count() |>
  ungroup()

# Average number of ASVs selected per configuration
avg_n_bm <- sel_asvs_bm_raw |>
  mutate(n.asv = map_int(sel.asv, length)) |>
  fgroup_by(species, ncomp, eta) |>
  fselect(n.asv) |>
  fmean()

# Stable-set candidate construction (mirrors full analysis logic)
p_ex    <- ncol(whales_ex |> select(starts_with('asv')))
pi.max  <- 0.8
EV.max  <- 1
q.max   <- sqrt((2 * pi.max - 1) * p_ex * EV.max)
q.min   <- 15

eta_intervals <- expand_grid(
  eta.min.ix = seq(1, length(eta_grid), by = 2),
  eta.max.ix = rev(seq(1, length(eta_grid), by = 2))
) |>
  filter(eta.max.ix - eta.min.ix > 4) |>
  mutate(eta.max = eta_grid[eta.max.ix],
         eta.min = eta_grid[eta.min.ix])

candidate_ranges_bm <- avg_n_bm |>
  nest(model.rslt = -eta) |>
  mutate(interval = list(eta_intervals)) |>
  unnest(cols = interval) |>
  filter(eta >= eta.min, eta <= eta.max) |>
  unnest(cols = model.rslt) |>
  fgroup_by(species, ncomp, eta.min, eta.max) |>
  fsummarize(q = fmean(n.asv), sd = fsd(n.asv)) |>
  filter(q < q.max, q > q.min, sd / q < 1) |>
  select(species, ncomp, eta.min, eta.max)

candidate_sets_bm <- lapply(seq_len(nrow(candidate_ranges_bm)), function(j) {
  .p <- slice(candidate_ranges_bm, j)
  sel_freq_bm |>
    fsubset(eta >= .p$eta.min & eta <= .p$eta.max &
              ncomp == .p$ncomp & n >= n_ex * pi.max &
              species == .p$species) |>
    pull(sel.asv) |>
    unique()
}) |>
  (\(x) mutate(candidate_ranges_bm, ss = x))() |>
  group_by(across(-c(eta.min, eta.max))) |>
  slice_min(eta.max - eta.min) |>
  ungroup() |>
  mutate(n.asv = map_int(ss, length)) |>
  filter(n.asv >= 2 * ncomp)

# Recovery check (informational)
bm_union   <- candidate_sets_bm$ss |> unlist() |> unique()
bm_stable  <- stable_sets |> filter(species == 'bm') |> pull(ss) |> unlist()
cat("\nBM recovery:", sum(bm_stable %in% bm_union), "/", length(bm_stable),
    "stable ASVs |", sum(!bm_union %in% bm_stable), "false positives\n")

## SAVE ------------------------------------------------------------------------

save(example_16s, example_density, example_asvs, all_stable,
     dens_means, dens_raw,
     sel_freq_bm, avg_n_bm, candidate_sets_bm,
     file = out_file, compress = 'xz')

cat("\nSaved:", out_file, "\n")
cat("File size:", round(file.size(out_file) / 1024), "KB\n")
