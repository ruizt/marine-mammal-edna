# Stability selection functions
#
# These functions implement the sparse PLS stability selection procedure used
# in the analysis. See analysis/stability-selection-*.R for how they are called,
# and vignettes/single-iteration.Rmd for a narrative walkthrough.

#' Fit sPLS models across leave-one-out partitions
#'
#' For a single (species, eta, ncomp) combination, fits a sparse PLS model to
#' each training partition and returns the selected ASVs per partition.
#'
#' @param species Character. Response variable name (e.g. "bm", "bp", "mn").
#' @param eta Numeric. sPLS sparsity parameter (0, 1).
#' @param ncomp Integer. Number of PLS components.
#' @param data Data frame. Full joined eDNA + density data.
#' @param partitions Tibble. LOO partition table with columns train, test, test.id.
#' @return Tibble with columns species, eta, ncomp, obs.id, sel.asv (list-col).
#' @export
fit_spls_partition <- function(species, eta, ncomp, data, partitions) {
  lapply(seq_len(nrow(partitions)), function(j) {
    train  <- partitions$train[[j]]
    x_tr   <- train |> dplyr::select(dplyr::starts_with('asv')) |> as.data.frame()
    y_tr   <- train |> dplyr::pull(!!species)
    fit    <- spls::spls(x_tr, y_tr, K = ncomp, eta = eta,
                         scale.x = FALSE, scale.y = FALSE)
    tibble::tibble(
      species = species, eta = eta, ncomp = ncomp,
      obs.id  = partitions$test.id[[j]],
      sel.asv = list(rownames(fit$projection))
    )
  }) |> dplyr::bind_rows()
}

#' Compute selection frequencies from raw sPLS partition fits
#'
#' @param sel_asvs Tibble of raw sPLS outputs from \code{fit_spls_partition()}.
#' @return Tibble with columns species, eta, ncomp, sel.asv, n (selection count).
#' @export
selection_frequencies <- function(sel_asvs) {
  sel_asvs |>
    tidyr::unnest(sel.asv) |>
    dplyr::group_by(species, eta, ncomp, sel.asv) |>
    dplyr::count() |>
    dplyr::ungroup()
}

#' Extract stable feature sets from selection frequencies
#'
#' Applies the Meinshausen-Bühlmann stability selection criterion to identify
#' features selected with probability >= pi.max across a range of eta values,
#' subject to a bound on expected false positives (EV.max).
#'
#' @param sel_freq Tibble from \code{selection_frequencies()}.
#' @param avg_n Tibble of mean selected features per (species, ncomp, eta).
#' @param eta_grid Numeric vector of eta values used in the grid.
#' @param n_obs Integer. Number of observations (LOO folds).
#' @param p Integer. Total number of features (ASVs).
#' @param pi.max Numeric. Stability threshold (default 0.8).
#' @param EV.max Numeric. Upper bound on expected false positives (default 1).
#' @param q.min Integer. Minimum stable set size (default 15).
#' @return Tibble of candidate stable sets with columns species, ncomp,
#'   eta.min, eta.max, ss (list-col of ASV character vectors), n.asv.
#' @export
extract_stable_sets <- function(sel_freq, avg_n, eta_grid, n_obs, p,
                                 pi.max = 0.8, EV.max = 1, q.min = 15) {
  q.max <- sqrt((2 * pi.max - 1) * p * EV.max)

  eta_intervals <- tidyr::expand_grid(
    eta.min.ix = seq(1, length(eta_grid), by = 2),
    eta.max.ix = rev(seq(1, length(eta_grid), by = 2))
  ) |>
    dplyr::filter(eta.max.ix - eta.min.ix > 4) |>
    dplyr::mutate(eta.max = eta_grid[eta.max.ix],
                  eta.min = eta_grid[eta.min.ix])

  candidate_ranges <- avg_n |>
    tidyr::nest(model.rslt = -eta) |>
    dplyr::mutate(interval = list(eta_intervals)) |>
    tidyr::unnest(cols = interval) |>
    dplyr::filter(eta >= eta.min, eta <= eta.max) |>
    tidyr::unnest(cols = model.rslt) |>
    collapse::fgroup_by(species, ncomp, eta.min, eta.max) |>
    collapse::fsummarize(q = collapse::fmean(n.asv), sd = collapse::fsd(n.asv)) |>
    dplyr::filter(q < q.max, q > q.min, sd / q < 1) |>
    dplyr::select(species, ncomp, eta.min, eta.max)

  candidate_sets <- lapply(seq_len(nrow(candidate_ranges)), function(j) {
    .p <- dplyr::slice(candidate_ranges, j)
    sel_freq |>
      collapse::fsubset(eta >= .p$eta.min & eta <= .p$eta.max &
                          ncomp == .p$ncomp & n >= n_obs * pi.max &
                          species == .p$species) |>
      dplyr::pull(sel.asv) |>
      unique()
  }) |>
    (\(x) dplyr::mutate(candidate_ranges, ss = x))() |>
    dplyr::group_by(dplyr::across(-c(eta.min, eta.max))) |>
    dplyr::slice_min(eta.max - eta.min) |>
    dplyr::ungroup() |>
    dplyr::mutate(n.asv = purrr::map_int(ss, length)) |>
    dplyr::filter(n.asv >= 2 * ncomp)

  candidate_sets
}
