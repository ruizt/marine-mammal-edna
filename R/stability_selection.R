# Stability selection functions
#
# These functions implement the sparse PLS stability selection procedure used
# in the analysis. See vignettes/single-iteration.Rmd for a narrative walkthrough.

#' Run sPLS stability selection across a hyperparameter grid
#'
#' Fits sparse PLS models across all (eta, ncomp) combinations and leave-one-out
#' subsamples, recording which ASVs are selected in each partition. Returns
#' per-ASV selection frequencies across the grid.
#'
#' @param data Data frame with a response variable and ASV predictor columns
#'   (column names starting with \code{"asv"}).
#' @param response Character. Name of the response variable.
#' @param eta_vals Numeric vector of sparsity values (from \code{\link{eta_grid}}).
#' @param ncomp_vals Integer vector of component counts to search.
#' @param partitions LOO partition tibble from \code{\link{loo_partitions}}.
#' @return Tibble with columns \code{eta}, \code{ncomp}, \code{asv}, \code{n}
#'   (selection count across partitions). The number of LOO folds is stored as
#'   attribute \code{n_obs}.
#' @export
run_stability_selection <- function(data, response, eta_vals, ncomp_vals, partitions) {
  n_obs <- nrow(partitions)
  model_grid <- tidyr::expand_grid(ncomp = ncomp_vals, eta = eta_vals)

  sel_raw <- purrr::pmap(model_grid, function(ncomp, eta) {
    lapply(seq_len(n_obs), function(j) {
      train <- partitions$train[[j]]
      x_tr  <- train |> dplyr::select(dplyr::starts_with("asv")) |> as.data.frame()
      y_tr  <- train[[response]]
      fit   <- spls::spls(x_tr, y_tr, K = ncomp, eta = eta,
                          scale.x = FALSE, scale.y = FALSE)
      tibble::tibble(eta = eta, ncomp = ncomp, asv = rownames(fit$projection))
    }) |> dplyr::bind_rows()
  }) |> dplyr::bind_rows()

  sel_freq <- sel_raw |>
    dplyr::group_by(eta, ncomp, asv) |>
    dplyr::count() |>
    dplyr::ungroup()

  attr(sel_freq, "n_obs") <- n_obs
  sel_freq
}

#' Extract stable feature sets from selection frequencies
#'
#' Applies the Meinshausen-Bühlmann stability selection criterion to identify
#' features selected with probability >= \code{pi.max} across a range of eta
#' values, subject to a bound on expected false positives (\code{EV.max}).
#' The mean selection count per partition and the eta hyperparameter grid are
#' inferred directly from \code{sel_freq}.
#'
#' @param sel_freq Tibble from \code{\link{run_stability_selection}}.
#' @param p Integer. Total number of ASV features (columns starting with
#'   \code{"asv"} in the original joined dataset).
#' @param n_obs Integer. Number of LOO folds. If \code{NULL}, read from the
#'   \code{n_obs} attribute of \code{sel_freq} (set automatically by
#'   \code{run_stability_selection}).
#' @param pi.max Numeric. Stability threshold (default 0.8).
#' @param EV.max Numeric. Upper bound on expected false positives (default 1).
#' @param q.min Integer. Minimum stable set size (default 15).
#' @return Tibble of candidate stable sets with columns \code{ncomp},
#'   \code{eta.min}, \code{eta.max}, \code{stable_asvs} (list-col of ASV
#'   character vectors), \code{n_stable}.
#' @export
extract_stable_sets <- function(sel_freq, p, n_obs = NULL,
                                pi.max = 0.8, EV.max = 1, q.min = 15) {
  if (is.null(n_obs)) n_obs <- attr(sel_freq, "n_obs")
  if (is.null(n_obs)) stop("n_obs could not be inferred; pass it explicitly.")

  q.max <- sqrt((2 * pi.max - 1) * p * EV.max)
  eta_grid_vals <- sort(unique(sel_freq$eta))

  # Mean ASVs selected per partition, per (ncomp, eta).
  # sum(n) across ASVs = total selection events = sum over partitions of
  # (# ASVs selected), so sum(n) / n_obs = mean selection count per partition.
  avg_n <- sel_freq |>
    dplyr::group_by(ncomp, eta) |>
    dplyr::summarize(n.asv = sum(n) / n_obs, .groups = "drop")

  ix <- seq(1, length(eta_grid_vals), by = 2)
  eta_intervals <- tidyr::expand_grid(
    eta.min.ix = ix, eta.max.ix = rev(ix)
  ) |>
    dplyr::filter(eta.max.ix - eta.min.ix > 4) |>
    dplyr::mutate(eta.max = eta_grid_vals[eta.max.ix],
                  eta.min = eta_grid_vals[eta.min.ix])

  candidate_ranges <- avg_n |>
    tidyr::nest(model.rslt = -eta) |>
    dplyr::mutate(interval = list(eta_intervals)) |>
    tidyr::unnest(cols = interval) |>
    dplyr::filter(eta >= eta.min, eta <= eta.max) |>
    tidyr::unnest(cols = model.rslt) |>
    collapse::fgroup_by(ncomp, eta.min, eta.max) |>
    collapse::fsummarize(q = collapse::fmean(n.asv), sd = collapse::fsd(n.asv)) |>
    dplyr::filter(q < q.max, q > q.min, sd / q < 1) |>
    dplyr::select(ncomp, eta.min, eta.max)

  candidate_sets <- lapply(seq_len(nrow(candidate_ranges)), function(j) {
    .p <- dplyr::slice(candidate_ranges, j)
    sel_freq |>
      collapse::fsubset(eta >= .p$eta.min & eta <= .p$eta.max &
                          ncomp == .p$ncomp & n >= n_obs * pi.max) |>
      dplyr::pull(asv) |>
      unique()
  }) |>
    (\(x) dplyr::mutate(candidate_ranges, stable_asvs = x))() |>
    dplyr::mutate(n_stable = purrr::map_int(stable_asvs, length)) |>
    dplyr::filter(n_stable >= 2 * ncomp)

  candidate_sets
}
