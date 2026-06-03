# Stability selection functions
#
# These functions implement the sparse PLS stability selection procedure
# described in Satterthwaite et al. See vignettes/single-iteration.qmd
# for a narrative walkthrough.

#' Run sPLS stability selection across a hyperparameter grid
#'
#' Fits sparse PLS (sPLS) models for every combination of sparsity (`eta`) and
#' component count (`ncomp`) on each leave-one-out (LOO) training subset, and
#' records which ASVs are selected in each fit. The output is a tally of how
#' many times each ASV was selected across all LOO folds for each
#' `(eta, ncomp)` configuration.
#'
#' @details
#' Stability selection (Meinshausen & Bühlmann, 2010) identifies features whose
#' inclusion in a model is robust to small perturbations of the data. This
#' implementation uses LOO subsampling: for each hyperparameter setting, sPLS
#' is fit to each of \eqn{n} training sets formed by leaving one observation
#' out, and the selected variables are recorded. Features that are selected in
#' a large fraction of subsamples (high selection frequency) are candidates for
#' inclusion in the final stable set.
#'
#' In the full analysis (Satterthwaite et al.), the hyperparameter grid uses
#' 55 eta values from \code{\link{eta_grid}()} spanning approximately 0.4--0.925,
#' and `ncomp` values 4--12, giving 495 configurations per LOO fold.
#' For 25 cruises, this is \eqn{25 \times 495 = 12{,}375} sPLS fits per
#' species.
#'
#' The sPLS models are fit via \code{\link[spls]{spls}} with
#' \code{scale.x = FALSE} and \code{scale.y = FALSE}. The input data should
#' already be appropriately transformed (e.g., CLR-transformed eDNA
#' abundances and log-ratio density estimates).
#'
#' @param data Data frame containing a response column and ASV predictor
#'   columns. Predictor columns must have names starting with \code{"asv"}.
#' @param response Character string. Name of the response column in
#'   \code{data} (e.g., \code{"bm"} for blue whale log-ratio density).
#' @param eta_vals Numeric vector of sparsity values to search over.
#'   Higher values produce sparser models (fewer selected features).
#'   Use \code{\link{eta_grid}()} to generate the default grid matching the
#'   paper analysis.
#' @param ncomp_vals Integer vector of PLS component counts to search
#'   (e.g., \code{4:12}).
#' @param partitions LOO partition tibble from \code{\link{loo_partitions}},
#'   with list-columns \code{train} and \code{test}.
#'
#' @return A tibble with four columns:
#'   \describe{
#'     \item{eta}{Numeric. Sparsity parameter value.}
#'     \item{ncomp}{Integer. Number of PLS components.}
#'     \item{asv}{Character. ASV identifier.}
#'     \item{n}{Integer. Number of LOO folds (out of \code{nrow(partitions)})
#'       in which this ASV was selected at this \code{(eta, ncomp)} setting.}
#'   }
#'   The total number of LOO folds is stored as attribute \code{"n_obs"}.
#'
#' @seealso \code{\link{extract_stable_sets}} for the next pipeline step,
#'   \code{\link{eta_grid}} and \code{\link{loo_partitions}} for helper
#'   functions, \code{\link[spls]{spls}} for the underlying sPLS implementation.
#'
#' @references
#' Meinshausen, N. & Bühlmann, P. (2010). Stability selection.
#' \emph{Journal of the Royal Statistical Society: Series B}, 72(4), 417--473.
#' \doi{10.1111/j.1467-9868.2010.00740.x}
#'
#' @examples
#' \dontrun{
#' library(satterthwaite2026)
#' library(dplyr)
#'
#' whales <- inner_join(example_density, example_16s, by = "cruise")
#' partitions <- loo_partitions(whales)
#'
#' sel_freq <- run_stability_selection(
#'   whales, response = "bm",
#'   eta_vals   = eta_grid(n = 55),
#'   ncomp_vals = 4:12,
#'   partitions = partitions
#' )
#' }
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
#' Scans over contiguous intervals of the eta grid and, for each number of
#' PLS components, identifies hyperparameter regions where the mean number of
#' selected features satisfies the Meinshausen--Bühlmann false-positive bound.
#' Within each qualifying region, features with selection frequency
#' \eqn{\geq \pi_{\max}}{>= pi.max} are declared stable. The result is a
#' table of candidate stable sets for downstream evaluation.
#'
#' @details
#' The false-positive control derives from Theorem 1 of Meinshausen &
#' Bühlmann (2010). The expected number of falsely selected variables is
#' bounded by:
#'
#' \deqn{E[V] \leq \frac{q^2}{(2\pi_{\max} - 1) \cdot p}}{
#'       E[V] <= q^2 / ((2 * pi.max - 1) * p)}
#'
#' where \eqn{q} is the mean number of features selected per subsample across
#' the eta interval, \eqn{p} is the total number of candidate features, and
#' \eqn{\pi_{\max}}{pi.max} is the selection frequency threshold. Requiring
#' \eqn{E[V] \leq}{E[V] <=} \code{EV.max} gives the upper bound:
#'
#' \deqn{q_{\max} = \sqrt{(2\pi_{\max} - 1) \cdot p \cdot E_{\max}[V]}}{
#'       q.max = sqrt((2 * pi.max - 1) * p * EV.max)}
#'
#' Only eta intervals where the mean selection count \eqn{q} satisfies
#' \eqn{q_{\min} < q < q_{\max}}{q.min < q < q.max} and the coefficient of
#' variation is less than 1 are retained. Within each retained interval,
#' ASVs appearing in \eqn{\geq \pi_{\max} \cdot n}{>= pi.max * n} of the
#' \eqn{n} LOO folds are declared stable, and the set is kept only if it
#' contains at least \eqn{2 \times}{2 *} \code{ncomp} features.
#'
#' @param sel_freq Tibble returned by \code{\link{run_stability_selection}}.
#' @param p Integer. Total number of ASV predictor columns in the original
#'   data. This is needed for the false-positive bound calculation. Compute
#'   as \code{sum(startsWith(names(data), "asv"))}.
#' @param n_obs Integer. Number of LOO folds. If \code{NULL} (the default),
#'   read from the \code{"n_obs"} attribute of \code{sel_freq} (set
#'   automatically by \code{\link{run_stability_selection}}).
#' @param pi.max Numeric in (0.5, 1]. Selection frequency threshold. An ASV
#'   must be selected in at least this fraction of LOO folds within the eta
#'   interval to be declared stable. Default 0.8.
#' @param EV.max Numeric > 0. Maximum expected number of false positives.
#'   Lower values give tighter false-positive control and smaller stable sets.
#'   Default 1.
#' @param q.min Integer. Minimum mean selection count per eta interval; avoids
#'   trivially sparse regions where too few features are selected to form a
#'   useful model. Default 15.
#'
#' @return A tibble of candidate stable sets with columns:
#'   \describe{
#'     \item{ncomp}{Integer. Number of PLS components.}
#'     \item{eta.min}{Numeric. Lower bound of the eta interval.}
#'     \item{eta.max}{Numeric. Upper bound of the eta interval.}
#'     \item{stable_asvs}{List of character vectors. ASV identifiers
#'       declared stable in this configuration.}
#'     \item{n_stable}{Integer. Number of stable ASVs.}
#'   }
#'   Rows are sorted by the filtering order. Use nested cross-validation
#'   to select among candidates; the vignette demonstrates a conservative
#'   heuristic (narrowest eta interval).
#'
#' @seealso \code{\link{run_stability_selection}} for the preceding step,
#'   \code{\link{fit_pls_stable}} for fitting a model to the chosen stable set.
#'
#' @references
#' Meinshausen, N. & Bühlmann, P. (2010). Stability selection.
#' \emph{Journal of the Royal Statistical Society: Series B}, 72(4), 417--473.
#' \doi{10.1111/j.1467-9868.2010.00740.x}
#'
#' @examples
#' \dontrun{
#' library(satterthwaite2026)
#' library(dplyr)
#'
#' whales <- inner_join(example_density, example_16s, by = "cruise")
#' p <- sum(startsWith(names(whales), "asv"))
#'
#' # Using pre-computed selection frequencies bundled with the package
#' candidates <- extract_stable_sets(sel_freq_bm, p = p)
#' candidates |> select(ncomp, eta.min, eta.max, n_stable)
#' }
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
