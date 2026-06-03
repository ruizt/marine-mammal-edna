# Utility functions
#
# Helper functions for setting up the stability selection pipeline.

#' Build an eta hyperparameter grid on a log scale
#'
#' Generates a sequence of sparsity (\code{eta}) values for use with
#' \code{\link{run_stability_selection}}. The values are spaced on a log
#' scale and transformed so that the resulting grid is denser near the upper
#' end (higher sparsity), where small changes in eta have a larger effect on
#' the number of selected features.
#'
#' @details
#' The grid is constructed as
#' \code{rev(1 - exp(seq(log(lo), log(hi), length = n)))}, which produces
#' \code{n} values ranging from approximately \code{1 - hi} to
#' \code{1 - lo}. With the defaults (\code{lo = 0.075}, \code{hi = 0.6}),
#' this gives a grid spanning roughly 0.4 to 0.925.
#'
#' In the full analysis (Satterthwaite et al.), all three markers (16S, 18S V4,
#' 18S V9) used the same grid: \code{eta_grid(n = 55, lo = 0.075, hi = 0.6)}
#' (see \code{analysis/2a-stability-selection-16s.R}, lines 32--33).
#'
#' The \code{eta} parameter controls the L1/L2 penalty balance in
#' \code{\link[spls]{spls}}: \code{eta = 1} is pure L1 (LASSO-like, maximal
#' sparsity), while lower values blend in an L2 (ridge) penalty and select
#' more features. The range 0.4--0.925 covers the regime from moderately
#' sparse to highly sparse models.
#'
#' @param n Integer. Number of grid points. Default 55, matching the paper.
#' @param lo Numeric. Lower bound of the log-scale sequence (maps to the
#'   upper end of the eta grid). Default 0.075.
#' @param hi Numeric. Upper bound of the log-scale sequence (maps to the
#'   lower end of the eta grid). Default 0.6.
#'
#' @return Numeric vector of length \code{n}, sorted in increasing order.
#'
#' @seealso \code{\link{run_stability_selection}}, \code{\link[spls]{spls}}
#'
#' @examples
#' vals <- eta_grid()
#' length(vals)        # 55
#' range(vals)         # approximately 0.40 to 0.93
#'
#' # Coarser grid for faster prototyping
#' eta_grid(n = 15)
#'
#' @export
eta_grid <- function(n = 55, lo = 0.075, hi = 0.6) {
  rev(1 - exp(seq(log(lo), log(hi), length = n)))
}

#' Create leave-one-out partitions
#'
#' Builds a tibble of leave-one-out (LOO) train/test splits from a data frame
#' with a \code{cruise} column. Each row holds one partition: the test set is
#' a single cruise, and the training set is all remaining cruises.
#'
#' @details
#' The number of partitions equals the number of rows in \code{data} (one per
#' cruise). Each training set has \eqn{n - 1} rows, and each test set has 1
#' row.
#'
#' The resulting tibble is used as the \code{partitions} argument to
#' \code{\link{run_stability_selection}}. It can also be used for manual LOO
#' cross-validation when evaluating candidate stable sets.
#'
#' @param data Data frame. Must contain a \code{cruise} column with unique
#'   identifiers for each observation.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{test.id}{Character. Cruise identifier for the held-out
#'       observation.}
#'     \item{train}{List of data frames. Each element is the training set
#'       (all rows except the held-out cruise).}
#'     \item{test}{List of data frames. Each element is a single-row data
#'       frame for the held-out cruise.}
#'   }
#'
#' @seealso \code{\link{run_stability_selection}}
#'
#' @examples
#' \dontrun{
#' library(satterthwaite2026)
#' library(dplyr)
#'
#' whales <- inner_join(example_density, example_16s, by = "cruise")
#' parts <- loo_partitions(whales)
#' nrow(parts)             # 25 partitions
#' nrow(parts$train[[1]])  # 24 training observations
#' }
#' @export
loo_partitions <- function(data) {
  tibble::tibble(
    test.id = data$cruise,
    train   = purrr::map(data$cruise, ~dplyr::filter(data, cruise != .x)),
    test    = purrr::map(data$cruise, ~dplyr::filter(data, cruise == .x))
  )
}
