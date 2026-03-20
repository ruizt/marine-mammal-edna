# Utility functions

#' Build an eta hyperparameter grid on a log scale
#'
#' @param n Integer. Number of grid points (default 55, matching the paper).
#' @param lo Numeric. Lower bound (default 0.075).
#' @param hi Numeric. Upper bound (default 0.6).
#' @return Numeric vector of eta values in increasing order.
#' @export
eta_grid <- function(n = 55, lo = 0.075, hi = 0.6) {
  rev(1 - exp(seq(log(lo), log(hi), length = n)))
}

#' Create leave-one-out partitions from a joined data frame
#'
#' @param data Data frame with a cruise column.
#' @return Tibble with columns test.id, train (list-col), test (list-col).
#' @export
loo_partitions <- function(data) {
  tibble::tibble(
    test.id = data$cruise,
    train   = purrr::map(data$cruise, ~dplyr::filter(data, cruise != .x)),
    test    = purrr::map(data$cruise, ~dplyr::filter(data, cruise == .x))
  )
}
