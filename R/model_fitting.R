# Model fitting functions
#
# Functions for fitting PLS models to stable feature sets.
# See vignettes/single-iteration.Rmd for a narrative walkthrough.

#' Fit a PLS model to a stable feature set
#'
#' Subsets the data to the stable ASV features and fits a PLS regression model.
#'
#' @param data Data frame with a response variable and ASV predictor columns.
#' @param response Character. Response variable name.
#' @param stable_asvs Character vector of ASV column names, typically
#'   \code{candidates$stable_asvs[[i]]} from \code{\link{extract_stable_sets}}.
#' @param ncomp Integer. Number of PLS components.
#' @return A \code{pls} model object.
#' @export
fit_pls_stable <- function(data, response, stable_asvs, ncomp) {
  df <- data |>
    dplyr::rename(y = !!response) |>
    dplyr::select(y, dplyr::all_of(stable_asvs))
  pls::plsr(y ~ ., data = df, ncomp = ncomp, scale = FALSE, center = TRUE)
}
