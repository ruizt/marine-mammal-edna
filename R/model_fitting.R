# Model fitting functions
#
# Functions for fitting PLS models to stable feature sets.
# See vignettes/single-iteration.qmd for a narrative walkthrough.

#' Fit a PLS model restricted to a stable feature set
#'
#' Subsets the data to the specified stable ASV features and fits a partial
#' least squares (PLS) regression model via \code{\link[pls]{plsr}}.
#'
#' @details
#' This function is the final step of the stability selection pipeline: once a
#' stable feature set has been identified by
#' \code{\link{extract_stable_sets}}, this function fits the actual predictive
#' model using only those features.
#'
#' The model is fit with \code{scale = FALSE} and \code{center = TRUE}
#' (predictors are mean-centered but not scaled to unit variance). The input
#' ASV abundances should already be on the CLR scale, and the response should
#' be a log-ratio density.
#'
#' The returned object is a standard \code{\link[pls]{mvr}} model, so all
#' \pkg{pls} methods apply:
#' \itemize{
#'   \item \code{\link[stats]{coef}(fit)[, , "K comps"]} extracts
#'     regression coefficients for \code{K} components
#'   \item \code{\link[stats]{fitted}(fit)[, , "K comps"]} extracts
#'     fitted values
#'   \item \code{\link[stats]{predict}(fit, newdata)} produces predictions
#'     on new data
#' }
#'
#' @param data Data frame containing the response variable and ASV predictor
#'   columns. Must include all columns named in \code{stable_asvs}.
#' @param response Character string. Name of the response column.
#' @param stable_asvs Character vector of ASV column names to include as
#'   predictors. Typically obtained from
#'   \code{candidates$stable_asvs[[i]]} via \code{\link{extract_stable_sets}}.
#' @param ncomp Integer. Number of PLS components to fit. Should match the
#'   \code{ncomp} value from the selected candidate row in
#'   \code{\link{extract_stable_sets}}.
#'
#' @return An object of class \code{\link[pls]{mvr}} (PLS regression fit).
#'
#' @seealso \code{\link{extract_stable_sets}} for obtaining the stable feature
#'   set, \code{\link[pls]{plsr}} for the underlying PLS implementation,
#'   \code{\link[pls]{mvr}} for methods available on the returned object.
#'
#' @examples
#' \dontrun{
#' library(satterthwaite2026)
#' library(dplyr)
#'
#' whales <- inner_join(example_density, example_16s, by = "cruise")
#' p <- sum(startsWith(names(whales), "asv"))
#' candidates <- extract_stable_sets(sel_freq_bm, p = p)
#'
#' # Take the most conservative candidate
#' best <- candidates |>
#'   mutate(interval_width = eta.max - eta.min) |>
#'   slice_min(interval_width, with_ties = FALSE)
#'
#' fit <- fit_pls_stable(
#'   whales, response = "bm",
#'   stable_asvs = best$stable_asvs[[1]],
#'   ncomp = best$ncomp
#' )
#'
#' # Regression coefficients
#' coef(fit)[, , paste(best$ncomp, "comps")]
#'
#' # Fitted values
#' fitted(fit)[, , paste(best$ncomp, "comps")]
#' }
#' @export
fit_pls_stable <- function(data, response, stable_asvs, ncomp) {
  df <- data |>
    dplyr::rename(y = !!response) |>
    dplyr::select(y, dplyr::all_of(stable_asvs))
  pls::plsr(y ~ ., data = df, ncomp = ncomp, scale = FALSE, center = TRUE)
}
