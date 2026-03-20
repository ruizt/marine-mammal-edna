# Model fitting functions
#
# Functions for fitting PLS models to stable feature sets and computing
# model metrics. See analysis/model-fitting-*.R and vignettes/single-iteration.Rmd.

#' Subset joined data to a stable feature set and response variable
#'
#' @param data Data frame. Full joined eDNA + density data.
#' @param stable_set Character vector of ASV column names.
#' @param species Character. Response variable name (e.g. "bm").
#' @return Data frame with columns y (response) and the stable ASV features.
#' @export
subset_to_stable <- function(data, stable_set, species) {
  data |>
    dplyr::rename(y = !!species) |>
    dplyr::select(y, dplyr::all_of(stable_set))
}

#' Fit a PLS regression model
#'
#' @param data Data frame with column y and predictor columns.
#' @param ncomp Integer. Number of PLS components.
#' @return A \code{pls} model object.
#' @export
fit_pls <- function(data, ncomp) {
  pls::plsr(y ~ ., data = data, ncomp = ncomp, scale = FALSE, center = TRUE)
}

#' Compute model metrics for fitted PLS models
#'
#' @param fit_df Tibble of observed and fitted values with columns
#'   species, lr.obs, lr.fit, dens.obs, dens.fit, lr.resid, dens.resid.
#' @param n_obs Integer. Number of observations.
#' @param stable_sets Tibble with columns species, ncomp, ss (list-col).
#' @return Tibble with adjusted R² on log-ratio and density scales per species.
#' @export
compute_model_metrics <- function(fit_df, n_obs, stable_sets) {
  fit_df |>
    dplyr::group_by(species) |>
    dplyr::summarize(
      adj.rsq.lr   = 1 - ((n_obs - 1) / (n_obs - unique(ncomp) - 1)) *
        stats::var(lr.resid) / stats::var(lr.obs),
      adj.rsq.dens = 1 - ((n_obs - 1) / (n_obs - unique(ncomp) - 1)) *
        stats::var(dens.resid) / stats::var(dens.obs),
      .groups = 'drop'
    ) |>
    dplyr::left_join(stable_sets, dplyr::join_by(species)) |>
    dplyr::mutate(n.asv = purrr::map_int(ss, length)) |>
    dplyr::select(species, n.asv, dplyr::ends_with('.lr'), dplyr::ends_with('.dens'))
}
