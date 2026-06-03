# Dataset documentation
#
# Primary example datasets and supporting objects bundled with the package
# for use in vignettes and interactive exploration.

#' Example 16S eDNA data (500-ASV subset)
#'
#' A 500-ASV subset of the full 16S rRNA eDNA dataset used in Satterthwaite
#' et al. Every ASV identified as stable in the full analysis is present.
#'
#' @format A tibble with 27 rows and 501 columns:
#' \describe{
#'   \item{cruise}{Character. CalCOFI cruise identifier.}
#'   \item{asv.*}{Numeric. Relative abundance of each ASV.}
#' }
#' @source Satterthwaite et al. (forthcoming). Processed from raw 16S
#'   sequencing data; full ASV table archived on Zenodo
#'   \doi{10.5281/zenodo.19139338}.
"example_16s"

#' Example whale density estimates
#'
#' Log-ratio density estimates for blue whale (\emph{Balaenoptera musculus},
#' \code{bm}), fin whale (\emph{B. physalus}, \code{bp}), and humpback whale
#' (\emph{Megaptera novaeangliae}, \code{mn}) across 25 CalCOFI cruises.
#' Each value is \eqn{\log(\hat{D} / \bar{D}_\text{season})}, where
#' \eqn{\hat{D}} is the LOO-cross-validated line-transect density estimate
#' and \eqn{\bar{D}_\text{season}} is the seasonal mean.
#'
#' @format A tibble with 25 rows and 4 columns:
#' \describe{
#'   \item{cruise}{Character. CalCOFI cruise identifier.}
#'   \item{bm}{Numeric. Log-ratio density for blue whale.}
#'   \item{bp}{Numeric. Log-ratio density for fin whale.}
#'   \item{mn}{Numeric. Log-ratio density for humpback whale.}
#' }
#' @source Derived from CalCOFI line-transect surveys via distance sampling;
#'   see \code{vignette("line-transect")}.
"example_density"

#' Pre-computed sPLS selection frequencies for blue whale (16S, example data)
#'
#' Output of \code{\link{run_stability_selection}} applied to the blue whale
#' (\code{bm}) response using \code{example_16s} and \code{example_density}.
#' Bundled to avoid the ~3-minute computation in the vignette.
#'
#' @format A tibble with columns \code{eta} (sparsity), \code{ncomp}
#'   (components), \code{asv} (ASV identifier), and \code{n} (selection count
#'   across LOO partitions). Attribute \code{n_obs} records the number of
#'   LOO folds (25).
#' @seealso \code{\link{run_stability_selection}}, \code{vignette("single-iteration")}
"sel_freq_bm"

#' Pre-computed candidate stable sets for blue whale (16S, example data)
#'
#' Output of \code{\link{extract_stable_sets}} applied to \code{sel_freq_bm}.
#'
#' @format A tibble with columns \code{ncomp}, \code{eta.min}, \code{eta.max},
#'   \code{stable_asvs} (list of character vectors), and \code{n_stable}.
#' @seealso \code{\link{extract_stable_sets}}, \code{vignette("single-iteration")}
"candidate_sets_bm"

#' Stable ASVs from the full 16S analysis
#'
#' Character vector of the 93 ASV identifiers judged stable for at least one
#' whale species in the full 16S analysis. Used in \code{vignette("single-iteration")}
#' to verify that features selected on the 500-ASV example subset are genuine.
#'
#' @format A character vector of length 93.
"all_stable"

#' Example ASVs
#'
#' Character vector of the 500 ASV identifiers included in \code{example_16s}.
#'
#' @format A character vector of length 500.
"example_asvs"

#' Seasonal mean log-densities
#'
#' Seasonal mean log-density for each whale species, used to back-transform
#' log-ratio fitted values to the original density scale.
#'
#' @format A tibble with 4 rows (seasons) and 4 columns:
#' \describe{
#'   \item{season}{Factor. Season (winter, spring, summer, fall).}
#'   \item{log.bm.imp.mean}{Numeric. Mean log-density for blue whale.}
#'   \item{log.bp.imp.mean}{Numeric. Mean log-density for fin whale.}
#'   \item{log.mn.imp.mean}{Numeric. Mean log-density for humpback whale.}
#' }
"dens_means"

#' Raw line-transect density estimates
#'
#' Cruise-level density estimates (individuals km\eqn{^{-2}}) for each whale
#' species from distance-sampling line-transect surveys.
#'
#' @format A tibble with 75 rows (25 cruises × 3 species) and 11 columns
#'   including \code{cruise}, \code{species} (\code{"bm"}, \code{"bp"},
#'   \code{"mn"}), \code{estimate}, standard error, CV, confidence limits,
#'   degrees of freedom, and \code{season}.
#' @source CalCOFI line-transect surveys; see \code{vignette("line-transect")}.
"dens_raw"
