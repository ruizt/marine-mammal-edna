# satterthwaite2026

Code repository for: *Microbial and small zooplankton communities predict density of baleen whales in the southern California Current Ecosystem.*

Citation:

> E.V. Satterthwaite, T.D. Ruiz, N.V. Patin, M.N. Alksne, L. Thomas, J. Dinasquet, R.H. Lampe, K.G. Chan, N.A. Patrick, A.E. Allen, S. Baumann-Pickering, B.X. Semmens. Microbial and small zooplankton communities predict density of baleen whales in the southern California Current Ecosystem. *Forthcoming in PLOS One.*

Repository contributors: T.D. Ruiz, N.A. Patrick, K.G. Chan, L. Thomas

------------------------------------------------------------------------

## About

This repository is structured as an R package but also, separately from package data and functions, contains scripts to reproduce the results of the paper. Thus it serves two distinct purposes.

-   Readers that wish to **implement the methodology** in a different context should install the package.

-   Readers that wish to **reproduce the paper results** should clone the repository.

Further instructions for each are below.

### For methods adopters

To use the methodology, install the package.

``` r
# install.packages("devtools")
devtools::install_github("ruizt/marine-mammal-edna")
```

Package functions provide a three-step pipeline for implementing stability selection adapted to sPLS: fit models on subsamples; extract stable sets; refit on the full data. An example is shown below.

``` r
library(satterthwaite2026)

# 1. Fit sPLS across the hyperparameter grid and LOO subsamples;
#    returns per-ASV selection frequencies
sel_freq <- run_stability_selection(data, response = "bm",
                                    eta_vals = eta_grid(),
                                    ncomp_vals = 4:12,
                                    partitions = loo_partitions(data))

# 2. Apply the Meinshausen-Bühlmann criterion to identify stable features
candidates <- extract_stable_sets(sel_freq, p = ncol_asvs)

# 3. Fit PLS restricted to the stable feature set
fit <- fit_pls_stable(data, response = "bm",
                      stable_asvs = candidates$stable_asvs[[1]],
                      ncomp = candidates$ncomp[[1]])
```

The [**single-iteration walkthrough**](https://htmlpreview.github.io/?https://github.com/ruizt/marine-mammal-edna/blob/main/doc/single-iteration.html) (`vignettes/single-iteration.Rmd`) demonstrates use for one eDNA marker (16S) and one whale species (blue whale) using a bundled 500-ASV example dataset.

### For paper reproducers

The full analysis from the paper implemented in a series of scripts stored in the `analysis` directory. Steps to reproduce the results of the paper in full are:

1.  Clone the repository.

``` bash
git clone https://github.com/ruizt/marine-mammal-edna
```

2.  Download processed data from Zenodo (v1.1.0, [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19139338.svg)](https://doi.org/10.5281/zenodo.19139338)):

``` r
source('analysis/setup-data.R')               # ~16 MB core files
# source('analysis/setup-data.R'); setup_data(partitions = TRUE)  # +5.5 GB partitions
```

3.  Run the full analysis pipeline (note: stability selection and nested validation take approximately 3 hours per marker):

``` r
source('analysis/run-all.R')
```

The **[full analysis vignette](https://htmlpreview.github.io/?https://github.com/ruizt/marine-mammal-edna/blob/main/doc/full-analysis.html)** (`vignettes/full-analysis.Rmd`) narrates each stage of `run-all.R` and shows summary tables of the key results.

Raw eDNA sequencing data (ASV count tables) are not included in this repository or the Zenodo archive but may be requested from the authors.

------------------------------------------------------------------------

## Repository structure

```         
├── R/                        # Package functions
│   ├── stability_selection.R # sPLS stability selection
│   ├── model_fitting.R       # PLS model fitting and metrics
│   └── utils.R               # Shared utilities (eta grid, LOO partitions)
│
├── analysis/                 # Executable analysis scripts
│   ├── run-all.R             # Top-level orchestration
│   ├── 0-setup-data.R        # Download data from Zenodo
│   ├── processing/           # Step 1: raw data processing (requires unpublished raw data)
│   │   └── 1a–1f-*.R
│   ├── 2a–2c-stability-selection-*.R   # Step 2: sPLS stability selection (a/b/c = 16S/18SV4/18SV9)
│   ├── 3a–3c-nested-validation-*.R     # Step 3: nested CV to select stable sets
│   ├── 4a–4c-model-fitting-*.R         # Step 4: PLS model fitting
│   ├── 5-naive-preds.R                 # Step 5: naive baseline predictions
│   ├── 6-line-transect.R               # Step 6: distance-sampling density estimation
│   ├── 7a-figures.R                    # Step 7: manuscript outputs
│   ├── 7b-tables.R
│   ├── 7c-data.R
│   └── create-example-data.R # Developer: generates bundled example dataset
│
├── vignettes/                # Narrative walkthroughs
│   ├── single-iteration.Rmd  # Package API walkthrough (one marker, one species)
│   ├── full-analysis.Rmd     # Paper reproduction guide (narrates run-all.R)
│   └── line-transect.qmd     # Supplement: distance-sampling analysis
│
├── tests/testthat/           # IO and unit tests
│
├── data/                     # Data directory (gitignored except example data)
│   ├── example-data-16s.RData  # Bundled 500-ASV example dataset
│   └── gshhg.RData             # Coastline shapefiles for maps
│
└── rslt/                     # Analysis outputs
    ├── stability-selection/
    ├── nested-validation/
    └── models/
```

------------------------------------------------------------------------

## Data

Processed data are archived on Zenodo (v1.1.0): [10.5281/zenodo.19139338](https://doi.org/10.5281/zenodo.19139338). These can be downloaded via

| File | Description |
|----|----|
| `ncog16s.RData` | Processed 16S eDNA ASV table |
| `ncog18sv4.RData` | Processed 18S V4 eDNA ASV table |
| `ncog18sv9.RData` | Processed 18S V9 eDNA ASV table |
| `density-estimates.RData` | Line-transect whale density estimates |
| `sightings.RData` | CalCOFI whale sighting records |
| `sample_table.Rds` | Survey effort data (line transect analysis) |
| `region_table.Rds` | Cruise/region metadata (line transect analysis) |
| `obs_table.Rds` | Individual sighting records (line transect analysis) |
| `_partitions.zip` | LOO cross-validation partitions (\~3.1 GB) |
| `_combined-partitions.zip` | Combined eDNA + density partitions (\~2.4 GB) |

Map figures use coastline shapefiles from GSHHG Release v2.3.7 (<https://www.soest.hawaii.edu/pwessel/gshhg/>): \> Wessel, P., and W. H. F. Smith (1996), A global, self-consistent, hierarchical, high-resolution shoreline database, J. Geophys. Res., 101(B4), 8741–8743, <doi:10.1029/96JB00104>.

------------------------------------------------------------------------

## Software

R version 4.4.3 (2025-02-28) · RStudio 2024.12.1+563 · macOS Sequoia 15.3.2

Dependencies are listed in `DESCRIPTION`. Key packages: `spls`, `pls`, `tidyverse`, `collapse`, `sf`, `patchwork`, `Distance`.

Package infrastructure assisted by [Posit AI](https://positai.com) (March 2026).
