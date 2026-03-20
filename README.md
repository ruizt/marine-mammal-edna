# predicting marine mammal density from eDNA

Code repository for: *Microbial and small zooplankton communities predict density of baleen whales in the southern California Current Ecosystem.*

Citation:

> E.V. Satterthwaite, T.D. Ruiz, N.V. Patin, M.N. Alksne, L. Thomas, J. Dinasquet, R.H. Lampe, K.G. Chan, N.A. Patrick, A.E. Allen, S. Baumann-Pickering, B.X. Semmens. Microbial and small zooplankton communities predict density of baleen whales in the southern California Current Ecosystem.

Repository contributors: T.D. Ruiz, N.A. Patrick, K.G. Chan, L. Thomas

---

## Quick start

### For methods adopters

The `vignettes/` directory contains self-contained walkthroughs of the analysis pipeline:

- **`vignettes/single-iteration.Rmd`** — demonstrates the full pipeline (stability selection → model fitting → interpretation) for one eDNA marker (16S) and one whale species (blue whale) using a bundled 500-ASV example dataset. No data download required.
- **`vignettes/full-analysis.Rmd`** — covers all three markers × three species combinations from the paper. Requires processed data (see below).
- **`vignettes/line-transect.qmd`** — documents the distance-sampling analysis used to estimate whale density from the CalCOFI survey data.

### For paper reproducers

1. Download processed data from Zenodo (v1.1.0, [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19139338.svg)](https://doi.org/10.5281/zenodo.19139338)):

```r
source('analysis/setup-data.R')               # ~16 MB core files
# source('analysis/setup-data.R'); setup_data(partitions = TRUE)  # +5.5 GB partitions
```

2. Run the full analysis pipeline (note: stability selection and nested validation take approximately 3 hours per marker):

```r
source('analysis/run-all.R')
```

Raw eDNA sequencing data (ASV count tables) are not included in this repository or the Zenodo archive but may be requested from the authors.

---

## Repository structure

```
├── R/                        # Package functions
│   ├── stability_selection.R # sPLS stability selection
│   ├── model_fitting.R       # PLS model fitting and metrics
│   └── utils.R               # Shared utilities (eta grid, LOO partitions)
│
├── analysis/                 # Executable analysis scripts
│   ├── run-all.R             # Top-level orchestration
│   ├── setup-data.R          # Download data from Zenodo
│   ├── processing/           # Raw data processing (requires unpublished raw data)
│   ├── stability-selection-*.R
│   ├── nested-validation-*.R
│   ├── model-fitting-*.R
│   ├── naive-preds.R
│   ├── line-transect.R       # Line transect density estimation
│   ├── figures.R
│   ├── tables.R
│   └── create-example-data.R # Backend: generates bundled example dataset
│
├── vignettes/                # Narrative walkthroughs
│   ├── single-iteration.Rmd  # Short: one marker, one species
│   ├── full-analysis.Rmd     # Long: all 3 markers × 3 species
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

---

## Data

Processed data are archived on Zenodo (v1.1.0): [10.5281/zenodo.19139338](https://doi.org/10.5281/zenodo.19139338).

| File | Description |
|------|-------------|
| `ncog16s.RData` | Processed 16S eDNA ASV table |
| `ncog18sv4.RData` | Processed 18S V4 eDNA ASV table |
| `ncog18sv9.RData` | Processed 18S V9 eDNA ASV table |
| `density-estimates.RData` | Line-transect whale density estimates |
| `sightings.RData` | CalCOFI whale sighting records |
| `sample_table.Rds` | Survey effort data (line transect analysis) |
| `region_table.Rds` | Cruise/region metadata (line transect analysis) |
| `obs_table.Rds` | Individual sighting records (line transect analysis) |
| `_partitions.zip` | LOO cross-validation partitions (~3.1 GB) |
| `_combined-partitions.zip` | Combined eDNA + density partitions (~2.4 GB) |

Map figures use coastline shapefiles from GSHHG Release v2.3.7 ([https://www.soest.hawaii.edu/pwessel/gshhg/](https://www.soest.hawaii.edu/pwessel/gshhg/)):
> Wessel, P., and W. H. F. Smith (1996), A global, self-consistent, hierarchical, high-resolution shoreline database, J. Geophys. Res., 101(B4), 8741–8743, doi:10.1029/96JB00104.

---

## Software

R version 4.4.3 (2025-02-28) · RStudio 2024.12.1+563 · macOS Sequoia 15.3.2

Dependencies are listed in `DESCRIPTION`. Key packages: `spls`, `pls`, `tidyverse`, `collapse`, `sf`, `patchwork`, `Distance`.
