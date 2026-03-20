## LINE TRANSECT AND DATA PROCESSING -------------------------------------------
#
# Produces density-estimates.RData and processed eDNA/sightings data.
# Processing scripts (0b-0g) require unpublished raw data; their outputs are
# available from https://doi.org/10.5281/zenodo.19139338 via 1-setup-data.R.

source('analysis/0_line-transect/0a-line-transect.R')
source('analysis/0_processing/0b-processing-16s.R')
source('analysis/0_processing/0c-processing-18sv4.R')
source('analysis/0_processing/0d-processing-18sv9.R')
source('analysis/0_processing/0e-processing-density.R')
source('analysis/0_processing/0f-processing-sightings.R')
source('analysis/0_processing/0g-validation-partitions.R')

## STATISTICAL MODELING --------------------------------------------------------

# perform stability selection to identify asvs of interest
source('analysis/2a-stability-selection-16s.R')
source('analysis/2b-stability-selection-18sv4.R')
source('analysis/2c-stability-selection-18sv9.R')

# nested validation procedure; requires nested partition files from zenodo
source('analysis/3a-nested-validation-16s.R')
source('analysis/3b-nested-validation-18sv4.R')
source('analysis/3c-nested-validation-18sv9.R')

# model fitting to stable sets
source('analysis/4a-model-fitting-16s.R')
source('analysis/4b-model-fitting-18sv4.R')
source('analysis/4c-model-fitting-18sv9.R')

## NAIVE PREDICTIONS -----------------------------------------------------------

source('analysis/5-naive-preds.R')

## RESULTS ---------------------------------------------------------------------

source('analysis/6a-figures.R')
source('analysis/6b-tables.R')
source('analysis/6c-data.R')
