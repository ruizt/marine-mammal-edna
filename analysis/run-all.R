## DATA PROCESSING -------------------------------------------------------------
#
# Requires unpublished raw data; outputs stored in data/
# Processed outputs are available from https://doi.org/10.5281/zenodo.19139338
# and can be downloaded automatically via analysis/setup-data.R

source('analysis/processing/processing-18sv9.R')
source('analysis/processing/processing-18sv4.R')
source('analysis/processing/processing-16s.R')
source('analysis/processing/processing-density.R')
source('analysis/processing/processing-sightings.R')
source('analysis/processing/validation-partitions.R')

## STATISTICAL MODELING --------------------------------------------------------

# perform stability selection to identify asvs of interest
source('analysis/stability-selection-18sv9.R')
source('analysis/stability-selection-18sv4.R')
source('analysis/stability-selection-16s.R')

# nested validation procedure; requires nested partition files from zenodo
source('analysis/nested-validation-18sv9.R')
source('analysis/nested-validation-18sv4.R')
source('analysis/nested-validation-16s.R')

# model fitting to stable sets
source('analysis/model-fitting-16s.R')
source('analysis/model-fitting-18sv4.R')
source('analysis/model-fitting-18sv9.R')

## NAIVE PREDICTIONS -----------------------------------------------------------

source('analysis/naive-preds.R')

## LINE TRANSECT ANALYSIS ------------------------------------------------------

source('analysis/line-transect.R')

## RESULTS ---------------------------------------------------------------------

source('analysis/figures.R')
source('analysis/tables.R')
source('analysis/data.R')
