## DATA PROCESSING -------------------------------------------------------------
#
# Requires unpublished raw data; outputs stored in data/
# Processed outputs are available from https://doi.org/10.5281/zenodo.19139338
# and can be downloaded automatically via analysis/0-setup-data.R

source('analysis/processing/1a-processing-16s.R')
source('analysis/processing/1b-processing-18sv4.R')
source('analysis/processing/1c-processing-18sv9.R')
source('analysis/processing/1d-processing-density.R')
source('analysis/processing/1e-processing-sightings.R')
source('analysis/processing/1f-validation-partitions.R')

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

## LINE TRANSECT ANALYSIS ------------------------------------------------------

source('analysis/6-line-transect.R')

## RESULTS ---------------------------------------------------------------------

source('analysis/7a-figures.R')
source('analysis/7b-tables.R')
source('analysis/7c-data.R')
