## DATA PROCESSING -------------------------------------------------------------

# requires unpublished raw data; outputs stored in ~/data/processed/
source('scripts/processing/processing-18sv9.R')
source('scripts/processing/processing-18sv4.R')
source('scripts/processing/processing-16s.R')
source('scripts/processing/processing-density.R')
source('scripts/processing/processing-sightings.R')
source('scripts/processing/validation-partitions.R')

## DENSITY MODELING ------------------------------------------------------------

# perform stability selection to identify asvs of interest
source('scripts/analysis/stability-selection-18sv9.R')
source('scripts/analysis/stability-selection-18sv4.R')
source('scripts/analysis/stability-selection-16s.R')

# nested validation procedure; requires nested partition files not on GH
source('scripts/analysis/nested-validation-18sv9.R')
source('scripts/analysis/nested-validation-18sv4.R')
source('scripts/analysis/nested-validation-16s.R')

# model fitting to stable sets
source('scripts/analysis/model-fitting-16s.R')
source('scripts/analysis/model-fitting-18sv4.R')
source('scripts/analysis/model-fitting-18sv9.R')

## NAIVE PREDICTIONS -----------------------------------------------------------

source('scripts/analysis/naive-preds.R')

## RESULTS ---------------------------------------------------------------------

source('scripts/figures.R')
source('scripts/tables.R')
source('scripts/data.R')