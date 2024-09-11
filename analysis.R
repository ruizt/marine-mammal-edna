## DATA PROCESSING -------------------------------------------------------------
source('scripts/processing/processing-18sv9.R')
source('scripts/processing/processing-18sv4.R')
source('scripts/processing/processing-16s.R')
source('scripts/processing/processing-mm.R')

## GENERATE DATA PARTITIONS FOR VALIDATION PROCEDURES --------------------------
source('scripts/processing/validation-partitions.R')

## SCALED SIGHTINGS ANALYSIS ---------------------------------------------------

# nested validation procedure to choose optimal hyperparameter settings
source('scripts/analysis-scaled-sightings/nested-validation-18sv9-ss.R')
source('scripts/analysis-scaled-sightings/nested-validation-18sv4-ss.R')
source('scripts/analysis-scaled-sightings/nested-validation-16s-ss.R')

# perform stability selection to identify asvs of interest
source('scripts/analysis-scaled-sightings/stability-selection-18sv9-ss.R')
source('scripts/analysis-scaled-sightings/stability-selection-18sv4-ss.R')
source('scripts/analysis-scaled-sightings/stability-selection-16s-ss.R')

# fit models on stable sets
source('scripts/analysis-scaled-sightings/model-fitting-18sv9-ss.R')
source('scripts/analysis-scaled-sightings/model-fitting-18sv4-ss.R')
source('scripts/analysis-scaled-sightings/model-fitting-16s-ss.R')
