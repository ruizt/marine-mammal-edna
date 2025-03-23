## DATA PROCESSING -------------------------------------------------------------
source('scripts/processing/processing-18sv9.R')
source('scripts/processing/processing-18sv4.R')
source('scripts/processing/processing-16s.R')
source('scripts/processing/processing-density.R')
source('scripts/processing/validation-partitions-dens.R')

## DENSITY MODELING ------------------------------------------------------------

# perform stability selection to identify asvs of interest
source('scripts/analysis-density-estimates/stability-selection-18sv9-dens.R')
source('scripts/analysis-density-estimates/stability-selection-18sv4-dens.R')
source('scripts/analysis-density-estimates/stability-selection-16s-dens.R')

# nested validation procedure to assess consistency of selection procedure
source('scripts/analysis-density-estimates/nested-validation-18sv9-dens.R')
source('scripts/analysis-density-estimates/nested-validation-18sv4-dens.R')
source('scripts/analysis-density-estimates/nested-validation-16s-ss.R')

# model fitting to stable sets
source('scripts/analysis-density-estimates/model-fitting-16s-dens.R')
source('scripts/analysis-density-estimates/model-fitting-18sv4-dens.R')
source('scripts/analysis-density-estimates/model-fitting-18sv9-dens.R')