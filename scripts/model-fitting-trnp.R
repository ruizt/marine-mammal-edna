library(tidyverse)
library(spls)

# read in selection frequencies from LOOCV
sel_freq <- read_rds('rslt/loocv/2024-07-15/selection-frequencies.rds')

# read in metrics from LOOCV
metrics <- read_rds('rslt/loocv/2024-07-15/metrics.rds')
