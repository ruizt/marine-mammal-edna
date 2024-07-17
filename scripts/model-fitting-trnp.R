library(tidyverse)
library(spls)

# read in selection frequencies from LOOCV
sel_freq <- read_rds('rslt/loocv/2024-07-15/selection-frequencies.rds')

# read in metrics from LOOCV
metrics <- read_rds('rslt/loocv/2024-07-15/metrics.rds')

## SPECIFYING NCOMP/ETA RANGE --------------------------------------------------
# chosen stability threshold (minimum max selection prob.)
pimax <- 0.8

# desired bound for expected number of false positives
EV <- 5

# total number of candidate ASVs
p <- 3248

# threshold for max average number of selected ASVs
qmax <- sqrt((2*pi_max - 1)*p*EV) 
