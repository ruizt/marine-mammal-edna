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
qmax <- sqrt((2*pimax - 1)*p*EV) 


metrics |> 
  select(species) |> 
  distinct()

sel_freq |> 
  select(species) |>  
  distinct()

## FILTER EXPLORATION - SPECIES = BP ---------------------------------------------

# grid search: ranges of ncomp, eta -> average # of asvs selected
ncomp_start_grid <- 1:9
ncomp_stop_grid <- 10:1
eta_min_grid <- seq(0,0.9, by = 0.1)
eta_max_grid <- seq(1,0.1, by = -0.1)

num_asvs_eta_ncomp_gs <- function(ncomp_start,ncomp_end, eta_min, eta_max){
  avg_asvs <- metrics |> 
    filter(species == "bp",
           ncomp >= ncomp_start,
           ncomp <= ncomp_end,
           eta >= eta_min,
           eta <= eta_max) |> 
    summarise(avg.n.asv = mean(n.asv))
  
  avg_asvs
}


gs_df <- data.frame(0,0,0,0,0)
names(gs_df) = c("ncomp.min", "ncomp.max", "eta.min","eta.max", "avg.n.csv")
for (ns in ncomp_start_grid){
  for (nf in ncomp_stop_grid){
    if (ns > nf){
      next
    }
    for (emin in eta_min_grid){
      for (emax in eta_max_grid){
        if (emin > emax){
          next
        }
        
        avg_asvs = num_asvs_eta_ncomp_gs(ns,nf,emin,emax)
        
        newrow = c(ns,nf,emin,emax,avg_asvs)
        
        print(newrow)
        
        gs_df <- rbind(gs_df, setNames(newrow, names(gs_df)))
        
      }
    }
  }
}

gs_df <- gs_df |> 
  slice(-1)


# Candidate ranges
candidate_ranges <- gs_df |> 
  filter(avg.n.csv <= qmax,
         avg.n.csv > 0)

candidate_ranges <- candidate_ranges |> 
  mutate(ncomp.range = ncomp.max - ncomp.min,
         eta.range = eta.max - eta.min)

candidate_ranges |> 
  filter(ncomp.range == max(ncomp.range) | eta.range == max(eta.range))


## STABLE SET -------------------------------------------------------------------
## INITIAL FILTERING: 1 <= ncomp <= 10 ; 0.4 <= eta <= 1
## => eta >= 0.4

# this is the widest range of ncomp and eta values that still obtains an 
# average number of asvs <= qmax (98.711)

# asvs that have a maximum selection prob of at least 0.8 (max n >= 20)

sel_freq |> 
  filter(species == "bp",
         eta >= 0.4,
         n >= 20) |> 
  select(asv) |> 
  distinct()
  
# asvs that have average selection prob of at least 0.8

sel_freq |> 
  filter(species == "bp",
         eta >= 0.4) |> 
  mutate(prop = n/25) |> 
  group_by(asv) |> 
  summarise(avg.prop = mean(prop)) |> 
  filter(avg.prop >= 0.8)


# asvs that have average selection prob of at least 0.5

sel_freq |> 
  filter(species == "bp",
         eta >= 0.4) |> 
  mutate(prop = n/25) |> 
  group_by(asv) |> 
  summarise(avg.prop = mean(prop)) |> 
  filter(avg.prop >= 0.5)
