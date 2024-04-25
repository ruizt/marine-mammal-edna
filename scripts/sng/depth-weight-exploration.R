library(tidyverse)
library(compositions)
library(zCompositions)
library(vegan)
load('data/ncog-18s.RData')

# dplyr::select columns of interest based on nonzero frequency across samples
cols_of_interest <- edna_samples |>
  dplyr::select(starts_with('asv')) |>
  as.matrix() |>
  apply(2, function(.x){mean(.x > 0)}) |>
  t() |>
  as_tibble() |>
  gather(col, prop.nz) |>
  filter(prop.nz > 0.3) |> # adjustable?
  pull(col)

# filter to samples with at least 25% nonzero reads
rows_of_interest <- edna_samples |> 
  dplyr::select(all_of(cols_of_interest)) |>
  mutate(across(where(is.numeric), ~ .x > 0)) |>
  mutate(total = rowMeans(pick(starts_with('asv'))),
         rownum = row_number()) |>
  dplyr::select(total, rownum) |>
  filter(total > 0.25) |>
  pull(rownum)

# zero imputation (bayesian multiplicative, martin-fernandez 2015)
imputation_out <- edna_samples |> 
  dplyr::select(all_of(cols_of_interest)) |>
  slice(rows_of_interest) |>
  zCompositions::cmultRepl(label = 0, 
                           method = 'GBM', 
                           output = 'prop',
                           z.warning = 0.75)

# bind imputed values to sample info
edna_imputed <- edna_samples |> 
  slice(rows_of_interest) |>
  dplyr::select(1:5) |>
  bind_cols(imputation_out)

# aggregate to cruise level
depth_weight_fn <- function(depth, alpha, beta){
  dgamma(depth, shape = alpha, scale = beta) ##CHANGE????
}

##Grid search over shape and scale to maximize alpha diversity
#full list of summary statistics (Mean, median, sd, quartiles, min, max, etc)
alpha.vec <- seq(2,11,1)
beta.vec <- seq(1,15,1)

gs1 <- data.frame(0,0,0,0,0,0,0,0,0)
names(gs1) = c("alpha", "beta", "avgDiv", "SD", "Min", "25", "50","75", "Max")


for (a in alpha.vec){
  for (b in beta.vec){
edna_agg <- edna_imputed |>
  mutate(depthm = as.numeric(depthm),
         weight = depth_weight_fn(depthm, a, b)) |>
  group_by(cruise, line, sta) |>
  summarize(across(starts_with('asv'), 
                   ~weighted.mean(log(.x), weight)), ## automatically normalizes weights
            .groups = 'drop') |>
  group_by(cruise, line) |>
  summarize(across(starts_with('asv'), ~mean(.x))) |>
  group_by(cruise) |>
  summarize(across(starts_with('asv'), ~mean(.x)))

edna_data <- edna_agg |> mutate(exp(pick(-cruise)))

# split into sample info and asvs
asv <- edna_data |> 
  dplyr::select(starts_with('asv'))

sampinfo <- edna_data |> 
  dplyr::select(-starts_with('asv'))

# alpha diversity measures
alpha_divs <- sampinfo %>%
  bind_cols(alpha.div.sh = diversity(asv, index = 'shannon')) 

alpha_div <- alpha_divs |> 
  slice(-14) |> 
  summarise(avgDiv = mean(alpha.div.sh),
            sd = sd(alpha.div.sh),
            min = min(alpha.div.sh),,
            twentyfifth = quantile(alpha.div.sh, .25),
            fifty = quantile(alpha.div.sh, .50),
            seventyfifth = quantile(alpha.div.sh, .75),
            max = max(alpha.div.sh))

newRow = c(a,b, alpha_div$avgDiv, alpha_div$sd, alpha_div$min, alpha_div$twentyfifth, alpha_div$fifty, alpha_div$seventyfifth, alpha_div$max)

print(newRow)

gs1 <- rbind(gs1, setNames(newRow,names(gs1)))
  }
}


gs1 <- gs1 |> 
  slice(-1)


save(gs1, file = "rslt/grid-search-result-wide.RData")


gs1 |> 
  filter(avgDiv == max(gs1$avgDiv))

library(rgl)

with(gs1, plot3d(x = alpha, y = beta, z = SD))

##########################################
## Other Distributions ###################
##########################################

# Normal
# Target certain depth ranges??


mu.vec = c(10,25,50,75,100)
sigma.vec = c(5,10,15)

depth_weight_fn2 <- function(depth, mu, sigma){
  dnorm(depth, mean = mu, sd= sigma) ##CHANGE????
}

normalDF <- data.frame(0,0,0,0,0,0,0,0,0)
names(normalDF) = c("mu", "sigma", "avgDiv", "SD", "Min", "25", "50","75", "Max")

for (a in mu.vec){
  for (b in sigma.vec){
    edna_agg <- edna_imputed |>
      mutate(depthm = as.numeric(depthm),
             weight = depth_weight_fn2(depthm, a, b)) |>
      group_by(cruise, line, sta) |>
      summarize(across(starts_with('asv'), 
                       ~weighted.mean(log(.x), weight)), ## automatically normalizes weights
                .groups = 'drop') |>
      group_by(cruise, line) |>
      summarize(across(starts_with('asv'), ~mean(.x))) |>
      group_by(cruise) |>
      summarize(across(starts_with('asv'), ~mean(.x)))
    
    edna_data <- edna_agg |> mutate(exp(pick(-cruise)))
    
    # split into sample info and asvs
    asv <- edna_data |> 
      dplyr::select(starts_with('asv'))
    
    sampinfo <- edna_data |> 
      dplyr::select(-starts_with('asv'))
    
    # alpha diversity measures
    alpha_divs <- sampinfo %>%
      bind_cols(alpha.div.sh = diversity(asv, index = 'shannon')) 
    
    alpha_div <- alpha_divs |> 
      slice(-14) |> 
      summarise(avgDiv = mean(alpha.div.sh),
                sd = sd(alpha.div.sh),
                min = min(alpha.div.sh),,
                twentyfifth = quantile(alpha.div.sh, .25),
                fifty = quantile(alpha.div.sh, .50),
                seventyfifth = quantile(alpha.div.sh, .75),
                max = max(alpha.div.sh))
    
    newRow = c(a,b, alpha_div$avgDiv, alpha_div$sd, alpha_div$min, alpha_div$twentyfifth, alpha_div$fifty, alpha_div$seventyfifth, alpha_div$max)
    
    print(newRow)
    
    normalDF <- rbind(normalDF, setNames(newRow,names(normalDF)))
  }
}


normalDF <- normalDF |> 
  slice(-1)



### Summarize sample sizes div by depth
## zoom in to upper corner alpha = 1 ->5, beta = 1 ->5
## is optimal better than raw average?
## different thing to optimize?

# calc raw mean
edna_agg <- edna_imputed |>
  mutate(depthm = as.numeric(depthm),
         weight = 1) |>
  group_by(cruise, line, sta) |>
  summarize(across(starts_with('asv'), 
                   ~weighted.mean(log(.x), weight)), ## automatically normalizes weights
            .groups = 'drop') |>
  group_by(cruise, line) |>
  summarize(across(starts_with('asv'), ~mean(.x))) |>
  group_by(cruise) |>
  summarize(across(starts_with('asv'), ~mean(.x)))

edna_data <- edna_agg |> mutate(exp(pick(-cruise)))

# split into sample info and asvs
asv <- edna_data |> 
  dplyr::select(starts_with('asv'))

sampinfo <- edna_data |> 
  dplyr::select(-starts_with('asv'))

# alpha diversity measures
alpha_divs_raw <- sampinfo %>%
  bind_cols(alpha.div.sh = diversity(asv, index = 'shannon'))

## Raw mean
alpha_divs_raw |> 
  summarise(meanDiv = mean(alpha.div.sh))


######################################################
# More Precise Grid Search alpha(2:5) beta(1:5) by 0.1
######################################################

depth_weight_fn <- function(depth, alpha, beta){
  dgamma(depth, shape = alpha, scale = beta) ##CHANGE????
}

alpha.vec <- seq(2,5,0.1)
beta.vec <- seq(1,5,0.1)

gs2 <- data.frame(0,0,0,0,0,0,0,0,0)
names(gs2) = c("alpha", "beta", "avgDiv", "SD", "Min", "25", "50","75", "Max")


for (a in alpha.vec){
  for (b in beta.vec){
    edna_agg <- edna_imputed |>
      mutate(depthm = as.numeric(depthm),
             weight = depth_weight_fn(depthm, a, b)) |>
      group_by(cruise, line, sta) |>
      summarize(across(starts_with('asv'), 
                       ~weighted.mean(log(.x), weight)), ## automatically normalizes weights
                .groups = 'drop') |>
      group_by(cruise, line) |>
      summarize(across(starts_with('asv'), ~mean(.x))) |>
      group_by(cruise) |>
      summarize(across(starts_with('asv'), ~mean(.x)))
    
    edna_data <- edna_agg |> mutate(exp(pick(-cruise)))
    
    # split into sample info and asvs
    asv <- edna_data |> 
      dplyr::select(starts_with('asv'))
    
    sampinfo <- edna_data |> 
      dplyr::select(-starts_with('asv'))
    
    # alpha diversity measures
    alpha_divs <- sampinfo %>%
      bind_cols(alpha.div.sh = diversity(asv, index = 'shannon')) 
    
    alpha_div <- alpha_divs |> 
      slice(-14) |> 
      summarise(avgDiv = mean(alpha.div.sh),
                sd = sd(alpha.div.sh),
                min = min(alpha.div.sh),,
                twentyfifth = quantile(alpha.div.sh, .25),
                fifty = quantile(alpha.div.sh, .50),
                seventyfifth = quantile(alpha.div.sh, .75),
                max = max(alpha.div.sh))
    
    newRow = c(a,b, alpha_div$avgDiv, alpha_div$sd, alpha_div$min, alpha_div$twentyfifth, alpha_div$fifty, alpha_div$seventyfifth, alpha_div$max)
    
    print(newRow)
    
    gs2 <- rbind(gs2, setNames(newRow,names(gs2)))
  }
}


gs2 <- gs2 |> 
  slice(-1)

save(gs2, file = "rslt/grid-search-narrow.RData")

gs2 |> 
  filter(avgDiv == max(gammaDF$avgDiv))

library(rgl)

with(gs2, plot3d(x = alpha, y = beta, z = avgDiv))

with(gs2, plot3d(x = alpha, y = beta, z = SD))


############################################################
## Binned depths ###########################################
############################################################

binnedDepth <- edna_imputed |> 
  mutate(depthRNG = cut_number(as.numeric(depthm),3)) 

binnedDepth |> 
  dplyr::select(depthRNG)

gs3 <- data.frame(0,0,0,0,0,0,0,0,0,0)
names(gs3) = c("w1", "w2", "w3", "avgDiv", "SD", "Min", "25", "50","75", "Max")


binned_weight_fn <- function(depthRNG, w1, w2, w3){
  return(ifelse(depthRNG == "[0,10]", w1, ifelse(depthRNG == "(10,36.3]", w2, w3)))
}

w1.vec = c(0.1,0.5,1)
w2.vec = c(0.1,0.5,1)
w3.vec = c(0.1,0.5,1)

for (w1 in w1.vec) {
  for (w2 in w2.vec) {
    for (w3 in w3.vec) {
      
      edna_agg <- binnedDepth |>
        mutate(depthm = as.numeric(depthm),
               weight = binned_weight_fn(depthRNG, w1, w2, w3)) |>
        group_by(cruise, line, sta) |>
        summarize(across(starts_with('asv'), 
                         ~weighted.mean(log(.x), weight)), ## automatically normalizes weights
                  .groups = 'drop') |>
        group_by(cruise, line) |>
        summarize(across(starts_with('asv'), ~mean(.x))) |>
        group_by(cruise) |>
        summarize(across(starts_with('asv'), ~mean(.x)))
      
      edna_data <- edna_agg |> mutate(exp(pick(-cruise))) 
      
      # split into sample info and asvs
      asv <- edna_data |> 
        dplyr::select(starts_with('asv'))
      
      sampinfo <- edna_data |> 
        dplyr::select(-starts_with('asv'))
      
      # alpha diversity measures
      alpha_divs <- sampinfo %>%
        bind_cols(alpha.div.sh = diversity(asv, index = 'shannon')) 
      
      alpha_div <- alpha_divs |> 
        #slice(-14) |> 
        summarise(avgDiv = mean(alpha.div.sh),
                  sd = sd(alpha.div.sh),
                  min = min(alpha.div.sh),,
                  twentyfifth = quantile(alpha.div.sh, .25),
                  fifty = quantile(alpha.div.sh, .50),
                  seventyfifth = quantile(alpha.div.sh, .75),
                  max = max(alpha.div.sh))
      
      newRow = c(w1,w2,w3, alpha_div$avgDiv, alpha_div$sd, alpha_div$min, alpha_div$twentyfifth, alpha_div$fifty, alpha_div$seventyfifth, alpha_div$max)
      
      print(newRow)
      
      gs3 <- rbind(gs3, setNames(newRow,names(gs3)))
      
    }
  }
}

gs3 <- gs3 |> 
  slice(-1)


  

          