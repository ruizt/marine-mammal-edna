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
  filter(prop.nz > 0.05, prop.nz < 0.9) |> 
  pull(col)

# filter to samples with at least 25% nonzero reads
rows_of_interest <- edna_samples |> 
  dplyr::select(all_of(cols_of_interest)) |>
  mutate(across(where(is.numeric), ~ .x > 0)) |>
  mutate(total = rowMeans(pick(starts_with('asv'))),
         rownum = row_number()) |>
  dplyr::select(total, rownum) |>
  filter(total > 0.1) |>
  pull(rownum)

# zero imputation (bayesian multiplicative, martin-fernandez 2015)
imputation_out <- edna_samples |> 
  dplyr::select(all_of(cols_of_interest)) |>
  slice(rows_of_interest) |>
  zCompositions::cmultRepl(label = 0, 
                           method = 'GBM', 
                           output = 'prop',
                           z.warning = 0.99)

# bind imputed values to sample info
edna_imputed <- edna_samples |> 
  slice(rows_of_interest) |>
  dplyr::select(1:5) |>
  bind_cols(imputation_out)




##############################################################
### Gamma Weight Function ####################################
##############################################################

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


##save(gs1, file = "rslt/grid-search-result-wide.RData")


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

##save(gs2, file = "rslt/grid-search-narrow.RData")

gs2 |> 
  filter(avgDiv == max(gs2$avgDiv))

gs2 |> 
  filter(SD == max(gs2$SD))
library(rgl)

with(gs2, plot3d(x = alpha, y = beta, z = avgDiv))

with(gs2, plot3d(x = alpha, y = beta, z = SD))


############################################################
## Binned depths ###########################################
############################################################

binnedDepth <- edna_imputed |> 
  mutate(depthRNG = case_when(as.numeric(depthm) <= 15 ~ "Surface",
                              as.numeric(depthm) > 15 ~ "Deep")) 

## Try 2 bins, maybe surface and below surface


gs3 <- data.frame(0,0,0,0,0,0,0,0,0)
names(gs3) = c("w1", "w2", "avgDiv", "SD", "Min", "25", "50","75", "Max")


binned_weight_fn <- function(depthRNG, weight1){
  return(ifelse(depthRNG == "Surface", w1, 1 - w1))
}

w1.vec = seq(0.05,0.99, by = 0.005)


for (w1 in w1.vec) {
      
      w2 = (1 - w1)
  
      edna_agg <- binnedDepth |>
        drop_na() |> 
        mutate(depthm = as.numeric(depthm),
               weight = binned_weight_fn(depthRNG, w1)) |>
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
      
      newRow = c(w1,w2, alpha_div$avgDiv, alpha_div$sd, alpha_div$min, alpha_div$twentyfifth, alpha_div$fifty, alpha_div$seventyfifth, alpha_div$max)
      
      print(newRow)
      
      gs3 <- rbind(gs3, setNames(newRow,names(gs3)))
      
    }
  

gs3 <- gs3 |> 
  slice(-1)


  
gs3 |> 
  filter(avgDiv == max(gs3$avgDiv))

gs3 |> 
  filter(SD == max(gs3$SD))

## Check: across all stations, shallowest depth for each station, greatest depth fr "surface level" observation
## Same thing for second shallowest depth
## find min of these
## group by cruise line, station or create a factor: is depth the minimum depth?

## Separate shallowest obs vs deeper obs

#### calculate avg div by depth range

######################################
## Create depth range dataframe 
######################################

min_depths <- edna_imputed |> 
  group_by(cruise, line, sta) |> 
  summarise(min.depth = min(as.numeric(depthm)), 
            md.binary = 1)

binnedDepth <- edna_imputed |> 
  mutate(depthm = as.numeric(depthm)) |> 
  left_join(min_depths, by = join_by(cruise,line,sta,depthm == min.depth))


binnedDepth <- binnedDepth |> 
  mutate(md.binary = replace_na(md.binary, 0))

binnedDepth <- binnedDepth |> 
  mutate(depth.range = case_when(md.binary == 1 ~ "Surface", 
                                 md.binary == 0 ~ "Deep"))

binnedDepth$depth.range <- factor(binnedDepth$depth.range, levels = c("Surface", "Deep"))

binnedDepth <- binnedDepth |> 
  mutate(depth.range = if_else(depth.range == "Surface" & depthm > 30, "Deep", depth.range))


####################################################
##### ALPHA DIV BY BIN EXPLORATION
####################################################

# split into sample info and asvs
asv <- binnedDepth |> 
  dplyr::select(starts_with('asv'))

sampinfo <- binnedDepth |> 
  dplyr::select(-starts_with('asv'))

# alpha diversity measures
alpha_divs <- sampinfo %>%
  bind_cols(alpha.div.sh = diversity(asv, index = 'shannon')) 

alpha_divs |> 
  filter(depth.range == "Surface",
         depthm > 50)

alpha_divs |> 
  group_by(depth.range) |> 
  ggplot((aes(x = depth.range, y = alpha.div.sh))) + 
  geom_boxplot() + 
  labs(y = "Alpha Diversity", x = "Depth Range", title = "Boxplot: Alpha Diversity By Depth")
  
#ggsave("rslt/plots/np_surface_vs_deep_BP.png")

alpha_divs |> 
  ggplot(aes(x = depthm, y = alpha.div.sh, col = depth.range)) +
  geom_point() + scale_color_manual(values = c("red2", "blue4")) + 
  labs(y = "Alpha Diversity", x = "Depth (Meters)", title = "Scatterplot: Alpha Diversity By Depth", color = "Depth")
ggsave("rslt/plots/np_surface_vs_deep_SP.png")

alpha_divs |> 
  filter(depth.range == "Surface",
         depthm > 30)


##############################################################################

# Wide Grid Search

gs4 <- data.frame(0,0,0,0,0,0,0,0,0)
names(gs4) = c("w1", "w2", "avgDiv", "SD", "Min", "25", "50","75", "Max")


binned_weight_fn <- function(depth.range, weight1){
  return(ifelse(depth.range == "Surface", w1, 1 - w1))
}

w1.vec = seq(0.05,0.99, by = 0.05)


for (w1 in w1.vec) {
  
  w2 = (1 - w1)
  
  edna_agg <- binnedDepth |>
    drop_na() |> 
    mutate(depthm = as.numeric(depthm),
           weight = binned_weight_fn(depth.range, w1)) |>
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
  
  newRow = c(w1,w2, alpha_div$avgDiv, alpha_div$sd, alpha_div$min, alpha_div$twentyfifth, alpha_div$fifty, alpha_div$seventyfifth, alpha_div$max)
  
  print(newRow)
  
  gs4 <- rbind(gs4, setNames(newRow,names(gs4)))
  
}


gs4 <- gs4 |> 
  slice(-1)



gs4 |> 
  filter(avgDiv == max(gs4$avgDiv))

gs4 |> 
  filter(SD == max(gs4$SD))


gs4 |> 
  ggplot(aes(x = w1, y = avgDiv)) +
  geom_point() + labs(x= "Weight of Surface Measurements")

# Results: Max alpha div with surface weight = 0.1, deep weight = 0.9
#          Max SD with surface weight = 0.95, deep weight = 0.05


#################################################################

# More precise grid search: High deep weight and low surface weight 


gs5 <- data.frame(0,0,0,0,0,0,0,0,0)
names(gs5) = c("w1", "w2", "avgDiv", "SD", "Min", "25", "50","75", "Max")


binned_weight_fn <- function(depth.range, weight1){
  return(ifelse(depth.range == "Surface", w1, 1 - w1))
}

w1.vec = seq(0.05,0.15, by = 0.001)


for (w1 in w1.vec) {
  
  w2 = (1 - w1)
  
  edna_agg <- binnedDepth |>
    drop_na() |> 
    mutate(depthm = as.numeric(depthm),
           weight = binned_weight_fn(depth.range, w1)) |>
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
  
  newRow = c(w1,w2, alpha_div$avgDiv, alpha_div$sd, alpha_div$min, alpha_div$twentyfifth, alpha_div$fifty, alpha_div$seventyfifth, alpha_div$max)
  
  print(newRow)
  
  gs5 <- rbind(gs5, setNames(newRow,names(gs5)))
  
}


gs5 <- gs5 |> 
  slice(-1)



gs5 |> 
  filter(avgDiv == max(gs5$avgDiv))

gs5 |> 
  filter(SD == max(gs5$SD))


save(gs4, file = "rslt/grid-search-wide.RData")
save(gs5, file = "rslt/grid-search-narrow.RData")


gs5 |> 
  ggplot(aes(x = w1, y = avgDiv)) +
  geom_point() + labs(x= "Weight of Surface Measurements")



### Implement this all into preprocessing.R
### Save plots somewhere nice
### Filter out "surface" measurements deeper than 30 meters