library(tidyverse)
library(spls)
library(modelr)

## PREDICTION OPTIMAL MODELS ---------------------------------------------------

# model fitting function
fit_fn <- function(.data, .var, .parms){
  x <- dplyr::select(.data, starts_with('asv'))
  y <- pull(.data, {{.var}})
  fit <- spls(x, y, 
              K = .parms['ncomp'], eta = .parms['eta'], 
              scale.x = F, scale.y = F)
}

# fitting 3 optimal models
fit_bm <- fit_fn(whales, bm, c(ncomp = 5, eta = 0.86))
fit_bp <- fit_fn(whales, bp, c(ncomp = 6, eta = 0.967))
fit_mn <- fit_fn(whales, mn, c(ncomp = 7, eta = 0.611))

## Model Results: BM, ncomp = 5, eta = 0.86

# R-squared
n <- length(fit_bm$y)
bm_fitted <- predict(fit_bm, type = 'fit')
bm_resid <- fit_bm$y - fitted
bm_r2 <- 1 - (var(resid)/var(fit_bm$y))
bm_adj_r2 <- 1 - ((1-r2)*(n-1))/(n-fit_bm$K-1)

bm_r2
bm_adj_r2

# loocv for mspe
loo_preds <- rep(NA, nrow(fit_bm$x))

for (i in 1:nrow(x)){
  x_train <- fit_bm$x[-i, ] # removes ith row
  y_train <- fit_bm$y[-i]  # removes ith element
  
  
  # fit cross validation model
  fit_cv <- spls(x_train, y_train, K = 5, eta = 0.86, 
                 scale.x = F, scale.y = F)
  
  # predict response variable
  loo_preds[i] <- predict(fit_cv, newx = x[i, ], type = "fit")
}

# mspe (mean squared prediction error)
bm_mspe <- mean((fit_bm$y - loo_preds)^2)

model_results <- data_frame("bm", 5, 0.86, bm_r2[[1]], bm_adj_r2[[1]], bm_mspe)

names(model_results) <- c("species", "ncomp", "eta", "r2", "adj.r2", "mspe")

# Model Results: BP, ncomp = 6, eta = 0.967

# R-squared
n <- length(fit_bp$y)
bp_fitted <- predict(fit_bp, type = 'fit')
bp_resid <- fit_bp$y - fitted
bp_r2 <- 1 - (var(resid)/var(fit_bp$y))
bp_adj_r2 <- 1 - ((1-r2)*(n-1))/(n-fit_bp$K-1)

bp_r2
bp_adj_r2

# loocv for mspe
loo_preds <- rep(NA, nrow(fit_bp$x))

for (i in 1:nrow(x)){
  x_train <- fit_bp$x[-i, ] # removes ith row
  y_train <- fit_bp$y[-i]  # removes ith element
  
  
  # fit cross validation model
  fit_cv <- spls(x_train, y_train, K = 5, eta = 0.86, 
                 scale.x = F, scale.y = F)
  
  # predict response variable
  loo_preds[i] <- predict(fit_cv, newx = x[i, ], type = "fit")
}

# mspe (mean squared prediction error)
bp_mspe <- mean((fit_bp$y - loo_preds)^2)

model_results <- rbind(model_results, c("bp", 6, 0.967, bp_r2[[1]], bp_adj_r2[[1]], bp_mspe))

# Model Results: MN, ncomp = 7, eta = 0.611

# R-squared
n <- length(fit_mn$y)
mn_fitted <- predict(fit_mn, type = 'fit')
mn_resid <- fit_mn$y - fitted
mn_r2 <- 1 - (var(resid)/var(fit_mn$y))
mn_adj_r2 <- 1 - ((1-r2)*(n-1))/(n-fit_mn$K-1)

mn_r2
mn_adj_r2

# loocv for mspe
loo_preds <- rep(NA, nrow(fit_bp$x))

for (i in 1:nrow(x)){
  x_train <- fit_mn$x[-i, ] # removes ith row
  y_train <- fit_mn$y[-i]  # removes ith element
  
  
  # fit cross validation model
  fit_cv <- spls(x_train, y_train, K = 5, eta = 0.86, 
                 scale.x = F, scale.y = F)
  
  # predict response variable
  loo_preds[i] <- predict(fit_cv, newx = x[i, ], type = "fit")
}

# mspe (mean squared prediction error)
mn_mspe <- mean((fit_mn$y - loo_preds)^2)

model_results <- rbind(model_results, c("mn", 7, 0.611, mn_r2[[1]], mn_adj_r2[[1]], mn_mspe))

## OPTIMAL MODEL RESULTS

model_results

#-------------------------------------------------------------------------------



## STABLE SET APPROACH ---------------------------------------------------------

# read in selection frequencies from LOOCV
sel_freq <- read_rds('rslt/loocv/2024-07-27/selection-frequencies.rds')

# read in metrics from LOOCV
metrics <- read_rds('rslt/loocv/2024-07-27/metrics.rds')

## SPECIFYING NCOMP/ETA RANGE --------------------------------------------------
# chosen stability threshold (minimum max selection prob.)
pimax <- 0.8

# desired bound for expected number of false positives
EV <- 5

# total number of candidate ASVs
p <- 6832

# threshold for max average number of selected ASVs
qmax <- sqrt((2*pimax - 1)*p*EV) 

## GRID SEARCH FOR ETA RANGES  -------------------------------------------------

# Modified grid search
# grid search: FIXED ncomp, eta range -> average # of asvs selected
ncomp_grid <- 5:8
eta_min_grid <- seq(0.4,0.88, by = 0.02)
eta_max_grid <- seq(0.88,0.4, by = -0.02)
species_grid <- c("bm", "bp", "mn")

num_asvs_eta_ncomp_gs <- function(ncomp_val, eta_min, eta_max){
  avg_asvs <- metrics |> 
    filter(ncomp == ncomp_val,
           eta >= eta_min,
           eta <= eta_max) |> 
    group_by(species) |> 
    summarise(avg.n.asv = mean(n.asv),
              sd.n.asv = sd(n.asv),
              .groups = 'drop')
  
  avg_asvs <- avg_asvs |> 
    mutate(ncomp = ncomp_val,
           eta.min = eta_min,
           eta.max = eta_max,
           exp.fp = avg.n.asv^2/(p*(2*pimax - 1)))
  
  col.order <- c("species", "ncomp", "eta.min", "eta.max", "avg.n.asv", "sd.n.asv", "exp.fp")
  
  avg_asvs <- avg_asvs[,col.order]
  
  return(avg_asvs)
}

num_asvs_eta_ncomp_gs(6,0.6,0.85)


gs_df <- data.frame(species = character(), 
                    ncomp = numeric(), 
                    eta.min = numeric(), 
                    eta.max = numeric(), 
                    avg.n.asv = numeric(), 
                    sd.n.asv = numeric(), 
                    exp.fp = numeric(),
                    stringsAsFactors = FALSE)

for (n in ncomp_grid) {
  for (emin in eta_min_grid){
    for (emax in eta_max_grid){
      if (emin > emax){
          next
      }
      
        newrows = num_asvs_eta_ncomp_gs(n,emin,emax)
        
        print(newrows)
        
        gs_df <- rbind(gs_df, newrows)
        
      }
        
      }
  }
 

#-------------------------------------------------------------------------------


# Candidate ranges: wide enough range with expected false positives < 1
candidate_ranges <- gs_df |> 
  filter(exp.fp < 1,
         eta.max - eta.min > 0.3)

# compute stable sets for all candidate ranges (582)
gs_stable_sets <- data.frame(species = character(), 
                          ncomp = numeric(), 
                          eta.min = numeric(), 
                          eta.max = numeric(), 
                          avg.n.asv = numeric(), 
                          sd.n.asv = numeric(), 
                          exp.fp = numeric(),
                          stable_sets = list(),
                          stringsAsFactors = FALSE)

for(i in 1:nrow(candidate_ranges)) {
  # get row
  row <- candidate_ranges[i,]
  
  # compute stable set
  ss <- sel_freq |> 
    filter(eta >= row$eta.min,
           eta <= row$eta.max,
           ncomp == row$ncomp,
           species == row$species, 
           n >= 20) |> 
    distinct(asv) |>
    nest(data = asv) |>
    transmute(stable.set = map(data, unlist))
  
  # add to df
  newrow <- cbind(row,ss)

  gs_stable_sets <- rbind(gs_stable_sets, newrow)
  
}

# STABLE SET EXPLORATION
# GOAL: FIND OPTIMAL NCOMP FOR EACH SPECIES (stable set that optimizes MSPE)
# TESTING NCOMP = 5, 6, 7, and 8

# size of stable sets
gs_stable_sets$ss.size <- sapply(gs_stable_sets$stable.set, length)


## find largest stable set for each species-ncomp combination (with expected FP < 1)

max_asvs_ss <- gs_stable_sets |> 
  group_by(species, ncomp) |> 
  slice_max(ss.size)

max_asvs_ss

# are the stable sets of same size identical?

identical(max_asvs_ss$stable.set[[2]], max_asvs_ss$stable.set[[3]])

identical(max_asvs_ss$stable.set[[4]], max_asvs_ss$stable.set[[5]])

identical(max_asvs_ss$stable.set[[6]], max_asvs_ss$stable.set[[7]])


# filter out one stable set per ncomp-species combo

stable_sets <- max_asvs_ss |> 
  group_by(species, ncomp) |> 
  slice_min(exp.fp) |> 
  slice_max(eta.max - eta.min)

# fit all 12 models to compare MSPE
load('data/processed/ncog18sv9-2024-07-20.RData')
load('data/processed/mm-sightings-2024-07-20.RData')
whales <- inner_join(sightings, edna, by = 'cruise') 


ss_fit_fn <- function(.data, .var, .ss, .ncomp){
  stable_asvs <- unlist(.ss, use.names = F)
  
  all_asvs <- dplyr::select(.data, starts_with('asv'))
  
  x <- all_asvs |> dplyr::select(all_of(stable_asvs))
  y <- pull(.data, {{.var}})
  fit <- spls(x, y, 
              K = .ncomp, eta = 0, 
              scale.x = F, scale.y = F)
}


model_outputs = data.frame(species = character(), 
                          ncomp = numeric(),
                          stable.set = list(),
                          r2 = numeric(),
                          adj.r2 = numeric(),
                          mspe = numeric()  
                          )


for (j in 1:nrow(stable_sets)) {
  row <- stable_sets[j,]
  print(row)
  stable_set = row$stable.set
  
  fit <- ss_fit_fn(whales, row$species, stable_set, row$ncomp)
  
  
  # R-squared
  n <- length(fit$y)
  fitted <- predict(fit, type = 'fit')
  resid <- fit$y - fitted
  r2 <- 1 - (var(resid)/var(fit$y))
  adj_r2 <- 1 - ((1-r2)*(n-1))/(n-fit$K-1)
  
  print(r2)
  print(adj_r2)
  
  
  # loocv for mspe
  loo_preds <- rep(NA, nrow(fit$x))
  
  for (i in 1:nrow(fit$x)){
    x_train <- fit$x[-i, ] # removes ith row
    y_train <- fit$y[-i]  # removes ith element
    
    
    # fit cross validation model
    fit_cv <- spls(x_train, y_train, K = row$ncomp, eta = 0, 
                   scale.x = F, scale.y = F)
    print(fit$x[i,])
    # predict response variable
    loo_preds[i] <- predict(fit_cv, newx = fit$x[i, , drop = F], type = "fit")
  }
  
  # mspe (mean squared prediction error)
  mspe <- mean((fit$y - loo_preds)^2)
  
  newrow <- data.frame(species = row$species,
                       ncomp = row$ncomp,
                       stable.set = I(list(row$stable.set)),
                       r2 = r2,
                       adj.r2 = adj_r2,
                       mspe = mspe)
  print(names(newrow))
  print(names(model_outputs))
  
  model_outputs <- rbind(model_outputs, newrow)
  
}


model_outputs |>
  group_by(species) |>
  arrange(mspe)

# find best ncomp / stable set for each species
model_outputs |> 
  group_by(species) |> 
  slice_min(mspe)


# in terms of adjusted r^2
model_outputs |> 
  group_by(species) |> 
  slice_max(adj.r2)


# FROM OLD GRID SEARCH
# candidate_ranges <- gs_df |> 
#   filter(avg.n.csv <= qmax,
#          avg.n.csv > 0)
# 
# candidate_ranges <- candidate_ranges |> 
#   mutate(ncomp.range = ncomp.max - ncomp.min,
#          eta.range = eta.max - eta.min)
# 
# candidate_ranges |> 
#   filter(ncomp.range == max(ncomp.range) | eta.range == max(eta.range))
# 
# # notice effect of range on sd, for example...
# candidate_ranges |> 
#   filter(ncomp.range == 9) |>
#   filter(avg.n.csv > 40, avg.n.csv < 70) |>
#   arrange(eta.range)
# 
# # just tinkering here ...
# candidate_ranges |> 
#   filter(avg.n.csv < 90,
#          avg.n.csv > 20,
#          eta.max <= 0.9, 
#          eta.min >= 0.4,
#          ncomp.min > 2,
#          ncomp.max == 8,
#          avg.n.csv > 15,
#          sd.n.csv < 50)
# 
# # recalculate EV for new range
# qmax_new <- candidate_ranges |>
#   filter(ncomp.min == 3, ncomp.max == 8, eta.min == 0.5, eta.max == 0.9) |>
#   pull(avg.n.csv)
# 
# qmax_new^2/(p*(2*pimax - 1))

## STABLE SET -------------------------------------------------------------------
## INITIAL FILTERING: 1 <= ncomp <= 10 ; 0.4 <= eta <= 1
## => eta >= 0.4

# this is the widest range of ncomp and eta values that still obtains an 
# average number of asvs <= qmax (98.711)

# asvs that have a maximum selection prob of at least 0.8 (max n >= 20)

sel_freq |> 
  filter(species == "bp",
         eta >= 0.6,
         eta <= 0.9,
         ncomp >= 8,
         ncomp <= 8,
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


## another approach -- bin and group
pal.grad <- colorRampPalette(c('red', 'blue'))
pal <- pal.grad(13)

# look at avg nsel by eta bin and ncomp
# find ranges where qmin < q < qmax for all ncomp
metrics |>
  select(species, ncomp, eta, n.asv) |>
  mutate(eta.bin = cut_interval(eta, n = 8)) |>
  group_by(species, ncomp, eta.bin) |>
  summarize(avg.n.asv = mean(n.asv),
            sd.n.asv = sd(n.asv)) |>
  ggplot(aes(x = eta.bin, y = avg.n.asv, color = factor(ncomp))) +
  facet_wrap(~species) +
  geom_point() +
  geom_path(aes(group = ncomp)) +
  scale_y_log10() +
  scale_color_manual(values = pal) +
  theme(axis.text.x = element_text(angle = 90)) 

# say, qmin = 10, qmax = 80 groupwise; find q across full range
metrics |>
  select(species, ncomp, eta, n.asv) |>
  filter(ncomp == 8, 
         eta >= 0.65, eta <= 0.85) |>
  group_by(species, ncomp) |>
  summarize(avg.n = mean(n.asv), 
            sd.n = sd(n.asv)) |>
  slice_max(avg.n)

# upper bound on expected no. false positives
qmax <- 57
qmax^2/(p*(2*pimax - 1))

# stable sets
stable_sets <- sel_freq |> 
  filter(eta >= 0.65,
         eta <= 0.88,
         ncomp == 8,
         n >= 20) |> 
  group_by(species) |> 
  distinct(asv) |>
  nest(data = asv) |>
  transmute(stable.set = map(data, unlist))



# tinkering with fitting models -- it works!
load('data/processed/ncog18sv9-2024-07-20.RData')
load('data/processed/mm-sightings-2024-07-20.RData')
whales <- inner_join(sightings, edna, by = 'cruise') 

loocv_fn <- function(.train, .var){
  .asv <- filter(stable_sets, species == .var) %>% pull(stable.set) %$% .[[1]]
  names(.asv) <- NULL
  x.train <- .train |> as.data.frame() |> dplyr::select(all_of(.asv))
  y.train <- .train |> as.data.frame() |> pull({{.var}})
  fit <- spls(x.train, y.train,
              K = 8, eta = 0,
              scale.x = F, scale.y = F)

  return(fit)
}

pred_fn <- function(.test, .fit){
  asv <- rownames(.fit$projection)
  newdata <- .test |> as.data.frame() |> dplyr::select(all_of(asv))
  out <- predict(.fit, newx = newdata, type = 'fit')[1, 1]
  return(out)
}

loocv <- crossv_loo(whales) |>
  expand_grid(species = c('bm', 'bp', 'mn')) |>
  mutate(fit = map2(train, species,
                    ~loocv_fn(.x, .y)),
         pred = map2(test, fit, ~pred_fn(.x, .y)),
         obs = map2(test, species, ~pull(as.data.frame(.x), .y))) |>
  unnest(c(pred, obs))

loocv |>
  group_by(species) |>
  summarize(cor = cor(pred, obs),
            mspe = mean((pred - obs)^2))

loocv |>
  mutate(fitted = map(fit, ~predict(.x, type = 'fit')[,1]),
         train.obs = map2(train, species, ~pull(as.data.frame(.x), .y))) |>
  select(species, fitted, train.obs) |>
  mutate(cor = map2(fitted, train.obs, ~cor(.x, .y))) |>
  unnest(cor) |>
  group_by(species) |>
  summarize(mean(cor))


metrics |>
  group_by(ncomp, eta, species) |>
  summarize(mspe = mean(sq.pred.err),
            cor = cor(pred, pred + pred.err),
            df = mean(n.asv)) |>
  group_by(species) |>
  slice_max(cor)
