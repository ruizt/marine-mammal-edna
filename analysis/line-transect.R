## LINE TRANSECT ANALYSIS
#
# Distance-sampling analysis of CalCOFI whale survey data.
# Produces density estimates for blue, fin, and humpback whales by cruise.
#
# Inputs (from data/, downloaded via setup-data.R):
#   sample_table.Rds, region_table.Rds, obs_table.Rds
#
# Executable R code extracted from vignettes/line-transect.qmd.

#Load required libraries
library(readxl)
library(Distance)
library(tidyverse)
library(kableExtra)
library(ggplot2)
set.seed(71246)

sample_table <- read_rds(here::here("data/sample_table.Rds"))
region_table <- read_rds(here::here("data/region_table.Rds"))
obs_table <- read_rds(here::here("data/obs_table.Rds"))

region_table %>% kbl(longtable =  TRUE) %>% 
  kable_classic(full_width = FALSE)

ggplot(region_table, aes(x = start_date)) +
  geom_line( aes(y = L), color = "red") + 
  geom_line( aes(y = k * 100), color = "blue", linetype = "dashed") +
  geom_point( aes(y = L), size = ifelse(region_table$season == "winter", 2, 0)) +
  scale_y_continuous(
    name = "Line length",
    sec.axis = sec_axis(~ . / 100, name="Days of surveying"),
    limits = c(0, 1600)
  ) + 
  xlab("Start date of cruise") +
  theme(
    axis.title.y = element_text(color = "red"),
    axis.title.y.right = element_text(color = "blue")
  )

#Left truncation
w <- 2400
obs_table <- obs_table %>% filter(distance <= w)

#Summary of observations
obs_summary <- obs_table %>%
  group_by(species) %>%
  summarise(n = n(), 
            sizeL = quantile(size, 0.05), sizeM = median(size), sizeU = quantile(size, 0.95), 
            ssL = quantile(ss, 0.05), ssM = median(ss), ssU = quantile(ss, 0.95), 
            swellL = quantile(swell, 0.05), swellM = median(swell), swellU = quantile(swell, 0.95))
obs_summary %>% kbl(col.names = c("species", "n", 
                                  "size 5", "size 50", "size 95", 
                                  "ss 5", "ss 50", "ss 95", 
                                  "swell 5", "swell 50", "swell 95")) %>% 
  kable_classic(full_width = FALSE)

#Define function to fit detection functions to data passed in as data frame/tibble dat. 
# w is the truncation distance.
# fits a set of key function and adjustments with forward selection of covariates, specified in the function body
# returns a data frame with one row for each fitted detection function, and an AIC column
#Note - season not passed in by default
fit_det_func <- function(dat, w,
  keys = c("unif", "hn", "hr"),
  adjustments = c("cos", "cos", "poly"),
  max_adjustments = c(3, 2, 1),
  covariates = c("ss", "swell", "size", "f_obs_height"),
  covariates_and_adjustments = FALSE){
  
  #Note only fits covariates for hn and hr key functions
  
  n_covariates <- length(covariates)
  model <- 0
  for(key_series in 1:length(keys)){
  
    key <- keys[key_series]
    adjustment <- adjustments[key_series]
    for(nadj in 0:max_adjustments[key_series]){
      
      #Fit model with no covariates - assumes no error will be generated here
      res <- ds(dat, w, key = key, adjustment = adjustment, nadj = nadj, quiet = TRUE)
      model <- model + 1
      min_AIC_within_keyfct <- AIC(res)$AIC
      if(model == 1){
        res_table <- data.frame(model = model, key = key, adjustment = adjustment,
                                nadj = nadj, covar = "", AIC = min_AIC_within_keyfct)
      } else {
        res_table <- rbind(res_table, 
                           data.frame(model = model, key = key, adjustment = adjustment,
                                nadj = nadj, covar = "", AIC = min_AIC_within_keyfct))
      }
      
      if((key %in% c("hn", "hr")) & ((nadj == 0) | covariates_and_adjustments)){
        #Forward selection for model covariates
        round <- 1
        covar_selected <- rep(FALSE, n_covariates)
        while(round <= n_covariates){
          res_AIC <- rep(Inf, n_covariates)
          for(cov in 1:n_covariates){
            #Skip this model if the covariate is already selected
            if(covar_selected[cov]) next
            #Build up the covariate string
            if(any(covar_selected)) {
              covar_text <- paste(paste(covariates[covar_selected], collapse = " + "), covariates[cov], sep = " + ")
            } else {
              covar_text <- covariates[cov]
            }
            #Fit the model
            res <- try(ds(dat, w, key = key, adjustment = adjustment, nadj = nadj,
                      formula = eval(parse(text = paste("~", covar_text))), quiet = TRUE))
            if(inherits(res, "try-error")){
              #There was an error during the fitting - do not save the result, and move on to the next model
              next()
            } else {
              #Save the result
              model <- model + 1
              res_AIC[cov] <- AIC(res)$AIC
              res_table <- rbind(res_table, 
                                 data.frame(model = model, key = key, adjustment = adjustment,
                                      nadj = nadj, covar = covar_text, AIC = res_AIC[cov]))
            }
          }
          #Select best model in this round
          best_model <- which.min(res_AIC)
          if(res_AIC[best_model] < min_AIC_within_keyfct){
            min_AIC_within_keyfct <- res_AIC[best_model]  
            covar_selected[best_model] <- TRUE
            round <- round + 1
          } else{
            break
          }
        }
      }
      
      #Determine whether to fit more adjustments
      if(nadj == 0){
        min_AIC <- min_AIC_within_keyfct
      } else {
        if(min_AIC_within_keyfct < min_AIC){
          #Go on to next adjustment term
          min_AIC <- min_AIC_within_keyfct
        } else {
          #No more adjustments
          break
        }
      }
      
    } #end adjustments
  } #end key_series

  res_order <- order(res_table$AIC)
  res_table <- res_table[res_order, ]
  res_table$Delta_AIC <- res_table$AIC - res_table$AIC[1]
  return(res_table)
}

#Fit each species in turn

species <- unique(obs_table$species)
n_species <- length(species)
convert_units <- 1/1000000
fits <- det_func_tables <- vector("list", n_species)
names(fits) <- names(det_func_tables) <- species

for(i in 1:n_species){
  
  dat <- obs_table[obs_table$species == species[i], ]
  
  #Fit the detection functions
  det_func_tables[[i]] <- fit_det_func(dat, w)
  
  #Refit with the AIC_best model
  if(det_func_tables[[i]]$covar[1] == "") {
    fits[[i]] <- ds(dat, w, key = det_func_tables[[i]]$key[1], 
                    adjustment = det_func_tables[[i]]$adjustment[1], nadj = det_func_tables[[i]]$nadj[1], 
                    region_table = region_table, sample_table = sample_table,
                    convert_units = convert_units, quiet = TRUE)
  } else {
    fits[[i]] <- ds(dat, w, key = det_func_tables[[i]]$key[1], 
                    adjustment = det_func_tables[[i]]$adjustment[1], 
                    nadj = det_func_tables[[i]]$nadj[1], region_table = region_table, 
                    sample_table = sample_table,
                    formula = eval(parse(text = paste("~", det_func_tables[[i]]$covar[1]))), 
                    convert_units = convert_units, quiet = TRUE)
  }

}

tmp <- det_func_tables[2][[1]][, -1] %>% 
  kbl(digits = 2, row.names = FALSE, longtable = TRUE) %>% 
  kable_classic(full_width = FALSE)
print(tmp)

summary_BM <- summary(fits[[2]]$ddf)
print(summary_BM)

plot(fits[[2]])

gof_BM <-gof_ds(fits[[2]])

tmp <- det_func_tables[1][[1]][, -1] %>% 
  kbl(digits = 2, row.names = FALSE, longtable = TRUE) %>% 
  kable_classic(full_width = FALSE)
print(tmp)

summary_BP <- summary(fits[[1]]$ddf)
print(summary_BP)

plot(fits[[1]])

gof_BP <-gof_ds(fits[[1]])

tmp <- det_func_tables[3][[1]][, -1] %>% 
  kbl(digits = 2, row.names = FALSE, longtable = TRUE) %>% 
  kable_classic(full_width = FALSE)
print(tmp)

summary_MN <- summary(fits[[3]]$ddf)
print(summary_MN)

plot(fits[[1]])

gof_MN <- gof_ds(fits[[1]])

tmp <- summary(fits[[2]])$dht$individuals$D
tmp2 <- tmp %>% 
  kbl(digits = 2, longtable = TRUE) %>%
  kable_classic(full_width = FALSE)
print(tmp2)

tmp <- summary(fits[[1]])$dht$individuals$D
tmp2 <- tmp %>% 
  kbl(digits = 2, longtable = TRUE) %>%
  kable_classic(full_width = FALSE)
print(tmp2)

tmp <- summary(fits[[3]])$dht$individuals$D
tmp2 <- tmp %>% 
  kbl(digits = 2, longtable = TRUE) %>%
  kable_classic(full_width = FALSE)
print(tmp2)

sessionInfo()

