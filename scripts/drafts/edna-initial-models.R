library(compositions)
library(pls)
library(spls)
library(tidyverse)
library(caret)

# load edna_data using EITHER 16s or 18s (note: edna_data is the same name in both...)
load("data/edna_18s_processed.RData")
load("data/edna_16s_processed.RData")
#load("data/ncog-18s.RData")
#ls() to see file names

head(edna_data)

# compute closure using centered logratio transform
cruises <- edna_data %>% select(cruise)
# transform log contrast
edna_data_lc <- edna_data %>% select(-cruise) %>% 
  clr() %>% 
  clrInv() %>% as.data.frame() 
# combine cruise and log ratio-transformed compositions
edna_data_lc <- cbind(cruises, edna_data_lc)


# add in response of sightings `Freq`
#sightings <- readxl::read_xlsx("data/CalCOFI_2004-2021_CombinedSightings.xlsx")
#sightings <- as.data.frame(table(sightings$Cruise))
#sightings$Var1 <- str_replace(sightings$Var1, "CC", "X20")
#head(sightings)
#asv_sightings <- left_join(edna_data_lc, sightings, by = c("cruise" = "Var1"))


# response variable: density of blue whales `Bm_scaled`
scaled_sightings <- read_csv("data/CC_on_effort_scaled_sightings_Ceta_SCB.csv")
scaled_sightings$cruiseID <- str_replace(scaled_sightings$cruiseID, "CC", "X20")
asv_sightings <- left_join(edna_data_lc, scaled_sightings, by = c("cruise" = "cruiseID")) 
asv_sightings <- asv_sightings %>% drop_na(Bm_scaled) 
Bm_scaled <- asv_sightings$Bm_scaled

asv_predictors <- asv_sightings %>% select(starts_with("asv"))
asv_sightings <- cbind(Bm_scaled, asv_predictors)


set.seed(1)


# run PCR
model <- pcr(Bm_scaled ~ ., data=asv_sightings, validation="CV")
summary(model)
model$scores

# run PLSR
model <- plsr(Bm_scaled ~ ., ncomp = 10, data=asv_sightings, validation = "CV")
summary(model)

var(model$residuals)/var(asv_sightings$Bm_scaled)

model$scores

#visualize cross-validation plots
validationplot(model)
validationplot(model, val.type="MSEP")
validationplot(model, val.type="R2")


# rough variance reduction
var(model$residuals)/var(asv_sightings$Bm_scaled)

cbind(model$fitted.values, Bm_scaled)

####################### sPLS
set.seed(4)
# cross validate using K=2 to 8 latent components and eta = 0.1 to 0.9
cv <- cv.spls(asv_predictors, Bm_scaled, 
              eta = seq(0.1,0.9,0.1), K = c(2:8))

# fit model with best K and eta
f <- spls(x = asv_predictors, y = asv_sightings$Bm_scaled, K = cv$K.opt, eta = cv$eta.opt)

pred <- predict.spls(f, type = 'fit')

var(Bm_scaled - pred)/var(Bm_scaled)
var(Bm_scaled - model$fitted.values)/var(Bm_scaled)

# show results of selected variables
print(f)

# estimate and show the coefficients
coef.f <- coef(f) %>% as.data.frame() %>% rename(estimate = V1)
coef.f %>% filter(estimate != 0)

# plot coefficients
plot.spls(f)


# Calculate bootstrapped confidence intervals of SPLS coefficients
ci.f <- ci.spls(f, plot.it=TRUE, plot.fix='x', plot.var=1)
cis <- ci.f$cibeta
cis[[1]]

cf <- correct.spls( ci.f ) %>% as.data.frame() %>% rename(estimate = V1)
cf %>% filter(estimate != 0)

# rough variance explained
pred <- predict.spls(f, type = 'fit')
var(Bm_scaled - pred)/var(Bm_scaled)

cbind(pred, Bm_scaled)

################### PLS with LASSO
model <- train(
  Bm_scaled ~ starts_with("asv"),
  data = asv_sightings,
  preProcess = c("center", "scale"),  
  method = 'pls'
)
model
plot(model)
