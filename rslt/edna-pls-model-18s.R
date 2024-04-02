###############################################################################
#
# Partial Least Squares (PLS) on 18s Data to Model Blue Whale Densities 
#
###############################################################################

# Set Up
###############################################################################
library(compositions)
library(pls)
library(spls)
library(tidyverse)
library(caret)

# load edna_data using 18s 
load("data/edna_18s_processed.RData")
#ls() to see file names
head(edna_data)

# read taxa data
edna_in <- read_tsv('data/NCOG_18sV9_asv_count_tax_S.tsv') %>%
  mutate(short.id = paste('asv', row_number(), sep = '.'))
taxa <- edna_in %>% 
  dplyr::select(where(is.character), silva_Confidence)


# compute closure using centered logratio transform
cruises <- edna_data %>% select(cruise)
# transform log contrast
edna_data_lc <- edna_data %>% select(-cruise) %>% 
  clr() %>% 
  clrInv() %>% as.data.frame() 
# combine cruise and log ratio-transformed compositions
edna_data_lc <- cbind(cruises, edna_data_lc)

# response variable: density of blue whales `Bm_scaled`
scaled_sightings <- read_csv("data/CC_on_effort_scaled_sightings_Ceta_SCB.csv")
scaled_sightings$cruiseID <- str_replace(scaled_sightings$cruiseID, "CC", "X20")

# adjust response for seasonality
scaled_sightings3 <- scaled_sightings %>% 
  mutate_all(~replace(., . == 0, 0.00000001)) %>% 
  group_by(season) %>% 
  mutate(Bm_scaled = log(Bm_scaled),
         Bm_season = mean(Bm_scaled, na.rm = TRUE),
         Bm_season = Bm_scaled - Bm_season) %>% 
  select(cruiseID, season, Bm_season, Bm_scaled)

# impute zeros with uniform random numbers between 0 and the minimum density

# scaled_sightings <- scaled_sightings %>% 
#   group_by(season) %>% 
#   mutate(Bm_season = mean(Bm_scaled, na.rm = TRUE),
#          Bm_season = Bm_scaled - Bm_season) %>% 
#   select(cruiseID, season, Bm_season, Bm_scaled) %>% 
#   # replace negatives w 0
#   rowwise() %>% 
#   mutate(Bm_season = max(c(0, Bm_season)))

asv_sightings <- left_join(edna_data_lc, scaled_sightings, by = c("cruise" = "cruiseID")) 
asv_sightings <- asv_sightings %>% drop_na(Bm_scaled) 
Bm_scaled <- asv_sightings$Bm_season   # orgininally Bm_scaled

# adjust ASV predictors for seasonality
asv_sightings <- asv_sightings %>% 
  group_by(season) %>% 
  # for each asv, subtract season avg
  mutate(across(starts_with("asv"),
                ~(.x-mean(.x)))) %>% 
  ungroup()

# dataset with only centered log transformed ASVs and blue whale density response
asv_predictors <- asv_sightings %>% select(starts_with("asv"))
asv_sightings <- cbind(Bm_scaled, asv_predictors)

# PLS Model
###############################################################################
# Run PLS model
model <- plsr(Bm_scaled ~ ., ncomp = 7, data=asv_sightings, validation = "CV")
summary(model) # note: 7 components captures 94.30% of variance in Bm_scaled
               # 6 explains 84.57% variance


# extract weights and estimates
l <- model$loadings 
loading_df <- data.frame(matrix(as.numeric(l), attributes(l)$dim, dimnames=attributes(l)$dimnames))
loading_df$asv = rownames(loading_df)

# estimates for 7 components
loading_df$coef_estimates <- model$coefficients %>% as.data.frame() %>% select(7) %>% pull()

weight_df <- loading_df  %>%
  rowwise() %>% 
  # for each component, max absolute loading as a weight for variable importance score
  mutate(weight = max(abs(c_across(Comp.1:Comp.7)))) %>%
  select(weight, coef_estimates) %>%
  bind_cols(short.id = loading_df$asv) %>%
  # incorporate taxa data
  left_join(taxa, by = "short.id") %>%
  select(short.id, silva_Taxon, weight, coef_estimates)

# clean dataframe
weight_df %>%
  separate(silva_Taxon, into = c('d', 'p', 'c', 'o', 'f', 'g'), sep = ';') %>%
  write_csv(file = 'rslt/taxa-weights-18s.csv')


#visualize cross-validation plots
validationplot(model)
validationplot(model, val.type="MSEP")
validationplot(model, val.type="R2")


# rough variance reduction
var(model$residuals)/var(asv_sightings$Bm_scaled)

# check by eye predicted vs actual
cbind(model$fitted.values, Bm_scaled)

# check residuals
plot(model$fitted.values, model$residuals)
plot(model)

# format taxa info in dataframe
clean_weight_df <- weight_df %>% 
  separate(silva_Taxon, into = c('d', 'p', 'c', 'o', 'f', 'g'), sep = ';') 

clean_weight_df %>% 
  arrange(desc(weight)) %>% 
  head(50) %>% 
  ggplot(aes(x = p, y = weight)) +
  geom_bar(aes(fill = short.id),
           position = "dodge",
           stat = "identity") + guides(fill="none") +
  theme_bw() +
  #scale_x_discrete(guide = guide_axis(n.dodge = 2)) + 
  theme(axis.text.x = element_text(angle = 45,hjust=1))

# Fit initial sPLS model
###############################################################################
set.seed(3)
# cross validate using K=2 to 8 latent components and eta = 0.1 to 0.9
cv <- cv.spls(asv_predictors, Bm_scaled, 
              eta = seq(0.1,0.9,0.1), K = c(2:8))

# fit model with best K and eta
#f <- spls(x = asv_predictors, y = asv_sightings$Bm_scaled, K = cv$K.opt, eta = cv$eta.opt)
f <- spls(x = asv_predictors, y = asv_sightings$Bm_scaled, K = 3, eta = cv$eta.opt)

# show results of selected variables
print(f)

# estimate and show the coefficients
coef.f <- coef(f) %>% as.data.frame() %>% rename(estimate = V1)
coef.f %>% filter(estimate != 0)

# plot coefficients
plot.spls(f)

# rough variance explained in spls model
pred <- predict.spls(f, type = 'fit')
var(Bm_scaled - pred)/var(Bm_scaled)
# compare to pls model
var(Bm_scaled - model$fitted.values)/var(Bm_scaled)


f$meanx[1:10]

f$projection

# predictions should add as 
# rep(f$mu, 25) + scale(as.matrix(asv_predictors)) %*% coef.f$estimates

## EXPERIMENTATION

tbl <- f$projection |>
  bind_cols(filter(coef.f, estimate != 0)) |>
  rename(reg.coef = estimate) |>
  rownames_to_column('short.id') |>
  left_join(dplyr::select(taxa, silva_Taxon, short.id), by = 'short.id') |>
  separate(silva_Taxon, into = c('d', 'p', 'c', 'o', 'f', 'g'), sep = ';') 

view(tbl)

coef_pred <- predict(f, type = 'coefficient') |> as.data.frame()
bind_cols(coef.f, coef_pred)
bind_cols(pred, Bm_scaled)

pred_manual <- rep(f$mu, 25) + scale(as.matrix(asv_predictors), center = T, scale = T) %*% coef.f$estimate
bind_cols(pred, pred_manual)

x_mx <- as.matrix(asv_predictors[, tbl$short.id]) %*% ß
fit <- lm(y ~ x, data = bind_cols(y = Bm_scaled,x = x_mx))
f$projection %*% coef(fit)[2:3]
filter(coef.f, estimate != 0)

tbl |>
  write_csv('rslt/taxa-weights-18s-spls.csv')
