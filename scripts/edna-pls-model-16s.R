###############################################################################
#
# Partial Least Squares (PLS) on 16s Data to Model Blue Whale Densities 
#
###############################################################################

# Set Up
###############################################################################
library(compositions)
library(pls)
library(spls)
library(tidyverse)
library(caret)

# load edna_data using 16s 
load("data/edna_16s_processed.RData")
# ls() to see file names

# read in taxa data
edna_in <- read_tsv('data/NCOG_16s_asv_count_tax_S.tsv') %>%
  mutate(short.id = paste('asv', row_number(), sep = '.'))
taxa <- edna_in %>% 
  dplyr::select(where(is.character), silva_Confidence)

# compute closure using centered log-ratio transform
cruises <- edna_data %>% select(cruise)
# transform log contrast
edna_data_lc <- edna_data %>% select(-cruise) %>% 
  clr() %>% 
  clrInv() %>% as.data.frame() 
# combine cruise and log ratio-transformed compositions
edna_data_lc <- cbind(cruises, edna_data_lc)

# get response variable: density of blue whales `Bm_scaled`
scaled_sightings <- read_csv("data/CC_on_effort_scaled_sightings_Ceta_SCB.csv")
scaled_sightings$cruiseID <- str_replace(scaled_sightings$cruiseID, "CC", "X20")
asv_sightings <- left_join(edna_data_lc, scaled_sightings, by = c("cruise" = "cruiseID")) 
asv_sightings <- asv_sightings %>% drop_na(Bm_scaled) 
Bm_scaled <- asv_sightings$Bm_scaled

# dataset with only centered log transformed ASVs and blue whale density response
asv_predictors <- asv_sightings %>% select(starts_with("asv"))
asv_sightings <- cbind(Bm_scaled, asv_predictors)

# PLS Model
###############################################################################
# Run PLS model
model <- plsr(Bm_scaled ~ ., ncomp = 5, data=asv_sightings, validation = "CV")
summary(model) # note: 5 components captures 90.05% of variance

# extract component weights as a dataframe
l <- model$loadings 
loading_df <- data.frame(matrix(as.numeric(l), attributes(l)$dim, dimnames=attributes(l)$dimnames))
loading_df$asv = rownames(loading_df)

# add in estimates for 5 components
loading_df$coef_estimates <- model$coefficients %>% as.data.frame() %>% select(5) %>% pull()

# calculate single weight for each ASV and combine with taxa info
weight_df <- loading_df  %>%
  rowwise() %>% 
  # for each component, max absolute loading as a weight for variable importance score
  mutate(weight = max(abs(c_across(Comp.1:Comp.5)))) %>%
  select(weight, coef_estimates) %>%
  bind_cols(short.id = loading_df$asv) %>%
  # incorporate taxa data
  left_join(taxa, by = "short.id") %>%
  select(short.id, silva_Taxon, weight, coef_estimates)

# format taxa info in dataframe
clean_weight_df <- weight_df %>% 
  separate(silva_Taxon, into = c('d', 'p', 'c', 'o', 'f', 'g'), sep = ';') 

# write weight/estimates df to a CSV file
clean_weight_df %>% 
  select(-short.id) %>% 
  write_csv(file = 'rslt/taxa-weights-16s.csv')

# Exploration of Model & Weights
###############################################################################

# average weight is 0.0495
clean_weight_df %>% summary()

# top 50 weighted ASVs
clean_weight_df %>% 
  arrange(desc(weight)) %>% 
  head(50) %>% 
  ggplot(aes(x = d, y = weight)) +
  geom_bar(aes(fill = short.id),
           position = "dodge",
           stat = "identity") + guides(fill="none") +
  theme_bw() 

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
abline(h=0, col="red")
plot(model)


