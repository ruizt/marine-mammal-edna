library(tidyverse)
library(compositions)
library(spls)
load('data/ncog-18s.RData')

## DATA FILTERING AND AGGREGATION ----------------------------------------------

# dplyr::select columns of interest based on nonzero frequency across samples
cols_of_interest <- edna_samples %>%
  dplyr::select(starts_with('asv')) %>%
  as.matrix() %>%
  apply(2, function(.x){mean(.x > 0)}) %>%
  t() %>%
  as_tibble() %>%
  gather(col, prop.nz) %>%
  filter(prop.nz > 0.5) %>% # experiment?
  pull(col)

# edna_samples %>%
#   dplyr::select(all_of(cols_of_interest)) %>%
#   mutate(across(everything(), ~(.x > 0))) %>%
#   rowSums() %>%
#   bind_cols(dplyr::select(edna_samples, 1:5)) %>%
#   mutate(prop.asv.nz = `...1`/length(cols_of_interest)) %>%
#   arrange(prop.asv.nz) %>%
#   mutate(q = row_number()/nrow(edna_samples)) %>%
#   ggplot(aes(x = prop.asv.nz, y = q)) +
#   geom_step() +
#   labs(y = 'proportion of samples', x = 'proportion of nonzero reads') +
#   geom_vline(xintercept = 0.5)

# zero imputation (bayesian multiplicative, martin-fernandez 2015, luz calle 2019)
imputation_out <- edna_samples %>% 
  dplyr::select(all_of(cols_of_interest)) %>%
  zCompositions::cmultRepl(label = 0, 
                           method = 'GBM', 
                           output = 'prop',
                           z.warning = 1)

# remove rows with all zeroes
edna_imputed <- edna_samples %>%
  slice(-c(215, 432)) %>%
  dplyr::select(1:5) %>%
  bind_cols(imputation_out)

# aggregate, first across depth, then across line and station
edna_agg <- edna_imputed %>%
  group_by(cruise, line, sta) %>% # depth filtering or depth weighting?
  summarize(across(starts_with('asv'), ~mean(log(.x))), .groups = 'drop') %>%
  group_by(cruise) %>%
  summarize(across(starts_with('asv'), ~mean(.x)))

# interpretation: average proportion of asv.XXX across cruise is ZZZ
edna_data <- edna_agg %>% mutate(exp(pick(-cruise)))


## MODEL FITTING ---------------------------------------------------------------

# read taxa data
taxa <- read_tsv('data/NCOG_18sV9_asv_count_tax_S.tsv') %>%
  mutate(short.id = paste('asv', row_number(), sep = '.')) %>% 
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
asv_sightings <- left_join(edna_data_lc, scaled_sightings, by = c("cruise" = "cruiseID")) 
asv_sightings <- asv_sightings %>% drop_na(Bm_scaled) 
Bm_scaled <- asv_sightings$Bm_scaled

# dataset with only centered log transformed ASVs and blue whale density response
asv_predictors <- asv_sightings %>% select(starts_with("asv"))
asv_sightings <- cbind(Bm_scaled, asv_predictors)

# cross validate using K=2 to 8 latent components and eta = 0.1 to 0.9
set.seed(30424)
cv <- cv.spls(asv_predictors, Bm_scaled, eta = seq(0.1,0.9,0.1), K = c(1:8))
cv$eta.opt

# fit model with best eta
f <- spls(x = asv_predictors, y = asv_sightings$Bm_scaled, K = 2, eta = 0.4)
print(f)

# variance explained (rough)
pred <- predict.spls(f, type = 'fit')
1 - var(Bm_scaled - pred)/var(Bm_scaled)

# estimate and show the coefficients
coef.f <- coef(f) %>% as.data.frame() %>% rename(estimate = V1)
coef.f %>% filter(estimate != 0)

# add taxon info
tbl <- f$projection |>
  bind_cols(filter(coef.f, estimate != 0)) |>
  rename(reg.coef = estimate) |>
  rownames_to_column('short.id') |>
  left_join(dplyr::select(taxa, silva_Taxon, short.id), by = 'short.id') |>
  separate(silva_Taxon, into = c('d', 'p', 'c', 'o', 'f', 'g'), sep = ';') 

view(tbl)
