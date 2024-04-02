library(tidyverse)
library(broom)
load('data/edna_16s_processed.RData')

# read in mammal data
mammal_raw <- read_csv('data/CC_on_effort_scaled_sightings_Ceta_SCB.csv')

# recode cruise IDs and select species of interest
mammal <- mammal_raw %>%
  mutate(cruise = str_replace(cruiseID, 'CC', 'X20')) %>%
  select(cruise, Bp_scaled, Bm_scaled, Mn_scaled) %>%
  rename(bp = Bp_scaled,
         bm = Bm_scaled,
         mn = Mn_scaled)

# project edna onto a few principal components
fit <- edna_data %>%
  mutate(across(starts_with('asv'), log)) %>%
  select(starts_with('asv')) %>%
  prcomp(scale = T)

plot(1:27, fit$sdev)

pc_data <- tidy(fit) %>%
  filter(PC <= 6) %>%
  pivot_wider(names_from = PC, values_from = value, id_cols = row, names_prefix = 'PC') %>%
  bind_cols(cruise = pull(edna_data, cruise)) %>%
  right_join(mammal, by = 'cruise') %>%
  select(-row)

# plot pcs against densities -- not very strong marginal relationships
pc_data %>%
  pivot_longer(cols = starts_with("PC"), names_to = 'component', values_to = 'score') %>%
  ggplot(aes(x = score, y = bp + 1)) +
  facet_wrap(~component, scales = 'free_x') +
  geom_point() +
  scale_y_log10()

# about 30-60% variance explained
fit_bp <- lm(bp ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data = pc_data) 
fit_bm <- lm(bm ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data = pc_data)
fit_mn <- lm(mn ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data = pc_data)

augment(fit_bp) %>%
  ggplot(aes(x = .fitted, y = bp)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

augment(fit_bm) %>%
  ggplot(aes(x = .fitted, y = bm)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)

augment(fit_mn) %>%
  ggplot(aes(x = .fitted, y = mn)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0)


# presence/absence

mammal %>%
  pivot_longer(cols = c(bp, bm, mn)) %>%
  group_by(name) %>%
  summarize(cruises_present = sum(value > 0, na.rm = T),
            cruises_absent = sum(value == 0, na.rm = T),
            cruises_missing = sum(is.na(value)),
            n = n())

fit_pa_bp <- glm(I(bp > 0) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data = pc_data, family = binomial()) 
summary(fit_pa_bm)
fit_pa_mn <- glm(I(mn > 0) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data = pc_data, family = binomial()) 
summary(fit_pa_mn)
fit_pa_bm <- glm(I(bm > 0) ~ PC1 + PC2 + PC3 + PC4 + PC5 + PC6, data = pc_data, family = binomial()) 
summary(fit_pa_mn)

pred <- predict(fit_pa_bm, type = 'response') > 0.5
obs <- fit_pa_bm$y

table(pred, obs)

##

# read in 16s reads
edna_in <- read_tsv('data/NCOG_16S_asv_count_tax_S.tsv') %>%
  mutate(short.id = paste('asv', row_number(), sep = '.'))
taxa <- edna_in %>% 
  dplyr::select(where(is.character), silva_Confidence)

sel_asv <- colnames(edna_data)

weight_df <- as_tibble(fit$rotation) %>%
  select(1:6) %>%
  rowwise() %>% 
  mutate(weight = max(abs(c_across(PC1:PC6)))) %>%
  # arrange(desc(weight)) %>%
  select(weight) %>%
  bind_cols(short.id = sel_asv[-1]) %>%
  left_join(taxa, by = "short.id") %>%
  select(silva_Taxon, weight)

weight_df %>%
  separate(silva_Taxon, into = c('d', 'p', 'c', 'o', 'f', 'g'), sep = ';') %>%
  write_csv(file = 'taxa-weights.csv')
