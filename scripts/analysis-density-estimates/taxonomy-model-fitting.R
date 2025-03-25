library(tidyverse)
library(spls)
library(readxl)
library(openxlsx)

# Data processing --------------------------------------------------------------

# directories
model_dir <- 'rslt/models/density/'
data_dir <- 'data/processed/'


# load in documented relationships table
doc.rel <- read.csv(here::here("data/whale-edna-relationships.csv"))

doc.rel <- doc.rel |> 
  mutate(Whale.species = str_replace(Whale.species, " whales", "s"))


# candidate asvs
# candidate asvs from 18sv9
paste(data_dir, 'ncog18sv9.RData', sep = '') |> load()
asv_taxa_18sv9 <- asv_taxa 

# candidate asvs from 18sv4
paste(data_dir, 'ncog18sv4.RData', sep = '') |> load()
asv_taxa_18sv4 <- asv_taxa

# candidate asvs from 16s
paste(data_dir, 'ncog16s.RData', sep = '') |> load()
asv_taxa_16s <- asv_taxa


# join candidate asv tables
asv_taxa_16s <- asv_taxa_16s |> 
  mutate(gene.seq = "16s") 

asv_taxa_18sv4 <- asv_taxa_18sv4|> 
  mutate(gene.seq = "18sv4")

asv_taxa_18sv9 <- asv_taxa_18sv9|> 
  mutate(gene.seq = "18sv9")

asv_taxa_all <- rbind(asv_taxa_16s,
                      asv_taxa_18sv4,
                      asv_taxa_18sv9)

asv_taxa_all


# 18sv9
paste(model_dir, 'fitted-models-18sv9.RData', sep = '') |> load()
# selected asvs
ss_asvs_18sv9 <- fitted_models |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback')),
         marker = '18SV9',
         coef.ss = map(coef, ~2^.x)) |>
  select(species, coef, coef.ss, ss) |>
  unnest(everything()) |>
  rename(asv = ss,
         coef.lr = coef) |>
  left_join(asv_taxa_18sv9, join_by(asv == short.id)) |>
  select(-feature.id) |>
  arrange(species, desc(coef.lr))

# 18sv4
paste(model_dir, 'fitted-models-18sv4.RData', sep = '') |> load()

# selected asvs
ss_asvs_18sv4 <- fitted_models |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback')),
         marker = '18SV4',
         coef.ss = map(coef, ~2^.x)) |>
  select(species, coef, coef.ss, ss) |>
  unnest(everything()) |>
  rename(asv = ss,
         coef.lr = coef) |>
  left_join(asv_taxa_18sv4, join_by(asv == short.id)) |>
  select(-feature.id) |>
  arrange(species, desc(coef.lr))

# 16s
paste(model_dir, 'fitted-models-16s.RData', sep = '') |> load()

# selected asvs
ss_asvs_16s <- fitted_models |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback')),
         marker = '16S',
         coef.ss = map(coef, ~2^.x)) |>
  select(species, coef, coef.ss, ss) |>
  unnest(everything()) |>
  rename(asv = ss,
         coef.lr = coef) |>
  left_join(asv_taxa_16s, join_by(asv == short.id)) |>
  select(-feature.id) |>
  arrange(species, desc(coef.lr))


# join selected asv tables
ss_asvs_16s <- ss_asvs_16s |> 
  mutate(gene.seq = "16s")

ss_asvs_18sv4 <- ss_asvs_18sv4 |> 
  mutate(gene.seq = "18sv4")

ss_asvs_18sv9 <- ss_asvs_18sv9 |> 
  mutate(gene.seq = "18sv9")

ss_asvs_all <- rbind(ss_asvs_16s,
                     ss_asvs_18sv4,
                     ss_asvs_18sv9)

ss_asvs_all <- ss_asvs_all |> 
  mutate(species = paste(tolower(species), "s", sep = ""))
# ------------------------------------------------------------------------------

## MODEL: STABLE SET ASV OVERLAP WITH DOCUMENTED RELATIONSHIPS

# ------------------------------------------------------------------------------

# all of the asvs 
ss_asvs_all |> 
  inner_join(doc.rel, join_by(species == Whale.species,
                              p == Phylum,
                              c == Class)) |> 
  select(gene.seq,asv,
         species,
         Connection,
         p,c,o,f,g) |> 
  distinct() #|> group_by(species) |> count()


# 16s ASVs to fit the model on
asvs_16s <- ss_asvs_all |> 
  filter(gene.seq == "16s") |> 
  inner_join(doc.rel, join_by(species == Whale.species,
                              p == Phylum,
                              c == Class)) |> 
  select(gene.seq,asv,
         species,
         Connection,
         p,c,o,f,g) |> 
  mutate(species = case_when(species == "blues" ~ "bm",
                             species == "humpbacks" ~ "mn")) |> 
  distinct()


# PLS MODELS -------------------------------------------------------------------

# NOTE: asv.30001 from 18sv9 corresponds to the same species as asv.15277 from 16s
# therefore, it was omitted from the model (only 16s data is being used)

#-------------------------------------------------------------------------------

# load edna data
data_dir <- 'data/processed/'
#16s
paste(data_dir, 'ncog16s.RData', sep = '') |> load()
edna_16s <- edna

paste(data_dir, 'density-estimates.RData', sep = '') |> load() 

# join sightings, and 16s 
whales <- inner_join(dens, edna_16s, by = 'cruise') 


# MODEL FITTING GRID SEARCH FOR NCOMP
species_grid <- c("bm","mn") # blue whales and humpbacks only
ncomp_grid <- 1:9
n <- 25

ss_taxonomy_results <- data.frame(
  Marker = character(),
  Model = character(),
  Species = character(),
  p = integer(),
  Ncomp = integer(),
  R2 = numeric(),
  Adj_R2 = numeric(),
  MSPE = numeric()
)

for(ncomp in ncomp_grid){
  for(spec in species_grid){
    # select relevant features
    species_asvs <- asvs_16s |> 
      filter(species == spec)
    
    print(species_asvs)
    
    p = nrow(species_asvs)
    
    x <- whales |> 
      select(all_of(species_asvs$asv)) 
    # target column
    y <- whales |> select(spec) |> pull()
    
    # fit model
    fit <- spls(x, y , eta = 0, K = ncomp)
    
    # r2 and adj r2
    fitted <- predict(fit, type = 'fit')
    resid <- y - fitted
    r2 <- 1 - (var(resid)/var(y))
    adj_r2 <- 1 - ((1-r2)*(n-1))/(n-ncomp-1)
    
    #loocv for mspe (this is why spls package is being used)
    bp_cv_results <- cv.spls(x, y, fold = 25, eta = 0, K = ncomp, plot.it = F)
    
    mspe <- bp_cv_results$mspemat[[1]]
    
    new_row <- data.frame(
      Marker = "16s",
      Model = "pls",
      Species = spec,
      P = p,
      Ncomp = ncomp,
      R2 = r2,
      Adj_R2 = adj_r2,
      MSPE = mspe
    )
    
    print(new_row)
    
    ss_taxonomy_results <- rbind(ss_taxonomy_results,new_row)
  }
}

# fix column names 
colnames(ss_taxonomy_results) <- c("marker", "model", "species", "p", "ncomp", "r2", "adj.r2", "mspe")

best_ss_models <- ss_taxonomy_results |> 
  group_by(species) |> 
  slice_min(mspe)


# LM MODELS --------------------------------------------------------------------

# BM MODEL
bm_16s_asvs_ss <- asvs_16s |> 
  filter(species == "bm")

bm_16s_ss <- whales |> 
  select(bm, all_of(bm_16s_asvs_ss$asv))



bm.fit <- lm(bm ~ ., data = bm_16s_ss)

bm_summary <- summary(bm.fit)

r2 <- bm_summary$r.squared
adj.r2 <- bm_summary$adj.r.squared

# loocv for mspe 
squared_errors <- numeric(nrow(whales))
for (i in 1:nrow(bm_16s_ss)) {
  loo <- bm_16s_ss[-i, ]
  
  bm.fit <- lm(bm ~ ., data = loo)
  
  prediction <- predict(bm.fit, newdata = bm_16s_ss[i, , drop = FALSE])
  
  squared_errors[i] <- (bm_16s_ss$bm[i] - prediction)^2
}

mspe <- mean(squared_errors)

mspe

bm_ss <- data.frame(marker = "16s", model = "lm",species = "bm", p = nrow(bm_16s_asvs_ss), ncomp = NA, r2 = r2, adj.r2 = adj.r2, mspe = mspe)
bm_ss

# MN MODEL 
mn_16s_asvs_ss <- asvs_16s |> 
  filter(species == "mn")

mn_16s_ss <- whales |> 
  select(mn, all_of(mn_16s_asvs_ss$asv))

mn.fit <- lm(mn ~ ., data = mn_16s_ss)

mn_summary <- summary(mn.fit)

r2 <- mn_summary$r.squared
adj.r2 <- mn_summary$adj.r.squared

# loocv for mspe 
squared_errors <- numeric(nrow(whales))
for (i in 1:nrow(mn_16s_ss)) {
  loo <- mn_16s_ss[-i, ]
  
  mn.fit <- lm(mn ~ ., data = loo)
  
  prediction <- predict(mn.fit, newdata = mn_16s_ss[i, , drop = FALSE])
  
  squared_errors[i] <- (mn_16s_ss$mn[i] - prediction)^2
}

mspe <- mean(squared_errors)

mspe

mn_ss <- data.frame(marker = "16s", model = "lm",species = "mn", p = nrow(mn_16s_asvs_ss), ncomp = NA, r2 = r2, adj.r2 = adj.r2, mspe = mspe)
mn_ss
# ------------------------------------------------------------------------------





# ------------------------------------------------------------------------------

## CANDIDATE ASV OVERLAP WITH DOCUMENTED RELATIONSHIPS MODELS

# ------------------------------------------------------------------------------

# all of the asvs (ORDER Level => 61 asvs)
asv_taxa_all |> 
  inner_join(doc.rel, join_by(p == Phylum,
                              c == Class,
                              o == Order)) |> 
  select(gene.seq,short.id,
         Whale.species,
         Connection,
         p,c,o,f,g) |> 
  distinct() |> group_by(Whale.species, gene.seq) |> count()


# 16s ASVs to fit the model on
asvs_16s <- asv_taxa_all |> 
  filter(gene.seq == "16s") |> 
  inner_join(doc.rel, join_by(p == Phylum,
                              c == Class,
                              o == Order)) |> 
  select(gene.seq,short.id,
         Whale.species,
         Connection,
         p,c,o,f,g) |> 
  mutate(Whale.species = case_when(Whale.species == "blues" ~ "bm",
                                   Whale.species == "humpbacks" ~ "mn")) |> 
  distinct() 

# 18sv9 ASVs to fit the model on
asvs_18sv9 <- asv_taxa_all |> 
  filter(gene.seq == "18sv9") |> 
  inner_join(doc.rel, join_by(p == Phylum,
                              c == Class,
                              o == Order)) |> 
  select(gene.seq,short.id,
         Whale.species,
         Connection,
         p,c,o,f,g) |> 
  distinct()

# 18sv4 ASVs to fit models on 
asvs_18sv4 <- asv_taxa_all |> 
  filter(gene.seq == "18sv4") |> 
  inner_join(doc.rel, join_by(p == Phylum,
                              c == Class,
                              o == Order)) |> 
  select(gene.seq,short.id,
         Whale.species,
         Connection,
         p,c,o,f,g) |> 
  distinct()



## PLS MODELS FOR 16s DATA ----------------------------------------------------

species_grid <- c("bm","mn") # blue whales and humpbacks only
ncomp_grid <- 1:10
n <- 25

candidate_taxonomy_results <- data.frame(
  Marker = character(),
  model = character(),
  Species = character(),
  P = integer(),
  Ncomp = integer(),
  R2 = numeric(),
  Adj_R2 = numeric(),
  MSPE = numeric()
)

for(ncomp in ncomp_grid){
  for(species in species_grid){
    # select relevant features
    species_asvs <- asvs_16s |> 
      filter(Whale.species == species)
    
    p <- species_asvs |> count()
    
    x <- whales |> 
      select(all_of(species_asvs$short.id)) 
    # target column
    y <- whales |> select(species) |> pull()
    
    # fit model
    fit <- spls(x, y , eta = 0, K = ncomp)
    
    # r2 and adj r2
    fitted <- predict(fit, type = 'fit')
    resid <- y - fitted
    r2 <- 1 - (var(resid)/var(y))
    adj_r2 <- 1 - ((1-r2)*(n-1))/(n-ncomp-1)
    
    #loocv for mspe (this is why spls package is being used)
    bp_cv_results <- cv.spls(x, y, fold = 25, eta = 0, K = ncomp, plot.it = F)
    
    mspe <- bp_cv_results$mspemat[[1]]
    
    new_row <- data.frame(
      Marker = "16s",
      Model = "pls",
      Species = species,
      P = p,
      Ncomp = ncomp,
      R2 = r2,
      Adj_R2 = adj_r2,
      MSPE = mspe
    )
    
    print(new_row)
    
    candidate_taxonomy_results <- rbind(candidate_taxonomy_results,new_row)
  }
}

# fix column names 
colnames(candidate_taxonomy_results) <- c("marker",
                                          "model",
                                          "species",
                                          "p",
                                          "ncomp",
                                          "r2",
                                          "adj.r2",
                                          "mspe")

best_16s_models <- candidate_taxonomy_results |> 
  group_by(species) |> 
  slice_min(mspe)

# ------------------------------------------------------------------------------


# 18sv9 models - lm() call -----------------------------------------------------


# load in 18sv9 data
paste(data_dir, 'ncog18sv9.RData', sep = '') |> load()
edna_18sv9 <- edna
# join asvs with whale sightings 
whales <- inner_join(dens, edna_18sv9, by = "cruise")
whales


## BM MODEL --------------------------------------------------------------------
bm_18sv9_asvs <- asvs_18sv9 |> 
  filter(Whale.species == "blues")

bm_18sv9 <- whales |> 
  select(bm, all_of(bm_18sv9_asvs$short.id))

whales

bm.fit <- lm(bm ~ ., data = bm_18sv9)

bm_summary <- summary(bm.fit)

r2 <- bm_summary$r.squared
adj.r2 <- bm_summary$adj.r.squared

# loocv for mspe 
squared_errors <- numeric(nrow(whales))
for (i in 1:nrow(bm_18sv9)) {
  loo <- bm_18sv9[-i, ]
  
  bm.fit <- lm(bm ~ ., data = loo)
  
  prediction <- predict(bm.fit, newdata = bm_18sv9[i, , drop = FALSE])
  
  squared_errors[i] <- (bm_18sv9$bm[i] - prediction)^2
}

mspe <- mean(squared_errors)

mspe

row1 <- data.frame(marker = "18sv9", model = "lm",species = "bm", p = nrow(bm_18sv9_asvs), ncomp = NA, r2 = r2, adj.r2 = adj.r2, mspe = mspe)


# MN MODEL ---------------------------------------------------------------------

mn_18sv9_asvs <- asvs_18sv9 |> 
  filter(Whale.species == "humpbacks")

mn_18sv9 <- whales |> 
  select(mn, all_of(mn_18sv9_asvs$short.id))

whales

mn.fit <- lm(mn ~ ., data = mn_18sv9)

mn_summary <- summary(mn.fit)

r2 <- mn_summary$r.squared
adj.r2 <- mn_summary$adj.r.squared

# loocv for mspe 
squared_errors <- numeric(nrow(whales))
for (i in 1:nrow(mn_18sv9)) {
  loo <- mn_18sv9[-i, ]
  
  mn.fit <- lm(mn ~ ., data = loo)
  
  prediction <- predict(mn.fit, newdata = mn_18sv9[i, , drop = FALSE])
  
  squared_errors[i] <- (mn_18sv9$mn[i] - prediction)^2
}

mspe <- mean(squared_errors)

mspe

row2 <- data.frame(marker = "18sv9",model = "lm", species = "mn", p = nrow(mn_18sv9_asvs),
                   ncomp = NA, r2 = r2, adj.r2 = adj.r2, mspe = mspe)


# 18sv4 model - ONLY HAS MATCHES FOR MN (Humpback Whales) ----------------------
# load in 18sv9 data
paste(data_dir, 'ncog18sv4.RData', sep = '') |> load()
edna_18sv4 <- edna
# join asvs with whale sightings 
whales <- inner_join(dens, edna_18sv4, by = "cruise")
whales

# matches
asvs_18sv4

# both asvs appear to be linked to the same species

mn_18sv4_asvs <- asvs_18sv4 

mn_18sv4 <- whales |> 
  select(mn, asv.1868) # keeping just one asv reduces mspe and increases adj.r2

mn_18sv4

mn.fit <- lm(mn ~ ., data = mn_18sv4)

mn_summary <- summary(mn.fit)

r2 <- mn_summary$r.squared
adj.r2 <- mn_summary$adj.r.squared

# loocv for mspe 
squared_errors <- numeric(nrow(whales))
for (i in 1:nrow(mn_18sv4)) {
  loo <- mn_18sv4[-i, ]
  
  mn.fit <- lm(mn ~ ., data = loo)
  
  prediction <- predict(mn.fit, newdata = mn_18sv4[i, , drop = FALSE])
  
  squared_errors[i] <- (mn_18sv4$mn[i] - prediction)^2
}

mspe <- mean(squared_errors)

mspe

row3 <- data.frame(marker = "18sv4", model = "lm", species = "mn", p = 1,
                   ncomp = NA, r2 = r2, adj.r2 = adj.r2, mspe = mspe)



# RESULTS ----------------------------------------------------------------------

# candidate asv models
candidate_asv_results <- rbind(best_16s_models,row1,row2,row3,)

candidate_asv_results


# stable set asv models

ss_asv_results <- rbind(best_ss_models, bm_ss, mn_ss)

ss_asv_results


# Save results for later
saveRDS(candidate_asv_results, "rslt/tbl/candidate-asv-model-res.rds")
saveRDS(ss_asv_results, "rslt/tbl/stable-set-asv-model-res.rds")
