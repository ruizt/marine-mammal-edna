library(tidyverse)
library(magrittr)
library(modelr)
library(spls)

load('data/ncog-18s-processed-2024-07-11.RData')
load('data/ceta-density-processed-2024-07-11.RData')

# combine seasonally adjusted density estimates and seasonally adjusted edna data
whales <- inner_join(log_density_estimates_adj, edna_clr_adj, by = 'cruise')

## LEAVE ONE OUT CROSS VALIDATION ----------------------------------------------

# function to fit spls model and compute evaluation metrics
loocv_fit_fn <- function(.train, .test, .var, .parms){
  # fit model
  x.train <- dplyr::select(.train, starts_with('asv'))
  y.train <- pull(.train, {{.var}})
  fit <- spls(x.train, y.train, 
              K = .parms['ncomp'], eta = .parms['eta'], 
              scale.x = F, scale.y = F)
  
  # fit summaries
  y.train.hat <- predict(fit, type = 'fit')
  train.resid <- y.train - y.train.hat
  n <- length(y.train)
  p <- ncol(fit$projection)
  df <- nrow(fit$projection)
  rsq <- as.numeric(1 - (n - p)*var(train.resid)/((n - 1)*var(y.train)))
  sel.asv <- fit$projection |> rownames()
  
  # squared prediction error
  x.test <- dplyr::select(.test, starts_with('asv'))
  y.test <- pull(.test, {{.var}})
  y.test.hat <- predict(fit, newx = x.test, type = 'fit')
  test.resid <- y.test - y.test.hat
  test.spe <- as.numeric(test.resid^2)
  
  # outputs
  suffix <- deparse(substitute(.var))
  out <- tibble(model = list(fit),
                sel.asv = list(sel.asv),
                metrics = data_frame(df = df,
                                     rsq = rsq,
                                     p = y.test.hat[,1],
                                     pe = test.resid[,1],
                                     spe = test.spe) |>
                  list()
                ) |>
    rename_with(~paste0(.x, '.', suffix))
  return(out)
}

# hyperparameter grids
eta_grid_res <- 50
eta_grid <- rev(1 - seq(0.01, 0.95, length = eta_grid_res)^2)
ncomp_grid <- 1:10
data_partitions <- crossv_loo(whales)
obs_grid <- 1:nrow(data_partitions)
species_grid <- list("bm", "bp", "mn")
# new ncomp grid to find global maximum for bm
ncomp_grid_for_bm <- 11:13


# leave one out cross validation, export 9 RDS files per component number
# 3 files per species, per component number
lapply(ncomp_grid_for_bm, function(.ncomp){
  cv_out <- lapply(eta_grid, function(.eta){
    loo_out <- lapply(obs_grid, function(.obs){
      .train <- data_partitions$train[.obs][[1]] %>% as.data.frame() 
      .test <- data_partitions$test[.obs][[1]] %>% as.data.frame() 
      
      
      row_out <- bind_cols(obs = .obs,
                           ncomp = .ncomp,
                           eta = .eta,
                           loocv_fit_fn(.train, .test, bm, c(eta = .eta, ncomp = .ncomp)),
                           loocv_fit_fn(.train, .test, bp, c(eta = .eta, ncomp = .ncomp)),
                           loocv_fit_fn(.train, .test, mn, c(eta = .eta, ncomp = .ncomp)))
      
      return(row_out)
    }) %>% Reduce(bind_rows, .)
    
    paste('K = ', .ncomp, ', eta = ', round(.eta, 4), sep = '') |> print()
    return(loo_out)
    }) %>% Reduce(bind_rows, .)
  
  # output 3 files per species, per ncomp iteration
  for(species in species_grid){
    
  # Model Output file
  model_output <-  cv_out |> 
    select(obs, ncomp, eta, starts_with(paste("model.", species, sep = "")))
  
  model_file <- paste("rslt/comp/", today(), "/", species, "/cv-", .ncomp, "comp-model.rds", sep = "")
  
  # Selected ASVs Output File
  selected_asvs <- cv_out |> 
    select(obs, ncomp, eta, starts_with(paste("sel.asv.", species, sep = "")))
  
  asv_file <- paste("rslt/comp/", today(), "/", species, "/cv-", .ncomp, "comp-asvs.rds", sep = "")
  
  # Metrics Output File
  model_metrics <- cv_out |> 
    select(obs, ncomp, eta, starts_with(paste("metrics.", species, sep = "")))
  
  metric_file <- paste("rslt/comp/",today(), "/", species, "/cv-", .ncomp, "comp-metrics.rds", sep = "")
  
  
  #filename <- paste('rslt/comp/', today(), '/cv-', .ncomp, 'comp', '.rds', sep = '')
  #write_rds(cv_out, file = filename)
  write_rds(model_output, file = model_file)
  write_rds(selected_asvs, file = asv_file)
  write_rds(model_metrics, file = metric_file)
  #print("model file:", model_file)
  #print("asv file:", asv_file)
  #print("metric file:", metric_file)
  }
})


## find eta that maximizes cor(pred, obs) for each species

# BM ------------------------------------------------------------------
max_cor_bm <- lapply(1:13, function(.ncomp){
  paste('rslt/comp/2024-07-11/bm/cv-', .ncomp, "comp-metrics.rds", sep = "") |> 
    read_rds() |> 
    unnest_wider("metrics.bm") |> 
    group_by(ncomp, eta) |> 
    summarize(cor = cor(p + pe, p),
              mspe = mean(spe),
              df = mean(df),
              rsp = mean(rsq)) |> 
    group_by(ncomp) |> 
    slice_max(cor)
})
max_cor_bm <- bind_rows(max_cor_bm)
max_cor_bm

# find row where correlation is maximized, then join with selected asv file
max_bm = max(max_cor_bm$cor)

max_cor_bm |> 
  filter(cor == max_bm)

# no global maximum, so using local max at k = 5
# Read in asv file for k = 5 (local max corr)
bm_asvs_5comp <- read_rds("rslt/comp/2024-07-11/bm/cv-5comp-asvs.rds")


# Join files (filter only the asvs for 'best' model)
bm_selected_asvs <- max_cor_bm |> 
  filter(ncomp == 5) |> 
  left_join(bm_asvs_5comp, join_by(ncomp, eta))

# Find stability of each asv (proportion of loocv runs where it appears)
unnested_bm <- bm_selected_asvs |> 
  unnest(sel.asv.bm)

asv_counts_bm <- unnested_bm |> 
  group_by(sel.asv.bm) |> 
  summarise(n = n(),
            prop = n()/25,
            ncomp = ncomp,
            eta = eta) |> 
  rename(asv = sel.asv.bm) |> 
  distinct()

asv_counts_bm |> 
  filter(prop >= 0.4)

# save results as csv
save(asv_counts_bp, file = "rslt/comp/2024-07-11/bm/asv-stability-bm.csv")



# BP --------------------------------------------------
max_cor_bp <- lapply(ncomp_grid, function(.ncomp){
  paste('rslt/comp/2024-07-11/bp/cv-', .ncomp, "comp-metrics.rds", sep = "") |> 
    read_rds() |> 
    unnest_wider("metrics.bp") |> 
    group_by(ncomp, eta) |> 
    summarize(cor = cor(p + pe, p),
              mspe = mean(spe),
              df = mean(df),
              rsp = mean(rsq)) |> 
    group_by(ncomp) |> 
    slice_max(cor)
})
max_cor_bp <- bind_rows(max_cor_bp)
max_cor_bp

# Find row where correlation is maximized
max_bp = max(max_cor_bp$cor)

max_cor_bp |> 
  filter(cor == max_bp)

# Read in asv file for k = 6 (max corr)
bp_asvs_6comp <- read_rds("rslt/comp/2024-07-11/bp/cv-6comp-asvs.rds")


# Join files (filter only the asvs for 'best' model)
bp_selected_asvs <- max_cor_bp |> 
  filter(cor == max_bp) |> 
  left_join(bp_asvs_6comp, join_by(ncomp, eta))

# Find stability of each asv (proportion of loocv runs where it appears)
unnested_bp <- bp_selected_asvs |> 
  unnest(sel.asv.bp)

asv_counts_bp <- unnested_bp |> 
  group_by(sel.asv.bp) |> 
  summarise(n = n(),
            prop = n()/25,
            ncomp = ncomp,
            eta = eta) |> 
  rename(asv = sel.asv.bp) |> 
  distinct()
  
asv_counts_bp |> 
  filter(prop >= 0.4)

# save results as csv
save(asv_counts_bp, file = "rslt/comp/2024-07-11/bp/asv-stability-bp.csv")

# MN ------------------------------------------------------------------------
max_cor_mn <- lapply(ncomp_grid, function(.ncomp){
  paste('rslt/comp/2024-07-11/mn/cv-', .ncomp, "comp-metrics.rds", sep = "") |> 
    read_rds() |> 
    unnest_wider("metrics.mn") |> 
    group_by(ncomp, eta) |> 
    summarize(cor = cor(p + pe, p),
              mspe = mean(spe),
              df = mean(df),
              rsp = mean(rsq)) |> 
    group_by(ncomp) |> 
    slice_max(cor)
})
max_cor_mn <- bind_rows(max_cor_mn)
max_cor_mn

# find row where correlation is maximized, then join with selected asv file
max_mn = max(max_cor_mn$cor)

max_cor_mn |> 
  filter(cor == max_mn)


# Read in asv file for k = 6 (max corr)
mn_asvs_7comp <- read_rds("rslt/comp/2024-07-11/mn/cv-7comp-asvs.rds")

# Join files (filter only the asvs for 'best' model)
mn_selected_asvs <- max_cor_mn |> 
  filter(cor == max_mn) |> 
  left_join(mn_asvs_7comp, join_by(ncomp, eta))

mn_selected_asvs
# Find stability of each asv (proportion of loocv runs where it appears)
unnested_mn <- mn_selected_asvs |> 
  unnest(sel.asv.mn)

asv_counts_mn <- unnested_mn |> 
  group_by(sel.asv.mn) |> 
  summarise(n = n(),
            prop = n()/25,
            ncomp = ncomp,
            eta = eta) |> 
  rename(asv = sel.asv.mn) |> 
  distinct()

asv_counts_mn |> 
  filter(prop >= 0.9)

# save results as csv
save(asv_counts_mn, file = "rslt/comp/2024-07-11/mn/asv-stability-mn.csv")


#---------------------------------------------------------------------------------
# for each species and each K, find eta that maximizes corr(pred, obs)
max_cors <- lapply(ncomp_grid, function(.ncomp){
  paste('rslt/comp/2024-07-11/cv-', .ncomp, 'comp', '.rds', sep = '') |>
    readRDS() |>
    select(obs, ncomp, eta, starts_with('metrics')) |>
    pivot_longer(starts_with('metrics')) |>
    separate(name, into = c('drop', 'species')) |>
    select(-drop) |>
    unnest(value) |>
    group_by(ncomp, eta, species) |>
    summarize(cor = cor(p + pe, p),
              mspe = mean(spe),
              df = mean(df),
              rsq = mean(rsq)) |>
    group_by(ncomp, species) |>
    slice_max(cor)}) |>  
  Reduce(bind_rows, .)

write_rds(max_cors, 'rslt/comp/2024-07-06/max-cors.rds')
max_cors <- read_rds('rslt/comp/2024-07-06/max-cors.rds')

# inspect correlations as function of K
max_cors |>
  slice_max(eta) |>
  ggplot(aes(x = ncomp, y = cor)) +
  facet_wrap(~species) +
  geom_path() +
  scale_x_continuous(breaks = seq(1, 10, 2))
  
# inspect other metrics at K = 6 (fairly uniformly good)
max_cors |> filter(ncomp == 6)

# fix K = 6 and retrieve best etas
eta_sel <- max_cors |> filter(ncomp == 6) |> ungroup() |> select(species, eta)
cv6 <- readRDS('rslt/comp/2024-07-06/cv-6comp.rds')

# cross validation runs for K = 6 at best eta values
cv6_best <- cv6 |>
  rename_with(~gsub('sel.', '', .x)) |>
  pivot_longer(c(starts_with('metrics'), 
               starts_with('model'), 
               starts_with('asv'))) |>
  separate(name, into = c('var', 'species')) |>
  left_join(eta_sel, by = 'species') |>
  group_by(species) |>
  slice_min(abs(eta.x - eta.y)) |>
  select(-starts_with('eta'), -ncomp) |>
  pivot_wider(names_from = var, values_from = value) 

# remove full loocv from memory
rm('cv6')


## MODEL FITTING ---------------------------------------------------------------

# fit spls models
fit_fn <- function(.data, .var, .parms){
  x <- dplyr::select(.data, starts_with('asv'))
  y <- pull(.data, {{.var}})
  fit <- spls(x, y, 
              K = .parms['ncomp'], eta = .parms['eta'], 
              scale.x = F, scale.y = F)
}

eta_sel_vec <- pull(eta_sel, eta)

fit_bm <- fit_fn(whales, bm, c(ncomp = 6, eta = eta_sel_vec[1]))
fit_bp <- fit_fn(whales, bp, c(ncomp = 6, eta = eta_sel_vec[2]))
fit_mn <- fit_fn(whales, mn, c(ncomp = 6, eta = eta_sel_vec[3]))

# extract responses
x <- dplyr::select(whales, starts_with('asv'))
bm <- pull(whales, bm)
bp <- pull(whales, bp)
mn <- pull(whales, mn)
n <- nrow(whales)

# extract selected asvs
asv_bp <- tibble(short.id = colnames(x)[fit_bp$A],
                 coef = fit_bp$betahat[fit_bp$A]) |>
  left_join(dplyr::select(taxa, silva_Taxon, short.id), by = 'short.id') |>
  separate(silva_Taxon, into = c('d', 'p', 'c', 'o', 'f', 'g'), sep = ';') |>
  arrange(coef)

asv_bm <- tibble(short.id = colnames(x)[fit_bm$A],
                 coef = fit_bm$betahat[fit_bm$A]) |>
  left_join(dplyr::select(taxa, silva_Taxon, short.id), by = 'short.id') |>
  separate(silva_Taxon, into = c('d', 'p', 'c', 'o', 'f', 'g'), sep = ';') |>
  arrange(coef)

asv_mn <- tibble(short.id = colnames(x)[fit_mn$A],
                 coef = fit_mn$betahat[fit_mn$A]) |>
  left_join(dplyr::select(taxa, silva_Taxon, short.id), by = 'short.id') |>
  separate(silva_Taxon, into = c('d', 'p', 'c', 'o', 'f', 'g'), sep = ';') |>
  arrange(coef)

# selection frequency of each asv across loocv runs
sel_freqs <- cv6_best |>
  select(obs, species, asv) |>
  mutate(asv = map(asv, ~str_remove_all(.x, '\\.') |> str_flatten(' '))) |>
  unnest(asv) |>
  group_by(obs, species) |>
  tidytext::unnest_tokens(output = asv, input = asv, token = 'words') |>
  ungroup() |>
  count(species, asv) |>
  mutate(freq = n/25) 

# find selection frequencies for fitted models

# obs vs pred
cv6_best |>
  unnest(metrics) |>
  ggplot(aes(x = p + pe, y = p)) +
  facet_wrap(~species) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = 'observation', y = 'prediction')


# # join asvs across species
# asv_sel <- full_join(asv_bm, asv_bp, 
#           by = c('short.id', 'd', 'p', 'c', 'o', 'f', 'g'),
#           suffix = c('.bm', '.bp')) |>
#   full_join(asv_mn, 
#             by = c('short.id', 'd', 'p', 'c', 'o', 'f', 'g')) |>
#   rename(coef.mn = coef) 
# 
# # export selected asvs and fit summaries as excel workbook
# xl_out <- list(blue = mutate(asv_bm, coef = 100*(2^coef - 1)),
#                fin = mutate(asv_bp, coef = 100*(2^coef - 1)),
#                humpback = mutate(asv_mn, coef = 100*(2^coef - 1)),
#                fit.summary = fit_summary)
# 
# openxlsx::write.xlsx(xl_out, 
#            paste('rslt/tbl/fit-tables-', today(), '.xlsx', sep = ''), 
#            rownames = F)
# 
# 
# # extract fitted values
# fitted_bp <- predict(fit_bp)
# fitted_mn <- predict(fit_mn)
# fitted_bm <- predict(fit_bm)
# 
# # summarize prediction errors from loocv
# pred_err_df <- loopreds_sel |>
#   filter(metric == 'pe') |>
#   group_by(species) |>
#   summarize(bias.ratio = mean(value) |> exp(),
#             mspe = mean(value^2))
# 
# # summary of model fits
# fit_summary <- tibble(species = c('bp', 'bm', 'mn'),
#                       n.asv = c(nrow(fit_bp$projection),
#                                 nrow(fit_bm$projection),
#                                 nrow(fit_mn$projection)),
#                       var.expl = c(1 - ((n - 2)*var(bp - fitted_bp)/((n - 1)*var(bp))),
#                                    1 - ((n - 2)*var(bm - fitted_bm)/((n - 1)*var(bm))),
#                                    1 - ((n - 2)*var(mn - fitted_mn)/((n - 1)*var(mn))))) |>
#   left_join(pred_err_df)
# 
# fit_summary |>
#   xtable::xtable() |>
#   print() |>
#   clipr::write_clip()


# ## DRAFT MATERIAL --------------------------------------------------------------
# 
# cv_out |>
#   filter(metric %in% c('df', 'rsq', 'spe')) |>
#   group_by(eta, species, metric) |>
#   summarize(mean = mean(value),
#             sd = sd(value)) |>
#   # mutate(sd = if_else(metric == 'spe', NA, sd)) |>
#   ggplot(aes(x = eta, y = mean)) +
#   facet_wrap(~metric + species, scales = 'free_y', nrow = 3) +
#   geom_path() +
#   geom_ribbon(aes(ymin = mean - sd,
#                   ymax = mean + sd),
#               alpha = 0.3) + 
#   theme_bw(base_size = 16) +
#   geom_vline(aes(xintercept = eta.approx), 
#              data = tibble(species = c('bm', 'bp', 'mn'),
#                            eta.approx = c(0.65, 0.55, 0.63)),
#              linetype = 'dotdash') +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linewidth = 0.1, 
#                                         color = "#04004c"),
#         panel.background = element_rect(fill = "#fffffeff"),
#         plot.background = element_rect(fill = "#fffffeff"),
#         text = element_text(color = "#49485a")) +
#   scale_x_continuous(n.breaks = 4) +
#   labs(y = '')
# 
# ggsave(filename = 'rslt/plots/51024/tuning.png',
#        height = 5, width = 7, dpi = 400)
# 
# 
# count(group_by(asv_bm, d), name = 'bm') |>
#   full_join(count(group_by(asv_bp, d), name = 'bp'), by = 'd') |>
#   full_join(count(group_by(asv_mn, d), name = 'mn'), by = 'd') |>
#   pivot_longer(-d) |>
#   mutate(phylum = str_remove(d, 'd__') |> str_trim()) |>
#   drop_na() |>
#   ggplot(aes(x = name, fill = phylum)) +
#   geom_bar(position = 'fill') +
#   scale_y_continuous(labels = scales::percent)
# 
# 
# asv_sel |>
#   pivot_longer(starts_with('coef')) |>
#   mutate(species = str_remove(name, 'coef.')) |>
#   drop_na(value) |>
#   group_by(species) |>
#   summarize(across(p:g, ~length(unique(.x)))) |>
#   write.csv('rslt/tbl/classification-counts.csv')
# 
# # # verify fitted model components
# # x_sel <- dplyr::select(x, rownames(fit$projection))
# # x_sel_proj <- as.matrix(x_sel) %*% fit$projection
# # 
# # lm(y ~ x_sel_proj)$fitted |> 
# #   bind_cols(fitted_vals)
# # 
# # (fit$projection %*% lm(y ~ x_sel_proj)$coef[2:3]) |>
# #   bind_cols(fit$betahat[fit$betahat != 0])
# 
# # visualize
# pal <- colorRampPalette(c('#649dfa', '#02204f'), bias = 1)
# 
# asv_sel |>
#   mutate(phylum = str_remove(p, 'p__') |> str_trim() |> replace_na('Unknown')) |> 
#   dplyr::select(short.id, phylum, starts_with('coef')) |>
#   rename_with(~str_remove(.x, 'coef.')) |>
#   rowwise() |>
#   mutate(n.na = sum(is.na(c(bm, bp, mn)))) |>
#   filter(n.na < 2, 
#          phylum != 'Unknown', 
#          phylum != 'uncultured') |>
#   pivot_longer(bm:mn, names_to = 'species', values_to = 'coef') |>
#   mutate(species = factor(species, levels = c('bm', 'bp', 'mn'), labels = c('blue', 'fin','humpback'))) |>
#   group_by(phylum) |>
#   mutate(count = n()) |>
#   ungroup() |>
#   mutate(phylum = fct_reorder(phylum, count),
#          coef.trans = 100*(2^coef - 1),
#          p.ix = as.numeric(phylum) %% 2) |>
#   ggplot(aes(y = phylum)) +
#   geom_col(aes(x = coef.trans, fill = factor(p.ix)),
#            position = position_dodge2(preserve = 'single')) +
#   facet_wrap(~species) +
#   scale_fill_manual(values = c("#04004c", "#8e87ea")) +
#   guides(fill = guide_none()) +
#   scale_x_continuous(n.breaks = 4) +
#   theme_bw(base_size = 16) +
#   geom_vline(xintercept = 0, linewidth = 0.3, color = "#04004c") +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major.y = element_blank(),
#         panel.grid.major.x = element_line(linewidth = 0.1, 
#                                           color = "#04004c"),
#         panel.background = element_rect(fill = "#fffffeff"),
#         plot.background = element_rect(fill = "#fffffeff"),
#         text = element_text(color = "#49485a")) +
#   labs(x = 'percent change in median scaled sighting', 
#        y = '')
# 
# ggsave(filename = 'rslt/plots/51024/asv-sel-all.png',
#        height = 4, width = 6, dpi = 400)
# 
# 
# asv_mn |>
#   mutate(phylum = str_remove(p, 'p__') |> str_trim() |> replace_na('Unknown')) |> 
#   mutate(phylum = fct_reorder(phylum, coef, ~max(abs(.x))),
#          x = 100*(2^coef - 1)) |>
#   ggplot(aes(y = phylum)) +
#   geom_col(aes(x = x, fill = phylum),
#            width = 1,
#            position = position_dodge2(preserve = 'single',
#                                       padding = 0.1)) +
#   scale_fill_manual(values = pal(19)) +
#   guides(fill = guide_none()) +
#   scale_x_continuous(n.breaks = 4) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         panel.grid.major.x = element_line(linewidth = 0.1, 
#                                           color = 'black'),
#         panel.grid.major.y = element_line(linewidth = 0.1, 
#                                           color = 'grey')) +
#   labs(x = 'percent change in density', 
#        y = '',
#        title = 'humpback')
# 
# ggsave(filename = 'rslt/plots/42524/asv-sel-mn.png',
#        height = 4, width = 5, dpi = 300)
# 
# 
# asv_bm |>
#   mutate(phylum = str_remove(p, 'p__') |> str_trim() |> replace_na('Unknown')) |> 
#   mutate(phylum = fct_reorder(phylum, coef, ~max(abs(.x))),
#          x = 100*(2^coef - 1)) |>
#   ggplot(aes(y = phylum)) +
#   geom_col(aes(x = x, fill = phylum),
#            width = 1,
#            position = position_dodge2(preserve = 'single',
#                                       padding = 0.1)) +
#   scale_fill_manual(values = pal(19)) +
#   guides(fill = guide_none()) +
#   scale_x_continuous(n.breaks = 4) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         panel.grid.major.x = element_line(linewidth = 0.1, 
#                                           color = 'black'),
#         panel.grid.major.y = element_line(linewidth = 0.1, 
#                                           color = 'grey')) +
#   labs(x = 'percent change in density', 
#        y = '',
#        title = 'blue')
# 
# ggsave(filename = 'rslt/plots/42524/asv-sel-bm.png',
#        height = 4, width = 5, dpi = 300)
# 
# asv_bp |>
#   mutate(phylum = str_remove(p, 'p__') |> str_trim() |> replace_na('Unknown')) |> 
#   mutate(phylum = fct_reorder(phylum, coef, ~max(abs(.x))),
#          x = 100*(2^coef - 1)) |>
#   ggplot(aes(y = phylum)) +
#   geom_col(aes(x = x, fill = phylum),
#            width = 1,
#            position = position_dodge2(preserve = 'single',
#                                       padding = 0.1)) +
#   scale_fill_manual(values = pal(19)) +
#   guides(fill = guide_none()) +
#   scale_x_continuous(n.breaks = 4) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         panel.grid.major.x = element_line(linewidth = 0.1, 
#                                           color = 'black'),
#         panel.grid.major.y = element_line(linewidth = 0.1, 
#                                           color = 'grey')) +
#   labs(x = 'percent change in density', 
#        y = '',
#        title = 'fin')
# 
# ggsave(filename = 'rslt/plots/42524/asv-sel-bp.png',
#        height = 4, width = 5, dpi = 300)
# 
# p1 <- loopreds_sel |>
#   filter(metric == 'p') |>
#   bind_cols(obs = c(bm, bp, mn)) |>
#   ggplot(aes(x = obs, y = value)) +
#   geom_point() +
#   facet_wrap(~species) +
#   geom_abline(slope = 1, intercept = 0) +
#   theme_bw() +
#   labs(x = '',
#        y = 'prediction')
# 
# tibble(fitted = c(fitted_bm, fitted_bp, fitted_mn),
#        observed = c(bm, bp, mn),
#        species = rep(c('blue', 'fin', 'humpback'), 
#                      times = c(length(bm), length(bp), length(mn)))) |>
#   ggplot(aes(x = observed, y = fitted)) +
#   facet_wrap(~species) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) +
#   theme_bw(base_size = 18) +
#   labs(x = 'scaled sighting logratio',
#        y = 'fitted value') +
#   theme(panel.grid.minor = element_blank(),
#         panel.grid.major = element_line(linewidth = 0.1, 
#                                         color = 'black'),
#         panel.background = element_rect(fill = "#fffffeff"),
#         plot.background = element_rect(fill = "#fffffeff"),
#         text = element_text(color = "#49485a"))
# 
# ggsave(filename = 'rslt/plots/51024/obs-fitted.png',
#        height = 2.5, width = 6, dpi = 400)
# 
# 
# library(patchwork)
# 
# fig <- p1 + p2 + plot_layout(nrow = 2)
# 
# ggsave(filename = 'rslt/plots/42524/pred-fitted.png',
#        height = 4, width = 5, dpi = 300)
