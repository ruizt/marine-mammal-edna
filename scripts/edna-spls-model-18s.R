library(compositions)
library(pls)
library(spls)
library(tidyverse)
library(RColorBrewer)


## DATA PREPARATION ------------------------------------------------------------

# load aggregated 18s relative abundances
load("data/edna_18s_processed.RData")
head(edna_data)

# read in taxonomic classifications
taxa <- read_tsv('data/NCOG_18sV9_asv_count_tax_S.tsv') |>
  mutate(short.id = paste('asv', row_number(), sep = '.')) |> 
  dplyr::select(where(is.character), silva_Confidence)

# import whale sighting data
scaled_sightings <- read_csv("data/CC_on_effort_scaled_sightings_Ceta_SCB.csv") |>
  rename_with(tolower) |>
  rename_with(~gsub('_.*', '', .x)) |>
  mutate(cruise = str_replace(cruiseid, "CC", "X20")) |>
  dplyr::select(cruise, season, bp, bm, mn)

# impute zeroes with uniform random numbers up to seasonal minima
set.seed(30924)
scaled_sightings_imp <- scaled_sightings |>
  group_by(season) |>
  summarize(across(where(is.numeric), 
                   .fns = list(min = ~min(.x[.x > 0], na.rm = T)),
                   .names = '{.col}.{.fn}')) |>
  right_join(scaled_sightings, by = 'season') |>
  mutate(bp.imp = if_else(bp == 0, 
                          runif(n = length(bp == 0), min = 0, max = bp.min),
                          bp),
         bm.imp = if_else(bm == 0,
                          runif(n = length(bm == 0), min = 0, max = bm.min),
                          bm),
         mn.imp = if_else(mn == 0,
                          runif(n = length(mn == 0), min = 0, max = mn.min),
                          mn))

# subtract seasonal averages
sightings <- scaled_sightings_imp |>
  dplyr::select(cruise, season, ends_with('imp')) |>
  group_by(season) |>
  summarize(across(ends_with('imp'), 
                   .fns = list(mean = ~mean(log(.x), na.rm = T)),
                   .names = 'log.{.col}.{.fn}')) |>
  right_join(scaled_sightings_imp) |>
  mutate(cruise = cruise,
            bp = log(bp.imp) - log.bp.imp.mean,
            bm = log(bm.imp) - log.bm.imp.mean,
            mn = log(mn.imp) - log.mn.imp.mean,
         .keep = "none")

# clr transformation
clr_out <- edna_data |> 
  dplyr::select(-cruise) |> 
  clr() |>
  as_tibble()

# center wrt seasonal (geometric) mean
clr_centered <- edna_data |>
  dplyr::select(cruise) |>
  bind_cols(clr_out) |>
  left_join(dplyr::select(scaled_sightings, cruise, season),
            by = 'cruise') |>
  group_by(season) |>
  mutate(across(starts_with('asv'), ~.x - mean(.x))) |>
  ungroup() |>
  dplyr::select(cruise, starts_with('asv'))

# combine sightings and edna
whales <- inner_join(sightings, clr_centered, by = 'cruise')

## MODEL FITTING ---------------------------------------------------------------

# data inputs
x <- dplyr::select(whales, starts_with('asv'))
mn <- pull(whales, mn) # humpback
bp <- pull(whales, bp) # fin
bm <- pull(whales, bm) # blue
n <- nrow(whales)

# hyperparameter grid
grid_res <- 200
eta_grid_logspace <- rev(1 - exp(seq(-4, -0.01, length = grid_res)))

# loocv to choose hyperparameters
cv_mn <- cv.spls(x, mn, 
                 K = 2, 
                 scale.x = F, scale.y = F,
                 eta = eta_grid_logspace, 
                 fold = n,
                 plot.it = F)

cv_bp <- cv.spls(x, bp, 
                 K = 2, 
                 scale.x = F, scale.y = F,
                 eta = eta_grid_logspace, 
                 fold = n,
                 plot.it = F)

cv_bm <- cv.spls(x, bm, 
                 K = 2, 
                 scale.x = F, scale.y = F,
                 eta = eta_grid_logspace, 
                 fold = n,
                 plot.it = F)

# examine mspe and manually pick eta
eta_df <- tibble(eta = eta_grid_logspace) |>
  bind_cols(bm = cv_bm$mspemat[, 1],
            bp = cv_bp$mspemat[, 1],
            mn = cv_mn$mspemat[, 1]) 

eta_df |>
  pivot_longer(-eta, values_to = 'mspe', names_to = 'species') |>
  ggplot(aes(x = (eta), y = mspe)) +
  geom_path() +
  facet_wrap(~species) +
  scale_y_log10() + 
  geom_vline(aes(xintercept = eta.approx), 
             data = tibble(species = c('bm', 'bp', 'mn'),
                           eta.approx = c(0.48, 0.61, 0.45)))

eta_sel <- tibble(species = c('bm', 'bp', 'mn'),
                  eta.approx = c(0.48, 0.61, 0.45))

# choose closest grid points to selected eta
mspe_df <- eta_df |>
  pivot_longer(-eta, values_to = 'mspe', names_to = 'species') |>
  left_join(eta_sel, by = 'species') |>
  group_by(species) |>
  slice_min(abs(eta - eta.approx)) |>
  dplyr::select(-eta.approx)

eta_bm <- filter(mspe_df, species == 'bm') |> pull(eta)
eta_bp <- filter(mspe_df, species == 'bp') |> pull(eta)
eta_mn <- filter(mspe_df, species == 'mn') |> pull(eta)

# fit spls models
fit_bm <- spls(x, bm, K = 2, eta = eta_bm, 
            scale.x = F, scale.y = F)
fit_bp <- spls(x, bp, K = 2, eta = eta_bp, 
               scale.x = F, scale.y = F)
fit_mn <- spls(x, mn, K = 2, eta = eta_mn, 
               scale.x = F, scale.y = F)

## MODEL OUTPUTS

# extract predictions
pred_bp <- predict(fit_bp)
pred_mn <- predict(fit_mn)
pred_bm <- predict(fit_bm)

# summary statistics for responses
whale_summaries <- whales |>
  dplyr::select(2:4) |>
  pivot_longer(everything(), names_to = 'species', values_to = 'density') |>
  group_by(species) |>
  summarize(mean.log.density = mean(density),
            var.log.density = var(density)) 

fit_summary <- tibble(species = c('bp', 'bm', 'mn'),
                      n.asv = c(nrow(fit_bp$projection),
                                nrow(fit_bm$projection),
                                nrow(fit_mn$projection)),
                      var.expl = c(1 - ((n - 2)*var(bp - pred_bp)/((n - 1)*var(bp))),
                                   1 - ((n - 2)*var(bm - pred_bm)/((n - 1)*var(bm))),
                                   1 - ((n - 2)*var(mn - pred_mn)/((n - 1)*var(mn))))) |>
  left_join(mspe_df, by = 'species') |>
  left_join(whale_summaries, by = 'species') |>
  dplyr::select(-eta)

fit_summary


# # verify fitted model components
# x_sel <- dplyr::select(x, rownames(fit$projection))
# x_sel_proj <- as.matrix(x_sel) %*% fit$projection
# 
# lm(y ~ x_sel_proj)$fitted |> 
#   bind_cols(fitted_vals)
# 
# (fit$projection %*% lm(y ~ x_sel_proj)$coef[2:3]) |>
#   bind_cols(fit$betahat[fit$betahat != 0])

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


# join asvs across species
asv_sel <- full_join(asv_bm, asv_bp, 
          by = c('short.id', 'd', 'p', 'c', 'o', 'f', 'g'),
          suffix = c('.bm', '.bp')) |>
  full_join(asv_mn, 
            by = c('short.id', 'd', 'p', 'c', 'o', 'f', 'g')) |>
  rename(coef.mn = coef) 

# visualize
pal <- colorRampPalette(c('#649dfa', '#02204f'), bias = 1)

asv_sel |>
  mutate(phylum = str_remove(p, 'p__') |> str_trim() |> replace_na('Unknown')) |> 
  dplyr::select(short.id, phylum, starts_with('coef')) |>
  rename_with(~str_remove(.x, 'coef.')) |>
  pivot_longer(bm:mn, names_to = 'species', values_to = 'coef') |>
  mutate(species = factor(species, levels = c('bm', 'bp', 'mn'), labels = c('blue', 'fin','humpback'))) |>
  mutate(phylum = fct_reorder(phylum, coef, ~max(abs(.x), na.rm = T), .na_rm = F),
         coef.trans = 100*(2^coef - 1)) |>
  ggplot(aes(y = phylum)) +
  geom_col(aes(x = coef.trans, fill = phylum),
           width = 1,
           position = position_dodge2(preserve = 'single',
                                      padding = 0.1)) +
  facet_wrap(~species) +
  scale_fill_manual(values = pal(26)) +
  guides(fill = guide_none()) +
  scale_x_continuous(n.breaks = 4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.1, 
                                          color = 'black'),
        panel.grid.major.y = element_line(linewidth = 0.1, 
                                          color = 'grey')) +
  labs(x = 'percent change in median density', 
       y = '',
       title = '')

ggsave(filename = 'rslt/plots/42524/asv-sel-all.png',
       height = 6, width = 8, dpi = 300)

asv_mn |>
  mutate(phylum = str_remove(p, 'p__') |> str_trim() |> replace_na('Unknown')) |> 
  mutate(phylum = fct_reorder(phylum, coef, ~max(abs(.x))),
         x = 100*(2^coef - 1)) |>
  ggplot(aes(y = phylum)) +
  geom_col(aes(x = x, fill = phylum),
           width = 1,
           position = position_dodge2(preserve = 'single',
                                      padding = 0.1)) +
  scale_fill_manual(values = pal(19)) +
  guides(fill = guide_none()) +
  scale_x_continuous(n.breaks = 4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.1, 
                                          color = 'black'),
        panel.grid.major.y = element_line(linewidth = 0.1, 
                                          color = 'grey')) +
  labs(x = 'percent change in density', 
       y = '',
       title = 'humpback')

ggsave(filename = 'rslt/plots/42524/asv-sel-mn.png',
       height = 4, width = 5, dpi = 300)


asv_bm |>
  mutate(phylum = str_remove(p, 'p__') |> str_trim() |> replace_na('Unknown')) |> 
  mutate(phylum = fct_reorder(phylum, coef, ~max(abs(.x))),
         x = 100*(2^coef - 1)) |>
  ggplot(aes(y = phylum)) +
  geom_col(aes(x = x, fill = phylum),
           width = 1,
           position = position_dodge2(preserve = 'single',
                                      padding = 0.1)) +
  scale_fill_manual(values = pal(19)) +
  guides(fill = guide_none()) +
  scale_x_continuous(n.breaks = 4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.1, 
                                          color = 'black'),
        panel.grid.major.y = element_line(linewidth = 0.1, 
                                          color = 'grey')) +
  labs(x = 'percent change in density', 
       y = '',
       title = 'blue')

ggsave(filename = 'rslt/plots/42524/asv-sel-bm.png',
       height = 4, width = 5, dpi = 300)

asv_bp |>
  mutate(phylum = str_remove(p, 'p__') |> str_trim() |> replace_na('Unknown')) |> 
  mutate(phylum = fct_reorder(phylum, coef, ~max(abs(.x))),
         x = 100*(2^coef - 1)) |>
  ggplot(aes(y = phylum)) +
  geom_col(aes(x = x, fill = phylum),
           width = 1,
           position = position_dodge2(preserve = 'single',
                                      padding = 0.1)) +
  scale_fill_manual(values = pal(19)) +
  guides(fill = guide_none()) +
  scale_x_continuous(n.breaks = 4) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.1, 
                                          color = 'black'),
        panel.grid.major.y = element_line(linewidth = 0.1, 
                                          color = 'grey')) +
  labs(x = 'percent change in density', 
       y = '',
       title = 'fin')

ggsave(filename = 'rslt/plots/42524/asv-sel-bp.png',
       height = 4, width = 5, dpi = 300)


# tables
library(openxlsx)

xl_out <- list(bm = mutate(asv_bm, coef = 100*(2^coef - 1)),
     bp = mutate(asv_bp, coef = 100*(2^coef - 1)),
     mn = mutate(asv_mn, coef = 100*(2^coef - 1)),
     fit_summary = fit_summary)

write.xlsx(xl_out, 'rslt/plots/42524/tables.xlsx', rownames = F)

count(group_by(asv_bm, p), name = 'bm') |>
  full_join(count(group_by(asv_bp, p), name = 'bp'), by = 'p') |>
  full_join(count(group_by(asv_mn, p), name = 'mn'), by = 'p')


fit_summary
          