library(tidyverse)
library(spls)
library(openxlsx)

load('data/ncog-18s-processed-2024-04-27.RData')
load('data/ceta-density-processed-2024-04-27.RData')

# combine seasonally adjusted density estimates and seasonally adjusted edna data
whales <- inner_join(log_density_estimates_adj, edna_clr_adj, by = 'cruise')

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
                           eta.approx = c(0.65, 0.58, 0.635)))

eta_sel <- tibble(species = c('bm', 'bp', 'mn'),
                  eta.approx = c(0.65, 0.58, 0.635))

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

## MODEL OUTPUTS ---------------------------------------------------------------

# extract predictions
pred_bp <- predict(fit_bp)
pred_mn <- predict(fit_mn)
pred_bm <- predict(fit_bm)

# summary statistics for responses
whale_summaries <- whales |>
  dplyr::select(2:4) |>
  pivot_longer(everything(), names_to = 'species', values_to = 'density') |>
  group_by(species) |>
  summarize(var.log.adj.density = var(density)) 

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

# export selected asvs and fit summaries as excel workbook
xl_out <- list(blue = mutate(asv_bm, coef = 100*(2^coef - 1)),
               fin = mutate(asv_bp, coef = 100*(2^coef - 1)),
               humpback = mutate(asv_mn, coef = 100*(2^coef - 1)),
               fit.summary = fit_summary)

openxlsx::write.xlsx(xl_out, 
           paste('rslt/tbl/fit-tables-', today(), '.xlsx', sep = ''), 
           rownames = F)



## DRAFT MATERIAL --------------------------------------------------------------

count(group_by(asv_bm, p), name = 'bm') |>
  full_join(count(group_by(asv_bp, p), name = 'bp'), by = 'p') |>
  full_join(count(group_by(asv_mn, p), name = 'mn'), by = 'p')

# # verify fitted model components
# x_sel <- dplyr::select(x, rownames(fit$projection))
# x_sel_proj <- as.matrix(x_sel) %*% fit$projection
# 
# lm(y ~ x_sel_proj)$fitted |> 
#   bind_cols(fitted_vals)
# 
# (fit$projection %*% lm(y ~ x_sel_proj)$coef[2:3]) |>
#   bind_cols(fit$betahat[fit$betahat != 0])

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
