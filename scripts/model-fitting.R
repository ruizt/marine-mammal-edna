library(tidyverse)
library(modelr)
library(spls)

load('data/ncog-18s-processed-2024-05-02.RData')
load('data/ceta-density-processed-2024-05-02.RData')

# combine seasonally adjusted density estimates and seasonally adjusted edna data
whales <- inner_join(log_density_estimates_adj, edna_clr_adj, by = 'cruise')

## MODEL FITTING ---------------------------------------------------------------

# function to fit spls model
fit_fn <- function(.data, .var, .eta){
  df <- as.data.frame(.data)
  x <- dplyr::select(df, starts_with('asv'))
  y <- pull(df, {{.var}})
  fit <- spls(x, y, K = 2, eta = .eta, scale.x = F, scale.y = F)
  return(fit)
}

# function to compute evaluation metrics
metrics <- function(.fit, .data, .var){
  # fit summaries
  y.train <- .fit$y
  y.train.hat <- predict(.fit, type = 'fit')
  train.resid <- y.train - y.train.hat
  n <- length(y.train)
  p <- ncol(.fit$projection)
  df <- nrow(.fit$projection)
  rsq <- as.numeric(1 - (n - p)*var(train.resid)/((n - 1)*var(y.train)))
  
  # squared prediction error
  test <- as.data.frame(.data)
  x.test <- dplyr::select(test, starts_with('asv'))
  y.test <- pull(test, {{.var}})
  y.test.hat <- predict(.fit, newx = x.test, type = 'fit')
  test.resid <- y.test - y.test.hat
  test.spe <- as.numeric(test.resid^2)
  
  # outputs
  suffix <- deparse(substitute(.var))
  out <- data.frame(df = df,
                    rsq = rsq,
                    p = y.test.hat,
                    pe = test.resid,
                    spe = test.spe) |>
    rename_with(~paste0(.x, '.', suffix))
  return(out)
}

# hyperparameter grid
grid_res <- 100
eta_grid <- rev(1 - seq(0.001, 0.999, length = grid_res)^2)

# # leave one out cross validation
# cv_out <- crossv_loo(whales) |>
#   rowwise() |>
#   mutate(eta = list(eta_grid)) |>
#   unnest(eta) |>
#   mutate(fit.bm = map2(train, eta, \(.x, .y) fit_fn(.x, bm, .y)),
#          fit.bp = map2(train, eta, \(.x, .y) fit_fn(.x, bp, .y)),
#          fit.mn = map2(train, eta, \(.x, .y) fit_fn(.x, mn, .y)),
#          bm = map2(fit.bm, test, \(.x, .y) metrics(.x, .y, bm)),
#          bp = map2(fit.bp, test, \(.x, .y) metrics(.x, .y, bp)),
#          mn = map2(fit.mn, test, \(.x, .y) metrics(.x, .y, mn))) |> 
#   dplyr::select(.id, eta, bp, bm, mn) |>
#   unnest(c(bm, bp, mn)) |>
#   pivot_longer(-(1:2)) |>
#   separate_wider_delim(name, delim = '.', names = c('metric', 'species'))

# save(cv_out, file = paste('rslt/comp/cv-eta-', today(), '.RData', sep = ''))

load('rslt/comp/cv-eta-2024-04-28.RData')

# inspect to manually choose eta (0.6 seems reasonable across the board)
cv_out |>
  filter(metric %in% c('df', 'rsq', 'spe')) |>
  group_by(eta, species, metric) |>
  summarize(mean = mean(value),
            sd = sd(value)) |>
  # mutate(sd = if_else(metric == 'spe', NA, sd)) |>
  ggplot(aes(x = eta, y = mean)) +
  facet_wrap(~metric + species, scales = 'free_y', nrow = 3) +
  geom_path() +
  geom_ribbon(aes(ymin = mean - sd,
                  ymax = mean + sd),
              alpha = 0.3) + 
  geom_vline(aes(xintercept = eta.approx), 
             data = tibble(species = c('bm', 'bp', 'mn'),
                           eta.approx = c(0.65, 0.55, 0.63)),
             linetype = 'dotdash')

eta_sel_approx <- tibble(species = c('bm', 'bp', 'mn'),
                  eta.approx = c(0.65, 0.55, 0.63))

# choose closest grid points to selected eta
loopreds_sel <- cv_out |>
  left_join(eta_sel_approx, by = 'species') |>
  group_by(species) |>
  slice_min(abs(eta - eta.approx)) |>
  dplyr::select(-eta.approx) |>
  filter(metric %in% c('p', 'pe', 'spe'))

eta_sel <- loopreds_sel |> distinct(eta, species)

eta_bm <- filter(eta_sel, species == 'bm') |> pull(eta)
eta_bp <- filter(eta_sel, species == 'bp') |> pull(eta)
eta_mn <- filter(eta_sel, species == 'mn') |> pull(eta)

cv_out |>
  left_join(eta_sel_approx, by = 'species') |>
  group_by(species) |>
  slice_min(abs(eta - eta.approx)) |>
  dplyr::select(-eta.approx) |>
  filter(metric %in% c('df', 'rsq', 'spe')) |>
  group_by(species, eta, metric) |>
  summarize(mean = mean(value)) |>
  spread(metric, mean) |>
  xtable::xtable() |>
  print() |>
  clipr::write_clip()

# fit spls models
fit_bm <- fit_fn(whales, bm, eta_bm)
fit_bp <- fit_fn(whales, bp, eta_bp)
fit_mn <- fit_fn(whales, mn, eta_mn)

## MODEL OUTPUTS ---------------------------------------------------------------

# extract responses
x <- dplyr::select(whales, starts_with('asv'))
bm <- pull(whales, bm)
bp <- pull(whales, bp)
mn <- pull(whales, mn)
n <- nrow(whales)

# extract fitted values
fitted_bp <- predict(fit_bp)
fitted_mn <- predict(fit_mn)
fitted_bm <- predict(fit_bm)

# summarize prediction errors from loocv
pred_err_df <- loopreds_sel |>
  filter(metric == 'pe') |>
  group_by(species) |>
  summarize(bias.ratio = mean(value) |> exp(),
            mspe = mean(value^2))

# summary of model fits
fit_summary <- tibble(species = c('bp', 'bm', 'mn'),
                      n.asv = c(nrow(fit_bp$projection),
                                nrow(fit_bm$projection),
                                nrow(fit_mn$projection)),
                      var.expl = c(1 - ((n - 2)*var(bp - fitted_bp)/((n - 1)*var(bp))),
                                   1 - ((n - 2)*var(bm - fitted_bm)/((n - 1)*var(bm))),
                                   1 - ((n - 2)*var(mn - fitted_mn)/((n - 1)*var(mn))))) |>
  left_join(pred_err_df)

fit_summary |>
  xtable::xtable() |>
  print() |>
  clipr::write_clip()

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

count(group_by(asv_bm, d), name = 'bm') |>
  full_join(count(group_by(asv_bp, d), name = 'bp'), by = 'd') |>
  full_join(count(group_by(asv_mn, d), name = 'mn'), by = 'd') |>
  pivot_longer(-d) |>
  mutate(phylum = str_remove(d, 'd__') |> str_trim()) |>
  drop_na() |>
  ggplot(aes(x = name, fill = phylum)) +
  geom_bar(position = 'fill') +
  scale_y_continuous(labels = scales::percent)

asv_sel |>
  pivot_longer(starts_with('coef')) |>
  mutate(species = str_remove(name, 'coef.')) |>
  drop_na(value) |>
  group_by(species) |>
  summarize(across(p:g, ~length(unique(.x)))) |>
  write.csv('rslt/tbl/classification-counts.csv')

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
  rowwise() |>
  mutate(n.na = sum(is.na(c(bm, bp, mn)))) |>
  filter(n.na < 2, 
         phylum != 'Unknown', 
         phylum != 'uncultured') |>
  pivot_longer(bm:mn, names_to = 'species', values_to = 'coef') |>
  mutate(species = factor(species, levels = c('bm', 'bp', 'mn'), labels = c('blue', 'fin','humpback'))) |>
  group_by(phylum) |>
  mutate(count = n()) |>
  ungroup() |>
  mutate(phylum = fct_reorder(phylum, count),
         coef.trans = 100*(2^coef - 1),
         p.ix = as.numeric(phylum) %% 2) |>
  ggplot(aes(y = phylum)) +
  geom_col(aes(x = coef.trans, fill = factor(p.ix)),
           position = position_dodge2(preserve = 'single')) +
  facet_wrap(~species) +
  scale_fill_manual(values = c("#04004c", "#8e87ea")) +
  guides(fill = guide_none()) +
  scale_x_continuous(n.breaks = 4) +
  theme_bw(base_size = 16) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "#04004c") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(linewidth = 0.1, 
                                          color = "#04004c"),
        panel.background = element_rect(fill = "#fffffeff"),
        plot.background = element_rect(fill = "#fffffeff"),
        text = element_text(color = "#49485a")) +
  labs(x = 'percent change in median scaled sighting', 
       y = '')

ggsave(filename = 'rslt/plots/51024/asv-sel-all.png',
       height = 4, width = 6, dpi = 400)


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

p1 <- loopreds_sel |>
  filter(metric == 'p') |>
  bind_cols(obs = c(bm, bp, mn)) |>
  ggplot(aes(x = obs, y = value)) +
  geom_point() +
  facet_wrap(~species) +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw() +
  labs(x = '',
       y = 'prediction')

tibble(fitted = c(fitted_bm, fitted_bp, fitted_mn),
       observed = c(bm, bp, mn),
       species = rep(c('blue', 'fin', 'humpback'), 
                     times = c(length(bm), length(bp), length(mn)))) |>
  ggplot(aes(x = observed, y = fitted)) +
  facet_wrap(~species) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  theme_bw(base_size = 18) +
  labs(x = 'scaled sighting logratio',
       y = 'fitted value') +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.1, 
                                        color = 'black'),
        panel.background = element_rect(fill = "#fffffeff"),
        plot.background = element_rect(fill = "#fffffeff"),
        text = element_text(color = "#49485a"))

ggsave(filename = 'rslt/plots/51024/obs-fitted.png',
       height = 2.5, width = 6, dpi = 400)


library(patchwork)

fig <- p1 + p2 + plot_layout(nrow = 2)

ggsave(filename = 'rslt/plots/42524/pred-fitted.png',
       height = 4, width = 5, dpi = 300)
