library(tidyverse)
library(lubridate)
library(patchwork)

## ASV TABLES ------------------------------------------------------------------
out_dir <- 'rslt/tbl/'
fs::dir_create(out_dir)

# candidate asvs from 18sv9
load('data/processed/ncog18sv9.RData')
asv_taxa_18sv9 <- asv_taxa 

# candidate asvs from 18sv4
load('data/processed/ncog18sv4.RData')
asv_taxa_18sv4 <- asv_taxa

# candidate asvs from 16s
load('data/processed/ncog16s.RData')
asv_taxa_16s <- asv_taxa

# selected asvs from 18sv9 & scaled sightings model 
load('rslt/models/scaled-sightings/fitted-models-18sv9-ss.RData')
load('data/processed/ncog18sv9.RData')
sel_asv_18sv9_ss <- fitted_models |>
  mutate(model.species = factor(species, 
                                levels = c('bm', 'bp', 'mn'),
                                labels = c('Blue whale', 'Fin whale', 'Humpback whale')),
         asv.id = map(data, ~select(.x, -y) |> colnames())) |>
  select(model.species, asv.id) |>
  unnest(everything()) |>
  left_join(asv_taxa, join_by(asv.id == short.id)) 

# selected asvs from 18sv4 & scaled sightings model 
load('rslt/models/scaled-sightings/fitted-models-18sv4-ss.RData')
load('data/processed/ncog18sv4.RData')
sel_asv_18sv4_ss <- fitted_models |>
  mutate(model.species = factor(species, 
                                levels = c('bm', 'bp', 'mn'),
                                labels = c('Blue whale', 'Fin whale', 'Humpback whale')),
         asv.id = map(data, ~select(.x, -y) |> colnames())) |>
  select(model.species, asv.id) |>
  unnest(everything()) |>
  left_join(asv_taxa, join_by(asv.id == short.id)) 

# selected asvs from 16s & scaled sightings model 
load('rslt/models/scaled-sightings/fitted-models-16s-ss.RData')
load('data/processed/ncog16s.RData')
sel_asv_16s_ss <- fitted_models |>
  mutate(model.species = factor(species, 
                                levels = c('bm', 'bp', 'mn'),
                                labels = c('Blue whale', 'Fin whale', 'Humpback whale')),
         asv.id = map(data, ~select(.x, -y) |> colnames())) |>
  select(model.species, asv.id) |>
  unnest(everything()) |>
  left_join(asv_taxa, join_by(asv.id == short.id)) 

# write as excel sheets
sheets <- list("18Sv9-candidates" = asv_taxa_18sv9, 
               "18Sv4-candidates" = asv_taxa_18sv4,
               "16S-candidates" = asv_taxa_16s,
               "18Sv9-selected" = sel_asv_18sv9_ss,
               "18Sv4-selected" = sel_asv_18sv4_ss,
               "16S-selected" = sel_asv_16s_ss)
writexl::write_xlsx(sheets, paste(out_dir, 'asv-tables.xlsx', sep = ''))

## TABLE: MODEL SUMMARIES ------------------------------------------------------

# model summary for 18sv9
load('rslt/models/scaled-sightings/fitted-models-18sv9-ss.RData')
summary_tbl_18sv9 <- pred_metrics |> 
  filter(scale == 'ss', metric == 'rmspe') |>
  select(-metric, -scale) |>
  pivot_longer(c(bm, bp, mn), names_to = 'species', values_to = 'rmspe') |>
  pivot_wider(names_from = 'model', values_from = 'rmspe') |>
  mutate(pct.reduction = (naive - pls)/naive,
         marker = '18SV9') |>
  rename(rmspe = pls) |>
  left_join(fit_metrics, join_by(species)) |>
  select(species, marker, n.asv, adj.rsq.ss, rmspe, pct.reduction)

# model summary for 18sv4
load('rslt/models/scaled-sightings/fitted-models-18sv4-ss.RData')
summary_tbl_18sv4 <- pred_metrics |> 
  filter(scale == 'ss', metric == 'rmspe') |>
  select(-metric, -scale) |>
  pivot_longer(c(bm, bp, mn), names_to = 'species', values_to = 'rmspe') |>
  pivot_wider(names_from = 'model', values_from = 'rmspe') |>
  mutate(pct.reduction = (naive - pls)/naive,
         marker = '18SV4') |>
  rename(rmspe = pls) |>
  left_join(fit_metrics, join_by(species)) |>
  select(species, marker, n.asv, adj.rsq.ss, rmspe, pct.reduction)

# model summary for 16s
load('rslt/models/scaled-sightings/increased-res-fitted-models-16s-ss.RData')
summary_tbl_16s <- pred_metrics |> 
  filter(scale == 'ss', metric == 'rmspe') |>
  select(-metric, -scale) |>
  pivot_longer(c(bm, bp, mn), names_to = 'species', values_to = 'rmspe') |>
  pivot_wider(names_from = 'model', values_from = 'rmspe') |>
  mutate(pct.reduction = (naive - pls)/naive,
         marker = '16S') |>
  rename(rmspe = pls) |>
  left_join(fit_metrics, join_by(species)) |>
  select(species, marker, n.asv, adj.rsq.ss, rmspe, pct.reduction)

# join
model_summaries <- bind_rows(summary_tbl_16s, summary_tbl_18sv4, summary_tbl_18sv9) |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales'))) |>
  arrange(species, marker)
model_summaries |>
  mutate(across(4:6, ~round_z(.x, 3))) |>
  write_csv('rslt/_draft/ss-model-summaries.csv')


# # counts of phyla of selected asvs
# sel_asv_18sv9_ss |>
#   group_by(model.species, d, p) |>
#   count() |>
#   arrange(model.species, desc(n)) |>
#   group_by(model.species) |>
#   mutate(prop = n/sum(n)) |>
#   print(n = 200)

# fit <- fitted_models |>
#   slice(1) |> pull(fit) %$% .[[1]]
# 
# coef(fit)[, , 1]
# 
# coefs <- fitted_models |>
#   mutate(coef = map(fit, ~coef(.x)[, , 1]),
#          asv = map(coef, names)) |>
#   unnest(c(coef, asv)) |>
#   select(species, coef, asv) |>
#   arrange(species, desc(coef))
# 
# coef_df <- loo_preds |>
#   left_join(best_settings, join_by(species), suffix = c('', '.best')) |>
#   filter(setting == setting.best) |>
#   select(test.cruise, test.season, species, fit) |>
#   mutate(coef = map(fit, ~coef(.x)[, , 1]),
#          asv = map(coef, names)) |>
#   unnest(c(coef, asv)) |>
#   group_by(species, asv) |>
#   summarize(coef.se = sd(coef)) |>
#   full_join(coefs, join_by(species, asv))
# 
# coef_df |>
#   left_join(asv_taxa, join_by(asv == short.id)) |>
#   # filter(species == '') |>
#   ggplot(aes(x = coef, 
#              y = fct_reorder(factor(asv), desc(coef)),
#              color = d)) +
#   geom_point() +
#   geom_errorbarh(aes(xmin = coef - 2*coef.se, xmax = coef + 2*coef.se)) +
#   geom_vline(xintercept = 0) +
#   labs(x = NULL, y = NULL) +
#   theme_minimal() +
#   theme(axis.text.y = element_blank()) +
#   facet_wrap(~species)

## FIGURE: PREDICTIONS ---------------------------------------------------------

# predictions from 18sv9
load('rslt/models/scaled-sightings/fitted-models-18sv9-ss.RData')
pred_pts_18sv9 <- loo_pred_df |>
  mutate(cruise = test.cruise,
         cruise.ym = ym(test.cruise),
         year = year(cruise.ym),
         quarter = factor(test.season, 
                          levels = c('winter', 'spring', 'summer', 'fall')) |>
           as.numeric(),
         cruise.yq = paste(year, quarter, sep = '-') |> yq(),
         species = factor(test.species, 
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales'))) |>
  select(cruise, cruise.ym, cruise.yq, species, ss.obs, ss.pred, ss.pred.q1, ss.pred.q3) |>
  arrange(species, cruise.yq) |>
  rename(observed = ss.obs,
         predicted = ss.pred) |>
  mutate(marker = "18S-V9")

# pred_pts_18sv9 |>
#   ggplot(aes(x = cruise.ym,
#              y = observed)) +
#   facet_wrap(~species, nrow = 3) +
#   geom_path() +
#   geom_path(aes(y = predicted), linetype = 2, color = 'blue') +
#   geom_ribbon(aes(ymin = ss.pred.q1,
#                   ymax = ss.pred.q3, 
#                   x = cruise.ym),
#               inherit.aes = F,
#               alpha = 0.1,
#               fill = 'blue')

# predictions from 18sv4
load('rslt/models/scaled-sightings/fitted-models-18sv4-ss.RData')
pred_pts_18sv4 <- loo_pred_df |>
  mutate(cruise = test.cruise,
         cruise.ym = ym(test.cruise),
         year = year(cruise.ym),
         quarter = factor(test.season, 
                          levels = c('winter', 'spring', 'summer', 'fall')) |>
           as.numeric(),
         cruise.yq = paste(year, quarter, sep = '-') |> yq(),
         species = factor(test.species, 
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales'))) |>
  select(cruise, cruise.ym, cruise.yq, species, ss.obs, ss.pred) |>
  arrange(species, cruise.yq) |>
  rename(observed = ss.obs,
         predicted = ss.pred) |>
  mutate(marker = "18S-V4")

# predictions from 16s
load('rslt/models/scaled-sightings/increased-res-fitted-models-16s-ss.RData')
pred_pts_16s <- loo_pred_df |>
  mutate(cruise = test.cruise,
         cruise.ym = ym(test.cruise),
         year = year(cruise.ym),
         quarter = factor(test.season, 
                          levels = c('winter', 'spring', 'summer', 'fall')) |>
           as.numeric(),
         cruise.yq = paste(year, quarter, sep = '-') |> yq(),
         species = factor(test.species, 
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales'))) |>
  select(cruise, cruise.ym, cruise.yq, species, ss.obs, ss.pred) |>
  arrange(species, cruise.yq) |>
  rename(observed = ss.obs,
         predicted = ss.pred) |>
  mutate(marker = "16S")



# combine predictions from each marker
pred_pts <- bind_rows(pred_pts_18sv9) |>
  mutate(key = paste('Predicted (', marker, ')', sep = '')) |>
  select(-observed, -marker) |>
  rename(value = predicted)

# observations
obs_pts <- pred_pts_18sv9 |>
  select(-predicted, -marker, -ss.pred.q1, -ss.pred.q3) |>
  rename(value = observed) |>
  mutate(key = 'Observed')

# color palette
pal <- c('#000000', RColorBrewer::brewer.pal(name = 'Dark2', n = 3))

# predictions from each marker overlaid on serial observations
p1 <- bind_rows(pred_pts, obs_pts) |>
  mutate(key = str_remove_all(key, '-')) |>
  ggplot(aes(x = cruise.ym, 
             y = value,
             color = key, 
             linetype = key)) +
  facet_wrap(~species, nrow = 3, scale = 'fixed') +
  geom_path(linewidth = 0.4) +
  geom_ribbon(inherit.aes = F,
              aes(x = cruise.ym,
                  ymin = ss.pred.q1,
                  ymax = ss.pred.q3,
                  fill = key),
              alpha = 0.1) +
  scale_x_date(date_breaks = '1 years', labels = ~year(.x)) +
  scale_color_manual(values = pal) +
  # scale_y_sqrt() +
  theme_bw() +
  # geom_blank(data = pivot_longer(fit_pts_df, c(Observed, Fitted)),
  #            aes(linetype = NULL)) +
  guides(linetype = guide_legend(title = NULL, position = 'top', nrow = 2),
         color = guide_legend(title = NULL, position = 'top', nrow = 2),
         fill = guide_legend(title = NULL, position = 'top', nrow = 2)) +
  labs(x = 'Year', y = 'Sightings per 1000km', title = NULL) +
  theme(panel.grid = element_line(linewidth = 0.1, color = 'black'),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 11))

p1
# predictions vs observations
p2 <- bind_rows(pred_pts_16s, pred_pts_18sv4, pred_pts_18sv9) |>
  mutate(marker = str_remove_all(marker, '-'),
         species = str_remove_all(species, ' whales')) |>
  ggplot(aes(x = observed, y = predicted)) +
  facet_grid(marker~species) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4) +
  theme_bw() +
  theme(panel.grid = element_line(linewidth = 0.1, color = 'black'),
        panel.grid.minor = element_blank(),
        text = element_text(size = 11)) +
  labs(x = 'Observed value', y = 'Predicted value')

# correlations between predictions and observations (for labeling)
pred_cors <- bind_rows(pred_pts_16s, pred_pts_18sv4, pred_pts_18sv9) |>
  mutate(marker = str_remove_all(marker, '-'),
         species = str_remove_all(species, ' whales')) |>
  group_by(species, marker) |>
  summarize(r = cor(observed, predicted), 
            observed = 30, 
            predicted = 0.05,
            .groups = 'drop') 

# function to retain trailing zeroes
round_z <- function(x, n){sprintf("%.3f", round(x,n))}
# round_z(0.910, 3)

# add correlation annotations to pred v. obs plots
p2_ann <- p2 + geom_label(data = pred_cors, 
                          aes(label = paste('r =', round_z(r, 3))),
                          size = 8,
                          size.unit = 'pt',
                          hjust = 1)

# panel layouts
# fig_long <- p1 + p2_ann +
#   plot_layout(nrow = 2, ncol = 1, heights = c(1.25, 1)) +
#   plot_annotation(tag_levels = 'A')
# ggsave(filename = 'rslt/_draft/plots/preds-ss-long.png',
#        width = 3.5, height = 6, units = 'in', dpi = 300)

fig_wide <- p1 + p2_ann + 
  plot_layout(nrow = 1, ncol = 2, widths = c(1, 1)) 
ggsave(fig_wide, filename = 'rslt/_draft/plots/preds-ss-wide.png',
       width = 6.5, height = 4, units = 'in', dpi = 400)


## FIGURE: MODEL DIAGNOSTICS ---------------------------------------------------

# predictions from 18sv9
load('rslt/models/scaled-sightings/fitted-models-18sv9-ss.RData')
fit_pts_18sv9 <- fit_df |>
  mutate(cruise.ym = ym(cruise),
         year = year(cruise.ym),
         quarter = factor(season,
                          levels = c('winter', 'spring', 'summer', 'fall')) |>
           as.numeric(),
         cruise.yq = paste(year, quarter, sep = '-') |> yq(),
         species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales'))) |>
  select(cruise, cruise.ym, cruise.yq, species, y, fitted, lr.resid) |>
  arrange(species, cruise.yq) |>
  rename(lr.fit = fitted, lr.obs = y) |>
  mutate(marker = '18SV9')

# predictions from 18sv4
load('rslt/models/scaled-sightings/fitted-models-18sv4-ss.RData')
fit_pts_18sv4 <- fit_df |>
  mutate(cruise.ym = ym(cruise),
         year = year(cruise.ym),
         quarter = factor(season,
                          levels = c('winter', 'spring', 'summer', 'fall')) |>
           as.numeric(),
         cruise.yq = paste(year, quarter, sep = '-') |> yq(),
         species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales'))) |>
  select(cruise, cruise.ym, cruise.yq, species, y, fitted, lr.resid) |>
  arrange(species, cruise.yq) |>
  rename(lr.fit = fitted, lr.obs = y) |>
  mutate(marker = '18SV4')

# predictions from 16s
load('rslt/models/scaled-sightings/increased-res-fitted-models-16s-ss.RData')
fit_pts_16s <- fit_df |>
  mutate(cruise.ym = ym(cruise),
         year = year(cruise.ym),
         quarter = factor(season,
                          levels = c('winter', 'spring', 'summer', 'fall')) |>
           as.numeric(),
         cruise.yq = paste(year, quarter, sep = '-') |> yq(),
         species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales'))) |>
  select(cruise, cruise.ym, cruise.yq, species, y, fitted, lr.resid) |>
  arrange(species, cruise.yq) |>
  rename(lr.fit = fitted, lr.obs = y) |>
  mutate(marker = '16S')

# merge above outputs
fit_pts <- bind_rows(fit_pts_16s, fit_pts_18sv4, fit_pts_18sv9) 

# for rescaling axes without <scales = 'free_x'>
zscore <- function(x){(x - mean(x))/sd(x)}

# annotations for adding adjusted rsq
ann_data <- fit_pts |>
  group_by(species, marker) |>
  mutate(across(c(lr.fit, lr.obs), zscore)) |>
  summarize(lr.obs = max(lr.obs),
            lr.fit = min(lr.fit)) |>
  left_join(model_summaries, join_by(species, marker)) |>
  ungroup() |>
  mutate(species = fct_relabel(species, ~str_remove_all(.x, ' whales')),
         lr.fit = min(lr.fit) + 0.5, 
         lr.obs = max(lr.obs))

# observed vs fitted
obs_fit <- fit_pts |>
  mutate(species = fct_relabel(species, ~str_remove_all(.x, ' whales'))) |>
  group_by(species, marker) |>
  mutate(across(c(lr.fit, lr.obs), zscore)) |>
  ggplot(aes(x = lr.obs, y = lr.fit)) +
  facet_grid(species ~ marker) +
  # facet_wrap(species ~ marker, scales = 'free') +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4) +
  # geom_smooth(span = 1.5) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = 'black', linewidth = 0.1),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = 'Observed values', y = 'Fitted values', title = 'A. Model fit') +
  geom_label(data = ann_data, 
             aes(label = paste("R^'2' == ", round(adj.rsq.ss, 2))),
             parse = T,
             size = 8,
             size.unit = 'pt',
             hjust = 1)

obs_fit

# plot residuals vs fit for each model
resid_fit <- bind_rows(fit_pts_16s, fit_pts_18sv4, fit_pts_18sv9) |>
  mutate(species = fct_relabel(species, ~str_remove_all(.x, ' whales'))) |>
  group_by(species, marker) |>
  mutate(across(c(lr.fit, lr.resid), zscore)) |>
  ggplot(aes(x = lr.fit, y = lr.resid)) +
  facet_grid(species ~ marker) +
  # facet_wrap(species ~ marker, scales = 'free') +
  geom_point() +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  # geom_smooth(span = 1.5) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = 'black', linewidth = 0.1),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = 'Fitted values', y = 'Residuals', title = 'B. Residual diagnostics')

resid_fit

# function for residual autocorrelation
pacf_fn <- function(x){
  pacf_out <- pacf(x, plot = F, lag.max = 10)
  out <- bind_cols(pacf = c(1, pacf_out$acf[, 1, 1]),
                   lag = c(0, pacf_out$lag[, , 1]),
                   se = 2/sqrt(pacf_out$n.used))
  return(out)
}

# plot pacf for each model
resid_pacf <- bind_rows(fit_pts_16s, fit_pts_18sv4, fit_pts_18sv9) |>
  select(species, marker, lr.resid) |>
  mutate(species = fct_relabel(species, ~str_remove_all(.x, ' whales'))) |>
  nest(resids = lr.resid, .by = c(species, marker)) |>
  mutate(pacf = map(resids, pacf_fn)) |>
  unnest(pacf) |>
  ggplot(aes(x = lag)) +
  facet_grid(species ~ marker) +
  scale_y_continuous(limits = c(-1, 1), n.breaks = 4) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  geom_linerange(aes(ymin = 0, ymax = pacf),
                 linewidth = 0.4) +
  geom_ribbon(aes(ymin = -se, ymax = se), 
              fill = 'blue', 
              alpha = 0.2) +
  theme_bw() +
  theme(panel.grid = element_line(color = 'black', linewidth = 0.1),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = 'Lag', title = 'C. Residual PACF')

resid_pacf

# panel layout
fig_resid <- obs_fit + resid_fit + resid_pacf + plot_layout(nrow = 1, widths = c(1, 1, 1))
ggsave(fig_resid, filename = 'rslt/_draft/plots/resids-ss-diagnostics.png',
       width = 9, height = 4, units = 'in', dpi = 400)


## TABLE: SELECTION CONSISTENCY ------------------------------------------------

# read in results for 18sv9-ss
outer_val_ss <- dir_ls('rslt/loocv/18sv9-ss/_outer-validation-ss/') |>
  lapply(read_rds) %>%
  Reduce(bind_rows, .) |>
  unnest(sel.asv) |>
  group_by(outer.id, species, eta, ncomp, sel.asv) |>
  count() |>
  ungroup() |>
  filter(n > pi.max*(n.obs - 1)) |>
  select(outer.id, species, sel.asv) |>
  group_by(outer.id, species) |>
  distinct(sel.asv) |>
  nest(ss = sel.asv) |>
  mutate(ss = map(ss, ~pull(.x, sel.asv))) |>
  ungroup()

# function to compute "soft intersection" of sets in <set_list>
intersect_fn <- function(set_list, thresh){
  out <- tibble(asv = Reduce(c, set_list)) |> 
    group_by(asv) |> 
    count() |>
    filter(n >= thresh*length(set_list)) |>
    pull(asv)
  return(out)
}

# function to compute no. of elements in union of sets in <set_list>
union_fn <- function(set_list){
  out <- Reduce(c, set_list) |> unique()
  return(out)
}

# asv-specific stability
ss_asv_stbl <- outer_val_ss |>
  group_by(species) |>
  summarize(int = intersect_fn(ss, 0.8) |> list(),
            un = union_fn(ss) |> list()) |>
  mutate(int = map(int, length),
         un = map(un, length),
         prop.stable = map2(int, un, ~.x/.y)) |>
  unnest(everything()) |>
  rename(intersection = int,
         union = un,
         asv.stability = prop.stable) |>
  mutate(method = 'ss', marker = '18SV9') |>
  select(method, species, marker, intersection, union, asv.stability)

# # family-specific stability
# outer_val_ss |>
#   unnest(ss) |>
#   left_join(asv_taxa, join_by(ss == short.id)) |>
#   unite(classification, c(d, p, c, o, g)) |>
#   select(species, ss, classification, outer.id) |>
#   group_by(species, outer.id, classification) |>
#   count() |>
#   pivot_wider(names_from = 'outer.id', values_from = 'n') |>
#   pivot_longer(-c(species, classification)) |>
#   mutate(selected = !is.na(value)) |>
#   group_by(species, classification) |>
#   summarize(sel.freq = mean(selected)) |>
#   summarize(prop.stable = mean(sel.freq > 0.8))

outer_val_spls <- dir_ls('rslt/loocv/18sv9-ss/_outer-validation-spls/') |>
  lapply(read_rds) %>%
  Reduce(bind_rows, .) |>
  select(.id, species, sel.asv)

# asv-specific stability for spls models
spls_asv_stbl <- best_spls_val |> 
  group_by(species) |>
  summarize(int = intersect_fn(sel.asv, 0.8) |> list(),
            un = union_fn(sel.asv) |> list()) |>
  mutate(int = map(int, length),
         un = map(un, length),
         prop.stable = map2(int, un, ~.x/.y)) |>
  unnest(everything()) |>
  rename(intersection = int,
         union = un,
         asv.stability = prop.stable) |>
  mutate(method = 'spls', marker = '18SV9') |>
  select(method, species, marker, intersection, union, asv.stability)

# # family-specific stability for spls models
# best_spls_val |>
#   unnest(sel.asv) |>
#   select(.id, species, sel.asv) |>
#   left_join(asv_taxa, join_by(sel.asv == short.id)) |>
#   unite(classification, c(d, p, c, o, f, g)) |>
#   group_by(.id, species, classification) |>
#   count() |>
#   pivot_wider(names_from = '.id', values_from = 'n') |>
#   pivot_longer(-c(species, classification), names_to = '.id') |>
#   mutate(selected = !is.na(value)) |>
#   group_by(species, classification) |>
#   summarize(sel.freq = mean(selected)) |>
#   summarize(prop.stable = mean(sel.freq > 0.8))

outerval_tbl_18sv9ss <- bind_rows(spls_asv_stbl, ss_asv_stbl) |>
  arrange(species, method)

outerval_tbl_18sv9ss
