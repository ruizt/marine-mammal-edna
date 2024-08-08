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


# remove data from environment
rm(list = setdiff(ls(), 'out_dir'))

# # counts of phyla of selected asvs
# sel_asv_18sv9_ss |>
#   group_by(model.species, d, p) |>
#   count() |>
#   arrange(model.species, desc(n)) |>
#   group_by(model.species) |>
#   mutate(prop = n/sum(n)) |>
#   print(n = 200)
# 
# fit <- fitted_models |>
#   slice(1) |> pull(fit) %$% .[[1]]
# 
# coef(fit)[, , 1]
# 
# fitted_models |>
#   mutate(coef = map(fit, ~coef(.x)[, , 1]),
#          asv = map(coef, names)) |>
#   unnest(c(coef, asv)) |>
#   select(species, coef, asv) |>
#   group_by(species) |>
#   arrange(desc(coef))

## FIGURES ---------------------------------------------------------------------

# # predictions (smooth)
# smooth_fn <- function(.data, .xout, .var){
#   .yout <- spline(x = .data$cruise.ym, 
#                   y = pull(.data, {{.var}}), 
#                   xout = .xout, 
#                   method = 'natural')$y
#   .yout[.yout < 0] <- 0
#   return(.yout)
# }

pred_pts_df <- loo_pred_df |>
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
  rename(Observed = ss.obs,
         Predicted = ss.pred)

fit_pts_df <- fit_df |>
  mutate(cruise.ym = ym(cruise),
         year = year(cruise.ym),
         quarter = factor(season, 
                          levels = c('winter', 'spring', 'summer', 'fall')) |>
           as.numeric(),
         cruise.yq = paste(year, quarter, sep = '-') |> yq(),
         species = factor(species, 
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales'))) |>
  select(cruise, cruise.ym, cruise.yq, species, ss.obs, ss.fit) |>
  arrange(species, cruise.yq) |>
  rename(Observed = ss.obs, Fitted = ss.fit)

# pred_interp_df <- pred_pts_df |>
#   nest(data = -species) |>
#   mutate(date = map(data, ~seq_range(.x$cruise.yq, n = 200)), 
#          Prediction = map2(data, date, ~smooth_fn(.x, .y, Prediction)),
#          Observation = map2(data, date, ~smooth_fn(.x, .y, Observation))) |>
#   select(-data) |>
#   unnest(everything())
#
# pred_interp_df |>
#   pivot_longer(c(Observation, Prediction)) |>
#   ggplot(aes(x = date, y = value, linetype = name)) +
#   facet_wrap(~species, nrow = 3, scale = 'free_y') +
#   geom_path(linewidth = 0.3) +
#   geom_point(data = pivot_longer(pred_pts_df, c(Observation, Prediction)),
#              aes(x = cruise.ym, shape = name)) +
#   theme_bw() +
#   guides(linetype = guide_legend(title = NULL, position = 'top'),
#          shape = guide_legend(title = NULL, position = 'top')) +
#   labs(x = NULL, y = 'Sightings per 1000km', title = NULL) +
#   theme(panel.grid = element_line(linewidth = 0.1, color = 'black'),
#         panel.grid.minor.y = element_blank())

p1 <- pred_pts_df |>
  pivot_longer(c(Observed, Predicted)) |>
  ggplot(aes(x = cruise.ym, y = value, linetype = name)) +
  facet_wrap(~species, nrow = 3, scale = 'free_y') +
  geom_path(linewidth = 0.3) +
  scale_x_date(date_breaks = '1 years', labels = ~year(.x)) +
  theme_bw() +
  geom_blank(data = pivot_longer(fit_pts_df, c(Observed, Fitted)),
             aes(linetype = NULL)) +
  guides(linetype = guide_legend(title = NULL, position = 'top'),
         shape = guide_legend(title = NULL, position = 'top')) +
  labs(x = NULL, y = 'Sightings per 1000km', title = NULL) +
  theme(panel.grid = element_line(linewidth = 0.1, color = 'black'),
        panel.grid.minor = element_blank())


p1

# fit_interp_df <- fit_pts_df |>
#   nest(data = -species) |>
#   mutate(date = map(data, ~seq_range(.x$cruise.yq, n = 200)), 
#          Fitted = map2(data, date, ~smooth_fn(.x, .y, Fitted)),
#          Observation = map2(data, date, ~smooth_fn(.x, .y, Observation))) |>
#   select(-data) |>
#   unnest(everything())
# 
# fit_interp_df |>
#   pivot_longer(c(Observation, Fitted)) |>
#   ggplot(aes(x = date, y = value, linetype = fct_rev(name))) +
#   facet_wrap(~species, nrow = 3, scale = 'free_y') +
#   geom_path(linewidth = 0.3) +
#   geom_point(data = pivot_longer(fit_pts_df, c(Observation, Fitted)),
#              aes(x = cruise.yq, shape = fct_rev(name))) +
#   theme_bw() +
#   guides(linetype = guide_legend(title = NULL, position = 'top'),
#          shape = guide_legend(title = NULL, position = 'top')) +
#   labs(x = NULL, y = 'Sightings per 1000km', title = NULL)

p2 <- fit_pts_df |>
  pivot_longer(c(Observed, Fitted)) |>
  ggplot(aes(x = cruise.ym, y = value, linetype = fct_rev(name))) +
  facet_wrap(~species, nrow = 3, scale = 'free_y') +
  geom_path(linewidth = 0.3) +
  theme_bw() +
  geom_blank(data = pivot_longer(pred_pts_df, c(Observed, Predicted)),
             aes(linetype = NULL)) +
  scale_x_date(date_breaks = '1 years', labels = ~year(.x)) +
  guides(linetype = guide_legend(title = NULL, position = 'top'),
         shape = guide_legend(title = NULL, position = 'top')) +
  labs(x = NULL, y = 'Sightings per 1000km', title = NULL) +
  theme(panel.grid = element_line(linewidth = 0.1, color = 'black'),
        panel.grid.minor = element_blank())


p3 <- fit_pts_df |>
  ggplot(aes(x = Observed, y = Fitted)) +
  facet_wrap(~species) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3) +
  theme_bw() +
  theme(panel.grid = element_line(linewidth = 0.1, color = 'black'),
        panel.grid.minor = element_blank()) +
  labs(x = 'Observed value', y = 'Fitted value')

p4 <- pred_pts_df |>
  ggplot(aes(x = Observed, y = Predicted)) +
  facet_wrap(~species) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3) +
  theme_bw() +
  theme(panel.grid = element_line(linewidth = 0.1, color = 'black'),
        panel.grid.minor = element_blank()) +
  labs(x = 'Observed value', y = 'Predicted value')


fig <- p1 + p2 + p4 + p3 + 
  plot_layout(ncol = 2, nrow = 2, heights = c(2, 1)) +
  plot_annotation(tag_levels = 'A')

fig

p1

ggsave(filename = 'rslt/_draft/plots/preds-18sv9-ss.png',
       width = 8, height = 6, units = 'in', dpi = 300)
