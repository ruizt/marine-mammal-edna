library(tidyverse)
library(lubridate)
library(patchwork)
library(ggmap)
library(mapdata)
library(sf)

## DIRECTORIES -----------------------------------------------------------------
data_dir <- 'data/processed/'
model_dir <- 'rslt/models/scaled-sightings/'
stbl_dir <- 'rslt/stability-selection/'
val_dir <- 'rslt/nested-validation/'

out_dir <- 'rslt/fig/'
fs::dir_create(out_dir)

## FIGURE 1: SAMPLING MAPS -----------------------------------------------------

# retrieve station metadata from NCOG
paste(data_dir, 'ncog16s.RData', sep = '') |> load()
ncog16s_meta <- sample_metadata |> 
  select(sta.id, cruise, lat.dec, lon.dec) |>
  mutate(marker = "16S")

paste(data_dir, 'ncog18sv4.RData', sep = '') |> load()
ncog18sv4_meta <- sample_metadata |> 
  select(sta.id, cruise, lat.dec, lon.dec) |>
  mutate(marker = "18SV4")

paste(data_dir, 'ncog18sv9.RData', sep = '') |> load()
ncog18sv9_meta <- sample_metadata |> 
  select(sta.id, cruise, lat.dec, lon.dec) |>
  mutate(marker = "18SV9")

ncog_meta <- bind_rows(ncog16s_meta, ncog18sv4_meta, ncog18sv9_meta) |>
  separate(sta.id, into = c('line', 'sta'), sep = ' ')

# station locations
stations <- ncog_meta |>
  group_by(line, sta) |>
  summarize(lat = mean(lat.dec), 
            long = mean(lon.dec), 
            .groups = 'drop') |>
  arrange(line, desc(lat)) |>
  filter(as.numeric(line) > 74)

# coordinates for labels
line_labels <- stations |>
  group_by(line) |>
  slice_min(order_by = lat, n = 1)

# bounding box
bbox <- c(xmin = min(stations$long) - 2, 
          xmax = max(stations$long) + 1, 
          ymin = min(stations$lat) - 1,
          ymax = max(stations$lat) + 1)

# read in shapefile and crop to coastline
gshhg <- read_sf("data/_gshhg-shp-2.3.7/GSHHS_shp/f/GSHHS_f_L1.shp") 
sf_use_s2(F)
gshhg_crop <- gshhg |> st_crop(bbox)

# plot ncog sampling locations
ncog_map <- ggplot(aes(x = long, y = lat), data = stations) +
  geom_path(aes(group = line), linewidth = 0.3) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_text(aes(label = line),
            hjust = 1, 
            angle = 29, 
            nudge_x = -0.25, 
            nudge_y = -0.25*0.29,
            data = line_labels) +
  geom_sf(color = 'black', linewidth = 0.25,
          data = gshhg_crop, inherit.aes = F) +
  scale_x_continuous(limits = bbox[1:2], expand = expansion(0, 0)) +
  scale_y_continuous(limits = bbox[3:4], expand = expansion(0, 0)) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = 'black', linewidth = 0.1),
        panel.grid.minor = element_blank()) +
  labs(x = NULL, y = NULL)

# read in sighting locations
visual_sightings <- read_csv('data/_raw/CalCOFI_2004-2021_CombinedSightings.csv') |>
  rename_with(~str_remove_all(.x, '[:punct:]') |> 
                str_squish() |> 
                str_replace_all(' ', '.') |> 
                tolower()) |>
  mutate(species = tolower(species.1)) |>
  filter(adjusted.both.on.effort.and.on.transect == 'ON',
         species %in% c('bm', 'bp', 'mn')) |>
  mutate(datetime = mdy_hm(datetime.local),
         year = year(datetime)) |>
  select(cruise, year, datetime, species, declat, declong) |>
  filter(year >= 2014) |>
  distinct(declat, declong) |>
  rename(lat = declat, long = declong)

# plot sighting locations
sighting_map <- ggplot(aes(x = long, y = lat), data = visual_sightings) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_path(aes(group = line), data = filter(stations, as.numeric(line) > 75)) +
  geom_text(aes(label = line),
            hjust = 1, 
            angle = 29, 
            nudge_x = -0.25, 
            nudge_y = -0.25*0.29,
            data = filter(line_labels, as.numeric(line) > 75)) +
  geom_sf(color = 'black', linewidth = 0.25,
          data = gshhg_crop, inherit.aes = F) +
  scale_x_continuous(limits = bbox[1:2], expand = expansion(0, 0)) +
  scale_y_continuous(limits = bbox[3:4], expand = expansion(0, 0)) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = 'black', linewidth = 0.1),
        panel.grid.minor = element_blank()) +
  labs(x = NULL, y = NULL)

fig_map <- ncog_map + sighting_map + 
  plot_layout(nrow = 1) + 
  plot_annotation(tag_levels = 'A')

ggsave(fig_map, file = paste(out_dir, 'fig1-map.png', sep = ''),
       width = 8, height = 4, units = 'in', dpi = 400)

ncog_map_slide <- ncog_map + labs(title = 'A. eDNA sampling')
sighting_map_slide <- sighting_map + labs(title = 'B. Whale sightings')
fig_map_slide <- ncog_map_slide + sighting_map_slide +
  plot_layout(nrow = 1)

ggsave(fig_map_slide, file = paste(out_dir, 'fig1-map-slide.png', sep = ''),
       width = 8, height = 4, units = 'in', dpi = 400)

## FIGURE 2: TIME SERIES -------------------------------------------------------

# load sighting data
paste(data_dir, 'mm-sightings.RData', sep = '') |> load()

sighting_data <- sightings_raw |>
  mutate(cruise.ym = ym(cruise)) |>
  pivot_longer(c(bp, bm, mn), names_to = 'species', values_to = 'ss') |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales')),
         season = factor(season, levels = c('winter', 'spring', 'summer', 'fall'))) 

seasonal_means <- ss_means |>
  mutate(across(-season, exp)) |>
  rename(bp = log.bp.imp.mean,
         bm = log.bm.imp.mean,
         mn = log.mn.imp.mean) |>
  pivot_longer(-season, names_to = 'species', values_to = 'ss') |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales')),
         season = factor(season, levels = c('winter', 'spring', 'summer', 'fall'))) |>
  arrange(species, season)

p1 <- sighting_data |>
  ggplot(aes(x = cruise.ym, y = ss)) +
  facet_wrap(~species, nrow = 3) +
  geom_path() +
  theme_bw() +
  labs(x = NULL, y = 'Sightings per 1000km', title = NULL) +
  theme(panel.grid = element_line(linewidth = 0.1, color = 'black'),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 11))


p2 <- sightings_raw |>
  mutate(cruise.ym = ym(cruise)) |>
  pivot_longer(c(bp, bm, mn), names_to = 'species', values_to = 'ss') |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales')),
         season = factor(season, levels = c('winter', 'spring', 'summer', 'fall'))) |>
  ggplot(aes(x = season, y = ss)) +
  facet_wrap(~species, nrow = 3) +
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) +
  geom_point(color = 'red', shape = 3, data = seasonal_means) +
  geom_path(aes(group = species), color = 'red', linewidth = 0.2, data = seasonal_means) +
  theme_bw() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(panel.grid = element_line(linewidth = 0.1, color = 'black'),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 11),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p1 + p2 + plot_layout(ncol = 2, widths = c(2, 1))

paste(out_dir, 'fig2-timeseries.png', sep = '') |>
  ggsave(width = 5, height = 4, dpi = 400, units = 'in')

paste(out_dir, 'fig2-timeseries-slide.png', sep = '') |>
  ggsave(width = 6, height = 3, dpi = 400, units = 'in')

## FIGURE 3: PREDICTIONS ---------------------------------------------------------

# predictions from 16s
pred_pts_16s <- paste(stbl_dir, '16s-ss/loo-preds.rds', sep = '') |>
  read_rds() |>
  mutate(cruise.ym = ym(obs.id),
         species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales')),
         marker = "16S") |>
  arrange(species, cruise.ym)

# predictions from 18sv4
pred_pts_18sv4 <- paste(stbl_dir, '18sv4-ss/loo-preds.rds', sep = '') |>
  read_rds() |>
  mutate(cruise.ym = ym(obs.id),
         species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales')),
         marker = "18SV4") |>
  arrange(species, cruise.ym)

# predictions from 18sv9
pred_pts_18sv9 <- paste(stbl_dir, '18sv9-ss/loo-preds.rds', sep = '') |>
  read_rds() |>
  mutate(cruise.ym = ym(obs.id),
         species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales')),
         marker = "18SV9") |>
  arrange(species, cruise.ym)

# combine predictions from each marker
pred_pts <- bind_rows(pred_pts_16s, pred_pts_18sv4, pred_pts_18sv9) |>
  mutate(key = paste('Predicted (', marker, ')', sep = '')) |>
  select(species, cruise.ym, pred.ss, key, pred.ss.qlo, pred.ss.qhi) |>
  rename(value = pred.ss)

# observations
obs_pts <- pred_pts_18sv9 |>
  select(species, cruise.ym, obs.ss) |>
  mutate(key = 'Observed') |>
  rename(value = obs.ss)

# color palette
pal <- c('#000000', RColorBrewer::brewer.pal(name = 'Dark2', n = 3))

# predictions from each marker overlaid on serial observations
p1 <- bind_rows(pred_pts, obs_pts) |>
  ggplot(aes(x = cruise.ym,
             y = value,
             color = key,
             linetype = key)) +
  facet_wrap(~species, nrow = 3, scale = 'fixed') +
  geom_path(linewidth = 0.4) +
  # geom_ribbon(inherit.aes = F,
  #             aes(x = cruise.ym,
  #                 ymin = pred.ss.qlo,
  #                 ymax = pred.ss.qhi,
  #                 fill = key),
  #             alpha = 0.2) +
  scale_x_date(date_breaks = '1 years', labels = ~year(.x)) +
  scale_color_manual(values = pal) +
  # scale_y_sqrt() +
  theme_bw() +
  # geom_blank(data = pivot_longer(fit_pts_df, c(Observed, Fitted)),
  #            aes(linetype = NULL)) +
  guides(linetype = guide_legend(title = NULL, position = 'top', nrow = 2),
         color = guide_legend(title = NULL, position = 'top', nrow = 2),
         fill = guide_none()) +
  labs(x = 'Year', y = 'Sightings per 1000km', title = NULL) +
  theme(panel.grid = element_line(linewidth = 0.1, color = 'black'),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 11))

# predictions vs observations
p2 <- bind_rows(pred_pts_16s, pred_pts_18sv4, pred_pts_18sv9) |>
  mutate(species = str_remove_all(species, ' whales')) |>
  ggplot(aes(x = obs.ss, y = pred.ss)) +
  facet_grid(species~marker) +
  geom_point(size = 0.8) +
  geom_linerange(aes(ymin = pred.ss.qlo, ymax = pred.ss.qhi),
                 linewidth = 0.4) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.2) +
  theme_bw() +
  theme(panel.grid = element_line(linewidth = 0.1, color = 'black'),
        panel.grid.minor = element_blank(),
        text = element_text(size = 11)) +
  labs(x = 'Observed value', y = 'Predicted value')

# correlations between predictions and observations (for labeling)
pred_cors <- bind_rows(pred_pts_16s, pred_pts_18sv4, pred_pts_18sv9) |>
  mutate(species = str_remove_all(species, ' whales')) |>
  group_by(species, marker) |>
  summarize(r = cor(obs.ss, pred.ss),
            obs.ss = 30,
            pred.ss = 0.05,
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

# arrange panel layout
fig_predictions <- p1 + p2_ann +
  plot_layout(nrow = 1, ncol = 2, widths = c(1, 1))

# export
ggsave(fig_predictions, filename = paste(out_dir, 'fig3-predictions.png', sep = ''),
       width = 6.5, height = 4, units = 'in', dpi = 400)


## SUPPLEMENTARY FIGURE 1: MODEL DIAGNOSTICS -----------------------------------

# predictions from 18sv9
load('rslt/models/scaled-sightings/fitted-models-18sv9.RData')
fit_pts_18sv9 <- fit_df |>
  mutate(cruise.ym = ym(cruise),
         year = year(cruise.ym),
         species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue', 'Fin', 'Humpback'))) |>
  select(cruise, cruise.ym, species, lr.obs, lr.fit, lr.resid) |>
  arrange(species, cruise.ym) |>
  mutate(marker = '18SV9')

# predictions from 18sv4
load('rslt/models/scaled-sightings/fitted-models-18sv4.RData')
fit_pts_18sv4 <- fit_df |>
  mutate(cruise.ym = ym(cruise),
         year = year(cruise.ym),
         species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue', 'Fin', 'Humpback'))) |>
  select(cruise, cruise.ym, species, lr.obs, lr.fit, lr.resid) |>
  arrange(species, cruise.ym) |>
  mutate(marker = '18SV4')

# predictions from 16s
load('rslt/models/scaled-sightings/fitted-models-16s.RData')
fit_pts_16s <- fit_df |>
  mutate(cruise.ym = ym(cruise),
         year = year(cruise.ym),
         species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue', 'Fin', 'Humpback'))) |>
  select(cruise, cruise.ym, species, lr.obs, lr.fit, lr.resid) |>
  arrange(species, cruise.ym) |>
  mutate(marker = '16S')

# merge above outputs
fit_pts <- bind_rows(fit_pts_16s, fit_pts_18sv4, fit_pts_18sv9)

# for rescaling axes without <scales = 'free_x'>
zscore <- function(x){(x - mean(x))/sd(x)}

# plot residuals vs fit for each model
resid_fit <- fit_pts |>
  group_by(species, marker) |>
  mutate(across(c(lr.fit, lr.resid), zscore)) |>
  ggplot(aes(x = lr.fit, y = lr.resid)) +
  facet_grid(species ~ marker) +
  # facet_wrap(species ~ marker, scales = 'free') +
  geom_point(alpha = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.2) +
  # geom_smooth(span = 1.5) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = 'black', linewidth = 0.1),
        panel.grid.minor = element_blank()) +
  labs(x = 'z(fitted)', y = 'z(residual)')

# function for residual autocorrelation
pacf_fn <- function(x){
  pacf_out <- pacf(x, plot = F, lag.max = 10)
  out <- bind_cols(pacf = c(1, pacf_out$acf[, 1, 1]),
                   lag = c(0, pacf_out$lag[, , 1]),
                   se = 2/sqrt(pacf_out$n.used))
  return(out)
}

# plot pacf for each model
resid_pacf <- fit_pts |>
  ungroup() |>
  select(species, marker, lr.resid) |>
  nest(resids = lr.resid, .by = c(species, marker)) |>
  mutate(pacf = purrr::map(resids, pacf_fn)) |>
  unnest(pacf) |>
  ggplot(aes(x = lag)) +
  facet_grid(species ~ marker) +
  scale_y_continuous(limits = c(-0.5, 1), n.breaks = 6) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  geom_linerange(aes(ymin = 0, ymax = pacf),
                 linewidth = 0.4) +
  geom_ribbon(aes(ymin = -se, ymax = se),
              fill = 'blue',
              alpha = 0.2) +
  theme_bw() +
  theme(panel.grid.major.y = element_line(color = 'black', linewidth = 0.1),
        panel.grid.minor.y = element_line(color = 'grey', linewidth = 0.1),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = 'Lag', y = 'Partial autocorrelation')

resid_pacf

# panel layout
resid_diagnostics <- resid_fit + resid_pacf + plot_layout(nrow = 1, widths = c(1, 1)) +
  plot_annotation(tag_levels = 'A')
ggsave(resid_diagnostics, filename = paste(out_dir, 'sfig1-resid.png', sep = ''),
       width = 8, height = 4, units = 'in', dpi = 400)



