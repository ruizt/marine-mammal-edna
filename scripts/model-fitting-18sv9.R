library(tidyverse)
library(spls)
library(pls)
library(modelr)
library(magrittr)
library(fs)
out_dir <- 'rslt/models/'
dir_create(out_dir)

## DATA INPUTS -----------------------------------------------------------------

# read in and combine 18sv9 data and scaled sightings
load('data/processed/ncog18sv9-2024-07-27.RData')
load('data/processed/mm-sightings-2024-07-27.RData')
whales <- inner_join(sightings, edna, by = 'cruise') 

# dimensions
p <- whales |> dplyr::select(starts_with('asv')) |> ncol()
n <- nrow(whales)

## STABILITY SELECTION ---------------------------------------------------------

# read in selection frequencies and fit metrics from LOOCV
sel_freq <- read_rds('rslt/loocv/2024-07-27/selection-frequencies.rds')
metrics <- read_rds('rslt/loocv/2024-07-27/metrics.rds')
loo_partitions <- read_rds('rslt/loocv/2024-07-27/partitions.rds')

# chosen stability threshold (minimax selection prob.)
pimax <- 0.8

# # plot expected selection error (EV) against average no. of selected asvs (q) 
# # across hyperparameter region
# curve(x^2/(p*(2*pimax - 1)), from = 1, to = 200,
#       xlab = 'q', ylab = 'EV')

# criterion: range of ncomp/eta to search should have an average no. of selected
# asvs not exceeding 100; controls EV at 1
EVmax <- 1
qmax <- sqrt((2*pimax - 1)*p*EVmax)
qmax

# criterion: range to search should not include regions having models with too
# few asv's selected
qmin <- 10

# using loocv results, plot q as a function of eta bin by ncomp
# visually search for ranges where qmin < q < qmax for a fixed ncomp
pal.grad <- colorRampPalette(c('red', 'blue'))
pal <- pal.grad(10)
metrics |>
  dplyr::select(species, ncomp, eta, n.asv) |>
  mutate(eta.bin = cut_interval(eta, n = 8)) |>
  group_by(species, ncomp, eta.bin) |>
  summarize(avg.n.asv = mean(n.asv),
            sd.n.asv = sd(n.asv)) |>
  ggplot(aes(x = eta.bin, y = avg.n.asv, color = factor(ncomp))) +
  facet_wrap(~species) +
  geom_point() +
  geom_path(aes(group = ncomp)) +
  geom_hline(yintercept = c(qmax, qmin), linetype = 3) +
  scale_y_log10() +
  scale_color_manual(values = pal) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90))

# check q for proposed region
metrics |>
  dplyr::select(species, ncomp, eta, n.asv) |>
  filter(ncomp == 6, 
         eta >= 0.625, eta <= 0.95) |>
  group_by(species, ncomp) |>
  summarize(avg.n = mean(n.asv), 
            sd.n = sd(n.asv)) |>
  slice_max(avg.n)

# upper bound on expected no. false positives: EV < 0.3
27^2/(p*(2*pimax - 1))

# find stable sets for each species
stable_sets <- sel_freq |> 
  filter(eta >= 0.625,
         eta <= 0.95,
         ncomp == 6,
         n >= 20) |> 
  group_by(species) |> 
  distinct(asv) |>
  nest(stable.set = asv) 

# check sizes
stable_sets

## FIT MODELS ON STABLE SETS ---------------------------------------------------

# fix number of components
ncomp <- 6

# subset data columns using stable sets and fit models for each species
fit <- whales |> 
  dplyr::select(-cruise) |>
  pivot_longer(c(bm, bp, mn), names_to = 'species', values_to = 'y') |>
  group_by(species) |>
  nest(data = -species) |>
  left_join(stable_sets) |>
  transmute(data = map2(data, stable.set, 
                        ~dplyr::select(.x, y, all_of(pull(.y, asv))))) |>
  mutate(fit = map(data, ~plsr(y ~ ., data = .x, ncomp = ncomp, 
                               scale = F, center = T)))

# compute fit metrics
ss_means_long <- ss_means |>
  pivot_longer(-season, names_to = 'species', values_to = 'seasonal.mean') |>
  mutate(species = str_remove(species, 'log.') |> str_remove('.imp.mean'))

sightings_raw_long <- sightings_raw |>
  pivot_longer(c(bp, bm, mn), names_to = 'species', values_to = 'ss.obs')

fit_df <- fit |>
  mutate(y = map(data, ~pull(.x, y)),
         fitted = map(fit, ~fitted(.x)[, , paste(ncomp, 'comps')]),
         cruise = pull(whales, cruise) |> list()) |>
  select(cruise, species, y, fitted) |>
  unnest(everything()) |>
  left_join(sightings_raw_long, join_by(species, cruise)) |>
  left_join(ss_means_long, join_by(season, species)) |>
  mutate(ss.fit = exp(fitted + seasonal.mean),
         lr.resid = y - fitted,
         ss.resid = ss.obs - ss.fit) 


fit_metrics <- fit_df |>
  summarize(adj.rsq.lr = (1 - ((n - 1)/(n - ncomp - 1))*var(lr.resid)/var(y)),
         adj.rsq.ss = (1 - ((n - 1)/(n - ncomp - 1))*var(ss.resid)/var(ss.obs))) |>
  left_join(stable_sets, join_by(species)) |>
  mutate(n.asv = map(stable.set, nrow)) |>
  unnest(n.asv) |>
  select(-stable.set)

fit_metrics

## LEAVE-ONE-OUT PREDICTIONS ---------------------------------------------------

# fit models on stable sets to leave one out partitions and compute predictions
loo_preds <- loo_partitions %>%
  expand_grid(species = c('bm', 'bp', 'mn')) |>
  mutate(train = map2(train, species, ~select(.x, {.y}, starts_with('asv')) |>
                        rename(y = {.y})),
         test = map2(test, species, ~select(.x, {.y}, starts_with('asv')) |>
                       rename(y = {.y}))) |>
  left_join(stable_sets, by = 'species') |>
  mutate(train = map2(train, stable.set, ~select(.x, y, any_of(.y$asv))),
         test = map2(test, stable.set, ~select(.x, y, any_of(.y$asv)))) |>
  select(test.cruise, species, train, test) |>
  mutate(fit = map(train, ~plsr(y ~ ., data = .x, ncomp = ncomp, scale = F, center = T)),
         pred = map2(fit, test, ~predict(.x, .y)[, , paste(ncomp, 'comps')]),
         y = map(test, ~pull(.x, y))) |>
  unnest(c(pred, y)) 

# back-transform predictions to original scale
loo_pred_df <- loo_sightings |> 
  select(test.cruise, test.season, seasonal.means) |>
  unnest(seasonal.means) |>
  filter(test.season == season) |>
  pivot_longer(starts_with('log'), 
               names_to = 'test.species',
               values_to = 'train.seasonal.mean') |>
  mutate(test.species = str_remove(test.species, 'log.') |> str_remove('.mean')) |>
  select(-season) |>
  inner_join(loo_preds, join_by(test.cruise, test.species == species)) |>
  mutate(ss.pred = exp(pred + train.seasonal.mean),
         ss.imp = exp(y + train.seasonal.mean)) |>
  inner_join(sightings_raw_long,
             join_by(test.cruise == cruise, test.species == species)) |>
  select(-c(train, test, fit))

# compute prediction metrics  
pred_metrics <- loo_pred_df |>
  group_by(test.species) |>
  summarize(rmspe.lr = mean((y - pred)^2) |> sqrt(),
            cor.lr = cor(y, pred),
            rmspe.ss = mean((ss.obs - ss.pred)^2) |> sqrt(),
            cor.ss = cor(ss.obs, ss.pred),
            rmspe.ss.naive = mean((ss.obs - exp(train.seasonal.mean))^2) |> sqrt(),
            cor.ss.naive = cor(ss.obs, exp(train.seasonal.mean))) |>
  pivot_longer(-test.species) |>
  pivot_wider(names_from = test.species, values_from = value) |>
  separate(name, into = c('metric', 'scale', 'model')) |>
  mutate(model = replace_na(model, 'pls')) |>
  arrange(scale, metric, model)


save(list = c('fit', 'loo_preds', 'loo_pred_df', 'fit_df', 'fit_metrics', 'pred_metrics'),
     file = paste(out_dir, 'fitted-models-18sv9-', today(), '.RData', sep = ''))

# 
# ## COMPARE WITH PREDICTION-OPTIMAL SPLS ----------------------------------------
# 
# # find prediction optimal spls models from loocv metrics
# loo_best <- metrics |>
#   group_by(ncomp, eta, species) |>
#   summarize(mspe = mean(sq.pred.err),
#             cor = cor(pred, pred + pred.err),
#             df = mean(n.asv)) |>
#   group_by(species) |>
#   slice_min(mspe)
# 
# # retrieve prediction optimal hyperparameter configurations
# parms <- loo_best |>
#   dplyr::select(species, ncomp, eta) |>
#   nest(parms = c(ncomp, eta))
# 
# # function to fit spls models
# fit_spls_fn <- function(.data, .parms){
#   x <- dplyr::select(.data, -y)
#   y <- pull(.data, y)
#   out <- spls(x, y, K = .parms$ncomp, eta = .parms$eta,
#               scale.x = F, scale.y = F)
#   return(out)
# }
# 
# # fit models
# fit_pos <- whales |> 
#   dplyr::select(-cruise) |>
#   pivot_longer(c(bm, bp, mn), names_to = 'species', values_to = 'y') |>
#   group_by(species) |>
#   nest(data = -species) |>
#   left_join(parms) |>
#   mutate(fit = map2(data, parms, fit_spls_fn),
#          df = map(fit, ~nrow(.x$projection)),
#          fitted = map(fit, ~predict(.x, type = 'fit')[, 1]),
#          y = map(data, ~pull(.x, y)),
#          adj.rsq = map2(y, fitted, ~(1 - (24/16)*var(.x - .y)/var(.y)))) 
# 
# # compute fit metrics
# loo_best |>
#   mutate(rmse = sqrt(mspe),
#          pred.corr = cor) |>
#   dplyr::select(species, rmse, pred.corr) |>
#   left_join(unnest(fit_pos, c(df, adj.rsq))) |>
#   dplyr::select(species, rmse, df, adj.rsq, pred.corr)
# 
# # compare with stability selection
# metrics_ss

