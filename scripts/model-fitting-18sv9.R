library(tidyverse)
library(spls)
library(pls)
library(modelr)
library(magrittr)
library(fs)
out_dir <- 'rslt/18sv9/'
dir_create(out_dir)

## DATA INPUTS -----------------------------------------------------------------

# read in and combine 18sv9 data and scaled sightings
load('data/processed/ncog18sv9-2024-07-20.RData')
load('data/processed/mm-sightings-2024-07-20.RData')
whales <- inner_join(sightings, edna, by = 'cruise') 

# dimensions
p <- whales |> dplyr::select(starts_with('asv')) |> ncol()
n <- nrow(whales)

## STABILITY SELECTION ---------------------------------------------------------

# read in selection frequencies and fit metrics from LOOCV
sel_freq <- read_rds('rslt/loocv/18sv9/selection-frequencies.rds')
metrics <- read_rds('rslt/loocv/18sv9/metrics.rds')

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
  filter(ncomp == 8, 
         eta >= 0.65, eta <= 0.95) |>
  group_by(species, ncomp) |>
  summarize(avg.n = mean(n.asv), 
            sd.n = sd(n.asv)) |>
  slice_max(avg.n)

# upper bound on expected no. false positives: EV < 0.3
35^2/(p*(2*pimax - 1))

# find stable sets for each species
stable_sets <- sel_freq |> 
  filter(eta >= 0.65,
         eta <= 0.95,
         ncomp == 8,
         n >= 20) |> 
  group_by(species) |> 
  distinct(asv) |>
  nest(stable.set = asv) 

# check sizes
stable_sets

## FIT MODELS ON STABLE SETS ---------------------------------------------------

# fit_bm <- stable_sets |>
#   unnest(stable.set) |>
#   filter(species == 'bm') |>
#   pull(asv) %>%
#   {dplyr::dplyr::select(whales, bm, all_of(.))} %>%
#   plsr(bm ~ ., data = ., ncomp = 8,
#        scale = F, center = T, validation = 'LOO')
# 
# # loocv root mspe
# rmse <- RMSEP(fit_bm)$val['CV', , paste(8, 'comps')] 
# 
# # loo predictions
# loo_preds <- fit_bm$validation$pred[, , paste(8, 'comps')]
# 
# # fitted values
# fitted <- fitted(fit_bm)[, , paste(8, 'comps')]
# 
# # coefficients
# coef(fit_bm, intercept = T)[, , paste(8, 'comps')]

# subset data columns using stable sets and fit models for each species
fit_ss <- whales |> 
  dplyr::select(-cruise) |>
  pivot_longer(c(bm, bp, mn), names_to = 'species', values_to = 'y') |>
  group_by(species) |>
  nest(data = -species) |>
  left_join(stable_sets) |>
  transmute(data = map2(data, stable.set, 
                        ~dplyr::select(.x, y, all_of(pull(.y, asv))))) |>
  mutate(fit = map(data, ~plsr(y ~ ., data = .x, ncomp = 8, 
                               scale = F, center = T, validation = 'LOO')))

# compute fit/prediction metrics
metrics_ss <- fit_ss |>
  mutate(rmse = map(fit, ~RMSEP(.x)$val['CV', , paste(8, 'comps')]),
         loo.preds = map(fit, ~.x$validation$pred[, , paste(8, 'comps')]),
         fitted = map(fit, ~fitted(.x)[, , paste(8, 'comps')]),
         y = map(data, ~pull(.x, y)),
         df = map(data, ~ncol(.x) - 1),
         adj.rsq = map2(y, fitted, ~(1 - (24/16)*var(.x - .y)/var(.y))),
         pred.corr = map2(y, loo.preds, ~cor(.x, .y))) |>
  unnest(c(rmse, df, adj.rsq, pred.corr)) |>
  dplyr::select(species, where(is.numeric))

save(list = c('fit_ss', 'metrics_ss'),
     file = paste(out_dir, 'fitted-models-18sv9-', today(), '.RData', sep = ''))

# # prediction vs observation
# fit_ss |>
#   transmute(loo.preds = map(fit, ~.x$validation$pred[, , paste(8, 'comps')]),
#             y = map(data, ~pull(.x, y))) |>
#   unnest(everything()) |>
#   ggplot(aes(x = y, y = loo.preds)) +
#   facet_wrap(~species) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0)
# 
# # prediction vs fitted
# fit_ss |>
#   transmute(fitted = map(fit, ~fitted(.x)[, , paste(8, 'comps')]),
#             y = map(data, ~pull(.x, y))) |>
#   unnest(everything()) |>
#   ggplot(aes(x = y, y = fitted)) +
#   facet_wrap(~species) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0)
# 
# # residual vs fit
# fit_ss |>
#   transmute(fitted = map(fit, ~fitted(.x)[, , paste(8, 'comps')]),
#             y = map(data, ~pull(.x, y))) |>
#   unnest(everything()) |>
#   mutate(resid = y - fitted) |>
#   ggplot(aes(y = resid, x = fitted)) +
#   facet_wrap(~species) +
#   geom_point() +
#   geom_hline(yintercept = 0)
# 
# # residual autocorrelation
# pacf_fn <- function(x){
#   pacf_out <- pacf(x, plot = F)
#   out <- bind_cols(pacf = pacf_out$acf[, 1, 1],
#                    lag = pacf_out$lag[, , 1],
#                    se = 2/sqrt(pacf_out$n.used))
#   return(out)
# }
# fit_ss |>
#   transmute(fitted = map(fit, ~fitted(.x)[, , paste(8, 'comps')]),
#             y = map(data, ~pull(.x, y))) |>
#   unnest(everything()) |>
#   transmute(resid = y - fitted) |>
#   nest(cols = resid) |>
#   mutate(pacf = map(cols, pacf_fn)) |>
#   unnest(pacf) %>%
#   ggplot(aes(x = lag)) +
#   facet_wrap(~species) +
#   geom_linerange(aes(ymin = 0, ymax = pacf)) +
#   geom_ribbon(aes(ymin = -se, ymax = se), 
#               fill = 'blue', 
#               alpha = 0.1)


## COMPARE WITH PREDICTION-OPTIMAL SPLS ----------------------------------------

# find prediction optimal spls models from loocv metrics
loo_best <- metrics |>
  group_by(ncomp, eta, species) |>
  summarize(mspe = mean(sq.pred.err),
            cor = cor(pred, pred + pred.err),
            df = mean(n.asv)) |>
  group_by(species) |>
  slice_min(mspe)

# retrieve prediction optimal hyperparameter configurations
parms <- loo_best |>
  dplyr::select(species, ncomp, eta) |>
  nest(parms = c(ncomp, eta))

# function to fit spls models
fit_spls_fn <- function(.data, .parms){
  x <- dplyr::select(.data, -y)
  y <- pull(.data, y)
  out <- spls(x, y, K = .parms$ncomp, eta = .parms$eta,
              scale.x = F, scale.y = F)
  return(out)
}

# fit models
fit_pos <- whales |> 
  dplyr::select(-cruise) |>
  pivot_longer(c(bm, bp, mn), names_to = 'species', values_to = 'y') |>
  group_by(species) |>
  nest(data = -species) |>
  left_join(parms) |>
  mutate(fit = map2(data, parms, fit_spls_fn),
         df = map(fit, ~nrow(.x$projection)),
         fitted = map(fit, ~predict(.x, type = 'fit')[, 1]),
         y = map(data, ~pull(.x, y)),
         adj.rsq = map2(y, fitted, ~(1 - (24/16)*var(.x - .y)/var(.y)))) 

# compute fit metrics
loo_best |>
  mutate(rmse = sqrt(mspe),
         pred.corr = cor) |>
  dplyr::select(species, rmse, pred.corr) |>
  left_join(unnest(fit_pos, c(df, adj.rsq))) |>
  dplyr::select(species, rmse, df, adj.rsq, pred.corr)

# compare with stability selection
metrics_ss

