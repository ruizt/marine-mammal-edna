library(compositions)
library(pls)
library(spls)
library(tidyverse)
library(viridis)
library(rgl)

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
y <- pull(whales, mn) 

# hyperparameter grid
grid_res <- 40
eta_grid_logspace <- seq(-2.5, -0.001, length = grid_res) |> exp()
eta_grid_linspace <- seq(0.001, 0.999, length = grid_res)

# loocv to choose hyperparameter -- initial pass
set.seed(32824)
cv_out_init <- cv.spls(x, y, 
                       K = 2, 
                       eta = eta_grid_logspace, 
                       fold = length(y), scale.x = F, scale.y = F,
                       plot.it = F)

plot(log(eta_grid_logspace), 
     cv_out_init$mspemat[,1], 
     type = 'l',
     xlab = expr(log(eta)),
     ylab = 'LOOCV MSPE') 
abline(v = log(cv_out_init$eta.opt), lty = 2)

# loocv to choose hyperparameter -- second pass
eta_init <- cv_out_init$eta.opt
mspe_grid_init <- cv_out_init$mspemat[, 1]
mspe_init <- mspe_grid_init[which(eta_grid_logspace == eta_init)]
mspe_init_se <- sd(mspe_grid_init)/sqrt(grid_res)
eta_range_ix <- which(mspe_grid_init <= mspe_init + 10*mspe_init_se)
eta_grid <- seq(from = eta_grid_logspace[min(eta_range_ix)],
                to = eta_grid_logspace[max(eta_range_ix)],
                length = grid_res)

cv_out <- cv.spls(x, y, 
                  K = 2, 
                  eta = eta_grid, 
                  fold = length(y), scale.x = F, scale.y = F,
                  plot.it = F)

eta_min <- cv_out$eta.opt
mspe_grid <- cv_out$mspemat[, 1]
mspe_min <- mspe_grid[which(eta_grid == eta_min)]
mspe_se <- sd(mspe_grid)/sqrt(grid_res)
eta_1se <- eta_grid[which(mspe_grid <= mspe_min + 5*mspe_se) |> max()]

plot(log(eta_grid), 
     cv_out$mspemat[,1], 
     type = 'l',
     xlab = expr(log(eta)),
     ylab = 'LOOCV MSPE') 
abline(h = mspe_min + 5*mspe_se, lty = 2)
abline(v = log(eta_1se), lty = 2)

# fit spls model
fit <- spls(x, y, K = 2, eta = eta_1se, scale.x = F, scale.y = F)
fit

# variability explained
fitted_vals <- predict(fit, x, type = 'fit')
1 - ((length(y) - 2)*var(y - fitted_vals))/((length(y) - 1)*var(y))

# estimate of prediction error
mspe <- mspe_grid[eta_grid == eta_1se]
sqrt(mspe) |> exp()

# # verify fitted model components
# x_sel <- dplyr::select(x, rownames(fit$projection))
# x_sel_proj <- as.matrix(x_sel) %*% fit$projection
# 
# lm(y ~ x_sel_proj)$fitted |> 
#   bind_cols(fitted_vals)
# 
# (fit$projection %*% lm(y ~ x_sel_proj)$coef[2:3]) |>
#   bind_cols(fit$betahat[fit$betahat != 0])

selected_asvs <- tibble(short.id = colnames(x)[fit$A],
                        coef = fit$betahat[fit$A]) |>
  left_join(dplyr::select(taxa, silva_Taxon, short.id), by = 'short.id') |>
  separate(silva_Taxon, into = c('d', 'p', 'c', 'o', 'f', 'g'), sep = ';') |>
  arrange(coef)

selected_asvs |>
  group_by(p) |>
  summarize(count = n(),
            mean.coef = mean(coef)) |>
  arrange(desc(count))

## MODEL SPARSITY ----------
x <- dplyr::select(whales, starts_with('asv'))
y <- pull(whales, mn) 

# num of components and sparsity grid (ncomp and eta)
ncomp_grid <- seq(1, 20, by=1)
eta_grid <- seq(0.1, 0.9, by=0.1)

# intialize vars to store results

ncomp_eta_gs <- data.frame(0,0,0,0,0)
names(ncomp_eta_gs) = c("ncomp", "eta", "r2", "adj_r2", "mspe")

# loop over diff vals of ncomp
for (ncomp in ncomp_grid) {
  for (eta in eta_grid) {
    # fit model
    # select eta param first before fitting model eta_1se
    
    model_fit <- spls(x, y, K = ncomp, eta = eta, scale.x = F, scale.y = F)

    
    # predict response variable
    pred <- predict(model_fit, newdata = x, type = "fit")
    #print(pred)
   
    
    # r^2
    r2 <- 1 - (var(pred)/var(y))
    
    n <- length(y)
    
    # adjusted r2
    adj_r2 <- 1 - ((1-r2)*(n-1))/(n-ncomp-1)
  
    
    # MSPE (mean squared prediction error)
    mspe <- mean((y - pred) ^2)
    
    # print(paste("ncomp:", ncomp, "eta:", eta))
    # print(paste("adj_r2:", adj_r2))
    # print(paste("mspe:", mspe))
    newRow <- c(ncomp, eta, r2, adj_r2 , mspe)
    #print(newRow)
    ncomp_eta_gs <- rbind(ncomp_eta_gs, setNames(newRow,names(ncomp_eta_gs)))
                                   
                                                      
  }
}

# get rid of row of 0s
ncomp_eta_gs <-  ncomp_eta_gs|> 
  slice(-1)

ncomp_eta_gs |> 
  filter(adj_r2 == min(adj_r2) | mspe == min(mspe))

## PLOTS ----------

# scatterplot of eta vs adjusted r^2
ncomp_eta_gs |>
  ggplot(aes(x=eta, y = adj_r2, color= factor(ncomp))) + geom_point() + 
  scale_color_viridis(discrete = TRUE) + theme_bw()

# scatterplot of eta vs mspe
ncomp_eta_gs |>
  ggplot(aes(x=eta, y = mspe, color= factor(ncomp))) + geom_point() + 
  scale_color_viridis(discrete = TRUE) + theme_bw()

# scatterplot of ncomp vs adjusted r^2
ncomp_eta_gs |>
  ggplot(aes(x=ncomp, y = adj_r2, color = eta)) + geom_point() +
  scale_color_viridis() + theme_bw()

# scatterplot of ncomp vs mspe
ncomp_eta_gs |>
  ggplot(aes(x=ncomp, y = mspe, color = eta)) + geom_point() +
  scale_color_viridis() + theme_bw()

# 3D graph
# adjusted r2
with(ncomp_eta_gs , plot3d(x = ncomp, y = eta, z = adj_r2))
# mspe
with(ncomp_eta_gs , plot3d(x = ncomp, y = eta, z = mspe))

