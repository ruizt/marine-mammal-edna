library(compositions)
library(pls)
library(spls)
library(tidyverse)
library(viridis)
library(rgl)

load('data/ncog-18s-processed-2024-05-19.RData')
load('data/ceta-density-processed-2024-05-19.RData')

# combine seasonally adjusted density estimates and seasonally adjusted edna data
whales <- inner_join(log_density_estimates_adj, edna_clr_adj, by = 'cruise')

x <- dplyr::select(whales, starts_with('asv'))
y <- pull(whales, bp) 
n <- length(y)

# num of components and sparsity grid (ncomp and eta)
ncomp_grid <- seq(1, 15, by=1)
eta_grid <- seq(0.01, 0.95, length = 50)

ncomp_eta_bp <- data.frame(0,0,0,0,0,0)
names(ncomp_eta_bp) = c("ncomp", "eta", "r2", "adj_r2", "mspe", "nsel")

# loop over diff vals of ncomp
for (ncomp in ncomp_grid) {
  for (eta in eta_grid) {
    
    # fit model to full dataset for given eta/ncomp    
    fit <- spls(x, y, K = ncomp, eta = eta, scale.x = F, scale.y = F)
    
    # r^2 and adjusted R2
    fitted <- predict(fit, type = 'fit')
    resid <- y - fitted
    r2 <- 1 - (var(resid)/var(y))
    adj_r2 <- 1 - ((1-r2)*(n-1))/(n-ncomp-1)
    
    # number of asvs selected
    nsel <- nrow(fit$projection)
    
    ## LEAVE ONE OUT PREDICTIONS
    
    # storage for leave one out predictions
    loo_preds <- rep(NA, nrow(x))
    
    for (i in 1:nrow(x)){
      x_train <- x[-i, ] # removes ith row
      y_train <- y[-i]  # removes ith element
      
      
      # fit cross validation model
      fit_cv <- spls(x_train, y_train, K = ncomp, eta = eta, 
                     scale.x = F, scale.y = F)
      
      # predict response variable
      loo_preds[i] <- predict(fit_cv, newx = x[i, ], type = "fit")
    }
    
    # mspe (mean squared prediction error)
    mspe <- mean((y - loo_preds)^2) 
    
    ## ALTERNATIVE USING CV.SPLS
    cv_out <- cv.spls(x, y, fold = n, K = ncomp, eta = eta, scale.x = F, scale.y = F, plot.it = F)
    mspe <- cv_out$mspemat[1]
    
    
    # print(paste("ncomp:", ncomp, "eta:", eta))
    # print(paste("adj_r2:", adj_r2))
    # print(paste("mspe:", mspe))
    
    newRow <- c(ncomp, eta, r2, adj_r2 , mspe, nsel)
    ncomp_eta_bp <- rbind(ncomp_eta_bp, setNames(newRow,names(ncomp_eta_bp)))
    
    print(c(ncomp, eta))
  }
}

# large variability w mspe

save(ncomp_eta_bp, file = paste('rslt/comp/ncomp_eta_gs_bp_NEW_', lubridate::today(), '.RData', sep = ''))

# get rid of row of 0s
ncomp_eta_bp <-  ncomp_eta_bp|> 
  slice(-1)

ncomp_eta_bp |> 
  filter(adj_r2 == max(adj_r2) | mspe == min(mspe))

# highest adj_r2: ncomp = 8, eta = 0.336
# lowest mspe: ncomp = 2, eta =0.259

# 3D graph
# adjusted r2
with(ncomp_eta_bp , plot3d(x = ncomp, y = eta, z = adj_r2))
# mspe
with(ncomp_eta_bp , plot3d(x = ncomp, y = eta, z = mspe))

# boxplots of ncomp vs adjusted r^2
# adjusted r2 starts to decline at around 8 => model gets less helpful after 8 components
# if you add in enough components, sparisty affects adjr2 less
ncomp_eta_bp |>
  ggplot(aes(x=ncomp, y = adj_r2)) + geom_boxplot(aes(group= ncomp)) +
  scale_color_viridis() + theme_bw() + 
  labs(title="Fin Whales - Adj R^2 vs ncomp")

# once above 5-6 comp, start overfitting
ncomp_eta_bp |>
  ggplot(aes(x=ncomp, y = mspe)) + geom_boxplot(aes(group=ncomp)) +
  scale_color_viridis() + theme_bw()+ 
  labs(title="Fin Whales - MSPE vs ncomp")

ncomp_eta_bp |>
  ggplot(aes(x=ncomp, y = nsel)) + geom_boxplot(aes(group=ncomp)) +
  scale_color_viridis() + theme_bw() + 
  labs(title="Fin Whales - # Selected ASVs vs ncomp")


# scatterplot of eta vs adjusted r^2
# ncomp around 9-11 for optimal adjusted r2
ncomp_eta_bp |>
  ggplot(aes(x=eta, y = adj_r2, color= factor(ncomp))) + geom_point() + 
  scale_color_viridis(discrete = TRUE) 

# eta around 0.25-0.75 for ncomp ~ 8
ncomp_eta_bp |>
  ggplot(aes(x=ncomp, y = adj_r2, color= eta)) + geom_point() + 
  scale_color_viridis(discrete = FALSE) 

# adj r2 line plot
# looks like eta in "green" range - 0.25-0.35 best for ncomp ~ 8 best
ggplot(ncomp_eta_bp, aes(x = ncomp, y = adj_r2, color = factor(eta))) +
  geom_line() +
  theme_bw() + labs(title="Fin Whale Adj R^2 Lineplot")

# mspe line plot
# looks like eta ~ 0.25 for ncomp ~ 2 best
# more variability w/ sparser model
ggplot(ncomp_eta_bp, aes(x = ncomp, y = mspe, color = factor(eta))) +
  geom_line() +
  theme_bw() + labs(title="Fin Whale MSPE Lineplot")

