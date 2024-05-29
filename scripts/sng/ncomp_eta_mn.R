library(compositions)
library(pls)
library(spls)
library(tidyverse)
library(viridis)
library(rgl)

load('data/ncog-18s-processed-2024-05-19.RData')
load('data/ceta-density-processed-2024-05-28.RData')

# combine seasonally adjusted density estimates and seasonally adjusted edna data
whales <- inner_join(log_density_estimates_adj, edna_clr_adj, by = 'cruise')

x <- dplyr::select(whales, starts_with('asv'))
y <- pull(whales, mn) 
n <- length(y)

# num of components and sparsity grid (ncomp and eta)
ncomp_grid <- seq(1, 15, by=1)
eta_grid <- seq(0.01, 0.95, length = 50)

ncomp_eta_gs <- data.frame(0,0,0,0,0,0)
names(ncomp_eta_gs) = c("ncomp", "eta", "r2", "adj_r2", "mspe", "nsel")

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
    ncomp_eta_gs <- rbind(ncomp_eta_gs, setNames(newRow,names(ncomp_eta_gs)))
    
    print(c(ncomp, eta))
  }
}

save(ncomp_eta_gs, file = paste('rslt/comp/ncomp_eta_gs_mn_NEW_', lubridate::today(), '.RData', sep = ''))

# get rid of row of 0s
ncomp_eta_gs <-  ncomp_eta_gs|> 
  slice(-1)

ncomp_eta_gs |> 
  filter(adj_r2 == max(adj_r2) | mspe == min(mspe))

# highest adj_r2: ncomp = 7, eta = 0.758
# lowest mspe: ncomp = 8, eta = 0.68

# 3D graph
# adjusted r2
with(ncomp_eta_gs , plot3d(x = ncomp, y = eta, z = adj_r2))
# mspe
with(ncomp_eta_gs , plot3d(x = ncomp, y = eta, z = mspe))

# boxplots of ncomp vs adjusted r^2
# adjusted r2 starts to decline at around 10 => model gets less helpful after 10 components
ncomp_eta_gs |>
  ggplot(aes(x=ncomp, y = adj_r2)) + geom_boxplot(aes(group= ncomp)) +
  scale_color_viridis() + theme_bw()

# less variation in r2 for higher ncomp
# plateus around 10
ncomp_eta_gs |>
  ggplot(aes(x=ncomp, y = r2)) + geom_boxplot(aes(group=ncomp)) +
  scale_color_viridis() + theme_bw()


# scatterplot of eta vs adjusted r^2
# ncomp around 12-16 for optimal adjusted r2
ncomp_eta_gs |>
  ggplot(aes(x=eta, y = adj_r2, color= factor(ncomp))) + geom_point() + 
  scale_color_viridis(discrete = TRUE) 
# look into lower number of components (try to fit with lower number)
# advantages to go to 3-4 instead of 2?

# lower etas for better r2
ncomp_eta_gs |>
  ggplot(aes(x=ncomp, y = adj_r2, color= eta)) + geom_point() + 
  scale_color_viridis(discrete = FALSE) 

# adj r2 line plot
# looks like eta in "blue" range - 0.45-0.64 best
ggplot(ncomp_eta_gs, aes(x = ncomp, y = adj_r2, color = factor(eta))) +
  geom_line() +
  theme_bw() 

# mspe line plot
# looks like eta ~ 0.63-0.68 for ncpmp = 7 has lowest mspe
ggplot(ncomp_eta_gs, aes(x = ncomp, y = mspe, color = factor(eta))) +
  geom_line() +
  theme_bw()

