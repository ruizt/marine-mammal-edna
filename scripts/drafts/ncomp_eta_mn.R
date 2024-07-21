library(compositions)
library(pls)
library(spls)
library(tidyverse)
library(viridis)
library(rgl)

load('data/ncog-18s-processed-2024-07-04.RData')
load('data/ceta-density-processed-2024-07-04.RData')

# combine seasonally adjusted density estimates and seasonally adjusted edna data
whales <- inner_join(log_density_estimates_adj, edna_clr_adj, by = 'cruise')

x <- dplyr::select(whales, starts_with('asv'))
y <- pull(whales, mn) 
n <- length(y)

# num of components and sparsity grid (ncomp and eta)
ncomp_grid <- seq(1, 15, by=1)
eta_grid <- seq(0.01, 0.95, length = 50)

ncomp_eta_mn <- data.frame(0,0,0,0,0,0)
names(ncomp_eta_mn) = c("ncomp", "eta", "r2", "adj_r2", "mspe", "nsel")

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
    
    #print(mspe)
    ## ALTERNATIVE USING CV.SPLS
    #cv_out <- cv.spls(x, y, fold = n, K = ncomp, eta = eta, scale.x = F, scale.y = F, plot.it = F)
    #mspe <- cv_out$mspemat[1]
    #print(mspe)
    
    # print(paste("ncomp:", ncomp, "eta:", eta))
    # print(paste("adj_r2:", adj_r2))
    # print(paste("mspe:", mspe))
    
    newRow <- c(ncomp, eta, r2, adj_r2 , mspe, nsel)
    ncomp_eta_mn <- rbind(ncomp_eta_mn, setNames(newRow,names(ncomp_eta_mn)))
    
    print(c(ncomp, eta))
  }
}

#save(ncomp_eta_mn, file = paste('rslt/comp/ncomp_eta_gs_mn_NEW_', lubridate::today(), '.RData', sep = ''))

# get rid of row of 0s
ncomp_eta_mn <-  ncomp_eta_mn|> 
  slice(-1)

ncomp_eta_mn |> 
  filter(adj_r2 == max(adj_r2) | mspe == min(mspe))


ncomp_eta_mn |> 
  filter(adj_r2 > 0.99 & mspe < 1)

# highest adj_r2: ncomp = 7, eta = 0.758
# lowest mspe: ncomp = 8, eta = 0.68

# 3D graph
# adjusted r2
with(ncomp_eta_mn , plot3d(x = ncomp, y = eta, z = adj_r2))
# mspe
with(ncomp_eta_mn , plot3d(x = ncomp, y = eta, z = mspe))

# boxplots of ncomp vs adjusted r^2
# adjusted r2 starts to decline at around 10 => model gets less helpful after 8 components
ncomp_eta_mn |>
  ggplot(aes(x=ncomp, y = adj_r2)) + geom_boxplot(aes(group= ncomp)) +
  scale_color_viridis() + theme_bw() + labs(title="Humpback Whales - Adj R^2 vs ncomp")


ncomp_eta_mn |>
  ggplot(aes(x=ncomp, y = mspe)) + geom_boxplot(aes(group=ncomp)) +
  scale_color_viridis() + theme_bw() + labs(title="Humpback Whales - MSPE vs ncomp")

ncomp_eta_mn |>
  ggplot(aes(x=ncomp, y = nsel)) + geom_boxplot(aes(group=ncomp)) +
  scale_color_viridis() + theme_bw() + labs(title="Humpback Whales - # ASVs selected vs ncomp")


# scatterplot of eta vs adjusted r^2
# ncomp around 12-16 for optimal adjusted r2
ncomp_eta_mn |>
  ggplot(aes(x=eta, y = adj_r2, color= factor(ncomp))) + geom_point() + 
  scale_color_viridis(discrete = TRUE) 
# look into lower number of components (try to fit with lower number)
# advantages to go to 3-4 instead of 2?

# lower etas for better r2
ncomp_eta_mn |>
  ggplot(aes(x=ncomp, y = adj_r2, color= eta)) + geom_point() + 
  scale_color_viridis(discrete = FALSE) 

# adj r2 line plot
# looks like eta in "green" range best (0.25-0.35)
# looks like ncomp around 8 is best
ggplot(ncomp_eta_mn, aes(x = ncomp, y = adj_r2, color = factor(eta))) +
  geom_line() +
  theme_bw() + labs(title="Humpback Adj R^2 Lineplot\n")

# mspe line plot
# looks like ncomp around 8 is best
# looks like in the blue/purple range is best (0.60-0.68)
ggplot(ncomp_eta_mn, aes(x = ncomp, y = mspe, color = factor(eta))) +
  geom_line() +
  theme_bw() + labs(title="Humpback MSPE Lineplot\n")

