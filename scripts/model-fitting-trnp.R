library(tidyverse)
library(spls)
library(modelr)

# read in selection frequencies from LOOCV
sel_freq <- read_rds('rslt/loocv/2024-07-20/selection-frequencies.rds')

# read in metrics from LOOCV
metrics <- read_rds('rslt/loocv/2024-07-20/metrics.rds')

## SPECIFYING NCOMP/ETA RANGE --------------------------------------------------
# chosen stability threshold (minimum max selection prob.)
pimax <- 0.8

# desired bound for expected number of false positives
EV <- 5

# total number of candidate ASVs
p <- 6832

# threshold for max average number of selected ASVs
qmax <- sqrt((2*pimax - 1)*p*EV) 

## FILTER EXPLORATION - SPECIES = BP ---------------------------------------------

# grid search: ranges of ncomp, eta -> average # of asvs selected
ncomp_start_grid <- 1:9
ncomp_stop_grid <- 10:1
eta_min_grid <- seq(0,0.9, by = 0.1)
eta_max_grid <- seq(1,0.1, by = -0.1)

num_asvs_eta_ncomp_gs <- function(ncomp_start, ncomp_end, eta_min, eta_max){
  avg_asvs <- metrics |> 
    filter(species == "bp",
           ncomp >= ncomp_start,
           ncomp <= ncomp_end,
           eta >= eta_min,
           eta <= eta_max) |> 
    summarise(avg.n.asv = mean(n.asv),
              sd.n.asv = sd(n.asv))
  
  return(avg_asvs)
}

num_asvs_eta_ncomp_gs(3, 8, 0.4, 0.9)

gs_df <- data.frame(0,0,0,0,0,0)
names(gs_df) = c("ncomp.min", "ncomp.max", "eta.min","eta.max", "avg.n.csv", 'sd.n.csv')
for (ns in ncomp_start_grid){
  for (nf in ncomp_stop_grid){
    if (ns > nf){
      next
    }
    for (emin in eta_min_grid){
      for (emax in eta_max_grid){
        if (emin > emax){
          next
        }
        
        avg_asvs = num_asvs_eta_ncomp_gs(ns,nf,emin,emax) |> select(avg.n.asv)
        sd_asvs = num_asvs_eta_ncomp_gs(ns,nf,emin,emax) |> select(sd.n.asv)
        
        newrow = c(ns,nf,emin,emax,avg_asvs,sd_asvs)
        
        print(newrow)
        
        gs_df <- rbind(gs_df, setNames(newrow, names(gs_df)))
        
      }
    }
  }
}

gs_df <- gs_df |> 
  slice(-1)


# Candidate ranges
candidate_ranges <- gs_df |> 
  filter(avg.n.csv <= qmax,
         avg.n.csv > 0)

candidate_ranges <- candidate_ranges |> 
  mutate(ncomp.range = ncomp.max - ncomp.min,
         eta.range = eta.max - eta.min)

candidate_ranges |> 
  filter(ncomp.range == max(ncomp.range) | eta.range == max(eta.range))

# notice effect of range on sd, for example...
candidate_ranges |> 
  filter(ncomp.range == 9) |>
  filter(avg.n.csv > 40, avg.n.csv < 70) |>
  arrange(eta.range)

# just tinkering here ...
candidate_ranges |> 
  filter(avg.n.csv < 90,
         avg.n.csv > 20,
         eta.max <= 0.9, 
         eta.min >= 0.4,
         ncomp.min > 2,
         ncomp.max == 8,
         avg.n.csv > 15,
         sd.n.csv < 50)

# recalculate EV for new range
qmax_new <- candidate_ranges |>
  filter(ncomp.min == 3, ncomp.max == 8, eta.min == 0.5, eta.max == 0.9) |>
  pull(avg.n.csv)

qmax_new^2/(p*(2*pimax - 1))

## STABLE SET -------------------------------------------------------------------
## INITIAL FILTERING: 1 <= ncomp <= 10 ; 0.4 <= eta <= 1
## => eta >= 0.4

# this is the widest range of ncomp and eta values that still obtains an 
# average number of asvs <= qmax (98.711)

# asvs that have a maximum selection prob of at least 0.8 (max n >= 20)

sel_freq |> 
  filter(species == "bp",
         eta >= 0.6,
         eta <= 0.9,
         ncomp >= 8,
         ncomp <= 8,
         n >= 20) |> 
  select(asv) |> 
  distinct()
  
# asvs that have average selection prob of at least 0.8

sel_freq |> 
  filter(species == "bp",
         eta >= 0.4) |> 
  mutate(prop = n/25) |> 
  group_by(asv) |> 
  summarise(avg.prop = mean(prop)) |> 
  filter(avg.prop >= 0.8)


# asvs that have average selection prob of at least 0.5

sel_freq |> 
  filter(species == "bp",
         eta >= 0.4) |> 
  mutate(prop = n/25) |> 
  group_by(asv) |> 
  summarise(avg.prop = mean(prop)) |> 
  filter(avg.prop >= 0.5)


## another approach -- bin and group
pal.grad <- colorRampPalette(c('red', 'blue'))
pal <- pal.grad(13)

# look at avg nsel by eta bin and ncomp
# find ranges where qmin < q < qmax for all ncomp
metrics |>
  select(species, ncomp, eta, n.asv) |>
  mutate(eta.bin = cut_interval(eta, n = 8)) |>
  group_by(species, ncomp, eta.bin) |>
  summarize(avg.n.asv = mean(n.asv),
            sd.n.asv = sd(n.asv)) |>
  ggplot(aes(x = eta.bin, y = avg.n.asv, color = factor(ncomp))) +
  facet_wrap(~species) +
  geom_point() +
  geom_path(aes(group = ncomp)) +
  scale_y_log10() +
  scale_color_manual(values = pal) +
  theme(axis.text.x = element_text(angle = 90)) 

# say, qmin = 10, qmax = 80 groupwise; find q across full range
metrics |>
  select(species, ncomp, eta, n.asv) |>
  filter(ncomp == 8, 
         eta >= 0.65, eta <= 0.85) |>
  group_by(species, ncomp) |>
  summarize(avg.n = mean(n.asv), 
            sd.n = sd(n.asv)) |>
  slice_max(avg.n)

# upper bound on expected no. false positives
qmax <- 57
qmax^2/(p*(2*pimax - 1))

# stable sets
stable_sets <- sel_freq |> 
  filter(eta >= 0.65,
         eta <= 0.88,
         ncomp == 8,
         n >= 20) |> 
  group_by(species) |> 
  distinct(asv) |>
  nest(data = asv) |>
  transmute(stable.set = map(data, unlist))

# tinkering with fitting models -- it works!

load('data/processed/ncog18sv9-2024-07-20.RData')
load('data/processed/mm-sightings-2024-07-20.RData')
whales <- inner_join(sightings, edna, by = 'cruise') 

loocv_fn <- function(.train, .var){
  .asv <- filter(stable_sets, species == .var) %>% pull(stable.set) %$% .[[1]]
  names(.asv) <- NULL
  x.train <- .train |> as.data.frame() |> dplyr::select(all_of(.asv))
  y.train <- .train |> as.data.frame() |> pull({{.var}})
  fit <- spls(x.train, y.train,
              K = 8, eta = 0,
              scale.x = F, scale.y = F)

  return(fit)
}

pred_fn <- function(.test, .fit){
  asv <- rownames(.fit$projection)
  newdata <- .test |> as.data.frame() |> dplyr::select(all_of(asv))
  out <- predict(.fit, newx = newdata, type = 'fit')[1, 1]
  return(out)
}

loocv <- crossv_loo(whales) |>
  expand_grid(species = c('bm', 'bp', 'mn')) |>
  mutate(fit = map2(train, species,
                    ~loocv_fn(.x, .y)),
         pred = map2(test, fit, ~pred_fn(.x, .y)),
         obs = map2(test, species, ~pull(as.data.frame(.x), .y))) |>
  unnest(c(pred, obs))

loocv |>
  group_by(species) |>
  summarize(cor = cor(pred, obs),
            mspe = mean((pred - obs)^2))

loocv |>
  mutate(fitted = map(fit, ~predict(.x, type = 'fit')[,1]),
         train.obs = map2(train, species, ~pull(as.data.frame(.x), .y))) |>
  select(species, fitted, train.obs) |>
  mutate(cor = map2(fitted, train.obs, ~cor(.x, .y))) |>
  unnest(cor) |>
  group_by(species) |>
  summarize(mean(cor))


metrics |>
  group_by(ncomp, eta, species) |>
  summarize(mspe = mean(sq.pred.err),
            cor = cor(pred, pred + pred.err),
            df = mean(n.asv)) |>
  group_by(species) |>
  slice_max(cor)
