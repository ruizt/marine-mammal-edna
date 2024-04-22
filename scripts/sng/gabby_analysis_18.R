# RUN ONCE, THEN COMMENT OUT
source('rslt/edna-preprocessing-18s.R')
rm(list = ls())

# script starts here
library(tidyverse)
library(vegan)
load('data/ncog-18s.RData')
library(janitor)
library(viridis)

## metadata summaries

# number of samples per cruise (27 cruises)

metadata18 <- metadata
metadata18 %>% count(Cruise, name = 'n.samples')

# number of samples by transect and cruise (15 transects)
metadata18 %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') %>%
  count(Cruise, line) %>%
  arrange(line, Cruise) %>%
  pivot_wider(names_from = line, values_from = n) %>% 
  print(n = 27)

# number of samples by station, transect, and cruise (91 stations)
metadata18 %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') %>%
  count(Cruise, line, sta)

# number of samples by depth and cruise
metadata18 %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') %>%
  mutate(depth.rng = cut_interval(Depthm, n = 2)) %>%
  count(Cruise, depth.rng) %>%
  pivot_wider(names_from = depth.rng, values_from = n)

# common transects sampled across cruises (filter to these for common area) 
sampling_counts18 <- metadata18 %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') %>%
  count(Cruise, line) %>%
  mutate(line = as.numeric(line)) %>%
  filter((76 < line) & (line < 93.4)) 

# transect IDs
transects18 <- sampling_counts18 %>% pull(line) %>% unique()

# number of observations per transect/cruise for common transects
sampling_counts18 %>%
  arrange(line, Cruise) %>%
  pivot_wider(names_from = line, values_from = n) %>% 
  print(n = 27)

## 18s edna summaries

cruises <- edna_samples %>% pull(cruise) %>% unique()

# filter to common transects 
edna18s_reads <- edna_samples %>%
  filter(as.numeric(line) %in% transects18)

# split into sample info and asvs
asv <- edna18s_reads %>%
  select(starts_with('asv'))
sampinfo <- edna18s_reads %>%
  select(-starts_with('asv'))


# asv counts by sample
sampinfo %>%
  bind_cols(num.asvs = specnumber(asv))

# alpha diversity measures
alpha_div18 <- sampinfo %>%
  bind_cols(alpha.div.sh = diversity(asv, index = 'shannon'),
            alpha.div.si = diversity(asv, index = 'simpson'),
            alpha.div.fi = fisher.alpha(asv))

# diversity by depth
alpha_div18 %>%
  mutate(depth = as.numeric(depthm),
         depth.fac = cut_number(depth, 3)) %>%
  ggplot(aes(x = alpha.div.sh, y = depth.fac)) + 
  geom_boxplot()

# diversity distribution by cruise
alpha_div18 %>%
  ggplot(aes(y = alpha.div.sh, x = cruise)) +
  geom_boxplot(coef = 2, outlier.shape = NA) +
  scale_y_continuous(limits = c(3, 6)) +
  theme(axis.text.x = element_text(angle = 90))

# diversity distribution by cruise
alpha_div18 %>%
  group_by(cruise) %>%
  summarize(alpha = mean(alpha.div.sh)) %>%
  ggplot(aes(x = cruise, y = alpha, group = 1)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90))

# from the whole dataset:
metadata18b<- metadata18 |>
  mutate(Sample.Name = paste("X", Sample.Name, sep = "")) |>
  rename("sample.id" = "Sample.Name") |>
  mutate(Lat_Dec = as.numeric(Lat_Dec)) 

# join alpha div 18S with metadata 18S
alpha_meta_18 <- alpha_div18|>
  inner_join(metadata18b,  by = "sample.id") |>
  clean_names()|>
  mutate(depthm = as.numeric(depthm)) |>
  rename("alpha.div.sh" = "alpha_div_sh")

alpha_meta_18 <- alpha_meta_18 |>
  mutate(depth_bins = cut_number(depthm, n = 3),
         depth_bins2 = cut_interval(depthm, n = 10))|>
  filter(!is.na(depth_bins)) |>
  filter(!is.na(depth_bins2))

# temp and alpha binned by depth
alpha_meta_18|>
  ggplot(aes(x=t_deg_c , y = alpha.div.sh, col="orangered", alpha=0.5)) +
  facet_wrap(~depth_bins) +
  geom_point()+
  guides(color=FALSE, alpha=FALSE)+
  labs(title="Temperature vs Alpha Div binned by Depth")

# color by depth
alpha_meta_18|>
  ggplot(aes(x=t_deg_c , y = alpha.div.sh, col=depthm, alpha=0.5)) +
  geom_point()+
  guides(alpha=FALSE)+
  scale_color_continuous(name = "Depth", trans = "reverse") +
  guides(alpha = FALSE)

# temp vs depth
alpha_meta_18|>
  group_by(cruise) |>
  mutate(mean_temp = mean(t_deg_c),
         mean_depth = mean(depthm))|>
  ggplot(aes(x=mean_depth, y = mean_temp, col=alpha.div.sh, alpha=0.5))+
  geom_point()+
  guides(alpha=FALSE)+
  #scale_color_continuous(name = "Alpha Div") +
  scale_color_viridis(name="Alpha Div") +
  guides(alpha = FALSE)

# salinity vs alpha div binned by depth
alpha_meta_18 |>
  ggplot(aes(x=salnty , y = alpha.div.sh, col="orangered", alpha=0.5)) +
  facet_wrap(~depth_bins) +
  geom_point()+
  guides(color=FALSE, alpha=FALSE)+
  labs(title="Salinity vs Alpha Div binned by Depth")
# salinity vs alpha div colored by depth
alpha_meta_18|>
  ggplot(aes(x=salnty , y = alpha.div.sh, col=depthm, alpha=0.5)) +
  geom_point()+
  guides(alpha=FALSE)+
  scale_color_continuous(name = "Depth", trans = "reverse") +
  guides(alpha = FALSE)
# salinity vs depth
alpha_meta_18|>
  group_by(cruise) |>
  mutate(mean_sal = mean(salnty),
         mean_depth = mean(depthm))|>
  ggplot(aes(x=mean_depth, y = mean_sal, col=alpha.div.sh, alpha=0.5))+
  geom_point()+
  guides(alpha=FALSE)+
  scale_color_viridis(name="Alpha Div") 


# oxygen vs alpha div binned by depth
alpha_meta_18 |>
  ggplot(aes(x=o2ml_l , y = alpha.div.sh, col="orangered", alpha=0.5)) +
  facet_wrap(~depth_bins) +
  geom_point()+
  guides(color=FALSE, alpha=FALSE)+
  labs(title="Oxygen vs Alpha Div binned by Depth")
# oxygen vs alpha div colored by depth
alpha_meta_18|>
  ggplot(aes(x=o2ml_l , y = alpha.div.sh, col=depthm, alpha=0.5)) +
  geom_point()+
  guides(alpha=FALSE)+
  scale_color_continuous(name = "Depth", trans = "reverse") +
  guides(alpha = FALSE)
# oxygen vs depth
alpha_meta_18|>
  group_by(cruise) |>
  mutate(mean_oxy = mean(o2ml_l),
         mean_depth = mean(depthm))|>
  ggplot(aes(x=mean_depth, y = mean_oxy, col=alpha.div.sh, alpha=0.5))+
  geom_point()+
  guides(alpha=FALSE)+
  scale_color_viridis(name="Alpha Div") 



# phosphate vs alpha div binned by depth
alpha_meta_18 |>
  ggplot(aes(x=po4ug , y = alpha.div.sh, col="orangered", alpha=0.5)) +
  facet_wrap(~depth_bins) +
  geom_point()+
  guides(color=FALSE, alpha=FALSE)+
  labs(title="Phosphate vs Alpha Div binned by Depth")
# phosphate vs alpha div colored by depth
alpha_meta_18|>
  ggplot(aes(x=po4ug , y = alpha.div.sh, col=depthm, alpha=0.5)) +
  geom_point()+
  guides(alpha=FALSE)+
  scale_color_continuous(name = "Depth", trans = "reverse") +
  guides(alpha = FALSE)
# phosphate vs depth
alpha_meta_18|>
  group_by(cruise) |>
  mutate(mean_po4 = mean(po4ug),
         mean_depth = mean(depthm))|>
  ggplot(aes(x=mean_depth, y = mean_po4, col=alpha.div.sh, alpha=0.5))+
  geom_point()+
  guides(alpha=FALSE)+
  scale_color_viridis(name="Alpha Div") 


# silicate vs alpha div binned by depth
alpha_meta_18 |>
  ggplot(aes(x=si_o3ug , y = alpha.div.sh, col="orangered", alpha=0.5)) +
  facet_wrap(~depth_bins) +
  geom_point()+
  guides(color=FALSE, alpha=FALSE)+
  labs(title="Silicate vs Alpha Div binned by Depth")
# silicate vs alpha div colored by depth
alpha_meta_18|>
  ggplot(aes(x=si_o3ug , y = alpha.div.sh, col=depthm, alpha=0.5)) +
  geom_point()+
  guides(alpha=FALSE)+
  scale_color_continuous(name = "Depth", trans = "reverse") +
  guides(alpha = FALSE)
# silicate vs depth
alpha_meta_18|>
  group_by(cruise) |>
  mutate(mean_sio3 = mean(si_o3ug),
         mean_depth = mean(depthm))|>
  ggplot(aes(x=mean_depth, y = mean_sio3, col=alpha.div.sh, alpha=0.5))+
  geom_point()+
  guides(alpha=FALSE)+
  scale_color_viridis(name="Alpha Div") 


# nitrate vs alpha div binned by depth
alpha_meta_18 |>
  ggplot(aes(x=no3ug , y = alpha.div.sh, col="orangered", alpha=0.5)) +
  facet_wrap(~depth_bins) +
  geom_point()+
  guides(color=FALSE, alpha=FALSE)+
  labs(title="Nitrate vs Alpha Div binned by Depth")
# nitrate vs alpha div colored by depth
alpha_meta_18|>
  ggplot(aes(x=no3ug , y = alpha.div.sh, col=depthm, alpha=0.5)) +
  geom_point()+
  guides(alpha=FALSE)+
  scale_color_continuous(name = "Depth", trans = "reverse") +
  guides(alpha = FALSE)
# nitrate vs depth
alpha_meta_18|>
  group_by(cruise) |>
  mutate(mean_no3 = mean(no3ug),
         mean_depth = mean(depthm))|>
  ggplot(aes(x=mean_depth, y = mean_no3, col=alpha.div.sh, alpha=0.5))+
  geom_point()+
  guides(alpha=FALSE)+
  scale_color_viridis(name="Alpha Div") 


# ammonia vs alpha div binned by depth
alpha_meta_18 |>
  ggplot(aes(x=nh3ug , y = alpha.div.sh, col="orangered", alpha=0.5)) +
  facet_wrap(~depth_bins) +
  geom_point()+
  guides(color=FALSE, alpha=FALSE)+
  labs(title="Ammonia vs Alpha Div binned by Depth")
# ammonia vs alpha div colored by depth
alpha_meta_18|>
  ggplot(aes(x=nh3ug , y = alpha.div.sh, col=depthm, alpha=0.5)) +
  geom_point()+
  guides(alpha=FALSE)+
  scale_color_continuous(name = "Depth", trans = "reverse") +
  guides(alpha = FALSE)
# ammonia vs depth
alpha_meta_18|>
  group_by(cruise) |>
  mutate(mean_nh3 = mean(nh3ug),
         mean_depth = mean(depthm))|>
  ggplot(aes(x=mean_depth, y = mean_nh3, col=alpha.div.sh, alpha=0.5))+
  geom_point()+
  guides(alpha=FALSE)+
  scale_color_viridis(name="Alpha Div") 

# chlorophyll vs alpha div binned by depth
alpha_meta_18|>
  ggplot(aes(x=chlor_a , y = alpha.div.sh, col="orangered", alpha=0.5)) +
  facet_wrap(~depth_bins) +
  geom_point()+
  guides(color=FALSE, alpha=FALSE)+
  labs(title="Chlorophyll vs Alpha Div binned by Depth")
# chlorophyll vs alpha div colored by depth
alpha_meta_18|>
  ggplot(aes(x=chlor_a , y = alpha.div.sh, col=depthm, alpha=0.5)) +
  geom_point()+
  guides(alpha=FALSE)+
  scale_color_continuous(name = "Depth", trans = "reverse") +
  guides(alpha = FALSE)
# chlorophyll vs depth
alpha_meta_18|>
  group_by(cruise) |>
  mutate(mean_chlor = mean(chlor_a),
         mean_depth = mean(depthm))|>
  ggplot(aes(x=mean_depth, y = mean_chlor, col=alpha.div.sh, alpha=0.5))+
  geom_point()+
  guides(alpha=FALSE)+
  scale_color_viridis(name="Alpha Div") 


# 9 taxa
taxa9 = read_csv("rslt/taxa-weights-18s-spls.csv")
head(taxa9)
taxa9$short.id

# get sample IDs of 9 taxa
reads_9 <- edna18s_reads|>
  select(sample.id, cruise, line, sta, asv.1638, asv.12465, asv.13584, asv.17696, 
         asv.22999, asv.28914, asv.31220, asv.32287, asv.37041) |>
  rename(sample_id = sample.id)

# merge with alpha meta data
alpha_meta_9 <- reads_9 |>
  inner_join(alpha_meta_18,  by = "sample_id") |>
  rename(cruise=cruise.x, line = line.x, sta=sta.x) |>
  select(-c(cruise.y, line.y, sta.y, cruise_2))
  

# Temperature vs asv.1638 binned by Depth
alpha_meta_9|>
  ggplot(aes(x=t_deg_c , y = asv.1638, col="orangered", alpha=0.5)) +
  facet_wrap(~depth_bins) +
  geom_point()+
  guides(color=FALSE, alpha=FALSE)+
  labs(title="Temperature vs asv.1638 binned by Depth")


# asv vs depth
asv_columns <- c("asv.1638", "asv.12465", "asv.13584", "asv.17696", 
                 "asv.22999", "asv.28914", "asv.31220", "asv.32287", "asv.37041")

sums <- list()
plots <- list()

for (asv_col in asv_columns) {
  # sums
  sum_val <- sum(alpha_meta_9[[asv_col]])
  sums[[asv_col]] <- sum_val
  
  # scatterplot
  plot_point <- alpha_meta_9 |>
    ggplot(aes(x=depthm , y = .data[[asv_col]], col="orangered", alpha=0.5)) +
    geom_point()+
    guides(color=FALSE, alpha=FALSE)+
    labs(title=paste("Depth vs", asv_col))
  
  # boxplot 
  # note: boxplots are cut_interval, not cut_number b/c R wouldn't let
  # me do cut_number
  plot_boxplot <- alpha_meta_9 |>
    ggplot(aes(x=depth_bins2 , y = .data[[asv_col]], col="orangered")) +
    geom_boxplot()+
    guides(color=FALSE, alpha=FALSE)+
    labs(title=paste("Depth vs", asv_col))
  
  # store plots
  plots[[paste("point", asv_col)]] <- plot_point
  plots[[paste("boxplot", asv_col)]] <- plot_boxplot
}

# asv.1638
plots[[paste("point", "asv.1638")]]
plots[[paste("boxplot", "asv.1638")]]

# asv.12465
plots[[paste("point", "asv.12465")]]
plots[[paste("boxplot", "asv.12465")]]

# asv.13584
plots[[paste("point", "asv.13584")]]
plots[[paste("boxplot", "asv.13584")]]

# asv.17696
plots[[paste("point", "asv.17696")]]
plots[[paste("boxplot", "asv.17696")]]

# asv.22999
plots[[paste("point", "asv.22999")]]
plots[[paste("boxplot", "asv.22999")]]

# asv.28914
plots[[paste("point", "asv.28914")]]
plots[[paste("boxplot", "asv.28914")]]

# asv.31220
plots[[paste("point", "asv.31220")]]
plots[[paste("boxplot", "asv.31220")]]

# asv.32287
plots[[paste("point", "asv.32287")]]
plots[[paste("boxplot", "asv.32287")]]

# asv.37041
plots[[paste("point", "asv.37041")]]
plots[[paste("boxplot", "asv.37041")]]

# df with the asv sums
asv_sums <- data.frame(
  asv = c("asv.1638", "asv.12465", "asv.13584", "asv.17696", "asv.22999", 
          "asv.28914", "asv.31220", "asv.32287", "asv.37041"), 
  asv_sum = c(sums[["asv.1638"]], sums[["asv.12465"]], sums[["asv.13584"]],
              sums[["asv.17696"]], sums[["asv.22999"]], sums[["asv.28914"]],
              sums[["asv.31220"]], sums[["asv.32287"]], sums[["asv.37041"]])
)

# divide all read counts by sum (find in preprocessing) -> relative abundances
# depth profile: density chart
# create more bins (10) for boxplots
# look at asvs vs oxygen and chlorophyll

plots_oxy <- list()
# Plot asvs vs. Oxygen binned by depth
for (asv_col in asv_columns) {
  oxy_plot <- alpha_meta_9 |>
    ggplot(aes(x = o2ml_l, y = .data[[asv_col]], col = "orangered", alpha = 0.5)) +
    facet_wrap(~depth_bins) +
    geom_point() +
    guides(color = FALSE, alpha = FALSE) +
    labs(title = paste("Oxygen vs", asv_col, "binned by Depth"))
  
  plots_oxy[[paste(asv_col)]] <- oxy_plot
}


# asv.1638
plots_oxy[[paste("asv.1638")]]

# asv.12465
plots_oxy[[paste("asv.12465")]]

# asv.13584
plots_oxy[[paste("asv.13584")]]

# asv.17696
plots_oxy[[paste("asv.17696")]]

# asv.22999
plots_oxy[[paste("asv.22999")]]

# asv.28914
plots_oxy[[paste("asv.28914")]]

# asv.31220
plots_oxy[[paste("asv.31220")]]

# asv.32287
plots_oxy[[paste("asv.32287")]]

# asv.37041
plots_oxy[[paste("asv.37041")]]



plots_chlor <- list()
# Plot asvs vs. chlorophyll binned by depth
for (asv_col in asv_columns) {
  chlor_plot <- alpha_meta_9 |>
    ggplot(aes(x = log(chlor_a), y = log(.data[[asv_col]]+1), col = "orangered", alpha = 0.5)) +
    facet_wrap(~depth_bins) +
    geom_point() +
    guides(color = FALSE, alpha = FALSE) +
    labs(title = paste("Chlorophyll vs", asv_col, "binned by Depth"))
  
  plots_chlor[[paste(asv_col)]] <- chlor_plot
}
# do log transformation
# asv.1638
plots_chlor[[paste("asv.1638")]]

# asv.12465 (abundant in almost every sample)
plots_chlor[[paste("asv.12465")]]

# asv.13584
plots_chlor[[paste("asv.13584")]]

# asv.17696
plots_chlor[[paste("asv.17696")]]

# asv.22999
plots_chlor[[paste("asv.22999")]]

# asv.28914
plots_chlor[[paste("asv.28914")]]

# asv.31220
plots_chlor[[paste("asv.31220")]]

# asv.32287
plots_chlor[[paste("asv.32287")]]

# asv.37041
plots_chlor[[paste("asv.37041")]]