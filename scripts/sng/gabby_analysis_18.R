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
transects18 <- sampling_counts %>% pull(line) %>% unique()

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
  mutate(depth_bins = cut_number(depthm, n = 3))|>
  filter(!is.na(depth_bins)) 

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


