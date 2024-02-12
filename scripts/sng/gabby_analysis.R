# RUN ONCE, THEN COMMENT OUT
#source('scripts/sng/edna-preprocessing.R')
#rm(list = ls())

# script starts here
library(readr)
library(tidyverse)
library(vegan)
load('data/ncog-16s.RData')
library(janitor)
library(viridis)

calcofi = read_csv("CalCOFIStationOrder.csv")

## metadata summaries

# number of samples per cruise (27 cruises)
metadata %>% count(Cruise, name = 'n.samples')

# number of samples by transect and cruise (15 transects)
metadata %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') %>%
  count(Cruise, line) %>%
  arrange(line, Cruise) %>%
  pivot_wider(names_from = line, values_from = n) %>% 
  print(n = 27)

# number of samples by station, transect, and cruise (91 stations)
metadata %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') %>%
  count(Cruise, line, sta)

# number of samples by depth and cruise
metadata %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') %>%
  mutate(depth.rng = cut_interval(Depthm, n = 2)) %>%
  count(Cruise, depth.rng) %>%
  pivot_wider(names_from = depth.rng, values_from = n)

# common transects sampled across cruises (filter to these for common area) 
sampling_counts <- metadata %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') %>%
  count(Cruise, line) %>%
  mutate(line = as.numeric(line)) %>%
  filter((76 < line) & (line < 93.4)) 

# transect IDs
transects <- sampling_counts %>% pull(line) %>% unique()

# number of observations per transect/cruise for common transects
sampling_counts %>%
  arrange(line, Cruise) %>%
  pivot_wider(names_from = line, values_from = n) %>% 
  print(n = 27)

## 16s edna summaries

cruises <- edna_samples %>% pull(cruise) %>% unique()

# filter to common transects 
edna16s_reads <- edna_samples %>%
  filter(as.numeric(line) %in% transects)

# split into sample info and asvs
asv <- edna16s_reads %>%
  select(starts_with('asv'))
sampinfo <- edna16s_reads %>%
  select(-starts_with('asv'))


# asv counts by sample
sampinfo %>%
  bind_cols(num.asvs = specnumber(asv))

# alpha diversity measures
alpha_div <- sampinfo %>%
  bind_cols(alpha.div.sh = diversity(asv, index = 'shannon'),
            alpha.div.si = diversity(asv, index = 'simpson'),
            alpha.div.fi = fisher.alpha(asv))

# diversity by depth
alpha_div %>%
  mutate(depth = as.numeric(depthm),
         depth.fac = cut_number(depth, 3)) %>%
  ggplot(aes(x = alpha.div.sh, y = depth.fac)) + 
  geom_boxplot()

# diversity distribution by cruise
alpha_div %>%
  ggplot(aes(y = alpha.div.sh, x = cruise)) +
  geom_boxplot(coef = 2, outlier.shape = NA) +
  scale_y_continuous(limits = c(3, 6)) +
  theme(axis.text.x = element_text(angle = 90))

# seasonal fluctuations?
alpha_div %>%
  group_by(cruise) %>%
  summarize(alpha = mean(alpha.div.sh)) %>%
  ggplot(aes(x = cruise, y = alpha, group = 1)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90))

# further explorations??


#alpha div grouping by line boxplot
alpha_div|>
  group_by(line)|>
  ggplot(aes(x=line, y=alpha.div.sh)) + 
  geom_boxplot() + 
  scale_y_continuous(breaks=seq(0,6,1)) +
  theme_bw()

# order by latitude (filter data by depth ranges to see changes from depth-to-depth) 

#alpha div vs line barchart
alpha_div|>
  group_by(line)|>
  summarize(alpha = mean(alpha.div.sh))|>
  ggplot(aes(x=line, y=alpha)) + 
  geom_bar(stat="identity", col="steelblue", fill="steelblue")

# alpha div vs station boxplot
alpha_div |>
  group_by(sta) |>
  ggplot(aes(x=sta, y=alpha.div.sh)) + geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90))

# alpha div vs station barchart
alpha_div|>
  group_by(sta)|>
  summarize(alpha = mean(alpha.div.sh))|>
  ggplot(aes(x=sta, y=alpha)) + 
  geom_bar(stat="identity", col="steelblue", fill="steelblue") +
  theme(axis.text.x = element_text(angle = 90))
# group by line

# alpha div vs depth by line
alpha_div|>
  ggplot(aes(x=depthm, y=alpha.div.sh, color=sta)) + geom_point() +
  facet_wrap(~line)


# pairwise beta diversity for one cruise
edna16s_reads %>%
  filter(cruise == cruises[1], as.numeric(depthm) < 30) %>%
  betadiver() %>%
  plot()

# -------------------------------------------

# join alpha_div data and metadata

metadata2 <- metadata |>
  mutate(Sample.Name = paste("X", Sample.Name, sep = "")) |>
  rename("sample.id" = "Sample.Name") |>
  mutate(Lat_Dec = as.numeric(Lat_Dec)) 
  

loc <- alpha_div|>
  inner_join(metadata2,  by = "sample.id") |>
  clean_names()|>
  mutate(depthm = as.numeric(depthm)) |>
  rename("alpha.div.sh" = "alpha_div_sh")

# alpha div vs station boxplot, ordered by latitude
loc |>
  mutate(lat_dec = as.numeric(lat_dec)) |>
  arrange(lat_dec) |>
  mutate(sta = factor(sta, levels = unique(sta))) |>
  ggplot(aes(x = sta, y = alpha.div.sh)) +  
  geom_boxplot() +
  xlab("Station (ordered by latitude)") + ylab("Alpha Diversity") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90)) 

# alpha div vs line boxplot, ordered by latitude
loc |>
  arrange(lat_dec) |>
  mutate(line = factor(line, levels = unique(line))) |>
  ggplot(aes(x = line, y = alpha.div.sh)) +  
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  xlab("Transect (ordered by latitude)") + ylab("Alpha Diversity") +
  theme_minimal()


# alpha diversity and depth (binned) (boxplots)
loc |>
  mutate(depth_bins = cut(depthm, breaks = seq(0, 150, by = 50)))|>
  filter(!is.na(depth_bins)) |>
  ggplot(aes(x = lat_dec, y = alpha.div.sh, col=line)) +
  geom_boxplot() +
  facet_wrap(~depth_bins)+
  labs(title = "Alpha Diversity by Latitude and Depth") +
  xlab("Latitude") + ylab("Alpha Diversity") +
  scale_color_discrete(name = "Transect") +
  theme_minimal()

# alpha diversity and depth (binned) (scatter)
loc |>
  mutate(depth_bins = cut(depthm, breaks = seq(0, 150, by = 50)))|>
  filter(!is.na(depth_bins)) |>
  ggplot(aes(x = lat_dec, y = alpha.div.sh, col=line)) +
  geom_point() +
  facet_wrap(~depth_bins)+
  labs(title = "Alpha Diversity by Latitude and Depth") +
  xlab("Latitude") + ylab("Alpha Diversity") +
  scale_color_discrete(name = "Transect") +
  theme_minimal()

# latitude vs Alpha Diversity (scatter) w/ depth as color map
loc |>
  ggplot(aes(x=lat_dec, y = alpha.div.sh, color=depthm, alpha=1/2)) + geom_point()+
  labs(title = "Latitude vs Depth") +
  xlab("Latitude") + ylab("Alpha Diversity") +
  scale_color_continuous(name = "Depth", trans = "reverse") +
  guides(alpha = FALSE)

# latitude vs alpha div faceted by season

loc2 <- loc |>
  mutate(date = mdy(date),  
         season = case_when(
           month(date) %in% 3:5 ~ "Spring",
           month(date) %in% 6:8 ~ "Summer",
           month(date) %in% 9:11 ~ "Fall",
           month(date) %in% c(12, 1, 2) ~ "Winter"
         )) 

loc2 |>
  filter(!is.na(season)) |>
  ggplot(aes(x=lat_dec, y = alpha.div.sh, color="coral", alpha=0.5)) + 
  geom_point() +
  facet_wrap(~season)+
  guides(alpha = FALSE, color = FALSE)

# read docs
?betadiver

# indices supported
betadiver(help = T)

# further explorations?

