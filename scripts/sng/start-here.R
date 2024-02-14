# RUN ONCE, THEN COMMENT OUT
#source('scripts/sng/edna-preprocessing.R')
 #rm(list = ls())

# script starts here
library(tidyverse)
library(vegan)
load('data/ncog-16s.RData')

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

#Alpha Diversity By Depth

#Boxplot
alpha_div |> 
   mutate( depth.rng = cut_number(depthm, n = 3)) |> 
  group_by(depth.rng) |>  
  ggplot(aes(x = depth.rng, y = alpha.div.sh, group = depth.rng)) +
  geom_boxplot() +
  theme_bw() + ylab("Alpha Diversity") + xlab("Depth (M)") + title("Alpha Diversity by Depth")

alpha_div |> 
  mutate( depth.rng = cut_number(depthm, n = 3)) |> 
  count(depth.rng)
#Scatter plot

alpha_div$depthm = as.numeric(alpha_div$depthm)


div.fit1 = lm(alpha.div.sh~depthm, data = alpha_div)
summary(div.fit1)

alpha_div |> 
  ggplot(aes(x= depthm, y = alpha.div.sh)) + 
  geom_point() + geom_smooth(method = "lm", se = FALSE)

# GROUP BY SEASON/QUARTER
#Log transform on X
div.fit2 = lm(alpha.div.sh~log(depthm+1), data = alpha_div)
summary(div.fit2)
alpha_div  |> 
  ggplot(aes(x= log(depthm+1), y = alpha.div.sh)) + 
  geom_point() + geom_smooth(method = "lm", se = FALSE)

#Log transform on X and Y^2 transformation
div.fit3 = lm(alpha.div.sh^2~log(depthm+1), data = alpha_div)
summary(div.fit3)

alpha_div  |>
  ggplot(aes(x= log(depthm+1), y = alpha.div.sh^2)) + 
  geom_point() + geom_smooth(method = "lm", se = FALSE)

# Cruise and depth
alpha_div |> 
  group_by(cruise) |>  
  ggplot(aes(x = cruise, y = depthm )) +
  geom_boxplot()  +
  theme_bw() + ylab("Depth") + xlab("Cruise") + title("Depth by Cruise")

# Deep alpha divs

alpha_div |> 
  filter(depthm >200) |> 
  #select(alpha.div.sh) |> 
  summarise(avgAlpha = mean(alpha.div.sh))

alpha_div |> 
  filter(depthm <25) |> 
  summarise(avgAlpha = mean(alpha.div.sh))

alpha_div |> 
  filter(depthm >25, depthm < 200) |> 
  summarise(avgAlpha = mean(alpha.div.sh))

alpha_div |> 
  filter(depthm >38, depthm <100) |> 
  summarise(avgAlpha = mean(alpha.div.sh))

alpha_div |> 
  filter(depthm >100, depthm <200) |> 
  summarise(avgAlpha = mean(alpha.div.sh))

# Average alpha div by depth

alpha_div |> 
  filter(depthm < 500) |> 
  group_by(depthm) |> 
  summarize(avgDiv = mean(alpha.div.sh)) |> 
  ggplot(aes(x=depthm, y = avgDiv)) + geom_point() + geom_smooth(method = "lm", se = FALSE)


alpha_div |> 
  filter(depthm < 500, depthm >10) |> 
  group_by(depthm) |> 
  summarize(avgDiv = mean(alpha.div.sh)) |> 
  ggplot(aes(x=log(depthm,2), y = avgDiv)) + geom_point() + geom_smooth(method = "lm", se = FALSE)


# By Month
alpha_div |> 
  mutate(month =substring(cruise, 6,8)) |> 
  group_by(month) |> 
  ggplot(aes(x = month, y = alpha.div.sh)) + geom_boxplot()

# By Year 

# Alpha diversity
alpha_div |> 
  mutate(year = substring(cruise, 2,5)) |> 
  group_by(year) |> 
  ggplot(aes(x = year, y = alpha.div.sh)) + geom_boxplot()

# Num samples
alpha_div |> 
  mutate(year = substring(cruise, 2,5)) |> 
  count(year)

#Average Depth By Year
alpha_div |> 
  mutate(year = substring(cruise, 2,5)) |> 
  group_by(year) |> 
  summarize(avgDepth = mean(depthm), avgAlphaDiv = mean(alpha.div.sh))


# pairwise beta diversity for one cruise
edna16s_reads %>%
  filter(cruise == cruises[1], as.numeric(depthm) < 30) %>%
  betadiver() %>%
  plot()

# read docs
?betadiver

# indices supported
betadiver(help = T)

# further explorations?
# Group by cruise, transect, station, depth
edna16s_reads %>%
  mutate(line = as.numeric(line)) |> 
  filter(line == transects[1]) %>%
  betadiver() %>%
  plot()


# CHECK OUT ALPHA DIV BY DEPTH AND SEASON
# CHECK OUT BETA DIV BY DEPTH AND SEASON

alpha_div <- alpha_div |> 
  mutate(month = as.numeric(substring(cruise, 6,8)))
  
#Create season variable
alpha_div_season <- 
  alpha_div |> 
  filter(depthm < 500) |> 
  mutate(season = case_when(
                    month %in% c(12,1,2) ~ "Winter",
                    month %in% c(3,4,5) ~ "Spring",
                    month %in% c(6,7,8) ~ "Summer",
                    month %in% c(9,10,11) ~ "Fall"))
alpha_div_season <- 
  alpha_div_season |> 
  mutate(season = as_factor(season),
         depth.rng = cut_number(depthm, n = 3))

# Sample sizes
alpha_div_season |> 
  count(season)

# Boxplot by season
alpha_div_season |> 
  ggplot(aes(x = season, y = alpha.div.sh)) + 
  geom_boxplot()

# scatterplot by season
alpha_div_season |> 
  ggplot(aes(x = depthm, y = alpha.div.sh)) + 
  geom_point() + 
  facet_wrap(~season)

# scatterplot with averages by season 
alpha_div_season |> 
  group_by(season, depthm) |> 
  summarize(avgDiv = mean(alpha.div.sh)) |> 
  ggplot(aes(x = depthm, y = avgDiv)) + 
  geom_point() + 
  facet_wrap(~season)

library(kableExtra)
library(knitr)
# Means by season

alpha_div_season |> 
  group_by(season) |>
  summarize(N = length(alpha.div.sh),
            avgDiv = round(mean(alpha.div.sh),2),
            sdDiv = round(sd(alpha.div.sh),2)
  ) |> 
  kbl(col.names = c("Season","Sample Size", "Mean Diversity", "SD"),
      caption = "Average Alpha Diversity By Season and Depth Range") |> 
  kable_classic_2(lightable_options = "striped") 

# Means by season and depth
alpha_div_season |> 
  group_by(depth.rng, season) |>
  summarize(N = length(alpha.div.sh),
            avgDiv = round(mean(alpha.div.sh),2),
            sdDiv = round(sd(alpha.div.sh),2)
            ) |> 
  kbl(col.names = c("Depth Range", "Season","Sample Size", "Mean Diversity", "SD"),
      caption = "Average Alpha Diversity By Season and Depth Range") |> 
  kable_classic_2(lightable_options = "striped") 

edna16s_reads |> 
  filter(as.numeric(depthm) < 500) |> 
  count()
# Beta Diversity 
edna16s_reads_season <- edna16s_reads |> 
  filter(as.numeric() < 500) |> 
  mutate(month = as.character(substring(cruise, 6,8)),
         depthm = as.character(depthm))

edna16s_reads_season <- edna16s_reads_season  |> 
  mutate(season = case_when(
    month %in% c("12","1","2") ~ "Winter",
    month %in% c("3","4","5") ~ "Spring",
    month %in% c("6","7","8") ~ "Summer",
    month %in% c("9","10","11") ~ "Fall"))

edna16s_reads_season |> 
  count(month)

edna16s_reads_season%>%
  filter(as.character(season) == "Winter",
         as.numeric(depthm) <= 10,
         as.numeric(depthm) >= 0) %>%
  betadiver() %>%
  plot()
  








