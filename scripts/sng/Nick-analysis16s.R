# Nick Patrick
# eDNA Analysis

library(tidyverse)
library(vegan)
load('data/ncog-16s.RData')

# alpha diversity measures
alpha_div <- sampinfo %>%
  bind_cols(alpha.div.sh = diversity(asv, index = 'shannon'),
            alpha.div.si = diversity(asv, index = 'simpson'),
            alpha.div.fi = fisher.alpha(asv))

# diversity by depth
alpha_div %>%
  mutate(depth = as.numeric(depthm),
         depth.fac = cut_number(depth, 3)) %>%
  ggplot(aes(x = depth.fac, y = alpha.div.sh)) + 
  geom_boxplot()


#Alpha Diversity By Depth

#Boxplot
alpha_div |> 
  mutate( depth.rng = cut_number(depthm, n = 3)) |> 
  group_by(depth.rng) |>  
  ggplot(aes(x = depth.rng, y = alpha.div.sh, group = depth.rng)) +
  geom_boxplot() +
  theme_bw() + ylab("Alpha Diversity") + xlab("Depth (M)") + title("Alpha Diversity by Depth")

alpha_div |> 
  mutate( depthm = as.numeric(depthm),
          depth.rng = cut_number(depthm, n = 3)
  ) |> 
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

# Average alpha div by depth
alpha_div |> 
  filter(depthm < 500) |> 
  group_by(depthm) |> 
  summarize(avgDiv = mean(alpha.div.sh)) |> 
  ggplot(aes(x=depthm, y = avgDiv)) + geom_point() + geom_smooth(method = "lm", se = FALSE)


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

## SEASONAL ALPHA DIV
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
         depth.rng = cut_number(as.numeric(depthm), n = 3))

# Sample sizes
alpha_div_season |> 
  count(season)

# Boxplot by season
alpha_div_season |> 
  ggplot(aes(x = season, y = alpha.div.sh)) + 
  geom_boxplot()

# scatterplot by season
alpha_div_season |> 
  mutate(depthm = as.numeric(depthm)) |> 
  ggplot(aes(x = depthm, y = alpha.div.sh)) + 
  geom_point() + 
  facet_wrap(~season) +
  scale_x_continuous(breaks = seq(0,200, by = 20))

# scatterplot with averages by season 
alpha_div_season |> 
  mutate(depthm = as.numeric(depthm)) |> 
  group_by(season, depthm) |> 
  summarize(avgDiv = mean(alpha.div.sh)) |> 
  ggplot(aes(x = depthm, y = avgDiv)) + 
  geom_point() + 
  facet_wrap(~season) +
  scale_x_continuous(breaks = seq(0,200, by = 20))


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

# Boxplot by Season, faceted on depth range
alpha_div_season |> 
  ggplot(aes(x = season, y = alpha.div.sh)) + 
  geom_boxplot() + facet_grid(~depth.rng)

# Boxplot grouped by depth range, faceted on season
alpha_div_season |> 
  ggplot(aes(x = depth.rng, y = alpha.div.sh)) + 
  geom_boxplot() + facet_grid(~season) 

# Beta Diversity 

edna16s_reads_season <- edna16s_reads |> 
  filter(as.numeric(depthm) < 500) |> 
  mutate(month = as.numeric(substring(cruise, 6,8)),
         depthm = as.numeric(depthm))

edna16s_reads_season <- edna16s_reads_season  |> 
  mutate(season = case_when(
    month %in% c(12,1,2) ~ "Winter",
    month %in% c(3,4,5) ~ "Spring",
    month %in% c(6,7,8) ~ "Summer",
    month %in% c(9,10,11) ~ "Fall"))

edna16s_reads_season |> 
  count(season)

edna16s_reads_season%>%
  filter(as.character(season) == "Winter") %>%
  betadiver() %>%
  plot()


edna16s_reads_season%>%
  filter(as.character(season) == "Winter",
         depthm >= 0,
         depthm <= 9) %>%
  betadiver() %>%
  plot()

edna16s_reads_season%>%
  filter(as.character(season) == "Winter",
         depthm ==10) %>%
  betadiver() %>%
  plot()

edna16s_reads_season%>%
  filter(as.character(season) == "Winter",
         depthm > 10,
         depthm <= 20) %>%
  betadiver() %>%
  plot()


##############################################################

#Mammal Exploration


##################################### 

# Create season dataframe
div_mammal_season <- alpha_div %>%
  filter(depthm < 150) %>%
  group_by(cruise) %>%
  summarize(alpha.mean = mean(alpha.div.sh),
            alpha.sd = sd(alpha.div.sh)) %>%
  left_join(mammal, by = 'cruise') %>%
  ungroup()

div_mammal_season <- div_mammal_season |> 
  mutate(month = as.numeric(substring(cruise,6,8)))

div_mammal_season <- div_mammal_season |> 
  mutate(season = case_when(
    month %in% c(12,1,2) ~ "Winter",
    month %in% c(3,4,5) ~ "Spring",
    month %in% c(6,7,8) ~ "Summer",
    month %in% c(9,10,11) ~ "Fall")
  ) 

div_mammal_season[is.na(div_mammal_season)] <- 0

## Mammal sightings by season
div_mammal_season |> 
  group_by(season) |> 
  count()


div_mammal_season |> 
  group_by(season) |> 
  summarise(sum_bp = sum(bp),
            sum_bm = sum(bm),
            sum_mn = sum(mn),
            avg_bp = mean(bp),
            avg_bm = mean(bm), 
            avg_mn = mean(mn))

library(ggbeeswarm)

# BP boxplot by season
div_mammal_season |> 
  group_by(season) |> 
  ggplot(aes(x=season, y = bp)) + 
  geom_boxplot() + 
  geom_beeswarm(color = "blue") +
  ggtitle("BP By Season")


# BM boxplot by season
div_mammal_season |> 
  group_by(season) |> 
  ggplot(aes(x=season, y = bm)) + 
  geom_boxplot() + 
  geom_beeswarm(color = "blue") +
  ggtitle("BM By Season")

#MN boxplot by season
div_mammal_season |> 
  group_by(season) |> 
  ggplot(aes(x=season, y = mn)) + 
  geom_boxplot() + 
  geom_beeswarm(color = "blue") +
  ggtitle("MN By Season")

## Scatter plots faceted by season
div_mammal_season |> 
  ggplot(aes(x=alpha.mean, y = bp)) +
  geom_point() + facet_wrap(~season)

div_mammal_season |> 
  ggplot(aes(x=alpha.mean, y = bm)) +
  geom_point() + facet_wrap(~season)

div_mammal_season |> 
  ggplot(aes(x=alpha.mean, y = mn)) +
  geom_point() + facet_wrap(~season)

# Aggregate metadata

ag_md <- metadata |> 
  mutate(Cruise = paste("X", Cruise, sep = "")) |> 
  group_by(Cruise) |> 
  summarize(avg_Temp = mean(T_degC),
            avg_Salnty = mean(Salnty),
            avg_Chlorophyll = mean(ChlorA),
            avg_Phosphate = mean(PO4ug),
            avg_Silicate = mean(SiO3ug),
            avg_Nitrile = mean(NO3ug),
            avg_Ammonia = mean(NH3ug))

ag_md <- ag_md |> 
  mutate(across(where(is.numeric), ~replace_na(., mean(., na.rm=TRUE))))

# Join metadata with div_mammal_season

div_mammal_meta <- div_mammal_season |> 
  left_join(ag_md, by = join_by(cruise == Cruise))

# Scatterplots BP~watercontents
div_mammal_meta |> 
  ggplot(aes(x=avg_Salnty, y = bp)) +
  geom_point()

div_mammal_meta |> 
  ggplot(aes(x=avg_Temp, y = bp)) +
  geom_point()

div_mammal_meta |> 
  ggplot(aes(x=avg_Nitrile, y = bp)) +
  geom_point()

div_mammal_meta |> 
  ggplot(aes(x=avg_Phosphate, y = bp)) +
  geom_point()

div_mammal_meta |> 
  ggplot(aes(x=avg_Chlorophyll, y = bp)) +
  geom_point()

div_mammal_meta |> 
  ggplot(aes(x=avg_Ammonia, y = bp)) +
  geom_point()

### Which physical variables correlated with depth
### => diversity
### Potential confounding
### non aggregated features
## 16s highly variable across bacteria
## 18s phytoplankton etc

metaSeason <- metadata |> 
  mutate(month = as.numeric(substring( Cruise,5,7)) )

metaSeason <- metaSeason |> 
  mutate(season = as.factor(case_when(
    month %in% c(12,1,2) ~ "Winter",
    month %in% c(3,4,5) ~ "Spring",
    month %in% c(6,7,8) ~ "Summer",
    month %in% c(9,10,11) ~ "Fall")))

# Temp compared by season

#Scatterplot temp by depth and season
metaSeason |> 
  ggplot(aes(x= Depthm, y = T_degC)) +
  geom_point() +
  facet_wrap(~fct_relevel(season, "Winter", "Spring", "Summer", "Fall"))

# Temp By Depth boxplot
metaSeason |> 
  filter(Depthm < 250) |> 
  mutate(depth.rng = cut_number(Depthm,3)) |> 
  ggplot(aes(x =depth.rng, y = T_degC )) +
  geom_boxplot() +
  geom_beeswarm() + labs(y = "Temp (C)", x = "Depth") +
  theme(axis.text = element_text(angle = 90))

# Temp by Depth and Season Boxplot
metaSeason |> 
  filter(Depthm < 250) |> 
  mutate(depth.rng = cut_number(Depthm,3)) |> 
  ggplot(aes(x =fct_relevel(season, "Winter", "Spring", "Summer", "Fall"), y = T_degC )) +
  geom_boxplot() +
  geom_beeswarm() +
  facet_wrap(~depth.rng) + labs(y = "Temp (C)", x = "Season") +
  theme(axis.text = element_text(angle = 90))

### METADATA BY DEPTH/SEASON EXPLORATION

# PO4ug (Phosphate)
# Scatter plot
metaSeason |> 
  ggplot(aes(x=Depthm, y = PO4ug)) +
  geom_point() +
  facet_wrap(~fct_relevel(season, "Winter", "Spring", "Summer", "Fall"))
# Boxplot
metaSeason |> 
  filter(Depthm < 250) |> 
  mutate(depth.rng = cut_number(Depthm,3)) |> 
  ggplot(aes(x=fct_relevel(season, "Winter", "Spring", "Summer", "Fall"), y = PO4ug)) +
  geom_boxplot() +
  facet_wrap(~depth.rng) +
  labs(x= "Season") +
  theme(axis.text = element_text(angle = 90))

# Spring and Summer show highest Phosphate levels across all depth
# No trend across depth

# NO3ug (Nitrate)
# Scatter plot
metaSeason |> 
  ggplot(aes(x=Depthm, y = NO3ug)) +
  geom_point() +
  facet_wrap(~fct_relevel(season, "Winter", "Spring", "Summer", "Fall"))
# Boxplot
metaSeason |> 
  filter(Depthm < 250) |> 
  mutate(depth.rng = cut_number(Depthm,3)) |> 
  ggplot(aes(x=fct_relevel(season, "Winter", "Spring", "Summer", "Fall"), y = NO3ug)) +
  geom_boxplot() +
  facet_wrap(~depth.rng) +
  labs(x= "Season", y = "Nitrate") +
  theme(axis.text = element_text(angle = 90))

# Similar pattern to Phosphate
# Highest in spring and summer
# No pattern across depth

# Si03 (Silicon)
# Scatter plot
metaSeason |> 
  ggplot(aes(x=Depthm, y = SiO3ug)) +
  geom_point() +
  facet_wrap(~fct_relevel(season, "Winter", "Spring", "Summer", "Fall"))

# Boxplot
metaSeason |> 
  filter(Depthm < 250) |> 
  mutate(depth.rng = cut_number(Depthm,3)) |> 
  ggplot(aes(x=fct_relevel(season, "Winter", "Spring", "Summer", "Fall"), y = SiO3ug)) +
  geom_boxplot() +
  facet_wrap(~depth.rng) +
  labs(x= "Season", y = "Silicon") +
  theme(axis.text = element_text(angle = 90))

# Similar pattern to Phosphate and Nitrate
# Highest in spring and summer
# No pattern across depth



# ChlorA (Chlorophyll)
# Scatter plot
metaSeason |> 
  ggplot(aes(x=Depthm, y = ChlorA)) +
  geom_point() +
  facet_wrap(~fct_relevel(season, "Winter", "Spring", "Summer", "Fall"))
# Boxplot
metaSeason |> 
  filter(Depthm < 250) |> 
  mutate(depth.rng = cut_number(Depthm,3)) |> 
  ggplot(aes(x=fct_relevel(season, "Winter", "Spring", "Summer", "Fall"), y = ChlorA)) +
  geom_boxplot() +
  facet_wrap(~depth.rng) +
  labs(x= "Season") +
  theme(axis.text = element_text(angle = 90))

# Spring and summer show highest Chlorophyll levels
# as depth increases, Chlorophyll levels decrease


# Phaeop (Phaeopigment)
# https://en.wiktionary.org/wiki/phaeopigment
#  A non-photosynthetic pigment which is the degradation product 
#  of algal chlorophyll pigments. It is commonly formed during and 
#  after marine phytoplankton blooms.

metaSeason |> 
  ggplot(aes(x=Depthm, y = Phaeop)) +
  geom_point() +
  facet_wrap(~fct_relevel(season, "Winter", "Spring", "Summer", "Fall"))

metaSeason |> 
  filter(Depthm < 250) |> 
  mutate(depth.rng = cut_number(Depthm,3)) |> 
  ggplot(aes(x=fct_relevel(season, "Winter", "Spring", "Summer", "Fall"), y = Phaeop)) +
  geom_boxplot() +
  facet_wrap(~depth.rng) +
  labs(x= "Season", y = "Phaeopigment") +
  theme(axis.text = element_text(angle = 90))

# Similar pattern to Chlorophyll
# Highest in spring and summer
# As depth increases, Phaeopigment decreases


# NH3 (Ammonia)
# Scatter plot
metaSeason |> 
  ggplot(aes(x=Depthm, y = NH3ug)) +
  geom_point() +
  facet_wrap(~fct_relevel(season, "Winter", "Spring", "Summer", "Fall"))
# Boxplot
metaSeason |> 
  filter(Depthm < 250) |> 
  mutate(depth.rng = cut_number(Depthm,3)) |> 
  ggplot(aes(x=fct_relevel(season, "Winter", "Spring", "Summer", "Fall"), y = NH3ug)) +
  geom_boxplot() +
  facet_wrap(~depth.rng) +
  labs(x= "Season", y = "Ammonia") +
  theme(axis.text = element_text(angle = 90))

# Lowest in Fall
# Not much of a pattern anywhere
# Maybe highest in spring/summer??




# ALPHA DIV COMPARED BY METADATA

# Join meta and alpha div
meta2 <- metaSeason %>%
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ')

alpha_meta <- alpha_div |> 
  mutate(cruise = as.double(substring(cruise,2)),
         depthm = as.numeric(depthm)) |> 
  left_join(meta2, by = join_by(cruise == Cruise, line, sta, depthm ==Depthm))


# TEMP 

#  Scatterplot of alpha div by temp, for each season
alpha_meta |> 
  filter(season %in% c("Winter", "Spring", "Summer", "Fall")) |> 
  ggplot(aes(x = T_degC, y = alpha.div.sh)) +
  geom_point() + 
  facet_wrap(~fct_relevel(season, "Winter", "Spring", "Summer", "Fall"))

alpha_meta |> 
  filter(season %in% c("Winter", "Spring", "Summer", "Fall")) |>
  group_by(T_degC, season) |> 
  summarise(avgDiv = mean(alpha.div.sh)) |> 
  ggplot(aes(x = T_degC, y = avgDiv)) +
  geom_point() + 
  facet_wrap(~fct_relevel(season, "Winter", "Spring", "Summer", "Fall"))


# Phosphate
alpha_meta |> 
  filter(season %in% c("Winter", "Spring", "Summer", "Fall")) |> 
  ggplot(aes(x = PO4ug, y = alpha.div.sh)) +
  geom_point() + 
  facet_wrap(~fct_relevel(season, "Winter", "Spring", "Summer", "Fall"))

# Nitrate
alpha_meta |> 
  filter(season %in% c("Winter", "Spring", "Summer", "Fall")) |> 
  ggplot(aes(x = NO3ug, y = alpha.div.sh)) +
  geom_point() + 
  facet_wrap(~fct_relevel(season, "Winter", "Spring", "Summer", "Fall"))

# Silicon

alpha_meta |> 
  filter(season %in% c("Winter", "Spring", "Summer", "Fall")) |> 
  ggplot(aes(x = SiO3ug, y = alpha.div.sh)) +
  geom_point() + 
  facet_wrap(~fct_relevel(season, "Winter", "Spring", "Summer", "Fall"))


# Chlorophyll
alpha_meta |> 
  filter(season %in% c("Winter", "Spring", "Summer", "Fall")) |> 
  ggplot(aes(x = ChlorA, y = alpha.div.sh)) +
  geom_point() + 
  facet_wrap(~fct_relevel(season, "Winter", "Spring", "Summer", "Fall"))

# Phaeopigment
alpha_meta |> 
  filter(season %in% c("Winter", "Spring", "Summer", "Fall")) |> 
  ggplot(aes(x = Phaeop, y = alpha.div.sh)) +
  geom_point() + 
  facet_wrap(~fct_relevel(season, "Winter", "Spring", "Summer", "Fall"))


