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
