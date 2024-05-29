##### 18s



source('~/Documents/marine-mammal-edna/scripts/edna-preprocessing-18s.R')
rm(list = ls())

# script starts here
library(tidyverse)
library(vegan)
load('~/Documents/marine-mammal-edna/data/ncog-18s.RData')
library(janitor)
library(viridis)

## metadata summaries
# number of samples per cruise (27 cruises)

metadata18 <- metadata
metadata18 %>% count(Cruise, name = 'n.samples')

curve(dgamma(x, shape = 10, scale = 15), from = 0, to = 200)

load('~/Documents/marine-mammal-edna/data/ncog-18s.RData')


## lines and temp

alpha_div <- alpha_div |>
  rename(Sample.Name = sample.id)

alpha_div2 <- alpha_div %>%
  select(Sample.Name, line, sta, depthm, alpha.div.sh) %>%
  mutate(Sample.Name = substr(Sample.Name, 2, nchar(alpha_div2$Sample.Name)))


chlor18 <- metadata18 %>%
  select(Sample.Name, Cruise, Depthm, CC_Depth, Time, Lat_Dec, Lon_Dec, Distance, NO3ug, NH3ug, ChlorA, T_degC, O2ml_L)

meta18 <- inner_join(alpha_div2, chlor18, by = "Sample.Name")

ggplot(data = meta18, aes(x = O2ml_L, y = alpha.div.sh)) + geom_point()

ggplot(data = meta18, aes(x = Depthm, y = log(ChlorA))) + geom_point()

ggplot(data = meta18, aes(x = depth)) + geom_bar()


### northern most line (temp, depth, o2, div)
landt1 <- meta18 |>
  select(Sample.Name, alpha.div.sh, Cruise, line, Depthm, Distance, T_degC, O2ml_L) |>
  filter(line == "076.7")

ggplot(data = landt1, aes(x = Depthm, y = alpha.div.sh)) + geom_point()

ggplot(data = landt1, aes(x = Distance, y = alpha.div.sh)) + geom_point()

ggplot(data = landt1, aes(x = T_degC, y = alpha.div.sh)) + geom_point()

ggplot(data = landt1, aes(x = Depthm, y = T_degC)) + geom_point()

ggplot(data = landt1, aes(x = T_degC, y = O2ml_L)) + geom_point()

### southern most line
landt2 <- meta18 |>
  select(Sample.Name, alpha.div.sh, Cruise, line, Distance, Depthm, T_degC) |>
  filter(line == "093.3")

ggplot(data = landt2, aes(x = Depthm, y = alpha.div.sh)) + geom_point()

ggplot(data = landt2, aes(x = T_degC, y = alpha.div.sh)) + geom_point()

ggplot(data = landt2, aes(x = Depthm, y = T_degC)) + geom_point()


## north line nearshore by diversity
meta3 <- meta2 %>%
  select(Sample.Name, Cruise, line, sta, Distance, Depthm, Lon_Dec, ChlorA, NO3ug, alpha.div.sh, T_degC) %>%
  filter(ChlorA != 0) %>%
  mutate(station = as.numeric(sta),
         nearshore = "na",
         nearshore = if_else(station <= 55 | line == "081.8", "close", nearshore),
         nearshore = if_else(station > 55 & station < 80 , "middle", nearshore),
         nearshore = if_else(station >= 80 , "not close", nearshore))

ggplot(data = meta3, aes(x = nearshore)) + geom_bar()



landt3 <- meta18  |>
  filter(line == "076.7")

ggplot(data = landt3, aes(x = nearshore, y = alpha.div.sh)) + geom_boxplot()

ggplot(data = landt3, aes(x = Distance, y = alpha.div.sh)) + geom_point()

#south line " "
landt4 <- meta3  |>
  filter(line == "093.3")

ggplot(data = landt4, aes(x = nearshore, y = alpha.div.sh)) + geom_boxplot()

#middle line ""
landt5 <- meta3  |>
  filter(line == "083.3")

ggplot(data = landt5, aes(x = nearshore, y = alpha.div.sh)) + geom_boxplot()

### ???
ggplot(data = landt4, aes(x = Depthm, y = nearshore)) + geom_boxplot()

ggplot(data = landt3, aes(x = Depthm, y = nearshore)) + geom_point()

ggplot(data = landt5, aes(x = Depthm, y = T_degC)) + geom_point()

### whale whale whale.....

whale <- density_estimates
  
ggplot(data = whale, aes(x = season, y = bp)) + geom_boxplot()

ggplot(data = whale, aes(x = season, y = bm)) + geom_boxplot()

ggplot(data = whale, aes(x = season, y = mn)) + geom_boxplot() + 
  scale_y_continuous(limits = c(0,40))


ggplot(data = alphaszn, aes(x = logdepth, y = alpha, color = season)) +
  geom_point() + facet_wrap(~ season, ncol = 4) +
  scale_color_manual(
    values = c("winter" = "#00B8E7", "spring" = "#39B600", 
               "summer" = "#D89000", "fall" = "#FF5733"))

ggplot(data = alphaszn, aes(x = logdepth, y = alpha, color = season)) +
  geom_point() + facet_wrap(~ season, ncol = 4) +
  scale_color_manual(
    values = c("winter" = "#00B8E7", "spring" = "#39B600", 
               "summer" = "#D89000", "fall" = "#FF5733"))

meta182 <- meta18 |>
  rename(cruise = Cruise)

inner_join(meta182, alphasam4, by = cruise)


#nitrate and chlor and div
ggplot(data = meta18, aes(x = NO3ug, y = ChlorA)) + 
  geom_point()

ggplot(data = meta18, aes(x = NO3ug, y = ChlorA)) + 
  geom_point() + 
  scale_y_continuous(limits = c(0,7.5)) + 
  scale_x_continuous(limits = c(0,7.5))

ggplot(data = meta18, aes(x = NO3ug, y = alpha.div.sh)) + 
  geom_point() + 
  scale_y_continuous(limits = c(3,7.5)) + 
  scale_x_continuous(limits = c(0,7.5)) 

nitrate.fit <- lm(NO3ug~ChlorA, data = meta18)

summary(nitrate.fit)


#### processing meta for whales #####

# 
meta_avg <- metadata18 |>
  select(Sample.Name, Cruise, Sta_ID, Depthm, Distance, ChlorA, O2ml_L, NO3ug) |>
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') |>
  group_by(Cruise) |>
  summarize(o2 = mean(O2ml_L),
            chlor = mean(ChlorA),
            dist = mean(Distance), 
            nitrate = mean(NO3ug))

meta_avg <- meta_avg %>%
  mutate(cruise = as.character(cruise))

unique_cruises <- unique(metadata18$Cruise)
num_cruises <- length(unique_cruises)
print(num_cruises)

whale <- whale |>
  mutate(cruise = substr(cruise, 2, nchar(cruise)))

metawhale <- left_join(meta_avg, whale, by = 'cruise')

metawhale <- metawhale |>
  select(cruise, season, dist, bp, bm, mn, o2, chlor, nitrate)


### same averaging for imputed values (or maybe non imputed idk)
whale_imp <- density_estimates_imputed |>
  mutate(cruise = substr(cruise, 2, nchar(cruise)))

metawhale_imp <- left_join(meta_avg, whale_imp, by = 'cruise')


# chlor and blue(?) whale densities 
ggplot(data = metawhale_imp, aes(x = chlor, y = bm)) + 
  geom_point()

# o2 and 
ggplot(data = metawhale_imp, aes(x = o2, y = bm.imp)) + 
  geom_point()

