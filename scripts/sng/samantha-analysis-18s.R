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
  select(Sample.Name, Cruise, Depthm, CC_Depth, Time, Lat_Dec, Lon_Dec, Distance, NO3ug, NH3ug, ChlorA, T_degC)

meta18 <- inner_join(alpha_div2, chlor18, by = "Sample.Name")

### northern most line
landt1 <- meta18 |>
  select(Sample.Name, alpha.div.sh, Cruise, line, Depthm, T_degC) |>
  filter(line == "076.7")

ggplot(data = landt1, aes(x = Depthm, y = alpha.div.sh)) + geom_point()

ggplot(data = landt1, aes(x = T_degC, y = alpha.div.sh)) + geom_point()

ggplot(data = landt1, aes(x = Depthm, y = T_degC)) + geom_point()

### southern most line
landt2 <- meta18 |>
  select(Sample.Name, alpha.div.sh, Cruise, line, Depthm, T_degC) |>
  filter(line == "093.3")

ggplot(data = landt2, aes(x = Depthm, y = alpha.div.sh)) + geom_point()

ggplot(data = landt2, aes(x = T_degC, y = alpha.div.sh)) + geom_point()

ggplot(data = landt2, aes(x = Depthm, y = T_degC)) + geom_point()

## north line nearshore by diversity
meta3 <- meta2 %>%
  select(Sample.Name, Cruise, line, sta, Depthm, Lon_Dec, ChlorA, NO3ug, alpha.div.sh, T_degC) %>%
  filter(ChlorA != 0) %>%
  mutate(station = as.numeric(sta),
         nearshore = "na",
         nearshore = if_else(station <= 55 | line == "081.8", "close", nearshore),
         nearshore = if_else(station > 55 & station < 80 , "middle", nearshore),
         nearshore = if_else(station >= 80 , "not close", nearshore))

landt3 <- meta3  |>
  filter(line == "093.3")

ggplot(data = landt3, aes(x = nearshore, y = alpha.div.sh)) + geom_boxplot()

ggplot(data = landt3, aes(x = nearshore, y = alpha.div.sh)) + geom_point()

#south line " "
landt3 <- meta3  |>
  filter(line == "093.3")

ggplot(data = landt3, aes(x = nearshore, y = alpha.div.sh)) + geom_boxplot()
