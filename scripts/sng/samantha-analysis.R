# further explorations by sam 
library(ggplot2)
library(tidyverse)
# grouping by station
alpha_div %>%
  group_by(depthm) %>%
  summarize(alpha = mean(alpha.div.sh)) %>%
  ggplot(aes(x = depthm, y = alpha, color = line)) +
  geom_line() +
  theme(axis.text.x = element_text(angle = 90))

#X cleaning for grouping by depth to compare depth and diversity 
alphasam <- alpha_div |>
  mutate(depthm == as.numeric(depthm)) |>
  group_by(depthm) |>
  summarise(alpha = mean(alpha.div.sh))

ggplot(data = alphasam, mapping = aes(x= alpha, y= depthm), by = 6) +
  geom_point() +
  scale_y_continuous(n.breaks = 17)

# clean
alphasam2 <- alpha_div %>%
  mutate(depthm = as.numeric(depthm)) %>%
  group_by(depthm, line) %>%
  summarise(alpha = mean(alpha.div.sh)) %>%
  filter(depthm <= 171) 

# diversity by depth by line
ggplot(data = alphasam2, aes(x = alpha, y = depthm, color = line)) +
  geom_point() +
  scale_y_continuous(breaks = seq(min(alphasam2$depthm), 
                                  max(alphasam2$depthm), 
                                  by = 17)) + scale_y_reverse() + facet_wrap(~ line)

# stacking, by quarter(color)

# diversity by depth with loess / lm line
ggplot(data = alphasam2, aes(x = depthm, y = alpha)) +
  geom_point() +
  scale_y_continuous(breaks = seq(min(alphasam2$depthm), 
                                  max(alphasam2$depthm), 
                                  by = 17)) +
  geom_smooth() #+ geom_smooth(method = lm)


# diversity by line boxplot
ggplot(data = alpha_div, aes(x = alpha.div.sh, y = line)) + 
  geom_boxplot()
# bin by depth 

# diversity v depth by cruise
alphacru <- alpha_div %>%
  mutate(depthm = as.numeric(depthm)) %>%
  group_by(depthm, cruise) %>%
  summarise(alpha = mean(alpha.div.sh)) %>%
  filter(depthm <= 171) 

ggplot(data = alphacru, aes(x = depthm, 
                            y = alpha, 
                            color = cruise)) +
  geom_point(size = 0.5) +
  scale_y_continuous(breaks = seq(min(alphacru$depthm), 
                                  max(alphacru$depthm), 
                                  by = 34)) +
  facet_wrap(~ cruise, nrow = 4, dir = 'v') + 
  geom_hline(aes(yintercept = mean(alpha), group = cruise), color = 'black')# + facet_grid()
# facet columns are by quarter, by transect
# try color by season or transect

# another diversity by cruise by alpha (not as good looking)
ggplot(data = alphacru, aes(x = cruise, 
                            y = alpha, 
                            size = depthm)) +
  geom_point(fill = 'white') +
  scale_y_continuous(breaks = seq(min(alphacru$depthm), 
                                  max(alphacru$depthm), 
                                  by = 24)) 

#sorting months and seasons
alphasam3 <- alpha_div %>%
  select(cruise, line, sta, depthm, alpha.div.sh) %>%
  mutate(cruise = as.numeric(str_sub(cruise, -6))) %>%
  mutate(logdepth = log(as.numeric(depthm))) %>%
  mutate(month = str_sub(cruise, -2), season = "fall") %>%
  mutate(season = if_else(month %in% c("01", "02"), "winter", season)) %>%
  mutate(season = if_else(month %in% c("03", "04", "05"), "spring", season)) %>%
  mutate(season = if_else(month %in% c("06", "07", "08"), "summer", season)) %>%
  mutate(year = str_sub(cruise, 1, 4))

alphasam4 <- alpha_div %>%
  select(cruise, line, sta, depthm, alpha.div.sh) %>%
  mutate(cruise = as.numeric(str_sub(cruise, -6))) %>%
  mutate(month = str_sub(cruise, -2), season = "fall") %>%
  mutate(season = if_else(month %in% c("01", "02"), "winter", season)) %>%
  mutate(season = if_else(month %in% c("03", "04", "05"), "spring", season)) %>%
  mutate(season = if_else(month %in% c("06", "07", "08"), "summer", season)) %>%
  mutate(year = str_sub(cruise, 1, 4))

ggplot(data = alphasam4, aes(x = season, y = alpha.div.sh)) + 
  geom_boxplot()

alphaszn <- alphasam3 %>%
  select(season, line, logdepth, alpha.div.sh) %>%
  group_by(logdepth, season, line) %>%
  summarise(alpha = mean(alpha.div.sh)) %>%
  filter(logdepth >= 0 & logdepth <= 7)

# diversity by season
ggplot(data = alphaszn, aes(x = logdepth, y = alpha, color = season)) +
  geom_point() + facet_wrap(~ line)

# only by season, looks cleaner
ggplot(data = alphaszn, aes(x = logdepth, y = alpha, color = season)) +
  geom_point() + facet_wrap(~ season, nrow = 1)

# 
levels(season$alphaszn)
season$alphaszn <- factor(season$alphaszn, levels = c("winter", 
                                                      "spring", 
                                                      "summer", 
                                                      "fall"))

ggplot(data = alphaszn, aes(x = logdepth, y = alpha, color = season)) +
  geom_point() + facet_wrap(~ season, ncol = 4) +
  scale_color_manual(
    values = c("winter" = "#00B8E7", "spring" = "#39B600", 
               "summer" = "#D89000", "fall" = "#FF5733"))

# new diversity by depth by cruise (months and years)
alphacru2 <- alphasam3 %>%
  select(cruise, logdepth, alpha.div.sh) %>%
  group_by(logdepth, season, line) %>%
  summarise(alpha = mean(alpha.div.sh)) %>%
  filter(logdepth >= 0 & logdepth <= 7)


ggplot(data = alphacru2, aes(x = depthm, 
                             y = alpha, 
                             color = cruise)) +
  geom_point(size = 0.5) +
  scale_y_continuous(breaks = seq(min(alphacru$depthm), 
                                  max(alphacru$depthm), 
                                  by = 34)) +
  facet_wrap(~ cruise) + facet_grid()

# pairwise beta diversity for one cruise
edna16s_reads %>%
  filter(cruise == cruises[1], as.numeric(depthm) < 30) %>%
  betadiver() %>%
  plot()


# further explorations?
# examining other variables 
alpha_div <- alpha_div |>
  rename(Sample.Name = sample.id)

alpha_div2 <- alpha_div |>
  select(Sample.Name, line, sta, depthm, alpha.div.sh) |>
  mutate(Sample.Name = substr(Sample.Name, 2, nchar(alpha_div2$Sample.Name)))


chlor <- metadata %>%
  select(Sample.Name, Cruise, Depthm, Time, Lat_Dec, Lon_Dec, Distance, NO3ug, NH3ug, ChlorA, T_degC)

meta2 <- inner_join(alpha_div2, chlor, by = "Sample.Name")



##### Investigating Nitrate, Distance from shore, Chlorophyll, and line #####

# Nitrate measurements compared to diversity scatter
ggplot(data = meta2, aes(x = NO3ug, y = alpha.div.sh)) + 
  geom_point() + 
  scale_y_continuous(limits = c(0,6))
# closer look 
meta2 %>%
  filter(NO3ug != 0) %>%
  ggplot(aes(x = NO3ug, y = alpha.div.sh)) + 
  geom_point() + 
  scale_y_continuous(limits = c(3,6)) +
  scale_x_continuous(limits = c(0, 10))

meta2 %>%
  filter(NO3ug != 0, NO3ug <10) %>%
  mutate(nitrate = as.numeric(NO3ug),
         nitrate.fac = cut_number(nitrate, 4)) %>%
  ggplot(aes(x = alpha.div.sh, y = nitrate.fac)) + 
  geom_boxplot()

meta2 %>%
  filter(NO3ug == 0) %>%
  ggplot(aes(x = alpha.div.sh)) + 
  geom_histogram()


# Nitrate by Chlorophyll scatter
ggplot(data = meta2, aes(x = NO3ug, y = ChlorA)) + 
  geom_point()
#  closer look
ggplot(data = meta2, aes(x = NO3ug, y = ChlorA)) + 
  geom_point() + 
  scale_y_continuous(limits = c(0,7.5)) + 
  scale_x_continuous(limits = c(0,7.5))

# boxplot distance from shore by diversity
meta2 %>%
  mutate(dist = as.numeric(Distance),
         dist.fac = cut_number(dist, 4)) %>%
  ggplot(aes(x = dist.fac, y = alpha.div.sh)) + 
  geom_boxplot()

meta2 %>%
  mutate(dist = as.numeric(Distance),
         dist.fac = cut_number(dist, 4)) %>%
  ggplot(aes(x = depthm, y = alpha.div.sh, color = dist.fac)) + 
  geom_point() + facet_wrap(~dist.fac)

## fix axes
# ggplot(data = meta2, aes(x = Distance, y = depthm, color = alpha.div.sh)) + geom_point()


# Nitrate and distance from shore
meta2 %>%
  mutate(dist = as.numeric(Distance),
         dist.fac = cut_number(dist, 2)) %>%
  filter(!is.na(dist.fac)) %>%
  ggplot(aes(x = dist.fac, y = NO3ug)) + 
  geom_boxplot()

meta2 %>%
  mutate(dist = as.numeric(Distance),
         dist.fac = cut_number(dist, 4)) %>%
  filter(!is.na(dist.fac)) %>%
  ggplot(aes(x = dist.fac, y = alpha.div.sh)) + 
  geom_boxplot()

# Chlorophyll and distance from shore
meta2 %>%
  mutate(dist = as.numeric(Distance),
         dist.fac = cut_number(dist, 4)) %>%
  ggplot(aes(x = ChlorA, y = dist.fac)) + 
  geom_boxplot()

# chlorophyll by line
ggplot(data = meta2, aes(x = line, y = ChlorA)) + geom_boxplot()


#### investigating closeness to shore with these factors

meta3 <- meta2 %>%
  select(Sample.Name, Cruise, line, sta, Depthm, Lon_Dec, ChlorA, NO3ug, alpha.div.sh, T_degC) %>%
  filter(ChlorA != 0) %>%
  mutate(station = as.numeric(sta),
         nearshore = "na",
         nearshore = if_else(station <= 55 | line == "081.8", "close", nearshore),
         nearshore = if_else(station > 55 & station < 80 , "middle", nearshore),
         nearshore = if_else(station >= 80 , "not close", nearshore))

# *chlorophyll based on closeness to shore
ggplot(data = meta3, aes(x = nearshore, y = ChlorA)) + geom_boxplot()

# *diversity based on closeness to shore
ggplot(data = meta3, aes(x = nearshore, y = alpha.div.sh)) + geom_boxplot()

# *nitrate based on closeness to shore
ggplot(data = meta3, aes(x = nearshore, y = NO3ug)) + geom_boxplot()

# *ammonia explorations
meta2 %>%
  mutate(dist = as.numeric(Distance),
         dist.fac = cut_number(dist, 4)) %>%
  filter(!is.na(dist.fac)) %>%
  ggplot(aes(x = dist.fac, y = NH3ug)) + 
  geom_boxplot()

ggplot(data = meta2, aes(x = NH3ug, y = alpha.div.sh)) + 
  geom_point() + 
  scale_y_continuous(limits = c(0,6))
# closer look 
meta2 %>%
  filter(NH3ug != 0) %>%
  ggplot(aes(x = NH3ug, y = alpha.div.sh)) + 
  geom_point() + 
  scale_y_continuous(limits = c(3,6)) +
  scale_x_continuous(limits = c(0, 2))

meta2 %>%
  filter(NH3ug != 0) %>%
  mutate(ammon = as.numeric(NO3ug),
         am.fac = cut_number(ammon, 4)) %>%
  ggplot(aes(x = alpha.div.sh, y = am.fac)) + 
  geom_boxplot()

ggplot(data = meta2, aes(x = NH3ug, y = NO3ug)) + geom_point()

ggplot(data = meta2, aes(x = NH3ug, y = alpha.div.sh)) + geom_point()



# temperature by depth by line 
ggplot(data = meta2, aes(x = Depthm, y = T_degC, color = line)) + 
  geom_point() + 
  facet_wrap(~line)

# water temp
meta2 %>%
  mutate(T_degC = as.numeric(T_degC),
         temp = cut_number(T_degC, 5)) %>%
  filter(!is.na(temp)) %>%
  ggplot(aes(x = temp, y = alpha.div.sh)) + 
  geom_boxplot()

ggplot(data = meta3, aes(x = T_degC, y = alpha.div.sh)) +
  geom_point() + 
  facet_wrap(~line)

# closeness to shore by temp and diversity, by line 
meta3 %>%
  filter(line != "081.8") %>%
  ggplot(aes(x = nearshore, y = T_degC)) + geom_boxplot() + facet_wrap(~line)

meta3 %>%
  filter(line != "081.8") %>%
  ggplot(aes(x = nearshore, y = alpha.div.sh)) + geom_boxplot() + facet_wrap(~line)

########## transition to 18s data ############

edna18s <- edna_samples




curve(dgamma(x, shape = 10, scale = 15), from = 0, to = 200)




