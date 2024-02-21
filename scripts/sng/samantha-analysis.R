# further explorations by sam 
library(ggplot2)
theme_set(theme_grey())
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
ggsave()
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

# read docs
?betadiver

# indices supported
betadiver(help = T)

# further explorations?
# examining other variables 
alpha_div <- alpha_div |>
  rename(sample.id = Sample.Name)


chlor <- metadata |>
  select(Sample.Name, Cruise, Depthm, Time, Distance, NO3ug, NH3ug, ChlorA)

full_join(alpha_div, chlor, by =)


