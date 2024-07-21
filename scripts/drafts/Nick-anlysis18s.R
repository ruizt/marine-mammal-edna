### Check out the 9 asvs by season
### Recompute alpha div, beta div
### Check to see same patterns over 18s as we saw in 16s


########################################################################

### 18s Data Exploration

load('data/ncog-18s.RData')

## BY SEASON

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
edna18s_reads <- edna_samples %>%
  filter(as.numeric(line) %in% transects)

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
  bind_cols(alpha.div.sh = diversity(asv, index = 'shannon'))



alpha_div_18_season <- alpha_div18 |> 
  mutate(month = as.numeric(substring(cruise,6,8)))

alpha_div_18_season <-  alpha_div_18_season |> 
  mutate(season = as.factor(case_when(
    month %in% c(12,1,2) ~ "Winter",
    month %in% c(3,4,5) ~ "Spring",
    month %in% c(6,7,8) ~ "Summer",
    month %in% c(9,10,11) ~ "Fall")))

alpha_div_18_season <-  alpha_div_18_season |> 
  filter(as.numeric(depthm) < 500) |> 
  mutate(season = fct_relevel(season, "Winter", "Spring", "Summer", "Fall"))


alpha_div_18_season |> 
  filter(as.numeric(depthm) < 500) |> 
  ggplot(aes(x = season, y = alpha.div.sh)) +
  geom_boxplot()


alpha_div_18_season |> 
  mutate( depth.rng = cut_number(as.numeric(depthm),3)) |> 
  ggplot(aes(x = depth.rng, y = alpha.div.sh)) +
  geom_boxplot()

alpha_div_18_season |> 
  mutate( depth.rng = cut_number(as.numeric(depthm),3)) |> 
  ggplot(aes(x = depth.rng, y = alpha.div.sh)) +
  geom_boxplot() + facet_wrap(~season)


alpha_div_18_season |> 
  mutate( depth.rng = cut_number(as.numeric(depthm),3)) |> 
  ggplot(aes(x = season, y = alpha.div.sh)) +
  geom_boxplot() + facet_wrap(~depth.rng)

# Overall Same pattern as 16s
# Winter [0,10] only difference in trend
# Much more variation in 18s alpha div 

alpha_div_18_season |> 
  ggplot(aes(x = as.numeric(depthm), y = alpha.div.sh)) +
  geom_point()

alpha_div_18_season |> 
  group_by(depthm) |> 
  summarise(avgDiv = mean(as.numeric(alpha.div.sh))) |> 
  ggplot(aes(x = as.numeric(depthm), y = avgDiv)) +
  geom_point()

# 9 taxa from spls model

spls = read_csv("/Users/nicholaspatrick/Desktop/eDNA Project/rslt/taxa-weights-18s-spls.csv")

# Beta Diversity in Winter
edna18s_reads |> 
  filter(as.numeric(substring(cruise, 6, 8)) %in% c(12,1,2)) |> 
  select(cruise, line, sta, depthm, all_of(spls$short.id)) |> 
  betadiver() |> 
  plot() |>  title(main = "Winter", xlab = NULL)

# Beta Diversity in Spring
edna18s_reads |> 
  filter(as.numeric(substring(cruise, 6, 8)) %in% c(3,4,5)) |> 
  select(cruise, line, sta, depthm, all_of(spls$short.id)) |> 
  betadiver() |> 
  plot() |> title(main = "Spring", xlab = NULL)

# Beta Diversity in Summer
edna18s_reads |> 
  filter(as.numeric(substring(cruise, 6, 8)) %in% c(6,7,8)) |> 
  select(cruise, line, sta, depthm, all_of(spls$short.id)) |> 
  betadiver() |> 
  plot() |>  title(main = "Summer", xlab = NULL)

# Beta Diversity in Fall
edna18s_reads |> 
  filter(as.numeric(substring(cruise, 6, 8)) %in% c(9,10,11)) |> 
  select(cruise, line, sta, depthm, all_of(spls$short.id)) |> 
  betadiver() |> 
  plot() |> title(main = "Fall", xlab = NULL)


# alpha diversity of those 9
filtered_asvs <- edna18s_reads  |> 
  select(all_of(spls$short.id))

alpha_div_9 <- sampinfo %>%
  bind_cols(alpha.div.sh = diversity(filtered_asvs, index = 'shannon'))

# add season and depth range
alpha_div_9 <- alpha_div_9 |> 
  mutate( month = as.numeric(substring(cruise, 6, 8)))

alpha_div_9 <- alpha_div_9 |> 
  mutate(season = as.factor(case_when(
    month %in% c(12,1,2) ~ "Winter",
    month %in% c(3,4,5) ~ "Spring",
    month %in% c(6,7,8) ~ "Summer",
    month %in% c(9,10,11) ~ "Fall")))

alpha_div_9 <- alpha_div_9 |> 
  filter(as.numeric(depthm) < 500) |> 
  mutate(season = fct_relevel(season, "Winter", "Spring", "Summer", "Fall"),
         depth.rng = cut_number(as.numeric(depthm), 3))


# Alpha diversity across season

alpha_div_9 |> 
  ggplot(aes(x=season, y = alpha.div.sh)) + 
  geom_boxplot() + labs(title = "9 asvs Alpha Div by Season")


# Alpha diversity across depth
alpha_div_9 |> 
  ggplot(aes(x = depth.rng, y = alpha.div.sh)) + 
  geom_boxplot()

# Alpha Div by Depth scatter
alpha_div_9 |>
  mutate(depthm = as.numeric(depthm)) |> 
  ggplot(aes(x = depthm, y = alpha.div.sh)) +
  geom_point() + facet_wrap(~season) +
  geom_hline(yintercept = 1, col = "blue")
# AVG alpha div by depth scatter

alpha_div_9 |> 
  mutate(depthm = as.numeric(depthm)) |> 
  group_by(depthm,season) |> 
  summarize(avgDiv = mean(alpha.div.sh)) |> 
  ggplot(aes(x = depthm, y = avgDiv)) + 
  geom_point() + facet_wrap(~season) +
  geom_hline(yintercept = 1, col = "blue")
# Alpha diversity across depth and season 

alpha_div_9 |> 
  ggplot(aes(x= depth.rng, y = alpha.div.sh)) + 
  geom_boxplot() + facet_grid(~season) + 
  geom_hline(yintercept = 1, col = "blue") +
  xlab("Depth") +
  ylab("Alpha Diversity") +
  labs(title = "9 ASVs Alpha Diversity", subtitle = "Season")

# Alpha diversity across season and depth 
alpha_div_9 |> 
  ggplot(aes(x= season, y = alpha.div.sh)) + 
  geom_boxplot() + facet_grid(~depth.rng) + 
  geom_hline(yintercept = 1, col = "blue") +
  xlab("Season") +
  ylab("Alpha Diversity") +
  labs(title = "9 ASVs Alpha Diversity", subtitle = "Depth")

# Join metadata with 9asv alpha div
meta18 <- metadata |> 
  separate(Sta_ID, into = c('line', 'sta'), sep = ' ') |> 
  mutate(Depthm = as.character(Depthm),
         Cruise = paste("X", Cruise, sep = "")) |> 
  left_join(alpha_div_9, by = join_by(Cruise == cruise, line, sta, Depthm == depthm) )


# Oxygen vs alpha div

meta18 |> 
  ggplot(aes(x = as.numeric(Depthm), y = O2ml_L)) +
  geom_point()
