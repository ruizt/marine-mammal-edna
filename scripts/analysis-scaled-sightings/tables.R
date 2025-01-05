library(tidyverse)
library(lubridate)

## DIRECTORIES -----------------------------------------------------------------
data_dir <- 'data/processed/'
model_dir <- 'rslt/models/scaled-sightings/'
stbl_dir <- 'rslt/stability-selection/'
val_dir <- 'rslt/nested-validation/'

out_dir <- 'rslt/tbl/'
fs::dir_create(out_dir)

## SUPPLEMENTARY TABLE 2: AMPLICON DATA SUMMARY (LONG) -------------------------

## 18SV9 ## 

paste(data_dir, 'ncog18sv9.RData', sep = '') |> load()

cruise_meta_18sv9 <- sample_metadata |>
  select(sample.id, cruise, sta.id, datetime, depthm) |>
  mutate(datetime = mdy_hm(datetime)) |>
  separate(sta.id, into = c('line', 'station'), sep = ' ') |>
  group_by(cruise) |>
  summarize(start.date = min(datetime) |> as_date(),
            end.date = max(datetime) |> as_date(),
            min.depth = min(depthm),
            max.depth = max(depthm),
            n.samples = n())

ncog_summary_18sv9 <- sample_metadata |>
  select(sample.id, cruise, sta.id, datetime, depthm) |>
  mutate(datetime = mdy_hm(datetime)) |>
  separate(sta.id, into = c('line', 'station'), sep = ' ') |>
  count(cruise, line) |>
  pivot_wider(names_from = line, values_from = n) |>
  left_join(cruise_meta_18sv9) |>
  mutate(cruise = ym(cruise), 
         year = year(cruise)) |>
  group_by(year) |> 
  mutate(qtr = paste('Q', row_number(), sep = ''),
         season = factor(qtr, levels = paste('Q', 1:4, sep = ''),
                         labels = c('winter', 'spring', 'summer', 'fall')),
         marker = '18SV9') |>
  select(marker, year, season, start.date, end.date, min.depth, max.depth, 
         where(is.integer), n.samples)

## 18SV4 ## 

paste(data_dir, 'ncog18sv4.RData', sep = '') |> load()

cruise_meta_18sv4 <- sample_metadata |>
  select(sample.id, cruise, sta.id, datetime, depthm) |>
  mutate(datetime = mdy_hm(datetime)) |>
  separate(sta.id, into = c('line', 'station'), sep = ' ') |>
  group_by(cruise) |>
  summarize(start.date = min(datetime) |> as_date(),
            end.date = max(datetime) |> as_date(),
            min.depth = min(depthm),
            max.depth = max(depthm),
            n.samples = n())

ncog_summary_18sv4 <- sample_metadata |>
  select(sample.id, cruise, sta.id, datetime, depthm) |>
  mutate(datetime = mdy_hm(datetime)) |>
  separate(sta.id, into = c('line', 'station'), sep = ' ') |>
  count(cruise, line) |>
  pivot_wider(names_from = line, values_from = n) |>
  left_join(cruise_meta_18sv4) |>
  mutate(cruise = ym(cruise), 
         year = year(cruise)) |>
  group_by(year) |> 
  mutate(qtr = paste('Q', row_number(), sep = ''),
         season = factor(qtr, levels = paste('Q', 1:4, sep = ''),
                         labels = c('winter', 'spring', 'summer', 'fall')),
         marker = '18SV4') |>
  select(marker, year, season, start.date, end.date, min.depth, max.depth, 
         where(is.integer), n.samples)

## 16S ## 

paste(data_dir, 'ncog16s.RData', sep = '') |> load()

cruise_meta_16s <- sample_metadata |>
  select(sample.id, cruise, sta.id, datetime, depthm) |>
  mutate(datetime = mdy_hm(datetime)) |>
  separate(sta.id, into = c('line', 'station'), sep = ' ') |>
  group_by(cruise) |>
  summarize(start.date = min(datetime) |> as_date(),
            end.date = max(datetime) |> as_date(),
            min.depth = min(depthm),
            max.depth = max(depthm),
            n.samples = n())

ncog_summary_16s <- sample_metadata |>
  select(sample.id, cruise, sta.id, datetime, depthm) |>
  mutate(datetime = mdy_hm(datetime)) |>
  separate(sta.id, into = c('line', 'station'), sep = ' ') |>
  count(cruise, line) |>
  pivot_wider(names_from = line, values_from = n) |>
  left_join(cruise_meta_16s) |>
  mutate(cruise = ym(cruise), 
         year = year(cruise)) |>
  group_by(year) |> 
  mutate(qtr = paste('Q', row_number(), sep = ''),
         season = factor(qtr, levels = paste('Q', 1:4, sep = ''),
                         labels = c('winter', 'spring', 'summer', 'fall')),
         marker = '16S') |>
  select(marker, year, season, start.date, end.date, min.depth, max.depth, 
         where(is.integer), n.samples)

# combine into long-form summary (number of samples per line)
ncog_summary <- bind_rows(ncog_summary_16s, 
                          ncog_summary_18sv4, 
                          ncog_summary_18sv9)

## TABLE 1a-b: SAMPLING SUMMARIES ----------------------------------------------

# aggregate previous table into short form summary
ncog_summary_short <- ncog_summary |>
  select(-n.samples) |>
  pivot_longer(where(is.integer)) |>
  drop_na() |>
  group_by(marker, year, season, start.date, end.date) |>
  summarize(n.samples = sum(value, na.rm = T),
            n.lines = n(),
            .groups = 'drop') |>
  pivot_wider(names_from = marker, values_from = starts_with('n')) |>
  select(year, season, ends_with('date'), ends_with('16S'), ends_with('18SV4'),
         ends_with('18SV9'))
  

## TABLE 1b: VISUAL SIGHTING SUMMARY

visual_sightings <- read_csv('data/_raw/CalCOFI_2004-2021_CombinedSightings.csv') |>
  rename_with(~str_remove_all(.x, '[:punct:]') |> 
                str_squish() |> 
                str_replace_all(' ', '.') |> 
                tolower()) |>
  mutate(species = tolower(species.1)) |>
  filter(adjusted.both.on.effort.and.on.transect == 'ON',
         species %in% c('bm', 'bp', 'mn')) |>
  mutate(datetime = mdy_hm(datetime.local),
         year = year(datetime),
         season = factor(season, levels = c('winter', 'spring', 'summer', 'fall'))) |>
  select(year, season, datetime, species, best) |>
  filter(year >= 2014) 

sighting_counts <- visual_sightings |>
  count(year, season, species) |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback'))) |>
  pivot_wider(names_from = species, values_from = n) |> 
  mutate(across(where(is.integer), ~replace_na(.x, replace = 0L)))

sighting_summary <- visual_sightings |>
  group_by(year, season) |>
  summarize(start.date = min(date(datetime)),
            end.date = max(date(datetime))) |>
  left_join(sighting_counts)

## SUPP TABLES 3a-c ------------------------------------------------------------

# candidate asvs from 18sv9
paste(data_dir, 'ncog18sv9.RData', sep = '') |> load()
asv_taxa_18sv9 <- asv_taxa 

# candidate asvs from 18sv4
paste(data_dir, 'ncog18sv4.RData', sep = '') |> load()
asv_taxa_18sv4 <- asv_taxa

# candidate asvs from 16s
paste(data_dir, 'ncog16s.RData', sep = '') |> load()
asv_taxa_16s <- asv_taxa

## TABLE 2: MODEL FIT ----------------------------------------------------------

## 18SV9 ##

# load model info
paste(model_dir, 'fitted-models-18sv9.RData', sep = '') |> load()

# selected asvs
sel_asv_18sv9 <- fitted_models |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback')),
         marker = '18SV9',
         coef.ss = map(coef, ~2^.x)) |>
  select(species, coef, coef.ss, ss) |>
  unnest(everything()) |>
  rename(asv = ss,
         coef.lr = coef) |>
  left_join(asv_taxa_18sv9, join_by(asv == short.id)) |>
  select(-asv, -feature.id) |>
  arrange(species, desc(coef.lr))

# number of unique families among selected asvs
nfam_18sv9 <- sel_asv_18sv9 |>
  group_by(species) |>
  distinct(f) |>
  count(name = 'n.family')

# number of unique orders among selected asvs
nord_18sv9 <- sel_asv_18sv9 |>
  group_by(species) |>
  distinct(o) |>
  count(name = 'n.order')

# number of unique classes among selected asvs
nclass_18sv9 <- sel_asv_18sv9 |>
  group_by(species) |>
  distinct(c) |>
  count(name = 'n.class')

# retreive hyperparameter info
parms_18sv9 <- paste(stbl_dir, '18sv9-ss/stable-sets.rds', sep = '') |> 
  read_rds() |>
  mutate(eta.range = paste('[', round(eta.min, 2), ', ', round(eta.max, 2), ']', sep = '')) |>
  select(species, ncomp, eta.range)

# join with model fit metrics
model_fit_18sv9 <- model_metrics |>
  left_join(parms_18sv9) |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback')),
         marker = '18SV9') |>
  left_join(nfam_18sv9) |>
  left_join(nord_18sv9) |>
  left_join(nclass_18sv9) |>
  select(species, marker, ncomp, eta.range, adj.rsq.lr, adj.rsq.ss, n.asv, n.family, n.order, n.class)

## 18SV4 ##

# load model info
paste(model_dir, 'fitted-models-18sv4.RData', sep = '') |> load()

# selected asvs
sel_asv_18sv4 <- fitted_models |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback')),
         marker = '18SV4',
         coef.ss = map(coef, ~2^.x)) |>
  select(species, coef, coef.ss, ss) |>
  unnest(everything()) |>
  rename(asv = ss,
         coef.lr = coef) |>
  left_join(asv_taxa_18sv4, join_by(asv == short.id)) |>
  select(-asv, -feature.id) |>
  arrange(species, desc(coef.lr))

# number of unique families among selected asvs
nfam_18sv4 <- sel_asv_18sv4 |>
  group_by(species) |>
  distinct(f) |>
  count(name = 'n.family')

# number of unique orders among selected asvs
nord_18sv4 <- sel_asv_18sv4 |>
  group_by(species) |>
  distinct(o) |>
  count(name = 'n.order')

# number of unique classes among selected asvs
nclass_18sv4 <- sel_asv_18sv4 |>
  group_by(species) |>
  distinct(c) |>
  count(name = 'n.class')

# retreive hyperparameter info
parms_18sv4 <- paste(stbl_dir, '18sv4-ss/stable-sets.rds', sep = '') |> 
  read_rds() |>
  mutate(eta.range = paste('[', round(eta.min, 2), ', ', round(eta.max, 2), ']', sep = '')) |>
  select(species, ncomp, eta.range)

# join with model fit metrics
model_fit_18sv4 <- model_metrics |>
  left_join(parms_18sv4) |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback')),
         marker = '18SV4') |>
  left_join(nfam_18sv4) |>
  left_join(nord_18sv4) |>
  left_join(nclass_18sv4) |>
  select(species, marker, ncomp, eta.range, adj.rsq.lr, adj.rsq.ss, n.asv, n.family, n.order, n.class)

## 16S ##

# load model info
paste(model_dir, 'fitted-models-16s.RData', sep = '') |> load()

# selected asvs
sel_asv_16s <- fitted_models |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback')),
         marker = '16S',
         coef.ss = map(coef, ~2^.x)) |>
  select(species, coef, coef.ss, ss) |>
  unnest(everything()) |>
  rename(asv = ss,
         coef.lr = coef) |>
  left_join(asv_taxa_16s, join_by(asv == short.id)) |>
  select(-asv, -feature.id) |>
  arrange(species, desc(coef.lr))

# number of unique families among selected asvs
nfam_16s <- sel_asv_16s |>
  group_by(species) |>
  distinct(f) |>
  count(name = 'n.family')

# number of unique orders among selected asvs
nord_16s <- sel_asv_16s |>
  group_by(species) |>
  distinct(o) |>
  count(name = 'n.order')

# number of unique classes among selected asvs
nclass_16s <- sel_asv_16s |>
  group_by(species) |>
  distinct(c) |>
  count(name = 'n.class')

# retreive hyperparameter info
parms_16s <- paste(stbl_dir, '16s-ss/stable-sets.rds', sep = '') |> 
  read_rds() |>
  mutate(eta.range = paste('[', round(eta.min, 2), ', ', round(eta.max, 2), ']', sep = '')) |>
  select(species, ncomp, eta.range)

# join with model fit metrics
model_fit_16s <- model_metrics |>
  left_join(parms_16s) |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback')),
         marker = '16S') |>
  left_join(nfam_16s) |>
  left_join(nord_16s) |>
  left_join(nclass_16s) |>
  select(species, marker, ncomp, eta.range, adj.rsq.lr, adj.rsq.ss, n.asv, n.family, n.order, n.class)

## COMBINE ##

model_fit <- bind_rows(model_fit_16s, model_fit_18sv4, model_fit_18sv9)

## TABLE 3: SELECTION CONSISTENCY ----------------------------------------------

## UNION/INTERSECTION FUNCTIONS ##

# function to compute "soft intersection" of sets in <set_list>
intersect_fn <- function(set_list, thresh){
  out <- tibble(asv = Reduce(c, set_list)) |>
    group_by(asv) |>
    count() |>
    filter(n >= thresh*length(set_list)) |>
    pull(asv)
  return(out)
}

# function to compute no. of elements in union of sets in <set_list>
union_fn <- function(set_list){
  out <- Reduce(c, set_list) |> unique()
  return(out)
}

# Function: find intersection based on taxonomic levels
# input: dataset containing validation id, species, asv, and taxonomic info
# output: 1 row per species containing intersection based on taxonomic level
intersect_fn_taxa <- function(data, tax.level, thresh){
  if(tolower(tax.level) == "c" | tolower(tax.level) == "class"){
    res <- data |> 
      group_by(species, outer.id, d,p,c) |> 
      drop_na() |> 
      unique() |> 
      group_by(outer.id, species) |> 
      summarize(class.list = list(c)) |> 
      mutate(class.list = map(class.list, unique)) |> 
      group_by(species) |> 
      summarise(intersect = list(intersect_fn(class.list, thresh)))
  }
  else if(tolower(tax.level) == "o" | tolower(tax.level) == "order"){
    res <- data |> 
      group_by(species, outer.id, d,p,c,o) |> 
      drop_na() |> 
      unique() |> 
      group_by(outer.id, species) |> 
      summarize(order.list = list(o)) |> 
      mutate(order.list = map(order.list, unique)) |> 
      group_by(species) |> 
      summarise(intersect = list(intersect_fn(order.list, thresh)))
  }
  else if(tolower(tax.level) == "f" | tolower(tax.level) == "family"){
    res <- data |> 
      group_by(species, outer.id, d,p,c,o,f) |> 
      drop_na() |> 
      unique() |> 
      group_by(outer.id, species) |> 
      summarize(family.list = list(f)) |> 
      mutate(family.list = map(family.list, unique)) |> 
      group_by(species) |> 
      summarise(intersect = list(intersect_fn(family.list, thresh)))
  }
  else{return("invalid taxonomic level")}
  return(res)
}

# Function: find union based on taxonomic level
# input: dataset containing validation id, species, asv, and taxonomic info
# output: 1 row per species containing intersection based on taxonomic level
union_fn_taxa <- function(data, tax.level){
  
  if(tolower(tax.level) == "c" | tolower(tax.level) == "class"){
    res <- data |> 
      group_by(species, outer.id, d,p,c) |> 
      drop_na() |> 
      unique() |> 
      group_by(outer.id, species) |> 
      summarize(class.list = list(c)) |> 
      mutate(class.list = map(class.list, unique)) |> 
      group_by(species) |> 
      summarise(union = list(union_fn(class.list)))
  }
  else if(tolower(tax.level) == "o" | tolower(tax.level) == "order"){
    res <- data |> 
      group_by(species, outer.id, d,p,c,o) |> 
      drop_na() |> 
      unique() |> 
      group_by(outer.id, species) |> 
      summarize(order.list = list(o)) |> 
      mutate(order.list = map(order.list, unique)) |> 
      group_by(species) |> 
      summarise(union = list(union_fn(order.list)))
  }
  else if(tolower(tax.level) == "f" | tolower(tax.level) == "family"){
    res <- data |> 
      group_by(species, outer.id, d,p,c,o,f) |> 
      drop_na() |> 
      unique() |> 
      group_by(outer.id, species) |> 
      summarize(family.list = list(f)) |> 
      mutate(family.list = map(family.list, unique)) |> 
      group_by(species) |> 
      summarise(union = list(union_fn(family.list)))
  }
  else{return("invalid taxonomic level")}
  return(res)
}

## 18SV9 ##

# merge validation stable sets with asv taxa
validation_stable_sets_18sv9 <- paste(val_dir, '18sv9-ss/validation-stable-sets.rds', sep = '') |>
  read_rds() |>
  mutate(asv = map(asv, unique)) |>
  unnest(asv)

validation_ss_taxa_18sv9 <- validation_stable_sets_18sv9 |>
  left_join(asv_taxa_18sv9, join_by(asv == short.id))

# measures of overlap at asv level
asv_jindex_18sv9 <- paste(val_dir, '18sv9-ss/validation-stable-sets.rds', sep = '') |>
  read_rds() |>
  mutate(asv = map(asv, unique)) |>
  select(species, asv)  |>
  group_by(species) |>
  summarize(intersect = intersect_fn(asv, 0.5) |> list(),
            union = union_fn(asv) |> list()) |>
  mutate(j.index = map2(intersect, union, ~length(.x)/length(.y))) |>
  mutate(n.intersect = map(intersect, length),
         n.union = map(union, length)) |>
  unnest(-c(intersect, union)) |>
  mutate(taxonomic.level = "asv") |>
  select(species, taxonomic.level, intersect, union, n.intersect, n.union, j.index)

# intersections and unions at class, order, and family levels (threshold = 0.5)
class_int_18sv9 <- intersect_fn_taxa(validation_ss_taxa_18sv9, "c", 0.5) |> 
  mutate(taxonomic.level = "class")

order_int_18sv9 <- intersect_fn_taxa(validation_ss_taxa_18sv9, "o", 0.5) |> 
  mutate(taxonomic.level = "order")

family_int_18sv9 <- intersect_fn_taxa(validation_ss_taxa_18sv9, "f", 0.5) |> 
  mutate(taxonomic.level = "family")

class_un_18sv9 <- union_fn_taxa(validation_ss_taxa_18sv9, "c") |> 
  mutate(taxonomic.level = "class")

order_un_18sv9 <- union_fn_taxa(validation_ss_taxa_18sv9, "o") |> 
  mutate(taxonomic.level = "order")

family_un_18sv9 <- union_fn_taxa(validation_ss_taxa_18sv9, "f") |> 
  mutate(taxonomic.level = "family")

# measures of overlap (j index) at class, order, and family level
intersection_18sv9_tax <- rbind(class_int_18sv9, order_int_18sv9, family_int_18sv9)
union_18sv9_tax <- rbind(class_un_18sv9, order_un_18sv9, family_un_18sv9)
tax_jindex_18sv9 <- intersection_18sv9_tax |> 
  inner_join(union_18sv9_tax, join_by(species, taxonomic.level)) |> 
  mutate(n.intersect = sapply(intersect, length),
         n.union = sapply(union, length),
         j.index = n.intersect/n.union) |>
  select(species, taxonomic.level, intersect, union, n.intersect, n.union, j.index)

# combine (retains specific asv lists)
selection_consistency_18sv9_long <- rbind(asv_jindex_18sv9, tax_jindex_18sv9) |>
  mutate(taxonomic.level = factor(taxonomic.level, levels = c('asv', 'family', 'class', 'order'))) |>
  arrange(species, taxonomic.level)

# render as table
selection_consistency_18sv9 <- selection_consistency_18sv9_long |>
  mutate(marker = '18SV9',
         n = n.union,
         j = j.index) |>
  select(species, marker, taxonomic.level, j, n) |>
  pivot_wider(names_from = taxonomic.level, values_from = c(j, n)) |>
  select(species, marker, ends_with('asv'), ends_with('family'), ends_with('class'), ends_with('order'))
  
## 18SV4 ##

# merge validation stable sets with asv taxa
validation_stable_sets_18sv4 <- paste(val_dir, '18sv4-ss/validation-stable-sets.rds', sep = '') |>
  read_rds() |>
  mutate(asv = map(asv, unique)) |>
  unnest(asv)

validation_ss_taxa_18sv4 <- validation_stable_sets_18sv4 |>
  left_join(asv_taxa_18sv4, join_by(asv == short.id))

# measures of overlap at asv level
asv_jindex_18sv4 <- paste(val_dir, '18sv4-ss/validation-stable-sets.rds', sep = '') |>
  read_rds() |>
  mutate(asv = map(asv, unique)) |>
  select(species, asv)  |>
  group_by(species) |>
  summarize(intersect = intersect_fn(asv, 0.5) |> list(),
            union = union_fn(asv) |> list()) |>
  mutate(j.index = map2(intersect, union, ~length(.x)/length(.y))) |>
  mutate(n.intersect = map(intersect, length),
         n.union = map(union, length)) |>
  unnest(-c(intersect, union)) |>
  mutate(taxonomic.level = "asv") |>
  select(species, taxonomic.level, intersect, union, n.intersect, n.union, j.index)

# intersections and unions at class, order, and family levels (threshold = 0.5)
class_int_18sv4 <- intersect_fn_taxa(validation_ss_taxa_18sv4, "c", 0.5) |> 
  mutate(taxonomic.level = "class")

order_int_18sv4 <- intersect_fn_taxa(validation_ss_taxa_18sv4, "o", 0.5) |> 
  mutate(taxonomic.level = "order")

family_int_18sv4 <- intersect_fn_taxa(validation_ss_taxa_18sv4, "f", 0.5) |> 
  mutate(taxonomic.level = "family")

class_un_18sv4 <- union_fn_taxa(validation_ss_taxa_18sv4, "c") |> 
  mutate(taxonomic.level = "class")

order_un_18sv4 <- union_fn_taxa(validation_ss_taxa_18sv4, "o") |> 
  mutate(taxonomic.level = "order")

family_un_18sv4 <- union_fn_taxa(validation_ss_taxa_18sv4, "f") |> 
  mutate(taxonomic.level = "family")

# measures of overlap (j index) at class, order, and family level
intersection_18sv4_tax <- rbind(class_int_18sv4, order_int_18sv4, family_int_18sv4)
union_18sv4_tax <- rbind(class_un_18sv4, order_un_18sv4, family_un_18sv4)
tax_jindex_18sv4 <- intersection_18sv4_tax |> 
  inner_join(union_18sv4_tax, join_by(species, taxonomic.level)) |> 
  mutate(n.intersect = sapply(intersect, length),
         n.union = sapply(union, length),
         j.index = n.intersect/n.union) |>
  select(species, taxonomic.level, intersect, union, n.intersect, n.union, j.index)

# combine (retains specific asv lists)
selection_consistency_18sv4_long <- rbind(asv_jindex_18sv4, tax_jindex_18sv4) |>
  mutate(taxonomic.level = factor(taxonomic.level, levels = c('asv', 'family', 'class', 'order'))) |>
  arrange(species, taxonomic.level)

# render as table
selection_consistency_18sv4 <- selection_consistency_18sv4_long |>
  mutate(marker = '18SV4',
         n = n.union,
         j = j.index) |>
  select(species, marker, taxonomic.level, j, n) |>
  pivot_wider(names_from = taxonomic.level, values_from = c(j, n)) |>
  select(species, marker, ends_with('asv'), ends_with('family'), ends_with('class'), ends_with('order'))

## 16S ##

# merge validation stable sets with asv taxa
validation_stable_sets_16s <- paste(val_dir, '16s-ss/validation-stable-sets.rds', sep = '') |>
  read_rds() |>
  mutate(asv = map(asv, unique)) |>
  unnest(asv)

validation_ss_taxa_16s <- validation_stable_sets_16s |>
  left_join(asv_taxa_16s, join_by(asv == short.id))

# measures of overlap at asv level
asv_jindex_16s <- paste(val_dir, '16s-ss/validation-stable-sets.rds', sep = '') |>
  read_rds() |>
  mutate(asv = map(asv, unique)) |>
  select(species, asv)  |>
  group_by(species) |>
  summarize(intersect = intersect_fn(asv, 0.5) |> list(),
            union = union_fn(asv) |> list()) |>
  mutate(j.index = map2(intersect, union, ~length(.x)/length(.y))) |>
  mutate(n.intersect = map(intersect, length),
         n.union = map(union, length)) |>
  unnest(-c(intersect, union)) |>
  mutate(taxonomic.level = "asv") |>
  select(species, taxonomic.level, intersect, union, n.intersect, n.union, j.index)

# intersections and unions at class, order, and family levels (threshold = 0.5)
class_int_16s <- intersect_fn_taxa(validation_ss_taxa_16s, "c", 0.5) |> 
  mutate(taxonomic.level = "class")

order_int_16s <- intersect_fn_taxa(validation_ss_taxa_16s, "o", 0.5) |> 
  mutate(taxonomic.level = "order")

family_int_16s <- intersect_fn_taxa(validation_ss_taxa_16s, "f", 0.5) |> 
  mutate(taxonomic.level = "family")

class_un_16s <- union_fn_taxa(validation_ss_taxa_16s, "c") |> 
  mutate(taxonomic.level = "class")

order_un_16s <- union_fn_taxa(validation_ss_taxa_16s, "o") |> 
  mutate(taxonomic.level = "order")

family_un_16s <- union_fn_taxa(validation_ss_taxa_16s, "f") |> 
  mutate(taxonomic.level = "family")

# measures of overlap (j index) at class, order, and family level
intersection_16s_tax <- rbind(class_int_16s, order_int_16s, family_int_16s)
union_16s_tax <- rbind(class_un_16s, order_un_16s, family_un_16s)
tax_jindex_16s <- intersection_16s_tax |> 
  inner_join(union_16s_tax, join_by(species, taxonomic.level)) |> 
  mutate(n.intersect = sapply(intersect, length),
         n.union = sapply(union, length),
         j.index = n.intersect/n.union) |>
  select(species, taxonomic.level, intersect, union, n.intersect, n.union, j.index)

# combine (retains specific asv lists)
selection_consistency_16s_long <- rbind(asv_jindex_16s, tax_jindex_16s) |>
  mutate(taxonomic.level = factor(taxonomic.level, levels = c('asv', 'family', 'class', 'order'))) |>
  arrange(species, taxonomic.level)

# render as table
selection_consistency_16s <- selection_consistency_16s_long |>
  mutate(marker = '16S',
         n = n.union,
         j = j.index) |>
  select(species, marker, taxonomic.level, j, n) |>
  pivot_wider(names_from = taxonomic.level, values_from = c(j, n)) |>
  select(species, marker, ends_with('asv'), ends_with('family'), ends_with('class'), ends_with('order'))

## COMBINE ALL ##

selection_consistency <- rbind(selection_consistency_16s, 
                                   selection_consistency_18sv4, 
                                   selection_consistency_18sv9) |>
  mutate(species = factor(species, 
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback')))

selection_consistency

## TABLES 4a-c: COEFFICIENTS ---------------------------------------------------

top_coef_16s <- sel_asv_16s |> 
  drop_na(o) |> 
  mutate(mag = abs(coef.lr),
         marker = '16S') |> 
  group_by(species) |> 
  slice_max(mag, n = 5) |> 
  select(species, marker, coef.lr, coef.ss, d, p, c, o, f, g) |>
  arrange(species, desc(coef.ss)) |>
  mutate(across(where(is.numeric), ~round(.x, 3)))

top_coef_18sv4 <- sel_asv_18sv4 |> 
  drop_na(o) |> 
  mutate(mag = abs(coef.lr),
         marker = '18SV4') |> 
  group_by(species) |> 
  slice_max(mag, n = 5) |> 
  select(species, marker, coef.lr, coef.ss, d, p, c, o, f, g) |>
  arrange(species, desc(coef.ss)) |>
  mutate(across(where(is.numeric), ~round(.x, 3)))

top_coef_18sv9 <- sel_asv_18sv9 |> 
  drop_na(o) |> 
  mutate(mag = abs(coef.lr),
         marker = '18SV9') |> 
  group_by(species) |> 
  slice_max(mag, n = 5) |> 
  select(species, marker, coef.lr, coef.ss, d, p, c, o, f, g) |>
  arrange(species, desc(coef.ss)) |>
  mutate(across(where(is.numeric), ~round(.x, 3)))

## TABLE 5: PREDICTIONS --------------------------------------------------------

# prediction metrics from 18sv9 model
pred_metrics_18sv9 <- paste(stbl_dir, '18sv9-ss/loo-preds.rds', sep = '') |>
  read_rds() |>
  group_by(species) |>
  summarize(cor.ss = cor(pred.ss, obs.ss),
            cor.lr = cor(pred.lr, obs.lr),
            rmspe.ss = (pred.ss - obs.ss)^2 |> mean() |> sqrt(),
            rmspe.lr = (pred.lr - obs.lr)^2 |> mean() |> sqrt()) |>
  mutate(marker = '18SV9')

# prediction metrics from 18sv4 model
pred_metrics_18sv4 <- paste(stbl_dir, '18sv4-ss/loo-preds.rds', sep = '') |> 
  read_rds() |>
  group_by(species) |>
  summarize(cor.ss = cor(pred.ss, obs.ss),
            cor.lr = cor(pred.lr, obs.lr),
            rmspe.ss = (pred.ss - obs.ss)^2 |> mean() |> sqrt(),
            rmspe.lr = (pred.lr - obs.lr)^2 |> mean() |> sqrt()) |>
  mutate(marker = '18SV4')

# prediction metrics from 16s model
pred_metrics_16s <- paste(stbl_dir, '16s-ss/loo-preds.rds', sep = '') |> 
  read_rds() |>
  group_by(species) |>
  summarize(cor.ss = cor(pred.ss, obs.ss),
            cor.lr = cor(pred.lr, obs.lr),
            rmspe.ss = (pred.ss - obs.ss)^2 |> mean() |> sqrt(),
            rmspe.lr = (pred.lr - obs.lr)^2 |> mean() |> sqrt()) |>
  mutate(marker = '16S')

# naive predictions
naive_preds <- paste(model_dir, 'naive-preds.rds', sep = '') |> read_rds()

# combine
pred_metrics <- rbind(pred_metrics_16s, pred_metrics_18sv4, pred_metrics_18sv9) |>
  select(species, marker, ends_with('lr'), ends_with('ss')) |>
  left_join(naive_preds) |>
  mutate(rel.reduction.lag = (rmspe.lag - rmspe.ss)/rmspe.lag,
         rel.reduction.mean = (rmspe.mean - rmspe.ss)/rmspe.mean,
         across(where(is.double), ~round(.x, 3)),
         across(starts_with('rel'), ~.x*100)) |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback'))) |>
  select(species, marker, ends_with('lr'), ends_with('ss'), starts_with('rel'))

## TABLE 6: LITERATURE OVERLAP -------------------------------------------------

# Read in lit review overlap tables (fixed overlap calculation)
class_overlap <- read_rds('rslt/tbl/overlap-table-class.rds') |>
  mutate(species = factor(species,
                          levels = c('blues', 'fins', 'humpbacks'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback')),
         taxonomic.level = 'class')
order_overlap <- read_rds('rslt/tbl/overlap-table-order.rds') |>
  mutate(species = factor(species,
                          levels = c('blues', 'fins', 'humpbacks'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback')),
         taxonomic.level = 'order') 

lit_overlap <- bind_rows(class_overlap, order_overlap) |>
  select(taxonomic.level, gene.seq, species, ends_with('count'), ends_with('overlap'))

## SUPP TABLE 5a-c: COMMONLY SELECTED TAXA -------------------------------------

#16s
tax_intersection_list_16s <- rbind(tax_jindex_16s, asv_jindex_16s) |> 
  select(species, taxonomic.level, intersect) |> 
  arrange(species) |> 
  rename(intersection = intersect) |> 
  mutate(species = factor(species, 
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback'))) |> 
  unnest(intersection)

#18sv4
tax_intersection_list_18sv4 <- rbind(tax_jindex_18sv4, asv_jindex_18sv4) |> 
  select(species, taxonomic.level, intersect) |> 
  arrange(species) |> 
  rename(intersection = intersect) |> 
  mutate(species = factor(species, 
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback'))) |> 
  unnest(intersection)

#18sv9
tax_intersection_list_18sv9 <- rbind(tax_jindex_18sv9, asv_jindex_18sv9) |> 
  select(species, taxonomic.level, intersect) |> 
  arrange(species) |> 
  rename(intersection = intersect) |> 
  mutate(species = factor(species, 
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback'))) |> 
  unnest(intersection)

## SUPP TABLE 6a-b: LIT/NCOG OVERLAP (ORDER) -----------------------------------

## 6a: all direct relationships
doc.rel <- read.csv(here::here("data/whale-edna-relationships.csv"))

## 6b: overlap with ncog data at order level

# FUNCTION: get_ncog_asvs()  ---------------------------------------------------------------------
# function takes input marker and taxonomic level
# function outputs all asvs (of that marker) that appear in the documented relationships file
get_ncog_asvs <- function(marker, taxonomic.level){
  if (!(tolower(marker) %in% c("16s", "18sv4", "18sv9"))){
    print("invalid marker")
  }
  if (tolower(taxonomic.level) == "genus" | tolower(taxonomic.level) == "g"){
    asv_taxa_all |> 
      filter(gene.seq == tolower(marker)) |> 
      inner_join(doc.rel, join_by(p == Phylum,
                                  c == Class,
                                  o == Order, 
                                  f == Family,
                                  g == Genus)) |> 
      select(gene.seq,short.id,
             Whale.species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
  else if (tolower(taxonomic.level) == "family" | tolower(taxonomic.level) == "f"){
    asv_taxa_all |> 
      filter(gene.seq == tolower(marker)) |> 
      inner_join(doc.rel, join_by(p == Phylum,
                                  c == Class,
                                  o == Order, 
                                  f == Family)) |> 
      select(gene.seq,short.id,
             Whale.species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
  else if (tolower(taxonomic.level) == "order" | tolower(taxonomic.level) == "o"){
    asv_taxa_all |> 
      filter(gene.seq == tolower(marker)) |> 
      inner_join(doc.rel, join_by(p == Phylum,
                                  c == Class,
                                  o == Order)) |> 
      select(gene.seq,short.id,
             Whale.species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
  else if (tolower(taxonomic.level) == "class" | tolower(taxonomic.level) == "c"){
    asv_taxa_all |> 
      filter(gene.seq == tolower(marker)) |> 
      inner_join(doc.rel, join_by(p == Phylum,
                                  c == Class),
                 relationship = "many-to-many") |> 
      select(gene.seq,short.id,
             Whale.species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
  else if (tolower(taxonomic.level) == "phylum" | tolower(taxonomic.level) == "p"){
    asv_taxa_all |> 
      filter(gene.seq == tolower(marker)) |> 
      inner_join(doc.rel, join_by(p == Phylum),
                 relationship = "many-to-many") |> 
      select(gene.seq,short.id,
             Whale.species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
  else{print("invalid taxonomic level")}
}

# join candidate asv tables
asv_taxa_16s <- asv_taxa_16s |> 
  mutate(gene.seq = "16s") 

asv_taxa_18sv4 <- asv_taxa_18sv4|> 
  mutate(gene.seq = "18sv4")

asv_taxa_18sv9 <- asv_taxa_18sv9|> 
  mutate(gene.seq = "18sv9")

asv_taxa_all <- rbind(asv_taxa_16s,
                      asv_taxa_18sv4,
                      asv_taxa_18sv9)


# create table
ncog_lit_overlap_order <- rbind(get_ncog_asvs("16s", "o"),
                                get_ncog_asvs("18sv4", "o"),
                                get_ncog_asvs("18sv9", "o")) |>
  rename_with(tolower) |>
  mutate(whale.species = factor(whale.species,
                          levels = c('blue whales', 'fin whales', 'humpback whales'),
                          labels = c('Blue',
                                     'Fin',
                                     'Humpback')))

## SUPP TABLE 7: MODELS FIT TO ASVS IN LIT REVIEW ------------------------------

# These models were fit on candidate asv overlap with documented relationships
# LEVEL OF OVERLAP: ORDER
candidate_model_results <- read_rds('rslt/tbl/candidate-asv-model-res.rds')

# These models were fit on selected asv overlap with documented relatioships
# LEVEL OF OVERLAP: CLASS
ss_model_results <- read_rds('rslt/tbl/stable-set-asv-model-res.rds')

## EXPORT TABLES ---------------------------------------------------------------

# named list of tables
sheets <- list("tbl1a-ednasampling" = ncog_summary_short,
               "tbl1b-visualsampling" = sighting_summary,
               "tbl2-modelfit" = model_fit,
               "tbl3-validation" = selection_consistency,
               "tbl4a-coef16s" = top_coef_16s,
               "tbl4b-coef18sv4" = top_coef_18sv4,
               "tbl4c-coef18sv9" = top_coef_18sv9,
               "tbl5-modelpred" = pred_metrics,
               "stbl2-ednasampling-long" = ncog_summary,
               "stbl3a-candidates16s" = asv_taxa_16s,
               "stbl3b-candidates18sv4" = asv_taxa_18sv4,
               "stbl3c-candidates18sv9" = asv_taxa_18sv9,
               "stbl4a-coef16s" = sel_asv_16s,
               "stbl4b-coef18sv4" = sel_asv_18sv4,
               "stbl4c-coef18sv9" = sel_asv_18sv9)

# export as excel workbook
writexl::write_xlsx(sheets, paste(out_dir, 'tables.xlsx', sep = ''))
