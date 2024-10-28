library(tidyverse)
library(lubridate)
library(patchwork)
data_dir <- 'data/processed/'
model_dir <- 'rslt/models/scaled-sightings/'
stbl_dir <- 'rslt/stability-selection/'
val_dir <- 'rslt/nested-validation/'
dirs <- c('data_dir', 'model_dir', 'stbl_dir', 'val_dir', 'dirs')

tbl_out_dir <- 'rslt/tbl/'
fs::dir_create(tbl_out_dir)
fig_out_dir <- 'rslt/fig/'
fs::dir_create(fig_out_dir)

## ASV TABLES ------------------------------------------------------------------

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

# candidate asvs from 18sv9
paste(data_dir, 'ncog18sv9.RData', sep = '') |> load()
asv_taxa_18sv9 <- asv_taxa 

# candidate asvs from 18sv4
paste(data_dir, 'ncog18sv4.RData', sep = '') |> load()
asv_taxa_18sv4 <- asv_taxa

# candidate asvs from 16s
paste(data_dir, 'ncog16s.RData', sep = '') |> load()
asv_taxa_16s <- asv_taxa

# selected asvs from 18sv9 model
paste(data_dir, 'ncog18sv9.RData', sep = '') |> load()
paste(model_dir, 'fitted-models-18sv9.RData', sep = '') |> load()
sel_asv_18sv9 <- fitted_models |>
  mutate(species = factor(species, 
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whale', 
                                     'Fin whale', 
                                     'Humpback whale')),
         coef.ss = map(coef, ~2^.x)) |>
  select(species, coef, coef.ss, ss) |>
  unnest(everything()) |>
  rename(asv = ss,
         coef.lr = coef) |>
  left_join(asv_taxa, join_by(asv == short.id))



# prediction metrics from 18sv9 model
pred_metrics <- paste(stbl_dir, '18sv9-ss/loo-preds.rds', sep = '') |> 
  read_rds() |>
  group_by(species) |>
  summarize(cor.ss = cor(pred.ss, obs.ss),
            cor.lr = cor(pred.lr, obs.lr),
            rmspe.ss = (pred.ss - obs.ss)^2 |> mean() |> sqrt(),
            rmspe.lr = (pred.lr - obs.lr)^2 |> mean() |> sqrt())

# merge validation stable sets with asv taxa
# why when run 71-73 theres 3 rows per outr.id?
validation_stable_sets_18sv9 <- paste(val_dir, '18sv9-ss/validation-stable-sets.rds', sep = '') |>
  read_rds() |>
  mutate(asv = map(asv, unique)) |> 
  unnest(asv)

# dataframe to pass into intersect/union functions
validation_ss_taxa_18sv9 <- validation_stable_sets_18sv9 |> 
  left_join(asv_taxa_18sv9, join_by(asv == short.id))

# Check counts are the same
validation_ss_taxa_18sv9 |> group_by(outer.id) |> count()

validation_stable_sets_18sv9 |> group_by(outer.id) |> count()

# CODE CHUNK FOR TESTING
# validation_ss_taxa |>
#   group_by(species, outer.id, d,p,c) |>
#   drop_na() |>
#   unique() |>
#   group_by(outer.id, species) |>
#   summarize(class.list = list(c)) |>
#   mutate(class.list = map(class.list, unique)) |>
#   group_by(species) |>
#   summarise(intersect = list(union_fn(class.list)))


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

# test function
intersect_fn_taxa(validation_ss_taxa_18sv9, "c", 0.5)

intersect_fn_taxa(validation_ss_taxa_18sv9, "o", 0.5)

intersect_fn_taxa(validation_ss_taxa_18sv9, "f", 0.5)

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

# test function
union_fn_taxa(validation_ss_taxa_18sv9, "c")

union_fn_taxa(validation_ss_taxa_18sv9, "o")

union_fn_taxa(validation_ss_taxa_18sv9, "f")

# Join all intersections (threshold = 0.5)
class_int_18sv9 <- intersect_fn_taxa(validation_ss_taxa_18sv9, "c", 0.5) |> 
  mutate(taxonomic.level = "class")

order_int_18sv9 <- intersect_fn_taxa(validation_ss_taxa_18sv9, "o", 0.5) |> 
  mutate(taxonomic.level = "order")

family_int_18sv9 <- intersect_fn_taxa(validation_ss_taxa_18sv9, "f", 0.5) |> 
  mutate(taxonomic.level = "family")

# final table with all intersections by species and taxonomic level
intersection_18sv9_tax <- rbind(class_int_18sv9, order_int_18sv9, family_int_18sv9)

# Join all unions
class_un_18sv9 <- union_fn_taxa(validation_ss_taxa_18sv9, "c") |> 
  mutate(taxonomic.level = "class")

order_un_18sv9 <- union_fn_taxa(validation_ss_taxa_18sv9, "o") |> 
  mutate(taxonomic.level = "order")

family_un_18sv9 <- union_fn_taxa(validation_ss_taxa_18sv9, "f") |> 
  mutate(taxonomic.level = "family")

# final table with all unions by species and taxonomic level
union_18sv9_tax <- rbind(class_un_18sv9, order_un_18sv9, family_un_18sv9)


# Final table containing j index for each species at class, order, and family level
j_index_18sv9 <- intersection_18sv9_tax |> 
  inner_join(union_18sv9_tax, join_by(species, taxonomic.level)) |> 
  mutate(n.intersect = sapply(intersect, length),
         n.union = sapply(union, length),
         j.index = n.intersect/n.union)

# consistency of selection procedure at asv level
selection_consistency <- paste(val_dir, '18sv9-ss/validation-stable-sets.rds', sep = '') |>
  read_rds() |>
  mutate(asv = map(asv, unique)) |>
  select(species, asv)  |>
  group_by(species) |>
  summarize(int = intersect_fn(asv, 0.5) |> list(),
            un = union_fn(asv) |> list()) |>
  mutate(j.index = map2(int, un, ~length(.x)/length(.y))) |>
  mutate(int = map(int, length),
         un = map(un, length)) |>
  unnest(everything()) |> 
  rename(n.intersect = int,
         n.union = un) |> 
  mutate(taxonomic.level = "asv",
         intersect = NA,
         union = NA)

j_index_18sv9 <- rbind(j_index_18sv9, selection_consistency)

# join metrics
model_summary_18sv9 <- model_metrics |>
  mutate(marker = '18sv9') |>
  select(species, marker, n.asv, adj.rsq.lr, adj.rsq.ss) |>
  left_join(pred_metrics) |>
  left_join(selection_consistency)

# selected asvs from 18sv4 model
paste(data_dir, 'ncog18sv4.RData', sep = '') |> load()
paste(model_dir, 'fitted-models-18sv4.RData', sep = '') |> load()
sel_asv_18sv4 <- fitted_models |>
  mutate(species = factor(species, 
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whale', 
                                     'Fin whale', 
                                     'Humpback whale')),
         coef.ss = map(coef, ~2^.x)) |>
  select(species, coef, coef.ss, ss) |>
  unnest(everything()) |>
  rename(asv = ss,
         coef.lr = coef) |>
  left_join(asv_taxa, join_by(asv == short.id))

# prediction metrics from 18sv4 model
pred_metrics <- paste(stbl_dir, '18sv4-ss/loo-preds.rds', sep = '') |> 
  read_rds() |>
  group_by(species) |>
  summarize(cor.ss = cor(pred.ss, obs.ss),
            cor.lr = cor(pred.lr, obs.lr),
            rmspe.ss = (pred.ss - obs.ss)^2 |> mean() |> sqrt(),
            rmspe.lr = (pred.lr - obs.lr)^2 |> mean() |> sqrt())

# merge validation stable sets with asv taxa
validation_stable_sets_18sv4 <- paste(val_dir, '18sv4-ss/validation-stable-sets.rds', sep = '') |>
  read_rds() |>
  mutate(asv = map(asv, unique)) |> 
  unnest(asv)

# dataframe to pass into intersect/union functions
validation_ss_taxa_18sv4 <- validation_stable_sets_18sv4 |> 
  left_join(asv_taxa_18sv4, join_by(asv == short.id))

# Check counts are the same
validation_ss_taxa_18sv4 |> group_by(outer.id) |> count()

validation_stable_sets_18sv4 |> group_by(outer.id) |> count()

# test intersect function
intersect_fn_taxa(validation_ss_taxa_18sv4, "c", 0.5)

intersect_fn_taxa(validation_ss_taxa_18sv4, "o", 0.5)

intersect_fn_taxa(validation_ss_taxa_18sv4, "f", 0.5)

# test union function
union_fn_taxa(validation_ss_taxa_18sv4, "c")

union_fn_taxa(validation_ss_taxa_18sv4, "o")

union_fn_taxa(validation_ss_taxa_18sv4, "f")

# Join all intersections (threshold = 0.5)
class_int_18sv4 <- intersect_fn_taxa(validation_ss_taxa_18sv4, "c", 0.5) |> 
  mutate(taxonomic.level = "class")

order_int_18sv4 <- intersect_fn_taxa(validation_ss_taxa_18sv4, "o", 0.5) |> 
  mutate(taxonomic.level = "order")

family_int_18sv4 <- intersect_fn_taxa(validation_ss_taxa_18sv4, "f", 0.5) |> 
  mutate(taxonomic.level = "family")

# final table with all intersections by species and taxonomic level
intersection_18sv4_tax <- rbind(class_int_18sv4, order_int_18sv4, family_int_18sv4)

# Join all unions
class_un_18sv4 <- union_fn_taxa(validation_ss_taxa_18sv4, "c") |> 
  mutate(taxonomic.level = "class")

order_un_18sv4 <- union_fn_taxa(validation_ss_taxa_18sv4, "o") |> 
  mutate(taxonomic.level = "order")

family_un_18sv4 <- union_fn_taxa(validation_ss_taxa_18sv4, "f") |> 
  mutate(taxonomic.level = "family")

# final table with all unions by species and taxonomic level
union_18sv4_tax <- rbind(class_un_18sv4, order_un_18sv4, family_un_18sv4)


# Final table containing j index for each species at class, order, and family level
j_index_18sv4 <- intersection_18sv4_tax |> 
  inner_join(union_18sv4_tax, join_by(species, taxonomic.level)) |> 
  mutate(n.intersect = sapply(intersect, length),
         n.union = sapply(union, length),
         j.index = n.intersect/n.union)

# consistency of selection procedure
selection_consistency <- paste(val_dir, '18sv4-ss/validation-stable-sets.rds', sep = '') |> 
  read_rds() |>
  select(species, asv)  |>
  group_by(species) |>
  summarize(int = intersect_fn(asv, 0.5) |> list(),
            un = union_fn(asv) |> list()) |>
  mutate(j.index = map2(int, un, ~length(.x)/length(.y))) |>
  mutate(int = map(int, length),
         un = map(un, length)) |>
  unnest(everything())  |> 
  rename(n.intersect = int,
         n.union = un) |> 
  mutate(taxonomic.level = "asv",
         intersect = NA,
         union = NA)

j_index_18sv4 <- rbind(j_index_18sv4, selection_consistency)

# join metrics
model_summary_18sv4 <- model_metrics |>
  mutate(marker = '18sv4') |>
  select(species, marker, n.asv, adj.rsq.lr, adj.rsq.ss) |>
  left_join(pred_metrics) |>
  left_join(selection_consistency)

# selected asvs from 18sv4 model
paste(data_dir, 'ncog16s.RData', sep = '') |> load()
paste(model_dir, 'fitted-models-16s.RData', sep = '') |> load()
sel_asv_16s <- fitted_models |>
  mutate(species = factor(species, 
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whale', 
                                     'Fin whale', 
                                     'Humpback whale')),
         coef.ss = map(coef, ~2^.x)) |>
  select(species, coef, coef.ss, ss) |>
  unnest(everything()) |>
  rename(asv = ss,
         coef.lr = coef) |>
  left_join(asv_taxa, join_by(asv == short.id))

# prediction metrics from 16s model
pred_metrics <- paste(stbl_dir, '16s-ss/loo-preds.rds', sep = '') |> 
  read_rds() |>
  group_by(species) |>
  summarize(cor.ss = cor(pred.ss, obs.ss),
            cor.lr = cor(pred.lr, obs.lr),
            rmspe.ss = (pred.ss - obs.ss)^2 |> mean() |> sqrt(),
            rmspe.lr = (pred.lr - obs.lr)^2 |> mean() |> sqrt())


# merge validation stable sets with asv taxa
validation_stable_sets_16s <- paste(val_dir, '16s-ss/validation-stable-sets.rds', sep = '') |>
  read_rds() |>
  mutate(asv = map(asv, unique)) |> 
  unnest(asv)

# dataframe to pass into intersect/union functions
validation_ss_taxa_16s <- validation_stable_sets_16s |> 
  left_join(asv_taxa_16s, join_by(asv == short.id))

# Check counts are the same
validation_ss_taxa_16s |> group_by(outer.id) |> count()

validation_stable_sets_16s |> group_by(outer.id) |> count()

# test intersect function
intersect_fn_taxa(validation_ss_taxa_16s, "c", 0.5)

intersect_fn_taxa(validation_ss_taxa_16s, "o", 0.5)

intersect_fn_taxa(validation_ss_taxa_16s, "f", 0.5)

# test union function
union_fn_taxa(validation_ss_taxa_16s, "c")

union_fn_taxa(validation_ss_taxa_16s, "o")

union_fn_taxa(validation_ss_taxa_16s, "f")

# Join all intersections (threshold = 0.5)
class_int_16s <- intersect_fn_taxa(validation_ss_taxa_16s, "c", 0.5) |> 
  mutate(taxonomic.level = "class")

order_int_16s <- intersect_fn_taxa(validation_ss_taxa_16s, "o", 0.5) |> 
  mutate(taxonomic.level = "order")

family_int_16s <- intersect_fn_taxa(validation_ss_taxa_16s, "f", 0.5) |> 
  mutate(taxonomic.level = "family")

# final table with all intersections by species and taxonomic level
intersection_16s_tax <- rbind(class_int_16s, order_int_16s, family_int_16s)

# Join all unions
class_un_16s <- union_fn_taxa(validation_ss_taxa_16s, "c") |> 
  mutate(taxonomic.level = "class")

order_un_16s <- union_fn_taxa(validation_ss_taxa_16s, "o") |> 
  mutate(taxonomic.level = "order")

family_un_16s <- union_fn_taxa(validation_ss_taxa_16s, "f") |> 
  mutate(taxonomic.level = "family")

# final table with all unions by species and taxonomic level
union_16s_tax <- rbind(class_un_16s, order_un_16s, family_un_16s)


# Final table containing j index for each species at class, order, and family level
j_index_16s <- intersection_16s_tax |> 
  inner_join(union_16s_tax, join_by(species, taxonomic.level)) |> 
  mutate(n.intersect = sapply(intersect, length),
         n.union = sapply(union, length),
         j.index = n.intersect/n.union)

# consistency of selection procedure
selection_consistency <- paste(val_dir, '16s-ss/validation-stable-sets.rds', sep = '') |> 
  read_rds() |>
  select(species, asv)  |>
  group_by(species) |>
  summarize(int = intersect_fn(asv, 0.5) |> list(),
            un = union_fn(asv) |> list()) |>
  mutate(j.index = map2(int, un, ~length(.x)/length(.y))) |>
  mutate(int = map(int, length),
         un = map(un, length)) |>
  unnest(everything())  |> 
  rename(n.intersect = int,
         n.union = un) |> 
  mutate(taxonomic.level = "asv",
         intersect = NA,
         union = NA)

j_index_16s <- rbind(j_index_16s, selection_consistency)

# join metrics
model_summary_16s <- model_metrics |>
  mutate(marker = '16s') |>
  select(species, marker, n.asv, adj.rsq.lr, adj.rsq.ss) |>
  left_join(pred_metrics) |>
  left_join(selection_consistency)

# join all model summaries
model_summaries <- bind_rows(model_summary_16s, model_summary_18sv4, model_summary_18sv9)

# Read in candidate set & stable set asv model results
candidate_model_results <- read_rds('rslt/tbl/candidate-asv-model-res.rds')
ss_model_results <- read_rds('rslt/tbl/stable-set-asv-model-res.rds')

# Read in lit review overlap tables
class_overlap <- read_rds('rslt/tbl/overlap-table-class.rds')
order_overlap <- read_rds('rslt/tbl/overlap-table-order.rds')

# create final j index table
j_index_16s <- j_index_16s |> 
  mutate(gene.seq = "16s")

j_index_18sv4 <- j_index_18sv4 |> 
  mutate(gene.seq = "18sv4")

j_index_18sv9 <- j_index_18sv9 |> 
  mutate(gene.seq = "18sv9")

j_index_final <- rbind(j_index_16s,j_index_18sv4, j_index_18sv9)

j_index_final <- j_index_final |> 
  select(species, taxonomic.level, gene.seq, j.index) |> 
  pivot_wider(names_from = taxonomic.level,
              values_from = j.index) |> 
  arrange(species, gene.seq) |> 
  rename(j.class = class,
         j.order= order,
         j.family = family,
         j.asv = asv)

# write as excel sheets
sheets <- list("18Sv9-candidates" = asv_taxa_18sv9, 
               "18Sv4-candidates" = asv_taxa_18sv4,
               "16S-candidates" = asv_taxa_16s,
               "18Sv9-selected" = sel_asv_18sv9,
               "18Sv4-selected" = sel_asv_18sv4,
               "16S-selected" = sel_asv_16s,
               "model-metrics" = model_summaries,
               "candidate-asv-overlap-results" = candidate_model_results,
               "ss-asv-overlap-results" = ss_model_results,
               "class-overlap-lit-review" = class_overlap,
               "order-overlap-lit-review" = order_overlap,
               "j-index" = j_index_final)
writexl::write_xlsx(sheets, paste(tbl_out_dir, 'summary-tables.xlsx', sep = ''))

rm(list = setdiff(ls(), dirs))

## FIGURE: TIME SERIES ---------------------------------------------------------

# load sighting data
paste(data_dir, 'mm-sightings.RData', sep = '') |> load()

sighting_data <- sightings_raw |>
  mutate(cruise.ym = ym(cruise)) |>
  pivot_longer(c(bp, bm, mn), names_to = 'species', values_to = 'ss') |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales')),
         season = factor(season, levels = c('winter', 'spring', 'summer', 'fall'))) 

seasonal_means <- ss_means |>
  mutate(across(-season, exp)) |>
  rename(bp = log.bp.imp.mean,
         bm = log.bm.imp.mean,
         mn = log.mn.imp.mean) |>
  pivot_longer(-season, names_to = 'species', values_to = 'ss') |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales')),
         season = factor(season, levels = c('winter', 'spring', 'summer', 'fall'))) |>
  arrange(species, season)

p1 <- sighting_data |>
  ggplot(aes(x = cruise.ym, y = ss)) +
  facet_wrap(~species, nrow = 3) +
  geom_path() +
  theme_bw() +
  labs(x = NULL, y = 'Sightings per 1000km', title = NULL) +
  theme(panel.grid = element_line(linewidth = 0.1, color = 'black'),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 11))


p2 <- sightings_raw |>
  mutate(cruise.ym = ym(cruise)) |>
  pivot_longer(c(bp, bm, mn), names_to = 'species', values_to = 'ss') |>
  mutate(species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales')),
         season = factor(season, levels = c('winter', 'spring', 'summer', 'fall'))) |>
  ggplot(aes(x = season, y = ss)) +
  facet_wrap(~species, nrow = 3) +
  geom_jitter(height = 0, width = 0.2, alpha = 0.5) +
  geom_point(color = 'red', shape = 3, data = seasonal_means) +
  geom_path(aes(group = species), color = 'red', linewidth = 0.2, data = seasonal_means) +
  theme_bw() +
  labs(x = NULL, y = NULL, title = NULL) +
  theme(panel.grid = element_line(linewidth = 0.1, color = 'black'),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 11),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p1 + p2 + plot_layout(ncol = 2, widths = c(2, 1))

paste(fig_out_dir, 'fig-timeseries.png', sep = '') |>
  ggsave(width = 5, height = 4, dpi = 400, units = 'in')

## FIGURE: PREDICTIONS ---------------------------------------------------------

# predictions from 16s
pred_pts_16s <- paste(stbl_dir, '16s-ss/loo-preds.rds', sep = '') |>
  read_rds() |>
  mutate(cruise.ym = ym(obs.id),
         species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales')),
         marker = "16S") |>
  arrange(species, cruise.ym)
  
# predictions from 18sv4
pred_pts_18sv4 <- paste(stbl_dir, '18sv4-ss/loo-preds.rds', sep = '') |>
  read_rds() |>
  mutate(cruise.ym = ym(obs.id),
         species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales')),
         marker = "18SV4") |>
  arrange(species, cruise.ym)

# predictions from 18sv9
pred_pts_18sv9 <- paste(stbl_dir, '18sv9-ss/loo-preds.rds', sep = '') |>
  read_rds() |>
  mutate(cruise.ym = ym(obs.id),
         species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales')),
         marker = "18SV9") |>
  arrange(species, cruise.ym)

# combine predictions from each marker
pred_pts <- bind_rows(pred_pts_16s, pred_pts_18sv4, pred_pts_18sv9) |>
  mutate(key = paste('Predicted (', marker, ')', sep = '')) |>
  select(species, cruise.ym, pred.ss, key, pred.ss.qlo, pred.ss.qhi) |>
  rename(value = pred.ss)

# observations
obs_pts <- pred_pts_18sv9 |>
  select(species, cruise.ym, obs.ss) |>
  mutate(key = 'Observed') |>
  rename(value = obs.ss)

# color palette
pal <- c('#000000', RColorBrewer::brewer.pal(name = 'Dark2', n = 3))

# predictions from each marker overlaid on serial observations
p1 <- bind_rows(pred_pts, obs_pts) |>
  ggplot(aes(x = cruise.ym,
             y = value,
             color = key,
             linetype = key)) +
  facet_wrap(~species, nrow = 3, scale = 'fixed') +
  geom_path(linewidth = 0.4) +
  # geom_ribbon(inherit.aes = F,
  #             aes(x = cruise.ym,
  #                 ymin = pred.ss.qlo,
  #                 ymax = pred.ss.qhi,
  #                 fill = key),
  #             alpha = 0.2) +
  scale_x_date(date_breaks = '1 years', labels = ~year(.x)) +
  scale_color_manual(values = pal) +
  # scale_y_sqrt() +
  theme_bw() +
  # geom_blank(data = pivot_longer(fit_pts_df, c(Observed, Fitted)),
  #            aes(linetype = NULL)) +
  guides(linetype = guide_legend(title = NULL, position = 'top', nrow = 2),
         color = guide_legend(title = NULL, position = 'top', nrow = 2),
         fill = guide_none()) +
  labs(x = 'Year', y = 'Sightings per 1000km', title = NULL) +
  theme(panel.grid = element_line(linewidth = 0.1, color = 'black'),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        text = element_text(size = 11))

# predictions vs observations
p2 <- bind_rows(pred_pts_16s, pred_pts_18sv4, pred_pts_18sv9) |>
  mutate(species = str_remove_all(species, ' whales')) |>
  ggplot(aes(x = obs.ss, y = pred.ss)) +
  facet_grid(species~marker) +
  geom_point(size = 0.8) +
  geom_linerange(aes(ymin = pred.ss.qlo, ymax = pred.ss.qhi),
                 linewidth = 0.4) +
  scale_x_log10() +
  scale_y_log10() +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.2) +
  theme_bw() +
  theme(panel.grid = element_line(linewidth = 0.1, color = 'black'),
        panel.grid.minor = element_blank(),
        text = element_text(size = 11)) +
  labs(x = 'Observed value', y = 'Predicted value')

# correlations between predictions and observations (for labeling)
pred_cors <- bind_rows(pred_pts_16s, pred_pts_18sv4, pred_pts_18sv9) |>
  mutate(species = str_remove_all(species, ' whales')) |>
  group_by(species, marker) |>
  summarize(r = cor(obs.ss, pred.ss),
            obs.ss = 30,
            pred.ss = 0.05,
            .groups = 'drop')

# function to retain trailing zeroes
round_z <- function(x, n){sprintf("%.3f", round(x,n))}
# round_z(0.910, 3)

# add correlation annotations to pred v. obs plots
p2_ann <- p2 + geom_label(data = pred_cors,
                          aes(label = paste('r =', round_z(r, 3))),
                          size = 8,
                          size.unit = 'pt',
                          hjust = 1)

# arrange panel layout
fig_predictions <- p1 + p2_ann +
  plot_layout(nrow = 1, ncol = 2, widths = c(1, 1))

# export
ggsave(fig_predictions, filename = paste(out_dir, 'fig-predictions.png', sep = ''),
       width = 6.5, height = 4, units = 'in', dpi = 400)


## FIGURE: MODEL DIAGNOSTICS ---------------------------------------------------

# predictions from 18sv9
load('rslt/models/scaled-sightings/fitted-models-18sv9-ss.RData')
fit_pts_18sv9 <- fit_df |>
  mutate(cruise.ym = ym(cruise),
         year = year(cruise.ym),
         quarter = factor(season,
                          levels = c('winter', 'spring', 'summer', 'fall')) |>
           as.numeric(),
         cruise.yq = paste(year, quarter, sep = '-') |> yq(),
         species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales'))) |>
  select(cruise, cruise.ym, cruise.yq, species, y, fitted, lr.resid) |>
  arrange(species, cruise.yq) |>
  rename(lr.fit = fitted, lr.obs = y) |>
  mutate(marker = '18SV9')

# predictions from 18sv4
load('rslt/models/scaled-sightings/fitted-models-18sv4-ss.RData')
fit_pts_18sv4 <- fit_df |>
  mutate(cruise.ym = ym(cruise),
         year = year(cruise.ym),
         quarter = factor(season,
                          levels = c('winter', 'spring', 'summer', 'fall')) |>
           as.numeric(),
         cruise.yq = paste(year, quarter, sep = '-') |> yq(),
         species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales'))) |>
  select(cruise, cruise.ym, cruise.yq, species, y, fitted, lr.resid) |>
  arrange(species, cruise.yq) |>
  rename(lr.fit = fitted, lr.obs = y) |>
  mutate(marker = '18SV4')

# predictions from 16s
load('rslt/models/scaled-sightings/increased-res-fitted-models-16s-ss.RData')
fit_pts_16s <- fit_df |>
  mutate(cruise.ym = ym(cruise),
         year = year(cruise.ym),
         quarter = factor(season,
                          levels = c('winter', 'spring', 'summer', 'fall')) |>
           as.numeric(),
         cruise.yq = paste(year, quarter, sep = '-') |> yq(),
         species = factor(species,
                          levels = c('bm', 'bp', 'mn'),
                          labels = c('Blue whales', 'Fin whales', 'Humpback whales'))) |>
  select(cruise, cruise.ym, cruise.yq, species, y, fitted, lr.resid) |>
  arrange(species, cruise.yq) |>
  rename(lr.fit = fitted, lr.obs = y) |>
  mutate(marker = '16S')

# merge above outputs
fit_pts <- bind_rows(fit_pts_16s, fit_pts_18sv4, fit_pts_18sv9)

# for rescaling axes without <scales = 'free_x'>
zscore <- function(x){(x - mean(x))/sd(x)}

# annotations for adding adjusted rsq
ann_data <- fit_pts |>
  group_by(species, marker) |>
  mutate(across(c(lr.fit, lr.obs), zscore)) |>
  summarize(lr.obs = max(lr.obs),
            lr.fit = min(lr.fit)) |>
  left_join(model_summaries, join_by(species, marker)) |>
  ungroup() |>
  mutate(species = fct_relabel(species, ~str_remove_all(.x, ' whales')),
         lr.fit = min(lr.fit) + 0.5,
         lr.obs = max(lr.obs))

# observed vs fitted
obs_fit <- fit_pts |>
  mutate(species = fct_relabel(species, ~str_remove_all(.x, ' whales'))) |>
  group_by(species, marker) |>
  mutate(across(c(lr.fit, lr.obs), zscore)) |>
  ggplot(aes(x = lr.obs, y = lr.fit)) +
  facet_grid(species ~ marker) +
  # facet_wrap(species ~ marker, scales = 'free') +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.4) +
  # geom_smooth(span = 1.5) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = 'black', linewidth = 0.1),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = 'Observed values', y = 'Fitted values', title = 'A. Model fit') +
  geom_label(data = ann_data,
             aes(label = paste("R^'2' == ", round(adj.rsq.ss, 2))),
             parse = T,
             size = 8,
             size.unit = 'pt',
             hjust = 1)

obs_fit

# plot residuals vs fit for each model
resid_fit <- bind_rows(fit_pts_16s, fit_pts_18sv4, fit_pts_18sv9) |>
  mutate(species = fct_relabel(species, ~str_remove_all(.x, ' whales'))) |>
  group_by(species, marker) |>
  mutate(across(c(lr.fit, lr.resid), zscore)) |>
  ggplot(aes(x = lr.fit, y = lr.resid)) +
  facet_grid(species ~ marker) +
  # facet_wrap(species ~ marker, scales = 'free') +
  geom_point() +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  # geom_smooth(span = 1.5) +
  theme_bw() +
  theme(panel.grid.major = element_line(color = 'black', linewidth = 0.1),
        panel.grid.minor = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(x = 'Fitted values', y = 'Residuals', title = 'B. Residual diagnostics')

resid_fit

# function for residual autocorrelation
pacf_fn <- function(x){
  pacf_out <- pacf(x, plot = F, lag.max = 10)
  out <- bind_cols(pacf = c(1, pacf_out$acf[, 1, 1]),
                   lag = c(0, pacf_out$lag[, , 1]),
                   se = 2/sqrt(pacf_out$n.used))
  return(out)
}

# plot pacf for each model
resid_pacf <- bind_rows(fit_pts_16s, fit_pts_18sv4, fit_pts_18sv9) |>
  select(species, marker, lr.resid) |>
  mutate(species = fct_relabel(species, ~str_remove_all(.x, ' whales'))) |>
  nest(resids = lr.resid, .by = c(species, marker)) |>
  mutate(pacf = map(resids, pacf_fn)) |>
  unnest(pacf) |>
  ggplot(aes(x = lag)) +
  facet_grid(species ~ marker) +
  scale_y_continuous(limits = c(-1, 1), n.breaks = 4) +
  scale_x_continuous(breaks = seq(0, 10, by = 2)) +
  geom_linerange(aes(ymin = 0, ymax = pacf),
                 linewidth = 0.4) +
  geom_ribbon(aes(ymin = -se, ymax = se),
              fill = 'blue',
              alpha = 0.2) +
  theme_bw() +
  theme(panel.grid = element_line(color = 'black', linewidth = 0.1),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  labs(x = 'Lag', title = 'C. Residual PACF')

resid_pacf

# panel layout
fig_resid <- obs_fit + resid_fit + resid_pacf + plot_layout(nrow = 1, widths = c(1, 1, 1))
ggsave(fig_resid, filename = 'rslt/_draft/plots/resids-ss-diagnostics.png',
       width = 9, height = 4, units = 'in', dpi = 400)



