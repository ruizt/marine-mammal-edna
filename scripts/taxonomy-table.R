# taxonomy table draft

library(tidyverse)
library(openxlsx)
library(readxl)
library(here)

# candidate taxa overlap with documented relationships --------------------------

# directories
model_dir <- 'rslt/models/scaled-sightings/'
data_dir <- 'data/processed/'


# load in documented relationships table
doc.rel <- read.csv(here::here("data/whale-edna-relationships.csv"))

doc.rel <- doc.rel |> 
  mutate(Whale.species = str_replace(Whale.species, " whales", "s"))


# candidate asvs
# candidate asvs from 18sv9
paste(data_dir, 'ncog18sv9.RData', sep = '') |> load()
asv_taxa_18sv9 <- asv_taxa 

# candidate asvs from 18sv4
paste(data_dir, 'ncog18sv4.RData', sep = '') |> load()
asv_taxa_18sv4 <- asv_taxa

# candidate asvs from 16s
paste(data_dir, 'ncog16s.RData', sep = '') |> load()
asv_taxa_16s <- asv_taxa


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

asv_taxa_all


# 18sv9
paste(model_dir, 'fitted-models-18sv9.RData', sep = '') |> load()
# selected asvs
ss_asvs_18sv9 <- fitted_models |>
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
  select(-feature.id) |>
  arrange(species, desc(coef.lr))

# 18sv4
paste(model_dir, 'fitted-models-18sv4.RData', sep = '') |> load()

# selected asvs
ss_asvs_18sv4 <- fitted_models |>
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
  select(-feature.id) |>
  arrange(species, desc(coef.lr))

# 16s
paste(model_dir, 'fitted-models-16s.RData', sep = '') |> load()

# selected asvs
ss_asvs_16s <- fitted_models |>
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
  select(-feature.id) |>
  arrange(species, desc(coef.lr))


# join selected asv tables
ss_asvs_16s <- ss_asvs_16s |> 
  mutate(gene.seq = "16s")

ss_asvs_18sv4 <- ss_asvs_18sv4 |> 
  mutate(gene.seq = "18sv4")

ss_asvs_18sv9 <- ss_asvs_18sv9 |> 
  mutate(gene.seq = "18sv9")

ss_asvs_all <- rbind(ss_asvs_16s,
                     ss_asvs_18sv4,
                     ss_asvs_18sv9)

ss_asvs_all <- ss_asvs_all |> 
  mutate(species = paste(tolower(species), "s", sep = ""))

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

# FUNCTION: get_model_asvs()------------------------------------------------------------------------

get_model_asvs <- function(marker, whale, taxonomic.level){
  if (!(tolower(marker) %in% c("16s", "18sv4", "18sv9"))){
    print("invalid marker")
  }
  if (!(tolower(whale) %in% c("blues", "fins", "humpbacks"))){
    print("invalid species")
  }
  
  if (tolower(taxonomic.level) == "genus" | tolower(taxonomic.level) == "g"){
    ss_asvs_all |> 
      filter(gene.seq == tolower(marker),
             species == tolower(whale)) |> 
      inner_join(doc.rel, join_by(species == Whale.species,
                                  p == Phylum,
                                  c == Class,
                                  o == Order, 
                                  f == Family,
                                  g == Genus)) |> 
      select(gene.seq,asv,
             species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
  else if (tolower(taxonomic.level) == "family" | tolower(taxonomic.level) == "f"){
    ss_asvs_all |> 
      filter(gene.seq == tolower(marker),
             species == tolower(whale)) |> 
      inner_join(doc.rel, join_by(species == Whale.species,
                                  p == Phylum,
                                  c == Class,
                                  o == Order,
                                  f == Family)) |> 
      select(gene.seq,asv,
             species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
  else if (tolower(taxonomic.level) == "order" | tolower(taxonomic.level) == "o"){
    ss_asvs_all |> 
      filter(gene.seq == tolower(marker),
             species == tolower(whale)) |> 
      inner_join(doc.rel, join_by(species == Whale.species,
                                  p == Phylum,
                                  c == Class,
                                  o == Order)) |> 
      select(gene.seq,asv,
             species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
  else if (tolower(taxonomic.level) == "class" | tolower(taxonomic.level) == "c"){
    ss_asvs_all |> 
      filter(gene.seq == tolower(marker),
             species == tolower(whale)) |> 
      inner_join(doc.rel, join_by(species == Whale.species,
                                  p == Phylum,
                                  c == Class)) |> 
      select(gene.seq,asv,
             species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
  else if (tolower(taxonomic.level) == "phylum" | tolower(taxonomic.level) == "p"){
    ss_asvs_all |> 
      filter(gene.seq == tolower(marker),
             species == tolower(whale)) |> 
      inner_join(doc.rel, join_by(species == Whale.species,
                                  p == Phylum)) |> 
      select(gene.seq,asv,
             species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
}

# ------------------------------------------------------------------------------
# Function to generate counts and overlap proportion summary by taxa level
# - Input: Taxonomic level
# - Output: A summary table with 7 columns
#     1) Gene Sequence (16s, 18sv4, 18sv9)
#     2) Whale Species (blue, fin, humpback)
#     3) Count of unique occurrences of the selected taxa level in the model
#     4) Count of unique occurrences of the selected taxa level in the lit. review ***CHECK/REDO***
#     5) Overlap of unique taxa occurrences between models and lit review (n) ***CHECK/REDO***
#     6) Proportion of taxa overlap in lit. review
#     7) Proportion of taxa overlap in selected model
# ------------------------------------------------------------------------------

get_summary_table <- function(taxonomic.level = "o"){
  
  # Check if the column exists in the dataframe
  if (!(taxonomic.level %in% names(ss_asvs_all))) {
    stop("Error: Valid taxonomic levels are 'd', 'p', 'c', 'o', 'f', 'g'.")
  }
  
  taxonomic.level = as.name(taxonomic.level)
  
  # 1) lit count (n) - need to split by whale species
  lit_count <- rbind(get_ncog_asvs("16s", taxonomic.level), 
                     get_ncog_asvs("18sv4", taxonomic.level), 
                     get_ncog_asvs("18sv9", taxonomic.level)) |> 
    group_by(Whale.species, gene.seq) |> 
    summarise(count = n_distinct(!!taxonomic.level),
              .groups = "keep")
  
  # 2) model count (n)
  model_count <- ss_asvs_all |> 
    group_by(gene.seq,
             species) |> 
    summarise(count = n_distinct(!!taxonomic.level),
              .groups = "keep")
  
  # Join the 2 tables of model and lit. review count
  merged_df <- model_count %>%
    left_join(lit_count, by = c("gene.seq" = "gene.seq", "species" = "Whale.species"))
  
  # 3) Overlap
  
  # distinct ncogs at specified taxonomic level
  ncogs <- rbind(get_ncog_asvs("16s", taxonomic.level), 
        get_ncog_asvs("18sv4", taxonomic.level), 
        get_ncog_asvs("18sv9", taxonomic.level)) |> 
    group_by(Whale.species, gene.seq) |> 
    distinct({{taxonomic.level}})
  
  # distinct models at specified taxonomic level
  models <- rbind(get_model_asvs("16s", "humpbacks", taxonomic.level),
                  get_model_asvs("16s", "blues", taxonomic.level),
                  get_model_asvs("16s", "fins", taxonomic.level),
                  get_model_asvs("18sv4", "humpbacks", taxonomic.level),
                  get_model_asvs("18sv4", "blues", taxonomic.level),
                  get_model_asvs("18sv4", "fins", taxonomic.level),
                  get_model_asvs("18sv9", "humpbacks", taxonomic.level),
                  get_model_asvs("18sv9", "blues", taxonomic.level),
                  get_model_asvs("18sv9", "fins", taxonomic.level)) |> 
    group_by(species,gene.seq) |> 
    distinct({{taxonomic.level}})
  
  # merge to find overlap
  overlap <- ncogs |> 
    inner_join(models, join_by(Whale.species == species, gene.seq, {{taxonomic.level}})) |> 
    count()
  
  overlap_prop_table <- merged_df |> 
    # rename count columns
    rename(model.count = count.x, lit.count = count.y) |> 
    mutate(
      # replace NA's with 0's
      across(c(model.count, lit.count), ~ replace_na(., 0)),
      
    )
  
  overlap_prop_table <- overlap_prop_table |> 
    left_join(overlap, join_by(species == Whale.species, gene.seq)) |> 
    rename(overlap = n) |> 
    mutate(
      overlap = ifelse(is.na(), 0, overlap),
      # 4) [Taxa overlap (n)] / [Lit count (n)]
      prop.lit.overlap = overlap/lit.count,
      
      # 5) [Taxa overlap (n)] / [Model count (n)]
      prop.model.overlap = overlap/model.count)
  
  overlap_prop_table
  
}

# MANUAL TESTING

# model count

get_ncog_asvs("16s", "o") |> 
  group_by(Whale.species) |> 
  summarize(count = n_distinct(o))

get_ncog_asvs("16s", "o") |> 
  group_by(Whale.species) |> 
 distinct(o)

# taxa
test1 <- rbind(get_ncog_asvs("16s", "o"), 
      get_ncog_asvs("18sv4", "o"), 
      get_ncog_asvs("18sv9", "o")) |> 
  group_by(Whale.species, gene.seq) |> 
  distinct(o)

# overlap
test2 <- rbind(get_model_asvs("16s", "humpbacks", "o"),
      get_model_asvs("16s", "blues", "o"),
      get_model_asvs("16s", "fins", "o"),
      get_model_asvs("18sv4", "humpbacks", "o"),
      get_model_asvs("18sv4", "blues", "o"),
      get_model_asvs("18sv4", "fins", "o"),
      get_model_asvs("18sv9", "humpbacks", "o"),
      get_model_asvs("18sv9", "blues", "o"),
      get_model_asvs("18sv9", "fins", "o")) |> 
  group_by(species, gene.seq) |> 
  distinct(o)

test3 <- test1 |> 
  inner_join(test2, join_by(Whale.species == species, gene.seq, o)) |> 
  count()


# Summary table for Order
order_table <- get_summary_table()

# Summary table for Class
class_table <- get_summary_table(taxonomic.level = "c")

order_table
class_table

# Save tables for later

saveRDS(order_table, 'rslt/tbl/overlap-table-order.rds')

saveRDS(class_table, 'rslt/tbl/overlap-table-class.rds')
