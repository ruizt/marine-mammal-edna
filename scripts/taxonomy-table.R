# taxonomy table draft

library(tidyverse)
library(openxlsx)
library(readxl)
library(here)

# candidate taxa overlap with documented relationships --------------------------

# load in documented relationships table
doc.rel <- read.csv(here::here("data/whale-edna-relationships.csv"))

# load asv tables 
table_names <- getSheetNames("rslt/tbl/summary-tables.xlsx")
all_sheets <- lapply(table_names, function(sheet) {
  read_excel("rslt/tbl/summary-tables.xlsx", sheet = sheet)
})

# candidate asvs
asv_taxa_18sv9 <- all_sheets[[1]] 
asv_taxa_18sv4 <- all_sheets[[2]] 
asv_taxa_16s <- all_sheets[[3]] 

# selected asvs
ss_asvs_18sv9 <- all_sheets[[4]] 
ss_asvs_18sv4 <- all_sheets[[5]] 
ss_asvs_16s <- all_sheets[[6]] 

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
# note: fit parameter is not used as of now

get_model_asvs <- function(marker, whale, taxonomic.level){
  if (!(tolower(marker) %in% c("16s", "18sv4", "18sv9"))){
    print("invalid marker")
  }
  if (!(tolower(whale) %in% c("blue whales", "fin whales", "humpback whales"))){
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
#     4) Count of unique occurrences of the selected taxa level in the lit. review
#     5) Overlap of unique taxa occurrences between models and lit review (n)
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
  
  overlap_prop_table <- merged_df |> 
    # rename count columns
    rename(model.count = count.x, lit.count = count.y) |> 
    mutate(
      # replace NA's with 0's
      across(c(model.count, lit.count), ~ replace_na(., 0)),
      
      # 3) Taxa overlap (n)
      overlap = pmap_int(list(gene.seq, species),
                         ~ nrow(get_model_asvs(..1, ..2, "o"))),
      
      # 4) [Taxa overlap (n)] / [Lit count (n)]
      prop.lit.overlap = overlap/lit.count,
      
      # 5) [Taxa overlap (n)] / [Model count (n)]
      prop.model.overlap = overlap/model.count
    )
  
  overlap_prop_table
  
}

# Summary table for Order
order_table <- get_summary_table()

# Summary table for Class
class_table <- get_summary_table(taxonomic.level = "c")

order_table
class_table

# Save tables for later

saveRDS(order_table, 'rslt/tbl/overlap-table-order.rds')

saveRDS(class_table, 'rslt/tbl/overlap-table-class.rds')