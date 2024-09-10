library(tidyverse)
library(openxlsx)
library(readxl)
library(here)

# candidate taxa overlap with documented relationships --------------------------

# load in documented relationships table
doc.rel <- read.csv(here::here("data/whale-edna-relationships.csv"))

# load asv tables 
table_names <- getSheetNames("rslt/tbl/asv-tables.xlsx")

all_sheets <- lapply(table_names, function(sheet) {
  read_excel("rslt/tbl/asv-tables.xlsx", sheet = sheet)
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
  mutate(model.species = paste(tolower(model.species), "s", sep = ""))


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
                                  c == Class)) |> 
      select(gene.seq,short.id,
             Whale.species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
  else if (tolower(taxonomic.level) == "phylum" | tolower(taxonomic.level) == "p"){
    asv_taxa_all |> 
      filter(gene.seq == tolower(marker)) |> 
      inner_join(doc.rel, join_by(p == Phylum)) |> 
      select(gene.seq,short.id,
             Whale.species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
  else{print("invalid taxonomic level")}
}


# 16s testing 
# genus level - 7 matches 
get_ncog_asvs("16s", "g")

# family level - 31 matches
get_ncog_asvs("16s", "f")

# order level - 47 matches
get_ncog_asvs("16s", "order")



# 18sv4 testing
# genus level - 0 matches 
get_ncog_asvs("18sV4", "g")

# family level - 0 matches
get_ncog_asvs("18sv4", "f")

# order level - 2 matches
get_ncog_asvs("18sv4", "order")



# 18sv9 testing
# genus level - 2 matches 
get_ncog_asvs("18sV9", "g")

# family level - 4 matches
get_ncog_asvs("18sv9", "f")

# order level - 12 matches
get_ncog_asvs("18sv9", "order")


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
             model.species == tolower(whale)) |> 
      inner_join(doc.rel, join_by(model.species == Whale.species,
                                  p == Phylum,
                                  c == Class,
                                  o == Order, 
                                  f == Family,
                                  g == Genus)) |> 
      select(gene.seq,asv.id,
             model.species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
  else if (tolower(taxonomic.level) == "family" | tolower(taxonomic.level) == "f"){
    ss_asvs_all |> 
      filter(gene.seq == tolower(marker),
             model.species == tolower(whale)) |> 
      inner_join(doc.rel, join_by(model.species == Whale.species,
                                  p == Phylum,
                                  c == Class,
                                  o == Order,
                                  f == Family)) |> 
      select(gene.seq,asv.id,
             model.species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
  else if (tolower(taxonomic.level) == "order" | tolower(taxonomic.level) == "o"){
  ss_asvs_all |> 
    filter(gene.seq == tolower(marker),
           model.species == tolower(whale)) |> 
    inner_join(doc.rel, join_by(model.species == Whale.species,
                                p == Phylum,
                                c == Class,
                                o == Order)) |> 
    select(gene.seq,asv.id,
           model.species,
           Connection,
           p,c,o,f,g) |> 
    distinct()
  }
  else if (tolower(taxonomic.level) == "class" | tolower(taxonomic.level) == "c"){
    ss_asvs_all |> 
      filter(gene.seq == tolower(marker),
             model.species == tolower(whale)) |> 
      inner_join(doc.rel, join_by(model.species == Whale.species,
                                  p == Phylum,
                                  c == Class)) |> 
      select(gene.seq,asv.id,
             model.species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
  else if (tolower(taxonomic.level) == "phylum" | tolower(taxonomic.level) == "p"){
    ss_asvs_all |> 
      filter(gene.seq == tolower(marker),
             model.species == tolower(whale)) |> 
      inner_join(doc.rel, join_by(model.species == Whale.species,
                                  p == Phylum)) |> 
      select(gene.seq,asv.id,
             model.species,
             Connection,
             p,c,o,f,g) |> 
      distinct()
  }
  
  
}

# only two at order level (across all markers/species)
get_model_asvs("16s", "humpback whales", "o")
get_model_asvs("18sv9", "humpback whales", "o")


ss_asvs_all |> 
  inner_join(doc.rel, join_by(model.species == Whale.species,
                                           p == Phylum,
                                           c == Class,
                                           o == Order)) |> 
  select(gene.seq,asv.id,
         model.species,
         Connection,
         p,c,o,f,g) |> 
  distinct()


# None at family level
ss_asvs_all |> 
  inner_join(doc.rel, join_by(model.species == Whale.species,
                              p == Phylum,
                              c == Class,
                              o == Order, 
                              f == Family)) |> 
  select(gene.seq,asv.id,
         model.species,
         Connection,
         p,c,o,f,g) |> 
  distinct()


# things get more interesting at class level

ss_asvs_all |> 
  inner_join(doc.rel, join_by(model.species == Whale.species,
                              p == Phylum,
                              c == Class)) |> 
  select(gene.seq,asv.id,
         model.species,
         Connection,
         p,c,o,f,g) |> 
  distinct()


# this is all matches at class level
get_model_asvs("16s", "blue whales",  "c")
get_model_asvs("16s", "humpback whales", "c")
get_model_asvs("18sv9", "humpback whales",  "c")

# ------------------------------------------------------------------------------
# Proportions Exploration
# ------------------------------------------------------------------------------

# Define possible markers and whales
markers <- c("16s", "18sv4", "18sv9")
whales <- c("blue whales", "fin whales", "humpback whales")

# OBJECTIVE 3 ------------------------------------------------------------------
# 3) proportion of lit review ASVs that could have been selected that in fact were selected 
#    (i.e., the number in (2) as a proportion of the number in (1)), for each model
# OBJECTIVE 3 ------------------------------------------------------------------

# Singular example for 16s humpbacks
# note: 1 model/lit review match over 47 ncog/lit review matches = 0.0212766
nrow(get_model_asvs("16s", "humpback whales", "p"))/nrow(get_ncog_asvs("16s", "p"))


# Now, summarize proportion for all marker/whale combinations

# Function to calculate proportion of possible lit review ASVs selection for each model
prop_asv_sel <- function(marker, whale, taxonomic.level = "order") {
  model_asvs <- get_model_asvs(marker, 
                               whale,
                               substring(taxonomic.level[1], 1, 1))
  ncog_asvs <- get_ncog_asvs(marker, 
                             taxonomic.level)
  ratio <- nrow(model_asvs) / nrow(ncog_asvs)
  return(ratio)
}

# Data frame with proportion for each marker/whale combo
prop_selected_df <- data.frame(
  Marker = rep(markers, each = length(whales)),
  Whale = rep(whales, length(markers)),
  # Use mapply to apply the function over all marker/whale combinations
  Proportion = mapply(prop_asv_sel, rep(markers, each = length(whales)), rep(whales, length(markers)))
)

# Show df (note: only 16s humpback and 189v9 humpback have a non-zero prop)
prop_selected_df

# OBJECTIVE 4 ------------------------------------------------------------------
# 4) proportion of ASVs selected in each model that were identified in the lit review
# OBJECTIVE 4 ------------------------------------------------------------------

# Singular example for 18sv9 humpbacks
# num ASVs selected in each model (marker/whale combo) included in lit review table
num_asvs_match <- get_model_asvs("18sv9", "humpback whales", "ss", "o") |> 
  nrow()

# num ASVs selected for the model (marker/whale combo)
num_asvs_ss <- ss_asvs_all |> 
  filter(gene.seq == "18sv9",
         model.species == "humpback whales") |> 
  nrow()

# proportion of ASVs selected in 18sv9 humpback model IDed in lit review
num_asvs_match/num_asvs_ss


# Now, summarize proportion for all marker/whale combinations
# Function to calculate proportion of selected ASVs that were IDed in lit revew
prop_sel_lit <- function(marker, whale, model = "ss", taxonomic.level = "order",
                         known = "n") {
  if (known == "n"){
    selected_asvs_df = ss_asvs_all
  }  
  else {
    selected_asvs_df = ss_asvs_all_known
  }
  num_asvs_match <- get_model_asvs(marker, 
                                   whale, 
                                   model, 
                                   substring(taxonomic.level[1], 1, 1)) |> 
    nrow()
  
  num_asvs_ss <- selected_asvs_df |> 
    filter(gene.seq == marker,
           model.species == whale) |> 
    nrow()
  
  ratio <- num_asvs_match/num_asvs_ss
  return(ratio)
}

# Data frame with proportion for each marker/whale combo
prop_asvs_ided_df <- data.frame(
  Marker = rep(markers, each = length(whales)),
  Whale = rep(whales, length(markers)),
  # Use mapply to apply the function over all marker/whale combinations
  Proportion_IDed = mapply(prop_sel_lit, rep(markers, each = length(whales)), rep(whales, length(markers)))
)

# Show df (note: only 16s humpback and 189v9 humpback have a non-zero prop)
prop_asvs_ided_df

# OBJECTIVE 5 ------------------------------------------------------------------
# 5) proportion of ASVs selected in each model that were not identified in the lit review
# OBJECTIVE 5 ------------------------------------------------------------------

prop_asvs_ided_df["Proportion_not_IDed"] = 1 - prop_asvs_ided_df["Proportion_IDed"]

prop_asvs_ided_df

# OBJECTIVE 6 ------------------------------------------------------------------
# 6) 4-5, but after removing unidentified ASVs
#    (i.e., missing taxonomic classification below domain level)
# OBJECTIVE 6 ------------------------------------------------------------------

# remove unidentified ASVs from selected model ASVs
ss_asvs_all_known <- ss_asvs_all |> 
  filter(!is.na(p), p != "uncultured") 

# Data frame with proportion for each marker/whale combo
prop_asvs_ided_df <- data.frame(
  Marker = rep(markers, each = length(whales)),
  Whale = rep(whales, length(markers)),
  # Use mapply to apply the function over all marker/whale combinations
  Proportion_IDed = mapply(prop_sel_lit, 
                           rep(markers, each = length(whales)), 
                           rep(whales, length(markers)),
                           MoreArgs = list(model = "ss", 
                                           taxonomic.level = "order", 
                                           known = "y"))
)

# Show df (note: only 16s humpback and 189v9 humpback have a non-zero prop)
prop_asvs_ided_df

prop_asvs_ided_df["Proportion_not_IDed"] = 1 - prop_asvs_ided_df["Proportion_IDed"]

prop_asvs_ided_df
  
# note:  prop_sel_lit currently works right if num_asvs_match <- get_model_asvs()
# contains known ASVs, so need to filter out if the ASVs are unknown

