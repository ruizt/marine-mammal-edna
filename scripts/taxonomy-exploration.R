library(tidyverse)
library(openxlsx)
library(readxl)


# candidate taxa overlap with documented relationships --------------------------

# load in documented relationships table
doc.rel <- read.csv("/Users/nicholaspatrick/Desktop/eDNA Project/data/whale-edna-relationships.csv")

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

get_model_asvs <- function(marker, whale, fit, taxonomic.level){
  if (!(tolower(marker) %in% c("16s", "18sv4", "18sv9"))){
    print("invalid marker")
  }
  if (!(tolower(whale) %in% c("blue whales", "fin whales", "humpback whales"))){
    print("invalid species")
  }
  
  if (tolower(taxonomic.level) == "order" | tolower(taxonomic.level) == "o"){
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
  
}

# only two at order level (across all markers/species)
get_model_asvs("16s", "humpback whales", "ss", "o")
get_model_asvs("18sv9", "humpback whales", "ss", "o")


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


# things get more interestig at class level

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
get_model_asvs("16s", "blue whales", "ss", "c")
get_model_asvs("16s", "humpback whales", "ss", "c")
get_model_asvs("18sv9", "humpback whales", "ss", "c")



