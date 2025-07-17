library(tidyverse)
out_dir <- '_tbl/'

load('data/processed/ncog16s.RData')
data_16s <- edna
metadata_16s <- sample_metadata
ann_16s <- asv_taxa |> 
  rename(domain = d,
         phylum = p,
         class = c,
         order = o,
         family = f,
         genus = g,
         species = s)

load('data/processed/ncog18sv4.RData')
data_18sv4 <- edna
metadata_18sv4 <- sample_metadata
ann_18sv4 <- asv_taxa |> 
  rename(domain = d,
         phylum = p,
         class = c,
         order = o,
         family = f,
         genus = g,
         species = s)

load('data/processed/ncog18sv9.RData')
data_18sv9 <- edna
metadata_18sv9 <- sample_metadata
ann_18sv9 <- asv_taxa |> 
  rename(domain = d,
         phylum = p,
         class = c,
         order = o,
         family = f,
         genus = g,
         species = s)

load('data/processed/density-estimates.RData')
load('data/processed/sightings.RData')

sheets <- list("annotations-16s" = ann_16s,
               "data-16s" = data_16s,
               "metadata-16s" = metadata_16s,
               "annotations-18sv4" = ann_18sv4,
               "data-18sv4" = data_18sv4,
               "metadata-18sv4" = metadata_18sv4,
               "annotations-18sv9" = ann_18sv9,
               "data-18sv9" = data_18sv9,
               "metadata-18sv9" = metadata_18sv9,
               "sightings" = sightings,
               "density-raw" = dens_raw,
               "density-processed" = dens)

# export as excel workbook
writexl::write_xlsx(sheets, paste(out_dir, 'data.xlsx', sep = ''))