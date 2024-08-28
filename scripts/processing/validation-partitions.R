library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(fs)
library(readr)
data_dir <- 'data/processed/'
out_dir <- paste(data_dir, '_combined-partitions/', sep = '')
dir_create(out_dir)

## 18Sv9-SS PARTITIONS ---------------------------------------------------------

# read in edna and sighting partitions
load('data/processed/_cv/ncog18sv9-partitions.RData')
load('data/processed/_cv/mm-sightings-partitions.RData')

# combine regular (non-nested) validation partitions
loo_partitions <- left_join(loo_sightings, 
                            loo_edna, 
                            join_by(test.id), 
                            suffix = c('.mm', '.edna')) |>
  mutate(train = map2(train.mm, train.edna, ~inner_join(.x, .y, join_by('cruise'))),
         test = map2(test.mm, test.edna, ~inner_join(.x, .y, join_by('cruise')))) |>
  select(test.id, test.season, train, test)

# store
write_rds(loo_partitions, 
          file = paste(out_dir, 'partitions-18sv9-ss.rds', sep = ''))

# read in nested partitions
load('data/processed/_cv/ncog18sv9-nested-partitions.RData')
load('data/processed/_cv/mm-sightings-nested-partitions.RData')

# combine nested validation partitions
loo_partitions_nested <- inner_join(loo_sightings_nested,
                                    loo_edna_nested,
                                    join_by(outer.id, inner.id),
                                    suffix = c('.mm', '.edna')) |>
  mutate(train = map2(train.mm, train.edna,
                      ~join(.x, .y, on = 'cruise', how = 'inner'))) |>
  select(outer.id, inner.id, train)

# store partitions
write_rds(loo_partitions_nested,
          file = paste(out_dir, '_nested-partitions-18sv9-ss.rds', sep = ''))


## 18Sv4-SS PARTITIONS ---------------------------------------------------------


## 16S-SS PARTITIONS ---------------------------------------------------------
