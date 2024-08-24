library(tibble)
library(dplyr)
library(tidyr)
library(purrr)
library(fs)
library(readr)

## 18Sv9-SS PARTITIONS ---------------------------------------------------------

part_dir <- 'data/processed/_cv'
out_dir <- 'rslt/loocv/18sv9-ss/validation'
dir_create(out_dir)

# read in and combine validation partitions
paste(part_dir, 'ncog18sv9-partitions.RData', sep = '/') |> load()
paste(part_dir, 'mm-sightings-partitions.RData', sep = '/') |> load()
loo_partitions <- inner_join(loo_sightings_nested,
                             loo_edna_nested,
                             join_by(test.outer.cruise, test.inner.cruise),
                             suffix = c('.mm', '.edna')) |>
  mutate(train.inner = map2(train.inner.mm, train.inner.edna,
                            ~inner_join(.x, .y, join_by('cruise'))),
         test.inner = map2(test.inner.mm, test.inner.edna,
                           ~inner_join(.x, .y, join_by('cruise')))) |>
  select(test.outer.cruise, test.inner.cruise, test.inner.season,
         train.inner, test.inner)

# store partitions
write_rds(loo_partitions,
          file = paste(out_dir, 'partitions.rds', sep = '/'))

## 18Sv4-SS PARTITIONS ---------------------------------------------------------

part_dir <- 'data/processed/_cv'
out_dir <- 'rslt/loocv/18sv4-ss/validation'
dir_create(out_dir)

# read in and combine validation partitions
paste(part_dir, 'ncog18sv4-partitions.RData', sep = '/') |> load()
paste(part_dir, 'mm-sightings-partitions.RData', sep = '/') |> load()
loo_partitions <- inner_join(loo_sightings_nested,
                             loo_edna_nested,
                             join_by(test.outer.cruise, test.inner.cruise),
                             suffix = c('.mm', '.edna')) |>
  mutate(train.inner = map2(train.inner.mm, train.inner.edna,
                            ~inner_join(.x, .y, join_by('cruise'))),
         test.inner = map2(test.inner.mm, test.inner.edna,
                           ~inner_join(.x, .y, join_by('cruise')))) |>
  select(test.outer.cruise, test.inner.cruise, test.inner.season,
         train.inner, test.inner)

# store partitions
write_rds(loo_partitions,
          file = paste(out_dir, 'partitions.rds', sep = '/'))

## 16S-SS PARTITIONS -----------------------------------------------------------

part_dir <- 'data/processed/_cv'
out_dir <- 'rslt/loocv/16s-ss/validation'
dir_create(out_dir)

# read in and combine validation partitions
paste(part_dir, 'ncog16s-partitions.RData', sep = '/') |> load()
paste(part_dir, 'mm-sightings-partitions.RData', sep = '/') |> load()
loo_partitions <- inner_join(loo_sightings_nested,
                             loo_edna_nested,
                             join_by(test.outer.cruise, test.inner.cruise),
                             suffix = c('.mm', '.edna')) |>
  mutate(train.inner = map2(train.inner.mm, train.inner.edna,
                            ~inner_join(.x, .y, join_by('cruise'))),
         test.inner = map2(test.inner.mm, test.inner.edna,
                           ~inner_join(.x, .y, join_by('cruise')))) |>
  select(test.outer.cruise, test.inner.cruise, test.inner.season,
         train.inner, test.inner)

# store partitions
write_rds(loo_partitions,
          file = paste(out_dir, 'partitions.rds', sep = '/'))