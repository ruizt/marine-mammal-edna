# test-paths.R
# Verify that expected data files exist in data/ after running setup-data.R.
# These tests catch path-reference regressions without running any analysis.
# Run after setup-data.R; will be skipped automatically if files are absent
# (so they don't break CI for users who haven't downloaded the data).

p <- function(...) here::here(...)  # anchor all paths to project root

skip_if_no_data <- function(path) {
  testthat::skip_if_not(file.exists(p(path)),
                        paste("Data file not present:", path))
}

test_that("core processed data files exist", {
  skip_if_no_data("data/ncog16s.RData")
  skip_if_no_data("data/ncog18sv4.RData")
  skip_if_no_data("data/ncog18sv9.RData")
  skip_if_no_data("data/density-estimates.RData")
  skip_if_no_data("data/sightings.RData")
  expect_true(file.exists(p("data/ncog16s.RData")))
  expect_true(file.exists(p("data/density-estimates.RData")))
})

test_that("line transect data files exist", {
  skip_if_no_data("data/sample_table.Rds")
  expect_true(file.exists(p("data/sample_table.Rds")))
  expect_true(file.exists(p("data/region_table.Rds")))
  expect_true(file.exists(p("data/obs_table.Rds")))
})

test_that("always-present repository files exist", {
  expect_true(file.exists(p("data/example-data-16s.RData")))
  expect_true(file.exists(p("data/gshhg.RData")))
})

test_that("partition directories exist (if downloaded)", {
  skip_if_no_data("data/_partitions")
  expect_true(dir.exists(p("data/_partitions")))
  expect_true(dir.exists(p("data/_combined-partitions")))
})
