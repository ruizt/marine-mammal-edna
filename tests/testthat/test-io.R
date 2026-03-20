# test-io.R
# Verify that key data files load without error and have expected structure.
# Skipped automatically when files are absent.

skip_if_no_data <- function(path) {
  testthat::skip_if_not(file.exists(path), paste("Data file not present:", path))
}

test_that("example data loads and has expected structure", {
  load("data/example-data-16s.RData")
  expect_s3_class(example_16s, "data.frame")
  expect_s3_class(example_density, "data.frame")
  expect_equal(nrow(example_16s), 27L)
  expect_equal(ncol(example_16s), 501L)   # cruise + 500 ASVs
  expect_true("cruise" %in% names(example_16s))
  expect_true(all(c("bm", "bp", "mn") %in% names(example_density)))
  expect_length(example_asvs, 500L)
  expect_length(all_stable, 93L)
  expect_true(all(all_stable %in% names(example_16s)))
})

test_that("16S eDNA data loads and has expected structure", {
  skip_if_no_data("data/ncog16s.RData")
  load("data/ncog16s.RData")
  expect_s3_class(edna, "data.frame")
  expect_true("cruise" %in% names(edna))
  expect_gt(sum(startsWith(names(edna), "asv")), 1000L)
})

test_that("density estimates load and have expected structure", {
  skip_if_no_data("data/density-estimates.RData")
  load("data/density-estimates.RData")
  expect_s3_class(dens, "data.frame")
  expect_true(all(c("cruise", "bm", "bp", "mn") %in% names(dens)))
  expect_s3_class(dens_raw, "data.frame")
  expect_true(all(c("cruise", "species", "estimate", "season") %in% names(dens_raw)))
  expect_s3_class(dens_means, "data.frame")
  expect_equal(nrow(dens_means), 4L)   # four seasons
})

test_that("stability selection results in example data are well-formed", {
  load("data/example-data-16s.RData")
  expect_s3_class(sel_freq_bm, "data.frame")
  expect_true(all(c("species", "eta", "ncomp", "sel.asv", "n") %in%
                    names(sel_freq_bm)))
  expect_s3_class(candidate_sets_bm, "data.frame")
  expect_true("ss" %in% names(candidate_sets_bm))
  expect_gt(nrow(candidate_sets_bm), 0L)
})
