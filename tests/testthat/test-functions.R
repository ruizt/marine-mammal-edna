# test-functions.R
# Lightweight unit tests for R/ package functions using a small synthetic
# dataset. These run fast and do not require any downloaded data.

library(spls)
library(pls)

make_toy_data <- function(n = 10, p = 30, seed = 5571) {
  set.seed(seed)
  asvs <- matrix(runif(n * p), nrow = n,
                 dimnames = list(NULL, paste0("asv.", seq_len(p))))
  tibble::tibble(
    cruise = paste0("cruise", seq_len(n)),
    bm     = rnorm(n),
    bp     = rnorm(n),
    mn     = rnorm(n),
    as.data.frame(asvs)
  )
}

test_that("eta_grid returns a numeric vector of the right length", {
  g <- eta_grid(n = 10)
  expect_length(g, 10L)
  expect_true(is.numeric(g))
  expect_true(all(diff(g) > 0))   # increasing
})

test_that("loo_partitions returns one fold per cruise", {
  toy <- make_toy_data()
  parts <- loo_partitions(toy)
  expect_equal(nrow(parts), nrow(toy))
  expect_equal(parts$test.id, toy$cruise)
  expect_equal(nrow(parts$train[[1]]), nrow(toy) - 1L)
  expect_equal(nrow(parts$test[[1]]),  1L)
})

test_that("fit_spls_partition returns expected columns and types", {
  toy   <- make_toy_data()
  parts <- loo_partitions(toy)
  result <- fit_spls_partition("bm", eta = 0.5, ncomp = 2, data = toy,
                               partitions = parts)
  expect_s3_class(result, "data.frame")
  expect_true(all(c("species", "eta", "ncomp", "obs.id", "sel.asv") %in%
                    names(result)))
  expect_equal(nrow(result), nrow(toy))
  expect_type(result$sel.asv, "list")
})

test_that("subset_to_stable and fit_pls run without error", {
  toy <- make_toy_data()
  ss  <- paste0("asv.", 1:5)
  d   <- subset_to_stable(toy, ss, "bm")
  expect_equal(names(d), c("y", ss))
  fit <- fit_pls(d, ncomp = 2)
  expect_s3_class(fit, "mvr")
})
