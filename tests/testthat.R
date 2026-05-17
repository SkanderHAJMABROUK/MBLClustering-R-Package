library(testthat)
library(MBLClustering)

test_that("compare_models works", {
  data <- matrix(rnorm(200), ncol = 2)
  comparison <- compare_models(data, K_range = 1:3, n_init = 2)
  
  expect_s3_class(comparison, "gmm_comparison")
  expect_true("K" %in% names(comparison))
  expect_true("criterion_value" %in% names(comparison))
})

test_that("plot_diagnostics works", {
  data <- matrix(rnorm(200), ncol = 2)
  result <- gmm_outliers(data, K = 2, n_init = 2)
  
  # Devrait fonctionner sans erreur
  expect_error(plot_diagnostics(result, data = data, type = "basic"), NA)
})

test_that("calculate_metrics works", {
  data <- matrix(rnorm(200), ncol = 2)
  result <- gmm_outliers(data, K = 2, n_init = 2)
  metrics <- calculate_metrics(result)
  
  expect_s3_class(metrics, "gmm_metrics")
  expect_true("basic" %in% names(metrics))
  expect_true("clustering" %in% names(metrics))
})
