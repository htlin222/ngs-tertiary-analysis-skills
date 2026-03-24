library(testthat)
library(here)
library(yaml)
library(fs)
library(withr)

# Source the utilities
source(here("R/utils.R"))

test_that("load_config reads YAML correctly", {
  config_path <- here("config/default.yaml")
  config <- load_config(config_path)

  # Check that key fields exist
  expect_true(is.list(config))
  expect_true("sample" %in% names(config))
  expect_true("qc" %in% names(config))
  expect_true("variant_calling" %in% names(config))
  expect_true("annotation" %in% names(config))

  # Check nested access
  expect_equal(config$sample$id, "SAMPLE_001")
  expect_equal(config$sample$tumor_type, "NSCLC")
  expect_equal(config$qc$min_mean_coverage, 200)
})

test_that("set_nested works for deep keys", {
  test_list <- list(a = list(b = list(c = 1)))

  # Test setting at leaf
  result <- set_nested(test_list, c("a", "b", "c"), 42)
  expect_equal(result$a$b$c, 42)

  # Test setting new path
  result2 <- set_nested(list(), c("x", "y", "z"), "value")
  expect_equal(result2$x$y$z, "value")

  # Test single level
  result3 <- set_nested(list(), c("key"), "val")
  expect_equal(result3$key, "val")
})

test_that("stage_output_dir creates directory", {
  temp_base <- local_tempdir()

  result_dir <- stage_output_dir("TEST_SAMPLE", "99-test", base_dir = temp_base)

  # Check that directory was created
  expect_true(dir.exists(result_dir))

  # Check path structure
  expect_true(grepl("TEST_SAMPLE", result_dir))
  expect_true(grepl("99-test", result_dir))

  # Second call should not error
  result_dir2 <- stage_output_dir("TEST_SAMPLE", "99-test", base_dir = temp_base)
  expect_equal(result_dir, result_dir2)
})

test_that("format_hgvs_protein formats correctly", {
  # Standard case
  result <- format_hgvs_protein("BRAF", "V600E")
  expect_equal(result, "BRAF p.V600E")

  # Other genes
  result2 <- format_hgvs_protein("TP53", "R175H")
  expect_equal(result2, "TP53 p.R175H")

  # Deletion notation
  result3 <- format_hgvs_protein("EGFR", "L858R")
  expect_equal(result3, "EGFR p.L858R")
})

test_that("adjust_vaf_for_purity handles edge cases", {
  # Normal case
  result <- adjust_vaf_for_purity(0.3, 0.8)
  expect_true(is.numeric(result))
  expect_true(result > 0)
  expect_true(result <= 1.0)

  # VAF of 0
  result_zero <- adjust_vaf_for_purity(0, 0.8)
  expect_equal(result_zero, 0)

  # Invalid purity (zero) - should return raw VAF
  result_invalid_zero <- adjust_vaf_for_purity(0.3, 0)
  expect_equal(result_invalid_zero, 0.3)

  # Invalid purity (>1) - should return raw VAF
  result_invalid_high <- adjust_vaf_for_purity(0.3, 1.5)
  expect_equal(result_invalid_high, 0.3)

  # High purity, high VAF
  result_high <- adjust_vaf_for_purity(0.5, 1.0)
  expect_true(result_high <= 1.0)
})
