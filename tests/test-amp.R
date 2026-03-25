library(testthat)
library(here)
library(dplyr)

source(here("R/amp_classification.R"))

test_that("classify_amp_tier Tier I Level A for FDA-approved", {
  result <- classify_amp_tier(
    oncokb_oncogenic = "Oncogenic", oncokb_level = "LEVEL_1",
    civic_amp_level = "TIER_I_LEVEL_A", vep_impact = "HIGH", clinvar = "pathogenic"
  )
  expect_equal(result$amp_tier, "Tier I")
  expect_equal(result$amp_level, "Level A")
  expect_true(grepl("LEVEL_1", result$amp_evidence))
})

test_that("classify_amp_tier Tier I Level A for LEVEL_2", {
  result <- classify_amp_tier(
    oncokb_oncogenic = "Oncogenic", oncokb_level = "LEVEL_2",
    civic_amp_level = NA, vep_impact = "HIGH", clinvar = "pathogenic"
  )
  expect_equal(result$amp_tier, "Tier I")
  expect_equal(result$amp_level, "Level A")
})

test_that("classify_amp_tier Tier I Level B for well-powered studies", {
  result <- classify_amp_tier(
    oncokb_oncogenic = "Oncogenic", oncokb_level = "LEVEL_3A",
    civic_amp_level = NA, vep_impact = "HIGH", clinvar = "pathogenic"
  )
  expect_equal(result$amp_tier, "Tier I")
  expect_equal(result$amp_level, "Level B")
  expect_true(grepl("LEVEL_3A", result$amp_evidence))
})

test_that("classify_amp_tier Tier II Level C for different tumor type", {
  result <- classify_amp_tier(
    oncokb_oncogenic = "Oncogenic", oncokb_level = "LEVEL_3B",
    civic_amp_level = NA, vep_impact = "HIGH", clinvar = NA
  )
  expect_equal(result$amp_tier, "Tier II")
  expect_equal(result$amp_level, "Level C")
})

test_that("classify_amp_tier Tier II Level D for LEVEL_4", {
  result <- classify_amp_tier(
    oncokb_oncogenic = "Oncogenic", oncokb_level = "LEVEL_4",
    civic_amp_level = NA, vep_impact = "MODERATE", clinvar = NA
  )
  expect_equal(result$amp_tier, "Tier II")
  expect_equal(result$amp_level, "Level D")
})

test_that("classify_amp_tier Tier II Level D for oncogenic + HIGH impact", {
  result <- classify_amp_tier(
    oncokb_oncogenic = "Oncogenic", oncokb_level = NA,
    civic_amp_level = NA, vep_impact = "HIGH", clinvar = NA
  )
  expect_equal(result$amp_tier, "Tier II")
  expect_equal(result$amp_level, "Level D")
  expect_true(grepl("HIGH VEP impact", result$amp_evidence))
})

test_that("classify_amp_tier Tier III for VUS", {
  result <- classify_amp_tier(
    oncokb_oncogenic = "Likely Neutral", oncokb_level = NA,
    civic_amp_level = NA, vep_impact = "MODERATE", clinvar = "uncertain_significance"
  )
  expect_equal(result$amp_tier, "Tier III")
  expect_equal(result$amp_level, "VUS")
})

test_that("classify_amp_tier Tier IV for benign", {
  result <- classify_amp_tier(
    oncokb_oncogenic = "Neutral", oncokb_level = NA,
    civic_amp_level = NA, vep_impact = "LOW", clinvar = "benign"
  )
  expect_equal(result$amp_tier, "Tier IV")
  expect_equal(result$amp_level, "Benign")
})

test_that("classify_amp_tier Tier IV for Neutral oncogenicity", {
  result <- classify_amp_tier(
    oncokb_oncogenic = "Neutral", oncokb_level = NA,
    civic_amp_level = NA, vep_impact = "LOW", clinvar = NA
  )
  expect_equal(result$amp_tier, "Tier IV")
  expect_equal(result$amp_level, "Benign")
})

test_that("classify_amp_tier handles all NA inputs gracefully", {
  result <- classify_amp_tier()
  expect_equal(result$amp_tier, "Tier III")
  expect_equal(result$amp_level, "VUS")
})

test_that("classify_all_amp returns tibble with amp columns", {
  mock_variants <- tibble(
    gene = c("BRAF", "TP53"),
    hgvsp = c("p.V600E", "p.R175H"),
    oncogenic = c("Oncogenic", "Oncogenic"),
    sensitive_level = c("LEVEL_1", "LEVEL_3A"),
    impact = c("HIGH", "HIGH"),
    clinvar_significance = c("pathogenic", "pathogenic")
  )
  mock_civic <- tibble(
    gene = "BRAF", variant_name = "V600E", amp_level = "TIER_I_LEVEL_A"
  )
  result <- classify_all_amp(mock_variants, mock_civic)
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("amp_tier", "amp_level", "amp_evidence") %in% names(result)))
  expect_equal(nrow(result), 2)
  # BRAF V600E should be Tier I Level A
  expect_equal(result$amp_tier[1], "Tier I")
  expect_equal(result$amp_level[1], "Level A")
  # TP53 should be Tier I Level B (LEVEL_3A)
  expect_equal(result$amp_tier[2], "Tier I")
  expect_equal(result$amp_level[2], "Level B")
})

test_that("classify_all_amp handles empty input", {
  empty <- tibble(
    gene = character(), hgvsp = character(), oncogenic = character(),
    sensitive_level = character(), impact = character(),
    clinvar_significance = character()
  )
  result <- classify_all_amp(empty)
  expect_equal(nrow(result), 0)
  expect_true(all(c("amp_tier", "amp_level", "amp_evidence") %in% names(result)))
})

test_that("classify_all_amp works without civic_assertions", {
  mock_variants <- tibble(
    gene = "EGFR", hgvsp = "p.L858R", oncogenic = "Oncogenic",
    sensitive_level = "LEVEL_1", impact = "HIGH",
    clinvar_significance = "pathogenic"
  )
  result <- classify_all_amp(mock_variants)
  expect_equal(nrow(result), 1)
  expect_equal(result$amp_tier[1], "Tier I")
  expect_equal(result$amp_level[1], "Level A")
})

test_that("classify_all_amp CiVIC matching strips p. prefix", {
  mock_variants <- tibble(
    gene = "BRAF", hgvsp = "p.V600E", oncogenic = "Likely Neutral",
    sensitive_level = NA_character_, impact = "MODERATE",
    clinvar_significance = NA_character_
  )
  # CiVIC has a Tier I Level A assertion that should override VUS

  mock_civic <- tibble(
    gene = "BRAF", variant_name = "V600E", amp_level = "TIER_I_LEVEL_A"
  )
  result <- classify_all_amp(mock_variants, mock_civic)
  expect_equal(result$amp_tier[1], "Tier I")
  expect_equal(result$amp_level[1], "Level A")
})
