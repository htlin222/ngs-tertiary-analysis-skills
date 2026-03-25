library(testthat)
library(here)
library(dplyr)

# Source the ESMO helpers
source(here("R/utils.R"))
source(here("R/esmo_helpers.R"))

test_that("oncokb_to_escat maps all levels correctly", {
  # FDA-approved levels
  expect_equal(oncokb_to_escat("LEVEL_1"), "I")
  expect_equal(oncokb_to_escat("LEVEL_2"), "I")

  # Investigational levels
  expect_equal(oncokb_to_escat("LEVEL_3A"), "II")
  expect_equal(oncokb_to_escat("LEVEL_3B"), "III")
  expect_equal(oncokb_to_escat("LEVEL_4"), "IV")

  # Resistance levels
  expect_equal(oncokb_to_escat("LEVEL_R1"), "I-R")
  expect_equal(oncokb_to_escat("LEVEL_R2"), "II-R")
})

test_that("oncokb_to_escat handles NA and NULL", {
  expect_equal(oncokb_to_escat(NA), "X")
  expect_equal(oncokb_to_escat(NULL), "X")
  expect_equal(oncokb_to_escat("UNKNOWN_LEVEL"), "X")
})

test_that("escat_description returns descriptions", {
  # Tier I
  desc_i <- escat_description("I")
  expect_true(grepl("routine clinical", desc_i, ignore.case = TRUE))

  # Tier X
  desc_x <- escat_description("X")
  expect_true(grepl("no evidence", desc_x, ignore.case = TRUE))

  # Tier II
  desc_ii <- escat_description("II")
  expect_true(grepl("investigational", desc_ii, ignore.case = TRUE))

  # Unknown tier
  desc_unknown <- escat_description("UNKNOWN")
  expect_equal(desc_unknown, "Unclassified")
})

test_that("escat_badge_html generates valid HTML", {
  html <- escat_badge_html("I")

  # Check HTML structure
  expect_true(grepl("span", html))
  expect_true(grepl("background-color", html))
  expect_true(grepl("ESCAT I", html))

  # Check color code
  expect_true(grepl("#d32f2f", html))  # Red for Tier I

  # Check that each tier has a color
  html_x <- escat_badge_html("X")
  expect_true(grepl("#9e9e9e", html_x))  # Grey for Tier X
})

test_that("format_variant_esmo formats correctly", {
  # Gene + NA protein (minimum args)
  result1 <- format_variant_esmo("BRAF", hgvsp = NA)
  expect_equal(result1, "BRAF")

  # Gene + protein change
  result2 <- format_variant_esmo("BRAF", "p.Val600Glu")
  expect_equal(result2, "BRAF p.Val600Glu")

  # Gene + protein + coding
  result3 <- format_variant_esmo("BRAF", "p.Val600Glu", "c.1799T>A")
  expect_equal(result3, "BRAF p.Val600Glu (c.1799T>A)")

  # Gene + protein + coding + VAF
  result4 <- format_variant_esmo("BRAF", "p.Val600Glu", "c.1799T>A", 0.35)
  expect_true(grepl("BRAF", result4))
  expect_true(grepl("p.Val600Glu", result4))
  expect_true(grepl("c.1799T>A", result4))
  expect_true(grepl("35", result4))  # VAF formatted

  # VAF only (hgvsp required arg)
  result5 <- format_variant_esmo("TP53", hgvsp = NA, vaf = 0.55)
  expect_true(grepl("TP53", result5))
  expect_true(grepl("55", result5))
})

test_that("format_fusion_esmo uses double colon", {
  # Basic fusion
  result1 <- format_fusion_esmo("ALK", "EML4")
  expect_equal(result1, "ALK::EML4")

  # Fusion with exon information
  result2 <- format_fusion_esmo("ALK", "EML4", exon_a = 20, exon_b = 1)
  expect_equal(result2, "ALK::EML4 (exon 20::exon 1)")

  # Partial exon info (should not add)
  result3 <- format_fusion_esmo("BRAF", "KIAA1549", exon_a = 8)
  expect_equal(result3, "BRAF::KIAA1549")
})

test_that("coverage_gaps identifies low-coverage genes", {
  # Create mock coverage data
  coverage_df <- dplyr::tibble(
    gene = c("BRCA1", "BRCA2", "TP53", "EGFR", "MYC"),
    mean_coverage = c(150, 250, 180, 350, 80)
  )

  # Find gaps with threshold 200
  gaps <- coverage_gaps(coverage_df, min_coverage = 200)

  # Should return 3 genes below 200
  expect_equal(nrow(gaps), 3)
  expect_true("BRCA1" %in% gaps$gene)
  expect_true("TP53" %in% gaps$gene)
  expect_true("MYC" %in% gaps$gene)

  # Verify sorted by coverage
  expect_true(gaps$mean_coverage[1] < gaps$mean_coverage[2])

  # No gaps scenario
  good_coverage <- dplyr::tibble(
    gene = c("GENE1", "GENE2", "GENE3"),
    mean_coverage = c(300, 400, 500)
  )
  gaps_none <- coverage_gaps(good_coverage, min_coverage = 200)
  expect_equal(nrow(gaps_none), 0)
})
