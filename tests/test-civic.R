library(testthat)
library(here)
library(dplyr)

# Source the CiVIC client
source(here("R/utils.R"))
source(here("R/civic_client.R"))

test_that("civic_search_variant returns tibble for known variant (BRAF V600E)", {
  skip_if_offline()

  result <- civic_search_variant("BRAF", "V600E")

  expect_s3_class(result, "tbl_df")
  expect_gt(nrow(result), 0)
  expect_true("variant_id" %in% names(result))
  expect_true("gene" %in% names(result))
  expect_true("variant_name" %in% names(result))
  expect_true("evidence_count" %in% names(result))
  expect_true("molecular_profile_id" %in% names(result))
})

test_that("civic_search_variant returns empty tibble for unknown variant", {
  skip_if_offline()

  result <- civic_search_variant("FAKEGENE123", "X999Z")

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
  expect_true("variant_id" %in% names(result))
  expect_true("gene" %in% names(result))
  expect_true("variant_name" %in% names(result))
  expect_true("evidence_count" %in% names(result))
  expect_true("molecular_profile_id" %in% names(result))
})

test_that("civic_get_evidence returns evidence items", {
  skip_if_offline()

  # First search for BRAF V600E to get IDs
  search <- civic_search_variant("BRAF", "V600E")
  skip_if(nrow(search) == 0, "BRAF V600E not found in CiVIC")

  result <- civic_get_evidence(search$molecular_profile_id[1])

  expect_s3_class(result, "tbl_df")
  expect_true("evidence_id" %in% names(result))
  expect_true("evidence_type" %in% names(result))
  expect_true("evidence_level" %in% names(result))
  expect_true("evidence_direction" %in% names(result))
  expect_true("significance" %in% names(result))
  expect_true("disease" %in% names(result))
  expect_true("therapies" %in% names(result))
  expect_true("source_citation" %in% names(result))
  expect_true("source_url" %in% names(result))
  expect_true("description" %in% names(result))
})

test_that("civic_get_gene_summary returns gene info", {
  skip_if_offline()

  result <- civic_get_gene_summary("BRAF")

  expect_type(result, "list")
  expect_false(is.na(result$gene_id))
  expect_equal(result$gene, "BRAF")
  expect_true("description" %in% names(result))
  expect_true("official_name" %in% names(result))
})

test_that("civic_get_assertions returns AMP assertions", {
  skip_if_offline()

  result <- civic_get_assertions("BRAF")

  expect_s3_class(result, "tbl_df")
  expect_true("assertion_id" %in% names(result))
  expect_true("amp_level" %in% names(result))
  expect_true("assertion_type" %in% names(result))
  expect_true("assertion_direction" %in% names(result))
  expect_true("significance" %in% names(result))
  expect_true("disease" %in% names(result))
  expect_true("therapies" %in% names(result))
  expect_true("variant_name" %in% names(result))
  expect_true("gene" %in% names(result))
  expect_true("description" %in% names(result))
  expect_true("nccn_guideline" %in% names(result))
  expect_true("regulatory_approval" %in% names(result))
  expect_true("fda_companion_test" %in% names(result))
})

test_that("civic_annotate_variants handles empty input", {
  result <- civic_annotate_variants(
    tibble(gene = character(), variant = character()),
    sample_id = "TEST_EMPTY"
  )

  expect_type(result, "list")
  expect_s3_class(result$variant_evidence, "tbl_df")
  expect_equal(nrow(result$variant_evidence), 0)
  expect_type(result$gene_summaries, "list")
  expect_equal(length(result$gene_summaries), 0)
  expect_s3_class(result$assertions, "tbl_df")
  expect_equal(nrow(result$assertions), 0)
})

# ── Batch GraphQL tests ──────────────────────────────────────────────────────

test_that("civic_batch_resolve_genes resolves known genes", {
  skip_if_offline()

  result <- civic_batch_resolve_genes(c("BRAF", "TP53", "FAKEGENE999"))

  expect_type(result, "list")
  expect_false(is.na(result[["BRAF"]]))
  expect_false(is.na(result[["TP53"]]))
  expect_true(is.na(result[["FAKEGENE999"]]))
})

test_that("civic_batch_resolve_genes handles empty input", {
  result <- civic_batch_resolve_genes(character(0))
  expect_type(result, "list")
  expect_equal(length(result), 0)
})

test_that("civic_batch_search_variants finds known variants", {
  skip_if_offline()

  gene_ids <- civic_batch_resolve_genes(c("BRAF"))
  skip_if(is.na(gene_ids[["BRAF"]]), "BRAF not resolved")

  queries <- tibble(
    gene = "BRAF",
    variant = "V600E",
    gene_id = as.integer(gene_ids[["BRAF"]])
  )
  result <- civic_batch_search_variants(queries)

  expect_s3_class(result, "tbl_df")
  expect_gt(nrow(result), 0)
  expect_true("molecular_profile_id" %in% names(result))
  expect_true("query_gene" %in% names(result))
  expect_true("query_variant" %in% names(result))
})

test_that("civic_batch_search_variants handles empty input", {
  result <- civic_batch_search_variants(
    tibble(gene = character(), variant = character(), gene_id = integer())
  )
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
})

test_that("civic_sequential_annotate_evidence works as fallback", {
  result <- civic_sequential_annotate_evidence(
    tibble(gene = character(), variant = character())
  )
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
})
