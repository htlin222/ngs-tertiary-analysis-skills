library(testthat)
library(here)
library(dplyr)
library(ggplot2)

source(here("R/plot_helpers.R"))

test_that("plot_vaf_distribution returns girafe for valid input", {
  mock <- tibble(gene = c("BRAF","TP53"), hgvsp = c("p.V600E","p.R175H"),
                 vaf = c(0.45, 0.30), consequence = c("missense_variant","missense_variant"),
                 pathogenicity = c("pathogenic","pathogenic"))
  expect_s3_class(plot_vaf_distribution(mock, 0.05), "girafe")
})

test_that("plot_vaf_distribution handles empty input", {
  empty <- tibble(gene=character(), hgvsp=character(), vaf=numeric(),
                  consequence=character(), pathogenicity=character())
  expect_s3_class(plot_vaf_distribution(empty, 0.05), "girafe")
})

test_that("plot_gene_coverage returns girafe", {
  mock <- tibble(gene=c("BRAF","TP53"), mean_coverage=c(450,150),
                 min_coverage=c(120,30), pct_above_200x=c(95,40))
  expect_s3_class(plot_gene_coverage(mock, 200), "girafe")
})

test_that("plot_cnv_genome_interactive returns plotly", {
  mock <- tibble(chromosome=c("7","17"), start=c(55e6L,75e5L), end=c(551e5L,77e5L),
                 gene=c("EGFR","TP53"), log2_ratio=c(2.1,-1.5),
                 copy_number=c(8,0.5), length=c(1e5L,2e5L),
                 type=c("AMPLIFICATION","DELETION"))
  expect_s3_class(plot_cnv_genome_interactive(mock), "plotly")
})

test_that("plot_cnv_genome_interactive handles empty input", {
  empty <- tibble(chromosome=character(), start=integer(), end=integer(),
                  gene=character(), log2_ratio=numeric(), copy_number=numeric(),
                  length=integer(), type=character())
  expect_s3_class(plot_cnv_genome_interactive(empty), "plotly")
})

test_that("plot_circos generates PNG file", {
  mock_cnv <- tibble(chromosome=c("chr7","chr17"), start=c(55e6L,75e5L),
                     end=c(552e5L,77e5L), gene=c("EGFR","TP53"),
                     log2_ratio=c(2.1,-1.5), type=c("AMPLIFICATION","DELETION"))
  mock_fusions <- tibble(gene_a="EML4", gene_b="ALK",
                         chr_a="chr2", pos_a=42500000L,
                         chr_b="chr2", pos_b=29400000L,
                         supporting_reads=15L)
  tmp <- tempfile(fileext=".png")
  plot_circos(mock_cnv, mock_fusions, output_file=tmp)
  expect_true(file.exists(tmp))
  expect_gt(file.size(tmp), 1000)
  unlink(tmp)
})

test_that("plot_fusion_arcs returns girafe", {
  mock <- tibble(gene_a="EML4", gene_b="ALK", supporting_reads=15L,
                 fusion_type="in-frame", known_fusion=TRUE)
  expect_s3_class(plot_fusion_arcs(mock), "girafe")
})

test_that("plot_fusion_arcs handles empty input", {
  empty <- tibble(gene_a=character(), gene_b=character(), supporting_reads=integer(),
                  fusion_type=character(), known_fusion=logical())
  expect_s3_class(plot_fusion_arcs(empty), "girafe")
})

test_that("plot_biomarker_gauges returns girafe", {
  tmb <- list(tmb_score=12.5, tmb_class="TMB-High", variant_count=24)
  msi <- list(msi_score=25, msi_status="MSI-H", unstable_sites=5, total_sites=20)
  hrd <- list(hrd_score=55, hrd_status="HRD-positive", loh_score=20, tai_score=18, lst_score=17)
  expect_s3_class(plot_biomarker_gauges(tmb, msi, hrd), "girafe")
})
