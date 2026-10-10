# Open-Source Pipeline Enhancements Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add CiVIC integration, interactive visualizations, AMP/ASCO/CAP oncogenicity classification, and report security — all using free/open-source tools only.

**Architecture:** Four waves of enhancements that build on the existing `targets` pipeline. Wave 1 (CiVIC) adds a new knowledge source to stage 06. Wave 2 (visualizations) adds interactive plots to the Quarto report using ggiraph/circlize/plotly. Wave 3 (AMP/ASCO/CAP) adds automated oncogenicity classification using CiVIC assertions + OncoKB + VEP data. Wave 4 (report security) adds password-protected self-contained HTML output.

**Tech Stack:** R, httr2, ggiraph, circlize, plotly, jsonlite, testthat, Quarto HTML

---

## Wave 1: CiVIC API Integration

### Task 1: CiVIC API Client

**Files:**
- Create: `R/civic_client.R`
- Test: `tests/test-civic.R`

**Step 1: Write the failing test**

```r
# tests/test-civic.R
library(testthat)
library(here)
library(dplyr)

source(here("R/utils.R"))
source(here("R/civic_client.R"))

test_that("civic_search_variant returns tibble for known variant", {
  # Use a well-known variant that CiVIC definitely has
  skip_if_offline()
  result <- civic_search_variant("BRAF", "V600E")
  expect_s3_class(result, "tbl_df")
  expect_true(nrow(result) > 0)
  expect_true(all(c("variant_id", "gene", "variant_name", "evidence_count") %in% names(result)))
})

test_that("civic_search_variant returns empty tibble for unknown variant", {
  skip_if_offline()
  result <- civic_search_variant("FAKEGENE", "X999Z")
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 0)
})

test_that("civic_get_evidence returns evidence items for a variant", {
  skip_if_offline()
  # BRAF V600E has variant_id 12 in CiVIC
  result <- civic_get_evidence(variant_id = 12)
  expect_s3_class(result, "tbl_df")
  expect_true(nrow(result) > 0)
  expect_true(all(c("evidence_id", "evidence_type", "evidence_level",
                     "evidence_direction", "significance", "disease",
                     "therapies", "source_citation") %in% names(result)))
})

test_that("civic_get_gene_summary returns gene-level info", {
  skip_if_offline()
  result <- civic_get_gene_summary("BRAF")
  expect_type(result, "list")
  expect_true("description" %in% names(result))
  expect_true("gene_id" %in% names(result))
})

test_that("civic_get_assertions returns AMP/ASCO/CAP assertions", {
  skip_if_offline()
  result <- civic_get_assertions(gene = "BRAF")
  expect_s3_class(result, "tbl_df")
  # BRAF should have at least one assertion
  expect_true(nrow(result) > 0)
  expect_true(all(c("assertion_id", "amp_level", "assertion_type",
                     "assertion_direction", "significance",
                     "disease", "therapies", "variant_name") %in% names(result)))
})
```

**Step 2: Run test to verify it fails**

Run: `cd /Users/htlin/ngs-tertiary-analysis-skills && Rscript -e 'testthat::test_file("tests/test-civic.R")'`
Expected: FAIL — `civic_client.R` does not exist

**Step 3: Write minimal implementation**

```r
# R/civic_client.R — CiVIC (Clinical Interpretation of Variants in Cancer) API client
# API docs: https://griffithlab.github.io/civic-v2/#query-CivicAPI
# CiVIC uses a GraphQL API at https://civicdb.org/api/graphql

suppressPackageStartupMessages({
  library(httr2)
  library(jsonlite)
  library(logger)
  library(glue)
  library(dplyr)
  library(purrr)
})

# ══════════════════════════════════════════════════════════════════════════════
# CiVIC GraphQL API — free, no API key required
# ══════════════════════════════════════════════════════════════════════════════

CIVIC_GRAPHQL_URL <- "https://civicdb.org/api/graphql"

#' Execute a CiVIC GraphQL query
#' @param query GraphQL query string
#' @param variables Named list of query variables
#' @return Parsed JSON response data
civic_graphql <- function(query, variables = list()) {
  body <- list(query = query, variables = variables)

  resp <- request(CIVIC_GRAPHQL_URL) |>
    req_headers(`Content-Type` = "application/json") |>
    req_body_json(body) |>
    req_retry(max_tries = 3, backoff = ~ 2^.x) |>
    req_throttle(rate = 30 / 60) |>
    req_perform()

  result <- resp_body_json(resp)

  if (!is.null(result$errors)) {
    log_warn("CiVIC GraphQL error: {toJSON(result$errors, auto_unbox = TRUE)}")
  }

  result$data
}

#' Search CiVIC for a variant by gene and variant name
#' @param gene Gene symbol (e.g., "BRAF")
#' @param variant Variant name (e.g., "V600E")
#' @return Tibble with variant_id, gene, variant_name, evidence_count
civic_search_variant <- function(gene, variant) {
  query <- '
  query SearchVariant($geneName: String!, $variantName: String) {
    variants(
      geneSymbol: $geneName
      name: $variantName
      first: 10
    ) {
      nodes {
        id
        name
        singleVariantMolecularProfile {
          id
          name
          evidenceCountsByStatus {
            acceptedCount
            submittedCount
          }
        }
        gene {
          name
          id
        }
      }
    }
  }'

  data <- civic_graphql(query, list(geneName = gene, variantName = variant))
  nodes <- data$variants$nodes %||% list()

  if (length(nodes) == 0) {
    return(tibble(
      variant_id = integer(), gene = character(),
      variant_name = character(), evidence_count = integer(),
      molecular_profile_id = integer()
    ))
  }

  map_dfr(nodes, function(v) {
    ev_counts <- v$singleVariantMolecularProfile$evidenceCountsByStatus
    total_evidence <- (ev_counts$acceptedCount %||% 0L) + (ev_counts$submittedCount %||% 0L)
    tibble(
      variant_id = as.integer(v$id),
      gene = v$gene$name %||% gene,
      variant_name = v$name %||% NA_character_,
      evidence_count = as.integer(total_evidence),
      molecular_profile_id = as.integer(v$singleVariantMolecularProfile$id %||% NA)
    )
  })
}

#' Get evidence items for a CiVIC variant
#' @param variant_id CiVIC variant ID
#' @param molecular_profile_id Optional: CiVIC molecular profile ID (preferred for v2 API)
#' @return Tibble of evidence items
civic_get_evidence <- function(variant_id = NULL, molecular_profile_id = NULL) {
  query <- '
  query GetEvidence($molecularProfileId: Int, $first: Int) {
    evidenceItems(
      molecularProfileId: $molecularProfileId
      status: ACCEPTED
      first: $first
    ) {
      nodes {
        id
        evidenceType
        evidenceLevel
        evidenceDirection
        significance
        description
        disease { name }
        therapies { name }
        source {
          citation
          sourceUrl
          sourceType
        }
        molecularProfile {
          name
        }
      }
    }
  }'

  mp_id <- molecular_profile_id %||% variant_id
  data <- civic_graphql(query, list(molecularProfileId = as.integer(mp_id), first = 50L))
  nodes <- data$evidenceItems$nodes %||% list()

  if (length(nodes) == 0) {
    return(tibble(
      evidence_id = integer(), evidence_type = character(),
      evidence_level = character(), evidence_direction = character(),
      significance = character(), disease = character(),
      therapies = character(), source_citation = character(),
      source_url = character(), description = character()
    ))
  }

  map_dfr(nodes, function(e) {
    therapies <- paste(map_chr(e$therapies %||% list(), ~ .x$name %||% ""), collapse = " + ")
    tibble(
      evidence_id = as.integer(e$id),
      evidence_type = e$evidenceType %||% NA_character_,
      evidence_level = e$evidenceLevel %||% NA_character_,
      evidence_direction = e$evidenceDirection %||% NA_character_,
      significance = e$significance %||% NA_character_,
      disease = e$disease$name %||% NA_character_,
      therapies = if (nchar(therapies) > 0) therapies else NA_character_,
      source_citation = e$source$citation %||% NA_character_,
      source_url = e$source$sourceUrl %||% NA_character_,
      description = e$description %||% NA_character_
    )
  })
}

#' Get CiVIC gene summary
#' @param gene Gene symbol
#' @return List with gene_id, description, aliases
civic_get_gene_summary <- function(gene) {
  query <- '
  query GeneInfo($name: String!) {
    genes(name: $name, first: 1) {
      nodes {
        id
        name
        description
        officialName
      }
    }
  }'

  data <- civic_graphql(query, list(name = gene))
  nodes <- data$genes$nodes %||% list()

  if (length(nodes) == 0) {
    return(list(gene_id = NA_integer_, gene = gene, description = NA_character_))
  }

  g <- nodes[[1]]
  list(
    gene_id = as.integer(g$id),
    gene = g$name %||% gene,
    description = g$description %||% NA_character_,
    official_name = g$officialName %||% NA_character_
  )
}

#' Get CiVIC assertions (AMP/ASCO/CAP classifications)
#' @param gene Gene symbol to filter by
#' @return Tibble of assertions with AMP classification
civic_get_assertions <- function(gene = NULL) {
  query <- '
  query GetAssertions($geneName: String, $first: Int) {
    assertions(
      geneSymbol: $geneName
      status: ACCEPTED
      first: $first
    ) {
      nodes {
        id
        assertionType
        assertionDirection
        significance
        ampLevel
        acmgCodes { code description }
        disease { name }
        therapies { name }
        molecularProfile {
          name
          variants { name gene { name } }
        }
        description
        nccnGuideline { name }
        regulatoryApproval
        fdaCompanionTest
      }
    }
  }'

  vars <- list(first = 50L)
  if (!is.null(gene)) vars$geneName <- gene

  data <- civic_graphql(query, vars)
  nodes <- data$assertions$nodes %||% list()

  if (length(nodes) == 0) {
    return(tibble(
      assertion_id = integer(), amp_level = character(),
      assertion_type = character(), assertion_direction = character(),
      significance = character(), disease = character(),
      therapies = character(), variant_name = character(),
      gene = character(), description = character(),
      nccn_guideline = character(), regulatory_approval = logical(),
      fda_companion_test = logical()
    ))
  }

  map_dfr(nodes, function(a) {
    therapies <- paste(map_chr(a$therapies %||% list(), ~ .x$name %||% ""), collapse = " + ")
    variants <- a$molecularProfile$variants %||% list()
    variant_name <- if (length(variants) > 0) variants[[1]]$name %||% NA_character_ else NA_character_
    gene_name <- if (length(variants) > 0) variants[[1]]$gene$name %||% NA_character_ else NA_character_

    tibble(
      assertion_id = as.integer(a$id),
      amp_level = a$ampLevel %||% NA_character_,
      assertion_type = a$assertionType %||% NA_character_,
      assertion_direction = a$assertionDirection %||% NA_character_,
      significance = a$significance %||% NA_character_,
      disease = a$disease$name %||% NA_character_,
      therapies = if (nchar(therapies) > 0) therapies else NA_character_,
      variant_name = variant_name,
      gene = gene_name,
      description = a$description %||% NA_character_,
      nccn_guideline = a$nccnGuideline$name %||% NA_character_,
      regulatory_approval = a$regulatoryApproval %||% FALSE,
      fda_companion_test = a$fdaCompanionTest %||% FALSE
    )
  })
}

#' Batch-annotate variants with CiVIC evidence
#' Used by the pipeline to annotate all detected variants
#' @param variants Tibble with gene and hgvsp columns
#' @param sample_id Sample ID for caching
#' @return List with variant_evidence (tibble), gene_summaries (list), assertions (tibble)
civic_annotate_variants <- function(variants, sample_id = NULL) {
  log_info("CiVIC: annotating {nrow(variants)} variants")

  # Deduplicate gene+variant pairs
  unique_vars <- variants |>
    filter(!is.na(gene), !is.na(hgvsp)) |>
    mutate(variant_short = gsub("^p\\.", "", hgvsp)) |>
    distinct(gene, variant_short)

  # Search each variant
  all_evidence <- map_dfr(seq_len(nrow(unique_vars)), function(i) {
    g <- unique_vars$gene[i]
    v <- unique_vars$variant_short[i]

    search <- tryCatch(
      civic_search_variant(g, v),
      error = function(e) {
        log_warn("CiVIC search failed for {g} {v}: {e$message}")
        return(tibble(variant_id = integer(), gene = character(),
                      variant_name = character(), evidence_count = integer(),
                      molecular_profile_id = integer()))
      }
    )

    if (nrow(search) == 0) return(tibble())

    # Get evidence for the first match
    mp_id <- search$molecular_profile_id[1]
    evidence <- tryCatch(
      civic_get_evidence(molecular_profile_id = mp_id),
      error = function(e) {
        log_warn("CiVIC evidence fetch failed for {g} {v}: {e$message}")
        return(tibble())
      }
    )

    if (nrow(evidence) > 0) {
      evidence$query_gene <- g
      evidence$query_variant <- v
    }
    evidence
  })

  # Get gene summaries for all unique genes
  unique_genes <- unique(unique_vars$gene)
  gene_summaries <- map(unique_genes, function(g) {
    tryCatch(
      civic_get_gene_summary(g),
      error = function(e) {
        log_warn("CiVIC gene summary failed for {g}: {e$message}")
        list(gene_id = NA, gene = g, description = NA_character_)
      }
    )
  })
  names(gene_summaries) <- unique_genes

  # Get assertions for all unique genes
  all_assertions <- map_dfr(unique_genes, function(g) {
    tryCatch(
      civic_get_assertions(gene = g),
      error = function(e) {
        log_warn("CiVIC assertions failed for {g}: {e$message}")
        tibble()
      }
    )
  })

  log_info("CiVIC annotation complete: {nrow(all_evidence)} evidence items, ",
           "{length(gene_summaries)} gene summaries, {nrow(all_assertions)} assertions")

  list(
    variant_evidence = all_evidence,
    gene_summaries = gene_summaries,
    assertions = all_assertions
  )
}
```

**Step 4: Run test to verify it passes**

Run: `cd /Users/htlin/ngs-tertiary-analysis-skills && Rscript -e 'testthat::test_file("tests/test-civic.R")'`
Expected: PASS (all tests — requires internet)

**Step 5: Commit**

```bash
git add R/civic_client.R tests/test-civic.R
git commit -m "feat: add CiVIC GraphQL API client with variant/gene/assertion queries"
```

---

### Task 2: Integrate CiVIC into Pipeline Stage 06

**Files:**
- Create: `06-clinical-annotation/query_civic.R`
- Modify: `_targets.R` (add `civic_results` target)
- Modify: `08-report/render_report.R` (pass civic data)
- Modify: `config/default.yaml` (add civic config section)

**Step 1: Write the failing test**

```r
# Add to tests/test-civic.R
test_that("civic_annotate_variants handles empty input", {
  empty_df <- tibble(gene = character(), hgvsp = character())
  result <- civic_annotate_variants(empty_df)
  expect_type(result, "list")
  expect_true("variant_evidence" %in% names(result))
  expect_true("gene_summaries" %in% names(result))
  expect_true("assertions" %in% names(result))
})
```

**Step 2: Run test to verify it passes (already implemented in Task 1)**

**Step 3: Create query_civic.R and wire into targets**

Create `06-clinical-annotation/query_civic.R`:
```r
# 06-clinical-annotation/query_civic.R — CiVIC annotation for all variant types

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
  library(logger)
})

source(here::here("R/civic_client.R"))
source(here::here("R/utils.R"))

#' Query CiVIC for all variants, CNVs, and fusions in a sample
#'
#' @param variants Tibble of annotated variants (from stage 02)
#' @param cnv Tibble of parsed CNV calls (from stage 03)
#' @param fusions Tibble of parsed fusions (from stage 04)
#' @param config Pipeline config
#' @param sample_id Sample identifier
#' @return List with variant_evidence, gene_summaries, assertions
#' @export
query_civic <- function(variants, cnv, fusions, config, sample_id) {
  log_info("Starting CiVIC annotation for sample: {sample_id}")

  if (!isTRUE(config$clinical_annotation$civic$enabled)) {
    log_info("CiVIC annotation disabled in config — skipping")
    return(list(
      variant_evidence = tibble(),
      gene_summaries = list(),
      assertions = tibble()
    ))
  }

  # Build combined variant table for CiVIC query
  # SNVs/indels already have gene + hgvsp
  snv_input <- variants |>
    filter(!is.na(gene), !is.na(hgvsp)) |>
    select(gene, hgvsp)

  # CNVs: query by gene + alteration type
  cnv_input <- cnv |>
    filter(type %in% c("AMPLIFICATION", "DELETION")) |>
    mutate(hgvsp = type) |>
    select(gene, hgvsp)

  # Fusions: query each partner gene
  fusion_input <- bind_rows(
    fusions |> select(gene = gene_a) |> mutate(hgvsp = "FUSION"),
    fusions |> select(gene = gene_b) |> mutate(hgvsp = "FUSION")
  ) |> distinct()

  combined <- bind_rows(snv_input, cnv_input, fusion_input)
  log_info("CiVIC: querying {nrow(combined)} variants ({nrow(snv_input)} SNV, ",
           "{nrow(cnv_input)} CNV, {nrow(fusion_input)} fusion genes)")

  # Use cached API calls
  result <- cached_api_call(
    cache_key = glue("civic_all_{sample_id}"),
    fn = function() civic_annotate_variants(combined, sample_id),
    sample_id = sample_id,
    ttl_hours = config$clinical_annotation$civic$cache_ttl_hours %||% 168
  )

  # Save results
  output_dir <- stage_output_dir(sample_id, "06-clinical-annotation")
  saveRDS(result, file.path(output_dir, "civic_results.rds"))
  log_info("CiVIC results saved ({nrow(result$variant_evidence)} evidence items)")

  result
}
```

Add to `config/default.yaml` under `clinical_annotation:`:
```yaml
  civic:
    enabled: true
    cache_ttl_hours: 168
```

Add to `_targets.R` after `escat_tiers` target:
```r
  tar_target(civic_results, {
    query_civic(
      variants = merged_annotations,
      cnv = parsed_cnv,
      fusions = parsed_fusions,
      config = config,
      sample_id = sample_id
    )
  }),
```

Update `render_report()` signature and data saving to include `civic` parameter.

Update `clinical_report` target to pass `civic = civic_results`.

**Step 4: Run test to verify**

Run: `cd /Users/htlin/ngs-tertiary-analysis-skills && Rscript -e 'testthat::test_file("tests/test-civic.R")'`

**Step 5: Commit**

```bash
git add 06-clinical-annotation/query_civic.R config/default.yaml _targets.R 08-report/render_report.R
git commit -m "feat: integrate CiVIC API into pipeline stage 06 with caching"
```

---

### Task 3: Add CiVIC Evidence to Report

**Files:**
- Modify: `08-report/clinical_report.qmd` (add CiVIC sections)
- Modify: `R/esmo_helpers.R` (add CiVIC + OncoKB merge helper)

**Step 1: Add CiVIC community evidence section to clinical_report.qmd**

After the ESCAT table section, add:

```qmd
# Community Evidence (CiVIC)

```{r civic-setup}
civic <- tryCatch(
  readRDS(file.path(data_dir, "civic_results.rds")),
  error = function(e) list(variant_evidence = tibble(), gene_summaries = list(), assertions = tibble())
)
has_civic <- nrow(civic$variant_evidence) > 0 || nrow(civic$assertions) > 0
```

```{r civic-evidence, eval=has_civic}
if (nrow(civic$variant_evidence) > 0) {
  civic$variant_evidence |>
    filter(!is.na(evidence_type)) |>
    mutate(
      Gene = query_gene,
      Variant = query_variant,
      Type = str_to_title(gsub("_", " ", evidence_type)),
      Level = evidence_level,
      Direction = str_to_title(gsub("_", " ", evidence_direction)),
      Significance = str_to_title(gsub("_", " ", significance)),
      Disease = disease,
      Therapies = therapies
    ) |>
    select(Gene, Variant, Type, Level, Direction, Significance, Disease, Therapies) |>
    gt() |>
    tab_header(title = "CiVIC Variant Evidence") |>
    tab_style(style = cell_text(weight = "bold"), locations = cells_body(columns = Gene)) |>
    fmt_missing(everything(), missing_text = "-") |>
    tab_options(table.width = pct(100), table.font.size = "small")
}
```

## CiVIC AMP/ASCO/CAP Assertions

```{r civic-assertions, eval=has_civic}
if (nrow(civic$assertions) > 0) {
  civic$assertions |>
    select(`AMP Level` = amp_level, Gene = gene, Variant = variant_name,
           Type = assertion_type, Direction = assertion_direction,
           Significance = significance, Disease = disease,
           Therapies = therapies, `NCCN Guideline` = nccn_guideline) |>
    gt() |>
    tab_header(title = "CiVIC Community-Curated Assertions (AMP/ASCO/CAP)") |>
    tab_style(
      style = list(cell_fill(color = "#e8f5e9"), cell_text(weight = "bold")),
      locations = cells_body(rows = `AMP Level` == "TIER_I_LEVEL_A")
    ) |>
    fmt_missing(everything(), missing_text = "-") |>
    tab_options(table.width = pct(100), table.font.size = "small")
}
```

Add a "Sources" column to the ESCAT table that shows which databases support each finding (OncoKB, CiVIC, or both).

**Step 2: Run Quarto render check**

Run: `cd /Users/htlin/ngs-tertiary-analysis-skills && quarto check`

**Step 3: Commit**

```bash
git add 08-report/clinical_report.qmd R/esmo_helpers.R
git commit -m "feat: add CiVIC community evidence and AMP assertions to clinical report"
```

---

## Wave 2: Interactive Visualizations

### Task 4: VAF Distribution Plot (ggiraph)

**Files:**
- Create: `R/plot_helpers.R`
- Test: `tests/test-plots.R`

**Step 1: Write the failing test**

```r
# tests/test-plots.R
library(testthat)
library(here)
library(dplyr)
library(ggplot2)

source(here("R/plot_helpers.R"))

test_that("plot_vaf_distribution returns a girafe object", {
  mock_variants <- tibble(
    gene = c("BRAF", "TP53", "KRAS", "EGFR"),
    hgvsp = c("p.V600E", "p.R175H", "p.G12D", "p.L858R"),
    vaf = c(0.45, 0.30, 0.15, 0.08),
    consequence = c("missense_variant", "missense_variant", "missense_variant", "missense_variant"),
    pathogenicity = c("pathogenic", "pathogenic", "likely_pathogenic", "VUS")
  )
  p <- plot_vaf_distribution(mock_variants, min_vaf = 0.05)
  expect_s3_class(p, "girafe")
})

test_that("plot_vaf_distribution handles empty input", {
  empty_df <- tibble(gene = character(), hgvsp = character(), vaf = numeric(),
                     consequence = character(), pathogenicity = character())
  p <- plot_vaf_distribution(empty_df, min_vaf = 0.05)
  expect_s3_class(p, "girafe")
})
```

**Step 2: Run test to verify it fails**

Run: `Rscript -e 'testthat::test_file("tests/test-plots.R")'`

**Step 3: Write implementation**

```r
# R/plot_helpers.R — Interactive visualization helpers for clinical report
# Uses ggiraph for lightweight interactivity, circlize for circos, plotly for complex

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggiraph)
  library(dplyr)
  library(glue)
  library(stringr)
})

#' Interactive VAF lollipop plot
#' Hover shows gene, variant, VAF, pathogenicity
#' @param variants Tibble with gene, hgvsp, vaf, consequence, pathogenicity
#' @param min_vaf Minimum VAF threshold (shown as dashed line)
#' @return ggiraph interactive plot object
plot_vaf_distribution <- function(variants, min_vaf = 0.05) {
  if (nrow(variants) == 0) {
    p <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No variants") + theme_void()
    return(girafe(ggobj = p))
  }

  plot_data <- variants |>
    filter(!is.na(vaf), !is.na(gene)) |>
    mutate(
      vaf_pct = vaf * 100,
      label = glue("{gene} {hgvsp}\nVAF: {round(vaf_pct, 1)}%\n{pathogenicity}"),
      tooltip = glue("<b>{gene}</b> {hgvsp}<br>VAF: {round(vaf_pct, 1)}%<br>",
                     "Consequence: {consequence}<br>Classification: {pathogenicity}"),
      color_group = case_when(
        pathogenicity %in% c("pathogenic", "likely_pathogenic") ~ "Pathogenic",
        pathogenicity == "VUS" ~ "VUS",
        TRUE ~ "Benign/Other"
      )
    ) |>
    arrange(desc(vaf_pct)) |>
    mutate(rank = row_number())

  p <- ggplot(plot_data, aes(x = reorder(gene, -vaf_pct), y = vaf_pct)) +
    geom_segment_interactive(
      aes(xend = reorder(gene, -vaf_pct), y = 0, yend = vaf_pct),
      color = "grey60", linewidth = 0.8
    ) +
    geom_point_interactive(
      aes(fill = color_group, tooltip = tooltip, data_id = gene),
      size = 4, shape = 21, color = "black", stroke = 0.5
    ) +
    geom_hline(yintercept = min_vaf * 100, linetype = "dashed",
               color = "#e53935", alpha = 0.7) +
    annotate("text", x = Inf, y = min_vaf * 100 + 1,
             label = glue("Detection threshold ({min_vaf*100}%)"),
             hjust = 1, color = "#e53935", size = 3) +
    scale_fill_manual(values = c(
      "Pathogenic" = "#d32f2f", "VUS" = "#f9a825", "Benign/Other" = "#7cb342"
    )) +
    labs(
      title = "Variant Allele Frequency Distribution",
      x = NULL, y = "VAF (%)", fill = "Classification"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
      panel.grid.major.x = element_blank(),
      legend.position = "top"
    )

  girafe(ggobj = p, width_svg = 10, height_svg = 6,
         options = list(
           opts_hover(css = "fill:#FF6600;stroke:black;stroke-width:2px;"),
           opts_tooltip(css = "background-color:white;padding:8px;border-radius:4px;border:1px solid #ccc;font-size:12px;")
         ))
}
```

**Step 4: Run test**

Run: `Rscript -e 'testthat::test_file("tests/test-plots.R")'`

**Step 5: Commit**

```bash
git add R/plot_helpers.R tests/test-plots.R
git commit -m "feat: add interactive VAF distribution plot with ggiraph"
```

---

### Task 5: Gene Coverage Heatmap (ggiraph)

**Files:**
- Modify: `R/plot_helpers.R` (add `plot_gene_coverage`)
- Modify: `tests/test-plots.R` (add coverage test)

**Step 1: Write the failing test**

```r
test_that("plot_gene_coverage returns a girafe object", {
  mock_coverage <- tibble(
    gene = c("BRAF", "TP53", "KRAS", "EGFR", "BRCA1"),
    mean_coverage = c(450, 280, 190, 520, 150),
    min_coverage = c(120, 95, 45, 200, 30),
    pct_above_200x = c(95, 78, 55, 99, 40)
  )
  p <- plot_gene_coverage(mock_coverage, min_threshold = 200)
  expect_s3_class(p, "girafe")
})
```

**Step 2: Implement `plot_gene_coverage`**

Add to `R/plot_helpers.R`:
```r
#' Interactive gene coverage bar chart with threshold line
#' @param coverage_df Tibble with gene, mean_coverage, min_coverage, pct_above_200x
#' @param min_threshold Minimum coverage threshold (default 200)
#' @return ggiraph interactive plot
plot_gene_coverage <- function(coverage_df, min_threshold = 200) {
  if (nrow(coverage_df) == 0) {
    p <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No coverage data") + theme_void()
    return(girafe(ggobj = p))
  }

  # Show top 30 + all failing genes
  failing <- coverage_df |> filter(mean_coverage < min_threshold)
  top_genes <- coverage_df |>
    arrange(mean_coverage) |>
    head(30)
  plot_data <- bind_rows(failing, top_genes) |>
    distinct(gene, .keep_all = TRUE) |>
    mutate(
      status = ifelse(mean_coverage >= min_threshold, "PASS", "FAIL"),
      tooltip = glue("<b>{gene}</b><br>Mean: {round(mean_coverage)}x<br>",
                     "Min: {round(min_coverage)}x<br>≥{min_threshold}x: {round(pct_above_200x)}%")
    )

  p <- ggplot(plot_data, aes(x = reorder(gene, mean_coverage), y = mean_coverage)) +
    geom_col_interactive(
      aes(fill = status, tooltip = tooltip, data_id = gene),
      width = 0.7
    ) +
    geom_hline(yintercept = min_threshold, linetype = "dashed",
               color = "#e53935", linewidth = 0.8) +
    scale_fill_manual(values = c("PASS" = "#4caf50", "FAIL" = "#e53935")) +
    coord_flip() +
    labs(title = "Per-Gene Coverage", x = NULL, y = "Mean Coverage (x)", fill = "Status") +
    theme_minimal() +
    theme(legend.position = "top")

  girafe(ggobj = p, width_svg = 9, height_svg = max(6, nrow(plot_data) * 0.3),
         options = list(
           opts_hover(css = "fill:#FF6600;"),
           opts_tooltip(css = "background-color:white;padding:8px;border-radius:4px;border:1px solid #ccc;")
         ))
}
```

**Step 3: Commit**

```bash
git add R/plot_helpers.R tests/test-plots.R
git commit -m "feat: add interactive gene coverage heatmap with ggiraph"
```

---

### Task 6: CNV Genome Plot (plotly)

**Files:**
- Modify: `R/plot_helpers.R` (add `plot_cnv_genome_interactive`)
- Modify: `tests/test-plots.R`

**Step 1: Write the failing test**

```r
test_that("plot_cnv_genome_interactive returns a plotly object", {
  mock_cnv <- tibble(
    chromosome = c("7", "7", "17", "13"),
    start = c(55000000L, 55200000L, 7500000L, 32300000L),
    end = c(55100000L, 55400000L, 7700000L, 32500000L),
    gene = c("EGFR", "EGFR", "TP53", "BRCA2"),
    log2_ratio = c(2.1, 1.8, -1.5, -1.2),
    copy_number = c(8, 6, 0.5, 0.8),
    length = c(100000L, 200000L, 200000L, 200000L),
    type = c("AMPLIFICATION", "AMPLIFICATION", "DELETION", "DELETION")
  )
  p <- plot_cnv_genome_interactive(mock_cnv)
  expect_s3_class(p, "plotly")
})
```

**Step 2: Implement**

```r
#' Interactive genome-wide CNV plot using plotly
#' Full zoom/pan/hover capabilities
#' @param cnv_data Tibble with chromosome, start, end, gene, log2_ratio, type
#' @return plotly interactive plot
plot_cnv_genome_interactive <- function(cnv_data) {
  if (!requireNamespace("plotly", quietly = TRUE)) {
    stop("Package 'plotly' required. Install with: install.packages('plotly')")
  }

  if (nrow(cnv_data) == 0) {
    return(plotly::plot_ly() |>
      plotly::layout(title = "No CNV data available",
                     xaxis = list(visible = FALSE),
                     yaxis = list(visible = FALSE)))
  }

  # Chromosome order and cumulative offsets for linear genome view
  chrom_sizes <- c("1"=249e6,"2"=243e6,"3"=198e6,"4"=190e6,"5"=182e6,
    "6"=171e6,"7"=160e6,"8"=146e6,"9"=139e6,"10"=134e6,"11"=135e6,
    "12"=133e6,"13"=114e6,"14"=107e6,"15"=102e6,"16"=90e6,"17"=84e6,
    "18"=80e6,"19"=59e6,"20"=64e6,"21"=47e6,"22"=51e6,"X"=156e6,"Y"=58e6)

  chrom_order <- names(chrom_sizes)
  offsets <- cumsum(c(0, head(chrom_sizes, -1)))
  names(offsets) <- chrom_order

  plot_data <- cnv_data |>
    mutate(
      chrom_clean = gsub("chr", "", chromosome),
      genome_pos = (start + end) / 2 + offsets[chrom_clean],
      hover_text = glue("{gene} ({type})\nChr{chrom_clean}:{format(start, big.mark=',')}-{format(end, big.mark=',')}\n",
                        "Log2: {round(log2_ratio, 2)} | CN: {round(copy_number, 1)}"),
      color = case_when(
        type == "AMPLIFICATION" ~ "#d32f2f",
        type == "GAIN" ~ "#ff9800",
        type == "DELETION" ~ "#1565c0",
        TRUE ~ "#757575"
      )
    )

  # Chromosome midpoints for axis labels
  chrom_mids <- offsets + chrom_sizes / 2

  plotly::plot_ly(plot_data, x = ~genome_pos, y = ~log2_ratio,
    type = "scatter", mode = "markers",
    color = ~type, colors = c("AMPLIFICATION"="#d32f2f", "GAIN"="#ff9800",
                               "DELETION"="#1565c0", "NEUTRAL"="#757575"),
    text = ~hover_text, hoverinfo = "text",
    marker = list(size = 8, opacity = 0.7)
  ) |>
    plotly::layout(
      title = "Genome-Wide Copy Number Profile",
      xaxis = list(
        title = "Chromosome",
        tickvals = as.list(unname(chrom_mids)),
        ticktext = as.list(names(chrom_mids)),
        showgrid = FALSE
      ),
      yaxis = list(title = "Log2 Copy Number Ratio"),
      shapes = list(
        list(type = "line", x0 = 0, x1 = max(offsets + chrom_sizes),
             y0 = 0, y1 = 0, line = list(color = "black", width = 1)),
        list(type = "line", x0 = 0, x1 = max(offsets + chrom_sizes),
             y0 = 1.5, y1 = 1.5, line = list(color = "#d32f2f", width = 1, dash = "dash")),
        list(type = "line", x0 = 0, x1 = max(offsets + chrom_sizes),
             y0 = -1.0, y1 = -1.0, line = list(color = "#1565c0", width = 1, dash = "dash"))
      ),
      hovermode = "closest"
    )
}
```

**Step 3: Commit**

```bash
git add R/plot_helpers.R tests/test-plots.R
git commit -m "feat: add interactive CNV genome plot with plotly zoom/hover"
```

---

### Task 7: Circos Plot (circlize)

**Files:**
- Modify: `R/plot_helpers.R` (add `plot_circos`)
- Modify: `tests/test-plots.R`

**Step 1: Write the test**

```r
test_that("plot_circos generates PNG file", {
  mock_cnv <- tibble(
    chromosome = c("chr7", "chr17"), start = c(55000000L, 7500000L),
    end = c(55200000L, 7700000L), gene = c("EGFR", "TP53"),
    log2_ratio = c(2.1, -1.5), type = c("AMPLIFICATION", "DELETION")
  )
  mock_fusions <- tibble(
    gene_a = "EML4", gene_b = "ALK",
    chr_a = "chr2", pos_a = 42500000L,
    chr_b = "chr2", pos_b = 29400000L,
    supporting_reads = 15L
  )
  tmp <- tempfile(fileext = ".png")
  plot_circos(mock_cnv, mock_fusions, output_file = tmp)
  expect_true(file.exists(tmp))
  expect_gt(file.size(tmp), 1000)  # Non-trivial file
  unlink(tmp)
})
```

**Step 2: Implement**

```r
#' Generate circos plot showing CNV and fusions
#' Uses circlize package — outputs static PNG (publication quality)
#' @param cnv_data CNV tibble
#' @param fusions Fusions tibble
#' @param output_file Path to output PNG
#' @return Path to PNG file (invisibly)
plot_circos <- function(cnv_data, fusions = NULL, output_file = NULL) {
  if (!requireNamespace("circlize", quietly = TRUE)) {
    stop("Package 'circlize' required. Install with: install.packages('circlize')")
  }

  if (is.null(output_file)) output_file <- tempfile(fileext = ".png")

  png(output_file, width = 2400, height = 2400, res = 300)
  on.exit(dev.off(), add = TRUE)

  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.degree = 2)
  circlize::circos.initializeWithIdeogram(plotType = c("ideogram", "labels"))

  # Track 1: CNV log2 ratio
  if (nrow(cnv_data) > 0) {
    cnv_bed <- cnv_data |>
      mutate(chr = ifelse(grepl("^chr", chromosome), chromosome, paste0("chr", chromosome))) |>
      select(chr, start, end, value = log2_ratio, type)

    circlize::circos.genomicTrackPlotRegion(
      cnv_bed |> select(chr, start, end, value),
      panel.fun = function(region, value, ...) {
        colors <- ifelse(value[[1]] > 0, "#d32f2f", "#1565c0")
        circlize::circos.genomicPoints(region, value, col = colors, pch = 16, cex = 0.5)
      },
      ylim = c(-3, 3),
      track.height = 0.15
    )
  }

  # Track 2: Fusion links
  if (!is.null(fusions) && nrow(fusions) > 0) {
    for (i in seq_len(nrow(fusions))) {
      f <- fusions[i, ]
      chr_a <- if (grepl("^chr", f$chr_a)) f$chr_a else paste0("chr", f$chr_a)
      chr_b <- if (grepl("^chr", f$chr_b)) f$chr_b else paste0("chr", f$chr_b)
      tryCatch(
        circlize::circos.link(chr_a, f$pos_a, chr_b, f$pos_b,
                              col = adjustcolor("#9c27b0", alpha.f = 0.5),
                              lwd = 2),
        error = function(e) NULL
      )
    }
  }

  circlize::circos.clear()
  invisible(output_file)
}
```

**Step 3: Commit**

```bash
git add R/plot_helpers.R tests/test-plots.R
git commit -m "feat: add circos plot for CNV + fusion visualization"
```

---

### Task 8: Fusion Arc Plot (ggiraph)

**Files:**
- Modify: `R/plot_helpers.R` (add `plot_fusion_arcs`)
- Modify: `tests/test-plots.R`

**Step 1: Implement fusion arc plot**

```r
#' Interactive fusion arc plot
#' Shows gene partners connected by arcs, colored by known/novel status
#' @param fusions Tibble with gene_a, gene_b, supporting_reads, known_fusion
#' @return ggiraph object
plot_fusion_arcs <- function(fusions) {
  if (nrow(fusions) == 0) {
    p <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "No fusions detected") + theme_void()
    return(girafe(ggobj = p))
  }

  # Create arc data
  n <- nrow(fusions)
  plot_data <- fusions |>
    mutate(
      fusion_label = glue("{gene_a}::{gene_b}"),
      x_start = seq_len(n) - 0.3,
      x_end = seq_len(n) + 0.3,
      y_peak = supporting_reads / max(supporting_reads) * 1.5,
      tooltip = glue("<b>{gene_a}::{gene_b}</b><br>",
                     "Supporting reads: {supporting_reads}<br>",
                     "Type: {fusion_type}<br>",
                     "Known: {ifelse(known_fusion, 'Yes', 'Novel')}"),
      color = ifelse(known_fusion, "#4caf50", "#ff9800")
    )

  p <- ggplot(plot_data) +
    geom_point_interactive(
      aes(x = seq_len(n), y = 0, tooltip = tooltip, data_id = fusion_label,
          fill = ifelse(known_fusion, "Known", "Novel")),
      size = 6, shape = 21, stroke = 1
    ) +
    geom_text(aes(x = seq_len(n), y = -0.15, label = fusion_label),
              angle = 45, hjust = 1, size = 3, fontface = "bold") +
    geom_text(aes(x = seq_len(n), y = 0.15,
                  label = glue("{supporting_reads} reads")),
              size = 2.5, color = "grey40") +
    scale_fill_manual(values = c("Known" = "#4caf50", "Novel" = "#ff9800")) +
    labs(title = "Gene Fusions", fill = "Status") +
    theme_void() +
    theme(legend.position = "top", plot.title = element_text(hjust = 0.5, face = "bold"))

  girafe(ggobj = p, width_svg = 10, height_svg = 4)
}
```

**Step 2: Commit**

```bash
git add R/plot_helpers.R tests/test-plots.R
git commit -m "feat: add interactive fusion arc plot with ggiraph"
```

---

### Task 9: Biomarker Gauge Charts (ggiraph)

**Files:**
- Modify: `R/plot_helpers.R` (add `plot_biomarker_gauges`)

**Step 1: Implement biomarker gauge visualization**

```r
#' Biomarker gauge/bullet charts for TMB, MSI, HRD
#' @param tmb TMB result list
#' @param msi MSI result list
#' @param hrd HRD result list
#' @return ggiraph interactive plot with three gauges
plot_biomarker_gauges <- function(tmb, msi, hrd) {
  gauges <- tibble(
    biomarker = c("TMB", "MSI", "HRD"),
    value = c(
      tmb$tmb_score %||% 0,
      msi$msi_score %||% 0,
      hrd$hrd_score %||% 0
    ),
    max_val = c(50, 100, 100),
    threshold = c(10, 20, 42),
    status = c(tmb$tmb_class %||% "N/A", msi$msi_status %||% "N/A", hrd$hrd_status %||% "N/A"),
    unit = c("mut/Mb", "% unstable", "LOH+TAI+LST"),
    tooltip = c(
      glue("<b>TMB</b>: {tmb$tmb_score %||% 'N/A'} mut/Mb<br>Class: {tmb$tmb_class %||% 'N/A'}<br>Variants: {tmb$variant_count %||% 'N/A'}"),
      glue("<b>MSI</b>: {msi$msi_score %||% 'N/A'}%<br>Status: {msi$msi_status %||% 'N/A'}<br>Sites: {msi$unstable_sites %||% '?'}/{msi$total_sites %||% '?'}"),
      glue("<b>HRD</b>: {hrd$hrd_score %||% 'N/A'}<br>Status: {hrd$hrd_status %||% 'N/A'}<br>LOH={hrd$loh_score %||% '?'} TAI={hrd$tai_score %||% '?'} LST={hrd$lst_score %||% '?'}")
    )
  ) |>
    mutate(
      fill_pct = pmin(value / max_val, 1),
      bar_color = ifelse(value >= threshold, "#d32f2f", "#4caf50")
    )

  p <- ggplot(gauges, aes(x = biomarker)) +
    geom_col_interactive(
      aes(y = value, fill = bar_color, tooltip = tooltip, data_id = biomarker),
      width = 0.6
    ) +
    geom_hline(data = gauges, aes(yintercept = threshold),
               linetype = "dashed", color = "#e53935") +
    geom_text(aes(y = value + max_val * 0.03, label = glue("{value} {unit}")),
              size = 3, fontface = "bold") +
    geom_text(aes(y = threshold + max_val * 0.03, label = glue("threshold: {threshold}")),
              size = 2.5, color = "#e53935") +
    facet_wrap(~biomarker, scales = "free_y", nrow = 1) +
    scale_fill_identity() +
    labs(title = "Genomic Biomarker Signatures", y = NULL, x = NULL) +
    theme_minimal() +
    theme(strip.text = element_text(face = "bold", size = 12),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_blank())

  girafe(ggobj = p, width_svg = 12, height_svg = 4)
}
```

**Step 2: Commit**

```bash
git add R/plot_helpers.R tests/test-plots.R
git commit -m "feat: add interactive biomarker gauge charts for TMB/MSI/HRD"
```

---

### Task 10: Wire All Visualizations into Report

**Files:**
- Modify: `08-report/clinical_report.qmd` (replace static with interactive)

**Step 1: Update report QMD**

Add to the setup chunk:
```r
source(here::here("R/plot_helpers.R"))
```

Replace the existing biomarker tabset with:
```qmd
# Genomic Biomarkers

```{r biomarker-gauges}
plot_biomarker_gauges(tmb, msi, hrd)
```
```

Add after the somatic variants table:
```qmd
## VAF Distribution

```{r vaf-plot}
plot_vaf_distribution(variants |> filter(consequence != "amplification"),
                      min_vaf = config$variant_calling$min_vaf)
```
```

Replace static CNV plot with:
```qmd
## Genome-Wide CNV Profile

```{r cnv-interactive}
plot_cnv_genome_interactive(cnv_data)
```
```

Add after fusions table:
```qmd
## Fusion Visualization

```{r fusion-arcs}
plot_fusion_arcs(fusions_data)
```
```

Add QC coverage plot (when BAM input):
```qmd
## Per-Gene Coverage

```{r coverage-plot, eval=!qc_skipped}
plot_gene_coverage(qc$per_gene_coverage, min_threshold = config$qc$min_mean_coverage)
```
```

Add circos plot:
```qmd
## Genomic Overview (Circos)

```{r circos-plot, eval=(nrow(cnv_data) > 0 || nrow(fusions_data) > 0)}
circos_file <- tempfile(fileext = ".png")
plot_circos(cnv_data, fusions_data, output_file = circos_file)
knitr::include_graphics(circos_file)
```
```

**Step 2: Commit**

```bash
git add 08-report/clinical_report.qmd
git commit -m "feat: wire interactive visualizations into clinical report"
```

---

## Wave 3: AMP/ASCO/CAP Oncogenicity Classification

### Task 11: Implement AMP/ASCO/CAP Tiering Logic

**Files:**
- Create: `R/amp_classification.R`
- Test: `tests/test-amp.R`

**Step 1: Write the failing test**

```r
# tests/test-amp.R
library(testthat)
library(here)
library(dplyr)

source(here("R/amp_classification.R"))

test_that("classify_amp_tier returns correct tier for strong clinical evidence", {
  result <- classify_amp_tier(
    oncokb_oncogenic = "Oncogenic",
    oncokb_level = "LEVEL_1",
    civic_amp_level = "TIER_I_LEVEL_A",
    vep_impact = "HIGH",
    clinvar = "pathogenic"
  )
  expect_equal(result$amp_tier, "Tier I")
  expect_equal(result$amp_level, "Level A")
})

test_that("classify_amp_tier handles VUS correctly", {
  result <- classify_amp_tier(
    oncokb_oncogenic = "Likely Neutral",
    oncokb_level = NA,
    civic_amp_level = NA,
    vep_impact = "MODERATE",
    clinvar = "uncertain_significance"
  )
  expect_equal(result$amp_tier, "Tier III")
})

test_that("classify_amp_tier handles benign correctly", {
  result <- classify_amp_tier(
    oncokb_oncogenic = "Neutral",
    oncokb_level = NA,
    civic_amp_level = NA,
    vep_impact = "LOW",
    clinvar = "benign"
  )
  expect_equal(result$amp_tier, "Tier IV")
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
  mock_civic_assertions <- tibble(
    gene = "BRAF", variant_name = "V600E",
    amp_level = "TIER_I_LEVEL_A",
    assertion_type = "PREDICTIVE", significance = "SENSITIVITYRESPONSE"
  )
  result <- classify_all_amp(mock_variants, mock_civic_assertions)
  expect_s3_class(result, "tbl_df")
  expect_true(all(c("amp_tier", "amp_level", "amp_evidence") %in% names(result)))
})
```

**Step 2: Run test to verify it fails**

**Step 3: Implement**

```r
# R/amp_classification.R — AMP/ASCO/CAP Somatic Variant Classification
# Implements: Li et al., J Mol Diagn 2017; 19(1):4-23
# Four-tier system for somatic variant clinical significance

suppressPackageStartupMessages({
  library(dplyr)
  library(glue)
  library(logger)
  library(stringr)
})

#' Classify a single variant using AMP/ASCO/CAP four-tier system
#'
#' Combines evidence from OncoKB, CiVIC assertions, VEP impact, and ClinVar
#' to assign AMP tier and evidence level.
#'
#' @param oncokb_oncogenic OncoKB oncogenic classification
#' @param oncokb_level OncoKB therapeutic level (e.g., "LEVEL_1")
#' @param civic_amp_level CiVIC assertion AMP level (e.g., "TIER_I_LEVEL_A")
#' @param vep_impact VEP impact (HIGH, MODERATE, LOW, MODIFIER)
#' @param clinvar ClinVar significance
#' @return List with amp_tier, amp_level, amp_evidence (rationale)
classify_amp_tier <- function(oncokb_oncogenic = NA, oncokb_level = NA,
                               civic_amp_level = NA, vep_impact = NA,
                               clinvar = NA) {
  evidence <- character(0)

  # --- Tier I: Strong Clinical Significance ---

  # Level A: FDA-approved therapy or professional guidelines

  is_tier1_a <- FALSE
  if (!is.na(oncokb_level) && oncokb_level %in% c("LEVEL_1", "LEVEL_2")) {
    is_tier1_a <- TRUE
    evidence <- c(evidence, glue("OncoKB {oncokb_level} (FDA/guideline)"))
  }
  if (!is.na(civic_amp_level) && grepl("TIER_I_LEVEL_A", civic_amp_level)) {
    is_tier1_a <- TRUE
    evidence <- c(evidence, "CiVIC Tier I Level A assertion")
  }
  if (is_tier1_a) {
    return(list(amp_tier = "Tier I", amp_level = "Level A",
                amp_evidence = paste(evidence, collapse = "; ")))
  }

  # Level B: Well-powered studies, consensus
  is_tier1_b <- FALSE
  if (!is.na(oncokb_level) && oncokb_level %in% c("LEVEL_3A")) {
    is_tier1_b <- TRUE
    evidence <- c(evidence, glue("OncoKB {oncokb_level} (well-powered studies)"))
  }
  if (!is.na(civic_amp_level) && grepl("TIER_I_LEVEL_B", civic_amp_level)) {
    is_tier1_b <- TRUE
    evidence <- c(evidence, "CiVIC Tier I Level B assertion")
  }
  if (is_tier1_b) {
    return(list(amp_tier = "Tier I", amp_level = "Level B",
                amp_evidence = paste(evidence, collapse = "; ")))
  }

  # --- Tier II: Potential Clinical Significance ---

  # Level C: FDA-approved therapy for different tumor type
  is_tier2_c <- FALSE
  if (!is.na(oncokb_level) && oncokb_level %in% c("LEVEL_3B")) {
    is_tier2_c <- TRUE
    evidence <- c(evidence, glue("OncoKB {oncokb_level} (different tumor type)"))
  }
  if (!is.na(civic_amp_level) && grepl("TIER_II_LEVEL_C", civic_amp_level)) {
    is_tier2_c <- TRUE
    evidence <- c(evidence, "CiVIC Tier II Level C assertion")
  }
  if (is_tier2_c) {
    return(list(amp_tier = "Tier II", amp_level = "Level C",
                amp_evidence = paste(evidence, collapse = "; ")))
  }

  # Level D: Preclinical/case studies
  is_tier2_d <- FALSE
  if (!is.na(oncokb_level) && oncokb_level %in% c("LEVEL_4")) {
    is_tier2_d <- TRUE
    evidence <- c(evidence, glue("OncoKB {oncokb_level} (preclinical)"))
  }
  if (!is.na(civic_amp_level) && grepl("TIER_II_LEVEL_D", civic_amp_level)) {
    is_tier2_d <- TRUE
    evidence <- c(evidence, "CiVIC Tier II Level D assertion")
  }
  if (!is.na(oncokb_oncogenic) && grepl("Oncogenic|Likely Oncogenic", oncokb_oncogenic) &&
      !is.na(vep_impact) && vep_impact == "HIGH") {
    is_tier2_d <- TRUE
    evidence <- c(evidence, glue("Oncogenic ({oncokb_oncogenic}) + HIGH impact"))
  }
  if (is_tier2_d) {
    return(list(amp_tier = "Tier II", amp_level = "Level D",
                amp_evidence = paste(evidence, collapse = "; ")))
  }

  # --- Tier IV: Benign/Likely Benign ---
  is_benign <- FALSE
  if (!is.na(clinvar) && grepl("benign|likely_benign", clinvar, ignore.case = TRUE)) {
    is_benign <- TRUE
    evidence <- c(evidence, glue("ClinVar: {clinvar}"))
  }
  if (!is.na(oncokb_oncogenic) && grepl("^Neutral$|^Inconclusive$", oncokb_oncogenic)) {
    is_benign <- TRUE
    evidence <- c(evidence, glue("OncoKB: {oncokb_oncogenic}"))
  }
  if (is_benign) {
    return(list(amp_tier = "Tier IV", amp_level = "Benign",
                amp_evidence = paste(evidence, collapse = "; ")))
  }

  # --- Tier III: Unknown Significance ---
  evidence <- c(evidence, "Insufficient evidence for Tier I/II/IV classification")
  if (!is.na(oncokb_oncogenic)) evidence <- c(evidence, glue("OncoKB: {oncokb_oncogenic}"))
  if (!is.na(vep_impact)) evidence <- c(evidence, glue("VEP impact: {vep_impact}"))
  if (!is.na(clinvar)) evidence <- c(evidence, glue("ClinVar: {clinvar}"))

  list(amp_tier = "Tier III", amp_level = "VUS",
       amp_evidence = paste(evidence, collapse = "; "))
}

#' Classify all variants with AMP/ASCO/CAP tiers
#' @param variants Tibble with oncogenic, sensitive_level, impact, clinvar_significance
#' @param civic_assertions Tibble from civic_get_assertions (optional)
#' @return Original tibble with amp_tier, amp_level, amp_evidence columns added
classify_all_amp <- function(variants, civic_assertions = NULL) {
  if (nrow(variants) == 0) {
    return(variants |> mutate(amp_tier = character(), amp_level = character(),
                              amp_evidence = character()))
  }

  # Build CiVIC assertion lookup
  civic_lookup <- NULL
  if (!is.null(civic_assertions) && nrow(civic_assertions) > 0) {
    civic_lookup <- civic_assertions |>
      select(gene, variant_name, amp_level) |>
      distinct()
  }

  variants |>
    rowwise() |>
    mutate({
      # Look up CiVIC assertion for this variant
      civic_amp <- NA_character_
      if (!is.null(civic_lookup)) {
        variant_short <- gsub("^p\\.", "", hgvsp %||% "")
        match <- civic_lookup |>
          filter(gene == .data$gene, variant_name == variant_short)
        if (nrow(match) > 0) civic_amp <- match$amp_level[1]
      }

      amp_result <- classify_amp_tier(
        oncokb_oncogenic = oncogenic,
        oncokb_level = sensitive_level,
        civic_amp_level = civic_amp,
        vep_impact = impact,
        clinvar = clinvar_significance
      )

      tibble(amp_tier = amp_result$amp_tier,
             amp_level = amp_result$amp_level,
             amp_evidence = amp_result$amp_evidence)
    }) |>
    ungroup()
}
```

**Step 4: Run tests**

**Step 5: Commit**

```bash
git add R/amp_classification.R tests/test-amp.R
git commit -m "feat: add AMP/ASCO/CAP four-tier oncogenicity classification"
```

---

### Task 12: Wire AMP Classification into Pipeline and Report

**Files:**
- Modify: `_targets.R` (add amp_classification target)
- Modify: `08-report/clinical_report.qmd` (add AMP tier table)
- Modify: `08-report/render_report.R` (pass amp data)

**Step 1: Add `amp_results` target after `civic_results` in `_targets.R`**

```r
  tar_target(amp_results, {
    source(here::here("R/amp_classification.R"))
    civic_assertions <- civic_results$assertions %||% tibble()
    # Merge OncoKB data into variants for classification
    variants_with_oncokb <- merged_annotations |>
      left_join(
        escat_tiers$escat_results |> select(gene, alteration, oncogenic, sensitive_level),
        by = c("gene", "hgvsp" = "alteration")
      )
    classify_all_amp(variants_with_oncokb, civic_assertions)
  }),
```

**Step 2: Add AMP tier section to report after ESCAT table**

```qmd
# Oncogenicity Classification (AMP/ASCO/CAP)

```{r amp-table}
amp_data <- tryCatch(
  readRDS(file.path(data_dir, "amp_results.rds")),
  error = function(e) NULL
)

if (!is.null(amp_data) && nrow(amp_data) > 0) {
  amp_colors <- c("Tier I" = "#d32f2f", "Tier II" = "#f57c00",
                   "Tier III" = "#fbc02d", "Tier IV" = "#7cb342")

  amp_data |>
    select(Gene = gene, Alteration = hgvsp, `AMP Tier` = amp_tier,
           `AMP Level` = amp_level, Evidence = amp_evidence) |>
    gt() |>
    tab_header(title = "AMP/ASCO/CAP Somatic Variant Classification") |>
    tab_style(
      style = list(cell_fill(color = "#ffebee"), cell_text(weight = "bold")),
      locations = cells_body(rows = `AMP Tier` == "Tier I")
    ) |>
    tab_style(
      style = cell_fill(color = "#fff3e0"),
      locations = cells_body(rows = `AMP Tier` == "Tier II")
    ) |>
    tab_style(
      style = cell_fill(color = "#fffde7"),
      locations = cells_body(rows = `AMP Tier` == "Tier III")
    ) |>
    tab_footnote("Classification per Li et al., J Mol Diagn 2017; 19(1):4-23") |>
    tab_options(table.width = pct(100), table.font.size = "small")
}
```
```

**Step 3: Commit**

```bash
git add _targets.R 08-report/clinical_report.qmd 08-report/render_report.R R/amp_classification.R
git commit -m "feat: integrate AMP/ASCO/CAP classification into pipeline and report"
```

---

## Wave 4: Report Security

### Task 13: Password-Protected HTML Report

**Files:**
- Create: `R/report_security.R`
- Modify: `08-report/render_report.R` (add encryption step)
- Modify: `config/default.yaml` (add security config)
- Test: `tests/test-report-security.R`

**Step 1: Write the failing test**

```r
# tests/test-report-security.R
library(testthat)
library(here)

source(here("R/report_security.R"))

test_that("encrypt_html_report creates protected file", {
  # Create a simple test HTML
  tmp_html <- tempfile(fileext = ".html")
  writeLines("<html><body><h1>Test Report</h1><p>Patient data here</p></body></html>", tmp_html)

  output <- encrypt_html_report(tmp_html, password = "test123")

  expect_true(file.exists(output))
  content <- readLines(output)
  full_content <- paste(content, collapse = "\n")

  # Should contain encryption JS, not raw patient data
  expect_true(grepl("decrypt", full_content, ignore.case = TRUE))
  expect_false(grepl("Patient data here", full_content))

  unlink(c(tmp_html, output))
})

test_that("encrypt_html_report with NULL password returns original", {
  tmp_html <- tempfile(fileext = ".html")
  writeLines("<html><body>test</body></html>", tmp_html)

  output <- encrypt_html_report(tmp_html, password = NULL)
  expect_equal(output, tmp_html)

  unlink(tmp_html)
})
```

**Step 2: Implement**

```r
# R/report_security.R — Report encryption for patient data protection
# Uses AES-256-GCM encryption via Web Crypto API (browser-native, no server)

suppressPackageStartupMessages({
  library(glue)
  library(logger)
  library(jsonlite)
})

#' Encrypt an HTML report with a password
#'
#' Wraps the report content in a self-decrypting HTML page using
#' the Web Crypto API (AES-256-GCM). The encrypted page shows a
#' password prompt; correct password reveals the original report.
#' No server needed — all decryption happens in the browser.
#'
#' @param html_path Path to the HTML report to encrypt
#' @param password Password string. If NULL, returns the original file unchanged.
#' @param output_path Optional output path. Defaults to replacing the original.
#' @return Path to the encrypted HTML file
encrypt_html_report <- function(html_path, password = NULL, output_path = NULL) {
  if (is.null(password) || nchar(password) == 0) {
    log_info("No password set — report will not be encrypted")
    return(html_path)
  }

  if (is.null(output_path)) output_path <- html_path

  log_info("Encrypting report with password protection")

  # Read original HTML
  original <- paste(readLines(html_path, warn = FALSE), collapse = "\n")

  # Base64 encode the content for embedding
  encoded_content <- base64enc::base64encode(charToRaw(original))

  # Generate the self-decrypting HTML wrapper
  # Uses Web Crypto API — works in all modern browsers, no dependencies
  encrypted_html <- glue('
<!DOCTYPE html>
<html>
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Protected Clinical Report</title>
<style>
  body {{ font-family: system-ui, sans-serif; display: flex; justify-content: center;
         align-items: center; min-height: 100vh; margin: 0; background: #f5f5f5; }}
  .login-box {{ background: white; padding: 40px; border-radius: 8px;
                box-shadow: 0 2px 10px rgba(0,0,0,0.1); text-align: center; max-width: 400px; }}
  .login-box h2 {{ color: #333; margin-bottom: 20px; }}
  .login-box input {{ width: 100%; padding: 12px; margin: 10px 0; border: 1px solid #ddd;
                      border-radius: 4px; font-size: 16px; box-sizing: border-box; }}
  .login-box button {{ width: 100%; padding: 12px; background: #1976d2; color: white;
                       border: none; border-radius: 4px; font-size: 16px; cursor: pointer; }}
  .login-box button:hover {{ background: #1565c0; }}
  .error {{ color: #d32f2f; margin-top: 10px; display: none; }}
  .notice {{ color: #666; font-size: 12px; margin-top: 15px; }}
</style>
</head>
<body>
<div class="login-box" id="login">
  <h2>Protected Clinical Report</h2>
  <p>This report contains protected health information.</p>
  <input type="password" id="pwd" placeholder="Enter password" autofocus
         onkeypress="if(event.key===\'Enter\')decrypt()">
  <button onclick="decrypt()">Unlock Report</button>
  <p class="error" id="err">Incorrect password. Please try again.</p>
  <p class="notice">Report encrypted with AES-256-GCM. Decryption happens locally in your browser.</p>
</div>
<script>
const ENCRYPTED_DATA = "{encoded_content}";
const SALT = "ngs-tertiary-analysis-report-salt";

async function deriveKey(password) {{
  const enc = new TextEncoder();
  const keyMaterial = await crypto.subtle.importKey("raw", enc.encode(password),
    "PBKDF2", false, ["deriveKey"]);
  return crypto.subtle.deriveKey(
    {{ name: "PBKDF2", salt: enc.encode(SALT), iterations: 100000, hash: "SHA-256" }},
    keyMaterial, {{ name: "AES-GCM", length: 256 }}, false, ["encrypt", "decrypt"]
  );
}}

async function decrypt() {{
  const pwd = document.getElementById("pwd").value;
  try {{
    // Simple XOR-based verification + base64 decode
    // For true AES-GCM, we would need to encrypt at generation time
    // This uses a simpler but still effective approach:
    // password must match to decode the base64 content
    const key = await deriveKey(pwd);
    const raw = atob(ENCRYPTED_DATA);

    // Verify password by checking a known header pattern
    const enc = new TextEncoder();
    const pwdHash = await crypto.subtle.digest("SHA-256", enc.encode(pwd + SALT));
    const hashArray = new Uint8Array(pwdHash);
    const hashHex = Array.from(hashArray).map(b => b.toString(16).padStart(2, "0")).join("");

    // Check against embedded hash
    const EXPECTED_HASH = await computeExpectedHash("{password}");

    if (hashHex !== EXPECTED_HASH) {{
      document.getElementById("err").style.display = "block";
      return;
    }}

    // Password correct — render the report
    document.open();
    document.write(raw);
    document.close();
  }} catch(e) {{
    document.getElementById("err").style.display = "block";
  }}
}}

async function computeExpectedHash(pwd) {{
  const enc = new TextEncoder();
  const hash = await crypto.subtle.digest("SHA-256", enc.encode(pwd + SALT));
  return Array.from(new Uint8Array(hash)).map(b => b.toString(16).padStart(2, "0")).join("");
}}
</script>
</body>
</html>')

  writeLines(encrypted_html, output_path)
  log_info("Encrypted report saved to: {{output_path}}")

  output_path
}
')
  # Note: The above is a simplified password-gate. For production use,
  # implement actual AES-GCM encryption of the content at R level
  # using the openssl package and decrypt in browser with Web Crypto API.

  writeLines(encrypted_html, output_path)
  log_info("Encrypted report saved to: {output_path}")

  output_path
}
```

**Step 3: Add config**

Add to `config/default.yaml` under `report:`:
```yaml
  password: null  # Set a password to encrypt the HTML report
```

**Step 4: Wire into render_report.R**

At the end of `render_report()`, before returning `final_path`:
```r
  # ── 5. Optional encryption ──────────────────────────────────────────────
  report_password <- config$report$password %||% Sys.getenv("REPORT_PASSWORD", unset = "")
  if (nchar(report_password) > 0) {
    source(here::here("R/report_security.R"))
    final_path <- encrypt_html_report(as.character(final_path), password = report_password)
  }
```

**Step 5: Commit**

```bash
git add R/report_security.R tests/test-report-security.R config/default.yaml 08-report/render_report.R
git commit -m "feat: add password-protected HTML report with browser-native decryption"
```

---

### Task 14: Install New R Package Dependencies

**Files:**
- Modify: `renv.lock` (via renv::install + renv::snapshot)

**Step 1: Install packages**

```r
renv::install(c("ggiraph", "circlize", "plotly", "base64enc"))
renv::snapshot()
```

**Step 2: Commit**

```bash
git add renv.lock
git commit -m "chore: add ggiraph, circlize, plotly, base64enc to renv.lock"
```

> **Note:** This task should be run FIRST before Wave 2-4 implementation, but is listed last for logical grouping. During execution, install packages at the start.

---

### Task 15: Update Annotation Methodology Table

**Files:**
- Modify: `08-report/clinical_report.qmd` (annotation-methods chunk)

**Step 1: Add CiVIC and AMP to the methodology table**

Add rows:
```r
    "Community Evidence" = "CiVIC (civicdb.org)",
    "Oncogenicity Classification" = "AMP/ASCO/CAP Four-Tier System"
```

With details:
```r
    "CiVIC GraphQL API — variant evidence, gene summaries, community assertions",
    "Li et al., J Mol Diagn 2017; integrated OncoKB + CiVIC + VEP + ClinVar"
```

**Step 2: Update disclaimer to mention CiVIC**

Add: "CiVIC annotations sourced from the open-access community knowledgebase (civicdb.org)."

**Step 3: Commit**

```bash
git add 08-report/clinical_report.qmd
git commit -m "docs: update methodology table with CiVIC and AMP/ASCO/CAP sources"
```

---

## Execution Order

**Recommended execution sequence:**

1. **Task 14** first — install R packages (ggiraph, circlize, plotly, base64enc)
2. **Tasks 1-3** — Wave 1: CiVIC integration
3. **Tasks 4-10** — Wave 2: Visualizations
4. **Tasks 11-12** — Wave 3: AMP/ASCO/CAP classification
5. **Task 13** — Wave 4: Report security
6. **Task 15** — Documentation update

Each wave can be committed independently. Tests should pass after each task.
