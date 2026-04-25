#!/usr/bin/env Rscript
# benchmarks/agent/03_prepare_stress_cohort.R
#
# Phase 2: build a cross-tumor STRESS cohort for the evidence_reconciler.
#
# Strategy: harvest CiVIC's ACCEPTED + AMP-level + single-gene assertions whose
# genes are OUTSIDE our existing 100-COSMIC cohort. This gives a naturally-
# diverse cohort drawn from CiVIC's curation (not constructed by us), exercises
# variants the deterministic baseline has not seen, and queries OncoKB live for
# each (cached). Outputs go to agents/inputs/evidence_reconciler/ with case_id
# prefix "stress_".

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(stringr)
  library(jsonlite)
  library(logger)
  library(glue)
})

source(here::here("R/utils.R"))
source(here::here("R/api_clients.R"))
source(here::here("R/civic_client.R"))
source(here::here("R/amp_classification.R"))
source(here::here("agents/io.R"))
source(here::here("agents/evidence_reconciler/format_input.R"))

load_env()
setup_logging("INFO")

# ── Pull all CiVIC assertions (paginated) ────────────────────────────────────

fetch_all_civic_assertions <- function(page_size = 50L, max_pages = 20L) {
  all_nodes <- list()
  cursor <- NULL
  page <- 0L
  repeat {
    page <- page + 1L
    after_clause <- if (is.null(cursor)) "" else paste0(', after: "', cursor, '"')
    q <- sprintf('
      query GetAssertions {
        assertions(first: %d%s) {
          nodes {
            id ampLevel status
            disease { name doid }
            molecularProfile {
              name
              variants { feature { name } }
            }
          }
          pageInfo { hasNextPage endCursor }
        }
      }', page_size, after_clause)
    res <- civic_graphql(q)
    nodes <- res$data$assertions$nodes %||% list()
    all_nodes <- c(all_nodes, nodes)
    pi <- res$data$assertions$pageInfo
    if (!isTRUE(pi$hasNextPage) || page >= max_pages) break
    cursor <- pi$endCursor
  }
  log_info("Fetched {length(all_nodes)} CiVIC assertions across {page} page(s)")
  all_nodes
}

raw_nodes <- fetch_all_civic_assertions()

civic_long <- map_dfr(raw_nodes, function(n) {
  if (toupper(n$status %||% "") != "ACCEPTED") return(NULL)
  if (is.null(n$ampLevel)) return(NULL)
  vs <- n$molecularProfile$variants %||% list()
  feats <- unique(unlist(map(vs, ~ .x$feature$name %||% NA_character_)))
  feats <- feats[!is.na(feats) & !grepl("::", feats, fixed = TRUE)]
  if (length(feats) == 0) return(NULL)
  mp_name <- n$molecularProfile$name %||% NA_character_
  map_dfr(feats, function(g) {
    tibble(
      gene         = g,
      variant_key  = str_remove(str_remove(mp_name, paste0("^", g, "\\s+")), "^p\\."),
      mp_name      = mp_name,
      amp_level    = n$ampLevel,
      disease      = n$disease$name %||% NA_character_,
      assertion_id = as.integer(n$id %||% NA_integer_)
    )
  })
})

log_info("ACCEPTED + has-ampLevel + single-gene assertions: {nrow(civic_long)}")

# ── Filter to genes outside existing cohort + map disease → OncoTree ─────────

cohort_path <- here::here("benchmarks/results/evidence_concordance.csv")
cohort_genes <- unique(read.csv(cohort_path, stringsAsFactors = FALSE)$gene)
log_info("Existing cohort genes: {length(cohort_genes)}")

# Map free-text CiVIC disease → OncoTree codes (best-effort manual map for the
# main cancer types). Only diseases we can reliably code go into the stress
# cohort; the rest get logged and dropped.
disease_to_oncotree <- function(d) {
  d_lower <- tolower(d %||% NA_character_)
  case_when(
    is.na(d_lower) ~ NA_character_,
    grepl("acute myeloid leukemia|aml", d_lower) ~ "AML",
    grepl("chronic myeloid leukemia|cml", d_lower) ~ "CML",
    grepl("acute lymphoblastic|acute lymphocytic|all", d_lower) ~ "ALL",
    grepl("chronic lymphocytic|cll", d_lower) ~ "CLL",
    grepl("chronic eosinophilic|hyper.?eosinophilic", d_lower) ~ "MPN",
    grepl("myelodysplastic|mds", d_lower) ~ "MDS",
    grepl("myeloproliferative|polycythemia|essential thrombocythemia", d_lower) ~ "MPN",
    grepl("langerhans cell histiocytosis|hairy cell leukemia|hcl", d_lower) ~ "HCL",
    grepl("erdheim.?chester", d_lower) ~ "ECD",
    grepl("non.?small cell lung|nsclc|lung adenocarcinoma|lung squamous", d_lower) ~ "NSCLC",
    grepl("small cell lung|sclc", d_lower) ~ "SCLC",
    grepl("colorectal|colon|rectal|coad", d_lower) ~ "COAD",
    grepl("melanoma|mel", d_lower) ~ "MEL",
    grepl("breast", d_lower) ~ "BRCA",
    grepl("ovarian|ovary|fallopian", d_lower) ~ "OVT",
    grepl("pancreatic|pancreas|pdac", d_lower) ~ "PAAD",
    grepl("prostate", d_lower) ~ "PRAD",
    grepl("medullary thyroid|mtc", d_lower) ~ "MTC",
    grepl("papillary thyroid|ptc", d_lower) ~ "THPA",
    grepl("anaplastic thyroid|atc", d_lower) ~ "THAP",
    grepl("thyroid", d_lower) ~ "THCA",
    grepl("cholangiocarcinoma|bile duct", d_lower) ~ "CHOL",
    grepl("hepatocellular|hcc", d_lower) ~ "HCC",
    grepl("gastrointestinal stromal|gist", d_lower) ~ "GIST",
    grepl("gastric|stomach|gastroesophageal", d_lower) ~ "STAD",
    grepl("esophageal|esca", d_lower) ~ "ESCA",
    grepl("urothelial|bladder|blca", d_lower) ~ "BLCA",
    grepl("renal cell|kidney|rcc", d_lower) ~ "RCC",
    grepl("glioblastoma|gbm", d_lower) ~ "GBM",
    grepl("glioma|astrocytoma|oligodendroglioma", d_lower) ~ "GLIOMA",
    grepl("endometrial|uterine corpus", d_lower) ~ "UCEC",
    grepl("cervical", d_lower) ~ "CESC",
    grepl("head and neck|hnsc", d_lower) ~ "HNSC",
    grepl("merkel cell", d_lower) ~ "MCC",
    grepl("salivary gland", d_lower) ~ "SACA",
    grepl("multiple myeloma|plasma cell", d_lower) ~ "MM",
    grepl("lymphoma|dlbcl|hodgkin", d_lower) ~ "LYMPH",
    TRUE ~ NA_character_
  )
}

stress_candidates <- civic_long |>
  filter(!gene %in% cohort_genes) |>
  mutate(tumor_type = disease_to_oncotree(disease)) |>
  filter(!is.na(tumor_type))

log_info("Stress candidates (outside cohort, OncoTree-mappable): {nrow(stress_candidates)}")

# Dedup by (gene, variant_key, tumor_type) keeping strongest level
stress_candidates <- stress_candidates |>
  mutate(level_rank = case_when(
    grepl("TIER_I_LEVEL_A",  amp_level, fixed = TRUE) ~ 1L,
    grepl("TIER_I_LEVEL_B",  amp_level, fixed = TRUE) ~ 2L,
    grepl("TIER_II_LEVEL_C", amp_level, fixed = TRUE) ~ 3L,
    grepl("TIER_II_LEVEL_D", amp_level, fixed = TRUE) ~ 4L,
    TRUE ~ 99L
  )) |>
  arrange(gene, variant_key, tumor_type, level_rank) |>
  distinct(gene, variant_key, tumor_type, .keep_all = TRUE)

log_info("Stress cohort after dedup: {nrow(stress_candidates)}")

# ── Query OncoKB for each stress variant (cached) ────────────────────────────

stress_with_oncokb <- pmap_dfr(
  list(stress_candidates$gene, stress_candidates$variant_key,
       stress_candidates$tumor_type),
  function(gene, variant, tumor_type) {
    tryCatch({
      result <- oncokb_annotate_mutation(
        hugo_symbol    = gene,
        protein_change = variant,
        tumor_type     = tumor_type
      )
      tibble(
        gene = gene, variant = variant, tumor_type = tumor_type,
        has_oncokb_evidence = !is.na(result$oncogenic) && result$oncogenic != "Unknown",
        oncokb_oncogenic    = result$oncogenic %||% NA_character_,
        oncokb_level        = result$highest_sensitive_level %||% NA_character_,
        oncokb_mutation_effect = result$mutation_effect %||% NA_character_
      )
    }, error = function(e) {
      log_debug("OncoKB failed for {gene} {variant} {tumor_type}: {e$message}")
      tibble(
        gene = gene, variant = variant, tumor_type = tumor_type,
        has_oncokb_evidence = FALSE,
        oncokb_oncogenic    = NA_character_,
        oncokb_level        = NA_character_,
        oncokb_mutation_effect = NA_character_
      )
    })
  }
)

stress_full <- stress_candidates |>
  rename(variant = variant_key, civic_amp_level = amp_level,
         civic_disease = disease) |>
  inner_join(stress_with_oncokb,
             by = c("gene", "variant", "tumor_type"))

log_info("Stress cohort with OncoKB queried: {nrow(stress_full)}")

# ── Compute baseline + write inputs ──────────────────────────────────────────

oncokb_to_amp <- function(level) {
  case_when(
    level %in% c("LEVEL_1", "LEVEL_2") ~ "Tier I Level A",
    level == "LEVEL_3A"                ~ "Tier I Level B",
    level == "LEVEL_3B"                ~ "Tier II Level C",
    level == "LEVEL_4"                 ~ "Tier II Level D",
    TRUE                               ~ NA_character_
  )
}

walk(seq_len(nrow(stress_full)), function(i) {
  row <- stress_full[i, ]

  case_id <- paste0("stress_",
                    str_replace_all(row$gene, "[^A-Za-z0-9]", ""), "_",
                    str_replace_all(row$variant, "[^A-Za-z0-9]", ""), "_",
                    row$tumor_type, "_",
                    row$assertion_id)

  oncokb_evidence <- list(
    level           = row$oncokb_level,
    oncogenic       = row$oncokb_oncogenic,
    mutation_effect = row$oncokb_mutation_effect,
    evidence_text   = glue::glue(
      "OncoKB level {row$oncokb_level %||% 'NA'}; ",
      "oncogenic={row$oncokb_oncogenic %||% 'NA'}; ",
      "mutation_effect={row$oncokb_mutation_effect %||% 'NA'}. ",
      "Tumor type queried: {row$tumor_type}."
    )
  )
  civic_evidence <- list(
    amp_level     = row$civic_amp_level,
    tumor_context = row$civic_disease,
    evidence_text = glue::glue(
      "CiVIC AMP-level assertion: {row$civic_amp_level} ",
      "in disease context '{row$civic_disease}' ",
      "(mapped to OncoTree code '{row$tumor_type}'). ",
      "Assertion ID: {row$assertion_id}."
    )
  )

  user_turn <- format_reconciler_input(
    gene       = row$gene,
    variant    = row$variant,
    tumor_type = row$tumor_type,
    oncokb     = oncokb_evidence,
    civic      = civic_evidence,
    vep_impact = NA_character_,
    clinvar    = NA_character_
  )

  baseline <- classify_amp_tier(
    oncokb_oncogenic = row$oncokb_oncogenic,
    oncokb_level     = row$oncokb_level,
    civic_amp_level  = row$civic_amp_level,
    vep_impact       = NA_character_,
    clinvar          = NA_character_
  )

  metadata <- list(
    gene             = row$gene,
    variant          = row$variant,
    tumor_type       = row$tumor_type,
    oncokb_level     = row$oncokb_level %||% NA_character_,
    civic_amp_level  = row$civic_amp_level,
    civic_disease    = row$civic_disease,
    civic_evidence_in_different_tumor = FALSE,  # by construction in stress
    oncokb_expected_amp_tier = oncokb_to_amp(row$oncokb_level),
    civic_mapped_amp_tier    = case_when(
      grepl("TIER_I_LEVEL_A",  row$civic_amp_level, fixed = TRUE) ~ "Tier I Level A",
      grepl("TIER_I_LEVEL_B",  row$civic_amp_level, fixed = TRUE) ~ "Tier I Level B",
      grepl("TIER_II_LEVEL_C", row$civic_amp_level, fixed = TRUE) ~ "Tier II Level C",
      grepl("TIER_II_LEVEL_D", row$civic_amp_level, fixed = TRUE) ~ "Tier II Level D",
      TRUE ~ NA_character_
    ),
    expected_stratum  = "stress",
    discordance_signal = "stress",
    baseline_tier     = baseline$amp_tier,
    baseline_level    = baseline$amp_level,
    baseline_evidence = baseline$amp_evidence
  )

  write_agent_input(
    agent_name = "evidence_reconciler",
    case_id    = case_id,
    user_turn  = user_turn,
    metadata   = metadata
  )
})

log_info("Wrote {nrow(stress_full)} stress-cohort input(s) under agents/inputs/evidence_reconciler/ (case_id prefix: stress_)")
