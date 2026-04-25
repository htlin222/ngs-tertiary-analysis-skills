#!/usr/bin/env Rscript
# benchmarks/agent/01_prepare_real_slice.R
#
# Phase 1 of Benchmark A real-data slice for the evidence_reconciler agent.
#
# Joins the existing 100-variant OncoKB cohort (benchmarks/results/evidence_concordance.csv)
# with fresh CiVIC AMP-level assertions, identifies cases of genuine discordance
# (OncoKB-mapped AMP tier ≠ CiVIC AMP-level assertion, or OncoKB-Neutral vs CiVIC
# Tier I/II), and writes one input JSON per discordant variant under
# agents/inputs/evidence_reconciler/.
#
# Phase 2: operator invokes the evidence_reconciler subagent on each input.
# Phase 3: 02_analyze_real_slice.R reads outputs, compares to baseline, writes summary.

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

# ── Inputs ────────────────────────────────────────────────────────────────────

cohort_path <- here::here("benchmarks/results/evidence_concordance.csv")
if (!file.exists(cohort_path)) {
  stop(glue("Missing {cohort_path}. Run benchmarks/scripts/05_evidence_concordance.R first."))
}
cohort <- read.csv(cohort_path, stringsAsFactors = FALSE) |> as_tibble()
log_info("Loaded cohort: {nrow(cohort)} variants from existing OncoKB benchmark")

# Restrict to rows where OncoKB has evidence (these are the variants worth reconciling)
have_oncokb <- cohort |> filter(has_oncokb_evidence)
log_info("Variants with OncoKB evidence: {nrow(have_oncokb)}")

# ── Pull all CiVIC assertions in one paginated walk ──────────────────────────
# civic_client.R::civic_get_assertions() fetches the first 100 assertions per
# call, so calling it once-per-gene burns 36 round-trips on the same first
# page. We do a single paginated GraphQL walk and filter in R instead.

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
            disease { name }
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

cohort_genes <- unique(have_oncokb$gene)
log_info("Cohort genes: {length(cohort_genes)}")

raw_nodes <- fetch_all_civic_assertions()

# Filter: ACCEPTED + has ampLevel + variant whose feature.name is a single gene
# (drop fusions like "v::ALK Fusion") + gene matches a cohort gene
civic_assertions_long <- map_dfr(raw_nodes, function(n) {
  if (toupper(n$status %||% "") != "ACCEPTED") return(NULL)
  if (is.null(n$ampLevel)) return(NULL)
  vs <- n$molecularProfile$variants %||% list()
  feats <- unique(unlist(map(vs, ~ .x$feature$name %||% NA_character_)))
  feats <- feats[!is.na(feats) & !grepl("::", feats, fixed = TRUE)]
  if (length(feats) == 0) return(NULL)
  mp_name <- n$molecularProfile$name %||% NA_character_
  map_dfr(feats, function(g) {
    tibble(
      gene           = g,
      mp_name        = mp_name,
      amp_level      = n$ampLevel,
      disease        = n$disease$name %||% NA_character_,
      assertion_id   = as.integer(n$id %||% NA_integer_)
    )
  })
})

log_info("ACCEPTED + has-ampLevel + single-gene assertions: {nrow(civic_assertions_long)}")

# NB: do NOT dedup by (gene, variant_key) here — the SAME variant can have
# different AMP-level assertions in different tumor types (e.g. BRAF V600E
# Tier I A in CRC AND in MEL). Dedup later, after picking the disease-matched
# entry per cohort row.
civic_long <- civic_assertions_long |>
  filter(gene %in% cohort_genes) |>
  mutate(variant_key = str_remove(mp_name, paste0("^", gene, "\\s+"))) |>
  mutate(variant_key = str_remove(variant_key, "^p\\.")) |>
  mutate(level_rank = case_when(
    grepl("TIER_I_LEVEL_A",  amp_level, fixed = TRUE) ~ 1L,
    grepl("TIER_I_LEVEL_B",  amp_level, fixed = TRUE) ~ 2L,
    grepl("TIER_II_LEVEL_C", amp_level, fixed = TRUE) ~ 3L,
    grepl("TIER_II_LEVEL_D", amp_level, fixed = TRUE) ~ 4L,
    TRUE ~ 99L
  )) |>
  select(gene, variant_key, amp_level, disease, level_rank)

log_info("Cohort-matching CiVIC assertions (long form, multi-disease): {nrow(civic_long)}")

# Map cohort OncoTree codes → CiVIC disease keywords for fuzzy matching.
tumor_keywords <- list(
  NSCLC = c("Lung", "Non-small", "NSCLC"),
  COAD  = c("Colorectal", "Colon", "COAD"),
  MEL   = c("Melanoma"),
  BRCA  = c("Breast"),
  OVT   = c("Ovarian")
)

# For each (gene, variant_key, tumor_type), pick the disease-matched CiVIC
# entry if any; otherwise pick the strongest cross-tumor entry and flag it.
civic_match_for <- function(gene_, vk_, tumor_) {
  cands <- civic_long |> filter(gene == gene_, variant_key == vk_)
  if (nrow(cands) == 0) {
    return(tibble(civic_amp_level = NA_character_,
                  civic_disease   = NA_character_,
                  civic_evidence_in_different_tumor = NA))
  }
  kws <- tumor_keywords[[tumor_]] %||% character()
  # Disease-matched candidates first
  if (length(kws) > 0) {
    pat <- paste0("(?i)", paste(kws, collapse = "|"))
    matched <- cands |> filter(grepl(pat, disease, perl = TRUE))
    if (nrow(matched) > 0) {
      best <- matched |> arrange(level_rank) |> slice(1)
      return(tibble(civic_amp_level = best$amp_level,
                    civic_disease   = best$disease,
                    civic_evidence_in_different_tumor = FALSE))
    }
  }
  # No matched disease — fall back to strongest cross-tumor entry (flagged)
  best <- cands |> arrange(level_rank) |> slice(1)
  tibble(civic_amp_level = best$amp_level,
         civic_disease   = best$disease,
         civic_evidence_in_different_tumor = TRUE)
}

log_info("CiVIC long-form rows ready for tumor-context join: {nrow(civic_long)}")

# ── Join cohort → CiVIC AMP-level (tumor-context-aware) ─────────────────────

joined <- have_oncokb |>
  mutate(variant_key = str_remove(variant, "^p\\.")) |>
  rowwise() |>
  mutate(civic = list(civic_match_for(gene, variant_key, tumor_type))) |>
  ungroup() |>
  unnest(civic)

# Diagnostics
n_with_civic_amp     <- sum(!is.na(joined$civic_amp_level))
n_with_civic_matched <- sum(!is.na(joined$civic_amp_level) &
                            !joined$civic_evidence_in_different_tumor)
n_with_civic_xtumor  <- sum(!is.na(joined$civic_amp_level) &
                            joined$civic_evidence_in_different_tumor)
log_info("Cohort variants with CiVIC AMP-level: {n_with_civic_amp}")
log_info("  - tumor-matched : {n_with_civic_matched}")
log_info("  - different tumor (cross-tumor evidence): {n_with_civic_xtumor}")

# ── Identify discordances ────────────────────────────────────────────────────
#
# Map OncoKB level → expected AMP tier (per amp_classification.R rules)
oncokb_to_amp <- function(level) {
  case_when(
    level %in% c("LEVEL_1", "LEVEL_2") ~ "Tier I Level A",
    level == "LEVEL_3A"                ~ "Tier I Level B",
    level == "LEVEL_3B"                ~ "Tier II Level C",
    level == "LEVEL_4"                 ~ "Tier II Level D",
    TRUE                               ~ NA_character_
  )
}

# Map CiVIC amp_level (TIER_I_LEVEL_A etc.) → AMP tier label
civic_to_amp <- function(amp_level) {
  case_when(
    grepl("TIER_I_LEVEL_A",  amp_level, fixed = TRUE) ~ "Tier I Level A",
    grepl("TIER_I_LEVEL_B",  amp_level, fixed = TRUE) ~ "Tier I Level B",
    grepl("TIER_II_LEVEL_C", amp_level, fixed = TRUE) ~ "Tier II Level C",
    grepl("TIER_II_LEVEL_D", amp_level, fixed = TRUE) ~ "Tier II Level D",
    TRUE                                              ~ NA_character_
  )
}

flagged <- joined |>
  mutate(
    oncokb_expected_amp = oncokb_to_amp(oncokb_level),
    civic_amp           = civic_to_amp(civic_amp_level),
    # NB: isTRUE() is scalar-only — use ` & !is.na(...)` for vector logic
    discordant = case_when(
      # CiVIC evidence is in a DIFFERENT tumor type than queried — real tumor-
      # specificity discordance regardless of nominal tier match
      !is.na(civic_amp_level) &
        !is.na(civic_evidence_in_different_tumor) &
        civic_evidence_in_different_tumor ~ TRUE,
      # Same-tumor: tiers differ
      !is.na(oncokb_expected_amp) & !is.na(civic_amp) & oncokb_expected_amp != civic_amp ~ TRUE,
      # OncoKB Neutral but CiVIC has Tier I/II assertion (oncogenicity mismatch)
      oncokb_oncogenic %in% c("Neutral", "Inconclusive") & !is.na(civic_amp) ~ TRUE,
      TRUE ~ FALSE
    )
  )

n_discordant <- sum(flagged$discordant)
log_info("Of {nrow(flagged)} cohort variants: {n_discordant} discordant, {nrow(flagged) - n_discordant} concordant")

# Write inputs for ALL OncoKB-positive variants (not just discordant). The
# concordant majority lets us measure the agent's false-positive rate (does it
# wrongly overrule a correct baseline?), which is the safety metric reviewers
# expect alongside discordance recall.
discordant <- flagged |>
  mutate(case_id = paste0("real_",
                          str_replace_all(gene, "[^A-Za-z0-9]", ""), "_",
                          str_replace_all(variant, "[^A-Za-z0-9]", ""), "_",
                          tumor_type))

cat("\n─── Discordant variants ───────────────────────────────────────\n")
print(discordant |>
        select(case_id, gene, variant, tumor_type,
               oncokb_level, oncokb_expected_amp,
               civic_amp_level, civic_amp,
               oncokb_oncogenic),
      n = 50)

# ── Write input JSONs ────────────────────────────────────────────────────────

if (n_discordant == 0) {
  log_warn("No discordant variants found — nothing to do.")
  quit(status = 0)
}

walk(seq_len(nrow(discordant)), function(i) {
  row <- discordant[i, ]

  # Synthesize the evidence dossier from cohort data + CiVIC assertion
  oncokb_evidence <- list(
    level           = row$oncokb_level,
    oncogenic       = row$oncokb_oncogenic,
    mutation_effect = row$oncokb_mutation_effect,
    evidence_text   = glue::glue(
      "OncoKB level {row$oncokb_level %||% 'NA'}; oncogenic={row$oncokb_oncogenic %||% 'NA'}; ",
      "mutation_effect={row$oncokb_mutation_effect %||% 'NA'}. ",
      "Tumor type queried: {row$tumor_type}."
    )
  )
  # Make tumor-context mismatch explicit in the evidence text so the agent does
  # not have to guess whether 'Colorectal Cancer' matches the queried 'NSCLC'.
  civic_text <- if (is.na(row$civic_amp_level)) {
    glue::glue("No tumor-matched CiVIC AMP-level assertion found for this variant ",
               "in the queried tumor context ({row$tumor_type}). ",
               "Total supporting evidence items in CiVIC (any disease): ",
               "{row$civic_evidence_count}.")
  } else if (isTRUE(row$civic_evidence_in_different_tumor)) {
    glue::glue("CiVIC AMP-level assertion: {row$civic_amp_level} ",
               "but in a DIFFERENT tumor context: '{row$civic_disease}'. ",
               "(Queried tumor: {row$tumor_type}.) ",
               "No same-tumor AMP-level assertion exists in CiVIC. ",
               "Total supporting evidence items: {row$civic_evidence_count}.")
  } else {
    glue::glue("CiVIC AMP-level assertion: {row$civic_amp_level} ",
               "in tumor-matched disease context '{row$civic_disease}'. ",
               "Total supporting evidence items: {row$civic_evidence_count}.")
  }
  civic_evidence <- list(
    amp_level       = row$civic_amp_level,
    evidence_count  = row$civic_evidence_count,
    tumor_context   = row$civic_disease,
    evidence_text   = civic_text
  )

  user_turn <- format_reconciler_input(
    gene       = row$gene,
    variant    = row$variant,
    tumor_type = row$tumor_type,
    oncokb     = oncokb_evidence,
    civic      = civic_evidence,
    vep_impact = NA_character_,         # not in this benchmark cohort
    clinvar    = NA_character_
  )

  # Deterministic baseline from the same inputs
  baseline <- classify_amp_tier(
    oncokb_oncogenic = row$oncokb_oncogenic,
    oncokb_level     = row$oncokb_level,
    civic_amp_level  = row$civic_amp_level,
    vep_impact       = NA_character_,
    clinvar          = NA_character_
  )

  expected_discordant <- isTRUE(row$discordant)
  metadata <- list(
    gene              = row$gene,
    variant           = row$variant,
    tumor_type        = row$tumor_type,
    oncokb_level      = row$oncokb_level %||% NA_character_,
    civic_amp_level   = row$civic_amp_level %||% NA_character_,
    civic_disease     = row$civic_disease %||% NA_character_,
    civic_evidence_in_different_tumor = row$civic_evidence_in_different_tumor %||% NA,
    oncokb_expected_amp_tier = row$oncokb_expected_amp %||% NA_character_,
    civic_mapped_amp_tier    = row$civic_amp %||% NA_character_,
    expected_stratum  = if (expected_discordant) "discordant" else "concordant",
    discordance_signal = if (!expected_discordant) {
      "concordant"
    } else if (row$oncokb_oncogenic %in% c("Neutral", "Inconclusive")) {
      "oncogenicity"
    } else {
      "tier_level"
    },
    baseline_tier     = baseline$amp_tier,
    baseline_level    = baseline$amp_level,
    baseline_evidence = baseline$amp_evidence
  )

  write_agent_input(
    agent_name = "evidence_reconciler",
    case_id    = row$case_id,
    user_turn  = user_turn,
    metadata   = metadata
  )
})

log_info("Wrote {nrow(discordant)} real-slice input(s) under agents/inputs/evidence_reconciler/")
log_info("Next: operator invokes the evidence_reconciler subagent once per case.")
