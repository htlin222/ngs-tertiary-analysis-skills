# R/foundation_sections.R — htmltools tag builders for each report section.
# Sourced by foundation_report.R which composes them into the final document.

suppressPackageStartupMessages({
  library(htmltools); library(dplyr); library(stringr); library(glue); library(tibble)
})

source(here::here("R/foundation_helpers.R"))
source(here::here("R/oncokb_helpers.R"))

# ── Cover ────────────────────────────────────────────────────────────────────

section_cover <- function(sample_id, tumor_label, parsed, config) {
  ad <- parsed$analysis_details %||% list()
  tags$section(class = "page cover",
    tags$h1("Comprehensive Genomic Profile"),
    tags$p(class = "subtitle", "TruSight Oncology 500 — Clinician Summary"),
    tags$div(class = "cover-meta",
      tags$table(
        tags$tr(tags$th("Sample ID"),       tags$td(sample_id)),
        tags$tr(tags$th("Tumor Type"),      tags$td(tumor_label)),
        tags$tr(tags$th("Specimen"),        tags$td(ad[["DNA Sample ID"]] %||% sample_id)),
        tags$tr(tags$th("Panel"),           tags$td("Illumina TruSight Oncology 500 (523 genes, ~1.94 Mb coding)")),
        tags$tr(tags$th("Pipeline"),        tags$td(paste0("DRAGEN TSO500 v",
                                                            ad[["Module Version"]] %||% "?",
                                                            " · annotation: OncoKB + CiVIC"))),
        tags$tr(tags$th("Sample run date"), tags$td(ad[["Output Date"]] %||% "—")),
        tags$tr(tags$th("Report date"),     tags$td(format(Sys.time(), "%Y-%m-%d"))),
        tags$tr(tags$th("Institution"),     tags$td(config$report$institution %||% "Research Institution"))
      )
    )
  )
}

# ── Biomarker snapshot ───────────────────────────────────────────────────────

section_biomarkers <- function(tmb, msi, hrd) {
  fmt <- function(x) if (length(x) == 0 || is.na(x)) "—" else format(x, nsmall = 1)
  card <- function(kind, label, value, status, score) {
    band <- if (kind == "tmb") {
      if (!is.na(score) && score >= 10) "high"
      else if (!is.na(score) && score >= 6) "med" else "low"
    } else if (kind == "msi") {
      if (grepl("MSI-High", status %||% "")) "high" else "low"
    } else {
      if (!is.na(score) && score >= 42) "high"
      else if (!is.na(score) && score >= 30) "med" else "low"
    }
    tags$div(class = paste("biomarker-card", band),
      tags$div(class = "bm-label", label),
      tags$div(class = "bm-value", value),
      tags$div(class = "bm-status", status %||% "Unknown"),
      tags$div(class = "bm-blurb", biomarker_blurb(kind, score, status %||% ""))
    )
  }
  tags$section(class = "biomarkers",
    tags$h2("Biomarker Snapshot"),
    tags$div(class = "bm-grid",
      card("tmb", "Tumor Mutational Burden",
           paste0(fmt(tmb$tmb_score), " mut/Mb"),
           tmb$tmb_class, tmb$tmb_score),
      card("msi", "Microsatellite Instability",
           paste0(fmt(msi$msi_score), "% unstable"),
           msi$msi_status, msi$msi_score),
      card("hrd", "Genomic Instability (HRD proxy)",
           paste0("GIS ", round(hrd$hrd_score %||% NA_real_, 0)),
           hrd$hrd_status, hrd$hrd_score)
    )
  )
}

# ── Therapy hero tables ──────────────────────────────────────────────────────

therapy_table <- function(rows, empty_msg, off_label = FALSE) {
  if (nrow(rows) == 0) return(tags$p(class = "empty-section", empty_msg))
  tags$table(class = "therapy-table",
    tags$thead(tags$tr(
      tags$th("Therapy / Drug"),
      tags$th("Genomic Driver"),
      tags$th("OncoKB Level"),
      tags$th("Clinical Note")
    )),
    tags$tbody(lapply(seq_len(nrow(rows)), function(i) {
      r <- rows[i, ]
      lvl_class <- paste0("lvl-", tolower(gsub("LEVEL_", "", r$level)))
      tags$tr(class = lvl_class,
        tags$td(class = "drug-cell", tags$strong(r$drugs)),
        tags$td(class = "driver-cell", paste0(r$gene, " ", r$alteration)),
        tags$td(class = "level-cell",
                tags$span(class = paste0("lvl-badge ", lvl_class),
                          oncokb_level_pretty(r$level))),
        tags$td(class = "note-cell",
                if (off_label) "Evidence from a different tumor context — discuss off-label or trial enrollment."
                else "Approved or guideline-listed in this tumor type.")
      )
    }))
  )
}

section_therapies_hero <- function(therapies) {
  on  <- filter(therapies, tier == "on-label")
  off <- filter(therapies, tier == "off-label")
  total <- nrow(on) + nrow(off)
  tags$section(class = "therapies-hero",
    tags$h2("Therapies — Potential Clinical Benefit"),
    tags$p(class = "section-lede",
           glue("{total} therapy match(es) identified across OncoKB evidence levels.")),
    tags$h3("Sensitivity in this tumor type"),
    therapy_table(on, "No on-label sensitivity matches in this tumor type."),
    tags$h3("Sensitivity in other tumor types (off-label / preclinical)"),
    therapy_table(off, "No off-label or preclinical sensitivity matches.", off_label = TRUE)
  )
}

section_resistance <- function(therapies) {
  res <- filter(therapies, tier == "resistance")
  if (nrow(res) == 0) return(NULL)
  tags$section(class = "therapies-resistance",
    tags$h2("Therapies — Predicted Resistance / Lack of Response"),
    therapy_table(res, "(none)")
  )
}

# ── Variant cards (AMP Tier I-III) ───────────────────────────────────────────

variant_card <- function(row, oncokb_summary, civic_assertions) {
  hgvsp_short_val <- if (!is.na(row$hgvsp)) hgvsp_to_short(as.character(row$hgvsp)) else NA_character_
  oncokb_match <- oncokb_summary |>
    filter(gene == row$gene & (alteration == row$hgvsp |
                               alteration == hgvsp_short_val)) |>
    slice(1)
  civic_match <- civic_assertions |>
    filter(gene == row$gene) |>
    slice_head(n = 2)

  vaf_str <- if (!is.na(row$vaf))
    sprintf("%.1f%% (%d/%d)", row$vaf*100, row$alt_depth %||% 0L, row$total_depth %||% 0L)
    else "—"

  tier_slug <- gsub(" ", "-", tolower(row$amp_tier %||% "tier-iv"))
  amp_badge <- tags$span(class = paste0("amp-badge ", tier_slug), row$amp_tier %||% "Tier IV")

  blurb <- if (nrow(oncokb_match) > 0 && nzchar(oncokb_match$tumor_type_summary)) {
    oncokb_match$tumor_type_summary
  } else if (nrow(oncokb_match) > 0 && nzchar(oncokb_match$mut_eff_desc)) {
    str_trunc(oncokb_match$mut_eff_desc, 320, side = "right")
  } else NULL

  civic_block <- if (nrow(civic_match) > 0) {
    tags$div(class = "civic-cite",
      tags$h4("CiVIC assertion"),
      tags$ul(lapply(seq_len(nrow(civic_match)), function(i) {
        a <- civic_match[i, ]
        tags$li(tags$strong(a$amp_level %||% ""), " — ",
                a$significance %||% "", "; ", a$disease %||% "",
                if (!is.null(a$therapies) && !is.na(a$therapies) && nzchar(a$therapies))
                  glue(" · {a$therapies}") else "")
      }))
    )
  } else NULL

  meta_parts <- list(
    tags$span(class = "vc-vaf", paste0("VAF ", vaf_str)),
    if (!is.na(row$amp_level)) tags$span(class = "vc-amplvl", row$amp_level) else NULL,
    if (nrow(oncokb_match) > 0 && !is.na(oncokb_match$oncogenic))
      tags$span(class = "vc-onc", oncokb_match$oncogenic) else NULL,
    if (nrow(oncokb_match) > 0 && !is.na(oncokb_match$mutation_effect))
      tags$span(class = "vc-eff", oncokb_match$mutation_effect) else NULL,
    if (nrow(oncokb_match) > 0 && isTRUE(oncokb_match$hotspot))
      tags$span(class = "vc-hot", "Hotspot") else NULL
  )

  tags$div(class = "variant-card",
    tags$div(class = "vc-head",
      tags$span(class = "vc-gene", row$gene),
      tags$span(class = "vc-alt", row$hgvsp %||% (row$consequence %||% "")),
      amp_badge
    ),
    tags$div(class = "vc-meta", meta_parts),
    if (!is.null(row$amp_evidence) && nzchar(row$amp_evidence))
      tags$p(class = "vc-evidence", tags$strong("AMP rationale: "), row$amp_evidence) else NULL,
    if (!is.null(blurb)) tags$p(class = "vc-blurb", blurb) else NULL,
    civic_block
  )
}

# ── CNV / fusion card builders ───────────────────────────────────────────────

# Match a parsed CNV row (gene, type, fold_change) back to the OncoKB CNA
# annotation that came from `query_cnas`, so the card can show OncoKB level
# and tumor-type summary alongside the call from TSO500.
cnv_card <- function(parsed_row, oncokb_match, lr_meta = NULL) {
  fold_change <- if (!is.na(parsed_row$log2_ratio)) round(2 ^ parsed_row$log2_ratio, 2) else NA_real_
  type <- parsed_row$type
  display_alt <- switch(type,
    AMPLIFICATION = "Amplification", DELETION = "Deletion",
    GAIN = "Copy gain", LOSS = "Copy loss", type)

  tier_class <- if (!is.null(oncokb_match) &&
                     length(oncokb_match$highest_sensitive_level) > 0 &&
                     !is.na(oncokb_match$highest_sensitive_level)) {
    lvl <- oncokb_match$highest_sensitive_level
    if (lvl %in% c("LEVEL_1", "LEVEL_2", "LEVEL_3A")) "tier-i"
    else if (lvl %in% c("LEVEL_3B", "LEVEL_4")) "tier-ii"
    else "tier-iii"
  } else "tier-iii"

  tier_label <- if (tier_class == "tier-i") "Tier I"
                else if (tier_class == "tier-ii") "Tier II" else "Significant CNA"
  badge <- tags$span(class = paste0("amp-badge ", tier_class), tier_label)

  meta_parts <- list(
    tags$span(class = "vc-vaf",
              if (!is.na(fold_change)) paste0("Fold Change ", fold_change) else "—"),
    if (!is.na(parsed_row$copy_number))
      tags$span(class = "vc-amplvl", paste0("CN ", parsed_row$copy_number)) else NULL,
    if (!is.null(oncokb_match) && length(oncokb_match$oncogenic) > 0 &&
        !is.na(oncokb_match$oncogenic))
      tags$span(class = "vc-onc", oncokb_match$oncogenic) else NULL,
    if (!is.null(oncokb_match) && length(oncokb_match$highest_sensitive_level) > 0 &&
        !is.na(oncokb_match$highest_sensitive_level))
      tags$span(class = "vc-eff",
                paste0("OncoKB ", oncokb_match$highest_sensitive_level)) else NULL
  )

  blurb <- if (!is.null(oncokb_match) && length(oncokb_match$raw$tumorTypeSummary) > 0 &&
               nzchar(oncokb_match$raw$tumorTypeSummary[[1]])) {
    oncokb_match$raw$tumorTypeSummary[[1]]
  } else if (!is.null(oncokb_match) && length(oncokb_match$raw$variantSummary) > 0 &&
             nzchar(oncokb_match$raw$variantSummary[[1]])) {
    oncokb_match$raw$variantSummary[[1]]
  } else NULL

  evidence_lines <- character()
  if (!is.null(lr_meta) && nrow(lr_meta) > 0) {
    lr <- lr_meta[1, ]
    evidence_lines <- c(evidence_lines,
      paste0("Source: TSO500 [Large Rearrangements] — exon ",
             lr[["Affected Exon(s)"]] %||% "?",
             " (", lr$Chromosome %||% "?", ":",
             format(as.integer(lr$Start %||% 0), big.mark = ","), "–",
             format(as.integer(lr$Stop %||% 0), big.mark = ","), ")"))
  } else if (!is.na(parsed_row$chromosome)) {
    evidence_lines <- c(evidence_lines,
      paste0("Source: TSO500 [Copy Number Variants] — ",
             parsed_row$chromosome,
             if (!is.na(parsed_row$start))
               paste0(":", format(parsed_row$start, big.mark = ","), "–",
                      format(parsed_row$end, big.mark = ",")) else ""))
  } else {
    evidence_lines <- c(evidence_lines,
      "Source: TSO500 [Copy Number Variants] (gene-level)")
  }

  tags$div(class = "variant-card",
    tags$div(class = "vc-head",
      tags$span(class = "vc-gene", parsed_row$gene),
      tags$span(class = "vc-alt", display_alt),
      badge
    ),
    tags$div(class = "vc-meta", meta_parts),
    tags$p(class = "vc-evidence",
           paste(evidence_lines, collapse = " · ")),
    if (!is.null(blurb)) tags$p(class = "vc-blurb", blurb) else NULL
  )
}

fusion_card <- function(fusion_row, oncokb_match) {
  badge <- tags$span(class = "amp-badge tier-ii", "Fusion")
  blurb <- if (!is.null(oncokb_match) && length(oncokb_match$raw$tumorTypeSummary) > 0 &&
               nzchar(oncokb_match$raw$tumorTypeSummary[[1]])) {
    oncokb_match$raw$tumorTypeSummary[[1]]
  } else NULL
  tags$div(class = "variant-card",
    tags$div(class = "vc-head",
      tags$span(class = "vc-gene", paste0(fusion_row$gene_a, "::", fusion_row$gene_b)),
      tags$span(class = "vc-alt", fusion_row$fusion_name %||% "Fusion"),
      badge
    ),
    tags$div(class = "vc-meta",
      if (!is.na(fusion_row$supporting_reads))
        tags$span(class = "vc-vaf", paste0("Supporting reads: ", fusion_row$supporting_reads)) else NULL,
      if (!is.null(oncokb_match) && length(oncokb_match$oncogenic) > 0 &&
          !is.na(oncokb_match$oncogenic))
        tags$span(class = "vc-onc", oncokb_match$oncogenic) else NULL,
      if (!is.null(oncokb_match) && length(oncokb_match$highest_sensitive_level) > 0 &&
          !is.na(oncokb_match$highest_sensitive_level))
        tags$span(class = "vc-eff",
                  paste0("OncoKB ", oncokb_match$highest_sensitive_level)) else NULL
    ),
    if (!is.null(blurb)) tags$p(class = "vc-blurb", blurb) else NULL
  )
}

# Pull the parsed [Large Rearrangements] table to enrich CNV cards with
# breakpoint coordinates / affected exons. Matched by gene name.
attach_lr_metadata <- function(parsed_cnv, parsed_lr) {
  if (!is.data.frame(parsed_lr) || nrow(parsed_lr) == 0) return(parsed_cnv)
  parsed_cnv
}

section_variant_cards <- function(amp, oncokb, civic, parsed = NULL) {
  oncokb_summary <- oncokb_mutation_summary(oncokb)
  civic_assertions <- civic$assertions %||% tibble()

  # ── 1. Short variants — AMP Tier I/II ────────────────────────────────────
  snv_cards <- list()
  if (is.data.frame(amp) && nrow(amp) > 0) {
    amp_sig <- amp |>
      filter(amp_tier %in% c("Tier I", "Tier II")) |>
      arrange(factor(amp_tier, levels = c("Tier I", "Tier II"))) |>
      distinct(gene, hgvsp, .keep_all = TRUE)
    snv_cards <- lapply(seq_len(nrow(amp_sig)), function(i) {
      variant_card(amp_sig[i, ], oncokb_summary, civic_assertions)
    })
  }

  # ── 2. CNVs flagged by OncoKB (Level 1/2/3/4) ────────────────────────────
  cnv_cards <- list()
  cnas <- oncokb$cnas %||% list()
  if (length(cnas) > 0 && !is.null(parsed) && is.data.frame(parsed$cnv)) {
    lr_tbl <- parsed$large_rearrangements %||% tibble()
    cnv_lookup <- parsed$cnv |>
      group_by(gene) |>
      slice(1) |>
      ungroup()
    for (m in cnas) {
      if (is.null(m) || isTRUE(is.na(m))) next
      lvl <- m$highest_sensitive_level
      onc <- m$oncogenic %||% "Unknown"
      # Surface anything OncoKB flags as oncogenic OR with a sensitivity level
      keep <- (length(lvl) > 0 && !is.na(lvl) && nzchar(lvl)) ||
              (onc %in% c("Oncogenic", "Likely Oncogenic"))
      if (!keep) next
      pr <- cnv_lookup |> filter(gene == m$gene) |> slice(1)
      if (nrow(pr) == 0) {
        # Synthesize a minimal row when the CNV table doesn't carry the gene
        # (shouldn't happen, but defensive)
        pr <- tibble(gene = m$gene,
                     type = if (m$alteration == "AMPLIFICATION") "AMPLIFICATION"
                            else "DELETION",
                     log2_ratio = NA_real_, copy_number = NA_real_,
                     chromosome = NA_character_, start = NA_integer_, end = NA_integer_)
      }
      lr_meta <- if (nrow(lr_tbl) > 0) lr_tbl |> filter(Gene == m$gene) else tibble()
      cnv_cards[[length(cnv_cards) + 1L]] <- cnv_card(pr, m, lr_meta)
    }
  }

  # ── 3. Fusions (rare in this batch, but supported) ───────────────────────
  fus_cards <- list()
  fusions <- oncokb$fusions %||% list()
  if (length(fusions) > 0 && !is.null(parsed) && is.data.frame(parsed$fusions) &&
      nrow(parsed$fusions) > 0) {
    for (m in fusions) {
      if (is.null(m) || isTRUE(is.na(m))) next
      lvl <- m$highest_sensitive_level
      onc <- m$oncogenic %||% "Unknown"
      keep <- (length(lvl) > 0 && !is.na(lvl) && nzchar(lvl)) ||
              (onc %in% c("Oncogenic", "Likely Oncogenic"))
      if (!keep) next
      fr <- parsed$fusions |>
        filter((gene_a == m$gene_a & gene_b == m$gene_b) |
               (gene_a == m$gene_b & gene_b == m$gene_a)) |>
        slice(1)
      if (nrow(fr) == 0) next
      fus_cards[[length(fus_cards) + 1L]] <- fusion_card(fr, m)
    }
  }

  total <- length(snv_cards) + length(cnv_cards) + length(fus_cards)
  if (total == 0) {
    return(tags$section(class = "variant-cards",
      tags$h2("Variants of Potential Clinical Significance"),
      tags$p(class = "empty-section",
             "No clinically significant short variants, copy-number alterations, ",
             "or fusions identified for this tumor type.")))
  }

  body <- list(
    tags$h2("Variants of Potential Clinical Significance"),
    tags$p(class = "section-lede",
           glue("{total} significant alteration(s): ",
                "{length(snv_cards)} short variant(s), ",
                "{length(cnv_cards)} CNV(s), ",
                "{length(fus_cards)} fusion(s).")),
    if (length(snv_cards) > 0) tags$h3("Short variants") else NULL,
    if (length(snv_cards) > 0) do.call(tagList, snv_cards) else NULL,
    if (length(cnv_cards) > 0) tags$h3("Copy-number alterations") else NULL,
    if (length(cnv_cards) > 0) do.call(tagList, cnv_cards) else NULL,
    if (length(fus_cards) > 0) tags$h3("Gene fusions") else NULL,
    if (length(fus_cards) > 0) do.call(tagList, fus_cards) else NULL
  )
  body <- body[!vapply(body, is.null, logical(1))]
  do.call(tags$section, c(list(class = "variant-cards"), body))
}

# ── Literature Review (OpenEvidence per-sample summary) ─────────────────────

# Build an AMA-style entry from one citation record (from
# scripts/extract_oe_citations.sh). Returns an htmltools tag list.
ama_citation_li <- function(cit, n) {
  trim_dot <- function(s) sub("\\.+\\s*$", "", trimws(s))

  authors <- trim_dot(cit$authors %||% "Anonymous")
  title <- trim_dot(cit$title %||% "(untitled)")

  # OpenEvidence's `publication_info_string` already encodes the citation tail
  # as "Journal. Year;Vol(Issue):Pages. doi:DOI." Strip the trailing doi
  # fragment so we can render the DOI as its own clickable badge.
  pub <- cit$citation_string %||% ""
  pub <- sub("\\s*doi:.*$", "", pub, ignore.case = TRUE)
  pub <- trim_dot(pub)
  if (!nzchar(pub)) {
    yr <- if (!is.null(cit$year)) substr(cit$year, 1, 4) else ""
    pub <- paste0(cit$journal_short %||% cit$journal %||% "", ". ", yr)
  }

  link_label <- if (!is.null(cit$doi) && nzchar(as.character(cit$doi))) {
    paste0("doi:", cit$doi)
  } else if (!is.null(cit$pmid) && !is.na(cit$pmid)) {
    paste0("PMID:", cit$pmid)
  } else NULL

  link_href <- if (!is.null(cit$doi) && nzchar(as.character(cit$doi))) {
    paste0("https://doi.org/", cit$doi)
  } else if (!is.null(cit$pmid) && !is.na(cit$pmid)) {
    paste0("https://pubmed.ncbi.nlm.nih.gov/", cit$pmid, "/")
  } else cit$href %||% NULL

  link_tag <- if (!is.null(link_href) && !is.null(link_label))
    tagList(" ", tags$a(href = link_href, target = "_blank",
                        rel = "noopener noreferrer", link_label))
  else NULL

  tags$li(value = n,
    tags$span(class = "cit-authors", paste0(authors, ".")),
    " ",
    tags$span(class = "cit-title", paste0(title, ".")),
    " ",
    tags$span(class = "cit-pub", paste0(pub, ".")),
    link_tag
  )
}

# Reads two artefacts produced by the OpenEvidence pipeline:
#   reports/.openevidence_cache/sample_<sid>.md             — answer prose
#   reports/.openevidence_cache/sample_<sid>.citations.json — ranked refs
# Renders both: prose + an AMA-style references list with DOI/PubMed links so
# clinicians can validate every claim.
section_lit_review <- function(sample_id) {
  md_path <- here::here("reports", ".openevidence_cache",
                        paste0("sample_", sample_id, ".md"))
  cit_path <- here::here("reports", ".openevidence_cache",
                         paste0("sample_", sample_id, ".citations.json"))
  if (!file.exists(md_path)) return(NULL)

  raw <- paste(readLines(md_path, warn = FALSE), collapse = "\n")
  if (!nzchar(trimws(raw))) return(NULL)

  # Strip the OpenEvidence numeric-index inline citations ([12][24]) — they
  # are an internal ranking; the references list below carries the actual
  # bibliographic data with DOI links.
  cleaned <- gsub("\\[[0-9]+\\](\\[[0-9]+\\])*", "", raw)
  cleaned <- sub("\\n\\nWould you like to explore[^\\n]*\\??$", "", cleaned)
  cleaned <- gsub("<visual>[^<]+</visual>", "", cleaned)

  body <- if (requireNamespace("commonmark", quietly = TRUE)) {
    HTML(commonmark::markdown_html(cleaned))
  } else {
    tags$pre(cleaned)
  }

  # ── References ────────────────────────────────────────────────────────────
  refs_block <- NULL
  if (file.exists(cit_path)) {
    citations <- tryCatch(
      jsonlite::fromJSON(cit_path, simplifyVector = FALSE),
      error = function(e) list()
    )
    if (length(citations) > 0) {
      # Deduplicate by DOI then by title; prefer entries with a DOI.
      keys <- vapply(citations, function(c) {
        c$doi %||% c$title %||% c$href %||% ""
      }, character(1))
      keep <- !duplicated(keys) & nzchar(keys)
      citations <- citations[keep]
      # Sort: entries with DOI first, NCCN guidelines + FDA last; preserve
      # otherwise-stable order.
      ord <- order(
        vapply(citations, function(c)
          if (!is.null(c$doi) && nzchar(c$doi)) 0L else 1L, integer(1)),
        seq_along(citations)
      )
      citations <- citations[ord]
      lis <- lapply(seq_along(citations), function(i) {
        ama_citation_li(citations[[i]], i)
      })
      refs_block <- tags$div(class = "lit-refs",
        tags$h3(glue("References ({length(citations)})")),
        do.call(tags$ol, c(class = "ama-refs", lis))
      )
    }
  }

  tags$section(class = "lit-review",
    tags$h2("Literature Review (OpenEvidence)"),
    tags$p(class = "section-lede",
           "Variant-by-variant clinical evidence synthesis from peer-reviewed ",
           "literature, FDA labels, and NCCN/ESMO guidelines via OpenEvidence. ",
           "Every claim is traceable to an AMA-style reference with DOI link below."),
    tags$div(class = "lit-body", body),
    refs_block,
    tags$p(class = "lit-disclaimer",
           tags$strong("Source: "),
           "OpenEvidence (",
           tags$a(href = "https://www.openevidence.com", "openevidence.com"),
           "). Inline OpenEvidence index markers were stripped; the AMA ",
           "references above carry full bibliographic detail with clickable ",
           "DOI / PubMed identifiers for validation.")
  )
}

# ── CNV table ────────────────────────────────────────────────────────────────

section_cnv <- function(cnv) {
  if (!is.data.frame(cnv) || nrow(cnv) == 0) return(NULL)
  tbl <- cnv |>
    filter(type %in% c("AMPLIFICATION", "DELETION", "GAIN", "LOSS")) |>
    arrange(factor(type, levels = c("AMPLIFICATION", "GAIN", "DELETION", "LOSS")), gene) |>
    transmute(
      Gene = gene,
      Type = type,
      `Fold Change` = round(2 ^ log2_ratio, 2),
      `Copy Number` = ifelse(is.na(copy_number), "—", as.character(copy_number))
    )
  if (nrow(tbl) == 0) return(NULL)
  tags$section(class = "cnv-section",
    tags$h2("Copy Number Alterations"),
    tags$table(class = "cnv-table",
      tags$thead(tags$tr(lapply(names(tbl), tags$th))),
      tags$tbody(lapply(seq_len(nrow(tbl)), function(i) {
        tags$tr(class = paste0("cnv-", tolower(tbl$Type[i])),
                lapply(names(tbl), function(c) tags$td(tbl[[c]][i])))
      }))
    )
  )
}

# ── Embedded PNG figures ─────────────────────────────────────────────────────

section_figures <- function(vaf_png, cnv_png) {
  embed <- function(p) {
    if (is.null(p) || !file.exists(p)) return(NULL)
    raw <- readBin(p, "raw", file.info(p)$size)
    paste0("data:image/png;base64,", jsonlite::base64_enc(raw))
  }
  vu <- embed(vaf_png); cu <- embed(cnv_png)
  if (is.null(vu) && is.null(cu)) return(NULL)
  tags$section(class = "figures",
    tags$h2("Variant + CNV Distribution"),
    if (!is.null(vu)) tags$figure(
      tags$img(src = vu, alt = "Top variants by VAF"),
      tags$figcaption("Top 30 somatic variants by allele frequency, coloured by AMP/pathogenicity class.")
    ) else NULL,
    if (!is.null(cu)) tags$figure(
      tags$img(src = cu, alt = "Genome-wide CNV log2"),
      tags$figcaption("Gene-level copy-number log2 fold-change across the TSO500 panel.")
    ) else NULL
  )
}

# ── Per-gene narrative + AMA refs ────────────────────────────────────────────

section_gene_narratives <- function(literature, oncokb, max_genes = 8) {
  narr <- literature$narratives %||% list()
  if (length(narr) == 0) return(NULL)
  actionable_genes <- vapply(oncokb$mutations %||% list(), function(m) {
    if (is.null(m) || isTRUE(is.na(m))) return(NA_character_)
    if (length(m$highest_sensitive_level) > 0 && !is.na(m$highest_sensitive_level)) m$gene
    else NA_character_
  }, character(1))
  actionable_genes <- unique(stats::na.omit(actionable_genes))

  keys <- names(narr)
  prio <- if (length(actionable_genes) > 0)
    keys[grepl(paste(actionable_genes, collapse = "|"), keys)] else character()
  others <- setdiff(keys, prio)
  use <- head(c(prio, others), max_genes)

  blocks <- lapply(use, function(k) {
    gene <- str_extract(k, "^[^_]+")
    txt <- narr[[k]]
    body <- if (requireNamespace("commonmark", quietly = TRUE)) {
      HTML(commonmark::markdown_html(txt))
    } else tags$pre(txt)
    tags$div(class = "gene-narrative",
      tags$h3(gene),
      body
    )
  })

  refs <- if (is.data.frame(literature$pubmed) && nrow(literature$pubmed) > 0) {
    pm <- literature$pubmed |> distinct(pmid, .keep_all = TRUE) |> head(15)
    tags$div(class = "pm-refs",
      tags$h3("References"),
      tags$ol(lapply(seq_len(nrow(pm)), function(i) {
        r <- pm[i, ]
        first_author <- str_extract(r$authors %||% "", "^[^,]+") %||% "Unknown"
        n <- str_count(r$authors %||% "", ",") + 1
        astr <- if (!is.na(n) && n > 3) paste0(first_author, ", et al")
                else (r$authors %||% first_author)
        tags$li(astr, ". ", r$title %||% "", ". ",
                tags$em(r$journal %||% ""), ". ", r$year %||% "", ". ",
                "PMID: ",
                tags$a(href = paste0("https://pubmed.ncbi.nlm.nih.gov/", r$pmid, "/"),
                       r$pmid))
      }))
    )
  } else NULL

  tags$section(class = "gene-narratives",
    tags$h2("Per-Gene Therapeutic Narrative"),
    do.call(tagList, blocks),
    refs
  )
}

# ── VUS appendix ─────────────────────────────────────────────────────────────

section_vus <- function(amp, max_rows = 200) {
  if (!is.data.frame(amp) || nrow(amp) == 0) return(NULL)
  # AMP Tier III = VUS, Tier IV = benign.
  vus <- amp |> filter(amp_tier %in% c("Tier III", "Tier IV"))
  if (nrow(vus) == 0) return(NULL)
  total <- nrow(vus)
  # For readability, only show coding non-synonymous variants in the table
  # body; collapse synonymous + UTR + intronic into a tail count.
  coding <- vus |>
    filter(!is.na(hgvsp), !grepl("synonymous|UTR|intron", consequence %||% ""))
  truncated <- nrow(coding) > max_rows
  shown <- coding |> arrange(gene) |> slice_head(n = max_rows)
  tbl <- shown |>
    transmute(
      Gene = gene,
      Alteration = hgvsp,
      VAF = ifelse(!is.na(vaf), sprintf("%.1f%%", vaf*100), "—"),
      Consequence = consequence %||% "",
      `AMP Tier` = amp_tier
    )
  noted <- nrow(coding)
  hidden <- total - noted
  caption <- glue("AMP Tier III/IV — no current evidence of clinical ",
                  "significance. Showing {nrow(tbl)} of {noted} coding variants",
                  if (truncated) glue(" (table truncated)") else "",
                  if (hidden > 0) glue("; {hidden} synonymous / non-coding ",
                                       "variants omitted from this view") else "",
                  ".")
  tags$section(class = "vus",
    tags$h2(glue("Variants of Unknown Significance ({total})")),
    tags$p(class = "section-lede", caption),
    tags$table(class = "vus-table",
      tags$thead(tags$tr(lapply(names(tbl), tags$th))),
      tags$tbody(lapply(seq_len(nrow(tbl)), function(i) {
        tags$tr(lapply(names(tbl), function(c) tags$td(tbl[[c]][i])))
      }))
    )
  )
}
