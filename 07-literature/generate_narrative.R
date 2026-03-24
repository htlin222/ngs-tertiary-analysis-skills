# 07-literature/generate_narrative.R — Orchestrate literature search and generate clinical narratives

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(glue)
  library(logger)
  library(jsonlite)
  library(tidyr)
})

source("R/api_clients.R")
source("07-literature/search_pubmed.R")
source("07-literature/search_scopus.R")

#' Generate clinical narratives for key alterations
#'
#' @param variants Tibble with gene, hgvsp, consequence, classification columns
#' @param cnv Tibble with copy number variants (gene, alteration, tier columns)
#' @param fusions Tibble with structural variants (gene_a, gene_b, tier columns)
#' @param config List with sample config and literature settings
#' @param sample_id Sample identifier
#'
#' @return List with: pubmed (tibble), scopus (tibble), narratives (named list of paragraphs)
#' @export
#'
#' @details
#' Orchestrates the full literature stage:
#'  1. Searches PubMed and Scopus for each significant variant
#'  2. For ESCAT tier I-III alterations or high-impact variants:
#'     - Generates narrative paragraph with clinical interpretation
#'     - Includes gene/alteration name, oncogenic classification, variant frequency,
#'       therapeutic relevance (from OncoKB), and key paper citations
#'  3. Saves combined results to reports/{sample_id}/07-literature/narratives.json
#'
#' Narrative template: "{gene} {alteration} is a {oncogenic_classification} alteration.
#' {frequency_statement}. {therapeutic_relevance}. {key_references}."
#'
#' @examples
#' \dontrun{
#'   variants <- tibble(
#'     gene = "BRAF", hgvsp = "p.V600E", consequence = "missense_variant",
#'     classification = "pathogenic"
#'   )
#'   cnv <- tibble(gene = "MYC", alteration = "AMPLIFICATION", tier = "Tier I")
#'   fusions <- tibble(gene_a = "EWSR1", gene_b = "FLI1", tier = "Tier I")
#'   config <- list(
#'     sample = list(tumor_type = "MEL"),
#'     literature = list(pubmed = list(max_results = 20), scopus = list(max_results = 20))
#'   )
#'   results <- generate_narrative(variants, cnv, fusions, config, "SAMPLE123")
#' }
generate_narrative <- function(variants, cnv, fusions, config, sample_id) {
  log_info("Starting literature narrative generation for sample: {sample_id}")

  # Validate inputs
  if (is.null(variants)) {
    variants <- tibble(
      gene = character(), hgvsp = character(), consequence = character(),
      classification = character()
    )
  }
  if (is.null(cnv)) {
    cnv <- tibble(gene = character(), alteration = character(), tier = character())
  }
  if (is.null(fusions)) {
    fusions <- tibble(gene_a = character(), gene_b = character(), tier = character())
  }

  if (is.null(config$sample$tumor_type)) {
    stop("config$sample$tumor_type is required")
  }

  tumor_type <- config$sample$tumor_type

  # Phase 1: Search PubMed
  log_info("Phase 1: Searching PubMed")
  pubmed_results <- tryCatch(
    {
      if (nrow(variants) > 0) {
        search_pubmed_for_variants(variants, config, sample_id)
      } else {
        tibble(
          pmid = character(), title = character(), authors = character(),
          journal = character(), year = character(), abstract = character(),
          query_gene = character(), query_variant = character()
        )
      }
    },
    error = function(e) {
      log_error("PubMed search failed: {conditionMessage(e)}")
      tibble(
        pmid = character(), title = character(), authors = character(),
        journal = character(), year = character(), abstract = character(),
        query_gene = character(), query_variant = character()
      )
    }
  )

  log_info("PubMed search returned {nrow(pubmed_results)} articles")

  # Phase 2: Search Scopus
  log_info("Phase 2: Searching Scopus")
  scopus_results <- tryCatch(
    {
      if (nrow(variants) > 0) {
        search_scopus_for_variants(variants, config, sample_id)
      } else {
        tibble(
          scopus_id = character(), title = character(), authors = character(),
          journal = character(), year = character(), doi = character(),
          cited_by = integer(), is_oa = logical(), oa_url = character(),
          oa_status = character(), query_gene = character(), query_variant = character()
        )
      }
    },
    error = function(e) {
      log_error("Scopus search failed: {conditionMessage(e)}")
      tibble(
        scopus_id = character(), title = character(), authors = character(),
        journal = character(), year = character(), doi = character(),
        cited_by = integer(), is_oa = logical(), oa_url = character(),
        oa_status = character(), query_gene = character(), query_variant = character()
      )
    }
  )

  log_info("Scopus search returned {nrow(scopus_results)} articles")

  # Phase 3: Identify key alterations for narrative generation
  log_info("Phase 3: Identifying key alterations for narratives")

  key_alterations <- list()

  # Add SNVs/Indels (pathogenic or likely pathogenic)
  if (nrow(variants) > 0) {
    snv_alts <- variants %>%
      filter(classification %in% c("pathogenic", "likely_pathogenic")) %>%
      select(gene, hgvsp) %>%
      distinct()

    if (nrow(snv_alts) > 0) {
      key_alterations[["snv"]] <- snv_alts %>%
        pmap(function(gene, hgvsp) {
          list(type = "SNV/Indel", gene = gene, alteration = hgvsp)
        })
    }
  }

  # Add CNVs (Tier I-III)
  if (nrow(cnv) > 0) {
    cnv_alts <- cnv %>%
      filter(tier %in% c("Tier I", "Tier II", "Tier III")) %>%
      select(gene, alteration) %>%
      distinct()

    if (nrow(cnv_alts) > 0) {
      key_alterations[["cnv"]] <- cnv_alts %>%
        pmap(function(gene, alteration) {
          list(type = "CNV", gene = gene, alteration = alteration)
        })
    }
  }

  # Add Fusions (Tier I-III)
  if (nrow(fusions) > 0) {
    fusion_alts <- fusions %>%
      filter(tier %in% c("Tier I", "Tier II", "Tier III")) %>%
      select(gene_a, gene_b) %>%
      distinct()

    if (nrow(fusion_alts) > 0) {
      key_alterations[["fusion"]] <- fusion_alts %>%
        pmap(function(gene_a, gene_b) {
          list(type = "Fusion", gene = gene_a, alteration = glue("{gene_a}::{gene_b}"))
        })
      }
  }

  all_key_alts <- unlist(key_alterations, recursive = FALSE)
  log_info("Found {length(all_key_alts)} key alterations for narrative")

  # Phase 4: Generate narratives for each key alteration
  log_info("Phase 4: Generating clinical narratives")

  narratives <- map_chr(all_key_alts, function(alt) {
    gene <- alt$gene
    alteration <- alt$alteration
    alt_type <- alt$type

    log_info("Generating narrative for {gene} {alteration}")

    # Get OncoKB annotation
    oncokb_data <- tryCatch(
      {
        if (alt_type == "SNV/Indel") {
          oncokb_annotate_mutation(gene, alteration, tumor_type)
        } else if (alt_type == "CNV") {
          # Extract CNA type from alteration
          cna_type <- sub(".*\\s", "", alteration)
          if (!cna_type %in% c("AMPLIFICATION", "DELETION")) {
            cna_type <- "AMPLIFICATION"
          }
          oncokb_annotate_cna(gene, cna_type, tumor_type)
        } else if (alt_type == "Fusion") {
          parts <- strsplit(alteration, "::")[[1]]
          if (length(parts) == 2) {
            oncokb_annotate_fusion(parts[1], parts[2], tumor_type)
          } else {
            list(oncogenic = "Unknown", treatments = tibble())
          }
        } else {
          list(oncogenic = "Unknown", treatments = tibble())
        }
      },
      error = function(e) {
        log_warn("OncoKB annotation failed for {gene} {alteration}: {conditionMessage(e)}")
        list(oncogenic = "Unknown", treatments = tibble())
      }
    )

    # Build oncogenic classification statement
    oncogenic <- oncokb_data$oncogenic %||% "Unknown"
    oncogenic_stmt <- switch(oncogenic,
      "Oncogenic" = "an oncogenic",
      "Likely Oncogenic" = "a likely oncogenic",
      "Likely Neutral" = "a likely neutral",
      "Neutral" = "a neutral",
      "Inconclusive" = "an inconclusive",
      "Unknown" = "an uncharacterized"
    )
    oncogenic_line <- glue("{gene} {alteration} is {oncogenic_stmt} alteration.")

    # Build frequency statement
    freq_statement <- ""
    pubmed_hits <- pubmed_results %>%
      filter(query_gene == gene | grepl(gene, title, ignore.case = TRUE))

    if (nrow(pubmed_hits) > 0) {
      freq_statement <- glue("This alteration has been reported in the literature with {nrow(pubmed_hits)} PubMed citations.")
    } else {
      freq_statement <- "Limited literature available for this specific alteration."
    }

    # Build therapeutic relevance statement
    therapy_statement <- ""
    if (nrow(oncokb_data$treatments %||% tibble()) > 0) {
      therapies <- oncokb_data$treatments %>%
        filter(fda_approved == TRUE) %>%
        pull(drugs) %>%
        paste(collapse = "; ")

      if (nchar(therapies) > 0) {
        therapy_statement <- glue("FDA-approved therapeutic options include {therapies}.")
      } else {
        therapies_all <- oncokb_data$treatments %>%
          pull(drugs) %>%
          paste(collapse = "; ")
        therapy_statement <- glue("Potential therapeutic options: {therapies_all}.")
      }
    } else {
      therapy_statement <- "No specific OncoKB-annotated therapies available for this alteration."
    }

    # Build key references section
    ref_statement <- ""

    # Combine PubMed and Scopus results for this alteration
    combined_hits <- bind_rows(
      pubmed_hits %>%
        select(title, authors, journal, year, pmid) %>%
        rename(article_id = pmid) %>%
        mutate(source = "PubMed"),
      if (nrow(scopus_results) > 0) {
        scopus_results %>%
          filter(query_gene == gene | grepl(gene, title, ignore.case = TRUE)) %>%
          select(title, authors, journal, year, scopus_id) %>%
          rename(article_id = scopus_id) %>%
          mutate(source = "Scopus")
      } else {
        tibble()
      }
    )

    if (nrow(combined_hits) > 0) {
      # Select top 3 most recent/relevant papers
      top_papers <- combined_hits %>%
        distinct(title, .keep_all = TRUE) %>%
        arrange(desc(year)) %>%
        slice_head(n = 3)

      refs <- top_papers %>%
        pmap_chr(function(title, authors, journal, year, ...) {
          # Extract first author
          first_author <- str_extract(authors, "^[^,]+") %||% "Anonymous"
          glue("{first_author} et al., {journal}, {year}")
        })

      ref_statement <- glue("Key references: {paste(refs, collapse = '; ')}.")
    } else {
      ref_statement <- "See PubMed and Scopus searches for comprehensive literature review."
    }

    # Combine all components
    narrative <- glue("{oncogenic_line} {freq_statement} {therapy_statement} {ref_statement}")

    log_info("Generated narrative ({nchar(narrative)} characters)")
    narrative
  })

  names(narratives) <- map_chr(all_key_alts, function(alt) {
    glue("{alt$gene}_{alt$alteration}")
  })

  log_info("Generated {length(narratives)} narratives")

  # Phase 5: Save combined results
  log_info("Phase 5: Saving results")

  output_dir <- glue("reports/{sample_id}/07-literature")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  combined_output <- list(
    sample_id = sample_id,
    tumor_type = tumor_type,
    generated_at = Sys.time(),
    pubmed_articles = if (nrow(pubmed_results) > 0) {
      pubmed_results %>%
        select(pmid, title, authors, journal, year, query_gene, query_variant)
    } else {
      tibble()
    },
    scopus_articles = if (nrow(scopus_results) > 0) {
      scopus_results %>%
        select(scopus_id, title, authors, journal, year, doi, cited_by, is_oa, oa_status)
    } else {
      tibble()
    },
    narratives = narratives
  )

  narratives_path <- file.path(output_dir, "narratives.json")
  tryCatch(
    {
      write_json(combined_output, narratives_path, pretty = TRUE)
      log_info("Saved narratives to {narratives_path}")
    },
    error = function(e) {
      log_error("Failed to save narratives: {conditionMessage(e)}")
    }
  )

  # Return structured result
  list(
    pubmed = pubmed_results,
    scopus = scopus_results,
    narratives = narratives,
    sample_id = sample_id,
    tumor_type = tumor_type
  )
}
