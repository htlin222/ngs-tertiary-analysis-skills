# R/civic_client.R — CiVIC GraphQL API client for community variant evidence
# CiVIC (Clinical Interpretation of Variants in Cancer)
# Endpoint: https://civicdb.org/api/graphql (free, no API key required)

suppressPackageStartupMessages({
  library(httr2)
  library(jsonlite)
  library(logger)
  library(glue)
  library(dplyr)
  library(purrr)
})

# ══════════════════════════════════════════════════════════════════════════════
# GraphQL Core
# ══════════════════════════════════════════════════════════════════════════════

#' Execute a GraphQL query against the CiVIC API
#'
#' @param query Character. GraphQL query string.
#' @param variables Named list of GraphQL variables (default: empty list).
#' @return Parsed JSON response body as a list.
#' @keywords internal
civic_graphql <- function(query, variables = list()) {
  body <- list(query = query)
  if (length(variables) > 0) body$variables <- variables

  resp <- request("https://civicdb.org/api/graphql") |>
    req_headers(
      `Content-Type` = "application/json",
      Accept = "application/json"
    ) |>
    req_body_json(body) |>
    req_retry(max_tries = 3, backoff = ~ 2^.x) |>
    req_throttle(rate = 30 / 60) |>
    req_perform()

  result <- resp_body_json(resp)

  if (!is.null(result$errors)) {
    error_msgs <- paste(map_chr(result$errors, ~ .x$message %||% "Unknown error"),
                        collapse = "; ")
    log_warn("CiVIC GraphQL errors: {error_msgs}")
  }

  result
}

# ══════════════════════════════════════════════════════════════════════════════
# Variant Search
# ══════════════════════════════════════════════════════════════════════════════

#' Search CiVIC for a variant by gene and variant name
#'
#' @param gene Character. Gene symbol (e.g., "BRAF").
#' @param variant Character. Variant name (e.g., "V600E").
#' @return Tibble with variant_id, gene, variant_name, evidence_count,
#'   molecular_profile_id. Empty tibble if no results.
#' @export
civic_search_variant <- function(gene, variant) {
  log_debug("CiVIC: searching variant {gene} {variant}")

  query <- '
    query SearchVariants($geneName: String!, $variantName: String) {
      variants(
        geneNames: [$geneName]
        variantName: $variantName
        first: 10
      ) {
        nodes {
          id
          name
          singleVariantMolecularProfile {
            id
            evidenceItems(first: 0) {
              totalCount
            }
          }
          feature {
            ... on Gene {
              name
            }
          }
        }
      }
    }
  '

  variables <- list(geneName = gene, variantName = variant)

  result <- tryCatch(
    civic_graphql(query, variables),
    error = function(e) {
      log_warn("CiVIC variant search failed for {gene} {variant}: {e$message}")
      return(NULL)
    }
  )

  empty_tbl <- tibble(
    variant_id = integer(), gene = character(),
    variant_name = character(), evidence_count = integer(),
    molecular_profile_id = integer()
  )

  if (is.null(result) || is.null(result$data$variants$nodes)) return(empty_tbl)
  nodes <- result$data$variants$nodes
  if (length(nodes) == 0) return(empty_tbl)

  map_dfr(nodes, function(node) {
    mp <- node$singleVariantMolecularProfile
    tibble(
      variant_id = as.integer(node$id %||% NA_integer_),
      gene = node$feature$name %||% gene,
      variant_name = node$name %||% NA_character_,
      evidence_count = as.integer(mp$evidenceItems$totalCount %||% 0L),
      molecular_profile_id = as.integer(mp$id %||% NA_integer_)
    )
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# Evidence Items
# ══════════════════════════════════════════════════════════════════════════════

#' Get evidence items for a variant from CiVIC
#'
#' @param variant_id Integer. CiVIC variant ID.
#' @param molecular_profile_id Integer. CiVIC molecular profile ID.
#' @return Tibble with evidence_id, evidence_type, evidence_level,
#'   evidence_direction, significance, disease, therapies, source_citation,
#'   source_url, description. Empty tibble if no results.
#' @export
civic_get_evidence <- function(variant_id, molecular_profile_id) {
  log_debug("CiVIC: fetching evidence for variant_id={variant_id}, mp_id={molecular_profile_id}")

  query <- '
    query GetEvidence($molecularProfileId: Int!) {
      molecularProfile(id: $molecularProfileId) {
        evidenceItems(first: 50) {
          nodes {
            id
            evidenceType
            evidenceLevel
            evidenceDirection
            significance
            disease { name }
            therapies { name }
            source { citation sourceUrl }
            description
            status
          }
        }
      }
    }
  '

  variables <- list(molecularProfileId = as.integer(molecular_profile_id))

  result <- tryCatch(
    civic_graphql(query, variables),
    error = function(e) {
      log_warn("CiVIC evidence fetch failed for mp_id={molecular_profile_id}: {e$message}")
      return(NULL)
    }
  )

  empty_tbl <- tibble(
    evidence_id = integer(), evidence_type = character(),
    evidence_level = character(), evidence_direction = character(),
    significance = character(), disease = character(),
    therapies = character(), source_citation = character(),
    source_url = character(), description = character()
  )

  if (is.null(result) || is.null(result$data$molecularProfile$evidenceItems$nodes)) {
    return(empty_tbl)
  }

  nodes <- result$data$molecularProfile$evidenceItems$nodes
  if (length(nodes) == 0) return(empty_tbl)

  # Filter to accepted evidence only
  nodes <- Filter(function(n) (n$status %||% "") == "accepted", nodes)
  if (length(nodes) == 0) return(empty_tbl)

  map_dfr(nodes, function(node) {
    therapy_names <- if (length(node$therapies) > 0) {
      paste(map_chr(node$therapies, ~ .x$name %||% ""), collapse = " + ")
    } else {
      NA_character_
    }
    tibble(
      evidence_id = as.integer(node$id %||% NA_integer_),
      evidence_type = node$evidenceType %||% NA_character_,
      evidence_level = node$evidenceLevel %||% NA_character_,
      evidence_direction = node$evidenceDirection %||% NA_character_,
      significance = node$significance %||% NA_character_,
      disease = node$disease$name %||% NA_character_,
      therapies = therapy_names,
      source_citation = node$source$citation %||% NA_character_,
      source_url = node$source$sourceUrl %||% NA_character_,
      description = node$description %||% NA_character_
    )
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# Gene Summary
# ══════════════════════════════════════════════════════════════════════════════

#' Get gene summary information from CiVIC
#'
#' @param gene Character. Gene symbol (e.g., "BRAF").
#' @return List with gene_id, gene, description, official_name.
#' @export
civic_get_gene_summary <- function(gene) {
  log_debug("CiVIC: fetching gene summary for {gene}")

  query <- '
    query GetGene($name: String!) {
      genes(name: $name) {
        nodes { id name description officialName }
      }
    }
  '

  result <- tryCatch(
    civic_graphql(query, list(name = gene)),
    error = function(e) {
      log_warn("CiVIC gene summary failed for {gene}: {e$message}")
      return(NULL)
    }
  )

  if (is.null(result) || is.null(result$data$genes$nodes) ||
      length(result$data$genes$nodes) == 0) {
    return(list(gene_id = NA_integer_, gene = gene,
                description = NA_character_, official_name = NA_character_))
  }

  node <- result$data$genes$nodes[[1]]
  list(
    gene_id = as.integer(node$id %||% NA_integer_),
    gene = node$name %||% gene,
    description = node$description %||% NA_character_,
    official_name = node$officialName %||% NA_character_
  )
}

# ══════════════════════════════════════════════════════════════════════════════
# AMP/ASCO/CAP Assertions
# ══════════════════════════════════════════════════════════════════════════════

#' Get AMP/ASCO/CAP assertions for a gene from CiVIC
#'
#' @param gene Character. Gene symbol (e.g., "BRAF").
#' @return Tibble with assertion_id, amp_level, assertion_type,
#'   assertion_direction, significance, disease, therapies, variant_name,
#'   gene, description, nccn_guideline, regulatory_approval, fda_companion_test.
#' @export
civic_get_assertions <- function(gene) {
  log_debug("CiVIC: fetching assertions for {gene}")

  query <- '
    query GetAssertions($geneName: String!) {
      assertions(geneNames: [$geneName], first: 50) {
        nodes {
          id
          ampLevel
          assertionType
          assertionDirection
          significance
          disease { name }
          therapies { name }
          molecularProfile { name }
          description
          nccnGuideline { name }
          regulatoryApproval
          fdaCompanionTest
          status
        }
      }
    }
  '

  result <- tryCatch(
    civic_graphql(query, list(geneName = gene)),
    error = function(e) {
      log_warn("CiVIC assertions fetch failed for {gene}: {e$message}")
      return(NULL)
    }
  )

  empty_tbl <- tibble(
    assertion_id = integer(), amp_level = character(),
    assertion_type = character(), assertion_direction = character(),
    significance = character(), disease = character(),
    therapies = character(), variant_name = character(),
    gene = character(), description = character(),
    nccn_guideline = character(), regulatory_approval = logical(),
    fda_companion_test = logical()
  )

  if (is.null(result) || is.null(result$data$assertions$nodes)) return(empty_tbl)
  nodes <- result$data$assertions$nodes
  if (length(nodes) == 0) return(empty_tbl)

  # Filter to accepted assertions only
  nodes <- Filter(function(n) (n$status %||% "") == "accepted", nodes)
  if (length(nodes) == 0) return(empty_tbl)

  map_dfr(nodes, function(node) {
    therapy_names <- if (length(node$therapies) > 0) {
      paste(map_chr(node$therapies, ~ .x$name %||% ""), collapse = " + ")
    } else {
      NA_character_
    }
    tibble(
      assertion_id = as.integer(node$id %||% NA_integer_),
      amp_level = node$ampLevel %||% NA_character_,
      assertion_type = node$assertionType %||% NA_character_,
      assertion_direction = node$assertionDirection %||% NA_character_,
      significance = node$significance %||% NA_character_,
      disease = node$disease$name %||% NA_character_,
      therapies = therapy_names,
      variant_name = node$molecularProfile$name %||% NA_character_,
      gene = gene,
      description = node$description %||% NA_character_,
      nccn_guideline = node$nccnGuideline$name %||% NA_character_,
      regulatory_approval = node$regulatoryApproval %||% FALSE,
      fda_companion_test = node$fdaCompanionTest %||% FALSE
    )
  })
}

# ══════════════════════════════════════════════════════════════════════════════
# Batch Annotation
# ══════════════════════════════════════════════════════════════════════════════

#' Batch annotate variants with CiVIC evidence
#'
#' Searches each variant in CiVIC, retrieves evidence items, gene summaries,
#' and AMP/ASCO/CAP assertions.
#'
#' @param variants Tibble with at least columns: gene, variant. Can be empty.
#' @param sample_id Character. Sample identifier for logging.
#' @return List with variant_evidence, gene_summaries, assertions.
#' @export
civic_annotate_variants <- function(variants, sample_id = "unknown") {
  log_info("CiVIC: annotating {nrow(variants)} variants for sample {sample_id}")

  empty_evidence <- tibble(
    evidence_id = integer(), evidence_type = character(),
    evidence_level = character(), evidence_direction = character(),
    significance = character(), disease = character(),
    therapies = character(), source_citation = character(),
    source_url = character(), description = character(),
    query_gene = character(), query_variant = character()
  )
  empty_assertions <- tibble(
    assertion_id = integer(), amp_level = character(),
    assertion_type = character(), assertion_direction = character(),
    significance = character(), disease = character(),
    therapies = character(), variant_name = character(),
    gene = character(), description = character(),
    nccn_guideline = character(), regulatory_approval = logical(),
    fda_companion_test = logical()
  )
  empty_result <- list(
    variant_evidence = empty_evidence,
    gene_summaries = list(),
    assertions = empty_assertions
  )

  if (nrow(variants) == 0) {
    log_info("CiVIC: no variants to annotate")
    return(empty_result)
  }

  unique_genes <- unique(variants$gene)
  unique_genes <- unique_genes[!is.na(unique_genes)]

  # Search each variant and collect evidence
  all_evidence <- map_dfr(seq_len(nrow(variants)), function(i) {
    gene <- variants$gene[i]
    variant <- variants$variant[i]
    if (is.na(gene) || is.na(variant)) return(tibble())

    search_result <- tryCatch(
      civic_search_variant(gene, variant),
      error = function(e) {
        log_warn("CiVIC search failed for {gene} {variant}: {e$message}")
        return(tibble())
      }
    )
    if (nrow(search_result) == 0) return(tibble())

    best <- search_result[1, ]
    if (is.na(best$molecular_profile_id)) return(tibble())

    evidence <- tryCatch(
      civic_get_evidence(best$variant_id, best$molecular_profile_id),
      error = function(e) {
        log_warn("CiVIC evidence fetch failed for {gene} {variant}: {e$message}")
        return(tibble())
      }
    )
    if (nrow(evidence) > 0) {
      evidence <- evidence |> mutate(query_gene = gene, query_variant = variant)
    }
    evidence
  })

  if (nrow(all_evidence) == 0) all_evidence <- empty_evidence

  # Gene summaries
  gene_summaries <- map(unique_genes, function(g) {
    tryCatch(civic_get_gene_summary(g), error = function(e) {
      log_warn("CiVIC gene summary failed for {g}: {e$message}")
      list(gene_id = NA_integer_, gene = g,
           description = NA_character_, official_name = NA_character_)
    })
  })
  names(gene_summaries) <- unique_genes

  # Assertions for all genes
  all_assertions <- map_dfr(unique_genes, function(g) {
    tryCatch(civic_get_assertions(g), error = function(e) {
      log_warn("CiVIC assertions failed for {g}: {e$message}")
      return(tibble())
    })
  })
  if (nrow(all_assertions) == 0) all_assertions <- empty_assertions

  log_info("CiVIC annotation complete: {nrow(all_evidence)} evidence items, ",
           "{length(gene_summaries)} gene summaries, {nrow(all_assertions)} assertions")

  list(
    variant_evidence = all_evidence,
    gene_summaries = gene_summaries,
    assertions = all_assertions
  )
}
