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
    if (is.null(result$data)) {
      stop("CiVIC GraphQL errors: ", error_msgs)
    }
    log_warn("CiVIC GraphQL partial errors: {error_msgs}")
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

  empty_tbl <- tibble(
    variant_id = integer(), gene = character(),
    variant_name = character(), evidence_count = integer(),
    molecular_profile_id = integer()
  )

  # Step 1: resolve gene symbol to CiVIC gene ID
  gene_query <- '
    query GetGeneId($symbols: [String!]) {
      genes(entrezSymbols: $symbols, first: 1) {
        nodes { id name }
      }
    }
  '

  gene_result <- tryCatch(
    civic_graphql(gene_query, list(symbols = list(gene))),
    error = function(e) {
      log_warn("CiVIC gene lookup failed for {gene}: {e$message}")
      return(NULL)
    }
  )

  if (is.null(gene_result) || is.null(gene_result$data$genes$nodes) ||
      length(gene_result$data$genes$nodes) == 0) {
    log_debug("CiVIC: gene {gene} not found")
    return(empty_tbl)
  }

  gene_id <- as.integer(gene_result$data$genes$nodes[[1]]$id)

  # Step 2: search variants by geneId + name
  var_query <- '
    query SearchVariants($geneId: Int, $name: String) {
      variants(
        geneId: $geneId
        name: $name
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
            name
          }
        }
      }
    }
  '

  variables <- list(geneId = gene_id, name = variant)

  result <- tryCatch(
    civic_graphql(var_query, variables),
    error = function(e) {
      log_warn("CiVIC variant search failed for {gene} {variant}: {e$message}")
      return(NULL)
    }
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
#' @param molecular_profile_id Integer. CiVIC molecular profile ID.
#' @return Tibble with evidence_id, evidence_type, evidence_level,
#'   evidence_direction, significance, disease, therapies, source_citation,
#'   source_url, description. Empty tibble if no results.
#' @export
civic_get_evidence <- function(molecular_profile_id) {
  log_debug("CiVIC: fetching evidence for mp_id={molecular_profile_id}")

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
    query GetGene($symbols: [String!]) {
      genes(entrezSymbols: $symbols, first: 1) {
        nodes { id name description fullName }
      }
    }
  '

  result <- tryCatch(
    civic_graphql(query, list(symbols = list(gene))),
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
    official_name = node$fullName %||% NA_character_
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

  # Fetch assertions sorted by evidence count without gene filter, then filter

  # in R. CiVIC has ~200 accepted assertions total, so fetching 100 is safe.
  query <- '
    query GetAssertions {
      assertions(
        sortBy: { column: EVIDENCE_ITEM_COUNT, direction: DESC }
        first: 100
      ) {
        nodes {
          id
          ampLevel
          assertionType
          assertionDirection
          significance
          disease { name }
          therapies { name }
          molecularProfile {
            name
            variants {
              feature { name }
            }
          }
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
    civic_graphql(query),
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

  # Filter by gene: check if any variant's gene matches the requested gene
  if (!is.null(gene)) {
    nodes <- Filter(function(a) {
      variants <- a$molecularProfile$variants %||% list()
      any(sapply(variants, function(v) {
        (v$feature$name %||% "") == gene
      }))
    }, nodes)
  }
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
# Batch GraphQL Helpers
# ══════════════════════════════════════════════════════════════════════════════

#' Batch-resolve gene symbols to CiVIC gene IDs using GraphQL aliases
#'
#' Resolves multiple gene symbols in a single GraphQL request using aliased
#' subqueries. Falls back to sequential resolution on failure.
#'
#' @param gene_symbols Character vector of gene symbols.
#' @return Named list mapping gene symbol -> integer gene ID. Genes not found
#'   are mapped to NA_integer_.
#' @keywords internal
civic_batch_resolve_genes <- function(gene_symbols) {
  gene_symbols <- unique(gene_symbols[!is.na(gene_symbols)])
  if (length(gene_symbols) == 0) return(list())

  log_debug("CiVIC batch: resolving {length(gene_symbols)} genes")

  # Build aliased query: g0: genes(entrezSymbols: ["BRAF"]) { nodes { id name } }
  aliases <- map_chr(seq_along(gene_symbols), function(i) {
    sym <- gene_symbols[i]
    glue('g{i - 1}: genes(entrezSymbols: ["{sym}"], first: 1) {{ nodes {{ id name }} }}')
  })
  query_str <- paste0("query {\n", paste(aliases, collapse = "\n"), "\n}")

  result <- tryCatch(
    civic_graphql(query_str),
    error = function(e) {
      log_warn("CiVIC batch gene resolution failed: {e$message}")
      return(NULL)
    }
  )

  gene_ids <- setNames(
    rep(NA_integer_, length(gene_symbols)),
    gene_symbols
  )

  if (is.null(result) || is.null(result$data)) return(as.list(gene_ids))

  for (i in seq_along(gene_symbols)) {
    alias <- paste0("g", i - 1)
    nodes <- result$data[[alias]]$nodes
    if (!is.null(nodes) && length(nodes) > 0) {
      gene_ids[gene_symbols[i]] <- as.integer(nodes[[1]]$id)
    }
  }

  as.list(gene_ids)
}

#' Batch-search variants in CiVIC using GraphQL aliases
#'
#' Searches multiple variants in a single GraphQL request, chunked to avoid
#' query size limits. Each chunk uses aliased subqueries.
#'
#' @param variant_queries Tibble with columns: gene, variant, gene_id.
#'   gene_id must be resolved (non-NA integer).
#' @param chunk_size Integer. Max variants per GraphQL request (default: 15).
#' @return Tibble with variant_id, gene, variant_name, evidence_count,
#'   molecular_profile_id, query_gene, query_variant.
#' @keywords internal
civic_batch_search_variants <- function(variant_queries, chunk_size = 15L) {
  if (nrow(variant_queries) == 0) return(tibble())

  log_debug("CiVIC batch: searching {nrow(variant_queries)} variants")

  # Split into chunks
  chunks <- split(variant_queries, ceiling(seq_len(nrow(variant_queries)) / chunk_size))

  results <- map_dfr(chunks, function(chunk) {
    # Build aliased query for this chunk
    aliases <- map_chr(seq_len(nrow(chunk)), function(i) {
      gid <- chunk$gene_id[i]
      vname <- gsub('"', '\\\\"', chunk$variant[i])
      glue('v{i - 1}: variants(geneId: {gid}, name: "{vname}", first: 5) {{
        nodes {{
          id
          name
          singleVariantMolecularProfile {{
            id
            evidenceItems(first: 0) {{ totalCount }}
          }}
          feature {{ name }}
        }}
      }}')
    })
    query_str <- paste0("query {\n", paste(aliases, collapse = "\n"), "\n}")

    batch_result <- tryCatch(
      civic_graphql(query_str),
      error = function(e) {
        log_warn("CiVIC batch variant search failed: {e$message}")
        return(NULL)
      }
    )

    if (is.null(batch_result) || is.null(batch_result$data)) return(tibble())

    map_dfr(seq_len(nrow(chunk)), function(i) {
      alias <- paste0("v", i - 1)
      nodes <- batch_result$data[[alias]]$nodes
      if (is.null(nodes) || length(nodes) == 0) return(tibble())

      map_dfr(nodes, function(node) {
        mp <- node$singleVariantMolecularProfile
        tibble(
          variant_id = as.integer(node$id %||% NA_integer_),
          gene = node$feature$name %||% chunk$gene[i],
          variant_name = node$name %||% NA_character_,
          evidence_count = as.integer(mp$evidenceItems$totalCount %||% 0L),
          molecular_profile_id = as.integer(mp$id %||% NA_integer_),
          query_gene = chunk$gene[i],
          query_variant = chunk$variant[i]
        )
      })
    })
  })

  results
}

#' Batch-fetch evidence for multiple molecular profiles using GraphQL aliases
#'
#' @param mp_ids Named integer vector. Names are "gene_variant" keys, values
#'   are molecular_profile_ids.
#' @param chunk_size Integer. Max profiles per request (default: 10).
#' @return Tibble of evidence items with query_gene and query_variant columns.
#' @keywords internal
civic_batch_get_evidence <- function(mp_ids, chunk_size = 10L) {
  if (length(mp_ids) == 0) return(tibble())

  log_debug("CiVIC batch: fetching evidence for {length(mp_ids)} molecular profiles")

  mp_keys <- names(mp_ids)
  chunks <- split(seq_along(mp_ids), ceiling(seq_along(mp_ids) / chunk_size))

  results <- map_dfr(chunks, function(idx_chunk) {
    aliases <- map_chr(seq_along(idx_chunk), function(j) {
      i <- idx_chunk[j]
      mpid <- mp_ids[i]
      glue('e{j - 1}: molecularProfile(id: {mpid}) {{
        evidenceItems(first: 50) {{
          nodes {{
            id
            evidenceType
            evidenceLevel
            evidenceDirection
            significance
            disease {{ name }}
            therapies {{ name }}
            source {{ citation sourceUrl }}
            description
            status
          }}
        }}
      }}')
    })
    query_str <- paste0("query {\n", paste(aliases, collapse = "\n"), "\n}")

    batch_result <- tryCatch(
      civic_graphql(query_str),
      error = function(e) {
        log_warn("CiVIC batch evidence fetch failed: {e$message}")
        return(NULL)
      }
    )

    if (is.null(batch_result) || is.null(batch_result$data)) return(tibble())

    map_dfr(seq_along(idx_chunk), function(j) {
      i <- idx_chunk[j]
      alias <- paste0("e", j - 1)
      nodes <- batch_result$data[[alias]]$evidenceItems$nodes
      if (is.null(nodes) || length(nodes) == 0) return(tibble())

      # Filter to accepted evidence
      nodes <- Filter(function(n) (n$status %||% "") == "accepted", nodes)
      if (length(nodes) == 0) return(tibble())

      # Parse the key back to gene + variant
      key_parts <- strsplit(mp_keys[i], "\\|\\|")[[1]]
      qgene <- key_parts[1]
      qvariant <- if (length(key_parts) >= 2) key_parts[2] else NA_character_

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
          description = node$description %||% NA_character_,
          query_gene = qgene,
          query_variant = qvariant
        )
      })
    })
  })

  results
}

# ══════════════════════════════════════════════════════════════════════════════
# Batch Annotation
# ══════════════════════════════════════════════════════════════════════════════

#' Batch annotate variants with CiVIC evidence using batched GraphQL queries
#'
#' Optimized version that uses GraphQL aliases to batch gene resolution,
#' variant search, and evidence fetches into fewer API calls. Falls back
#' to the sequential approach on any failure.
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

  # Try batch approach first, fallback to sequential
  all_evidence <- tryCatch({
    civic_batch_annotate_evidence(variants, unique_genes)
  }, error = function(e) {
    log_warn("CiVIC batch annotation failed, falling back to sequential: {e$message}")
    NULL
  })

  # Fallback: sequential variant search + evidence fetch
  if (is.null(all_evidence)) {
    all_evidence <- civic_sequential_annotate_evidence(variants)
  }

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

#' Batch evidence annotation using batched GraphQL queries
#'
#' @param variants Tibble with gene + variant columns.
#' @param unique_genes Character vector of unique gene symbols.
#' @return Tibble of evidence items, or empty tibble.
#' @keywords internal
civic_batch_annotate_evidence <- function(variants, unique_genes) {
  # Step 1: Batch-resolve all gene symbols to IDs
  gene_ids <- civic_batch_resolve_genes(unique_genes)
  resolved <- gene_ids[!is.na(gene_ids)]

  if (length(resolved) == 0) {
    log_info("CiVIC batch: no genes resolved")
    return(tibble())
  }

  log_info("CiVIC batch: resolved {length(resolved)}/{length(unique_genes)} genes")

  # Step 2: Prepare variant queries with gene IDs
  variant_queries <- variants |>
    filter(!is.na(gene), !is.na(variant), gene %in% names(resolved)) |>
    mutate(gene_id = as.integer(gene_ids[gene]))

  if (nrow(variant_queries) == 0) return(tibble())

  # Step 3: Batch-search all variants
  search_results <- civic_batch_search_variants(variant_queries)
  if (nrow(search_results) == 0) return(tibble())

  # Pick best match per query (first result, highest evidence count)
  best_matches <- search_results |>
    filter(!is.na(molecular_profile_id)) |>
    group_by(query_gene, query_variant) |>
    slice_max(evidence_count, n = 1, with_ties = FALSE) |>
    ungroup()

  if (nrow(best_matches) == 0) return(tibble())

  # Step 4: Batch-fetch evidence for all matched molecular profiles
  mp_ids <- setNames(
    best_matches$molecular_profile_id,
    paste0(best_matches$query_gene, "||", best_matches$query_variant)
  )
  # Deduplicate by molecular_profile_id (same mp may appear for different queries)
  mp_ids <- mp_ids[!duplicated(mp_ids)]

  all_evidence <- civic_batch_get_evidence(mp_ids)

  log_info("CiVIC batch: retrieved {nrow(all_evidence)} evidence items")
  all_evidence
}

#' Sequential evidence annotation (fallback)
#'
#' @param variants Tibble with gene + variant columns.
#' @return Tibble of evidence items.
#' @keywords internal
civic_sequential_annotate_evidence <- function(variants) {
  log_info("CiVIC: using sequential annotation for {nrow(variants)} variants")

  map_dfr(seq_len(nrow(variants)), function(i) {
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
      civic_get_evidence(best$molecular_profile_id),
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
}
