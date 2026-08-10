# Helper functions for the internal FASH/FASH-CL versus Strober comparison.

select_cumulative_lfdr_calls <- function(lfdr, alpha = 0.05) {
  lfdr <- as.numeric(lfdr)
  if (!length(lfdr) || any(!is.finite(lfdr)) ||
      any(lfdr < 0 | lfdr > 1) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("lfdr and alpha must be valid probabilities.")
  }
  ordering <- order(lfdr, method = "radix")
  accepted <- which(cumsum(lfdr[ordering]) / seq_along(ordering) <= alpha)
  if (!length(accepted)) {
    return(integer())
  }
  ordering[seq_len(max(accepted))]
}

parse_pair_keys <- function(pair_keys) {
  pair_keys <- as.character(pair_keys)
  if (!length(pair_keys) || anyNA(pair_keys) ||
      any(!grepl("^[^_]+_.+$", pair_keys)) || anyDuplicated(pair_keys)) {
    stop("Pair keys must be unique and use the gene_variant format.")
  }
  data.frame(
    pair_key = pair_keys,
    gene_id = sub("_.*$", "", pair_keys),
    variant_id = sub("^[^_]+_", "", pair_keys),
    stringsAsFactors = FALSE
  )
}

derive_ranked_discovery_sets <- function(pair_table, score) {
  required_columns <- c("pair_key", "gene_id", "variant_id")
  score <- as.numeric(score)
  if (!all(required_columns %in% names(pair_table)) ||
      nrow(pair_table) != length(score) || anyDuplicated(pair_table$pair_key) ||
      any(!is.finite(score))) {
    stop("The discovered pair table or ranking score is invalid.")
  }
  ranked_pairs <- pair_table
  ranked_pairs$score <- score
  ranked_pairs <- ranked_pairs[
    order(
      ranked_pairs$score,
      ranked_pairs$variant_id,
      ranked_pairs$pair_key,
      method = "radix"
    ),
    ,
    drop = FALSE
  ]
  row.names(ranked_pairs) <- NULL
  lead_pairs <- ranked_pairs[!duplicated(ranked_pairs$gene_id), , drop = FALSE]
  list(
    all_pairs = ranked_pairs,
    lead_pairs = lead_pairs,
    all_variants = unique(ranked_pairs$variant_id),
    lead_variants = unique(lead_pairs$variant_id)
  )
}

derive_fash_discovery_sets <- function(pair_keys, lfdr, alpha = 0.05) {
  pair_table <- parse_pair_keys(pair_keys)
  if (length(lfdr) != nrow(pair_table)) {
    stop("FASH pair keys and lfdr values are not aligned.")
  }
  discovered <- select_cumulative_lfdr_calls(lfdr, alpha = alpha)
  derive_ranked_discovery_sets(
    pair_table[discovered, , drop = FALSE],
    lfdr[discovered]
  )
}

derive_strober_discovery_sets <- function(result,
                                          pair_keys,
                                          alpha = 0.05) {
  required_columns <- c("rs_id", "ensamble_id", "pvalue", "eFDR")
  if (!identical(names(result), required_columns) ||
      any(!is.finite(result$pvalue)) || any(!is.finite(result$eFDR)) ||
      any(result$eFDR < 0 | result$eFDR > 1)) {
    stop("The Strober result table has an unexpected schema or values.")
  }
  result_keys <- paste(result$ensamble_id, result$rs_id, sep = "_")
  if (!setequal(result_keys, pair_keys) || anyDuplicated(result_keys)) {
    stop("The Strober result table does not match the tested pair universe.")
  }
  discovered <- result[result$eFDR <= alpha, , drop = FALSE]
  pair_table <- data.frame(
    pair_key = paste(discovered$ensamble_id, discovered$rs_id, sep = "_"),
    gene_id = as.character(discovered$ensamble_id),
    variant_id = as.character(discovered$rs_id),
    stringsAsFactors = FALSE
  )
  derive_ranked_discovery_sets(pair_table, discovered$pvalue)
}

compute_unmatched_enrichment <- function(annotation_matrix,
                                         selected_sets,
                                         annotations) {
  if (!"variant_id" %in% names(annotation_matrix) ||
      !all(annotations %in% names(annotation_matrix)) ||
      !length(selected_sets) || is.null(names(selected_sets))) {
    stop("The annotation matrix or selected sets are invalid.")
  }
  background_rate <- vapply(
    annotation_matrix[annotations],
    function(column) mean(as.logical(column)),
    numeric(1)
  )
  rows <- lapply(names(selected_sets), function(set_name) {
    covered <- intersect(
      unique(as.character(selected_sets[[set_name]])),
      annotation_matrix$variant_id
    )
    indices <- match(covered, annotation_matrix$variant_id)
    do.call(rbind, lapply(annotations, function(annotation) {
      overlap <- sum(as.logical(annotation_matrix[[annotation]][indices]))
      selected_rate <- if (length(indices)) {
        overlap / length(indices)
      } else {
        NA_real_
      }
      data.frame(
        discovery_set = set_name,
        annotation = annotation,
        unmatched_selected_total = length(indices),
        unmatched_background_rate = unname(background_rate[[annotation]]),
        unmatched_selected_rate = selected_rate,
        unmatched_enrichment =
          selected_rate / unname(background_rate[[annotation]]),
        stringsAsFactors = FALSE
      )
    }))
  })
  output <- do.call(rbind, rows)
  row.names(output) <- NULL
  output
}

summarize_discovery_set <- function(method,
                                    view,
                                    selection_strategy,
                                    pairs,
                                    variants) {
  data.frame(
    method = method,
    view = view,
    selection_strategy = selection_strategy,
    pair_count = nrow(pairs),
    unique_variant_count = length(unique(variants)),
    unique_gene_count = length(unique(pairs$gene_id)),
    stringsAsFactors = FALSE
  )
}

build_discovery_overlap <- function(selected_sets, set_metadata) {
  ordinary <- set_metadata$view == "Ordinary"
  set_names <- set_metadata$discovery_set[ordinary]
  strategies <- unique(set_metadata$selection_strategy[ordinary])
  rows <- list()
  output_index <- 0L
  for (strategy in strategies) {
    strategy_names <- set_metadata$discovery_set[
      ordinary & set_metadata$selection_strategy == strategy
    ]
    for (first_index in seq_along(strategy_names)) {
      for (second_index in seq_along(strategy_names)) {
        first <- unique(selected_sets[[strategy_names[first_index]]])
        second <- unique(selected_sets[[strategy_names[second_index]]])
        intersection_count <- length(intersect(first, second))
        union_count <- length(union(first, second))
        output_index <- output_index + 1L
        rows[[output_index]] <- data.frame(
          selection_strategy = strategy,
          first_set = strategy_names[first_index],
          second_set = strategy_names[second_index],
          intersection_count = intersection_count,
          union_count = union_count,
          jaccard = if (union_count) intersection_count / union_count else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

summarize_enhancer_panel <- function(enhancer_results) {
  required_columns <- c(
    "annotation_system", "discovery_set", "method", "view",
    "selection_strategy", "log2_enrichment", "enrichment", "p_value",
    "q_value_within_set"
  )
  if (!all(required_columns %in% names(enhancer_results))) {
    stop("Enhancer results lack required columns.")
  }
  groups <- split(
    seq_len(nrow(enhancer_results)),
    interaction(
      enhancer_results$annotation_system,
      enhancer_results$discovery_set,
      drop = TRUE,
      lex.order = TRUE
    )
  )
  rows <- lapply(groups, function(indices) {
    group <- enhancer_results[indices, , drop = FALSE]
    finite <- is.finite(group$log2_enrichment)
    data.frame(
      annotation_system = group$annotation_system[1L],
      discovery_set = group$discovery_set[1L],
      method = group$method[1L],
      view = group$view[1L],
      selection_strategy = group$selection_strategy[1L],
      annotation_count = nrow(group),
      finite_annotation_count = sum(finite),
      positive_annotation_count = sum(
        group$log2_enrichment > 0,
        na.rm = TRUE
      ),
      median_log2_enrichment = if (any(finite)) {
        stats::median(group$log2_enrichment[finite])
      } else {
        NA_real_
      },
      maximum_enrichment = if (any(finite)) {
        max(group$enrichment[finite])
      } else {
        NA_real_
      },
      minimum_p_value = suppressWarnings(min(group$p_value, na.rm = TRUE)),
      minimum_q_value = suppressWarnings(min(
        group$q_value_within_set,
        na.rm = TRUE
      )),
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, rows)
  for (column in c("minimum_p_value", "minimum_q_value")) {
    output[[column]][!is.finite(output[[column]])] <- NA_real_
  }
  row.names(output) <- NULL
  output
}
