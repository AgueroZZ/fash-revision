# Helper functions for the internal FASH versus FASH-CL enrichment comparison.

require_comparison_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The ", package, " package is required.")
  }
  invisible(TRUE)
}

row_log_sum_exp <- function(values) {
  values <- as.matrix(values)
  if (!is.numeric(values) || !nrow(values) || !ncol(values) ||
      any(is.na(values)) || any(is.nan(values))) {
    stop("values must be a non-empty numeric matrix without NA or NaN.")
  }
  row_maximum <- apply(values, 1L, max)
  if (any(!is.finite(row_maximum))) {
    stop("Every row must contain at least one finite value.")
  }
  row_maximum + log(rowSums(exp(values - row_maximum)))
}

compute_log_bayes_factor <- function(log_likelihood_matrix,
                                     alternative_weights,
                                     chunk_size = 50000L) {
  log_likelihood_matrix <- as.matrix(log_likelihood_matrix)
  chunk_size <- as.integer(chunk_size)
  if (!is.numeric(log_likelihood_matrix) ||
      nrow(log_likelihood_matrix) < 1L ||
      ncol(log_likelihood_matrix) < 2L ||
      any(!is.finite(log_likelihood_matrix)) ||
      length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L) {
    stop("The log-likelihood matrix or chunk size is invalid.")
  }
  alternative_weights <- as.numeric(alternative_weights)
  if (length(alternative_weights) != ncol(log_likelihood_matrix) - 1L ||
      any(!is.finite(alternative_weights)) ||
      any(alternative_weights < 0) || sum(alternative_weights) <= 0) {
    stop("Alternative weights do not match the log-likelihood matrix.")
  }
  alternative_weights <- alternative_weights / sum(alternative_weights)
  log_alternative_weights <- log(alternative_weights)
  log_bayes_factor <- numeric(nrow(log_likelihood_matrix))
  starts <- seq.int(1L, nrow(log_likelihood_matrix), by = chunk_size)
  for (start in starts) {
    rows <- seq.int(start, min(nrow(log_likelihood_matrix),
                              start + chunk_size - 1L))
    weighted_alternative <- sweep(
      log_likelihood_matrix[rows, -1L, drop = FALSE],
      2L,
      log_alternative_weights,
      "+"
    )
    log_bayes_factor[rows] <-
      row_log_sum_exp(weighted_alternative) - log_likelihood_matrix[rows, 1L]
  }
  log_bayes_factor
}

compute_bf_adjusted_lfdr <- function(log_likelihood_matrix,
                                     chunk_size = 50000L,
                                     verbose = FALSE) {
  require_comparison_namespace("mixsqp")
  log_likelihood_matrix <- as.matrix(log_likelihood_matrix)
  if (!is.numeric(log_likelihood_matrix) ||
      nrow(log_likelihood_matrix) < 1L ||
      ncol(log_likelihood_matrix) < 2L ||
      any(!is.finite(log_likelihood_matrix))) {
    stop("The log-likelihood matrix is invalid.")
  }

  started <- proc.time()[["elapsed"]]
  likelihood_matrix <- exp(log_likelihood_matrix)
  if (any(!is.finite(likelihood_matrix))) {
    stop("Exponentiated likelihood values are not finite.")
  }
  bf_mixture_fit <- mixsqp::mixsqp(
    L = likelihood_matrix,
    log = FALSE,
    control = list(verbose = isTRUE(verbose))
  )
  posterior_mixture_fit <- mixsqp::mixsqp(
    L = likelihood_matrix,
    log = FALSE,
    control = list(verbose = isTRUE(verbose))
  )
  rm(likelihood_matrix)
  invisible(gc())
  bf_alternative_weights <- as.numeric(bf_mixture_fit$x)[-1L]
  bf_alternative_weights <-
    bf_alternative_weights / sum(bf_alternative_weights)
  alternative_weights <- as.numeric(posterior_mixture_fit$x)[-1L]
  if (any(!is.finite(alternative_weights)) ||
      sum(alternative_weights) <= 0) {
    stop("The conditional alternative mixture weights are invalid.")
  }
  alternative_weights <- alternative_weights / sum(alternative_weights)
  bf_log_bayes_factor <- compute_log_bayes_factor(
    log_likelihood_matrix,
    bf_alternative_weights,
    chunk_size = chunk_size
  )
  posterior_log_bayes_factor <- compute_log_bayes_factor(
    log_likelihood_matrix,
    alternative_weights,
    chunk_size = chunk_size
  )

  bayes_factor <- exp(bf_log_bayes_factor)
  ordering <- order(bayes_factor, method = "radix")
  cumulative_bf_mean <- cumsum(bayes_factor[ordering]) / seq_along(ordering)
  crossing <- which(cumulative_bf_mean >= 1)[1L]
  null_weight <- if (is.na(crossing)) {
    1
  } else {
    crossing / length(bayes_factor)
  }
  adjusted_lfdr <- if (null_weight >= 1) {
    rep(1, length(posterior_log_bayes_factor))
  } else {
    stats::plogis(
      stats::qlogis(null_weight) - posterior_log_bayes_factor
    )
  }

  list(
    lfdr = adjusted_lfdr,
    null_weight = null_weight,
    bf_alternative_weights = bf_alternative_weights,
    alternative_weights = alternative_weights,
    log_bayes_factor = bf_log_bayes_factor,
    posterior_log_bayes_factor = posterior_log_bayes_factor,
    bf_mixture_status = as.character(bf_mixture_fit$status),
    posterior_mixture_status = as.character(posterior_mixture_fit$status),
    elapsed_seconds = proc.time()[["elapsed"]] - started
  )
}

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

parse_comparison_pair_keys <- function(pair_keys) {
  pair_keys <- as.character(pair_keys)
  if (!length(pair_keys) || anyNA(pair_keys) ||
      any(!grepl("^[^_]+_.+$", pair_keys))) {
    stop("Pair keys must use the gene_variant format.")
  }
  data.frame(
    pair_key = pair_keys,
    gene_id = sub("_.*$", "", pair_keys),
    variant_id = sub("^[^_]+_", "", pair_keys),
    stringsAsFactors = FALSE
  )
}

derive_all_and_lead_sets <- function(pair_keys, lfdr, alpha = 0.05) {
  pair_table <- parse_comparison_pair_keys(pair_keys)
  if (length(lfdr) != nrow(pair_table)) {
    stop("Pair keys and lfdr values are not aligned.")
  }
  discovered_indices <- select_cumulative_lfdr_calls(lfdr, alpha = alpha)
  discovered_pairs <- pair_table[discovered_indices, , drop = FALSE]
  discovered_pairs$lfdr <- as.numeric(lfdr[discovered_indices])
  discovered_pairs <- discovered_pairs[
    order(
      discovered_pairs$lfdr,
      discovered_pairs$variant_id,
      discovered_pairs$pair_key,
      method = "radix"
    ),
    ,
    drop = FALSE
  ]
  row.names(discovered_pairs) <- NULL
  lead_pairs <- discovered_pairs[
    !duplicated(discovered_pairs$gene_id),
    ,
    drop = FALSE
  ]
  list(
    all_variants = unique(discovered_pairs$variant_id),
    lead_variants = unique(lead_pairs$variant_id),
    discovered_pairs = discovered_pairs,
    lead_pairs = lead_pairs,
    summary = data.frame(
      selection_strategy = c("all", "one_lead_per_gene"),
      pair_count = c(nrow(discovered_pairs), nrow(lead_pairs)),
      unique_variant_count = c(
        length(unique(discovered_pairs$variant_id)),
        length(unique(lead_pairs$variant_id))
      ),
      unique_gene_count = c(
        length(unique(discovered_pairs$gene_id)),
        length(unique(lead_pairs$gene_id))
      ),
      stringsAsFactors = FALSE
    )
  )
}

summarize_method_overlap <- function(current_sets, fash_cl_sets) {
  required_names <- c("all", "one_lead_per_gene")
  if (!identical(names(current_sets), required_names) ||
      !identical(names(fash_cl_sets), required_names)) {
    stop("Method sets must contain all and one_lead_per_gene in order.")
  }
  do.call(rbind, lapply(required_names, function(strategy) {
    current <- unique(as.character(current_sets[[strategy]]))
    fash_cl <- unique(as.character(fash_cl_sets[[strategy]]))
    intersection_count <- length(intersect(current, fash_cl))
    union_count <- length(union(current, fash_cl))
    data.frame(
      selection_strategy = strategy,
      current_variant_count = length(current),
      fash_cl_variant_count = length(fash_cl),
      intersection_count = intersection_count,
      union_count = union_count,
      jaccard = intersection_count / union_count,
      stringsAsFactors = FALSE
    )
  }))
}

summarize_enhancer_panel <- function(enhancer_results) {
  required_columns <- c(
    "annotation_system", "method", "selection_strategy", "annotation",
    "log2_enrichment", "enrichment", "p_value", "q_value_within_set"
  )
  if (!all(required_columns %in% names(enhancer_results))) {
    stop("Enhancer results lack required columns.")
  }
  groups <- split(
    seq_len(nrow(enhancer_results)),
    interaction(
      enhancer_results$annotation_system,
      enhancer_results$method,
      enhancer_results$selection_strategy,
      drop = TRUE,
      lex.order = TRUE
    )
  )
  rows <- lapply(groups, function(indices) {
    group <- enhancer_results[indices, , drop = FALSE]
    finite <- is.finite(group$log2_enrichment)
    data.frame(
      annotation_system = group$annotation_system[1L],
      method = group$method[1L],
      selection_strategy = group$selection_strategy[1L],
      annotation_count = nrow(group),
      finite_annotation_count = sum(finite),
      median_log2_enrichment = if (any(finite)) {
        stats::median(group$log2_enrichment[finite])
      } else {
        NA_real_
      },
      positive_annotation_count = sum(group$log2_enrichment > 0, na.rm = TRUE),
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
  numeric_columns <- c("minimum_p_value", "minimum_q_value")
  for (column in numeric_columns) {
    output[[column]][!is.finite(output[[column]])] <- NA_real_
  }
  row.names(output) <- NULL
  output
}
