# Helpers for the R5 balanced-variant thinning sensitivity analysis.

validate_pair_keys <- function(pair_keys) {
  pair_keys <- as.character(pair_keys)
  if (length(pair_keys) < 2L || anyNA(pair_keys) || any(!nzchar(pair_keys)) ||
      anyDuplicated(pair_keys) || any(!grepl("_", pair_keys, fixed = TRUE))) {
    stop("pair_keys must contain at least two unique gene-variant keys.")
  }
  gene_id <- sub("_.*$", "", pair_keys)
  variant_id <- sub("^[^_]+_", "", pair_keys)
  if (any(!nzchar(gene_id)) || any(!nzchar(variant_id))) {
    stop("Every pair key must contain a non-empty gene and variant ID.")
  }
  list(pair_keys = pair_keys, gene_id = gene_id, variant_id = variant_id)
}

make_gene_index <- function(pair_keys) {
  parsed <- validate_pair_keys(pair_keys)
  split(
    seq_along(parsed$pair_keys),
    factor(
      parsed$gene_id,
      levels = sort(unique(parsed$gene_id), method = "radix")
    )
  )
}

validate_gene_index <- function(gene_index, pair_keys) {
  parsed <- validate_pair_keys(pair_keys)
  flattened <- unlist(gene_index, use.names = FALSE)
  valid <-
    is.list(gene_index) &&
    length(gene_index) >= 2L &&
    all(nzchar(names(gene_index))) &&
    anyDuplicated(names(gene_index)) == 0L &&
    all(lengths(gene_index) >= 1L) &&
    length(flattened) == length(parsed$pair_keys) &&
    !anyNA(flattened) &&
    all(flattened >= 1L & flattened <= length(parsed$pair_keys)) &&
    anyDuplicated(flattened) == 0L &&
    identical(sort(as.integer(flattened)), seq_along(parsed$pair_keys)) &&
    all(vapply(names(gene_index), function(gene) {
      indices <- gene_index[[gene]]
      all(parsed$gene_id[indices] == gene)
    }, logical(1)))
  if (!isTRUE(valid)) {
    stop("gene_index must partition every pair key into its named gene group.")
  }
  invisible(TRUE)
}

select_balanced_variants_per_gene <- function(pair_keys,
                                              seed,
                                              target_per_gene,
                                              gene_index = NULL) {
  parsed <- validate_pair_keys(pair_keys)
  seed <- as.integer(seed)
  target_per_gene <- as.integer(target_per_gene)
  if (length(seed) != 1L || is.na(seed) || seed < 0L ||
      length(target_per_gene) != 1L || is.na(target_per_gene) ||
      target_per_gene < 1L) {
    stop("seed and target_per_gene must be non-negative and positive integers.")
  }
  if (is.null(gene_index)) {
    gene_index <- make_gene_index(parsed$pair_keys)
  }
  validate_gene_index(gene_index, parsed$pair_keys)
  eligible_gene_index <- gene_index[lengths(gene_index) >= target_per_gene]
  if (length(eligible_gene_index) < 2L) {
    stop("At least two genes must meet target_per_gene.")
  }

  set.seed(seed)
  selected_by_gene <- lapply(eligible_gene_index, function(indices) {
    sort(sample(indices, size = target_per_gene, replace = FALSE))
  })
  selected_indices <- as.integer(unlist(selected_by_gene, use.names = FALSE))
  selected_gene_id <- rep(
    names(selected_by_gene),
    times = lengths(selected_by_gene)
  )
  expected_rows <- length(eligible_gene_index) * target_per_gene
  if (length(selected_indices) != expected_rows ||
      anyDuplicated(selected_indices) ||
      any(parsed$gene_id[selected_indices] != selected_gene_id) ||
      any(table(selected_gene_id) != target_per_gene)) {
    stop("Balanced sampling did not return the requested pairs per gene.")
  }

  data.frame(
    seed = seed,
    target_per_gene = target_per_gene,
    fash_index = selected_indices,
    pair_key = parsed$pair_keys[selected_indices],
    gene_id = selected_gene_id,
    variant_id = parsed$variant_id[selected_indices],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

validate_probability_vector <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) < 2L || any(!is.finite(x)) || any(x < 0 | x > 1)) {
    stop(name, " must contain at least two finite probabilities.")
  }
  x
}

cumulative_fdr_calls <- function(lfdr, alpha = 0.05) {
  lfdr <- validate_probability_vector(lfdr, "lfdr")
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be one finite number strictly between zero and one.")
  }
  ordering <- order(lfdr, method = "radix")
  cumulative_fdr <- cumsum(lfdr[ordering]) / seq_along(ordering)
  selected_rank <- which(cumulative_fdr <= alpha)
  if (length(selected_rank) == 0L) {
    return(integer(0))
  }
  ordering[selected_rank]
}

validate_prior_weights <- function(prior_weights, name) {
  prior_weights <- as.data.frame(prior_weights, stringsAsFactors = FALSE)
  if (!all(c("psd", "prior_weight") %in% names(prior_weights))) {
    stop(name, " must contain psd and prior_weight columns.")
  }
  prior_weights <- prior_weights[, c("psd", "prior_weight")]
  prior_weights$psd <- as.numeric(prior_weights$psd)
  prior_weights$prior_weight <- as.numeric(prior_weights$prior_weight)
  if (nrow(prior_weights) < 2L ||
      any(!is.finite(as.matrix(prior_weights))) ||
      any(prior_weights$psd < 0) ||
      any(prior_weights$prior_weight < 0) ||
      anyDuplicated(prior_weights$psd) ||
      abs(sum(prior_weights$prior_weight) - 1) > 1e-6 ||
      sum(prior_weights$psd == 0) != 1L) {
    stop(name, " is not a valid normalized FASH prior-weight table.")
  }
  prior_weights[order(prior_weights$psd), , drop = FALSE]
}

compare_prior_weights <- function(full_prior, thinned_prior) {
  full_prior <- validate_prior_weights(full_prior, "full_prior")
  thinned_prior <- validate_prior_weights(thinned_prior, "thinned_prior")
  comparison <- merge(
    full_prior,
    thinned_prior,
    by = "psd",
    all = TRUE,
    suffixes = c("_full", "_thinned"),
    sort = TRUE
  )
  names(comparison)[names(comparison) == "prior_weight_full"] <- "full_weight"
  names(comparison)[names(comparison) == "prior_weight_thinned"] <-
    "thinned_weight"
  comparison$full_weight[is.na(comparison$full_weight)] <- 0
  comparison$thinned_weight[is.na(comparison$thinned_weight)] <- 0
  comparison$difference <- comparison$thinned_weight - comparison$full_weight
  comparison$absolute_difference <- abs(comparison$difference)
  null_row <- which(comparison$psd == 0)
  list(
    table = comparison,
    summary = data.frame(
      full_pi0 = comparison$full_weight[null_row],
      thinned_pi0 = comparison$thinned_weight[null_row],
      pi0_difference = comparison$difference[null_row],
      prior_total_variation = 0.5 * sum(comparison$absolute_difference),
      stringsAsFactors = FALSE
    )
  )
}

safe_correlation <- function(x, y, method) {
  if (stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }
  unname(stats::cor(x, y, method = method))
}

compare_paired_lfdr <- function(full_lfdr,
                                thinned_lfdr,
                                pair_keys,
                                alpha = 0.05) {
  full_lfdr <- validate_probability_vector(full_lfdr, "full_lfdr")
  thinned_lfdr <- validate_probability_vector(thinned_lfdr, "thinned_lfdr")
  pair_keys <- as.character(pair_keys)
  if (length(full_lfdr) != length(thinned_lfdr) ||
      length(pair_keys) != length(full_lfdr) ||
      any(!nzchar(pair_keys)) || anyDuplicated(pair_keys)) {
    stop("Paired lfdr vectors must align with unique pair keys.")
  }
  full_calls <- cumulative_fdr_calls(full_lfdr, alpha = alpha)
  thinned_calls <- cumulative_fdr_calls(thinned_lfdr, alpha = alpha)
  difference <- thinned_lfdr - full_lfdr
  call_union <- union(full_calls, thinned_calls)
  call_intersection <- intersect(full_calls, thinned_calls)
  jaccard <- if (length(call_union) == 0L) {
    1
  } else {
    length(call_intersection) / length(call_union)
  }
  list(
    table = data.frame(
      pair_key = pair_keys,
      full_lfdr = full_lfdr,
      thinned_lfdr = thinned_lfdr,
      lfdr_difference = difference,
      absolute_difference = abs(difference),
      full_fdr_call = seq_along(full_lfdr) %in% full_calls,
      thinned_fdr_call = seq_along(thinned_lfdr) %in% thinned_calls,
      stringsAsFactors = FALSE
    ),
    summary = data.frame(
      n_units = length(full_lfdr),
      pearson_lfdr = safe_correlation(full_lfdr, thinned_lfdr, "pearson"),
      spearman_lfdr = safe_correlation(full_lfdr, thinned_lfdr, "spearman"),
      mean_absolute_lfdr_difference = mean(abs(difference)),
      median_absolute_lfdr_difference = stats::median(abs(difference)),
      rmse_lfdr = sqrt(mean(difference^2)),
      full_fdr_calls = length(full_calls),
      thinned_fdr_calls = length(thinned_calls),
      fdr_call_intersection = length(call_intersection),
      fdr_call_union = length(call_union),
      fdr_call_jaccard = jaccard,
      alpha = alpha,
      stringsAsFactors = FALSE
    ),
    full_call_indices = full_calls,
    thinned_call_indices = thinned_calls
  )
}

refit_bf_from_cached_likelihood <- function(full_likelihood,
                                             selected_indices,
                                             selected_pair_keys,
                                             psd_grid,
                                             penalty) {
  if (!requireNamespace("fashr", quietly = TRUE)) {
    stop("The fashr package is required.")
  }
  selected_indices <- as.integer(selected_indices)
  selected_pair_keys <- as.character(selected_pair_keys)
  psd_grid <- as.numeric(psd_grid)
  penalty <- as.numeric(penalty)
  if (!is.matrix(full_likelihood) || nrow(full_likelihood) < 2L ||
      ncol(full_likelihood) != length(psd_grid) ||
      anyNA(selected_indices) || anyDuplicated(selected_indices) ||
      any(selected_indices < 1L | selected_indices > nrow(full_likelihood)) ||
      length(selected_pair_keys) != length(selected_indices) ||
      any(!nzchar(selected_pair_keys)) || anyDuplicated(selected_pair_keys) ||
      length(penalty) != 1L || !is.finite(penalty) || penalty < 1) {
    stop("Cached-likelihood BF refit inputs are invalid.")
  }
  selected_likelihood <- full_likelihood[selected_indices, , drop = FALSE]
  if (anyNA(selected_likelihood) || any(is.nan(selected_likelihood)) ||
      any(selected_likelihood == Inf)) {
    stop("The selected likelihood matrix contains invalid values.")
  }
  rownames(selected_likelihood) <- selected_pair_keys
  eb_result <- fashr::fash_eb_est(
    L_matrix = selected_likelihood,
    grid = psd_grid,
    penalty = penalty
  )
  rownames(eb_result$posterior_weight) <- selected_pair_keys
  null_column <- which(eb_result$prior_weight$psd == 0)
  if (length(null_column) != 1L) {
    stop("The raw thinned refit does not contain exactly one null component.")
  }
  minimal_fit <- structure(
    list(
      prior_weights = eb_result$prior_weight,
      posterior_weights = eb_result$posterior_weight,
      psd_grid = psd_grid,
      lfdr = eb_result$posterior_weight[, null_column],
      L_matrix = selected_likelihood
    ),
    class = "fash"
  )
  fashr::BF_update(minimal_fit, plot = FALSE)
}

cumulative_gene_union <- function(discovered_gene_sets,
                                  seeds,
                                  full_discovered_genes) {
  if (!is.list(discovered_gene_sets) ||
      length(discovered_gene_sets) != length(seeds) ||
      length(seeds) < 2L || anyDuplicated(seeds)) {
    stop("discovered_gene_sets and seeds must define aligned unique replicates.")
  }
  full_discovered_genes <- unique(as.character(full_discovered_genes))
  running_union <- character()
  rows <- vector("list", length(seeds))
  for (index in seq_along(seeds)) {
    genes <- unique(as.character(discovered_gene_sets[[index]]))
    if (anyNA(genes) || any(!nzchar(genes))) {
      stop("Every discovered-gene set must contain non-empty identifiers.")
    }
    running_union <- union(running_union, genes)
    rows[[index]] <- data.frame(
      n_seeds = index,
      added_seed = seeds[index],
      seed_discovered_genes = length(genes),
      cumulative_unique_genes = length(running_union),
      cumulative_full_genes_recovered = length(intersect(
        running_union,
        full_discovered_genes
      )),
      stringsAsFactors = FALSE
    )
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  if (any(diff(result$cumulative_unique_genes) < 0L) ||
      any(diff(result$cumulative_full_genes_recovered) < 0L)) {
    stop("Cumulative discovered-gene counts must be nondecreasing.")
  }
  result
}

seed_subset_union_envelope <- function(discovered_gene_sets,
                                       seeds,
                                       full_discovered_genes) {
  if (!is.list(discovered_gene_sets) ||
      length(discovered_gene_sets) != length(seeds) ||
      length(seeds) < 2L || anyDuplicated(seeds)) {
    stop("discovered_gene_sets and seeds must define aligned unique replicates.")
  }
  discovered_gene_sets <- lapply(discovered_gene_sets, function(genes) {
    unique(as.character(genes))
  })
  full_discovered_genes <- unique(as.character(full_discovered_genes))
  rows <- lapply(seq_along(seeds), function(k) {
    subsets <- utils::combn(seq_along(seeds), k, simplify = FALSE)
    union_counts <- vapply(subsets, function(indices) {
      length(unique(unlist(discovered_gene_sets[indices], use.names = FALSE)))
    }, integer(1))
    recovered_counts <- vapply(subsets, function(indices) {
      genes <- unique(unlist(discovered_gene_sets[indices], use.names = FALSE))
      length(intersect(genes, full_discovered_genes))
    }, integer(1))
    data.frame(
      n_seeds = k,
      n_subsets = length(subsets),
      unique_min = min(union_counts),
      unique_median = stats::median(union_counts),
      unique_max = max(union_counts),
      recovered_min = min(recovered_counts),
      recovered_median = stats::median(recovered_counts),
      recovered_max = max(recovered_counts),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
