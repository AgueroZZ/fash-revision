# Helpers for one-variant-per-gene FASH refit sensitivity analyses.

validate_probability_vector <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) < 2L || any(!is.finite(x)) || any(x < 0 | x > 1)) {
    stop(name, " must contain at least two finite probabilities.")
  }
  x
}

select_minimum_lfdr_variant_per_gene <- function(pair_keys,
                                                 lfdr,
                                                 gene_index = NULL) {
  pair_keys <- as.character(pair_keys)
  lfdr <- validate_probability_vector(lfdr, "lfdr")
  if (length(pair_keys) != length(lfdr) || any(!nzchar(pair_keys)) ||
      anyDuplicated(pair_keys)) {
    stop("pair_keys and lfdr must be aligned unique full-data vectors.")
  }

  gene_id <- sub("_.*$", "", pair_keys)
  if (is.null(gene_index)) {
    gene_index <- split(
      seq_along(pair_keys),
      factor(gene_id, levels = sort(unique(gene_id)))
    )
  }
  flattened_indices <- unlist(gene_index, use.names = FALSE)
  if (!is.list(gene_index) || length(gene_index) < 2L ||
      any(lengths(gene_index) < 1L) ||
      length(flattened_indices) != length(pair_keys) ||
      anyNA(flattened_indices) ||
      any(flattened_indices < 1L | flattened_indices > length(pair_keys)) ||
      anyDuplicated(flattened_indices) ||
      !identical(sort(as.integer(flattened_indices)), seq_along(pair_keys)) ||
      any(vapply(gene_index, function(indices) {
        length(unique(gene_id[indices])) != 1L
      }, logical(1)))) {
    stop("gene_index must partition all pair keys into non-empty gene groups.")
  }

  selected_indices <- vapply(gene_index, function(indices) {
    minimum_lfdr <- min(lfdr[indices])
    candidates <- indices[lfdr[indices] == minimum_lfdr]
    candidates[order(pair_keys[candidates], method = "radix")][1]
  }, integer(1))
  selected_gene <- gene_id[selected_indices]
  if (length(selected_indices) != length(gene_index) ||
      anyDuplicated(selected_indices) || anyDuplicated(selected_gene)) {
    stop("Minimum-lfdr selection did not produce exactly one pair per gene.")
  }

  data.frame(
    selection_rule = "minimum_full_bf_lfdr",
    fash_index = selected_indices,
    pair_key = pair_keys[selected_indices],
    gene_id = selected_gene,
    variant_id = sub("^[^_]+_", "", pair_keys[selected_indices]),
    selection_lfdr = lfdr[selected_indices],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

make_fash_refit_datasets <- function(fash_fit, selected_indices) {
  if (is.null(fash_fit$fash_data$data_list) || is.null(fash_fit$fash_data$S)) {
    stop("The FASH object must contain fash_data$data_list and fash_data$S.")
  }
  data_list <- fash_fit$fash_data$data_list
  se_list <- fash_fit$fash_data$S
  if (length(data_list) != length(se_list) || is.null(names(data_list)) ||
      any(!nzchar(names(data_list))) || anyDuplicated(names(data_list))) {
    stop("The FASH data and SE lists must have matching unique pair keys.")
  }
  selected_indices <- as.integer(selected_indices)
  if (length(selected_indices) < 2L || anyNA(selected_indices) ||
      any(selected_indices < 1L) || any(selected_indices > length(data_list)) ||
      anyDuplicated(selected_indices)) {
    stop("selected_indices must contain at least two unique in-range indices.")
  }

  selected_data <- data_list[selected_indices]
  selected_se <- se_list[selected_indices]
  expected_time <- NULL
  refit_datasets <- Map(function(dataset, se, pair_key) {
    if (!is.data.frame(dataset) || !all(c("y", "x") %in% names(dataset))) {
      stop("Every selected FASH dataset must contain y and x columns.")
    }
    beta <- as.numeric(dataset$y)
    time <- as.numeric(dataset$x)
    se <- as.numeric(se)
    if (length(beta) < 2L || length(time) != length(beta) ||
        length(se) != length(beta) || any(!is.finite(beta)) ||
        any(!is.finite(time)) || any(!is.finite(se)) || any(se <= 0) ||
        anyDuplicated(time)) {
      stop("Invalid beta, time, or SE values for selected pair ", pair_key, ".")
    }
    if (is.null(expected_time)) {
      expected_time <<- time
    } else if (!identical(time, expected_time)) {
      stop("All selected FASH datasets must use an identical time grid.")
    }
    data.frame(beta = beta, time = time, SE = se)
  }, selected_data, selected_se, names(selected_data))
  names(refit_datasets) <- names(selected_data)
  refit_datasets
}

refit_fash_from_cached_likelihood <- function(
    fash_fit,
    selected_indices,
    penalty = fash_fit$settings$penalty) {
  if (!requireNamespace("fashr", quietly = TRUE)) {
    stop("The fashr package is required for cached-likelihood refitting.")
  }
  required_fields <- c(
    "prior_weights", "posterior_weights", "psd_grid", "lfdr",
    "settings", "fash_data", "L_matrix"
  )
  required_settings <- c(
    "num_basis", "order", "betaprec", "pred_step", "likelihood", "penalty"
  )
  if (!all(required_fields %in% names(fash_fit)) ||
      !all(required_settings %in% names(fash_fit$settings))) {
    stop("The source FASH object or its settings are incomplete.")
  }
  penalty <- as.numeric(penalty)
  if (length(penalty) != 1L || !is.finite(penalty) || penalty < 1 ||
      abs(penalty - round(penalty)) > sqrt(.Machine$double.eps)) {
    stop("penalty must be one finite integer greater than or equal to one.")
  }
  selected_indices <- as.integer(selected_indices)
  refit_datasets <- make_fash_refit_datasets(fash_fit, selected_indices)
  pair_keys <- names(refit_datasets)
  selected_likelihood <- fash_fit$L_matrix[selected_indices, , drop = FALSE]
  if (nrow(selected_likelihood) != length(selected_indices) ||
      ncol(selected_likelihood) != length(fash_fit$psd_grid) ||
      anyNA(selected_likelihood) || any(is.nan(selected_likelihood)) ||
      any(selected_likelihood == Inf)) {
    stop("The selected cached likelihood matrix is invalid.")
  }
  rownames(selected_likelihood) <- pair_keys

  eb_result <- fashr::fash_eb_est(
    L_matrix = selected_likelihood,
    grid = fash_fit$psd_grid,
    penalty = penalty
  )
  rownames(eb_result$posterior_weight) <- pair_keys
  null_column <- which(eb_result$prior_weight$psd == 0)
  lfdr <- if (length(null_column) == 1L) {
    eb_result$posterior_weight[, null_column]
  } else {
    rep(0, length(pair_keys))
  }
  names(lfdr) <- pair_keys

  source_data <- fash_fit$fash_data
  selected_omega <- source_data$Omega
  if (!is.null(selected_omega) && is.list(selected_omega) &&
      length(selected_omega) == length(source_data$data_list)) {
    selected_omega <- selected_omega[selected_indices]
  }
  selected_fash_data <- list(
    data_list = source_data$data_list[selected_indices],
    S = source_data$S[selected_indices],
    Omega = selected_omega
  )
  refit_settings <- fash_fit$settings
  refit_settings$penalty <- penalty
  structure(
    list(
      prior_weights = eb_result$prior_weight,
      posterior_weights = eb_result$posterior_weight,
      psd_grid = fash_fit$psd_grid,
      lfdr = lfdr,
      settings = refit_settings,
      fash_data = selected_fash_data,
      L_matrix = selected_likelihood,
      eb_result = eb_result
    ),
    class = "fash"
  )
}

validate_prior_weights <- function(prior_weights, name) {
  prior_weights <- as.data.frame(prior_weights, stringsAsFactors = FALSE)
  if (!all(c("psd", "prior_weight") %in% names(prior_weights))) {
    stop(name, " must contain psd and prior_weight columns.")
  }
  prior_weights <- prior_weights[, c("psd", "prior_weight")]
  prior_weights$psd <- as.numeric(prior_weights$psd)
  prior_weights$prior_weight <- as.numeric(prior_weights$prior_weight)
  if (nrow(prior_weights) < 1L || any(!is.finite(as.matrix(prior_weights))) ||
      any(prior_weights$psd < 0) || any(prior_weights$prior_weight < 0) ||
      anyDuplicated(prior_weights$psd) ||
      abs(sum(prior_weights$prior_weight) - 1) > 1e-6 ||
      !any(prior_weights$psd == 0)) {
    stop(name, " is not a valid normalized FASH prior-weight table.")
  }
  prior_weights[order(prior_weights$psd), , drop = FALSE]
}

compare_prior_weights <- function(full_prior, thinned_prior, fit_stage) {
  full_prior <- validate_prior_weights(full_prior, "full_prior")
  thinned_prior <- validate_prior_weights(thinned_prior, "thinned_prior")
  fit_stage <- as.character(fit_stage)
  if (length(fit_stage) != 1L || is.na(fit_stage) || !nzchar(fit_stage)) {
    stop("fit_stage must be one non-empty label.")
  }
  support <- sort(unique(c(full_prior$psd, thinned_prior$psd)))
  full_weight <- full_prior$prior_weight[match(support, full_prior$psd)]
  thinned_weight <- thinned_prior$prior_weight[
    match(support, thinned_prior$psd)
  ]
  full_weight[is.na(full_weight)] <- 0
  thinned_weight[is.na(thinned_weight)] <- 0
  comparison <- data.frame(
    fit_stage = fit_stage,
    psd = support,
    full_weight = full_weight,
    thinned_weight = thinned_weight,
    absolute_difference = abs(full_weight - thinned_weight),
    stringsAsFactors = FALSE
  )
  null_row <- which(comparison$psd == 0)
  summary <- data.frame(
    fit_stage = fit_stage,
    full_pi0 = comparison$full_weight[null_row],
    thinned_pi0 = comparison$thinned_weight[null_row],
    pi0_difference = comparison$thinned_weight[null_row] -
      comparison$full_weight[null_row],
    prior_total_variation = 0.5 * sum(comparison$absolute_difference),
    stringsAsFactors = FALSE
  )
  list(table = comparison, summary = summary)
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

safe_correlation <- function(x, y, method) {
  if (stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }
  unname(stats::cor(x, y, method = method))
}

compare_paired_lfdr <- function(full_lfdr,
                                thinned_lfdr,
                                pair_keys,
                                fit_stage,
                                alpha = 0.05) {
  full_lfdr <- validate_probability_vector(full_lfdr, "full_lfdr")
  thinned_lfdr <- validate_probability_vector(thinned_lfdr, "thinned_lfdr")
  pair_keys <- as.character(pair_keys)
  fit_stage <- as.character(fit_stage)
  if (length(full_lfdr) != length(thinned_lfdr) ||
      length(pair_keys) != length(full_lfdr) || any(!nzchar(pair_keys)) ||
      anyDuplicated(pair_keys) || length(fit_stage) != 1L ||
      is.na(fit_stage) || !nzchar(fit_stage)) {
    stop("Paired lfdr inputs, pair keys, or fit_stage are invalid.")
  }
  full_calls <- cumulative_fdr_calls(full_lfdr, alpha = alpha)
  thinned_calls <- cumulative_fdr_calls(thinned_lfdr, alpha = alpha)
  call_union <- union(full_calls, thinned_calls)
  call_intersection <- intersect(full_calls, thinned_calls)
  jaccard <- if (length(call_union) == 0L) {
    1
  } else {
    length(call_intersection) / length(call_union)
  }
  difference <- thinned_lfdr - full_lfdr
  comparison <- data.frame(
    fit_stage = fit_stage,
    pair_key = pair_keys,
    full_lfdr = full_lfdr,
    thinned_lfdr = thinned_lfdr,
    lfdr_difference = difference,
    absolute_difference = abs(difference),
    full_fdr_call = seq_along(full_lfdr) %in% full_calls,
    thinned_fdr_call = seq_along(thinned_lfdr) %in% thinned_calls,
    stringsAsFactors = FALSE
  )
  summary <- data.frame(
    fit_stage = fit_stage,
    n_units = length(full_lfdr),
    pearson_lfdr = safe_correlation(full_lfdr, thinned_lfdr, "pearson"),
    spearman_lfdr = safe_correlation(full_lfdr, thinned_lfdr, "spearman"),
    mean_absolute_lfdr_difference = mean(abs(difference)),
    median_absolute_lfdr_difference = stats::median(abs(difference)),
    rmse_lfdr = sqrt(mean(difference^2)),
    full_mean_lfdr = mean(full_lfdr),
    thinned_mean_lfdr = mean(thinned_lfdr),
    full_fdr_calls = length(full_calls),
    thinned_fdr_calls = length(thinned_calls),
    fdr_call_intersection = length(call_intersection),
    fdr_call_union = length(call_union),
    fdr_call_jaccard = jaccard,
    alpha = alpha,
    stringsAsFactors = FALSE
  )
  list(table = comparison, summary = summary)
}
