# Helpers for pairwise-difference correlation estimates across lfdr thresholds.

select_gene_representatives_above_lfdr <- function(pair_keys,
                                                   lfdr,
                                                   threshold,
                                                   minimum_selected = 20L) {
  pair_keys <- as.character(pair_keys)
  lfdr <- as.numeric(lfdr)
  threshold <- as.numeric(threshold)
  minimum_selected <- as.integer(minimum_selected)
  if (length(pair_keys) != length(lfdr) || length(pair_keys) == 0L ||
      any(!nzchar(pair_keys)) || any(!is.finite(lfdr)) ||
      any(lfdr < 0 | lfdr > 1) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold <= 0 || threshold >= 1 ||
      is.na(minimum_selected) || minimum_selected < 2L) {
    stop("Invalid pair keys, lfdr values, threshold, or minimum selection size.")
  }

  gene_id <- parse_gene_ids(pair_keys)
  within_gene_order <- order(
    gene_id,
    -lfdr,
    pair_keys,
    method = "radix"
  )
  representative_indices <- within_gene_order[
    !duplicated(gene_id[within_gene_order])
  ]
  representative_order <- order(
    -lfdr[representative_indices],
    pair_keys[representative_indices],
    method = "radix"
  )
  representative_indices <- representative_indices[representative_order]
  selected_indices <- representative_indices[
    lfdr[representative_indices] > threshold
  ]
  if (length(selected_indices) < minimum_selected) {
    stop(
      "Only ",
      length(selected_indices),
      " gene representatives satisfy lfdr > ",
      threshold,
      "."
    )
  }

  selected <- data.frame(
    rank = seq_along(selected_indices),
    fash_index = selected_indices,
    pair_key = pair_keys[selected_indices],
    gene_id = gene_id[selected_indices],
    variant_id = sub("^[^_]+_", "", pair_keys[selected_indices]),
    lfdr = lfdr[selected_indices],
    threshold = threshold,
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(selected$gene_id) || anyDuplicated(selected$pair_key) ||
      any(selected$lfdr <= threshold) || any(diff(selected$lfdr) > 0)) {
    stop("The lfdr-threshold selection failed its invariants.")
  }

  list(
    selected = selected,
    selected_indices = selected_indices,
    representative_indices = representative_indices,
    n_gene_representatives = length(representative_indices),
    threshold = threshold
  )
}

bootstrap_pairwise_lag_variogram <- function(beta_hat,
                                             se,
                                             n_bootstrap = 1000L,
                                             seed = 20260821L) {
  beta_hat <- validate_finite_matrix(beta_hat, "beta_hat")
  se <- validate_finite_matrix(se, "se", positive = TRUE)
  n_bootstrap <- as.integer(n_bootstrap)
  seed <- as.integer(seed)
  if (!identical(dim(beta_hat), dim(se)) ||
      is.na(n_bootstrap) || n_bootstrap < 20L ||
      is.na(seed) || nrow(beta_hat) < 20L || ncol(beta_hat) < 2L) {
    stop("Invalid inputs for the pairwise gene bootstrap.")
  }

  observed_matrix <- estimate_pairwise_difference_correlation(beta_hat, se)
  observed_correlation <- lag_average_correlation(observed_matrix)
  n_units <- nrow(beta_hat)
  draws <- matrix(
    NA_real_,
    nrow = n_bootstrap,
    ncol = length(observed_correlation)
  )
  set.seed(seed)
  for (replication in seq_len(n_bootstrap)) {
    index <- sample.int(n_units, size = n_units, replace = TRUE)
    draws[replication, ] <- lag_average_correlation(
      estimate_pairwise_difference_correlation(
        beta_hat[index, , drop = FALSE],
        se[index, , drop = FALSE]
      )
    )
  }
  semivariogram_draws <- 1 - draws

  summary <- data.frame(
    lag = seq_along(observed_correlation),
    observed_correlation = observed_correlation,
    observed_semivariogram = 1 - observed_correlation,
    bootstrap_mean_correlation = colMeans(draws),
    correlation_ci_lower = apply(
      draws,
      2L,
      stats::quantile,
      probs = 0.025,
      names = FALSE
    ),
    correlation_ci_upper = apply(
      draws,
      2L,
      stats::quantile,
      probs = 0.975,
      names = FALSE
    ),
    bootstrap_mean_semivariogram = colMeans(semivariogram_draws),
    semivariogram_ci_lower = apply(
      semivariogram_draws,
      2L,
      stats::quantile,
      probs = 0.025,
      names = FALSE
    ),
    semivariogram_ci_upper = apply(
      semivariogram_draws,
      2L,
      stats::quantile,
      probs = 0.975,
      names = FALSE
    ),
    n_bootstrap = n_bootstrap,
    stringsAsFactors = FALSE
  )

  list(
    observed_matrix = observed_matrix,
    draws = draws,
    summary = summary,
    seed = seed
  )
}
