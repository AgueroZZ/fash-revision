# Helpers for observationally null-enriched residual-covariance pilots.

discover_by_cumulative_lfdr <- function(lfdr, alpha = 0.05) {
  lfdr <- as.numeric(lfdr)
  alpha <- as.numeric(alpha)
  if (length(lfdr) == 0L || any(!is.finite(lfdr)) ||
      any(lfdr < 0 | lfdr > 1) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Invalid lfdr values or alpha.")
  }

  ordering <- order(lfdr, method = "radix")
  cumulative_fdr <- cumsum(lfdr[ordering]) / seq_along(ordering)
  passing_ranks <- which(cumulative_fdr <= alpha)
  discovered <- rep(FALSE, length(lfdr))
  cutoff_rank <- 0L
  cutoff_lfdr <- NA_real_
  if (length(passing_ranks) > 0L) {
    cutoff_rank <- max(passing_ranks)
    discovered[ordering[seq_len(cutoff_rank)]] <- TRUE
    cutoff_lfdr <- lfdr[ordering[cutoff_rank]]
  }

  list(
    discovered = discovered,
    ordering = ordering,
    cumulative_fdr = cumulative_fdr,
    cutoff_rank = cutoff_rank,
    cutoff_lfdr = cutoff_lfdr,
    alpha = alpha
  )
}

make_randomized_pair_catalog <- function(pair_keys, seed = 20260806L) {
  pair_keys <- as.character(pair_keys)
  seed <- as.integer(seed)
  if (length(pair_keys) == 0L || any(!nzchar(pair_keys)) ||
      anyDuplicated(pair_keys) || any(!grepl("_", pair_keys, fixed = TRUE)) ||
      is.na(seed)) {
    stop("Invalid pair keys or random seed.")
  }

  gene_id <- sub("_.*$", "", pair_keys)
  variant_id <- sub("^[^_]+_", "", pair_keys)
  if (any(!nzchar(gene_id)) || any(!nzchar(variant_id))) {
    stop("At least one pair key could not be parsed.")
  }

  pair_order <- order(pair_keys, method = "radix")
  pair_score <- numeric(length(pair_keys))
  set.seed(seed)
  pair_score[pair_order] <- stats::runif(length(pair_keys))

  sorted_genes <- sort(unique(gene_id), method = "radix")
  gene_score_lookup <- stats::setNames(
    stats::runif(length(sorted_genes)),
    sorted_genes
  )
  variant_gene_count <- table(variant_id)

  data.frame(
    fash_index = seq_along(pair_keys),
    pair_key = pair_keys,
    gene_id = gene_id,
    variant_id = variant_id,
    variant_gene_count = as.integer(variant_gene_count[variant_id]),
    pair_random_score = pair_score,
    gene_random_score = unname(gene_score_lookup[gene_id]),
    stringsAsFactors = FALSE
  )
}

select_random_unique_variant_per_gene <- function(catalog,
                                                  eligible_genes,
                                                  n_select = 200L,
                                                  set_id = "null_set") {
  required_columns <- c(
    "fash_index", "pair_key", "gene_id", "variant_id",
    "variant_gene_count", "pair_random_score", "gene_random_score"
  )
  eligible_genes <- unique(as.character(eligible_genes))
  n_select <- as.integer(n_select)
  if (!is.data.frame(catalog) ||
      !all(required_columns %in% names(catalog)) ||
      length(eligible_genes) == 0L || any(!nzchar(eligible_genes)) ||
      is.na(n_select) || n_select < 2L || !nzchar(set_id)) {
    stop("Invalid catalog, eligible genes, selection size, or set ID.")
  }

  candidates <- catalog[
    catalog$gene_id %in% eligible_genes & catalog$variant_gene_count == 1L,
    required_columns,
    drop = FALSE
  ]
  if (nrow(candidates) == 0L) {
    stop("No globally gene-unique variants are available for eligible genes.")
  }

  candidate_order <- order(
    candidates$gene_id,
    candidates$pair_random_score,
    candidates$pair_key,
    method = "radix"
  )
  representatives <- candidates[candidate_order, , drop = FALSE]
  representatives <- representatives[
    !duplicated(representatives$gene_id),
    ,
    drop = FALSE
  ]
  representative_order <- order(
    representatives$gene_random_score,
    representatives$gene_id,
    method = "radix"
  )
  representatives <- representatives[representative_order, , drop = FALSE]
  realized_n <- min(n_select, nrow(representatives))
  selected <- representatives[seq_len(realized_n), , drop = FALSE]
  selected$set_id <- set_id
  selected$selection_rank <- seq_len(realized_n)
  selected <- selected[, c(
    "set_id", "selection_rank", required_columns
  )]
  rownames(selected) <- NULL

  if (nrow(selected) < 2L || anyDuplicated(selected$gene_id) ||
      anyDuplicated(selected$variant_id) ||
      any(selected$variant_gene_count != 1L)) {
    stop("The one-pair-per-gene selection failed its invariants.")
  }

  list(
    selected = selected,
    selected_indices = selected$fash_index,
    n_requested = n_select,
    n_selected = nrow(selected),
    n_eligible_genes = length(eligible_genes),
    n_eligible_genes_with_unique_variant = nrow(representatives),
    n_candidate_pairs = nrow(candidates),
    set_id = set_id
  )
}

estimate_direct_zero_mean_covariance <- function(beta_hat, se) {
  beta_hat <- validate_finite_matrix(beta_hat, "beta_hat")
  se <- validate_finite_matrix(se, "se", positive = TRUE)
  if (!identical(dim(beta_hat), dim(se))) {
    stop("beta_hat and se must have the same dimensions.")
  }
  z_score <- beta_hat / se
  covariance <- crossprod(z_score) / nrow(z_score)
  covariance <- (covariance + t(covariance)) / 2
  dimnames(covariance) <- list(colnames(beta_hat), colnames(beta_hat))
  covariance
}

estimate_within_unit_centered_matrices <- function(beta_hat, se) {
  beta_hat <- validate_finite_matrix(beta_hat, "beta_hat")
  se <- validate_finite_matrix(se, "se", positive = TRUE)
  if (!identical(dim(beta_hat), dim(se))) {
    stop("beta_hat and se must have the same dimensions.")
  }

  centered <- weighted_center_standardize(beta_hat, se)
  standardized_residual <- centered$standardized_residual
  covariance <- crossprod(standardized_residual) /
    nrow(standardized_residual)
  covariance <- (covariance + t(covariance)) / 2
  correlation <- stats::cor(standardized_residual)
  correlation <- (correlation + t(correlation)) / 2
  dimnames(covariance) <- dimnames(correlation) <- list(
    colnames(beta_hat),
    colnames(beta_hat)
  )
  if (any(!is.finite(covariance)) || any(!is.finite(correlation))) {
    stop("The within-unit centered matrices contain non-finite values.")
  }

  list(
    covariance = covariance,
    correlation = correlation,
    standardized_residual = standardized_residual,
    constant_estimate = centered$constant_estimate
  )
}

sample_lfdr_quartile_gallery <- function(metadata,
                                         n_per_quartile = 25L,
                                         seed = 20260891L) {
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  n_per_quartile <- as.integer(n_per_quartile)
  seed <- as.integer(seed)
  required_columns <- c("pair_key", "bf_adjusted_lfdr")
  if (!all(required_columns %in% names(metadata)) ||
      nrow(metadata) < 4L || any(!nzchar(metadata$pair_key)) ||
      anyDuplicated(metadata$pair_key) ||
      any(!is.finite(metadata$bf_adjusted_lfdr)) ||
      any(metadata$bf_adjusted_lfdr < 0 |
            metadata$bf_adjusted_lfdr > 1) ||
      is.na(n_per_quartile) || n_per_quartile < 1L || is.na(seed)) {
    stop("Invalid metadata, gallery size, or seed.")
  }

  ordering <- order(
    -metadata$bf_adjusted_lfdr,
    metadata$pair_key,
    method = "radix"
  )
  ranked <- metadata[ordering, , drop = FALSE]
  ranked$set_a_lfdr_rank <- seq_len(nrow(ranked))
  ranked$lfdr_quartile <- pmin(
    4L,
    ceiling(4 * ranked$set_a_lfdr_rank / nrow(ranked))
  )
  quartile_sizes <- tabulate(ranked$lfdr_quartile, nbins = 4L)
  if (any(quartile_sizes < n_per_quartile)) {
    stop("Every lfdr quartile must contain at least n_per_quartile units.")
  }

  set.seed(seed)
  sampled_rows <- unlist(lapply(seq_len(4L), function(quartile) {
    candidates <- which(ranked$lfdr_quartile == quartile)
    sample(candidates, n_per_quartile, replace = FALSE)
  }))
  selected <- ranked[sampled_rows, , drop = FALSE]
  selected <- selected[order(
    selected$lfdr_quartile,
    selected$set_a_lfdr_rank,
    method = "radix"
  ), , drop = FALSE]
  selected$gallery_page_position <- ave(
    selected$set_a_lfdr_rank,
    selected$lfdr_quartile,
    FUN = seq_along
  )
  rownames(selected) <- NULL
  if (nrow(selected) != 4L * n_per_quartile ||
      anyDuplicated(selected$pair_key) ||
      !identical(
        as.integer(table(selected$lfdr_quartile)),
        rep(n_per_quartile, 4L)
      )) {
    stop("The stratified gallery sample failed its invariants.")
  }
  selected
}

lag_covariance_variogram <- function(covariance, se = NULL) {
  covariance <- validate_finite_matrix(covariance, "covariance")
  if (nrow(covariance) != ncol(covariance)) {
    stop("covariance must be square.")
  }
  if (!is.null(se)) {
    se <- validate_finite_matrix(se, "se", positive = TRUE)
    if (ncol(se) != ncol(covariance)) {
      stop("se and covariance use different time dimensions.")
    }
  }

  n_time <- ncol(covariance)
  output <- lapply(seq_len(n_time - 1L), function(lag) {
    time_a <- seq_len(n_time - lag)
    time_b <- time_a + lag
    covariance_values <- covariance[cbind(time_a, time_b)]
    variogram_values <- 0.5 * (
      diag(covariance)[time_a] + diag(covariance)[time_b] -
        2 * covariance_values
    )
    result <- data.frame(
      lag = lag,
      mean_standardized_covariance = mean(covariance_values),
      mean_semivariogram = mean(variogram_values),
      stringsAsFactors = FALSE
    )

    if (!is.null(se)) {
      beta_covariance <- vapply(seq_along(time_a), function(index) {
        stats::median(
          se[, time_a[index]] * se[, time_b[index]] *
            covariance[time_a[index], time_b[index]]
        )
      }, numeric(1))
      beta_semivariance <- vapply(seq_along(time_a), function(index) {
        stats::median(0.5 * (
          se[, time_a[index]]^2 * covariance[time_a[index], time_a[index]] +
            se[, time_b[index]]^2 * covariance[time_b[index], time_b[index]] -
            2 * se[, time_a[index]] * se[, time_b[index]] *
              covariance[time_a[index], time_b[index]]
        ))
      }, numeric(1))
      result$mean_median_beta_scale_covariance <- mean(beta_covariance)
      result$mean_median_beta_scale_semivariance <- mean(beta_semivariance)
    }
    result
  })
  do.call(rbind, output)
}

scale_covariance_by_median_se <- function(covariance, se) {
  covariance <- validate_finite_matrix(covariance, "covariance")
  se <- validate_finite_matrix(se, "se", positive = TRUE)
  if (nrow(covariance) != ncol(covariance) ||
      ncol(se) != ncol(covariance)) {
    stop("Invalid covariance or se dimensions.")
  }
  n_time <- ncol(covariance)
  scaled <- matrix(NA_real_, nrow = n_time, ncol = n_time)
  for (time_a in seq_len(n_time)) {
    for (time_b in seq_len(n_time)) {
      scaled[time_a, time_b] <- covariance[time_a, time_b] *
        stats::median(se[, time_a] * se[, time_b])
    }
  }
  dimnames(scaled) <- dimnames(covariance)
  (scaled + t(scaled)) / 2
}

covariance_matrix_diagnostics <- function(covariance,
                                          set_id,
                                          estimator) {
  covariance <- validate_finite_matrix(covariance, "covariance")
  if (nrow(covariance) != ncol(covariance)) {
    stop("covariance must be square.")
  }
  eigenvalues <- eigen(
    (covariance + t(covariance)) / 2,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  off_diagonal <- covariance[upper.tri(covariance)]
  data.frame(
    set_id = set_id,
    estimator = estimator,
    minimum_eigenvalue = min(eigenvalues),
    maximum_eigenvalue = max(eigenvalues),
    n_negative_eigenvalues = sum(eigenvalues < -1e-10),
    minimum_diagonal = min(diag(covariance)),
    maximum_diagonal = max(diag(covariance)),
    mean_off_diagonal = mean(off_diagonal),
    minimum_off_diagonal = min(off_diagonal),
    maximum_off_diagonal = max(off_diagonal),
    stringsAsFactors = FALSE
  )
}

covariance_matrix_to_long <- function(covariance,
                                      set_id,
                                      estimator,
                                      scale) {
  covariance <- validate_finite_matrix(covariance, "covariance")
  if (nrow(covariance) != ncol(covariance) || !nzchar(scale)) {
    stop("Invalid covariance matrix or scale label.")
  }
  time_labels <- colnames(covariance)
  if (is.null(time_labels)) {
    time_labels <- as.character(seq_len(ncol(covariance)) - 1L)
  } else {
    time_labels <- sub("^time_", "", time_labels)
  }
  grid <- expand.grid(
    time_a = as.numeric(time_labels),
    time_b = as.numeric(time_labels),
    stringsAsFactors = FALSE
  )
  grid$covariance <- as.vector(covariance)
  grid$set_id <- set_id
  grid$estimator <- estimator
  grid$scale <- scale
  grid[, c(
    "set_id", "estimator", "scale", "time_a", "time_b", "covariance"
  )]
}

bootstrap_covariance_lags <- function(beta_hat,
                                      se,
                                      estimators,
                                      n_bootstrap = 200L,
                                      seed = 20260816L) {
  beta_hat <- validate_finite_matrix(beta_hat, "beta_hat")
  se <- validate_finite_matrix(se, "se", positive = TRUE)
  n_bootstrap <- as.integer(n_bootstrap)
  seed <- as.integer(seed)
  valid_estimators <- c(
    "Pairwise difference",
    "Direct zero mean",
    "Within-unit centered covariance"
  )
  estimators <- unique(as.character(estimators))
  if (!identical(dim(beta_hat), dim(se)) ||
      length(estimators) == 0L || !all(estimators %in% valid_estimators) ||
      is.na(n_bootstrap) || n_bootstrap < 20L || is.na(seed)) {
    stop("Invalid bootstrap inputs.")
  }

  estimate_one <- function(estimator, index) {
    if (identical(estimator, "Pairwise difference")) {
      estimate_pairwise_difference_correlation(
        beta_hat[index, , drop = FALSE],
        se[index, , drop = FALSE]
      )
    } else if (identical(estimator, "Direct zero mean")) {
      estimate_direct_zero_mean_covariance(
        beta_hat[index, , drop = FALSE],
        se[index, , drop = FALSE]
      )
    } else {
      estimate_within_unit_centered_matrices(
        beta_hat[index, , drop = FALSE],
        se[index, , drop = FALSE]
      )$covariance
    }
  }

  all_indices <- seq_len(nrow(beta_hat))
  observed_matrices <- lapply(estimators, estimate_one, index = all_indices)
  names(observed_matrices) <- estimators
  observed_lags <- lapply(observed_matrices, lag_covariance_variogram, se = se)
  n_lag <- ncol(beta_hat) - 1L
  draws <- lapply(estimators, function(unused) {
    list(
      covariance = matrix(NA_real_, nrow = n_bootstrap, ncol = n_lag),
      semivariogram = matrix(NA_real_, nrow = n_bootstrap, ncol = n_lag)
    )
  })
  names(draws) <- estimators

  set.seed(seed)
  for (replication in seq_len(n_bootstrap)) {
    index <- sample.int(nrow(beta_hat), nrow(beta_hat), replace = TRUE)
    for (estimator in estimators) {
      matrix_estimate <- estimate_one(estimator, index)
      lag_estimate <- lag_covariance_variogram(
        matrix_estimate,
        se[index, , drop = FALSE]
      )
      draws[[estimator]]$covariance[replication, ] <-
        lag_estimate$mean_standardized_covariance
      draws[[estimator]]$semivariogram[replication, ] <-
        lag_estimate$mean_semivariogram
    }
  }

  summary_rows <- lapply(estimators, function(estimator) {
    covariance_draws <- draws[[estimator]]$covariance
    variogram_draws <- draws[[estimator]]$semivariogram
    observed <- observed_lags[[estimator]]
    data.frame(
      estimator = estimator,
      lag = seq_len(n_lag),
      observed_covariance = observed$mean_standardized_covariance,
      covariance_bootstrap_mean = colMeans(covariance_draws),
      covariance_ci_lower = apply(
        covariance_draws, 2L, stats::quantile,
        probs = 0.025, names = FALSE
      ),
      covariance_ci_upper = apply(
        covariance_draws, 2L, stats::quantile,
        probs = 0.975, names = FALSE
      ),
      observed_semivariogram = observed$mean_semivariogram,
      semivariogram_bootstrap_mean = colMeans(variogram_draws),
      semivariogram_ci_lower = apply(
        variogram_draws, 2L, stats::quantile,
        probs = 0.025, names = FALSE
      ),
      semivariogram_ci_upper = apply(
        variogram_draws, 2L, stats::quantile,
        probs = 0.975, names = FALSE
      ),
      n_bootstrap = n_bootstrap,
      bootstrap_seed = seed,
      stringsAsFactors = FALSE
    )
  })

  list(
    observed_matrices = observed_matrices,
    observed_lags = observed_lags,
    draws = draws,
    summary = do.call(rbind, summary_rows),
    n_bootstrap = n_bootstrap,
    seed = seed
  )
}
