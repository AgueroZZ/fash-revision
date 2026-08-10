# Helpers for estimating full time-correlation patterns from null-like eQTLs.

validate_finite_matrix <- function(x,
                                   name,
                                   minimum_rows = 2L,
                                   minimum_columns = 2L,
                                   positive = FALSE) {
  x <- as.matrix(x)
  if (nrow(x) < minimum_rows || ncol(x) < minimum_columns ||
      any(!is.finite(x)) || (positive && any(x <= 0))) {
    stop(name, " must be a finite matrix with the required dimensions",
         if (positive) " and positive entries." else ".")
  }
  x
}

parse_gene_ids <- function(pair_keys) {
  pair_keys <- as.character(pair_keys)
  if (length(pair_keys) == 0L || any(!nzchar(pair_keys)) ||
      any(!grepl("_", pair_keys, fixed = TRUE))) {
    stop("Every FASH dataset key must contain a gene and variant separated by an underscore.")
  }
  sub("_.*$", "", pair_keys)
}

select_highest_lfdr_per_gene <- function(pair_keys, lfdr, top_n = 500L) {
  pair_keys <- as.character(pair_keys)
  lfdr <- as.numeric(lfdr)
  top_n <- as.integer(top_n)
  if (length(pair_keys) != length(lfdr) || length(lfdr) == 0L ||
      any(!is.finite(lfdr)) || any(lfdr < 0 | lfdr > 1) ||
      is.na(top_n) || top_n < 2L) {
    stop("Invalid pair keys, lfdr values, or top_n.")
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
  if (length(representative_indices) < top_n) {
    stop("There are fewer gene representatives than top_n.")
  }
  representative_order <- order(
    -lfdr[representative_indices],
    pair_keys[representative_indices],
    method = "radix"
  )
  representative_indices <- representative_indices[representative_order]
  selected_indices <- representative_indices[seq_len(top_n)]

  selected <- data.frame(
    rank = seq_len(top_n),
    fash_index = selected_indices,
    pair_key = pair_keys[selected_indices],
    gene_id = gene_id[selected_indices],
    variant_id = sub("^[^_]+_", "", pair_keys[selected_indices]),
    lfdr = lfdr[selected_indices],
    stringsAsFactors = FALSE
  )
  if (anyDuplicated(selected$gene_id) || anyDuplicated(selected$pair_key)) {
    stop("The null-like selection is not one-pair-per-gene.")
  }

  list(
    selected = selected,
    selected_indices = selected_indices,
    representative_indices = representative_indices,
    n_genes = length(representative_indices)
  )
}

extract_fash_beta_se <- function(fash_fit, indices) {
  if (is.null(fash_fit$fash_data$data_list) || is.null(fash_fit$fash_data$S)) {
    stop("The fitted FASH object does not contain data_list and S.")
  }
  indices <- as.integer(indices)
  if (length(indices) < 2L || anyNA(indices) || any(indices < 1L) ||
      any(indices > length(fash_fit$fash_data$data_list))) {
    stop("Invalid FASH dataset indices.")
  }
  data_subset <- fash_fit$fash_data$data_list[indices]
  se_subset <- fash_fit$fash_data$S[indices]
  beta_hat <- do.call(rbind, lapply(data_subset, function(x) x$y))
  se <- do.call(rbind, se_subset)
  beta_hat <- validate_finite_matrix(beta_hat, "beta_hat")
  se <- validate_finite_matrix(se, "se", positive = TRUE)
  if (!identical(dim(beta_hat), dim(se))) {
    stop("Extracted beta_hat and se matrices have different dimensions.")
  }

  time_grid <- as.numeric(data_subset[[1]]$x)
  same_grid <- vapply(data_subset, function(x) {
    isTRUE(all.equal(as.numeric(x$x), time_grid))
  }, logical(1))
  if (length(time_grid) != ncol(beta_hat) || any(!same_grid)) {
    stop("The selected FASH datasets do not share one time grid.")
  }
  colnames(beta_hat) <- colnames(se) <- paste0("time_", time_grid)
  rownames(beta_hat) <- rownames(se) <- names(data_subset)

  list(beta_hat = beta_hat, se = se, time_grid = time_grid)
}

weighted_center_standardize <- function(beta_hat, se) {
  beta_hat <- validate_finite_matrix(beta_hat, "beta_hat")
  se <- validate_finite_matrix(se, "se", positive = TRUE)
  if (!identical(dim(beta_hat), dim(se))) {
    stop("beta_hat and se must have the same dimensions.")
  }
  weights <- 1 / se^2
  constant_estimate <- rowSums(beta_hat * weights) / rowSums(weights)
  standardized_residual <- sweep(
    beta_hat,
    1L,
    constant_estimate,
    `-`
  ) / se
  list(
    constant_estimate = constant_estimate,
    standardized_residual = standardized_residual,
    residual_score = rowSums(standardized_residual^2)
  )
}

estimate_direct_centered_correlation <- function(beta_hat, se) {
  centered <- weighted_center_standardize(beta_hat, se)
  correlation <- stats::cor(centered$standardized_residual)
  if (any(!is.finite(correlation))) {
    stop("The direct centered correlation matrix contains non-finite values.")
  }
  correlation
}

estimate_pairwise_difference_correlation <- function(beta_hat, se) {
  beta_hat <- validate_finite_matrix(beta_hat, "beta_hat")
  se <- validate_finite_matrix(se, "se", positive = TRUE)
  if (!identical(dim(beta_hat), dim(se))) {
    stop("beta_hat and se must have the same dimensions.")
  }
  n_time <- ncol(beta_hat)
  correlation <- diag(n_time)
  dimnames(correlation) <- list(colnames(beta_hat), colnames(beta_hat))

  for (time_a in seq_len(n_time - 1L)) {
    for (time_b in (time_a + 1L):n_time) {
      variance_sum <- se[, time_a]^2 + se[, time_b]^2
      correlation_coefficient <-
        2 * se[, time_a] * se[, time_b] / variance_sum
      squared_standardized_difference <-
        (beta_hat[, time_a] - beta_hat[, time_b])^2 / variance_sum
      denominator <- sum(correlation_coefficient^2)
      estimate <- sum(
        correlation_coefficient * (1 - squared_standardized_difference)
      ) / denominator
      correlation[time_a, time_b] <- correlation[time_b, time_a] <- estimate
    }
  }
  correlation
}

invert_t_to_normal_se <- function(beta_hat, adjusted_se, df) {
  beta_hat <- validate_finite_matrix(beta_hat, "beta_hat")
  adjusted_se <- validate_finite_matrix(
    adjusted_se,
    "adjusted_se",
    positive = TRUE
  )
  df <- as.numeric(df)
  if (!identical(dim(beta_hat), dim(adjusted_se)) ||
      length(df) != ncol(beta_hat) || any(!is.finite(df)) || any(df <= 0)) {
    stop("Invalid beta_hat, adjusted_se, or df for SE inversion.")
  }

  raw_se <- adjusted_se
  for (time_index in seq_len(ncol(beta_hat))) {
    z_value <- abs(beta_hat[, time_index] / adjusted_se[, time_index])
    probability <- pmin(stats::pnorm(z_value), 1 - 1e-15)
    t_value <- stats::qt(probability, df = df[time_index])
    ratio <- z_value / t_value
    zero_like <- z_value < sqrt(.Machine$double.eps)
    if (any(zero_like)) {
      density_ratio <- stats::dt(0, df = df[time_index]) / stats::dnorm(0)
      ratio[zero_like] <- density_ratio
    }
    if (any(!is.finite(ratio)) || any(ratio <= 0)) {
      stop("Could not invert the t-to-normal SE correction.")
    }
    raw_se[, time_index] <- adjusted_se[, time_index] * ratio
  }
  raw_se
}

project_to_positive_definite_correlation <- function(correlation,
                                                     eigen_tolerance = 1e-8) {
  correlation <- validate_finite_matrix(correlation, "correlation")
  if (nrow(correlation) != ncol(correlation)) {
    stop("correlation must be square.")
  }
  correlation <- (correlation + t(correlation)) / 2
  diag(correlation) <- 1
  raw_eigenvalues <- eigen(
    correlation,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("The Matrix package is required for nearest-PD projection.")
  }
  projection <- Matrix::nearPD(
    correlation,
    corr = TRUE,
    keepDiag = TRUE,
    eig.tol = eigen_tolerance,
    posd.tol = eigen_tolerance,
    base.matrix = TRUE
  )
  projected <- as.matrix(projection$mat)
  projected <- (projected + t(projected)) / 2
  projected_eigenvalues <- eigen(
    projected,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  difference <- projected - correlation
  list(
    raw = correlation,
    projected = projected,
    diagnostics = data.frame(
      raw_minimum_eigenvalue = min(raw_eigenvalues),
      projected_minimum_eigenvalue = min(projected_eigenvalues),
      maximum_absolute_change = max(abs(difference)),
      frobenius_change = sqrt(sum(difference^2)),
      converged = isTRUE(projection$converged),
      iterations = projection$iterations,
      stringsAsFactors = FALSE
    )
  )
}

lag_average_correlation <- function(correlation) {
  correlation <- validate_finite_matrix(correlation, "correlation")
  if (nrow(correlation) != ncol(correlation)) {
    stop("correlation must be square.")
  }
  n_time <- ncol(correlation)
  vapply(seq_len(n_time - 1L), function(lag) {
    mean(correlation[cbind(
      seq_len(n_time - lag),
      (lag + 1L):n_time
    )])
  }, numeric(1))
}

make_lag_summary <- function(correlation,
                             estimator,
                             matrix_version,
                             se_scale = "t-adjusted") {
  lag_correlation <- lag_average_correlation(correlation)
  data.frame(
    estimator = estimator,
    matrix_version = matrix_version,
    se_scale = se_scale,
    lag = seq_along(lag_correlation),
    mean_correlation = lag_correlation,
    semivariogram = 1 - lag_correlation,
    stringsAsFactors = FALSE
  )
}

summarize_resampled_lags <- function(correlation_array,
                                     estimator,
                                     benchmark = NA_character_) {
  if (length(dim(correlation_array)) != 2L ||
      ncol(correlation_array) < 1L || any(!is.finite(correlation_array))) {
    stop("correlation_array must be a finite replication-by-lag matrix.")
  }
  semivariogram_array <- 1 - correlation_array
  data.frame(
    estimator = estimator,
    benchmark = benchmark,
    lag = seq_len(ncol(correlation_array)),
    mean_correlation = colMeans(correlation_array),
    correlation_ci_lower = apply(
      correlation_array,
      2L,
      stats::quantile,
      probs = 0.025,
      names = FALSE
    ),
    correlation_ci_upper = apply(
      correlation_array,
      2L,
      stats::quantile,
      probs = 0.975,
      names = FALSE
    ),
    mean_semivariogram = colMeans(semivariogram_array),
    semivariogram_ci_lower = apply(
      semivariogram_array,
      2L,
      stats::quantile,
      probs = 0.025,
      names = FALSE
    ),
    semivariogram_ci_upper = apply(
      semivariogram_array,
      2L,
      stats::quantile,
      probs = 0.975,
      names = FALSE
    ),
    n_replications = nrow(correlation_array),
    stringsAsFactors = FALSE
  )
}

bootstrap_correlation_lags <- function(beta_hat,
                                       se,
                                       n_bootstrap = 1000L,
                                       seed = 20260805L) {
  beta_hat <- validate_finite_matrix(beta_hat, "beta_hat")
  se <- validate_finite_matrix(se, "se", positive = TRUE)
  n_bootstrap <- as.integer(n_bootstrap)
  if (!identical(dim(beta_hat), dim(se)) || is.na(n_bootstrap) ||
      n_bootstrap < 20L) {
    stop("Invalid inputs for the gene bootstrap.")
  }
  n_units <- nrow(beta_hat)
  n_lags <- ncol(beta_hat) - 1L
  direct_lags <- matrix(NA_real_, nrow = n_bootstrap, ncol = n_lags)
  pairwise_lags <- matrix(NA_real_, nrow = n_bootstrap, ncol = n_lags)
  set.seed(seed)
  for (replication in seq_len(n_bootstrap)) {
    index <- sample.int(n_units, size = n_units, replace = TRUE)
    direct_lags[replication, ] <- lag_average_correlation(
      estimate_direct_centered_correlation(
        beta_hat[index, , drop = FALSE],
        se[index, , drop = FALSE]
      )
    )
    pairwise_lags[replication, ] <- lag_average_correlation(
      estimate_pairwise_difference_correlation(
        beta_hat[index, , drop = FALSE],
        se[index, , drop = FALSE]
      )
    )
  }
  rbind(
    summarize_resampled_lags(direct_lags, "Direct centered"),
    summarize_resampled_lags(pairwise_lags, "Pairwise difference")
  )
}

run_independence_benchmark <- function(candidate_se,
                                       fixed_selected_se,
                                       top_n = 500L,
                                       n_replications = 1000L,
                                       seed = 20260806L) {
  candidate_se <- validate_finite_matrix(
    candidate_se,
    "candidate_se",
    positive = TRUE
  )
  fixed_selected_se <- validate_finite_matrix(
    fixed_selected_se,
    "fixed_selected_se",
    positive = TRUE
  )
  top_n <- as.integer(top_n)
  n_replications <- as.integer(n_replications)
  if (ncol(candidate_se) != ncol(fixed_selected_se) ||
      nrow(candidate_se) < top_n || nrow(fixed_selected_se) != top_n ||
      is.na(n_replications) || n_replications < 20L) {
    stop("Invalid inputs for the independence benchmark.")
  }

  n_lags <- ncol(candidate_se) - 1L
  benchmark_names <- c(
    "Fixed selected SE patterns",
    "Flatness-selected gene representatives"
  )
  estimator_names <- c("Direct centered", "Pairwise difference")
  arrays <- lapply(benchmark_names, function(unused) {
    lapply(estimator_names, function(unused_estimator) {
      matrix(NA_real_, nrow = n_replications, ncol = n_lags)
    })
  })
  set.seed(seed)

  for (replication in seq_len(n_replications)) {
    fixed_beta <- fixed_selected_se * matrix(
      stats::rnorm(length(fixed_selected_se)),
      nrow = nrow(fixed_selected_se)
    )
    arrays[[1]][[1]][replication, ] <- lag_average_correlation(
      estimate_direct_centered_correlation(fixed_beta, fixed_selected_se)
    )
    arrays[[1]][[2]][replication, ] <- lag_average_correlation(
      estimate_pairwise_difference_correlation(fixed_beta, fixed_selected_se)
    )

    candidate_beta <- candidate_se * matrix(
      stats::rnorm(length(candidate_se)),
      nrow = nrow(candidate_se)
    )
    candidate_centered <- weighted_center_standardize(
      candidate_beta,
      candidate_se
    )
    selected <- order(candidate_centered$residual_score)[seq_len(top_n)]
    selected_beta <- candidate_beta[selected, , drop = FALSE]
    selected_se <- candidate_se[selected, , drop = FALSE]
    arrays[[2]][[1]][replication, ] <- lag_average_correlation(
      estimate_direct_centered_correlation(selected_beta, selected_se)
    )
    arrays[[2]][[2]][replication, ] <- lag_average_correlation(
      estimate_pairwise_difference_correlation(selected_beta, selected_se)
    )
  }

  summaries <- list()
  output_index <- 1L
  for (benchmark_index in seq_along(benchmark_names)) {
    for (estimator_index in seq_along(estimator_names)) {
      summaries[[output_index]] <- summarize_resampled_lags(
        arrays[[benchmark_index]][[estimator_index]],
        estimator = estimator_names[estimator_index],
        benchmark = benchmark_names[benchmark_index]
      )
      output_index <- output_index + 1L
    }
  }
  do.call(rbind, summaries)
}

correlation_matrix_to_long <- function(correlation,
                                       estimator,
                                       matrix_version,
                                       se_scale = "t-adjusted") {
  correlation <- validate_finite_matrix(correlation, "correlation")
  if (nrow(correlation) != ncol(correlation)) {
    stop("correlation must be square.")
  }
  time_labels <- colnames(correlation)
  if (is.null(time_labels)) {
    time_labels <- as.character(seq_len(ncol(correlation)) - 1L)
  } else {
    time_labels <- sub("^time_", "", time_labels)
  }
  grid <- expand.grid(
    time_a = time_labels,
    time_b = time_labels,
    stringsAsFactors = FALSE
  )
  grid$correlation <- as.vector(correlation)
  grid$estimator <- estimator
  grid$matrix_version <- matrix_version
  grid$se_scale <- se_scale
  grid[, c(
    "estimator", "matrix_version", "se_scale",
    "time_a", "time_b", "correlation"
  )]
}
