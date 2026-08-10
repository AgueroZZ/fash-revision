# Helpers for the internal zero-intercept z-score correlation exploration.

load_bf_adjusted_real_fash <- function(path) {
  if (!file.exists(path)) {
    stop("The BF-adjusted real-data FASH fit is missing: ", path)
  }
  fit_environment <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = fit_environment)
  if (!identical(loaded_names, "fash_fit1_update")) {
    stop("The fit file must contain only fash_fit1_update.")
  }
  fit <- fit_environment$fash_fit1_update
  required_fields <- c("fash_data", "lfdr")
  if (!all(required_fields %in% names(fit)) ||
      is.null(fit$fash_data$data_list) || is.null(fit$fash_data$S) ||
      length(fit$lfdr) != length(fit$fash_data$data_list)) {
    stop("The FASH object does not contain the required data and lfdr fields.")
  }
  fit
}

make_gene_index <- function(pair_keys) {
  gene_id <- parse_gene_ids(pair_keys)
  split(seq_along(pair_keys), factor(gene_id, levels = sort(unique(gene_id))))
}

select_random_variant_per_gene <- function(pair_keys,
                                           seed,
                                           gene_index = NULL) {
  pair_keys <- as.character(pair_keys)
  seed <- as.integer(seed)
  if (length(pair_keys) < 2L || any(!nzchar(pair_keys)) ||
      length(seed) != 1L || is.na(seed)) {
    stop("Invalid pair keys or seed for random variant selection.")
  }
  if (is.null(gene_index)) {
    gene_index <- make_gene_index(pair_keys)
  }
  if (!is.list(gene_index) || length(gene_index) < 2L ||
      any(lengths(gene_index) < 1L)) {
    stop("gene_index must contain at least two non-empty gene groups.")
  }
  set.seed(seed)
  selected_indices <- vapply(
    gene_index,
    function(indices) {
      if (length(indices) == 1L) indices else sample(indices, size = 1L)
    },
    integer(1)
  )
  selected_gene <- parse_gene_ids(pair_keys[selected_indices])
  if (anyDuplicated(selected_gene) || length(selected_indices) != length(gene_index)) {
    stop("Random selection did not produce exactly one variant per gene.")
  }
  data.frame(
    seed = seed,
    fash_index = selected_indices,
    pair_key = pair_keys[selected_indices],
    gene_id = selected_gene,
    variant_id = sub("^[^_]+_", "", pair_keys[selected_indices]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

stratified_lfdr_gallery_sample <- function(metadata,
                                           n_per_quartile = 25L,
                                           seed = 20260891L) {
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  n_per_quartile <- as.integer(n_per_quartile)
  seed <- as.integer(seed)
  required_columns <- c("pair_key", "lfdr")
  if (!all(required_columns %in% names(metadata)) ||
      nrow(metadata) < 4L || any(!nzchar(metadata$pair_key)) ||
      anyDuplicated(metadata$pair_key) || any(!is.finite(metadata$lfdr)) ||
      any(metadata$lfdr < 0 | metadata$lfdr > 1) ||
      length(n_per_quartile) != 1L || is.na(n_per_quartile) ||
      n_per_quartile < 1L || length(seed) != 1L || is.na(seed)) {
    stop("Invalid metadata, n_per_quartile, or seed for the lfdr gallery.")
  }

  retained_order <- order(
    -metadata$lfdr,
    metadata$pair_key,
    method = "radix"
  )
  ranked <- metadata[retained_order, , drop = FALSE]
  ranked$retained_lfdr_rank <- seq_len(nrow(ranked))
  ranked$lfdr_quartile <- pmin(
    4L,
    ceiling(4 * ranked$retained_lfdr_rank / nrow(ranked))
  )
  quartile_sizes <- tabulate(ranked$lfdr_quartile, nbins = 4L)
  if (any(quartile_sizes < n_per_quartile)) {
    stop("Every lfdr quartile must contain at least n_per_quartile units.")
  }

  set.seed(seed)
  selected_rows <- unlist(lapply(1:4, function(quartile) {
    candidates <- which(ranked$lfdr_quartile == quartile)
    sample(candidates, size = n_per_quartile, replace = FALSE)
  }))
  selected <- ranked[selected_rows, , drop = FALSE]
  selected <- selected[order(
    selected$lfdr_quartile,
    selected$retained_lfdr_rank,
    method = "radix"
  ), , drop = FALSE]
  selected$gallery_page_position <- ave(
    selected$retained_lfdr_rank,
    selected$lfdr_quartile,
    FUN = seq_along
  )
  rownames(selected) <- NULL
  selected
}

extract_z_matrix <- function(fash_fit, indices) {
  extracted <- extract_fash_beta_se(fash_fit, indices)
  z <- extracted$beta_hat / extracted$se
  if (any(!is.finite(z))) {
    stop("The extracted z-score matrix contains non-finite values.")
  }
  list(
    beta_hat = extracted$beta_hat,
    adjusted_se = extracted$se,
    z = z,
    time_grid = extracted$time_grid
  )
}

filter_zero_intercept_z <- function(z, threshold = 2) {
  z <- validate_finite_matrix(z, "z")
  threshold <- as.numeric(threshold)
  if (length(threshold) != 1L || !is.finite(threshold) || threshold <= 0) {
    stop("threshold must be one positive finite number.")
  }
  max_absolute_z <- apply(abs(z), 1L, max)
  keep <- max_absolute_z < threshold
  list(
    threshold = threshold,
    keep = keep,
    selected_indices = which(keep),
    max_absolute_z = max_absolute_z,
    n_candidates = nrow(z),
    n_selected = sum(keep),
    selected_fraction = mean(keep)
  )
}

filter_zero_intercept_z_with_mean <- function(z,
                                              max_threshold = 2,
                                              mean_z_threshold = 2) {
  z <- validate_finite_matrix(z, "z")
  max_threshold <- as.numeric(max_threshold)
  mean_z_threshold <- as.numeric(mean_z_threshold)
  if (length(max_threshold) != 1L || !is.finite(max_threshold) ||
      max_threshold <= 0 || length(mean_z_threshold) != 1L ||
      !is.finite(mean_z_threshold) || mean_z_threshold <= 0) {
    stop("Both thresholds must be positive finite numbers.")
  }

  max_absolute_z <- apply(abs(z), 1L, max)
  mean_z_score <- sqrt(ncol(z)) * rowMeans(z)
  keep_maximum <- max_absolute_z < max_threshold
  keep_mean_z <- abs(mean_z_score) < mean_z_threshold
  keep <- keep_maximum & keep_mean_z

  list(
    max_threshold = max_threshold,
    mean_z_threshold = mean_z_threshold,
    keep_maximum = keep_maximum,
    keep_mean_z = keep_mean_z,
    keep = keep,
    selected_indices = which(keep),
    max_absolute_z = max_absolute_z,
    mean_z_score = mean_z_score,
    n_candidates = nrow(z),
    n_selected_by_maximum = sum(keep_maximum),
    n_selected = sum(keep),
    selected_fraction_by_maximum = mean(keep_maximum),
    selected_fraction = mean(keep)
  )
}

validate_correlation_matrix <- function(correlation,
                                        name = "correlation",
                                        tolerance = 1e-8) {
  correlation <- validate_finite_matrix(correlation, name)
  if (nrow(correlation) != ncol(correlation) ||
      max(abs(correlation - t(correlation))) > tolerance ||
      max(abs(diag(correlation) - 1)) > tolerance ||
      max(abs(correlation)) > 1 + tolerance) {
    stop(name, " is not a valid finite correlation matrix.")
  }
  correlation
}

estimate_zero_intercept_correlation <- function(z_nullish) {
  z_nullish <- validate_finite_matrix(z_nullish, "z_nullish")
  if (nrow(z_nullish) <= ncol(z_nullish)) {
    stop("z_nullish must contain more rows than time points.")
  }
  sample_correlation <- stats::cor(z_nullish)
  second_moment <- crossprod(z_nullish) / nrow(z_nullish)
  normalized_second_moment <- stats::cov2cor(second_moment)
  sample_correlation <- validate_correlation_matrix(
    sample_correlation,
    "sample_correlation"
  )
  normalized_second_moment <- validate_correlation_matrix(
    normalized_second_moment,
    "normalized_second_moment"
  )
  sample_projection <- project_to_positive_definite_correlation(
    sample_correlation
  )
  normalized_projection <- project_to_positive_definite_correlation(
    normalized_second_moment
  )
  list(
    sample_correlation = sample_correlation,
    second_moment = second_moment,
    normalized_second_moment = normalized_second_moment,
    sample_projection = sample_projection,
    normalized_projection = normalized_projection,
    column_means = colMeans(z_nullish),
    second_moment_diagonal = diag(second_moment),
    n_selected = nrow(z_nullish),
    maximum_correlation_difference = max(abs(
      sample_correlation - normalized_second_moment
    ))
  )
}

summarize_correlation_lags <- function(correlation,
                                       design,
                                       threshold,
                                       estimator = "mashr cor(z)",
                                       seed = NA_integer_) {
  correlation <- validate_correlation_matrix(correlation)
  lag_correlation <- lag_average_correlation(correlation)
  data.frame(
    design = as.character(design),
    threshold = as.numeric(threshold),
    estimator = as.character(estimator),
    seed = as.integer(seed),
    lag = seq_along(lag_correlation),
    mean_correlation = as.numeric(lag_correlation),
    semivariogram = 1 - as.numeric(lag_correlation),
    stringsAsFactors = FALSE
  )
}

summarize_adjacent_correlations <- function(correlation,
                                            design,
                                            threshold,
                                            estimator = "mashr cor(z)",
                                            seed = NA_integer_) {
  correlation <- validate_correlation_matrix(correlation)
  n_time <- ncol(correlation)
  data.frame(
    design = as.character(design),
    threshold = as.numeric(threshold),
    estimator = as.character(estimator),
    seed = as.integer(seed),
    time_a = seq_len(n_time - 1L) - 1L,
    time_b = seq_len(n_time - 1L),
    correlation = correlation[cbind(seq_len(n_time - 1L), 2:n_time)],
    stringsAsFactors = FALSE
  )
}

correlation_to_long <- function(correlation,
                                design,
                                threshold,
                                estimator,
                                seed = NA_integer_) {
  correlation <- validate_correlation_matrix(correlation)
  grid <- expand.grid(
    time_a = seq_len(nrow(correlation)) - 1L,
    time_b = seq_len(ncol(correlation)) - 1L
  )
  grid$design <- as.character(design)
  grid$threshold <- as.numeric(threshold)
  grid$estimator <- as.character(estimator)
  grid$seed <- as.integer(seed)
  grid$correlation <- correlation[cbind(grid$time_a + 1L, grid$time_b + 1L)]
  grid[, c(
    "design", "threshold", "estimator", "seed",
    "time_a", "time_b", "correlation"
  )]
}

finite_square_matrix_to_long <- function(matrix,
                                         design,
                                         threshold,
                                         estimator,
                                         seed = NA_integer_) {
  matrix <- validate_finite_matrix(matrix, "matrix")
  if (nrow(matrix) != ncol(matrix)) {
    stop("matrix must be square.")
  }
  grid <- expand.grid(
    time_a = seq_len(nrow(matrix)) - 1L,
    time_b = seq_len(ncol(matrix)) - 1L
  )
  grid$design <- as.character(design)
  grid$threshold <- as.numeric(threshold)
  grid$estimator <- as.character(estimator)
  grid$seed <- as.integer(seed)
  grid$value <- matrix[cbind(grid$time_a + 1L, grid$time_b + 1L)]
  grid[, c(
    "design", "threshold", "estimator", "seed",
    "time_a", "time_b", "value"
  )]
}

summarize_resampled_values <- function(values,
                                       labels,
                                       probability = c(0.025, 0.975)) {
  values <- as.matrix(values)
  if (nrow(values) < 20L || any(!is.finite(values)) ||
      length(labels) != ncol(values)) {
    stop("Invalid resampled values or labels.")
  }
  data.frame(
    label = labels,
    mean = colMeans(values),
    median = apply(values, 2L, stats::median),
    lower = apply(values, 2L, stats::quantile, probs = probability[1]),
    upper = apply(values, 2L, stats::quantile, probs = probability[2]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

bootstrap_zero_intercept_correlations <- function(z_nullish,
                                                  n_bootstrap = 1000L,
                                                  seed = 20260810L) {
  z_nullish <- validate_finite_matrix(z_nullish, "z_nullish")
  n_bootstrap <- as.integer(n_bootstrap)
  seed <- as.integer(seed)
  if (nrow(z_nullish) <= ncol(z_nullish) || is.na(n_bootstrap) ||
      n_bootstrap < 20L || is.na(seed)) {
    stop("Invalid bootstrap inputs.")
  }
  n_time <- ncol(z_nullish)
  adjacent <- matrix(NA_real_, nrow = n_bootstrap, ncol = n_time - 1L)
  lags <- matrix(NA_real_, nrow = n_bootstrap, ncol = n_time - 1L)
  set.seed(seed)
  for (replication in seq_len(n_bootstrap)) {
    index <- sample.int(nrow(z_nullish), nrow(z_nullish), replace = TRUE)
    correlation <- stats::cor(z_nullish[index, , drop = FALSE])
    adjacent[replication, ] <- correlation[cbind(
      seq_len(n_time - 1L),
      2:n_time
    )]
    lags[replication, ] <- lag_average_correlation(correlation)
  }
  adjacent_summary <- summarize_resampled_values(
    adjacent,
    paste0(seq_len(n_time - 1L) - 1L, "-", seq_len(n_time - 1L))
  )
  adjacent_summary$summary_type <- "adjacent_pair"
  adjacent_summary$index <- seq_len(n_time - 1L)
  lag_summary <- summarize_resampled_values(
    lags,
    as.character(seq_len(n_time - 1L))
  )
  lag_summary$summary_type <- "lag_average"
  lag_summary$index <- seq_len(n_time - 1L)
  rbind(adjacent_summary, lag_summary)[, c(
    "summary_type", "index", "label", "mean", "median", "lower", "upper"
  )]
}

make_lag1_only_correlation <- function(n_time, rho) {
  n_time <- as.integer(n_time)
  rho <- as.numeric(rho)
  if (is.na(n_time) || n_time < 2L || length(rho) != 1L ||
      !is.finite(rho)) {
    stop("Invalid n_time or rho.")
  }
  correlation <- diag(n_time)
  correlation[cbind(seq_len(n_time - 1L), 2:n_time)] <- rho
  correlation[cbind(2:n_time, seq_len(n_time - 1L))] <- rho
  if (min(eigen(correlation, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("The requested lag-1-only matrix is not positive definite.")
  }
  correlation
}

simulate_truncated_correlation <- function(n_candidates,
                                           threshold,
                                           correlation,
                                           seed) {
  n_candidates <- as.integer(n_candidates)
  seed <- as.integer(seed)
  correlation <- validate_correlation_matrix(correlation)
  if (is.na(n_candidates) || n_candidates <= ncol(correlation) || is.na(seed)) {
    stop("Invalid truncated-correlation simulation inputs.")
  }
  set.seed(seed)
  independent <- matrix(
    stats::rnorm(n_candidates * ncol(correlation)),
    nrow = n_candidates,
    ncol = ncol(correlation)
  )
  z <- independent %*% chol(correlation)
  filtered <- filter_zero_intercept_z(z, threshold)
  if (filtered$n_selected <= ncol(z)) {
    stop("The truncation simulation retained too few rows.")
  }
  estimate <- estimate_zero_intercept_correlation(
    z[filtered$keep, , drop = FALSE]
  )
  list(
    n_selected = filtered$n_selected,
    sample_correlation = estimate$sample_correlation,
    normalized_second_moment = estimate$normalized_second_moment
  )
}
