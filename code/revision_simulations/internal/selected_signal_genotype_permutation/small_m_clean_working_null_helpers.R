# Helpers for the small-M clean working-model null sensitivity analysis.

simulate_clean_working_null <- function(target_data_list,
                                        target_se_list,
                                        source_unit_keys,
                                        seed) {
  n_units <- length(target_data_list)
  source_unit_keys <- as.character(source_unit_keys)
  seed <- as.integer(seed)
  if (n_units < 1L || length(target_se_list) != n_units ||
      length(source_unit_keys) != n_units || any(!nzchar(source_unit_keys)) ||
      anyDuplicated(source_unit_keys) || length(seed) != 1L || is.na(seed)) {
    stop("Invalid clean working-model null inputs.")
  }

  data_valid <- vapply(target_data_list, function(value) {
    is.data.frame(value) && all(c("y", "x", "offset") %in% names(value)) &&
      nrow(value) >= 2L && all(is.finite(value$y)) &&
      all(is.finite(value$x)) && all(is.finite(value$offset))
  }, logical(1))
  data_lengths <- vapply(target_data_list, nrow, integer(1))
  se_valid <- vapply(seq_len(n_units), function(index) {
    value <- target_se_list[[index]]
    is.numeric(value) && length(value) == data_lengths[index] &&
      all(is.finite(value)) && all(value > 0)
  }, logical(1))
  if (!all(data_valid) || !all(se_valid) ||
      length(unique(data_lengths)) != 1L) {
    stop("Clean working-model null data or SE vectors are invalid.")
  }

  n_time <- unique(data_lengths)
  set.seed(seed)
  z_matrix <- matrix(
    stats::rnorm(n_units * n_time),
    nrow = n_units,
    ncol = n_time,
    byrow = TRUE,
    dimnames = list(source_unit_keys, paste0("time_", seq_len(n_time)))
  )
  null_unit_keys <- paste0("clean_null::", source_unit_keys)
  data_list <- lapply(seq_len(n_units), function(index) {
    result <- target_data_list[[index]]
    result$y <- as.numeric(target_se_list[[index]]) * z_matrix[index, ]
    result
  })
  se_list <- lapply(target_se_list, as.numeric)
  names(data_list) <- null_unit_keys
  names(se_list) <- null_unit_keys

  if (any(!is.finite(z_matrix)) ||
      any(!vapply(data_list, function(value) all(is.finite(value$y)),
                  logical(1)))) {
    stop("Clean working-model null simulation produced invalid values.")
  }

  list(
    data_list = data_list,
    se_list = se_list,
    z_matrix = z_matrix,
    source_unit_keys = source_unit_keys,
    null_unit_keys = null_unit_keys,
    seed = seed
  )
}

summarize_clean_null_z <- function(z_matrix) {
  z_matrix <- as.matrix(z_matrix)
  if (nrow(z_matrix) < 2L || ncol(z_matrix) < 2L ||
      any(!is.finite(z_matrix))) {
    stop("Invalid clean-null z matrix.")
  }
  correlation <- stats::cor(z_matrix)
  if (any(!is.finite(correlation))) {
    stop("Clean-null z correlations are not finite.")
  }
  time_labels <- colnames(z_matrix)
  if (is.null(time_labels)) {
    time_labels <- paste0("time_", seq_len(ncol(z_matrix)))
  }
  off_diagonal <- correlation[row(correlation) != col(correlation)]
  correlation_long <- expand.grid(
    time_1 = time_labels,
    time_2 = time_labels,
    stringsAsFactors = FALSE
  )
  correlation_long$correlation <- as.numeric(correlation)

  list(
    overall = data.frame(
      n_units = nrow(z_matrix),
      n_time = ncol(z_matrix),
      overall_mean = mean(z_matrix),
      overall_sd = stats::sd(as.numeric(z_matrix)),
      mean_off_diagonal_correlation = mean(off_diagonal),
      mean_absolute_off_diagonal_correlation = mean(abs(off_diagonal)),
      max_absolute_off_diagonal_correlation = max(abs(off_diagonal)),
      stringsAsFactors = FALSE
    ),
    time_summary = data.frame(
      time = time_labels,
      mean_z = colMeans(z_matrix),
      sd_z = apply(z_matrix, 2L, stats::sd),
      stringsAsFactors = FALSE
    ),
    correlation_long = correlation_long
  )
}
