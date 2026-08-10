#!/usr/bin/env Rscript

# Estimate raw pairwise-difference residual correlations across lfdr thresholds.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(args, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(args[equals_hit[1]], nchar(equals_prefix) + 1L))
  }
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1] == length(args)) {
    return(default)
  }
  args[hit[1] + 1L]
}

parse_numeric_list <- function(value, name) {
  parsed <- as.numeric(strsplit(value, ",", fixed = TRUE)[[1]])
  if (length(parsed) == 0L || any(!is.finite(parsed))) {
    stop("Invalid ", name, ".")
  }
  parsed
}

threshold_id <- function(threshold) {
  sub("\\.", "p", format(threshold, nsmall = 2, trim = TRUE))
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "r4_correlated_errors",
  "real_data_correlation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "covariance_estimation",
  "pairwise_lfdr_threshold_helpers.R"
))

thresholds <- parse_numeric_list(
  get_arg("--lfdr-thresholds", "0.97,0.96,0.95"),
  "lfdr thresholds"
)
thresholds <- sort(unique(thresholds), decreasing = TRUE)
n_bootstrap <- as.integer(get_arg("--bootstrap-reps", "1000"))
bootstrap_seed_start <- as.integer(
  get_arg("--bootstrap-seed-start", "20260821")
)
output_id <- get_arg(
  "--output-id",
  "pairwise_lfdr_threshold_correlation"
)
if (length(thresholds) < 1L || any(thresholds <= 0 | thresholds >= 1) ||
    is.na(n_bootstrap) || n_bootstrap < 20L ||
    is.na(bootstrap_seed_start) || !nzchar(output_id)) {
  stop("Invalid pairwise threshold-exploration arguments.")
}

fit_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "fash_fit1_update.RData"
)
if (!file.exists(fit_path)) {
  stop("The BF-adjusted real-data FASH fit is missing: ", fit_path)
}
fit_info <- file.info(fit_path)
fit_environment <- new.env(parent = emptyenv())
loaded_names <- load(fit_path, envir = fit_environment)
if (!identical(loaded_names, "fash_fit1_update")) {
  stop("The real-data fit file did not contain only fash_fit1_update.")
}
fash_fit <- fit_environment$fash_fit1_update
if (is.null(fash_fit$fash_data$data_list) ||
    length(fash_fit$lfdr) != length(fash_fit$fash_data$data_list) ||
    any(!is.finite(fash_fit$lfdr))) {
  stop("The fitted FASH object has incomplete data or lfdr values.")
}
pair_keys <- names(fash_fit$fash_data$data_list)

selections <- list()
matrices <- list()
bootstrap_results <- list()
selection_rows <- list()
count_rows <- list()
diagnostic_rows <- list()
matrix_rows <- list()
lag_rows <- list()
bootstrap_interval_rows <- list()

for (threshold_index in seq_along(thresholds)) {
  threshold <- thresholds[threshold_index]
  id <- paste0("lfdr_", threshold_id(threshold))
  message("Selecting one pair per gene with lfdr > ", threshold, ".")
  selection <- select_gene_representatives_above_lfdr(
    pair_keys,
    fash_fit$lfdr,
    threshold
  )
  extracted <- extract_fash_beta_se(
    fash_fit,
    selection$selected_indices
  )
  if (!isTRUE(all.equal(extracted$time_grid, 0:15))) {
    stop("The selected FASH data do not use the expected time grid 0:15.")
  }

  bootstrap_seed <- bootstrap_seed_start + threshold_index - 1L
  message(
    "Running ",
    n_bootstrap,
    " gene-bootstrap replications for lfdr > ",
    threshold,
    "."
  )
  bootstrap <- bootstrap_pairwise_lag_variogram(
    extracted$beta_hat,
    extracted$se,
    n_bootstrap = n_bootstrap,
    seed = bootstrap_seed
  )
  correlation <- bootstrap$observed_matrix
  eigenvalues <- eigen(
    correlation,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  off_diagonal <- correlation[upper.tri(correlation)]
  lag_correlation <- lag_average_correlation(correlation)

  selected_metadata <- selection$selected
  selected_metadata$bootstrap_seed <- bootstrap_seed
  selection_rows[[id]] <- selected_metadata
  count_rows[[id]] <- data.frame(
    threshold = threshold,
    n_selected = nrow(selected_metadata),
    n_gene_representatives = selection$n_gene_representatives,
    minimum_selected_lfdr = min(selected_metadata$lfdr),
    maximum_selected_lfdr = max(selected_metadata$lfdr),
    stringsAsFactors = FALSE
  )
  diagnostic_rows[[id]] <- data.frame(
    threshold = threshold,
    n_selected = nrow(selected_metadata),
    minimum_eigenvalue = min(eigenvalues),
    maximum_eigenvalue = max(eigenvalues),
    n_negative_eigenvalues = sum(eigenvalues < -1e-10),
    mean_off_diagonal_correlation = mean(off_diagonal),
    minimum_off_diagonal_correlation = min(off_diagonal),
    maximum_off_diagonal_correlation = max(off_diagonal),
    mean_lag1_correlation = lag_correlation[1],
    mean_lag15_correlation = lag_correlation[15],
    stringsAsFactors = FALSE
  )
  matrix_long <- correlation_matrix_to_long(
    correlation,
    estimator = "Pairwise difference",
    matrix_version = "Raw unprojected estimate",
    se_scale = "t-adjusted"
  )
  matrix_long$threshold <- threshold
  matrix_rows[[id]] <- matrix_long[, c(
    "threshold",
    "estimator",
    "matrix_version",
    "se_scale",
    "time_a",
    "time_b",
    "correlation"
  )]
  lag_rows[[id]] <- data.frame(
    threshold = threshold,
    n_selected = nrow(selected_metadata),
    lag = seq_along(lag_correlation),
    mean_correlation = lag_correlation,
    semivariogram = 1 - lag_correlation,
    stringsAsFactors = FALSE
  )
  bootstrap_summary <- bootstrap$summary
  bootstrap_summary$threshold <- threshold
  bootstrap_summary$n_selected <- nrow(selected_metadata)
  bootstrap_summary$bootstrap_seed <- bootstrap_seed
  bootstrap_interval_rows[[id]] <- bootstrap_summary[, c(
    "threshold",
    "n_selected",
    "bootstrap_seed",
    "lag",
    "observed_correlation",
    "observed_semivariogram",
    "bootstrap_mean_correlation",
    "correlation_ci_lower",
    "correlation_ci_upper",
    "bootstrap_mean_semivariogram",
    "semivariogram_ci_lower",
    "semivariogram_ci_upper",
    "n_bootstrap"
  )]

  selections[[id]] <- selection
  matrices[[id]] <- correlation
  bootstrap_results[[id]] <- list(
    seed = bootstrap_seed,
    draws = bootstrap$draws
  )
}

for (threshold_index in seq_len(length(thresholds) - 1L)) {
  tighter <- selections[[threshold_index]]$selected$pair_key
  looser <- selections[[threshold_index + 1L]]$selected$pair_key
  if (!all(tighter %in% looser)) {
    stop("The lfdr-threshold selections are not nested.")
  }
}

selected_units <- do.call(rbind, selection_rows)
selection_counts <- do.call(rbind, count_rows)
matrix_diagnostics <- do.call(rbind, diagnostic_rows)
correlation_matrices_long <- do.call(rbind, matrix_rows)
lag_variogram_summaries <- do.call(rbind, lag_rows)
bootstrap_lag_variogram_intervals <- do.call(
  rbind,
  bootstrap_interval_rows
)
rownames(selected_units) <- NULL
rownames(selection_counts) <- NULL
rownames(matrix_diagnostics) <- NULL
rownames(correlation_matrices_long) <- NULL
rownames(lag_variogram_summaries) <- NULL
rownames(bootstrap_lag_variogram_intervals) <- NULL

configuration <- list(
  analysis_id = output_id,
  created_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  fit_path = fit_path,
  fit_object = "fash_fit1_update",
  fit_size_bytes = unname(fit_info$size),
  fit_mtime = format(fit_info$mtime, tz = "UTC", usetz = TRUE),
  n_fash_pairs = length(pair_keys),
  n_gene_representatives = selections[[1]]$n_gene_representatives,
  time_grid = 0:15,
  thresholds = thresholds,
  threshold_rule = "Strictly greater than threshold",
  selection = paste(
    "Highest BF-adjusted FASH lfdr variant per gene, then all",
    "gene representatives with lfdr above the threshold"
  ),
  estimator = paste(
    "OLS through the origin for 1 - squared beta-hat difference divided by",
    "the sum of squared t-adjusted SEs, against",
    "2*s_t*s_s/(s_t^2+s_s^2)"
  ),
  matrix_version = "Raw pairwise-difference estimate; no PD projection",
  se_scale = "t-adjusted",
  n_bootstrap = n_bootstrap,
  bootstrap_seed_start = bootstrap_seed_start,
  r_version = R.version.string,
  fashr_version = as.character(utils::packageVersion("fashr"))
)

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  output_id
)
summary_dir <- file.path(output_dir, "summary")
invisible(lapply(
  c(output_dir, summary_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

analysis_result <- list(
  configuration = configuration,
  selections = selections,
  matrices = matrices,
  bootstrap_results = bootstrap_results,
  selected_units = selected_units,
  selection_counts = selection_counts,
  matrix_diagnostics = matrix_diagnostics,
  correlation_matrices_long = correlation_matrices_long,
  lag_variogram_summaries = lag_variogram_summaries,
  bootstrap_lag_variogram_intervals = bootstrap_lag_variogram_intervals
)
saveRDS(
  configuration,
  file.path(output_dir, "configuration.rds")
)
saveRDS(
  analysis_result,
  file.path(output_dir, "pairwise_lfdr_threshold_analysis.rds")
)
saveRDS(
  matrices,
  file.path(output_dir, "raw_pairwise_correlation_matrices.rds")
)
write_csv(selected_units, file.path(summary_dir, "selected_units.csv"))
write_csv(selection_counts, file.path(summary_dir, "selection_counts.csv"))
write_csv(matrix_diagnostics, file.path(summary_dir, "matrix_diagnostics.csv"))
write_csv(
  correlation_matrices_long,
  file.path(summary_dir, "correlation_matrices_long.csv")
)
write_csv(
  lag_variogram_summaries,
  file.path(summary_dir, "lag_variogram_summaries.csv")
)
write_csv(
  bootstrap_lag_variogram_intervals,
  file.path(summary_dir, "bootstrap_lag_variogram_intervals.csv")
)

cat("\nSelection counts:\n")
print(selection_counts)
cat("\nRaw matrix diagnostics:\n")
print(matrix_diagnostics)
cat("\nPairwise lfdr-threshold exploration completed: ", output_dir, "\n", sep = "")
