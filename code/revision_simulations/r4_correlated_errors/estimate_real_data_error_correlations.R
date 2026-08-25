#!/usr/bin/env Rscript

# Estimate full time-correlation patterns from the most null-like real-data eQTLs.

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

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

as_flag <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

sha256_file <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  command <- if (identical(Sys.info()[["sysname"]], "Darwin")) {
    "shasum"
  } else {
    "sha256sum"
  }
  arguments <- if (identical(command, "shasum")) {
    c("-a", "256", path)
  } else {
    path
  }
  output <- system2(command, arguments, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Unable to compute SHA-256 for ", path, ".")
  }
  sub("[[:space:]].*$", "", output[[1L]])
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "r4_correlated_errors",
  "real_data_correlation_helpers.R"
))

top_n <- as.integer(get_arg("--top-n", "500"))
n_bootstrap <- as.integer(get_arg("--bootstrap-reps", "1000"))
n_benchmark <- as.integer(get_arg("--benchmark-reps", "1000"))
bootstrap_seed <- as.integer(get_arg("--bootstrap-seed", "20260805"))
benchmark_seed <- as.integer(get_arg("--benchmark-seed", "20260806"))
output_id <- get_arg(
  "--output-id",
  "r4_null_like_top500_full_correlations"
)
fit_path_argument <- get_arg("--fit-path", "")
output_dir_argument <- get_arg("--output-dir", "")
expected_fit_sha256 <- get_arg("--expected-fit-sha256", "")
expected_fashr_version <- get_arg("--expected-fashr-version", "")
expected_fashr_remote_sha <- get_arg("--expected-fashr-remote-sha", "")
overwrite <- as_flag(get_arg("--overwrite", "false"))
if (is.na(top_n) || top_n < 20L || is.na(n_bootstrap) || n_bootstrap < 20L ||
    is.na(n_benchmark) || n_benchmark < 20L || is.na(bootstrap_seed) ||
    is.na(benchmark_seed) || !nzchar(output_id)) {
  stop("Invalid real-data correlation-estimation arguments.")
}

fit_path <- if (nzchar(fit_path_argument)) {
  fit_path_argument
} else {
  file.path(
    workflowr_root,
    "output",
    "dynamic_eQTL_real",
    "fash_fit1_update.RData"
  )
}
if (!file.exists(fit_path)) {
  stop("The BF-adjusted real-data FASH fit is missing: ", fit_path)
}
fit_path <- normalizePath(fit_path, winslash = "/", mustWork = TRUE)
fit_sha256 <- sha256_file(fit_path)
if (nzchar(expected_fit_sha256) && !identical(fit_sha256, expected_fit_sha256)) {
  stop(
    "Expected real-data fit SHA-256 ", expected_fit_sha256,
    "; found ", fit_sha256, "."
  )
}
if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required to load the real-data fit.")
}
fashr_description <- utils::packageDescription("fashr")
package_provenance <- list(
  package = "fashr",
  version = as.character(utils::packageVersion("fashr")),
  remote_sha = if (is.null(fashr_description$RemoteSha)) {
    NA_character_
  } else {
    as.character(fashr_description$RemoteSha)
  },
  library_path = normalizePath(
    find.package("fashr"), winslash = "/", mustWork = TRUE
  ),
  r_version = R.version.string,
  platform = R.version$platform
)
if (nzchar(expected_fashr_version) &&
    !identical(package_provenance$version, expected_fashr_version)) {
  stop(
    "Expected fashr ", expected_fashr_version,
    "; found ", package_provenance$version, "."
  )
}
if (nzchar(expected_fashr_remote_sha) &&
    !identical(package_provenance$remote_sha, expected_fashr_remote_sha)) {
  stop(
    "Expected fashr RemoteSha ", expected_fashr_remote_sha,
    "; found ", package_provenance$remote_sha, "."
  )
}
fit_environment <- new.env(parent = emptyenv())
loaded_names <- load(fit_path, envir = fit_environment)
if (!identical(loaded_names, "fash_fit1_update")) {
  stop("The real-data fit file did not contain only fash_fit1_update.")
}
fash_fit <- fit_environment$fash_fit1_update
pair_keys <- names(fash_fit$fash_data$data_list)
selection <- select_highest_lfdr_per_gene(
  pair_keys = pair_keys,
  lfdr = fash_fit$lfdr,
  top_n = top_n
)
selected_data <- extract_fash_beta_se(
  fash_fit,
  selection$selected_indices
)
representative_data <- extract_fash_beta_se(
  fash_fit,
  selection$representative_indices
)
if (!isTRUE(all.equal(selected_data$time_grid, 0:15))) {
  stop("The real-data time grid is not 0:15.")
}

adjusted_se <- selected_data$se
beta_hat <- selected_data$beta_hat
residual_df <- c(19, 19, 16, 19, 16, 19, 19, 19,
                 19, 19, 19, 19, 19, 18, 19, 19) - 7
raw_regression_se <- invert_t_to_normal_se(
  beta_hat,
  adjusted_se,
  df = residual_df
)

matrix_inputs <- list(
  direct_adjusted = list(
    estimator = "Direct centered",
    se_scale = "t-adjusted",
    matrix = estimate_direct_centered_correlation(beta_hat, adjusted_se)
  ),
  pairwise_adjusted = list(
    estimator = "Pairwise difference",
    se_scale = "t-adjusted",
    matrix = estimate_pairwise_difference_correlation(beta_hat, adjusted_se)
  ),
  direct_raw_se = list(
    estimator = "Direct centered",
    se_scale = "raw regression",
    matrix = estimate_direct_centered_correlation(beta_hat, raw_regression_se)
  ),
  pairwise_raw_se = list(
    estimator = "Pairwise difference",
    se_scale = "raw regression",
    matrix = estimate_pairwise_difference_correlation(beta_hat, raw_regression_se)
  )
)

matrix_results <- lapply(matrix_inputs, function(input) {
  projection <- project_to_positive_definite_correlation(input$matrix)
  projection$estimator <- input$estimator
  projection$se_scale <- input$se_scale
  projection
})

projection_diagnostics <- do.call(rbind, lapply(
  names(matrix_results),
  function(name) {
    result <- matrix_results[[name]]
    data.frame(
      matrix_id = name,
      estimator = result$estimator,
      se_scale = result$se_scale,
      result$diagnostics,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
))
rownames(projection_diagnostics) <- NULL

lag_summaries <- do.call(rbind, lapply(matrix_results, function(result) {
  rbind(
    make_lag_summary(
      result$raw,
      estimator = result$estimator,
      matrix_version = "Raw estimate",
      se_scale = result$se_scale
    ),
    make_lag_summary(
      result$projected,
      estimator = result$estimator,
      matrix_version = "Nearest-PD",
      se_scale = result$se_scale
    )
  )
}))
rownames(lag_summaries) <- NULL

message("Running ", n_bootstrap, " selected-gene bootstrap replications.")
bootstrap_lags <- bootstrap_correlation_lags(
  beta_hat,
  adjusted_se,
  n_bootstrap = n_bootstrap,
  seed = bootstrap_seed
)
observed_lags <- rbind(
  make_lag_summary(
    matrix_results$direct_adjusted$raw,
    "Direct centered",
    "Raw estimate"
  ),
  make_lag_summary(
    matrix_results$pairwise_adjusted$raw,
    "Pairwise difference",
    "Raw estimate"
  )
)
bootstrap_lags <- merge(
  bootstrap_lags,
  observed_lags[, c(
    "estimator", "lag", "mean_correlation", "semivariogram"
  )],
  by = c("estimator", "lag"),
  suffixes = c("_bootstrap", "_observed"),
  sort = FALSE
)
bootstrap_lags <- bootstrap_lags[
  order(bootstrap_lags$estimator, bootstrap_lags$lag),
]
rownames(bootstrap_lags) <- NULL

message("Running ", n_benchmark, " matched independence benchmark replications.")
independence_benchmark <- run_independence_benchmark(
  candidate_se = representative_data$se,
  fixed_selected_se = adjusted_se,
  top_n = top_n,
  n_replications = n_benchmark,
  seed = benchmark_seed
)

matrix_long <- do.call(rbind, lapply(matrix_results, function(result) {
  rbind(
    correlation_matrix_to_long(
      result$raw,
      result$estimator,
      "Raw estimate",
      result$se_scale
    ),
    correlation_matrix_to_long(
      result$projected,
      result$estimator,
      "Nearest-PD",
      result$se_scale
    )
  )
}))
rownames(matrix_long) <- NULL

configuration <- list(
  output_id = output_id,
  fit_path = fit_path,
  fit_sha256 = fit_sha256,
  fit_object = "fash_fit1_update",
  package_provenance = package_provenance,
  selection = "Highest BF-adjusted lfdr variant per gene, then top lfdr genes",
  top_n = top_n,
  n_gene_representatives = selection$n_genes,
  time_grid = selected_data$time_grid,
  residual_df = residual_df,
  n_bootstrap = n_bootstrap,
  bootstrap_seed = bootstrap_seed,
  n_independence_benchmark = n_benchmark,
  benchmark_seed = benchmark_seed,
  pairwise_estimator = paste(
    "OLS through the origin for 1 - standardized squared difference",
    "against 2*s_t*s_s/(s_t^2+s_s^2)"
  ),
  projection = "Matrix::nearPD(corr = TRUE, keepDiag = TRUE)"
)

output_dir <- if (nzchar(output_dir_argument)) {
  output_dir_argument
} else {
  file.path(
    workflowr_root,
    "output",
    "revision_simulations",
    "real_data",
    output_id
  )
}
if ((dir.exists(output_dir) || file.exists(output_dir)) && !overwrite) {
  stop("Refusing to overwrite real-data R4 output: ", output_dir)
}
summary_dir <- file.path(output_dir, "summary")
invisible(lapply(
  c(output_dir, summary_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

analysis_result <- list(
  configuration = configuration,
  selected_units = selection$selected,
  beta_hat = beta_hat,
  adjusted_se = adjusted_se,
  raw_regression_se = raw_regression_se,
  centered_standardized_residual = weighted_center_standardize(
    beta_hat,
    adjusted_se
  )$standardized_residual,
  matrices = matrix_results,
  projection_diagnostics = projection_diagnostics,
  lag_summaries = lag_summaries,
  bootstrap_lags = bootstrap_lags,
  independence_benchmark = independence_benchmark
)
saveRDS(analysis_result, file.path(output_dir, "real_data_correlation_analysis.rds"))
saveRDS(configuration, file.path(output_dir, "configuration.rds"))
saveRDS(
  list(
    direct_centered = matrix_results$direct_adjusted$projected,
    pairwise_difference = matrix_results$pairwise_adjusted$projected
  ),
  file.path(output_dir, "simulation_correlation_matrices.rds")
)
write_csv(selection$selected, file.path(summary_dir, "selected_units.csv"))
write_csv(projection_diagnostics, file.path(
  summary_dir,
  "projection_diagnostics.csv"
))
write_csv(lag_summaries, file.path(summary_dir, "lag_summaries.csv"))
write_csv(bootstrap_lags, file.path(summary_dir, "bootstrap_lag_intervals.csv"))
write_csv(independence_benchmark, file.path(
  summary_dir,
  "independence_selection_benchmark.csv"
))
write_csv(matrix_long, file.path(summary_dir, "correlation_matrices_long.csv"))

cat("\nSelected-unit lfdr range:\n")
print(range(selection$selected$lfdr))
cat("\nProjection diagnostics:\n")
print(projection_diagnostics)
cat("\nAdjusted-SE lag-1 summaries:\n")
print(lag_summaries[
  lag_summaries$se_scale == "t-adjusted" & lag_summaries$lag == 1L,
])
cat("\nIndependence-benchmark lag-1 summaries:\n")
print(independence_benchmark[independence_benchmark$lag == 1L, ])
