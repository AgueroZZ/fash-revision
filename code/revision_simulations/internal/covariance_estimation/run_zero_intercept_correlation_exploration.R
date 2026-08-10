#!/usr/bin/env Rscript

# Run the internal mash-style zero-intercept correlation exploration.

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

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

summarize_interval <- function(x) {
  if (length(x) < 20L || any(!is.finite(x))) {
    stop("At least 20 finite values are required for an interval summary.")
  }
  c(
    mean = mean(x),
    median = stats::median(x),
    lower = stats::quantile(x, 0.025, names = FALSE),
    upper = stats::quantile(x, 0.975, names = FALSE)
  )
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
  "zero_intercept_correlation_helpers.R"
))

top_n <- as.integer(get_arg("--top-n", "1000"))
thresholds <- parse_numeric_list(
  get_arg("--threshold-list", "1.5,2,2.5"),
  "threshold list"
)
random_reps <- as.integer(get_arg("--random-reps", "50"))
random_seed_start <- as.integer(get_arg("--random-seed-start", "20260811"))
bootstrap_reps <- as.integer(get_arg("--bootstrap-reps", "1000"))
bootstrap_seed <- as.integer(get_arg("--bootstrap-seed", "20260810"))
benchmark_reps <- as.integer(get_arg("--benchmark-reps", "1000"))
benchmark_seed <- as.integer(get_arg("--benchmark-seed", "20260870"))
calibration_n <- as.integer(get_arg("--calibration-n", "200000"))
output_id <- get_arg("--output-id", "zero_intercept_correlation")

if (is.na(top_n) || top_n < 100L || any(thresholds <= 0) ||
    anyDuplicated(thresholds) || is.na(random_reps) || random_reps < 2L ||
    is.na(random_seed_start) || is.na(bootstrap_reps) || bootstrap_reps < 20L ||
    is.na(bootstrap_seed) || is.na(benchmark_reps) || benchmark_reps < 20L ||
    is.na(benchmark_seed) || is.na(calibration_n) || calibration_n < 1000L ||
    !nzchar(output_id)) {
  stop("Invalid exploration arguments.")
}

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  output_id
)
summary_dir <- file.path(output_dir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

fit_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "fash_fit1_update.RData"
)
message("Loading BF-adjusted real-data FASH fit.")
fash_fit <- load_bf_adjusted_real_fash(fit_path)
pair_keys <- names(fash_fit$fash_data$data_list)
gene_ids <- parse_gene_ids(pair_keys)
gene_index <- make_gene_index(pair_keys)
random_seeds <- random_seed_start + seq_len(random_reps) - 1L

message("Running the top-lfdr one-pair-per-gene design.")
top_selection <- select_highest_lfdr_per_gene(
  pair_keys,
  fash_fit$lfdr,
  top_n = top_n
)
top_data <- extract_z_matrix(fash_fit, top_selection$selected_indices)
if (!isTRUE(all.equal(top_data$time_grid, 0:15))) {
  stop("The real-data time grid is not 0:15.")
}

top_results <- list()
selection_counts <- list()
diagnostics <- list()
matrix_long <- list()
second_moment_long <- list()
adjacent_summary <- list()
lag_summary <- list()
output_index <- 1L

for (threshold in thresholds) {
  filtered <- filter_zero_intercept_z(top_data$z, threshold)
  estimate <- estimate_zero_intercept_correlation(
    top_data$z[filtered$keep, , drop = FALSE]
  )
  threshold_id <- paste0("threshold_", format(threshold, trim = TRUE))
  top_results[[threshold_id]] <- list(
    filter = filtered,
    estimate = estimate,
    z_nullish = top_data$z[filtered$keep, , drop = FALSE]
  )
  selection_counts[[length(selection_counts) + 1L]] <- data.frame(
    design = "Top-1000 highest lfdr per gene",
    threshold = threshold,
    seed = NA_integer_,
    n_candidates = filtered$n_candidates,
    n_selected = filtered$n_selected,
    selected_fraction = filtered$selected_fraction,
    stringsAsFactors = FALSE
  )
  diagnostics[[length(diagnostics) + 1L]] <- data.frame(
    design = "Top-1000 highest lfdr per gene",
    threshold = threshold,
    seed = NA_integer_,
    n_selected = estimate$n_selected,
    maximum_absolute_column_mean = max(abs(estimate$column_means)),
    minimum_second_moment_diagonal = min(estimate$second_moment_diagonal),
    maximum_second_moment_diagonal = max(estimate$second_moment_diagonal),
    maximum_correlation_difference = estimate$maximum_correlation_difference,
    sample_minimum_eigenvalue = estimate$sample_projection$diagnostics$raw_minimum_eigenvalue,
    sample_projection_maximum_change = estimate$sample_projection$diagnostics$maximum_absolute_change,
    normalized_minimum_eigenvalue = estimate$normalized_projection$diagnostics$raw_minimum_eigenvalue,
    normalized_projection_maximum_change = estimate$normalized_projection$diagnostics$maximum_absolute_change,
    stringsAsFactors = FALSE
  )
  second_moment_long[[length(second_moment_long) + 1L]] <-
    finite_square_matrix_to_long(
      estimate$second_moment,
      design = "Top-1000 highest lfdr per gene",
      threshold = threshold,
      estimator = "literal nullish-z second moment"
    )
  for (estimator_name in c("mashr cor(z)", "normalized second moment")) {
    correlation <- if (estimator_name == "mashr cor(z)") {
      estimate$sample_correlation
    } else {
      estimate$normalized_second_moment
    }
    matrix_long[[length(matrix_long) + 1L]] <- correlation_to_long(
      correlation,
      design = "Top-1000 highest lfdr per gene",
      threshold = threshold,
      estimator = estimator_name
    )
    adjacent_summary[[length(adjacent_summary) + 1L]] <-
      summarize_adjacent_correlations(
        correlation,
        design = "Top-1000 highest lfdr per gene",
        threshold = threshold,
        estimator = estimator_name
      )
    lag_summary[[length(lag_summary) + 1L]] <- summarize_correlation_lags(
      correlation,
      design = "Top-1000 highest lfdr per gene",
      threshold = threshold,
      estimator = estimator_name
    )
  }
  output_index <- output_index + 1L
}

message("Running pre-z random-one-variant-per-gene selections.")
random_results <- vector("list", length(random_seeds))
names(random_results) <- as.character(random_seeds)
random_stability <- list()
random_lag_stability <- list()

for (seed_index in seq_along(random_seeds)) {
  seed <- random_seeds[seed_index]
  if (seed_index == 1L || seed_index %% 5L == 0L ||
      seed_index == length(random_seeds)) {
    message(
      "Random selection ", seed_index, "/", length(random_seeds),
      " (seed ", seed, ")."
    )
  }
  selected <- select_random_variant_per_gene(
    pair_keys,
    seed = seed,
    gene_index = gene_index
  )
  extracted <- extract_z_matrix(fash_fit, selected$fash_index)
  seed_results <- list(selection = selected, thresholds = list())
  for (threshold in thresholds) {
    filtered <- filter_zero_intercept_z(extracted$z, threshold)
    estimate <- estimate_zero_intercept_correlation(
      extracted$z[filtered$keep, , drop = FALSE]
    )
    threshold_id <- paste0("threshold_", format(threshold, trim = TRUE))
    seed_results$thresholds[[threshold_id]] <- list(
      filter = filtered,
      sample_correlation = estimate$sample_correlation,
      second_moment = estimate$second_moment,
      normalized_second_moment = estimate$normalized_second_moment,
      column_means = estimate$column_means,
      second_moment_diagonal = estimate$second_moment_diagonal
    )
    selection_counts[[length(selection_counts) + 1L]] <- data.frame(
      design = "Random variant per gene before z",
      threshold = threshold,
      seed = seed,
      n_candidates = filtered$n_candidates,
      n_selected = filtered$n_selected,
      selected_fraction = filtered$selected_fraction,
      stringsAsFactors = FALSE
    )
    diagnostics[[length(diagnostics) + 1L]] <- data.frame(
      design = "Random variant per gene before z",
      threshold = threshold,
      seed = seed,
      n_selected = estimate$n_selected,
      maximum_absolute_column_mean = max(abs(estimate$column_means)),
      minimum_second_moment_diagonal = min(estimate$second_moment_diagonal),
      maximum_second_moment_diagonal = max(estimate$second_moment_diagonal),
      maximum_correlation_difference = estimate$maximum_correlation_difference,
      sample_minimum_eigenvalue = estimate$sample_projection$diagnostics$raw_minimum_eigenvalue,
      sample_projection_maximum_change = estimate$sample_projection$diagnostics$maximum_absolute_change,
      normalized_minimum_eigenvalue = estimate$normalized_projection$diagnostics$raw_minimum_eigenvalue,
      normalized_projection_maximum_change = estimate$normalized_projection$diagnostics$maximum_absolute_change,
      stringsAsFactors = FALSE
    )
    sample_correlation <- estimate$sample_correlation
    normalized_correlation <- estimate$normalized_second_moment
    random_stability[[length(random_stability) + 1L]] <- data.frame(
      seed = seed,
      threshold = threshold,
      n_selected = filtered$n_selected,
      selected_fraction = filtered$selected_fraction,
      day0_day1_correlation = sample_correlation[1, 2],
      mean_lag1_correlation = mean(sample_correlation[cbind(1:15, 2:16)]),
      normalized_day0_day1_correlation = normalized_correlation[1, 2],
      normalized_mean_lag1_correlation = mean(normalized_correlation[cbind(
        1:15,
        2:16
      )]),
      maximum_absolute_column_mean = max(abs(estimate$column_means)),
      minimum_second_moment_diagonal = min(estimate$second_moment_diagonal),
      maximum_second_moment_diagonal = max(estimate$second_moment_diagonal),
      stringsAsFactors = FALSE
    )
    random_lag_stability[[length(random_lag_stability) + 1L]] <-
      summarize_correlation_lags(
        sample_correlation,
        design = "Random variant per gene before z",
        threshold = threshold,
        estimator = "mashr cor(z)",
        seed = seed
      )
  }
  random_results[[seed_index]] <- seed_results
  rm(extracted)
  invisible(gc(FALSE))
}

selection_counts <- do.call(rbind, selection_counts)
diagnostics <- do.call(rbind, diagnostics)
random_stability <- do.call(rbind, random_stability)
random_lag_stability <- do.call(rbind, random_lag_stability)

message("Summarizing random-selection matrices.")
random_mean_results <- list()
random_mean_design <- paste0(
  "Mean across ",
  random_reps,
  " pre-z random selections"
)
for (threshold in thresholds) {
  threshold_id <- paste0("threshold_", format(threshold, trim = TRUE))
  sample_array <- simplify2array(lapply(
    random_results,
    function(x) x$thresholds[[threshold_id]]$sample_correlation
  ))
  normalized_array <- simplify2array(lapply(
    random_results,
    function(x) x$thresholds[[threshold_id]]$normalized_second_moment
  ))
  second_moment_array <- simplify2array(lapply(
    random_results,
    function(x) x$thresholds[[threshold_id]]$second_moment
  ))
  mean_sample <- apply(sample_array, c(1L, 2L), mean)
  mean_normalized <- apply(normalized_array, c(1L, 2L), mean)
  mean_second_moment <- apply(second_moment_array, c(1L, 2L), mean)
  random_mean_results[[threshold_id]] <- list(
    sample_correlation = mean_sample,
    second_moment = mean_second_moment,
    normalized_second_moment = mean_normalized,
    sample_array = sample_array,
    second_moment_array = second_moment_array,
    normalized_array = normalized_array
  )
  second_moment_long[[length(second_moment_long) + 1L]] <-
    finite_square_matrix_to_long(
      mean_second_moment,
      design = random_mean_design,
      threshold = threshold,
      estimator = "literal nullish-z second moment"
    )
  for (estimator_name in c("mashr cor(z)", "normalized second moment")) {
    correlation <- if (estimator_name == "mashr cor(z)") {
      mean_sample
    } else {
      mean_normalized
    }
    matrix_long[[length(matrix_long) + 1L]] <- correlation_to_long(
      correlation,
      design = random_mean_design,
      threshold = threshold,
      estimator = estimator_name
    )
    adjacent_summary[[length(adjacent_summary) + 1L]] <-
      summarize_adjacent_correlations(
        correlation,
        design = random_mean_design,
        threshold = threshold,
        estimator = estimator_name
      )
    lag_summary[[length(lag_summary) + 1L]] <- summarize_correlation_lags(
      correlation,
      design = random_mean_design,
      threshold = threshold,
      estimator = estimator_name
    )
  }
}

matrix_long <- do.call(rbind, matrix_long)
second_moment_long <- do.call(rbind, second_moment_long)
adjacent_summary <- do.call(rbind, adjacent_summary)
lag_summary <- do.call(rbind, lag_summary)

message("Running gene bootstraps for the threshold-2 primary comparisons.")
if (!2 %in% thresholds) {
  stop("The threshold list must contain 2 for the primary bootstrap.")
}
top_threshold2 <- top_results[["threshold_2"]]
top_bootstrap <- bootstrap_zero_intercept_correlations(
  top_threshold2$z_nullish,
  n_bootstrap = bootstrap_reps,
  seed = bootstrap_seed
)
top_bootstrap$design <- "Top-1000 highest lfdr per gene"

primary_random <- random_results[[1L]]
primary_random_threshold2 <- primary_random$thresholds[["threshold_2"]]
primary_random_z <- extract_z_matrix(
  fash_fit,
  primary_random$selection$fash_index
)$z
primary_random_z <- primary_random_z[
  primary_random_threshold2$filter$keep,
  ,
  drop = FALSE
]
random_bootstrap <- bootstrap_zero_intercept_correlations(
  primary_random_z,
  n_bootstrap = bootstrap_reps,
  seed = bootstrap_seed + 1L
)
random_bootstrap$design <- paste0(
  "Random variant per gene before z, seed ",
  random_seeds[1]
)
bootstrap_intervals <- rbind(top_bootstrap, random_bootstrap)

message("Running independent rectangular-truncation benchmarks.")
benchmark_candidates <- c(
  `Top-1000 candidate count` = top_n,
  `One-candidate-per-gene count` = length(gene_index)
)
benchmark_rows <- list()
set.seed(benchmark_seed)
for (candidate_name in names(benchmark_candidates)) {
  n_candidates <- benchmark_candidates[[candidate_name]]
  day01 <- matrix(NA_real_, nrow = benchmark_reps, ncol = length(thresholds))
  lag1 <- matrix(NA_real_, nrow = benchmark_reps, ncol = length(thresholds))
  selected_n <- matrix(NA_real_, nrow = benchmark_reps, ncol = length(thresholds))
  for (replication in seq_len(benchmark_reps)) {
    z <- matrix(
      stats::rnorm(n_candidates * 16L),
      nrow = n_candidates,
      ncol = 16L
    )
    for (threshold_index in seq_along(thresholds)) {
      threshold <- thresholds[threshold_index]
      filtered <- filter_zero_intercept_z(z, threshold)
      correlation <- stats::cor(z[filtered$keep, , drop = FALSE])
      selected_n[replication, threshold_index] <- filtered$n_selected
      day01[replication, threshold_index] <- correlation[1, 2]
      lag1[replication, threshold_index] <- mean(correlation[cbind(
        1:15,
        2:16
      )])
    }
  }
  for (threshold_index in seq_along(thresholds)) {
    for (statistic in c("n_selected", "day0_day1", "mean_lag1")) {
      values <- switch(
        statistic,
        n_selected = selected_n[, threshold_index],
        day0_day1 = day01[, threshold_index],
        mean_lag1 = lag1[, threshold_index]
      )
      interval <- summarize_interval(values)
      benchmark_rows[[length(benchmark_rows) + 1L]] <- data.frame(
        candidate_design = candidate_name,
        n_candidates = n_candidates,
        threshold = thresholds[threshold_index],
        statistic = statistic,
        mean = interval["mean"],
        median = interval["median"],
        lower = interval["lower"],
        upper = interval["upper"],
        n_replications = benchmark_reps,
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    }
  }
}
independence_benchmark <- do.call(rbind, benchmark_rows)

message("Calibrating correlation attenuation caused by z truncation.")
rho_grid <- seq(-0.3, 0.3, by = 0.1)
calibration_rows <- list()
for (rho_index in seq_along(rho_grid)) {
  rho <- rho_grid[rho_index]
  target <- make_lag1_only_correlation(16L, rho)
  set.seed(benchmark_seed + 1000L + rho_index)
  independent <- matrix(
    stats::rnorm(calibration_n * 16L),
    nrow = calibration_n,
    ncol = 16L
  )
  z <- independent %*% chol(target)
  for (threshold in thresholds) {
    filtered <- filter_zero_intercept_z(z, threshold)
    correlation <- stats::cor(z[filtered$keep, , drop = FALSE])
    calibration_rows[[length(calibration_rows) + 1L]] <- data.frame(
      generating_rho = rho,
      threshold = threshold,
      n_candidates = calibration_n,
      n_selected = filtered$n_selected,
      selected_fraction = filtered$selected_fraction,
      day0_day1_correlation = correlation[1, 2],
      mean_lag1_correlation = mean(correlation[cbind(1:15, 2:16)]),
      stringsAsFactors = FALSE
    )
  }
}
truncation_calibration <- do.call(rbind, calibration_rows)

fit_info <- file.info(fit_path)
configuration <- list(
  analysis_id = output_id,
  created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  fit_path = normalizePath(fit_path),
  fit_mtime = as.character(fit_info$mtime),
  fit_size = unname(fit_info$size),
  fit_md5 = unname(tools::md5sum(fit_path)),
  n_fash_pairs = length(pair_keys),
  n_genes = length(gene_index),
  top_n = top_n,
  thresholds = thresholds,
  random_seeds = random_seeds,
  bootstrap_reps = bootstrap_reps,
  bootstrap_seed = bootstrap_seed,
  benchmark_reps = benchmark_reps,
  benchmark_seed = benchmark_seed,
  calibration_n = calibration_n,
  rho_grid = rho_grid,
  primary_estimator = "stats::cor(z[apply(abs(z) < threshold, 1L, all), ])",
  second_moment_estimator = "crossprod(z_nullish) / nrow(z_nullish)",
  projection = "Matrix::nearPD(corr = TRUE, keepDiag = TRUE)",
  r_version = R.version.string,
  mashr_version = if (requireNamespace("mashr", quietly = TRUE)) {
    as.character(utils::packageVersion("mashr"))
  } else {
    NA_character_
  },
  matrix_version = as.character(utils::packageVersion("Matrix"))
)

top_selected_metadata <- top_selection$selected
top_selected_metadata$max_absolute_z <- apply(abs(top_data$z), 1L, max)
top_selected_metadata$passes_threshold_2 <-
  top_selected_metadata$max_absolute_z < 2
primary_random_metadata <- primary_random$selection
primary_random_metadata$max_absolute_z <- apply(
  abs(extract_z_matrix(fash_fit, primary_random$selection$fash_index)$z),
  1L,
  max
)
primary_random_metadata$passes_threshold_2 <-
  primary_random_metadata$max_absolute_z < 2

analysis <- list(
  configuration = configuration,
  top_selection = top_selected_metadata,
  primary_random_selection = primary_random_metadata,
  top_results = top_results,
  random_results = random_results,
  random_mean_results = random_mean_results,
  selection_counts = selection_counts,
  estimator_diagnostics = diagnostics,
  correlation_matrices_long = matrix_long,
  second_moments_long = second_moment_long,
  adjacent_correlations = adjacent_summary,
  lag_summaries = lag_summary,
  bootstrap_intervals = bootstrap_intervals,
  random_variant_stability = random_stability,
  random_lag_stability = random_lag_stability,
  independence_benchmark = independence_benchmark,
  truncation_calibration = truncation_calibration
)

saveRDS(configuration, file.path(output_dir, "configuration.rds"))
saveRDS(analysis, file.path(output_dir, "zero_intercept_correlation_analysis.rds"))
write_csv(selection_counts, file.path(summary_dir, "selection_counts.csv"))
write_csv(diagnostics, file.path(summary_dir, "estimator_diagnostics.csv"))
write_csv(matrix_long, file.path(summary_dir, "correlation_matrices_long.csv"))
write_csv(second_moment_long, file.path(summary_dir, "second_moments_long.csv"))
write_csv(adjacent_summary, file.path(summary_dir, "adjacent_correlations.csv"))
write_csv(lag_summary, file.path(summary_dir, "lag_summaries.csv"))
write_csv(bootstrap_intervals, file.path(summary_dir, "bootstrap_intervals.csv"))
write_csv(random_stability, file.path(summary_dir, "random_variant_stability.csv"))
write_csv(random_lag_stability, file.path(
  summary_dir,
  "random_lag_stability.csv"
))
write_csv(independence_benchmark, file.path(summary_dir, "independence_benchmark.csv"))
write_csv(truncation_calibration, file.path(summary_dir, "truncation_calibration.csv"))
write_csv(top_selected_metadata, file.path(summary_dir, "top_lfdr_selected_units.csv"))
write_csv(primary_random_metadata, file.path(
  summary_dir,
  "primary_random_selected_units.csv"
))

message("Zero-intercept correlation exploration completed: ", output_dir)
message(
  "Top-lfdr threshold-2 retained ",
  top_results$threshold_2$filter$n_selected,
  "/", top_n,
  "; mean lag-1 correlation = ",
  format(
    mean(top_results$threshold_2$estimate$sample_correlation[cbind(
      1:15,
      2:16
    )]),
    digits = 5
  ),
  "."
)
random_threshold2 <- random_stability[random_stability$threshold == 2, ]
message(
  "Random-pre-z threshold-2 mean lag-1 correlation across selections = ",
  format(mean(random_threshold2$mean_lag1_correlation), digits = 5),
  " [range ",
  format(min(random_threshold2$mean_lag1_correlation), digits = 5),
  ", ",
  format(max(random_threshold2$mean_lag1_correlation), digits = 5),
  "]."
)
