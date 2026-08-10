# Load, validate, and format cached results for the R4 workflowr report.

source("code/revision_simulations/r4_correlated_errors/plotting.R")

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required R4 cache file is missing: ", path)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

expected_seeds <- c(12345L, 22345L, 32345L, 42345L, 52345L)
expected_rho <- seq(-0.3, 0.3, by = 0.1)
expected_alpha_methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "FASH-linear-Raw",
  "FASH-linear-BF"
)
expected_iwp_methods <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
expected_class_probs <- c(
  dynamic_bspline = 0.20,
  constant = 0.40,
  zero = 0.40
)
expected_full_conditions <- c(
  "Independent",
  "Direct centered full matrix",
  "Pairwise-difference full matrix"
)

# Real-data full-correlation analysis.
real_output_dir <- file.path(
  "output",
  "revision_simulations",
  "real_data",
  "r4_null_like_top500_full_correlations"
)
real_analysis_path <- file.path(
  real_output_dir,
  "real_data_correlation_analysis.rds"
)
if (!file.exists(real_analysis_path)) {
  stop("The real-data full-correlation cache is missing.")
}
real_analysis <- readRDS(real_analysis_path)
real_configuration <- real_analysis$configuration

if (!identical(real_configuration$fit_object, "fash_fit1_update") ||
    !isTRUE(all.equal(real_configuration$top_n, 500L)) ||
    !isTRUE(all.equal(real_configuration$n_gene_representatives, 6362L)) ||
    !isTRUE(all.equal(real_configuration$time_grid, 0:15)) ||
    !isTRUE(all.equal(real_configuration$n_bootstrap, 1000L)) ||
    !isTRUE(all.equal(real_configuration$n_independence_benchmark, 1000L)) ||
    nrow(real_analysis$selected_units) != 500L ||
    anyDuplicated(real_analysis$selected_units$gene_id) ||
    anyDuplicated(real_analysis$selected_units$pair_key) ||
    nrow(real_analysis$projection_diagnostics) != 4L ||
    nrow(real_analysis$bootstrap_lags) != 30L ||
    nrow(real_analysis$independence_benchmark) != 60L ||
    any(real_analysis$bootstrap_lags$n_replications != 1000L) ||
    any(real_analysis$independence_benchmark$n_replications != 1000L) ||
    any(!real_analysis$projection_diagnostics$converged)) {
  stop("The real-data full-correlation cache is incomplete or mismatched.")
}

primary_matrix_ids <- c("direct_adjusted", "pairwise_adjusted")
if (!all(primary_matrix_ids %in% names(real_analysis$matrices))) {
  stop("The real-data cache does not contain both primary full matrices.")
}
for (matrix_id in primary_matrix_ids) {
  result <- real_analysis$matrices[[matrix_id]]
  if (!identical(dim(result$raw), c(16L, 16L)) ||
      !identical(dim(result$projected), c(16L, 16L)) ||
      any(!is.finite(result$raw)) || any(!is.finite(result$projected)) ||
      max(abs(diag(result$projected) - 1)) > 1e-10 ||
      min(eigen(
        result$projected,
        symmetric = TRUE,
        only.values = TRUE
      )$values) <= 0) {
    stop("A projected real-data full matrix is invalid: ", matrix_id)
  }
}

# Paired R1 simulations using the two full matrices.
full_output_dir <- file.path(
  "output",
  "revision_simulations",
  "mc",
  "r4_full_empirical_correlations_top500_pilot5"
)
full_summary_dir <- file.path(full_output_dir, "summary")
full_configuration <- readRDS(file.path(full_output_dir, "configuration.rds"))
full_condition_alpha <- read_required_csv(file.path(
  full_summary_dir,
  "condition_mc_alpha005_summary.csv"
))
full_paired_alpha <- read_required_csv(file.path(
  full_summary_dir,
  "paired_alpha005_differences.csv"
))
full_condition_pi0 <- read_required_csv(file.path(
  full_summary_dir,
  "condition_mc_pi0_summary.csv"
))
full_paired_pi0 <- read_required_csv(file.path(
  full_summary_dir,
  "paired_pi0_differences.csv"
))
full_matrix_summary <- read_required_csv(file.path(
  full_summary_dir,
  "realized_correlation_matrix_summary.csv"
))
full_lag_summary <- read_required_csv(file.path(
  full_summary_dir,
  "realized_lag_summary.csv"
))
full_pairing_check <- read_required_csv(file.path(
  full_summary_dir,
  "pairing_check.csv"
))
full_target_checks <- read_required_csv(file.path(
  full_summary_dir,
  "target_matrix_checks.csv"
))
full_r1_reference_check <- read_required_csv(file.path(
  full_summary_dir,
  "r1_reference_check.csv"
))

if (!identical(
      full_configuration$scenario,
      "genotype_random_bspline_main_effect_dynamic_eqtl"
    ) ||
    !isTRUE(all.equal(full_configuration$J, 1000L)) ||
    !isTRUE(all.equal(full_configuration$n_donors, 19L)) ||
    !isTRUE(all.equal(full_configuration$n_covariates, 5L)) ||
    !isTRUE(all.equal(full_configuration$time_grid, 0:15)) ||
    !isTRUE(all.equal(full_configuration$class_probs, expected_class_probs)) ||
    !isTRUE(all.equal(full_configuration$expression_noise_sd, 1)) ||
    !isTRUE(all.equal(full_configuration$dynamic_main_effect_sd, 1)) ||
    !isTRUE(all.equal(full_configuration$num_basis, 20L)) ||
    !isTRUE(all.equal(full_configuration$true_pi0, 0.8)) ||
    !identical(full_configuration$seed_list, expected_seeds) ||
    !identical(sort(full_configuration$conditions),
               sort(expected_full_conditions)) ||
    !identical(sort(full_configuration$methods), sort(expected_iwp_methods)) ||
    nrow(full_condition_alpha) != 6L ||
    nrow(full_paired_alpha) != 4L ||
    nrow(full_condition_pi0) != 6L ||
    nrow(full_paired_pi0) != 4L ||
    nrow(full_matrix_summary) != 2304L ||
    nrow(full_lag_summary) != 135L ||
    nrow(full_pairing_check) != 15L ||
    nrow(full_target_checks) != 3L ||
    nrow(full_r1_reference_check) != 2L ||
    any(full_condition_alpha$n_replications != 5L) ||
    any(full_condition_pi0$n_replications != 5L) ||
    any(full_lag_summary$n_replications != 5L) ||
    any(!full_pairing_check$passed) ||
    any(full_pairing_check$maximum_absolute_error_difference >
        full_pairing_check$tolerance) ||
    any(!full_target_checks$positive_definite) ||
    any(full_target_checks$minimum_eigenvalue <= 0) ||
    any(!full_r1_reference_check$passed) ||
    any(full_r1_reference_check$maximum_absolute_difference >
        full_r1_reference_check$tolerance)) {
  stop("The paired full-matrix R1 cache is incomplete or mismatched.")
}

# Lag-1 strength sweep with IWP1 pi0 summaries.
sweep_output_dir <- file.path(
  "output",
  "revision_simulations",
  "mc",
  "r4_lag1_correlation_sweep_m0p3_to_p0p3_pi0_pilot5"
)
sweep_summary_dir <- file.path(sweep_output_dir, "summary")
sweep_configuration <- readRDS(file.path(sweep_output_dir, "configuration.rds"))
sweep_alpha_005 <- read_required_csv(file.path(
  sweep_summary_dir,
  "mc_alpha005_summary.csv"
))
sweep_paired_005 <- read_required_csv(file.path(
  sweep_summary_dir,
  "paired_vs_zero_alpha005_summary.csv"
))
sweep_pi0 <- read_required_csv(file.path(
  sweep_summary_dir,
  "mc_pi0_summary.csv"
))
sweep_paired_pi0 <- read_required_csv(file.path(
  sweep_summary_dir,
  "paired_vs_zero_pi0_summary.csv"
))
sweep_lag1 <- read_required_csv(file.path(
  sweep_summary_dir,
  "lag1_correlation_summary.csv"
))
sweep_pairing_check <- read_required_csv(file.path(
  sweep_summary_dir,
  "pairing_check.csv"
))
sweep_matrix_check <- read_required_csv(file.path(
  sweep_summary_dir,
  "correlation_matrix_check.csv"
))
sweep_r1_reference_check <- read_required_csv(file.path(
  sweep_summary_dir,
  "r1_reference_check.csv"
))
sweep_r1_pi0_reference_check <- read_required_csv(file.path(
  sweep_summary_dir,
  "r1_pi0_reference_check.csv"
))

if (!identical(
      sweep_configuration$scenario,
      "genotype_random_bspline_main_effect_dynamic_eqtl"
    ) ||
    !isTRUE(all.equal(sweep_configuration$J, 1000L)) ||
    !isTRUE(all.equal(sweep_configuration$n_donors, 19L)) ||
    !isTRUE(all.equal(sweep_configuration$n_covariates, 5L)) ||
    !isTRUE(all.equal(sweep_configuration$time_grid, 0:15)) ||
    !isTRUE(all.equal(sweep_configuration$class_probs, expected_class_probs)) ||
    !isTRUE(all.equal(sweep_configuration$expression_noise_sd, 1)) ||
    !isTRUE(all.equal(sweep_configuration$dynamic_main_effect_sd, 1)) ||
    !isTRUE(all.equal(sweep_configuration$num_basis, 20L)) ||
    !isTRUE(all.equal(sweep_configuration$true_pi0, 0.8)) ||
    !identical(sweep_configuration$seed_list, expected_seeds) ||
    !isTRUE(all.equal(sweep_configuration$rho_grid, expected_rho)) ||
    !identical(sweep_configuration$correlation_structure, "lag1_only") ||
    !identical(sort(sweep_configuration$pi0_methods),
               sort(expected_iwp_methods)) ||
    nrow(sweep_alpha_005) != 28L ||
    nrow(sweep_paired_005) != 28L ||
    nrow(sweep_pi0) != 14L ||
    nrow(sweep_paired_pi0) != 14L ||
    nrow(sweep_lag1) != 14L ||
    nrow(sweep_pairing_check) != 35L ||
    nrow(sweep_matrix_check) != 7L ||
    !identical(sort(unique(sweep_alpha_005$method)),
               sort(expected_alpha_methods)) ||
    !identical(sort(unique(sweep_pi0$method)), sort(expected_iwp_methods)) ||
    !isTRUE(all.equal(sort(unique(sweep_alpha_005$rho)), expected_rho)) ||
    !isTRUE(all.equal(sort(unique(sweep_pi0$rho)), expected_rho)) ||
    any(sweep_alpha_005$n_replications != 5L) ||
    any(sweep_pi0$n_replications != 5L) ||
    any(!sweep_pairing_check$passed) ||
    any(sweep_pairing_check$maximum_absolute_error_difference >
        sweep_pairing_check$tolerance) ||
    any(!sweep_matrix_check$positive_definite) ||
    any(sweep_matrix_check$minimum_eigenvalue <= 0) ||
    any(!sweep_r1_reference_check$passed) ||
    any(!sweep_r1_pi0_reference_check$passed) ||
    any(sweep_r1_reference_check$maximum_absolute_difference >
        sweep_r1_reference_check$tolerance) ||
    any(sweep_r1_pi0_reference_check$maximum_absolute_difference >
        sweep_r1_pi0_reference_check$tolerance)) {
  stop("The lag-1 correlation-sweep cache is incomplete or mismatched.")
}

format_decimal <- function(x, digits = 3) {
  formatC(x, format = "f", digits = digits)
}

format_interval <- function(mean, lower, upper, digits = 3) {
  paste0(
    format_decimal(mean, digits),
    " [",
    format_decimal(lower, digits),
    ", ",
    format_decimal(upper, digits),
    "]"
  )
}

format_signed_interval <- function(mean, lower, upper, digits = 3) {
  paste0(
    ifelse(mean >= 0, "+", ""),
    format_decimal(mean, digits),
    " [",
    ifelse(lower >= 0, "+", ""),
    format_decimal(lower, digits),
    ", ",
    ifelse(upper >= 0, "+", ""),
    format_decimal(upper, digits),
    "]"
  )
}

format_rho <- function(x, digits = 1) {
  value <- formatC(x, format = "f", digits = digits)
  ifelse(x > 0, paste0("+", value), value)
}

render_scrollable_table <- function(x,
                                    caption = NULL,
                                    align = NULL,
                                    minimum_width = "720px") {
  table_html <- knitr::kable(
    x,
    format = "html",
    caption = caption,
    align = align,
    table.attr = paste0(
      "class=\"table table-condensed\" style=\"min-width:",
      minimum_width,
      ";\""
    )
  )
  knitr::asis_output(paste0(
    "<div class=\"r4-table-scroll\">\n",
    table_html,
    "\n</div>"
  ))
}

projection_table <- do.call(rbind, lapply(primary_matrix_ids, function(matrix_id) {
  result <- real_analysis$matrices[[matrix_id]]
  diagnostic <- real_analysis$projection_diagnostics[
    real_analysis$projection_diagnostics$matrix_id == matrix_id,
    ,
    drop = FALSE
  ]
  raw_lag1 <- mean(result$raw[cbind(1:15, 2:16)])
  projected_lag1 <- mean(result$projected[cbind(1:15, 2:16)])
  data.frame(
    Estimator = result$estimator,
    `Raw minimum eigenvalue` = format_decimal(
      diagnostic$raw_minimum_eigenvalue,
      4
    ),
    `Projected minimum eigenvalue` = format_decimal(
      diagnostic$projected_minimum_eigenvalue,
      7
    ),
    `Maximum entry change` = format_decimal(
      diagnostic$maximum_absolute_change,
      6
    ),
    `Frobenius change` = format_decimal(diagnostic$frobenius_change, 5),
    `Raw lag-1 correlation` = format_decimal(raw_lag1),
    `Projected lag-1 correlation` = format_decimal(projected_lag1),
    check.names = FALSE
  )
}))
rownames(projection_table) <- NULL

se_sensitivity_rows <- real_analysis$lag_summaries[
  real_analysis$lag_summaries$matrix_version == "Raw estimate" &
    real_analysis$lag_summaries$lag == 1L,
  ,
  drop = FALSE
]
se_sensitivity_table <- data.frame(
  Estimator = se_sensitivity_rows$estimator,
  `SE scale` = se_sensitivity_rows$se_scale,
  `Lag-1 correlation` = format_decimal(se_sensitivity_rows$mean_correlation),
  `Lag-1 semivariogram` = format_decimal(se_sensitivity_rows$semivariogram),
  check.names = FALSE
)

observed_lag1 <- real_analysis$bootstrap_lags[
  real_analysis$bootstrap_lags$lag == 1L,
  ,
  drop = FALSE
]
fixed_benchmark_lag1 <- real_analysis$independence_benchmark[
  real_analysis$independence_benchmark$lag == 1L &
    real_analysis$independence_benchmark$benchmark ==
      "Fixed selected SE patterns",
  ,
  drop = FALSE
]
flat_benchmark_lag1 <- real_analysis$independence_benchmark[
  real_analysis$independence_benchmark$lag == 1L &
    real_analysis$independence_benchmark$benchmark ==
      "Flatness-selected gene representatives",
  ,
  drop = FALSE
]
benchmark_lag1 <- merge(
  observed_lag1,
  fixed_benchmark_lag1,
  by = "estimator",
  suffixes = c("_observed", "_fixed"),
  sort = FALSE
)
names(benchmark_lag1)[names(benchmark_lag1) == "mean_correlation"] <-
  "mean_correlation_fixed"
benchmark_lag1 <- merge(
  benchmark_lag1,
  flat_benchmark_lag1,
  by = "estimator",
  sort = FALSE
)
benchmark_lag1_table <- data.frame(
  Estimator = benchmark_lag1$estimator,
  `Observed correlation (95% gene-bootstrap CI)` = format_interval(
    benchmark_lag1$mean_correlation_observed,
    benchmark_lag1$correlation_ci_lower_observed,
    benchmark_lag1$correlation_ci_upper_observed
  ),
  `Independent, fixed SE (95% interval)` = format_interval(
    benchmark_lag1$mean_correlation_fixed,
    benchmark_lag1$correlation_ci_lower_fixed,
    benchmark_lag1$correlation_ci_upper_fixed
  ),
  `Independent, flatness-selected (95% interval)` = format_interval(
    benchmark_lag1$mean_correlation,
    benchmark_lag1$correlation_ci_lower,
    benchmark_lag1$correlation_ci_upper
  ),
  check.names = FALSE
)

full_results <- merge(
  full_condition_alpha,
  full_condition_pi0,
  by = c("condition", "method", "n_replications"),
  sort = FALSE
)
full_results <- merge(
  full_results,
  full_paired_alpha[, c(
    "condition", "method", "mean_delta_power", "power_delta_ci_lower",
    "power_delta_ci_upper", "mean_delta_fdr", "fdr_delta_ci_lower",
    "fdr_delta_ci_upper"
  )],
  by = c("condition", "method"),
  all.x = TRUE,
  sort = FALSE
)
full_results <- merge(
  full_results,
  full_paired_pi0[, c(
    "condition", "method", "mean_delta_pi0", "pi0_delta_ci_lower",
    "pi0_delta_ci_upper"
  )],
  by = c("condition", "method"),
  all.x = TRUE,
  sort = FALSE
)
full_results$condition <- factor(
  full_results$condition,
  levels = expected_full_conditions
)
full_results$method <- factor(full_results$method, levels = expected_iwp_methods)
full_results <- full_results[order(full_results$condition, full_results$method), ]
format_delta_or_reference <- function(mean, lower, upper) {
  ifelse(
    is.na(mean),
    "Reference",
    format_signed_interval(mean, lower, upper)
  )
}
full_results_table <- data.frame(
  Condition = as.character(full_results$condition),
  Fit = ifelse(full_results$method == "FASH-IWP1-Raw", "Raw", "BF-corrected"),
  `Power (95% MC CI)` = format_interval(
    full_results$mean_power,
    full_results$power_ci_lower,
    full_results$power_ci_upper
  ),
  `E[FDP] (95% MC CI)` = format_interval(
    full_results$mean_fdr,
    full_results$fdr_ci_lower,
    full_results$fdr_ci_upper
  ),
  `Estimated pi0 (95% MC CI)` = format_interval(
    full_results$mean_estimated_pi0,
    full_results$pi0_ci_lower,
    full_results$pi0_ci_upper
  ),
  `Power change vs independent` = format_delta_or_reference(
    full_results$mean_delta_power,
    full_results$power_delta_ci_lower,
    full_results$power_delta_ci_upper
  ),
  `FDR change vs independent` = format_delta_or_reference(
    full_results$mean_delta_fdr,
    full_results$fdr_delta_ci_lower,
    full_results$fdr_delta_ci_upper
  ),
  `pi0 change vs independent` = format_delta_or_reference(
    full_results$mean_delta_pi0,
    full_results$pi0_delta_ci_lower,
    full_results$pi0_delta_ci_upper
  ),
  check.names = FALSE
)

full_lag1 <- full_lag_summary[full_lag_summary$lag == 1L, , drop = FALSE]
lag1_wide <- reshape(
  full_lag1[, c(
    "condition", "diagnostic", "mean_correlation",
    "correlation_ci_lower", "correlation_ci_upper"
  )],
  idvar = "condition",
  timevar = "diagnostic",
  direction = "wide"
)
lag1_wide <- merge(
  full_target_checks[, c("condition", "mean_lag1_correlation")],
  lag1_wide,
  by = "condition",
  sort = FALSE
)
lag1_wide$condition <- factor(lag1_wide$condition, levels = expected_full_conditions)
lag1_wide <- lag1_wide[order(lag1_wide$condition), ]
lag_interval_from_wide <- function(data, diagnostic) {
  format_interval(
    data[[paste0("mean_correlation.", diagnostic)]],
    data[[paste0("correlation_ci_lower.", diagnostic)]],
    data[[paste0("correlation_ci_upper.", diagnostic)]]
  )
}
full_lag1_table <- data.frame(
  Condition = as.character(lag1_wide$condition),
  `Target lag-1 correlation` = format_decimal(lag1_wide$mean_lag1_correlation),
  `Expression error (95% MC CI)` = lag_interval_from_wide(
    lag1_wide,
    "Expression error"
  ),
  `Truth-known beta-hat error (95% MC CI)` = lag_interval_from_wide(
    lag1_wide,
    "Truth-known standardized beta-hat error"
  ),
  `Centered dynamic-null residual (95% MC CI)` = lag_interval_from_wide(
    lag1_wide,
    "Centered dynamic-null beta-hat residual"
  ),
  check.names = FALSE
)

primary_sweep_method <- "FASH-IWP1-BF"
sweep_primary <- sweep_alpha_005[
  sweep_alpha_005$method == primary_sweep_method,
  ,
  drop = FALSE
]
sweep_primary <- merge(
  sweep_primary,
  sweep_paired_005[sweep_paired_005$method == primary_sweep_method, ],
  by = c("rho", "scenario", "target", "method", "alpha"),
  suffixes = c("", "_paired"),
  sort = FALSE
)
sweep_primary <- merge(
  sweep_primary,
  sweep_pi0[sweep_pi0$method == primary_sweep_method, ],
  by = c("rho", "method"),
  sort = FALSE
)
sweep_primary <- merge(
  sweep_primary,
  sweep_lag1[
    sweep_lag1$diagnostic == "Truth-known standardized error",
    ,
    drop = FALSE
  ],
  by = "rho",
  sort = FALSE
)
sweep_primary <- sweep_primary[order(sweep_primary$rho), ]
sweep_primary_table <- data.frame(
  `Generating rho` = format_rho(sweep_primary$rho),
  `Realized beta-hat correlation (95% MC CI)` = format_interval(
    sweep_primary$mean_lag1_correlation,
    sweep_primary$lag1_ci_lower,
    sweep_primary$lag1_ci_upper
  ),
  `Power (95% MC CI)` = format_interval(
    sweep_primary$mean_power,
    sweep_primary$power_ci_lower,
    sweep_primary$power_ci_upper
  ),
  `E[FDP] (95% MC CI)` = format_interval(
    sweep_primary$mean_fdr,
    sweep_primary$fdr_ci_lower,
    sweep_primary$fdr_ci_upper
  ),
  `Estimated pi0 (95% MC CI)` = format_interval(
    sweep_primary$mean_estimated_pi0,
    sweep_primary$pi0_ci_lower,
    sweep_primary$pi0_ci_upper
  ),
  `Power change vs rho=0` = format_signed_interval(
    sweep_primary$mean_delta_power,
    sweep_primary$power_delta_ci_lower,
    sweep_primary$power_delta_ci_upper
  ),
  `FDR change vs rho=0` = format_signed_interval(
    sweep_primary$mean_delta_fdr,
    sweep_primary$fdr_delta_ci_lower,
    sweep_primary$fdr_delta_ci_upper
  ),
  check.names = FALSE
)

get_real_lag_metric <- function(estimator, metric, lag = 1L) {
  row <- real_analysis$bootstrap_lags[
    real_analysis$bootstrap_lags$estimator == estimator &
      real_analysis$bootstrap_lags$lag == lag,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L || !metric %in% names(row)) {
    stop("Could not extract a unique real-data lag metric.")
  }
  row[[metric]]
}

get_benchmark_metric <- function(estimator, benchmark, metric, lag = 1L) {
  row <- real_analysis$independence_benchmark[
    real_analysis$independence_benchmark$estimator == estimator &
      real_analysis$independence_benchmark$benchmark == benchmark &
      real_analysis$independence_benchmark$lag == lag,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L || !metric %in% names(row)) {
    stop("Could not extract a unique independence-benchmark metric.")
  }
  row[[metric]]
}

get_full_metric <- function(condition, method, metric) {
  row <- full_results[
    as.character(full_results$condition) == condition &
      as.character(full_results$method) == method,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L || !metric %in% names(row)) {
    stop("Could not extract a unique full-matrix metric.")
  }
  row[[metric]]
}

get_sweep_metric <- function(rho,
                             metric,
                             method = "FASH-IWP1-BF") {
  row <- sweep_alpha_005[
    abs(sweep_alpha_005$rho - rho) < 1e-12 &
      sweep_alpha_005$method == method,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L || !metric %in% names(row)) {
    stop("Could not extract a unique sweep alpha metric.")
  }
  row[[metric]]
}

get_sweep_paired_metric <- function(rho,
                                    metric,
                                    method = "FASH-IWP1-BF") {
  row <- sweep_paired_005[
    abs(sweep_paired_005$rho - rho) < 1e-12 &
      sweep_paired_005$method == method,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L || !metric %in% names(row)) {
    stop("Could not extract a unique paired sweep metric.")
  }
  row[[metric]]
}

get_sweep_pi0_metric <- function(rho,
                                 metric = "mean_estimated_pi0",
                                 method = "FASH-IWP1-BF") {
  row <- sweep_pi0[
    abs(sweep_pi0$rho - rho) < 1e-12 & sweep_pi0$method == method,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L || !metric %in% names(row)) {
    stop("Could not extract a unique sweep pi0 metric.")
  }
  row[[metric]]
}
