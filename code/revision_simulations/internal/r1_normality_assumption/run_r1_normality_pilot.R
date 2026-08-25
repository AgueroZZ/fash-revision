#!/usr/bin/env Rscript

# Run the paired seed-12345 R1 finite-sample normality pilot.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return(normalizePath("coderepo-local", winslash = "/", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(args, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(args[equals_hit[1L]], nchar(equals_prefix) + 1L))
  }
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1L] == length(args)) {
    return(default)
  }
  args[hit[1L] + 1L]
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

capture_warnings <- function(expression) {
  warnings <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "r1_normality_assumption", "r1_normality_helpers.R"
))

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required.")
}

num_cores <- as.integer(get_arg("--num-cores", "8"))
t_df <- as.numeric(get_arg("--t-df", "5"))
output_id <- get_arg(
  "--output-id",
  "r1_normality_assumption_seed12345_t5_v1"
)
if (is.na(num_cores) || num_cores < 1L ||
    length(t_df) != 1L || !is.finite(t_df) || t_df <= 2 ||
    length(output_id) != 1L || !nzchar(output_id)) {
  stop("Invalid pilot arguments.")
}

source_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "mc",
  paste0(
    "r1_real_genotype_one_per_gene_J6362_random_bspline_main_effect_",
    "linear_mixture_predstep1_penalty10_pilot5"
  ),
  "full_fits", "seed_12345.rds"
)
if (!file.exists(source_path)) {
  stop("The formal R1 seed-12345 source object is missing: ", source_path)
}

output_root <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal"
)
final_directory <- file.path(output_root, output_id)
checkpoint_directory <- file.path(
  output_root,
  paste0(output_id, "_checkpoints")
)
staging_directory <- paste0(
  final_directory,
  ".staging_",
  Sys.getpid()
)
if (dir.exists(final_directory) || dir.exists(staging_directory)) {
  stop("The requested final or staging output directory already exists.")
}
dir.create(staging_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_directory, recursive = TRUE, showWarnings = FALSE)
completed <- FALSE
on.exit({
  if (!completed && dir.exists(staging_directory)) {
    unlink(staging_directory, recursive = TRUE, force = TRUE)
  }
}, add = TRUE)

message("Loading and validating the formal R1 seed-12345 source object.")
source_object <- readRDS(source_path)
required_source_fields <- c(
  "unit_info", "genotype", "covariates", "true_beta", "expression",
  "expression_simulation", "eqtl_summary", "fash_fits", "settings"
)
if (!is.list(source_object) ||
    !all(required_source_fields %in% names(source_object))) {
  stop("The formal R1 source object is incomplete.")
}
expected_classes <- c(
  constant = 2545L,
  dynamic_bspline = 1272L,
  zero = 2545L
)
observed_classes <- table(factor(
  source_object$unit_info$effect_class,
  levels = names(expected_classes)
))
source_fit_raw <- source_object$fash_fits$fash_iwp1_raw
source_fit_bf <- source_object$fash_fits$fash_iwp1_bf
if (!identical(dim(source_object$genotype), c(19L, 6362L)) ||
    !identical(dim(source_object$true_beta), c(6362L, 16L)) ||
    !identical(dim(source_object$expression), c(19L, 6362L, 16L)) ||
    !identical(dim(source_object$covariates), c(19L, 5L)) ||
    !identical(as.integer(observed_classes), unname(expected_classes)) ||
    !identical(as.integer(source_object$settings$seed), 12345L) ||
    !identical(as.integer(source_object$settings$n_donors), 19L) ||
    !identical(as.integer(source_object$settings$n_variants), 6362L) ||
    !identical(as.integer(source_object$settings$n_covariates), 5L) ||
    !all(as.integer(source_object$eqtl_summary$df) == 12L) ||
    !identical(as.integer(unique(source_object$eqtl_summary$rank)), 7L) ||
    length(source_fit_raw$psd_grid) != 52L ||
    !identical(as.numeric(source_fit_raw$settings$penalty), 10) ||
    !identical(as.numeric(source_fit_raw$settings$pred_step), 1) ||
    !identical(as.integer(source_fit_raw$settings$num_basis), 20L)) {
  stop("The formal R1 source object does not match the fixed pilot design.")
}

time_grid <- as.numeric(source_object$settings$time_grid)
unit_info <- source_object$unit_info
true_beta <- source_object$true_beta
noise_sd <- as.numeric(source_object$settings$expression_noise_sd)
expression_mean <- reconstruct_r1_expression_mean(
  G = source_object$genotype,
  beta_matrix = true_beta,
  covariates = source_object$covariates,
  covariate_effects = source_object$expression_simulation$covariate_effects,
  intercepts = source_object$expression_simulation$intercepts
)
gaussian_error <- source_object$expression - expression_mean
t_error <- gaussian_to_standardized_t(gaussian_error, df = t_df)
t_expression <- expression_mean + t_error
if (!identical(
      order(as.vector(gaussian_error)),
      order(as.vector(t_error))
    ) ||
    !identical(
      sign(as.vector(gaussian_error)),
      sign(as.vector(t_error))
    )) {
  stop("The paired Student-t transform did not preserve error ranks and signs.")
}

message("Re-estimating Gaussian and variance-matched Student-t summaries.")
gaussian_summary <- estimate_eqtl_summaries_from_genotypes(
  G = source_object$genotype,
  expression = source_object$expression,
  covariates = source_object$covariates,
  apply_t_se_correction = TRUE
)
t_summary <- estimate_eqtl_summaries_from_genotypes(
  G = source_object$genotype,
  expression = t_expression,
  covariates = source_object$covariates,
  apply_t_se_correction = TRUE
)

gaussian_reconstruction_error <- c(
  beta_hat = max(abs(
    gaussian_summary$beta_hat - source_object$eqtl_summary$beta_hat
  )),
  raw_se = max(abs(
    gaussian_summary$se_uncorrected -
      source_object$eqtl_summary$se_uncorrected
  )),
  adjusted_se = max(abs(
    gaussian_summary$se - source_object$eqtl_summary$se
  ))
)
if (any(!is.finite(gaussian_reconstruction_error)) ||
    any(gaussian_reconstruction_error > 1e-10)) {
  stop("The formal Gaussian summaries were not reproduced exactly enough.")
}

oracle_se <- compute_oracle_regression_se(
  G = source_object$genotype,
  covariates = source_object$covariates,
  noise_sd = noise_sd,
  n_time = length(time_grid)
)
dimnames(oracle_se) <- dimnames(gaussian_summary$beta_hat)

message("Computing marginal and SE-conditional distribution diagnostics.")
gaussian_error_table <- make_standardized_error_table(
  beta_hat = gaussian_summary$beta_hat,
  true_beta = true_beta,
  se_uncorrected = gaussian_summary$se_uncorrected,
  se_adjusted = gaussian_summary$se,
  oracle_se = oracle_se,
  residual_df = gaussian_summary$df,
  unit_info = unit_info,
  error_distribution = "Gaussian"
)
t_error_table <- make_standardized_error_table(
  beta_hat = t_summary$beta_hat,
  true_beta = true_beta,
  se_uncorrected = t_summary$se_uncorrected,
  se_adjusted = t_summary$se,
  oracle_se = oracle_se,
  residual_df = t_summary$df,
  unit_info = unit_info,
  error_distribution = "Standardized t5"
)
error_table <- rbind(gaussian_error_table, t_error_table)
distribution_summary <- summarize_standardized_errors(error_table)
qq_probabilities <- unique(c(
  seq(0.001, 0.01, length.out = 10L),
  seq(0.015, 0.985, length.out = 195L),
  seq(0.99, 0.999, length.out = 10L)
))
qq_quantiles <- make_qq_quantiles(error_table, qq_probabilities)
se_uncertainty_summary <- summarize_se_uncertainty_bins(
  error_table,
  n_bins = 4L
)

tail_specs <- data.frame(
  nominal_alpha = c(0.05, 0.01, 0.001),
  column = c("tail_rate_0_05", "tail_rate_0_01", "tail_rate_0_001"),
  stringsAsFactors = FALSE
)
tail_calibration <- do.call(rbind, lapply(seq_len(nrow(tail_specs)), function(i) {
  data.frame(
    distribution_summary[, c(
      "error_distribution", "population", "effect_class", "se_scale",
      "reference_distribution", "residual_df", "n"
    )],
    nominal_alpha = tail_specs$nominal_alpha[i],
    empirical_tail_rate = distribution_summary[[tail_specs$column[i]]],
    tail_rate_difference =
      distribution_summary[[tail_specs$column[i]]] -
      tail_specs$nominal_alpha[i],
    stringsAsFactors = FALSE
  )
}))
rownames(tail_calibration) <- NULL

rm(gaussian_error_table, t_error_table, error_table)
invisible(gc())

make_datasets <- function(summary_object, se, scenario) {
  make_fash_datasets_from_eqtl_summary(
    beta_hat = summary_object$beta_hat,
    se = se,
    true_beta = true_beta,
    time_grid = time_grid,
    unit_info = unit_info,
    scenario = scenario
  )
}

summarize_fit_pair <- function(raw_fit,
                               bf_fit,
                               error_distribution,
                               se_scale) {
  stage_fits <- list(Raw = raw_fit, BF = bf_fit)
  result_rows <- list()
  curve_rows <- list()
  alpha005_rows <- list()
  prior_rows <- list()
  for (fit_stage in names(stage_fits)) {
    fit <- stage_fits[[fit_stage]]
    method_label <- paste(
      error_distribution,
      se_scale,
      fit_stage,
      sep = " / "
    )
    result <- evaluate_lfdr_method(
      lfdr = get_fash_lfdr(fit),
      unit_info = unit_info,
      method = method_label,
      target = "dynamic",
      alpha = 0.05
    )
    curve <- compute_alpha_curve(
      result,
      alpha_grid = seq(0, 0.20, by = 0.005)
    )
    curve$error_distribution <- error_distribution
    curve$se_scale <- se_scale
    curve$fit_stage <- fit_stage
    names(curve)[names(curve) == "empirical_fdr"] <- "realized_fdp"
    alpha005 <- curve[abs(curve$alpha - 0.05) < 1e-12, , drop = FALSE]
    prior_rows[[fit_stage]] <- data.frame(
      error_distribution = error_distribution,
      se_scale = se_scale,
      fit_stage = fit_stage,
      estimated_dynamic_null_weight = constant_component_prior_weight(fit),
      stringsAsFactors = FALSE
    )
    result_rows[[fit_stage]] <- result
    curve_rows[[fit_stage]] <- curve
    alpha005_rows[[fit_stage]] <- alpha005
  }
  list(
    unit_results = do.call(rbind, result_rows),
    alpha_curve = do.call(rbind, curve_rows),
    alpha005 = do.call(rbind, alpha005_rows),
    prior = do.call(rbind, prior_rows)
  )
}

fit_missing_arm <- function(summary_object,
                            se,
                            error_distribution,
                            se_scale) {
  scenario <- paste(
    "r1_normality_assumption",
    gsub("[^A-Za-z0-9]+", "_", tolower(error_distribution)),
    gsub("[^A-Za-z0-9]+", "_", tolower(se_scale)),
    sep = "_"
  )
  datasets <- make_datasets(summary_object, se, scenario)
  capture <- capture_warnings(fit_fash_for_revision(
    datasets = datasets,
    orders = 1,
    grid = source_fit_raw$psd_grid,
    num_basis = source_fit_raw$settings$num_basis,
    penalty = source_fit_raw$settings$penalty,
    pred_step = source_fit_raw$settings$pred_step,
    num_cores = num_cores,
    apply_bf = TRUE,
    verbose = FALSE
  ))
  fits <- capture$value
  summary <- summarize_fit_pair(
    raw_fit = fits$fash_iwp1_raw,
    bf_fit = fits$fash_iwp1_bf,
    error_distribution = error_distribution,
    se_scale = se_scale
  )
  summary$warnings <- capture$warnings
  rm(datasets, fits)
  invisible(gc())
  summary
}

load_or_compute_checkpoint <- function(key, computation) {
  checkpoint_path <- file.path(
    checkpoint_directory,
    paste0(key, ".rds")
  )
  if (file.exists(checkpoint_path)) {
    message("Loading completed fit checkpoint: ", key)
    checkpoint <- readRDS(checkpoint_path)
    if (!is.list(checkpoint) ||
        !all(c(
          "unit_results", "alpha_curve", "alpha005", "prior", "warnings"
        ) %in% names(checkpoint))) {
      stop("The fit checkpoint is invalid: ", checkpoint_path)
    }
    return(checkpoint)
  }
  checkpoint <- force(computation)
  saveRDS(checkpoint, checkpoint_path)
  checkpoint
}

message("Summarizing the retained Gaussian/t-adjusted FASH arm.")
fit_summaries <- list()
fit_summaries[["gaussian_adjusted"]] <- summarize_fit_pair(
  raw_fit = source_fit_raw,
  bf_fit = source_fit_bf,
  error_distribution = "Gaussian",
  se_scale = "t-adjusted SE"
)
fit_summaries[["gaussian_adjusted"]]$warnings <- character()

message("Fitting the missing Gaussian/raw-SE arm.")
fit_summaries[["gaussian_raw"]] <- load_or_compute_checkpoint(
  "gaussian_raw",
  fit_missing_arm(
  summary_object = gaussian_summary,
  se = gaussian_summary$se_uncorrected,
  error_distribution = "Gaussian",
  se_scale = "raw regression SE"
  )
)

message("Fitting the Student-t/raw-SE arm.")
fit_summaries[["t_raw"]] <- load_or_compute_checkpoint(
  "t_raw",
  fit_missing_arm(
  summary_object = t_summary,
  se = t_summary$se_uncorrected,
  error_distribution = "Standardized t5",
  se_scale = "raw regression SE"
  )
)

message("Fitting the Student-t/t-adjusted-SE arm.")
fit_summaries[["t_adjusted"]] <- load_or_compute_checkpoint(
  "t_adjusted",
  fit_missing_arm(
  summary_object = t_summary,
  se = t_summary$se,
  error_distribution = "Standardized t5",
  se_scale = "t-adjusted SE"
  )
)

performance_alpha005 <- do.call(rbind, lapply(
  fit_summaries,
  `[[`,
  "alpha005"
))
performance_alpha_curve <- do.call(rbind, lapply(
  fit_summaries,
  `[[`,
  "alpha_curve"
))
prior_summary <- do.call(rbind, lapply(
  fit_summaries,
  `[[`,
  "prior"
))
rownames(performance_alpha005) <- NULL
rownames(performance_alpha_curve) <- NULL
rownames(prior_summary) <- NULL

warning_rows <- do.call(rbind, lapply(names(fit_summaries), function(arm) {
  warnings <- fit_summaries[[arm]]$warnings
  if (length(warnings) == 0L) {
    return(data.frame(
      arm = arm,
      warning = NA_character_,
      stringsAsFactors = FALSE
    ))
  }
  data.frame(arm = arm, warning = warnings, stringsAsFactors = FALSE)
}))

input_error_summary <- data.frame(
  error_distribution = c("Gaussian", "Standardized t5"),
  n = c(length(gaussian_error), length(t_error)),
  mean = c(mean(gaussian_error), mean(t_error)),
  sd = c(stats::sd(as.vector(gaussian_error)), stats::sd(as.vector(t_error))),
  skewness = c(
    mean(((gaussian_error - mean(gaussian_error)) /
      stats::sd(as.vector(gaussian_error)))^3),
    mean(((t_error - mean(t_error)) /
      stats::sd(as.vector(t_error)))^3)
  ),
  excess_kurtosis = c(
    mean(((gaussian_error - mean(gaussian_error)) /
      stats::sd(as.vector(gaussian_error)))^4) - 3,
    mean(((t_error - mean(t_error)) /
      stats::sd(as.vector(t_error)))^4) - 3
  ),
  stringsAsFactors = FALSE
)

configuration <- list(
  experiment = paste(
    "Internal R1 finite-sample normality and standardized Student-t error",
    "pilot"
  ),
  seed = 12345L,
  n_units = 6362L,
  n_donors = 19L,
  n_time = 16L,
  n_covariates = 5L,
  residual_df = 12L,
  t_df = t_df,
  t_scale = sqrt((t_df - 2) / t_df),
  pairing = paste(
    "Common probability-integral transform of the retained Gaussian errors;",
    "signs and ranks are preserved."
  ),
  se_scales = c("oracle known sigma", "raw regression", "t adjusted"),
  performance_scope = paste(
    "One seed; performance quantities are realized FDP and power, not",
    "Monte Carlo empirical FDR."
  ),
  source_path = normalizePath(source_path, winslash = "/", mustWork = TRUE),
  source_fingerprint = artifact_fingerprint(source_path),
  source_scenario = source_object$settings$scenario,
  psd_grid = source_fit_raw$psd_grid,
  num_basis = source_fit_raw$settings$num_basis,
  order = source_fit_raw$settings$order,
  pred_step = source_fit_raw$settings$pred_step,
  penalty = source_fit_raw$settings$penalty,
  num_cores = num_cores,
  r_version = R.version.string,
  package_versions = c(fashr = as.character(utils::packageVersion("fashr"))),
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

validation <- data.frame(
  check = c(
    "gaussian_beta_reconstruction",
    "gaussian_raw_se_reconstruction",
    "gaussian_adjusted_se_reconstruction",
    "student_t_rank_pairing",
    "student_t_sign_pairing",
    "residual_df",
    "fit_warning_count",
    "distribution_summary_finite",
    "qq_quantiles_finite",
    "performance_rows_complete",
    "prior_weights_valid"
  ),
  value = c(
    gaussian_reconstruction_error["beta_hat"],
    gaussian_reconstruction_error["raw_se"],
    gaussian_reconstruction_error["adjusted_se"],
    as.numeric(identical(
      order(as.vector(gaussian_error)),
      order(as.vector(t_error))
    )),
    as.numeric(identical(
      sign(as.vector(gaussian_error)),
      sign(as.vector(t_error))
    )),
    unique(as.vector(gaussian_summary$df)),
    sum(!is.na(warning_rows$warning)),
    as.numeric(all(is.finite(distribution_summary$n)) &&
      all(is.finite(distribution_summary$sd))),
    as.numeric(all(is.finite(qq_quantiles$reference_quantile)) &&
      all(is.finite(qq_quantiles$empirical_quantile))),
    nrow(performance_alpha005),
    as.numeric(all(is.finite(prior_summary$estimated_dynamic_null_weight)) &&
      all(prior_summary$estimated_dynamic_null_weight >= 0) &&
      all(prior_summary$estimated_dynamic_null_weight <= 1))
  ),
  expected = c(
    "<= 1e-10", "<= 1e-10", "<= 1e-10", "1", "1", "12", "0",
    "1", "1", "8", "1"
  ),
  pass = c(
    gaussian_reconstruction_error["beta_hat"] <= 1e-10,
    gaussian_reconstruction_error["raw_se"] <= 1e-10,
    gaussian_reconstruction_error["adjusted_se"] <= 1e-10,
    identical(order(as.vector(gaussian_error)), order(as.vector(t_error))),
    identical(sign(as.vector(gaussian_error)), sign(as.vector(t_error))),
    all(as.integer(gaussian_summary$df) == 12L),
    sum(!is.na(warning_rows$warning)) == 0L,
    all(is.finite(distribution_summary$n)) &&
      all(is.finite(distribution_summary$sd)),
    all(is.finite(qq_quantiles$reference_quantile)) &&
      all(is.finite(qq_quantiles$empirical_quantile)),
    nrow(performance_alpha005) == 8L,
    all(is.finite(prior_summary$estimated_dynamic_null_weight)) &&
      all(prior_summary$estimated_dynamic_null_weight >= 0) &&
      all(prior_summary$estimated_dynamic_null_weight <= 1)
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$pass)) {
  stop(
    "Pilot validation failed: ",
    paste(validation$check[!validation$pass], collapse = ", ")
  )
}

saveRDS(configuration, file.path(staging_directory, "configuration.rds"))
write_csv(input_error_summary, file.path(
  staging_directory,
  "input_error_summary.csv"
))
write_csv(distribution_summary, file.path(
  staging_directory,
  "distribution_summary.csv"
))
write_csv(tail_calibration, file.path(
  staging_directory,
  "tail_calibration.csv"
))
write_csv(qq_quantiles, file.path(staging_directory, "qq_quantiles.csv"))
write_csv(se_uncertainty_summary, file.path(
  staging_directory,
  "se_uncertainty_summary.csv"
))
write_csv(performance_alpha005, file.path(
  staging_directory,
  "performance_alpha005.csv"
))
write_csv(performance_alpha_curve, file.path(
  staging_directory,
  "performance_alpha_curve.csv"
))
write_csv(prior_summary, file.path(staging_directory, "prior_summary.csv"))
write_csv(warning_rows, file.path(staging_directory, "fit_warnings.csv"))
write_csv(validation, file.path(staging_directory, "validation.csv"))

if (!file.rename(staging_directory, final_directory)) {
  stop("Could not finalize the pilot cache directory.")
}
unlink(checkpoint_directory, recursive = TRUE, force = TRUE)
completed <- TRUE

cat("\nR1 normality-assumption seed-12345 pilot completed.\n")
cat("Output:", final_directory, "\n")
cat("Gaussian summary reconstruction max errors:\n")
print(gaussian_reconstruction_error)
cat("\nInput error summaries:\n")
print(input_error_summary)
cat("\nZero-null adjusted-z summaries:\n")
print(distribution_summary[
  distribution_summary$population == "zero_null" &
    distribution_summary$se_scale == "t_adjusted",
  c(
    "error_distribution", "n", "mean", "sd", "skewness",
    "excess_kurtosis", "tail_rate_0_05", "tail_rate_0_01"
  )
])
cat("\nFASH realized outcomes at alpha 0.05:\n")
print(performance_alpha005[, c(
  "error_distribution", "se_scale", "fit_stage", "n_discoveries",
  "false_discoveries", "realized_fdp", "power"
)])
