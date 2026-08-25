#!/usr/bin/env Rscript

# Run a known-truth R1 shared-genotype-permutation calibration pilot.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
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

capture_warnings <- function(expression) {
  warning_messages <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(warning) {
      warning_messages <<- c(warning_messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warning_messages))
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

file_metadata <- function(paths, roles) {
  normalized_paths <- normalizePath(paths, mustWork = TRUE)
  information <- file.info(normalized_paths)
  data.frame(
    role = roles,
    path = normalized_paths,
    size_bytes = unname(information$size),
    modification_time = format(information$mtime, tz = "UTC", usetz = TRUE),
    md5 = unname(tools::md5sum(normalized_paths)),
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation",
  "selected_signal_genotype_permutation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "r1_known_truth_genotype_permutation",
  "r1_known_truth_genotype_permutation_helpers.R"
))

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required.")
}

source_seed <- 12345L
permutation_seed <- as.integer(get_arg("--permutation-seed", "20260811"))
n_alternatives <- as.integer(get_arg("--n-alternatives", "200"))
num_cores <- as.integer(get_arg("--num-cores", "8"))
if (length(permutation_seed) != 1L || is.na(permutation_seed) ||
    length(n_alternatives) != 1L || is.na(n_alternatives) ||
    n_alternatives < 2L || n_alternatives > 200L ||
    length(num_cores) != 1L || is.na(num_cores) || num_cores < 1L) {
  stop("Invalid permutation seed, alternative count, or core count.")
}

output_id <- paste0(
  "r1_known_truth_genotype_permutation_seed", source_seed,
  "_perm", permutation_seed,
  "_J", n_alternatives
)
source_fit_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "mc",
  "r1_random_bspline_main_effect_profile_sigma_pilot5",
  "full_fits", paste0("seed_", source_seed, ".rds")
)
simulation_functions_path <- file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
)
selected_helper_path <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation",
  "selected_signal_genotype_permutation_helpers.R"
)
known_truth_helper_path <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "r1_known_truth_genotype_permutation",
  "r1_known_truth_genotype_permutation_helpers.R"
)
required_paths <- c(
  source_fit_path,
  simulation_functions_path,
  selected_helper_path,
  known_truth_helper_path
)
required_roles <- c(
  "cached_r1_full_fit",
  "shared_simulation_functions",
  "shared_genotype_permutation_helpers",
  "known_truth_pilot_helpers"
)
if (any(!file.exists(required_paths))) {
  stop(
    "At least one required input is missing: ",
    paste(required_paths[!file.exists(required_paths)], collapse = ", ")
  )
}

output_parent <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal"
)
final_directory <- file.path(output_parent, output_id)
staging_directory <- file.path(
  output_parent,
  paste0(".", output_id, "_staging_", Sys.getpid())
)
if (file.exists(final_directory) || file.exists(staging_directory)) {
  stop("Refusing to overwrite an existing output or staging directory.")
}
dir.create(staging_directory, recursive = FALSE)
staging_complete <- FALSE
on.exit({
  if (!staging_complete && dir.exists(staging_directory)) {
    unlink(staging_directory, recursive = TRUE)
  }
}, add = TRUE)

analysis_start <- proc.time()[["elapsed"]]
message("Loading and validating cached R1 seed 12345.")
source_information <- file_metadata(required_paths, required_roles)
source_object <- readRDS(source_fit_path)
source_fit <- source_object$fash_fits$fash_iwp1_raw
expected_grid <- default_revision_grid()
effect_class <- as.character(source_object$unit_info$effect_class)
alternative_indices_all <- which(effect_class == "dynamic_bspline")
zero_indices_all <- which(effect_class == "zero")
expression_correlation <- source_object$settings$expression_error_correlation
if (!identical(dim(source_object$genotype), c(19L, 1000L)) ||
    !identical(dim(source_object$expression), c(19L, 1000L, 16L)) ||
    !identical(dim(source_object$covariates), c(19L, 5L)) ||
    length(alternative_indices_all) != 200L ||
    length(zero_indices_all) != 400L ||
    !identical(source_fit$settings$order, 1) ||
    !identical(source_fit$settings$num_basis, 20L) ||
    !isTRUE(all.equal(source_fit$psd_grid, expected_grid)) ||
    !identical(source_fit$settings$pred_step, 1) ||
    !identical(source_fit$settings$penalty, 10) ||
    !isTRUE(all.equal(expression_correlation, diag(16), tolerance = 0)) ||
    nrow(source_fit$L_matrix) != 1000L ||
    ncol(source_fit$L_matrix) != length(expected_grid) ||
    length(source_fit$fash_data$data_list) != 1000L ||
    length(source_fit$fash_data$S) != 1000L) {
  stop("The cached R1 source object does not match the expected fixed version.")
}

alternative_indices <- alternative_indices_all[seq_len(n_alternatives)]
genuine_null_indices <- zero_indices_all[seq_len(n_alternatives)]
alternative_keys <- source_object$unit_info$unit_id[alternative_indices]
genuine_null_keys <- source_object$unit_info$unit_id[genuine_null_indices]
if (!identical(
      names(source_fit$fash_data$data_list),
      source_object$unit_info$unit_id
    ) || anyDuplicated(c(alternative_keys, genuine_null_keys))) {
  stop("The cached R1 FASH rows are not aligned with source unit identifiers.")
}

message("Reproducing cached alternative beta and SE values from individual data.")
alternative_genotype <- source_object$genotype[
  , alternative_indices, drop = FALSE
]
colnames(alternative_genotype) <- alternative_keys
alternative_expression <- source_object$expression[
  , alternative_indices, , drop = FALSE
]
original_eqtl <- estimate_eqtl_summaries_from_genotypes(
  G = alternative_genotype,
  expression = alternative_expression,
  covariates = source_object$covariates,
  apply_t_se_correction = TRUE
)
cached_beta <- source_object$eqtl_summary$beta_hat[
  alternative_indices, , drop = FALSE
]
cached_se <- source_object$eqtl_summary$se[
  alternative_indices, , drop = FALSE
]
maximum_original_beta_difference <- max(abs(
  original_eqtl$beta_hat - cached_beta
))
maximum_original_se_difference <- max(abs(original_eqtl$se - cached_se))
if (maximum_original_beta_difference > 1e-12 ||
    maximum_original_se_difference > 1e-12) {
  stop("The individual-level R1 refit did not reproduce cached beta/SE values.")
}

message("Applying one shared donor permutation to all 200 alternative genotypes.")
shared_permutation <- make_shared_genotype_permutation(
  alternative_genotype,
  seed = permutation_seed
)
permuted_genotype <- shared_permutation$genotype
donor_map <- shared_permutation$donor_map
if (all(donor_map$fixed_point)) {
  stop("The fixed-seed donor permutation is the identity permutation.")
}
alignment <- residualized_genotype_alignment(
  alternative_genotype,
  permuted_genotype,
  source_object$covariates
)

message("Estimating eQTL summaries for genotype-permuted copies.")
permuted_eqtl <- estimate_eqtl_summaries_from_genotypes(
  G = permuted_genotype,
  expression = alternative_expression,
  covariates = source_object$covariates,
  apply_t_se_correction = TRUE
)
permuted_unit_info <- source_object$unit_info[
  alternative_indices, , drop = FALSE
]
permuted_unit_info$unit_id <- paste0(
  alternative_keys,
  "__permuted_null_seed", permutation_seed
)
permuted_unit_info$effect_class <- "permuted_null"
permuted_unit_info$genetic_main_effect <- NA_real_
permuted_unit_info$scenario <- "r1_known_truth_shared_genotype_permutation"
permuted_true_beta <- matrix(
  0,
  nrow = n_alternatives,
  ncol = 16L,
  dimnames = dimnames(permuted_eqtl$beta_hat)
)
permuted_datasets <- make_fash_datasets_from_eqtl_summary(
  beta_hat = permuted_eqtl$beta_hat,
  se = permuted_eqtl$se,
  true_beta = permuted_true_beta,
  time_grid = source_object$settings$time_grid,
  unit_info = permuted_unit_info,
  scenario = "r1_known_truth_shared_genotype_permutation"
)
permuted_keys <- permuted_unit_info$unit_id
names(permuted_datasets) <- permuted_keys

message("Computing 200 new FASH likelihood rows for permuted copies.")
permuted_likelihood_start <- proc.time()[["elapsed"]]
permuted_capture <- capture_warnings(fashr::fash(
  Y = "y",
  smooth_var = "x",
  S = "sd",
  data_list = permuted_datasets,
  num_basis = source_fit$settings$num_basis,
  order = source_fit$settings$order,
  betaprec = source_fit$settings$betaprec,
  pred_step = source_fit$settings$pred_step,
  penalty = source_fit$settings$penalty,
  grid = source_fit$psd_grid,
  num_cores = num_cores,
  verbose = TRUE
))
permuted_likelihood_elapsed <-
  proc.time()[["elapsed"]] - permuted_likelihood_start
permuted_fit <- permuted_capture$value
if (nrow(permuted_fit$L_matrix) != n_alternatives ||
    ncol(permuted_fit$L_matrix) != length(source_fit$psd_grid) ||
    !isTRUE(all.equal(permuted_fit$psd_grid, source_fit$psd_grid))) {
  stop("The permuted-copy FASH likelihood matrix is invalid.")
}
names(permuted_fit$fash_data$data_list) <- permuted_keys
names(permuted_fit$fash_data$S) <- permuted_keys
rownames(permuted_fit$L_matrix) <- permuted_keys

alternative_data <- source_fit$fash_data$data_list[alternative_indices]
alternative_se <- source_fit$fash_data$S[alternative_indices]
alternative_likelihood <- source_fit$L_matrix[
  alternative_indices, , drop = FALSE
]
names(alternative_data) <- alternative_keys
names(alternative_se) <- alternative_keys
rownames(alternative_likelihood) <- alternative_keys

genuine_null_data <- source_fit$fash_data$data_list[genuine_null_indices]
genuine_null_se <- source_fit$fash_data$S[genuine_null_indices]
genuine_null_likelihood <- source_fit$L_matrix[
  genuine_null_indices, , drop = FALSE
]
names(genuine_null_data) <- genuine_null_keys
names(genuine_null_se) <- genuine_null_keys
rownames(genuine_null_likelihood) <- genuine_null_keys

arm_inputs <- list(
  genuine_null_baseline = list(
    unit_keys = c(alternative_keys, genuine_null_keys),
    source_unit_id = c(alternative_keys, genuine_null_keys),
    group = rep(c("known_alternative", "genuine_zero_null"),
                each = n_alternatives),
    data_list = c(alternative_data, genuine_null_data),
    se_list = c(alternative_se, genuine_null_se),
    likelihood = rbind(alternative_likelihood, genuine_null_likelihood)
  ),
  shared_genotype_permutation = list(
    unit_keys = c(alternative_keys, permuted_keys),
    source_unit_id = rep(alternative_keys, 2L),
    group = rep(c("known_alternative", "permuted_null"),
                each = n_alternatives),
    data_list = c(alternative_data, permuted_fit$fash_data$data_list),
    se_list = c(alternative_se, permuted_fit$fash_data$S),
    likelihood = rbind(alternative_likelihood, permuted_fit$L_matrix)
  )
)
true_null <- c(
  rep(FALSE, n_alternatives),
  rep(TRUE, n_alternatives)
)
alpha_grid <- c(0, seq(0.001, 0.20, by = 0.001))

message("Refitting raw and BF-adjusted mixtures for both known-truth arms.")
arm_results <- list()
for (arm in names(arm_inputs)) {
  input <- arm_inputs[[arm]]
  raw_start <- proc.time()[["elapsed"]]
  raw_capture <- capture_warnings(refit_fash_from_likelihood(
    source_fit = source_fit,
    data_list = input$data_list,
    se_list = input$se_list,
    likelihood_matrix = input$likelihood,
    unit_keys = input$unit_keys,
    penalty = source_fit$settings$penalty
  ))
  raw_elapsed <- proc.time()[["elapsed"]] - raw_start
  bf_start <- proc.time()[["elapsed"]]
  bf_capture <- capture_warnings(fashr::BF_update(
    raw_capture$value,
    plot = FALSE
  ))
  bf_elapsed <- proc.time()[["elapsed"]] - bf_start
  names(bf_capture$value$lfdr) <- input$unit_keys
  arm_results[[arm]] <- list(
    raw_fit = raw_capture$value,
    bf_fit = bf_capture$value,
    true_null = true_null,
    unit_keys = input$unit_keys,
    source_unit_id = input$source_unit_id,
    group = input$group,
    raw_warnings = raw_capture$warnings,
    bf_warnings = bf_capture$warnings,
    raw_elapsed_seconds = unname(raw_elapsed),
    bf_elapsed_seconds = unname(bf_elapsed)
  )
}

pi0_rows <- list()
alpha_rows <- list()
truth_rows <- list()
bf_summary_rows <- list()
bf_unit_rows <- list()
row_index <- 0L
for (arm in names(arm_results)) {
  result <- arm_results[[arm]]
  truth_rows[[arm]] <- data.frame(
    arm = arm,
    unit_key = result$unit_keys,
    source_unit_id = result$source_unit_id,
    group = result$group,
    true_null = result$true_null,
    stringsAsFactors = FALSE
  )
  for (fit_stage in c("Raw", "BF")) {
    fit <- if (fit_stage == "Raw") result$raw_fit else result$bf_fit
    row_index <- row_index + 1L
    pi0_rows[[row_index]] <- data.frame(
      arm = arm,
      fit_stage = fit_stage,
      estimated_pi0 = extract_pi0(fit),
      true_pi0 = mean(result$true_null),
      difference_from_true_pi0 = extract_pi0(fit) - mean(result$true_null),
      stringsAsFactors = FALSE
    )
    alpha_rows[[row_index]] <- known_truth_alpha_curve(
      lfdr = fit$lfdr,
      true_null = result$true_null,
      alpha_grid = alpha_grid,
      arm = arm,
      fit_stage = fit_stage
    )
  }
  bf_summary_rows[[arm]] <- summarize_null_bf(
    bf = result$bf_fit$BF,
    true_null = result$true_null,
    arm = arm
  )
  bf_unit_rows[[arm]] <- data.frame(
    arm = arm,
    unit_key = result$unit_keys[result$true_null],
    source_unit_id = result$source_unit_id[result$true_null],
    group = result$group[result$true_null],
    bayes_factor = result$bf_fit$BF[result$true_null],
    stringsAsFactors = FALSE
  )
}
pi0_summary <- do.call(rbind, pi0_rows)
alpha_curve <- do.call(rbind, alpha_rows)
alpha005_summary <- alpha_curve[
  abs(alpha_curve$alpha - 0.05) < 1e-12, , drop = FALSE
]
unit_truth <- do.call(rbind, truth_rows)
null_bf_summary <- do.call(rbind, bf_summary_rows)
null_bf_units <- do.call(rbind, bf_unit_rows)
rownames(pi0_summary) <- NULL
rownames(alpha_curve) <- NULL
rownames(alpha005_summary) <- NULL
rownames(unit_truth) <- NULL
rownames(null_bf_summary) <- NULL
rownames(null_bf_units) <- NULL

permuted_result <- arm_results$shared_genotype_permutation
permuted_null_bf <- permuted_result$bf_fit$BF[permuted_result$true_null]
true_beta <- source_object$true_beta[alternative_indices, , drop = FALSE]
predicted_leakage_z <- sweep(
  true_beta,
  1L,
  alignment$leakage_coefficient,
  `*`
) / permuted_eqtl$se
actual_permuted_z <- permuted_eqtl$beta_hat / permuted_eqtl$se
alignment$source_unit_id <- alternative_keys
alignment$permuted_unit_key <- permuted_keys
alignment$permuted_bayes_factor <- permuted_null_bf
alignment$log10_permuted_bayes_factor <- log10(permuted_null_bf)
alignment$maximum_absolute_predicted_leakage_z <- apply(
  abs(predicted_leakage_z), 1L, max
)
alignment$rms_predicted_leakage_z <- sqrt(rowMeans(predicted_leakage_z^2))
alignment$maximum_absolute_actual_permuted_z <- apply(
  abs(actual_permuted_z), 1L, max
)
alignment$rms_actual_permuted_z <- sqrt(rowMeans(actual_permuted_z^2))

analysis_elapsed <- proc.time()[["elapsed"]] - analysis_start
all_warning_messages <- unique(c(
  permuted_capture$warnings,
  unlist(lapply(arm_results, `[[`, "raw_warnings")),
  unlist(lapply(arm_results, `[[`, "bf_warnings"))
))
configuration <- list(
  experiment = "R1 known-truth shared-genotype-permutation calibration pilot",
  source_seed = source_seed,
  permutation_seed = permutation_seed,
  n_alternatives = n_alternatives,
  n_null_per_arm = n_alternatives,
  n_units_per_arm = 2L * n_alternatives,
  true_pi0 = 0.5,
  source_scenario = source_object$settings$scenario,
  null_arms = c(
    "genuine_zero_units",
    "shared_genotype_permutation_of_known_alternative_units"
  ),
  donor_permutation = paste(
    "One donor-row permutation of the complete alternative genotype matrix,",
    "shared across all units and 16 time points."
  ),
  expression_error_correlation = expression_correlation,
  source_fash_settings = source_fit$settings,
  psd_grid = source_fit$psd_grid,
  num_cores = num_cores,
  permuted_likelihood_elapsed_seconds = unname(permuted_likelihood_elapsed),
  arm_refit_elapsed_seconds = vapply(
    arm_results,
    function(x) x$raw_elapsed_seconds + x$bf_elapsed_seconds,
    numeric(1)
  ),
  total_elapsed_seconds = unname(analysis_elapsed),
  warnings = all_warning_messages,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  interpretation_boundary = paste(
    "One source simulation seed and one donor-permutation seed provide an",
    "internal mechanism diagnostic, not repeated-simulation FDR evidence."
  )
)

validation <- data.frame(
  component = c(
    "R1 genotype dimensions",
    "R1 expression dimensions",
    "Known alternatives",
    "Known nulls per arm",
    "True pi0 per arm",
    "Original beta reproduction",
    "Original SE reproduction",
    "Permutation is non-identity",
    "Permutation genotype cross-products",
    "Arm fits",
    "Alpha 0.05 rows",
    "Positive finite null Bayes factors",
    "Independent expression-error setting"
  ),
  observed = c(
    paste(dim(source_object$genotype), collapse = " x "),
    paste(dim(source_object$expression), collapse = " x "),
    n_alternatives,
    n_alternatives,
    paste(unique(pi0_summary$true_pi0), collapse = ","),
    format(maximum_original_beta_difference, digits = 16),
    format(maximum_original_se_difference, digits = 16),
    sum(!donor_map$fixed_point),
    format(max(abs(
      crossprod(alternative_genotype) - crossprod(permuted_genotype)
    )), digits = 16),
    length(arm_results),
    nrow(alpha005_summary),
    sum(is.finite(null_bf_units$bayes_factor) &
          null_bf_units$bayes_factor > 0),
    max(abs(expression_correlation - diag(16)))
  ),
  expected = c(
    "19 x 1000",
    "19 x 1000 x 16",
    as.character(n_alternatives),
    as.character(n_alternatives),
    "0.5",
    "<= 1e-12",
    "<= 1e-12",
    "> 0",
    "<= 1e-12",
    "2",
    "4",
    as.character(2L * n_alternatives),
    "0"
  ),
  pass = c(
    identical(dim(source_object$genotype), c(19L, 1000L)),
    identical(dim(source_object$expression), c(19L, 1000L, 16L)),
    length(alternative_indices) == n_alternatives,
    length(genuine_null_indices) == n_alternatives,
    all(abs(pi0_summary$true_pi0 - 0.5) < 1e-12),
    maximum_original_beta_difference <= 1e-12,
    maximum_original_se_difference <= 1e-12,
    any(!donor_map$fixed_point),
    max(abs(crossprod(alternative_genotype) -
              crossprod(permuted_genotype))) <= 1e-12,
    length(arm_results) == 2L,
    nrow(alpha005_summary) == 4L,
    all(is.finite(null_bf_units$bayes_factor) &
          null_bf_units$bayes_factor > 0),
    max(abs(expression_correlation - diag(16))) <= 1e-12
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$pass)) {
  stop(
    "Internal validation failed: ",
    paste(validation$component[!validation$pass], collapse = ", ")
  )
}

effect_estimates <- list(
  source_unit_id = alternative_keys,
  true_beta = true_beta,
  cached_original_beta = cached_beta,
  cached_original_se = cached_se,
  reproduced_original_beta = original_eqtl$beta_hat,
  reproduced_original_se = original_eqtl$se,
  permuted_beta = permuted_eqtl$beta_hat,
  permuted_se_uncorrected = permuted_eqtl$se_uncorrected,
  permuted_se = permuted_eqtl$se,
  predicted_leakage_z = predicted_leakage_z,
  actual_permuted_z = actual_permuted_z
)
arm_fit_bundle <- lapply(arm_results, function(x) {
  x[c(
    "raw_fit", "bf_fit", "true_null", "unit_keys", "source_unit_id",
    "group", "raw_warnings", "bf_warnings", "raw_elapsed_seconds",
    "bf_elapsed_seconds"
  )]
})

saveRDS(configuration, file.path(staging_directory, "configuration.rds"))
saveRDS(effect_estimates, file.path(staging_directory, "effect_estimates.rds"))
saveRDS(arm_fit_bundle, file.path(staging_directory, "arm_fits.rds"))
write_csv(donor_map, file.path(staging_directory, "donor_permutation.csv"))
write_csv(unit_truth, file.path(staging_directory, "unit_truth.csv"))
write_csv(pi0_summary, file.path(staging_directory, "pi0_summary.csv"))
write_csv(alpha_curve, file.path(staging_directory, "alpha_curve.csv"))
write_csv(
  alpha005_summary,
  file.path(staging_directory, "alpha005_summary.csv")
)
write_csv(
  null_bf_summary,
  file.path(staging_directory, "null_bf_summary.csv")
)
write_csv(null_bf_units, file.path(staging_directory, "null_bf_units.csv"))
write_csv(alignment, file.path(staging_directory, "genotype_alignment.csv"))
write_csv(validation, file.path(staging_directory, "validation.csv"))
write_csv(
  source_information,
  file.path(staging_directory, "source_information.csv")
)

if (!file.rename(staging_directory, final_directory)) {
  stop("Failed to atomically finalize the known-truth pilot output.")
}
staging_complete <- TRUE

cat("\nR1 known-truth genotype-permutation pilot completed.\n")
cat("Output: ", final_directory, "\n", sep = "")
cat("\nEstimated pi0:\n")
print(pi0_summary, row.names = FALSE)
cat("\nAlpha 0.05 calibration:\n")
print(alpha005_summary, row.names = FALSE)
cat("\nKnown-null BF summaries:\n")
print(null_bf_summary, row.names = FALSE)
