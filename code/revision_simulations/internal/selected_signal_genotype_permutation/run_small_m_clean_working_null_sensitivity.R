#!/usr/bin/env Rscript

# Run the single-draw small-M clean working-model null sensitivity analysis.

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

write_csv <- function(value, path) {
  utils::write.csv(value, path, row.names = FALSE, na = "")
}

capture_warnings <- function(expression) {
  warning_messages <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(condition) {
      warning_messages <<- c(warning_messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warning_messages))
}

make_input_provenance <- function(paths, roles) {
  information <- file.info(paths)
  data.frame(
    role = roles,
    path = normalizePath(paths, mustWork = TRUE),
    size_bytes = as.numeric(information$size),
    modified_at = format(information$mtime, tz = "UTC", usetz = TRUE),
    md5 = unname(tools::md5sum(paths)),
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
script_directory <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation"
)
source(file.path(
  script_directory,
  "selected_signal_genotype_permutation_helpers.R"
))
source(file.path(script_directory, "small_m_eb_sensitivity_helpers.R"))
source(file.path(
  script_directory,
  "small_m_clean_working_null_helpers.R"
))

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required.")
}

source_output_id <- paste0(
  "all_gene_random_variant_signal_stripped_residual_block_permutation_",
  "selection20260817_seed20260811"
)
source_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  source_output_id
)
source_fit_path <- file.path(source_directory, "merged_fash_fit.rds")
source_selection_path <- file.path(source_directory, "selection.csv")
input_paths <- c(source_fit_path, source_selection_path)
input_roles <- c(
  "source_target_likelihood_and_settings",
  "one_random_variant_per_gene_selection"
)
if (any(!file.exists(input_paths))) {
  stop("At least one required source input is missing.")
}

output_id <- paste0(
  "all_gene_random_variant_small_m_clean_working_null_",
  "selection20260817_nullseed20260819"
)
output_parent <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal"
)
final_directory <- file.path(output_parent, output_id)
checkpoint_path <- file.path(
  output_parent,
  paste0(output_id, "_likelihood_checkpoint.rds")
)
if (dir.exists(final_directory)) {
  stop("The completed clean-null output already exists: ", final_directory)
}
staging_directory <- paste0(final_directory, ".staging_", Sys.getpid())
if (dir.exists(staging_directory)) {
  stop("The clean-null staging directory already exists: ", staging_directory)
}
dir.create(staging_directory, recursive = TRUE, showWarnings = FALSE)
staging_complete <- FALSE
on.exit({
  if (!staging_complete && dir.exists(staging_directory)) {
    unlink(staging_directory, recursive = TRUE)
  }
}, add = TRUE)

analysis_start <- proc.time()[["elapsed"]]
input_provenance <- make_input_provenance(input_paths, input_roles)
source_bundle <- readRDS(source_fit_path)
selection <- utils::read.csv(source_selection_path, stringsAsFactors = FALSE)

J <- as.integer(source_bundle$configuration$n_target_units)
expected_J <- 6362L
source_raw <- source_bundle$raw_fit
all_source_keys <- rownames(source_raw$L_matrix)
target_indices <- seq_len(J)
target_keys <- all_source_keys[target_indices]
target_data <- source_raw$fash_data$data_list[target_indices]
target_se <- source_raw$fash_data$S[target_indices]
target_likelihood <- source_raw$L_matrix[target_indices, , drop = FALSE]
names(target_data) <- target_keys
names(target_se) <- target_keys
rownames(target_likelihood) <- target_keys

required_settings <- c(
  "num_basis", "order", "betaprec", "pred_step", "likelihood", "penalty"
)
target_lengths <- vapply(target_data, nrow, integer(1))
target_se_lengths <- vapply(target_se, length, integer(1))
if (!identical(source_bundle$configuration$target_selection_method,
               "random_all_genes") ||
    J != expected_J || nrow(selection) != J ||
    !all(required_settings %in% names(source_raw$settings)) ||
    nrow(source_raw$L_matrix) != 2L * J ||
    ncol(source_raw$L_matrix) != 52L ||
    length(all_source_keys) != 2L * J || anyDuplicated(all_source_keys) ||
    !identical(target_keys, selection$pair_key) ||
    !identical(names(target_data), target_keys) ||
    !identical(names(target_se), target_keys) ||
    !all(target_lengths == 16L) ||
    !identical(target_lengths, target_se_lengths) ||
    any(!vapply(target_se, function(value) {
      all(is.finite(value)) && all(value > 0)
    }, logical(1))) ||
    anyNA(target_likelihood) || any(is.nan(target_likelihood)) ||
    any(target_likelihood == Inf)) {
  stop("The cached one-random-variant-per-gene target universe is invalid.")
}

null_seed <- 20260819L
subset_seed <- 2026081901L
num_cores <- 8L
alpha_grid <- seq(0.001, 0.200, by = 0.001)
ratio_to_size <- c(
  "0.05" = as.integer(round(0.05 * J)),
  "0.10" = as.integer(round(0.10 * J)),
  "0.20" = as.integer(round(0.20 * J)),
  "1.00" = J
)

message("Simulating one iid Gaussian working-model null pool.")
clean_null <- simulate_clean_working_null(
  target_data_list = target_data,
  target_se_list = target_se,
  source_unit_keys = target_keys,
  seed = null_seed
)
z_diagnostics <- summarize_clean_null_z(clean_null$z_matrix)
membership <- make_nested_null_subsets(
  n_null = J,
  ratio_to_size = ratio_to_size,
  n_replicates = 1L,
  seed = subset_seed
)
membership$null_unit_key <- clean_null$null_unit_keys[
  membership$null_index
]

source_fit_md5 <- unname(tools::md5sum(source_fit_path))
checkpoint_reused <- file.exists(checkpoint_path)
if (checkpoint_reused) {
  message("Reusing the validated clean-null likelihood checkpoint.")
  null_cache <- readRDS(checkpoint_path)
  if (!identical(null_cache$cache_version,
                 "small_m_clean_working_null_v1") ||
      !identical(null_cache$source_fit_md5, source_fit_md5) ||
      !identical(null_cache$null_seed, null_seed) ||
      !identical(null_cache$subset_seed, subset_seed) ||
      !identical(null_cache$source_unit_keys, target_keys) ||
      !identical(null_cache$null_unit_keys, clean_null$null_unit_keys) ||
      !isTRUE(all.equal(
        null_cache$z_matrix,
        clean_null$z_matrix,
        tolerance = 0
      ))) {
    stop("The clean-null likelihood checkpoint does not match this run.")
  }
  null_fit <- null_cache$fit
  null_likelihood_seconds <- null_cache$elapsed_seconds
  null_capture <- list(
    value = null_fit,
    warnings = null_cache$warnings
  )
} else {
  message("Computing 6,362 clean-null likelihood rows.")
  null_likelihood_start <- proc.time()[["elapsed"]]
  null_capture <- capture_warnings(fashr::fash(
    Y = "y",
    smooth_var = "x",
    S = clean_null$se_list,
    data_list = clean_null$data_list,
    num_basis = source_raw$settings$num_basis,
    order = source_raw$settings$order,
    betaprec = source_raw$settings$betaprec,
    pred_step = source_raw$settings$pred_step,
    likelihood = source_raw$settings$likelihood,
    penalty = source_raw$settings$penalty,
    grid = source_raw$psd_grid,
    num_cores = num_cores,
    verbose = TRUE
  ))
  null_likelihood_seconds <-
    proc.time()[["elapsed"]] - null_likelihood_start
  null_fit <- null_capture$value
}
names(null_fit$fash_data$data_list) <- clean_null$null_unit_keys
names(null_fit$fash_data$S) <- clean_null$null_unit_keys
rownames(null_fit$L_matrix) <- clean_null$null_unit_keys
names(null_fit$lfdr) <- clean_null$null_unit_keys
if (nrow(null_fit$L_matrix) != J || ncol(null_fit$L_matrix) != 52L ||
    length(null_fit$fash_data$data_list) != J ||
    length(null_fit$fash_data$S) != J ||
    !identical(rownames(null_fit$L_matrix), clean_null$null_unit_keys) ||
    !isTRUE(all.equal(null_fit$psd_grid, source_raw$psd_grid, tolerance = 0)) ||
    anyNA(null_fit$L_matrix) || any(is.nan(null_fit$L_matrix)) ||
    any(null_fit$L_matrix == Inf)) {
  stop("The clean-null likelihood fit failed validation.")
}

null_cache <- list(
  cache_version = "small_m_clean_working_null_v1",
  source_fit_path = normalizePath(source_fit_path, mustWork = TRUE),
  source_fit_md5 = source_fit_md5,
  null_seed = null_seed,
  subset_seed = subset_seed,
  source_unit_keys = target_keys,
  null_unit_keys = clean_null$null_unit_keys,
  z_matrix = clean_null$z_matrix,
  fit = null_fit,
  elapsed_seconds = null_likelihood_seconds,
  warnings = null_capture$warnings
)
if (!file.exists(checkpoint_path)) {
  temporary_checkpoint <- paste0(checkpoint_path, ".tmp_", Sys.getpid())
  saveRDS(null_cache, temporary_checkpoint)
  if (!file.rename(temporary_checkpoint, checkpoint_path)) {
    unlink(temporary_checkpoint)
    stop("Could not atomically save the clean-null likelihood checkpoint.")
  }
  message("Saved the validated clean-null likelihood checkpoint.")
}

fit_one_ratio <- function(m_ratio, null_indices) {
  null_indices <- as.integer(null_indices)
  unit_keys <- c(target_keys, clean_null$null_unit_keys[null_indices])
  merged_data <- c(target_data, null_fit$fash_data$data_list[null_indices])
  merged_se <- c(target_se, null_fit$fash_data$S[null_indices])
  merged_likelihood <- rbind(
    target_likelihood,
    null_fit$L_matrix[null_indices, , drop = FALSE]
  )
  rownames(merged_likelihood) <- unit_keys
  names(merged_data) <- unit_keys
  names(merged_se) <- unit_keys

  raw_start <- proc.time()[["elapsed"]]
  raw_capture <- capture_warnings(refit_fash_from_likelihood(
    source_fit = source_raw,
    data_list = merged_data,
    se_list = merged_se,
    likelihood_matrix = merged_likelihood,
    unit_keys = unit_keys,
    penalty = source_raw$settings$penalty
  ))
  raw_seconds <- proc.time()[["elapsed"]] - raw_start
  bf_start <- proc.time()[["elapsed"]]
  bf_capture <- capture_warnings(fashr::BF_update(
    raw_capture$value,
    plot = FALSE
  ))
  bf_seconds <- proc.time()[["elapsed"]] - bf_start
  names(bf_capture$value$lfdr) <- unit_keys
  list(
    m_ratio = m_ratio,
    m_size = length(null_indices),
    null_indices = null_indices,
    raw_fit = raw_capture$value,
    bf_fit = bf_capture$value,
    raw_seconds = raw_seconds,
    bf_seconds = bf_seconds,
    raw_warnings = raw_capture$warnings,
    bf_warnings = bf_capture$warnings
  )
}

message("Refitting the target-only reference and four nested merged fits.")
fit_results <- list()
fit_results[["0.00"]] <- fit_one_ratio(0, integer())
for (ratio_name in names(ratio_to_size)) {
  current_ratio <- as.numeric(ratio_name)
  current_indices <- membership$null_index[
    abs(membership$m_ratio - current_ratio) < 1e-12
  ]
  fit_results[[sprintf("%.2f", current_ratio)]] <- fit_one_ratio(
    current_ratio,
    current_indices
  )
}

make_pi0_row <- function(fit, m_ratio, m_size, fit_stage) {
  pi0_merged <- extract_pi0(fit)
  pi0_target <- if (m_size == 0L) {
    pi0_merged
  } else {
    ((J + m_size) * pi0_merged - m_size) / J
  }
  data.frame(
    m_ratio = m_ratio,
    m_size = m_size,
    fit_stage = fit_stage,
    pi0_merged = pi0_merged,
    pi0_target_unbounded = pi0_target,
    pi0_target_valid = pi0_target >= 0 && pi0_target <= 1,
    pi0_target_bounded = min(1, max(0, pi0_target)),
    stringsAsFactors = FALSE
  )
}

make_prior_rows <- function(fit, m_ratio, m_size, fit_stage) {
  data.frame(
    m_ratio = m_ratio,
    m_size = m_size,
    fit_stage = fit_stage,
    psd = as.numeric(fit$prior_weights$psd),
    prior_weight = as.numeric(fit$prior_weights$prior_weight),
    component_type = ifelse(
      fit$prior_weights$psd == 0,
      "Exact null",
      "Alternative"
    ),
    stringsAsFactors = FALSE
  )
}

pi0_rows <- list()
prior_rows <- list()
curve_rows <- list()
runtime_rows <- list()
warning_rows <- list()
fit_summary_bundle <- list()
row_index <- 0L
for (ratio_name in names(fit_results)) {
  result <- fit_results[[ratio_name]]
  fits <- list("Raw" = result$raw_fit, "BF-adjusted" = result$bf_fit)
  row_index <- row_index + 1L
  runtime_rows[[row_index]] <- data.frame(
    m_ratio = result$m_ratio,
    m_size = result$m_size,
    raw_seconds = result$raw_seconds,
    bf_seconds = result$bf_seconds,
    total_seconds = result$raw_seconds + result$bf_seconds,
    stringsAsFactors = FALSE
  )
  current_warnings <- c(
    if (length(result$raw_warnings) > 0L) {
      paste0("Raw: ", result$raw_warnings)
    } else {
      character()
    },
    if (length(result$bf_warnings) > 0L) {
      paste0("BF-adjusted: ", result$bf_warnings)
    } else {
      character()
    }
  )
  warning_rows[[row_index]] <- if (length(current_warnings) == 0L) {
    data.frame(
      m_ratio = numeric(),
      m_size = integer(),
      warning = character(),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      m_ratio = result$m_ratio,
      m_size = result$m_size,
      warning = current_warnings,
      stringsAsFactors = FALSE
    )
  }
  fit_summary_bundle[[ratio_name]] <- list(
    m_ratio = result$m_ratio,
    m_size = result$m_size,
    null_indices = result$null_indices,
    raw = list(
      prior_weights = result$raw_fit$prior_weights,
      lfdr = result$raw_fit$lfdr
    ),
    bf_adjusted = list(
      prior_weights = result$bf_fit$prior_weights,
      lfdr = result$bf_fit$lfdr
    )
  )

  for (fit_stage in names(fits)) {
    current_fit <- fits[[fit_stage]]
    pi0_rows[[length(pi0_rows) + 1L]] <- make_pi0_row(
      current_fit,
      result$m_ratio,
      result$m_size,
      fit_stage
    )
    prior_rows[[length(prior_rows) + 1L]] <- make_prior_rows(
      current_fit,
      result$m_ratio,
      result$m_size,
      fit_stage
    )
    if (result$m_size > 0L) {
      group <- c(rep("target", J), rep("permuted_null", result$m_size))
      curve_rows[[length(curve_rows) + 1L]] <- summarize_small_m_curve(
        lfdr = current_fit$lfdr,
        group = group,
        pi0_merged = extract_pi0(current_fit),
        fit_stage = fit_stage,
        m_ratio = result$m_ratio,
        replicate_id = 1L,
        alpha_grid = alpha_grid
      )
    }
  }
}

pi0_summary <- do.call(rbind, pi0_rows)
prior_weights <- do.call(rbind, prior_rows)
alpha_curves <- do.call(rbind, curve_rows)
alpha_curves$is_single_draw <- TRUE
alpha005 <- alpha_curves[abs(alpha_curves$nominal_alpha - 0.05) < 1e-12, ]
runtime <- do.call(rbind, runtime_rows)
warnings <- do.call(rbind, warning_rows)

prior_checks <- aggregate(
  prior_weight ~ m_ratio + fit_stage,
  data = prior_weights,
  FUN = sum
)
nested_checks <- all(
  membership$null_index[membership$m_ratio == 0.05] %in%
    membership$null_index[membership$m_ratio == 0.10]
) && all(
  membership$null_index[membership$m_ratio == 0.10] %in%
    membership$null_index[membership$m_ratio == 0.20]
) && all(
  membership$null_index[membership$m_ratio == 0.20] %in%
    membership$null_index[membership$m_ratio == 1.00]
)
all_warning_messages <- c(
  null_capture$warnings,
  warnings$warning
)
all_warning_messages <- all_warning_messages[nzchar(all_warning_messages)]
validation <- data.frame(
  check = c(
    "clean_null_likelihood_dimensions",
    "clean_null_likelihood_finite",
    "expected_pi0_rows",
    "expected_alpha_curve_rows",
    "expected_alpha005_rows",
    "nested_membership",
    "prior_weights_sum_to_one",
    "merged_pi0_in_unit_interval",
    "overall_null_z_mean",
    "overall_null_z_sd",
    "null_z_off_diagonal_correlation",
    "zero_fit_warnings"
  ),
  passed = c(
    nrow(null_fit$L_matrix) == J && ncol(null_fit$L_matrix) == 52L,
    all(is.finite(null_fit$L_matrix)),
    nrow(pi0_summary) == 10L,
    nrow(alpha_curves) == 4L * 2L * length(alpha_grid),
    nrow(alpha005) == 8L,
    nested_checks,
    all(abs(prior_checks$prior_weight - 1) < 1e-8),
    all(pi0_summary$pi0_merged >= 0 & pi0_summary$pi0_merged <= 1),
    abs(z_diagnostics$overall$overall_mean) < 0.03,
    abs(z_diagnostics$overall$overall_sd - 1) < 0.03,
    z_diagnostics$overall$mean_absolute_off_diagonal_correlation < 0.03,
    length(all_warning_messages) == 0L
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  print(validation)
  if (length(all_warning_messages) > 0L) {
    message("Captured warning messages:")
    message(paste0("- ", all_warning_messages, collapse = "\n"))
  }
  stop("At least one clean working-model null validation failed.")
}

current_process_seconds <- proc.time()[["elapsed"]] - analysis_start
total_seconds <- if (checkpoint_reused) {
  current_process_seconds + null_likelihood_seconds
} else {
  current_process_seconds
}
configuration <- list(
  experiment = "small-M clean working-model null sensitivity",
  output_id = output_id,
  J = J,
  m_ratios = as.numeric(names(ratio_to_size)),
  m_sizes = unname(ratio_to_size),
  null_definition = "beta_hat_jt = SE_jt * Z_jt, iid Z_jt ~ N(0, 1)",
  null_seed = null_seed,
  subset_seed = subset_seed,
  n_simulation_draws = 1L,
  target_selection_method = "random_all_genes",
  target_selection_seed = 20260817L,
  alpha_grid = alpha_grid,
  num_cores = num_cores,
  null_likelihood_seconds = null_likelihood_seconds,
  checkpoint_reused = checkpoint_reused,
  current_process_seconds = unname(current_process_seconds),
  total_seconds = unname(total_seconds),
  r_version = R.version.string,
  package_versions = c(fashr = as.character(utils::packageVersion("fashr"))),
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

saveRDS(configuration, file.path(staging_directory, "configuration.rds"))
saveRDS(null_cache, file.path(staging_directory, "clean_null_fit.rds"))
saveRDS(fit_summary_bundle, file.path(staging_directory, "fit_summaries.rds"))
write_csv(input_provenance, file.path(staging_directory, "input_provenance.csv"))
write_csv(membership, file.path(staging_directory, "subset_membership.csv"))
write_csv(z_diagnostics$overall,
          file.path(staging_directory, "null_z_overall.csv"))
write_csv(z_diagnostics$time_summary,
          file.path(staging_directory, "null_z_by_time.csv"))
write_csv(z_diagnostics$correlation_long,
          file.path(staging_directory, "null_z_correlations.csv"))
write_csv(pi0_summary, file.path(staging_directory, "pi0_summary.csv"))
write_csv(prior_weights, file.path(staging_directory, "prior_weights.csv"))
write_csv(alpha005, file.path(staging_directory, "alpha005.csv"))
write_csv(alpha_curves, file.path(staging_directory, "alpha_curves.csv"))
write_csv(runtime, file.path(staging_directory, "runtime.csv"))
write_csv(warnings, file.path(staging_directory, "warnings.csv"))
write_csv(validation, file.path(staging_directory, "validation.csv"))

if (!file.rename(staging_directory, final_directory)) {
  stop("Could not atomically finalize the clean working-model null output.")
}
staging_complete <- TRUE
if (file.exists(checkpoint_path)) {
  unlink(checkpoint_path)
}

cat("\nSmall-M clean working-model null sensitivity completed.\n")
cat("Output: ", final_directory, "\n", sep = "")
cat("Null likelihood seconds: ",
    format(null_likelihood_seconds, digits = 8), "\n", sep = "")
cat("Total seconds: ", format(total_seconds, digits = 8), "\n", sep = "")
cat("Warnings: ", length(all_warning_messages), "\n", sep = "")
