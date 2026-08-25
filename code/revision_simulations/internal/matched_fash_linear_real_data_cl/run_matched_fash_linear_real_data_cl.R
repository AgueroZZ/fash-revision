#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

workflowr_root <- if (file.exists("_workflowr.yml")) {
  "."
} else if (file.exists("coderepo-local/_workflowr.yml")) {
  "coderepo-local"
} else {
  stop("Run this script from the workflowr root or its parent.")
}

shared_path <- file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
)
historical_helper_path <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "fash_linear_real_data_ablation", "fash_linear_real_data_helpers.R"
)
matched_helper_path <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "matched_fash_linear_real_data", "matched_fash_linear_helpers.R"
)
overlap_helper_path <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "matched_fash_linear_real_data", "iwp_overlap_proportion_helpers.R"
)
runner_path <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "matched_fash_linear_real_data_cl",
  "run_matched_fash_linear_real_data_cl.R"
)
source(shared_path)
source(historical_helper_path)
source(matched_helper_path)
source(overlap_helper_path)

analysis_id <- "matched_fash_linear_real_data_cl_fashr_0_1_43"
output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal", analysis_id
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

fit_directory <- file.path(workflowr_root, "output", "dynamic_eQTL_real")
cl_iwp1_raw_path <- file.path(fit_directory, "fash_fit1_all_CL.RData")
cl_iwp2_raw_path <- file.path(fit_directory, "fash_fit2_all_CL.RData")
cl_statistics_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "fash_linear_real_data_ablation_cl", "sufficient_statistics.rds"
)
strober_linear_path <- file.path(
  workflowr_root,
  "data", "dynamic_eQTL_real", "strober_linear",
  "linear_dynamic_eqtls_5_pc.txt"
)
strober_quadratic_path <- file.path(
  workflowr_root,
  "data", "dynamic_eQTL_real", "strober_nonlinear",
  "non_linear_dynamic_eqtls_5_pc.txt"
)
input_paths <- c(
  cl_iwp1_raw = cl_iwp1_raw_path,
  cl_iwp2_raw = cl_iwp2_raw_path,
  cl_sufficient_statistics = cl_statistics_path,
  strober_linear = strober_linear_path,
  strober_quadratic = strober_quadratic_path,
  shared_functions = shared_path,
  historical_linear_helpers = historical_helper_path,
  matched_helpers = matched_helper_path,
  overlap_helpers = overlap_helper_path,
  runner = runner_path
)
if (any(!file.exists(input_paths))) {
  stop("At least one required CL-PC input is missing.")
}

analysis_cache_path <- file.path(output_directory, "analysis_cache.rds")
run_status_path <- file.path(output_directory, "run_status.rds")
linear_raw_path <- file.path(output_directory, "linear_fit_raw.rds")
linear_bf_path <- file.path(output_directory, "linear_fit_bf.rds")
iwp1_bf_path <- file.path(output_directory, "cl_iwp1_bf_adjustment.rds")
iwp2_bf_path <- file.path(output_directory, "cl_iwp2_bf_adjustment.rds")

analysis_started <- proc.time()[["elapsed"]]
atomic_save_rds_matched(list(
  analysis_id = analysis_id,
  status = "in_progress",
  started_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
), run_status_path)

compact_iwp_adjustment <- function(fit, method, package_provenance) {
  list(
    method = method,
    package_version = package_provenance$version,
    package_sha = package_provenance$remote_sha,
    settings = fit$settings,
    psd_grid = as.numeric(fit$psd_grid),
    prior_weights = fit$prior_weights,
    estimated_pi0 = extract_null_weight_matched(fit),
    lfdr = as.numeric(fit$lfdr)
  )
}

validate_iwp_fit <- function(fit,
                             expected_order,
                             expected_keys,
                             expected_grid,
                             expected_settings) {
  if (!inherits(fit, "fash") ||
      length(fit$lfdr) != length(expected_keys) ||
      !identical(names(fit$fash_data$data_list), expected_keys) ||
      !identical(as.numeric(fit$psd_grid), expected_grid) ||
      !identical(as.integer(fit$settings$num_basis), 20L) ||
      !identical(as.numeric(fit$settings$betaprec), 0) ||
      !identical(as.integer(fit$settings$order), as.integer(expected_order)) ||
      !identical(as.numeric(fit$settings$pred_step), 1) ||
      !identical(as.character(fit$settings$likelihood), "gaussian") ||
      !identical(as.integer(fit$settings$penalty), 10L) ||
      !identical(
        as.numeric(fit$settings$pred_step),
        expected_settings$pred_step
      )) {
    stop("A CL-PC IWP fit is not aligned with the approved settings.")
  }
  invisible(TRUE)
}

conditional_alternative_difference <- function(raw_fit, bf_fit, grid) {
  raw_weights <- expand_grid_prior_weights(raw_fit$prior_weights, grid)
  bf_weights <- expand_grid_prior_weights(bf_fit$prior_weights, grid)
  raw_conditional <- raw_weights[-1L] / sum(raw_weights[-1L])
  bf_conditional <- bf_weights[-1L] / sum(bf_weights[-1L])
  max(abs(raw_conditional - bf_conditional))
}

message("[1/8] Validating package, CL-PC statistics, and immutable inputs.")
package_provenance <- extract_package_provenance_matched()
expected_package_sha <- "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
if (!identical(package_provenance$version, "0.1.43") ||
    !identical(package_provenance$remote_sha, expected_package_sha)) {
  stop("The installed fashr package is not the approved 0.1.43 build.")
}
input_provenance <- make_file_provenance_matched(input_paths)
input_md5_before <- stats::setNames(
  input_provenance$md5,
  input_provenance$label
)

statistics_cache <- readRDS(cl_statistics_path)
required_statistics_fields <- c(
  "analysis_id", "cl_raw_md5", "expected_time", "scale_time", "ridge",
  "statistics"
)
if (!is.list(statistics_cache) ||
    !all(required_statistics_fields %in% names(statistics_cache)) ||
    !identical(statistics_cache$expected_time, 0:15) ||
    !isTRUE(statistics_cache$scale_time) ||
    !identical(statistics_cache$ridge, 1e-10) ||
    !identical(
      as.character(statistics_cache$cl_raw_md5),
      unname(tools::md5sum(cl_iwp1_raw_path))
    ) ||
    !is.data.frame(statistics_cache$statistics) ||
    nrow(statistics_cache$statistics) != 1009173L) {
  stop("The retained CL-PC sufficient-statistic cache is invalid.")
}
statistics <- statistics_cache$statistics

message("[2/8] Loading CL-PC IWP1 and validating the matched design.")
cl_iwp1_raw <- load_exact_object_matched(cl_iwp1_raw_path, "fash_fit1")
pair_keys <- names(cl_iwp1_raw$fash_data$data_list)
matched_grid <- as.numeric(cl_iwp1_raw$psd_grid)
matched_settings <- list(
  likelihood = as.character(cl_iwp1_raw$settings$likelihood),
  num_basis_iwp1 = as.integer(cl_iwp1_raw$settings$num_basis),
  betaprec = as.numeric(cl_iwp1_raw$settings$betaprec),
  order_iwp1 = as.integer(cl_iwp1_raw$settings$order),
  pred_step = as.numeric(cl_iwp1_raw$settings$pred_step),
  penalty = as.integer(cl_iwp1_raw$settings$penalty),
  grid = matched_grid,
  time_grid = 0:15,
  statistic_time_span = 15,
  alpha = 0.05
)
if (!identical(matched_settings$likelihood, "gaussian") ||
    !identical(matched_settings$num_basis_iwp1, 20L) ||
    !identical(matched_settings$betaprec, 0) ||
    !identical(matched_settings$order_iwp1, 1L) ||
    !identical(matched_settings$pred_step, 1) ||
    !identical(matched_settings$penalty, 10L) ||
    length(matched_grid) != 52L ||
    !identical(matched_grid, default_revision_grid()) ||
    length(pair_keys) != 1009173L ||
    anyDuplicated(pair_keys) ||
    !identical(statistics$unit_id, pair_keys)) {
  stop("The CL-PC IWP1 settings, grid, or pair order are invalid.")
}
validate_iwp_fit(
  cl_iwp1_raw,
  expected_order = 1L,
  expected_keys = pair_keys,
  expected_grid = matched_grid,
  expected_settings = matched_settings
)

validation_indices <- c(1L, 12345L, 250000L, 500000L, 750000L, 1009173L)
validation_data <- lapply(validation_indices, function(original_index) {
  dataset <- cl_iwp1_raw$fash_data$data_list[[original_index]]
  data.frame(
    time = as.numeric(dataset$x),
    beta = as.numeric(dataset$y),
    SE = as.numeric(cl_iwp1_raw$fash_data$S[[original_index]]),
    stringsAsFactors = FALSE
  )
})
names(validation_data) <- pair_keys[validation_indices]
recomputed_statistics <- compute_linear_sufficient_statistics_block(
  validation_data,
  expected_time = matched_settings$time_grid,
  scale_time = TRUE,
  ridge = statistics_cache$ridge
)
cached_statistics <- statistics[validation_indices, , drop = FALSE]
statistic_columns <- setdiff(names(cached_statistics), "unit_id")
statistics_match <- identical(
  recomputed_statistics$unit_id,
  cached_statistics$unit_id
) && max(abs(
  as.matrix(recomputed_statistics[statistic_columns]) -
    as.matrix(cached_statistics[statistic_columns])
)) < 1e-8

direct_validation_data <- lapply(validation_data, function(dataset) {
  list(x = dataset$time, y = dataset$beta, sd = dataset$SE)
})
validation_grid <- matched_grid[c(1L, 5L, 15L, 25L, 35L, 45L, 52L)]
direct_validation_likelihood <- compute_linear_mixture_log_likelihood(
  datasets = direct_validation_data,
  grid = validation_grid,
  pred_step = matched_settings$pred_step
)
cached_validation_likelihood <- compute_linear_mixture_log_likelihood_from_stats(
  statistics = cached_statistics,
  grid = validation_grid,
  pred_step = matched_settings$pred_step,
  statistic_time_span = matched_settings$statistic_time_span,
  n_time = length(matched_settings$time_grid)
)
likelihood_match <- identical(
  dimnames(direct_validation_likelihood),
  dimnames(cached_validation_likelihood)
) && max(abs(
  direct_validation_likelihood - cached_validation_likelihood
)) < 1e-8

pair_table <- parse_linear_pair_keys(pair_keys)
iwp1_raw_lfdr <- as.numeric(cl_iwp1_raw$lfdr)
iwp1_raw_pi0 <- extract_null_weight_matched(cl_iwp1_raw)

message("[3/8] Applying corrected fashr 0.1.43 BF adjustment to CL-PC IWP1.")
iwp1_bf_started <- proc.time()[["elapsed"]]
cl_iwp1_bf <- fashr::BF_update(cl_iwp1_raw)
validate_iwp_fit(
  cl_iwp1_bf,
  expected_order = 1L,
  expected_keys = pair_keys,
  expected_grid = matched_grid,
  expected_settings = matched_settings
)
iwp1_conditional_difference <- conditional_alternative_difference(
  cl_iwp1_raw,
  cl_iwp1_bf,
  matched_grid
)
iwp1_bf_lfdr <- as.numeric(cl_iwp1_bf$lfdr)
iwp1_bf_pi0 <- extract_null_weight_matched(cl_iwp1_bf)
atomic_save_rds_matched(
  compact_iwp_adjustment(cl_iwp1_bf, "FASH-IWP1", package_provenance),
  iwp1_bf_path
)
iwp1_bf_elapsed <- proc.time()[["elapsed"]] - iwp1_bf_started
rm(
  cl_iwp1_raw,
  cl_iwp1_bf,
  validation_data,
  direct_validation_data,
  recomputed_statistics,
  direct_validation_likelihood,
  cached_validation_likelihood
)
invisible(gc())

message("[4/8] Fitting the matched CL-PC FASH-linear mixture.")
linear_fit_started <- proc.time()[["elapsed"]]
linear_raw <- fit_linear_mixture_fash_from_stats(
  statistics = statistics,
  grid = matched_grid,
  pred_step = matched_settings$pred_step,
  penalty = matched_settings$penalty,
  statistic_time_span = matched_settings$statistic_time_span,
  n_time = length(matched_settings$time_grid)
)
linear_raw$settings$betaprec <- matched_settings$betaprec
linear_raw$settings$likelihood <- matched_settings$likelihood
linear_raw$settings$reference_iwp_num_basis <- matched_settings$num_basis_iwp1
linear_raw$settings$reference_iwp_order <- matched_settings$order_iwp1
validate_linear_mixture_fash(
  linear_raw,
  expected_grid = matched_grid,
  expected_pred_step = matched_settings$pred_step,
  expected_penalty = matched_settings$penalty
)
atomic_save_rds_matched(linear_raw, linear_raw_path)

linear_bf <- BF_update_linear_mixture_fash(linear_raw)
validate_linear_mixture_fash(
  linear_bf,
  expected_grid = matched_grid,
  expected_pred_step = matched_settings$pred_step,
  expected_penalty = matched_settings$penalty
)
atomic_save_rds_matched(linear_bf, linear_bf_path)
linear_fit_elapsed <- proc.time()[["elapsed"]] - linear_fit_started

raw_linear_weights <- expand_grid_prior_weights(
  linear_raw$prior_weights,
  matched_grid
)
bf_linear_weights <- expand_grid_prior_weights(
  linear_bf$prior_weights,
  matched_grid
)
linear_conditional_difference <- conditional_alternative_difference(
  linear_raw,
  linear_bf,
  matched_grid
)
linear_raw_lfdr <- as.numeric(linear_raw$lfdr)
linear_bf_lfdr <- as.numeric(linear_bf$lfdr)
linear_raw_pi0 <- extract_null_weight_matched(linear_raw)
linear_bf_pi0 <- extract_null_weight_matched(linear_bf)

message("[5/8] Applying corrected fashr 0.1.43 BF adjustment to CL-PC IWP2.")
iwp2_bf_started <- proc.time()[["elapsed"]]
cl_iwp2_raw <- load_exact_object_matched(cl_iwp2_raw_path, "fash_fit2")
validate_iwp_fit(
  cl_iwp2_raw,
  expected_order = 2L,
  expected_keys = pair_keys,
  expected_grid = matched_grid,
  expected_settings = matched_settings
)
cl_iwp2_bf <- fashr::BF_update(cl_iwp2_raw)
validate_iwp_fit(
  cl_iwp2_bf,
  expected_order = 2L,
  expected_keys = pair_keys,
  expected_grid = matched_grid,
  expected_settings = matched_settings
)
iwp2_conditional_difference <- conditional_alternative_difference(
  cl_iwp2_raw,
  cl_iwp2_bf,
  matched_grid
)
iwp2_bf_lfdr <- as.numeric(cl_iwp2_bf$lfdr)
iwp2_bf_pi0 <- extract_null_weight_matched(cl_iwp2_bf)
atomic_save_rds_matched(
  compact_iwp_adjustment(cl_iwp2_bf, "FASH-IWP2", package_provenance),
  iwp2_bf_path
)
iwp2_bf_elapsed <- proc.time()[["elapsed"]] - iwp2_bf_started
rm(cl_iwp2_raw, cl_iwp2_bf)
invisible(gc())

message("[6/8] Constructing discoveries, directional summaries, and Venn sets.")
discoveries <- list(
  iwp1_raw = summarize_discoveries_matched(
    "FASH-IWP1", "Raw", iwp1_raw_lfdr, pair_table,
    iwp1_raw_pi0, matched_settings$alpha
  ),
  iwp1_bf = summarize_discoveries_matched(
    "FASH-IWP1", "BF-adjusted", iwp1_bf_lfdr, pair_table,
    iwp1_bf_pi0, matched_settings$alpha
  ),
  linear_raw = summarize_discoveries_matched(
    "FASH-linear", "Raw", linear_raw_lfdr, pair_table,
    linear_raw_pi0, matched_settings$alpha
  ),
  linear_bf = summarize_discoveries_matched(
    "FASH-linear", "BF-adjusted", linear_bf_lfdr, pair_table,
    linear_bf_pi0, matched_settings$alpha
  )
)
discovery_counts <- do.call(rbind, lapply(discoveries, `[[`, "summary"))
rownames(discovery_counts) <- NULL

strober_linear <- utils::read.delim(
  strober_linear_path,
  stringsAsFactors = FALSE
)
strober_quadratic <- utils::read.delim(
  strober_quadratic_path,
  stringsAsFactors = FALSE
)
for (table_name in c("strober_linear", "strober_quadratic")) {
  table <- get(table_name)
  required_columns <- c("ensamble_id", "rs_id", "eFDR")
  if (!all(required_columns %in% names(table))) {
    stop(table_name, " is missing required columns.")
  }
  table$key <- paste0(table$ensamble_id, "_", table$rs_id)
  table <- table[
    is.finite(table$eFDR) & table$eFDR <= matched_settings$alpha,
    ,
    drop = FALSE
  ]
  assign(table_name, table)
}

venn_sets <- build_four_method_venn_sets_matched(
  strober_quadratic = strober_quadratic,
  strober_linear = strober_linear,
  fash_iwp1 = discoveries$iwp1_bf$table,
  fash_linear = discoveries$linear_bf$table
)
validate_four_method_venn_sets_matched(venn_sets)
venn_region_counts <- do.call(rbind, lapply(names(venn_sets), function(unit) {
  exclusive_set_regions_matched(venn_sets[[unit]], unit)
}))
pairwise_overlap <- do.call(rbind, lapply(names(venn_sets), function(unit) {
  pairwise_overlap_matched(venn_sets[[unit]], unit)
}))
rownames(venn_region_counts) <- NULL
rownames(pairwise_overlap) <- NULL

linear_sets <- lapply(venn_sets, function(sets) sets[["FASH-linear BF"]])
iwp1_sets <- lapply(venn_sets, function(sets) sets[["FASH-IWP1 BF"]])
iwp1_linear <- directional_overlap_summary_matched(
  reference_sets = linear_sets,
  comparison_sets = iwp1_sets
)
iwp2_indices <- select_cumulative_lfdr_calls_linear(
  iwp2_bf_lfdr,
  matched_settings$alpha
)
iwp2_table <- pair_table[iwp2_indices, , drop = FALSE]
iwp2_sets <- list(
  `Gene-variant pairs` = unique(iwp2_table$key),
  Genes = unique(iwp2_table$gene_id),
  Variants = unique(iwp2_table$variant_id)
)
iwp2_linear <- directional_overlap_summary_matched(
  reference_sets = linear_sets,
  comparison_sets = iwp2_sets
)

message("[7/8] Building validation, provenance, and reporting summaries.")
matched_settings_table <- data.frame(
  setting = c(
    "Tested units", "PC covariates", "Time grid", "Likelihood",
    "Null trajectory", "Alternative family", "Scale parameter",
    "Grid components", "Positive grid range", "pred_step", "betaprec",
    "penalty", "IWP basis functions", "FDR rule"
  ),
  fash_iwp1 = c(
    format(length(pair_keys), big.mark = ","),
    "Five cell-line-collapsed PCs", "0:15", "Gaussian",
    "Unrestricted constant", "IWP1 Gaussian process mixture", "psd",
    length(matched_grid),
    sprintf("%.6f to %.1f", min(matched_grid[-1L]), max(matched_grid)),
    matched_settings$pred_step, matched_settings$betaprec,
    matched_settings$penalty, matched_settings$num_basis_iwp1,
    "Cumulative lfdr at 0.05"
  ),
  fash_linear = c(
    format(length(pair_keys), big.mark = ","),
    "Five cell-line-collapsed PCs", "0:15", "Gaussian",
    "Unrestricted constant", "Zero-centered Gaussian slope mixture",
    "SD(beta * pred_step)", length(matched_grid),
    sprintf("%.6f to %.1f", min(matched_grid[-1L]), max(matched_grid)),
    matched_settings$pred_step, matched_settings$betaprec,
    matched_settings$penalty, "Not used: exact rank-one linear kernel",
    "Cumulative lfdr at 0.05"
  ),
  alignment = c(
    "Identical", "Identical", "Identical", "Identical", "Identical",
    "Matched scale; different covariance geometry",
    "Equivalent one-step SD", "Identical", "Identical", "Identical",
    "Matched; not used by linear kernel", "Identical", "Not applicable",
    "Identical"
  ),
  stringsAsFactors = FALSE
)

prior_weights <- rbind(
  data.frame(
    method = "FASH-linear", adjustment = "Raw",
    scale = matched_grid, prior_weight = raw_linear_weights,
    stringsAsFactors = FALSE
  ),
  data.frame(
    method = "FASH-linear", adjustment = "BF-adjusted",
    scale = matched_grid, prior_weight = bf_linear_weights,
    stringsAsFactors = FALSE
  )
)

immutable_input_labels <- c(
  "cl_iwp1_raw", "cl_iwp2_raw", "cl_sufficient_statistics",
  "strober_linear", "strober_quadratic"
)
input_md5_after <- unname(tools::md5sum(
  input_paths[immutable_input_labels]
))
input_immutable <- identical(
  input_md5_after,
  unname(input_md5_before[immutable_input_labels])
)
region_invariants <- all(vapply(names(venn_sets), function(unit) {
  sets <- venn_sets[[unit]]
  regions <- venn_region_counts[venn_region_counts$unit == unit, ]
  sum(regions$count) == length(unique(unlist(sets, use.names = FALSE)))
}, logical(1)))

validation <- data.frame(
  check = c(
    "approved_fashr_version_and_sha",
    "cl_pc_iwp1_settings_match_requested_hyperparameters",
    "grid_values_and_size_identical",
    "pair_count_and_order_identical",
    "retained_statistics_match_cl_pc_iwp_data",
    "direct_and_accelerated_likelihoods_match",
    "raw_linear_fit_valid",
    "bf_linear_fit_valid",
    "linear_bf_preserves_conditional_alternative_weights",
    "iwp1_bf_preserves_conditional_alternative_weights",
    "iwp2_bf_preserves_conditional_alternative_weights",
    "authoritative_inputs_immutable",
    "four_method_venn_sets_complete",
    "exclusive_venn_regions_partition_each_union",
    "single_process_single_thread_configuration"
  ),
  passed = c(
    TRUE,
    TRUE,
    identical(matched_grid, default_revision_grid()),
    identical(statistics$unit_id, pair_keys),
    statistics_match,
    likelihood_match,
    isTRUE(tryCatch({
      validate_linear_mixture_fash(linear_raw)
      TRUE
    }, error = function(condition) FALSE)),
    isTRUE(tryCatch({
      validate_linear_mixture_fash(linear_bf)
      TRUE
    }, error = function(condition) FALSE)),
    linear_conditional_difference < 1e-10,
    iwp1_conditional_difference < 1e-10,
    iwp2_conditional_difference < 1e-10,
    input_immutable,
    isTRUE(tryCatch({
      validate_four_method_venn_sets_matched(venn_sets)
      TRUE
    }, error = function(condition) FALSE)),
    region_invariants,
    all(Sys.getenv(c(
      "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
      "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"
    )) == "1")
  ),
  stringsAsFactors = FALSE
)
if (any(!validation$passed)) {
  print(validation)
  stop("At least one matched CL-PC validation failed.")
}

fit_provenance <- make_file_provenance_matched(c(
  linear_raw = linear_raw_path,
  linear_bf = linear_bf_path,
  cl_iwp1_bf = iwp1_bf_path,
  cl_iwp2_bf = iwp2_bf_path
))
runtime_summary <- data.frame(
  stage = c(
    "CL-PC IWP1 BF update",
    "Full CL-PC linear fit plus BF update",
    "CL-PC IWP2 BF update",
    "Complete runner"
  ),
  elapsed_seconds = c(
    iwp1_bf_elapsed,
    linear_fit_elapsed,
    iwp2_bf_elapsed,
    proc.time()[["elapsed"]] - analysis_started
  ),
  stringsAsFactors = FALSE
)
configuration <- list(
  analysis_id = analysis_id,
  scope = "Internal full real-data matched FASH-linear CL-PC rerun",
  package_version = package_provenance$version,
  package_sha = package_provenance$remote_sha,
  pc_correction = "Cell-line-collapsed PCs repeated across time",
  expected_pairs = length(pair_keys),
  time_grid = matched_settings$time_grid,
  likelihood = matched_settings$likelihood,
  betaprec = matched_settings$betaprec,
  reference_iwp_order = matched_settings$order_iwp1,
  reference_iwp_num_basis = matched_settings$num_basis_iwp1,
  pred_step = matched_settings$pred_step,
  penalty = matched_settings$penalty,
  grid = matched_grid,
  linear_scale_definition = "sd_linear_step = SD(beta * pred_step)",
  linear_prior = paste(
    "exact constant null plus a finite mixture of zero-centered",
    "Gaussian slopes on the matched one-step SD grid"
  ),
  alpha = matched_settings$alpha,
  strober_rule = "eFDR <= 0.05"
)

analysis_cache <- list(
  configuration = configuration,
  package_provenance = package_provenance,
  input_provenance = input_provenance,
  fit_provenance = fit_provenance,
  matched_settings_table = matched_settings_table,
  discovery_counts = discovery_counts,
  prior_weights = prior_weights,
  venn_sets = venn_sets,
  venn_region_counts = venn_region_counts,
  pairwise_overlap = pairwise_overlap,
  iwp1_linear = iwp1_linear,
  iwp2_linear = iwp2_linear,
  validation = validation,
  runtime_summary = runtime_summary
)

message("[8/8] Saving the versioned CL-PC reporting cache and exports.")
atomic_save_rds_matched(analysis_cache, analysis_cache_path)
exports <- list(
  matched_settings = matched_settings_table,
  discovery_counts = discovery_counts,
  prior_weights = prior_weights,
  venn_region_counts = venn_region_counts,
  pairwise_overlap = pairwise_overlap,
  iwp1_linear = iwp1_linear,
  iwp2_linear = iwp2_linear,
  validation = validation,
  input_provenance = input_provenance,
  fit_provenance = fit_provenance,
  runtime_summary = runtime_summary
)
for (label in names(exports)) {
  utils::write.csv(
    exports[[label]],
    file.path(output_directory, paste0(label, ".csv")),
    row.names = FALSE
  )
}

atomic_save_rds_matched(list(
  analysis_id = analysis_id,
  status = "complete",
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  analysis_cache_md5 = unname(tools::md5sum(analysis_cache_path)),
  fit_md5 = stats::setNames(fit_provenance$md5, fit_provenance$label),
  package_sha = package_provenance$remote_sha
), run_status_path)

message("Matched CL-PC FASH-linear rerun completed.")
print(discovery_counts)
print(iwp1_linear)
print(iwp2_linear)
print(validation)
