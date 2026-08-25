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
helper_path <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "matched_fash_linear_real_data", "matched_fash_linear_helpers.R"
)
runner_path <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "matched_fash_linear_real_data", "run_matched_fash_linear_real_data.R"
)
source(shared_path)
source(historical_helper_path)
source(helper_path)

analysis_id <- "matched_fash_linear_real_data_fashr_0_1_43"
output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal", analysis_id
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

fit_directory <- file.path(workflowr_root, "output", "dynamic_eQTL_real")
current_raw_path <- file.path(fit_directory, "fash_fit1_all.RData")
current_bf_path <- file.path(fit_directory, "fash_fit1_update.RData")
sufficient_statistics_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "fash_linear_real_data_ablation_mixture_predstep1_penalty10",
  "sufficient_statistics.rds"
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
  current_iwp1_raw = current_raw_path,
  current_iwp1_bf = current_bf_path,
  sufficient_statistics = sufficient_statistics_path,
  strober_linear = strober_linear_path,
  strober_quadratic = strober_quadratic_path,
  shared_functions = shared_path,
  historical_linear_helpers = historical_helper_path,
  matched_helpers = helper_path,
  runner = runner_path
)
if (any(!file.exists(input_paths))) {
  stop("At least one required input is missing.")
}

analysis_started <- proc.time()[["elapsed"]]
run_status_path <- file.path(output_directory, "run_status.rds")
analysis_cache_path <- file.path(output_directory, "analysis_cache.rds")
linear_raw_path <- file.path(output_directory, "linear_fit_raw.rds")
linear_bf_path <- file.path(output_directory, "linear_fit_bf.rds")
atomic_save_rds_matched(list(
  analysis_id = analysis_id,
  status = "in_progress",
  started_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
), run_status_path)

message("[1/7] Validating package and retained inputs.")
package_provenance <- extract_package_provenance_matched()
expected_package_sha <- "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
if (!identical(package_provenance$version, "0.1.43") ||
    !identical(package_provenance$remote_sha, expected_package_sha)) {
  stop("The installed fashr package is not the approved 0.1.43 build.")
}
input_provenance <- make_file_provenance_matched(input_paths)
input_md5_before <- stats::setNames(input_provenance$md5, input_provenance$label)

statistics_cache <- readRDS(sufficient_statistics_path)
required_statistics_fields <- c(
  "analysis_id", "dataset_md5", "expected_time", "scale_time", "ridge",
  "statistics"
)
if (!is.list(statistics_cache) ||
    !all(required_statistics_fields %in% names(statistics_cache)) ||
    !identical(statistics_cache$expected_time, 0:15) ||
    !isTRUE(statistics_cache$scale_time) ||
    !identical(statistics_cache$ridge, 1e-10) ||
    !is.data.frame(statistics_cache$statistics) ||
    nrow(statistics_cache$statistics) != 1009173L) {
  stop("The retained sufficient-statistic cache is invalid.")
}
statistics <- statistics_cache$statistics

message("[2/7] Loading authoritative IWP1 settings and validating statistics.")
current_raw <- load_exact_object_matched(current_raw_path, "fash_fit1")
pair_keys <- names(current_raw$fash_data$data_list)
matched_grid <- as.numeric(current_raw$psd_grid)
matched_settings <- list(
  likelihood = as.character(current_raw$settings$likelihood),
  num_basis_iwp1 = as.integer(current_raw$settings$num_basis),
  betaprec = as.numeric(current_raw$settings$betaprec),
  order_iwp1 = as.integer(current_raw$settings$order),
  pred_step = as.numeric(current_raw$settings$pred_step),
  penalty = as.integer(current_raw$settings$penalty),
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
  stop("The IWP1 settings, grid, or unit order do not match the approved design.")
}

validation_indices <- c(1L, 12345L, 250000L, 500000L, 750000L, 1009173L)
validation_data <- lapply(seq_along(validation_indices), function(index) {
  original_index <- validation_indices[index]
  dataset <- current_raw$fash_data$data_list[[original_index]]
  data.frame(
    time = as.numeric(dataset$x),
    beta = as.numeric(dataset$y),
    SE = as.numeric(current_raw$fash_data$S[[original_index]]),
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
current_raw_lfdr <- as.numeric(current_raw$lfdr)
current_raw_pi0 <- extract_null_weight_matched(current_raw)
rm(current_raw, validation_data, direct_validation_data)
invisible(gc())

message("[3/7] Recomputing the full matched FASH-linear raw and BF fits.")
fit_started <- proc.time()[["elapsed"]]
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
fit_elapsed <- proc.time()[["elapsed"]] - fit_started

raw_linear_weights <- expand_grid_prior_weights(linear_raw$prior_weights, matched_grid)
bf_linear_weights <- expand_grid_prior_weights(linear_bf$prior_weights, matched_grid)
raw_conditional_alternative <- raw_linear_weights[-1L] /
  sum(raw_linear_weights[-1L])
bf_conditional_alternative <- bf_linear_weights[-1L] /
  sum(bf_linear_weights[-1L])
conditional_alternative_difference <- max(abs(
  raw_conditional_alternative - bf_conditional_alternative
))

linear_raw_lfdr <- as.numeric(linear_raw$lfdr)
linear_bf_lfdr <- as.numeric(linear_bf$lfdr)
linear_raw_pi0 <- extract_null_weight_matched(linear_raw)
linear_bf_pi0 <- extract_null_weight_matched(linear_bf)

message("[4/7] Loading BF-adjusted IWP1 and defining all discovery sets.")
current_bf <- load_exact_object_matched(current_bf_path, "fash_fit1_update")
if (!identical(names(current_bf$fash_data$data_list), pair_keys) ||
    !identical(as.numeric(current_bf$psd_grid), matched_grid) ||
    !identical(
      as.numeric(current_bf$settings$betaprec),
      matched_settings$betaprec
    ) ||
    !identical(
      as.numeric(current_bf$settings$pred_step),
      matched_settings$pred_step
    ) ||
    !identical(
      as.integer(current_bf$settings$penalty),
      matched_settings$penalty
    )) {
  stop("The BF-adjusted IWP1 fit is not aligned with the raw reference fit.")
}
current_bf_lfdr <- as.numeric(current_bf$lfdr)
current_bf_pi0 <- extract_null_weight_matched(current_bf)
rm(current_bf)
invisible(gc())

discoveries <- list(
  iwp1_raw = summarize_discoveries_matched(
    "FASH-IWP1", "Raw", current_raw_lfdr, pair_table,
    current_raw_pi0, matched_settings$alpha
  ),
  iwp1_bf = summarize_discoveries_matched(
    "FASH-IWP1", "BF-adjusted", current_bf_lfdr, pair_table,
    current_bf_pi0, matched_settings$alpha
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

strober_linear <- utils::read.delim(strober_linear_path, stringsAsFactors = FALSE)
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
  table <- table[is.finite(table$eFDR) & table$eFDR <= matched_settings$alpha, ]
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

message("[5/7] Building matched-setting, prior, and validation summaries.")
matched_settings_table <- data.frame(
  setting = c(
    "Tested units", "Time grid", "Likelihood", "Null trajectory",
    "Alternative family", "Scale parameter", "Grid components",
    "Positive grid range", "pred_step", "betaprec", "penalty",
    "IWP basis functions", "FDR rule"
  ),
  fash_iwp1 = c(
    format(length(pair_keys), big.mark = ","), "0:15", "Gaussian",
    "Unrestricted constant", "IWP1 Gaussian process mixture", "psd",
    length(matched_grid),
    sprintf("%.6f to %.1f", min(matched_grid[-1L]), max(matched_grid)),
    matched_settings$pred_step, matched_settings$betaprec,
    matched_settings$penalty, matched_settings$num_basis_iwp1,
    "Cumulative lfdr at 0.05"
  ),
  fash_linear = c(
    format(length(pair_keys), big.mark = ","), "0:15", "Gaussian",
    "Unrestricted constant", "Zero-centered Gaussian slope mixture",
    "SD(beta * pred_step)", length(matched_grid),
    sprintf("%.6f to %.1f", min(matched_grid[-1L]), max(matched_grid)),
    matched_settings$pred_step, matched_settings$betaprec,
    matched_settings$penalty, "Not used: exact rank-one linear kernel",
    "Cumulative lfdr at 0.05"
  ),
  alignment = c(
    rep("Identical", 4L), "Matched scale; different covariance geometry",
    "Equivalent one-step SD", rep("Identical", 5L), "Not applicable",
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
  "current_iwp1_raw", "current_iwp1_bf", "sufficient_statistics",
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
    "iwp1_settings_match_requested_hyperparameters",
    "grid_values_and_size_identical",
    "pair_count_and_order_identical",
    "retained_statistics_match_current_iwp_data",
    "direct_and_accelerated_likelihoods_match",
    "raw_linear_fit_valid",
    "bf_linear_fit_valid",
    "bf_preserves_conditional_alternative_weights",
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
      validate_linear_mixture_fash(linear_raw); TRUE
    }, error = function(condition) FALSE)),
    isTRUE(tryCatch({
      validate_linear_mixture_fash(linear_bf); TRUE
    }, error = function(condition) FALSE)),
    conditional_alternative_difference < 1e-10,
    input_immutable,
    isTRUE(tryCatch({
      validate_four_method_venn_sets_matched(venn_sets); TRUE
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
  stop("At least one matched FASH-linear validation failed.")
}

fit_provenance <- make_file_provenance_matched(c(
  linear_raw = linear_raw_path,
  linear_bf = linear_bf_path
))
runtime_summary <- data.frame(
  stage = c("Full linear fit plus BF update", "Complete runner"),
  elapsed_seconds = c(
    fit_elapsed,
    proc.time()[["elapsed"]] - analysis_started
  ),
  stringsAsFactors = FALSE
)
configuration <- list(
  analysis_id = analysis_id,
  scope = "Internal full real-data matched FASH-linear rerun",
  package_version = package_provenance$version,
  package_sha = package_provenance$remote_sha,
  pc_correction = "Time-specific PCs",
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
  validation = validation,
  runtime_summary = runtime_summary
)

message("[6/7] Saving the versioned reporting cache and CSV exports.")
atomic_save_rds_matched(analysis_cache, analysis_cache_path)
utils::write.csv(
  matched_settings_table,
  file.path(output_directory, "matched_settings.csv"),
  row.names = FALSE
)
utils::write.csv(
  discovery_counts,
  file.path(output_directory, "discovery_counts.csv"),
  row.names = FALSE
)
utils::write.csv(
  prior_weights,
  file.path(output_directory, "prior_weights.csv"),
  row.names = FALSE
)
utils::write.csv(
  venn_region_counts,
  file.path(output_directory, "venn_region_counts.csv"),
  row.names = FALSE
)
utils::write.csv(
  pairwise_overlap,
  file.path(output_directory, "pairwise_overlap.csv"),
  row.names = FALSE
)
utils::write.csv(
  validation,
  file.path(output_directory, "validation.csv"),
  row.names = FALSE
)
utils::write.csv(
  input_provenance,
  file.path(output_directory, "input_provenance.csv"),
  row.names = FALSE
)
utils::write.csv(
  fit_provenance,
  file.path(output_directory, "fit_provenance.csv"),
  row.names = FALSE
)
utils::write.csv(
  runtime_summary,
  file.path(output_directory, "runtime_summary.csv"),
  row.names = FALSE
)

atomic_save_rds_matched(list(
  analysis_id = analysis_id,
  status = "complete",
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  analysis_cache_md5 = unname(tools::md5sum(analysis_cache_path)),
  fit_md5 = stats::setNames(fit_provenance$md5, fit_provenance$label),
  package_sha = package_provenance$remote_sha
), run_status_path)

message("[7/7] Matched FASH-linear rerun completed.")
print(discovery_counts)
print(validation)
