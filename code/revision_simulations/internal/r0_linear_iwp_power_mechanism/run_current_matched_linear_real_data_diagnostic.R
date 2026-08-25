#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr project root.")
}

summarize_calls <- function(selected, keys, method, adjustment) {
  selected_keys <- keys[selected]
  genes <- sub("_.*$", "", selected_keys)
  variants <- sub("^[^_]*_", "", selected_keys)
  data.frame(
    method = method,
    adjustment = adjustment,
    pair_count = length(selected_keys),
    gene_count = length(unique(genes)),
    variant_count = length(unique(variants)),
    stringsAsFactors = FALSE
  )
}

summarize_overlap <- function(first_selected,
                              second_selected,
                              first_label,
                              second_label) {
  data.frame(
    first_method = first_label,
    second_method = second_label,
    first_count = sum(first_selected),
    second_count = sum(second_selected),
    intersection_count = sum(first_selected & second_selected),
    first_only_count = sum(first_selected & !second_selected),
    second_only_count = sum(!first_selected & second_selected),
    union_count = sum(first_selected | second_selected),
    jaccard = sum(first_selected & second_selected) /
      sum(first_selected | second_selected),
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
  "fash_linear_real_data_ablation", "fash_linear_real_data_helpers.R"
))

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

input_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "fash_linear_real_data_ablation"
)
statistics_path <- file.path(input_directory, "sufficient_statistics.rds")
retained_analysis_path <- file.path(input_directory, "analysis_cache.rds")
retained_linear_path <- file.path(input_directory, "linear_fit_bf.rds")

statistics_cache <- readRDS(statistics_path)
retained_analysis <- readRDS(retained_analysis_path)
retained_linear <- readRDS(retained_linear_path)
statistics <- statistics_cache$statistics
keys <- as.character(statistics$unit_id)

if (!inherits(retained_linear, "profiled_linear_fash") ||
    !identical(keys, as.character(retained_linear$log_marginal$unit_id)) ||
    nrow(retained_analysis$lfdr_scatter_all) != length(keys)) {
  stop("The retained artifacts are not aligned with the sufficient statistics.")
}

grid <- default_revision_grid()
pred_step <- 1
penalty <- 10L
alpha <- 0.05

start_time <- proc.time()[["elapsed"]]
matched_raw <- fit_linear_mixture_fash_from_stats(
  statistics = statistics,
  grid = grid,
  pred_step = pred_step,
  penalty = penalty,
  statistic_time_span = 15,
  n_time = 16L
)
raw_elapsed_seconds <- proc.time()[["elapsed"]] - start_time
validate_linear_mixture_fash(
  matched_raw,
  expected_grid = grid,
  expected_pred_step = pred_step,
  expected_penalty = penalty
)

bf_start_time <- proc.time()[["elapsed"]]
matched_bf <- BF_update_linear_mixture_fash(matched_raw)
bf_elapsed_seconds <- proc.time()[["elapsed"]] - bf_start_time
validate_linear_mixture_fash(
  matched_bf,
  expected_grid = grid,
  expected_pred_step = pred_step,
  expected_penalty = penalty
)

matched_raw_selected <- rep(FALSE, length(keys))
matched_raw_selected[select_cumulative_lfdr_calls_linear(
  matched_raw$lfdr,
  alpha
)] <- TRUE
matched_bf_selected <- rep(FALSE, length(keys))
matched_bf_selected[select_cumulative_lfdr_calls_linear(
  matched_bf$lfdr,
  alpha
)] <- TRUE
retained_profiled_selected <-
  retained_analysis$lfdr_scatter_all$discovery_status %in%
    c("FASH-linear only", "Both")
retained_iwp_selected <-
  retained_analysis$lfdr_scatter_all$discovery_status %in%
    c("Current FASH only", "Both")

discovery_summary <- rbind(
  summarize_calls(
    matched_raw_selected,
    keys,
    "Current matched-mixture FASH-linear",
    "Raw"
  ),
  summarize_calls(
    matched_bf_selected,
    keys,
    "Current matched-mixture FASH-linear",
    "BF-adjusted"
  ),
  summarize_calls(
    retained_profiled_selected,
    keys,
    "Retained profiled-slab FASH-linear",
    "BF-adjusted"
  ),
  summarize_calls(
    retained_iwp_selected,
    keys,
    "Retained FASH-IWP1",
    "BF-adjusted"
  )
)

overlap_summary <- rbind(
  summarize_overlap(
    matched_bf_selected,
    retained_iwp_selected,
    "Current matched-mixture FASH-linear",
    "Retained FASH-IWP1"
  ),
  summarize_overlap(
    matched_bf_selected,
    retained_profiled_selected,
    "Current matched-mixture FASH-linear",
    "Retained profiled-slab FASH-linear"
  )
)

prior_summary <- rbind(
  transform(
    extract_linear_mixture_prior_table(matched_raw, 0L, "Raw"),
    adjustment = "Raw"
  ),
  transform(
    extract_linear_mixture_prior_table(matched_bf, 0L, "BF"),
    adjustment = "BF-adjusted"
  )
)
null_weight_summary <- data.frame(
  method = c(
    "Current matched-mixture FASH-linear",
    "Current matched-mixture FASH-linear",
    "Retained profiled-slab FASH-linear",
    "Retained FASH-IWP1"
  ),
  adjustment = c("Raw", "BF-adjusted", "BF-adjusted", "BF-adjusted"),
  null_weight = c(
    constant_component_prior_weight(matched_raw),
    constant_component_prior_weight(matched_bf),
    retained_linear$prior_weights$prior_weight[[1L]],
    0.938153319599315
  ),
  stringsAsFactors = FALSE
)

compact_raw <- compact_linear_mixture_fash(matched_raw)
compact_bf <- compact_linear_mixture_fash(matched_bf)
rm(matched_raw, matched_bf, statistics)
invisible(gc())

validation <- data.frame(
  check = c(
    "retained comparator provenance is profiled slab",
    "current comparator uses matched grid and penalty",
    "unit alignment is exact",
    "retained discovery counts reproduce",
    "current summaries are finite",
    "one-job thread cap"
  ),
  passed = c(
    inherits(retained_linear, "profiled_linear_fash"),
    identical(compact_bf$psd_grid, grid) &&
      identical(compact_bf$settings$pred_step, pred_step) &&
      identical(as.integer(compact_bf$settings$penalty), penalty),
    identical(compact_bf$unit_ids, keys),
    sum(retained_profiled_selected) == 15865L &&
      sum(retained_iwp_selected) == 9205L,
    all(is.finite(discovery_summary$pair_count)) &&
      all(is.finite(null_weight_summary$null_weight)) &&
      all(is.finite(overlap_summary$jaccard)),
    Sys.getenv("VECLIB_MAXIMUM_THREADS") == "1"
  ),
  stringsAsFactors = FALSE
)
if (any(!validation$passed)) {
  print(validation)
  stop("The current matched-linear real-data diagnostic failed validation.")
}

output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r0_linear_iwp_power_mechanism_current_matched_real_data"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
analysis_cache <- list(
  configuration = list(
    scope = "non-overwriting current-definition real-data diagnostic",
    grid = grid,
    pred_step = pred_step,
    penalty = penalty,
    alpha = alpha,
    sufficient_statistics_path = normalizePath(statistics_path),
    retained_analysis_path = normalizePath(retained_analysis_path),
    retained_linear_path = normalizePath(retained_linear_path)
  ),
  discovery_summary = discovery_summary,
  overlap_summary = overlap_summary,
  null_weight_summary = null_weight_summary,
  prior_summary = prior_summary,
  compact_raw = compact_raw,
  compact_bf = compact_bf,
  runtime = data.frame(
    stage = c("matched raw fit", "BF update"),
    elapsed_seconds = c(raw_elapsed_seconds, bf_elapsed_seconds),
    stringsAsFactors = FALSE
  ),
  validation = validation
)
saveRDS(analysis_cache, file.path(output_directory, "analysis_cache.rds"))
utils::write.csv(
  discovery_summary,
  file.path(output_directory, "discovery_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  overlap_summary,
  file.path(output_directory, "overlap_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  null_weight_summary,
  file.path(output_directory, "null_weight_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  prior_summary,
  file.path(output_directory, "prior_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  analysis_cache$runtime,
  file.path(output_directory, "runtime.csv"),
  row.names = FALSE
)
utils::write.csv(
  validation,
  file.path(output_directory, "validation.csv"),
  row.names = FALSE
)

print(discovery_summary)
print(overlap_summary)
print(null_weight_summary)
print(analysis_cache$runtime)
print(validation)
message("Saved output to: ", normalizePath(output_directory))
