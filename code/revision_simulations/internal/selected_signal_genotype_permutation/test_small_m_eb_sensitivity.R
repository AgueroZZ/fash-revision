#!/usr/bin/env Rscript

# Deterministic tests for the small-M matched-null EB sensitivity helpers.

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

ratio_to_size <- c("0.2" = 4L, "0.4" = 8L, "0.6" = 12L)
subsets_a <- make_nested_null_subsets(
  n_null = 20L,
  ratio_to_size = ratio_to_size,
  n_replicates = 3L,
  seed = 17L
)
subsets_b <- make_nested_null_subsets(
  n_null = 20L,
  ratio_to_size = ratio_to_size,
  n_replicates = 3L,
  seed = 17L
)
stopifnot(identical(subsets_a, subsets_b))
for (replicate_id in 1:3) {
  replicate_rows <- subsets_a[subsets_a$replicate_id == replicate_id, ]
  set_02 <- replicate_rows$null_index[replicate_rows$m_ratio == 0.2]
  set_04 <- replicate_rows$null_index[replicate_rows$m_ratio == 0.4]
  set_06 <- replicate_rows$null_index[replicate_rows$m_ratio == 0.6]
  stopifnot(
    length(set_02) == 4L,
    length(set_04) == 8L,
    length(set_06) == 12L,
    all(set_02 %in% set_04),
    all(set_04 %in% set_06)
  )
}

synthetic_lfdr <- c(0.01, 0.04, 0.90, 0.95, 0.02, 0.80)
synthetic_group <- c(rep("target", 4L), rep("permuted_null", 2L))
synthetic_summary <- summarize_small_m_calibration(
  lfdr = synthetic_lfdr,
  group = synthetic_group,
  pi0_merged = 0.8,
  fit_stage = "Synthetic",
  m_ratio = 0.5,
  replicate_id = 1L,
  alpha = 0.05
)
stopifnot(
  synthetic_summary$target_calls == 2L,
  synthetic_summary$permuted_null_calls == 1L,
  synthetic_summary$merged_calls == 3L,
  abs(synthetic_summary$permuted_null_call_rate - 0.5) < 1e-15,
  abs(synthetic_summary$pi0_target_unbounded - 0.7) < 1e-15,
  abs(synthetic_summary$target_pi0_plugin_fdr - 0.7) < 1e-15,
  abs(synthetic_summary$merged_pi0_plugin_fdr - 0.8) < 1e-15,
  abs(synthetic_summary$known_null_discovery_fraction - 1 / 3) < 1e-15
)

zero_summary <- summarize_small_m_calibration(
  lfdr = synthetic_lfdr,
  group = synthetic_group,
  pi0_merged = 0.8,
  fit_stage = "Synthetic",
  m_ratio = 0.5,
  replicate_id = 1L,
  alpha = 0.001
)
stopifnot(
  zero_summary$merged_calls == 0L,
  zero_summary$permuted_null_call_rate == 0,
  is.na(zero_summary$target_pi0_plugin_fdr),
  is.na(zero_summary$merged_pi0_plugin_fdr),
  is.na(zero_summary$known_null_discovery_fraction)
)

full_output_id <- paste0(
  "all_gene_random_variant_signal_stripped_residual_block_permutation_",
  "selection20260817_seed20260811"
)
full_output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  full_output_id
)
full_bundle <- readRDS(file.path(
  full_output_directory,
  "merged_fash_fit.rds"
))
saved_calibration <- utils::read.csv(file.path(
  full_output_directory,
  "calibration_diagnostics.csv"
), stringsAsFactors = FALSE)
bf_group <- rep(
  c("target", "permuted_null"),
  each = full_bundle$configuration$n_target_units
)
bf_summary <- summarize_small_m_calibration(
  lfdr = full_bundle$bf_adjusted_fit$lfdr,
  group = bf_group,
  pi0_merged = extract_pi0(full_bundle$bf_adjusted_fit),
  fit_stage = "BF-adjusted",
  m_ratio = 1,
  replicate_id = 0L,
  alpha = 0.05
)
saved_bf <- saved_calibration[
  saved_calibration$fit_stage == "BF-adjusted",
  ,
  drop = FALSE
]
stopifnot(
  nrow(saved_bf) == 1L,
  bf_summary$target_calls == saved_bf$target_calls,
  bf_summary$permuted_null_calls == saved_bf$permuted_null_calls,
  bf_summary$merged_calls == saved_bf$total_calls,
  abs(
    bf_summary$permuted_null_call_rate - saved_bf$permuted_null_call_rate
  ) < 1e-15,
  abs(
    bf_summary$target_pi0_plugin_fdr -
      saved_bf$post_selection_fdr_target_from_pi0
  ) < 1e-14,
  abs(
    bf_summary$merged_pi0_plugin_fdr -
      saved_bf$scaled_fdr_merged_from_estimated_pi0
  ) < 1e-14
)

cat("Small-M matched-null EB sensitivity tests passed.\n")
