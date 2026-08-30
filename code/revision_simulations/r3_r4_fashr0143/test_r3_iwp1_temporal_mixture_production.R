#!/usr/bin/env Rscript

# Regression tests for the canonical IWP1 temporal-mixture production contract.

find_workflowr_root <- function() {
  if (file.exists(file.path(
    "code", "revision_simulations", "shared", "simulation_functions.R"
  ))) return(".")
  if (file.exists(file.path(
    "coderepo-local", "code", "revision_simulations", "shared",
    "simulation_functions.R"
  ))) return("coderepo-local")
  stop("Could not locate the workflowr repository root.", call. = FALSE)
}

expect_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared",
  "r3_iwp1_temporal_mixture_contract.R"
))

contract <- r3_iwp1_temporal_mixture_contract(1272L)
expected_probabilities <- c(early = 0.29, middle = 0.42, late = 0.29)
expected_location_counts <- c(early = 369L, middle = 534L, late = 369L)
expected_truth_group_counts <- c(
  `early / switch` = 185L,
  `early / non-switch` = 184L,
  `middle / switch` = 267L,
  `middle / non-switch` = 267L,
  `late / switch` = 185L,
  `late / non-switch` = 184L
)
expected_switch_counts <- c(switch = 637L, `non-switch` = 635L)

expect_true(
  identical(contract$temporal_category_probs, expected_probabilities),
  "The canonical temporal-category probabilities changed."
)
expect_true(
  identical(contract$temporal_category_counts, expected_location_counts),
  "The canonical temporal-category integer counts changed."
)
expect_true(
  identical(contract$truth_group_counts, expected_truth_group_counts),
  "The canonical six-cell integer counts changed."
)
expect_true(
  identical(contract$switch_status_counts, expected_switch_counts),
  "The canonical marginal Switch counts changed."
)
expect_true(
  sum(contract$truth_group_counts) == 1272L,
  "The canonical truth-group counts do not sum to 1,272."
)

production_dir <- file.path(
  workflowr_root,
  "code", "revision_simulations", "r3_r4_fashr0143"
)
formal_entry_text <- paste(readLines(file.path(
  production_dir, "24_r3_full_universe_iwp1_mixture_fashr0143.R"
)), collapse = "\n")
formal_core_text <- paste(readLines(file.path(
  production_dir, "24_r3_full_universe_iwp1_mixture_core.R"
)), collapse = "\n")
formal_slurm_text <- paste(readLines(file.path(
  production_dir, "24_r3_full_universe_iwp1_mixture_fashr0143.sbatch"
)), collapse = "\n")
expected_formal_result_id <- paste0(
  "r3_real_genotype_one_per_gene_J6362_",
  "matched_functional_open_middle_3_12_center_aligned_",
  "iwp1_geometry_mixture_relative_location_clearance_",
  "full_universe_paired_posterior_fashr0143_pilot5"
)
expect_true(
  grepl("0.29,0.42,0.29", formal_entry_text, fixed = TRUE),
  "The formal entry point does not request the canonical mixture."
)
expect_true(
  grepl("iwp1_geometry_mixture_", formal_entry_text, fixed = TRUE) &&
    grepl(expected_formal_result_id, formal_slurm_text, fixed = TRUE),
  "The formal production result identity is inconsistent."
)
expect_true(
  grepl(
    "user-specified temporal-category probabilities",
    formal_core_text,
    fixed = TRUE
  ) &&
    grepl(
      "r3_iwp1_temporal_mixture_contract.R",
      formal_core_text,
      fixed = TRUE
    ),
  "The formal core does not validate the canonical mixture contract."
)
expect_true(
  !grepl("equal_cells", formal_entry_text, fixed = TRUE) &&
    !grepl("equal_cells", formal_slurm_text, fixed = TRUE),
  "The new formal production chain still advertises equal cells."
)

ideal_dir <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "r3_ideal_gaussian_measurement"
)
ideal_runner_text <- paste(readLines(file.path(
  ideal_dir, "run_r3_ideal_gaussian_measurement_iwp1_mixture.R"
)), collapse = "\n")
ideal_slurm_text <- paste(readLines(file.path(
  ideal_dir, "25_r3_ideal_gaussian_measurement_iwp1_mixture.sbatch"
)), collapse = "\n")
expected_ideal_result_id <- paste0(
  "r3_ideal_gaussian_known_t_adjusted_se_matched_truth_",
  "open_middle_3_12_center_aligned_iwp1_geometry_mixture_",
  "full_universe_fashr0143_pilot5"
)
expect_true(
  grepl(expected_ideal_result_id, ideal_slurm_text, fixed = TRUE) &&
    grepl("iwp1_geometry_mixture_", ideal_runner_text, fixed = TRUE),
  "The ideal-Gaussian production result identity is inconsistent."
)
expect_true(
  grepl(
    "r3_iwp1_temporal_mixture_contract.R",
    ideal_runner_text,
    fixed = TRUE
  ) &&
    grepl(
      "TEMPORAL_MIXTURE_CONTRACT$truth_group_counts",
      ideal_runner_text,
      fixed = TRUE
    ),
  "The ideal-Gaussian runner does not use the common mixture contract."
)
expect_true(
  grepl("all_alpha$alpha >= 0.05 - 1e-12", ideal_runner_text, fixed = TRUE) &&
    grepl("nrow(middle_curve) != 62L", ideal_runner_text, fixed = TRUE),
  "The ideal-Gaussian runner lacks the production alpha-boundary contract."
)
expect_true(
  !grepl("equal_cells", ideal_runner_text, fixed = TRUE) &&
    !grepl("equal_cells", ideal_slurm_text, fixed = TRUE),
  "The new ideal-Gaussian production chain still advertises equal cells."
)

cat("Canonical IWP1 temporal-mixture production contract tests passed.\n")
