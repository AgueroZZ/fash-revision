#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath("."))
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath("coderepo-local"))
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
cache_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "one_variant_per_gene_penalty_sensitivity"
)
required_files <- c(
  "configuration.rds",
  "full_fit_summary.csv",
  "full_penalty_comparison.csv",
  "seed_summary.csv",
  "seed12345_unpenalized_lfdr_comparison.csv"
)
stopifnot(
  dir.exists(cache_directory),
  all(file.exists(file.path(cache_directory, required_files)))
)

configuration <- readRDS(file.path(cache_directory, "configuration.rds"))
full_fit_summary <- utils::read.csv(
  file.path(cache_directory, "full_fit_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
full_penalty_comparison <- utils::read.csv(
  file.path(cache_directory, "full_penalty_comparison.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
seed_summary <- utils::read.csv(
  file.path(cache_directory, "seed_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
seed_specific <- utils::read.csv(
  file.path(cache_directory, "seed12345_unpenalized_lfdr_comparison.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  identical(configuration$seeds, c(12345L, 22345L, 32345L, 42345L, 52345L)),
  identical(configuration$n_full_units, 1009173L),
  identical(configuration$n_genes, 6362L),
  identical(configuration$n_audited_likelihood_rows, 31327L),
  identical(configuration$penalty_10_null_pseudo_units, 9L),
  configuration$likelihood_audit_max_absolute_difference == 0,
  identical(configuration$full_penalized_settings$penalty, 10),
  identical(configuration$full_unpenalized_settings$penalty, 1),
  nrow(full_fit_summary) == 3L,
  identical(
    full_fit_summary$Fit,
    c("Raw, penalty = 10", "Raw, penalty = 1", "BF-adjusted")
  ),
  identical(full_fit_summary$pair_discoveries, c(43860L, 43950L, 9205L)),
  identical(full_fit_summary$gene_discoveries, c(3258L, 3262L, 1177L)),
  nrow(full_penalty_comparison) == 1L,
  full_penalty_comparison$mean_absolute_lfdr_difference < 0.00125,
  full_penalty_comparison$mean_absolute_lfdr_difference > 0.00124,
  nrow(seed_summary) == 5L,
  identical(seed_summary$seed, configuration$seeds),
  all(seed_summary$n_units == 6362L),
  range(seed_summary$thinned_penalty1_pi0)[1] < 0.375,
  range(seed_summary$thinned_penalty1_pi0)[2] > 0.539,
  range(seed_summary$thinning_penalty1_lfdr_mae)[1] > 0.052,
  range(seed_summary$thinning_penalty1_lfdr_mae)[2] > 0.111,
  max(seed_summary$bf_penalty_prior_total_variation) == 0,
  max(seed_summary$bf_penalty_lfdr_max_absolute_difference) == 0,
  nrow(seed_specific) == 6362L,
  all(seed_specific$seed == 12345L),
  all(seed_specific$fit_stage == "Raw, penalty = 1"),
  !anyDuplicated(seed_specific$pair_key)
)

cat("One-variant-per-gene penalty-sensitivity tests passed.\n")
