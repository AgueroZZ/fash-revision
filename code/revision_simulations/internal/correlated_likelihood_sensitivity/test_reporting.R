#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath("."))
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return(normalizePath("coderepo-local"))
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
original_working_directory <- getwd()
on.exit(setwd(original_working_directory), add = TRUE)
setwd(workflowr_root)

source(file.path(
  "code",
  "revision_simulations",
  "internal",
  "correlated_likelihood_sensitivity",
  "reporting.R"
))

stopifnot(
  identical(configuration$analysis_id, "correlated_likelihood_sensitivity"),
  identical(configuration$thinning_seed, 20260811L),
  identical(configuration$n_units, 6362L),
  identical(configuration$n_genes, 6362L),
  identical(configuration$n_time, 16L),
  nrow(pair_metadata) == 6362L,
  anyDuplicated(pair_metadata$pair_key) == 0L,
  anyDuplicated(pair_metadata$gene_id) == 0L,
  nrow(method_stage_summary) == 6L,
  isTRUE(all.equal(raw_diagonal_pi0, 0.482599786139328)),
  isTRUE(all.equal(raw_c1_pi0, 0.990269970125036)),
  identical(raw_c2_pi0, 1),
  isTRUE(all.equal(bf_diagonal_pi0, 0.936183590066017)),
  isTRUE(all.equal(bf_c1_pi0, 0.999371266897202)),
  identical(raw_diagonal_calls, 268L),
  identical(bf_diagonal_calls, 57L),
  identical(c1_estimating_pairs, 3436L),
  identical(c2_estimating_pairs, 2522L),
  identical(bf_update_status$bf_update_available, c(TRUE, TRUE, FALSE)),
  identical(bf_update_status$raw_alternative_prior_mass[3], 0),
  all(is.na(lfdr_wide$C2_bf)),
  nrow(lfdr_long) == 31810L,
  nrow(lfdr_pairwise_metrics) == 4L,
  nrow(prior_pairwise_metrics) == 4L,
  nrow(top_discrepancy_display) == 20L,
  max(identity_validation$row_centered_likelihood_maximum_difference) < 3e-13,
  max(bf_rebuild_validation$lfdr_maximum_difference) < 1e-10,
  min(vapply(correlations, function(C) {
    min(eigen(C, symmetric = TRUE, only.values = TRUE)$values)
  }, numeric(1))) > 0.38,
  inherits(plot_correlation_heatmaps(), "ggplot"),
  inherits(plot_prior_weights(), "ggplot"),
  inherits(plot_lfdr_comparisons(), "ggplot"),
  inherits(plot_lfdr_difference_distributions(), "ggplot"),
  inherits(plot_discovery_counts(), "ggplot")
)

cat("Correlated-likelihood reporting tests passed.\n")
