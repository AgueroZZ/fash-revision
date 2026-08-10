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
original_working_directory <- getwd()
on.exit(setwd(original_working_directory), add = TRUE)
setwd(workflowr_root)

source(file.path(
  "code",
  "revision_simulations",
  "internal",
  "one_variant_per_gene_refit",
  "reporting.R"
))

stopifnot(
  nrow(aggregate_summary) == 20L,
  nrow(selection_audit) == 10L,
  all(selection_audit$`Unique genes` == 6362L),
  identical(unique_variant_range, c(6333L, 6351L)),
  identical(repeated_assignment_range, c(11L, 29L)),
  nrow(validation_summary) == 3L,
  max(validation_summary$likelihood_max_absolute_difference) < 1.5e-13,
  isTRUE(all.equal(full_raw_pi0, 0.428990301465861)),
  isTRUE(all.equal(full_bf_pi0, 0.938153319599315)),
  isTRUE(all.equal(raw_pi0_range, c(0.470505310346353, 0.559068881644855))),
  isTRUE(all.equal(bf_pi0_range, c(0.934768940584722, 0.946243319710783))),
  nrow(fit_range_table) == 2L,
  identical(fit_range_table$Fit, c("Raw", "BF-adjusted")),
  nrow(full_bf_discovery_summary) == 1L,
  identical(full_bf_pair_discoveries, 9205L),
  identical(full_bf_gene_discoveries, 1177L),
  nrow(variant_count_by_gene) == 6362L,
  sum(variant_count_by_gene$n_tested_variants) == 1009173L,
  identical(variant_count_min, 2L),
  isTRUE(all.equal(variant_count_q1, 111)),
  isTRUE(all.equal(variant_count_median, 150)),
  isTRUE(all.equal(variant_count_mean, 1009173 / 6362)),
  isTRUE(all.equal(variant_count_q3, 196)),
  identical(variant_count_max, 1859L),
  nrow(variant_count_summary_table) == 6L,
  identical(bf_thinned_trained_gene_range, c(42L, 70L)),
  identical(bf_thinned_trained_gene_range_width, 28L),
  identical(bf_thinned_gene_median, 58.5),
  identical(bf_thinned_gene_iqr, c(50.5, 65.25)),
  nrow(thinned_gene_discovery_table) == 10L,
  nrow(specific_seed_table) == 2L,
  min(bf_spearman_range) > 0.9978,
  max(bf_lfdr_mae_range) < 0.009,
  min(bf_jaccard_range) > 0.86,
  inherits(plot_variants_per_gene_histogram(), "ggplot"),
  inherits(plot_specific_seed_lfdr(), "ggplot"),
  inherits(plot_pi0_stability(), "ggplot"),
  inherits(plot_lfdr_mae(), "ggplot"),
  inherits(plot_discovery_jaccard(), "ggplot"),
  inherits(plot_thinned_gene_discoveries(), "ggplot")
)

cat("One-variant-per-gene reporting tests passed.\n")
