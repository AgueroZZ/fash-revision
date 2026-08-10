#!/usr/bin/env Rscript

# Validate the retained R5 reporting layer and its derived summaries.

find_workflowr_root <- function() {
  if (file.exists("analysis/index.Rmd")) {
    return(normalizePath("."))
  }
  if (file.exists("coderepo-local/analysis/index.Rmd")) {
    return(normalizePath("coderepo-local"))
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
setwd(workflowr_root)
source("code/revision_simulations/r5_balanced_thinning/reporting.R")

stopifnot(
  identical(expected_target, 20L),
  identical(nrow(variant_count_by_gene), 6362L),
  sum(variant_count_values) == 1009173L,
  identical(nrow(excluded_genes), 10L),
  !any(excluded_genes$full_data_discovered_gene),
  identical(nrow(seed_summary), 10L),
  all(seed_summary$warning_count == 0L),
  all(seed_summary$fit_stage == "BF-updated"),
  all(diff(cumulative_discoveries$cumulative_unique_genes) >= 0L),
  all(diff(cumulative_discoveries$cumulative_full_genes_recovered) >= 0L),
  identical(final_cumulative_unique, 1113L),
  identical(final_cumulative_recovered, 1075L),
  identical(final_additional_unique, 38L),
  isTRUE(all.equal(full_pi0, 0.938153319599315)),
  min(seed_summary$spearman_lfdr) > 0.9997,
  max(seed_summary$mean_absolute_lfdr_difference) < 0.004
)

plot_functions <- list(
  plot_variants_per_gene_histogram,
  plot_seed_12345_lfdr,
  plot_null_weight_stability,
  plot_paired_lfdr_agreement,
  plot_cumulative_discovered_genes
)
stopifnot(all(vapply(plot_functions, function(plot_function) {
  inherits(plot_function(), "ggplot")
}, logical(1))))

message("R5 reporting tests passed.")
