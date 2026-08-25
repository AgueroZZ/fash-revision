#!/usr/bin/env Rscript

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
source("code/revision_simulations/r5_one_variant_per_gene/reporting.R")

stopifnot(
  identical(expected_seeds, seq(12345L, by = 10000L, length.out = 100L)),
  identical(nrow(variant_count_by_gene), 6362L),
  sum(variant_count_values) == 1009173L,
  identical(nrow(seed_summary), 100L),
  all(seed_summary$warning_count == 0L),
  all(seed_summary$fit_stage == "BF-updated"),
  all(seed_summary$n_selected_units == 6362L),
  identical(nrow(gene_discovery_frequency), 6362L),
  isTRUE(all.equal(unname(thinned_discovery_range), c(45, 86))),
  identical(ever_discovered, 1041L),
  identical(ever_discovered_full, 943L),
  identical(ever_discovered_additional, 98L),
  identical(repeated_variant_assignment_range, c(11L, 29L)),
  identical(maximum_variant_multiplicity, 3L),
  identical(shared_configuration$cache_id, expected_cache_id),
  identical(
    shared_configuration$package_provenance,
    expected_fashr_provenance
  ),
  isTRUE(all.equal(full_pi0, 0.938159265061590)),
  min(seed_summary$thinned_pi0) > 0.929,
  max(seed_summary$thinned_pi0) < 0.947,
  min(seed_summary$spearman_lfdr) > 0.996,
  max(seed_summary$mean_absolute_lfdr_difference) < 0.012,
  min(seed_summary$fdr_call_jaccard) > 0.81
)

table_output <- render_scrollable_table(
  data.frame(Item = "Example", Value = "1"),
  caption = "Example table",
  align = "ll",
  minimum_width = "400px"
)
stopifnot(
  inherits(table_output, "knit_asis"),
  grepl("<table", as.character(table_output), fixed = TRUE),
  grepl('class="r5-table-scroll"', as.character(table_output), fixed = TRUE),
  grepl('style="min-width:400px;"', as.character(table_output), fixed = TRUE),
  !grepl("[1]", as.character(table_output), fixed = TRUE)
)

plot_functions <- list(
  plot_variants_per_gene_histogram,
  plot_seed_12345_lfdr,
  plot_null_weight_stability,
  plot_paired_lfdr_agreement,
  plot_discovery_robustness
)
stopifnot(all(vapply(plot_functions, function(plot_function) {
  inherits(plot_function(), "ggplot")
}, logical(1))))

message("R5 one-variant-per-gene reporting tests passed.")
