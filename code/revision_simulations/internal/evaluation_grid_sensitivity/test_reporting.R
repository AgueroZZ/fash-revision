#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  reporting_file <- file.path(
    "code", "revision_simulations", "internal",
    "evaluation_grid_sensitivity", "reporting.R"
  )
  if (file.exists(reporting_file)) return(".")
  if (file.exists(file.path("coderepo-local", reporting_file))) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
old_working_directory <- getwd()
on.exit(setwd(old_working_directory), add = TRUE)
setwd(workflowr_root)
source(file.path(
  "code", "revision_simulations", "internal",
  "evaluation_grid_sensitivity", "reporting.R"
))

stopifnot(
  identical(population_summary$category, category_order),
  identical(population_summary$n_pairs, unname(expected_pair_counts)),
  identical(population_summary$n_genes, unname(expected_gene_counts)),
  nrow(comparison_summary) == 16L,
  nrow(pair_lfsr) == 5930L,
  max(grid_comparisons$maximum_absolute_change) == 0.01,
  max(mc_3000_vs_5000$maximum_absolute_change) < 0.0122,
  min(mc_3000_vs_5000$fraction_within_0p01) > 0.9989,
  abs(sampling_runtime_seconds - 170.8) < 0.1
)

message("Current R4 numerical-sensitivity reporting tests passed.")
