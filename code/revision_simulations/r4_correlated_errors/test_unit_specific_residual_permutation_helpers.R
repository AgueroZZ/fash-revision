#!/usr/bin/env Rscript

# Focused contracts for R4 unit-specific residual-permutation summaries.

source("code/revision_simulations/internal/covariance_estimation/donor_null_permutation_helpers.R")
source("code/revision_simulations/r4_correlated_errors/unit_specific_residual_permutation_helpers.R")

pair_keys <- sprintf("ENSG%06d_rs%03d", 1:10, 1:10)
selected <- select_random_r4_units(pair_keys, n_units = 4L, seed = 20260824L)
selected_repeat <- select_random_r4_units(pair_keys, n_units = 4L, seed = 20260824L)
if (!identical(selected, selected_repeat) || nrow(selected) != 4L ||
    anyDuplicated(selected$pair_key)) {
  stop("Random unit selection is not deterministic and unique.")
}

donor_ids <- paste0("donor", 1:8)
patterns <- setNames(c(rep("1111", 6), rep("1011", 2)), donor_ids)
map <- make_r4_donor_map(donor_ids, patterns, seed = 20260825L)
if (nrow(map) != length(donor_ids) || anyDuplicated(map$source_donor) ||
    !identical(unname(patterns[map$source_donor]),
               unname(patterns[map$target_donor]))) {
  stop("The R4 donor map does not preserve observation patterns.")
}

set.seed(7)
draws <- matrix(stats::rnorm(400L * 16L), nrow = 400L, ncol = 16L)
summary <- summarize_r4_null_draws(draws)
if (!identical(dim(summary$correlation), c(16L, 16L)) ||
    !identical(summary$variogram$lag, 1:15) ||
    any(!is.finite(summary$variogram$semivariogram))) {
  stop("The null correlation summary is malformed.")
}
identity_summary <- summarize_r4_correlation(diag(4L), time_grid = 0:3)
if (!identical(identity_summary$variogram$lag, 1:3) ||
    any(identity_summary$variogram$semivariogram != 1)) {
  stop("The direct correlation summary is malformed.")
}
intervals <- bootstrap_r4_null_variogram_intervals(
  draws, n_bootstrap = 100L, seed = 20260826L
)
if (nrow(intervals) != 15L || any(intervals$lower > intervals$upper) ||
    any(!is.finite(as.matrix(intervals[, c("lower", "upper")]))) ||
    !all(intervals$n_bootstrap == 100L)) {
  stop("The bootstrap variogram intervals are malformed.")
}
long <- correlation_matrix_to_long_r4(summary$correlation, "Unit 1")
if (nrow(long) != 256L || !identical(unique(long$unit_label), "Unit 1")) {
  stop("The heatmap long form is malformed.")
}

message("PASS: unit-specific residual-permutation helper contracts.")
