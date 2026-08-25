#!/usr/bin/env Rscript

# Static and generative tests for the formal R3/R4 fashr 0.1.43 contracts.
# This file intentionally does not fit a formal simulation replicate.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
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
  "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "r4_correlated_errors",
  "real_genotype_r1_helpers.R"
))

expected_version <- "0.1.43"
expected_sha <- "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
expect_true(requireNamespace("fashr", quietly = TRUE), "fashr is unavailable.")
description <- utils::packageDescription("fashr")
observed_sha <- if (is.null(description$RemoteSha)) {
  NA_character_
} else {
  as.character(description$RemoteSha)
}
expect_true(
  identical(as.character(utils::packageVersion("fashr")), expected_version),
  "Unexpected fashr version."
)
expect_true(identical(observed_sha, expected_sha), "Unexpected fashr RemoteSha.")

evaluation_grid <- seq(0, 15, by = 0.1)
open_middle <- temporal_middle_membership(
  evaluation_grid,
  middle_window = c(3, 12),
  middle_boundary = "open"
)
legacy_middle <- temporal_middle_membership(
  evaluation_grid,
  middle_window = c(4, 11),
  middle_boundary = "closed"
)
expect_true(
  isTRUE(all.equal(
    evaluation_grid[open_middle],
    seq(3.1, 11.9, by = 0.1),
    tolerance = 1e-12,
    check.attributes = FALSE
  )),
  "The open-Middle grid is not exactly 3.1 through 11.9."
)
expect_true(
  isTRUE(all.equal(
    evaluation_grid[legacy_middle],
    seq(4, 11, by = 0.1),
    tolerance = 1e-12,
    check.attributes = FALSE
  )),
  "The legacy closed-Middle grid is not exactly 4.0 through 11.0."
)
expect_true(
  sum(open_middle) == 89L && sum(legacy_middle) == 71L,
  "The open and legacy Middle grid sizes are unexpected."
)

open_middle_effects <- simulate_matched_functional_effect_set(
  n_variants = 60L,
  truth_mechanism = "raised_cosine",
  evaluation_grid = evaluation_grid,
  non_switch_min_range_fraction = 0.10,
  seed = 20260823L,
  scenario = "r3_open_middle_contract_test",
  middle_window = c(3, 12),
  middle_boundary = "open"
)
recomputed_open_middle_functionals <- evaluate_temporal_functionals(
  open_middle_effects$beta_evaluation,
  smooth_var = evaluation_grid,
  switch_threshold = 0.25,
  middle_window = c(3, 12),
  middle_boundary = "open"
)
expect_true(
  isTRUE(all.equal(
    open_middle_effects$true_functionals,
    recomputed_open_middle_functionals,
    tolerance = 1e-12
  )),
  "Open-Middle truth generation and functional evaluation are inconsistent."
)
expect_true(
  identical(open_middle_effects$settings$middle_window, c(3, 12)) &&
    identical(open_middle_effects$settings$middle_boundary, "open"),
  "The matched truth artifact does not record the open-Middle definition."
)
expected_center_ranges <- list(
  early = c(0.5, 1.5),
  middle = c(4.5, 10.5),
  late = c(13.5, 14.5)
)
expect_true(
  identical(
    open_middle_effects$settings$raised_cosine_center_ranges,
    expected_center_ranges
  ),
  paste(
    "The raised-cosine truth artifact does not record",
    "support-contained center ranges."
  )
)
dynamic_indices <- which(
  open_middle_effects$unit_info$effect_class == "dynamic_bspline"
)
primary_centers_are_aligned <- vapply(
  dynamic_indices,
  function(index) {
    time_group <- open_middle_effects$unit_info$time_group[[index]]
    center_range <- expected_center_ranges[[time_group]]
    primary_center <- open_middle_effects$unit_info$peak_centers[[index]][[1L]]
    is.finite(primary_center) &&
      primary_center >= center_range[[1L]] &&
      primary_center <= center_range[[2L]]
  },
  logical(1)
)
expect_true(
  all(primary_centers_are_aligned),
  paste(
    "A raised-cosine primary peak center is outside its assigned",
    "support-contained range."
  )
)

genotype_cache_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "shared",
  "real_genotype_one_per_gene_J6362_pilot5", "genotype_samples.rds"
)
r1_cache_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "mc", "r1_r2_fashr0143"
)
matrix_cache_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "real_data",
  "r4_null_like_top500_full_correlations_fashr0143"
)
for (path in c(
  genotype_cache_path,
  file.path(r1_cache_dir, "complete.flag"),
  file.path(r1_cache_dir, "replicates", "r1", "seed_12345.rds"),
  file.path(matrix_cache_dir, "real_data_correlation_analysis.rds"),
  file.path(matrix_cache_dir, "simulation_correlation_matrices.rds")
)) {
  expect_true(file.exists(path), paste("Required test input is missing:", path))
}

seed_list <- c(12345L, 22345L, 32345L, 42345L, 52345L)
genotype_cache <- validate_r4_real_genotype_cache(
  readRDS(genotype_cache_path),
  seed_list = seed_list,
  J = 6362L,
  n_donors = 19L
)
seed_inputs <- prepare_r4_real_genotype_r1_seed(
  genotype_cache,
  seed = 12345L
)
r1_reference <- readRDS(file.path(
  r1_cache_dir, "replicates", "r1", "seed_12345.rds"
))
expect_true(
  identical(
    seed_inputs$genotype_sample$selection$pair_key,
    r1_reference$selected_pair_keys
  ),
  "The R4 seed constructor does not use the formal R1 pair ordering."
)
expect_true(
  identical(
    seed_inputs$effect_sim$unit_info$effect_class,
    r1_reference$unit_info$effect_class
  ),
  "The R4 seed constructor does not reproduce the formal R1 truth classes."
)
expect_true(
  isTRUE(all.equal(
    seed_inputs$effect_sim$unit_info$genetic_main_effect,
    r1_reference$unit_info$genetic_main_effect,
    tolerance = 1e-12
  )),
  "The R4 seed constructor does not reproduce the formal R1 main effects."
)
expect_true(
  identical(seed_inputs$scenario, r1_reference$scenario),
  "The R4 scenario label does not match formal R1."
)

expected_counts <- exact_proportional_counts(
  6362L,
  c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
)
observed_counts <- table(factor(
  seed_inputs$effect_sim$unit_info$effect_class,
  levels = names(expected_counts)
))
expect_true(
  identical(as.integer(observed_counts), as.integer(expected_counts)),
  "The R4 truth-class counts are invalid."
)

test_correlation <- make_lag1_correlation(16L, 0.1)
candidate_expression <- make_r4_correlated_expression(
  seed_inputs,
  expression_error_correlation = test_correlation
)
maximum_pairing_difference <- validate_r4_condition_pair(
  seed_inputs,
  reference_expression = seed_inputs$independent_expression,
  candidate_expression = candidate_expression,
  correlation = test_correlation
)
expect_true(
  maximum_pairing_difference <= 1e-10,
  "The R4 correlated-error transformation is not exactly paired."
)

real_analysis <- readRDS(file.path(
  matrix_cache_dir,
  "real_data_correlation_analysis.rds"
))
expect_true(
  identical(
    real_analysis$configuration$fit_sha256,
    "7f0ca9ab0fbeab89a13c83d2a0fb7c24195f7b5a5835f209399cf0e359001f50"
  ),
  "The R4 real-data cache does not use the retained 0.1.43 fit."
)
expect_true(
  identical(
    real_analysis$configuration$package_provenance$version,
    expected_version
  ) && identical(
    real_analysis$configuration$package_provenance$remote_sha,
    expected_sha
  ),
  "The R4 real-data cache has unexpected package provenance."
)
expect_true(
  nrow(real_analysis$selected_units) == 500L &&
    length(unique(real_analysis$selected_units$gene_id)) == 500L &&
    real_analysis$configuration$n_gene_representatives == 6362L,
  "The R4 real-data selected-unit universe is invalid."
)
matrices <- readRDS(file.path(
  matrix_cache_dir,
  "simulation_correlation_matrices.rds"
))
expect_true(
  identical(sort(names(matrices)), sort(c(
    "direct_centered", "pairwise_difference"
  ))),
  "The R4 simulation matrix cache has unexpected names."
)
for (matrix in matrices) {
  matrix <- validate_time_correlation(matrix, n_time = 16L)
  expect_true(
    min(eigen(matrix, symmetric = TRUE, only.values = TRUE)$values) > 0,
    "An R4 simulation matrix is not positive definite."
  )
}

cat("R3/R4 fashr 0.1.43 static and generative tests passed.\n")
