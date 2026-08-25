#!/usr/bin/env Rscript

# Regression test for the formal R3 full-support raised-cosine truth geometry.

find_workflowr_root <- function() {
  if (file.exists(file.path(
    "code", "revision_simulations", "r3_r4_fashr0143",
    "source_snapshots", "r3_prior_geometry_simulation_functions.R"
  ))) return(".")
  if (file.exists(file.path(
    "coderepo-local", "code", "revision_simulations", "r3_r4_fashr0143",
    "source_snapshots", "r3_prior_geometry_simulation_functions.R"
  ))) return("coderepo-local")
  stop("Could not find the workflowr repository root.", call. = FALSE)
}

expect_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root, "code", "revision_simulations", "r3_r4_fashr0143",
  "source_snapshots", "r3_prior_geometry_simulation_functions.R"
))

temporal_probs <- c(early = 0.29, middle = 0.42, late = 0.29)
simulation <- simulate_matched_functional_effect_set(
  n_variants = 120L,
  truth_mechanism = "raised_cosine",
  temporal_category_probs = temporal_probs,
  middle_window = c(3, 12),
  middle_boundary = "open",
  seed = 812345L
)

dynamic <- simulation$unit_info$effect_class == "dynamic_bspline"
unit_info <- simulation$unit_info[dynamic, , drop = FALSE]
primary_centers <- vapply(
  unit_info$peak_centers,
  function(x) x[[1L]],
  numeric(1)
)
support_left <- primary_centers - simulation$settings$cosine_width_half
support_right <- primary_centers + simulation$settings$cosine_width_half

expect_true(
  identical(
    simulation$settings$raised_cosine_center_ranges,
    list(
      early = c(1.5, 1.5),
      middle = c(4.5, 10.5),
      late = c(13.5, 13.5)
    )
  ),
  "Unexpected frozen raised-cosine center ranges."
)
expect_true(
  all(primary_centers[unit_info$time_group == "early"] == 1.5) &&
    all(primary_centers[unit_info$time_group == "late"] == 13.5),
  "Early or Late primary centers are not fixed at the full-support locations."
)
expect_true(
  all(support_left >= 0) && all(support_right <= 15),
  "At least one primary peak extends beyond the observed time domain."
)
expect_true(
  all(support_left[unit_info$time_group == "early"] >= 0) &&
    all(support_right[unit_info$time_group == "early"] <= 3) &&
    all(support_left[unit_info$time_group == "middle"] >= 3) &&
    all(support_right[unit_info$time_group == "middle"] <= 12) &&
    all(support_left[unit_info$time_group == "late"] >= 12) &&
    all(support_right[unit_info$time_group == "late"] <= 15),
  "At least one primary peak extends outside its target time region."
)

target_columns <- match(
  unit_info$time_group,
  colnames(simulation$true_functionals)
)
target_values <- simulation$true_functionals[
  cbind(which(dynamic), target_columns)
]
expect_true(
  all(target_values >= simulation$settings$location_truth_margin),
  "At least one generated curve fails the target-functional truth margin."
)

expected_counts <- exact_temporal_truth_group_counts(
  sum(dynamic), temporal_category_probs = temporal_probs
)
observed_counts <- table(factor(
  unit_info$truth_group, levels = names(expected_counts)
))
expect_true(
  identical(as.integer(observed_counts), as.integer(expected_counts)),
  "The generated truth groups do not match the frozen temporal mixture."
)

cat("R3 full-support raised-cosine geometry regression test passed.\n")
cat("Dynamic truth-group counts:\n")
print(expected_counts)
