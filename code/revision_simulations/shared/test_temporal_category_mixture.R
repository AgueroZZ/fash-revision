#!/usr/bin/env Rscript

# Regression tests for configurable R3 temporal truth-category mixtures.

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

expect_error <- function(expression, pattern) {
  observed <- tryCatch(
    {
      force(expression)
      NULL
    },
    error = function(error) conditionMessage(error)
  )
  expect_true(
    !is.null(observed) && grepl(pattern, observed, fixed = TRUE),
    paste("Expected error containing:", pattern)
  )
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code/revision_simulations/shared/simulation_functions.R"
))

truth_group_levels <- c(
  "early / switch",
  "early / non-switch",
  "middle / switch",
  "middle / non-switch",
  "late / switch",
  "late / non-switch"
)

balanced <- exact_temporal_truth_group_counts(1272L)
expect_true(
  identical(
    as.integer(balanced[truth_group_levels]),
    rep(212L, length(truth_group_levels))
  ),
  "The default balanced truth counts changed."
)

reference_probs <- c(early = 0.29, middle = 0.42, late = 0.29)
reference <- exact_temporal_truth_group_counts(
  1272L,
  temporal_category_probs = reference_probs
)
expect_true(
  identical(
    as.integer(reference[truth_group_levels]),
    c(185L, 184L, 267L, 267L, 185L, 184L)
  ),
  "The canonical IWP1 truth-group counts are incorrect."
)
expect_true(sum(reference) == 1272L, "Truth-group counts do not sum to 1272.")

expect_error(
  exact_temporal_truth_group_counts(
    100L,
    c(early = 0.30, middle = 0.40, late = 0.20)
  ),
  "sum to one"
)
expect_error(
  exact_temporal_truth_group_counts(
    100L,
    c(middle = 0.42, early = 0.29, late = 0.29)
  ),
  "early, middle, and late"
)

small <- simulate_matched_functional_effect_set(
  n_variants = 60L,
  truth_mechanism = "random_bspline",
  time_grid = make_time_grid(),
  evaluation_grid = seq(0, 15, by = 0.1),
  class_probs = c(
    dynamic_bspline = 0.20,
    constant = 0.40,
    zero = 0.40
  ),
  temporal_category_probs = reference_probs,
  middle_window = c(3, 12),
  middle_boundary = "open",
  seed = 9917L
)
dynamic <- small$unit_info$effect_class == "dynamic_bspline"
observed_counts <- table(factor(
  small$unit_info$truth_group[dynamic],
  levels = truth_group_levels
))
expected_counts <- exact_temporal_truth_group_counts(
  sum(dynamic),
  temporal_category_probs = reference_probs
)[truth_group_levels]
expect_true(
  identical(as.integer(observed_counts), as.integer(expected_counts)),
  "The generated truth groups do not match the requested mixture."
)
target_columns <- match(
  small$unit_info$time_group[dynamic],
  colnames(small$true_functionals)
)
expect_true(
  all(small$true_functionals[cbind(which(dynamic), target_columns)] >= 0.10),
  "At least one generated target functional violates the truth margin."
)
expect_true(
  identical(small$settings$temporal_category_probs, reference_probs),
  "The requested temporal mixture is missing from simulation provenance."
)

cat("Temporal truth-category mixture regression tests passed.\n")
