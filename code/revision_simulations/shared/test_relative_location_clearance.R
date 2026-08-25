#!/usr/bin/env Rscript

# Regression tests for scale-adaptive temporal location truth clearance.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.", call. = FALSE)
}

expect_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

expect_equal <- function(observed, expected, tolerance = 1e-12) {
  expect_true(
    isTRUE(all.equal(observed, expected, tolerance = tolerance)),
    paste("Expected", expected, "but observed", observed)
  )
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

evaluation_grid <- seq(0, 15, by = 0.1)
middle_window <- c(3, 12)

make_location_curve <- function(time_group) {
  curve <- rep(1, length(evaluation_grid))
  target_time <- switch(
    time_group,
    early = 1.5,
    middle = 7.5,
    late = 13.5
  )
  curve[which.min(abs(evaluation_grid - target_time))] <- 3
  curve
}

for (time_group in c("early", "middle", "late")) {
  curve <- make_location_curve(time_group)
  historical <- assess_matched_functional_truth(
    beta_evaluation = curve,
    time_group = time_group,
    switch_status = "non-switch",
    evaluation_grid = evaluation_grid,
    location_truth_margin = 0.10,
    location_truth_min_range_fraction = 0,
    non_switch_min_abs = 0.10,
    non_switch_min_range_fraction = 0.10,
    middle_window = middle_window,
    middle_boundary = "open"
  )
  relative <- assess_matched_functional_truth(
    beta_evaluation = curve,
    time_group = time_group,
    switch_status = "non-switch",
    evaluation_grid = evaluation_grid,
    location_truth_margin = 0.10,
    location_truth_min_range_fraction = 0.75,
    non_switch_min_abs = 0.10,
    non_switch_min_range_fraction = 0.10,
    middle_window = middle_window,
    middle_boundary = "open"
  )

  expect_true(historical$valid, paste(time_group, "historical truth failed."))
  expect_true(relative$valid, paste(time_group, "relative truth failed."))
  expect_equal(historical$required_location_clearance, 0.10)
  expect_equal(relative$required_location_clearance, 1.50)
  expect_equal(relative$effect_range, 2)
  expect_true(
    relative$contrasts[[time_group]] >=
      relative$required_location_clearance,
    paste(time_group, "target contrast missed its relative clearance.")
  )
  competing_groups <- setdiff(c("early", "middle", "late"), time_group)
  expect_true(
    all(relative$contrasts[competing_groups] <=
      -relative$required_location_clearance),
    paste(time_group, "competing contrast missed its relative clearance.")
  )
}

expect_error(
  assess_matched_functional_truth(
    beta_evaluation = make_location_curve("middle"),
    time_group = "middle",
    switch_status = "non-switch",
    evaluation_grid = evaluation_grid,
    location_truth_min_range_fraction = 1,
    middle_window = middle_window,
    middle_boundary = "open"
  ),
  "Invalid matched functional-truth assessment settings."
)

driver_path <- file.path(
  workflowr_root,
  "code/revision_simulations/r3_functional_testing",
  "run_matched_functional_testing_mc_replication.R"
)
driver_text <- readLines(driver_path, warn = FALSE)
expect_true(
  sum(grepl("seed = functional_posterior_seed,", driver_text, fixed = TRUE)) ==
    2L,
  "Raw and BF must both use the named functional posterior seed."
)
expect_true(
  !any(grepl("functional_posterior.*100", driver_text)),
  "The legacy offset functional-posterior seed reappeared."
)

cat("Relative location-clearance regression tests passed.\n")
