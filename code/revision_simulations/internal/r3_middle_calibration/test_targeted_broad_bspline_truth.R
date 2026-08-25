#!/usr/bin/env Rscript

# Verify the labelled broad-B-spline R3A candidate before a real-genotype
# calibration run. The construction remains non-IWP: it draws smooth cubic
# B-spline lobes, then accepts only curves with a prespecified functional
# separation and switch/non-switch status on the dense evaluation grid.

options(stringsAsFactors = FALSE)

if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
  workflowr_root <- "."
} else if (file.exists(
  "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
)) {
  workflowr_root <- "coderepo-local"
} else {
  stop("Could not locate simulation_functions.R.")
}
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "shared",
  "simulation_functions.R"
))

minimum_location_margin <- 0.60
minimum_location_ratio <- 1.40
middle_window <- c(3, 12)
simulation <- simulate_targeted_local_bspline_effect_set(
  n_variants = 300L,
  time_grid = 0:15,
  evaluation_grid = seq(0, 15, by = 0.1),
  dynamic_amplitude = 2,
  switch_threshold = 0.25,
  minimum_location_margin = minimum_location_margin,
  minimum_location_ratio = minimum_location_ratio,
  non_switch_baseline_fraction = 0.75,
  non_switch_background_fraction = 0.05,
  profile = "broad",
  seed = 98765L,
  middle_window = middle_window,
  middle_boundary = "open",
  scenario = "r3a_targeted_broad_truth_test"
)

stopifnot(identical(
  targeted_time_window(
    "middle", c(3, 3.1, 11, 11.9, 12),
    middle_window = middle_window,
    middle_boundary = "open"
  ),
  c(FALSE, TRUE, TRUE, TRUE, FALSE)
))

recomputed_functionals <- evaluate_temporal_functionals(
  simulation$beta_evaluation,
  smooth_var = simulation$evaluation_grid,
  switch_threshold = 0.25,
  middle_window = middle_window,
  middle_boundary = "open"
)
legacy_functionals <- evaluate_temporal_functionals(
  simulation$beta_evaluation,
  smooth_var = simulation$evaluation_grid,
  switch_threshold = 0.25
)
stopifnot(
  isTRUE(all.equal(simulation$true_functionals, recomputed_functionals)),
  !isTRUE(all.equal(simulation$true_functionals, legacy_functionals)),
  identical(simulation$settings$middle_window, middle_window),
  identical(simulation$settings$middle_boundary, "open")
)

dynamic <- simulation$unit_info[
  simulation$unit_info$effect_class == "dynamic_bspline",
  ,
  drop = FALSE
]
target_column <- match(
  dynamic$time_group,
  colnames(simulation$true_functionals)
)
target_contrast <- simulation$true_functionals[
  cbind(which(simulation$unit_info$effect_class == "dynamic_bspline"), target_column)
]
switch_rows <- dynamic$switch_status == "switch"

if (any(target_contrast < minimum_location_margin) ||
    any(dynamic$target_to_outside_ratio < minimum_location_ratio) ||
    any(simulation$true_functionals[
      which(simulation$unit_info$effect_class == "dynamic_bspline")[switch_rows],
      "switch"
    ] <= 0) ||
    any(simulation$true_functionals[
      which(simulation$unit_info$effect_class == "dynamic_bspline")[!switch_rows],
      "switch"
    ] > 0) ||
    any(dynamic$effective_sign_transitions[switch_rows] != 1L) ||
    any(dynamic$effective_sign_transitions[!switch_rows] != 0L)) {
  stop("The targeted broad B-spline candidate failed its truth contract.")
}

print(aggregate(
  cbind(target_contrast, target_to_outside_ratio, effective_sign_transitions) ~
    truth_group,
  data = transform(
    dynamic,
    target_contrast = target_contrast
  ),
  FUN = min
))
cat(
  "Verified functional separation and switch labels for ",
  nrow(dynamic),
  " broad-B-spline dynamic truths.\n",
  sep = ""
)
