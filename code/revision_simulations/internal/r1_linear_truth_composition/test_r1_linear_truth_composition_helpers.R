#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr project root.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "r1_linear_truth_composition", "r1_linear_truth_composition_helpers.R"
))

time_grid <- 0:15
positive <- make_centered_linear_deviation(time_grid, amplitude = 2, direction = 1)
negative <- make_centered_linear_deviation(time_grid, amplitude = 2, direction = -1)
stopifnot(
  abs(mean(positive)) < 1e-12,
  abs(max(abs(positive)) - 2) < 1e-12,
  identical(negative, -positive)
)

base_effect_sim <- simulate_variant_effect_curves(
  n_variants = 100,
  time_grid = time_grid,
  class_probs = c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40),
  scenario = "toy_r1_random_bspline",
  dynamic_amplitude = 2,
  bspline_df = 6,
  bspline_coefficient_sd = 1,
  constant_sd = 1,
  dynamic_main_effect_sd = 1,
  exact_class_counts = TRUE,
  seed = 101
)

truths <- make_r1_linear_truth_scenarios(
  base_effect_sim = base_effect_sim,
  time_grid = time_grid,
  linear_amplitude = 2,
  linear_sign_seed = 202,
  mixture_seed = 303
)
truths_repeat <- make_r1_linear_truth_scenarios(
  base_effect_sim = base_effect_sim,
  time_grid = time_grid,
  linear_amplitude = 2,
  linear_sign_seed = 202,
  mixture_seed = 303
)

all_linear <- truths$scenarios$all_linear
mixed <- truths$scenarios$linear90_bspline10
base_classes <- base_effect_sim$unit_info$effect_class
dynamic_indices <- which(base_classes == "dynamic_bspline")
non_dynamic_indices <- which(base_classes != "dynamic_bspline")
mixed_bspline_indices <- which(mixed$unit_info$effect_class == "dynamic_bspline")
mixed_linear_indices <- which(mixed$unit_info$effect_class == "dynamic_linear")

stopifnot(
  sum(all_linear$unit_info$effect_class == "dynamic_linear") == 20L,
  sum(all_linear$unit_info$effect_class == "dynamic_bspline") == 0L,
  sum(mixed$unit_info$effect_class == "dynamic_linear") == 18L,
  sum(mixed$unit_info$effect_class == "dynamic_bspline") == 2L,
  identical(truths$mixture_counts, c(dynamic_linear = 18L, dynamic_bspline = 2L)),
  identical(
    all_linear$beta_matrix[non_dynamic_indices, , drop = FALSE],
    base_effect_sim$beta_matrix[non_dynamic_indices, , drop = FALSE]
  ),
  identical(
    mixed$beta_matrix[non_dynamic_indices, , drop = FALSE],
    base_effect_sim$beta_matrix[non_dynamic_indices, , drop = FALSE]
  ),
  identical(
    mixed$beta_matrix[mixed_bspline_indices, , drop = FALSE],
    base_effect_sim$beta_matrix[mixed_bspline_indices, , drop = FALSE]
  ),
  identical(
    mixed$beta_matrix[mixed_linear_indices, , drop = FALSE],
    all_linear$beta_matrix[mixed_linear_indices, , drop = FALSE]
  ),
  identical(truths$membership, truths_repeat$membership),
  identical(
    truths$scenarios$linear90_bspline10$beta_matrix,
    truths_repeat$scenarios$linear90_bspline10$beta_matrix
  )
)

for (index in dynamic_indices) {
  deviation <-
    all_linear$beta_matrix[index, ] -
    all_linear$unit_info$genetic_main_effect[index]
  stopifnot(
    abs(mean(deviation)) < 1e-12,
    abs(max(abs(deviation)) - 2) < 1e-12
  )
}

stopifnot(
  nrow(truths$membership) == 200L,
  !any(truths$membership$retained_bspline[
    truths$membership$scenario == "r1_all_linear_dynamic_truth"
  ]),
  all(
    truths$membership$selected_for_mixed_bspline[
      truths$membership$scenario == "r1_all_linear_dynamic_truth"
    ] ==
      truths$membership$unit_index[
        truths$membership$scenario == "r1_all_linear_dynamic_truth"
      ] %in% mixed_bspline_indices
  )
)

cat("R1 linear-truth composition helper tests passed.\n")
