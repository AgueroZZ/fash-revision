# Regression checks for the center-aligned, relative-clearance, paired R3 run.

simulation_functions_path <- file.path(
  "code", "revision_simulations", "r3_r4_fashr0143", "source_snapshots",
  "r3_center_aligned_relative_location_paired_simulation_functions.R"
)
driver_path <- file.path(
  "code", "revision_simulations", "r3_r4_fashr0143", "source_snapshots",
  "r3_center_aligned_relative_location_paired_driver.R"
)
core_path <- file.path(
  "code", "revision_simulations", "r3_r4_fashr0143",
  "21_r3_center_aligned_relative_location_paired_core.R"
)
wrapper_path <- file.path(
  "code", "revision_simulations", "r3_r4_fashr0143",
  "21_r3_center_aligned_relative_location_paired_fashr0143.R"
)
sbatch_path <- file.path(
  "code", "revision_simulations", "r3_r4_fashr0143",
  "21_r3_center_aligned_relative_location_paired_fashr0143.sbatch"
)

source(simulation_functions_path)

center_ranges <- list(
  early = c(1.5, 2.5),
  middle = c(4.5, 10.5),
  late = c(12.5, 13.5)
)
temporal_category_probs <- stats::setNames(
  rep(1 / 3, 3),
  c("early", "middle", "late")
)
evaluation_grid <- seq(0, 15, by = 0.1)

truth <- simulate_matched_functional_effect_set(
  n_variants = 60L,
  truth_mechanism = "raised_cosine",
  time_grid = 0:15,
  evaluation_grid = evaluation_grid,
  class_probs = c(
    dynamic_bspline = 0.60,
    constant = 0.20,
    zero = 0.20
  ),
  dynamic_main_effect_sd = 1,
  cosine_width_half = 1.5,
  cosine_center_ranges = center_ranges,
  location_truth_margin = 0.10,
  location_truth_min_range_fraction = 0.10,
  switch_truth_margin = 0.10,
  non_switch_min_abs = 0.10,
  non_switch_min_range_fraction = 0.10,
  temporal_category_probs = temporal_category_probs,
  seed = 12345L,
  middle_window = c(3, 12),
  middle_boundary = "open"
)

dynamic <- truth$unit_info$effect_class == "dynamic_bspline"
dynamic_info <- truth$unit_info[dynamic, , drop = FALSE]
expected_truth_groups <- c(
  "early / switch",
  "early / non-switch",
  "middle / switch",
  "middle / non-switch",
  "late / switch",
  "late / non-switch"
)
observed_counts <- table(factor(
  dynamic_info$truth_group,
  levels = expected_truth_groups
))
stopifnot(
  identical(as.integer(observed_counts), rep(6L, 6L)),
  identical(truth$settings$raised_cosine_center_ranges, center_ranges),
  identical(truth$settings$middle_window, c(3, 12)),
  identical(truth$settings$middle_boundary, "open"),
  isTRUE(all.equal(
    evaluation_grid[
      temporal_middle_membership(
        evaluation_grid,
        middle_window = c(3, 12),
        middle_boundary = "open"
      )
    ],
    seq(3.1, 11.9, by = 0.1),
    tolerance = 1e-12,
    check.attributes = FALSE
  ))
)

primary_centers <- vapply(
  truth$unit_info$peak_centers[dynamic],
  function(x) x[[1L]],
  numeric(1)
)
for (time_group in names(center_ranges)) {
  group_rows <- dynamic_info$time_group == time_group
  group_range <- center_ranges[[time_group]]
  stopifnot(
    all(primary_centers[group_rows] >= group_range[[1L]]),
    all(primary_centers[group_rows] <= group_range[[2L]])
  )
}

location_targets <- c("early", "middle", "late")
dynamic_indices <- which(dynamic)
for (index in dynamic_indices) {
  target <- truth$unit_info$time_group[[index]]
  competing_targets <- setdiff(location_targets, target)
  required_margin <- max(
    0.10,
    0.10 * diff(range(truth$beta_evaluation[index, ]))
  )
  stopifnot(
    isTRUE(all.equal(
      truth$unit_info$required_location_clearance[[index]],
      required_margin,
      tolerance = 1e-12
    )),
    truth$true_functionals[index, target] >= required_margin - 1e-12,
    all(
      truth$true_functionals[index, competing_targets] <=
        -required_margin + 1e-12
    )
  )
}

switch_rows <- dynamic_indices[
  truth$unit_info$switch_status[dynamic_indices] == "switch"
]
non_switch_rows <- dynamic_indices[
  truth$unit_info$switch_status[dynamic_indices] == "non-switch"
]
stopifnot(
  all(truth$true_functionals[switch_rows, "switch"] >= 0.10 - 1e-12),
  all(vapply(
    non_switch_rows,
    function(index) {
      curve <- truth$beta_evaluation[index, ]
      same_sign <- all(curve > 0) || all(curve < 0)
      required_clearance <- max(0.10, 0.10 * diff(range(curve)))
      same_sign && min(abs(curve)) >= required_clearance - 1e-12
    },
    logical(1)
  ))
)

driver_text <- paste(readLines(driver_path, warn = FALSE), collapse = "\n")
core_text <- paste(readLines(core_path, warn = FALSE), collapse = "\n")
wrapper_text <- paste(readLines(wrapper_path, warn = FALSE), collapse = "\n")
sbatch_text <- paste(readLines(sbatch_path, warn = FALSE), collapse = "\n")
expected_result_id <- paste0(
  "r3_real_genotype_one_per_gene_J6362_",
  "matched_functional_open_middle_3_12_center_aligned_equal_cells_",
  "relative_location_clearance_paired_posterior_fashr0143_pilot5"
)

seed_argument_count <- lengths(regmatches(
  driver_text,
  gregexpr("seed = functional_posterior_seed", driver_text, fixed = TRUE)
))
stopifnot(
  identical(seed_argument_count, 2L),
  !grepl("functional_posterior.*\\+ *100", driver_text),
  grepl(
    "cosine_center_ranges = configuration$raised_cosine$center_ranges",
    driver_text,
    fixed = TRUE
  ),
  grepl("equal temporal categories", driver_text, fixed = TRUE),
  grepl("center_aligned_equal_cells_", core_text, fixed = TRUE),
  grepl("center_aligned_equal_cells_", wrapper_text, fixed = TRUE),
  grepl(expected_result_id, sbatch_text, fixed = TRUE),
  grepl("common_random_seed_raw_bf", core_text, fixed = TRUE),
  grepl("temporal_category_probs_match", core_text, fixed = TRUE),
  grepl("tolerance = 1e-12", core_text, fixed = TRUE),
  grepl('sprintf("%.17g", TEMPORAL_CATEGORY_PROBS)', core_text, fixed = TRUE)
)

message(
  "Center-aligned relative-clearance paired-posterior R3 tests passed."
)
