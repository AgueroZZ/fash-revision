#!/usr/bin/env Rscript

# Audit formal R3B truth labels on the 0.1 evaluation grid without fitting FASH.

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
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))

seed_list <- c(12345L, 22345L, 32345L, 42345L, 52345L)
evaluation_grid <- seq(0, 15, by = 0.1)
timing_targets <- c("early", "middle", "late")
location_truth_margin <- 0.10
switch_truth_margin <- 0.10
non_switch_min_abs <- 0.10
non_switch_min_range_fraction <- 0.10
expected_center_ranges <- list(
  early = c(0.5, 1.5),
  middle = c(4.5, 10.5),
  late = c(13.5, 14.5)
)
expected_group_counts <- stats::setNames(
  rep(212L, 6L),
  c(
    "early / switch",
    "early / non-switch",
    "middle / switch",
    "middle / non-switch",
    "late / switch",
    "late / non-switch"
  )
)
tolerance <- 1e-12

audit_rows <- lapply(seed_list, function(seed) {
  effect_sim <- simulate_matched_functional_effect_set(
    n_variants = 6362L,
    truth_mechanism = "raised_cosine",
    evaluation_grid = evaluation_grid,
    class_probs = c(
      dynamic_bspline = 0.20,
      constant = 0.40,
      zero = 0.40
    ),
    dynamic_main_effect_sd = 1,
    switch_threshold = 0.25,
    location_truth_margin = location_truth_margin,
    switch_truth_margin = switch_truth_margin,
    non_switch_min_abs = non_switch_min_abs,
    non_switch_min_range_fraction = non_switch_min_range_fraction,
    seed = seed,
    scenario = "r3b_support_contained_truth_audit",
    middle_window = c(3, 12),
    middle_boundary = "open"
  )

  expect_true(
    identical(
      effect_sim$settings$raised_cosine_center_ranges,
      expected_center_ranges
    ),
    paste("Unexpected center ranges for seed", seed)
  )
  recomputed <- evaluate_temporal_functionals(
    effect_sim$beta_evaluation,
    smooth_var = evaluation_grid,
    switch_threshold = 0.25,
    middle_window = c(3, 12),
    middle_boundary = "open"
  )
  expect_true(
    isTRUE(all.equal(
      effect_sim$true_functionals,
      recomputed,
      tolerance = tolerance
    )),
    paste("Stored and recomputed functionals differ for seed", seed)
  )

  dynamic_index <- which(
    effect_sim$unit_info$effect_class == "dynamic_bspline"
  )
  dynamic_info <- effect_sim$unit_info[dynamic_index, , drop = FALSE]
  observed_group_counts <- table(factor(
    dynamic_info$truth_group,
    levels = names(expected_group_counts)
  ))
  expect_true(
    identical(
      stats::setNames(as.integer(observed_group_counts), names(observed_group_counts)),
      expected_group_counts
    ),
    paste("Unexpected functional-cell counts for seed", seed)
  )

  target_columns <- match(dynamic_info$time_group, colnames(recomputed))
  target_values <- recomputed[cbind(dynamic_index, target_columns)]
  non_target_values <- t(vapply(
    seq_along(dynamic_index),
    function(position) {
      index <- dynamic_index[[position]]
      other_targets <- setdiff(
        timing_targets,
        dynamic_info$time_group[[position]]
      )
      recomputed[index, other_targets]
    },
    numeric(2)
  ))
  target_margin_failures <- sum(
    target_values < location_truth_margin - tolerance
  )
  non_target_margin_failures <- sum(
    apply(
      non_target_values,
      1,
      function(values) any(values > -location_truth_margin + tolerance)
    )
  )

  switch_positions <- which(dynamic_info$switch_status == "switch")
  non_switch_positions <- which(dynamic_info$switch_status == "non-switch")
  switch_margin_failures <- sum(
    recomputed[dynamic_index[switch_positions], "switch"] <
      switch_truth_margin - tolerance
  )
  non_switch_curves <- effect_sim$beta_evaluation[
    dynamic_index[non_switch_positions],
    ,
    drop = FALSE
  ]
  non_switch_minimum_absolute_effect <- apply(
    abs(non_switch_curves),
    1,
    min
  )
  non_switch_effect_range <- apply(non_switch_curves, 1, function(curve) {
    diff(range(curve))
  })
  required_non_switch_clearance <- pmax(
    non_switch_min_abs,
    non_switch_min_range_fraction * non_switch_effect_range
  )
  non_switch_sign_failures <- sum(!apply(
    non_switch_curves,
    1,
    function(curve) all(curve > 0) || all(curve < 0)
  ))
  non_switch_clearance_failures <- sum(
    non_switch_minimum_absolute_effect <
      required_non_switch_clearance - tolerance
  )

  support_containment_failures <- sum(vapply(
    seq_along(dynamic_index),
    function(position) {
      time_group <- dynamic_info$time_group[[position]]
      primary_center <- dynamic_info$peak_centers[[position]][[1L]]
      if (time_group == "early") {
        return(primary_center + 1.5 > 3 + tolerance)
      }
      if (time_group == "middle") {
        return(
          primary_center - 1.5 < 3 - tolerance ||
            primary_center + 1.5 > 12 + tolerance
        )
      }
      primary_center - 1.5 < 12 - tolerance
    },
    logical(1)
  ))

  data.frame(
    seed = seed,
    dynamic_units = length(dynamic_index),
    minimum_target_functional = min(target_values),
    maximum_non_target_functional = max(non_target_values),
    target_margin_failures = target_margin_failures,
    non_target_margin_failures = non_target_margin_failures,
    switch_margin_failures = switch_margin_failures,
    non_switch_sign_failures = non_switch_sign_failures,
    non_switch_clearance_failures = non_switch_clearance_failures,
    support_containment_failures = support_containment_failures,
    maximum_generation_attempt = max(
      dynamic_info$generation_attempt,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
})

audit <- do.call(rbind, audit_rows)
failure_columns <- grep("_failures$", names(audit), value = TRUE)
expect_true(
  all(audit$dynamic_units == 1272L),
  "The truth audit did not generate 1,272 dynamic units per seed."
)
expect_true(
  all(as.matrix(audit[, failure_columns, drop = FALSE]) == 0L),
  "At least one support-contained truth contract failed."
)

output_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "diagnostics",
  "r3_support_contained_truth_audit.csv"
)
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(audit, output_path, row.names = FALSE)

print(audit, row.names = FALSE)
message("R3B support-contained five-seed truth audit passed: ", output_path)
