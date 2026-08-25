#!/usr/bin/env Rscript

# Deterministic tests for clean working-model null helpers.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
helper_directory <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation"
)
source(file.path(
  helper_directory,
  "selected_signal_genotype_permutation_helpers.R"
))
source(file.path(helper_directory, "small_m_eb_sensitivity_helpers.R"))
source(file.path(helper_directory, "small_m_clean_working_null_helpers.R"))

target_data <- list(
  pair_a = data.frame(
    y = c(10, 20, 30),
    x = c(0, 1, 2),
    offset = c(0, 0, 0)
  ),
  pair_b = data.frame(
    y = c(-2, -1, 0),
    x = c(0, 1, 2),
    offset = c(0.5, 0.5, 0.5)
  )
)
target_se <- list(
  pair_a = c(0.5, 1.0, 2.0),
  pair_b = c(0.25, 0.75, 1.25)
)

first <- simulate_clean_working_null(
  target_data_list = target_data,
  target_se_list = target_se,
  source_unit_keys = names(target_data),
  seed = 101L
)
second <- simulate_clean_working_null(
  target_data_list = target_data,
  target_se_list = target_se,
  source_unit_keys = names(target_data),
  seed = 101L
)
third <- simulate_clean_working_null(
  target_data_list = target_data,
  target_se_list = target_se,
  source_unit_keys = names(target_data),
  seed = 102L
)

stopifnot(
  identical(first, second),
  !isTRUE(all.equal(first$z_matrix, third$z_matrix)),
  identical(dim(first$z_matrix), c(2L, 3L)),
  identical(first$source_unit_keys, names(target_data)),
  identical(first$null_unit_keys, paste0("clean_null::", names(target_data))),
  identical(names(first$data_list), first$null_unit_keys),
  identical(names(first$se_list), first$null_unit_keys),
  identical(first$data_list[[1]]$x, target_data[[1]]$x),
  identical(first$data_list[[2]]$offset, target_data[[2]]$offset),
  identical(unname(first$se_list[[1]]), target_se[[1]]),
  isTRUE(all.equal(
    first$data_list[[1]]$y / first$se_list[[1]],
    unname(first$z_matrix[1, ]),
    tolerance = 1e-15
  )),
  isTRUE(all.equal(
    first$data_list[[2]]$y / first$se_list[[2]],
    unname(first$z_matrix[2, ]),
    tolerance = 1e-15
  ))
)

diagnostics <- summarize_clean_null_z(first$z_matrix)
expected_correlation <- stats::cor(first$z_matrix)
stopifnot(
  nrow(diagnostics$time_summary) == 3L,
  nrow(diagnostics$correlation_long) == 9L,
  isTRUE(all.equal(
    diagnostics$overall$overall_mean,
    mean(first$z_matrix),
    tolerance = 1e-15
  )),
  isTRUE(all.equal(
    diagnostics$overall$overall_sd,
    stats::sd(as.numeric(first$z_matrix)),
    tolerance = 1e-15
  )),
  isTRUE(all.equal(
    diagnostics$overall$mean_off_diagonal_correlation,
    mean(expected_correlation[row(expected_correlation) !=
                                col(expected_correlation)]),
    tolerance = 1e-15
  ))
)

membership <- make_nested_null_subsets(
  n_null = 20L,
  ratio_to_size = c("0.25" = 5L, "0.50" = 10L, "1.00" = 20L),
  n_replicates = 1L,
  seed = 202L
)
indices_025 <- membership$null_index[membership$m_ratio == 0.25]
indices_050 <- membership$null_index[membership$m_ratio == 0.50]
indices_100 <- membership$null_index[membership$m_ratio == 1.00]
stopifnot(
  length(indices_025) == 5L,
  length(indices_050) == 10L,
  length(indices_100) == 20L,
  all(indices_025 %in% indices_050),
  all(indices_050 %in% indices_100),
  !anyDuplicated(indices_100)
)

synthetic_summary <- summarize_small_m_calibration(
  lfdr = c(0.01, 0.04, 0.90, 0.95, 0.02, 0.80),
  group = c(rep("target", 4L), rep("permuted_null", 2L)),
  pi0_merged = 0.8,
  fit_stage = "Synthetic",
  m_ratio = 0.5,
  replicate_id = 1L,
  alpha = 0.05
)
stopifnot(
  synthetic_summary$target_calls == 2L,
  synthetic_summary$permuted_null_calls == 1L,
  synthetic_summary$merged_calls == 3L,
  abs(synthetic_summary$known_null_discovery_fraction - 1 / 3) < 1e-15,
  abs(synthetic_summary$target_pi0_plugin_fdr - 0.7) < 1e-15,
  abs(synthetic_summary$merged_pi0_plugin_fdr - 0.8) < 1e-15
)

cat("Small-M clean working-model null helper tests passed.\n")
