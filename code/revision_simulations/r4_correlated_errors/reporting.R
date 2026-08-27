# Load, validate, and plot the R4 unit-specific residual-permutation cache.

r4_output_dir <- file.path(
  "output", "revision_simulations", "real_data",
  "r4_unit_specific_residual_permutation_z"
)
r4_result_path <- file.path(r4_output_dir, "unit_specific_null_correlation.rds")
r4_complete_flag_path <- file.path(r4_output_dir, "complete.flag")
if (!file.exists(r4_result_path) || !file.exists(r4_complete_flag_path)) {
  stop("The completed R4 unit-specific residual-permutation cache is missing.")
}

r4_analysis <- readRDS(r4_result_path)
r4_configuration <- r4_analysis$configuration
required_configuration <- list(
  schema_version = "r4-unit-specific-residual-permutation-z-v1",
  output_id = "r4_unit_specific_residual_permutation_z",
  n_units = 6L,
  n_permutations = 400L,
  selection_seed = 20260824L,
  permutation_seed = 20260825L,
  time_grid = 0:15
)
for (field in names(required_configuration)) {
  if (!identical(r4_configuration[[field]], required_configuration[[field]])) {
    stop("The R4 cache has an unexpected ", field, ".")
  }
}
completion <- readLines(r4_complete_flag_path, warn = FALSE)
if (!all(c(
  "output_id=r4_unit_specific_residual_permutation_z",
  "n_units=6",
  "n_permutations=400",
  "residual_method=signal_stripped_full_donor_trajectory_permutation"
) %in% completion)) {
  stop("The R4 unit-specific complete.flag is invalid.")
}
if (nrow(r4_analysis$selected_units) != 6L ||
    length(r4_analysis$unit_results) != 6L ||
    anyDuplicated(r4_analysis$selected_units$pair_key) ||
    nrow(r4_analysis$donor_maps) != 6L * 400L * 19L) {
  stop("The R4 unit-specific cache is incomplete.")
}
for (unit_result in r4_analysis$unit_results) {
  if (!identical(dim(unit_result$null_beta_draws), c(400L, 16L)) ||
      !identical(dim(unit_result$null_t_adjusted_se_draws), c(400L, 16L)) ||
      !identical(dim(unit_result$null_z_draws), c(400L, 16L)) ||
      !identical(dim(unit_result$correlation), c(16L, 16L)) ||
      nrow(unit_result$variogram) != 15L ||
      any(!is.finite(unit_result$null_beta_draws)) ||
      any(!is.finite(unit_result$null_t_adjusted_se_draws)) ||
      any(unit_result$null_t_adjusted_se_draws <= 0) ||
      any(!is.finite(unit_result$null_z_draws)) ||
      any(!is.finite(unit_result$correlation)) ||
      unit_result$maximum_residual_genotype_correlation > 1e-10) {
    stop("An R4 unit-specific null result failed validation.")
  }
}

r4_residual_output_dir <- file.path(
  "output", "revision_simulations", "real_data",
  "r4_full_model_residual_expression_correlation"
)
r4_residual_result_path <- file.path(
  r4_residual_output_dir, "full_model_residual_correlation.rds"
)
r4_residual_complete_flag_path <- file.path(r4_residual_output_dir, "complete.flag")
if (!file.exists(r4_residual_result_path) ||
    !file.exists(r4_residual_complete_flag_path)) {
  stop("The R4 full-model residual correlation cache is missing.")
}
r4_residual_analysis <- readRDS(r4_residual_result_path)
r4_residual_configuration <- r4_residual_analysis$configuration
if (!identical(
  r4_residual_configuration$schema_version,
  "r4-full-model-residual-correlation-v1"
) || !identical(r4_residual_configuration$n_units, 6L) ||
    !identical(r4_residual_configuration$time_grid, 0:15) ||
    !identical(
      r4_residual_configuration$correlation_definition,
      "stats::cor(..., use = 'pairwise.complete.obs')"
    ) || !identical(
      r4_residual_configuration$model,
      "Y ~ 1 + PC1 + ... + PC5 + G"
    ) || !identical(r4_residual_configuration$genotype_in_residualization, TRUE) ||
    !identical(
      r4_residual_analysis$selected_units$pair_key,
      r4_analysis$selected_units$pair_key
    )) {
  stop("The R4 full-model residual cache has an unexpected design.")
}
residual_completion <- readLines(r4_residual_complete_flag_path, warn = FALSE)
if (!all(c(
  "result_id=r4_full_model_residual_expression_correlation",
  "n_units=6",
  "model=Y~1+PC1+...+PC5+G",
  "genotype_in_residualization=true",
  "correlation=pairwise_complete_donors",
  "minimum_pairwise_donors=13"
) %in% residual_completion)) {
  stop("The R4 full-model residual complete.flag is invalid.")
}
if (length(r4_residual_analysis$unit_results) != 6L ||
    any(vapply(r4_residual_analysis$unit_results, function(unit_result) {
      !identical(dim(unit_result$correlation), c(16L, 16L)) ||
      nrow(unit_result$variogram) != 15L ||
      min(unit_result$pairwise_donor_count) < 13L ||
        unit_result$maximum_observed_beta_difference > 1e-10 ||
        any(!is.finite(unit_result$correlation))
    }, logical(1)))) {
  stop("An R4 full-model residual result failed validation.")
}

r4_selected_units <- r4_analysis$selected_units
r4_selected_units$unit <- paste("Unit", r4_selected_units$selected_order)
r4_selected_units$source_bf_lfdr <- round(r4_selected_units$source_bf_lfdr, 4)
r4_selected_units$maximum_observed_beta_difference <- format(
  r4_selected_units$maximum_observed_beta_difference,
  scientific = TRUE,
  digits = 2
)
r4_unit_label <- function(unit_result) {
  paste0(
    "Unit ", unit_result$selected_unit$selected_order,
    ": ", unit_result$selected_unit$pair_key
  )
}

r4_short_unit_label <- function(unit_result) {
  paste("Unit", unit_result$selected_unit$selected_order)
}

r4_variogram_bootstrap_reps <- 1000L
r4_variogram_intervals <- do.call(rbind, lapply(
  seq_along(r4_analysis$unit_results), function(unit_index) {
    unit_result <- r4_analysis$unit_results[[unit_index]]
    transform(
      bootstrap_r4_null_variogram_intervals(
        beta_draws = unit_result$null_z_draws,
        n_bootstrap = r4_variogram_bootstrap_reps,
        seed = 20260826L + unit_index
      ),
      unit_label = r4_short_unit_label(unit_result)
    )
  }
))

plot_r4_unit_null_heatmaps <- function(unit_results) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to plot R4 null correlations.")
  }
  heatmap_data <- do.call(rbind, lapply(unit_results, function(unit_result) {
    correlation_matrix_to_long_r4(
      unit_result$correlation,
      unit_label = r4_unit_label(unit_result)
    )
  }))
  heatmap_data$unit_label <- factor(
    heatmap_data$unit_label,
    levels = unique(heatmap_data$unit_label)
  )
  ggplot2::ggplot(
    heatmap_data,
    ggplot2::aes(x = time_column, y = time_row, fill = correlation)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.12) +
    ggplot2::coord_equal() +
    ggplot2::scale_y_reverse(breaks = seq(0, 15, by = 3)) +
    ggplot2::scale_x_continuous(breaks = seq(0, 15, by = 3)) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 0, limits = c(-1, 1), name = "Correlation"
    ) +
    ggplot2::facet_wrap(~unit_label, nrow = 1) +
    ggplot2::labs(x = "Time", y = "Time") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 8),
      legend.position = "bottom"
    )
}

plot_r4_unit_null_variograms <- function(unit_results) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to plot R4 variograms.")
  }
  selected_labels <- vapply(unit_results, r4_short_unit_label, character(1))
  variogram_data <- r4_variogram_intervals[
    r4_variogram_intervals$unit_label %in% selected_labels,
    , drop = FALSE
  ]
  variogram_data$unit_label <- factor(
    variogram_data$unit_label,
    levels = unique(variogram_data$unit_label)
  )
  ggplot2::ggplot(
    variogram_data,
    ggplot2::aes(x = lag, y = semivariogram, color = unit_label)
  ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey55") +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper, fill = unit_label),
      alpha = 0.14, linewidth = 0, show.legend = FALSE
    ) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_x_continuous(breaks = seq(1, 15, by = 2)) +
    ggplot2::labs(
      x = "Time lag", y = "Estimated semivariogram",
      color = "Randomly selected unit"
    ) +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

plot_r4_full_model_residual_variograms <- function(residual_results) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required to plot R4 full-model residual variograms.")
  }
  variogram_data <- do.call(rbind, lapply(residual_results, function(unit_result) {
    transform(
      unit_result$variogram,
    unit_label = r4_short_unit_label(unit_result)
    )
  }))
  variogram_data$unit_label <- factor(
    variogram_data$unit_label,
    levels = paste("Unit", seq_along(residual_results))
  )
  ggplot2::ggplot(
    variogram_data,
    ggplot2::aes(x = lag, y = semivariogram, color = unit_label)
  ) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey55") +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_x_continuous(breaks = seq(1, 15, by = 2)) +
    ggplot2::labs(
      x = "Time lag",
      y = "Full-model residual semivariogram",
      color = "Randomly selected unit"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

r4_variogram_summary <- do.call(rbind, lapply(r4_analysis$unit_results, function(unit_result) {
  first_lag <- unit_result$variogram[unit_result$variogram$lag == 1L, , drop = FALSE]
  last_lags <- unit_result$variogram[unit_result$variogram$lag >= 9L, , drop = FALSE]
  data.frame(
    unit = paste("Unit", unit_result$selected_unit$selected_order),
    pair_key = unit_result$selected_unit$pair_key,
    lag_1_correlation = first_lag$mean_correlation,
    mean_lag_9_to_15_correlation = mean(last_lags$mean_correlation),
    stringsAsFactors = FALSE
  )
}))
r4_variogram_summary[, c("lag_1_correlation", "mean_lag_9_to_15_correlation")] <-
  lapply(r4_variogram_summary[, c("lag_1_correlation", "mean_lag_9_to_15_correlation")], round, digits = 3)
