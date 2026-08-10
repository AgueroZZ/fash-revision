# Reporting helpers for the matched donor-permutation variogram analysis.

permuted_workflowr_root <- if (
  file.exists("code/revision_simulations/shared/simulation_functions.R")
) {
  "."
} else if (file.exists(
  "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
)) {
  "coderepo-local"
} else {
  stop("Could not find the workflowr repository root.")
}

permuted_cache_dir <- file.path(
  permuted_workflowr_root,
  "output", "revision_simulations", "internal",
  "permuted_variogram_covariance_summary"
)
permuted_analysis_path <- file.path(
  permuted_cache_dir,
  "permuted_variogram_covariance_summary.rds"
)
if (!file.exists(permuted_analysis_path)) {
  stop("The permuted variogram/covariance cache is missing.")
}
permuted_analysis <- readRDS(permuted_analysis_path)

permuted_required_names <- c(
  "configuration",
  "aggregate_lag_summaries",
  "aggregate_method_contrasts",
  "mean_covariance_matrices_long",
  "mean_matrix_diagnostics",
  "matrix_agreement",
  "independent_covariance_draws",
  "covariance_signal_draws"
)
if (!all(permuted_required_names %in% names(permuted_analysis)) ||
    !identical(
      permuted_analysis$configuration$analysis_id,
      "permuted_variogram_covariance_summary"
    ) || permuted_analysis$configuration$n_replications != 100L ||
    permuted_analysis$configuration$primary_threshold != 0.95 ||
    !identical(
      permuted_analysis$configuration$primary_se_scale,
      "t_adjusted"
    )) {
  stop("The permuted variogram/covariance cache is incompatible.")
}

permuted_configuration <- permuted_analysis$configuration
permuted_lag_summaries <- permuted_analysis$aggregate_lag_summaries
permuted_matrix_diagnostics <- permuted_analysis$mean_matrix_diagnostics
permuted_matrix_agreement <- permuted_analysis$matrix_agreement
permuted_method_labels <- c(
  pairwise_ols = "Pairwise OLS",
  within_unit_centered = "Within-unit centered",
  direct_zero_null = "Direct zero-null"
)
permuted_method_order <- names(permuted_method_labels)
permuted_primary_lags <- c(1L, 9L, 15L)

permuted_primary_lag_data <- permuted_lag_summaries[
  permuted_lag_summaries$threshold == 0.95 &
    permuted_lag_summaries$se_scale == "t_adjusted" &
    permuted_lag_summaries$metric == "covariance" &
    permuted_lag_summaries$lag %in% permuted_primary_lags,
]
permuted_primary_lag_table <- do.call(rbind, lapply(
  permuted_method_order,
  function(estimator) {
    estimator_data <- permuted_primary_lag_data[
      permuted_primary_lag_data$estimator == estimator,
    ]
    estimator_data <- estimator_data[
      match(permuted_primary_lags, estimator_data$lag),
    ]
    data.frame(
      Estimator = unname(permuted_method_labels[estimator]),
      `Lag 1 covariance (95% interval)` = format_interval(
        estimator_data$mean[1L],
        estimator_data$lower[1L],
        estimator_data$upper[1L]
      ),
      `Lag 9 covariance (95% interval)` = format_interval(
        estimator_data$mean[2L],
        estimator_data$lower[2L],
        estimator_data$upper[2L]
      ),
      `Lag 15 covariance (95% interval)` = format_interval(
        estimator_data$mean[3L],
        estimator_data$lower[3L],
        estimator_data$upper[3L]
      ),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
))

permuted_primary_agreement <- permuted_matrix_agreement[
  permuted_matrix_agreement$threshold == 0.95 &
    permuted_matrix_agreement$se_scale == "t_adjusted",
]
permuted_primary_agreement <- permuted_primary_agreement[match(
  c("pairwise_ols", "within_unit_centered"),
  permuted_primary_agreement$estimator
), ]
permuted_method_agreement_table <- data.frame(
  Estimator = unname(permuted_method_labels[
    permuted_primary_agreement$estimator
  ]),
  `Off-diagonal matrix correlation with direct` =
    permuted_primary_agreement$off_diagonal_matrix_correlation,
  `Off-diagonal RMSE` = permuted_primary_agreement$off_diagonal_rmse,
  `Mean covariance difference` =
    permuted_primary_agreement$off_diagonal_mean_difference,
  `Maximum absolute difference` =
    permuted_primary_agreement$off_diagonal_maximum_absolute_difference,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

permuted_primary_diagnostics <- permuted_matrix_diagnostics[
  permuted_matrix_diagnostics$threshold == 0.95 &
    permuted_matrix_diagnostics$se_scale == "t_adjusted",
]
permuted_primary_diagnostics <- permuted_primary_diagnostics[match(
  permuted_method_order,
  permuted_primary_diagnostics$estimator
), ]
permuted_matrix_diagnostics_table <- data.frame(
  Estimator = unname(permuted_method_labels[
    permuted_primary_diagnostics$estimator
  ]),
  `Mean diagonal` = permuted_primary_diagnostics$mean_diagonal,
  `Mean off-diagonal` = permuted_primary_diagnostics$mean_off_diagonal,
  `Off-diagonal minimum` =
    permuted_primary_diagnostics$minimum_off_diagonal,
  `Off-diagonal maximum` =
    permuted_primary_diagnostics$maximum_off_diagonal,
  `Minimum eigenvalue` = permuted_primary_diagnostics$minimum_eigenvalue,
  `Leading eigenvalue fraction` =
    permuted_primary_diagnostics$leading_eigenvalue_fraction,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

permuted_threshold_data <- permuted_lag_summaries[
  permuted_lag_summaries$se_scale == "t_adjusted" &
    permuted_lag_summaries$estimator == "direct_zero_null" &
    permuted_lag_summaries$metric == "covariance" &
    permuted_lag_summaries$lag %in% permuted_primary_lags,
]
permuted_threshold_sensitivity_table <- do.call(rbind, lapply(
  sort(unique(permuted_threshold_data$threshold), decreasing = TRUE),
  function(threshold) {
    threshold_data <- permuted_threshold_data[
      permuted_threshold_data$threshold == threshold,
    ]
    threshold_data <- threshold_data[
      match(permuted_primary_lags, threshold_data$lag),
    ]
    data.frame(
      `Strict lfdr rule` = paste0(
        "lfdr > ",
        format_decimal(threshold, 2L)
      ),
      Units = threshold_data$n_units[1L],
      `Lag 1 covariance` = threshold_data$mean[1L],
      `Lag 9 covariance` = threshold_data$mean[2L],
      `Lag 15 covariance` = threshold_data$mean[3L],
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
))

permuted_primary_matrix_data <-
  permuted_analysis$mean_covariance_matrices_long[
    permuted_analysis$mean_covariance_matrices_long$threshold == 0.95 &
      permuted_analysis$mean_covariance_matrices_long$se_scale ==
        "t_adjusted" &
      permuted_analysis$mean_covariance_matrices_long$estimator ==
        "direct_zero_null",
  ]
permuted_time_grid <- sort(unique(permuted_primary_matrix_data$time_a))
permuted_primary_matrix <- matrix(
  permuted_primary_matrix_data$covariance,
  nrow = length(permuted_time_grid),
  ncol = length(permuted_time_grid),
  dimnames = list(
    paste0("Time ", permuted_time_grid),
    paste0("Time ", permuted_time_grid)
  )
)
permuted_primary_matrix_table <- data.frame(
  Time = permuted_time_grid,
  round(permuted_primary_matrix, 3L),
  check.names = FALSE
)
colnames(permuted_primary_matrix_table)[-1L] <- paste0(
  "Time ",
  permuted_time_grid
)

permuted_direct_diagnostics <- permuted_primary_diagnostics[
  permuted_primary_diagnostics$estimator == "direct_zero_null",
]
permuted_primary_key <- "lfdr_0p95__t_adjusted__direct_zero_null"
permuted_independent_mean <- apply(
  permuted_analysis$independent_covariance_draws[[permuted_primary_key]],
  c(1L, 2L),
  mean
)
permuted_independent_off_diagonal_mean <- mean(
  permuted_independent_mean[upper.tri(permuted_independent_mean)]
)
permuted_independent_off_diagonal_range <- range(
  permuted_independent_mean[upper.tri(permuted_independent_mean)]
)

permuted_figure_paths <- c(
  variogram_covariance = file.path(
    permuted_cache_dir,
    "figure",
    "matched_variogram_and_covariance_t_adjusted.png"
  ),
  estimator_matrices = file.path(
    permuted_cache_dir,
    "figure",
    "matched_empirical_covariance_matrices_t_adjusted.png"
  ),
  direct_matrix = file.path(
    permuted_cache_dir,
    "figure",
    "direct_empirical_covariance_matrix_t_adjusted.png"
  )
)
if (any(!file.exists(permuted_figure_paths))) {
  stop("At least one permuted variogram/covariance figure is missing.")
}

permuted_colors <- c(
  pairwise_ols = "#0072B2",
  within_unit_centered = "#D55E00",
  direct_zero_null = "#009E73"
)

plot_permuted_variogram_covariance <- function() {
  plot_data <- permuted_lag_summaries[
    permuted_lag_summaries$threshold == 0.95 &
      permuted_lag_summaries$se_scale == "t_adjusted",
  ]
  synchronized_data <- plot_data[
    plot_data$metric == "synchronized_variogram",
  ]
  independent_data <- plot_data[
    plot_data$metric == "independent_variogram",
  ]
  covariance_data <- plot_data[plot_data$metric == "covariance", ]
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::layout(matrix(1:2, ncol = 1L))

  graphics::par(mar = c(3.0, 5.1, 3.4, 1.2))
  top_range <- range(c(
    synchronized_data$lower,
    synchronized_data$upper,
    independent_data$lower,
    independent_data$upper
  ))
  top_range[2L] <- max(top_range[2L], 1.22)
  graphics::plot(
    NA,
    xlim = c(1, 15),
    ylim = top_range,
    xlab = "",
    ylab = "Covariance-scale semivariogram",
    main = "Synchronized variogram and matched independence benchmark",
    xaxt = "n"
  )
  graphics::axis(1, at = 1:15, labels = FALSE)
  for (estimator in permuted_method_order) {
    synchronized_estimator <- synchronized_data[
      synchronized_data$estimator == estimator,
    ]
    synchronized_estimator <- synchronized_estimator[
      order(synchronized_estimator$lag),
    ]
    independent_estimator <- independent_data[
      independent_data$estimator == estimator,
    ]
    independent_estimator <- independent_estimator[
      order(independent_estimator$lag),
    ]
    graphics::polygon(
      c(
        synchronized_estimator$lag,
        rev(synchronized_estimator$lag)
      ),
      c(
        synchronized_estimator$lower,
        rev(synchronized_estimator$upper)
      ),
      col = grDevices::adjustcolor(
        permuted_colors[estimator],
        alpha.f = 0.12
      ),
      border = NA
    )
    graphics::lines(
      synchronized_estimator$lag,
      synchronized_estimator$mean,
      col = permuted_colors[estimator],
      lwd = 2.5
    )
    graphics::lines(
      independent_estimator$lag,
      independent_estimator$mean,
      col = permuted_colors[estimator],
      lwd = 2.1,
      lty = 2
    )
  }
  graphics::legend(
    "topleft",
    legend = c(
      permuted_method_labels[permuted_method_order],
      "Solid: synchronized",
      "Dashed: independent-time benchmark"
    ),
    col = c(
      permuted_colors[permuted_method_order],
      "black",
      "black"
    ),
    lwd = c(rep(2.5, 3L), 2.2, 2.2),
    lty = c(rep(1, 4L), 2),
    bty = "n",
    cex = 0.82
  )

  graphics::par(mar = c(4.8, 5.1, 3.4, 1.2))
  bottom_range <- range(c(covariance_data$lower, covariance_data$upper, 0))
  graphics::plot(
    NA,
    xlim = c(1, 15),
    ylim = bottom_range,
    xlab = "Lag",
    ylab = "Estimated covariance",
    main = "Independence benchmark minus synchronized variogram",
    xaxt = "n"
  )
  graphics::axis(1, at = 1:15)
  graphics::abline(h = 0, lty = 3, col = "grey45")
  for (estimator in permuted_method_order) {
    estimator_data <- covariance_data[covariance_data$estimator == estimator, ]
    estimator_data <- estimator_data[order(estimator_data$lag), ]
    graphics::polygon(
      c(estimator_data$lag, rev(estimator_data$lag)),
      c(estimator_data$lower, rev(estimator_data$upper)),
      col = grDevices::adjustcolor(
        permuted_colors[estimator],
        alpha.f = 0.14
      ),
      border = NA
    )
    graphics::lines(
      estimator_data$lag,
      estimator_data$mean,
      col = permuted_colors[estimator],
      lwd = 2.5
    )
  }
  graphics::legend(
    "topright",
    legend = permuted_method_labels[permuted_method_order],
    col = permuted_colors[permuted_method_order],
    lwd = 2.5,
    bty = "n"
  )
}

get_permuted_mean_matrix <- function(estimator) {
  matrix_data <- permuted_analysis$mean_covariance_matrices_long[
    permuted_analysis$mean_covariance_matrices_long$threshold == 0.95 &
      permuted_analysis$mean_covariance_matrices_long$se_scale ==
        "t_adjusted" &
      permuted_analysis$mean_covariance_matrices_long$estimator == estimator,
  ]
  matrix(
    matrix_data$covariance,
    nrow = length(permuted_time_grid),
    ncol = length(permuted_time_grid)
  )
}

plot_permuted_covariance_matrices <- function() {
  matrices <- lapply(permuted_method_order, get_permuted_mean_matrix)
  names(matrices) <- permuted_method_order
  z_range <- range(c(0, unlist(matrices)))
  palette <- grDevices::hcl.colors(201L, "Reds 3", rev = TRUE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::layout(matrix(1:4, nrow = 1L), widths = c(1, 1, 1, 0.20))
  graphics::par(oma = c(0, 0, 2.7, 0))
  for (estimator in permuted_method_order) {
    graphics::par(mar = c(4.2, 4.2, 3.2, 1.0))
    graphics::image(
      x = permuted_time_grid,
      y = permuted_time_grid,
      z = matrices[[estimator]],
      col = palette,
      zlim = z_range,
      xlab = "Time",
      ylab = "Time",
      main = permuted_method_labels[estimator],
      axes = FALSE,
      useRaster = TRUE
    )
    graphics::axis(1, at = c(0, 5, 10, 15))
    graphics::axis(2, at = c(0, 5, 10, 15))
    graphics::box()
  }
  graphics::par(mar = c(4.2, 0.2, 3.2, 2.5))
  color_values <- seq(
    z_range[1L],
    z_range[2L],
    length.out = length(palette)
  )
  graphics::image(
    x = c(0, 1),
    y = color_values,
    z = matrix(rep(color_values, each = 2L), nrow = 2L),
    col = palette,
    xlab = "",
    ylab = "",
    axes = FALSE,
    useRaster = TRUE
  )
  graphics::axis(4, las = 1)
  graphics::box()
  graphics::mtext(
    paste(
      "Mean empirical covariance matrices from matched permutations;",
      "lfdr > 0.95, t-adjusted SE"
    ),
    outer = TRUE,
    line = 0.7,
    cex = 1.05
  )
}

plot_permuted_direct_covariance_matrix <- function() {
  palette <- grDevices::colorRampPalette(
    c("white", "#FEE8C8", "#F46D43", "#A50026")
  )(201L)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(5.0, 5.0, 4.2, 1.5))
  graphics::image(
    x = permuted_time_grid,
    y = permuted_time_grid,
    z = permuted_primary_matrix,
    col = palette,
    zlim = c(0, max(permuted_primary_matrix)),
    xlab = "Time",
    ylab = "Time",
    main = paste(
      "Direct zero-null empirical covariance matrix\n",
      "matched donor permutations; lfdr > 0.95, t-adjusted SE"
    ),
    axes = FALSE,
    useRaster = TRUE
  )
  graphics::axis(1, at = permuted_time_grid)
  graphics::axis(2, at = permuted_time_grid)
  graphics::box()
  for (row_index in seq_along(permuted_time_grid)) {
    for (column_index in seq_along(permuted_time_grid)) {
      value <- permuted_primary_matrix[column_index, row_index]
      text_color <- if (
        value > 0.55 * max(permuted_primary_matrix)
      ) {
        "white"
      } else {
        "black"
      }
      graphics::text(
        permuted_time_grid[row_index],
        permuted_time_grid[column_index],
        labels = formatC(value, format = "f", digits = 2L),
        cex = 0.66,
        col = text_color
      )
    }
  }
}
