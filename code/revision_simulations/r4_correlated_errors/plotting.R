# Plotting helpers for the R4 paired correlated-error analysis.

r4_method_colors <- function() {
  c(
    "FASH-IWP1-Raw" = "#1B9E77",
    "FASH-IWP1-BF" = "#0072B2",
    "FASH-linear-Raw" = "#D55E00",
    "FASH-linear-BF" = "#CC79A7"
  )
}

r4_condition_line_types <- function() {
  c("Independent" = 2, "Correlated" = 1)
}

plot_r4_condition_curves <- function(mc_alpha,
                                     metric = c("power", "fdr"),
                                     file = NULL,
                                     title = NULL,
                                     width = 1900,
                                     height = 850,
                                     res = 170) {
  metric <- match.arg(metric)
  required_columns <- c(
    "condition", "method", "alpha",
    if (metric == "power") "mean_power" else "mean_fdr"
  )
  missing_columns <- setdiff(required_columns, names(mc_alpha))
  if (length(missing_columns) > 0) {
    stop(
      "mc_alpha is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  methods <- names(r4_method_colors())
  mc_alpha <- mc_alpha[
    mc_alpha$method %in% methods &
      mc_alpha$condition %in% names(r4_condition_line_types()),
    ,
    drop = FALSE
  ]
  if (nrow(mc_alpha) == 0) {
    stop("No R4 condition curves are available to plot.")
  }

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(1, 2),
    mar = c(4.4, 4.5, 3.2, 1.0),
    oma = c(0, 0, 2.4, 0)
  )

  method_groups <- list(
    "IWP1 alternative" = c("FASH-IWP1-Raw", "FASH-IWP1-BF"),
    "Linear alternative" = c("FASH-linear-Raw", "FASH-linear-BF")
  )
  value_column <- if (metric == "power") "mean_power" else "mean_fdr"
  y_label <- if (metric == "power") "Mean power" else "Mean empirical FDR"
  colors <- r4_method_colors()
  line_types <- r4_condition_line_types()

  for (panel_name in names(method_groups)) {
    panel_methods <- method_groups[[panel_name]]
    panel <- mc_alpha[mc_alpha$method %in% panel_methods, , drop = FALSE]
    y_max <- if (metric == "power") {
      1
    } else {
      max(0.12, panel[[value_column]], panel$alpha, na.rm = TRUE)
    }
    graphics::plot(
      NA,
      xlim = range(panel$alpha, na.rm = TRUE),
      ylim = c(0, y_max),
      xlab = "Nominal FDR level alpha",
      ylab = y_label,
      main = panel_name
    )
    graphics::grid(col = "gray90")
    if (metric == "fdr") {
      graphics::abline(a = 0, b = 1, col = "gray45", lty = 3, lwd = 1.5)
    }
    graphics::abline(v = 0.05, col = "gray65", lty = 3)

    for (method in panel_methods) {
      for (condition in names(line_types)) {
        curve <- panel[
          panel$method == method & panel$condition == condition,
          ,
          drop = FALSE
        ]
        curve <- curve[order(curve$alpha), , drop = FALSE]
        graphics::lines(
          curve$alpha,
          curve[[value_column]],
          col = colors[[method]],
          lty = line_types[[condition]],
          lwd = 2.5
        )
      }
    }

    graphics::legend(
      if (metric == "power") "bottomright" else "topleft",
      legend = c(panel_methods, names(line_types)),
      col = c(colors[panel_methods], rep("gray25", length(line_types))),
      lty = c(rep(1, length(panel_methods)), unname(line_types)),
      lwd = 2.5,
      bty = "n",
      cex = 0.78
    )
  }

  if (is.null(title)) {
    title <- if (metric == "power") {
      "R4 power under independent and correlated errors"
    } else {
      "R4 empirical FDR under independent and correlated errors"
    }
  }
  graphics::mtext(title, outer = TRUE, cex = 1.15, font = 2)
  invisible(mc_alpha)
}

plot_r4_lag_diagnostics <- function(lag_summary,
                                    file = NULL,
                                    title = "Realized time dependence in R4",
                                    empirical_centered_target = -0.152,
                                    width = 1500,
                                    height = 1000,
                                    res = 170) {
  required_columns <- c(
    "condition", "diagnostic", "lag", "mean_correlation",
    "ci_lower", "ci_upper"
  )
  missing_columns <- setdiff(required_columns, names(lag_summary))
  if (length(missing_columns) > 0) {
    stop(
      "lag_summary is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  diagnostic_order <- c(
    "Truth-known standardized error",
    "Centered dynamic-null residual"
  )
  colors <- c(
    "Truth-known standardized error" = "#0072B2",
    "Centered dynamic-null residual" = "#D55E00"
  )
  line_types <- r4_condition_line_types()
  y_range <- range(
    lag_summary$ci_lower,
    lag_summary$ci_upper,
    empirical_centered_target,
    0,
    na.rm = TRUE
  )
  y_padding <- 0.05 * diff(y_range)

  graphics::plot(
    NA,
    xlim = range(lag_summary$lag),
    ylim = y_range + c(-y_padding, y_padding),
    xlab = "Time lag",
    ylab = "Mean correlation",
    main = title
  )
  graphics::grid(col = "gray90")
  graphics::abline(h = 0, col = "gray55", lty = 3)
  graphics::abline(
    h = empirical_centered_target,
    col = "#B2182B",
    lty = 3,
    lwd = 1.6
  )

  for (diagnostic in diagnostic_order) {
    for (condition in names(line_types)) {
      curve <- lag_summary[
        lag_summary$diagnostic == diagnostic &
          lag_summary$condition == condition,
        ,
        drop = FALSE
      ]
      curve <- curve[order(curve$lag), , drop = FALSE]
      graphics::polygon(
        c(curve$lag, rev(curve$lag)),
        c(curve$ci_lower, rev(curve$ci_upper)),
        col = grDevices::adjustcolor(colors[[diagnostic]], alpha.f = 0.08),
        border = NA
      )
      graphics::lines(
        curve$lag,
        curve$mean_correlation,
        col = colors[[diagnostic]],
        lty = line_types[[condition]],
        lwd = 2.4
      )
    }
  }

  graphics::legend(
    "topright",
    legend = c(
      diagnostic_order,
      names(line_types),
      "Real-data centered lag-1 target"
    ),
    col = c(
      colors[diagnostic_order],
      rep("gray25", length(line_types)),
      "#B2182B"
    ),
    lty = c(
      rep(1, length(diagnostic_order)),
      unname(line_types),
      3
    ),
    lwd = 2.4,
    bty = "n",
    cex = 0.78
  )
  invisible(lag_summary)
}

plot_r4_correlation_sweep <- function(sweep_summary,
                                      file = NULL,
                                      title = paste(
                                        "R4 lag-1 correlation sweep",
                                        "at alpha = 0.05"
                                      ),
                                      empirical_rho = -0.09,
                                      nominal_alpha = 0.05,
                                      methods = names(r4_method_colors()),
                                      width = 1900,
                                      height = 1450,
                                      res = 170) {
  required_columns <- c(
    "rho", "method", "mean_power", "power_ci_lower", "power_ci_upper",
    "mean_fdr", "fdr_ci_lower", "fdr_ci_upper"
  )
  missing_columns <- setdiff(required_columns, names(sweep_summary))
  if (length(missing_columns) > 0L) {
    stop(
      "sweep_summary is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  methods <- as.character(methods)
  available_methods <- names(r4_method_colors())
  if (length(methods) == 0L || anyDuplicated(methods) ||
      any(!methods %in% available_methods)) {
    stop("The requested R4 sweep methods are invalid.")
  }
  sweep_summary <- sweep_summary[
    sweep_summary$method %in% methods,
    ,
    drop = FALSE
  ]
  expected_rho <- sort(unique(sweep_summary$rho))
  if (length(expected_rho) < 3L ||
      any(!is.finite(expected_rho)) ||
      !any(abs(expected_rho) < 1e-12)) {
    stop("The R4 correlation sweep must contain at least three finite rho values including zero.")
  }
  method_rho_counts <- table(sweep_summary$method, sweep_summary$rho)
  if (!identical(sort(unique(sweep_summary$method)), sort(methods)) ||
      any(method_rho_counts != 1L)) {
    stop("The R4 correlation sweep is incomplete across methods and rho values.")
  }

  all_method_families <- list(
    "IWP1 alternative" = c("FASH-IWP1-Raw", "FASH-IWP1-BF"),
    "Linear alternative" = c("FASH-linear-Raw", "FASH-linear-BF")
  )
  method_families <- lapply(all_method_families, function(family_methods) {
    intersect(family_methods, methods)
  })
  method_families <- method_families[lengths(method_families) > 0L]

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    plot_height <- if (length(method_families) == 1L) min(height, 850) else height
    grDevices::png(file, width = width, height = plot_height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(length(method_families), 2),
    mar = c(4.2, 4.6, 3.0, 1.0),
    oma = c(0, 0, 2.8, 0)
  )

  metric_specs <- list(
    power = list(
      mean = "mean_power",
      lower = "power_ci_lower",
      upper = "power_ci_upper",
      label = "Mean power"
    ),
    fdr = list(
      mean = "mean_fdr",
      lower = "fdr_ci_lower",
      upper = "fdr_ci_upper",
      label = "Mean empirical FDR"
    )
  )
  colors <- r4_method_colors()

  for (family_name in names(method_families)) {
    family_methods <- method_families[[family_name]]
    family_rows <- sweep_summary[
      sweep_summary$method %in% family_methods,
      ,
      drop = FALSE
    ]

    for (metric_name in names(metric_specs)) {
      spec <- metric_specs[[metric_name]]
      if (metric_name == "power") {
        observed_range <- range(
          family_rows[[spec$lower]],
          family_rows[[spec$upper]],
          na.rm = TRUE
        )
        padding <- max(0.025, 0.12 * diff(observed_range))
        y_limits <- c(
          max(0, observed_range[1] - padding),
          min(1, observed_range[2] + padding)
        )
      } else {
        y_limits <- c(
          0,
          max(
            0.10,
            nominal_alpha * 1.35,
            family_rows[[spec$upper]],
            na.rm = TRUE
          )
        )
      }
      if (!all(is.finite(y_limits)) || diff(y_limits) <= 0) {
        stop("Could not construct finite plotting limits for the R4 sweep.")
      }

      graphics::plot(
        NA,
        xlim = range(expected_rho),
        ylim = y_limits,
        xlab = "Generating lag-1 correlation (rho)",
        ylab = spec$label,
        main = paste(family_name, spec$label, sep = ": ")
      )
      graphics::grid(col = "gray90")

      for (method in family_methods) {
        curve <- family_rows[
          family_rows$method == method,
          ,
          drop = FALSE
        ]
        curve <- curve[order(curve$rho), , drop = FALSE]
        graphics::polygon(
          c(curve$rho, rev(curve$rho)),
          c(curve[[spec$lower]], rev(curve[[spec$upper]])),
          col = grDevices::adjustcolor(colors[[method]], alpha.f = 0.12),
          border = NA
        )
        graphics::lines(
          curve$rho,
          curve[[spec$mean]],
          col = colors[[method]],
          lwd = 2.5,
          type = "o",
          pch = if (grepl("-BF$", method)) 16 else 1
        )
      }

      graphics::abline(v = 0, col = "gray35", lty = 2, lwd = 1.4)
      graphics::abline(
        v = empirical_rho,
        col = "#B2182B",
        lty = 3,
        lwd = 1.7
      )
      if (metric_name == "fdr") {
        graphics::abline(
          h = nominal_alpha,
          col = "gray20",
          lty = 3,
          lwd = 1.5
        )
      }
      graphics::legend(
        if (metric_name == "power") "bottomleft" else "topleft",
        legend = c(
          sub("FASH-", "", family_methods, fixed = TRUE),
          "rho = 0",
          paste0("Empirical rho = ", format(empirical_rho, nsmall = 2)),
          if (metric_name == "fdr") "Nominal FDR = 0.05" else NULL
        ),
        col = c(
          colors[family_methods],
          "gray35",
          "#B2182B",
          if (metric_name == "fdr") "gray20" else NULL
        ),
        lty = c(
          rep(1, length(family_methods)),
          2,
          3,
          if (metric_name == "fdr") 3 else NULL
        ),
        pch = c(
          ifelse(grepl("-BF$", family_methods), 16, 1),
          NA,
          NA,
          if (metric_name == "fdr") NA else NULL
        ),
        lwd = c(
          rep(2.5, length(family_methods)),
          1.4,
          1.7,
          if (metric_name == "fdr") 1.5 else NULL
        ),
        bty = "n",
        cex = 0.72
      )
    }
  }

  graphics::mtext(title, outer = TRUE, cex = 1.15, font = 2)
  invisible(sweep_summary)
}

plot_r4_correlation_sweep_pi0 <- function(pi0_summary,
                                          file = NULL,
                                          title = paste(
                                            "Estimated dynamic-null proportion",
                                            "across lag-1 correlations"
                                          ),
                                          empirical_rho = -0.09,
                                          true_pi0 = 0.8,
                                          width = 1500,
                                          height = 900,
                                          res = 170) {
  required_columns <- c(
    "rho", "method", "mean_estimated_pi0", "pi0_ci_lower", "pi0_ci_upper"
  )
  missing_columns <- setdiff(required_columns, names(pi0_summary))
  if (length(missing_columns) > 0L) {
    stop(
      "pi0_summary is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  methods <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
  pi0_summary <- pi0_summary[pi0_summary$method %in% methods, , drop = FALSE]
  method_rho_counts <- table(pi0_summary$method, pi0_summary$rho)
  if (!identical(sort(unique(pi0_summary$method)), sort(methods)) ||
      any(method_rho_counts != 1L) || any(!is.finite(pi0_summary$rho))) {
    stop("The IWP1 pi0 sweep is incomplete.")
  }

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mar = c(4.4, 4.8, 3.4, 1.0))

  y_range <- range(
    pi0_summary$pi0_ci_lower,
    pi0_summary$pi0_ci_upper,
    true_pi0,
    na.rm = TRUE
  )
  padding <- max(0.015, 0.12 * diff(y_range))
  y_limits <- c(max(0, y_range[1] - padding), min(1, y_range[2] + padding))
  graphics::plot(
    NA,
    xlim = range(pi0_summary$rho),
    ylim = y_limits,
    xlab = "Generating lag-1 correlation (rho)",
    ylab = "Estimated dynamic-null proportion",
    main = title
  )
  graphics::grid(col = "gray90")
  colors <- r4_method_colors()
  for (method in methods) {
    curve <- pi0_summary[pi0_summary$method == method, , drop = FALSE]
    curve <- curve[order(curve$rho), , drop = FALSE]
    graphics::polygon(
      c(curve$rho, rev(curve$rho)),
      c(curve$pi0_ci_lower, rev(curve$pi0_ci_upper)),
      col = grDevices::adjustcolor(colors[[method]], alpha.f = 0.12),
      border = NA
    )
    graphics::lines(
      curve$rho,
      curve$mean_estimated_pi0,
      col = colors[[method]],
      lwd = 2.5,
      type = "o",
      pch = if (grepl("-BF$", method)) 16 else 1
    )
  }
  graphics::abline(h = true_pi0, col = "gray20", lty = 3, lwd = 1.6)
  graphics::abline(v = 0, col = "gray35", lty = 2, lwd = 1.4)
  graphics::abline(v = empirical_rho, col = "#B2182B", lty = 3, lwd = 1.7)
  graphics::legend(
    "topright",
    legend = c(
      "IWP1 Raw",
      "IWP1 BF-corrected",
      paste0("True pi0 = ", format(true_pi0, nsmall = 1)),
      "rho = 0",
      paste0("Empirical rho = ", format(empirical_rho, nsmall = 2))
    ),
    col = c(colors[methods], "gray20", "gray35", "#B2182B"),
    lty = c(1, 1, 3, 2, 3),
    pch = c(1, 16, NA, NA, NA),
    lwd = c(2.5, 2.5, 1.6, 1.4, 1.7),
    bty = "n",
    cex = 0.82
  )
  invisible(pi0_summary)
}

r4_correlation_palette <- function(n = 201L) {
  grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(n)
}

draw_r4_correlation_heatmap <- function(correlation,
                                        title,
                                        show_legend = FALSE,
                                        zlim = c(-1, 1)) {
  correlation <- as.matrix(correlation)
  if (nrow(correlation) != ncol(correlation) || any(!is.finite(correlation))) {
    stop("A correlation heatmap requires one finite square matrix.")
  }
  n_time <- ncol(correlation)
  time_grid <- seq_len(n_time) - 1L
  colors <- r4_correlation_palette()
  graphics::image(
    x = time_grid,
    y = time_grid,
    z = correlation,
    zlim = zlim,
    col = colors,
    xlab = "Time",
    ylab = "Time",
    main = title,
    axes = FALSE,
    useRaster = TRUE
  )
  graphics::axis(1, at = time_grid, labels = time_grid, cex.axis = 0.72)
  graphics::axis(2, at = time_grid, labels = time_grid, las = 1, cex.axis = 0.72)
  graphics::box()
  if (show_legend) {
    legend_values <- seq(zlim[1], zlim[2], length.out = 101L)
    legend_x <- max(time_grid) + c(1.2, 2.0)
    legend_y <- seq(min(time_grid), max(time_grid), length.out = 101L)
    graphics::image(
      x = legend_x,
      y = legend_y,
      z = matrix(legend_values, nrow = 2L, ncol = length(legend_values), byrow = TRUE),
      zlim = zlim,
      col = colors,
      add = TRUE
    )
  }
  invisible(correlation)
}

plot_r4_real_data_correlation_heatmaps <- function(real_analysis,
                                                   file = NULL,
                                                   width = 1800,
                                                   height = 1600,
                                                   res = 180) {
  required <- c("direct_adjusted", "pairwise_adjusted")
  if (is.null(real_analysis$matrices) ||
      !all(required %in% names(real_analysis$matrices))) {
    stop("The real-data analysis does not contain both correlation estimators.")
  }
  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(2, 2),
    mar = c(4.0, 4.0, 3.2, 1.0),
    oma = c(1.8, 0, 2.7, 0)
  )
  draw_r4_correlation_heatmap(
    real_analysis$matrices$direct_adjusted$raw,
    "Direct centered: raw estimate"
  )
  draw_r4_correlation_heatmap(
    real_analysis$matrices$direct_adjusted$projected,
    "Direct centered: nearest-PD"
  )
  draw_r4_correlation_heatmap(
    real_analysis$matrices$pairwise_adjusted$raw,
    "Pairwise difference: raw estimate"
  )
  draw_r4_correlation_heatmap(
    real_analysis$matrices$pairwise_adjusted$projected,
    "Pairwise difference: nearest-PD"
  )
  graphics::mtext(
    "Full 16 by 16 correlation estimates from 500 null-like gene-variant pairs",
    outer = TRUE,
    cex = 1.1,
    font = 2
  )
  graphics::mtext(
    "Common color scale: blue = -1, white = 0, red = +1",
    side = 1,
    outer = TRUE,
    line = 0.2,
    cex = 0.82
  )
  invisible(real_analysis)
}

plot_r4_real_data_variograms <- function(bootstrap_lags,
                                         independence_benchmark,
                                         file = NULL,
                                         width = 1900,
                                         height = 850,
                                         res = 180) {
  observed_required <- c(
    "estimator", "lag", "semivariogram", "semivariogram_ci_lower",
    "semivariogram_ci_upper"
  )
  benchmark_required <- c(
    "estimator", "benchmark", "lag", "mean_semivariogram",
    "semivariogram_ci_lower", "semivariogram_ci_upper"
  )
  if (length(setdiff(observed_required, names(bootstrap_lags))) > 0L ||
      length(setdiff(benchmark_required, names(independence_benchmark))) > 0L) {
    stop("The real-data variogram inputs are incomplete.")
  }
  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(1, 2),
    mar = c(4.4, 4.7, 3.3, 1.0),
    oma = c(0, 0, 2.7, 0)
  )
  estimators <- c("Direct centered", "Pairwise difference")
  benchmark_types <- c(
    "Fixed selected SE patterns",
    "Flatness-selected gene representatives"
  )
  benchmark_colors <- c(
    "Fixed selected SE patterns" = "#4D4D4D",
    "Flatness-selected gene representatives" = "#E69F00"
  )
  benchmark_line_types <- c(
    "Fixed selected SE patterns" = 2,
    "Flatness-selected gene representatives" = 3
  )

  for (estimator in estimators) {
    observed <- bootstrap_lags[
      bootstrap_lags$estimator == estimator,
      ,
      drop = FALSE
    ]
    observed <- observed[order(observed$lag), , drop = FALSE]
    benchmarks <- independence_benchmark[
      independence_benchmark$estimator == estimator,
      ,
      drop = FALSE
    ]
    y_range <- range(
      observed$semivariogram_ci_lower,
      observed$semivariogram_ci_upper,
      benchmarks$semivariogram_ci_lower,
      benchmarks$semivariogram_ci_upper,
      na.rm = TRUE
    )
    padding <- max(0.025, 0.08 * diff(y_range))
    graphics::plot(
      NA,
      xlim = range(observed$lag),
      ylim = y_range + c(-padding, padding),
      xlab = "Time lag",
      ylab = "Standardized semivariogram",
      main = estimator
    )
    graphics::grid(col = "gray90")
    graphics::polygon(
      c(observed$lag, rev(observed$lag)),
      c(
        observed$semivariogram_ci_lower,
        rev(observed$semivariogram_ci_upper)
      ),
      col = grDevices::adjustcolor("#0072B2", alpha.f = 0.15),
      border = NA
    )
    graphics::lines(
      observed$lag,
      observed$semivariogram,
      col = "#0072B2",
      lwd = 2.7,
      type = "o",
      pch = 16
    )
    for (benchmark in benchmark_types) {
      curve <- benchmarks[benchmarks$benchmark == benchmark, , drop = FALSE]
      curve <- curve[order(curve$lag), , drop = FALSE]
      graphics::polygon(
        c(curve$lag, rev(curve$lag)),
        c(
          curve$semivariogram_ci_lower,
          rev(curve$semivariogram_ci_upper)
        ),
        col = grDevices::adjustcolor(
          benchmark_colors[[benchmark]],
          alpha.f = 0.07
        ),
        border = NA
      )
      graphics::lines(
        curve$lag,
        curve$mean_semivariogram,
        col = benchmark_colors[[benchmark]],
        lty = benchmark_line_types[[benchmark]],
        lwd = 2.0
      )
    }
    graphics::legend(
      if (estimator == "Direct centered") "bottomright" else "topright",
      legend = c(
        "Observed top-500 (95% gene bootstrap)",
        "Independent, fixed SE patterns",
        "Independent, flatness-selected"
      ),
      col = c("#0072B2", benchmark_colors[benchmark_types]),
      lty = c(1, unname(benchmark_line_types[benchmark_types])),
      lwd = c(2.7, 2.0, 2.0),
      pch = c(16, NA, NA),
      bty = "n",
      cex = 0.76
    )
  }
  graphics::mtext(
    "Real-data variograms and matched independence benchmarks",
    outer = TRUE,
    cex = 1.1,
    font = 2
  )
  invisible(bootstrap_lags)
}

plot_r4_full_matrix_target_realized <- function(configuration,
                                                matrix_summary,
                                                file = NULL,
                                                width = 1800,
                                                height = 1600,
                                                res = 180) {
  conditions <- c(
    "Direct centered full matrix",
    "Pairwise-difference full matrix"
  )
  if (is.null(configuration$condition_matrices) ||
      !all(conditions %in% names(configuration$condition_matrices))) {
    stop("The full-matrix configuration is incomplete.")
  }
  realized_rows <- matrix_summary[
    matrix_summary$condition %in% conditions &
      matrix_summary$diagnostic == "Expression error",
    ,
    drop = FALSE
  ]
  if (nrow(realized_rows) != length(conditions) * 16^2) {
    stop("The mean realized expression-error matrices are incomplete.")
  }
  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(2, 2),
    mar = c(4.0, 4.0, 3.2, 1.0),
    oma = c(1.8, 0, 2.7, 0)
  )
  for (condition in conditions) {
    short_name <- if (grepl("Direct", condition)) {
      "Direct centered"
    } else {
      "Pairwise difference"
    }
    target <- configuration$condition_matrices[[condition]]
    rows <- realized_rows[realized_rows$condition == condition, , drop = FALSE]
    rows <- rows[order(rows$time_a, rows$time_b), , drop = FALSE]
    realized <- matrix(NA_real_, nrow = 16L, ncol = 16L)
    realized[cbind(rows$time_a + 1L, rows$time_b + 1L)] <-
      rows$mean_correlation
    draw_r4_correlation_heatmap(target, paste(short_name, "target"))
    draw_r4_correlation_heatmap(
      realized,
      paste(short_name, "mean realized expression error")
    )
  }
  graphics::mtext(
    "Target and realized full error-correlation matrices",
    outer = TRUE,
    cex = 1.1,
    font = 2
  )
  graphics::mtext(
    "Common color scale: blue = -1, white = 0, red = +1",
    side = 1,
    outer = TRUE,
    line = 0.2,
    cex = 0.82
  )
  invisible(matrix_summary)
}

plot_r4_full_matrix_metrics <- function(condition_alpha,
                                        condition_pi0,
                                        file = NULL,
                                        true_pi0 = 0.8,
                                        nominal_alpha = 0.05,
                                        width = 2100,
                                        height = 850,
                                        res = 180) {
  methods <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
  conditions <- c(
    "Independent",
    "Direct centered full matrix",
    "Pairwise-difference full matrix"
  )
  alpha_rows <- condition_alpha[
    condition_alpha$method %in% methods &
      condition_alpha$condition %in% conditions,
    ,
    drop = FALSE
  ]
  pi0_rows <- condition_pi0[
    condition_pi0$method %in% methods &
      condition_pi0$condition %in% conditions,
    ,
    drop = FALSE
  ]
  if (nrow(alpha_rows) != 6L || nrow(pi0_rows) != 6L) {
    stop("The full-matrix IWP1 metric summaries are incomplete.")
  }
  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(file, width = width, height = height, res = res)
    on.exit(grDevices::dev.off(), add = TRUE)
  }
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(1, 3),
    mar = c(7.0, 4.7, 3.3, 1.0),
    oma = c(0, 0, 2.6, 0)
  )
  colors <- r4_method_colors()
  condition_labels <- c(
    "Independent" = "Independent",
    "Direct centered full matrix" = "Direct full",
    "Pairwise-difference full matrix" = "Pairwise full"
  )
  metric_specs <- list(
    power = list(
      data = alpha_rows,
      mean = "mean_power",
      lower = "power_ci_lower",
      upper = "power_ci_upper",
      label = "Power",
      reference = NA_real_
    ),
    fdr = list(
      data = alpha_rows,
      mean = "mean_fdr",
      lower = "fdr_ci_lower",
      upper = "fdr_ci_upper",
      label = "Empirical FDR",
      reference = nominal_alpha
    ),
    pi0 = list(
      data = pi0_rows,
      mean = "mean_estimated_pi0",
      lower = "pi0_ci_lower",
      upper = "pi0_ci_upper",
      label = "Estimated pi0",
      reference = true_pi0
    )
  )
  for (metric_name in names(metric_specs)) {
    spec <- metric_specs[[metric_name]]
    panel <- spec$data
    x_positions <- seq_along(conditions)
    y_range <- range(
      panel[[spec$lower]],
      panel[[spec$upper]],
      spec$reference,
      na.rm = TRUE
    )
    padding <- max(0.02, 0.12 * diff(y_range))
    lower_limit <- if (metric_name == "fdr") {
      y_range[1] - padding
    } else {
      max(0, y_range[1] - padding)
    }
    graphics::plot(
      NA,
      xlim = c(0.6, length(conditions) + 0.4),
      ylim = c(lower_limit, min(1, y_range[2] + padding)),
      xaxt = "n",
      xlab = "",
      ylab = spec$label,
      main = spec$label
    )
    graphics::axis(
      1,
      at = x_positions,
      labels = unname(condition_labels[conditions]),
      las = 2,
      cex.axis = 0.82
    )
    graphics::grid(nx = NA, ny = NULL, col = "gray90")
    for (method_index in seq_along(methods)) {
      method <- methods[method_index]
      offset <- if (method_index == 1L) -0.08 else 0.08
      rows <- panel[panel$method == method, , drop = FALSE]
      rows$condition <- factor(rows$condition, levels = conditions)
      rows <- rows[order(rows$condition), , drop = FALSE]
      x <- x_positions + offset
      graphics::segments(
        x,
        rows[[spec$lower]],
        x,
        rows[[spec$upper]],
        col = colors[[method]],
        lwd = 2
      )
      graphics::points(
        x,
        rows[[spec$mean]],
        col = colors[[method]],
        pch = if (grepl("-BF$", method)) 16 else 1,
        cex = 1.15,
        lwd = 2
      )
    }
    if (is.finite(spec$reference)) {
      graphics::abline(h = spec$reference, col = "gray25", lty = 3, lwd = 1.5)
    }
    graphics::legend(
      "bottomleft",
      legend = c("Raw", "BF-corrected"),
      col = colors[methods],
      pch = c(1, 16),
      bty = "n",
      cex = 0.78
    )
  }
  graphics::mtext(
    "IWP1 performance under complete empirical correlation patterns",
    outer = TRUE,
    cex = 1.1,
    font = 2
  )
  invisible(condition_alpha)
}
