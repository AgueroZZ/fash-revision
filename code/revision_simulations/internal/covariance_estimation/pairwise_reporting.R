# Reporting helpers for the pairwise lfdr-threshold correlation exploration.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

format_decimal <- function(x, digits = 3L) {
  formatC(x, format = "f", digits = digits)
}

format_interval <- function(estimate, lower, upper, digits = 3L) {
  paste0(
    format_decimal(estimate, digits),
    " [",
    format_decimal(lower, digits),
    ", ",
    format_decimal(upper, digits),
    "]"
  )
}

render_scrollable_table <- function(x,
                                    align = NULL,
                                    minimum_width = "1000px",
                                    digits = 3L) {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("The knitr package is required to render tables.")
  }
  table_html <- knitr::kable(
    x,
    format = "html",
    align = align,
    digits = digits,
    escape = TRUE,
    row.names = FALSE
  )
  cat(
    paste0(
      '<div class="internal-table-scroll" style="min-width:0;">',
      '<div style="min-width:', minimum_width, ';">',
      table_html,
      "</div></div>\n"
    )
  )
  invisible(x)
}

threshold_key <- function(threshold) {
  paste0(
    "lfdr_",
    sub("\\.", "p", format(threshold, nsmall = 2, trim = TRUE))
  )
}

workflowr_root <- find_workflowr_root()
cache_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "pairwise_lfdr_threshold_correlation"
)
analysis_path <- file.path(
  cache_dir,
  "pairwise_lfdr_threshold_analysis.rds"
)
if (!file.exists(analysis_path)) {
  stop("The pairwise lfdr-threshold analysis cache is missing: ", analysis_path)
}
pairwise_analysis <- readRDS(analysis_path)
configuration <- pairwise_analysis$configuration
expected_thresholds <- c(0.97, 0.96, 0.95)
expected_matrix_names <- vapply(
  expected_thresholds,
  threshold_key,
  character(1)
)

required_names <- c(
  "configuration",
  "selections",
  "matrices",
  "bootstrap_results",
  "selected_units",
  "selection_counts",
  "matrix_diagnostics",
  "correlation_matrices_long",
  "lag_variogram_summaries",
  "bootstrap_lag_variogram_intervals"
)
if (!all(required_names %in% names(pairwise_analysis)) ||
    !identical(configuration$analysis_id, "pairwise_lfdr_threshold_correlation") ||
    !isTRUE(all.equal(configuration$thresholds, expected_thresholds)) ||
    !identical(names(pairwise_analysis$matrices), expected_matrix_names) ||
    !identical(names(pairwise_analysis$selections), expected_matrix_names) ||
    configuration$n_bootstrap != 1000L ||
    !identical(configuration$matrix_version,
               "Raw pairwise-difference estimate; no PD projection")) {
  stop("The pairwise lfdr-threshold cache configuration is incomplete.")
}

fit_path <- file.path(workflowr_root, configuration$fit_path)
if (!file.exists(fit_path)) {
  stop("The source FASH fit recorded by the cache is missing.")
}
fit_info <- file.info(fit_path)
if (!isTRUE(all.equal(unname(fit_info$size), configuration$fit_size_bytes))) {
  stop("The source FASH fit size no longer matches the cache metadata.")
}

for (threshold_index in seq_along(expected_thresholds)) {
  name <- expected_matrix_names[threshold_index]
  threshold <- expected_thresholds[threshold_index]
  correlation <- pairwise_analysis$matrices[[name]]
  selection <- pairwise_analysis$selections[[name]]$selected
  if (!is.matrix(correlation) || !identical(dim(correlation), c(16L, 16L)) ||
      any(!is.finite(correlation)) ||
      max(abs(correlation - t(correlation))) > 1e-12 ||
      max(abs(diag(correlation) - 1)) > 1e-12 ||
      nrow(selection) < 20L || any(selection$lfdr <= threshold) ||
      anyDuplicated(selection$gene_id) || anyDuplicated(selection$pair_key)) {
    stop("A threshold-specific raw matrix or selection is invalid.")
  }
}
for (threshold_index in seq_len(length(expected_thresholds) - 1L)) {
  tighter <- pairwise_analysis$selections[[threshold_index]]$selected$pair_key
  looser <- pairwise_analysis$selections[[threshold_index + 1L]]$selected$pair_key
  if (!all(tighter %in% looser)) {
    stop("The cached lfdr-threshold selections are not nested.")
  }
}

selection_counts <- pairwise_analysis$selection_counts
matrix_diagnostics <- pairwise_analysis$matrix_diagnostics
lag_summaries <- pairwise_analysis$lag_variogram_summaries
bootstrap_intervals <- pairwise_analysis$bootstrap_lag_variogram_intervals
if (!isTRUE(all.equal(selection_counts$threshold, expected_thresholds)) ||
    !isTRUE(all.equal(matrix_diagnostics$threshold, expected_thresholds)) ||
    nrow(lag_summaries) != 45L ||
    nrow(bootstrap_intervals) != 45L ||
    any(table(lag_summaries$threshold) != 15L) ||
    any(table(bootstrap_intervals$threshold) != 15L) ||
    any(bootstrap_intervals$correlation_ci_lower >
        bootstrap_intervals$correlation_ci_upper) ||
    any(bootstrap_intervals$semivariogram_ci_lower >
        bootstrap_intervals$semivariogram_ci_upper)) {
  stop("The cached threshold summaries or bootstrap intervals are invalid.")
}

get_lag_interval <- function(threshold, lag) {
  hit <- bootstrap_intervals[
    bootstrap_intervals$threshold == threshold &
      bootstrap_intervals$lag == lag,
    ,
    drop = FALSE
  ]
  if (nrow(hit) != 1L) {
    stop("Could not identify one threshold-specific lag interval.")
  }
  hit
}

summary_rows <- lapply(expected_thresholds, function(threshold) {
  count <- selection_counts[selection_counts$threshold == threshold, ]
  diagnostics <- matrix_diagnostics[
    matrix_diagnostics$threshold == threshold,
  ]
  lag1 <- get_lag_interval(threshold, 1L)
  lag15 <- get_lag_interval(threshold, 15L)
  data.frame(
    `lfdr rule` = paste0("lfdr > ", format_decimal(threshold, 2)),
    `Selected genes` = count$n_selected,
    `Selected lfdr range` = paste0(
      "[",
      format_decimal(count$minimum_selected_lfdr, 4),
      ", ",
      format_decimal(count$maximum_selected_lfdr, 4),
      "]"
    ),
    `Mean lag-1 correlation (95% CI)` = format_interval(
      lag1$observed_correlation,
      lag1$correlation_ci_lower,
      lag1$correlation_ci_upper
    ),
    `Mean lag-15 correlation (95% CI)` = format_interval(
      lag15$observed_correlation,
      lag15$correlation_ci_lower,
      lag15$correlation_ci_upper
    ),
    `Mean off-diagonal correlation` = format_decimal(
      diagnostics$mean_off_diagonal_correlation
    ),
    `Off-diagonal range` = paste0(
      "[",
      format_decimal(diagnostics$minimum_off_diagonal_correlation),
      ", ",
      format_decimal(diagnostics$maximum_off_diagonal_correlation),
      "]"
    ),
    `Minimum eigenvalue` = format_decimal(
      diagnostics$minimum_eigenvalue,
      4L
    ),
    `Negative eigenvalues` = diagnostics$n_negative_eigenvalues,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})
threshold_summary_table <- do.call(rbind, summary_rows)

matrix_difference_rows <- list()
row_index <- 1L
for (first_index in seq_len(length(expected_thresholds) - 1L)) {
  for (second_index in (first_index + 1L):length(expected_thresholds)) {
    first_threshold <- expected_thresholds[first_index]
    second_threshold <- expected_thresholds[second_index]
    difference <- pairwise_analysis$matrices[[first_index]] -
      pairwise_analysis$matrices[[second_index]]
    off_diagonal <- abs(difference[upper.tri(difference)])
    matrix_difference_rows[[row_index]] <- data.frame(
      Comparison = paste0(
        "lfdr > ",
        format_decimal(first_threshold, 2),
        " vs > ",
        format_decimal(second_threshold, 2)
      ),
      `Maximum absolute entry difference` = max(abs(difference)),
      `Mean absolute off-diagonal difference` = mean(off_diagonal),
      `Frobenius difference` = sqrt(sum(difference^2)),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L
  }
}
matrix_difference_table <- do.call(rbind, matrix_difference_rows)

lag1_values <- matrix_diagnostics$mean_lag1_correlation
lag15_values <- matrix_diagnostics$mean_lag15_correlation
off_diagonal_values <- matrix_diagnostics$mean_off_diagonal_correlation
lag1_range <- range(lag1_values)
lag15_range <- range(lag15_values)
off_diagonal_range <- range(off_diagonal_values)
extreme_matrix_difference <- matrix_difference_table[
  matrix_difference_table$Comparison == "lfdr > 0.97 vs > 0.95",
  ,
  drop = FALSE
]

plot_raw_matrix_panel <- function(correlation, threshold, diagnostics, z_limit) {
  palette <- grDevices::colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(201)
  graphics::image(
    x = 0:15,
    y = 0:15,
    z = correlation,
    zlim = c(-z_limit, z_limit),
    col = palette,
    axes = FALSE,
    xlab = "Time",
    ylab = "Time",
    main = paste0(
      "lfdr > ",
      format_decimal(threshold, 2),
      " (n = ",
      diagnostics$n_selected,
      ")\nmin eigenvalue = ",
      format_decimal(diagnostics$minimum_eigenvalue, 4)
    ),
    cex.main = 0.9,
    asp = 1
  )
  graphics::axis(1, at = 0:15, labels = 0:15, cex.axis = 0.66)
  graphics::axis(2, at = 0:15, labels = 0:15, las = 1, cex.axis = 0.66)
  graphics::box()
}

plot_pairwise_raw_matrices <- function() {
  z_limit <- max(abs(unlist(pairwise_analysis$matrices)))
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(1, 3),
    mar = c(4.1, 4.0, 4.5, 0.8),
    oma = c(2.0, 0.4, 3.0, 0.4)
  )
  for (threshold_index in seq_along(expected_thresholds)) {
    threshold <- expected_thresholds[threshold_index]
    diagnostics <- matrix_diagnostics[
      matrix_diagnostics$threshold == threshold,
    ]
    plot_raw_matrix_panel(
      pairwise_analysis$matrices[[threshold_index]],
      threshold,
      diagnostics,
      z_limit
    )
  }
  graphics::mtext(
    "Raw pairwise-difference residual correlation estimates",
    side = 3,
    outer = TRUE,
    line = 1.0,
    cex = 1.35,
    font = 2
  )
  graphics::mtext(
    "Common scale: blue = negative, white = zero, red = positive; no PD projection",
    side = 1,
    outer = TRUE,
    line = 0.3,
    cex = 0.82
  )
}

plot_pairwise_variograms <- function() {
  colors <- c("0.97" = "#0072B2", "0.96" = "#D55E00", "0.95" = "#009E73")
  y_limits <- range(c(
    bootstrap_intervals$semivariogram_ci_lower,
    bootstrap_intervals$semivariogram_ci_upper
  ))
  graphics::plot(
    1:15,
    rep(NA_real_, 15L),
    type = "n",
    xlab = "Time lag",
    ylab = "Standardized semivariogram",
    ylim = y_limits,
    main = "Raw pairwise-difference variogram"
  )
  for (threshold in rev(expected_thresholds)) {
    subset <- bootstrap_intervals[
      bootstrap_intervals$threshold == threshold,
    ]
    color <- colors[format_decimal(threshold, 2)]
    graphics::polygon(
      c(subset$lag, rev(subset$lag)),
      c(
        subset$semivariogram_ci_lower,
        rev(subset$semivariogram_ci_upper)
      ),
      col = grDevices::adjustcolor(color, alpha.f = 0.12),
      border = NA
    )
  }
  for (threshold in expected_thresholds) {
    subset <- bootstrap_intervals[
      bootstrap_intervals$threshold == threshold,
    ]
    color <- colors[format_decimal(threshold, 2)]
    graphics::lines(
      subset$lag,
      subset$observed_semivariogram,
      type = "b",
      pch = 19,
      lwd = 2,
      col = color
    )
  }
  graphics::legend(
    "topright",
    legend = paste0(
      "lfdr > ",
      format_decimal(expected_thresholds, 2),
      " (95% gene bootstrap)"
    ),
    col = colors[format_decimal(expected_thresholds, 2)],
    pch = 19,
    lwd = 2,
    bty = "n"
  )
}

plot_pairwise_threshold_sensitivity <- function() {
  colors <- c("#0072B2", "#D55E00")
  lag1 <- bootstrap_intervals[bootstrap_intervals$lag == 1L, ]
  lag15 <- bootstrap_intervals[bootstrap_intervals$lag == 15L, ]
  x <- expected_thresholds
  y_limits <- range(c(
    lag1$correlation_ci_lower,
    lag1$correlation_ci_upper,
    lag15$correlation_ci_lower,
    lag15$correlation_ci_upper
  ))
  graphics::plot(
    x,
    lag1$observed_correlation,
    type = "b",
    pch = 19,
    lwd = 2,
    col = colors[1],
    xlim = rev(range(x)),
    ylim = y_limits,
    xlab = "FASH lfdr threshold",
    ylab = "Mean correlation",
    main = "Threshold sensitivity"
  )
  graphics::arrows(
    x,
    lag1$correlation_ci_lower,
    x,
    lag1$correlation_ci_upper,
    angle = 90,
    code = 3,
    length = 0.04,
    col = colors[1]
  )
  graphics::lines(
    x,
    lag15$observed_correlation,
    type = "b",
    pch = 17,
    lwd = 2,
    col = colors[2]
  )
  graphics::arrows(
    x,
    lag15$correlation_ci_lower,
    x,
    lag15$correlation_ci_upper,
    angle = 90,
    code = 3,
    length = 0.04,
    col = colors[2]
  )
  graphics::legend(
    "bottomleft",
    legend = c("Mean lag-1 correlation", "Mean lag-15 correlation"),
    col = colors,
    pch = c(19, 17),
    lwd = 2,
    bty = "n"
  )
}
