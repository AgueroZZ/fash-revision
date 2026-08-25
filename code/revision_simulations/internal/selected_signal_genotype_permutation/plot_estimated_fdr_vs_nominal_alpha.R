#!/usr/bin/env Rscript

# Plot permutation-estimated target FDR across nominal alpha thresholds.

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

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required.")
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("The scales package is required.")
}

cumulative_lfdr_calls <- function(lfdr, alpha) {
  ordering <- order(lfdr, method = "radix")
  cumulative_mean <- cumsum(lfdr[ordering]) / seq_along(ordering)
  accepted <- which(cumulative_mean <= alpha)
  if (length(accepted) == 0L) {
    return(integer())
  }
  ordering[accepted]
}

summarize_curve <- function(directory, selection_method, alpha_grid) {
  lfdr_path <- file.path(directory, "unit_lfdr.csv")
  calibration_path <- file.path(directory, "calibration_diagnostics.csv")
  if (!file.exists(lfdr_path) || !file.exists(calibration_path)) {
    stop("A required completed-pilot artifact is missing: ", directory)
  }
  unit_lfdr <- utils::read.csv(lfdr_path, stringsAsFactors = FALSE)
  calibration <- utils::read.csv(
    calibration_path,
    stringsAsFactors = FALSE
  )
  stage_data <- unit_lfdr[
    unit_lfdr$fit_stage == "BF-adjusted",
    ,
    drop = FALSE
  ]
  calibration_row <- calibration[
    calibration$fit_stage == "BF-adjusted",
    ,
    drop = FALSE
  ]
  if (nrow(stage_data) != 2354L || nrow(calibration_row) != 1L ||
      !setequal(unique(stage_data$group), c("target", "permuted_null")) ||
      any(!is.finite(stage_data$lfdr)) ||
      any(stage_data$lfdr < 0 | stage_data$lfdr > 1)) {
    stop("The saved BF-adjusted lfdr data are invalid: ", selection_method)
  }
  target <- stage_data$group == "target"
  permuted_null <- stage_data$group == "permuted_null"
  n_target <- sum(target)
  n_null <- sum(permuted_null)
  pi0_target <- calibration_row$pi0_target_unbounded
  pi0_target_valid <- isTRUE(calibration_row$pi0_target_valid)

  rows <- lapply(alpha_grid, function(alpha) {
    selected_indices <- cumulative_lfdr_calls(stage_data$lfdr, alpha)
    selected <- seq_len(nrow(stage_data)) %in% selected_indices
    target_calls <- sum(selected & target)
    null_calls <- sum(selected & permuted_null)
    total_calls <- target_calls + null_calls
    permutation_fdr <- if (target_calls > 0L) {
      (null_calls / n_null) / (target_calls / n_target)
    } else {
      NA_real_
    }
    data.frame(
      target_selection_method = selection_method,
      nominal_alpha = alpha,
      n_target = n_target,
      n_permuted_null = n_null,
      target_calls = target_calls,
      permuted_null_calls = null_calls,
      total_calls = total_calls,
      target_call_rate = target_calls / n_target,
      permuted_null_call_rate = null_calls / n_null,
      permutation_estimated_target_fdr = permutation_fdr,
      known_null_fraction_merged = if (total_calls > 0L) {
        null_calls / total_calls
      } else {
        NA_real_
      },
      fitted_pi0_target = pi0_target,
      fitted_pi0_target_valid = pi0_target_valid,
      fitted_pi0_scaled_target_fdr = if (pi0_target_valid) {
        pi0_target * permutation_fdr
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

workflowr_root <- find_workflowr_root()
output_parent <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal"
)
pilot_directories <- c(
  "Top variant" = file.path(
    output_parent,
    "selected_signal_residual_block_permutation_seed20260811"
  ),
  "Random variant" = file.path(
    output_parent,
    paste0(
      "selected_signal_random_all_tested_residual_block_permutation_",
      "selection20260817_seed20260811"
    )
  )
)
random_directory <- unname(pilot_directories[["Random variant"]])
figure_directory <- file.path(random_directory, "figures")
dir.create(figure_directory, showWarnings = FALSE)

alpha_grid <- seq(0.001, 0.20, by = 0.001)
curve_data <- do.call(rbind, lapply(names(pilot_directories), function(method) {
  summarize_curve(
    directory = pilot_directories[[method]],
    selection_method = method,
    alpha_grid = alpha_grid
  )
}))
rownames(curve_data) <- NULL
curve_data$target_selection_method <- factor(
  curve_data$target_selection_method,
  levels = c("Top variant", "Random variant")
)

alpha_005 <- curve_data[
  abs(curve_data$nominal_alpha - 0.05) < 1e-12,
  ,
  drop = FALSE
]
top_005 <- alpha_005[
  alpha_005$target_selection_method == "Top variant",
  ,
  drop = FALSE
]
random_005 <- alpha_005[
  alpha_005$target_selection_method == "Random variant",
  ,
  drop = FALSE
]
if (nrow(alpha_005) != 2L || top_005$target_calls != 1177L ||
    top_005$permuted_null_calls != 385L ||
    random_005$target_calls != 123L ||
    random_005$permuted_null_calls != 24L ||
    abs(top_005$permutation_estimated_target_fdr - 385 / 1177) > 1e-14 ||
    abs(random_005$permutation_estimated_target_fdr - 24 / 123) > 1e-14) {
  stop("The alpha-0.05 curve values do not reproduce the validated results.")
}

curve_path <- file.path(
  random_directory,
  "estimated_target_fdr_vs_nominal_alpha.csv"
)
utils::write.csv(curve_data, curve_path, row.names = FALSE)

plot_data <- curve_data[
  is.finite(curve_data$permutation_estimated_target_fdr),
  ,
  drop = FALSE
]
curve_maximum <- max(plot_data$permutation_estimated_target_fdr)
y_upper <- max(0.40, 1.08 * curve_maximum)
method_colors <- c(
  "Top variant" = "#D55E00",
  "Random variant" = "#0072B2"
)

figure <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = nominal_alpha,
    y = permutation_estimated_target_fdr,
    color = target_selection_method
  )
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    color = "#6B6B6B",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  ggplot2::geom_step(linewidth = 1.15, direction = "hv") +
  ggplot2::geom_vline(
    xintercept = 0.05,
    color = "#9A9A9A",
    linewidth = 0.55,
    linetype = "dotted"
  ) +
  ggplot2::geom_point(
    data = alpha_005,
    size = 2.8,
    fill = "white",
    shape = 21,
    stroke = 1
  ) +
  ggplot2::annotate(
    "label",
    x = 0.054,
    y = top_005$permutation_estimated_target_fdr,
    label = sprintf(
      "Top: %.1f%% (%d/%d)",
      100 * top_005$permutation_estimated_target_fdr,
      top_005$permuted_null_calls,
      top_005$target_calls
    ),
    hjust = 0,
    vjust = -0.25,
    size = 3.4,
    color = method_colors[["Top variant"]],
    fill = "white",
    linewidth = 0.2
  ) +
  ggplot2::annotate(
    "label",
    x = 0.054,
    y = random_005$permutation_estimated_target_fdr,
    label = sprintf(
      "Random: %.1f%% (%d/%d)",
      100 * random_005$permutation_estimated_target_fdr,
      random_005$permuted_null_calls,
      random_005$target_calls
    ),
    hjust = 0,
    vjust = 1.25,
    size = 3.4,
    color = method_colors[["Random variant"]],
    fill = "white",
    linewidth = 0.2
  ) +
  ggplot2::annotate(
    "text",
    x = 0.175,
    y = 0.175,
    label = "Nominal reference",
    angle = 31,
    hjust = 0.5,
    vjust = -0.6,
    size = 3.2,
    color = "#5F5F5F"
  ) +
  ggplot2::scale_color_manual(values = method_colors, name = NULL) +
  ggplot2::scale_x_continuous(
    limits = c(0, 0.20),
    breaks = seq(0, 0.20, by = 0.025),
    labels = scales::label_percent(accuracy = 0.1),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, y_upper),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    title = "Random variant selection lowers permutation-estimated target FDR",
    subtitle = paste(
      "BF-adjusted cumulative-lfdr calls across nominal thresholds;",
      "the dashed line marks estimated FDR = nominal alpha"
    ),
    x = "Nominal alpha",
    y = "Permutation-estimated target FDR",
    caption = paste(
      paste(
        "Estimator: (V / M) / (R[target] / N[target]), equivalent to",
        "a conservative target pi0 of 1."
      ),
      paste(
        "One shared residual-block donor map; curves are conditional",
        "diagnostics, not formal repeated-sampling FDR."
      ),
      sep = "\n"
    )
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 14,
      color = "#111111",
      margin = ggplot2::margin(b = 5)
    ),
    plot.subtitle = ggplot2::element_text(
      size = 10.5,
      color = "#4D4D4D",
      margin = ggplot2::margin(b = 11)
    ),
    plot.caption = ggplot2::element_text(
      size = 8.8,
      color = "#5A5A5A",
      hjust = 0,
      lineheight = 1.05,
      margin = ggplot2::margin(t = 10)
    ),
    axis.title = ggplot2::element_text(size = 11.5, color = "#222222"),
    axis.text = ggplot2::element_text(size = 10, color = "#333333"),
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(
      color = "#E4E4E4",
      linewidth = 0.4
    ),
    legend.position = "top",
    legend.justification = "left",
    legend.text = ggplot2::element_text(size = 10.5),
    plot.margin = ggplot2::margin(14, 18, 12, 14)
  )

png_path <- file.path(
  figure_directory,
  "estimated_target_fdr_vs_nominal_alpha.png"
)
pdf_path <- file.path(
  figure_directory,
  "estimated_target_fdr_vs_nominal_alpha.pdf"
)
ggplot2::ggsave(
  png_path,
  figure,
  width = 9.5,
  height = 6.4,
  dpi = 320,
  bg = "white"
)
ggplot2::ggsave(
  pdf_path,
  figure,
  width = 9.5,
  height = 6.4,
  device = grDevices::cairo_pdf,
  bg = "white"
)

cat("Estimated-FDR curve plot created.\n")
cat("Curve data: ", curve_path, "\n", sep = "")
cat("PNG: ", png_path, "\n", sep = "")
cat("PDF: ", pdf_path, "\n", sep = "")
print(alpha_005, row.names = FALSE)
