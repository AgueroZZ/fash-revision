#!/usr/bin/env Rscript

# Plot single-seed empirical FDR curves for the R1 known-truth pilot.

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

if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("scales", quietly = TRUE)) {
  stop("The ggplot2 and scales packages are required.")
}

workflowr_root <- find_workflowr_root()
output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r1_known_truth_genotype_permutation_seed12345_perm20260811_J200"
)
curve_path <- file.path(output_directory, "alpha_curve.csv")
figure_directory <- file.path(output_directory, "figures")
if (!file.exists(curve_path)) {
  stop("The validated fixed-seed alpha curve is missing.")
}
dir.create(figure_directory, showWarnings = FALSE)

curve <- utils::read.csv(curve_path, stringsAsFactors = FALSE)
required_columns <- c(
  "arm", "fit_stage", "alpha", "n_units", "n_true_null",
  "n_true_alternative", "n_discoveries", "false_discoveries",
  "true_positives", "realized_fdp"
)
missing_columns <- setdiff(required_columns, names(curve))
if (length(missing_columns) > 0L || nrow(curve) != 804L ||
    !setequal(curve$arm, c(
      "genuine_null_baseline",
      "shared_genotype_permutation"
    )) || !setequal(curve$fit_stage, c("Raw", "BF")) ||
    any(curve$n_units != 400L) || any(curve$n_true_null != 200L) ||
    any(curve$n_true_alternative != 200L) ||
    any(!is.finite(curve$alpha)) || any(curve$alpha < 0 | curve$alpha > 0.20) ||
    any(!is.finite(curve$realized_fdp)) ||
    any(curve$realized_fdp < 0 | curve$realized_fdp > 1)) {
  stop("The saved known-truth alpha curve is invalid.")
}

arm_labels <- c(
  genuine_null_baseline = "Genuine simulated nulls",
  shared_genotype_permutation = "Genotype-permuted copies"
)
fit_labels <- c(
  Raw = "Raw empirical-Bayes fit",
  BF = "BF-adjusted fit"
)
curve$null_construction <- factor(
  arm_labels[curve$arm],
  levels = unname(arm_labels)
)
curve$fit <- factor(
  fit_labels[curve$fit_stage],
  levels = unname(fit_labels)
)
alpha005 <- curve[abs(curve$alpha - 0.05) < 1e-12, , drop = FALSE]
if (nrow(alpha005) != 4L) {
  stop("The alpha-0.05 calibration rows are incomplete.")
}
alpha005$point_label <- scales::percent(
  alpha005$realized_fdp,
  accuracy = 0.1
)

arm_colors <- c(
  `Genuine simulated nulls` = "#0072B2",
  `Genotype-permuted copies` = "#D55E00"
)
arm_linetypes <- c(
  `Genuine simulated nulls` = "solid",
  `Genotype-permuted copies` = "longdash"
)

figure <- ggplot2::ggplot(
  curve,
  ggplot2::aes(
    x = alpha,
    y = realized_fdp,
    color = null_construction,
    linetype = null_construction
  )
) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    color = "#555555",
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  ggplot2::geom_line(linewidth = 1.05) +
  ggplot2::geom_point(
    data = alpha005,
    size = 2.6,
    stroke = 0.4
  ) +
  ggplot2::geom_text(
    data = alpha005,
    ggplot2::aes(label = point_label),
    hjust = -0.20,
    vjust = c(1.35, -0.55, 1.35, -0.55),
    size = 3.25,
    show.legend = FALSE
  ) +
  ggplot2::facet_wrap(~fit, nrow = 1L) +
  ggplot2::scale_color_manual(
    values = arm_colors,
    name = "Null construction"
  ) +
  ggplot2::scale_linetype_manual(
    values = arm_linetypes,
    name = "Null construction"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 0.20),
    breaks = seq(0, 0.20, by = 0.05),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 0.32),
    breaks = seq(0, 0.30, by = 0.05),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0.01, 0.03))
  ) +
  ggplot2::labs(
    title = "Genotype-permuted nulls show higher realized FDP across alpha",
    subtitle = paste(
      "R1 seed 12345; each arm has 200 known alternatives and 200 known nulls.",
      "Labels mark alpha = 5%."
    ),
    x = expression(paste("Nominal cumulative-FDR level, ", alpha)),
    y = "Single-seed empirical FDR (realized FDP)",
    caption = paste(
      "The dotted diagonal is empirical FDR = nominal alpha.",
      "One simulation and permutation seed diagnose this realization, not expected FDR."
    )
  ) +
  ggplot2::theme_classic(base_size = 12) +
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
      margin = ggplot2::margin(b = 10)
    ),
    plot.caption = ggplot2::element_text(
      size = 9,
      color = "#5A5A5A",
      hjust = 0,
      margin = ggplot2::margin(t = 10)
    ),
    strip.background = ggplot2::element_rect(
      fill = "#F3F3F3",
      color = NA
    ),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 11.5,
      color = "#222222",
      margin = ggplot2::margin(7, 7, 7, 7)
    ),
    axis.title = ggplot2::element_text(size = 11.5, color = "#222222"),
    axis.text = ggplot2::element_text(size = 10, color = "#333333"),
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
    panel.grid.major.y = ggplot2::element_line(
      color = "#E6E6E6",
      linewidth = 0.4
    ),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(1.2, "lines"),
    legend.position = "top",
    legend.title = ggplot2::element_text(face = "bold", size = 10.5),
    legend.text = ggplot2::element_text(size = 10),
    plot.margin = ggplot2::margin(14, 18, 12, 14)
  )

ggplot2::ggsave(
  file.path(figure_directory, "known_truth_empirical_fdr_vs_alpha.png"),
  figure,
  width = 11.5,
  height = 6.6,
  dpi = 320,
  bg = "white"
)
ggplot2::ggsave(
  file.path(figure_directory, "known_truth_empirical_fdr_vs_alpha.pdf"),
  figure,
  width = 11.5,
  height = 6.6,
  device = grDevices::cairo_pdf,
  bg = "white"
)

cat("R1 known-truth empirical-FDR calibration figure created.\n")
