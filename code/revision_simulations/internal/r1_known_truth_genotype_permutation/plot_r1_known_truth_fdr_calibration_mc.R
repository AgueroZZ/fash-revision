#!/usr/bin/env Rscript

# Plot five-seed Monte Carlo FDR curves for the R1 known-truth permutation
# diagnostic, with the formal R1 calibration curve as a reference.

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
experiment_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r1_known_truth_genotype_permutation_mc5_J200_v1"
)
experiment_curve_path <- file.path(
  experiment_directory, "summary", "mc_alpha_curve.csv"
)
formal_curve_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "mc",
  "r1_random_bspline_main_effect_profile_sigma_pilot5",
  "summary", "iwp_vs_linear_fash_mc_alpha_curve.csv"
)
figure_directory <- file.path(experiment_directory, "figures")
if (!file.exists(experiment_curve_path) || !file.exists(formal_curve_path)) {
  stop("A five-seed experiment or formal R1 calibration summary is missing.")
}
dir.create(figure_directory, showWarnings = FALSE)

experiment <- utils::read.csv(
  experiment_curve_path, stringsAsFactors = FALSE
)
formal <- utils::read.csv(formal_curve_path, stringsAsFactors = FALSE)
required_curve_columns <- c(
  "alpha", "n_replications", "mean_fdr", "fdr_ci_lower", "fdr_ci_upper"
)
if (length(setdiff(
      c(required_curve_columns, "arm", "fit_stage"), names(experiment)
    )) > 0L ||
    length(setdiff(
      c(required_curve_columns, "method"), names(formal)
    )) > 0L) {
  stop("A calibration summary is missing required columns.")
}

experiment <- experiment[
  experiment$alpha <= 0.20 + 1e-12, , drop = FALSE
]
formal <- formal[
  formal$method %in% c("FASH-IWP1-Raw", "FASH-IWP1-BF") &
    formal$alpha <= 0.20 + 1e-12,
  , drop = FALSE
]
if (nrow(experiment) != 804L || nrow(formal) != 82L ||
    any(experiment$n_replications != 5L) ||
    any(formal$n_replications != 5L)) {
  stop("The five-seed calibration curves are incomplete.")
}

series_labels <- c(
  formal_r1 = "Formal R1: all 1,000 units",
  genuine_null_baseline = "200 alternatives + 200 genuine nulls",
  shared_genotype_permutation =
    "200 alternatives + 200 genotype-permuted copies"
)
fit_labels <- c(
  Raw = "Raw empirical-Bayes fit",
  BF = "BF-adjusted fit"
)
experiment$series_key <- experiment$arm
experiment$fit_key <- experiment$fit_stage
formal$series_key <- "formal_r1"
formal$fit_key <- ifelse(
  formal$method == "FASH-IWP1-Raw", "Raw", "BF"
)
curve <- rbind(
  experiment[, c(
    "series_key", "fit_key", "alpha", "n_replications", "mean_fdr",
    "fdr_ci_lower", "fdr_ci_upper"
  )],
  formal[, c(
    "series_key", "fit_key", "alpha", "n_replications", "mean_fdr",
    "fdr_ci_lower", "fdr_ci_upper"
  )]
)
curve$series <- factor(
  series_labels[curve$series_key], levels = unname(series_labels)
)
curve$fit <- factor(
  fit_labels[curve$fit_key], levels = unname(fit_labels)
)
curve$fdr_ci_lower <- pmax(0, curve$fdr_ci_lower)
curve$fdr_ci_upper <- pmin(1, curve$fdr_ci_upper)

alpha005 <- curve[abs(curve$alpha - 0.05) < 1e-12, , drop = FALSE]
if (nrow(alpha005) != 6L) {
  stop("The alpha-0.05 calibration values are incomplete.")
}
annotation_rows <- lapply(levels(curve$fit), function(fit_name) {
  rows <- alpha005[alpha005$fit == fit_name, , drop = FALSE]
  rows <- rows[match(unname(series_labels), as.character(rows$series)), ]
  short_names <- c("Formal R1", "Genuine null", "Permuted copy")
  data.frame(
    fit = factor(fit_name, levels = levels(curve$fit)),
    x = 0.006,
    y = 0.358,
    label = paste0(
      expression(alpha),
      " = 5% mean empirical FDR\n",
      paste0(
        short_names, ": ",
        scales::percent(rows$mean_fdr, accuracy = 0.1),
        collapse = "   |   "
      )
    ),
    stringsAsFactors = FALSE
  )
})
annotations <- do.call(rbind, annotation_rows)

series_colors <- c(
  `Formal R1: all 1,000 units` = "#4D4D4D",
  `200 alternatives + 200 genuine nulls` = "#0072B2",
  `200 alternatives + 200 genotype-permuted copies` = "#D55E00"
)
series_linetypes <- c(
  `Formal R1: all 1,000 units` = "dotdash",
  `200 alternatives + 200 genuine nulls` = "solid",
  `200 alternatives + 200 genotype-permuted copies` = "longdash"
)

figure <- ggplot2::ggplot(
  curve,
  ggplot2::aes(
    x = alpha,
    y = mean_fdr,
    color = series,
    fill = series,
    linetype = series,
    group = series
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(ymin = fdr_ci_lower, ymax = fdr_ci_upper),
    alpha = 0.09,
    color = NA,
    show.legend = FALSE
  ) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    color = "#777777",
    linewidth = 0.75,
    linetype = "dotted"
  ) +
  ggplot2::geom_vline(
    xintercept = 0.05,
    color = "#AAAAAA",
    linewidth = 0.45,
    linetype = "dashed"
  ) +
  ggplot2::geom_line(linewidth = 1.05) +
  ggplot2::geom_point(
    data = alpha005,
    size = 2.5,
    stroke = 0.4,
    show.legend = FALSE
  ) +
  ggplot2::geom_label(
    data = annotations,
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3.05,
    lineheight = 1.15,
    linewidth = 0.2,
    label.padding = grid::unit(0.14, "lines"),
    color = "#333333",
    fill = "white"
  ) +
  ggplot2::facet_wrap(~fit, nrow = 1L) +
  ggplot2::scale_color_manual(
    values = series_colors,
    name = "Analysis"
  ) +
  ggplot2::scale_fill_manual(values = series_colors, guide = "none") +
  ggplot2::scale_linetype_manual(
    values = series_linetypes,
    name = "Analysis"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 0.20),
    breaks = seq(0, 0.20, by = 0.05),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 0.37),
    breaks = seq(0, 0.35, by = 0.05),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::labs(
    title = "Five-seed R1 diagnostic isolates permutation-induced miscalibration",
    subtitle = paste(
      "Genuine simulated nulls recover the original R1 calibration pattern;",
      "genotype-permuted copies of known alternatives do not."
    ),
    x = expression(paste("Nominal cumulative-FDR level, ", alpha)),
    y = "Mean empirical FDR across five seeds",
    caption = paste(
      "Curves are means of realized FDP over the five formal R1 source seeds;",
      "shading gives 95% t-based Monte Carlo intervals.",
      "The dotted diagonal is empirical FDR = nominal alpha."
    )
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold", size = 14, color = "#111111",
      margin = ggplot2::margin(b = 5)
    ),
    plot.subtitle = ggplot2::element_text(
      size = 10.5, color = "#4D4D4D",
      margin = ggplot2::margin(b = 10)
    ),
    plot.caption = ggplot2::element_text(
      size = 8.8, color = "#5A5A5A", hjust = 0,
      margin = ggplot2::margin(t = 10)
    ),
    strip.background = ggplot2::element_rect(fill = "#F3F3F3", color = NA),
    strip.text = ggplot2::element_text(
      face = "bold", size = 11.5, color = "#222222",
      margin = ggplot2::margin(7, 7, 7, 7)
    ),
    axis.title = ggplot2::element_text(size = 11.5, color = "#222222"),
    axis.text = ggplot2::element_text(size = 10, color = "#333333"),
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
    panel.grid.major.y = ggplot2::element_line(
      color = "#E6E6E6", linewidth = 0.4
    ),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(1.2, "lines"),
    legend.position = "top",
    legend.title = ggplot2::element_text(face = "bold", size = 10.5),
    legend.text = ggplot2::element_text(size = 9.5),
    legend.key.width = grid::unit(2.1, "lines"),
    plot.margin = ggplot2::margin(14, 18, 12, 14)
  )

png_path <- file.path(
  figure_directory, "known_truth_mc_fdr_vs_alpha_with_formal_r1.png"
)
pdf_path <- file.path(
  figure_directory, "known_truth_mc_fdr_vs_alpha_with_formal_r1.pdf"
)
ggplot2::ggsave(
  png_path,
  figure,
  width = 12.5,
  height = 7.0,
  dpi = 320,
  bg = "white"
)
ggplot2::ggsave(
  pdf_path,
  figure,
  width = 12.5,
  height = 7.0,
  device = grDevices::cairo_pdf,
  bg = "white"
)

cat("R1 known-truth five-seed MC calibration figure created.\n")
