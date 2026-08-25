#!/usr/bin/env Rscript

# Compare formal R1, genuine nulls, naive genotype permutations, and
# signal-stripped donor-residual permutations across five source seeds.

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
control_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r1_known_truth_genotype_permutation_mc5_J200_v1"
)
residual_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r1_signal_stripped_unadjusted_residual_permutation_mc5_J200_v1"
)
formal_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "mc",
  "r1_random_bspline_main_effect_profile_sigma_pilot5",
  "summary", "iwp_vs_linear_fash_mc_alpha_curve.csv"
)
control_path <- file.path(
  control_directory, "summary", "mc_alpha_curve.csv"
)
residual_path <- file.path(
  residual_directory, "summary", "mc_alpha_curve.csv"
)
figure_directory <- file.path(residual_directory, "figures")
if (any(!file.exists(c(formal_path, control_path, residual_path)))) {
  stop("At least one five-seed calibration summary is missing.")
}
dir.create(figure_directory, showWarnings = FALSE)

formal <- utils::read.csv(formal_path, stringsAsFactors = FALSE)
control <- utils::read.csv(control_path, stringsAsFactors = FALSE)
residual <- utils::read.csv(residual_path, stringsAsFactors = FALSE)
required_columns <- c(
  "alpha", "n_replications", "mean_fdr", "fdr_ci_lower", "fdr_ci_upper"
)
if (length(setdiff(c(required_columns, "method"), names(formal))) > 0L ||
    length(setdiff(c(required_columns, "arm", "fit_stage"), names(control))) > 0L ||
    length(setdiff(c(required_columns, "arm", "fit_stage"), names(residual))) > 0L) {
  stop("A calibration table is missing required columns.")
}
formal <- formal[
  formal$method %in% c("FASH-IWP1-Raw", "FASH-IWP1-BF") &
    formal$alpha <= 0.20 + 1e-12,
  , drop = FALSE
]
control <- control[control$alpha <= 0.20 + 1e-12, , drop = FALSE]
residual <- residual[residual$alpha <= 0.20 + 1e-12, , drop = FALSE]
if (nrow(formal) != 82L || nrow(control) != 804L || nrow(residual) != 402L ||
    any(formal$n_replications != 5L) ||
    any(control$n_replications != 5L) ||
    any(residual$n_replications != 5L)) {
  stop("The five-seed calibration tables are incomplete.")
}

series_labels <- c(
  formal_r1 = "Formal R1: all 1,000 units",
  genuine_null_baseline = "Genuine simulated nulls",
  shared_genotype_permutation = "Naive genotype-permuted copies",
  signal_stripped_residual_permutation =
    "Signal-stripped residual permutations"
)
fit_labels <- c(
  Raw = "Raw empirical-Bayes fit",
  BF = "BF-adjusted fit"
)
formal$series_key <- "formal_r1"
formal$fit_key <- ifelse(
  formal$method == "FASH-IWP1-Raw", "Raw", "BF"
)
control$series_key <- control$arm
control$fit_key <- control$fit_stage
residual$series_key <- residual$arm
residual$fit_key <- residual$fit_stage
curve <- rbind(
  formal[, c(
    "series_key", "fit_key", "alpha", "n_replications", "mean_fdr",
    "fdr_ci_lower", "fdr_ci_upper"
  )],
  control[, c(
    "series_key", "fit_key", "alpha", "n_replications", "mean_fdr",
    "fdr_ci_lower", "fdr_ci_upper"
  )],
  residual[, c(
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
if (nrow(alpha005) != 8L) {
  stop("The alpha-0.05 values are incomplete.")
}
annotation_rows <- lapply(levels(curve$fit), function(fit_name) {
  rows <- alpha005[alpha005$fit == fit_name, , drop = FALSE]
  rows <- rows[match(unname(series_labels), as.character(rows$series)), ]
  short_names <- c("Formal", "Genuine", "Naive", "Signal-stripped")
  data.frame(
    fit = factor(fit_name, levels = levels(curve$fit)),
    x = 0.006,
    y = 0.358,
    label = paste0(
      "alpha = 5% mean empirical FDR\n",
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
  `Genuine simulated nulls` = "#0072B2",
  `Naive genotype-permuted copies` = "#D55E00",
  `Signal-stripped residual permutations` = "#009E73"
)
series_linetypes <- c(
  `Formal R1: all 1,000 units` = "dotdash",
  `Genuine simulated nulls` = "solid",
  `Naive genotype-permuted copies` = "longdash",
  `Signal-stripped residual permutations` = "twodash"
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
    alpha = 0.075,
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
    size = 2.4,
    stroke = 0.4,
    show.legend = FALSE
  ) +
  ggplot2::geom_label(
    data = annotations,
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 2.95,
    lineheight = 1.15,
    linewidth = 0.2,
    label.padding = grid::unit(0.14, "lines"),
    color = "#333333",
    fill = "white"
  ) +
  ggplot2::facet_wrap(~fit, nrow = 1L) +
  ggplot2::scale_color_manual(values = series_colors, name = "Analysis") +
  ggplot2::scale_fill_manual(values = series_colors, guide = "none") +
  ggplot2::scale_linetype_manual(
    values = series_linetypes, name = "Analysis"
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
    title = "Removing the fitted genotype effect restores R1 null calibration",
    subtitle = paste(
      "The signal-stripped residual-permutation curve tracks genuine simulated",
      "nulls; naive genotype permutation retains a large leakage bias."
    ),
    x = expression(paste("Nominal cumulative-FDR level, ", alpha)),
    y = "Mean empirical FDR across five seeds",
    caption = paste(
      "Each 400-unit arm contains the same 200 known alternatives plus 200 nulls.",
      "Signal-stripped nulls use unadjusted full-model residuals because four",
      "seeds contain one exact-leverage-one unit, for which HC2 is undefined.\n",
      "Shading gives 95% t-based Monte Carlo intervals."
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
      size = 8.6, color = "#5A5A5A", hjust = 0,
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
    legend.text = ggplot2::element_text(size = 9.2),
    legend.key.width = grid::unit(1.8, "lines"),
    plot.margin = ggplot2::margin(14, 18, 12, 14)
  )

png_path <- file.path(
  figure_directory, "signal_stripped_residual_permutation_mc_fdr.png"
)
pdf_path <- file.path(
  figure_directory, "signal_stripped_residual_permutation_mc_fdr.pdf"
)
ggplot2::ggsave(
  png_path,
  figure,
  width = 13.0,
  height = 7.1,
  dpi = 320,
  bg = "white"
)
ggplot2::ggsave(
  pdf_path,
  figure,
  width = 13.0,
  height = 7.1,
  device = grDevices::cairo_pdf,
  bg = "white"
)

cat("Signal-stripped residual-permutation calibration figure created.\n")
