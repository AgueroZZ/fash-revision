#!/usr/bin/env Rscript

# Plot paired R3 Middle calibration curves across immutable designs.

options(stringsAsFactors = FALSE)

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[[1L]] == length(args)) return(default)
  args[[hit[[1L]] + 1L]]
}

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }
  if (file.exists(file.path(
    "coderepo-local", "code", "revision_simulations", "shared",
    "simulation_functions.R"
  ))) {
    return(normalizePath(
      "coderepo-local", winslash = "/", mustWork = TRUE
    ))
  }
  stop("Could not find the workflowr repository root.", call. = FALSE)
}

read_design <- function(summary_path, design) {
  if (!file.exists(summary_path)) {
    stop("Missing summary: ", summary_path, call. = FALSE)
  }
  data <- utils::read.csv(summary_path)
  data <- data[
    data$target == "middle" & data$method == "FASH-IWP1-BF",
    , drop = FALSE
  ]
  data$truth_mechanism <- ifelse(
    grepl("^r3a_", data$scenario),
    "R3A: broad random B-spline",
    "R3B: compact raised cosine"
  )
  data$design <- design
  data
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The plot requires the existing ggplot2 installation.", call. = FALSE)
}

workflowr_root <- find_workflowr_root()
full_result_argument <- get_arg("--full-result-dir", "")
if (!nzchar(full_result_argument)) {
  stop("--full-result-dir is required.", call. = FALSE)
}
full_result_dir <- normalizePath(
  full_result_argument, winslash = "/", mustWork = TRUE
)

cache_root <- file.path(
  workflowr_root, "output", "revision_simulations", "mc"
)
design_dirs <- c(
  "Historical closed Middle" = file.path(
    cache_root,
    paste0(
      "r3_real_genotype_one_per_gene_J6362_matched_functional_",
      "relative_clearance_main_effect_fashr0143_pilot5"
    )
  ),
  "Open Middle, balanced boundary design" = file.path(
    cache_root,
    paste0(
      "r3_real_genotype_one_per_gene_J6362_matched_functional_",
      "open_middle_3_12_support_contained_relative_clearance_",
      "main_effect_fashr0143_pilot5"
    )
  ),
  "Full-support IWP1-mixture candidate" = full_result_dir
)

plot_rows <- Map(
  function(result_dir, design) {
    read_design(
      file.path(
        result_dir, "summary", "functional_testing_mc_alpha_curve.csv"
      ),
      design
    )
  },
  result_dir = unname(design_dirs),
  design = names(design_dirs)
)
plot_data <- do.call(rbind, plot_rows)
rownames(plot_data) <- NULL
plot_data$design <- factor(plot_data$design, levels = names(design_dirs))
plot_data$truth_mechanism <- factor(
  plot_data$truth_mechanism,
  levels = c(
    "R3A: broad random B-spline",
    "R3B: compact raised cosine"
  )
)

output_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "diagnostics",
  "r3_middle_calibration"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
source_columns <- c(
  "design", "truth_mechanism", "alpha", "n_replications",
  "mean_discoveries", "mean_estimated_fsr", "mean_empirical_fsr",
  "empirical_fsr_mc_se", "empirical_fsr_ci_lower",
  "empirical_fsr_ci_upper"
)
utils::write.csv(
  plot_data[, source_columns, drop = FALSE],
  file.path(output_dir, "r3_middle_calibration_design_plot_data.csv"),
  row.names = FALSE
)

candidate <- plot_data[
  plot_data$design == "Full-support IWP1-mixture candidate",
  , drop = FALSE
]
r3a_annotation <- candidate[
  candidate$truth_mechanism == "R3A: broad random B-spline" &
    abs(candidate$alpha - 0.20) < 1e-12,
  , drop = FALSE
]
r3b_candidate <- candidate[
  candidate$truth_mechanism == "R3B: compact raised cosine",
  , drop = FALSE
]
r3b_maximum_index <- which.max(
  r3b_candidate$mean_empirical_fsr - r3b_candidate$alpha
)
r3b_annotation <- r3b_candidate[r3b_maximum_index, , drop = FALSE]

annotation_points <- rbind(r3a_annotation, r3b_annotation)
annotation_points$label <- c(
  "0.249 at alpha = 0.20\n(excess = 0.049)",
  "Maximum excess = 0.035\nat alpha = 0.075"
)
annotation_points$label_x <- c(0.115, 0.128)
annotation_points$label_y <- c(0.292, 0.175)

colors <- c(
  "Historical closed Middle" = "#777777",
  "Open Middle, balanced boundary design" = "#D55E00",
  "Full-support IWP1-mixture candidate" = "#0072B2"
)
line_types <- c(
  "Historical closed Middle" = "longdash",
  "Open Middle, balanced boundary design" = "dotted",
  "Full-support IWP1-mixture candidate" = "solid"
)

ggplot2::theme_set(ggplot2::theme_minimal(base_size = 11))
plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = alpha,
    y = mean_empirical_fsr,
    color = design,
    fill = design,
    linetype = design,
    group = design
  )
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = empirical_fsr_ci_lower,
      ymax = empirical_fsr_ci_upper
    ),
    alpha = 0.11,
    linewidth = 0,
    show.legend = FALSE
  ) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    color = "#222222",
    linetype = "dotdash",
    linewidth = 0.45,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_point(
    data = annotation_points,
    ggplot2::aes(x = alpha, y = mean_empirical_fsr),
    color = "#0072B2",
    size = 2.2,
    inherit.aes = FALSE
  ) +
  ggplot2::geom_segment(
    data = annotation_points,
    ggplot2::aes(
      x = label_x,
      y = label_y,
      xend = alpha,
      yend = mean_empirical_fsr
    ),
    color = "#0072B2",
    linewidth = 0.45,
    arrow = grid::arrow(length = grid::unit(0.08, "inches")),
    inherit.aes = FALSE
  ) +
  ggplot2::geom_text(
    data = annotation_points,
    ggplot2::aes(x = label_x, y = label_y, label = label),
    color = "#005A8D",
    size = 3.2,
    hjust = 0,
    lineheight = 1.05,
    inherit.aes = FALSE
  ) +
  ggplot2::facet_wrap(~truth_mechanism, nrow = 1) +
  ggplot2::scale_color_manual(values = colors, name = NULL) +
  ggplot2::scale_fill_manual(values = colors, name = NULL) +
  ggplot2::scale_linetype_manual(values = line_types, name = NULL) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 0.20, by = 0.05),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(0, 0.30, by = 0.05),
    expand = ggplot2::expansion(mult = c(0, 0.015))
  ) +
  ggplot2::coord_cartesian(
    xlim = c(0, 0.20),
    ylim = c(0, 0.33),
    expand = FALSE
  ) +
  ggplot2::labs(
    title = paste(
      "The full-support mixture substantially reduces, but does not",
      "eliminate, Middle FSR inflation"
    ),
    subtitle = paste(
      "Curves show mean empirical FSR; shaded bands are 95% t intervals",
      "across five Monte Carlo replicates"
    ),
    x = "Nominal FSR threshold (alpha)",
    y = "Mean empirical FSR",
    caption = paste(
      "The candidate misses the prespecified maximum-excess <= 0.03 gate;",
      "the magnitude and uncertainty should be interpreted separately."
    )
  ) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(
      color = "#D9D9D9", linewidth = 0.35
    ),
    strip.text = ggplot2::element_text(face = "bold", size = 11.5),
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(color = "#444444", size = 10.5),
    plot.caption = ggplot2::element_text(
      color = "#444444", hjust = 0, size = 9
    ),
    legend.position = "top",
    legend.box = "vertical",
    legend.text = ggplot2::element_text(size = 9.5),
    axis.title = ggplot2::element_text(size = 10.5),
    panel.spacing = grid::unit(1.2, "lines")
  ) +
  ggplot2::guides(
    color = ggplot2::guide_legend(nrow = 2, byrow = TRUE),
    linetype = ggplot2::guide_legend(nrow = 2, byrow = TRUE)
  )

png_path <- file.path(
  output_dir, "r3_middle_calibration_design_comparison.png"
)
pdf_path <- file.path(
  output_dir, "r3_middle_calibration_design_comparison.pdf"
)
ggplot2::ggsave(
  png_path, plot = plot, width = 12.4, height = 5.8, dpi = 220,
  bg = "white"
)
ggplot2::ggsave(
  pdf_path, plot = plot, width = 12.4, height = 5.8,
  device = grDevices::cairo_pdf, bg = "white"
)

cat(png_path, "\n")
cat(pdf_path, "\n")
