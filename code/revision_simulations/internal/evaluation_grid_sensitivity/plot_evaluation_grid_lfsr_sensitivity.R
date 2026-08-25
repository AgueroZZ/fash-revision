#!/usr/bin/env Rscript

# Plot paired functional-LFSR agreement and numerical errors across nested
# evaluation-grid resolutions.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

suppressPackageStartupMessages(library(ggplot2))

workflowr_root <- find_workflowr_root()
output_dir <- get_arg(
  "--output-dir",
  file.path(
    workflowr_root,
    "output",
    "revision_simulations",
    "internal",
    "evaluation_grid_lfsr_sensitivity_pilot"
  )
)
results_path <- file.path(output_dir, "pair_category_lfsr_by_grid.csv")
summary_path <- file.path(output_dir, "grid_sensitivity_summary.csv")
if (!file.exists(results_path) || !file.exists(summary_path)) {
  stop("The evaluation-grid sensitivity results are incomplete.")
}

results <- read.csv(results_path, stringsAsFactors = FALSE)
summary <- read.csv(summary_path, stringsAsFactors = FALSE)
category_order <- c("early", "middle", "late", "switch")
category_labels <- c(
  early = "Early",
  middle = "Middle",
  late = "Late",
  switch = "Switch"
)
results$category <- factor(results$category, levels = category_order, labels = category_labels)
summary$category <- factor(summary$category, levels = category_order, labels = category_labels)

current_step <- 0.10
fine_step <- 0.05
scatter_data <- results[abs(results$grid_step - current_step) < 1e-12, , drop = FALSE]
scatter_summary <- summary[abs(summary$grid_step - current_step) < 1e-12, , drop = FALSE]
scatter_summary$label <- paste0(
  "Mean |diff| = ",
  formatC(
    scatter_summary$mean_absolute_difference_vs_fine,
    format = "f",
    digits = 4
  ),
  "\n90th pct = ",
  formatC(
    scatter_summary$q90_absolute_difference_vs_fine,
    format = "f",
    digits = 4
  ),
  "\nSpearman rho = ",
  formatC(scatter_summary$spearman_vs_fine, format = "f", digits = 3)
)

observed_max <- max(scatter_data$lfsr, scatter_data$reference_lfsr)
common_limit <- min(1, max(0.12, ceiling(observed_max * 20) / 20))
scatter_summary$x <- 0.97 * common_limit
scatter_summary$y <- 0.97 * common_limit
category_colors <- c(
  Early = "#0072B2",
  Middle = "#009E73",
  Late = "#CC79A7",
  Switch = "#D55E00"
)

scatter_plot <- ggplot(
  scatter_data,
  aes(x = lfsr, y = reference_lfsr, color = category)
) +
  geom_abline(slope = 1, intercept = 0, color = "#666666", linewidth = 0.6) +
  geom_point(size = 1.8, alpha = 0.62, stroke = 0) +
  geom_label(
    data = scatter_summary,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1,
    size = 2.8,
    lineheight = 1.05,
    color = "#1F2937",
    fill = scales::alpha("white", 0.92),
    linewidth = 0.25,
    label.padding = grid::unit(0.20, "lines")
  ) +
  facet_wrap(~category, ncol = 2) +
  scale_color_manual(values = category_colors, guide = "none") +
  scale_x_continuous(
    limits = c(0, common_limit),
    breaks = scales::breaks_pretty(n = 5),
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(0, common_limit),
    breaks = scales::breaks_pretty(n = 5),
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = 0)
  ) +
  coord_equal() +
  labs(
    title = "Functional LFSR agreement under evaluation-grid refinement",
    subtitle = "Current 0.10-day grid versus the finer 0.05-day reference grid",
    x = "Functional LFSR: 0.10-day grid",
    y = "Functional LFSR: 0.05-day grid",
    caption = paste(
      "Pairs are restricted to the primary reported category sets; both grids",
      "use common posterior draws."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#4B5563"),
    plot.caption = element_text(color = "#4B5563", hjust = 0),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.35),
    panel.spacing = grid::unit(1.0, "lines")
  )

scatter_path <- file.path(output_dir, "grid_resolution_lfsr_scatter.png")
ggsave(
  filename = scatter_path,
  plot = scatter_plot,
  width = 9.5,
  height = 8.2,
  units = "in",
  dpi = 200,
  bg = "white"
)

difference_data <- results[results$grid_step > fine_step + 1e-12, , drop = FALSE]
difference_data$grid_label <- factor(
  format(difference_data$grid_step, nsmall = 2),
  levels = c("0.15", "0.10"),
  labels = c("0.15-day grid", "0.10-day grid")
)

difference_plot <- ggplot(
  difference_data,
  aes(x = grid_label, y = absolute_difference_vs_fine, fill = category)
) +
  geom_violin(
    scale = "width",
    trim = TRUE,
    alpha = 0.25,
    color = NA
  ) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    alpha = 0.60,
    linewidth = 0.35
  ) +
  geom_jitter(
    aes(color = category),
    width = 0.08,
    height = 0,
    size = 0.7,
    alpha = 0.16,
    show.legend = FALSE
  ) +
  facet_wrap(~category, nrow = 1) +
  scale_fill_manual(values = category_colors, guide = "none") +
  scale_color_manual(values = category_colors, guide = "none") +
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.001),
    expand = expansion(mult = c(0.03, 0.10))
  ) +
  labs(
    title = "Absolute functional-LFSR error relative to the 0.05-day grid",
    subtitle = "Smaller values indicate numerical stability under grid refinement",
    x = NULL,
    y = "Absolute LFSR difference versus 0.05-day grid",
    caption = "Distributions are across primary reported pairs within each category."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#4B5563"),
    plot.caption = element_text(color = "#4B5563", hjust = 0),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E7EB", linewidth = 0.35),
    axis.text.x = element_text(angle = 18, hjust = 1)
  )

difference_path <- file.path(output_dir, "grid_resolution_absolute_difference.png")
ggsave(
  filename = difference_path,
  plot = difference_plot,
  width = 12,
  height = 5.4,
  units = "in",
  dpi = 200,
  bg = "white"
)

cat(scatter_path, "\n")
cat(difference_path, "\n")
