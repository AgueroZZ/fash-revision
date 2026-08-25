#!/usr/bin/env Rscript

# Plot paired middle-functional LFSR agreement and absolute numerical changes
# for the fixed most-significant dynamic candidate pair per gene.

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
    "middle_gene_representative_grid_sensitivity_pilot"
  )
)
results_path <- file.path(output_dir, "middle_lfsr_grid_sensitivity_by_gene.csv")
summary_path <- file.path(output_dir, "middle_lfsr_grid_sensitivity_summary.csv")
if (!file.exists(results_path) || !file.exists(summary_path)) {
  stop("The middle gene-representative grid-sensitivity results are incomplete.")
}

results <- read.csv(results_path, stringsAsFactors = FALSE)
summary <- read.csv(summary_path, stringsAsFactors = FALSE)
overall <- summary[summary$stratum == "All genes", , drop = FALSE]
if (nrow(overall) != 1L) {
  stop("The overall sensitivity summary is missing or duplicated.")
}

results$saved_lfsr_stratum <- factor(
  results$saved_lfsr_stratum,
  levels = c("[0, 0.05]", "(0.05, 0.10]", "(0.10, 0.25]", "(0.25, 0.50]", "(0.50, 1]")
)
stratum_colors <- c(
  "[0, 0.05]" = "#D55E00",
  "(0.05, 0.10]" = "#E69F00",
  "(0.10, 0.25]" = "#009E73",
  "(0.25, 0.50]" = "#0072B2",
  "(0.50, 1]" = "#9CA3AF"
)

annotation <- paste0(
  "n = ", format(overall$n_genes, big.mark = ","), " genes",
  "\nMean increase = ", formatC(overall$mean_absolute_difference, format = "f", digits = 4),
  "\n90th percentile = ", formatC(overall$q90_absolute_difference, format = "f", digits = 4),
  "\nSpearman rho = ", formatC(overall$spearman_0p10_vs_0p05, format = "f", digits = 3)
)

scatter_plot <- ggplot(
  results,
  aes(x = lfsr_0p10, y = lfsr_0p05, color = saved_lfsr_stratum)
) +
  geom_abline(slope = 1, intercept = 0, color = "#4B5563", linewidth = 0.65) +
  geom_point(size = 1.65, alpha = 0.50, stroke = 0) +
  annotate(
    "label",
    x = 0.98,
    y = 0.02,
    label = annotation,
    hjust = 1,
    vjust = 0,
    size = 3.2,
    lineheight = 1.08,
    color = "#1F2937",
    fill = scales::alpha("white", 0.94),
    linewidth = 0.25,
    label.padding = grid::unit(0.25, "lines")
  ) +
  scale_color_manual(
    values = stratum_colors,
    drop = FALSE,
    name = "Saved 0.10-grid LFSR"
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = 0)
  ) +
  coord_equal() +
  labs(
    title = "Middle-functional LFSR shifts upward under grid refinement",
    subtitle = "One fixed, most-significant dynamic candidate pair per gene",
    x = "Middle-functional LFSR: 0.10-day grid",
    y = "Middle-functional LFSR: 0.05-day grid",
    caption = paste(
      "Pair selection uses saved 0.10-grid LFSR; paired grid estimates use",
      "common posterior draws."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#4B5563"),
    plot.caption = element_text(color = "#4B5563", hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.35),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.key.height = grid::unit(0.55, "lines")
  )

scatter_path <- file.path(output_dir, "middle_gene_representative_lfsr_scatter.png")
ggsave(
  filename = scatter_path,
  plot = scatter_plot,
  width = 8.8,
  height = 7.2,
  units = "in",
  dpi = 220,
  bg = "white"
)

difference_plot <- ggplot(
  results,
  aes(
    x = saved_lfsr_stratum,
    y = lfsr_0p05 - lfsr_0p10,
    fill = saved_lfsr_stratum
  )
) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.28, color = NA) +
  geom_boxplot(
    width = 0.18,
    outlier.shape = NA,
    alpha = 0.70,
    linewidth = 0.35
  ) +
  geom_jitter(
    width = 0.08,
    height = 0,
    size = 0.65,
    alpha = 0.16,
    color = "#374151"
  ) +
  scale_fill_manual(values = stratum_colors, guide = "none", drop = FALSE) +
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.005),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  labs(
    title = "LFSR increase from the 0.10-day to 0.05-day grid",
    subtitle = "Differences are shown by the saved LFSR used for pair selection",
    x = "Saved middle-functional LFSR on the 0.10-day grid",
    y = "Paired LFSR increase (0.05-day minus 0.10-day)",
    caption = "Each point represents the fixed representative pair for one gene."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#4B5563"),
    plot.caption = element_text(color = "#4B5563", hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "#E5E7EB", linewidth = 0.35),
    axis.text.x = element_text(angle = 18, hjust = 1)
  )

difference_path <- file.path(output_dir, "middle_gene_representative_lfsr_increase.png")
ggsave(
  filename = difference_path,
  plot = difference_plot,
  width = 9.5,
  height = 5.8,
  units = "in",
  dpi = 220,
  bg = "white"
)

cat(scatter_path, "\n")
cat(difference_path, "\n")
