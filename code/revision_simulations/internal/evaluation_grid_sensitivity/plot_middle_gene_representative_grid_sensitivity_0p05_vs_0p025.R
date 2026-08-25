#!/usr/bin/env Rscript

# Plot further middle-functional LFSR grid refinement and compare its change
# with the preceding 0.10-day to 0.05-day refinement.

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
    "middle_gene_representative_grid_sensitivity_0p05_vs_0p025_pilot"
  )
)
results_path <- file.path(output_dir, "middle_lfsr_grid_sensitivity_by_gene.csv")
summary_path <- file.path(output_dir, "middle_lfsr_grid_sensitivity_summary.csv")
if (!file.exists(results_path) || !file.exists(summary_path)) {
  stop("The 0.05 versus 0.025 grid-sensitivity results are incomplete.")
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
  "\nMean |difference| = ", formatC(overall$mean_absolute_difference, format = "f", digits = 4),
  "\n90th percentile = ", formatC(overall$q90_absolute_difference, format = "f", digits = 4),
  "\nSpearman rho = ", formatC(overall$spearman_0p05_vs_0p025, format = "f", digits = 5)
)

scatter_plot <- ggplot(
  results,
  aes(x = lfsr_0p05, y = lfsr_0p025, color = saved_lfsr_stratum)
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
    title = "Middle-functional LFSR is unchanged below 0.05-day resolution",
    subtitle = "The same fixed gene-representative pairs are compared at both resolutions",
    x = "Middle-functional LFSR: 0.05-day grid",
    y = "Middle-functional LFSR: 0.025-day grid",
    caption = "Both grid estimates use common posterior draws on the nested 0.025-day grid."
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

observed_max <- max(
  results$previous_absolute_difference_0p10_to_0p05,
  results$absolute_difference
)
comparison_limit <- max(0.05, ceiling(observed_max * 100) / 100)
convergence_annotation <- paste0(
  "Mean |difference| ratio = ",
  formatC(overall$ratio_of_mean_absolute_differences, format = "f", digits = 3),
  "\nNew difference smaller for ",
  scales::percent(
    overall$fraction_new_difference_smaller_than_previous,
    accuracy = 0.1
  ),
  " of genes"
)

convergence_plot <- ggplot(
  results,
  aes(
    x = previous_absolute_difference_0p10_to_0p05,
    y = absolute_difference,
    color = saved_lfsr_stratum
  )
) +
  geom_abline(slope = 1, intercept = 0, color = "#4B5563", linewidth = 0.65) +
  geom_point(size = 1.65, alpha = 0.50, stroke = 0) +
  annotate(
    "label",
    x = 0.98 * comparison_limit,
    y = 0.02 * comparison_limit,
    label = convergence_annotation,
    hjust = 1,
    vjust = 0,
    size = 3.2,
    lineheight = 1.08,
    color = "#1F2937",
    fill = scales::alpha("white", 0.94),
    linewidth = 0.25,
    label.padding = grid::unit(0.25, "lines")
  ) +
  scale_color_manual(values = stratum_colors, drop = FALSE, guide = "none") +
  scale_x_continuous(
    limits = c(0, comparison_limit),
    breaks = scales::breaks_pretty(n = 5),
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = 0)
  ) +
  scale_y_continuous(
    limits = c(0, comparison_limit),
    breaks = scales::breaks_pretty(n = 5),
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = 0)
  ) +
  coord_equal() +
  labs(
    title = "Further grid refinement produces much smaller LFSR changes",
    subtitle = "Points below the diagonal have a smaller change in the second refinement",
    x = "Absolute LFSR change: 0.10-day to 0.05-day grid",
    y = "Absolute LFSR change: 0.05-day to 0.025-day grid",
    caption = "Each point represents the same fixed pair for one gene."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#4B5563"),
    plot.caption = element_text(color = "#4B5563", hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.35)
  )

convergence_path <- file.path(output_dir, "middle_gene_representative_convergence_comparison.png")
ggsave(
  filename = convergence_path,
  plot = convergence_plot,
  width = 8.0,
  height = 7.2,
  units = "in",
  dpi = 220,
  bg = "white"
)

cat(scatter_path, "\n")
cat(convergence_path, "\n")
