#!/usr/bin/env Rscript

# Plot paired LFSR sensitivity for alternative early-window cutoffs.

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
    "early_window_lfsr_sensitivity_pilot"
  )
)
input_path <- file.path(output_dir, "early_window_pair_lfsr.csv")
if (!file.exists(input_path)) {
  stop("Missing early-window sensitivity results: ", input_path)
}

results <- read.csv(input_path, stringsAsFactors = FALSE)
required_columns <- c("lfsr_cutoff_2p5", "lfsr_cutoff_3", "lfsr_cutoff_3p5")
if (!all(required_columns %in% names(results))) {
  stop("The sensitivity results are missing required LFSR columns.")
}

plot_data <- rbind(
  data.frame(
    pair_id = results$pair_id,
    primary_lfsr = results$lfsr_cutoff_3,
    alternative_lfsr = results$lfsr_cutoff_2p5,
    scenario = "Narrower window: days 0-2.5",
    stringsAsFactors = FALSE
  ),
  data.frame(
    pair_id = results$pair_id,
    primary_lfsr = results$lfsr_cutoff_3,
    alternative_lfsr = results$lfsr_cutoff_3p5,
    scenario = "Wider window: days 0-3.5",
    stringsAsFactors = FALSE
  )
)
plot_data$scenario <- factor(
  plot_data$scenario,
  levels = c("Narrower window: days 0-2.5", "Wider window: days 0-3.5")
)

if (any(!is.finite(plot_data$primary_lfsr)) ||
    any(!is.finite(plot_data$alternative_lfsr)) ||
    any(plot_data$primary_lfsr < 0 | plot_data$primary_lfsr > 1) ||
    any(plot_data$alternative_lfsr < 0 | plot_data$alternative_lfsr > 1)) {
  stop("LFSR values must be finite and lie in [0, 1].")
}

annotation <- do.call(rbind, lapply(levels(plot_data$scenario), function(label) {
  current <- plot_data[plot_data$scenario == label, , drop = FALSE]
  data.frame(
    scenario = factor(label, levels = levels(plot_data$scenario)),
    label = paste0(
      "Spearman rho = ",
      formatC(cor(
        current$primary_lfsr,
        current$alternative_lfsr,
        method = "spearman"
      ), format = "f", digits = 3),
      "\nMean delta LFSR = ",
      sprintf("%+.3f", mean(current$alternative_lfsr - current$primary_lfsr)),
      "\nMean alternative LFSR = ",
      formatC(mean(current$alternative_lfsr), format = "f", digits = 3)
    ),
    stringsAsFactors = FALSE
  )
}))

observed_max <- max(plot_data$primary_lfsr, plot_data$alternative_lfsr)
common_limit <- min(1, max(0.12, ceiling(observed_max * 20) / 20))
annotation$x <- 0.97 * common_limit
annotation$y <- 0.97 * common_limit

scenario_colors <- c(
  "Narrower window: days 0-2.5" = "#0072B2",
  "Wider window: days 0-3.5" = "#D55E00"
)

plot_object <- ggplot(
  plot_data,
  aes(x = primary_lfsr, y = alternative_lfsr, color = scenario)
) +
  geom_abline(slope = 1, intercept = 0, color = "#666666", linewidth = 0.6) +
  geom_point(size = 2.1, alpha = 0.72, stroke = 0) +
  geom_label(
    data = annotation,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 1,
    size = 3.2,
    lineheight = 1.08,
    linewidth = 0.25,
    label.padding = grid::unit(0.25, "lines"),
    color = "#1F2937",
    fill = scales::alpha("white", 0.92)
  ) +
  facet_wrap(~scenario, nrow = 1) +
  scale_color_manual(values = scenario_colors, guide = "none") +
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
    title = "Early-category support under half-day window changes",
    subtitle = paste0(
      "Primary early-category pairs only (n = ",
      nrow(results),
      "); both axes use common posterior draws"
    ),
    x = "LFSR: primary early window (days 0-3)",
    y = "LFSR: alternative early window",
    caption = paste(
      "This conditional sensitivity analysis assesses reported-pair stability;",
      "it does not identify new pairs under alternative definitions."
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
    panel.spacing = grid::unit(1.1, "lines")
  )

output_path <- file.path(output_dir, "early_window_lfsr_scatter.png")
ggsave(
  filename = output_path,
  plot = plot_object,
  width = 11.5,
  height = 5.3,
  units = "in",
  dpi = 200,
  bg = "white"
)
cat(output_path, "\n")
