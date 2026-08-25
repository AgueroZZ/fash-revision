#!/usr/bin/env Rscript

# Plot the local functional-LFSR profile around the primary day-3 cutoff.

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
    "early_window_lfsr_local_profile_pilot"
  )
)
profile_path <- file.path(output_dir, "early_window_pair_lfsr_profile.csv")
summary_path <- file.path(output_dir, "early_window_profile_summary.csv")
if (!file.exists(profile_path) || !file.exists(summary_path)) {
  stop("The early-window local-profile results are incomplete.")
}

profile <- read.csv(profile_path, stringsAsFactors = FALSE)
summary <- read.csv(summary_path, stringsAsFactors = FALSE)
profile_required <- c("pair_id", "cutoff_day", "lfsr", "delta_vs_primary")
summary_required <- c(
  "cutoff_day", "mean_lfsr", "mean_delta", "q10_delta", "q25_delta",
  "median_delta", "q75_delta", "q90_delta"
)
if (!all(profile_required %in% names(profile)) ||
    !all(summary_required %in% names(summary))) {
  stop("The early-window local-profile tables are missing required columns.")
}
if (any(!is.finite(profile$delta_vs_primary)) ||
    any(!is.finite(as.matrix(summary[, summary_required[-1L]])))) {
  stop("The early-window local-profile tables contain non-finite values.")
}

primary_cutoff <- 3
endpoint_labels <- summary[
  summary$cutoff_day %in% c(min(summary$cutoff_day), max(summary$cutoff_day)),
  c("cutoff_day", "mean_delta"),
  drop = FALSE
]
endpoint_labels$label <- paste0(
  "Mean delta = ",
  sprintf("%+.3f", endpoint_labels$mean_delta)
)
endpoint_labels$hjust <- ifelse(endpoint_labels$cutoff_day < primary_cutoff, 0, 1)
endpoint_labels$x_offset <- ifelse(
  endpoint_labels$cutoff_day < primary_cutoff,
  0.015,
  -0.015
)

plot_object <- ggplot() +
  geom_hline(yintercept = 0, color = "#4B5563", linewidth = 0.55) +
  geom_vline(
    xintercept = primary_cutoff,
    color = "#4B5563",
    linewidth = 0.55,
    linetype = "dashed"
  ) +
  geom_line(
    data = profile,
    aes(x = cutoff_day, y = delta_vs_primary, group = pair_id),
    color = "#9CA3AF",
    linewidth = 0.24,
    alpha = 0.16
  ) +
  geom_ribbon(
    data = summary,
    aes(x = cutoff_day, ymin = q10_delta, ymax = q90_delta),
    fill = "#56B4E9",
    alpha = 0.18
  ) +
  geom_ribbon(
    data = summary,
    aes(x = cutoff_day, ymin = q25_delta, ymax = q75_delta),
    fill = "#0072B2",
    alpha = 0.24
  ) +
  geom_line(
    data = summary,
    aes(x = cutoff_day, y = mean_delta),
    color = "#0072B2",
    linewidth = 1.15
  ) +
  geom_point(
    data = summary,
    aes(x = cutoff_day, y = mean_delta),
    color = "#0072B2",
    size = 2.2
  ) +
  geom_label(
    data = endpoint_labels,
    aes(
      x = cutoff_day + x_offset,
      y = mean_delta,
      label = label,
      hjust = hjust
    ),
    size = 3.3,
    color = "#1F2937",
    fill = scales::alpha("white", 0.92),
    linewidth = 0.25,
    label.padding = grid::unit(0.22, "lines")
  ) +
  scale_x_continuous(
    breaks = seq(min(summary$cutoff_day), max(summary$cutoff_day), by = 0.1),
    labels = scales::label_number(accuracy = 0.1),
    expand = expansion(mult = c(0.025, 0.025))
  ) +
  scale_y_continuous(
    labels = scales::label_number(accuracy = 0.01),
    expand = expansion(mult = c(0.10, 0.16))
  ) +
  labs(
    title = "Local sensitivity of early-category LFSR to the time-window cutoff",
    subtitle = paste0(
      "Conditional on ",
      length(unique(profile$pair_id)),
      " primary early pairs; thin lines are pair-level profiles"
    ),
    x = "Early-window cutoff day",
    y = "Change in functional LFSR relative to cutoff day 3",
    caption = paste(
      "Dark and light ribbons show the pairwise 25th-75th and 10th-90th",
      "percentiles, respectively; they are not confidence intervals."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "#4B5563"),
    plot.caption = element_text(color = "#4B5563", hjust = 0),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "#E5E7EB", linewidth = 0.35),
    axis.text.x = element_text(size = 9)
  )

output_path <- file.path(output_dir, "early_window_lfsr_local_profile.png")
ggsave(
  filename = output_path,
  plot = plot_object,
  width = 10,
  height = 6.2,
  units = "in",
  dpi = 200,
  bg = "white"
)
cat(output_path, "\n")

