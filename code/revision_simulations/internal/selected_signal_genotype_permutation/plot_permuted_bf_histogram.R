#!/usr/bin/env Rscript

# Plot the Bayes-factor distribution for permuted negative-control units.

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

workflowr_root <- find_workflowr_root()
output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "selected_signal_genotype_permutation_seed20260811"
)
fit_path <- file.path(output_directory, "merged_fash_fit.rds")
figure_directory <- file.path(output_directory, "figures")
if (!file.exists(fit_path)) {
  stop("The fixed-seed merged FASH fit is missing.")
}
dir.create(figure_directory, showWarnings = FALSE)

fit_bundle <- readRDS(fit_path)
n_target <- fit_bundle$configuration$n_target_units
bayes_factor <- fit_bundle$bf_adjusted_fit$BF
if (n_target != 1177L || length(bayes_factor) != 2L * n_target ||
    any(!is.finite(bayes_factor)) || any(bayes_factor <= 0)) {
  stop("The saved Bayes factors are not the expected fixed-seed results.")
}

permuted_bf <- bayes_factor[n_target + seq_len(n_target)]
permuted_data <- data.frame(
  permuted_unit_index = seq_len(n_target),
  source_pair_key = fit_bundle$selection$pair_key,
  bayes_factor = permuted_bf,
  log10_bayes_factor = log10(permuted_bf),
  stringsAsFactors = FALSE
)

descending_order <- order(permuted_bf, decreasing = TRUE)
total_bf_mass <- sum(permuted_bf)
top_mass_share <- function(proportion) {
  n_top <- ceiling(length(permuted_bf) * proportion)
  sum(permuted_bf[descending_order[seq_len(n_top)]]) / total_bf_mass
}
trimmed_mean <- function(proportion) {
  n_top <- ceiling(length(permuted_bf) * proportion)
  mean(permuted_bf[-descending_order[seq_len(n_top)]])
}

summary_data <- data.frame(
  n_permuted = length(permuted_bf),
  mean_bf = mean(permuted_bf),
  median_bf = stats::median(permuted_bf),
  sd_bf = stats::sd(permuted_bf),
  minimum_bf = min(permuted_bf),
  maximum_bf = max(permuted_bf),
  proportion_bf_greater_than_one = mean(permuted_bf > 1),
  top_1_percent_bf_mass_share = top_mass_share(0.01),
  top_5_percent_bf_mass_share = top_mass_share(0.05),
  mean_after_removing_top_1_percent = trimmed_mean(0.01),
  mean_after_removing_top_5_percent = trimmed_mean(0.05),
  stringsAsFactors = FALSE
)

utils::write.csv(
  permuted_data,
  file.path(output_directory, "permuted_bayes_factors.csv"),
  row.names = FALSE
)
utils::write.csv(
  summary_data,
  file.path(output_directory, "permuted_bayes_factor_summary.csv"),
  row.names = FALSE
)

reference_data <- data.frame(
  value = log10(c(stats::median(permuted_bf), 1, mean(permuted_bf))),
  label = c(
    sprintf("Median = %.3f", stats::median(permuted_bf)),
    "Null reference = 1",
    sprintf("Arithmetic mean = %.1f", mean(permuted_bf))
  ),
  line_type = factor(
    c("Median", "Null reference", "Arithmetic mean"),
    levels = c("Median", "Null reference", "Arithmetic mean")
  ),
  stringsAsFactors = FALSE
)

histogram_counts <- hist(
  permuted_data$log10_bayes_factor,
  breaks = seq(-1.6, 4.8, by = 0.16),
  plot = FALSE
)$counts
maximum_count <- max(histogram_counts)

figure <- ggplot2::ggplot(
  permuted_data,
  ggplot2::aes(x = log10_bayes_factor)
) +
  ggplot2::geom_histogram(
    breaks = seq(-1.6, 4.8, by = 0.16),
    fill = "#4C78A8",
    color = "white",
    linewidth = 0.25
  ) +
  ggplot2::geom_vline(
    data = reference_data,
    ggplot2::aes(xintercept = value, linetype = line_type),
    color = "#252525",
    linewidth = 0.8
  ) +
  ggplot2::annotate(
    "text",
    x = reference_data$value,
    y = c(0.94, 0.78, 0.62) * maximum_count,
    label = reference_data$label,
    angle = 90,
    hjust = 1,
    vjust = c(-0.35, -0.35, -0.35),
    size = 3.5,
    color = "#252525"
  ) +
  ggplot2::annotate(
    "label",
    x = 4.55,
    y = 0.92 * maximum_count,
    hjust = 1,
    vjust = 1,
    label = sprintf(
      paste0(
        "Top 1%% of units contribute %.1f%% of total BF\n",
        "Top 5%% contribute %.1f%%\n",
        "%.1f%% of units have BF > 1"
      ),
      100 * summary_data$top_1_percent_bf_mass_share,
      100 * summary_data$top_5_percent_bf_mass_share,
      100 * summary_data$proportion_bf_greater_than_one
    ),
    size = 3.6,
    lineheight = 1.15,
    color = "#252525",
    fill = "white",
    linewidth = 0.25
  ) +
  ggplot2::scale_linetype_manual(
    values = c(
      `Median` = "dashed",
      `Null reference` = "dotted",
      `Arithmetic mean` = "dotdash"
    ),
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(
    breaks = -1:4,
    labels = c("0.1", "1", "10", "100", "1,000", "10,000"),
    expand = ggplot2::expansion(mult = c(0.01, 0.01))
  ) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(mult = c(0, 0.05))
  ) +
  ggplot2::labs(
    title = "A small right tail drives the permuted-control mean BF",
    subtitle = paste(
      "1,177 genotype-permuted units from the fixed seed 20260811;",
      "the x-axis is logarithmic"
    ),
    x = "Bayes factor for IWP alternative vs exact null (log10 scale)",
    y = "Number of permuted units",
    caption = paste(
      "BF uses the alternative mixture learned in the merged original-plus-permuted fit.",
      "Vertical lines mark the median, BF = 1, and arithmetic mean."
    )
  ) +
  ggplot2::coord_cartesian(xlim = c(-1.6, 4.8)) +
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
    axis.title = ggplot2::element_text(size = 11.5, color = "#222222"),
    axis.text = ggplot2::element_text(size = 10, color = "#333333"),
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
    panel.grid.major.y = ggplot2::element_line(
      color = "#E6E6E6",
      linewidth = 0.4
    ),
    panel.grid.minor = ggplot2::element_blank(),
    plot.margin = ggplot2::margin(14, 18, 12, 14)
  )

ggplot2::ggsave(
  file.path(figure_directory, "permuted_bayes_factor_histogram.png"),
  figure,
  width = 10.5,
  height = 6.5,
  dpi = 320,
  bg = "white"
)
ggplot2::ggsave(
  file.path(figure_directory, "permuted_bayes_factor_histogram.pdf"),
  figure,
  width = 10.5,
  height = 6.5,
  device = grDevices::cairo_pdf,
  bg = "white"
)

print(summary_data)
