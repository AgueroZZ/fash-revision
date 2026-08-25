#!/usr/bin/env Rscript

# Histograms of the merged-fit lfdr for target and decoy units.
#
# Both sets of lfdr values come from the SAME merged empirical-Bayes fit, so
# they are directly comparable: any difference between the two histograms is a
# difference in evidence, not in calibration.

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

for (package in c("ggplot2", "patchwork")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The ", package, " package is required.")
  }
}
library(ggplot2)

workflowr_root <- find_workflowr_root()
statistics_path <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real",
  "fash_knockoff_full_seed20260823", "pair_statistics.rds"
)
figure_directory <- file.path(
  workflowr_root, "output", "revision_simulations", "internal",
  "fash_knockoff_full_seed20260823_from_midway3", "figures"
)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

statistics <- readRDS(statistics_path)
n_pair <- nrow(statistics)
lfdr <- data.frame(
  arm = rep(c("target (observed genotype)", "decoy (permuted genotype)"),
            each = n_pair),
  lfdr = c(statistics$target_lfdr, statistics$decoy_lfdr)
)
lfdr$arm <- factor(lfdr$arm, levels = c(
  "target (observed genotype)", "decoy (permuted genotype)"
))
arm_colours <- c(
  "target (observed genotype)" = "#1f6f8b",
  "decoy (permuted genotype)" = "#c1462f"
)

summary_table <- do.call(rbind, lapply(levels(lfdr$arm), function(arm) {
  values <- lfdr$lfdr[lfdr$arm == arm]
  data.frame(
    arm = arm,
    n = length(values),
    median = median(values),
    below_0.05 = sum(values < 0.05),
    below_1e_6 = sum(values < 1e-6),
    stringsAsFactors = FALSE
  )
}))
print(summary_table, row.names = FALSE)
utils::write.csv(
  summary_table,
  file.path(figure_directory, "lfdr_histogram_summary.csv"),
  row.names = FALSE
)

base_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11)
  )

linear_panel <- ggplot(lfdr, aes(x = lfdr, fill = arm)) +
  geom_histogram(
    binwidth = 0.02, boundary = 0, position = "identity", alpha = 0.5,
    colour = NA
  ) +
  scale_fill_manual(values = arm_colours) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
  labs(
    x = "lfdr (merged fit)", y = "gene-variant pairs",
    title = "a  Full lfdr range \u2014 the two arms are indistinguishable",
    subtitle = paste0(
      format(n_pair, big.mark = ","),
      " pairs per arm; merged pi0 = 0.942; medians 0.94477 vs 0.94497"
    )
  ) +
  base_theme

# A log y-axis is essential: the pile at lfdr near 1 is four orders of
# magnitude taller than the tail that contains every discovery.
log_panel <- ggplot(lfdr, aes(x = log10(lfdr), fill = arm)) +
  geom_histogram(
    binwidth = 0.25, position = "identity", alpha = 0.5, colour = NA
  ) +
  geom_vline(xintercept = log10(0.05), linetype = "dashed",
             colour = "grey30", linewidth = 0.4) +
  annotate("text", x = log10(0.05) - 0.3, y = 3e5, label = "lfdr = 0.05",
           hjust = 1, size = 3, colour = "grey30") +
  scale_fill_manual(values = arm_colours) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
  labs(
    x = expression(log[10] ~ "lfdr (merged fit)"),
    y = "gene-variant pairs (log scale)",
    title = "b  Same data, log counts \u2014 the tail becomes visible"
  ) +
  base_theme

# The decision-relevant region, counted rather than eyeballed.
tail_breaks <- c(-Inf, -12, -9, -6, -4, -3, -2, log10(0.05), -1, 0)
tail_labels <- c("<1e-12", "1e-12..1e-9", "1e-9..1e-6", "1e-6..1e-4",
                 "1e-4..1e-3", "1e-3..0.01", "0.01..0.05", "0.05..0.1",
                 "0.1..1")
tail <- lfdr
tail$bin <- cut(log10(tail$lfdr), breaks = tail_breaks, labels = tail_labels)
tail <- subset(tail, !is.na(bin) & bin != "0.1..1")
tail$bin <- droplevels(tail$bin)
tail_counts <- as.data.frame(table(arm = tail$arm, bin = tail$bin))
band_ratio <- reshape(tail_counts, idvar = "bin", timevar = "arm",
                      direction = "wide")
names(band_ratio) <- c("bin", "target", "decoy")
band_ratio$decoy_over_target <- round(band_ratio$decoy / band_ratio$target, 3)

tail_panel <- ggplot(tail_counts, aes(x = bin, y = Freq, fill = arm)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7,
           alpha = 0.5) +
  geom_text(aes(label = Freq), position = position_dodge(width = 0.75),
            vjust = -0.35, size = 2.7) +
  scale_fill_manual(values = arm_colours) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "lfdr band", y = "gene-variant pairs",
    title = "c  Only the left tail, counted",
    subtitle = paste0(
      "lfdr < 0.05: ", summary_table$below_0.05[1], " target vs ",
      summary_table$below_0.05[2], " decoy  (ratio ",
      round(summary_table$below_0.05[2] / summary_table$below_0.05[1], 3), ")"
    )
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

combined <- (linear_panel / log_panel / tail_panel) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "top")
ggsave(
  file.path(figure_directory, "lfdr_histograms_target_vs_decoy.png"),
  combined, width = 8.5, height = 11, dpi = 200
)
ggsave(
  file.path(figure_directory, "lfdr_histograms_target_vs_decoy.pdf"),
  combined, width = 8.5, height = 11
)
utils::write.csv(
  tail_counts, file.path(figure_directory, "lfdr_tail_counts.csv"),
  row.names = FALSE
)
print(band_ratio, row.names = FALSE)
utils::write.csv(band_ratio,
                 file.path(figure_directory, "lfdr_band_ratios.csv"),
                 row.names = FALSE)
message("Wrote figures to ", figure_directory)
