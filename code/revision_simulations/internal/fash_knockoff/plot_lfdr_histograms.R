#!/usr/bin/env Rscript

# Writes the target-versus-decoy lfdr figures and band table to disk. The
# figure is defined in reporting.R; this script only renders it.

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

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "internal",
                 "fash_knockoff", "fash_knockoff_helpers.R"))
source(file.path(workflowr_root, "code", "revision_simulations", "internal",
                 "fash_knockoff", "reporting.R"))
library(ggplot2)

statistics <- load_fash_knockoff_statistics(workflowr_root)
figure_directory <- file.path(
  fash_knockoff_paths(workflowr_root)$summary, "figures"
)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

bands <- fash_knockoff_lfdr_bands(statistics)
print(bands, row.names = FALSE)
utils::write.csv(bands, file.path(figure_directory, "lfdr_band_ratios.csv"),
                 row.names = FALSE)

# Full-range and log-count views, kept for diagnosis; the band panel is the one
# that carries the result.
long <- fash_knockoff_lfdr_long(statistics)
linear_panel <- ggplot(long, aes(x = lfdr, fill = arm)) +
  geom_histogram(binwidth = 0.02, boundary = 0, position = "identity",
                 alpha = 0.5, colour = NA) +
  scale_fill_manual(values = FASH_KNOCKOFF_ARM_COLOURS) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
  labs(x = "lfdr (merged fit)", y = "gene-variant pairs",
       title = "a  Full lfdr range: the two arms are indistinguishable",
       subtitle = sprintf("%s pairs per arm",
                          format(nrow(statistics), big.mark = ","))) +
  fash_knockoff_theme() + theme(legend.position = "top")
log_panel <- ggplot(long, aes(x = log10(lfdr), fill = arm)) +
  geom_histogram(binwidth = 0.25, position = "identity", alpha = 0.5,
                 colour = NA) +
  geom_vline(xintercept = log10(0.05), linetype = "dashed",
             colour = "grey30", linewidth = 0.4) +
  scale_fill_manual(values = FASH_KNOCKOFF_ARM_COLOURS) +
  scale_y_log10(labels = function(x) format(x, big.mark = ",",
                                            scientific = FALSE)) +
  labs(x = expression(log[10] ~ "lfdr (merged fit)"),
       y = "gene-variant pairs (log scale)",
       title = "b  Log counts: the tail becomes visible") +
  fash_knockoff_theme() + theme(legend.position = "none")
band_panel <- plot_lfdr_bands(statistics) +
  labs(title = "c  Only the left tail, counted") +
  theme(legend.position = "none")

combined <- (linear_panel / log_panel / band_panel)
ggsave(file.path(figure_directory, "lfdr_histograms_target_vs_decoy.png"),
       combined, width = 8.5, height = 11, dpi = 200)
ggsave(file.path(figure_directory, "lfdr_histograms_target_vs_decoy.pdf"),
       combined, width = 8.5, height = 11)
message("Wrote figures to ", figure_directory)
