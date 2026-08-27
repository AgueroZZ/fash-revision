#!/usr/bin/env Rscript

# Writes the W-by-discovery figure to disk. The figure itself is defined in
# reporting.R so that the workflowr page and this script cannot drift apart.

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

summary_table <- fash_knockoff_group_summary(statistics)
print(summary_table, row.names = FALSE, digits = 4)
utils::write.csv(
  summary_table, file.path(figure_directory, "w_by_discovery_summary.csv"),
  row.names = FALSE
)

figure <- plot_w_by_discovery(statistics) +
  patchwork::plot_annotation(
    caption = paste(
      "Boxes are quartiles with 1.5 IQR whiskers. Points in panel a are jittered at alpha 0.04.",
      "The grouping conditions on the target: the published rule selects pairs with a high",
      "\ntarget Bayes factor, and a high target BF mechanically forces W > 0 unless the decoy is also high.",
      "A null pair that reached a high BF by chance would sit in the same box.",
      "\nThe figure locates W by discovery status; it does not validate the calls."
    ),
    theme = theme(plot.caption = element_text(hjust = 0, size = 8,
                                             colour = "grey25"))
  )
ggsave(file.path(figure_directory, "w_boxplot_by_discovery.png"),
       figure, width = 9.5, height = 5.8, dpi = 200)
ggsave(file.path(figure_directory, "w_boxplot_by_discovery.pdf"),
       figure, width = 9.5, height = 5.8)
message("Wrote figures to ", figure_directory)
