#!/usr/bin/env Rscript

# Render internal plots for the raw-expression minimum-lfdr covariance analysis.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(args, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(args[equals_hit[1L]], nchar(equals_prefix) + 1L))
  }
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1L] == length(args)) {
    return(default)
  }
  args[hit[1L] + 1L]
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required to render the diagnostic plots.")
}

workflowr_root <- find_workflowr_root()
output_id <- get_arg("--output-id", "raw_expression_min_lfdr_covariance")
analysis_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "internal", output_id
)
analysis_path <- file.path(analysis_dir, "raw_expression_min_lfdr_covariance.rds")
summary_dir <- file.path(analysis_dir, "summary")
figure_dir <- file.path(analysis_dir, "figures")
if (!file.exists(analysis_path)) {
  stop("The completed raw-expression analysis is missing: ", analysis_path)
}
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

analysis <- readRDS(analysis_path)
selection_counts <- utils::read.csv(
  file.path(summary_dir, "selection_counts.csv"),
  stringsAsFactors = FALSE
)
lag_profiles <- utils::read.csv(
  file.path(summary_dir, "lag_profiles.csv"),
  stringsAsFactors = FALSE
)
bootstrap_intervals <- utils::read.csv(
  file.path(summary_dir, "bootstrap_intervals.csv"),
  stringsAsFactors = FALSE
)

heatmap_rows <- do.call(rbind, lapply(names(analysis$matrices), function(set_id) {
  matrix <- analysis$matrices[[set_id]]$raw_expression_matrices[[
    "correlation_of_mean_covariance"
  ]]
  grid <- expand.grid(
    time_row = seq_len(nrow(matrix)),
    time_column = seq_len(ncol(matrix)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  threshold <- selection_counts$threshold[
    match(set_id, selection_counts$set_id)
  ]
  n_genes <- selection_counts$n_genes_passing_min_lfdr[
    match(set_id, selection_counts$set_id)
  ]
  data.frame(
    set_id = set_id,
    cutoff_label = paste0("m_g > ", format(threshold, nsmall = 3L, trim = TRUE),
                          "   (n = ", n_genes, ")"),
    time_row = grid$time_row - 1L,
    time_column = grid$time_column - 1L,
    value = as.vector(matrix),
    stringsAsFactors = FALSE
  )
}))
heatmap_rows$cutoff_label <- factor(
  heatmap_rows$cutoff_label,
  levels = unique(heatmap_rows$cutoff_label)
)

heatmap_plot <- ggplot2::ggplot(
  heatmap_rows,
  ggplot2::aes(x = time_column, y = time_row, fill = value)
) +
  ggplot2::geom_tile() +
  ggplot2::facet_wrap(~ cutoff_label, nrow = 1L) +
  ggplot2::scale_y_reverse(breaks = 0:15) +
  ggplot2::scale_x_continuous(breaks = 0:15) +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-0.10, 0.75),
    oob = scales::squish,
    name = "Correlation"
  ) +
  ggplot2::coord_fixed() +
  ggplot2::labs(
    title = "Raw-expression covariance-derived C",
    subtitle = "Mean within-gene covariance across 13 complete repeated cell lines, then cov2cor()",
    x = "Time (column)",
    y = "Time (row)"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    plot.title = ggplot2::element_text(face = "bold")
  )

lag_data <- lag_profiles[
  lag_profiles$expression_representation == "raw_expression" &
    lag_profiles$matrix_type %in% c(
      "correlation_of_mean_covariance",
      "mean_gene_correlation"
    ) &
    grepl("^lag_", lag_profiles$statistic),
  ,
  drop = FALSE
]
lag_data$lag <- as.integer(sub("^lag_", "", lag_data$statistic))
lag_data$cutoff_label <- paste0(
  "m_g > ",
  format(lag_data$threshold, nsmall = 3L, trim = TRUE)
)
lag_data$matrix_label <- ifelse(
  lag_data$matrix_type == "correlation_of_mean_covariance",
  "cov2cor(mean raw covariance)",
  "mean gene-level correlation"
)

lag_intervals <- bootstrap_intervals[
  bootstrap_intervals$expression_representation == "raw_expression" &
    bootstrap_intervals$matrix_type == "correlation_of_mean_covariance" &
    grepl("^lag_", bootstrap_intervals$statistic),
  ,
  drop = FALSE
]
lag_intervals$lag <- as.integer(sub("^lag_", "", lag_intervals$statistic))
lag_intervals$cutoff_label <- paste0(
  "m_g > ",
  format(lag_intervals$threshold, nsmall = 3L, trim = TRUE)
)

lag_plot <- ggplot2::ggplot() +
  ggplot2::geom_ribbon(
    data = lag_intervals,
    ggplot2::aes(
      x = lag,
      ymin = bootstrap_lower_95,
      ymax = bootstrap_upper_95
    ),
    alpha = 0.16,
    colour = NA,
    fill = "#0072B2"
  ) +
  ggplot2::geom_line(
    data = lag_data,
    ggplot2::aes(
      x = lag,
      y = value,
      colour = matrix_label,
      linetype = matrix_label
    ),
    linewidth = 0.85
  ) +
  ggplot2::geom_point(
    data = subset(
      lag_data,
      matrix_type == "correlation_of_mean_covariance"
    ),
    ggplot2::aes(x = lag, y = value, colour = matrix_label),
    size = 1.6
  ) +
  ggplot2::geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.35) +
  ggplot2::facet_wrap(~ cutoff_label, nrow = 1L) +
  ggplot2::scale_x_continuous(breaks = c(1L, 5L, 10L, 15L)) +
  ggplot2::scale_colour_manual(
    values = c(
      "cov2cor(mean raw covariance)" = "#0072B2",
      "mean gene-level correlation" = "#D55E00"
    ),
    name = "Aggregation"
  ) +
  ggplot2::scale_linetype_manual(
    values = c(
      "cov2cor(mean raw covariance)" = "solid",
      "mean gene-level correlation" = "dashed"
    ),
    name = "Aggregation"
  ) +
  ggplot2::labs(
    title = "Lag decay in repeated-cell-line raw expression",
    subtitle = "Shaded bands are 95% gene-bootstrap intervals for covariance-derived C",
    x = "Time lag",
    y = "Estimated correlation"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  filename = file.path(figure_dir, "raw_expression_covariance_C_heatmaps.png"),
  plot = heatmap_plot,
  width = 10,
  height = 4.8,
  dpi = 180
)
ggplot2::ggsave(
  filename = file.path(figure_dir, "raw_expression_covariance_lag_profiles.png"),
  plot = lag_plot,
  width = 8.8,
  height = 5.2,
  dpi = 180
)

message("Saved raw-expression covariance figures to ", figure_dir)
