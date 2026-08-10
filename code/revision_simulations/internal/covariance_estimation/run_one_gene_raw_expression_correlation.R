#!/usr/bin/env Rscript

# Inspect one gene's raw repeated-cell-line expression correlation matrix.

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

matrix_to_long <- function(matrix, value_name) {
  grid <- expand.grid(
    row_index = seq_len(nrow(matrix)),
    column_index = seq_len(ncol(matrix)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  result <- data.frame(
    time_row = sub("^time_", "", rownames(matrix)[grid$row_index]),
    time_column = sub("^time_", "", colnames(matrix)[grid$column_index]),
    stringsAsFactors = FALSE
  )
  result[[value_name]] <- as.vector(matrix)
  result
}

workflowr_root <- find_workflowr_root()
project_root <- normalizePath(file.path(workflowr_root, ".."), mustWork = TRUE)

threshold <- as.numeric(get_arg("--threshold", "0.925"))
requested_gene_id <- get_arg("--gene-id", "")
output_id <- get_arg("--output-id", "one_gene_raw_expression_correlation")
if (!is.finite(threshold) || threshold <= 0 || threshold >= 1 ||
    !nzchar(output_id)) {
  stop("Invalid threshold or output ID.")
}

fit_path <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
)
expression_path <- file.path(
  project_root, "iPSC-data", "expression-data",
  "quantile_normalized_no_projection.txt"
)
if (!file.exists(fit_path) || !file.exists(expression_path)) {
  stop("The FASH fit or raw expression matrix is missing.")
}

output_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "internal", output_id
)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading the FASH fit to select one null-enriched gene.")
fit_environment <- new.env(parent = emptyenv())
loaded_names <- load(fit_path, envir = fit_environment)
if (!identical(loaded_names, "fash_fit1_update")) {
  stop("The fitted-object file must contain only fash_fit1_update.")
}
fash_fit <- fit_environment$fash_fit1_update
pair_keys <- names(fash_fit$fash_data$data_list)
lfdr <- as.numeric(fash_fit$lfdr)
if (length(pair_keys) != length(lfdr) || any(!is.finite(lfdr)) ||
    any(lfdr < 0 | lfdr > 1)) {
  stop("The FASH pair keys or lfdr values are invalid.")
}
gene_ids <- sub("_.*$", "", pair_keys)
minimum_lfdr <- tapply(lfdr, gene_ids, min)
eligible_genes <- data.frame(
  gene_id = names(minimum_lfdr),
  min_lfdr = as.numeric(minimum_lfdr),
  stringsAsFactors = FALSE
)
eligible_genes <- eligible_genes[eligible_genes$min_lfdr > threshold, ]
eligible_genes <- eligible_genes[order(
  eligible_genes$min_lfdr,
  eligible_genes$gene_id,
  method = "radix"
), ]
rownames(eligible_genes) <- NULL
if (nrow(eligible_genes) < 1L) {
  stop("No genes satisfy the minimum-lfdr threshold.")
}

if (nzchar(requested_gene_id)) {
  selected_row <- match(requested_gene_id, eligible_genes$gene_id)
  if (is.na(selected_row)) {
    stop("The requested gene is not in the selected minimum-lfdr set.")
  }
  selected_gene <- eligible_genes[selected_row, , drop = FALSE]
  selection_rule <- "User-specified gene in the minimum-lfdr set"
} else {
  selected_row <- ceiling(nrow(eligible_genes) / 2L)
  selected_gene <- eligible_genes[selected_row, , drop = FALSE]
  selection_rule <- paste(
    "Deterministic median gene after sorting eligible genes by minimum lfdr",
    "and then gene ID"
  )
}

message("Reading raw expression for ", selected_gene$gene_id, ".")
expression_data <- utils::read.csv(
  expression_path,
  sep = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!identical(names(expression_data)[1L], "Gene_id") ||
    anyDuplicated(expression_data$Gene_id)) {
  stop("Raw expression gene identifiers are invalid.")
}
gene_row <- match(selected_gene$gene_id, expression_data$Gene_id)
if (is.na(gene_row)) {
  stop("The selected FASH gene is absent from raw expression.")
}

time_grid <- 0:15
sample_ids <- names(expression_data)[-1L]
sample_time <- suppressWarnings(as.integer(sub(
  "^.*_([0-9]+)$", "\\1", sample_ids
)))
cell_lines <- sort(unique(sub("_[0-9]+$", "", sample_ids)))
if (anyNA(sample_time) || any(!sample_time %in% time_grid) ||
    length(cell_lines) != 19L) {
  stop("The raw expression sample identifiers do not form the expected 19-cell-line panel.")
}

expression_matrix <- matrix(
  NA_real_,
  nrow = length(cell_lines),
  ncol = length(time_grid),
  dimnames = list(cell_lines, paste0("time_", time_grid))
)
for (time_index in seq_along(time_grid)) {
  time_value <- time_grid[time_index]
  expected_ids <- paste0(cell_lines, "_", time_value)
  columns <- match(expected_ids, names(expression_data))
  observed <- !is.na(columns)
  expression_matrix[observed, time_index] <- as.numeric(
    expression_data[gene_row, columns[observed], drop = TRUE]
  )
}
if (any(colSums(!is.na(expression_matrix)) < 2L) ||
    any(!is.finite(expression_matrix[!is.na(expression_matrix)]))) {
  stop("The selected gene cannot support pairwise expression correlations.")
}

pairwise_n_cell_lines <- crossprod(!is.na(expression_matrix))
storage.mode(pairwise_n_cell_lines) <- "integer"
dimnames(pairwise_n_cell_lines) <- list(
  colnames(expression_matrix),
  colnames(expression_matrix)
)
correlation_matrix <- stats::cor(
  expression_matrix,
  use = "pairwise.complete.obs",
  method = "pearson"
)
if (any(!is.finite(correlation_matrix)) ||
    max(abs(correlation_matrix - t(correlation_matrix))) > 1e-12 ||
    max(abs(diag(correlation_matrix) - 1)) > 1e-12) {
  stop("The selected gene has an invalid pairwise-complete correlation matrix.")
}

observed_index <- which(!is.na(expression_matrix), arr.ind = TRUE)
linear_time_data <- data.frame(
  expression = expression_matrix[observed_index],
  time = as.integer(sub(
    "^time_",
    "",
    colnames(expression_matrix)[observed_index[, 2L]]
  )),
  stringsAsFactors = FALSE
)
linear_time_fit <- stats::lm(expression ~ time, data = linear_time_data)
linear_time_residual_matrix <- matrix(
  NA_real_,
  nrow = nrow(expression_matrix),
  ncol = ncol(expression_matrix),
  dimnames = dimnames(expression_matrix)
)
linear_time_residual_matrix[observed_index] <- stats::residuals(linear_time_fit)
linear_time_correlation <- stats::cor(
  linear_time_residual_matrix,
  use = "pairwise.complete.obs",
  method = "pearson"
)
correlation_invariance_difference <- max(abs(
  correlation_matrix - linear_time_correlation
))
if (correlation_invariance_difference > 1e-10) {
  stop("Removing the pooled linear time effect unexpectedly changed correlation.")
}

time_mean <- colMeans(expression_matrix, na.rm = TRUE)
time_sd <- apply(expression_matrix, 2L, stats::sd, na.rm = TRUE)
standardized_expression <- sweep(
  sweep(expression_matrix, 2L, time_mean, "-"),
  2L,
  time_sd,
  "/"
)
time_summary <- data.frame(
  time = time_grid,
  n_cell_lines = colSums(!is.na(expression_matrix)),
  mean_expression = time_mean,
  sd_expression = time_sd,
  stringsAsFactors = FALSE
)

expression_long <- do.call(rbind, lapply(seq_along(time_grid), function(index) {
  data.frame(
    cell_line = rownames(expression_matrix),
    time = time_grid[index],
    raw_expression = expression_matrix[, index],
    standardized_expression = standardized_expression[, index],
    stringsAsFactors = FALSE
  )
}))
expression_long <- expression_long[!is.na(expression_long$raw_expression), ]
residual_long <- do.call(rbind, lapply(seq_along(time_grid), function(index) {
  data.frame(
    cell_line = rownames(linear_time_residual_matrix),
    time = time_grid[index],
    linear_time_residual = linear_time_residual_matrix[, index],
    stringsAsFactors = FALSE
  )
}))
residual_long <- residual_long[!is.na(residual_long$linear_time_residual), ]
cell_line_coverage <- data.frame(
  cell_line = rownames(expression_matrix),
  n_time_points_observed = rowSums(!is.na(expression_matrix)),
  stringsAsFactors = FALSE
)
cell_line_coverage$trajectory_status <- ifelse(
  cell_line_coverage$n_time_points_observed == length(time_grid),
  "Complete trajectory",
  "Partial trajectory"
)
expression_long$trajectory_status <- cell_line_coverage$trajectory_status[
  match(expression_long$cell_line, cell_line_coverage$cell_line)
]
expression_long$cell_line <- factor(
  expression_long$cell_line,
  levels = cell_line_coverage$cell_line
)
residual_long$trajectory_status <- cell_line_coverage$trajectory_status[
  match(residual_long$cell_line, cell_line_coverage$cell_line)
]
residual_long$cell_line <- factor(
  residual_long$cell_line,
  levels = cell_line_coverage$cell_line
)

utils::write.csv(
  selected_gene,
  file.path(output_dir, "selected_gene.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(cell_line = rownames(expression_matrix), expression_matrix,
             check.names = FALSE),
  file.path(output_dir, "cell_line_by_time_raw_expression.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(cell_line = rownames(standardized_expression),
             standardized_expression, check.names = FALSE),
  file.path(output_dir, "cell_line_by_time_standardized_expression.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(cell_line = rownames(linear_time_residual_matrix),
             linear_time_residual_matrix, check.names = FALSE),
  file.path(output_dir, "cell_line_by_time_linear_time_residual.csv"),
  row.names = FALSE
)
utils::write.csv(
  matrix_to_long(correlation_matrix, "correlation"),
  file.path(output_dir, "pairwise_correlation.csv"),
  row.names = FALSE
)
utils::write.csv(
  matrix_to_long(pairwise_n_cell_lines, "n_cell_lines"),
  file.path(output_dir, "pairwise_n_cell_lines.csv"),
  row.names = FALSE
)
utils::write.csv(
  time_summary,
  file.path(output_dir, "time_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  cell_line_coverage,
  file.path(output_dir, "cell_line_coverage.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(
    term = names(stats::coef(linear_time_fit)),
    estimate = as.numeric(stats::coef(linear_time_fit)),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "linear_time_fit_coefficients.csv"),
  row.names = FALSE
)

configuration <- list(
  created_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  fit_path = fit_path,
  expression_path = expression_path,
  threshold = threshold,
  selection_rule = selection_rule,
  n_eligible_genes = nrow(eligible_genes),
  selected_gene = selected_gene,
  n_cell_lines = length(cell_lines),
  time_grid = time_grid,
  correlation_definition = paste(
    "Pearson correlation across cell lines for each time-pair using",
    "pairwise-complete observed cell lines."
  ),
  linear_time_residualization = "Pooled expression ~ 1 + numeric time",
  max_abs_correlation_change_after_linear_time_residualization =
    correlation_invariance_difference
)
saveRDS(
  list(
    configuration = configuration,
    expression_matrix = expression_matrix,
    standardized_expression = standardized_expression,
    correlation_matrix = correlation_matrix,
    linear_time_residual_matrix = linear_time_residual_matrix,
    linear_time_correlation = linear_time_correlation,
    correlation_invariance_difference = correlation_invariance_difference,
    pairwise_n_cell_lines = pairwise_n_cell_lines,
    time_summary = time_summary
  ),
  file.path(output_dir, "one_gene_raw_expression_correlation.rds")
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  trajectory_plot <- ggplot2::ggplot(
    expression_long,
    ggplot2::aes(x = time, y = raw_expression, group = cell_line)
  ) +
    ggplot2::geom_line(colour = "grey55", alpha = 0.60, linewidth = 0.45) +
    ggplot2::geom_point(colour = "grey40", alpha = 0.70, size = 1.1) +
    ggplot2::geom_line(
      data = time_summary,
      ggplot2::aes(x = time, y = mean_expression),
      inherit.aes = FALSE,
      colour = "#0072B2",
      linewidth = 1.05
    ) +
    ggplot2::geom_point(
      data = time_summary,
      ggplot2::aes(x = time, y = mean_expression),
      inherit.aes = FALSE,
      colour = "#0072B2",
      size = 2.1
    ) +
    ggplot2::scale_x_continuous(breaks = time_grid) +
    ggplot2::labs(
      title = paste0("Raw expression trajectories: ", selected_gene$gene_id),
      subtitle = "Grey: individual cell lines; blue: mean among observed cell lines",
      x = "Time",
      y = "Quantile-normalized, unprojected expression"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

  trajectory_grid_plot <- ggplot2::ggplot(
    expression_long,
    ggplot2::aes(x = time, y = raw_expression, colour = trajectory_status)
  ) +
    ggplot2::geom_line(linewidth = 0.65) +
    ggplot2::geom_point(size = 1.25) +
    ggplot2::facet_wrap(~ cell_line, ncol = 5L) +
    ggplot2::scale_x_continuous(breaks = c(0L, 5L, 10L, 15L)) +
    ggplot2::scale_colour_manual(
      values = c(
        "Complete trajectory" = "#0072B2",
        "Partial trajectory" = "#D55E00"
      ),
      name = NULL
    ) +
    ggplot2::labs(
      title = paste0(
        "Cell-line-specific expression trajectories: ",
        selected_gene$gene_id
      ),
      subtitle = "Each panel is one cell line; 13 complete trajectories and 6 partial trajectories",
      x = "Time",
      y = "Quantile-normalized, unprojected expression"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  residual_grid_plot <- ggplot2::ggplot(
    residual_long,
    ggplot2::aes(
      x = time,
      y = linear_time_residual,
      colour = trajectory_status
    )
  ) +
    ggplot2::geom_hline(
      yintercept = 0,
      colour = "grey65",
      linewidth = 0.35
    ) +
    ggplot2::geom_line(linewidth = 0.65) +
    ggplot2::geom_point(size = 1.25) +
    ggplot2::facet_wrap(~ cell_line, ncol = 5L) +
    ggplot2::scale_x_continuous(breaks = c(0L, 5L, 10L, 15L)) +
    ggplot2::scale_colour_manual(
      values = c(
        "Complete trajectory" = "#0072B2",
        "Partial trajectory" = "#D55E00"
      ),
      name = NULL
    ) +
    ggplot2::labs(
      title = paste0(
        "Residual trajectories after pooled expression ~ 1 + time: ",
        selected_gene$gene_id
      ),
      subtitle = "The correlation matrix is unchanged because Pearson correlation centers each time pair",
      x = "Time",
      y = "Residual expression"
    ) +
    ggplot2::theme_minimal(base_size = 10) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )

  correlation_long <- matrix_to_long(correlation_matrix, "correlation")
  correlation_long$time_row <- as.integer(correlation_long$time_row)
  correlation_long$time_column <- as.integer(correlation_long$time_column)
  correlation_plot <- ggplot2::ggplot(
    correlation_long,
    ggplot2::aes(x = time_column, y = time_row, fill = correlation)
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", correlation)),
      size = 2.4
    ) +
    ggplot2::scale_x_continuous(breaks = time_grid) +
    ggplot2::scale_y_reverse(breaks = time_grid) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-1, 1),
      name = "Pearson r"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = paste0("One-gene cross-time correlation: ", selected_gene$gene_id),
      subtitle = "Each entry correlates matched cell lines at the two time points",
      x = "Time (column)",
      y = "Time (row)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )

  n_long <- matrix_to_long(pairwise_n_cell_lines, "n_cell_lines")
  n_long$time_row <- as.integer(n_long$time_row)
  n_long$time_column <- as.integer(n_long$time_column)
  n_plot <- ggplot2::ggplot(
    n_long,
    ggplot2::aes(x = time_column, y = time_row, fill = n_cell_lines)
  ) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = n_cell_lines), size = 3.0) +
    ggplot2::scale_x_continuous(breaks = time_grid) +
    ggplot2::scale_y_reverse(breaks = time_grid) +
    ggplot2::scale_fill_gradient(
      low = "#FEE8C8",
      high = "#E34A33",
      limits = c(13, 19),
      name = "Matched lines"
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = "Matched cell-line count for each correlation entry",
      subtitle = "The raw data contain 19 cell lines overall, with missing time points",
      x = "Time (column)",
      y = "Time (row)"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )

  ggplot2::ggsave(
    file.path(figure_dir, "selected_gene_raw_expression_trajectories.png"),
    trajectory_plot,
    width = 9,
    height = 4.8,
    dpi = 180
  )
  ggplot2::ggsave(
    file.path(figure_dir, "selected_gene_cell_line_trajectory_grid.png"),
    trajectory_grid_plot,
    width = 13,
    height = 10,
    dpi = 180
  )
  ggplot2::ggsave(
    file.path(figure_dir, "selected_gene_linear_time_residual_trajectory_grid.png"),
    residual_grid_plot,
    width = 13,
    height = 10,
    dpi = 180
  )
  ggplot2::ggsave(
    file.path(figure_dir, "selected_gene_pairwise_correlation.png"),
    correlation_plot,
    width = 10.5,
    height = 9.5,
    dpi = 180
  )
  ggplot2::ggsave(
    file.path(figure_dir, "selected_gene_pairwise_n_cell_lines.png"),
    n_plot,
    width = 10.5,
    height = 9.5,
    dpi = 180
  )
}

message(
  "Completed one-gene raw-expression correlation for ",
  selected_gene$gene_id,
  ": ",
  output_dir
)
