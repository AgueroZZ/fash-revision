#!/usr/bin/env Rscript

# Verify the global PCA convention of eQTL covariates and plot its factors.

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

build_global_pca <- function(expression_by_sample, transform_name, n_components) {
  centered <- sweep(
    expression_by_sample,
    2L,
    colMeans(expression_by_sample),
    "-"
  )
  if (identical(transform_name, "gene_centered_and_scaled")) {
    gene_sd <- apply(centered, 2L, stats::sd)
    if (any(!is.finite(gene_sd)) || any(gene_sd <= 0)) {
      stop("Gene-wise scaling is not defined for the expression matrix.")
    }
    transformed <- sweep(centered, 2L, gene_sd, "/")
  } else if (identical(transform_name, "gene_centered")) {
    transformed <- centered
  } else {
    stop("Unknown global PCA transformation.")
  }

  gram <- tcrossprod(transformed)
  eig <- eigen(gram, symmetric = TRUE)
  if (sum(eig$values > 0) < n_components) {
    stop("The global expression matrix has insufficient positive singular values.")
  }
  component_index <- seq_len(n_components)
  left_factor <- eig$vectors[, component_index, drop = FALSE]
  singular_values <- sqrt(pmax(eig$values[component_index], 0))
  gene_factor <- crossprod(transformed, left_factor)
  colnames(left_factor) <- paste0("PC", component_index)
  colnames(gene_factor) <- paste0("PC", component_index)

  list(
    transform_name = transform_name,
    left_factor = left_factor,
    gene_factor = gene_factor,
    singular_values = singular_values
  )
}

align_to_saved_pcs <- function(reconstruction, saved_left_factor) {
  reconstructed_left_factor <- reconstruction$left_factor
  signs <- sign(colSums(reconstructed_left_factor * saved_left_factor))
  signs[signs == 0] <- 1
  aligned_left_factor <- sweep(reconstructed_left_factor, 2L, signs, "*")
  aligned_gene_factor <- sweep(reconstruction$gene_factor, 2L, signs, "*")
  component_correlations <- vapply(
    seq_len(ncol(saved_left_factor)),
    function(index) stats::cor(
      saved_left_factor[, index],
      aligned_left_factor[, index]
    ),
    numeric(1L)
  )

  list(
    transform_name = reconstruction$transform_name,
    signs = signs,
    left_factor = aligned_left_factor,
    gene_factor = aligned_gene_factor,
    singular_values = reconstruction$singular_values,
    component_correlations = component_correlations,
    max_abs_difference = max(abs(saved_left_factor - aligned_left_factor))
  )
}

make_long_loadings <- function(left_factor, sample_info) {
  pieces <- lapply(seq_len(ncol(left_factor)), function(index) {
    data.frame(
      sample_id = sample_info$sample_id,
      cell_line = sample_info$cell_line,
      time = sample_info$time,
      pc = colnames(left_factor)[index],
      loading = left_factor[, index],
      stringsAsFactors = FALSE
    )
  })
  loading_long <- do.call(rbind, pieces)
  loading_long <- loading_long[order(
    loading_long$pc,
    loading_long$cell_line,
    loading_long$time
  ), ]
  trajectory_group <- interaction(
    loading_long$pc,
    loading_long$cell_line,
    drop = TRUE
  )
  loading_long$segment <- ave(
    loading_long$time,
    trajectory_group,
    FUN = function(time_values) cumsum(c(TRUE, diff(time_values) > 1L))
  )
  loading_long$line_group <- interaction(
    loading_long$pc,
    loading_long$cell_line,
    loading_long$segment,
    drop = TRUE
  )
  loading_long
}

make_cell_line_by_time_matrix <- function(sample_info, values) {
  if (nrow(sample_info) != length(values)) {
    stop("Sample metadata and values have incompatible lengths.")
  }
  cell_lines <- levels(sample_info$cell_line)
  time_grid <- sort(unique(sample_info$time))
  row_index <- match(as.character(sample_info$cell_line), cell_lines)
  column_index <- match(sample_info$time, time_grid)
  if (anyNA(row_index) || anyNA(column_index) || anyDuplicated(
    paste(row_index, column_index, sep = "_"))
  ) {
    stop("The sample metadata cannot be converted to a cell-line-by-time matrix.")
  }
  result <- matrix(
    NA_real_,
    nrow = length(cell_lines),
    ncol = length(time_grid),
    dimnames = list(cell_lines, paste0("time_", time_grid))
  )
  result[cbind(row_index, column_index)] <- values
  result
}

matrix_to_time_long <- function(input_matrix, value_name) {
  grid <- expand.grid(
    row_index = seq_len(nrow(input_matrix)),
    column_index = seq_len(ncol(input_matrix)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  result <- data.frame(
    time_row = as.integer(sub(
      "^time_", "", rownames(input_matrix)[grid$row_index]
    )),
    time_column = as.integer(sub(
      "^time_", "", colnames(input_matrix)[grid$column_index]
    )),
    stringsAsFactors = FALSE
  )
  result[[value_name]] <- input_matrix[cbind(
    grid$row_index,
    grid$column_index
  )]
  result
}

summarize_correlation_lags <- function(correlation_matrix, matrix_label) {
  time_values <- as.integer(sub("^time_", "", colnames(correlation_matrix)))
  do.call(rbind, lapply(0:max(time_values), function(lag_value) {
    selected <- which(
      upper.tri(correlation_matrix, diag = TRUE) &
        abs(outer(time_values, time_values, "-")) == lag_value,
      arr.ind = TRUE
    )
    correlations <- correlation_matrix[selected]
    data.frame(
      matrix = matrix_label,
      lag = lag_value,
      n_time_pairs = length(correlations),
      mean_correlation = mean(correlations),
      median_correlation = stats::median(correlations),
      min_correlation = min(correlations),
      max_correlation = max(correlations),
      stringsAsFactors = FALSE
    )
  }))
}

workflowr_root <- find_workflowr_root()
project_root <- normalizePath(file.path(workflowr_root, ".."), mustWork = TRUE)

selected_gene_id <- get_arg("--gene-id", "ENSG00000067704")
output_id <- get_arg("--output-id", "global_pca_factor_visualization")
n_components <- as.integer(get_arg("--n-components", "10"))
if (!nzchar(selected_gene_id) || !nzchar(output_id) ||
    !is.finite(n_components) || n_components < 5L) {
  stop("Invalid gene ID, output ID, or number of components.")
}

pca_expression_path_arg <- get_arg("--pca-expression-path", "")
pca_expression_candidates <- if (nzchar(pca_expression_path_arg)) {
  pca_expression_path_arg
} else {
  c(
    file.path(
      project_root,
      "iPSC-data",
      "expression-data",
      "log_quantile_normalized_no_projection.txt"
    ),
    file.path(
      dirname(dirname(project_root)),
      "eQTL_analysis",
      "data",
      "log_quantile_normalized_no_projection.txt"
    )
  )
}
pca_expression_path <- pca_expression_candidates[
  file.exists(pca_expression_candidates)
][1L]
pc_path <- file.path(
  project_root,
  "iPSC-data",
  "pc-data",
  "principal_components_10.txt"
)
if (is.na(pca_expression_path) || !file.exists(pc_path)) {
  stop(paste(
    "The log quantile-normalized PCA input or saved PC covariate file is missing.",
    "Supply --pca-expression-path to specify the former."
  ))
}

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  output_id
)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading saved eQTL PC covariates.")
saved_pc_data <- utils::read.table(
  pc_path,
  header = TRUE,
  sep = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
expected_pc_names <- paste0("PC", seq_len(n_components))
if (!identical(names(saved_pc_data)[1L], "Sample_id") ||
    !all(expected_pc_names %in% names(saved_pc_data)) ||
    anyDuplicated(saved_pc_data$Sample_id)) {
  stop("The saved PC covariate file is malformed.")
}
saved_left_factor <- as.matrix(saved_pc_data[, expected_pc_names, drop = FALSE])
storage.mode(saved_left_factor) <- "double"
if (any(!is.finite(saved_left_factor))) {
  stop("The saved PC covariate file contains non-finite values.")
}

message("Loading the global sample-by-gene PCA input: ", pca_expression_path)
pca_expression_data <- utils::read.csv(
  pca_expression_path,
  sep = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!identical(names(pca_expression_data)[1L], "Gene_id") ||
    anyDuplicated(pca_expression_data$Gene_id)) {
  stop("PCA-input gene identifiers are invalid.")
}
sample_ids <- names(pca_expression_data)[-1L]
if (!identical(saved_pc_data$Sample_id, sample_ids)) {
  stop("Saved PC rows do not exactly match PCA-input sample columns.")
}
selected_gene_index <- match(selected_gene_id, pca_expression_data$Gene_id)
if (is.na(selected_gene_index)) {
  stop("The selected gene is absent from the PCA input.")
}
expression_by_sample <- t(as.matrix(pca_expression_data[, -1L, drop = FALSE]))
storage.mode(expression_by_sample) <- "double"
if (any(!is.finite(expression_by_sample))) {
  stop("The PCA input contains non-finite values.")
}

sample_info <- data.frame(
  sample_id = sample_ids,
  cell_line = sub("_[0-9]+$", "", sample_ids),
  time = suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", sample_ids))),
  stringsAsFactors = FALSE
)
if (anyNA(sample_info$time) || any(sample_info$time < 0L) ||
    any(sample_info$time > 15L) || length(unique(sample_info$cell_line)) != 19L) {
  stop("Sample identifiers do not form the expected 19-cell-line, 16-time-point panel.")
}
sample_info$cell_line <- factor(
  sample_info$cell_line,
  levels = sort(unique(sample_info$cell_line))
)

message("Reconstructing candidate global PCA factorizations.")
candidate_names <- c("gene_centered", "gene_centered_and_scaled")
candidate_results <- lapply(candidate_names, function(transform_name) {
  alignment <- align_to_saved_pcs(
    build_global_pca(expression_by_sample, transform_name, n_components),
    saved_left_factor
  )
  data.frame(
    transform = alignment$transform_name,
    pca_expression_path = pca_expression_path,
    min_component_correlation = min(alignment$component_correlations),
    mean_component_correlation = mean(alignment$component_correlations),
    max_abs_difference = alignment$max_abs_difference,
    stringsAsFactors = FALSE
  )
})
validation_table <- do.call(rbind, candidate_results)
utils::write.csv(
  validation_table,
  file.path(output_dir, "pca_reconstruction_validation.csv"),
  row.names = FALSE
)

match_index <- which(
  validation_table$min_component_correlation >= 0.999999 &
    validation_table$max_abs_difference <= 1e-5
)
if (length(match_index) != 1L) {
  stop(paste(
    "No unique candidate global PCA reconstruction matches the saved PC covariates.",
    "See pca_reconstruction_validation.csv before interpreting factor scores."
  ))
}
matched_transform <- validation_table$transform[match_index]
matched_alignment <- align_to_saved_pcs(
  build_global_pca(expression_by_sample, matched_transform, n_components),
  saved_left_factor
)

message("The saved PC covariates exactly match: ", matched_transform, ".")
loading_long <- make_long_loadings(
  saved_left_factor[, seq_len(5L), drop = FALSE],
  sample_info
)
factor_scores <- data.frame(
  gene_id = selected_gene_id,
  pc = expected_pc_names[seq_len(5L)],
  gene_factor_score = matched_alignment$gene_factor[
    selected_gene_index,
    seq_len(5L)
  ],
  stringsAsFactors = FALSE
)
factor_scores$direction <- ifelse(
  factor_scores$gene_factor_score >= 0,
  "Positive",
  "Negative"
)
selected_gene_standardized_expression <- as.numeric(
  scale(expression_by_sample[, selected_gene_index], center = TRUE, scale = TRUE)
)
top_five_fitted_expression <- as.numeric(
  saved_left_factor[, seq_len(5L), drop = FALSE] %*%
    matched_alignment$gene_factor[selected_gene_index, seq_len(5L)]
)
top_five_residual_expression <-
  selected_gene_standardized_expression - top_five_fitted_expression
residual_trajectory <- data.frame(
  sample_id = sample_info$sample_id,
  cell_line = sample_info$cell_line,
  time = sample_info$time,
  standardized_log_expression = selected_gene_standardized_expression,
  top_five_pc_fitted_expression = top_five_fitted_expression,
  top_five_pc_residual = top_five_residual_expression,
  stringsAsFactors = FALSE
)
residual_trajectory <- residual_trajectory[order(
  residual_trajectory$cell_line,
  residual_trajectory$time
), ]
residual_trajectory$trajectory_segment <- ave(
  residual_trajectory$time,
  residual_trajectory$cell_line,
  FUN = function(time_values) cumsum(c(TRUE, diff(time_values) > 1L))
)
residual_trajectory$line_group <- interaction(
  residual_trajectory$cell_line,
  residual_trajectory$trajectory_segment,
  drop = TRUE
)
total_sum_squares <- sum(selected_gene_standardized_expression^2)
fitted_sum_squares <- sum(top_five_fitted_expression^2)
residual_sum_squares <- sum(top_five_residual_expression^2)
residual_summary <- data.frame(
  gene_id = selected_gene_id,
  n_sample_time_observations = nrow(residual_trajectory),
  total_sum_squares = total_sum_squares,
  top_five_pc_fitted_sum_squares = fitted_sum_squares,
  top_five_pc_residual_sum_squares = residual_sum_squares,
  fraction_sum_squares_explained_by_top_five =
    fitted_sum_squares / total_sum_squares,
  max_abs_inner_product_residual_with_top_five_pcs = max(abs(
    crossprod(
      saved_left_factor[, seq_len(5L), drop = FALSE],
      top_five_residual_expression
    )
  )),
  stringsAsFactors = FALSE
)
pca_scale_expression_matrix <- make_cell_line_by_time_matrix(
  sample_info,
  selected_gene_standardized_expression
)
top_five_residual_matrix <- make_cell_line_by_time_matrix(
  sample_info,
  top_five_residual_expression
)
pairwise_n_cell_lines <- crossprod(!is.na(top_five_residual_matrix))
storage.mode(pairwise_n_cell_lines) <- "integer"
dimnames(pairwise_n_cell_lines) <- list(
  colnames(top_five_residual_matrix),
  colnames(top_five_residual_matrix)
)
pca_scale_correlation <- stats::cor(
  pca_scale_expression_matrix,
  use = "pairwise.complete.obs",
  method = "pearson"
)
top_five_residual_correlation <- stats::cor(
  top_five_residual_matrix,
  use = "pairwise.complete.obs",
  method = "pearson"
)
if (any(!is.finite(pca_scale_correlation)) ||
    any(!is.finite(top_five_residual_correlation)) ||
    max(abs(diag(top_five_residual_correlation) - 1)) > 1e-12) {
  stop("The selected gene cannot support valid pairwise residual correlations.")
}
top_five_residual_correlation_change <-
  top_five_residual_correlation - pca_scale_correlation
correlation_lag_summary <- rbind(
  summarize_correlation_lags(
    pca_scale_correlation,
    "PCA-input expression"
  ),
  summarize_correlation_lags(
    top_five_residual_correlation,
    "Top-five-PC residual"
  )
)
correlation_summary <- data.frame(
  matrix = c("PCA-input expression", "Top-five-PC residual"),
  mean_off_diagonal_correlation = c(
    mean(pca_scale_correlation[upper.tri(pca_scale_correlation)]),
    mean(top_five_residual_correlation[upper.tri(top_five_residual_correlation)])
  ),
  median_off_diagonal_correlation = c(
    stats::median(pca_scale_correlation[upper.tri(pca_scale_correlation)]),
    stats::median(
      top_five_residual_correlation[upper.tri(top_five_residual_correlation)]
    )
  ),
  min_off_diagonal_correlation = c(
    min(pca_scale_correlation[upper.tri(pca_scale_correlation)]),
    min(
      top_five_residual_correlation[upper.tri(top_five_residual_correlation)]
    )
  ),
  max_off_diagonal_correlation = c(
    max(pca_scale_correlation[upper.tri(pca_scale_correlation)]),
    max(
      top_five_residual_correlation[upper.tri(top_five_residual_correlation)]
    )
  ),
  minimum_eigenvalue = c(
    min(eigen(pca_scale_correlation, symmetric = TRUE, only.values = TRUE)$values),
    min(
      eigen(
        top_five_residual_correlation,
        symmetric = TRUE,
        only.values = TRUE
      )$values
    )
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  sample_info,
  file.path(output_dir, "sample_info.csv"),
  row.names = FALSE
)
utils::write.csv(
  loading_long,
  file.path(output_dir, "pc1_to_pc5_loading_trajectories.csv"),
  row.names = FALSE
)
utils::write.csv(
  factor_scores,
  file.path(output_dir, "selected_gene_factor_scores.csv"),
  row.names = FALSE
)
utils::write.csv(
  residual_trajectory,
  file.path(output_dir, "selected_gene_top5_pc_residual_trajectories.csv"),
  row.names = FALSE
)
utils::write.csv(
  residual_summary,
  file.path(output_dir, "selected_gene_top5_pc_residual_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(
    cell_line = rownames(pca_scale_expression_matrix),
    pca_scale_expression_matrix,
    check.names = FALSE
  ),
  file.path(output_dir, "selected_gene_pca_scale_expression_matrix.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(
    cell_line = rownames(top_five_residual_matrix),
    top_five_residual_matrix,
    check.names = FALSE
  ),
  file.path(output_dir, "selected_gene_top5_pc_residual_matrix.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(
    time = rownames(pca_scale_correlation),
    pca_scale_correlation,
    check.names = FALSE
  ),
  file.path(output_dir, "selected_gene_pca_scale_pairwise_correlation.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(
    time = rownames(top_five_residual_correlation),
    top_five_residual_correlation,
    check.names = FALSE
  ),
  file.path(
    output_dir,
    "selected_gene_top5_pc_residual_pairwise_correlation.csv"
  ),
  row.names = FALSE
)
utils::write.csv(
  data.frame(
    time = rownames(top_five_residual_correlation_change),
    top_five_residual_correlation_change,
    check.names = FALSE
  ),
  file.path(
    output_dir,
    "selected_gene_top5_pc_residual_correlation_change.csv"
  ),
  row.names = FALSE
)
utils::write.csv(
  data.frame(
    time = rownames(pairwise_n_cell_lines),
    pairwise_n_cell_lines,
    check.names = FALSE
  ),
  file.path(output_dir, "selected_gene_pairwise_n_cell_lines.csv"),
  row.names = FALSE
)
utils::write.csv(
  correlation_lag_summary,
  file.path(output_dir, "selected_gene_correlation_lag_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  correlation_summary,
  file.path(output_dir, "selected_gene_correlation_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  data.frame(
    pc = expected_pc_names,
    sign_to_saved_pc = matched_alignment$signs,
    component_correlation = matched_alignment$component_correlations,
    singular_value = matched_alignment$singular_values,
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "matched_component_details.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    selected_gene_id = selected_gene_id,
    pca_convention = paste(
      "Global PCA of a gene-wise centered and scaled log quantile-normalized",
      "sample-by-gene matrix, with rows indexed by cell-line--time pairs;",
      "the saved PC columns are the unit-norm left factors."
    ),
    pca_expression_path = pca_expression_path,
    matched_transform = matched_transform,
    validation_table = validation_table,
    sample_info = sample_info,
    saved_left_factor = saved_left_factor,
    aligned_gene_factor_for_selected_gene = factor_scores,
    selected_gene_top_five_pc_residual_trajectory = residual_trajectory,
    selected_gene_top_five_pc_residual_summary = residual_summary,
    selected_gene_pca_scale_correlation = pca_scale_correlation,
    selected_gene_top_five_pc_residual_correlation =
      top_five_residual_correlation,
    selected_gene_top_five_pc_residual_correlation_change =
      top_five_residual_correlation_change,
    selected_gene_pairwise_n_cell_lines = pairwise_n_cell_lines,
    selected_gene_correlation_lag_summary = correlation_lag_summary,
    selected_gene_correlation_summary = correlation_summary
  ),
  file.path(output_dir, "global_pca_factor_visualization.rds")
)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("ggplot2 is required to render the requested factor visualizations.")
}

cell_line_colours <- grDevices::hcl.colors(
  length(levels(sample_info$cell_line)),
  palette = "Dynamic"
)
names(cell_line_colours) <- levels(sample_info$cell_line)
loading_plot <- ggplot2::ggplot(
  loading_long,
  ggplot2::aes(
    x = time,
    y = loading,
    colour = cell_line,
    group = line_group
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey75",
    linewidth = 0.35
  ) +
  ggplot2::geom_line(linewidth = 0.50, alpha = 0.82) +
  ggplot2::geom_point(size = 1.15, alpha = 0.90) +
  ggplot2::facet_wrap(~ pc, ncol = 1L) +
  ggplot2::scale_colour_manual(values = cell_line_colours, name = "Cell line") +
  ggplot2::scale_x_continuous(breaks = 0:15) +
  ggplot2::labs(
    title = "Global PCA left factors by cell line and time",
    subtitle = paste(
      "Each line is L[(cell line, time), k] from a global standardized-log PCA;",
      "breaks denote absent sample-time observations."
    ),
    x = "Time",
    y = "Unit-norm left-factor value"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    strip.text = ggplot2::element_text(face = "bold"),
    legend.position = "bottom",
    legend.text = ggplot2::element_text(size = 7),
    legend.key.height = grid::unit(0.30, "cm"),
    legend.key.width = grid::unit(0.55, "cm")
  ) +
  ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2L, byrow = TRUE))

factor_score_plot <- ggplot2::ggplot(
  factor_scores,
  ggplot2::aes(x = pc, y = gene_factor_score, fill = direction)
) +
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey55",
    linewidth = 0.45
  ) +
  ggplot2::geom_col(width = 0.68) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", gene_factor_score)),
    vjust = ifelse(factor_scores$gene_factor_score >= 0, -0.40, 1.30),
    size = 3.5
  ) +
  ggplot2::scale_fill_manual(
    values = c("Positive" = "#0072B2", "Negative" = "#D55E00"),
    guide = "none"
  ) +
  ggplot2::labs(
    title = paste0("Static global-PCA gene factor scores: ", selected_gene_id),
    subtitle = paste(
      "F[gene, k] in Y_standardized = L F^T, after sign alignment to the saved PCs;",
      "these are not time-specific scores."
    ),
    x = NULL,
    y = "Gene factor score"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    panel.grid.minor = ggplot2::element_blank()
  )

residual_trajectory_plot <- ggplot2::ggplot(
  residual_trajectory,
  ggplot2::aes(x = time, y = top_five_pc_residual, group = line_group)
) +
  ggplot2::geom_hline(
    yintercept = 0,
    colour = "grey65",
    linewidth = 0.35
  ) +
  ggplot2::geom_line(
    colour = "#0072B2",
    linewidth = 0.65
  ) +
  ggplot2::geom_point(
    colour = "#0072B2",
    size = 1.25
  ) +
  ggplot2::facet_wrap(~ cell_line, ncol = 5L) +
  ggplot2::scale_x_continuous(breaks = c(0L, 5L, 10L, 15L)) +
  ggplot2::labs(
    title = paste0(
      "Residual trajectories after top-five global PCs: ", selected_gene_id
    ),
    subtitle = sprintf(
      paste(
        "Each panel is one cell line; PC1--PC5 explain %.1f%% of this gene's",
        "global standardized log-expression sum of squares."
      ),
      100 * residual_summary$fraction_sum_squares_explained_by_top_five
    ),
    x = "Time",
    y = "Standardized log-expression residual"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    strip.text = ggplot2::element_text(face = "bold")
  )

correlation_comparison_long <- rbind(
  transform(
    matrix_to_time_long(pca_scale_correlation, "correlation"),
    matrix = "PCA-input expression"
  ),
  transform(
    matrix_to_time_long(top_five_residual_correlation, "correlation"),
    matrix = "Top-five-PC residual"
  )
)
correlation_comparison_long$matrix <- factor(
  correlation_comparison_long$matrix,
  levels = c("PCA-input expression", "Top-five-PC residual")
)
correlation_comparison_plot <- ggplot2::ggplot(
  correlation_comparison_long,
  ggplot2::aes(x = time_column, y = time_row, fill = correlation)
) +
  ggplot2::geom_tile() +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.2f", correlation)),
    size = 2.25
  ) +
  ggplot2::facet_wrap(~ matrix, ncol = 2L) +
  ggplot2::scale_x_continuous(breaks = 0:15) +
  ggplot2::scale_y_reverse(breaks = 0:15) +
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
    title = paste0(
      "Cross-time correlation after top-five global-PC removal: ",
      selected_gene_id
    ),
    subtitle = paste(
      "Each entry is Pearson correlation across pairwise matched cell lines;",
      "the same standardized-log PCA scale is used in both panels."
    ),
    x = "Time (column)",
    y = "Time (row)"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold"),
    strip.text = ggplot2::element_text(face = "bold")
  )

ggplot2::ggsave(
  file.path(figure_dir, "global_pca_cell_line_loadings_over_time.png"),
  loading_plot,
  width = 13,
  height = 15,
  dpi = 180
)
ggplot2::ggsave(
  file.path(figure_dir, "global_pca_selected_gene_factor_scores.png"),
  factor_score_plot,
  width = 8.5,
  height = 5.5,
  dpi = 180
)
ggplot2::ggsave(
  file.path(
    figure_dir,
    "selected_gene_top5_pc_residual_trajectory_grid.png"
  ),
  residual_trajectory_plot,
  width = 13,
  height = 10,
  dpi = 180
)
ggplot2::ggsave(
  file.path(
    figure_dir,
    "selected_gene_top5_pc_residual_correlation_comparison.png"
  ),
  correlation_comparison_plot,
  width = 15.5,
  height = 8.5,
  dpi = 180
)

message(
  "Completed global PCA factor visualization for ",
  selected_gene_id,
  ": ", output_dir
)
