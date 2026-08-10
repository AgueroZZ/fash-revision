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
    aligned_gene_factor_for_selected_gene = factor_scores
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

message(
  "Completed global PCA factor visualization for ",
  selected_gene_id,
  ": ", output_dir
)
