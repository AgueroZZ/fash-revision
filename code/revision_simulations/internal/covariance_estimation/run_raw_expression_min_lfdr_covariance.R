#!/usr/bin/env Rscript

# Estimate repeated-cell-line raw-expression covariance from gene-level
# minimum-lfdr null-enriched sets.

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

parse_numeric_list <- function(value, name) {
  parsed <- suppressWarnings(as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (length(parsed) == 0L || any(!is.finite(parsed))) {
    stop("Invalid ", name, ".")
  }
  parsed
}

threshold_id <- function(threshold) {
  sub("\\.", "p", formatC(threshold, format = "f", digits = 3L))
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

parse_gene_ids <- function(pair_keys) {
  pair_keys <- as.character(pair_keys)
  if (length(pair_keys) == 0L || any(!nzchar(pair_keys)) ||
      any(!grepl("_", pair_keys, fixed = TRUE))) {
    stop("Every FASH key must include a gene and a variant.")
  }
  sub("_.*$", "", pair_keys)
}

mean_matrix_array <- function(matrices) {
  if (length(dim(matrices)) != 3L || any(dim(matrices) < 1L)) {
    stop("Expected a non-empty three-dimensional matrix array.")
  }
  result <- apply(matrices, c(2L, 3L), mean)
  dimnames(result) <- dimnames(matrices)[2:3]
  result
}

lag_summary_rows <- function(matrix) {
  if (!is.matrix(matrix) || nrow(matrix) != ncol(matrix) ||
      any(!is.finite(matrix))) {
    stop("Expected a finite square matrix.")
  }
  n_time <- ncol(matrix)
  lags <- seq_len(n_time - 1L)
  lag_values <- vapply(lags, function(lag) {
    mean(matrix[cbind(seq_len(n_time - lag), (lag + 1L):n_time)])
  }, numeric(1L))
  data.frame(
    statistic = c("mean_off_diagonal", paste0("lag_", lags)),
    value = c(mean(matrix[row(matrix) != col(matrix)]), lag_values),
    stringsAsFactors = FALSE
  )
}

matrix_diagnostics <- function(matrix, set_id, threshold, representation,
                               matrix_type) {
  if (!is.matrix(matrix) || nrow(matrix) != ncol(matrix) ||
      any(!is.finite(matrix))) {
    stop("Expected a finite square matrix for diagnostics.")
  }
  matrix <- (matrix + t(matrix)) / 2
  eigenvalues <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  off_diagonal <- matrix[row(matrix) != col(matrix)]
  data.frame(
    set_id = set_id,
    threshold = threshold,
    expression_representation = representation,
    matrix_type = matrix_type,
    minimum_eigenvalue = min(eigenvalues),
    maximum_eigenvalue = max(eigenvalues),
    n_negative_eigenvalues = sum(eigenvalues < -1e-10),
    minimum_diagonal = min(diag(matrix)),
    maximum_diagonal = max(diag(matrix)),
    mean_off_diagonal = mean(off_diagonal),
    minimum_off_diagonal = min(off_diagonal),
    maximum_off_diagonal = max(off_diagonal),
    stringsAsFactors = FALSE
  )
}

matrix_to_long <- function(matrix, set_id, threshold, representation,
                           matrix_type) {
  if (!is.matrix(matrix) || nrow(matrix) != ncol(matrix)) {
    stop("Expected a square matrix.")
  }
  grid <- expand.grid(
    row_index = seq_len(nrow(matrix)),
    column_index = seq_len(ncol(matrix)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  data.frame(
    set_id = set_id,
    threshold = threshold,
    expression_representation = representation,
    matrix_type = matrix_type,
    time_row = sub("^time_", "", rownames(matrix)[grid$row_index]),
    time_column = sub("^time_", "", colnames(matrix)[grid$column_index]),
    value = as.vector(matrix),
    stringsAsFactors = FALSE
  )
}

make_raw_expression_array <- function(expression_data, gene_rows,
                                      gene_ids, complete_cell_lines,
                                      time_grid) {
  time_labels <- paste0("time_", time_grid)
  result <- array(
    NA_real_,
    dim = c(length(gene_rows), length(complete_cell_lines), length(time_grid)),
    dimnames = list(gene_ids, complete_cell_lines, time_labels)
  )
  for (time_index in seq_along(time_grid)) {
    sample_ids <- paste0(complete_cell_lines, "_", time_grid[time_index])
    columns <- match(sample_ids, names(expression_data))
    if (anyNA(columns)) {
      stop("A complete-panel expression sample is missing at time ",
           time_grid[time_index], ".")
    }
    values <- as.matrix(expression_data[gene_rows, columns, drop = FALSE])
    storage.mode(values) <- "double"
    if (any(!is.finite(values))) {
      stop("Raw expression contains a non-finite selected value at time ",
           time_grid[time_index], ".")
    }
    result[, , time_index] <- values
  }
  result
}

make_pc_residual_expression_array <- function(expression_data, pc_data,
                                              gene_rows, gene_ids,
                                              complete_cell_lines, time_grid) {
  time_labels <- paste0("time_", time_grid)
  result <- array(
    NA_real_,
    dim = c(length(gene_rows), length(complete_cell_lines), length(time_grid)),
    dimnames = list(gene_ids, complete_cell_lines, time_labels)
  )
  expression_sample_ids <- names(expression_data)[-1L]
  for (time_index in seq_along(time_grid)) {
    time_value <- time_grid[time_index]
    sample_ids <- grep(
      paste0("_", time_value, "$"),
      expression_sample_ids,
      value = TRUE
    )
    expression_columns <- match(sample_ids, names(expression_data))
    pc_rows <- match(sample_ids, pc_data$Sample_id)
    if (anyNA(expression_columns) || anyNA(pc_rows)) {
      stop("Expression or PC sample matching failed at time ", time_value, ".")
    }
    expression_matrix <- t(as.matrix(
      expression_data[gene_rows, expression_columns, drop = FALSE]
    ))
    storage.mode(expression_matrix) <- "double"
    pc_matrix <- as.matrix(
      pc_data[pc_rows, paste0("PC", 1:5), drop = FALSE]
    )
    storage.mode(pc_matrix) <- "double"
    if (any(!is.finite(expression_matrix)) || any(!is.finite(pc_matrix))) {
      stop("Non-finite expression or PC value at time ", time_value, ".")
    }
    design_matrix <- cbind(intercept = 1, pc_matrix)
    design_qr <- qr(design_matrix)
    if (design_qr$rank != ncol(design_matrix)) {
      stop("The intercept-plus-PC design is rank deficient at time ",
           time_value, ".")
    }
    orthonormal_basis <- qr.Q(design_qr, complete = FALSE)
    residual_matrix <- expression_matrix -
      orthonormal_basis %*% crossprod(orthonormal_basis, expression_matrix)
    complete_rows <- match(
      paste0(complete_cell_lines, "_", time_value),
      sample_ids
    )
    if (anyNA(complete_rows)) {
      stop("A complete-panel PC-residual sample is missing at time ",
           time_value, ".")
    }
    result[, , time_index] <- t(residual_matrix[complete_rows, , drop = FALSE])
  }
  result
}

extract_per_gene_matrices <- function(expression_array) {
  dimensions <- dim(expression_array)
  if (length(dimensions) != 3L || dimensions[1L] < 2L ||
      dimensions[2L] < 2L || dimensions[3L] < 2L ||
      any(!is.finite(expression_array))) {
    stop("Expression array has invalid dimensions or values.")
  }
  n_gene <- dimensions[1L]
  n_time <- dimensions[3L]
  time_labels <- dimnames(expression_array)[[3L]]
  covariance <- array(
    NA_real_,
    dim = c(n_gene, n_time, n_time),
    dimnames = list(dimnames(expression_array)[[1L]], time_labels, time_labels)
  )
  correlation <- covariance
  for (gene_index in seq_len(n_gene)) {
    gene_expression <- expression_array[gene_index, , , drop = TRUE]
    if (!is.matrix(gene_expression)) {
      stop("Could not recover a cell-line by time expression matrix.")
    }
    current_covariance <- stats::cov(gene_expression)
    if (any(!is.finite(current_covariance)) ||
        any(diag(current_covariance) <= 0)) {
      stop("A selected gene has invalid raw-expression variance.")
    }
    current_correlation <- stats::cov2cor(current_covariance)
    if (any(!is.finite(current_correlation)) ||
        max(abs(diag(current_correlation) - 1)) > 1e-12) {
      stop("A selected gene has invalid raw-expression correlation.")
    }
    covariance[gene_index, , ] <- current_covariance
    correlation[gene_index, , ] <- current_correlation
  }
  list(covariance = covariance, correlation = correlation)
}

gene_scale_summary <- function(per_gene_covariance, set_id, threshold,
                               representation) {
  n_time <- dim(per_gene_covariance)[2L]
  time_labels <- dimnames(per_gene_covariance)[[2L]]
  data.frame(
    set_id = set_id,
    threshold = threshold,
    expression_representation = representation,
    time = sub("^time_", "", time_labels),
    minimum_gene_sd = vapply(seq_len(n_time), function(index) {
      min(sqrt(per_gene_covariance[, index, index]))
    }, numeric(1L)),
    median_gene_sd = vapply(seq_len(n_time), function(index) {
      stats::median(sqrt(per_gene_covariance[, index, index]))
    }, numeric(1L)),
    mean_gene_sd = vapply(seq_len(n_time), function(index) {
      mean(sqrt(per_gene_covariance[, index, index]))
    }, numeric(1L)),
    maximum_gene_sd = vapply(seq_len(n_time), function(index) {
      max(sqrt(per_gene_covariance[, index, index]))
    }, numeric(1L)),
    stringsAsFactors = FALSE
  )
}

bootstrap_correlation_summaries <- function(per_gene_covariance,
                                            per_gene_correlation,
                                            n_bootstrap, seed) {
  if (!identical(dim(per_gene_covariance), dim(per_gene_correlation))) {
    stop("Per-gene covariance and correlation arrays have different dimensions.")
  }
  n_gene <- dim(per_gene_covariance)[1L]
  result <- vector("list", 2L * n_bootstrap)
  output_index <- 0L
  set.seed(seed)
  for (replication in seq_len(n_bootstrap)) {
    selected_genes <- sample.int(n_gene, n_gene, replace = TRUE)
    mean_covariance <- mean_matrix_array(
      per_gene_covariance[selected_genes, , , drop = FALSE]
    )
    covariance_derived_correlation <- stats::cov2cor(mean_covariance)
    mean_gene_correlation <- mean_matrix_array(
      per_gene_correlation[selected_genes, , , drop = FALSE]
    )
    for (matrix_type in c(
      "correlation_of_mean_covariance",
      "mean_gene_correlation"
    )) {
      matrix <- if (identical(matrix_type, "correlation_of_mean_covariance")) {
        covariance_derived_correlation
      } else {
        mean_gene_correlation
      }
      output_index <- output_index + 1L
      current <- lag_summary_rows(matrix)
      current$replicate <- replication
      current$matrix_type <- matrix_type
      result[[output_index]] <- current
    }
  }
  do.call(rbind, result)
}

summarize_bootstrap_intervals <- function(bootstrap_draws, observed_rows) {
  grouping_columns <- c(
    "set_id", "threshold", "expression_representation",
    "matrix_type", "statistic"
  )
  split_key <- do.call(
    interaction,
    c(bootstrap_draws[, grouping_columns, drop = FALSE], list(
      drop = TRUE,
      lex.order = TRUE
    ))
  )
  intervals <- do.call(rbind, lapply(split(bootstrap_draws, split_key), function(x) {
    data.frame(
      x[1L, grouping_columns, drop = FALSE],
      bootstrap_median = stats::median(x$value),
      bootstrap_lower_95 = as.numeric(stats::quantile(x$value, 0.025)),
      bootstrap_upper_95 = as.numeric(stats::quantile(x$value, 0.975)),
      stringsAsFactors = FALSE
    )
  }))
  observed_rows <- observed_rows[
    observed_rows$matrix_type != "mean_covariance",
    c(grouping_columns, "value"),
    drop = FALSE
  ]
  names(observed_rows)[names(observed_rows) == "value"] <- "observed_value"
  merge(observed_rows, intervals, by = grouping_columns, sort = FALSE)
}

analyze_expression_array <- function(expression_array, set_id, threshold,
                                     representation, n_bootstrap,
                                     bootstrap_seed) {
  per_gene <- extract_per_gene_matrices(expression_array)
  mean_covariance <- mean_matrix_array(per_gene$covariance)
  covariance_derived_correlation <- stats::cov2cor(mean_covariance)
  mean_gene_correlation <- mean_matrix_array(per_gene$correlation)
  if (max(abs(diag(covariance_derived_correlation) - 1)) > 1e-12 ||
      max(abs(diag(mean_gene_correlation) - 1)) > 1e-12) {
    stop("An average expression correlation does not have a unit diagonal.")
  }
  list(
    matrices = list(
      mean_covariance = mean_covariance,
      correlation_of_mean_covariance = covariance_derived_correlation,
      mean_gene_correlation = mean_gene_correlation
    ),
    diagnostics = do.call(rbind, lapply(
      c("mean_covariance", "correlation_of_mean_covariance",
        "mean_gene_correlation"),
      function(matrix_type) {
        matrix_diagnostics(
          matrix = switch(
            matrix_type,
            mean_covariance = mean_covariance,
            correlation_of_mean_covariance = covariance_derived_correlation,
            mean_gene_correlation = mean_gene_correlation
          ),
          set_id = set_id,
          threshold = threshold,
          representation = representation,
          matrix_type = matrix_type
        )
      }
    )),
    lag_profiles = do.call(rbind, lapply(
      c("mean_covariance", "correlation_of_mean_covariance",
        "mean_gene_correlation"),
      function(matrix_type) {
        current <- lag_summary_rows(switch(
          matrix_type,
          mean_covariance = mean_covariance,
          correlation_of_mean_covariance = covariance_derived_correlation,
          mean_gene_correlation = mean_gene_correlation
        ))
        current$set_id <- set_id
        current$threshold <- threshold
        current$expression_representation <- representation
        current$matrix_type <- matrix_type
        current[, c(
          "set_id", "threshold", "expression_representation",
          "matrix_type", "statistic", "value"
        )]
      }
    )),
    bootstrap_draws = bootstrap_correlation_summaries(
      per_gene$covariance,
      per_gene$correlation,
      n_bootstrap = n_bootstrap,
      seed = bootstrap_seed
    ),
    gene_scale_summary = gene_scale_summary(
      per_gene$covariance,
      set_id = set_id,
      threshold = threshold,
      representation = representation
    )
  )
}

workflowr_root <- find_workflowr_root()
project_root <- normalizePath(file.path(workflowr_root, ".."), mustWork = TRUE)

thresholds <- sort(unique(parse_numeric_list(
  get_arg("--thresholds", "0.90,0.925"),
  "thresholds"
)))
n_bootstrap <- as.integer(get_arg("--bootstrap-reps", "500"))
bootstrap_seed <- as.integer(get_arg("--bootstrap-seed", "20260829"))
output_id <- get_arg("--output-id", "raw_expression_min_lfdr_covariance")
if (length(thresholds) != 2L || any(thresholds <= 0 | thresholds >= 1) ||
    is.na(n_bootstrap) || n_bootstrap < 20L || is.na(bootstrap_seed) ||
    !nzchar(output_id)) {
  stop("Invalid thresholds, bootstrap settings, or output ID.")
}

fit_path <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
)
expression_path <- file.path(
  project_root, "iPSC-data", "expression-data",
  "quantile_normalized_no_projection.txt"
)
pc_path <- file.path(
  project_root, "iPSC-data", "pc-data", "principal_components_10.txt"
)
required_inputs <- c(fit_path, expression_path, pc_path)
if (any(!file.exists(required_inputs))) {
  stop("Required input is missing: ",
       paste(required_inputs[!file.exists(required_inputs)], collapse = ", "))
}

output_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "internal", output_id
)
summary_dir <- file.path(output_dir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading the BF-adjusted FASH fit.")
fit_environment <- new.env(parent = emptyenv())
loaded_names <- load(fit_path, envir = fit_environment)
if (!identical(loaded_names, "fash_fit1_update")) {
  stop("The fitted-object file must contain only fash_fit1_update.")
}
fash_fit <- fit_environment$fash_fit1_update
pair_keys <- names(fash_fit$fash_data$data_list)
lfdr <- as.numeric(fash_fit$lfdr)
if (length(pair_keys) != length(lfdr) || length(pair_keys) < 2L ||
    any(!is.finite(lfdr)) || any(lfdr < 0 | lfdr > 1) ||
    anyDuplicated(pair_keys)) {
  stop("The fitted FASH object is incomplete or invalid.")
}
gene_id <- parse_gene_ids(pair_keys)
min_lfdr <- tapply(lfdr, gene_id, min)
gene_summary <- data.frame(
  gene_id = names(min_lfdr),
  min_lfdr = as.numeric(min_lfdr),
  stringsAsFactors = FALSE
)
gene_summary <- gene_summary[order(gene_summary$min_lfdr, decreasing = TRUE), ]
rownames(gene_summary) <- NULL

message("Reading raw expression and time-specific PCs.")
expression_data <- utils::read.csv(
  expression_path,
  sep = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
pc_data <- utils::read.delim(
  pc_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!identical(names(expression_data)[1L], "Gene_id") ||
    anyDuplicated(expression_data$Gene_id) ||
    anyDuplicated(pc_data$Sample_id) ||
    !all(c("Sample_id", paste0("PC", 1:5)) %in% names(pc_data))) {
  stop("The expression or PC input has invalid identifiers.")
}
expression_sample_ids <- names(expression_data)[-1L]
if (length(expression_sample_ids) != 297L ||
    !setequal(expression_sample_ids, pc_data$Sample_id)) {
  stop("Expression and PC sample IDs do not match exactly.")
}
time_grid <- 0:15
sample_time <- suppressWarnings(as.integer(sub(
  "^.*_([0-9]+)$", "\\1", expression_sample_ids
)))
sample_cell_line <- sub("_[0-9]+$", "", expression_sample_ids)
if (anyNA(sample_time) || any(!sample_time %in% time_grid) ||
    any(!nzchar(sample_cell_line))) {
  stop("Could not parse raw-expression cell-line time sample IDs.")
}
time_sample_ids <- lapply(time_grid, function(time_value) {
  expression_sample_ids[sample_time == time_value]
})
names(time_sample_ids) <- paste0("time_", time_grid)
time_cell_lines <- lapply(time_sample_ids, function(sample_ids) {
  sub("_[0-9]+$", "", sample_ids)
})
if (any(vapply(time_cell_lines, function(cell_lines) {
  anyDuplicated(cell_lines) > 0L
}, logical(1L)))) {
  stop("At least one time point contains duplicate cell-line observations.")
}
complete_cell_lines <- sort(Reduce(intersect, time_cell_lines))
if (length(complete_cell_lines) < 2L) {
  stop("Fewer than two cell lines are observed at every time point.")
}
complete_sample_ids <- unlist(lapply(time_grid, function(time_value) {
  paste0(complete_cell_lines, "_", time_value)
}), use.names = FALSE)
if (!all(complete_sample_ids %in% expression_sample_ids)) {
  stop("The complete repeated-measurement panel could not be reconstructed.")
}

sample_coverage <- data.frame(
  time = time_grid,
  n_cell_lines_observed = vapply(time_cell_lines, length, integer(1L)),
  n_complete_panel_cell_lines = length(complete_cell_lines),
  stringsAsFactors = FALSE
)
pairwise_complete_counts <- vapply(seq_along(time_grid), function(row_index) {
  vapply(seq_along(time_grid), function(column_index) {
    length(intersect(
      time_cell_lines[[row_index]],
      time_cell_lines[[column_index]]
    ))
  }, integer(1L))
}, integer(length(time_grid)))
dimnames(pairwise_complete_counts) <- list(
  paste0("time_", time_grid), paste0("time_", time_grid)
)
pairwise_coverage <- matrix_to_long(
  pairwise_complete_counts,
  set_id = "sample_coverage",
  threshold = NA_real_,
  representation = "raw_expression",
  matrix_type = "n_pairwise_cell_lines"
)

selection_counts <- list()
matrix_diagnostics_rows <- list()
lag_profile_rows <- list()
bootstrap_rows <- list()
gene_scale_rows <- list()
matrix_rows <- list()
selected_gene_tables <- list()
analysis_by_threshold <- list()

for (threshold_index in seq_along(thresholds)) {
  threshold <- thresholds[threshold_index]
  set_id <- paste0("min_lfdr_gt_", threshold_id(threshold))
  selected_genes <- gene_summary[
    gene_summary$min_lfdr > threshold,
    ,
    drop = FALSE
  ]
  if (nrow(selected_genes) < 2L) {
    stop("Fewer than two genes satisfy the minimum-lfdr threshold.")
  }
  gene_rows <- match(selected_genes$gene_id, expression_data$Gene_id)
  if (anyNA(gene_rows)) {
    stop("At least one selected gene is missing from raw expression.")
  }
  message(
    "Calculating raw-expression covariance for ", nrow(selected_genes),
    " genes with m_g > ", threshold, "."
  )
  raw_expression <- make_raw_expression_array(
    expression_data = expression_data,
    gene_rows = gene_rows,
    gene_ids = selected_genes$gene_id,
    complete_cell_lines = complete_cell_lines,
    time_grid = time_grid
  )
  raw_analysis <- analyze_expression_array(
    expression_array = raw_expression,
    set_id = set_id,
    threshold = threshold,
    representation = "raw_expression",
    n_bootstrap = n_bootstrap,
    bootstrap_seed = bootstrap_seed + threshold_index - 1L
  )
  message("Calculating the PC-residualized expression control for ", set_id, ".")
  pc_residual_expression <- make_pc_residual_expression_array(
    expression_data = expression_data,
    pc_data = pc_data,
    gene_rows = gene_rows,
    gene_ids = selected_genes$gene_id,
    complete_cell_lines = complete_cell_lines,
    time_grid = time_grid
  )
  pc_residual_analysis <- analyze_expression_array(
    expression_array = pc_residual_expression,
    set_id = set_id,
    threshold = threshold,
    representation = "pc_residual_expression",
    n_bootstrap = n_bootstrap,
    bootstrap_seed = bootstrap_seed + 100L + threshold_index - 1L
  )

  current_analyses <- list(
    raw_expression = raw_analysis,
    pc_residual_expression = pc_residual_analysis
  )
  for (representation in names(current_analyses)) {
    current <- current_analyses[[representation]]
    matrix_diagnostics_rows[[length(matrix_diagnostics_rows) + 1L]] <-
      current$diagnostics
    lag_profile_rows[[length(lag_profile_rows) + 1L]] <- current$lag_profiles
    gene_scale_rows[[length(gene_scale_rows) + 1L]] <- current$gene_scale_summary
    current_bootstrap <- current$bootstrap_draws
    current_bootstrap$set_id <- set_id
    current_bootstrap$threshold <- threshold
    current_bootstrap$expression_representation <- representation
    current_bootstrap <- current_bootstrap[, c(
      "set_id", "threshold", "expression_representation", "matrix_type",
      "replicate", "statistic", "value"
    )]
    bootstrap_rows[[length(bootstrap_rows) + 1L]] <- current_bootstrap
    for (matrix_type in names(current$matrices)) {
      matrix_rows[[length(matrix_rows) + 1L]] <- matrix_to_long(
        current$matrices[[matrix_type]],
        set_id = set_id,
        threshold = threshold,
        representation = representation,
        matrix_type = matrix_type
      )
    }
  }
  selection_counts[[set_id]] <- data.frame(
    set_id = set_id,
    threshold = threshold,
    n_genes_passing_min_lfdr = nrow(selected_genes),
    n_genes_in_raw_expression = length(gene_rows),
    n_complete_panel_cell_lines = length(complete_cell_lines),
    minimum_selected_min_lfdr = min(selected_genes$min_lfdr),
    maximum_selected_min_lfdr = max(selected_genes$min_lfdr),
    stringsAsFactors = FALSE
  )
  selected_gene_tables[[set_id]] <- selected_genes
  analysis_by_threshold[[set_id]] <- list(
    raw_expression_matrices = raw_analysis$matrices,
    pc_residual_expression_matrices = pc_residual_analysis$matrices
  )
}

selected_gene_id_sets <- lapply(selected_gene_tables, function(x) x$gene_id)
if (!all(selected_gene_id_sets[[2L]] %in% selected_gene_id_sets[[1L]])) {
  stop("The stricter minimum-lfdr gene set is not nested in the broader set.")
}

selection_counts <- do.call(rbind, selection_counts)
matrix_diagnostics_table <- do.call(rbind, matrix_diagnostics_rows)
lag_profiles <- do.call(rbind, lag_profile_rows)
bootstrap_draws <- do.call(rbind, bootstrap_rows)
gene_scale_table <- do.call(rbind, gene_scale_rows)
matrix_table <- do.call(rbind, matrix_rows)
bootstrap_intervals <- summarize_bootstrap_intervals(
  bootstrap_draws = bootstrap_draws,
  observed_rows = lag_profiles
)

configuration <- list(
  analysis_id = output_id,
  created_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  fit_path = fit_path,
  fit_object = "fash_fit1_update",
  expression_path = expression_path,
  pc_path = pc_path,
  n_fash_pairs = length(pair_keys),
  n_fash_genes = nrow(gene_summary),
  thresholds = thresholds,
  threshold_rule = "Strictly greater than the gene-level minimum BF-adjusted FASH lfdr",
  primary_estimand = paste(
    "Equal-gene average of within-gene raw-expression covariance matrices",
    "across the complete 13-cell-line by 16-time-point panel; its cov2cor",
    "standardization is the primary candidate C matrix."
  ),
  secondary_estimand = paste(
    "Equal-gene average of per-gene raw-expression correlation matrices;",
    "this reduces the influence of gene-specific expression scale."
  ),
  pc_control = paste(
    "At each time point, expression was residualized on an intercept and",
    "PCs 1 through 5 before applying the same complete-panel covariance."
  ),
  n_complete_panel_cell_lines = length(complete_cell_lines),
  complete_cell_lines = complete_cell_lines,
  time_grid = time_grid,
  n_bootstrap = n_bootstrap,
  bootstrap_seed = bootstrap_seed,
  r_version = R.version.string
)
source_information <- data.frame(
  role = c("BF-adjusted_FASH_fit", "raw_expression", "time_specific_PCs"),
  path = required_inputs,
  size_bytes = unname(file.info(required_inputs)$size),
  mtime = format(file.info(required_inputs)$mtime, tz = "UTC", usetz = TRUE),
  md5 = unname(tools::md5sum(required_inputs)),
  stringsAsFactors = FALSE
)
analysis_result <- list(
  configuration = configuration,
  gene_min_lfdr = gene_summary,
  selection_counts = selection_counts,
  sample_coverage = sample_coverage,
  pairwise_cell_line_coverage = pairwise_complete_counts,
  matrix_diagnostics = matrix_diagnostics_table,
  lag_profiles = lag_profiles,
  bootstrap_intervals = bootstrap_intervals,
  matrices = analysis_by_threshold
)

saveRDS(configuration, file.path(output_dir, "configuration.rds"))
saveRDS(
  analysis_result,
  file.path(output_dir, "raw_expression_min_lfdr_covariance.rds")
)
write_csv(gene_summary, file.path(summary_dir, "gene_min_lfdr.csv"))
write_csv(selection_counts, file.path(summary_dir, "selection_counts.csv"))
write_csv(sample_coverage, file.path(summary_dir, "sample_coverage.csv"))
write_csv(pairwise_coverage, file.path(summary_dir, "pairwise_cell_line_coverage.csv"))
write_csv(source_information, file.path(summary_dir, "source_information.csv"))
write_csv(matrix_diagnostics_table, file.path(summary_dir, "matrix_diagnostics.csv"))
write_csv(lag_profiles, file.path(summary_dir, "lag_profiles.csv"))
write_csv(bootstrap_draws, file.path(summary_dir, "bootstrap_draws.csv"))
write_csv(bootstrap_intervals, file.path(summary_dir, "bootstrap_intervals.csv"))
write_csv(gene_scale_table, file.path(summary_dir, "gene_scale_summary.csv"))
write_csv(matrix_table, file.path(summary_dir, "matrices_long.csv"))
for (set_id in names(selected_gene_tables)) {
  write_csv(
    selected_gene_tables[[set_id]],
    file.path(summary_dir, paste0("selected_genes_", set_id, ".csv"))
  )
}

message("Completed raw-expression covariance analysis: ", output_dir)
