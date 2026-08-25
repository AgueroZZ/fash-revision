#!/usr/bin/env Rscript

# Propagate matched-donor expression residual correlation through eQTL OLS weights.

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

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
}

matrix_to_long <- function(correlation, selection, method) {
  labels <- sub("^time_", "", colnames(correlation))
  grid <- expand.grid(time_a = labels, time_b = labels,
                      stringsAsFactors = FALSE)
  grid$correlation <- as.vector(correlation)
  grid$selection <- selection
  grid$method <- method
  grid
}

summarize_unit_correlations <- function(unit_correlations,
                                        selection,
                                        method,
                                        n_bootstrap = 500L,
                                        seed = 20260810L) {
  n_units <- dim(unit_correlations)[1L]
  n_time <- dim(unit_correlations)[2L]
  if (length(dim(unit_correlations)) != 3L || n_units < 2L ||
      n_time < 2L || any(!is.finite(unit_correlations))) {
    stop("Invalid unit-level correlation array.")
  }
  mean_matrix <- apply(unit_correlations, c(2L, 3L), mean)
  diag(mean_matrix) <- 1
  set.seed(seed)
  bootstrap_matrices <- array(
    NA_real_,
    dim = c(n_time, n_time, n_bootstrap),
    dimnames = list(
      dimnames(unit_correlations)[[2L]],
      dimnames(unit_correlations)[[3L]],
      paste0("bootstrap_", seq_len(n_bootstrap))
    )
  )
  bootstrap_lags <- matrix(NA_real_, nrow = n_bootstrap,
                           ncol = n_time - 1L)
  for (bootstrap_index in seq_len(n_bootstrap)) {
    sampled_units <- sample.int(n_units, n_units, replace = TRUE)
    matrix_draw <- apply(
      unit_correlations[sampled_units, , , drop = FALSE],
      c(2L, 3L),
      mean
    )
    diag(matrix_draw) <- 1
    bootstrap_matrices[, , bootstrap_index] <- matrix_draw
    bootstrap_lags[bootstrap_index, ] <- lag_average_correlation(matrix_draw)
  }
  list(
    mean_matrix = mean_matrix,
    lower_matrix = apply(bootstrap_matrices, c(1L, 2L), stats::quantile,
                         probs = 0.025, names = FALSE),
    upper_matrix = apply(bootstrap_matrices, c(1L, 2L), stats::quantile,
                         probs = 0.975, names = FALSE),
    lag_summary = data.frame(
      selection = selection,
      method = method,
      lag = seq_len(ncol(bootstrap_lags)),
      mean = lag_average_correlation(mean_matrix),
      lower = apply(bootstrap_lags, 2L, stats::quantile,
                    probs = 0.025, names = FALSE),
      upper = apply(bootstrap_lags, 2L, stats::quantile,
                    probs = 0.975, names = FALSE),
      n_units = n_units,
      n_bootstrap = n_bootstrap,
      stringsAsFactors = FALSE
    ),
    bootstrap_matrices = bootstrap_matrices
  )
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "r4_correlated_errors",
  "real_data_correlation_helpers.R"
))

output_id <- get_arg("--output-id", "multigene_null_beta_covariance")
n_bootstrap <- as.integer(get_arg("--bootstrap-reps", "500"))
seed <- as.integer(get_arg("--seed", "20260810"))
if (!nzchar(output_id) || is.na(n_bootstrap) || n_bootstrap < 100L ||
    is.na(seed)) {
  stop("Invalid output ID, bootstrap count, or seed.")
}

output_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "internal", output_id
)
raw_input_path <- file.path(output_dir, "input", "selected_raw_data.rds")
analysis_path <- file.path(output_dir, "multigene_null_beta_covariance.rds")
if (!file.exists(raw_input_path) || !file.exists(analysis_path)) {
  stop("The multigene null-covariance output is incomplete.")
}
summary_dir <- file.path(output_dir, "summary")
figure_dir <- file.path(output_dir, "figures")

raw_input <- readRDS(raw_input_path)
analysis <- readRDS(analysis_path)
required_input_fields <- c(
  "unit_table", "selection_indices", "selection_counts", "time_grid",
  "unit_dosage", "time_inputs"
)
if (!all(required_input_fields %in% names(raw_input))) {
  stop("The saved raw input does not contain all required fields.")
}
time_grid <- raw_input$time_grid
n_time <- length(time_grid)
unit_table <- raw_input$unit_table
if (n_time != 16L || nrow(unit_table) < 2L ||
    !identical(colnames(raw_input$unit_dosage), unit_table$pair_key)) {
  stop("The saved multigene input has invalid time or unit dimensions.")
}

make_ols_weights <- function(input, genotype) {
  genotype_residual <- input$projection$residualizer %*% genotype
  denominator <- colSums(genotype_residual^2)
  tolerance <- 1e-12 * pmax(1, colSums(genotype^2))
  if (any(!is.finite(denominator)) || any(denominator <= tolerance)) {
    stop("At least one observed genotype has zero residual information.")
  }
  weights <- sweep(genotype_residual, 2L, denominator, "/")
  rownames(weights) <- input$donors
  colnames(weights) <- colnames(genotype)
  weights
}

message("Constructing observed eQTL OLS weights at each time point.")
weights_by_time <- lapply(raw_input$time_inputs, function(input) {
  genotype <- raw_input$unit_dosage[input$donors, , drop = FALSE]
  make_ols_weights(input, genotype)
})

selection_summaries <- list()
comparison_rows <- list()
matrix_rows <- list()
lag_rows <- list()
matched_counts <- matrix(NA_integer_, nrow = n_time, ncol = n_time,
                         dimnames = list(paste0("time_", time_grid),
                                         paste0("time_", time_grid)))
for (time_one in seq_len(n_time)) {
  for (time_two in seq_len(n_time)) {
    matched_counts[time_one, time_two] <- length(intersect(
      raw_input$time_inputs[[time_one]]$donors,
      raw_input$time_inputs[[time_two]]$donors
    ))
  }
}

for (selection_index in seq_along(raw_input$selection_indices)) {
  selection <- names(raw_input$selection_indices)[selection_index]
  selected <- raw_input$selection_indices[[selection]]
  unit_correlations <- array(
    NA_real_,
    dim = c(length(selected), n_time, n_time),
    dimnames = list(
      unit_table$pair_key[selected],
      paste0("time_", time_grid),
      paste0("time_", time_grid)
    )
  )
  for (unit_position in seq_along(selected)) {
    unit_index <- selected[unit_position]
    covariance <- matrix(NA_real_, nrow = n_time, ncol = n_time,
                         dimnames = list(paste0("time_", time_grid),
                                         paste0("time_", time_grid)))
    for (time_one in seq_len(n_time)) {
      input_one <- raw_input$time_inputs[[time_one]]
      residual_one <- input_one$expression_residual[, unit_index]
      weight_one <- weights_by_time[[time_one]][, unit_index]
      residual_sd_one <- input_one$observed_residual_sd[unit_index]
      for (time_two in time_one:n_time) {
        input_two <- raw_input$time_inputs[[time_two]]
        residual_two <- input_two$expression_residual[, unit_index]
        weight_two <- weights_by_time[[time_two]][, unit_index]
        shared <- intersect(names(weight_one), names(weight_two))
        residual_correlation <- if (time_one == time_two) {
          1
        } else {
          stats::cor(residual_one[shared], residual_two[shared])
        }
        if (!is.finite(residual_correlation)) {
          stop("A matched-donor residual correlation is non-finite.")
        }
        residual_sd_two <- input_two$observed_residual_sd[unit_index]
        covariance_value <- residual_correlation * residual_sd_one *
          residual_sd_two * sum(weight_one[shared] * weight_two[shared])
        covariance[time_one, time_two] <- covariance[time_two, time_one] <-
          covariance_value
      }
    }
    if (any(!is.finite(covariance)) || any(diag(covariance) <= 0)) {
      stop("A propagated covariance matrix is invalid.")
    }
    correlation <- stats::cov2cor(covariance)
    correlation <- (correlation + t(correlation)) / 2
    diag(correlation) <- 1
    unit_correlations[unit_position, , ] <- correlation
  }

  summary <- summarize_unit_correlations(
    unit_correlations,
    selection = selection,
    method = "design_propagated_matched_donor_residual",
    n_bootstrap = n_bootstrap,
    seed = seed + selection_index
  )
  selection_summaries[[selection]] <- list(
    unit_correlations = unit_correlations,
    summary = summary
  )
  matrix_rows[[selection]] <- matrix_to_long(
    summary$mean_matrix,
    selection = selection,
    method = "design_propagated_matched_donor_residual"
  )
  lag_rows[[selection]] <- summary$lag_summary

  mc_key <- paste(
    selection,
    "donor_residual_block_permutation",
    "direct_mc_unit_mean",
    sep = "__"
  )
  if (!mc_key %in% names(analysis$mc_summaries)) {
    stop("The primary donor-residual-block Monte Carlo matrix is missing.")
  }
  permutation_matrix <- analysis$mc_summaries[[mc_key]]$mean_matrix
  propagated_matrix <- summary$mean_matrix
  off_diagonal <- row(propagated_matrix) != col(propagated_matrix)
  comparison_rows[[selection]] <- data.frame(
    selection = selection,
    matrix_correlation = stats::cor(
      propagated_matrix[off_diagonal], permutation_matrix[off_diagonal]
    ),
    mean_absolute_difference = mean(abs(
      propagated_matrix[off_diagonal] - permutation_matrix[off_diagonal]
    )),
    maximum_absolute_difference = max(abs(
      propagated_matrix[off_diagonal] - permutation_matrix[off_diagonal]
    )),
    propagated_minimum_eigenvalue = min(eigen(
      propagated_matrix,
      symmetric = TRUE,
      only.values = TRUE
    )$values),
    permutation_minimum_eigenvalue = min(eigen(
      permutation_matrix,
      symmetric = TRUE,
      only.values = TRUE
    )$values),
    stringsAsFactors = FALSE
  )
  matrix_rows[[paste0(selection, "_permutation")]] <- matrix_to_long(
    permutation_matrix,
    selection = selection,
    method = "donor_residual_block_permutation_direct_mc"
  )
  write_csv(
    data.frame(time = sub("^time_", "", rownames(permutation_matrix)),
               permutation_matrix, check.names = FALSE),
    file.path(
      summary_dir,
      paste0("primary_donor_residual_block_C_", selection, ".csv")
    )
  )
  write_csv(
    data.frame(time = sub("^time_", "", rownames(propagated_matrix)),
               propagated_matrix, check.names = FALSE),
    file.path(
      summary_dir,
      paste0("design_propagated_C_", selection, ".csv")
    )
  )
}

matrix_rows <- do.call(rbind, matrix_rows)
lag_rows <- do.call(rbind, lag_rows)
comparison <- do.call(rbind, comparison_rows)
rownames(matrix_rows) <- rownames(lag_rows) <- rownames(comparison) <- NULL

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required for figures.")
}
plot_data <- matrix_rows
plot_data$selection_label <- factor(
  plot_data$selection,
  levels = c("m_g_gt_0p900", "m_g_gt_0p925"),
  labels = c("m[g] > 0.90", "m[g] > 0.925")
)
plot_data$method_label <- factor(
  plot_data$method,
  levels = c(
    "design_propagated_matched_donor_residual",
    "donor_residual_block_permutation_direct_mc"
  ),
  labels = c("Design-propagated residual model", "Residual-block permutation")
)
plot_data$time_a <- factor(plot_data$time_a, levels = as.character(time_grid))
plot_data$time_b <- factor(plot_data$time_b,
                           levels = rev(as.character(time_grid)))
comparison_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = time_a, y = time_b, fill = correlation)
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-0.1, 0.8), oob = scales::squish, name = "Correlation"
  ) +
  ggplot2::facet_grid(selection_label ~ method_label) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    x = "Time point", y = "Time point",
    title = "Design propagation versus residual-block permutation",
    subtitle = "The propagated model assumes independent residuals across cell lines"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold")
  )
ggplot2::ggsave(
  file.path(figure_dir, "design_propagated_vs_permutation_heatmaps.png"),
  comparison_plot,
  width = 9.4,
  height = 6.4,
  dpi = 220
)

configuration <- list(
  output_id = output_id,
  created_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  model = paste(
    "For each gene-variant unit, matched-donor correlation of the expression",
    "residuals after the original intercept plus PC1--PC5 adjustment is",
    "propagated through the exact residualized-genotype OLS weights."
  ),
  cross_cell_line_assumption = paste(
    "Residuals from distinct cell lines are independent after conditioning",
    "on time-specific covariates; only same-cell-line covariance contributes."
  ),
  n_bootstrap = n_bootstrap,
  seed = seed
)
result <- list(
  configuration = configuration,
  matched_donor_counts = matched_counts,
  selection_summaries = selection_summaries,
  comparison = comparison
)
saveRDS(result, file.path(output_dir, "design_propagated_null_correlation.rds"),
        compress = "xz")
write_csv(as.data.frame(matched_counts),
          file.path(summary_dir, "matched_donor_counts.csv"))
write_csv(matrix_rows,
          file.path(summary_dir, "design_propagated_and_permutation_matrices_long.csv"))
write_csv(lag_rows,
          file.path(summary_dir, "design_propagated_lag_summaries.csv"))
write_csv(comparison,
          file.path(summary_dir, "design_propagated_vs_permutation_agreement.csv"))
message("Saved design-propagated correlation sensitivity to ", output_dir)
