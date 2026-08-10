#!/usr/bin/env Rscript

# Validate the matched donor-permutation variogram/covariance summary cache.

workflowr_root <- if (
  file.exists("code/revision_simulations/shared/simulation_functions.R")
) {
  normalizePath(".", mustWork = TRUE)
} else if (file.exists(
  "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
)) {
  normalizePath("coderepo-local", mustWork = TRUE)
} else {
  stop("Could not find the workflowr repository root.")
}

cache_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "permuted_variogram_covariance_summary"
)
result_path <- file.path(
  cache_dir,
  "permuted_variogram_covariance_summary.rds"
)
if (!file.exists(result_path)) {
  stop("The matched permutation summary cache is missing.")
}
result <- readRDS(result_path)

required_names <- c(
  "configuration",
  "replicate_lag_summaries",
  "aggregate_lag_summaries",
  "aggregate_method_contrasts",
  "mean_covariance_matrices_long",
  "mean_matrix_diagnostics",
  "matrix_agreement",
  "independent_covariance_draws",
  "covariance_signal_draws"
)
stopifnot(
  all(required_names %in% names(result)),
  identical(
    result$configuration$analysis_id,
    "permuted_variogram_covariance_summary"
  ),
  result$configuration$n_replications == 100L,
  identical(result$configuration$primary_se_scale, "t_adjusted"),
  nrow(result$aggregate_lag_summaries) == 810L,
  nrow(result$aggregate_method_contrasts) == 180L,
  nrow(result$mean_covariance_matrices_long) == 4608L,
  all(is.finite(result$aggregate_lag_summaries$mean)),
  all(is.finite(result$aggregate_method_contrasts$mean)),
  all(is.finite(result$mean_covariance_matrices_long$covariance))
)

validate_draw_list <- function(draw_list) {
  stopifnot(length(draw_list) == 18L)
  for (draws in draw_list) {
    stopifnot(
      identical(dim(draws), c(16L, 16L, 100L)),
      all(is.finite(draws)),
      max(abs(draws - aperm(draws, c(2L, 1L, 3L)))) < 1e-10
    )
  }
}
validate_draw_list(result$independent_covariance_draws)
validate_draw_list(result$covariance_signal_draws)

primary_key <- "lfdr_0p95__t_adjusted__direct_zero_null"
independent_mean <- apply(
  result$independent_covariance_draws[[primary_key]],
  c(1L, 2L),
  mean
)
independent_off_diagonal <- independent_mean[upper.tri(independent_mean)]
stopifnot(
  abs(mean(independent_off_diagonal)) < 0.005,
  max(abs(independent_off_diagonal)) < 0.025
)

primary_matrix <- apply(
  result$covariance_signal_draws[[primary_key]],
  c(1L, 2L),
  mean
)
lag_average <- function(matrix, lag) {
  mean(matrix[cbind(
    seq_len(ncol(matrix) - lag),
    (lag + 1L):ncol(matrix)
  )])
}
primary_eigenvalues <- eigen(
  primary_matrix,
  symmetric = TRUE,
  only.values = TRUE
)$values
stopifnot(
  max(abs(primary_matrix - t(primary_matrix))) < 1e-10,
  mean(diag(primary_matrix)) > 0.9,
  mean(diag(primary_matrix)) < 1.1,
  mean(primary_matrix[upper.tri(primary_matrix)]) > 0.15,
  lag_average(primary_matrix, 1L) > 0.35,
  lag_average(primary_matrix, 1L) > lag_average(primary_matrix, 15L),
  min(primary_eigenvalues) > 0
)

primary_agreement <- result$matrix_agreement[
  result$matrix_agreement$threshold == 0.95 &
    result$matrix_agreement$se_scale == "t_adjusted",
]
pairwise_agreement <- primary_agreement[
  primary_agreement$estimator == "pairwise_ols",
]
centered_agreement <- primary_agreement[
  primary_agreement$estimator == "within_unit_centered",
]
stopifnot(
  pairwise_agreement$off_diagonal_matrix_correlation > 0.995,
  pairwise_agreement$off_diagonal_rmse < 0.02,
  centered_agreement$off_diagonal_matrix_correlation > 0.98,
  centered_agreement$off_diagonal_rmse < 0.05
)

figure_paths <- file.path(
  cache_dir,
  "figure",
  c(
    "matched_variogram_and_covariance_t_adjusted.png",
    "matched_empirical_covariance_matrices_t_adjusted.png",
    "direct_empirical_covariance_matrix_t_adjusted.png"
  )
)
stopifnot(all(file.exists(figure_paths)), all(file.info(figure_paths)$size > 0))

cat("Matched donor-permutation variogram/covariance summary checks passed.\n")
