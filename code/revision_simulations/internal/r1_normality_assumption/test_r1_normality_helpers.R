#!/usr/bin/env Rscript

# Focused tests for the internal R1 normality-assumption helpers.

helper_path <- file.path(
  "code", "revision_simulations", "internal",
  "r1_normality_assumption", "r1_normality_helpers.R"
)
if (!file.exists(helper_path)) {
  stop("The R1 normality helper file is missing: ", helper_path)
}
source(helper_path)

probability_grid <- (seq_len(200000L) - 0.5) / 200000L
gaussian_grid <- stats::qnorm(probability_grid)
t5_grid <- gaussian_to_standardized_t(gaussian_grid, df = 5)
stopifnot(
  all(is.finite(t5_grid)),
  all(sign(t5_grid) == sign(gaussian_grid)),
  identical(order(t5_grid), order(gaussian_grid)),
  abs(mean(t5_grid)) < 1e-12,
  abs(mean(t5_grid^2) - 1) < 0.01,
  identical(
    t5_grid,
    gaussian_to_standardized_t(gaussian_grid, df = 5)
  )
)

invalid_df_failed <- tryCatch(
  {
    gaussian_to_standardized_t(gaussian_grid, df = 2)
    FALSE
  },
  error = function(error) TRUE
)
stopifnot(invalid_df_failed)

G <- matrix(
  c(0, 1, 2, 0, 2, 1, 0, 2),
  nrow = 4L,
  ncol = 2L,
  dimnames = list(paste0("donor_", 1:4), paste0("unit_", 1:2))
)
covariates <- matrix(
  c(-1.5, -0.5, 0.5, 1.5),
  nrow = 4L,
  ncol = 1L,
  dimnames = list(rownames(G), "PC1")
)
beta_matrix <- matrix(
  c(0.2, -0.1, 0.4, 0.3),
  nrow = 2L,
  ncol = 2L,
  byrow = TRUE
)
covariate_effects <- array(
  c(0.5, -0.25, 0.2, 0.4),
  dim = c(1L, 2L, 2L)
)
intercepts <- matrix(
  c(1, 2, -1, 0.5),
  nrow = 2L,
  ncol = 2L
)

reconstructed_mean <- reconstruct_r1_expression_mean(
  G = G,
  beta_matrix = beta_matrix,
  covariates = covariates,
  covariate_effects = covariate_effects,
  intercepts = intercepts
)
manual_mean <- array(NA_real_, dim = c(4L, 2L, 2L))
for (time_index in 1:2) {
  manual_mean[, , time_index] <-
    sweep(G, 2L, beta_matrix[, time_index], `*`) +
    covariates %*% covariate_effects[, , time_index] +
    matrix(
      intercepts[, time_index],
      nrow = nrow(G),
      ncol = ncol(G),
      byrow = TRUE
    )
}
stopifnot(isTRUE(all.equal(
  unname(reconstructed_mean),
  unname(manual_mean),
  tolerance = 0
)))

oracle_se <- compute_oracle_regression_se(
  G = G,
  covariates = covariates,
  noise_sd = 2,
  n_time = 2L
)
manual_oracle <- vapply(seq_len(ncol(G)), function(unit_index) {
  design <- cbind(intercept = 1, G = G[, unit_index], covariates)
  2 * sqrt(solve(crossprod(design))["G", "G"])
}, numeric(1))
stopifnot(
  identical(dim(oracle_se), c(2L, 2L)),
  max(abs(oracle_se[, 1L] - manual_oracle)) < 1e-12,
  max(abs(oracle_se[, 2L] - manual_oracle)) < 1e-12
)

unit_info <- data.frame(
  unit_index = 1:4,
  unit_id = paste0("unit_", 1:4),
  effect_class = c("zero", "constant", "dynamic_bspline", "zero"),
  scenario = "test",
  stringsAsFactors = FALSE
)
true_beta <- matrix(
  c(
    0, 0,
    1, 1,
    -0.5, 0.5,
    0, 0
  ),
  nrow = 4L,
  ncol = 2L,
  byrow = TRUE
)
centered_error <- matrix(
  c(
    -0.5, 0.5,
    -1, 1,
    -1.5, 1.5,
    -2, 2
  ),
  nrow = 4L,
  ncol = 2L,
  byrow = TRUE
)
beta_hat <- true_beta + centered_error
se_uncorrected <- matrix(1, nrow = 4L, ncol = 2L)
se_adjusted <- matrix(1.25, nrow = 4L, ncol = 2L)
oracle_se_small <- matrix(0.8, nrow = 4L, ncol = 2L)
residual_df <- matrix(12, nrow = 4L, ncol = 2L)

error_table <- make_standardized_error_table(
  beta_hat = beta_hat,
  true_beta = true_beta,
  se_uncorrected = se_uncorrected,
  se_adjusted = se_adjusted,
  oracle_se = oracle_se_small,
  residual_df = residual_df,
  unit_info = unit_info,
  error_distribution = "gaussian"
)
expected_rows <- 3L * (
  sum(unit_info$effect_class == "zero") * ncol(beta_hat) +
    nrow(beta_hat) * ncol(beta_hat)
)
stopifnot(
  nrow(error_table) == expected_rows,
  setequal(error_table$population, c("zero_null", "all_centered")),
  setequal(
    error_table$se_scale,
    c("oracle_known_sigma", "raw_regression", "t_adjusted")
  ),
  all(error_table$residual_df == 12),
  all(is.finite(error_table$statistic)),
  all(error_table$se_ratio > 0)
)

distribution_summary <- summarize_standardized_errors(error_table)
stopifnot(
  nrow(distribution_summary) == 12L,
  all(c(
    "n", "mean", "sd", "median", "iqr", "skewness",
    "excess_kurtosis", "tail_rate_0_05", "tail_rate_0_01",
    "tail_rate_0_001"
  ) %in% names(distribution_summary)),
  all(is.finite(as.matrix(distribution_summary[, c(
    "mean", "sd", "median", "iqr", "skewness", "excess_kurtosis"
  )])))
)

qq_quantiles <- make_qq_quantiles(
  error_table,
  probabilities = c(0.1, 0.5, 0.9)
)
stopifnot(
  nrow(qq_quantiles) == 36L,
  setequal(qq_quantiles$reference_distribution, c("normal", "t")),
  all(qq_quantiles$reference_df[
    qq_quantiles$reference_distribution == "t"
  ] == 12),
  all(is.finite(qq_quantiles$empirical_quantile)),
  all(is.finite(qq_quantiles$reference_quantile))
)

uncertainty_summary <- summarize_se_uncertainty_bins(
  error_table,
  n_bins = 2L
)
stopifnot(
  nrow(uncertainty_summary) == 4L,
  setequal(uncertainty_summary$population, c("zero_null", "all_centered")),
  identical(
    as.integer(uncertainty_summary$se_ratio_bin),
    c(1L, 2L, 1L, 2L)
  ),
  all(uncertainty_summary$n[
    uncertainty_summary$population == "all_centered"
  ] == 4L),
  all(uncertainty_summary$n[
    uncertainty_summary$population == "zero_null"
  ] == 2L),
  all(is.finite(uncertainty_summary$sd)),
  all(uncertainty_summary$tail_rate_0_05 >= 0 &
        uncertainty_summary$tail_rate_0_05 <= 1)
)

cat("R1 normality helper tests passed.\n")
