#!/usr/bin/env Rscript

# Test the traditional standardized-residual calculation used by Figure 1.

helper_path <- file.path(
  "code", "revision_simulations", "internal",
  "r1_normality_assumption",
  "real_data_residual_normality_helpers.R"
)
if (file.exists(helper_path)) {
  source(helper_path)
}

set.seed(20260820)
n_donors <- 19L
test_data <- data.frame(
  genotype = rep(c(0, 1, 2), length.out = n_donors),
  PC1 = stats::rnorm(n_donors),
  PC2 = stats::rnorm(n_donors),
  PC3 = stats::rnorm(n_donors),
  PC4 = stats::rnorm(n_donors),
  PC5 = stats::rnorm(n_donors)
)
test_data$expression <- with(
  test_data,
  0.7 + 0.4 * genotype - 0.2 * PC1 + 0.1 * PC3 +
    stats::rnorm(n_donors, sd = 0.6)
)
reference_fit <- stats::lm(
  expression ~ genotype + PC1 + PC2 + PC3 + PC4 + PC5,
  data = test_data
)
test_design <- stats::model.matrix(reference_fit)

tested_fit <- fit_standardized_residuals(
  expression = test_data$expression,
  design = test_design
)

stopifnot(
  identical(tested_fit$n_observations, 19L),
  identical(tested_fit$n_parameters, 7L),
  identical(tested_fit$residual_df, 12L),
  max(abs(tested_fit$fitted - stats::fitted(reference_fit))) < 1e-10,
  max(abs(tested_fit$residual - stats::residuals(reference_fit))) < 1e-10,
  max(abs(tested_fit$leverage - stats::hatvalues(reference_fit))) < 1e-10,
  abs(
    tested_fit$residual_standard_error -
      summary(reference_fit)$sigma
  ) < 1e-10,
  max(abs(
    tested_fit$standardized_residual -
      stats::rstandard(reference_fit)
  )) < 1e-10
)

rank_deficient_design <- test_design
rank_deficient_design[, "PC5"] <- rank_deficient_design[, "PC4"]
stopifnot(is.null(fit_standardized_residuals(
  expression = test_data$expression,
  design = rank_deficient_design
)))

nonfinite_expression <- test_data$expression
nonfinite_expression[1L] <- NA_real_
stopifnot(is.null(fit_standardized_residuals(
  expression = nonfinite_expression,
  design = test_design
)))

cat("Real-data standardized-residual helper tests passed.\n")
