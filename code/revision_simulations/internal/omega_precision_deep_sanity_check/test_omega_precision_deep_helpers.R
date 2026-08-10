#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/internal/omega_precision_deep_sanity_check/omega_precision_deep_helpers.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/internal/omega_precision_deep_sanity_check/omega_precision_deep_helpers.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "omega_precision_deep_sanity_check",
  "omega_precision_deep_helpers.R"
))

suppressPackageStartupMessages({
  library(fashr)
  library(Matrix)
})

expect_true <- function(value, message) {
  if (!isTRUE(value)) {
    stop(message)
  }
}

expect_error <- function(expression, pattern) {
  result <- try(force(expression), silent = TRUE)
  if (!inherits(result, "try-error") || !grepl(pattern, as.character(result))) {
    stop("Expected an error matching: ", pattern)
  }
}

x <- 0:7
grid <- c(0, 0.45, 0.9)
true_weights <- c(0.5, 0.3, 0.2)
num_basis <- 4L
betaprec <- 0.8
order <- 1L
pred_step <- 1
error_sd <- 0.65
correlated_covariance <- error_sd^2 * ar1_correlation(length(x), 0.65)
diagonal_covariance <- diag(diag(correlated_covariance))
precision <- solve(correlated_covariance)

prior_covariances <- make_component_prior_covariances(
  x = x,
  grid = grid,
  num_basis = num_basis,
  betaprec = betaprec,
  order = order,
  pred_step = pred_step
)
full_covariances <- make_component_marginal_covariances(
  prior_covariances,
  correlated_covariance
)
diagonal_covariances <- make_component_marginal_covariances(
  prior_covariances,
  diagonal_covariance
)
simulated <- simulate_marginal_mixture(
  seed = 20260807L,
  n_units = 30L,
  true_weights = true_weights,
  component_covariances = full_covariances
)
full_log_likelihood <- component_log_likelihood_matrix(
  simulated$response,
  full_covariances
)
diagonal_log_likelihood <- component_log_likelihood_matrix(
  simulated$response,
  diagonal_covariances
)

package_full <- matrix(NA_real_, nrow = 3L, ncol = length(grid))
package_diagonal <- matrix(NA_real_, nrow = 3L, ncol = length(grid))
for (unit in seq_len(nrow(package_full))) {
  data_i <- data.frame(y = simulated$response[unit, ], x = x, offset = 0)
  for (component_index in seq_along(grid)) {
    package_full[unit, component_index] <- fashr:::compute_L_gaussian_helper(
      data_i = data_i,
      Si = NULL,
      Omegai = precision,
      psd_iwp = grid[[component_index]],
      num_basis = num_basis,
      betaprec = betaprec,
      order = order,
      pred_step = pred_step
    )
    package_diagonal[unit, component_index] <- fashr:::compute_L_gaussian_helper(
      data_i = data_i,
      Si = rep(error_sd, length(x)),
      Omegai = NULL,
      psd_iwp = grid[[component_index]],
      num_basis = num_basis,
      betaprec = betaprec,
      order = order,
      pred_step = pred_step
    )
  }
}
expect_true(
  max(abs(package_full - full_log_likelihood[seq_len(3L), ])) < 1e-8,
  "Vectorized full-Omega likelihood does not match the package."
)
expect_true(
  max(abs(package_diagonal - diagonal_log_likelihood[seq_len(3L), ])) < 1e-8,
  "Vectorized diagonal likelihood does not match the package."
)

independent_covariance <- error_sd^2 * diag(length(x))
independent_full <- make_component_marginal_covariances(
  prior_covariances,
  independent_covariance
)
independent_diagonal <- make_component_marginal_covariances(
  prior_covariances,
  diag(diag(independent_covariance))
)
independent_simulated <- simulate_marginal_mixture(
  seed = 20260808L,
  n_units = 30L,
  true_weights = true_weights,
  component_covariances = independent_full
)
independent_log_full <- component_log_likelihood_matrix(
  independent_simulated$response,
  independent_full
)
independent_log_diagonal <- component_log_likelihood_matrix(
  independent_simulated$response,
  independent_diagonal
)
expect_true(
  max(abs(independent_log_full - independent_log_diagonal)) < 1e-12,
  "rho = 0 did not produce identical likelihoods."
)

raw_full <- fit_raw_from_likelihood(full_log_likelihood, grid, penalty = 1)
bf_full <- fit_bf_from_raw(raw_full, grid)
oracle_full <- fit_oracle_from_likelihood(full_log_likelihood, grid, true_weights)
for (fit in list(raw_full, bf_full, oracle_full)) {
  expect_true(
    all(is.finite(expand_prior_weights(fit, grid))),
    "A fit contains non-finite prior weights."
  )
  expect_true(
    all(is.finite(expand_posterior_weights(fit, grid))),
    "A fit contains non-finite posterior weights."
  )
}

fdr_rows <- evaluate_fdr_fit(
  fit = oracle_full,
  truth_component = simulated$component,
  alpha_grid = c(0.01, 0.05, 0.10),
  seed = 20260807L,
  n_units = nrow(simulated$response),
  observation_model = "Omega",
  correction = "Oracle"
)
expect_true(nrow(fdr_rows) == 3L, "FDR evaluation returned the wrong row count.")
expect_true(
  all(fdr_rows$fdp >= 0 & fdr_rows$fdp <= 1),
  "FDP values are outside [0, 1]."
)
for (alpha in c(0.01, 0.05, 0.10)) {
  helper_selection <- select_lfdr_prefix(oracle_full$lfdr, alpha)
  package_selection <- capture.output(
    package_result <- fashr::fdr_control(
      oracle_full,
      alpha = alpha,
      plot = FALSE,
      sort = FALSE
    )
  )
  package_selected <- if (is.character(package_result$significant_units)) {
    match(package_result$significant_units, names(oracle_full$lfdr))
  } else {
    as.integer(package_result$significant_units)
  }
  expect_true(
    identical(
      sort(as.integer(helper_selection$selected)),
      sort(as.integer(package_selected))
    ),
    "The helper discovery rule does not match fashr::fdr_control()."
  )
}

prior_rows <- prior_weight_rows(
  fit = raw_full,
  truth_component = simulated$component,
  true_weights = true_weights,
  grid = grid,
  seed = 20260807L,
  n_units = nrow(simulated$response),
  observation_model = "Omega",
  correction = "Raw"
)
prior_summary <- summarize_prior_accuracy(prior_rows)
expect_true(
  nrow(prior_summary$summary) == length(grid) + 1L,
  "Prior summary returned the wrong row count."
)

expect_error(
  validate_symmetric_positive_definite(matrix(c(1, 2, 0, 1), 2L)),
  "symmetric"
)
expect_error(
  validate_symmetric_positive_definite(matrix(c(1, 2, 2, 1), 2L)),
  "positive definite"
)
expect_error(
  make_component_prior_covariances(
    x = x,
    grid = c(0.45, 0, 0.9),
    num_basis = num_basis,
    betaprec = betaprec,
    order = order,
    pred_step = pred_step
  ),
  "start with zero"
)
expect_error(
  evaluate_fdr_fit(
    fit = oracle_full,
    truth_component = simulated$component[-1L],
    alpha_grid = 0.05,
    seed = 1L,
    n_units = 29L,
    observation_model = "Omega",
    correction = "Oracle"
  ),
  "align"
)
underflow_raw <- raw_full
underflow_raw$L_matrix <- underflow_raw$L_matrix - 1e6
expect_error(
  fit_bf_from_raw(underflow_raw, grid),
  "non-negative, finite"
)

cat("All omega_precision_deep_helpers tests passed.\n")
