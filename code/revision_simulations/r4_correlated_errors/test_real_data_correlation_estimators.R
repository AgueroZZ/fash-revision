#!/usr/bin/env Rscript

# Deterministic tests for the real-data full-correlation estimators.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/r4_correlated_errors/real_data_correlation_helpers.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/r4_correlated_errors/real_data_correlation_helpers.R")) {
    return("coderepo-local")
  }
  stop("Could not find the R4 correlation helpers.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "r4_correlated_errors",
  "real_data_correlation_helpers.R"
))

set.seed(91023)
n_units <- 20000L
n_time <- 6L
true_correlation <- matrix(0.20, nrow = n_time, ncol = n_time)
diag(true_correlation) <- 1
se <- matrix(
  stats::runif(n_units * n_time, min = 0.7, max = 1.3),
  nrow = n_units,
  ncol = n_time
)
constant_effect <- stats::rnorm(n_units)
errors <- matrix(
  stats::rnorm(n_units * n_time),
  nrow = n_units,
  ncol = n_time
) %*% chol(true_correlation)
beta_hat <- constant_effect + se * errors
pairwise_estimate <- estimate_pairwise_difference_correlation(beta_hat, se)
if (max(abs(pairwise_estimate - true_correlation)) > 0.035) {
  stop("The pairwise-difference estimator did not recover a known correlation matrix.")
}

direct_estimate <- estimate_direct_centered_correlation(beta_hat, se)
if (!isTRUE(all.equal(direct_estimate, t(direct_estimate), tolerance = 1e-12)) ||
    max(abs(diag(direct_estimate) - 1)) > 1e-12) {
  stop("The direct centered estimator is not a symmetric correlation matrix.")
}

non_pd <- true_correlation
non_pd[1, 2] <- non_pd[2, 1] <- 0.99
non_pd[1, 3] <- non_pd[3, 1] <- -0.99
projection <- project_to_positive_definite_correlation(non_pd)
if (min(eigen(
      projection$projected,
      symmetric = TRUE,
      only.values = TRUE
    )$values) <= 0 || max(abs(diag(projection$projected) - 1)) > 1e-10) {
  stop("Nearest-PD projection did not produce a positive-definite correlation matrix.")
}

raw_se <- matrix(
  stats::runif(200 * 6, min = 0.05, max = 0.3),
  nrow = 200,
  ncol = 6
)
beta_for_se <- matrix(stats::rnorm(200 * 6), nrow = 200) * raw_se
df <- c(12, 12, 9, 12, 9, 12)
t_value <- beta_for_se / raw_se
p_value <- 2 * stats::pt(abs(t_value), df = matrix(
  df,
  nrow = nrow(raw_se),
  ncol = ncol(raw_se),
  byrow = TRUE
), lower.tail = FALSE)
z_value <- stats::qnorm(1 - p_value / 2) * sign(t_value)
adjusted_se <- abs(beta_for_se) / abs(z_value)
recovered_raw_se <- invert_t_to_normal_se(beta_for_se, adjusted_se, df)
if (max(abs(recovered_raw_se - raw_se), na.rm = TRUE) > 1e-10) {
  stop("The inverse t-to-normal SE transformation is inaccurate.")
}

set.seed(12004)
n_candidates <- 6000L
selection_se <- matrix(
  stats::runif(n_candidates * n_time, min = 0.8, max = 1.2),
  nrow = n_candidates,
  ncol = n_time
)
independent_beta <- selection_se * matrix(
  stats::rnorm(n_candidates * n_time),
  nrow = n_candidates,
  ncol = n_time
)
unselected_pairwise <- estimate_pairwise_difference_correlation(
  independent_beta,
  selection_se
)
scores <- weighted_center_standardize(
  independent_beta,
  selection_se
)$residual_score
selected <- order(scores)[seq_len(500L)]
selected_pairwise <- estimate_pairwise_difference_correlation(
  independent_beta[selected, , drop = FALSE],
  selection_se[selected, , drop = FALSE]
)
if (abs(mean(unselected_pairwise[upper.tri(unselected_pairwise)])) > 0.05 ||
    mean(selected_pairwise[upper.tri(selected_pairwise)]) < 0.35) {
  stop("The flatness-selection test did not expose the expected pairwise-estimator bias.")
}

cat("Real-data correlation estimator tests passed.\n")
