#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
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
  "covariance_estimation",
  "mashr_screen_helpers.R"
))

set.seed(20260844)
n_null <- 500L
n_signal <- 100L
n_time <- 4L
adjusted_se <- matrix(
  stats::runif((n_null + n_signal) * n_time, 0.8, 1.2),
  nrow = n_null + n_signal,
  ncol = n_time
)
true_beta <- matrix(0, nrow = n_null + n_signal, ncol = n_time)
true_beta[(n_null + 1L):(n_null + n_signal), ] <- 3
beta_hat <- true_beta + matrix(
  stats::rnorm((n_null + n_signal) * n_time),
  nrow = n_null + n_signal,
  ncol = n_time
) * adjusted_se
pair_names <- sprintf("pair_%04d", seq_len(nrow(beta_hat)))
time_names <- paste0("time_", 0:(n_time - 1L))
rownames(beta_hat) <- rownames(adjusted_se) <- pair_names
colnames(beta_hat) <- colnames(adjusted_se) <- time_names

first_fit <- fit_mashr_pair_screen(
  beta_hat,
  adjusted_se,
  z_threshold = 2,
  seed = 20260845L,
  verbose = FALSE
)
second_fit <- fit_mashr_pair_screen(
  beta_hat,
  adjusted_se,
  z_threshold = 2,
  seed = 20260845L,
  verbose = FALSE
)

z <- beta_hat / adjusted_se
maximum_keep <- apply(abs(z), 1L, max) < 2
expected_null_correlation <- stats::cor(z[maximum_keep, , drop = FALSE])
null_lfdr <- first_fit$pair_lfdr[seq_len(n_null)]
signal_lfdr <- first_fit$pair_lfdr[(n_null + 1L):(n_null + n_signal)]

stopifnot(
  identical(length(first_fit$pair_lfdr), nrow(beta_hat)),
  identical(dim(first_fit$condition_lfdr), dim(beta_hat)),
  identical(names(first_fit$pair_lfdr), pair_names),
  max(abs(first_fit$null_correlation - expected_null_correlation)) < 1e-12,
  first_fit$diagnostics$null_correlation_maximum_difference < 1e-12,
  first_fit$diagnostics$minimum_null_correlation_eigenvalue > 0,
  is.finite(first_fit$fitted_pi0),
  first_fit$fitted_pi0 >= 0,
  first_fit$fitted_pi0 <= 1,
  "null" %in% colnames(first_fit$fit$posterior_weights),
  all(first_fit$pair_lfdr >= 0 & first_fit$pair_lfdr <= 1),
  all(first_fit$condition_lfdr >= 0 & first_fit$condition_lfdr <= 1),
  stats::median(signal_lfdr) < stats::median(null_lfdr),
  mean(signal_lfdr <= 0.05) > 0.8,
  isTRUE(all.equal(
    first_fit$null_correlation,
    second_fit$null_correlation,
    tolerance = 0
  )),
  isTRUE(all.equal(
    first_fit$pair_lfdr,
    second_fit$pair_lfdr,
    tolerance = 1e-12
  ))
)

cat("Mashr Screen 3 helper tests passed.\n")
