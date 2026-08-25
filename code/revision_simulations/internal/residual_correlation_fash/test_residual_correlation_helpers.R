#!/usr/bin/env Rscript

# Tests for the residual-based common-correlation helpers.
#
# The load-bearing test is `test_closed_form_matches_tmb`: the closed-form
# posterior mean must reproduce the optimum that `fashr` itself finds by TMB
# plus nlminb, for both the diagonal and the dependent likelihood, at psd = 0
# and at nonzero psd.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "residual_correlation_fash", "residual_correlation_helpers.R"
))

suppressPackageStartupMessages({
  library(fashr)
  library(Matrix)
})

expect_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop("FAILED: ", message)
  }
  invisible(TRUE)
}

expect_error <- function(expression, message) {
  outcome <- tryCatch({
    force(expression)
    FALSE
  }, error = function(error) TRUE)
  expect_true(outcome, message)
}

expect_close <- function(actual, expected, tolerance, message) {
  difference <- max(abs(as.numeric(actual) - as.numeric(expected)))
  expect_true(
    is.finite(difference) && difference <= tolerance,
    paste0(message, " (maximum absolute difference ", format(difference), ")")
  )
  difference
}

# Reference implementation: the package's own TMB optimum for one component.
fashr_component_mean <- function(y, x, precision, psd, num_basis, order,
                                 betaprec, pred_step) {
  namespace <- asNamespace("fashr")
  diagonal <- !is.matrix(precision)
  data_i <- list(y = y, x = x, offset = rep(0, length(y)))
  standard_errors <- if (diagonal) 1 / sqrt(precision) else NULL
  omega_matrix <- if (diagonal) NULL else precision
  tmbdat <- namespace$fash_set_tmbdat(
    data_i,
    Si = standard_errors,
    Omegai = omega_matrix,
    num_basis = num_basis,
    betaprec = betaprec,
    order = order
  )
  knots <- seq(min(x), max(x), length.out = num_basis)
  global_design <- as.matrix(namespace$global_poly_helper(x, p = order))
  if (psd != 0) {
    tmbdat$sigmaIWP <- psd / sqrt(
      (pred_step^((2 * order) - 1)) /
        (((2 * order) - 1) * (factorial(order - 1)^2))
    )
    spline_design <- as.matrix(
      namespace$local_poly_helper(knots = knots, refined_x = x, p = order)
    )
    design <- cbind(spline_design, global_design)
    dll <- if (diagonal) "Gaussian_ind" else "Gaussian_dep"
    n_parameter <- ncol(tmbdat$X) + ncol(tmbdat$B)
  } else {
    design <- global_design
    dll <- if (diagonal) "Gaussian_ind_fixed" else "Gaussian_dep_fixed"
    n_parameter <- ncol(tmbdat$X)
  }
  ff <- TMB::MakeADFun(
    data = tmbdat,
    parameters = list(W = rep(0, n_parameter)),
    DLL = dll,
    silent = TRUE
  )
  ff$he <- function(w) numDeriv::jacobian(ff$gr, w)
  opt <- stats::nlminb(
    start = ff$par, objective = ff$fn, gradient = ff$gr, hessian = ff$he,
    control = list(eval.max = 20000, iter.max = 20000)
  )
  as.numeric(design %*% opt$par)
}

set.seed(20260820)
n_time <- 16L
time_grid <- 0:15
num_basis <- 20L
pred_step <- 1
betaprec <- 0

test_closed_form_matches_tmb <- function() {
  worst <- 0
  for (order in c(1L, 2L)) {
    shared_design <- build_shared_design(
      time_grid = time_grid,
      num_basis = num_basis,
      order = order,
      betaprec = betaprec
    )
    psd_values <- c(0, 0.05, 0.2, 1)
    scales <- psd_to_prior_scale(psd_values, order, pred_step)$prior_scale
    y <- as.numeric(sin(time_grid / 4) + stats::rnorm(n_time, sd = 0.4))
    standard_errors <- stats::runif(n_time, 0.2, 0.7)

    diagonal_means <- unit_component_posterior_means(
      y = y,
      precision = 1 / standard_errors^2,
      shared_design = shared_design,
      prior_scale = scales
    )
    for (index in seq_along(psd_values)) {
      reference <- fashr_component_mean(
        y, time_grid, 1 / standard_errors^2, psd_values[index],
        num_basis, order, betaprec, pred_step
      )
      worst <- max(worst, expect_close(
        diagonal_means[, index], reference, 1e-7,
        paste0("diagonal closed form matches TMB, order ", order,
               ", psd ", psd_values[index])
      ))
    }

    random <- matrix(stats::rnorm(n_time * n_time), n_time, n_time)
    correlation <- stats::cov2cor(
      crossprod(random) / n_time + diag(0.8, n_time)
    )
    omega <- diag(1 / standard_errors) %*% solve(correlation) %*%
      diag(1 / standard_errors)
    omega <- (omega + t(omega)) / 2
    dependent_means <- unit_component_posterior_means(
      y = y,
      precision = omega,
      shared_design = shared_design,
      prior_scale = scales
    )
    for (index in seq_along(psd_values)) {
      reference <- fashr_component_mean(
        y, time_grid, omega, psd_values[index],
        num_basis, order, betaprec, pred_step
      )
      worst <- max(worst, expect_close(
        dependent_means[, index], reference, 1e-7,
        paste0("dependent closed form matches TMB, order ", order,
               ", psd ", psd_values[index])
      ))
    }
  }
  cat("  closed form vs TMB, worst absolute difference:",
      format(worst, digits = 3), "\n")
}

test_bma_matches_predict_fash <- function() {
  order <- 1L
  n_unit <- 6L
  psd_grid <- sort(c(0, exp(-0.5 * seq(0, 10, by = 2))))
  data_list <- list()
  standard_errors <- matrix(NA_real_, n_unit, n_time)
  beta_hat <- matrix(NA_real_, n_unit, n_time)
  for (index in seq_len(n_unit)) {
    se <- stats::runif(n_time, 0.2, 0.6)
    signal <- if (index %% 2 == 0) sin(time_grid / 3) else rep(0.1, n_time)
    y <- signal + stats::rnorm(n_time, sd = se)
    data_list[[paste0("unit", index)]] <- data.frame(beta = y, time = time_grid,
                                                     SE = se)
    standard_errors[index, ] <- se
    beta_hat[index, ] <- y
  }
  rownames(beta_hat) <- names(data_list)
  rownames(standard_errors) <- names(data_list)

  fit <- fashr::fash(
    Y = "beta", smooth_var = "time", S = "SE", data_list = data_list,
    num_basis = num_basis, order = order, betaprec = betaprec,
    pred_step = pred_step, penalty = 1, grid = psd_grid,
    num_cores = 1, verbose = FALSE
  )
  shared_design <- build_shared_design(time_grid, num_basis, order, betaprec)
  retained <- retained_psd_values(fit)
  expect_true(
    length(retained) <= length(fit$psd_grid),
    "the retained PSD values are a subset of the grid"
  )
  scales <- psd_to_prior_scale(retained, order, pred_step)$prior_scale
  closed_form <- bma_posterior_means(
    beta_hat = beta_hat,
    standard_errors = standard_errors,
    posterior_weights = fit$posterior_weights,
    shared_design = shared_design,
    prior_scale = scales
  )
  set.seed(99)
  monte_carlo <- t(vapply(seq_len(n_unit), function(index) {
    stats::predict(fit, index = index, M = 8000)$mean
  }, numeric(n_time)))
  standard_error_of_mean <- max(
    apply(stats::predict(fit, index = 1, only.samples = TRUE, M = 8000), 1, sd)
  ) / sqrt(8000)
  difference <- max(abs(closed_form - monte_carlo))
  cat("  closed form vs predict.fash Monte Carlo, maximum difference:",
      format(difference, digits = 3),
      "(Monte Carlo standard error approximately",
      format(standard_error_of_mean, digits = 2), ")\n")
  expect_true(
    difference < 12 * standard_error_of_mean,
    "BMA posterior mean agrees with predict.fash within Monte Carlo error"
  )
}

test_null_component_is_weighted_mean <- function() {
  # At order 1 the psd = 0 component is a constant, and with betaprec = 0 that
  # constant is the precision-weighted mean. This is the source of the
  # -1/(R-1) floor on the residual correlation.
  shared_design <- build_shared_design(time_grid, num_basis, 1L, 0)
  scales <- psd_to_prior_scale(c(0, 0.5), 1L, pred_step)$prior_scale
  y <- stats::rnorm(n_time)
  se <- stats::runif(n_time, 0.2, 0.6)
  means <- unit_component_posterior_means(y, 1 / se^2, shared_design, scales)
  weighted_mean <- sum(y / se^2) / sum(1 / se^2)
  expect_close(means[, 1], rep(weighted_mean, n_time), 1e-10,
               "psd = 0 component equals the precision-weighted mean")
  expect_true(
    stats::sd(means[, 2]) > 1e-8,
    "a nonzero psd component is not constant"
  )
}

test_correlation_estimators <- function() {
  z <- matrix(stats::rnorm(4000 * 4), 4000, 4)
  colnames(z) <- paste0("t", 1:4)
  for (method in c("pearson", "uncentred", "spearman")) {
    correlation <- estimate_common_correlation(z, method)
    expect_true(identical(dim(correlation), c(4L, 4L)),
                paste(method, "correlation has the right dimension"))
    expect_close(diag(correlation), rep(1, 4), 1e-12,
                 paste(method, "correlation has unit diagonal"))
    expect_close(correlation, t(correlation), 1e-12,
                 paste(method, "correlation is symmetric"))
    expect_true(max(abs(correlation[upper.tri(correlation)])) < 0.08,
                paste(method, "recovers near-independence on independent z"))
  }
  expect_error(estimate_common_correlation(matrix(1, 2, 2)),
               "estimate_common_correlation rejects a tiny matrix")
}

test_lag_and_diagnostics <- function() {
  n <- 6L
  correlation <- 0.5^abs(outer(seq_len(n), seq_len(n), "-"))
  profile <- correlation_lag_profile(correlation, "ar1")
  expect_close(profile$mean_correlation, 0.5^(1:(n - 1)), 1e-12,
               "lag profile recovers the AR(1) autocorrelation")
  diagnostics <- correlation_diagnostics(correlation, "ar1")
  expect_true(diagnostics$positive_definite,
              "AR(1) correlation is positive definite")
  expect_true(diagnostics$condition_number > 1,
              "condition number exceeds one")
}

test_simulation_is_reproducible <- function() {
  means <- matrix(stats::rnorm(50 * 4), 50, 4)
  se <- matrix(stats::runif(50 * 4, 0.2, 0.5), 50, 4)
  first <- simulate_independent_replicate(means, se, 4242)
  second <- simulate_independent_replicate(means, se, 4242)
  third <- simulate_independent_replicate(means, se, 4243)
  expect_close(first, second, 0,
               "the independence replicate is reproducible for a fixed seed")
  expect_true(max(abs(first - third)) > 0,
              "different seeds give different replicates")
  standardised <- (first - means) / se
  expect_true(abs(stats::sd(standardised) - 1) < 0.1,
              "the replicate noise has the requested scale")
}

test_analytic_second_moment_matches_monte_carlo <- function() {
  # With the mixture weights held fixed the residual map is linear, so the
  # analytic second moment must match a direct simulation of independent
  # errors through the same posterior-mean machinery.
  order <- 1L
  shared_design <- build_shared_design(time_grid, num_basis, order, betaprec)
  psd_values <- c(0, 0.02, 0.1, 0.4)
  scales <- psd_to_prior_scale(psd_values, order, pred_step)$prior_scale
  weights <- c(0.5, 0.2, 0.2, 0.1)
  se <- stats::runif(n_time, 0.2, 0.7)
  truth <- sin(time_grid / 5)

  analytic <- unit_residual_second_moment(
    standard_errors = se, weights = weights,
    shared_design = shared_design, prior_scale = scales
  )

  n_draw <- 40000L
  set.seed(31337)
  z <- matrix(NA_real_, n_draw, n_time)
  for (draw in seq_len(n_draw)) {
    y <- truth + stats::rnorm(n_time, sd = se)
    component_means <- unit_component_posterior_means(
      y = y, precision = 1 / se^2, shared_design = shared_design,
      prior_scale = scales
    )
    z[draw, ] <- (y - as.numeric(component_means %*% weights)) / se
  }
  empirical <- crossprod(z) / n_draw
  # The analytic moment is centred on the residual of `truth`, which the
  # simulation also carries, so compare the centred second moments.
  offset <- as.numeric(
    (truth - as.numeric(unit_component_posterior_means(
      y = truth, precision = 1 / se^2, shared_design = shared_design,
      prior_scale = scales
    ) %*% weights)) / se
  )
  empirical <- empirical - tcrossprod(offset)
  difference <- max(abs(empirical - analytic))
  cat("  analytic vs simulated residual second moment, maximum difference:",
      format(difference, digits = 3), "\n")
  expect_true(difference < 0.03,
              "analytic residual second moment matches simulation")

  correlation_analytic <- analytic / tcrossprod(sqrt(diag(analytic)))
  expect_true(
    max(abs(correlation_analytic[upper.tri(correlation_analytic)])) > 0.02,
    "shrinkage alone induces non-negligible residual correlation"
  )
  expect_true(
    all(sqrt(diag(analytic)) < 1),
    "shrinkage deflates the residual z scale below one"
  )
}

test_transfer_operator_and_deconvolution <- function() {
  order <- 1L
  n_unit <- 40L
  shared_design <- build_shared_design(time_grid, num_basis, order, betaprec)
  psd_values <- c(0, 0.02, 0.1, 0.4)
  scales <- psd_to_prior_scale(psd_values, order, pred_step)$prior_scale
  standard_errors <- matrix(stats::runif(n_unit * n_time, 0.2, 0.8),
                            n_unit, n_time)
  weights <- matrix(stats::runif(n_unit * length(scales)), n_unit,
                    length(scales))
  weights <- weights / rowSums(weights)

  # The transfer operator must reproduce the second moment it was derived from.
  transfer <- unit_residual_transfer(standard_errors[1, ], weights[1, ],
                                     shared_design, scales)
  expect_close(
    tcrossprod(transfer),
    unit_residual_second_moment(standard_errors[1, ], weights[1, ],
                                shared_design, scales),
    1e-12,
    "the transfer operator reproduces the residual second moment"
  )

  kernel <- residual_transfer_kronecker(standard_errors, weights,
                                        shared_design, scales)

  # Forward then inverse must be the identity map on a known correlation.
  truth_correlation <- 0.6^abs(outer(seq_len(n_time), seq_len(n_time), "-"))
  predicted <- matrix(kernel %*% as.numeric(truth_correlation), n_time)
  recovered <- deconvolve_error_correlation(kernel, predicted)
  expect_close(recovered$solution, truth_correlation, 1e-4,
               "deconvolution recovers a known AR(1) error correlation")
  identity_predicted <- matrix(
    kernel %*% as.numeric(diag(1, n_time)), n_time
  )
  expect_close(
    deconvolve_error_correlation(kernel, identity_predicted)$solution,
    diag(1, n_time), 1e-4,
    "deconvolution recovers the identity"
  )
  cat("  transfer-operator condition number:",
      format(recovered$condition_number, digits = 4),
      "with effective rank", recovered$effective_rank, "of",
      n_time * (n_time + 1) / 2, "\n")
  expect_true(recovered$condition_number > 1,
              "the transfer operator is not perfectly conditioned")
}

test_symmetric_basis <- function() {
  basis <- symmetric_basis(4L)
  expect_true(ncol(basis$basis) == 10L,
              "the symmetric basis for dimension four has ten columns")
  reconstructed <- matrix(basis$basis %*% rep(1, 10), 4, 4)
  expect_close(reconstructed, matrix(1, 4, 4), 1e-12,
               "the symmetric basis sums to the all-ones matrix")
  expect_error(symmetric_basis(1L), "symmetric_basis rejects dimension one")
}

cat("Running residual-correlation helper tests.\n")
test_closed_form_matches_tmb()
test_bma_matches_predict_fash()
test_null_component_is_weighted_mean()
test_correlation_estimators()
test_lag_and_diagnostics()
test_simulation_is_reproducible()
test_analytic_second_moment_matches_monte_carlo()
test_transfer_operator_and_deconvolution()
test_symmetric_basis()
cat("Residual-correlation helper tests passed.\n")
