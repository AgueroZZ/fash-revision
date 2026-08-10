#!/usr/bin/env Rscript

parse_arguments <- function(args) {
  settings <- list(
    output_dir = file.path(
      "output",
      "revision_simulations",
      "internal",
      "omega_precision_sanity_check"
    ),
    label = "installed",
    posterior_samples = 30000L,
    simulation_units = 120L,
    simulation_seeds = c(20260807L, 20260808L, 20260809L)
  )

  index <- 1L
  while (index <= length(args)) {
    argument <- args[[index]]
    if (argument == "--output-dir") {
      index <- index + 1L
      settings$output_dir <- args[[index]]
    } else if (argument == "--label") {
      index <- index + 1L
      settings$label <- args[[index]]
    } else if (argument == "--posterior-samples") {
      index <- index + 1L
      settings$posterior_samples <- as.integer(args[[index]])
    } else if (argument == "--simulation-units") {
      index <- index + 1L
      settings$simulation_units <- as.integer(args[[index]])
    } else if (argument == "--simulation-seeds") {
      index <- index + 1L
      settings$simulation_seeds <- as.integer(strsplit(args[[index]], ",", fixed = TRUE)[[1]])
    } else {
      stop("Unknown argument: ", argument)
    }
    index <- index + 1L
  }

  if (!nzchar(settings$label)) {
    stop("--label must be non-empty.")
  }
  if (is.na(settings$posterior_samples) || settings$posterior_samples < 5000L) {
    stop("--posterior-samples must be at least 5000.")
  }
  if (is.na(settings$simulation_units) || settings$simulation_units < 30L) {
    stop("--simulation-units must be at least 30.")
  }
  if (anyNA(settings$simulation_seeds)) {
    stop("--simulation-seeds must be comma-separated integers.")
  }

  settings
}

log_mvn_zero <- function(residual, covariance) {
  chol_covariance <- chol(covariance)
  whitened <- forwardsolve(t(chol_covariance), residual)
  -0.5 * (
    length(residual) * log(2 * pi) +
      2 * sum(log(diag(chol_covariance))) +
      sum(whitened^2)
  )
}

sigma_iwp_from_psd <- function(psd, order, pred_step) {
  psd / sqrt(
    pred_step^(2 * order - 1) /
      ((2 * order - 1) * factorial(order - 1)^2)
  )
}

model_matrices <- function(data_i, precision, psd, num_basis,
                           betaprec, order, pred_step) {
  tmb_data <- fashr::fash_set_tmbdat(
    data_i = data_i,
    Si = NULL,
    Omegai = precision,
    num_basis = num_basis,
    betaprec = betaprec,
    order = order
  )

  if (psd == 0) {
    design <- as.matrix(tmb_data$X)
    prior_precision <- diag(betaprec, ncol(design))
  } else {
    sigma_iwp <- sigma_iwp_from_psd(psd, order, pred_step)
    design <- cbind(as.matrix(tmb_data$B), as.matrix(tmb_data$X))
    prior_precision <- as.matrix(Matrix::bdiag(
      as.matrix(tmb_data$P) / sigma_iwp^2,
      diag(betaprec, ncol(tmb_data$X))
    ))
  }

  list(
    tmb_data = tmb_data,
    design = design,
    prior_precision = prior_precision
  )
}

analytic_log_marginal <- function(data_i, precision, psd, num_basis,
                                  betaprec, order, pred_step) {
  matrices <- model_matrices(
    data_i = data_i,
    precision = precision,
    psd = psd,
    num_basis = num_basis,
    betaprec = betaprec,
    order = order,
    pred_step = pred_step
  )
  observation_covariance <- solve(precision)
  prior_covariance <- solve(matrices$prior_precision)
  marginal_covariance <- observation_covariance +
    matrices$design %*% prior_covariance %*% t(matrices$design)
  log_mvn_zero(data_i$y - data_i$offset, marginal_covariance)
}

analytic_posterior_paths <- function(data_i, precision, psd, num_basis,
                                     betaprec, order, pred_step) {
  matrices <- model_matrices(
    data_i = data_i,
    precision = precision,
    psd = psd,
    num_basis = num_basis,
    betaprec = betaprec,
    order = order,
    pred_step = pred_step
  )
  posterior_precision <- matrices$prior_precision +
    crossprod(matrices$design, precision %*% matrices$design)
  posterior_covariance <- solve(posterior_precision)
  posterior_mean <- posterior_covariance %*%
    crossprod(matrices$design, precision %*% (data_i$y - data_i$offset))

  list(
    mean = as.vector(matrices$design %*% posterior_mean),
    covariance = matrices$design %*% posterior_covariance %*% t(matrices$design),
    coefficient_mean = as.vector(posterior_mean),
    coefficient_covariance = posterior_covariance
  )
}

ar1_correlation <- function(size, rho) {
  outer(seq_len(size), seq_len(size), function(i, j) rho^abs(i - j))
}

expand_prior_weights <- function(fit, grid) {
  result <- setNames(numeric(length(grid)), as.character(grid))
  indices <- match(as.character(fit$prior_weights$psd), names(result))
  result[indices] <- fit$prior_weights$prior_weight
  result
}

expand_posterior_weights <- function(fit, grid) {
  result <- matrix(
    0,
    nrow = nrow(fit$posterior_weights),
    ncol = length(grid),
    dimnames = list(NULL, as.character(grid))
  )
  result[, colnames(fit$posterior_weights)] <- fit$posterior_weights
  result
}

simulate_model_aligned_data <- function(seed, n_units, grid, true_weights,
                                        x, observation_covariance,
                                        num_basis, betaprec, order,
                                        pred_step) {
  set.seed(seed)
  component <- sample(
    seq_along(grid),
    size = n_units,
    replace = TRUE,
    prob = true_weights
  )
  precision <- solve(observation_covariance)
  template_data <- data.frame(y = numeric(length(x)), x = x, offset = 0)
  error_factor <- chol(observation_covariance)

  response <- matrix(NA_real_, nrow = n_units, ncol = length(x))
  latent_mean <- matrix(NA_real_, nrow = n_units, ncol = length(x))

  for (unit in seq_len(n_units)) {
    matrices <- model_matrices(
      data_i = template_data,
      precision = precision,
      psd = grid[[component[[unit]]]],
      num_basis = num_basis,
      betaprec = betaprec,
      order = order,
      pred_step = pred_step
    )
    coefficient <- as.vector(
      rnorm(ncol(matrices$design)) %*% chol(solve(matrices$prior_precision))
    )
    latent_mean[unit, ] <- as.vector(matrices$design %*% coefficient)
    response[unit, ] <- latent_mean[unit, ] +
      as.vector(rnorm(length(x)) %*% error_factor)
  }

  list(
    response = response,
    latent_mean = latent_mean,
    component = component
  )
}

summarize_simulation_fit <- function(fit, method, seed, grid,
                                     true_weights, true_component) {
  prior_weights <- expand_prior_weights(fit, grid)
  posterior_weights <- expand_posterior_weights(fit, grid)
  fitted_component <- max.col(posterior_weights, ties.method = "first")
  true_probability <- posterior_weights[
    cbind(seq_along(true_component), true_component)
  ]

  data.frame(
    seed = seed,
    method = method,
    n_units = length(true_component),
    prior_l1_error = sum(abs(prior_weights - true_weights)),
    prior_total_variation = 0.5 * sum(abs(prior_weights - true_weights)),
    classification_accuracy = mean(fitted_component == true_component),
    mean_true_component_probability = mean(true_probability),
    mean_true_component_log_score = mean(log(pmax(true_probability, 1e-300))),
    true_weight_0 = true_weights[[1]],
    estimated_weight_0 = prior_weights[[1]],
    true_weight_middle = true_weights[[2]],
    estimated_weight_middle = prior_weights[[2]],
    true_weight_high = true_weights[[3]],
    estimated_weight_high = prior_weights[[3]],
    stringsAsFactors = FALSE
  )
}

settings <- parse_arguments(commandArgs(trailingOnly = TRUE))
dir.create(settings$output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(fashr)
  library(Matrix)
})

package_description <- utils::packageDescription("fashr")
package_info <- data.frame(
  label = settings$label,
  version = as.character(utils::packageVersion("fashr")),
  package_path = find.package("fashr"),
  built = if (is.null(package_description$Built)) NA_character_ else package_description$Built,
  stringsAsFactors = FALSE
)

num_basis <- 6L
betaprec <- 0.8
order <- 1L
pred_step <- 1
x <- seq(0, 7, length.out = 8)
offset <- 0.15 * cos(x / 2)
y <- c(-0.45, 0.10, 0.62, 0.91, 0.25, -0.18, 0.36, 0.82)
data_i <- data.frame(y = y, x = x, offset = offset)
standard_errors <- seq(0.55, 0.95, length.out = length(x))
diagonal_covariance <- diag(standard_errors^2)
correlated_covariance <- diag(standard_errors) %*%
  ar1_correlation(length(x), rho = 0.65) %*%
  diag(standard_errors)
precision_list <- list(
  diagonal = solve(diagonal_covariance),
  correlated = solve(correlated_covariance)
)
grid <- c(0, 0.35, 0.9)

exact_rows <- list()
row_index <- 1L
for (precision_name in names(precision_list)) {
  precision <- precision_list[[precision_name]]
  for (psd in grid) {
    package_value <- fashr:::compute_L_gaussian_helper(
      data_i = data_i,
      Si = NULL,
      Omegai = precision,
      psd_iwp = psd,
      num_basis = num_basis,
      betaprec = betaprec,
      order = order,
      pred_step = pred_step
    )
    analytic_value <- analytic_log_marginal(
      data_i = data_i,
      precision = precision,
      psd = psd,
      num_basis = num_basis,
      betaprec = betaprec,
      order = order,
      pred_step = pred_step
    )
    absolute_error <- abs(package_value - analytic_value)
    exact_rows[[row_index]] <- data.frame(
      check = "component_log_marginal",
      scenario = precision_name,
      psd = psd,
      package_value = package_value,
      reference_value = analytic_value,
      absolute_error = absolute_error,
      tolerance = 1e-6,
      passed = absolute_error < 1e-6,
      stringsAsFactors = FALSE
    )
    row_index <- row_index + 1L
  }
}

for (psd in grid) {
  omega_value <- fashr:::compute_L_gaussian_helper(
    data_i = data_i,
    Si = NULL,
    Omegai = precision_list$diagonal,
    psd_iwp = psd,
    num_basis = num_basis,
    betaprec = betaprec,
    order = order,
    pred_step = pred_step
  )
  standard_error_value <- fashr:::compute_L_gaussian_helper(
    data_i = data_i,
    Si = standard_errors,
    Omegai = NULL,
    psd_iwp = psd,
    num_basis = num_basis,
    betaprec = betaprec,
    order = order,
    pred_step = pred_step
  )
  absolute_error <- abs(omega_value - standard_error_value)
  exact_rows[[row_index]] <- data.frame(
    check = "diagonal_omega_equals_standard_errors",
    scenario = "diagonal",
    psd = psd,
    package_value = omega_value,
    reference_value = standard_error_value,
    absolute_error = absolute_error,
    tolerance = 1e-6,
    passed = absolute_error < 1e-6,
    stringsAsFactors = FALSE
  )
  row_index <- row_index + 1L
}

public_response <- rbind(
  y,
  y + c(0.12, -0.04, 0.08, -0.05, 0.03, 0.06, -0.09, 0.04),
  y + c(-0.10, 0.07, -0.06, 0.11, -0.05, 0.02, 0.08, -0.03),
  y + c(0.04, 0.09, -0.11, 0.06, 0.10, -0.08, 0.05, -0.07)
)
public_fit <- fashr::fash(
  Y = public_response,
  smooth_var = x,
  offset = offset,
  S = NULL,
  Omega = precision_list$correlated,
  grid = grid,
  likelihood = "gaussian",
  num_basis = num_basis,
  betaprec = betaprec,
  order = order,
  pred_step = pred_step,
  penalty = 1,
  num_cores = 1,
  verbose = FALSE
)
public_reference <- vapply(
  seq_len(nrow(public_response)),
  function(unit) {
    vapply(
      grid,
      function(psd) {
        analytic_log_marginal(
          data_i = public_fit$fash_data$data_list[[unit]],
          precision = precision_list$correlated,
          psd = psd,
          num_basis = num_basis,
          betaprec = betaprec,
          order = order,
          pred_step = pred_step
        )
      },
      numeric(1)
    )
  },
  numeric(length(grid))
)
public_reference <- t(public_reference)
public_absolute_error <- max(abs(public_fit$L_matrix - public_reference))
exact_rows[[row_index]] <- data.frame(
  check = "public_fash_l_matrix",
  scenario = "correlated",
  psd = NA_real_,
  package_value = max(public_fit$L_matrix),
  reference_value = max(public_reference),
  absolute_error = public_absolute_error,
  tolerance = 1e-6,
  passed = public_absolute_error < 1e-6,
  stringsAsFactors = FALSE
)
row_index <- row_index + 1L

posterior_psd <- 0.9
posterior_reference <- analytic_posterior_paths(
  data_i = data_i,
  precision = precision_list$correlated,
  psd = posterior_psd,
  num_basis = num_basis,
  betaprec = betaprec,
  order = order,
  pred_step = pred_step
)
set.seed(20260807L)
posterior_samples <- fashr:::fash_fit_once(
  data_i = data_i,
  refined_x = x,
  M = settings$posterior_samples,
  psd_iwp = posterior_psd,
  Si = NULL,
  Omegai = precision_list$correlated,
  num_basis = num_basis,
  betaprec = betaprec,
  order = order,
  pred_step = pred_step,
  likelihood = "gaussian",
  deriv = 0
)
sample_mean <- rowMeans(posterior_samples)
sample_variance <- apply(posterior_samples, 1, stats::var)
reference_variance <- diag(posterior_reference$covariance)
mean_monte_carlo_z <- max(
  abs(sample_mean - posterior_reference$mean) /
    sqrt(reference_variance / settings$posterior_samples)
)
variance_relative_error <- max(
  abs(sample_variance - reference_variance) / reference_variance
)
exact_rows[[row_index]] <- data.frame(
  check = "posterior_path_mean_monte_carlo_z",
  scenario = "correlated",
  psd = posterior_psd,
  package_value = mean_monte_carlo_z,
  reference_value = 0,
  absolute_error = mean_monte_carlo_z,
  tolerance = 6,
  passed = mean_monte_carlo_z < 6,
  stringsAsFactors = FALSE
)
row_index <- row_index + 1L
exact_rows[[row_index]] <- data.frame(
  check = "posterior_path_variance_relative_error",
  scenario = "correlated",
  psd = posterior_psd,
  package_value = variance_relative_error,
  reference_value = 0,
  absolute_error = variance_relative_error,
  tolerance = 0.08,
  passed = variance_relative_error < 0.08,
  stringsAsFactors = FALSE
)

exact_checks <- do.call(rbind, exact_rows)

simulation_grid <- c(0, 0.45, 0.9)
true_weights <- c(0.50, 0.30, 0.20)
simulation_x <- 0:9
simulation_standard_error <- 0.65
simulation_covariance <- simulation_standard_error^2 *
  ar1_correlation(length(simulation_x), rho = 0.75)
simulation_precision <- solve(simulation_covariance)
simulation_rows <- list()
simulation_index <- 1L

for (simulation_seed in settings$simulation_seeds) {
  simulated <- simulate_model_aligned_data(
    seed = simulation_seed,
    n_units = settings$simulation_units,
    grid = simulation_grid,
    true_weights = true_weights,
    x = simulation_x,
    observation_covariance = simulation_covariance,
    num_basis = num_basis,
    betaprec = betaprec,
    order = order,
    pred_step = pred_step
  )

  correct_fit <- fashr::fash(
    Y = simulated$response,
    smooth_var = simulation_x,
    offset = 0,
    S = NULL,
    Omega = simulation_precision,
    grid = simulation_grid,
    likelihood = "gaussian",
    num_basis = num_basis,
    betaprec = betaprec,
    order = order,
    pred_step = pred_step,
    penalty = 1,
    num_cores = 1,
    verbose = FALSE
  )
  diagonal_fit <- fashr::fash(
    Y = simulated$response,
    smooth_var = simulation_x,
    offset = 0,
    S = rep(simulation_standard_error, length(simulation_x)),
    Omega = NULL,
    grid = simulation_grid,
    likelihood = "gaussian",
    num_basis = num_basis,
    betaprec = betaprec,
    order = order,
    pred_step = pred_step,
    penalty = 1,
    num_cores = 1,
    verbose = FALSE
  )

  simulation_rows[[simulation_index]] <- summarize_simulation_fit(
    fit = correct_fit,
    method = "correct_full_precision",
    seed = simulation_seed,
    grid = simulation_grid,
    true_weights = true_weights,
    true_component = simulated$component
  )
  simulation_index <- simulation_index + 1L
  simulation_rows[[simulation_index]] <- summarize_simulation_fit(
    fit = diagonal_fit,
    method = "diagonal_independence",
    seed = simulation_seed,
    grid = simulation_grid,
    true_weights = true_weights,
    true_component = simulated$component
  )
  simulation_index <- simulation_index + 1L
}

simulation_summary <- do.call(rbind, simulation_rows)

artifact_prefix <- file.path(settings$output_dir, settings$label)
utils::write.csv(
  package_info,
  paste0(artifact_prefix, "_package_info.csv"),
  row.names = FALSE
)
utils::write.csv(
  exact_checks,
  paste0(artifact_prefix, "_exact_checks.csv"),
  row.names = FALSE
)
utils::write.csv(
  simulation_summary,
  paste0(artifact_prefix, "_simulation_summary.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    package_info = package_info,
    settings = settings,
    exact_checks = exact_checks,
    simulation_summary = simulation_summary,
    exact_inputs = list(
      data = data_i,
      standard_errors = standard_errors,
      diagonal_covariance = diagonal_covariance,
      correlated_covariance = correlated_covariance,
      grid = grid,
      num_basis = num_basis,
      betaprec = betaprec,
      order = order,
      pred_step = pred_step
    ),
    simulation_inputs = list(
      grid = simulation_grid,
      true_weights = true_weights,
      x = simulation_x,
      standard_error = simulation_standard_error,
      rho = 0.75
    )
  ),
  paste0(artifact_prefix, "_audit_bundle.rds")
)
writeLines(
  capture.output(utils::sessionInfo()),
  paste0(artifact_prefix, "_session_info.txt")
)

cat("Package label:", settings$label, "\n")
cat("Package version:", package_info$version, "\n")
cat("Package path:", package_info$package_path, "\n")
cat("Exact checks passed:", sum(exact_checks$passed), "/", nrow(exact_checks), "\n")
cat("Maximum analytic log-marginal error:",
    max(exact_checks$absolute_error[exact_checks$check == "component_log_marginal"]),
    "\n")
cat("Maximum public L-matrix error:", public_absolute_error, "\n")
cat("Posterior mean maximum Monte Carlo z-score:", mean_monte_carlo_z, "\n")
cat("Posterior variance maximum relative error:", variance_relative_error, "\n")
print(simulation_summary)

if (!all(exact_checks$passed)) {
  failed <- exact_checks[!exact_checks$passed, , drop = FALSE]
  print(failed)
  stop("One or more exact Omega checks failed.")
}
