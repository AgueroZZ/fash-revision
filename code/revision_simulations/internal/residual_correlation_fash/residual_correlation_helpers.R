# Helpers for the residual-based common-correlation FASH experiment.
#
# The pipeline is:
#   1. ordinary conditional-independence FASH on a fixed subset;
#   2. closed-form Bayesian-model-averaged posterior means;
#   3. residual z-scores r / se;
#   4. one common cross-time correlation matrix;
#   5. a dependent-likelihood refit built from that matrix.
#
# The posterior mean is computed in closed form rather than by
# `fashr:::predict.fash` Monte Carlo. For a fixed PSD value the FASH model is
# linear-Gaussian, so the mode equals the mean and
#
#   w(c) = (D' Omega D + c * Pfull)^{-1} D' Omega (y - offset),
#   fitted(c) = D w(c),
#
# with D = [B, X] the O-spline basis stacked with the global polynomial,
# Pfull = bdiag(P, betaprec * I) and c = 1 / sigmaIWP^2. This reproduces the
# package's own TMB optimum to machine precision; see
# `test_residual_correlation_helpers.R`.
#
# `fashr:::compute_marginal_mean_var` is deliberately NOT used: it returns only
# the O-spline block, drops the global polynomial, and returns exactly zero at
# psd = 0.

build_shared_design <- function(time_grid,
                                num_basis,
                                order,
                                betaprec) {
  if (!requireNamespace("fashr", quietly = TRUE)) {
    stop("The fashr package is required to build the shared design.")
  }
  time_grid <- as.numeric(time_grid)
  num_basis <- as.integer(num_basis)
  order <- as.integer(order)
  betaprec <- as.numeric(betaprec)
  if (length(time_grid) < 3L || any(!is.finite(time_grid)) ||
      anyDuplicated(time_grid) || is.unsorted(time_grid) ||
      length(num_basis) != 1L || is.na(num_basis) || num_basis < 2L ||
      length(order) != 1L || is.na(order) || order < 1L ||
      length(betaprec) != 1L || !is.finite(betaprec) || betaprec < 0) {
    stop("The shared-design inputs are invalid.")
  }
  namespace <- asNamespace("fashr")
  knots <- seq(min(time_grid), max(time_grid), length.out = num_basis)
  global_design <- as.matrix(
    namespace$global_poly_helper(time_grid, p = order)
  )
  spline_design <- as.matrix(
    namespace$local_poly_helper(
      knots = knots,
      refined_x = time_grid,
      p = order
    )
  )
  penalty_matrix <- as.matrix(
    namespace$compute_weights_precision_helper(knots)
  )
  if (nrow(global_design) != length(time_grid) ||
      nrow(spline_design) != length(time_grid) ||
      !identical(dim(penalty_matrix), rep(ncol(spline_design), 2L))) {
    stop("The fashr design helpers returned inconsistent dimensions.")
  }
  design <- cbind(spline_design, global_design)
  spline_index <- seq_len(ncol(spline_design))
  global_index <- ncol(spline_design) + seq_len(ncol(global_design))

  # The PSD-dependent scale multiplies the O-spline penalty only. `betaprec`
  # is a fixed ridge on the global-polynomial block and must not be rescaled.
  scaled_prior <- matrix(0, nrow = ncol(design), ncol = ncol(design))
  scaled_prior[spline_index, spline_index] <- penalty_matrix
  fixed_prior <- matrix(0, nrow = ncol(design), ncol = ncol(design))
  diag(fixed_prior)[global_index] <- betaprec

  list(
    time_grid = time_grid,
    knots = knots,
    spline_design = spline_design,
    global_design = global_design,
    design = design,
    penalty_matrix = penalty_matrix,
    scaled_prior_precision = scaled_prior,
    fixed_prior_precision = fixed_prior,
    spline_index = spline_index,
    global_index = global_index,
    num_basis = num_basis,
    order = order,
    betaprec = betaprec
  )
}

# `fash_eb_est` prunes grid components whose estimated prior weight is zero, so
# the columns of `posterior_weights` index `prior_weights$psd`, NOT `psd_grid`.
# Everything downstream must use these retained values.
retained_psd_values <- function(fit) {
  if (!is.list(fit) || is.null(fit$prior_weights) ||
      is.null(fit$posterior_weights)) {
    stop("retained_psd_values requires a fitted FASH object.")
  }
  psd <- as.numeric(fit$prior_weights$psd)
  if (length(psd) < 1L || any(!is.finite(psd)) || any(psd < 0) ||
      anyDuplicated(psd) || is.unsorted(psd) ||
      ncol(as.matrix(fit$posterior_weights)) != length(psd)) {
    stop("The retained PSD values are not aligned with the posterior weights.")
  }
  psd
}

psd_to_prior_scale <- function(psd_values, order, pred_step) {
  psd_values <- as.numeric(psd_values)
  order <- as.integer(order)
  pred_step <- as.numeric(pred_step)
  if (length(psd_values) < 2L || any(!is.finite(psd_values)) ||
      any(psd_values < 0) ||
      length(order) != 1L || is.na(order) || order < 1L ||
      length(pred_step) != 1L || !is.finite(pred_step) || pred_step <= 0) {
    stop("psd_to_prior_scale received invalid inputs.")
  }
  conversion <- sqrt(
    (pred_step^((2 * order) - 1)) /
      (((2 * order) - 1) * (factorial(order - 1)^2))
  )
  sigma_iwp <- psd_values / conversion
  scale <- rep(NA_real_, length(psd_values))
  nonzero <- psd_values != 0
  scale[nonzero] <- 1 / (sigma_iwp[nonzero]^2)
  if (!any(nonzero) || any(is.finite(scale) & scale <= 0)) {
    stop("The PSD grid did not yield a usable set of positive prior scales.")
  }
  list(sigma_iwp = sigma_iwp, prior_scale = scale, conversion = conversion)
}

# Component-wise posterior means for one unit, evaluated at the observed
# times. Returns a length(time) x length(psd_values) matrix.
#
# All nonzero-PSD components share one Cholesky factor and one symmetric
# eigen-decomposition. Anchoring at the smallest prior scale `a` present in the
# grid,
#
#   H(c) = H(a) + (c - a) Pscaled = R' U (I + (c - a) Lambda) U' R,
#
# so every component is a diagonal rescaling of the same n_time x p
# projection. Lambda is nonnegative because Pscaled is positive semidefinite,
# and c >= a by construction, so H(c) stays positive definite throughout.
unit_component_posterior_means <- function(y,
                                           precision,
                                           shared_design,
                                           prior_scale,
                                           offset = 0) {
  y <- as.numeric(y)
  n_time <- length(y)
  offset <- if (length(offset) == 1L) rep(offset, n_time) else as.numeric(offset)
  prior_scale <- as.numeric(prior_scale)
  design <- shared_design$design
  global_design <- shared_design$global_design
  scaled_prior <- shared_design$scaled_prior_precision
  fixed_prior <- shared_design$fixed_prior_precision
  if (n_time != nrow(design) || length(offset) != n_time ||
      any(!is.finite(y)) || any(!is.finite(offset)) ||
      length(prior_scale) < 2L || !any(is.finite(prior_scale))) {
    stop("unit_component_posterior_means received invalid unit inputs.")
  }
  anchor_scale <- min(prior_scale[is.finite(prior_scale)])
  if (is.matrix(precision)) {
    if (!identical(dim(precision), c(n_time, n_time)) ||
        any(!is.finite(precision))) {
      stop("The unit precision matrix is invalid.")
    }
    omega <- precision
  } else {
    precision <- as.numeric(precision)
    if (length(precision) != n_time || any(!is.finite(precision)) ||
        any(precision <= 0)) {
      stop("The unit precision diagonal is invalid.")
    }
    omega <- diag(precision, n_time)
  }
  centred <- y - offset
  cross_product <- crossprod(design, omega) %*% design
  cross_product <- (cross_product + t(cross_product)) / 2
  moment <- as.numeric(crossprod(design, omega) %*% centred)

  anchor <- cross_product + fixed_prior + anchor_scale * scaled_prior
  anchor <- (anchor + t(anchor)) / 2
  root <- chol(anchor)
  whitened_prior <- backsolve(root, scaled_prior, transpose = TRUE)
  whitened_prior <- t(backsolve(root, t(whitened_prior), transpose = TRUE))
  whitened_prior <- (whitened_prior + t(whitened_prior)) / 2
  spectrum <- eigen(whitened_prior, symmetric = TRUE)
  projection <- design %*% backsolve(root, spectrum$vectors)
  rotated_moment <- as.numeric(
    crossprod(spectrum$vectors, backsolve(root, moment, transpose = TRUE))
  )

  null_precision <- crossprod(global_design, omega) %*% global_design
  null_precision <- (null_precision + t(null_precision)) / 2
  diag(null_precision) <- diag(null_precision) + shared_design$betaprec
  null_moment <- as.numeric(crossprod(global_design, omega) %*% centred)
  null_fitted <- as.numeric(
    global_design %*% solve(null_precision, null_moment)
  )

  means <- matrix(0, nrow = n_time, ncol = length(prior_scale))
  for (index in seq_along(prior_scale)) {
    scale <- prior_scale[index]
    if (is.na(scale)) {
      means[, index] <- null_fitted + offset
    } else {
      shrunk <- rotated_moment /
        (1 + (scale - anchor_scale) * spectrum$values)
      means[, index] <- as.numeric(projection %*% shrunk) + offset
    }
  }
  means
}

# Bayesian-model-averaged posterior means for every unit.
#
# `posterior_weights` is units x components and must be row-normalised.
bma_posterior_means <- function(beta_hat,
                                standard_errors,
                                posterior_weights,
                                shared_design,
                                prior_scale,
                                offset = 0,
                                precision_matrices = NULL,
                                num_cores = 1L) {
  beta_hat <- as.matrix(beta_hat)
  standard_errors <- as.matrix(standard_errors)
  posterior_weights <- as.matrix(posterior_weights)
  num_cores <- as.integer(num_cores)
  n_unit <- nrow(beta_hat)
  n_time <- ncol(beta_hat)
  if (!identical(dim(beta_hat), dim(standard_errors)) ||
      nrow(posterior_weights) != n_unit ||
      ncol(posterior_weights) != length(prior_scale) ||
      n_time != nrow(shared_design$design) ||
      any(!is.finite(beta_hat)) || any(!is.finite(standard_errors)) ||
      any(standard_errors <= 0) || any(!is.finite(posterior_weights)) ||
      any(posterior_weights < 0) ||
      max(abs(rowSums(posterior_weights) - 1)) > 1e-8 ||
      length(num_cores) != 1L || is.na(num_cores) || num_cores < 1L) {
    stop("bma_posterior_means received invalid or misaligned inputs.")
  }
  if (!is.null(precision_matrices) &&
      length(precision_matrices) != n_unit) {
    stop("precision_matrices must have one entry per unit.")
  }
  one_unit <- function(index) {
    precision <- if (is.null(precision_matrices)) {
      1 / (standard_errors[index, ]^2)
    } else {
      precision_matrices[[index]]
    }
    component_means <- unit_component_posterior_means(
      y = beta_hat[index, ],
      precision = precision,
      shared_design = shared_design,
      prior_scale = prior_scale,
      offset = offset
    )
    as.numeric(component_means %*% posterior_weights[index, ])
  }
  results <- if (num_cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(n_unit), one_unit, mc.cores = num_cores)
  } else {
    lapply(seq_len(n_unit), one_unit)
  }
  failed <- !vapply(results, function(x) {
    is.numeric(x) && length(x) == n_time && all(is.finite(x))
  }, logical(1))
  if (any(failed)) {
    stop(sum(failed), " units returned an invalid posterior mean.")
  }
  means <- do.call(rbind, results)
  dimnames(means) <- dimnames(beta_hat)
  means
}

residual_z_scores <- function(beta_hat, posterior_means, standard_errors) {
  beta_hat <- as.matrix(beta_hat)
  posterior_means <- as.matrix(posterior_means)
  standard_errors <- as.matrix(standard_errors)
  if (!identical(dim(beta_hat), dim(posterior_means)) ||
      !identical(dim(beta_hat), dim(standard_errors)) ||
      any(!is.finite(beta_hat)) || any(!is.finite(posterior_means)) ||
      any(!is.finite(standard_errors)) || any(standard_errors <= 0)) {
    stop("residual_z_scores received invalid or misaligned inputs.")
  }
  residual <- beta_hat - posterior_means
  z <- residual / standard_errors
  list(residual = residual, z = z)
}

# Primary estimator is stats::cor. The uncentred second-moment version is
# reported alongside it because the residual z-scores are not exactly mean
# zero, and a Spearman variant is reported because they are heavy tailed.
estimate_common_correlation <- function(z, method = c("pearson",
                                                      "uncentred",
                                                      "spearman")) {
  method <- match.arg(method)
  z <- as.matrix(z)
  if (nrow(z) < 10L || ncol(z) < 2L || any(!is.finite(z))) {
    stop("estimate_common_correlation received an invalid z matrix.")
  }
  correlation <- if (method == "uncentred") {
    second_moment <- crossprod(z) / nrow(z)
    scaling <- sqrt(diag(second_moment))
    second_moment / tcrossprod(scaling)
  } else {
    stats::cor(z, method = method)
  }
  correlation <- (correlation + t(correlation)) / 2
  diag(correlation) <- 1
  dimnames(correlation) <- list(colnames(z), colnames(z))
  correlation
}

correlation_diagnostics <- function(correlation, name = "correlation") {
  correlation <- as.matrix(correlation)
  if (nrow(correlation) != ncol(correlation) ||
      any(!is.finite(correlation))) {
    stop(name, " must be a finite square matrix.")
  }
  eigenvalues <- eigen(correlation, symmetric = TRUE, only.values = TRUE)$values
  offdiag <- correlation[upper.tri(correlation)]
  data.frame(
    name = name,
    dimension = nrow(correlation),
    minimum_eigenvalue = min(eigenvalues),
    maximum_eigenvalue = max(eigenvalues),
    condition_number = max(eigenvalues) / min(eigenvalues),
    positive_definite = min(eigenvalues) > 1e-10,
    mean_offdiagonal = mean(offdiag),
    minimum_offdiagonal = min(offdiag),
    maximum_offdiagonal = max(offdiag),
    stringsAsFactors = FALSE
  )
}

correlation_lag_profile <- function(correlation, name = "correlation") {
  correlation <- as.matrix(correlation)
  n_time <- nrow(correlation)
  if (n_time != ncol(correlation) || n_time < 2L ||
      any(!is.finite(correlation))) {
    stop(name, " must be a finite square matrix.")
  }
  lags <- seq_len(n_time - 1L)
  profile <- vapply(lags, function(lag) {
    index <- seq_len(n_time - lag)
    mean(correlation[cbind(index, index + lag)])
  }, numeric(1))
  data.frame(
    name = name,
    lag = lags,
    mean_correlation = profile,
    stringsAsFactors = FALSE
  )
}

# Exact analytic companion to the bootstrap calibration.
#
# Holding the mixture weights fixed, the residual map of one unit is linear:
#
#   r_j = (I - M_j) beta_hat_j,
#   M_j = w_j0 S_j0 + Proj_j Dbar_j Proj_j' Omega_j,
#
# where `Proj_j` and the eigenvalues come from the same factorisation used for
# the posterior mean, so `M_j` costs nothing extra. Under independent errors
# with Sigma_j = diag(se_j^2) the residual z second moment is
#
#   D_j^{-1} (I - M_j) Sigma_j (I - M_j)' D_j^{-1}.
#
# Averaging that over units and converting to a correlation gives what the
# estimator must be expected to return when there is no error correlation at
# all. It is the honest null for `estimate_common_correlation`, and it is not
# the identity.
unit_residual_transfer <- function(standard_errors,
                                   weights,
                                   shared_design,
                                   prior_scale) {
  standard_errors <- as.numeric(standard_errors)
  weights <- as.numeric(weights)
  prior_scale <- as.numeric(prior_scale)
  n_time <- length(standard_errors)
  design <- shared_design$design
  global_design <- shared_design$global_design
  scaled_prior <- shared_design$scaled_prior_precision
  fixed_prior <- shared_design$fixed_prior_precision
  if (n_time != nrow(design) || any(!is.finite(standard_errors)) ||
      any(standard_errors <= 0) || length(weights) != length(prior_scale) ||
      any(!is.finite(weights)) || abs(sum(weights) - 1) > 1e-8) {
    stop("unit_residual_second_moment received invalid unit inputs.")
  }
  omega <- diag(1 / standard_errors^2, n_time)
  anchor_scale <- min(prior_scale[is.finite(prior_scale)])
  anchor <- crossprod(design, omega) %*% design + fixed_prior +
    anchor_scale * scaled_prior
  anchor <- (anchor + t(anchor)) / 2
  root <- chol(anchor)
  whitened_prior <- backsolve(root, scaled_prior, transpose = TRUE)
  whitened_prior <- t(backsolve(root, t(whitened_prior), transpose = TRUE))
  spectrum <- eigen((whitened_prior + t(whitened_prior)) / 2, symmetric = TRUE)
  projection <- design %*% backsolve(root, spectrum$vectors)

  finite <- is.finite(prior_scale)
  averaged <- if (any(finite)) {
    colSums(
      outer(weights[finite], rep(1, length(spectrum$values))) /
        (1 + outer(prior_scale[finite] - anchor_scale, spectrum$values))
    )
  } else {
    rep(0, length(spectrum$values))
  }
  smoother <- projection %*% (averaged * t(projection)) %*% omega

  null_weight <- sum(weights[!finite])
  if (null_weight > 0) {
    null_precision <- crossprod(global_design, omega) %*% global_design
    diag(null_precision) <- diag(null_precision) + shared_design$betaprec
    null_smoother <- global_design %*%
      solve(null_precision, crossprod(global_design, omega))
    smoother <- smoother + null_weight * null_smoother
  }
  residual_operator <- diag(1, n_time) - smoother
  # T_j = D_j^{-1} (I - M_j) D_j maps an error correlation on the z scale to
  # the residual-z scale, so this unit's residual-z second moment is
  # T_j P T_j' for any true error correlation P.
  residual_operator * outer(1 / standard_errors, standard_errors)
}

unit_residual_second_moment <- function(standard_errors,
                                        weights,
                                        shared_design,
                                        prior_scale) {
  tcrossprod(unit_residual_transfer(
    standard_errors = standard_errors,
    weights = weights,
    shared_design = shared_design,
    prior_scale = prior_scale
  ))
}

expected_residual_z_moment <- function(beta_hat,
                                       standard_errors,
                                       posterior_weights,
                                       shared_design,
                                       prior_scale,
                                       num_cores = 1L) {
  beta_hat <- as.matrix(beta_hat)
  standard_errors <- as.matrix(standard_errors)
  posterior_weights <- as.matrix(posterior_weights)
  n_unit <- nrow(beta_hat)
  if (!identical(dim(beta_hat), dim(standard_errors)) ||
      nrow(posterior_weights) != n_unit ||
      ncol(posterior_weights) != length(prior_scale)) {
    stop("expected_residual_z_moment received misaligned inputs.")
  }
  one_unit <- function(index) {
    unit_residual_second_moment(
      standard_errors = standard_errors[index, ],
      weights = posterior_weights[index, ],
      shared_design = shared_design,
      prior_scale = prior_scale
    )
  }
  moments <- if (num_cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(n_unit), one_unit, mc.cores = num_cores)
  } else {
    lapply(seq_len(n_unit), one_unit)
  }
  valid <- vapply(moments, function(x) {
    is.matrix(x) && all(is.finite(x))
  }, logical(1))
  if (!all(valid)) {
    stop(sum(!valid), " units returned an invalid residual second moment.")
  }
  average <- Reduce(`+`, moments) / n_unit
  average <- (average + t(average)) / 2
  scaling <- sqrt(diag(average))
  correlation <- average / tcrossprod(scaling)
  diag(correlation) <- 1
  dimnames(average) <- dimnames(correlation) <-
    list(colnames(beta_hat), colnames(beta_hat))
  list(
    second_moment = average,
    correlation = correlation,
    expected_z_sd = scaling
  )
}

# The reference distribution for the estimator: what does this exact pipeline
# return when the errors are independent by construction?
simulate_independent_replicate <- function(posterior_means,
                                           standard_errors,
                                           seed) {
  posterior_means <- as.matrix(posterior_means)
  standard_errors <- as.matrix(standard_errors)
  seed <- as.integer(seed)
  if (!identical(dim(posterior_means), dim(standard_errors)) ||
      any(!is.finite(posterior_means)) || any(!is.finite(standard_errors)) ||
      any(standard_errors <= 0) || length(seed) != 1L || is.na(seed)) {
    stop("simulate_independent_replicate received invalid inputs.")
  }
  state <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  } else {
    NULL
  }
  on.exit({
    if (!is.null(state)) {
      assign(".Random.seed", state, envir = globalenv())
    }
  }, add = TRUE)
  set.seed(seed)
  noise <- matrix(
    stats::rnorm(length(posterior_means)),
    nrow = nrow(posterior_means),
    ncol = ncol(posterior_means)
  )
  simulated <- posterior_means + standard_errors * noise
  dimnames(simulated) <- dimnames(posterior_means)
  simulated
}

compare_correlation_matrices <- function(reference,
                                         comparison,
                                         reference_name,
                                         comparison_name) {
  reference <- as.matrix(reference)
  comparison <- as.matrix(comparison)
  if (!identical(dim(reference), dim(comparison))) {
    stop("Correlation matrices must have identical dimensions.")
  }
  upper <- upper.tri(reference)
  data.frame(
    reference = reference_name,
    comparison = comparison_name,
    mean_absolute_difference = mean(abs(reference[upper] - comparison[upper])),
    maximum_absolute_difference = max(abs(reference[upper] - comparison[upper])),
    offdiagonal_pearson = stats::cor(reference[upper], comparison[upper]),
    offdiagonal_spearman = stats::cor(
      reference[upper], comparison[upper], method = "spearman"
    ),
    stringsAsFactors = FALSE
  )
}

# Discovery counting deliberately reuses `cumulative_lfdr_calls()` from
# `correlated_likelihood_sensitivity/correlated_likelihood_helpers.R`, which
# returns the selected indices, so callers here wrap it in `length()`. Do not
# redefine it: this file is sourced alongside that one.

# ------------------------------------------------------------------------
# Can this estimator identify the error correlation at all?
#
# The residual-z second moment is linear in the true error correlation P:
#
#   M(P) = (1/n) sum_j T_j P T_j',   vec(M) = K vec(P),
#   K    = (1/n) sum_j (T_j kron T_j).
#
# So the observed moment can in principle be deconvolved back to P. Whether
# that inverse is usable is an empirical question about the spectrum of K
# restricted to symmetric matrices: directions of P that the shrinkage
# operator annihilates cannot be recovered from residuals at any sample size.
# ------------------------------------------------------------------------

symmetric_basis <- function(n_time) {
  n_time <- as.integer(n_time)
  if (length(n_time) != 1L || is.na(n_time) || n_time < 2L) {
    stop("symmetric_basis requires a dimension of at least two.")
  }
  index <- which(upper.tri(diag(n_time), diag = TRUE), arr.ind = TRUE)
  basis <- matrix(0, nrow = n_time * n_time, ncol = nrow(index))
  for (column in seq_len(nrow(index))) {
    element <- matrix(0, n_time, n_time)
    element[index[column, 1], index[column, 2]] <- 1
    element[index[column, 2], index[column, 1]] <- 1
    basis[, column] <- as.numeric(element)
  }
  list(basis = basis, row = index[, 1], column = index[, 2])
}

residual_transfer_kronecker <- function(standard_errors,
                                        posterior_weights,
                                        shared_design,
                                        prior_scale) {
  standard_errors <- as.matrix(standard_errors)
  posterior_weights <- as.matrix(posterior_weights)
  n_unit <- nrow(standard_errors)
  n_time <- ncol(standard_errors)
  if (nrow(posterior_weights) != n_unit ||
      ncol(posterior_weights) != length(prior_scale)) {
    stop("residual_transfer_kronecker received misaligned inputs.")
  }
  kernel <- matrix(0, n_time * n_time, n_time * n_time)
  for (index in seq_len(n_unit)) {
    transfer <- unit_residual_transfer(
      standard_errors = standard_errors[index, ],
      weights = posterior_weights[index, ],
      shared_design = shared_design,
      prior_scale = prior_scale
    )
    kernel <- kernel + kronecker(transfer, transfer)
  }
  kernel / n_unit
}

# Least-squares inverse on the symmetric subspace, with optional ridge
# regularisation. Reports the spectrum so the caller can judge whether the
# inverse means anything.
deconvolve_error_correlation <- function(kernel,
                                        observed_second_moment,
                                        ridge = 0) {
  observed_second_moment <- as.matrix(observed_second_moment)
  n_time <- nrow(observed_second_moment)
  ridge <- as.numeric(ridge)
  if (!identical(dim(kernel), rep(n_time * n_time, 2L)) ||
      ncol(observed_second_moment) != n_time ||
      length(ridge) != 1L || !is.finite(ridge) || ridge < 0) {
    stop("deconvolve_error_correlation received invalid inputs.")
  }
  basis <- symmetric_basis(n_time)
  design <- kernel %*% basis$basis
  spectrum <- svd(design)
  scale_reference <- max(spectrum$d)
  coefficients <- if (ridge > 0) {
    filtered <- spectrum$d / (spectrum$d^2 + (ridge * scale_reference)^2)
    as.numeric(spectrum$v %*% (filtered *
      crossprod(spectrum$u, as.numeric(observed_second_moment))))
  } else {
    as.numeric(qr.solve(design, as.numeric(observed_second_moment)))
  }
  solution <- matrix(0, n_time, n_time)
  solution[cbind(basis$row, basis$column)] <- coefficients
  solution[cbind(basis$column, basis$row)] <- coefficients
  eigenvalues <- eigen(solution, symmetric = TRUE, only.values = TRUE)$values
  scaling <- sqrt(pmax(diag(solution), .Machine$double.eps))
  correlation <- solution / tcrossprod(scaling)
  diag(correlation) <- 1
  list(
    solution = solution,
    correlation = correlation,
    implied_diagonal = diag(solution),
    solution_minimum_eigenvalue = min(eigenvalues),
    singular_values = spectrum$d,
    condition_number = max(spectrum$d) / min(spectrum$d),
    effective_rank = sum(spectrum$d > 1e-8 * scale_reference),
    least_identified_directions = lapply(
      seq(ncol(spectrum$v), max(1L, ncol(spectrum$v) - 2L)),
      function(column) {
        direction <- matrix(0, n_time, n_time)
        direction[cbind(basis$row, basis$column)] <- spectrum$v[, column]
        direction[cbind(basis$column, basis$row)] <- spectrum$v[, column]
        list(singular_value = spectrum$d[column], direction = direction)
      }
    )
  )
}
