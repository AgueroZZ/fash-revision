# Helpers for the mashr-based pair-level null screen.

fit_mashr_pair_screen <- function(beta_hat,
                                  adjusted_se,
                                  z_threshold = 2,
                                  cov_methods = c(
                                    "identity",
                                    "singletons",
                                    "equal_effects",
                                    "simple_het"
                                  ),
                                  prior = "nullbiased",
                                  nullweight = 10,
                                  optmethod = "mixSQP",
                                  seed = 123L,
                                  verbose = TRUE) {
  if (!requireNamespace("mashr", quietly = TRUE)) {
    stop("The mashr package is required for Screen 3.")
  }
  beta_hat <- as.matrix(beta_hat)
  adjusted_se <- as.matrix(adjusted_se)
  z_threshold <- as.numeric(z_threshold)
  nullweight <- as.numeric(nullweight)
  seed <- as.integer(seed)
  if (!identical(dim(beta_hat), dim(adjusted_se)) ||
      nrow(beta_hat) <= ncol(beta_hat) || ncol(beta_hat) < 2L ||
      any(!is.finite(beta_hat)) || any(!is.finite(adjusted_se)) ||
      any(adjusted_se <= 0) || length(z_threshold) != 1L ||
      !is.finite(z_threshold) || z_threshold <= 0 ||
      length(nullweight) != 1L || !is.finite(nullweight) ||
      nullweight <= 0 || length(seed) != 1L || is.na(seed) ||
      length(cov_methods) < 1L || any(!nzchar(cov_methods))) {
    stop("Invalid beta estimates, standard errors, or mashr settings.")
  }
  if (!is.null(rownames(beta_hat)) &&
      !identical(rownames(beta_hat), rownames(adjusted_se))) {
    stop("beta_hat and adjusted_se row names must match.")
  }

  z <- beta_hat / adjusted_se
  maximum_absolute_z <- apply(abs(z), 1L, max)
  nullish_keep <- maximum_absolute_z < z_threshold
  if (sum(nullish_keep) <= ncol(z)) {
    stop("Too few maximum-z null-like rows to estimate null correlation.")
  }

  initial_data <- mashr::mash_set_data(beta_hat, adjusted_se)
  null_correlation <- mashr::estimate_null_correlation_simple(
    initial_data,
    z_thresh = z_threshold,
    est_cor = TRUE
  )
  direct_null_correlation <- stats::cor(z[nullish_keep, , drop = FALSE])
  null_correlation_difference <- max(abs(
    null_correlation - direct_null_correlation
  ))
  if (!is.finite(null_correlation_difference) ||
      null_correlation_difference > 1e-12) {
    stop("mashr's null correlation does not match the documented rule.")
  }

  updated_data <- mashr::mash_update_data(
    initial_data,
    V = null_correlation
  )
  canonical_covariances <- mashr::cov_canonical(
    updated_data,
    cov_methods = cov_methods
  )
  fit <- mashr::mash(
    updated_data,
    Ulist = canonical_covariances,
    usepointmass = TRUE,
    prior = prior,
    nullweight = nullweight,
    optmethod = optmethod,
    verbose = verbose,
    seed = seed,
    outputlevel = 2,
    output_lfdr = TRUE
  )

  if (is.null(fit$posterior_weights) ||
      !"null" %in% colnames(fit$posterior_weights) ||
      nrow(fit$posterior_weights) != nrow(beta_hat) ||
      is.null(fit$result$lfdr) ||
      !identical(dim(fit$result$lfdr), dim(beta_hat))) {
    stop("The mashr fit does not contain the required null probabilities.")
  }
  pair_lfdr <- as.numeric(fit$posterior_weights[, "null"])
  condition_lfdr <- as.matrix(fit$result$lfdr)
  if (any(!is.finite(pair_lfdr)) || any(pair_lfdr < 0 | pair_lfdr > 1) ||
      any(!is.finite(condition_lfdr)) ||
      any(condition_lfdr < 0 | condition_lfdr > 1)) {
    stop("The mashr null probabilities are invalid.")
  }
  names(pair_lfdr) <- rownames(beta_hat)
  rownames(condition_lfdr) <- rownames(beta_hat)
  colnames(condition_lfdr) <- colnames(beta_hat)

  fitted_pi <- fit$fitted_g$pi
  fitted_pi0 <- if (!is.null(names(fitted_pi)) && "null" %in% names(fitted_pi)) {
    unname(fitted_pi["null"])
  } else {
    NA_real_
  }

  list(
    fit = fit,
    null_correlation = null_correlation,
    nullish_keep = nullish_keep,
    maximum_absolute_z = maximum_absolute_z,
    pair_lfdr = pair_lfdr,
    condition_lfdr = condition_lfdr,
    fitted_pi0 = fitted_pi0,
    canonical_covariance_names = names(canonical_covariances),
    diagnostics = list(
      n_pairs = nrow(beta_hat),
      n_conditions = ncol(beta_hat),
      n_nullish_for_correlation = sum(nullish_keep),
      z_threshold = z_threshold,
      null_correlation_maximum_difference = null_correlation_difference,
      minimum_null_correlation_eigenvalue = min(eigen(
        null_correlation,
        symmetric = TRUE,
        only.values = TRUE
      )$values),
      prior = prior,
      nullweight = nullweight,
      optmethod = optmethod,
      seed = seed,
      cov_methods = cov_methods
    )
  )
}
