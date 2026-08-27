# Helpers for per-unit cross-time expression correlation by discovery status.

# Residual maker: returns the residualizing projection for a covariate matrix.
make_residualizer <- function(covariates) {
  if (!is.matrix(covariates) || nrow(covariates) <= ncol(covariates) ||
      any(!is.finite(covariates))) {
    stop("Expected a finite covariate matrix with more rows than columns.")
  }
  qr_decomposition <- qr(covariates)
  if (qr_decomposition$rank != ncol(covariates)) {
    stop("The covariate matrix is rank deficient.")
  }
  diag(nrow(covariates)) - qr.fitted(qr_decomposition, diag(nrow(covariates)))
}

# Residualized-genotype OLS weight vector: beta_hat = sum_i weight_i * y_i.
make_genotype_weights <- function(residualizer, genotype) {
  genotype <- as.numeric(genotype)
  if (length(genotype) != nrow(residualizer)) {
    stop("Genotype length does not match the residualizer dimension.")
  }
  genotype_residual <- as.numeric(residualizer %*% genotype)
  denominator <- sum(genotype_residual^2)
  if (!is.finite(denominator) ||
      denominator <= 1e-12 * max(1, sum(genotype^2))) {
    stop("The observed genotype has no residual information.")
  }
  genotype_residual / denominator
}

# Pairwise-complete correlation of per-time donor residual vectors.
# residuals: list of named numeric vectors, one per time.
matched_donor_correlation <- function(residuals) {
  n_time <- length(residuals)
  if (n_time < 2L || any(vapply(residuals, function(r) is.null(names(r)),
                                logical(1L)))) {
    stop("Expected at least two named per-time residual vectors.")
  }
  correlation <- diag(n_time)
  for (one in seq_len(n_time - 1L)) {
    for (two in (one + 1L):n_time) {
      shared <- intersect(names(residuals[[one]]), names(residuals[[two]]))
      if (length(shared) < 4L) {
        stop("Fewer than four shared donors between two times.")
      }
      value <- stats::cor(residuals[[one]][shared], residuals[[two]][shared])
      if (!is.finite(value)) {
        stop("A matched-donor correlation is non-finite.")
      }
      correlation[one, two] <- correlation[two, one] <- value
    }
  }
  correlation
}

# Design factor: cosine similarity of the OLS weight vectors, with the
# numerator over shared donors and each denominator over that time's own
# donors, matching run_design_propagated_null_correlation.R.
weight_design_factor <- function(weights) {
  n_time <- length(weights)
  if (n_time < 2L) {
    stop("Expected at least two per-time weight vectors.")
  }
  norms <- vapply(weights, function(w) sqrt(sum(w^2)), numeric(1L))
  if (any(!is.finite(norms)) || any(norms <= 0)) {
    stop("A weight vector has zero or non-finite norm.")
  }
  factor_matrix <- diag(n_time)
  for (one in seq_len(n_time - 1L)) {
    for (two in (one + 1L):n_time) {
      shared <- intersect(names(weights[[one]]), names(weights[[two]]))
      value <- sum(weights[[one]][shared] * weights[[two]][shared]) /
        (norms[one] * norms[two])
      factor_matrix[one, two] <- factor_matrix[two, one] <- value
    }
  }
  factor_matrix
}

# Mean correlation at each lag of a square matrix.
lag_profile <- function(matrix) {
  if (!is.matrix(matrix) || nrow(matrix) != ncol(matrix) ||
      any(!is.finite(matrix))) {
    stop("Expected a finite square matrix.")
  }
  n_time <- ncol(matrix)
  vapply(seq_len(n_time - 1L), function(lag) {
    mean(matrix[cbind(seq_len(n_time - lag), (lag + 1L):n_time)])
  }, numeric(1L))
}

mean_off_diagonal <- function(matrix) {
  mean(matrix[row(matrix) != col(matrix)])
}
