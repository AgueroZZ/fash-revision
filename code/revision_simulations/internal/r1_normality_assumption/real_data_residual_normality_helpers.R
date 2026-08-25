# Calculate the conventional standardized residuals used in regression diagnostics.

fit_standardized_residuals <- function(expression, design) {
  expression <- as.numeric(expression)
  design <- as.matrix(design)
  storage.mode(design) <- "double"

  n_observations <- length(expression)
  n_parameters <- ncol(design)
  if (nrow(design) != n_observations ||
      n_observations <= n_parameters ||
      any(!is.finite(expression)) ||
      any(!is.finite(design))) {
    return(NULL)
  }

  fit <- stats::lm.fit(x = design, y = expression)
  if (fit$rank != n_parameters) {
    return(NULL)
  }

  residual_df <- n_observations - n_parameters
  residual <- as.numeric(fit$residuals)
  residual_standard_error <- sqrt(sum(residual^2) / residual_df)
  orthonormal_design <- qr.Q(fit$qr, complete = FALSE)
  leverage <- rowSums(orthonormal_design^2)
  residual_scale <- residual_standard_error * sqrt(1 - leverage)
  if (!is.finite(residual_standard_error) ||
      residual_standard_error <= 0 ||
      any(!is.finite(residual_scale)) ||
      any(residual_scale <= 0)) {
    return(NULL)
  }

  list(
    n_observations = as.integer(n_observations),
    n_parameters = as.integer(n_parameters),
    residual_df = as.integer(residual_df),
    fitted = as.numeric(fit$fitted.values),
    residual = residual,
    leverage = as.numeric(leverage),
    residual_standard_error = residual_standard_error,
    standardized_residual = residual / residual_scale
  )
}

calculate_skewness <- function(x) {
  centered <- x - mean(x)
  mean(centered^3) / mean(centered^2)^(3 / 2)
}

calculate_excess_kurtosis <- function(x) {
  centered <- x - mean(x)
  mean(centered^4) / mean(centered^2)^2 - 3
}
