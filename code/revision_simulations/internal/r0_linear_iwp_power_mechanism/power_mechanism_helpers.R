compute_iwp1_kernel <- function(time_grid, num_basis = 20L) {
  template <- data.frame(
    x = time_grid,
    y = rep(0, length(time_grid)),
    offset = rep(0, length(time_grid))
  )
  tmb_data <- fashr:::fash_set_tmbdat(
    template,
    Si = rep(1, length(time_grid)),
    num_basis = num_basis,
    order = 1
  )
  basis <- as.matrix(tmb_data$B)
  precision <- as.matrix(tmb_data$P)
  basis %*% solve(precision, t(basis))
}

make_clean_standard_error <- function(n_time, pattern) {
  if (pattern == "homoskedastic") {
    standard_error <- rep(1, n_time)
  } else if (pattern == "endpoint_mild") {
    standard_error <- rep(1.05, n_time)
    standard_error[c(1L, n_time)] <- 0.60
  } else if (pattern == "endpoint_strong") {
    standard_error <- rep(1.06, n_time)
    standard_error[c(1L, n_time)] <- 0.30
  } else if (pattern == "edge_quartet") {
    standard_error <- rep(1.10, n_time)
    standard_error[c(1L, 2L, n_time - 1L, n_time)] <- 0.60
  } else if (pattern == "middle_precise") {
    standardized_time <- seq(-1, 1, length.out = n_time)
    standard_error <- 0.55 + 0.90 * abs(standardized_time)
  } else {
    stop("Unknown standard-error pattern.")
  }
  standard_error / sqrt(mean(standard_error^2))
}

compute_common_intercept_log_likelihood <- function(y_matrix,
                                                    standard_error,
                                                    kernel,
                                                    grid) {
  y_matrix <- as.matrix(y_matrix)
  n_time <- ncol(y_matrix)
  standard_error <- as.numeric(standard_error)
  if (!is.matrix(kernel) ||
      any(dim(kernel) != n_time) ||
      !length(standard_error) %in% c(1L, n_time) ||
      any(!is.finite(standard_error)) ||
      any(standard_error <= 0)) {
    stop("The likelihood inputs are invalid.")
  }
  standard_error <- rep(standard_error, length.out = n_time)
  intercept <- rep(1, n_time)
  likelihood <- vapply(grid, function(scale) {
    covariance <- diag(standard_error^2, n_time) + scale^2 * kernel
    chol_covariance <- chol(covariance)
    inverse_covariance <- chol2inv(chol_covariance)
    inverse_intercept <- as.numeric(inverse_covariance %*% intercept)
    intercept_precision <- sum(inverse_intercept)
    residual_precision <- inverse_covariance -
      tcrossprod(inverse_intercept) / intercept_precision
    quadratic <- rowSums((y_matrix %*% residual_precision) * y_matrix)
    log_determinant <- 2 * sum(log(diag(chol_covariance)))
    -0.5 * (
      (n_time - 1) * log(2 * pi) +
        log_determinant +
        log(intercept_precision) +
        quadratic
    )
  }, numeric(nrow(y_matrix)))
  dimnames(likelihood) <- list(rownames(y_matrix), as.character(grid))
  likelihood
}

select_cumulative_lfdr_power_mechanism <- function(lfdr, alpha = 0.05) {
  ordering <- order(lfdr, seq_along(lfdr))
  cumulative_fdr <- cumsum(lfdr[ordering]) / seq_along(ordering)
  last_selected <- max(c(0L, which(cumulative_fdr <= alpha)))
  selected <- rep(FALSE, length(lfdr))
  if (last_selected > 0L) {
    selected[ordering[seq_len(last_selected)]] <- TRUE
  }
  selected
}

summarize_power_mechanism_fit <- function(fit,
                                          true_null,
                                          family,
                                          adjustment,
                                          alpha = 0.05) {
  selected <- select_cumulative_lfdr_power_mechanism(fit$lfdr, alpha)
  discoveries <- sum(selected)
  false_discoveries <- sum(selected & true_null)
  true_positives <- sum(selected & !true_null)
  bf_available <- isTRUE(fit$bf_adjusted) &&
    !is.null(fit$BF) &&
    all(is.finite(fit$BF)) &&
    all(fit$BF > 0)
  data.frame(
    family = family,
    adjustment = adjustment,
    n_discoveries = discoveries,
    false_discoveries = false_discoveries,
    true_positives = true_positives,
    power = true_positives / sum(!true_null),
    realized_fdp = if (discoveries == 0L) 0 else {
      false_discoveries / discoveries
    },
    estimated_pi0 = constant_component_prior_weight(fit),
    bf_available = bf_available,
    median_log_bf_alternative = if (bf_available) {
      stats::median(log(fit$BF[!true_null]))
    } else {
      NA_real_
    },
    median_log_bf_null = if (bf_available) {
      stats::median(log(fit$BF[true_null]))
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

fit_power_mechanism_family <- function(likelihood,
                                       grid,
                                       true_null,
                                       family,
                                       alpha = 0.05) {
  raw <- fit_linear_mixture_fash_from_log_likelihood(
    L_matrix = likelihood,
    grid = grid,
    pred_step = 1,
    penalty = 10
  )
  bf <- tryCatch(
    suppressWarnings(BF_update_linear_mixture_fash(raw)),
    error = function(condition) NULL
  )
  if (is.null(bf) || is.null(bf$BF)) {
    bf <- list(
      lfdr = rep(1, nrow(likelihood)),
      prior_weights = data.frame(psd = 0, prior_weight = 1),
      BF = rep(NA_real_, nrow(likelihood)),
      bf_adjusted = TRUE
    )
  }
  rbind(
    summarize_power_mechanism_fit(
      raw,
      true_null,
      family,
      "Raw",
      alpha
    ),
    summarize_power_mechanism_fit(
      bf,
      true_null,
      family,
      "BF",
      alpha
    )
  )
}

simulate_linear_power_mechanism <- function(J,
                                            pi0,
                                            time_grid,
                                            endpoint_scale,
                                            standard_error,
                                            alternative_distribution = c(
                                              "fixed",
                                              "narrow_normal",
                                              "gaussian"
                                            ),
                                            equicorrelation = 0,
                                            seed) {
  alternative_distribution <- match.arg(alternative_distribution)
  standard_error <- as.numeric(standard_error)
  if (!length(standard_error) %in% c(1L, length(time_grid)) ||
      any(!is.finite(standard_error)) ||
      any(standard_error <= 0)) {
    stop("standard_error must be positive and scalar or time-specific.")
  }
  standard_error <- rep(standard_error, length.out = length(time_grid))
  if (length(equicorrelation) != 1L ||
      !is.finite(equicorrelation) ||
      equicorrelation < 0 ||
      equicorrelation >= 1) {
    stop("equicorrelation must be one value in [0, 1).")
  }
  set.seed(seed)
  n_null <- as.integer(round(J * pi0))
  true_null <- c(rep(TRUE, n_null), rep(FALSE, J - n_null))
  n_alternative <- sum(!true_null)
  intercept <- stats::rnorm(J)
  endpoint <- numeric(J)
  if (alternative_distribution == "fixed") {
    endpoint[!true_null] <- endpoint_scale *
      sample(c(-1, 1), n_alternative, replace = TRUE)
  } else if (alternative_distribution == "narrow_normal") {
    magnitude <- pmax(
      stats::rnorm(n_alternative, mean = endpoint_scale, sd = 0.10),
      0
    )
    endpoint[!true_null] <- magnitude *
      sample(c(-1, 1), n_alternative, replace = TRUE)
  } else {
    endpoint[!true_null] <- stats::rnorm(
      n_alternative,
      mean = 0,
      sd = endpoint_scale
    )
  }
  mean_matrix <- intercept + outer(
    endpoint,
    (time_grid - min(time_grid)) / diff(range(time_grid))
  )
  independent_noise <- matrix(
    stats::rnorm(J * length(time_grid)),
    nrow = J,
    ncol = length(time_grid)
  )
  common_noise <- stats::rnorm(J)
  standardized_noise <- (
    sqrt(1 - equicorrelation) * independent_noise +
      sqrt(equicorrelation) * common_noise
  )
  noise <- sweep(standardized_noise, 2L, standard_error, `*`)
  y_matrix <- mean_matrix + noise
  rownames(y_matrix) <- sprintf("unit_%05d", seq_len(J))
  list(
    y_matrix = y_matrix,
    true_null = true_null,
    endpoint = endpoint
  )
}
