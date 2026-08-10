validate_symmetric_positive_definite <- function(matrix_value,
                                                 name = "matrix_value",
                                                 tolerance = 1e-10) {
  matrix_value <- as.matrix(matrix_value)
  if (!is.numeric(matrix_value) || nrow(matrix_value) != ncol(matrix_value) ||
      any(!is.finite(matrix_value))) {
    stop(name, " must be a finite numeric square matrix.")
  }
  symmetry_error <- max(abs(matrix_value - t(matrix_value)))
  if (symmetry_error > tolerance) {
    stop(name, " must be symmetric; maximum asymmetry = ", symmetry_error)
  }
  chol_result <- try(chol((matrix_value + t(matrix_value)) / 2), silent = TRUE)
  if (inherits(chol_result, "try-error")) {
    stop(name, " must be positive definite.")
  }
  invisible(TRUE)
}

ar1_correlation <- function(size, rho) {
  if (length(size) != 1L || !is.finite(size) || size < 1 || size != floor(size)) {
    stop("size must be a positive integer.")
  }
  if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 1) {
    stop("rho must be a finite scalar strictly between -1 and 1.")
  }
  outer(seq_len(size), seq_len(size), function(i, j) rho^abs(i - j))
}

sigma_iwp_from_psd <- function(psd, order, pred_step) {
  psd / sqrt(
    pred_step^(2 * order - 1) /
      ((2 * order - 1) * factorial(order - 1)^2)
  )
}

model_matrices <- function(data_i, precision, psd, num_basis,
                           betaprec, order, pred_step) {
  if (betaprec <= 0) {
    stop("betaprec must be positive for the proper model-aligned simulation.")
  }
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

make_component_prior_covariances <- function(x, grid, num_basis,
                                             betaprec, order, pred_step) {
  if (length(grid) < 2L || grid[[1]] != 0 || any(!is.finite(grid)) ||
      any(grid[-1] <= 0) || anyDuplicated(grid)) {
    stop("grid must start with zero followed by unique positive PSD values.")
  }
  template <- data.frame(y = numeric(length(x)), x = x, offset = 0)
  identity_precision <- diag(length(x))
  lapply(grid, function(psd) {
    matrices <- model_matrices(
      data_i = template,
      precision = identity_precision,
      psd = psd,
      num_basis = num_basis,
      betaprec = betaprec,
      order = order,
      pred_step = pred_step
    )
    covariance <- matrices$design %*% solve(matrices$prior_precision) %*%
      t(matrices$design)
    covariance <- (covariance + t(covariance)) / 2
    validate_symmetric_positive_definite(
      covariance + diag(1e-12, nrow(covariance)),
      "component prior covariance"
    )
    covariance
  })
}

make_component_marginal_covariances <- function(prior_covariances,
                                                observation_covariance) {
  validate_symmetric_positive_definite(
    observation_covariance,
    "observation_covariance"
  )
  lapply(prior_covariances, function(prior_covariance) {
    covariance <- observation_covariance + prior_covariance
    covariance <- (covariance + t(covariance)) / 2
    validate_symmetric_positive_definite(covariance, "marginal covariance")
    covariance
  })
}

log_mvn_rows <- function(response, covariance) {
  response <- as.matrix(response)
  covariance <- as.matrix(covariance)
  if (ncol(response) != nrow(covariance)) {
    stop("response columns must match covariance dimensions.")
  }
  validate_symmetric_positive_definite(covariance, "covariance")
  chol_covariance <- chol(covariance)
  whitened <- t(forwardsolve(t(chol_covariance), t(response)))
  -0.5 * (
    ncol(response) * log(2 * pi) +
      2 * sum(log(diag(chol_covariance))) +
      rowSums(whitened^2)
  )
}

component_log_likelihood_matrix <- function(response, component_covariances) {
  response <- as.matrix(response)
  if (!is.list(component_covariances) || length(component_covariances) < 2L) {
    stop("component_covariances must be a list with at least two matrices.")
  }
  result <- vapply(
    component_covariances,
    function(covariance) log_mvn_rows(response, covariance),
    numeric(nrow(response))
  )
  if (is.null(dim(result))) {
    result <- matrix(result, ncol = length(component_covariances))
  }
  if (any(!is.finite(result))) {
    stop("The component log-likelihood matrix contains non-finite values.")
  }
  result
}

simulate_marginal_mixture <- function(seed, n_units, true_weights,
                                      component_covariances) {
  if (length(true_weights) != length(component_covariances) ||
      any(!is.finite(true_weights)) || any(true_weights <= 0) ||
      abs(sum(true_weights) - 1) > 1e-10) {
    stop("true_weights must be positive, sum to one, and match the components.")
  }
  if (length(n_units) != 1L || n_units < 1L || n_units != floor(n_units)) {
    stop("n_units must be a positive integer.")
  }
  set.seed(seed)
  component <- sample.int(
    length(true_weights),
    size = n_units,
    replace = TRUE,
    prob = true_weights
  )
  n_time <- nrow(component_covariances[[1]])
  response <- matrix(NA_real_, nrow = n_units, ncol = n_time)
  for (component_index in seq_along(component_covariances)) {
    selected <- which(component == component_index)
    if (length(selected) == 0L) {
      next
    }
    factor <- chol(component_covariances[[component_index]])
    response[selected, ] <- matrix(
      rnorm(length(selected) * n_time),
      nrow = length(selected),
      ncol = n_time
    ) %*% factor
  }
  if (any(!is.finite(response))) {
    stop("Simulation produced non-finite responses.")
  }
  list(response = response, component = component)
}

expand_prior_weights <- function(fit, grid) {
  result <- setNames(numeric(length(grid)), as.character(grid))
  if (is.null(fit$prior_weights) ||
      !all(c("psd", "prior_weight") %in% names(fit$prior_weights))) {
    stop("fit does not contain valid prior_weights.")
  }
  indices <- match(as.character(fit$prior_weights$psd), names(result))
  if (anyNA(indices)) {
    stop("fit prior weights contain PSD values outside the requested grid.")
  }
  result[indices] <- fit$prior_weights$prior_weight
  if (any(!is.finite(result)) || any(result < 0) ||
      abs(sum(result) - 1) > 1e-6) {
    stop("Expanded prior weights are invalid.")
  }
  result
}

expand_posterior_weights <- function(fit, grid) {
  if (is.null(fit$posterior_weights) || is.null(colnames(fit$posterior_weights))) {
    stop("fit does not contain named posterior weights.")
  }
  result <- matrix(
    0,
    nrow = nrow(fit$posterior_weights),
    ncol = length(grid),
    dimnames = list(rownames(fit$posterior_weights), as.character(grid))
  )
  indices <- match(colnames(fit$posterior_weights), colnames(result))
  if (anyNA(indices)) {
    stop("fit posterior weights contain PSD values outside the requested grid.")
  }
  result[, indices] <- fit$posterior_weights
  if (any(!is.finite(result)) || any(result < -1e-12) ||
      max(abs(rowSums(result) - 1)) > 1e-6) {
    stop("Expanded posterior weights are invalid.")
  }
  result
}

fit_raw_from_likelihood <- function(log_likelihood, grid, penalty = 1) {
  log_likelihood <- as.matrix(log_likelihood)
  if (ncol(log_likelihood) != length(grid) || any(!is.finite(log_likelihood))) {
    stop("log_likelihood must be finite and have one column per grid value.")
  }
  rownames(log_likelihood) <- sprintf("unit_%06d", seq_len(nrow(log_likelihood)))
  empirical_bayes <- fashr::fash_eb_est(
    L_matrix = log_likelihood,
    penalty = penalty,
    grid = grid
  )
  rownames(empirical_bayes$posterior_weight) <- rownames(log_likelihood)
  lfdr <- if ("0" %in% colnames(empirical_bayes$posterior_weight)) {
    empirical_bayes$posterior_weight[, "0"]
  } else {
    rep(0, nrow(log_likelihood))
  }
  names(lfdr) <- rownames(log_likelihood)
  fit <- structure(
    list(
      prior_weights = empirical_bayes$prior_weight,
      posterior_weights = empirical_bayes$posterior_weight,
      psd_grid = grid,
      lfdr = lfdr,
      L_matrix = log_likelihood,
      eb_result = empirical_bayes
    ),
    class = "fash"
  )
  expand_prior_weights(fit, grid)
  expand_posterior_weights(fit, grid)
  fit
}

fit_oracle_from_likelihood <- function(log_likelihood, grid, true_weights) {
  log_likelihood <- as.matrix(log_likelihood)
  if (length(true_weights) != length(grid) || any(true_weights <= 0) ||
      abs(sum(true_weights) - 1) > 1e-10) {
    stop("true_weights must be positive, sum to one, and match grid.")
  }
  log_weighted <- sweep(log_likelihood, 2L, log(true_weights), `+`)
  row_maximum <- apply(log_weighted, 1L, max)
  posterior <- exp(log_weighted - row_maximum)
  posterior <- posterior / rowSums(posterior)
  colnames(posterior) <- as.character(grid)
  rownames(posterior) <- rownames(log_likelihood)
  lfdr <- posterior[, 1L]
  names(lfdr) <- rownames(log_likelihood)
  structure(
    list(
      prior_weights = data.frame(psd = grid, prior_weight = true_weights),
      posterior_weights = posterior,
      psd_grid = grid,
      lfdr = lfdr,
      L_matrix = log_likelihood,
      oracle = TRUE
    ),
    class = "fash"
  )
}

fit_bf_from_raw <- function(raw_fit, grid) {
  warning_messages <- character()
  bf_fit <- withCallingHandlers(
    fashr::BF_update(raw_fit, plot = FALSE),
    warning = function(warning_condition) {
      warning_messages <<- c(warning_messages, conditionMessage(warning_condition))
      invokeRestart("muffleWarning")
    }
  )
  if (is.null(bf_fit$BF) || any(!is.finite(bf_fit$BF))) {
    stop(
      "BF_update did not produce finite Bayes factors.",
      if (length(warning_messages) > 0L) paste0(" Warnings: ", paste(warning_messages, collapse = " | ")) else ""
    )
  }
  expand_prior_weights(bf_fit, grid)
  posterior <- expand_posterior_weights(bf_fit, grid)
  if (max(abs(bf_fit$lfdr - posterior[, 1L])) > 1e-10) {
    stop("BF-adjusted lfdr does not match the null posterior weight.")
  }
  attr(bf_fit, "bf_warnings") <- warning_messages
  bf_fit
}

select_lfdr_prefix <- function(lfdr, alpha) {
  if (any(!is.finite(lfdr)) || any(lfdr < 0) || any(lfdr > 1)) {
    stop("lfdr must contain finite probabilities.")
  }
  if (length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("alpha must be a finite scalar in [0, 1].")
  }
  ordering <- order(lfdr)
  cumulative_fdr <- cumsum(lfdr[ordering]) / seq_along(lfdr)
  eligible <- which(cumulative_fdr <= alpha)
  selected_count <- if (length(eligible) == 0L) 0L else max(eligible)
  selected <- if (selected_count == 0L) integer() else ordering[seq_len(selected_count)]
  list(
    selected = selected,
    selected_count = selected_count,
    reported_fdr = if (selected_count == 0L) 0 else cumulative_fdr[[selected_count]],
    ordering = ordering,
    cumulative_fdr = cumulative_fdr
  )
}

evaluate_fdr_fit <- function(fit, truth_component, alpha_grid, seed, n_units,
                             observation_model, correction) {
  if (length(truth_component) != length(fit$lfdr)) {
    stop("truth_component must align with fit$lfdr.")
  }
  true_null <- truth_component == 1L
  true_alternative_count <- sum(!true_null)
  rows <- lapply(alpha_grid, function(alpha) {
    selection <- select_lfdr_prefix(fit$lfdr, alpha)
    selected <- selection$selected
    false_discoveries <- if (length(selected) == 0L) 0L else sum(true_null[selected])
    true_positives <- if (length(selected) == 0L) 0L else sum(!true_null[selected])
    discoveries <- length(selected)
    data.frame(
      seed = seed,
      n_units = n_units,
      observation_model = observation_model,
      correction = correction,
      method = paste(observation_model, correction, sep = "-"),
      alpha = alpha,
      discoveries = discoveries,
      false_discoveries = false_discoveries,
      true_positives = true_positives,
      fdp = if (discoveries == 0L) 0 else false_discoveries / discoveries,
      power = if (true_alternative_count == 0L) NA_real_ else true_positives / true_alternative_count,
      reported_fdr = selection$reported_fdr,
      calibration_gap = if (discoveries == 0L) 0 else
        false_discoveries / discoveries - selection$reported_fdr,
      no_discoveries = discoveries == 0L,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

prior_weight_rows <- function(fit, truth_component, true_weights, grid,
                              seed, n_units, observation_model, correction) {
  estimated <- expand_prior_weights(fit, grid)
  realized <- tabulate(truth_component, nbins = length(grid)) / length(truth_component)
  data.frame(
    seed = seed,
    n_units = n_units,
    observation_model = observation_model,
    correction = correction,
    method = paste(observation_model, correction, sep = "-"),
    psd = grid,
    true_weight = true_weights,
    realized_weight = realized,
    estimated_weight = as.numeric(estimated),
    population_error = as.numeric(estimated) - true_weights,
    realized_error = as.numeric(estimated) - realized,
    stringsAsFactors = FALSE
  )
}

bf_diagnostic_row <- function(fit, true_weights, grid, seed, n_units,
                              observation_model) {
  estimated <- expand_prior_weights(fit, grid)
  estimated_alt_total <- sum(estimated[-1L])
  conditional_alt <- if (estimated_alt_total > 0) {
    estimated[-1L] / estimated_alt_total
  } else {
    rep(NA_real_, length(grid) - 1L)
  }
  true_conditional_alt <- true_weights[-1L] / sum(true_weights[-1L])
  data.frame(
    seed = seed,
    n_units = n_units,
    observation_model = observation_model,
    true_pi0 = true_weights[[1]],
    estimated_pi0 = estimated[[1]],
    pi0_error = estimated[[1]] - true_weights[[1]],
    conservative_pi0 = estimated[[1]] >= true_weights[[1]],
    conditional_alt_l1 = if (anyNA(conditional_alt)) NA_real_ else
      sum(abs(conditional_alt - true_conditional_alt)),
    min_bf = min(fit$BF),
    median_bf = stats::median(fit$BF),
    max_bf = max(fit$BF),
    stringsAsFactors = FALSE
  )
}

summarize_numeric <- function(values, confidence_level = 0.95) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  n <- length(values)
  if (n == 0L) {
    return(c(n = 0, mean = NA, sd = NA, se = NA, lower = NA, upper = NA))
  }
  mean_value <- mean(values)
  if (n == 1L) {
    return(c(n = 1, mean = mean_value, sd = NA, se = NA, lower = NA, upper = NA))
  }
  sd_value <- stats::sd(values)
  se_value <- sd_value / sqrt(n)
  critical <- stats::qt((1 + confidence_level) / 2, df = n - 1L)
  c(
    n = n,
    mean = mean_value,
    sd = sd_value,
    se = se_value,
    lower = mean_value - critical * se_value,
    upper = mean_value + critical * se_value
  )
}

summarize_prior_accuracy <- function(prior_rows) {
  key <- interaction(
    prior_rows$observation_model,
    prior_rows$correction,
    prior_rows$n_units,
    prior_rows$psd,
    drop = TRUE
  )
  component_rows <- lapply(split(prior_rows, key), function(group) {
    error <- group$population_error
    data.frame(
      summary_level = "component",
      observation_model = group$observation_model[[1]],
      correction = group$correction[[1]],
      method = group$method[[1]],
      n_units = group$n_units[[1]],
      psd = group$psd[[1]],
      n_replications = length(unique(group$seed)),
      bias = mean(error),
      sd = stats::sd(error),
      rmse = sqrt(mean(error^2)),
      mean_absolute_error = mean(abs(error)),
      mean_tv = NA_real_,
      tv_mc_se = NA_real_,
      tv_ci_lower = NA_real_,
      tv_ci_upper = NA_real_,
      stringsAsFactors = FALSE
    )
  })

  seed_key <- interaction(
    prior_rows$seed,
    prior_rows$observation_model,
    prior_rows$correction,
    prior_rows$n_units,
    drop = TRUE
  )
  vector_by_seed <- lapply(split(prior_rows, seed_key), function(group) {
    data.frame(
      seed = group$seed[[1]],
      observation_model = group$observation_model[[1]],
      correction = group$correction[[1]],
      method = group$method[[1]],
      n_units = group$n_units[[1]],
      tv = 0.5 * sum(abs(group$population_error)),
      max_absolute_error = max(abs(group$population_error)),
      stringsAsFactors = FALSE
    )
  })
  vector_by_seed <- do.call(rbind, vector_by_seed)
  vector_key <- interaction(
    vector_by_seed$observation_model,
    vector_by_seed$correction,
    vector_by_seed$n_units,
    drop = TRUE
  )
  vector_rows <- lapply(split(vector_by_seed, vector_key), function(group) {
    tv_summary <- summarize_numeric(group$tv)
    data.frame(
      summary_level = "vector",
      observation_model = group$observation_model[[1]],
      correction = group$correction[[1]],
      method = group$method[[1]],
      n_units = group$n_units[[1]],
      psd = NA_real_,
      n_replications = nrow(group),
      bias = mean(group$max_absolute_error),
      sd = stats::sd(group$max_absolute_error),
      rmse = sqrt(mean(group$max_absolute_error^2)),
      mean_absolute_error = mean(group$max_absolute_error),
      mean_tv = tv_summary[["mean"]],
      tv_mc_se = tv_summary[["se"]],
      tv_ci_lower = pmax(0, tv_summary[["lower"]]),
      tv_ci_upper = tv_summary[["upper"]],
      stringsAsFactors = FALSE
    )
  })
  list(
    summary = do.call(rbind, c(component_rows, vector_rows)),
    vector_by_seed = vector_by_seed
  )
}

summarize_fdr_replications <- function(fdr_rows) {
  key <- interaction(
    fdr_rows$observation_model,
    fdr_rows$correction,
    fdr_rows$n_units,
    fdr_rows$alpha,
    drop = TRUE
  )
  rows <- lapply(split(fdr_rows, key), function(group) {
    fdr <- summarize_numeric(group$fdp)
    power <- summarize_numeric(group$power)
    gap <- summarize_numeric(group$calibration_gap)
    data.frame(
      observation_model = group$observation_model[[1]],
      correction = group$correction[[1]],
      method = group$method[[1]],
      n_units = group$n_units[[1]],
      alpha = group$alpha[[1]],
      n_replications = nrow(group),
      mean_fdr = fdr[["mean"]],
      fdr_mc_se = fdr[["se"]],
      fdr_ci_lower = pmax(0, fdr[["lower"]]),
      fdr_ci_upper = pmin(1, fdr[["upper"]]),
      mean_power = power[["mean"]],
      power_mc_se = power[["se"]],
      power_ci_lower = pmax(0, power[["lower"]]),
      power_ci_upper = pmin(1, power[["upper"]]),
      mean_discoveries = mean(group$discoveries),
      mean_false_discoveries = mean(group$false_discoveries),
      probability_no_discoveries = mean(group$no_discoveries),
      mean_reported_fdr = mean(group$reported_fdr),
      mean_calibration_gap = gap[["mean"]],
      calibration_gap_mc_se = gap[["se"]],
      calibration_gap_ci_lower = gap[["lower"]],
      calibration_gap_ci_upper = gap[["upper"]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_bf_diagnostics <- function(bf_rows) {
  key <- interaction(
    bf_rows$observation_model,
    bf_rows$n_units,
    drop = TRUE
  )
  rows <- lapply(split(bf_rows, key), function(group) {
    pi0 <- summarize_numeric(group$estimated_pi0)
    error <- summarize_numeric(group$pi0_error)
    data.frame(
      observation_model = group$observation_model[[1]],
      n_units = group$n_units[[1]],
      n_replications = nrow(group),
      true_pi0 = group$true_pi0[[1]],
      mean_estimated_pi0 = pi0[["mean"]],
      pi0_mc_se = pi0[["se"]],
      pi0_ci_lower = pmax(0, pi0[["lower"]]),
      pi0_ci_upper = pmin(1, pi0[["upper"]]),
      mean_pi0_error = error[["mean"]],
      probability_conservative_pi0 = mean(group$conservative_pi0),
      mean_conditional_alt_l1 = mean(group$conditional_alt_l1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

bootstrap_mean_ci <- function(values, n_bootstrap = 2000L, seed = 1L,
                              confidence_level = 0.95) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  if (length(values) < 2L) {
    return(c(mean = mean(values), lower = NA_real_, upper = NA_real_))
  }
  set.seed(seed)
  bootstrap_means <- replicate(
    n_bootstrap,
    mean(sample(values, size = length(values), replace = TRUE))
  )
  tail_probability <- (1 - confidence_level) / 2
  c(
    mean = mean(values),
    lower = unname(stats::quantile(bootstrap_means, tail_probability)),
    upper = unname(stats::quantile(bootstrap_means, 1 - tail_probability))
  )
}
