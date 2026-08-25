# Scientific helpers for the profiled FASH-linear real-data ablation.

linear_log_sum_exp_pair <- function(a, b) {
  maximum <- pmax(a, b)
  maximum + log(exp(a - maximum) + exp(b - maximum))
}

linear_fit_two_component_pi <- function(log_m0, log_m1, eps = 1e-6) {
  if (length(log_m0) != length(log_m1) || length(log_m0) == 0L) {
    stop("log_m0 and log_m1 must have the same positive length.")
  }
  if (any(!is.finite(log_m0)) || any(!is.finite(log_m1))) {
    stop("The marginal log likelihoods must be finite.")
  }

  loglik_at <- function(pi0) {
    sum(linear_log_sum_exp_pair(
      log(pi0) + log_m0,
      log1p(-pi0) + log_m1
    ))
  }
  optimum <- stats::optimize(
    function(pi0) -loglik_at(pi0),
    interval = c(eps, 1 - eps)
  )
  candidates <- c(eps, optimum$minimum, 1 - eps)
  loglik <- vapply(candidates, loglik_at, numeric(1))
  best <- which.max(loglik)
  list(pi0 = candidates[best], loglik = loglik[best])
}

compute_linear_sufficient_statistics_block <- function(datasets,
                                                       expected_time = 0:15,
                                                       time_column = "time",
                                                       response_column = "beta",
                                                       se_column = "SE",
                                                       scale_time = TRUE,
                                                       ridge = 1e-10) {
  if (!is.list(datasets) || length(datasets) == 0L) {
    stop("datasets must be a nonempty list.")
  }
  if (is.null(names(datasets)) || any(names(datasets) == "")) {
    stop("Every dataset must have a nonempty unit name.")
  }
  expected_time <- as.numeric(expected_time)
  if (length(expected_time) < 2L || any(!is.finite(expected_time))) {
    stop("expected_time must contain at least two finite values.")
  }
  scaled_time <- if (scale_time) {
    (expected_time - min(expected_time)) / diff(range(expected_time))
  } else {
    expected_time
  }

  statistics <- vapply(datasets, function(dataset) {
    required_columns <- c(time_column, response_column, se_column)
    if (!is.data.frame(dataset) ||
        !all(required_columns %in% names(dataset)) ||
        nrow(dataset) != length(expected_time)) {
      stop("A dataset has invalid columns or time dimension.")
    }
    observed_time <- as.numeric(dataset[[time_column]])
    y <- as.numeric(dataset[[response_column]])
    se <- as.numeric(dataset[[se_column]])
    if (!isTRUE(all.equal(observed_time, expected_time, tolerance = 0)) ||
        any(!is.finite(y)) || any(!is.finite(se)) || any(se <= 0)) {
      stop("A dataset has invalid time, response, or standard-error values.")
    }
    variance <- se^2 + ridge
    weight <- 1 / variance
    c(
      sum_w = sum(weight),
      sum_wx = sum(weight * scaled_time),
      sum_wxx = sum(weight * scaled_time^2),
      sum_wy = sum(weight * y),
      sum_wxy = sum(weight * scaled_time * y),
      sum_wyy = sum(weight * y^2),
      logdet_d = sum(log(variance))
    )
  }, numeric(7))

  result <- data.frame(
    unit_id = names(datasets),
    t(statistics),
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  numeric_columns <- setdiff(names(result), "unit_id")
  result[numeric_columns] <- lapply(result[numeric_columns], as.numeric)
  result
}

linear_log_marginal_from_stats <- function(statistics,
                                           slope_sd,
                                           n_time = 16L) {
  required_columns <- c(
    "sum_w", "sum_wx", "sum_wxx", "sum_wy", "sum_wxy",
    "sum_wyy", "logdet_d"
  )
  if (!is.data.frame(statistics) ||
      !all(required_columns %in% names(statistics)) ||
      nrow(statistics) == 0L) {
    stop("statistics is missing required sufficient-statistic columns.")
  }
  if (length(slope_sd) != 1L || !is.finite(slope_sd) || slope_sd < 0) {
    stop("slope_sd must be one finite nonnegative value.")
  }
  tau2 <- slope_sd^2
  denominator <- 1 + tau2 * statistics$sum_wxx
  intercept_precision <- statistics$sum_w -
    tau2 * statistics$sum_wx^2 / denominator
  intercept_response <- statistics$sum_wy -
    tau2 * statistics$sum_wx * statistics$sum_wxy / denominator
  response_precision <- statistics$sum_wyy -
    tau2 * statistics$sum_wxy^2 / denominator
  if (any(!is.finite(denominator)) || any(denominator <= 0) ||
      any(!is.finite(intercept_precision)) ||
      any(intercept_precision <= 0)) {
    stop("The rank-one covariance update produced invalid precision terms.")
  }
  logdet_v <- statistics$logdet_d + log(denominator)
  quadratic <- response_precision -
    intercept_response^2 / intercept_precision
  result <- -0.5 * (
    (as.integer(n_time) - 1L) * log(2 * pi) +
      logdet_v + log(intercept_precision) + quadratic
  )
  if (any(!is.finite(result))) {
    stop("The vectorized marginal log likelihood is not finite.")
  }
  result
}

compute_linear_mixture_log_likelihood_from_stats <- function(
    statistics,
    grid = default_revision_grid(),
    pred_step = 1,
    statistic_time_span = 15,
    n_time = 16L) {
  grid <- validate_linear_mixture_grid(grid)
  if (length(pred_step) != 1L ||
      !is.finite(pred_step) ||
      pred_step <= 0) {
    stop("pred_step must be one positive finite value.")
  }
  if (length(statistic_time_span) != 1L ||
      !is.finite(statistic_time_span) ||
      statistic_time_span <= 0) {
    stop("statistic_time_span must be one positive finite value.")
  }
  if (length(n_time) != 1L ||
      !is.finite(n_time) ||
      n_time < 2L ||
      n_time != as.integer(n_time)) {
    stop("n_time must be one integer greater than one.")
  }
  required_columns <- c(
    "unit_id", "sum_w", "sum_wx", "sum_wxx", "sum_wy", "sum_wxy",
    "sum_wyy", "logdet_d"
  )
  if (!is.data.frame(statistics) ||
      !all(required_columns %in% names(statistics)) ||
      nrow(statistics) < 1L ||
      anyNA(statistics$unit_id) ||
      any(statistics$unit_id == "") ||
      anyDuplicated(statistics$unit_id)) {
    stop("statistics must contain unique units and all required columns.")
  }
  numeric_columns <- setdiff(required_columns, "unit_id")
  if (any(!vapply(statistics[numeric_columns], is.numeric, logical(1))) ||
      any(!is.finite(as.matrix(statistics[numeric_columns])))) {
    stop("The sufficient-statistic columns must be finite and numeric.")
  }

  likelihood <- vapply(grid, function(predstep_sd) {
    # The cached moments use endpoint-scaled time. Convert the predictive-step
    # slope SD back to that coefficient scale before evaluating the likelihood.
    cached_slope_sd <- predstep_sd * statistic_time_span / pred_step
    linear_log_marginal_from_stats(
      statistics = statistics,
      slope_sd = cached_slope_sd,
      n_time = as.integer(n_time)
    )
  }, numeric(nrow(statistics)))
  dimnames(likelihood) <- list(
    as.character(statistics$unit_id),
    as.character(grid)
  )
  if (any(!is.finite(likelihood))) {
    stop("The sufficient-statistic mixture likelihood is not finite.")
  }
  likelihood
}

fit_linear_mixture_fash_from_stats <- function(
    statistics,
    grid = default_revision_grid(),
    pred_step = 1,
    penalty = 10,
    statistic_time_span = 15,
    n_time = 16L) {
  L_matrix <- compute_linear_mixture_log_likelihood_from_stats(
    statistics = statistics,
    grid = grid,
    pred_step = pred_step,
    statistic_time_span = statistic_time_span,
    n_time = n_time
  )
  fit_linear_mixture_fash_from_log_likelihood(
    L_matrix = L_matrix,
    grid = grid,
    pred_step = pred_step,
    penalty = penalty
  )
}

fit_profiled_linear_fash_from_stats <- function(statistics,
                                                sigma_beta_grid = exp(seq(
                                                  log(0.05),
                                                  log(5),
                                                  length.out = 25
                                                )),
                                                n_time = 16L) {
  sigma_beta_grid <- sort(unique(as.numeric(sigma_beta_grid)))
  if (length(sigma_beta_grid) == 0L ||
      any(!is.finite(sigma_beta_grid)) ||
      any(sigma_beta_grid <= 0)) {
    stop("sigma_beta_grid must contain finite positive values.")
  }
  log_m0 <- linear_log_marginal_from_stats(
    statistics,
    slope_sd = 0,
    n_time = n_time
  )
  best <- NULL
  profile_rows <- vector("list", length(sigma_beta_grid))
  for (index in seq_along(sigma_beta_grid)) {
    sigma_beta <- sigma_beta_grid[index]
    log_m1 <- linear_log_marginal_from_stats(
      statistics,
      slope_sd = sigma_beta,
      n_time = n_time
    )
    pi_fit <- linear_fit_two_component_pi(log_m0, log_m1)
    profile_rows[[index]] <- data.frame(
      sigma_beta = sigma_beta,
      estimated_pi0 = pi_fit$pi0,
      loglik = pi_fit$loglik,
      stringsAsFactors = FALSE
    )
    if (is.null(best) || pi_fit$loglik > best$loglik) {
      best <- list(
        index = index,
        sigma_beta = sigma_beta,
        pi0 = pi_fit$pi0,
        loglik = pi_fit$loglik,
        log_m1 = log_m1
      )
    }
  }
  profile <- do.call(rbind, profile_rows)
  profile$selected <- seq_len(nrow(profile)) == best$index
  profile$grid_boundary <- seq_len(nrow(profile)) %in%
    c(1L, nrow(profile))
  log_denominator <- linear_log_sum_exp_pair(
    log(best$pi0) + log_m0,
    log1p(-best$pi0) + best$log_m1
  )
  lfdr <- exp(log(best$pi0) + log_m0 - log_denominator)

  structure(list(
    sigma_beta = best$sigma_beta,
    sigma_beta_candidates = sigma_beta_grid,
    sigma_profile = profile,
    selected_sigma_on_boundary = profile$grid_boundary[best$index],
    prior_weights = data.frame(
      component = c("constant", "linear_slope"),
      prior_weight = c(best$pi0, 1 - best$pi0),
      stringsAsFactors = FALSE
    ),
    log_marginal = data.frame(
      unit_id = statistics$unit_id,
      log_m0 = log_m0,
      log_m1 = best$log_m1,
      stringsAsFactors = FALSE
    ),
    lfdr = lfdr,
    posterior_weights = cbind(null = lfdr, alternative = 1 - lfdr),
    loglik = best$loglik,
    bf_adjusted = FALSE
  ), class = "profiled_linear_fash")
}

bf_update_profiled_linear_fash <- function(fit) {
  if (!inherits(fit, "profiled_linear_fash")) {
    stop("fit must inherit from profiled_linear_fash.")
  }
  log_bf <- fit$log_marginal$log_m1 - fit$log_marginal$log_m0
  bayes_factor <- exp(log_bf)
  ordering <- order(bayes_factor, na.last = NA)
  sorted_bf <- bayes_factor[ordering]
  cumulative_mean <- cumsum(sorted_bf) / seq_along(sorted_bf)
  estimated_pi0 <- seq_along(sorted_bf) / length(sorted_bf)
  pi0_bf <- if (max(cumulative_mean, na.rm = TRUE) < 1) {
    1
  } else {
    estimated_pi0[which(cumulative_mean >= 1)[1]]
  }
  if (!is.finite(pi0_bf) || pi0_bf < 0 || pi0_bf > 1) {
    stop("The BF-adjusted null weight is invalid.")
  }
  log_denominator <- linear_log_sum_exp_pair(
    log(pi0_bf) + fit$log_marginal$log_m0,
    log1p(-pi0_bf) + fit$log_marginal$log_m1
  )
  lfdr <- exp(
    log(pi0_bf) + fit$log_marginal$log_m0 - log_denominator
  )
  fit$prior_weights <- data.frame(
    component = c("constant", "linear_slope"),
    prior_weight = c(pi0_bf, 1 - pi0_bf),
    stringsAsFactors = FALSE
  )
  fit$posterior_weights <- cbind(null = lfdr, alternative = 1 - lfdr)
  fit$lfdr <- lfdr
  fit$BF <- bayes_factor
  fit$bf_adjusted <- TRUE
  fit
}

select_cumulative_lfdr_calls_linear <- function(lfdr, alpha = 0.05) {
  lfdr <- as.numeric(lfdr)
  if (length(lfdr) == 0L || any(!is.finite(lfdr)) ||
      any(lfdr < 0 | lfdr > 1)) {
    stop("lfdr must contain finite values in [0, 1].")
  }
  ordering <- order(lfdr, seq_along(lfdr))
  cumulative_fdr <- cumsum(lfdr[ordering]) / seq_along(ordering)
  last_call <- max(c(0L, which(cumulative_fdr <= alpha)))
  if (last_call == 0L) integer() else sort(ordering[seq_len(last_call)])
}

parse_linear_pair_keys <- function(keys) {
  separator <- regexpr("_", keys, fixed = TRUE)
  if (any(separator < 2L)) {
    stop("Every pair key must contain a gene and variant separator.")
  }
  data.frame(
    key = keys,
    gene_id = substr(keys, 1L, separator - 1L),
    variant_id = substr(keys, separator + 1L, nchar(keys)),
    stringsAsFactors = FALSE
  )
}

select_current_top_pair_per_gene <- function(pair_table, current_lfdr) {
  if (nrow(pair_table) != length(current_lfdr)) {
    stop("pair_table and current_lfdr are not aligned.")
  }
  ordering <- order(
    pair_table$gene_id,
    current_lfdr,
    pair_table$variant_id,
    pair_table$key
  )
  ordered_genes <- pair_table$gene_id[ordering]
  ordering[!duplicated(ordered_genes)]
}

extract_profiled_linear_posterior_plot_data <- function(dataset,
                                                        standard_error,
                                                        slope_sd,
                                                        lfdr,
                                                        grid = seq(
                                                          0,
                                                          15,
                                                          by = 0.1
                                                        ),
                                                        sample_size = 10000L,
                                                        seed = 20260810L,
                                                        ridge = 1e-10) {
  if (!is.data.frame(dataset) || !all(c("x", "y") %in% names(dataset))) {
    stop("dataset must contain x and y columns.")
  }
  x <- as.numeric(dataset$x)
  y <- as.numeric(dataset$y)
  standard_error <- as.numeric(standard_error)
  if (length(x) != length(y) || length(y) != length(standard_error) ||
      any(!is.finite(x)) || any(!is.finite(y)) ||
      any(!is.finite(standard_error)) || any(standard_error <= 0)) {
    stop("The observed trajectory or standard errors are invalid.")
  }
  if (!is.finite(slope_sd) || slope_sd <= 0 ||
      !is.finite(lfdr) || lfdr < 0 || lfdr > 1 ||
      sample_size < 100L) {
    stop("The slope prior or posterior sampling settings are invalid.")
  }
  scaled_x <- (x - min(x)) / diff(range(x))
  scaled_grid <- (grid - min(x)) / diff(range(x))
  weight <- 1 / (standard_error^2 + ridge)

  null_variance <- 1 / sum(weight)
  null_mean <- sum(weight * y) * null_variance

  design <- cbind(1, scaled_x)
  posterior_precision <- crossprod(design, weight * design) +
    diag(c(0, 1 / slope_sd^2), nrow = 2L)
  posterior_covariance <- solve(posterior_precision)
  posterior_mean <- as.numeric(
    posterior_covariance %*% crossprod(design, weight * y)
  )

  set.seed(seed)
  alternative <- stats::runif(sample_size) > lfdr
  intercept_draw <- stats::rnorm(
    sample_size,
    mean = null_mean,
    sd = sqrt(null_variance)
  )
  slope_draw <- numeric(sample_size)
  alternative_count <- sum(alternative)
  if (alternative_count > 0L) {
    standard_draws <- matrix(
      stats::rnorm(2L * alternative_count),
      ncol = 2L
    )
    coefficient_draws <- sweep(
      standard_draws %*% chol(posterior_covariance),
      2L,
      posterior_mean,
      `+`
    )
    intercept_draw[alternative] <- coefficient_draws[, 1L]
    slope_draw[alternative] <- coefficient_draws[, 2L]
  }
  trajectory_draws <- outer(intercept_draw, rep(1, length(grid))) +
    outer(slope_draw, scaled_grid)
  intervals <- apply(
    trajectory_draws,
    2L,
    stats::quantile,
    probs = c(0.025, 0.975),
    names = FALSE
  )
  data.frame(
    time = as.numeric(grid),
    posterior_mean = colMeans(trajectory_draws),
    lower = intervals[1L, ],
    upper = intervals[2L, ],
    stringsAsFactors = FALSE
  )
}

extract_linear_mixture_posterior_plot_data <- function(
    dataset,
    standard_error,
    fit,
    unit_index,
    grid = seq(0, 15, by = 0.1),
    sample_size = 10000L,
    seed = 20260810L) {
  if (inherits(fit, "linear_mixture_fash")) {
    validate_linear_mixture_fash(fit)
  } else if (inherits(fit, "compact_linear_mixture_fash")) {
    validate_compact_linear_mixture_fash(fit)
  } else {
    stop("fit must be a full or compact linear-mixture FASH object.")
  }
  if (!is.data.frame(dataset) || !all(c("x", "y") %in% names(dataset))) {
    stop("dataset must contain x and y columns.")
  }
  x <- as.numeric(dataset$x)
  y <- as.numeric(dataset$y)
  standard_error <- as.numeric(standard_error)
  grid <- as.numeric(grid)
  if (length(x) < 2L ||
      length(y) != length(x) ||
      length(standard_error) != length(x) ||
      any(!is.finite(x)) ||
      any(!is.finite(y)) ||
      any(!is.finite(standard_error)) ||
      any(standard_error <= 0) ||
      diff(range(x)) <= 0) {
    stop("The observed trajectory or standard errors are invalid.")
  }
  if (length(grid) < 1L || any(!is.finite(grid))) {
    stop("grid must contain finite evaluation times.")
  }
  if (length(unit_index) != 1L ||
      !is.finite(unit_index) ||
      unit_index != as.integer(unit_index) ||
      unit_index < 1L ||
      unit_index > nrow(fit$posterior_weights)) {
    stop("unit_index must identify one row of fit$posterior_weights.")
  }
  if (length(sample_size) != 1L ||
      !is.finite(sample_size) ||
      sample_size < 100L ||
      sample_size != as.integer(sample_size)) {
    stop("sample_size must be one integer of at least 100.")
  }
  if (length(seed) != 1L || !is.finite(seed)) {
    stop("seed must be one finite value.")
  }
  unit_index <- as.integer(unit_index)
  sample_size <- as.integer(sample_size)
  pred_step <- fit$settings$pred_step
  normalized_time <- (x - min(x)) / pred_step
  normalized_grid <- (grid - min(x)) / pred_step
  weight <- 1 / (standard_error^2 + 1e-10)

  component_sd <- as.numeric(colnames(fit$posterior_weights))
  component_probability <- as.numeric(
    fit$posterior_weights[unit_index, , drop = TRUE]
  )
  if (any(!is.finite(component_sd)) ||
      component_sd[1] != 0 ||
      any(component_probability < -1e-12)) {
    stop("The unit-specific mixture components are invalid.")
  }
  component_probability <- pmax(component_probability, 0)
  component_probability <- component_probability /
    sum(component_probability)

  set.seed(seed)
  component_draw <- sample.int(
    length(component_sd),
    size = sample_size,
    replace = TRUE,
    prob = component_probability
  )
  intercept_draw <- numeric(sample_size)
  slope_draw <- numeric(sample_size)
  design <- cbind(1, normalized_time)

  for (component_index in seq_along(component_sd)) {
    selected <- which(component_draw == component_index)
    selected_count <- length(selected)
    if (selected_count == 0L) {
      next
    }
    slope_sd <- component_sd[component_index]
    if (slope_sd == 0) {
      intercept_variance <- 1 / sum(weight)
      intercept_mean <- sum(weight * y) * intercept_variance
      intercept_draw[selected] <- stats::rnorm(
        selected_count,
        mean = intercept_mean,
        sd = sqrt(intercept_variance)
      )
      next
    }

    prior_precision <- matrix(0, nrow = 2L, ncol = 2L)
    prior_precision[2L, 2L] <- 1 / slope_sd^2
    posterior_precision <- crossprod(design, weight * design) +
      prior_precision
    posterior_covariance <- solve(posterior_precision)
    posterior_mean <- as.numeric(
      posterior_covariance %*% crossprod(design, weight * y)
    )
    standard_draw <- matrix(
      stats::rnorm(2L * selected_count),
      ncol = 2L
    )
    coefficient_draw <- sweep(
      standard_draw %*% chol(posterior_covariance),
      2L,
      posterior_mean,
      `+`
    )
    intercept_draw[selected] <- coefficient_draw[, 1L]
    slope_draw[selected] <- coefficient_draw[, 2L]
  }

  trajectory_draws <- outer(intercept_draw, rep(1, length(grid))) +
    outer(slope_draw, normalized_grid)
  intervals <- apply(
    trajectory_draws,
    2L,
    stats::quantile,
    probs = c(0.025, 0.975),
    names = FALSE
  )
  data.frame(
    time = grid,
    posterior_mean = colMeans(trajectory_draws),
    lower = intervals[1L, ],
    upper = intervals[2L, ],
    stringsAsFactors = FALSE
  )
}
