# Helpers for the internal R1 finite-sample normality diagnostic.

validate_finite_numeric <- function(x, name, positive = FALSE) {
  if (!is.numeric(x) || any(!is.finite(x))) {
    stop(name, " must contain only finite numeric values.")
  }
  if (positive && any(x <= 0)) {
    stop(name, " must contain only positive values.")
  }
  invisible(x)
}

gaussian_to_standardized_t <- function(error, df) {
  validate_finite_numeric(error, "error")
  df <- as.numeric(df)
  if (length(df) != 1L || !is.finite(df) || df <= 2) {
    stop("df must be one finite value greater than two.")
  }
  probability <- pmin(
    pmax(stats::pnorm(error), .Machine$double.eps),
    1 - .Machine$double.eps
  )
  stats::qt(probability, df = df) * sqrt((df - 2) / df)
}

reconstruct_r1_expression_mean <- function(G,
                                           beta_matrix,
                                           covariates,
                                           covariate_effects,
                                           intercepts) {
  G <- as.matrix(G)
  beta_matrix <- as.matrix(beta_matrix)
  covariates <- as.matrix(covariates)
  covariate_effects <- as.array(covariate_effects)
  intercepts <- as.matrix(intercepts)
  storage.mode(G) <- "double"
  storage.mode(beta_matrix) <- "double"
  storage.mode(covariates) <- "double"
  storage.mode(covariate_effects) <- "double"
  storage.mode(intercepts) <- "double"

  n_donors <- nrow(G)
  n_units <- ncol(G)
  n_time <- ncol(beta_matrix)
  expected_covariate_dimensions <- c(ncol(covariates), n_units, n_time)
  if (nrow(beta_matrix) != n_units ||
      nrow(covariates) != n_donors ||
      !identical(dim(covariate_effects), expected_covariate_dimensions) ||
      !identical(dim(intercepts), c(n_units, n_time)) ||
      any(!is.finite(G)) || any(!is.finite(beta_matrix)) ||
      any(!is.finite(covariates)) || any(!is.finite(covariate_effects)) ||
      any(!is.finite(intercepts))) {
    stop("Invalid R1 inputs for expression-mean reconstruction.")
  }

  expression_mean <- array(
    NA_real_,
    dim = c(n_donors, n_units, n_time),
    dimnames = list(rownames(G), colnames(G), colnames(beta_matrix))
  )
  for (time_index in seq_len(n_time)) {
    genetic_mean <- sweep(
      G,
      2L,
      beta_matrix[, time_index],
      `*`
    )
    covariate_mean <- covariates %*%
      covariate_effects[, , time_index, drop = FALSE][, , 1L]
    intercept_mean <- matrix(
      intercepts[, time_index],
      nrow = n_donors,
      ncol = n_units,
      byrow = TRUE
    )
    expression_mean[, , time_index] <-
      genetic_mean + covariate_mean + intercept_mean
  }
  expression_mean
}

compute_oracle_regression_se <- function(G,
                                         covariates,
                                         noise_sd,
                                         n_time) {
  G <- as.matrix(G)
  covariates <- as.matrix(covariates)
  storage.mode(G) <- "double"
  storage.mode(covariates) <- "double"
  noise_sd <- as.numeric(noise_sd)
  n_time <- as.integer(n_time)
  if (nrow(covariates) != nrow(G) ||
      any(!is.finite(G)) || any(!is.finite(covariates)) ||
      length(noise_sd) != 1L || !is.finite(noise_sd) || noise_sd <= 0 ||
      length(n_time) != 1L || is.na(n_time) || n_time < 1L) {
    stop("Invalid inputs for oracle regression standard errors.")
  }

  oracle_by_unit <- vapply(seq_len(ncol(G)), function(unit_index) {
    design <- cbind(intercept = 1, G = G[, unit_index], covariates)
    if (qr(design)$rank != ncol(design)) {
      stop("The oracle-SE design is rank deficient for unit ", unit_index, ".")
    }
    noise_sd * sqrt(solve(crossprod(design))["G", "G"])
  }, numeric(1))
  oracle_se <- matrix(
    oracle_by_unit,
    nrow = ncol(G),
    ncol = n_time,
    dimnames = list(colnames(G), paste0("time_", seq_len(n_time)))
  )
  oracle_se
}

validate_standardized_error_inputs <- function(beta_hat,
                                                true_beta,
                                                se_uncorrected,
                                                se_adjusted,
                                                oracle_se,
                                                residual_df,
                                                unit_info) {
  matrices <- lapply(list(
    beta_hat = beta_hat,
    true_beta = true_beta,
    se_uncorrected = se_uncorrected,
    se_adjusted = se_adjusted,
    oracle_se = oracle_se,
    residual_df = residual_df
  ), as.matrix)
  dimensions <- lapply(matrices, dim)
  if (!all(vapply(dimensions, identical, logical(1), dimensions[[1L]]))) {
    stop("All standardized-error matrices must have identical dimensions.")
  }
  if (any(!is.finite(matrices$beta_hat)) ||
      any(!is.finite(matrices$true_beta)) ||
      any(!is.finite(matrices$se_uncorrected)) ||
      any(!is.finite(matrices$se_adjusted)) ||
      any(!is.finite(matrices$oracle_se)) ||
      any(!is.finite(matrices$residual_df)) ||
      any(matrices$se_uncorrected <= 0) ||
      any(matrices$se_adjusted <= 0) ||
      any(matrices$oracle_se <= 0) ||
      any(matrices$residual_df <= 0)) {
    stop("Standardized-error matrices contain invalid values.")
  }
  required_unit_columns <- c(
    "unit_index", "unit_id", "effect_class", "scenario"
  )
  if (!is.data.frame(unit_info) ||
      !all(required_unit_columns %in% names(unit_info)) ||
      nrow(unit_info) != nrow(matrices$beta_hat)) {
    stop("unit_info does not align with the standardized-error matrices.")
  }
  matrices
}

make_standardized_error_table <- function(beta_hat,
                                          true_beta,
                                          se_uncorrected,
                                          se_adjusted,
                                          oracle_se,
                                          residual_df,
                                          unit_info,
                                          error_distribution) {
  matrices <- validate_standardized_error_inputs(
    beta_hat = beta_hat,
    true_beta = true_beta,
    se_uncorrected = se_uncorrected,
    se_adjusted = se_adjusted,
    oracle_se = oracle_se,
    residual_df = residual_df,
    unit_info = unit_info
  )
  error_distribution <- as.character(error_distribution)
  if (length(error_distribution) != 1L || !nzchar(error_distribution)) {
    stop("error_distribution must be one non-empty label.")
  }

  n_units <- nrow(matrices$beta_hat)
  n_time <- ncol(matrices$beta_hat)
  unit_index <- rep(seq_len(n_units), times = n_time)
  time_index <- rep(seq_len(n_time), each = n_units)
  beta_vector <- as.vector(matrices$beta_hat)
  truth_vector <- as.vector(matrices$true_beta)
  oracle_vector <- as.vector(matrices$oracle_se)
  df_vector <- as.vector(matrices$residual_df)
  se_matrices <- list(
    oracle_known_sigma = matrices$oracle_se,
    raw_regression = matrices$se_uncorrected,
    t_adjusted = matrices$se_adjusted
  )
  zero_unit <- unit_info$effect_class == "zero"
  zero_observation <- zero_unit[unit_index]

  rows <- list()
  row_index <- 0L
  for (se_scale in names(se_matrices)) {
    se_vector <- as.vector(se_matrices[[se_scale]])
    reference_distribution <- if (se_scale == "raw_regression") {
      "t"
    } else {
      "normal"
    }
    centered_statistic <- (beta_vector - truth_vector) / se_vector
    common <- data.frame(
      error_distribution = error_distribution,
      se_scale = se_scale,
      reference_distribution = reference_distribution,
      unit_index = unit_info$unit_index[unit_index],
      unit_id = unit_info$unit_id[unit_index],
      effect_class = as.character(unit_info$effect_class[unit_index]),
      time_index = time_index,
      beta_hat = beta_vector,
      true_beta = truth_vector,
      se = se_vector,
      oracle_se = oracle_vector,
      se_ratio = se_vector / oracle_vector,
      residual_df = df_vector,
      stringsAsFactors = FALSE
    )

    row_index <- row_index + 1L
    zero_rows <- common[zero_observation, , drop = FALSE]
    zero_rows$population <- "zero_null"
    zero_rows$statistic <- beta_vector[zero_observation] /
      se_vector[zero_observation]
    rows[[row_index]] <- zero_rows

    row_index <- row_index + 1L
    centered_rows <- common
    centered_rows$population <- "all_centered"
    centered_rows$statistic <- centered_statistic
    rows[[row_index]] <- centered_rows
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[, c(
    "error_distribution", "population", "effect_class", "se_scale",
    "reference_distribution", "unit_index", "unit_id", "time_index",
    "beta_hat", "true_beta", "se", "oracle_se", "se_ratio",
    "residual_df", "statistic"
  )]
}

distribution_moments <- function(values) {
  validate_finite_numeric(values, "values")
  if (length(values) < 2L) {
    stop("At least two values are required for distribution summaries.")
  }
  center <- mean(values)
  scale <- stats::sd(values)
  if (!is.finite(scale) || scale <= 0) {
    stop("The distribution summary has zero or invalid spread.")
  }
  standardized <- (values - center) / scale
  c(
    n = length(values),
    mean = center,
    sd = scale,
    median = stats::median(values),
    iqr = stats::IQR(values),
    skewness = mean(standardized^3),
    excess_kurtosis = mean(standardized^4) - 3
  )
}

reference_tail_rate <- function(values,
                                reference_distribution,
                                residual_df,
                                alpha) {
  if (reference_distribution == "normal") {
    threshold <- stats::qnorm(1 - alpha / 2)
  } else if (reference_distribution == "t") {
    unique_df <- unique(as.numeric(residual_df))
    if (length(unique_df) != 1L) {
      stop("A raw-regression summary must have one residual df value.")
    }
    threshold <- stats::qt(1 - alpha / 2, df = unique_df)
  } else {
    stop("Unknown reference distribution: ", reference_distribution)
  }
  mean(abs(values) > threshold)
}

split_error_groups <- function(error_table) {
  group_columns <- c(
    "error_distribution", "population", "effect_class", "se_scale",
    "reference_distribution"
  )
  missing_columns <- setdiff(
    c(group_columns, "residual_df", "statistic"),
    names(error_table)
  )
  if (length(missing_columns) > 0L) {
    stop(
      "The error table is missing columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  interaction_key <- interaction(
    error_table[, group_columns],
    drop = TRUE,
    lex.order = TRUE,
    sep = "\r"
  )
  split(error_table, interaction_key, drop = TRUE)
}

summarize_standardized_errors <- function(error_table) {
  groups <- split_error_groups(error_table)
  rows <- lapply(groups, function(group) {
    moments <- distribution_moments(group$statistic)
    data.frame(
      error_distribution = group$error_distribution[1L],
      population = group$population[1L],
      effect_class = group$effect_class[1L],
      se_scale = group$se_scale[1L],
      reference_distribution = group$reference_distribution[1L],
      residual_df = if (group$reference_distribution[1L] == "t") {
        unique(group$residual_df)
      } else {
        NA_real_
      },
      n = unname(moments["n"]),
      mean = unname(moments["mean"]),
      sd = unname(moments["sd"]),
      median = unname(moments["median"]),
      iqr = unname(moments["iqr"]),
      skewness = unname(moments["skewness"]),
      excess_kurtosis = unname(moments["excess_kurtosis"]),
      tail_rate_0_05 = reference_tail_rate(
        group$statistic,
        group$reference_distribution[1L],
        group$residual_df,
        alpha = 0.05
      ),
      tail_rate_0_01 = reference_tail_rate(
        group$statistic,
        group$reference_distribution[1L],
        group$residual_df,
        alpha = 0.01
      ),
      tail_rate_0_001 = reference_tail_rate(
        group$statistic,
        group$reference_distribution[1L],
        group$residual_df,
        alpha = 0.001
      ),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(
    out$error_distribution,
    out$population,
    out$effect_class,
    out$se_scale
  ), , drop = FALSE]
}

make_qq_quantiles <- function(error_table,
                              probabilities = stats::ppoints(199L)) {
  probabilities <- as.numeric(probabilities)
  if (length(probabilities) < 2L ||
      any(!is.finite(probabilities)) ||
      any(probabilities <= 0 | probabilities >= 1) ||
      anyDuplicated(probabilities)) {
    stop("probabilities must be unique finite values strictly inside (0, 1).")
  }
  groups <- split_error_groups(error_table)
  rows <- lapply(groups, function(group) {
    reference_distribution <- group$reference_distribution[1L]
    reference_df <- if (reference_distribution == "t") {
      unique_df <- unique(group$residual_df)
      if (length(unique_df) != 1L) {
        stop("A raw-regression QQ group must have one residual df value.")
      }
      unique_df
    } else {
      NA_real_
    }
    reference_quantile <- if (reference_distribution == "normal") {
      stats::qnorm(probabilities)
    } else {
      stats::qt(probabilities, df = reference_df)
    }
    data.frame(
      error_distribution = group$error_distribution[1L],
      population = group$population[1L],
      effect_class = group$effect_class[1L],
      se_scale = group$se_scale[1L],
      reference_distribution = reference_distribution,
      reference_df = reference_df,
      probability = probabilities,
      reference_quantile = reference_quantile,
      empirical_quantile = as.numeric(stats::quantile(
        group$statistic,
        probs = probabilities,
        names = FALSE,
        type = 8
      )),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

summarize_se_uncertainty_bins <- function(error_table, n_bins = 4L) {
  n_bins <- as.integer(n_bins)
  if (length(n_bins) != 1L || is.na(n_bins) || n_bins < 2L) {
    stop("n_bins must be one integer of at least two.")
  }
  selected <- error_table[
    error_table$population %in% c("zero_null", "all_centered") &
      error_table$se_scale == "t_adjusted",
    ,
    drop = FALSE
  ]
  if (nrow(selected) < n_bins * 4L) {
    stop("The adjusted centered-error table is too small for binning.")
  }
  group_key <- interaction(
    selected$error_distribution,
    selected$population,
    drop = TRUE,
    lex.order = TRUE,
    sep = "\r"
  )
  groups <- split(selected, group_key, drop = TRUE)
  rows <- lapply(groups, function(group) {
    ordering <- order(group$se_ratio, group$unit_index, group$time_index)
    bin <- integer(nrow(group))
    bin[ordering] <- ceiling(seq_along(ordering) * n_bins / length(ordering))
    bin <- pmin(bin, n_bins)
    group$se_ratio_bin <- bin
    do.call(rbind, lapply(split(group, group$se_ratio_bin), function(bin_group) {
      moments <- distribution_moments(bin_group$statistic)
      data.frame(
        error_distribution = bin_group$error_distribution[1L],
        population = bin_group$population[1L],
        se_ratio_bin = bin_group$se_ratio_bin[1L],
        n = unname(moments["n"]),
        minimum_se_ratio = min(bin_group$se_ratio),
        median_se_ratio = stats::median(bin_group$se_ratio),
        maximum_se_ratio = max(bin_group$se_ratio),
        mean = unname(moments["mean"]),
        sd = unname(moments["sd"]),
        skewness = unname(moments["skewness"]),
        excess_kurtosis = unname(moments["excess_kurtosis"]),
        tail_rate_0_05 = mean(
          abs(bin_group$statistic) > stats::qnorm(0.975)
        ),
        tail_rate_0_01 = mean(
          abs(bin_group$statistic) > stats::qnorm(0.995)
        ),
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(
    out$error_distribution,
    out$population,
    out$se_ratio_bin
  ), , drop = FALSE]
}
