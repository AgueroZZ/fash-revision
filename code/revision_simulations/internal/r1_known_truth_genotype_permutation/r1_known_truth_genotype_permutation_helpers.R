# Helpers for the R1 known-truth shared-genotype-permutation pilot.

known_truth_alpha_curve <- function(lfdr,
                                    true_null,
                                    alpha_grid,
                                    arm,
                                    fit_stage) {
  lfdr <- as.numeric(lfdr)
  true_null <- as.logical(true_null)
  alpha_grid <- sort(unique(as.numeric(alpha_grid)))
  arm <- as.character(arm)
  fit_stage <- as.character(fit_stage)
  if (length(lfdr) < 1L || length(true_null) != length(lfdr) ||
      any(!is.finite(lfdr)) || any(lfdr < 0 | lfdr > 1) ||
      anyNA(true_null) || length(alpha_grid) < 1L ||
      any(!is.finite(alpha_grid)) || any(alpha_grid < 0 | alpha_grid >= 1) ||
      length(arm) != 1L || !nzchar(arm) || length(fit_stage) != 1L ||
      !nzchar(fit_stage)) {
    stop("Invalid known-truth alpha-curve inputs.")
  }

  true_alternatives <- sum(!true_null)
  rows <- lapply(alpha_grid, function(alpha) {
    selected_indices <- if (alpha == 0) {
      integer()
    } else {
      cumulative_lfdr_calls(lfdr, alpha = alpha)
    }
    selected <- seq_along(lfdr) %in% selected_indices
    discoveries <- sum(selected)
    false_discoveries <- sum(selected & true_null)
    true_positives <- sum(selected & !true_null)
    data.frame(
      arm = arm,
      fit_stage = fit_stage,
      alpha = alpha,
      n_units = length(lfdr),
      n_true_null = sum(true_null),
      n_true_alternative = true_alternatives,
      n_discoveries = discoveries,
      false_discoveries = false_discoveries,
      true_positives = true_positives,
      realized_fdp = if (discoveries == 0L) 0 else {
        false_discoveries / discoveries
      },
      power = if (true_alternatives == 0L) NA_real_ else {
        true_positives / true_alternatives
      },
      mean_selected_lfdr = if (discoveries == 0L) NA_real_ else {
        mean(lfdr[selected])
      },
      maximum_selected_lfdr = if (discoveries == 0L) NA_real_ else {
        max(lfdr[selected])
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_null_bf <- function(bf, true_null, arm) {
  bf <- as.numeric(bf)
  true_null <- as.logical(true_null)
  arm <- as.character(arm)
  if (length(bf) < 1L || length(true_null) != length(bf) ||
      any(!is.finite(bf)) || any(bf <= 0) || anyNA(true_null) ||
      length(arm) != 1L || !nzchar(arm) || !any(true_null)) {
    stop("Invalid known-null Bayes-factor inputs.")
  }
  null_bf <- bf[true_null]
  descending <- order(null_bf, decreasing = TRUE)
  total_mass <- sum(null_bf)
  mass_share <- function(proportion) {
    n_top <- ceiling(length(null_bf) * proportion)
    sum(null_bf[descending[seq_len(n_top)]]) / total_mass
  }
  quantiles <- stats::quantile(
    null_bf,
    probs = c(0.90, 0.95, 0.99),
    names = FALSE,
    type = 8
  )
  data.frame(
    arm = arm,
    n_null = length(null_bf),
    mean_bf = mean(null_bf),
    median_bf = stats::median(null_bf),
    sd_bf = stats::sd(null_bf),
    minimum_bf = min(null_bf),
    q90_bf = quantiles[1L],
    q95_bf = quantiles[2L],
    q99_bf = quantiles[3L],
    maximum_bf = max(null_bf),
    proportion_bf_greater_than_one = mean(null_bf > 1),
    top_1_percent_bf_mass_share = mass_share(0.01),
    top_5_percent_bf_mass_share = mass_share(0.05),
    stringsAsFactors = FALSE
  )
}

residualized_genotype_alignment <- function(original_genotype,
                                             permuted_genotype,
                                             covariates) {
  original_genotype <- as.matrix(original_genotype)
  permuted_genotype <- as.matrix(permuted_genotype)
  covariates <- as.matrix(covariates)
  storage.mode(original_genotype) <- "double"
  storage.mode(permuted_genotype) <- "double"
  storage.mode(covariates) <- "double"
  if (!identical(dim(original_genotype), dim(permuted_genotype)) ||
      nrow(covariates) != nrow(original_genotype) ||
      nrow(original_genotype) < 4L || ncol(original_genotype) < 1L ||
      any(!is.finite(original_genotype)) ||
      any(!is.finite(permuted_genotype)) || any(!is.finite(covariates))) {
    stop("Invalid genotype-alignment inputs.")
  }
  design <- cbind(intercept = 1, covariates)
  design_qr <- qr(design)
  if (design_qr$rank != ncol(design)) {
    stop("The genotype-alignment covariate design is rank deficient.")
  }
  residualizer <- diag(nrow(design)) - tcrossprod(qr.Q(design_qr))
  original_residual <- residualizer %*% original_genotype
  permuted_residual <- residualizer %*% permuted_genotype
  original_ss <- colSums(original_residual^2)
  permuted_ss <- colSums(permuted_residual^2)
  cross_product <- colSums(original_residual * permuted_residual)
  tolerance <- 1e-12 * pmax(1, colSums(original_genotype^2))
  if (any(original_ss <= tolerance) || any(permuted_ss <= tolerance)) {
    stop("At least one residualized genotype vector is degenerate.")
  }
  unit_key <- colnames(original_genotype)
  if (is.null(unit_key)) {
    unit_key <- paste0("unit_", seq_len(ncol(original_genotype)))
  }
  data.frame(
    unit_key = unit_key,
    residualized_correlation = cross_product / sqrt(original_ss * permuted_ss),
    absolute_residualized_correlation = abs(
      cross_product / sqrt(original_ss * permuted_ss)
    ),
    leakage_coefficient = cross_product / permuted_ss,
    original_residual_sum_squares = original_ss,
    permuted_residual_sum_squares = permuted_ss,
    stringsAsFactors = FALSE
  )
}

make_signal_stripped_residual_null <- function(
    genotype,
    expression,
    covariates,
    donor_map,
    leverage_adjustment = c("HC2", "none"),
    tolerance = 1e-10) {
  leverage_adjustment <- match.arg(leverage_adjustment)
  genotype <- as.matrix(genotype)
  covariates <- as.matrix(covariates)
  storage.mode(genotype) <- "double"
  storage.mode(covariates) <- "double"
  expression <- as.array(expression)
  storage.mode(expression) <- "double"
  tolerance <- as.numeric(tolerance)
  donor_ids <- rownames(genotype)
  unit_ids <- colnames(genotype)
  expression_dimnames <- dimnames(expression)
  required_map_columns <- c("target_donor", "source_donor")
  if (length(dim(expression)) != 3L ||
      !identical(dim(expression)[1:2], dim(genotype)) ||
      nrow(covariates) != nrow(genotype) || nrow(genotype) < 4L ||
      ncol(genotype) < 1L || ncol(covariates) < 1L ||
      nrow(genotype) <= ncol(covariates) + 2L ||
      any(!is.finite(genotype)) || any(!is.finite(expression)) ||
      any(!is.finite(covariates)) ||
      is.null(donor_ids) || is.null(unit_ids) ||
      any(!nzchar(donor_ids)) || any(!nzchar(unit_ids)) ||
      anyDuplicated(donor_ids) || anyDuplicated(unit_ids) ||
      is.null(rownames(covariates)) ||
      !identical(rownames(covariates), donor_ids) ||
      is.null(expression_dimnames) ||
      !identical(expression_dimnames[[1L]], donor_ids) ||
      !identical(expression_dimnames[[2L]], unit_ids) ||
      !is.data.frame(donor_map) ||
      !all(required_map_columns %in% names(donor_map)) ||
      nrow(donor_map) != nrow(genotype) ||
      length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) {
    stop("Invalid signal-stripped residual-null inputs.")
  }
  donor_map$target_donor <- as.character(donor_map$target_donor)
  donor_map$source_donor <- as.character(donor_map$source_donor)
  target_order <- match(donor_ids, donor_map$target_donor)
  if (anyNA(target_order) || anyDuplicated(donor_map$target_donor) ||
      anyDuplicated(donor_map$source_donor) ||
      !setequal(donor_map$target_donor, donor_ids) ||
      !setequal(donor_map$source_donor, donor_ids)) {
    stop("donor_map must contain a one-to-one permutation of donor IDs.")
  }
  donor_map <- donor_map[target_order, , drop = FALSE]
  rownames(donor_map) <- NULL
  source_index <- match(donor_map$source_donor, donor_ids)
  if (anyNA(source_index)) {
    stop("Could not align donor_map source donors to the genotype matrix.")
  }

  n_donors <- nrow(genotype)
  n_units <- ncol(genotype)
  n_times <- dim(expression)[3L]
  time_ids <- expression_dimnames[[3L]]
  if (is.null(time_ids)) {
    time_ids <- paste0("time_", seq_len(n_times))
  }
  array_dim <- c(n_donors, n_units, n_times)
  array_dimnames <- list(donor_ids, unit_ids, time_ids)
  full_model_residual <- array(
    NA_real_, dim = array_dim, dimnames = array_dimnames
  )
  adjusted_residual <- array(
    NA_real_, dim = array_dim, dimnames = array_dimnames
  )
  nuisance_fitted <- array(
    NA_real_, dim = array_dim, dimnames = array_dimnames
  )
  source_genotype_fitted <- array(
    NA_real_, dim = array_dim, dimnames = array_dimnames
  )
  null_expression <- array(
    NA_real_, dim = array_dim, dimnames = array_dimnames
  )
  leverage <- matrix(
    NA_real_,
    nrow = n_donors,
    ncol = n_units,
    dimnames = list(donor_ids, unit_ids)
  )
  genotype_coefficient <- matrix(
    NA_real_,
    nrow = n_units,
    ncol = n_times,
    dimnames = list(unit_ids, time_ids)
  )
  maximum_design_cross_product <- 0
  maximum_nuisance_partial_genotype <- 0

  for (unit_index in seq_len(n_units)) {
    design <- cbind(
      intercept = 1,
      genotype = genotype[, unit_index],
      covariates
    )
    design_qr <- qr(design)
    if (design_qr$rank != ncol(design)) {
      stop("At least one full-model design is rank deficient.")
    }
    response <- expression[, unit_index, , drop = FALSE][, 1L, ]
    if (n_times == 1L) {
      response <- matrix(response, ncol = 1L)
    }
    coefficients <- qr.coef(design_qr, response)
    fitted_full <- design %*% coefficients
    residual <- response - fitted_full
    unit_leverage <- rowSums(qr.Q(design_qr)^2)
    if (any(!is.finite(unit_leverage)) ||
        any(unit_leverage < -tolerance) || any(unit_leverage > 1 + tolerance) ||
        (leverage_adjustment == "HC2" &&
          any(unit_leverage >= 1 - tolerance))) {
      stop("At least one full-model leverage is invalid for the requested adjustment.")
    }
    nuisance_columns <- c("intercept", colnames(covariates))
    unit_nuisance_fitted <- design[, nuisance_columns, drop = FALSE] %*%
      coefficients[nuisance_columns, , drop = FALSE]
    unit_genotype_fitted <- outer(
      genotype[, unit_index], coefficients["genotype", ]
    )
    unit_adjusted_residual <- if (leverage_adjustment == "HC2") {
      sweep(residual, 1L, sqrt(1 - unit_leverage), `/`)
    } else {
      residual
    }
    permuted_adjusted_residual <- unit_adjusted_residual[
      source_index, , drop = FALSE
    ]
    unit_null_expression <-
      unit_nuisance_fitted + permuted_adjusted_residual

    residual_cross_product <- crossprod(design, residual)
    nuisance_partial_genotype <- qr.coef(
      design_qr, unit_nuisance_fitted
    )["genotype", ]
    maximum_design_cross_product <- max(
      maximum_design_cross_product,
      abs(residual_cross_product)
    )
    maximum_nuisance_partial_genotype <- max(
      maximum_nuisance_partial_genotype,
      abs(nuisance_partial_genotype)
    )

    full_model_residual[, unit_index, ] <- residual
    adjusted_residual[, unit_index, ] <- unit_adjusted_residual
    nuisance_fitted[, unit_index, ] <- unit_nuisance_fitted
    source_genotype_fitted[, unit_index, ] <- unit_genotype_fitted
    null_expression[, unit_index, ] <- unit_null_expression
    leverage[, unit_index] <- unit_leverage
    genotype_coefficient[unit_index, ] <- coefficients["genotype", ]
  }

  response_scale <- max(1, max(abs(expression)))
  if (maximum_design_cross_product > tolerance * response_scale * n_donors ||
      maximum_nuisance_partial_genotype > tolerance * response_scale) {
    stop("Signal-stripped residual-null OLS invariants failed.")
  }
  list(
    null_expression = null_expression,
    nuisance_fitted = nuisance_fitted,
    source_genotype_fitted = source_genotype_fitted,
    full_model_residual = full_model_residual,
    adjusted_residual = adjusted_residual,
    leverage = leverage,
    genotype_coefficient = genotype_coefficient,
    donor_map = donor_map,
    diagnostics = data.frame(
      n_donors = n_donors,
      n_units = n_units,
      n_time_points = n_times,
      leverage_adjustment = leverage_adjustment,
      minimum_leverage = min(leverage),
      maximum_leverage = max(leverage),
      maximum_full_residual_design_cross_product =
        maximum_design_cross_product,
      maximum_nuisance_partial_genotype_coefficient =
        maximum_nuisance_partial_genotype,
      source_genotype_fitted_rms = sqrt(mean(source_genotype_fitted^2)),
      full_model_residual_rms = sqrt(mean(full_model_residual^2)),
      adjusted_residual_rms = sqrt(mean(adjusted_residual^2)),
      null_expression_rms = sqrt(mean(null_expression^2)),
      stringsAsFactors = FALSE
    )
  )
}
