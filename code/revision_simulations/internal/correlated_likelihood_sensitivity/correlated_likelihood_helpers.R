# Helpers for FASH sensitivity analyses with shared time correlation.

validate_shared_correlation <- function(correlation,
                                        n_time = nrow(correlation),
                                        name = "correlation",
                                        tolerance = 1e-10) {
  correlation <- as.matrix(correlation)
  n_time <- as.integer(n_time)
  tolerance <- as.numeric(tolerance)
  if (length(n_time) != 1L || is.na(n_time) || n_time < 2L ||
      length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0 ||
      !identical(dim(correlation), c(n_time, n_time)) ||
      any(!is.finite(correlation)) ||
      max(abs(correlation - t(correlation))) > tolerance ||
      max(abs(diag(correlation) - 1)) > 1e-8) {
    stop(name, " must be a finite symmetric unit-diagonal matrix.")
  }
  correlation <- (correlation + t(correlation)) / 2
  eigenvalues <- eigen(
    correlation,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  if (min(eigenvalues) <= tolerance) {
    stop(name, " must be strictly positive definite.")
  }
  determinant_result <- determinant(correlation, logarithm = TRUE)
  if (determinant_result$sign <= 0) {
    stop(name, " must have a positive determinant.")
  }
  list(
    matrix = correlation,
    minimum_eigenvalue = min(eigenvalues),
    maximum_eigenvalue = max(eigenvalues),
    eigenvalue_condition_number = max(eigenvalues) / min(eigenvalues),
    log_determinant = as.numeric(determinant_result$modulus)
  )
}

construct_unit_precision_matrices <- function(adjusted_se,
                                              correlation,
                                              validate_all = TRUE) {
  adjusted_se <- as.matrix(adjusted_se)
  validate_all <- isTRUE(validate_all)
  if (nrow(adjusted_se) < 1L || ncol(adjusted_se) < 2L ||
      any(!is.finite(adjusted_se)) || any(adjusted_se <= 0)) {
    stop("adjusted_se must be a finite positive matrix.")
  }
  correlation_diagnostics <- validate_shared_correlation(
    correlation,
    n_time = ncol(adjusted_se)
  )
  correlation <- correlation_diagnostics$matrix
  correlation_precision <- solve(correlation)
  correlation_precision <- (correlation_precision +
    t(correlation_precision)) / 2
  time_names <- colnames(adjusted_se)
  pair_names <- rownames(adjusted_se)
  precision_matrices <- lapply(seq_len(nrow(adjusted_se)), function(index) {
    inverse_se <- 1 / adjusted_se[index, ]
    precision <- correlation_precision * tcrossprod(inverse_se)
    precision <- (precision + t(precision)) / 2
    dimnames(precision) <- list(time_names, time_names)
    precision
  })
  names(precision_matrices) <- pair_names
  if (validate_all) {
    valid_precision <- vapply(precision_matrices, function(precision) {
      isTRUE(tryCatch({
        chol(precision)
        TRUE
      }, error = function(error) FALSE))
    }, logical(1))
    if (!all(valid_precision)) {
      stop("At least one constructed precision matrix is not positive definite.")
    }
  }
  list(
    precision_matrices = precision_matrices,
    correlation_precision = correlation_precision,
    correlation_diagnostics = correlation_diagnostics
  )
}

make_fash_data_list <- function(beta_hat, time_grid) {
  beta_hat <- as.matrix(beta_hat)
  time_grid <- as.numeric(time_grid)
  pair_keys <- rownames(beta_hat)
  if (nrow(beta_hat) < 2L || ncol(beta_hat) < 2L ||
      any(!is.finite(beta_hat)) || length(time_grid) != ncol(beta_hat) ||
      any(!is.finite(time_grid)) || anyDuplicated(time_grid) ||
      is.null(pair_keys) || any(!nzchar(pair_keys)) ||
      anyDuplicated(pair_keys)) {
    stop("beta_hat and time_grid are invalid or beta_hat lacks unique row names.")
  }
  data_list <- lapply(seq_len(nrow(beta_hat)), function(index) {
    data.frame(
      beta = as.numeric(beta_hat[index, ]),
      time = time_grid
    )
  })
  names(data_list) <- pair_keys
  data_list
}

validate_fash_settings <- function(settings) {
  required <- c(
    "num_basis", "order", "betaprec", "pred_step", "likelihood", "penalty"
  )
  if (!is.list(settings) || !all(required %in% names(settings)) ||
      !identical(as.character(settings$likelihood), "gaussian") ||
      !is.finite(settings$num_basis) || settings$num_basis < 2 ||
      !is.finite(settings$order) || settings$order < 1 ||
      !is.finite(settings$betaprec) || settings$betaprec < 0 ||
      !is.finite(settings$pred_step) || settings$pred_step <= 0 ||
      !is.finite(settings$penalty) || settings$penalty < 1) {
    stop("The source FASH settings are incomplete or invalid.")
  }
  settings
}

fit_fash_with_shared_correlation <- function(beta_hat,
                                             adjusted_se,
                                             time_grid,
                                             correlation,
                                             settings,
                                             psd_grid,
                                             num_cores = 1L,
                                             verbose = TRUE) {
  if (!requireNamespace("fashr", quietly = TRUE)) {
    stop("The fashr package is required for correlated-likelihood fitting.")
  }
  beta_hat <- as.matrix(beta_hat)
  adjusted_se <- as.matrix(adjusted_se)
  settings <- validate_fash_settings(settings)
  psd_grid <- as.numeric(psd_grid)
  num_cores <- as.integer(num_cores)
  if (!identical(dim(beta_hat), dim(adjusted_se)) ||
      !identical(rownames(beta_hat), rownames(adjusted_se)) ||
      length(psd_grid) < 2L || any(!is.finite(psd_grid)) ||
      any(psd_grid < 0) || anyDuplicated(psd_grid) ||
      !any(psd_grid == 0) || length(num_cores) != 1L ||
      is.na(num_cores) || num_cores < 1L) {
    stop("The correlated FASH inputs are invalid or misaligned.")
  }
  data_list <- make_fash_data_list(beta_hat, time_grid)
  precision <- construct_unit_precision_matrices(
    adjusted_se,
    correlation,
    validate_all = TRUE
  )
  fit <- fashr::fash(
    Y = "beta",
    smooth_var = "time",
    Omega = precision$precision_matrices,
    data_list = data_list,
    grid = psd_grid,
    likelihood = settings$likelihood,
    num_basis = settings$num_basis,
    betaprec = settings$betaprec,
    order = settings$order,
    pred_step = settings$pred_step,
    penalty = settings$penalty,
    num_cores = num_cores,
    verbose = verbose
  )
  names(fit$lfdr) <- rownames(beta_hat)
  rownames(fit$posterior_weights) <- rownames(beta_hat)
  list(
    fit = fit,
    correlation_diagnostics = precision$correlation_diagnostics,
    correlation_precision = precision$correlation_precision
  )
}

validate_probability_vector <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) < 2L || any(!is.finite(x)) || any(x < 0 | x > 1)) {
    stop(name, " must contain finite probabilities in [0, 1].")
  }
  x
}

validate_prior_weight_table <- function(prior_weights, name) {
  prior_weights <- as.data.frame(prior_weights, stringsAsFactors = FALSE)
  if (!all(c("psd", "prior_weight") %in% names(prior_weights))) {
    stop(name, " must contain psd and prior_weight columns.")
  }
  prior_weights <- prior_weights[, c("psd", "prior_weight")]
  prior_weights$psd <- as.numeric(prior_weights$psd)
  prior_weights$prior_weight <- as.numeric(prior_weights$prior_weight)
  if (nrow(prior_weights) < 1L || any(!is.finite(as.matrix(prior_weights))) ||
      any(prior_weights$psd < 0) || any(prior_weights$prior_weight < 0) ||
      anyDuplicated(prior_weights$psd) ||
      abs(sum(prior_weights$prior_weight) - 1) > 1e-6) {
    stop(name, " is not a normalized FASH prior-weight table.")
  }
  prior_weights[order(prior_weights$psd), , drop = FALSE]
}

extract_fit_prior_weights <- function(fit, name) {
  prior <- validate_prior_weight_table(fit$prior_weights, name)
  psd_grid <- as.numeric(fit$psd_grid)
  if (length(psd_grid) < 2L || any(!is.finite(psd_grid)) ||
      any(psd_grid < 0) || anyDuplicated(psd_grid) ||
      sum(psd_grid == 0) != 1L ||
      any(!prior$psd %in% psd_grid)) {
    stop(name, " is not aligned to a valid complete PSD grid.")
  }
  complete <- data.frame(
    psd = psd_grid,
    prior_weight = 0,
    stringsAsFactors = FALSE
  )
  matched <- match(prior$psd, complete$psd)
  complete$prior_weight[matched] <- prior$prior_weight
  complete
}

cumulative_lfdr_calls <- function(lfdr, alpha = 0.05) {
  lfdr <- validate_probability_vector(lfdr, "lfdr")
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must be strictly between zero and one.")
  }
  ordering <- order(lfdr, method = "radix")
  selected <- which(cumsum(lfdr[ordering]) / seq_along(ordering) <= alpha)
  if (length(selected) == 0L) integer(0) else ordering[selected]
}

safe_correlation <- function(x, y, method) {
  if (stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }
  unname(stats::cor(x, y, method = method))
}

compare_prior_weight_fits <- function(fits, fit_stage, method_labels = names(fits)) {
  if (!is.list(fits) || length(fits) < 2L || is.null(names(fits)) ||
      any(!nzchar(names(fits))) || anyDuplicated(names(fits)) ||
      length(method_labels) != length(fits)) {
    stop("fits must be a named list with at least two methods.")
  }
  prior_tables <- Map(function(fit, method_id, method_label) {
    prior <- extract_fit_prior_weights(
      fit,
      paste0(method_id, " prior weights")
    )
    data.frame(
      fit_stage = fit_stage,
      method_id = method_id,
      method_label = method_label,
      psd = prior$psd,
      prior_weight = prior$prior_weight,
      stringsAsFactors = FALSE
    )
  }, fits, names(fits), method_labels)
  reference_grid <- prior_tables[[1]]$psd
  if (!all(vapply(prior_tables, function(table) {
    identical(table$psd, reference_grid)
  }, logical(1)))) {
    stop("All fits must use the same ordered PSD grid.")
  }
  comparisons <- utils::combn(seq_along(fits), 2L)
  metrics <- lapply(seq_len(ncol(comparisons)), function(index) {
    reference_index <- comparisons[1L, index]
    comparison_index <- comparisons[2L, index]
    reference <- prior_tables[[reference_index]]
    comparison <- prior_tables[[comparison_index]]
    difference <- comparison$prior_weight - reference$prior_weight
    null_index <- which(reference$psd == 0)
    data.frame(
      fit_stage = fit_stage,
      reference_method_id = names(fits)[reference_index],
      comparison_method_id = names(fits)[comparison_index],
      reference_method_label = method_labels[reference_index],
      comparison_method_label = method_labels[comparison_index],
      reference_pi0 = reference$prior_weight[null_index],
      comparison_pi0 = comparison$prior_weight[null_index],
      pi0_difference = difference[null_index],
      prior_total_variation = 0.5 * sum(abs(difference)),
      maximum_absolute_weight_difference = max(abs(difference)),
      stringsAsFactors = FALSE
    )
  })
  list(
    prior_weights = do.call(rbind, prior_tables),
    pairwise_metrics = do.call(rbind, metrics)
  )
}

compare_lfdr_fits <- function(fits,
                              pair_metadata,
                              fit_stage,
                              method_labels = names(fits),
                              alpha = 0.05,
                              top_n = 100L) {
  if (!is.list(fits) || length(fits) < 2L || is.null(names(fits)) ||
      any(!nzchar(names(fits))) || anyDuplicated(names(fits)) ||
      length(method_labels) != length(fits) ||
      !is.data.frame(pair_metadata) ||
      !all(c("pair_key", "gene_id", "variant_id") %in% names(pair_metadata)) ||
      anyDuplicated(pair_metadata$pair_key)) {
    stop("The fit list or pair metadata is invalid.")
  }
  pair_keys <- as.character(pair_metadata$pair_key)
  lfdr_by_method <- lapply(seq_along(fits), function(index) {
    lfdr <- validate_probability_vector(
      fits[[index]]$lfdr,
      paste0(names(fits)[index], " lfdr")
    )
    fit_keys <- names(fits[[index]]$lfdr)
    if (is.null(fit_keys)) {
      fit_keys <- rownames(fits[[index]]$posterior_weights)
    }
    if (!identical(as.character(fit_keys), pair_keys)) {
      stop("Every fit must have lfdr values aligned to pair_metadata.")
    }
    lfdr
  })
  names(lfdr_by_method) <- names(fits)
  lfdr_wide <- pair_metadata[, c("pair_key", "gene_id", "variant_id")]
  for (method_id in names(fits)) {
    lfdr_wide[[method_id]] <- lfdr_by_method[[method_id]]
  }
  lfdr_long <- do.call(rbind, lapply(seq_along(fits), function(index) {
    data.frame(
      fit_stage = fit_stage,
      method_id = names(fits)[index],
      method_label = method_labels[index],
      pair_key = pair_keys,
      gene_id = pair_metadata$gene_id,
      variant_id = pair_metadata$variant_id,
      lfdr = lfdr_by_method[[index]],
      stringsAsFactors = FALSE
    )
  }))
  calls <- lapply(lfdr_by_method, cumulative_lfdr_calls, alpha = alpha)
  discovery_summary <- do.call(rbind, lapply(seq_along(fits), function(index) {
    call_indices <- calls[[index]]
    data.frame(
      fit_stage = fit_stage,
      method_id = names(fits)[index],
      method_label = method_labels[index],
      n_units = length(pair_keys),
      mean_lfdr = mean(lfdr_by_method[[index]]),
      median_lfdr = stats::median(lfdr_by_method[[index]]),
      discovered_units = length(call_indices),
      discovered_genes = length(unique(pair_metadata$gene_id[call_indices])),
      alpha = alpha,
      stringsAsFactors = FALSE
    )
  }))
  comparisons <- utils::combn(seq_along(fits), 2L)
  metric_rows <- vector("list", ncol(comparisons))
  top_rows <- vector("list", ncol(comparisons))
  for (index in seq_len(ncol(comparisons))) {
    reference_index <- comparisons[1L, index]
    comparison_index <- comparisons[2L, index]
    reference_id <- names(fits)[reference_index]
    comparison_id <- names(fits)[comparison_index]
    reference_lfdr <- lfdr_by_method[[reference_id]]
    comparison_lfdr <- lfdr_by_method[[comparison_id]]
    difference <- comparison_lfdr - reference_lfdr
    reference_calls <- calls[[reference_id]]
    comparison_calls <- calls[[comparison_id]]
    call_union <- union(reference_calls, comparison_calls)
    call_intersection <- intersect(reference_calls, comparison_calls)
    jaccard <- if (length(call_union) == 0L) {
      1
    } else {
      length(call_intersection) / length(call_union)
    }
    metric_rows[[index]] <- data.frame(
      fit_stage = fit_stage,
      reference_method_id = reference_id,
      comparison_method_id = comparison_id,
      reference_method_label = method_labels[reference_index],
      comparison_method_label = method_labels[comparison_index],
      pearson_lfdr = safe_correlation(reference_lfdr, comparison_lfdr, "pearson"),
      spearman_lfdr = safe_correlation(reference_lfdr, comparison_lfdr, "spearman"),
      mean_signed_lfdr_difference = mean(difference),
      mean_absolute_lfdr_difference = mean(abs(difference)),
      median_absolute_lfdr_difference = stats::median(abs(difference)),
      rmse_lfdr = sqrt(mean(difference^2)),
      maximum_absolute_lfdr_difference = max(abs(difference)),
      reference_discoveries = length(reference_calls),
      comparison_discoveries = length(comparison_calls),
      discovery_intersection = length(call_intersection),
      discovery_union = length(call_union),
      discovery_jaccard = jaccard,
      alpha = alpha,
      stringsAsFactors = FALSE
    )
    ordering <- order(abs(difference), decreasing = TRUE, method = "radix")
    ordering <- ordering[seq_len(min(as.integer(top_n), length(ordering)))]
    top_rows[[index]] <- data.frame(
      fit_stage = fit_stage,
      reference_method_id = reference_id,
      comparison_method_id = comparison_id,
      discrepancy_rank = seq_along(ordering),
      pair_key = pair_keys[ordering],
      gene_id = pair_metadata$gene_id[ordering],
      variant_id = pair_metadata$variant_id[ordering],
      reference_lfdr = reference_lfdr[ordering],
      comparison_lfdr = comparison_lfdr[ordering],
      lfdr_difference = difference[ordering],
      absolute_lfdr_difference = abs(difference[ordering]),
      stringsAsFactors = FALSE
    )
  }
  list(
    lfdr_wide = lfdr_wide,
    lfdr_long = lfdr_long,
    discovery_summary = discovery_summary,
    pairwise_metrics = do.call(rbind, metric_rows),
    top_discrepancies = do.call(rbind, top_rows)
  )
}

validate_identity_path_equivalence <- function(standard_error_fit,
                                               identity_precision_fit,
                                               likelihood_tolerance = 1e-7,
                                               prior_tolerance = 1e-6,
                                               lfdr_tolerance = 1e-6) {
  standard_likelihood <- as.matrix(standard_error_fit$L_matrix)
  identity_likelihood <- as.matrix(identity_precision_fit$L_matrix)
  if (!identical(dim(standard_likelihood), dim(identity_likelihood)) ||
      ncol(standard_likelihood) < 2L) {
    stop("Identity-path fits have incompatible likelihood matrices.")
  }
  standard_centered <- standard_likelihood - standard_likelihood[, 1L]
  identity_centered <- identity_likelihood - identity_likelihood[, 1L]
  likelihood_difference <- max(abs(standard_centered - identity_centered))
  standard_prior <- extract_fit_prior_weights(
    standard_error_fit,
    "standard_error_fit prior"
  )
  identity_prior <- extract_fit_prior_weights(
    identity_precision_fit,
    "identity_precision_fit prior"
  )
  prior_difference <- max(abs(
    standard_prior$prior_weight - identity_prior$prior_weight
  ))
  lfdr_difference <- max(abs(
    as.numeric(standard_error_fit$lfdr) -
      as.numeric(identity_precision_fit$lfdr)
  ))
  result <- list(
    row_centered_likelihood_maximum_difference = likelihood_difference,
    prior_weight_maximum_difference = prior_difference,
    lfdr_maximum_difference = lfdr_difference,
    tolerances = c(
      likelihood = likelihood_tolerance,
      prior_weight = prior_tolerance,
      lfdr = lfdr_tolerance
    )
  )
  if (likelihood_difference > likelihood_tolerance ||
      prior_difference > prior_tolerance ||
      lfdr_difference > lfdr_tolerance) {
    stop("The identity-precision and diagonal-SE likelihood paths disagree.")
  }
  result
}

run_bf_updates_checked <- function(raw_fits,
                                   method_labels = names(raw_fits),
                                   pair_keys = names(raw_fits[[1L]]$lfdr)) {
  if (!requireNamespace("fashr", quietly = TRUE)) {
    stop("The fashr package is required for BF updates.")
  }
  if (!is.list(raw_fits) || length(raw_fits) < 1L ||
      is.null(names(raw_fits)) || any(!nzchar(names(raw_fits))) ||
      anyDuplicated(names(raw_fits)) ||
      length(method_labels) != length(raw_fits) ||
      length(pair_keys) < 2L || any(!nzchar(pair_keys)) ||
      anyDuplicated(pair_keys)) {
    stop("The raw fit list, method labels, or pair keys are invalid.")
  }
  successful_fits <- list()
  warning_messages <- list()
  status_rows <- vector("list", length(raw_fits))
  elapsed_seconds <- numeric(length(raw_fits))
  names(elapsed_seconds) <- names(raw_fits)
  for (index in seq_along(raw_fits)) {
    method_id <- names(raw_fits)[index]
    raw_fit <- raw_fits[[index]]
    fit_pair_keys <- names(raw_fit$lfdr)
    if (is.null(fit_pair_keys)) {
      fit_pair_keys <- rownames(raw_fit$posterior_weights)
    }
    if (length(raw_fit$lfdr) != length(pair_keys) ||
        nrow(raw_fit$L_matrix) != length(pair_keys) ||
        !identical(as.character(fit_pair_keys), pair_keys)) {
      stop("Raw fit alignment failed for ", method_id, ".")
    }
    names(raw_fit$lfdr) <- pair_keys
    rownames(raw_fit$posterior_weights) <- pair_keys
    rownames(raw_fit$L_matrix) <- pair_keys
    complete_prior <- extract_fit_prior_weights(
      raw_fit,
      paste0(method_id, " raw prior")
    )
    alternative_mass <- sum(
      complete_prior$prior_weight[complete_prior$psd > 0]
    )
    centered_input <- raw_fit
    row_offsets <- apply(centered_input$L_matrix, 1L, max)
    centered_input$L_matrix <- centered_input$L_matrix - row_offsets
    warnings <- character()
    update_start <- proc.time()[["elapsed"]]
    updated <- withCallingHandlers(
      fashr::BF_update(centered_input, plot = FALSE),
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    )
    elapsed_seconds[method_id] <- proc.time()[["elapsed"]] - update_start
    warnings <- unique(warnings)
    success <- "BF" %in% names(updated) &&
      length(updated$BF) == length(pair_keys) &&
      any(is.finite(updated$BF)) &&
      length(updated$lfdr) == length(pair_keys) &&
      all(is.finite(updated$lfdr)) &&
      all(updated$lfdr >= 0 & updated$lfdr <= 1)
    if (success) {
      updated$L_matrix <- raw_fit$L_matrix
      names(updated$lfdr) <- pair_keys
      rownames(updated$posterior_weights) <- pair_keys
      rownames(updated$L_matrix) <- pair_keys
      successful_fits[[method_id]] <- updated
    }
    warning_messages[[method_id]] <- warnings
    status_rows[[index]] <- data.frame(
      method_id = method_id,
      method_label = method_labels[index],
      raw_alternative_prior_mass = alternative_mass,
      bf_update_available = success,
      bf_update_status = if (success) {
        "Available"
      } else {
        "Unavailable: conditional alternative mixture is undefined"
      },
      warning_message = if (length(warnings) == 0L) {
        ""
      } else {
        paste(warnings, collapse = " | ")
      },
      elapsed_seconds = elapsed_seconds[method_id],
      row_offset_minimum = min(row_offsets),
      row_offset_maximum = max(row_offsets),
      stringsAsFactors = FALSE
    )
  }
  list(
    successful_fits = successful_fits,
    status = do.call(rbind, status_rows),
    warnings = warning_messages,
    elapsed_seconds = elapsed_seconds,
    strategy = paste(
      "Subtract each row maximum from L_matrix before BF_update;",
      "restore the original L_matrix after a successful update."
    )
  )
}
