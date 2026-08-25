# Helpers for the selected-signal matched-null permutation pilots.

cumulative_lfdr_calls <- function(lfdr, alpha = 0.05) {
  lfdr <- as.numeric(lfdr)
  alpha <- as.numeric(alpha)
  if (length(lfdr) < 1L || any(!is.finite(lfdr)) ||
      any(lfdr < 0 | lfdr > 1) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("Invalid lfdr vector or alpha.")
  }
  ordering <- order(lfdr, method = "radix")
  accepted_ranks <- which(cumsum(lfdr[ordering]) / seq_along(ordering) <= alpha)
  if (length(accepted_ranks) == 0L) {
    return(integer())
  }
  ordering[accepted_ranks]
}

select_discovered_gene_top_pairs <- function(pair_keys,
                                              bf_lfdr,
                                              alpha = 0.05) {
  pair_keys <- as.character(pair_keys)
  bf_lfdr <- as.numeric(bf_lfdr)
  if (length(pair_keys) != length(bf_lfdr) || length(pair_keys) < 1L ||
      any(!nzchar(pair_keys)) || anyDuplicated(pair_keys) ||
      any(!is.finite(bf_lfdr)) || any(bf_lfdr < 0 | bf_lfdr > 1)) {
    stop("pair_keys and bf_lfdr must be aligned finite unique vectors.")
  }
  gene_id <- sub("_.*$", "", pair_keys)
  variant_id <- sub("^[^_]+_", "", pair_keys)
  discovered_indices <- cumulative_lfdr_calls(bf_lfdr, alpha = alpha)
  if (length(discovered_indices) == 0L) {
    stop("The BF-adjusted fit has no discoveries at the requested alpha.")
  }
  discovered_gene <- sort(unique(gene_id[discovered_indices]))
  gene_indices <- split(
    discovered_indices,
    factor(gene_id[discovered_indices], levels = discovered_gene)
  )
  selected_indices <- vapply(gene_indices, function(indices) {
    minimum_lfdr <- min(bf_lfdr[indices])
    candidates <- indices[bf_lfdr[indices] == minimum_lfdr]
    candidates[order(pair_keys[candidates], method = "radix")][1L]
  }, integer(1))
  selection <- data.frame(
    target_index = seq_along(selected_indices),
    source_fash_index = selected_indices,
    pair_key = pair_keys[selected_indices],
    gene_id = gene_id[selected_indices],
    variant_id = variant_id[selected_indices],
    source_bf_lfdr = bf_lfdr[selected_indices],
    source_pair_level_discovery = TRUE,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  if (nrow(selection) != length(discovered_gene) ||
      anyDuplicated(selection$pair_key) || anyDuplicated(selection$gene_id) ||
      !setequal(selection$gene_id, discovered_gene) ||
      any(!selection$source_fash_index %in% discovered_indices)) {
    stop("The discovered-gene top-pair selection failed validation.")
  }
  attr(selection, "n_pair_level_discoveries") <- length(discovered_indices)
  attr(selection, "alpha") <- alpha
  selection
}

select_random_tested_pair_per_discovered_gene <- function(pair_keys,
                                                           bf_lfdr,
                                                           alpha = 0.05,
                                                           selection_seed) {
  pair_keys <- as.character(pair_keys)
  bf_lfdr <- as.numeric(bf_lfdr)
  selection_seed <- as.integer(selection_seed)
  if (length(pair_keys) != length(bf_lfdr) || length(pair_keys) < 1L ||
      any(!nzchar(pair_keys)) || anyDuplicated(pair_keys) ||
      any(!is.finite(bf_lfdr)) || any(bf_lfdr < 0 | bf_lfdr > 1) ||
      length(selection_seed) != 1L || is.na(selection_seed)) {
    stop("Invalid aligned pairs, BF-adjusted lfdr values, or selection seed.")
  }
  discovered_indices <- cumulative_lfdr_calls(bf_lfdr, alpha = alpha)
  gene_id <- sub("_.*$", "", pair_keys)
  variant_id <- sub("^[^_]+_", "", pair_keys)
  discovered_genes <- unique(gene_id[discovered_indices])
  eligible <- gene_id %in% discovered_genes
  candidate_indices <- split(
    which(eligible),
    factor(gene_id[eligible], levels = discovered_genes),
    drop = FALSE
  )
  if (length(candidate_indices) != length(discovered_genes) ||
      any(lengths(candidate_indices) == 0L)) {
    stop("Could not construct complete tested-variant pools for discovered genes.")
  }

  had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_random_seed) {
    previous_random_seed <- get(".Random.seed", envir = .GlobalEnv)
  }
  on.exit({
    if (had_random_seed) {
      assign(".Random.seed", previous_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(selection_seed)
  selected_indices <- vapply(candidate_indices, function(indices) {
    indices[sample.int(length(indices), size = 1L)]
  }, integer(1))

  selection <- data.frame(
    target_index = seq_along(selected_indices),
    source_fash_index = unname(selected_indices),
    pair_key = pair_keys[selected_indices],
    gene_id = gene_id[selected_indices],
    variant_id = variant_id[selected_indices],
    candidate_variant_count = unname(lengths(candidate_indices)),
    source_bf_lfdr = bf_lfdr[selected_indices],
    source_pair_level_discovery = selected_indices %in% discovered_indices,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  if (nrow(selection) != length(discovered_genes) ||
      anyDuplicated(selection$pair_key) || anyDuplicated(selection$gene_id) ||
      !identical(selection$gene_id, discovered_genes) ||
      any(selection$candidate_variant_count < 1L) ||
      any(gene_id[selection$source_fash_index] != selection$gene_id)) {
    stop("The random tested-pair selection failed validation.")
  }
  attr(selection, "n_pair_level_discoveries") <- length(discovered_indices)
  attr(selection, "alpha") <- alpha
  attr(selection, "selection_seed") <- selection_seed
  attr(selection, "candidate_pool") <- "all_tested_variants_within_discovered_gene"
  selection
}

select_random_tested_pair_per_gene <- function(pair_keys,
                                                bf_lfdr,
                                                alpha = 0.05,
                                                selection_seed) {
  pair_keys <- as.character(pair_keys)
  bf_lfdr <- as.numeric(bf_lfdr)
  selection_seed <- as.integer(selection_seed)
  if (length(pair_keys) != length(bf_lfdr) || length(pair_keys) < 1L ||
      any(!nzchar(pair_keys)) || anyDuplicated(pair_keys) ||
      any(!is.finite(bf_lfdr)) || any(bf_lfdr < 0 | bf_lfdr > 1) ||
      length(selection_seed) != 1L || is.na(selection_seed)) {
    stop("Invalid aligned pairs, BF-adjusted lfdr values, or selection seed.")
  }
  discovered_indices <- cumulative_lfdr_calls(bf_lfdr, alpha = alpha)
  gene_id <- sub("_.*$", "", pair_keys)
  variant_id <- sub("^[^_]+_", "", pair_keys)
  tested_genes <- unique(gene_id)
  candidate_indices <- split(
    seq_along(pair_keys),
    factor(gene_id, levels = tested_genes),
    drop = FALSE
  )
  if (length(candidate_indices) != length(tested_genes) ||
      any(lengths(candidate_indices) == 0L)) {
    stop("Could not construct complete tested-variant pools for all genes.")
  }

  had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_random_seed) {
    previous_random_seed <- get(".Random.seed", envir = .GlobalEnv)
  }
  on.exit({
    if (had_random_seed) {
      assign(".Random.seed", previous_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(selection_seed)
  selected_indices <- vapply(candidate_indices, function(indices) {
    indices[sample.int(length(indices), size = 1L)]
  }, integer(1))

  selection <- data.frame(
    target_index = seq_along(selected_indices),
    source_fash_index = unname(selected_indices),
    pair_key = pair_keys[selected_indices],
    gene_id = gene_id[selected_indices],
    variant_id = variant_id[selected_indices],
    candidate_variant_count = unname(lengths(candidate_indices)),
    source_bf_lfdr = bf_lfdr[selected_indices],
    source_pair_level_discovery = selected_indices %in% discovered_indices,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  if (nrow(selection) != length(tested_genes) ||
      anyDuplicated(selection$pair_key) || anyDuplicated(selection$gene_id) ||
      !identical(selection$gene_id, tested_genes) ||
      any(selection$candidate_variant_count < 1L) ||
      any(gene_id[selection$source_fash_index] != selection$gene_id)) {
    stop("The all-gene random tested-pair selection failed validation.")
  }
  attr(selection, "n_pair_level_discoveries") <- length(discovered_indices)
  attr(selection, "alpha") <- alpha
  attr(selection, "selection_seed") <- selection_seed
  attr(selection, "candidate_pool") <- "all_tested_variants_within_all_genes"
  selection
}

make_shared_genotype_permutation <- function(genotype, seed) {
  genotype <- as.matrix(genotype)
  storage.mode(genotype) <- "double"
  seed <- as.integer(seed)
  donor_ids <- rownames(genotype)
  if (nrow(genotype) < 2L || ncol(genotype) < 1L ||
      is.null(donor_ids) || any(!nzchar(donor_ids)) ||
      anyDuplicated(donor_ids) || any(!is.finite(genotype)) ||
      length(seed) != 1L || is.na(seed)) {
    stop("Invalid genotype matrix or permutation seed.")
  }
  set.seed(seed)
  source_donor <- sample(donor_ids, length(donor_ids), replace = FALSE)
  names(source_donor) <- donor_ids
  permuted <- genotype[source_donor, , drop = FALSE]
  rownames(permuted) <- donor_ids
  if (anyDuplicated(source_donor) || !setequal(source_donor, donor_ids) ||
      !isTRUE(all.equal(
        unname(crossprod(genotype)),
        unname(crossprod(permuted)),
        tolerance = 1e-12
      ))) {
    stop("The shared genotype permutation failed its invariants.")
  }
  list(
    genotype = permuted,
    donor_map = data.frame(
      target_donor = donor_ids,
      source_donor = unname(source_donor),
      fixed_point = donor_ids == unname(source_donor),
      stringsAsFactors = FALSE
    )
  )
}

apply_donor_map_to_genotype <- function(
    genotype,
    donor_map,
    target_donors) {
  genotype <- as.matrix(genotype)
  target_donors <- as.character(target_donors)
  donor_ids <- rownames(genotype)
  required_map_columns <- c("target_donor", "source_donor")
  if (nrow(genotype) < 2L || ncol(genotype) < 1L ||
      is.null(donor_ids) || any(!nzchar(donor_ids)) ||
      anyDuplicated(donor_ids) || any(!is.finite(genotype)) ||
      is.null(colnames(genotype)) || any(!nzchar(colnames(genotype))) ||
      anyDuplicated(colnames(genotype)) || length(target_donors) < 2L ||
      any(!nzchar(target_donors)) || anyDuplicated(target_donors) ||
      !all(required_map_columns %in% names(donor_map))) {
    stop("Invalid genotype matrix, donor map, or target donors.")
  }
  map_target <- as.character(donor_map$target_donor)
  map_source <- as.character(donor_map$source_donor)
  if (nrow(donor_map) != length(target_donors) ||
      any(!nzchar(map_target)) || any(!nzchar(map_source)) ||
      anyDuplicated(map_target) || anyDuplicated(map_source) ||
      !setequal(map_target, target_donors) ||
      !setequal(map_source, target_donors) ||
      any(!target_donors %in% donor_ids)) {
    stop("The donor map is not a bijection over the target donors.")
  }
  source_donors <- map_source[match(target_donors, map_target)]
  permuted <- genotype[source_donors, , drop = FALSE]
  rownames(permuted) <- target_donors
  if (!isTRUE(all.equal(
    unname(crossprod(genotype[target_donors, , drop = FALSE])),
    unname(crossprod(permuted)),
    tolerance = 1e-12
  ))) {
    stop("The mapped genotype matrix failed its invariants.")
  }
  permuted
}

make_independent_time_donor_permutations <- function(
    donor_observation_matrix,
    time_grid,
    seed) {
  donor_observation_matrix <- as.matrix(donor_observation_matrix)
  time_grid <- as.numeric(time_grid)
  seed <- as.integer(seed)
  donor_ids <- rownames(donor_observation_matrix)
  if (!is.logical(donor_observation_matrix) ||
      nrow(donor_observation_matrix) < 2L ||
      ncol(donor_observation_matrix) != length(time_grid) ||
      is.null(donor_ids) || any(!nzchar(donor_ids)) ||
      anyDuplicated(donor_ids) || anyNA(donor_observation_matrix) ||
      any(colSums(donor_observation_matrix) < 2L) ||
      any(!is.finite(time_grid)) || anyDuplicated(time_grid) ||
      length(seed) != 1L || is.na(seed)) {
    stop("Invalid donor-observation matrix, time grid, or permutation seed.")
  }

  had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_random_seed) {
    previous_random_seed <- get(".Random.seed", envir = .GlobalEnv)
  }
  on.exit({
    if (had_random_seed) {
      assign(".Random.seed", previous_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  rows <- lapply(seq_along(time_grid), function(time_index) {
    observed_donors <- donor_ids[donor_observation_matrix[, time_index]]
    source_donors <- observed_donors
    for (attempt in seq_len(100L)) {
      source_donors <- sample(
        observed_donors,
        length(observed_donors),
        replace = FALSE
      )
      if (!identical(source_donors, observed_donors)) {
        break
      }
    }
    if (identical(source_donors, observed_donors)) {
      source_donors <- c(observed_donors[-1L], observed_donors[1L])
    }
    data.frame(
      time_index = time_index,
      time_value = time_grid[time_index],
      target_donor = observed_donors,
      source_donor = source_donors,
      fixed_point = observed_donors == source_donors,
      stringsAsFactors = FALSE
    )
  })
  donor_map <- do.call(rbind, rows)
  rownames(donor_map) <- NULL

  valid_time_maps <- vapply(seq_along(time_grid), function(time_index) {
    time_map <- donor_map[donor_map$time_index == time_index, , drop = FALSE]
    observed_donors <- donor_ids[donor_observation_matrix[, time_index]]
    nrow(time_map) == length(observed_donors) &&
      !anyDuplicated(time_map$target_donor) &&
      !anyDuplicated(time_map$source_donor) &&
      setequal(time_map$target_donor, observed_donors) &&
      setequal(time_map$source_donor, observed_donors) &&
      any(!time_map$fixed_point) &&
      all(time_map$fixed_point ==
            (time_map$target_donor == time_map$source_donor))
  }, logical(1))
  if (!all(valid_time_maps)) {
    stop("An independent time-specific donor permutation failed validation.")
  }
  donor_map
}

make_unit_specific_donor_block_permutations <- function(
    donor_ids,
    observation_patterns,
    unit_keys,
    seed) {
  donor_ids <- as.character(donor_ids)
  unit_keys <- as.character(unit_keys)
  pattern_names <- names(observation_patterns)
  observation_patterns <- setNames(
    as.character(observation_patterns),
    pattern_names
  )
  seed <- as.integer(seed)
  if (length(donor_ids) < 2L || any(!nzchar(donor_ids)) ||
      anyDuplicated(donor_ids) || length(unit_keys) < 1L ||
      any(!nzchar(unit_keys)) || anyDuplicated(unit_keys) ||
      is.null(names(observation_patterns)) ||
      any(!nzchar(names(observation_patterns))) ||
      anyDuplicated(names(observation_patterns)) ||
      !setequal(donor_ids, names(observation_patterns)) ||
      any(!nzchar(observation_patterns)) ||
      length(seed) != 1L || is.na(seed)) {
    stop("Invalid donors, observation patterns, unit keys, or seed.")
  }
  observation_patterns <- observation_patterns[donor_ids]
  pattern_groups <- split(donor_ids, observation_patterns)
  movable_groups <- pattern_groups[lengths(pattern_groups) > 1L]
  if (length(movable_groups) == 0L) {
    stop("The observation-pattern strata do not permit a non-identity map.")
  }

  had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_random_seed) {
    previous_random_seed <- get(".Random.seed", envir = .GlobalEnv)
  }
  on.exit({
    if (had_random_seed) {
      assign(".Random.seed", previous_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)

  unit_maps <- lapply(seq_along(unit_keys), function(unit_index) {
    source_donor <- setNames(donor_ids, donor_ids)
    for (attempt in seq_len(100L)) {
      for (group in pattern_groups) {
        source_donor[group] <- sample(group, length(group), replace = FALSE)
      }
      if (any(source_donor[donor_ids] != donor_ids)) {
        break
      }
    }
    if (all(source_donor[donor_ids] == donor_ids)) {
      group <- movable_groups[[1L]]
      source_donor[group] <- c(group[-1L], group[1L])
    }
    data.frame(
      unit_index = unit_index,
      unit_key = unit_keys[unit_index],
      target_donor = donor_ids,
      source_donor = unname(source_donor[donor_ids]),
      fixed_point = donor_ids == unname(source_donor[donor_ids]),
      observation_pattern = unname(observation_patterns[donor_ids]),
      stringsAsFactors = FALSE
    )
  })
  donor_map <- do.call(rbind, unit_maps)
  rownames(donor_map) <- NULL

  source_patterns <- unname(
    observation_patterns[donor_map$source_donor]
  )
  valid_unit_maps <- vapply(seq_along(unit_keys), function(unit_index) {
    unit_map <- donor_map[
      donor_map$unit_index == unit_index,
      ,
      drop = FALSE
    ]
    nrow(unit_map) == length(donor_ids) &&
      identical(unit_map$unit_key, rep(unit_keys[unit_index], length(donor_ids))) &&
      identical(unit_map$target_donor, donor_ids) &&
      !anyDuplicated(unit_map$source_donor) &&
      setequal(unit_map$source_donor, donor_ids) &&
      any(!unit_map$fixed_point)
  }, logical(1))
  if (!all(valid_unit_maps) ||
      !identical(source_patterns, donor_map$observation_pattern) ||
      any(donor_map$fixed_point !=
            (donor_map$target_donor == donor_map$source_donor))) {
    stop("A unit-specific donor-block permutation failed validation.")
  }
  donor_map
}

make_signal_stripped_residual_block_null <- function(
    expression,
    genotype,
    covariates,
    source_rows,
    tolerance = 1e-10) {
  expression <- as.matrix(expression)
  genotype <- as.matrix(genotype)
  covariates <- as.matrix(covariates)
  storage.mode(expression) <- "double"
  storage.mode(genotype) <- "double"
  storage.mode(covariates) <- "double"
  source_row_matrix <- if (is.null(dim(source_rows))) {
    matrix(
      rep(as.integer(source_rows), ncol(expression)),
      nrow = nrow(expression),
      ncol = ncol(expression)
    )
  } else {
    source_rows <- as.matrix(source_rows)
    if (any(!is.finite(source_rows)) ||
        any(source_rows != as.integer(source_rows))) {
      stop("The source-row matrix contains invalid indices.")
    }
    matrix(
      as.integer(source_rows),
      nrow = nrow(source_rows),
      ncol = ncol(source_rows)
    )
  }
  tolerance <- as.numeric(tolerance)
  valid_source_columns <- if (identical(
    dim(source_row_matrix),
    dim(expression)
  )) {
    vapply(seq_len(ncol(source_row_matrix)), function(unit_index) {
      identical(
        sort(source_row_matrix[, unit_index]),
        seq_len(nrow(expression))
      )
    }, logical(1))
  } else {
    FALSE
  }
  if (!identical(dim(expression), dim(genotype)) ||
      nrow(expression) < 4L || ncol(expression) < 1L ||
      nrow(covariates) != nrow(expression) || ncol(covariates) < 1L ||
      any(!is.finite(expression)) || any(!is.finite(genotype)) ||
      any(!is.finite(covariates)) ||
      anyNA(source_row_matrix) || !all(valid_source_columns) ||
      length(tolerance) != 1L || !is.finite(tolerance) || tolerance <= 0) {
    stop("Invalid signal-stripped residual-block inputs.")
  }

  observed_fit <- fit_many_genotype_regressions(
    expression,
    genotype,
    covariates
  )
  fitted_genotype <- sweep(genotype, 2L, observed_fit$beta, `*`)
  nuisance_fitted <- expression - observed_fit$residual - fitted_genotype
  source_full_model_residual <- observed_fit$residual
  residual_genotype_cross_product <- colSums(
    observed_fit$genotype_residual * source_full_model_residual
  )
  residual_genotype_scale <- sqrt(
    colSums(observed_fit$genotype_residual^2) *
      colSums(source_full_model_residual^2)
  )
  if (any(!is.finite(residual_genotype_scale)) ||
      any(residual_genotype_scale <= 0)) {
    stop("The source full-model residual has an invalid scale.")
  }
  residual_genotype_correlation <- abs(
    residual_genotype_cross_product / residual_genotype_scale
  )
  maximum_residual_genotype_correlation <- max(
    residual_genotype_correlation
  )
  if (!is.finite(maximum_residual_genotype_correlation) ||
      maximum_residual_genotype_correlation > tolerance) {
    stop("The source full-model residual is not orthogonal to genotype.")
  }

  source_linear_indices <- as.vector(source_row_matrix) + rep(
    (seq_len(ncol(expression)) - 1L) * nrow(expression),
    each = nrow(expression)
  )
  permuted_full_model_residual <- matrix(
    source_full_model_residual[source_linear_indices],
    nrow = nrow(expression),
    ncol = ncol(expression),
    dimnames = dimnames(source_full_model_residual)
  )
  null_expression <- nuisance_fitted + permuted_full_model_residual
  null_fit <- fit_many_genotype_regressions(
    null_expression,
    genotype,
    covariates
  )
  if (observed_fit$residual_df != null_fit$residual_df) {
    stop("Observed and signal-stripped null residual degrees of freedom differ.")
  }

  list(
    observed_fit = observed_fit,
    null_fit = null_fit,
    null_expression = null_expression,
    nuisance_fitted = nuisance_fitted,
    fitted_genotype = fitted_genotype,
    source_full_model_residual = source_full_model_residual,
    permuted_full_model_residual = permuted_full_model_residual,
    residual_genotype_correlation = residual_genotype_correlation,
    maximum_residual_genotype_correlation =
      maximum_residual_genotype_correlation
  )
}

convert_raw_to_original_t_adjusted_se <- function(beta_hat,
                                                   raw_se,
                                                   residual_df) {
  beta_hat <- as.matrix(beta_hat)
  raw_se <- as.matrix(raw_se)
  residual_df <- as.numeric(residual_df)
  if (!identical(dim(beta_hat), dim(raw_se)) ||
      length(residual_df) != ncol(beta_hat) ||
      any(!is.finite(beta_hat)) || any(!is.finite(raw_se)) ||
      any(raw_se <= 0) || any(!is.finite(residual_df)) ||
      any(residual_df <= 0)) {
    stop("Invalid beta, raw-SE, or residual-df inputs.")
  }
  adjusted <- raw_se
  for (time_index in seq_len(ncol(beta_hat))) {
    t_value <- beta_hat[, time_index] / raw_se[, time_index]
    two_sided_p <- 2 * stats::pt(
      abs(t_value),
      df = residual_df[time_index],
      lower.tail = FALSE
    )
    z_value <- stats::qnorm(1 - two_sided_p / 2)
    adjusted[, time_index] <- abs(beta_hat[, time_index]) / abs(z_value)
    near_zero <- abs(t_value) < 1e-8 |
      !is.finite(adjusted[, time_index]) |
      adjusted[, time_index] <= 0
    if (any(near_zero)) {
      adjusted[near_zero, time_index] <- raw_se[near_zero, time_index] *
        stats::dnorm(0) / stats::dt(0, df = residual_df[time_index])
    }
  }
  if (any(!is.finite(adjusted)) || any(adjusted <= 0)) {
    stop("The original t-based SE adjustment produced invalid values.")
  }
  adjusted
}

safe_divide <- function(numerator, denominator) {
  if (length(numerator) != 1L || length(denominator) != 1L ||
      !is.finite(numerator) || !is.finite(denominator) || denominator == 0) {
    return(NA_real_)
  }
  numerator / denominator
}

summarize_matched_null_calibration <- function(lfdr,
                                                group,
                                                pi0_merged,
                                                fit_stage,
                                                alpha = 0.05) {
  lfdr <- as.numeric(lfdr)
  group <- as.character(group)
  pi0_merged <- as.numeric(pi0_merged)
  fit_stage <- as.character(fit_stage)
  expected_groups <- c("target", "permuted_null")
  if (length(lfdr) != length(group) || length(lfdr) < 2L ||
      any(!is.finite(lfdr)) || any(lfdr < 0 | lfdr > 1) ||
      !setequal(unique(group), expected_groups) ||
      length(pi0_merged) != 1L || !is.finite(pi0_merged) ||
      pi0_merged < 0 || pi0_merged > 1 || length(fit_stage) != 1L ||
      is.na(fit_stage) || !nzchar(fit_stage)) {
    stop("Invalid matched-null calibration inputs.")
  }
  calls <- seq_along(lfdr) %in% cumulative_lfdr_calls(lfdr, alpha = alpha)
  target <- group == "target"
  permuted_null <- group == "permuted_null"
  n_target <- sum(target)
  n_null <- sum(permuted_null)
  n_total <- length(lfdr)
  target_calls <- sum(calls & target)
  null_calls <- sum(calls & permuted_null)
  total_calls <- sum(calls)
  target_call_rate <- target_calls / n_target
  null_call_rate <- null_calls / n_null
  merged_call_rate <- total_calls / n_total
  design_pi0_lower_bound <- n_null / n_total
  known_null_discovery_fraction <- safe_divide(null_calls, total_calls)
  pi0_target_unbounded <- (n_total * pi0_merged - n_null) / n_target
  pi0_target_bounded <- min(1, max(0, pi0_target_unbounded))
  pi0_target_valid <- pi0_target_unbounded >= 0 && pi0_target_unbounded <= 1
  data.frame(
    fit_stage = fit_stage,
    alpha = alpha,
    n_target = n_target,
    n_permuted_null = n_null,
    n_total = n_total,
    target_calls = target_calls,
    permuted_null_calls = null_calls,
    total_calls = total_calls,
    target_call_rate = target_call_rate,
    permuted_null_call_rate = null_call_rate,
    merged_call_rate = merged_call_rate,
    pi0_merged = pi0_merged,
    design_pi0_lower_bound = design_pi0_lower_bound,
    pi0_merged_below_design_lower_bound =
      pi0_merged < design_pi0_lower_bound,
    pi0_target_unbounded = pi0_target_unbounded,
    pi0_target_bounded = pi0_target_bounded,
    pi0_target_valid = pi0_target_valid,
    known_null_discovery_fraction = known_null_discovery_fraction,
    scaled_fdr_merged_from_estimated_pi0 = safe_divide(
      pi0_merged * null_call_rate,
      merged_call_rate
    ),
    scaled_fdr_merged_from_design_lower_bound = safe_divide(
      design_pi0_lower_bound * null_call_rate,
      merged_call_rate
    ),
    post_selection_fdr_target_from_pi0 = if (pi0_target_valid) {
      safe_divide(
        pi0_target_unbounded * null_call_rate,
        target_call_rate
      )
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

summarize_group_lfdr <- function(lfdr, group, fit_stage) {
  lfdr <- as.numeric(lfdr)
  group <- as.character(group)
  probabilities <- c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1)
  rows <- lapply(sort(unique(group)), function(group_name) {
    values <- lfdr[group == group_name]
    quantiles <- stats::quantile(
      values,
      probs = probabilities,
      names = FALSE,
      type = 8
    )
    data.frame(
      fit_stage = fit_stage,
      group = group_name,
      probability = probabilities,
      lfdr_quantile = quantiles,
      group_mean_lfdr = mean(values),
      group_sd_lfdr = stats::sd(values),
      n_units = length(values),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

refit_fash_from_likelihood <- function(source_fit,
                                       data_list,
                                       se_list,
                                       likelihood_matrix,
                                       unit_keys,
                                       penalty = source_fit$settings$penalty) {
  if (!requireNamespace("fashr", quietly = TRUE)) {
    stop("The fashr package is required.")
  }
  unit_keys <- as.character(unit_keys)
  likelihood_matrix <- as.matrix(likelihood_matrix)
  if (length(data_list) != length(se_list) ||
      length(data_list) != nrow(likelihood_matrix) ||
      length(unit_keys) != length(data_list) || any(!nzchar(unit_keys)) ||
      anyDuplicated(unit_keys) || ncol(likelihood_matrix) !=
        length(source_fit$psd_grid) || anyNA(likelihood_matrix) ||
      any(is.nan(likelihood_matrix)) || any(likelihood_matrix == Inf)) {
    stop("Invalid merged FASH likelihood inputs.")
  }
  names(data_list) <- unit_keys
  names(se_list) <- unit_keys
  rownames(likelihood_matrix) <- unit_keys
  empirical_bayes <- fashr::fash_eb_est(
    L_matrix = likelihood_matrix,
    grid = source_fit$psd_grid,
    penalty = penalty
  )
  rownames(empirical_bayes$posterior_weight) <- unit_keys
  null_column <- which(empirical_bayes$prior_weight$psd == 0)
  if (length(null_column) != 1L) {
    stop("The merged FASH fit does not have exactly one exact-null component.")
  }
  lfdr <- empirical_bayes$posterior_weight[, null_column]
  names(lfdr) <- unit_keys
  structure(
    list(
      prior_weights = empirical_bayes$prior_weight,
      posterior_weights = empirical_bayes$posterior_weight,
      psd_grid = source_fit$psd_grid,
      lfdr = lfdr,
      settings = source_fit$settings,
      fash_data = list(data_list = data_list, S = se_list, Omega = NULL),
      L_matrix = likelihood_matrix,
      eb_result = empirical_bayes
    ),
    class = "fash"
  )
}

extract_pi0 <- function(fash_fit) {
  null_row <- which(fash_fit$prior_weights$psd == 0)
  if (length(null_row) != 1L) {
    stop("The FASH fit does not contain exactly one null prior weight.")
  }
  as.numeric(fash_fit$prior_weights$prior_weight[null_row])
}
