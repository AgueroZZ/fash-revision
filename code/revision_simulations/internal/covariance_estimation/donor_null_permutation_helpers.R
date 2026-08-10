# Helpers for the donor-level null-permutation correlation pilot.

read_selected_vcf_dosages <- function(vcf_path,
                                      variant_ids,
                                      chunk_size = 100000L) {
  variant_ids <- as.character(variant_ids)
  chunk_size <- as.integer(chunk_size)
  if (!file.exists(vcf_path) || length(variant_ids) == 0L ||
      any(!nzchar(variant_ids)) || anyDuplicated(variant_ids) ||
      is.na(chunk_size) || chunk_size < 100L) {
    stop("Invalid VCF path, variant IDs, or chunk size.")
  }

  connection <- gzfile(vcf_path, open = "rt")
  on.exit(close(connection), add = TRUE)
  sample_ids <- NULL
  found <- vector("list", length(variant_ids))
  names(found) <- variant_ids
  requested <- setNames(seq_along(variant_ids), variant_ids)

  repeat {
    lines <- readLines(connection, n = chunk_size, warn = FALSE)
    if (length(lines) == 0L) {
      break
    }

    chrom_header <- lines[startsWith(lines, "#CHROM\t")]
    if (length(chrom_header) > 1L) {
      stop("The VCF contains more than one #CHROM header.")
    }
    if (length(chrom_header) == 1L) {
      header_fields <- strsplit(chrom_header, "\t", fixed = TRUE)[[1L]]
      if (length(header_fields) < 10L) {
        stop("The VCF header does not contain sample columns.")
      }
      sample_ids <- header_fields[10:length(header_fields)]
      if (any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
        stop("The VCF sample IDs are empty or duplicated.")
      }
    }

    data_lines <- lines[!startsWith(lines, "#")]
    if (length(data_lines) == 0L) {
      next
    }
    line_ids <- sub(
      "^[^\t]*\t[^\t]*\t([^\t]*)\t.*$",
      "\\1",
      data_lines
    )
    keep <- line_ids %in% variant_ids
    if (!any(keep)) {
      next
    }

    matched_lines <- data_lines[keep]
    for (line in matched_lines) {
      fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      variant_id <- fields[3L]
      output_index <- requested[[variant_id]]
      if (!is.null(found[[output_index]])) {
        stop("The requested VCF variant appears more than once: ", variant_id)
      }
      if (is.null(sample_ids) || length(fields) != 9L + length(sample_ids)) {
        stop("A requested VCF row has an unexpected number of fields: ", variant_id)
      }
      format_names <- strsplit(fields[9L], ":", fixed = TRUE)[[1L]]
      dosage_index <- match("DS", format_names)
      if (is.na(dosage_index)) {
        stop("A requested VCF row does not contain the DS field: ", variant_id)
      }
      sample_fields <- strsplit(fields[10:length(fields)], ":", fixed = TRUE)
      dosage <- vapply(sample_fields, function(values) {
        if (length(values) < dosage_index) {
          return(NA_real_)
        }
        suppressWarnings(as.numeric(values[dosage_index]))
      }, numeric(1))
      if (any(!is.finite(dosage))) {
        stop("A requested VCF row contains a missing or non-finite dosage: ",
             variant_id)
      }
      found[[output_index]] <- dosage
    }
  }

  missing <- variant_ids[vapply(found, is.null, logical(1))]
  if (length(missing) > 0L) {
    stop(
      "The VCF is missing ",
      length(missing),
      " requested variants; first missing ID: ",
      missing[1L]
    )
  }
  dosage <- do.call(cbind, found)
  rownames(dosage) <- sample_ids
  colnames(dosage) <- variant_ids
  storage.mode(dosage) <- "double"
  dosage
}

make_covariate_projection <- function(covariates) {
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  if (nrow(covariates) < 4L || ncol(covariates) < 1L ||
      any(!is.finite(covariates))) {
    stop("The covariate matrix is invalid.")
  }
  design <- cbind(intercept = 1, covariates)
  design_qr <- qr(design)
  if (design_qr$rank != ncol(design)) {
    stop("The covariate design is rank deficient.")
  }
  q_matrix <- qr.Q(design_qr)
  residualizer <- diag(nrow(design)) - tcrossprod(q_matrix)
  list(
    design = design,
    rank = design_qr$rank,
    residualizer = residualizer
  )
}

make_shared_donor_block_permutation <- function(donor_ids,
                                                 observation_patterns) {
  donor_ids <- as.character(donor_ids)
  pattern_names <- names(observation_patterns)
  observation_patterns <- setNames(
    as.character(observation_patterns),
    pattern_names
  )
  if (length(donor_ids) < 2L || any(!nzchar(donor_ids)) ||
      anyDuplicated(donor_ids) || is.null(names(observation_patterns)) ||
      any(!nzchar(names(observation_patterns))) ||
      anyDuplicated(names(observation_patterns)) ||
      !setequal(donor_ids, names(observation_patterns)) ||
      any(!nzchar(observation_patterns))) {
    stop("Invalid donor IDs or observation patterns.")
  }
  observation_patterns <- observation_patterns[donor_ids]
  pattern_groups <- split(donor_ids, observation_patterns)
  source_donor <- setNames(donor_ids, donor_ids)
  for (group in pattern_groups) {
    source_donor[group] <- sample(group, length(group), replace = FALSE)
  }
  if (any(observation_patterns[source_donor] != observation_patterns)) {
    stop("A donor block was mapped across incompatible missingness patterns.")
  }
  source_donor
}

fit_many_genotype_regressions <- function(expression,
                                          genotype,
                                          covariates) {
  expression <- as.matrix(expression)
  genotype <- as.matrix(genotype)
  storage.mode(expression) <- "double"
  storage.mode(genotype) <- "double"
  if (!identical(dim(expression), dim(genotype)) ||
      nrow(expression) < 4L || ncol(expression) < 1L ||
      any(!is.finite(expression)) || any(!is.finite(genotype))) {
    stop("Expression and genotype must be matching finite matrices.")
  }
  projection <- make_covariate_projection(covariates)
  residual_df <- nrow(expression) - projection$rank - 1L
  if (residual_df < 2L) {
    stop("The genotype regression has fewer than two residual degrees of freedom.")
  }
  expression_residual <- projection$residualizer %*% expression
  genotype_residual <- projection$residualizer %*% genotype
  denominator <- colSums(genotype_residual^2)
  tolerance <- 1e-12 * pmax(1, colSums(genotype^2))
  if (any(!is.finite(denominator)) || any(denominator <= tolerance)) {
    stop("At least one genotype column is constant or explained by the covariates.")
  }
  beta <- colSums(genotype_residual * expression_residual) / denominator
  fitted_genotype <- sweep(genotype_residual, 2L, beta, `*`)
  residual <- expression_residual - fitted_genotype
  residual_variance <- colSums(residual^2) / residual_df
  standard_error <- sqrt(residual_variance / denominator)
  if (any(!is.finite(beta)) || any(!is.finite(standard_error)) ||
      any(standard_error <= 0)) {
    stop("The vectorized genotype regression produced invalid estimates.")
  }
  list(
    beta = beta,
    standard_error = standard_error,
    residual_df = residual_df,
    expression_residual = expression_residual,
    genotype_residual = genotype_residual,
    residual = residual,
    projection = projection
  )
}

fit_residualized_genotype_regressions <- function(expression_residual,
                                                  genotype,
                                                  residualizer,
                                                  covariate_rank) {
  expression_residual <- as.matrix(expression_residual)
  genotype <- as.matrix(genotype)
  residualizer <- as.matrix(residualizer)
  covariate_rank <- as.integer(covariate_rank)
  storage.mode(expression_residual) <- "double"
  storage.mode(genotype) <- "double"
  if (!identical(dim(expression_residual), dim(genotype)) ||
      nrow(residualizer) != nrow(genotype) ||
      ncol(residualizer) != nrow(genotype) ||
      any(!is.finite(expression_residual)) || any(!is.finite(genotype)) ||
      any(!is.finite(residualizer)) || length(covariate_rank) != 1L ||
      is.na(covariate_rank)) {
    stop("Invalid residualized-regression inputs.")
  }
  residual_df <- nrow(genotype) - covariate_rank - 1L
  if (residual_df < 2L) {
    stop("The residualized genotype regression has insufficient degrees of freedom.")
  }
  genotype_residual <- residualizer %*% genotype
  denominator <- colSums(genotype_residual^2)
  tolerance <- 1e-12 * pmax(1, colSums(genotype^2))
  if (any(!is.finite(denominator)) || any(denominator <= tolerance)) {
    stop("At least one residualized genotype column has zero information.")
  }
  beta <- colSums(genotype_residual * expression_residual) / denominator
  residual <- expression_residual - sweep(genotype_residual, 2L, beta, `*`)
  standard_error <- sqrt(
    (colSums(residual^2) / residual_df) / denominator
  )
  if (any(!is.finite(beta)) || any(!is.finite(standard_error)) ||
      any(standard_error <= 0)) {
    stop("The residualized genotype regression produced invalid estimates.")
  }
  list(
    beta = beta,
    standard_error = standard_error,
    residual_df = residual_df
  )
}

convert_raw_to_t_adjusted_se <- function(beta_hat,
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
    absolute_t <- abs(t_value)
    probability <- stats::pt(
      absolute_t,
      df = residual_df[time_index],
      lower.tail = TRUE
    )
    z_value <- stats::qnorm(probability)
    ratio <- absolute_t / z_value
    near_zero <- absolute_t < 1e-8 | !is.finite(ratio)
    if (any(near_zero)) {
      ratio[near_zero] <- stats::dnorm(0) /
        stats::dt(0, df = residual_df[time_index])
    }
    adjusted[, time_index] <- raw_se[, time_index] * ratio
  }
  if (any(!is.finite(adjusted)) || any(adjusted <= 0)) {
    stop("The t-based SE adjustment produced invalid values.")
  }
  adjusted
}

estimate_ordinary_pairwise_correlation <- function(beta_hat, se) {
  beta_hat <- validate_finite_matrix(beta_hat, "beta_hat")
  se <- validate_finite_matrix(se, "se", positive = TRUE)
  if (!identical(dim(beta_hat), dim(se))) {
    stop("beta_hat and se must have the same dimensions.")
  }
  n_time <- ncol(beta_hat)
  correlation <- diag(n_time)
  dimnames(correlation) <- list(colnames(beta_hat), colnames(beta_hat))
  for (time_a in seq_len(n_time - 1L)) {
    for (time_b in (time_a + 1L):n_time) {
      variance_sum <- se[, time_a]^2 + se[, time_b]^2
      covariance_scale <- 2 * se[, time_a] * se[, time_b]
      unit_estimate <- (
        variance_sum -
          (beta_hat[, time_a] - beta_hat[, time_b])^2
      ) / covariance_scale
      estimate <- mean(unit_estimate)
      correlation[time_a, time_b] <- correlation[time_b, time_a] <- estimate
    }
  }
  correlation
}

summarize_raw_correlation_matrix <- function(correlation) {
  correlation <- validate_finite_matrix(correlation, "correlation")
  if (nrow(correlation) != ncol(correlation) ||
      !isTRUE(all.equal(correlation, t(correlation), tolerance = 1e-10))) {
    stop("correlation must be a finite symmetric square matrix.")
  }
  lag_correlation <- lag_average_correlation(correlation)
  off_diagonal <- correlation[row(correlation) != col(correlation)]
  eigenvalues <- eigen(correlation, symmetric = TRUE, only.values = TRUE)$values
  data.frame(
    lag1 = lag_correlation[1L],
    mean_lags_9_15 = if (length(lag_correlation) >= 15L) {
      mean(lag_correlation[9:15])
    } else {
      NA_real_
    },
    lag15 = if (length(lag_correlation) >= 15L) {
      lag_correlation[15L]
    } else {
      NA_real_
    },
    mean_off_diagonal = mean(off_diagonal),
    minimum_off_diagonal = min(off_diagonal),
    maximum_off_diagonal = max(off_diagonal),
    minimum_eigenvalue = min(eigenvalues),
    n_negative_eigenvalues = sum(eigenvalues < -1e-8),
    stringsAsFactors = FALSE
  )
}
