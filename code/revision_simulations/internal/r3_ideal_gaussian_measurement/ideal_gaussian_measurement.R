# Pure helpers for the R3 ideal-Gaussian measurement isolation experiment.

validate_numeric_matrix <- function(x, name) {
  x <- as.matrix(x)
  if (!is.numeric(x) || length(dim(x)) != 2L || any(dim(x) < 1L)) {
    stop(name, " must be a non-empty numeric matrix.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop(name, " must contain only finite values.", call. = FALSE)
  }
  x
}

validate_single_seed <- function(seed, name = "seed") {
  if (length(seed) != 1L || !is.finite(seed) || seed != round(seed)) {
    stop(name, " must be one finite integer.", call. = FALSE)
  }
  max_seed <- .Machine$integer.max - 1
  seed <- as.double(seed)
  as.integer(((seed - 1) %% max_seed) + 1)
}

ideal_gaussian_measurement_seed <- function(seed) {
  seed <- validate_single_seed(seed)
  max_seed <- .Machine$integer.max - 1
  offset <- 23003
  as.integer(((as.double(seed) + offset - 1) %% max_seed) + 1)
}

with_preserved_seed <- function(seed, code) {
  seed <- validate_single_seed(seed)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    previous_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_seed) {
      assign(".Random.seed", previous_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(code)
}

validate_regression_summary_for_adjusted_se <- function(regression_summary,
                                                        adjusted_se) {
  if (is.null(regression_summary)) return(invisible(TRUE))
  required <- c(
    "beta_hat", "se", "se_uncorrected", "df", "rank",
    "apply_t_se_correction"
  )
  if (!is.list(regression_summary) ||
      !all(required %in% names(regression_summary))) {
    stop(
      "regression_summary does not contain the required R3 regression fields.",
      call. = FALSE
    )
  }
  if (!isTRUE(regression_summary$apply_t_se_correction)) {
    stop("regression_summary must contain t-adjusted standard errors.", call. = FALSE)
  }
  if (!identical(as.matrix(regression_summary$se), adjusted_se)) {
    stop("regression_summary$se does not match adjusted_se.", call. = FALSE)
  }
  invisible(TRUE)
}

simulate_ideal_gaussian_eqtl_summary <- function(true_beta,
                                                 adjusted_se,
                                                 seed,
                                                 regression_summary = NULL) {
  true_beta <- validate_numeric_matrix(true_beta, "true_beta")
  adjusted_se <- validate_numeric_matrix(adjusted_se, "adjusted_se")
  if (!identical(dim(true_beta), dim(adjusted_se))) {
    stop("true_beta and adjusted_se must have identical dimensions.", call. = FALSE)
  }
  if (any(adjusted_se <= 0)) {
    stop("adjusted_se must be strictly positive.", call. = FALSE)
  }
  if (!is.null(rownames(true_beta)) && !is.null(rownames(adjusted_se)) &&
      !identical(rownames(true_beta), rownames(adjusted_se))) {
    stop("true_beta and adjusted_se row names must match.", call. = FALSE)
  }
  validate_regression_summary_for_adjusted_se(
    regression_summary = regression_summary,
    adjusted_se = adjusted_se
  )

  seed <- validate_single_seed(seed)
  standard_normal_draw <- with_preserved_seed(seed, {
    matrix(
      stats::rnorm(length(true_beta)),
      nrow = nrow(true_beta),
      ncol = ncol(true_beta),
      dimnames = dimnames(true_beta)
    )
  })
  dimnames(adjusted_se) <- dimnames(true_beta)
  beta_hat <- true_beta + adjusted_se * standard_normal_draw
  standardized_error <- (beta_hat - true_beta) / adjusted_se

  if (!isTRUE(all.equal(
    standardized_error,
    standard_normal_draw,
    tolerance = 1e-13,
    check.attributes = TRUE
  ))) {
    stop("The ideal-Gaussian standardized-error identity failed.", call. = FALSE)
  }

  list(
    beta_hat = beta_hat,
    se = adjusted_se,
    standardized_error = standardized_error,
    standard_normal_draw = standard_normal_draw,
    measurement_seed = seed,
    apply_t_se_correction = FALSE,
    se_source = "R3 final t-adjusted SE treated as fixed"
  )
}

validate_standardized_error_matrix <- function(z) {
  z <- validate_numeric_matrix(z, "z")
  if (nrow(z) < 2L || ncol(z) < 2L) {
    stop("z must contain at least two units and two time points.", call. = FALSE)
  }
  z
}

standardized_error_metadata <- function(seed,
                                        truth_mechanism,
                                        observation_model) {
  seed <- validate_single_seed(seed)
  if (length(truth_mechanism) != 1L || !nzchar(truth_mechanism)) {
    stop("truth_mechanism must be one non-empty string.", call. = FALSE)
  }
  if (length(observation_model) != 1L || !nzchar(observation_model)) {
    stop("observation_model must be one non-empty string.", call. = FALSE)
  }
  list(
    seed = seed,
    truth_mechanism = as.character(truth_mechanism),
    observation_model = as.character(observation_model)
  )
}

summarize_standardized_error_matrix <- function(z,
                                                seed,
                                                truth_mechanism,
                                                observation_model) {
  z <- validate_standardized_error_matrix(z)
  metadata <- standardized_error_metadata(
    seed, truth_mechanism, observation_model
  )
  values <- as.numeric(z)
  center <- mean(values)
  scale <- stats::sd(values)
  if (!is.finite(scale) || scale <= 0) {
    stop("z must have positive finite standard deviation.", call. = FALSE)
  }
  standardized <- (values - center) / scale
  quantiles <- stats::quantile(
    values,
    probs = c(0.005, 0.025, 0.25, 0.5, 0.75, 0.975, 0.995),
    names = FALSE,
    type = 8
  )
  data.frame(
    seed = metadata$seed,
    truth_mechanism = metadata$truth_mechanism,
    observation_model = metadata$observation_model,
    n = length(values),
    mean = center,
    sd = scale,
    skewness = mean(standardized^3),
    excess_kurtosis = mean(standardized^4) - 3,
    q005 = quantiles[[1L]],
    q025 = quantiles[[2L]],
    q25 = quantiles[[3L]],
    median = quantiles[[4L]],
    q75 = quantiles[[5L]],
    q975 = quantiles[[6L]],
    q995 = quantiles[[7L]],
    tail_abs_1_96 = mean(abs(values) > 1.96),
    tail_abs_3 = mean(abs(values) > 3),
    stringsAsFactors = FALSE
  )
}

summarize_standardized_error_histogram <- function(
    z,
    seed,
    truth_mechanism,
    observation_model,
    breaks = seq(-6, 6, by = 0.1)) {
  z <- validate_standardized_error_matrix(z)
  metadata <- standardized_error_metadata(
    seed, truth_mechanism, observation_model
  )
  breaks <- as.numeric(breaks)
  if (length(breaks) < 2L || any(!is.finite(breaks)) ||
      any(diff(breaks) <= 0)) {
    stop("breaks must be increasing finite values.", call. = FALSE)
  }
  complete_breaks <- c(-Inf, breaks, Inf)
  bin_index <- cut(
    as.numeric(z),
    breaks = complete_breaks,
    include.lowest = TRUE,
    right = FALSE,
    labels = FALSE
  )
  count <- tabulate(bin_index, nbins = length(complete_breaks) - 1L)
  left <- head(complete_breaks, -1L)
  right <- tail(complete_breaks, -1L)
  bin_label <- ifelse(
    is.infinite(left),
    paste0("< ", format(right, trim = TRUE)),
    ifelse(
      is.infinite(right),
      paste0(">= ", format(left, trim = TRUE)),
      paste0(
        "[", format(left, trim = TRUE), ", ",
        format(right, trim = TRUE), ")"
      )
    )
  )
  data.frame(
    seed = metadata$seed,
    truth_mechanism = metadata$truth_mechanism,
    observation_model = metadata$observation_model,
    bin_index = seq_along(count),
    bin_left = left,
    bin_right = right,
    bin_label = bin_label,
    count = count,
    probability = count / sum(count),
    stringsAsFactors = FALSE
  )
}

summarize_standardized_error_lag_correlation <- function(
    z,
    seed,
    truth_mechanism,
    observation_model) {
  z <- validate_standardized_error_matrix(z)
  metadata <- standardized_error_metadata(
    seed, truth_mechanism, observation_model
  )
  correlation <- stats::cor(z)
  if (any(!is.finite(correlation))) {
    stop("z produced non-finite time-point correlations.", call. = FALSE)
  }
  output <- lapply(seq_len(ncol(z) - 1L), function(lag) {
    row_index <- seq_len(ncol(z) - lag)
    values <- correlation[cbind(row_index, row_index + lag)]
    data.frame(
      seed = metadata$seed,
      truth_mechanism = metadata$truth_mechanism,
      observation_model = metadata$observation_model,
      lag = lag,
      n_pairs = length(values),
      mean_correlation = mean(values),
      min_correlation = min(values),
      max_correlation = max(values),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, output)
}

serialized_object_md5 <- function(object) {
  path <- tempfile("fash-r3-ideal-md5-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(object, path, version = 3)
  unname(as.character(tools::md5sum(path))[[1L]])
}
