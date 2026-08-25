# Reusable simulation and inference functions for the FASH JASA revision.
#
# This script is designed to be sourced from the workflowr project root:
#   source("code/revision_simulations/shared/simulation_functions.R")
#
# It keeps the Appendix B simulation structure, but exposes reusable functions
# for reviewer-driven scenarios: FASH calibration, Strober-style parametric
# comparators, real-data-like effect shapes, correlated errors, and compact
# summary tables.

default_revision_grid <- function(log_prec = seq(0, 10, by = 0.2)) {
  sort(c(0, exp(-0.5 * log_prec)))
}

revision_output_dirs <- function(base_dir = "output/revision_simulations") {
  dirs <- list(
    base = base_dir,
    raw = file.path(base_dir, "raw"),
    summary = file.path(base_dir, "summary"),
    figures = file.path(base_dir, "figures")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
  dirs
}

numeric_vector_label <- function(x) {
  paste(gsub("\\.", "p", format(x, trim = TRUE, scientific = FALSE)), collapse = "-")
}

make_time_grid <- function(n_time = 16, start = 0, end = 15) {
  seq(start, end, length.out = n_time)
}

make_ar1_correlation <- function(n, rho) {
  idx <- seq_len(n)
  rho ^ abs(outer(idx, idx, "-"))
}

make_lag1_correlation <- function(n, rho) {
  if (length(n) != 1L || !is.finite(n) || n < 2 || n != as.integer(n)) {
    stop("n must be an integer greater than one.")
  }
  if (length(rho) != 1L || !is.finite(rho) || abs(rho) >= 0.5) {
    stop("rho must be finite and satisfy abs(rho) < 0.5.")
  }

  n <- as.integer(n)
  correlation <- diag(n)
  correlation[row(correlation) == col(correlation) + 1L] <- rho
  correlation[col(correlation) == row(correlation) + 1L] <- rho
  correlation
}

validate_time_correlation <- function(correlation,
                                      n_time,
                                      tolerance = 1e-8) {
  if (length(n_time) != 1L || !is.finite(n_time) ||
      n_time < 2 || n_time != as.integer(n_time)) {
    stop("n_time must be an integer greater than one.")
  }
  n_time <- as.integer(n_time)
  correlation <- as.matrix(correlation)
  storage.mode(correlation) <- "numeric"

  if (!identical(dim(correlation), c(n_time, n_time))) {
    stop("The time-correlation matrix dimensions must match n_time.")
  }
  if (any(!is.finite(correlation))) {
    stop("The time-correlation matrix must contain only finite values.")
  }
  if (max(abs(correlation - t(correlation))) > tolerance) {
    stop("The time-correlation matrix must be symmetric.")
  }
  if (max(abs(diag(correlation) - 1)) > tolerance) {
    stop("The time-correlation matrix must have a unit diagonal.")
  }

  correlation <- (correlation + t(correlation)) / 2
  diag(correlation) <- 1
  minimum_eigenvalue <- min(eigen(
    correlation,
    symmetric = TRUE,
    only.values = TRUE
  )$values)
  if (minimum_eigenvalue <= tolerance) {
    stop("The time-correlation matrix must be positive definite.")
  }
  correlation
}

make_exchangeable_correlation <- function(n, rho) {
  mat <- matrix(rho, nrow = n, ncol = n)
  diag(mat) <- 1
  mat
}

sample_effect_classes <- function(J, class_probs, exact_class_counts = TRUE) {
  if (!exact_class_counts) {
    return(sample(names(class_probs), size = J, replace = TRUE, prob = class_probs))
  }

  raw_counts <- J * class_probs
  class_counts <- floor(raw_counts)
  remainder <- J - sum(class_counts)
  if (remainder > 0) {
    fractional_order <- order(raw_counts - class_counts, decreasing = TRUE)
    add_to <- fractional_order[seq_len(remainder)]
    class_counts[add_to] <- class_counts[add_to] + 1
  }

  sample(rep(names(class_counts), times = class_counts), size = J, replace = FALSE)
}

draw_correlated_noise <- function(sd_vec, correlation = NULL) {
  n <- length(sd_vec)
  if (is.null(correlation)) {
    return(rnorm(n, mean = 0, sd = sd_vec))
  }
  cov_mat <- diag(sd_vec, nrow = n) %*% correlation %*% diag(sd_vec, nrow = n)
  drop(t(chol(cov_mat + diag(1e-10, n))) %*% rnorm(n))
}

effect_function <- function(effect_class, x, amplitude = 1) {
  x_scaled <- (x - min(x)) / diff(range(x))
  midpoint <- 0.5

  switch(
    effect_class,
    constant = rep(rnorm(1, 0, 1), length(x)),
    linear = rnorm(1, 0, 1) + amplitude * rnorm(1, 0.5, 0.15) * x_scaled,
    quadratic = rnorm(1, 0, 1) +
      amplitude * (x_scaled - midpoint)^2 * sample(c(-1, 1), 1),
    early = amplitude * exp(-0.5 * ((x_scaled - 0.15) / 0.12)^2),
    middle = amplitude * exp(-0.5 * ((x_scaled - 0.50) / 0.12)^2),
    late = amplitude * exp(-0.5 * ((x_scaled - 0.85) / 0.12)^2),
    switch = amplitude * tanh((x_scaled - midpoint) / 0.12),
    transient = amplitude * exp(-0.5 * ((x_scaled - 0.50) / 0.08)^2),
    abrupt = amplitude * as.numeric(x_scaled >= midpoint),
    stop("Unknown effect class: ", effect_class)
  )
}

simulate_bspline_effect <- function(x,
                                    amplitude = 2,
                                    df = 6,
                                    coefficient_sd = 1) {
  x_scaled <- (x - min(x)) / diff(range(x))
  basis <- splines::bs(x_scaled, df = df, degree = 3, intercept = TRUE)
  coefs <- rnorm(ncol(basis), mean = 0, sd = coefficient_sd)
  f <- drop(basis %*% coefs)
  f <- f - mean(f)
  max_abs <- max(abs(f))
  if (max_abs > 0) {
    f <- amplitude * f / max_abs
  }
  f
}

simulate_bspline_effect_on_grids <- function(time_grid,
                                             evaluation_grid,
                                             amplitude = 2,
                                             df = 6,
                                             coefficient_sd = 1) {
  time_grid <- as.numeric(time_grid)
  evaluation_grid <- as.numeric(evaluation_grid)
  if (length(time_grid) < 4 ||
      any(!is.finite(time_grid)) ||
      any(!is.finite(evaluation_grid)) ||
      min(evaluation_grid) < min(time_grid) ||
      max(evaluation_grid) > max(time_grid) ||
      !is.finite(amplitude) ||
      amplitude <= 0 ||
      !is.finite(df) ||
      df <= 3 ||
      !is.finite(coefficient_sd) ||
      coefficient_sd <= 0) {
    stop("Invalid paired B-spline effect settings.")
  }

  time_scaled <- (time_grid - min(time_grid)) / diff(range(time_grid))
  evaluation_scaled <- (
    evaluation_grid - min(time_grid)
  ) / diff(range(time_grid))
  basis_observed <- splines::bs(
    time_scaled,
    df = df,
    degree = 3,
    intercept = TRUE
  )
  basis_evaluation <- predict(basis_observed, newx = evaluation_scaled)
  coefficients <- stats::rnorm(
    ncol(basis_observed),
    mean = 0,
    sd = coefficient_sd
  )
  raw_observed <- drop(basis_observed %*% coefficients)
  raw_evaluation <- drop(basis_evaluation %*% coefficients)
  observed_mean <- mean(raw_observed)
  deviation_observed <- raw_observed - observed_mean
  deviation_evaluation <- raw_evaluation - observed_mean
  max_abs <- max(abs(deviation_observed))
  if (!is.finite(max_abs) || max_abs <= .Machine$double.eps) {
    stop("The sampled B-spline deviation is numerically constant.")
  }
  scale_factor <- amplitude / max_abs

  list(
    deviation_observed = scale_factor * deviation_observed,
    deviation_evaluation = scale_factor * deviation_evaluation,
    coefficients = coefficients,
    observed_mean = observed_mean,
    scale_factor = scale_factor
  )
}

make_temporal_functionals <- function(smooth_var,
                                      switch_threshold = 0.25,
                                      switch_minimum_duration = 0) {
  smooth_var <- as.numeric(smooth_var)
  if (length(smooth_var) < 2 || any(!is.finite(smooth_var))) {
    stop("smooth_var must contain at least two finite time points.")
  }
  if (!is.finite(switch_threshold) || switch_threshold <= 0) {
    stop("switch_threshold must be positive and finite.")
  }
  if (!is.finite(switch_minimum_duration) || switch_minimum_duration < 0) {
    stop("switch_minimum_duration must be non-negative and finite.")
  }
  if (switch_minimum_duration > 0 && any(diff(smooth_var) <= 0)) {
    stop("smooth_var must be strictly increasing for persistent-switch testing.")
  }

  grid_cell_widths <- if (switch_minimum_duration > 0) {
    midpoints <- (smooth_var[-1] + smooth_var[-length(smooth_var)]) / 2
    boundaries <- c(smooth_var[1], midpoints, smooth_var[length(smooth_var)])
    diff(boundaries)
  } else {
    NULL
  }

  validate_curve <- function(x) {
    x <- as.numeric(x)
    if (length(x) != length(smooth_var) || any(!is.finite(x))) {
      stop("Each curve must be finite and have the same length as smooth_var.")
    }
    x
  }

  longest_contiguous_duration <- function(indicator) {
    runs <- rle(as.logical(indicator))
    if (!any(runs$values)) {
      return(0)
    }
    run_ends <- cumsum(runs$lengths)
    run_starts <- run_ends - runs$lengths + 1L
    positive_runs <- which(runs$values)
    max(vapply(
      positive_runs,
      function(i) sum(grid_cell_widths[run_starts[i]:run_ends[i]]),
      numeric(1)
    ))
  }

  list(
    early = function(x) {
      x <- validate_curve(x)
      max(abs(x[smooth_var <= 3])) - max(abs(x[smooth_var > 3]))
    },
    middle = function(x) {
      x <- validate_curve(x)
      middle <- smooth_var >= 4 & smooth_var <= 11
      outside <- !middle
      max(abs(x[middle])) - max(abs(x[outside]))
    },
    late = function(x) {
      x <- validate_curve(x)
      max(abs(x[smooth_var >= 12])) - max(abs(x[smooth_var < 12]))
    },
    switch = function(x) {
      x <- validate_curve(x)
      if (switch_minimum_duration > 0) {
        positive_duration <- longest_contiguous_duration(x > switch_threshold)
        negative_duration <- longest_contiguous_duration(x < -switch_threshold)
        return(
          min(positive_duration, negative_duration) -
            switch_minimum_duration
        )
      }
      x_positive <- x[x > 0]
      x_negative <- x[x < 0]
      if (length(x_positive) == 0 || length(x_negative) == 0) {
        return(-switch_threshold)
      }
      min(max(abs(x_positive)), max(abs(x_negative))) - switch_threshold
    }
  )
}

evaluate_temporal_functionals <- function(curves,
                                          smooth_var,
                                          switch_threshold = 0.25,
                                          switch_minimum_duration = 0) {
  curves <- as.matrix(curves)
  storage.mode(curves) <- "numeric"
  if (ncol(curves) != length(smooth_var)) {
    stop("curves must have one column for each smooth_var value.")
  }
  functionals <- make_temporal_functionals(
    smooth_var = smooth_var,
    switch_threshold = switch_threshold,
    switch_minimum_duration = switch_minimum_duration
  )
  out <- vapply(
    functionals,
    function(functional) apply(curves, 1, functional),
    numeric(nrow(curves))
  )
  if (is.null(dim(out))) {
    out <- matrix(out, nrow = nrow(curves), dimnames = list(NULL, names(functionals)))
  }
  rownames(out) <- rownames(curves)
  out
}

random_bspline_truth_curve <- function(time_grid,
                                       evaluation_grid,
                                       amplitude = 2,
                                       df = 6,
                                       degree = 3,
                                       coefficient_sd = 1,
                                       coefficients = NULL,
                                       baseline = 0) {
  time_grid <- as.numeric(time_grid)
  evaluation_grid <- as.numeric(evaluation_grid)
  if (length(time_grid) < 4 || any(!is.finite(time_grid)) ||
      any(!is.finite(evaluation_grid)) || diff(range(time_grid)) <= 0 ||
      any(evaluation_grid < min(time_grid)) || any(evaluation_grid > max(time_grid))) {
    stop("time_grid and evaluation_grid must be finite, with evaluation_grid inside time_grid bounds.")
  }
  if (!is.finite(amplitude) || amplitude <= 0 || !is.finite(coefficient_sd) ||
      coefficient_sd <= 0 || !is.finite(baseline)) {
    stop("amplitude, coefficient_sd, and baseline must be finite; amplitude and coefficient_sd must be positive.")
  }

  scale_time <- function(x) (x - min(time_grid)) / diff(range(time_grid))
  basis_observed <- splines::bs(
    scale_time(time_grid),
    df = df,
    degree = degree,
    intercept = TRUE
  )
  if (is.null(coefficients)) {
    coefficients <- rnorm(ncol(basis_observed), mean = 0, sd = coefficient_sd)
  }
  coefficients <- as.numeric(coefficients)
  if (length(coefficients) != ncol(basis_observed) || any(!is.finite(coefficients))) {
    stop("coefficients must be finite and match the random B-spline basis dimension.")
  }

  basis_evaluation <- splines::bs(
    scale_time(evaluation_grid),
    knots = attr(basis_observed, "knots"),
    degree = attr(basis_observed, "degree"),
    Boundary.knots = attr(basis_observed, "Boundary.knots"),
    intercept = TRUE
  )
  raw_observed <- drop(basis_observed %*% coefficients)
  raw_evaluation <- drop(basis_evaluation %*% coefficients)
  centered_mean <- mean(raw_observed)
  centered_observed <- raw_observed - centered_mean
  scale_factor <- max(abs(centered_observed))
  if (!is.finite(scale_factor) || scale_factor <= 0) {
    stop("The random B-spline draw has no dynamic variation.")
  }
  scale_factor <- amplitude / scale_factor
  deviation_observed <- centered_observed * scale_factor
  deviation_evaluation <- (raw_evaluation - centered_mean) * scale_factor

  list(
    beta_observed = baseline + deviation_observed,
    beta_evaluation = baseline + deviation_evaluation,
    deviation_observed = deviation_observed,
    deviation_evaluation = deviation_evaluation,
    coefficients = coefficients,
    baseline = baseline,
    knots = attr(basis_observed, "knots"),
    boundary_knots = attr(basis_observed, "Boundary.knots"),
    degree = attr(basis_observed, "degree")
  )
}

targeted_time_window <- function(time_group, x) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) {
    stop("x must contain finite time values.")
  }
  switch(
    time_group,
    early = x <= 3,
    middle = x >= 4 & x <= 11,
    late = x >= 12,
    stop("Unsupported time group: ", time_group)
  )
}

targeted_local_bspline_parameter_ranges <- function(time_group,
                                                     profile = c("narrow", "broad", "spiky")) {
  profile <- match.arg(profile)
  if (profile == "spiky") {
    profile <- "narrow"
  }
  ranges <- if (profile == "narrow") {
    data.frame(
      time_group = c("early", "middle", "late"),
      center_min = c(1.20, 6.50, 13.20),
      center_max = c(1.80, 8.50, 13.80),
      scale_min = c(0.55, 0.70, 0.55),
      scale_max = c(0.70, 1.00, 0.70),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      time_group = c("early", "middle", "late"),
      center_min = c(1.40, 7.25, 13.40),
      center_max = c(1.60, 7.75, 13.60),
      scale_min = c(1.85, 2.50, 1.85),
      scale_max = c(2.05, 2.80, 2.05),
      stringsAsFactors = FALSE
    )
  }
  row <- ranges[ranges$time_group == time_group, , drop = FALSE]
  if (nrow(row) != 1) {
    stop("Unsupported time group: ", time_group)
  }
  row
}

local_cubic_bspline_lobe <- function(x, center, scale) {
  x <- as.numeric(x)
  if (any(!is.finite(x)) || !is.finite(center) || !is.finite(scale) || scale <= 0) {
    stop("x, center, and scale must be finite, and scale must be positive.")
  }
  u <- abs((x - center) / scale)
  out <- numeric(length(x))
  near_center <- u < 1
  outer_support <- u >= 1 & u < 2
  out[near_center] <- (4 - 6 * u[near_center]^2 + 3 * u[near_center]^3) / 6
  out[outer_support] <- (2 - u[outer_support])^3 / 6
  out / (2 / 3)
}

sample_targeted_local_bspline_lobe <- function(x,
                                               time_group,
                                               profile = c("narrow", "broad", "spiky")) {
  profile <- match.arg(profile)
  ranges <- targeted_local_bspline_parameter_ranges(time_group, profile = profile)
  center <- stats::runif(1, ranges$center_min, ranges$center_max)
  scale <- stats::runif(1, ranges$scale_min, ranges$scale_max)
  list(
    values = local_cubic_bspline_lobe(x, center = center, scale = scale),
    center = center,
    scale = scale,
    time_group = time_group
  )
}

count_effective_sign_transitions <- function(x, tolerance = 1e-8) {
  x <- as.numeric(x)
  if (any(!is.finite(x)) || !is.finite(tolerance) || tolerance < 0) {
    stop("x must be finite and tolerance must be nonnegative and finite.")
  }
  signs <- sign(x)
  signs[abs(x) <= tolerance] <- 0
  nonzero_signs <- signs[signs != 0]
  if (length(nonzero_signs) < 2) {
    return(0L)
  }
  as.integer(sum(diff(nonzero_signs) != 0))
}

summarize_targeted_truth_curve <- function(curve,
                                           evaluation_grid,
                                           time_group,
                                           switch_threshold = 0.25) {
  curve <- as.numeric(curve)
  evaluation_grid <- as.numeric(evaluation_grid)
  if (length(curve) != length(evaluation_grid) || any(!is.finite(curve)) ||
      any(!is.finite(evaluation_grid))) {
    stop("curve and evaluation_grid must have the same finite length.")
  }
  target_mask <- targeted_time_window(time_group, evaluation_grid)
  if (!any(target_mask) || all(target_mask)) {
    stop("The target time window must contain, but not exhaust, evaluation_grid.")
  }
  target_peak <- max(abs(curve[target_mask]))
  outside_peak <- max(abs(curve[!target_mask]))
  positive_peak <- if (any(curve > 0)) max(curve[curve > 0]) else 0
  negative_peak <- if (any(curve < 0)) max(abs(curve[curve < 0])) else 0
  functionals <- evaluate_temporal_functionals(
    matrix(curve, nrow = 1),
    smooth_var = evaluation_grid,
    switch_threshold = switch_threshold
  )
  list(
    target_peak = target_peak,
    outside_peak = outside_peak,
    target_to_outside_ratio = if (outside_peak > 0) target_peak / outside_peak else Inf,
    positive_peak = positive_peak,
    negative_peak = negative_peak,
    effective_sign_transitions = count_effective_sign_transitions(
      curve,
      tolerance = max(abs(curve)) * 0.02
    ),
    contrasts = drop(functionals[1, ])
  )
}

sample_spiky_local_bspline_truth <- function(time_group,
                                             switch_status,
                                             time_grid,
                                             evaluation_grid,
                                             df = 16,
                                             degree = 3,
                                             switch_threshold = 0.25,
                                             minimum_location_margin = 0.60,
                                             minimum_location_ratio = 2,
                                             non_switch_baseline_fraction = 0.75,
                                             target_centered_rms = 0.70,
                                             max_attempts = 1000) {
  allowed_time_groups <- c("early", "middle", "late")
  allowed_switch_status <- c("switch", "non-switch")
  if (!time_group %in% allowed_time_groups ||
      !switch_status %in% allowed_switch_status ||
      df <= degree || degree < 1 || !is.finite(target_centered_rms) ||
      target_centered_rms <= 0 || max_attempts < 1) {
    stop("Invalid spiky local B-spline truth settings.")
  }

  boundary_knots <- range(evaluation_grid)
  evaluation_scaled <- (evaluation_grid - boundary_knots[1]) / diff(boundary_knots)
  observed_scaled <- (time_grid - boundary_knots[1]) / diff(boundary_knots)
  basis_evaluation <- splines::bs(
    evaluation_scaled,
    df = df,
    degree = degree,
    intercept = TRUE,
    Boundary.knots = c(0, 1)
  )
  basis_observed <- splines::bs(
    observed_scaled,
    knots = attr(basis_evaluation, "knots"),
    degree = degree,
    intercept = TRUE,
    Boundary.knots = c(0, 1)
  )
  interior <- seq_len(ncol(basis_evaluation))
  if (length(interior) > 4) {
    interior <- interior[-c(1, 2, length(interior) - 1, length(interior))]
  }
  peak_times <- evaluation_grid[
    apply(basis_evaluation[, interior, drop = FALSE], 2, which.max)
  ]
  eligible_basis <- interior[targeted_time_window(time_group, peak_times)]
  if (length(eligible_basis) == 0) {
    stop("The spiky B-spline basis has no support centered in the target window.")
  }

  for (attempt in seq_len(max_attempts)) {
    selected_basis <- sample(eligible_basis, size = 1)
    direction <- sample(c(-1, 1), size = 1)
    deviation_observed <- direction * basis_observed[, selected_basis]
    deviation_evaluation <- direction * basis_evaluation[, selected_basis]
    centering_constant <- mean(deviation_observed)
    deviation_observed <- deviation_observed - centering_constant
    deviation_evaluation <- deviation_evaluation - centering_constant
    observed_rms <- sqrt(mean(deviation_observed^2))
    if (!is.finite(observed_rms) || observed_rms <= .Machine$double.eps) {
      next
    }
    rms_scale_factor <- target_centered_rms / observed_rms
    deviation_observed <- deviation_observed * rms_scale_factor
    deviation_evaluation <- deviation_evaluation * rms_scale_factor

    if (switch_status == "switch") {
      baseline <- 0
    } else {
      opposite_background <- if (direction > 0) {
        max(0, -min(deviation_evaluation))
      } else {
        max(0, max(deviation_evaluation))
      }
      baseline <- direction * (
        opposite_background +
          switch_threshold +
          non_switch_baseline_fraction * target_centered_rms
      )
    }
    beta_observed <- deviation_observed + baseline
    beta_evaluation <- deviation_evaluation + baseline
    diagnostics <- summarize_targeted_truth_curve(
      curve = beta_evaluation,
      evaluation_grid = evaluation_grid,
      time_group = time_group,
      switch_threshold = switch_threshold
    )
    location_ok <- diagnostics$contrasts[[time_group]] >= minimum_location_margin &&
      diagnostics$target_to_outside_ratio >= minimum_location_ratio
    switch_ok <- if (switch_status == "switch") {
      diagnostics$contrasts[["switch"]] > 0 &&
        diagnostics$effective_sign_transitions >= 1L
    } else {
      diagnostics$contrasts[["switch"]] <= 0 &&
        diagnostics$effective_sign_transitions == 0L
    }
    if (location_ok && switch_ok) {
      return(list(
        beta_observed = beta_observed,
        beta_evaluation = beta_evaluation,
        contrasts = diagnostics$contrasts,
        diagnostics = data.frame(
          time_group = time_group,
          switch_status = switch_status,
          profile = "spiky",
          primary_center = evaluation_grid[which.max(abs(deviation_evaluation))],
          primary_scale = NA_real_,
          secondary_group = NA_character_,
          secondary_center = NA_real_,
          secondary_scale = NA_real_,
          primary_amplitude = max(abs(deviation_evaluation)),
          secondary_amplitude = NA_real_,
          baseline = baseline,
          centered_rms = sqrt(mean((beta_observed - mean(beta_observed))^2)),
          rms_scale_factor = rms_scale_factor,
          selected_basis = selected_basis,
          target_peak = diagnostics$target_peak,
          outside_peak = diagnostics$outside_peak,
          target_to_outside_ratio = diagnostics$target_to_outside_ratio,
          positive_peak = diagnostics$positive_peak,
          negative_peak = diagnostics$negative_peak,
          effective_sign_transitions = diagnostics$effective_sign_transitions,
          stringsAsFactors = FALSE
        ),
        attempt = attempt
      ))
    }
  }
  stop(
    "Could not generate a spiky local B-spline truth for ",
    time_group,
    " / ",
    switch_status,
    " within max_attempts."
  )
}

sample_multispike_local_bspline_truth <- function(time_group,
                                                  switch_status,
                                                  time_grid,
                                                  evaluation_grid,
                                                  spike_count = NULL,
                                                  df = 16,
                                                  degree = 3,
                                                  secondary_fraction = c(0.40, 0.65),
                                                  switch_threshold = 0.25,
                                                  minimum_location_margin = 0.60,
                                                  minimum_location_ratio = 1.4,
                                                  target_centered_rms = 0.90,
                                                  minimum_peak_separation = 3,
                                                  secondary_selection = c("random", "nearest"),
                                                  non_switch_baseline_fraction = 0,
                                                  direction = NULL,
                                                  max_attempts = 1000) {
  secondary_selection <- match.arg(secondary_selection)
  allowed_time_groups <- c("early", "middle", "late")
  allowed_switch_status <- c("switch", "non-switch")
  if (!time_group %in% allowed_time_groups ||
      !switch_status %in% allowed_switch_status ||
      df <= degree || degree < 1 ||
      length(secondary_fraction) != 2 ||
      any(!is.finite(secondary_fraction)) ||
      any(secondary_fraction <= 0) ||
      secondary_fraction[1] > secondary_fraction[2] ||
      secondary_fraction[2] >= 1 ||
      !is.finite(switch_threshold) || switch_threshold <= 0 ||
      !is.finite(minimum_location_margin) || minimum_location_margin <= 0 ||
      !is.finite(minimum_location_ratio) || minimum_location_ratio <= 1 ||
      !is.finite(target_centered_rms) || target_centered_rms <= 0 ||
      !is.finite(minimum_peak_separation) || minimum_peak_separation < 0 ||
      !is.finite(non_switch_baseline_fraction) ||
      non_switch_baseline_fraction < 0 ||
      (!is.null(direction) &&
        (length(direction) != 1 || !is.finite(direction) ||
          !direction %in% c(-1, 1))) ||
      max_attempts < 1) {
    stop("Invalid multi-spike local B-spline truth settings.")
  }
  if (is.null(spike_count)) {
    spike_count <- if (switch_status == "switch") 2L else sample(1:2, size = 1)
  }
  spike_count <- as.integer(spike_count)
  if (!spike_count %in% 1:2 ||
      (switch_status == "switch" && spike_count != 2L)) {
    stop("Switch truths require two spikes; non-switch truths allow one or two spikes.")
  }

  boundary_knots <- range(evaluation_grid)
  if (length(time_grid) < 4 || any(!is.finite(time_grid)) ||
      any(!is.finite(evaluation_grid)) || diff(boundary_knots) <= 0 ||
      any(time_grid < boundary_knots[1]) || any(time_grid > boundary_knots[2])) {
    stop("time_grid and evaluation_grid must be finite and cover the same time range.")
  }
  evaluation_scaled <- (evaluation_grid - boundary_knots[1]) / diff(boundary_knots)
  observed_scaled <- (time_grid - boundary_knots[1]) / diff(boundary_knots)
  basis_evaluation <- splines::bs(
    evaluation_scaled,
    df = df,
    degree = degree,
    intercept = TRUE,
    Boundary.knots = c(0, 1)
  )
  basis_observed <- splines::bs(
    observed_scaled,
    knots = attr(basis_evaluation, "knots"),
    degree = degree,
    intercept = TRUE,
    Boundary.knots = c(0, 1)
  )
  basis_evaluation[basis_evaluation < 0] <- 0
  basis_observed[basis_observed < 0] <- 0

  interior <- seq_len(ncol(basis_evaluation))
  if (length(interior) > 4) {
    interior <- interior[-c(1, 2, length(interior) - 1, length(interior))]
  }
  peak_times <- vapply(
    interior,
    function(index) evaluation_grid[which.max(basis_evaluation[, index])],
    numeric(1)
  )
  names(peak_times) <- as.character(interior)
  target_basis <- interior[targeted_time_window(time_group, peak_times)]
  if (length(target_basis) == 0) {
    stop("The multi-spike B-spline basis has no peak in the target window.")
  }

  for (attempt in seq_len(max_attempts)) {
    primary_basis <- target_basis[sample.int(length(target_basis), size = 1)]
    primary_peak_time <- unname(peak_times[as.character(primary_basis)])
    secondary_basis <- NA_integer_
    secondary_peak_time <- NA_real_
    secondary_amplitude_fraction <- NA_real_

    if (spike_count == 2L) {
      secondary_candidates <- interior[
        !targeted_time_window(time_group, peak_times) &
          abs(peak_times - primary_peak_time) >= minimum_peak_separation
      ]
      if (length(secondary_candidates) == 0) {
        next
      }
      if (secondary_selection == "nearest") {
        candidate_distance <- abs(
          peak_times[as.character(secondary_candidates)] - primary_peak_time
        )
        nearest_candidates <- secondary_candidates[
          candidate_distance == min(candidate_distance)
        ]
        secondary_basis <- nearest_candidates[
          sample.int(length(nearest_candidates), size = 1)
        ]
      } else {
        secondary_basis <- secondary_candidates[
          sample.int(length(secondary_candidates), size = 1)
        ]
      }
      secondary_peak_time <- unname(peak_times[as.character(secondary_basis)])
      secondary_amplitude_fraction <- stats::runif(
        1,
        secondary_fraction[1],
        secondary_fraction[2]
      )
    }

    curve_direction <- if (is.null(direction)) {
      sample(c(-1, 1), size = 1)
    } else {
      direction
    }
    beta_observed <- curve_direction * basis_observed[, primary_basis]
    beta_evaluation <- curve_direction * basis_evaluation[, primary_basis]
    if (spike_count == 2L) {
      secondary_direction <- if (switch_status == "switch") {
        -curve_direction
      } else {
        curve_direction
      }
      beta_observed <- beta_observed +
        secondary_direction * secondary_amplitude_fraction *
          basis_observed[, secondary_basis]
      beta_evaluation <- beta_evaluation +
        secondary_direction * secondary_amplitude_fraction *
          basis_evaluation[, secondary_basis]
    }

    baseline <- if (switch_status == "non-switch") {
      curve_direction * non_switch_baseline_fraction
    } else {
      0
    }
    beta_observed <- beta_observed + baseline
    beta_evaluation <- beta_evaluation + baseline
    observed_centered_rms <- sqrt(mean((beta_observed - mean(beta_observed))^2))
    if (!is.finite(observed_centered_rms) ||
        observed_centered_rms <= .Machine$double.eps) {
      next
    }
    rms_scale_factor <- target_centered_rms / observed_centered_rms
    beta_observed <- beta_observed * rms_scale_factor
    beta_evaluation <- beta_evaluation * rms_scale_factor
    baseline <- baseline * rms_scale_factor
    diagnostics <- summarize_targeted_truth_curve(
      curve = beta_evaluation,
      evaluation_grid = evaluation_grid,
      time_group = time_group,
      switch_threshold = switch_threshold
    )
    location_ok <- diagnostics$contrasts[[time_group]] >= minimum_location_margin &&
      diagnostics$target_to_outside_ratio >= minimum_location_ratio
    switch_ok <- if (switch_status == "switch") {
      diagnostics$contrasts[["switch"]] > 0 &&
        diagnostics$effective_sign_transitions >= 1L
    } else {
      diagnostics$contrasts[["switch"]] <= 0 &&
        diagnostics$effective_sign_transitions == 0L
    }
    if (location_ok && switch_ok) {
      spike_pattern <- if (switch_status == "switch") {
        "opposite-sign double"
      } else if (spike_count == 1L) {
        "single"
      } else {
        "same-sign double"
      }
      return(list(
        beta_observed = beta_observed,
        beta_evaluation = beta_evaluation,
        contrasts = diagnostics$contrasts,
        diagnostics = data.frame(
          time_group = time_group,
          switch_status = switch_status,
          profile = "spiky_multispike_v2",
          spike_count = spike_count,
          spike_pattern = spike_pattern,
          primary_center = primary_peak_time,
          primary_scale = NA_real_,
          secondary_group = if (spike_count == 2L) {
            if (secondary_peak_time <= 3) {
              "early"
            } else if (secondary_peak_time >= 12) {
              "late"
            } else {
              "middle"
            }
          } else {
            NA_character_
          },
          secondary_center = secondary_peak_time,
          secondary_scale = NA_real_,
          primary_amplitude = max(abs(
            rms_scale_factor * basis_evaluation[, primary_basis]
          )),
          secondary_amplitude = if (spike_count == 2L) {
            max(abs(
              rms_scale_factor * secondary_amplitude_fraction *
                basis_evaluation[, secondary_basis]
            ))
          } else {
            NA_real_
          },
          secondary_amplitude_fraction = secondary_amplitude_fraction,
          secondary_selection = secondary_selection,
          baseline = baseline,
          centered_rms = sqrt(mean((beta_observed - mean(beta_observed))^2)),
          rms_scale_factor = rms_scale_factor,
          selected_basis = primary_basis,
          secondary_basis = secondary_basis,
          target_peak = diagnostics$target_peak,
          outside_peak = diagnostics$outside_peak,
          target_to_outside_ratio = diagnostics$target_to_outside_ratio,
          positive_peak = diagnostics$positive_peak,
          negative_peak = diagnostics$negative_peak,
          effective_sign_transitions = diagnostics$effective_sign_transitions,
          stringsAsFactors = FALSE
        ),
        attempt = attempt
      ))
    }
  }
  stop(
    "Could not generate a multi-spike local B-spline truth for ",
    time_group,
    " / ",
    switch_status,
    " / ",
    spike_count,
    " spike(s) within max_attempts."
  )
}

raised_cosine_peak <- function(x, center, half_width) {
  if (length(center) != 1 || !is.finite(center) ||
      length(half_width) != 1 || !is.finite(half_width) ||
      half_width <= 0 || any(!is.finite(x))) {
    stop("Invalid raised-cosine peak settings.")
  }
  distance <- abs(x - center)
  out <- numeric(length(x))
  inside <- distance <= half_width
  out[inside] <- (
    1 + cos(pi * (x[inside] - center) / half_width)
  ) / 2
  out
}

simulate_raised_cosine_multipeak_effect_set <- function(
    n_variants = 2400,
    time_grid = make_time_grid(),
    evaluation_grid = seq(0, 15, by = 0.1),
    class_probs = c(
      dynamic_bspline = 0.20,
      constant = 0.40,
      zero = 0.40
    ),
    width_levels = c(
      not_spiky = 3.00,
      mildly_spiky = 2.25,
      spiky = 1.50,
      very_spiky = 0.90
    ),
    spike_counts = 1:3,
    shape_cell_probs = NULL,
    primary_time_groups = NULL,
    center_by_observed_mean = TRUE,
    validate_functional_labels = TRUE,
    switch_threshold = 0.25,
    relative_amplitude_range = c(0.75, 1.00),
    target_centered_rms = 0.90,
    baseline_sd = 1,
    constant_sd = NULL,
    dynamic_baseline_sd = NULL,
    exact_class_counts = TRUE,
    seed = 12345,
    class_seed = NULL,
    constant_seed = NULL,
    shape_seed = NULL,
    scenario = "internal_raised_cosine_multipeak_factorial") {
  if (is.null(constant_sd)) constant_sd <- baseline_sd
  if (is.null(dynamic_baseline_sd)) dynamic_baseline_sd <- baseline_sd
  required_classes <- c("dynamic_bspline", "constant", "zero")
  if (!identical(sort(names(class_probs)), sort(required_classes)) ||
      abs(sum(class_probs) - 1) > 1e-8 ||
      n_variants < 20 ||
      length(width_levels) < 1 ||
      is.null(names(width_levels)) ||
      any(!nzchar(names(width_levels))) ||
      any(!is.finite(width_levels)) ||
      any(width_levels <= 0) ||
      length(spike_counts) < 1 ||
      any(!is.finite(spike_counts)) ||
      any(as.integer(spike_counts) != spike_counts) ||
      any(!(as.integer(spike_counts) %in% 1:3)) ||
      anyDuplicated(as.integer(spike_counts)) ||
      (!is.null(shape_cell_probs) &&
        (length(shape_cell_probs) < 1 ||
          is.null(names(shape_cell_probs)) ||
          any(!nzchar(names(shape_cell_probs))) ||
          anyDuplicated(names(shape_cell_probs)) ||
          any(!is.finite(shape_cell_probs)) ||
          any(shape_cell_probs <= 0) ||
          abs(sum(shape_cell_probs) - 1) > 1e-8)) ||
      (!is.null(primary_time_groups) &&
        (length(primary_time_groups) < 1 ||
          anyDuplicated(primary_time_groups) ||
          any(!primary_time_groups %in% c("early", "middle", "late")))) ||
      length(center_by_observed_mean) != 1 ||
      is.na(center_by_observed_mean) ||
      !is.logical(center_by_observed_mean) ||
      length(validate_functional_labels) != 1 ||
      is.na(validate_functional_labels) ||
      !is.logical(validate_functional_labels) ||
      !is.finite(switch_threshold) ||
      switch_threshold <= 0 ||
      length(relative_amplitude_range) != 2 ||
      any(!is.finite(relative_amplitude_range)) ||
      any(relative_amplitude_range <= 0) ||
      relative_amplitude_range[1] > relative_amplitude_range[2] ||
      !is.finite(target_centered_rms) ||
      target_centered_rms <= 0 ||
      !is.finite(baseline_sd) ||
      baseline_sd < 0 ||
      !is.finite(constant_sd) ||
      constant_sd < 0 ||
      !is.finite(dynamic_baseline_sd) ||
      dynamic_baseline_sd < 0 ||
      any(vapply(
        list(class_seed, constant_seed, shape_seed),
        function(x) !is.null(x) && (length(x) != 1 || !is.finite(x)),
        logical(1)
      )) ||
      any(!is.finite(time_grid)) ||
      any(!is.finite(evaluation_grid)) ||
      min(evaluation_grid) > min(time_grid) ||
      max(evaluation_grid) < max(time_grid) ||
      !nzchar(scenario)) {
    stop("Invalid raised-cosine multi-peak effect settings.")
  }
  if (!is.null(primary_time_groups) &&
      max(as.integer(spike_counts)) > length(primary_time_groups)) {
    stop("Timed multi-peak effects require at least one time group per peak.")
  }

  split_seed_streams <- any(vapply(
    list(class_seed, constant_seed, shape_seed),
    Negate(is.null),
    logical(1)
  ))
  set.seed(if (is.null(class_seed)) seed else class_seed)
  effect_class <- sample_effect_classes(
    J = n_variants,
    class_probs = class_probs,
    exact_class_counts = exact_class_counts
  )
  dynamic_index <- which(effect_class == "dynamic_bspline")
  constant_index <- which(effect_class == "constant")
  if (split_seed_streams) {
    set.seed(if (is.null(shape_seed)) seed else shape_seed)
  }

  spike_counts <- sort(as.integer(spike_counts))
  factorial_cells <- do.call(rbind, lapply(names(width_levels), function(label) {
    cells <- list()
    if (1L %in% spike_counts) {
      cells[[length(cells) + 1L]] <- data.frame(
        width_label = label,
        width_half = unname(width_levels[label]),
        spike_count = 1L,
        sign_pattern = "single",
        stringsAsFactors = FALSE
      )
    }
    multi_counts <- spike_counts[spike_counts > 1L]
    if (length(multi_counts) > 0) {
      multi_cells <- expand.grid(
        spike_count = multi_counts,
        sign_pattern = c("same-sign", "alternating-sign"),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
      )
      multi_cells$width_label <- label
      multi_cells$width_half <- unname(width_levels[label])
      multi_cells <- multi_cells[
        ,
        c("width_label", "width_half", "spike_count", "sign_pattern")
      ]
      cells[[length(cells) + 1L]] <- multi_cells
    }
    do.call(rbind, cells)
  }))
  factorial_cells$cell_id <- paste(
    paste0("k", factorial_cells$spike_count),
    factorial_cells$width_label,
    factorial_cells$sign_pattern,
    sep = "__"
  )
  if (is.null(shape_cell_probs)) {
    cell_counts <- exact_balanced_counts(
      length(dynamic_index),
      factorial_cells$cell_id
    )
  } else {
    if (!setequal(names(shape_cell_probs), factorial_cells$cell_id)) {
      stop(
        "shape_cell_probs names must match the generated shape-cell IDs: ",
        paste(factorial_cells$cell_id, collapse = ", ")
      )
    }
    cell_counts <- exact_proportional_counts(
      length(dynamic_index),
      shape_cell_probs[factorial_cells$cell_id]
    )
  }
  assignment <- do.call(rbind, lapply(
    factorial_cells$cell_id,
    function(cell_id) {
      cell <- factorial_cells[
        factorial_cells$cell_id == cell_id,
        ,
        drop = FALSE
      ]
      count <- unname(cell_counts[cell_id])
      cell_assignment <- data.frame(
        cell[rep(1, count), , drop = FALSE],
        latent_id = seq_len(count),
        stringsAsFactors = FALSE
      )
      if (!is.null(primary_time_groups)) {
        cell_assignment$time_group <- rep(
          primary_time_groups,
          length.out = count
        )
      }
      cell_assignment
    }
  ))
  assignment <- assignment[
    sample.int(nrow(assignment), size = nrow(assignment)),
    ,
    drop = FALSE
  ]
  rownames(assignment) <- NULL

  if (is.null(primary_time_groups)) {
    latent_groups <- unique(assignment[, c("spike_count", "sign_pattern")])
    latent_specs <- do.call(rbind, lapply(seq_len(nrow(latent_groups)), function(i) {
      spike_count <- latent_groups$spike_count[i]
      sign_pattern <- latent_groups$sign_pattern[i]
      n_latent <- max(
        assignment$latent_id[
          assignment$spike_count == spike_count &
            assignment$sign_pattern == sign_pattern
        ]
      )
      do.call(rbind, lapply(seq_len(n_latent), function(latent_id) {
        center_shift <- if (spike_count == 3L) {
          stats::runif(1, -0.25, 0.25)
        } else {
          stats::runif(1, -0.35, 0.35)
        }
        centers <- switch(
          as.character(spike_count),
          "1" = c(3, 7.5, 12)[(latent_id - 1L) %% 3L + 1L] + center_shift,
          "2" = c(4, 11) + center_shift,
          "3" = c(3, 7.5, 12) + center_shift
        )
        direction <- sample(c(-1, 1), size = 1)
        signs <- if (sign_pattern == "alternating-sign") {
          direction * (-1)^(seq_len(spike_count) - 1L)
        } else {
          rep(direction, spike_count)
        }
        relative_amplitudes <- c(
          1,
          if (spike_count > 1L) {
            stats::runif(
              spike_count - 1L,
              relative_amplitude_range[1],
              relative_amplitude_range[2]
            )
          } else {
            numeric(0)
          }
        )
        data.frame(
          spike_count = spike_count,
          sign_pattern = sign_pattern,
          latent_id = latent_id,
          centers = I(list(centers)),
          signs = I(list(signs)),
          relative_amplitudes = I(list(relative_amplitudes)),
          baseline = stats::rnorm(1, mean = 0, sd = dynamic_baseline_sd),
          stringsAsFactors = FALSE
        )
      }))
    }))
    latent_specs$time_group <- NA_character_
    latent_specs$secondary_group <- NA_character_
  } else {
    timed_assignments <- unique(assignment[
      ,
      c(
        "cell_id", "spike_count", "sign_pattern",
        "latent_id", "time_group"
      )
    ])
    time_group_center_ranges <- list(
      early = c(1.50, 3.50),
      middle = c(5.50, 9.50),
      late = c(11.50, 13.50)
    )
    latent_specs <- do.call(rbind, lapply(seq_len(nrow(timed_assignments)), function(i) {
      specification <- timed_assignments[i, , drop = FALSE]
      spike_count <- specification$spike_count
      sign_pattern <- specification$sign_pattern
      time_group <- specification$time_group
      primary_range <- time_group_center_ranges[[time_group]]
      primary_center <- stats::runif(
        1,
        min = primary_range[1],
        max = primary_range[2]
      )
      secondary_group <- NA_character_
      centers <- primary_center
      if (spike_count > 1L) {
        secondary_groups <- sample(
          setdiff(primary_time_groups, time_group),
          size = spike_count - 1L,
          replace = FALSE
        )
        secondary_centers <- vapply(
          secondary_groups,
          function(group) {
            center_range <- time_group_center_ranges[[group]]
            stats::runif(
              1,
              min = center_range[1],
              max = center_range[2]
            )
          },
          numeric(1)
        )
        secondary_group <- paste(secondary_groups, collapse = ",")
        centers <- c(primary_center, secondary_centers)
      }
      direction <- sample(c(-1, 1), size = 1)
      signs <- if (sign_pattern == "alternating-sign") {
        direction * (-1)^(seq_len(spike_count) - 1L)
      } else {
        rep(direction, spike_count)
      }
      relative_amplitudes <- c(
        1,
        if (spike_count > 1L) {
          stats::runif(
            spike_count - 1L,
            relative_amplitude_range[1],
            relative_amplitude_range[2]
          )
        } else {
          numeric(0)
        }
      )
      data.frame(
        cell_id = specification$cell_id,
        spike_count = spike_count,
        sign_pattern = sign_pattern,
        latent_id = specification$latent_id,
        time_group = time_group,
        secondary_group = secondary_group,
        centers = I(list(centers)),
        signs = I(list(signs)),
        relative_amplitudes = I(list(relative_amplitudes)),
        baseline = stats::rnorm(1, mean = 0, sd = dynamic_baseline_sd),
        stringsAsFactors = FALSE
      )
    }))
  }

  beta_matrix <- matrix(
    0,
    nrow = n_variants,
    ncol = length(time_grid),
    dimnames = list(
      sprintf("variant_%04d", seq_len(n_variants)),
      sprintf("time_%02d", seq_along(time_grid))
    )
  )
  beta_evaluation <- matrix(
    0,
    nrow = n_variants,
    ncol = length(evaluation_grid),
    dimnames = list(
      rownames(beta_matrix),
      sprintf("evaluation_%03d", seq_along(evaluation_grid))
    )
  )

  set.seed(if (is.null(constant_seed)) seed + 1L else constant_seed)
  if (length(constant_index) > 0) {
    constant_effects <- stats::rnorm(
      length(constant_index),
      mean = 0,
      sd = constant_sd
    )
    beta_matrix[constant_index, ] <- constant_effects
    beta_evaluation[constant_index, ] <- constant_effects
  }

  width_label <- rep(NA_character_, n_variants)
  width_half <- rep(NA_real_, n_variants)
  spike_count <- rep(NA_integer_, n_variants)
  sign_pattern <- rep(NA_character_, n_variants)
  cell_id <- rep("dynamic-null", n_variants)
  latent_id <- rep(NA_integer_, n_variants)
  time_group <- rep(NA_character_, n_variants)
  secondary_group <- rep(NA_character_, n_variants)
  switch_status <- rep(NA_character_, n_variants)
  baseline <- rep(NA_real_, n_variants)
  centered_rms <- rep(NA_real_, n_variants)
  peak_centers <- vector("list", n_variants)
  peak_signs <- vector("list", n_variants)
  peak_relative_amplitudes <- vector("list", n_variants)

  for (position in seq_along(dynamic_index)) {
    index <- dynamic_index[position]
    cell <- assignment[position, , drop = FALSE]
    latent <- if (is.null(primary_time_groups)) {
      latent_specs[
        latent_specs$spike_count == cell$spike_count &
          latent_specs$sign_pattern == cell$sign_pattern &
          latent_specs$latent_id == cell$latent_id,
        ,
        drop = FALSE
      ]
    } else {
      latent_specs[
        latent_specs$cell_id == cell$cell_id &
          latent_specs$latent_id == cell$latent_id,
        ,
        drop = FALSE
      ]
    }
    centers <- latent$centers[[1]]
    signs <- latent$signs[[1]]
    relative_amplitudes <- latent$relative_amplitudes[[1]]
    deviation_observed <- numeric(length(time_grid))
    deviation_evaluation <- numeric(length(evaluation_grid))
    for (peak_index in seq_along(centers)) {
      weight <- signs[peak_index] * relative_amplitudes[peak_index]
      deviation_observed <- deviation_observed + weight * raised_cosine_peak(
        time_grid,
        center = centers[peak_index],
        half_width = cell$width_half
      )
      deviation_evaluation <- deviation_evaluation + weight * raised_cosine_peak(
        evaluation_grid,
        center = centers[peak_index],
        half_width = cell$width_half
      )
    }
    observed_mean <- mean(deviation_observed)
    if (center_by_observed_mean) {
      deviation_observed <- deviation_observed - observed_mean
      deviation_evaluation <- deviation_evaluation - observed_mean
      observed_rms <- sqrt(mean(deviation_observed^2))
    } else {
      observed_rms <- sqrt(mean(
        (deviation_observed - mean(deviation_observed))^2
      ))
    }
    if (!is.finite(observed_rms) ||
        observed_rms <= .Machine$double.eps) {
      stop("A raised-cosine multi-peak curve has zero observed RMS.")
    }
    scale_factor <- target_centered_rms / observed_rms
    beta_matrix[index, ] <-
      latent$baseline + scale_factor * deviation_observed
    beta_evaluation[index, ] <-
      latent$baseline + scale_factor * deviation_evaluation

    width_label[index] <- cell$width_label
    width_half[index] <- cell$width_half
    spike_count[index] <- cell$spike_count
    sign_pattern[index] <- cell$sign_pattern
    cell_id[index] <- cell$cell_id
    latent_id[index] <- cell$latent_id
    time_group[index] <- latent$time_group
    secondary_group[index] <- latent$secondary_group
    switch_status[index] <- if (cell$sign_pattern == "alternating-sign") {
      "switch"
    } else {
      "non-switch"
    }
    baseline[index] <- latent$baseline
    centered_rms[index] <- sqrt(mean(
      (beta_matrix[index, ] - mean(beta_matrix[index, ]))^2
    ))
    peak_centers[[index]] <- centers
    peak_signs[[index]] <- signs
    peak_relative_amplitudes[[index]] <- relative_amplitudes
  }

  unit_info <- data.frame(
    unit_index = seq_len(n_variants),
    unit_id = sprintf("%s_%04d", effect_class, seq_len(n_variants)),
    variant_id = rownames(beta_matrix),
    effect_class = effect_class,
    width_label = width_label,
    width_half = width_half,
    spike_count = spike_count,
    sign_pattern = sign_pattern,
    cell_id = cell_id,
    latent_id = latent_id,
    time_group = time_group,
    secondary_group = secondary_group,
    switch_status = switch_status,
    baseline = baseline,
    genetic_main_effect = baseline,
    centered_rms = centered_rms,
    scenario = scenario,
    stringsAsFactors = FALSE
  )
  unit_info$peak_centers <- I(peak_centers)
  unit_info$peak_signs <- I(peak_signs)
  unit_info$peak_relative_amplitudes <- I(peak_relative_amplitudes)

  true_functionals <- evaluate_temporal_functionals(
    beta_evaluation,
    smooth_var = evaluation_grid,
    switch_threshold = switch_threshold
  )
  if (!is.null(primary_time_groups) && length(dynamic_index) > 0) {
    target_index <- match(
      unit_info$time_group[dynamic_index],
      colnames(true_functionals)
    )
    target_values <- true_functionals[
      cbind(dynamic_index, target_index)
    ]
    switch_values <- true_functionals[dynamic_index, "switch"]
    expected_switch <-
      unit_info$sign_pattern[dynamic_index] == "alternating-sign"
    observed_switch <- switch_values > 0
    unit_info$switch_status[dynamic_index] <- ifelse(
      observed_switch,
      "switch",
      "non-switch"
    )
    switch_pattern_invalid <- validate_functional_labels &&
      !center_by_observed_mean &&
      any(observed_switch != expected_switch)
    location_invalid <- validate_functional_labels &&
      any(target_values <= 0)
    if (location_invalid || switch_pattern_invalid) {
      stop("The timed raised-cosine truth failed functional-label validation.")
    }
  }

  list(
    beta_matrix = beta_matrix,
    beta_evaluation = beta_evaluation,
    evaluation_grid = evaluation_grid,
    unit_info = unit_info,
    true_functionals = true_functionals,
    group_counts = cell_counts,
    settings = list(
      truth_mechanism = "raised_cosine_multipeak_factorial",
      width_levels = width_levels,
      spike_counts = spike_counts,
      shape_cell_probs = shape_cell_probs,
      primary_time_groups = primary_time_groups,
      center_by_observed_mean = center_by_observed_mean,
      validate_functional_labels = validate_functional_labels,
      switch_threshold = switch_threshold,
      relative_amplitude_range = relative_amplitude_range,
      target_centered_rms = target_centered_rms,
      baseline_sd = baseline_sd,
      constant_sd = constant_sd,
      dynamic_baseline_sd = dynamic_baseline_sd,
      exact_class_counts = exact_class_counts,
      seed = seed,
      class_seed = class_seed,
      constant_seed = constant_seed,
      shape_seed = shape_seed
    )
  )
}

sample_constrained_broad_bspline_truth <- function(time_group,
                                                   switch_status,
                                                   time_grid,
                                                   evaluation_grid,
                                                   df = 6,
                                                   degree = 3,
                                                   coefficient_sd = 1,
                                                   switch_threshold = 0.25,
                                                   minimum_location_margin = 0.60,
                                                   minimum_location_ratio = 1.4,
                                                   non_switch_baseline_fraction = 0.75,
                                                   target_centered_rms = 0.90,
                                                   max_attempts = 10000) {
  allowed_time_groups <- c("early", "middle", "late")
  allowed_switch_status <- c("switch", "non-switch")
  if (!time_group %in% allowed_time_groups ||
      !switch_status %in% allowed_switch_status ||
      df <= degree || degree < 1 || coefficient_sd <= 0 ||
      !is.finite(target_centered_rms) || target_centered_rms <= 0 ||
      max_attempts < 1) {
    stop("Invalid constrained broad B-spline truth settings.")
  }

  for (attempt in seq_len(max_attempts)) {
    curve <- random_bspline_truth_curve(
      time_grid = time_grid,
      evaluation_grid = evaluation_grid,
      amplitude = 1,
      df = df,
      degree = degree,
      coefficient_sd = coefficient_sd,
      baseline = 0
    )
    observed_rms <- sqrt(mean(curve$deviation_observed^2))
    if (!is.finite(observed_rms) || observed_rms <= .Machine$double.eps) {
      next
    }
    rms_scale_factor <- target_centered_rms / observed_rms
    deviation_observed <- curve$deviation_observed * rms_scale_factor
    deviation_evaluation <- curve$deviation_evaluation * rms_scale_factor

    target_mask <- targeted_time_window(time_group, evaluation_grid)
    target_peak_index <- which(target_mask)[
      which.max(abs(deviation_evaluation[target_mask]))
    ]
    target_direction <- sign(deviation_evaluation[target_peak_index])
    if (target_direction == 0) {
      next
    }
    if (switch_status == "switch") {
      baseline <- 0
    } else {
      opposite_background <- if (target_direction > 0) {
        max(0, -min(deviation_evaluation))
      } else {
        max(0, max(deviation_evaluation))
      }
      baseline <- target_direction * (
        opposite_background +
          switch_threshold +
          non_switch_baseline_fraction * target_centered_rms
      )
    }
    beta_observed <- deviation_observed + baseline
    beta_evaluation <- deviation_evaluation + baseline
    diagnostics <- summarize_targeted_truth_curve(
      curve = beta_evaluation,
      evaluation_grid = evaluation_grid,
      time_group = time_group,
      switch_threshold = switch_threshold
    )
    location_ok <- diagnostics$contrasts[[time_group]] >= minimum_location_margin &&
      diagnostics$target_to_outside_ratio >= minimum_location_ratio
    switch_ok <- if (switch_status == "switch") {
      diagnostics$contrasts[["switch"]] > 0 &&
        diagnostics$effective_sign_transitions >= 1L
    } else {
      diagnostics$contrasts[["switch"]] <= 0 &&
        diagnostics$effective_sign_transitions == 0L
    }
    if (location_ok && switch_ok) {
      return(list(
        beta_observed = beta_observed,
        beta_evaluation = beta_evaluation,
        contrasts = diagnostics$contrasts,
        diagnostics = data.frame(
          time_group = time_group,
          switch_status = switch_status,
          profile = "broad",
          primary_center = evaluation_grid[target_peak_index],
          primary_scale = NA_real_,
          secondary_group = NA_character_,
          secondary_center = NA_real_,
          secondary_scale = NA_real_,
          primary_amplitude = max(abs(deviation_evaluation)),
          secondary_amplitude = NA_real_,
          baseline = baseline,
          centered_rms = sqrt(mean((beta_observed - mean(beta_observed))^2)),
          rms_scale_factor = rms_scale_factor,
          selected_basis = NA_integer_,
          target_peak = diagnostics$target_peak,
          outside_peak = diagnostics$outside_peak,
          target_to_outside_ratio = diagnostics$target_to_outside_ratio,
          positive_peak = diagnostics$positive_peak,
          negative_peak = diagnostics$negative_peak,
          effective_sign_transitions = diagnostics$effective_sign_transitions,
          stringsAsFactors = FALSE
        ),
        attempt = attempt
      ))
    }
  }
  stop(
    "Could not generate a constrained broad B-spline truth for ",
    time_group,
    " / ",
    switch_status,
    " within max_attempts."
  )
}

sample_targeted_local_bspline_truth <- function(time_group,
                                                switch_status,
                                                time_grid,
                                                evaluation_grid,
                                                amplitude = 2,
                                                switch_threshold = 0.25,
                                                minimum_location_margin = 0.60,
                                                minimum_location_ratio = 2,
                                                non_switch_baseline_fraction = 0.15,
                                                non_switch_background_fraction = 0.05,
                                                switch_secondary_fraction = c(0.38, 0.50),
                                                profile = c(
                                                  "narrow",
                                                  "broad",
                                                  "spiky",
                                                  "random_broad"
                                                ),
                                                target_centered_rms = NULL,
                                                max_attempts = 1000) {
  profile <- match.arg(profile)
  allowed_time_groups <- c("early", "middle", "late")
  allowed_switch_status <- c("switch", "non-switch")
  if (!time_group %in% allowed_time_groups || !switch_status %in% allowed_switch_status) {
    stop("Unsupported targeted local B-spline truth group.")
  }
  if (length(time_grid) < 4 || any(!is.finite(time_grid)) ||
      any(!is.finite(evaluation_grid)) || any(evaluation_grid < min(time_grid)) ||
      any(evaluation_grid > max(time_grid)) || !is.finite(amplitude) || amplitude <= 0 ||
      !is.finite(switch_threshold) || switch_threshold <= 0 ||
      !is.finite(minimum_location_margin) || minimum_location_margin <= 0 ||
      !is.finite(minimum_location_ratio) || minimum_location_ratio <= 1 ||
      !is.finite(non_switch_baseline_fraction) || non_switch_baseline_fraction <= 0 ||
      !is.finite(non_switch_background_fraction) || non_switch_background_fraction < 0 ||
      length(switch_secondary_fraction) != 2 || any(!is.finite(switch_secondary_fraction)) ||
      any(switch_secondary_fraction <= 0) ||
      switch_secondary_fraction[1] > switch_secondary_fraction[2] ||
      (!is.null(target_centered_rms) &&
        (length(target_centered_rms) != 1 || !is.finite(target_centered_rms) ||
          target_centered_rms <= 0)) ||
      max_attempts < 1) {
    stop("Invalid targeted local B-spline truth settings.")
  }
  if (profile == "spiky") {
    if (is.null(target_centered_rms)) {
      stop("Spiky targeted local B-spline truths require target_centered_rms.")
    }
    return(sample_spiky_local_bspline_truth(
      time_group = time_group,
      switch_status = switch_status,
      time_grid = time_grid,
      evaluation_grid = evaluation_grid,
      switch_threshold = switch_threshold,
      minimum_location_margin = minimum_location_margin,
      minimum_location_ratio = minimum_location_ratio,
      non_switch_baseline_fraction = non_switch_baseline_fraction,
      target_centered_rms = target_centered_rms,
      max_attempts = max_attempts
    ))
  }
  if (profile == "random_broad") {
    if (is.null(target_centered_rms)) {
      stop("Constrained broad B-spline truths require target_centered_rms.")
    }
    return(sample_constrained_broad_bspline_truth(
      time_group = time_group,
      switch_status = switch_status,
      time_grid = time_grid,
      evaluation_grid = evaluation_grid,
      switch_threshold = switch_threshold,
      minimum_location_margin = minimum_location_margin,
      minimum_location_ratio = minimum_location_ratio,
      non_switch_baseline_fraction = non_switch_baseline_fraction,
      target_centered_rms = target_centered_rms,
      max_attempts = max_attempts * 10
    ))
  }

  other_time_groups <- setdiff(allowed_time_groups, time_group)
  for (attempt in seq_len(max_attempts)) {
    primary_observed <- sample_targeted_local_bspline_lobe(
      time_grid,
      time_group,
      profile = profile
    )
    primary_evaluation <- local_cubic_bspline_lobe(
      evaluation_grid,
      center = primary_observed$center,
      scale = primary_observed$scale
    )
    primary_amplitude <- amplitude * stats::runif(1, 0.95, 1.05)
    direction <- sample(c(-1, 1), size = 1)
    secondary_group <- sample(other_time_groups, size = 1)
    secondary_observed <- sample_targeted_local_bspline_lobe(
      time_grid,
      secondary_group,
      profile = profile
    )
    secondary_evaluation <- local_cubic_bspline_lobe(
      evaluation_grid,
      center = secondary_observed$center,
      scale = secondary_observed$scale
    )

    if (switch_status == "switch") {
      secondary_amplitude <- primary_amplitude * stats::runif(
        1,
        switch_secondary_fraction[1],
        switch_secondary_fraction[2]
      )
      beta_observed <- direction * (
        primary_amplitude * primary_observed$values -
          secondary_amplitude * secondary_observed$values
      )
      beta_evaluation <- direction * (
        primary_amplitude * primary_evaluation -
          secondary_amplitude * secondary_evaluation
      )
      baseline <- 0
      background_amplitude <- 0
    } else {
      baseline <- primary_amplitude * non_switch_baseline_fraction
      background_amplitude <- primary_amplitude * stats::runif(
        1,
        0,
        non_switch_background_fraction
      )
      beta_observed <- direction * (
        baseline +
          primary_amplitude * primary_observed$values +
          background_amplitude * secondary_observed$values
      )
      beta_evaluation <- direction * (
        baseline +
          primary_amplitude * primary_evaluation +
          background_amplitude * secondary_evaluation
      )
      secondary_amplitude <- background_amplitude
    }

    centered_rms <- sqrt(mean((beta_observed - mean(beta_observed))^2))
    rms_scale_factor <- 1
    if (!is.null(target_centered_rms)) {
      if (!is.finite(centered_rms) || centered_rms <= .Machine$double.eps) {
        next
      }
      rms_scale_factor <- target_centered_rms / centered_rms
      beta_observed <- beta_observed * rms_scale_factor
      beta_evaluation <- beta_evaluation * rms_scale_factor
      primary_amplitude <- primary_amplitude * rms_scale_factor
      secondary_amplitude <- secondary_amplitude * rms_scale_factor
      baseline <- baseline * rms_scale_factor
      background_amplitude <- background_amplitude * rms_scale_factor
      centered_rms <- sqrt(mean((beta_observed - mean(beta_observed))^2))
    }

    diagnostics <- summarize_targeted_truth_curve(
      curve = beta_evaluation,
      evaluation_grid = evaluation_grid,
      time_group = time_group,
      switch_threshold = switch_threshold
    )
    location_ok <- diagnostics$contrasts[[time_group]] >= minimum_location_margin &&
      diagnostics$target_to_outside_ratio >= minimum_location_ratio
    switch_ok <- if (switch_status == "switch") {
      diagnostics$contrasts[["switch"]] > 0 &&
        diagnostics$effective_sign_transitions == 1L
    } else {
      diagnostics$contrasts[["switch"]] <= 0 &&
        diagnostics$effective_sign_transitions == 0L
    }
    if (location_ok && switch_ok) {
      return(list(
        beta_observed = beta_observed,
        beta_evaluation = beta_evaluation,
        contrasts = diagnostics$contrasts,
        diagnostics = data.frame(
          time_group = time_group,
          switch_status = switch_status,
          profile = profile,
          primary_center = primary_observed$center,
          primary_scale = primary_observed$scale,
          secondary_group = secondary_group,
          secondary_center = secondary_observed$center,
          secondary_scale = secondary_observed$scale,
          primary_amplitude = primary_amplitude,
          secondary_amplitude = secondary_amplitude,
          baseline = baseline,
          centered_rms = centered_rms,
          rms_scale_factor = rms_scale_factor,
          target_peak = diagnostics$target_peak,
          outside_peak = diagnostics$outside_peak,
          target_to_outside_ratio = diagnostics$target_to_outside_ratio,
          positive_peak = diagnostics$positive_peak,
          negative_peak = diagnostics$negative_peak,
          effective_sign_transitions = diagnostics$effective_sign_transitions,
          stringsAsFactors = FALSE
        ),
        attempt = attempt
      ))
    }
  }
  stop(
    "Could not generate a targeted local B-spline truth for ",
    time_group,
    " / ",
    switch_status,
    " within max_attempts."
  )
}

exact_balanced_counts <- function(n, labels) {
  if (length(labels) == 0 || anyDuplicated(labels) || n < length(labels)) {
    stop("n must be at least the number of unique labels.")
  }
  base_count <- rep(floor(n / length(labels)), length(labels))
  remainder <- n - sum(base_count)
  if (remainder > 0) {
    base_count[seq_len(remainder)] <- base_count[seq_len(remainder)] + 1L
  }
  stats::setNames(base_count, labels)
}

exact_proportional_counts <- function(n, probabilities) {
  if (length(probabilities) == 0 ||
      is.null(names(probabilities)) ||
      any(!nzchar(names(probabilities))) ||
      anyDuplicated(names(probabilities)) ||
      any(!is.finite(probabilities)) ||
      any(probabilities < 0) ||
      abs(sum(probabilities) - 1) > 1e-8 ||
      n < sum(probabilities > 0)) {
    stop("Invalid named probabilities for exact proportional counts.")
  }
  raw_counts <- n * probabilities
  counts <- floor(raw_counts)
  remainder <- n - sum(counts)
  if (remainder > 0) {
    fractional_order <- order(raw_counts - counts, decreasing = TRUE)
    counts[fractional_order[seq_len(remainder)]] <-
      counts[fractional_order[seq_len(remainder)]] + 1L
  }
  stats::setNames(as.integer(counts), names(probabilities))
}

assess_matched_functional_truth <- function(beta_evaluation,
                                            time_group,
                                            switch_status,
                                            evaluation_grid,
                                            switch_threshold = 0.25,
                                            location_truth_margin = 0.10,
                                            switch_truth_margin = 0.10,
                                            non_switch_min_abs = 0.10,
                                            non_switch_min_range_fraction = 0) {
  allowed_time_groups <- c("early", "middle", "late")
  allowed_switch_status <- c("switch", "non-switch")
  beta_evaluation <- as.numeric(beta_evaluation)
  if (!time_group %in% allowed_time_groups ||
      !switch_status %in% allowed_switch_status ||
      length(beta_evaluation) != length(evaluation_grid) ||
      any(!is.finite(beta_evaluation)) ||
      any(!is.finite(evaluation_grid)) ||
      !is.finite(switch_threshold) ||
      switch_threshold <= 0 ||
      !is.finite(location_truth_margin) ||
      location_truth_margin <= 0 ||
      !is.finite(switch_truth_margin) ||
      switch_truth_margin <= 0 ||
      !is.finite(non_switch_min_abs) ||
      non_switch_min_abs <= 0 ||
      !is.finite(non_switch_min_range_fraction) ||
      non_switch_min_range_fraction < 0 ||
      non_switch_min_range_fraction >= 1) {
    stop("Invalid matched functional-truth assessment settings.")
  }

  contrasts <- evaluate_temporal_functionals(
    matrix(beta_evaluation, nrow = 1),
    smooth_var = evaluation_grid,
    switch_threshold = switch_threshold
  )[1, ]
  other_time_groups <- setdiff(allowed_time_groups, time_group)
  timing_valid <- contrasts[[time_group]] >= location_truth_margin &&
    all(contrasts[other_time_groups] <= -location_truth_margin)
  positive_peak <- max(beta_evaluation)
  negative_peak <- min(beta_evaluation)
  same_sign <- all(beta_evaluation > 0) || all(beta_evaluation < 0)
  minimum_absolute_effect <- min(abs(beta_evaluation))
  effect_range <- diff(range(beta_evaluation))
  relative_clearance <- if (effect_range > .Machine$double.eps) {
    minimum_absolute_effect / effect_range
  } else {
    Inf
  }
  required_non_switch_clearance <- max(
    non_switch_min_abs,
    non_switch_min_range_fraction * effect_range
  )
  switch_valid <- if (switch_status == "switch") {
    contrasts[["switch"]] >= switch_truth_margin
  } else {
    same_sign &&
      minimum_absolute_effect >= required_non_switch_clearance
  }

  list(
    valid = timing_valid && switch_valid,
    timing_valid = timing_valid,
    switch_valid = switch_valid,
    contrasts = contrasts,
    positive_peak = positive_peak,
    negative_peak = negative_peak,
    same_sign = same_sign,
    minimum_absolute_effect = minimum_absolute_effect,
    effect_range = effect_range,
    relative_clearance = relative_clearance,
    required_non_switch_clearance = required_non_switch_clearance
  )
}

simulate_matched_functional_effect_set <- function(
    n_variants = 1000,
    truth_mechanism = c("random_bspline", "raised_cosine"),
    time_grid = make_time_grid(),
    evaluation_grid = seq(0, 15, by = 0.1),
    class_probs = c(
      dynamic_bspline = 0.20,
      constant = 0.40,
      zero = 0.40
    ),
    dynamic_main_effect_sd = 1,
    constant_sd = 1,
    bspline_amplitude = 2,
    bspline_df = 6,
    bspline_coefficient_sd = 1,
    cosine_width_half = 1.5,
    cosine_spike_counts = 1:3,
    cosine_relative_amplitude_range = c(0.35, 0.75),
    cosine_target_centered_rms = 0.90,
    switch_threshold = 0.25,
    location_truth_margin = 0.10,
    switch_truth_margin = 0.10,
    non_switch_min_abs = 0.10,
    non_switch_min_range_fraction = 0,
    exact_class_counts = TRUE,
    seed = 12345,
    class_seed = NULL,
    constant_seed = NULL,
    shape_seed = NULL,
    max_attempts = 20000,
    scenario = NULL) {
  truth_mechanism <- match.arg(truth_mechanism)
  required_classes <- c("dynamic_bspline", "constant", "zero")
  cosine_spike_counts <- sort(as.integer(cosine_spike_counts))
  if (!identical(sort(names(class_probs)), sort(required_classes)) ||
      abs(sum(class_probs) - 1) > 1e-8 ||
      n_variants < 30 ||
      !is.finite(dynamic_main_effect_sd) ||
      dynamic_main_effect_sd <= 0 ||
      !is.finite(constant_sd) ||
      constant_sd < 0 ||
      !is.finite(bspline_amplitude) ||
      bspline_amplitude <= 0 ||
      !is.finite(bspline_df) ||
      bspline_df <= 3 ||
      !is.finite(bspline_coefficient_sd) ||
      bspline_coefficient_sd <= 0 ||
      !is.finite(cosine_width_half) ||
      cosine_width_half <= 0 ||
      length(cosine_spike_counts) < 1 ||
      any(!cosine_spike_counts %in% 1:3) ||
      anyDuplicated(cosine_spike_counts) ||
      length(cosine_relative_amplitude_range) != 2 ||
      any(!is.finite(cosine_relative_amplitude_range)) ||
      any(cosine_relative_amplitude_range <= 0) ||
      cosine_relative_amplitude_range[1] >
        cosine_relative_amplitude_range[2] ||
      cosine_relative_amplitude_range[2] >= 1 ||
      !is.finite(cosine_target_centered_rms) ||
      cosine_target_centered_rms <= 0 ||
      !is.finite(switch_threshold) ||
      switch_threshold <= 0 ||
      !is.finite(location_truth_margin) ||
      location_truth_margin <= 0 ||
      !is.finite(switch_truth_margin) ||
      switch_truth_margin <= 0 ||
      !is.finite(non_switch_min_abs) ||
      non_switch_min_abs <= 0 ||
      !is.finite(non_switch_min_range_fraction) ||
      non_switch_min_range_fraction < 0 ||
      non_switch_min_range_fraction >= 1 ||
      !is.finite(max_attempts) ||
      max_attempts < 1 ||
      any(!is.finite(time_grid)) ||
      any(!is.finite(evaluation_grid)) ||
      min(evaluation_grid) < min(time_grid) ||
      max(evaluation_grid) > max(time_grid)) {
    stop("Invalid matched functional-effect simulation settings.")
  }
  if (is.null(scenario)) {
    scenario <- paste0("r3_matched_functional_", truth_mechanism)
  }

  component_seeds <- revision_component_seeds(seed)
  class_seed <- if (is.null(class_seed)) {
    component_seeds[["classes"]]
  } else {
    class_seed
  }
  constant_seed <- if (is.null(constant_seed)) {
    component_seeds[["constant_effects"]]
  } else {
    constant_seed
  }
  shape_seed <- if (is.null(shape_seed)) {
    component_seeds[["functional_truth"]]
  } else {
    shape_seed
  }

  set.seed(class_seed)
  effect_class <- sample_effect_classes(
    J = n_variants,
    class_probs = class_probs,
    exact_class_counts = exact_class_counts
  )
  dynamic_index <- which(effect_class == "dynamic_bspline")
  constant_index <- which(effect_class == "constant")
  if (length(dynamic_index) < 6) {
    stop("The matched functional simulation requires at least six dynamic variants.")
  }

  truth_cells <- expand.grid(
    time_group = c("early", "middle", "late"),
    switch_status = c("switch", "non-switch"),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  truth_cells$truth_group <- paste(
    truth_cells$time_group,
    truth_cells$switch_status,
    sep = " / "
  )
  truth_group_counts <- exact_balanced_counts(
    length(dynamic_index),
    truth_cells$truth_group
  )
  dynamic_truth_group <- sample(
    rep(
      names(truth_group_counts),
      times = unname(truth_group_counts)
    ),
    size = length(dynamic_index),
    replace = FALSE
  )
  dynamic_assignment <- truth_cells[
    match(dynamic_truth_group, truth_cells$truth_group),
    ,
    drop = FALSE
  ]
  dynamic_assignment$spike_count <- NA_integer_
  if (truth_mechanism == "raised_cosine") {
    for (truth_group in truth_cells$truth_group) {
      group_positions <- which(
        dynamic_assignment$truth_group == truth_group
      )
      dynamic_assignment$spike_count[group_positions] <- sample(
        rep(
          cosine_spike_counts,
          length.out = length(group_positions)
        ),
        size = length(group_positions),
        replace = FALSE
      )
    }
  }

  beta_matrix <- matrix(
    0,
    nrow = n_variants,
    ncol = length(time_grid),
    dimnames = list(
      sprintf("variant_%04d", seq_len(n_variants)),
      sprintf("time_%02d", seq_along(time_grid))
    )
  )
  beta_evaluation <- matrix(
    0,
    nrow = n_variants,
    ncol = length(evaluation_grid),
    dimnames = list(
      rownames(beta_matrix),
      sprintf("evaluation_%03d", seq_along(evaluation_grid))
    )
  )

  set.seed(constant_seed)
  if (length(constant_index) > 0) {
    constant_effects <- stats::rnorm(
      length(constant_index),
      mean = 0,
      sd = constant_sd
    )
    beta_matrix[constant_index, ] <- constant_effects
    beta_evaluation[constant_index, ] <- constant_effects
  }

  time_group <- rep(NA_character_, n_variants)
  switch_status <- rep(NA_character_, n_variants)
  truth_group <- rep("dynamic-null", n_variants)
  genetic_main_effect <- rep(NA_real_, n_variants)
  generation_attempt <- rep(NA_integer_, n_variants)
  spike_count <- rep(NA_integer_, n_variants)
  sign_pattern <- rep(NA_character_, n_variants)
  centered_rms <- rep(NA_real_, n_variants)
  minimum_absolute_effect <- rep(NA_real_, n_variants)
  effect_range <- rep(NA_real_, n_variants)
  relative_clearance <- rep(NA_real_, n_variants)
  required_non_switch_clearance <- rep(NA_real_, n_variants)
  peak_centers <- vector("list", n_variants)
  peak_signs <- vector("list", n_variants)
  peak_relative_amplitudes <- vector("list", n_variants)
  bspline_coefficients <- vector("list", n_variants)
  center_ranges <- list(
    early = c(1.50, 3.50),
    middle = c(5.50, 9.50),
    late = c(11.50, 13.50)
  )

  set.seed(shape_seed)
  for (position in seq_along(dynamic_index)) {
    index <- dynamic_index[position]
    assigned_time_group <- dynamic_assignment$time_group[position]
    assigned_switch_status <- dynamic_assignment$switch_status[position]
    assigned_spike_count <- dynamic_assignment$spike_count[position]

    accepted <- FALSE
    for (attempt in seq_len(max_attempts)) {
      main_effect <- stats::rnorm(
        1,
        mean = 0,
        sd = dynamic_main_effect_sd
      )
      centers <- numeric(0)
      signs <- numeric(0)
      relative_amplitudes <- numeric(0)
      coefficients <- numeric(0)
      if (truth_mechanism == "random_bspline") {
        deviation <- simulate_bspline_effect_on_grids(
          time_grid = time_grid,
          evaluation_grid = evaluation_grid,
          amplitude = bspline_amplitude,
          df = bspline_df,
          coefficient_sd = bspline_coefficient_sd
        )
        deviation_observed <- deviation$deviation_observed
        deviation_evaluation <- deviation$deviation_evaluation
        coefficients <- deviation$coefficients
        observed_centered_rms <- sqrt(mean(
          (deviation_observed - mean(deviation_observed))^2
        ))
      } else {
        other_groups <- setdiff(
          c("early", "middle", "late"),
          assigned_time_group
        )
        peak_groups <- c(
          assigned_time_group,
          if (assigned_spike_count > 1L) {
            sample(
              other_groups,
              size = assigned_spike_count - 1L,
              replace = FALSE
            )
          } else {
            character(0)
          }
        )
        centers <- vapply(
          peak_groups,
          function(group) {
            range <- center_ranges[[group]]
            stats::runif(1, min = range[1], max = range[2])
          },
          numeric(1)
        )
        primary_direction <- sample(c(-1, 1), size = 1)
        signs <- if (assigned_switch_status == "non-switch") {
          direction <- if (main_effect == 0) {
            primary_direction
          } else {
            sign(main_effect)
          }
          rep(direction, assigned_spike_count)
        } else if (assigned_spike_count == 1L) {
          direction <- if (main_effect == 0) {
            primary_direction
          } else {
            -sign(main_effect)
          }
          direction
        } else {
          primary_direction * (-1)^(seq_len(assigned_spike_count) - 1L)
        }
        relative_amplitudes <- c(
          1,
          if (assigned_spike_count > 1L) {
            stats::runif(
              assigned_spike_count - 1L,
              min = cosine_relative_amplitude_range[1],
              max = cosine_relative_amplitude_range[2]
            )
          } else {
            numeric(0)
          }
        )
        deviation_observed <- numeric(length(time_grid))
        deviation_evaluation <- numeric(length(evaluation_grid))
        for (peak_index in seq_len(assigned_spike_count)) {
          peak_weight <- signs[peak_index] *
            relative_amplitudes[peak_index]
          deviation_observed <- deviation_observed +
            peak_weight * raised_cosine_peak(
              time_grid,
              center = centers[peak_index],
              half_width = cosine_width_half
            )
          deviation_evaluation <- deviation_evaluation +
            peak_weight * raised_cosine_peak(
              evaluation_grid,
              center = centers[peak_index],
              half_width = cosine_width_half
            )
        }
        observed_centered_rms <- sqrt(mean(
          (
            deviation_observed -
              mean(deviation_observed)
          )^2
        ))
        if (!is.finite(observed_centered_rms) ||
            observed_centered_rms <= .Machine$double.eps) {
          next
        }
        scale_factor <- cosine_target_centered_rms /
          observed_centered_rms
        deviation_observed <- scale_factor * deviation_observed
        deviation_evaluation <- scale_factor * deviation_evaluation
        observed_centered_rms <- sqrt(mean(
          (
            deviation_observed -
              mean(deviation_observed)
          )^2
        ))
      }

      candidate_observed <- main_effect + deviation_observed
      candidate_evaluation <- main_effect + deviation_evaluation
      assessment <- assess_matched_functional_truth(
        beta_evaluation = candidate_evaluation,
        time_group = assigned_time_group,
        switch_status = assigned_switch_status,
        evaluation_grid = evaluation_grid,
        switch_threshold = switch_threshold,
        location_truth_margin = location_truth_margin,
        switch_truth_margin = switch_truth_margin,
        non_switch_min_abs = non_switch_min_abs,
        non_switch_min_range_fraction =
          non_switch_min_range_fraction
      )
      if (!assessment$valid) {
        next
      }

      beta_matrix[index, ] <- candidate_observed
      beta_evaluation[index, ] <- candidate_evaluation
      time_group[index] <- assigned_time_group
      switch_status[index] <- assigned_switch_status
      truth_group[index] <- dynamic_assignment$truth_group[position]
      genetic_main_effect[index] <- main_effect
      generation_attempt[index] <- attempt
      spike_count[index] <- assigned_spike_count
      sign_pattern[index] <- if (truth_mechanism == "random_bspline") {
        "random B-spline"
      } else if (assigned_spike_count == 1L) {
        "single"
      } else if (length(unique(sign(signs))) == 1L) {
        "same-sign"
      } else {
        "alternating-sign"
      }
      centered_rms[index] <- observed_centered_rms
      minimum_absolute_effect[index] <-
        assessment$minimum_absolute_effect
      effect_range[index] <- assessment$effect_range
      relative_clearance[index] <- assessment$relative_clearance
      required_non_switch_clearance[index] <-
        assessment$required_non_switch_clearance
      peak_centers[[index]] <- centers
      peak_signs[[index]] <- signs
      peak_relative_amplitudes[[index]] <- relative_amplitudes
      bspline_coefficients[[index]] <- coefficients
      accepted <- TRUE
      break
    }
    if (!accepted) {
      stop(
        "Could not generate matched functional truth for ",
        assigned_time_group,
        " / ",
        assigned_switch_status,
        " with mechanism ",
        truth_mechanism,
        " within max_attempts."
      )
    }
  }

  true_functionals <- evaluate_temporal_functionals(
    beta_evaluation,
    smooth_var = evaluation_grid,
    switch_threshold = switch_threshold
  )
  dynamic_target_columns <- match(
    time_group[dynamic_index],
    colnames(true_functionals)
  )
  dynamic_target_values <- true_functionals[
    cbind(dynamic_index, dynamic_target_columns)
  ]
  switch_dynamic <- dynamic_index[
    switch_status[dynamic_index] == "switch"
  ]
  non_switch_dynamic <- dynamic_index[
    switch_status[dynamic_index] == "non-switch"
  ]
  non_switch_same_sign <- vapply(
    non_switch_dynamic,
    function(index) {
      all(beta_evaluation[index, ] > 0) ||
        all(beta_evaluation[index, ] < 0)
    },
    logical(1)
  )
  if (any(dynamic_target_values < location_truth_margin) ||
      any(true_functionals[switch_dynamic, "switch"] <
        switch_truth_margin) ||
      any(!non_switch_same_sign) ||
      any(minimum_absolute_effect[non_switch_dynamic] <
        required_non_switch_clearance[non_switch_dynamic])) {
    stop("Matched functional truth failed post-generation validation.")
  }

  unit_info <- data.frame(
    unit_index = seq_len(n_variants),
    unit_id = sprintf("%s_%04d", effect_class, seq_len(n_variants)),
    variant_id = rownames(beta_matrix),
    effect_class = effect_class,
    time_group = time_group,
    switch_status = switch_status,
    truth_group = truth_group,
    genetic_main_effect = genetic_main_effect,
    generation_attempt = generation_attempt,
    truth_mechanism = ifelse(
      effect_class == "dynamic_bspline",
      truth_mechanism,
      "dynamic-null"
    ),
    spike_count = spike_count,
    sign_pattern = sign_pattern,
    centered_rms = centered_rms,
    minimum_absolute_effect = minimum_absolute_effect,
    effect_range = effect_range,
    relative_clearance = relative_clearance,
    required_non_switch_clearance =
      required_non_switch_clearance,
    scenario = scenario,
    stringsAsFactors = FALSE
  )
  unit_info$peak_centers <- I(peak_centers)
  unit_info$peak_signs <- I(peak_signs)
  unit_info$peak_relative_amplitudes <- I(
    peak_relative_amplitudes
  )
  unit_info$bspline_coefficients <- I(bspline_coefficients)

  list(
    beta_matrix = beta_matrix,
    beta_evaluation = beta_evaluation,
    evaluation_grid = evaluation_grid,
    true_functionals = true_functionals,
    unit_info = unit_info,
    group_counts = truth_group_counts,
    settings = list(
      truth_mechanism = truth_mechanism,
      dynamic_main_effect_sd = dynamic_main_effect_sd,
      constant_sd = constant_sd,
      bspline_amplitude = bspline_amplitude,
      bspline_df = bspline_df,
      bspline_coefficient_sd = bspline_coefficient_sd,
      cosine_width_half = cosine_width_half,
      cosine_spike_counts = cosine_spike_counts,
      cosine_relative_amplitude_range =
        cosine_relative_amplitude_range,
      cosine_target_centered_rms = cosine_target_centered_rms,
      switch_threshold = switch_threshold,
      location_truth_margin = location_truth_margin,
      switch_truth_margin = switch_truth_margin,
      non_switch_min_abs = non_switch_min_abs,
      non_switch_min_range_fraction =
        non_switch_min_range_fraction,
      exact_class_counts = exact_class_counts,
      seed = seed,
      class_seed = class_seed,
      constant_seed = constant_seed,
      shape_seed = shape_seed,
      max_attempts = max_attempts
    )
  )
}

sample_labeled_random_bspline_truth <- function(time_group,
                                                switch_status,
                                                time_grid,
                                                evaluation_grid,
                                                amplitude = 2,
                                                df = 6,
                                                degree = 3,
                                                coefficient_sd = 1,
                                                switch_threshold = 0.25,
                                                truth_margin = 0.10,
                                                non_switch_baseline_margin = 0.35,
                                                max_attempts = 10000) {
  allowed_time_groups <- c("early", "middle", "late")
  allowed_switch_status <- c("switch", "non-switch")
  if (!time_group %in% allowed_time_groups || !switch_status %in% allowed_switch_status) {
    stop("Unsupported labelled random B-spline truth group.")
  }
  if (!is.finite(truth_margin) || truth_margin <= 0 ||
      !is.finite(non_switch_baseline_margin) || non_switch_baseline_margin <= 0) {
    stop("truth_margin and non_switch_baseline_margin must be positive and finite.")
  }
  functionals <- make_temporal_functionals(
    smooth_var = evaluation_grid,
    switch_threshold = switch_threshold
  )
  other_time_groups <- setdiff(allowed_time_groups, time_group)

  for (attempt in seq_len(max_attempts)) {
    provisional <- random_bspline_truth_curve(
      time_grid = time_grid,
      evaluation_grid = evaluation_grid,
      amplitude = amplitude,
      df = df,
      degree = degree,
      coefficient_sd = coefficient_sd,
      baseline = 0
    )
    baseline <- if (switch_status == "switch") {
      0
    } else {
      sample(c(-1, 1), size = 1) * (
        max(abs(provisional$deviation_evaluation)) +
          switch_threshold + non_switch_baseline_margin
      )
    }
    curve <- random_bspline_truth_curve(
      time_grid = time_grid,
      evaluation_grid = evaluation_grid,
      amplitude = amplitude,
      df = df,
      degree = degree,
      coefficient_sd = coefficient_sd,
      coefficients = provisional$coefficients,
      baseline = baseline
    )
    contrasts <- vapply(functionals, function(functional) {
      functional(curve$beta_evaluation)
    }, numeric(1))
    time_ok <- contrasts[[time_group]] >= truth_margin &&
      all(contrasts[other_time_groups] <= -truth_margin)
    switch_ok <- if (switch_status == "switch") {
      contrasts[["switch"]] >= truth_margin
    } else {
      contrasts[["switch"]] <= -truth_margin
    }
    if (time_ok && switch_ok) {
      curve$contrasts <- contrasts
      curve$attempt <- attempt
      return(curve)
    }
  }
  stop(
    "Could not generate a labelled ", time_group, " / ", switch_status,
    " random B-spline curve within max_attempts."
  )
}

simulate_labeled_random_bspline_effect_set <- function(n_variants = 1000,
                                                        time_grid = make_time_grid(),
                                                        evaluation_grid = seq(0, 15, by = 0.1),
                                                        class_probs = c(
                                                          dynamic_bspline = 0.20,
                                                          constant = 0.40,
                                                          zero = 0.40
                                                        ),
                                                        dynamic_amplitude = 2,
                                                        bspline_df = 6,
                                                        bspline_degree = 3,
                                                        bspline_coefficient_sd = 1,
                                                        constant_sd = 1,
                                                        switch_threshold = 0.25,
                                                        truth_margin = 0.10,
                                                        non_switch_baseline_margin = 0.35,
                                                        exact_class_counts = TRUE,
                                                        seed = 12345,
                                                        scenario = "genotype_functional_bspline") {
  required_classes <- c("dynamic_bspline", "constant", "zero")
  if (!identical(sort(names(class_probs)), sort(required_classes)) ||
      abs(sum(class_probs) - 1) > 1e-8) {
    stop("class_probs must contain dynamic_bspline, constant, and zero and sum to one.")
  }
  if (!is.finite(seed)) {
    stop("seed must be finite.")
  }
  set.seed(seed)
  effect_class <- sample_effect_classes(
    J = n_variants,
    class_probs = class_probs,
    exact_class_counts = exact_class_counts
  )
  dynamic_index <- which(effect_class == "dynamic_bspline")
  constant_index <- which(effect_class == "constant")
  if (length(dynamic_index) < 6) {
    stop("The labelled functional simulation requires at least six dynamic variants.")
  }

  group_table <- expand.grid(
    time_group = c("early", "middle", "late"),
    switch_status = c("switch", "non-switch"),
    stringsAsFactors = FALSE
  )
  group_table$group_id <- paste(group_table$time_group, group_table$switch_status, sep = "_")
  group_counts <- exact_balanced_counts(length(dynamic_index), group_table$group_id)
  dynamic_group_id <- sample(
    rep(names(group_counts), times = unname(group_counts)),
    size = length(dynamic_index),
    replace = FALSE
  )

  beta_matrix <- matrix(
    0,
    nrow = n_variants,
    ncol = length(time_grid),
    dimnames = list(
      sprintf("variant_%04d", seq_len(n_variants)),
      sprintf("time_%02d", seq_along(time_grid))
    )
  )
  beta_evaluation <- matrix(
    0,
    nrow = n_variants,
    ncol = length(evaluation_grid),
    dimnames = list(
      rownames(beta_matrix),
      sprintf("evaluation_%03d", seq_along(evaluation_grid))
    )
  )
  time_group <- rep(NA_character_, n_variants)
  switch_status <- rep(NA_character_, n_variants)
  generation_attempt <- rep(NA_integer_, n_variants)

  set.seed(seed + 1L)
  if (length(constant_index) > 0) {
    constants <- rnorm(length(constant_index), mean = 0, sd = constant_sd)
    beta_matrix[constant_index, ] <- matrix(
      rep(constants, times = length(time_grid)),
      nrow = length(constant_index),
      ncol = length(time_grid)
    )
    beta_evaluation[constant_index, ] <- matrix(
      rep(constants, times = length(evaluation_grid)),
      nrow = length(constant_index),
      ncol = length(evaluation_grid)
    )
  }

  set.seed(seed + 2L)
  for (position in seq_along(dynamic_index)) {
    group <- group_table[match(dynamic_group_id[position], group_table$group_id), , drop = FALSE]
    truth <- sample_labeled_random_bspline_truth(
      time_group = group$time_group,
      switch_status = group$switch_status,
      time_grid = time_grid,
      evaluation_grid = evaluation_grid,
      amplitude = dynamic_amplitude,
      df = bspline_df,
      degree = bspline_degree,
      coefficient_sd = bspline_coefficient_sd,
      switch_threshold = switch_threshold,
      truth_margin = truth_margin,
      non_switch_baseline_margin = non_switch_baseline_margin
    )
    index <- dynamic_index[position]
    beta_matrix[index, ] <- truth$beta_observed
    beta_evaluation[index, ] <- truth$beta_evaluation
    time_group[index] <- group$time_group
    switch_status[index] <- group$switch_status
    generation_attempt[index] <- truth$attempt
  }

  true_functionals <- evaluate_temporal_functionals(
    beta_evaluation,
    smooth_var = evaluation_grid,
    switch_threshold = switch_threshold
  )
  target_lookup <- match(time_group[dynamic_index], colnames(true_functionals))
  if (any(true_functionals[cbind(dynamic_index, target_lookup)] <= 0) ||
      any(true_functionals[dynamic_index[switch_status[dynamic_index] == "switch"], "switch"] <= 0) ||
      any(true_functionals[dynamic_index[switch_status[dynamic_index] == "non-switch"], "switch"] > 0)) {
    stop("Generated labelled B-spline truth does not satisfy its functional labels.")
  }

  unit_id <- sprintf("%s_%04d", effect_class, seq_len(n_variants))
  unit_info <- data.frame(
    unit_index = seq_len(n_variants),
    unit_id = unit_id,
    variant_id = rownames(beta_matrix),
    effect_class = effect_class,
    time_group = time_group,
    switch_status = switch_status,
    truth_group = ifelse(
      is.na(time_group),
      "dynamic-null",
      paste(time_group, switch_status, sep = " / ")
    ),
    generation_attempt = generation_attempt,
    scenario = scenario,
    stringsAsFactors = FALSE
  )

  list(
    beta_matrix = beta_matrix,
    beta_evaluation = beta_evaluation,
    evaluation_grid = evaluation_grid,
    true_functionals = true_functionals,
    unit_info = unit_info,
    group_counts = group_counts,
    settings = list(
      dynamic_amplitude = dynamic_amplitude,
      bspline_df = bspline_df,
      bspline_degree = bspline_degree,
      bspline_coefficient_sd = bspline_coefficient_sd,
      switch_threshold = switch_threshold,
      truth_margin = truth_margin,
      non_switch_baseline_margin = non_switch_baseline_margin,
      seed = seed
    )
  )
}

simulate_targeted_local_bspline_effect_set <- function(n_variants = 1000,
                                                       time_grid = make_time_grid(),
                                                       evaluation_grid = seq(0, 15, by = 0.1),
                                                       class_probs = c(
                                                         dynamic_bspline = 0.20,
                                                         constant = 0.40,
                                                         zero = 0.40
                                                       ),
                                                       dynamic_amplitude = 2,
                                                       constant_sd = 1,
                                                       switch_threshold = 0.25,
                                                       minimum_location_margin = 0.60,
                                                       minimum_location_ratio = 2,
                                                       non_switch_baseline_fraction = 0.15,
                                                       non_switch_background_fraction = 0.05,
                                                       switch_secondary_fraction = c(0.38, 0.50),
                                                       profile = c("narrow", "broad", "mixed"),
                                                       spiky_truth_version = c(
                                                         "centered_single_v1",
                                                         "mixed_single_double_v2"
                                                       ),
                                                       spiky_secondary_fraction = c(0.40, 0.65),
                                                       spiky_bspline_df = 16,
                                                       spiky_secondary_selection = c(
                                                         "random",
                                                         "nearest"
                                                       ),
                                                       spiky_minimum_peak_separation = 3,
                                                       spiky_non_switch_baseline_fraction = 0,
                                                       target_centered_rms = NULL,
                                                       exact_class_counts = TRUE,
                                                       seed = 12345,
                                                       scenario = "genotype_functional_targeted_local_bspline") {
  profile <- match.arg(profile)
  spiky_truth_version <- match.arg(spiky_truth_version)
  spiky_secondary_selection <- match.arg(spiky_secondary_selection)
  required_classes <- c("dynamic_bspline", "constant", "zero")
  if (!identical(sort(names(class_probs)), sort(required_classes)) ||
      abs(sum(class_probs) - 1) > 1e-8) {
    stop("class_probs must contain dynamic_bspline, constant, and zero and sum to one.")
  }
  if (!is.finite(seed) || !is.finite(constant_sd) || constant_sd <= 0 ||
      length(spiky_secondary_fraction) != 2 ||
      any(!is.finite(spiky_secondary_fraction)) ||
      any(spiky_secondary_fraction <= 0) ||
      spiky_secondary_fraction[1] > spiky_secondary_fraction[2] ||
      spiky_secondary_fraction[2] >= 1 ||
      !is.finite(spiky_bspline_df) ||
      spiky_bspline_df <= 3 ||
      !is.finite(spiky_minimum_peak_separation) ||
      spiky_minimum_peak_separation < 0 ||
      !is.finite(spiky_non_switch_baseline_fraction) ||
      spiky_non_switch_baseline_fraction < 0 ||
      (!is.null(target_centered_rms) &&
        (length(target_centered_rms) != 1 || !is.finite(target_centered_rms) ||
          target_centered_rms <= 0))) {
    stop("seed must be finite and constant_sd must be positive and finite.")
  }
  set.seed(seed)
  effect_class <- sample_effect_classes(
    J = n_variants,
    class_probs = class_probs,
    exact_class_counts = exact_class_counts
  )
  dynamic_index <- which(effect_class == "dynamic_bspline")
  constant_index <- which(effect_class == "constant")
  if (length(dynamic_index) < 6) {
    stop("The labelled functional simulation requires at least six dynamic variants.")
  }

  group_table <- expand.grid(
    time_group = c("early", "middle", "late"),
    switch_status = c("switch", "non-switch"),
    shape_profile = if (profile == "mixed") c("broad", "spiky") else profile,
    stringsAsFactors = FALSE
  )
  group_table$group_id <- paste(
    group_table$time_group,
    group_table$switch_status,
    group_table$shape_profile,
    sep = "_"
  )
  group_counts <- exact_balanced_counts(length(dynamic_index), group_table$group_id)
  dynamic_group_id <- sample(
    rep(names(group_counts), times = unname(group_counts)),
    size = length(dynamic_index),
    replace = FALSE
  )
  planned_spike_count <- rep(NA_integer_, length(dynamic_index))
  if (profile == "mixed" &&
      spiky_truth_version == "mixed_single_double_v2") {
    spiky_switch_groups <- group_table$group_id[
      group_table$shape_profile == "spiky" &
        group_table$switch_status == "switch"
    ]
    planned_spike_count[dynamic_group_id %in% spiky_switch_groups] <- 2L

    spiky_non_switch_groups <- group_table$group_id[
      group_table$shape_profile == "spiky" &
        group_table$switch_status == "non-switch"
    ]
    alternation_offset <- 0L
    for (group_id in spiky_non_switch_groups) {
      positions <- which(dynamic_group_id == group_id)
      alternating_counts <- rep(c(1L, 2L), length.out = length(positions) + 1L)
      if (alternation_offset %% 2L == 1L) {
        alternating_counts <- rep(c(2L, 1L), length.out = length(positions) + 1L)
      }
      alternating_counts <- alternating_counts[seq_along(positions)]
      planned_spike_count[positions] <- sample(
        alternating_counts,
        size = length(positions),
        replace = FALSE
      )
      alternation_offset <- alternation_offset + length(positions)
    }
  }

  beta_matrix <- matrix(
    0,
    nrow = n_variants,
    ncol = length(time_grid),
    dimnames = list(
      sprintf("variant_%04d", seq_len(n_variants)),
      sprintf("time_%02d", seq_along(time_grid))
    )
  )
  beta_evaluation <- matrix(
    0,
    nrow = n_variants,
    ncol = length(evaluation_grid),
    dimnames = list(
      rownames(beta_matrix),
      sprintf("evaluation_%03d", seq_along(evaluation_grid))
    )
  )
  time_group <- rep(NA_character_, n_variants)
  switch_status <- rep(NA_character_, n_variants)
  shape_profile <- rep(NA_character_, n_variants)
  generation_attempt <- rep(NA_integer_, n_variants)
  target_to_outside_ratio <- rep(NA_real_, n_variants)
  effective_sign_transitions <- rep(NA_integer_, n_variants)
  centered_rms <- rep(NA_real_, n_variants)
  spike_count <- rep(NA_integer_, n_variants)
  spike_pattern <- rep(NA_character_, n_variants)

  set.seed(seed + 1L)
  if (length(constant_index) > 0) {
    constants <- rnorm(length(constant_index), mean = 0, sd = constant_sd)
    beta_matrix[constant_index, ] <- matrix(
      rep(constants, times = length(time_grid)),
      nrow = length(constant_index),
      ncol = length(time_grid)
    )
    beta_evaluation[constant_index, ] <- matrix(
      rep(constants, times = length(evaluation_grid)),
      nrow = length(constant_index),
      ncol = length(evaluation_grid)
    )
  }

  set.seed(seed + 2L)
  for (position in seq_along(dynamic_index)) {
    group <- group_table[match(dynamic_group_id[position], group_table$group_id), , drop = FALSE]
    generator_profile <- if (profile == "mixed" && group$shape_profile == "broad") {
      "random_broad"
    } else {
      group$shape_profile
    }
    truth <- if (profile == "mixed" &&
      group$shape_profile == "spiky" &&
      spiky_truth_version == "mixed_single_double_v2") {
      sample_multispike_local_bspline_truth(
        time_group = group$time_group,
        switch_status = group$switch_status,
        time_grid = time_grid,
        evaluation_grid = evaluation_grid,
        spike_count = planned_spike_count[position],
        df = spiky_bspline_df,
        secondary_fraction = spiky_secondary_fraction,
        switch_threshold = switch_threshold,
        minimum_location_margin = minimum_location_margin,
        minimum_location_ratio = minimum_location_ratio,
        target_centered_rms = target_centered_rms,
        minimum_peak_separation = spiky_minimum_peak_separation,
        secondary_selection = spiky_secondary_selection,
        non_switch_baseline_fraction = spiky_non_switch_baseline_fraction
      )
    } else {
      sample_targeted_local_bspline_truth(
        time_group = group$time_group,
        switch_status = group$switch_status,
        time_grid = time_grid,
        evaluation_grid = evaluation_grid,
        amplitude = dynamic_amplitude,
        switch_threshold = switch_threshold,
        minimum_location_margin = minimum_location_margin,
        minimum_location_ratio = minimum_location_ratio,
        non_switch_baseline_fraction = non_switch_baseline_fraction,
        non_switch_background_fraction = non_switch_background_fraction,
        switch_secondary_fraction = switch_secondary_fraction,
        profile = generator_profile,
        target_centered_rms = target_centered_rms
      )
    }
    index <- dynamic_index[position]
    beta_matrix[index, ] <- truth$beta_observed
    beta_evaluation[index, ] <- truth$beta_evaluation
    time_group[index] <- group$time_group
    switch_status[index] <- group$switch_status
    shape_profile[index] <- group$shape_profile
    generation_attempt[index] <- truth$attempt
    target_to_outside_ratio[index] <- truth$diagnostics$target_to_outside_ratio
    effective_sign_transitions[index] <- truth$diagnostics$effective_sign_transitions
    centered_rms[index] <- truth$diagnostics$centered_rms
    if ("spike_count" %in% names(truth$diagnostics)) {
      spike_count[index] <- truth$diagnostics$spike_count
      spike_pattern[index] <- truth$diagnostics$spike_pattern
    } else if (group$shape_profile == "spiky") {
      spike_count[index] <- 1L
      spike_pattern[index] <- "centered single"
    }
  }

  true_functionals <- evaluate_temporal_functionals(
    beta_evaluation,
    smooth_var = evaluation_grid,
    switch_threshold = switch_threshold
  )
  target_lookup <- match(time_group[dynamic_index], colnames(true_functionals))
  target_contrasts <- true_functionals[cbind(dynamic_index, target_lookup)]
  switch_index <- dynamic_index[switch_status[dynamic_index] == "switch"]
  non_switch_index <- dynamic_index[switch_status[dynamic_index] == "non-switch"]
  if (any(target_contrasts < minimum_location_margin) ||
      any(target_to_outside_ratio[dynamic_index] < minimum_location_ratio) ||
      any(true_functionals[switch_index, "switch"] <= 0) ||
      if (profile == "mixed") {
        any(effective_sign_transitions[switch_index] < 1L)
      } else {
        any(effective_sign_transitions[switch_index] != 1L)
      } ||
      any(true_functionals[non_switch_index, "switch"] > 0) ||
      any(effective_sign_transitions[non_switch_index] != 0L) ||
      (!is.null(target_centered_rms) &&
        any(abs(centered_rms[dynamic_index] - target_centered_rms) > 1e-8))) {
    stop("Generated targeted local B-spline truth does not satisfy its constraints.")
  }
  if (profile == "mixed" &&
      spiky_truth_version == "mixed_single_double_v2") {
    spiky_switch_index <- dynamic_index[
      shape_profile[dynamic_index] == "spiky" &
        switch_status[dynamic_index] == "switch"
    ]
    spiky_non_switch_index <- dynamic_index[
      shape_profile[dynamic_index] == "spiky" &
        switch_status[dynamic_index] == "non-switch"
    ]
    if (any(spike_count[spiky_switch_index] != 2L) ||
        any(spike_pattern[spiky_switch_index] != "opposite-sign double") ||
        any(!spike_count[spiky_non_switch_index] %in% 1:2) ||
        any(!spike_pattern[spiky_non_switch_index] %in%
          c("single", "same-sign double")) ||
        abs(sum(spike_count[spiky_non_switch_index] == 1L) -
          sum(spike_count[spiky_non_switch_index] == 2L)) > 1L) {
      stop("The multi-spike truth allocation failed its pattern constraints.")
    }
    for (target_group in c("early", "middle", "late")) {
      target_index <- spiky_non_switch_index[
        time_group[spiky_non_switch_index] == target_group
      ]
      if (abs(sum(spike_count[target_index] == 1L) -
          sum(spike_count[target_index] == 2L)) > 1L) {
        stop("Single and double spikes are not balanced within ", target_group, ".")
      }
    }
  }

  unit_id <- sprintf("%s_%04d", effect_class, seq_len(n_variants))
  unit_info <- data.frame(
    unit_index = seq_len(n_variants),
    unit_id = unit_id,
    variant_id = rownames(beta_matrix),
    effect_class = effect_class,
    time_group = time_group,
    switch_status = switch_status,
    shape_profile = shape_profile,
    truth_group = ifelse(
      is.na(time_group),
      "dynamic-null",
      paste(time_group, switch_status, sep = " / ")
    ),
    truth_stratum = ifelse(
      is.na(time_group),
      "dynamic-null",
      paste(time_group, switch_status, shape_profile, sep = " / ")
    ),
    generation_attempt = generation_attempt,
    target_to_outside_ratio = target_to_outside_ratio,
    effective_sign_transitions = effective_sign_transitions,
    centered_rms = centered_rms,
    spike_count = spike_count,
    spike_pattern = spike_pattern,
    scenario = scenario,
    stringsAsFactors = FALSE
  )

  list(
    beta_matrix = beta_matrix,
    beta_evaluation = beta_evaluation,
    evaluation_grid = evaluation_grid,
    true_functionals = true_functionals,
    unit_info = unit_info,
    group_counts = group_counts,
    settings = list(
      truth_mechanism = "targeted_local_bspline",
      profile = profile,
      dynamic_amplitude = dynamic_amplitude,
      switch_threshold = switch_threshold,
      minimum_location_margin = minimum_location_margin,
      minimum_location_ratio = minimum_location_ratio,
      non_switch_baseline_fraction = non_switch_baseline_fraction,
      non_switch_background_fraction = non_switch_background_fraction,
      switch_secondary_fraction = switch_secondary_fraction,
      spiky_truth_version = spiky_truth_version,
      spiky_secondary_fraction = spiky_secondary_fraction,
      spiky_bspline_df = spiky_bspline_df,
      spiky_secondary_selection = spiky_secondary_selection,
      spiky_minimum_peak_separation = spiky_minimum_peak_separation,
      spiky_non_switch_baseline_fraction = spiky_non_switch_baseline_fraction,
      target_centered_rms = target_centered_rms,
      seed = seed
    )
  )
}

simulate_smooth_abrupt_effect <- function(x,
                                          amplitude = 2,
                                          tau_values = c(5, 8, 11),
                                          width = 0.8) {
  tau <- sample(tau_values, 1)
  direction <- sample(c(-1, 1), 1)
  f <- amplitude * direction * tanh((x - tau) / width)
  f - mean(f)
}

simulate_multi_smooth_abrupt_effect <- function(x,
                                                amplitude = 2,
                                                n_jumps_values = 1:3,
                                                tau_values = 3:12,
                                                width = 0.8) {
  n_jumps <- sample(n_jumps_values, 1)
  n_jumps <- min(n_jumps, length(tau_values))
  taus <- sort(sample(tau_values, n_jumps))
  jump_sizes <- rnorm(n_jumps)
  jump_sizes <- jump_sizes / sqrt(sum(jump_sizes^2))
  f <- rep(0, length(x))
  for (k in seq_along(taus)) {
    f <- f + jump_sizes[k] * tanh((x - taus[k]) / width)
  }
  f <- f - mean(f)
  max_abs <- max(abs(f))
  if (max_abs > 0) {
    f <- amplitude * f / max_abs
  }
  f
}

simulate_local_bspline_transient_effect <- function(x,
                                                    amplitude = 2,
                                                    df = 10,
                                                    degree = 3,
                                                    selected_basis = NULL,
                                                    location_quantile = NULL,
                                                    direction = NULL,
                                                    normalization = c(
                                                      "legacy",
                                                      "center_then_scale"
                                                    )) {
  normalization <- match.arg(normalization)
  x_scaled <- (x - min(x)) / diff(range(x))
  basis <- splines::bs(x_scaled, df = df, degree = degree, intercept = TRUE)
  interior <- seq_len(ncol(basis))
  if (length(interior) > 4) {
    interior <- interior[-c(1, 2, length(interior) - 1, length(interior))]
  }

  if (!is.null(selected_basis) && !is.null(location_quantile)) {
    stop("Specify at most one of selected_basis and location_quantile.")
  }
  if (is.null(selected_basis)) {
    if (is.null(location_quantile)) {
      selected_basis <- sample(interior, 1)
    } else {
      if (length(location_quantile) != 1 ||
          !is.finite(location_quantile) ||
          location_quantile < 0 ||
          location_quantile > 1) {
        stop("location_quantile must be a finite scalar between 0 and 1.")
      }
      location_index <- min(
        length(interior),
        floor(location_quantile * length(interior)) + 1L
      )
      selected_basis <- interior[location_index]
    }
  }
  if (length(selected_basis) != 1 || !selected_basis %in% interior) {
    stop("selected_basis must identify one of the retained interior basis functions.")
  }

  if (is.null(direction)) {
    direction <- sample(c(-1, 1), 1)
  }
  if (length(direction) != 1 || !direction %in% c(-1, 1)) {
    stop("direction must be either -1 or 1.")
  }

  f <- direction * basis[, selected_basis]
  if (normalization == "center_then_scale") {
    f <- f - mean(f)
    max_abs <- max(abs(f))
    if (max_abs > 0) {
      f <- amplitude * f / max_abs
    }
  } else {
    max_abs <- max(abs(f))
    if (max_abs > 0) {
      f <- amplitude * f / max_abs
    }
    f <- f - mean(f)
  }
  f
}

revision_component_seeds <- function(seed) {
  if (length(seed) != 1 || !is.finite(seed)) {
    stop("seed must be a finite scalar.")
  }
  seed <- as.double(seed)
  offsets <- c(
    genotype = 101,
    covariates = 211,
    classes = 307,
    constant_effects = 401,
    transient_locations = 503,
    transient_directions = 601,
    expression = 701,
    functional_truth = 1301,
    functional_posterior = 1901,
    permutations = 10007
  )
  max_seed <- .Machine$integer.max - 1
  out <- ((seed + offsets - 1) %% max_seed) + 1
  stats::setNames(as.integer(out), names(offsets))
}

simulate_paired_local_bspline_effect_sets <- function(candidate_settings,
                                                       n_variants = 1000,
                                                       time_grid = make_time_grid(),
                                                       class_probs = c(
                                                         dynamic_local_bspline_transient = 0.20,
                                                         constant = 0.40,
                                                         zero = 0.40
                                                       ),
                                                       constant_sd = 1,
                                                       transient_bspline_degree = 3,
                                                       exact_class_counts = TRUE,
                                                       seed = 12345,
                                                       scenario_prefix = "genotype_spiky_signal_pilot") {
  required_columns <- c(
    "setting_id",
    "transient_bspline_df",
    "dynamic_amplitude",
    "normalization"
  )
  missing_columns <- setdiff(required_columns, colnames(candidate_settings))
  if (length(missing_columns) > 0) {
    stop(
      "candidate_settings is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  if (anyDuplicated(candidate_settings$setting_id)) {
    stop("candidate_settings$setting_id must be unique.")
  }
  allowed_classes <- c(
    "dynamic_local_bspline_transient",
    "constant",
    "zero"
  )
  if (!identical(sort(names(class_probs)), sort(allowed_classes)) ||
      abs(sum(class_probs) - 1) > 1e-8) {
    stop("class_probs must define dynamic transient, constant, and zero classes and sum to 1.")
  }
  if (any(!candidate_settings$normalization %in%
          c("legacy", "center_then_scale"))) {
    stop("Unsupported local B-spline normalization.")
  }
  if (any(candidate_settings$transient_bspline_df <= transient_bspline_degree) ||
      any(candidate_settings$dynamic_amplitude <= 0)) {
    stop("Candidate degrees of freedom and amplitudes must be positive and valid.")
  }

  seeds <- revision_component_seeds(seed)
  set.seed(seeds[["classes"]])
  effect_class <- sample_effect_classes(
    J = n_variants,
    class_probs = class_probs,
    exact_class_counts = exact_class_counts
  )
  dynamic_index <- which(effect_class == "dynamic_local_bspline_transient")
  constant_index <- which(effect_class == "constant")

  set.seed(seeds[["constant_effects"]])
  constant_effects <- rnorm(length(constant_index), mean = 0, sd = constant_sd)
  set.seed(seeds[["transient_locations"]])
  transient_locations <- runif(length(dynamic_index))
  set.seed(seeds[["transient_directions"]])
  transient_directions <- sample(c(-1, 1), length(dynamic_index), replace = TRUE)

  variant_id <- sprintf("variant_%04d", seq_len(n_variants))
  unit_id <- sprintf("%s_%04d", effect_class, seq_len(n_variants))
  effect_sets <- vector("list", nrow(candidate_settings))
  names(effect_sets) <- candidate_settings$setting_id

  for (setting_row in seq_len(nrow(candidate_settings))) {
    setting <- candidate_settings[setting_row, , drop = FALSE]
    beta_matrix <- matrix(
      0,
      nrow = n_variants,
      ncol = length(time_grid),
      dimnames = list(variant_id, sprintf("time_%02d", seq_along(time_grid)))
    )
    if (length(constant_index) > 0) {
      beta_matrix[constant_index, ] <- matrix(
        rep(constant_effects, times = length(time_grid)),
        nrow = length(constant_index),
        ncol = length(time_grid)
      )
    }
    for (dynamic_row in seq_along(dynamic_index)) {
      beta_matrix[dynamic_index[dynamic_row], ] <-
        simulate_local_bspline_transient_effect(
          x = time_grid,
          amplitude = setting$dynamic_amplitude,
          df = setting$transient_bspline_df,
          degree = transient_bspline_degree,
          location_quantile = transient_locations[dynamic_row],
          direction = transient_directions[dynamic_row],
          normalization = setting$normalization
        )
    }

    scenario <- paste(scenario_prefix, setting$setting_id, sep = "_")
    unit_info <- data.frame(
      unit_index = seq_len(n_variants),
      unit_id = unit_id,
      variant_id = variant_id,
      effect_class = effect_class,
      scenario = scenario,
      stringsAsFactors = FALSE
    )
    effect_sets[[setting$setting_id]] <- list(
      beta_matrix = beta_matrix,
      unit_info = unit_info,
      settings = as.list(setting)
    )
  }

  list(
    effect_sets = effect_sets,
    candidate_settings = candidate_settings,
    seeds = seeds,
    latent = list(
      effect_class = effect_class,
      constant_index = constant_index,
      constant_effects = constant_effects,
      dynamic_index = dynamic_index,
      transient_locations = transient_locations,
      transient_directions = transient_directions
    )
  )
}

summarize_dynamic_effect_geometry <- function(effect_sim,
                                              dynamic_class =
                                                "dynamic_local_bspline_transient") {
  dynamic_index <- which(effect_sim$unit_info$effect_class == dynamic_class)
  if (length(dynamic_index) == 0) {
    stop("effect_sim contains no requested dynamic effects.")
  }
  beta <- effect_sim$beta_matrix[dynamic_index, , drop = FALSE]
  max_abs <- apply(abs(beta), 1, max)
  rms <- sqrt(rowMeans(beta^2))
  half_support <- rowSums(abs(beta) >= max_abs / 2)
  quarter_support <- rowSums(abs(beta) >= max_abs / 4)
  total_variation <- apply(beta, 1, function(x) sum(abs(diff(x))))

  data.frame(
    n_dynamic = length(dynamic_index),
    max_abs_min = min(max_abs),
    max_abs_median = stats::median(max_abs),
    max_abs_max = max(max_abs),
    rms_min = min(rms),
    rms_median = stats::median(rms),
    rms_max = max(rms),
    half_support_median = stats::median(half_support),
    half_support_max = max(half_support),
    quarter_support_median = stats::median(quarter_support),
    total_variation_median = stats::median(total_variation),
    stringsAsFactors = FALSE
  )
}

simulate_iwp_dataset <- function(effect_class,
                                 sd_vec,
                                 time_grid,
                                 correlation = NULL) {
  if (!requireNamespace("fashr", quietly = TRUE)) {
    stop("The fashr package is required for IWP simulation.")
  }

  type <- switch(
    effect_class,
    constant = "nondynamic",
    linear = "linear",
    quadratic = "quadratic",
    nonlinear_iwp2 = "nonlinear",
    stop("Unsupported IWP effect class: ", effect_class)
  )

  sim <- fashr::simulate_process(
    x = time_grid,
    sd_poly = if (effect_class == "nonlinear_iwp2") 0 else 1,
    type = type,
    sd = sd_vec,
    sd_fun = 5,
    p = if (effect_class == "nonlinear_iwp2") 2 else 1,
    normalize = FALSE
  )

  if (!is.null(correlation)) {
    sim$y <- sim$truef + draw_correlated_noise(sim$sd, correlation)
  }

  data.frame(
    x = sim$x,
    y = sim$y,
    sd = sim$sd,
    truef = sim$truef,
    effect_class = effect_class,
    stringsAsFactors = FALSE
  )
}

simulate_shape_dataset <- function(effect_class,
                                   sd_vec,
                                   time_grid,
                                   amplitude = 1,
                                   bspline_df = 6,
                                   bspline_coefficient_sd = 1,
                                   smooth_abrupt_width = 0.8,
                                   smooth_abrupt_tau_values = c(5, 8, 11),
                                   multi_smooth_abrupt_n_jumps_values = 1:3,
                                   multi_smooth_abrupt_tau_values = 3:12,
                                   transient_bspline_df = 10,
                                   transient_bspline_degree = 3,
                                   correlation = NULL) {
  truef <- if (effect_class == "bspline") {
    simulate_bspline_effect(
      x = time_grid,
      amplitude = amplitude,
      df = bspline_df,
      coefficient_sd = bspline_coefficient_sd
    )
  } else if (effect_class == "smooth_abrupt") {
    simulate_smooth_abrupt_effect(
      x = time_grid,
      amplitude = amplitude,
      tau_values = smooth_abrupt_tau_values,
      width = smooth_abrupt_width
    )
  } else if (effect_class == "multi_smooth_abrupt") {
    simulate_multi_smooth_abrupt_effect(
      x = time_grid,
      amplitude = amplitude,
      n_jumps_values = multi_smooth_abrupt_n_jumps_values,
      tau_values = multi_smooth_abrupt_tau_values,
      width = smooth_abrupt_width
    )
  } else if (effect_class == "local_bspline_transient") {
    simulate_local_bspline_transient_effect(
      x = time_grid,
      amplitude = amplitude,
      df = transient_bspline_df,
      degree = transient_bspline_degree
    )
  } else {
    effect_function(effect_class, time_grid, amplitude = amplitude)
  }
  y <- truef + draw_correlated_noise(sd_vec, correlation)
  data.frame(
    x = time_grid,
    y = y,
    sd = sd_vec,
    truef = truef,
    effect_class = effect_class,
    stringsAsFactors = FALSE
  )
}

simulate_revision_datasets <- function(J,
                                       class_probs = c(
                                         constant = 0.80,
                                         linear = 0.10,
                                         nonlinear_iwp2 = 0.10
                                       ),
                                       time_grid = make_time_grid(),
                                       sd_values = c(0.1, 0.3, 0.5),
                                       scenario = "appendixB_like",
                                       correlation = NULL,
                                       amplitude = 1,
                                       bspline_df = 6,
                                       bspline_coefficient_sd = 1,
                                       smooth_abrupt_width = 0.8,
                                       smooth_abrupt_tau_values = c(5, 8, 11),
                                       multi_smooth_abrupt_n_jumps_values = 1:3,
                                       multi_smooth_abrupt_tau_values = 3:12,
                                       transient_bspline_df = 10,
                                       transient_bspline_degree = 3,
                                       exact_class_counts = TRUE,
                                       seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (abs(sum(class_probs) - 1) > 1e-8) {
    stop("class_probs must sum to 1.")
  }

  classes <- sample_effect_classes(
    J = J,
    class_probs = class_probs,
    exact_class_counts = exact_class_counts
  )
  datasets <- vector("list", J)

  for (j in seq_len(J)) {
    sd_vec <- sample(sd_values, size = length(time_grid), replace = TRUE)
    effect_class <- classes[j]
    if (effect_class %in% c("constant", "linear", "quadratic", "nonlinear_iwp2")) {
      dat <- simulate_iwp_dataset(effect_class, sd_vec, time_grid, correlation)
    } else {
      dat <- simulate_shape_dataset(
        effect_class = effect_class,
        sd_vec = sd_vec,
        time_grid = time_grid,
        amplitude = amplitude,
        bspline_df = bspline_df,
        bspline_coefficient_sd = bspline_coefficient_sd,
        smooth_abrupt_width = smooth_abrupt_width,
        smooth_abrupt_tau_values = smooth_abrupt_tau_values,
        multi_smooth_abrupt_n_jumps_values = multi_smooth_abrupt_n_jumps_values,
        multi_smooth_abrupt_tau_values = multi_smooth_abrupt_tau_values,
        transient_bspline_df = transient_bspline_df,
        transient_bspline_degree = transient_bspline_degree,
        correlation = correlation
      )
    }
    dat$unit_id <- sprintf("%s_%04d", effect_class, j)
    dat$scenario <- scenario
    datasets[[j]] <- dat
  }

  names(datasets) <- vapply(datasets, function(dat) dat$unit_id[1], character(1))
  attr(datasets, "unit_info") <- data.frame(
    unit_index = seq_len(J),
    unit_id = names(datasets),
    effect_class = classes,
    scenario = scenario,
    stringsAsFactors = FALSE
  )
  datasets
}

is_null_class <- function(effect_class, target = c("dynamic", "nonlinear")) {
  target <- match.arg(target)
  if (target == "dynamic") {
    return(effect_class %in% c("constant", "zero"))
  }
  effect_class %in% c("constant", "linear", "zero")
}

fit_fash_for_revision <- function(datasets,
                                  orders = c(1, 2),
                                  grid = default_revision_grid(),
                                  num_basis = 20,
                                  penalty = 10,
                                  pred_step = 1,
                                  num_cores = 1,
                                  apply_bf = TRUE,
                                  verbose = FALSE) {
  if (!requireNamespace("fashr", quietly = TRUE)) {
    stop("The fashr package is required to fit FASH.")
  }
  if (length(pred_step) != 1L ||
      !is.finite(pred_step) ||
      pred_step <= 0) {
    stop("pred_step must be one positive finite value.")
  }

  fits <- list()
  for (order in orders) {
    fit <- fashr::fash(
      Y = "y",
      smooth_var = "x",
      S = "sd",
      data_list = datasets,
      order = order,
      grid = grid,
      num_basis = num_basis,
      penalty = penalty,
      pred_step = pred_step,
      num_cores = num_cores,
      verbose = verbose
    )
    fits[[paste0("fash_iwp", order, "_raw")]] <- fit
    if (apply_bf) {
      fits[[paste0("fash_iwp", order, "_bf")]] <- fashr::BF_update(fit, plot = FALSE)
    }
  }
  fits
}

get_fash_lfdr <- function(fit) {
  if (!is.null(fit$lfdr)) {
    return(fit$lfdr)
  }
  if (!is.null(fit$posterior_weights)) {
    return(fit$posterior_weights[, 1])
  }
  stop("Could not find lfdr values in the FASH fit.")
}

get_fash_fdr_table <- function(fit) {
  result <- NULL
  invisible(utils::capture.output(
    result <- fashr::fdr_control(fit, alpha = 1, plot = FALSE)
  ))
  fdr_table <- result$fdr_results
  if (is.null(fdr_table) || !all(c("index", "FDR") %in% names(fdr_table))) {
    stop("fashr::fdr_control did not return index and FDR columns.")
  }
  fdr_table
}

functional_cfsr_table <- function(indices, lfsr) {
  indices <- as.integer(indices)
  lfsr <- as.numeric(lfsr)
  if (length(indices) != length(lfsr) || any(!is.finite(lfsr)) ||
      any(lfsr < 0) || any(lfsr > 1)) {
    stop("indices and lfsr must have the same length, with lfsr in [0, 1].")
  }
  if (length(indices) == 0) {
    return(data.frame(index = integer(), lfsr = numeric(), cfsr = numeric()))
  }
  out <- data.frame(index = indices, lfsr = lfsr)
  out <- out[order(out$lfsr, out$index), , drop = FALSE]
  out$cfsr <- cumsum(out$lfsr) / seq_len(nrow(out))
  out
}

compute_functional_lfsr <- function(fit,
                                    functionals,
                                    indices,
                                    smooth_var,
                                    num_cores = 1,
                                    seed = NULL) {
  if (!is.list(functionals) || length(functionals) == 0 ||
      is.null(names(functionals)) || any(!nzchar(names(functionals)))) {
    stop("functionals must be a named non-empty list.")
  }
  indices <- sort(unique(as.integer(indices)))
  if (length(indices) == 0) {
    return(stats::setNames(vector("list", length(functionals)), names(functionals)))
  }
  if (any(indices < 1)) {
    stop("indices must be positive integers.")
  }

  output <- vector("list", length(functionals))
  names(output) <- names(functionals)
  for (k in seq_along(functionals)) {
    if (!is.null(seed)) {
      set.seed(as.integer(seed + k - 1L))
    }
    result <- fashr::testing_functional(
      functional = functionals[[k]],
      lfsr_cal = function(x) mean(x <= 0),
      fash = fit,
      indices = indices,
      smooth_var = smooth_var,
      num_cores = num_cores
    )
    if (!all(c("indices", "lfsr") %in% names(result)) ||
        nrow(result) != length(indices)) {
      stop("testing_functional returned an incomplete result.")
    }
    output[[k]] <- stats::setNames(result$lfsr, as.character(result$indices))
  }
  output
}

evaluate_fash_functional_testing <- function(fit,
                                             true_functionals,
                                             evaluation_grid,
                                             alpha_grid = seq(0.005, 0.20, by = 0.005),
                                             method = "FASH-IWP1-BF",
                                             scenario = "genotype_functional_bspline",
                                             switch_threshold = 0.25,
                                             switch_minimum_duration = 0,
                                             true_dynamic = NULL,
                                             num_cores = 1,
                                             seed = NULL) {
  true_functionals <- as.matrix(true_functionals)
  storage.mode(true_functionals) <- "numeric"
  required_targets <- c("early", "middle", "late", "switch")
  if (!all(required_targets %in% colnames(true_functionals))) {
    stop("true_functionals must contain early, middle, late, and switch columns.")
  }
  if (is.null(true_dynamic)) {
    true_dynamic <- rep(TRUE, nrow(true_functionals))
  }
  true_dynamic <- as.logical(true_dynamic)
  if (length(true_dynamic) != nrow(true_functionals) || anyNA(true_dynamic)) {
    stop("true_dynamic must contain one non-missing logical value per unit.")
  }
  alpha_grid <- sort(unique(as.numeric(alpha_grid)))
  if (length(alpha_grid) == 0 || any(!is.finite(alpha_grid)) ||
      any(alpha_grid <= 0) || any(alpha_grid > 1)) {
    stop("alpha_grid must contain finite values in (0, 1].")
  }

  fdr_table <- get_fash_fdr_table(fit)
  maximum_dynamic_indices <- sort(unique(fdr_table$index[fdr_table$FDR <= max(alpha_grid)]))
  functionals <- make_temporal_functionals(
    smooth_var = evaluation_grid,
    switch_threshold = switch_threshold,
    switch_minimum_duration = switch_minimum_duration
  )
  lfsr_by_target <- compute_functional_lfsr(
    fit = fit,
    functionals = functionals,
    indices = maximum_dynamic_indices,
    smooth_var = evaluation_grid,
    num_cores = num_cores,
    seed = seed
  )

  rows <- vector("list", length(alpha_grid) * length(required_targets))
  row_index <- 1L
  for (alpha in alpha_grid) {
    dynamic_indices <- sort(unique(fdr_table$index[fdr_table$FDR <= alpha]))
    if (!all(dynamic_indices %in% maximum_dynamic_indices)) {
      stop("Dynamic selections are not nested across the requested alpha grid.")
    }
    for (target in required_targets) {
      lfsr_map <- lfsr_by_target[[target]]
      candidate_lfsr <- if (length(dynamic_indices) == 0) {
        numeric()
      } else {
        unname(lfsr_map[as.character(dynamic_indices)])
      }
      if (anyNA(candidate_lfsr)) {
        stop("Missing functional lfsr values for dynamically selected variants.")
      }
      cfsr_table <- functional_cfsr_table(dynamic_indices, candidate_lfsr)
      selected_indices <- cfsr_table$index[cfsr_table$cfsr <= alpha]
      true_null <- true_functionals[, target] <= 0
      discoveries <- length(selected_indices)
      false_discoveries <- if (discoveries == 0) 0L else {
        sum(true_null[selected_indices])
      }
      true_positives <- if (discoveries == 0) 0L else {
        sum(!true_null[selected_indices])
      }
      conditional_indices <- selected_indices[true_dynamic[selected_indices]]
      conditional_discoveries <- length(conditional_indices)
      conditional_false_discoveries <- if (conditional_discoveries == 0) 0L else {
        sum(true_null[conditional_indices])
      }
      first_stage_null_calls <- if (discoveries == 0) 0L else {
        sum(!true_dynamic[selected_indices])
      }
      true_alternatives <- sum(!true_null)
      selected_lfsr <- if (discoveries == 0) numeric() else {
        cfsr_table$lfsr[cfsr_table$index %in% selected_indices]
      }
      rows[[row_index]] <- data.frame(
        scenario = scenario,
        target = target,
        method = method,
        alpha = alpha,
        dynamic_discoveries = length(dynamic_indices),
        n_discoveries = discoveries,
        false_discoveries = false_discoveries,
        conditional_discoveries = conditional_discoveries,
        conditional_false_discoveries = conditional_false_discoveries,
        first_stage_null_calls = first_stage_null_calls,
        true_positives = true_positives,
        estimated_fsr = if (discoveries == 0) 0 else mean(selected_lfsr),
        empirical_fsr = if (discoveries == 0) 0 else false_discoveries / discoveries,
        conditional_empirical_fsr = if (conditional_discoveries == 0) {
          0
        } else {
          conditional_false_discoveries / conditional_discoveries
        },
        power = if (true_alternatives == 0) NA_real_ else true_positives / true_alternatives,
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
  alpha_curve <- do.call(rbind, rows)
  rownames(alpha_curve) <- NULL
  list(
    alpha_curve = alpha_curve,
    fdr_table = fdr_table,
    maximum_dynamic_indices = maximum_dynamic_indices,
    lfsr_by_target = lfsr_by_target
  )
}

solve_from_chol <- function(chol_mat, b) {
  backsolve(chol_mat, forwardsolve(t(chol_mat), b))
}

log_marginal_common_intercept <- function(y, se, x, slope_sd = 0) {
  n <- length(y)
  if (length(se) != n || length(x) != n) {
    stop("y, se, and x must have the same length.")
  }
  if (any(se <= 0)) {
    stop("All standard errors must be positive.")
  }

  cov_mat <- diag(se^2, nrow = n)
  if (slope_sd > 0) {
    cov_mat <- cov_mat + slope_sd^2 * tcrossprod(x)
  }

  chol_cov <- chol(cov_mat + diag(1e-10, n))
  logdet_cov <- 2 * sum(log(diag(chol_cov)))

  y_mat <- matrix(y, ncol = 1)
  intercept <- matrix(1, nrow = n, ncol = 1)
  cov_inv_y <- solve_from_chol(chol_cov, y_mat)
  cov_inv_intercept <- solve_from_chol(chol_cov, intercept)

  xt_cinv_x <- crossprod(intercept, cov_inv_intercept)
  xt_cinv_y <- crossprod(intercept, cov_inv_y)
  fitted_quad <- crossprod(xt_cinv_y, solve(xt_cinv_x, xt_cinv_y))
  quad <- as.numeric(crossprod(y_mat, cov_inv_y) - fitted_quad)
  logdet_xt_cinv_x <- as.numeric(determinant(xt_cinv_x, logarithm = TRUE)$modulus)

  -0.5 * (
    (n - 1) * log(2 * pi) +
      logdet_cov +
      logdet_xt_cinv_x +
      quad
  )
}

log_sum_exp_pair <- function(a, b) {
  m <- pmax(a, b)
  m + log(exp(a - m) + exp(b - m))
}

fit_two_component_pi <- function(log_m0, log_m1, eps = 1e-6) {
  if (length(log_m0) != length(log_m1)) {
    stop("log_m0 and log_m1 must have the same length.")
  }

  loglik_at <- function(pi0) {
    sum(log_sum_exp_pair(log(pi0) + log_m0, log1p(-pi0) + log_m1))
  }

  opt <- optimize(
    f = function(pi0) -loglik_at(pi0),
    interval = c(eps, 1 - eps)
  )

  candidates <- c(eps, opt$minimum, 1 - eps)
  loglik <- vapply(candidates, loglik_at, numeric(1))
  best <- which.max(loglik)

  list(pi0 = candidates[best], loglik = loglik[best])
}

fit_simplified_fash <- function(datasets,
                                sigma_beta = 1,
                                estimate_sigma = FALSE,
                                sigma_beta_grid = exp(seq(log(0.05), log(5), length.out = 25)),
                                scale_time = TRUE) {
  candidate_sigmas <- if (estimate_sigma) sigma_beta_grid else sigma_beta
  candidate_sigmas <- sort(unique(candidate_sigmas))
  if (any(candidate_sigmas <= 0)) {
    stop("All sigma_beta values must be positive.")
  }

  get_scaled_time <- function(dat) {
    x <- dat$x
    if (scale_time) {
      x <- (x - min(x)) / diff(range(x))
    }
    x
  }

  log_m0 <- vapply(
    datasets,
    function(dat) {
      log_marginal_common_intercept(
        y = dat$y,
        se = dat$sd,
        x = get_scaled_time(dat),
        slope_sd = 0
      )
    },
    numeric(1)
  )

  sigma_fits <- lapply(candidate_sigmas, function(sigma) {
    log_m1 <- vapply(
      datasets,
      function(dat) {
        log_marginal_common_intercept(
          y = dat$y,
          se = dat$sd,
          x = get_scaled_time(dat),
          slope_sd = sigma
        )
      },
      numeric(1)
    )
    pi_fit <- fit_two_component_pi(log_m0, log_m1)
    list(
      sigma_beta = sigma,
      log_m0 = log_m0,
      log_m1 = log_m1,
      pi0 = pi_fit$pi0,
      loglik = pi_fit$loglik
    )
  })

  best_idx <- which.max(vapply(sigma_fits, `[[`, numeric(1), "loglik"))
  best <- sigma_fits[[best_idx]]
  sigma_profile <- data.frame(
    sigma_beta = vapply(sigma_fits, `[[`, numeric(1), "sigma_beta"),
    estimated_pi0 = vapply(sigma_fits, `[[`, numeric(1), "pi0"),
    loglik = vapply(sigma_fits, `[[`, numeric(1), "loglik"),
    selected = seq_along(sigma_fits) == best_idx,
    grid_boundary = seq_along(sigma_fits) %in%
      c(1L, length(sigma_fits)),
    stringsAsFactors = FALSE
  )
  log_denom <- log_sum_exp_pair(
    log(best$pi0) + best$log_m0,
    log1p(-best$pi0) + best$log_m1
  )
  lfdr <- exp(log(best$pi0) + best$log_m0 - log_denom)

  fit <- list(
    prior_weights = data.frame(
      component = c("constant", "linear_slope"),
      prior_weight = c(best$pi0, 1 - best$pi0),
      stringsAsFactors = FALSE
    ),
    sigma_beta = best$sigma_beta,
    sigma_beta_candidates = candidate_sigmas,
    estimate_sigma = estimate_sigma,
    scale_time = scale_time,
    sigma_profile = sigma_profile,
    selected_sigma_on_boundary = sigma_profile$grid_boundary[best_idx],
    log_marginal = data.frame(
      unit_id = names(datasets),
      log_m0 = best$log_m0,
      log_m1 = best$log_m1,
      stringsAsFactors = FALSE
    ),
    lfdr = lfdr,
    posterior_weights = cbind(null = lfdr, alternative = 1 - lfdr),
    loglik = best$loglik
  )
  class(fit) <- "simplified_fash"
  fit
}

validate_simplified_sigma_profile <- function(fit,
                                              require_interior = FALSE) {
  if (!inherits(fit, "simplified_fash")) {
    stop("fit must be a simplified_fash object.")
  }
  required_fit_fields <- c(
    "sigma_beta", "sigma_beta_candidates", "sigma_profile",
    "selected_sigma_on_boundary", "loglik"
  )
  if (!all(required_fit_fields %in% names(fit))) {
    stop("The simplified-FASH fit is missing sigma-profile fields.")
  }

  profile <- fit$sigma_profile
  required_profile_columns <- c(
    "sigma_beta", "estimated_pi0", "loglik", "selected",
    "grid_boundary"
  )
  if (!is.data.frame(profile) ||
      !all(required_profile_columns %in% names(profile)) ||
      nrow(profile) != length(fit$sigma_beta_candidates)) {
    stop("The simplified-FASH sigma profile has invalid dimensions.")
  }
  if (any(!is.finite(profile$sigma_beta)) ||
      any(profile$sigma_beta <= 0) ||
      any(!is.finite(profile$estimated_pi0)) ||
      any(profile$estimated_pi0 < 0 | profile$estimated_pi0 > 1) ||
      any(!is.finite(profile$loglik)) ||
      anyNA(profile$selected) ||
      anyNA(profile$grid_boundary)) {
    stop("The simplified-FASH sigma profile contains invalid values.")
  }
  if (!isTRUE(all.equal(
    profile$sigma_beta,
    fit$sigma_beta_candidates,
    tolerance = 0
  ))) {
    stop("The sigma profile does not match the recorded candidate grid.")
  }

  selected_index <- which(profile$selected)
  if (length(selected_index) != 1L ||
      !isTRUE(all.equal(
        fit$sigma_beta,
        profile$sigma_beta[selected_index],
        tolerance = 1e-12
      )) ||
      !isTRUE(all.equal(
        fit$loglik,
        profile$loglik[selected_index],
        tolerance = 1e-8
      )) ||
      !isTRUE(all.equal(
        profile$loglik[selected_index],
        max(profile$loglik),
        tolerance = 1e-8
      ))) {
    stop("The selected sigma-profile row does not match the fitted model.")
  }

  expected_boundary <- selected_index %in% c(1L, nrow(profile))
  if (!identical(
        as.logical(profile$grid_boundary[selected_index]),
        as.logical(expected_boundary)
      ) ||
      !identical(
        as.logical(fit$selected_sigma_on_boundary),
        as.logical(expected_boundary)
      )) {
    stop("The selected sigma boundary indicator is inconsistent.")
  }
  if (require_interior && expected_boundary) {
    stop("The selected slope scale lies on the candidate-grid boundary.")
  }

  invisible(TRUE)
}

bf_control_from_bf <- function(BF) {
  if (all(is.na(BF)) || all(is.nan(BF))) {
    stop("Bayes factors contain only NA or NaN values.")
  }
  if (any(is.na(BF)) || any(is.nan(BF))) {
    BF <- BF[!is.na(BF) & !is.nan(BF)]
  }

  BF_sorted <- sort(BF, decreasing = FALSE)
  mu <- cumsum(BF_sorted) / seq_along(BF_sorted)
  pi0_hat <- seq_along(BF_sorted) / length(BF_sorted)
  pi0_hat_star <- if (max(mu, na.rm = TRUE) < 1) {
    1
  } else {
    pi0_hat[which(mu >= 1)[1]]
  }

  list(mu = mu, pi0_hat = pi0_hat, pi0_hat_star = pi0_hat_star)
}

BF_update_simplified_fash <- function(fit) {
  if (!inherits(fit, "simplified_fash")) {
    stop("fit must be a simplified_fash object.")
  }

  BF <- exp(fit$log_marginal$log_m1 - fit$log_marginal$log_m0)
  pi0_bf <- bf_control_from_bf(BF)$pi0_hat_star
  log_denom <- log_sum_exp_pair(
    log(pi0_bf) + fit$log_marginal$log_m0,
    log1p(-pi0_bf) + fit$log_marginal$log_m1
  )
  lfdr <- exp(log(pi0_bf) + fit$log_marginal$log_m0 - log_denom)

  fit$prior_weights <- data.frame(
    component = c("constant", "linear_slope"),
    prior_weight = c(pi0_bf, 1 - pi0_bf),
    stringsAsFactors = FALSE
  )
  fit$posterior_weights <- cbind(null = lfdr, alternative = 1 - lfdr)
  fit$lfdr <- lfdr
  fit$BF <- BF
  fit$bf_adjusted <- TRUE
  fit
}

validate_linear_mixture_grid <- function(grid) {
  grid <- as.numeric(grid)
  if (length(grid) < 2L ||
      any(!is.finite(grid)) ||
      any(grid < 0) ||
      anyDuplicated(grid) ||
      grid[1] != 0 ||
      sum(grid == 0) != 1L ||
      any(diff(grid) <= 0)) {
    stop(
      "grid must be strictly increasing with one exact zero as its first value."
    )
  }
  grid
}

expand_grid_prior_weights <- function(prior_weights, grid) {
  grid <- validate_linear_mixture_grid(grid)
  if (!is.data.frame(prior_weights) ||
      !all(c("psd", "prior_weight") %in% names(prior_weights)) ||
      nrow(prior_weights) < 1L) {
    stop("prior_weights must contain psd and prior_weight columns.")
  }
  psd <- as.numeric(prior_weights$psd)
  weight <- as.numeric(prior_weights$prior_weight)
  if (any(!is.finite(psd)) ||
      any(!is.finite(weight)) ||
      any(weight < 0) ||
      anyDuplicated(psd) ||
      abs(sum(weight) - 1) > 1e-6) {
    stop("prior_weights contains invalid scales or weights.")
  }

  result <- numeric(length(grid))
  names(result) <- as.character(grid)
  index <- match(psd, grid)
  if (anyNA(index)) {
    stop("prior_weights contains scales outside the grid.")
  }
  result[index] <- weight
  result
}

compute_linear_mixture_log_likelihood <- function(datasets,
                                                  grid = default_revision_grid(),
                                                  pred_step = 1) {
  grid <- validate_linear_mixture_grid(grid)
  if (length(pred_step) != 1L ||
      !is.finite(pred_step) ||
      pred_step <= 0) {
    stop("pred_step must be one positive finite value.")
  }
  if (!is.list(datasets) || length(datasets) == 0L) {
    stop("datasets must be a nonempty list.")
  }
  if (is.null(names(datasets)) ||
      any(names(datasets) == "") ||
      anyDuplicated(names(datasets))) {
    stop("datasets must have unique nonempty names.")
  }

  normalized_time <- lapply(datasets, function(dat) {
    if (is.null(dat$x) || is.null(dat$y) || is.null(dat$sd)) {
      stop("Every dataset must contain x, y, and sd.")
    }
    x <- as.numeric(dat$x)
    y <- as.numeric(dat$y)
    se <- as.numeric(dat$sd)
    if (length(x) < 2L ||
        length(y) != length(x) ||
        length(se) != length(x) ||
        any(!is.finite(x)) ||
        any(!is.finite(y)) ||
        any(!is.finite(se)) ||
        any(se <= 0) ||
        diff(range(x)) <= 0) {
      stop("A dataset contains invalid time, response, or standard errors.")
    }
    (x - min(x)) / pred_step
  })

  likelihood <- vapply(grid, function(predstep_sd) {
    vapply(seq_along(datasets), function(index) {
      dat <- datasets[[index]]
      log_marginal_common_intercept(
        y = dat$y,
        se = dat$sd,
        x = normalized_time[[index]],
        slope_sd = predstep_sd
      )
    }, numeric(1))
  }, numeric(length(datasets)))
  dimnames(likelihood) <- list(names(datasets), as.character(grid))
  if (any(!is.finite(likelihood))) {
    stop("The linear-mixture marginal log likelihood is not finite.")
  }
  likelihood
}

fit_linear_mixture_fash_from_log_likelihood <- function(
    L_matrix,
    grid = default_revision_grid(),
    pred_step = 1,
    penalty = 10) {
  if (!requireNamespace("fashr", quietly = TRUE)) {
    stop("The fashr package is required to fit the linear mixture.")
  }
  grid <- validate_linear_mixture_grid(grid)
  if (length(pred_step) != 1L ||
      !is.finite(pred_step) ||
      pred_step <= 0) {
    stop("pred_step must be one positive finite value.")
  }
  if (length(penalty) != 1L ||
      !is.finite(penalty) ||
      penalty < 1 ||
      penalty != as.integer(penalty)) {
    stop("penalty must be one positive integer.")
  }
  penalty <- as.integer(penalty)
  if (!is.matrix(L_matrix) ||
      !is.numeric(L_matrix) ||
      nrow(L_matrix) < 1L ||
      ncol(L_matrix) != length(grid) ||
      any(!is.finite(L_matrix)) ||
      is.null(rownames(L_matrix)) ||
      anyNA(rownames(L_matrix)) ||
      any(rownames(L_matrix) == "") ||
      anyDuplicated(rownames(L_matrix)) ||
      !identical(colnames(L_matrix), as.character(grid))) {
    stop(
      "L_matrix must be a finite named matrix with one column per grid value."
    )
  }
  likelihood_scale <- exp(L_matrix)
  if (any(!is.finite(likelihood_scale)) ||
      any(likelihood_scale[, 1] <= 0) ||
      any(rowSums(likelihood_scale) <= 0)) {
    stop(
      "The likelihood scale is unsafe for fashr empirical-Bayes fitting."
    )
  }
  rm(likelihood_scale)
  eb_result <- fashr::fash_eb_est(
    L_matrix = L_matrix,
    penalty = penalty,
    grid = grid
  )
  rownames(eb_result$posterior_weight) <- rownames(L_matrix)
  if (nrow(eb_result$prior_weight) < 1L ||
      eb_result$prior_weight$psd[1] != 0 ||
      ncol(eb_result$posterior_weight) < 1L ||
      colnames(eb_result$posterior_weight)[1] != "0") {
    stop(
      "The null-favoring linear-mixture fit did not retain the exact null."
    )
  }
  lfdr <- eb_result$posterior_weight[, 1]
  names(lfdr) <- rownames(L_matrix)

  fit <- structure(list(
    prior_weights = eb_result$prior_weight,
    posterior_weights = eb_result$posterior_weight,
    psd_grid = grid,
    lfdr = lfdr,
    settings = list(
      prior_family = "finite_mixture_gaussian_linear_slope",
      scale_definition = "sd_linear_departure_at_pred_step",
      time_origin = "per_dataset_minimum",
      pred_step = pred_step,
      penalty = penalty
    ),
    L_matrix = L_matrix,
    eb_result = eb_result,
    bf_adjusted = FALSE
  ), class = "linear_mixture_fash")
  validate_linear_mixture_fash(fit)
  fit
}

fit_linear_mixture_fash <- function(datasets,
                                    grid = default_revision_grid(),
                                    pred_step = 1,
                                    penalty = 10) {
  L_matrix <- compute_linear_mixture_log_likelihood(
    datasets = datasets,
    grid = grid,
    pred_step = pred_step
  )
  fit_linear_mixture_fash_from_log_likelihood(
    L_matrix = L_matrix,
    grid = grid,
    pred_step = pred_step,
    penalty = penalty
  )
}

validate_linear_mixture_fash <- function(fit,
                                         expected_grid = NULL,
                                         expected_pred_step = NULL,
                                         expected_penalty = NULL) {
  if (!inherits(fit, "linear_mixture_fash")) {
    stop("fit must inherit from linear_mixture_fash.")
  }
  required_fields <- c(
    "prior_weights", "posterior_weights", "psd_grid", "lfdr",
    "settings", "L_matrix", "bf_adjusted"
  )
  if (!all(required_fields %in% names(fit))) {
    stop("The linear-mixture fit is missing required fields.")
  }

  grid <- validate_linear_mixture_grid(fit$psd_grid)
  if (!is.null(expected_grid) &&
      !isTRUE(all.equal(
        grid,
        validate_linear_mixture_grid(expected_grid),
        tolerance = 0
      ))) {
    stop("The linear-mixture grid does not match the expected grid.")
  }
  if (!is.matrix(fit$L_matrix) ||
      nrow(fit$L_matrix) < 1L ||
      ncol(fit$L_matrix) != length(grid) ||
      any(!is.finite(fit$L_matrix)) ||
      is.null(rownames(fit$L_matrix)) ||
      anyNA(rownames(fit$L_matrix)) ||
      any(rownames(fit$L_matrix) == "") ||
      anyDuplicated(rownames(fit$L_matrix)) ||
      !identical(colnames(fit$L_matrix), as.character(grid))) {
    stop("The linear-mixture likelihood matrix is invalid.")
  }
  likelihood_scale <- exp(fit$L_matrix)
  if (any(!is.finite(likelihood_scale)) ||
      any(likelihood_scale[, 1] <= 0) ||
      any(rowSums(likelihood_scale) <= 0)) {
    stop("The linear-mixture likelihood scale is numerically unsafe.")
  }
  rm(likelihood_scale)

  expanded_prior <- expand_grid_prior_weights(fit$prior_weights, grid)
  posterior <- fit$posterior_weights
  retained_scales <- as.character(fit$prior_weights$psd)
  if (fit$prior_weights$psd[1] != 0 ||
      retained_scales[1] != "0" ||
      !is.matrix(posterior) ||
      nrow(posterior) != nrow(fit$L_matrix) ||
      ncol(posterior) != length(retained_scales) ||
      !identical(rownames(posterior), rownames(fit$L_matrix)) ||
      !identical(colnames(posterior), retained_scales) ||
      any(!is.finite(posterior)) ||
      any(posterior < -1e-12) ||
      any(abs(rowSums(posterior) - 1) > 1e-8)) {
    stop("The linear-mixture posterior weights are invalid.")
  }

  expected_lfdr <- posterior[, 1]
  if (length(fit$lfdr) != nrow(posterior) ||
      any(!is.finite(fit$lfdr)) ||
      any(fit$lfdr < -1e-12 | fit$lfdr > 1 + 1e-12) ||
      !isTRUE(all.equal(
        as.numeric(fit$lfdr),
        as.numeric(expected_lfdr),
        tolerance = 1e-10
      ))) {
    stop("The linear-mixture lfdr is inconsistent with the exact null.")
  }

  settings <- fit$settings
  if (!is.list(settings) ||
      !identical(
        settings$prior_family,
        "finite_mixture_gaussian_linear_slope"
      ) ||
      !identical(
        settings$scale_definition,
        "sd_linear_departure_at_pred_step"
      ) ||
      length(settings$pred_step) != 1L ||
      !is.finite(settings$pred_step) ||
      settings$pred_step <= 0 ||
      length(settings$penalty) != 1L ||
      !is.finite(settings$penalty) ||
      settings$penalty < 1 ||
      settings$penalty != as.integer(settings$penalty)) {
    stop("The linear-mixture settings are invalid.")
  }
  if (!is.null(expected_pred_step) &&
      !isTRUE(all.equal(
        settings$pred_step,
        expected_pred_step,
        tolerance = 0
      ))) {
    stop("pred_step does not match the expected value.")
  }
  if (!is.null(expected_penalty) &&
      !identical(as.integer(settings$penalty), as.integer(expected_penalty))) {
    stop("penalty does not match the expected value.")
  }
  if (!is.logical(fit$bf_adjusted) ||
      length(fit$bf_adjusted) != 1L ||
      is.na(fit$bf_adjusted)) {
    stop("The BF-adjustment indicator is invalid.")
  }
  if (!isTRUE(fit$bf_adjusted) && !is.null(fit$BF)) {
    stop("An unadjusted linear-mixture fit must not contain Bayes factors.")
  }
  if (isTRUE(fit$bf_adjusted) &&
      (is.null(fit$BF) ||
        length(fit$BF) != nrow(fit$L_matrix) ||
        any(!is.finite(fit$BF)) ||
        any(fit$BF < 0))) {
    stop("The BF-adjusted fit does not contain valid Bayes factors.")
  }
  if (abs(sum(expanded_prior) - 1) > 1e-6) {
    stop("The expanded prior weights do not sum to one.")
  }
  invisible(TRUE)
}

validate_compact_linear_mixture_fash <- function(
    fit,
    expected_grid = NULL,
    expected_pred_step = NULL,
    expected_penalty = NULL) {
  if (!inherits(fit, "compact_linear_mixture_fash")) {
    stop("fit must inherit from compact_linear_mixture_fash.")
  }
  required_fields <- c(
    "unit_ids", "prior_weights", "posterior_weights", "psd_grid", "lfdr",
    "settings", "bf_adjusted"
  )
  if (!all(required_fields %in% names(fit))) {
    stop("The compact linear-mixture fit is missing required fields.")
  }

  grid <- validate_linear_mixture_grid(fit$psd_grid)
  if (!is.null(expected_grid) &&
      !isTRUE(all.equal(
        grid,
        validate_linear_mixture_grid(expected_grid),
        tolerance = 0
      ))) {
    stop("The compact linear-mixture grid does not match the expected grid.")
  }
  unit_ids <- as.character(fit$unit_ids)
  if (length(unit_ids) < 1L ||
      anyNA(unit_ids) ||
      any(unit_ids == "") ||
      anyDuplicated(unit_ids)) {
    stop("The compact linear-mixture unit IDs are invalid.")
  }

  expanded_prior <- expand_grid_prior_weights(fit$prior_weights, grid)
  retained_scales <- as.character(fit$prior_weights$psd)
  posterior <- fit$posterior_weights
  if (fit$prior_weights$psd[1] != 0 ||
      retained_scales[1] != "0" ||
      !is.matrix(posterior) ||
      nrow(posterior) != length(unit_ids) ||
      ncol(posterior) != length(retained_scales) ||
      !identical(rownames(posterior), unit_ids) ||
      !identical(colnames(posterior), retained_scales) ||
      any(!is.finite(posterior)) ||
      any(posterior < -1e-12) ||
      any(abs(rowSums(posterior) - 1) > 1e-8)) {
    stop("The compact linear-mixture posterior weights are invalid.")
  }
  if (length(fit$lfdr) != length(unit_ids) ||
      !identical(names(fit$lfdr), unit_ids) ||
      any(!is.finite(fit$lfdr)) ||
      any(fit$lfdr < -1e-12 | fit$lfdr > 1 + 1e-12) ||
      !isTRUE(all.equal(
        as.numeric(fit$lfdr),
        as.numeric(posterior[, 1]),
        tolerance = 1e-10
      ))) {
    stop("The compact linear-mixture lfdr is inconsistent with the null.")
  }

  settings <- fit$settings
  if (!is.list(settings) ||
      !identical(
        settings$prior_family,
        "finite_mixture_gaussian_linear_slope"
      ) ||
      !identical(
        settings$scale_definition,
        "sd_linear_departure_at_pred_step"
      ) ||
      length(settings$pred_step) != 1L ||
      !is.finite(settings$pred_step) ||
      settings$pred_step <= 0 ||
      length(settings$penalty) != 1L ||
      !is.finite(settings$penalty) ||
      settings$penalty < 1 ||
      settings$penalty != as.integer(settings$penalty)) {
    stop("The compact linear-mixture settings are invalid.")
  }
  if (!is.null(expected_pred_step) &&
      !isTRUE(all.equal(
        settings$pred_step,
        expected_pred_step,
        tolerance = 0
      ))) {
    stop("pred_step does not match the expected value.")
  }
  if (!is.null(expected_penalty) &&
      !identical(as.integer(settings$penalty), as.integer(expected_penalty))) {
    stop("penalty does not match the expected value.")
  }
  if (!is.logical(fit$bf_adjusted) ||
      length(fit$bf_adjusted) != 1L ||
      is.na(fit$bf_adjusted)) {
    stop("The compact BF-adjustment indicator is invalid.")
  }
  if (!isTRUE(fit$bf_adjusted) && !is.null(fit$BF)) {
    stop("An unadjusted compact fit must not contain Bayes factors.")
  }
  if (isTRUE(fit$bf_adjusted) &&
      (is.null(fit$BF) ||
        length(fit$BF) != length(unit_ids) ||
        any(!is.finite(fit$BF)) ||
        any(fit$BF < 0))) {
    stop("The compact BF-adjusted fit has invalid Bayes factors.")
  }
  if (abs(sum(expanded_prior) - 1) > 1e-6) {
    stop("The compact expanded prior weights do not sum to one.")
  }
  invisible(TRUE)
}

compact_linear_mixture_fash <- function(fit) {
  validate_linear_mixture_fash(fit)
  compact <- list(
    unit_ids = rownames(fit$posterior_weights),
    prior_weights = fit$prior_weights,
    posterior_weights = fit$posterior_weights,
    psd_grid = fit$psd_grid,
    lfdr = fit$lfdr,
    settings = fit$settings,
    bf_adjusted = fit$bf_adjusted
  )
  if (!is.null(fit$BF)) {
    compact$BF <- fit$BF
  }
  class(compact) <- "compact_linear_mixture_fash"
  validate_compact_linear_mixture_fash(compact)
  compact
}

BF_update_linear_mixture_fash <- function(fit) {
  validate_linear_mixture_fash(fit)
  updated <- fashr::BF_update(fit, plot = FALSE)
  updated$bf_adjusted <- TRUE
  class(updated) <- "linear_mixture_fash"
  validate_linear_mixture_fash(updated)
  updated
}

extract_linear_mixture_prior_table <- function(fit,
                                               seed,
                                               fit_label) {
  validate_linear_mixture_fash(fit)
  if (length(seed) != 1L || !is.finite(seed)) {
    stop("seed must be one finite value.")
  }
  if (length(fit_label) != 1L ||
      is.na(fit_label) ||
      !nzchar(fit_label)) {
    stop("fit_label must be one nonempty value.")
  }
  grid <- fit$psd_grid
  weight <- expand_grid_prior_weights(fit$prior_weights, grid)
  data.frame(
    seed = seed,
    fit = fit_label,
    predstep_sd = grid,
    prior_weight = unname(weight),
    is_null = grid == 0,
    active = unname(weight) > 0,
    stringsAsFactors = FALSE
  )
}

summarize_linear_mixture_prior_fit <- function(fit,
                                               seed,
                                               fit_label) {
  prior <- extract_linear_mixture_prior_table(
    fit = fit,
    seed = seed,
    fit_label = fit_label
  )
  null_weight <- prior$prior_weight[prior$is_null]
  alternative <- prior[!prior$is_null, , drop = FALSE]
  alternative_weight <- sum(alternative$prior_weight)
  alternative_rms_predstep_sd <- if (alternative_weight > 0) {
    sqrt(sum(
      alternative$prior_weight / alternative_weight *
        alternative$predstep_sd^2
    ))
  } else {
    NA_real_
  }
  data.frame(
    seed = seed,
    fit = fit_label,
    estimated_pi0 = null_weight,
    active_nonnull_components = sum(alternative$active),
    alternative_rms_predstep_sd = alternative_rms_predstep_sd,
    stringsAsFactors = FALSE
  )
}

evaluate_simplified_fash_fit <- function(fit,
                                         unit_info,
                                         alpha = 0.05,
                                         method = "FASH-linear-Raw") {
  evaluate_lfdr_method(
    lfdr = fit$lfdr,
    unit_info = unit_info,
    method = method,
    target = "dynamic",
    alpha = alpha
  )
}

evaluate_lfdr_method <- function(lfdr,
                                 unit_info,
                                 method,
                                 target = c("dynamic", "nonlinear"),
                                 alpha = 0.05) {
  target <- match.arg(target)
  stopifnot(length(lfdr) == nrow(unit_info))

  ord <- order(lfdr)
  cumulative_fdr <- cumsum(lfdr[ord]) / seq_along(lfdr)
  selected_ordered <- cumulative_fdr <= alpha
  selected <- rep(FALSE, length(lfdr))
  selected[ord[selected_ordered]] <- TRUE

  data.frame(
    unit_index = unit_info$unit_index,
    unit_id = unit_info$unit_id,
    effect_class = unit_info$effect_class,
    scenario = unit_info$scenario,
    target = target,
    method = method,
    score_type = "lfdr",
    score = lfdr,
    adjusted_score = cumulative_fdr[match(seq_along(lfdr), ord)],
    alpha = alpha,
    selected = selected,
    true_null = is_null_class(unit_info$effect_class, target),
    stringsAsFactors = FALSE
  )
}

weighted_lm_pvalue <- function(dat, target = c("dynamic", "nonlinear")) {
  target <- match.arg(target)
  dat$x_centered <- as.numeric(scale(dat$x, center = TRUE, scale = TRUE))
  weights <- 1 / pmax(dat$sd, .Machine$double.eps)^2

  if (target == "dynamic") {
    fit <- lm(y ~ x_centered, data = dat, weights = weights)
    coef_table <- summary(fit)$coefficients
    return(coef_table["x_centered", "Pr(>|t|)"])
  }

  fit <- lm(y ~ x_centered + I(x_centered^2), data = dat, weights = weights)
  coef_table <- summary(fit)$coefficients
  coef_table["I(x_centered^2)", "Pr(>|t|)"]
}

evaluate_parametric_method <- function(datasets,
                                       unit_info,
                                       method = c("linear_summary", "quadratic_summary"),
                                       target = c("dynamic", "nonlinear"),
                                       alpha = 0.05,
                                       adjust_method = "BH") {
  method <- match.arg(method)
  target <- match.arg(target)

  pvalues <- vapply(
    datasets,
    weighted_lm_pvalue,
    numeric(1),
    target = target
  )
  qvalues <- p.adjust(pvalues, method = adjust_method)

  data.frame(
    unit_index = unit_info$unit_index,
    unit_id = unit_info$unit_id,
    effect_class = unit_info$effect_class,
    scenario = unit_info$scenario,
    target = target,
    method = method,
    score_type = "pvalue",
    score = pvalues,
    adjusted_score = qvalues,
    alpha = alpha,
    selected = qvalues <= alpha,
    true_null = is_null_class(unit_info$effect_class, target),
    stringsAsFactors = FALSE
  )
}

evaluate_fash_fits <- function(fits, unit_info, alpha = 0.05) {
  out <- list()
  for (nm in names(fits)) {
    target <- if (grepl("iwp1", nm)) "dynamic" else "nonlinear"
    out[[nm]] <- evaluate_lfdr_method(
      lfdr = get_fash_lfdr(fits[[nm]]),
      unit_info = unit_info,
      method = nm,
      target = target,
      alpha = alpha
    )
  }
  do.call(rbind, out)
}

revision_method_order <- function() {
  c(
    "FASH-IWP1-Raw",
    "FASH-IWP1-BF",
    "FASH-linear-Raw",
    "FASH-linear-BF",
    "Direct-linear-LRT",
    "Direct-linear-LRT-eFDR",
    "Direct-linear-LRT-eFDR-true-pi0",
    "Direct-quadratic-LRT",
    "Direct-quadratic-LRT-eFDR",
    "Direct-quadratic-LRT-eFDR-true-pi0",
    "fash_iwp1_raw",
    "fash_iwp1_bf",
    "fash_iwp2_raw",
    "fash_iwp2_bf",
    "linear_summary",
    "quadratic_summary"
  )
}

revision_final_summary_methods <- function() {
  c(
    "FASH-IWP1-Raw",
    "FASH-IWP1-BF",
    "FASH-linear-Raw",
    "FASH-linear-BF",
    "Direct-linear-LRT-eFDR",
    "Direct-linear-LRT-eFDR-true-pi0",
    "Direct-quadratic-LRT-eFDR",
    "Direct-quadratic-LRT-eFDR-true-pi0"
  )
}

rank_revision_methods <- function(method) {
  method_rank <- match(method, revision_method_order())
  method_rank[is.na(method_rank)] <- length(revision_method_order()) + 1
  method_rank
}

summarize_method_results <- function(result_table) {
  split_results <- split(
    result_table,
    list(result_table$scenario, result_table$target, result_table$method),
    drop = TRUE
  )

  summaries <- lapply(split_results, function(x) {
    discoveries <- sum(x$selected)
    false_discoveries <- sum(x$selected & x$true_null)
    true_alternatives <- sum(!x$true_null)
    true_positives <- sum(x$selected & !x$true_null)
    empirical_fdr <- if (discoveries == 0) 0 else false_discoveries / discoveries
    power <- if (true_alternatives == 0) NA_real_ else true_positives / true_alternatives

    data.frame(
      scenario = x$scenario[1],
      target = x$target[1],
      method = x$method[1],
      alpha = x$alpha[1],
      n_units = nrow(x),
      n_discoveries = discoveries,
      false_discoveries = false_discoveries,
      true_positives = true_positives,
      empirical_fdr = empirical_fdr,
      power = power,
      stringsAsFactors = FALSE
    )
  })

  summary_table <- do.call(rbind, summaries)
  summary_table$method_rank <- rank_revision_methods(summary_table$method)
  summary_table <- summary_table[order(
    summary_table$scenario,
    summary_table$target,
    summary_table$method_rank,
    summary_table$method
  ), ]
  summary_table$method_rank <- NULL
  rownames(summary_table) <- NULL
  summary_table
}

compute_alpha_curve <- function(result_table,
                                alpha_grid = seq(0, 0.20, by = 0.005)) {
  required_columns <- c(
    "scenario",
    "target",
    "method",
    "adjusted_score",
    "true_null"
  )
  missing_columns <- setdiff(required_columns, colnames(result_table))
  if (length(missing_columns) > 0) {
    stop("result_table is missing required columns: ", paste(missing_columns, collapse = ", "))
  }

  alpha_grid <- sort(unique(alpha_grid))
  if (any(is.na(alpha_grid)) || any(alpha_grid < 0) || any(alpha_grid > 1)) {
    stop("alpha_grid must contain finite values between 0 and 1.")
  }

  split_results <- split(
    result_table,
    list(result_table$scenario, result_table$target, result_table$method),
    drop = TRUE
  )

  curve_rows <- lapply(split_results, function(x) {
    do.call(rbind, lapply(alpha_grid, function(alpha) {
      selected <- x$adjusted_score <= alpha
      discoveries <- sum(selected, na.rm = TRUE)
      false_discoveries <- sum(selected & x$true_null, na.rm = TRUE)
      true_alternatives <- sum(!x$true_null, na.rm = TRUE)
      true_positives <- sum(selected & !x$true_null, na.rm = TRUE)
      empirical_fdr <- if (discoveries == 0) 0 else false_discoveries / discoveries
      power <- if (true_alternatives == 0) NA_real_ else true_positives / true_alternatives

      data.frame(
        scenario = x$scenario[1],
        target = x$target[1],
        method = x$method[1],
        alpha = alpha,
        n_units = nrow(x),
        n_discoveries = discoveries,
        false_discoveries = false_discoveries,
        true_positives = true_positives,
        empirical_fdr = empirical_fdr,
        power = power,
        stringsAsFactors = FALSE
      )
    }))
  })

  alpha_curve <- do.call(rbind, curve_rows)
  alpha_curve$method_rank <- rank_revision_methods(alpha_curve$method)
  alpha_curve <- alpha_curve[order(
    alpha_curve$scenario,
    alpha_curve$target,
    alpha_curve$method_rank,
    alpha_curve$method,
    alpha_curve$alpha
  ), ]
  alpha_curve$method_rank <- NULL
  rownames(alpha_curve) <- NULL
  alpha_curve
}

constant_component_prior_weight <- function(fit) {
  prior_weights <- fit$prior_weights
  if (is.null(prior_weights)) {
    stop("fit does not contain prior_weights.")
  }
  if ("psd" %in% names(prior_weights)) {
    weight <- prior_weights$prior_weight[prior_weights$psd == 0]
  } else if ("component" %in% names(prior_weights)) {
    weight <- prior_weights$prior_weight[prior_weights$component == "constant"]
  } else {
    stop("Could not identify the constant component in prior_weights.")
  }
  if (length(weight) != 1 || !is.finite(weight)) {
    stop("Expected exactly one finite constant-component prior weight.")
  }
  unname(weight)
}

summarize_mc_values <- function(values, confidence_level = 0.95) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]
  n <- length(values)
  if (n == 0) {
    return(c(mean = NA_real_, sd = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  mean_value <- mean(values)
  if (n == 1) {
    return(c(mean = mean_value, sd = NA_real_, se = NA_real_, lower = NA_real_, upper = NA_real_))
  }
  sd_value <- stats::sd(values)
  se_value <- sd_value / sqrt(n)
  critical_value <- stats::qt((1 + confidence_level) / 2, df = n - 1)
  c(
    mean = mean_value,
    sd = sd_value,
    se = se_value,
    lower = mean_value - critical_value * se_value,
    upper = mean_value + critical_value * se_value
  )
}

summarize_mc_alpha_curves <- function(alpha_rows,
                                      confidence_level = 0.95) {
  required_columns <- c(
    "seed", "scenario", "target", "method", "alpha", "n_discoveries",
    "false_discoveries", "true_positives", "empirical_fdr", "power"
  )
  missing_columns <- setdiff(required_columns, names(alpha_rows))
  if (length(missing_columns) > 0) {
    stop("alpha_rows is missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  if (anyDuplicated(alpha_rows[c("seed", "scenario", "target", "method", "alpha")])) {
    stop("alpha_rows must contain at most one row per seed, method, and alpha.")
  }

  groups <- split(
    alpha_rows,
    list(alpha_rows$scenario, alpha_rows$target, alpha_rows$method, alpha_rows$alpha),
    drop = TRUE
  )
  summaries <- lapply(groups, function(x) {
    metrics <- list(
      discoveries = summarize_mc_values(x$n_discoveries, confidence_level),
      false_discoveries = summarize_mc_values(x$false_discoveries, confidence_level),
      true_positives = summarize_mc_values(x$true_positives, confidence_level),
      power = summarize_mc_values(x$power, confidence_level),
      fdr = summarize_mc_values(x$empirical_fdr, confidence_level)
    )
    data.frame(
      scenario = x$scenario[1],
      target = x$target[1],
      method = x$method[1],
      alpha = x$alpha[1],
      n_replications = length(unique(x$seed)),
      mean_discoveries = metrics$discoveries[["mean"]],
      mean_false_discoveries = metrics$false_discoveries[["mean"]],
      mean_true_positives = metrics$true_positives[["mean"]],
      mean_power = metrics$power[["mean"]],
      power_sd = metrics$power[["sd"]],
      power_mc_se = metrics$power[["se"]],
      power_ci_lower = pmax(0, metrics$power[["lower"]]),
      power_ci_upper = pmin(1, metrics$power[["upper"]]),
      mean_fdr = metrics$fdr[["mean"]],
      fdr_sd = metrics$fdr[["sd"]],
      fdr_mc_se = metrics$fdr[["se"]],
      fdr_ci_lower = pmax(0, metrics$fdr[["lower"]]),
      fdr_ci_upper = pmin(1, metrics$fdr[["upper"]]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, summaries)
  out$method_rank <- rank_revision_methods(out$method)
  out <- out[order(out$scenario, out$target, out$method_rank, out$method, out$alpha), ]
  out$method_rank <- NULL
  rownames(out) <- NULL
  out
}

summarize_mc_functional_alpha_curves <- function(alpha_rows,
                                                 confidence_level = 0.95) {
  required_columns <- c(
    "seed", "scenario", "target", "method", "alpha", "dynamic_discoveries",
    "n_discoveries", "false_discoveries", "conditional_discoveries",
    "conditional_false_discoveries", "first_stage_null_calls", "true_positives",
    "estimated_fsr", "empirical_fsr", "conditional_empirical_fsr", "power"
  )
  missing_columns <- setdiff(required_columns, names(alpha_rows))
  if (length(missing_columns) > 0) {
    stop("alpha_rows is missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  if (anyDuplicated(alpha_rows[c("seed", "scenario", "target", "method", "alpha")])) {
    stop("alpha_rows must contain at most one row per seed, target, method, and alpha.")
  }

  groups <- split(
    alpha_rows,
    list(alpha_rows$scenario, alpha_rows$target, alpha_rows$method, alpha_rows$alpha),
    drop = TRUE
  )
  summaries <- lapply(groups, function(x) {
    metrics <- list(
      dynamic_discoveries = summarize_mc_values(x$dynamic_discoveries, confidence_level),
      discoveries = summarize_mc_values(x$n_discoveries, confidence_level),
      false_discoveries = summarize_mc_values(x$false_discoveries, confidence_level),
      conditional_discoveries = summarize_mc_values(
        x$conditional_discoveries,
        confidence_level
      ),
      conditional_false_discoveries = summarize_mc_values(
        x$conditional_false_discoveries,
        confidence_level
      ),
      first_stage_null_calls = summarize_mc_values(
        x$first_stage_null_calls,
        confidence_level
      ),
      true_positives = summarize_mc_values(x$true_positives, confidence_level),
      power = summarize_mc_values(x$power, confidence_level),
      estimated_fsr = summarize_mc_values(x$estimated_fsr, confidence_level),
      empirical_fsr = summarize_mc_values(x$empirical_fsr, confidence_level),
      conditional_empirical_fsr = summarize_mc_values(
        x$conditional_empirical_fsr,
        confidence_level
      )
    )
    data.frame(
      scenario = x$scenario[1],
      target = x$target[1],
      method = x$method[1],
      alpha = x$alpha[1],
      n_replications = length(unique(x$seed)),
      mean_dynamic_discoveries = metrics$dynamic_discoveries[["mean"]],
      mean_discoveries = metrics$discoveries[["mean"]],
      mean_false_discoveries = metrics$false_discoveries[["mean"]],
      mean_conditional_discoveries = metrics$conditional_discoveries[["mean"]],
      mean_conditional_false_discoveries =
        metrics$conditional_false_discoveries[["mean"]],
      mean_first_stage_null_calls = metrics$first_stage_null_calls[["mean"]],
      mean_true_positives = metrics$true_positives[["mean"]],
      mean_power = metrics$power[["mean"]],
      power_mc_se = metrics$power[["se"]],
      power_ci_lower = pmax(0, metrics$power[["lower"]]),
      power_ci_upper = pmin(1, metrics$power[["upper"]]),
      mean_estimated_fsr = metrics$estimated_fsr[["mean"]],
      estimated_fsr_mc_se = metrics$estimated_fsr[["se"]],
      estimated_fsr_ci_lower = pmax(0, metrics$estimated_fsr[["lower"]]),
      estimated_fsr_ci_upper = pmin(1, metrics$estimated_fsr[["upper"]]),
      mean_empirical_fsr = metrics$empirical_fsr[["mean"]],
      empirical_fsr_mc_se = metrics$empirical_fsr[["se"]],
      empirical_fsr_ci_lower = pmax(0, metrics$empirical_fsr[["lower"]]),
      empirical_fsr_ci_upper = pmin(1, metrics$empirical_fsr[["upper"]]),
      mean_conditional_empirical_fsr =
        metrics$conditional_empirical_fsr[["mean"]],
      conditional_empirical_fsr_mc_se =
        metrics$conditional_empirical_fsr[["se"]],
      conditional_empirical_fsr_ci_lower = pmax(
        0,
        metrics$conditional_empirical_fsr[["lower"]]
      ),
      conditional_empirical_fsr_ci_upper = pmin(
        1,
        metrics$conditional_empirical_fsr[["upper"]]
      ),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, summaries)
  out$method_rank <- rank_revision_methods(out$method)
  out <- out[order(out$scenario, out$target, out$method_rank, out$method, out$alpha), ]
  out$method_rank <- NULL
  rownames(out) <- NULL
  out
}

plot_mc_functional_curve_grid <- function(mc_curve,
                                          metric = c(
                                            "power",
                                            "estimated_fsr",
                                            "empirical_fsr",
                                            "conditional_empirical_fsr"
                                          ),
                                          target_order = c("early", "middle", "late", "switch"),
                                          file = NULL,
                                          title = NULL,
                                          alpha_reference = 0.05,
                                          ncol = 2,
                                          width = 2200,
                                          height = 1600,
                                          res = 180) {
  metric <- match.arg(metric)
  metric_columns <- switch(
    metric,
    power = c("mean_power", "power_ci_lower", "power_ci_upper"),
    estimated_fsr = c(
      "mean_estimated_fsr", "estimated_fsr_ci_lower", "estimated_fsr_ci_upper"
    ),
    empirical_fsr = c(
      "mean_empirical_fsr", "empirical_fsr_ci_lower", "empirical_fsr_ci_upper"
    ),
    conditional_empirical_fsr = c(
      "mean_conditional_empirical_fsr",
      "conditional_empirical_fsr_ci_lower",
      "conditional_empirical_fsr_ci_upper"
    )
  )
  required_columns <- c("target", "method", "alpha", metric_columns)
  missing_columns <- setdiff(required_columns, names(mc_curve))
  if (length(missing_columns) > 0) {
    stop("mc_curve is missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  target_order <- target_order[target_order %in% unique(mc_curve$target)]
  if (length(target_order) == 0) {
    stop("mc_curve contains no requested targets.")
  }

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    png(file, width = width, height = height, res = res)
    on.exit(dev.off(), add = TRUE)
  }
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  nrow <- ceiling(length(target_order) / ncol)
  par(mfrow = c(nrow, ncol), mar = c(4.1, 4.2, 3.2, 1.0), oma = c(0, 0, 3.0, 0))

  styles <- revision_method_styles(unique(mc_curve$method), style_profile = "combined")
  x_limits <- range(mc_curve$alpha, na.rm = TRUE)
  value_column <- metric_columns[1]
  lower_column <- metric_columns[2]
  upper_column <- metric_columns[3]
  y_limits <- if (metric == "power") {
    c(0, 1)
  } else {
    upper <- max(c(x_limits[2], mc_curve[[upper_column]]), na.rm = TRUE)
    c(0, min(1, upper * 1.08))
  }
  y_label <- switch(
    metric,
    power = "Mean power across replications",
    estimated_fsr = "Mean estimated FSR",
    empirical_fsr = "Mean empirical FSR across replications",
    conditional_empirical_fsr = "Mean conditional empirical FSR"
  )
  if (is.null(title)) {
    title <- switch(
      metric,
      power = "Functional-testing power",
      estimated_fsr = "Estimated FSR from posterior cFSR calls",
      empirical_fsr = "End-to-end false-call proportion from simulation truth",
      conditional_empirical_fsr =
        "Conditional functional FSR among true dynamic eQTLs"
    )
  }

  for (target in target_order) {
    panel_curve <- mc_curve[mc_curve$target == target, , drop = FALSE]
    plot(
      NA,
      xlim = x_limits,
      ylim = y_limits,
      xlab = "Nominal level alpha",
      ylab = y_label,
      main = ""
    )
    mtext(
      tools::toTitleCase(target),
      side = 3,
      line = 1,
      font = 2,
      cex = 1.05
    )
    grid(col = "gray90")
    if (metric != "power") {
      abline(a = 0, b = 1, col = "gray35", lty = 3)
    }
    if (!is.null(alpha_reference)) {
      abline(v = alpha_reference, col = "gray55", lty = 3)
    }
    for (i in seq_along(styles$methods)) {
      method_curve <- panel_curve[panel_curve$method == styles$methods[i], , drop = FALSE]
      method_curve <- method_curve[order(method_curve$alpha), , drop = FALSE]
      if (nrow(method_curve) == 0) next
      lower <- pmax(y_limits[1], method_curve[[lower_column]])
      upper <- pmin(y_limits[2], method_curve[[upper_column]])
      if (all(is.finite(c(lower, upper)))) {
        polygon(
          x = c(method_curve$alpha, rev(method_curve$alpha)),
          y = c(lower, rev(upper)),
          col = grDevices::adjustcolor(styles$col[i], alpha.f = 0.14),
          border = NA
        )
      }
      lines(
        method_curve$alpha,
        method_curve[[value_column]],
        col = styles$col[i],
        lty = styles$lty[i],
        lwd = 2.4
      )
    }
    legend(
      "topleft",
      legend = c(styles$methods, if (metric != "power") "Nominal alpha"),
      col = c(styles$col, if (metric != "power") "gray35"),
      lty = c(styles$lty, if (metric != "power") 3),
      lwd = c(rep(2.4, length(styles$methods)), if (metric != "power") 1.8),
      bty = "n",
      cex = 0.78
    )
  }
  mtext(title, outer = TRUE, cex = 1.2, font = 2)
  invisible(mc_curve)
}

summarize_mc_pi0 <- function(pi0_rows, confidence_level = 0.95) {
  required_columns <- c("seed", "method", "fit", "estimated_pi0")
  missing_columns <- setdiff(required_columns, names(pi0_rows))
  if (length(missing_columns) > 0) {
    stop("pi0_rows is missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  if (anyDuplicated(pi0_rows[c("seed", "method", "fit")])) {
    stop("pi0_rows must contain at most one row per seed, method, and fit.")
  }

  groups <- split(pi0_rows, list(pi0_rows$method, pi0_rows$fit), drop = TRUE)
  summaries <- lapply(groups, function(x) {
    metrics <- summarize_mc_values(x$estimated_pi0, confidence_level)
    data.frame(
      method = x$method[1],
      fit = x$fit[1],
      n_replications = length(unique(x$seed)),
      mean_estimated_pi0 = metrics[["mean"]],
      pi0_sd = metrics[["sd"]],
      pi0_mc_se = metrics[["se"]],
      pi0_ci_lower = pmax(0, metrics[["lower"]]),
      pi0_ci_upper = pmin(1, metrics[["upper"]]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, summaries)
  rownames(out) <- NULL
  out
}

revision_method_styles <- function(methods,
                                   style_profile = c("default", "combined")) {
  style_profile <- match.arg(style_profile)
  ordered_methods <- revision_method_order()
  ordered_methods <- ordered_methods[ordered_methods %in% methods]
  remaining_methods <- sort(setdiff(methods, ordered_methods))
  methods <- c(ordered_methods, remaining_methods)

  colors <- c(
    "FASH-IWP1-Raw" = "#1B9E77",
    "FASH-IWP1-BF" = "#0072B2",
    "FASH-linear-Raw" = "#D55E00",
    "FASH-linear-BF" = "#CC79A7",
    "Direct-linear-LRT" = "#4D4D4D",
    "Direct-linear-LRT-eFDR" = "#4D4D4D",
    "Direct-linear-LRT-eFDR-true-pi0" = "#4D4D4D",
    "Direct-quadratic-LRT" = "#984EA3",
    "Direct-quadratic-LRT-eFDR" = "#984EA3",
    "Direct-quadratic-LRT-eFDR-true-pi0" = "#984EA3",
    "fash_iwp1_raw" = "#1B9E77",
    "fash_iwp1_bf" = "#0072B2",
    "fash_iwp2_raw" = "#7570B3",
    "fash_iwp2_bf" = "#E7298A",
    "linear_summary" = "#D55E00",
    "quadratic_summary" = "#CC79A7"
  )
  line_types <- c(
    "FASH-IWP1-Raw" = 1,
    "FASH-IWP1-BF" = 2,
    "FASH-linear-Raw" = 1,
    "FASH-linear-BF" = 2,
    "Direct-linear-LRT" = 3,
    "Direct-linear-LRT-eFDR" = 2,
    "Direct-linear-LRT-eFDR-true-pi0" = 5,
    "Direct-quadratic-LRT" = 4,
    "Direct-quadratic-LRT-eFDR" = 2,
    "Direct-quadratic-LRT-eFDR-true-pi0" = 5,
    "fash_iwp1_raw" = 1,
    "fash_iwp1_bf" = 2,
    "fash_iwp2_raw" = 1,
    "fash_iwp2_bf" = 2,
    "linear_summary" = 1,
    "quadratic_summary" = 2
  )
  points <- c(
    "FASH-IWP1-Raw" = 16,
    "FASH-IWP1-BF" = 17,
    "FASH-linear-Raw" = 15,
    "FASH-linear-BF" = 18,
    "Direct-linear-LRT" = 3,
    "Direct-linear-LRT-eFDR" = 8,
    "Direct-linear-LRT-eFDR-true-pi0" = 1,
    "Direct-quadratic-LRT" = 4,
    "Direct-quadratic-LRT-eFDR" = 8,
    "Direct-quadratic-LRT-eFDR-true-pi0" = 1,
    "fash_iwp1_raw" = 16,
    "fash_iwp1_bf" = 17,
    "fash_iwp2_raw" = 15,
    "fash_iwp2_bf" = 18,
    "linear_summary" = 15,
    "quadratic_summary" = 18
  )

  if (style_profile == "combined") {
    # Use colour for model family and line type for the FASH BF correction.
    colors["FASH-IWP1-BF"] <- colors["FASH-IWP1-Raw"]
    colors["FASH-linear-BF"] <- colors["FASH-linear-Raw"]
    line_types[c("FASH-IWP1-Raw", "FASH-linear-Raw")] <- 1
    line_types[c("FASH-IWP1-BF", "FASH-linear-BF")] <- 2
    line_types[c(
      "Direct-linear-LRT-eFDR-true-pi0",
      "Direct-quadratic-LRT-eFDR-true-pi0"
    )] <- 1
  }

  list(
    methods = methods,
    col = unname(colors[methods]),
    lty = unname(line_types[methods]),
    pch = unname(points[methods])
  )
}

plot_power_alpha_curves <- function(alpha_curve,
                                    file = NULL,
                                    title = "Power as a function of alpha",
                                    subtitle = NULL,
                                    alpha_reference = 0.05,
                                    legend_position = "bottomright",
                                    style_profile = c("default", "combined"),
                                    width = 1500,
                                    height = 1050,
                                    res = 170) {
  style_profile <- match.arg(style_profile)
  if (!all(c("method", "alpha", "power") %in% colnames(alpha_curve))) {
    stop("alpha_curve must contain method, alpha, and power columns.")
  }

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    png(file, width = width, height = height, res = res)
    on.exit(dev.off(), add = TRUE)
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(4.4, 4.4, 3.6, 1.2))

  styles <- revision_method_styles(
    unique(alpha_curve$method),
    style_profile = style_profile
  )
  x_limits <- range(alpha_curve$alpha, na.rm = TRUE)
  y_limits <- c(0, 1)
  plot(
    NA,
    xlim = x_limits,
    ylim = y_limits,
    xlab = "Nominal FDR level alpha",
    ylab = "Power",
    main = title
  )
  grid(col = "gray90")
  if (!is.null(alpha_reference)) {
    abline(v = alpha_reference, col = "gray55", lty = 3)
  }
  for (i in seq_along(styles$methods)) {
    method <- styles$methods[i]
    method_curve <- alpha_curve[alpha_curve$method == method, ]
    method_curve <- method_curve[order(method_curve$alpha), ]
    lines(
      method_curve$alpha,
      method_curve$power,
      col = styles$col[i],
      lty = styles$lty[i],
      lwd = 2.4
    )
  }
  legend(
    legend_position,
    legend = styles$methods,
    col = styles$col,
    lty = styles$lty,
    lwd = 2.4,
    bty = "n",
    cex = 0.9
  )
  if (!is.null(subtitle)) {
    mtext(subtitle, side = 3, line = 0.25, cex = 0.85)
  }

  invisible(alpha_curve)
}

plot_power_alpha_curve_grid <- function(alpha_curve,
                                        panel_var = "panel_label",
                                        file = NULL,
                                        title = "Power as a function of alpha",
                                        alpha_reference = 0.05,
                                        ncol = 2,
                                        width = 2200,
                                        height = 1600,
                                        res = 180) {
  if (!all(c(panel_var, "method", "alpha", "power") %in% colnames(alpha_curve))) {
    stop("alpha_curve must contain panel_var, method, alpha, and power columns.")
  }

  panels <- unique(alpha_curve[[panel_var]])
  nrow <- ceiling(length(panels) / ncol)

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    png(file, width = width, height = height, res = res)
    on.exit(dev.off(), add = TRUE)
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(nrow, ncol), mar = c(4.1, 4.2, 3.0, 1.0), oma = c(0, 0, 3.0, 0))

  styles <- revision_method_styles(unique(alpha_curve$method))
  x_limits <- range(alpha_curve$alpha, na.rm = TRUE)

  for (panel in panels) {
    panel_curve <- alpha_curve[alpha_curve[[panel_var]] == panel, ]
    plot(
      NA,
      xlim = x_limits,
      ylim = c(0, 1),
      xlab = "Nominal FDR level alpha",
      ylab = "Power",
      main = panel
    )
    grid(col = "gray90")
    if (!is.null(alpha_reference)) {
      abline(v = alpha_reference, col = "gray55", lty = 3)
    }
    for (i in seq_along(styles$methods)) {
      method <- styles$methods[i]
      method_curve <- panel_curve[panel_curve$method == method, ]
      method_curve <- method_curve[order(method_curve$alpha), ]
      lines(
        method_curve$alpha,
        method_curve$power,
        col = styles$col[i],
        lty = styles$lty[i],
        lwd = 2.4
      )
    }
    legend(
      "bottomright",
      legend = styles$methods,
      col = styles$col,
      lty = styles$lty,
      lwd = 2.4,
      bty = "n",
      cex = 0.78
    )
  }

  mtext(title, outer = TRUE, cex = 1.2, font = 2)
  invisible(alpha_curve)
}

plot_empirical_fdr_alpha_curves <- function(alpha_curve,
                                            file = NULL,
                                            title = "Empirical FDR as a function of alpha",
                                            subtitle = NULL,
                                            y_label = "Empirical FDR",
                                            alpha_reference = 0.05,
                                            legend_position = "topleft",
                                            style_profile = c("default", "combined"),
                                            width = 1500,
                                            height = 1050,
                                            res = 170) {
  style_profile <- match.arg(style_profile)
  if (!all(c("method", "alpha", "empirical_fdr") %in% colnames(alpha_curve))) {
    stop("alpha_curve must contain method, alpha, and empirical_fdr columns.")
  }

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    png(file, width = width, height = height, res = res)
    on.exit(dev.off(), add = TRUE)
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(4.4, 4.4, 3.6, 1.2))

  styles <- revision_method_styles(
    unique(alpha_curve$method),
    style_profile = style_profile
  )
  x_limits <- range(alpha_curve$alpha, na.rm = TRUE)
  y_upper <- max(c(x_limits[2], alpha_curve$empirical_fdr), na.rm = TRUE)
  y_limits <- c(0, min(1, y_upper * 1.08))
  plot(
    NA,
    xlim = x_limits,
    ylim = y_limits,
    xlab = "Nominal FDR level alpha",
    ylab = y_label,
    main = title
  )
  grid(col = "gray90")
  abline(a = 0, b = 1, col = "gray35", lty = 3)
  if (!is.null(alpha_reference)) {
    abline(v = alpha_reference, col = "gray55", lty = 3)
  }
  for (i in seq_along(styles$methods)) {
    method <- styles$methods[i]
    method_curve <- alpha_curve[alpha_curve$method == method, ]
    method_curve <- method_curve[order(method_curve$alpha), ]
    lines(
      method_curve$alpha,
      method_curve$empirical_fdr,
      col = styles$col[i],
      lty = styles$lty[i],
      lwd = 2.4
    )
  }
  legend(
    legend_position,
    legend = c(styles$methods, "Nominal alpha"),
    col = c(styles$col, "gray35"),
    lty = c(styles$lty, 3),
    lwd = c(rep(2.4, length(styles$methods)), 1.8),
    bty = "n",
    cex = 0.9
  )
  if (!is.null(subtitle)) {
    mtext(subtitle, side = 3, line = 0.25, cex = 0.85)
  }

  invisible(alpha_curve)
}

plot_mc_alpha_curves <- function(mc_curve,
                                 metric = c("power", "fdr"),
                                 file = NULL,
                                 title = NULL,
                                 subtitle = NULL,
                                 alpha_reference = 0.05,
                                 legend_position = "bottomright",
                                 style_profile = c("default", "combined"),
                                 width = 1500,
                                 height = 1050,
                                 res = 170) {
  metric <- match.arg(metric)
  style_profile <- match.arg(style_profile)
  value_column <- if (metric == "power") "mean_power" else "mean_fdr"
  lower_column <- if (metric == "power") "power_ci_lower" else "fdr_ci_lower"
  upper_column <- if (metric == "power") "power_ci_upper" else "fdr_ci_upper"
  required_columns <- c("method", "alpha", value_column, lower_column, upper_column)
  missing_columns <- setdiff(required_columns, names(mc_curve))
  if (length(missing_columns) > 0) {
    stop("mc_curve is missing required columns: ", paste(missing_columns, collapse = ", "))
  }

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    png(file, width = width, height = height, res = res)
    on.exit(dev.off(), add = TRUE)
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(4.4, 4.4, 4.8, 1.2))

  styles <- revision_method_styles(
    unique(mc_curve$method),
    style_profile = style_profile
  )
  x_limits <- range(mc_curve$alpha, na.rm = TRUE)
  if (metric == "power") {
    y_limits <- c(0, 1)
    y_label <- "Mean power across replications"
    if (is.null(title)) title <- "Monte Carlo power as a function of alpha"
  } else {
    y_upper <- max(c(x_limits[2], mc_curve[[upper_column]]), na.rm = TRUE)
    y_limits <- c(0, min(1, y_upper * 1.08))
    y_label <- "Monte Carlo estimate of FDR: E[FDP]"
    if (is.null(title)) title <- "Monte Carlo FDR as a function of alpha"
  }
  plot(
    NA,
    xlim = x_limits,
    ylim = y_limits,
    xlab = "Nominal FDR level alpha",
    ylab = y_label,
    main = title
  )
  grid(col = "gray90")
  if (metric == "fdr") {
    abline(a = 0, b = 1, col = "gray35", lty = 3)
  }
  if (!is.null(alpha_reference)) {
    abline(v = alpha_reference, col = "gray55", lty = 3)
  }
  for (i in seq_along(styles$methods)) {
    method <- styles$methods[i]
    method_curve <- mc_curve[mc_curve$method == method, ]
    method_curve <- method_curve[order(method_curve$alpha), ]
    lower <- pmax(y_limits[1], method_curve[[lower_column]])
    upper <- pmin(y_limits[2], method_curve[[upper_column]])
    if (all(is.finite(c(lower, upper)))) {
      polygon(
        x = c(method_curve$alpha, rev(method_curve$alpha)),
        y = c(lower, rev(upper)),
        col = grDevices::adjustcolor(styles$col[i], alpha.f = 0.14),
        border = NA
      )
    }
    lines(
      method_curve$alpha,
      method_curve[[value_column]],
      col = styles$col[i],
      lty = styles$lty[i],
      lwd = 2.4
    )
  }
  legend_labels <- styles$methods
  legend_col <- styles$col
  legend_lty <- styles$lty
  if (metric == "fdr") {
    legend_labels <- c(legend_labels, "Nominal alpha")
    legend_col <- c(legend_col, "gray35")
    legend_lty <- c(legend_lty, 3)
  }
  legend(
    legend_position,
    legend = legend_labels,
    col = legend_col,
    lty = legend_lty,
    lwd = c(rep(2.4, length(styles$methods)), if (metric == "fdr") 1.8 else numeric()),
    bty = "n",
    cex = 0.9
  )
  if (!is.null(subtitle)) {
    mtext(subtitle, side = 3, line = 0.25, cex = 0.85)
  }

  invisible(mc_curve)
}

plot_mc_shape_power_grid <- function(mc_curve,
                                     shape_order = c("broad", "spiky"),
                                     file = NULL,
                                     title = "Shape-stratified Monte Carlo power",
                                     subtitle = NULL,
                                     alpha_reference = 0.05,
                                     style_profile = c("default", "combined"),
                                     width = 2300,
                                     height = 1050,
                                     res = 180) {
  style_profile <- match.arg(style_profile)
  required_columns <- c(
    "shape_profile", "method", "alpha", "mean_power",
    "power_ci_lower", "power_ci_upper"
  )
  missing_columns <- setdiff(required_columns, names(mc_curve))
  if (length(missing_columns) > 0) {
    stop("mc_curve is missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  shape_order <- shape_order[shape_order %in% unique(mc_curve$shape_profile)]
  if (length(shape_order) == 0) {
    stop("mc_curve contains no requested shape profiles.")
  }

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    png(file, width = width, height = height, res = res)
    on.exit(dev.off(), add = TRUE)
  }
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(
    mfrow = c(1, length(shape_order)),
    mar = c(4.4, 4.4, 3.3, 1.0),
    oma = c(0, 0, if (is.null(subtitle)) 2.5 else 4.0, 0)
  )
  styles <- revision_method_styles(
    unique(mc_curve$method),
    style_profile = style_profile
  )
  x_limits <- range(mc_curve$alpha, na.rm = TRUE)

  for (shape in shape_order) {
    panel <- mc_curve[mc_curve$shape_profile == shape, , drop = FALSE]
    plot(
      NA,
      xlim = x_limits,
      ylim = c(0, 1),
      xlab = "Nominal FDR level alpha",
      ylab = "Mean power across replications",
      main = paste(tools::toTitleCase(shape), "effects")
    )
    grid(col = "gray90")
    if (!is.null(alpha_reference)) {
      abline(v = alpha_reference, col = "gray55", lty = 3)
    }
    for (i in seq_along(styles$methods)) {
      method <- styles$methods[i]
      method_curve <- panel[panel$method == method, , drop = FALSE]
      method_curve <- method_curve[order(method_curve$alpha), , drop = FALSE]
      if (nrow(method_curve) == 0) {
        next
      }
      polygon(
        x = c(method_curve$alpha, rev(method_curve$alpha)),
        y = c(
          pmax(0, method_curve$power_ci_lower),
          rev(pmin(1, method_curve$power_ci_upper))
        ),
        col = grDevices::adjustcolor(styles$col[i], alpha.f = 0.14),
        border = NA
      )
      lines(
        method_curve$alpha,
        method_curve$mean_power,
        col = styles$col[i],
        lty = styles$lty[i],
        lwd = 2.4
      )
    }
    legend(
      "bottomright",
      legend = styles$methods,
      col = styles$col,
      lty = styles$lty,
      lwd = 2.4,
      bty = "n",
      cex = 0.78
    )
  }
  mtext(title, outer = TRUE, line = if (is.null(subtitle)) 0.8 else 2.1, cex = 1.2)
  if (!is.null(subtitle)) {
    mtext(subtitle, outer = TRUE, line = 0.5, cex = 0.85)
  }
  invisible(mc_curve)
}

plot_empirical_fdr_alpha_curve_grid <- function(alpha_curve,
                                                panel_var = "panel_label",
                                                file = NULL,
                                                title = "Empirical FDR as a function of alpha",
                                                alpha_reference = 0.05,
                                                ncol = 2,
                                                width = 2200,
                                                height = 1600,
                                                res = 180) {
  if (!all(c(panel_var, "method", "alpha", "empirical_fdr") %in% colnames(alpha_curve))) {
    stop("alpha_curve must contain panel_var, method, alpha, and empirical_fdr columns.")
  }

  panels <- unique(alpha_curve[[panel_var]])
  nrow <- ceiling(length(panels) / ncol)

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    png(file, width = width, height = height, res = res)
    on.exit(dev.off(), add = TRUE)
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(nrow, ncol), mar = c(4.1, 4.2, 3.0, 1.0), oma = c(0, 0, 3.0, 0))

  styles <- revision_method_styles(unique(alpha_curve$method))
  x_limits <- range(alpha_curve$alpha, na.rm = TRUE)
  y_upper <- max(c(x_limits[2], alpha_curve$empirical_fdr), na.rm = TRUE)
  y_limits <- c(0, min(1, y_upper * 1.08))

  for (panel in panels) {
    panel_curve <- alpha_curve[alpha_curve[[panel_var]] == panel, ]
    plot(
      NA,
      xlim = x_limits,
      ylim = y_limits,
      xlab = "Nominal FDR level alpha",
      ylab = "Empirical FDR",
      main = panel
    )
    grid(col = "gray90")
    abline(a = 0, b = 1, col = "gray35", lty = 3)
    if (!is.null(alpha_reference)) {
      abline(v = alpha_reference, col = "gray55", lty = 3)
    }
    for (i in seq_along(styles$methods)) {
      method <- styles$methods[i]
      method_curve <- panel_curve[panel_curve$method == method, ]
      method_curve <- method_curve[order(method_curve$alpha), ]
      lines(
        method_curve$alpha,
        method_curve$empirical_fdr,
        col = styles$col[i],
        lty = styles$lty[i],
        lwd = 2.4
      )
    }
    legend(
      "topleft",
      legend = c(styles$methods, "Nominal alpha"),
      col = c(styles$col, "gray35"),
      lty = c(styles$lty, 3),
      lwd = c(rep(2.4, length(styles$methods)), 1.8),
      bty = "n",
      cex = 0.72
    )
  }

  mtext(title, outer = TRUE, cex = 1.2, font = 2)
  invisible(alpha_curve)
}

run_revision_simulation <- function(J = 300,
                                    class_probs = c(
                                      constant = 0.80,
                                      linear = 0.10,
                                      nonlinear_iwp2 = 0.10
                                    ),
                                    scenario = "appendixB_like",
                                    alpha = 0.05,
                                    seed = 12345,
                                    num_cores = 1,
                                    num_basis = 20,
                                    grid = default_revision_grid(),
                                    penalty = 10,
                                    sd_values = c(0.1, 0.3, 0.5),
                                    correlation = NULL,
                                    exact_class_counts = TRUE,
                                    output_dir = "output/revision_simulations",
                                    save_outputs = TRUE,
                                    verbose = FALSE) {
  dirs <- revision_output_dirs(output_dir)
  datasets <- simulate_revision_datasets(
    J = J,
    class_probs = class_probs,
    sd_values = sd_values,
    scenario = scenario,
    correlation = correlation,
    exact_class_counts = exact_class_counts,
    seed = seed
  )
  unit_info <- attr(datasets, "unit_info")

  fits <- fit_fash_for_revision(
    datasets = datasets,
    orders = c(1, 2),
    grid = grid,
    num_basis = num_basis,
    penalty = penalty,
    num_cores = num_cores,
    apply_bf = TRUE,
    verbose = verbose
  )

  fash_results <- evaluate_fash_fits(fits, unit_info, alpha = alpha)
  parametric_results <- rbind(
    evaluate_parametric_method(
      datasets = datasets,
      unit_info = unit_info,
      method = "linear_summary",
      target = "dynamic",
      alpha = alpha
    ),
    evaluate_parametric_method(
      datasets = datasets,
      unit_info = unit_info,
      method = "quadratic_summary",
      target = "nonlinear",
      alpha = alpha
    )
  )

  result_table <- rbind(fash_results, parametric_results)
  summary_table <- summarize_method_results(result_table)

  out <- list(
    datasets = datasets,
    unit_info = unit_info,
    fits = fits,
    result_table = result_table,
    summary_table = summary_table,
    settings = list(
      J = J,
      class_probs = class_probs,
      scenario = scenario,
      alpha = alpha,
      seed = seed,
      num_cores = num_cores,
      num_basis = num_basis,
      penalty = penalty,
      sd_values = sd_values,
      exact_class_counts = exact_class_counts,
      correlation = correlation
    )
  )

  if (save_outputs) {
    stem <- paste0(
      scenario,
      "_sd", numeric_vector_label(sd_values),
      "_seed", seed,
      "_J", J
    )
    saveRDS(out, file = file.path(dirs$raw, paste0(stem, ".rds")))
    write.csv(
      result_table,
      file = file.path(dirs$summary, paste0(stem, "_unit_results.csv")),
      row.names = FALSE
    )
    write.csv(
      summary_table,
      file = file.path(dirs$summary, paste0(stem, "_method_summary.csv")),
      row.names = FALSE
    )
  }

  out
}

run_constant_iwp2_simplified_comparison <- function(J = 300,
                                                    pi_dynamic = 0.20,
                                                    dynamic_class = "nonlinear_iwp2",
                                                    alpha = 0.05,
                                                    seed = 12345,
                                                    sigma_beta = 1,
                                                    estimate_sigma = FALSE,
                                                    sigma_beta_grid = exp(seq(log(0.05), log(5), length.out = 25)),
                                                    num_cores = 1,
                                                    num_basis = 20,
                                                    grid = default_revision_grid(),
                                                    penalty = 10,
                                                    sd_values = c(0.1, 0.3, 0.5),
                                                    dynamic_amplitude = 2,
                                                    bspline_df = 6,
                                                    bspline_coefficient_sd = 1,
                                                    smooth_abrupt_width = 0.8,
                                                    smooth_abrupt_tau_values = c(5, 8, 11),
                                                    multi_smooth_abrupt_n_jumps_values = 1:3,
                                                    multi_smooth_abrupt_tau_values = 3:12,
                                                    transient_bspline_df = 10,
                                                    transient_bspline_degree = 3,
                                                    exact_class_counts = TRUE,
                                                    output_dir = "output/revision_simulations",
                                                    save_outputs = TRUE,
                                                    verbose = FALSE) {
  dirs <- revision_output_dirs(output_dir)
  scenario <- paste0("constant_vs_", dynamic_class, "_simplified")
  class_probs <- c(constant = 1 - pi_dynamic, pi_dynamic)
  names(class_probs)[2] <- dynamic_class
  datasets <- simulate_revision_datasets(
    J = J,
    class_probs = class_probs,
    sd_values = sd_values,
    scenario = scenario,
    amplitude = dynamic_amplitude,
    bspline_df = bspline_df,
    bspline_coefficient_sd = bspline_coefficient_sd,
    smooth_abrupt_width = smooth_abrupt_width,
    smooth_abrupt_tau_values = smooth_abrupt_tau_values,
    multi_smooth_abrupt_n_jumps_values = multi_smooth_abrupt_n_jumps_values,
    multi_smooth_abrupt_tau_values = multi_smooth_abrupt_tau_values,
    transient_bspline_df = transient_bspline_df,
    transient_bspline_degree = transient_bspline_degree,
    exact_class_counts = exact_class_counts,
    seed = seed
  )
  unit_info <- attr(datasets, "unit_info")

  fash_fits <- fit_fash_for_revision(
    datasets = datasets,
    orders = 1,
    grid = grid,
    num_basis = num_basis,
    penalty = penalty,
    num_cores = num_cores,
    apply_bf = TRUE,
    verbose = verbose
  )
  simplified_fit <- fit_simplified_fash(
    datasets = datasets,
    sigma_beta = sigma_beta,
    estimate_sigma = estimate_sigma,
    sigma_beta_grid = sigma_beta_grid,
    scale_time = TRUE
  )
  simplified_fit_bf <- BF_update_simplified_fash(simplified_fit)

  result_table <- rbind(
    evaluate_lfdr_method(
      lfdr = get_fash_lfdr(fash_fits$fash_iwp1_raw),
      unit_info = unit_info,
      method = "FASH-IWP1-Raw",
      target = "dynamic",
      alpha = alpha
    ),
    evaluate_lfdr_method(
      lfdr = get_fash_lfdr(fash_fits$fash_iwp1_bf),
      unit_info = unit_info,
      method = "FASH-IWP1-BF",
      target = "dynamic",
      alpha = alpha
    ),
    evaluate_simplified_fash_fit(
      fit = simplified_fit,
      unit_info = unit_info,
      alpha = alpha,
      method = "FASH-linear-Raw"
    ),
    evaluate_simplified_fash_fit(
      fit = simplified_fit_bf,
      unit_info = unit_info,
      alpha = alpha,
      method = "FASH-linear-BF"
    )
  )
  summary_table <- summarize_method_results(result_table)

  out <- list(
    datasets = datasets,
    unit_info = unit_info,
    fash_fits = fash_fits,
    simplified_fit = simplified_fit,
    simplified_fit_bf = simplified_fit_bf,
    result_table = result_table,
    summary_table = summary_table,
    settings = list(
      J = J,
      pi_dynamic = pi_dynamic,
      alpha = alpha,
      seed = seed,
      sigma_beta = sigma_beta,
      estimate_sigma = estimate_sigma,
      sigma_beta_grid = sigma_beta_grid,
      num_cores = num_cores,
      num_basis = num_basis,
      penalty = penalty,
      sd_values = sd_values,
      dynamic_class = dynamic_class,
      dynamic_amplitude = dynamic_amplitude,
      bspline_df = bspline_df,
      bspline_coefficient_sd = bspline_coefficient_sd,
      smooth_abrupt_width = smooth_abrupt_width,
      smooth_abrupt_tau_values = smooth_abrupt_tau_values,
      multi_smooth_abrupt_n_jumps_values = multi_smooth_abrupt_n_jumps_values,
      multi_smooth_abrupt_tau_values = multi_smooth_abrupt_tau_values,
      transient_bspline_df = transient_bspline_df,
      transient_bspline_degree = transient_bspline_degree,
      exact_class_counts = exact_class_counts
    )
  )

  if (save_outputs) {
    stem <- paste0(
      scenario,
      "_sd", numeric_vector_label(sd_values),
      "_sigmabeta",
      if (estimate_sigma) "estimated" else numeric_vector_label(sigma_beta),
      "_amp", numeric_vector_label(dynamic_amplitude),
      "_bsdf", bspline_df,
      "_tdf", transient_bspline_df,
      "_width", numeric_vector_label(smooth_abrupt_width),
      "_seed", seed,
      "_J", J
    )
    saveRDS(out, file = file.path(dirs$raw, paste0(stem, ".rds")))
    write.csv(
      result_table,
      file = file.path(dirs$summary, paste0(stem, "_unit_results.csv")),
      row.names = FALSE
    )
    write.csv(
      summary_table,
      file = file.path(dirs$summary, paste0(stem, "_method_summary.csv")),
      row.names = FALSE
    )
  }

  out
}

validate_genotype_matrix <- function(G) {
  G <- as.matrix(G)
  storage.mode(G) <- "numeric"
  if (any(!is.finite(G))) {
    stop("G must contain only finite genotype values.")
  }
  if (is.null(rownames(G))) {
    rownames(G) <- sprintf("donor_%02d", seq_len(nrow(G)))
  }
  if (is.null(colnames(G))) {
    colnames(G) <- sprintf("variant_%04d", seq_len(ncol(G)))
  }
  genotype_sd <- apply(G, 2, sd)
  if (any(genotype_sd <= 0)) {
    bad <- colnames(G)[genotype_sd <= 0]
    stop(
      "G contains monomorphic variants; remove them before fitting. First affected variant: ",
      bad[1]
    )
  }
  G
}

validate_covariate_matrix <- function(covariates, n_donors) {
  if (is.null(covariates)) {
    return(NULL)
  }
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "numeric"
  if (nrow(covariates) != n_donors) {
    stop("covariates must have the same number of rows as G.")
  }
  if (any(!is.finite(covariates))) {
    stop("covariates must contain only finite values.")
  }
  if (is.null(colnames(covariates))) {
    colnames(covariates) <- sprintf("PC%d", seq_len(ncol(covariates)))
  }
  if (is.null(rownames(covariates))) {
    rownames(covariates) <- sprintf("donor_%02d", seq_len(nrow(covariates)))
  }
  covariates
}

simulate_genotype_matrix <- function(n_donors = 19,
                                     n_variants = 1000,
                                     maf_range = c(0.1, 0.5),
                                     seed = NULL,
                                     max_attempts = 100) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (length(maf_range) != 2 || any(maf_range <= 0) || any(maf_range >= 1)) {
    stop("maf_range must contain two values strictly between 0 and 1.")
  }
  maf_range <- sort(maf_range)
  maf <- runif(n_variants, min = maf_range[1], max = maf_range[2])
  G <- matrix(NA_real_, nrow = n_donors, ncol = n_variants)

  for (j in seq_len(n_variants)) {
    for (attempt in seq_len(max_attempts)) {
      g <- rbinom(n_donors, size = 2, prob = maf[j])
      if (sd(g) > 0) {
        G[, j] <- g
        break
      }
    }
    if (anyNA(G[, j])) {
      stop("Failed to simulate a polymorphic genotype for variant ", j, ".")
    }
  }

  rownames(G) <- sprintf("donor_%02d", seq_len(n_donors))
  colnames(G) <- sprintf("variant_%04d", seq_len(n_variants))
  variant_info <- data.frame(
    variant_index = seq_len(n_variants),
    variant_id = colnames(G),
    maf = maf,
    observed_maf = colMeans(G) / 2,
    genotype_sd = apply(G, 2, sd),
    stringsAsFactors = FALSE
  )

  list(G = G, variant_info = variant_info)
}

simulate_covariate_matrix <- function(n_donors = 19,
                                      n_covariates = 5,
                                      seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (n_covariates <= 0) {
    return(NULL)
  }
  covariates <- matrix(rnorm(n_donors * n_covariates), nrow = n_donors)
  covariates <- scale(covariates, center = TRUE, scale = TRUE)
  covariates <- as.matrix(covariates)
  rownames(covariates) <- sprintf("donor_%02d", seq_len(n_donors))
  colnames(covariates) <- sprintf("PC%d", seq_len(n_covariates))
  covariates
}

simulate_variant_effect_curves <- function(n_variants,
                                           time_grid = make_time_grid(),
                                           class_probs = c(
                                             dynamic_bspline = 0.20,
                                             constant = 0.40,
                                             zero = 0.40
                                           ),
                                           scenario = "genotype_bspline_dynamic_eqtl",
                                           dynamic_amplitude = 2,
                                           bspline_df = 6,
                                           bspline_coefficient_sd = 1,
                                           transient_bspline_df = 16,
                                           transient_bspline_degree = 3,
                                           constant_sd = 1,
                                           dynamic_main_effect_sd = 0,
                                           exact_class_counts = TRUE,
                                           seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (abs(sum(class_probs) - 1) > 1e-8) {
    stop("class_probs must sum to 1.")
  }
  allowed_classes <- c(
    "dynamic_bspline",
    "dynamic_local_bspline_transient",
    "constant",
    "zero"
  )
  if (!all(names(class_probs) %in% allowed_classes)) {
    stop("Unsupported class in class_probs.")
  }
  if (!is.finite(dynamic_main_effect_sd) || dynamic_main_effect_sd < 0) {
    stop("dynamic_main_effect_sd must be finite and nonnegative.")
  }

  effect_class <- sample_effect_classes(
    J = n_variants,
    class_probs = class_probs,
    exact_class_counts = exact_class_counts
  )
  beta_matrix <- matrix(0, nrow = n_variants, ncol = length(time_grid))

  for (j in seq_len(n_variants)) {
    beta_matrix[j, ] <- switch(
      effect_class[j],
      dynamic_bspline = simulate_bspline_effect(
        x = time_grid,
        amplitude = dynamic_amplitude,
        df = bspline_df,
        coefficient_sd = bspline_coefficient_sd
      ),
      dynamic_local_bspline_transient = simulate_local_bspline_transient_effect(
        x = time_grid,
        amplitude = dynamic_amplitude,
        df = transient_bspline_df,
        degree = transient_bspline_degree
      ),
      constant = rep(rnorm(1, mean = 0, sd = constant_sd), length(time_grid)),
      zero = rep(0, length(time_grid)),
      stop("Unsupported effect class: ", effect_class[j])
    )
  }
  dynamic <- effect_class %in% c(
    "dynamic_bspline",
    "dynamic_local_bspline_transient"
  )
  genetic_main_effect <- rep(NA_real_, n_variants)
  if (any(dynamic)) {
    genetic_main_effect[dynamic] <- if (dynamic_main_effect_sd > 0) {
      stats::rnorm(sum(dynamic), mean = 0, sd = dynamic_main_effect_sd)
    } else {
      0
    }
    beta_matrix[dynamic, ] <- sweep(
      beta_matrix[dynamic, , drop = FALSE],
      1,
      genetic_main_effect[dynamic],
      `+`
    )
  }

  variant_id <- sprintf("variant_%04d", seq_len(n_variants))
  unit_id <- sprintf("%s_%04d", effect_class, seq_len(n_variants))
  rownames(beta_matrix) <- variant_id
  colnames(beta_matrix) <- sprintf("time_%02d", seq_along(time_grid))

  unit_info <- data.frame(
    unit_index = seq_len(n_variants),
    unit_id = unit_id,
    variant_id = variant_id,
    effect_class = effect_class,
    genetic_main_effect = genetic_main_effect,
    scenario = scenario,
    stringsAsFactors = FALSE
  )

  list(
    beta_matrix = beta_matrix,
    unit_info = unit_info,
    settings = list(
      dynamic_main_effect_sd = dynamic_main_effect_sd
    )
  )
}

simulate_eqtl_expression_from_genotypes <- function(G,
                                                    beta_matrix,
                                                    time_grid = make_time_grid(),
                                                    covariates = NULL,
                                                    expression_noise_sd = 1,
                                                    covariate_effect_sd = 0.5,
                                                    intercept_sd = 0,
                                                    seed = NULL,
                                                    expression_error_correlation = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  G <- validate_genotype_matrix(G)
  covariates <- validate_covariate_matrix(covariates, nrow(G))
  beta_matrix <- as.matrix(beta_matrix)
  storage.mode(beta_matrix) <- "numeric"

  if (nrow(beta_matrix) != ncol(G)) {
    stop("beta_matrix must have one row per variant in G.")
  }
  if (ncol(beta_matrix) != length(time_grid)) {
    stop("beta_matrix must have one column per time point.")
  }
  if (expression_noise_sd <= 0) {
    stop("expression_noise_sd must be positive.")
  }

  n_donors <- nrow(G)
  n_variants <- ncol(G)
  n_time <- length(time_grid)
  n_covariates <- if (is.null(covariates)) 0 else ncol(covariates)
  time_correlation <- if (is.null(expression_error_correlation)) {
    NULL
  } else {
    validate_time_correlation(expression_error_correlation, n_time = n_time)
  }
  use_identity_correlation <- !is.null(time_correlation) &&
    max(abs(time_correlation - diag(n_time))) <= 1e-12

  expression <- array(
    NA_real_,
    dim = c(n_donors, n_variants, n_time),
    dimnames = list(rownames(G), colnames(G), colnames(beta_matrix))
  )
  intercepts <- matrix(
    rnorm(n_variants * n_time, mean = 0, sd = intercept_sd),
    nrow = n_variants,
    ncol = n_time
  )
  covariate_effects <- if (n_covariates > 0) {
    array(
      NA_real_,
      dim = c(n_covariates, n_variants, n_time),
      dimnames = list(colnames(covariates), colnames(G), colnames(beta_matrix))
    )
  } else {
    NULL
  }
  expression_errors <- array(
    NA_real_,
    dim = c(n_donors, n_variants, n_time),
    dimnames = dimnames(expression)
  )

  for (tt in seq_len(n_time)) {
    genetic_mean <- sweep(G, 2, beta_matrix[, tt], `*`)
    covariate_mean <- matrix(0, nrow = n_donors, ncol = n_variants)
    if (n_covariates > 0) {
      gamma <- matrix(
        rnorm(n_covariates * n_variants, mean = 0, sd = covariate_effect_sd),
        nrow = n_covariates,
        ncol = n_variants
      )
      covariate_effects[, , tt] <- gamma
      covariate_mean <- covariates %*% gamma
    }
    expression_errors[, , tt] <- matrix(
      rnorm(n_donors * n_variants, mean = 0, sd = expression_noise_sd),
      nrow = n_donors,
      ncol = n_variants
    )
    expression[, , tt] <- genetic_mean + covariate_mean
  }

  if (!is.null(time_correlation) && !use_identity_correlation) {
    error_matrix <- matrix(
      expression_errors,
      nrow = n_donors * n_variants,
      ncol = n_time
    )
    error_matrix <- error_matrix %*% chol(time_correlation)
    expression_errors <- array(
      error_matrix,
      dim = c(n_donors, n_variants, n_time),
      dimnames = dimnames(expression)
    )
  }
  for (tt in seq_len(n_time)) {
    expression[, , tt] <- sweep(
      expression[, , tt] + expression_errors[, , tt],
      2,
      intercepts[, tt],
      `+`
    )
  }

  list(
    expression = expression,
    intercepts = intercepts,
    covariate_effects = covariate_effects,
    settings = list(
      expression_noise_sd = expression_noise_sd,
      covariate_effect_sd = covariate_effect_sd,
      intercept_sd = intercept_sd,
      expression_error_correlation = if (is.null(time_correlation)) {
        diag(n_time)
      } else {
        time_correlation
      },
      expression_error_correlation_supplied = !is.null(time_correlation)
    )
  )
}

correct_se_from_t <- function(beta_hat, se, df) {
  beta_hat <- as.matrix(beta_hat)
  se <- as.matrix(se)
  if (!all(dim(beta_hat) == dim(se))) {
    stop("beta_hat and se must have the same dimensions.")
  }

  df_matrix <- matrix(df, nrow = nrow(beta_hat), ncol = ncol(beta_hat))
  corrected <- se
  valid <- is.finite(beta_hat) &
    is.finite(se) &
    se > 0 &
    is.finite(df_matrix) &
    df_matrix > 0

  t_val <- beta_hat[valid] / se[valid]
  p_val <- 2 * stats::pt(abs(t_val), df = df_matrix[valid], lower.tail = FALSE)
  z_val <- stats::qnorm(1 - p_val / 2) * sign(t_val)
  corrected_values <- abs(beta_hat[valid]) / abs(z_val)
  corrected[valid] <- corrected_values

  unstable <- !is.finite(corrected) | corrected <= 0
  corrected[unstable] <- se[unstable]
  corrected
}

estimate_eqtl_summaries_from_genotypes <- function(G,
                                                   expression,
                                                   covariates = NULL,
                                                   apply_t_se_correction = TRUE) {
  G <- validate_genotype_matrix(G)
  covariates <- validate_covariate_matrix(covariates, nrow(G))

  if (is.list(expression) && !is.array(expression)) {
    n_time <- length(expression)
    expression_array <- array(NA_real_, dim = c(nrow(G), ncol(G), n_time))
    for (tt in seq_len(n_time)) {
      expression_array[, , tt] <- as.matrix(expression[[tt]])
    }
    expression <- expression_array
  }
  if (length(dim(expression)) != 3) {
    stop("expression must be a donor-by-variant-by-time array or a list of matrices.")
  }
  if (dim(expression)[1] != nrow(G) || dim(expression)[2] != ncol(G)) {
    stop("expression dimensions must match G.")
  }

  n_donors <- nrow(G)
  n_variants <- ncol(G)
  n_time <- dim(expression)[3]
  beta_hat <- matrix(NA_real_, nrow = n_variants, ncol = n_time)
  se_uncorrected <- matrix(NA_real_, nrow = n_variants, ncol = n_time)
  df <- rep(NA_integer_, n_variants)
  rank <- rep(NA_integer_, n_variants)

  rownames(beta_hat) <- colnames(G)
  rownames(se_uncorrected) <- colnames(G)
  colnames(beta_hat) <- paste0("time_", seq_len(n_time))
  colnames(se_uncorrected) <- paste0("time_", seq_len(n_time))

  for (j in seq_len(n_variants)) {
    design <- cbind(intercept = 1, G = G[, j])
    if (!is.null(covariates)) {
      design <- cbind(design, covariates)
    }
    design_rank <- qr(design)$rank
    if (design_rank < ncol(design)) {
      stop("The regression design is rank deficient for variant ", colnames(G)[j], ".")
    }
    df[j] <- n_donors - ncol(design)
    if (df[j] <= 0) {
      stop("The regression degrees of freedom must be positive.")
    }
    rank[j] <- design_rank

    xtx_inv <- solve(crossprod(design))
    hat_solver <- xtx_inv %*% t(design)
    y_j <- matrix(expression[, j, ], nrow = n_donors, ncol = n_time)
    coef_mat <- hat_solver %*% y_j
    fitted <- design %*% coef_mat
    residual <- y_j - fitted
    sigma2 <- colSums(residual^2) / df[j]

    beta_hat[j, ] <- coef_mat["G", ]
    se_uncorrected[j, ] <- sqrt(pmax(sigma2, 0) * xtx_inv["G", "G"])
  }

  df_matrix <- matrix(df, nrow = n_variants, ncol = n_time)
  se <- if (apply_t_se_correction) {
    correct_se_from_t(beta_hat, se_uncorrected, df_matrix)
  } else {
    se_uncorrected
  }
  colnames(se) <- colnames(se_uncorrected)
  rownames(se) <- rownames(se_uncorrected)

  list(
    beta_hat = beta_hat,
    se = se,
    se_uncorrected = se_uncorrected,
    df = df_matrix,
    rank = rank,
    apply_t_se_correction = apply_t_se_correction
  )
}

scaled_time_polynomial <- function(time_grid, degree = 1) {
  if (degree < 1) {
    stop("degree must be at least 1.")
  }
  t_scaled <- as.numeric(scale(time_grid, center = TRUE, scale = TRUE))
  out <- sapply(seq_len(degree), function(d) t_scaled^d)
  out <- as.matrix(out)
  colnames(out) <- paste0("time_degree_", seq_len(degree))
  out
}

make_direct_interaction_base_design <- function(n_donors,
                                                time_grid,
                                                covariates = NULL) {
  covariates <- validate_covariate_matrix(covariates, n_donors)
  n_time <- length(time_grid)
  time_index <- rep(seq_len(n_time), each = n_donors)
  time_design <- stats::model.matrix(~ 0 + factor(time_index))
  colnames(time_design) <- paste0("time_", seq_len(n_time))

  base_design <- time_design
  covariate_repeated <- NULL
  if (!is.null(covariates)) {
    covariate_repeated <- covariates[rep(seq_len(n_donors), times = n_time), , drop = FALSE]
    covariate_time_designs <- lapply(seq_len(ncol(covariates)), function(k) {
      out <- sweep(time_design, 1, covariate_repeated[, k], `*`)
      colnames(out) <- paste0(colnames(covariates)[k], "_", colnames(time_design))
      out
    })
    base_design <- do.call(cbind, c(list(base_design), covariate_time_designs))
  }

  list(
    base_design = base_design,
    time_index = time_index,
    time_basis = scaled_time_polynomial(time_grid, degree = 2),
    covariate_repeated = covariate_repeated
  )
}

linear_model_residual_summary <- function(y, design) {
  qr_design <- qr(design)
  fitted <- qr.fitted(qr_design, y)
  residual <- y - fitted
  list(
    rss = sum(residual^2),
    rank = qr_design$rank
  )
}

direct_interaction_lrt_pvalues <- function(G,
                                           expression,
                                           covariates = NULL,
                                           time_grid = make_time_grid(),
                                           interaction_degree = c(1, 2)) {
  interaction_degree <- match.arg(as.character(interaction_degree), c("1", "2"))
  interaction_degree <- as.integer(interaction_degree)
  G <- validate_genotype_matrix(G)
  covariates <- validate_covariate_matrix(covariates, nrow(G))
  if (length(dim(expression)) != 3) {
    stop("expression must be a donor-by-variant-by-time array.")
  }
  if (dim(expression)[1] != nrow(G) || dim(expression)[2] != ncol(G)) {
    stop("expression dimensions must match G.")
  }
  if (dim(expression)[3] != length(time_grid)) {
    stop("time_grid length must match the expression time dimension.")
  }

  n_donors <- nrow(G)
  n_variants <- ncol(G)
  n_time <- length(time_grid)
  n_observations <- n_donors * n_time
  base_info <- make_direct_interaction_base_design(
    n_donors = n_donors,
    time_grid = time_grid,
    covariates = covariates
  )
  base_design <- base_info$base_design
  time_basis <- base_info$time_basis[base_info$time_index, seq_len(interaction_degree), drop = FALSE]

  out <- data.frame(
    variant_id = colnames(G),
    lrt_statistic = NA_real_,
    lrt_df = NA_integer_,
    pvalue = NA_real_,
    rss_null = NA_real_,
    rss_full = NA_real_,
    rank_null = NA_integer_,
    rank_full = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (j in seq_len(n_variants)) {
    genotype_repeated <- rep(G[, j], times = n_time)
    null_design <- cbind(base_design, G = genotype_repeated)
    interaction_design <- sweep(time_basis, 1, genotype_repeated, `*`)
    colnames(interaction_design) <- paste0("G_x_", colnames(time_basis))
    full_design <- cbind(null_design, interaction_design)

    y <- as.vector(expression[, j, ])
    null_fit <- linear_model_residual_summary(y, null_design)
    full_fit <- linear_model_residual_summary(y, full_design)
    lrt_df <- full_fit$rank - null_fit$rank
    lrt_statistic <- if (lrt_df > 0 && full_fit$rss > 0 && null_fit$rss >= full_fit$rss) {
      n_observations * log(null_fit$rss / full_fit$rss)
    } else {
      0
    }

    out$lrt_statistic[j] <- lrt_statistic
    out$lrt_df[j] <- lrt_df
    out$pvalue[j] <- if (lrt_df > 0) {
      stats::pchisq(lrt_statistic, df = lrt_df, lower.tail = FALSE)
    } else {
      NA_real_
    }
    out$rss_null[j] <- null_fit$rss
    out$rss_full[j] <- full_fit$rss
    out$rank_null[j] <- null_fit$rank
    out$rank_full[j] <- full_fit$rank
  }

  out
}

evaluate_direct_interaction_lrt <- function(G,
                                            expression,
                                            covariates,
                                            unit_info,
                                            time_grid = make_time_grid(),
                                            interaction_degree = c(1, 2),
                                            method = NULL,
                                            alpha = 0.05,
                                            adjust_method = "BH") {
  interaction_degree <- match.arg(as.character(interaction_degree), c("1", "2"))
  interaction_degree <- as.integer(interaction_degree)
  if (is.null(method)) {
    method <- if (interaction_degree == 1) "Direct-linear-LRT" else "Direct-quadratic-LRT"
  }

  lrt <- direct_interaction_lrt_pvalues(
    G = G,
    expression = expression,
    covariates = covariates,
    time_grid = time_grid,
    interaction_degree = interaction_degree
  )
  if (nrow(lrt) != nrow(unit_info)) {
    stop("LRT output and unit_info must have the same number of rows.")
  }
  qvalues <- stats::p.adjust(lrt$pvalue, method = adjust_method)

  result <- data.frame(
    unit_index = unit_info$unit_index,
    unit_id = unit_info$unit_id,
    effect_class = unit_info$effect_class,
    scenario = unit_info$scenario,
    target = "dynamic",
    method = method,
    score_type = "pvalue",
    score = lrt$pvalue,
    adjusted_score = qvalues,
    alpha = alpha,
    selected = qvalues <= alpha,
    true_null = is_null_class(unit_info$effect_class, "dynamic"),
    stringsAsFactors = FALSE
  )
  list(result = result, lrt = lrt)
}

add_direct_interaction_results_to_genotype_output <- function(out,
                                                             alpha = out$settings$alpha,
                                                             adjust_method = "BH") {
  required_fields <- c("genotype", "expression", "covariates", "unit_info", "settings")
  missing_fields <- setdiff(required_fields, names(out))
  if (length(missing_fields) > 0) {
    stop("out is missing required fields: ", paste(missing_fields, collapse = ", "))
  }
  time_grid <- out$settings$time_grid
  direct_linear <- evaluate_direct_interaction_lrt(
    G = out$genotype,
    expression = out$expression,
    covariates = out$covariates,
    unit_info = out$unit_info,
    time_grid = time_grid,
    interaction_degree = 1,
    method = "Direct-linear-LRT",
    alpha = alpha,
    adjust_method = adjust_method
  )
  direct_quadratic <- evaluate_direct_interaction_lrt(
    G = out$genotype,
    expression = out$expression,
    covariates = out$covariates,
    unit_info = out$unit_info,
    time_grid = time_grid,
    interaction_degree = 2,
    method = "Direct-quadratic-LRT",
    alpha = alpha,
    adjust_method = adjust_method
  )

  keep_existing <- !out$result_table$method %in% c("Direct-linear-LRT", "Direct-quadratic-LRT")
  out$result_table <- rbind(
    out$result_table[keep_existing, ],
    direct_linear$result,
    direct_quadratic$result
  )
  out$direct_interaction_lrt <- list(
    linear = direct_linear$lrt,
    quadratic = direct_quadratic$lrt,
    adjust_method = adjust_method
  )
  out$summary_table <- summarize_method_results(out$result_table)
  out$alpha_curve <- compute_alpha_curve(
    result_table = out$result_table,
    alpha_grid = seq(0, 0.20, by = 0.005)
  )
  out
}

empirical_fdr_qvalues <- function(observed_pvalues,
                                  null_pvalues,
                                  lambda = 0.5,
                                  pi0_method = c("conservative", "storey", "fixed"),
                                  fixed_pi0 = NULL) {
  pi0_method <- match.arg(pi0_method)
  if (pi0_method == "fixed") {
    if (is.null(fixed_pi0) ||
        length(fixed_pi0) != 1 ||
        !is.finite(fixed_pi0) ||
        fixed_pi0 < 0 ||
        fixed_pi0 > 1) {
      stop("fixed_pi0 must be a finite scalar between 0 and 1 when pi0_method = 'fixed'.")
    }
  }
  qvalues <- rep(NA_real_, length(observed_pvalues))
  valid_observed <- is.finite(observed_pvalues) &
    observed_pvalues >= 0 &
    observed_pvalues <= 1
  valid_null <- is.finite(null_pvalues) &
    null_pvalues >= 0 &
    null_pvalues <= 1

  if (!any(valid_observed) || !any(valid_null)) {
    return(list(
      qvalue = qvalues,
      threshold_table = data.frame(),
      pi0_multiplier = NA_real_
    ))
  }

  observed <- observed_pvalues[valid_observed]
  null <- null_pvalues[valid_null]
  thresholds <- sort(unique(observed))
  observed_count <- vapply(thresholds, function(t) sum(observed <= t), integer(1))
  null_count <- vapply(thresholds, function(t) sum(null <= t), integer(1))
  observed_cdf <- observed_count / length(observed)
  null_cdf <- null_count / length(null)

  pi0_multiplier <- 1
  if (pi0_method == "storey") {
    observed_tail <- mean(observed > lambda)
    null_tail <- mean(null > lambda)
    pi0_multiplier <- if (null_tail > 0) observed_tail / null_tail else 1
  } else if (pi0_method == "fixed") {
    pi0_multiplier <- fixed_pi0
  }

  efdr <- pi0_multiplier * null_cdf / observed_cdf
  efdr <- pmin(pmax(efdr, 0), 1)
  threshold_qvalue <- rev(cummin(rev(efdr)))
  qvalues[valid_observed] <- threshold_qvalue[match(observed, thresholds)]

  threshold_table <- data.frame(
    threshold = thresholds,
    observed_count = observed_count,
    null_count = null_count,
    observed_cdf = observed_cdf,
    null_cdf = null_cdf,
    efdr = efdr,
    qvalue = threshold_qvalue,
    pi0_method = pi0_method,
    pi0_multiplier = pi0_multiplier,
    stringsAsFactors = FALSE
  )

  list(
    qvalue = qvalues,
    threshold_table = threshold_table,
    pi0_multiplier = pi0_multiplier
  )
}

direct_interaction_degree_label <- function(interaction_degree) {
  interaction_degree <- as.integer(interaction_degree)
  if (!interaction_degree %in% c(1L, 2L)) {
    stop("interaction_degree must be 1 or 2.")
  }
  if (interaction_degree == 1L) "linear" else "quadratic"
}

direct_interaction_method_name <- function(interaction_degree,
                                           adjustment = c("BH", "eFDR", "eFDR_true_pi0")) {
  adjustment <- match.arg(adjustment)
  degree_label <- direct_interaction_degree_label(interaction_degree)
  base <- if (degree_label == "linear") {
    "Direct-linear-LRT"
  } else {
    "Direct-quadratic-LRT"
  }
  if (adjustment == "eFDR") {
    paste0(base, "-eFDR")
  } else if (adjustment == "eFDR_true_pi0") {
    paste0(base, "-eFDR-true-pi0")
  } else {
    base
  }
}

dynamic_null_proportion <- function(unit_info, target = "dynamic") {
  mean(is_null_class(unit_info$effect_class, target))
}

get_observed_direct_interaction_lrt <- function(out, interaction_degree) {
  degree_label <- direct_interaction_degree_label(interaction_degree)
  if (!is.null(out$direct_interaction_lrt[[degree_label]])) {
    return(out$direct_interaction_lrt[[degree_label]])
  }
  direct_interaction_lrt_pvalues(
    G = out$genotype,
    expression = out$expression,
    covariates = out$covariates,
    time_grid = out$settings$time_grid,
    interaction_degree = interaction_degree
  )
}

compute_direct_interaction_permutation_null <- function(out,
                                                        n_permutations = 100,
                                                        interaction_degrees = c(1, 2),
                                                        seed = 22345,
                                                        permute_covariates_with_expression = TRUE,
                                                        num_cores = 1,
                                                        verbose = TRUE) {
  required_fields <- c("genotype", "expression", "covariates", "settings")
  missing_fields <- setdiff(required_fields, names(out))
  if (length(missing_fields) > 0) {
    stop("out is missing required fields: ", paste(missing_fields, collapse = ", "))
  }
  if (!is.numeric(n_permutations) || length(n_permutations) != 1 || n_permutations < 1) {
    stop("n_permutations must be a positive integer.")
  }
  if (!is.numeric(num_cores) || length(num_cores) != 1 || num_cores < 1) {
    stop("num_cores must be a positive integer.")
  }
  n_permutations <- as.integer(n_permutations)
  num_cores <- as.integer(num_cores)
  interaction_degrees <- sort(unique(as.integer(interaction_degrees)))
  if (!all(interaction_degrees %in% c(1L, 2L))) {
    stop("interaction_degrees must contain only 1 and/or 2.")
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }
  n_donors <- nrow(out$genotype)
  n_variants <- ncol(out$genotype)
  null_pvalues <- lapply(interaction_degrees, function(degree) {
    matrix(NA_real_, nrow = n_permutations, ncol = n_variants)
  })
  names(null_pvalues) <- vapply(
    interaction_degrees,
    direct_interaction_degree_label,
    character(1)
  )
  permutation_index <- t(replicate(n_permutations, sample.int(n_donors)))
  colnames(permutation_index) <- rownames(out$genotype)

  evaluate_permutation <- function(b) {
    perm <- permutation_index[b, ]
    expression_perm <- out$expression[perm, , , drop = FALSE]
    covariates_perm <- out$covariates
    if (permute_covariates_with_expression && !is.null(out$covariates)) {
      covariates_perm <- out$covariates[perm, , drop = FALSE]
    }

    pvalues <- lapply(interaction_degrees, function(degree) {
      degree_label <- direct_interaction_degree_label(degree)
      lrt <- direct_interaction_lrt_pvalues(
        G = out$genotype,
        expression = expression_perm,
        covariates = covariates_perm,
        time_grid = out$settings$time_grid,
        interaction_degree = degree
      )
      lrt$pvalue
    })
    names(pvalues) <- vapply(interaction_degrees, direct_interaction_degree_label, character(1))
    pvalues
  }

  permutation_results <- if (num_cores > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(
      seq_len(n_permutations),
      evaluate_permutation,
      mc.cores = min(num_cores, n_permutations),
      mc.preschedule = FALSE
    )
  } else {
    lapply(seq_len(n_permutations), evaluate_permutation)
  }
  for (b in seq_len(n_permutations)) {
    for (degree_label in names(null_pvalues)) {
      null_pvalues[[degree_label]][b, ] <- permutation_results[[b]][[degree_label]]
    }
  }
  if (verbose) {
    message(
      "Completed direct-interaction eFDR permutations: ",
      n_permutations,
      " using ",
      min(num_cores, n_permutations),
      " core(s)."
    )
  }

  list(
    null_pvalues = null_pvalues,
    permutation_index = permutation_index,
    settings = list(
      n_permutations = n_permutations,
      interaction_degrees = interaction_degrees,
      seed = seed,
      num_cores = num_cores,
      permute_covariates_with_expression = permute_covariates_with_expression,
      permutation_unit = "donor",
      permutation_scheme = if (permute_covariates_with_expression) {
        "donor-level expression and covariate row permutation relative to genotype"
      } else {
        "donor-level expression permutation relative to genotype"
      }
    )
  )
}

evaluate_direct_interaction_efdr <- function(out,
                                             permutation_null,
                                             interaction_degree = c(1, 2),
                                             alpha = out$settings$alpha,
                                             lambda = 0.5,
                                             pi0_method = c("conservative", "storey", "fixed"),
                                             fixed_pi0 = NULL,
                                             method_adjustment = NULL) {
  interaction_degree <- match.arg(as.character(interaction_degree), c("1", "2"))
  interaction_degree <- as.integer(interaction_degree)
  pi0_method <- match.arg(pi0_method)
  if (is.null(method_adjustment)) {
    method_adjustment <- if (pi0_method == "fixed") "eFDR_true_pi0" else "eFDR"
  }
  degree_label <- direct_interaction_degree_label(interaction_degree)
  if (is.null(permutation_null$null_pvalues[[degree_label]])) {
    stop("permutation_null does not contain null p-values for degree ", interaction_degree, ".")
  }

  lrt <- get_observed_direct_interaction_lrt(out, interaction_degree)
  efdr <- empirical_fdr_qvalues(
    observed_pvalues = lrt$pvalue,
    null_pvalues = as.vector(permutation_null$null_pvalues[[degree_label]]),
    lambda = lambda,
    pi0_method = pi0_method,
    fixed_pi0 = fixed_pi0
  )

  result <- data.frame(
    unit_index = out$unit_info$unit_index,
    unit_id = out$unit_info$unit_id,
    effect_class = out$unit_info$effect_class,
    scenario = out$unit_info$scenario,
    target = "dynamic",
    method = direct_interaction_method_name(interaction_degree, adjustment = method_adjustment),
    score_type = "pvalue",
    score = lrt$pvalue,
    adjusted_score = efdr$qvalue,
    alpha = alpha,
    selected = efdr$qvalue <= alpha,
    true_null = is_null_class(out$unit_info$effect_class, "dynamic"),
    stringsAsFactors = FALSE
  )

  list(result = result, lrt = lrt, efdr = efdr)
}

add_direct_interaction_efdr_results_to_genotype_output <- function(out,
                                                                   n_permutations = 100,
                                                                   alpha = out$settings$alpha,
                                                                   seed = 22345,
                                                                   lambda = 0.5,
                                                                   pi0_method = c("conservative", "storey", "fixed"),
                                                                   true_pi0 = NULL,
                                                                   include_true_pi0 = TRUE,
                                                                   permute_covariates_with_expression = TRUE,
                                                                   num_cores = 1,
                                                                   overwrite = FALSE,
                                                                   verbose = TRUE) {
  pi0_method <- match.arg(pi0_method)
  if (is.null(true_pi0)) {
    true_pi0 <- dynamic_null_proportion(out$unit_info, target = "dynamic")
  }
  existing <- out$direct_interaction_efdr
  can_reuse <- !overwrite &&
    !is.null(existing) &&
    !is.null(existing$null_pvalues) &&
    !is.null(existing$settings$n_permutations) &&
    existing$settings$n_permutations >= n_permutations &&
    identical(
      existing$settings$permute_covariates_with_expression,
      permute_covariates_with_expression
    )

  permutation_null <- if (can_reuse) {
    existing
  } else {
    compute_direct_interaction_permutation_null(
      out = out,
      n_permutations = n_permutations,
      interaction_degrees = c(1, 2),
      seed = seed,
      permute_covariates_with_expression = permute_covariates_with_expression,
      num_cores = num_cores,
      verbose = verbose
    )
  }

  direct_linear_efdr <- evaluate_direct_interaction_efdr(
    out = out,
    permutation_null = permutation_null,
    interaction_degree = 1,
    alpha = alpha,
    lambda = lambda,
    pi0_method = pi0_method
  )
  direct_quadratic_efdr <- evaluate_direct_interaction_efdr(
    out = out,
    permutation_null = permutation_null,
    interaction_degree = 2,
    alpha = alpha,
    lambda = lambda,
    pi0_method = pi0_method
  )
  direct_linear_efdr_true_pi0 <- NULL
  direct_quadratic_efdr_true_pi0 <- NULL
  if (include_true_pi0) {
    direct_linear_efdr_true_pi0 <- evaluate_direct_interaction_efdr(
      out = out,
      permutation_null = permutation_null,
      interaction_degree = 1,
      alpha = alpha,
      lambda = lambda,
      pi0_method = "fixed",
      fixed_pi0 = true_pi0,
      method_adjustment = "eFDR_true_pi0"
    )
    direct_quadratic_efdr_true_pi0 <- evaluate_direct_interaction_efdr(
      out = out,
      permutation_null = permutation_null,
      interaction_degree = 2,
      alpha = alpha,
      lambda = lambda,
      pi0_method = "fixed",
      fixed_pi0 = true_pi0,
      method_adjustment = "eFDR_true_pi0"
    )
  }

  efdr_methods <- c(
    "Direct-linear-LRT-eFDR",
    "Direct-quadratic-LRT-eFDR",
    "Direct-linear-LRT-eFDR-true-pi0",
    "Direct-quadratic-LRT-eFDR-true-pi0"
  )
  keep_existing <- !out$result_table$method %in% efdr_methods
  result_pieces <- list(
    out$result_table[keep_existing, ],
    direct_linear_efdr$result,
    direct_quadratic_efdr$result
  )
  if (include_true_pi0) {
    result_pieces <- c(
      result_pieces,
      list(
        direct_linear_efdr_true_pi0$result,
        direct_quadratic_efdr_true_pi0$result
      )
    )
  }
  out$result_table <- do.call(rbind, result_pieces)
  base_settings <- permutation_null$settings
  base_settings <- base_settings[
    !names(base_settings) %in% c("lambda", "pi0_method", "true_pi0", "include_true_pi0")
  ]
  threshold_tables <- list(
    linear = direct_linear_efdr$efdr$threshold_table,
    quadratic = direct_quadratic_efdr$efdr$threshold_table
  )
  pi0_multipliers <- c(
    linear = direct_linear_efdr$efdr$pi0_multiplier,
    quadratic = direct_quadratic_efdr$efdr$pi0_multiplier
  )
  if (include_true_pi0) {
    threshold_tables$linear_true_pi0 <- direct_linear_efdr_true_pi0$efdr$threshold_table
    threshold_tables$quadratic_true_pi0 <- direct_quadratic_efdr_true_pi0$efdr$threshold_table
    pi0_multipliers <- c(
      pi0_multipliers,
      linear_true_pi0 = direct_linear_efdr_true_pi0$efdr$pi0_multiplier,
      quadratic_true_pi0 = direct_quadratic_efdr_true_pi0$efdr$pi0_multiplier
    )
  }
  out$direct_interaction_efdr <- list(
    null_pvalues = permutation_null$null_pvalues,
    permutation_index = permutation_null$permutation_index,
    settings = c(
      base_settings,
      list(
        lambda = lambda,
        pi0_method = pi0_method,
        true_pi0 = true_pi0,
        include_true_pi0 = include_true_pi0
      )
    ),
    threshold_tables = threshold_tables,
    pi0_multipliers = pi0_multipliers
  )
  out$summary_table <- summarize_method_results(out$result_table)
  out$alpha_curve <- compute_alpha_curve(
    result_table = out$result_table,
    alpha_grid = seq(0, 0.20, by = 0.005)
  )
  out
}

make_fash_datasets_from_eqtl_summary <- function(beta_hat,
                                                 se,
                                                 true_beta,
                                                 time_grid,
                                                 unit_info,
                                                 scenario = "genotype_bspline_dynamic_eqtl") {
  beta_hat <- as.matrix(beta_hat)
  se <- as.matrix(se)
  true_beta <- as.matrix(true_beta)
  if (!all(dim(beta_hat) == dim(se)) || !all(dim(beta_hat) == dim(true_beta))) {
    stop("beta_hat, se, and true_beta must have the same dimensions.")
  }
  if (ncol(beta_hat) != length(time_grid)) {
    stop("time_grid length must match the number of summary columns.")
  }
  if (nrow(unit_info) != nrow(beta_hat)) {
    stop("unit_info must have one row per variant.")
  }

  unit_info$scenario <- scenario
  datasets <- vector("list", nrow(beta_hat))
  for (j in seq_len(nrow(beta_hat))) {
    datasets[[j]] <- data.frame(
      x = time_grid,
      y = beta_hat[j, ],
      sd = se[j, ],
      truef = true_beta[j, ],
      effect_class = unit_info$effect_class[j],
      unit_id = unit_info$unit_id[j],
      variant_id = unit_info$variant_id[j],
      scenario = scenario,
      stringsAsFactors = FALSE
    )
  }
  names(datasets) <- unit_info$unit_id
  attr(datasets, "unit_info") <- unit_info
  datasets
}

genotype_dynamic_eqtl_output_stem <- function(n_donors,
                                              n_variants,
                                              time_grid,
                                              n_covariates,
                                              expression_noise_sd,
                                              class_probs = c(
                                                dynamic_bspline = 0.20,
                                                constant = 0.40,
                                                zero = 0.40
                                              ),
                                              seed,
                                              scenario = "genotype_bspline_dynamic_eqtl") {
  paste0(
    scenario,
    "_N", n_donors,
    "_J", n_variants,
    "_T", length(time_grid),
    "_pc", n_covariates,
    "_mix", numeric_vector_label(unname(class_probs)),
    "_noise", numeric_vector_label(expression_noise_sd),
    "_seed", seed
  )
}

genotype_bspline_eqtl_output_stem <- genotype_dynamic_eqtl_output_stem

genotype_cosine_multipeak_output_stem <- function(
    n_donors,
    n_variants,
    time_grid,
    n_covariates,
    expression_noise_sd,
    width_half,
    target_centered_rms,
    shape_cell_probs = NULL,
    class_probs = c(
      dynamic_bspline = 0.20,
      constant = 0.40,
      zero = 0.40
    ),
    seed,
    scenario = "genotype_cosine_multipeak_dynamic_eqtl") {
  cell_label <- if (is.null(shape_cell_probs)) {
    "cells5"
  } else {
    paste0(
      "cells",
      length(shape_cell_probs),
      "_weights",
      numeric_vector_label(unname(shape_cell_probs))
    )
  }
  paste0(
    scenario,
    "_N", n_donors,
    "_J", n_variants,
    "_T", length(time_grid),
    "_pc", n_covariates,
    "_mix", numeric_vector_label(unname(class_probs)),
    "_noise", numeric_vector_label(expression_noise_sd),
    "_cosinev1",
    "_w", numeric_vector_label(width_half),
    "_rms", numeric_vector_label(target_centered_rms),
    "_", cell_label,
    "_seed", seed
  )
}

genotype_spiky_transient_output_stem <- function(n_donors,
                                                  n_variants,
                                                  time_grid,
                                                  n_covariates,
                                                  expression_noise_sd,
                                                  transient_bspline_df,
                                                  transient_bspline_degree,
                                                  dynamic_amplitude,
                                                  normalization = c(
                                                    "center_then_scale",
                                                    "legacy"
                                                  ),
                                                  class_probs = c(
                                                    dynamic_local_bspline_transient = 0.20,
                                                    constant = 0.40,
                                                    zero = 0.40
                                                  ),
                                                  seed,
                                                  scenario =
                                                    "genotype_spiky_transient_dynamic_eqtl") {
  normalization <- match.arg(normalization)
  normalization_label <- if (normalization == "center_then_scale") {
    "exact"
  } else {
    "legacy"
  }
  paste0(
    scenario,
    "_N", n_donors,
    "_J", n_variants,
    "_T", length(time_grid),
    "_pc", n_covariates,
    "_mix", numeric_vector_label(unname(class_probs)),
    "_noise", numeric_vector_label(expression_noise_sd),
    "_pairedv1",
    "_tdf", transient_bspline_df,
    "_deg", transient_bspline_degree,
    "_amp", numeric_vector_label(dynamic_amplitude),
    "_", normalization_label,
    "_seed", seed
  )
}

summarize_se_correction <- function(eqtl_summary) {
  ratio <- as.numeric(eqtl_summary$se / eqtl_summary$se_uncorrected)
  ratio <- ratio[is.finite(ratio)]
  data.frame(
    n = length(ratio),
    min = min(ratio),
    q25 = unname(stats::quantile(ratio, 0.25)),
    median = stats::median(ratio),
    mean = mean(ratio),
    q75 = unname(stats::quantile(ratio, 0.75)),
    max = max(ratio),
    stringsAsFactors = FALSE
  )
}

run_genotype_level_dynamic_eqtl_simulation <- function(G = NULL,
                                                       n_donors = 19,
                                                       n_variants = 1000,
                                                       time_grid = make_time_grid(),
                                                       n_covariates = 5,
                                                       covariates = NULL,
                                                       class_probs = c(
                                                         dynamic_bspline = 0.20,
                                                         constant = 0.40,
                                                         zero = 0.40
                                                       ),
                                                       maf_range = c(0.1, 0.5),
                                                       expression_noise_sd = 1,
                                                       covariate_effect_sd = 0.5,
                                                       intercept_sd = 0,
                                                       dynamic_amplitude = 2,
                                                       bspline_df = 6,
                                                       bspline_coefficient_sd = 1,
                                                       transient_bspline_df = 16,
                                                       transient_bspline_degree = 3,
                                                       constant_sd = 1,
                                                       dynamic_main_effect_sd = 0,
                                                       alpha = 0.05,
                                                       seed = 12345,
                                                       sigma_beta = 1,
                                                       estimate_sigma = FALSE,
                                                       sigma_beta_grid = exp(seq(log(0.05), log(5), length.out = 25)),
                                                       num_cores = 1,
                                                       num_basis = 20,
                                                       grid = default_revision_grid(),
                                                       penalty = 10,
                                                       exact_class_counts = TRUE,
                                                       apply_t_se_correction = TRUE,
                                                       scenario = "genotype_bspline_dynamic_eqtl",
                                                       output_dir = "output/revision_simulations",
                                                       save_outputs = TRUE,
                                                       verbose = FALSE,
                                                       effect_sim = NULL,
                                                       expression_sim = NULL,
                                                       expression_error_correlation = NULL,
                                                       linear_prior_mode = c(
                                                         "profile_single",
                                                         "mixture_grid"
                                                       ),
                                                       pred_step = 1) {
  linear_prior_mode <- match.arg(linear_prior_mode)
  if (!is.null(seed)) {
    set.seed(seed)
  }
  dirs <- revision_output_dirs(output_dir)

  if (is.null(G)) {
    genotype_sim <- simulate_genotype_matrix(
      n_donors = n_donors,
      n_variants = n_variants,
      maf_range = maf_range
    )
    G <- genotype_sim$G
    variant_info <- genotype_sim$variant_info
  } else {
    G <- validate_genotype_matrix(G)
    n_donors <- nrow(G)
    n_variants <- ncol(G)
    variant_info <- data.frame(
      variant_index = seq_len(n_variants),
      variant_id = colnames(G),
      maf = NA_real_,
      observed_maf = colMeans(G) / 2,
      genotype_sd = apply(G, 2, sd),
      stringsAsFactors = FALSE
    )
  }

  if (is.null(covariates)) {
    covariates <- simulate_covariate_matrix(
      n_donors = nrow(G),
      n_covariates = n_covariates
    )
  } else {
    covariates <- validate_covariate_matrix(covariates, nrow(G))
  }
  n_covariates_actual <- if (is.null(covariates)) 0 else ncol(covariates)

  if (is.null(effect_sim)) {
    effect_sim <- simulate_variant_effect_curves(
      n_variants = ncol(G),
      time_grid = time_grid,
      class_probs = class_probs,
      scenario = scenario,
      dynamic_amplitude = dynamic_amplitude,
      bspline_df = bspline_df,
      bspline_coefficient_sd = bspline_coefficient_sd,
      transient_bspline_df = transient_bspline_df,
      transient_bspline_degree = transient_bspline_degree,
      constant_sd = constant_sd,
      dynamic_main_effect_sd = dynamic_main_effect_sd,
      exact_class_counts = exact_class_counts,
      seed = if (is.null(seed)) NULL else
        revision_component_seeds(seed)[["functional_truth"]]
    )
  } else {
    if (!is.list(effect_sim) ||
        is.null(effect_sim$beta_matrix) ||
        is.null(effect_sim$unit_info)) {
      stop("effect_sim must contain beta_matrix and unit_info.")
    }
    effect_sim$beta_matrix <- as.matrix(effect_sim$beta_matrix)
    if (!identical(dim(effect_sim$beta_matrix), c(ncol(G), length(time_grid)))) {
      stop("effect_sim$beta_matrix dimensions must match variants and time points.")
    }
    if (nrow(effect_sim$unit_info) != ncol(G) ||
        !"effect_class" %in% colnames(effect_sim$unit_info)) {
      stop("effect_sim$unit_info must contain one effect_class row per variant.")
    }
  }
  rownames(effect_sim$beta_matrix) <- colnames(G)
  effect_sim$unit_info$variant_id <- colnames(G)

  if (is.null(expression_sim)) {
    expression_sim <- simulate_eqtl_expression_from_genotypes(
      G = G,
      beta_matrix = effect_sim$beta_matrix,
      time_grid = time_grid,
      covariates = covariates,
      expression_noise_sd = expression_noise_sd,
      covariate_effect_sd = covariate_effect_sd,
      intercept_sd = intercept_sd,
      expression_error_correlation = expression_error_correlation
    )
  } else {
    if (is.array(expression_sim)) {
      expression_sim <- list(expression = expression_sim)
    }
    if (!is.list(expression_sim) || is.null(expression_sim$expression)) {
      stop("expression_sim must be an expression array or a list containing expression.")
    }
    expected_expression_dim <- c(nrow(G), ncol(G), length(time_grid))
    if (!identical(dim(expression_sim$expression), expected_expression_dim)) {
      stop("expression_sim$expression dimensions must match donors, variants, and time points.")
    }
  }
  eqtl_summary <- estimate_eqtl_summaries_from_genotypes(
    G = G,
    expression = expression_sim$expression,
    covariates = covariates,
    apply_t_se_correction = apply_t_se_correction
  )

  datasets <- make_fash_datasets_from_eqtl_summary(
    beta_hat = eqtl_summary$beta_hat,
    se = eqtl_summary$se,
    true_beta = effect_sim$beta_matrix,
    time_grid = time_grid,
    unit_info = effect_sim$unit_info,
    scenario = scenario
  )
  unit_info <- attr(datasets, "unit_info")

  fash_fits <- fit_fash_for_revision(
    datasets = datasets,
    orders = 1,
    grid = grid,
    num_basis = num_basis,
    penalty = penalty,
    pred_step = pred_step,
    num_cores = num_cores,
    apply_bf = TRUE,
    verbose = verbose
  )
  if (linear_prior_mode == "mixture_grid") {
    simplified_fit <- fit_linear_mixture_fash(
      datasets = datasets,
      grid = grid,
      pred_step = pred_step,
      penalty = penalty
    )
    simplified_fit_bf <- BF_update_linear_mixture_fash(simplified_fit)
  } else {
    simplified_fit <- fit_simplified_fash(
      datasets = datasets,
      sigma_beta = sigma_beta,
      estimate_sigma = estimate_sigma,
      sigma_beta_grid = sigma_beta_grid,
      scale_time = TRUE
    )
    simplified_fit_bf <- BF_update_simplified_fash(simplified_fit)
  }

  result_table <- rbind(
    evaluate_lfdr_method(
      lfdr = get_fash_lfdr(fash_fits$fash_iwp1_raw),
      unit_info = unit_info,
      method = "FASH-IWP1-Raw",
      target = "dynamic",
      alpha = alpha
    ),
    evaluate_lfdr_method(
      lfdr = get_fash_lfdr(fash_fits$fash_iwp1_bf),
      unit_info = unit_info,
      method = "FASH-IWP1-BF",
      target = "dynamic",
      alpha = alpha
    ),
    evaluate_simplified_fash_fit(
      fit = simplified_fit,
      unit_info = unit_info,
      alpha = alpha,
      method = "FASH-linear-Raw"
    ),
    evaluate_simplified_fash_fit(
      fit = simplified_fit_bf,
      unit_info = unit_info,
      alpha = alpha,
      method = "FASH-linear-BF"
    )
  )
  direct_linear <- evaluate_direct_interaction_lrt(
    G = G,
    expression = expression_sim$expression,
    covariates = covariates,
    unit_info = unit_info,
    time_grid = time_grid,
    interaction_degree = 1,
    method = "Direct-linear-LRT",
    alpha = alpha
  )
  direct_quadratic <- evaluate_direct_interaction_lrt(
    G = G,
    expression = expression_sim$expression,
    covariates = covariates,
    unit_info = unit_info,
    time_grid = time_grid,
    interaction_degree = 2,
    method = "Direct-quadratic-LRT",
    alpha = alpha
  )
  result_table <- rbind(
    result_table,
    direct_linear$result,
    direct_quadratic$result
  )
  summary_table <- summarize_method_results(result_table)
  alpha_curve <- compute_alpha_curve(
    result_table = result_table,
    alpha_grid = seq(0, 0.20, by = 0.005)
  )
  se_correction_summary <- summarize_se_correction(eqtl_summary)

  settings <- list(
    n_donors = nrow(G),
    n_variants = ncol(G),
    time_grid = time_grid,
    n_covariates = n_covariates_actual,
    class_probs = class_probs,
    maf_range = maf_range,
    expression_noise_sd = expression_noise_sd,
    covariate_effect_sd = covariate_effect_sd,
    intercept_sd = intercept_sd,
    dynamic_amplitude = dynamic_amplitude,
    bspline_df = bspline_df,
    bspline_coefficient_sd = bspline_coefficient_sd,
    transient_bspline_df = transient_bspline_df,
    transient_bspline_degree = transient_bspline_degree,
    constant_sd = constant_sd,
    dynamic_main_effect_sd = dynamic_main_effect_sd,
    alpha = alpha,
    seed = seed,
    sigma_beta = sigma_beta,
    estimate_sigma = estimate_sigma,
    sigma_beta_grid = sigma_beta_grid,
    linear_prior_mode = linear_prior_mode,
    pred_step = pred_step,
    num_cores = num_cores,
    num_basis = num_basis,
    penalty = penalty,
    exact_class_counts = exact_class_counts,
    apply_t_se_correction = apply_t_se_correction,
    scenario = scenario,
    expression_error_correlation = if (is.null(expression_error_correlation)) {
      diag(length(time_grid))
    } else {
      validate_time_correlation(
        expression_error_correlation,
        n_time = length(time_grid)
      )
    }
  )

  out <- list(
    datasets = datasets,
    unit_info = unit_info,
    genotype = G,
    variant_info = variant_info,
    covariates = covariates,
    true_beta = effect_sim$beta_matrix,
    expression = expression_sim$expression,
    expression_simulation = expression_sim,
    eqtl_summary = eqtl_summary,
    se_correction_summary = se_correction_summary,
    fash_fits = fash_fits,
    simplified_fit = simplified_fit,
    simplified_fit_bf = simplified_fit_bf,
    direct_interaction_lrt = list(
      linear = direct_linear$lrt,
      quadratic = direct_quadratic$lrt,
      adjust_method = "BH"
    ),
    result_table = result_table,
    summary_table = summary_table,
    alpha_curve = alpha_curve,
    settings = settings
  )

  if (save_outputs) {
    stem <- genotype_dynamic_eqtl_output_stem(
      n_donors = nrow(G),
      n_variants = ncol(G),
      time_grid = time_grid,
      n_covariates = n_covariates_actual,
      expression_noise_sd = expression_noise_sd,
      class_probs = class_probs,
      seed = seed,
      scenario = scenario
    )
    saveRDS(out, file = file.path(dirs$raw, paste0(stem, ".rds")))
    write.csv(
      result_table,
      file = file.path(dirs$summary, paste0(stem, "_unit_results.csv")),
      row.names = FALSE
    )
    write.csv(
      summary_table,
      file = file.path(dirs$summary, paste0(stem, "_method_summary.csv")),
      row.names = FALSE
    )
    write.csv(
      alpha_curve,
      file = file.path(dirs$summary, paste0(stem, "_alpha_curve.csv")),
      row.names = FALSE
    )
    write.csv(
      se_correction_summary,
      file = file.path(dirs$summary, paste0(stem, "_se_correction_summary.csv")),
      row.names = FALSE
    )
  }

  out
}

run_genotype_level_bspline_eqtl_simulation <-
  run_genotype_level_dynamic_eqtl_simulation

interpolate_true_curve <- function(dat, n = 300) {
  x <- dat$x
  y <- dat$truef
  if (length(unique(x)) < 4 || n <= length(x)) {
    return(data.frame(x = x, truef = y))
  }

  interpolated <- stats::spline(
    x = x,
    y = y,
    n = n,
    method = "natural"
  )
  data.frame(x = interpolated$x, truef = interpolated$y)
}

plot_constant_dynamic_examples <- function(out,
                                           n_per_class = 3,
                                           file = NULL,
                                           seed = 1,
                                           true_curve_n = 300) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  datasets <- out$datasets
  unit_info <- out$unit_info
  null_ids <- unit_info$unit_id[unit_info$effect_class == "constant"]
  dynamic_class <- setdiff(unique(unit_info$effect_class), "constant")[1]
  dynamic_ids <- unit_info$unit_id[unit_info$effect_class == dynamic_class]
  selected_ids <- c(
    sample(null_ids, min(n_per_class, length(null_ids))),
    sample(dynamic_ids, min(n_per_class, length(dynamic_ids)))
  )

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    png(file, width = 1800, height = 950, res = 160)
    on.exit(dev.off(), add = TRUE)
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(2, n_per_class), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  for (unit_id in selected_ids) {
    dat <- datasets[[unit_id]]
    true_curve <- interpolate_true_curve(dat, n = true_curve_n)
    y_limits <- range(c(
      dat$y - 2 * dat$sd,
      dat$y + 2 * dat$sd,
      dat$truef,
      true_curve$truef
    ))
    plot(
      dat$x,
      dat$y,
      pch = 19,
      col = "black",
      ylim = y_limits,
      xlab = "Time",
      ylab = "Effect estimate",
      main = paste(dat$effect_class[1], unit_id, sep = "\n")
    )
    arrows(
      dat$x,
      dat$y - 2 * dat$sd,
      dat$x,
      dat$y + 2 * dat$sd,
      length = 0.03,
      angle = 90,
      code = 3,
      col = "gray40"
    )
    lines(true_curve$x, true_curve$truef, col = "red", lwd = 2)
    legend(
      "topleft",
      legend = c("Observed", "True effect"),
      pch = c(19, NA),
      lty = c(NA, 1),
      col = c("black", "red"),
      bty = "n",
      cex = 0.8
    )
  }
  mtext(
    paste0(
      if (!is.na(dynamic_class)) paste0("Constant null and ", dynamic_class, " dynamic examples; SE values = ") else
        "Constant null and dynamic examples; SE values = ",
      paste(out$settings$sd_values, collapse = ", ")
    ),
    outer = TRUE,
    cex = 1.1
  )
  invisible(selected_ids)
}

plot_constant_iwp2_examples <- plot_constant_dynamic_examples

plot_genotype_eqtl_examples <- function(out,
                                        classes = c("zero", "constant", "dynamic_bspline"),
                                        n_per_class = 3,
                                        file = NULL,
                                        seed = 1,
                                        true_curve_n = 300) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  datasets <- out$datasets
  unit_info <- out$unit_info
  classes <- classes[classes %in% unique(unit_info$effect_class)]
  if (length(classes) == 0) {
    stop("No requested classes are present in out$unit_info.")
  }

  selected_by_class <- lapply(classes, function(effect_class) {
    ids <- unit_info$unit_id[unit_info$effect_class == effect_class]
    sample(ids, min(n_per_class, length(ids)))
  })
  names(selected_by_class) <- classes

  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    png(
      file,
      width = 1800,
      height = 380 * length(classes),
      res = 160
    )
    on.exit(dev.off(), add = TRUE)
  }

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(
    mfrow = c(length(classes), n_per_class),
    mar = c(4, 4, 3, 1),
    oma = c(0, 0, 2, 0)
  )

  for (effect_class in classes) {
    selected_ids <- selected_by_class[[effect_class]]
    for (panel_index in seq_len(n_per_class)) {
      if (panel_index > length(selected_ids)) {
        plot.new()
        title(main = paste(effect_class, "not available"))
        next
      }
      unit_id <- selected_ids[panel_index]
      dat <- datasets[[unit_id]]
      true_curve <- interpolate_true_curve(dat, n = true_curve_n)
      y_limits <- range(c(
        dat$y - 2 * dat$sd,
        dat$y + 2 * dat$sd,
        dat$truef,
        true_curve$truef
      ))
      plot(
        dat$x,
        dat$y,
        pch = 19,
        col = "black",
        ylim = y_limits,
        xlab = "Time",
        ylab = "Estimated effect",
        main = paste(effect_class, unit_id, sep = "\n")
      )
      arrows(
        dat$x,
        dat$y - 2 * dat$sd,
        dat$x,
        dat$y + 2 * dat$sd,
        length = 0.03,
        angle = 90,
        code = 3,
        col = "gray40"
      )
      lines(true_curve$x, true_curve$truef, col = "red", lwd = 2)
      legend(
        "topleft",
        legend = c("Estimated", "True effect"),
        pch = c(19, NA),
        lty = c(NA, 1),
        col = c("black", "red"),
        bty = "n",
        cex = 0.75
      )
    }
  }

  mtext(
    paste0(
      "Genotype-level eQTL simulation examples; N = ",
      out$settings$n_donors,
      ", J = ",
      out$settings$n_variants,
      ", PCs = ",
      out$settings$n_covariates
    ),
    outer = TRUE,
    cex = 1.1
  )
  invisible(selected_by_class)
}

compute_dynamic_subgroup_alpha_curve <- function(
    result_table,
    unit_info,
    subgroup_var = "spike_count",
    alpha_grid = seq(0, 0.20, by = 0.005)) {
  required_result_columns <- c("unit_index", "method", "adjusted_score")
  if (!all(required_result_columns %in% names(result_table)) ||
      !all(c("unit_index", "effect_class", subgroup_var) %in% names(unit_info))) {
    stop("The result table or unit metadata is missing subgroup-curve columns.")
  }
  dynamic <- !is_null_class(unit_info$effect_class, target = "dynamic")
  subgroup_values <- sort(unique(unit_info[[subgroup_var]][dynamic]))
  subgroup_values <- subgroup_values[!is.na(subgroup_values)]
  methods <- unique(result_table$method)
  out <- do.call(rbind, lapply(methods, function(method) {
    method_rows <- result_table[result_table$method == method, ]
    do.call(rbind, lapply(alpha_grid, function(alpha) {
      selected <- method_rows$unit_index[
        is.finite(method_rows$adjusted_score) &
          method_rows$adjusted_score <= alpha
      ]
      do.call(rbind, lapply(subgroup_values, function(subgroup_value) {
        subgroup <- dynamic & unit_info[[subgroup_var]] == subgroup_value
        n_dynamic <- sum(subgroup)
        true_positives <- sum(unit_info$unit_index[subgroup] %in% selected)
        data.frame(
          scenario = unit_info$scenario[which(subgroup)[1]],
          target = "dynamic",
          method = method,
          alpha = alpha,
          subgroup_var = subgroup_var,
          subgroup_value = as.character(subgroup_value),
          n_dynamic = n_dynamic,
          true_positives = true_positives,
          power = true_positives / n_dynamic,
          stringsAsFactors = FALSE
        )
      }))
    }))
  }))
  rownames(out) <- NULL
  out
}

plot_cosine_peak_cell_examples <- function(out,
                                           file = NULL,
                                           seed = 1,
                                           ncol = 3,
                                           include_time_groups = FALSE,
                                           true_color = "#D55E00") {
  required_fields <- c(
    "datasets", "unit_info", "true_beta_evaluation", "evaluation_grid"
  )
  if (!all(required_fields %in% names(out))) {
    stop("out is missing high-resolution cosine example fields.")
  }
  dynamic <- !is_null_class(out$unit_info$effect_class, target = "dynamic")
  metadata_columns <- c("cell_id", "spike_count", "sign_pattern")
  if (include_time_groups) {
    if (!"time_group" %in% names(out$unit_info) ||
        any(is.na(out$unit_info$time_group[dynamic]))) {
      stop("Timed cosine examples require non-missing time_group metadata.")
    }
    metadata_columns <- c(metadata_columns, "time_group")
  }
  cell_metadata <- unique(out$unit_info[
    dynamic,
    metadata_columns,
    drop = FALSE
  ])
  cell_metadata <- cell_metadata[!is.na(cell_metadata$cell_id), , drop = FALSE]
  sign_order <- match(
    cell_metadata$sign_pattern,
    c("single", "same-sign", "alternating-sign")
  )
  if (include_time_groups) {
    time_order <- match(
      cell_metadata$time_group,
      c("early", "middle", "late")
    )
    cell_metadata <- cell_metadata[
      order(time_order, cell_metadata$spike_count, sign_order),
      ,
      drop = FALSE
    ]
  } else {
    cell_metadata <- cell_metadata[
      order(cell_metadata$spike_count, sign_order),
      ,
      drop = FALSE
    ]
  }
  if (nrow(cell_metadata) == 0) {
    stop("No dynamic cosine shape cells are available.")
  }
  if (!is.null(seed)) set.seed(seed)
  selected <- vapply(seq_len(nrow(cell_metadata)), function(i) {
    candidates <- which(
      dynamic & out$unit_info$cell_id == cell_metadata$cell_id[i]
    )
    if (include_time_groups) {
      candidates <- candidates[
        out$unit_info$time_group[candidates] == cell_metadata$time_group[i]
      ]
    }
    candidates[sample.int(length(candidates), size = 1)]
  }, integer(1))

  nrow <- ceiling(length(selected) / ncol)
  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    png(file, width = 650 * ncol, height = 470 * nrow, res = 170)
    on.exit(dev.off(), add = TRUE)
  }
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mfrow = c(nrow, ncol), mar = c(3.8, 4.0, 3.1, 0.8))

  for (index in selected) {
    dat <- out$datasets[[out$unit_info$unit_id[index]]]
    true_curve <- out$true_beta_evaluation[index, ]
    y_limits <- range(
      dat$y - 2 * dat$sd,
      dat$y + 2 * dat$sd,
      true_curve
    )
    shape_title <- paste0(
      out$unit_info$spike_count[index],
      if (out$unit_info$spike_count[index] == 1L) " peak" else " peaks",
      if (out$unit_info$spike_count[index] > 1L) {
        paste0("; ", gsub("-", " ", out$unit_info$sign_pattern[index]))
      } else {
        ""
      }
    )
    title <- if (include_time_groups) {
      paste0(
        tools::toTitleCase(out$unit_info$time_group[index]),
        ": ",
        shape_title
      )
    } else {
      shape_title
    }
    plot(
      dat$x,
      dat$y,
      pch = 19,
      col = "black",
      ylim = y_limits,
      xlab = "Time",
      ylab = "Estimated effect",
      main = title
    )
    arrows(
      dat$x,
      dat$y - 2 * dat$sd,
      dat$x,
      dat$y + 2 * dat$sd,
      length = 0.03,
      angle = 90,
      code = 3,
      col = "gray45"
    )
    lines(
      out$evaluation_grid,
      true_curve,
      col = true_color,
      lwd = 2.2
    )
    legend(
      "topleft",
      legend = c("Estimated", "True effect"),
      pch = c(19, NA),
      lty = c(NA, 1),
      col = c("black", true_color),
      bty = "n",
      cex = 0.75
    )
  }
  if (length(selected) < nrow * ncol) {
    for (i in seq_len(nrow * ncol - length(selected))) plot.new()
  }
  invisible(out$unit_info[selected, ])
}

run_revision_smoke_test <- function() {
  run_revision_simulation(
    J = 30,
    class_probs = c(constant = 0.60, linear = 0.20, nonlinear_iwp2 = 0.20),
    scenario = "smoke_test",
    alpha = 0.05,
    seed = 1,
    num_cores = 1,
    num_basis = 8,
    grid = sort(c(0, exp(-0.5 * seq(0, 4, by = 1)))),
    penalty = 1,
    save_outputs = FALSE,
    verbose = FALSE
  )
}

if (identical(Sys.getenv("FASH_REVISION_RUN_SMOKE_TEST"), "true")) {
  smoke <- run_revision_smoke_test()
  print(smoke$summary_table)
}
