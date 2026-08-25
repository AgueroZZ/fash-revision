# Shared helpers for the open-Middle-window numerical sensitivity analysis.

middle_open_window_mask <- function(time,
                                    lower_boundary = 3,
                                    upper_boundary = 12,
                                    tolerance = sqrt(.Machine$double.eps)) {
  if (!is.numeric(time) || any(!is.finite(time)) ||
      !is.finite(lower_boundary) || !is.finite(upper_boundary) ||
      lower_boundary >= upper_boundary ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("Invalid time grid or Middle-window boundary.")
  }

  time > lower_boundary + tolerance & time < upper_boundary - tolerance
}

middle_open_window_lfsr <- function(samples, evaluation_grid) {
  if (!is.matrix(samples) || nrow(samples) != length(evaluation_grid) ||
      any(!is.finite(samples))) {
    stop("Posterior samples do not match the evaluation grid.")
  }

  inside <- middle_open_window_mask(evaluation_grid)
  if (!any(inside) || !any(!inside)) {
    stop("The evaluation grid does not support the open Middle functional.")
  }

  absolute_samples <- abs(samples)
  statistic <-
    matrixStats::colMaxs(absolute_samples, rows = which(inside)) -
    matrixStats::colMaxs(absolute_samples, rows = which(!inside))
  mean(statistic <= 0)
}

deterministic_cumulative_fsr <- function(lfsr, pair_id) {
  if (!is.numeric(lfsr) || length(lfsr) != length(pair_id) ||
      any(!is.finite(lfsr)) || any(lfsr < 0 | lfsr > 1) ||
      anyNA(pair_id) || anyDuplicated(pair_id)) {
    stop("Invalid LFSR vector or pair identifiers.")
  }

  ordering <- order(lfsr, pair_id)
  ranked_cumulative_fsr <- cumsum(lfsr[ordering]) / seq_along(ordering)
  cumulative_fsr <- numeric(length(lfsr))
  cumulative_fsr[ordering] <- ranked_cumulative_fsr
  cumulative_fsr
}

set_jaccard <- function(left, right) {
  left <- unique(left)
  right <- unique(right)
  union_size <- length(union(left, right))
  if (union_size == 0L) return(1)
  length(intersect(left, right)) / union_size
}
