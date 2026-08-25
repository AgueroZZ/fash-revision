# Helpers for the small-M matched-null empirical-Bayes sensitivity analysis.

safe_scalar_ratio <- function(numerator, denominator) {
  if (length(numerator) != 1L || length(denominator) != 1L ||
      !is.finite(numerator) || !is.finite(denominator) || denominator == 0) {
    return(NA_real_)
  }
  numerator / denominator
}

make_nested_null_subsets <- function(n_null,
                                     ratio_to_size,
                                     n_replicates,
                                     seed) {
  n_null <- as.integer(n_null)
  ratio_names <- names(ratio_to_size)
  ratio_to_size <- as.integer(ratio_to_size)
  ratios <- suppressWarnings(as.numeric(ratio_names))
  n_replicates <- as.integer(n_replicates)
  seed <- as.integer(seed)
  if (length(n_null) != 1L || is.na(n_null) || n_null < 1L ||
      length(ratio_to_size) < 1L || anyNA(ratio_to_size) ||
      any(ratio_to_size < 1L) || any(ratio_to_size > n_null) ||
      is.unsorted(ratio_to_size, strictly = TRUE) ||
      length(ratios) != length(ratio_to_size) || any(!is.finite(ratios)) ||
      any(ratios <= 0) || is.unsorted(ratios, strictly = TRUE) ||
      length(n_replicates) != 1L || is.na(n_replicates) ||
      n_replicates < 1L || length(seed) != 1L || is.na(seed)) {
    stop("Invalid nested null-subset specification.")
  }

  set.seed(seed)
  rows <- vector("list", n_replicates * length(ratio_to_size))
  row_index <- 0L
  for (replicate_id in seq_len(n_replicates)) {
    ordering <- sample.int(n_null, size = n_null, replace = FALSE)
    for (ratio_index in seq_along(ratio_to_size)) {
      row_index <- row_index + 1L
      subset_size <- ratio_to_size[ratio_index]
      rows[[row_index]] <- data.frame(
        replicate_id = replicate_id,
        m_ratio = ratios[ratio_index],
        m_size = subset_size,
        subset_position = seq_len(subset_size),
        null_index = ordering[seq_len(subset_size)],
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

summarize_small_m_calibration <- function(lfdr,
                                          group,
                                          pi0_merged,
                                          fit_stage,
                                          m_ratio,
                                          replicate_id,
                                          alpha,
                                          selected_indices = NULL) {
  lfdr <- as.numeric(lfdr)
  group <- as.character(group)
  pi0_merged <- as.numeric(pi0_merged)
  alpha <- as.numeric(alpha)
  if (length(lfdr) != length(group) || length(lfdr) < 2L ||
      any(!is.finite(lfdr)) || any(lfdr < 0 | lfdr > 1) ||
      !setequal(unique(group), c("target", "permuted_null")) ||
      length(pi0_merged) != 1L || !is.finite(pi0_merged) ||
      pi0_merged < 0 || pi0_merged > 1 ||
      length(alpha) != 1L || !is.finite(alpha) || alpha < 0 || alpha > 1) {
    stop("Invalid small-M calibration inputs.")
  }

  if (is.null(selected_indices)) {
    selected_indices <- cumulative_lfdr_calls(lfdr, alpha = alpha)
  } else {
    selected_indices <- as.integer(selected_indices)
    if (anyNA(selected_indices) || any(selected_indices < 1L) ||
        any(selected_indices > length(lfdr)) ||
        anyDuplicated(selected_indices)) {
      stop("Invalid precomputed selected indices.")
    }
  }
  selected <- seq_along(lfdr) %in% selected_indices
  target <- group == "target"
  permuted_null <- group == "permuted_null"
  n_target <- sum(target)
  n_null <- sum(permuted_null)
  n_merged <- length(lfdr)
  target_calls <- sum(selected & target)
  null_calls <- sum(selected & permuted_null)
  merged_calls <- sum(selected)
  target_call_rate <- target_calls / n_target
  null_call_rate <- null_calls / n_null
  merged_call_rate <- merged_calls / n_merged
  pi0_target_unbounded <- (
    n_merged * pi0_merged - n_null
  ) / n_target
  pi0_target_valid <- pi0_target_unbounded >= 0 &&
    pi0_target_unbounded <= 1

  data.frame(
    replicate_id = as.integer(replicate_id),
    m_ratio = as.numeric(m_ratio),
    m_size = n_null,
    fit_stage = as.character(fit_stage),
    nominal_alpha = alpha,
    n_target = n_target,
    n_merged = n_merged,
    target_calls = target_calls,
    permuted_null_calls = null_calls,
    merged_calls = merged_calls,
    target_call_rate = target_call_rate,
    permuted_null_call_rate = null_call_rate,
    merged_call_rate = merged_call_rate,
    pi0_merged = pi0_merged,
    pi0_target_unbounded = pi0_target_unbounded,
    pi0_target_valid = pi0_target_valid,
    pi0_target_bounded = min(1, max(0, pi0_target_unbounded)),
    target_pi0_plugin_fdr = if (pi0_target_valid) {
      safe_scalar_ratio(
        pi0_target_unbounded * null_call_rate,
        target_call_rate
      )
    } else {
      NA_real_
    },
    merged_pi0_plugin_fdr = safe_scalar_ratio(
      pi0_merged * null_call_rate,
      merged_call_rate
    ),
    known_null_discovery_fraction = safe_scalar_ratio(
      null_calls,
      merged_calls
    ),
    stringsAsFactors = FALSE
  )
}

summarize_small_m_curve <- function(lfdr,
                                    group,
                                    pi0_merged,
                                    fit_stage,
                                    m_ratio,
                                    replicate_id,
                                    alpha_grid) {
  ordering <- order(lfdr, method = "radix")
  cumulative_mean <- cumsum(lfdr[ordering]) / seq_along(ordering)
  rows <- lapply(alpha_grid, function(alpha) {
    accepted <- which(cumulative_mean <= alpha)
    selected_indices <- if (length(accepted) == 0L) {
      integer()
    } else {
      ordering[seq_len(max(accepted))]
    }
    summarize_small_m_calibration(
      lfdr = lfdr,
      group = group,
      pi0_merged = pi0_merged,
      fit_stage = fit_stage,
      m_ratio = m_ratio,
      replicate_id = replicate_id,
      alpha = alpha,
      selected_indices = selected_indices
    )
  })
  do.call(rbind, rows)
}
