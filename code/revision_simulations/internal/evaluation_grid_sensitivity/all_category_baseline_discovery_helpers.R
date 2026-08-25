# Shared helpers for the all-category baseline-discovery sensitivity analysis.

functional_window_mask <- function(evaluation_grid,
                                   category,
                                   middle_definition = c(
                                     "open_3_12",
                                     "closed_4_11"
                                   ),
                                   tolerance = sqrt(.Machine$double.eps)) {
  middle_definition <- match.arg(middle_definition)
  if (!is.numeric(evaluation_grid) || length(evaluation_grid) < 2L ||
      any(!is.finite(evaluation_grid)) || is.unsorted(evaluation_grid) ||
      !is.character(category) || length(category) != 1L ||
      !category %in% c("early", "middle", "late") ||
      !is.finite(tolerance) || tolerance < 0) {
    stop("Invalid evaluation grid or location category.")
  }

  switch(
    category,
    early = evaluation_grid <= 3 + tolerance,
    middle = if (middle_definition == "open_3_12") {
      evaluation_grid > 3 + tolerance &
        evaluation_grid < 12 - tolerance
    } else {
      evaluation_grid >= 4 - tolerance &
        evaluation_grid <= 11 + tolerance
    },
    late = evaluation_grid >= 12 - tolerance
  )
}

functional_failure_indicator <- function(samples,
                                         evaluation_grid,
                                         category,
                                         switch_threshold = 0.25,
                                         middle_definition = c(
                                           "open_3_12",
                                           "closed_4_11"
                                         )) {
  middle_definition <- match.arg(middle_definition)
  if (!is.matrix(samples) || nrow(samples) != length(evaluation_grid) ||
      ncol(samples) < 1L || any(!is.finite(samples)) ||
      !is.character(category) || length(category) != 1L ||
      !category %in% c("early", "middle", "late", "switch") ||
      !is.finite(switch_threshold) || switch_threshold <= 0) {
    stop("Invalid posterior samples or functional specification.")
  }

  if (category == "switch") {
    positive_maximum <- pmax(matrixStats::colMaxs(samples), 0)
    negative_maximum <- pmax(-matrixStats::colMins(samples), 0)
    switch_statistic <-
      pmin(positive_maximum, negative_maximum) - switch_threshold
    return(as.integer(switch_statistic <= 0))
  }

  inside_window <- functional_window_mask(
    evaluation_grid,
    category,
    middle_definition = middle_definition
  )
  if (!any(inside_window) || !any(!inside_window)) {
    stop("The evaluation grid does not support the requested functional.")
  }

  absolute_samples <- abs(samples)
  location_statistic <-
    matrixStats::colMaxs(
      absolute_samples,
      rows = which(inside_window)
    ) -
    matrixStats::colMaxs(
      absolute_samples,
      rows = which(!inside_window)
    )
  as.integer(location_statistic <= 0)
}

functional_lfsr <- function(samples,
                            evaluation_grid,
                            category,
                            switch_threshold = 0.25,
                            middle_definition = c(
                              "open_3_12",
                              "closed_4_11"
                            )) {
  middle_definition <- match.arg(middle_definition)
  mean(functional_failure_indicator(
    samples = samples,
    evaluation_grid = evaluation_grid,
    category = category,
    switch_threshold = switch_threshold,
    middle_definition = middle_definition
  ))
}

shuffle_posterior_draws <- function(samples, seed) {
  if (!is.matrix(samples) || ncol(samples) < 2L ||
      any(!is.finite(samples)) || length(seed) != 1L ||
      is.na(seed) || !is.finite(seed)) {
    stop("Invalid posterior sample matrix or shuffle seed.")
  }
  set.seed(as.integer(seed))
  samples[, sample.int(ncol(samples)), drop = FALSE]
}

set_jaccard <- function(left, right) {
  left <- unique(left)
  right <- unique(right)
  union_size <- length(union(left, right))
  if (union_size == 0L) return(1)
  length(intersect(left, right)) / union_size
}

summarize_reported_pair_comparison <- function(reference,
                                               comparison,
                                               category,
                                               comparison_label) {
  required_columns <- c("pair_id", "gene_id", "lfsr")
  if (!is.data.frame(reference) || !is.data.frame(comparison) ||
      !all(required_columns %in% names(reference)) ||
      !all(required_columns %in% names(comparison)) ||
      anyDuplicated(reference$pair_id) || anyDuplicated(comparison$pair_id) ||
      !setequal(reference$pair_id, comparison$pair_id)) {
    stop("The reported-pair comparison is incomplete or duplicated.")
  }

  comparison <- comparison[
    match(reference$pair_id, comparison$pair_id),
    ,
    drop = FALSE
  ]
  if (anyNA(comparison$pair_id) ||
      any(!is.finite(reference$lfsr)) ||
      any(!is.finite(comparison$lfsr)) ||
      any(reference$lfsr < 0 | reference$lfsr > 1) ||
      any(comparison$lfsr < 0 | comparison$lfsr > 1)) {
    stop("The paired LFSR values are invalid.")
  }

  difference <- comparison$lfsr - reference$lfsr
  absolute_difference <- abs(difference)
  reference_below <- reference$lfsr <= 0.05
  comparison_below <- comparison$lfsr <= 0.05
  numerical_tolerance <- sqrt(.Machine$double.eps)

  data.frame(
    category = category,
    comparison = comparison_label,
    n_pairs = nrow(reference),
    n_genes = length(unique(reference$gene_id)),
    mean_signed_change = mean(difference),
    mean_absolute_change = mean(absolute_difference),
    median_absolute_change = median(absolute_difference),
    q90_absolute_change = unname(quantile(absolute_difference, 0.90)),
    maximum_absolute_change = max(absolute_difference),
    spearman_correlation = suppressWarnings(cor(
      reference$lfsr,
      comparison$lfsr,
      method = "spearman"
    )),
    fraction_within_0p005 = mean(
      absolute_difference <= 0.005 + numerical_tolerance
    ),
    fraction_within_0p01 = mean(
      absolute_difference <= 0.01 + numerical_tolerance
    ),
    n_lfsr_below_0p05_reference = sum(reference_below),
    n_lfsr_below_0p05_comparison = sum(comparison_below),
    lfsr_0p05_classification_agreement = mean(
      reference_below == comparison_below
    ),
    lfsr_0p05_jaccard = set_jaccard(
      reference$pair_id[reference_below],
      comparison$pair_id[comparison_below]
    ),
    stringsAsFactors = FALSE
  )
}
