# Define the canonical IWP1 temporal truth-category mixture used by R3.

r3_iwp1_temporal_category_probs <- function() {
  c(early = 0.29, middle = 0.42, late = 0.29)
}

r3_iwp1_temporal_mixture_contract <- function(n_dynamic = 1272L) {
  if (length(n_dynamic) != 1L || is.na(n_dynamic) ||
      n_dynamic != round(n_dynamic) || n_dynamic < 6L) {
    stop("n_dynamic must be one integer of at least six.", call. = FALSE)
  }
  if (!exists("exact_proportional_counts", mode = "function", inherits = TRUE) ||
      !exists(
        "exact_temporal_truth_group_counts",
        mode = "function",
        inherits = TRUE
      )) {
    stop(
      "Source simulation_functions.R before the R3 temporal-mixture contract.",
      call. = FALSE
    )
  }

  n_dynamic <- as.integer(n_dynamic)
  probabilities <- r3_iwp1_temporal_category_probs()
  truth_group_levels <- c(
    "early / switch",
    "early / non-switch",
    "middle / switch",
    "middle / non-switch",
    "late / switch",
    "late / non-switch"
  )
  temporal_counts <- exact_proportional_counts(n_dynamic, probabilities)
  truth_group_counts <- exact_temporal_truth_group_counts(
    n_dynamic,
    temporal_category_probs = probabilities
  )[truth_group_levels]
  truth_group_counts <- stats::setNames(
    as.integer(truth_group_counts), truth_group_levels
  )
  switch_status_counts <- c(
    switch = sum(truth_group_counts[grepl(" / switch$", truth_group_levels)]),
    `non-switch` = sum(
      truth_group_counts[grepl(" / non-switch$", truth_group_levels)]
    )
  )
  switch_status_counts <- as.integer(switch_status_counts) |>
    stats::setNames(c("switch", "non-switch"))

  list(
    name = "canonical IWP1 temporal-category mixture",
    temporal_category_probs = probabilities,
    temporal_category_counts = temporal_counts,
    truth_group_levels = truth_group_levels,
    truth_group_counts = truth_group_counts,
    switch_status_counts = switch_status_counts,
    n_dynamic = n_dynamic
  )
}
