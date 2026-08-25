#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  helper <- file.path(
    "code", "revision_simulations", "internal",
    "evaluation_grid_sensitivity",
    "all_category_baseline_discovery_helpers.R"
  )
  if (file.exists(helper)) return(".")
  if (file.exists(file.path("coderepo-local", helper))) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "evaluation_grid_sensitivity",
  "all_category_baseline_discovery_helpers.R"
))

if (!requireNamespace("matrixStats", quietly = TRUE)) {
  stop("The matrixStats package is required.")
}

evaluation_grid <- c(0, 3, 6, 12, 15)
stopifnot(
  identical(
    functional_window_mask(evaluation_grid, "early"),
    c(TRUE, TRUE, FALSE, FALSE, FALSE)
  ),
  identical(
    functional_window_mask(evaluation_grid, "middle"),
    c(FALSE, FALSE, TRUE, FALSE, FALSE)
  ),
  identical(
    functional_window_mask(evaluation_grid, "late"),
    c(FALSE, FALSE, FALSE, TRUE, TRUE)
  )
)

middle_boundary_grid <- c(3, 3.5, 4, 6, 11, 11.5, 12)
stopifnot(
  identical(
    functional_window_mask(
      middle_boundary_grid,
      "middle",
      middle_definition = "open_3_12"
    ),
    c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE)
  ),
  identical(
    functional_window_mask(
      middle_boundary_grid,
      "middle",
      middle_definition = "closed_4_11"
    ),
    c(FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE)
  )
)

samples <- cbind(
  early = c(5, 4, 1, 1, 0),
  middle = c(1, 0, 5, 1, 0),
  late = c(0, 1, 1, 5, 4),
  switch = c(2, 0, -2, 0, 0),
  no_switch = c(0.1, 0, 0, 0, 0)
)
stopifnot(
  functional_failure_indicator(
    samples[, "early", drop = FALSE],
    evaluation_grid,
    "early"
  ) == 0L,
  functional_failure_indicator(
    samples[, "middle", drop = FALSE],
    evaluation_grid,
    "middle"
  ) == 0L,
  functional_failure_indicator(
    samples[, "late", drop = FALSE],
    evaluation_grid,
    "late"
  ) == 0L,
  functional_failure_indicator(
    samples[, "switch", drop = FALSE],
    evaluation_grid,
    "switch"
  ) == 0L,
  functional_failure_indicator(
    samples[, "no_switch", drop = FALSE],
    evaluation_grid,
    "switch"
  ) == 1L
)

ordered_draws <- matrix(seq_len(30L), nrow = 5L, ncol = 6L)
shuffled_draws <- shuffle_posterior_draws(ordered_draws, seed = 20260820L)
stopifnot(
  identical(
    shuffled_draws,
    shuffle_posterior_draws(ordered_draws, seed = 20260820L)
  ),
  !identical(shuffled_draws, ordered_draws),
  identical(
    sort(apply(shuffled_draws, 2L, paste, collapse = ":")),
    sort(apply(ordered_draws, 2L, paste, collapse = ":"))
  )
)

reference <- data.frame(
  pair_id = c("gene1_rs1", "gene2_rs2"),
  gene_id = c("gene1", "gene2"),
  lfsr = c(0.01, 0.04),
  stringsAsFactors = FALSE
)
comparison <- data.frame(
  pair_id = rev(reference$pair_id),
  gene_id = rev(reference$gene_id),
  lfsr = c(0.03, 0.02),
  stringsAsFactors = FALSE
)
summary_row <- summarize_reported_pair_comparison(
  reference,
  comparison,
  category = "early",
  comparison_label = "synthetic"
)
stopifnot(
  summary_row$n_pairs == 2L,
  summary_row$n_genes == 2L,
  abs(summary_row$mean_signed_change) < 1e-12,
  abs(summary_row$mean_absolute_change - 0.01) < 1e-12,
  abs(summary_row$maximum_absolute_change - 0.01) < 1e-12,
  summary_row$fraction_within_0p01 == 1,
  summary_row$lfsr_0p05_classification_agreement == 1
)

expect_error <- function(expression) {
  inherits(try(force(expression), silent = TRUE), "try-error")
}
stopifnot(
  expect_error(functional_window_mask(c(0, NA_real_), "early")),
  expect_error(functional_window_mask(evaluation_grid, "unknown")),
  expect_error(functional_window_mask(
    evaluation_grid,
    "middle",
    middle_definition = "unknown"
  )),
  expect_error(functional_lfsr(samples[-1, , drop = FALSE], evaluation_grid, "early")),
  expect_error(shuffle_posterior_draws(ordered_draws, seed = NA_integer_)),
  expect_error(summarize_reported_pair_comparison(
    reference,
    comparison[-1, , drop = FALSE],
    "early",
    "incomplete"
  ))
)

message("All all-category baseline-discovery sensitivity helper tests passed.")
