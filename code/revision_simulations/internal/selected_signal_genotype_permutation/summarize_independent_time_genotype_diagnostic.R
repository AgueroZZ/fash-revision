#!/usr/bin/env Rscript

# Summarize temporal dependence and signal leakage in one completed diagnostic.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(args, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(args[equals_hit[1L]], nchar(equals_prefix) + 1L))
  }
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1L] == length(args)) {
    return(default)
  }
  args[hit[1L] + 1L]
}

cumulative_lfdr_calls <- function(lfdr, alpha) {
  ordering <- order(lfdr, method = "radix")
  cumulative_mean <- cumsum(lfdr[ordering]) / seq_along(ordering)
  accepted <- which(cumulative_mean <= alpha)
  if (length(accepted) == 0L) {
    return(integer())
  }
  ordering[accepted]
}

workflowr_root <- find_workflowr_root()
output_id <- get_arg("--output-id", "")
if (length(output_id) != 1L || !nzchar(output_id) ||
    grepl("/", output_id, fixed = TRUE)) {
  stop("The requested output ID is invalid.")
}
output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  output_id
)
required_paths <- file.path(
  output_directory,
  c(
    "effect_estimates.rds",
    "merged_fash_fit.rds",
    "selection.csv",
    "unit_lfdr.csv"
  )
)
if (any(!file.exists(required_paths))) {
  stop("The completed diagnostic bundle is incomplete.")
}

effect_estimates <- readRDS(required_paths[1L])
fit_bundle <- readRDS(required_paths[2L])
selection <- utils::read.csv(required_paths[3L], stringsAsFactors = FALSE)
unit_lfdr <- utils::read.csv(required_paths[4L], stringsAsFactors = FALSE)
configuration <- fit_bundle$configuration
if (!identical(
      configuration$permutation_method,
      "genotype_label_independent_time"
    ) || nrow(selection) != configuration$n_target_units) {
  stop("The requested bundle is not an independent-time genotype diagnostic.")
}

observed_z <- effect_estimates$target_beta /
  effect_estimates$target_adjusted_se
null_z <- effect_estimates$permuted_beta /
  effect_estimates$permuted_adjusted_se
if (!identical(dim(observed_z), dim(null_z)) ||
    any(!is.finite(observed_z)) || any(!is.finite(null_z))) {
  stop("The observed or permuted z-score matrix is invalid.")
}

null_correlation <- stats::cor(null_z)
off_diagonal <- row(null_correlation) != col(null_correlation)
adjacent <- abs(row(null_correlation) - col(null_correlation)) == 1L
temporal_summary <- data.frame(
  metric = c(
    "mean_off_diagonal_null_z_correlation",
    "mean_adjacent_null_z_correlation",
    "maximum_absolute_off_diagonal_null_z_correlation"
  ),
  value = c(
    mean(null_correlation[off_diagonal]),
    mean(null_correlation[adjacent]),
    max(abs(null_correlation[off_diagonal]))
  ),
  stringsAsFactors = FALSE
)

bf_rows <- unit_lfdr$fit_stage == "BF-adjusted"
bf_data <- unit_lfdr[bf_rows, , drop = FALSE]
null_rows <- bf_data$group == "permuted_null"
null_data <- bf_data[null_rows, , drop = FALSE]
null_order <- match(selection$pair_key, null_data$source_pair_key)
if (anyNA(null_order) || anyDuplicated(null_order)) {
  stop("The permuted-null lfdr values do not align with the selection.")
}
null_lfdr <- null_data$lfdr[null_order]
observed_max_abs_z <- apply(abs(observed_z), 1L, max)
observed_rms_z <- sqrt(rowMeans(observed_z^2))
null_max_abs_z <- apply(abs(null_z), 1L, max)
null_evidence <- -log10(pmax(null_lfdr, .Machine$double.xmin))
source_evidence <- -log10(pmax(
  selection$source_bf_lfdr,
  .Machine$double.xmin
))
time_specific_observed_null_correlations <- vapply(
  seq_len(ncol(observed_z)),
  function(time_index) {
    stats::cor(
      observed_z[, time_index],
      null_z[, time_index],
      method = "spearman"
    )
  },
  numeric(1)
)

selected_indices <- cumulative_lfdr_calls(bf_data$lfdr, alpha = 0.05)
called_null_keys <- bf_data$unit_key[
  selected_indices[bf_data$group[selected_indices] == "permuted_null"]
]
null_called <- null_data$unit_key[null_order] %in% called_null_keys
safe_median <- function(x) {
  if (length(x) == 0L) NA_real_ else stats::median(x)
}
leakage_summary <- data.frame(
  metric = c(
    "spearman_null_evidence_vs_source_evidence",
    "spearman_null_evidence_vs_observed_max_abs_z",
    "spearman_null_evidence_vs_observed_rms_z",
    "spearman_null_evidence_vs_null_max_abs_z",
    "mean_time_specific_observed_vs_null_z_spearman",
    "maximum_absolute_time_specific_observed_vs_null_z_spearman",
    "called_null_units_at_alpha_0.05",
    "called_null_median_source_lfdr",
    "uncalled_null_median_source_lfdr",
    "called_null_median_observed_max_abs_z",
    "uncalled_null_median_observed_max_abs_z"
  ),
  value = c(
    stats::cor(null_evidence, source_evidence, method = "spearman"),
    stats::cor(null_evidence, observed_max_abs_z, method = "spearman"),
    stats::cor(null_evidence, observed_rms_z, method = "spearman"),
    stats::cor(null_evidence, null_max_abs_z, method = "spearman"),
    mean(time_specific_observed_null_correlations),
    max(abs(time_specific_observed_null_correlations)),
    sum(null_called),
    safe_median(selection$source_bf_lfdr[null_called]),
    safe_median(selection$source_bf_lfdr[!null_called]),
    safe_median(observed_max_abs_z[null_called]),
    safe_median(observed_max_abs_z[!null_called])
  ),
  stringsAsFactors = FALSE
)

utils::write.csv(
  temporal_summary,
  file.path(output_directory, "temporal_correlation_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  leakage_summary,
  file.path(output_directory, "signal_leakage_diagnostics.csv"),
  row.names = FALSE
)

cat("Independent-time genotype diagnostic summary created.\n")
print(temporal_summary, row.names = FALSE)
print(leakage_summary, row.names = FALSE)
