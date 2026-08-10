#!/usr/bin/env Rscript

# Separate likelihood changes from empirical-Bayes prior re-estimation.

find_workflowr_root <- function() {
  if (file.exists("output/revision_simulations/internal/correlated_likelihood_sensitivity/fit_bundle.rds")) {
    return(normalizePath("."))
  }
  if (file.exists(
    "coderepo-local/output/revision_simulations/internal/correlated_likelihood_sensitivity/fit_bundle.rds"
  )) {
    return(normalizePath("coderepo-local"))
  }
  stop("Could not find the workflowr repository root.")
}

log_sum_exp <- function(x) {
  maximum <- max(x)
  maximum + log(sum(exp(x - maximum)))
}

count_cumulative_lfdr_calls <- function(lfdr, alpha = 0.05) {
  ordered_lfdr <- sort(lfdr)
  eligible <- which(cumsum(ordered_lfdr) / seq_along(ordered_lfdr) <= alpha)
  if (length(eligible) == 0L) {
    return(0L)
  }
  max(eligible)
}

workflowr_root <- find_workflowr_root()
analysis_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "correlated_likelihood_sensitivity"
)
output_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fixed_diagonal_prior_likelihood_diagnostic"
)
if (file.exists(output_directory)) {
  stop("Refusing to overwrite the existing diagnostic output directory.")
}

fit_bundle <- readRDS(file.path(analysis_directory, "fit_bundle.rds"))
method_ids <- c("diagonal", "C1", "C2")
if (!identical(names(fit_bundle$raw_fits), method_ids)) {
  stop("The expected diagonal, C1, and C2 raw fits are unavailable.")
}

diagonal_fit <- fit_bundle$raw_fits$diagonal
psd_grid <- as.numeric(diagonal_fit$psd_grid)
diagonal_prior <- diagonal_fit$prior_weights
support_columns <- vapply(diagonal_prior$psd, function(psd) {
  which.min(abs(psd_grid - psd))
}, integer(1))
if (any(abs(psd_grid[support_columns] - diagonal_prior$psd) > 1e-10) ||
    anyDuplicated(support_columns)) {
  stop("The diagonal prior support cannot be mapped uniquely to the PSD grid.")
}

null_support_row <- which(diagonal_prior$psd == 0)
alternative_support_rows <- which(diagonal_prior$psd != 0)
if (length(null_support_row) != 1L ||
    length(alternative_support_rows) == 0L) {
  stop("The diagonal prior must contain one null and at least one alternative.")
}
log_prior_weights <- log(diagonal_prior$prior_weight)
conditional_alternative_weights <- diagonal_prior$prior_weight[
  alternative_support_rows
]
conditional_alternative_weights <- conditional_alternative_weights /
  sum(conditional_alternative_weights)

unit_results <- data.frame(
  pair_key = fit_bundle$pair_metadata$pair_key,
  stringsAsFactors = FALSE
)
method_summaries <- lapply(method_ids, function(method_id) {
  likelihood <- as.matrix(
    fit_bundle$raw_fits[[method_id]]$L_matrix[
      , support_columns, drop = FALSE
    ]
  )
  fixed_prior_lfdr <- vapply(seq_len(nrow(likelihood)), function(row_index) {
    log_denominator <- log_sum_exp(
      log_prior_weights + likelihood[row_index, ]
    )
    exp(
      log_prior_weights[null_support_row] +
        likelihood[row_index, null_support_row] -
        log_denominator
    )
  }, numeric(1))
  fixed_mixture_log_bf <- vapply(seq_len(nrow(likelihood)), function(row_index) {
    log_sum_exp(
      log(conditional_alternative_weights) +
        likelihood[row_index, alternative_support_rows] -
        likelihood[row_index, null_support_row]
    )
  }, numeric(1))

  unit_results[[paste0(method_id, "_fixed_prior_lfdr")]] <<-
    fixed_prior_lfdr
  unit_results[[paste0(method_id, "_fixed_mixture_log_bf")]] <<-
    fixed_mixture_log_bf

  data.frame(
    method_id = method_id,
    fixed_prior_pi0 = diagonal_prior$prior_weight[null_support_row],
    cumulative_lfdr_calls_at_005 = count_cumulative_lfdr_calls(
      fixed_prior_lfdr,
      alpha = 0.05
    ),
    units_with_lfdr_at_most_005 = sum(fixed_prior_lfdr <= 0.05),
    units_with_bf_above_one = sum(fixed_mixture_log_bf > 0),
    mean_fixed_mixture_log_bf = mean(fixed_mixture_log_bf),
    median_fixed_mixture_log_bf = stats::median(fixed_mixture_log_bf),
    boundary_score_sum_bf_minus_one = sum(
      exp(pmin(fixed_mixture_log_bf, 700)) - 1
    ),
    stringsAsFactors = FALSE
  )
})
method_summary <- do.call(rbind, method_summaries)

diagonal_reconstruction_error <- max(abs(
  unit_results$diagonal_fixed_prior_lfdr - diagonal_fit$lfdr
))
if (!is.finite(diagonal_reconstruction_error) ||
    diagonal_reconstruction_error > 1e-10) {
  stop("The fixed-prior calculation does not reconstruct diagonal lfdr.")
}

for (method_id in c("C1", "C2")) {
  unit_results[[paste0(method_id, "_minus_diagonal_log_bf")]] <-
    unit_results[[paste0(method_id, "_fixed_mixture_log_bf")]] -
    unit_results$diagonal_fixed_mixture_log_bf
}

configuration <- list(
  experiment = paste(
    "Apply the fitted diagonal raw prior weights to diagonal, C1, and C2",
    "likelihood matrices to isolate likelihood changes from prior refitting"
  ),
  n_units = nrow(unit_results),
  alpha = 0.05,
  diagonal_prior = diagonal_prior,
  diagonal_lfdr_reconstruction_maximum_error =
    diagonal_reconstruction_error,
  generated_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE)
)

dir.create(output_directory, recursive = FALSE)
saveRDS(configuration, file.path(output_directory, "configuration.rds"))
utils::write.csv(
  method_summary,
  file.path(output_directory, "method_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  unit_results,
  file.path(output_directory, "unit_results.csv"),
  row.names = FALSE
)

cat("\nFixed-diagonal-prior likelihood diagnostic completed.\n")
cat(
  "Maximum diagonal lfdr reconstruction error:",
  diagonal_reconstruction_error,
  "\n\n"
)
print(method_summary)
