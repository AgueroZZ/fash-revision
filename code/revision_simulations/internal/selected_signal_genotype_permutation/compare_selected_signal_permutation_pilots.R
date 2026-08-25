#!/usr/bin/env Rscript

# Compare the fixed-seed genotype-label and residual-block permutation pilots.

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

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation",
  "selected_signal_genotype_permutation_helpers.R"
))
output_parent <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal"
)
pilot_directories <- c(
  genotype_label = file.path(
    output_parent,
    "selected_signal_genotype_permutation_seed20260811"
  ),
  residual_block = file.path(
    output_parent,
    "selected_signal_residual_block_permutation_seed20260811"
  )
)
required_files <- unlist(lapply(pilot_directories, function(directory) {
  file.path(
    directory,
    c("calibration_diagnostics.csv", "lfdr_quantiles.csv", "configuration.rds")
  )
}), use.names = FALSE)
if (any(!file.exists(required_files))) {
  stop("At least one completed permutation-pilot cache is missing.")
}

calibration <- do.call(rbind, lapply(names(pilot_directories), function(method) {
  result <- utils::read.csv(
    file.path(pilot_directories[[method]], "calibration_diagnostics.csv"),
    stringsAsFactors = FALSE
  )
  result$permutation_method <- method
  result
}))
rownames(calibration) <- NULL
calibration <- calibration[, c(
  "permutation_method",
  "fit_stage",
  "alpha",
  "n_target",
  "n_permuted_null",
  "target_calls",
  "permuted_null_calls",
  "total_calls",
  "permuted_null_call_rate",
  "known_null_discovery_fraction",
  "pi0_merged",
  "pi0_merged_below_design_lower_bound",
  "scaled_fdr_merged_from_estimated_pi0",
  "scaled_fdr_merged_from_design_lower_bound",
  "pi0_target_unbounded",
  "pi0_target_valid",
  "post_selection_fdr_target_from_pi0"
)]

lfdr_quantiles <- do.call(rbind, lapply(
  names(pilot_directories),
  function(method) {
    result <- utils::read.csv(
      file.path(pilot_directories[[method]], "lfdr_quantiles.csv"),
      stringsAsFactors = FALSE
    )
    result$permutation_method <- method
    result
  }
))
rownames(lfdr_quantiles) <- NULL
lfdr_quantiles <- lfdr_quantiles[, c(
  "permutation_method",
  setdiff(names(lfdr_quantiles), "permutation_method")
)]

runtime <- do.call(rbind, lapply(names(pilot_directories), function(method) {
  configuration <- readRDS(
    file.path(pilot_directories[[method]], "configuration.rds")
  )
  data.frame(
    permutation_method = method,
    null_likelihood_elapsed_seconds =
      configuration$null_likelihood_elapsed_seconds,
    raw_refit_elapsed_seconds = configuration$raw_refit_elapsed_seconds,
    bf_update_elapsed_seconds = configuration$bf_update_elapsed_seconds,
    total_elapsed_seconds = configuration$total_elapsed_seconds,
    null_likelihood_warning_count =
      length(configuration$null_likelihood_warnings),
    raw_refit_warning_count = length(configuration$raw_refit_warnings),
    bf_update_warning_count = length(configuration$bf_update_warnings),
    stringsAsFactors = FALSE
  )
}))
rownames(runtime) <- NULL

bf_summary <- do.call(rbind, lapply(names(pilot_directories), function(method) {
  fit_bundle <- readRDS(
    file.path(pilot_directories[[method]], "merged_fash_fit.rds")
  )
  n_target <- fit_bundle$configuration$n_target_units
  null_bf <- as.numeric(
    fit_bundle$bf_adjusted_fit$BF[n_target + seq_len(n_target)]
  )
  descending <- order(null_bf, decreasing = TRUE)
  n_top_1 <- ceiling(0.01 * length(null_bf))
  n_top_5 <- ceiling(0.05 * length(null_bf))
  data.frame(
    permutation_method = method,
    n_permuted_null = length(null_bf),
    minimum_bf = min(null_bf),
    q05_bf = unname(stats::quantile(null_bf, 0.05)),
    q25_bf = unname(stats::quantile(null_bf, 0.25)),
    median_bf = stats::median(null_bf),
    mean_bf = mean(null_bf),
    q75_bf = unname(stats::quantile(null_bf, 0.75)),
    q95_bf = unname(stats::quantile(null_bf, 0.95)),
    q99_bf = unname(stats::quantile(null_bf, 0.99)),
    maximum_bf = max(null_bf),
    proportion_bf_greater_than_one = mean(null_bf > 1),
    top_1_percent_bf_mass_share =
      sum(null_bf[descending[seq_len(n_top_1)]]) / sum(null_bf),
    top_5_percent_bf_mass_share =
      sum(null_bf[descending[seq_len(n_top_5)]]) / sum(null_bf),
    stringsAsFactors = FALSE
  )
}))
rownames(bf_summary) <- NULL

selection_behavior <- do.call(rbind, lapply(
  names(pilot_directories),
  function(method) {
    unit_lfdr <- utils::read.csv(
      file.path(pilot_directories[[method]], "unit_lfdr.csv"),
      stringsAsFactors = FALSE
    )
    do.call(rbind, lapply(unique(unit_lfdr$fit_stage), function(stage) {
      stage_data <- unit_lfdr[unit_lfdr$fit_stage == stage, , drop = FALSE]
      calls <- seq_len(nrow(stage_data)) %in% cumulative_lfdr_calls(
        stage_data$lfdr,
        alpha = 0.05
      )
      null <- stage_data$group == "permuted_null"
      data.frame(
        permutation_method = method,
        fit_stage = stage,
        individual_null_lfdr_le_0_05 = sum(null & stage_data$lfdr <= 0.05),
        cumulative_rule_null_calls = sum(null & calls),
        maximum_selected_lfdr = if (any(calls)) {
          max(stage_data$lfdr[calls])
        } else {
          NA_real_
        },
        mean_selected_lfdr = if (any(calls)) {
          mean(stage_data$lfdr[calls])
        } else {
          NA_real_
        },
        stringsAsFactors = FALSE
      )
    }))
  }
))
rownames(selection_behavior) <- NULL

residual_directory <- pilot_directories[["residual_block"]]
utils::write.csv(
  calibration,
  file.path(residual_directory, "comparison_calibration.csv"),
  row.names = FALSE
)
utils::write.csv(
  lfdr_quantiles,
  file.path(residual_directory, "comparison_lfdr_quantiles.csv"),
  row.names = FALSE
)
utils::write.csv(
  runtime,
  file.path(residual_directory, "comparison_runtime.csv"),
  row.names = FALSE
)
utils::write.csv(
  bf_summary,
  file.path(residual_directory, "comparison_null_bf_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  selection_behavior,
  file.path(residual_directory, "comparison_selection_behavior.csv"),
  row.names = FALSE
)

cat("Selected-signal permutation-pilot comparison completed.\n")
print(calibration, row.names = FALSE)
