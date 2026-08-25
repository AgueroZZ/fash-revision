#!/usr/bin/env Rscript

# Compare top-variant and random-variant residual-block permutation pilots.

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
  top_discovered = file.path(
    output_parent,
    "selected_signal_residual_block_permutation_seed20260811"
  ),
  random_all_tested = file.path(
    output_parent,
    paste0(
      "selected_signal_random_all_tested_residual_block_permutation_",
      "selection20260817_seed20260811"
    )
  )
)
required_filenames <- c(
  "calibration_diagnostics.csv", "lfdr_quantiles.csv", "configuration.rds",
  "selection.csv", "unit_lfdr.csv", "merged_fash_fit.rds",
  "donor_permutation.csv"
)
required_files <- unlist(lapply(pilot_directories, function(directory) {
  file.path(directory, required_filenames)
}), use.names = FALSE)
if (any(!file.exists(required_files))) {
  stop("At least one completed residual-block pilot cache is missing.")
}

read_pilot <- function(directory) {
  list(
    configuration = readRDS(file.path(directory, "configuration.rds")),
    selection = utils::read.csv(
      file.path(directory, "selection.csv"),
      stringsAsFactors = FALSE
    ),
    calibration = utils::read.csv(
      file.path(directory, "calibration_diagnostics.csv"),
      stringsAsFactors = FALSE
    ),
    lfdr_quantiles = utils::read.csv(
      file.path(directory, "lfdr_quantiles.csv"),
      stringsAsFactors = FALSE
    ),
    unit_lfdr = utils::read.csv(
      file.path(directory, "unit_lfdr.csv"),
      stringsAsFactors = FALSE
    ),
    fit_bundle = readRDS(file.path(directory, "merged_fash_fit.rds")),
    donor_map_md5 = unname(tools::md5sum(file.path(
      directory,
      "donor_permutation.csv"
    )))
  )
}
pilots <- lapply(pilot_directories, read_pilot)
if (pilots$top_discovered$donor_map_md5 !=
    pilots$random_all_tested$donor_map_md5) {
  stop("The top and random pilots did not use the same donor map.")
}
if (!setequal(
  pilots$top_discovered$selection$gene_id,
  pilots$random_all_tested$selection$gene_id
)) {
  stop("The top and random pilots do not contain the same target genes.")
}

calibration <- do.call(rbind, lapply(names(pilots), function(selection_method) {
  result <- pilots[[selection_method]]$calibration
  result$target_selection_method <- selection_method
  result
}))
rownames(calibration) <- NULL
calibration <- calibration[, c(
  "target_selection_method",
  "fit_stage",
  "alpha",
  "n_target",
  "n_permuted_null",
  "target_calls",
  "permuted_null_calls",
  "total_calls",
  "target_call_rate",
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

difference_columns <- c(
  "target_calls", "permuted_null_calls", "total_calls", "target_call_rate",
  "permuted_null_call_rate", "known_null_discovery_fraction", "pi0_merged",
  "scaled_fdr_merged_from_estimated_pi0", "pi0_target_unbounded"
)
calibration_difference <- do.call(rbind, lapply(
  unique(calibration$fit_stage),
  function(stage) {
    stage_rows <- calibration[calibration$fit_stage == stage, , drop = FALSE]
    top_row <- stage_rows[
      stage_rows$target_selection_method == "top_discovered",
      ,
      drop = FALSE
    ]
    random_row <- stage_rows[
      stage_rows$target_selection_method == "random_all_tested",
      ,
      drop = FALSE
    ]
    data.frame(
      fit_stage = stage,
      metric = difference_columns,
      top_discovered = as.numeric(top_row[1L, difference_columns]),
      random_all_tested = as.numeric(random_row[1L, difference_columns]),
      random_minus_top = as.numeric(random_row[1L, difference_columns]) -
        as.numeric(top_row[1L, difference_columns]),
      stringsAsFactors = FALSE
    )
  }
))
rownames(calibration_difference) <- NULL

random_candidate_counts <- pilots$random_all_tested$selection[, c(
  "gene_id", "candidate_variant_count"
)]
selection_summary <- do.call(rbind, lapply(
  names(pilots),
  function(selection_method) {
    selection <- pilots[[selection_method]]$selection
    candidate_rows <- match(selection$gene_id, random_candidate_counts$gene_id)
    candidate_counts <- random_candidate_counts$candidate_variant_count[
      candidate_rows
    ]
    data.frame(
      target_selection_method = selection_method,
      n_target_genes = nrow(selection),
      n_original_pair_level_discoveries =
        sum(selection$source_pair_level_discovery),
      proportion_original_pair_level_discoveries =
        mean(selection$source_pair_level_discovery),
      minimum_candidate_variants = min(candidate_counts),
      median_candidate_variants = stats::median(candidate_counts),
      mean_candidate_variants = mean(candidate_counts),
      maximum_candidate_variants = max(candidate_counts),
      minimum_source_bf_lfdr = min(selection$source_bf_lfdr),
      q25_source_bf_lfdr = unname(stats::quantile(
        selection$source_bf_lfdr,
        0.25
      )),
      median_source_bf_lfdr = stats::median(selection$source_bf_lfdr),
      mean_source_bf_lfdr = mean(selection$source_bf_lfdr),
      q75_source_bf_lfdr = unname(stats::quantile(
        selection$source_bf_lfdr,
        0.75
      )),
      maximum_source_bf_lfdr = max(selection$source_bf_lfdr),
      stringsAsFactors = FALSE
    )
  }
))
rownames(selection_summary) <- NULL

lfdr_quantiles <- do.call(rbind, lapply(
  names(pilots),
  function(selection_method) {
    result <- pilots[[selection_method]]$lfdr_quantiles
    result$target_selection_method <- selection_method
    result
  }
))
rownames(lfdr_quantiles) <- NULL
lfdr_quantiles <- lfdr_quantiles[, c(
  "target_selection_method",
  setdiff(names(lfdr_quantiles), "target_selection_method")
)]

runtime <- do.call(rbind, lapply(names(pilots), function(selection_method) {
  configuration <- pilots[[selection_method]]$configuration
  data.frame(
    target_selection_method = selection_method,
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

bf_summary <- do.call(rbind, lapply(names(pilots), function(selection_method) {
  fit_bundle <- pilots[[selection_method]]$fit_bundle
  n_target <- fit_bundle$configuration$n_target_units
  bayes_factor <- as.numeric(fit_bundle$bf_adjusted_fit$BF)
  if (length(bayes_factor) != 2L * n_target ||
      any(!is.finite(bayes_factor)) || any(bayes_factor <= 0)) {
    stop("Invalid Bayes factors for selection method: ", selection_method)
  }
  group <- rep(c("target", "permuted_null"), each = n_target)
  do.call(rbind, lapply(unique(group), function(group_name) {
    values <- bayes_factor[group == group_name]
    descending <- order(values, decreasing = TRUE)
    n_top_1 <- ceiling(0.01 * length(values))
    n_top_5 <- ceiling(0.05 * length(values))
    data.frame(
      target_selection_method = selection_method,
      group = group_name,
      n_units = length(values),
      minimum_bf = min(values),
      q05_bf = unname(stats::quantile(values, 0.05)),
      q25_bf = unname(stats::quantile(values, 0.25)),
      median_bf = stats::median(values),
      mean_bf = mean(values),
      q75_bf = unname(stats::quantile(values, 0.75)),
      q95_bf = unname(stats::quantile(values, 0.95)),
      q99_bf = unname(stats::quantile(values, 0.99)),
      maximum_bf = max(values),
      proportion_bf_greater_than_one = mean(values > 1),
      top_1_percent_bf_mass_share =
        sum(values[descending[seq_len(n_top_1)]]) / sum(values),
      top_5_percent_bf_mass_share =
        sum(values[descending[seq_len(n_top_5)]]) / sum(values),
      stringsAsFactors = FALSE
    )
  }))
}))
rownames(bf_summary) <- NULL

selection_behavior <- do.call(rbind, lapply(
  names(pilots),
  function(selection_method) {
    unit_lfdr <- pilots[[selection_method]]$unit_lfdr
    do.call(rbind, lapply(unique(unit_lfdr$fit_stage), function(stage) {
      stage_data <- unit_lfdr[unit_lfdr$fit_stage == stage, , drop = FALSE]
      calls <- seq_len(nrow(stage_data)) %in% cumulative_lfdr_calls(
        stage_data$lfdr,
        alpha = 0.05
      )
      null <- stage_data$group == "permuted_null"
      target <- stage_data$group == "target"
      data.frame(
        target_selection_method = selection_method,
        fit_stage = stage,
        individual_target_lfdr_le_0_05 =
          sum(target & stage_data$lfdr <= 0.05),
        individual_null_lfdr_le_0_05 =
          sum(null & stage_data$lfdr <= 0.05),
        cumulative_rule_target_calls = sum(target & calls),
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

random_directory <- pilot_directories[["random_all_tested"]]
write_comparison <- function(x, filename) {
  utils::write.csv(
    x,
    file.path(random_directory, filename),
    row.names = FALSE
  )
}
write_comparison(calibration, "comparison_top_random_calibration.csv")
write_comparison(
  calibration_difference,
  "comparison_top_random_calibration_differences.csv"
)
write_comparison(selection_summary, "comparison_top_random_selection.csv")
write_comparison(lfdr_quantiles, "comparison_top_random_lfdr_quantiles.csv")
write_comparison(runtime, "comparison_top_random_runtime.csv")
write_comparison(bf_summary, "comparison_top_random_bf_summary.csv")
write_comparison(
  selection_behavior,
  "comparison_top_random_selection_behavior.csv"
)

cat("Top-versus-random residual-block pilot comparison completed.\n")
print(calibration, row.names = FALSE)
