#!/usr/bin/env Rscript

# Evaluate further middle-functional LFSR convergence from a 0.05-day grid to
# a nested 0.025-day grid. The same one-pair-per-gene representatives selected
# by saved 0.10-day middle LFSR are fixed throughout all grid comparisons.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

as_flag <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

middle_indicator <- function(samples, evaluation_grid) {
  tolerance <- sqrt(.Machine$double.eps)
  inside <- evaluation_grid >= 4 - tolerance & evaluation_grid <= 11 + tolerance
  if (!any(inside) || !any(!inside)) {
    stop("The evaluation grid does not support the middle functional.")
  }
  absolute_samples <- abs(samples)
  statistic <- matrixStats::colMaxs(
    absolute_samples,
    rows = which(inside)
  ) - matrixStats::colMaxs(
    absolute_samples,
    rows = which(!inside)
  )
  as.integer(statistic <= 0)
}

summarize_sensitivity <- function(data, stratum) {
  valid_ratio <- is.finite(data$convergence_ratio)
  data.frame(
    stratum = stratum,
    n_genes = nrow(data),
    mean_lfsr_0p05 = mean(data$lfsr_0p05),
    mean_lfsr_0p025 = mean(data$lfsr_0p025),
    mean_increase_0p05_to_0p025 = mean(data$increase_0p05_to_0p025),
    mean_absolute_difference = mean(data$absolute_difference),
    median_absolute_difference = median(data$absolute_difference),
    q90_absolute_difference = unname(quantile(data$absolute_difference, 0.90)),
    q95_absolute_difference = unname(quantile(data$absolute_difference, 0.95)),
    maximum_absolute_difference = max(data$absolute_difference),
    spearman_0p05_vs_0p025 = suppressWarnings(cor(
      data$lfsr_0p05,
      data$lfsr_0p025,
      method = "spearman"
    )),
    fraction_absolute_difference_gt_0p005 = mean(data$absolute_difference > 0.005),
    fraction_absolute_difference_gt_0p01 = mean(data$absolute_difference > 0.01),
    fraction_absolute_difference_gt_0p05 = mean(data$absolute_difference > 0.05),
    mean_previous_absolute_difference_0p10_to_0p05 = mean(
      data$previous_absolute_difference_0p10_to_0p05
    ),
    ratio_of_mean_absolute_differences = mean(data$absolute_difference) / mean(
      data$previous_absolute_difference_0p10_to_0p05
    ),
    median_pairwise_convergence_ratio = if (any(valid_ratio)) {
      median(data$convergence_ratio[valid_ratio])
    } else {
      NA_real_
    },
    fraction_new_difference_smaller_than_previous = mean(
      data$absolute_difference < data$previous_absolute_difference_0p10_to_0p05
    ),
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
posterior_draws <- as.integer(get_arg("--posterior-draws", "3000"))
num_cores <- as.integer(get_arg("--num-cores", "8"))
seed <- as.integer(get_arg("--seed", "20260819"))
output_id <- get_arg(
  "--output-id",
  "middle_gene_representative_grid_sensitivity_0p05_vs_0p025_pilot"
)
overwrite <- as_flag(get_arg("--overwrite", "false"))
fine_step <- 0.025
coarse_step <- 0.05
fine_grid <- seq(0, 15, by = fine_step)
coarse_rows <- which(abs(fine_grid / coarse_step - round(fine_grid / coarse_step)) < 1e-8)
coarse_grid <- fine_grid[coarse_rows]

if (posterior_draws < 100L || num_cores < 1L || is.na(seed) || !nzchar(output_id)) {
  stop("Invalid middle grid-sensitivity arguments.")
}
if (length(fine_grid) != 601L || length(coarse_grid) != 301L ||
    tail(fine_grid, 1L) != 15 || tail(coarse_grid, 1L) != 15) {
  stop("The 0.05-day grid is not nested correctly in the 0.025-day grid.")
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

suppressPackageStartupMessages(library(fashr))
if (!requireNamespace("matrixStats", quietly = TRUE)) {
  stop("The matrixStats package is required for efficient functional evaluation.")
}

fit_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "fash_fit1_update.RData"
)
previous_output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "middle_gene_representative_grid_sensitivity_pilot"
)
selection_path <- file.path(previous_output_dir, "selected_gene_representative_pairs.csv")
previous_results_path <- file.path(
  previous_output_dir,
  "middle_lfsr_grid_sensitivity_by_gene.csv"
)
output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  output_id
)

if (!file.exists(fit_path) || !file.exists(selection_path) ||
    !file.exists(previous_results_path)) {
  stop("The fitted model or previous fixed-pair grid-sensitivity output is missing.")
}
if (dir.exists(output_dir) && !overwrite) {
  stop("The output directory already exists. Use --overwrite true to replace its files.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

selected_pairs <- read.csv(selection_path, stringsAsFactors = FALSE)
previous_results <- read.csv(previous_results_path, stringsAsFactors = FALSE)
selection_columns <- c(
  "gene_id",
  "pair_id",
  "index",
  "n_candidate_pairs_for_gene",
  "saved_lfsr_0p10"
)
previous_columns <- c(
  "gene_id",
  "pair_id",
  "index",
  "lfsr_0p10",
  "lfsr_0p05",
  "absolute_difference"
)
if (!all(selection_columns %in% names(selected_pairs)) ||
    !all(previous_columns %in% names(previous_results)) ||
    nrow(selected_pairs) == 0L || anyDuplicated(selected_pairs$gene_id) ||
    anyDuplicated(selected_pairs$pair_id) || anyDuplicated(selected_pairs$index) ||
    anyDuplicated(previous_results$gene_id)) {
  stop("The previous fixed-pair outputs are incomplete or invalid.")
}
previous_rows <- match(selected_pairs$gene_id, previous_results$gene_id)
if (anyNA(previous_rows) ||
    !identical(selected_pairs$pair_id, previous_results$pair_id[previous_rows]) ||
    !identical(as.integer(selected_pairs$index), as.integer(previous_results$index[previous_rows]))) {
  stop("The previous selection and result files do not contain identical fixed pairs.")
}
previous_results <- previous_results[previous_rows, , drop = FALSE]

total_start <- proc.time()[["elapsed"]]
load_start <- proc.time()[["elapsed"]]
load(fit_path)
fit_load_seconds <- proc.time()[["elapsed"]] - load_start
if (!exists("fash_fit1_update") || !inherits(fash_fit1_update, "fash")) {
  stop("fash_fit1_update.RData did not contain the expected FASH fit.")
}

sample_one_pair <- function(row_index) {
  selection <- selected_pairs[row_index, , drop = FALSE]
  pair_index <- as.integer(selection$index[1])
  set.seed(seed + pair_index)
  samples <- predict(
    fash_fit1_update,
    index = pair_index,
    smooth_var = fine_grid,
    only.samples = TRUE,
    M = posterior_draws
  )
  if (!is.matrix(samples) || nrow(samples) != length(fine_grid) ||
      ncol(samples) != posterior_draws || any(!is.finite(samples))) {
    stop("Posterior sampling returned an invalid matrix for index ", pair_index, ".")
  }

  indicator_0p025 <- middle_indicator(samples, fine_grid)
  indicator_0p05 <- middle_indicator(samples[coarse_rows, , drop = FALSE], coarse_grid)
  paired_indicator_difference <- indicator_0p05 - indicator_0p025
  data.frame(
    gene_id = selection$gene_id,
    pair_id = selection$pair_id,
    index = pair_index,
    n_candidate_pairs_for_gene = selection$n_candidate_pairs_for_gene,
    saved_lfsr_0p10 = selection$saved_lfsr_0p10,
    lfsr_0p05 = mean(indicator_0p05),
    lfsr_0p025 = mean(indicator_0p025),
    difference_0p05_minus_0p025 = mean(paired_indicator_difference),
    paired_difference_mcse = sd(paired_indicator_difference) / sqrt(posterior_draws),
    discordant_posterior_draw_fraction = mean(indicator_0p05 != indicator_0p025),
    stringsAsFactors = FALSE
  )
}

sampling_start <- proc.time()[["elapsed"]]
row_indices <- seq_len(nrow(selected_pairs))
if (num_cores > 1L) {
  result_list <- parallel::mclapply(
    row_indices,
    sample_one_pair,
    mc.cores = num_cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  result_list <- lapply(row_indices, sample_one_pair)
}
sampling_seconds <- proc.time()[["elapsed"]] - sampling_start
if (any(vapply(result_list, inherits, logical(1), "try-error"))) {
  stop("At least one parallel posterior-sampling task failed.")
}

results <- do.call(rbind, result_list)
previous_rows <- match(results$gene_id, previous_results$gene_id)
results$previous_lfsr_0p10 <- previous_results$lfsr_0p10[previous_rows]
results$previous_lfsr_0p05 <- previous_results$lfsr_0p05[previous_rows]
results$previous_absolute_difference_0p10_to_0p05 <- previous_results$absolute_difference[
  previous_rows
]
results$increase_0p05_to_0p025 <- results$lfsr_0p025 - results$lfsr_0p05
results$absolute_difference <- abs(results$difference_0p05_minus_0p025)
results$convergence_ratio <- results$absolute_difference /
  results$previous_absolute_difference_0p10_to_0p05
results$convergence_ratio[
  results$previous_absolute_difference_0p10_to_0p05 == 0
] <- NA_real_
results$saved_lfsr_stratum <- cut(
  results$saved_lfsr_0p10,
  breaks = c(-Inf, 0.05, 0.10, 0.25, 0.50, Inf),
  labels = c("[0, 0.05]", "(0.05, 0.10]", "(0.10, 0.25]", "(0.25, 0.50]", "(0.50, 1]"),
  right = TRUE
)
results <- results[order(results$saved_lfsr_0p10, results$gene_id), , drop = FALSE]
rownames(results) <- NULL
if (nrow(results) != nrow(selected_pairs) || anyDuplicated(results$gene_id) ||
    any(!is.finite(results$lfsr_0p05)) || any(!is.finite(results$lfsr_0p025)) ||
    any(results$lfsr_0p05 < 0) || any(results$lfsr_0p05 > 1) ||
    any(results$lfsr_0p025 < 0) || any(results$lfsr_0p025 > 1)) {
  stop("The paired grid-sensitivity results failed validation.")
}

summary_rows <- list(summarize_sensitivity(results, "All genes"))
for (stratum in levels(results$saved_lfsr_stratum)) {
  stratum_data <- results[results$saved_lfsr_stratum == stratum, , drop = FALSE]
  if (nrow(stratum_data) > 0L) {
    summary_rows[[length(summary_rows) + 1L]] <- summarize_sensitivity(
      stratum_data,
      paste0("Saved LFSR ", stratum)
    )
  }
}
for (cutoff in c(0.10, 0.25, 0.50)) {
  cutoff_data <- results[results$saved_lfsr_0p10 <= cutoff, , drop = FALSE]
  summary_rows[[length(summary_rows) + 1L]] <- summarize_sensitivity(
    cutoff_data,
    paste0("Saved LFSR <= ", format(cutoff, nsmall = 2))
  )
}
sensitivity_summary <- do.call(rbind, summary_rows)

coarse_grid_qc <- data.frame(
  n_genes = nrow(results),
  mean_previous_lfsr_0p05 = mean(results$previous_lfsr_0p05),
  mean_recomputed_lfsr_0p05 = mean(results$lfsr_0p05),
  mean_absolute_difference = mean(abs(
    results$lfsr_0p05 - results$previous_lfsr_0p05
  )),
  q90_absolute_difference = unname(quantile(
    abs(results$lfsr_0p05 - results$previous_lfsr_0p05),
    0.90
  )),
  maximum_absolute_difference = max(abs(
    results$lfsr_0p05 - results$previous_lfsr_0p05
  )),
  spearman_correlation = suppressWarnings(cor(
    results$lfsr_0p05,
    results$previous_lfsr_0p05,
    method = "spearman"
  )),
  stringsAsFactors = FALSE
)

total_seconds <- proc.time()[["elapsed"]] - total_start
runtime <- data.frame(
  stage = c("fit_load", "posterior_sampling_and_grid_evaluation", "total"),
  elapsed_seconds = c(fit_load_seconds, sampling_seconds, total_seconds)
)
configuration <- list(
  fit_path = normalizePath(fit_path),
  selection_path = normalizePath(selection_path),
  previous_results_path = normalizePath(previous_results_path),
  output_id = output_id,
  posterior_draws = posterior_draws,
  num_cores = num_cores,
  seed = seed,
  coarse_step = coarse_step,
  fine_step = fine_step,
  coarse_grid = coarse_grid,
  fine_grid = fine_grid,
  n_genes = nrow(selected_pairs),
  estimand = paste(
    "Further middle-functional LFSR numerical convergence from a 0.05-day",
    "to a nested 0.025-day evaluation grid, conditional on the same pair",
    "with minimum saved 0.10-day middle LFSR within each gene."
  )
)

write_csv(selected_pairs, file.path(output_dir, "selected_gene_representative_pairs.csv"))
write_csv(results, file.path(output_dir, "middle_lfsr_grid_sensitivity_by_gene.csv"))
write_csv(sensitivity_summary, file.path(output_dir, "middle_lfsr_grid_sensitivity_summary.csv"))
write_csv(coarse_grid_qc, file.path(output_dir, "previous_0p05_grid_qc.csv"))
write_csv(runtime, file.path(output_dir, "runtime.csv"))
saveRDS(configuration, file.path(output_dir, "configuration.rds"))

cat("Middle 0.05 versus 0.025 grid sensitivity completed.\n")
cat("Fixed gene-representative pairs:", nrow(selected_pairs), "\n")
cat("Sampling and evaluation seconds:", round(sampling_seconds, 3), "\n")
cat("Total seconds:", round(total_seconds, 3), "\n\n")
print(sensitivity_summary, row.names = FALSE)
cat("\nPrevious 0.05-grid Monte Carlo QC:\n")
print(coarse_grid_qc, row.names = FALSE)
