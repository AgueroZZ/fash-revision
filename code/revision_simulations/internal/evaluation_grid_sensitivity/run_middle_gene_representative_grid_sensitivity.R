#!/usr/bin/env Rscript

# Evaluate middle-functional LFSR sensitivity to grid refinement after fixing
# the most significant dynamic candidate pair per gene on the saved 0.10-day
# grid. The selected pair is held fixed for the paired 0.10 versus 0.05
# comparison, and both resolutions use the same posterior draws.

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

parse_gene_id <- function(pair_id) {
  gene_id <- sub("_(rs[^_]+)$", "", pair_id)
  if (any(gene_id == pair_id) || any(!nzchar(gene_id))) {
    stop("At least one pair identifier could not be parsed.")
  }
  gene_id
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
  data.frame(
    stratum = stratum,
    n_genes = nrow(data),
    mean_saved_lfsr_0p10 = mean(data$saved_lfsr_0p10),
    mean_recomputed_lfsr_0p10 = mean(data$lfsr_0p10),
    mean_lfsr_0p05 = mean(data$lfsr_0p05),
    mean_difference_0p10_minus_0p05 = mean(data$difference_0p10_minus_0p05),
    mean_absolute_difference = mean(data$absolute_difference),
    median_absolute_difference = median(data$absolute_difference),
    q90_absolute_difference = unname(quantile(data$absolute_difference, 0.90)),
    q95_absolute_difference = unname(quantile(data$absolute_difference, 0.95)),
    maximum_absolute_difference = max(data$absolute_difference),
    spearman_0p10_vs_0p05 = suppressWarnings(cor(
      data$lfsr_0p10,
      data$lfsr_0p05,
      method = "spearman"
    )),
    fraction_absolute_difference_gt_0p005 = mean(data$absolute_difference > 0.005),
    fraction_absolute_difference_gt_0p01 = mean(data$absolute_difference > 0.01),
    fraction_absolute_difference_gt_0p05 = mean(data$absolute_difference > 0.05),
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
posterior_draws <- as.integer(get_arg("--posterior-draws", "3000"))
num_cores <- as.integer(get_arg("--num-cores", "8"))
seed <- as.integer(get_arg("--seed", "20260819"))
output_id <- get_arg(
  "--output-id",
  "middle_gene_representative_grid_sensitivity_pilot"
)
overwrite <- as_flag(get_arg("--overwrite", "false"))
fine_step <- 0.05
current_step <- 0.10
fine_grid <- seq(0, 15, by = fine_step)
current_rows <- which(abs(fine_grid / current_step - round(fine_grid / current_step)) < 1e-8)
current_grid <- fine_grid[current_rows]

if (posterior_draws < 100L || num_cores < 1L || is.na(seed) || !nzchar(output_id)) {
  stop("Invalid middle grid-sensitivity arguments.")
}
if (length(current_grid) != 151L || tail(current_grid, 1L) != 15) {
  stop("The current 0.10-day grid is not nested correctly in the fine grid.")
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
middle_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "classify_dyn_eQTLs_middle.RData"
)
output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  output_id
)

if (!file.exists(fit_path) || !file.exists(middle_path)) {
  stop("The fitted model or middle-functional cache is missing.")
}
if (dir.exists(output_dir) && !overwrite) {
  stop("The output directory already exists. Use --overwrite true to replace its files.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

middle_environment <- new.env(parent = emptyenv())
load(middle_path, envir = middle_environment)
if (!exists("testing_middle_dyn", envir = middle_environment, inherits = FALSE)) {
  stop("The middle-functional cache does not contain testing_middle_dyn.")
}
testing_middle_dyn <- middle_environment$testing_middle_dyn
required_columns <- c("indices", "lfsr")
pair_ids <- rownames(testing_middle_dyn)
if (!all(required_columns %in% names(testing_middle_dyn)) ||
    is.null(pair_ids) || any(!nzchar(pair_ids)) || anyDuplicated(pair_ids) ||
    anyNA(testing_middle_dyn$lfsr) || any(!is.finite(testing_middle_dyn$lfsr))) {
  stop("The middle-functional cache is incomplete or invalid.")
}

candidate_pairs <- data.frame(
  gene_id = parse_gene_id(pair_ids),
  pair_id = pair_ids,
  index = as.integer(testing_middle_dyn$indices),
  saved_lfsr_0p10 = as.numeric(testing_middle_dyn$lfsr),
  stringsAsFactors = FALSE
)
candidate_pairs <- candidate_pairs[
  order(
    candidate_pairs$gene_id,
    candidate_pairs$saved_lfsr_0p10,
    candidate_pairs$pair_id
  ),
  ,
  drop = FALSE
]
candidate_pairs$n_candidate_pairs_for_gene <- ave(
  candidate_pairs$index,
  candidate_pairs$gene_id,
  FUN = length
)
selected_pairs <- candidate_pairs[!duplicated(candidate_pairs$gene_id), , drop = FALSE]
rownames(selected_pairs) <- NULL
if (anyDuplicated(selected_pairs$gene_id) || anyDuplicated(selected_pairs$index) ||
    nrow(selected_pairs) != length(unique(candidate_pairs$gene_id))) {
  stop("The one-pair-per-gene selection failed validation.")
}

total_start <- proc.time()[["elapsed"]]
load_start <- proc.time()[["elapsed"]]
load(fit_path)
fit_load_seconds <- proc.time()[["elapsed"]] - load_start
if (!exists("fash_fit1_update") || !inherits(fash_fit1_update, "fash")) {
  stop("fash_fit1_update.RData did not contain the expected FASH fit.")
}

sample_one_pair <- function(row_index) {
  selection <- selected_pairs[row_index, , drop = FALSE]
  pair_index <- selection$index[1]
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

  indicator_0p05 <- middle_indicator(samples, fine_grid)
  indicator_0p10 <- middle_indicator(samples[current_rows, , drop = FALSE], current_grid)
  paired_indicator_difference <- indicator_0p10 - indicator_0p05
  data.frame(
    gene_id = selection$gene_id,
    pair_id = selection$pair_id,
    index = pair_index,
    n_candidate_pairs_for_gene = selection$n_candidate_pairs_for_gene,
    saved_lfsr_0p10 = selection$saved_lfsr_0p10,
    lfsr_0p10 = mean(indicator_0p10),
    lfsr_0p05 = mean(indicator_0p05),
    difference_0p10_minus_0p05 = mean(paired_indicator_difference),
    paired_difference_mcse = sd(paired_indicator_difference) / sqrt(posterior_draws),
    discordant_posterior_draw_fraction = mean(indicator_0p10 != indicator_0p05),
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
results$absolute_difference <- abs(results$difference_0p10_minus_0p05)
results$saved_lfsr_stratum <- cut(
  results$saved_lfsr_0p10,
  breaks = c(-Inf, 0.05, 0.10, 0.25, 0.50, Inf),
  labels = c("[0, 0.05]", "(0.05, 0.10]", "(0.10, 0.25]", "(0.25, 0.50]", "(0.50, 1]"),
  right = TRUE
)
results <- results[order(results$saved_lfsr_0p10, results$gene_id), , drop = FALSE]
rownames(results) <- NULL
if (nrow(results) != nrow(selected_pairs) || anyDuplicated(results$gene_id) ||
    any(!is.finite(results$lfsr_0p10)) || any(!is.finite(results$lfsr_0p05)) ||
    any(results$lfsr_0p10 < 0) || any(results$lfsr_0p10 > 1) ||
    any(results$lfsr_0p05 < 0) || any(results$lfsr_0p05 > 1)) {
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
sensitivity_summary <- do.call(rbind, summary_rows)

saved_qc <- data.frame(
  n_genes = nrow(results),
  mean_saved_lfsr_0p10 = mean(results$saved_lfsr_0p10),
  mean_recomputed_lfsr_0p10 = mean(results$lfsr_0p10),
  mean_absolute_difference = mean(abs(results$lfsr_0p10 - results$saved_lfsr_0p10)),
  q90_absolute_difference = unname(quantile(
    abs(results$lfsr_0p10 - results$saved_lfsr_0p10),
    0.90
  )),
  maximum_absolute_difference = max(abs(results$lfsr_0p10 - results$saved_lfsr_0p10)),
  spearman_correlation = suppressWarnings(cor(
    results$lfsr_0p10,
    results$saved_lfsr_0p10,
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
  middle_path = normalizePath(middle_path),
  output_id = output_id,
  posterior_draws = posterior_draws,
  num_cores = num_cores,
  seed = seed,
  current_step = current_step,
  fine_step = fine_step,
  current_grid = current_grid,
  fine_grid = fine_grid,
  n_dynamic_candidate_pairs = nrow(candidate_pairs),
  n_genes = nrow(selected_pairs),
  estimand = paste(
    "Middle-functional LFSR numerical sensitivity to refinement from a",
    "0.10-day to a 0.05-day evaluation grid, conditional on the pair with",
    "the minimum saved 0.10-day middle LFSR within each dynamic-candidate gene."
  )
)

write_csv(selected_pairs, file.path(output_dir, "selected_gene_representative_pairs.csv"))
write_csv(results, file.path(output_dir, "middle_lfsr_grid_sensitivity_by_gene.csv"))
write_csv(sensitivity_summary, file.path(output_dir, "middle_lfsr_grid_sensitivity_summary.csv"))
write_csv(saved_qc, file.path(output_dir, "saved_current_grid_qc.csv"))
write_csv(runtime, file.path(output_dir, "runtime.csv"))
saveRDS(configuration, file.path(output_dir, "configuration.rds"))

cat("Middle one-pair-per-gene grid sensitivity completed.\n")
cat("Dynamic candidate pairs:", nrow(candidate_pairs), "\n")
cat("Genes and selected pairs:", nrow(selected_pairs), "\n")
cat("Sampling and evaluation seconds:", round(sampling_seconds, 3), "\n")
cat("Total seconds:", round(total_seconds, 3), "\n\n")
print(sensitivity_summary, row.names = FALSE)
cat("\nSaved-grid Monte Carlo QC:\n")
print(saved_qc, row.names = FALSE)
