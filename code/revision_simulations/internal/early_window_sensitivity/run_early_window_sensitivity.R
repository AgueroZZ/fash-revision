#!/usr/bin/env Rscript

# Run a conditional sensitivity analysis for pairs assigned to the early
# category under the primary time-window definition. The fitted FASH model and
# the primary dynamic discovery set are held fixed.

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

column_max <- function(x) {
  apply(x, 2L, max)
}

early_lfsr <- function(samples, evaluation_grid, cutoff) {
  inside <- evaluation_grid <= cutoff
  if (!any(inside) || !any(!inside)) {
    stop("Each cutoff must leave evaluation points inside and outside the early window.")
  }
  statistic <- column_max(abs(samples[inside, , drop = FALSE])) -
    column_max(abs(samples[!inside, , drop = FALSE]))
  mean(statistic <= 0)
}

workflowr_root <- find_workflowr_root()
posterior_draws <- as.integer(get_arg("--posterior-draws", "3000"))
num_cores <- as.integer(get_arg("--num-cores", "8"))
seed <- as.integer(get_arg("--seed", "20260819"))
grid_step <- as.numeric(get_arg("--grid-step", "0.1"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
output_id <- get_arg("--output-id", "early_window_lfsr_sensitivity_pilot")
overwrite <- as_flag(get_arg("--overwrite", "false"))
cutoffs <- c(narrower = 2.5, primary = 3, wider = 3.5)

if (posterior_draws < 100L || num_cores < 1L || is.na(seed) ||
    !is.finite(grid_step) || grid_step <= 0 ||
    !is.finite(alpha) || alpha <= 0 || alpha >= 1 || !nzchar(output_id)) {
  stop("Invalid early-window sensitivity arguments.")
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

suppressPackageStartupMessages(library(fashr))

fit_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "fash_fit1_update.RData"
)
classification_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "classify_dyn_eQTLs_early.RData"
)
output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  output_id
)

if (!file.exists(fit_path) || !file.exists(classification_path)) {
  stop("The fitted model or primary early-classification cache is missing.")
}
if (dir.exists(output_dir) && !overwrite) {
  stop("The output directory already exists. Use --overwrite true to replace its files.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

total_start <- proc.time()[["elapsed"]]
load_start <- proc.time()[["elapsed"]]
load(fit_path)
fit_load_seconds <- proc.time()[["elapsed"]] - load_start
if (!exists("fash_fit1_update") || !inherits(fash_fit1_update, "fash")) {
  stop("fash_fit1_update.RData did not contain the expected FASH fit.")
}

classification_environment <- new.env(parent = emptyenv())
load(classification_path, envir = classification_environment)
if (!exists("testing_early_dyn", envir = classification_environment, inherits = FALSE)) {
  stop("The primary early-classification cache is incomplete.")
}
testing_early_dyn <- classification_environment$testing_early_dyn
required_columns <- c("indices", "lfsr", "cfsr")
if (!all(required_columns %in% names(testing_early_dyn))) {
  stop("The primary early-classification table is incomplete.")
}

primary <- testing_early_dyn[testing_early_dyn$cfsr <= alpha, , drop = FALSE]
primary$pair_id <- rownames(primary)
primary <- primary[order(primary$indices), , drop = FALSE]
rownames(primary) <- NULL
if (nrow(primary) == 0L || anyDuplicated(primary$indices) || any(!nzchar(primary$pair_id))) {
  stop("The selected primary early pairs are invalid.")
}

evaluation_grid <- seq(0, 15, by = grid_step)
if (!isTRUE(all.equal(tail(evaluation_grid, 1L), 15))) {
  stop("The evaluation grid must include day 15 exactly.")
}

sample_one_pair <- function(row_index) {
  pair_index <- primary$indices[row_index]
  set.seed(seed + as.integer(pair_index))
  samples <- predict(
    fash_fit1_update,
    index = pair_index,
    smooth_var = evaluation_grid,
    only.samples = TRUE,
    M = posterior_draws
  )
  if (!is.matrix(samples) || nrow(samples) != length(evaluation_grid) ||
      ncol(samples) != posterior_draws || any(!is.finite(samples))) {
    stop("Posterior sampling returned an invalid matrix for index ", pair_index, ".")
  }
  lfsr <- vapply(
    cutoffs,
    function(cutoff) early_lfsr(samples, evaluation_grid, cutoff),
    numeric(1)
  )
  data.frame(
    index = pair_index,
    pair_id = primary$pair_id[row_index],
    saved_primary_lfsr = primary$lfsr[row_index],
    lfsr_cutoff_2p5 = unname(lfsr["narrower"]),
    lfsr_cutoff_3 = unname(lfsr["primary"]),
    lfsr_cutoff_3p5 = unname(lfsr["wider"]),
    stringsAsFactors = FALSE
  )
}

sampling_start <- proc.time()[["elapsed"]]
if (num_cores > 1L) {
  result_list <- parallel::mclapply(
    seq_len(nrow(primary)),
    sample_one_pair,
    mc.cores = num_cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  result_list <- lapply(seq_len(nrow(primary)), sample_one_pair)
}
sampling_seconds <- proc.time()[["elapsed"]] - sampling_start
if (any(vapply(result_list, inherits, logical(1), "try-error"))) {
  stop("At least one parallel posterior-sampling task failed.")
}
results <- do.call(rbind, result_list)
results <- results[order(results$index), , drop = FALSE]
rownames(results) <- NULL

lfsr_columns <- c("lfsr_cutoff_2p5", "lfsr_cutoff_3", "lfsr_cutoff_3p5")
if (nrow(results) != nrow(primary) || anyDuplicated(results$index) ||
    any(!is.finite(as.matrix(results[, lfsr_columns]))) ||
    any(as.matrix(results[, lfsr_columns]) < 0) ||
    any(as.matrix(results[, lfsr_columns]) > 1)) {
  stop("The early-window sensitivity results failed validation.")
}

primary_lfsr <- results$lfsr_cutoff_3
definition_summary <- do.call(rbind, lapply(seq_along(cutoffs), function(i) {
  current <- results[[lfsr_columns[i]]]
  delta <- current - primary_lfsr
  data.frame(
    definition = names(cutoffs)[i],
    early_cutoff_day = unname(cutoffs[i]),
    n_primary_pairs = nrow(results),
    mean_lfsr = mean(current),
    median_lfsr = median(current),
    q10_lfsr = unname(quantile(current, 0.10)),
    q90_lfsr = unname(quantile(current, 0.90)),
    mean_delta_vs_primary = mean(delta),
    median_delta_vs_primary = median(delta),
    q10_delta_vs_primary = unname(quantile(delta, 0.10)),
    q90_delta_vs_primary = unname(quantile(delta, 0.90)),
    spearman_vs_primary = suppressWarnings(cor(current, primary_lfsr, method = "spearman")),
    stringsAsFactors = FALSE
  )
}))

saved_primary_qc <- data.frame(
  n_primary_pairs = nrow(results),
  saved_mean_lfsr = mean(results$saved_primary_lfsr),
  recomputed_mean_lfsr = mean(results$lfsr_cutoff_3),
  mean_absolute_difference = mean(abs(results$lfsr_cutoff_3 - results$saved_primary_lfsr)),
  median_absolute_difference = median(abs(results$lfsr_cutoff_3 - results$saved_primary_lfsr)),
  maximum_absolute_difference = max(abs(results$lfsr_cutoff_3 - results$saved_primary_lfsr)),
  spearman_correlation = cor(
    results$lfsr_cutoff_3,
    results$saved_primary_lfsr,
    method = "spearman"
  )
)

total_seconds <- proc.time()[["elapsed"]] - total_start
runtime <- data.frame(
  stage = c("fit_load", "posterior_sampling_and_functional_evaluation", "total"),
  elapsed_seconds = c(fit_load_seconds, sampling_seconds, total_seconds)
)
configuration <- list(
  fit_path = normalizePath(fit_path),
  classification_path = normalizePath(classification_path),
  output_id = output_id,
  alpha = alpha,
  posterior_draws = posterior_draws,
  num_cores = num_cores,
  seed = seed,
  evaluation_grid = evaluation_grid,
  cutoffs = cutoffs,
  estimand = paste(
    "Posterior non-support probabilities under alternative early-window",
    "definitions, conditional on the primary early pair set."
  )
)

write_csv(results, file.path(output_dir, "early_window_pair_lfsr.csv"))
write_csv(definition_summary, file.path(output_dir, "early_window_summary.csv"))
write_csv(saved_primary_qc, file.path(output_dir, "saved_primary_lfsr_qc.csv"))
write_csv(runtime, file.path(output_dir, "runtime.csv"))
saveRDS(configuration, file.path(output_dir, "configuration.rds"))

cat("Early-window sensitivity completed.\n")
cat("Primary pairs:", nrow(results), "\n")
cat("Posterior draws per pair:", posterior_draws, "\n")
cat("Parallel workers:", num_cores, "\n")
cat("Sampling and functional evaluation seconds:", round(sampling_seconds, 3), "\n")
cat("Total seconds:", round(total_seconds, 3), "\n\n")
print(definition_summary, row.names = FALSE)
cat("\nSaved-versus-recomputed primary QC:\n")
print(saved_primary_qc, row.names = FALSE)
