#!/usr/bin/env Rscript

# Profile local changes in functional LFSR as the early-window cutoff moves
# around day 3. The fitted FASH model and primary early pair set are fixed.

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

profile_early_lfsr <- function(samples, evaluation_grid, cutoffs) {
  absolute_samples <- abs(samples)
  tolerance <- sqrt(.Machine$double.eps)
  inside_rows <- lapply(
    cutoffs,
    function(cutoff) which(evaluation_grid <= cutoff + tolerance)
  )
  if (any(lengths(inside_rows) < 1L) ||
      any(lengths(inside_rows) >= nrow(absolute_samples))) {
    stop("Every cutoff must leave evaluation points inside and outside the early window.")
  }

  vapply(inside_rows, function(inside) {
    outside <- setdiff(seq_len(nrow(absolute_samples)), inside)
    statistic <- matrixStats::colMaxs(absolute_samples, rows = inside) -
      matrixStats::colMaxs(absolute_samples, rows = outside)
    mean(statistic <= 0)
  }, numeric(1))
}

workflowr_root <- find_workflowr_root()
posterior_draws <- as.integer(get_arg("--posterior-draws", "3000"))
num_cores <- as.integer(get_arg("--num-cores", "8"))
seed <- as.integer(get_arg("--seed", "20260819"))
grid_step <- as.numeric(get_arg("--grid-step", "0.1"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
output_id <- get_arg("--output-id", "early_window_lfsr_local_profile_pilot")
overwrite <- as_flag(get_arg("--overwrite", "false"))
cutoffs <- seq(2.5, 3.5, by = 0.1)
primary_cutoff <- 3

if (posterior_draws < 100L || num_cores < 1L || is.na(seed) ||
    !is.finite(grid_step) || grid_step <= 0 ||
    !is.finite(alpha) || alpha <= 0 || alpha >= 1 || !nzchar(output_id)) {
  stop("Invalid early-window local-profile arguments.")
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

suppressPackageStartupMessages(library(fashr))
if (!requireNamespace("matrixStats", quietly = TRUE)) {
  stop("The matrixStats package is required for efficient profile evaluation.")
}

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
if (!any(abs(cutoffs - primary_cutoff) < 1e-12)) {
  stop("The cutoff profile must include the primary cutoff.")
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
  lfsr <- profile_early_lfsr(samples, evaluation_grid, cutoffs)
  data.frame(
    index = pair_index,
    pair_id = primary$pair_id[row_index],
    saved_primary_lfsr = primary$lfsr[row_index],
    cutoff_day = cutoffs,
    lfsr = lfsr,
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
profile <- do.call(rbind, result_list)
profile <- profile[order(profile$index, profile$cutoff_day), , drop = FALSE]
rownames(profile) <- NULL

expected_rows <- nrow(primary) * length(cutoffs)
if (nrow(profile) != expected_rows ||
    any(!is.finite(profile$lfsr)) || any(profile$lfsr < 0) || any(profile$lfsr > 1) ||
    any(table(profile$index) != length(cutoffs))) {
  stop("The early-window LFSR profile failed validation.")
}

primary_rows <- abs(profile$cutoff_day - primary_cutoff) < 1e-12
primary_lfsr <- profile$lfsr[primary_rows]
names(primary_lfsr) <- profile$pair_id[primary_rows]
profile$primary_lfsr <- unname(primary_lfsr[profile$pair_id])
profile$delta_vs_primary <- profile$lfsr - profile$primary_lfsr
if (anyNA(profile$primary_lfsr) ||
    max(abs(profile$delta_vs_primary[primary_rows])) > 1e-12) {
  stop("Could not align the primary LFSR values to the cutoff profile.")
}

profile_summary <- do.call(rbind, lapply(cutoffs, function(cutoff) {
  current <- profile[abs(profile$cutoff_day - cutoff) < 1e-12, , drop = FALSE]
  data.frame(
    cutoff_day = cutoff,
    n_primary_pairs = nrow(current),
    mean_lfsr = mean(current$lfsr),
    median_lfsr = median(current$lfsr),
    mean_delta = mean(current$delta_vs_primary),
    q10_delta = unname(quantile(current$delta_vs_primary, 0.10)),
    q25_delta = unname(quantile(current$delta_vs_primary, 0.25)),
    median_delta = median(current$delta_vs_primary),
    q75_delta = unname(quantile(current$delta_vs_primary, 0.75)),
    q90_delta = unname(quantile(current$delta_vs_primary, 0.90)),
    fraction_changed = mean(abs(current$delta_vs_primary) > 1e-12),
    stringsAsFactors = FALSE
  )
}))

finite_difference_steps <- c(0.1, 0.2, 0.5)
finite_difference <- do.call(rbind, lapply(finite_difference_steps, function(h) {
  lower <- profile_summary$mean_lfsr[
    abs(profile_summary$cutoff_day - (primary_cutoff - h)) < 1e-12
  ]
  upper <- profile_summary$mean_lfsr[
    abs(profile_summary$cutoff_day - (primary_cutoff + h)) < 1e-12
  ]
  data.frame(
    half_width_days = h,
    lower_cutoff_day = primary_cutoff - h,
    upper_cutoff_day = primary_cutoff + h,
    mean_lfsr_lower = lower,
    mean_lfsr_upper = upper,
    central_finite_difference_per_day = (upper - lower) / (2 * h)
  )
}))

saved_primary_qc <- data.frame(
  n_primary_pairs = nrow(primary),
  saved_mean_lfsr = mean(primary$lfsr),
  recomputed_mean_lfsr = mean(primary_lfsr),
  mean_absolute_difference = mean(abs(primary_lfsr - primary$lfsr)),
  median_absolute_difference = median(abs(primary_lfsr - primary$lfsr)),
  maximum_absolute_difference = max(abs(primary_lfsr - primary$lfsr)),
  spearman_correlation = cor(primary_lfsr, primary$lfsr, method = "spearman")
)

total_seconds <- proc.time()[["elapsed"]] - total_start
runtime <- data.frame(
  stage = c("fit_load", "posterior_sampling_and_profile_evaluation", "total"),
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
  primary_cutoff = primary_cutoff,
  estimand = paste(
    "Local functional-LFSR profile under alternative early-window cutoffs,",
    "conditional on the primary early pair set."
  )
)

write_csv(profile, file.path(output_dir, "early_window_pair_lfsr_profile.csv"))
write_csv(profile_summary, file.path(output_dir, "early_window_profile_summary.csv"))
write_csv(finite_difference, file.path(output_dir, "local_finite_difference.csv"))
write_csv(saved_primary_qc, file.path(output_dir, "saved_primary_lfsr_qc.csv"))
write_csv(runtime, file.path(output_dir, "runtime.csv"))
saveRDS(configuration, file.path(output_dir, "configuration.rds"))

cat("Early-window local LFSR profile completed.\n")
cat("Primary pairs:", nrow(primary), "\n")
cat("Cutoffs:", paste(format(cutoffs, nsmall = 1), collapse = ", "), "\n")
cat("Sampling and profile-evaluation seconds:", round(sampling_seconds, 3), "\n")
cat("Total seconds:", round(total_seconds, 3), "\n\n")
print(profile_summary, row.names = FALSE)
cat("\nLocal central finite differences:\n")
print(finite_difference, row.names = FALSE)
