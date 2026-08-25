#!/usr/bin/env Rscript

# Evaluate functional-LFSR sensitivity to nested evaluation grids for the
# currently reported early, middle, late, and switch pair sets. The fitted
# FASH model and primary category assignments are held fixed.

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

grid_rows <- function(fine_grid, step) {
  scaled <- fine_grid / step
  which(abs(scaled - round(scaled)) < 1e-8)
}

functional_lfsr <- function(samples, evaluation_grid, category, switch_threshold) {
  tolerance <- sqrt(.Machine$double.eps)
  if (category == "switch") {
    positive_max <- pmax(matrixStats::colMaxs(samples), 0)
    negative_max <- pmax(-matrixStats::colMins(samples), 0)
    statistic <- pmin(positive_max, negative_max) - switch_threshold
    return(mean(statistic <= 0))
  }

  absolute_samples <- abs(samples)
  inside <- switch(
    category,
    early = evaluation_grid <= 3 + tolerance,
    middle = evaluation_grid >= 4 - tolerance & evaluation_grid <= 11 + tolerance,
    late = evaluation_grid >= 12 - tolerance,
    stop("Unknown functional category: ", category)
  )
  if (!any(inside) || !any(!inside)) {
    stop("The evaluation grid does not support the ", category, " functional.")
  }
  statistic <- matrixStats::colMaxs(absolute_samples, rows = which(inside)) -
    matrixStats::colMaxs(absolute_samples, rows = which(!inside))
  mean(statistic <= 0)
}

load_primary_category <- function(path, object_name, category, alpha) {
  environment <- new.env(parent = emptyenv())
  load(path, envir = environment)
  if (!exists(object_name, envir = environment, inherits = FALSE)) {
    stop("Missing ", object_name, " in ", path, ".")
  }
  object <- environment[[object_name]]
  required <- c("indices", "lfsr", "cfsr")
  if (!all(required %in% names(object))) {
    stop("The category cache is incomplete: ", path)
  }
  selected <- object[object$cfsr <= alpha, , drop = FALSE]
  data.frame(
    category = category,
    index = as.integer(selected$indices),
    pair_id = rownames(selected),
    saved_primary_lfsr = as.numeric(selected$lfsr),
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
posterior_draws <- as.integer(get_arg("--posterior-draws", "3000"))
num_cores <- as.integer(get_arg("--num-cores", "8"))
seed <- as.integer(get_arg("--seed", "20260819"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
switch_threshold <- as.numeric(get_arg("--switch-threshold", "0.25"))
output_id <- get_arg("--output-id", "evaluation_grid_lfsr_sensitivity_pilot")
overwrite <- as_flag(get_arg("--overwrite", "false"))
grid_steps <- c(0.15, 0.10, 0.05)
reference_step <- 0.05
current_step <- 0.10
fine_grid <- seq(0, 15, by = reference_step)
category_order <- c("early", "middle", "late", "switch")

if (posterior_draws < 100L || num_cores < 1L || is.na(seed) ||
    !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
    !is.finite(switch_threshold) || switch_threshold <= 0 ||
    !nzchar(output_id)) {
  stop("Invalid evaluation-grid sensitivity arguments.")
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
category_specs <- data.frame(
  category = category_order,
  object_name = paste0("testing_", category_order, "_dyn"),
  file_name = paste0("classify_dyn_eQTLs_", category_order, ".RData"),
  stringsAsFactors = FALSE
)
category_specs$path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  category_specs$file_name
)
output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  output_id
)

if (!file.exists(fit_path) || any(!file.exists(category_specs$path))) {
  stop("The fitted model or at least one primary category cache is missing.")
}
if (dir.exists(output_dir) && !overwrite) {
  stop("The output directory already exists. Use --overwrite true to replace its files.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

assignments <- do.call(rbind, lapply(seq_len(nrow(category_specs)), function(i) {
  load_primary_category(
    category_specs$path[i],
    category_specs$object_name[i],
    category_specs$category[i],
    alpha
  )
}))
assignments$category <- factor(assignments$category, levels = category_order)
assignments <- assignments[order(assignments$index, assignments$category), , drop = FALSE]
rownames(assignments) <- NULL
if (nrow(assignments) == 0L ||
    anyDuplicated(assignments[, c("category", "index")]) ||
    any(!nzchar(assignments$pair_id))) {
  stop("The primary category assignments are invalid.")
}

unique_indices <- sort(unique(assignments$index))
category_counts <- as.data.frame(table(assignments$category), stringsAsFactors = FALSE)
names(category_counts) <- c("category", "n_primary_assignments")
category_counts$category <- as.character(category_counts$category)
category_counts$n_unique_pairs <- vapply(category_order, function(category) {
  length(unique(assignments$index[assignments$category == category]))
}, integer(1))

total_start <- proc.time()[["elapsed"]]
load_start <- proc.time()[["elapsed"]]
load(fit_path)
fit_load_seconds <- proc.time()[["elapsed"]] - load_start
if (!exists("fash_fit1_update") || !inherits(fash_fit1_update, "fash")) {
  stop("fash_fit1_update.RData did not contain the expected FASH fit.")
}

rows_by_step <- lapply(grid_steps, function(step) grid_rows(fine_grid, step))
names(rows_by_step) <- format(grid_steps, nsmall = 2)
expected_grid_sizes <- as.integer(round(15 / grid_steps)) + 1L
if (!all(lengths(rows_by_step) == expected_grid_sizes) ||
    any(vapply(rows_by_step, function(rows) tail(fine_grid[rows], 1L), numeric(1)) != 15)) {
  stop("The requested grid steps are not valid nested subgrids of the fine grid.")
}

sample_one_pair <- function(pair_index) {
  current_assignments <- assignments[assignments$index == pair_index, , drop = FALSE]
  categories <- as.character(current_assignments$category)
  set.seed(seed + as.integer(pair_index))
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

  rows <- vector("list", length(categories) * length(grid_steps))
  result_index <- 0L
  for (step_index in seq_along(grid_steps)) {
    grid_step <- grid_steps[step_index]
    fine_rows <- rows_by_step[[step_index]]
    grid <- fine_grid[fine_rows]
    grid_samples <- samples[fine_rows, , drop = FALSE]
    for (category in categories) {
      result_index <- result_index + 1L
      assignment <- current_assignments[current_assignments$category == category, , drop = FALSE]
      rows[[result_index]] <- data.frame(
        category = category,
        index = pair_index,
        pair_id = assignment$pair_id[1],
        saved_primary_lfsr = assignment$saved_primary_lfsr[1],
        grid_step = grid_step,
        grid_points = length(grid),
        lfsr = functional_lfsr(
          samples = grid_samples,
          evaluation_grid = grid,
          category = category,
          switch_threshold = switch_threshold
        ),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

sampling_start <- proc.time()[["elapsed"]]
if (num_cores > 1L) {
  result_list <- parallel::mclapply(
    unique_indices,
    sample_one_pair,
    mc.cores = num_cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  result_list <- lapply(unique_indices, sample_one_pair)
}
sampling_seconds <- proc.time()[["elapsed"]] - sampling_start
if (any(vapply(result_list, inherits, logical(1), "try-error"))) {
  stop("At least one parallel posterior-sampling task failed.")
}
results <- do.call(rbind, result_list)
results$category <- factor(results$category, levels = category_order)
results <- results[order(results$category, results$index, -results$grid_step), , drop = FALSE]
rownames(results) <- NULL

expected_rows <- nrow(assignments) * length(grid_steps)
if (nrow(results) != expected_rows ||
    any(!is.finite(results$lfsr)) || any(results$lfsr < 0) || any(results$lfsr > 1) ||
    any(table(interaction(results$category, results$index, drop = TRUE)) != length(grid_steps))) {
  stop("The evaluation-grid sensitivity results failed validation.")
}

reference <- results[abs(results$grid_step - reference_step) < 1e-12, , drop = FALSE]
reference_key <- paste(reference$category, reference$index, sep = "__")
reference_lfsr <- reference$lfsr
names(reference_lfsr) <- reference_key
results$key <- paste(results$category, results$index, sep = "__")
results$reference_lfsr <- unname(reference_lfsr[results$key])
results$difference_vs_fine <- results$lfsr - results$reference_lfsr
results$absolute_difference_vs_fine <- abs(results$difference_vs_fine)
results$key <- NULL
if (anyNA(results$reference_lfsr)) {
  stop("Could not align the finest-grid LFSR values.")
}

summary <- do.call(rbind, lapply(category_order, function(category) {
  do.call(rbind, lapply(grid_steps, function(step) {
    current <- results[
      results$category == category & abs(results$grid_step - step) < 1e-12,
      ,
      drop = FALSE
    ]
    data.frame(
      category = category,
      grid_step = step,
      grid_points = unique(current$grid_points),
      n_primary_pairs = nrow(current),
      mean_lfsr = mean(current$lfsr),
      median_lfsr = median(current$lfsr),
      mean_difference_vs_fine = mean(current$difference_vs_fine),
      mean_absolute_difference_vs_fine = mean(current$absolute_difference_vs_fine),
      q90_absolute_difference_vs_fine = unname(quantile(
        current$absolute_difference_vs_fine,
        0.90
      )),
      q95_absolute_difference_vs_fine = unname(quantile(
        current$absolute_difference_vs_fine,
        0.95
      )),
      maximum_absolute_difference_vs_fine = max(current$absolute_difference_vs_fine),
      spearman_vs_fine = suppressWarnings(cor(
        current$lfsr,
        current$reference_lfsr,
        method = "spearman"
      )),
      fraction_absolute_difference_gt_0p005 = mean(
        current$absolute_difference_vs_fine > 0.005
      ),
      fraction_absolute_difference_gt_0p01 = mean(
        current$absolute_difference_vs_fine > 0.01
      ),
      stringsAsFactors = FALSE
    )
  }))
}))

current <- results[abs(results$grid_step - current_step) < 1e-12, , drop = FALSE]
current_qc <- do.call(rbind, lapply(category_order, function(category) {
  category_data <- current[current$category == category, , drop = FALSE]
  data.frame(
    category = category,
    n_primary_pairs = nrow(category_data),
    saved_mean_lfsr = mean(category_data$saved_primary_lfsr),
    recomputed_mean_lfsr = mean(category_data$lfsr),
    mean_absolute_difference = mean(abs(
      category_data$lfsr - category_data$saved_primary_lfsr
    )),
    maximum_absolute_difference = max(abs(
      category_data$lfsr - category_data$saved_primary_lfsr
    )),
    spearman_correlation = suppressWarnings(cor(
      category_data$lfsr,
      category_data$saved_primary_lfsr,
      method = "spearman"
    )),
    stringsAsFactors = FALSE
  )
}))

total_seconds <- proc.time()[["elapsed"]] - total_start
runtime <- data.frame(
  stage = c("fit_load", "posterior_sampling_and_grid_evaluation", "total"),
  elapsed_seconds = c(fit_load_seconds, sampling_seconds, total_seconds)
)
configuration <- list(
  fit_path = normalizePath(fit_path),
  category_paths = normalizePath(category_specs$path),
  output_id = output_id,
  alpha = alpha,
  posterior_draws = posterior_draws,
  num_cores = num_cores,
  seed = seed,
  switch_threshold = switch_threshold,
  grid_steps = grid_steps,
  reference_step = reference_step,
  current_step = current_step,
  fine_grid = fine_grid,
  category_order = category_order,
  n_primary_assignments = nrow(assignments),
  n_unique_indices = length(unique_indices),
  estimand = paste(
    "Functional-LFSR numerical sensitivity to nested evaluation grids,",
    "conditional on the primary reported category assignments."
  )
)

write_csv(assignments, file.path(output_dir, "primary_category_assignments.csv"))
write_csv(category_counts, file.path(output_dir, "category_counts.csv"))
write_csv(results, file.path(output_dir, "pair_category_lfsr_by_grid.csv"))
write_csv(summary, file.path(output_dir, "grid_sensitivity_summary.csv"))
write_csv(current_qc, file.path(output_dir, "saved_current_grid_qc.csv"))
write_csv(runtime, file.path(output_dir, "runtime.csv"))
saveRDS(configuration, file.path(output_dir, "configuration.rds"))

cat("Evaluation-grid LFSR sensitivity completed.\n")
cat("Primary category assignments:", nrow(assignments), "\n")
cat("Unique pair indices sampled:", length(unique_indices), "\n")
cat("Grid steps:", paste(format(grid_steps, nsmall = 2), collapse = ", "), "\n")
cat("Sampling and evaluation seconds:", round(sampling_seconds, 3), "\n")
cat("Total seconds:", round(total_seconds, 3), "\n\n")
print(category_counts, row.names = FALSE)
cat("\nSensitivity summary:\n")
print(summary, row.names = FALSE)
