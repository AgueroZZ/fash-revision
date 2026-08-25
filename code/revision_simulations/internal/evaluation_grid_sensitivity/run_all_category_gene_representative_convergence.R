#!/usr/bin/env Rscript

# Evaluate functional-LFSR convergence for one representative pair per gene
# and category. The representative is the pair with the smallest saved primary
# LFSR within each gene. Separate fixed-seed posterior calls are used for each
# Monte Carlo sample size, and every call is evaluated on nested 0.10, 0.05,
# and 0.025 time grids.

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

grid_rows <- function(fine_grid, step) {
  scaled <- fine_grid / step
  which(abs(scaled - round(scaled)) < 1e-8)
}

functional_indicator <- function(samples,
                                 evaluation_grid,
                                 category,
                                 switch_threshold) {
  tolerance <- sqrt(.Machine$double.eps)

  if (category == "switch") {
    positive_maximum <- pmax(matrixStats::colMaxs(samples), 0)
    negative_maximum <- pmax(-matrixStats::colMins(samples), 0)
    switch_statistic <- pmin(positive_maximum, negative_maximum) - switch_threshold
    return(as.integer(switch_statistic <= 0))
  }

  inside_window <- switch(
    category,
    early = evaluation_grid <= 3 + tolerance,
    middle = evaluation_grid >= 4 - tolerance & evaluation_grid <= 11 + tolerance,
    late = evaluation_grid >= 12 - tolerance,
    stop("Unknown functional category: ", category)
  )
  if (!any(inside_window) || !any(!inside_window)) {
    stop("The evaluation grid does not support the ", category, " functional.")
  }

  absolute_samples <- abs(samples)
  location_statistic <-
    matrixStats::colMaxs(absolute_samples, rows = which(inside_window)) -
    matrixStats::colMaxs(absolute_samples, rows = which(!inside_window))
  as.integer(location_statistic <= 0)
}

load_primary_assignments <- function(path, object_name, category, alpha) {
  input_environment <- new.env(parent = emptyenv())
  load(path, envir = input_environment)
  if (!exists(object_name, envir = input_environment, inherits = FALSE)) {
    stop("Missing ", object_name, " in ", path, ".")
  }

  category_table <- input_environment[[object_name]]
  required_columns <- c("indices", "lfsr", "cfsr")
  if (!all(required_columns %in% names(category_table))) {
    stop("The primary category cache is incomplete: ", path)
  }

  selected <- category_table[category_table$cfsr <= alpha, , drop = FALSE]
  data.frame(
    category = category,
    index = as.integer(selected$indices),
    pair_id = rownames(selected),
    saved_primary_lfsr = as.numeric(selected$lfsr),
    stringsAsFactors = FALSE
  )
}

summarize_pairwise_comparison <- function(current,
                                          category,
                                          comparison,
                                          estimate_column,
                                          reference_column) {
  estimate <- current[[estimate_column]]
  reference <- current[[reference_column]]
  absolute_difference <- abs(estimate - reference)
  estimate_selected <- estimate <= 0.05
  reference_selected <- reference <= 0.05
  selected_union <- sum(estimate_selected | reference_selected)

  data.frame(
    category = category,
    comparison = comparison,
    n_pairs = nrow(current),
    mean_signed_change = mean(estimate - reference),
    mean_absolute_change = mean(absolute_difference),
    median_absolute_change = median(absolute_difference),
    q90_absolute_change = unname(quantile(absolute_difference, 0.90)),
    maximum_absolute_change = max(absolute_difference),
    spearman_correlation = suppressWarnings(cor(
      estimate,
      reference,
      method = "spearman"
    )),
    fraction_within_0p005 = mean(absolute_difference <= 0.005),
    fraction_within_0p01 = mean(absolute_difference <= 0.01),
    n_lfsr_below_0p05_estimate = sum(estimate_selected),
    n_lfsr_below_0p05_reference = sum(reference_selected),
    lfsr_0p05_classification_agreement = mean(estimate_selected == reference_selected),
    lfsr_0p05_jaccard = if (selected_union == 0L) {
      1
    } else {
      sum(estimate_selected & reference_selected) / selected_union
    },
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
num_cores <- as.integer(get_arg("--num-cores", "2"))
seed <- as.integer(get_arg("--seed", "20260820"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
switch_threshold <- as.numeric(get_arg("--switch-threshold", "0.25"))
output_id <- get_arg(
  "--output-id",
  "evaluation_grid_mc_all_category_gene_representatives"
)

category_order <- c("early", "middle", "late", "switch")
expected_representative_counts <- c(
  early = 8L,
  middle = 5L,
  late = 12L,
  switch = 250L
)
grid_steps <- c(0.10, 0.05, 0.025)
grid_labels <- c("0.10", "0.05", "0.025")
posterior_draw_values <- c(3000L, 10000L, 30000L)
fine_grid <- seq(0, 15, by = 0.025)

if (is.na(num_cores) || num_cores < 1L || num_cores > 2L ||
    is.na(seed) || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
    !is.finite(switch_threshold) || switch_threshold <= 0 ||
    !nzchar(output_id)) {
  stop("Invalid convergence-analysis arguments.")
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

suppressPackageStartupMessages(library(fashr))
if (!requireNamespace("matrixStats", quietly = TRUE)) {
  stop("The matrixStats package is required.")
}

fit_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "fash_fit1_update.RData"
)
gene_map_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "cache_gene_map.rds"
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
staging_dir <- paste0(output_dir, ".staging-", Sys.getpid())

if (!file.exists(fit_path) || !file.exists(gene_map_path) ||
    any(!file.exists(category_specs$path))) {
  stop("The fitted model, gene map, or at least one category cache is missing.")
}
if (dir.exists(output_dir) || dir.exists(staging_dir)) {
  stop("The requested output or staging directory already exists: ", output_dir)
}
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

all_assignments <- do.call(rbind, lapply(seq_len(nrow(category_specs)), function(i) {
  load_primary_assignments(
    category_specs$path[i],
    category_specs$object_name[i],
    category_specs$category[i],
    alpha
  )
}))
all_assignments$gene_id <- sub("_(rs[^_]+)$", "", all_assignments$pair_id)
all_assignments$variant_id <- sub("^.*_(rs[^_]+)$", "\\1", all_assignments$pair_id)
all_assignments <- all_assignments[
  order(
    match(all_assignments$category, category_order),
    all_assignments$gene_id,
    all_assignments$saved_primary_lfsr,
    all_assignments$pair_id
  ),
  ,
  drop = FALSE
]

# The first row within each category-gene block has the smallest saved LFSR,
# with pair_id providing a deterministic tie-breaker.
representative_key <- paste(all_assignments$category, all_assignments$gene_id)
representatives <- all_assignments[!duplicated(representative_key), , drop = FALSE]
rownames(representatives) <- NULL

gene_map <- readRDS(gene_map_path)
representatives$gene_symbol <- gene_map$hgnc_symbol[
  match(representatives$gene_id, gene_map$ensembl_gene_id)
]
missing_symbol <- is.na(representatives$gene_symbol) |
  !nzchar(representatives$gene_symbol)
representatives$gene_symbol[missing_symbol] <- representatives$gene_id[missing_symbol]

observed_counts <- table(factor(representatives$category, levels = category_order))
if (!identical(as.integer(observed_counts), unname(expected_representative_counts)) ||
    anyDuplicated(representatives[c("category", "gene_id")]) ||
    any(!is.finite(representatives$saved_primary_lfsr)) ||
    any(representatives$saved_primary_lfsr < 0 |
        representatives$saved_primary_lfsr > 1)) {
  stop("The representative-pair universe failed validation.")
}

rows_by_grid <- lapply(grid_steps, function(step) grid_rows(fine_grid, step))
names(rows_by_grid) <- grid_labels
expected_grid_points <- as.integer(round(15 / grid_steps)) + 1L
if (!identical(unname(lengths(rows_by_grid)), expected_grid_points)) {
  stop("The requested grids are not valid nested subgrids of the 0.025 grid.")
}

load_start <- proc.time()[["elapsed"]]
load(fit_path)
fit_load_seconds <- proc.time()[["elapsed"]] - load_start
if (!exists("fash_fit1_update") || !inherits(fash_fit1_update, "fash")) {
  stop("The fitted-model file did not contain fash_fit1_update.")
}

representatives_by_index <- split(
  representatives,
  representatives$index,
  drop = TRUE
)
representative_indices <- as.integer(names(representatives_by_index))

sample_one_index <- function(pair_index) {
  current_representatives <- representatives_by_index[[as.character(pair_index)]]
  result_rows <- vector(
    "list",
    nrow(current_representatives) * length(posterior_draw_values) * length(grid_steps)
  )
  result_position <- 0L

  for (posterior_draws in posterior_draw_values) {
    set.seed(seed + pair_index)
    posterior_samples <- predict(
      fash_fit1_update,
      index = pair_index,
      smooth_var = fine_grid,
      only.samples = TRUE,
      M = posterior_draws
    )
    if (!is.matrix(posterior_samples) ||
        nrow(posterior_samples) != length(fine_grid) ||
        ncol(posterior_samples) != posterior_draws ||
        any(!is.finite(posterior_samples))) {
      stop("Posterior sampling failed for fitted pair index ", pair_index, ".")
    }

    for (grid_label in grid_labels) {
      keep_rows <- rows_by_grid[[grid_label]]
      evaluation_grid <- fine_grid[keep_rows]
      grid_samples <- posterior_samples[keep_rows, , drop = FALSE]

      for (representative_row in seq_len(nrow(current_representatives))) {
        current <- current_representatives[representative_row, , drop = FALSE]
        indicators <- functional_indicator(
          samples = grid_samples,
          evaluation_grid = evaluation_grid,
          category = current$category,
          switch_threshold = switch_threshold
        )
        lfsr <- mean(indicators)

        result_position <- result_position + 1L
        result_rows[[result_position]] <- data.frame(
          category = current$category,
          gene_id = current$gene_id,
          gene_symbol = current$gene_symbol,
          index = pair_index,
          pair_id = current$pair_id,
          variant_id = current$variant_id,
          saved_primary_lfsr = current$saved_primary_lfsr,
          posterior_draws = posterior_draws,
          grid_step = grid_label,
          grid_points = length(evaluation_grid),
          lfsr = lfsr,
          mcse = sqrt(lfsr * (1 - lfsr) / posterior_draws),
          stringsAsFactors = FALSE
        )
      }
    }
    rm(posterior_samples)
    gc()
  }

  do.call(rbind, result_rows)
}

sampling_start <- proc.time()[["elapsed"]]
if (num_cores > 1L) {
  sampled <- parallel::mclapply(
    representative_indices,
    sample_one_index,
    mc.cores = num_cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  sampled <- lapply(representative_indices, sample_one_index)
}
sampling_seconds <- proc.time()[["elapsed"]] - sampling_start

if (any(vapply(sampled, inherits, logical(1), "try-error"))) {
  stop("At least one posterior-sampling task failed.")
}

pair_lfsr <- do.call(rbind, sampled)
pair_lfsr <- pair_lfsr[
  order(
    match(pair_lfsr$category, category_order),
    pair_lfsr$gene_id,
    pair_lfsr$posterior_draws,
    match(pair_lfsr$grid_step, grid_labels)
  ),
  ,
  drop = FALSE
]
rownames(pair_lfsr) <- NULL

expected_rows <- nrow(representatives) *
  length(posterior_draw_values) *
  length(grid_steps)
if (nrow(pair_lfsr) != expected_rows ||
    anyDuplicated(pair_lfsr[c("category", "gene_id", "posterior_draws", "grid_step")]) ||
    any(!is.finite(pair_lfsr$lfsr)) ||
    any(pair_lfsr$lfsr < 0 | pair_lfsr$lfsr > 1)) {
  stop("The pair-level convergence results failed validation.")
}

# Grid-resolution summaries hold the Monte Carlo size at 30,000 and compare
# adjacent evaluation-grid steps within the same representative pairs.
grid_reference_rows <- pair_lfsr[
  pair_lfsr$posterior_draws == 30000L,
  c("category", "pair_id", "grid_step", "lfsr"),
  drop = FALSE
]
grid_wide <- reshape(
  grid_reference_rows,
  idvar = c("category", "pair_id"),
  timevar = "grid_step",
  direction = "wide"
)
names(grid_wide)[names(grid_wide) == "lfsr.0.10"] <- "lfsr_0p10"
names(grid_wide)[names(grid_wide) == "lfsr.0.05"] <- "lfsr_0p05"
names(grid_wide)[names(grid_wide) == "lfsr.0.025"] <- "lfsr_0p025"

grid_summary_rows <- list()
grid_summary_position <- 0L
for (category in category_order) {
  current <- grid_wide[grid_wide$category == category, , drop = FALSE]
  for (comparison in c("0.10 vs 0.05", "0.05 vs 0.025")) {
    comparison_columns <- if (comparison == "0.10 vs 0.05") {
      c("lfsr_0p10", "lfsr_0p05")
    } else {
      c("lfsr_0p05", "lfsr_0p025")
    }
    grid_summary_position <- grid_summary_position + 1L
    grid_summary_rows[[grid_summary_position]] <- summarize_pairwise_comparison(
      current = current,
      category = category,
      comparison = comparison,
      estimate_column = comparison_columns[2],
      reference_column = comparison_columns[1]
    )
  }
}
grid_summary <- do.call(rbind, grid_summary_rows)

# Monte Carlo summaries hold the grid step at 0.05 and compare each smaller
# run with the separate fixed-seed 30,000-sample result.
mc_reference_rows <- pair_lfsr[
  pair_lfsr$grid_step == "0.05",
  c("category", "pair_id", "posterior_draws", "lfsr"),
  drop = FALSE
]
mc_wide <- reshape(
  mc_reference_rows,
  idvar = c("category", "pair_id"),
  timevar = "posterior_draws",
  direction = "wide"
)
names(mc_wide)[names(mc_wide) == "lfsr.3000"] <- "lfsr_3000"
names(mc_wide)[names(mc_wide) == "lfsr.10000"] <- "lfsr_10000"
names(mc_wide)[names(mc_wide) == "lfsr.30000"] <- "lfsr_30000"

mc_summary_rows <- list()
mc_summary_position <- 0L
for (category in category_order) {
  current <- mc_wide[mc_wide$category == category, , drop = FALSE]
  for (comparison in c("3,000 vs 30,000", "10,000 vs 30,000")) {
    estimate_column <- if (comparison == "3,000 vs 30,000") {
      "lfsr_3000"
    } else {
      "lfsr_10000"
    }
    mc_summary_position <- mc_summary_position + 1L
    mc_summary_rows[[mc_summary_position]] <- summarize_pairwise_comparison(
      current = current,
      category = category,
      comparison = comparison,
      estimate_column = estimate_column,
      reference_column = "lfsr_30000"
    )
  }
}
mc_summary <- do.call(rbind, mc_summary_rows)

representative_counts <- data.frame(
  category = category_order,
  n_primary_pairs = as.integer(table(factor(
    all_assignments$category,
    levels = category_order
  ))),
  n_representative_genes = as.integer(observed_counts),
  stringsAsFactors = FALSE
)
runtime <- data.frame(
  stage = c("fit_load", "posterior_sampling_and_evaluation"),
  elapsed_seconds = c(fit_load_seconds, sampling_seconds),
  stringsAsFactors = FALSE
)
validation <- data.frame(
  check = c(
    "representative_counts_match",
    "expected_pair_level_rows",
    "unique_category_gene_setting_rows",
    "finite_bounded_lfsr",
    "complete_grid_summary",
    "complete_mc_summary"
  ),
  passed = c(
    identical(as.integer(observed_counts), unname(expected_representative_counts)),
    nrow(pair_lfsr) == expected_rows,
    !anyDuplicated(pair_lfsr[c("category", "gene_id", "posterior_draws", "grid_step")]),
    all(is.finite(pair_lfsr$lfsr) & pair_lfsr$lfsr >= 0 & pair_lfsr$lfsr <= 1),
    nrow(grid_summary) == length(category_order) * 2L,
    nrow(mc_summary) == length(category_order) * 2L
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  stop("At least one final cache validation failed.")
}

configuration <- list(
  fit_path = normalizePath(fit_path),
  category_paths = normalizePath(category_specs$path),
  gene_map_path = normalizePath(gene_map_path),
  output_id = output_id,
  seed = seed,
  num_cores = num_cores,
  alpha = alpha,
  switch_threshold = switch_threshold,
  category_order = category_order,
  expected_representative_counts = expected_representative_counts,
  grid_steps = grid_steps,
  posterior_draw_values = posterior_draw_values,
  fine_grid = fine_grid,
  sampling_construction = "separate_fixed_seed_calls",
  representative_rule = "minimum saved primary LFSR per category and gene",
  estimand = paste(
    "Functional-LFSR numerical convergence for one minimum-LFSR",
    "primary pair per category and gene."
  )
)

write.csv(representatives, file.path(staging_dir, "gene_representatives.csv"), row.names = FALSE)
write.csv(representative_counts, file.path(staging_dir, "representative_counts.csv"), row.names = FALSE)
write.csv(pair_lfsr, file.path(staging_dir, "pair_lfsr_by_grid_and_mc_size.csv"), row.names = FALSE)
write.csv(grid_summary, file.path(staging_dir, "grid_resolution_summary.csv"), row.names = FALSE)
write.csv(mc_summary, file.path(staging_dir, "mc_sample_summary.csv"), row.names = FALSE)
write.csv(runtime, file.path(staging_dir, "runtime.csv"), row.names = FALSE)
write.csv(validation, file.path(staging_dir, "validation.csv"), row.names = FALSE)
saveRDS(configuration, file.path(staging_dir, "configuration.rds"))

if (!file.rename(staging_dir, output_dir)) {
  stop("Failed to promote the validated staging directory to the final cache.")
}

cat("All-category representative convergence analysis completed.\n")
cat("Output: ", normalizePath(output_dir), "\n", sep = "")
cat("Representatives: ", nrow(representatives), "\n", sep = "")
cat("Sampling seconds: ", formatC(sampling_seconds, format = "f", digits = 3), "\n", sep = "")
print(representative_counts)
print(grid_summary)
print(mc_summary)
