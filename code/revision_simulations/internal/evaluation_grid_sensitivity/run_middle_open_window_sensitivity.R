#!/usr/bin/env Rscript

# Recompute the Middle functional with the presentation-oriented open window
# 3 < t < 12 over the complete saved manuscript-classification universe.
# The current setting is grid step 0.10 with 3,000 posterior draws. Grid and
# Monte Carlo alternatives are evaluated with pair-specific fixed seeds.

find_workflowr_root <- function() {
  helper_path <- file.path(
    "code", "revision_simulations", "internal",
    "evaluation_grid_sensitivity", "middle_open_window_helpers.R"
  )
  if (file.exists(helper_path)) return(".")
  if (file.exists(file.path("coderepo-local", helper_path))) {
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

match_subgrid_rows <- function(fine_grid, requested_grid) {
  rows <- match(round(requested_grid, 10), round(fine_grid, 10))
  if (anyNA(rows)) stop("A requested grid is not contained in the fine grid.")
  rows
}

summarize_comparison <- function(pair_results,
                                 candidate_table,
                                 baseline_id,
                                 comparison_id,
                                 comparison_label) {
  baseline <- pair_results[pair_results$setting_id == baseline_id, , drop = FALSE]
  comparison <- pair_results[
    pair_results$setting_id == comparison_id,
    ,
    drop = FALSE
  ]
  comparison <- comparison[match(baseline$pair_id, comparison$pair_id), , drop = FALSE]
  if (anyNA(comparison$pair_id) || !identical(baseline$pair_id, comparison$pair_id)) {
    stop("A paired comparison is incomplete: ", comparison_label)
  }

  difference <- comparison$lfsr - baseline$lfsr
  absolute_difference <- abs(difference)
  baseline_pairs <- baseline$pair_id[baseline$selected]
  comparison_pairs <- comparison$pair_id[comparison$selected]
  baseline_genes <- unique(baseline$gene_id[baseline$selected])
  comparison_genes <- unique(comparison$gene_id[comparison$selected])
  old_pairs <- candidate_table$pair_id[candidate_table$saved_middle_selected]
  old_genes <- unique(candidate_table$gene_id[candidate_table$saved_middle_selected])

  data.frame(
    comparison = comparison_label,
    baseline_setting = baseline_id,
    comparison_setting = comparison_id,
    n_pairs = nrow(baseline),
    mean_signed_difference = mean(difference),
    mean_absolute_difference = mean(absolute_difference),
    median_absolute_difference = median(absolute_difference),
    q90_absolute_difference = unname(quantile(absolute_difference, 0.90)),
    maximum_absolute_difference = max(absolute_difference),
    spearman_correlation = suppressWarnings(cor(
      baseline$lfsr,
      comparison$lfsr,
      method = "spearman"
    )),
    fraction_within_0p005 = mean(absolute_difference <= 0.005),
    fraction_within_0p01 = mean(absolute_difference <= 0.01),
    selected_pair_agreement = mean(baseline$selected == comparison$selected),
    selected_pair_jaccard = set_jaccard(baseline_pairs, comparison_pairs),
    selected_gene_jaccard = set_jaccard(baseline_genes, comparison_genes),
    old_pair_retention_baseline = mean(old_pairs %in% baseline_pairs),
    old_pair_retention_comparison = mean(old_pairs %in% comparison_pairs),
    old_gene_retention_baseline = mean(old_genes %in% baseline_genes),
    old_gene_retention_comparison = mean(old_genes %in% comparison_genes),
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
helper_path <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "evaluation_grid_sensitivity", "middle_open_window_helpers.R"
)
source(helper_path)

num_cores <- as.integer(get_arg("--num-cores", "2"))
seed <- as.integer(get_arg("--seed", "20260820"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
output_id <- get_arg(
  "--output-id",
  "middle_open_3_12_full_grid_mc_sensitivity"
)

if (is.na(num_cores) || num_cores < 1L || num_cores > 4L ||
    is.na(seed) || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
    !nzchar(output_id)) {
  stop("Invalid sensitivity-analysis arguments.")
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
  workflowr_root, "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
)
saved_middle_path <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real",
  "classify_dyn_eQTLs_middle.RData"
)
gene_map_path <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real", "cache_gene_map.rds"
)
output_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "internal", output_id
)
staging_dir <- paste0(output_dir, ".staging-", Sys.getpid())

if (!all(file.exists(c(fit_path, saved_middle_path, gene_map_path)))) {
  stop("A required fitted-model or manuscript-classification input is missing.")
}
if (dir.exists(output_dir) || dir.exists(staging_dir)) {
  stop("The requested output or staging directory already exists: ", output_dir)
}
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

load_start <- proc.time()[["elapsed"]]
load(fit_path)
load(saved_middle_path)
gene_map <- readRDS(gene_map_path)
load_seconds <- proc.time()[["elapsed"]] - load_start

if (!exists("fash_fit1_update") || !inherits(fash_fit1_update, "fash") ||
    !exists("testing_middle_dyn") || !is.data.frame(testing_middle_dyn)) {
  stop("The expected fitted model or saved Middle table was not loaded.")
}
required_saved_columns <- c("indices", "lfsr", "cfsr")
if (!all(required_saved_columns %in% names(testing_middle_dyn))) {
  stop("The saved Middle table is incomplete.")
}

candidate_table <- data.frame(
  pair_id = rownames(testing_middle_dyn),
  index = as.integer(testing_middle_dyn$indices),
  saved_middle_4_11_lfsr = as.numeric(testing_middle_dyn$lfsr),
  saved_middle_4_11_cfsr = as.numeric(testing_middle_dyn$cfsr),
  stringsAsFactors = FALSE
)
candidate_table$saved_middle_selected <-
  candidate_table$saved_middle_4_11_cfsr <= alpha
candidate_table$gene_id <- sub("_(rs[^_]+)$", "", candidate_table$pair_id)
candidate_table$variant_id <- sub("^.*_(rs[^_]+)$", "\\1", candidate_table$pair_id)
candidate_table$gene_symbol <- gene_map$hgnc_symbol[
  match(candidate_table$gene_id, gene_map$ensembl_gene_id)
]
missing_symbol <- is.na(candidate_table$gene_symbol) |
  !nzchar(candidate_table$gene_symbol)
candidate_table$gene_symbol[missing_symbol] <-
  candidate_table$gene_id[missing_symbol]

expected_candidate_count <- 9205L
expected_saved_pair_count <- 24L
expected_saved_gene_count <- 5L
saved_genes <- unique(candidate_table$gene_id[candidate_table$saved_middle_selected])
if (nrow(candidate_table) != expected_candidate_count ||
    anyDuplicated(candidate_table$pair_id) ||
    anyDuplicated(candidate_table$index) ||
    sum(candidate_table$saved_middle_selected) != expected_saved_pair_count ||
    length(saved_genes) != expected_saved_gene_count ||
    any(!is.finite(candidate_table$saved_middle_4_11_lfsr)) ||
    any(!is.finite(candidate_table$saved_middle_4_11_cfsr))) {
  stop("The saved Middle-classification universe failed validation.")
}

settings <- data.frame(
  setting_id = c(
    "grid_0p15_M3000",
    "grid_0p10_M3000",
    "grid_0p05_M3000",
    "grid_0p10_M10000",
    "grid_0p10_M30000"
  ),
  grid_step = c(0.15, 0.10, 0.05, 0.10, 0.10),
  posterior_draws = c(3000L, 3000L, 3000L, 10000L, 30000L),
  role = c(
    "coarser_grid",
    "current_baseline",
    "finer_grid",
    "larger_mc",
    "largest_mc"
  ),
  stringsAsFactors = FALSE
)
baseline_id <- "grid_0p10_M3000"

fine_grid <- seq(0, 15, by = 0.05)
grid_values <- list(
  `0.15` = seq(0, 15, by = 0.15),
  `0.10` = seq(0, 15, by = 0.10),
  `0.05` = fine_grid
)
grid_rows <- lapply(grid_values, function(x) match_subgrid_rows(fine_grid, x))
if (!identical(unname(lengths(grid_values)), c(101L, 151L, 301L))) {
  stop("The evaluation grids have unexpected sizes.")
}

# Verify that sampling on the 0.05 grid and retaining every second row is
# numerically equivalent to sampling directly on the 0.10 grid with the same
# pair-specific seed. This permits exact reuse of the current 3,000-draw result.
audit_start <- proc.time()[["elapsed"]]
audit_indices <- head(candidate_table$index, 3L)
audit_maximum_difference <- 0
for (pair_index in audit_indices) {
  set.seed(seed + pair_index)
  audit_fine <- predict(
    fash_fit1_update,
    index = pair_index,
    smooth_var = fine_grid,
    only.samples = TRUE,
    M = 3000L
  )
  set.seed(seed + pair_index)
  audit_direct <- predict(
    fash_fit1_update,
    index = pair_index,
    smooth_var = grid_values[["0.10"]],
    only.samples = TRUE,
    M = 3000L
  )
  audit_maximum_difference <- max(
    audit_maximum_difference,
    abs(audit_fine[grid_rows[["0.10"]], , drop = FALSE] - audit_direct)
  )
  rm(audit_fine, audit_direct)
}
audit_seconds <- proc.time()[["elapsed"]] - audit_start
if (!is.finite(audit_maximum_difference) || audit_maximum_difference > 1e-12) {
  stop("The fine-grid reuse audit failed.")
}

sample_one_pair <- function(pair_index) {
  output <- matrix(NA_real_, nrow = nrow(settings), ncol = 3L)
  colnames(output) <- c("index", "setting_position", "lfsr")
  output[, "index"] <- pair_index
  output[, "setting_position"] <- seq_len(nrow(settings))

  set.seed(seed + pair_index)
  samples_3000 <- predict(
    fash_fit1_update,
    index = pair_index,
    smooth_var = fine_grid,
    only.samples = TRUE,
    M = 3000L
  )
  if (!is.matrix(samples_3000) || nrow(samples_3000) != length(fine_grid) ||
      ncol(samples_3000) != 3000L || any(!is.finite(samples_3000))) {
    stop("The 3,000-draw posterior sample failed for pair index ", pair_index)
  }
  for (setting_position in 1:3) {
    grid_label <- c("0.15", "0.10", "0.05")[setting_position]
    rows <- grid_rows[[grid_label]]
    output[setting_position, "lfsr"] <- middle_open_window_lfsr(
      samples_3000[rows, , drop = FALSE],
      grid_values[[grid_label]]
    )
  }
  rm(samples_3000)

  for (setting_position in 4:5) {
    posterior_draws <- settings$posterior_draws[setting_position]
    set.seed(seed + pair_index)
    samples <- predict(
      fash_fit1_update,
      index = pair_index,
      smooth_var = grid_values[["0.10"]],
      only.samples = TRUE,
      M = posterior_draws
    )
    if (!is.matrix(samples) || nrow(samples) != length(grid_values[["0.10"]]) ||
        ncol(samples) != posterior_draws || any(!is.finite(samples))) {
      stop(
        "The ", posterior_draws,
        "-draw posterior sample failed for pair index ", pair_index
      )
    }
    output[setting_position, "lfsr"] <- middle_open_window_lfsr(
      samples,
      grid_values[["0.10"]]
    )
    rm(samples)
  }

  output
}

sampling_start <- proc.time()[["elapsed"]]
if (num_cores > 1L) {
  sampled <- parallel::mclapply(
    candidate_table$index,
    sample_one_pair,
    mc.cores = num_cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  sampled <- lapply(candidate_table$index, sample_one_pair)
}
sampling_seconds <- proc.time()[["elapsed"]] - sampling_start
if (any(vapply(sampled, inherits, logical(1), "try-error"))) {
  stop("At least one pair-level posterior-sampling task failed.")
}

sample_matrix <- do.call(rbind, sampled)
pair_results <- data.frame(
  index = as.integer(sample_matrix[, "index"]),
  setting_position = as.integer(sample_matrix[, "setting_position"]),
  lfsr = as.numeric(sample_matrix[, "lfsr"]),
  stringsAsFactors = FALSE
)
pair_results$setting_id <- settings$setting_id[pair_results$setting_position]
candidate_match <- match(pair_results$index, candidate_table$index)
pair_results <- cbind(
  candidate_table[
    candidate_match,
    c(
      "pair_id", "gene_id", "gene_symbol", "variant_id",
      "saved_middle_4_11_lfsr", "saved_middle_4_11_cfsr",
      "saved_middle_selected"
    ),
    drop = FALSE
  ],
  pair_results,
  stringsAsFactors = FALSE
)
pair_results$grid_step <- settings$grid_step[pair_results$setting_position]
pair_results$posterior_draws <-
  settings$posterior_draws[pair_results$setting_position]
pair_results$cfsr <- NA_real_
pair_results$selected <- FALSE

for (setting_id in settings$setting_id) {
  rows <- which(pair_results$setting_id == setting_id)
  pair_results$cfsr[rows] <- deterministic_cumulative_fsr(
    pair_results$lfsr[rows],
    pair_results$pair_id[rows]
  )
  pair_results$selected[rows] <- pair_results$cfsr[rows] <= alpha
}
pair_results <- pair_results[
  order(pair_results$setting_position, pair_results$pair_id),
  ,
  drop = FALSE
]
rownames(pair_results) <- NULL

setting_summary_rows <- lapply(seq_len(nrow(settings)), function(i) {
  current <- pair_results[
    pair_results$setting_id == settings$setting_id[i],
    ,
    drop = FALSE
  ]
  selected_pairs <- current$pair_id[current$selected]
  selected_genes <- unique(current$gene_id[current$selected])
  old_pairs <- candidate_table$pair_id[candidate_table$saved_middle_selected]
  data.frame(
    settings[i, , drop = FALSE],
    candidate_pair_count = nrow(current),
    selected_pair_count = length(selected_pairs),
    selected_gene_count = length(selected_genes),
    old_selected_pair_count = length(old_pairs),
    old_selected_gene_count = length(saved_genes),
    old_pair_retained_count = sum(old_pairs %in% selected_pairs),
    old_pair_retention_fraction = mean(old_pairs %in% selected_pairs),
    old_gene_retained_count = sum(saved_genes %in% selected_genes),
    old_gene_retention_fraction = mean(saved_genes %in% selected_genes),
    selected_pair_jaccard_with_old = set_jaccard(selected_pairs, old_pairs),
    selected_gene_jaccard_with_old = set_jaccard(selected_genes, saved_genes),
    stringsAsFactors = FALSE
  )
})
setting_summary <- do.call(rbind, setting_summary_rows)

comparison_spec <- data.frame(
  comparison_id = c(
    "grid_0p15_M3000",
    "grid_0p05_M3000",
    "grid_0p10_M10000",
    "grid_0p10_M30000"
  ),
  comparison_label = c(
    "Grid 0.10 vs 0.15 at M = 3,000",
    "Grid 0.10 vs 0.05 at M = 3,000",
    "M = 3,000 vs 10,000 at grid 0.10",
    "M = 3,000 vs 30,000 at grid 0.10"
  ),
  stringsAsFactors = FALSE
)
comparison_summary <- do.call(rbind, lapply(seq_len(nrow(comparison_spec)), function(i) {
  summarize_comparison(
    pair_results = pair_results,
    candidate_table = candidate_table,
    baseline_id = baseline_id,
    comparison_id = comparison_spec$comparison_id[i],
    comparison_label = comparison_spec$comparison_label[i]
  )
}))

current_discovery_retention <- pair_results[
  pair_results$saved_middle_selected,
  c(
    "setting_id", "grid_step", "posterior_draws", "pair_id", "gene_id",
    "gene_symbol", "variant_id", "saved_middle_4_11_lfsr",
    "saved_middle_4_11_cfsr", "lfsr", "cfsr", "selected"
  ),
  drop = FALSE
]

expected_result_rows <- expected_candidate_count * nrow(settings)
validation <- data.frame(
  check = c(
    "expected_candidate_universe",
    "expected_saved_discoveries",
    "expected_pair_setting_rows",
    "unique_pair_setting_rows",
    "finite_bounded_lfsr",
    "finite_bounded_cfsr",
    "complete_setting_summary",
    "complete_comparison_summary",
    "fine_grid_reuse_matches_direct_grid"
  ),
  passed = c(
    nrow(candidate_table) == expected_candidate_count,
    sum(candidate_table$saved_middle_selected) == expected_saved_pair_count &
      length(saved_genes) == expected_saved_gene_count,
    nrow(pair_results) == expected_result_rows,
    !anyDuplicated(pair_results[c("pair_id", "setting_id")]),
    all(is.finite(pair_results$lfsr) &
          pair_results$lfsr >= 0 & pair_results$lfsr <= 1),
    all(is.finite(pair_results$cfsr) &
          pair_results$cfsr >= 0 & pair_results$cfsr <= 1),
    nrow(setting_summary) == nrow(settings),
    nrow(comparison_summary) == nrow(comparison_spec),
    audit_maximum_difference <= 1e-12
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  stop("At least one final cache validation failed.")
}

runtime <- data.frame(
  stage = c("input_load", "fine_grid_reuse_audit", "posterior_sampling"),
  elapsed_seconds = c(load_seconds, audit_seconds, sampling_seconds),
  stringsAsFactors = FALSE
)
configuration <- list(
  fit_path = normalizePath(fit_path),
  saved_middle_path = normalizePath(saved_middle_path),
  gene_map_path = normalizePath(gene_map_path),
  output_id = output_id,
  seed = seed,
  num_cores = num_cores,
  alpha = alpha,
  candidate_count = expected_candidate_count,
  saved_middle_pair_count = expected_saved_pair_count,
  saved_middle_gene_count = expected_saved_gene_count,
  middle_window = c(lower_open = 3, upper_open = 12),
  settings = settings,
  baseline_setting = baseline_id,
  sampling_construction = "separate_fixed_seed_calls_by_M",
  grid_construction = paste(
    "One M=3000 call on grid 0.05; exact row subsets for grid 0.15 and 0.10."
  ),
  estimand = paste(
    "Full-universe numerical sensitivity and current-discovery retention for",
    "the presentation-oriented Middle functional with 3 < t < 12."
  )
)

write.csv(
  pair_results,
  file.path(staging_dir, "pair_lfsr_by_setting.csv"),
  row.names = FALSE
)
write.csv(
  setting_summary,
  file.path(staging_dir, "setting_selection_summary.csv"),
  row.names = FALSE
)
write.csv(
  comparison_summary,
  file.path(staging_dir, "comparison_summary.csv"),
  row.names = FALSE
)
write.csv(
  current_discovery_retention,
  file.path(staging_dir, "current_discovery_retention.csv"),
  row.names = FALSE
)
write.csv(runtime, file.path(staging_dir, "runtime.csv"), row.names = FALSE)
write.csv(validation, file.path(staging_dir, "validation.csv"), row.names = FALSE)
saveRDS(configuration, file.path(staging_dir, "configuration.rds"))

if (!file.rename(staging_dir, output_dir)) {
  stop("Failed to promote the validated staging directory to the final cache.")
}

cat("Open-Middle-window sensitivity analysis completed.\n")
cat("Output: ", normalizePath(output_dir), "\n", sep = "")
cat("Candidates: ", nrow(candidate_table), "\n", sep = "")
cat("Sampling seconds: ", formatC(sampling_seconds, format = "f", digits = 3), "\n", sep = "")
print(setting_summary)
print(comparison_summary)
