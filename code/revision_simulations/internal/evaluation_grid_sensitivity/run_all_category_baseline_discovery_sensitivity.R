#!/usr/bin/env Rscript

# Evaluate numerical sensitivity for every pair selected at the baseline
# setting in Early, Middle, Late, or Switch. The Middle definition is selected
# explicitly as open 3-to-12 or original closed 4-to-11. Each fitted pair is
# sampled once at the finest grid and largest Monte Carlo size. Nested grid
# rows and draw prefixes isolate the requested numerical changes.

find_workflowr_root <- function() {
  helper <- file.path(
    "code", "revision_simulations", "internal",
    "evaluation_grid_sensitivity",
    "all_category_baseline_discovery_helpers.R"
  )
  if (file.exists(helper)) return(".")
  if (file.exists(file.path("coderepo-local", helper))) {
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

load_saved_category_discoveries <- function(path,
                                            object_name,
                                            category,
                                            alpha) {
  input_environment <- new.env(parent = emptyenv())
  load(path, envir = input_environment)
  if (!exists(object_name, envir = input_environment, inherits = FALSE)) {
    stop("Missing ", object_name, " in ", path, ".")
  }

  category_table <- input_environment[[object_name]]
  required_columns <- c("indices", "lfsr", "cfsr")
  if (!is.data.frame(category_table) ||
      !all(required_columns %in% names(category_table))) {
    stop("The saved category cache is incomplete: ", path)
  }

  selected <- category_table[category_table$cfsr <= alpha, , drop = FALSE]
  data.frame(
    category = category,
    pair_id = rownames(selected),
    index = as.integer(selected$indices),
    selection_lfsr = as.numeric(selected$lfsr),
    selection_cfsr = as.numeric(selected$cfsr),
    selection_source = paste0("saved_", category, "_classification"),
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "evaluation_grid_sensitivity",
  "all_category_baseline_discovery_helpers.R"
))

num_cores <- as.integer(get_arg("--num-cores", "2"))
seed <- as.integer(get_arg("--seed", "20260820"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
switch_threshold <- as.numeric(get_arg("--switch-threshold", "0.25"))
shuffle_seed_offset <- as.integer(get_arg("--shuffle-seed-offset", "1000000"))
middle_definition <- get_arg("--middle-definition", "open_3_12")
output_id <- get_arg(
  "--output-id",
  if (identical(middle_definition, "closed_4_11")) {
    "evaluation_grid_mc_all_category_baseline_discoveries_middle_4_11"
  } else {
    "evaluation_grid_mc_all_category_baseline_discoveries_open_middle"
  }
)

if (is.na(num_cores) || num_cores < 1L || num_cores > 4L ||
    is.na(seed) || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
    !is.finite(switch_threshold) || switch_threshold <= 0 ||
    is.na(shuffle_seed_offset) || shuffle_seed_offset < 1L ||
    !middle_definition %in% c("open_3_12", "closed_4_11") ||
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

category_order <- c("early", "middle", "late", "switch")
expected_pair_counts <- c(
  early = 124L,
  middle = if (middle_definition == "open_3_12") 60L else 24L,
  late = 20L,
  switch = 984L
)
expected_gene_counts <- c(
  early = 8L,
  middle = if (middle_definition == "open_3_12") 15L else 5L,
  late = 12L,
  switch = 250L
)

fit_path <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
)
gene_map_path <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real", "cache_gene_map.rds"
)
saved_category_specs <- data.frame(
  category = c("early", "late", "switch"),
  object_name = c(
    "testing_early_dyn",
    "testing_late_dyn",
    "testing_switch_dyn"
  ),
  file_name = c(
    "classify_dyn_eQTLs_early.RData",
    "classify_dyn_eQTLs_late.RData",
    "classify_dyn_eQTLs_switch.RData"
  ),
  stringsAsFactors = FALSE
)
saved_category_specs$path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  saved_category_specs$file_name
)
open_middle_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "middle_open_3_12_full_grid_mc_sensitivity",
  "pair_lfsr_by_setting.csv"
)
saved_middle_path <- file.path(
  workflowr_root,
  "output", "dynamic_eQTL_real",
  "classify_dyn_eQTLs_middle.RData"
)
output_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  output_id
)
staging_dir <- paste0(output_dir, ".staging-", Sys.getpid())

required_paths <- c(
  fit_path,
  gene_map_path,
  saved_category_specs$path,
  if (middle_definition == "open_3_12") {
    open_middle_path
  } else {
    saved_middle_path
  }
)
if (any(!file.exists(required_paths))) {
  stop("A required fitted-model or category-classification input is missing.")
}
if (dir.exists(output_dir) || dir.exists(staging_dir)) {
  stop("The requested output or staging directory already exists: ", output_dir)
}
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

saved_assignments <- do.call(rbind, lapply(
  seq_len(nrow(saved_category_specs)),
  function(i) {
    load_saved_category_discoveries(
      path = saved_category_specs$path[i],
      object_name = saved_category_specs$object_name[i],
      category = saved_category_specs$category[i],
      alpha = alpha
    )
  }
))

if (middle_definition == "open_3_12") {
  open_middle_all <- read.csv(
    open_middle_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  required_middle_columns <- c(
    "pair_id", "index", "lfsr", "cfsr", "setting_id", "selected"
  )
  if (!all(required_middle_columns %in% names(open_middle_all))) {
    stop("The open-Middle full-universe cache is incomplete.")
  }
  open_middle <- open_middle_all[
    open_middle_all$setting_id == "grid_0p10_M3000" &
      open_middle_all$selected,
    ,
    drop = FALSE
  ]
  middle_assignments <- data.frame(
    category = "middle",
    pair_id = open_middle$pair_id,
    index = as.integer(open_middle$index),
    selection_lfsr = as.numeric(open_middle$lfsr),
    selection_cfsr = as.numeric(open_middle$cfsr),
    selection_source = "open_middle_grid_0p10_M3000",
    stringsAsFactors = FALSE
  )
} else {
  middle_assignments <- load_saved_category_discoveries(
    path = saved_middle_path,
    object_name = "testing_middle_dyn",
    category = "middle",
    alpha = alpha
  )
  middle_assignments$selection_source <-
    "saved_middle_classification_4_11"
}

assignments <- rbind(saved_assignments, middle_assignments)
assignments$gene_id <- sub("_(rs[^_]+)$", "", assignments$pair_id)
assignments$variant_id <- sub("^.*_(rs[^_]+)$", "\\1", assignments$pair_id)
gene_map <- readRDS(gene_map_path)
assignments$gene_symbol <- gene_map$hgnc_symbol[
  match(assignments$gene_id, gene_map$ensembl_gene_id)
]
missing_symbol <- is.na(assignments$gene_symbol) |
  !nzchar(assignments$gene_symbol)
assignments$gene_symbol[missing_symbol] <-
  assignments$gene_id[missing_symbol]
assignments <- assignments[
  order(
    match(assignments$category, category_order),
    assignments$gene_id,
    assignments$selection_lfsr,
    assignments$pair_id
  ),
  ,
  drop = FALSE
]
rownames(assignments) <- NULL

observed_pair_counts <- table(factor(
  assignments$category,
  levels = category_order
))
observed_gene_counts <- vapply(
  category_order,
  function(category) {
    length(unique(assignments$gene_id[assignments$category == category]))
  },
  integer(1)
)
if (!identical(as.integer(observed_pair_counts), unname(expected_pair_counts)) ||
    !identical(unname(observed_gene_counts), unname(expected_gene_counts)) ||
    anyDuplicated(assignments[c("category", "pair_id")]) ||
    any(!is.finite(assignments$selection_lfsr)) ||
    any(!is.finite(assignments$selection_cfsr)) ||
    any(assignments$selection_lfsr < 0 | assignments$selection_lfsr > 1) ||
    any(assignments$selection_cfsr < 0 | assignments$selection_cfsr > alpha)) {
  stop("The baseline-discovery analysis population failed validation.")
}

pair_index_map <- unique(assignments[c("pair_id", "index")])
if (anyDuplicated(pair_index_map$pair_id) ||
    anyDuplicated(pair_index_map$index)) {
  stop("Pair identifiers and fitted-model indices are not one-to-one.")
}

settings <- data.frame(
  setting_id = c(
    "grid_0p15_M3000",
    "grid_0p10_M3000",
    "grid_0p05_M3000",
    "grid_0p10_M2000",
    "grid_0p10_M5000"
  ),
  grid_step = c(0.15, 0.10, 0.05, 0.10, 0.10),
  posterior_draws = c(3000L, 3000L, 3000L, 2000L, 5000L),
  role = c(
    "coarser_grid",
    "baseline",
    "finer_grid",
    "smaller_mc",
    "larger_mc"
  ),
  stringsAsFactors = FALSE
)
comparisons <- data.frame(
  comparison_id = c(
    "grid_0p10_vs_0p15_M3000",
    "grid_0p10_vs_0p05_M3000",
    "M2000_vs_M3000_grid_0p10",
    "M3000_vs_M5000_grid_0p10"
  ),
  comparison = c(
    "Grid 0.10 vs 0.15 at M = 3,000",
    "Grid 0.10 vs 0.05 at M = 3,000",
    "M = 2,000 vs 3,000 at grid 0.10",
    "M = 3,000 vs 5,000 at grid 0.10"
  ),
  x_setting_id = c(
    "grid_0p10_M3000",
    "grid_0p10_M3000",
    "grid_0p10_M2000",
    "grid_0p10_M3000"
  ),
  y_setting_id = c(
    "grid_0p15_M3000",
    "grid_0p05_M3000",
    "grid_0p10_M3000",
    "grid_0p10_M5000"
  ),
  stringsAsFactors = FALSE
)

fine_grid <- seq(0, 15, by = 0.05)
grid_values <- list(
  `0.15` = seq(0, 15, by = 0.15),
  `0.10` = seq(0, 15, by = 0.10),
  `0.05` = fine_grid
)
grid_rows <- lapply(
  grid_values,
  function(grid) match_subgrid_rows(fine_grid, grid)
)
if (!identical(unname(lengths(grid_values)), c(101L, 151L, 301L))) {
  stop("The evaluation grids have unexpected sizes.")
}

load_start <- proc.time()[["elapsed"]]
load(fit_path)
fit_load_seconds <- proc.time()[["elapsed"]] - load_start
if (!exists("fash_fit1_update") || !inherits(fash_fit1_update, "fash")) {
  stop("The fitted-model file did not contain fash_fit1_update.")
}

# Confirm exact grid reuse and document why Monte Carlo prefixes must be taken
# from one common 5,000-draw call. predict.fash() changes its random-number
# consumption when M changes, so a direct M = 3,000 call is not the prefix of
# a direct M = 5,000 call under the same seed.
audit_start <- proc.time()[["elapsed"]]
audit_indices <- head(pair_index_map$index, 3L)
direct_vs_nested_trajectory_maximum_difference <- 0
direct_vs_nested_lfsr_maximum_difference <- 0
subgrid_maximum_difference <- 0
for (pair_index in audit_indices) {
  set.seed(seed + pair_index)
  audit_5000 <- predict(
    fash_fit1_update,
    index = pair_index,
    smooth_var = fine_grid,
    only.samples = TRUE,
    M = 5000L
  )
  audit_5000 <- shuffle_posterior_draws(
    audit_5000,
    seed = seed + shuffle_seed_offset + pair_index
  )
  set.seed(seed + pair_index)
  audit_3000 <- predict(
    fash_fit1_update,
    index = pair_index,
    smooth_var = fine_grid,
    only.samples = TRUE,
    M = 3000L
  )
  set.seed(seed + pair_index)
  audit_direct_0p10 <- predict(
    fash_fit1_update,
    index = pair_index,
    smooth_var = grid_values[["0.10"]],
    only.samples = TRUE,
    M = 3000L
  )
  direct_vs_nested_trajectory_maximum_difference <- max(
    direct_vs_nested_trajectory_maximum_difference,
    abs(audit_5000[, seq_len(3000L), drop = FALSE] - audit_3000)
  )
  subgrid_maximum_difference <- max(
    subgrid_maximum_difference,
    abs(
      audit_3000[grid_rows[["0.10"]], , drop = FALSE] -
      audit_direct_0p10
    )
  )

  audit_categories <- assignments$category[assignments$index == pair_index]
  nested_0p10 <- audit_5000[
    grid_rows[["0.10"]],
    seq_len(3000L),
    drop = FALSE
  ]
  for (category in audit_categories) {
    direct_lfsr <- functional_lfsr(
      samples = audit_direct_0p10,
      evaluation_grid = grid_values[["0.10"]],
      category = category,
      switch_threshold = switch_threshold,
      middle_definition = middle_definition
    )
    nested_lfsr <- functional_lfsr(
      samples = nested_0p10,
      evaluation_grid = grid_values[["0.10"]],
      category = category,
      switch_threshold = switch_threshold,
      middle_definition = middle_definition
    )
    direct_vs_nested_lfsr_maximum_difference <- max(
      direct_vs_nested_lfsr_maximum_difference,
      abs(nested_lfsr - direct_lfsr)
    )
  }
  rm(audit_5000, audit_3000, audit_direct_0p10)
}
audit_seconds <- proc.time()[["elapsed"]] - audit_start
if (!is.finite(direct_vs_nested_trajectory_maximum_difference) ||
    !is.finite(direct_vs_nested_lfsr_maximum_difference) ||
    !is.finite(subgrid_maximum_difference) ||
    subgrid_maximum_difference > 1e-12) {
  stop("The grid-subsetting reuse audit failed.")
}

assignments_by_index <- split(assignments, assignments$index, drop = TRUE)
unique_pair_indices <- as.integer(names(assignments_by_index))

sample_one_pair <- function(pair_index) {
  current_assignments <- assignments_by_index[[as.character(pair_index)]]
  set.seed(seed + pair_index)
  posterior_samples <- predict(
    fash_fit1_update,
    index = pair_index,
    smooth_var = fine_grid,
    only.samples = TRUE,
    M = 5000L
  )
  if (!is.matrix(posterior_samples) ||
      nrow(posterior_samples) != length(fine_grid) ||
      ncol(posterior_samples) != 5000L ||
      any(!is.finite(posterior_samples))) {
    stop("Posterior sampling failed for fitted pair index ", pair_index, ".")
  }
  posterior_samples <- shuffle_posterior_draws(
    posterior_samples,
    seed = seed + shuffle_seed_offset + pair_index
  )

  result_rows <- vector(
    "list",
    nrow(current_assignments) * nrow(settings)
  )
  result_position <- 0L
  for (setting_position in seq_len(nrow(settings))) {
    grid_label <- sprintf("%.2f", settings$grid_step[setting_position])
    keep_rows <- grid_rows[[grid_label]]
    keep_draws <- seq_len(settings$posterior_draws[setting_position])
    evaluation_grid <- fine_grid[keep_rows]
    setting_samples <- posterior_samples[
      keep_rows,
      keep_draws,
      drop = FALSE
    ]

    for (assignment_position in seq_len(nrow(current_assignments))) {
      current <- current_assignments[assignment_position, , drop = FALSE]
      lfsr <- functional_lfsr(
        samples = setting_samples,
        evaluation_grid = evaluation_grid,
        category = current$category,
        switch_threshold = switch_threshold,
        middle_definition = middle_definition
      )
      result_position <- result_position + 1L
      result_rows[[result_position]] <- data.frame(
        category = current$category,
        pair_id = current$pair_id,
        gene_id = current$gene_id,
        gene_symbol = current$gene_symbol,
        variant_id = current$variant_id,
        index = pair_index,
        selection_lfsr = current$selection_lfsr,
        selection_cfsr = current$selection_cfsr,
        selection_source = current$selection_source,
        setting_id = settings$setting_id[setting_position],
        grid_step = settings$grid_step[setting_position],
        grid_points = length(evaluation_grid),
        posterior_draws = settings$posterior_draws[setting_position],
        lfsr = lfsr,
        mcse = sqrt(lfsr * (1 - lfsr) /
          settings$posterior_draws[setting_position]),
        stringsAsFactors = FALSE
      )
    }
  }
  result_rows
}

sampling_start <- proc.time()[["elapsed"]]
if (num_cores > 1L) {
  sampled <- parallel::mclapply(
    unique_pair_indices,
    sample_one_pair,
    mc.cores = num_cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  sampled <- lapply(unique_pair_indices, sample_one_pair)
}
sampling_seconds <- proc.time()[["elapsed"]] - sampling_start
if (any(vapply(sampled, inherits, logical(1), "try-error"))) {
  stop("At least one posterior-sampling task failed.")
}

pair_lfsr <- do.call(rbind, unlist(sampled, recursive = FALSE))
pair_lfsr <- pair_lfsr[
  order(
    match(pair_lfsr$category, category_order),
    pair_lfsr$gene_id,
    pair_lfsr$pair_id,
    match(pair_lfsr$setting_id, settings$setting_id)
  ),
  ,
  drop = FALSE
]
rownames(pair_lfsr) <- NULL

comparison_summaries <- list()
comparison_position <- 0L
for (comparison_row in seq_len(nrow(comparisons))) {
  current_comparison <- comparisons[comparison_row, , drop = FALSE]
  for (category in category_order) {
    x_rows <- pair_lfsr[
      pair_lfsr$category == category &
        pair_lfsr$setting_id == current_comparison$x_setting_id,
      ,
      drop = FALSE
    ]
    y_rows <- pair_lfsr[
      pair_lfsr$category == category &
        pair_lfsr$setting_id == current_comparison$y_setting_id,
      ,
      drop = FALSE
    ]
    comparison_position <- comparison_position + 1L
    current_summary <- summarize_reported_pair_comparison(
      reference = x_rows,
      comparison = y_rows,
      category = category,
      comparison_label = current_comparison$comparison
    )
    current_summary$comparison_id <- current_comparison$comparison_id
    current_summary$x_setting_id <- current_comparison$x_setting_id
    current_summary$y_setting_id <- current_comparison$y_setting_id
    comparison_summaries[[comparison_position]] <- current_summary
  }
}
comparison_summary <- do.call(rbind, comparison_summaries)
comparison_summary <- comparison_summary[
  order(
    match(comparison_summary$comparison_id, comparisons$comparison_id),
    match(comparison_summary$category, category_order)
  ),
  ,
  drop = FALSE
]
rownames(comparison_summary) <- NULL

baseline_rows <- pair_lfsr[
  pair_lfsr$setting_id == "grid_0p10_M3000",
  ,
  drop = FALSE
]
middle_baseline <- baseline_rows[baseline_rows$category == "middle", , drop = FALSE]
middle_selection_vs_nested_baseline_maximum_difference <- max(abs(
  middle_baseline$lfsr - middle_baseline$selection_lfsr
))

population_summary <- do.call(rbind, lapply(category_order, function(category) {
  current <- assignments[assignments$category == category, , drop = FALSE]
  data.frame(
    category = category,
    n_pairs = nrow(current),
    n_genes = length(unique(current$gene_id)),
    n_unique_fitted_indices = length(unique(current$index)),
    selection_source = paste(unique(current$selection_source), collapse = "; "),
    stringsAsFactors = FALSE
  )
}))

validation <- data.frame(
  check = c(
    "expected_category_pair_counts",
    "expected_category_gene_counts",
    "unique_category_pair_assignments",
    "expected_pair_setting_rows",
    "unique_category_pair_setting_rows",
    "finite_bounded_lfsr",
    "complete_category_comparison_summary",
    "grid_subsetting_reuse_audit",
    "nested_draw_prefixes_share_one_common_sample",
    "middle_selection_population_fixed"
  ),
  passed = c(
    identical(as.integer(observed_pair_counts), unname(expected_pair_counts)),
    identical(unname(observed_gene_counts), unname(expected_gene_counts)),
    anyDuplicated(assignments[c("category", "pair_id")]) == 0L,
    nrow(pair_lfsr) == nrow(assignments) * nrow(settings),
    anyDuplicated(pair_lfsr[c("category", "pair_id", "setting_id")]) == 0L,
    all(is.finite(pair_lfsr$lfsr)) &&
      all(pair_lfsr$lfsr >= 0 & pair_lfsr$lfsr <= 1),
    nrow(comparison_summary) == length(category_order) * nrow(comparisons),
    subgrid_maximum_difference <= 1e-12,
    all(settings$posterior_draws <= 5000L),
    identical(
      sort(middle_baseline$pair_id),
      sort(middle_assignments$pair_id)
    )
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$passed)) {
  failed_checks <- validation$check[!validation$passed]
  stop("Production validation failed: ", paste(failed_checks, collapse = ", "))
}

runtime <- data.frame(
  stage = c("fit_load", "reuse_audit", "posterior_sampling_and_evaluation"),
  elapsed_seconds = c(
    fit_load_seconds,
    audit_seconds,
    sampling_seconds
  ),
  stringsAsFactors = FALSE
)

configuration <- list(
  generated_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  seed = seed,
  num_cores = num_cores,
  alpha = alpha,
  switch_threshold = switch_threshold,
  middle_definition = middle_definition,
  category_order = category_order,
  category_definitions = c(
    early = "t <= 3",
    middle = if (middle_definition == "open_3_12") {
      "3 < t < 12"
    } else {
      "4 <= t <= 11"
    },
    late = "t >= 12",
    switch = "min(max positive, max negative) > 0.25"
  ),
  settings = settings,
  comparisons = comparisons,
  baseline_setting = "grid_0p10_M3000",
  fine_grid = fine_grid,
  expected_pair_counts = expected_pair_counts,
  expected_gene_counts = expected_gene_counts,
  category_pair_count = nrow(assignments),
  unique_pair_count = length(unique_pair_indices),
  posterior_draw_construction = paste(
    "One 5,000-draw call per fitted pair; columns are deterministically",
    "shuffled because predict.fash groups them by PSD component; M = 2,000",
    "and M = 3,000 use nested prefixes of the shuffled matrix."
  ),
  shuffle_seed_offset = shuffle_seed_offset,
  direct_vs_nested_trajectory_maximum_difference =
    direct_vs_nested_trajectory_maximum_difference,
  direct_vs_nested_lfsr_maximum_difference =
    direct_vs_nested_lfsr_maximum_difference,
  subgrid_maximum_difference = subgrid_maximum_difference,
  middle_selection_vs_nested_baseline_maximum_difference =
    middle_selection_vs_nested_baseline_maximum_difference,
  input_paths = required_paths,
  input_md5 = unname(tools::md5sum(required_paths))
)

write.csv(
  assignments,
  file.path(staging_dir, "baseline_discovery_population.csv"),
  row.names = FALSE
)
write.csv(
  population_summary,
  file.path(staging_dir, "population_summary.csv"),
  row.names = FALSE
)
write.csv(
  pair_lfsr,
  file.path(staging_dir, "pair_lfsr_by_setting.csv"),
  row.names = FALSE
)
write.csv(
  comparison_summary,
  file.path(staging_dir, "comparison_summary.csv"),
  row.names = FALSE
)
write.csv(
  runtime,
  file.path(staging_dir, "runtime.csv"),
  row.names = FALSE
)
write.csv(
  validation,
  file.path(staging_dir, "validation.csv"),
  row.names = FALSE
)
saveRDS(configuration, file.path(staging_dir, "configuration.rds"))

if (!file.rename(staging_dir, output_dir)) {
  stop("Could not atomically promote the validated cache directory.")
}

message("Wrote all-category baseline-discovery sensitivity cache to ", output_dir)
message(
  "Category-pair rows: ", nrow(assignments),
  "; unique fitted pairs: ", length(unique_pair_indices),
  "; sampling seconds: ", sprintf("%.1f", sampling_seconds)
)
