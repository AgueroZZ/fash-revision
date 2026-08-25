#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

find_workflowr_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "_workflowr.yml"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the workflowr project root.")
    }
    current <- parent
  }
}

atomic_save_rds <- function(object, path) {
  temporary_path <- paste0(path, ".tmp")
  saveRDS(object, temporary_path)
  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically replace ", path)
  }
  invisible(path)
}

load_exact_object <- function(path, object_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  if (!identical(loaded, object_name)) {
    stop("Expected only object ", object_name, " in ", path)
  }
  environment[[object_name]]
}

extract_null_weight <- function(fit) {
  prior <- fit$prior_weights
  if (!is.data.frame(prior) || !"prior_weight" %in% names(prior)) {
    stop("The retained fit has no valid prior-weight table.")
  }
  if ("psd" %in% names(prior)) {
    index <- which(prior$psd == 0)
  } else if ("component" %in% names(prior)) {
    index <- which(prior$component == "constant")
  } else {
    index <- 1L
  }
  if (length(index) != 1L || !is.finite(prior$prior_weight[index])) {
    stop("The null prior weight could not be identified uniquely.")
  }
  as.numeric(prior$prior_weight[index])
}

make_input_provenance <- function(paths) {
  information <- file.info(paths)
  data.frame(
    label = names(paths),
    path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(information$size),
    modified = format(information$mtime, tz = "UTC", usetz = TRUE),
    md5 = unname(tools::md5sum(paths)),
    stringsAsFactors = FALSE
  )
}

make_discovery_row <- function(method, adjustment, lfdr, pair_table, alpha) {
  indices <- select_cumulative_lfdr_calls_linear(lfdr, alpha)
  selected <- pair_table[indices, , drop = FALSE]
  list(
    row = data.frame(
      method = method,
      adjustment = adjustment,
      pair_count = nrow(selected),
      gene_count = length(unique(selected$gene_id)),
      variant_count = length(unique(selected$variant_id)),
      stringsAsFactors = FALSE
    ),
    indices = indices,
    table = selected
  )
}

make_overlap_row <- function(unit, current_set, linear_set) {
  intersection_count <- length(intersect(current_set, linear_set))
  union_count <- length(union(current_set, linear_set))
  data.frame(
    unit = unit,
    current_count = length(current_set),
    linear_count = length(linear_set),
    intersection_count = intersection_count,
    current_only_count = length(setdiff(current_set, linear_set)),
    linear_only_count = length(setdiff(linear_set, current_set)),
    union_count = union_count,
    jaccard = if (union_count == 0L) NA_real_ else {
      intersection_count / union_count
    },
    stringsAsFactors = FALSE
  )
}

read_valid_statistics_cache <- function(path,
                                        dataset_md5,
                                        expected_time,
                                        expected_pairs) {
  if (!file.exists(path)) {
    return(NULL)
  }
  cache <- readRDS(path)
  required_columns <- c(
    "unit_id", "sum_w", "sum_wx", "sum_wxx", "sum_wy", "sum_wxy",
    "sum_wyy", "logdet_d"
  )
  valid <- is.list(cache) &&
    identical(cache$dataset_md5, dataset_md5) &&
    identical(cache$expected_time, expected_time) &&
    identical(cache$scale_time, TRUE) &&
    identical(cache$ridge, 1e-10) &&
    is.data.frame(cache$statistics) &&
    nrow(cache$statistics) == expected_pairs &&
    all(required_columns %in% names(cache$statistics)) &&
    !anyDuplicated(cache$statistics$unit_id) &&
    all(nzchar(cache$statistics$unit_id))
  if (!valid) NULL else cache
}

validate_matched_iwp_fit <- function(fit,
                                     pair_keys,
                                     grid,
                                     pred_step,
                                     penalty,
                                     label) {
  if (!inherits(fit, "fash") ||
      !identical(names(fit$fash_data$data_list), pair_keys) ||
      !isTRUE(all.equal(fit$psd_grid, grid, tolerance = 0)) ||
      !is.list(fit$settings) ||
      !isTRUE(all.equal(fit$settings$pred_step, pred_step, tolerance = 0)) ||
      !identical(as.integer(fit$settings$penalty), as.integer(penalty)) ||
      !is.data.frame(fit$prior_weights) ||
      fit$prior_weights$psd[1] != 0 ||
      length(fit$lfdr) != length(pair_keys) ||
      any(!is.finite(fit$lfdr)) ||
      any(fit$lfdr < 0 | fit$lfdr > 1)) {
    stop(label, " does not satisfy the matched IWP1 contract.")
  }
  invisible(TRUE)
}

make_prior_weight_rows <- function(fit, method, adjustment, grid) {
  weight <- expand_grid_prior_weights(fit$prior_weights, grid)
  data.frame(
    method = method,
    adjustment = adjustment,
    predstep_sd = grid,
    prior_weight = unname(weight),
    is_null = grid == 0,
    active = unname(weight) > 0,
    stringsAsFactors = FALSE
  )
}

make_prior_summary_row <- function(prior_rows) {
  null_weight <- prior_rows$prior_weight[prior_rows$is_null]
  alternative <- prior_rows[!prior_rows$is_null, , drop = FALSE]
  alternative_weight <- sum(alternative$prior_weight)
  data.frame(
    method = prior_rows$method[1],
    adjustment = prior_rows$adjustment[1],
    null_weight = null_weight,
    active_nonnull_components = sum(alternative$active),
    alternative_rms_predstep_sd = if (alternative_weight > 0) {
      sqrt(sum(
        alternative$prior_weight / alternative_weight *
          alternative$predstep_sd^2
      ))
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

compare_linear_versions <- function(adjustment,
                                    profile_lfdr,
                                    mixture_lfdr,
                                    pair_table,
                                    alpha) {
  profile_indices <- select_cumulative_lfdr_calls_linear(profile_lfdr, alpha)
  mixture_indices <- select_cumulative_lfdr_calls_linear(mixture_lfdr, alpha)
  profile_table <- pair_table[profile_indices, , drop = FALSE]
  mixture_table <- pair_table[mixture_indices, , drop = FALSE]
  intersection_count <- length(intersect(profile_indices, mixture_indices))
  union_count <- length(union(profile_indices, mixture_indices))
  data.frame(
    adjustment = adjustment,
    profile_pair_count = nrow(profile_table),
    mixture_pair_count = nrow(mixture_table),
    pair_count_difference = nrow(mixture_table) - nrow(profile_table),
    profile_gene_count = length(unique(profile_table$gene_id)),
    mixture_gene_count = length(unique(mixture_table$gene_id)),
    profile_variant_count = length(unique(profile_table$variant_id)),
    mixture_variant_count = length(unique(mixture_table$variant_id)),
    intersection_pair_count = intersection_count,
    profile_only_pair_count = length(setdiff(profile_indices, mixture_indices)),
    mixture_only_pair_count = length(setdiff(mixture_indices, profile_indices)),
    pair_jaccard = if (union_count == 0L) NA_real_ else {
      intersection_count / union_count
    },
    pearson_lfdr = stats::cor(profile_lfdr, mixture_lfdr),
    spearman_lfdr = stats::cor(profile_lfdr, mixture_lfdr, method = "spearman"),
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
analysis_directory <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "fash_linear_real_data_ablation"
)
helper_path <- file.path(
  analysis_directory,
  "fash_linear_real_data_helpers.R"
)
runner_path <- file.path(
  analysis_directory,
  "run_fash_linear_real_data_ablation.R"
)
shared_path <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "shared",
  "simulation_functions.R"
)
source(shared_path)
source(helper_path)
if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required for posterior trajectory prediction.")
}

input_directory <- file.path(workflowr_root, "output", "dynamic_eQTL_real")
datasets_path <- file.path(input_directory, "datasets_corrected.RData")
current_raw_path <- file.path(input_directory, "fash_fit1_all.RData")
current_bf_path <- file.path(input_directory, "fash_fit1_update.RData")
gene_map_path <- file.path(input_directory, "cache_gene_map.rds")
historical_output_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fash_linear_real_data_ablation"
)
historical_statistics_path <- file.path(
  historical_output_directory,
  "sufficient_statistics.rds"
)
historical_profile_raw_path <- file.path(
  historical_output_directory,
  "linear_fit_raw.rds"
)
historical_profile_bf_path <- file.path(
  historical_output_directory,
  "linear_fit_bf.rds"
)
historical_cache_path <- file.path(
  historical_output_directory,
  "analysis_cache.rds"
)
required_input_paths <- c(
  datasets = datasets_path,
  current_iwp1_raw = current_raw_path,
  current_iwp1_bf = current_bf_path,
  gene_map = gene_map_path,
  historical_profile_raw = historical_profile_raw_path,
  historical_profile_bf = historical_profile_bf_path,
  historical_profile_cache = historical_cache_path,
  helper = helper_path,
  shared_functions = shared_path,
  runner = runner_path
)
if (any(!file.exists(required_input_paths))) {
  stop("At least one required input is missing.")
}
input_paths <- required_input_paths
if (file.exists(historical_statistics_path)) {
  input_paths <- c(
    input_paths,
    historical_sufficient_statistics = historical_statistics_path
  )
}
input_provenance <- make_input_provenance(input_paths)

output_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fash_linear_real_data_ablation_mixture_predstep1_penalty10"
)
block_directory <- file.path(output_directory, "sufficient_statistic_blocks")
dir.create(block_directory, recursive = TRUE, showWarnings = FALSE)

analysis_started <- proc.time()[["elapsed"]]
analysis_id <- paste0(
  "revision_internal_fash_linear_real_data_ablation_",
  "mixture_predstep1_penalty10"
)
alpha <- 0.05
expected_pairs <- 1009173L
expected_time <- 0:15
block_size <- 25000L
grid <- default_revision_grid()
pred_step <- 1
penalty <- 10L
statistic_time_span <- diff(range(expected_time))
dataset_md5 <- input_provenance$md5[input_provenance$label == "datasets"]
sufficient_statistics_path <- file.path(
  output_directory,
  "sufficient_statistics.rds"
)
run_status_path <- file.path(output_directory, "run_status.rds")
atomic_save_rds(list(
  analysis_id = analysis_id,
  status = "in_progress",
  started_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
), run_status_path)

message("[1/7] Loading validated weighted sufficient statistics.")
statistics_cache <- NULL
statistics_cache_source <- NA_character_
for (candidate in c(sufficient_statistics_path, historical_statistics_path)) {
  candidate_cache <- read_valid_statistics_cache(
    path = candidate,
    dataset_md5 = dataset_md5,
    expected_time = expected_time,
    expected_pairs = expected_pairs
  )
  if (!is.null(candidate_cache)) {
    statistics_cache <- candidate_cache
    statistics_cache_source <- normalizePath(
      candidate,
      winslash = "/",
      mustWork = TRUE
    )
    break
  }
}

if (is.null(statistics_cache)) {
  message("  No valid retained cache found; extracting from the datasets.")
  dataset_environment <- new.env(parent = emptyenv())
  loaded <- load(datasets_path, envir = dataset_environment)
  if (!identical(loaded, "datasets")) {
    stop("The retained dataset file must contain only datasets.")
  }
  datasets <- dataset_environment$datasets
  rm(dataset_environment)
  if (length(datasets) != expected_pairs ||
      is.null(names(datasets)) || any(names(datasets) == "") ||
      anyDuplicated(names(datasets))) {
    stop("The retained datasets failed pair-count or key validation.")
  }
  pair_keys <- names(datasets)
  block_starts <- seq.int(1L, expected_pairs, by = block_size)
  block_paths <- character(length(block_starts))
  for (block_index in seq_along(block_starts)) {
    start <- block_starts[block_index]
    end <- min(start + block_size - 1L, expected_pairs)
    block_path <- file.path(
      block_directory,
      sprintf("block_%03d.rds", block_index)
    )
    block_paths[block_index] <- block_path
    valid_existing <- FALSE
    if (file.exists(block_path)) {
      block <- readRDS(block_path)
      valid_existing <- is.list(block) &&
        identical(block$dataset_md5, dataset_md5) &&
        identical(block$start, start) && identical(block$end, end) &&
        identical(block$expected_time, expected_time) &&
        identical(block$scale_time, TRUE) &&
        identical(block$ridge, 1e-10) &&
        is.data.frame(block$statistics) &&
        nrow(block$statistics) == end - start + 1L &&
        identical(block$statistics$unit_id, pair_keys[start:end]) &&
        all(c(
          "sum_w", "sum_wx", "sum_wxx", "sum_wy", "sum_wxy",
          "sum_wyy", "logdet_d"
        ) %in% names(block$statistics)) &&
        all(is.finite(as.matrix(block$statistics[, setdiff(
          names(block$statistics),
          "unit_id"
        ), drop = FALSE])))
    }
    if (!valid_existing) {
      block <- list(
        dataset_md5 = dataset_md5,
        start = start,
        end = end,
        expected_time = expected_time,
        scale_time = TRUE,
        ridge = 1e-10,
        statistics = compute_linear_sufficient_statistics_block(
          datasets[start:end],
          expected_time = expected_time,
          scale_time = TRUE
        )
      )
      atomic_save_rds(block, block_path)
    }
    message(
      "  sufficient statistics: ", end, "/", expected_pairs,
      if (valid_existing) " (cached)" else ""
    )
  }
  statistics <- do.call(rbind, lapply(block_paths, function(path) {
    readRDS(path)$statistics
  }))
  rownames(statistics) <- NULL
  rm(datasets)
  invisible(gc())
  statistics_cache_source <- "authoritative datasets"
} else {
  statistics <- statistics_cache$statistics
  pair_keys <- statistics$unit_id
  message("  Reused: ", statistics_cache_source)
}
if (nrow(statistics) != expected_pairs ||
    !identical(statistics$unit_id, pair_keys)) {
  stop("The assembled sufficient statistics are not aligned.")
}
statistics_cache <- list(
  analysis_id = analysis_id,
  dataset_md5 = dataset_md5,
  expected_time = expected_time,
  scale_time = TRUE,
  ridge = 1e-10,
  statistics = statistics
)
atomic_save_rds(statistics_cache, sufficient_statistics_path)

message("[2/7] Fitting the matched predictive-SD linear mixture.")
linear_fit_raw <- fit_linear_mixture_fash_from_stats(
  statistics = statistics,
  grid = grid,
  pred_step = pred_step,
  penalty = penalty,
  statistic_time_span = statistic_time_span,
  n_time = length(expected_time)
)
validate_linear_mixture_fash(
  linear_fit_raw,
  expected_grid = grid,
  expected_pred_step = pred_step,
  expected_penalty = penalty
)
linear_raw_path <- file.path(output_directory, "linear_fit_raw.rds")
linear_bf_path <- file.path(output_directory, "linear_fit_bf.rds")
linear_fit_raw_compact <- compact_linear_mixture_fash(linear_fit_raw)
validate_compact_linear_mixture_fash(
  linear_fit_raw_compact,
  expected_grid = grid,
  expected_pred_step = pred_step,
  expected_penalty = penalty
)
atomic_save_rds(
  linear_fit_raw_compact,
  linear_raw_path
)
linear_raw_lfdr <- as.numeric(linear_fit_raw$lfdr)
linear_prior_raw <- make_prior_weight_rows(
  linear_fit_raw,
  "FASH-linear",
  "Raw",
  grid
)

linear_fit_bf <- BF_update_linear_mixture_fash(linear_fit_raw)
rm(linear_fit_raw, linear_fit_raw_compact)
invisible(gc())

linear_bf_lfdr <- as.numeric(linear_fit_bf$lfdr)
linear_prior_bf <- make_prior_weight_rows(
  linear_fit_bf,
  "FASH-linear",
  "BF-adjusted",
  grid
)
linear_fit_bf_compact <- compact_linear_mixture_fash(linear_fit_bf)
rm(linear_fit_bf)
invisible(gc())
validate_compact_linear_mixture_fash(
  linear_fit_bf_compact,
  expected_grid = grid,
  expected_pred_step = pred_step,
  expected_penalty = penalty
)
atomic_save_rds(linear_fit_bf_compact, linear_bf_path)
linear_prior_rows <- rbind(linear_prior_raw, linear_prior_bf)
fit_provenance <- make_input_provenance(c(
  linear_raw = linear_raw_path,
  linear_bf = linear_bf_path
))
linear_fit_bf <- linear_fit_bf_compact
rm(
  statistics,
  linear_fit_bf_compact
)
invisible(gc())

message("[3/7] Loading and validating the retained matched IWP1 fits.")
current_raw <- load_exact_object(current_raw_path, "fash_fit1")
validate_matched_iwp_fit(
  current_raw,
  pair_keys = pair_keys,
  grid = grid,
  pred_step = pred_step,
  penalty = penalty,
  label = "Raw current-PC IWP1 fit"
)
current_raw_lfdr <- as.numeric(current_raw$lfdr)
current_prior_raw <- make_prior_weight_rows(
  current_raw,
  "Current FASH",
  "Raw",
  grid
)
rm(current_raw)
invisible(gc())

current_bf <- load_exact_object(current_bf_path, "fash_fit1_update")
validate_matched_iwp_fit(
  current_bf,
  pair_keys = pair_keys,
  grid = grid,
  pred_step = pred_step,
  penalty = penalty,
  label = "BF-adjusted current-PC IWP1 fit"
)
current_bf_lfdr <- as.numeric(current_bf$lfdr)
current_prior_bf <- make_prior_weight_rows(
  current_bf,
  "Current FASH",
  "BF-adjusted",
  grid
)
prior_weights <- rbind(
  current_prior_raw,
  current_prior_bf,
  linear_prior_rows
)
prior_groups <- split(
  prior_weights,
  interaction(
    prior_weights$method,
    prior_weights$adjustment,
    drop = TRUE,
    lex.order = TRUE
  )
)
prior_summary <- do.call(rbind, lapply(prior_groups, make_prior_summary_row))
rownames(prior_summary) <- NULL
prior_summary$method <- factor(
  prior_summary$method,
  levels = c("Current FASH", "FASH-linear")
)
prior_summary$adjustment <- factor(
  prior_summary$adjustment,
  levels = c("Raw", "BF-adjusted")
)
prior_summary <- prior_summary[order(
  prior_summary$method,
  prior_summary$adjustment
), , drop = FALSE]
prior_summary$method <- as.character(prior_summary$method)
prior_summary$adjustment <- as.character(prior_summary$adjustment)
rownames(prior_summary) <- NULL

message("[4/7] Building matched discovery and overlap summaries.")
pair_table <- parse_linear_pair_keys(pair_keys)
discoveries <- list(
  current_raw = make_discovery_row(
    "Current FASH", "Raw", current_raw_lfdr, pair_table, alpha
  ),
  current_bf = make_discovery_row(
    "Current FASH", "BF-adjusted", current_bf_lfdr, pair_table, alpha
  ),
  linear_raw = make_discovery_row(
    "FASH-linear", "Raw", linear_raw_lfdr, pair_table, alpha
  ),
  linear_bf = make_discovery_row(
    "FASH-linear", "BF-adjusted", linear_bf_lfdr, pair_table, alpha
  )
)
discovery_counts <- do.call(rbind, lapply(discoveries, `[[`, "row"))
rownames(discovery_counts) <- NULL
expected_current_counts <- data.frame(
  method = rep("Current FASH", 2L),
  adjustment = c("Raw", "BF-adjusted"),
  pair_count = c(43860L, 9205L),
  gene_count = c(3258L, 1177L),
  variant_count = c(42893L, 9139L),
  stringsAsFactors = FALSE
)
observed_current_counts <- discovery_counts[
  discovery_counts$method == "Current FASH",
  names(expected_current_counts),
  drop = FALSE
]
rownames(observed_current_counts) <- NULL
if (!identical(observed_current_counts, expected_current_counts)) {
  stop("The authoritative current-PC IWP1 discovery invariants changed.")
}

historical_cache <- readRDS(historical_cache_path)
expected_input_md5 <- c(
  datasets = "689d8c8e63844b1a9d7814d468533dd8",
  current_iwp1_raw = "e8aa07c1ebbd600b5dd192e5e91cb794",
  current_iwp1_bf = "af46b39a04b241cf1e6116a645dbcc3e",
  historical_profile_raw = "50ebfabca31b4b6aad93428768212fd5",
  historical_profile_bf = "54488e51ad6b737ee3c1741d9daafa35",
  historical_profile_cache = "f2fede0d8590e35c3524dd07906ce6f4"
)
observed_input_md5 <- stats::setNames(
  input_provenance$md5[match(names(expected_input_md5), input_provenance$label)],
  names(expected_input_md5)
)
unchanged_iwp_validation <- data.frame(
  check = c(
    paste0("input_md5_", names(expected_input_md5)),
    "pair_count_and_order",
    "matched_grid_predstep_penalty",
    "raw_discovery_counts",
    "bf_discovery_counts",
    "historical_bf_lfdr"
  ),
  pass = c(
    observed_input_md5 == expected_input_md5,
    length(pair_keys) == expected_pairs && !anyDuplicated(pair_keys),
    TRUE,
    identical(
      observed_current_counts[1L, , drop = FALSE],
      expected_current_counts[1L, , drop = FALSE]
    ),
    identical(
      observed_current_counts[2L, , drop = FALSE],
      expected_current_counts[2L, , drop = FALSE]
    ),
    isTRUE(all.equal(
      historical_cache$lfdr_scatter_all$current_lfdr,
      current_bf_lfdr,
      tolerance = 0
    ))
  ),
  stringsAsFactors = FALSE
)
if (any(!unchanged_iwp_validation$pass)) {
  stop("At least one unchanged-IWP validation check failed.")
}

current_bf_table <- discoveries$current_bf$table
linear_bf_table <- discoveries$linear_bf$table
venn_sets <- list(
  `Gene-variant pairs` = list(
    `Current FASH` = current_bf_table$key,
    `FASH-linear` = linear_bf_table$key
  ),
  Genes = list(
    `Current FASH` = unique(current_bf_table$gene_id),
    `FASH-linear` = unique(linear_bf_table$gene_id)
  ),
  Variants = list(
    `Current FASH` = unique(current_bf_table$variant_id),
    `FASH-linear` = unique(linear_bf_table$variant_id)
  )
)
overlap_summary <- do.call(rbind, Map(
  make_overlap_row,
  names(venn_sets),
  lapply(venn_sets, `[[`, "Current FASH"),
  lapply(venn_sets, `[[`, "FASH-linear")
))
rownames(overlap_summary) <- NULL

current_called <- logical(expected_pairs)
current_called[discoveries$current_bf$indices] <- TRUE
linear_called <- logical(expected_pairs)
linear_called[discoveries$linear_bf$indices] <- TRUE
discovery_status <- ifelse(
  current_called & linear_called,
  "Both",
  ifelse(
    current_called,
    "Current FASH only",
    ifelse(linear_called, "FASH-linear only", "Neither")
  )
)
discovery_status <- factor(
  discovery_status,
  levels = c("Neither", "Current FASH only", "FASH-linear only", "Both")
)
lfdr_scatter_all <- data.frame(
  current_lfdr = current_bf_lfdr,
  linear_lfdr = linear_bf_lfdr,
  discovery_status = discovery_status
)
top_indices <- select_current_top_pair_per_gene(pair_table, current_bf_lfdr)
lfdr_scatter_top <- data.frame(
  key = pair_table$key[top_indices],
  gene_id = pair_table$gene_id[top_indices],
  variant_id = pair_table$variant_id[top_indices],
  current_lfdr = current_bf_lfdr[top_indices],
  linear_lfdr = linear_fit_bf$lfdr[top_indices],
  discovery_status = discovery_status[top_indices],
  stringsAsFactors = FALSE
)
lfdr_summary <- data.frame(
  selection = c(
    "All tested gene-variant pairs",
    "Current-FASH top pair per gene"
  ),
  pair_count = c(expected_pairs, nrow(lfdr_scatter_top)),
  gene_count = c(
    length(unique(pair_table$gene_id)),
    nrow(lfdr_scatter_top)
  ),
  pearson_correlation = c(
    stats::cor(current_bf_lfdr, linear_fit_bf$lfdr),
    stats::cor(lfdr_scatter_top$current_lfdr, lfdr_scatter_top$linear_lfdr)
  ),
  spearman_correlation = c(
    stats::cor(current_bf_lfdr, linear_fit_bf$lfdr, method = "spearman"),
    stats::cor(
      lfdr_scatter_top$current_lfdr,
      lfdr_scatter_top$linear_lfdr,
      method = "spearman"
    )
  ),
  stringsAsFactors = FALSE
)

message("[5/7] Comparing the mixture with the historical profiled slab.")
historical_profile_raw <- readRDS(historical_profile_raw_path)
historical_profile_bf <- readRDS(historical_profile_bf_path)
if (!inherits(historical_profile_raw, "profiled_linear_fash") ||
    !inherits(historical_profile_bf, "profiled_linear_fash") ||
    !identical(historical_profile_raw$log_marginal$unit_id, pair_keys) ||
    !identical(historical_profile_bf$log_marginal$unit_id, pair_keys)) {
  stop("The historical profiled linear fits are not pair-key aligned.")
}
linear_version_comparison <- rbind(
  compare_linear_versions(
    "Raw",
    historical_profile_raw$lfdr,
    linear_raw_lfdr,
    pair_table,
    alpha
  ),
  compare_linear_versions(
    "BF-adjusted",
    historical_profile_bf$lfdr,
    linear_bf_lfdr,
    pair_table,
    alpha
  )
)
rm(historical_profile_raw, historical_profile_bf, historical_cache)
invisible(gc())

linear_only_indices <- setdiff(
  discoveries$linear_bf$indices,
  discoveries$current_bf$indices
)
if (length(linear_only_indices) == 0L) {
  stop("No BF-adjusted FASH-linear-only pair is available for the example.")
}
example_order <- order(
  linear_bf_lfdr[linear_only_indices],
  -current_bf_lfdr[linear_only_indices],
  pair_keys[linear_only_indices]
)
example_index <- linear_only_indices[example_order[1L]]
example_key <- pair_keys[example_index]
example_pair <- pair_table[example_index, , drop = FALSE]
gene_map <- readRDS(gene_map_path)
example_symbol <- gene_map$hgnc_symbol[match(
  example_pair$gene_id,
  gene_map$ensembl_gene_id
)]
if (is.na(example_symbol) || example_symbol == "") {
  example_symbol <- example_pair$gene_id
}
example_grid <- seq(0, 15, by = 0.1)
example_data <- current_bf$fash_data$data_list[[example_index]]
example_standard_error <- as.numeric(current_bf$fash_data$S[[example_index]])
example_iwp_prediction <- predict(
  current_bf,
  index = example_index,
  smooth_var = example_grid
)
example_iwp <- data.frame(
  method = "IWP1 FASH",
  time = example_iwp_prediction$x,
  posterior_mean = example_iwp_prediction$mean,
  lower = example_iwp_prediction$lower,
  upper = example_iwp_prediction$upper,
  stringsAsFactors = FALSE
)
example_linear <- extract_linear_mixture_posterior_plot_data(
  dataset = example_data,
  standard_error = example_standard_error,
  fit = linear_fit_bf,
  unit_index = example_index,
  grid = example_grid,
  sample_size = 10000L,
  seed = 20260810L + example_index
)
example_linear$method <- "FASH-linear"
example_observed <- data.frame(
  time = as.numeric(example_data$x),
  beta = as.numeric(example_data$y),
  standard_error = example_standard_error,
  stringsAsFactors = FALSE
)
example_summary <- data.frame(
  selection_rule = paste(
    "Smallest BF-adjusted FASH-linear lfdr among FASH-linear-only calls;",
    "then largest IWP1 lfdr and pair key"
  ),
  original_index = example_index,
  key = example_key,
  gene_id = example_pair$gene_id,
  gene_symbol = example_symbol,
  variant_id = example_pair$variant_id,
  current_lfdr = current_bf_lfdr[example_index],
  linear_lfdr = linear_bf_lfdr[example_index],
  stringsAsFactors = FALSE
)
rm(current_bf, linear_fit_bf)
invisible(gc())

analysis_elapsed <- proc.time()[["elapsed"]] - analysis_started
runtime_summary <- data.frame(
  task = paste(
    "Full current-PC predictive-SD mixture FASH-linear",
    "real-data ablation"
  ),
  elapsed_seconds = analysis_elapsed,
  pair_count = expected_pairs,
  predictive_sd_grid_size = length(grid),
  stringsAsFactors = FALSE
)
configuration <- list(
  analysis_id = analysis_id,
  output_id = basename(output_directory),
  alpha = alpha,
  pc_correction = "Time-specific PCs",
  current_model = "IWP1 FASH",
  comparator_model = paste(
    "Unrestricted constant plus a finite Gaussian mixture over linear",
    "one-step departures"
  ),
  linear_prior_mode = "mixture_grid",
  scale_definition = "sd_linear_departure_at_pred_step",
  time_grid = expected_time,
  psd_grid = grid,
  pred_step = pred_step,
  penalty = penalty,
  sufficient_statistic_time_parameterization = "(time - minimum) / 15",
  sufficient_statistic_time_span = statistic_time_span,
  sufficient_statistic_sd_multiplier = statistic_time_span / pred_step,
  sufficient_statistics_source = statistics_cache_source,
  sufficient_statistic_block_size = block_size,
  ridge = 1e-10,
  historical_comparator = "profiled single Gaussian slope slab",
  historical_output_id = basename(historical_output_directory)
)

analysis_cache <- list(
  configuration = configuration,
  input_provenance = input_provenance,
  fit_provenance = fit_provenance,
  unchanged_iwp_validation = unchanged_iwp_validation,
  discovery_counts = discovery_counts,
  overlap_summary = overlap_summary,
  venn_sets = venn_sets,
  lfdr_scatter_all = lfdr_scatter_all,
  lfdr_scatter_top = lfdr_scatter_top,
  lfdr_summary = lfdr_summary,
  prior_summary = prior_summary,
  prior_weights = prior_weights,
  linear_version_comparison = linear_version_comparison,
  example_summary = example_summary,
  example_observed = example_observed,
  example_iwp = example_iwp,
  example_linear = example_linear,
  runtime_summary = runtime_summary
)

message("[6/7] Saving the retained cache and CSV exports.")
atomic_save_rds(
  analysis_cache,
  file.path(output_directory, "analysis_cache.rds")
)
utils::write.csv(
  prior_weights,
  file.path(output_directory, "prior_weights.csv"),
  row.names = FALSE
)
utils::write.csv(
  discovery_counts,
  file.path(output_directory, "discovery_counts.csv"),
  row.names = FALSE
)
utils::write.csv(
  overlap_summary,
  file.path(output_directory, "overlap_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  lfdr_summary,
  file.path(output_directory, "lfdr_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  lfdr_scatter_top,
  file.path(output_directory, "top_pair_comparison.csv"),
  row.names = FALSE
)
utils::write.csv(
  prior_summary,
  file.path(output_directory, "prior_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  linear_version_comparison,
  file.path(output_directory, "linear_version_comparison.csv"),
  row.names = FALSE
)
utils::write.csv(
  unchanged_iwp_validation,
  file.path(output_directory, "unchanged_iwp_validation.csv"),
  row.names = FALSE
)
utils::write.csv(
  runtime_summary,
  file.path(output_directory, "runtime_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  input_provenance,
  file.path(output_directory, "input_provenance.csv"),
  row.names = FALSE
)
utils::write.csv(
  fit_provenance,
  file.path(output_directory, "fit_provenance.csv"),
  row.names = FALSE
)

analysis_cache_path <- file.path(output_directory, "analysis_cache.rds")
atomic_save_rds(list(
  analysis_id = analysis_id,
  status = "complete",
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  analysis_cache_md5 = unname(tools::md5sum(analysis_cache_path)),
  fit_md5 = stats::setNames(fit_provenance$md5, fit_provenance$label)
), run_status_path)

message("[7/7] Completed all matched-mixture validation checks.")
message(
  "Completed in ", format(round(analysis_elapsed, 1), nsmall = 1),
  " seconds with ", length(grid), " predictive-SD components."
)
print(discovery_counts)
print(overlap_summary)
print(linear_version_comparison)
