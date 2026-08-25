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

make_cl_dataset <- function(data, standard_error) {
  data.frame(
    time = as.numeric(data$x),
    beta = as.numeric(data$y),
    SE = as.numeric(standard_error),
    stringsAsFactors = FALSE
  )
}

parse_pair_table <- function(pair_keys) {
  table <- parse_linear_pair_keys(pair_keys)
  names(table)[names(table) == "key"] <- "pair_key"
  table
}

select_top_pairs_per_gene <- function(pair_table, lfdr, selected_indices) {
  selected <- pair_table[selected_indices, , drop = FALSE]
  selected$lfdr <- as.numeric(lfdr[selected_indices])
  selected <- selected[
    order(selected$gene_id, selected$lfdr, selected$variant_id,
          selected$pair_key, method = "radix"),
    ,
    drop = FALSE
  ]
  selected[!duplicated(selected$gene_id), , drop = FALSE]
}

make_method_set <- function(pair_keys, selected_indices) {
  selected <- pair_keys[selected_indices]
  parsed <- parse_pair_table(selected)
  list(
    pairs = selected,
    genes = unique(parsed$gene_id),
    variants = unique(parsed$variant_id)
  )
}

make_method_summary <- function(method, adjustment, method_set) {
  data.frame(
    method = method,
    adjustment = adjustment,
    pair_count = length(method_set$pairs),
    gene_count = length(method_set$genes),
    variant_count = length(method_set$variants),
    stringsAsFactors = FALSE
  )
}

make_overlap_rows <- function(first_method, first_set, second_method, second_set) {
  make_row <- function(unit, first_values, second_values) {
    intersection_count <- length(intersect(first_values, second_values))
    union_count <- length(union(first_values, second_values))
    data.frame(
      first_method = first_method,
      second_method = second_method,
      unit = unit,
      first_count = length(first_values),
      second_count = length(second_values),
      intersection_count = intersection_count,
      first_only_count = length(setdiff(first_values, second_values)),
      second_only_count = length(setdiff(second_values, first_values)),
      union_count = union_count,
      jaccard = if (union_count) intersection_count / union_count else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, list(
    make_row("Gene-variant pairs", first_set$pairs, second_set$pairs),
    make_row("Genes", first_set$genes, second_set$genes),
    make_row("Unique variants", first_set$variants, second_set$variants)
  ))
}

make_top_category_table <- function(top_pairs, strober_keys, linear_keys) {
  top_pairs$strober_found <- top_pairs$pair_key %in% strober_keys
  top_pairs$linear_found <- top_pairs$pair_key %in% linear_keys
  top_pairs$category <- ifelse(
    top_pairs$strober_found & top_pairs$linear_found,
    "Strober found / CL-fash-linear found",
    ifelse(
      top_pairs$strober_found & !top_pairs$linear_found,
      "Strober found / CL-fash-linear missed",
      ifelse(
        !top_pairs$strober_found & top_pairs$linear_found,
        "Strober missed / CL-fash-linear found",
        "Strober missed / CL-fash-linear missed"
      )
    )
  )
  category_levels <- c(
    "Strober found / CL-fash-linear found",
    "Strober found / CL-fash-linear missed",
    "Strober missed / CL-fash-linear found",
    "Strober missed / CL-fash-linear missed"
  )
  counts <- table(factor(top_pairs$category, levels = category_levels))
  data.frame(
    category = names(counts),
    pair_count = as.integer(counts),
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
helper_path <- file.path(analysis_directory, "fash_linear_real_data_helpers.R")
runner_path <- file.path(
  analysis_directory,
  "run_fash_linear_real_data_ablation_cl.R"
)
source(helper_path)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.")
}

input_directory <- file.path(workflowr_root, "output", "dynamic_eQTL_real")
cl_raw_path <- file.path(input_directory, "fash_fit1_all_CL.RData")
cl_bf_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fash_cl_variant_enrichment_comparison",
  "fash_cl_bf_adjustment.rds"
)
strober_linear_path <- file.path(
  workflowr_root,
  "data",
  "dynamic_eQTL_real",
  "strober_linear",
  "linear_dynamic_eqtls_5_pc.txt"
)
gene_map_path <- file.path(input_directory, "cache_gene_map.rds")
helper_current_fit_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fash_linear_real_data_ablation",
  "linear_fit_bf.rds"
)
fash_strober_cache_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fash_strober_enhancer_comparison",
  "analysis_cache.rds"
)
shared_path <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "shared",
  "simulation_functions.R"
)
input_paths <- c(
  cl_raw = cl_raw_path,
  cl_bf_adjustment = cl_bf_path,
  strober_linear = strober_linear_path,
  gene_map = gene_map_path,
  current_linear_bf = helper_current_fit_path,
  fash_strober_cache = fash_strober_cache_path,
  helper = helper_path,
  shared_functions = shared_path,
  runner = runner_path
)
if (any(!file.exists(input_paths))) {
  stop("At least one required CL comparison input is missing.")
}

output_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fash_linear_real_data_ablation_cl"
)
block_directory <- file.path(output_directory, "sufficient_statistic_blocks")
dir.create(block_directory, recursive = TRUE, showWarnings = FALSE)

input_provenance <- make_input_provenance(input_paths)
alpha <- 0.05
expected_pairs <- 1009173L
expected_time <- 0:15
block_size <- 25000L
sigma_beta_grid <- exp(seq(log(0.05), log(5), length.out = 25))
analysis_started <- proc.time()[["elapsed"]]

message("[1/7] Loading and validating the retained FASH-CL fit.")
cl_fit <- load_exact_object(cl_raw_path, "fash_fit1")
cl_pair_keys <- names(cl_fit$fash_data$data_list)
if (length(cl_pair_keys) != expected_pairs ||
    anyDuplicated(cl_pair_keys) ||
    length(cl_fit$fash_data$S) != expected_pairs) {
  stop("The retained FASH-CL fit failed pair-key or standard-error validation.")
}
if (!all(vapply(
  cl_fit$fash_data$data_list,
  function(data) is.data.frame(data) && nrow(data) == length(expected_time),
  logical(1)
))) {
  stop("At least one retained FASH-CL trajectory has the wrong time dimension.")
}
if (!all(vapply(
  cl_fit$fash_data$S,
  function(standard_error) {
    length(standard_error) == length(expected_time) &&
      all(is.finite(standard_error)) && all(standard_error > 0)
  },
  logical(1)
))) {
  stop("At least one retained FASH-CL standard-error vector is invalid.")
}

message("[2/7] Extracting or resuming CL sufficient statistics.")
block_starts <- seq.int(1L, expected_pairs, by = block_size)
block_paths <- character(length(block_starts))
cl_raw_md5 <- input_provenance$md5[input_provenance$label == "cl_raw"]
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
    existing <- readRDS(block_path)
    valid_existing <- is.list(existing) &&
      identical(existing$cl_raw_md5, cl_raw_md5) &&
      identical(existing$start, start) &&
      identical(existing$end, end) &&
      identical(existing$expected_time, expected_time) &&
      identical(existing$scale_time, TRUE) &&
      is.data.frame(existing$statistics) &&
      nrow(existing$statistics) == end - start + 1L &&
      identical(existing$statistics$unit_id, cl_pair_keys[start:end])
  }
  if (!valid_existing) {
    cl_data_block <- Map(
      make_cl_dataset,
      cl_fit$fash_data$data_list[start:end],
      cl_fit$fash_data$S[start:end]
    )
    statistics <- compute_linear_sufficient_statistics_block(
      cl_data_block,
      expected_time = expected_time,
      scale_time = TRUE
    )
    existing <- list(
      cl_raw_md5 = cl_raw_md5,
      start = start,
      end = end,
      expected_time = expected_time,
      scale_time = TRUE,
      statistics = statistics
    )
    atomic_save_rds(existing, block_path)
    rm(cl_data_block, statistics)
  }
  rm(existing)
  if (block_index %% 4L == 0L || end == expected_pairs) {
    invisible(gc())
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
if (nrow(statistics) != expected_pairs ||
    !identical(statistics$unit_id, cl_pair_keys)) {
  stop("The assembled CL sufficient statistics are not aligned.")
}
sufficient_statistics_path <- file.path(output_directory, "sufficient_statistics.rds")
atomic_save_rds(
  list(
    analysis_id = "revision_internal_fash_linear_real_data_ablation_cl",
    cl_raw_md5 = cl_raw_md5,
    expected_time = expected_time,
    scale_time = TRUE,
    ridge = 1e-10,
    statistics = statistics
  ),
  sufficient_statistics_path
)
rm(cl_fit)
invisible(gc())

message("[3/7] Profiling CL FASH-linear and applying the BF update.")
linear_fit_raw <- fit_profiled_linear_fash_from_stats(
  statistics,
  sigma_beta_grid = sigma_beta_grid,
  n_time = length(expected_time)
)
linear_fit_bf <- bf_update_profiled_linear_fash(linear_fit_raw)
atomic_save_rds(
  linear_fit_raw,
  file.path(output_directory, "linear_fit_raw.rds")
)
atomic_save_rds(
  linear_fit_bf,
  file.path(output_directory, "linear_fit_bf.rds")
)

message("[4/7] Loading CL IWP1 and Strober-linear discovery inputs.")
cl_bf_adjustment <- readRDS(cl_bf_path)
if (cl_bf_adjustment$pair_count != expected_pairs ||
    length(cl_bf_adjustment$lfdr) != expected_pairs ||
    any(!is.finite(cl_bf_adjustment$lfdr)) ||
    any(cl_bf_adjustment$lfdr < 0 | cl_bf_adjustment$lfdr > 1)) {
  stop("The retained FASH-CL BF-adjustment cache failed validation.")
}
strober_linear <- data.table::fread(
  strober_linear_path,
  data.table = FALSE,
  showProgress = FALSE
)
required_strober_columns <- c("rs_id", "ensamble_id", "pvalue", "eFDR")
if (!identical(names(strober_linear), required_strober_columns) ||
    nrow(strober_linear) != expected_pairs ||
    any(!is.finite(strober_linear$eFDR)) ||
    any(strober_linear$eFDR < 0 | strober_linear$eFDR > 1)) {
  stop("The Strober-linear table failed schema or value validation.")
}
strober_linear$pair_key <- paste(
  strober_linear$ensamble_id,
  strober_linear$rs_id,
  sep = "_"
)
if (anyDuplicated(strober_linear$pair_key) ||
    !setequal(strober_linear$pair_key, cl_pair_keys)) {
  stop("The Strober-linear table does not match the CL pair universe.")
}
strober_keys <- strober_linear$pair_key[strober_linear$eFDR <= alpha]

pair_table <- parse_pair_table(cl_pair_keys)
cl_iwp1_calls <- select_cumulative_lfdr_calls_linear(
  cl_bf_adjustment$lfdr,
  alpha = alpha
)
cl_linear_calls <- select_cumulative_lfdr_calls_linear(
  linear_fit_bf$lfdr,
  alpha = alpha
)
strober_calls <- match(strober_keys, cl_pair_keys)
method_sets <- list(
  `FASH-CL IWP1` = make_method_set(cl_pair_keys, cl_iwp1_calls),
  `FASH-CL-linear` = make_method_set(cl_pair_keys, cl_linear_calls),
  `Strober linear` = make_method_set(cl_pair_keys, strober_calls)
)
if (length(method_sets[["FASH-CL IWP1"]]$pairs) != 5393L ||
    length(method_sets[["FASH-CL IWP1"]]$genes) != 686L ||
    length(method_sets[["Strober linear"]]$pairs) != 5404L) {
  stop("The retained CL or Strober discovery invariants changed.")
}

discovery_summary <- do.call(rbind, list(
  make_method_summary("FASH-CL IWP1", "BF-adjusted", method_sets[["FASH-CL IWP1"]]),
  make_method_summary("FASH-CL-linear", "BF-adjusted", method_sets[["FASH-CL-linear"]]),
  make_method_summary("Strober linear", "eFDR", method_sets[["Strober linear"]])
))
row.names(discovery_summary) <- NULL

message("[5/7] Building CL overlap and top-pair summaries.")
overlap_pairs <- list(
  c("FASH-CL IWP1", "Strober linear"),
  c("FASH-CL-linear", "Strober linear"),
  c("FASH-CL IWP1", "FASH-CL-linear")
)
overlap_summary <- do.call(rbind, lapply(overlap_pairs, function(pair) {
  make_overlap_rows(
    pair[[1L]],
    method_sets[[pair[[1L]]]],
    pair[[2L]],
    method_sets[[pair[[2L]]]]
  )
}))
row.names(overlap_summary) <- NULL

top_pairs <- select_top_pairs_per_gene(
  pair_table,
  cl_bf_adjustment$lfdr,
  cl_iwp1_calls
)
top_pairs$strober_found <- top_pairs$pair_key %in% strober_keys
top_pairs$linear_found <- top_pairs$pair_key %in%
  method_sets[["FASH-CL-linear"]]$pairs
top_pairs$category <- ifelse(
  top_pairs$strober_found & top_pairs$linear_found,
  "Strober found / CL-fash-linear found",
  ifelse(
    top_pairs$strober_found & !top_pairs$linear_found,
    "Strober found / CL-fash-linear missed",
    ifelse(
      !top_pairs$strober_found & top_pairs$linear_found,
      "Strober missed / CL-fash-linear found",
      "Strober missed / CL-fash-linear missed"
    )
  )
)
top_category_counts <- make_top_category_table(
  top_pairs,
  strober_keys,
  method_sets[["FASH-CL-linear"]]$pairs
)
gene_map <- readRDS(gene_map_path)
if (!all(c("ensembl_gene_id", "hgnc_symbol") %in% names(gene_map))) {
  stop("The gene map is missing required columns.")
}
top_pairs$gene_symbol <- gene_map$hgnc_symbol[
  match(top_pairs$gene_id, gene_map$ensembl_gene_id)
]
missing_symbol <- is.na(top_pairs$gene_symbol) | top_pairs$gene_symbol == ""
top_pairs$gene_symbol[missing_symbol] <- top_pairs$gene_id[missing_symbol]
category_examples <- do.call(rbind, lapply(
  unique(top_category_counts$category),
  function(category) {
    selected <- top_pairs[top_pairs$category == category, , drop = FALSE]
    selected <- selected[
      order(selected$lfdr, selected$pair_key, method = "radix"),
      ,
      drop = FALSE
    ]
    selected <- selected[seq_len(min(5L, nrow(selected))), , drop = FALSE]
    data.frame(
      category = category,
      gene_id = selected$gene_id,
      gene_symbol = selected$gene_symbol,
      variant_id = selected$variant_id,
      pair_key = selected$pair_key,
      cl_iwp1_lfdr = selected$lfdr,
      stringsAsFactors = FALSE
    )
  }
))
row.names(category_examples) <- NULL

message("[6/7] Adding current-PC FASH-linear baseline closeness.")
current_cache <- readRDS(fash_strober_cache_path)
current_export <- current_cache$discovery_export
current_export$pair_key <- paste(
  current_export$gene_id,
  current_export$variant_id,
  sep = "_"
)
current_fash_set <- list(
  pairs = current_export$pair_key[
    current_export$discovery_set == "current_all"
  ],
  genes = unique(current_export$gene_id[
    current_export$discovery_set == "current_all"
  ]),
  variants = unique(current_export$variant_id[
    current_export$discovery_set == "current_all"
  ])
)
current_linear_fit <- readRDS(helper_current_fit_path)
current_linear_calls <- select_cumulative_lfdr_calls_linear(
  current_linear_fit$lfdr,
  alpha = alpha
)
current_linear_keys <- current_linear_fit$log_marginal$unit_id
current_linear_set <- make_method_set(current_linear_keys, current_linear_calls)
current_strober_set <- list(
  pairs = strober_keys,
  genes = method_sets[["Strober linear"]]$genes,
  variants = method_sets[["Strober linear"]]$variants
)
baseline_overlap <- do.call(rbind, list(
  make_overlap_rows(
    "Current FASH IWP1",
    current_fash_set,
    "Strober linear",
    current_strober_set
  ),
  make_overlap_rows(
    "Current FASH-linear",
    current_linear_set,
    "Strober linear",
    current_strober_set
  )
))
baseline_overlap$comparison_scope <- "Current-PC baseline"
overlap_summary$comparison_scope <- "FASH-CL experiment"
overlap_summary <- rbind(overlap_summary, baseline_overlap)
row.names(overlap_summary) <- NULL

make_rank_correlation_row <- function(method, pair_keys, lfdr) {
  strober_index <- match(pair_keys, strober_linear$pair_key)
  data.frame(
    method = method,
    spearman_lfdr_pvalue = cor(
      as.numeric(lfdr),
      strober_linear$pvalue[strober_index],
      method = "spearman"
    ),
    spearman_lfdr_eFDR = cor(
      as.numeric(lfdr),
      strober_linear$eFDR[strober_index],
      method = "spearman"
    ),
    stringsAsFactors = FALSE
  )
}
rank_correlation_summary <- rbind(
  make_rank_correlation_row(
    "FASH-CL-linear",
    cl_pair_keys,
    linear_fit_bf$lfdr
  ),
  make_rank_correlation_row(
    "Current FASH-linear",
    current_linear_keys,
    current_linear_fit$lfdr
  )
)

prior_summary <- data.frame(
  method = c("FASH-CL IWP1", "FASH-CL-linear", "Current FASH-linear"),
  adjustment = c("BF-adjusted", "BF-adjusted", "BF-adjusted"),
  null_weight = c(
    cl_bf_adjustment$null_weight,
    linear_fit_bf$prior_weights$prior_weight[
      linear_fit_bf$prior_weights$component == "constant"
    ],
    current_linear_fit$prior_weights$prior_weight[
      current_linear_fit$prior_weights$component == "constant"
    ]
  ),
  selected_sigma_beta = c(
    NA_real_,
    linear_fit_bf$sigma_beta,
    current_linear_fit$sigma_beta
  ),
  stringsAsFactors = FALSE
)

runtime_summary <- data.frame(
  stage = "CL sufficient-statistic extraction and profiled FASH-linear fit",
  elapsed_seconds = proc.time()[["elapsed"]] - analysis_started,
  stringsAsFactors = FALSE
)
configuration <- list(
  analysis_id = "revision_internal_fash_linear_real_data_ablation_cl",
  alpha = alpha,
  expected_pairs = expected_pairs,
  expected_time = expected_time,
  block_size = block_size,
  sigma_beta_grid = sigma_beta_grid,
  scale_time = TRUE,
  pc_correction = "Cell-line-specific PCs repeated across time",
  strober_pc_correction = "Cell-line-specific PCs in the published 5-PC Strober result",
  pc_correction_matched_to_strober = TRUE,
  cl_iwp1_bf_cache_null_weight = cl_bf_adjustment$null_weight,
  cl_linear_selected_sigma_beta = linear_fit_bf$sigma_beta
)

analysis_cache <- list(
  configuration = configuration,
  input_provenance = input_provenance,
  discovery_summary = discovery_summary,
  overlap_summary = overlap_summary,
  top_pair_categories = top_category_counts,
  top_pair_table = top_pairs,
  category_examples = category_examples,
  rank_correlation_summary = rank_correlation_summary,
  prior_summary = prior_summary,
  linear_sigma_profile = linear_fit_bf$sigma_profile,
  runtime_summary = runtime_summary
)

message("[7/7] Saving CL comparison summaries.")
atomic_save_rds(analysis_cache, file.path(output_directory, "analysis_cache.rds"))
write.csv(
  discovery_summary,
  file.path(output_directory, "discovery_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  overlap_summary,
  file.path(output_directory, "overlap_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  top_category_counts,
  file.path(output_directory, "top_pair_categories.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  category_examples,
  file.path(output_directory, "category_examples.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  rank_correlation_summary,
  file.path(output_directory, "rank_correlation_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  prior_summary,
  file.path(output_directory, "prior_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  input_provenance,
  file.path(output_directory, "input_provenance.csv"),
  row.names = FALSE,
  quote = TRUE
)
write.csv(
  runtime_summary,
  file.path(output_directory, "runtime_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)

cat("FASH-linear CL comparison completed.\n")
cat("Selected CL slope SD:", format(linear_fit_bf$sigma_beta, digits = 12), "\n")
cat("CL FASH-linear BF-adjusted pairs:", length(method_sets[["FASH-CL-linear"]]$pairs), "\n")
