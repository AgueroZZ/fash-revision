#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

find_workflowr_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "_workflowr.yml"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) stop("Could not locate the workflowr project root.")
    current <- parent
  }
}

load_exact_object <- function(path, object_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  if (!identical(loaded, object_name)) {
    stop("Expected only object ", object_name, " in ", path)
  }
  environment[[object_name]]
}

make_provenance <- function(paths) {
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

make_call_vector <- function(lfdr, alpha = 0.05) {
  result <- logical(length(lfdr))
  result[select_cumulative_lfdr_calls_linear(lfdr, alpha)] <- TRUE
  result
}

read_strober_table <- function(path, pair_keys, label) {
  table <- data.table::fread(
    path,
    select = c("rs_id", "ensamble_id", "pvalue", "eFDR"),
    data.table = FALSE,
    showProgress = FALSE
  )
  keys <- paste(table$ensamble_id, table$rs_id, sep = "_")
  if (nrow(table) != length(pair_keys) || anyDuplicated(keys)) {
    stop(label, " Strober results failed row-count or uniqueness validation.")
  }
  alignment <- match(pair_keys, keys)
  if (anyNA(alignment) || any(!is.finite(table$pvalue[alignment])) ||
      any(!is.finite(table$eFDR[alignment])) ||
      any(table$pvalue[alignment] < 0 | table$pvalue[alignment] > 1) ||
      any(table$eFDR[alignment] < 0 | table$eFDR[alignment] > 1)) {
    stop(label, " Strober results are incomplete or invalid after key matching.")
  }
  data.frame(
    pvalue = as.numeric(table$pvalue[alignment]),
    efdr = as.numeric(table$eFDR[alignment]),
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
app_directory <- file.path(workflowr_root, "apps", "fash_trajectory_explorer")
analysis_directory <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "fash_linear_real_data_ablation"
)
shared_path <- file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
)
source(shared_path)
source(file.path(analysis_directory, "fash_linear_real_data_helpers.R"))
source(file.path(app_directory, "explorer_helpers.R"))
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required to build the explorer cache.")
}

dynamic_directory <- file.path(workflowr_root, "output", "dynamic_eQTL_real")
linear_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "fash_linear_real_data_ablation_mixture_predstep1_penalty10"
)
paths <- c(
  current_raw = file.path(dynamic_directory, "fash_fit1_all.RData"),
  current_bf = file.path(dynamic_directory, "fash_fit1_update.RData"),
  linear_raw = file.path(linear_directory, "linear_fit_raw.rds"),
  linear_bf = file.path(linear_directory, "linear_fit_bf.rds"),
  linear_run_status = file.path(linear_directory, "run_status.rds"),
  gene_map = file.path(dynamic_directory, "cache_gene_map.rds"),
  strober_linear = file.path(
    workflowr_root,
    "data", "dynamic_eQTL_real", "strober_linear", "linear_dynamic_eqtls_5_pc.txt"
  ),
  strober_nonlinear = file.path(
    workflowr_root,
    "data", "dynamic_eQTL_real", "strober_nonlinear", "non_linear_dynamic_eqtls_5_pc.txt"
  ),
  shared_functions = shared_path,
  linear_helpers = file.path(analysis_directory, "fash_linear_real_data_helpers.R"),
  explorer_helpers = file.path(app_directory, "explorer_helpers.R"),
  builder = file.path(app_directory, "build_explorer_cache.R")
)
if (any(!file.exists(paths))) stop("At least one required explorer input is missing.")

linear_run_status <- readRDS(paths[["linear_run_status"]])
current_linear_md5 <- unname(tools::md5sum(paths[c("linear_raw", "linear_bf")]))
if (!identical(
      linear_run_status$analysis_id,
      paste0(
        "revision_internal_fash_linear_real_data_ablation_",
        "mixture_predstep1_penalty10"
      )
    ) ||
    !identical(linear_run_status$status, "complete") ||
    !identical(
      current_linear_md5,
      unname(linear_run_status$fit_md5[c("linear_raw", "linear_bf")])
    )) {
  stop("The current-PC mixture run is incomplete or its fits changed.")
}

alpha <- 0.05
expected_pairs <- 1009173L
started <- proc.time()[["elapsed"]]

message("[1/5] Loading raw and BF-adjusted IWP1 lfdr values.")
current_raw_fit <- load_exact_object(paths[["current_raw"]], "fash_fit1")
pair_keys <- names(current_raw_fit$fash_data$data_list)
current_raw_lfdr <- as.numeric(current_raw_fit$lfdr)
rm(current_raw_fit)
invisible(gc())

current_bf_fit <- load_exact_object(paths[["current_bf"]], "fash_fit1_update")
if (!identical(names(current_bf_fit$fash_data$data_list), pair_keys)) {
  stop("Raw and BF-adjusted IWP1 fits are not pair-key aligned.")
}
current_bf_lfdr <- as.numeric(current_bf_fit$lfdr)
rm(current_bf_fit)
invisible(gc())

if (length(pair_keys) != expected_pairs || anyDuplicated(pair_keys) ||
    length(current_raw_lfdr) != expected_pairs ||
    length(current_bf_lfdr) != expected_pairs) {
  stop("The IWP1 fits failed pair-count, key, or lfdr validation.")
}

message("[2/5] Loading predictive-SD mixture FASH-linear lfdr values and calls.")
linear_raw_fit <- readRDS(paths[["linear_raw"]])
linear_bf_fit <- readRDS(paths[["linear_bf"]])
validate_compact_linear_mixture_fash(
  linear_raw_fit,
  expected_grid = default_revision_grid(),
  expected_pred_step = 1,
  expected_penalty = 10L
)
validate_compact_linear_mixture_fash(
  linear_bf_fit,
  expected_grid = default_revision_grid(),
  expected_pred_step = 1,
  expected_penalty = 10L
)
linear_raw_lfdr <- as.numeric(linear_raw_fit$lfdr)
linear_bf_lfdr <- as.numeric(linear_bf_fit$lfdr)
if (length(linear_raw_lfdr) != expected_pairs ||
    length(linear_bf_lfdr) != expected_pairs ||
    isTRUE(linear_raw_fit$bf_adjusted) ||
    !isTRUE(linear_bf_fit$bf_adjusted) ||
    !identical(linear_raw_fit$unit_ids, pair_keys) ||
    !identical(linear_bf_fit$unit_ids, pair_keys)) {
  stop("The mixture FASH-linear fits are not aligned with the IWP1 fits.")
}
rm(linear_raw_fit, linear_bf_fit)
invisible(gc())

current_raw_called <- make_call_vector(current_raw_lfdr, alpha)
current_bf_called <- make_call_vector(current_bf_lfdr, alpha)
linear_raw_called <- make_call_vector(linear_raw_lfdr, alpha)
linear_bf_called <- make_call_vector(linear_bf_lfdr, alpha)

message("[3/5] Joining gene symbols and both Strober analyses by pair key.")
pair_table <- parse_linear_pair_keys(pair_keys)
gene_map <- readRDS(paths[["gene_map"]])
gene_symbol <- gene_map$hgnc_symbol[match(
  pair_table$gene_id,
  gene_map$ensembl_gene_id
)]
missing_symbol <- is.na(gene_symbol) | gene_symbol == ""
gene_symbol[missing_symbol] <- pair_table$gene_id[missing_symbol]
strober_linear <- read_strober_table(paths[["strober_linear"]], pair_keys, "Linear")
strober_nonlinear <- read_strober_table(paths[["strober_nonlinear"]], pair_keys, "Nonlinear")

status <- ifelse(
  current_bf_called & linear_bf_called,
  "Both",
  ifelse(
    current_bf_called,
    "IWP1 only",
    ifelse(linear_bf_called, "FASH-linear only", "Neither")
  )
)
index <- data.frame(
  original_index = seq_len(expected_pairs),
  key = pair_table$key,
  gene_symbol = gene_symbol,
  gene_id = pair_table$gene_id,
  variant_id = pair_table$variant_id,
  discovery_status_bf = status,
  iwp1_lfdr_raw = current_raw_lfdr,
  iwp1_lfdr_bf = current_bf_lfdr,
  linear_lfdr_raw = linear_raw_lfdr,
  linear_lfdr_bf = linear_bf_lfdr,
  iwp1_called_raw = current_raw_called,
  iwp1_called_bf = current_bf_called,
  linear_called_raw = linear_raw_called,
  linear_called_bf = linear_bf_called,
  strober_linear_pvalue = strober_linear$pvalue,
  strober_linear_efdr = strober_linear$efdr,
  strober_nonlinear_pvalue = strober_nonlinear$pvalue,
  strober_nonlinear_efdr = strober_nonlinear$efdr,
  stringsAsFactors = FALSE
)

message("[4/5] Selecting deterministic preset examples.")
linear_only <- which(linear_bf_called & !current_bf_called)
current_only <- which(current_bf_called & !linear_bf_called)
both <- which(current_bf_called & linear_bf_called)
linear_examples <- select_distinct_gene_examples(
  linear_only, pair_table, linear_bf_lfdr, current_bf_lfdr, n = 5L
)
current_examples <- select_distinct_gene_examples(
  current_only, pair_table, current_bf_lfdr, linear_bf_lfdr, n = 5L
)
shared_examples <- select_distinct_gene_examples(
  both, pair_table, pmax(current_bf_lfdr, linear_bf_lfdr),
  pmin(current_bf_lfdr, linear_bf_lfdr), n = 5L,
  secondary_decreasing = FALSE
)
discordant_examples <- select_distinct_gene_examples(
  seq_len(expected_pairs), pair_table,
  abs(current_bf_lfdr - linear_bf_lfdr),
  pmin(current_bf_lfdr, linear_bf_lfdr), n = 5L,
  primary_decreasing = TRUE, secondary_decreasing = FALSE
)
preset_groups <- list(
  `Strong FASH-linear-only` = linear_examples,
  `Strong IWP1-only` = current_examples,
  `Strong shared calls` = shared_examples,
  `Largest lfdr disagreements` = discordant_examples
)
presets <- do.call(rbind, lapply(names(preset_groups), function(group) {
  indices <- preset_groups[[group]]
  data.frame(
    group = group,
    rank = seq_along(indices),
    original_index = indices,
    label = paste0(
      group, " #", seq_along(indices), ": ",
      index$gene_symbol[indices], " / ", index$variant_id[indices]
    ),
    stringsAsFactors = FALSE
  )
}))
rownames(presets) <- NULL

if (sum(current_raw_called) != 43860L ||
    sum(current_bf_called) != 9205L ||
    length(linear_examples) != 5L ||
    length(current_examples) != 5L ||
    length(shared_examples) != 5L ||
    length(discordant_examples) != 5L ||
    nrow(presets) != 20L ||
    nrow(index) != expected_pairs || anyDuplicated(index$key) ||
    any(!is.finite(index$iwp1_lfdr_raw)) ||
    any(!is.finite(index$iwp1_lfdr_bf)) ||
    any(!is.finite(index$linear_lfdr_raw)) ||
    any(!is.finite(index$linear_lfdr_bf))) {
  stop("The assembled explorer index failed scientific invariant checks.")
}

message("[5/5] Saving the compact explorer index and provenance.")
output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "fash_trajectory_explorer_mixture_predstep1_penalty10"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
cache <- list(
  configuration = list(
    analysis_id = "fash_trajectory_explorer_mixture_predstep1_penalty10",
    pair_count = expected_pairs,
    alpha = alpha,
    pc_correction = "Time-specific PCs",
    iwp_model = "BF-adjusted IWP1 FASH",
    linear_model = "BF-adjusted predictive-SD mixture FASH-linear",
    linear_output_id = basename(linear_directory),
    linear_prior_mode = "mixture_grid",
    linear_grid = default_revision_grid(),
    pred_step = 1,
    penalty = 10L,
    strober_tests = c("linear", "nonlinear")
  ),
  provenance = make_provenance(paths),
  index = index,
  presets = presets,
  runtime_seconds = proc.time()[["elapsed"]] - started
)
temporary_path <- file.path(output_directory, "explorer_index.rds.tmp")
saveRDS(cache, temporary_path)
output_path <- file.path(output_directory, "explorer_index.rds")
if (!file.rename(temporary_path, output_path)) {
  stop("Could not atomically replace the explorer index cache.")
}

message(
  "Explorer cache saved: ", output_path, " (",
  format(round(cache$runtime_seconds, 1), nsmall = 1), " seconds)."
)
print(table(index$discovery_status_bf))
print(presets)
