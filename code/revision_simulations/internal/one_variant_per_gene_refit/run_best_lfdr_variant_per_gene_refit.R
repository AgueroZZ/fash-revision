#!/usr/bin/env Rscript

# Refit FASH(1) after keeping the minimum full-data BF-adjusted lfdr pair
# for each gene. This is an outcome-informed sensitivity analysis.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(args, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(args[equals_hit[1]], nchar(equals_prefix) + 1L))
  }
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1] == length(args)) {
    return(default)
  }
  args[hit[1] + 1L]
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

file_metadata <- function(path) {
  info <- file.info(path)
  if (nrow(info) != 1L || is.na(info$size)) {
    stop("Could not read source-file metadata: ", path)
  }
  list(
    path = normalizePath(path, mustWork = TRUE),
    size_bytes = unname(info$size),
    modification_time = format(info$mtime, tz = "UTC", usetz = TRUE)
  )
}

load_exact_object <- function(path, expected_name) {
  object_environment <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = object_environment)
  if (!identical(loaded_names, expected_name)) {
    stop(path, " must contain only ", expected_name, ".")
  }
  object_environment[[expected_name]]
}

validate_full_fit <- function(fit, name) {
  required_fields <- c(
    "prior_weights", "posterior_weights", "psd_grid", "lfdr",
    "settings", "fash_data", "L_matrix"
  )
  n_units <- length(fit$fash_data$data_list)
  pair_keys <- names(fit$fash_data$data_list)
  if (!all(required_fields %in% names(fit)) ||
      is.null(fit$fash_data$data_list) || is.null(fit$fash_data$S) ||
      n_units < 2L || length(fit$fash_data$S) != n_units ||
      length(fit$lfdr) != n_units || length(pair_keys) != n_units ||
      any(!nzchar(pair_keys)) || anyDuplicated(pair_keys) ||
      any(!is.finite(fit$lfdr)) || any(fit$lfdr < 0 | fit$lfdr > 1) ||
      length(fit$psd_grid) < 2L || any(!is.finite(fit$psd_grid)) ||
      !any(fit$psd_grid == 0)) {
    stop(name, " does not contain a valid full FASH fit.")
  }
  invisible(TRUE)
}

validate_thinned_fit <- function(fit, pair_keys, name) {
  if (length(fit$lfdr) != length(pair_keys) ||
      any(!is.finite(fit$lfdr)) || any(fit$lfdr < 0 | fit$lfdr > 1) ||
      nrow(fit$posterior_weights) != length(pair_keys) ||
      nrow(fit$L_matrix) != length(pair_keys)) {
    stop(name, " has invalid posterior dimensions or lfdr values.")
  }
  names(fit$lfdr) <- pair_keys
  fit
}

capture_fit_warnings <- function(expression) {
  warning_messages <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(warning) {
      warning_messages <<- c(warning_messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warning_messages))
}

safe_ratio <- function(numerator, denominator) {
  if (denominator == 0L) {
    return(NA_real_)
  }
  numerator / denominator
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "r4_correlated_errors",
  "real_data_correlation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "covariance_estimation",
  "zero_intercept_correlation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "one_variant_per_gene_refit",
  "one_variant_per_gene_refit_helpers.R"
))

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required.")
}

alpha <- as.numeric(get_arg("--alpha", "0.05"))
output_id <- get_arg(
  "--output-id",
  "one_variant_per_gene_best_lfdr"
)
if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
    !nzchar(output_id) || grepl("/", output_id, fixed = TRUE)) {
  stop("Invalid alpha or output_id.")
}

raw_fit_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "fash_fit1_all.RData"
)
bf_fit_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "fash_fit1_update.RData"
)
output_parent <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal"
)
final_directory <- file.path(output_parent, output_id)
staging_directory <- file.path(
  output_parent,
  paste0(".", output_id, "_staging_", Sys.getpid())
)
if (!file.exists(raw_fit_path) || !file.exists(bf_fit_path) ||
    !dir.exists(output_parent)) {
  stop("The full fits or internal output directory are missing.")
}
if (file.exists(final_directory) || file.exists(staging_directory)) {
  stop("Refusing to overwrite an existing output or staging directory.")
}

source_files <- list(
  raw_fit = file_metadata(raw_fit_path),
  bf_adjusted_fit = file_metadata(bf_fit_path)
)

message("Loading the full BF-adjusted fit for deterministic selection.")
full_bf <- load_exact_object(bf_fit_path, "fash_fit1_update")
validate_full_fit(full_bf, "fash_fit1_update")
pair_keys <- names(full_bf$fash_data$data_list)
gene_index <- make_gene_index(pair_keys)
n_full_units <- length(pair_keys)
n_genes <- length(gene_index)
selection <- select_minimum_lfdr_variant_per_gene(
  pair_keys = pair_keys,
  lfdr = full_bf$lfdr,
  gene_index = gene_index
)
if (nrow(selection) != n_genes || anyDuplicated(selection$gene_id) ||
    anyDuplicated(selection$pair_key)) {
  stop("The minimum-lfdr selection is not exactly one pair per gene.")
}

full_bf_lfdr <- full_bf$lfdr
full_bf_prior_weights <- full_bf$prior_weights
full_psd_grid <- full_bf$psd_grid
full_call_indices <- cumulative_fdr_calls(full_bf_lfdr, alpha = alpha)
full_discovered_genes <- sort(unique(parse_gene_ids(
  pair_keys[full_call_indices]
)))
selection$selected_pair_full_data_fdr_call <-
  selection$fash_index %in% full_call_indices
if (!setequal(
  selection$gene_id[selection$selected_pair_full_data_fdr_call],
  full_discovered_genes
)) {
  stop("The selected minimum-lfdr pairs do not represent the full gene set.")
}
rm(full_bf, gene_index)
gc(verbose = FALSE)

message("Loading the full raw fit and refitting the selected likelihood rows.")
full_raw <- load_exact_object(raw_fit_path, "fash_fit1")
validate_full_fit(full_raw, "fash_fit1")
if (!identical(names(full_raw$fash_data$data_list), pair_keys) ||
    !isTRUE(all.equal(full_raw$psd_grid, full_psd_grid)) ||
    !identical(full_raw$settings$penalty, 10)) {
  stop("The raw and BF-adjusted fits are not aligned or penalty is not 10.")
}

raw_start <- proc.time()[["elapsed"]]
raw_capture <- capture_fit_warnings(
  refit_fash_from_cached_likelihood(full_raw, selection$fash_index)
)
raw_elapsed <- proc.time()[["elapsed"]] - raw_start
thinned_raw <- validate_thinned_fit(
  raw_capture$value,
  selection$pair_key,
  "thinned_raw"
)

message("Applying the BF prior update to the best-lfdr thinned refit.")
bf_start <- proc.time()[["elapsed"]]
bf_capture <- capture_fit_warnings(
  fashr::BF_update(thinned_raw, plot = FALSE)
)
bf_elapsed <- proc.time()[["elapsed"]] - bf_start
thinned_bf <- validate_thinned_fit(
  bf_capture$value,
  selection$pair_key,
  "thinned_bf"
)

raw_prior <- compare_prior_weights(
  full_raw$prior_weights,
  thinned_raw$prior_weights,
  fit_stage = "Raw"
)
raw_lfdr <- compare_paired_lfdr(
  full_raw$lfdr[selection$fash_index],
  thinned_raw$lfdr,
  selection$pair_key,
  fit_stage = "Raw",
  alpha = alpha
)
full_raw_settings <- full_raw$settings
full_raw_prior_weights <- full_raw$prior_weights
rm(full_raw)
gc(verbose = FALSE)

bf_prior <- compare_prior_weights(
  full_bf_prior_weights,
  thinned_bf$prior_weights,
  fit_stage = "BF-adjusted"
)
bf_lfdr <- compare_paired_lfdr(
  full_bf_lfdr[selection$fash_index],
  thinned_bf$lfdr,
  selection$pair_key,
  fit_stage = "BF-adjusted",
  alpha = alpha
)

prior_table <- rbind(raw_prior$table, bf_prior$table)
lfdr_table <- rbind(raw_lfdr$table, bf_lfdr$table)
comparison_summary <- merge(
  rbind(raw_prior$summary, bf_prior$summary),
  rbind(raw_lfdr$summary, bf_lfdr$summary),
  by = "fit_stage",
  all = TRUE,
  sort = FALSE
)
comparison_summary$selection_rule <- "minimum_full_bf_lfdr"
comparison_summary$n_full_units <- n_full_units
comparison_summary$n_genes <- n_genes
comparison_summary$raw_eb_refit_elapsed_seconds <- unname(raw_elapsed)
comparison_summary$bf_update_elapsed_seconds <- unname(bf_elapsed)
comparison_summary <- comparison_summary[, c(
  "selection_rule", "fit_stage", "n_full_units", "n_genes", "n_units",
  "full_pi0", "thinned_pi0", "pi0_difference",
  "prior_total_variation", "pearson_lfdr", "spearman_lfdr",
  "mean_absolute_lfdr_difference", "median_absolute_lfdr_difference",
  "rmse_lfdr", "full_mean_lfdr", "thinned_mean_lfdr",
  "full_fdr_calls", "thinned_fdr_calls", "fdr_call_intersection",
  "fdr_call_union", "fdr_call_jaccard", "alpha",
  "raw_eb_refit_elapsed_seconds", "bf_update_elapsed_seconds"
)]

thinned_call_indices <- cumulative_fdr_calls(thinned_bf$lfdr, alpha = alpha)
thinned_discovered_genes <- sort(selection$gene_id[thinned_call_indices])
overlap_genes <- intersect(full_discovered_genes, thinned_discovered_genes)
missed_genes <- setdiff(full_discovered_genes, thinned_discovered_genes)
additional_genes <- setdiff(thinned_discovered_genes, full_discovered_genes)
gene_union <- union(full_discovered_genes, thinned_discovered_genes)

gene_discovery_comparison <- data.frame(
  gene_id = selection$gene_id,
  selected_pair_key = selection$pair_key,
  selection_lfdr = selection$selection_lfdr,
  selected_pair_full_data_fdr_call =
    selection$selected_pair_full_data_fdr_call,
  full_data_discovered_gene = selection$gene_id %in% full_discovered_genes,
  thinned_refit_discovered_gene =
    seq_len(nrow(selection)) %in% thinned_call_indices,
  stringsAsFactors = FALSE
)
gene_discovery_comparison$discovery_category <- with(
  gene_discovery_comparison,
  ifelse(
    full_data_discovered_gene & thinned_refit_discovered_gene,
    "Recovered full-data discovery",
    ifelse(
      full_data_discovered_gene & !thinned_refit_discovered_gene,
      "Missed full-data discovery",
      ifelse(
        !full_data_discovered_gene & thinned_refit_discovered_gene,
        "Additional thinned-refit discovery",
        "Neither"
      )
    )
  )
)

gene_discovery_summary <- data.frame(
  selection_rule = "minimum_full_bf_lfdr",
  n_selected_pairs = nrow(selection),
  full_discovered_pairs = length(full_call_indices),
  full_discovered_genes = length(full_discovered_genes),
  selected_full_data_fdr_call_pairs =
    sum(selection$selected_pair_full_data_fdr_call),
  thinned_discovered_genes = length(thinned_discovered_genes),
  recovered_full_discovered_genes = length(overlap_genes),
  missed_full_discovered_genes = length(missed_genes),
  additional_thinned_discovered_genes = length(additional_genes),
  gene_recall = safe_ratio(length(overlap_genes), length(full_discovered_genes)),
  gene_overlap_fraction = safe_ratio(
    length(overlap_genes),
    length(thinned_discovered_genes)
  ),
  gene_jaccard = safe_ratio(length(overlap_genes), length(gene_union)),
  alpha = alpha,
  stringsAsFactors = FALSE
)

configuration <- list(
  experiment = paste(
    "One minimum full-data BF-adjusted lfdr variant per gene FASH(1) refit"
  ),
  selection_rule = "minimum_full_bf_lfdr",
  selection_source = "Full-data BF-adjusted lfdr",
  outcome_informed = TRUE,
  n_full_units = n_full_units,
  n_genes = n_genes,
  n_selected_units = nrow(selection),
  alpha = alpha,
  refit_strategy = paste(
    "Subset fixed per-unit likelihood rows from the full raw fit and rerun",
    "fash_eb_est with penalty 10, followed by BF_update."
  ),
  full_raw_settings = full_raw_settings,
  full_psd_grid = full_psd_grid,
  source_files = source_files,
  raw_fit_warnings = raw_capture$warnings,
  bf_update_warnings = bf_capture$warnings,
  raw_eb_refit_elapsed_seconds = unname(raw_elapsed),
  bf_update_elapsed_seconds = unname(bf_elapsed),
  r_version = R.version.string,
  package_versions = c(fashr = as.character(utils::packageVersion("fashr"))),
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
fit_bundle <- list(
  configuration = configuration,
  selection = selection,
  raw_fit = thinned_raw,
  bf_adjusted_fit = thinned_bf,
  full_reference = list(
    raw_prior_weights = full_raw_prior_weights,
    bf_adjusted_prior_weights = full_bf_prior_weights,
    raw_selected_lfdr = raw_lfdr$table[, c("pair_key", "full_lfdr")],
    bf_adjusted_selected_lfdr = bf_lfdr$table[, c(
      "pair_key",
      "full_lfdr"
    )]
  )
)

dir.create(staging_directory, recursive = FALSE)

saveRDS(configuration, file.path(staging_directory, "configuration.rds"))
saveRDS(fit_bundle, file.path(staging_directory, "thinned_fash_fit.rds"))
write_csv(selection, file.path(staging_directory, "selection.csv"))
write_csv(
  prior_table,
  file.path(staging_directory, "prior_weight_comparison.csv")
)
write_csv(
  lfdr_table,
  file.path(staging_directory, "lfdr_comparison.csv")
)
write_csv(
  comparison_summary,
  file.path(staging_directory, "comparison_summary.csv")
)
write_csv(
  gene_discovery_comparison,
  file.path(staging_directory, "gene_discovery_comparison.csv")
)
write_csv(
  gene_discovery_summary,
  file.path(staging_directory, "gene_discovery_summary.csv")
)

if (!file.rename(staging_directory, final_directory)) {
  stop("Could not finalize output directory: ", final_directory)
}

cat(
  "\nBest-lfdr one-variant-per-gene FASH refit completed: ",
  final_directory,
  "\n",
  sep = ""
)
