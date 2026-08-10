#!/usr/bin/env Rscript

# Estimate pairwise-difference residual correlations from gene-level minimum-lfdr sets.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(args, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(args[equals_hit[1L]], nchar(equals_prefix) + 1L))
  }
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1L] == length(args)) {
    return(default)
  }
  args[hit[1L] + 1L]
}

parse_numeric_list <- function(value, name) {
  parsed <- as.numeric(strsplit(value, ",", fixed = TRUE)[[1L]])
  if (length(parsed) == 0L || any(!is.finite(parsed))) {
    stop("Invalid ", name, ".")
  }
  parsed
}

threshold_id <- function(threshold) {
  sub("\\.", "p", format(threshold, nsmall = 3L, trim = TRUE))
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "r4_correlated_errors",
  "real_data_correlation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "covariance_estimation",
  "observational_null_set_helpers.R"
))

thresholds <- sort(unique(parse_numeric_list(
  get_arg("--thresholds", "0.90,0.925"),
  "thresholds"
)))
n_bootstrap <- as.integer(get_arg("--bootstrap-reps", "500"))
selection_seed <- as.integer(get_arg("--selection-seed", "20260809"))
bootstrap_seed <- as.integer(get_arg("--bootstrap-seed", "20260819"))
output_id <- get_arg("--output-id", "gene_min_lfdr_covariance")
if (length(thresholds) != 2L || any(thresholds <= 0 | thresholds >= 1) ||
    is.na(n_bootstrap) || n_bootstrap < 20L || is.na(selection_seed) ||
    is.na(bootstrap_seed) || !nzchar(output_id)) {
  stop("Invalid pilot arguments.")
}

fit_path <- file.path(
  workflowr_root,
  "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
)
if (!file.exists(fit_path)) {
  stop("The BF-adjusted FASH fit is missing: ", fit_path)
}

output_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "internal", output_id
)
summary_dir <- file.path(output_dir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading the BF-adjusted FASH fit.")
fit_environment <- new.env(parent = emptyenv())
loaded_names <- load(fit_path, envir = fit_environment)
if (!identical(loaded_names, "fash_fit1_update")) {
  stop("The fitted-object file must contain only fash_fit1_update.")
}
fash_fit <- fit_environment$fash_fit1_update
pair_keys <- names(fash_fit$fash_data$data_list)
lfdr <- as.numeric(fash_fit$lfdr)
if (length(pair_keys) != length(lfdr) ||
    length(pair_keys) != length(fash_fit$fash_data$S) ||
    any(!is.finite(lfdr)) || any(lfdr < 0 | lfdr > 1) ||
    anyDuplicated(pair_keys)) {
  stop("The fitted FASH object is incomplete or invalid.")
}

gene_id <- parse_gene_ids(pair_keys)
min_lfdr <- tapply(lfdr, gene_id, min)
gene_summary <- data.frame(
  gene_id = names(min_lfdr),
  min_lfdr = as.numeric(min_lfdr),
  stringsAsFactors = FALSE
)

message("Constructing the outcome-independent variant catalog.")
catalog <- make_randomized_pair_catalog(pair_keys, seed = selection_seed)

matrices <- list()
lag_summaries <- list()
bootstrap_summaries <- list()
diagnostics <- list()
selection_counts <- list()
selected_units <- list()
bootstrap_results <- list()

for (threshold_index in seq_along(thresholds)) {
  threshold <- thresholds[threshold_index]
  id <- paste0("min_lfdr_gt_", threshold_id(threshold))
  eligible_genes <- names(min_lfdr)[min_lfdr > threshold]
  message(
    "Selecting globally gene-unique variants for ", length(eligible_genes),
    " genes with m_g > ", threshold, "."
  )
  selection <- select_random_unique_variant_per_gene(
    catalog,
    eligible_genes = eligible_genes,
    n_select = length(eligible_genes),
    set_id = id
  )
  selected <- selection$selected
  selected$gene_min_lfdr <- unname(min_lfdr[selected$gene_id])
  if (any(selected$gene_min_lfdr <= threshold) ||
      anyDuplicated(selected$gene_id) || anyDuplicated(selected$variant_id)) {
    stop("The gene-level minimum-lfdr selection failed its invariants.")
  }

  extracted <- extract_fash_beta_se(fash_fit, selected$fash_index)
  if (!isTRUE(all.equal(extracted$time_grid, 0:15))) {
    stop("The selected FASH data do not use the expected time grid 0:15.")
  }

  current_bootstrap_seed <- bootstrap_seed + threshold_index - 1L
  message(
    "Running ", n_bootstrap,
    " gene-bootstrap replications for ", id, "."
  )
  bootstrap <- bootstrap_covariance_lags(
    extracted$beta_hat,
    extracted$se,
    estimators = "Pairwise difference",
    n_bootstrap = n_bootstrap,
    seed = current_bootstrap_seed
  )
  covariance <- bootstrap$observed_matrices[["Pairwise difference"]]
  lag_summary <- lag_covariance_variogram(covariance, extracted$se)
  lag_summary$set_id <- id
  lag_summary$threshold <- threshold
  lag_summary$n_selected <- nrow(selected)
  lag_summary <- lag_summary[, c(
    "set_id", "threshold", "n_selected", "lag",
    "mean_standardized_covariance", "mean_semivariogram",
    "mean_median_beta_scale_covariance",
    "mean_median_beta_scale_semivariance"
  )]

  diagnostics[[id]] <- cbind(
    data.frame(
      set_id = id,
      threshold = threshold,
      n_genes_passing_min_lfdr = length(eligible_genes),
      n_genes_with_unique_variant = selection$n_eligible_genes_with_unique_variant,
      n_selected = nrow(selected),
      minimum_selected_gene_min_lfdr = min(selected$gene_min_lfdr),
      maximum_selected_gene_min_lfdr = max(selected$gene_min_lfdr),
      stringsAsFactors = FALSE
    ),
    covariance_matrix_diagnostics(
      covariance,
      set_id = id,
      estimator = "Pairwise difference"
    )[, c(
      "minimum_eigenvalue", "maximum_eigenvalue",
      "n_negative_eigenvalues", "minimum_diagonal", "maximum_diagonal",
      "mean_off_diagonal", "minimum_off_diagonal", "maximum_off_diagonal"
    )]
  )
  diagnostics[[id]]$lag_1 <- lag_summary$mean_standardized_covariance[1L]
  diagnostics[[id]]$lag_15 <- lag_summary$mean_standardized_covariance[15L]

  selection_counts[[id]] <- data.frame(
    set_id = id,
    threshold = threshold,
    n_genes_passing_min_lfdr = length(eligible_genes),
    n_genes_with_unique_variant = selection$n_eligible_genes_with_unique_variant,
    n_selected = nrow(selected),
    stringsAsFactors = FALSE
  )
  selected_units[[id]] <- selected
  matrices[[id]] <- covariance
  lag_summaries[[id]] <- lag_summary
  bootstrap_summaries[[id]] <- bootstrap$summary
  bootstrap_summaries[[id]]$set_id <- id
  bootstrap_summaries[[id]]$threshold <- threshold
  bootstrap_summaries[[id]]$n_selected <- nrow(selected)
  bootstrap_results[[id]] <- list(
    seed = current_bootstrap_seed,
    draws = bootstrap$draws
  )
}

if (!all(selected_units[[2L]]$gene_id %in% selected_units[[1L]]$gene_id)) {
  stop("The stricter gene-level set is not nested in the broader set.")
}

configuration <- list(
  analysis_id = output_id,
  created_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  fit_path = fit_path,
  fit_object = "fash_fit1_update",
  n_fash_pairs = length(pair_keys),
  n_genes = length(min_lfdr),
  thresholds = thresholds,
  threshold_rule = "Strictly greater than gene-level minimum lfdr",
  gene_statistic = "Minimum pair-level BF-adjusted FASH lfdr within gene",
  variant_selection = paste(
    "All available globally gene-unique variants were cataloged; one",
    "variant per retained gene was selected by a deterministic random score",
    "that does not use lfdr values."
  ),
  estimator = paste(
    "Existing OLS-through-origin pairwise-difference estimator on the",
    "standardized beta-hat/t-adjusted-SE scale."
  ),
  matrix_version = "Raw unprojected pairwise-difference estimate",
  se_scale = "t-adjusted",
  n_bootstrap = n_bootstrap,
  selection_seed = selection_seed,
  bootstrap_seed = bootstrap_seed,
  time_grid = 0:15,
  r_version = R.version.string,
  fashr_version = as.character(utils::packageVersion("fashr"))
)

analysis_result <- list(
  configuration = configuration,
  gene_min_lfdr = gene_summary,
  selection_counts = do.call(rbind, selection_counts),
  selected_units = selected_units,
  matrices = matrices,
  lag_summaries = do.call(rbind, lag_summaries),
  diagnostics = do.call(rbind, diagnostics),
  bootstrap_summaries = do.call(rbind, bootstrap_summaries),
  bootstrap_results = bootstrap_results
)

saveRDS(configuration, file.path(output_dir, "configuration.rds"))
saveRDS(analysis_result, file.path(output_dir, "gene_min_lfdr_covariance.rds"))
saveRDS(matrices, file.path(output_dir, "raw_pairwise_correlation_matrices.rds"))
write_csv(gene_summary, file.path(summary_dir, "gene_min_lfdr.csv"))
write_csv(analysis_result$selection_counts, file.path(summary_dir, "selection_counts.csv"))
write_csv(analysis_result$diagnostics, file.path(summary_dir, "matrix_diagnostics.csv"))
write_csv(analysis_result$lag_summaries, file.path(summary_dir, "lag_summaries.csv"))
write_csv(analysis_result$bootstrap_summaries, file.path(summary_dir, "bootstrap_summaries.csv"))
write_csv(do.call(rbind, selected_units), file.path(summary_dir, "selected_units.csv"))

cat("\nSelection counts:\n")
print(analysis_result$selection_counts)
cat("\nMatrix diagnostics:\n")
print(analysis_result$diagnostics)
cat("\nGene-minimum-lfdr covariance pilot completed: ", output_dir, "\n", sep = "")
