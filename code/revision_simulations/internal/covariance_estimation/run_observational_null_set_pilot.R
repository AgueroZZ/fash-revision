#!/usr/bin/env Rscript

# Run a small residual-covariance pilot on two observational null-enriched sets.

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

pilot_size <- as.integer(get_arg("--pilot-size", "200"))
n_bootstrap <- as.integer(get_arg("--bootstrap-reps", "200"))
selection_seed <- as.integer(get_arg("--selection-seed", "20260806"))
bootstrap_seed <- as.integer(get_arg("--bootstrap-seed", "20260816"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
output_id <- get_arg("--output-id", "observational_null_set_pilot")
if (is.na(pilot_size) || pilot_size < 20L ||
    is.na(n_bootstrap) || n_bootstrap < 20L ||
    is.na(selection_seed) || is.na(bootstrap_seed) ||
    !is.finite(alpha) || alpha <= 0 || alpha >= 1 || !nzchar(output_id)) {
  stop("Invalid pilot arguments.")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required for the time-specific screens.")
}

fit_path <- file.path(
  workflowr_root,
  "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
)
time_specific_dir <- file.path(
  workflowr_root,
  "data", "dynamic_eQTL_real", "strober_nondyn"
)
time_grid <- 0:15
time_specific_paths <- file.path(
  time_specific_dir,
  paste0(
    "non_dynamic_time_", time_grid,
    "_eqtl_results_3_pc.txt"
  )
)
required_inputs <- c(fit_path, time_specific_paths)
if (any(!file.exists(required_inputs))) {
  stop(
    "At least one required pilot input is missing: ",
    paste(required_inputs[!file.exists(required_inputs)], collapse = ", ")
  )
}

output_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal", output_id
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
if (length(pair_keys) != 1009173L ||
    length(fash_fit$lfdr) != length(pair_keys) ||
    length(fash_fit$fash_data$S) != length(pair_keys) ||
    any(!is.finite(fash_fit$lfdr)) || anyDuplicated(pair_keys)) {
  stop("The fitted FASH object is not the expected complete fit.")
}

message("Constructing the pair catalog and FASH dynamic-discovery screen.")
catalog <- make_randomized_pair_catalog(pair_keys, seed = selection_seed)
all_genes <- sort(unique(catalog$gene_id), method = "radix")
dynamic_discovery <- discover_by_cumulative_lfdr(fash_fit$lfdr, alpha)
dynamic_discovery_genes <- sort(unique(
  catalog$gene_id[dynamic_discovery$discovered]
), method = "radix")

message("Reading 16 time-specific eQTL discovery screens.")
time_specific_discovery_genes <- character()
time_specific_screen_rows <- vector("list", length(time_grid))
for (time_index in seq_along(time_grid)) {
  time_value <- time_grid[time_index]
  screen <- data.table::fread(
    time_specific_paths[time_index],
    select = c("ensamble_id", "eFDR"),
    showProgress = FALSE
  )
  if (nrow(screen) != length(pair_keys) ||
      !all(c("ensamble_id", "eFDR") %in% names(screen))) {
    stop("Unexpected time-specific screen structure at day ", time_value, ".")
  }
  significant <- !is.na(screen$eFDR) & screen$eFDR <= alpha
  significant_genes <- unique(as.character(screen$ensamble_id[significant]))
  time_specific_discovery_genes <- union(
    time_specific_discovery_genes,
    significant_genes
  )
  time_specific_screen_rows[[time_index]] <- data.frame(
    time = time_value,
    n_pairs = nrow(screen),
    n_missing_efdr = sum(is.na(screen$eFDR)),
    n_discovered_pairs = sum(significant),
    n_discovered_genes = length(significant_genes),
    stringsAsFactors = FALSE
  )
  rm(screen)
}
time_specific_discovery_genes <- sort(
  unique(time_specific_discovery_genes),
  method = "radix"
)
time_specific_screen_summary <- do.call(rbind, time_specific_screen_rows)

eligible_set_a <- setdiff(all_genes, dynamic_discovery_genes)
eligible_set_b <- setdiff(eligible_set_a, time_specific_discovery_genes)
selection_a <- select_random_unique_variant_per_gene(
  catalog,
  eligible_set_a,
  n_select = pilot_size,
  set_id = "A_dynamic_null"
)
selection_b <- select_random_unique_variant_per_gene(
  catalog,
  eligible_set_b,
  n_select = pilot_size,
  set_id = "B_zero_eqtl_null"
)
selected_units <- rbind(selection_a$selected, selection_b$selected)
selected_units$bf_adjusted_lfdr <- fash_fit$lfdr[selected_units$fash_index]
selected_units$gene_has_dynamic_discovery <-
  selected_units$gene_id %in% dynamic_discovery_genes
selected_units$gene_has_time_specific_discovery <-
  selected_units$gene_id %in% time_specific_discovery_genes

if (any(selected_units$gene_has_dynamic_discovery) ||
    any(selected_units$set_id == "B_zero_eqtl_null" &
          selected_units$gene_has_time_specific_discovery)) {
  stop("At least one selected gene violates its null-enrichment screen.")
}

message("Extracting selected beta estimates and t-adjusted standard errors.")
extracted_a <- extract_fash_beta_se(fash_fit, selection_a$selected_indices)
extracted_b <- extract_fash_beta_se(fash_fit, selection_b$selected_indices)
if (!isTRUE(all.equal(extracted_a$time_grid, time_grid)) ||
    !isTRUE(all.equal(extracted_b$time_grid, time_grid))) {
  stop("The selected FASH units do not use the expected 0:15 time grid.")
}

message("Estimating raw covariance patterns and gene-bootstrap intervals.")
bootstrap_a <- bootstrap_covariance_lags(
  extracted_a$beta_hat,
  extracted_a$se,
  estimators = c(
    "Pairwise difference",
    "Within-unit centered covariance"
  ),
  n_bootstrap = n_bootstrap,
  seed = bootstrap_seed
)
bootstrap_b <- bootstrap_covariance_lags(
  extracted_b$beta_hat,
  extracted_b$se,
  estimators = c("Pairwise difference", "Direct zero mean"),
  n_bootstrap = n_bootstrap,
  seed = bootstrap_seed + 1L
)
centered_a <- estimate_within_unit_centered_matrices(
  extracted_a$beta_hat,
  extracted_a$se
)

matrices <- list(
  A_dynamic_null__pairwise =
    bootstrap_a$observed_matrices[["Pairwise difference"]],
  A_dynamic_null__within_unit_centered_covariance =
    bootstrap_a$observed_matrices[["Within-unit centered covariance"]],
  B_zero_eqtl_null__pairwise =
    bootstrap_b$observed_matrices[["Pairwise difference"]],
  B_zero_eqtl_null__direct =
    bootstrap_b$observed_matrices[["Direct zero mean"]]
)
matrix_metadata <- data.frame(
  matrix_id = names(matrices),
  set_id = c(
    "A_dynamic_null", "A_dynamic_null",
    "B_zero_eqtl_null", "B_zero_eqtl_null"
  ),
  estimator = c(
    "Pairwise difference", "Within-unit centered covariance",
    "Pairwise difference", "Direct zero mean"
  ),
  stringsAsFactors = FALSE
)

matrix_diagnostics <- vector("list", length(matrices))
matrix_long <- vector("list", 2L * length(matrices))
lag_summaries <- vector("list", length(matrices))
beta_scale_matrices <- vector("list", length(matrices))
names(beta_scale_matrices) <- names(matrices)
for (matrix_index in seq_along(matrices)) {
  matrix_id <- names(matrices)[matrix_index]
  metadata <- matrix_metadata[matrix_index, ]
  selected_se <- if (metadata$set_id == "A_dynamic_null") {
    extracted_a$se
  } else {
    extracted_b$se
  }
  covariance <- matrices[[matrix_index]]
  beta_scale <- scale_covariance_by_median_se(covariance, selected_se)
  beta_scale_matrices[[matrix_id]] <- beta_scale

  matrix_diagnostics[[matrix_index]] <- covariance_matrix_diagnostics(
    covariance,
    metadata$set_id,
    metadata$estimator
  )
  standardized_long <- covariance_matrix_to_long(
    covariance,
    metadata$set_id,
    metadata$estimator,
    scale = "Standardized residual"
  )
  standardized_long$matrix_id <- matrix_id
  beta_long <- covariance_matrix_to_long(
    beta_scale,
    metadata$set_id,
    metadata$estimator,
    scale = "Median beta-hat scale"
  )
  beta_long$matrix_id <- matrix_id
  matrix_long[[2L * matrix_index - 1L]] <- standardized_long
  matrix_long[[2L * matrix_index]] <- beta_long

  lag_summary <- lag_covariance_variogram(covariance, selected_se)
  lag_summary$matrix_id <- matrix_id
  lag_summary$set_id <- metadata$set_id
  lag_summary$estimator <- metadata$estimator
  lag_summaries[[matrix_index]] <- lag_summary
}
matrix_diagnostics <- do.call(rbind, matrix_diagnostics)
covariance_matrices_long <- do.call(rbind, matrix_long)
lag_summaries <- do.call(rbind, lag_summaries)
rownames(matrix_diagnostics) <- NULL
rownames(covariance_matrices_long) <- NULL
rownames(lag_summaries) <- NULL

bootstrap_summary_a <- bootstrap_a$summary
bootstrap_summary_a$set_id <- "A_dynamic_null"
bootstrap_summary_b <- bootstrap_b$summary
bootstrap_summary_b$set_id <- "B_zero_eqtl_null"
bootstrap_lag_intervals <- rbind(bootstrap_summary_a, bootstrap_summary_b)
rownames(bootstrap_lag_intervals) <- NULL

selection_counts <- data.frame(
  set_id = c("A_dynamic_null", "B_zero_eqtl_null"),
  definition = c(
    "No gene-level BF-adjusted FASH dynamic discovery at cumulative-lfdr FDR 0.05",
    paste(
      "Set A plus no gene-level time-specific eQTL discovery at eFDR 0.05",
      "on any of 16 days"
    )
  ),
  n_eligible_genes = c(length(eligible_set_a), length(eligible_set_b)),
  n_eligible_genes_with_gene_unique_variant = c(
    selection_a$n_eligible_genes_with_unique_variant,
    selection_b$n_eligible_genes_with_unique_variant
  ),
  n_selected = c(selection_a$n_selected, selection_b$n_selected),
  n_candidate_pairs = c(
    selection_a$n_candidate_pairs,
    selection_b$n_candidate_pairs
  ),
  stringsAsFactors = FALSE
)

fit_info <- file.info(fit_path)
screen_info <- file.info(time_specific_paths)
configuration <- list(
  analysis_id = output_id,
  created_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  status = "Internal small-scale pilot",
  fit_path = fit_path,
  fit_object = "fash_fit1_update",
  fit_size_bytes = unname(fit_info$size),
  fit_mtime_utc = format(fit_info$mtime, tz = "UTC", usetz = TRUE),
  time_specific_paths = time_specific_paths,
  time_specific_sizes_bytes = unname(screen_info$size),
  time_specific_mtimes_utc = format(screen_info$mtime, tz = "UTC", usetz = TRUE),
  n_fash_pairs = length(pair_keys),
  n_genes = length(all_genes),
  n_gene_unique_variants = sum(catalog$variant_gene_count == 1L),
  alpha = alpha,
  dynamic_rule = paste(
    "BF-adjusted FASH cumulative-lfdr FDR at alpha",
    format(alpha)
  ),
  time_specific_rule = paste(
    "Published time-specific eFDR at alpha", format(alpha),
    "on each of 16 days"
  ),
  selection_rule = paste(
    "One globally gene-unique variant per eligible gene; fixed random scores",
    "generated without using beta estimates, SE, lfdr, p-value, or eFDR"
  ),
  set_a_definition = selection_counts$definition[1L],
  set_b_definition = selection_counts$definition[2L],
  pilot_size = pilot_size,
  selection_seed = selection_seed,
  n_bootstrap = n_bootstrap,
  bootstrap_seed_set_a = bootstrap_seed,
  bootstrap_seed_set_b = bootstrap_seed + 1L,
  time_grid = time_grid,
  se_scale = "FASH t-adjusted standard error",
  pairwise_estimator = paste(
    "OLS through the origin using unequal-SE pairwise differences;",
    "diagonal fixed at one"
  ),
  direct_estimator = paste(
    "Uncentered empirical second moment of z = beta_hat / SE;",
    "used only in the stricter zero-eQTL-null set"
  ),
  centered_estimator = paste(
    "Inverse-SE-squared within-unit constant removal, followed by",
    "crossprod of standardized residuals divided by the number of units"
  ),
  matrix_policy = "Raw unprojected estimates; no nearest-PD repair",
  r_version = R.version.string,
  fashr_version = as.character(utils::packageVersion("fashr")),
  data_table_version = as.character(utils::packageVersion("data.table"))
)

analysis_result <- list(
  configuration = configuration,
  selection_counts = selection_counts,
  selected_units = selected_units,
  time_specific_screen_summary = time_specific_screen_summary,
  discovery_summary = data.frame(
    n_dynamic_discovered_pairs = sum(dynamic_discovery$discovered),
    n_dynamic_discovered_genes = length(dynamic_discovery_genes),
    dynamic_lfdr_cutoff = dynamic_discovery$cutoff_lfdr,
    n_time_specific_discovered_genes_union =
      length(time_specific_discovery_genes),
    stringsAsFactors = FALSE
  ),
  selected_data = list(
    set_a = list(
      beta_hat = extracted_a$beta_hat,
      adjusted_se = extracted_a$se,
      time_grid = extracted_a$time_grid
    ),
    set_b = list(
      beta_hat = extracted_b$beta_hat,
      adjusted_se = extracted_b$se,
      time_grid = extracted_b$time_grid
    )
  ),
  matrices = matrices,
  centered_correlations = list(
    A_dynamic_null__within_unit_centered_correlation =
      centered_a$correlation
  ),
  beta_scale_matrices = beta_scale_matrices,
  matrix_metadata = matrix_metadata,
  matrix_diagnostics = matrix_diagnostics,
  covariance_matrices_long = covariance_matrices_long,
  lag_summaries = lag_summaries,
  bootstrap_lag_intervals = bootstrap_lag_intervals,
  bootstrap_draws = list(set_a = bootstrap_a$draws, set_b = bootstrap_b$draws)
)

saveRDS(configuration, file.path(output_dir, "configuration.rds"))
saveRDS(
  analysis_result,
  file.path(output_dir, "observational_null_set_pilot.rds")
)
saveRDS(matrices, file.path(output_dir, "raw_standardized_covariance_matrices.rds"))
saveRDS(
  beta_scale_matrices,
  file.path(output_dir, "median_beta_scale_covariance_matrices.rds")
)
write_csv(selection_counts, file.path(summary_dir, "selection_counts.csv"))
write_csv(selected_units, file.path(summary_dir, "selected_units.csv"))
write_csv(
  time_specific_screen_summary,
  file.path(summary_dir, "time_specific_screen_summary.csv")
)
write_csv(
  analysis_result$discovery_summary,
  file.path(summary_dir, "discovery_summary.csv")
)
write_csv(
  matrix_diagnostics,
  file.path(summary_dir, "matrix_diagnostics.csv")
)
write_csv(
  covariance_matrices_long,
  file.path(summary_dir, "covariance_matrices_long.csv")
)
write_csv(lag_summaries, file.path(summary_dir, "lag_summaries.csv"))
write_csv(
  bootstrap_lag_intervals,
  file.path(summary_dir, "bootstrap_lag_intervals.csv")
)

cat("\nSelection counts:\n")
print(selection_counts)
cat("\nDiscovery summary:\n")
print(analysis_result$discovery_summary)
cat("\nMatrix diagnostics:\n")
print(matrix_diagnostics)
cat("\nLag summaries at lags 1, 5, 10, and 15:\n")
print(lag_summaries[lag_summaries$lag %in% c(1L, 5L, 10L, 15L), ])
cat("\nPilot completed: ", output_dir, "\n", sep = "")
