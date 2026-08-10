#!/usr/bin/env Rscript

# Run the fixed-subset FASH sensitivity analysis for diagonal, C1, and C2 errors.

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
  arguments <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(arguments, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(arguments[equals_hit[1L]], nchar(equals_prefix) + 1L))
  }
  hit <- which(arguments == name)
  if (length(hit) == 0L || hit[1L] == length(arguments)) {
    return(default)
  }
  arguments[hit[1L] + 1L]
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

capture_warnings <- function(expression) {
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

file_metadata <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required for SHA-256 provenance.")
  }
  information <- file.info(path)
  if (nrow(information) != 1L || is.na(information$size)) {
    stop("Could not read source-file metadata: ", path)
  }
  list(
    path = normalizePath(path, mustWork = TRUE),
    size_bytes = unname(information$size),
    modification_time = format(information$mtime, tz = "UTC", usetz = TRUE),
    sha256 = digest::digest(path, algo = "sha256", file = TRUE)
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

validate_fash_fit <- function(fit, pair_keys, psd_grid, name) {
  required <- c(
    "prior_weights", "posterior_weights", "psd_grid", "lfdr",
    "settings", "fash_data", "L_matrix"
  )
  if (!inherits(fit, "fash") || !all(required %in% names(fit)) ||
      !identical(as.numeric(fit$psd_grid), as.numeric(psd_grid)) ||
      length(fit$lfdr) != length(pair_keys) ||
      nrow(fit$posterior_weights) != length(pair_keys) ||
      nrow(fit$L_matrix) != length(pair_keys) ||
      ncol(fit$L_matrix) != length(psd_grid) ||
      any(!is.finite(fit$lfdr)) || any(fit$lfdr < 0 | fit$lfdr > 1) ||
      anyNA(fit$L_matrix) || any(is.nan(fit$L_matrix)) ||
      any(fit$L_matrix == Inf)) {
    stop(name, " is not a valid aligned FASH fit.")
  }
  names(fit$lfdr) <- pair_keys
  rownames(fit$posterior_weights) <- pair_keys
  rownames(fit$L_matrix) <- pair_keys
  extract_fit_prior_weights(fit, paste0(name, " prior weights"))
  fit
}

matrix_to_long <- function(matrix, matrix_id, matrix_label) {
  matrix <- as.matrix(matrix)
  grid <- expand.grid(
    time_a = seq_len(nrow(matrix)) - 1L,
    time_b = seq_len(ncol(matrix)) - 1L
  )
  grid$correlation <- matrix[cbind(grid$time_a + 1L, grid$time_b + 1L)]
  grid$matrix_id <- matrix_id
  grid$matrix_label <- matrix_label
  grid[, c("matrix_id", "matrix_label", "time_a", "time_b", "correlation")]
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "one_variant_per_gene_refit",
  "one_variant_per_gene_refit_helpers.R"
))
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "correlated_likelihood_sensitivity",
  "correlated_likelihood_helpers.R"
))

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required.")
}

num_cores <- as.integer(get_arg("--num-cores", "8"))
validation_units <- as.integer(get_arg("--validation-units", "24"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
expected_thinning_seed <- as.integer(get_arg("--thinning-seed", "20260811"))
output_id <- get_arg(
  "--output-id",
  "correlated_likelihood_sensitivity"
)
if (anyNA(c(num_cores, validation_units, alpha, expected_thinning_seed)) ||
    num_cores < 1L || validation_units < 2L || alpha <= 0 || alpha >= 1 ||
    !nzchar(output_id) || grepl("/", output_id, fixed = TRUE)) {
  stop("Invalid analysis arguments.")
}

output_parent <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal"
)
final_directory <- file.path(output_parent, output_id)
if (file.exists(final_directory)) {
  stop("Refusing to overwrite the existing output directory: ", final_directory)
}
staging_directory <- file.path(
  output_parent,
  paste0(".", output_id, "_staging_", Sys.getpid())
)
if (file.exists(staging_directory)) {
  stop("Unexpected staging-directory collision: ", staging_directory)
}
dir.create(staging_directory, recursive = FALSE)
summary_directory <- file.path(staging_directory, "summary")
dir.create(summary_directory, recursive = FALSE)

covariance_directory <- file.path(
  output_parent,
  "mashr_mean_z_null_correlation"
)
covariance_configuration_path <- file.path(
  covariance_directory,
  "configuration.rds"
)
covariance_analysis_path <- file.path(
  covariance_directory,
  "mashr_mean_z_null_correlation.rds"
)
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
required_paths <- c(
  covariance_configuration_path,
  covariance_analysis_path,
  raw_fit_path,
  bf_fit_path
)
if (any(!file.exists(required_paths))) {
  stop("One or more required input caches are missing.")
}
source_files <- lapply(required_paths, file_metadata)
names(source_files) <- c(
  "covariance_configuration",
  "covariance_analysis",
  "raw_fash_fit",
  "bf_adjusted_fash_fit"
)

message("Loading the fixed thinning and C1/C2 correlation estimates.")
covariance_configuration <- readRDS(covariance_configuration_path)
covariance_analysis <- readRDS(covariance_analysis_path)
if (!identical(as.integer(covariance_configuration$thinning_seed),
               expected_thinning_seed) ||
    !identical(as.integer(covariance_configuration$n_genes), 6362L) ||
    !identical(as.integer(covariance_configuration$n_time), 16L)) {
  stop("The covariance cache does not use the confirmed thinning configuration.")
}
pair_metadata <- covariance_analysis$candidate_metadata
beta_hat <- as.matrix(covariance_analysis$candidate_beta_hat)
adjusted_se <- as.matrix(covariance_analysis$candidate_adjusted_se)
time_grid <- as.numeric(covariance_configuration$time_grid)
required_metadata <- c(
  "fash_index", "pair_key", "gene_id", "variant_id"
)
if (!is.data.frame(pair_metadata) ||
    !all(required_metadata %in% names(pair_metadata)) ||
    nrow(pair_metadata) != 6362L || anyDuplicated(pair_metadata$pair_key) ||
    anyDuplicated(pair_metadata$gene_id) ||
    !identical(dim(beta_hat), c(6362L, 16L)) ||
    !identical(dim(adjusted_se), dim(beta_hat)) ||
    !identical(rownames(beta_hat), pair_metadata$pair_key) ||
    !identical(rownames(adjusted_se), pair_metadata$pair_key) ||
    any(!is.finite(beta_hat)) || any(!is.finite(adjusted_se)) ||
    any(adjusted_se <= 0) || !identical(time_grid, as.numeric(0:15))) {
  stop("The fixed candidate matrices or metadata are invalid.")
}

C1 <- covariance_analysis$estimates$maximum_z$estimate$sample_correlation
C2 <- covariance_analysis$estimates$maximum_z_and_mean_z$estimate$sample_correlation
C1_diagnostics <- validate_shared_correlation(C1, 16L, "C1")
C2_diagnostics <- validate_shared_correlation(C2, 16L, "C2")
correlation_diagnostics <- data.frame(
  matrix_id = c("C1", "C2"),
  matrix_label = c("C1: Screen 1", "C2: Screen 2"),
  n_estimating_pairs = c(
    covariance_analysis$selection_counts$n_selected[
      covariance_analysis$selection_counts$filter_id == "maximum_z"
    ],
    covariance_analysis$selection_counts$n_selected[
      covariance_analysis$selection_counts$filter_id ==
        "maximum_z_and_mean_z"
    ]
  ),
  minimum_eigenvalue = c(
    C1_diagnostics$minimum_eigenvalue,
    C2_diagnostics$minimum_eigenvalue
  ),
  maximum_eigenvalue = c(
    C1_diagnostics$maximum_eigenvalue,
    C2_diagnostics$maximum_eigenvalue
  ),
  eigenvalue_condition_number = c(
    C1_diagnostics$eigenvalue_condition_number,
    C2_diagnostics$eigenvalue_condition_number
  ),
  log_determinant = c(
    C1_diagnostics$log_determinant,
    C2_diagnostics$log_determinant
  ),
  stringsAsFactors = FALSE
)
correlation_matrices_long <- rbind(
  matrix_to_long(C1, "C1", "C1: Screen 1"),
  matrix_to_long(C2, "C2", "C2: Screen 2")
)

message("Loading and validating the original raw FASH fit.")
full_raw <- load_exact_object(raw_fit_path, "fash_fit1")
required_fit_fields <- c(
  "prior_weights", "posterior_weights", "psd_grid", "lfdr",
  "settings", "fash_data", "L_matrix"
)
if (!all(required_fit_fields %in% names(full_raw)) ||
    is.null(full_raw$fash_data$data_list) || is.null(full_raw$fash_data$S) ||
    length(full_raw$fash_data$data_list) != nrow(full_raw$L_matrix) ||
    length(full_raw$fash_data$S) != nrow(full_raw$L_matrix)) {
  stop("The original raw FASH fit is incomplete.")
}
settings <- validate_fash_settings(full_raw$settings)
if (!identical(as.numeric(settings$penalty), 10)) {
  stop("The confirmed sensitivity analysis requires the original penalty 10.")
}
psd_grid <- as.numeric(full_raw$psd_grid)
full_pair_keys <- names(full_raw$fash_data$data_list)
selected_indices <- as.integer(pair_metadata$fash_index)
if (anyNA(selected_indices) || anyDuplicated(selected_indices) ||
    any(selected_indices < 1L | selected_indices > length(full_pair_keys)) ||
    !identical(full_pair_keys[selected_indices], pair_metadata$pair_key)) {
  stop("The covariance-cache selection is not aligned to the original FASH fit.")
}
n_time <- ncol(beta_hat)
selected_beta_from_raw <- t(vapply(
  full_raw$fash_data$data_list[selected_indices],
  function(dataset) as.numeric(dataset$y),
  numeric(n_time)
))
selected_time_from_raw <- t(vapply(
  full_raw$fash_data$data_list[selected_indices],
  function(dataset) as.numeric(dataset$x),
  numeric(n_time)
))
selected_se_from_raw <- t(vapply(
  full_raw$fash_data$S[selected_indices],
  as.numeric,
  numeric(n_time)
))
if (max(abs(selected_beta_from_raw - beta_hat)) > 1e-12 ||
    max(abs(selected_se_from_raw - adjusted_se)) > 1e-12 ||
    max(abs(selected_time_from_raw -
      matrix(time_grid, nrow = nrow(beta_hat), ncol = n_time, byrow = TRUE))) >
      1e-12) {
  stop("The cached beta estimates, adjusted SEs, or time grid do not match raw FASH.")
}
rm(selected_beta_from_raw, selected_se_from_raw, selected_time_from_raw)

method_ids <- c("diagonal", "C1", "C2")
method_labels <- c("Diagonal SE", "C1: Screen 1", "C2: Screen 2")
names(method_labels) <- method_ids
raw_fits <- list()
raw_warnings <- list()
raw_elapsed <- numeric(length(method_ids))
names(raw_elapsed) <- method_ids

message("Refitting the diagonal empirical-Bayes mixture on the fixed subset.")
diagonal_start <- proc.time()[["elapsed"]]
diagonal_capture <- capture_warnings(
  refit_fash_from_cached_likelihood(
    full_raw,
    selected_indices,
    penalty = settings$penalty
  )
)
raw_elapsed["diagonal"] <- proc.time()[["elapsed"]] - diagonal_start
raw_fits$diagonal <- validate_fash_fit(
  diagonal_capture$value,
  pair_metadata$pair_key,
  psd_grid,
  "diagonal raw fit"
)
raw_warnings$diagonal <- diagonal_capture$warnings

validation_count <- min(validation_units, nrow(pair_metadata))
validation_rows <- seq_len(validation_count)
validation_beta <- beta_hat[validation_rows, , drop = FALSE]
validation_se <- adjusted_se[validation_rows, , drop = FALSE]
validation_data <- make_fash_data_list(validation_beta, time_grid)
validation_data_with_se <- Map(function(dataset, se) {
  dataset$SE <- as.numeric(se)
  dataset
}, validation_data, split(validation_se, row(validation_se)))
names(validation_data_with_se) <- rownames(validation_beta)

message(
  "Validating the identity-precision path on ", validation_count,
  " real-data units."
)
identity_validation_start <- proc.time()[["elapsed"]]
validation_diagonal_capture <- capture_warnings(fashr::fash(
  Y = "beta",
  smooth_var = "time",
  S = "SE",
  data_list = validation_data_with_se,
  grid = psd_grid,
  likelihood = settings$likelihood,
  num_basis = settings$num_basis,
  betaprec = settings$betaprec,
  order = settings$order,
  pred_step = settings$pred_step,
  penalty = settings$penalty,
  num_cores = num_cores,
  verbose = FALSE
))
validation_identity_capture <- capture_warnings(
  fit_fash_with_shared_correlation(
    beta_hat = validation_beta,
    adjusted_se = validation_se,
    time_grid = time_grid,
    correlation = diag(n_time),
    settings = settings,
    psd_grid = psd_grid,
    num_cores = num_cores,
    verbose = FALSE
  )
)
identity_validation <- validate_identity_path_equivalence(
  validation_diagonal_capture$value,
  validation_identity_capture$value$fit,
  likelihood_tolerance = 1e-7,
  prior_tolerance = 1e-6,
  lfdr_tolerance = 1e-6
)
identity_validation$n_units <- validation_count
identity_validation$pair_keys <- rownames(validation_beta)
identity_validation$elapsed_seconds <-
  proc.time()[["elapsed"]] - identity_validation_start
identity_validation$diagonal_warnings <- validation_diagonal_capture$warnings
identity_validation$identity_precision_warnings <-
  validation_identity_capture$warnings
rm(
  validation_data,
  validation_data_with_se,
  validation_diagonal_capture,
  validation_identity_capture,
  validation_beta,
  validation_se
)
gc(verbose = FALSE)

full_raw_settings <- settings
rm(full_raw)
gc(verbose = FALSE)

correlation_list <- list(C1 = C1, C2 = C2)
correlated_diagnostics <- list()
for (method_id in names(correlation_list)) {
  message("Fitting the full correlated likelihood for ", method_id, ".")
  fit_start <- proc.time()[["elapsed"]]
  fit_capture <- capture_warnings(fit_fash_with_shared_correlation(
    beta_hat = beta_hat,
    adjusted_se = adjusted_se,
    time_grid = time_grid,
    correlation = correlation_list[[method_id]],
    settings = settings,
    psd_grid = psd_grid,
    num_cores = num_cores,
    verbose = TRUE
  ))
  raw_elapsed[method_id] <- proc.time()[["elapsed"]] - fit_start
  raw_fits[[method_id]] <- validate_fash_fit(
    fit_capture$value$fit,
    pair_metadata$pair_key,
    psd_grid,
    paste0(method_id, " raw fit")
  )
  raw_warnings[[method_id]] <- fit_capture$warnings
  correlated_diagnostics[[method_id]] <- list(
    correlation_diagnostics = fit_capture$value$correlation_diagnostics,
    correlation_precision = fit_capture$value$correlation_precision
  )
  rm(fit_capture)
  gc(verbose = FALSE)
}

message("Applying checked BF updates independently to all three raw fits.")
bf_updates <- run_bf_updates_checked(
  raw_fits,
  method_labels = unname(method_labels),
  pair_keys = pair_metadata$pair_key
)
bf_fits <- bf_updates$successful_fits
for (method_id in names(bf_fits)) {
  bf_fits[[method_id]] <- validate_fash_fit(
    bf_fits[[method_id]],
    pair_metadata$pair_key,
    psd_grid,
    paste0(method_id, " BF-adjusted fit")
  )
}
bf_update_status <- bf_updates$status
bf_warnings <- bf_updates$warnings
bf_elapsed <- bf_updates$elapsed_seconds

raw_prior_comparison <- compare_prior_weight_fits(
  raw_fits,
  fit_stage = "Raw",
  method_labels = unname(method_labels)
)
bf_prior_comparison <- compare_prior_weight_fits(
  bf_fits,
  fit_stage = "BF-adjusted",
  method_labels = unname(method_labels[names(bf_fits)])
)
raw_lfdr_comparison <- compare_lfdr_fits(
  raw_fits,
  pair_metadata = pair_metadata,
  fit_stage = "Raw",
  method_labels = unname(method_labels),
  alpha = alpha,
  top_n = 100L
)
bf_lfdr_comparison <- compare_lfdr_fits(
  bf_fits,
  pair_metadata = pair_metadata,
  fit_stage = "BF-adjusted",
  method_labels = unname(method_labels[names(bf_fits)]),
  alpha = alpha,
  top_n = 100L
)

prior_weights <- rbind(
  raw_prior_comparison$prior_weights,
  bf_prior_comparison$prior_weights
)
prior_pairwise_metrics <- rbind(
  raw_prior_comparison$pairwise_metrics,
  bf_prior_comparison$pairwise_metrics
)
lfdr_long <- rbind(
  raw_lfdr_comparison$lfdr_long,
  bf_lfdr_comparison$lfdr_long
)
lfdr_pairwise_metrics <- rbind(
  raw_lfdr_comparison$pairwise_metrics,
  bf_lfdr_comparison$pairwise_metrics
)
unavailable_bf_ids <- bf_update_status$method_id[
  !bf_update_status$bf_update_available
]
unavailable_bf_summary <- do.call(rbind, lapply(unavailable_bf_ids, function(
    method_id) {
  data.frame(
    fit_stage = "BF-adjusted",
    method_id = method_id,
    method_label = unname(method_labels[method_id]),
    n_units = nrow(pair_metadata),
    mean_lfdr = NA_real_,
    median_lfdr = NA_real_,
    discovered_units = NA_integer_,
    discovered_genes = NA_integer_,
    alpha = alpha,
    stringsAsFactors = FALSE
  )
}))
discovery_summary <- rbind(
  raw_lfdr_comparison$discovery_summary,
  bf_lfdr_comparison$discovery_summary,
  unavailable_bf_summary
)
top_lfdr_discrepancies <- rbind(
  raw_lfdr_comparison$top_discrepancies,
  bf_lfdr_comparison$top_discrepancies
)
lfdr_wide <- pair_metadata[, c("pair_key", "gene_id", "variant_id")]
for (method_id in method_ids) {
  lfdr_wide[[paste0(method_id, "_raw")]] <-
    raw_lfdr_comparison$lfdr_wide[[method_id]]
  lfdr_wide[[paste0(method_id, "_bf")]] <- if (method_id %in% names(bf_fits)) {
    bf_lfdr_comparison$lfdr_wide[[method_id]]
  } else {
    NA_real_
  }
}

pi0_rows <- prior_weights$psd == 0
pi0_table <- prior_weights[pi0_rows, c(
  "fit_stage", "method_id", "prior_weight"
)]
names(pi0_table)[3L] <- "pi0"
method_stage_summary <- merge(
  discovery_summary,
  pi0_table,
  by = c("fit_stage", "method_id"),
  all.x = TRUE,
  sort = FALSE
)
method_stage_summary$elapsed_seconds <- ifelse(
  method_stage_summary$fit_stage == "Raw",
  raw_elapsed[method_stage_summary$method_id],
  bf_elapsed[method_stage_summary$method_id]
)
method_stage_summary$result_status <- ifelse(
  method_stage_summary$fit_stage == "Raw",
  "Available",
  bf_update_status$bf_update_status[match(
    method_stage_summary$method_id,
    bf_update_status$method_id
  )]
)
method_stage_summary$likelihood_source <- ifelse(
  method_stage_summary$method_id == "diagonal",
  "Selected cached diagonal-SE likelihood rows",
  "Recomputed full-precision likelihood"
)
method_stage_summary <- method_stage_summary[order(
  match(method_stage_summary$fit_stage, c("Raw", "BF-adjusted")),
  match(method_stage_summary$method_id, method_ids)
), ]
rownames(method_stage_summary) <- NULL

identity_validation_table <- data.frame(
  n_units = identity_validation$n_units,
  row_centered_likelihood_maximum_difference =
    identity_validation$row_centered_likelihood_maximum_difference,
  prior_weight_maximum_difference =
    identity_validation$prior_weight_maximum_difference,
  lfdr_maximum_difference = identity_validation$lfdr_maximum_difference,
  elapsed_seconds = identity_validation$elapsed_seconds,
  stringsAsFactors = FALSE
)

configuration <- list(
  analysis_id = output_id,
  experiment = paste(
    "Fixed one-random-variant-per-gene FASH likelihood sensitivity for",
    "diagonal, Screen-1 C1, and Screen-2 C2 error structures"
  ),
  thinning_seed = expected_thinning_seed,
  n_units = nrow(pair_metadata),
  n_genes = length(unique(pair_metadata$gene_id)),
  n_time = n_time,
  time_grid = time_grid,
  alpha = alpha,
  num_cores = num_cores,
  validation_units = validation_count,
  method_ids = method_ids,
  method_labels = method_labels,
  precision_formula = "diag(1 / se_j) %*% solve(C) %*% diag(1 / se_j)",
  correlation_sources = c(
    C1 = "Screen 1: max_t |z_j(t)| < 2",
    C2 = paste(
      "Screen 2: max_t |z_j(t)| < 2 and",
      "|sqrt(16) * mean_t z_j(t)| < 2"
    )
  ),
  raw_settings = full_raw_settings,
  psd_grid = psd_grid,
  raw_penalty = settings$penalty,
  bf_update = bf_updates$strategy,
  bf_update_status = bf_update_status,
  source_files = source_files,
  raw_warnings = raw_warnings,
  bf_warnings = bf_warnings,
  raw_elapsed_seconds = raw_elapsed,
  bf_elapsed_seconds = bf_elapsed,
  identity_validation = identity_validation,
  package_versions = c(
    fashr = as.character(utils::packageVersion("fashr")),
    workflowr = as.character(utils::packageVersion("workflowr")),
    ggplot2 = as.character(utils::packageVersion("ggplot2"))
  ),
  r_version = R.version.string,
  generated_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE)
)

analysis <- list(
  configuration = configuration,
  pair_metadata = pair_metadata,
  correlations = list(C1 = C1, C2 = C2),
  correlation_diagnostics = correlation_diagnostics,
  correlation_matrices_long = correlation_matrices_long,
  identity_validation_table = identity_validation_table,
  bf_update_status = bf_update_status,
  method_stage_summary = method_stage_summary,
  prior_weights = prior_weights,
  prior_pairwise_metrics = prior_pairwise_metrics,
  lfdr_wide = lfdr_wide,
  lfdr_long = lfdr_long,
  lfdr_pairwise_metrics = lfdr_pairwise_metrics,
  discovery_summary = discovery_summary,
  top_lfdr_discrepancies = top_lfdr_discrepancies,
  correlated_diagnostics = correlated_diagnostics
)
fit_bundle <- list(
  configuration = configuration,
  pair_metadata = pair_metadata,
  correlations = list(C1 = C1, C2 = C2),
  raw_fits = raw_fits,
  bf_adjusted_fits = bf_fits,
  bf_update_status = bf_update_status
)

saveRDS(configuration, file.path(staging_directory, "configuration.rds"))
saveRDS(analysis, file.path(staging_directory, "analysis.rds"))
saveRDS(fit_bundle, file.path(staging_directory, "fit_bundle.rds"))
write_csv(pair_metadata, file.path(summary_directory, "pair_metadata.csv"))
write_csv(
  correlation_diagnostics,
  file.path(summary_directory, "correlation_diagnostics.csv")
)
write_csv(
  correlation_matrices_long,
  file.path(summary_directory, "correlation_matrices_long.csv")
)
write_csv(
  identity_validation_table,
  file.path(summary_directory, "identity_path_validation.csv")
)
write_csv(
  bf_update_status,
  file.path(summary_directory, "bf_update_status.csv")
)
write_csv(
  method_stage_summary,
  file.path(summary_directory, "method_stage_summary.csv")
)
write_csv(prior_weights, file.path(summary_directory, "prior_weights.csv"))
write_csv(
  prior_pairwise_metrics,
  file.path(summary_directory, "prior_pairwise_metrics.csv")
)
write_csv(lfdr_wide, file.path(summary_directory, "unit_lfdr_wide.csv"))
write_csv(lfdr_long, file.path(summary_directory, "unit_lfdr_long.csv"))
write_csv(
  lfdr_pairwise_metrics,
  file.path(summary_directory, "lfdr_pairwise_metrics.csv")
)
write_csv(
  discovery_summary,
  file.path(summary_directory, "discovery_summary.csv")
)
write_csv(
  top_lfdr_discrepancies,
  file.path(summary_directory, "top_lfdr_discrepancies.csv")
)

if (!file.rename(staging_directory, final_directory)) {
  stop("Could not atomically finalize the sensitivity cache.")
}

cat("\nCorrelated-likelihood sensitivity completed.\n")
print(method_stage_summary[, c(
  "fit_stage", "method_label", "pi0", "discovered_units",
  "discovered_genes", "elapsed_seconds"
)])
print(lfdr_pairwise_metrics[, c(
  "fit_stage", "reference_method_label", "comparison_method_label",
  "spearman_lfdr", "mean_absolute_lfdr_difference",
  "maximum_absolute_lfdr_difference", "discovery_jaccard"
)])
