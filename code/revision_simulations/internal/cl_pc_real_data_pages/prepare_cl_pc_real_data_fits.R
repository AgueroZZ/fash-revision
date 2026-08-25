#!/usr/bin/env Rscript

# Prepare corrected full CL-PC FASH fits for the dynamic and nonlinear pages.

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

find_workflowr_root <- function(start = getwd()) {
  nested_root <- file.path(start, "coderepo-local")
  if (file.exists(file.path(nested_root, "_workflowr.yml"))) {
    return(normalizePath(nested_root, winslash = "/", mustWork = TRUE))
  }
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

sha256_file <- function(path) {
  output <- system2(
    "shasum",
    args = c("-a", "256", normalizePath(path, mustWork = TRUE)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || length(output) != 1L) {
    stop("Could not compute SHA-256 for ", path, ".")
  }
  hash <- strsplit(output, "[[:space:]]+")[[1L]][1L]
  if (!grepl("^[0-9a-f]{64}$", hash)) {
    stop("Unexpected SHA-256 output for ", path, ".")
  }
  hash
}

file_provenance <- function(label, path) {
  information <- file.info(path)
  data.frame(
    label = label,
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    byte_size = unname(information$size),
    md5 = unname(tools::md5sum(path)),
    sha256 = sha256_file(path),
    modified_at = format(information$mtime, tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}

package_provenance <- function(package) {
  description <- utils::packageDescription(package)
  remote_sha <- description[["RemoteSha"]]
  if (is.null(remote_sha) || length(remote_sha) != 1L || is.na(remote_sha)) {
    remote_sha <- ""
  }
  data.frame(
    package = package,
    version = as.character(description[["Version"]]),
    remote_sha = as.character(remote_sha),
    library_path = normalizePath(
      find.package(package),
      winslash = "/",
      mustWork = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

load_exact_object <- function(path, expected_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  if (!identical(loaded, expected_name)) {
    stop("Unexpected object in ", path, ": ", paste(loaded, collapse = ", "))
  }
  environment[[expected_name]]
}

select_cumulative_lfdr <- function(lfdr, alpha = 0.05) {
  lfdr <- as.numeric(lfdr)
  if (!length(lfdr) || any(!is.finite(lfdr)) || any(lfdr < 0 | lfdr > 1)) {
    stop("lfdr must contain finite probabilities.")
  }
  ordering <- order(lfdr, seq_along(lfdr), method = "radix")
  accepted <- which(cumsum(lfdr[ordering]) / seq_along(ordering) <= alpha)
  if (!length(accepted)) integer() else ordering[seq_len(max(accepted))]
}

validate_raw_fit <- function(fit, order, pair_keys, expected_grid) {
  required_settings <- c(
    "num_basis", "betaprec", "order", "pred_step", "likelihood", "penalty"
  )
  if (!inherits(fit, "fash") ||
      !all(required_settings %in% names(fit$settings)) ||
      !identical(names(fit$fash_data$data_list), pair_keys) ||
      !identical(as.numeric(fit$psd_grid), expected_grid) ||
      as.integer(fit$settings$order) != as.integer(order) ||
      as.numeric(fit$settings$betaprec) != 0 ||
      as.numeric(fit$settings$pred_step) != 1 ||
      as.integer(fit$settings$penalty) != 10L ||
      !identical(as.character(fit$settings$likelihood), "gaussian")) {
    stop("A raw CL-PC fit failed the matched settings contract.")
  }
  invisible(TRUE)
}

validate_adjusted_fit <- function(fit,
                                  compact,
                                  expected_method,
                                  expected_pairs,
                                  expected_genes,
                                  pair_keys,
                                  alpha = 0.05) {
  prior <- fit$prior_weights[, c("psd", "prior_weight"), drop = FALSE]
  compact_prior <- compact$prior_weights[, c("psd", "prior_weight"), drop = FALSE]
  selected <- select_cumulative_lfdr(fit$lfdr, alpha)
  genes <- sub("_.*$", "", pair_keys[selected])
  if (!identical(compact$method, expected_method) ||
      !identical(as.numeric(fit$lfdr), as.numeric(compact$lfdr)) ||
      !identical(prior, compact_prior) ||
      nrow(fit$posterior_weights) != length(pair_keys) ||
      length(selected) != expected_pairs ||
      length(unique(genes)) != expected_genes) {
    stop("A corrected full CL-PC fit does not match the compact reference.")
  }
  list(
    selected_indices = selected,
    pair_count = length(selected),
    gene_count = length(unique(genes)),
    variant_count = length(unique(sub("^[^_]+_", "", pair_keys[selected]))),
    estimated_pi0 = compact$estimated_pi0
  )
}

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required.")
}
suppressPackageStartupMessages(library(fashr))

workflowr_root <- find_workflowr_root()
raw_directory <- file.path(workflowr_root, "output", "dynamic_eQTL_real")
matched_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "matched_fash_linear_real_data_cl_fashr_0_1_43"
)
statistics_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fash_linear_real_data_ablation_cl",
  "sufficient_statistics.rds"
)
raw_paths <- c(
  iwp1 = file.path(raw_directory, "fash_fit1_all_CL.RData"),
  iwp2 = file.path(raw_directory, "fash_fit2_all_CL.RData")
)
compact_paths <- c(
  iwp1 = file.path(matched_directory, "cl_iwp1_bf_adjustment.rds"),
  iwp2 = file.path(matched_directory, "cl_iwp2_bf_adjustment.rds")
)
output_parent <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal"
)
output_id <- "cl_pc_real_data_pages_fashr0143"
final_output_directory <- file.path(output_parent, output_id)
staging_directory <- file.path(
  output_parent,
  paste0(".", output_id, "_staging_", Sys.getpid())
)
required_paths <- c(raw_paths, compact_paths, statistics_path)
if (any(!file.exists(required_paths))) {
  stop("At least one required CL-PC input is missing.")
}
if (file.exists(final_output_directory)) {
  stop("Refusing to overwrite the retained CL-PC page cache.")
}
if (file.exists(staging_directory)) {
  stop("Unexpected staging-directory collision.")
}

expected_package <- data.frame(
  package = "fashr",
  version = "0.1.43",
  remote_sha = "bf223df75da6e41ae48607a56b4cd12d7c3b24e7",
  stringsAsFactors = FALSE
)
observed_package <- package_provenance("fashr")
if (!identical(
  observed_package[, names(expected_package), drop = FALSE],
  expected_package
)) {
  stop("The approved fashr 0.1.43 build is not installed.")
}
expected_raw_md5 <- c(
  iwp1 = "57383f6068304370930f2410021808cb",
  iwp2 = "550fd3d6bb5cf7a5a6d80288da7974a2"
)
if (!identical(unname(tools::md5sum(raw_paths)), unname(expected_raw_md5))) {
  stop("A raw CL-PC fit changed from the approved input.")
}

statistics_cache <- readRDS(statistics_path)
pair_keys <- as.character(statistics_cache$statistics$unit_id)
expected_grid <- sort(c(0, exp(-0.5 * seq(0, 10, by = 0.2))))
if (length(pair_keys) != 1009173L || anyDuplicated(pair_keys)) {
  stop("The retained CL-PC pair-key cache is invalid.")
}
rm(statistics_cache)
invisible(gc())

input_provenance <- do.call(rbind, Map(
  file_provenance,
  label = c("cl_iwp1_raw", "cl_iwp2_raw", "cl_iwp1_compact_bf",
            "cl_iwp2_compact_bf", "cl_pair_statistics"),
  path = c(raw_paths, compact_paths, statistics_path)
))

dir.create(staging_directory, recursive = FALSE)
completed <- FALSE
on.exit({
  if (!completed && dir.exists(staging_directory)) {
    unlink(staging_directory, recursive = TRUE, force = TRUE)
  }
}, add = TRUE)

analysis_started <- proc.time()[["elapsed"]]
fit_results <- list()
prior_tables <- list()
runtime_rows <- list()
validation_rows <- list()

for (order_label in names(raw_paths)) {
  order <- if (identical(order_label, "iwp1")) 1L else 2L
  expected_method <- if (order == 1L) "FASH-IWP1" else "FASH-IWP2"
  expected_pairs <- if (order == 1L) 5395L else 60L
  expected_genes <- if (order == 1L) 686L else 6L
  raw_object <- if (order == 1L) "fash_fit1" else "fash_fit2"
  adjusted_object <- if (order == 1L) "fash_fit1_update" else "fash_fit2_update"
  adjusted_file <- if (order == 1L) {
    "fash_fit1_update_CL.RData"
  } else {
    "fash_fit2_update_CL.RData"
  }

  message("Loading and validating ", order_label, " raw CL-PC fit.")
  raw_fit <- load_exact_object(raw_paths[[order_label]], raw_object)
  validate_raw_fit(raw_fit, order, pair_keys, expected_grid)
  compact <- readRDS(compact_paths[[order_label]])

  update_started <- proc.time()[["elapsed"]]
  adjusted_fit <- fashr::BF_update(raw_fit, plot = FALSE)
  update_elapsed <- proc.time()[["elapsed"]] - update_started
  validation <- validate_adjusted_fit(
    adjusted_fit,
    compact,
    expected_method,
    expected_pairs,
    expected_genes,
    pair_keys
  )

  output_environment <- new.env(parent = emptyenv())
  assign(adjusted_object, adjusted_fit, envir = output_environment)
  save(
    list = adjusted_object,
    file = file.path(staging_directory, adjusted_file),
    envir = output_environment,
    compress = "gzip"
  )

  fit_results[[order_label]] <- data.frame(
    method = expected_method,
    adjustment = "BF-adjusted",
    pair_count = validation$pair_count,
    gene_count = validation$gene_count,
    variant_count = validation$variant_count,
    estimated_pi0 = validation$estimated_pi0,
    alpha = 0.05,
    stringsAsFactors = FALSE
  )
  prior <- adjusted_fit$prior_weights
  prior$method <- expected_method
  prior_tables[[order_label]] <- prior[, c("method", "psd", "prior_weight")]
  runtime_rows[[order_label]] <- data.frame(
    stage = paste(expected_method, "BF update"),
    elapsed_seconds = update_elapsed,
    stringsAsFactors = FALSE
  )
  validation_rows[[order_label]] <- data.frame(
    check = c(
      paste0(order_label, "_settings_and_pair_order"),
      paste0(order_label, "_full_fit_matches_compact_reference"),
      paste0(order_label, "_discovery_counts")
    ),
    passed = TRUE,
    stringsAsFactors = FALSE
  )

  rm(raw_fit, adjusted_fit, compact, output_environment)
  invisible(gc())
}

discovery_counts <- do.call(rbind, fit_results)
prior_weights <- do.call(rbind, prior_tables)
runtime_summary <- do.call(rbind, runtime_rows)
runtime_summary <- rbind(
  runtime_summary,
  data.frame(
    stage = "Complete fit preparation",
    elapsed_seconds = proc.time()[["elapsed"]] - analysis_started,
    stringsAsFactors = FALSE
  )
)
validation <- rbind(
  data.frame(
    check = c(
      "approved_fashr_version_and_sha",
      "raw_cl_fit_md5_values",
      "pair_keys_and_grid"
    ),
    passed = TRUE,
    stringsAsFactors = FALSE
  ),
  do.call(rbind, validation_rows)
)

adjusted_paths <- c(
  iwp1 = file.path(staging_directory, "fash_fit1_update_CL.RData"),
  iwp2 = file.path(staging_directory, "fash_fit2_update_CL.RData")
)
adjusted_provenance <- do.call(rbind, Map(
  file_provenance,
  label = c("cl_iwp1_adjusted_full", "cl_iwp2_adjusted_full"),
  path = adjusted_paths
))
adjusted_provenance$path <- file.path(
  normalizePath(final_output_directory, winslash = "/", mustWork = FALSE),
  basename(adjusted_provenance$path)
)

configuration <- list(
  analysis_id = output_id,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  package_version = observed_package$version,
  package_sha = observed_package$remote_sha,
  pc_correction = "Five cell-line-collapsed PCs repeated across time",
  pair_count = length(pair_keys),
  time_grid = 0:15,
  evaluation_grid = seq(0, 15, by = 0.1),
  posterior_draws = 3000L,
  classification_seed = 20260820L,
  alpha = 0.05,
  switch_threshold = 0.25,
  category_definitions = c(
    early = "t <= 3 versus t > 3",
    middle = "3 < t < 12 versus t <= 3 or t >= 12",
    late = "t >= 12 versus t < 12",
    switch = "signed excursion on both sides of zero exceeds 0.25"
  )
)
analysis_cache <- list(
  configuration = configuration,
  package_provenance = observed_package,
  input_provenance = input_provenance,
  adjusted_provenance = adjusted_provenance,
  discovery_counts = discovery_counts,
  prior_weights = prior_weights,
  validation = validation,
  runtime_summary = runtime_summary,
  session_info = utils::sessionInfo()
)
saveRDS(
  analysis_cache,
  file.path(staging_directory, "analysis_cache.rds"),
  compress = "gzip"
)
utils::write.csv(
  input_provenance,
  file.path(staging_directory, "input_provenance.csv"),
  row.names = FALSE
)
utils::write.csv(
  adjusted_provenance,
  file.path(staging_directory, "adjusted_fit_provenance.csv"),
  row.names = FALSE
)
utils::write.csv(
  discovery_counts,
  file.path(staging_directory, "discovery_counts.csv"),
  row.names = FALSE
)
utils::write.csv(
  prior_weights,
  file.path(staging_directory, "prior_weights.csv"),
  row.names = FALSE
)
utils::write.csv(
  validation,
  file.path(staging_directory, "validation.csv"),
  row.names = FALSE
)
utils::write.csv(
  runtime_summary,
  file.path(staging_directory, "runtime_summary.csv"),
  row.names = FALSE
)

if (!file.rename(staging_directory, final_output_directory)) {
  stop("Could not finalize the retained CL-PC page cache.")
}
completed <- TRUE
message("Completed corrected CL-PC fit preparation: ", final_output_directory)
print(discovery_counts, row.names = FALSE)
