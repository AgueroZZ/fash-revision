#!/usr/bin/env Rscript

# Build retained caches for the internal FASH versus FASH-CL enrichment page.

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

file_provenance <- function(path) {
  information <- file.info(path)
  data.frame(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    byte_size = unname(information$size),
    md5 = unname(tools::md5sum(path)),
    modified_at = format(information$mtime, tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}

message_step <- function(index, total, text) {
  message("[", index, "/", total, "] ", text)
}

load_raw_fit_components <- function(path, expected_object) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  if (!identical(loaded, expected_object)) {
    stop("Unexpected object in ", path, ".")
  }
  fit <- environment[[expected_object]]
  output <- list(
    pair_keys = names(fit$fash_data$data_list),
    log_likelihood_matrix = fit$L_matrix,
    psd_grid = fit$psd_grid
  )
  rm(fit, environment)
  invisible(gc())
  output
}

workflowr_root <- find_workflowr_root()
analysis_directory <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "fash_cl_variant_enrichment_comparison"
)
source(file.path(analysis_directory, "fash_cl_variant_enrichment_helpers.R"))
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "variant_annotation_enrichment",
  "variant_annotation_enrichment_helpers.R"
))

current_raw_path <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real", "fash_fit1_all.RData"
)
current_adjusted_path <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
)
fash_cl_raw_path <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real", "fash_fit1_all_CL.RData"
)
matching_covariate_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "variant_annotation_enrichment",
  "variant_matching_covariates.rds"
)
custom_annotation_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "variant_annotation_enrichment",
  "variant_annotation_matrix.rds"
)
baseline_annotation_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "baseline_ld_variant_enrichment",
  "baseline_ld_binary_annotation_matrix.rds"
)
output_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fash_cl_variant_enrichment_comparison"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

input_paths <- c(
  current_raw_path,
  current_adjusted_path,
  fash_cl_raw_path,
  matching_covariate_path,
  custom_annotation_path,
  baseline_annotation_path,
  file.path(analysis_directory, "fash_cl_variant_enrichment_helpers.R"),
  file.path(
    analysis_directory,
    "run_fash_cl_variant_enrichment_comparison.R"
  ),
  file.path(
    workflowr_root,
    "code",
    "revision_simulations",
    "internal",
    "variant_annotation_enrichment",
    "variant_annotation_enrichment_helpers.R"
  )
)
if (any(!file.exists(input_paths))) {
  stop("At least one retained input is missing.")
}

analysis_started <- proc.time()[["elapsed"]]
alpha <- 0.05
matching_seeds <- seq.int(20260807L, length.out = 100L)
controls_per_variant <- 5L

message_step(1, 8, "Validating memory-bounded BF adjustment on current FASH.")
current_raw <- load_raw_fit_components(current_raw_path, "fash_fit1")
current_validation_started <- proc.time()[["elapsed"]]
current_recomputed <- compute_bf_adjusted_lfdr(
  current_raw$log_likelihood_matrix,
  chunk_size = 50000L,
  verbose = FALSE
)
current_validation_elapsed <-
  proc.time()[["elapsed"]] - current_validation_started
current_pair_keys <- current_raw$pair_keys
rm(current_raw)
invisible(gc())

adjusted_environment <- new.env(parent = emptyenv())
loaded <- load(current_adjusted_path, envir = adjusted_environment)
if (!identical(loaded, "fash_fit1_update")) {
  stop("The retained current adjusted fit has an unexpected object name.")
}
current_reference <- adjusted_environment$fash_fit1_update
current_reference_keys <- names(current_reference$fash_data$data_list)
current_reference_lfdr <- as.numeric(current_reference$lfdr)
current_reference_null_weight <- current_reference$prior_weights$prior_weight[
  current_reference$prior_weights$psd == 0
]
current_reference_alternative_weights <- numeric(
  length(current_reference$psd_grid) - 1L
)
reference_nonnull <- current_reference$prior_weights[
  current_reference$prior_weights$psd != 0,
  ,
  drop = FALSE
]
reference_weight_indices <- match(
  reference_nonnull$psd,
  current_reference$psd_grid[-1L]
)
if (anyNA(reference_weight_indices)) {
  stop("Retained current alternative weights do not align to the PSD grid.")
}
current_reference_alternative_weights[reference_weight_indices] <-
  reference_nonnull$prior_weight / (1 - current_reference_null_weight)
current_reference_log_bf <- compute_log_bayes_factor(
  current_reference$L_matrix,
  current_reference_alternative_weights,
  chunk_size = 50000L
)
current_reference_reconstructed_lfdr <- stats::plogis(
  stats::qlogis(current_reference_null_weight) - current_reference_log_bf
)
rm(current_reference, adjusted_environment)
invisible(gc())

current_recomputed_calls <- select_cumulative_lfdr_calls(
  current_recomputed$lfdr,
  alpha
)
current_reference_calls <- select_cumulative_lfdr_calls(
  current_reference_lfdr,
  alpha
)
current_call_intersection <- length(intersect(
  current_recomputed_calls,
  current_reference_calls
))
current_call_union <- length(union(
  current_recomputed_calls,
  current_reference_calls
))
bf_validation <- data.frame(
  pair_keys_identical = identical(current_pair_keys, current_reference_keys),
  stored_lfdr_reconstruction_maximum_absolute_difference = max(abs(
    current_reference_reconstructed_lfdr - current_reference_lfdr
  )),
  null_weight_recomputed = current_recomputed$null_weight,
  null_weight_reference = current_reference_null_weight,
  null_weight_absolute_difference = abs(
    current_recomputed$null_weight - current_reference_null_weight
  ),
  maximum_lfdr_absolute_difference = max(abs(
    current_recomputed$lfdr - current_reference_lfdr
  )),
  mean_lfdr_absolute_difference = mean(abs(
    current_recomputed$lfdr - current_reference_lfdr
  )),
  lfdr_correlation = stats::cor(
    current_recomputed$lfdr,
    current_reference_lfdr
  ),
  recomputed_discovery_call_count = length(current_recomputed_calls),
  reference_discovery_call_count = length(current_reference_calls),
  discovery_call_intersection = current_call_intersection,
  discovery_call_union = current_call_union,
  discovery_call_jaccard = current_call_intersection / current_call_union,
  discovery_calls_identical = identical(
    sort(current_recomputed_calls),
    sort(current_reference_calls)
  ),
  stringsAsFactors = FALSE
)
message(
  "  BF validation: stored-lfdr reconstruction difference = ",
  format(
    bf_validation$stored_lfdr_reconstruction_maximum_absolute_difference,
    scientific = TRUE
  ),
  ", null-weight rerun difference = ",
  format(bf_validation$null_weight_absolute_difference, scientific = TRUE),
  ", mean lfdr rerun difference = ",
  format(bf_validation$mean_lfdr_absolute_difference, scientific = TRUE),
  ", call Jaccard = ",
  format(bf_validation$discovery_call_jaccard, digits = 6)
)
if (!bf_validation$pair_keys_identical ||
    bf_validation$stored_lfdr_reconstruction_maximum_absolute_difference >=
      1e-12 ||
    bf_validation$null_weight_absolute_difference >= 1e-7 ||
    bf_validation$mean_lfdr_absolute_difference >= 1e-3 ||
    bf_validation$lfdr_correlation <= 0.999) {
  stop("The memory-bounded BF adjustment failed retained-current validation.")
}

current_sets <- derive_all_and_lead_sets(
  current_reference_keys,
  current_reference_lfdr,
  alpha = alpha
)
if (!identical(
  unname(unlist(current_sets$summary[, c(
    "pair_count", "unique_variant_count", "unique_gene_count"
  )])),
  c(9205L, 1177L, 9139L, 1170L, 1177L, 1177L)
)) {
  stop("Current FASH discovery counts no longer match retained invariants.")
}
rm(
  current_reference_lfdr,
  current_reference_alternative_weights,
  current_reference_log_bf,
  current_reference_reconstructed_lfdr,
  current_recomputed,
  current_recomputed_calls,
  current_reference_calls
)
invisible(gc())

message_step(2, 8, "Applying the validated BF adjustment to FASH-CL.")
fash_cl_raw <- load_raw_fit_components(fash_cl_raw_path, "fash_fit1")
if (!identical(current_pair_keys, fash_cl_raw$pair_keys)) {
  stop("Current FASH and FASH-CL pair keys are not identically aligned.")
}
fash_cl_started <- proc.time()[["elapsed"]]
fash_cl_adjustment <- compute_bf_adjusted_lfdr(
  fash_cl_raw$log_likelihood_matrix,
  chunk_size = 50000L,
  verbose = FALSE
)
fash_cl_adjustment_elapsed <- proc.time()[["elapsed"]] - fash_cl_started
fash_cl_sets <- derive_all_and_lead_sets(
  fash_cl_raw$pair_keys,
  fash_cl_adjustment$lfdr,
  alpha = alpha
)
rm(fash_cl_raw)
invisible(gc())

current_method_sets <- list(
  all = current_sets$all_variants,
  one_lead_per_gene = current_sets$lead_variants
)
fash_cl_method_sets <- list(
  all = fash_cl_sets$all_variants,
  one_lead_per_gene = fash_cl_sets$lead_variants
)
method_overlap <- summarize_method_overlap(
  current_method_sets,
  fash_cl_method_sets
)
discovery_set_summary <- rbind(
  transform(current_sets$summary, method = "Current FASH"),
  transform(fash_cl_sets$summary, method = "FASH-CL")
)
discovery_set_summary <- discovery_set_summary[, c(
  "method", "selection_strategy", "pair_count", "unique_variant_count",
  "unique_gene_count"
)]

selected_sets <- list(
  current_all = current_sets$all_variants,
  current_one_lead = current_sets$lead_variants,
  fash_cl_all = fash_cl_sets$all_variants,
  fash_cl_one_lead = fash_cl_sets$lead_variants
)
set_metadata <- data.frame(
  discovery_set = names(selected_sets),
  method = c("Current FASH", "Current FASH", "FASH-CL", "FASH-CL"),
  selection_strategy = c(
    "All discoveries", "One lead variant per gene",
    "All discoveries", "One lead variant per gene"
  ),
  original_variant_count = lengths(selected_sets),
  stringsAsFactors = FALSE
)

message_step(3, 8, "Loading and validating retained annotations and covariates.")
variant_table <- readRDS(matching_covariate_path)
custom_matrix <- readRDS(custom_annotation_path)
baseline_matrix <- readRDS(baseline_annotation_path)
if (nrow(variant_table) != 745867L || anyDuplicated(variant_table$variant_id) ||
    !identical(variant_table$variant_id, custom_matrix$variant_id) ||
    anyDuplicated(baseline_matrix$variant_id)) {
  stop("Retained matching or annotation inputs are not aligned.")
}
for (matrix in list(custom_matrix, baseline_matrix)) {
  annotation_columns <- setdiff(names(matrix), c("variant_id", "chromosome"))
  if (!all(vapply(matrix[annotation_columns], is.logical, logical(1)))) {
    stop("Retained annotation columns must be logical.")
  }
}
custom_enhancers <- c(
  "ENCODE cCRE enhancer-like",
  "Roadmap E020 iPS-20b: Enhancer",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer",
  "Roadmap E095 left ventricle: Enhancer"
)
baseline_enhancers <- c(
  "Enhancer_Andersson",
  "Enhancer_Hoffman",
  "WeakEnhancer_Hoffman",
  "SuperEnhancer_Hnisz",
  "Human_Enhancer_Villar"
)
if (!all(custom_enhancers %in% names(custom_matrix)) ||
    !all(baseline_enhancers %in% names(baseline_matrix))) {
  stop("A pre-specified enhancer annotation is missing.")
}

build_coverage_rows <- function(annotation_system, covered_ids) {
  do.call(rbind, lapply(names(selected_sets), function(set_name) {
    ids <- unique(selected_sets[[set_name]])
    data.frame(
      annotation_system = annotation_system,
      discovery_set = set_name,
      original_variant_count = length(ids),
      covered_variant_count = sum(ids %in% covered_ids),
      coverage_proportion = mean(ids %in% covered_ids),
      stringsAsFactors = FALSE
    )
  }))
}
coverage_by_annotation_system <- rbind(
  build_coverage_rows("Custom regulatory", custom_matrix$variant_id),
  build_coverage_rows("baselineLD v2.2", baseline_matrix$variant_id)
)

add_result_metadata <- function(analysis, annotation_system) {
  results <- merge(
    analysis$results,
    set_metadata,
    by = "discovery_set",
    all.x = TRUE,
    sort = FALSE
  )
  results$annotation_system <- annotation_system
  results$q_value_within_set <- ave(
    results$p_value,
    results$discovery_set,
    FUN = function(values) stats::p.adjust(values, method = "BH")
  )
  results
}

message_step(4, 8, "Running matched enrichment for custom regulatory annotations.")
custom_sets <- lapply(selected_sets, function(ids) {
  intersect(ids, custom_matrix$variant_id)
})
custom_started <- proc.time()[["elapsed"]]
custom_analysis <- run_streaming_matched_enrichment(
  variant_table = variant_table,
  annotation_matrix = custom_matrix,
  selected_sets = custom_sets,
  seeds = matching_seeds,
  controls_per_variant = controls_per_variant,
  chromosomes = as.character(1:22),
  minimum_overlap_count = 10L,
  verbose = TRUE
)
custom_elapsed <- proc.time()[["elapsed"]] - custom_started
custom_results <- add_result_metadata(custom_analysis, "Custom regulatory")
custom_balance <- custom_analysis$matching_balance
custom_balance$annotation_system <- "Custom regulatory"
custom_relaxation <- custom_analysis$relaxation_audit
custom_relaxation$annotation_system <- "Custom regulatory"

message_step(5, 8, "Running matched enrichment for baselineLD v2.2.")
baseline_variant_table <- variant_table[
  match(baseline_matrix$variant_id, variant_table$variant_id),
  ,
  drop = FALSE
]
if (anyNA(baseline_variant_table$variant_id) ||
    !identical(baseline_variant_table$variant_id, baseline_matrix$variant_id)) {
  stop("baselineLD annotations and matching covariates are not aligned.")
}
baseline_sets <- lapply(selected_sets, function(ids) {
  intersect(ids, baseline_matrix$variant_id)
})
baseline_started <- proc.time()[["elapsed"]]
baseline_analysis <- run_streaming_matched_enrichment(
  variant_table = baseline_variant_table,
  annotation_matrix = baseline_matrix,
  selected_sets = baseline_sets,
  seeds = matching_seeds,
  controls_per_variant = controls_per_variant,
  chromosomes = as.character(1:22),
  minimum_overlap_count = 10L,
  verbose = TRUE
)
baseline_elapsed <- proc.time()[["elapsed"]] - baseline_started
baseline_results <- add_result_metadata(baseline_analysis, "baselineLD v2.2")
baseline_balance <- baseline_analysis$matching_balance
baseline_balance$annotation_system <- "baselineLD v2.2"
baseline_relaxation <- baseline_analysis$relaxation_audit
baseline_relaxation$annotation_system <- "baselineLD v2.2"

message_step(6, 8, "Summarizing the pre-specified enhancer panels.")
enrichment_results <- rbind(custom_results, baseline_results)
enhancer_results <- rbind(
  custom_results[custom_results$annotation %in% custom_enhancers, ],
  baseline_results[baseline_results$annotation %in% baseline_enhancers, ]
)
enhancer_summary <- summarize_enhancer_panel(enhancer_results)
matching_balance <- rbind(custom_balance, baseline_balance)
matching_relaxation <- rbind(custom_relaxation, baseline_relaxation)
maximum_absolute_smd <- max(
  abs(matching_balance$standardized_mean_difference),
  na.rm = TRUE
)
if (!is.finite(maximum_absolute_smd) || maximum_absolute_smd >= 0.10) {
  stop("Matched controls failed the balance threshold.")
}

message_step(7, 8, "Building provenance and retained cache metadata.")
total_elapsed <- proc.time()[["elapsed"]] - analysis_started
runtime_summary <- data.frame(
  stage = c(
    "current_bf_validation",
    "fash_cl_bf_adjustment",
    "custom_matched_enrichment",
    "baseline_ld_matched_enrichment",
    "analysis_through_cache_assembly"
  ),
  elapsed_seconds = c(
    current_validation_elapsed,
    fash_cl_adjustment_elapsed,
    custom_elapsed,
    baseline_elapsed,
    total_elapsed
  ),
  stringsAsFactors = FALSE
)
input_provenance <- do.call(rbind, lapply(input_paths, file_provenance))
configuration <- list(
  analysis_id = "revision_internal_fash_cl_variant_enrichment_comparison",
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  current_pc_definition = "Five time-point-specific PCs",
  fash_cl_pc_definition = "Five cell-line-collapsed PCs",
  fash_order = 1L,
  fdr_threshold = alpha,
  controls_per_variant = controls_per_variant,
  matching_seed_count = length(matching_seeds),
  matching_seeds = matching_seeds,
  jackknife_blocks = as.character(1:22),
  annotation_systems = c("Custom regulatory", "baselineLD v2.2"),
  custom_enhancers = custom_enhancers,
  baseline_enhancers = baseline_enhancers,
  comparison_scope = paste(
    "Descriptive comparison of enrichment point estimates; discovery-set",
    "sizes and membership differ, so this is not a formal paired method test."
  )
)
cache <- list(
  configuration = configuration,
  input_provenance = input_provenance,
  bf_validation = bf_validation,
  discovery_set_summary = discovery_set_summary,
  method_overlap = method_overlap,
  set_metadata = set_metadata,
  coverage_by_annotation_system = coverage_by_annotation_system,
  enrichment_results = enrichment_results,
  enhancer_results = enhancer_results,
  enhancer_summary = enhancer_summary,
  matching_balance = matching_balance,
  matching_relaxation = matching_relaxation,
  runtime_summary = runtime_summary,
  maximum_absolute_smd = maximum_absolute_smd,
  session_info = utils::sessionInfo()
)

message_step(8, 8, "Writing reproducibility artifacts.")
saveRDS(cache, file.path(output_directory, "analysis_cache.rds"),
        compress = "gzip")
saveRDS(
  list(
    pair_count = length(current_pair_keys),
    null_weight = fash_cl_adjustment$null_weight,
    alternative_weights = fash_cl_adjustment$alternative_weights,
    lfdr = fash_cl_adjustment$lfdr,
    log_bayes_factor = fash_cl_adjustment$log_bayes_factor,
    bf_mixture_status = fash_cl_adjustment$bf_mixture_status,
    posterior_mixture_status = fash_cl_adjustment$posterior_mixture_status,
    elapsed_seconds = fash_cl_adjustment$elapsed_seconds,
    current_validation = bf_validation
  ),
  file.path(output_directory, "fash_cl_bf_adjustment.rds"),
  compress = "gzip"
)
saveRDS(
  selected_sets,
  file.path(output_directory, "discovery_sets.rds"),
  compress = "gzip"
)
write_output <- function(value, filename) {
  utils::write.csv(
    value,
    file.path(output_directory, filename),
    row.names = FALSE,
    quote = TRUE
  )
}
write_output(discovery_set_summary, "discovery_set_summary.csv")
write_output(method_overlap, "method_overlap.csv")
write_output(coverage_by_annotation_system,
             "coverage_by_annotation_system.csv")
write_output(enrichment_results, "enrichment_results.csv")
write_output(enhancer_results, "enhancer_results.csv")
write_output(enhancer_summary, "enhancer_summary.csv")
write_output(matching_balance, "matching_balance.csv")
write_output(runtime_summary, "runtime_summary.csv")
write_output(input_provenance, "input_provenance.csv")
message(
  "Completed FASH versus FASH-CL enrichment comparison in ",
  round(total_elapsed, 1),
  " seconds: ",
  file.path(output_directory, "analysis_cache.rds")
)
