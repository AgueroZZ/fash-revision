#!/usr/bin/env Rscript

# Rebuild BF-adjusted summaries from completed raw fits without recomputing likelihoods.

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

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "correlated_likelihood_sensitivity",
  "correlated_likelihood_helpers.R"
))

source_output_id <- get_arg(
  "--source-output-id",
  "correlated_likelihood_sensitivity"
)
output_id <- get_arg(
  "--output-id",
  "correlated_likelihood_sensitivity_finalized"
)
output_parent <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal"
)
source_directory <- file.path(output_parent, source_output_id)
final_directory <- file.path(output_parent, output_id)
if (!dir.exists(source_directory) ||
    !all(file.exists(file.path(
      source_directory,
      c("configuration.rds", "analysis.rds", "fit_bundle.rds")
    )))) {
  stop("The source raw-fit cache is incomplete.")
}
if (file.exists(final_directory)) {
  stop("Refusing to overwrite the finalized output directory.")
}
staging_directory <- file.path(
  output_parent,
  paste0(".", output_id, "_staging_", Sys.getpid())
)
dir.create(staging_directory, recursive = FALSE)
summary_directory <- file.path(staging_directory, "summary")
dir.create(summary_directory, recursive = FALSE)

old_configuration <- readRDS(file.path(source_directory, "configuration.rds"))
old_analysis <- readRDS(file.path(source_directory, "analysis.rds"))
old_fit_bundle <- readRDS(file.path(source_directory, "fit_bundle.rds"))
raw_fits <- old_fit_bundle$raw_fits
pair_metadata <- old_fit_bundle$pair_metadata
method_ids <- old_configuration$method_ids
method_labels <- old_configuration$method_labels
alpha <- old_configuration$alpha
if (!identical(names(raw_fits), method_ids) ||
    !identical(names(method_labels), method_ids) ||
    !is.data.frame(pair_metadata) || nrow(pair_metadata) != 6362L ||
    anyDuplicated(pair_metadata$pair_key)) {
  stop("The source raw fits, labels, or pair metadata are invalid.")
}
for (method_id in method_ids) {
  fit <- raw_fits[[method_id]]
  fit_keys <- names(fit$lfdr)
  if (is.null(fit_keys)) {
    fit_keys <- rownames(fit$posterior_weights)
  }
  if (!identical(as.character(fit_keys), pair_metadata$pair_key) ||
      any(!is.finite(fit$lfdr)) || any(fit$lfdr < 0 | fit$lfdr > 1)) {
    stop("The raw fit is invalid for ", method_id, ".")
  }
  names(raw_fits[[method_id]]$lfdr) <- pair_metadata$pair_key
  rownames(raw_fits[[method_id]]$posterior_weights) <- pair_metadata$pair_key
  rownames(raw_fits[[method_id]]$L_matrix) <- pair_metadata$pair_key
}

bf_updates <- run_bf_updates_checked(
  raw_fits,
  method_labels = unname(method_labels),
  pair_keys = pair_metadata$pair_key
)
bf_fits <- bf_updates$successful_fits
bf_update_status <- bf_updates$status

old_bf_fits <- old_fit_bundle$bf_adjusted_fits
common_valid_ids <- intersect(names(bf_fits), names(old_bf_fits))
bf_rebuild_validation <- do.call(rbind, lapply(common_valid_ids, function(
    method_id) {
  old_prior <- extract_fit_prior_weights(
    old_bf_fits[[method_id]],
    paste0(method_id, " previous BF prior")
  )
  new_prior <- extract_fit_prior_weights(
    bf_fits[[method_id]],
    paste0(method_id, " rebuilt BF prior")
  )
  data.frame(
    method_id = method_id,
    method_label = unname(method_labels[method_id]),
    lfdr_maximum_difference = max(abs(
      old_bf_fits[[method_id]]$lfdr - bf_fits[[method_id]]$lfdr
    )),
    prior_weight_maximum_difference = max(abs(
      old_prior$prior_weight - new_prior$prior_weight
    )),
    stringsAsFactors = FALSE
  )
}))
if (any(bf_rebuild_validation$lfdr_maximum_difference > 1e-8) ||
    any(bf_rebuild_validation$prior_weight_maximum_difference > 1e-8)) {
  stop("The stable BF rebuild changed a previously valid BF result.")
}

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

pi0_table <- prior_weights[prior_weights$psd == 0, c(
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
raw_elapsed <- old_configuration$raw_elapsed_seconds
bf_elapsed <- bf_updates$elapsed_seconds
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

configuration <- old_configuration
configuration$bf_update <- bf_updates$strategy
configuration$bf_update_status <- bf_update_status
configuration$bf_warnings <- bf_updates$warnings
configuration$bf_elapsed_seconds <- bf_updates$elapsed_seconds
configuration$bf_rebuild_validation <- bf_rebuild_validation
configuration$bf_rebuild_source <- normalizePath(
  source_directory,
  mustWork = TRUE
)
configuration$generated_at <- format(
  Sys.time(),
  tz = "America/Chicago",
  usetz = TRUE
)

analysis <- old_analysis
analysis$configuration <- configuration
analysis$bf_update_status <- bf_update_status
analysis$bf_rebuild_validation <- bf_rebuild_validation
analysis$method_stage_summary <- method_stage_summary
analysis$prior_weights <- prior_weights
analysis$prior_pairwise_metrics <- prior_pairwise_metrics
analysis$lfdr_wide <- lfdr_wide
analysis$lfdr_long <- lfdr_long
analysis$lfdr_pairwise_metrics <- lfdr_pairwise_metrics
analysis$discovery_summary <- discovery_summary
analysis$top_lfdr_discrepancies <- top_lfdr_discrepancies

fit_bundle <- old_fit_bundle
fit_bundle$configuration <- configuration
fit_bundle$raw_fits <- raw_fits
fit_bundle$bf_adjusted_fits <- bf_fits
fit_bundle$bf_update_status <- bf_update_status

saveRDS(configuration, file.path(staging_directory, "configuration.rds"))
saveRDS(analysis, file.path(staging_directory, "analysis.rds"))
saveRDS(fit_bundle, file.path(staging_directory, "fit_bundle.rds"))
write_csv(pair_metadata, file.path(summary_directory, "pair_metadata.csv"))
write_csv(
  analysis$correlation_diagnostics,
  file.path(summary_directory, "correlation_diagnostics.csv")
)
write_csv(
  analysis$correlation_matrices_long,
  file.path(summary_directory, "correlation_matrices_long.csv")
)
write_csv(
  analysis$identity_validation_table,
  file.path(summary_directory, "identity_path_validation.csv")
)
write_csv(
  bf_update_status,
  file.path(summary_directory, "bf_update_status.csv")
)
write_csv(
  bf_rebuild_validation,
  file.path(summary_directory, "bf_rebuild_validation.csv")
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
  stop("Could not finalize the rebuilt sensitivity cache.")
}

cat("Rebuilt correlated-likelihood summaries from raw fits.\n")
print(bf_update_status)
print(method_stage_summary[, c(
  "fit_stage", "method_label", "pi0", "discovered_units",
  "result_status"
)])
