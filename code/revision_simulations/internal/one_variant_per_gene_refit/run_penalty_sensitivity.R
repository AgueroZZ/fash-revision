#!/usr/bin/env Rscript

# Compare penalized and unpenalized empirical-Bayes refits after thinning.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath("."))
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath("coderepo-local"))
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

parse_integer_list <- function(value, name) {
  parsed <- suppressWarnings(as.integer(strsplit(value, ",", fixed = TRUE)[[1]]))
  if (length(parsed) < 1L || anyNA(parsed) || anyDuplicated(parsed)) {
    stop("Invalid ", name, ".")
  }
  parsed
}

load_exact_object <- function(path, expected_name) {
  object_environment <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = object_environment)
  if (!identical(loaded_names, expected_name)) {
    stop(path, " must contain only ", expected_name, ".")
  }
  object_environment[[expected_name]]
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

validate_fit <- function(fit, expected_penalty = NULL, name = "fit") {
  required_fields <- c(
    "prior_weights", "posterior_weights", "psd_grid", "lfdr",
    "settings", "fash_data", "L_matrix"
  )
  if (!all(required_fields %in% names(fit)) ||
      is.null(fit$fash_data$data_list) ||
      length(fit$lfdr) != length(fit$fash_data$data_list) ||
      nrow(fit$L_matrix) != length(fit$lfdr) ||
      ncol(fit$L_matrix) != length(fit$psd_grid) ||
      any(!is.finite(fit$lfdr)) || any(fit$lfdr < 0 | fit$lfdr > 1)) {
    stop(name, " does not contain a valid FASH fit.")
  }
  if (!is.null(expected_penalty) &&
      !isTRUE(all.equal(fit$settings$penalty, expected_penalty))) {
    stop(name, " has an unexpected penalty setting.")
  }
  invisible(TRUE)
}

summarize_discoveries <- function(lfdr, pair_keys, alpha) {
  calls <- cumulative_fdr_calls(lfdr, alpha = alpha)
  genes <- sub("_.*$", "", pair_keys[calls])
  data.frame(
    pair_discoveries = length(calls),
    gene_discoveries = length(unique(genes)),
    stringsAsFactors = FALSE
  )
}

summarize_paired_lfdr_compact <- function(reference_lfdr,
                                           comparison_lfdr,
                                           pair_keys,
                                           alpha) {
  reference_lfdr <- validate_probability_vector(
    reference_lfdr,
    "reference_lfdr"
  )
  comparison_lfdr <- validate_probability_vector(
    comparison_lfdr,
    "comparison_lfdr"
  )
  if (length(reference_lfdr) != length(comparison_lfdr) ||
      length(pair_keys) != length(reference_lfdr)) {
    stop("Paired lfdr inputs must have identical lengths.")
  }
  reference_calls <- cumulative_fdr_calls(reference_lfdr, alpha = alpha)
  comparison_calls <- cumulative_fdr_calls(comparison_lfdr, alpha = alpha)
  call_union <- union(reference_calls, comparison_calls)
  call_intersection <- intersect(reference_calls, comparison_calls)
  difference <- comparison_lfdr - reference_lfdr
  data.frame(
    n_units = length(reference_lfdr),
    pearson_lfdr = safe_correlation(
      reference_lfdr,
      comparison_lfdr,
      "pearson"
    ),
    spearman_lfdr = safe_correlation(
      reference_lfdr,
      comparison_lfdr,
      "spearman"
    ),
    mean_absolute_lfdr_difference = mean(abs(difference)),
    maximum_absolute_lfdr_difference = max(abs(difference)),
    rmse_lfdr = sqrt(mean(difference^2)),
    reference_fdr_calls = length(reference_calls),
    comparison_fdr_calls = length(comparison_calls),
    fdr_call_intersection = length(call_intersection),
    fdr_call_union = length(call_union),
    fdr_call_jaccard = if (length(call_union) == 0L) {
      1
    } else {
      length(call_intersection) / length(call_union)
    },
    reference_gene_calls = length(unique(sub(
      "_.*$",
      "",
      pair_keys[reference_calls]
    ))),
    comparison_gene_calls = length(unique(sub(
      "_.*$",
      "",
      pair_keys[comparison_calls]
    ))),
    alpha = alpha,
    stringsAsFactors = FALSE
  )
}

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
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

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required.")
}

seeds <- parse_integer_list(
  get_arg("--seeds", "12345,22345,32345,42345,52345"),
  "seeds"
)
alpha <- as.numeric(get_arg("--alpha", "0.05"))
output_id <- get_arg(
  "--output-id",
  "one_variant_per_gene_penalty_sensitivity"
)
if (!is.finite(alpha) || alpha <= 0 || alpha >= 1 || !nzchar(output_id) ||
    grepl("/", output_id, fixed = TRUE)) {
  stop("Invalid alpha or output_id.")
}

output_parent <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal"
)
output_directory <- file.path(output_parent, output_id)
if (file.exists(output_directory)) {
  stop("Refusing to overwrite existing output directory: ", output_directory)
}
staging_directory <- file.path(
  output_parent,
  paste0(".", output_id, "_staging_", Sys.getpid())
)
if (file.exists(staging_directory)) {
  stop("Unexpected staging-directory collision: ", staging_directory)
}
dir.create(staging_directory, recursive = FALSE)

penalized_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "fash_fit1_all.RData"
)
unpenalized_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "fash_fit1_unpenalized_all.RData"
)
bf_adjusted_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "fash_fit1_update.RData"
)
source_paths <- c(
  penalized_raw = penalized_path,
  unpenalized_raw = unpenalized_path,
  bf_adjusted = bf_adjusted_path
)
if (any(!file.exists(source_paths))) {
  stop("One or more full-data source fits are missing.")
}
source_files <- lapply(source_paths, file_metadata)

seed_cache_directories <- file.path(
  output_parent,
  paste0("one_variant_per_gene_refit_seed", seeds)
)
seed_bundle_paths <- file.path(seed_cache_directories, "thinned_fash_fit.rds")
seed_selection_paths <- file.path(seed_cache_directories, "selection.csv")
if (any(!file.exists(seed_bundle_paths)) || any(!file.exists(seed_selection_paths))) {
  stop("One or more retained seed caches are missing.")
}
selections <- lapply(seed_selection_paths, function(path) {
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
})
audit_indices <- sort(unique(unlist(lapply(selections, `[[`, "fash_index"))))

message("Loading the full penalty-10 fit.")
full_penalized <- load_exact_object(penalized_path, "fash_fit1")
validate_fit(full_penalized, expected_penalty = 10, name = "fash_fit1")
full_pair_keys <- names(full_penalized$fash_data$data_list)
if (is.null(full_pair_keys) || any(!nzchar(full_pair_keys)) ||
    anyDuplicated(full_pair_keys)) {
  stop("The full penalty-10 fit has invalid pair keys.")
}
penalized_settings <- full_penalized$settings
penalized_grid <- full_penalized$psd_grid
penalized_prior <- full_penalized$prior_weights
penalized_full_lfdr <- as.numeric(full_penalized$lfdr)
penalized_audit_likelihood <- full_penalized$L_matrix[
  audit_indices,
  ,
  drop = FALSE
]
penalized_discoveries <- summarize_discoveries(
  penalized_full_lfdr,
  full_pair_keys,
  alpha
)
rm(full_penalized)
gc(verbose = FALSE)

message("Loading the full penalty-1 fit.")
full_unpenalized <- load_exact_object(
  unpenalized_path,
  "fash_fit1_unpenalized"
)
validate_fit(
  full_unpenalized,
  expected_penalty = 1,
  name = "fash_fit1_unpenalized"
)
unpenalized_settings <- full_unpenalized$settings
settings_without_penalty <- setdiff(names(penalized_settings), "penalty")
if (!identical(names(full_unpenalized$fash_data$data_list), full_pair_keys) ||
    !isTRUE(all.equal(full_unpenalized$psd_grid, penalized_grid)) ||
    !isTRUE(all.equal(
      full_unpenalized$settings[settings_without_penalty],
      penalized_settings[settings_without_penalty]
    ))) {
  stop("The full penalty-1 and penalty-10 fits are not otherwise aligned.")
}
likelihood_audit_max_difference <- max(abs(
  full_unpenalized$L_matrix[audit_indices, , drop = FALSE] -
    penalized_audit_likelihood
))
if (!is.finite(likelihood_audit_max_difference) ||
    likelihood_audit_max_difference > 1e-12) {
  stop("The audited full likelihood rows differ between penalty settings.")
}
unpenalized_prior <- full_unpenalized$prior_weights
unpenalized_full_lfdr <- as.numeric(full_unpenalized$lfdr)
unpenalized_discoveries <- summarize_discoveries(
  unpenalized_full_lfdr,
  full_pair_keys,
  alpha
)
full_penalty_lfdr_comparison <- summarize_paired_lfdr_compact(
  reference_lfdr = penalized_full_lfdr,
  comparison_lfdr = unpenalized_full_lfdr,
  pair_keys = full_pair_keys,
  alpha = alpha
)
full_penalty_prior_comparison <- compare_prior_weights(
  penalized_prior,
  unpenalized_prior,
  fit_stage = "Full raw penalty 10 versus penalty 1"
)$summary
rm(full_unpenalized, penalized_audit_likelihood)
gc(verbose = FALSE)

message("Loading the full BF-adjusted fit and counting discoveries.")
full_bf <- load_exact_object(bf_adjusted_path, "fash_fit1_update")
validate_fit(full_bf, name = "fash_fit1_update")
if (!identical(names(full_bf$fash_data$data_list), full_pair_keys) ||
    !isTRUE(all.equal(full_bf$psd_grid, penalized_grid))) {
  stop("The BF-adjusted full fit is not aligned with the raw fits.")
}
bf_discoveries <- summarize_discoveries(full_bf$lfdr, full_pair_keys, alpha)
bf_pi0 <- full_bf$prior_weights$prior_weight[
  full_bf$prior_weights$psd == 0
]
rm(full_bf)
gc(verbose = FALSE)

full_fit_summary <- rbind(
  data.frame(
    Fit = "Raw, penalty = 10",
    Penalty = 10,
    pi0 = penalized_prior$prior_weight[penalized_prior$psd == 0],
    penalized_discoveries,
    check.names = FALSE
  ),
  data.frame(
    Fit = "Raw, penalty = 1",
    Penalty = 1,
    pi0 = unpenalized_prior$prior_weight[unpenalized_prior$psd == 0],
    unpenalized_discoveries,
    check.names = FALSE
  ),
  data.frame(
    Fit = "BF-adjusted",
    Penalty = NA_real_,
    pi0 = bf_pi0,
    bf_discoveries,
    check.names = FALSE
  )
)

seed_summary_rows <- vector("list", length(seeds))
seed_specific_unpenalized_lfdr <- NULL
for (seed_index in seq_along(seeds)) {
  seed <- seeds[seed_index]
  message("Running the penalty sensitivity for seed ", seed, ".")
  selection <- selections[[seed_index]]
  bundle <- readRDS(seed_bundle_paths[seed_index])
  required_bundle_fields <- c("selection", "raw_fit", "bf_adjusted_fit")
  if (!all(required_bundle_fields %in% names(bundle)) ||
      !identical(selection$pair_key, bundle$selection$pair_key) ||
      !identical(selection$fash_index, bundle$selection$fash_index)) {
    stop("The seed ", seed, " fit bundle and selection do not agree.")
  }
  validate_fit(bundle$raw_fit, expected_penalty = 10, name = "thinned raw fit")
  validate_fit(bundle$bf_adjusted_fit, name = "thinned BF-adjusted fit")
  if (!identical(names(bundle$raw_fit$fash_data$data_list), selection$pair_key)) {
    stop("The seed ", seed, " thinned fit has unexpected pair keys.")
  }

  thinned_unpenalized <- refit_fash_from_cached_likelihood(
    bundle$raw_fit,
    seq_len(nrow(selection)),
    penalty = 1
  )
  validate_fit(
    thinned_unpenalized,
    expected_penalty = 1,
    name = "thinned unpenalized fit"
  )
  unpenalized_thinning_prior <- compare_prior_weights(
    unpenalized_prior,
    thinned_unpenalized$prior_weights,
    fit_stage = "Raw, penalty = 1"
  )$summary
  unpenalized_thinning_lfdr <- compare_paired_lfdr(
    unpenalized_full_lfdr[selection$fash_index],
    thinned_unpenalized$lfdr,
    selection$pair_key,
    fit_stage = "Raw, penalty = 1",
    alpha = alpha
  )
  penalty_effect_prior <- compare_prior_weights(
    bundle$raw_fit$prior_weights,
    thinned_unpenalized$prior_weights,
    fit_stage = "Thinned penalty 10 versus penalty 1"
  )$summary
  penalty_effect_lfdr <- compare_paired_lfdr(
    bundle$raw_fit$lfdr,
    thinned_unpenalized$lfdr,
    selection$pair_key,
    fit_stage = "Thinned penalty 10 versus penalty 1",
    alpha = alpha
  )$summary

  bf_from_unpenalized <- fashr::BF_update(thinned_unpenalized, plot = FALSE)
  bf_prior_invariance <- compare_prior_weights(
    bundle$bf_adjusted_fit$prior_weights,
    bf_from_unpenalized$prior_weights,
    fit_stage = "BF penalty invariance"
  )$summary
  bf_lfdr_max_difference <- max(abs(
    bundle$bf_adjusted_fit$lfdr - bf_from_unpenalized$lfdr
  ))

  seed_summary_rows[[seed_index]] <- data.frame(
    seed = seed,
    n_units = nrow(selection),
    full_penalty1_pi0 = unpenalized_thinning_prior$full_pi0,
    thinned_penalty1_pi0 = unpenalized_thinning_prior$thinned_pi0,
    thinning_penalty1_pi0_difference =
      unpenalized_thinning_prior$pi0_difference,
    thinning_penalty1_prior_total_variation =
      unpenalized_thinning_prior$prior_total_variation,
    thinning_penalty1_spearman_lfdr =
      unpenalized_thinning_lfdr$summary$spearman_lfdr,
    thinning_penalty1_lfdr_mae =
      unpenalized_thinning_lfdr$summary$mean_absolute_lfdr_difference,
    thinning_penalty1_full_calls =
      unpenalized_thinning_lfdr$summary$full_fdr_calls,
    thinning_penalty1_thinned_calls =
      unpenalized_thinning_lfdr$summary$thinned_fdr_calls,
    thinning_penalty1_call_jaccard =
      unpenalized_thinning_lfdr$summary$fdr_call_jaccard,
    thinned_penalty10_pi0 = penalty_effect_prior$full_pi0,
    thinned_penalty1_pi0_for_penalty_comparison =
      penalty_effect_prior$thinned_pi0,
    thinned_penalty1_minus_penalty10_pi0 =
      penalty_effect_prior$pi0_difference,
    thinned_penalty_prior_total_variation =
      penalty_effect_prior$prior_total_variation,
    thinned_penalty_lfdr_spearman =
      penalty_effect_lfdr$spearman_lfdr,
    thinned_penalty_lfdr_mae =
      penalty_effect_lfdr$mean_absolute_lfdr_difference,
    thinned_penalty10_calls = penalty_effect_lfdr$full_fdr_calls,
    thinned_penalty1_calls = penalty_effect_lfdr$thinned_fdr_calls,
    bf_penalty_prior_total_variation =
      bf_prior_invariance$prior_total_variation,
    bf_penalty_lfdr_max_absolute_difference = bf_lfdr_max_difference,
    alpha = alpha,
    stringsAsFactors = FALSE
  )
  if (seed_index == 1L) {
    seed_specific_unpenalized_lfdr <- unpenalized_thinning_lfdr$table
    seed_specific_unpenalized_lfdr$seed <- seed
  }
  rm(
    bundle,
    thinned_unpenalized,
    unpenalized_thinning_lfdr,
    bf_from_unpenalized
  )
  gc(verbose = FALSE)
}
seed_summary <- do.call(rbind, seed_summary_rows)
rownames(seed_summary) <- NULL

configuration <- list(
  experiment = paste(
    "Penalty sensitivity for one uniformly selected tested variant per gene"
  ),
  seeds = seeds,
  alpha = alpha,
  n_full_units = length(full_pair_keys),
  n_genes = length(unique(sub("_.*$", "", full_pair_keys))),
  n_audited_likelihood_rows = length(audit_indices),
  likelihood_audit_max_absolute_difference = likelihood_audit_max_difference,
  penalty_10_null_pseudo_units = 9L,
  full_penalized_settings = penalized_settings,
  full_unpenalized_settings = unpenalized_settings,
  full_psd_grid = penalized_grid,
  source_files = source_files,
  bf_update_dependency = paste(
    "BF_update uses only L_matrix and psd_grid from the input FASH fit;",
    "the raw penalty and raw prior weights are not used."
  ),
  r_version = R.version.string,
  package_versions = c(
    fashr = as.character(utils::packageVersion("fashr")),
    mixsqp = as.character(utils::packageVersion("mixsqp"))
  ),
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

write_csv(full_fit_summary, file.path(staging_directory, "full_fit_summary.csv"))
write_csv(
  cbind(full_penalty_prior_comparison, full_penalty_lfdr_comparison),
  file.path(staging_directory, "full_penalty_comparison.csv")
)
write_csv(seed_summary, file.path(staging_directory, "seed_summary.csv"))
write_csv(
  seed_specific_unpenalized_lfdr,
  file.path(staging_directory, "seed12345_unpenalized_lfdr_comparison.csv")
)
saveRDS(configuration, file.path(staging_directory, "configuration.rds"))

if (!file.rename(staging_directory, output_directory)) {
  stop("Could not finalize output directory: ", output_directory)
}

cat("\nPenalty sensitivity completed: ", output_directory, "\n", sep = "")
