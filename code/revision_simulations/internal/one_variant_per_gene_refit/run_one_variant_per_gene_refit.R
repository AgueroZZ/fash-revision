#!/usr/bin/env Rscript

# Refit FASH(1) after uniformly sampling one tested variant per gene.

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
  if (!all(required_fields %in% names(fit)) ||
      is.null(fit$fash_data$data_list) || is.null(fit$fash_data$S)) {
    stop(name, " does not contain the required FASH fields.")
  }
  n_units <- length(fit$fash_data$data_list)
  pair_keys <- names(fit$fash_data$data_list)
  if (n_units < 2L || length(fit$fash_data$S) != n_units ||
      length(fit$lfdr) != n_units || length(pair_keys) != n_units ||
      any(!nzchar(pair_keys)) || anyDuplicated(pair_keys) ||
      any(!is.finite(fit$lfdr)) || any(fit$lfdr < 0 | fit$lfdr > 1) ||
      length(fit$psd_grid) < 2L || any(!is.finite(fit$psd_grid)) ||
      !any(fit$psd_grid == 0)) {
    stop(name, " has invalid dimensions, pair keys, lfdr values, or PSD grid.")
  }
  invisible(TRUE)
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

fit_direct_fash <- function(datasets, full_fit, num_cores) {
  settings <- full_fit$settings
  required_settings <- c(
    "num_basis", "order", "betaprec", "pred_step", "likelihood", "penalty"
  )
  if (!all(required_settings %in% names(settings))) {
    stop("The full FASH settings are incomplete.")
  }
  capture_fit_warnings(fashr::fash(
    Y = "beta",
    smooth_var = "time",
    S = "SE",
    data_list = datasets,
    num_basis = settings$num_basis,
    order = settings$order,
    betaprec = settings$betaprec,
    pred_step = settings$pred_step,
    likelihood = settings$likelihood,
    penalty = settings$penalty,
    grid = full_fit$psd_grid,
    num_cores = num_cores,
    verbose = TRUE
  ))
}

validate_thinned_fit <- function(fit, pair_keys, name) {
  if (length(fit$lfdr) != length(pair_keys) || any(!is.finite(fit$lfdr)) ||
      any(fit$lfdr < 0 | fit$lfdr > 1) ||
      nrow(fit$posterior_weights) != length(pair_keys) ||
      nrow(fit$L_matrix) != length(pair_keys)) {
    stop(name, " has invalid posterior dimensions or lfdr values.")
  }
  names(fit$lfdr) <- pair_keys
  fit
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

seeds <- parse_integer_list(get_arg("--seeds", "12345"), "seeds")
num_cores <- as.integer(get_arg("--num-cores", "8"))
validation_units <- as.integer(get_arg("--validation-units", "24"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
output_prefix <- get_arg(
  "--output-prefix",
  "one_variant_per_gene_refit_seed"
)
full_summary_file <- get_arg(
  "--full-summary-file",
  "one_variant_per_gene_full_bf_discovery_summary.csv"
)
variant_count_file <- get_arg(
  "--variant-count-file",
  "one_variant_per_gene_variant_counts.csv"
)
if (is.na(num_cores) || num_cores < 1L || is.na(validation_units) ||
    validation_units < 2L || !is.finite(alpha) ||
    alpha <= 0 || alpha >= 1 || !nzchar(output_prefix) ||
    grepl("/", output_prefix, fixed = TRUE) || !nzchar(full_summary_file) ||
    grepl("/", full_summary_file, fixed = TRUE) ||
    !nzchar(variant_count_file) ||
    grepl("/", variant_count_file, fixed = TRUE)) {
  stop(
    paste(
      "Invalid num_cores, alpha, output_prefix, full_summary_file,",
      "or variant_count_file."
    )
  )
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
if (!file.exists(raw_fit_path) || !file.exists(bf_fit_path) ||
    !dir.exists(output_parent)) {
  stop("The full fits or internal output directory are missing.")
}
final_directories <- file.path(output_parent, paste0(output_prefix, seeds))
if (any(file.exists(final_directories))) {
  stop(
    "Refusing to overwrite existing output directories: ",
    paste(final_directories[file.exists(final_directories)], collapse = ", ")
  )
}

source_files <- list(
  raw_fit = file_metadata(raw_fit_path),
  bf_adjusted_fit = file_metadata(bf_fit_path)
)
message("Loading and validating the full raw FASH fit.")
full_raw <- load_exact_object(raw_fit_path, "fash_fit1")
validate_full_fit(full_raw, "fash_fit1")
pair_keys <- names(full_raw$fash_data$data_list)
gene_index <- make_gene_index(pair_keys)
n_genes <- length(gene_index)
n_full_units <- length(pair_keys)
variant_count_by_gene <- data.frame(
  gene_id = names(gene_index),
  n_tested_variants = as.integer(lengths(gene_index)),
  stringsAsFactors = FALSE
)
if (nrow(variant_count_by_gene) != n_genes ||
    any(!nzchar(variant_count_by_gene$gene_id)) ||
    anyDuplicated(variant_count_by_gene$gene_id) ||
    any(variant_count_by_gene$n_tested_variants < 1L) ||
    sum(variant_count_by_gene$n_tested_variants) != n_full_units) {
  stop("The full-data per-gene variant counts failed validation.")
}
message(
  "Validated ", format(n_full_units, big.mark = ","),
  " full units across ", format(n_genes, big.mark = ","), " genes."
)

staging_directories <- character(length(seeds))
names(staging_directories) <- as.character(seeds)
raw_stage_results <- vector("list", length(seeds))
names(raw_stage_results) <- as.character(seeds)
likelihood_validation <- NULL

for (seed_index in seq_along(seeds)) {
  seed <- seeds[seed_index]
  output_id <- paste0(output_prefix, seed)
  staging_directory <- file.path(
    output_parent,
    paste0(".", output_id, "_staging_", Sys.getpid())
  )
  if (file.exists(staging_directory)) {
    stop("Unexpected staging-directory collision: ", staging_directory)
  }
  dir.create(staging_directory, recursive = FALSE)
  staging_directories[seed_index] <- staging_directory

  message("Selecting one random variant per gene for seed ", seed, ".")
  selection <- select_random_variant_per_gene(
    pair_keys,
    seed = seed,
    gene_index = gene_index
  )
  if (nrow(selection) != n_genes || anyDuplicated(selection$gene_id) ||
      anyDuplicated(selection$pair_key)) {
    stop("The seed ", seed, " selection is not exactly one pair per gene.")
  }
  if (seed_index == 1L) {
    validation_count <- min(validation_units, nrow(selection))
    validation_indices <- selection$fash_index[seq_len(validation_count)]
    validation_pair_keys <- selection$pair_key[seq_len(validation_count)]
    validation_datasets <- make_fash_refit_datasets(
      full_raw,
      validation_indices
    )
    message(
      "Validating ", validation_count,
      " cached likelihood rows against a direct beta/SE recomputation."
    )
    validation_start <- proc.time()[["elapsed"]]
    direct_capture <- fit_direct_fash(
      validation_datasets,
      full_raw,
      num_cores
    )
    direct_validation_fit <- validate_thinned_fit(
      direct_capture$value,
      validation_pair_keys,
      "direct_validation_fit"
    )
    cached_validation_fit <- validate_thinned_fit(
      refit_fash_from_cached_likelihood(full_raw, validation_indices),
      validation_pair_keys,
      "cached_validation_fit"
    )
    likelihood_max_abs_difference <- max(abs(
      direct_validation_fit$L_matrix - cached_validation_fit$L_matrix
    ))
    lfdr_max_abs_difference <- max(abs(
      direct_validation_fit$lfdr - cached_validation_fit$lfdr
    ))
    prior_validation <- compare_prior_weights(
      direct_validation_fit$prior_weights,
      cached_validation_fit$prior_weights,
      fit_stage = "Cached-likelihood validation"
    )
    validation_elapsed <- proc.time()[["elapsed"]] - validation_start
    likelihood_validation <- list(
      selection_seed = seed,
      n_units = validation_count,
      pair_keys = validation_pair_keys,
      likelihood_max_absolute_difference = likelihood_max_abs_difference,
      raw_lfdr_max_absolute_difference = lfdr_max_abs_difference,
      prior_total_variation =
        prior_validation$summary$prior_total_variation,
      elapsed_seconds = unname(validation_elapsed),
      warnings = direct_capture$warnings,
      tolerances = c(likelihood = 1e-8, lfdr = 1e-6, prior_tv = 1e-6)
    )
    if (likelihood_max_abs_difference > 1e-8 ||
        lfdr_max_abs_difference > 1e-6 ||
        prior_validation$summary$prior_total_variation > 1e-6) {
      stop("Cached likelihood validation failed numerical tolerances.")
    }
    rm(
      validation_datasets,
      direct_capture,
      direct_validation_fit,
      cached_validation_fit,
      prior_validation
    )
    gc(verbose = FALSE)
  }

  message("Refitting the thinned empirical-Bayes mixture for seed ", seed, ".")
  fit_start <- proc.time()[["elapsed"]]
  raw_capture <- capture_fit_warnings(
    refit_fash_from_cached_likelihood(full_raw, selection$fash_index)
  )
  raw_elapsed <- proc.time()[["elapsed"]] - fit_start
  thinned_raw <- validate_thinned_fit(
    raw_capture$value,
    selection$pair_key,
    "thinned_raw"
  )

  message("Applying the BF prior update for seed ", seed, ".")
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
  raw_stage_results[[seed_index]] <- list(
    seed = seed,
    selection = selection,
    raw_prior = raw_prior,
    raw_lfdr = raw_lfdr,
    thinned_raw = thinned_raw,
    thinned_bf = thinned_bf,
    raw_fit_warnings = raw_capture$warnings,
    bf_update_warnings = bf_capture$warnings,
    raw_eb_refit_elapsed_seconds = unname(raw_elapsed),
    bf_update_elapsed_seconds = unname(bf_elapsed)
  )
  rm(raw_capture, bf_capture)
  gc(verbose = FALSE)
}

full_raw_settings <- full_raw$settings
full_psd_grid <- full_raw$psd_grid
full_raw_prior_weights <- full_raw$prior_weights
rm(full_raw, gene_index)
gc(verbose = FALSE)

message("Loading and validating the full BF-adjusted FASH fit.")
full_bf <- load_exact_object(bf_fit_path, "fash_fit1_update")
validate_full_fit(full_bf, "fash_fit1_update")
if (!identical(names(full_bf$fash_data$data_list), pair_keys)) {
  stop("The raw and BF-adjusted full fits do not use identical pair keys.")
}
full_bf_calls <- cumulative_fdr_calls(full_bf$lfdr, alpha = alpha)
full_bf_discovery_summary <- data.frame(
  fit = "Full-data BF-adjusted FASH(1)",
  n_tested_pairs = length(full_bf$lfdr),
  n_tested_genes = length(unique(sub("_.*$", "", pair_keys))),
  discovered_pairs = length(full_bf_calls),
  discovered_genes = length(unique(sub(
    "_.*$",
    "",
    pair_keys[full_bf_calls]
  ))),
  alpha = alpha,
  stringsAsFactors = FALSE
)

completed_summaries <- vector("list", length(seeds))
for (seed_index in seq_along(seeds)) {
  stage_result <- raw_stage_results[[seed_index]]
  seed <- stage_result$seed
  selection <- stage_result$selection
  bf_prior <- compare_prior_weights(
    full_bf$prior_weights,
    stage_result$thinned_bf$prior_weights,
    fit_stage = "BF-adjusted"
  )
  bf_lfdr <- compare_paired_lfdr(
    full_bf$lfdr[selection$fash_index],
    stage_result$thinned_bf$lfdr,
    selection$pair_key,
    fit_stage = "BF-adjusted",
    alpha = alpha
  )
  prior_table <- rbind(stage_result$raw_prior$table, bf_prior$table)
  lfdr_table <- rbind(stage_result$raw_lfdr$table, bf_lfdr$table)
  prior_summary <- rbind(stage_result$raw_prior$summary, bf_prior$summary)
  lfdr_summary <- rbind(stage_result$raw_lfdr$summary, bf_lfdr$summary)
  comparison_summary <- merge(
    prior_summary,
    lfdr_summary,
    by = "fit_stage",
    all = TRUE,
    sort = FALSE
  )
  comparison_summary$seed <- seed
  comparison_summary$n_full_units <- n_full_units
  comparison_summary$n_genes <- n_genes
  comparison_summary$raw_eb_refit_elapsed_seconds <-
    stage_result$raw_eb_refit_elapsed_seconds
  comparison_summary$bf_update_elapsed_seconds <-
    stage_result$bf_update_elapsed_seconds
  comparison_summary <- comparison_summary[, c(
    "seed", "fit_stage", "n_full_units", "n_genes", "n_units",
    "full_pi0", "thinned_pi0", "pi0_difference",
    "prior_total_variation", "pearson_lfdr", "spearman_lfdr",
    "mean_absolute_lfdr_difference", "median_absolute_lfdr_difference",
    "rmse_lfdr", "full_mean_lfdr", "thinned_mean_lfdr",
    "full_fdr_calls", "thinned_fdr_calls", "fdr_call_intersection",
    "fdr_call_union", "fdr_call_jaccard", "alpha",
    "raw_eb_refit_elapsed_seconds", "bf_update_elapsed_seconds"
  )]

  configuration <- list(
    experiment = "One uniformly random tested variant per gene FASH(1) refit",
    seed = seed,
    n_full_units = n_full_units,
    n_genes = n_genes,
    n_selected_units = nrow(selection),
    alpha = alpha,
    num_cores = num_cores,
    refit_strategy = paste(
      "Subset fixed per-unit likelihood rows from the full fit and rerun",
      "only fash_eb_est; validated against a direct beta/SE recomputation."
    ),
    likelihood_validation = likelihood_validation,
    full_raw_settings = full_raw_settings,
    full_psd_grid = full_psd_grid,
    source_files = source_files,
    raw_fit_warnings = stage_result$raw_fit_warnings,
    bf_update_warnings = stage_result$bf_update_warnings,
    raw_eb_refit_elapsed_seconds = stage_result$raw_eb_refit_elapsed_seconds,
    bf_update_elapsed_seconds = stage_result$bf_update_elapsed_seconds,
    r_version = R.version.string,
    package_versions = c(fashr = as.character(utils::packageVersion("fashr"))),
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  fit_bundle <- list(
    configuration = configuration,
    selection = selection,
    raw_fit = stage_result$thinned_raw,
    bf_adjusted_fit = stage_result$thinned_bf,
    full_reference = list(
      raw_prior_weights = full_raw_prior_weights,
      bf_adjusted_prior_weights = full_bf$prior_weights,
      raw_selected_lfdr = stage_result$raw_lfdr$table[, c(
        "pair_key", "full_lfdr"
      )],
      bf_adjusted_selected_lfdr = bf_lfdr$table[, c(
        "pair_key", "full_lfdr"
      )]
    )
  )

  staging_directory <- staging_directories[seed_index]
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

  final_directory <- final_directories[seed_index]
  if (!file.rename(staging_directory, final_directory)) {
    stop("Could not finalize output directory: ", final_directory)
  }
  completed_summaries[[seed_index]] <- comparison_summary
  message("Completed seed ", seed, ": ", final_directory)
}

write_csv(
  full_bf_discovery_summary,
  file.path(
    output_parent,
    full_summary_file
  )
)
write_csv(
  variant_count_by_gene,
  file.path(
    output_parent,
    variant_count_file
  )
)

rm(full_bf, raw_stage_results)
gc(verbose = FALSE)

candidate_directories <- list.dirs(
  output_parent,
  recursive = FALSE,
  full.names = TRUE
)
candidate_directories <- candidate_directories[grepl(
  paste0("^", output_prefix, "[0-9]+$"),
  basename(candidate_directories)
)]
summary_paths <- file.path(candidate_directories, "comparison_summary.csv")
summary_paths <- summary_paths[file.exists(summary_paths)]
if (length(summary_paths) > 0L) {
  aggregate_summary <- do.call(rbind, lapply(summary_paths, utils::read.csv))
  aggregate_summary <- aggregate_summary[order(
    aggregate_summary$seed,
    match(aggregate_summary$fit_stage, c("Raw", "BF-adjusted"))
  ), , drop = FALSE]
  rownames(aggregate_summary) <- NULL
  write_csv(
    aggregate_summary,
    file.path(output_parent, paste0(output_prefix, "_summary.csv"))
  )
}

cat(
  "\nOne-variant-per-gene FASH refit completed for seeds: ",
  paste(seeds, collapse = ", "),
  "\n",
  sep = ""
)
