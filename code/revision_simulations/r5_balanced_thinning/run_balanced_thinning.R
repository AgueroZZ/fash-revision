#!/usr/bin/env Rscript

# Run the R5 BF-updated balanced-variant thinning sensitivity analysis.

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
  if (length(parsed) < 1L || anyNA(parsed) || any(parsed < 0L) ||
      anyDuplicated(parsed)) {
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
  hash <- strsplit(output, "[[:space:]]+")[[1]][1]
  if (!grepl("^[0-9a-f]{64}$", hash)) {
    stop("Unexpected SHA-256 output for ", path, ".")
  }
  hash
}

file_metadata <- function(path) {
  info <- file.info(path)
  if (nrow(info) != 1L || is.na(info$size)) {
    stop("Could not read source-file metadata: ", path)
  }
  list(
    path = normalizePath(path, mustWork = TRUE),
    size_bytes = unname(info$size),
    modification_time = format(info$mtime, tz = "UTC", usetz = TRUE),
    sha256 = sha256_file(path)
  )
}

validate_full_fit <- function(fit, name) {
  required_fields <- c(
    "prior_weights", "posterior_weights", "psd_grid", "lfdr",
    "settings", "fash_data", "L_matrix"
  )
  pair_keys <- names(fit$fash_data$data_list)
  n_units <- length(pair_keys)
  valid <-
    all(required_fields %in% names(fit)) &&
    n_units >= 2L &&
    all(nzchar(pair_keys)) &&
    anyDuplicated(pair_keys) == 0L &&
    length(fit$lfdr) == n_units &&
    all(is.finite(fit$lfdr)) &&
    all(fit$lfdr >= 0 & fit$lfdr <= 1) &&
    is.matrix(fit$L_matrix) &&
    nrow(fit$L_matrix) == n_units &&
    ncol(fit$L_matrix) == length(fit$psd_grid) &&
    length(fit$psd_grid) >= 2L &&
    all(is.finite(fit$psd_grid)) &&
    sum(fit$psd_grid == 0) == 1L
  if (!isTRUE(valid)) {
    stop(name, " does not contain a valid aligned FASH fit.")
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

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "r5_balanced_thinning",
  "balanced_thinning_helpers.R"
))

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required.")
}

seeds <- parse_integer_list(
  get_arg(
    "--seeds",
    "12345,22345,32345,42345,52345,62345,72345,82345,92345,102345"
  ),
  "seeds"
)
target_per_gene <- as.integer(get_arg("--target-per-gene", "20"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
output_id <- get_arg("--output-id", "r5_balanced_thinning")
if (is.na(target_per_gene) || target_per_gene < 1L ||
    !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
    !nzchar(output_id) || grepl("/", output_id, fixed = TRUE)) {
  stop("Invalid target_per_gene, alpha, or output_id.")
}

expected_seeds <- seq(12345L, 102345L, by = 10000L)
expected_full_units <- 1009173L
expected_full_genes <- 6362L
expected_target <- 20L
expected_eligible_genes <- 6352L
expected_excluded_genes <- 10L
expected_selected_units <- 127040L
expected_full_discovered_pairs <- 9205L
expected_full_discovered_genes <- 1177L
expected_penalty <- 10
expected_grid <- sort(c(0, exp(-0.5 * seq(0, 10, by = 0.2))))
expected_source_sizes <- c(
  raw_fit = 576110361,
  bf_adjusted_fit = 581122554
)
expected_source_hashes <- c(
  raw_fit = "9e9aaf8a405f7ca83439990656666035fc0a74bb1ae2858ace79944fac6ec929",
  bf_adjusted_fit = "3e3d6735b9da734a3ab16ed53713908c02d0af4e15a659760fc5f892ef64023b"
)
if (!identical(seeds, expected_seeds) || target_per_gene != expected_target) {
  stop("The retained R5 cache requires the ten confirmed seeds and target 20.")
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
  "revision_simulations"
)
output_directory <- file.path(output_parent, output_id)
if (!file.exists(raw_fit_path) || !file.exists(bf_fit_path) ||
    !dir.exists(output_parent)) {
  stop("The full FASH fits or revision output directory are missing.")
}
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = FALSE)
}
seed_directories <- file.path(output_directory, paste0("seed", seeds))
if (any(file.exists(seed_directories))) {
  stop(
    "Refusing to overwrite existing R5 seed directories: ",
    paste(seed_directories[file.exists(seed_directories)], collapse = ", ")
  )
}
shared_files <- file.path(
  output_directory,
  c(
    "configuration.rds", "full_variant_counts.csv", "excluded_genes.csv",
    "full_discovery_summary.csv", "full_discovered_genes.csv",
    "seed_summary.csv"
  )
)
if (any(file.exists(shared_files))) {
  stop("Refusing to overwrite existing shared R5 cache files.")
}

source_files <- list(
  raw_fit = file_metadata(raw_fit_path),
  bf_adjusted_fit = file_metadata(bf_fit_path)
)
if (!identical(
  unname(vapply(source_files, `[[`, numeric(1), "size_bytes")),
  unname(expected_source_sizes)
) || !identical(
  unname(vapply(source_files, `[[`, character(1), "sha256")),
  unname(expected_source_hashes)
)) {
  stop("The full-fit source files do not match the confirmed R5 inputs.")
}

message("Loading and validating the full BF-adjusted reference fit.")
full_bf <- load_exact_object(bf_fit_path, "fash_fit1_update")
validate_full_fit(full_bf, "fash_fit1_update")
pair_keys <- names(full_bf$fash_data$data_list)
full_bf_lfdr <- as.numeric(full_bf$lfdr)
names(full_bf_lfdr) <- pair_keys
full_bf_prior <- full_bf$prior_weights
full_psd_grid <- full_bf$psd_grid
full_call_indices <- cumulative_fdr_calls(full_bf_lfdr, alpha = alpha)
full_discovered_gene_set <- sort(unique(sub(
  "_.*$",
  "",
  pair_keys[full_call_indices]
)), method = "radix")
if (length(pair_keys) != expected_full_units ||
    length(unique(sub("_.*$", "", pair_keys))) != expected_full_genes ||
    length(full_call_indices) != expected_full_discovered_pairs ||
    length(full_discovered_gene_set) != expected_full_discovered_genes ||
    !isTRUE(all.equal(full_psd_grid, expected_grid))) {
  stop("The full BF-adjusted reference does not match confirmed dimensions.")
}
rm(full_bf)
gc(verbose = FALSE)

message("Loading and validating the full raw likelihood source.")
full_raw <- load_exact_object(raw_fit_path, "fash_fit1")
validate_full_fit(full_raw, "fash_fit1")
if (!identical(names(full_raw$fash_data$data_list), pair_keys) ||
    !isTRUE(all.equal(full_raw$psd_grid, full_psd_grid)) ||
    !identical(full_raw$settings$penalty, expected_penalty)) {
  stop("The raw and BF-adjusted full fits are not aligned.")
}
full_raw_likelihood <- full_raw$L_matrix
full_raw_settings <- full_raw$settings
rm(full_raw)
gc(verbose = FALSE)

gene_index <- make_gene_index(pair_keys)
variant_counts <- data.frame(
  gene_id = names(gene_index),
  n_tested_variants = as.integer(lengths(gene_index)),
  stringsAsFactors = FALSE
)
eligible_gene_ids <- variant_counts$gene_id[
  variant_counts$n_tested_variants >= target_per_gene
]
excluded_genes <- variant_counts[
  variant_counts$n_tested_variants < target_per_gene,
  ,
  drop = FALSE
]
excluded_genes$full_data_discovered_gene <-
  excluded_genes$gene_id %in% full_discovered_gene_set
if (nrow(variant_counts) != expected_full_genes ||
    sum(variant_counts$n_tested_variants) != expected_full_units ||
    length(eligible_gene_ids) != expected_eligible_genes ||
    nrow(excluded_genes) != expected_excluded_genes ||
    any(excluded_genes$full_data_discovered_gene) ||
    !all(full_discovered_gene_set %in% eligible_gene_ids)) {
  stop("The confirmed 20-variant eligible-gene universe failed validation.")
}

full_discovery_summary <- data.frame(
  fit = "Full-data BF-adjusted FASH(1)",
  n_tested_pairs = expected_full_units,
  n_tested_genes = expected_full_genes,
  discovered_pairs = length(full_call_indices),
  discovered_genes = length(full_discovered_gene_set),
  eligible_genes = length(eligible_gene_ids),
  eligible_full_discovered_genes = sum(
    full_discovered_gene_set %in% eligible_gene_ids
  ),
  excluded_genes = nrow(excluded_genes),
  excluded_full_discovered_genes = sum(
    excluded_genes$full_data_discovered_gene
  ),
  alpha = alpha,
  stringsAsFactors = FALSE
)

staging_directories <- character(length(seeds))
completed_summaries <- vector("list", length(seeds))
names(staging_directories) <- names(completed_summaries) <- as.character(seeds)
on.exit({
  remaining <- staging_directories[
    nzchar(staging_directories) & dir.exists(staging_directories)
  ]
  if (length(remaining) > 0L) {
    unlink(remaining, recursive = TRUE, force = TRUE)
  }
}, add = TRUE)

for (seed_index in seq_along(seeds)) {
  seed <- seeds[seed_index]
  final_directory <- seed_directories[seed_index]
  staging_directory <- file.path(
    output_directory,
    paste0(".seed", seed, "_staging_", Sys.getpid())
  )
  if (file.exists(staging_directory)) {
    stop("Unexpected R5 staging-directory collision: ", staging_directory)
  }
  dir.create(staging_directory, recursive = FALSE)
  staging_directories[seed_index] <- staging_directory

  message("Selecting 20 variants per eligible gene for seed ", seed, ".")
  selection <- select_balanced_variants_per_gene(
    pair_keys,
    seed = seed,
    target_per_gene = target_per_gene,
    gene_index = gene_index
  )
  if (nrow(selection) != expected_selected_units ||
      length(unique(selection$gene_id)) != expected_eligible_genes ||
      any(table(selection$gene_id) != target_per_gene)) {
    stop("Seed ", seed, " did not produce the confirmed balanced selection.")
  }

  message("Refitting BF-updated mixture weights for seed ", seed, ".")
  elapsed_start <- proc.time()[["elapsed"]]
  fit_capture <- capture_fit_warnings(refit_bf_from_cached_likelihood(
    full_likelihood = full_raw_likelihood,
    selected_indices = selection$fash_index,
    selected_pair_keys = selection$pair_key,
    psd_grid = full_psd_grid,
    penalty = expected_penalty
  ))
  elapsed_seconds <- proc.time()[["elapsed"]] - elapsed_start
  thinned_bf <- fit_capture$value
  if (length(thinned_bf$lfdr) != expected_selected_units ||
      any(!is.finite(thinned_bf$lfdr)) ||
      any(thinned_bf$lfdr < 0 | thinned_bf$lfdr > 1) ||
      !identical(rownames(thinned_bf$posterior_weights), selection$pair_key)) {
    stop("Seed ", seed, " produced an invalid BF-updated refit.")
  }

  prior_comparison <- compare_prior_weights(
    full_bf_prior,
    thinned_bf$prior_weights
  )
  lfdr_comparison <- compare_paired_lfdr(
    full_lfdr = full_bf_lfdr[selection$fash_index],
    thinned_lfdr = thinned_bf$lfdr,
    pair_keys = selection$pair_key,
    alpha = alpha
  )
  lfdr_table <- lfdr_comparison$table
  lfdr_table$gene_id <- selection$gene_id
  lfdr_table$variant_id <- selection$variant_id
  lfdr_table <- lfdr_table[, c(
    "pair_key", "gene_id", "variant_id", "full_lfdr", "thinned_lfdr",
    "lfdr_difference", "absolute_difference", "full_fdr_call",
    "thinned_fdr_call"
  )]
  called_gene_ids <- selection$gene_id[lfdr_comparison$thinned_call_indices]
  discovered_gene_counts <- sort(table(called_gene_ids), decreasing = TRUE)
  discovered_genes <- data.frame(
    gene_id = names(discovered_gene_counts),
    discovered_pairs = as.integer(discovered_gene_counts),
    full_data_discovered_gene =
      names(discovered_gene_counts) %in% full_discovered_gene_set,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  comparison_summary <- cbind(
    data.frame(
      seed = seed,
      fit_stage = "BF-updated",
      target_per_gene = target_per_gene,
      n_full_units = expected_full_units,
      n_full_genes = expected_full_genes,
      n_eligible_genes = expected_eligible_genes,
      n_selected_units = expected_selected_units,
      stringsAsFactors = FALSE
    ),
    prior_comparison$summary,
    lfdr_comparison$summary,
    data.frame(
      thinned_discovered_genes = nrow(discovered_genes),
      thinned_discovered_full_genes = sum(
        discovered_genes$full_data_discovered_gene
      ),
      thinned_discovered_additional_genes = sum(
        !discovered_genes$full_data_discovered_gene
      ),
      elapsed_seconds = unname(elapsed_seconds),
      warning_count = length(fit_capture$warnings),
      stringsAsFactors = FALSE
    )
  )
  configuration <- list(
    experiment = paste(
      "R5 BF-updated FASH(1) refit after uniformly sampling exactly",
      "20 tested variants per eligible gene"
    ),
    seed = seed,
    target_per_gene = target_per_gene,
    selection_rule = paste(
      "Uniform sampling without replacement within every gene with at least",
      "20 tested variants; independent of beta, SE, likelihood, lfdr, and calls."
    ),
    n_full_units = expected_full_units,
    n_full_genes = expected_full_genes,
    n_eligible_genes = expected_eligible_genes,
    n_excluded_genes = expected_excluded_genes,
    n_selected_units = expected_selected_units,
    alpha = alpha,
    fit_stage = "BF-updated",
    raw_penalty = expected_penalty,
    refit_strategy = paste(
      "Subset fixed full-data likelihood rows, re-estimate empirical-Bayes",
      "weights with penalty 10, then apply BF_update."
    ),
    source_files = source_files,
    full_raw_settings = full_raw_settings,
    full_psd_grid = full_psd_grid,
    warnings = fit_capture$warnings,
    elapsed_seconds = unname(elapsed_seconds),
    r_version = R.version.string,
    package_versions = c(
      fashr = as.character(utils::packageVersion("fashr"))
    ),
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )

  saveRDS(configuration, file.path(staging_directory, "configuration.rds"))
  saveRDS(selection, file.path(staging_directory, "selection.rds"))
  saveRDS(
    lfdr_table,
    file.path(staging_directory, "lfdr_comparison.rds"),
    compress = "gzip"
  )
  write_csv(
    prior_comparison$table,
    file.path(staging_directory, "prior_weight_comparison.csv")
  )
  write_csv(
    comparison_summary,
    file.path(staging_directory, "comparison_summary.csv")
  )
  write_csv(
    discovered_genes,
    file.path(staging_directory, "discovered_genes.csv")
  )
  if (!file.rename(staging_directory, final_directory)) {
    stop("Could not finalize R5 seed directory: ", final_directory)
  }
  staging_directories[seed_index] <- ""
  completed_summaries[[seed_index]] <- comparison_summary
  message(
    "Completed seed ", seed, ": ",
    nrow(discovered_genes), " unique discovered genes in ",
    format(unname(elapsed_seconds), digits = 4), " seconds."
  )
  rm(
    selection, fit_capture, thinned_bf, prior_comparison, lfdr_comparison,
    lfdr_table, discovered_gene_counts, discovered_genes, configuration
  )
  gc(verbose = FALSE)
}

seed_summary <- do.call(rbind, completed_summaries)
rownames(seed_summary) <- NULL
shared_configuration <- list(
  experiment = "Revision Real-data Sensitivity R5: balanced variant thinning",
  seeds = seeds,
  target_per_gene = target_per_gene,
  alpha = alpha,
  fit_stage = "BF-updated",
  n_full_units = expected_full_units,
  n_full_genes = expected_full_genes,
  n_eligible_genes = expected_eligible_genes,
  n_excluded_genes = expected_excluded_genes,
  n_selected_units_per_seed = expected_selected_units,
  source_files = source_files,
  full_raw_settings = full_raw_settings,
  full_psd_grid = full_psd_grid,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
saveRDS(shared_configuration, file.path(output_directory, "configuration.rds"))
write_csv(
  variant_counts,
  file.path(output_directory, "full_variant_counts.csv")
)
write_csv(
  excluded_genes,
  file.path(output_directory, "excluded_genes.csv")
)
write_csv(
  full_discovery_summary,
  file.path(output_directory, "full_discovery_summary.csv")
)
write_csv(
  data.frame(
    gene_id = full_discovered_gene_set,
    stringsAsFactors = FALSE
  ),
  file.path(output_directory, "full_discovered_genes.csv")
)
write_csv(
  seed_summary,
  file.path(output_directory, "seed_summary.csv")
)

cat(
  "\nR5 balanced-thinning BF-updated refits completed for seeds: ",
  paste(seeds, collapse = ", "),
  "\n",
  sep = ""
)
