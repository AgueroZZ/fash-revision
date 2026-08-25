#!/usr/bin/env Rscript

# Run the R5 BF-updated one-variant-per-gene sensitivity analysis.

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

package_provenance <- function(package) {
  description <- utils::packageDescription(package)
  remote_sha <- description[["RemoteSha"]]
  if (is.null(remote_sha) || length(remote_sha) != 1L || is.na(remote_sha)) {
    remote_sha <- ""
  }
  list(
    package = package,
    version = as.character(description[["Version"]]),
    remote_sha = as.character(remote_sha)
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
  "r5_one_variant_per_gene",
  "one_variant_per_gene_helpers.R"
))

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required.")
}
expected_fashr_provenance <- list(
  package = "fashr",
  version = "0.1.43",
  remote_sha = "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
)
fashr_provenance <- package_provenance("fashr")
if (!identical(fashr_provenance, expected_fashr_provenance)) {
  stop(
    "The retained R5 analysis requires fashr 0.1.43 at RemoteSha ",
    expected_fashr_provenance$remote_sha, "."
  )
}

expected_seeds <- seq(12345L, by = 10000L, length.out = 100L)
seeds <- parse_integer_list(
  get_arg(
    "--seeds",
    paste(expected_seeds, collapse = ",")
  ),
  "seeds"
)
alpha <- as.numeric(get_arg("--alpha", "0.05"))
output_id <- get_arg(
  "--output-id",
  "r5_one_variant_per_gene_100seed_fashr0143"
)
if (!is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
    !nzchar(output_id) || grepl("/", output_id, fixed = TRUE)) {
  stop("Invalid alpha or output_id.")
}

expected_full_units <- 1009173L
expected_full_genes <- 6362L
expected_target <- 1L
expected_selected_units <- expected_full_genes
expected_full_discovered_pairs <- 9214L
expected_full_discovered_genes <- 1176L
expected_full_pi0 <- 0.938159265061590
expected_penalty <- 10
expected_grid <- sort(c(0, exp(-0.5 * seq(0, 10, by = 0.2))))
expected_source_sizes <- c(
  raw_fit = 576110361,
  bf_adjusted_fit = 581123504
)
expected_source_hashes <- c(
  raw_fit = "9e9aaf8a405f7ca83439990656666035fc0a74bb1ae2858ace79944fac6ec929",
  bf_adjusted_fit = "7f0ca9ab0fbeab89a13c83d2a0fb7c24195f7b5a5835f209399cf0e359001f50"
)
if (!identical(seeds, expected_seeds) ||
    !identical(output_id, "r5_one_variant_per_gene_100seed_fashr0143")) {
  stop(
    paste(
      "The retained R5 cache requires the 100 prespecified seeds and",
      "output_id r5_one_variant_per_gene_100seed_fashr0143."
    )
  )
}

total_elapsed_start <- proc.time()[["elapsed"]]

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
    "configuration.rds", "full_variant_counts.csv",
    "full_discovery_summary.csv", "full_discovered_genes.csv",
    "seed_summary.csv", "gene_discovery_frequency.csv"
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
full_null_weight <- full_bf_prior$prior_weight[full_bf_prior$psd == 0]
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
    length(full_null_weight) != 1L ||
    !isTRUE(all.equal(full_null_weight, expected_full_pi0)) ||
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
all_gene_ids <- variant_counts$gene_id
if (nrow(variant_counts) != expected_full_genes ||
    sum(variant_counts$n_tested_variants) != expected_full_units ||
    length(all_gene_ids) != expected_full_genes ||
    anyDuplicated(all_gene_ids) ||
    !all(full_discovered_gene_set %in% all_gene_ids)) {
  stop("The confirmed one-variant-per-gene universe failed validation.")
}

full_discovery_summary <- data.frame(
  fit = "Full-data BF-adjusted FASH(1)",
  n_tested_pairs = expected_full_units,
  n_tested_genes = expected_full_genes,
  discovered_pairs = length(full_call_indices),
  discovered_genes = length(full_discovered_gene_set),
  representative_genes_per_seed = expected_full_genes,
  alpha = alpha,
  stringsAsFactors = FALSE
)

staging_directories <- character(length(seeds))
completed_summaries <- vector("list", length(seeds))
completed_gene_sets <- vector("list", length(seeds))
names(staging_directories) <- names(completed_summaries) <-
  names(completed_gene_sets) <- as.character(seeds)
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

  message("Selecting one variant per gene for seed ", seed, ".")
  selection <- select_one_variant_per_gene(
    pair_keys,
    seed = seed,
    gene_index = gene_index
  )
  if (nrow(selection) != expected_selected_units ||
      length(unique(selection$gene_id)) != expected_full_genes ||
      any(table(selection$gene_id) != 1L)) {
    stop("Seed ", seed, " did not produce one pair for every tested gene.")
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
      target_per_gene = expected_target,
      n_full_units = expected_full_units,
      n_full_genes = expected_full_genes,
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
    cache_id = output_id,
    experiment = paste(
      "R5 BF-updated FASH(1) refit after uniformly sampling exactly",
      "one tested variant per gene"
    ),
    seed = seed,
    target_per_gene = expected_target,
    selection_rule = paste(
      "One uniform draw from every tested gene; independent of beta, SE,",
      "likelihood, lfdr, and calls."
    ),
    n_full_units = expected_full_units,
    n_full_genes = expected_full_genes,
    n_selected_units = expected_selected_units,
    alpha = alpha,
    fit_stage = "BF-updated",
    raw_penalty = expected_penalty,
    refit_strategy = paste(
      "Subset fixed full-data likelihood rows, re-estimate empirical-Bayes",
      "weights with penalty 10, then apply BF_update."
    ),
    source_files = source_files,
    package_provenance = fashr_provenance,
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
  completed_gene_sets[[seed_index]] <- discovered_genes$gene_id
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
if (nrow(seed_summary) != length(expected_seeds) ||
    !identical(seed_summary$seed, expected_seeds) ||
    any(seed_summary$warning_count != 0L) ||
    any(seed_summary$n_selected_units != expected_selected_units)) {
  stop("The completed seed summaries failed the retained-cache contract.")
}
gene_discovery_frequency <- summarize_gene_discovery_frequency(
  discovered_gene_sets = completed_gene_sets,
  seeds = seeds,
  all_genes = all_gene_ids,
  full_discovered_genes = full_discovered_gene_set
)
total_elapsed_seconds <- proc.time()[["elapsed"]] - total_elapsed_start
shared_configuration <- list(
  cache_id = output_id,
  experiment = paste(
    "Revision Real-data Sensitivity R5: one uniformly sampled variant",
    "per gene across 100 prespecified seeds"
  ),
  seeds = seeds,
  target_per_gene = expected_target,
  alpha = alpha,
  fit_stage = "BF-updated",
  n_full_units = expected_full_units,
  n_full_genes = expected_full_genes,
  n_selected_units_per_seed = expected_selected_units,
  selection_rule = paste(
    "Uniform outcome-independent sampling of one tested pair per gene",
    "within each seed."
  ),
  monte_carlo_interpretation = paste(
    "Seed variation quantifies sensitivity to representative-variant choice;",
    "it is not a cell-line bootstrap or a formal confidence interval."
  ),
  source_files = source_files,
  package_provenance = fashr_provenance,
  full_raw_settings = full_raw_settings,
  full_psd_grid = full_psd_grid,
  total_elapsed_seconds = unname(total_elapsed_seconds),
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
saveRDS(shared_configuration, file.path(output_directory, "configuration.rds"))
write_csv(
  variant_counts,
  file.path(output_directory, "full_variant_counts.csv")
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
write_csv(
  gene_discovery_frequency,
  file.path(output_directory, "gene_discovery_frequency.csv")
)

cat(
  "\nR5 one-variant-per-gene BF-updated refits completed for ",
  length(seeds),
  " seeds in ",
  format(unname(total_elapsed_seconds), digits = 6),
  " seconds.\n",
  sep = ""
)
