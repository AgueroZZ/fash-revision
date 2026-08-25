#!/usr/bin/env Rscript

# Estimate unit-specific null correlations from signal-stripped donor residual blocks.

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
  arguments <- commandArgs(trailingOnly = TRUE)
  prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(arguments, prefix))
  if (length(equals_hit) > 0L) {
    return(substring(arguments[equals_hit[1L]], nchar(prefix) + 1L))
  }
  flag_hit <- which(arguments == name)
  if (length(flag_hit) == 0L || flag_hit[1L] == length(arguments)) {
    return(default)
  }
  arguments[flag_hit[1L] + 1L]
}

file_metadata <- function(paths, roles) {
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  information <- file.info(paths)
  data.frame(
    role = roles,
    path = paths,
    size_bytes = unname(information$size),
    modification_time = format(information$mtime, tz = "UTC", usetz = TRUE),
    md5 = unname(tools::md5sum(paths)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

workflowr_root <- find_workflowr_root()
project_root <- normalizePath(file.path(workflowr_root, ".."), mustWork = TRUE)
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "covariance_estimation", "donor_null_permutation_helpers.R"
))
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation",
  "selected_signal_genotype_permutation_helpers.R"
))
source(file.path(
  workflowr_root, "code", "revision_simulations", "r4_correlated_errors",
  "unit_specific_residual_permutation_helpers.R"
))

n_units <- as.integer(get_arg("--n-units", "6"))
n_permutations <- as.integer(get_arg("--n-permutations", "400"))
selection_seed <- as.integer(get_arg("--selection-seed", "20260824"))
permutation_seed <- as.integer(get_arg("--permutation-seed", "20260825"))
output_id <- get_arg("--output-id", "r4_unit_specific_residual_permutation_z")
if (length(n_units) != 1L || is.na(n_units) || n_units < 2L ||
    length(n_permutations) != 1L || is.na(n_permutations) ||
    n_permutations < 20L || length(selection_seed) != 1L ||
    is.na(selection_seed) || length(permutation_seed) != 1L ||
    is.na(permutation_seed) || !nzchar(output_id) ||
    grepl("/", output_id, fixed = TRUE)) {
  stop("Invalid R4 unit count, permutation count, seeds, or output ID.")
}

bf_fit_path <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
)
vcf_path <- file.path(project_root, "iPSC-data", "genotype-data", "YRI_genotype.vcf.gz")
expression_path <- file.path(
  project_root, "iPSC-data", "expression-data", "quantile_normalized_no_projection.txt"
)
pc_path <- file.path(
  workflowr_root, "data", "dynamic_eQTL_real", "principal_components_10.txt"
)
source_paths <- c(bf_fit_path, vcf_path, expression_path, pc_path)
source_roles <- c("bf_adjusted_fash", "genotype_vcf", "expression_matrix", "time_specific_pc_data")
if (any(!file.exists(source_paths))) {
  stop("A required input is missing: ", paste(source_paths[!file.exists(source_paths)], collapse = ", "))
}

output_parent <- file.path(workflowr_root, "output", "revision_simulations", "real_data")
output_dir <- file.path(output_parent, output_id)
if (file.exists(output_dir)) {
  stop("Refusing to overwrite existing R4 unit-specific output: ", output_dir)
}
dir.create(output_dir, recursive = TRUE)
summary_dir <- file.path(output_dir, "summary")
dir.create(summary_dir)
completed <- FALSE
on.exit({
  if (!completed) {
    message("R4 unit-specific run did not complete; no complete.flag was written.")
  }
}, add = TRUE)

message("Loading the BF-adjusted FASH fit and randomly selecting units.")
fit_environment <- new.env(parent = emptyenv())
loaded_names <- load(bf_fit_path, envir = fit_environment)
if (!identical(loaded_names, "fash_fit1_update")) {
  stop("The BF-adjusted fit must contain only fash_fit1_update.")
}
fash_fit <- fit_environment$fash_fit1_update
pair_keys <- names(fash_fit$fash_data$data_list)
if (length(pair_keys) < n_units || length(fash_fit$lfdr) != length(pair_keys) ||
    !identical(names(fash_fit$lfdr), pair_keys) || anyDuplicated(pair_keys)) {
  stop("The BF-adjusted FASH fit is not aligned by unique pair key.")
}
selected_units <- select_random_r4_units(pair_keys, n_units, selection_seed)
selected_units$source_bf_lfdr <- unname(fash_fit$lfdr[selected_units$source_fash_index])
selected_data <- fash_fit$fash_data$data_list[selected_units$source_fash_index]
saved_beta <- t(vapply(selected_data, function(unit_data) {
  if (!identical(as.numeric(unit_data$x), as.numeric(0:15))) {
    stop("A selected FASH unit does not use the expected time grid.")
  }
  as.numeric(unit_data$y)
}, numeric(16L)))
rownames(saved_beta) <- selected_units$pair_key
rm(fit_environment, fash_fit)
gc(verbose = FALSE)

message("Reading expression, covariate, and selected dosage inputs.")
expression_data <- utils::read.csv(
  expression_path, sep = "", check.names = FALSE, stringsAsFactors = FALSE
)
pc_data <- utils::read.delim(pc_path, check.names = FALSE, stringsAsFactors = FALSE)
if (!identical(names(expression_data)[1L], "Gene_id") ||
    !all(c("Sample_id", paste0("PC", 1:5)) %in% names(pc_data)) ||
    anyDuplicated(expression_data$Gene_id) || anyDuplicated(pc_data$Sample_id)) {
  stop("The expression or PC input has invalid identifiers.")
}
sample_ids <- names(expression_data)[-1L]
if (!setequal(sample_ids, pc_data$Sample_id)) {
  stop("Expression and PC sample IDs do not match exactly.")
}
gene_rows <- match(selected_units$gene_id, expression_data$Gene_id)
if (anyNA(gene_rows)) {
  stop("At least one randomly selected gene is missing from the expression matrix.")
}
unique_variant_ids <- unique(selected_units$variant_id)
unique_dosage <- read_selected_vcf_dosages(vcf_path, unique_variant_ids)
unit_dosage <- unique_dosage[, match(selected_units$variant_id, unique_variant_ids), drop = FALSE]
colnames(unit_dosage) <- selected_units$pair_key
donor_ids <- rownames(unit_dosage)
if (nrow(unit_dosage) != 19L || anyDuplicated(donor_ids)) {
  stop("Selected VCF dosage must contain 19 unique donors.")
}

time_grid <- 0:15
expected_sample_counts <- c(19L, 19L, 16L, 19L, 16L, 19L, 19L, 19L,
                            19L, 19L, 19L, 19L, 19L, 18L, 19L, 19L)
donor_observation_matrix <- vapply(time_grid, function(time_value) {
  paste0(donor_ids, "_", time_value) %in% sample_ids
}, logical(length(donor_ids)))
rownames(donor_observation_matrix) <- donor_ids
observation_patterns <- apply(donor_observation_matrix, 1L, paste0, collapse = "")

make_unit_input <- function(unit_index) {
  time_inputs <- lapply(seq_along(time_grid), function(time_index) {
    time_value <- time_grid[time_index]
    current_sample_ids <- grep(paste0("_", time_value, "$"), sample_ids, value = TRUE)
    donors <- sub(paste0("_", time_value, "$"), "", current_sample_ids)
    if (length(current_sample_ids) != expected_sample_counts[time_index] ||
        anyDuplicated(donors) || any(!donors %in% donor_ids)) {
      stop("Unexpected donor coverage at time ", time_value, ".")
    }
    pc_rows <- match(current_sample_ids, pc_data$Sample_id)
    list(
      donors = donors,
      expression = matrix(
        as.numeric(expression_data[gene_rows[unit_index], match(current_sample_ids, names(expression_data))]),
        ncol = 1L,
        dimnames = list(donors, selected_units$pair_key[unit_index])
      ),
      genotype = unit_dosage[donors, unit_index, drop = FALSE],
      covariates = as.matrix(pc_data[pc_rows, paste0("PC", 1:5), drop = FALSE])
    )
  })
  lapply(time_inputs, function(input) {
    storage.mode(input$genotype) <- "double"
    storage.mode(input$covariates) <- "double"
    rownames(input$covariates) <- input$donors
    input
  })
}
unit_inputs <- lapply(seq_len(n_units), make_unit_input)

message("Checking that reconstructed observed regressions reproduce FASH betas.")
observed_beta <- matrix(NA_real_, nrow = n_units, ncol = length(time_grid))
for (unit_index in seq_len(n_units)) {
  for (time_index in seq_along(time_grid)) {
    input <- unit_inputs[[unit_index]][[time_index]]
    observed_beta[unit_index, time_index] <- fit_many_genotype_regressions(
      input$expression, input$genotype, input$covariates
    )$beta
  }
}
rownames(observed_beta) <- selected_units$pair_key
maximum_observed_beta_difference <- apply(abs(observed_beta - saved_beta), 1L, max)
if (any(maximum_observed_beta_difference > 1e-10)) {
  stop("A reconstructed observed beta does not reproduce the saved FASH input.")
}
selected_units$maximum_observed_beta_difference <- maximum_observed_beta_difference

message("Generating ", n_permutations, " full donor-trajectory residual permutations per unit.")
unit_results <- vector("list", n_units)
all_map_rows <- vector("list", n_units * n_permutations)
map_row_index <- 0L
for (unit_index in seq_len(n_units)) {
  beta_draws <- matrix(NA_real_, nrow = n_permutations, ncol = length(time_grid))
  raw_se_draws <- beta_draws
  residual_correlation <- matrix(NA_real_, nrow = n_permutations, ncol = length(time_grid))
  residual_df <- integer(length(time_grid))
  for (draw_index in seq_len(n_permutations)) {
    draw_seed <- permutation_seed + 100000L * unit_index + draw_index
    donor_map <- make_r4_donor_map(donor_ids, observation_patterns, draw_seed)
    map_row_index <- map_row_index + 1L
    donor_map$selected_order <- unit_index
    donor_map$pair_key <- selected_units$pair_key[unit_index]
    donor_map$draw_index <- draw_index
    donor_map$draw_seed <- draw_seed
    all_map_rows[[map_row_index]] <- donor_map
    for (time_index in seq_along(time_grid)) {
      input <- unit_inputs[[unit_index]][[time_index]]
      source_donors <- donor_map$source_donor[match(input$donors, donor_map$target_donor)]
      source_rows <- match(source_donors, input$donors)
      if (anyNA(source_rows)) {
        stop("A donor map is incompatible with observed donors at time ", time_grid[time_index], ".")
      }
      null_fit <- make_signal_stripped_residual_block_null(
        expression = input$expression,
        genotype = input$genotype,
        covariates = input$covariates,
        source_rows = source_rows
      )
      beta_draws[draw_index, time_index] <- null_fit$null_fit$beta
      raw_se_draws[draw_index, time_index] <- null_fit$null_fit$standard_error
      residual_correlation[draw_index, time_index] <-
        null_fit$maximum_residual_genotype_correlation
      if (draw_index == 1L) {
        residual_df[time_index] <- null_fit$null_fit$residual_df
      } else if (residual_df[time_index] != null_fit$null_fit$residual_df) {
        stop("Residual degrees of freedom changed across permutation draws.")
      }
    }
  }
  colnames(beta_draws) <- paste0("time_", time_grid)
  colnames(raw_se_draws) <- colnames(beta_draws)
  adjusted_se_draws <- convert_raw_to_original_t_adjusted_se(
    beta_draws, raw_se_draws, residual_df
  )
  z_draws <- beta_draws / adjusted_se_draws
  if (any(!is.finite(z_draws))) {
    stop("The t-adjusted null z-score draws are invalid.")
  }
  summary <- summarize_r4_null_draws(z_draws, time_grid)
  unit_results[[unit_index]] <- list(
    selected_unit = selected_units[unit_index, , drop = FALSE],
    observed_beta = observed_beta[unit_index, ],
    saved_beta = saved_beta[unit_index, ],
    null_beta_draws = beta_draws,
    null_raw_se_draws = raw_se_draws,
    null_t_adjusted_se_draws = adjusted_se_draws,
    null_z_draws = z_draws,
    residual_df = residual_df,
    maximum_residual_genotype_correlation = max(residual_correlation),
    correlation = summary$correlation,
    variogram = summary$variogram
  )
}
names(unit_results) <- selected_units$pair_key
donor_maps <- do.call(rbind, all_map_rows)
rownames(donor_maps) <- NULL
if (any(!is.finite(vapply(unit_results, `[[`, numeric(1), "maximum_residual_genotype_correlation"))) ||
    any(vapply(unit_results, `[[`, numeric(1), "maximum_residual_genotype_correlation") > 1e-10)) {
  stop("A source full-model residual was not orthogonal to the genotype.")
}

configuration <- list(
  schema_version = "r4-unit-specific-residual-permutation-z-v1",
  output_id = output_id,
  n_units = n_units,
  n_permutations = n_permutations,
  selection_seed = selection_seed,
  permutation_seed = permutation_seed,
  time_grid = time_grid,
  estimand = paste(
    "Per-unit correlation of t-adjusted signal-stripped residual-permutation z scores",
    "conditional on observed genotype, covariates, donor coverage, and expression residuals"
  ),
  t_adjustment = paste(
    "Per-draw raw regression SEs are converted using",
    "convert_raw_to_original_t_adjusted_se() before z-score construction."
  ),
  residual_randomization = paste(
    "Full-model residual blocks are permuted by one donor map per draw and unit,",
    "shared across all time points and restricted to donor missingness-pattern strata."
  ),
  source_information = file_metadata(source_paths, source_roles)
)
result <- list(
  configuration = configuration,
  selected_units = selected_units,
  donor_maps = donor_maps,
  unit_results = unit_results
)
saveRDS(result, file.path(output_dir, "unit_specific_null_correlation.rds"))
utils::write.csv(selected_units, file.path(summary_dir, "selected_units.csv"), row.names = FALSE)
utils::write.csv(donor_maps, file.path(summary_dir, "donor_maps.csv"), row.names = FALSE)
variogram_table <- do.call(rbind, lapply(unit_results, function(unit_result) {
  transform(unit_result$variogram, pair_key = unit_result$selected_unit$pair_key)
}))
utils::write.csv(variogram_table, file.path(summary_dir, "variograms.csv"), row.names = FALSE)
writeLines(c(
  paste0("output_id=", output_id),
  paste0("n_units=", n_units),
  paste0("n_permutations=", n_permutations),
  "residual_method=signal_stripped_full_donor_trajectory_permutation"
), file.path(output_dir, "complete.flag"))
completed <- TRUE
message("Completed R4 unit-specific residual-permutation output: ", output_dir)
