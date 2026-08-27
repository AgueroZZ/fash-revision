#!/usr/bin/env Rscript

# Compute full FASH-model residual correlations for the R4 selected units.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "covariance_estimation", "donor_null_permutation_helpers.R"
))
source(file.path(
  workflowr_root, "code", "revision_simulations", "r4_correlated_errors",
  "unit_specific_residual_permutation_helpers.R"
))

z_output_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "real_data",
  "r4_unit_specific_residual_permutation_z"
)
z_result_path <- file.path(z_output_dir, "unit_specific_null_correlation.rds")
z_complete_flag_path <- file.path(z_output_dir, "complete.flag")
if (!file.exists(z_result_path) || !file.exists(z_complete_flag_path)) {
  stop("The completed z-scale R4 cache is missing.")
}
z_analysis <- readRDS(z_result_path)
if (!identical(
  z_analysis$configuration$schema_version,
  "r4-unit-specific-residual-permutation-z-v1"
) || nrow(z_analysis$selected_units) != 6L ||
    anyDuplicated(z_analysis$selected_units$pair_key)) {
  stop("The z-scale R4 cache does not contain the expected six selected units.")
}

source_information <- z_analysis$configuration$source_information
source_paths <- setNames(source_information$path, source_information$role)
required_roles <- c("expression_matrix", "time_specific_pc_data", "genotype_vcf")
if (!all(required_roles %in% names(source_paths)) ||
    any(!file.exists(source_paths[required_roles]))) {
  stop("The z-scale source provenance does not provide readable full-model inputs.")
}

output_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "real_data",
  "r4_full_model_residual_expression_correlation"
)
if (file.exists(output_dir)) {
  stop("Refusing to overwrite existing full-model-residual R4 output: ", output_dir)
}
dir.create(output_dir, recursive = TRUE)
summary_dir <- file.path(output_dir, "summary")
dir.create(summary_dir)
completed <- FALSE
on.exit({
  if (!completed) {
    message("Full-model-residual R4 computation did not complete; no complete.flag was written.")
  }
}, add = TRUE)

message("Loading original expression, genotype, and time-specific PC data.")
expression_data <- utils::read.csv(
  source_paths[["expression_matrix"]],
  sep = "", check.names = FALSE, stringsAsFactors = FALSE
)
pc_data <- utils::read.delim(
  source_paths[["time_specific_pc_data"]],
  check.names = FALSE, stringsAsFactors = FALSE
)
if (!identical(names(expression_data)[1L], "Gene_id") ||
    !all(c("Sample_id", paste0("PC", 1:5)) %in% names(pc_data)) ||
    anyDuplicated(expression_data$Gene_id) || anyDuplicated(pc_data$Sample_id)) {
  stop("The expression or PC input has invalid identifiers.")
}
sample_ids <- names(expression_data)[-1L]
if (!setequal(sample_ids, pc_data$Sample_id)) {
  stop("Expression and PC sample IDs do not match exactly.")
}

selected_units <- z_analysis$selected_units
gene_rows <- match(selected_units$gene_id, expression_data$Gene_id)
if (anyNA(gene_rows)) {
  stop("At least one selected gene is missing from the expression matrix.")
}
unique_variant_ids <- unique(selected_units$variant_id)
unique_dosage <- read_selected_vcf_dosages(
  source_paths[["genotype_vcf"]], unique_variant_ids
)
unit_dosage <- unique_dosage[, match(selected_units$variant_id, unique_variant_ids), drop = FALSE]
colnames(unit_dosage) <- selected_units$pair_key
donor_ids <- rownames(unit_dosage)
if (nrow(unit_dosage) != 19L || anyDuplicated(donor_ids)) {
  stop("Selected VCF dosage must contain 19 unique donors.")
}
time_grid <- 0:15
expected_sample_counts <- c(19L, 19L, 16L, 19L, 16L, 19L, 19L, 19L,
                            19L, 19L, 19L, 19L, 19L, 18L, 19L, 19L)
saved_beta <- do.call(rbind, lapply(z_analysis$unit_results, `[[`, "saved_beta"))
rownames(saved_beta) <- selected_units$pair_key

message("Refitting Y ~ 1 + PC1 + ... + PC5 + G and retaining full-model residuals.")
unit_results <- lapply(seq_len(nrow(selected_units)), function(unit_index) {
  full_model_residual <- matrix(
    NA_real_, nrow = length(donor_ids), ncol = length(time_grid),
    dimnames = list(donor_ids, paste0("time_", time_grid))
  )
  observed_beta <- numeric(length(time_grid))
  residual_df <- integer(length(time_grid))
  for (time_index in seq_along(time_grid)) {
    time_value <- time_grid[time_index]
    current_sample_ids <- grep(paste0("_", time_value, "$"), sample_ids, value = TRUE)
    donors <- sub(paste0("_", time_value, "$"), "", current_sample_ids)
    if (length(current_sample_ids) != expected_sample_counts[time_index] ||
        anyDuplicated(donors) || any(!donors %in% donor_ids)) {
      stop("Unexpected donor coverage at time ", time_value, ".")
    }
    pc_rows <- match(current_sample_ids, pc_data$Sample_id)
    covariates <- as.matrix(pc_data[pc_rows, paste0("PC", 1:5), drop = FALSE])
    storage.mode(covariates) <- "double"
    expression_vector <- matrix(
      as.numeric(expression_data[
        gene_rows[unit_index], match(current_sample_ids, names(expression_data))
      ]),
      ncol = 1L,
      dimnames = list(donors, selected_units$pair_key[unit_index])
    )
    fit <- fit_many_genotype_regressions(
      expression = expression_vector,
      genotype = unit_dosage[donors, unit_index, drop = FALSE],
      covariates = covariates
    )
    full_model_residual[donors, time_index] <- fit$residual[, 1L]
    observed_beta[time_index] <- fit$beta
    residual_df[time_index] <- fit$residual_df
  }
  maximum_observed_beta_difference <- max(abs(observed_beta - saved_beta[unit_index, ]))
  if (maximum_observed_beta_difference > 1e-10) {
    stop("The full-model refit does not reproduce a saved FASH beta estimate.")
  }
  pairwise_n <- crossprod(!is.na(full_model_residual))
  correlation <- stats::cor(full_model_residual, use = "pairwise.complete.obs")
  summary <- summarize_r4_correlation(correlation, time_grid)
  list(
    selected_unit = selected_units[unit_index, , drop = FALSE],
    full_model_residual = full_model_residual,
    observed_beta = observed_beta,
    residual_df = residual_df,
    maximum_observed_beta_difference = maximum_observed_beta_difference,
    pairwise_donor_count = pairwise_n,
    correlation = summary$correlation,
    variogram = summary$variogram
  )
})
names(unit_results) <- selected_units$pair_key
if (any(vapply(unit_results, function(unit_result) {
  any(unit_result$pairwise_donor_count < 13L) ||
    any(!is.finite(unit_result$correlation))
}, logical(1)))) {
  stop("The full-model residual correlations do not have adequate pairwise donor coverage.")
}

configuration <- list(
  schema_version = "r4-full-model-residual-correlation-v1",
  source_z_cache = normalizePath(z_result_path, winslash = "/", mustWork = TRUE),
  n_units = nrow(selected_units),
  time_grid = time_grid,
  model = "Y ~ 1 + PC1 + ... + PC5 + G",
  genotype_in_residualization = TRUE,
  residual_definition = "Full-model cell-line residual after fitted PC and genotype components are removed.",
  correlation_definition = "stats::cor(..., use = 'pairwise.complete.obs')",
  source_information = source_information
)
result <- list(
  configuration = configuration,
  selected_units = selected_units,
  unit_results = unit_results
)
saveRDS(result, file.path(output_dir, "full_model_residual_correlation.rds"))
utils::write.csv(selected_units, file.path(summary_dir, "selected_units.csv"), row.names = FALSE)
variogram_table <- do.call(rbind, lapply(unit_results, function(unit_result) {
  transform(unit_result$variogram, pair_key = unit_result$selected_unit$pair_key)
}))
utils::write.csv(variogram_table, file.path(summary_dir, "variograms.csv"), row.names = FALSE)
pairwise_n_table <- do.call(rbind, lapply(unit_results, function(unit_result) {
  data.frame(
    pair_key = unit_result$selected_unit$pair_key,
    minimum_pairwise_donors = min(unit_result$pairwise_donor_count),
    maximum_pairwise_donors = max(unit_result$pairwise_donor_count),
    maximum_observed_beta_difference = unit_result$maximum_observed_beta_difference,
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(pairwise_n_table, file.path(summary_dir, "pairwise_donor_counts.csv"), row.names = FALSE)
writeLines(c(
  "result_id=r4_full_model_residual_expression_correlation",
  "n_units=6",
  "model=Y~1+PC1+...+PC5+G",
  "genotype_in_residualization=true",
  "correlation=pairwise_complete_donors",
  "minimum_pairwise_donors=13"
), file.path(output_dir, "complete.flag"))
completed <- TRUE
message("Completed full-model residual correlation output: ", output_dir)
