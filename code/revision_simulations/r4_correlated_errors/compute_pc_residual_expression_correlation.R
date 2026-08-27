#!/usr/bin/env Rscript

# Compute direct PC-residualized expression correlations for the R4 units.

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
required_roles <- c("expression_matrix", "time_specific_pc_data")
if (!all(required_roles %in% names(source_paths)) ||
    any(!file.exists(source_paths[required_roles]))) {
  stop("The z-scale source provenance does not provide readable expression and PC inputs.")
}

output_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "real_data",
  "r4_pc_residual_expression_correlation"
)
if (file.exists(output_dir)) {
  stop("Refusing to overwrite existing direct-expression R4 output: ", output_dir)
}
dir.create(output_dir, recursive = TRUE)
summary_dir <- file.path(output_dir, "summary")
dir.create(summary_dir)
completed <- FALSE
on.exit({
  if (!completed) {
    message("Direct-expression R4 computation did not complete; no complete.flag was written.")
  }
}, add = TRUE)

message("Loading original expression and time-specific PC data.")
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
donor_ids <- sort(unique(sub("_[0-9]+$", "", sample_ids)))
if (length(donor_ids) != 19L || anyDuplicated(donor_ids)) {
  stop("The expression data do not encode the expected 19 donors.")
}
time_grid <- 0:15
expected_sample_counts <- c(19L, 19L, 16L, 19L, 16L, 19L, 19L, 19L,
                            19L, 19L, 19L, 19L, 19L, 18L, 19L, 19L)

message("Residualizing each unit's expression on intercept plus PC1--PC5.")
unit_results <- lapply(seq_len(nrow(selected_units)), function(unit_index) {
  residual_expression <- matrix(
    NA_real_, nrow = length(donor_ids), ncol = length(time_grid),
    dimnames = list(donor_ids, paste0("time_", time_grid))
  )
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
    expression_vector <- as.numeric(expression_data[
      gene_rows[unit_index], match(current_sample_ids, names(expression_data))
    ])
    projection <- make_covariate_projection(covariates)
    residual_expression[donors, time_index] <- as.numeric(
      projection$residualizer %*% expression_vector
    )
  }
  pairwise_n <- crossprod(!is.na(residual_expression))
  correlation <- stats::cor(residual_expression, use = "pairwise.complete.obs")
  summary <- summarize_r4_correlation(correlation, time_grid)
  list(
    selected_unit = selected_units[unit_index, , drop = FALSE],
    pc_residual_expression = residual_expression,
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
  stop("The direct expression correlations do not have adequate pairwise donor coverage.")
}

configuration <- list(
  schema_version = "r4-pc-residual-expression-correlation-v1",
  source_z_cache = normalizePath(z_result_path, winslash = "/", mustWork = TRUE),
  n_units = nrow(selected_units),
  time_grid = time_grid,
  expression_definition = paste(
    "Quantile-normalized expression residualized separately at each time point",
    "on intercept plus PC1--PC5; genotype is deliberately not included."
  ),
  correlation_definition = "stats::cor(..., use = 'pairwise.complete.obs')",
  source_information = source_information
)
result <- list(
  configuration = configuration,
  selected_units = selected_units,
  unit_results = unit_results
)
saveRDS(result, file.path(output_dir, "pc_residual_expression_correlation.rds"))
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
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(pairwise_n_table, file.path(summary_dir, "pairwise_donor_counts.csv"), row.names = FALSE)
writeLines(c(
  "result_id=r4_pc_residual_expression_correlation",
  "n_units=6",
  "residualization=intercept_plus_PC1_to_PC5",
  "genotype_in_residualization=false",
  "correlation=pairwise_complete_donors",
  "minimum_pairwise_donors=13"
), file.path(output_dir, "complete.flag"))
completed <- TRUE
message("Completed PC-residual expression correlation output: ", output_dir)
