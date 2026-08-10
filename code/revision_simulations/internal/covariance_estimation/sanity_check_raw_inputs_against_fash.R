#!/usr/bin/env Rscript

# Reproduce a deterministic random subset of processed FASH inputs from raw data.

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
  args <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(args, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(args[equals_hit[1L]], nchar(equals_prefix) + 1L))
  }
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1L] == length(args)) {
    return(default)
  }
  args[hit[1L] + 1L]
}

workflowr_root <- find_workflowr_root()
project_root <- normalizePath(file.path(workflowr_root, ".."), mustWork = TRUE)
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "covariance_estimation",
  "donor_null_permutation_helpers.R"
))

n_units <- as.integer(get_arg("--n-units", "50"))
random_seed <- as.integer(get_arg("--seed", "20260806"))
output_id <- get_arg("--output-id", "raw_input_fash_sanity_check")
if (is.na(n_units) || n_units < 1L || is.na(random_seed) || !nzchar(output_id)) {
  stop("Invalid subset size, seed, or output ID.")
}

fit_path <- file.path(
  workflowr_root,
  "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
)
vcf_path <- file.path(
  project_root,
  "iPSC-data", "genotype-data", "YRI_genotype.vcf.gz"
)
expression_path <- file.path(
  project_root,
  "iPSC-data", "expression-data", "quantile_normalized_no_projection.txt"
)
pc_path <- file.path(
  project_root,
  "iPSC-data", "pc-data", "principal_components_10.txt"
)
required_inputs <- c(fit_path, vcf_path, expression_path, pc_path)
if (any(!file.exists(required_inputs))) {
  stop(
    "At least one required input is missing: ",
    paste(required_inputs[!file.exists(required_inputs)], collapse = ", ")
  )
}

output_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal", output_id
)
summary_dir <- file.path(output_dir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading the fitted FASH object and selecting random units.")
fit_environment <- new.env(parent = emptyenv())
loaded_names <- load(fit_path, envir = fit_environment)
if (!identical(loaded_names, "fash_fit1_update")) {
  stop("The fitted-object file must contain only fash_fit1_update.")
}
fash_fit <- fit_environment$fash_fit1_update
fit_keys <- names(fash_fit$fash_data$data_list)
if (length(fit_keys) != 1009173L ||
    length(fash_fit$fash_data$S) != length(fit_keys) ||
    any(!nzchar(fit_keys)) || anyDuplicated(fit_keys) ||
    n_units > length(fit_keys)) {
  stop("The fitted FASH input structure is not the expected version.")
}

set.seed(random_seed)
selected_fit_indices <- sort(sample.int(length(fit_keys), n_units))
selected_keys <- fit_keys[selected_fit_indices]
selected_gene_ids <- sub("_[^_]+$", "", selected_keys)
selected_variant_ids <- sub("^.*_", "", selected_keys)
if (any(!nzchar(selected_gene_ids)) || any(!nzchar(selected_variant_ids))) {
  stop("At least one selected FASH key could not be parsed.")
}

message("Reading selected genes, copied PC data, and selected VCF dosages.")
expression_data <- utils::read.csv(
  expression_path,
  sep = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
pc_data <- utils::read.delim(
  pc_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!identical(names(expression_data)[1L], "Gene_id") ||
    anyDuplicated(expression_data$Gene_id) ||
    anyDuplicated(pc_data$Sample_id) ||
    !all(c("Sample_id", paste0("PC", 1:5)) %in% names(pc_data))) {
  stop("The expression or PC input has invalid identifiers.")
}
expression_sample_ids <- names(expression_data)[-1L]
if (length(expression_sample_ids) != 297L ||
    !setequal(expression_sample_ids, pc_data$Sample_id)) {
  stop("The copied PC file does not match the 297 expression sample IDs.")
}
gene_rows <- match(selected_gene_ids, expression_data$Gene_id)
if (anyNA(gene_rows)) {
  stop("At least one randomly selected gene is missing from expression data.")
}

unique_variant_ids <- unique(selected_variant_ids)
unique_dosage <- read_selected_vcf_dosages(
  vcf_path,
  unique_variant_ids,
  chunk_size = 100000L
)
selected_dosage <- unique_dosage[
  ,
  match(selected_variant_ids, unique_variant_ids),
  drop = FALSE
]
colnames(selected_dosage) <- selected_keys

time_grid <- 0:15
refit_beta <- matrix(
  NA_real_,
  nrow = n_units,
  ncol = length(time_grid),
  dimnames = list(selected_keys, paste0("time_", time_grid))
)
raw_se <- refit_beta
residual_df <- integer(length(time_grid))
sample_counts <- integer(length(time_grid))

message("Refitting 800 time-specific regressions with lm().")
for (time_index in seq_along(time_grid)) {
  time_value <- time_grid[time_index]
  sample_ids <- grep(
    paste0("_", time_value, "$"),
    expression_sample_ids,
    value = TRUE
  )
  donors <- sub(paste0("_", time_value, "$"), "", sample_ids)
  pc_rows <- match(sample_ids, pc_data$Sample_id)
  expression_columns <- match(sample_ids, names(expression_data))
  expression_matrix <- as.matrix(
    expression_data[gene_rows, expression_columns, drop = FALSE]
  )
  storage.mode(expression_matrix) <- "double"
  genotype_matrix <- selected_dosage[donors, , drop = FALSE]
  pc_matrix <- as.matrix(
    pc_data[pc_rows, paste0("PC", 1:5), drop = FALSE]
  )
  storage.mode(pc_matrix) <- "double"
  sample_counts[time_index] <- length(sample_ids)

  for (unit_index in seq_len(n_units)) {
    regression_data <- data.frame(
      y = as.numeric(expression_matrix[unit_index, ]),
      g = as.numeric(genotype_matrix[, unit_index]),
      PC1 = pc_matrix[, 1L],
      PC2 = pc_matrix[, 2L],
      PC3 = pc_matrix[, 3L],
      PC4 = pc_matrix[, 4L],
      PC5 = pc_matrix[, 5L]
    )
    fit <- stats::lm(y ~ g + PC1 + PC2 + PC3 + PC4 + PC5,
                     data = regression_data)
    coefficient_table <- summary(fit)$coefficients
    if (!"g" %in% rownames(coefficient_table)) {
      stop("A selected genotype was singular at time ", time_value, ".")
    }
    refit_beta[unit_index, time_index] <- coefficient_table["g", "Estimate"]
    raw_se[unit_index, time_index] <- coefficient_table["g", "Std. Error"]
    if (unit_index == 1L) {
      residual_df[time_index] <- fit$df.residual
    } else if (fit$df.residual != residual_df[time_index]) {
      stop("Residual degrees of freedom vary across units at one time point.")
    }
  }
}

adjusted_se <- convert_raw_to_t_adjusted_se(refit_beta, raw_se, residual_df)
fit_beta <- t(vapply(selected_fit_indices, function(index) {
  unit_data <- fash_fit$fash_data$data_list[[index]]
  if (!identical(as.numeric(unit_data$x), as.numeric(time_grid))) {
    stop("A selected fitted unit does not use the expected time grid.")
  }
  as.numeric(unit_data$y)
}, numeric(length(time_grid))))
fit_adjusted_se <- t(vapply(selected_fit_indices, function(index) {
  as.numeric(fash_fit$fash_data$S[[index]])
}, numeric(length(time_grid))))
dimnames(fit_beta) <- dimnames(refit_beta)
dimnames(fit_adjusted_se) <- dimnames(refit_beta)

comparison_rows <- data.frame(
  fit_index = rep(selected_fit_indices, each = length(time_grid)),
  pair_key = rep(selected_keys, each = length(time_grid)),
  gene_id = rep(selected_gene_ids, each = length(time_grid)),
  variant_id = rep(selected_variant_ids, each = length(time_grid)),
  time = rep(time_grid, times = n_units),
  n_donors = rep(sample_counts, times = n_units),
  residual_df = rep(residual_df, times = n_units),
  refit_beta = as.vector(t(refit_beta)),
  fash_beta = as.vector(t(fit_beta)),
  beta_difference = as.vector(t(refit_beta - fit_beta)),
  raw_regression_se = as.vector(t(raw_se)),
  refit_adjusted_se = as.vector(t(adjusted_se)),
  fash_adjusted_se = as.vector(t(fit_adjusted_se)),
  adjusted_se_difference = as.vector(t(adjusted_se - fit_adjusted_se)),
  stringsAsFactors = FALSE
)

absolute_beta_difference <- abs(comparison_rows$beta_difference)
absolute_se_difference <- abs(comparison_rows$adjusted_se_difference)
tolerance <- 1e-10
numerical_summary <- data.frame(
  quantity = c(
    "Number of random units",
    "Number of unit-time comparisons",
    "Maximum absolute beta difference",
    "Mean absolute beta difference",
    "Maximum absolute adjusted-SE difference",
    "Mean absolute adjusted-SE difference",
    "Comparisons exceeding tolerance"
  ),
  value = c(
    n_units,
    nrow(comparison_rows),
    max(absolute_beta_difference),
    mean(absolute_beta_difference),
    max(absolute_se_difference),
    mean(absolute_se_difference),
    sum(absolute_beta_difference > tolerance |
          absolute_se_difference > tolerance)
  ),
  stringsAsFactors = FALSE
)
selected_units <- data.frame(
  fit_index = selected_fit_indices,
  pair_key = selected_keys,
  gene_id = selected_gene_ids,
  variant_id = selected_variant_ids,
  stringsAsFactors = FALSE
)
source_information <- data.frame(
  role = c("fitted_fash", "genotype_vcf", "expression", "copied_pc_data"),
  path = required_inputs,
  size_bytes = unname(file.info(required_inputs)$size),
  mtime = format(file.info(required_inputs)$mtime, tz = "UTC", usetz = TRUE),
  md5 = unname(tools::md5sum(required_inputs)),
  stringsAsFactors = FALSE
)
configuration <- list(
  analysis_id = output_id,
  created_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  random_seed = random_seed,
  n_units = n_units,
  n_time_points = length(time_grid),
  regression = "Time-specific Y ~ 1 + G + PC1 + ... + PC5 via lm()",
  se_processing = paste(
    "Raw regression t statistics converted to two-sided normal z statistics",
    "using the time-specific residual degrees of freedom"
  ),
  tolerance = tolerance,
  source_information = source_information,
  r_version = R.version.string
)
result <- list(
  configuration = configuration,
  selected_units = selected_units,
  sample_counts = sample_counts,
  residual_df = residual_df,
  comparison_rows = comparison_rows,
  numerical_summary = numerical_summary
)

saveRDS(configuration, file.path(output_dir, "configuration.rds"))
saveRDS(result, file.path(output_dir, "sanity_check_results.rds"), compress = "xz")
utils::write.csv(
  selected_units,
  file.path(summary_dir, "selected_units.csv"),
  row.names = FALSE
)
utils::write.csv(
  comparison_rows,
  file.path(summary_dir, "comparison_rows.csv"),
  row.names = FALSE
)
utils::write.csv(
  numerical_summary,
  file.path(summary_dir, "numerical_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  source_information,
  file.path(summary_dir, "source_information.csv"),
  row.names = FALSE
)

print(numerical_summary, row.names = FALSE)
if (any(absolute_beta_difference > tolerance) ||
    any(absolute_se_difference > tolerance)) {
  stop("At least one raw-data reconstruction discrepancy exceeds tolerance.")
}
cat("Raw-input FASH sanity check passed.\n")
