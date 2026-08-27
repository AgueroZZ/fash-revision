#!/usr/bin/env Rscript

# Per-unit cross-time correlation of a gene's own expression, by discovery
# status, compared with the single common C used in the full-data correlated
# FASH refit.
#
# Three estimators per unit, all 16x16:
#   expr_pc_only      intercept + PC1-5 removed, genotype retained
#   expr_pc_genotype  intercept + PC1-5 + genotype removed
#   beta_scale        expr_pc_genotype propagated through the residualized
#                     genotype OLS weights; comparable to the common C

suppressWarnings(suppressPackageStartupMessages(library(stats)))

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
project_root <- normalizePath(file.path(workflowr_root, ".."), mustWork = TRUE)
source(file.path(workflowr_root, "code", "revision_simulations", "internal",
                 "per_unit_expression_correlation",
                 "per_unit_expression_correlation_helpers.R"))
source(file.path(workflowr_root, "code", "revision_simulations", "internal",
                 "covariance_estimation", "donor_null_permutation_helpers.R"))

internal_output <- file.path(workflowr_root, "output", "revision_simulations",
                             "internal")
output_dir <- file.path(internal_output, "per_unit_expression_correlation")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

panel_dir <- file.path(internal_output, "correlated_permutation_discoveries")
common_c_path <- file.path(internal_output, "permutation_correlation_export",
                           "permutation_null_correlation_mg090.rds")
expression_path <- file.path(project_root, "iPSC-data", "expression-data",
                             "quantile_normalized_no_projection.txt")
pc_path <- file.path(workflowr_root, "data", "dynamic_eQTL_real",
                     "principal_components_10.txt")
vcf_path <- file.path(project_root, "iPSC-data", "genotype-data",
                      "YRI_genotype.vcf.gz")
required <- c(common_c_path, expression_path, pc_path, vcf_path,
              file.path(panel_dir, "panel_topgene_metadata.csv"),
              file.path(panel_dir, "panel_original874_keys.csv"),
              file.path(panel_dir, "panel_discovery43_keys.csv"))
if (any(!file.exists(required))) {
  stop("Missing input: ", paste(required[!file.exists(required)], collapse = ", "))
}

# ------------------------------------------------------------------ panels ---

read_keys <- function(path) {
  as.character(utils::read.csv(path, stringsAsFactors = FALSE)$pair_key)
}
topgene <- utils::read.csv(file.path(panel_dir, "panel_topgene_metadata.csv"),
                          stringsAsFactors = FALSE)
panels <- list(
  discovery1169 = topgene$pair_key,
  null874 = read_keys(file.path(panel_dir, "panel_original874_keys.csv")),
  survivor43 = read_keys(file.path(panel_dir, "panel_discovery43_keys.csv"))
)
unit_table <- data.frame(
  pair_key = unique(unlist(panels, use.names = FALSE)),
  stringsAsFactors = FALSE
)
unit_table$gene_id <- sub("_.*$", "", unit_table$pair_key)
unit_table$variant_id <- sub("^[^_]*_", "", unit_table$pair_key)
if (anyDuplicated(unit_table$pair_key) || any(!nzchar(unit_table$variant_id))) {
  stop("The pooled unit table is invalid.")
}
message(nrow(unit_table), " unique units across ", length(panels), " panels.")

# -------------------------------------------------------------------- data ---

message("Reading expression, PC, and genotype inputs.")
expression_data <- utils::read.csv(expression_path, sep = "",
                                   check.names = FALSE,
                                   stringsAsFactors = FALSE)
pc_data <- utils::read.delim(pc_path, check.names = FALSE,
                             stringsAsFactors = FALSE)
if (!identical(names(expression_data)[1L], "Gene_id") ||
    !all(c("Sample_id", paste0("PC", 1:5)) %in% names(pc_data))) {
  stop("The expression or PC input has invalid identifiers.")
}
sample_ids <- names(expression_data)[-1L]
if (!setequal(sample_ids, pc_data$Sample_id)) {
  stop("Expression and PC sample IDs do not match exactly.")
}
gene_rows <- match(unit_table$gene_id, expression_data$Gene_id)
if (anyNA(gene_rows)) stop("A selected gene is missing from the expression matrix.")

unique_variants <- unique(unit_table$variant_id)
message("Extracting ", length(unique_variants), " variant dosages from the VCF.")
dosage <- read_selected_vcf_dosages(vcf_path, unique_variants)
donor_ids <- rownames(dosage)
if (nrow(dosage) != 19L) stop("Expected 19 donors in the VCF dosage matrix.")

time_grid <- 0:15
n_time <- length(time_grid)
time_donors <- lapply(time_grid, function(time_value) {
  current <- grep(paste0("_", time_value, "$"), sample_ids, value = TRUE)
  donors <- sub(paste0("_", time_value, "$"), "", current)
  if (anyDuplicated(donors) || any(!donors %in% donor_ids)) {
    stop("Unexpected donor coverage at time ", time_value, ".")
  }
  list(donors = donors, sample_ids = current)
})
message("Donor counts per time: ",
        paste(vapply(time_donors, function(x) length(x$donors), integer(1L)),
              collapse = " "))

# Per-time expression block (genes x donors) and covariate residualizers.
time_blocks <- lapply(seq_len(n_time), function(time_index) {
  info <- time_donors[[time_index]]
  columns <- match(info$sample_ids, names(expression_data))
  block <- as.matrix(expression_data[gene_rows, columns, drop = FALSE])
  storage.mode(block) <- "double"
  colnames(block) <- info$donors
  pc_rows <- match(info$sample_ids, pc_data$Sample_id)
  covariates <- cbind(1, as.matrix(pc_data[pc_rows, paste0("PC", 1:5),
                                           drop = FALSE]))
  storage.mode(covariates) <- "double"
  rownames(covariates) <- info$donors
  list(donors = info$donors, expression = block,
       covariates = covariates,
       residualizer = make_residualizer(covariates))
})
rm(expression_data)
invisible(gc(verbose = FALSE))

# ------------------------------------------------------------- estimation ---

n_unit <- nrow(unit_table)
estimators <- c("expr_pc_only", "expr_pc_genotype", "beta_scale")
matrices <- lapply(estimators, function(nm) {
  array(NA_real_, dim = c(n_unit, n_time, n_time),
        dimnames = list(unit_table$pair_key, paste0("day_", time_grid),
                        paste0("day_", time_grid)))
})
names(matrices) <- estimators
design_lag1 <- numeric(n_unit)

message("Estimating per-unit correlations.")
start_time <- Sys.time()
for (unit_index in seq_len(n_unit)) {
  variant_column <- match(unit_table$variant_id[unit_index], unique_variants)
  residual_pc <- vector("list", n_time)
  residual_pcg <- vector("list", n_time)
  weights <- vector("list", n_time)
  for (time_index in seq_len(n_time)) {
    block <- time_blocks[[time_index]]
    donors <- block$donors
    outcome <- block$expression[unit_index, ]
    genotype <- dosage[donors, variant_column]
    residual_pc_values <- as.numeric(block$residualizer %*% outcome)
    weight <- make_genotype_weights(block$residualizer, genotype)
    genotype_residual <- weight * sum((block$residualizer %*% genotype)^2)
    effect <- sum(weight * outcome)
    residual_pcg_values <- residual_pc_values - effect * genotype_residual
    residual_pc[[time_index]] <- setNames(residual_pc_values, donors)
    residual_pcg[[time_index]] <- setNames(residual_pcg_values, donors)
    weights[[time_index]] <- setNames(weight, donors)
  }
  correlation_pc <- matched_donor_correlation(residual_pc)
  correlation_pcg <- matched_donor_correlation(residual_pcg)
  design <- weight_design_factor(weights)
  beta_correlation <- correlation_pcg * design
  diag(beta_correlation) <- 1
  matrices$expr_pc_only[unit_index, , ] <- correlation_pc
  matrices$expr_pc_genotype[unit_index, , ] <- correlation_pcg
  matrices$beta_scale[unit_index, , ] <- beta_correlation
  design_lag1[unit_index] <- mean(design[cbind(1:(n_time - 1L), 2:n_time)])
  if (unit_index %% 250L == 0L) {
    message("  ", unit_index, " / ", n_unit, " units")
  }
}
elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
message("Estimation took ", round(elapsed, 1), " seconds.")

# ---------------------------------------------------------------- summary ---

common_c <- readRDS(common_c_path)
if (!is.matrix(common_c) || any(dim(common_c) != c(n_time, n_time))) {
  stop("The common correlation matrix has unexpected dimensions.")
}

per_unit <- do.call(rbind, lapply(estimators, function(nm) {
  data.frame(
    estimator = nm,
    pair_key = unit_table$pair_key,
    lag1 = apply(matrices[[nm]], 1L, function(m) lag_profile(m)[1L]),
    mean_off_diagonal = apply(matrices[[nm]], 1L, mean_off_diagonal),
    stringsAsFactors = FALSE
  )
}))
per_unit$gene_id <- sub("_.*$", "", per_unit$pair_key)

panel_membership <- do.call(rbind, lapply(names(panels), function(nm) {
  data.frame(panel = nm, pair_key = panels[[nm]], stringsAsFactors = FALSE)
}))

panel_lag_rows <- list()
panel_matrices <- list()
for (panel_name in names(panels)) {
  index <- match(panels[[panel_name]], unit_table$pair_key)
  for (nm in estimators) {
    mean_matrix <- apply(matrices[[nm]][index, , , drop = FALSE], c(2L, 3L), mean)
    panel_matrices[[paste(panel_name, nm, sep = "__")]] <- mean_matrix
    panel_lag_rows[[paste(panel_name, nm, sep = "__")]] <- data.frame(
      panel = panel_name, estimator = nm, n_units = length(index),
      lag = seq_len(n_time - 1L), correlation = lag_profile(mean_matrix),
      stringsAsFactors = FALSE
    )
  }
}
panel_lag_rows[["common_c"]] <- data.frame(
  panel = "common_C", estimator = "permutation_common_C", n_units = 874L,
  lag = seq_len(n_time - 1L), correlation = lag_profile(common_c),
  stringsAsFactors = FALSE
)
lag_table <- do.call(rbind, panel_lag_rows)
rownames(lag_table) <- NULL

result <- list(
  unit_table = unit_table,
  panels = panels,
  panel_membership = panel_membership,
  matrices = matrices,
  panel_matrices = panel_matrices,
  common_c = common_c,
  per_unit = per_unit,
  design_lag1 = setNames(design_lag1, unit_table$pair_key),
  lag_table = lag_table,
  topgene = topgene,
  elapsed_seconds = elapsed
)
saveRDS(result, file.path(output_dir, "per_unit_expression_correlation.rds"))
utils::write.csv(lag_table, file.path(output_dir, "panel_lag_profiles.csv"),
                 row.names = FALSE)
utils::write.csv(merge(panel_membership, per_unit, by = "pair_key"),
                 file.path(output_dir, "per_unit_summary.csv"), row.names = FALSE)
writeLines(format(Sys.time()), file.path(output_dir, "complete.flag"))

summary_rows <- do.call(rbind, lapply(names(panels), function(panel_name) {
  keys <- panels[[panel_name]]
  do.call(rbind, lapply(estimators, function(nm) {
    values <- per_unit$lag1[per_unit$estimator == nm &
                              per_unit$pair_key %in% keys]
    data.frame(panel = panel_name, estimator = nm, n = length(values),
               mean_lag1 = mean(values), median_lag1 = median(values),
               sd_lag1 = sd(values), stringsAsFactors = FALSE)
  }))
}))
print(summary_rows, row.names = FALSE)
message("common C: lag-1 ", round(lag_profile(common_c)[1L], 4),
        ", mean off-diagonal ", round(mean_off_diagonal(common_c), 4))
message("Design factor lag-1: mean ", round(mean(design_lag1), 4),
        ", range [", round(min(design_lag1), 3), ", ",
        round(max(design_lag1), 3), "]")
