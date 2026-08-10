#!/usr/bin/env Rscript

# Run exact-design null controls and synchronized donor-label permutations.

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

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

workflowr_root <- find_workflowr_root()
project_root <- normalizePath(file.path(workflowr_root, ".."), mustWork = TRUE)
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "r4_correlated_errors",
  "real_data_correlation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "covariance_estimation",
  "donor_null_permutation_helpers.R"
))

n_replications <- as.integer(get_arg("--n-reps", "100"))
seed_start <- as.integer(get_arg("--seed-start", "20260831"))
output_id <- get_arg("--output-id", "donor_null_permutation_pilot")
if (is.na(n_replications) || n_replications < 2L ||
    is.na(seed_start) || !nzchar(output_id)) {
  stop("Invalid replication count, seed, or output ID.")
}

selection_cache_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "pairwise_lfdr_threshold_correlation",
  "pairwise_lfdr_threshold_analysis.rds"
)
vcf_path <- file.path(
  project_root,
  "iPSC-data", "genotype-data", "YRI_genotype.vcf.gz"
)
expression_path <- file.path(
  project_root,
  "iPSC-data", "expression-data",
  "quantile_normalized_no_projection.txt"
)
pc_path <- file.path(
  workflowr_root,
  "data", "dynamic_eQTL_real", "principal_components_10.txt"
)
required_inputs <- c(selection_cache_path, vcf_path, expression_path, pc_path)
if (any(!file.exists(required_inputs))) {
  stop("At least one required pilot input is missing: ",
       paste(required_inputs[!file.exists(required_inputs)], collapse = ", "))
}

output_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal", output_id
)
input_dir <- file.path(output_dir, "input")
summary_dir <- file.path(output_dir, "summary")
invisible(lapply(
  c(output_dir, input_dir, summary_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

message("Loading and validating the fixed high-lfdr selections.")
selection_cache <- readRDS(selection_cache_path)
expected_selection_names <- c("lfdr_0p97", "lfdr_0p96", "lfdr_0p95")
expected_thresholds <- c(0.97, 0.96, 0.95)
expected_counts <- c(1046L, 2323L, 4025L)
if (!identical(names(selection_cache$selections), expected_selection_names)) {
  stop("The selection cache does not contain the expected threshold sets.")
}
for (selection_index in seq_along(expected_selection_names)) {
  selection <- selection_cache$selections[[selection_index]]
  if (!isTRUE(all.equal(selection$threshold, expected_thresholds[selection_index])) ||
      nrow(selection$selected) != expected_counts[selection_index] ||
      anyDuplicated(selection$selected$gene_id) ||
      any(selection$selected$lfdr <= expected_thresholds[selection_index])) {
    stop("A fixed lfdr selection failed validation.")
  }
}
unit_table <- selection_cache$selections$lfdr_0p95$selected
unit_table$unit_index <- seq_len(nrow(unit_table))
threshold_indices <- lapply(selection_cache$selections, function(selection) {
  index <- match(selection$selected$pair_key, unit_table$pair_key)
  if (anyNA(index)) {
    stop("A tighter threshold set is not nested in the 0.95 set.")
  }
  index
})

message("Reading selected expression rows and time-specific PCs.")
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
    !all(c("Sample_id", paste0("PC", 1:5)) %in% names(pc_data)) ||
    anyDuplicated(expression_data$Gene_id) || anyDuplicated(pc_data$Sample_id)) {
  stop("The expression or PC input has invalid identifiers.")
}
expression_sample_ids <- names(expression_data)[-1L]
if (length(expression_sample_ids) != 297L || nrow(pc_data) != 297L ||
    !setequal(expression_sample_ids, pc_data$Sample_id)) {
  stop("Expression and PC donor-time IDs do not match exactly.")
}
gene_row <- match(unit_table$gene_id, expression_data$Gene_id)
if (anyNA(gene_row)) {
  stop("At least one selected gene is missing from the expression matrix.")
}

message("Streaming selected dosages from the raw VCF.")
unique_variant_ids <- unique(unit_table$variant_id)
unique_dosage <- read_selected_vcf_dosages(
  vcf_path,
  unique_variant_ids,
  chunk_size = 100000L
)
unit_dosage <- unique_dosage[, match(unit_table$variant_id, unique_variant_ids),
                              drop = FALSE]
colnames(unit_dosage) <- unit_table$pair_key
vcf_donors <- rownames(unit_dosage)
if (nrow(unit_dosage) != 19L || anyDuplicated(vcf_donors)) {
  stop("The selected dosage matrix does not contain 19 unique donors.")
}

time_grid <- 0:15
expected_sample_counts <- c(19L, 19L, 16L, 19L, 16L, 19L, 19L, 19L,
                            19L, 19L, 19L, 19L, 19L, 18L, 19L, 19L)
time_inputs <- vector("list", length(time_grid))
names(time_inputs) <- paste0("time_", time_grid)
for (time_index in seq_along(time_grid)) {
  time_value <- time_grid[time_index]
  sample_ids <- grep(
    paste0("_", time_value, "$"),
    expression_sample_ids,
    value = TRUE
  )
  donors <- sub(paste0("_", time_value, "$"), "", sample_ids)
  if (length(sample_ids) != expected_sample_counts[time_index] ||
      anyDuplicated(donors) || any(!donors %in% vcf_donors)) {
    stop("Unexpected donor coverage at time ", time_value, ".")
  }
  pc_rows <- match(sample_ids, pc_data$Sample_id)
  expression_columns <- match(sample_ids, names(expression_data))
  expression_matrix <- t(as.matrix(
    expression_data[gene_row, expression_columns, drop = FALSE]
  ))
  storage.mode(expression_matrix) <- "double"
  rownames(expression_matrix) <- donors
  colnames(expression_matrix) <- unit_table$pair_key
  covariates <- as.matrix(pc_data[pc_rows, paste0("PC", 1:5), drop = FALSE])
  storage.mode(covariates) <- "double"
  rownames(covariates) <- donors
  projection <- make_covariate_projection(covariates)
  expression_residual <- projection$residualizer %*% expression_matrix
  rownames(expression_residual) <- donors
  colnames(expression_residual) <- unit_table$pair_key
  covariate_residual_df <- nrow(expression_matrix) - projection$rank
  observed_residual_sd <- sqrt(
    colSums(expression_residual^2) / covariate_residual_df
  )
  if (any(!is.finite(observed_residual_sd)) ||
      any(observed_residual_sd <= 0)) {
    stop("The observed covariate residual SD is invalid at time ", time_value, ".")
  }
  time_inputs[[time_index]] <- list(
    time = time_value,
    sample_ids = sample_ids,
    donors = donors,
    expression = expression_matrix,
    expression_residual = expression_residual,
    covariates = covariates,
    projection = projection,
    covariate_residual_df = covariate_residual_df,
    observed_residual_sd = observed_residual_sd
  )
}

donor_observation_matrix <- vapply(time_inputs, function(input) {
  vcf_donors %in% input$donors
}, logical(length(vcf_donors)))
rownames(donor_observation_matrix) <- vcf_donors
colnames(donor_observation_matrix) <- paste0("time_", time_grid)
donor_observation_patterns <- apply(
  donor_observation_matrix,
  1L,
  paste0,
  collapse = ""
)

input_coverage <- data.frame(
  input = c(
    "Selected units at lfdr > 0.97",
    "Selected units at lfdr > 0.96",
    "Selected units at lfdr > 0.95",
    "Unique selected genes",
    "Unique selected variants",
    "Expression donor-time columns",
    "PC donor-time rows",
    "VCF donors"
  ),
  count = c(
    expected_counts,
    length(unique(unit_table$gene_id)),
    length(unique_variant_ids),
    length(expression_sample_ids),
    nrow(pc_data),
    length(vcf_donors)
  ),
  stringsAsFactors = FALSE
)

source_information <- data.frame(
  role = c("selection_cache", "genotype_vcf", "expression", "pc_data"),
  path = required_inputs,
  size_bytes = unname(file.info(required_inputs)$size),
  mtime = format(file.info(required_inputs)$mtime, tz = "UTC", usetz = TRUE),
  md5 = unname(tools::md5sum(required_inputs)),
  stringsAsFactors = FALSE
)

selected_raw_data <- list(
  unit_table = unit_table,
  threshold_indices = threshold_indices,
  threshold_values = expected_thresholds,
  time_grid = time_grid,
  expected_sample_counts = expected_sample_counts,
  vcf_donors = vcf_donors,
  donor_observation_matrix = donor_observation_matrix,
  donor_observation_patterns = donor_observation_patterns,
  unit_dosage = unit_dosage,
  time_inputs = time_inputs,
  source_information = source_information,
  input_coverage = input_coverage
)
saveRDS(
  selected_raw_data,
  file.path(input_dir, "selected_raw_data.rds"),
  compress = "xz"
)
write_csv(input_coverage, file.path(summary_dir, "input_coverage.csv"))
write_csv(source_information, file.path(summary_dir, "source_information.csv"))

fit_one_replication <- function(generator,
                                replication_seed) {
  set.seed(replication_seed)
  beta_hat <- matrix(
    NA_real_,
    nrow = nrow(unit_table),
    ncol = length(time_grid),
    dimnames = list(unit_table$pair_key, paste0("time_", time_grid))
  )
  raw_se <- beta_hat
  oracle_se <- beta_hat
  residual_df <- integer(length(time_grid))
  if (identical(generator, "donor_label_permutation")) {
    donor_permutation <- sample(vcf_donors, length(vcf_donors), replace = FALSE)
    names(donor_permutation) <- vcf_donors
  } else if (identical(generator, "donor_residual_block_permutation")) {
    donor_permutation <- make_shared_donor_block_permutation(
      vcf_donors,
      donor_observation_patterns
    )
  } else if (!identical(generator, "independent_gaussian")) {
    stop("Unknown null generator: ", generator)
  }

  for (time_index in seq_along(time_grid)) {
    input <- time_inputs[[time_index]]
    if (identical(generator, "independent_gaussian")) {
      simulated_expression <- matrix(
        stats::rnorm(length(input$expression)),
        nrow = nrow(input$expression),
        ncol = ncol(input$expression)
      )
      simulated_expression <- sweep(
        simulated_expression,
        2L,
        input$observed_residual_sd,
        `*`
      )
      expression_residual <- input$projection$residualizer %*%
        simulated_expression
      genotype <- unit_dosage[input$donors, , drop = FALSE]
    } else if (identical(generator, "donor_label_permutation")) {
      expression_residual <- input$expression_residual
      permuted_donors <- unname(donor_permutation[input$donors])
      genotype <- unit_dosage[permuted_donors, , drop = FALSE]
    } else {
      source_donors <- unname(donor_permutation[input$donors])
      source_rows <- match(source_donors, input$donors)
      if (anyNA(source_rows)) {
        stop("A residual donor block is unavailable at time ", input$time, ".")
      }
      permuted_residual <- input$expression_residual[
        source_rows,
        ,
        drop = FALSE
      ]
      expression_residual <- input$projection$residualizer %*%
        permuted_residual
      genotype <- unit_dosage[input$donors, , drop = FALSE]
    }
    if (identical(generator, "independent_gaussian")) {
      genotype_residual <- input$projection$residualizer %*% genotype
      oracle_se[, time_index] <- input$observed_residual_sd /
        sqrt(colSums(genotype_residual^2))
    }
    fit <- fit_residualized_genotype_regressions(
      expression_residual,
      genotype,
      input$projection$residualizer,
      input$projection$rank
    )
    beta_hat[, time_index] <- fit$beta
    raw_se[, time_index] <- fit$standard_error
    residual_df[time_index] <- fit$residual_df
  }
  list(
    beta_hat = beta_hat,
    raw_se = raw_se,
    adjusted_se = convert_raw_to_t_adjusted_se(
      beta_hat,
      raw_se,
      residual_df
    ),
    oracle_se = if (identical(generator, "independent_gaussian")) {
      oracle_se
    } else {
      NULL
    },
    residual_df = residual_df
  )
}

estimator_functions <- list(
  ordinary_mean = estimate_ordinary_pairwise_correlation,
  ols_weighted = estimate_pairwise_difference_correlation
)
generator_names <- c(
  "independent_gaussian",
  "donor_residual_block_permutation",
  "donor_label_permutation"
)
se_scales_by_generator <- list(
  independent_gaussian = c(
    "oracle_sampling", "raw_regression", "t_adjusted"
  ),
  donor_residual_block_permutation = c("raw_regression", "t_adjusted"),
  donor_label_permutation = c("raw_regression", "t_adjusted")
)
matrix_key_grid <- do.call(rbind, lapply(generator_names, function(generator) {
  expand.grid(
    generator = generator,
    selection = expected_selection_names,
    se_scale = se_scales_by_generator[[generator]],
    estimator = names(estimator_functions),
    stringsAsFactors = FALSE
  )
}))
matrix_keys <- do.call(
  paste,
  c(matrix_key_grid, list(sep = "__"))
)
matrix_draws <- setNames(lapply(matrix_keys, function(unused) {
  array(
    NA_real_,
    dim = c(length(time_grid), length(time_grid), n_replications),
    dimnames = list(
      paste0("time_", time_grid),
      paste0("time_", time_grid),
      paste0("rep_", seq_len(n_replications))
    )
  )
}), matrix_keys)
lag_rows <- vector(
  "list",
  n_replications * length(expected_thresholds) *
    length(estimator_functions) *
    sum(lengths(se_scales_by_generator))
)
diagnostic_rows <- lag_rows
output_index <- 1L

iid_seeds <- seed_start + seq_len(n_replications) - 1L
permutation_seeds <- seed_start + 100000L + seq_len(n_replications) - 1L
residual_permutation_seeds <- seed_start + 200000L +
  seq_len(n_replications) - 1L
seed_sets <- list(
  independent_gaussian = iid_seeds,
  donor_residual_block_permutation = residual_permutation_seeds,
  donor_label_permutation = permutation_seeds
)

message(
  "Running ",
  n_replications,
  " independent-Gaussian controls and ",
  n_replications,
  " synchronized residual-block permutations, and ",
  n_replications,
  " donor-label permutations."
)
for (generator in generator_names) {
  for (replication in seq_len(n_replications)) {
    fit <- fit_one_replication(
      generator,
      seed_sets[[generator]][replication]
    )
    for (threshold_index in seq_along(expected_thresholds)) {
      selected_index <- threshold_indices[[threshold_index]]
      for (se_scale in se_scales_by_generator[[generator]]) {
        selected_se <- switch(
          se_scale,
          oracle_sampling = fit$oracle_se[selected_index, , drop = FALSE],
          raw_regression = fit$raw_se[selected_index, , drop = FALSE],
          t_adjusted = fit$adjusted_se[selected_index, , drop = FALSE],
          stop("Unknown SE scale: ", se_scale)
        )
        selected_beta <- fit$beta_hat[selected_index, , drop = FALSE]
        for (estimator_name in names(estimator_functions)) {
          correlation <- estimator_functions[[estimator_name]](
            selected_beta,
            selected_se
          )
          key <- paste(
            generator,
            expected_selection_names[threshold_index],
            se_scale,
            estimator_name,
            sep = "__"
          )
          matrix_draws[[key]][, , replication] <- correlation
          lag_correlation <- lag_average_correlation(correlation)
          lag_rows[[output_index]] <- data.frame(
            generator = generator,
            replication = replication,
            seed = seed_sets[[generator]][replication],
            threshold = expected_thresholds[threshold_index],
            n_units = length(selected_index),
            se_scale = se_scale,
            estimator = estimator_name,
            lag = seq_along(lag_correlation),
            correlation = lag_correlation,
            semivariogram = 1 - lag_correlation,
            stringsAsFactors = FALSE
          )
          diagnostics <- summarize_raw_correlation_matrix(correlation)
          diagnostic_rows[[output_index]] <- cbind(
            data.frame(
              generator = generator,
              replication = replication,
              seed = seed_sets[[generator]][replication],
              threshold = expected_thresholds[threshold_index],
              n_units = length(selected_index),
              se_scale = se_scale,
              estimator = estimator_name,
              stringsAsFactors = FALSE
            ),
            diagnostics
          )
          output_index <- output_index + 1L
        }
      }
    }
    if (replication %% max(1L, floor(n_replications / 10L)) == 0L) {
      message(
        generator,
        ": completed ",
        replication,
        " of ",
        n_replications,
        " replications."
      )
    }
  }
}

replicate_lag_summaries <- do.call(rbind, lag_rows)
replicate_matrix_diagnostics <- do.call(rbind, diagnostic_rows)
rownames(replicate_lag_summaries) <- NULL
rownames(replicate_matrix_diagnostics) <- NULL

aggregate_values <- function(data,
                             value_name,
                             grouping_names) {
  split_key <- interaction(data[grouping_names], drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(data)), split_key)
  rows <- lapply(groups, function(index) {
    values <- data[[value_name]][index]
    identifiers <- data[index[1L], grouping_names, drop = FALSE]
    cbind(
      identifiers,
      data.frame(
        metric = value_name,
        mean = mean(values),
        median = stats::median(values),
        lower = unname(stats::quantile(values, 0.025, names = FALSE)),
        upper = unname(stats::quantile(values, 0.975, names = FALSE)),
        sd = stats::sd(values),
        n_replications = length(values),
        stringsAsFactors = FALSE
      )
    )
  })
  output <- do.call(rbind, rows)
  rownames(output) <- NULL
  output
}

lag_grouping <- c(
  "generator", "threshold", "n_units", "se_scale", "estimator", "lag"
)
aggregate_correlation <- aggregate_values(
  replicate_lag_summaries,
  "correlation",
  lag_grouping
)
aggregate_semivariogram <- aggregate_values(
  replicate_lag_summaries,
  "semivariogram",
  lag_grouping
)
aggregate_lag_summaries <- rbind(
  aggregate_correlation,
  aggregate_semivariogram
)
rownames(aggregate_lag_summaries) <- NULL

diagnostic_metrics <- c(
  "lag1", "mean_lags_9_15", "lag15", "mean_off_diagonal",
  "minimum_off_diagonal", "maximum_off_diagonal", "minimum_eigenvalue",
  "n_negative_eigenvalues"
)
diagnostic_grouping <- c(
  "generator", "threshold", "n_units", "se_scale", "estimator"
)
aggregate_matrix_diagnostics <- do.call(rbind, lapply(
  diagnostic_metrics,
  function(metric) {
    aggregate_values(
      replicate_matrix_diagnostics,
      metric,
      diagnostic_grouping
    )
  }
))
rownames(aggregate_matrix_diagnostics) <- NULL

key_contrasts <- aggregate_matrix_diagnostics[
  aggregate_matrix_diagnostics$metric %in%
    c("lag1", "mean_lags_9_15", "lag15", "mean_off_diagonal"),
]
key_contrasts <- key_contrasts[order(
  key_contrasts$generator,
  -key_contrasts$threshold,
  key_contrasts$se_scale,
  key_contrasts$estimator,
  match(key_contrasts$metric,
        c("lag1", "mean_lags_9_15", "lag15", "mean_off_diagonal"))
), ]
rownames(key_contrasts) <- NULL

configuration <- list(
  analysis_id = output_id,
  created_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  pilot_status = "Internal synchronized donor-null permutation pilot",
  thresholds = expected_thresholds,
  selected_counts = expected_counts,
  selection = paste(
    "Fixed one-highest-BF-adjusted-lfdr variant per gene before",
    "null generation; use all units with lfdr strictly above threshold"
  ),
  regression = "Time-specific Y ~ 1 + G + PC1 + ... + PC5",
  pc_source = pc_path,
  pc_version = "Timepoint-specific PCs used by code/00_eQTLs.R",
  null_generators = c(
    independent_gaussian = paste(
      "Independent Gaussian donor-by-unit errors at every time point,",
      "scaled by observed covariate-only residual SD"
    ),
    donor_label_permutation = paste(
      "One global permutation of 19 genotype donor labels per replicate,",
      "shared by every unit and all time points; expression, PCs, and",
      "missingness remain fixed"
    ),
    donor_residual_block_permutation = paste(
      "Covariate-only residual donor blocks are permuted within the four",
      "observed missingness-pattern strata using one shared mapping for",
      "every unit and time point; residuals are reprojected against each",
      "time-specific PC design before the genotype regression"
    )
  ),
  n_replications = n_replications,
  iid_seeds = iid_seeds,
  permutation_seeds = permutation_seeds,
  residual_permutation_seeds = residual_permutation_seeds,
  donor_observation_pattern_counts = table(donor_observation_patterns),
  se_scales = c(
    "Oracle sampling SD available only in the independent-Gaussian control",
    "Estimated raw regression SE",
    "Conservative t-to-normal adjusted SE used by FASH"
  ),
  estimators = c(
    "Ordinary mean of unit-specific pairwise moment estimates",
    "OLS-through-origin x-squared-weighted pairwise moment estimate"
  ),
  matrix_version = "Raw unprojected pairwise estimates",
  source_information = source_information,
  r_version = R.version.string
)

analysis_result <- list(
  configuration = configuration,
  input_coverage = input_coverage,
  replicate_lag_summaries = replicate_lag_summaries,
  replicate_matrix_diagnostics = replicate_matrix_diagnostics,
  aggregate_lag_summaries = aggregate_lag_summaries,
  aggregate_matrix_diagnostics = aggregate_matrix_diagnostics,
  key_contrasts = key_contrasts,
  matrix_draws = matrix_draws
)
saveRDS(configuration, file.path(output_dir, "configuration.rds"))
saveRDS(
  analysis_result,
  file.path(output_dir, "donor_null_permutation_analysis.rds"),
  compress = "xz"
)
write_csv(
  replicate_lag_summaries,
  file.path(summary_dir, "replicate_lag_summaries.csv")
)
write_csv(
  aggregate_lag_summaries,
  file.path(summary_dir, "aggregate_lag_summaries.csv")
)
write_csv(
  replicate_matrix_diagnostics,
  file.path(summary_dir, "replicate_matrix_diagnostics.csv")
)
write_csv(
  aggregate_matrix_diagnostics,
  file.path(summary_dir, "aggregate_matrix_diagnostics.csv")
)
write_csv(key_contrasts, file.path(summary_dir, "key_contrasts.csv"))

cat("\nKey pilot contrasts:\n")
print(
  key_contrasts[
    key_contrasts$threshold == 0.95 &
      key_contrasts$estimator == "ols_weighted" &
      key_contrasts$metric %in% c("lag1", "mean_lags_9_15"),
  ],
  row.names = FALSE
)
cat("\nSaved donor null-permutation pilot to:\n", output_dir, "\n", sep = "")
