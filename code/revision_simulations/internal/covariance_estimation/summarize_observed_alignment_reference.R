#!/usr/bin/env Rscript

# Refit the observed genotype-expression alignment for pilot comparison.

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
  workflowr_root,
  "code", "revision_simulations", "r4_correlated_errors",
  "real_data_correlation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "covariance_estimation",
  "donor_null_permutation_helpers.R"
))

output_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "donor_null_permutation_pilot"
)
input_path <- file.path(output_dir, "input", "selected_raw_data.rds")
if (!file.exists(input_path)) {
  stop("The selected raw-data pilot cache does not exist: ", input_path)
}
input <- readRDS(input_path)

n_units <- nrow(input$unit_table)
n_times <- length(input$time_grid)
beta_hat <- matrix(NA_real_, nrow = n_units, ncol = n_times)
raw_se <- matrix(NA_real_, nrow = n_units, ncol = n_times)
residual_df <- integer(n_times)
for (time_index in seq_along(input$time_inputs)) {
  time_input <- input$time_inputs[[time_index]]
  genotype <- input$unit_dosage[
    time_input$donors,
    ,
    drop = FALSE
  ]
  fit <- fit_residualized_genotype_regressions(
    time_input$expression_residual,
    genotype,
    time_input$projection$residualizer,
    time_input$projection$rank
  )
  beta_hat[, time_index] <- fit$beta
  raw_se[, time_index] <- fit$standard_error
  residual_df[time_index] <- fit$residual_df
}
adjusted_se <- convert_raw_to_t_adjusted_se(beta_hat, raw_se, residual_df)

estimator_functions <- list(
  ordinary_mean = estimate_ordinary_pairwise_correlation,
  ols_weighted = estimate_pairwise_difference_correlation
)
se_scales <- list(
  raw_regression = raw_se,
  t_adjusted = adjusted_se
)
matrices <- list()
diagnostic_rows <- list()
lag_rows <- list()
output_index <- 1L
for (threshold_index in seq_along(input$threshold_values)) {
  selected_index <- input$threshold_indices[[threshold_index]]
  for (se_scale in names(se_scales)) {
    for (estimator in names(estimator_functions)) {
      correlation <- estimator_functions[[estimator]](
        beta_hat[selected_index, , drop = FALSE],
        se_scales[[se_scale]][selected_index, , drop = FALSE]
      )
      key <- paste(
        input$threshold_values[threshold_index],
        se_scale,
        estimator,
        sep = "__"
      )
      matrices[[key]] <- correlation
      identifiers <- data.frame(
        reference = "observed_alignment",
        threshold = input$threshold_values[threshold_index],
        n_units = length(selected_index),
        se_scale = se_scale,
        estimator = estimator,
        stringsAsFactors = FALSE
      )
      diagnostic_rows[[output_index]] <- cbind(
        identifiers,
        summarize_raw_correlation_matrix(correlation)
      )
      lag_correlation <- lag_average_correlation(correlation)
      lag_rows[[output_index]] <- cbind(
        identifiers[rep(1L, length(lag_correlation)), , drop = FALSE],
        data.frame(
          lag = seq_along(lag_correlation),
          correlation = lag_correlation,
          semivariogram = 1 - lag_correlation
        )
      )
      output_index <- output_index + 1L
    }
  }
}
diagnostics <- do.call(rbind, diagnostic_rows)
lag_summaries <- do.call(rbind, lag_rows)
rownames(diagnostics) <- NULL
rownames(lag_summaries) <- NULL

reference <- list(
  description = paste(
    "Observed, unpermuted time-specific Y ~ 1 + G + PC1 + ... + PC5",
    "refit on the same selected units and raw inputs used by the pilot"
  ),
  input_cache = input_path,
  input_cache_md5 = unname(tools::md5sum(input_path)),
  beta_hat = beta_hat,
  raw_se = raw_se,
  adjusted_se = adjusted_se,
  residual_df = residual_df,
  matrices = matrices,
  diagnostics = diagnostics,
  lag_summaries = lag_summaries,
  r_version = R.version.string
)
saveRDS(
  reference,
  file.path(output_dir, "observed_alignment_reference.rds"),
  compress = "xz"
)
utils::write.csv(
  diagnostics,
  file.path(output_dir, "summary", "observed_alignment_key_contrasts.csv"),
  row.names = FALSE
)
utils::write.csv(
  lag_summaries,
  file.path(output_dir, "summary", "observed_alignment_lag_summaries.csv"),
  row.names = FALSE
)

print(diagnostics, row.names = FALSE)
