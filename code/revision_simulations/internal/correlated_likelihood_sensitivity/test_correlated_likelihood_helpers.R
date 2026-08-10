#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "correlated_likelihood_sensitivity",
  "correlated_likelihood_helpers.R"
))

set.seed(20260851)
n_units <- 24L
n_time <- 6L
pair_keys <- sprintf("gene_%03d_variant", seq_len(n_units))
time_names <- paste0("time_", 0:(n_time - 1L))
adjusted_se <- matrix(
  stats::runif(n_units * n_time, 0.55, 1.45),
  nrow = n_units,
  dimnames = list(pair_keys, time_names)
)
correlation <- matrix(0.25, nrow = n_time, ncol = n_time)
diag(correlation) <- 1
dimnames(correlation) <- list(time_names, time_names)
precision <- construct_unit_precision_matrices(adjusted_se, correlation)

for (index in c(1L, 7L, n_units)) {
  marginal_scale <- diag(adjusted_se[index, ])
  expected_covariance <- marginal_scale %*% correlation %*% marginal_scale
  observed_covariance <- solve(precision$precision_matrices[[index]])
  stopifnot(
    max(abs(observed_covariance - expected_covariance)) < 1e-10,
    max(abs(diag(observed_covariance) - adjusted_se[index, ]^2)) < 1e-10,
    max(abs(stats::cov2cor(observed_covariance) - correlation)) < 1e-10
  )
}

non_pd <- correlation
non_pd[1, 2] <- non_pd[2, 1] <- 1.2
non_pd_error <- try(
  validate_shared_correlation(non_pd, n_time = n_time),
  silent = TRUE
)
stopifnot(inherits(non_pd_error, "try-error"))

beta_hat <- matrix(
  stats::rnorm(n_units * n_time, sd = 0.7),
  nrow = n_units,
  dimnames = list(pair_keys, time_names)
)
time_grid <- 0:(n_time - 1L)
data_list <- make_fash_data_list(beta_hat, time_grid)
data_list_with_se <- Map(function(dataset, se) {
  dataset$SE <- as.numeric(se)
  dataset
}, data_list, split(adjusted_se, row(adjusted_se)))
names(data_list_with_se) <- pair_keys
settings <- list(
  num_basis = 8L,
  order = 1L,
  betaprec = 1e-6,
  pred_step = 1,
  likelihood = "gaussian",
  penalty = 1L
)
psd_grid <- c(0, 0.05, 0.15)

standard_error_fit <- fashr::fash(
  Y = "beta",
  smooth_var = "time",
  S = "SE",
  data_list = data_list_with_se,
  grid = psd_grid,
  likelihood = settings$likelihood,
  num_basis = settings$num_basis,
  betaprec = settings$betaprec,
  order = settings$order,
  pred_step = settings$pred_step,
  penalty = settings$penalty,
  num_cores = 1L,
  verbose = FALSE
)
identity_result <- fit_fash_with_shared_correlation(
  beta_hat = beta_hat,
  adjusted_se = adjusted_se,
  time_grid = time_grid,
  correlation = diag(n_time),
  settings = settings,
  psd_grid = psd_grid,
  num_cores = 1L,
  verbose = FALSE
)
identity_validation <- validate_identity_path_equivalence(
  standard_error_fit,
  identity_result$fit,
  likelihood_tolerance = 1e-7,
  prior_tolerance = 1e-6,
  lfdr_tolerance = 1e-6
)

method_fits <- list(
  diagonal = standard_error_fit,
  identity_precision = identity_result$fit
)
method_labels <- c("Diagonal SE", "Identity precision")
prior_comparison <- compare_prior_weight_fits(
  method_fits,
  fit_stage = "Raw",
  method_labels = method_labels
)
pair_metadata <- data.frame(
  pair_key = pair_keys,
  gene_id = sub("_variant$", "", pair_keys),
  variant_id = "variant",
  stringsAsFactors = FALSE
)
lfdr_comparison <- compare_lfdr_fits(
  method_fits,
  pair_metadata = pair_metadata,
  fit_stage = "Raw",
  method_labels = method_labels,
  top_n = 5L
)
bf_updates <- run_bf_updates_checked(
  method_fits,
  method_labels = method_labels,
  pair_keys = pair_keys
)

all_null_fit <- standard_error_fit
all_null_fit$L_matrix[, 1L] <- 0
all_null_fit$L_matrix[, -1L] <- -100
all_null_eb <- fashr::fash_eb_est(
  all_null_fit$L_matrix,
  grid = all_null_fit$psd_grid,
  penalty = 1L
)
all_null_fit$prior_weights <- all_null_eb$prior_weight
all_null_fit$posterior_weights <- all_null_eb$posterior_weight
rownames(all_null_fit$posterior_weights) <- pair_keys
all_null_fit$lfdr <- rep(1, n_units)
names(all_null_fit$lfdr) <- pair_keys
all_null_update <- run_bf_updates_checked(
  list(all_null = all_null_fit),
  method_labels = "All-null test",
  pair_keys = pair_keys
)

stopifnot(
  identical(names(precision$precision_matrices), pair_keys),
  precision$correlation_diagnostics$minimum_eigenvalue > 0,
  identity_validation$row_centered_likelihood_maximum_difference < 1e-7,
  nrow(prior_comparison$pairwise_metrics) == 1L,
  prior_comparison$pairwise_metrics$prior_total_variation < 1e-6,
  nrow(lfdr_comparison$lfdr_wide) == n_units,
  nrow(lfdr_comparison$lfdr_long) == 2L * n_units,
  nrow(lfdr_comparison$pairwise_metrics) == 1L,
  lfdr_comparison$pairwise_metrics$maximum_absolute_lfdr_difference < 1e-6,
  nrow(lfdr_comparison$top_discrepancies) == 5L,
  all(bf_updates$status$bf_update_available),
  length(bf_updates$successful_fits) == 2L,
  !all_null_update$status$bf_update_available,
  length(all_null_update$successful_fits) == 0L,
  grepl("conditional alternative mixture", all_null_update$status$bf_update_status)
)

cat("Correlated-likelihood helper tests passed.\n")
