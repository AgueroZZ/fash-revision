#!/usr/bin/env Rscript

# Focused tests for the FASH-style finite mixture of Gaussian linear slopes.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the shared revision simulation functions.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "shared",
  "simulation_functions.R"
))

time_grid <- 0:4
test_grid <- c(0, 0.05, 0.1, 0.25, 0.5, 1)
test_penalty <- 10L
test_pred_step <- 1
n_units <- 18L
intercepts <- seq(-0.45, 0.4, length.out = n_units)
slopes <- c(rep(0, 12L), 0.08, -0.12, 0.25, -0.35, 0.7, -0.9)

datasets <- lapply(seq_len(n_units), function(index) {
  standard_error <- 0.12 + 0.01 * ((seq_along(time_grid) + index) %% 4)
  deterministic_residual <- 0.015 * sin(
    time_grid + index * seq_along(time_grid) / 5
  )
  list(
    x = time_grid,
    y = intercepts[index] + slopes[index] * time_grid +
      deterministic_residual,
    sd = standard_error
  )
})
names(datasets) <- sprintf("unit_%02d", seq_len(n_units))

fit <- fit_linear_mixture_fash(
  datasets = datasets,
  grid = test_grid,
  pred_step = test_pred_step,
  penalty = test_penalty
)

expected_likelihood <- vapply(test_grid, function(tau) {
  vapply(datasets, function(dataset) {
    normalized_time <- (
      dataset$x - min(dataset$x)
    ) / test_pred_step
    log_marginal_common_intercept(
      y = dataset$y,
      se = dataset$sd,
      x = normalized_time,
      slope_sd = tau
    )
  }, numeric(1))
}, numeric(length(datasets)))
dimnames(expected_likelihood) <- list(names(datasets), as.character(test_grid))

precomputed_fit <- fit_linear_mixture_fash_from_log_likelihood(
  L_matrix = expected_likelihood,
  grid = test_grid,
  pred_step = test_pred_step,
  penalty = test_penalty
)

stopifnot(
  inherits(fit, "linear_mixture_fash"),
  identical(fit$psd_grid, test_grid),
  identical(fit$settings$pred_step, test_pred_step),
  identical(fit$settings$penalty, test_penalty),
  identical(fit$bf_adjusted, FALSE),
  identical(dim(fit$L_matrix), c(n_units, length(test_grid))),
  identical(dimnames(fit$L_matrix), dimnames(expected_likelihood)),
  max(abs(fit$L_matrix - expected_likelihood)) < 1e-10,
  any(fit$prior_weights$psd == 0),
  fit$prior_weights$psd[1] == 0,
  colnames(fit$posterior_weights)[1] == "0",
  identical(
    colnames(fit$posterior_weights),
    as.character(fit$prior_weights$psd)
  ),
  abs(sum(fit$prior_weights$prior_weight) - 1) < 1e-10,
  all(abs(rowSums(fit$posterior_weights) - 1) < 1e-10),
  all(is.finite(fit$lfdr)),
  all(fit$lfdr >= 0 & fit$lfdr <= 1),
  identical(precomputed_fit$L_matrix, fit$L_matrix),
  identical(precomputed_fit$settings, fit$settings),
  max(abs(
    expand_grid_prior_weights(precomputed_fit$prior_weights, test_grid) -
      expand_grid_prior_weights(fit$prior_weights, test_grid)
  )) < 1e-8,
  max(abs(precomputed_fit$posterior_weights - fit$posterior_weights)) < 1e-8,
  max(abs(precomputed_fit$lfdr - fit$lfdr)) < 1e-8,
  isTRUE(validate_linear_mixture_fash(
    fit,
    expected_grid = test_grid,
    expected_pred_step = test_pred_step,
    expected_penalty = test_penalty
  ))
)

rescaled_datasets <- lapply(datasets, function(dataset) {
  dataset$x <- 3 * dataset$x
  dataset
})
rescaled_fit <- fit_linear_mixture_fash(
  datasets = rescaled_datasets,
  grid = test_grid,
  pred_step = 3 * test_pred_step,
  penalty = test_penalty
)
stopifnot(
  max(abs(rescaled_fit$L_matrix - fit$L_matrix)) < 1e-10,
  max(abs(
    expand_grid_prior_weights(rescaled_fit$prior_weights, test_grid) -
      expand_grid_prior_weights(fit$prior_weights, test_grid)
  )) < 1e-8,
  max(abs(rescaled_fit$lfdr - fit$lfdr)) < 1e-8
)

direct_eb <- fashr::fash_eb_est(
  L_matrix = fit$L_matrix,
  grid = test_grid,
  penalty = test_penalty
)
stopifnot(
  max(abs(
    expand_grid_prior_weights(fit$prior_weights, test_grid) -
      expand_grid_prior_weights(direct_eb$prior_weight, test_grid)
  )) < 1e-8,
  max(abs(fit$posterior_weights - direct_eb$posterior_weight)) < 1e-8
)

fit_bf <- BF_update_linear_mixture_fash(fit)
precomputed_fit_bf <- BF_update_linear_mixture_fash(precomputed_fit)
direct_bf <- fashr::BF_update(fit, plot = FALSE)
stopifnot(
  inherits(fit_bf, "linear_mixture_fash"),
  isTRUE(fit_bf$bf_adjusted),
  identical(fit_bf$psd_grid, fit$psd_grid),
  identical(fit_bf$L_matrix, fit$L_matrix),
  identical(fit_bf$settings, fit$settings),
  identical(fit_bf$eb_result, fit$eb_result),
  max(abs(
    expand_grid_prior_weights(fit_bf$prior_weights, test_grid) -
      expand_grid_prior_weights(direct_bf$prior_weights, test_grid)
  )) < 1e-8,
  max(abs(fit_bf$posterior_weights - direct_bf$posterior_weights)) < 1e-8,
  max(abs(fit_bf$lfdr - direct_bf$lfdr)) < 1e-8,
  max(abs(
    expand_grid_prior_weights(precomputed_fit_bf$prior_weights, test_grid) -
      expand_grid_prior_weights(fit_bf$prior_weights, test_grid)
  )) < 1e-8,
  max(abs(
    precomputed_fit_bf$posterior_weights - fit_bf$posterior_weights
  )) < 1e-8,
  max(abs(precomputed_fit_bf$lfdr - fit_bf$lfdr)) < 1e-8,
  isTRUE(all.equal(
    precomputed_fit_bf$BF,
    fit_bf$BF,
    tolerance = 1e-8
  )),
  isTRUE(all.equal(fit_bf$BF, direct_bf$BF, tolerance = 1e-8)),
  isTRUE(validate_linear_mixture_fash(
    fit_bf,
    expected_grid = test_grid,
    expected_pred_step = test_pred_step,
    expected_penalty = test_penalty
  ))
)

invalid_likelihood <- expected_likelihood
colnames(invalid_likelihood)[2] <- "mismatched"
invalid_likelihood_error <- tryCatch(
  {
    fit_linear_mixture_fash_from_log_likelihood(
      invalid_likelihood,
      grid = test_grid,
      pred_step = test_pred_step,
      penalty = test_penalty
    )
    NULL
  },
  error = identity
)
stopifnot(
  inherits(invalid_likelihood_error, "error"),
  grepl("L_matrix must be", conditionMessage(invalid_likelihood_error))
)

compact_fit <- compact_linear_mixture_fash(fit)
compact_fit_bf <- compact_linear_mixture_fash(fit_bf)
stopifnot(
  inherits(compact_fit, "compact_linear_mixture_fash"),
  inherits(compact_fit_bf, "compact_linear_mixture_fash"),
  !any(c("L_matrix", "eb_result") %in% names(compact_fit)),
  !any(c("L_matrix", "eb_result") %in% names(compact_fit_bf)),
  identical(compact_fit$unit_ids, rownames(fit$L_matrix)),
  identical(compact_fit$prior_weights, fit$prior_weights),
  identical(compact_fit$posterior_weights, fit$posterior_weights),
  identical(compact_fit$psd_grid, fit$psd_grid),
  identical(compact_fit$lfdr, fit$lfdr),
  identical(compact_fit$settings, fit$settings),
  identical(compact_fit$bf_adjusted, fit$bf_adjusted),
  is.null(compact_fit$BF),
  identical(compact_fit_bf$unit_ids, rownames(fit_bf$L_matrix)),
  identical(compact_fit_bf$prior_weights, fit_bf$prior_weights),
  identical(compact_fit_bf$posterior_weights, fit_bf$posterior_weights),
  identical(compact_fit_bf$psd_grid, fit_bf$psd_grid),
  identical(compact_fit_bf$lfdr, fit_bf$lfdr),
  identical(compact_fit_bf$settings, fit_bf$settings),
  identical(compact_fit_bf$bf_adjusted, fit_bf$bf_adjusted),
  identical(compact_fit_bf$BF, fit_bf$BF),
  isTRUE(validate_compact_linear_mixture_fash(
    compact_fit,
    expected_grid = test_grid,
    expected_pred_step = test_pred_step,
    expected_penalty = test_penalty
  )),
  isTRUE(validate_compact_linear_mixture_fash(
    compact_fit_bf,
    expected_grid = test_grid,
    expected_pred_step = test_pred_step,
    expected_penalty = test_penalty
  ))
)

fit_without_null <- fit
fit_without_null$prior_weights <- fit_without_null$prior_weights[-1, , drop = FALSE]
fit_without_null$prior_weights$prior_weight <-
  fit_without_null$prior_weights$prior_weight /
  sum(fit_without_null$prior_weights$prior_weight)
fit_without_null$posterior_weights <-
  fit_without_null$posterior_weights[, -1, drop = FALSE]
fit_without_null$posterior_weights <-
  fit_without_null$posterior_weights /
  rowSums(fit_without_null$posterior_weights)
fit_without_null$lfdr <- stats::setNames(
  numeric(nrow(fit_without_null$posterior_weights)),
  rownames(fit_without_null$posterior_weights)
)
stopifnot(
  abs(sum(fit_without_null$prior_weights$prior_weight) - 1) < 1e-10,
  all(abs(rowSums(fit_without_null$posterior_weights) - 1) < 1e-10),
  all(fit_without_null$lfdr == 0)
)
missing_null_error <- tryCatch(
  {
    validate_linear_mixture_fash(fit_without_null)
    NULL
  },
  error = identity
)
stopifnot(
  inherits(missing_null_error, "error"),
  grepl(
    "prior_weights contains invalid|posterior weights are invalid",
    conditionMessage(missing_null_error)
  )
)

raw_prior_table <- extract_linear_mixture_prior_table(
  fit,
  seed = 20260811L,
  fit_label = "Raw"
)
bf_prior_summary <- summarize_linear_mixture_prior_fit(
  fit_bf,
  seed = 20260811L,
  fit_label = "BF-corrected"
)
expanded_bf_prior <- expand_grid_prior_weights(
  fit_bf$prior_weights,
  test_grid
)
conditional_bf_prior <- expanded_bf_prior[-1] /
  sum(expanded_bf_prior[-1])
stopifnot(
  nrow(raw_prior_table) == length(test_grid),
  identical(raw_prior_table$predstep_sd, test_grid),
  sum(raw_prior_table$is_null) == 1L,
  abs(sum(raw_prior_table$prior_weight) - 1) < 1e-10,
  bf_prior_summary$estimated_pi0 == expanded_bf_prior[1],
  bf_prior_summary$active_nonnull_components ==
    sum(expanded_bf_prior[-1] > 0),
  isTRUE(all.equal(
    bf_prior_summary$alternative_rms_predstep_sd,
    sqrt(sum(conditional_bf_prior * test_grid[-1]^2)),
    tolerance = 1e-10
  ))
)

cat("Linear-mixture FASH tests passed.\n")
