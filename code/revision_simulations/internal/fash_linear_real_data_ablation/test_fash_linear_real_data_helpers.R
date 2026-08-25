#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

find_workflowr_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "_workflowr.yml"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the workflowr project root.")
    }
    current <- parent
  }
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "shared",
  "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "fash_linear_real_data_ablation",
  "fash_linear_real_data_helpers.R"
))

set.seed(20260810)
toy_datasets <- lapply(seq_len(30L), function(index) {
  x <- 0:15
  scaled_x <- x / 15
  se <- 0.05 + 0.02 * ((x + index) %% 5)
  signal <- 0.15 * index / 30 +
    if (index %% 3L == 0L) 0.8 * scaled_x else 0 +
    if (index %% 5L == 0L) 0.4 * sin(2 * pi * scaled_x) else 0
  data.frame(
    x = x,
    y = signal + stats::rnorm(length(x), sd = se),
    sd = se,
    stringsAsFactors = FALSE
  )
})
names(toy_datasets) <- paste0("gene", seq_along(toy_datasets), "_rs", 1000 + seq_along(toy_datasets))

real_data_format <- lapply(toy_datasets, function(dataset) {
  data.frame(
    time = dataset$x,
    beta = dataset$y,
    SE = dataset$sd,
    stringsAsFactors = FALSE
  )
})
statistics <- compute_linear_sufficient_statistics_block(real_data_format)

for (slope_sd in c(0, 0.05, 0.5, 2)) {
  reference <- vapply(toy_datasets, function(dataset) {
    log_marginal_common_intercept(
      y = dataset$y,
      se = dataset$sd,
      x = dataset$x / 15,
      slope_sd = slope_sd
    )
  }, numeric(1))
  accelerated <- linear_log_marginal_from_stats(
    statistics,
    slope_sd = slope_sd
  )
  stopifnot(max(abs(reference - accelerated)) < 1e-8)
}

sigma_grid <- c(0.05, 0.2, 0.8, 2)
reference_fit <- fit_simplified_fash(
  toy_datasets,
  estimate_sigma = TRUE,
  sigma_beta_grid = sigma_grid,
  scale_time = TRUE
)
accelerated_fit <- fit_profiled_linear_fash_from_stats(
  statistics,
  sigma_beta_grid = sigma_grid
)
stopifnot(
  identical(reference_fit$sigma_beta, accelerated_fit$sigma_beta),
  max(abs(
    reference_fit$prior_weights$prior_weight -
      accelerated_fit$prior_weights$prior_weight
  )) < 1e-8,
  max(abs(reference_fit$lfdr - accelerated_fit$lfdr)) < 1e-8,
  max(abs(
    reference_fit$sigma_profile$loglik -
      accelerated_fit$sigma_profile$loglik
  )) < 1e-7
)

reference_bf <- BF_update_simplified_fash(reference_fit)
accelerated_bf <- bf_update_profiled_linear_fash(accelerated_fit)
stopifnot(
  max(abs(
    reference_bf$prior_weights$prior_weight -
      accelerated_bf$prior_weights$prior_weight
  )) < 1e-12,
  max(abs(reference_bf$lfdr - accelerated_bf$lfdr)) < 1e-8
)

mixture_grid <- c(0, 0.01, 0.025, 0.05, 0.1, 0.25)
mixture_penalty <- 10L
mixture_pred_step <- 1
reference_mixture_likelihood <- compute_linear_mixture_log_likelihood(
  datasets = toy_datasets,
  grid = mixture_grid,
  pred_step = mixture_pred_step
)
accelerated_mixture_likelihood <-
  compute_linear_mixture_log_likelihood_from_stats(
    statistics = statistics,
    grid = mixture_grid,
    pred_step = mixture_pred_step,
    statistic_time_span = 15,
    n_time = 16L
  )
stopifnot(
  identical(
    dimnames(accelerated_mixture_likelihood),
    dimnames(reference_mixture_likelihood)
  ),
  max(abs(
    accelerated_mixture_likelihood - reference_mixture_likelihood
  )) < 1e-8,
  max(abs(
    accelerated_mixture_likelihood[, "0.1"] -
      linear_log_marginal_from_stats(
        statistics,
        slope_sd = 0.1 * 15,
        n_time = 16L
      )
  )) < 1e-10
)

reference_mixture_fit <- fit_linear_mixture_fash(
  datasets = toy_datasets,
  grid = mixture_grid,
  pred_step = mixture_pred_step,
  penalty = mixture_penalty
)
accelerated_mixture_fit <- fit_linear_mixture_fash_from_stats(
  statistics = statistics,
  grid = mixture_grid,
  pred_step = mixture_pred_step,
  penalty = mixture_penalty,
  statistic_time_span = 15,
  n_time = 16L
)
stopifnot(
  max(abs(
    expand_grid_prior_weights(
      reference_mixture_fit$prior_weights,
      mixture_grid
    ) -
      expand_grid_prior_weights(
        accelerated_mixture_fit$prior_weights,
        mixture_grid
      )
  )) < 1e-8,
  max(abs(
    reference_mixture_fit$posterior_weights -
      accelerated_mixture_fit$posterior_weights
  )) < 1e-8,
  max(abs(
    reference_mixture_fit$lfdr - accelerated_mixture_fit$lfdr
  )) < 1e-8
)

reference_mixture_bf <- BF_update_linear_mixture_fash(
  reference_mixture_fit
)
accelerated_mixture_bf <- BF_update_linear_mixture_fash(
  accelerated_mixture_fit
)
stopifnot(
  max(abs(
    expand_grid_prior_weights(
      reference_mixture_bf$prior_weights,
      mixture_grid
    ) -
      expand_grid_prior_weights(
        accelerated_mixture_bf$prior_weights,
        mixture_grid
      )
  )) < 1e-8,
  max(abs(
    reference_mixture_bf$posterior_weights -
      accelerated_mixture_bf$posterior_weights
  )) < 1e-8,
  max(abs(reference_mixture_bf$lfdr - accelerated_mixture_bf$lfdr)) < 1e-8,
  isTRUE(all.equal(
    reference_mixture_bf$BF,
    accelerated_mixture_bf$BF,
    tolerance = 1e-8
  ))
)

mixture_unit_index <- which.min(reference_mixture_bf$lfdr)
reference_mixture_plot <- extract_linear_mixture_posterior_plot_data(
  dataset = toy_datasets[[mixture_unit_index]],
  standard_error = toy_datasets[[mixture_unit_index]]$sd,
  fit = reference_mixture_bf,
  unit_index = mixture_unit_index,
  sample_size = 5000L,
  seed = 20260811L
)
accelerated_mixture_plot <- extract_linear_mixture_posterior_plot_data(
  dataset = toy_datasets[[mixture_unit_index]],
  standard_error = toy_datasets[[mixture_unit_index]]$sd,
  fit = accelerated_mixture_bf,
  unit_index = mixture_unit_index,
  sample_size = 5000L,
  seed = 20260811L
)
compact_mixture_bf <- compact_linear_mixture_fash(
  accelerated_mixture_bf
)
compact_mixture_plot <- extract_linear_mixture_posterior_plot_data(
  dataset = toy_datasets[[mixture_unit_index]],
  standard_error = toy_datasets[[mixture_unit_index]]$sd,
  fit = compact_mixture_bf,
  unit_index = mixture_unit_index,
  sample_size = 5000L,
  seed = 20260811L
)
stopifnot(
  nrow(reference_mixture_plot) == 151L,
  all(is.finite(as.matrix(reference_mixture_plot))),
  all(reference_mixture_plot$lower <= reference_mixture_plot$upper),
  max(abs(
    as.matrix(reference_mixture_plot) -
      as.matrix(accelerated_mixture_plot)
  )) < 1e-8,
  identical(compact_mixture_plot, accelerated_mixture_plot),
  identical(compact_mixture_bf$unit_ids, names(toy_datasets)),
  identical(compact_mixture_bf$prior_weights,
            accelerated_mixture_bf$prior_weights),
  identical(compact_mixture_bf$posterior_weights,
            accelerated_mixture_bf$posterior_weights),
  identical(compact_mixture_bf$lfdr, accelerated_mixture_bf$lfdr),
  identical(compact_mixture_bf$BF, accelerated_mixture_bf$BF)
)

calls <- select_cumulative_lfdr_calls_linear(c(0.001, 0.01, 0.1, 0.9))
stopifnot(identical(calls, 1:3))

pair_table <- parse_linear_pair_keys(names(toy_datasets))
top_indices <- select_current_top_pair_per_gene(
  pair_table,
  rev(seq_along(toy_datasets)) / length(toy_datasets)
)
stopifnot(length(top_indices) == nrow(pair_table))

posterior_plot_data <- extract_profiled_linear_posterior_plot_data(
  dataset = data.frame(
    x = toy_datasets[[1]]$x,
    y = toy_datasets[[1]]$y
  ),
  standard_error = toy_datasets[[1]]$sd,
  slope_sd = accelerated_bf$sigma_beta,
  lfdr = accelerated_bf$lfdr[1],
  sample_size = 1000L,
  seed = 1L
)
stopifnot(
  nrow(posterior_plot_data) == 151L,
  all(is.finite(as.matrix(posterior_plot_data))),
  all(posterior_plot_data$lower <= posterior_plot_data$upper)
)

cat("FASH-linear real-data helper tests passed.\n")
