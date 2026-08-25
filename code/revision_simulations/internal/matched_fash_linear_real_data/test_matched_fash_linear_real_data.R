#!/usr/bin/env Rscript

workflowr_root <- if (file.exists("_workflowr.yml")) {
  "."
} else if (file.exists("coderepo-local/_workflowr.yml")) {
  "coderepo-local"
} else {
  stop("Run this test from the workflowr root or its parent.")
}

source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "fash_linear_real_data_ablation", "fash_linear_real_data_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "matched_fash_linear_real_data", "matched_fash_linear_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "matched_fash_linear_real_data", "iwp_overlap_proportion_helpers.R"
))

package_provenance <- extract_package_provenance_matched()
stopifnot(
  identical(package_provenance$version, "0.1.43"),
  identical(
    package_provenance$remote_sha,
    "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
  )
)

set.seed(20260821)
toy_datasets <- lapply(seq_len(40L), function(index) {
  time <- 0:15
  standard_error <- 0.08 + 0.015 * ((time + index) %% 4L)
  endpoint <- if (index %% 3L == 0L) 1.5 else 0
  data.frame(
    time = time,
    beta = 0.1 * index / 40 + endpoint * time / max(time) +
      stats::rnorm(length(time), sd = standard_error),
    SE = standard_error,
    stringsAsFactors = FALSE
  )
})
names(toy_datasets) <- paste0(
  "ENSG", sprintf("%011d", seq_along(toy_datasets)),
  "_rs", 100000 + seq_along(toy_datasets)
)

statistics <- compute_linear_sufficient_statistics_block(toy_datasets)
direct_datasets <- lapply(toy_datasets, function(dataset) {
  list(x = dataset$time, y = dataset$beta, sd = dataset$SE)
})
grid <- c(0, 0.02, 0.05, 0.1, 0.25, 0.5)
direct_likelihood <- compute_linear_mixture_log_likelihood(
  datasets = direct_datasets,
  grid = grid,
  pred_step = 1
)
accelerated_likelihood <- compute_linear_mixture_log_likelihood_from_stats(
  statistics = statistics,
  grid = grid,
  pred_step = 1,
  statistic_time_span = 15,
  n_time = 16L
)
stopifnot(
  identical(dimnames(direct_likelihood), dimnames(accelerated_likelihood)),
  max(abs(direct_likelihood - accelerated_likelihood)) < 1e-8
)

raw_fit <- fit_linear_mixture_fash_from_stats(
  statistics = statistics,
  grid = grid,
  pred_step = 1,
  penalty = 10L,
  statistic_time_span = 15,
  n_time = 16L
)
raw_fit$settings$betaprec <- 0
bf_fit <- BF_update_linear_mixture_fash(raw_fit)
raw_weights <- expand_grid_prior_weights(raw_fit$prior_weights, grid)
bf_weights <- expand_grid_prior_weights(bf_fit$prior_weights, grid)
raw_conditional_alternative <- raw_weights[-1L] / sum(raw_weights[-1L])
bf_conditional_alternative <- bf_weights[-1L] / sum(bf_weights[-1L])
stopifnot(
  max(abs(raw_conditional_alternative - bf_conditional_alternative)) < 1e-10,
  all(is.finite(bf_fit$BF)),
  all(bf_fit$lfdr >= 0 & bf_fit$lfdr <= 1),
  identical(bf_fit$settings$betaprec, 0)
)

toy_sets <- list(
  `Strober quadratic` = c("a", "b", "c", "d"),
  `Strober linear` = c("b", "c", "e"),
  `FASH-IWP1 BF` = c("c", "d", "f"),
  `FASH-linear BF` = c("c", "e", "f", "g")
)
regions <- exclusive_set_regions_matched(toy_sets, "Toy units")
pairwise <- pairwise_overlap_matched(toy_sets, "Toy units")
toy_reference <- list(
  `Gene-variant pairs` = c("g1_v1", "g1_v2", "g2_v3"),
  Genes = c("g1", "g2"),
  Variants = c("v1", "v2", "v3")
)
toy_comparison <- list(
  `Gene-variant pairs` = c("g1_v1", "g3_v4"),
  Genes = c("g1", "g3"),
  Variants = c("v1", "v4")
)
directional <- directional_overlap_summary_matched(
  toy_reference,
  toy_comparison
)
stopifnot(
  nrow(regions) == 15L,
  sum(regions$count) == length(unique(unlist(toy_sets))),
  nrow(pairwise) == 6L,
  regions$count[regions$region == paste(names(toy_sets), collapse = " & ")] == 1L,
  identical(directional$intersection_count, c(1L, 1L, 1L)),
  identical(directional$comparison_covered_by_reference, rep(0.5, 3L))
)

cat("Matched FASH-linear helper tests passed.\n")
