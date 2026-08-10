#!/usr/bin/env Rscript

# Focused tests for the observational null-set covariance pilot.

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
  "code", "revision_simulations", "r4_correlated_errors",
  "real_data_correlation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "covariance_estimation",
  "observational_null_set_helpers.R"
))

discovery <- discover_by_cumulative_lfdr(
  c(0.001, 0.01, 0.04, 0.30, 0.90),
  alpha = 0.05
)
stopifnot(
  identical(discovery$cutoff_rank, 3L),
  identical(which(discovery$discovered), 1:3)
)

pair_keys <- c(
  "gene_a_variant_1",
  "gene_a_variant_shared",
  "gene_b_variant_2",
  "gene_b_variant_shared",
  "gene_c_variant_3",
  "gene_d_variant_4"
)
# The production keys contain one underscore. Use an equivalent compact catalog
# here so that the shared-variant exclusion can be tested directly.
catalog <- data.frame(
  fash_index = seq_along(pair_keys),
  pair_key = pair_keys,
  gene_id = c("gene_a", "gene_a", "gene_b", "gene_b", "gene_c", "gene_d"),
  variant_id = c(
    "variant_1", "variant_shared", "variant_2", "variant_shared",
    "variant_3", "variant_4"
  ),
  variant_gene_count = c(1L, 2L, 1L, 2L, 1L, 1L),
  pair_random_score = c(0.2, 0.1, 0.2, 0.1, 0.3, 0.4),
  gene_random_score = c(0.1, 0.1, 0.2, 0.2, 0.3, 0.4),
  stringsAsFactors = FALSE
)
selection <- select_random_unique_variant_per_gene(
  catalog,
  eligible_genes = c("gene_a", "gene_b", "gene_c", "gene_d"),
  n_select = 3L,
  set_id = "test"
)
stopifnot(
  nrow(selection$selected) == 3L,
  !any(selection$selected$variant_id == "variant_shared"),
  !anyDuplicated(selection$selected$gene_id),
  !anyDuplicated(selection$selected$variant_id)
)

randomized_a <- make_randomized_pair_catalog(
  c("geneA_rs1", "geneA_rs2", "geneB_rs3"),
  seed = 31L
)
randomized_b <- make_randomized_pair_catalog(
  c("geneA_rs1", "geneA_rs2", "geneB_rs3"),
  seed = 31L
)
stopifnot(identical(randomized_a, randomized_b))

gallery_metadata <- data.frame(
  pair_key = paste0("gene", seq_len(40L), "_rs", seq_len(40L)),
  bf_adjusted_lfdr = seq(0.51, 0.99, length.out = 40L),
  stringsAsFactors = FALSE
)
gallery_a <- sample_lfdr_quartile_gallery(
  gallery_metadata,
  n_per_quartile = 5L,
  seed = 41L
)
gallery_b <- sample_lfdr_quartile_gallery(
  gallery_metadata,
  n_per_quartile = 5L,
  seed = 41L
)
stopifnot(
  identical(gallery_a, gallery_b),
  nrow(gallery_a) == 20L,
  !anyDuplicated(gallery_a$pair_key),
  identical(as.integer(table(gallery_a$lfdr_quartile)), rep(5L, 4L))
)

set.seed(20260806)
n_units <- 5000L
n_time <- 4L
rho <- 0.25
truth <- outer(seq_len(n_time), seq_len(n_time), function(x, y) {
  rho^abs(x - y)
})
standardized_error <- matrix(
  stats::rnorm(n_units * n_time),
  nrow = n_units
) %*% chol(truth)
se <- matrix(
  stats::runif(n_units * n_time, min = 0.06, max = 0.24),
  nrow = n_units
)
colnames(se) <- paste0("time_", 0:(n_time - 1L))
constant_effect <- stats::rnorm(n_units, sd = 0.4)
beta_constant <- constant_effect + se * standardized_error
beta_zero <- se * standardized_error
colnames(beta_constant) <- colnames(beta_zero) <- colnames(se)

pairwise_estimate <- estimate_pairwise_difference_correlation(
  beta_constant,
  se
)
direct_estimate <- estimate_direct_zero_mean_covariance(beta_zero, se)
centered_estimate <- estimate_within_unit_centered_matrices(
  beta_constant,
  se
)
stopifnot(
  max(abs(pairwise_estimate - truth)) < 0.06,
  max(abs(direct_estimate - truth)) < 0.06,
  identical(
    centered_estimate$covariance,
    crossprod(centered_estimate$standardized_residual) / n_units
  ),
  max(abs(diag(centered_estimate$correlation) - 1)) < 1e-12
)

lag_summary <- lag_covariance_variogram(pairwise_estimate, se)
scaled <- scale_covariance_by_median_se(pairwise_estimate, se)
stopifnot(
  nrow(lag_summary) == n_time - 1L,
  all(is.finite(as.matrix(lag_summary[, -1L]))),
  identical(dim(scaled), c(n_time, n_time))
)

bootstrap <- bootstrap_covariance_lags(
  beta_zero[seq_len(100L), , drop = FALSE],
  se[seq_len(100L), , drop = FALSE],
  estimators = c(
    "Pairwise difference",
    "Direct zero mean",
    "Within-unit centered covariance"
  ),
  n_bootstrap = 20L,
  seed = 17L
)
stopifnot(
  nrow(bootstrap$summary) == 3L * (n_time - 1L),
  all(is.finite(bootstrap$summary$observed_semivariogram)),
  all(bootstrap$summary$n_bootstrap == 20L)
)

cat("All observational null-set helper tests passed.\n")
