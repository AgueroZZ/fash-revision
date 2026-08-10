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
  "r4_correlated_errors",
  "real_data_correlation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "covariance_estimation",
  "zero_intercept_correlation_helpers.R"
))

pair_keys <- c(
  "geneA_variant1",
  "geneA_variant2",
  "geneB_variant3",
  "geneC_variant4",
  "geneC_variant5"
)
first_selection <- select_random_variant_per_gene(pair_keys, seed = 11L)
second_selection <- select_random_variant_per_gene(pair_keys, seed = 11L)
stopifnot(
  identical(first_selection, second_selection),
  nrow(first_selection) == 3L,
  !anyDuplicated(first_selection$gene_id),
  first_selection$pair_key[first_selection$gene_id == "geneB"] ==
    "geneB_variant3"
)

gallery_metadata <- data.frame(
  pair_key = sprintf("gene%03d_variant%03d", 1:80, 1:80),
  lfdr = seq(0.999, 0.900, length.out = 80),
  stringsAsFactors = FALSE
)
first_gallery_sample <- stratified_lfdr_gallery_sample(
  gallery_metadata,
  n_per_quartile = 5L,
  seed = 105L
)
second_gallery_sample <- stratified_lfdr_gallery_sample(
  gallery_metadata,
  n_per_quartile = 5L,
  seed = 105L
)
quartile_bounds <- cbind(
  lower = c(1L, 21L, 41L, 61L),
  upper = c(20L, 40L, 60L, 80L)
)
stopifnot(
  identical(first_gallery_sample, second_gallery_sample),
  nrow(first_gallery_sample) == 20L,
  !anyDuplicated(first_gallery_sample$pair_key),
  identical(
    as.integer(table(first_gallery_sample$lfdr_quartile)),
    rep(5L, 4L)
  ),
  all(vapply(1:4, function(quartile) {
    ranks <- first_gallery_sample$retained_lfdr_rank[
      first_gallery_sample$lfdr_quartile == quartile
    ]
    all(ranks >= quartile_bounds[quartile, "lower"]) &&
      all(ranks <= quartile_bounds[quartile, "upper"]) &&
      !is.unsorted(ranks)
  }, logical(1)))
)

threshold_matrix <- rbind(
  c(0, 1.99, rep(0, 14)),
  c(0, 2.00, rep(0, 14)),
  c(0, 2.01, rep(0, 14))
)
threshold_result <- filter_zero_intercept_z(threshold_matrix, threshold = 2)
stopifnot(identical(threshold_result$keep, c(TRUE, FALSE, FALSE)))

mean_z_matrix <- rbind(
  rep(0.49, 16L),
  rep(0.50, 16L),
  c(rep(1.90, 8L), rep(-1.90, 8L)),
  seq(-0.75, 0.75, length.out = 16L)
)
mean_z_result <- filter_zero_intercept_z_with_mean(
  mean_z_matrix,
  max_threshold = 2,
  mean_z_threshold = 2
)
stopifnot(
  isTRUE(all.equal(
    mean_z_result$mean_z_score,
    sqrt(ncol(mean_z_matrix)) * rowMeans(mean_z_matrix)
  )),
  identical(mean_z_result$keep_maximum, rep(TRUE, 4L)),
  identical(mean_z_result$keep_mean_z, c(TRUE, FALSE, TRUE, TRUE)),
  identical(mean_z_result$keep, c(TRUE, FALSE, TRUE, TRUE)),
  all(!mean_z_result$keep | mean_z_result$keep_maximum),
  mean_z_result$n_selected_by_maximum == 4L,
  mean_z_result$n_selected == 3L
)

set.seed(101)
n_independent <- 200000L
z_independent <- matrix(
  stats::rnorm(n_independent * 16L),
  nrow = n_independent,
  ncol = 16L
)
independent_filter <- filter_zero_intercept_z(z_independent, threshold = 2)
independent_estimate <- estimate_zero_intercept_correlation(
  z_independent[independent_filter$keep, , drop = FALSE]
)
independent_upper <- independent_estimate$sample_correlation[
  upper.tri(independent_estimate$sample_correlation)
]
stopifnot(
  independent_filter$n_selected > 80000L,
  max(abs(independent_upper)) < 0.03,
  min(independent_estimate$second_moment_diagonal) > 0.70,
  max(independent_estimate$second_moment_diagonal) < 0.85
)

set.seed(102)
intercept <- c(rep(0, 50000L), rep(4, 50000L))
z_intercept_mixture <- matrix(
  stats::rnorm(length(intercept) * 16L),
  nrow = length(intercept),
  ncol = 16L
) + intercept
intercept_filter <- filter_zero_intercept_z(z_intercept_mixture, threshold = 2)
stopifnot(
  mean(intercept_filter$keep[intercept == 4]) < 0.001,
  max(abs(colMeans(
    z_intercept_mixture[intercept_filter$keep, , drop = FALSE]
  ))) < 0.03
)

positive_target <- make_lag1_only_correlation(16L, rho = 0.25)
positive_simulation <- simulate_truncated_correlation(
  n_candidates = 200000L,
  threshold = 2,
  correlation = positive_target,
  seed = 103L
)
positive_lag1 <- mean(positive_simulation$sample_correlation[cbind(
  1:15,
  2:16
)])
stopifnot(positive_lag1 > 0.15, positive_lag1 < 0.25)

bootstrap_input <- z_independent[independent_filter$selected_indices[1:500], ]
first_bootstrap <- bootstrap_zero_intercept_correlations(
  bootstrap_input,
  n_bootstrap = 30L,
  seed = 104L
)
second_bootstrap <- bootstrap_zero_intercept_correlations(
  bootstrap_input,
  n_bootstrap = 30L,
  seed = 104L
)
stopifnot(
  identical(first_bootstrap, second_bootstrap),
  nrow(first_bootstrap) == 30L,
  all(is.finite(unlist(first_bootstrap[c("mean", "lower", "upper")]))),
  all(first_bootstrap$lower <= first_bootstrap$upper)
)

projected <- project_to_positive_definite_correlation(
  positive_simulation$sample_correlation
)
stopifnot(
  projected$diagnostics$projected_minimum_eigenvalue > 0,
  max(abs(diag(projected$projected) - 1)) < 1e-10
)

cat("Zero-intercept correlation helper tests passed.\n")
