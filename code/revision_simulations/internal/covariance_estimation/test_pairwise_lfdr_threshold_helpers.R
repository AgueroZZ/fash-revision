#!/usr/bin/env Rscript

# Deterministic tests for the pairwise lfdr-threshold exploration helpers.

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
  "pairwise_lfdr_threshold_helpers.R"
))

n_genes <- 60L
gene_id <- sprintf("gene%03d", seq_len(n_genes))
pair_keys <- as.vector(rbind(
  paste0(gene_id, "_variant_b"),
  paste0(gene_id, "_variant_a")
))
representative_lfdr <- seq(0.999, 0.951, length.out = n_genes)
lfdr <- as.vector(rbind(
  representative_lfdr,
  representative_lfdr - 0.02
))

# Force one exact tie and require lexical pair-key tie breaking.
lfdr[1:2] <- 0.999

thresholds <- c(0.97, 0.96, 0.95)
selections <- lapply(thresholds, function(threshold) {
  select_gene_representatives_above_lfdr(
    pair_keys,
    lfdr,
    threshold
  )
})

if (!identical(selections[[1]]$selected$pair_key[1], "gene001_variant_a") ||
    anyDuplicated(selections[[3]]$selected$gene_id) ||
    any(diff(selections[[3]]$selected$lfdr) > 0)) {
  stop("The representative selection is not deterministic and one-per-gene.")
}
if (!all(selections[[1]]$selected$lfdr > thresholds[1]) ||
    !all(selections[[2]]$selected$lfdr > thresholds[2]) ||
    !all(selections[[3]]$selected$lfdr > thresholds[3])) {
  stop("The strict lfdr threshold rule was not applied.")
}
if (!all(selections[[1]]$selected$pair_key %in%
         selections[[2]]$selected$pair_key) ||
    !all(selections[[2]]$selected$pair_key %in%
         selections[[3]]$selected$pair_key)) {
  stop("The threshold-selected sets are not nested.")
}

set.seed(38102)
n_units <- 800L
n_time <- 6L
true_correlation <- outer(seq_len(n_time), seq_len(n_time), function(a, b) {
  0.25^abs(a - b)
})
se <- matrix(
  stats::runif(n_units * n_time, min = 0.7, max = 1.3),
  nrow = n_units,
  ncol = n_time
)
constant_effect <- stats::rnorm(n_units, sd = 0.8)
errors <- matrix(
  stats::rnorm(n_units * n_time),
  nrow = n_units,
  ncol = n_time
) %*% chol(true_correlation)
beta_hat <- constant_effect + se * errors

first_bootstrap <- bootstrap_pairwise_lag_variogram(
  beta_hat,
  se,
  n_bootstrap = 30L,
  seed = 20260821L
)
second_bootstrap <- bootstrap_pairwise_lag_variogram(
  beta_hat,
  se,
  n_bootstrap = 30L,
  seed = 20260821L
)
if (!identical(first_bootstrap, second_bootstrap) ||
    nrow(first_bootstrap$summary) != n_time - 1L ||
    any(!is.finite(first_bootstrap$draws)) ||
    any(first_bootstrap$summary$correlation_ci_lower >
        first_bootstrap$summary$correlation_ci_upper) ||
    any(first_bootstrap$summary$semivariogram_ci_lower >
        first_bootstrap$summary$semivariogram_ci_upper)) {
  stop("The pairwise bootstrap summary is invalid or non-deterministic.")
}

cat("Pairwise lfdr-threshold helper tests passed.\n")
