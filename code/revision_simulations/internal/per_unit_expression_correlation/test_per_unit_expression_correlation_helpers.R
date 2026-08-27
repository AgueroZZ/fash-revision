#!/usr/bin/env Rscript
source(file.path(
  "code", "revision_simulations", "internal", "per_unit_expression_correlation",
  "per_unit_expression_correlation_helpers.R"
))

stopifnot_message <- function(condition, label) {
  if (!isTRUE(condition)) stop("FAILED: ", label)
  message("ok  ", label)
}

set.seed(20260825)

# 1. Residualizer is a projection and kills the covariate span.
covariates <- cbind(1, matrix(rnorm(40), nrow = 10))
residualizer <- make_residualizer(covariates)
stopifnot_message(
  max(abs(residualizer %*% residualizer - residualizer)) < 1e-10,
  "residualizer is idempotent"
)
stopifnot_message(
  max(abs(residualizer %*% covariates)) < 1e-10,
  "residualizer annihilates the covariate span"
)

# 2. Genotype weights reproduce the OLS genotype coefficient exactly.
genotype <- round(runif(10, 0, 2), 2)
weights <- make_genotype_weights(residualizer, genotype)
outcome <- rnorm(10)
coefficient <- unname(coef(lm(outcome ~ covariates[, -1] + genotype))["genotype"])
stopifnot_message(
  abs(sum(weights * outcome) - coefficient) < 1e-10,
  "OLS weights reproduce the genotype coefficient"
)

# 3. Matched-donor correlation reproduces stats::cor on complete data.
donors <- paste0("d", 1:12)
truth <- matrix(rnorm(12 * 3), nrow = 12, dimnames = list(donors, NULL))
residuals <- lapply(1:3, function(k) setNames(truth[, k], donors))
stopifnot_message(
  max(abs(matched_donor_correlation(residuals) - cor(truth))) < 1e-12,
  "matched-donor correlation matches stats::cor when donors are complete"
)

# 4. Missing donors are dropped pairwise, not globally.
ragged <- residuals
ragged[[2]] <- ragged[[2]][-1]
observed <- matched_donor_correlation(ragged)
expected13 <- cor(truth[, 1], truth[, 3])
expected12 <- cor(truth[-1, 1], truth[-1, 2])
stopifnot_message(
  abs(observed[1, 3] - expected13) < 1e-12 &&
    abs(observed[1, 2] - expected12) < 1e-12,
  "pairwise-complete handling keeps the full pair where both times are observed"
)

# 5. Identical weight vectors give a design factor of one.
constant <- lapply(1:3, function(k) setNames(rep(0.5, 12), donors))
stopifnot_message(
  max(abs(weight_design_factor(constant) - 1)) < 1e-12,
  "identical weight vectors give design factor one"
)

# 6. Design factor is bounded by one and shrinks when donors are missing.
ragged_weights <- constant
ragged_weights[[2]] <- ragged_weights[[2]][-(1:4)]
ragged_factor <- weight_design_factor(ragged_weights)
stopifnot_message(
  max(ragged_factor) <= 1 + 1e-12 && ragged_factor[1, 2] < 0.95,
  "design factor is at most one and shrinks under missing donors"
)

# 7. Lag profile on a known AR(1).
rho <- 0.6
ar1 <- rho^abs(outer(1:16, 1:16, "-"))
stopifnot_message(
  max(abs(lag_profile(ar1) - rho^(1:15))) < 1e-12,
  "lag profile recovers AR(1) lags"
)
stopifnot_message(
  abs(mean_off_diagonal(diag(4)) - 0) < 1e-12,
  "mean off-diagonal of the identity is zero"
)

message("All per-unit expression correlation helper tests passed.")
