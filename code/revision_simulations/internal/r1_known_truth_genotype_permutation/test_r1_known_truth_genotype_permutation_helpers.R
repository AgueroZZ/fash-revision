#!/usr/bin/env Rscript

# Deterministic tests for the R1 known-truth permutation helpers.

helper_directory <- file.path(
  "code", "revision_simulations", "internal",
  "r1_known_truth_genotype_permutation"
)
source(file.path(
  "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation",
  "selected_signal_genotype_permutation_helpers.R"
))
source(file.path(
  helper_directory,
  "r1_known_truth_genotype_permutation_helpers.R"
))

lfdr <- c(0.01, 0.02, 0.20, 0.90)
true_null <- c(FALSE, TRUE, FALSE, TRUE)
curve <- known_truth_alpha_curve(
  lfdr = lfdr,
  true_null = true_null,
  alpha_grid = c(0, 0.01, 0.02, 0.10),
  arm = "toy",
  fit_stage = "BF"
)
alpha_002 <- curve[curve$alpha == 0.02, ]
if (nrow(curve) != 4L || alpha_002$n_discoveries != 2L ||
    alpha_002$false_discoveries != 1L ||
    abs(alpha_002$realized_fdp - 0.5) > 1e-12 ||
    curve$n_discoveries[curve$alpha == 0] != 0L) {
  stop("The known-truth alpha curve returned incorrect toy results.")
}

bf_summary <- summarize_null_bf(
  bf = c(100, 0.5, 25, 2),
  true_null = true_null,
  arm = "toy"
)
if (bf_summary$n_null != 2L ||
    abs(bf_summary$mean_bf - 1.25) > 1e-12 ||
    abs(bf_summary$median_bf - 1.25) > 1e-12 ||
    abs(bf_summary$proportion_bf_greater_than_one - 0.5) > 1e-12) {
  stop("The known-null BF summary returned incorrect toy results.")
}

original_genotype <- matrix(
  c(
    0, 0, 1, 1, 2, 2,
    0, 1, 0, 2, 1, 2
  ),
  nrow = 6,
  ncol = 2
)
colnames(original_genotype) <- c("u1", "u2")
permuted_genotype <- original_genotype[c(2, 1, 4, 3, 6, 5), , drop = FALSE]
covariates <- matrix(c(-2, -1, 0, 0, 1, 2), ncol = 1)
alignment <- residualized_genotype_alignment(
  original_genotype,
  permuted_genotype,
  covariates
)
if (nrow(alignment) != 2L || !identical(alignment$unit_key, c("u1", "u2")) ||
    any(!is.finite(as.matrix(alignment[, -1L]))) ||
    any(abs(alignment$residualized_correlation) > 1 + 1e-12)) {
  stop("The residualized genotype alignment is invalid.")
}

toy_donors <- paste0("d", seq_len(6L))
toy_units <- c("u1", "u2")
toy_times <- paste0("t", seq_len(3L))
toy_genotype <- matrix(
  c(
    0, 1, 0, 2, 1, 2,
    2, 0, 1, 0, 2, 1
  ),
  nrow = 6L,
  ncol = 2L,
  dimnames = list(toy_donors, toy_units)
)
toy_covariates <- matrix(
  c(-1.4, -0.8, -0.1, 0.3, 0.9, 1.7),
  ncol = 1L,
  dimnames = list(toy_donors, "PC1")
)
toy_expression <- array(
  NA_real_,
  dim = c(6L, 2L, 3L),
  dimnames = list(toy_donors, toy_units, toy_times)
)
for (unit_index in seq_len(2L)) {
  for (time_index in seq_len(3L)) {
    toy_expression[, unit_index, time_index] <-
      0.4 * unit_index - 0.2 * time_index +
      (0.7 + 0.1 * time_index) * toy_genotype[, unit_index] +
      (0.3 - 0.05 * unit_index) * toy_covariates[, 1L] +
      c(-0.4, 0.2, 0.1, -0.1, 0.35, -0.15) * (unit_index + time_index)
  }
}
toy_donor_map <- data.frame(
  target_donor = toy_donors,
  source_donor = toy_donors[c(2, 1, 4, 3, 6, 5)],
  fixed_point = FALSE,
  stringsAsFactors = FALSE
)
toy_null <- make_signal_stripped_residual_null(
  genotype = toy_genotype,
  expression = toy_expression,
  covariates = toy_covariates,
  donor_map = toy_donor_map,
  leverage_adjustment = "HC2"
)
for (unit_index in seq_len(2L)) {
  design <- cbind(
    intercept = 1,
    genotype = toy_genotype[, unit_index],
    toy_covariates
  )
  design_qr <- qr(design)
  coefficients <- qr.coef(
    design_qr,
    toy_expression[, unit_index, , drop = FALSE][, 1L, ]
  )
  fitted_full <- design %*% coefficients
  residual <- toy_expression[, unit_index, , drop = FALSE][, 1L, ] -
    fitted_full
  leverage <- rowSums(qr.Q(design_qr)^2)
  adjusted_residual <- residual / sqrt(1 - leverage)
  nuisance_fitted <- cbind(intercept = 1, toy_covariates) %*%
    coefficients[c("intercept", "PC1"), , drop = FALSE]
  expected_null <- nuisance_fitted + adjusted_residual[
    match(toy_donor_map$source_donor, toy_donors), , drop = FALSE
  ]
  if (max(abs(
        toy_null$full_model_residual[, unit_index, ] - residual
      )) > 1e-12 ||
      max(abs(
        toy_null$adjusted_residual[, unit_index, ] - adjusted_residual
      )) > 1e-12 ||
      max(abs(
        toy_null$nuisance_fitted[, unit_index, ] - nuisance_fitted
      )) > 1e-12 ||
      max(abs(
        toy_null$null_expression[, unit_index, ] - expected_null
      )) > 1e-12 ||
      max(abs(toy_null$leverage[, unit_index] - leverage)) > 1e-12) {
    stop("The signal-stripped residual null did not reproduce direct OLS.")
  }
}
if (toy_null$diagnostics$maximum_full_residual_design_cross_product > 1e-10 ||
    toy_null$diagnostics$maximum_nuisance_partial_genotype_coefficient > 1e-10 ||
    !identical(toy_null$donor_map, toy_donor_map)) {
  stop("The signal-stripped residual-null invariants failed.")
}

cat("R1 known-truth genotype-permutation helper tests passed.\n")
