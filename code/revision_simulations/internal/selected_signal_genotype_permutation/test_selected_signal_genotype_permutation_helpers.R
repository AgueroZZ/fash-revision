#!/usr/bin/env Rscript

# Deterministic tests for the selected-signal matched-null permutation helpers.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation",
  "selected_signal_genotype_permutation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "covariance_estimation",
  "donor_null_permutation_helpers.R"
))

pair_keys <- c("g1_v1", "g1_v2", "g2_v1", "g3_v1", "g3_v2")
bf_lfdr <- c(0.01, 0.04, 0.02, 0.5, 0.6)
selection <- select_discovered_gene_top_pairs(
  pair_keys,
  bf_lfdr,
  alpha = 0.05
)
if (!identical(selection$pair_key, c("g1_v1", "g2_v1")) ||
    attr(selection, "n_pair_level_discoveries") != 3L ||
    anyDuplicated(selection$gene_id)) {
  stop("The discovered-gene top-pair selector returned the wrong units.")
}

random_pair_keys <- c(
  "g1_v1", "g1_v2", "g1_v3", "g2_v1", "g2_v2", "g3_v1", "g3_v2"
)
random_bf_lfdr <- c(0.01, 0.8, 0.9, 0.02, 0.7, 0.8, 0.9)
set.seed(20260817)
random_state_before <- .Random.seed
random_selection <- select_random_tested_pair_per_discovered_gene(
  random_pair_keys,
  random_bf_lfdr,
  alpha = 0.05,
  selection_seed = 8L
)
random_selection_repeat <- select_random_tested_pair_per_discovered_gene(
  random_pair_keys,
  random_bf_lfdr,
  alpha = 0.05,
  selection_seed = 8L
)
if (!identical(random_selection, random_selection_repeat) ||
    !identical(.Random.seed, random_state_before) ||
    !identical(random_selection$gene_id, c("g1", "g2")) ||
    !identical(random_selection$candidate_variant_count, c(3L, 2L)) ||
    any(random_selection$source_pair_level_discovery) ||
    anyDuplicated(random_selection$gene_id) ||
    any(random_selection$pair_key %in% c("g3_v1", "g3_v2"))) {
  stop("The random tested-pair selector failed its deterministic invariants.")
}

set.seed(20260818)
all_gene_state_before <- .Random.seed
all_gene_selection <- select_random_tested_pair_per_gene(
  random_pair_keys,
  random_bf_lfdr,
  alpha = 0.05,
  selection_seed = 8L
)
all_gene_selection_repeat <- select_random_tested_pair_per_gene(
  random_pair_keys,
  random_bf_lfdr,
  alpha = 0.05,
  selection_seed = 8L
)
if (!identical(all_gene_selection, all_gene_selection_repeat) ||
    !identical(.Random.seed, all_gene_state_before) ||
    !identical(all_gene_selection$gene_id, c("g1", "g2", "g3")) ||
    !identical(
      all_gene_selection$candidate_variant_count,
      c(3L, 2L, 2L)
    ) ||
    !identical(
      attr(all_gene_selection, "candidate_pool"),
      "all_tested_variants_within_all_genes"
    ) ||
    all(all_gene_selection$source_pair_level_discovery) ||
    anyDuplicated(all_gene_selection$gene_id)) {
  stop("The all-gene random tested-pair selector failed its invariants.")
}

genotype <- matrix(
  seq_len(24),
  nrow = 6L,
  dimnames = list(paste0("d", 1:6), paste0("v", 1:4))
)
permutation <- make_shared_genotype_permutation(genotype, seed = 20260811)
expected <- genotype[permutation$donor_map$source_donor, , drop = FALSE]
rownames(expected) <- rownames(genotype)
if (!isTRUE(all.equal(permutation$genotype, expected, tolerance = 0)) ||
    anyDuplicated(permutation$donor_map$source_donor) ||
    !isTRUE(all.equal(
      crossprod(genotype),
      crossprod(permutation$genotype),
      tolerance = 0
    ))) {
  stop("The shared genotype permutation did not preserve matrix structure.")
}

time_donors <- paste0("td", 1:6)
time_observation_matrix <- matrix(
  c(
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
    TRUE, TRUE, TRUE, TRUE, FALSE, FALSE,
    TRUE, FALSE, TRUE, TRUE, TRUE, FALSE
  ),
  nrow = length(time_donors),
  dimnames = list(time_donors, paste0("time_", 1:3))
)
set.seed(20260818)
time_state_before <- .Random.seed
time_donor_map <- make_independent_time_donor_permutations(
  time_observation_matrix,
  time_grid = c(0, 1, 2),
  seed = 20260811L
)
time_donor_map_repeat <- make_independent_time_donor_permutations(
  time_observation_matrix,
  time_grid = c(0, 1, 2),
  seed = 20260811L
)
time_map_valid <- vapply(seq_len(ncol(time_observation_matrix)), function(j) {
  observed <- time_donors[time_observation_matrix[, j]]
  current <- time_donor_map[time_donor_map$time_index == j, , drop = FALSE]
  setequal(current$target_donor, observed) &&
    setequal(current$source_donor, observed) &&
    !anyDuplicated(current$target_donor) &&
    !anyDuplicated(current$source_donor) &&
    any(!current$fixed_point)
}, logical(1))
if (!identical(time_donor_map, time_donor_map_repeat) ||
    !identical(.Random.seed, time_state_before) ||
    !all(time_map_valid) ||
    any(time_donor_map$fixed_point !=
          (time_donor_map$target_donor == time_donor_map$source_donor))) {
  stop("The independent time-specific donor permutations failed invariants.")
}

time_genotype <- matrix(
  seq_len(length(time_donors) * 4L),
  nrow = length(time_donors),
  dimnames = list(time_donors, paste0("variant_", 1:4))
)
for (time_index in seq_len(ncol(time_observation_matrix))) {
  observed <- time_donors[time_observation_matrix[, time_index]]
  current_map <- time_donor_map[
    time_donor_map$time_index == time_index,
    ,
    drop = FALSE
  ]
  permuted_time_genotype <- apply_donor_map_to_genotype(
    genotype = time_genotype,
    donor_map = current_map,
    target_donors = observed
  )
  expected_time_genotype <- time_genotype[
    current_map$source_donor[match(observed, current_map$target_donor)],
    ,
    drop = FALSE
  ]
  rownames(expected_time_genotype) <- observed
  if (!identical(permuted_time_genotype, expected_time_genotype) ||
      !identical(dimnames(permuted_time_genotype), list(
        observed,
        colnames(time_genotype)
      )) ||
      !isTRUE(all.equal(
        unname(crossprod(time_genotype[observed, , drop = FALSE])),
        unname(crossprod(permuted_time_genotype)),
        tolerance = 0
      ))) {
    stop("The time-specific genotype permutation failed invariants.")
  }
}

set.seed(20260816)
residual_donors <- paste0("rd", 1:8)
observation_patterns <- setNames(
  c("111", "111", "111", "111", "110", "110", "101", "101"),
  residual_donors
)
source_donor <- make_shared_donor_block_permutation(
  residual_donors,
  observation_patterns
)
if (anyDuplicated(source_donor) ||
    any(observation_patterns[source_donor] != observation_patterns)) {
  stop("The residual donor-block map crossed an observation-pattern stratum.")
}

set.seed(20260819)
unit_map_state_before <- .Random.seed
unit_keys <- paste0("unit_", 1:5)
unit_donor_map <- make_unit_specific_donor_block_permutations(
  donor_ids = residual_donors,
  observation_patterns = observation_patterns,
  unit_keys = unit_keys,
  seed = 20260811L
)
unit_donor_map_repeat <- make_unit_specific_donor_block_permutations(
  donor_ids = residual_donors,
  observation_patterns = observation_patterns,
  unit_keys = unit_keys,
  seed = 20260811L
)
unit_map_signatures <- vapply(unit_keys, function(unit_key) {
  paste(
    unit_donor_map$source_donor[unit_donor_map$unit_key == unit_key],
    collapse = ":"
  )
}, character(1))
valid_unit_maps <- vapply(seq_along(unit_keys), function(unit_index) {
  current <- unit_donor_map[
    unit_donor_map$unit_index == unit_index,
    ,
    drop = FALSE
  ]
  identical(current$target_donor, residual_donors) &&
    setequal(current$source_donor, residual_donors) &&
    !anyDuplicated(current$source_donor) &&
    any(!current$fixed_point) &&
    all(
      observation_patterns[current$source_donor] ==
        observation_patterns[current$target_donor]
    )
}, logical(1))
if (!identical(unit_donor_map, unit_donor_map_repeat) ||
    !identical(.Random.seed, unit_map_state_before) ||
    !all(valid_unit_maps) || length(unique(unit_map_signatures)) < 2L) {
  stop("The unit-specific donor-block maps failed deterministic invariants.")
}

n_residual_units <- 3L
residual_covariates <- matrix(
  stats::rnorm(length(residual_donors)),
  ncol = 1L,
  dimnames = list(residual_donors, "PC1")
)
residual_genotype <- matrix(
  stats::rnorm(length(residual_donors) * n_residual_units),
  nrow = length(residual_donors),
  dimnames = list(residual_donors, paste0("unit_", 1:n_residual_units))
)
residual_expression <- cbind(1, residual_covariates) %*%
  matrix(stats::rnorm(2L * n_residual_units), nrow = 2L) +
  matrix(
    stats::rnorm(length(residual_donors) * n_residual_units),
    nrow = length(residual_donors)
  )
projection <- make_covariate_projection(residual_covariates)
expression_residual <- projection$residualizer %*% residual_expression
source_rows <- match(source_donor[residual_donors], residual_donors)
permuted_expression_residual <- projection$residualizer %*%
  expression_residual[source_rows, , drop = FALSE]
residual_fit <- fit_residualized_genotype_regressions(
  expression_residual = permuted_expression_residual,
  genotype = residual_genotype,
  residualizer = projection$residualizer,
  covariate_rank = projection$rank
)
for (unit_index in seq_len(n_residual_units)) {
  reference_data <- data.frame(
    y = permuted_expression_residual[, unit_index],
    genotype = residual_genotype[, unit_index],
    PC1 = residual_covariates[, 1L]
  )
  reference_fit <- summary(
    stats::lm(y ~ genotype + PC1, data = reference_data)
  )$coefficients
  if (!isTRUE(all.equal(
    unname(residual_fit$beta[unit_index]),
    unname(reference_fit["genotype", "Estimate"]),
    tolerance = 1e-10
  )) || !isTRUE(all.equal(
    unname(residual_fit$standard_error[unit_index]),
    unname(reference_fit["genotype", "Std. Error"]),
    tolerance = 1e-10
  ))) {
    stop("The residual-block vectorized regression disagrees with lm().")
  }
}

signal_beta <- c(1.75, -1.25, 0.8)
signal_nuisance_coefficients <- matrix(
  c(0.5, -0.2, 0.9, 0.4, -0.6, 0.3),
  nrow = 2L
)
signal_noise <- matrix(
  c(
    -0.30, 0.10, 0.25, -0.05, 0.20, -0.15, 0.05, -0.10,
    0.15, -0.25, 0.05, 0.30, -0.20, 0.10, -0.05, 0.12,
    -0.18, 0.22, -0.08, 0.16, 0.04, -0.12, 0.20, -0.06
  ),
  nrow = length(residual_donors),
  ncol = n_residual_units
)
signal_expression <- cbind(1, residual_covariates) %*%
  signal_nuisance_coefficients +
  sweep(residual_genotype, 2L, signal_beta, `*`) +
  signal_noise
signal_stripped <- make_signal_stripped_residual_block_null(
  expression = signal_expression,
  genotype = residual_genotype,
  covariates = residual_covariates,
  source_rows = source_rows
)
for (unit_index in seq_len(n_residual_units)) {
  observed_data <- data.frame(
    y = signal_expression[, unit_index],
    genotype = residual_genotype[, unit_index],
    PC1 = residual_covariates[, 1L]
  )
  observed_reference <- stats::lm(
    y ~ genotype + PC1,
    data = observed_data
  )
  observed_coefficients <- stats::coef(observed_reference)
  nuisance_fitted <- observed_coefficients[["(Intercept)"]] +
    observed_coefficients[["PC1"]] * residual_covariates[, 1L]
  null_expression <- nuisance_fitted +
    stats::residuals(observed_reference)[source_rows]
  null_reference <- summary(stats::lm(
    null_expression ~ residual_genotype[, unit_index] +
      residual_covariates[, 1L]
  ))$coefficients
  if (!isTRUE(all.equal(
    unname(signal_stripped$null_fit$beta[unit_index]),
    unname(null_reference[2L, "Estimate"]),
    tolerance = 1e-10
  )) || !isTRUE(all.equal(
    unname(signal_stripped$null_fit$standard_error[unit_index]),
    unname(null_reference[2L, "Std. Error"]),
    tolerance = 1e-10
  ))) {
    stop("The signal-stripped residual null disagrees with direct lm().")
  }
}

unit_specific_source_rows <- vapply(seq_len(n_residual_units), function(j) {
  unit_map <- unit_donor_map[
    unit_donor_map$unit_index == j,
    ,
    drop = FALSE
  ]
  match(unit_map$source_donor, residual_donors)
}, integer(length(residual_donors)))
unit_specific_signal_stripped <- make_signal_stripped_residual_block_null(
  expression = signal_expression,
  genotype = residual_genotype,
  covariates = residual_covariates,
  source_rows = unit_specific_source_rows
)
for (unit_index in seq_len(n_residual_units)) {
  observed_data <- data.frame(
    y = signal_expression[, unit_index],
    genotype = residual_genotype[, unit_index],
    PC1 = residual_covariates[, 1L]
  )
  observed_reference <- stats::lm(y ~ genotype + PC1, data = observed_data)
  observed_coefficients <- stats::coef(observed_reference)
  nuisance_fitted <- observed_coefficients[["(Intercept)"]] +
    observed_coefficients[["PC1"]] * residual_covariates[, 1L]
  null_expression <- nuisance_fitted + stats::residuals(observed_reference)[
    unit_specific_source_rows[, unit_index]
  ]
  null_reference <- summary(stats::lm(
    null_expression ~ residual_genotype[, unit_index] +
      residual_covariates[, 1L]
  ))$coefficients
  if (!isTRUE(all.equal(
    unname(unit_specific_signal_stripped$null_fit$beta[unit_index]),
    unname(null_reference[2L, "Estimate"]),
    tolerance = 1e-10
  )) || !isTRUE(all.equal(
    unname(unit_specific_signal_stripped$null_fit$standard_error[unit_index]),
    unname(null_reference[2L, "Std. Error"]),
    tolerance = 1e-10
  ))) {
    stop("The unit-specific signal-stripped null disagrees with direct lm().")
  }
}

augmented_signal_expression <- signal_expression + sweep(
  residual_genotype,
  2L,
  c(4.0, -3.0, 2.5),
  `*`
)
augmented_signal_stripped <- make_signal_stripped_residual_block_null(
  expression = augmented_signal_expression,
  genotype = residual_genotype,
  covariates = residual_covariates,
  source_rows = source_rows
)
if (!isTRUE(all.equal(
  signal_stripped$null_expression,
  augmented_signal_stripped$null_expression,
  tolerance = 1e-10
)) || !isTRUE(all.equal(
  signal_stripped$null_fit$beta,
  augmented_signal_stripped$null_fit$beta,
  tolerance = 1e-10
)) || !isTRUE(all.equal(
  signal_stripped$null_fit$standard_error,
  augmented_signal_stripped$null_fit$standard_error,
  tolerance = 1e-10
)) || signal_stripped$maximum_residual_genotype_correlation > 1e-10) {
  stop("The signal-stripped null retained an added genotype component.")
}

raw_se <- matrix(c(0.2, 0.3, 0.4, 0.5), nrow = 2L)
beta_hat <- matrix(c(0.8, -0.6, 0.2, -0.1), nrow = 2L)
residual_df <- c(12, 9)
adjusted_se <- convert_raw_to_original_t_adjusted_se(
  beta_hat,
  raw_se,
  residual_df
)
reference_adjusted_se <- raw_se
for (time_index in seq_len(ncol(beta_hat))) {
  t_value <- beta_hat[, time_index] / raw_se[, time_index]
  p_value <- 2 * stats::pt(
    abs(t_value),
    df = residual_df[time_index],
    lower.tail = FALSE
  )
  z_value <- stats::qnorm(1 - p_value / 2)
  reference_adjusted_se[, time_index] <-
    abs(beta_hat[, time_index]) / abs(z_value)
}
if (!isTRUE(all.equal(adjusted_se, reference_adjusted_se, tolerance = 0))) {
  stop("The original t-adjusted SE calculation was not reproduced exactly.")
}

lfdr <- c(0.01, 0.02, 0.03, 0.4)
group <- c("target", "target", "permuted_null", "permuted_null")
calibration <- summarize_matched_null_calibration(
  lfdr = lfdr,
  group = group,
  pi0_merged = 0.75,
  fit_stage = "Test",
  alpha = 0.05
)
if (calibration$target_calls != 2L ||
    calibration$permuted_null_calls != 1L ||
    calibration$pi0_target_unbounded != 0.5 ||
    !calibration$pi0_target_valid ||
    abs(calibration$known_null_discovery_fraction - 1 / 3) > 1e-12 ||
    abs(calibration$scaled_fdr_merged_from_estimated_pi0 - 0.5) > 1e-12 ||
    abs(calibration$post_selection_fdr_target_from_pi0 - 0.25) > 1e-12) {
  stop("The matched-null calibration formulas returned the wrong values.")
}

no_calls <- summarize_matched_null_calibration(
  lfdr = rep(0.9, 4),
  group = group,
  pi0_merged = 0.5,
  fit_stage = "No calls",
  alpha = 0.05
)
if (no_calls$total_calls != 0L ||
    !is.na(no_calls$known_null_discovery_fraction) ||
    !is.na(no_calls$scaled_fdr_merged_from_estimated_pi0) ||
    !is.na(no_calls$post_selection_fdr_target_from_pi0)) {
  stop("Zero-call calibration outputs were not handled safely.")
}

invalid_target_pi0 <- summarize_matched_null_calibration(
  lfdr = lfdr,
  group = group,
  pi0_merged = 0.4,
  fit_stage = "Invalid target pi0",
  alpha = 0.05
)
if (!invalid_target_pi0$pi0_merged_below_design_lower_bound ||
    invalid_target_pi0$pi0_target_valid ||
    !is.na(invalid_target_pi0$post_selection_fdr_target_from_pi0)) {
  stop("An invalid implied target pi0 was not flagged correctly.")
}

cat("Selected-signal matched-null permutation helper tests passed.\n")
