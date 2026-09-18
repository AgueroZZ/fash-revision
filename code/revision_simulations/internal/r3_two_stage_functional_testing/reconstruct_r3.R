# Reconstruct inputs with the frozen formal R3 functions, without changing truth.

r3ts_validate_configuration <- function(config) {
  contract <- r3ts_contract()
  r3ts_assert(identical(config$output_id, contract$formal_id), "Incorrect base R3 cache.")
  checks <- list(
    J = 6362L, n_donors = 19L, n_covariates = 5L,
    time_grid = 0:15, evaluation_grid = seq(0, 15, by = 0.1),
    middle_window = c(3, 12), middle_boundary = "open", switch_threshold = 0.25,
    seed_list = contract$seeds, truth_mechanisms = contract$mechanisms,
    class_probs = c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40),
    expected_class_counts = c(dynamic_bspline = 1272L, constant = 2545L, zero = 2545L),
    temporal_category_probs = c(early = 0.29, middle = 0.42, late = 0.29),
    location_truth_margin = 0.10, location_truth_min_range_fraction = 0.10,
    switch_truth_margin = 0.10, non_switch_min_abs = 0.10,
    non_switch_min_range_fraction = 0.10, dynamic_main_effect_sd = 1,
    expression_noise_sd = 1, covariate_effect_sd = 0.5, intercept_sd = 0,
    num_basis = 20L
  )
  for (field in names(checks)) {
    r3ts_assert(isTRUE(all.equal(config[[field]], checks[[field]])),
                 paste("Unexpected formal setting:", field))
  }
  r3ts_assert(identical(unname(config$expected_truth_group_counts),
                         c(185L, 184L, 267L, 267L, 185L, 184L)), "Incorrect truth-cell counts.")
  r3ts_assert(identical(config$raised_cosine$center_ranges,
                         list(early = c(1.5, 2.5), middle = c(4.5, 10.5), late = c(12.5, 13.5))),
               "Incorrect raised-cosine center ranges.")
  r3ts_assert(identical(config$random_bspline, list(amplitude = 2, df = 6, coefficient_sd = 1)) &&
                 isTRUE(all.equal(config$raised_cosine[c("width_half", "spike_counts",
                    "relative_amplitude_range", "target_centered_rms")],
                    list(width_half = 1.5, spike_counts = 1:3,
                         relative_amplitude_range = c(0.35, 0.75), target_centered_rms = 0.9))),
               "Formal curve scales or peak settings changed.")
  invisible(TRUE)
}

r3ts_reconstruct <- function(config, genotype_sample, seed, mechanism) {
  cs <- revision_component_seeds(seed)
  n <- ncol(genotype_sample$G)
  # Preserve the formal driver's historical scenario metadata as well as its inputs.
  scenario <- paste0(
    if (mechanism == "random_bspline") "r3a_" else "r3b_",
    "real_genotype_one_per_gene_matched_functional_", mechanism,
    "_open_middle_3_12_center_aligned_equal_cells_relative_location_clearance_",
    "full_universe_paired_posterior_main_effect"
  )
  covariates <- simulate_covariate_matrix(
    n_donors = config$n_donors, n_covariates = config$n_covariates,
    seed = cs[["covariates"]]
  )
  effects <- simulate_matched_functional_effect_set(
    n_variants = n, truth_mechanism = mechanism,
    time_grid = config$time_grid, evaluation_grid = config$evaluation_grid,
    class_probs = config$class_probs, dynamic_main_effect_sd = config$dynamic_main_effect_sd,
    cosine_center_ranges = config$raised_cosine$center_ranges,
    switch_threshold = config$switch_threshold,
    location_truth_margin = config$location_truth_margin,
    location_truth_min_range_fraction = config$location_truth_min_range_fraction,
    switch_truth_margin = config$switch_truth_margin,
    non_switch_min_abs = config$non_switch_min_abs,
    non_switch_min_range_fraction = config$non_switch_min_range_fraction,
    temporal_category_probs = config$temporal_category_probs,
    seed = seed, class_seed = cs[["classes"]], constant_seed = cs[["constant_effects"]],
    shape_seed = cs[["functional_truth"]], scenario = scenario,
    middle_window = config$middle_window, middle_boundary = config$middle_boundary
  )
  effects <- reassign_effect_simulation_by_maf(
    effect_sim = effects, maf = genotype_sample$variant_info$observed_maf,
    class_probs = config$class_probs, seed = cs[["classes"]], n_strata = 10L
  )
  for (field in c("beta_matrix", "beta_evaluation", "true_functionals")) {
    rownames(effects[[field]]) <- genotype_sample$selection$pair_key
  }
  effects$unit_info$variant_id <- genotype_sample$selection$pair_key
  expression <- simulate_eqtl_expression_from_genotypes(
    G = genotype_sample$G, beta_matrix = effects$beta_matrix,
    time_grid = config$time_grid, covariates = covariates,
    expression_noise_sd = config$expression_noise_sd,
    covariate_effect_sd = config$covariate_effect_sd, intercept_sd = config$intercept_sd,
    seed = cs[["expression"]]
  )
  summaries <- estimate_eqtl_summaries_from_genotypes(
    G = genotype_sample$G, expression = expression$expression,
    covariates = covariates, apply_t_se_correction = TRUE
  )
  digests <- c(
    genotype_content_md5 = genotype_content_md5(
      genotype_sample$selection$pair_key, rownames(genotype_sample$G), genotype_sample$G
    ),
    true_beta_md5 = serialized_object_md5(effects$beta_matrix),
    adjusted_se_md5 = serialized_object_md5(summaries$se),
    regression_beta_hat_md5 = serialized_object_md5(summaries$beta_hat)
  )
  list(seed = seed, truth_mechanism = mechanism, component_seeds = cs,
       scenario = scenario, effects = effects, covariates = covariates,
       expression = expression, regression = summaries, input_digests = digests)
}

r3ts_verify_reconstruction <- function(inputs, reference_digests, config, formal_replicate) {
  ref <- reference_digests[reference_digests$seed == inputs$seed &
                            reference_digests$truth_mechanism == inputs$truth_mechanism, , drop = FALSE]
  r3ts_assert(nrow(ref) == 1L, "Missing or duplicated reference input digests.")
  for (field in names(inputs$input_digests)) {
    r3ts_assert(identical(unname(inputs$input_digests[[field]]), ref[[field]][[1L]]),
                 paste("Reconstruction does not match retained R3 inputs:", field))
  }
  effects <- inputs$effects
  dynamic <- effects$unit_info$effect_class == "dynamic_bspline"
  counts <- table(factor(effects$unit_info$effect_class, levels = names(config$expected_class_counts)))
  cells <- table(factor(effects$unit_info$truth_group[dynamic],
                        levels = names(config$expected_truth_group_counts)))
  r3ts_assert(identical(as.integer(counts), unname(config$expected_class_counts)) &&
                 identical(as.integer(cells), unname(config$expected_truth_group_counts)),
               "Reconstructed truth counts differ from formal R3.")
  r3ts_assert(identical(effects$unit_info$variant_id, formal_replicate$selected_pair_keys),
               "Reconstructed pair order differs from formal R3.")
  calls <- formal_replicate$call_diagnostics_alpha005
  for (target in r3ts_contract()$targets) {
    x <- calls[calls$target == target, , drop = FALSE]
    r3ts_assert(isTRUE(all.equal(unname(effects$true_functionals[x$unit_index, target]),
                                  x$true_functional, tolerance = 1e-12)),
                 paste("Retained functional truth mismatch:", target))
  }
  invisible(TRUE)
}

r3ts_fit_all <- function(inputs, config, num_cores = 1L) {
  datasets <- make_fash_datasets_from_eqtl_summary(
    beta_hat = inputs$regression$beta_hat, se = inputs$regression$se,
    true_beta = inputs$effects$beta_matrix, time_grid = config$time_grid,
    unit_info = inputs$effects$unit_info, scenario = inputs$scenario
  )
  # The formal wrapper resets this seed before fitting the complete data list.
  set.seed(inputs$seed)
  fit_fash_for_revision(
    datasets = datasets, orders = 1, grid = default_revision_grid(),
    num_basis = config$num_basis, penalty = 10, pred_step = 1,
    num_cores = num_cores, apply_bf = TRUE, verbose = FALSE
  )
}

r3ts_infer <- function(fit, inputs, config, num_cores = 1L) {
  contract <- r3ts_contract()
  n <- nrow(inputs$effects$true_functionals)
  r3ts_assert(length(fit$fash_data$data_list) == n && length(get_fash_lfdr(fit)) == n,
               "FASH was not fitted to every simulated unit.")
  fdr_table <- get_fash_fdr_table(fit)
  candidates <- r3ts_stage1(fdr_table, n, contract$q_dyn)
  functionals <- make_temporal_functionals(
    config$evaluation_grid, switch_threshold = config$switch_threshold,
    middle_window = config$middle_window, middle_boundary = config$middle_boundary
  )
  lfsr <- compute_functional_lfsr(
    fit = fit, functionals = functionals, indices = candidates,
    smooth_var = config$evaluation_grid, num_cores = num_cores,
    seed = inputs$component_seeds[["functional_posterior"]]
  )
  stage2 <- r3ts_stage2(candidates, lfsr, n, contract$q_func)
  truth <- inputs$effects$true_functionals[, contract$targets, drop = FALSE]
  dynamic <- inputs$effects$unit_info$effect_class == "dynamic_bspline"
  metrics <- r3ts_metrics(truth, dynamic, candidates, stage2,
                          inputs$seed, inputs$truth_mechanism)
  list(candidates = candidates, fdr_table = fdr_table, lfsr_by_target = lfsr,
       stage2 = stage2, metrics = metrics, true_functionals = truth,
       true_dynamic = dynamic)
}
