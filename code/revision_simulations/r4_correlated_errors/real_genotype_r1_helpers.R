# Shared construction for R4 conditions aligned with the formal R1 design.

validate_r4_real_genotype_cache <- function(genotype_cache,
                                            seed_list,
                                            J = 6362L,
                                            n_donors = 19L) {
  if (!is.list(genotype_cache) ||
      !all(c("configuration", "sample_ids", "samples") %in%
        names(genotype_cache)) ||
      !identical(genotype_cache$configuration$n_genes, as.integer(J)) ||
      !identical(
        genotype_cache$configuration$n_donors,
        as.integer(n_donors)
      ) ||
      !all(seed_list %in% genotype_cache$configuration$seed_list) ||
      !all(as.character(seed_list) %in% names(genotype_cache$samples))) {
    stop("The R4 real-genotype cache has an unexpected schema.")
  }
  for (seed in seed_list) {
    validate_real_genotype_sample(
      genotype_cache$samples[[as.character(seed)]],
      expected_genes = as.integer(J),
      expected_donors = as.integer(n_donors),
      maf_min = genotype_cache$configuration$maf_min
    )
  }
  genotype_cache
}

prepare_r4_real_genotype_r1_seed <- function(genotype_cache,
                                              seed,
                                              J = 6362L,
                                              n_donors = 19L,
                                              n_covariates = 5L,
                                              time_grid = 0:15,
                                              class_probs = c(
                                                dynamic_bspline = 0.20,
                                                constant = 0.40,
                                                zero = 0.40
                                              ),
                                              expression_noise_sd = 1,
                                              dynamic_main_effect_sd = 1) {
  seed <- as.integer(seed)
  J <- as.integer(J)
  n_donors <- as.integer(n_donors)
  n_covariates <- as.integer(n_covariates)
  if (length(seed) != 1L || is.na(seed) ||
      length(time_grid) < 2L || any(!is.finite(time_grid)) ||
      !is.finite(expression_noise_sd) || expression_noise_sd <= 0 ||
      !is.finite(dynamic_main_effect_sd) || dynamic_main_effect_sd <= 0) {
    stop("Invalid formal R1 seed-construction arguments.")
  }
  validate_r4_real_genotype_cache(
    genotype_cache,
    seed_list = seed,
    J = J,
    n_donors = n_donors
  )
  genotype_sample <- validate_real_genotype_sample(
    genotype_cache$samples[[as.character(seed)]],
    expected_genes = J,
    expected_donors = n_donors,
    maf_min = genotype_cache$configuration$maf_min
  )
  component_seeds <- revision_component_seeds(seed)
  covariates <- simulate_covariate_matrix(
    n_donors = n_donors,
    n_covariates = n_covariates,
    seed = component_seeds[["covariates"]]
  )
  scenario <- paste0(
    "r1_real_genotype_one_per_gene_",
    "random_bspline_main_effect_dynamic_eqtl"
  )
  effect_sim <- simulate_variant_effect_curves(
    n_variants = J,
    time_grid = time_grid,
    class_probs = class_probs,
    scenario = scenario,
    dynamic_amplitude = 2,
    bspline_df = 6,
    bspline_coefficient_sd = 1,
    constant_sd = 1,
    dynamic_main_effect_sd = dynamic_main_effect_sd,
    exact_class_counts = TRUE,
    seed = component_seeds[["functional_truth"]]
  )
  effect_sim <- reassign_effect_simulation_by_maf(
    effect_sim = effect_sim,
    maf = genotype_sample$variant_info$observed_maf,
    class_probs = class_probs,
    seed = component_seeds[["classes"]],
    n_strata = 10L
  )
  rownames(effect_sim$beta_matrix) <- genotype_sample$selection$pair_key
  effect_sim$unit_info$variant_id <- genotype_sample$selection$pair_key
  independent_expression <- simulate_eqtl_expression_from_genotypes(
    G = genotype_sample$G,
    beta_matrix = effect_sim$beta_matrix,
    time_grid = time_grid,
    covariates = covariates,
    expression_noise_sd = expression_noise_sd,
    covariate_effect_sd = 0.5,
    intercept_sd = 0,
    seed = component_seeds[["expression"]]
  )
  expected_counts <- exact_proportional_counts(J, class_probs)
  observed_counts <- table(factor(
    effect_sim$unit_info$effect_class,
    levels = names(class_probs)
  ))
  if (!identical(as.integer(observed_counts), as.integer(expected_counts))) {
    stop("The formal R1 truth-class counts are invalid.")
  }
  list(
    seed = seed,
    scenario = scenario,
    component_seeds = component_seeds,
    genotype_sample = genotype_sample,
    covariates = covariates,
    effect_sim = effect_sim,
    independent_expression = independent_expression,
    settings = list(
      J = J,
      n_donors = n_donors,
      n_covariates = n_covariates,
      time_grid = time_grid,
      class_probs = class_probs,
      expression_noise_sd = expression_noise_sd,
      dynamic_main_effect_sd = dynamic_main_effect_sd,
      covariate_effect_sd = 0.5,
      intercept_sd = 0,
      dynamic_amplitude = 2,
      bspline_df = 6,
      bspline_coefficient_sd = 1,
      common_sd_grid = default_revision_grid(),
      pred_step = 1,
      penalty = 10L,
      linear_prior_mode = "mixture_grid",
      genotype_selection_rule = genotype_cache$configuration$selection_rule
    )
  )
}

make_r4_correlated_expression <- function(seed_inputs,
                                          expression_error_correlation = NULL) {
  if (!is.list(seed_inputs) ||
      !all(c(
        "component_seeds", "genotype_sample", "covariates", "effect_sim",
        "settings"
      ) %in% names(seed_inputs))) {
    stop("seed_inputs is not a formal R4 real-genotype seed object.")
  }
  settings <- seed_inputs$settings
  correlation <- if (is.null(expression_error_correlation)) {
    NULL
  } else {
    validate_time_correlation(
      expression_error_correlation,
      n_time = length(settings$time_grid)
    )
  }
  simulate_eqtl_expression_from_genotypes(
    G = seed_inputs$genotype_sample$G,
    beta_matrix = seed_inputs$effect_sim$beta_matrix,
    time_grid = settings$time_grid,
    covariates = seed_inputs$covariates,
    expression_noise_sd = settings$expression_noise_sd,
    covariate_effect_sd = settings$covariate_effect_sd,
    intercept_sd = settings$intercept_sd,
    seed = seed_inputs$component_seeds[["expression"]],
    expression_error_correlation = correlation
  )
}

run_r4_real_genotype_r1_condition <- function(
    seed_inputs,
    expression_error_correlation = NULL,
    num_cores = 1L,
    num_basis = 20L,
    nominal_alpha = 0.05,
    output_dir = tempdir(),
    verbose = FALSE) {
  if (!is.list(seed_inputs) || is.null(seed_inputs$settings)) {
    stop("seed_inputs is not a formal R4 real-genotype seed object.")
  }
  settings <- seed_inputs$settings
  expression_sim <- make_r4_correlated_expression(
    seed_inputs,
    expression_error_correlation = expression_error_correlation
  )
  out <- run_genotype_level_dynamic_eqtl_simulation(
    G = seed_inputs$genotype_sample$G,
    time_grid = settings$time_grid,
    covariates = seed_inputs$covariates,
    class_probs = settings$class_probs,
    expression_noise_sd = settings$expression_noise_sd,
    covariate_effect_sd = settings$covariate_effect_sd,
    intercept_sd = settings$intercept_sd,
    dynamic_amplitude = settings$dynamic_amplitude,
    bspline_df = settings$bspline_df,
    bspline_coefficient_sd = settings$bspline_coefficient_sd,
    dynamic_main_effect_sd = settings$dynamic_main_effect_sd,
    alpha = nominal_alpha,
    seed = seed_inputs$seed,
    num_cores = as.integer(num_cores),
    num_basis = as.integer(num_basis),
    grid = settings$common_sd_grid,
    penalty = settings$penalty,
    pred_step = settings$pred_step,
    linear_prior_mode = settings$linear_prior_mode,
    scenario = seed_inputs$scenario,
    output_dir = output_dir,
    save_outputs = FALSE,
    verbose = verbose,
    effect_sim = seed_inputs$effect_sim,
    expression_sim = expression_sim
  )
  out$settings$genotype_source <- "paper-derived YRI DS dosage"
  out$settings$genotype_selection_rule <-
    settings$genotype_selection_rule
  out$settings$genotype_digest <-
    seed_inputs$genotype_sample$genotype_digest
  out$settings$real_genotype_r1_aligned <- TRUE
  out
}

validate_r4_condition_pair <- function(seed_inputs,
                                       reference_expression,
                                       candidate_expression,
                                       correlation,
                                       tolerance = 1e-10) {
  required_fields <- c("expression", "intercepts", "covariate_effects")
  if (!all(required_fields %in% names(reference_expression)) ||
      !all(required_fields %in% names(candidate_expression))) {
    stop("An R4 expression object is incomplete.")
  }
  if (!identical(
        reference_expression$intercepts,
        candidate_expression$intercepts
      ) ||
      !identical(
        reference_expression$covariate_effects,
        candidate_expression$covariate_effects
      )) {
    stop("Paired R4 conditions differ in nuisance effects.")
  }
  extract_errors <- function(expression_sim) {
    genotype <- seed_inputs$genotype_sample$G
    beta <- seed_inputs$effect_sim$beta_matrix
    covariates <- seed_inputs$covariates
    n_time <- ncol(beta)
    errors <- array(NA_real_, dim = dim(expression_sim$expression))
    for (time_index in seq_len(n_time)) {
      genetic_mean <- sweep(
        genotype,
        2L,
        beta[, time_index],
        `*`
      )
      covariate_mean <- covariates %*%
        expression_sim$covariate_effects[, , time_index]
      systematic_mean <- sweep(
        genetic_mean + covariate_mean,
        2L,
        expression_sim$intercepts[, time_index],
        `+`
      )
      errors[, , time_index] <-
        expression_sim$expression[, , time_index] - systematic_mean
    }
    errors
  }
  reference_errors <- extract_errors(reference_expression)
  candidate_errors <- extract_errors(candidate_expression)
  n_series <- dim(reference_errors)[1L] * dim(reference_errors)[2L]
  n_time <- dim(reference_errors)[3L]
  expected <- matrix(reference_errors, nrow = n_series, ncol = n_time) %*%
    chol(validate_time_correlation(correlation, n_time = n_time))
  observed <- matrix(candidate_errors, nrow = n_series, ncol = n_time)
  maximum_difference <- max(abs(expected - observed))
  if (!is.finite(maximum_difference) || maximum_difference > tolerance) {
    stop(
      "The paired R4 Gaussian innovations are invalid; maximum difference = ",
      format(maximum_difference, scientific = TRUE),
      "."
    )
  }
  maximum_difference
}
