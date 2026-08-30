# Unit and static contract tests for the R3 ideal-Gaussian measurement experiment.

find_workflowr_root <- function() {
  candidates <- c(".", "coderepo-local", "../..", "../../..", "../../../..")
  marker <- file.path(
    candidates,
    "code", "revision_simulations", "shared", "simulation_functions.R"
  )
  hits <- candidates[file.exists(marker)]
  if (length(hits) == 0L) {
    stop("Could not locate the workflowr project root.", call. = FALSE)
  }
  normalizePath(hits[[1L]], winslash = "/", mustWork = TRUE)
}

expect_error <- function(expression, pattern) {
  message <- tryCatch(
    {
      force(expression)
      NA_character_
    },
    error = function(condition) conditionMessage(condition)
  )
  stopifnot(!is.na(message), grepl(pattern, message, fixed = TRUE))
}

project_root <- find_workflowr_root()
source(file.path(
  project_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
helper_path <- file.path(
  project_root,
  "code", "revision_simulations", "internal",
  "r3_ideal_gaussian_measurement", "ideal_gaussian_measurement.R"
)
source(helper_path)

seed <- 12345L
ideal_seed <- ideal_gaussian_measurement_seed(seed)
stopifnot(
  identical(ideal_seed, ideal_gaussian_measurement_seed(seed)),
  length(ideal_seed) == 1L,
  is.integer(ideal_seed),
  ideal_seed > 0L,
  !ideal_seed %in% unname(revision_component_seeds(seed)),
  ideal_seed != ideal_gaussian_measurement_seed(seed + 1L)
)

true_beta <- matrix(
  seq(-1.1, 1.2, length.out = 24L),
  nrow = 3L,
  ncol = 8L,
  dimnames = list(
    paste0("unit_", 1:3),
    paste0("time_", 1:8)
  )
)
adjusted_se <- matrix(
  seq(0.1, 0.4, length.out = length(true_beta)),
  nrow = nrow(true_beta),
  ncol = ncol(true_beta),
  dimnames = dimnames(true_beta)
)
regression_summary <- list(
  beta_hat = true_beta + 0.1,
  se = adjusted_se,
  se_uncorrected = adjusted_se * 0.9,
  df = matrix(12L, nrow(true_beta), ncol(true_beta)),
  rank = rep(7L, nrow(true_beta)),
  apply_t_se_correction = TRUE
)

ideal_a <- simulate_ideal_gaussian_eqtl_summary(
  true_beta = true_beta,
  adjusted_se = adjusted_se,
  seed = ideal_seed,
  regression_summary = regression_summary
)
ideal_b <- simulate_ideal_gaussian_eqtl_summary(
  true_beta = true_beta,
  adjusted_se = adjusted_se,
  seed = ideal_seed,
  regression_summary = regression_summary
)
ideal_c <- simulate_ideal_gaussian_eqtl_summary(
  true_beta = true_beta,
  adjusted_se = adjusted_se,
  seed = ideal_gaussian_measurement_seed(seed + 1L),
  regression_summary = regression_summary
)

stopifnot(
  identical(ideal_a, ideal_b),
  !identical(ideal_a$beta_hat, ideal_c$beta_hat),
  identical(dim(ideal_a$beta_hat), dim(true_beta)),
  identical(dimnames(ideal_a$beta_hat), dimnames(true_beta)),
  identical(ideal_a$se, adjusted_se),
  isTRUE(all.equal(
    (ideal_a$beta_hat - true_beta) / adjusted_se,
    ideal_a$standard_normal_draw,
    tolerance = 1e-14,
    check.attributes = TRUE
  )),
  identical(ideal_a$apply_t_se_correction, FALSE),
  identical(ideal_a$se_source, "R3 final t-adjusted SE treated as fixed")
)

expect_error(
  simulate_ideal_gaussian_eqtl_summary(
    true_beta = true_beta,
    adjusted_se = adjusted_se[-1, , drop = FALSE],
    seed = ideal_seed
  ),
  "must have identical dimensions"
)
invalid_se <- adjusted_se
invalid_se[1, 1] <- 0
expect_error(
  simulate_ideal_gaussian_eqtl_summary(
    true_beta = true_beta,
    adjusted_se = invalid_se,
    seed = ideal_seed
  ),
  "strictly positive"
)
bad_regression_summary <- regression_summary
bad_regression_summary$se[1, 1] <- adjusted_se[1, 1] + 1
expect_error(
  simulate_ideal_gaussian_eqtl_summary(
    true_beta = true_beta,
    adjusted_se = adjusted_se,
    seed = ideal_seed,
    regression_summary = bad_regression_summary
  ),
  "does not match adjusted_se"
)

moment_summary <- summarize_standardized_error_matrix(
  z = ideal_a$standardized_error,
  seed = seed,
  truth_mechanism = "random_bspline",
  observation_model = "ideal_gaussian"
)
stopifnot(
  nrow(moment_summary) == 1L,
  moment_summary$n == length(true_beta),
  identical(moment_summary$seed, seed),
  identical(moment_summary$truth_mechanism, "random_bspline"),
  identical(moment_summary$observation_model, "ideal_gaussian"),
  all(is.finite(unlist(moment_summary[setdiff(
    names(moment_summary),
    c("truth_mechanism", "observation_model")
  )])))
)

histogram_summary <- summarize_standardized_error_histogram(
  z = ideal_a$standardized_error,
  seed = seed,
  truth_mechanism = "random_bspline",
  observation_model = "ideal_gaussian",
  breaks = seq(-3, 3, by = 0.5)
)
stopifnot(
  sum(histogram_summary$count) == length(true_beta),
  abs(sum(histogram_summary$probability) - 1) < 1e-14,
  nrow(histogram_summary) == length(seq(-3, 3, by = 0.5)) + 1L,
  all(!is.na(histogram_summary$bin_label))
)

lag_summary <- summarize_standardized_error_lag_correlation(
  z = ideal_a$standardized_error,
  seed = seed,
  truth_mechanism = "random_bspline",
  observation_model = "ideal_gaussian"
)
stopifnot(
  nrow(lag_summary) == ncol(true_beta) - 1L,
  identical(lag_summary$lag, seq_len(ncol(true_beta) - 1L)),
  identical(lag_summary$n_pairs, ncol(true_beta) - lag_summary$lag),
  all(is.finite(lag_summary$mean_correlation))
)

digest_a <- serialized_object_md5(adjusted_se)
digest_b <- serialized_object_md5(adjusted_se)
stopifnot(
  identical(digest_a, digest_b),
  grepl("^[[:xdigit:]]{32}$", digest_a)
)

if (tolower(Sys.getenv("FASH_RUN_R3_IDEAL_SMOKE", unset = "false")) %in%
    c("1", "true", "t", "yes", "y")) {
  if (!requireNamespace("fashr", quietly = TRUE)) {
    stop("The fashr package is required for the integration smoke test.")
  }
  source(file.path(
    project_root,
    "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
  ))
  genotype_cache <- readRDS(file.path(
    project_root,
    "output", "revision_simulations", "shared",
    "real_genotype_one_per_gene_J6362_pilot5", "genotype_samples.rds"
  ))
  smoke_seed <- 12345L
  component_seeds <- revision_component_seeds(smoke_seed)
  genotype_sample <- validate_real_genotype_sample(
    genotype_cache$samples[[as.character(smoke_seed)]],
    expected_genes = 6362L,
    expected_donors = 19L,
    maf_min = genotype_cache$configuration$maf_min
  )
  smoke_index <- seq_len(60L)
  smoke_G <- genotype_sample$G[, smoke_index, drop = FALSE]
  smoke_variant_info <- genotype_sample$variant_info[smoke_index, , drop = FALSE]
  smoke_covariates <- simulate_covariate_matrix(
    n_donors = 19L,
    n_covariates = 5L,
    seed = component_seeds[["covariates"]]
  )
  smoke_effect <- simulate_matched_functional_effect_set(
    n_variants = 60L,
    truth_mechanism = "random_bspline",
    time_grid = 0:15,
    evaluation_grid = seq(0, 15, by = 0.1),
    class_probs = c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40),
    dynamic_main_effect_sd = 1,
    cosine_center_ranges = list(
      early = c(1.5, 2.5),
      middle = c(4.5, 10.5),
      late = c(12.5, 13.5)
    ),
    switch_threshold = 0.25,
    location_truth_margin = 0.10,
    location_truth_min_range_fraction = 0.10,
    switch_truth_margin = 0.10,
    non_switch_min_abs = 0.10,
    non_switch_min_range_fraction = 0.10,
    temporal_category_probs = stats::setNames(
      rep(1 / 3, 3), c("early", "middle", "late")
    ),
    seed = smoke_seed,
    class_seed = component_seeds[["classes"]],
    constant_seed = component_seeds[["constant_effects"]],
    shape_seed = component_seeds[["functional_truth"]],
    scenario = "r3_ideal_gaussian_smoke",
    middle_window = c(3, 12),
    middle_boundary = "open"
  )
  smoke_effect <- reassign_effect_simulation_by_maf(
    effect_sim = smoke_effect,
    maf = smoke_variant_info$observed_maf,
    class_probs = c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40),
    seed = component_seeds[["classes"]],
    n_strata = 10L
  )
  for (field in c("beta_matrix", "beta_evaluation", "true_functionals")) {
    rownames(smoke_effect[[field]]) <- colnames(smoke_G)
  }
  smoke_effect$unit_info$variant_id <- colnames(smoke_G)
  smoke_expression <- simulate_eqtl_expression_from_genotypes(
    G = smoke_G,
    beta_matrix = smoke_effect$beta_matrix,
    time_grid = 0:15,
    covariates = smoke_covariates,
    expression_noise_sd = 1,
    covariate_effect_sd = 0.5,
    intercept_sd = 0,
    seed = component_seeds[["expression"]]
  )
  smoke_regression <- estimate_eqtl_summaries_from_genotypes(
    G = smoke_G,
    expression = smoke_expression$expression,
    covariates = smoke_covariates,
    apply_t_se_correction = TRUE
  )
  smoke_ideal <- simulate_ideal_gaussian_eqtl_summary(
    true_beta = smoke_effect$beta_matrix,
    adjusted_se = smoke_regression$se,
    seed = ideal_gaussian_measurement_seed(smoke_seed),
    regression_summary = smoke_regression
  )
  smoke_datasets <- make_fash_datasets_from_eqtl_summary(
    beta_hat = smoke_ideal$beta_hat,
    se = smoke_ideal$se,
    true_beta = smoke_effect$beta_matrix,
    time_grid = 0:15,
    unit_info = smoke_effect$unit_info,
    scenario = "r3_ideal_gaussian_smoke"
  )
  smoke_fits <- fit_fash_for_revision(
    datasets = smoke_datasets,
    orders = 1,
    grid = default_revision_grid(),
    num_basis = 20,
    penalty = 10,
    pred_step = 1,
    num_cores = min(2L, as.integer(Sys.getenv(
      "SLURM_CPUS_PER_TASK", unset = "1"
    ))),
    apply_bf = TRUE,
    verbose = FALSE
  )
  smoke_true_functionals <- smoke_effect$true_functionals[
    , c("early", "middle", "late", "switch"), drop = FALSE
  ]
  smoke_true_dynamic <-
    smoke_effect$unit_info$effect_class == "dynamic_bspline"
  smoke_raw <- evaluate_fash_functional_testing(
    fit = smoke_fits$fash_iwp1_raw,
    true_functionals = smoke_true_functionals,
    evaluation_grid = seq(0, 15, by = 0.1),
    alpha_grid = c(0.05, 0.10),
    method = "FASH-IWP1-Raw",
    scenario = "r3_ideal_gaussian_smoke",
    switch_threshold = 0.25,
    candidate_scope = "full_universe",
    true_dynamic = smoke_true_dynamic,
    num_cores = 2L,
    seed = component_seeds[["functional_posterior"]],
    middle_window = c(3, 12),
    middle_boundary = "open"
  )
  smoke_bf <- evaluate_fash_functional_testing(
    fit = smoke_fits$fash_iwp1_bf,
    true_functionals = smoke_true_functionals,
    evaluation_grid = seq(0, 15, by = 0.1),
    alpha_grid = c(0.05, 0.10),
    method = "FASH-IWP1-BF",
    scenario = "r3_ideal_gaussian_smoke",
    switch_threshold = 0.25,
    candidate_scope = "full_universe",
    true_dynamic = smoke_true_dynamic,
    num_cores = 2L,
    seed = component_seeds[["functional_posterior"]],
    middle_window = c(3, 12),
    middle_boundary = "open"
  )
  smoke_alpha <- rbind(smoke_raw$alpha_curve, smoke_bf$alpha_curve)
  stopifnot(
    nrow(smoke_alpha) == 16L,
    all(smoke_alpha$candidate_scope == "full_universe"),
    all(smoke_alpha$candidate_count == 60L),
    all(smoke_alpha$first_stage_null_calls == 0L),
    all(is.finite(smoke_alpha$power)),
    all(is.finite(smoke_alpha$empirical_fsr))
  )
  message("Ideal-Gaussian 60-unit integration smoke test passed.")
}

production_files <- c(
  "run_r3_ideal_gaussian_measurement.R",
  "23_r3_ideal_gaussian_measurement.sbatch",
  "reporting.R"
)
production_paths <- file.path(dirname(helper_path), production_files)
if (all(file.exists(production_paths))) {
  runner_text <- paste(readLines(production_paths[[1L]], warn = FALSE), collapse = "\n")
  sbatch_text <- paste(readLines(production_paths[[2L]], warn = FALSE), collapse = "\n")
  reporting_text <- paste(readLines(production_paths[[3L]], warn = FALSE), collapse = "\n")
  stopifnot(
    grepl("ideal_gaussian_known_t_adjusted_se", runner_text, fixed = TRUE),
    grepl('candidate_scope = "full_universe"', runner_text, fixed = TRUE),
    grepl("functional_posterior", runner_text, fixed = TRUE),
    grepl("complete.flag", runner_text, fixed = TRUE),
    grepl("--preflight-only", sbatch_text, fixed = TRUE),
    grepl("#SBATCH --cpus-per-task=16", sbatch_text, fixed = TRUE),
    grepl("FASH-IWP1-BF", reporting_text, fixed = TRUE)
  )
}

message("Ideal-Gaussian measurement unit tests passed.")
