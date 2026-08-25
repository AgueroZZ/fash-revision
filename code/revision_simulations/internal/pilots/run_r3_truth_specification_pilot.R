#!/usr/bin/env Rscript

# Run a deliberately exploratory real-genotype R3 truth-design pilot. This is
# not a formal replacement for the frozen R3 cache: its role is to compare two
# prespecified, more sharply labelled truth geometries in the same FASH
# pipeline before any full five-seed refresh is proposed.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not locate simulation_functions.R.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  position <- match(name, args)
  if (is.na(position) || position == length(args)) return(default)
  args[[position + 1L]]
}

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE)

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
))

truth_mechanism <- match.arg(
  get_arg("--truth-mechanism", "r3b_interior_single_lobe"),
  c(
    "r3b_interior_single_lobe",
    "r3a_targeted_broad_bspline",
    "r3a_constrained_global_bspline"
  )
)
seed <- as.integer(get_arg("--seed", "12345"))
J <- as.integer(get_arg("--J", "6362"))
num_cores <- as.integer(get_arg("--num-cores", "8"))
output_id <- get_arg(
  "--output-id",
  paste0("r3_truth_specification_", truth_mechanism, "_seed", seed)
)
genotype_cache_path <- get_arg(
  "--genotype-cache",
  file.path(
    workflowr_root,
    "..", "workspace", "inputs", "r1_r2_fashr0143", "genotype_samples.rds"
  )
)

if (!is.finite(seed) || J < 30L || !is.finite(num_cores) || num_cores < 1L ||
    !nzchar(output_id) || !file.exists(genotype_cache_path)) {
  stop("Invalid pilot arguments or missing genotype cache.")
}
if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The pinned fashr package is required for this pilot.")
}

fashr_description <- utils::packageDescription("fashr")
if (!identical(as.character(utils::packageVersion("fashr")), "0.1.43") ||
    !identical(
      as.character(fashr_description$RemoteSha),
      "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
    )) {
  stop("The pilot requires fashr 0.1.43 at the pinned revision SHA.")
}

genotype_cache <- readRDS(genotype_cache_path)
if (!all(c("configuration", "samples") %in% names(genotype_cache)) ||
    !as.character(seed) %in% names(genotype_cache$samples)) {
  stop("The requested seed is absent from the genotype cache.")
}
genotype_sample <- validate_real_genotype_sample(
  genotype_cache$samples[[as.character(seed)]],
  expected_genes = genotype_cache$configuration$n_genes,
  expected_donors = genotype_cache$configuration$n_donors,
  maf_min = genotype_cache$configuration$maf_min
)
if (J > ncol(genotype_sample$G)) {
  stop("J exceeds the available number of real-genotype gene--variant pairs.")
}
if (J < ncol(genotype_sample$G)) {
  keep <- seq_len(J)
  genotype_sample$G <- genotype_sample$G[, keep, drop = FALSE]
  genotype_sample$variant_info <- genotype_sample$variant_info[keep, , drop = FALSE]
  genotype_sample$selection <- genotype_sample$selection[keep, , drop = FALSE]
}

component_seeds <- revision_component_seeds(seed)
time_grid <- 0:15
evaluation_grid <- seq(0, 15, by = 0.1)
middle_window <- c(3, 12)
class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
targets <- c("early", "middle", "late", "switch")
covariates <- simulate_covariate_matrix(
  n_donors = nrow(genotype_sample$G),
  n_covariates = 5L,
  seed = component_seeds[["covariates"]]
)

truth_settings <- if (truth_mechanism == "r3b_interior_single_lobe") {
  list(
    truth_mechanism = "raised_cosine",
    primary_center_ranges = list(
      early = c(1.5, 1.5),
      middle = c(4.5, 10.5),
      late = c(13.5, 13.5)
    ),
    cosine_half_width = 1.2,
    cosine_spike_counts = 1L,
    label = "strictly target-contained single-lobe raised cosine"
  )
} else if (truth_mechanism == "r3a_targeted_broad_bspline") {
  list(
    truth_mechanism = "targeted_broad_bspline",
    minimum_location_margin = 0.60,
    minimum_location_ratio = 1.40,
    non_switch_baseline_fraction = 0.75,
    non_switch_background_fraction = 0.05,
    label = "functionally separated broad cubic B-spline"
  )
} else {
  list(
    truth_mechanism = "constrained_global_bspline",
    minimum_location_margin = 0.60,
    minimum_location_ratio = 1.40,
    non_switch_baseline_fraction = 0.75,
    target_centered_rms = 0.90,
    label = "functionally separated global random cubic B-spline"
  )
}

effect_sim <- if (truth_mechanism == "r3b_interior_single_lobe") {
  simulate_matched_functional_effect_set(
    n_variants = J,
    truth_mechanism = "raised_cosine",
    time_grid = time_grid,
    evaluation_grid = evaluation_grid,
    class_probs = class_probs,
    dynamic_main_effect_sd = 1,
    cosine_width_half = truth_settings$cosine_half_width,
    cosine_spike_counts = truth_settings$cosine_spike_counts,
    cosine_center_ranges = truth_settings$primary_center_ranges,
    switch_threshold = 0.25,
    location_truth_margin = 0.10,
    location_truth_min_range_fraction = 0.10,
    switch_truth_margin = 0.10,
    non_switch_min_abs = 0.10,
    non_switch_min_range_fraction = 0.10,
    seed = seed,
    class_seed = component_seeds[["classes"]],
    constant_seed = component_seeds[["constant_effects"]],
    shape_seed = component_seeds[["functional_truth"]],
    middle_window = middle_window,
    middle_boundary = "open"
  )
} else if (truth_mechanism == "r3a_targeted_broad_bspline") {
  simulate_targeted_local_bspline_effect_set(
    n_variants = J,
    time_grid = time_grid,
    evaluation_grid = evaluation_grid,
    class_probs = class_probs,
    dynamic_amplitude = 2,
    switch_threshold = 0.25,
    minimum_location_margin = truth_settings$minimum_location_margin,
    minimum_location_ratio = truth_settings$minimum_location_ratio,
    non_switch_baseline_fraction = truth_settings$non_switch_baseline_fraction,
    non_switch_background_fraction = truth_settings$non_switch_background_fraction,
    profile = "broad",
    seed = component_seeds[["functional_truth"]],
    scenario = "r3a_targeted_broad_real_genotype_pilot",
    middle_window = middle_window,
    middle_boundary = "open"
  )
} else {
  simulate_targeted_local_bspline_effect_set(
    n_variants = J,
    time_grid = time_grid,
    evaluation_grid = evaluation_grid,
    class_probs = class_probs,
    dynamic_amplitude = 2,
    switch_threshold = 0.25,
    minimum_location_margin = truth_settings$minimum_location_margin,
    minimum_location_ratio = truth_settings$minimum_location_ratio,
    non_switch_baseline_fraction = truth_settings$non_switch_baseline_fraction,
    target_centered_rms = truth_settings$target_centered_rms,
    profile = "random_broad",
    seed = component_seeds[["functional_truth"]],
    scenario = "r3a_constrained_global_bspline_real_genotype_pilot",
    middle_window = middle_window,
    middle_boundary = "open"
  )
}

if (!identical(effect_sim$settings$middle_window, middle_window) ||
    !identical(effect_sim$settings$middle_boundary, "open")) {
  stop("Truth generation and posterior evaluation use different Middle estimands.")
}

recomputed_true_functionals <- evaluate_temporal_functionals(
  effect_sim$beta_evaluation,
  smooth_var = effect_sim$evaluation_grid,
  switch_threshold = 0.25,
  middle_window = middle_window,
  middle_boundary = "open"
)
if (!isTRUE(all.equal(
  effect_sim$true_functionals,
  recomputed_true_functionals,
  check.attributes = TRUE
))) {
  stop("Saved truth functionals do not match the requested Middle estimand.")
}

effect_sim <- reassign_effect_simulation_by_maf(
  effect_sim = effect_sim,
  maf = genotype_sample$variant_info$observed_maf,
  class_probs = class_probs,
  seed = component_seeds[["classes"]],
  n_strata = 10L
)
for (field in c("beta_matrix", "beta_evaluation", "true_functionals")) {
  rownames(effect_sim[[field]]) <- genotype_sample$selection$pair_key
}
effect_sim$unit_info$variant_id <- genotype_sample$selection$pair_key

expression_sim <- simulate_eqtl_expression_from_genotypes(
  G = genotype_sample$G,
  beta_matrix = effect_sim$beta_matrix,
  time_grid = time_grid,
  covariates = covariates,
  expression_noise_sd = 1,
  covariate_effect_sd = 0.5,
  intercept_sd = 0,
  seed = component_seeds[["expression"]]
)
output_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal", output_id
)
if (dir.exists(output_dir) && length(list.files(
  output_dir,
  all.files = TRUE,
  no.. = TRUE
)) > 0L) {
  stop("Refusing to overwrite a non-empty pilot output directory: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
out <- run_genotype_level_dynamic_eqtl_simulation(
  G = genotype_sample$G,
  time_grid = time_grid,
  covariates = covariates,
  class_probs = class_probs,
  expression_noise_sd = 1,
  covariate_effect_sd = 0.5,
  intercept_sd = 0,
  dynamic_main_effect_sd = 1,
  alpha = 0.05,
  seed = seed,
  num_cores = num_cores,
  num_basis = 20L,
  scenario = output_id,
  output_dir = output_dir,
  save_outputs = FALSE,
  verbose = FALSE,
  effect_sim = effect_sim,
  expression_sim = expression_sim
)

true_functionals <- effect_sim$true_functionals[, targets, drop = FALSE]
true_dynamic <- out$unit_info$effect_class == "dynamic_bspline"
functional_seed <- component_seeds[["functional_posterior"]]
raw <- evaluate_fash_functional_testing(
  fit = out$fash_fits$fash_iwp1_raw,
  true_functionals = true_functionals,
  evaluation_grid = evaluation_grid,
  method = "FASH-IWP1-Raw",
  scenario = output_id,
  switch_threshold = 0.25,
  true_dynamic = true_dynamic,
  num_cores = num_cores,
  seed = functional_seed,
  middle_window = middle_window,
  middle_boundary = "open"
)
bf <- evaluate_fash_functional_testing(
  fit = out$fash_fits$fash_iwp1_bf,
  true_functionals = true_functionals,
  evaluation_grid = evaluation_grid,
  method = "FASH-IWP1-BF",
  scenario = output_id,
  switch_threshold = 0.25,
  true_dynamic = true_dynamic,
  num_cores = num_cores,
  seed = functional_seed,
  middle_window = middle_window,
  middle_boundary = "open"
)

configuration <- list(
  output_id = output_id,
  exploratory = TRUE,
  seed = seed,
  J = J,
  genotype_source = "paper-derived YRI DS dosage",
  genotype_selection_rule = genotype_cache$configuration$selection_rule,
  truth_settings = truth_settings,
  truth_generation_settings = effect_sim$settings,
  middle_window = middle_window,
  middle_boundary = "open",
  middle_definition = "3 < t < 12",
  package_provenance = list(
    version = as.character(utils::packageVersion("fashr")),
    remote_sha = as.character(fashr_description$RemoteSha)
  )
)
saveRDS(configuration, file.path(output_dir, "configuration.rds"))
saveRDS(
  list(
    configuration = configuration,
    functional_alpha = rbind(raw$alpha_curve, bf$alpha_curve),
    truth_group_counts = as.data.frame(table(effect_sim$unit_info$truth_group)),
    truth_maf_balance = summarize_truth_maf_balance(
      genotype_sample$variant_info,
      out$unit_info,
      seed = seed
    )
  ),
  file.path(output_dir, "pilot_result.rds")
)
write_csv(
  rbind(raw$alpha_curve, bf$alpha_curve),
  file.path(output_dir, "functional_alpha_curve.csv")
)
message("Saved exploratory real-genotype R3 truth-design pilot to: ", output_dir)
