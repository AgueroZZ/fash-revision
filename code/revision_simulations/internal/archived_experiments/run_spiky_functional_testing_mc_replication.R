#!/usr/bin/env Rscript

# Run functional-testing replications for the reviewer-facing sparse
# raised-cosine one-/two-peak genotype simulation.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) return(".")
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find code/revision_simulations/shared/simulation_functions.R.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0 || hit[1] == length(args)) return(default)
  args[hit[1] + 1]
}

as_flag <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

parse_seed_list <- function(x) {
  seeds <- suppressWarnings(as.integer(trimws(strsplit(x, ",", fixed = TRUE)[[1]])))
  if (length(seeds) == 0 || anyNA(seeds) || anyDuplicated(seeds)) {
    stop("--seed-list must contain unique comma-separated integer seeds.")
  }
  seeds
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

J <- as.integer(get_arg("--J", "1000"))
n_donors <- as.integer(get_arg("--n-donors", "19"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
expression_noise_sd <- as.numeric(get_arg("--noise-sd", "1"))
width_half <- as.numeric(get_arg("--width-half", "1.5"))
target_centered_rms <- as.numeric(get_arg("--target-centered-rms", "0.9"))
dynamic_main_effect_sd <- as.numeric(get_arg(
  "--dynamic-main-effect-sd",
  "0"
))
non_switch_min_abs_effect <- as.numeric(get_arg(
  "--non-switch-min-abs-effect",
  "0"
))
switch_effect_threshold <- as.numeric(get_arg("--switch-effect-threshold", "0.25"))
switch_minimum_duration <- as.numeric(get_arg("--switch-minimum-duration", "0"))
num_basis <- as.integer(get_arg("--num-basis", "20"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
seed_list <- parse_seed_list(get_arg(
  "--seed-list",
  "12345,22345,32345,42345,52345"
))
output_id <- get_arg(
  "--output-id",
  "sparse_timed_cosine_functional_pilot5"
)
overwrite <- as_flag(get_arg("--overwrite", "false"))
cache_full_fits <- as_flag(get_arg("--cache-full-fits", "true"))

if (J < 10 || n_donors < n_covariates + 3 || n_covariates < 0 ||
    !is.finite(expression_noise_sd) || expression_noise_sd <= 0 ||
    !is.finite(width_half) || width_half <= 0 ||
    !is.finite(target_centered_rms) || target_centered_rms <= 0 ||
    !is.finite(dynamic_main_effect_sd) || dynamic_main_effect_sd < 0 ||
    !is.finite(non_switch_min_abs_effect) ||
    non_switch_min_abs_effect < 0 ||
    (dynamic_main_effect_sd > 0 && non_switch_min_abs_effect > 0) ||
    !is.finite(switch_effect_threshold) || switch_effect_threshold <= 0 ||
    !is.finite(switch_minimum_duration) || switch_minimum_duration < 0 ||
    num_basis < 2 || num_cores < 1 || !nzchar(output_id)) {
  stop("Invalid sparse-cosine functional-testing arguments.")
}

time_grid <- make_time_grid()
evaluation_grid <- seq(min(time_grid), max(time_grid), by = 0.1)
alpha_grid <- seq(0.005, 0.20, by = 0.005)
class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
shape_cell_probs <- c(
  k1__spiky__single = 0.50,
  `k2__spiky__same-sign` = 0.25,
  `k2__spiky__alternating-sign` = 0.25
)
primary_time_groups <- c("early", "middle", "late")
relative_amplitude_range <- c(0.35, 0.60)
center_by_observed_mean <- FALSE
truth_generation_switch_threshold <- 0.25
scenario <- "genotype_sparse_timed_cosine_one_two_peak_dynamic_eqtl"
methods <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
targets <- c("early", "middle", "late", "switch")
true_pi0 <- unname(class_probs["constant"] + class_probs["zero"])

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "mc",
  output_id
)
replicate_dir <- file.path(output_dir, "replicates")
summary_dir <- file.path(output_dir, "summary")
figure_dir <- file.path(output_dir, "figures")
invisible(lapply(
  c(output_dir, replicate_dir, summary_dir, figure_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

configuration <- list(
  output_id = output_id,
  scenario = scenario,
  J = J,
  n_donors = n_donors,
  n_covariates = n_covariates,
  time_grid = time_grid,
  evaluation_grid = evaluation_grid,
  expression_noise_sd = expression_noise_sd,
  width_half = width_half,
  target_centered_rms = target_centered_rms,
  num_basis = num_basis,
  class_probs = class_probs,
  spike_counts = 1:2,
  shape_cell_probs = shape_cell_probs,
  shape_cell_counts = J * class_probs[["dynamic_bspline"]] * shape_cell_probs,
  primary_time_groups = primary_time_groups,
  relative_amplitude_range = relative_amplitude_range,
  center_by_observed_mean = center_by_observed_mean,
  truth_generation_switch_threshold = truth_generation_switch_threshold,
  switch_effect_threshold = switch_effect_threshold,
  switch_minimum_duration = switch_minimum_duration,
  constant_sd = 1,
  dynamic_baseline_sd = 0,
  true_pi0 = true_pi0,
  methods = methods,
  targets = targets,
  alpha_grid = alpha_grid,
  seed_list = seed_list
)
if (non_switch_min_abs_effect > 0) {
  configuration$non_switch_min_abs_effect <- non_switch_min_abs_effect
}
if (dynamic_main_effect_sd > 0) {
  configuration$dynamic_main_effect_sd <- dynamic_main_effect_sd
}

configuration_matches <- function(candidate) {
  if (isTRUE(all.equal(candidate, configuration))) {
    return(TRUE)
  }
  if (non_switch_min_abs_effect > 0 ||
      dynamic_main_effect_sd > 0 ||
      !isTRUE(all.equal(
        truth_generation_switch_threshold,
        switch_effect_threshold
      )) ||
      switch_minimum_duration != 0) {
    return(FALSE)
  }
  legacy_configuration <- configuration[
    setdiff(
      names(configuration),
      c(
        "truth_generation_switch_threshold",
        "switch_effect_threshold",
        "switch_minimum_duration"
      )
    )
  ]
  legacy_configuration$switch_threshold <- switch_effect_threshold
  legacy_order <- c(
    names(legacy_configuration)[
      seq_len(match("center_by_observed_mean", names(legacy_configuration)))
    ],
    "switch_threshold",
    names(legacy_configuration)[
      (match("center_by_observed_mean", names(legacy_configuration)) + 1):
        (length(legacy_configuration) - 1)
    ]
  )
  legacy_configuration <- legacy_configuration[legacy_order]
  isTRUE(all.equal(candidate, legacy_configuration))
}

configuration_path <- file.path(output_dir, "configuration.rds")
if (file.exists(configuration_path) && !overwrite) {
  cached_configuration <- readRDS(configuration_path)
  if (!configuration_matches(cached_configuration)) {
    stop("The existing output id has different settings.")
  }
} else {
  saveRDS(configuration, configuration_path)
}

validate_truth <- function(effect_sim, seed) {
  dynamic <- effect_sim$unit_info$effect_class == "dynamic_bspline"
  pattern_switch <- dynamic &
    effect_sim$unit_info$sign_pattern == "alternating-sign"
  pattern_non_switch <- dynamic & !pattern_switch
  realized_switch <- effect_sim$true_functionals[dynamic, "switch"] > 0
  observed_cell_counts <- table(effect_sim$unit_info$cell_id[dynamic])
  timing_counts <- table(
    effect_sim$unit_info$cell_id[dynamic],
    effect_sim$unit_info$time_group[dynamic]
  )
  zero_fraction <- rowMeans(
    abs(effect_sim$beta_evaluation[dynamic, , drop = FALSE]) < 1e-12
  )
  minimum_zero_fraction <- max(
    0,
    1 - 4 * width_half / diff(range(evaluation_grid)) - 0.05
  )
  centered_rms <- sqrt(rowMeans(
    (
      effect_sim$beta_matrix[dynamic, , drop = FALSE] -
        rowMeans(effect_sim$beta_matrix[dynamic, , drop = FALSE])
    )^2
  ))
  non_switch_indices <- which(pattern_non_switch)
  non_switch_direction <- vapply(
    effect_sim$unit_info$peak_signs[non_switch_indices],
    function(x) sign(x[1]),
    numeric(1)
  )
  non_switch_min_signed_effect <- vapply(
    seq_along(non_switch_indices),
    function(position) {
      index <- non_switch_indices[position]
      min(non_switch_direction[position] * effect_sim$beta_evaluation[index, ])
    },
    numeric(1)
  )
  support_valid <- if (dynamic_main_effect_sd > 0) {
    any(realized_switch) &&
      any(!realized_switch) &&
      all(is.finite(effect_sim$unit_info$baseline[dynamic]))
  } else if (non_switch_min_abs_effect > 0) {
    all(zero_fraction[pattern_switch[dynamic]] > minimum_zero_fraction) &&
      all(
        non_switch_min_signed_effect >=
          non_switch_min_abs_effect - 1e-10
      )
  } else {
    all(zero_fraction > minimum_zero_fraction)
  }
  valid <- sum(dynamic) == J * class_probs[["dynamic_bspline"]] &&
    identical(
      as.integer(observed_cell_counts[names(shape_cell_probs)]),
      as.integer(configuration$shape_cell_counts)
    ) &&
    setequal(colnames(timing_counts), primary_time_groups) &&
    all(apply(timing_counts, 1, function(x) diff(range(x))) <= 1L) &&
    support_valid &&
    all(abs(centered_rms - target_centered_rms) <= 1e-10)
  if (!valid) {
    stop("Invalid sparse timed cosine truth for seed ", seed, ".")
  }
  invisible(TRUE)
}

add_dynamic_main_effects <- function(effect_sim, seed) {
  if (dynamic_main_effect_sd <= 0) {
    return(effect_sim)
  }
  dynamic <- effect_sim$unit_info$effect_class == "dynamic_bspline"
  max_seed <- .Machine$integer.max - 1
  main_effect_seed <- as.integer(((seed + 2309 - 1) %% max_seed) + 1)
  set.seed(main_effect_seed)
  main_effect <- stats::rnorm(
    sum(dynamic),
    mean = 0,
    sd = dynamic_main_effect_sd
  )
  effect_sim$unit_info$dynamic_main_effect <- 0
  effect_sim$unit_info$dynamic_main_effect[dynamic] <- main_effect
  effect_sim$beta_matrix[dynamic, ] <-
    effect_sim$beta_matrix[dynamic, , drop = FALSE] + main_effect
  effect_sim$beta_evaluation[dynamic, ] <-
    effect_sim$beta_evaluation[dynamic, , drop = FALSE] + main_effect
  effect_sim$unit_info$baseline[dynamic] <-
    effect_sim$unit_info$baseline[dynamic] + main_effect
  effect_sim$true_functionals <- evaluate_temporal_functionals(
    effect_sim$beta_evaluation,
    smooth_var = evaluation_grid,
    switch_threshold = truth_generation_switch_threshold
  )
  effect_sim$unit_info$switch_status[dynamic] <- ifelse(
    effect_sim$true_functionals[dynamic, "switch"] > 0,
    "switch",
    "non-switch"
  )
  effect_sim$settings$dynamic_main_effect_sd <- dynamic_main_effect_sd
  effect_sim$settings$dynamic_main_effect_seed <- main_effect_seed
  effect_sim
}

separate_non_switch_curves <- function(effect_sim) {
  if (non_switch_min_abs_effect <= 0) {
    return(effect_sim)
  }
  dynamic_non_switch <- effect_sim$unit_info$effect_class == "dynamic_bspline" &
    effect_sim$unit_info$sign_pattern != "alternating-sign"
  for (index in which(dynamic_non_switch)) {
    direction <- sign(effect_sim$unit_info$peak_signs[[index]][1])
    current_minimum <- min(direction * effect_sim$beta_evaluation[index, ])
    offset <- max(0, non_switch_min_abs_effect - current_minimum)
    effect_sim$beta_matrix[index, ] <-
      effect_sim$beta_matrix[index, ] + direction * offset
    effect_sim$beta_evaluation[index, ] <-
      effect_sim$beta_evaluation[index, ] + direction * offset
    effect_sim$unit_info$baseline[index] <-
      effect_sim$unit_info$baseline[index] + direction * offset
  }
  effect_sim$true_functionals <- evaluate_temporal_functionals(
    effect_sim$beta_evaluation,
    smooth_var = evaluation_grid,
    switch_threshold = truth_generation_switch_threshold
  )
  effect_sim
}

simulate_fash_output <- function(seed) {
  component_seeds <- revision_component_seeds(seed)
  genotype_sim <- simulate_genotype_matrix(
    n_donors = n_donors,
    n_variants = J,
    seed = component_seeds[["genotype"]]
  )
  covariates <- simulate_covariate_matrix(
    n_donors = n_donors,
    n_covariates = n_covariates,
    seed = component_seeds[["covariates"]]
  )
  effect_sim <- simulate_raised_cosine_multipeak_effect_set(
    n_variants = J,
    time_grid = time_grid,
    evaluation_grid = evaluation_grid,
    class_probs = class_probs,
    width_levels = c(spiky = width_half),
    spike_counts = 1:2,
    shape_cell_probs = shape_cell_probs,
    primary_time_groups = primary_time_groups,
    center_by_observed_mean = center_by_observed_mean,
    switch_threshold = truth_generation_switch_threshold,
    relative_amplitude_range = relative_amplitude_range,
    target_centered_rms = target_centered_rms,
    baseline_sd = 1,
    constant_sd = 1,
    dynamic_baseline_sd = 0,
    exact_class_counts = TRUE,
    seed = seed,
    class_seed = component_seeds[["classes"]],
    constant_seed = component_seeds[["constant_effects"]],
    shape_seed = component_seeds[["functional_truth"]],
    scenario = scenario
  )
  effect_sim <- add_dynamic_main_effects(effect_sim, seed)
  effect_sim <- separate_non_switch_curves(effect_sim)
  validate_truth(effect_sim, seed)
  expression_sim <- simulate_eqtl_expression_from_genotypes(
    G = genotype_sim$G,
    beta_matrix = effect_sim$beta_matrix,
    time_grid = time_grid,
    covariates = covariates,
    expression_noise_sd = expression_noise_sd,
    seed = component_seeds[["expression"]]
  )
  out <- run_genotype_level_dynamic_eqtl_simulation(
    G = genotype_sim$G,
    time_grid = time_grid,
    covariates = covariates,
    class_probs = class_probs,
    expression_noise_sd = expression_noise_sd,
    alpha = 0.05,
    seed = seed,
    num_cores = num_cores,
    num_basis = num_basis,
    scenario = scenario,
    output_dir = output_dir,
    save_outputs = FALSE,
    verbose = FALSE,
    effect_sim = effect_sim,
    expression_sim = expression_sim
  )
  out$true_functionals <- effect_sim$true_functionals
  out$true_beta_evaluation <- effect_sim$beta_evaluation
  out$evaluation_grid <- evaluation_grid
  out$component_seeds <- component_seeds
  out
}

load_or_simulate_fash_output <- function(seed) {
  raw_path <- if (
    non_switch_min_abs_effect > 0 ||
      dynamic_main_effect_sd > 0
  ) {
    file.path(output_dir, "full_fits", paste0("seed_", seed, ".rds"))
  } else {
    stem <- genotype_cosine_multipeak_output_stem(
      n_donors = n_donors,
      n_variants = J,
      time_grid = time_grid,
      n_covariates = n_covariates,
      expression_noise_sd = expression_noise_sd,
      width_half = width_half,
      target_centered_rms = target_centered_rms,
      shape_cell_probs = shape_cell_probs,
      class_probs = class_probs,
      seed = seed,
      scenario = scenario
    )
    file.path(
      workflowr_root,
      "output",
      "revision_simulations",
      "raw",
      paste0(stem, ".rds")
    )
  }
  if (file.exists(raw_path)) {
    out <- readRDS(raw_path)
    expected_fields <- c(
      "unit_info", "fash_fits", "true_functionals", "true_beta_evaluation",
      "evaluation_grid", "component_seeds"
    )
    if (all(expected_fields %in% names(out)) &&
        identical(out$settings$scenario, scenario) &&
        isTRUE(all.equal(out$settings$n_variants, J)) &&
        isTRUE(all.equal(out$settings$n_donors, n_donors)) &&
        isTRUE(all.equal(out$settings$n_covariates, n_covariates)) &&
        isTRUE(all.equal(out$evaluation_grid, evaluation_grid))) {
      message("Reusing full fitted output: ", raw_path)
      return(out)
    }
    stop("The existing full fitted output does not match the requested setting.")
  }
  message("Simulating and fitting sparse-cosine replicate with seed ", seed, ".")
  out <- simulate_fash_output(seed)
  if (cache_full_fits) {
    dir.create(dirname(raw_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(out, raw_path)
    message("Cached full fitted output: ", raw_path)
  }
  out
}

make_replicate <- function(seed) {
  out <- load_or_simulate_fash_output(seed)
  true_curves <- out$true_beta_evaluation[
    out$unit_info$variant_id,
    ,
    drop = FALSE
  ]
  true_functionals <- evaluate_temporal_functionals(
    curves = true_curves,
    smooth_var = evaluation_grid,
    switch_threshold = switch_effect_threshold,
    switch_minimum_duration = switch_minimum_duration
  )[, targets, drop = FALSE]
  true_dynamic <- out$unit_info$effect_class == "dynamic_bspline"
  component_seeds <- revision_component_seeds(seed)
  raw_results <- evaluate_fash_functional_testing(
    fit = out$fash_fits$fash_iwp1_raw,
    true_functionals = true_functionals,
    evaluation_grid = evaluation_grid,
    alpha_grid = alpha_grid,
    method = "FASH-IWP1-Raw",
    scenario = scenario,
    switch_threshold = switch_effect_threshold,
    switch_minimum_duration = switch_minimum_duration,
    true_dynamic = true_dynamic,
    num_cores = num_cores,
    seed = component_seeds[["functional_posterior"]]
  )
  bf_results <- evaluate_fash_functional_testing(
    fit = out$fash_fits$fash_iwp1_bf,
    true_functionals = true_functionals,
    evaluation_grid = evaluation_grid,
    alpha_grid = alpha_grid,
    method = "FASH-IWP1-BF",
    scenario = scenario,
    switch_threshold = switch_effect_threshold,
    switch_minimum_duration = switch_minimum_duration,
    true_dynamic = true_dynamic,
    num_cores = num_cores,
    seed = component_seeds[["functional_posterior"]] + 100L
  )
  functional_alpha <- rbind(raw_results$alpha_curve, bf_results$alpha_curve)
  functional_alpha$seed <- seed
  truth_counts <- data.frame(
    seed = seed,
    target = targets,
    n_true_alternatives = colSums(true_functionals[, targets, drop = FALSE] > 0),
    stringsAsFactors = FALSE
  )
  list(
    configuration = configuration,
    seed = seed,
    component_seeds = component_seeds,
    functional_alpha = functional_alpha,
    functional_alpha_005 = functional_alpha[
      abs(functional_alpha$alpha - 0.05) < 1e-12,
      ,
      drop = FALSE
    ],
    truth_counts = truth_counts,
    estimated_pi0 = data.frame(
      seed = seed,
      method = c("FASH-IWP1", "FASH-IWP1"),
      fit = c("Raw", "BF-corrected"),
      estimated_pi0 = c(
        constant_component_prior_weight(out$fash_fits$fash_iwp1_raw),
        constant_component_prior_weight(out$fash_fits$fash_iwp1_bf)
      ),
      stringsAsFactors = FALSE
    )
  )
}

validate_replicate <- function(x, seed) {
  required <- c(
    "configuration", "seed", "component_seeds", "functional_alpha",
    "functional_alpha_005", "truth_counts", "estimated_pi0"
  )
  all(required %in% names(x)) &&
    identical(x$seed, seed) &&
    configuration_matches(x$configuration) &&
    nrow(x$functional_alpha) == length(methods) * length(targets) * length(alpha_grid) &&
    nrow(x$functional_alpha_005) == length(methods) * length(targets) &&
    all(methods %in% unique(x$functional_alpha$method)) &&
    all(targets %in% unique(x$functional_alpha$target)) &&
    nrow(x$truth_counts) == length(targets) &&
    nrow(x$estimated_pi0) == length(methods)
}

replicates <- lapply(seed_list, function(seed) {
  path <- file.path(replicate_dir, paste0("seed_", seed, ".rds"))
  if (file.exists(path) && !overwrite) {
    cached <- readRDS(path)
    if (validate_replicate(cached, seed)) {
      message("Reusing sparse-cosine functional replicate: ", path)
      return(cached)
    }
    stop("Cached functional replicate does not match: ", path)
  }
  replicate <- make_replicate(seed)
  saveRDS(replicate, path)
  replicate
})

all_alpha <- do.call(rbind, lapply(replicates, `[[`, "functional_alpha"))
all_alpha_005 <- do.call(rbind, lapply(replicates, `[[`, "functional_alpha_005"))
all_truth_counts <- do.call(rbind, lapply(replicates, `[[`, "truth_counts"))
all_pi0 <- do.call(rbind, lapply(replicates, `[[`, "estimated_pi0"))
mc_alpha <- summarize_mc_functional_alpha_curves(all_alpha)
mc_alpha_005 <- mc_alpha[abs(mc_alpha$alpha - 0.05) < 1e-12, , drop = FALSE]
mc_pi0 <- summarize_mc_pi0(all_pi0)

if (any(mc_alpha$n_replications != length(seed_list)) ||
    any(mc_alpha$power_ci_lower < 0, na.rm = TRUE) ||
    any(mc_alpha$power_ci_upper > 1, na.rm = TRUE) ||
    any(mc_alpha$empirical_fsr_ci_lower < 0, na.rm = TRUE) ||
    any(mc_alpha$empirical_fsr_ci_upper > 1, na.rm = TRUE)) {
  stop("Sparse-cosine functional Monte Carlo summaries failed validation.")
}

write_csv(
  all_alpha,
  file.path(summary_dir, "all_replicate_functional_alpha_curves.csv")
)
write_csv(
  all_alpha_005,
  file.path(summary_dir, "all_replicate_functional_alpha005.csv")
)
write_csv(
  all_truth_counts,
  file.path(summary_dir, "all_replicate_truth_counts.csv")
)
write_csv(all_pi0, file.path(summary_dir, "all_replicate_pi0.csv"))
write_csv(
  mc_alpha,
  file.path(summary_dir, "functional_testing_mc_alpha_curve.csv")
)
write_csv(
  mc_alpha_005,
  file.path(summary_dir, "functional_testing_mc_alpha005_summary.csv")
)
write_csv(
  mc_pi0,
  file.path(summary_dir, "functional_testing_mc_pi0_summary.csv")
)

plot_mc_functional_curve_grid(
  mc_curve = mc_alpha,
  metric = "power",
  file = file.path(figure_dir, "functional_testing_mc_power.png"),
  title = "Functional-testing power for sparse transient effects"
)
plot_mc_functional_curve_grid(
  mc_curve = mc_alpha,
  metric = "empirical_fsr",
  file = file.path(figure_dir, "functional_testing_mc_empirical_fsr.png"),
  title = "Empirical FSR for sparse transient effects"
)

print(mc_alpha_005)
print(mc_pi0)
