#!/usr/bin/env Rscript

# Run a paired internal sensitivity analysis for the sparse raised-cosine
# genotype simulation. This script does not build or modify workflowr pages.

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
source_centered_rms <- as.numeric(get_arg("--source-centered-rms", "0.9"))
target_max_abs <- as.numeric(get_arg("--target-max-abs", "1"))
current_switch_threshold <- as.numeric(get_arg("--current-switch-threshold", "0.25"))
relative_switch_threshold <- as.numeric(get_arg("--relative-switch-threshold", "0.75"))
num_basis <- as.integer(get_arg("--num-basis", "20"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
efdr_permutations <- as.integer(get_arg("--efdr-permutations", "100"))
posterior_samples <- as.integer(get_arg("--posterior-samples", "3000"))
seed_list <- parse_seed_list(get_arg(
  "--seed-list",
  "12345,22345,32345,42345,52345"
))
output_id <- get_arg("--output-id", "spiky_maxabs_switch_sensitivity")
overwrite <- as_flag(get_arg("--overwrite", "false"))

if (J < 10 || n_donors < n_covariates + 3 || n_covariates < 0 ||
    !is.finite(expression_noise_sd) || expression_noise_sd <= 0 ||
    !is.finite(width_half) || width_half <= 0 ||
    !is.finite(source_centered_rms) || source_centered_rms <= 0 ||
    !is.finite(target_max_abs) || target_max_abs <= 0 ||
    !is.finite(current_switch_threshold) || current_switch_threshold <= 0 ||
    !is.finite(relative_switch_threshold) || relative_switch_threshold <= 0 ||
    num_basis < 2 || num_cores < 1 || efdr_permutations < 1 ||
    posterior_samples < 500 || !nzchar(output_id)) {
  stop("Invalid max-absolute-effect sensitivity arguments.")
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
source_scenario <- "genotype_sparse_timed_cosine_one_two_peak_dynamic_eqtl"
maxabs_scenario <- "internal_genotype_sparse_timed_cosine_maxabs1_dynamic_eqtl"
true_pi0 <- unname(class_probs["constant"] + class_probs["zero"])
global_methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "Direct-linear-LRT-eFDR-true-pi0",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
functional_methods <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")

condition_labels <- c(
  current = "current_scale_threshold_0.25",
  relative_threshold = "current_scale_threshold_0.75",
  maxabs = "max_abs_1_threshold_0.25"
)

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  output_id
)
fit_dir <- file.path(output_dir, "full_fits")
replicate_dir <- file.path(output_dir, "replicates")
summary_dir <- file.path(output_dir, "summary")
figure_dir <- file.path(output_dir, "figures")
invisible(lapply(
  c(output_dir, fit_dir, replicate_dir, summary_dir, figure_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

configuration <- list(
  output_id = output_id,
  J = J,
  n_donors = n_donors,
  n_covariates = n_covariates,
  time_grid = time_grid,
  evaluation_grid = evaluation_grid,
  expression_noise_sd = expression_noise_sd,
  width_half = width_half,
  source_centered_rms = source_centered_rms,
  target_max_abs = target_max_abs,
  current_switch_threshold = current_switch_threshold,
  relative_switch_threshold = relative_switch_threshold,
  num_basis = num_basis,
  efdr_permutations = efdr_permutations,
  posterior_samples = posterior_samples,
  class_probs = class_probs,
  shape_cell_probs = shape_cell_probs,
  primary_time_groups = primary_time_groups,
  relative_amplitude_range = relative_amplitude_range,
  center_by_observed_mean = FALSE,
  dynamic_baseline_sd = 0,
  constant_sd = 1,
  true_pi0 = true_pi0,
  alpha_grid = alpha_grid,
  seed_list = seed_list
)

configuration_path <- file.path(output_dir, "configuration.rds")
if (file.exists(configuration_path) && !overwrite) {
  cached_configuration <- readRDS(configuration_path)
  if (!isTRUE(all.equal(cached_configuration, configuration))) {
    stop("The existing output id has different settings.")
  }
} else {
  saveRDS(configuration, configuration_path)
}

source_output_path <- function(seed) {
  stem <- genotype_cosine_multipeak_output_stem(
    n_donors = n_donors,
    n_variants = J,
    time_grid = time_grid,
    n_covariates = n_covariates,
    expression_noise_sd = expression_noise_sd,
    width_half = width_half,
    target_centered_rms = source_centered_rms,
    shape_cell_probs = shape_cell_probs,
    class_probs = class_probs,
    seed = seed,
    scenario = source_scenario
  )
  file.path(
    workflowr_root,
    "output",
    "revision_simulations",
    "raw",
    paste0(stem, ".rds")
  )
}

load_source_output <- function(seed) {
  path <- source_output_path(seed)
  if (!file.exists(path)) {
    stop(
      "The paired current-scale fitted output is missing for seed ",
      seed,
      ": ",
      path
    )
  }
  out <- readRDS(path)
  required <- c(
    "unit_info", "true_beta_evaluation", "evaluation_grid", "fash_fits",
    "result_table", "alpha_curve", "eqtl_summary", "settings"
  )
  if (!all(required %in% names(out)) ||
      !isTRUE(all.equal(out$evaluation_grid, evaluation_grid)) ||
      !isTRUE(all.equal(out$settings$n_variants, J)) ||
      !isTRUE(all.equal(out$settings$n_donors, n_donors)) ||
      !isTRUE(all.equal(out$settings$n_covariates, n_covariates))) {
    stop("The current-scale fitted output failed validation for seed ", seed, ".")
  }
  out
}

load_current_mc_replicate <- function(seed) {
  path <- file.path(
    workflowr_root,
    "output",
    "revision_simulations",
    "mc",
    "sparse_timed_cosine_one_two_peak_pilot5",
    "replicates",
    paste0("seed_", seed, ".rds")
  )
  if (!file.exists(path)) {
    stop("The current-scale Monte Carlo replicate is missing: ", path)
  }
  replicate <- readRDS(path)
  required <- c("seed", "alpha_curve", "peak_alpha_curve")
  if (!all(required %in% names(replicate)) ||
      !identical(replicate$seed, seed) ||
      !all(global_methods %in% unique(replicate$alpha_curve$method)) ||
      !all(global_methods %in% unique(replicate$peak_alpha_curve$method))) {
    stop("The current-scale Monte Carlo replicate failed validation for seed ", seed, ".")
  }
  replicate
}

simulate_source_effects <- function(seed) {
  component_seeds <- revision_component_seeds(seed)
  simulate_raised_cosine_multipeak_effect_set(
    n_variants = J,
    time_grid = time_grid,
    evaluation_grid = evaluation_grid,
    class_probs = class_probs,
    width_levels = c(spiky = width_half),
    spike_counts = 1:2,
    shape_cell_probs = shape_cell_probs,
    primary_time_groups = primary_time_groups,
    center_by_observed_mean = FALSE,
    switch_threshold = current_switch_threshold,
    relative_amplitude_range = relative_amplitude_range,
    target_centered_rms = source_centered_rms,
    baseline_sd = 1,
    constant_sd = 1,
    dynamic_baseline_sd = 0,
    exact_class_counts = TRUE,
    seed = seed,
    class_seed = component_seeds[["classes"]],
    constant_seed = component_seeds[["constant_effects"]],
    shape_seed = component_seeds[["functional_truth"]],
    scenario = source_scenario
  )
}

rescale_dynamic_effects_to_max_abs <- function(effect_sim, target) {
  dynamic <- effect_sim$unit_info$effect_class == "dynamic_bspline"
  current_max_abs <- apply(
    abs(effect_sim$beta_evaluation[dynamic, , drop = FALSE]),
    1,
    max
  )
  if (any(!is.finite(current_max_abs)) ||
      any(current_max_abs <= .Machine$double.eps)) {
    stop("Cannot rescale a dynamic curve with zero or invalid maximum magnitude.")
  }
  scale_factor <- target / current_max_abs
  effect_sim$beta_matrix[dynamic, ] <-
    effect_sim$beta_matrix[dynamic, , drop = FALSE] * scale_factor
  effect_sim$beta_evaluation[dynamic, ] <-
    effect_sim$beta_evaluation[dynamic, , drop = FALSE] * scale_factor
  effect_sim$unit_info$baseline[dynamic] <-
    effect_sim$unit_info$baseline[dynamic] * scale_factor
  effect_sim$unit_info$centered_rms[dynamic] <- sqrt(rowMeans(
    (
      effect_sim$beta_matrix[dynamic, , drop = FALSE] -
        rowMeans(effect_sim$beta_matrix[dynamic, , drop = FALSE])
    )^2
  ))
  effect_sim$true_functionals <- evaluate_temporal_functionals(
    curves = effect_sim$beta_evaluation,
    smooth_var = evaluation_grid,
    switch_threshold = current_switch_threshold
  )
  effect_sim$unit_info$scenario <- maxabs_scenario

  rescaled_max_abs <- apply(
    abs(effect_sim$beta_evaluation[dynamic, , drop = FALSE]),
    1,
    max
  )
  if (any(abs(rescaled_max_abs - target) > 1e-10)) {
    stop("The max-absolute-effect normalization failed.")
  }

  geometry <- data.frame(
    variant_id = effect_sim$unit_info$variant_id[dynamic],
    cell_id = effect_sim$unit_info$cell_id[dynamic],
    spike_count = effect_sim$unit_info$spike_count[dynamic],
    sign_pattern = effect_sim$unit_info$sign_pattern[dynamic],
    current_max_abs = current_max_abs,
    scale_factor = scale_factor,
    rescaled_max_abs = rescaled_max_abs,
    rescaled_centered_rms = effect_sim$unit_info$centered_rms[dynamic],
    stringsAsFactors = FALSE
  )
  list(effect_sim = effect_sim, geometry = geometry)
}

build_maxabs_output <- function(seed) {
  fit_path <- file.path(fit_dir, paste0("seed_", seed, ".rds"))
  if (file.exists(fit_path) && !overwrite) {
    cached <- readRDS(fit_path)
    if (isTRUE(all.equal(cached$internal_configuration, configuration)) &&
        all(global_methods %in% unique(cached$result_table$method))) {
      message("Reusing max-absolute-effect fitted output: ", fit_path)
      return(cached)
    }
    stop("The cached max-absolute-effect fit does not match seed ", seed, ".")
  }

  message("Fitting max-absolute-effect condition for seed ", seed, ".")
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
  rescaled <- rescale_dynamic_effects_to_max_abs(
    simulate_source_effects(seed),
    target = target_max_abs
  )
  effect_sim <- rescaled$effect_sim
  dynamic <- effect_sim$unit_info$effect_class == "dynamic_bspline"
  zero_fraction <- rowMeans(
    abs(effect_sim$beta_evaluation[dynamic, , drop = FALSE]) < 1e-12
  )
  observed_cell_counts <- table(effect_sim$unit_info$cell_id[dynamic])
  if (sum(dynamic) != J * class_probs[["dynamic_bspline"]] ||
      !identical(
        as.integer(observed_cell_counts[names(shape_cell_probs)]),
        as.integer(J * class_probs[["dynamic_bspline"]] * shape_cell_probs)
      ) ||
      any(zero_fraction <= 0.55)) {
    stop("The rescaled sparse truth allocation failed validation.")
  }

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
    scenario = maxabs_scenario,
    output_dir = output_dir,
    save_outputs = FALSE,
    verbose = FALSE,
    effect_sim = effect_sim,
    expression_sim = expression_sim
  )
  out <- add_direct_interaction_efdr_results_to_genotype_output(
    out = out,
    n_permutations = efdr_permutations,
    alpha = 0.05,
    seed = component_seeds[["permutations"]],
    lambda = 0.5,
    pi0_method = "conservative",
    true_pi0 = true_pi0,
    include_true_pi0 = TRUE,
    permute_covariates_with_expression = TRUE,
    num_cores = num_cores,
    overwrite = TRUE,
    verbose = FALSE
  )
  if (!all(global_methods %in% unique(out$result_table$method))) {
    stop("The max-absolute-effect output is missing required methods.")
  }
  out$true_beta_evaluation <- effect_sim$beta_evaluation
  out$evaluation_grid <- evaluation_grid
  out$true_functionals <- effect_sim$true_functionals
  out$internal_geometry <- rescaled$geometry
  out$internal_configuration <- configuration
  out$component_seeds <- component_seeds
  saveRDS(out, fit_path)
  out
}

extract_global_rows <- function(out, seed, condition) {
  rows <- out$alpha_curve[out$alpha_curve$method %in% global_methods, , drop = FALSE]
  rows$scenario <- condition
  rows$seed <- seed
  rows
}

extract_peak_rows <- function(out, seed, condition) {
  result_table <- out$result_table[
    out$result_table$method %in% global_methods,
    ,
    drop = FALSE
  ]
  rows <- compute_dynamic_subgroup_alpha_curve(
    result_table = result_table,
    unit_info = out$unit_info,
    subgroup_var = "spike_count",
    alpha_grid = alpha_grid
  )
  rows$scenario <- condition
  rows$seed <- seed
  rows
}

evaluate_switch_curve <- function(out,
                                  fit,
                                  method,
                                  seed,
                                  condition,
                                  threshold,
                                  posterior_seed) {
  true_curves <- out$true_beta_evaluation[
    out$unit_info$variant_id,
    ,
    drop = FALSE
  ]
  true_switch <- evaluate_temporal_functionals(
    curves = true_curves,
    smooth_var = evaluation_grid,
    switch_threshold = threshold
  )[, "switch"]
  true_dynamic <- out$unit_info$effect_class == "dynamic_bspline"
  fdr_table <- get_fash_fdr_table(fit)
  maximum_dynamic_indices <- sort(unique(
    fdr_table$index[fdr_table$FDR <= max(alpha_grid)]
  ))
  switch_functional <- make_temporal_functionals(
    smooth_var = evaluation_grid,
    switch_threshold = threshold
  )["switch"]
  lfsr_map <- compute_functional_lfsr(
    fit = fit,
    functionals = switch_functional,
    indices = maximum_dynamic_indices,
    smooth_var = evaluation_grid,
    num_cores = num_cores,
    seed = posterior_seed
  )$switch

  alpha_rows <- lapply(alpha_grid, function(alpha) {
    dynamic_indices <- sort(unique(fdr_table$index[fdr_table$FDR <= alpha]))
    candidate_lfsr <- if (length(dynamic_indices) == 0) {
      numeric()
    } else {
      unname(lfsr_map[as.character(dynamic_indices)])
    }
    if (anyNA(candidate_lfsr)) {
      stop("Missing switch lfsr values for dynamically selected variants.")
    }
    cfsr_table <- functional_cfsr_table(dynamic_indices, candidate_lfsr)
    selected_indices <- cfsr_table$index[cfsr_table$cfsr <= alpha]
    true_null <- true_switch <= 0
    false_discoveries <- sum(true_null[selected_indices])
    true_positives <- sum(!true_null[selected_indices])
    conditional_indices <- selected_indices[true_dynamic[selected_indices]]
    conditional_false <- sum(true_null[conditional_indices])
    data.frame(
      scenario = condition,
      target = "switch",
      method = method,
      alpha = alpha,
      dynamic_discoveries = length(dynamic_indices),
      n_discoveries = length(selected_indices),
      false_discoveries = false_discoveries,
      conditional_discoveries = length(conditional_indices),
      conditional_false_discoveries = conditional_false,
      first_stage_null_calls = sum(!true_dynamic[selected_indices]),
      true_positives = true_positives,
      estimated_fsr = if (length(selected_indices) == 0) 0 else {
        mean(cfsr_table$lfsr[cfsr_table$index %in% selected_indices])
      },
      empirical_fsr = if (length(selected_indices) == 0) 0 else {
        false_discoveries / length(selected_indices)
      },
      conditional_empirical_fsr = if (length(conditional_indices) == 0) 0 else {
        conditional_false / length(conditional_indices)
      },
      power = if (sum(!true_null) == 0) NA_real_ else {
        true_positives / sum(!true_null)
      },
      seed = seed,
      stringsAsFactors = FALSE
    )
  })
  alpha_rows <- do.call(rbind, alpha_rows)

  alpha_005 <- 0.05
  dynamic_005 <- sort(unique(fdr_table$index[fdr_table$FDR <= alpha_005]))
  lfsr_005 <- if (length(dynamic_005) == 0) {
    numeric()
  } else {
    unname(lfsr_map[as.character(dynamic_005)])
  }
  cfsr_005 <- functional_cfsr_table(dynamic_005, lfsr_005)
  cfsr_005$functional_call <- cfsr_005$cfsr <= alpha_005
  cfsr_005$true_switch <- true_switch[cfsr_005$index] > 0
  cfsr_005$false_functional_call <-
    cfsr_005$functional_call & !cfsr_005$true_switch
  cfsr_005$effect_class <- out$unit_info$effect_class[cfsr_005$index]
  cfsr_005$cell_id <- out$unit_info$cell_id[cfsr_005$index]
  cfsr_005$spike_count <- out$unit_info$spike_count[cfsr_005$index]
  cfsr_005$sign_pattern <- out$unit_info$sign_pattern[cfsr_005$index]
  cfsr_005$method <- rep(method, nrow(cfsr_005))
  cfsr_005$condition <- rep(condition, nrow(cfsr_005))
  cfsr_005$seed <- rep(seed, nrow(cfsr_005))

  list(
    alpha_curve = alpha_rows,
    per_unit_alpha005 = cfsr_005,
    n_true_switch = sum(true_switch > 0)
  )
}

make_replicate <- function(seed) {
  current_out <- load_source_output(seed)
  current_mc <- load_current_mc_replicate(seed)
  maxabs_out <- build_maxabs_output(seed)
  component_seeds <- revision_component_seeds(seed)

  current_global <- current_mc$alpha_curve[
    current_mc$alpha_curve$method %in% global_methods,
    ,
    drop = FALSE
  ]
  current_global$scenario <- condition_labels[["current"]]
  current_global$seed <- seed
  global_alpha <- rbind(
    current_global,
    extract_global_rows(maxabs_out, seed, condition_labels[["maxabs"]])
  )
  current_peak <- current_mc$peak_alpha_curve[
    current_mc$peak_alpha_curve$method %in% global_methods,
    ,
    drop = FALSE
  ]
  current_peak$scenario <- condition_labels[["current"]]
  current_peak$seed <- seed
  peak_alpha <- rbind(
    current_peak,
    extract_peak_rows(maxabs_out, seed, condition_labels[["maxabs"]])
  )

  switch_results <- list()
  result_index <- 1L
  for (method_index in seq_along(functional_methods)) {
    method <- functional_methods[method_index]
    fit_name <- if (method == "FASH-IWP1-Raw") {
      "fash_iwp1_raw"
    } else {
      "fash_iwp1_bf"
    }
    switch_results[[result_index]] <- evaluate_switch_curve(
      out = current_out,
      fit = current_out$fash_fits[[fit_name]],
      method = method,
      seed = seed,
      condition = condition_labels[["current"]],
      threshold = current_switch_threshold,
      posterior_seed = component_seeds[["functional_posterior"]] +
        1000L * method_index
    )
    result_index <- result_index + 1L
    switch_results[[result_index]] <- evaluate_switch_curve(
      out = current_out,
      fit = current_out$fash_fits[[fit_name]],
      method = method,
      seed = seed,
      condition = condition_labels[["relative_threshold"]],
      threshold = relative_switch_threshold,
      posterior_seed = component_seeds[["functional_posterior"]] +
        10000L + 1000L * method_index
    )
    result_index <- result_index + 1L
    switch_results[[result_index]] <- evaluate_switch_curve(
      out = maxabs_out,
      fit = maxabs_out$fash_fits[[fit_name]],
      method = method,
      seed = seed,
      condition = condition_labels[["maxabs"]],
      threshold = current_switch_threshold,
      posterior_seed = component_seeds[["functional_posterior"]] +
        20000L + 1000L * method_index
    )
    result_index <- result_index + 1L
  }

  current_dynamic <- current_out$unit_info$effect_class == "dynamic_bspline"
  current_max_abs <- apply(
    abs(current_out$true_beta_evaluation[current_dynamic, , drop = FALSE]),
    1,
    max
  )
  maxabs_dynamic <- maxabs_out$unit_info$effect_class == "dynamic_bspline"
  maxabs_max_abs <- apply(
    abs(maxabs_out$true_beta_evaluation[maxabs_dynamic, , drop = FALSE]),
    1,
    max
  )
  geometry_summary <- data.frame(
    seed = seed,
    condition = c(condition_labels[["current"]], condition_labels[["maxabs"]]),
    max_abs_min = c(min(current_max_abs), min(maxabs_max_abs)),
    max_abs_median = c(stats::median(current_max_abs), stats::median(maxabs_max_abs)),
    max_abs_max = c(max(current_max_abs), max(maxabs_max_abs)),
    centered_rms_median = c(
      stats::median(sqrt(rowMeans(
        (
          current_out$true_beta[current_dynamic, , drop = FALSE] -
            rowMeans(current_out$true_beta[current_dynamic, , drop = FALSE])
        )^2
      ))),
      stats::median(sqrt(rowMeans(
        (
          maxabs_out$true_beta[maxabs_dynamic, , drop = FALSE] -
            rowMeans(maxabs_out$true_beta[maxabs_dynamic, , drop = FALSE])
        )^2
      )))
    ),
    stringsAsFactors = FALSE
  )

  list(
    configuration = configuration,
    seed = seed,
    global_alpha = global_alpha,
    peak_alpha = peak_alpha,
    switch_alpha = do.call(rbind, lapply(switch_results, `[[`, "alpha_curve")),
    switch_per_unit_alpha005 = do.call(
      rbind,
      lapply(switch_results, `[[`, "per_unit_alpha005")
    ),
    geometry_summary = geometry_summary,
    maxabs_unit_geometry = maxabs_out$internal_geometry
  )
}

validate_replicate <- function(x, seed) {
  required <- c(
    "configuration", "seed", "global_alpha", "peak_alpha", "switch_alpha",
    "switch_per_unit_alpha005", "geometry_summary", "maxabs_unit_geometry"
  )
  all(required %in% names(x)) &&
    identical(x$seed, seed) &&
    isTRUE(all.equal(x$configuration, configuration)) &&
    nrow(x$geometry_summary) == 2L &&
    all(global_methods %in% unique(x$global_alpha$method)) &&
    all(functional_methods %in% unique(x$switch_alpha$method)) &&
    all(condition_labels %in% unique(x$switch_alpha$scenario))
}

replicates <- lapply(seed_list, function(seed) {
  path <- file.path(replicate_dir, paste0("seed_", seed, ".rds"))
  if (file.exists(path) && !overwrite) {
    cached <- readRDS(path)
    if (validate_replicate(cached, seed)) {
      message("Reusing internal sensitivity replicate: ", path)
      return(cached)
    }
    stop("The cached sensitivity replicate does not match seed ", seed, ".")
  }
  message("Running paired internal sensitivity replicate with seed ", seed, ".")
  result <- make_replicate(seed)
  saveRDS(result, path)
  result
})

all_global <- do.call(rbind, lapply(replicates, `[[`, "global_alpha"))
all_peak <- do.call(rbind, lapply(replicates, `[[`, "peak_alpha"))
all_switch <- do.call(rbind, lapply(replicates, `[[`, "switch_alpha"))
all_switch_per_unit <- do.call(
  rbind,
  lapply(replicates, `[[`, "switch_per_unit_alpha005")
)
all_geometry <- do.call(rbind, lapply(replicates, `[[`, "geometry_summary"))
all_maxabs_unit_geometry <- do.call(rbind, lapply(replicates, function(x) {
  rows <- x$maxabs_unit_geometry
  rows$seed <- x$seed
  rows
}))

mc_global <- summarize_mc_alpha_curves(all_global)
mc_switch <- summarize_mc_functional_alpha_curves(all_switch)

summarize_peak_power <- function(rows) {
  groups <- split(
    rows,
    list(rows$scenario, rows$subgroup_value, rows$method, rows$alpha),
    drop = TRUE
  )
  out <- lapply(groups, function(x) {
    power <- summarize_mc_values(x$power)
    data.frame(
      scenario = x$scenario[1],
      spike_count = as.integer(x$subgroup_value[1]),
      method = x$method[1],
      alpha = x$alpha[1],
      n_replications = length(unique(x$seed)),
      mean_power = power[["mean"]],
      power_mc_se = power[["se"]],
      power_ci_lower = pmax(0, power[["lower"]]),
      power_ci_upper = pmin(1, power[["upper"]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

mc_peak <- summarize_peak_power(all_peak)
global_alpha005 <- mc_global[abs(mc_global$alpha - 0.05) < 1e-12, , drop = FALSE]
peak_alpha005 <- mc_peak[abs(mc_peak$alpha - 0.05) < 1e-12, , drop = FALSE]
switch_alpha005 <- mc_switch[abs(mc_switch$alpha - 0.05) < 1e-12, , drop = FALSE]

write_csv(all_global, file.path(summary_dir, "all_replicate_global_alpha.csv"))
write_csv(mc_global, file.path(summary_dir, "mc_global_alpha.csv"))
write_csv(global_alpha005, file.path(summary_dir, "mc_global_alpha005.csv"))
write_csv(all_peak, file.path(summary_dir, "all_replicate_peak_alpha.csv"))
write_csv(mc_peak, file.path(summary_dir, "mc_peak_alpha.csv"))
write_csv(peak_alpha005, file.path(summary_dir, "mc_peak_alpha005.csv"))
write_csv(all_switch, file.path(summary_dir, "all_replicate_switch_alpha.csv"))
write_csv(mc_switch, file.path(summary_dir, "mc_switch_alpha.csv"))
write_csv(switch_alpha005, file.path(summary_dir, "mc_switch_alpha005.csv"))
write_csv(
  all_switch_per_unit,
  file.path(summary_dir, "all_replicate_switch_per_unit_alpha005.csv")
)
write_csv(all_geometry, file.path(summary_dir, "geometry_summary.csv"))
write_csv(
  all_maxabs_unit_geometry,
  file.path(summary_dir, "maxabs_unit_geometry.csv")
)

plot_global_power <- function() {
  png(
    file.path(figure_dir, "global_power_by_effect_scale.png"),
    width = 2400,
    height = 1100,
    res = 180
  )
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.8, 1.0))
  styles <- revision_method_styles(global_methods, style_profile = "combined")
  for (condition in c(condition_labels[["current"]], condition_labels[["maxabs"]])) {
    rows <- mc_global[mc_global$scenario == condition, , drop = FALSE]
    plot(
      NA,
      xlim = range(alpha_grid),
      ylim = c(0, 1),
      xlab = "Nominal FDR level alpha",
      ylab = "Mean power",
      main = condition
    )
    grid(col = "gray90")
    abline(v = 0.05, col = "gray45", lty = 3)
    for (i in seq_along(styles$methods)) {
      method_rows <- rows[rows$method == styles$methods[i], , drop = FALSE]
      method_rows <- method_rows[order(method_rows$alpha), , drop = FALSE]
      lines(
        method_rows$alpha,
        method_rows$mean_power,
        col = styles$col[i],
        lty = styles$lty[i],
        lwd = 2.5
      )
    }
    if (condition == condition_labels[["current"]]) {
      legend(
        "bottomright",
        legend = styles$methods,
        col = styles$col,
        lty = styles$lty,
        lwd = 2.5,
        bty = "n",
        cex = 0.78
      )
    }
  }
}

plot_switch_sensitivity <- function() {
  png(
    file.path(figure_dir, "switch_power_and_empirical_fsr.png"),
    width = 2400,
    height = 1100,
    res = 180
  )
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.8, 1.0))
  colors <- c(
    current_scale_threshold_0.25 = "#D55E00",
    current_scale_threshold_0.75 = "#0072B2",
    max_abs_1_threshold_0.25 = "#009E73"
  )
  line_types <- c(`FASH-IWP1-Raw` = 1, `FASH-IWP1-BF` = 2)
  for (metric in c("mean_power", "mean_empirical_fsr")) {
    plot(
      NA,
      xlim = range(alpha_grid),
      ylim = c(0, 1),
      xlab = "Nominal FSR level alpha",
      ylab = if (metric == "mean_power") "Mean switch power" else {
        "Mean empirical switch FSR"
      },
      main = if (metric == "mean_power") {
        "Switch power"
      } else {
        "Switch empirical FSR"
      }
    )
    grid(col = "gray90")
    abline(v = 0.05, col = "gray45", lty = 3)
    if (metric == "mean_empirical_fsr") {
      abline(a = 0, b = 1, col = "gray35", lty = 3)
    }
    for (condition in names(colors)) {
      for (method in names(line_types)) {
        rows <- mc_switch[
          mc_switch$scenario == condition & mc_switch$method == method,
          ,
          drop = FALSE
        ]
        rows <- rows[order(rows$alpha), , drop = FALSE]
        lines(
          rows$alpha,
          rows[[metric]],
          col = colors[[condition]],
          lty = line_types[[method]],
          lwd = 2.5
        )
      }
    }
    if (metric == "mean_power") {
      legend(
        "bottomright",
        legend = c(
          "Current scale, c = 0.25",
          "Current scale, c = 0.75",
          "Max abs = 1, c = 0.25",
          "Raw",
          "BF-updated"
        ),
        col = c(unname(colors), "gray25", "gray25"),
        lty = c(1, 1, 1, 1, 2),
        lwd = 2.5,
        bty = "n",
        cex = 0.76
      )
    }
  }
}

plot_global_power()
plot_switch_sensitivity()

plot_false_switch_examples <- function(seed = seed_list[1], num_per_group = 3) {
  out <- readRDS(file.path(fit_dir, paste0("seed_", seed, ".rds")))
  rows <- all_switch_per_unit[
    all_switch_per_unit$seed == seed &
      all_switch_per_unit$condition == condition_labels[["maxabs"]] &
      all_switch_per_unit$method == "FASH-IWP1-BF" &
      all_switch_per_unit$false_functional_call,
    ,
    drop = FALSE
  ]
  rows$truth_group <- ifelse(
    rows$spike_count == 1,
    "one-peak",
    ifelse(rows$sign_pattern == "same-sign", "same-sign two-peak", "other")
  )
  rows <- rows[rows$truth_group %in% c("one-peak", "same-sign two-peak"), ]
  select_group <- function(group) {
    candidates <- rows[rows$truth_group == group, , drop = FALSE]
    candidates <- candidates[order(candidates$lfsr, candidates$index), ]
    head(candidates, num_per_group)
  }
  selected <- rbind(
    select_group("one-peak"),
    select_group("same-sign two-peak")
  )
  csv_path <- file.path(
    figure_dir,
    paste0("seed_", seed, "_bf_false_switch_examples.csv")
  )
  write_csv(selected, csv_path)
  if (nrow(selected) == 0) {
    message("No BF false switch discoveries were available for plotting.")
    return(invisible(NULL))
  }

  fit <- out$fash_fits$fash_iwp1_bf
  posterior <- lapply(selected$index, function(index) {
    set.seed(990000L + seed + index)
    samples <- predict(
      fit,
      index = index,
      smooth_var = evaluation_grid,
      only.samples = TRUE,
      M = posterior_samples
    )
    list(
      mean = rowMeans(samples),
      lower = apply(samples, 1, stats::quantile, probs = 0.025),
      upper = apply(samples, 1, stats::quantile, probs = 0.975)
    )
  })
  plot_path <- file.path(
    figure_dir,
    paste0("seed_", seed, "_bf_false_switch_examples.png")
  )
  panel_columns <- min(3L, nrow(selected))
  panel_rows <- ceiling(nrow(selected) / panel_columns)
  png(
    plot_path,
    width = 850 * panel_columns,
    height = 760 * panel_rows,
    res = 180
  )
  old_par <- par(no.readonly = TRUE)
  par(
    mfrow = c(panel_rows, panel_columns),
    mar = c(4.2, 4.3, 4.2, 1),
    oma = c(0, 0, 2.4, 0)
  )
  for (i in seq_len(nrow(selected))) {
    index <- selected$index[i]
    estimate <- out$eqtl_summary$beta_hat[index, ]
    standard_error <- out$eqtl_summary$se[index, ]
    truth <- out$true_beta_evaluation[index, ]
    posterior_i <- posterior[[i]]
    ylim <- range(
      truth,
      estimate - 1.96 * standard_error,
      estimate + 1.96 * standard_error,
      posterior_i$lower,
      posterior_i$upper
    )
    plot(
      evaluation_grid,
      posterior_i$mean,
      type = "n",
      ylim = ylim,
      xlab = "Time",
      ylab = "Genetic effect",
      main = paste0(
        selected$truth_group[i],
        "\nvariant ",
        index,
        "; lfsr = ",
        formatC(selected$lfsr[i], format = "f", digits = 3)
      )
    )
    polygon(
      c(evaluation_grid, rev(evaluation_grid)),
      c(posterior_i$lower, rev(posterior_i$upper)),
      border = NA,
      col = grDevices::adjustcolor("#4C78A8", alpha.f = 0.18)
    )
    abline(
      h = c(-current_switch_threshold, 0, current_switch_threshold),
      lty = c(3, 2, 3),
      col = "gray50"
    )
    lines(evaluation_grid, truth, col = "#E45756", lwd = 2.8)
    lines(evaluation_grid, posterior_i$mean, col = "#4C78A8", lwd = 2.4)
    arrows(
      time_grid,
      estimate - 1.96 * standard_error,
      time_grid,
      estimate + 1.96 * standard_error,
      angle = 90,
      code = 3,
      length = 0.025,
      col = "gray45"
    )
    points(time_grid, estimate, pch = 19, cex = 0.65)
  }
  if (nrow(selected) < panel_rows * panel_columns) {
    for (i in seq_len(panel_rows * panel_columns - nrow(selected))) plot.new()
  }
  mtext(
    "BF false switch discoveries after per-curve max-absolute-effect normalization",
    outer = TRUE,
    side = 3,
    line = 0.5,
    cex = 1.15,
    font = 2
  )
  par(old_par)
  dev.off()
  invisible(plot_path)
}

plot_false_switch_examples()

cat("\nGlobal results at alpha = 0.05:\n")
print(global_alpha005[
  global_alpha005$method %in% global_methods,
  c("scenario", "method", "mean_power", "mean_fdr")
])
cat("\nPeak-stratified results at alpha = 0.05:\n")
print(peak_alpha005[
  peak_alpha005$method %in% c(
    "FASH-IWP1-BF",
    "Direct-linear-LRT-eFDR-true-pi0",
    "Direct-quadratic-LRT-eFDR-true-pi0"
  ),
  c("scenario", "spike_count", "method", "mean_power")
])
cat("\nSwitch results at alpha = 0.05:\n")
print(switch_alpha005[
  ,
  c(
    "scenario", "method", "mean_discoveries", "mean_false_discoveries",
    "mean_power", "mean_empirical_fsr", "mean_conditional_empirical_fsr"
  )
])
cat("\nOutputs written to: ", normalizePath(output_dir), "\n", sep = "")
