#!/usr/bin/env Rscript

# Run compact five-seed replications for the reviewer-facing compact
# raised-cosine one-/two-/three-peak genotype simulation.

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

as_flag <- function(x) tolower(x) %in% c("1", "true", "t", "yes", "y")

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

summarize_peak_power <- function(rows, confidence_level = 0.95) {
  groups <- split(
    rows,
    list(rows$subgroup_value, rows$method, rows$alpha),
    drop = TRUE
  )
  out <- lapply(groups, function(x) {
    power <- summarize_mc_values(x$power, confidence_level)
    true_positives <- summarize_mc_values(x$true_positives, confidence_level)
    data.frame(
      spike_count = as.integer(x$subgroup_value[1]),
      shape_profile = paste0(x$subgroup_value[1], "-peak"),
      method = x$method[1],
      alpha = x$alpha[1],
      n_dynamic = x$n_dynamic[1],
      n_replications = length(unique(x$seed)),
      mean_true_positives = true_positives[["mean"]],
      mean_power = power[["mean"]],
      power_sd = power[["sd"]],
      power_mc_se = power[["se"]],
      power_ci_lower = pmax(0, power[["lower"]]),
      power_ci_upper = pmin(1, power[["upper"]]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

J <- as.integer(get_arg("--J", "1000"))
n_donors <- as.integer(get_arg("--n-donors", "19"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
expression_noise_sd <- as.numeric(get_arg("--noise-sd", "1"))
width_half <- as.numeric(get_arg("--width-half", "1.5"))
target_centered_rms <- as.numeric(get_arg("--target-centered-rms", "0.9"))
dynamic_main_effect_sd <- as.numeric(get_arg("--dynamic-main-effect-sd", "1"))
num_basis <- as.integer(get_arg("--num-basis", "20"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
efdr_permutations <- as.integer(get_arg("--efdr-permutations", "100"))
seed_list <- parse_seed_list(get_arg(
  "--seed-list",
  "12345,22345,32345,42345,52345"
))
center_by_observed_mean <- as_flag(get_arg(
  "--center-by-observed-mean",
  "false"
))
default_output_id <- if (center_by_observed_mean) {
  paste0(
    "r2_centered_timed_cosine_one_two_three_peak_",
    "main_effect_profile_sigma_pilot5"
  )
} else {
  paste0(
    "r2_timed_cosine_one_two_three_peak_",
    "main_effect_profile_sigma_pilot5"
  )
}
output_id <- get_arg("--output-id", default_output_id)
overwrite <- as_flag(get_arg("--overwrite", "false"))

if (J < 10 || n_donors < n_covariates + 3 || n_covariates < 0 ||
    !is.finite(expression_noise_sd) || expression_noise_sd <= 0 ||
    !is.finite(width_half) || width_half <= 0 ||
    !is.finite(target_centered_rms) || target_centered_rms <= 0 ||
    !is.finite(dynamic_main_effect_sd) || dynamic_main_effect_sd <= 0 ||
    num_basis < 2 || num_cores < 1 || efdr_permutations < 1 ||
    !nzchar(output_id)) {
  stop("Invalid cosine one-/two-/three-peak Monte Carlo arguments.")
}

time_grid <- make_time_grid()
evaluation_grid <- seq(min(time_grid), max(time_grid), by = 0.1)
class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
shape_cell_probs <- c(
  k1__spiky__single = 1 / 3,
  `k2__spiky__same-sign` = 1 / 6,
  `k2__spiky__alternating-sign` = 1 / 6,
  `k3__spiky__same-sign` = 1 / 6,
  `k3__spiky__alternating-sign` = 1 / 6
)
primary_time_groups <- c("early", "middle", "late")
relative_amplitude_range <- c(0.35, 0.75)
switch_threshold <- 0.25
scenario <- if (center_by_observed_mean) {
  "r2_genotype_centered_timed_cosine_one_two_three_peak_main_effect_dynamic_eqtl"
} else {
  "r2_genotype_timed_cosine_one_two_three_peak_main_effect_dynamic_eqtl"
}
true_pi0 <- unname(class_probs["constant"] + class_probs["zero"])
linear_sigma_estimation <- "profile_grid"
linear_sigma_grid <- exp(seq(log(0.05), log(5), length.out = 25))
fash_methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "FASH-linear-Raw",
  "FASH-linear-BF"
)
direct_methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "Direct-linear-LRT-eFDR-true-pi0",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
all_methods <- unique(c(fash_methods, direct_methods))

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
full_fit_dir <- file.path(output_dir, "full_fits")
invisible(lapply(
  c(output_dir, replicate_dir, summary_dir, figure_dir, full_fit_dir),
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
  dynamic_main_effect_sd = dynamic_main_effect_sd,
  num_basis = num_basis,
  class_probs = class_probs,
  spike_counts = 1:3,
  shape_cell_probs = shape_cell_probs,
  shape_cell_counts = exact_proportional_counts(
    exact_proportional_counts(J, class_probs)[["dynamic_bspline"]],
    shape_cell_probs
  ),
  primary_time_groups = primary_time_groups,
  relative_amplitude_range = relative_amplitude_range,
  center_by_observed_mean = center_by_observed_mean,
  switch_threshold = switch_threshold,
  constant_sd = 1,
  dynamic_baseline_sd = dynamic_main_effect_sd,
  efdr_permutations = efdr_permutations,
  true_pi0 = true_pi0,
  linear_sigma_estimation = linear_sigma_estimation,
  linear_sigma_grid = linear_sigma_grid,
  full_fit_seed = seed_list[1],
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

make_replicate <- function(seed) {
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
    spike_counts = 1:3,
    shape_cell_probs = shape_cell_probs,
    primary_time_groups = primary_time_groups,
    center_by_observed_mean = center_by_observed_mean,
    validate_functional_labels = FALSE,
    switch_threshold = switch_threshold,
    relative_amplitude_range = relative_amplitude_range,
    target_centered_rms = target_centered_rms,
    baseline_sd = 1,
    constant_sd = 1,
    dynamic_baseline_sd = dynamic_main_effect_sd,
    exact_class_counts = TRUE,
    seed = seed,
    class_seed = component_seeds[["classes"]],
    constant_seed = component_seeds[["constant_effects"]],
    shape_seed = component_seeds[["functional_truth"]],
    scenario = scenario
  )
  dynamic <- effect_sim$unit_info$effect_class == "dynamic_bspline"
  observed_cell_counts <- table(effect_sim$unit_info$cell_id[dynamic])
  timing_counts <- table(
    effect_sim$unit_info$cell_id[dynamic],
    effect_sim$unit_info$time_group[dynamic]
  )
  observed_mean <- rowMeans(
    effect_sim$beta_matrix[dynamic, , drop = FALSE]
  )
  centered_rms <- sqrt(rowMeans(
    (
      effect_sim$beta_matrix[dynamic, , drop = FALSE] -
        rowMeans(effect_sim$beta_matrix[dynamic, , drop = FALSE])
    )^2
  ))
  peak_lengths <- vapply(
    effect_sim$unit_info$peak_centers[dynamic],
    length,
    integer(1)
  )
  expected_class_counts <- exact_proportional_counts(J, class_probs)
  if (sum(dynamic) != expected_class_counts[["dynamic_bspline"]] ||
      !identical(
        as.integer(observed_cell_counts[names(shape_cell_probs)]),
        as.integer(configuration$shape_cell_counts)
      ) ||
      !setequal(colnames(timing_counts), primary_time_groups) ||
      any(apply(timing_counts, 1, function(x) diff(range(x))) > 1L) ||
      (
        center_by_observed_mean &&
          any(abs(
            observed_mean -
              effect_sim$unit_info$genetic_main_effect[dynamic]
          ) > 1e-10)
      ) ||
      any(peak_lengths != effect_sim$unit_info$spike_count[dynamic]) ||
      any(!is.finite(effect_sim$unit_info$genetic_main_effect[dynamic])) ||
      any(abs(centered_rms - target_centered_rms) > 1e-10)) {
    stop("Invalid timed cosine truth for seed ", seed, ".")
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
    estimate_sigma = TRUE,
    sigma_beta_grid = linear_sigma_grid,
    num_cores = num_cores,
    num_basis = num_basis,
    scenario = scenario,
    output_dir = output_dir,
    save_outputs = FALSE,
    verbose = FALSE,
    effect_sim = effect_sim,
    expression_sim = expression_sim
  )
  validate_simplified_sigma_profile(
    out$simplified_fit,
    require_interior = TRUE
  )
  validate_simplified_sigma_profile(
    out$simplified_fit_bf,
    require_interior = TRUE
  )
  if (!isTRUE(all.equal(
        out$simplified_fit$sigma_beta,
        out$simplified_fit_bf$sigma_beta,
        tolerance = 0
      ))) {
    stop("The BF update changed the selected linear slope scale.")
  }
  linear_sigma_profile <- out$simplified_fit$sigma_profile
  linear_sigma_profile$seed <- seed
  linear_sigma_profile <- linear_sigma_profile[, c(
    "seed", "sigma_beta", "estimated_pi0", "loglik", "selected",
    "grid_boundary"
  )]
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
  out$true_beta_evaluation <- effect_sim$beta_evaluation
  out$evaluation_grid <- evaluation_grid
  out$true_functionals <- effect_sim$true_functionals
  if (identical(seed, seed_list[1])) {
    saveRDS(
      out,
      file.path(full_fit_dir, paste0("seed_", seed, ".rds"))
    )
  }
  missing_methods <- setdiff(all_methods, unique(out$result_table$method))
  if (length(missing_methods) > 0) {
    stop("Missing methods for seed ", seed, ": ", paste(missing_methods, collapse = ", "))
  }
  result_table <- out$result_table[out$result_table$method %in% all_methods, ]
  alpha_curve <- out$alpha_curve[out$alpha_curve$method %in% all_methods, ]
  alpha_curve$seed <- seed
  peak_alpha_curve <- compute_dynamic_subgroup_alpha_curve(
    result_table = result_table,
    unit_info = out$unit_info,
    subgroup_var = "spike_count",
    alpha_grid = seq(0, 0.20, by = 0.005)
  )
  peak_alpha_curve$seed <- seed
  pi0 <- data.frame(
    seed = seed,
    method = c("FASH-IWP1", "FASH-IWP1", "FASH-linear", "FASH-linear"),
    fit = c("Raw", "BF-corrected", "Raw", "BF-corrected"),
    estimated_pi0 = c(
      constant_component_prior_weight(out$fash_fits$fash_iwp1_raw),
      constant_component_prior_weight(out$fash_fits$fash_iwp1_bf),
      constant_component_prior_weight(out$simplified_fit),
      constant_component_prior_weight(out$simplified_fit_bf)
    ),
    stringsAsFactors = FALSE
  )
  geometry <- cbind(
    seed = seed,
    summarize_dynamic_effect_geometry(
      effect_sim,
      dynamic_class = "dynamic_bspline"
    )
  )
  list(
    configuration = configuration,
    seed = seed,
    component_seeds = component_seeds,
    alpha_curve = alpha_curve,
    alpha_005 = alpha_curve[abs(alpha_curve$alpha - 0.05) < 1e-12, ],
    peak_alpha_curve = peak_alpha_curve,
    peak_alpha_005 = peak_alpha_curve[
      abs(peak_alpha_curve$alpha - 0.05) < 1e-12,
    ],
    pi0 = pi0,
    geometry = geometry,
    linear_sigma_profile = linear_sigma_profile
  )
}

validate_replicate <- function(x, seed) {
  required <- c(
    "configuration", "seed", "component_seeds", "alpha_curve", "alpha_005",
    "peak_alpha_curve", "peak_alpha_005", "pi0", "geometry",
    "linear_sigma_profile"
  )
  profile <- x$linear_sigma_profile
  all(required %in% names(x)) &&
    identical(x$seed, seed) &&
    isTRUE(all.equal(x$configuration, configuration)) &&
    all(all_methods %in% unique(x$alpha_curve$method)) &&
    all(all_methods %in% unique(x$peak_alpha_curve$method)) &&
    identical(sort(unique(x$peak_alpha_curve$subgroup_value)), c("1", "2", "3")) &&
    is.data.frame(profile) &&
    nrow(profile) == length(linear_sigma_grid) &&
    identical(profile$seed, rep(seed, length(linear_sigma_grid))) &&
    isTRUE(all.equal(profile$sigma_beta, linear_sigma_grid, tolerance = 0)) &&
    all(is.finite(profile$estimated_pi0)) &&
    all(is.finite(profile$loglik)) &&
    sum(profile$selected) == 1L &&
    !any(profile$selected & profile$grid_boundary)
}

replicates <- lapply(seed_list, function(seed) {
  path <- file.path(replicate_dir, paste0("seed_", seed, ".rds"))
  if (file.exists(path) && !overwrite) {
    cached <- readRDS(path)
    if (validate_replicate(cached, seed)) {
      message("Reusing cosine one-/two-/three-peak replicate: ", path)
      return(cached)
    }
    stop("Cached replicate does not match: ", path)
  }
  message("Running cosine one-/two-/three-peak replicate with seed ", seed, ".")
  replicate <- make_replicate(seed)
  saveRDS(replicate, path)
  replicate
})

all_alpha <- do.call(rbind, lapply(replicates, `[[`, "alpha_curve"))
all_alpha_005 <- do.call(rbind, lapply(replicates, `[[`, "alpha_005"))
all_peak_alpha <- do.call(rbind, lapply(replicates, `[[`, "peak_alpha_curve"))
all_peak_alpha_005 <- do.call(rbind, lapply(replicates, `[[`, "peak_alpha_005"))
all_pi0 <- do.call(rbind, lapply(replicates, `[[`, "pi0"))
all_geometry <- do.call(rbind, lapply(replicates, `[[`, "geometry"))
all_linear_sigma_profiles <- do.call(
  rbind,
  lapply(replicates, `[[`, "linear_sigma_profile")
)
selected_linear_sigma <- all_linear_sigma_profiles[
  all_linear_sigma_profiles$selected,
  ,
  drop = FALSE
]
if (nrow(selected_linear_sigma) != length(seed_list) ||
    any(selected_linear_sigma$grid_boundary) ||
    !identical(sort(selected_linear_sigma$seed), sort(seed_list))) {
  stop("The selected linear slope scales are incomplete or on a grid boundary.")
}
linear_sigma_values <- summarize_mc_values(selected_linear_sigma$sigma_beta)
linear_sigma_summary <- data.frame(
  estimation = linear_sigma_estimation,
  n_replications = length(seed_list),
  mean_selected_sigma = linear_sigma_values[["mean"]],
  selected_sigma_sd = linear_sigma_values[["sd"]],
  selected_sigma_mc_se = linear_sigma_values[["se"]],
  selected_sigma_ci_lower = pmax(0, linear_sigma_values[["lower"]]),
  selected_sigma_ci_upper = linear_sigma_values[["upper"]],
  min_selected_sigma = min(selected_linear_sigma$sigma_beta),
  max_selected_sigma = max(selected_linear_sigma$sigma_beta),
  stringsAsFactors = FALSE
)

mc_alpha <- summarize_mc_alpha_curves(all_alpha)
mc_alpha_005 <- mc_alpha[abs(mc_alpha$alpha - 0.05) < 1e-12, ]
mc_peak_alpha <- summarize_peak_power(all_peak_alpha)
mc_peak_alpha_005 <- mc_peak_alpha[abs(mc_peak_alpha$alpha - 0.05) < 1e-12, ]
mc_pi0 <- summarize_mc_pi0(all_pi0)

fash_all_alpha <- all_alpha[all_alpha$method %in% fash_methods, ]
direct_all_alpha <- all_alpha[all_alpha$method %in% direct_methods, ]
fash_mc_alpha <- mc_alpha[mc_alpha$method %in% fash_methods, ]
direct_mc_alpha <- mc_alpha[mc_alpha$method %in% direct_methods, ]
fash_mc_alpha_005 <- mc_alpha_005[mc_alpha_005$method %in% fash_methods, ]
direct_mc_alpha_005 <- mc_alpha_005[mc_alpha_005$method %in% direct_methods, ]

write_csv(all_alpha, file.path(summary_dir, "all_replicate_alpha_curves.csv"))
write_csv(all_alpha_005, file.path(summary_dir, "all_replicate_alpha005.csv"))
write_csv(all_peak_alpha, file.path(summary_dir, "all_replicate_peak_alpha_curves.csv"))
write_csv(all_peak_alpha_005, file.path(summary_dir, "all_replicate_peak_alpha005.csv"))
write_csv(all_pi0, file.path(summary_dir, "all_replicate_pi0.csv"))
write_csv(all_geometry, file.path(summary_dir, "all_replicate_geometry.csv"))
write_csv(
  all_linear_sigma_profiles,
  file.path(summary_dir, "all_replicate_linear_sigma_profiles.csv")
)
write_csv(
  linear_sigma_summary,
  file.path(summary_dir, "linear_sigma_summary.csv")
)
write_csv(mc_alpha, file.path(summary_dir, "mc_alpha_curve.csv"))
write_csv(mc_alpha_005, file.path(summary_dir, "mc_alpha005_summary.csv"))
write_csv(mc_peak_alpha, file.path(summary_dir, "mc_peak_alpha_curve.csv"))
write_csv(mc_peak_alpha_005, file.path(summary_dir, "mc_peak_alpha005_summary.csv"))
write_csv(mc_pi0, file.path(summary_dir, "mc_pi0_summary.csv"))
write_csv(fash_all_alpha, file.path(summary_dir, "iwp_vs_linear_fash_replicate_alpha_curves.csv"))
write_csv(fash_mc_alpha, file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha_curve.csv"))
write_csv(fash_mc_alpha_005, file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha005_summary.csv"))
write_csv(direct_all_alpha, file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_replicate_alpha_curves.csv"))
write_csv(direct_mc_alpha, file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_mc_alpha_curve.csv"))
write_csv(direct_mc_alpha_005, file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_mc_alpha005_summary.csv"))

subtitle <- paste0(
  length(seed_list),
  " seeds; ",
  if (center_by_observed_mean) "centered" else "sparse",
  " timed raised-cosine one-, two-, and three-peak effects; N = ",
  n_donors,
  ", J = ",
  J
)
plot_mc_alpha_curves(
  fash_mc_alpha,
  metric = "power",
  file = file.path(figure_dir, "iwp_vs_linear_fash_mc_power.png"),
  title = "IWP versus linear FASH: Monte Carlo power",
  subtitle = subtitle,
  style_profile = "combined"
)
plot_mc_alpha_curves(
  fash_mc_alpha,
  metric = "fdr",
  file = file.path(figure_dir, "iwp_vs_linear_fash_mc_fdr.png"),
  title = "IWP versus linear FASH: Monte Carlo FDR",
  subtitle = subtitle,
  legend_position = "topleft",
  style_profile = "combined"
)
plot_mc_alpha_curves(
  direct_mc_alpha,
  metric = "power",
  file = file.path(figure_dir, "iwp_fash_vs_direct_true_pi0_mc_power.png"),
  title = "IWP FASH versus direct interactions: Monte Carlo power",
  subtitle = paste0(subtitle, "; direct eFDR uses true pi0 = 0.8"),
  style_profile = "combined"
)
plot_mc_alpha_curves(
  direct_mc_alpha,
  metric = "fdr",
  file = file.path(figure_dir, "iwp_fash_vs_direct_true_pi0_mc_fdr.png"),
  title = "IWP FASH versus direct interactions: Monte Carlo FDR",
  subtitle = paste0(subtitle, "; direct eFDR uses true pi0 = 0.8"),
  legend_position = "bottomleft",
  style_profile = "combined"
)
peak_plot <- mc_peak_alpha[
  mc_peak_alpha$method %in% c(
    "FASH-IWP1-BF",
    "Direct-linear-LRT-eFDR-true-pi0",
    "Direct-quadratic-LRT-eFDR-true-pi0"
  ),
]
plot_mc_shape_power_grid(
  peak_plot,
  shape_order = c("1-peak", "2-peak", "3-peak"),
  file = file.path(figure_dir, "peak_count_power_mc.png"),
  title = "Power stratified by the number of transient peaks",
  subtitle = "Global method thresholds are applied before subgroup power is calculated",
  style_profile = "combined",
  width = 2600
)

print(fash_mc_alpha_005)
print(direct_mc_alpha_005)
print(mc_peak_alpha_005)
print(mc_pi0)
print(linear_sigma_summary)
