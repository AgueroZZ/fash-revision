#!/usr/bin/env Rscript

# Run the reviewer-facing genotype-level simulation with compact
# raised-cosine one-, two-, and three-peak dynamic effects.

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

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

J <- as.integer(get_arg("--J", "1000"))
n_donors <- as.integer(get_arg("--n-donors", "19"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
expression_noise_sd <- as.numeric(get_arg("--noise-sd", "1"))
width_half <- as.numeric(get_arg("--width-half", "1.5"))
target_centered_rms <- as.numeric(get_arg("--target-centered-rms", "0.9"))
dynamic_main_effect_sd <- as.numeric(get_arg("--dynamic-main-effect-sd", "1"))
seed <- as.integer(get_arg("--seed", "12345"))
num_basis <- as.integer(get_arg("--num-basis", "20"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
efdr_permutations <- as.integer(get_arg("--efdr-permutations", "100"))
overwrite <- as_flag(get_arg("--overwrite", "false"))

if (J < 10 || n_donors < n_covariates + 3 || n_covariates < 0 ||
    !is.finite(expression_noise_sd) || expression_noise_sd <= 0 ||
    !is.finite(width_half) || width_half <= 0 ||
    !is.finite(target_centered_rms) || target_centered_rms <= 0 ||
    !is.finite(dynamic_main_effect_sd) || dynamic_main_effect_sd <= 0 ||
    !is.finite(seed) || num_basis < 2 || num_cores < 1 ||
    efdr_permutations < 1) {
  stop("Invalid cosine one-/two-/three-peak simulation arguments.")
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
center_by_observed_mean <- FALSE
switch_threshold <- 0.25
scenario <- "r2_genotype_timed_cosine_one_two_three_peak_main_effect_dynamic_eqtl"
component_seeds <- revision_component_seeds(seed)
true_pi0 <- unname(class_probs["constant"] + class_probs["zero"])
linear_prior_mode <- "mixture_grid"
common_sd_grid <- default_revision_grid()
common_pred_step <- 1
common_penalty <- 10L
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
stem <- paste0(stem, "_linear_mixture_predstep1_penalty10")
raw_dir <- file.path(workflowr_root, "output", "revision_simulations", "raw")
summary_dir <- file.path(workflowr_root, "output", "revision_simulations", "summary")
figure_dir <- file.path(workflowr_root, "output", "revision_simulations", "figures")
invisible(lapply(
  c(raw_dir, summary_dir, figure_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))
raw_path <- file.path(raw_dir, paste0(stem, ".rds"))

configuration <- list(
  stem = stem,
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
  class_probs = class_probs,
  spike_counts = 1:3,
  shape_cell_probs = shape_cell_probs,
  primary_time_groups = primary_time_groups,
  relative_amplitude_range = relative_amplitude_range,
  center_by_observed_mean = center_by_observed_mean,
  switch_threshold = switch_threshold,
  constant_sd = 1,
  dynamic_baseline_sd = dynamic_main_effect_sd,
  seed = seed,
  num_basis = num_basis,
  efdr_permutations = efdr_permutations,
  true_pi0 = true_pi0,
  linear_prior_mode = linear_prior_mode,
  common_sd_grid = common_sd_grid,
  common_pred_step = common_pred_step,
  common_penalty = common_penalty
)

if (file.exists(raw_path) && !overwrite) {
  message("Reusing existing cosine one-/two-/three-peak output: ", raw_path)
  out <- readRDS(raw_path)
  if (!isTRUE(all.equal(out$cosine_configuration, configuration))) {
    stop("The cached output does not match the requested cosine configuration.")
  }
} else {
  message("Simulating cosine one-/two-/three-peak genotype data: ", stem)
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
  constant <- effect_sim$unit_info$effect_class == "constant"
  zero <- effect_sim$unit_info$effect_class == "zero"
  observed_cell_counts <- table(effect_sim$unit_info$cell_id[dynamic])
  expected_cell_counts <- exact_proportional_counts(
    sum(dynamic),
    shape_cell_probs
  )
  timing_counts <- table(
    effect_sim$unit_info$cell_id[dynamic],
    effect_sim$unit_info$time_group[dynamic]
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
  if (sum(dynamic) != exact_proportional_counts(J, class_probs)["dynamic_bspline"] ||
      sum(constant) != exact_proportional_counts(J, class_probs)["constant"] ||
      sum(zero) != exact_proportional_counts(J, class_probs)["zero"] ||
      !identical(
        as.integer(observed_cell_counts[names(shape_cell_probs)]),
        as.integer(expected_cell_counts)
      ) ||
      !setequal(colnames(timing_counts), primary_time_groups) ||
      any(apply(timing_counts, 1, function(x) diff(range(x))) > 1L) ||
      any(peak_lengths != effect_sim$unit_info$spike_count[dynamic]) ||
      any(!is.finite(effect_sim$unit_info$genetic_main_effect[dynamic])) ||
      any(abs(centered_rms - target_centered_rms) > 1e-10)) {
    stop("The sparse timed cosine truth allocation is invalid.")
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
    grid = common_sd_grid,
    penalty = common_penalty,
    pred_step = common_pred_step,
    linear_prior_mode = linear_prior_mode,
    num_cores = num_cores,
    num_basis = num_basis,
    scenario = scenario,
    output_dir = raw_dir,
    save_outputs = FALSE,
    verbose = FALSE,
    effect_sim = effect_sim,
    expression_sim = expression_sim
  )
  validate_linear_mixture_fash(
    out$simplified_fit,
    expected_grid = common_sd_grid,
    expected_pred_step = common_pred_step,
    expected_penalty = common_penalty
  )
  validate_linear_mixture_fash(
    out$simplified_fit_bf,
    expected_grid = common_sd_grid,
    expected_pred_step = common_pred_step,
    expected_penalty = common_penalty
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
  out$true_beta_evaluation <- effect_sim$beta_evaluation
  out$evaluation_grid <- evaluation_grid
  out$true_functionals <- effect_sim$true_functionals
  out$cosine_configuration <- configuration
  out$component_seeds <- component_seeds
  saveRDS(out, raw_path)
}

missing_methods <- setdiff(all_methods, unique(out$result_table$method))
if (length(missing_methods) > 0) {
  stop("Missing reviewer-facing methods: ", paste(missing_methods, collapse = ", "))
}
validate_linear_mixture_fash(
  out$simplified_fit,
  expected_grid = common_sd_grid,
  expected_pred_step = common_pred_step,
  expected_penalty = common_penalty
)
validate_linear_mixture_fash(
  out$simplified_fit_bf,
  expected_grid = common_sd_grid,
  expected_pred_step = common_pred_step,
  expected_penalty = common_penalty
)
for (fit_name in c("fash_iwp1_raw", "fash_iwp1_bf")) {
  iwp_fit <- out$fash_fits[[fit_name]]
  if (!isTRUE(all.equal(iwp_fit$psd_grid, common_sd_grid, tolerance = 0)) ||
      !isTRUE(all.equal(
        iwp_fit$settings$pred_step,
        common_pred_step,
        tolerance = 0
      )) ||
      !identical(
        as.integer(iwp_fit$settings$penalty),
        as.integer(common_penalty)
      )) {
    stop("The IWP and linear fits do not share grid, pred_step, and penalty.")
  }
}
selected_results <- out$result_table[out$result_table$method %in% all_methods, ]
alpha_curve <- out$alpha_curve[out$alpha_curve$method %in% all_methods, ]
alpha_005 <- alpha_curve[abs(alpha_curve$alpha - 0.05) < 1e-12, ]
peak_alpha_curve <- compute_dynamic_subgroup_alpha_curve(
  result_table = selected_results,
  unit_info = out$unit_info,
  subgroup_var = "spike_count",
  alpha_grid = seq(0, 0.20, by = 0.005)
)
peak_alpha_005 <- peak_alpha_curve[abs(peak_alpha_curve$alpha - 0.05) < 1e-12, ]
pi0_table <- data.frame(
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
geometry <- summarize_dynamic_effect_geometry(
  list(beta_matrix = out$true_beta, unit_info = out$unit_info),
  dynamic_class = "dynamic_bspline"
)
linear_prior_weights <- rbind(
  extract_linear_mixture_prior_table(
    out$simplified_fit,
    seed = seed,
    fit_label = "Raw"
  ),
  extract_linear_mixture_prior_table(
    out$simplified_fit_bf,
    seed = seed,
    fit_label = "BF-corrected"
  )
)
linear_prior_summary <- rbind(
  summarize_linear_mixture_prior_fit(
    out$simplified_fit,
    seed = seed,
    fit_label = "Raw"
  ),
  summarize_linear_mixture_prior_fit(
    out$simplified_fit_bf,
    seed = seed,
    fit_label = "BF-corrected"
  )
)

saveRDS(configuration, file.path(summary_dir, paste0(stem, "_configuration.rds")))
write.csv(selected_results, file.path(summary_dir, paste0(stem, "_unit_results.csv")),
          row.names = FALSE)
write.csv(alpha_curve, file.path(summary_dir, paste0(stem, "_alpha_curve.csv")),
          row.names = FALSE)
write.csv(alpha_005, file.path(summary_dir, paste0(stem, "_alpha005.csv")),
          row.names = FALSE)
write.csv(peak_alpha_curve, file.path(summary_dir, paste0(stem, "_peak_alpha_curve.csv")),
          row.names = FALSE)
write.csv(peak_alpha_005, file.path(summary_dir, paste0(stem, "_peak_alpha005.csv")),
          row.names = FALSE)
write.csv(pi0_table, file.path(summary_dir, paste0(stem, "_pi0.csv")),
          row.names = FALSE)
write.csv(
  linear_prior_weights,
  file.path(summary_dir, paste0(stem, "_linear_prior_weights.csv")),
  row.names = FALSE
)
write.csv(
  linear_prior_summary,
  file.path(summary_dir, paste0(stem, "_linear_prior_summary.csv")),
  row.names = FALSE
)
write.csv(geometry, file.path(summary_dir, paste0(stem, "_geometry.csv")),
          row.names = FALSE)

plot_genotype_eqtl_examples(
  out,
  classes = c("zero", "constant"),
  n_per_class = 3,
  file = file.path(figure_dir, paste0(stem, "_null_examples.png")),
  seed = 1
)
plot_cosine_peak_cell_examples(
  out,
  file = file.path(figure_dir, paste0(stem, "_dynamic_cell_examples.png")),
  seed = 1,
  include_time_groups = TRUE
)
plot_power_alpha_curves(
  alpha_curve[alpha_curve$method %in% fash_methods, ],
  file = file.path(figure_dir, paste0(stem, "_fash_power.png")),
  title = "IWP versus linear FASH: single-replicate power",
  subtitle = "Timed raised-cosine one-, two-, and three-peak effects",
  style_profile = "combined"
)
plot_power_alpha_curves(
  alpha_curve[alpha_curve$method %in% direct_methods, ],
  file = file.path(figure_dir, paste0(stem, "_direct_power.png")),
  title = "IWP FASH versus direct interactions: single-replicate power",
  subtitle = "Direct eFDR uses true pi0 = 0.8 and B = 100",
  style_profile = "combined"
)

message("Saved reviewer-facing sparse timed cosine output: ", raw_path)
