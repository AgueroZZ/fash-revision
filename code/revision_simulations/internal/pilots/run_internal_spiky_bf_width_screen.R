#!/usr/bin/env Rscript

# Diagnose and screen the temporal support width of localized multi-spike
# effects while keeping the broad effects and total centered signal fixed.

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

parse_integer_list <- function(x) {
  out <- suppressWarnings(as.integer(trimws(strsplit(x, ",", fixed = TRUE)[[1]])))
  if (length(out) == 0 || anyNA(out) || anyDuplicated(out)) {
    stop("Integer list arguments must contain unique comma-separated integers.")
  }
  out
}

parse_logical <- function(x, name) {
  value <- tolower(trimws(x))
  if (value %in% c("true", "t", "1", "yes", "y")) return(TRUE)
  if (value %in% c("false", "f", "0", "no", "n")) return(FALSE)
  stop(name, " must be true or false.")
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

cumulative_lfdr_score <- function(lfdr) {
  ordering <- order(lfdr)
  cumulative <- cumsum(lfdr[ordering]) / seq_along(lfdr)
  cumulative[match(seq_along(lfdr), ordering)]
}

truth_geometry <- function(beta_matrix, time_grid) {
  beta_matrix <- as.matrix(beta_matrix)
  quadratic_basis <- scaled_time_polynomial(time_grid, degree = 2)
  quadratic_qr <- qr(quadratic_basis)
  centered <- t(apply(beta_matrix, 1, function(x) x - mean(x)))
  maximum <- apply(abs(centered), 1, max)
  energy <- rowSums(centered^2)
  data.frame(
    centered_rms = sqrt(rowMeans(centered^2)),
    maximum_absolute_centered_effect = maximum,
    support_25 = rowSums(abs(centered) >= 0.25 * maximum),
    support_50 = rowSums(abs(centered) >= 0.50 * maximum),
    first_difference_roughness = apply(
      beta_matrix,
      1,
      function(x) sum(diff(x)^2) / sum((x - mean(x))^2)
    ),
    second_difference_roughness = apply(
      beta_matrix,
      1,
      function(x) sum(diff(x, differences = 2)^2) / sum((x - mean(x))^2)
    ),
    quadratic_projection = apply(centered, 1, function(x) {
      if (sum(x^2) <= .Machine$double.eps) return(0)
      fitted <- qr.fitted(quadratic_qr, x)
      sum(fitted^2) / sum(x^2)
    }),
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

J <- as.integer(get_arg("--J", "1000"))
n_donors <- as.integer(get_arg("--n-donors", "19"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
num_basis <- as.integer(get_arg("--num-basis", "20"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
efdr_permutations <- as.integer(get_arg("--efdr-permutations", "30"))
seed <- as.integer(get_arg("--seed", "12345"))
spiky_dfs <- parse_integer_list(get_arg("--spiky-dfs", "12,14,16"))
base_dgp <- get_arg("--base-dgp", "internal_broad075")
secondary_selection <- get_arg("--secondary-selection", "random")
secondary_fraction <- c(
  as.numeric(get_arg("--secondary-fraction-min", "0.40")),
  as.numeric(get_arg("--secondary-fraction-max", "0.65"))
)
pair_truth_rng_across_widths <- parse_logical(
  get_arg("--pair-truth-rng-across-widths", "false"),
  "--pair-truth-rng-across-widths"
)
output_id <- get_arg("--output-id", paste0("spiky_bf_width_seed", seed, "_B", efdr_permutations))

if (J != 1000 || n_donors != 19 || n_covariates != 5 ||
    num_basis < 2 || num_cores < 1 || efdr_permutations < 1 ||
    any(spiky_dfs <= 3) ||
    !base_dgp %in% c("internal_broad075", "reviewer_mixed") ||
    !secondary_selection %in% c("random", "nearest") ||
    length(secondary_fraction) != 2 ||
    any(!is.finite(secondary_fraction)) ||
    any(secondary_fraction <= 0) ||
    secondary_fraction[1] > secondary_fraction[2] ||
    secondary_fraction[2] >= 1 ||
    !nzchar(output_id)) {
  stop("Invalid spiky-width screen arguments.")
}

time_grid <- make_time_grid()
evaluation_grid <- seq(min(time_grid), max(time_grid), by = 0.1)
class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
true_pi0 <- 0.8
component_seeds <- revision_component_seeds(seed)
base_effects <- NULL
if (base_dgp == "internal_broad075") {
  base_candidate_path <- file.path(
    workflowr_root,
    "output",
    "revision_simulations",
    "internal",
    paste0("confirm_broad075_df10_seed", seed, "_B100"),
    "candidates",
    "broad075_df10.rds"
  )
  if (!file.exists(base_candidate_path)) {
    stop("Missing locked broad/spiky candidate cache: ", base_candidate_path)
  }
  base_candidate <- readRDS(base_candidate_path)
  base_effects <- base_candidate$effect_sim
}

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  output_id
)
candidate_dir <- file.path(output_dir, "candidates")
dir.create(candidate_dir, recursive = TRUE, showWarnings = FALSE)

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

configuration <- list(
  output_id = output_id,
  internal_only = TRUE,
  seed = seed,
  spiky_dfs = spiky_dfs,
  J = J,
  n_donors = n_donors,
  n_covariates = n_covariates,
  num_basis = num_basis,
  efdr_permutations = efdr_permutations,
  base_dgp = base_dgp,
  secondary_selection = secondary_selection,
  secondary_fraction = secondary_fraction,
  pair_truth_rng_across_widths = pair_truth_rng_across_widths,
  true_pi0 = true_pi0,
  spiky_centered_rms = 0.9,
  broad_component = if (base_dgp == "reviewer_mixed") {
    "50% constrained random cubic B-spline, df=6"
  } else {
    "75% random cubic B-spline, df=10, amplitude=2"
  },
  paired_genotype_covariates_noise_and_permutations = TRUE
)
saveRDS(configuration, file.path(output_dir, "configuration.rds"))

candidate_results <- vector("list", length(spiky_dfs))
for (candidate_index in seq_along(spiky_dfs)) {
  spiky_df <- spiky_dfs[candidate_index]
  if (base_dgp == "reviewer_mixed") {
    effect_sim <- simulate_targeted_local_bspline_effect_set(
      n_variants = J,
      time_grid = time_grid,
      evaluation_grid = evaluation_grid,
      class_probs = class_probs,
      dynamic_amplitude = 2,
      switch_threshold = 0.25,
      minimum_location_margin = 0.60,
      minimum_location_ratio = 1.4,
      non_switch_baseline_fraction = 0.75,
      non_switch_background_fraction = 0.05,
      profile = "mixed",
      spiky_truth_version = "mixed_single_double_v2",
      spiky_secondary_fraction = secondary_fraction,
      spiky_bspline_df = spiky_df,
      spiky_secondary_selection = secondary_selection,
      spiky_minimum_peak_separation = 3,
      spiky_non_switch_baseline_fraction = 0.20,
      target_centered_rms = 0.90,
      exact_class_counts = TRUE,
      seed = component_seeds[["functional_truth"]],
      scenario = "internal_reviewer_mixed_spiky_width_screen"
    )
  } else {
    effect_sim <- base_effects
    spiky_index <- which(
      effect_sim$unit_info$effect_class == "dynamic_bspline" &
        effect_sim$unit_info$shape_profile == "spiky"
    )
    if (length(spiky_index) != 50) {
      stop("The locked candidate must contain exactly 50 spiky dynamic effects.")
    }
    if (!pair_truth_rng_across_widths) {
      set.seed(component_seeds[["functional_truth"]] + 5000L + spiky_df)
    }
    for (j in spiky_index) {
      if (pair_truth_rng_across_widths) {
        set.seed(component_seeds[["functional_truth"]] + 5000L + j)
      }
      truth <- sample_multispike_local_bspline_truth(
        time_group = effect_sim$unit_info$time_group[j],
        switch_status = effect_sim$unit_info$switch_status[j],
        time_grid = time_grid,
        evaluation_grid = evaluation_grid,
        spike_count = effect_sim$unit_info$spike_count[j],
        df = spiky_df,
        degree = 3,
        secondary_fraction = secondary_fraction,
        switch_threshold = 0.25,
        minimum_location_margin = 0.60,
        minimum_location_ratio = 1.4,
        target_centered_rms = 0.90,
        minimum_peak_separation = 3,
        secondary_selection = secondary_selection,
        non_switch_baseline_fraction = 0.20
      )
      effect_sim$beta_matrix[j, ] <- truth$beta_observed
      effect_sim$beta_evaluation[j, ] <- truth$beta_evaluation
      effect_sim$unit_info$generation_attempt[j] <- truth$attempt
      effect_sim$unit_info$target_to_outside_ratio[j] <-
        truth$diagnostics$target_to_outside_ratio
      effect_sim$unit_info$effective_sign_transitions[j] <-
        truth$diagnostics$effective_sign_transitions
      effect_sim$unit_info$centered_rms[j] <- truth$diagnostics$centered_rms
      effect_sim$unit_info$spike_pattern[j] <- truth$diagnostics$spike_pattern
    }
    effect_sim$true_functionals <- evaluate_temporal_functionals(
      curves = effect_sim$beta_evaluation,
      smooth_var = evaluation_grid,
      switch_threshold = 0.25
    )
    effect_sim$settings$spiky_bspline_df <- spiky_df
  }
  spiky_index <- which(
    effect_sim$unit_info$effect_class == "dynamic_bspline" &
      effect_sim$unit_info$shape_profile == "spiky"
  )
  expected_spiky <- if (base_dgp == "reviewer_mixed") 98L else 50L
  if (length(spiky_index) != expected_spiky) {
    stop("Unexpected number of spiky dynamic effects.")
  }

  expression_sim <- simulate_eqtl_expression_from_genotypes(
    G = genotype_sim$G,
    beta_matrix = effect_sim$beta_matrix,
    time_grid = time_grid,
    covariates = covariates,
    expression_noise_sd = 1,
    seed = component_seeds[["expression"]]
  )
  message(
    "Running spiky-width candidate ",
    candidate_index,
    "/",
    length(spiky_dfs),
    ": df=",
    spiky_df
  )
  fit <- run_genotype_level_dynamic_eqtl_simulation(
    G = genotype_sim$G,
    time_grid = time_grid,
    covariates = covariates,
    class_probs = class_probs,
    expression_noise_sd = 1,
    alpha = 0.05,
    seed = seed,
    num_cores = num_cores,
    num_basis = num_basis,
    scenario = "internal_spiky_bf_width_screen",
    output_dir = output_dir,
    save_outputs = FALSE,
    verbose = FALSE,
    effect_sim = effect_sim,
    expression_sim = expression_sim
  )
  fit <- add_direct_interaction_efdr_results_to_genotype_output(
    out = fit,
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

  methods <- c(
    "FASH-IWP1-Raw",
    "FASH-IWP1-BF",
    "Direct-quadratic-LRT-eFDR-true-pi0"
  )
  result <- fit$result_table[fit$result_table$method %in% methods, ]
  geometry <- truth_geometry(effect_sim$beta_matrix, time_grid)
  unit_diagnostics <- cbind(
    effect_sim$unit_info,
    geometry,
    raw_lfdr = fit$fash_fits$fash_iwp1_raw$lfdr,
    bf_lfdr = fit$fash_fits$fash_iwp1_bf$lfdr,
    log_bayes_factor = log(fit$fash_fits$fash_iwp1_bf$BF)
  )
  unit_diagnostics$raw_cumulative_fdr <- cumulative_lfdr_score(
    unit_diagnostics$raw_lfdr
  )
  unit_diagnostics$bf_cumulative_fdr <- cumulative_lfdr_score(
    unit_diagnostics$bf_lfdr
  )
  unit_diagnostics$raw_selected <- unit_diagnostics$raw_cumulative_fdr <= 0.05
  unit_diagnostics$bf_selected <- unit_diagnostics$bf_cumulative_fdr <= 0.05
  direct_result <- result[
    result$method == "Direct-quadratic-LRT-eFDR-true-pi0",
  ]
  unit_diagnostics$direct_quadratic_qvalue <- direct_result$adjusted_score[
    match(unit_diagnostics$unit_index, direct_result$unit_index)
  ]
  unit_diagnostics$direct_quadratic_selected <-
    unit_diagnostics$direct_quadratic_qvalue <= 0.05

  spiky_diagnostics <- unit_diagnostics[spiky_index, ]
  summary <- data.frame(
    spiky_df = spiky_df,
    mean_support_25 = mean(spiky_diagnostics$support_25),
    mean_support_50 = mean(spiky_diagnostics$support_50),
    mean_second_difference_roughness =
      mean(spiky_diagnostics$second_difference_roughness),
    mean_quadratic_projection = mean(spiky_diagnostics$quadratic_projection),
    median_log_bayes_factor = median(spiky_diagnostics$log_bayes_factor),
    mean_log_bayes_factor = mean(spiky_diagnostics$log_bayes_factor),
    raw_spiky_power = mean(spiky_diagnostics$raw_selected),
    bf_spiky_power = mean(spiky_diagnostics$bf_selected),
    direct_quadratic_spiky_power =
      mean(spiky_diagnostics$direct_quadratic_selected),
    raw_pi0 = constant_component_prior_weight(fit$fash_fits$fash_iwp1_raw),
    bf_pi0 = constant_component_prior_weight(fit$fash_fits$fash_iwp1_bf),
    stringsAsFactors = FALSE
  )
  compact <- list(
    configuration = configuration,
    spiky_df = spiky_df,
    effect_sim = effect_sim,
    unit_diagnostics = unit_diagnostics,
    result_table = result,
    alpha_curve = fit$alpha_curve[fit$alpha_curve$method %in% methods, ],
    summary = summary,
    raw_prior_weights = fit$fash_fits$fash_iwp1_raw$prior_weights,
    bf_prior_weights = fit$fash_fits$fash_iwp1_bf$prior_weights
  )
  saveRDS(
    compact,
    file.path(candidate_dir, paste0("spiky_df", spiky_df, ".rds"))
  )
  candidate_results[[candidate_index]] <- summary
  write_csv(
    do.call(rbind, candidate_results[seq_len(candidate_index)]),
    file.path(output_dir, "width_screen_summary.csv")
  )
}

message("Saved internal spiky BF width screen to: ", output_dir)
