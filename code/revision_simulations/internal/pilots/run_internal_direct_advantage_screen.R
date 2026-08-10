#!/usr/bin/env Rscript

# Screen scientifically interpretable mixed-shape settings for a reproducible
# IWP FASH power advantage over low-degree direct interaction tests.

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

parse_numeric_list <- function(x) {
  out <- suppressWarnings(as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1]])))
  if (length(out) == 0 || any(!is.finite(out))) {
    stop("Numeric list arguments must contain finite comma-separated values.")
  }
  out
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

balanced_group_subset <- function(index, groups, size, seed) {
  if (size < 0 || size > length(index)) {
    stop("Requested balanced subset size is invalid.")
  }
  if (size == 0) return(integer())
  if (size == length(index)) return(sort(index))
  groups <- as.character(groups)
  group_levels <- sort(unique(groups))
  capacity <- vapply(group_levels, function(group) sum(groups == group), integer(1))
  expected <- size * capacity / sum(capacity)
  target <- floor(expected)
  remainder <- size - sum(target)
  if (remainder > 0) {
    allocation_order <- order(expected - target, decreasing = TRUE)
    for (k in allocation_order) {
      if (remainder == 0) break
      if (target[k] < capacity[k]) {
        target[k] <- target[k] + 1L
        remainder <- remainder - 1L
      }
    }
  }
  if (sum(target) != size || any(target > capacity)) {
    stop("Could not allocate the requested subset across shape strata.")
  }

  set.seed(seed)
  selected <- integer()
  for (k in seq_along(group_levels)) {
    candidates <- index[groups == group_levels[k]]
    if (length(candidates) < target[k]) {
      stop("A shape stratum does not contain enough candidates.")
    }
    selected <- c(selected, sample(candidates, target[k]))
  }
  sort(selected)
}

quadratic_projection_fraction <- function(curves, time_grid) {
  curves <- as.matrix(curves)
  basis <- scaled_time_polynomial(time_grid, degree = 2)
  basis_qr <- qr(basis)
  apply(curves, 1, function(curve) {
    centered <- curve - mean(curve)
    denominator <- sum(centered^2)
    if (denominator <= .Machine$double.eps) return(0)
    fitted <- qr.fitted(basis_qr, centered)
    sum(fitted^2) / denominator
  })
}

compute_internal_shape_power <- function(result_table,
                                         unit_info,
                                         methods,
                                         alpha_grid,
                                         seed) {
  dynamic <- unit_info$effect_class == "dynamic_bspline"
  shapes <- sort(unique(unit_info$shape_profile[dynamic]))
  rows <- vector("list", length(methods) * length(alpha_grid) * length(shapes))
  row_index <- 0L
  for (method in methods) {
    method_result <- result_table[result_table$method == method, , drop = FALSE]
    matched_unit <- match(method_result$unit_index, unit_info$unit_index)
    if (anyNA(matched_unit)) stop("Could not align result rows with unit metadata.")
    for (alpha in alpha_grid) {
      selected <- is.finite(method_result$adjusted_score) &
        method_result$adjusted_score <= alpha
      for (shape in shapes) {
        row_index <- row_index + 1L
        truth <- dynamic[matched_unit] &
          unit_info$shape_profile[matched_unit] == shape
        rows[[row_index]] <- data.frame(
          seed = seed,
          shape_profile = shape,
          method = method,
          alpha = alpha,
          n_dynamic = sum(truth),
          true_positives = sum(selected & truth),
          power = if (sum(truth) > 0) sum(selected & truth) / sum(truth) else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

J <- as.integer(get_arg("--J", "1000"))
n_donors <- as.integer(get_arg("--n-donors", "19"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
num_basis <- as.integer(get_arg("--num-basis", "20"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
expression_noise_sd <- as.numeric(get_arg("--noise-sd", "1"))
efdr_permutations <- as.integer(get_arg("--efdr-permutations", "30"))
seed <- as.integer(get_arg("--seed", "12345"))
broad_fractions <- parse_numeric_list(get_arg("--broad-fractions", "0.50,0.75"))
broad_dfs <- as.integer(parse_numeric_list(get_arg("--broad-dfs", "6,8,10")))
broad_amplitude <- as.numeric(get_arg("--broad-amplitude", "2"))
spiky_centered_rms <- as.numeric(get_arg("--spiky-centered-rms", "0.90"))
output_id <- get_arg("--output-id", "mixed_random_broad_internal_screen")

if (J < 100 || n_donors < n_covariates + 3 || n_covariates < 0 ||
    num_basis < 2 || num_cores < 1 || expression_noise_sd <= 0 ||
    efdr_permutations < 1 || any(broad_fractions < 0.5) ||
    any(broad_fractions >= 1) || any(broad_dfs <= 3) ||
    broad_amplitude <= 0 || spiky_centered_rms <= 0 || !nzchar(output_id)) {
  stop("Invalid internal-screen arguments.")
}

time_grid <- make_time_grid()
evaluation_grid <- seq(min(time_grid), max(time_grid), by = 0.1)
alpha_grid <- seq(0.005, 0.20, by = 0.005)
class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
true_pi0 <- unname(class_probs["constant"] + class_probs["zero"])
methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "Direct-linear-LRT-eFDR-true-pi0",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
component_seeds <- revision_component_seeds(seed)

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
base_effects <- simulate_targeted_local_bspline_effect_set(
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
  spiky_secondary_fraction = c(0.40, 0.65),
  spiky_minimum_peak_separation = 3,
  spiky_non_switch_baseline_fraction = 0.20,
  target_centered_rms = spiky_centered_rms,
  seed = component_seeds[["functional_truth"]],
  scenario = "internal_mixed_random_broad_spiky"
)

dynamic_index <- which(base_effects$unit_info$effect_class == "dynamic_bspline")
available_spiky <- dynamic_index[
  base_effects$unit_info$shape_profile[dynamic_index] == "spiky"
]
candidate_grid <- expand.grid(
  broad_fraction_requested = broad_fractions,
  broad_df = broad_dfs,
  stringsAsFactors = FALSE
)

configuration <- list(
  output_id = output_id,
  internal_only = TRUE,
  J = J,
  n_donors = n_donors,
  n_covariates = n_covariates,
  time_grid = time_grid,
  expression_noise_sd = expression_noise_sd,
  class_probs = class_probs,
  seed = seed,
  efdr_permutations = efdr_permutations,
  broad_fractions = broad_fractions,
  broad_dfs = broad_dfs,
  broad_amplitude = broad_amplitude,
  spiky_centered_rms = spiky_centered_rms,
  true_pi0 = true_pi0,
  paired_genotype_covariates_noise_and_permutations = TRUE
)
saveRDS(configuration, file.path(output_dir, "configuration.rds"))

candidate_results <- vector("list", nrow(candidate_grid))
for (candidate_index in seq_len(nrow(candidate_grid))) {
  setting <- candidate_grid[candidate_index, ]
  n_spiky_target <- min(
    length(available_spiky),
    round(length(dynamic_index) * (1 - setting$broad_fraction_requested))
  )
  keep_spiky <- balanced_group_subset(
    index = available_spiky,
    groups = base_effects$unit_info$truth_stratum[available_spiky],
    size = n_spiky_target,
    seed = component_seeds[["functional_truth"]] + candidate_index
  )
  broad_index <- setdiff(dynamic_index, keep_spiky)

  effect_sim <- base_effects
  set.seed(component_seeds[["functional_truth"]] + 1000L + setting$broad_df)
  for (j in broad_index) {
    effect_sim$beta_matrix[j, ] <- simulate_bspline_effect(
      x = time_grid,
      amplitude = broad_amplitude,
      df = setting$broad_df,
      coefficient_sd = 1
    )
  }
  effect_sim$unit_info$shape_profile[dynamic_index] <- "broad_random"
  effect_sim$unit_info$shape_profile[keep_spiky] <- "spiky"
  effect_sim$unit_info$spike_count[broad_index] <- NA_integer_
  effect_sim$unit_info$spike_pattern[broad_index] <- NA_character_
  effect_sim$unit_info$truth_stratum[broad_index] <- "broad_random"
  effect_sim$true_functionals <- evaluate_temporal_functionals(
    curves = effect_sim$beta_matrix,
    smooth_var = time_grid,
    switch_threshold = 0.25
  )

  expression_sim <- simulate_eqtl_expression_from_genotypes(
    G = genotype_sim$G,
    beta_matrix = effect_sim$beta_matrix,
    time_grid = time_grid,
    covariates = covariates,
    expression_noise_sd = expression_noise_sd,
    seed = component_seeds[["expression"]]
  )
  message(
    "Running candidate ",
    candidate_index,
    "/",
    nrow(candidate_grid),
    ": broad fraction ",
    round(length(broad_index) / length(dynamic_index), 3),
    ", df ",
    setting$broad_df
  )
  fit <- run_genotype_level_dynamic_eqtl_simulation(
    G = genotype_sim$G,
    time_grid = time_grid,
    covariates = covariates,
    class_probs = class_probs,
    expression_noise_sd = expression_noise_sd,
    alpha = 0.05,
    seed = seed,
    num_cores = num_cores,
    num_basis = num_basis,
    scenario = "internal_mixed_random_broad_spiky",
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

  alpha_curve <- fit$alpha_curve[fit$alpha_curve$method %in% methods, ]
  alpha_005 <- alpha_curve[abs(alpha_curve$alpha - 0.05) < 1e-12, ]
  projection <- quadratic_projection_fraction(
    effect_sim$beta_matrix[dynamic_index, , drop = FALSE],
    time_grid = time_grid
  )
  shape <- effect_sim$unit_info$shape_profile[dynamic_index]
  projection_summary <- data.frame(
    shape_profile = c("all", "broad_random", "spiky"),
    n_dynamic = c(
      length(dynamic_index),
      sum(shape == "broad_random"),
      sum(shape == "spiky")
    ),
    mean_quadratic_projection = c(
      mean(projection),
      mean(projection[shape == "broad_random"]),
      mean(projection[shape == "spiky"])
    ),
    median_quadratic_projection = c(
      median(projection),
      median(projection[shape == "broad_random"]),
      median(projection[shape == "spiky"])
    ),
    stringsAsFactors = FALSE
  )
  shape_power <- compute_internal_shape_power(
    result_table = fit$result_table,
    unit_info = effect_sim$unit_info,
    methods = methods,
    alpha_grid = alpha_grid,
    seed = seed
  )

  candidate_id <- sprintf(
    "broad%03d_df%02d",
    round(100 * length(broad_index) / length(dynamic_index)),
    setting$broad_df
  )
  compact <- list(
    candidate_id = candidate_id,
    setting = setting,
    broad_fraction_realized = length(broad_index) / length(dynamic_index),
    broad_index = broad_index,
    spiky_index = keep_spiky,
    effect_sim = effect_sim,
    alpha_curve = alpha_curve,
    alpha_005 = alpha_005,
    shape_power = shape_power,
    projection_summary = projection_summary,
    estimated_pi0 = data.frame(
      fit = c("raw", "BF-corrected"),
      estimated_pi0 = c(
        constant_component_prior_weight(fit$fash_fits$fash_iwp1_raw),
        constant_component_prior_weight(fit$fash_fits$fash_iwp1_bf)
      ),
      stringsAsFactors = FALSE
    )
  )
  saveRDS(compact, file.path(candidate_dir, paste0(candidate_id, ".rds")))

  comparison <- alpha_005[, c(
    "method", "n_discoveries", "false_discoveries", "power", "empirical_fdr"
  )]
  comparison$candidate_id <- candidate_id
  comparison$broad_fraction <- compact$broad_fraction_realized
  comparison$broad_df <- setting$broad_df
  comparison$mean_quadratic_projection <- projection_summary$mean_quadratic_projection[1]
  candidate_results[[candidate_index]] <- comparison[, c(
    "candidate_id", "broad_fraction", "broad_df",
    "mean_quadratic_projection", "method", "n_discoveries",
    "false_discoveries", "power", "empirical_fdr"
  )]
  write_csv(
    do.call(rbind, candidate_results[seq_len(candidate_index)]),
    file.path(output_dir, "screen_results_alpha005.csv")
  )
}

screen_results <- do.call(rbind, candidate_results)
wide_power <- reshape(
  screen_results[, c("candidate_id", "method", "power")],
  idvar = "candidate_id",
  timevar = "method",
  direction = "wide"
)
names(wide_power) <- sub("^power\\.", "", names(wide_power))
wide_power$raw_iwp_minus_direct_quadratic <-
  wide_power[["FASH-IWP1-Raw"]] -
  wide_power[["Direct-quadratic-LRT-eFDR-true-pi0"]]
wide_power$bf_iwp_minus_direct_quadratic <-
  wide_power[["FASH-IWP1-BF"]] -
  wide_power[["Direct-quadratic-LRT-eFDR-true-pi0"]]
wide_power <- wide_power[
  order(
    -wide_power$bf_iwp_minus_direct_quadratic,
    -wide_power$raw_iwp_minus_direct_quadratic
  ),
]
write_csv(wide_power, file.path(output_dir, "candidate_ranking.csv"))

message("Saved internal direct-interaction screen to: ", output_dir)
