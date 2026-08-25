#!/usr/bin/env Rscript

# Run compact five-seed replications for the formal real-genotype compact
# raised-cosine one-/two-/three-peak simulation.

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

resolve_input_path <- function(path, workflowr_root, argument_name) {
  candidates <- unique(c(path, file.path(workflowr_root, path)))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop(argument_name, " does not exist: ", path)
  }
  normalizePath(existing[[1L]], winslash = "/", mustWork = TRUE)
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
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
))

J <- as.integer(get_arg("--J", "6362"))
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
    "main_effect_linear_mixture_predstep1_penalty10_pilot5"
  )
} else {
  paste0(
    "r2_real_genotype_one_per_gene_J6362_",
    "timed_cosine_one_two_three_peak_",
    "main_effect_linear_mixture_predstep1_penalty10_pilot5"
  )
}
output_id <- get_arg("--output-id", default_output_id)
default_genotype_cache_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "shared",
  "real_genotype_one_per_gene_J6362_pilot5",
  "genotype_samples.rds"
)
genotype_cache_path <- resolve_input_path(
  get_arg("--genotype-cache", default_genotype_cache_path),
  workflowr_root,
  "--genotype-cache"
)
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

genotype_cache <- readRDS(genotype_cache_path)
if (!is.list(genotype_cache) ||
    !all(c("configuration", "sample_ids", "samples") %in% names(genotype_cache)) ||
    !identical(genotype_cache$configuration$n_genes, J) ||
    !identical(genotype_cache$configuration$n_donors, n_donors) ||
    !all(seed_list %in% genotype_cache$configuration$seed_list) ||
    !all(as.character(seed_list) %in% names(genotype_cache$samples))) {
  stop("The shared real-genotype cache does not match J, donors, or seeds.")
}
for (seed in seed_list) {
  sample <- validate_real_genotype_sample(
    genotype_cache$samples[[as.character(seed)]],
    expected_genes = J,
    expected_donors = n_donors,
    maf_min = genotype_cache$configuration$maf_min
  )
  expected_digest <- object_md5(list(
    pair_key = sample$selection$pair_key,
    sample_ids = rownames(sample$G),
    G = sample$G
  ))
  if (!identical(sample$genotype_digest, expected_digest)) {
    stop("The shared genotype digest is invalid for seed ", seed, ".")
  }
}
genotype_cache_fingerprint <- artifact_fingerprint(genotype_cache_path)

time_grid <- make_time_grid()
evaluation_grid <- seq(min(time_grid), max(time_grid), by = 0.1)
class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
expected_class_counts <- exact_proportional_counts(J, class_probs)
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
  paste0(
    "r2_real_genotype_one_per_gene_centered_",
    "timed_cosine_one_two_three_peak_main_effect_dynamic_eqtl"
  )
} else {
  paste0(
    "r2_real_genotype_one_per_gene_",
    "timed_cosine_one_two_three_peak_main_effect_dynamic_eqtl"
  )
}
true_pi0 <- unname(sum(expected_class_counts[c("constant", "zero")]) / J)
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
  expected_class_counts = expected_class_counts,
  spike_counts = 1:3,
  shape_cell_probs = shape_cell_probs,
  shape_cell_counts = exact_proportional_counts(
    expected_class_counts[["dynamic_bspline"]],
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
  linear_prior_mode = linear_prior_mode,
  common_sd_grid = common_sd_grid,
  common_pred_step = common_pred_step,
  common_penalty = common_penalty,
  genotype_source = "paper-derived YRI DS dosage",
  genotype_selection_rule = genotype_cache$configuration$selection_rule,
  genotype_repeated_variant_rule = genotype_cache$configuration$repeated_variant_rule,
  genotype_maf_min = genotype_cache$configuration$maf_min,
  genotype_sample_ids = genotype_cache$sample_ids,
  genotype_cache_fingerprint = genotype_cache_fingerprint,
  genotype_source_configuration = genotype_cache$configuration,
  maf_truth_balance_method = paste(
    "exact global class counts with within-MAF-decile permutation"
  ),
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
  genotype_sample <- validate_real_genotype_sample(
    genotype_cache$samples[[as.character(seed)]],
    expected_genes = J,
    expected_donors = n_donors,
    maf_min = genotype_cache$configuration$maf_min
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
  effect_sim <- reassign_effect_simulation_by_maf(
    effect_sim = effect_sim,
    maf = genotype_sample$variant_info$observed_maf,
    class_probs = class_probs,
    seed = component_seeds[["classes"]],
    n_strata = 10L
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
    G = genotype_sample$G,
    beta_matrix = effect_sim$beta_matrix,
    time_grid = time_grid,
    covariates = covariates,
    expression_noise_sd = expression_noise_sd,
    seed = component_seeds[["expression"]]
  )
  out <- run_genotype_level_dynamic_eqtl_simulation(
    G = genotype_sample$G,
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
    output_dir = output_dir,
    save_outputs = FALSE,
    verbose = FALSE,
    effect_sim = effect_sim,
    expression_sim = expression_sim
  )
  out$settings$genotype_source <- configuration$genotype_source
  out$settings$genotype_selection_rule <- configuration$genotype_selection_rule
  out$settings$genotype_digest <- genotype_sample$genotype_digest
  truth_maf_balance <- summarize_truth_maf_balance(
    variant_info = genotype_sample$variant_info,
    unit_info = out$unit_info,
    seed = seed
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
    out$genotype_metadata <- genotype_sample$variant_info
    out$genotype_selection <- genotype_sample$selection
    out$genotype_source_configuration <- genotype_cache$configuration
    out$genotype_digest <- genotype_sample$genotype_digest
    out$truth_maf_balance <- truth_maf_balance
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
    genotype_digest = genotype_sample$genotype_digest,
    selected_pair_keys = genotype_sample$selection$pair_key,
    genotype_selection_summary = genotype_sample$selection_summary,
    truth_maf_balance = truth_maf_balance,
    alpha_curve = alpha_curve,
    alpha_005 = alpha_curve[abs(alpha_curve$alpha - 0.05) < 1e-12, ],
    peak_alpha_curve = peak_alpha_curve,
    peak_alpha_005 = peak_alpha_curve[
      abs(peak_alpha_curve$alpha - 0.05) < 1e-12,
    ],
    pi0 = pi0,
    geometry = geometry,
    linear_prior_weights = linear_prior_weights,
    linear_prior_summary = linear_prior_summary
  )
}

validate_replicate <- function(x, seed) {
  required <- c(
    "configuration", "seed", "component_seeds", "alpha_curve", "alpha_005",
    "genotype_digest", "selected_pair_keys", "genotype_selection_summary",
    "truth_maf_balance",
    "peak_alpha_curve", "peak_alpha_005", "pi0", "geometry",
    "linear_prior_weights", "linear_prior_summary"
  )
  prior_weights <- x$linear_prior_weights
  prior_summary <- x$linear_prior_summary
  weight_groups <- split(prior_weights, prior_weights$fit)
  weights_valid <- length(weight_groups) == 2L &&
    all(vapply(weight_groups, function(rows) {
      identical(rows$seed, rep(seed, length(common_sd_grid))) &&
        isTRUE(all.equal(
          rows$predstep_sd,
          common_sd_grid,
          tolerance = 0
        )) &&
        sum(rows$is_null) == 1L &&
        rows$is_null[1] &&
        all(is.finite(rows$prior_weight)) &&
        all(rows$prior_weight >= 0) &&
        abs(sum(rows$prior_weight) - 1) < 1e-6 &&
        identical(rows$active, rows$prior_weight > 0)
    }, logical(1)))
  all(required %in% names(x)) &&
    identical(x$seed, seed) &&
    isTRUE(all.equal(x$configuration, configuration)) &&
    identical(
      x$genotype_digest,
      genotype_cache$samples[[as.character(seed)]]$genotype_digest
    ) &&
    identical(
      x$selected_pair_keys,
      genotype_cache$samples[[as.character(seed)]]$selection$pair_key
    ) &&
    is.data.frame(x$genotype_selection_summary) &&
    nrow(x$genotype_selection_summary) == 1L &&
    x$genotype_selection_summary$genes == J &&
    is.data.frame(x$truth_maf_balance) &&
    nrow(x$truth_maf_balance) == length(class_probs) &&
    all(
      x$truth_maf_balance$n ==
        exact_proportional_counts(J, class_probs)[
          x$truth_maf_balance$effect_class
        ]
    ) &&
    all(all_methods %in% unique(x$alpha_curve$method)) &&
    all(all_methods %in% unique(x$peak_alpha_curve$method)) &&
    identical(sort(unique(x$peak_alpha_curve$subgroup_value)), c("1", "2", "3")) &&
    is.data.frame(prior_weights) &&
    nrow(prior_weights) == 2L * length(common_sd_grid) &&
    weights_valid &&
    is.data.frame(prior_summary) &&
    nrow(prior_summary) == 2L &&
    identical(prior_summary$seed, rep(seed, 2L)) &&
    setequal(prior_summary$fit, c("Raw", "BF-corrected")) &&
    all(is.finite(prior_summary$estimated_pi0)) &&
    all(prior_summary$estimated_pi0 >= 0) &&
    all(prior_summary$estimated_pi0 <= 1) &&
    all(prior_summary$active_nonnull_components >= 0) &&
    all(
      is.finite(prior_summary$alternative_rms_predstep_sd) |
        (
          is.na(prior_summary$alternative_rms_predstep_sd) &
            prior_summary$active_nonnull_components == 0
        )
    )
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
all_linear_prior_weights <- do.call(
  rbind,
  lapply(replicates, `[[`, "linear_prior_weights")
)
all_linear_prior_summary <- do.call(
  rbind,
  lapply(replicates, `[[`, "linear_prior_summary")
)
all_genotype_selection_summary <- do.call(
  rbind,
  lapply(replicates, `[[`, "genotype_selection_summary")
)
all_truth_maf_balance <- do.call(
  rbind,
  lapply(replicates, `[[`, "truth_maf_balance")
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
  all_linear_prior_weights,
  file.path(summary_dir, "all_replicate_linear_prior_weights.csv")
)
write_csv(
  all_linear_prior_summary,
  file.path(summary_dir, "all_replicate_linear_prior_summary.csv")
)
write_csv(
  all_genotype_selection_summary,
  file.path(summary_dir, "genotype_selection_summary.csv")
)
write_csv(
  all_truth_maf_balance,
  file.path(summary_dir, "truth_maf_balance.csv")
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
print(all_linear_prior_summary)
