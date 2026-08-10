#!/usr/bin/env Rscript

# Run compact multi-seed replications for early, middle, late, and switch
# functional testing under balanced labelled B-spline truth mechanisms.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) return(".")
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find code/revision_simulations/shared/simulation_functions.R from the current working directory.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0 || hit[1] == length(args)) return(default)
  args[hit[1] + 1]
}

as_flag <- function(x) tolower(x) %in% c("1", "true", "t", "yes", "y")

parse_seed_list <- function(x) {
  pieces <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  if (length(pieces) == 0 || any(!nzchar(pieces))) {
    stop("--seed-list must contain one or more comma-separated integer seeds.")
  }
  seeds <- suppressWarnings(as.integer(pieces))
  if (any(is.na(seeds)) || anyDuplicated(seeds)) {
    stop("--seed-list must contain unique integer seeds.")
  }
  seeds
}

write_csv <- function(x, path) write.csv(x, file = path, row.names = FALSE)

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

J <- as.integer(get_arg("--J", "1000"))
n_donors <- as.integer(get_arg("--n-donors", "19"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
expression_noise_sd <- as.numeric(get_arg("--noise-sd", "1"))
num_basis <- as.integer(get_arg("--num-basis", "20"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
include_global_comparisons <- as_flag(get_arg("--include-global-comparisons", "false"))
efdr_permutations <- as.integer(get_arg("--efdr-permutations", "100"))
dynamic_amplitude <- as.numeric(get_arg("--dynamic-amplitude", "2"))
bspline_df <- as.integer(get_arg("--bspline-df", "6"))
switch_threshold <- as.numeric(get_arg("--switch-threshold", "0.25"))
truth_margin <- as.numeric(get_arg("--truth-margin", "0.10"))
non_switch_baseline_margin <- as.numeric(get_arg("--non-switch-baseline-margin", "0.35"))
truth_mechanism <- get_arg("--truth-mechanism", "random_bspline")
minimum_location_margin <- as.numeric(get_arg("--minimum-location-margin", "0.60"))
minimum_location_ratio <- as.numeric(get_arg("--minimum-location-ratio", "2"))
targeted_profile <- get_arg("--targeted-profile", "narrow")
spiky_truth_version <- get_arg("--spiky-truth-version", "centered_single_v1")
spiky_secondary_fraction <- c(
  as.numeric(get_arg("--spiky-secondary-fraction-min", "0.40")),
  as.numeric(get_arg("--spiky-secondary-fraction-max", "0.65"))
)
spiky_minimum_peak_separation <- as.numeric(get_arg(
  "--spiky-minimum-peak-separation",
  "3"
))
spiky_non_switch_baseline_fraction <- as.numeric(get_arg(
  "--spiky-non-switch-baseline-fraction",
  "0"
))
target_centered_rms_arg <- get_arg("--target-centered-rms", "none")
target_centered_rms <- if (tolower(target_centered_rms_arg) %in% c("none", "null")) {
  NULL
} else {
  as.numeric(target_centered_rms_arg)
}
non_switch_baseline_fraction <- as.numeric(get_arg("--non-switch-baseline-fraction", "0.15"))
non_switch_background_fraction <- as.numeric(get_arg("--non-switch-background-fraction", "0.05"))
n_examples_per_group <- as.integer(get_arg("--n-examples-per-group", "2"))
seed_list <- parse_seed_list(get_arg(
  "--seed-list",
  "12345,22345,32345,42345,52345"
))
output_id <- get_arg("--output-id", "functional_bspline_pilot5")
overwrite <- as_flag(get_arg("--overwrite", "false"))
cache_full_fits <- as_flag(get_arg("--cache-full-fits", "false"))

if (!truth_mechanism %in% c("random_bspline", "targeted_local_bspline") ||
    !targeted_profile %in% c("narrow", "broad", "mixed") ||
    !spiky_truth_version %in% c("centered_single_v1", "mixed_single_double_v2") ||
    (spiky_truth_version == "mixed_single_double_v2" &&
      (truth_mechanism != "targeted_local_bspline" ||
        targeted_profile != "mixed")) ||
    J < 30 || n_donors < n_covariates + 3 || n_covariates < 0 ||
    !is.finite(expression_noise_sd) || expression_noise_sd <= 0 ||
    num_basis < 2 || num_cores < 1 || dynamic_amplitude <= 0 ||
    !is.finite(efdr_permutations) || efdr_permutations < 1 ||
    bspline_df <= 3 || switch_threshold <= 0 || truth_margin <= 0 ||
    non_switch_baseline_margin <= 0 || minimum_location_margin <= 0 ||
    minimum_location_ratio <= 1 || non_switch_baseline_fraction <= 0 ||
    non_switch_background_fraction < 0 ||
    length(spiky_secondary_fraction) != 2 ||
    any(!is.finite(spiky_secondary_fraction)) ||
    any(spiky_secondary_fraction <= 0) ||
    spiky_secondary_fraction[1] > spiky_secondary_fraction[2] ||
    spiky_secondary_fraction[2] >= 1 ||
    !is.finite(spiky_minimum_peak_separation) ||
    spiky_minimum_peak_separation < 0 ||
    !is.finite(spiky_non_switch_baseline_fraction) ||
    spiky_non_switch_baseline_fraction < 0 ||
    (!is.null(target_centered_rms) &&
      (length(target_centered_rms) != 1 || !is.finite(target_centered_rms) ||
        target_centered_rms <= 0)) ||
    n_examples_per_group < 1 || !nzchar(output_id)) {
  stop("Invalid simulation arguments.")
}

time_grid <- make_time_grid()
evaluation_grid <- seq(min(time_grid), max(time_grid), by = 0.1)
alpha_grid <- seq(0.005, 0.20, by = 0.005)
class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
scenario <- if (truth_mechanism == "random_bspline") {
  "genotype_functional_random_bspline"
} else {
  "genotype_functional_targeted_local_bspline"
}
true_pi0_expected <- unname(class_probs["constant"] + class_probs["zero"])
methods <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
targets <- c("early", "middle", "late", "switch")
global_fash_methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "FASH-linear-Raw",
  "FASH-linear-BF"
)
global_direct_methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "Direct-linear-LRT-eFDR-true-pi0",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
global_methods <- unique(c(global_fash_methods, global_direct_methods))

output_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "mc", output_id
)
replicate_dir <- file.path(output_dir, "replicates")
summary_dir <- file.path(output_dir, "summary")
figure_dir <- file.path(output_dir, "figures")
full_fit_dir <- file.path(output_dir, "full_fits")
invisible(lapply(
  c(
    output_dir,
    replicate_dir,
    summary_dir,
    figure_dir,
    if (cache_full_fits) full_fit_dir
  ),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

configuration <- list(
  output_id = output_id,
  scenario = scenario,
  J = J,
  n_donors = n_donors,
  time_grid = time_grid,
  evaluation_grid = evaluation_grid,
  n_covariates = n_covariates,
  expression_noise_sd = expression_noise_sd,
  num_basis = num_basis,
  class_probs = class_probs,
  truth_mechanism = truth_mechanism,
  targeted_profile = targeted_profile,
  spiky_truth_version = spiky_truth_version,
  spiky_secondary_fraction = spiky_secondary_fraction,
  spiky_minimum_peak_separation = spiky_minimum_peak_separation,
  spiky_non_switch_baseline_fraction = spiky_non_switch_baseline_fraction,
  dynamic_amplitude = dynamic_amplitude,
  bspline_df = bspline_df,
  switch_threshold = switch_threshold,
  truth_margin = truth_margin,
  non_switch_baseline_margin = non_switch_baseline_margin,
  minimum_location_margin = minimum_location_margin,
  minimum_location_ratio = minimum_location_ratio,
  non_switch_baseline_fraction = non_switch_baseline_fraction,
  non_switch_background_fraction = non_switch_background_fraction,
  target_centered_rms = target_centered_rms,
  alpha_grid = alpha_grid,
  n_examples_per_group = n_examples_per_group,
  true_pi0 = true_pi0_expected,
  seed_list = seed_list
)
if (include_global_comparisons) {
  configuration$global_comparisons <- list(
    efdr_permutations = efdr_permutations,
    direct_efdr_uses_true_pi0 = TRUE
  )
}
configuration_matches <- function(candidate_configuration) {
  isTRUE(all.equal(candidate_configuration, configuration)) ||
    (!include_global_comparisons && truth_mechanism == "random_bspline" &&
      isTRUE(all.equal(
        candidate_configuration,
        configuration[names(candidate_configuration)]
      ))) ||
    (!include_global_comparisons &&
      isTRUE(all.equal(
        candidate_configuration,
        configuration[setdiff(names(configuration), "global_comparisons")]
      )))
}
configuration_path <- file.path(output_dir, "configuration.rds")
if (file.exists(configuration_path) && !overwrite) {
  cached_configuration <- readRDS(configuration_path)
  if (!configuration_matches(cached_configuration)) {
    stop("The existing output id has different settings. Choose a new --output-id or use --overwrite true.")
  }
} else {
  saveRDS(configuration, configuration_path)
}

extract_truth_examples <- function(out, effect_sim) {
  mixed_shape <- identical(targeted_profile, "mixed")
  group_order <- if (mixed_shape) {
    as.vector(outer(
      c(
        "early / switch", "early / non-switch",
        "middle / switch", "middle / non-switch",
        "late / switch", "late / non-switch"
      ),
      c("broad", "spiky"),
      paste,
      sep = " / "
    ))
  } else {
    c(
      "early / switch", "early / non-switch",
      "middle / switch", "middle / non-switch",
      "late / switch", "late / non-switch"
    )
  }
  examples <- lapply(group_order, function(group) {
    group_column <- if (mixed_shape) "truth_stratum" else "truth_group"
    available <- which(out$unit_info[[group_column]] == group)
    multispike_non_switch <- mixed_shape &&
      spiky_truth_version == "mixed_single_double_v2" &&
      grepl(" / non-switch / spiky$", group)
    selected <- if (multispike_non_switch) {
      unlist(lapply(
        c("single", "same-sign double"),
        function(pattern) {
          head(
            available[out$unit_info$spike_pattern[available] == pattern],
            n_examples_per_group
          )
        }
      ))
    } else {
      head(available, n_examples_per_group)
    }
    expected_examples <- n_examples_per_group * if (multispike_non_switch) 2L else 1L
    if (length(selected) != expected_examples) {
      stop("Could not extract the requested number of examples for ", group, ".")
    }
    lapply(selected, function(index) {
      list(
        group = group,
        shape_profile = out$unit_info$shape_profile[index],
        spike_count = out$unit_info$spike_count[index],
        spike_pattern = out$unit_info$spike_pattern[index],
        variant_id = out$unit_info$variant_id[index],
        observed = data.frame(
          time = time_grid,
          estimate = out$eqtl_summary$beta_hat[index, ],
          se = out$eqtl_summary$se[index, ],
          true_effect = effect_sim$beta_matrix[index, ],
          stringsAsFactors = FALSE
        ),
        true_curve = data.frame(
          time = evaluation_grid,
          true_effect = effect_sim$beta_evaluation[index, ],
          stringsAsFactors = FALSE
        ),
        true_functionals = effect_sim$true_functionals[index, ]
      )
    })
  })
  names(examples) <- group_order
  examples
}

compute_shape_stratified_power <- function(result_table,
                                           unit_info,
                                           method_order,
                                           alpha_values,
                                           seed) {
  if (!"shape_profile" %in% names(unit_info)) {
    return(NULL)
  }
  shape_levels <- c("broad", "spiky")
  dynamic_by_shape <- lapply(shape_levels, function(shape) {
    unit_info$unit_index[
      unit_info$effect_class == "dynamic_bspline" &
        unit_info$shape_profile == shape
    ]
  })
  names(dynamic_by_shape) <- shape_levels
  if (any(lengths(dynamic_by_shape) == 0)) {
    stop("Mixed-shape truth is missing a broad or spiky dynamic stratum.")
  }

  rows <- list()
  row_index <- 1L
  for (method in method_order) {
    method_results <- result_table[result_table$method == method, , drop = FALSE]
    if (nrow(method_results) != nrow(unit_info)) {
      stop("Every shape-stratified method must contain one row per variant.")
    }
    for (shape in shape_levels) {
      alternative_indices <- dynamic_by_shape[[shape]]
      for (alpha in alpha_values) {
        selected_indices <- method_results$unit_index[
          method_results$adjusted_score <= alpha
        ]
        true_positives <- sum(alternative_indices %in% selected_indices)
        rows[[row_index]] <- data.frame(
          seed = seed,
          scenario = method_results$scenario[1],
          target = "dynamic",
          shape_profile = shape,
          method = method,
          alpha = alpha,
          n_true_alternatives = length(alternative_indices),
          true_positives = true_positives,
          power = true_positives / length(alternative_indices),
          stringsAsFactors = FALSE
        )
        row_index <- row_index + 1L
      }
    }
  }
  do.call(rbind, rows)
}

summarize_mc_shape_power <- function(shape_rows,
                                     confidence_level = 0.95) {
  split_rows <- split(
    shape_rows,
    list(shape_rows$shape_profile, shape_rows$method, shape_rows$alpha),
    drop = TRUE
  )
  summaries <- lapply(split_rows, function(x) {
    power_summary <- summarize_mc_values(x$power, confidence_level)
    tp_summary <- summarize_mc_values(x$true_positives, confidence_level)
    data.frame(
      scenario = x$scenario[1],
      target = x$target[1],
      shape_profile = x$shape_profile[1],
      method = x$method[1],
      alpha = x$alpha[1],
      n_true_alternatives = x$n_true_alternatives[1],
      n_replications = nrow(x),
      mean_true_positives = tp_summary[["mean"]],
      mean_power = power_summary[["mean"]],
      power_ci_lower = pmax(0, power_summary[["lower"]]),
      power_ci_upper = pmin(1, power_summary[["upper"]]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, summaries)
  out$method_rank <- rank_revision_methods(out$method)
  out <- out[order(
    match(out$shape_profile, c("broad", "spiky")),
    out$method_rank,
    out$method,
    out$alpha
  ), , drop = FALSE]
  out$method_rank <- NULL
  rownames(out) <- NULL
  out
}

extract_functional_call_diagnostics <- function(functional_result,
                                                method,
                                                unit_info,
                                                true_functionals,
                                                alpha = 0.05) {
  targets_to_audit <- colnames(true_functionals)
  rows <- vector("list", length(targets_to_audit))
  for (target_index in seq_along(targets_to_audit)) {
    target <- targets_to_audit[target_index]
    dynamic_indices <- sort(unique(
      functional_result$fdr_table$index[
        functional_result$fdr_table$FDR <= alpha
      ]
    ))
    lfsr_map <- functional_result$lfsr_by_target[[target]]
    candidate_lfsr <- unname(lfsr_map[as.character(dynamic_indices)])
    cfsr_table <- functional_cfsr_table(dynamic_indices, candidate_lfsr)
    selected <- cfsr_table[cfsr_table$cfsr <= alpha, , drop = FALSE]
    if (nrow(selected) == 0) {
      rows[[target_index]] <- NULL
      next
    }
    metadata <- unit_info[selected$index, c(
      "unit_index", "variant_id", "effect_class", "time_group",
      "switch_status", "shape_profile", "spike_count", "spike_pattern"
    ), drop = FALSE]
    rows[[target_index]] <- data.frame(
      metadata,
      target = target,
      method = method,
      alpha = alpha,
      lfsr = selected$lfsr,
      cfsr = selected$cfsr,
      true_functional = true_functionals[selected$index, target],
      true_null = true_functionals[selected$index, target] <= 0,
      stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    return(data.frame())
  }
  rownames(out) <- NULL
  out
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
  effect_sim <- if (truth_mechanism == "random_bspline") {
    simulate_labeled_random_bspline_effect_set(
      n_variants = J,
      time_grid = time_grid,
      evaluation_grid = evaluation_grid,
      class_probs = class_probs,
      dynamic_amplitude = dynamic_amplitude,
      bspline_df = bspline_df,
      switch_threshold = switch_threshold,
      truth_margin = truth_margin,
      non_switch_baseline_margin = non_switch_baseline_margin,
      seed = component_seeds[["functional_truth"]],
      scenario = scenario
    )
  } else {
    simulate_targeted_local_bspline_effect_set(
      n_variants = J,
      time_grid = time_grid,
      evaluation_grid = evaluation_grid,
      class_probs = class_probs,
      dynamic_amplitude = dynamic_amplitude,
      switch_threshold = switch_threshold,
      minimum_location_margin = minimum_location_margin,
      minimum_location_ratio = minimum_location_ratio,
      non_switch_baseline_fraction = non_switch_baseline_fraction,
      non_switch_background_fraction = non_switch_background_fraction,
      profile = targeted_profile,
      spiky_truth_version = spiky_truth_version,
      spiky_secondary_fraction = spiky_secondary_fraction,
      spiky_minimum_peak_separation = spiky_minimum_peak_separation,
      spiky_non_switch_baseline_fraction = spiky_non_switch_baseline_fraction,
      target_centered_rms = target_centered_rms,
      seed = component_seeds[["functional_truth"]],
      scenario = scenario
    )
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
    dynamic_amplitude = dynamic_amplitude,
    bspline_df = bspline_df,
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
  true_pi0 <- dynamic_null_proportion(out$unit_info, target = "dynamic")
  if (!isTRUE(all.equal(true_pi0, true_pi0_expected))) {
    stop("The simulated dynamic-null proportion does not match the configured mixture.")
  }

  global_results <- NULL
  if (include_global_comparisons) {
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
    missing_global_methods <- setdiff(global_methods, unique(out$result_table$method))
    if (length(missing_global_methods) > 0) {
      stop(
        "Global genotype-level output is missing reviewer-facing methods: ",
        paste(missing_global_methods, collapse = ", ")
      )
    }
    global_alpha <- out$alpha_curve[out$alpha_curve$method %in% global_methods, , drop = FALSE]
    global_alpha$seed <- seed
    shape_power <- if (identical(targeted_profile, "mixed")) {
      compute_shape_stratified_power(
        result_table = out$result_table,
        unit_info = out$unit_info,
        method_order = global_methods,
        alpha_values = alpha_grid,
        seed = seed
      )
    } else {
      NULL
    }
    global_results <- list(
      alpha_curve = global_alpha,
      alpha_005 = global_alpha[abs(global_alpha$alpha - 0.05) < 1e-12, , drop = FALSE],
      shape_power = shape_power,
      pi0 = data.frame(
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
    )
    if (nrow(global_results$alpha_005) != length(global_methods)) {
      stop("Global genotype-level output is missing alpha = 0.05 results.")
    }
  }

  true_functionals <- effect_sim$true_functionals[out$unit_info$variant_id, targets, drop = FALSE]
  dynamic_truth <- out$unit_info$effect_class == "dynamic_bspline"
  raw_results <- evaluate_fash_functional_testing(
    fit = out$fash_fits$fash_iwp1_raw,
    true_functionals = true_functionals,
    evaluation_grid = evaluation_grid,
    alpha_grid = alpha_grid,
    method = "FASH-IWP1-Raw",
    scenario = scenario,
    switch_threshold = switch_threshold,
    true_dynamic = dynamic_truth,
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
    switch_threshold = switch_threshold,
    true_dynamic = dynamic_truth,
    num_cores = num_cores,
    seed = component_seeds[["functional_posterior"]] + 100L
  )
  functional_alpha <- rbind(raw_results$alpha_curve, bf_results$alpha_curve)
  functional_alpha$seed <- seed
  if (!all(methods %in% unique(functional_alpha$method)) ||
      !all(targets %in% unique(functional_alpha$target))) {
    stop("Functional-testing output is missing a reviewer-facing method or target.")
  }
  functional_call_diagnostics <- rbind(
    extract_functional_call_diagnostics(
      functional_result = raw_results,
      method = "FASH-IWP1-Raw",
      unit_info = out$unit_info,
      true_functionals = true_functionals,
      alpha = 0.05
    ),
    extract_functional_call_diagnostics(
      functional_result = bf_results,
      method = "FASH-IWP1-BF",
      unit_info = out$unit_info,
      true_functionals = true_functionals,
      alpha = 0.05
    )
  )
  functional_call_diagnostics$seed <- seed

  truth_group_column <- if (identical(targeted_profile, "mixed")) {
    "truth_stratum"
  } else {
    "truth_group"
  }
  truth_group_counts <- as.data.frame(table(out$unit_info[[truth_group_column]][dynamic_truth]))
  colnames(truth_group_counts) <- c("truth_stratum", "n_dynamic")
  truth_group_counts$seed <- seed
  truth_group_counts <- truth_group_counts[, c("seed", "truth_stratum", "n_dynamic")]
  spiky_dynamic <- dynamic_truth & out$unit_info$shape_profile == "spiky"
  spike_pattern_counts <- if (any(spiky_dynamic)) {
    counts <- as.data.frame(table(
      time_group = out$unit_info$time_group[spiky_dynamic],
      switch_status = out$unit_info$switch_status[spiky_dynamic],
      spike_pattern = out$unit_info$spike_pattern[spiky_dynamic]
    ))
    counts <- counts[counts$Freq > 0, ]
    names(counts)[names(counts) == "Freq"] <- "n_dynamic"
    counts$seed <- seed
    counts[, c(
      "seed", "time_group", "switch_status", "spike_pattern", "n_dynamic"
    )]
  } else {
    data.frame(
      seed = integer(),
      time_group = character(),
      switch_status = character(),
      spike_pattern = character(),
      n_dynamic = integer(),
      stringsAsFactors = FALSE
    )
  }

  if (cache_full_fits) {
    full_fit_path <- file.path(full_fit_dir, paste0("seed_", seed, ".rds"))
    saveRDS(out, full_fit_path)
    message("Cached full fitted output: ", full_fit_path)
  }

  list(
    configuration = configuration,
    seed = seed,
    component_seeds = component_seeds,
    true_pi0 = true_pi0,
    functional_alpha = functional_alpha,
    functional_alpha_005 = functional_alpha[abs(functional_alpha$alpha - 0.05) < 1e-12, ],
    functional_call_diagnostics = functional_call_diagnostics,
    truth_group_counts = truth_group_counts,
    spike_pattern_counts = spike_pattern_counts,
    example_curves = extract_truth_examples(out, effect_sim),
    estimated_pi0 = data.frame(
      seed = seed,
      method = c("FASH-IWP1", "FASH-IWP1"),
      fit = c("Raw", "BF-corrected"),
      estimated_pi0 = c(
        constant_component_prior_weight(out$fash_fits$fash_iwp1_raw),
        constant_component_prior_weight(out$fash_fits$fash_iwp1_bf)
      ),
      stringsAsFactors = FALSE
    ),
    global_alpha = if (is.null(global_results)) NULL else global_results$alpha_curve,
    global_alpha_005 = if (is.null(global_results)) NULL else global_results$alpha_005,
    global_shape_power = if (is.null(global_results)) NULL else global_results$shape_power,
    global_pi0 = if (is.null(global_results)) NULL else global_results$pi0
  )
}

validate_replicate <- function(replicate, seed) {
  required_fields <- c(
    "configuration", "seed", "component_seeds", "true_pi0", "functional_alpha",
    "functional_alpha_005", "functional_call_diagnostics", "truth_group_counts",
    "spike_pattern_counts",
    "example_curves", "estimated_pi0"
  )
  if (!all(required_fields %in% names(replicate)) ||
      !identical(replicate$seed, seed) ||
      !configuration_matches(replicate$configuration) ||
      !isTRUE(all.equal(replicate$true_pi0, true_pi0_expected))) {
    return(FALSE)
  }
  expected_alpha_rows <- length(methods) * length(targets) * length(alpha_grid)
  expected_truth_strata <- if (identical(targeted_profile, "mixed")) 12L else 6L
  functional_valid <- nrow(replicate$functional_alpha) == expected_alpha_rows &&
    nrow(replicate$functional_alpha_005) == length(methods) * length(targets) &&
    nrow(replicate$truth_group_counts) == expected_truth_strata &&
    nrow(replicate$estimated_pi0) == 2
  if (spiky_truth_version == "mixed_single_double_v2") {
    functional_valid <- functional_valid &&
      nrow(replicate$spike_pattern_counts) == 9L &&
      all(c("single", "same-sign double", "opposite-sign double") %in%
        replicate$spike_pattern_counts$spike_pattern)
  }
  if (!functional_valid) {
    return(FALSE)
  }
  if (!include_global_comparisons) {
    return(TRUE)
  }
  global_fields <- c("global_alpha", "global_alpha_005", "global_pi0")
  if (!all(global_fields %in% names(replicate)) ||
      any(vapply(replicate[global_fields], is.null, logical(1)))) {
    return(FALSE)
  }
  global_valid <- all(global_methods %in% unique(replicate$global_alpha$method)) &&
    nrow(replicate$global_alpha_005) == length(global_methods) &&
    nrow(replicate$global_pi0) == 4
  if (!global_valid) {
    return(FALSE)
  }
  if (identical(targeted_profile, "mixed")) {
    expected_shape_rows <- length(global_methods) * 2L * length(alpha_grid)
    return(
      !is.null(replicate$global_shape_power) &&
        nrow(replicate$global_shape_power) == expected_shape_rows
    )
  }
  TRUE
}

replicates <- lapply(seed_list, function(seed) {
  replicate_path <- file.path(replicate_dir, paste0("seed_", seed, ".rds"))
  full_fit_path <- file.path(full_fit_dir, paste0("seed_", seed, ".rds"))
  if (file.exists(replicate_path) && !overwrite) {
    cached <- readRDS(replicate_path)
    cache_is_complete <- !cache_full_fits || file.exists(full_fit_path)
    if (validate_replicate(cached, seed) && cache_is_complete) {
      message("Reusing compact replicate cache: ", replicate_path)
      return(cached)
    }
    if (validate_replicate(cached, seed) && !cache_is_complete) {
      message(
        "Compact cache exists but full fitted output is missing; rerunning seed ",
        seed,
        "."
      )
    } else {
      stop("Cached replicate does not match the requested settings: ", replicate_path)
    }
  }
  message("Running functional-testing ", truth_mechanism, " replicate with seed ", seed, ".")
  replicate <- make_replicate(seed)
  saveRDS(replicate, replicate_path)
  replicate
})

all_alpha <- do.call(rbind, lapply(replicates, `[[`, "functional_alpha"))
all_alpha_005 <- do.call(rbind, lapply(replicates, `[[`, "functional_alpha_005"))
all_functional_calls <- do.call(
  rbind,
  lapply(replicates, `[[`, "functional_call_diagnostics")
)
all_truth_groups <- do.call(rbind, lapply(replicates, `[[`, "truth_group_counts"))
all_spike_patterns <- do.call(rbind, lapply(replicates, `[[`, "spike_pattern_counts"))
all_pi0 <- do.call(rbind, lapply(replicates, `[[`, "estimated_pi0"))
mc_alpha <- summarize_mc_functional_alpha_curves(all_alpha)
mc_alpha_005 <- mc_alpha[abs(mc_alpha$alpha - 0.05) < 1e-12, ]
mc_pi0 <- summarize_mc_pi0(all_pi0)
example_curves <- replicates[[1]]$example_curves

if (any(mc_alpha$n_replications != length(seed_list)) ||
    any(mc_alpha$power_ci_lower < 0, na.rm = TRUE) ||
    any(mc_alpha$power_ci_upper > 1, na.rm = TRUE) ||
    any(mc_alpha$estimated_fsr_ci_lower < 0, na.rm = TRUE) ||
    any(mc_alpha$estimated_fsr_ci_upper > 1, na.rm = TRUE) ||
    any(mc_alpha$empirical_fsr_ci_lower < 0, na.rm = TRUE) ||
    any(mc_alpha$empirical_fsr_ci_upper > 1, na.rm = TRUE) ||
    any(mc_alpha$conditional_empirical_fsr_ci_lower < 0, na.rm = TRUE) ||
    any(mc_alpha$conditional_empirical_fsr_ci_upper > 1, na.rm = TRUE)) {
  stop("Functional-testing Monte Carlo summaries failed validation.")
}

write_csv(all_alpha, file.path(summary_dir, "all_replicate_functional_alpha_curves.csv"))
write_csv(all_alpha_005, file.path(summary_dir, "all_replicate_functional_alpha005.csv"))
write_csv(
  all_functional_calls,
  file.path(summary_dir, "all_replicate_functional_call_diagnostics_alpha005.csv")
)
write_csv(all_truth_groups, file.path(summary_dir, "all_replicate_truth_group_counts.csv"))
write_csv(
  all_spike_patterns,
  file.path(summary_dir, "all_replicate_spike_pattern_counts.csv")
)
write_csv(all_pi0, file.path(summary_dir, "all_replicate_pi0.csv"))
write_csv(mc_alpha, file.path(summary_dir, "functional_testing_mc_alpha_curve.csv"))
write_csv(mc_alpha_005, file.path(summary_dir, "functional_testing_mc_alpha005_summary.csv"))
write_csv(mc_pi0, file.path(summary_dir, "functional_testing_mc_pi0_summary.csv"))
saveRDS(example_curves, file.path(output_dir, "example_curves.rds"))

if (include_global_comparisons) {
  all_global_alpha <- do.call(rbind, lapply(replicates, `[[`, "global_alpha"))
  all_global_alpha_005 <- do.call(rbind, lapply(replicates, `[[`, "global_alpha_005"))
  all_global_pi0 <- do.call(rbind, lapply(replicates, `[[`, "global_pi0"))
  all_global_shape_power <- if (identical(targeted_profile, "mixed")) {
    do.call(rbind, lapply(replicates, `[[`, "global_shape_power"))
  } else {
    NULL
  }
  mc_global_alpha <- summarize_mc_alpha_curves(all_global_alpha)
  mc_global_alpha_005 <- mc_global_alpha[
    abs(mc_global_alpha$alpha - 0.05) < 1e-12,
    ,
    drop = FALSE
  ]
  mc_global_pi0 <- summarize_mc_pi0(all_global_pi0)

  if (any(mc_global_alpha$n_replications != length(seed_list)) ||
      any(mc_global_alpha$power_ci_lower < 0, na.rm = TRUE) ||
      any(mc_global_alpha$power_ci_upper > 1, na.rm = TRUE) ||
      any(mc_global_alpha$fdr_ci_lower < 0, na.rm = TRUE) ||
      any(mc_global_alpha$fdr_ci_upper > 1, na.rm = TRUE)) {
    stop("Global genotype-level Monte Carlo summaries failed validation.")
  }

  mc_global_fash_alpha <- mc_global_alpha[
    mc_global_alpha$method %in% global_fash_methods,
    ,
    drop = FALSE
  ]
  mc_global_direct_alpha <- mc_global_alpha[
    mc_global_alpha$method %in% global_direct_methods,
    ,
    drop = FALSE
  ]
  mc_global_fash_alpha_005 <- mc_global_alpha_005[
    mc_global_alpha_005$method %in% global_fash_methods,
    ,
    drop = FALSE
  ]
  mc_global_direct_alpha_005 <- mc_global_alpha_005[
    mc_global_alpha_005$method %in% global_direct_methods,
    ,
    drop = FALSE
  ]

  write_csv(all_global_alpha, file.path(summary_dir, "all_replicate_global_alpha_curves.csv"))
  write_csv(all_global_alpha_005, file.path(summary_dir, "all_replicate_global_alpha005.csv"))
  write_csv(all_global_pi0, file.path(summary_dir, "all_replicate_global_pi0.csv"))
  write_csv(mc_global_alpha, file.path(summary_dir, "global_mc_alpha_curve.csv"))
  write_csv(mc_global_alpha_005, file.path(summary_dir, "global_mc_alpha005_summary.csv"))
  write_csv(mc_global_pi0, file.path(summary_dir, "global_mc_pi0_summary.csv"))
  if (!is.null(all_global_shape_power)) {
    mc_global_shape_power <- summarize_mc_shape_power(all_global_shape_power)
    if (any(mc_global_shape_power$n_replications != length(seed_list)) ||
        any(mc_global_shape_power$power_ci_lower < 0, na.rm = TRUE) ||
        any(mc_global_shape_power$power_ci_upper > 1, na.rm = TRUE)) {
      stop("Shape-stratified Monte Carlo power summaries failed validation.")
    }
    write_csv(
      all_global_shape_power,
      file.path(summary_dir, "all_replicate_global_shape_power.csv")
    )
    write_csv(
      mc_global_shape_power,
      file.path(summary_dir, "global_mc_shape_power_curve.csv")
    )
    write_csv(
      mc_global_shape_power[abs(mc_global_shape_power$alpha - 0.05) < 1e-12, ],
      file.path(summary_dir, "global_mc_shape_power_alpha005_summary.csv")
    )
  }
  write_csv(
    mc_global_fash_alpha,
    file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha_curve.csv")
  )
  write_csv(
    mc_global_fash_alpha_005,
    file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha005_summary.csv")
  )
  write_csv(
    mc_global_direct_alpha,
    file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_mc_alpha_curve.csv")
  )
  write_csv(
    mc_global_direct_alpha_005,
    file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_mc_alpha005_summary.csv")
  )

  global_plot_subtitle <- paste0(
    length(seed_list), " seeds; ",
    if (identical(targeted_profile, "mixed")) {
      "balanced broad and spiky local cubic B-spline effects"
    } else {
      "targeted local cubic B-spline effects"
    },
    "; N = ",
    n_donors, ", J = ", J, ", five covariates"
  )
  plot_mc_alpha_curves(
    mc_global_fash_alpha,
    metric = "power",
    file = file.path(figure_dir, "iwp_vs_linear_fash_mc_power.png"),
    title = "IWP versus linear FASH: Monte Carlo power",
    subtitle = global_plot_subtitle,
    style_profile = "combined"
  )
  plot_mc_alpha_curves(
    mc_global_fash_alpha,
    metric = "fdr",
    file = file.path(figure_dir, "iwp_vs_linear_fash_mc_fdr.png"),
    title = "IWP versus linear FASH: Monte Carlo FDR estimate",
    subtitle = global_plot_subtitle,
    legend_position = "topleft",
    style_profile = "combined"
  )
  plot_mc_alpha_curves(
    mc_global_direct_alpha,
    metric = "power",
    file = file.path(figure_dir, "iwp_fash_vs_direct_true_pi0_mc_power.png"),
    title = "IWP FASH versus direct interaction: Monte Carlo power",
    subtitle = paste0(global_plot_subtitle, "; direct eFDR uses true pi0 = ", true_pi0_expected),
    style_profile = "combined"
  )
  plot_mc_alpha_curves(
    mc_global_direct_alpha,
    metric = "fdr",
    file = file.path(figure_dir, "iwp_fash_vs_direct_true_pi0_mc_fdr.png"),
    title = "IWP FASH versus direct interaction: Monte Carlo FDR estimate",
    subtitle = paste0(global_plot_subtitle, "; direct eFDR uses true pi0 = ", true_pi0_expected),
    legend_position = "bottomright",
    style_profile = "combined"
  )
}

truth_title <- if (truth_mechanism == "random_bspline") {
  "labelled random B-spline effects"
} else if (identical(targeted_profile, "mixed")) {
  "balanced broad and spiky local cubic B-spline effects"
} else {
  "targeted local cubic B-spline effects"
}

plot_mc_functional_curve_grid(
  mc_curve = mc_alpha,
  metric = "power",
  file = file.path(figure_dir, "functional_testing_mc_power.png"),
  title = paste("Functional-testing power for", truth_title)
)
plot_mc_functional_curve_grid(
  mc_curve = mc_alpha,
  metric = "estimated_fsr",
  file = file.path(figure_dir, "functional_testing_mc_estimated_fsr.png"),
  title = paste("Posterior estimated FSR for", truth_title)
)
plot_mc_functional_curve_grid(
  mc_curve = mc_alpha,
  metric = "empirical_fsr",
  file = file.path(figure_dir, "functional_testing_mc_empirical_fsr.png"),
  title = paste("End-to-end false-call proportion for", truth_title)
)
plot_mc_functional_curve_grid(
  mc_curve = mc_alpha,
  metric = "conditional_empirical_fsr",
  file = file.path(figure_dir, "functional_testing_mc_conditional_empirical_fsr.png"),
  title = paste("Conditional empirical FSR for", truth_title)
)

message("Saved functional-testing cache to: ", output_dir)
