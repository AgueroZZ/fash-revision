#!/usr/bin/env Rscript

# Run paired functional-testing simulations under broad random B-spline
# and compact raised-cosine truth mechanisms.

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

parse_integer_list <- function(x, name) {
  values <- suppressWarnings(as.integer(
    trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  ))
  if (length(values) == 0 || anyNA(values) || anyDuplicated(values)) {
    stop(name, " must contain unique comma-separated integers.")
  }
  values
}

parse_character_list <- function(x, name) {
  values <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  if (length(values) == 0 || any(!nzchar(values)) || anyDuplicated(values)) {
    stop(name, " must contain unique comma-separated values.")
  }
  values
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "shared",
  "simulation_functions.R"
))

J <- as.integer(get_arg("--J", "1000"))
n_donors <- as.integer(get_arg("--n-donors", "19"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
expression_noise_sd <- as.numeric(get_arg("--noise-sd", "1"))
dynamic_main_effect_sd <- as.numeric(get_arg(
  "--dynamic-main-effect-sd",
  "1"
))
num_basis <- as.integer(get_arg("--num-basis", "20"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
seed_list <- parse_integer_list(
  get_arg("--seed-list", "12345,22345,32345,42345,52345"),
  "--seed-list"
)
truth_mechanisms <- parse_character_list(
  get_arg("--truth-mechanisms", "random_bspline,raised_cosine"),
  "--truth-mechanisms"
)
output_id <- get_arg(
  "--output-id",
  "r3_matched_functional_relative_clearance_main_effect_pilot5"
)
n_examples_per_group <- as.integer(get_arg(
  "--n-examples-per-group",
  "2"
))
overwrite <- as_flag(get_arg("--overwrite", "false"))

allowed_mechanisms <- c("random_bspline", "raised_cosine")
if (!all(truth_mechanisms %in% allowed_mechanisms) ||
    J < 30 ||
    n_donors < n_covariates + 3 ||
    n_covariates < 0 ||
    !is.finite(expression_noise_sd) ||
    expression_noise_sd <= 0 ||
    !is.finite(dynamic_main_effect_sd) ||
    dynamic_main_effect_sd <= 0 ||
    num_basis < 2 ||
    num_cores < 1 ||
    n_examples_per_group < 1 ||
    !nzchar(output_id)) {
  stop("Invalid matched functional-testing arguments.")
}

time_grid <- make_time_grid()
evaluation_grid <- seq(min(time_grid), max(time_grid), by = 0.1)
alpha_grid <- seq(0.005, 0.20, by = 0.005)
class_probs <- c(
  dynamic_bspline = 0.20,
  constant = 0.40,
  zero = 0.40
)
methods <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
targets <- c("early", "middle", "late", "switch")
switch_threshold <- 0.25
location_truth_margin <- 0.10
switch_truth_margin <- 0.10
non_switch_min_abs <- 0.10
non_switch_min_range_fraction <- 0.10
maf_range <- c(0.10, 0.50)
covariate_effect_sd <- 0.50
intercept_sd <- 0

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
  J = J,
  n_donors = n_donors,
  n_covariates = n_covariates,
  time_grid = time_grid,
  evaluation_grid = evaluation_grid,
  expression_noise_sd = expression_noise_sd,
  maf_range = maf_range,
  covariate_effect_sd = covariate_effect_sd,
  intercept_sd = intercept_sd,
  dynamic_main_effect_sd = dynamic_main_effect_sd,
  num_basis = num_basis,
  class_probs = class_probs,
  truth_mechanisms = truth_mechanisms,
  random_bspline = list(
    amplitude = 2,
    df = 6,
    coefficient_sd = 1
  ),
  raised_cosine = list(
    width_half = 1.5,
    spike_counts = 1:3,
    relative_amplitude_range = c(0.35, 0.75),
    target_centered_rms = 0.90
  ),
  switch_threshold = switch_threshold,
  location_truth_margin = location_truth_margin,
  switch_truth_margin = switch_truth_margin,
  non_switch_min_abs = non_switch_min_abs,
  non_switch_min_range_fraction =
    non_switch_min_range_fraction,
  alpha_grid = alpha_grid,
  n_examples_per_group = n_examples_per_group,
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

extract_examples <- function(out, effect_sim, truth_mechanism) {
  dynamic <- out$unit_info$effect_class == "dynamic_bspline"
  group_order <- c(
    "early / switch",
    "early / non-switch",
    "middle / switch",
    "middle / non-switch",
    "late / switch",
    "late / non-switch"
  )
  examples <- lapply(group_order, function(group) {
    candidates <- which(dynamic & out$unit_info$truth_group == group)
    if (truth_mechanism == "raised_cosine") {
      spike_levels <- sort(unique(
        out$unit_info$spike_count[candidates]
      ))
      selected_levels <- spike_levels[unique(round(seq(
        1,
        length(spike_levels),
        length.out = min(n_examples_per_group, length(spike_levels))
      )))]
      diverse <- vapply(
        selected_levels,
        function(level) {
          candidates[
            which(out$unit_info$spike_count[candidates] == level)[1]
          ]
        },
        integer(1)
      )
      candidates <- c(diverse, setdiff(candidates, diverse))
    }
    selected <- head(candidates, n_examples_per_group)
    if (length(selected) != n_examples_per_group) {
      stop("Could not extract enough examples for ", group, ".")
    }
    lapply(selected, function(index) {
      list(
        mechanism = truth_mechanism,
        truth_group = group,
        time_group = out$unit_info$time_group[index],
        switch_status = out$unit_info$switch_status[index],
        spike_count = out$unit_info$spike_count[index],
        variant_id = out$unit_info$variant_id[index],
        genetic_main_effect =
          out$unit_info$genetic_main_effect[index],
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

extract_call_diagnostics <- function(functional_result,
                                     method,
                                     out,
                                     true_functionals,
                                     truth_mechanism,
                                     seed,
                                     alpha = 0.05) {
  dynamic_indices <- sort(unique(
    functional_result$fdr_table$index[
      functional_result$fdr_table$FDR <= alpha
    ]
  ))
  rows <- lapply(targets, function(target) {
    lfsr_map <- functional_result$lfsr_by_target[[target]]
    candidate_lfsr <- if (length(dynamic_indices) == 0) {
      numeric()
    } else {
      unname(lfsr_map[as.character(dynamic_indices)])
    }
    cfsr_table <- functional_cfsr_table(
      dynamic_indices,
      candidate_lfsr
    )
    selected <- cfsr_table[cfsr_table$cfsr <= alpha, , drop = FALSE]
    if (nrow(selected) == 0) return(NULL)
    indices <- selected$index
    data.frame(
      seed = seed,
      truth_mechanism = truth_mechanism,
      method = method,
      target = target,
      alpha = alpha,
      unit_index = indices,
      variant_id = out$unit_info$variant_id[indices],
      effect_class = out$unit_info$effect_class[indices],
      truth_group = out$unit_info$truth_group[indices],
      genetic_main_effect =
        out$unit_info$genetic_main_effect[indices],
      true_functional = true_functionals[indices, target],
      false_discovery = true_functionals[indices, target] <= 0,
      lfsr = selected$lfsr,
      cfsr = selected$cfsr,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    return(data.frame(
      seed = integer(),
      truth_mechanism = character(),
      method = character(),
      target = character(),
      alpha = numeric(),
      unit_index = integer(),
      variant_id = character(),
      effect_class = character(),
      truth_group = character(),
      genetic_main_effect = numeric(),
      true_functional = numeric(),
      false_discovery = logical(),
      lfsr = numeric(),
      cfsr = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  rownames(out) <- NULL
  out
}

make_replicate <- function(seed, truth_mechanism) {
  component_seeds <- revision_component_seeds(seed)
  scenario <- paste0(
    if (truth_mechanism == "random_bspline") "r3a_" else "r3b_",
    "matched_functional_",
    truth_mechanism,
    "_relative_clearance_main_effect"
  )
  genotype_sim <- simulate_genotype_matrix(
    n_donors = n_donors,
    n_variants = J,
    maf_range = maf_range,
    seed = component_seeds[["genotype"]]
  )
  covariates <- simulate_covariate_matrix(
    n_donors = n_donors,
    n_covariates = n_covariates,
    seed = component_seeds[["covariates"]]
  )
  effect_sim <- simulate_matched_functional_effect_set(
    n_variants = J,
    truth_mechanism = truth_mechanism,
    time_grid = time_grid,
    evaluation_grid = evaluation_grid,
    class_probs = class_probs,
    dynamic_main_effect_sd = dynamic_main_effect_sd,
    switch_threshold = switch_threshold,
    location_truth_margin = location_truth_margin,
    switch_truth_margin = switch_truth_margin,
    non_switch_min_abs = non_switch_min_abs,
    non_switch_min_range_fraction =
      non_switch_min_range_fraction,
    seed = seed,
    class_seed = component_seeds[["classes"]],
    constant_seed = component_seeds[["constant_effects"]],
    shape_seed = component_seeds[["functional_truth"]],
    scenario = scenario
  )
  expression_sim <- simulate_eqtl_expression_from_genotypes(
    G = genotype_sim$G,
    beta_matrix = effect_sim$beta_matrix,
    time_grid = time_grid,
    covariates = covariates,
    expression_noise_sd = expression_noise_sd,
    covariate_effect_sd = covariate_effect_sd,
    intercept_sd = intercept_sd,
    seed = component_seeds[["expression"]]
  )
  out <- run_genotype_level_dynamic_eqtl_simulation(
    G = genotype_sim$G,
    time_grid = time_grid,
    covariates = covariates,
    class_probs = class_probs,
    maf_range = maf_range,
    expression_noise_sd = expression_noise_sd,
    covariate_effect_sd = covariate_effect_sd,
    intercept_sd = intercept_sd,
    dynamic_main_effect_sd = dynamic_main_effect_sd,
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

  true_functionals <- effect_sim$true_functionals[
    out$unit_info$variant_id,
    targets,
    drop = FALSE
  ]
  true_dynamic <- out$unit_info$effect_class == "dynamic_bspline"
  raw_results <- evaluate_fash_functional_testing(
    fit = out$fash_fits$fash_iwp1_raw,
    true_functionals = true_functionals,
    evaluation_grid = evaluation_grid,
    alpha_grid = alpha_grid,
    method = "FASH-IWP1-Raw",
    scenario = scenario,
    switch_threshold = switch_threshold,
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
    switch_threshold = switch_threshold,
    true_dynamic = true_dynamic,
    num_cores = num_cores,
    seed = component_seeds[["functional_posterior"]] + 100L
  )
  functional_alpha <- rbind(
    raw_results$alpha_curve,
    bf_results$alpha_curve
  )
  functional_alpha$seed <- seed
  functional_alpha$truth_mechanism <- truth_mechanism

  call_diagnostics <- rbind(
    extract_call_diagnostics(
      raw_results,
      "FASH-IWP1-Raw",
      out,
      true_functionals,
      truth_mechanism,
      seed
    ),
    extract_call_diagnostics(
      bf_results,
      "FASH-IWP1-BF",
      out,
      true_functionals,
      truth_mechanism,
      seed
    )
  )
  truth_group_counts <- as.data.frame(table(
    out$unit_info$truth_group[true_dynamic]
  ))
  names(truth_group_counts) <- c("truth_group", "n_dynamic")
  truth_group_counts$seed <- seed
  truth_group_counts$truth_mechanism <- truth_mechanism

  list(
    configuration = configuration,
    seed = seed,
    truth_mechanism = truth_mechanism,
    component_seeds = component_seeds,
    functional_alpha = functional_alpha,
    functional_alpha_005 = functional_alpha[
      abs(functional_alpha$alpha - 0.05) < 1e-12,
      ,
      drop = FALSE
    ],
    call_diagnostics_alpha005 = call_diagnostics,
    truth_group_counts = truth_group_counts,
    estimated_pi0 = data.frame(
      seed = seed,
      truth_mechanism = truth_mechanism,
      method = c("FASH-IWP1", "FASH-IWP1"),
      fit = c("Raw", "BF-corrected"),
      estimated_pi0 = c(
        constant_component_prior_weight(
          out$fash_fits$fash_iwp1_raw
        ),
        constant_component_prior_weight(
          out$fash_fits$fash_iwp1_bf
        )
      ),
      stringsAsFactors = FALSE
    ),
    example_curves = extract_examples(
      out,
      effect_sim,
      truth_mechanism
    )
  )
}

validate_replicate <- function(x, seed, truth_mechanism) {
  required <- c(
    "configuration",
    "seed",
    "truth_mechanism",
    "component_seeds",
    "functional_alpha",
    "functional_alpha_005",
    "call_diagnostics_alpha005",
    "truth_group_counts",
    "estimated_pi0",
    "example_curves"
  )
  expected_rows <- length(methods) * length(targets) *
    length(alpha_grid)
  all(required %in% names(x)) &&
    identical(x$seed, seed) &&
    identical(x$truth_mechanism, truth_mechanism) &&
    isTRUE(all.equal(x$configuration, configuration)) &&
    nrow(x$functional_alpha) == expected_rows &&
    nrow(x$functional_alpha_005) == length(methods) * length(targets) &&
    all(methods %in% unique(x$functional_alpha$method)) &&
    all(targets %in% unique(x$functional_alpha$target)) &&
    nrow(x$truth_group_counts) == 6
}

replicates <- list()
replicate_index <- 1L
for (truth_mechanism in truth_mechanisms) {
  for (seed in seed_list) {
    path <- file.path(
      replicate_dir,
      paste0(truth_mechanism, "_seed_", seed, ".rds")
    )
    if (file.exists(path) && !overwrite) {
      replicate <- readRDS(path)
      if (!validate_replicate(replicate, seed, truth_mechanism)) {
        stop("Cached replicate does not match: ", path)
      }
      message("Reusing matched functional replicate: ", path)
    } else {
      message(
        "Running ",
        truth_mechanism,
        " functional replicate with seed ",
        seed,
        "."
      )
      replicate <- make_replicate(seed, truth_mechanism)
      saveRDS(replicate, path)
    }
    replicates[[replicate_index]] <- replicate
    replicate_index <- replicate_index + 1L
  }
}

all_alpha <- do.call(rbind, lapply(
  replicates,
  `[[`,
  "functional_alpha"
))
all_alpha_005 <- do.call(rbind, lapply(
  replicates,
  `[[`,
  "functional_alpha_005"
))
all_calls <- do.call(rbind, lapply(
  replicates,
  `[[`,
  "call_diagnostics_alpha005"
))
all_truth_groups <- do.call(rbind, lapply(
  replicates,
  `[[`,
  "truth_group_counts"
))
all_pi0 <- do.call(rbind, lapply(replicates, `[[`, "estimated_pi0"))
mc_alpha <- summarize_mc_functional_alpha_curves(all_alpha)
mc_alpha_005 <- mc_alpha[
  abs(mc_alpha$alpha - 0.05) < 1e-12,
  ,
  drop = FALSE
]
mc_pi0 <- do.call(rbind, lapply(
  names(split(all_pi0, all_pi0$truth_mechanism)),
  function(mechanism) {
    summary <- summarize_mc_pi0(
      all_pi0[all_pi0$truth_mechanism == mechanism, , drop = FALSE]
    )
    summary$truth_mechanism <- mechanism
    summary
  }
))
example_curves <- lapply(truth_mechanisms, function(mechanism) {
  matching <- which(vapply(
    replicates,
    function(x) identical(x$truth_mechanism, mechanism),
    logical(1)
  ))
  replicates[[matching[1]]]$example_curves
})
names(example_curves) <- truth_mechanisms

write_csv(
  all_alpha,
  file.path(summary_dir, "all_replicate_functional_alpha_curves.csv")
)
write_csv(
  all_alpha_005,
  file.path(summary_dir, "all_replicate_functional_alpha005.csv")
)
write_csv(
  all_calls,
  file.path(summary_dir, "all_functional_calls_alpha005.csv")
)
write_csv(
  all_truth_groups,
  file.path(summary_dir, "all_truth_group_counts.csv")
)
write_csv(
  all_pi0,
  file.path(summary_dir, "all_replicate_pi0.csv")
)
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
saveRDS(
  example_curves,
  file.path(output_dir, "example_curves.rds")
)

for (truth_mechanism in truth_mechanisms) {
  scenario_prefix <- if (truth_mechanism == "random_bspline") {
    "r3a_"
  } else {
    "r3b_"
  }
  mechanism_curve <- mc_alpha[
    grepl(paste0("^", scenario_prefix), mc_alpha$scenario),
    ,
    drop = FALSE
  ]
  truth_title <- if (truth_mechanism == "random_bspline") {
    "broad random B-spline effects"
  } else {
    "compact raised-cosine effects"
  }
  plot_mc_functional_curve_grid(
    mc_curve = mechanism_curve,
    metric = "power",
    file = file.path(
      figure_dir,
      paste0(truth_mechanism, "_functional_power.png")
    ),
    title = paste("Functional-testing power:", truth_title)
  )
  plot_mc_functional_curve_grid(
    mc_curve = mechanism_curve,
    metric = "empirical_fsr",
    file = file.path(
      figure_dir,
      paste0(truth_mechanism, "_functional_empirical_fsr.png")
    ),
    title = paste("Empirical functional FSR:", truth_title)
  )
}

print(mc_alpha_005)
print(mc_pi0)
message("Saved matched functional-testing cache to: ", output_dir)
