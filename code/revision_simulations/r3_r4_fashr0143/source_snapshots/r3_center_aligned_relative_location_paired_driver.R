#!/usr/bin/env Rscript

# Run the center-aligned, relative-clearance, paired-posterior R3 simulation
# under broad random B-spline and compact raised-cosine truth mechanisms. Each
# seed uses one paper-tested YRI dosage variant per gene from the shared cache.

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

parse_numeric_pair <- function(x, name) {
  values <- suppressWarnings(as.numeric(
    trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  ))
  if (length(values) != 2L ||
      any(!is.finite(values)) ||
      values[[1L]] >= values[[2L]]) {
    stop(name, " must contain two increasing finite values.")
  }
  values
}

parse_temporal_category_probs <- function(x, name) {
  if (!nzchar(x)) return(NULL)
  values <- suppressWarnings(as.numeric(
    trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  ))
  if (length(values) != 3L ||
      any(!is.finite(values)) ||
      any(values <= 0) ||
      abs(sum(values) - 1) > 1e-8) {
    stop(paste(
      name,
      "must contain positive Early,Middle,Late probabilities that sum to one."
    ))
  }
  stats::setNames(values, c("early", "middle", "late"))
}

parse_character_list <- function(x, name) {
  values <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  if (length(values) == 0 || any(!nzchar(values)) || anyDuplicated(values)) {
    stop(name, " must contain unique comma-separated values.")
  }
  values
}

parse_md5_list <- function(x, seed_list, name) {
  if (!nzchar(x)) return(NULL)
  values <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  if (length(values) != length(seed_list) ||
      any(!grepl("^[[:xdigit:]]{32}$", values))) {
    stop(
      name,
      " must contain one comma-separated 32-character MD5 per configured seed."
    )
  }
  stats::setNames(tolower(values), as.character(seed_list))
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

workflowr_root <- find_workflowr_root()
simulation_functions_path <- resolve_input_path(
  Sys.getenv(
    "FASH_R3_SIMULATION_FUNCTIONS",
    unset = file.path(
      workflowr_root,
      "code", "revision_simulations", "shared", "simulation_functions.R"
    )
  ),
  workflowr_root,
  "FASH_R3_SIMULATION_FUNCTIONS"
)
real_genotype_helper_path <- resolve_input_path(
  Sys.getenv(
    "FASH_R3_REAL_GENOTYPE_HELPER",
    unset = file.path(
      workflowr_root,
      "code", "revision_simulations", "shared",
      "real_genotype_one_per_gene.R"
    )
  ),
  workflowr_root,
  "FASH_R3_REAL_GENOTYPE_HELPER"
)
source(simulation_functions_path)
source(real_genotype_helper_path)

J <- as.integer(get_arg("--J", "6362"))
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
run_seed_list <- parse_integer_list(
  get_arg("--run-seed-list", paste(seed_list, collapse = ",")),
  "--run-seed-list"
)
run_truth_mechanisms <- parse_character_list(
  get_arg(
    "--run-truth-mechanisms",
    paste(truth_mechanisms, collapse = ",")
  ),
  "--run-truth-mechanisms"
)
output_id <- get_arg(
  "--output-id",
  paste0(
    "r3_real_genotype_one_per_gene_J6362_",
    "matched_functional_relative_clearance_main_effect_pilot5"
  )
)
output_dir_argument <- get_arg("--output-dir", "")
expected_fashr_version <- get_arg("--expected-fashr-version", "")
expected_fashr_remote_sha <- get_arg("--expected-fashr-remote-sha", "")
expected_genotype_content_md5 <- parse_md5_list(
  get_arg("--expected-genotype-content-md5", ""),
  seed_list,
  "--expected-genotype-content-md5"
)
middle_window <- parse_numeric_pair(
  get_arg("--middle-window", "4,11"),
  "--middle-window"
)
middle_boundary <- match.arg(
  get_arg("--middle-boundary", "closed"),
  c("closed", "open")
)
temporal_category_probs <- parse_temporal_category_probs(
  get_arg("--temporal-category-probs", ""),
  "--temporal-category-probs"
)
location_truth_min_range_fraction <- as.numeric(get_arg(
  "--location-truth-min-range-fraction",
  "0"
))
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
n_examples_per_group <- as.integer(get_arg(
  "--n-examples-per-group",
  "2"
))
replicate_only <- as_flag(get_arg("--replicate-only", "false"))
preflight_only <- as_flag(get_arg("--preflight-only", "false"))
overwrite <- as_flag(get_arg("--overwrite", "false"))

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required for the R3 simulation.")
}
fashr_description <- utils::packageDescription("fashr")
package_provenance <- list(
  package = "fashr",
  version = as.character(utils::packageVersion("fashr")),
  remote_sha = if (is.null(fashr_description$RemoteSha)) {
    NA_character_
  } else {
    as.character(fashr_description$RemoteSha)
  },
  library_path = normalizePath(
    find.package("fashr"),
    winslash = "/",
    mustWork = TRUE
  ),
  r_version = R.version.string,
  platform = R.version$platform
)
if (nzchar(expected_fashr_version) &&
    !identical(package_provenance$version, expected_fashr_version)) {
  stop(
    "Expected fashr ", expected_fashr_version,
    "; found ", package_provenance$version, "."
  )
}
if (nzchar(expected_fashr_remote_sha) &&
    !identical(package_provenance$remote_sha, expected_fashr_remote_sha)) {
  stop(
    "Expected fashr RemoteSha ", expected_fashr_remote_sha,
    "; found ", package_provenance$remote_sha, "."
  )
}

allowed_mechanisms <- c("random_bspline", "raised_cosine")
if (!all(truth_mechanisms %in% allowed_mechanisms) ||
    !all(run_truth_mechanisms %in% truth_mechanisms) ||
    !all(run_seed_list %in% seed_list) ||
    J < 30 ||
    n_donors < n_covariates + 3 ||
    n_covariates < 0 ||
    !is.finite(expression_noise_sd) ||
    expression_noise_sd <= 0 ||
    !is.finite(dynamic_main_effect_sd) ||
    dynamic_main_effect_sd <= 0 ||
    num_basis < 2 ||
    num_cores < 1 ||
    !is.finite(location_truth_min_range_fraction) ||
    location_truth_min_range_fraction < 0 ||
    location_truth_min_range_fraction >= 1 ||
    n_examples_per_group < 1 ||
    !nzchar(output_id)) {
  stop("Invalid matched functional-testing arguments.")
}

time_grid <- make_time_grid()
evaluation_grid <- seq(min(time_grid), max(time_grid), by = 0.1)
middle_membership <- temporal_middle_membership(
  smooth_var = evaluation_grid,
  middle_window = middle_window,
  middle_boundary = middle_boundary
)
middle_grid <- evaluation_grid[middle_membership]
middle_expression <- if (identical(middle_boundary, "open")) {
  sprintf("%g < t < %g", middle_window[[1L]], middle_window[[2L]])
} else {
  sprintf("%g <= t <= %g", middle_window[[1L]], middle_window[[2L]])
}
format_scenario_endpoint <- function(x) {
  gsub("[.]", "p", format(x, trim = TRUE, scientific = FALSE))
}
middle_scenario_label <- paste(
  middle_boundary,
  "middle",
  format_scenario_endpoint(middle_window[[1L]]),
  format_scenario_endpoint(middle_window[[2L]]),
  sep = "_"
)

genotype_cache <- readRDS(genotype_cache_path)
if (!is.list(genotype_cache) ||
    !all(c("configuration", "sample_ids", "samples") %in% names(genotype_cache)) ||
    !identical(genotype_cache$configuration$n_genes, J) ||
    !identical(genotype_cache$configuration$n_donors, n_donors) ||
    !all(seed_list %in% genotype_cache$configuration$seed_list) ||
    !all(as.character(seed_list) %in% names(genotype_cache$samples))) {
  stop("The shared real-genotype cache does not match J, donors, or seeds.")
}
genotype_content_digests <- stats::setNames(
  character(length(seed_list)),
  as.character(seed_list)
)
for (seed in seed_list) {
  sample <- validate_real_genotype_sample(
    genotype_cache$samples[[as.character(seed)]],
    expected_genes = J,
    expected_donors = n_donors,
    maf_min = genotype_cache$configuration$maf_min
  )
  observed_digest <- genotype_content_md5(
    pair_key = sample$selection$pair_key,
    sample_ids = rownames(sample$G),
    G = sample$G
  )
  if (!is.null(expected_genotype_content_md5) && !identical(
    observed_digest,
    expected_genotype_content_md5[[as.character(seed)]]
  )) {
    stop(
      "Unexpected canonical genotype-content digest for seed ", seed,
      ": expected ",
      expected_genotype_content_md5[[as.character(seed)]],
      "; found ", observed_digest, "."
    )
  }
  genotype_content_digests[[as.character(seed)]] <- observed_digest
}
genotype_cache_fingerprint <- artifact_fingerprint(genotype_cache_path)

if (preflight_only) {
  message(
    "Validated the R3 package, arguments, real-genotype cache, and Middle definition: ",
    middle_expression,
    "."
  )
  quit(save = "no", status = 0L)
}

alpha_grid <- seq(0.005, 0.20, by = 0.005)
class_probs <- c(
  dynamic_bspline = 0.20,
  constant = 0.40,
  zero = 0.40
)
expected_class_counts <- exact_proportional_counts(J, class_probs)
truth_group_levels <- c(
  "early / switch",
  "early / non-switch",
  "middle / switch",
  "middle / non-switch",
  "late / switch",
  "late / non-switch"
)
expected_truth_group_counts <- exact_temporal_truth_group_counts(
  unname(expected_class_counts[["dynamic_bspline"]]),
  temporal_category_probs = temporal_category_probs
)
expected_truth_group_counts <- expected_truth_group_counts[truth_group_levels]
expected_truth_group_counts <- setNames(
  as.integer(expected_truth_group_counts),
  names(expected_truth_group_counts)
)
methods <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
targets <- c("early", "middle", "late", "switch")
switch_threshold <- 0.25
location_truth_margin <- 0.10
switch_truth_margin <- 0.10
non_switch_min_abs <- 0.10
non_switch_min_range_fraction <- 0.10
covariate_effect_sd <- 0.50
intercept_sd <- 0

output_dir <- if (nzchar(output_dir_argument)) {
  output_dir_argument
} else {
  file.path(
    workflowr_root,
    "output",
    "revision_simulations",
    "mc",
    output_id
  )
}
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
  covariate_effect_sd = covariate_effect_sd,
  intercept_sd = intercept_sd,
  dynamic_main_effect_sd = dynamic_main_effect_sd,
  num_basis = num_basis,
  class_probs = class_probs,
  expected_class_counts = expected_class_counts,
  expected_truth_group_counts = expected_truth_group_counts,
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
    target_centered_rms = 0.90,
    center_ranges = list(
      early = c(1.5, 2.5),
      middle = c(4.5, 10.5),
      late = c(12.5, 13.5)
    )
  ),
  middle_window = middle_window,
  middle_boundary = middle_boundary,
  middle_grid = middle_grid,
  middle_expression = middle_expression,
  switch_threshold = switch_threshold,
  location_truth_margin = location_truth_margin,
  location_truth_min_range_fraction =
    location_truth_min_range_fraction,
  switch_truth_margin = switch_truth_margin,
  non_switch_min_abs = non_switch_min_abs,
  non_switch_min_range_fraction =
    non_switch_min_range_fraction,
  functional_posterior_pairing = "common_random_seed_raw_bf",
  alpha_grid = alpha_grid,
  n_examples_per_group = n_examples_per_group,
  genotype_source = "paper-derived YRI DS dosage",
  genotype_selection_rule = genotype_cache$configuration$selection_rule,
  genotype_repeated_variant_rule =
    genotype_cache$configuration$repeated_variant_rule,
  genotype_dosage_field = genotype_cache$configuration$dosage_field,
  genotype_maf_min = genotype_cache$configuration$maf_min,
  genotype_sample_ids = genotype_cache$sample_ids,
  genotype_digest_method = "fash-genotype-content-md5-v1",
  genotype_content_digests = genotype_content_digests,
  genotype_cache_fingerprint = genotype_cache_fingerprint,
  genotype_source_configuration = genotype_cache$configuration,
  maf_truth_balance_method = paste(
    "exact global class counts with within-MAF-decile permutation"
  ),
  seed_list = seed_list,
  package_provenance = package_provenance
)
if (!is.null(temporal_category_probs)) {
  equal_temporal_probabilities <- stats::setNames(
    rep(1 / 3, 3),
    c("early", "middle", "late")
  )
  configuration$temporal_category_design <- if (isTRUE(all.equal(
    temporal_category_probs,
    equal_temporal_probabilities,
    tolerance = 1e-8,
    check.attributes = TRUE
  ))) {
    "equal temporal categories"
  } else {
    "user-specified temporal-category probabilities"
  }
  configuration$temporal_category_probs <- temporal_category_probs
}
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
    "real_genotype_one_per_gene_matched_functional_",
    truth_mechanism,
    "_",
    middle_scenario_label,
    "_center_aligned_equal_cells_relative_location_clearance_",
    "paired_posterior_main_effect"
  )
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
  effect_sim <- simulate_matched_functional_effect_set(
    n_variants = J,
    truth_mechanism = truth_mechanism,
    time_grid = time_grid,
    evaluation_grid = evaluation_grid,
    class_probs = class_probs,
    dynamic_main_effect_sd = dynamic_main_effect_sd,
    cosine_center_ranges = configuration$raised_cosine$center_ranges,
    switch_threshold = switch_threshold,
    location_truth_margin = location_truth_margin,
    location_truth_min_range_fraction =
      location_truth_min_range_fraction,
    switch_truth_margin = switch_truth_margin,
    non_switch_min_abs = non_switch_min_abs,
    non_switch_min_range_fraction =
      non_switch_min_range_fraction,
    temporal_category_probs = temporal_category_probs,
    seed = seed,
    class_seed = component_seeds[["classes"]],
    constant_seed = component_seeds[["constant_effects"]],
    shape_seed = component_seeds[["functional_truth"]],
    scenario = scenario,
    middle_window = middle_window,
    middle_boundary = middle_boundary
  )
  if (!identical(
    effect_sim$settings$raised_cosine_center_ranges,
    configuration$raised_cosine$center_ranges
  )) {
    stop("The generated truth does not use the configured center ranges.")
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
    expression_noise_sd = expression_noise_sd,
    covariate_effect_sd = covariate_effect_sd,
    intercept_sd = intercept_sd,
    seed = component_seeds[["expression"]]
  )
  out <- run_genotype_level_dynamic_eqtl_simulation(
    G = genotype_sample$G,
    time_grid = time_grid,
    covariates = covariates,
    class_probs = class_probs,
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
  out$settings$genotype_source <- configuration$genotype_source
  out$settings$genotype_selection_rule <- configuration$genotype_selection_rule
  out$settings$genotype_digest <-
    genotype_content_digests[[as.character(seed)]]
  truth_maf_balance <- summarize_truth_maf_balance(
    variant_info = genotype_sample$variant_info,
    unit_info = out$unit_info,
    seed = seed
  )

  true_functionals <- effect_sim$true_functionals[, targets, drop = FALSE]
  true_dynamic <- out$unit_info$effect_class == "dynamic_bspline"
  functional_posterior_seed <-
    component_seeds[["functional_posterior"]]
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
    seed = functional_posterior_seed,
    middle_window = middle_window,
    middle_boundary = middle_boundary
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
    seed = functional_posterior_seed,
    middle_window = middle_window,
    middle_boundary = middle_boundary
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
    factor(
      out$unit_info$truth_group[true_dynamic],
      levels = truth_group_levels
    )
  ))
  names(truth_group_counts) <- c("truth_group", "n_dynamic")
  truth_group_counts$seed <- seed
  truth_group_counts$truth_mechanism <- truth_mechanism
  if (!identical(
    setNames(truth_group_counts$n_dynamic, truth_group_counts$truth_group),
    expected_truth_group_counts
  )) {
    stop("The dynamic functional truth groups do not match the configured counts.")
  }

  list(
    configuration = configuration,
    seed = seed,
    truth_mechanism = truth_mechanism,
    component_seeds = component_seeds,
    genotype_digest = genotype_content_digests[[as.character(seed)]],
    selected_pair_keys = genotype_sample$selection$pair_key,
    genotype_selection_summary = genotype_sample$selection_summary,
    truth_maf_balance = truth_maf_balance,
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
    "genotype_digest",
    "selected_pair_keys",
    "genotype_selection_summary",
    "truth_maf_balance",
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
    identical(
      x$genotype_digest,
      genotype_content_digests[[as.character(seed)]]
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
        expected_class_counts[x$truth_maf_balance$effect_class]
    ) &&
    nrow(x$functional_alpha) == expected_rows &&
    nrow(x$functional_alpha_005) == length(methods) * length(targets) &&
    all(methods %in% unique(x$functional_alpha$method)) &&
    all(targets %in% unique(x$functional_alpha$target)) &&
    all(is.finite(x$functional_alpha$power)) &&
    all(is.finite(x$functional_alpha$empirical_fsr)) &&
    nrow(x$truth_group_counts) == length(truth_group_levels) &&
    identical(
      setNames(
        x$truth_group_counts$n_dynamic,
        as.character(x$truth_group_counts$truth_group)
      ),
      expected_truth_group_counts
    )
}

replicates <- list()
replicate_index <- 1L
for (truth_mechanism in run_truth_mechanisms) {
  for (seed in run_seed_list) {
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

if (length(run_truth_mechanisms) > 1L) {
  for (seed in run_seed_list) {
  seed_replicates <- replicates[vapply(
    replicates,
    function(x) identical(x$seed, seed),
    logical(1)
  )]
  if (length(seed_replicates) != length(run_truth_mechanisms) ||
      length(unique(vapply(
        seed_replicates,
        `[[`,
        character(1),
        "genotype_digest"
      ))) != 1L ||
      !all(vapply(
        seed_replicates,
        function(x) identical(
          x$selected_pair_keys,
          genotype_cache$samples[[as.character(seed)]]$selection$pair_key
        ),
        logical(1)
      ))) {
    stop("R3A and R3B do not share the genotype selection for seed ", seed, ".")
  }
  }
}

if (replicate_only) {
  message(
    "Saved requested matched functional-testing replicates to: ",
    replicate_dir
  )
  quit(save = "no", status = 0L)
}

if (!identical(run_seed_list, seed_list) ||
    !identical(run_truth_mechanisms, truth_mechanisms)) {
  stop(
    "Summary generation requires the full configured seed and mechanism lists. ",
    "Use --replicate-only true for subset runs."
  )
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
all_genotype_selection_summary <- do.call(
  rbind,
  lapply(replicates, `[[`, "genotype_selection_summary")
)
all_genotype_selection_summary$truth_mechanism <- rep(
  truth_mechanisms,
  each = length(seed_list)
)
all_truth_maf_balance <- do.call(
  rbind,
  Map(function(replicate, mechanism) {
    out <- replicate$truth_maf_balance
    out$truth_mechanism <- mechanism
    out
  }, replicates, vapply(replicates, `[[`, character(1), "truth_mechanism"))
)
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
  all_genotype_selection_summary,
  file.path(summary_dir, "genotype_selection_summary.csv")
)
write_csv(
  all_truth_maf_balance,
  file.path(summary_dir, "truth_maf_balance.csv")
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
