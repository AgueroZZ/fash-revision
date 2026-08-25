#!/usr/bin/env Rscript

# Run compact multi-seed replications for the formal real-genotype random
# B-spline revision simulation. Each seed uses one paper-tested YRI dosage
# variant per gene from the shared R1/R2 genotype cache.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find code/revision_simulations/shared/simulation_functions.R from the current working directory.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0 || hit[1] == length(args)) {
    return(default)
  }
  args[hit[1] + 1]
}

as_flag <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

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
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
))

J <- as.integer(get_arg("--J", "6362"))
n_donors <- as.integer(get_arg("--n-donors", "19"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
expression_noise_sd <- as.numeric(get_arg("--noise-sd", "1"))
dynamic_main_effect_sd <- as.numeric(get_arg("--dynamic-main-effect-sd", "1"))
num_basis <- as.integer(get_arg("--num-basis", "20"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
efdr_permutations <- as.integer(get_arg("--efdr-permutations", "100"))
seed_list <- parse_seed_list(get_arg(
  "--seed-list",
  "12345,22345,32345,42345,52345"
))
output_id <- get_arg(
  "--output-id",
  paste0(
    "r1_real_genotype_one_per_gene_J6362_",
    "random_bspline_main_effect_",
    "linear_mixture_predstep1_penalty10_pilot5"
  )
)
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
    !is.finite(dynamic_main_effect_sd) || dynamic_main_effect_sd <= 0 ||
    num_basis < 2 || num_cores < 1 || efdr_permutations < 1 || !nzchar(output_id)) {
  stop("Invalid simulation arguments.")
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
class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
expected_class_counts <- exact_proportional_counts(J, class_probs)
scenario <- paste0(
  "r1_real_genotype_one_per_gene_",
  "random_bspline_main_effect_dynamic_eqtl"
)
true_pi0_expected <- unname(
  sum(expected_class_counts[c("constant", "zero")]) / J
)
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
invisible(lapply(c(
  output_dir,
  replicate_dir,
  summary_dir,
  figure_dir,
  full_fit_dir
), dir.create,
  recursive = TRUE, showWarnings = FALSE
))

configuration <- list(
  output_id = output_id,
  scenario = scenario,
  J = J,
  n_donors = n_donors,
  time_grid = time_grid,
  n_covariates = n_covariates,
  expression_noise_sd = expression_noise_sd,
  dynamic_main_effect_sd = dynamic_main_effect_sd,
  num_basis = num_basis,
  class_probs = class_probs,
  expected_class_counts = expected_class_counts,
  dynamic_amplitude = 2,
  bspline_df = 6,
  bspline_coefficient_sd = 1,
  efdr_permutations = efdr_permutations,
  true_pi0 = true_pi0_expected,
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
    stop(
      "The existing output id has different settings. Choose a new --output-id or use --overwrite true."
    )
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
  effect_sim <- simulate_variant_effect_curves(
    n_variants = J,
    time_grid = time_grid,
    class_probs = class_probs,
    scenario = scenario,
    dynamic_amplitude = 2,
    bspline_df = 6,
    bspline_coefficient_sd = 1,
    constant_sd = 1,
    dynamic_main_effect_sd = dynamic_main_effect_sd,
    exact_class_counts = TRUE,
    seed = component_seeds[["functional_truth"]]
  )
  effect_sim <- reassign_effect_simulation_by_maf(
    effect_sim = effect_sim,
    maf = genotype_sample$variant_info$observed_maf,
    class_probs = class_probs,
    seed = component_seeds[["classes"]],
    n_strata = 10L
  )
  expression_sim <- simulate_eqtl_expression_from_genotypes(
    G = genotype_sample$G,
    beta_matrix = effect_sim$beta_matrix,
    time_grid = time_grid,
    covariates = covariates,
    expression_noise_sd = expression_noise_sd,
    covariate_effect_sd = 0.5,
    intercept_sd = 0,
    seed = component_seeds[["expression"]]
  )
  out <- run_genotype_level_bspline_eqtl_simulation(
    G = genotype_sample$G,
    time_grid = time_grid,
    covariates = covariates,
    class_probs = class_probs,
    expression_noise_sd = expression_noise_sd,
    dynamic_main_effect_sd = dynamic_main_effect_sd,
    scenario = scenario,
    alpha = 0.05,
    seed = seed,
    grid = common_sd_grid,
    penalty = common_penalty,
    pred_step = common_pred_step,
    linear_prior_mode = linear_prior_mode,
    num_cores = num_cores,
    num_basis = num_basis,
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
  true_pi0 <- dynamic_null_proportion(out$unit_info, target = "dynamic")
  if (!isTRUE(all.equal(true_pi0, true_pi0_expected))) {
    stop("The simulated dynamic-null proportion does not match the configured mixture.")
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
    out,
    n_permutations = efdr_permutations,
    alpha = 0.05,
    seed = seed + 10000L,
    lambda = 0.5,
    pi0_method = "conservative",
    true_pi0 = true_pi0,
    include_true_pi0 = TRUE,
    permute_covariates_with_expression = TRUE,
    num_cores = num_cores,
    overwrite = TRUE,
    verbose = FALSE
  )
  available_methods <- unique(out$result_table$method)
  missing_methods <- setdiff(all_methods, available_methods)
  if (length(missing_methods) > 0) {
    stop("Missing reviewer-facing methods: ", paste(missing_methods, collapse = ", "))
  }
  alpha_curve <- out$alpha_curve[out$alpha_curve$method %in% all_methods, ]
  alpha_curve$seed <- seed
  alpha_005 <- alpha_curve[abs(alpha_curve$alpha - 0.05) < 1e-12, ]
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
  list(
    configuration = configuration,
    seed = seed,
    component_seeds = component_seeds,
    permutation_seed = seed + 10000L,
    true_pi0 = true_pi0,
    genotype_digest = genotype_sample$genotype_digest,
    selected_pair_keys = genotype_sample$selection$pair_key,
    genotype_selection_summary = genotype_sample$selection_summary,
    truth_maf_balance = truth_maf_balance,
    alpha_curve = alpha_curve,
    alpha_005 = alpha_005,
    pi0 = pi0,
    linear_prior_weights = linear_prior_weights,
    linear_prior_summary = linear_prior_summary
  )
}

validate_replicate <- function(replicate, seed) {
  required_fields <- c(
    "configuration", "seed", "component_seeds", "permutation_seed", "true_pi0",
    "genotype_digest", "selected_pair_keys", "genotype_selection_summary",
    "truth_maf_balance",
    "alpha_curve", "alpha_005", "pi0", "linear_prior_weights",
    "linear_prior_summary"
  )
  if (!all(required_fields %in% names(replicate)) ||
      !identical(replicate$seed, seed) ||
      !isTRUE(all.equal(replicate$configuration, configuration)) ||
      !isTRUE(all.equal(replicate$true_pi0, true_pi0_expected)) ||
      !identical(
        replicate$genotype_digest,
        genotype_cache$samples[[as.character(seed)]]$genotype_digest
      ) ||
      !identical(
        replicate$selected_pair_keys,
        genotype_cache$samples[[as.character(seed)]]$selection$pair_key
      )) {
    return(FALSE)
  }
  observed_methods <- unique(replicate$alpha_curve$method)
  prior_weights <- replicate$linear_prior_weights
  prior_summary <- replicate$linear_prior_summary
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
  all(all_methods %in% observed_methods) &&
    is.data.frame(replicate$genotype_selection_summary) &&
    nrow(replicate$genotype_selection_summary) == 1L &&
    replicate$genotype_selection_summary$genes == J &&
    is.data.frame(replicate$truth_maf_balance) &&
    nrow(replicate$truth_maf_balance) == length(class_probs) &&
    all(
      replicate$truth_maf_balance$n ==
        exact_proportional_counts(J, class_probs)[
          replicate$truth_maf_balance$effect_class
        ]
    ) &&
    nrow(replicate$alpha_005) == length(all_methods) &&
    nrow(replicate$pi0) == 4 &&
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
  replicate_path <- file.path(replicate_dir, paste0("seed_", seed, ".rds"))
  if (file.exists(replicate_path) && !overwrite) {
    cached <- readRDS(replicate_path)
    if (validate_replicate(cached, seed)) {
      message("Reusing compact replicate cache: ", replicate_path)
      return(cached)
    }
    stop("Cached replicate does not match the requested settings: ", replicate_path)
  }
  message("Running random B-spline replicate with seed ", seed, ".")
  replicate <- make_replicate(seed)
  saveRDS(replicate, replicate_path)
  replicate
})

all_alpha <- do.call(rbind, lapply(replicates, `[[`, "alpha_curve"))
all_alpha_005 <- do.call(rbind, lapply(replicates, `[[`, "alpha_005"))
all_pi0 <- do.call(rbind, lapply(replicates, `[[`, "pi0"))
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
mc_pi0 <- summarize_mc_pi0(all_pi0)

fash_all_alpha <- all_alpha[all_alpha$method %in% fash_methods, ]
direct_all_alpha <- all_alpha[all_alpha$method %in% direct_methods, ]
fash_mc_alpha <- mc_alpha[mc_alpha$method %in% fash_methods, ]
direct_mc_alpha <- mc_alpha[mc_alpha$method %in% direct_methods, ]
fash_mc_alpha_005 <- mc_alpha_005[mc_alpha_005$method %in% fash_methods, ]
direct_mc_alpha_005 <- mc_alpha_005[mc_alpha_005$method %in% direct_methods, ]

write_csv(all_alpha, file.path(summary_dir, "all_replicate_alpha_curves.csv"))
write_csv(all_alpha_005, file.path(summary_dir, "all_replicate_alpha005.csv"))
write_csv(all_pi0, file.path(summary_dir, "all_replicate_pi0.csv"))
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
write_csv(mc_pi0, file.path(summary_dir, "mc_pi0_summary.csv"))
write_csv(fash_all_alpha, file.path(summary_dir, "iwp_vs_linear_fash_replicate_alpha_curves.csv"))
write_csv(fash_mc_alpha, file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha_curve.csv"))
write_csv(fash_mc_alpha_005, file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha005_summary.csv"))
write_csv(direct_all_alpha, file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_replicate_alpha_curves.csv"))
write_csv(direct_mc_alpha, file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_mc_alpha_curve.csv"))
write_csv(direct_mc_alpha_005, file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_mc_alpha005_summary.csv"))

plot_subtitle <- paste0(
  length(seed_list), " independent seeds; N = ", n_donors,
  ", J = ", J, ", five covariates"
)
plot_mc_alpha_curves(
  fash_mc_alpha,
  metric = "power",
  file = file.path(figure_dir, "iwp_vs_linear_fash_mc_power.png"),
  title = "IWP versus linear FASH: Monte Carlo power",
  subtitle = plot_subtitle,
  style_profile = "combined"
)
plot_mc_alpha_curves(
  fash_mc_alpha,
  metric = "fdr",
  file = file.path(figure_dir, "iwp_vs_linear_fash_mc_fdr.png"),
  title = "IWP versus linear FASH: Monte Carlo FDR estimate",
  subtitle = plot_subtitle,
  legend_position = "topleft",
  style_profile = "combined"
)
plot_mc_alpha_curves(
  direct_mc_alpha,
  metric = "power",
  file = file.path(figure_dir, "iwp_fash_vs_direct_true_pi0_mc_power.png"),
  title = "IWP FASH versus direct interaction: Monte Carlo power",
  subtitle = paste0(plot_subtitle, "; direct eFDR uses true pi0 = ", true_pi0_expected),
  style_profile = "combined"
)
plot_mc_alpha_curves(
  direct_mc_alpha,
  metric = "fdr",
  file = file.path(figure_dir, "iwp_fash_vs_direct_true_pi0_mc_fdr.png"),
  title = "IWP FASH versus direct interaction: Monte Carlo FDR estimate",
  subtitle = paste0(plot_subtitle, "; direct eFDR uses true pi0 = ", true_pi0_expected),
  legend_position = "bottomright",
  style_profile = "combined"
)

print(fash_mc_alpha_005)
print(direct_mc_alpha_005)
print(mc_pi0)
print(all_linear_prior_summary)
