#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit)) {
    return(default)
  }
  if (hit == length(args)) {
    stop("Missing value after ", flag, ".")
  }
  args[[hit + 1L]]
}

parse_numeric_grid <- function(x) {
  as.numeric(strsplit(x, ",", fixed = TRUE)[[1]])
}

workflowr_root <- normalizePath(
  get_arg("--workflowr-root", "."),
  mustWork = TRUE
)
seed <- as.integer(get_arg("--seed", "12345"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
baseline_grid <- parse_numeric_grid(
  get_arg("--baseline-grid", "0,0.25,0.5,0.75")
)
num_cores <- as.integer(get_arg("--num-cores", "4"))
overwrite <- identical(tolower(get_arg("--overwrite", "false")), "true")
se_mode <- get_arg("--se-mode", "pipeline")
output_dir <- get_arg(
  "--output-dir",
  file.path(
    "output",
    "revision_simulations",
    "diagnostics",
    "switch_non_switch_baseline_ablation"
  )
)

if (!se_mode %in% c("pipeline", "fixed-reference") ||
    !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
    length(baseline_grid) == 0 || any(!is.finite(baseline_grid)) ||
    any(baseline_grid < 0) || !is.finite(num_cores) || num_cores < 1) {
  stop("Invalid ablation arguments.")
}

setwd(workflowr_root)
source("code/revision_simulations/shared/simulation_functions.R")

raw_candidates <- list.files(
  file.path("output", "revision_simulations", "raw"),
  pattern = paste0(
    "^genotype_sparse_timed_cosine_one_two_peak.*seed",
    seed,
    "[.]rds$"
  ),
  full.names = TRUE
)
if (length(raw_candidates) != 1L) {
  stop("Expected exactly one fitted sparse timed cosine object for seed ", seed, ".")
}
reference <- readRDS(raw_candidates)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

switch_functional <- make_temporal_functionals(
  smooth_var = reference$evaluation_grid,
  switch_threshold = 0.25
)[["switch"]]

make_baseline_truth <- function(fraction) {
  beta_observed <- reference$true_beta
  beta_evaluation <- reference$true_beta_evaluation
  non_switch <- which(
    reference$unit_info$effect_class == "dynamic_bspline" &
      reference$unit_info$switch_status == "non-switch"
  )
  for (index in non_switch) {
    peak_sign <- sign(reference$unit_info$peak_signs[[index]][1])
    peak_size <- max(abs(beta_evaluation[index, ]))
    baseline <- peak_sign * fraction * peak_size
    beta_observed[index, ] <- beta_observed[index, ] + baseline
    beta_evaluation[index, ] <- beta_evaluation[index, ] + baseline
  }
  true_functionals <- evaluate_temporal_functionals(
    curves = beta_evaluation,
    smooth_var = reference$evaluation_grid,
    switch_threshold = 0.25
  )
  if (any(true_functionals[non_switch, "switch"] > 0)) {
    stop("A non-switch curve became a switch after adding its same-sign baseline.")
  }
  list(
    beta_observed = beta_observed,
    beta_evaluation = beta_evaluation,
    true_functionals = true_functionals
  )
}

evaluate_fit <- function(fit, truth, baseline_fraction, posterior_seed) {
  true_switch <- truth$true_functionals[, "switch"] > 0
  fdr_table <- get_fash_fdr_table(fit)
  dynamic_indices <- sort(unique(
    fdr_table$index[fdr_table$FDR <= alpha]
  ))
  lfsr_map <- compute_functional_lfsr(
    fit = fit,
    functionals = list(switch = switch_functional),
    indices = dynamic_indices,
    smooth_var = reference$evaluation_grid,
    num_cores = num_cores,
    seed = posterior_seed
  )$switch
  lfsr <- unname(lfsr_map[as.character(dynamic_indices)])
  cfsr_table <- functional_cfsr_table(dynamic_indices, lfsr)
  calls <- cfsr_table$index[cfsr_table$cfsr <= alpha]
  data.frame(
    non_switch_baseline_fraction = baseline_fraction,
    se_mode = se_mode,
    estimated_pi0 = constant_component_prior_weight(fit),
    dynamic_discoveries = length(dynamic_indices),
    switch_calls = length(calls),
    false_switch_calls = sum(!true_switch[calls]),
    true_switch_calls = sum(true_switch[calls]),
    power = sum(true_switch[calls]) / sum(true_switch),
    empirical_fsr = if (length(calls) == 0) 0 else {
      mean(!true_switch[calls])
    },
    estimated_fsr = if (length(calls) == 0) 0 else {
      mean(lfsr[dynamic_indices %in% calls])
    },
    stringsAsFactors = FALSE
  )
}

fit_one_baseline <- function(fraction) {
  truth <- make_baseline_truth(fraction)
  if (fraction == 0) {
    return(list(
      fit = reference$fash_fits$fash_iwp1_bf,
      truth = truth
    ))
  }
  fraction_stem <- gsub("[.]", "p", format(fraction, trim = TRUE))
  cache_path <- file.path(
    output_dir,
    paste0(
      "seed",
      seed,
      "_baseline",
      fraction_stem,
      "_",
      gsub("-", "_", se_mode),
      "_fash_fit.rds"
    )
  )
  if (file.exists(cache_path) && !overwrite) {
    cached <- readRDS(cache_path)
    cached$truth <- truth
    return(cached)
  }
  if (se_mode == "pipeline") {
    expression_sim <- simulate_eqtl_expression_from_genotypes(
      G = reference$genotype,
      beta_matrix = truth$beta_observed,
      time_grid = reference$settings$time_grid,
      covariates = reference$covariates,
      expression_noise_sd = reference$settings$expression_noise_sd,
      covariate_effect_sd = reference$settings$covariate_effect_sd,
      intercept_sd = reference$settings$intercept_sd,
      seed = reference$component_seeds[["expression"]]
    )
    eqtl_summary <- estimate_eqtl_summaries_from_genotypes(
      G = reference$genotype,
      expression = expression_sim$expression,
      covariates = reference$covariates,
      apply_t_se_correction = TRUE
    )
  } else {
    baseline_shift <- truth$beta_observed - reference$true_beta
    eqtl_summary <- reference$eqtl_summary
    eqtl_summary$beta_hat <-
      reference$eqtl_summary$beta_hat + baseline_shift
  }
  datasets <- make_fash_datasets_from_eqtl_summary(
    beta_hat = eqtl_summary$beta_hat,
    se = eqtl_summary$se,
    true_beta = truth$beta_observed,
    time_grid = reference$settings$time_grid,
    unit_info = reference$unit_info,
    scenario = paste0(
      reference$settings$scenario,
      "_non_switch_baseline_",
      fraction_stem
    )
  )
  fits <- fit_fash_for_revision(
    datasets = datasets,
    orders = 1,
    grid = default_revision_grid(),
    num_basis = reference$settings$num_basis,
    penalty = reference$settings$penalty,
    num_cores = num_cores,
    apply_bf = TRUE,
    verbose = FALSE
  )
  result <- list(
    fit = fits$fash_iwp1_bf,
    eqtl_summary = eqtl_summary,
    truth = truth
  )
  saveRDS(result, cache_path)
  result
}

results <- do.call(rbind, lapply(seq_along(baseline_grid), function(i) {
  fitted <- fit_one_baseline(baseline_grid[i])
  evaluate_fit(
    fit = fitted$fit,
    truth = fitted$truth,
    baseline_fraction = baseline_grid[i],
    posterior_seed = 940000L + seed + i * 1000L
  )
}))
results <- results[
  order(results$non_switch_baseline_fraction),
  ,
  drop = FALSE
]

write.csv(
  results,
  file.path(
    output_dir,
    paste0(
      "seed",
      seed,
      "_alpha005_",
      gsub("-", "_", se_mode),
      "_baseline_comparison.csv"
    )
  ),
  row.names = FALSE
)
print(results)
