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
noise_grid <- parse_numeric_grid(get_arg("--noise-grid", "1,0.75,0.5,0.25"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
overwrite <- identical(tolower(get_arg("--overwrite", "false")), "true")
output_dir <- get_arg(
  "--output-dir",
  file.path(
    "output",
    "revision_simulations",
    "diagnostics",
    "switch_noise_ablation"
  )
)

if (!is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
    length(noise_grid) == 0 || any(!is.finite(noise_grid)) ||
    any(noise_grid <= 0) || !is.finite(num_cores) || num_cores < 1) {
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
true_switch <- reference$true_functionals[, "switch"] > 0

evaluate_fit <- function(fit, eqtl_summary, noise_sd, posterior_seed) {
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
  non_switch_dynamic <- which(
    reference$unit_info$effect_class == "dynamic_bspline" & !true_switch
  )
  data.frame(
    expression_noise_sd = noise_sd,
    median_se_all = stats::median(eqtl_summary$se),
    median_se_non_switch_dynamic =
      stats::median(eqtl_summary$se[non_switch_dynamic, ]),
    threshold_over_median_se_non_switch =
      0.25 / stats::median(eqtl_summary$se[non_switch_dynamic, ]),
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

fit_one_noise <- function(noise_sd, noise_index) {
  if (isTRUE(all.equal(noise_sd, reference$settings$expression_noise_sd))) {
    return(list(
      fit = reference$fash_fits$fash_iwp1_bf,
      eqtl_summary = reference$eqtl_summary
    ))
  }
  noise_stem <- gsub("[.]", "p", format(noise_sd, trim = TRUE))
  cache_path <- file.path(
    output_dir,
    paste0("seed", seed, "_noise", noise_stem, "_fash_fit.rds")
  )
  if (file.exists(cache_path) && !overwrite) {
    return(readRDS(cache_path))
  }
  expression_sim <- simulate_eqtl_expression_from_genotypes(
    G = reference$genotype,
    beta_matrix = reference$true_beta,
    time_grid = reference$settings$time_grid,
    covariates = reference$covariates,
    expression_noise_sd = noise_sd,
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
  datasets <- make_fash_datasets_from_eqtl_summary(
    beta_hat = eqtl_summary$beta_hat,
    se = eqtl_summary$se,
    true_beta = reference$true_beta,
    time_grid = reference$settings$time_grid,
    unit_info = reference$unit_info,
    scenario = paste0(reference$settings$scenario, "_noise_", noise_stem)
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
    eqtl_summary = eqtl_summary
  )
  saveRDS(result, cache_path)
  result
}

results <- do.call(rbind, lapply(seq_along(noise_grid), function(i) {
  fitted <- fit_one_noise(noise_grid[i], i)
  evaluate_fit(
    fit = fitted$fit,
    eqtl_summary = fitted$eqtl_summary,
    noise_sd = noise_grid[i],
    posterior_seed = 930000L + seed + i * 1000L
  )
}))
results <- results[order(results$expression_noise_sd, decreasing = TRUE), ]

write.csv(
  results,
  file.path(output_dir, paste0("seed", seed, "_alpha005_noise_comparison.csv")),
  row.names = FALSE
)
print(results)
