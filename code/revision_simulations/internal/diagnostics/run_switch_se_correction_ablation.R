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

workflowr_root <- normalizePath(
  get_arg("--workflowr-root", "."),
  mustWork = TRUE
)
seed <- as.integer(get_arg("--seed", "12345"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
overwrite <- identical(tolower(get_arg("--overwrite", "false")), "true")
output_dir <- get_arg(
  "--output-dir",
  file.path(
    "output",
    "revision_simulations",
    "diagnostics",
    "switch_se_correction"
  )
)

if (!is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
    !is.finite(num_cores) || num_cores < 1) {
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
out <- readRDS(raw_candidates)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

uncorrected_fit_path <- file.path(
  output_dir,
  paste0("seed", seed, "_uncorrected_se_fash_fit.rds")
)
if (file.exists(uncorrected_fit_path) && !overwrite) {
  uncorrected_fits <- readRDS(uncorrected_fit_path)
} else {
  uncorrected_datasets <- make_fash_datasets_from_eqtl_summary(
    beta_hat = out$eqtl_summary$beta_hat,
    se = out$eqtl_summary$se_uncorrected,
    true_beta = out$true_beta,
    time_grid = out$settings$time_grid,
    unit_info = out$unit_info,
    scenario = paste0(out$settings$scenario, "_uncorrected_se")
  )
  uncorrected_fits <- fit_fash_for_revision(
    datasets = uncorrected_datasets,
    orders = 1,
    grid = default_revision_grid(),
    num_basis = out$settings$num_basis,
    penalty = out$settings$penalty,
    num_cores = num_cores,
    apply_bf = TRUE,
    verbose = FALSE
  )
  saveRDS(uncorrected_fits, uncorrected_fit_path)
}

switch_functional <- make_temporal_functionals(
  smooth_var = out$evaluation_grid,
  switch_threshold = 0.25
)[["switch"]]
true_switch <- out$true_functionals[, "switch"] > 0

evaluate_fit <- function(fit, se_setting, posterior_seed) {
  fdr_table <- get_fash_fdr_table(fit)
  dynamic_indices <- sort(unique(
    fdr_table$index[fdr_table$FDR <= alpha]
  ))
  set.seed(posterior_seed)
  lfsr_map <- compute_functional_lfsr(
    fit = fit,
    functionals = list(switch = switch_functional),
    indices = dynamic_indices,
    smooth_var = out$evaluation_grid,
    num_cores = num_cores,
    seed = posterior_seed
  )$switch
  lfsr <- unname(lfsr_map[as.character(dynamic_indices)])
  cfsr_table <- functional_cfsr_table(dynamic_indices, lfsr)
  calls <- cfsr_table$index[cfsr_table$cfsr <= alpha]
  data.frame(
    se_setting = se_setting,
    method = "FASH-IWP1-BF",
    alpha = alpha,
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

results <- rbind(
  evaluate_fit(
    out$fash_fits$fash_iwp1_bf,
    "t-corrected",
    910000L + seed
  ),
  evaluate_fit(
    uncorrected_fits$fash_iwp1_bf,
    "uncorrected",
    920000L + seed
  )
)

write.csv(
  results,
  file.path(output_dir, paste0("seed", seed, "_alpha005_comparison.csv")),
  row.names = FALSE
)

se_ratio <- as.numeric(out$eqtl_summary$se / out$eqtl_summary$se_uncorrected)
se_summary <- data.frame(
  minimum = min(se_ratio),
  q025 = unname(stats::quantile(se_ratio, 0.025)),
  median = stats::median(se_ratio),
  q975 = unname(stats::quantile(se_ratio, 0.975)),
  maximum = max(se_ratio)
)
write.csv(
  se_summary,
  file.path(output_dir, paste0("seed", seed, "_se_ratio_summary.csv")),
  row.names = FALSE
)

print(results)
print(se_summary)
