#!/usr/bin/env Rscript

# Evaluate whether the R3 Middle calibration shift is explained by the mismatch
# between the balanced truth-category mixture and the temporal-category mixture
# induced by the fitted dynamic IWP1 prior. This script only reweights retained
# Monte Carlo results; it does not refit FASH or overwrite formal results.

options(stringsAsFactors = FALSE)

project_root <- normalizePath(getwd(), mustWork = TRUE)
diagnostic_dir <- file.path(
  project_root,
  "output/revision_simulations/diagnostics/r3_middle_calibration"
)
cache_dir <- file.path(
  project_root,
  "output/revision_simulations/mc",
  paste0(
    "r3_real_genotype_one_per_gene_J6362_",
    "matched_functional_open_middle_3_12_center_aligned_",
    "relative_clearance_main_effect_fashr0143_pilot5"
  )
)
replicate_dir <- file.path(cache_dir, "replicates")

prior_path <- file.path(
  diagnostic_dir,
  "fitted_dynamic_prior_category_probabilities.csv"
)
canonical_path <- file.path(
  diagnostic_dir,
  "canonical_iwp1_reference_mixture.csv"
)
if (!file.exists(prior_path) || !file.exists(canonical_path)) {
  stop(paste(
    "Run audit_fitted_prior_geometry.R and",
    "derive_canonical_iwp1_reference_mixture.R before this script."
  ))
}

prior_draws <- read.csv(prior_path, check.names = FALSE)
prior_open <- prior_draws[
  prior_draws$component == "all" &
    prior_draws$definition == "open",
  ,
  drop = FALSE
]
fitted_middle_probability <- mean(
  prior_open$probability[prior_open$category == "middle"]
)
canonical_reference <- read.csv(canonical_path, check.names = FALSE)
reference_probability <- setNames(
  canonical_reference$frozen_formal_probability,
  canonical_reference$category
)
if (!identical(names(reference_probability), c("early", "middle", "late")) ||
    abs(sum(reference_probability) - 1) > 1e-12) {
  stop("The frozen canonical reference mixture is invalid.")
}
balanced_probability <- c(early = 1 / 3, middle = 1 / 3, late = 1 / 3)
importance_weight <- reference_probability / balanced_probability
prior_odds_ratio <-
  (balanced_probability[["middle"]] /
     (1 - balanced_probability[["middle"]])) /
  (reference_probability[["middle"]] /
     (1 - reference_probability[["middle"]]))

reference_table <- data.frame(
  category = names(reference_probability),
  balanced_truth_probability = as.numeric(balanced_probability),
  canonical_reference_probability = as.numeric(reference_probability),
  retained_fit_mean_probability = c(
    (1 - fitted_middle_probability) / 2,
    fitted_middle_probability,
    (1 - fitted_middle_probability) / 2
  ),
  importance_weight = as.numeric(importance_weight),
  reference_source = "canonical IWP1 prior, frozen before formal rerun",
  stringsAsFactors = FALSE
)
write.csv(
  reference_table,
  file.path(diagnostic_dir, "reference_category_mixture.csv"),
  row.names = FALSE
)

replicate_paths <- sort(list.files(
  replicate_dir,
  pattern = "^(random_bspline|raised_cosine)_seed_[0-9]+[.]rds$",
  full.names = TRUE
))
if (length(replicate_paths) != 10L) {
  stop("Expected ten retained R3 replicate files.")
}
replicates <- lapply(replicate_paths, readRDS)

middle_rows <- do.call(rbind, lapply(replicates, function(replicate) {
  rows <- replicate$functional_alpha
  rows <- rows[rows$target == "middle", , drop = FALSE]
  false_dynamic_weight <- importance_weight[["early"]]
  true_middle_weight <- importance_weight[["middle"]]
  weighted_false <-
    false_dynamic_weight * rows$conditional_false_discoveries +
    rows$first_stage_null_calls
  weighted_true <- true_middle_weight * rows$true_positives
  rows$matched_mixture_empirical_fsr <-
    weighted_false / (weighted_false + weighted_true)
  rows$prior_odds_mapped_estimated_fsr <-
    rows$estimated_fsr /
    (
      rows$estimated_fsr +
        prior_odds_ratio * (1 - rows$estimated_fsr)
    )
  rows$reference_middle_probability <-
    reference_probability[["middle"]]
  rows$balanced_middle_probability <-
    balanced_probability[["middle"]]
  rows$middle_prior_odds_ratio <- prior_odds_ratio
  rows
}))
write.csv(
  middle_rows,
  file.path(diagnostic_dir, "middle_mixture_sensitivity_by_replicate.csv"),
  row.names = FALSE
)

summarize_mean_ci <- function(x) {
  n <- length(x)
  estimate <- mean(x)
  standard_error <- if (n > 1L) stats::sd(x) / sqrt(n) else NA_real_
  critical_value <- if (n > 1L) stats::qt(0.975, df = n - 1L) else NA_real_
  c(
    mean = estimate,
    mc_se = standard_error,
    ci_lower = max(0, estimate - critical_value * standard_error),
    ci_upper = min(1, estimate + critical_value * standard_error)
  )
}

split_middle <- split(
  middle_rows,
  interaction(
    middle_rows$truth_mechanism,
    middle_rows$method,
    middle_rows$alpha,
    drop = TRUE
  )
)
middle_summary <- do.call(rbind, lapply(split_middle, function(rows) {
  observed <- summarize_mean_ci(rows$empirical_fsr)
  matched <- summarize_mean_ci(rows$matched_mixture_empirical_fsr)
  mapped <- summarize_mean_ci(rows$prior_odds_mapped_estimated_fsr)
  data.frame(
    truth_mechanism = rows$truth_mechanism[[1L]],
    method = rows$method[[1L]],
    alpha = rows$alpha[[1L]],
    replications = nrow(rows),
    observed_mean = observed[["mean"]],
    observed_mc_se = observed[["mc_se"]],
    observed_ci_lower = observed[["ci_lower"]],
    observed_ci_upper = observed[["ci_upper"]],
    matched_mixture_mean = matched[["mean"]],
    matched_mixture_mc_se = matched[["mc_se"]],
    matched_mixture_ci_lower = matched[["ci_lower"]],
    matched_mixture_ci_upper = matched[["ci_upper"]],
    prior_odds_map_mean = mapped[["mean"]],
    prior_odds_map_mc_se = mapped[["mc_se"]],
    stringsAsFactors = FALSE
  )
}))
middle_summary <- middle_summary[
  order(
    middle_summary$truth_mechanism,
    middle_summary$method,
    middle_summary$alpha
  ),
  ,
  drop = FALSE
]
write.csv(
  middle_summary,
  file.path(diagnostic_dir, "middle_mixture_sensitivity_summary.csv"),
  row.names = FALSE
)

call_rows <- do.call(rbind, lapply(replicates, function(replicate) {
  rows <- replicate$call_diagnostics_alpha005
  time_group <- sub(" /.*$", "", rows$truth_group)
  rows$importance_weight <- 1
  dynamic_group <- time_group %in% names(importance_weight)
  rows$importance_weight[dynamic_group] <-
    importance_weight[time_group[dynamic_group]]
  rows
}))
split_calls <- split(
  call_rows,
  interaction(
    call_rows$truth_mechanism,
    call_rows$method,
    call_rows$target,
    call_rows$seed,
    drop = TRUE
  )
)
weighted_calls <- do.call(rbind, lapply(split_calls, function(rows) {
  data.frame(
    seed = rows$seed[[1L]],
    truth_mechanism = rows$truth_mechanism[[1L]],
    method = rows$method[[1L]],
    target = rows$target[[1L]],
    alpha = rows$alpha[[1L]],
    discoveries = nrow(rows),
    observed_empirical_fsr = mean(rows$false_discovery),
    matched_mixture_empirical_fsr =
      sum(rows$importance_weight * rows$false_discovery) /
      sum(rows$importance_weight),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  weighted_calls,
  file.path(diagnostic_dir, "alpha005_all_target_mixture_sensitivity.csv"),
  row.names = FALSE
)

middle_calls <- call_rows[call_rows$target == "middle", , drop = FALSE]
split_middle_calls <- split(
  middle_calls,
  interaction(
    middle_calls$truth_mechanism,
    middle_calls$method,
    middle_calls$seed,
    drop = TRUE
  )
)
odds_corrected_calls <- do.call(rbind, lapply(
  split_middle_calls,
  function(rows) {
    rows$odds_corrected_lfsr <- rows$lfsr /
      (rows$lfsr + prior_odds_ratio * (1 - rows$lfsr))
    rows <- rows[order(rows$odds_corrected_lfsr, rows$unit_index), ]
    rows$odds_corrected_cfsr <-
      cumsum(rows$odds_corrected_lfsr) / seq_len(nrow(rows))
    selected <- rows$odds_corrected_cfsr <= 0.05
    data.frame(
      seed = rows$seed[[1L]],
      truth_mechanism = rows$truth_mechanism[[1L]],
      method = rows$method[[1L]],
      original_discoveries = nrow(rows),
      odds_corrected_discoveries = sum(selected),
      original_empirical_fsr = mean(rows$false_discovery),
      odds_corrected_empirical_fsr = if (any(selected)) {
        mean(rows$false_discovery[selected])
      } else {
        NA_real_
      },
      odds_corrected_estimated_fsr = if (any(selected)) {
        mean(rows$odds_corrected_lfsr[selected])
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }
))
write.csv(
  odds_corrected_calls,
  file.path(
    diagnostic_dir,
    "alpha005_middle_inference_odds_correction.csv"
  ),
  row.names = FALSE
)

cat("Reference open-region category probabilities:\n")
print(reference_table, row.names = FALSE)
cat("\nWrote retained-result mixture sensitivity diagnostics to:\n")
cat(diagnostic_dir, "\n")
