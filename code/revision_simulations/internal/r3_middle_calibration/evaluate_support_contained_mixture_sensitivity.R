#!/usr/bin/env Rscript

# Reweight the immutable balanced support-contained R3 replicates to the frozen
# canonical IWP1 temporal-category mixture. This is a counterfactual diagnostic,
# not a substitute for the prespecified formal rerun.

options(stringsAsFactors = FALSE)

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }
  if (file.exists(file.path(
    "coderepo-local", "code", "revision_simulations", "shared",
    "simulation_functions.R"
  ))) {
    return(normalizePath(
      "coderepo-local", winslash = "/", mustWork = TRUE
    ))
  }
  stop("Could not find the workflowr repository root.", call. = FALSE)
}

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

workflowr_root <- find_workflowr_root()
diagnostic_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "diagnostics",
  "r3_middle_calibration"
)
canonical <- utils::read.csv(file.path(
  diagnostic_dir, "canonical_iwp1_reference_mixture.csv"
))
reference_probability <- stats::setNames(
  canonical$frozen_formal_probability, canonical$category
)
if (!identical(names(reference_probability), c("early", "middle", "late")) ||
    abs(sum(reference_probability) - 1) > 1e-12) {
  stop("The frozen canonical reference mixture is invalid.", call. = FALSE)
}
importance_weight <- reference_probability / (1 / 3)

cache_id <- paste0(
  "r3_real_genotype_one_per_gene_J6362_matched_functional_",
  "open_middle_3_12_support_contained_relative_clearance_",
  "main_effect_fashr0143_pilot5"
)
replicate_paths <- sort(list.files(
  file.path(
    workflowr_root, "output", "revision_simulations", "mc", cache_id,
    "replicates"
  ),
  pattern = "^(random_bspline|raised_cosine)_seed_[0-9]+[.]rds$",
  full.names = TRUE
))
if (length(replicate_paths) != 10L) {
  stop("Expected ten retained formal replicate files.", call. = FALSE)
}

rows <- do.call(rbind, lapply(replicate_paths, function(path) {
  result <- readRDS(path)$functional_alpha
  result <- result[result$target == "middle", , drop = FALSE]
  result$matched_mixture_false_discoveries <-
    importance_weight[["early"]] * result$conditional_false_discoveries +
    result$first_stage_null_calls
  result$matched_mixture_true_positives <-
    importance_weight[["middle"]] * result$true_positives
  result$matched_mixture_empirical_fsr <- with(
    result,
    matched_mixture_false_discoveries /
      (matched_mixture_false_discoveries + matched_mixture_true_positives)
  )
  result
}))
utils::write.csv(
  rows,
  file.path(
    diagnostic_dir,
    "support_contained_mixture_sensitivity_by_replicate.csv"
  ),
  row.names = FALSE
)

split_rows <- split(
  rows,
  interaction(rows$truth_mechanism, rows$method, rows$alpha, drop = TRUE)
)
summary <- do.call(rbind, lapply(split_rows, function(group) {
  observed <- summarize_mean_ci(group$empirical_fsr)
  matched <- summarize_mean_ci(group$matched_mixture_empirical_fsr)
  data.frame(
    truth_mechanism = group$truth_mechanism[[1L]],
    method = group$method[[1L]],
    alpha = group$alpha[[1L]],
    replications = nrow(group),
    mean_discoveries = mean(group$n_discoveries),
    observed_mean = observed[["mean"]],
    observed_mc_se = observed[["mc_se"]],
    observed_ci_lower = observed[["ci_lower"]],
    observed_ci_upper = observed[["ci_upper"]],
    matched_mixture_mean = matched[["mean"]],
    matched_mixture_mc_se = matched[["mc_se"]],
    matched_mixture_ci_lower = matched[["ci_lower"]],
    matched_mixture_ci_upper = matched[["ci_upper"]],
    stringsAsFactors = FALSE
  )
}))
summary <- summary[
  order(summary$truth_mechanism, summary$method, summary$alpha),
  , drop = FALSE
]
utils::write.csv(
  summary,
  file.path(
    diagnostic_dir, "support_contained_mixture_sensitivity_summary.csv"
  ),
  row.names = FALSE
)

bf <- summary[summary$method == "FASH-IWP1-BF", , drop = FALSE]
gate <- do.call(rbind, lapply(split(bf, bf$truth_mechanism), function(group) {
  excess <- group$matched_mixture_mean - group$alpha
  index <- which.max(excess)
  data.frame(
    truth_mechanism = group$truth_mechanism[[1L]],
    maximum_excess = excess[[index]],
    alpha_at_maximum = group$alpha[[index]],
    matched_mixture_fsr_at_maximum =
      group$matched_mixture_mean[[index]],
    prespecified_tolerance = 0.03,
    predicted_gate_pass = excess[[index]] <= 0.03,
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  gate,
  file.path(
    diagnostic_dir, "support_contained_mixture_predicted_gate.csv"
  ),
  row.names = FALSE
)

selected <- bf[vapply(
  bf$alpha,
  function(alpha) any(abs(alpha - c(0.05, 0.10, 0.15, 0.20)) < 1e-12),
  logical(1)
), c(
  "truth_mechanism", "alpha", "mean_discoveries", "observed_mean",
  "matched_mixture_mean", "matched_mixture_ci_lower",
  "matched_mixture_ci_upper"
), drop = FALSE]
cat("Frozen mixture weights:\n")
print(importance_weight)
cat("\nSelected BF Middle sensitivity rows:\n")
print(selected, row.names = FALSE)
cat("\nPredicted formal gate:\n")
print(gate, row.names = FALSE)
