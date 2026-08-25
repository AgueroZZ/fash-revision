#!/usr/bin/env Rscript

# Independently validate and refresh summaries from a completed pilot cache.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation",
  "selected_signal_genotype_permutation_helpers.R"
))

output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "selected_signal_genotype_permutation_seed20260811"
)
fit_path <- file.path(output_directory, "merged_fash_fit.rds")
unit_lfdr_path <- file.path(output_directory, "unit_lfdr.csv")
diagnostic_path <- file.path(output_directory, "calibration_diagnostics.csv")
if (!file.exists(fit_path) || !file.exists(unit_lfdr_path)) {
  stop("The completed fixed-seed pilot cache is missing.")
}

fit_bundle <- readRDS(fit_path)
unit_lfdr <- utils::read.csv(unit_lfdr_path, stringsAsFactors = FALSE)
expected_stages <- c("Raw", "BF-adjusted")
if (!identical(unique(unit_lfdr$fit_stage), expected_stages) ||
    nrow(unit_lfdr) != 4708L ||
    anyDuplicated(unit_lfdr[, c("fit_stage", "unit_key")]) ||
    any(!is.finite(unit_lfdr$lfdr)) ||
    any(unit_lfdr$lfdr < 0 | unit_lfdr$lfdr > 1)) {
  stop("The saved unit-level lfdr table is invalid.")
}

fits <- list(
  Raw = fit_bundle$raw_fit,
  `BF-adjusted` = fit_bundle$bf_adjusted_fit
)
diagnostics <- do.call(rbind, lapply(expected_stages, function(stage) {
  rows <- unit_lfdr$fit_stage == stage
  fit <- fits[[stage]]
  saved_lfdr <- unit_lfdr$lfdr[rows]
  maximum_lfdr_roundtrip_difference <- max(abs(
    saved_lfdr - as.numeric(fit$lfdr)
  ))
  if (maximum_lfdr_roundtrip_difference > 1e-14) {
    stop("Saved and fitted lfdr values disagree for stage: ", stage)
  }
  summarize_matched_null_calibration(
    lfdr = saved_lfdr,
    group = unit_lfdr$group[rows],
    pi0_merged = extract_pi0(fit),
    fit_stage = stage,
    alpha = fit_bundle$configuration$alpha
  )
}))
rownames(diagnostics) <- NULL

if (any(abs(
  diagnostics$known_null_discovery_fraction -
    diagnostics$permuted_null_calls / diagnostics$total_calls
) > 1e-15) ||
    any(abs(
      diagnostics$scaled_fdr_merged_from_design_lower_bound -
        diagnostics$known_null_discovery_fraction
    ) > 1e-15) ||
    any(!diagnostics$pi0_merged_below_design_lower_bound) ||
    any(diagnostics$pi0_target_valid) ||
    any(!is.na(diagnostics$post_selection_fdr_target_from_pi0))) {
  stop("The independently recomputed diagnostics failed validation.")
}

utils::write.csv(diagnostics, diagnostic_path, row.names = FALSE)
cat("Selected-signal genotype-permutation pilot validation passed.\n")
print(diagnostics, row.names = FALSE)
