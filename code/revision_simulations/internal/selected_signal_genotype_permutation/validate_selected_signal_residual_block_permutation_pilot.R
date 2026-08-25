#!/usr/bin/env Rscript

# Independently validate the selected-signal residual-block permutation cache.

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
  "selected_signal_residual_block_permutation_seed20260811"
)
fit_path <- file.path(output_directory, "merged_fash_fit.rds")
effect_path <- file.path(output_directory, "effect_estimates.rds")
unit_lfdr_path <- file.path(output_directory, "unit_lfdr.csv")
diagnostic_path <- file.path(output_directory, "calibration_diagnostics.csv")
donor_map_path <- file.path(output_directory, "donor_permutation.csv")
source_path <- file.path(output_directory, "source_information.csv")
required_paths <- c(
  fit_path,
  effect_path,
  unit_lfdr_path,
  diagnostic_path,
  donor_map_path,
  source_path
)
if (any(!file.exists(required_paths))) {
  stop("The completed residual-block pilot cache is incomplete.")
}

fit_bundle <- readRDS(fit_path)
effect_estimates <- readRDS(effect_path)
unit_lfdr <- utils::read.csv(unit_lfdr_path, stringsAsFactors = FALSE)
saved_diagnostics <- utils::read.csv(
  diagnostic_path,
  stringsAsFactors = FALSE
)
donor_map <- utils::read.csv(donor_map_path, stringsAsFactors = FALSE)
source_information <- utils::read.csv(source_path, stringsAsFactors = FALSE)

configuration <- fit_bundle$configuration
if (!identical(configuration$permutation_method, "residual_block") ||
    configuration$seed != 20260811L || configuration$alpha != 0.05 ||
    configuration$n_target_units != 1177L ||
    configuration$n_permuted_null_units != 1177L ||
    !identical(as.integer(configuration$sample_counts),
               c(19L, 19L, 16L, 19L, 16L, rep(19L, 8L), 18L, 19L, 19L)) ||
    !identical(as.integer(configuration$residual_df),
               as.integer(configuration$sample_counts) - 7L)) {
  stop("The residual-block pilot configuration is invalid.")
}

expected_stages <- c("Raw", "BF-adjusted")
if (!identical(unique(unit_lfdr$fit_stage), expected_stages) ||
    nrow(unit_lfdr) != 4708L ||
    anyDuplicated(unit_lfdr[, c("fit_stage", "unit_key")]) ||
    !setequal(unique(unit_lfdr$group), c("target", "permuted_null")) ||
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
  if (max(abs(saved_lfdr - as.numeric(fit$lfdr))) > 1e-14) {
    stop("Saved and fitted lfdr values disagree for stage: ", stage)
  }
  summarize_matched_null_calibration(
    lfdr = saved_lfdr,
    group = unit_lfdr$group[rows],
    pi0_merged = extract_pi0(fit),
    fit_stage = stage,
    alpha = configuration$alpha
  )
}))
rownames(diagnostics) <- NULL
numeric_diagnostic_columns <- names(diagnostics)[vapply(
  diagnostics,
  is.numeric,
  logical(1)
)]
for (column in numeric_diagnostic_columns) {
  saved_diagnostics[[column]] <- as.numeric(saved_diagnostics[[column]])
}

if (!isTRUE(all.equal(
  diagnostics,
  saved_diagnostics,
  tolerance = 1e-14,
  check.attributes = FALSE
))) {
  stop("The independently recomputed calibration diagnostics disagree.")
}
if (any(abs(
  diagnostics$known_null_discovery_fraction -
    diagnostics$permuted_null_calls / diagnostics$total_calls
) > 1e-15, na.rm = TRUE) ||
    any(abs(
      diagnostics$scaled_fdr_merged_from_design_lower_bound -
        diagnostics$known_null_discovery_fraction
    ) > 1e-15, na.rm = TRUE) ||
    any(
      diagnostics$pi0_target_valid !=
        is.finite(diagnostics$post_selection_fdr_target_from_pi0)
    )) {
  stop("The recomputed calibration identities failed validation.")
}

required_donor_columns <- c(
  "target_donor", "source_donor", "fixed_point", "observation_pattern"
)
if (!all(required_donor_columns %in% names(donor_map)) ||
    nrow(donor_map) != 19L || anyDuplicated(donor_map$target_donor) ||
    anyDuplicated(donor_map$source_donor) ||
    !setequal(donor_map$target_donor, donor_map$source_donor) ||
    all(donor_map$fixed_point)) {
  stop("The residual-block donor map is invalid.")
}
source_pattern <- donor_map$observation_pattern[
  match(donor_map$source_donor, donor_map$target_donor)
]
if (!identical(source_pattern, donor_map$observation_pattern)) {
  stop("The residual-block donor map crossed missingness-pattern strata.")
}

effect_matrices <- effect_estimates[c(
  "target_beta",
  "target_adjusted_se",
  "observed_refit_beta",
  "observed_refit_raw_se",
  "observed_refit_adjusted_se",
  "permuted_beta",
  "permuted_raw_se",
  "permuted_adjusted_se"
)]
if (any(vapply(effect_matrices, function(x) {
  !identical(dim(x), c(1177L, 16L)) || any(!is.finite(x))
}, logical(1))) ||
    any(effect_estimates$target_adjusted_se <= 0) ||
    any(effect_estimates$permuted_raw_se <= 0) ||
    any(effect_estimates$permuted_adjusted_se <= 0)) {
  stop("The saved effect-estimate matrices are invalid.")
}

for (stage in expected_stages) {
  weights <- fits[[stage]]$prior_weights$prior_weight
  if (any(!is.finite(weights)) || any(weights < 0) ||
      abs(sum(weights) - 1) > 1e-8) {
    stop("The fitted mixture weights are invalid for stage: ", stage)
  }
}

if (any(!file.exists(source_information$path)) ||
    !identical(
      unname(tools::md5sum(source_information$path)),
      source_information$md5
    )) {
  stop("At least one immutable source file no longer matches its saved hash.")
}

utils::write.csv(diagnostics, diagnostic_path, row.names = FALSE)
cat("Selected-signal residual-block permutation pilot validation passed.\n")
print(diagnostics, row.names = FALSE)
