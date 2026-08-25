#!/usr/bin/env Rscript

# Validate the fixed-seed R1 known-truth genotype-permutation pilot.

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
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "r1_known_truth_genotype_permutation",
  "r1_known_truth_genotype_permutation_helpers.R"
))

output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r1_known_truth_genotype_permutation_seed12345_perm20260811_J200"
)
required_files <- c(
  "configuration.rds",
  "donor_permutation.csv",
  "unit_truth.csv",
  "effect_estimates.rds",
  "arm_fits.rds",
  "pi0_summary.csv",
  "alpha_curve.csv",
  "alpha005_summary.csv",
  "null_bf_summary.csv",
  "null_bf_units.csv",
  "genotype_alignment.csv",
  "validation.csv",
  "source_information.csv"
)
required_paths <- file.path(output_directory, required_files)
if (!dir.exists(output_directory) || any(!file.exists(required_paths))) {
  stop("The completed fixed-seed known-truth pilot cache is missing.")
}

configuration <- readRDS(file.path(output_directory, "configuration.rds"))
arm_fits <- readRDS(file.path(output_directory, "arm_fits.rds"))
unit_truth <- utils::read.csv(
  file.path(output_directory, "unit_truth.csv"),
  stringsAsFactors = FALSE
)
pi0_summary <- utils::read.csv(
  file.path(output_directory, "pi0_summary.csv"),
  stringsAsFactors = FALSE
)
alpha005_summary <- utils::read.csv(
  file.path(output_directory, "alpha005_summary.csv"),
  stringsAsFactors = FALSE
)
null_bf_summary <- utils::read.csv(
  file.path(output_directory, "null_bf_summary.csv"),
  stringsAsFactors = FALSE
)
null_bf_units <- utils::read.csv(
  file.path(output_directory, "null_bf_units.csv"),
  stringsAsFactors = FALSE
)
alignment <- utils::read.csv(
  file.path(output_directory, "genotype_alignment.csv"),
  stringsAsFactors = FALSE
)
validation <- utils::read.csv(
  file.path(output_directory, "validation.csv"),
  stringsAsFactors = FALSE
)

expected_arms <- c(
  "genuine_null_baseline",
  "shared_genotype_permutation"
)
if (!identical(sort(names(arm_fits)), sort(expected_arms)) ||
    configuration$n_alternatives != 200L ||
    configuration$n_units_per_arm != 400L ||
    abs(configuration$true_pi0 - 0.5) > 1e-12 ||
    nrow(unit_truth) != 800L ||
    any(table(unit_truth$arm) != 400L) ||
    any(tapply(unit_truth$true_null, unit_truth$arm, sum) != 200L) ||
    nrow(alignment) != 200L || anyDuplicated(alignment$unit_key) ||
    !all(validation$pass)) {
  stop("The known-truth pilot dimensions, truth labels, or validations failed.")
}

numeric_curve_columns <- c(
  "n_units", "n_true_null", "n_true_alternative", "n_discoveries",
  "false_discoveries", "true_positives", "realized_fdp", "power",
  "mean_selected_lfdr", "maximum_selected_lfdr"
)
numeric_bf_columns <- setdiff(names(null_bf_summary), "arm")
recomputed_alpha005 <- list()
recomputed_pi0 <- list()
recomputed_bf_summary <- list()
recomputed_bf_units <- list()
row_index <- 0L
for (arm in expected_arms) {
  bundle <- arm_fits[[arm]]
  if (length(bundle$unit_keys) != 400L ||
      length(bundle$true_null) != 400L ||
      sum(bundle$true_null) != 200L ||
      anyDuplicated(bundle$unit_keys) ||
      nrow(bundle$raw_fit$L_matrix) != 400L ||
      nrow(bundle$bf_fit$L_matrix) != 400L ||
      length(bundle$bf_fit$BF) != 400L ||
      any(!is.finite(bundle$bf_fit$BF)) || any(bundle$bf_fit$BF <= 0)) {
    stop("A saved arm fit is incomplete or misaligned: ", arm)
  }
  for (fit_stage in c("Raw", "BF")) {
    fit <- if (fit_stage == "Raw") bundle$raw_fit else bundle$bf_fit
    row_index <- row_index + 1L
    recomputed_alpha005[[row_index]] <- known_truth_alpha_curve(
      lfdr = fit$lfdr,
      true_null = bundle$true_null,
      alpha_grid = 0.05,
      arm = arm,
      fit_stage = fit_stage
    )
    recomputed_pi0[[row_index]] <- data.frame(
      arm = arm,
      fit_stage = fit_stage,
      estimated_pi0 = extract_pi0(fit),
      true_pi0 = mean(bundle$true_null),
      difference_from_true_pi0 = extract_pi0(fit) - mean(bundle$true_null),
      stringsAsFactors = FALSE
    )
  }
  recomputed_bf_summary[[arm]] <- summarize_null_bf(
    bundle$bf_fit$BF,
    bundle$true_null,
    arm
  )
  recomputed_bf_units[[arm]] <- data.frame(
    arm = arm,
    unit_key = bundle$unit_keys[bundle$true_null],
    source_unit_id = bundle$source_unit_id[bundle$true_null],
    group = bundle$group[bundle$true_null],
    bayes_factor = bundle$bf_fit$BF[bundle$true_null],
    stringsAsFactors = FALSE
  )
}

recomputed_alpha005 <- do.call(rbind, recomputed_alpha005)
recomputed_pi0 <- do.call(rbind, recomputed_pi0)
recomputed_bf_summary <- do.call(rbind, recomputed_bf_summary)
recomputed_bf_units <- do.call(rbind, recomputed_bf_units)
rownames(recomputed_alpha005) <- NULL
rownames(recomputed_pi0) <- NULL
rownames(recomputed_bf_summary) <- NULL
rownames(recomputed_bf_units) <- NULL

sort_rows <- function(x, keys) {
  x[do.call(order, unname(x[keys])), , drop = FALSE]
}
alpha005_summary <- sort_rows(alpha005_summary, c("arm", "fit_stage"))
recomputed_alpha005 <- sort_rows(
  recomputed_alpha005,
  c("arm", "fit_stage")
)
pi0_summary <- sort_rows(pi0_summary, c("arm", "fit_stage"))
recomputed_pi0 <- sort_rows(recomputed_pi0, c("arm", "fit_stage"))
null_bf_summary <- sort_rows(null_bf_summary, "arm")
recomputed_bf_summary <- sort_rows(recomputed_bf_summary, "arm")
null_bf_units <- sort_rows(null_bf_units, c("arm", "unit_key"))
recomputed_bf_units <- sort_rows(
  recomputed_bf_units,
  c("arm", "unit_key")
)
rownames(alpha005_summary) <- NULL
rownames(recomputed_alpha005) <- NULL
rownames(pi0_summary) <- NULL
rownames(recomputed_pi0) <- NULL
rownames(null_bf_summary) <- NULL
rownames(recomputed_bf_summary) <- NULL
rownames(null_bf_units) <- NULL
rownames(recomputed_bf_units) <- NULL

if (!identical(
      alpha005_summary[c("arm", "fit_stage")],
      recomputed_alpha005[c("arm", "fit_stage")]
    ) || max(abs(
      as.matrix(alpha005_summary[numeric_curve_columns]) -
        as.matrix(recomputed_alpha005[numeric_curve_columns])
    ), na.rm = TRUE) > 1e-12 ||
    !identical(
      pi0_summary[c("arm", "fit_stage")],
      recomputed_pi0[c("arm", "fit_stage")]
    ) || max(abs(
      as.matrix(pi0_summary[c(
        "estimated_pi0", "true_pi0", "difference_from_true_pi0"
      )]) -
        as.matrix(recomputed_pi0[c(
          "estimated_pi0", "true_pi0", "difference_from_true_pi0"
        )])
    )) > 1e-12 ||
    !identical(null_bf_summary$arm, recomputed_bf_summary$arm) ||
    max(abs(
      as.matrix(null_bf_summary[numeric_bf_columns]) -
        as.matrix(recomputed_bf_summary[numeric_bf_columns])
    )) > 1e-10 ||
    !identical(
      null_bf_units[c("arm", "unit_key", "source_unit_id", "group")],
      recomputed_bf_units[c("arm", "unit_key", "source_unit_id", "group")]
    ) || max(abs(
      null_bf_units$bayes_factor - recomputed_bf_units$bayes_factor
    )) > 1e-10) {
  stop("Saved known-truth summaries do not reproduce from the saved arm fits.")
}

cat("R1 known-truth genotype-permutation pilot validation passed.\n")
