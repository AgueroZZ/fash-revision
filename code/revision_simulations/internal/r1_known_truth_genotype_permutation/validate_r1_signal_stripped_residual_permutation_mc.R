#!/usr/bin/env Rscript

# Independently validate the five-seed R1 signal-stripped residual permutation.

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
  "r1_signal_stripped_unadjusted_residual_permutation_mc5_J200_v1"
)
replicate_directory <- file.path(output_directory, "replicates")
summary_directory <- file.path(output_directory, "summary")
source_seeds <- c(12345L, 22345L, 32345L, 42345L, 52345L)
replicate_paths <- file.path(
  replicate_directory, paste0("seed_", source_seeds, ".rds")
)
summary_paths <- file.path(
  summary_directory,
  c(
    "replicate_alpha005.csv", "mc_alpha005_summary.csv",
    "source_reconstruction_validation.csv", "construction_diagnostics.csv"
  )
)
if (any(!file.exists(c(replicate_paths, summary_paths)))) {
  stop("A required residual-permutation artifact is missing.")
}

replicates <- lapply(replicate_paths, readRDS)
saved_alpha005 <- utils::read.csv(
  summary_paths[1L], stringsAsFactors = FALSE
)
saved_mc_alpha005 <- utils::read.csv(
  summary_paths[2L], stringsAsFactors = FALSE
)
source_validation <- utils::read.csv(
  summary_paths[3L], stringsAsFactors = FALSE
)
construction_diagnostics <- utils::read.csv(
  summary_paths[4L], stringsAsFactors = FALSE
)
if (nrow(saved_alpha005) != 10L || nrow(saved_mc_alpha005) != 2L ||
    nrow(source_validation) != 7L || nrow(construction_diagnostics) != 5L ||
    any(source_validation$maximum_absolute_difference >
          source_validation$tolerance) ||
    any(construction_diagnostics$leverage_adjustment != "none") ||
    any(construction_diagnostics$
      maximum_full_residual_design_cross_product >= 1e-10) ||
    any(construction_diagnostics$
      maximum_nuisance_partial_genotype_coefficient >= 1e-10)) {
  stop("A saved residual-permutation summary failed structural validation.")
}

recomputed_rows <- list()
row_index <- 0L
for (seed_index in seq_along(replicates)) {
  replicate <- replicates[[seed_index]]
  if (!identical(replicate$seed, source_seeds[seed_index]) ||
      !all(replicate$validation$pass) || length(replicate$warnings) != 0L) {
    stop("A replicate cache failed provenance or internal validation.")
  }
  unit_groups <- split(
    replicate$unit_results,
    replicate$unit_results$fit_stage,
    drop = TRUE
  )
  if (length(unit_groups) != 2L) {
    stop("A replicate does not contain two fit-stage unit groups.")
  }
  for (units in unit_groups) {
    row_index <- row_index + 1L
    selected_indices <- cumulative_lfdr_calls(units$lfdr, alpha = 0.05)
    selected <- seq_len(nrow(units)) %in% selected_indices
    discoveries <- sum(selected)
    false_discoveries <- sum(selected & units$true_null)
    true_positives <- sum(selected & !units$true_null)
    recomputed_rows[[row_index]] <- data.frame(
      seed = units$seed[1L],
      permutation_seed = units$permutation_seed[1L],
      arm = units$arm[1L],
      fit_stage = units$fit_stage[1L],
      alpha = 0.05,
      n_discoveries = discoveries,
      false_discoveries = false_discoveries,
      true_positives = true_positives,
      realized_fdp = if (discoveries == 0L) 0 else {
        false_discoveries / discoveries
      },
      power = true_positives / sum(!units$true_null),
      stringsAsFactors = FALSE
    )
  }
}
recomputed <- do.call(rbind, recomputed_rows)
key <- c("seed", "permutation_seed", "arm", "fit_stage", "alpha")
metric <- c(
  "n_discoveries", "false_discoveries", "true_positives",
  "realized_fdp", "power"
)
recomputed <- recomputed[do.call(order, recomputed[key]), c(key, metric)]
saved_comparison <- saved_alpha005[
  do.call(order, saved_alpha005[key]), c(key, metric)
]
rownames(recomputed) <- NULL
rownames(saved_comparison) <- NULL
if (!isTRUE(all.equal(recomputed, saved_comparison, tolerance = 1e-12))) {
  stop("Saved alpha-0.05 results do not reproduce from unit lfdr.")
}

means <- aggregate(
  realized_fdp ~ arm + fit_stage,
  data = recomputed,
  FUN = mean
)
saved_means <- saved_mc_alpha005[, c("arm", "fit_stage", "mean_fdr")]
means <- means[order(means$arm, means$fit_stage), ]
saved_means <- saved_means[order(saved_means$arm, saved_means$fit_stage), ]
rownames(means) <- NULL
rownames(saved_means) <- NULL
if (!identical(means[c("arm", "fit_stage")],
               saved_means[c("arm", "fit_stage")]) ||
    max(abs(means$realized_fdp - saved_means$mean_fdr)) > 1e-12) {
  stop("Saved alpha-0.05 Monte Carlo means do not reproduce.")
}

cat("R1 signal-stripped residual-permutation MC validation passed.\n")
