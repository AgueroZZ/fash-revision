#!/usr/bin/env Rscript

# Compare the immutable seed-12345 R3B designs with the frozen full-support
# pilot. This diagnostic reports paired retained results and does not select or
# modify any simulation setting.

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

workflowr_root <- find_workflowr_root()
cache_ids <- c(
  closed = paste0(
    "r3_real_genotype_one_per_gene_J6362_matched_functional_",
    "relative_clearance_main_effect_fashr0143_pilot5"
  ),
  open_original = paste0(
    "r3_real_genotype_one_per_gene_J6362_matched_functional_",
    "open_middle_3_12_relative_clearance_main_effect_fashr0143_pilot5"
  ),
  center_aligned = paste0(
    "r3_real_genotype_one_per_gene_J6362_matched_functional_",
    "open_middle_3_12_center_aligned_relative_clearance_",
    "main_effect_fashr0143_pilot5"
  ),
  boundary_shifted = paste0(
    "r3_real_genotype_one_per_gene_J6362_matched_functional_",
    "open_middle_3_12_support_contained_relative_clearance_",
    "main_effect_fashr0143_pilot5"
  ),
  full_support_mixture = paste0(
    "r3b_real_genotype_one_per_gene_J6362_open_middle_3_12_",
    "full_support_iwp1_geometry_mixture_fashr0143_seed12345_pilot"
  )
)

read_seed_rows <- function(design, cache_id) {
  path <- file.path(
    workflowr_root, "output", "revision_simulations", "mc", cache_id,
    "replicates", "raised_cosine_seed_12345.rds"
  )
  if (!file.exists(path)) {
    stop("Missing paired R3B replicate: ", path, call. = FALSE)
  }
  rows <- readRDS(path)$functional_alpha
  rows <- rows[
    rows$method == "FASH-IWP1-BF" & rows$target == "middle",
    , drop = FALSE
  ]
  false_weight <- 0.29 / (1 / 3)
  true_weight <- 0.42 / (1 / 3)
  weighted_false <-
    false_weight * rows$conditional_false_discoveries +
    rows$first_stage_null_calls
  weighted_true <- true_weight * rows$true_positives
  data.frame(
    design = design,
    seed = 12345L,
    alpha = rows$alpha,
    discoveries = rows$n_discoveries,
    false_discoveries = rows$false_discoveries,
    empirical_fsr = rows$empirical_fsr,
    balanced_to_iwp_counterfactual = if (
      identical(design, "full_support_mixture")
    ) {
      NA_real_
    } else {
      weighted_false / (weighted_false + weighted_true)
    },
    stringsAsFactors = FALSE
  )
}

comparison <- do.call(rbind, Map(
  read_seed_rows,
  design = names(cache_ids),
  cache_id = unname(cache_ids)
))
comparison <- comparison[
  order(match(comparison$design, names(cache_ids)), comparison$alpha),
  , drop = FALSE
]

diagnostic_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "diagnostics",
  "r3_middle_calibration"
)
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  comparison,
  file.path(diagnostic_dir, "r3b_seed12345_design_comparison.csv"),
  row.names = FALSE
)

selected_alpha <- c(0.05, 0.10, 0.15, 0.20)
selected <- comparison[vapply(
  comparison$alpha,
  function(alpha) any(abs(alpha - selected_alpha) < 1e-12),
  logical(1)
), , drop = FALSE]
cat("Paired seed-12345 BF Middle comparison:\n")
print(selected, row.names = FALSE)
