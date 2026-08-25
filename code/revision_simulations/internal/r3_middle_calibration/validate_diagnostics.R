#!/usr/bin/env Rscript

# Validate the internal R3 Middle diagnostic artifacts using independent
# probability, count, and cross-representation consistency checks.

options(stringsAsFactors = FALSE)

project_root <- normalizePath(getwd(), mustWork = TRUE)
diagnostic_dir <- file.path(
  project_root,
  "output/revision_simulations/diagnostics/r3_middle_calibration"
)

checks <- list()
check_index <- 0L
add_check <- function(name, passed, detail) {
  check_index <<- check_index + 1L
  checks[[check_index]] <<- data.frame(
    check = name,
    passed = isTRUE(passed),
    detail = as.character(detail),
    stringsAsFactors = FALSE
  )
}

proposal <- read.csv(file.path(
  diagnostic_dir,
  "proposal_category_probabilities.csv"
))
proposal_sums <- aggregate(
  proposal$probability,
  proposal["definition"],
  sum
)
add_check(
  "proposal probabilities sum to one",
  max(abs(proposal_sums$x - 1)) < 1e-12,
  paste(proposal_sums$definition, sprintf("%.12f", proposal_sums$x), collapse = "; ")
)

transition <- read.csv(file.path(
  diagnostic_dir,
  "proposal_transition_matrix.csv"
))
transition_total <- sum(transition$count)
gap_to_middle <- transition$count[
  transition$closed_category == "unassigned" &
    transition$open_category == "middle"
]
closed_gap <- proposal$count[
  proposal$definition == "closed" & proposal$category == "unassigned"
]
add_check(
  "paired proposal transition total",
  transition_total == unique(proposal$proposals),
  paste("transition total =", transition_total)
)
add_check(
  "all closed-gap proposals become open Middle",
  length(gap_to_middle) == 1L && gap_to_middle == closed_gap,
  paste("gap-to-Middle =", gap_to_middle, "; closed gap =", closed_gap)
)

accepted <- readRDS(file.path(
  diagnostic_dir,
  "accepted_truth_unit_metrics.rds"
))
accepted_groups <- aggregate(
  rep(1L, nrow(accepted)),
  accepted[c("definition", "seed")],
  sum
)
add_check(
  "accepted truth count per definition and seed",
  all(accepted_groups$x == 1272L),
  paste("range =", paste(range(accepted_groups$x), collapse = "--"))
)
add_check(
  "accepted truth labels satisfy positive target functionals",
  all(accepted$target_functional >= 0.10 - 1e-10),
  paste("minimum target functional =", min(accepted$target_functional))
)

prior <- read.csv(file.path(
  diagnostic_dir,
  "fitted_dynamic_prior_category_probabilities.csv"
))
prior_key <- interaction(
  prior$fit_set,
  prior$seed,
  prior$method,
  prior$definition,
  prior$component,
  drop = TRUE
)
prior_sums <- vapply(
  split(prior$probability, prior_key),
  sum,
  numeric(1)
)
add_check(
  "fitted dynamic-prior probabilities sum to one",
  max(abs(prior_sums - 1)) < 1e-12,
  paste("maximum absolute deviation =", max(abs(prior_sums - 1)))
)

reference <- read.csv(file.path(
  diagnostic_dir,
  "reference_category_mixture.csv"
))
add_check(
  "reference category mixture is a probability distribution",
  abs(sum(reference$canonical_reference_probability) - 1) < 1e-12 &&
    all(reference$canonical_reference_probability > 0),
  paste(
    reference$category,
    sprintf("%.6f", reference$canonical_reference_probability),
    collapse = "; "
  )
)

by_replicate <- read.csv(file.path(
  diagnostic_dir,
  "middle_mixture_sensitivity_by_replicate.csv"
))
alpha005_calls <- read.csv(file.path(
  diagnostic_dir,
  "alpha005_all_target_mixture_sensitivity.csv"
))
alpha005_aggregate <- by_replicate[
  by_replicate$target == "middle" &
    abs(by_replicate$alpha - 0.05) < 1e-12,
  c(
    "seed", "truth_mechanism", "method",
    "matched_mixture_empirical_fsr"
  ),
  drop = FALSE
]
alpha005_calls <- alpha005_calls[
  alpha005_calls$target == "middle",
  c(
    "seed", "truth_mechanism", "method",
    "matched_mixture_empirical_fsr"
  ),
  drop = FALSE
]
comparison <- merge(
  alpha005_aggregate,
  alpha005_calls,
  by = c("seed", "truth_mechanism", "method"),
  suffixes = c("_aggregate", "_calls"),
  sort = TRUE
)
path_difference <- abs(
  comparison$matched_mixture_empirical_fsr_aggregate -
    comparison$matched_mixture_empirical_fsr_calls
)
add_check(
  "alpha-0.05 mixture sensitivity agrees across aggregate and call paths",
  nrow(comparison) == 20L && max(path_difference) < 1e-12,
  paste("maximum absolute difference =", max(path_difference))
)

gate <- read.csv(file.path(
  diagnostic_dir,
  "candidate_calibration_gate_metrics.csv"
))
add_check(
  "matched-mixture Middle excess satisfies the prespecified gate",
  all(gate$matched_mixture_max_excess <= 0.03),
  paste(
    gate$truth_mechanism,
    sprintf("%.6f", gate$matched_mixture_max_excess),
    collapse = "; "
  )
)

formal_complete <- file.path(
  project_root,
  "output/revision_simulations/mc",
  paste0(
    "r3_real_genotype_one_per_gene_J6362_",
    "matched_functional_open_middle_3_12_center_aligned_",
    "relative_clearance_main_effect_fashr0143_pilot5/complete.flag"
  )
)
add_check(
  "historical center-aligned formal cache remains complete",
  file.exists(formal_complete),
  formal_complete
)

validation <- do.call(rbind, checks)
write.csv(
  validation,
  file.path(diagnostic_dir, "validation_checks.csv"),
  row.names = FALSE
)

if (!all(validation$passed)) {
  print(validation[!validation$passed, , drop = FALSE], row.names = FALSE)
  stop("At least one R3 Middle diagnostic validation failed.")
}

cat("All ", nrow(validation), " diagnostic validation checks passed.\n", sep = "")
