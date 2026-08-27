#!/usr/bin/env Rscript

# Focused tests for the fash-knockoff reporting helpers.

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
source(file.path(workflowr_root, "code", "revision_simulations", "internal",
                 "fash_knockoff", "fash_knockoff_helpers.R"))
source(file.path(workflowr_root, "code", "revision_simulations", "internal",
                 "fash_knockoff", "reporting.R"))

failures <- 0L
check <- function(label, condition) {
  if (isTRUE(condition)) {
    message("PASS: ", label)
  } else {
    failures <<- failures + 1L
    message("FAIL: ", label)
  }
}

statistics <- load_fash_knockoff_statistics(workflowr_root)

check(
  "the pair set matches the completed run",
  nrow(statistics) == 1009136L && sum(statistics$discovered) == 9214L &&
    length(unique(statistics$gene_id[statistics$discovered])) == 1176L
)
check(
  "the competition statistic reproduces from the stored Bayes factors",
  isTRUE(all.equal(
    statistics$W_updated_weights,
    log(statistics$target_bf_updated_weights) -
      log(statistics$decoy_bf_updated_weights),
    tolerance = 1e-10
  ))
)

group_summary <- fash_knockoff_group_summary(statistics)
check(
  "the non-discovered group is centred on zero",
  abs(group_summary$median[group_summary$group == "not discovered"]) < 1e-3 &&
    abs(group_summary$frac_positive[
      group_summary$group == "not discovered"] - 0.5) < 0.005
)

evidence <- fash_knockoff_sign_evidence(statistics)
naive <- evidence$value[evidence$quantity ==
                          "naive sign-test z (independent pairs)"]
robust <- evidence$value[evidence$quantity ==
                           "cluster-robust z (genes as clusters)"]
check(
  "the naive sign test overstates the evidence relative to cluster-robust",
  naive > 4 * robust && robust > 1.5 && robust < 3
)
check(
  "the implied within-gene sign correlation is positive and well below one",
  { rho <- evidence$value[evidence$quantity ==
                            "implied within-gene sign correlation"]
    rho > 0.05 && rho < 0.5 }
)

bands <- fash_knockoff_lfdr_bands(statistics)
check(
  "the lfdr band table totals the pairs below 0.05 on both arms",
  { row <- bands[bands$lfdr_band == "all < 0.05", ]
    row$target == sum(statistics$target_lfdr < 0.05) &&
      row$decoy == sum(statistics$decoy_lfdr < 0.05) }
)

tails <- fash_knockoff_tail_concentration(statistics, c(10, 25))
check(
  "collapsing the tail to genes changes the ratio materially",
  tails$gene_ratio[tails$threshold == 10] >
    tails$pair_ratio[tails$threshold == 10]
)

selected <- fash_knockoff_selected_pairs(statistics, 0.15)
check(
  "the q = 0.15 selection is a small set inside few genes",
  nrow(selected) == 17L && length(unique(selected$gene_id)) == 6L &&
    all(selected$log10_target_BF > 9) && all(selected$log10_decoy_BF < 1)
)
check(
  "an unreachable q returns an empty selection rather than an error",
  nrow(fash_knockoff_selected_pairs(statistics, 0.01)) == 0L
)

comparison <- fash_knockoff_gene_statistic_comparison(statistics, c(0.05, 0.20))
check(
  "every gene-level statistic agrees that nothing survives q = 0.05",
  nrow(comparison) == 5L &&
    all(comparison$min_estimated_fdr > 0.05)
)

if (failures > 0L) stop(failures, " reporting test(s) failed.")
message("All fash-knockoff reporting tests passed.")
