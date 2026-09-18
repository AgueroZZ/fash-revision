#!/usr/bin/env Rscript

# Hand-calculated inference fixtures and checks against the frozen R3 baseline.
experiment <- "code/revision_simulations/internal/r3_stage1_method_comparison"
source(file.path(experiment, "comparison_helpers.R"))
source(file.path(experiment, "reporting.R"))
must_fail <- function(code) stopifnot(inherits(try(force(code), silent = TRUE), "try-error"))

# At 0.05 the two selected units are one true dynamic and one null; at 0.01 none are selected.
scores <- list(fixture = data.frame(unit_index = 1:4, score = c(0.02, 0.06, 0.10, 0.40),
                                   adjusted_score = c(0.02, 0.04, 0.06, 0.145)))
truth <- c(TRUE, FALSE, TRUE, FALSE)
x <- r3c_curves(scores, truth, 1L, "fixture", c(0.01, 0.05, 0.10))
stopifnot(identical(x$calls, c(0L, 2L, 3L)), identical(x$power, c(0, 0.5, 1)),
          isTRUE(all.equal(x$empirical_fdr, c(0, 0.5, 1 / 3))))
bad <- scores
bad[[1]] <- bad[[1]][c(2, 1, 3, 4), ]
must_fail(r3c_curves(bad, truth, 1L, "fixture", 0.05))
bad[[1]] <- scores[[1]]
bad[[1]]$adjusted_score[1] <- NA_real_
must_fail(r3c_curves(bad, truth, 1L, "fixture", 0.05))

# Rate averages are unweighted means across replicate-specific rates.
second <- x
second$seed <- 2L
second$power <- c(0, 1, 1)
second$empirical_fdr <- c(0, 0, 0.5)
fixture_contract <- list(seeds = 1:2, mechanisms = "fixture", methods = "fixture", alpha_grid = c(0.01, 0.05, 0.10))
fixture <- rbind(x, second)
summary <- r3c_summarize(fixture, fixture_contract)
stopifnot(summary$mean_power[2] == 0.75, summary$mean_empirical_fdr[2] == 0.25,
          summary$power_min[2] == 0.5, summary$power_max[2] == 1)
must_fail(r3c_summarize(fixture[-1, ], fixture_contract))
must_fail(r3c_summarize(rbind(fixture, fixture[1, ]), fixture_contract))

base <- r3c_base(".")
context <- r3c_context(".", base)
# With true pi0 = 1/2 and this null CDF, monotone eFDR q-values are 0, 1/6, 1/6, 1/2.
efdr <- context$cmp$empirical_fdr_qvalues(c(0.01, 0.02, 0.5, 0.9),
  c(0.02, 0.8, 0.85, 0.9), pi0_method = "fixed", fixed_pi0 = 0.5)
stopifnot(isTRUE(all.equal(efdr$qvalue, c(0, 1 / 6, 1 / 6, 0.5))))
fixture_scores <- list(direct = data.frame(unit_index = 1:4, score = c(0.01, 0.02, 0.5, 0.9),
                                          adjusted_score = efdr$qvalue))
diagnostics <- list(direct = list(threshold_table = efdr$threshold_table))
r3c_validate_direct(fixture_scores, diagnostics, 4L, 0.5)
diagnostics$direct$threshold_table$null_count[1] <- 1L
must_fail(r3c_validate_direct(fixture_scores, diagnostics, 4L, 0.5))

# Reproduce every existing IWP1 BF decision at every alpha from the retained all-unit table.
iwp_curves <- do.call(rbind, lapply(base$records, function(record) {
  table <- record$inference$fdr_table
  n <- length(record$inference$true_dynamic)
  aligned <- r3c_scores_from_fdr(table[nrow(table):1, ], n)
  r3c_curves(list("FASH-IWP1-BF" = aligned), record$inference$true_dynamic,
              record$seed, record$truth_mechanism)
}))
iwp_curves <- iwp_curves[, names(base$tables$stage1_alpha_by_seed)]
rownames(iwp_curves) <- NULL
stopifnot(isTRUE(all.equal(iwp_curves, base$tables$stage1_alpha_by_seed)))
reference <- base$records[[1]]$inference$fdr_table
reference$index[1] <- reference$index[2]
must_fail(r3c_scores_from_fdr(reference, nrow(reference)))
cat("Selection, per-seed aggregation, eFDR reconstruction, rejection checks, and all 800 baseline IWP1 points passed.\n")
