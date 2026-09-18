#!/usr/bin/env Rscript

# Hand-calculated tests of selection, denominators, and Monte Carlo aggregation.
experiment <- "code/revision_simulations/internal/r3_two_stage_functional_testing"
source("code/revision_simulations/r3_r4_fashr0143/source_snapshots/r3_full_universe_functional_simulation_functions.R")
source(file.path(experiment, "two_stage_helpers.R"))

must_fail <- function(expression) {
  stopifnot(inherits(tryCatch(force(expression), error = identity), "error"))
}
targets <- r3ts_contract()$targets
lfdr <- c(0.01, 0.08, 0.20, 0.40, 0.80, 0.90)
fdr <- data.frame(index = 1:6, lfdr = lfdr, FDR = cumsum(lfdr) / seq_along(lfdr))
candidates <- r3ts_stage1(fdr, 6L)
stopifnot(identical(candidates, 1:2), lfdr[[2L]] > 0.05)
must_fail(r3ts_stage1(fdr[-1, ], 6L))
lfsr <- list(early = c(0.01, 0.08), middle = c(0.06, 0.10),
             late = c(0.02, 0.20), switch = c(0.01, 0.01))
lfsr <- lapply(lfsr, setNames, as.character(candidates))
selected <- r3ts_stage2(candidates, lfsr, 6L)
stopifnot(identical(selected$index[selected$target == "early" & selected$selected], 1:2))
bad <- lfsr
names(bad$early) <- c("1", "3")
must_fail(r3ts_stage2(candidates, bad, 6L))
truth <- matrix(-1, 6, 4, dimnames = list(NULL, targets))
truth[c(1, 3), "early"] <- 1
truth[3, "middle"] <- 1
truth[2, "late"] <- 1
truth[1:3, "switch"] <- 1
truth[2, "early"] <- 0 # An exact functional null is a false call.
dynamic <- c(TRUE, FALSE, TRUE, FALSE, FALSE, FALSE)
m <- r3ts_metrics(truth, dynamic, candidates, selected, 12345L, "random_bspline")
early <- m$functional[m$functional$target == "early", ]
middle <- m$functional[m$functional$target == "middle", ]
stopifnot(m$stage1$dynamic_fdr == 0.5, m$stage1$dynamic_power == 0.5,
          early$empirical_fsr == 0.5, early$empirical_power == 0.5,
          early$conditional_classification_power == 1,
          middle$calls == 0, middle$empirical_fsr == 0,
          middle$empirical_power == 0, is.na(middle$conditional_classification_power))

# Changing truth changes evaluation, but cannot change either selection interface.
flipped <- r3ts_metrics(-truth, !dynamic, candidates, selected, 12345L, "random_bspline")
stopifnot(identical(candidates, r3ts_stage1(fdr, 6L)),
          identical(selected, r3ts_stage2(candidates, lfsr, 6L)),
          flipped$stage1$dynamic_power != m$stage1$dynamic_power)

empty_lfsr <- setNames(vector("list", 4L), targets)
empty_stage2 <- r3ts_stage2(integer(), empty_lfsr, 6L)
empty <- r3ts_metrics(truth, dynamic, integer(), empty_stage2, 12345L, "random_bspline")
stopifnot(empty$stage1$dynamic_fdr == 0, empty$stage1$dynamic_power == 0,
          all(empty$functional$calls == 0), all(empty$functional$empirical_fsr == 0),
          all(is.na(empty$functional$conditional_classification_power)))

# Subset re-ranking is not intersection with calls selected in a larger universe.
whole <- functional_cfsr_table(1:3, c(0.001, 0.08, 0.07))
subset <- functional_cfsr_table(2:3, c(0.08, 0.07))
stopifnot(identical(intersect(whole$index[whole$cfsr <= 0.05], 2:3), 3L),
          length(subset$index[subset$cfsr <= 0.05]) == 0L)

stage1_rows <- list()
functional_rows <- list()
for (mechanism in r3ts_contract()$mechanisms) {
  for (i in seq_along(r3ts_contract()$seeds)) {
    seed <- r3ts_contract()$seeds[[i]]
    s <- m$stage1
    f <- m$functional
    s$seed <- f$seed <- seed
    s$truth_mechanism <- f$truth_mechanism <- mechanism
    f$empirical_fsr <- if (i == 1L) 1 else 0
    f$calls <- if (i == 1L) 1L else 100L
    f$false_calls <- if (i == 1L) 1L else 0L
    stage1_rows[[length(stage1_rows) + 1L]] <- s
    functional_rows[[length(functional_rows) + 1L]] <- f
  }
}
s <- do.call(rbind, stage1_rows)
f <- do.call(rbind, functional_rows)
mc <- r3ts_summaries(s, f)
stopifnot(all(mc$primary$mean_empirical_fsr == 0.2),
          all(mc$primary$mean_calls == 80.2),
          sum(f$false_calls) / sum(f$calls) != 0.2)
must_fail(r3ts_summaries(s[-1, ], f))
must_fail(r3ts_summaries(rbind(s, s[1, ]), f))
must_fail(r3ts_summaries(s, f[-1, ]))
cat("All two-stage selection, denominator, and five-seed aggregation tests passed.\n")
