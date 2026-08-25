#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/internal/evaluation_grid_sensitivity/middle_open_window_helpers.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/internal/evaluation_grid_sensitivity/middle_open_window_helpers.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "evaluation_grid_sensitivity",
  "middle_open_window_helpers.R"
))

for (step in c(0.15, 0.10, 0.05)) {
  evaluation_grid <- seq(0, 15, by = step)
  inside <- middle_open_window_mask(evaluation_grid)
  early <- evaluation_grid <= 3 + sqrt(.Machine$double.eps)
  late <- evaluation_grid >= 12 - sqrt(.Machine$double.eps)

  stopifnot(
    any(inside),
    any(!inside),
    !any(inside & early),
    !any(inside & late),
    !inside[which.min(abs(evaluation_grid - 3))],
    !inside[which.min(abs(evaluation_grid - 12))]
  )
}

test_lfsr <- c(0.04, 0.01, 0.01, 0.30)
test_pair_id <- c("pair_b", "pair_c", "pair_a", "pair_d")
test_cfsr <- deterministic_cumulative_fsr(test_lfsr, test_pair_id)
test_order <- order(test_lfsr, test_pair_id)
stopifnot(
  isTRUE(all.equal(
    test_cfsr[test_order],
    cumsum(test_lfsr[test_order]) / seq_along(test_order)
  )),
  set_jaccard(character(), character()) == 1,
  set_jaccard(c("a", "b"), c("b", "c")) == 1 / 3
)

cat("Middle open-window helper tests passed.\n")
