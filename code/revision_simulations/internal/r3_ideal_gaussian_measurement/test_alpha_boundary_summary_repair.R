# Tests for the R3 ideal-Gaussian alpha-boundary summary repair.

find_workflowr_root <- function() {
  candidates <- c(".", "coderepo-local", "../..", "../../..", "../../../..")
  marker <- file.path(
    candidates,
    "code", "revision_simulations", "shared", "simulation_functions.R"
  )
  hits <- candidates[file.exists(marker)]
  if (length(hits) == 0L) {
    stop("Could not locate the workflowr project root.", call. = FALSE)
  }
  normalizePath(hits[[1L]], winslash = "/", mustWork = TRUE)
}

project_root <- find_workflowr_root()
repair_path <- file.path(
  project_root,
  "code", "revision_simulations", "internal",
  "r3_ideal_gaussian_measurement", "repair_alpha_boundary_summary.R"
)
source(repair_path)

alpha_near_005 <- 0.05 - 5e-17
stopifnot(
  alpha_near_005 < 0.05,
  isTRUE(select_alpha_interval(alpha_near_005, 0.05, 0.20, 1e-12)),
  !isTRUE(select_alpha_interval(alpha_near_005, 0.05, 0.20, 1e-18)),
  identical(canonicalize_alpha(alpha_near_005), 0.05)
)

mechanisms <- c("random_bspline", "raised_cosine")
seeds <- c(101L, 202L)
alpha_grid <- c(alpha_near_005, 0.055, 0.20)
synthetic <- do.call(rbind, lapply(mechanisms, function(mechanism) {
  do.call(rbind, lapply(seeds, function(seed) {
    data.frame(
      method = "FASH-IWP1-BF",
      target = "middle",
      alpha = alpha_grid,
      empirical_fsr = if (mechanism == "random_bspline") {
        c(0.04, 0.06, 0.31) + (seed == seeds[[2L]]) * 0.01
      } else {
        c(0.03, 0.04, 0.18) + (seed == seeds[[2L]]) * 0.01
      },
      power = c(0.20, 0.25, 0.50),
      seed = seed,
      truth_mechanism = mechanism,
      stringsAsFactors = FALSE
    )
  }))
}))
irrelevant <- synthetic
irrelevant$method <- "FASH-IWP1-Raw"
irrelevant$target <- "early"
synthetic <- rbind(synthetic, irrelevant)

rebuilt <- rebuild_middle_summaries(
  all_alpha = synthetic,
  mechanisms = mechanisms,
  expected_seeds = seeds,
  lower = 0.05,
  upper = 0.20,
  tolerance = 1e-12
)
stopifnot(
  nrow(rebuilt$curve) == length(mechanisms) * length(alpha_grid),
  nrow(rebuilt$primary) == length(mechanisms),
  identical(sort(unique(rebuilt$curve$alpha)), c(0.05, 0.055, 0.20)),
  min(rebuilt$curve$alpha) == 0.05,
  all(rebuilt$primary$alpha_min == 0.05),
  rebuilt$primary$alpha_at_maximum[
    rebuilt$primary$truth_mechanism == "random_bspline"
  ] == 0.20
)

temporary_csv <- tempfile(fileext = ".csv")
on.exit(unlink(temporary_csv), add = TRUE)
utils::write.csv(rebuilt$curve, temporary_csv, row.names = FALSE)
round_trip <- utils::read.csv(temporary_csv, stringsAsFactors = FALSE)
stopifnot(
  nrow(round_trip) == nrow(rebuilt$curve),
  min(round_trip$alpha) == 0.05,
  identical(sort(unique(round_trip$alpha)), c(0.05, 0.055, 0.20))
)

message("Alpha-boundary summary repair tests passed.")
