#!/usr/bin/env Rscript

# Focused checks for the ten-seed summaries and pilot-5 record migration.

project_root <- if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
  "."
} else {
  "coderepo-local"
}
experiment <- file.path(
  project_root,
  "code/revision_simulations/internal/r3_two_stage_functional_testing_pilot10"
)
base_experiment <- file.path(
  project_root,
  "code/revision_simulations/internal/r3_two_stage_functional_testing"
)
source(file.path(experiment, "pilot10_helpers.R"))
source(file.path(experiment, "reporting.R"))
base_environment <- new.env(parent = globalenv())
sys.source(file.path(base_experiment, "two_stage_helpers.R"), envir = base_environment)

expect_error <- function(expression) {
  inherits(try(force(expression), silent = TRUE), "try-error")
}

contract <- r3ts10_contract()
stopifnot(
  length(contract$seeds) == 10L,
  length(contract$added_seeds) == 5L,
  length(contract$alpha_grid) == 40L,
  identical(contract$q_dyn, 0.05),
  identical(contract$q_func, 0.05),
  length(r3ts10_expected_keys()) == 20L,
  !anyDuplicated(r3ts10_expected_keys())
)

indices <- c(8L, 2L, 5L)
lfsr <- c(0.03, 0.01, 0.02)
expected_cfsr <- data.frame(
  index = c(2L, 5L, 8L),
  lfsr = c(0.01, 0.02, 0.03),
  cfsr = c(0.01, 0.015, 0.02)
)
attr(expected_cfsr, "row.names") <- c(2L, 3L, 1L)
stopifnot(identical(r3ts10_functional_cfsr_table(indices, lfsr), expected_cfsr))

formal_dir <- file.path(
  project_root,
  "output/revision_simulations/mc",
  contract$formal_id
)
formal_configuration <- readRDS(file.path(formal_dir, "configuration.rds"))
r3ts10_validate_formal_configuration(formal_configuration)

base_result_dir <- file.path(
  project_root,
  "output/revision_simulations/internal",
  contract$base_result_id
)
base_record_path <- file.path(
  base_result_dir,
  "replicates/random_bspline_seed_12345.rds"
)
base_record <- readRDS(base_record_path)
migrated <- r3ts10_migrate_base_record(
  base_record,
  configuration_token = "fixture-token",
  source_sha256 = r3ts10_sha256(base_record_path),
  base_environment = base_environment,
  cfsr_function = r3ts10_functional_cfsr_table
)
stopifnot(
  identical(migrated$inference, base_record$inference),
  identical(migrated$selected_pair_keys, base_record$selected_pair_keys),
  identical(migrated$validation$kind, "validated_pilot5_reuse"),
  identical(migrated$configuration_token, "fixture-token")
)

actual_curve <- r3ts10_stage1_alpha_by_seed(migrated, base_environment)
actual_alpha005 <- actual_curve[abs(actual_curve$alpha - 0.05) < 1e-12, ]
stopifnot(
  nrow(actual_curve) == 40L,
  nrow(actual_alpha005) == 1L,
  identical(actual_alpha005$calls, base_record$inference$metrics$stage1$dynamic_calls),
  isTRUE(all.equal(
    actual_alpha005$empirical_fdr,
    base_record$inference$metrics$stage1$dynamic_fdr
  )),
  isTRUE(all.equal(
    actual_alpha005$power,
    base_record$inference$metrics$stage1$dynamic_power
  )),
  all(diff(actual_curve$calls) >= 0L),
  all(diff(actual_curve$power) >= -1e-12)
)

base_manifest <- readRDS(file.path(base_result_dir, "manifest.rds"))
base_replicate_names <- as.vector(outer(
  contract$mechanisms,
  contract$base_seeds,
  function(mechanism, seed) paste0(mechanism, "_seed_", seed, ".rds")
))
all_migrated <- lapply(base_replicate_names, function(file_name) {
  relative_path <- file.path("replicates", file_name)
  source_path <- file.path(base_result_dir, relative_path)
  r3ts10_migrate_base_record(
    readRDS(source_path),
    configuration_token = "all-base-fixture-token",
    source_sha256 = base_manifest$artifact_sha256[[relative_path]],
    base_environment = base_environment,
    cfsr_function = r3ts10_functional_cfsr_table
  )
})
all_migrated_keys <- vapply(
  all_migrated,
  function(record) paste(record$truth_mechanism, record$seed, sep = ":"),
  character(1)
)
all_base_curves <- do.call(rbind, lapply(
  all_migrated,
  r3ts10_stage1_alpha_by_seed,
  base_environment = base_environment
))
stopifnot(
  length(all_migrated) == 10L,
  !anyDuplicated(all_migrated_keys),
  setequal(all_migrated_keys, base_environment$r3ts_expected_keys()),
  nrow(all_base_curves) == 400L
)

stage1_fixture <- expand.grid(
  truth_mechanism = contract$mechanisms,
  seed = contract$seeds,
  alpha = contract$alpha_grid,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
seed_index <- match(stage1_fixture$seed, contract$seeds)
mechanism_index <- match(stage1_fixture$truth_mechanism, contract$mechanisms)
stage1_fixture$calls <- 100L + seed_index + round(100 * stage1_fixture$alpha)
stage1_fixture$true_calls <- stage1_fixture$calls - mechanism_index
stage1_fixture$false_calls <- mechanism_index
stage1_fixture$empirical_fdr <- stage1_fixture$false_calls / stage1_fixture$calls
stage1_fixture$power <- 0.20 + 0.01 * seed_index +
  0.50 * stage1_fixture$alpha - 0.02 * (mechanism_index - 1L)

stage2_fixture <- expand.grid(
  truth_mechanism = contract$mechanisms,
  seed = contract$seeds,
  target = contract$targets,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
seed_index2 <- match(stage2_fixture$seed, contract$seeds)
target_index <- match(stage2_fixture$target, contract$targets)
mechanism_index2 <- match(stage2_fixture$truth_mechanism, contract$mechanisms)
stage2_fixture$calls <- 20L + seed_index2 + target_index
stage2_fixture$empirical_fsr <- 0.001 * seed_index2 + 0.002 * target_index
stage2_fixture$empirical_power <- 0.10 + 0.01 * seed_index2
stage2_fixture$conditional_classification_power <-
  0.30 + 0.01 * seed_index2 + 0.02 * target_index -
  0.01 * (mechanism_index2 - 1L)

summaries <- r3ts10_summarize(stage1_fixture, stage2_fixture)
one_group <- stage1_fixture[
  stage1_fixture$truth_mechanism == "random_bspline" &
    abs(stage1_fixture$alpha - 0.05) < 1e-12,
]
saved_group <- summaries$stage1[
  summaries$stage1$truth_mechanism == "random_bspline" &
    abs(summaries$stage1$alpha - 0.05) < 1e-12,
]
stopifnot(
  nrow(summaries$stage1) == 80L,
  nrow(summaries$stage2) == 8L,
  all(summaries$stage1$n_seeds == 10L),
  all(summaries$stage2$n_seeds == 10L),
  isTRUE(all.equal(saved_group$mean_power, mean(one_group$power))),
  isTRUE(all.equal(saved_group$power_min, min(one_group$power))),
  isTRUE(all.equal(saved_group$power_max, max(one_group$power)))
)

display <- r3ts10_stage2_display_values(summaries$stage2, "random_bspline")
stopifnot(
  identical(
    display$Metric,
    c("Empirical FSR", "Conditional classification power", "Mean number of calls")
  ),
  identical(names(display), c("Metric", "Early", "Middle", "Late", "Switch")),
  !any(grepl("unconditional|end-to-end", display$Metric, ignore.case = TRUE))
)
formatted <- r3ts10_format_stage2_table(summaries$stage2, "random_bspline")
stopifnot(is.character(formatted$Early), nrow(formatted) == 3L)

plot_path <- tempfile("r3-stage1-plot-", fileext = ".png")
grDevices::png(plot_path, width = 1200, height = 800, res = 150)
r3ts10_plot_stage1(summaries$stage1, "power")
grDevices::dev.off()
stopifnot(file.exists(plot_path), file.info(plot_path)$size > 1000L)
unlink(plot_path)
plot_path <- tempfile("r3-stage1-fdr-plot-", fileext = ".png")
grDevices::png(plot_path, width = 1200, height = 800, res = 150)
r3ts10_plot_stage1(summaries$stage1, "fdr")
grDevices::dev.off()
stopifnot(file.exists(plot_path), file.info(plot_path)$size > 1000L)
unlink(plot_path)

stopifnot(
  expect_error(r3ts10_summarize(stage1_fixture[-1L, ], stage2_fixture)),
  expect_error(r3ts10_summarize(stage1_fixture, stage2_fixture[-1L, ])),
  expect_error(r3ts10_stage1_at_alpha(summaries$stage1, 0.051))
)

cat("Pilot-10 helper and migration checks passed.\n")
