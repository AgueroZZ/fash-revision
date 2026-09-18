#!/usr/bin/env Rscript

# Exercise production reporting with a clearly synthetic, temporary cache.
experiment <- "code/revision_simulations/internal/r3_two_stage_functional_testing"
source(file.path(experiment, "test_two_stage.R"))
source(file.path(experiment, "reporting.R"))
contract <- r3ts_contract()
n <- contract$n_units
pair_keys <- paste0("synthetic_test_unit_", seq_len(n))
truth_large <- matrix(-1, n, 4, dimnames = list(pair_keys, contract$targets))
truth_large[1:6, ] <- truth
dynamic_large <- c(dynamic, rep(FALSE, n - 6L))
lfdr_large <- c(lfdr, rep(1, n - 6L))
fdr_large <- data.frame(index = seq_len(n), lfdr = lfdr_large,
                        FDR = cumsum(lfdr_large) / seq_len(n))
records <- list()
for (mechanism in contract$mechanisms) {
  for (seed in contract$seeds) {
    metrics <- r3ts_metrics(truth_large, dynamic_large, candidates, selected, seed, mechanism)
    record <- list(schema = contract$schema, seed = seed, truth_mechanism = mechanism,
      fitted_n_units = n, selected_pair_keys = pair_keys,
      inference = list(fdr_table = fdr_large, candidates = candidates,
        lfsr_by_target = lfsr, stage2 = selected, metrics = metrics,
        true_functionals = truth_large, true_dynamic = dynamic_large))
    r3ts_validate_replicate(record)
    records[[paste0(mechanism, "_seed_", seed, ".rds")]] <- record
  }
}
wrong <- records[[1L]]
wrong$inference$metrics$functional$empirical_power[[1L]] <- 1
must_fail(r3ts_validate_replicate(wrong))

temporary_root <- tempfile("r3-two-stage-reporting-test-")
dir.create(temporary_root)
stopifnot(!r3ts_load_report(temporary_root)$ready)
frozen_relative <- "code/revision_simulations/r3_r4_fashr0143/source_snapshots/r3_full_universe_functional_simulation_functions.R"
dir.create(dirname(file.path(temporary_root, frozen_relative)), recursive = TRUE)
stopifnot(file.copy(frozen_relative, file.path(temporary_root, frozen_relative)))
directory <- file.path(temporary_root, "output/revision_simulations/internal", contract$result_id)
dir.create(file.path(directory, "replicates"), recursive = TRUE)
for (name in names(records)) r3ts_write_rds(records[[name]], file.path(directory, "replicates", name))
ordered_names <- as.vector(outer(contract$mechanisms, contract$seeds,
                                 function(m, s) paste0(m, "_seed_", s, ".rds")))
records <- records[ordered_names]
stage1 <- do.call(rbind, lapply(unname(records), function(x) x$inference$metrics$stage1))
functional <- do.call(rbind, lapply(unname(records), function(x) x$inference$metrics$functional))
summaries <- r3ts_summaries(stage1, functional)
tables <- list(stage1_by_seed = stage1, functional_by_seed = functional,
  stage1_mc_summary = summaries$stage1, functional_mc_summary = summaries$primary,
  input_validation = data.frame(fitted_n_units = rep(n, 10), stage1_calls = 2,
    original_dynamic_calls = 2, fitted_bf_pi0 = 0.8, original_bf_pi0 = 0.8))
for (mechanism in contract$mechanisms) {
  tables[[paste0("primary_", mechanism)]] <- r3ts_primary_values(summaries$primary, mechanism)
  shown <- r3ts_primary_table(summaries$primary, mechanism)
  stopifnot(identical(names(shown), c("Metric", "Early", "Middle", "Late", "Switch")),
    identical(shown$Metric, c("Empirical FSR", "empirical power", "Mean number of calls")),
    nrow(r3ts_stage1_table(stage1, mechanism)) == 6L,
    nrow(r3ts_stage2_table(summaries$primary, mechanism)) == 4L)
}
configuration <- list(contract = contract, q_dyn = 0.05, q_func = 0.05,
  package = list(version = contract$package_version, remote_sha = contract$package_sha),
  fit_universe = "all_units", prior_refit_after_screening = FALSE)
r3ts_write_rds(configuration, file.path(directory, "configuration.rds"))
r3ts_write_rds(tables, file.path(directory, "summary_tables.rds"))
artifacts <- list.files(directory, recursive = TRUE)
manifest <- list(schema = contract$schema, result_id = contract$result_id,
  configuration = configuration,
  artifact_sha256 = setNames(vapply(file.path(directory, artifacts), r3ts_sha256, character(1)), artifacts))
r3ts_write_rds(manifest, file.path(directory, "manifest.rds"))
writeLines(c(paste0("result_id=", contract$result_id), "replicates=10", "q_dyn=0.05", "q_func=0.05",
  "input_digest_checks=passed", "fit_universe=6362", "functional_candidate_scope=dynamic_fdr_screen"),
  file.path(directory, "complete.flag"))
stopifnot(r3ts_load_report(temporary_root)$ready)
tables$functional_mc_summary$mean_empirical_power[[1L]] <- 0.999
r3ts_write_rds(tables, file.path(directory, "summary_tables.rds"))
must_fail(r3ts_load_report(temporary_root))
cat("Synthetic reporting cache, exact primary-table layout, and corruption gates passed.\n")
cat("Synthetic fixture retained only in temporary test directory:", temporary_root, "\n")
