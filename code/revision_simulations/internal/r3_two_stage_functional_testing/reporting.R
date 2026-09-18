# Load only completed, validated production results for the internal page.

r3ts_show_table <- function(x, caption) {
  cat('<div class="r3-two-stage-table">\n\n')
  print(knitr::kable(x, caption = caption, row.names = FALSE))
  cat('\n\n</div>\n\n')
}

r3ts_primary_table <- function(summary, mechanism) {
  out <- r3ts_primary_values(summary, mechanism)
  for (column in names(out)[-1L]) {
    out[[column]] <- c(sprintf("%.3f", out[[column]][1:2]), sprintf("%.1f", out[[column]][3]))
  }
  out
}

r3ts_stage1_table <- function(stage1, mechanism) {
  x <- stage1[stage1$truth_mechanism == mechanism, , drop = FALSE]
  x <- x[order(x$seed), , drop = FALSE]
  x <- rbind(x, transform(x[1, , drop = FALSE],
    seed = NA_integer_, dynamic_calls = mean(x$dynamic_calls),
    dynamic_fdr = mean(x$dynamic_fdr), dynamic_power = mean(x$dynamic_power)))
  data.frame(
    Seed = ifelse(is.na(x$seed), "Mean", as.character(x$seed)),
    `Dynamic discoveries` = sprintf("%.1f", x$dynamic_calls),
    `Empirical dynamic FDR` = sprintf("%.3f", x$dynamic_fdr),
    `Dynamic power` = sprintf("%.3f", x$dynamic_power), check.names = FALSE
  )
}

r3ts_stage2_table <- function(summary, mechanism) {
  x <- summary[summary$truth_mechanism == mechanism, , drop = FALSE]
  x <- x[match(r3ts_contract()$targets, x$target), , drop = FALSE]
  data.frame(
    Target = tools::toTitleCase(x$target),
    `Mean calls` = sprintf("%.1f", x$mean_calls),
    `Calls range` = paste(x$min_calls, x$max_calls, sep = " to "),
    `Calls SD` = sprintf("%.1f", x$sd_calls),
    `Conditional classification power` = ifelse(
      is.na(x$mean_conditional_classification_power), "NA",
      sprintf("%.3f", x$mean_conditional_classification_power)),
    `Valid seeds` = x$conditional_valid_seeds, check.names = FALSE
  )
}

r3ts_load_report <- function(project_root) {
  contract <- r3ts_contract()
  directory <- file.path(project_root, "output/revision_simulations/internal", contract$result_id)
  if (!dir.exists(directory)) {
    return(list(ready = FALSE, partial = dir.exists(paste0(directory, ".partial")),
                directory = directory, contract = contract))
  }
  r3ts_assert(file.exists(file.path(directory, "complete.flag")),
               "An incomplete production cache cannot be reported.")
  flag <- readLines(file.path(directory, "complete.flag"), warn = FALSE)
  r3ts_assert(all(c(paste0("result_id=", contract$result_id), "replicates=10",
    "q_dyn=0.05", "q_func=0.05", "input_digest_checks=passed", "fit_universe=6362",
    "functional_candidate_scope=dynamic_fdr_screen") %in% flag),
    "The completion marker does not satisfy the production contract.")
  manifest <- readRDS(file.path(directory, "manifest.rds"))
  r3ts_assert(identical(manifest$schema, contract$schema) &&
                 identical(manifest$result_id, contract$result_id), "Incorrect internal result cache.")
  configuration <- readRDS(file.path(directory, "configuration.rds"))
  r3ts_assert(identical(configuration, manifest$configuration) &&
                 identical(configuration$contract, contract) &&
                 configuration$q_dyn == 0.05 && configuration$q_func == 0.05 &&
                 identical(configuration$package$version, contract$package_version) &&
                 identical(configuration$package$remote_sha, contract$package_sha) &&
                 identical(configuration$fit_universe, "all_units") &&
                 identical(configuration$prior_refit_after_screening, FALSE),
               "The cache does not implement the fixed two-stage contract.")
  actual <- vapply(file.path(directory, names(manifest$artifact_sha256)), r3ts_sha256, character(1))
  r3ts_assert(identical(unname(actual), unname(manifest$artifact_sha256)),
               "Internal result artifact hash mismatch.")
  frozen <- file.path(project_root,
    "code/revision_simulations/r3_r4_fashr0143/source_snapshots/r3_full_universe_functional_simulation_functions.R")
  r3ts_assert(identical(r3ts_sha256(frozen), contract$source_sha[["simulation"]]),
               "The reporting selection helper is not the frozen R3 version.")
  selection_environment <- new.env(parent = environment())
  sys.source(frozen, envir = selection_environment)
  paths <- file.path(directory, "replicates", as.vector(outer(contract$mechanisms, contract$seeds,
    function(m, s) paste0(m, "_seed_", s, ".rds"))))
  records <- lapply(paths, readRDS)
  invisible(lapply(records, r3ts_validate_replicate,
                   cfsr_function = selection_environment$functional_cfsr_table))
  tables <- readRDS(file.path(directory, "summary_tables.rds"))
  r3ts_assert(isTRUE(all.equal(tables$stage1_by_seed,
    do.call(rbind, lapply(records, function(x) x$inference$metrics$stage1)))) &&
    isTRUE(all.equal(tables$functional_by_seed,
    do.call(rbind, lapply(records, function(x) x$inference$metrics$functional)))),
    "Summary tables disagree with auditable replicate records.")
  computed <- r3ts_summaries(tables$stage1_by_seed, tables$functional_by_seed)
  r3ts_assert(isTRUE(all.equal(computed$primary, tables$functional_mc_summary)) &&
                 isTRUE(all.equal(computed$stage1, tables$stage1_mc_summary)),
               "Reported means disagree with replicate-specific metrics.")
  for (mechanism in contract$mechanisms) {
    r3ts_assert(identical(tables[[paste0("primary_", mechanism)]],
      r3ts_primary_values(computed$primary, mechanism)), "Primary CSV table mismatch.")
  }
  validation <- tables$input_validation
  r3ts_assert(nrow(validation) == 10L && all(validation$fitted_n_units == contract$n_units) &&
                 all(validation$stage1_calls == validation$original_dynamic_calls) &&
                 all(abs(validation$fitted_bf_pi0 - validation$original_bf_pi0) < 1e-9),
               "Full-data fit reconstruction checks did not pass.")
  list(ready = TRUE, directory = directory, contract = contract,
       manifest = manifest, configuration = configuration, tables = tables)
}

r3ts_flag_text <- function(summary, mechanism) {
  x <- summary[summary$truth_mechanism == mechanism, , drop = FALSE]
  describe <- function(flag, label) {
    targets <- tools::toTitleCase(x$target[flag])
    if (!length(targets)) return(character())
    paste0(label, ": ", paste(targets, collapse = ", "), ".")
  }
  flags <- c(describe(x$fsr_above_nominal, "Mean empirical FSR above 0.05"),
             describe(x$small_call_count, "Fewer than 10 calls in at least one seed"),
             describe(x$large_power_gap, "Conditional minus end-to-end power above 0.10"))
  if (!length(flags)) "No prespecified descriptive flag was triggered." else paste(flags, collapse = " ")
}
