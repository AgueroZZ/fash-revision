# Selection and evaluation for the internal two-stage R3 experiment.

r3ts_contract <- function() {
  list(
    schema = "r3-two-stage-functional-testing-v1",
    result_id = paste0(
      "r3_two_stage_functional_testing_iwp1_mixture_",
      "qdyn005_qfunc005_fashr0143_pilot5"
    ),
    formal_id = paste0(
      "r3_real_genotype_one_per_gene_J6362_matched_functional_",
      "open_middle_3_12_center_aligned_iwp1_geometry_mixture_",
      "relative_location_clearance_full_universe_paired_posterior_",
      "fashr0143_pilot5"
    ),
    ideal_id = paste0(
      "r3_ideal_gaussian_known_t_adjusted_se_matched_truth_",
      "open_middle_3_12_center_aligned_iwp1_geometry_mixture_",
      "full_universe_fashr0143_pilot5"
    ),
    seeds = c(12345L, 22345L, 32345L, 42345L, 52345L),
    mechanisms = c("random_bspline", "raised_cosine"),
    targets = c("early", "middle", "late", "switch"),
    q_dyn = 0.05, q_func = 0.05, n_units = 6362L,
    posterior_draws = 3000L,
    package_version = "0.1.43",
    package_sha = "bf223df75da6e41ae48607a56b4cd12d7c3b24e7",
    source_sha = c(
      simulation = "45267b0884168e5ae33cc4f14e3f05b711d961b65bf9c6fbd880e748de064a6e",
      genotype_helper = "c03c01a188503336a77793c96f4e2d3ac7e0cbd56f4028b552da4b2f88e6b9d7",
      digest_helper = "84bc06be91531f1587b0f298bea82f0bd937444663419c0b77d589b48d3fe84e",
      genotypes = "81bbef5f323a0bab2ca993c782d8a9b7c63518b83c2cdb46ef7ed1d46f65af24"
    )
  )
}

r3ts_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

r3ts_sha256 <- function(path) {
  r3ts_assert(file.exists(path) && !dir.exists(path), paste("Missing file:", path))
  mac <- identical(Sys.info()[["sysname"]], "Darwin")
  output <- system2(
    if (mac) "shasum" else "sha256sum",
    c(if (mac) c("-a", "256"), shQuote(path)), stdout = TRUE
  )
  r3ts_assert(is.null(attr(output, "status")), paste("Hash failed:", path))
  sub("[[:space:]].*$", "", output[[1L]])
}

r3ts_write_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3)
  r3ts_assert(file.rename(temporary, path), paste("Atomic write failed:", path))
}

# Truth is deliberately absent from both selection interfaces.
r3ts_stage1 <- function(fdr_table, n_units, q = 0.05) {
  r3ts_assert(all(c("index", "FDR") %in% names(fdr_table)), "Incomplete FDR table.")
  r3ts_assert(
    nrow(fdr_table) == n_units && !anyDuplicated(fdr_table$index) &&
      setequal(fdr_table$index, seq_len(n_units)),
    "Stage 1 requires the full-data FDR table."
  )
  r3ts_assert(all(is.finite(fdr_table$FDR) & fdr_table$FDR >= 0 &
                     fdr_table$FDR <= 1), "Invalid cumulative FDR values.")
  sort(as.integer(fdr_table$index[fdr_table$FDR <= q]))
}

r3ts_stage2 <- function(candidates, lfsr_by_target, n_units, q = 0.05,
                        cfsr_function = functional_cfsr_table) {
  targets <- r3ts_contract()$targets
  candidates <- as.integer(candidates)
  r3ts_assert(!anyNA(candidates) && !anyDuplicated(candidates) &&
                 all(candidates %in% seq_len(n_units)), "Invalid Stage-1 candidates.")
  r3ts_assert(identical(names(lfsr_by_target), targets), "Incorrect functional targets.")
  rows <- lapply(targets, function(target) {
    values <- lfsr_by_target[[target]]
    r3ts_assert(length(values) == length(candidates) &&
                   !anyDuplicated(names(values)) &&
                   setequal(names(values), as.character(candidates)),
                 "Stage-2 candidates must exactly equal Stage-1 discoveries.")
    values <- unname(values[as.character(candidates)])
    tab <- cfsr_function(candidates, values)
    data.frame(target = rep(target, nrow(tab)), tab,
               selected = tab$cfsr <= q, stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

r3ts_metrics <- function(true_functionals, true_dynamic, candidates,
                         stage2, seed, mechanism) {
  targets <- r3ts_contract()$targets
  truth <- as.matrix(true_functionals)
  n <- nrow(truth)
  r3ts_assert(identical(colnames(truth), targets) && all(is.finite(truth)),
               "Incomplete or non-finite functional truth.")
  r3ts_assert(is.logical(true_dynamic) && length(true_dynamic) == n &&
                 !anyNA(true_dynamic), "Invalid dynamic truth.")
  r3ts_assert(!anyDuplicated(candidates) && all(candidates %in% seq_len(n)),
               "Invalid candidates in metric evaluation.")
  r3ts_assert(all(stage2$target %in% targets), "Unexpected Stage-2 target.")
  dynamic_calls <- length(candidates)
  dynamic_tp <- sum(true_dynamic[candidates])
  dynamic_total <- sum(true_dynamic)
  stage1 <- data.frame(
    seed = seed, truth_mechanism = mechanism, n_units = n,
    dynamic_calls = dynamic_calls, dynamic_true_positives = dynamic_tp,
    dynamic_false_calls = dynamic_calls - dynamic_tp,
    dynamic_true_total = dynamic_total,
    dynamic_fdr = (dynamic_calls - dynamic_tp) / max(dynamic_calls, 1L),
    dynamic_power = if (dynamic_total) dynamic_tp / dynamic_total else NA_real_
  )
  functional <- do.call(rbind, lapply(targets, function(target) {
    tab <- stage2[stage2$target == target, , drop = FALSE]
    r3ts_assert(nrow(tab) == length(candidates) && !anyDuplicated(tab$index) &&
                   setequal(tab$index, candidates), "Stage-2 candidate mismatch.")
    selected <- tab$index[tab$selected]
    positive <- truth[, target] > 0
    total <- sum(positive)
    retained <- sum(positive[candidates])
    tp <- sum(positive[selected])
    calls <- length(selected)
    data.frame(
      seed = seed, truth_mechanism = mechanism, target = target,
      candidate_count = length(candidates), calls = calls,
      true_calls = tp, false_calls = calls - tp,
      true_target_total = total, true_target_in_stage1 = retained,
      empirical_fsr = (calls - tp) / max(calls, 1L),
      empirical_power = if (total) tp / total else NA_real_,
      conditional_classification_power = if (retained) tp / retained else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  list(stage1 = stage1, functional = functional)
}

r3ts_expected_keys <- function() {
  contract <- r3ts_contract()
  as.vector(outer(contract$mechanisms, contract$seeds, paste, sep = ":"))
}

r3ts_validate_replicate <- function(record, cfsr_function = functional_cfsr_table) {
  contract <- r3ts_contract()
  r3ts_assert(identical(record$schema, contract$schema) &&
                 record$seed %in% contract$seeds &&
                 record$truth_mechanism %in% contract$mechanisms &&
                 record$fitted_n_units == contract$n_units &&
                 length(record$selected_pair_keys) == contract$n_units,
               "Invalid production replicate identity or fit universe.")
  x <- record$inference
  r3ts_assert(nrow(x$true_functionals) == contract$n_units &&
                 identical(rownames(x$true_functionals), record$selected_pair_keys),
               "Per-unit functional truth does not match the fitted pair order.")
  candidates <- r3ts_stage1(x$fdr_table, contract$n_units, contract$q_dyn)
  stage2 <- r3ts_stage2(candidates, x$lfsr_by_target, contract$n_units, contract$q_func,
                        cfsr_function = cfsr_function)
  metrics <- r3ts_metrics(x$true_functionals, x$true_dynamic, candidates, stage2,
                          record$seed, record$truth_mechanism)
  r3ts_assert(identical(candidates, x$candidates) &&
                 isTRUE(all.equal(stage2, x$stage2)) &&
                 isTRUE(all.equal(metrics, x$metrics)),
               "Saved selections or metric denominators cannot be reproduced.")
  invisible(TRUE)
}

r3ts_primary_values <- function(summary, mechanism) {
  targets <- r3ts_contract()$targets
  x <- summary[summary$truth_mechanism == mechanism, , drop = FALSE]
  x <- x[match(targets, x$target), , drop = FALSE]
  r3ts_assert(nrow(x) == 4L && !anyNA(x$target), "Missing primary-table targets.")
  out <- data.frame(Metric = c("Empirical FSR", "empirical power", "Mean number of calls"))
  for (target in targets) {
    row <- x[x$target == target, ]
    out[[tools::toTitleCase(target)]] <- c(
      row$mean_empirical_fsr, row$mean_empirical_power, row$mean_calls
    )
  }
  out
}

r3ts_summaries <- function(stage1, functional) {
  contract <- r3ts_contract()
  keys <- paste(stage1$truth_mechanism, stage1$seed, sep = ":")
  r3ts_assert(!anyDuplicated(keys) && setequal(keys, r3ts_expected_keys()),
               "Reporting requires all five seeds for both mechanisms.")
  fkeys <- paste(functional$truth_mechanism, functional$seed, functional$target, sep = ":")
  expected <- as.vector(outer(r3ts_expected_keys(), contract$targets, paste, sep = ":"))
  r3ts_assert(!anyDuplicated(fkeys) && setequal(fkeys, expected),
               "Reporting requires all four targets for all ten replicates.")
  summary_rows <- lapply(contract$mechanisms, function(mechanism) {
    do.call(rbind, lapply(contract$targets, function(target) {
      x <- functional[functional$truth_mechanism == mechanism &
                        functional$target == target, , drop = FALSE]
      cp <- x$conditional_classification_power
      r3ts_assert(all(is.finite(x$empirical_fsr)) &&
                     all(is.finite(x$empirical_power)), "Undefined primary metric.")
      cp_mean <- if (all(is.na(cp))) NA_real_ else mean(cp, na.rm = TRUE)
      data.frame(
        truth_mechanism = mechanism, target = target, n_seeds = nrow(x),
        mean_empirical_fsr = mean(x$empirical_fsr),
        mean_empirical_power = mean(x$empirical_power), mean_calls = mean(x$calls),
        min_calls = min(x$calls), max_calls = max(x$calls), sd_calls = sd(x$calls),
        mean_conditional_classification_power = cp_mean,
        conditional_valid_seeds = sum(!is.na(cp)),
        min_conditional_classification_power = if (all(is.na(cp))) NA_real_ else min(cp, na.rm = TRUE),
        max_conditional_classification_power = if (all(is.na(cp))) NA_real_ else max(cp, na.rm = TRUE),
        mean_power_gap = cp_mean - mean(x$empirical_power),
        fsr_above_nominal = mean(x$empirical_fsr) > contract$q_func,
        small_call_count = any(x$calls < 10L),
        large_power_gap = is.finite(cp_mean) && cp_mean - mean(x$empirical_power) > 0.10
      )
    }))
  })
  stage1_summary <- do.call(rbind, lapply(contract$mechanisms, function(mechanism) {
    x <- stage1[stage1$truth_mechanism == mechanism, , drop = FALSE]
    do.call(rbind, lapply(c("dynamic_calls", "dynamic_fdr", "dynamic_power"), function(metric) {
      values <- x[[metric]]
      data.frame(truth_mechanism = mechanism, metric = metric,
                 mean = mean(values), min = min(values), max = max(values), sd = sd(values))
    }))
  }))
  list(primary = do.call(rbind, summary_rows), stage1 = stage1_summary)
}
