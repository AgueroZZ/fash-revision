# Ten-seed contracts and summaries for the formal two-stage R3 replacement.

r3ts10_contract <- function() {
  list(
    schema = "r3-two-stage-functional-testing-v2-pilot10",
    result_id = paste0(
      "r3_two_stage_functional_testing_iwp1_mixture_alpha_grid_",
      "qdyn005_qfunc005_fashr0143_pilot10"
    ),
    genotype_id = "real_genotype_one_per_gene_J6362_pilot10",
    base_result_id = paste0(
      "r3_two_stage_functional_testing_iwp1_mixture_",
      "qdyn005_qfunc005_fashr0143_pilot5"
    ),
    formal_id = paste0(
      "r3_real_genotype_one_per_gene_J6362_matched_functional_",
      "open_middle_3_12_center_aligned_iwp1_geometry_mixture_",
      "relative_location_clearance_full_universe_paired_posterior_",
      "fashr0143_pilot5"
    ),
    base_seeds = c(12345L, 22345L, 32345L, 42345L, 52345L),
    added_seeds = c(62345L, 72345L, 82345L, 92345L, 102345L),
    seeds = c(
      12345L, 22345L, 32345L, 42345L, 52345L,
      62345L, 72345L, 82345L, 92345L, 102345L
    ),
    mechanisms = c("random_bspline", "raised_cosine"),
    mechanism_labels = c(
      random_bspline = "Broad random B-spline",
      raised_cosine = "Compact raised cosine"
    ),
    targets = c("early", "middle", "late", "switch"),
    alpha_grid = seq(0.005, 0.20, by = 0.005),
    q_dyn = 0.05,
    q_func = 0.05,
    n_units = 6362L,
    posterior_draws = 3000L,
    package_version = "0.1.43",
    package_sha = "bf223df75da6e41ae48607a56b4cd12d7c3b24e7",
    base_source_sha = c(
      simulation = "45267b0884168e5ae33cc4f14e3f05b711d961b65bf9c6fbd880e748de064a6e",
      genotype_helper = "c03c01a188503336a77793c96f4e2d3ac7e0cbd56f4028b552da4b2f88e6b9d7",
      digest_helper = "84bc06be91531f1587b0f298bea82f0bd937444663419c0b77d589b48d3fe84e",
      base_genotypes = "81bbef5f323a0bab2ca993c782d8a9b7c63518b83c2cdb46ef7ed1d46f65af24",
      base_two_stage_helper = "d46d3a4686a58c9f81958fb060e012fbd3c2b2ac59fa35c7198d1b5935070793",
      base_reconstruction = "c8a185fa1cf7a36c5d330119bf56087c39b5ae877097ebf8da6130ea22683835",
      base_result_manifest = "b1a7bb3306f9a972004ccf4f88869541b66cb1546c6908e03397b9d6871785ae",
      formal_result_manifest = "a7e296eef2f2ca24e9a802e8f78e4c1d8e629b0585460068fb7fed657d2ce99b"
    )
  )
}

r3ts10_assert <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

r3ts10_sha256 <- function(path) {
  r3ts10_assert(
    length(path) == 1L && file.exists(path) && !dir.exists(path),
    paste("Missing file:", path)
  )
  mac <- identical(Sys.info()[["sysname"]], "Darwin")
  output <- system2(
    if (mac) "shasum" else "sha256sum",
    c(if (mac) c("-a", "256"), shQuote(path)),
    stdout = TRUE
  )
  r3ts10_assert(is.null(attr(output, "status")), paste("Hash failed:", path))
  sub("[[:space:]].*$", "", output[[1L]])
}

r3ts10_write_rds <- function(object, path, version = 3L) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  saveRDS(object, temporary, version = version)
  r3ts10_assert(file.rename(temporary, path), paste("Atomic write failed:", path))
  invisible(path)
}

r3ts10_serialized_object_md5 <- function(object) {
  path <- tempfile("r3-two-stage-pilot10-md5-", fileext = ".rds")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  saveRDS(object, path, version = 3)
  checksum <- unname(as.character(tools::md5sum(path))[[1L]])
  r3ts10_assert(!is.na(checksum), "Could not compute a serialized-object MD5.")
  checksum
}

r3ts10_functional_cfsr_table <- function(indices, lfsr) {
  indices <- as.integer(indices)
  lfsr <- as.numeric(lfsr)
  r3ts10_assert(
    length(indices) == length(lfsr) && all(is.finite(lfsr)) &&
      all(lfsr >= 0) && all(lfsr <= 1),
    "indices and lfsr must have the same length, with lfsr in [0, 1]."
  )
  if (!length(indices)) {
    return(data.frame(index = integer(), lfsr = numeric(), cfsr = numeric()))
  }
  out <- data.frame(index = indices, lfsr = lfsr)
  out <- out[order(out$lfsr, out$index), , drop = FALSE]
  out$cfsr <- cumsum(out$lfsr) / seq_len(nrow(out))
  out
}

r3ts10_expected_keys <- function() {
  contract <- r3ts10_contract()
  as.vector(outer(contract$mechanisms, contract$seeds, paste, sep = ":"))
}

r3ts10_replicate_names <- function() {
  contract <- r3ts10_contract()
  as.vector(outer(
    contract$mechanisms,
    contract$seeds,
    function(mechanism, seed) paste0(mechanism, "_seed_", seed, ".rds")
  ))
}

r3ts10_require_base_environment <- function(base_environment) {
  required <- c(
    "r3ts_stage1", "r3ts_stage2", "r3ts_metrics",
    "r3ts_validate_replicate", "r3ts_write_rds", "r3ts_sha256"
  )
  r3ts10_assert(
    is.environment(base_environment) &&
      all(vapply(required, exists, logical(1), envir = base_environment, inherits = FALSE)),
    "The validated pilot-5 helper environment is incomplete."
  )
}

r3ts10_validate_formal_configuration <- function(config) {
  contract <- r3ts10_contract()
  expected <- list(
    output_id = contract$formal_id,
    J = 6362L,
    n_donors = 19L,
    n_covariates = 5L,
    time_grid = 0:15,
    evaluation_grid = seq(0, 15, by = 0.1),
    middle_window = c(3, 12),
    middle_boundary = "open",
    switch_threshold = 0.25,
    seed_list = contract$base_seeds,
    truth_mechanisms = contract$mechanisms,
    class_probs = c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40),
    expected_class_counts = c(dynamic_bspline = 1272L, constant = 2545L, zero = 2545L),
    temporal_category_probs = c(early = 0.29, middle = 0.42, late = 0.29),
    location_truth_margin = 0.10,
    location_truth_min_range_fraction = 0.10,
    switch_truth_margin = 0.10,
    non_switch_min_abs = 0.10,
    non_switch_min_range_fraction = 0.10,
    dynamic_main_effect_sd = 1,
    expression_noise_sd = 1,
    covariate_effect_sd = 0.5,
    intercept_sd = 0,
    num_basis = 20L
  )
  for (field in names(expected)) {
    r3ts10_assert(
      isTRUE(all.equal(config[[field]], expected[[field]])),
      paste("Unexpected frozen R3 setting:", field)
    )
  }
  r3ts10_assert(
    identical(
      unname(config$expected_truth_group_counts),
      c(185L, 184L, 267L, 267L, 185L, 184L)
    ),
    "The frozen R3 functional-cell counts changed."
  )
  r3ts10_assert(
    identical(
      config$raised_cosine$center_ranges,
      list(early = c(1.5, 2.5), middle = c(4.5, 10.5), late = c(12.5, 13.5))
    ) &&
      identical(
        config$random_bspline,
        list(amplitude = 2, df = 6, coefficient_sd = 1)
      ),
    "The frozen R3 truth geometry changed."
  )
  invisible(TRUE)
}

r3ts10_validate_added_inputs <- function(inputs, genotype_sample, config) {
  contract <- r3ts10_contract()
  r3ts10_assert(
    inputs$seed %in% contract$added_seeds &&
      inputs$truth_mechanism %in% contract$mechanisms,
    "An added-seed input has the wrong identity."
  )
  effects <- inputs$effects
  dynamic <- effects$unit_info$effect_class == "dynamic_bspline"
  class_counts <- table(factor(
    effects$unit_info$effect_class,
    levels = names(config$expected_class_counts)
  ))
  cell_counts <- table(factor(
    effects$unit_info$truth_group[dynamic],
    levels = names(config$expected_truth_group_counts)
  ))
  r3ts10_assert(
    identical(as.integer(class_counts), unname(config$expected_class_counts)) &&
      identical(as.integer(cell_counts), unname(config$expected_truth_group_counts)),
    "An added seed does not reproduce the fixed truth counts."
  )
  pair_keys <- genotype_sample$selection$pair_key
  r3ts10_assert(
    identical(effects$unit_info$variant_id, pair_keys) &&
      identical(rownames(effects$true_functionals), pair_keys) &&
      identical(colnames(effects$true_functionals), contract$targets),
    "Added-seed truth is not aligned to the sampled genotype pairs."
  )
  regression <- inputs$regression
  r3ts10_assert(
    identical(dim(regression$beta_hat), c(contract$n_units, 16L)) &&
      identical(dim(regression$se), c(contract$n_units, 16L)) &&
      all(is.finite(regression$beta_hat)) &&
      all(is.finite(regression$se)) &&
      all(regression$se > 0),
    "Added-seed regression summaries are incomplete."
  )
  expected_digests <- c(
    genotype_content_md5 = genotype_content_md5(
      pair_keys, rownames(genotype_sample$G), genotype_sample$G
    ),
    true_beta_md5 = r3ts10_serialized_object_md5(effects$beta_matrix),
    adjusted_se_md5 = r3ts10_serialized_object_md5(regression$se),
    regression_beta_hat_md5 = r3ts10_serialized_object_md5(regression$beta_hat)
  )
  r3ts10_assert(
    identical(inputs$input_digests, expected_digests),
    "Added-seed input digests cannot be reproduced."
  )
  invisible(TRUE)
}

r3ts10_validate_record <- function(record, base_environment, cfsr_function) {
  r3ts10_require_base_environment(base_environment)
  contract <- r3ts10_contract()
  r3ts10_assert(
    identical(record$schema, contract$schema) &&
      is.character(record$configuration_token) && length(record$configuration_token) == 1L &&
      nzchar(record$configuration_token) &&
      record$seed %in% contract$seeds &&
      record$truth_mechanism %in% contract$mechanisms &&
      record$fitted_n_units == contract$n_units &&
      length(record$selected_pair_keys) == contract$n_units &&
      !anyDuplicated(record$selected_pair_keys),
    "Invalid pilot-10 replicate identity or fit universe."
  )
  inference <- record$inference
  r3ts10_assert(
    nrow(inference$true_functionals) == contract$n_units &&
      identical(colnames(inference$true_functionals), contract$targets) &&
      identical(rownames(inference$true_functionals), record$selected_pair_keys),
    "Pilot-10 functional truth is incomplete or misaligned."
  )
  candidates <- base_environment$r3ts_stage1(
    inference$fdr_table,
    contract$n_units,
    contract$q_dyn
  )
  stage2 <- base_environment$r3ts_stage2(
    candidates,
    inference$lfsr_by_target,
    contract$n_units,
    contract$q_func,
    cfsr_function = cfsr_function
  )
  metrics <- base_environment$r3ts_metrics(
    inference$true_functionals,
    inference$true_dynamic,
    candidates,
    stage2,
    record$seed,
    record$truth_mechanism
  )
  r3ts10_assert(
    identical(candidates, inference$candidates) &&
      isTRUE(all.equal(stage2, inference$stage2)) &&
      isTRUE(all.equal(metrics, inference$metrics)),
    "Pilot-10 selections or metrics cannot be reproduced."
  )
  r3ts10_assert(
    is.list(record$validation) &&
      identical(record$validation$kind, if (record$seed %in% contract$base_seeds) {
        "validated_pilot5_reuse"
      } else {
        "new_seed_frozen_generation"
      }),
    "Pilot-10 replicate provenance is incomplete."
  )
  if (record$seed %in% contract$base_seeds) {
    r3ts10_assert(
      identical(record$validation$base_result_id, contract$base_result_id) &&
        is.character(record$validation$source_replicate_sha256) &&
        length(record$validation$source_replicate_sha256) == 1L &&
        nzchar(record$validation$source_replicate_sha256),
      "A migrated pilot-5 replicate lacks source provenance."
    )
  } else {
    r3ts10_assert(
      identical(record$validation$input_validation, "passed") &&
        identical(record$validation$fit_universe, contract$n_units) &&
        identical(record$validation$functional_candidate_scope, "dynamic_fdr_screen"),
      "A new-seed replicate lacks generation or fit validation."
    )
  }
  invisible(TRUE)
}

r3ts10_migrate_base_record <- function(record, configuration_token,
                                        source_sha256, base_environment,
                                        cfsr_function) {
  r3ts10_require_base_environment(base_environment)
  contract <- r3ts10_contract()
  base_environment$r3ts_validate_replicate(
    record,
    cfsr_function = cfsr_function
  )
  r3ts10_assert(
    record$seed %in% contract$base_seeds &&
      record$truth_mechanism %in% contract$mechanisms,
    "Only the validated pilot-5 seeds can be migrated."
  )
  migrated <- record
  migrated$schema <- contract$schema
  migrated$configuration_token <- configuration_token
  migrated$reuse <- c(migrated$reuse, validated_pilot5_replicate = TRUE)
  migrated$validation <- list(
    kind = "validated_pilot5_reuse",
    base_result_id = contract$base_result_id,
    source_replicate_sha256 = unname(source_sha256)
  )
  r3ts10_validate_record(
    migrated,
    base_environment = base_environment,
    cfsr_function = cfsr_function
  )
  migrated
}

r3ts10_stage1_alpha_by_seed <- function(record, base_environment) {
  r3ts10_require_base_environment(base_environment)
  contract <- r3ts10_contract()
  truth <- record$inference$true_dynamic
  r3ts10_assert(
    is.logical(truth) && length(truth) == contract$n_units && !anyNA(truth),
    "Dynamic truth is invalid for the Stage-1 alpha curve."
  )
  do.call(rbind, lapply(contract$alpha_grid, function(alpha) {
    selected <- base_environment$r3ts_stage1(
      record$inference$fdr_table,
      contract$n_units,
      alpha
    )
    calls <- length(selected)
    true_calls <- sum(truth[selected])
    data.frame(
      seed = record$seed,
      truth_mechanism = record$truth_mechanism,
      alpha = alpha,
      calls = calls,
      true_calls = true_calls,
      false_calls = calls - true_calls,
      empirical_fdr = (calls - true_calls) / max(calls, 1L),
      power = true_calls / sum(truth),
      stringsAsFactors = FALSE
    )
  }))
}

r3ts10_summarize <- function(stage1_alpha_by_seed, stage2_by_seed) {
  contract <- r3ts10_contract()
  stage1_keys <- paste(
    stage1_alpha_by_seed$truth_mechanism,
    stage1_alpha_by_seed$seed,
    sprintf("%.3f", stage1_alpha_by_seed$alpha),
    sep = ":"
  )
  expected_stage1 <- as.vector(outer(
    r3ts10_expected_keys(),
    sprintf("%.3f", contract$alpha_grid),
    paste,
    sep = ":"
  ))
  r3ts10_assert(
    !anyDuplicated(stage1_keys) && setequal(stage1_keys, expected_stage1),
    "Stage-1 curves require all 20 replicates and all 40 alpha values."
  )
  stage2_keys <- paste(
    stage2_by_seed$truth_mechanism,
    stage2_by_seed$seed,
    stage2_by_seed$target,
    sep = ":"
  )
  expected_stage2 <- as.vector(outer(
    r3ts10_expected_keys(),
    contract$targets,
    paste,
    sep = ":"
  ))
  r3ts10_assert(
    !anyDuplicated(stage2_keys) && setequal(stage2_keys, expected_stage2),
    "Stage-2 summaries require all targets in all 20 replicates."
  )
  summarize_values <- function(values) {
    r3ts10_assert(all(is.finite(values)), "A Monte Carlo summary contains a non-finite value.")
    c(
      mean = mean(values),
      min = min(values),
      max = max(values),
      sd = stats::sd(values),
      mc_se = stats::sd(values) / sqrt(length(values))
    )
  }
  stage1 <- do.call(rbind, lapply(contract$mechanisms, function(mechanism) {
    do.call(rbind, lapply(contract$alpha_grid, function(alpha) {
      x <- stage1_alpha_by_seed[
        stage1_alpha_by_seed$truth_mechanism == mechanism &
          abs(stage1_alpha_by_seed$alpha - alpha) < 1e-12,
        , drop = FALSE
      ]
      power <- summarize_values(x$power)
      fdr <- summarize_values(x$empirical_fdr)
      calls <- summarize_values(x$calls)
      r3ts10_assert(nrow(x) == length(contract$seeds),
                    "Each Stage-1 point must contain ten seeds.")
      data.frame(
        truth_mechanism = mechanism,
        alpha = alpha,
        n_seeds = nrow(x),
        mean_power = power[["mean"]],
        power_min = power[["min"]],
        power_max = power[["max"]],
        power_sd = power[["sd"]],
        power_mc_se = power[["mc_se"]],
        mean_empirical_fdr = fdr[["mean"]],
        empirical_fdr_min = fdr[["min"]],
        empirical_fdr_max = fdr[["max"]],
        empirical_fdr_sd = fdr[["sd"]],
        empirical_fdr_mc_se = fdr[["mc_se"]],
        mean_calls = calls[["mean"]],
        min_calls = calls[["min"]],
        max_calls = calls[["max"]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  stage2 <- do.call(rbind, lapply(contract$mechanisms, function(mechanism) {
    do.call(rbind, lapply(contract$targets, function(target) {
      x <- stage2_by_seed[
        stage2_by_seed$truth_mechanism == mechanism & stage2_by_seed$target == target,
        , drop = FALSE
      ]
      conditional <- x$conditional_classification_power
      r3ts10_assert(
        all(is.finite(x$empirical_fsr)) && all(is.finite(x$empirical_power)) &&
          all(is.finite(conditional)),
        "A Stage-2 rate is undefined."
      )
      r3ts10_assert(nrow(x) == length(contract$seeds),
                    "Each Stage-2 target must contain ten seeds.")
      fsr <- summarize_values(x$empirical_fsr)
      cp <- summarize_values(conditional)
      calls <- summarize_values(x$calls)
      unconditional <- summarize_values(x$empirical_power)
      data.frame(
        truth_mechanism = mechanism,
        target = target,
        n_seeds = nrow(x),
        mean_empirical_fsr = fsr[["mean"]],
        empirical_fsr_min = fsr[["min"]],
        empirical_fsr_max = fsr[["max"]],
        empirical_fsr_sd = fsr[["sd"]],
        empirical_fsr_mc_se = fsr[["mc_se"]],
        mean_conditional_classification_power = cp[["mean"]],
        conditional_power_min = cp[["min"]],
        conditional_power_max = cp[["max"]],
        conditional_power_sd = cp[["sd"]],
        conditional_power_mc_se = cp[["mc_se"]],
        mean_calls = calls[["mean"]],
        min_calls = calls[["min"]],
        max_calls = calls[["max"]],
        calls_sd = calls[["sd"]],
        calls_mc_se = calls[["mc_se"]],
        mean_unconditional_power_audit = unconditional[["mean"]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  list(stage1 = stage1, stage2 = stage2)
}

r3ts10_stage2_display_values <- function(stage2_summary, mechanism) {
  contract <- r3ts10_contract()
  x <- stage2_summary[stage2_summary$truth_mechanism == mechanism, , drop = FALSE]
  x <- x[match(contract$targets, x$target), , drop = FALSE]
  r3ts10_assert(
    nrow(x) == length(contract$targets) && !anyNA(x$target),
    "The Stage-2 display table is missing a target."
  )
  out <- data.frame(
    Metric = c(
      "Empirical FSR",
      "Conditional classification power",
      "Mean number of calls"
    ),
    stringsAsFactors = FALSE
  )
  for (target in contract$targets) {
    row <- x[x$target == target, , drop = FALSE]
    out[[tools::toTitleCase(target)]] <- c(
      row$mean_empirical_fsr,
      row$mean_conditional_classification_power,
      row$mean_calls
    )
  }
  r3ts10_assert(
    !any(grepl("unconditional|end-to-end", out$Metric, ignore.case = TRUE)),
    "Unconditional power must not enter the formal Stage-2 table."
  )
  out
}

r3ts10_stage1_at_alpha <- function(stage1_summary, alpha = 0.05) {
  contract <- r3ts10_contract()
  r3ts10_assert(
    length(alpha) == 1L && is.finite(alpha) &&
      any(abs(contract$alpha_grid - alpha) < 1e-12),
    "The requested Stage-1 alpha is outside the saved grid."
  )
  out <- stage1_summary[abs(stage1_summary$alpha - alpha) < 1e-12, , drop = FALSE]
  out <- out[match(contract$mechanisms, out$truth_mechanism), , drop = FALSE]
  r3ts10_assert(
    nrow(out) == length(contract$mechanisms) && !anyNA(out$truth_mechanism),
    "The Stage-1 summary is incomplete at the requested alpha."
  )
  out
}
