# Validate completed paired comparisons before rendering the consolidated R3 page.

r3c_load_report <- function(root) {
  contract <- r3c_contract()
  base <- r3c_base(root)
  directory <- file.path(root, "output/revision_simulations/mc", contract$result_id)
  manifest <- r3c_manifest(directory)
  configuration <- readRDS(file.path(directory, "configuration.rds"))
  r3c_assert(identical(manifest$schema, contract$schema) &&
               identical(manifest$result_id, contract$result_id) &&
               identical(manifest$configuration, configuration) &&
               identical(configuration$contract, contract) &&
               identical(configuration$baseline_manifest_sha256, contract$base_manifest_sha256) &&
               identical(configuration$genotype_sha256, base$config$genotype_cache_sha256) &&
               identical(configuration$scientific_source_sha256[["comparator"]], contract$comparator_sha256) &&
               identical(configuration$package_version, contract$package_version) &&
               identical(configuration$package_sha, contract$package_sha),
             "Comparison cache configuration or provenance changed.")
  producer <- file.path(root, "code/revision_simulations/internal/r3_stage1_method_comparison")
  r3c_assert(identical(r3c_sha(file.path(producer, "comparison_helpers.R")),
                       configuration$experiment_source_sha256[["helpers"]]),
             "Comparison validation helpers differ from the production source.")
  complete <- readLines(file.path(directory, "complete.flag"), warn = FALSE)
  r3c_assert(all(c("replicates=20", "seeds=10", "retained_methods=6", "display_methods=4",
                    "stage1_alpha_rows=4800", "raw_and_bf_retained=IWP1,FASH-linear",
                    "stage2=unchanged", "donor_permutations=100") %in% complete),
             "The completion marker does not satisfy the paired comparison contract.")
  records <- lapply(base$records, function(reference) {
    name <- paste0(reference$truth_mechanism, "_seed_", reference$seed, ".rds")
    relative <- file.path("replicates", name)
    r3c_assert(relative %in% names(manifest$artifact_sha256), "The manifest omits a replicate.")
    record <- readRDS(file.path(directory, relative))
    r3c_assert(identical(record$configuration_token, manifest$configuration_token),
               "A comparison replicate uses a different configuration token.")
    r3c_validate_record(record, reference)
    record
  })
  curves <- do.call(rbind, lapply(records, function(x) r3c_curves(
    x$scores, x$true_dynamic, x$seed, x$truth_mechanism)))
  rownames(curves) <- NULL
  summary <- r3c_summarize(curves)
  tables <- readRDS(file.path(directory, "summary_tables.rds"))
  r3c_assert(isTRUE(all.equal(curves, tables$stage1_alpha_by_seed)) &&
               isTRUE(all.equal(summary, tables$stage1_mc_summary)) &&
               isTRUE(all.equal(summary[abs(summary$alpha - 0.05) < 1e-12, , drop = FALSE],
                                tables$stage1_alpha005_summary)),
             "Saved comparison curves do not reproduce from the all-unit scores.")
  for (name in c("functional_by_seed", "functional_mc_summary", "stage2_random_bspline", "stage2_raised_cosine")) {
    r3c_assert(identical(tables[[name]], base$tables[[name]]),
               paste("Stage 2 changed:", name))
  }
  iwp <- curves[curves$method == "FASH-IWP1-BF", names(base$tables$stage1_alpha_by_seed), drop = FALSE]
  rownames(iwp) <- NULL
  r3c_assert(isTRUE(all.equal(iwp, base$tables$stage1_alpha_by_seed)),
             "IWP1 BF curves differ from the original R3 page.")
  list(directory = directory, contract = contract, manifest = manifest, tables = tables,
       baseline = base, records = records)
}

r3c_load_truth_examples <- function(root, report) {
  base <- report$baseline
  directory <- file.path(root, "output/revision_simulations/mc", base$config$contract$formal_id)
  manifest <- r3c_manifest(directory, base$config$formal_manifest_sha256)
  r3c_assert("example_curves.rds" %in% names(manifest$artifact_sha256),
             "The pinned truth cache does not include example curves.")
  examples <- readRDS(file.path(directory, "example_curves.rds"))
  source_path <- file.path(root, "code/revision_simulations/r3_r4_fashr0143/source_snapshots",
                           "r3_full_universe_functional_simulation_functions.R")
  r3c_assert(identical(r3c_sha(source_path), base$config$contract$base_source_sha[["simulation"]]),
             "The frozen functional definitions changed.")
  definitions <- new.env(parent = globalenv())
  sys.source(source_path, envir = definitions)
  config <- base$config$formal_configuration
  result <- list()
  observed <- list()
  for (mechanism in report$contract$mechanisms) {
    reference <- base$records[[paste(mechanism, 12345L, sep = ":")]]
    original <- readRDS(file.path(directory, "replicates", paste0(mechanism, "_seed_12345.rds")))
    r3c_assert(identical(examples[[mechanism]], original$example_curves) &&
                 identical(original$selected_pair_keys, reference$selected_pair_keys),
               "Truth examples do not match the original seed and current pair order.")
    pool <- unlist(examples[[mechanism]], recursive = FALSE, use.names = FALSE)
    ids <- vapply(pool, function(x) x$variant_id, character(1))
    r3c_assert(length(pool) == 12L && !anyDuplicated(ids), "Incomplete or duplicated truth examples.")
    if (mechanism == "raised_cosine") {
      recovered_directory <- file.path(root,
        "output/revision_simulations/internal/r3_six_examples_20260918")
      recovered_path <- file.path(recovered_directory, "two_peak_example.rds")
      r3c_assert(identical(r3c_sha(recovered_path),
        readLines(file.path(recovered_directory, "two_peak_example.sha256"))),
        "Recovered two-peak illustration checksum changed.")
      recovered <- readRDS(recovered_path)
      expected_sources <- c(base$config$contract$base_source_sha[
        c("simulation", "genotype_helper", "digest_helper", "base_reconstruction")],
        comparator = report$contract$comparator_sha256)
      r3c_assert(identical(recovered$schema, "r3-two-peak-illustration-v1") &&
        identical(recovered$seed, 12345L) &&
        identical(recovered$baseline_manifest_sha256, report$contract$base_manifest_sha256) &&
        identical(recovered$scientific_source_sha256, expected_sources) &&
        identical(recovered$genotype_sha256, base$config$genotype_cache_sha256) &&
        identical(recovered$original_input_digests, reference$input_digests) &&
        identical(recovered$producer_sha256, r3c_sha(file.path(root,
          "code/revision_simulations/internal/r3_stage1_method_comparison/prepare_two_peak_example.R"))) &&
        is.finite(recovered$max_functional_error) && recovered$max_functional_error <= 1e-12 &&
        identical(dim(recovered$example_validation_errors), c(3L, 12L)) &&
        all(is.finite(recovered$example_validation_errors)) &&
        max(recovered$example_validation_errors) <= 1e-12 &&
        identical(recovered$example$spike_count, 2L),
        "Recovered two-peak illustration failed its provenance or numerical checks.")
      pool <- c(pool, list(recovered$example))
      r3c_assert(!anyDuplicated(vapply(pool, function(x) x$variant_id, character(1))),
        "The recovered illustration duplicates a stored example.")
    }
    for (example in pool) {
      index <- match(example$variant_id, reference$selected_pair_keys)
      curve <- example$true_curve
      measurements <- example$observed
      r3c_assert(!is.na(index) && isTRUE(reference$inference$true_dynamic[index]) &&
                   identical(curve$time, config$evaluation_grid) &&
                   all(is.finite(curve$true_effect)), "Invalid example unit or dense curve.")
      r3c_assert(nrow(measurements) == 16L &&
                   identical(measurements$time, as.numeric(0:15)) &&
                   all(is.finite(measurements$estimate)) &&
                   all(is.finite(measurements$se) & measurements$se > 0) &&
                   (mechanism != "raised_cosine" || example$spike_count %in% 1:3),
                 "Invalid stored effect estimates, standard errors, or generating peak count.")
      observed[[length(observed) + 1L]] <- data.frame(
        truth_mechanism = mechanism, seed = 12345L, unit_index = index,
        pair_key = example$variant_id, measurements
      )
      truth <- definitions$evaluate_temporal_functionals(
        matrix(curve$true_effect, nrow = 1L), smooth_var = curve$time,
        switch_threshold = config$switch_threshold, middle_window = config$middle_window,
        middle_boundary = config$middle_boundary
      )[1L, ]
      r3c_assert(isTRUE(all.equal(truth, example$true_functionals, tolerance = 1e-12)) &&
                   isTRUE(all.equal(unname(truth),
                     unname(reference$inference$true_functionals[index, ]), tolerance = 1e-12)) &&
                   isTRUE(all.equal(curve$true_effect[match(example$observed$time, curve$time)],
                     example$observed$true_effect, tolerance = 1e-12)),
                 "Dense example truth differs from current R3 functional truth or observed-time truth.")
      for (target in names(truth)) {
        result[[length(result) + 1L]] <- data.frame(
          truth_mechanism = mechanism, seed = 12345L, unit_index = index,
          pair_key = example$variant_id, target = target, true_functional = truth[[target]],
          peak_count = example$spike_count,
          time = curve$time, true_effect = curve$true_effect
        )
      }
    }
  }
  result <- do.call(rbind, result)
  attr(result, "observed") <- do.call(rbind, observed)
  result
}

r3c_peak_power <- function(root, report) {
  directory <- file.path(root, "output/revision_simulations/internal/r3_peak_stratified_power_20260918")
  path <- file.path(directory, "peak_labels.rds")
  r3c_assert(identical(r3c_sha(path), readLines(file.path(directory, "peak_labels.sha256"))),
             "Peak-label cache checksum changed.")
  cache <- readRDS(path)
  base <- report$baseline
  contract <- report$contract
  producer <- file.path(root, "code/revision_simulations/internal/r3_stage1_method_comparison",
                         "prepare_peak_stratification.R")
  expected_sources <- c(base$config$contract$base_source_sha[
    c("simulation", "genotype_helper", "digest_helper", "base_reconstruction")],
    comparator = contract$comparator_sha256)
  r3c_assert(identical(cache$schema, "r3-generating-peak-labels-v1") &&
    identical(cache$baseline_manifest_sha256, contract$base_manifest_sha256) &&
    identical(cache$genotype_sha256, base$config$genotype_cache_sha256) &&
    identical(cache$scientific_source_sha256, expected_sources) &&
    identical(cache$producer_sha256, r3c_sha(producer)) &&
    identical(names(cache$records), as.character(contract$seeds)),
    "Peak-label provenance or seed coverage differs from the paired experiment.")
  curves <- list()
  counts <- list()
  for (record in report$records) {
    if (record$truth_mechanism != "raised_cosine") next
    metadata <- cache$records[[as.character(record$seed)]]
    units <- metadata$units
    reference_truth <- base$records[[paste("raised_cosine", record$seed, sep = ":")]]$inference$true_functionals
    r3c_assert(identical(metadata$seed, record$seed) &&
      identical(metadata$true_beta_md5, unname(record$input_digests[["true_beta_md5"]])) &&
      identical(units$unit_index, seq_len(contract$n_units)) &&
      identical(units$pair_key, record$selected_pair_keys) &&
      identical(units$true_dynamic, record$true_dynamic) &&
      identical(dimnames(metadata$true_functionals), dimnames(reference_truth)) &&
      max(abs(metadata$true_functionals - reference_truth)) <= 1e-12 &&
      is.finite(metadata$max_functional_error) && metadata$max_functional_error <= 1e-12 &&
      all(units$peak_count[units$true_dynamic] %in% 1:3) &&
      all(is.na(units$peak_count[!units$true_dynamic])),
      "Peak labels are not aligned to this replicate's original truth and scores.")
    for (peaks in 1:3) {
      indices <- which(units$true_dynamic & units$peak_count == peaks)
      r3c_assert(length(indices) > 0L, "A peak stratum has no true dynamic units.")
      counts[[length(counts) + 1L]] <- data.frame(seed = record$seed,
        peak_count = peaks, n_true = length(indices))
      for (method in contract$methods) {
        # Keep the original global adjusted scores: no within-stratum reranking.
        scores <- record$scores[[method]]$adjusted_score[indices]
        detected <- vapply(contract$alpha_grid, function(alpha) sum(scores <= alpha), integer(1))
        curves[[length(curves) + 1L]] <- data.frame(seed = record$seed,
          peak_count = peaks, method = method, alpha = contract$alpha_grid,
          n_true = length(indices), n_detected = detected, power = detected / length(indices))
      }
    }
  }
  by_seed <- do.call(rbind, curves)
  counts <- do.call(rbind, counts)
  r3c_assert(nrow(by_seed) == 7200L && nrow(counts) == 30L,
             "Peak-stratified reporting requires all seeds, methods, and alpha values.")
  r3c_assert(all(counts$n_true == c(426L, 424L, 422L)[counts$peak_count]),
             "Peak-stratum sizes differ from the values stated on the page.")
  combined <- stats::aggregate(cbind(n_true, n_detected) ~ seed + method + alpha,
                               data = by_seed, FUN = sum)
  original <- report$tables$stage1_alpha_by_seed
  original <- original[original$truth_mechanism == "raised_cosine", ]
  key <- function(x) paste(x$seed, x$method, sprintf("%.3f", x$alpha), sep = ":")
  matched <- match(key(combined), key(original))
  r3c_assert(!anyNA(matched) && all(combined$n_true == 1272L) &&
    isTRUE(all.equal(combined$n_detected / combined$n_true, original$power[matched], tolerance = 1e-12)),
    "Count-weighted subgroup powers do not reproduce the original compact-truth curves.")
  summary <- do.call(rbind, lapply(split(by_seed, interaction(by_seed$peak_count,
    by_seed$method, by_seed$alpha, drop = TRUE)), function(x) {
      r3c_assert(nrow(x) == 10L && setequal(x$seed, contract$seeds), "Incomplete subgroup seed coverage.")
      data.frame(peak_count = x$peak_count[1L], method = x$method[1L], alpha = x$alpha[1L],
        mean_power = mean(x$power), power_min = min(x$power), power_max = max(x$power))
    }))
  rownames(summary) <- NULL
  list(by_seed = by_seed, summary = summary, counts = counts)
}

r3c_location_power <- function(report) {
  targets <- c("early", "middle", "late")
  expected_counts <- c(early = 369L, middle = 534L, late = 369L)
  contract <- report$contract
  curves <- list()
  counts <- list()
  for (record in report$records) {
    reference <- report$baseline$records[[paste(record$truth_mechanism, record$seed, sep = ":")]]
    truth <- reference$inference$true_functionals[, targets, drop = FALSE]
    membership <- truth > 0
    r3c_assert(identical(record$selected_pair_keys, reference$selected_pair_keys) &&
      identical(record$true_dynamic, reference$inference$true_dynamic) &&
      !anyNA(membership) &&
      all(rowSums(membership)[record$true_dynamic] == 1L) &&
      all(rowSums(membership)[!record$true_dynamic] == 0L),
      "Temporal truth does not partition the original dynamic units.")
    for (target in targets) {
      indices <- which(record$true_dynamic & membership[, target])
      r3c_assert(length(indices) == expected_counts[[target]],
                 "Temporal stratum size differs from the value stated on the page.")
      counts[[length(counts) + 1L]] <- data.frame(truth_mechanism = record$truth_mechanism, seed = record$seed,
        target = target, n_true = length(indices))
      for (method in contract$methods) {
        # Original global discovery rules; both Switch statuses remain included.
        scores <- record$scores[[method]]$adjusted_score[indices]
        detected <- vapply(contract$alpha_grid, function(alpha) sum(scores <= alpha), integer(1))
        curves[[length(curves) + 1L]] <- data.frame(truth_mechanism = record$truth_mechanism, seed = record$seed,
          target = target, method = method, alpha = contract$alpha_grid,
          n_true = length(indices), n_detected = detected, power = detected / length(indices))
      }
    }
  }
  by_seed <- do.call(rbind, curves)
  counts <- do.call(rbind, counts)
  r3c_assert(nrow(by_seed) == 14400L && nrow(counts) == 60L,
             "Temporal-stratified reporting lacks complete seed/method/alpha coverage.")
  combined <- stats::aggregate(cbind(n_true, n_detected) ~ truth_mechanism + seed + method + alpha,
                               data = by_seed, FUN = sum)
  original <- report$tables$stage1_alpha_by_seed
  key <- function(x) paste(x$truth_mechanism, x$seed, x$method, sprintf("%.3f", x$alpha), sep = ":")
  matched <- match(key(combined), key(original))
  r3c_assert(!anyNA(matched) && all(combined$n_true == 1272L) &&
    isTRUE(all.equal(combined$n_detected / combined$n_true, original$power[matched], tolerance = 1e-12)),
    "Count-weighted temporal powers do not reproduce the original curves for each truth mechanism.")
  summary <- do.call(rbind, lapply(split(by_seed, interaction(by_seed$truth_mechanism, by_seed$target,
    by_seed$method, by_seed$alpha, drop = TRUE)), function(x) {
      r3c_assert(nrow(x) == 10L && setequal(x$seed, contract$seeds), "Incomplete temporal-stratum seed coverage.")
      data.frame(truth_mechanism = x$truth_mechanism[1L], target = x$target[1L], method = x$method[1L], alpha = x$alpha[1L],
        mean_power = mean(x$power), power_min = min(x$power), power_max = max(x$power))
    }))
  rownames(summary) <- NULL
  list(by_seed = by_seed, summary = summary, counts = counts)
}

r3c_plot_stage1 <- function(curves, metric = c("power", "fdr")) {
  metric <- match.arg(metric)
  contract <- r3c_contract()
  r3c_assert(setequal(unique(curves$method), contract$display_methods) && nrow(curves) == 320L,
             "A main figure must contain four methods under both truth mechanisms.")
  labels <- c("FASH-IWP1 (BF)", "FASH-linear (BF)", "Direct linear", "Direct quadratic")
  names(labels) <- contract$display_methods
  curves$method_label <- factor(labels[curves$method], levels = unname(labels))
  curves$truth_label <- factor(curves$truth_mechanism, levels = contract$mechanisms,
    labels = c("Broad random B-spline", "Compact raised cosine"))
  r3c_assert(!anyNA(curves$method_label) && !anyNA(curves$truth_label), "Unmapped plot labels.")
  colors <- setNames(c("#0072B2", "#D55E00", "#009E73", "#CC79A7"), unname(labels))
  types <- setNames(c("solid", "longdash", "dotted", "dotdash"), unname(labels))
  columns <- if (metric == "power") c("mean_power", "power_min", "power_max") else
    c("mean_empirical_fdr", "empirical_fdr_min", "empirical_fdr_max")
  curves$mean <- curves[[columns[[1L]]]]
  curves$lower <- curves[[columns[[2L]]]]
  curves$upper <- curves[[columns[[3L]]]]
  r3c_assert(all(is.finite(as.matrix(curves[, c("mean", "lower", "upper")]))),
             "The figure contains non-finite rates.")
  upper <- if (metric == "power") 1 else min(1, 1.08 * max(0.20, curves$upper))
  plot <- ggplot2::ggplot(curves, ggplot2::aes(x = alpha, y = mean,
      color = method_label, fill = method_label, linetype = method_label)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.09,
                          color = NA, show.legend = FALSE) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = 0.05, color = "grey50", linetype = "dashed", linewidth = 0.4) +
    ggplot2::facet_wrap(~truth_label, nrow = 1) +
    ggplot2::scale_color_manual(values = colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = types, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = c(0.005, 0.05, 0.10, 0.15, 0.20)) +
    ggplot2::coord_cartesian(ylim = c(0, upper)) +
    ggplot2::labs(x = "Nominal Stage-1 FDR level", y = if (metric == "power")
      "Mean dynamic discovery power" else "Mean empirical dynamic FDR",
      color = NULL, fill = NULL, linetype = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank(),
                    strip.background = ggplot2::element_rect(fill = "grey95")) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE))
  if (metric == "fdr") plot <- plot +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "grey35", linetype = "dotted", linewidth = 0.4)
  plot
}
