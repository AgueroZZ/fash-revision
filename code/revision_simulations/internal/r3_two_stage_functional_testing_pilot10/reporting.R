# Validate and display the completed ten-seed two-stage R3 cache.

r3ts10_show_table <- function(x, caption) {
  cat('<div class="r3-two-stage-table">\n\n')
  print(knitr::kable(x, caption = caption, row.names = FALSE, align = "l"))
  cat('\n\n</div>\n\n')
}

r3ts10_format_stage2_table <- function(stage2_summary, mechanism) {
  out <- r3ts10_stage2_display_values(stage2_summary, mechanism)
  for (column in names(out)[-1L]) {
    out[[column]] <- c(
      sprintf("%.3f", out[[column]][1:2]),
      sprintf("%.1f", out[[column]][3L])
    )
  }
  out
}

r3ts10_plot_stage1 <- function(stage1_summary,
                               metric = c("power", "fdr"),
                               alpha_reference = 0.05) {
  metric <- match.arg(metric)
  contract <- r3ts10_contract()
  columns <- if (metric == "power") {
    c(mean = "mean_power", lower = "power_min", upper = "power_max")
  } else {
    c(
      mean = "mean_empirical_fdr",
      lower = "empirical_fdr_min",
      upper = "empirical_fdr_max"
    )
  }
  required <- c("truth_mechanism", "alpha", unname(columns))
  r3ts10_assert(
    all(required %in% names(stage1_summary)) &&
      setequal(stage1_summary$truth_mechanism, contract$mechanisms),
    "The Stage-1 plotting summary is incomplete."
  )
  colors <- c(random_bspline = "#0072B2", raised_cosine = "#D55E00")
  line_types <- c(random_bspline = 1L, raised_cosine = 2L)
  x_limits <- range(contract$alpha_grid)
  if (metric == "power") {
    y_limits <- c(0, 1)
    y_label <- "Mean dynamic discovery power"
    title <- "Dynamic discovery power"
  } else {
    y_upper <- max(c(x_limits[[2L]], stage1_summary[[columns[["upper"]]]]))
    y_limits <- c(0, min(1, y_upper * 1.08))
    y_label <- "Mean empirical dynamic FDR"
    title <- "Empirical dynamic FDR"
  }
  old_parameters <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_parameters), add = TRUE)
  graphics::par(mar = c(4.4, 4.6, 3.4, 1.2))
  graphics::plot(
    NA,
    xlim = x_limits,
    ylim = y_limits,
    xlab = "Nominal Stage-1 FDR level alpha",
    ylab = y_label,
    main = title
  )
  graphics::grid(col = "gray90")
  if (metric == "fdr") {
    graphics::abline(a = 0, b = 1, col = "gray35", lty = 3, lwd = 1.5)
  }
  graphics::abline(v = alpha_reference, col = "gray55", lty = 3, lwd = 1.3)
  for (mechanism in contract$mechanisms) {
    curve <- stage1_summary[
      stage1_summary$truth_mechanism == mechanism,
      ,
      drop = FALSE
    ]
    curve <- curve[order(curve$alpha), , drop = FALSE]
    lower <- pmax(y_limits[[1L]], curve[[columns[["lower"]]]])
    upper <- pmin(y_limits[[2L]], curve[[columns[["upper"]]]])
    graphics::polygon(
      x = c(curve$alpha, rev(curve$alpha)),
      y = c(lower, rev(upper)),
      col = grDevices::adjustcolor(colors[[mechanism]], alpha.f = 0.14),
      border = NA
    )
    graphics::lines(
      curve$alpha,
      curve[[columns[["mean"]]]],
      col = colors[[mechanism]],
      lty = line_types[[mechanism]],
      lwd = 2.4
    )
  }
  legend_labels <- unname(contract$mechanism_labels[contract$mechanisms])
  legend_colors <- unname(colors[contract$mechanisms])
  legend_types <- unname(line_types[contract$mechanisms])
  if (metric == "fdr") {
    legend_labels <- c(legend_labels, "Nominal alpha")
    legend_colors <- c(legend_colors, "gray35")
    legend_types <- c(legend_types, 3L)
  }
  graphics::legend(
    "topleft",
    legend = legend_labels,
    col = legend_colors,
    lty = legend_types,
    lwd = c(rep(2.4, length(contract$mechanisms)), if (metric == "fdr") 1.5),
    bty = "n",
    cex = 0.88
  )
  invisible(stage1_summary)
}

r3ts10_load_report <- function(project_root) {
  contract <- r3ts10_contract()
  experiment <- file.path(
    project_root,
    "code/revision_simulations/internal/r3_two_stage_functional_testing_pilot10"
  )
  base_experiment <- file.path(
    project_root,
    "code/revision_simulations/internal/r3_two_stage_functional_testing"
  )
  directory <- file.path(
    project_root,
    "output/revision_simulations/mc",
    contract$result_id
  )
  required <- c(
    "complete.flag", "manifest.rds", "configuration.rds", "summary_tables.rds",
    "sessionInfo.txt", file.path("replicates", r3ts10_replicate_names())
  )
  r3ts10_assert(
    dir.exists(directory) && all(file.exists(file.path(directory, required))),
    "The completed ten-seed two-stage R3 cache is unavailable."
  )
  completion <- readLines(file.path(directory, "complete.flag"), warn = FALSE)
  required_completion <- c(
    paste0("result_id=", contract$result_id),
    "replicates=20",
    "seeds=10",
    "base_replicates_reused=10",
    "new_replicates_computed=10",
    "stage1_method=FASH-IWP1-BF",
    "stage1_alpha_grid=0.005:0.005:0.200",
    "q_dyn=0.05",
    "q_func=0.05",
    "fit_universe=6362",
    "functional_candidate_scope=dynamic_fdr_screen",
    "stage2_primary_power=conditional_classification_power",
    "unconditional_power=retained_for_audit_only"
  )
  r3ts10_assert(
    all(required_completion %in% completion),
    "The completion marker does not satisfy the formal two-stage contract."
  )
  manifest <- readRDS(file.path(directory, "manifest.rds"))
  configuration <- readRDS(file.path(directory, "configuration.rds"))
  r3ts10_assert(
    identical(manifest$schema, contract$schema) &&
      identical(manifest$result_id, contract$result_id) &&
      identical(manifest$configuration, configuration) &&
      is.character(manifest$configuration_token) &&
      length(manifest$configuration_token) == 1L &&
      nzchar(manifest$configuration_token) &&
      identical(configuration$contract, contract) &&
      identical(configuration$fit_universe, "all_units") &&
      identical(configuration$prior_refit_after_screening, FALSE) &&
      identical(configuration$functional_candidate_scope, "dynamic_fdr_screen") &&
      identical(configuration$stage1_method, "FASH-IWP1-BF") &&
      identical(configuration$stage1_alpha_grid, contract$alpha_grid) &&
      identical(configuration$stage2_q_dyn, contract$q_dyn) &&
      identical(configuration$stage2_q_func, contract$q_func) &&
      identical(configuration$package$version, contract$package_version) &&
      identical(configuration$package$remote_sha, contract$package_sha),
    "The result cache does not implement the fixed pilot-10 contract."
  )
  # The configuration file hash and record tokens validate producer identity.
  # Re-serializing configuration on another R version changes its header.
  r3ts10_assert(
    all(c("configuration.rds", "summary_tables.rds",
          file.path("replicates", r3ts10_replicate_names())) %in%
        names(manifest$artifact_sha256)),
    "The manifest omits required configuration or result hashes."
  )
  actual_artifact_sha256 <- vapply(
    file.path(directory, names(manifest$artifact_sha256)),
    r3ts10_sha256,
    character(1)
  )
  r3ts10_assert(
    identical(unname(actual_artifact_sha256), unname(manifest$artifact_sha256)),
    "A result artifact no longer matches the production manifest."
  )
  current_sources <- c(
    pilot10_helpers = file.path(experiment, "pilot10_helpers.R"),
    base_two_stage_helper = file.path(base_experiment, "two_stage_helpers.R")
  )
  current_hashes <- vapply(current_sources, r3ts10_sha256, character(1))
  r3ts10_assert(
    identical(
      current_hashes[["pilot10_helpers"]],
      configuration$experiment_source_sha256[["pilot10_helpers"]]
    ) &&
      identical(
        current_hashes[["base_two_stage_helper"]],
        contract$base_source_sha[["base_two_stage_helper"]]
      ),
    "Current reporting helpers differ from the production sources."
  )
  base_environment <- new.env(parent = environment())
  sys.source(current_sources[["base_two_stage_helper"]], envir = base_environment)
  records <- lapply(
    file.path(directory, "replicates", r3ts10_replicate_names()),
    readRDS
  )
  invisible(lapply(records, function(record) {
    r3ts10_assert(
      identical(record$configuration_token, manifest$configuration_token),
      "A replicate has a different configuration token."
    )
    r3ts10_validate_record(
      record,
      base_environment = base_environment,
      cfsr_function = r3ts10_functional_cfsr_table
    )
  }))
  record_keys <- vapply(
    records,
    function(record) paste(record$truth_mechanism, record$seed, sep = ":"),
    character(1)
  )
  r3ts10_assert(
    !anyDuplicated(record_keys) && setequal(record_keys, r3ts10_expected_keys()),
    "The cache lacks one or more mechanism-by-seed records."
  )
  tables <- readRDS(file.path(directory, "summary_tables.rds"))
  expected_tables <- c(
    "stage1_alpha_by_seed", "stage1_mc_summary", "stage1_alpha005_summary",
    "functional_by_seed", "functional_mc_summary", "genotype_validation",
    "stage2_random_bspline", "stage2_raised_cosine",
    "unconditional_power_audit", "input_validation", "all_stage1_units",
    "all_stage2_candidates"
  )
  r3ts10_assert(
    identical(names(tables), expected_tables),
    "The saved summary-table schema changed."
  )
  reconstructed_stage1 <- do.call(rbind, lapply(
    records,
    r3ts10_stage1_alpha_by_seed,
    base_environment = base_environment
  ))
  rownames(reconstructed_stage1) <- NULL
  reconstructed_functional <- do.call(rbind, lapply(
    records,
    function(record) record$inference$metrics$functional
  ))
  rownames(reconstructed_functional) <- NULL
  summaries <- r3ts10_summarize(reconstructed_stage1, reconstructed_functional)
  r3ts10_assert(
    isTRUE(all.equal(tables$stage1_alpha_by_seed, reconstructed_stage1)) &&
      isTRUE(all.equal(tables$functional_by_seed, reconstructed_functional)) &&
      isTRUE(all.equal(tables$stage1_mc_summary, summaries$stage1)) &&
      isTRUE(all.equal(tables$functional_mc_summary, summaries$stage2)) &&
      isTRUE(all.equal(
        tables$stage1_alpha005_summary,
        r3ts10_stage1_at_alpha(summaries$stage1, contract$q_dyn)
      )),
    "The saved Monte Carlo summaries do not reproduce from replicate records."
  )
  for (mechanism in contract$mechanisms) {
    table_name <- paste0("stage2_", mechanism)
    r3ts10_assert(
      # Summation may differ at machine precision across R platforms.
      isTRUE(all.equal(
        tables[[table_name]],
        r3ts10_stage2_display_values(summaries$stage2, mechanism),
        tolerance = 1e-12,
        scale = 1
      )),
      paste("The Stage-2 display table changed for", mechanism, ".")
    )
  }
  input_validation <- tables$input_validation
  r3ts10_assert(
    nrow(input_validation) == 20L &&
      sum(input_validation$record_kind == "validated_pilot5_reuse") == 10L &&
      sum(input_validation$record_kind == "new_seed_frozen_generation") == 10L &&
      all(input_validation$fitted_n_units == contract$n_units) &&
      all(is.finite(input_validation$fitted_bf_pi0)) &&
      all(input_validation$fitted_bf_pi0 >= 0) &&
      all(input_validation$fitted_bf_pi0 <= 1),
    "The replicate-level provenance or full-fit validation is incomplete."
  )
  r3ts10_assert(
    nrow(tables$genotype_validation) == 10L &&
      setequal(tables$genotype_validation$seed, contract$seeds) &&
      all(tables$genotype_validation$genes == contract$n_units) &&
      all(tables$genotype_validation$donors == 19L) &&
      all(tables$genotype_validation$maf_min >= 0.10),
    "The ten-seed genotype validation is incomplete."
  )
  genotype_path <- file.path(
    project_root,
    "output/revision_simulations/shared",
    contract$genotype_id,
    "genotype_samples.rds"
  )
  r3ts10_assert(
    file.exists(genotype_path) &&
      identical(
        r3ts10_sha256(genotype_path),
        configuration$genotype_cache_sha256
      ),
    "The local ten-seed genotype cache differs from the fitted cache."
  )
  list(
    directory = directory,
    contract = contract,
    manifest = manifest,
    configuration = configuration,
    tables = tables
  )
}
