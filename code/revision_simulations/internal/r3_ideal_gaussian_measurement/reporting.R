# Validate and format the current-versus-ideal R3 measurement comparison.

if (!exists("project_root", inherits = TRUE)) {
  project_root <- if (file.exists(
    "code/revision_simulations/shared/simulation_functions.R"
  )) {
    "."
  } else if (file.exists(
    "../code/revision_simulations/shared/simulation_functions.R"
  )) {
    ".."
  } else {
    stop("Could not locate the workflowr project root.", call. = FALSE)
  }
}

# The formal reporting layer is the authoritative validator for the current R3
# cache. Preserve its validated objects before loading the ideal experiment.
source(file.path(
  project_root,
  "code", "revision_simulations", "r3_functional_testing", "reporting.R"
))
formal_r3_configuration <- configuration
formal_r3_manifest <- r3_manifest
formal_r3_alpha_replicates <- mc_alpha_replicates
formal_r3_complete_flag <- r3_completion

ideal_result_id <- paste0(
  "r3_ideal_gaussian_known_t_adjusted_se_",
  paste0(
    "matched_truth_open_middle_3_12_center_aligned_",
    "iwp1_geometry_mixture_"
  ),
  "full_universe_",
  "fashr0143_pilot5"
)
ideal_schema_version <- "r3-ideal-gaussian-measurement-iwp1-mixture-v2"
ideal_output_dir <- file.path(
  project_root,
  "output", "revision_simulations", "internal", ideal_result_id
)
ideal_summary_dir <- file.path(ideal_output_dir, "summary")
ideal_manifest_path <- file.path(ideal_output_dir, "manifest.rds")
ideal_configuration_path <- file.path(ideal_output_dir, "configuration.rds")
ideal_complete_flag_path <- file.path(ideal_output_dir, "complete.flag")
ideal_example_path <- file.path(ideal_output_dir, "example_observations.rds")

ideal_summary_paths <- c(
  alpha_replicates = file.path(
    ideal_summary_dir, "all_replicate_functional_alpha_curves.csv"
  ),
  alpha005 = file.path(
    ideal_summary_dir, "all_replicate_functional_alpha005.csv"
  ),
  alpha_mc = file.path(
    ideal_summary_dir, "functional_testing_mc_alpha_curve.csv"
  ),
  pi0 = file.path(ideal_summary_dir, "all_replicate_pi0.csv"),
  truth_groups = file.path(
    ideal_summary_dir, "all_truth_group_counts.csv"
  ),
  truth_maf = file.path(ideal_summary_dir, "truth_maf_balance.csv"),
  moments = file.path(
    ideal_summary_dir, "standardized_error_moments.csv"
  ),
  histogram = file.path(
    ideal_summary_dir, "standardized_error_histogram.csv"
  ),
  lag_correlation = file.path(
    ideal_summary_dir, "standardized_error_lag_correlation.csv"
  ),
  digests = file.path(
    ideal_summary_dir, "replicate_input_digests.csv"
  ),
  middle_curve = file.path(ideal_summary_dir, "ideal_middle_curve.csv"),
  primary = file.path(ideal_summary_dir, "primary_middle_summary.csv")
)
ideal_replicate_paths <- unlist(lapply(
  c("random_bspline", "raised_cosine"),
  function(truth_mechanism) {
    file.path(
      ideal_output_dir,
      "replicates",
      paste0(
        truth_mechanism,
        "_seed_",
        c(12345L, 22345L, 32345L, 42345L, 52345L),
        ".rds"
      )
    )
  }
), use.names = FALSE)
ideal_required_paths <- c(
  ideal_manifest_path,
  ideal_configuration_path,
  ideal_complete_flag_path,
  ideal_example_path,
  unname(ideal_summary_paths),
  ideal_replicate_paths
)
if (any(!file.exists(ideal_required_paths))) {
  stop(
    "The ideal-Gaussian R3 cache is incomplete. Expected ten replicates, ",
    "all summaries, manifest.rds, and complete.flag under: ",
    ideal_output_dir,
    call. = FALSE
  )
}

report_sha256_file <- function(path) {
  command <- if (identical(Sys.info()[["sysname"]], "Darwin")) {
    "shasum"
  } else {
    "sha256sum"
  }
  arguments <- if (identical(command, "shasum")) c("-a", "256", path) else path
  output <- system2(command, arguments, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Unable to compute SHA-256 for ", path, ".", call. = FALSE)
  }
  sub("[[:space:]].*$", "", output[[1L]])
}

ideal_manifest <- readRDS(ideal_manifest_path)
ideal_configuration <- readRDS(ideal_configuration_path)
ideal_completion <- readLines(ideal_complete_flag_path, warn = FALSE)
ideal_alpha_replicates <- utils::read.csv(
  ideal_summary_paths[["alpha_replicates"]], stringsAsFactors = FALSE
)
ideal_alpha_005 <- utils::read.csv(
  ideal_summary_paths[["alpha005"]], stringsAsFactors = FALSE
)
ideal_alpha_mc <- utils::read.csv(
  ideal_summary_paths[["alpha_mc"]], stringsAsFactors = FALSE
)
ideal_pi0 <- utils::read.csv(
  ideal_summary_paths[["pi0"]], stringsAsFactors = FALSE
)
ideal_truth_groups <- utils::read.csv(
  ideal_summary_paths[["truth_groups"]], stringsAsFactors = FALSE
)
ideal_truth_maf <- utils::read.csv(
  ideal_summary_paths[["truth_maf"]], stringsAsFactors = FALSE
)
ideal_moments <- utils::read.csv(
  ideal_summary_paths[["moments"]], stringsAsFactors = FALSE
)
ideal_histogram <- utils::read.csv(
  ideal_summary_paths[["histogram"]], stringsAsFactors = FALSE
)
ideal_lag_correlation <- utils::read.csv(
  ideal_summary_paths[["lag_correlation"]], stringsAsFactors = FALSE
)
ideal_digests <- utils::read.csv(
  ideal_summary_paths[["digests"]], stringsAsFactors = FALSE
)
ideal_middle_curve_saved <- utils::read.csv(
  ideal_summary_paths[["middle_curve"]], stringsAsFactors = FALSE
)
ideal_primary_saved <- utils::read.csv(
  ideal_summary_paths[["primary"]], stringsAsFactors = FALSE
)
ideal_examples <- readRDS(ideal_example_path)

expected_seed_list <- c(12345L, 22345L, 32345L, 42345L, 52345L)
expected_mechanisms <- c("random_bspline", "raised_cosine")
expected_methods <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
expected_targets <- c("early", "middle", "late", "switch")
expected_alpha_grid <- seq(0.005, 0.20, by = 0.005)
expected_source_sha256 <- c(
  helper = "84bc06be91531f1587b0f298bea82f0bd937444663419c0b77d589b48d3fe84e",
  simulation_functions =
    "45267b0884168e5ae33cc4f14e3f05b711d961b65bf9c6fbd880e748de064a6e",
  real_genotype_helper =
    "c03c01a188503336a77793c96f4e2d3ac7e0cbd56f4028b552da4b2f88e6b9d7",
  genotype_cache =
    "81bbef5f323a0bab2ca993c782d8a9b7c63518b83c2cdb46ef7ed1d46f65af24",
  temporal_mixture_contract =
    "8e2b4c527b6c17fd2645b9d46e47a7b99410d2e9be1db1a38d541380d91ad723",
  runner =
    "3c59198b276ad841480e288c549ee50f684536b8c29e4c156608ff08dc191b6c"
)
expected_temporal_category_probs <- stats::setNames(
  c(0.29, 0.42, 0.29),
  c("early", "middle", "late")
)
expected_temporal_category_counts <- c(
  early = 369L,
  middle = 534L,
  late = 369L
)
expected_switch_status_counts <- c(
  switch = 637L,
  "non-switch" = 635L
)
expected_genotype_content_md5 <- c(
  `12345` = "526a7318aa2af901e09252f5a6ca3c46",
  `22345` = "517faa30d5218a956f1be84f2567369c",
  `32345` = "9b9be3205d7db54dac31763492bcb2eb",
  `42345` = "9ef73aa94a061df868b1a951fe495d9f",
  `52345` = "7dab159b8453e2f66188ae313bfbd611"
)
expected_truth_group_counts <- stats::setNames(
  c(185L, 184L, 267L, 267L, 185L, 184L),
  c(
    "early / switch",
    "early / non-switch",
    "middle / switch",
    "middle / non-switch",
    "late / switch",
    "late / non-switch"
  )
)

expected_replicate_keys <- expand.grid(
  seed = expected_seed_list,
  truth_mechanism = expected_mechanisms,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
expected_replicate_keys <- expected_replicate_keys[order(
  match(expected_replicate_keys$truth_mechanism, expected_mechanisms),
  match(expected_replicate_keys$seed, expected_seed_list)
), ]
rownames(expected_replicate_keys) <- NULL

if (!identical(ideal_manifest$schema_version, ideal_schema_version) ||
    !identical(ideal_manifest$result_id, ideal_result_id) ||
    !identical(ideal_configuration$schema_version, ideal_schema_version) ||
    !identical(ideal_configuration$output_id, ideal_result_id) ||
    !identical(ideal_manifest$configuration, ideal_configuration) ||
    !identical(
      ideal_configuration$base_r3_result_id,
      formal_r3_configuration$output_id
    ) ||
    !identical(ideal_configuration$J, 6362L) ||
    !identical(ideal_configuration$n_donors, 19L) ||
    !identical(ideal_configuration$n_covariates, 5L) ||
    !identical(ideal_configuration$time_grid, 0:15) ||
    !isTRUE(all.equal(
      ideal_configuration$evaluation_grid,
      seq(0, 15, by = 0.1)
    )) ||
    !identical(ideal_configuration$middle_window, c(3, 12)) ||
    !identical(ideal_configuration$middle_boundary, "open") ||
    !identical(ideal_configuration$middle_expression, "3 < t < 12") ||
    !identical(
      ideal_configuration$temporal_category_design,
      "canonical IWP1 temporal-category mixture"
    ) ||
    !isTRUE(all.equal(
      ideal_configuration$temporal_category_probs,
      expected_temporal_category_probs,
      tolerance = 1e-12,
      check.attributes = TRUE
    )) ||
    !identical(
      ideal_configuration$temporal_category_counts,
      expected_temporal_category_counts
    ) ||
    !identical(
      ideal_configuration$switch_status_counts,
      expected_switch_status_counts
    ) ||
    !identical(
      ideal_configuration$temporal_mixture_contract_sha256,
      expected_source_sha256[["temporal_mixture_contract"]]
    ) ||
    !identical(
      ideal_configuration$expected_truth_group_counts,
      expected_truth_group_counts
    ) ||
    !identical(
      ideal_configuration$observation_model,
      "ideal_gaussian_known_t_adjusted_se"
    ) ||
    !identical(
      ideal_configuration$standard_error_source,
      "R3 final t-adjusted SE regenerated from frozen seeds"
    ) ||
    !identical(ideal_configuration$ideal_measurement_seed_offset, 23003L) ||
    !identical(
      ideal_configuration$functional_candidate_scope,
      "full_universe"
    ) ||
    !identical(ideal_configuration$functional_candidate_universe_size, 6362L) ||
    !identical(
      ideal_configuration$functional_posterior_pairing,
      "common_random_seed_raw_bf"
    ) ||
    !identical(ideal_configuration$seed_list, expected_seed_list) ||
    !identical(ideal_configuration$truth_mechanisms, expected_mechanisms) ||
    !identical(
      ideal_configuration$genotype_content_digests,
      expected_genotype_content_md5
    ) ||
    !identical(
      ideal_manifest$source_provenance$sha256[
        names(expected_source_sha256)
      ],
      expected_source_sha256
    ) ||
    !identical(
      ideal_configuration$source_sha256[names(expected_source_sha256)],
      expected_source_sha256
    ) ||
    !identical(ideal_manifest$replicate_keys, expected_replicate_keys) ||
    !identical(ideal_manifest$package_provenance$version, "0.1.43") ||
    !identical(
      ideal_manifest$package_provenance$remote_sha,
      "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
    )) {
  stop("The ideal-Gaussian manifest or configuration is invalid.", call. = FALSE)
}

required_completion_lines <- c(
  paste0("result_id=", ideal_result_id),
  paste0("schema_version=", ideal_schema_version),
  paste0("base_r3_result_id=", formal_r3_configuration$output_id),
  "observation_model=ideal_gaussian_known_t_adjusted_se",
  "middle_definition=3 < t < 12",
  "temporal_category_probs=early:0.29;middle:0.42;late:0.29",
  paste0(
    "truth_group_counts=early / switch:185;",
    "early / non-switch:184;middle / switch:267;",
    "middle / non-switch:267;late / switch:185;",
    "late / non-switch:184"
  ),
  paste0(
    "temporal_mixture_contract_sha256=",
    expected_source_sha256[["temporal_mixture_contract"]]
  ),
  "summary_alpha_min=0.05",
  "summary_alpha_max=0.20",
  "summary_alpha_tolerance=1e-12",
  "summary_middle_rows=62",
  "functional_candidate_scope=full_universe",
  "functional_candidate_universe_size=6362",
  "functional_posterior_pairing=common_random_seed_raw_bf",
  "replicates=10"
)
if (!all(required_completion_lines %in% ideal_completion)) {
  stop("The ideal-Gaussian completion flag is invalid.", call. = FALSE)
}

relative_to_ideal <- function(path) {
  substring(
    normalizePath(path, winslash = "/", mustWork = TRUE),
    nchar(normalizePath(ideal_output_dir, winslash = "/", mustWork = TRUE)) + 2L
  )
}
artifact_paths_to_check <- c(
  ideal_configuration_path,
  ideal_example_path,
  unname(ideal_summary_paths),
  ideal_replicate_paths
)
artifact_relative_paths <- vapply(
  artifact_paths_to_check, relative_to_ideal, character(1)
)
expected_artifact_hashes <- ideal_manifest$artifact_sha256[
  artifact_relative_paths
]
observed_artifact_hashes <- stats::setNames(
  vapply(artifact_paths_to_check, report_sha256_file, character(1)),
  artifact_relative_paths
)
if (anyNA(expected_artifact_hashes) ||
    !identical(unname(observed_artifact_hashes), unname(expected_artifact_hashes))) {
  stop("One or more ideal-Gaussian cache artifacts failed SHA-256 validation.", call. = FALSE)
}

required_alpha_columns <- c(
  "scenario", "target", "method", "alpha", "candidate_scope",
  "candidate_count", "first_stage_null_calls", "n_discoveries",
  "false_discoveries", "power", "empirical_fsr", "seed",
  "truth_mechanism"
)
expected_alpha_rows <-
  length(expected_seed_list) *
  length(expected_mechanisms) *
  length(expected_methods) *
  length(expected_targets) *
  length(expected_alpha_grid)
if (!all(required_alpha_columns %in% names(ideal_alpha_replicates)) ||
    nrow(ideal_alpha_replicates) != expected_alpha_rows ||
    anyDuplicated(ideal_alpha_replicates[c(
      "seed", "truth_mechanism", "method", "target", "alpha"
    )]) ||
    !setequal(ideal_alpha_replicates$seed, expected_seed_list) ||
    !setequal(ideal_alpha_replicates$truth_mechanism, expected_mechanisms) ||
    !setequal(ideal_alpha_replicates$method, expected_methods) ||
    !setequal(ideal_alpha_replicates$target, expected_targets) ||
    !isTRUE(all.equal(
      sort(unique(ideal_alpha_replicates$alpha)),
      expected_alpha_grid
    )) ||
    any(ideal_alpha_replicates$candidate_scope != "full_universe") ||
    any(ideal_alpha_replicates$candidate_count != 6362L) ||
    any(ideal_alpha_replicates$first_stage_null_calls != 0L) ||
    any(!is.finite(ideal_alpha_replicates$power)) ||
    any(!is.finite(ideal_alpha_replicates$empirical_fsr)) ||
    nrow(ideal_alpha_005) != 80L ||
    nrow(ideal_alpha_mc) != 640L ||
    nrow(ideal_pi0) != 20L ||
    nrow(ideal_truth_groups) != 60L ||
    nrow(ideal_truth_maf) != 30L ||
    nrow(ideal_moments) != 20L ||
    nrow(ideal_histogram) != 2440L ||
    nrow(ideal_lag_correlation) != 300L ||
    nrow(ideal_digests) != 10L ||
    nrow(ideal_middle_curve_saved) != 62L ||
    nrow(ideal_primary_saved) != 2L ||
    nrow(ideal_examples) != 192L) {
  stop("The ideal-Gaussian summaries have an invalid shape or value contract.", call. = FALSE)
}

for (truth_mechanism in expected_mechanisms) {
  for (seed in expected_seed_list) {
    rows <- ideal_truth_groups[
      ideal_truth_groups$truth_mechanism == truth_mechanism &
        ideal_truth_groups$seed == seed,
      ,
      drop = FALSE
    ]
    observed_counts <- stats::setNames(
      as.integer(rows$n_dynamic), rows$truth_group
    )[names(expected_truth_group_counts)]
    if (!identical(observed_counts, expected_truth_group_counts)) {
      stop("An ideal-Gaussian truth-group count is invalid.", call. = FALSE)
    }
  }
}
histogram_totals <- stats::aggregate(
  count ~ seed + truth_mechanism + observation_model,
  data = ideal_histogram,
  FUN = sum
)
if (nrow(histogram_totals) != 20L ||
    any(histogram_totals$count != 6362L * 16L) ||
    !setequal(
      ideal_moments$observation_model,
      c("regression_t_adjusted", "ideal_gaussian")
    ) ||
    any(!is.finite(ideal_moments$mean)) ||
    any(!is.finite(ideal_moments$sd)) ||
    any(!is.finite(ideal_lag_correlation$mean_correlation))) {
  stop("The standardized-error diagnostics are invalid.", call. = FALSE)
}
md5_columns <- grep("md5$", names(ideal_digests), value = TRUE)
if (length(md5_columns) != 5L ||
    any(!vapply(
      ideal_digests[md5_columns],
      function(values) all(grepl("^[[:xdigit:]]{32}$", values)),
      logical(1)
    ))) {
  stop("The replicate input digests are invalid.", call. = FALSE)
}

strict_factor <- function(values, levels, name) {
  output <- factor(values, levels = levels)
  if (anyNA(output)) {
    stop(name, " contains a value outside its display levels.", call. = FALSE)
  }
  output
}

observation_model_levels <- c(
  "Regression + t-adjusted SE",
  "Ideal Gaussian, fixed adjusted SE"
)
mechanism_label_map <- c(
  random_bspline = "R3A: broad random B-spline",
  raised_cosine = "R3B: compact raised cosine"
)
target_label_map <- c(
  early = "Early", middle = "Middle", late = "Late", switch = "Switch"
)
method_label_map <- c(
  `FASH-IWP1-Raw` = "Raw",
  `FASH-IWP1-BF` = "BF-adjusted"
)

formal_alpha_for_comparison <- formal_r3_alpha_replicates
formal_alpha_for_comparison$observation_model <- observation_model_levels[[1L]]
ideal_alpha_for_comparison <- ideal_alpha_replicates
ideal_alpha_for_comparison$observation_model <- observation_model_levels[[2L]]
combined_alpha_replicates <- rbind(
  formal_alpha_for_comparison,
  ideal_alpha_for_comparison
)
combined_alpha_replicates$observation_model <- strict_factor(
  combined_alpha_replicates$observation_model,
  observation_model_levels,
  "observation_model"
)
combined_alpha_replicates$truth_mechanism_label <- strict_factor(
  unname(mechanism_label_map[combined_alpha_replicates$truth_mechanism]),
  unname(mechanism_label_map),
  "truth_mechanism"
)
combined_alpha_replicates$target_label <- strict_factor(
  unname(target_label_map[combined_alpha_replicates$target]),
  unname(target_label_map),
  "target"
)
combined_alpha_replicates$method_label <- strict_factor(
  unname(method_label_map[combined_alpha_replicates$method]),
  unname(method_label_map),
  "method"
)

summarize_replicate_curve <- function(data, metric) {
  if (!metric %in% c("power", "empirical_fsr")) {
    stop("metric must be power or empirical_fsr.", call. = FALSE)
  }
  key <- interaction(
    data$observation_model,
    data$truth_mechanism,
    data$target,
    data$method,
    data$alpha,
    drop = TRUE,
    lex.order = TRUE
  )
  output <- do.call(rbind, lapply(split(data, key), function(rows) {
    values <- rows[[metric]]
    data.frame(
      observation_model = as.character(rows$observation_model[[1L]]),
      truth_mechanism = rows$truth_mechanism[[1L]],
      truth_mechanism_label = as.character(rows$truth_mechanism_label[[1L]]),
      target = rows$target[[1L]],
      target_label = as.character(rows$target_label[[1L]]),
      method = rows$method[[1L]],
      method_label = as.character(rows$method_label[[1L]]),
      alpha = rows$alpha[[1L]],
      mean = mean(values),
      minimum = min(values),
      maximum = max(values),
      n_replications = nrow(rows),
      stringsAsFactors = FALSE
    )
  }))
  output$observation_model <- strict_factor(
    output$observation_model, observation_model_levels, "observation_model"
  )
  output$truth_mechanism_label <- strict_factor(
    output$truth_mechanism_label,
    unname(mechanism_label_map),
    "truth_mechanism_label"
  )
  output$target_label <- strict_factor(
    output$target_label, unname(target_label_map), "target_label"
  )
  output$method_label <- strict_factor(
    output$method_label, unname(method_label_map), "method_label"
  )
  output[order(
    output$truth_mechanism_label,
    output$target_label,
    output$method_label,
    output$observation_model,
    output$alpha
  ), ]
}

fsr_comparison_curve <- summarize_replicate_curve(
  combined_alpha_replicates[
    combined_alpha_replicates$method == "FASH-IWP1-BF",
    ,
    drop = FALSE
  ],
  "empirical_fsr"
)
power_comparison_curve <- summarize_replicate_curve(
  combined_alpha_replicates,
  "power"
)
if (nrow(fsr_comparison_curve) != 640L ||
    nrow(power_comparison_curve) != 1280L ||
    any(fsr_comparison_curve$n_replications != 5L) ||
    any(power_comparison_curve$n_replications != 5L)) {
  stop("The comparison curves have unexpected row counts.", call. = FALSE)
}

summarize_primary_middle <- function(data, observation_model) {
  rows <- data[
    data$method == "FASH-IWP1-BF" &
      data$target == "middle" &
      is.finite(data$alpha) &
      data$alpha >= 0.05 - 1e-12 &
      data$alpha <= 0.20 + 1e-12,
    ,
    drop = FALSE
  ]
  curve <- do.call(rbind, lapply(
    split(rows, list(rows$truth_mechanism, rows$alpha)),
    function(group) {
      data.frame(
        truth_mechanism = group$truth_mechanism[[1L]],
        alpha = group$alpha[[1L]],
        mean_empirical_fsr = mean(group$empirical_fsr),
        stringsAsFactors = FALSE
      )
    }
  ))
  do.call(rbind, lapply(split(curve, curve$truth_mechanism), function(group) {
    excess <- group$mean_empirical_fsr - group$alpha
    index <- which.max(excess)
    data.frame(
      truth_mechanism = group$truth_mechanism[[1L]],
      observation_model = observation_model,
      maximum_mean_fsr_excess = excess[[index]],
      alpha_at_maximum = group$alpha[[index]],
      mean_empirical_fsr_at_maximum = group$mean_empirical_fsr[[index]],
      stringsAsFactors = FALSE
    )
  }))
}

primary_current <- summarize_primary_middle(
  formal_alpha_for_comparison, observation_model_levels[[1L]]
)
primary_ideal <- summarize_primary_middle(
  ideal_alpha_for_comparison, observation_model_levels[[2L]]
)
primary_comparison <- merge(
  primary_current,
  primary_ideal,
  by = "truth_mechanism",
  suffixes = c("_current", "_ideal"),
  sort = FALSE
)
primary_comparison$current_minus_ideal_excess <-
  primary_comparison$maximum_mean_fsr_excess_current -
  primary_comparison$maximum_mean_fsr_excess_ideal
primary_comparison$truth_mechanism_label <- unname(
  mechanism_label_map[primary_comparison$truth_mechanism]
)
primary_comparison <- primary_comparison[match(
  expected_mechanisms, primary_comparison$truth_mechanism
), ]
rownames(primary_comparison) <- NULL
if (!isTRUE(all.equal(
  primary_ideal$maximum_mean_fsr_excess[match(
    ideal_primary_saved$truth_mechanism,
    primary_ideal$truth_mechanism
  )],
  ideal_primary_saved$maximum_mean_fsr_excess,
  tolerance = 1e-12,
  check.attributes = FALSE
))) {
  stop("The saved ideal primary summary does not reproduce from replicates.", call. = FALSE)
}

selected_alpha_levels <- c(0.05, 0.10, 0.15, 0.20)
middle_selected_rows <- combined_alpha_replicates[
  combined_alpha_replicates$method == "FASH-IWP1-BF" &
    combined_alpha_replicates$target == "middle" &
    combined_alpha_replicates$alpha %in% selected_alpha_levels,
  ,
  drop = FALSE
]
middle_selected_summary <- do.call(rbind, lapply(
  split(
    middle_selected_rows,
    interaction(
      middle_selected_rows$truth_mechanism,
      middle_selected_rows$observation_model,
      middle_selected_rows$alpha,
      drop = TRUE
    )
  ),
  function(rows) {
    data.frame(
      truth_mechanism = rows$truth_mechanism[[1L]],
      truth_mechanism_label = as.character(rows$truth_mechanism_label[[1L]]),
      observation_model = as.character(rows$observation_model[[1L]]),
      alpha = rows$alpha[[1L]],
      mean_empirical_fsr = mean(rows$empirical_fsr),
      minimum_seed_fsr = min(rows$empirical_fsr),
      maximum_seed_fsr = max(rows$empirical_fsr),
      mean_power = mean(rows$power),
      mean_calls = mean(rows$n_discoveries),
      mean_false_calls = mean(rows$false_discoveries),
      stringsAsFactors = FALSE
    )
  }
))
middle_selected_summary <- middle_selected_summary[order(
  match(middle_selected_summary$truth_mechanism, expected_mechanisms),
  match(middle_selected_summary$observation_model, observation_model_levels),
  middle_selected_summary$alpha
), ]
rownames(middle_selected_summary) <- NULL
if (nrow(middle_selected_summary) != 16L) {
  stop("The selected-alpha Middle table is incomplete.", call. = FALSE)
}

diagnostic_model_map <- c(
  regression_t_adjusted = observation_model_levels[[1L]],
  ideal_gaussian = observation_model_levels[[2L]]
)
ideal_moments$observation_model_label <- unname(
  diagnostic_model_map[ideal_moments$observation_model]
)
diagnostic_moment_columns <- c(
  "mean", "sd", "skewness", "excess_kurtosis", "tail_abs_1_96",
  "tail_abs_3"
)
standardized_error_summary <- do.call(rbind, lapply(
  split(
    ideal_moments,
    list(ideal_moments$truth_mechanism, ideal_moments$observation_model)
  ),
  function(rows) {
    values <- vapply(
      diagnostic_moment_columns,
      function(column) mean(rows[[column]]),
      numeric(1)
    )
    data.frame(
      truth_mechanism = rows$truth_mechanism[[1L]],
      truth_mechanism_label = unname(
        mechanism_label_map[rows$truth_mechanism[[1L]]]
      ),
      observation_model = rows$observation_model_label[[1L]],
      as.list(values),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
))
standardized_error_summary <- standardized_error_summary[order(
  match(standardized_error_summary$truth_mechanism, expected_mechanisms),
  match(standardized_error_summary$observation_model, observation_model_levels)
), ]
rownames(standardized_error_summary) <- NULL

finite_histogram <- ideal_histogram[
  is.finite(ideal_histogram$bin_left) &
    is.finite(ideal_histogram$bin_right),
  ,
  drop = FALSE
]
pooled_histogram <- stats::aggregate(
  count ~ truth_mechanism + observation_model + bin_left + bin_right,
  data = finite_histogram,
  FUN = sum
)
pooled_histogram$bin_midpoint <-
  (pooled_histogram$bin_left + pooled_histogram$bin_right) / 2
pooled_histogram$bin_width <-
  pooled_histogram$bin_right - pooled_histogram$bin_left
group_total <- stats::aggregate(
  count ~ truth_mechanism + observation_model,
  data = ideal_histogram,
  FUN = sum
)
names(group_total)[[3L]] <- "total_count"
pooled_histogram <- merge(
  pooled_histogram,
  group_total,
  by = c("truth_mechanism", "observation_model"),
  sort = FALSE
)
pooled_histogram$density <-
  pooled_histogram$count /
  pooled_histogram$total_count /
  pooled_histogram$bin_width
pooled_histogram$truth_mechanism_label <- strict_factor(
  unname(mechanism_label_map[pooled_histogram$truth_mechanism]),
  unname(mechanism_label_map),
  "histogram truth mechanism"
)
pooled_histogram$observation_model_label <- strict_factor(
  unname(diagnostic_model_map[pooled_histogram$observation_model]),
  observation_model_levels,
  "histogram observation model"
)

lag_correlation_summary <- stats::aggregate(
  mean_correlation ~ truth_mechanism + observation_model + lag,
  data = ideal_lag_correlation,
  FUN = mean
)
lag_correlation_summary$truth_mechanism_label <- strict_factor(
  unname(mechanism_label_map[lag_correlation_summary$truth_mechanism]),
  unname(mechanism_label_map),
  "lag truth mechanism"
)
lag_correlation_summary$observation_model_label <- strict_factor(
  unname(diagnostic_model_map[lag_correlation_summary$observation_model]),
  observation_model_levels,
  "lag observation model"
)

comparison_colors <- stats::setNames(
  c("#4D4D4D", "#0072B2"), observation_model_levels
)
comparison_linetypes <- stats::setNames(c("solid", "22"), observation_model_levels)

plot_fsr_comparison <- function(truth_mechanism) {
  data <- fsr_comparison_curve[
    fsr_comparison_curve$truth_mechanism == truth_mechanism,
    ,
    drop = FALSE
  ]
  if (nrow(data) != 320L) {
    stop("The FSR panel has an unexpected row count.", call. = FALSE)
  }
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = alpha,
      y = mean,
      color = observation_model,
      linetype = observation_model
    )
  ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      color = "grey55",
      linewidth = 0.45,
      linetype = "dotted"
    ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::facet_wrap(~target_label, ncol = 2) +
    ggplot2::scale_color_manual(values = comparison_colors, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = comparison_linetypes, drop = FALSE) +
    ggplot2::coord_cartesian(xlim = c(0, 0.20), ylim = c(0, NA)) +
    ggplot2::labs(
      x = "Nominal alpha",
      y = "Mean empirical FSR across five seeds",
      color = "Observation model",
      linetype = "Observation model"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}

plot_power_comparison <- function(truth_mechanism) {
  data <- power_comparison_curve[
    power_comparison_curve$truth_mechanism == truth_mechanism,
    ,
    drop = FALSE
  ]
  if (nrow(data) != 640L) {
    stop("The power panel has an unexpected row count.", call. = FALSE)
  }
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = alpha,
      y = mean,
      color = observation_model,
      linetype = observation_model
    )
  ) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::facet_grid(method_label ~ target_label) +
    ggplot2::scale_color_manual(values = comparison_colors, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = comparison_linetypes, drop = FALSE) +
    ggplot2::coord_cartesian(xlim = c(0, 0.20), ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Nominal alpha",
      y = "Mean power across five seeds",
      color = "Observation model",
      linetype = "Observation model"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")
}

plot_standardized_error_histogram <- function() {
  ggplot2::ggplot(
    pooled_histogram,
    ggplot2::aes(x = bin_midpoint, y = density)
  ) +
    ggplot2::geom_col(width = 0.1, fill = "#56B4E9", alpha = 0.75) +
    ggplot2::stat_function(
      fun = stats::dnorm,
      color = "#D55E00",
      linewidth = 0.7
    ) +
    ggplot2::facet_grid(
      truth_mechanism_label ~ observation_model_label
    ) +
    ggplot2::coord_cartesian(xlim = c(-4, 4)) +
    ggplot2::labs(
      x = "Standardized estimation error",
      y = "Pooled density"
    ) +
    ggplot2::theme_bw(base_size = 10)
}

plot_standardized_error_lag_correlation <- function() {
  ggplot2::ggplot(
    lag_correlation_summary,
    ggplot2::aes(
      x = lag,
      y = mean_correlation,
      color = observation_model_label,
      linetype = observation_model_label
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey65", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = 0.75) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::facet_wrap(~truth_mechanism_label, ncol = 1) +
    ggplot2::scale_color_manual(values = comparison_colors, drop = FALSE) +
    ggplot2::scale_linetype_manual(values = comparison_linetypes, drop = FALSE) +
    ggplot2::labs(
      x = "Time lag",
      y = "Mean correlation across time pairs and seeds",
      color = "Observation model",
      linetype = "Observation model"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "bottom")
}
