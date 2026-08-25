# Validate and expose the targeted fashr 0.1.43 R1/R2 reporting cache.

R1_R2_FASHR0143_CACHE_DIR <- file.path(
  "output", "revision_simulations", "mc", "r1_r2_fashr0143"
)
R1_R2_FASHR0143_VERSION <- "0.1.43"
R1_R2_FASHR0143_REMOTE_SHA <-
  "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
R1_R2_FASHR0143_SEEDS <- c(12345L, 22345L, 32345L, 42345L, 52345L)
R1_R2_FASHR0143_SOURCE_SHA256 <- c(
  simulation_functions =
    "93f9a2c5606ae74763fb53e189af62c4d3b7c973aaf38d5b3eb3bf32f97a6487",
  real_genotype_helper =
    "c03c01a188503336a77793c96f4e2d3ac7e0cbd56f4028b552da4b2f88e6b9d7",
  genotype_cache =
    "81bbef5f323a0bab2ca993c782d8a9b7c63518b83c2cdb46ef7ed1d46f65af24",
  runner =
    "0301a2e9499c0a5ee4262ea8679e91926ca5e6f7425256050d15a4a027c9481d"
)

r1_r2_fashr0143_sha256_file <- function(path) {
  if (!file.exists(path) || dir.exists(path)) {
    stop("Required file does not exist: ", path, call. = FALSE)
  }
  command <- if (identical(Sys.info()[["sysname"]], "Darwin")) {
    "shasum"
  } else {
    "sha256sum"
  }
  arguments <- if (identical(command, "shasum")) {
    c("-a", "256", path)
  } else {
    path
  }
  output <- system2(command, arguments, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Unable to compute SHA-256 for ", path, ".", call. = FALSE)
  }
  sub("[[:space:]].*$", "", output[[1L]])
}

r1_r2_fashr0143_require_equal <- function(
    observed,
    expected,
    label,
    tolerance = 0,
    check_attributes = TRUE) {
  comparison <- all.equal(
    observed,
    expected,
    tolerance = tolerance,
    check.attributes = check_attributes
  )
  if (!isTRUE(comparison)) {
    stop(label, " failed: ", paste(comparison, collapse = "; "), call. = FALSE)
  }
  invisible(observed)
}

r1_r2_fashr0143_parse_flag <- function(path) {
  lines <- readLines(path, warn = FALSE)
  pieces <- strsplit(lines, "=", fixed = TRUE)
  if (length(lines) == 0L || any(lengths(pieces) != 2L)) {
    stop("The R1/R2 complete.flag has an invalid format.", call. = FALSE)
  }
  values <- vapply(pieces, `[[`, character(1), 2L)
  names(values) <- vapply(pieces, `[[`, character(1), 1L)
  values
}

r1_r2_fashr0143_read_csv <- function(cache_dir, file_name) {
  utils::read.csv(
    file.path(cache_dir, "summary", file_name),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

r1_r2_fashr0143_validate_manifest <- function(cache_dir) {
  expected_artifacts <- unlist(lapply(
    c("r1", "r2"),
    function(scenario_id) {
      file.path(
        "replicates",
        scenario_id,
        paste0("seed_", R1_R2_FASHR0143_SEEDS, ".rds")
      )
    }
  ))
  expected_summaries <- file.path("summary", c(
    "r1_all_replicate_fash_alpha_curves.csv",
    "r1_all_replicate_linear_prior_summary.csv",
    "r1_all_replicate_linear_prior_weights.csv",
    "r1_all_replicate_pi0.csv",
    "r1_fash_mc_alpha005_summary.csv",
    "r1_fash_mc_alpha_curve.csv",
    "r1_genotype_selection_summary.csv",
    "r1_mc_pi0_summary.csv",
    "r1_truth_maf_balance.csv",
    "r2_all_replicate_fash_alpha_curves.csv",
    "r2_all_replicate_geometry.csv",
    "r2_all_replicate_linear_prior_summary.csv",
    "r2_all_replicate_linear_prior_weights.csv",
    "r2_all_replicate_peak_alpha_curves.csv",
    "r2_all_replicate_pi0.csv",
    "r2_fash_mc_alpha005_summary.csv",
    "r2_fash_mc_alpha_curve.csv",
    "r2_genotype_selection_summary.csv",
    "r2_mc_pi0_summary.csv",
    "r2_peak_mc_alpha005_summary.csv",
    "r2_peak_mc_alpha_curve.csv",
    "r2_truth_maf_balance.csv"
  ))
  required_paths <- file.path(
    cache_dir,
    c("complete.flag", "manifest.rds", expected_artifacts, expected_summaries)
  )
  missing_paths <- required_paths[!file.exists(required_paths)]
  if (length(missing_paths) > 0L) {
    stop(
      "The targeted R1/R2 fashr 0.1.43 cache is incomplete. Missing: ",
      paste(missing_paths, collapse = ", "),
      call. = FALSE
    )
  }

  completion <- r1_r2_fashr0143_parse_flag(
    file.path(cache_dir, "complete.flag")
  )
  manifest <- readRDS(file.path(cache_dir, "manifest.rds"))
  stopifnot(
    identical(manifest$schema_version, "r1-r2-fashr0143-manifest-v1"),
    identical(manifest$result_id, "r1_r2_fashr0143"),
    identical(completion[["result_id"]], manifest$result_id),
    identical(completion[["fashr_version"]], R1_R2_FASHR0143_VERSION),
    identical(
      completion[["fashr_remote_sha"]],
      R1_R2_FASHR0143_REMOTE_SHA
    ),
    identical(completion[["r1_replicates"]], "5"),
    identical(completion[["r2_replicates"]], "5"),
    identical(
      completion[["direct_interaction_tests"]],
      "reused_from_historical_formal_caches"
    ),
    identical(
      manifest$package_provenance$version,
      R1_R2_FASHR0143_VERSION
    ),
    identical(
      manifest$package_provenance$remote_sha,
      R1_R2_FASHR0143_REMOTE_SHA
    ),
    identical(
      manifest$configuration$direct_interaction_tests,
      "excluded; reuse formal historical direct caches"
    ),
    identical(
      as.integer(manifest$configuration$seed_list),
      R1_R2_FASHR0143_SEEDS
    )
  )

  recorded_source_sha256 <- unlist(
    manifest$source_provenance$sha256,
    use.names = TRUE
  )
  r1_r2_fashr0143_require_equal(
    recorded_source_sha256[names(R1_R2_FASHR0143_SOURCE_SHA256)],
    R1_R2_FASHR0143_SOURCE_SHA256,
    "R1/R2 recorded source provenance"
  )

  expected_hashes <- c(
    unlist(manifest$artifact_sha256, use.names = TRUE),
    unlist(manifest$summary_sha256, use.names = TRUE)
  )
  if (!setequal(names(manifest$artifact_sha256), expected_artifacts) ||
      !setequal(names(manifest$summary_sha256), expected_summaries) ||
      anyDuplicated(names(expected_hashes))) {
    stop("The R1/R2 manifest file inventory is invalid.", call. = FALSE)
  }
  observed_hashes <- vapply(
    file.path(cache_dir, names(expected_hashes)),
    r1_r2_fashr0143_sha256_file,
    character(1)
  )
  names(observed_hashes) <- names(expected_hashes)
  r1_r2_fashr0143_require_equal(
    observed_hashes,
    expected_hashes,
    "R1/R2 retained output SHA-256 validation"
  )

  local_package <- utils::packageDescription("fashr")
  stopifnot(
    identical(as.character(utils::packageVersion("fashr")),
              R1_R2_FASHR0143_VERSION),
    identical(local_package$RemoteSha, R1_R2_FASHR0143_REMOTE_SHA)
  )

  local_sources <- c(
    simulation_functions = file.path(
      "code", "revision_simulations", "shared", "simulation_functions.R"
    ),
    real_genotype_helper = file.path(
      "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
    ),
    genotype_cache = file.path(
      "output", "revision_simulations", "shared",
      "real_genotype_one_per_gene_J6362_pilot5", "genotype_samples.rds"
    ),
    runner = file.path(
      "code", "revision_simulations", "r1_r2_fashr0143",
      "11_r1_r2_fashr0143.R"
    )
  )
  observed_source_sha256 <- vapply(
    local_sources,
    r1_r2_fashr0143_sha256_file,
    character(1)
  )
  r1_r2_fashr0143_require_equal(
    observed_source_sha256,
    R1_R2_FASHR0143_SOURCE_SHA256,
    "R1/R2 local source provenance"
  )
  manifest
}

r1_r2_fashr0143_require_numeric_pairing <- function(
    current,
    historical,
    label,
    tolerance = 1e-12) {
  if (!identical(dim(current), dim(historical)) ||
      length(current) != length(historical)) {
    stop(label, " has incompatible dimensions.", call. = FALSE)
  }
  r1_r2_fashr0143_require_equal(
    as.numeric(current),
    as.numeric(historical),
    label,
    tolerance = tolerance,
    check_attributes = FALSE
  )
}

r1_r2_fashr0143_validate_payload <- function(
    artifact,
    historical_out,
    scenario_id) {
  payload <- artifact$example_payload
  historical_pair_keys <- historical_out$genotype_selection$pair_key
  stopifnot(
    identical(artifact$schema_version, "r1-r2-fashr0143-seed-v2"),
    identical(artifact$result_id, "r1_r2_fashr0143"),
    identical(artifact$scenario_id, scenario_id),
    identical(artifact$seed, 12345L),
    identical(artifact$selected_pair_keys, historical_pair_keys),
    identical(payload$variant_info$unit_id, historical_pair_keys),
    identical(
      payload$unit_info$unit_id,
      historical_out$unit_info$unit_id
    ),
    identical(
      payload$unit_info$effect_class,
      historical_out$unit_info$effect_class
    ),
    length(artifact$warnings) == 0L,
    identical(
      artifact$package_provenance$version,
      R1_R2_FASHR0143_VERSION
    ),
    identical(
      artifact$package_provenance$remote_sha,
      R1_R2_FASHR0143_REMOTE_SHA
    )
  )

  if (scenario_id == "r2") {
    stopifnot(
      identical(
        payload$unit_info$cell_id,
        historical_out$unit_info$cell_id
      ),
      identical(
        payload$unit_info$spike_count,
        historical_out$unit_info$spike_count
      ),
      identical(
        payload$unit_info$time_group,
        historical_out$unit_info$time_group
      )
    )
  }

  r1_r2_fashr0143_require_numeric_pairing(
    payload$genotype,
    historical_out$genotype,
    paste0(toupper(scenario_id), " genotype pairing"),
    tolerance = 0
  )
  r1_r2_fashr0143_require_numeric_pairing(
    payload$covariates,
    historical_out$covariates,
    paste0(toupper(scenario_id), " covariate pairing")
  )
  r1_r2_fashr0143_require_numeric_pairing(
    payload$true_beta,
    historical_out$true_beta,
    paste0(toupper(scenario_id), " truth pairing")
  )
  r1_r2_fashr0143_require_numeric_pairing(
    payload$expression,
    historical_out$expression,
    paste0(toupper(scenario_id), " expression pairing")
  )
  for (field in c("beta_hat", "se", "se_uncorrected")) {
    r1_r2_fashr0143_require_numeric_pairing(
      payload$eqtl_summary[[field]],
      historical_out$eqtl_summary[[field]],
      paste0(toupper(scenario_id), " ", field, " pairing")
    )
  }
  r1_r2_fashr0143_require_numeric_pairing(
    payload$variant_info$observed_maf,
    pmin(
      historical_out$variant_info$observed_maf,
      1 - historical_out$variant_info$observed_maf
    ),
    paste0(toupper(scenario_id), " observed-MAF pairing")
  )
  if (scenario_id == "r2") {
    r1_r2_fashr0143_require_numeric_pairing(
      payload$evaluation_grid,
      historical_out$evaluation_grid,
      "R2 evaluation-grid pairing",
      tolerance = 0
    )
    r1_r2_fashr0143_require_numeric_pairing(
      payload$true_beta_evaluation,
      historical_out$true_beta_evaluation,
      "R2 evaluation-truth pairing"
    )
  }
  invisible(payload)
}

r1_r2_fashr0143_sort_curve <- function(table) {
  ordering_columns <- intersect(
    c("method", "spike_count", "shape_profile", "alpha"),
    names(table)
  )
  ordering <- do.call(order, unname(table[ordering_columns]))
  table[ordering, , drop = FALSE]
}

r1_r2_fashr0143_validate_curve <- function(
    table,
    methods,
    label,
    expected_rows_per_method = 41L) {
  if (!setequal(unique(table$method), methods) ||
      anyDuplicated(table[intersect(
        c("method", "spike_count", "shape_profile", "alpha"),
        names(table)
      )]) ||
      anyNA(table)) {
    stop(label, " has invalid method-alpha keys or missing values.", call. = FALSE)
  }
  if ("spike_count" %in% names(table)) {
    expected_rows <- length(methods) * 3L * expected_rows_per_method
  } else {
    expected_rows <- length(methods) * expected_rows_per_method
  }
  if (nrow(table) != expected_rows) {
    stop(label, " has an unexpected number of rows.", call. = FALSE)
  }
  invisible(table)
}

load_r1_r2_fashr0143_overlay <- function(
    scenario_id,
    historical_out,
    historical_configuration,
    historical_fash_alpha,
    historical_direct_alpha,
    historical_direct_alpha_005,
    historical_peak_alpha = NULL,
    historical_peak_alpha_005 = NULL,
    cache_dir = R1_R2_FASHR0143_CACHE_DIR) {
  if (!scenario_id %in% c("r1", "r2")) {
    stop("scenario_id must be r1 or r2.", call. = FALSE)
  }
  manifest <- r1_r2_fashr0143_validate_manifest(cache_dir)
  artifact <- readRDS(file.path(
    cache_dir, "replicates", scenario_id, "seed_12345.rds"
  ))
  r1_r2_fashr0143_validate_payload(artifact, historical_out, scenario_id)

  stopifnot(
    identical(historical_configuration$J, manifest$configuration$J),
    identical(
      as.integer(historical_configuration$seed_list),
      R1_R2_FASHR0143_SEEDS
    ),
    identical(
      historical_configuration$common_sd_grid,
      manifest$configuration$common_sd_grid
    ),
    identical(
      historical_configuration$common_pred_step,
      manifest$configuration$pred_step
    ),
    identical(
      as.integer(historical_configuration$common_penalty),
      manifest$configuration$penalty
    ),
    identical(
      artifact$configuration$scenario,
      artifact$scenario
    )
  )

  prefix <- paste0(scenario_id, "_")
  fash_alpha <- r1_r2_fashr0143_read_csv(
    cache_dir, paste0(prefix, "fash_mc_alpha_curve.csv")
  )
  fash_alpha_005 <- r1_r2_fashr0143_read_csv(
    cache_dir, paste0(prefix, "fash_mc_alpha005_summary.csv")
  )
  pi0_summary <- r1_r2_fashr0143_read_csv(
    cache_dir, paste0(prefix, "mc_pi0_summary.csv")
  )
  linear_prior_weights <- r1_r2_fashr0143_read_csv(
    cache_dir, paste0(prefix, "all_replicate_linear_prior_weights.csv")
  )
  linear_prior_summary <- r1_r2_fashr0143_read_csv(
    cache_dir, paste0(prefix, "all_replicate_linear_prior_summary.csv")
  )
  genotype_selection_summary <- r1_r2_fashr0143_read_csv(
    cache_dir, paste0(prefix, "genotype_selection_summary.csv")
  )
  truth_maf_balance <- r1_r2_fashr0143_read_csv(
    cache_dir, paste0(prefix, "truth_maf_balance.csv")
  )

  fash_methods <- c(
    "FASH-IWP1-Raw", "FASH-IWP1-BF",
    "FASH-linear-Raw", "FASH-linear-BF"
  )
  r1_r2_fashr0143_validate_curve(
    fash_alpha, fash_methods, paste0(toupper(scenario_id), " FASH curve")
  )
  r1_r2_fashr0143_validate_curve(
    fash_alpha_005,
    fash_methods,
    paste0(toupper(scenario_id), " FASH alpha-0.05 summary"),
    expected_rows_per_method = 1L
  )

  current_raw <- r1_r2_fashr0143_sort_curve(
    fash_alpha[grepl("-Raw$", fash_alpha$method), , drop = FALSE]
  )
  historical_raw <- r1_r2_fashr0143_sort_curve(
    historical_fash_alpha[
      grepl("-Raw$", historical_fash_alpha$method),
      ,
      drop = FALSE
    ]
  )
  r1_r2_fashr0143_require_equal(
    current_raw,
    historical_raw,
    paste0(toupper(scenario_id), " raw-summary pairing"),
    tolerance = 1e-12,
    check_attributes = FALSE
  )

  direct_rows <- grepl("^Direct-", historical_direct_alpha$method)
  direct_rows_005 <- grepl("^Direct-", historical_direct_alpha_005$method)
  direct_alpha <- rbind(
    fash_alpha[fash_alpha$method %in% c("FASH-IWP1-Raw", "FASH-IWP1-BF"), ],
    historical_direct_alpha[direct_rows, ]
  )
  direct_alpha_005 <- rbind(
    fash_alpha_005[
      fash_alpha_005$method %in% c("FASH-IWP1-Raw", "FASH-IWP1-BF"),
    ],
    historical_direct_alpha_005[direct_rows_005, ]
  )
  direct_methods <- c(
    "FASH-IWP1-Raw",
    "FASH-IWP1-BF",
    "Direct-linear-LRT-eFDR-true-pi0",
    "Direct-quadratic-LRT-eFDR-true-pi0"
  )
  r1_r2_fashr0143_validate_curve(
    direct_alpha,
    direct_methods,
    paste0(toupper(scenario_id), " combined direct curve")
  )
  r1_r2_fashr0143_validate_curve(
    direct_alpha_005,
    direct_methods,
    paste0(toupper(scenario_id), " combined direct alpha-0.05 summary"),
    expected_rows_per_method = 1L
  )

  peak_alpha <- peak_alpha_005 <- NULL
  if (scenario_id == "r2") {
    peak_fash <- r1_r2_fashr0143_read_csv(
      cache_dir, "r2_peak_mc_alpha_curve.csv"
    )
    peak_fash_005 <- r1_r2_fashr0143_read_csv(
      cache_dir, "r2_peak_mc_alpha005_summary.csv"
    )
    peak_direct <- grepl("^Direct-", historical_peak_alpha$method)
    peak_direct_005 <- grepl("^Direct-", historical_peak_alpha_005$method)
    peak_alpha <- rbind(
      peak_fash,
      historical_peak_alpha[peak_direct, ]
    )
    peak_alpha_005 <- rbind(
      peak_fash_005,
      historical_peak_alpha_005[peak_direct_005, ]
    )
    peak_methods <- c(
      fash_methods,
      "Direct-linear-LRT-eFDR-true-pi0",
      "Direct-quadratic-LRT-eFDR-true-pi0"
    )
    r1_r2_fashr0143_validate_curve(
      peak_alpha, peak_methods, "R2 combined peak curve"
    )
    r1_r2_fashr0143_validate_curve(
      peak_alpha_005,
      peak_methods,
      "R2 combined peak alpha-0.05 summary",
      expected_rows_per_method = 1L
    )
  }

  historical_source_label <- paste(
    "Historical", toupper(scenario_id), "formal cache"
  )
  provenance_table <- data.frame(
    `Result component` = c(
      "FASH-IWP1 and FASH-linear, Raw and BF",
      "Direct linear and quadratic interaction eFDR",
      "Representative seed-12345 simulation"
    ),
    `Retained source` = c(
      "Targeted R1/R2 0.1.43 cache",
      historical_source_label,
      historical_source_label
    ),
    `Software provenance` = c(
      paste0(
        "fashr ", R1_R2_FASHR0143_VERSION, " @ ",
        substr(R1_R2_FASHR0143_REMOTE_SHA, 1L, 12L)
      ),
      "Historical direct LRT/eFDR cache; unchanged",
      "Numerically paired to the retained 0.1.43 payload"
    ),
    Status = c(
      "Recomputed on Midway3",
      "Reused; not rerun",
      "Strict pairing check passed"
    ),
    check.names = FALSE
  )

  list(
    manifest = manifest,
    artifact = artifact,
    out = historical_out,
    fash_alpha = fash_alpha,
    fash_alpha_005 = fash_alpha_005,
    direct_alpha = direct_alpha,
    direct_alpha_005 = direct_alpha_005,
    peak_alpha = peak_alpha,
    peak_alpha_005 = peak_alpha_005,
    pi0_summary = pi0_summary,
    linear_prior_weights = linear_prior_weights,
    linear_prior_summary = linear_prior_summary,
    genotype_selection_summary = genotype_selection_summary,
    truth_maf_balance = truth_maf_balance,
    provenance_table = provenance_table
  )
}
