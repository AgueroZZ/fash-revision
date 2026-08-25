#!/usr/bin/env Rscript

# Validate the downloaded formal R3 and R4 fashr 0.1.43 result families.
# This script is read-only: it verifies manifests, provenance, row contracts,
# real-genotype pairing, and numerical diagnostics without modifying caches.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not locate the workflowr project root.", call. = FALSE)
}

expect_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

expect_equal <- function(observed, expected, message, tolerance = 0) {
  if (!isTRUE(all.equal(
    observed,
    expected,
    tolerance = tolerance,
    check.attributes = TRUE
  ))) {
    stop(message, call. = FALSE)
  }
  invisible(TRUE)
}

require_file <- function(path, label) {
  if (!file.exists(path) || dir.exists(path)) {
    stop(label, " is missing: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

sha256_file <- function(path) {
  path <- require_file(path, "SHA-256 input")
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

read_csv_rows <- function(path, expected_rows) {
  path <- require_file(path, "Required summary")
  object <- utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  expect_true(
    identical(nrow(object), as.integer(expected_rows)),
    paste0(
      "Unexpected row count for ", path, ": expected ", expected_rows,
      "; found ", nrow(object), "."
    )
  )
  object
}

expect_finite_numeric_columns <- function(object, label) {
  numeric_columns <- vapply(object, is.numeric, logical(1))
  expect_true(any(numeric_columns), paste0(label, " has no numeric columns."))
  expect_true(
    all(is.finite(as.matrix(object[numeric_columns]))),
    paste0(label, " contains non-finite numeric values.")
  )
  invisible(TRUE)
}

validate_package_provenance <- function(provenance, label) {
  expect_true(is.list(provenance), paste0(label, " provenance is missing."))
  expect_true(
    identical(provenance$package, "fashr") &&
      identical(provenance$version, "0.1.43") &&
      identical(
        provenance$remote_sha,
        "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
      ),
    paste0(label, " does not use the pinned fashr 0.1.43 package.")
  )
  invisible(TRUE)
}

validate_manifest <- function(result_dir, expected_schema, expected_result_id) {
  result_dir <- normalizePath(result_dir, winslash = "/", mustWork = TRUE)
  manifest_path <- require_file(
    file.path(result_dir, "manifest.rds"),
    "Result manifest"
  )
  flag_path <- require_file(
    file.path(result_dir, "complete.flag"),
    "Result completion flag"
  )
  manifest <- readRDS(manifest_path)
  expect_true(
    identical(manifest$schema_version, expected_schema) &&
      identical(manifest$result_id, expected_result_id),
    paste0("Unexpected manifest identity in ", result_dir, ".")
  )
  validate_package_provenance(
    manifest$package_provenance,
    paste0(expected_result_id, " manifest")
  )

  artifact_sha256 <- unlist(manifest$artifact_sha256, use.names = TRUE)
  expect_true(
    is.character(artifact_sha256) &&
      length(artifact_sha256) > 0L &&
      !is.null(names(artifact_sha256)) &&
      !anyDuplicated(names(artifact_sha256)) &&
      all(grepl("^[[:xdigit:]]{64}$", artifact_sha256)),
    paste0(expected_result_id, " has an invalid artifact hash manifest.")
  )
  actual_files <- list.files(
    result_dir,
    recursive = TRUE,
    full.names = FALSE,
    include.dirs = FALSE
  )
  actual_artifacts <- setdiff(actual_files, c("manifest.rds", "complete.flag"))
  expect_true(
    identical(sort(actual_artifacts), sort(names(artifact_sha256))),
    paste0(expected_result_id, " has missing or unmanifested artifact files.")
  )
  observed_sha256 <- vapply(
    names(artifact_sha256),
    function(relative_path) sha256_file(file.path(result_dir, relative_path)),
    character(1)
  )
  expect_true(
    identical(unname(observed_sha256), unname(artifact_sha256)),
    paste0(expected_result_id, " contains an artifact SHA-256 mismatch.")
  )

  list(
    result_dir = result_dir,
    manifest = manifest,
    completion = readLines(flag_path, warn = FALSE)
  )
}

validate_local_sources <- function(manifest, source_paths, label) {
  source_sha256 <- unlist(manifest$source_provenance$sha256, use.names = TRUE)
  expect_true(
    all(names(source_paths) %in% names(source_sha256)),
    paste0(label, " source manifest is incomplete.")
  )
  for (name in names(source_paths)) {
    expect_true(
      identical(sha256_file(source_paths[[name]]), source_sha256[[name]]),
      paste0(label, " source hash mismatch for ", name, ".")
    )
  }
  invisible(TRUE)
}

workflowr_root <- normalizePath(
  find_workflowr_root(),
  winslash = "/",
  mustWork = TRUE
)
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
))

EXPECTED_SEEDS <- c(12345L, 22345L, 32345L, 42345L, 52345L)
EXPECTED_CLASS_PROBS <- c(
  dynamic_bspline = 0.20,
  constant = 0.40,
  zero = 0.40
)
EXPECTED_CLASS_COUNTS <- c(
  dynamic_bspline = 1272L,
  constant = 2545L,
  zero = 2545L
)
EXPECTED_TRUE_PI0 <- sum(EXPECTED_CLASS_COUNTS[c("constant", "zero")]) / 6362L
EXPECTED_GENOTYPE_CONTENT_MD5 <- c(
  `12345` = "526a7318aa2af901e09252f5a6ca3c46",
  `22345` = "517faa30d5218a956f1be84f2567369c",
  `32345` = "9b9be3205d7db54dac31763492bcb2eb",
  `42345` = "9ef73aa94a061df868b1a951fe495d9f",
  `52345` = "7dab159b8453e2f66188ae313bfbd611"
)
EXPECTED_GENOTYPE_CACHE_SHA256 <-
  "81bbef5f323a0bab2ca993c782d8a9b7c63518b83c2cdb46ef7ed1d46f65af24"
EXPECTED_MATRIX_CACHE_SHA256 <-
  "3f63144baa1489533ef13238270c361632d417a7506c9a905aa26e979ad4e1aa"

genotype_cache_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "shared",
  "real_genotype_one_per_gene_J6362_pilot5", "genotype_samples.rds"
)
r1_r2_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "mc", "r1_r2_fashr0143"
)
r1_r2_manifest <- readRDS(require_file(
  file.path(r1_r2_dir, "manifest.rds"),
  "Formal R1/R2 manifest"
))
validate_package_provenance(r1_r2_manifest$package_provenance, "R1/R2 manifest")
expect_true(
  identical(
    r1_r2_manifest$source_provenance$sha256$genotype_cache,
    EXPECTED_GENOTYPE_CACHE_SHA256
  ) && identical(r1_r2_manifest$configuration$J, 6362L),
  "The local formal R1/R2 reference does not match the pinned genotype design."
)
expect_true(
  identical(sha256_file(genotype_cache_path), EXPECTED_GENOTYPE_CACHE_SHA256),
  "The local formal genotype cache SHA-256 is unexpected."
)
genotype_cache <- readRDS(genotype_cache_path)

observed_genotype_digests <- stats::setNames(
  vapply(EXPECTED_SEEDS, function(seed) {
    genotype_sample <- validate_real_genotype_sample(
      genotype_cache$samples[[as.character(seed)]],
      expected_genes = 6362L,
      expected_donors = 19L,
      maf_min = genotype_cache$configuration$maf_min
    )
    genotype_content_md5(
      pair_key = genotype_sample$selection$pair_key,
      sample_ids = rownames(genotype_sample$G),
      G = genotype_sample$G
    )
  }, character(1)),
  as.character(EXPECTED_SEEDS)
)
expect_true(
  identical(observed_genotype_digests, EXPECTED_GENOTYPE_CONTENT_MD5),
  "The canonical genotype-content digests are unexpected."
)

# R3 validation.
R3_RESULT_ID <- paste0(
  "r3_real_genotype_one_per_gene_J6362_",
  paste0(
    "matched_functional_open_middle_3_12_",
    "center_aligned_relative_clearance_main_effect_fashr0143_pilot5"
  )
)
r3_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "mc", R3_RESULT_ID
)
r3_validation <- validate_manifest(
  r3_dir,
  expected_schema = "r3-fashr0143-manifest-v3",
  expected_result_id = R3_RESULT_ID
)
expect_true(
  all(c(
    paste0("result_id=", R3_RESULT_ID),
    "fashr_version=0.1.43",
    "fashr_remote_sha=bf223df75da6e41ae48607a56b4cd12d7c3b24e7",
    "middle_definition=3 < t < 12",
    paste0(
      "raised_cosine_center_ranges=",
      "early:1.5,2.5;middle:4.5,10.5;late:12.5,13.5"
    ),
    "replicates=10"
  ) %in% r3_validation$completion),
  "The R3 completion flag is invalid."
)
validate_local_sources(
  r3_validation$manifest,
  c(
    r3_driver = file.path(
      workflowr_root,
      "code", "revision_simulations", "r3_r4_fashr0143",
      "source_snapshots", "r3_center_aligned_driver.R"
    ),
    simulation_functions = file.path(
      workflowr_root,
      "code", "revision_simulations", "r3_r4_fashr0143",
      "source_snapshots", "r3_center_aligned_simulation_functions.R"
    ),
    real_genotype_helper = file.path(
      workflowr_root,
      "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
    ),
    genotype_cache = genotype_cache_path,
    wrapper_core = file.path(
      workflowr_root,
      "code", "revision_simulations", "r3_r4_fashr0143",
      "16_r3_center_aligned_core.R"
    ),
    runner = file.path(
      workflowr_root,
      "code", "revision_simulations", "r3_r4_fashr0143",
      "16_r3_center_aligned_fashr0143.R"
    )
  ),
  "R3"
)

r3_configuration <- readRDS(file.path(r3_dir, "configuration.rds"))
validate_package_provenance(r3_configuration$package_provenance, "R3 configuration")
expect_true(
  identical(r3_configuration$output_id, R3_RESULT_ID) &&
    identical(r3_configuration$J, 6362L) &&
    identical(r3_configuration$n_donors, 19L) &&
    identical(r3_configuration$n_covariates, 5L) &&
    identical(r3_configuration$middle_window, c(3, 12)) &&
    identical(r3_configuration$middle_boundary, "open") &&
    isTRUE(all.equal(
      r3_configuration$middle_grid,
      seq(3.1, 11.9, by = 0.1)
    )) &&
    identical(r3_configuration$middle_expression, "3 < t < 12") &&
    identical(r3_configuration$seed_list, EXPECTED_SEEDS) &&
    identical(
      r3_configuration$truth_mechanisms,
      c("random_bspline", "raised_cosine")
    ) &&
    identical(
      r3_configuration$raised_cosine$center_ranges,
      list(
        early = c(1.5, 2.5),
        middle = c(4.5, 10.5),
        late = c(12.5, 13.5)
      )
    ) &&
    identical(
      r3_configuration$genotype_digest_method,
      "fash-genotype-content-md5-v1"
    ) &&
    identical(
      r3_configuration$genotype_content_digests,
      EXPECTED_GENOTYPE_CONTENT_MD5
    ),
  "The R3 configuration identity or genotype provenance is invalid."
)
expect_equal(
  r3_configuration$class_probs,
  EXPECTED_CLASS_PROBS,
  "The R3 class probabilities are invalid."
)
expect_equal(
  r3_configuration$expected_class_counts,
  EXPECTED_CLASS_COUNTS,
  "The R3 class counts are invalid."
)
expect_true(
  identical(
    r3_configuration$expected_truth_group_counts,
    stats::setNames(
      rep(212L, 6L),
      c(
        "early / switch", "early / non-switch",
        "middle / switch", "middle / non-switch",
        "late / switch", "late / non-switch"
      )
    )
  ),
  "The R3 functional truth-group counts are invalid."
)

r3_replicate_paths <- unlist(lapply(
  c("random_bspline", "raised_cosine"),
  function(mechanism) file.path(
    r3_dir,
    "replicates",
    paste0(mechanism, "_seed_", EXPECTED_SEEDS, ".rds")
  )
), use.names = FALSE)
expect_true(
  length(r3_replicate_paths) == 10L && all(file.exists(r3_replicate_paths)),
  "The R3 replicate set is incomplete."
)
for (path in r3_replicate_paths) {
  replicate <- readRDS(path)
  seed <- replicate$seed
  mechanism <- replicate$truth_mechanism
  genotype_sample <- genotype_cache$samples[[as.character(seed)]]
  r1_reference <- readRDS(file.path(
    r1_r2_dir,
    "replicates", "r1", paste0("seed_", seed, ".rds")
  ))
  expect_true(
    seed %in% EXPECTED_SEEDS &&
      mechanism %in% c("random_bspline", "raised_cosine") &&
      isTRUE(all.equal(replicate$configuration, r3_configuration)) &&
      identical(
        replicate$genotype_digest,
        EXPECTED_GENOTYPE_CONTENT_MD5[[as.character(seed)]]
      ) &&
      identical(replicate$selected_pair_keys, genotype_sample$selection$pair_key) &&
      identical(replicate$selected_pair_keys, r1_reference$selected_pair_keys) &&
      identical(replicate$genotype_digest, r1_reference$genotype_digest) &&
      nrow(replicate$functional_alpha) == 320L &&
      nrow(replicate$functional_alpha_005) == 8L &&
      nrow(replicate$truth_group_counts) == 6L &&
      all(replicate$truth_group_counts$n_dynamic == 212L),
    paste0("An R3 replicate failed its contract: ", path)
  )
  expect_finite_numeric_columns(
    replicate$functional_alpha,
    paste0("R3 functional alpha replicate ", basename(path))
  )
}

r3_summary_rows <- c(
  all_functional_calls_alpha005.csv = 15512L,
  all_replicate_functional_alpha_curves.csv = 3200L,
  all_replicate_functional_alpha005.csv = 80L,
  all_replicate_pi0.csv = 20L,
  all_truth_group_counts.csv = 60L,
  functional_testing_mc_alpha_curve.csv = 640L,
  functional_testing_mc_alpha005_summary.csv = 16L,
  functional_testing_mc_pi0_summary.csv = 4L,
  genotype_selection_summary.csv = 10L,
  truth_maf_balance.csv = 30L
)
r3_summaries <- lapply(names(r3_summary_rows), function(file_name) {
  read_csv_rows(
    file.path(r3_dir, "summary", file_name),
    r3_summary_rows[[file_name]]
  )
})
names(r3_summaries) <- names(r3_summary_rows)
invisible(lapply(names(r3_summaries), function(file_name) {
  object <- r3_summaries[[file_name]]
  if (identical(file_name, "all_functional_calls_alpha005.csv")) {
    missing_main_effect <- is.na(object$genetic_main_effect)
    expect_true(
      all(
        object$truth_group[missing_main_effect] == "dynamic-null" &
          object$effect_class[missing_main_effect] %in% c("constant", "zero")
      ) &&
        all(is.finite(object$genetic_main_effect[!missing_main_effect])),
      paste0(
        "R3 ", file_name,
        " has unexpected missing genetic-main-effect annotations."
      )
    )
    object$genetic_main_effect <- NULL
  }
  expect_finite_numeric_columns(object, paste0("R3 ", file_name))
}))
cat("Validated downloaded R3 fashr 0.1.43 results.\n")

# R4 validation.
R4_RESULT_ID <- "r4_fashr0143"
r4_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "mc", R4_RESULT_ID
)
r4_validation <- validate_manifest(
  r4_dir,
  expected_schema = "r4-fashr0143-manifest-v1",
  expected_result_id = R4_RESULT_ID
)
expect_true(
  all(c(
    paste0("result_id=", R4_RESULT_ID),
    "fashr_version=0.1.43",
    "fashr_remote_sha=bf223df75da6e41ae48607a56b4cd12d7c3b24e7",
    "full_matrix_replicates=5",
    "lag1_sweep_replicates=5"
  ) %in% r4_validation$completion),
  "The R4 completion flag is invalid."
)

r4_real_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "real_data",
  "r4_null_like_top500_full_correlations_fashr0143"
)
r4_matrix_cache_path <- file.path(
  r4_real_dir,
  "simulation_correlation_matrices.rds"
)
validate_local_sources(
  r4_validation$manifest,
  c(
    simulation_functions = file.path(
      workflowr_root,
      "code", "revision_simulations", "r3_r4_fashr0143",
      "source_snapshots", "r4_simulation_functions.R"
    ),
    real_genotype_helper = file.path(
      workflowr_root,
      "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
    ),
    r4_real_genotype_helper = file.path(
      workflowr_root,
      "code", "revision_simulations", "r4_correlated_errors",
      "real_genotype_r1_helpers.R"
    ),
    full_driver = file.path(
      workflowr_root,
      "code", "revision_simulations", "r4_correlated_errors",
      "run_full_empirical_correlation_mc.R"
    ),
    sweep_driver = file.path(
      workflowr_root,
      "code", "revision_simulations", "r4_correlated_errors",
      "run_lag1_correlation_sweep.R"
    ),
    plotting = file.path(
      workflowr_root,
      "code", "revision_simulations", "r4_correlated_errors", "plotting.R"
    ),
    genotype_cache = genotype_cache_path,
    matrix_cache = r4_matrix_cache_path,
    runner = file.path(
      workflowr_root,
      "code", "revision_simulations", "r3_r4_fashr0143", "14_r4_fashr0143.R"
    )
  ),
  "R4"
)
expect_true(
  identical(sha256_file(r4_matrix_cache_path), EXPECTED_MATRIX_CACHE_SHA256),
  "The local R4 correlation-matrix cache SHA-256 is unexpected."
)

real_analysis <- readRDS(require_file(
  file.path(r4_real_dir, "real_data_correlation_analysis.rds"),
  "Versioned R4 real-data correlation analysis"
))
validate_package_provenance(
  real_analysis$configuration$package_provenance,
  "R4 real-data correlation input"
)
expect_true(
  identical(
    real_analysis$configuration$output_id,
    "r4_null_like_top500_full_correlations_fashr0143"
  ) &&
    identical(
      real_analysis$configuration$fit_sha256,
      "7f0ca9ab0fbeab89a13c83d2a0fb7c24195f7b5a5835f209399cf0e359001f50"
    ) &&
    identical(real_analysis$configuration$top_n, 500L) &&
    identical(real_analysis$configuration$n_gene_representatives, 6362L) &&
    identical(real_analysis$configuration$n_bootstrap, 1000L) &&
    identical(real_analysis$configuration$n_independence_benchmark, 1000L),
  "The versioned R4 real-data correlation input is invalid."
)
for (matrix_id in c("direct_adjusted", "pairwise_adjusted")) {
  matrix <- real_analysis$matrices[[matrix_id]]$projected
  expect_true(
    identical(dim(matrix), c(16L, 16L)) &&
      all(is.finite(matrix)) &&
      max(abs(diag(matrix) - 1)) <= 1e-10 &&
      min(eigen(matrix, symmetric = TRUE, only.values = TRUE)$values) > 0,
    paste0("The R4 real-data matrix is invalid: ", matrix_id)
  )
}

validate_r4_configuration <- function(configuration, output_id) {
  validate_package_provenance(
    configuration$package_provenance,
    paste0(output_id, " configuration")
  )
  expect_true(
    identical(configuration$output_id, output_id) &&
      identical(
        configuration$scenario,
        "r1_real_genotype_one_per_gene_random_bspline_main_effect_dynamic_eqtl"
      ) &&
      identical(configuration$J, 6362L) &&
      identical(configuration$n_donors, 19L) &&
      identical(configuration$n_covariates, 5L) &&
      identical(configuration$seed_list, EXPECTED_SEEDS) &&
      identical(configuration$r1_reference_output_id, "r1_r2_fashr0143") &&
      identical(configuration$genotype_source, "paper-derived YRI DS dosage") &&
      identical(
        configuration$genotype_selection_rule,
        "one uniformly sampled tested variant per gene"
      ),
    paste0(output_id, " has an invalid formal R1 design contract.")
  )
  expect_equal(
    configuration$class_probs,
    EXPECTED_CLASS_PROBS,
    paste0(output_id, " has invalid class probabilities.")
  )
  expect_equal(
    configuration$expected_class_counts,
    EXPECTED_CLASS_COUNTS,
    paste0(output_id, " has invalid class counts.")
  )
  expect_equal(
    configuration$true_pi0,
    EXPECTED_TRUE_PI0,
    paste0(output_id, " has an invalid exact true pi0."),
    tolerance = 1e-15
  )
  invisible(TRUE)
}

full_dir <- file.path(r4_dir, "full_matrix")
sweep_dir <- file.path(r4_dir, "lag1_sweep")
full_configuration <- readRDS(file.path(full_dir, "configuration.rds"))
sweep_configuration <- readRDS(file.path(sweep_dir, "configuration.rds"))
validate_r4_configuration(
  full_configuration,
  "r4_full_empirical_correlations_top500_J6362_fashr0143_pilot5"
)
validate_r4_configuration(
  sweep_configuration,
  "r4_lag1_correlation_sweep_J6362_m0p3_to_p0p3_fashr0143_pilot5"
)
expect_true(
  identical(
    full_configuration$conditions,
    c(
      "Independent",
      "Direct centered full matrix",
      "Pairwise-difference full matrix"
    )
  ) &&
    identical(
      full_configuration$methods,
      c("FASH-IWP1-Raw", "FASH-IWP1-BF")
    ),
  "The R4 full-matrix condition or method contract is invalid."
)
for (condition in names(full_configuration$condition_matrices)) {
  matrix <- full_configuration$condition_matrices[[condition]]
  expect_true(
    identical(dim(matrix), c(16L, 16L)) &&
      all(is.finite(matrix)) &&
      max(abs(diag(matrix) - 1)) <= 1e-10 &&
      min(eigen(matrix, symmetric = TRUE, only.values = TRUE)$values) > 0,
    paste0("The R4 target matrix is invalid: ", condition)
  )
}
expect_equal(
  sweep_configuration$rho_grid,
  seq(-0.3, 0.3, by = 0.1),
  "The R4 lag-1 grid is invalid.",
  tolerance = 1e-12
)
expect_true(
  identical(sweep_configuration$correlation_structure, "lag1_only") &&
    identical(
      sweep_configuration$pi0_methods,
      c("FASH-IWP1-Raw", "FASH-IWP1-BF")
    ),
  "The R4 lag-1 structure or pi0 method contract is invalid."
)

full_summary_rows <- c(
  all_condition_alpha005.csv = 30L,
  all_condition_pi0.csv = 30L,
  all_realized_correlation_matrices.csv = 11520L,
  all_realized_lag_summaries.csv = 675L,
  condition_mc_alpha005_summary.csv = 6L,
  condition_mc_pi0_summary.csv = 6L,
  matrix_match_diagnostics.csv = 45L,
  paired_alpha005_differences.csv = 4L,
  paired_pi0_differences.csv = 4L,
  pairing_check.csv = 15L,
  r1_reference_check.csv = 2L,
  realized_correlation_matrix_summary.csv = 2304L,
  realized_lag_summary.csv = 135L,
  target_matrix_checks.csv = 3L
)
full_summaries <- lapply(names(full_summary_rows), function(file_name) {
  read_csv_rows(
    file.path(full_dir, "summary", file_name),
    full_summary_rows[[file_name]]
  )
})
names(full_summaries) <- names(full_summary_rows)
invisible(lapply(names(full_summaries), function(file_name) {
  expect_finite_numeric_columns(
    full_summaries[[file_name]],
    paste0("R4 full-matrix ", file_name)
  )
}))
expect_true(
  all(full_summaries[["pairing_check.csv"]]$passed) &&
    all(full_summaries[["r1_reference_check.csv"]]$passed) &&
    all(full_summaries[["target_matrix_checks.csv"]]$positive_definite) &&
    all(full_summaries[["target_matrix_checks.csv"]]$minimum_eigenvalue > 0),
  "The R4 full-matrix pairing, R1 reference, or matrix checks failed."
)

sweep_summary_rows <- c(
  all_alpha005.csv = 140L,
  all_lag1_correlations.csv = 70L,
  all_pi0.csv = 70L,
  correlation_matrix_check.csv = 7L,
  lag1_correlation_summary.csv = 14L,
  mc_alpha005_summary.csv = 28L,
  mc_pi0_summary.csv = 14L,
  paired_vs_zero_alpha005_summary.csv = 28L,
  paired_vs_zero_pi0_summary.csv = 14L,
  pairing_check.csv = 35L,
  r1_pi0_reference_check.csv = 1L,
  r1_reference_check.csv = 1L
)
sweep_summaries <- lapply(names(sweep_summary_rows), function(file_name) {
  read_csv_rows(
    file.path(sweep_dir, "summary", file_name),
    sweep_summary_rows[[file_name]]
  )
})
names(sweep_summaries) <- names(sweep_summary_rows)
invisible(lapply(names(sweep_summaries), function(file_name) {
  expect_finite_numeric_columns(
    sweep_summaries[[file_name]],
    paste0("R4 lag-1 sweep ", file_name)
  )
}))
expect_true(
  all(sweep_summaries[["pairing_check.csv"]]$passed) &&
    all(sweep_summaries[["r1_reference_check.csv"]]$passed) &&
    all(sweep_summaries[["r1_pi0_reference_check.csv"]]$passed) &&
    all(sweep_summaries[["correlation_matrix_check.csv"]]$positive_definite) &&
    all(sweep_summaries[["correlation_matrix_check.csv"]]$minimum_eigenvalue > 0),
  "The R4 sweep pairing, R1 reference, or matrix checks failed."
)

for (seed in EXPECTED_SEEDS) {
  genotype_sample <- genotype_cache$samples[[as.character(seed)]]
  full_replicate <- readRDS(file.path(
    full_dir, "replicates", paste0("seed_", seed, ".rds")
  ))
  sweep_replicate <- readRDS(file.path(
    sweep_dir, "replicates", paste0("seed_", seed, ".rds")
  ))
  expect_true(
    identical(full_replicate$seed, seed) &&
      identical(sweep_replicate$seed, seed) &&
      identical(full_replicate$selected_pair_keys, genotype_sample$selection$pair_key) &&
      identical(sweep_replicate$selected_pair_keys, genotype_sample$selection$pair_key) &&
      isTRUE(all.equal(full_replicate$configuration, full_configuration)) &&
      isTRUE(all.equal(sweep_replicate$configuration, sweep_configuration)) &&
      nrow(full_replicate$alpha_005) == 6L &&
      nrow(full_replicate$pi0) == 6L &&
      nrow(full_replicate$pairing_check) == 3L &&
      all(full_replicate$pairing_check$passed) &&
      nrow(sweep_replicate$alpha_005) == 28L &&
      nrow(sweep_replicate$pi0) == 14L &&
      nrow(sweep_replicate$lag1_diagnostics) == 14L &&
      nrow(sweep_replicate$pairing_check) == 7L &&
      all(sweep_replicate$pairing_check$passed),
    paste0("The R4 replicate contract failed for seed ", seed, ".")
  )
  expect_finite_numeric_columns(
    full_replicate$alpha_005,
    paste0("R4 full-matrix alpha replicate ", seed)
  )
  expect_finite_numeric_columns(
    sweep_replicate$alpha_005,
    paste0("R4 sweep alpha replicate ", seed)
  )
}

cat("Validated downloaded R4 fashr 0.1.43 results.\n")
cat("All downloaded R3/R4 manifests and scientific contracts passed.\n")
