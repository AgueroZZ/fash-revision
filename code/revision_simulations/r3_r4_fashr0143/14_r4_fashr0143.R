#!/usr/bin/env Rscript

# Run and validate the formal R4 temporal-correlation sensitivity aligned with
# the current J=6362 real-genotype R1 design under fashr 0.1.43.

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[[1L]] == length(args)) return(default)
  args[[hit[[1L]] + 1L]]
}

as_flag <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

log_message <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    ...,
    "\n",
    sep = ""
  )
  flush.console()
}

require_file <- function(path, label) {
  if (!file.exists(path) || dir.exists(path)) {
    stop(label, " is missing: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

ensure_directory <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) {
    stop("Unable to create directory: ", path, call. = FALSE)
  }
  invisible(path)
}

sha256_file <- function(path) {
  path <- require_file(path, "SHA-256 input")
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

atomic_save_rds <- function(object, path) {
  temporary_path <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)
  saveRDS(object, temporary_path, version = 3)
  if (!file.rename(temporary_path, path)) {
    stop("Unable to atomically replace ", path, ".", call. = FALSE)
  }
  invisible(path)
}

atomic_write_lines <- function(text, path) {
  temporary_path <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)
  writeLines(text, temporary_path, useBytes = TRUE)
  if (!file.rename(temporary_path, path)) {
    stop("Unable to atomically replace ", path, ".", call. = FALSE)
  }
  invisible(path)
}

read_required_csv <- function(path, expected_rows) {
  path <- require_file(path, "R4 summary")
  object <- utils::read.csv(path, stringsAsFactors = FALSE)
  if (!identical(nrow(object), as.integer(expected_rows))) {
    stop(
      "Unexpected row count for ", path, ": expected ", expected_rows,
      "; found ", nrow(object), ".",
      call. = FALSE
    )
  }
  object
}

EXPECTED_FASHR_VERSION <- "0.1.43"
EXPECTED_FASHR_REMOTE_SHA <-
  "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
EXPECTED_SIMULATION_FUNCTIONS_SHA256 <-
  "93f9a2c5606ae74763fb53e189af62c4d3b7c973aaf38d5b3eb3bf32f97a6487"
EXPECTED_REAL_GENOTYPE_HELPER_SHA256 <-
  "c03c01a188503336a77793c96f4e2d3ac7e0cbd56f4028b552da4b2f88e6b9d7"
EXPECTED_R4_REAL_GENOTYPE_HELPER_SHA256 <-
  "3f39a3c90578552d8d8a9af097748e2a95f4c60db43e01d04816732726685dd8"
EXPECTED_FULL_DRIVER_SHA256 <-
  "8dc4d9c4df52326efddca255b337b46021419d6ce6414f7435eced710246e75e"
EXPECTED_SWEEP_DRIVER_SHA256 <-
  "d7124291f1929ad050cf97e64f9e5f754427973a3226e0a2cbd4be54e0934d85"
EXPECTED_PLOTTING_SHA256 <-
  "8455eb799a976ab4de3e14dcc0a3cac05a1a65a04f245ff27688807d8deaac5d"
EXPECTED_GENOTYPE_CACHE_SHA256 <-
  "81bbef5f323a0bab2ca993c782d8a9b7c63518b83c2cdb46ef7ed1d46f65af24"
EXPECTED_MATRIX_CACHE_SHA256 <-
  "3f63144baa1489533ef13238270c361632d417a7506c9a905aa26e979ad4e1aa"
EXPECTED_GENOTYPE_CONTENT_MD5 <- c(
  `12345` = "526a7318aa2af901e09252f5a6ca3c46",
  `22345` = "517faa30d5218a956f1be84f2567369c",
  `32345` = "9b9be3205d7db54dac31763492bcb2eb",
  `42345` = "9ef73aa94a061df868b1a951fe495d9f",
  `52345` = "7dab159b8453e2f66188ae313bfbd611"
)

RESULT_ID <- "r4_fashr0143"
FULL_OUTPUT_ID <-
  "r4_full_empirical_correlations_top500_J6362_fashr0143_pilot5"
SWEEP_OUTPUT_ID <-
  "r4_lag1_correlation_sweep_J6362_m0p3_to_p0p3_fashr0143_pilot5"
SEED_LIST <- c(12345L, 22345L, 32345L, 42345L, 52345L)
RHO_GRID <- seq(-0.3, 0.3, by = 0.1)
J <- 6362L
TRUE_PI0 <- 5090 / 6362
NUM_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
PREFLIGHT_ONLY <- as_flag(get_arg("--preflight-only", "false"))

if (is.na(NUM_CORES) || NUM_CORES < 1L) {
  stop("SLURM_CPUS_PER_TASK must be a positive integer.", call. = FALSE)
}

WORKSPACE_DIR <- Sys.getenv(
  "FASH_WORKSPACE",
  unset = "/project/mstephens/ziangzhang/fash/workspace"
)
WORKSPACE_DIR <- normalizePath(WORKSPACE_DIR, winslash = "/", mustWork = TRUE)
RESULT_PARENT <- Sys.getenv(
  "FASH_R4_RESULT_PARENT",
  unset = file.path(WORKSPACE_DIR, "results", "revision_simulations")
)
FINAL_DIR <- file.path(RESULT_PARENT, RESULT_ID)
PARTIAL_DIR <- file.path(RESULT_PARENT, paste0(RESULT_ID, "_partial"))
FULL_OUTPUT_DIR <- file.path(PARTIAL_DIR, "full_matrix")
SWEEP_OUTPUT_DIR <- file.path(PARTIAL_DIR, "lag1_sweep")

SIMULATION_FUNCTIONS_PATH <- require_file(
  Sys.getenv(
    "FASH_R4_SIMULATION_FUNCTIONS",
    unset = file.path(
      WORKSPACE_DIR,
      "code", "revision_simulations", "shared", "simulation_functions.R"
    )
  ),
  "Shared simulation functions"
)
REAL_GENOTYPE_HELPER_PATH <- require_file(
  Sys.getenv(
    "FASH_R4_REAL_GENOTYPE_HELPER",
    unset = file.path(
      WORKSPACE_DIR,
      "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
    )
  ),
  "Real-genotype helper"
)
R4_REAL_GENOTYPE_HELPER_PATH <- require_file(
  Sys.getenv(
    "FASH_R4_R1_HELPER",
    unset = file.path(
      WORKSPACE_DIR,
      "code", "revision_simulations", "r4_correlated_errors",
      "real_genotype_r1_helpers.R"
    )
  ),
  "R4 formal-R1 helper"
)
FULL_DRIVER_PATH <- require_file(
  Sys.getenv(
    "FASH_R4_FULL_DRIVER",
    unset = file.path(
      WORKSPACE_DIR,
      "code", "revision_simulations", "r4_correlated_errors",
      "run_full_empirical_correlation_mc.R"
    )
  ),
  "R4 full-matrix driver"
)
SWEEP_DRIVER_PATH <- require_file(
  Sys.getenv(
    "FASH_R4_SWEEP_DRIVER",
    unset = file.path(
      WORKSPACE_DIR,
      "code", "revision_simulations", "r4_correlated_errors",
      "run_lag1_correlation_sweep.R"
    )
  ),
  "R4 lag-1 sweep driver"
)
PLOTTING_PATH <- require_file(
  Sys.getenv(
    "FASH_R4_PLOTTING",
    unset = file.path(
      WORKSPACE_DIR,
      "code", "revision_simulations", "r4_correlated_errors", "plotting.R"
    )
  ),
  "R4 plotting helper"
)
GENOTYPE_CACHE_PATH <- require_file(
  Sys.getenv(
    "FASH_R4_GENOTYPE_CACHE",
    unset = file.path(
      WORKSPACE_DIR,
      "inputs", "r1_r2_fashr0143", "genotype_samples.rds"
    )
  ),
  "Formal genotype cache"
)
MATRIX_CACHE_PATH <- require_file(
  Sys.getenv(
    "FASH_R4_MATRIX_CACHE",
    unset = file.path(
      WORKSPACE_DIR,
      "inputs", "r3_r4_fashr0143", "simulation_correlation_matrices.rds"
    )
  ),
  "R4 full-correlation matrix cache"
)
R1_CACHE_DIR <- Sys.getenv(
  "FASH_R4_R1_CACHE_DIR",
  unset = file.path(
    WORKSPACE_DIR,
    "results", "revision_simulations", "r1_r2_fashr0143"
  )
)
R1_CACHE_DIR <- normalizePath(R1_CACHE_DIR, winslash = "/", mustWork = TRUE)
R1_MANIFEST_PATH <- require_file(
  file.path(R1_CACHE_DIR, "manifest.rds"),
  "R1 0.1.43 manifest"
)
invisible(require_file(
  file.path(R1_CACHE_DIR, "complete.flag"),
  "R1 completion flag"
))
RUNNER_PATH <- normalizePath(
  sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][[1L]]),
  winslash = "/",
  mustWork = TRUE
)

observed_source_hashes <- list(
  simulation_functions = sha256_file(SIMULATION_FUNCTIONS_PATH),
  real_genotype_helper = sha256_file(REAL_GENOTYPE_HELPER_PATH),
  r4_real_genotype_helper = sha256_file(R4_REAL_GENOTYPE_HELPER_PATH),
  full_driver = sha256_file(FULL_DRIVER_PATH),
  sweep_driver = sha256_file(SWEEP_DRIVER_PATH),
  plotting = sha256_file(PLOTTING_PATH),
  genotype_cache = sha256_file(GENOTYPE_CACHE_PATH),
  matrix_cache = sha256_file(MATRIX_CACHE_PATH),
  runner = sha256_file(RUNNER_PATH)
)
expected_source_hashes <- list(
  simulation_functions = EXPECTED_SIMULATION_FUNCTIONS_SHA256,
  real_genotype_helper = EXPECTED_REAL_GENOTYPE_HELPER_SHA256,
  r4_real_genotype_helper = EXPECTED_R4_REAL_GENOTYPE_HELPER_SHA256,
  full_driver = EXPECTED_FULL_DRIVER_SHA256,
  sweep_driver = EXPECTED_SWEEP_DRIVER_SHA256,
  plotting = EXPECTED_PLOTTING_SHA256,
  genotype_cache = EXPECTED_GENOTYPE_CACHE_SHA256,
  matrix_cache = EXPECTED_MATRIX_CACHE_SHA256
)
for (name in names(expected_source_hashes)) {
  if (!identical(observed_source_hashes[[name]], expected_source_hashes[[name]])) {
    stop(
      "Unexpected SHA-256 for ", name, ": expected ",
      expected_source_hashes[[name]], "; found ", observed_source_hashes[[name]],
      ".",
      call. = FALSE
    )
  }
}

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is not installed.", call. = FALSE)
}
fashr_description <- utils::packageDescription("fashr")
package_provenance <- list(
  package = "fashr",
  version = as.character(utils::packageVersion("fashr")),
  remote_sha = if (is.null(fashr_description$RemoteSha)) {
    NA_character_
  } else {
    as.character(fashr_description$RemoteSha)
  },
  library_path = normalizePath(
    find.package("fashr"), winslash = "/", mustWork = TRUE
  ),
  r_version = R.version.string,
  platform = R.version$platform
)
if (!identical(package_provenance$version, EXPECTED_FASHR_VERSION) ||
    !identical(package_provenance$remote_sha, EXPECTED_FASHR_REMOTE_SHA)) {
  stop(
    "Expected fashr ", EXPECTED_FASHR_VERSION, " at ",
    EXPECTED_FASHR_REMOTE_SHA, "; found ", package_provenance$version,
    " at ", package_provenance$remote_sha, ".",
    call. = FALSE
  )
}

source(SIMULATION_FUNCTIONS_PATH)
source(REAL_GENOTYPE_HELPER_PATH)
source(R4_REAL_GENOTYPE_HELPER_PATH)
genotype_cache <- validate_r4_real_genotype_cache(
  readRDS(GENOTYPE_CACHE_PATH),
  seed_list = SEED_LIST,
  J = J,
  n_donors = 19L
)
for (seed in SEED_LIST) {
  sample <- genotype_cache$samples[[as.character(seed)]]
  observed_digest <- genotype_content_md5(
    pair_key = sample$selection$pair_key,
    sample_ids = rownames(sample$G),
    G = sample$G
  )
  if (!identical(
    observed_digest,
    EXPECTED_GENOTYPE_CONTENT_MD5[[as.character(seed)]]
  )) {
    stop("Unexpected genotype-content digest for seed ", seed, ".", call. = FALSE)
  }
}

matrix_cache <- readRDS(MATRIX_CACHE_PATH)
expected_matrix_names <- c("direct_centered", "pairwise_difference")
if (!identical(sort(names(matrix_cache)), sort(expected_matrix_names))) {
  stop("The R4 matrix cache has unexpected entries.", call. = FALSE)
}
for (name in expected_matrix_names) {
  matrix <- validate_time_correlation(matrix_cache[[name]], n_time = 16L)
  if (min(eigen(matrix, symmetric = TRUE, only.values = TRUE)$values) <= 0) {
    stop("The R4 matrix is not positive definite: ", name, call. = FALSE)
  }
}

r1_manifest <- readRDS(R1_MANIFEST_PATH)
if (!identical(r1_manifest$result_id, "r1_r2_fashr0143") ||
    !identical(
      r1_manifest$package_provenance$version,
      EXPECTED_FASHR_VERSION
    ) ||
    !identical(
      r1_manifest$package_provenance$remote_sha,
      EXPECTED_FASHR_REMOTE_SHA
    ) ||
    !identical(r1_manifest$configuration$J, J)) {
  stop("The R1 reference manifest is invalid.", call. = FALSE)
}
for (relative_path in c(
  "summary/r1_all_replicate_fash_alpha_curves.csv",
  "summary/r1_all_replicate_pi0.csv"
)) {
  expected_hash <- r1_manifest$summary_sha256[[relative_path]]
  observed_hash <- sha256_file(file.path(R1_CACHE_DIR, relative_path))
  if (is.null(expected_hash) || !identical(observed_hash, expected_hash)) {
    stop("The R1 reference summary failed its manifest hash: ", relative_path)
  }
}

log_message("R4 result: ", FINAL_DIR)
log_message("fashr ", package_provenance$version, " at ", package_provenance$remote_sha)
log_message("Validated R4 sources, matrices, genotype cache, and R1 reference.")

if (PREFLIGHT_ONLY) {
  log_message("R4 preflight completed successfully.")
  quit(save = "no", status = 0L)
}
if (dir.exists(FINAL_DIR) || file.exists(FINAL_DIR)) {
  stop("Refusing to overwrite completed R4 output: ", FINAL_DIR, call. = FALSE)
}
ensure_directory(RESULT_PARENT)
ensure_directory(PARTIAL_DIR)

old_working_directory <- setwd(WORKSPACE_DIR)
on.exit(setwd(old_working_directory), add = TRUE)
rscript <- file.path(R.home("bin"), "Rscript")
common_arguments <- c(
  "--J", as.character(J),
  "--n-donors", "19",
  "--n-covariates", "5",
  "--noise-sd", "1",
  "--seed-list", paste(SEED_LIST, collapse = ","),
  "--num-cores", as.character(NUM_CORES),
  "--num-basis", "20",
  "--genotype-cache", GENOTYPE_CACHE_PATH,
  "--r1-cache-dir", R1_CACHE_DIR,
  "--expected-fashr-version", EXPECTED_FASHR_VERSION,
  "--expected-fashr-remote-sha", EXPECTED_FASHR_REMOTE_SHA
)

log_message("Starting or resuming the R4 complete-matrix experiment.")
full_status <- system2(
  rscript,
  c(
    "--vanilla", FULL_DRIVER_PATH, common_arguments,
    "--matrix-cache", MATRIX_CACHE_PATH,
    "--output-id", FULL_OUTPUT_ID,
    "--output-dir", FULL_OUTPUT_DIR
  )
)
if (!identical(full_status, 0L)) {
  stop("The R4 full-matrix driver failed with status ", full_status, ".")
}

log_message("Starting or resuming the R4 lag-1 sweep.")
sweep_status <- system2(
  rscript,
  c(
    "--vanilla", SWEEP_DRIVER_PATH, common_arguments,
    "--rho-list", paste(RHO_GRID, collapse = ","),
    "--output-id", SWEEP_OUTPUT_ID,
    "--output-dir", SWEEP_OUTPUT_DIR
  )
)
if (!identical(sweep_status, 0L)) {
  stop("The R4 lag-1 driver failed with status ", sweep_status, ".")
}

full_configuration <- readRDS(require_file(
  file.path(FULL_OUTPUT_DIR, "configuration.rds"),
  "R4 full-matrix configuration"
))
sweep_configuration <- readRDS(require_file(
  file.path(SWEEP_OUTPUT_DIR, "configuration.rds"),
  "R4 sweep configuration"
))
for (configuration in list(full_configuration, sweep_configuration)) {
  if (!identical(configuration$J, J) ||
      !identical(configuration$seed_list, SEED_LIST) ||
      !isTRUE(all.equal(configuration$true_pi0, TRUE_PI0, tolerance = 0)) ||
      !identical(
        configuration$package_provenance$version,
        EXPECTED_FASHR_VERSION
      ) ||
      !identical(
        configuration$package_provenance$remote_sha,
        EXPECTED_FASHR_REMOTE_SHA
      ) ||
      !identical(configuration$r1_reference_output_id, "r1_r2_fashr0143")) {
    stop("An R4 configuration failed the formal contract.", call. = FALSE)
  }
}

for (output_dir in c(FULL_OUTPUT_DIR, SWEEP_OUTPUT_DIR)) {
  replicate_paths <- file.path(
    output_dir,
    "replicates",
    paste0("seed_", SEED_LIST, ".rds")
  )
  for (path in replicate_paths) require_file(path, "R4 replicate")
  replicates <- lapply(replicate_paths, readRDS)
  for (index in seq_along(replicates)) {
    replicate <- replicates[[index]]
    seed <- SEED_LIST[[index]]
    if (!identical(replicate$seed, seed) ||
        !identical(
          replicate$selected_pair_keys,
          genotype_cache$samples[[as.character(seed)]]$selection$pair_key
        ) ||
        !identical(
          replicate$genotype_digest,
          genotype_cache$samples[[as.character(seed)]]$genotype_digest
        ) ||
        any(!replicate$pairing_check$passed) ||
        any(replicate$pairing_check$maximum_absolute_error_difference >
          replicate$pairing_check$tolerance)) {
      stop("An R4 replicate failed pairing or genotype validation: ", path)
    }
  }
}

full_summary <- file.path(FULL_OUTPUT_DIR, "summary")
full_alpha <- read_required_csv(
  file.path(full_summary, "all_condition_alpha005.csv"), 30L
)
read_required_csv(
  file.path(full_summary, "condition_mc_alpha005_summary.csv"), 6L
)
read_required_csv(
  file.path(full_summary, "paired_alpha005_differences.csv"), 4L
)
full_pi0 <- read_required_csv(
  file.path(full_summary, "all_condition_pi0.csv"), 30L
)
read_required_csv(
  file.path(full_summary, "condition_mc_pi0_summary.csv"), 6L
)
read_required_csv(
  file.path(full_summary, "paired_pi0_differences.csv"), 4L
)
read_required_csv(
  file.path(full_summary, "all_realized_correlation_matrices.csv"), 11520L
)
read_required_csv(
  file.path(full_summary, "realized_correlation_matrix_summary.csv"), 2304L
)
read_required_csv(
  file.path(full_summary, "all_realized_lag_summaries.csv"), 675L
)
read_required_csv(
  file.path(full_summary, "realized_lag_summary.csv"), 135L
)
read_required_csv(
  file.path(full_summary, "matrix_match_diagnostics.csv"), 45L
)
full_pairing <- read_required_csv(
  file.path(full_summary, "pairing_check.csv"), 15L
)
full_matrix_checks <- read_required_csv(
  file.path(full_summary, "target_matrix_checks.csv"), 3L
)
full_reference <- read_required_csv(
  file.path(full_summary, "r1_reference_check.csv"), 2L
)
if (any(!is.finite(full_alpha$power)) ||
    any(!is.finite(full_alpha$empirical_fdr)) ||
    any(!is.finite(full_pi0$estimated_pi0)) ||
    any(!full_pairing$passed) ||
    any(!full_matrix_checks$positive_definite) ||
    any(!full_reference$passed)) {
  stop("The R4 full-matrix summaries failed numerical validation.")
}

sweep_summary <- file.path(SWEEP_OUTPUT_DIR, "summary")
sweep_alpha <- read_required_csv(
  file.path(sweep_summary, "all_alpha005.csv"), 140L
)
read_required_csv(file.path(sweep_summary, "mc_alpha005_summary.csv"), 28L)
read_required_csv(
  file.path(sweep_summary, "paired_vs_zero_alpha005_summary.csv"), 28L
)
sweep_pi0 <- read_required_csv(file.path(sweep_summary, "all_pi0.csv"), 70L)
read_required_csv(file.path(sweep_summary, "mc_pi0_summary.csv"), 14L)
read_required_csv(
  file.path(sweep_summary, "paired_vs_zero_pi0_summary.csv"), 14L
)
read_required_csv(
  file.path(sweep_summary, "all_lag1_correlations.csv"), 70L
)
read_required_csv(
  file.path(sweep_summary, "lag1_correlation_summary.csv"), 14L
)
sweep_pairing <- read_required_csv(
  file.path(sweep_summary, "pairing_check.csv"), 35L
)
sweep_matrix_checks <- read_required_csv(
  file.path(sweep_summary, "correlation_matrix_check.csv"), 7L
)
sweep_reference <- read_required_csv(
  file.path(sweep_summary, "r1_reference_check.csv"), 1L
)
sweep_pi0_reference <- read_required_csv(
  file.path(sweep_summary, "r1_pi0_reference_check.csv"), 1L
)
if (any(!is.finite(sweep_alpha$power)) ||
    any(!is.finite(sweep_alpha$empirical_fdr)) ||
    any(!is.finite(sweep_pi0$estimated_pi0)) ||
    any(!sweep_pairing$passed) ||
    any(!sweep_matrix_checks$positive_definite) ||
    any(!sweep_reference$passed) ||
    any(!sweep_pi0_reference$passed)) {
  stop("The R4 sweep summaries failed numerical validation.")
}

artifact_paths <- list.files(
  PARTIAL_DIR,
  recursive = TRUE,
  full.names = TRUE,
  all.files = FALSE
)
artifact_paths <- artifact_paths[file.info(artifact_paths)$isdir %in% FALSE]
partial_normalized <- normalizePath(PARTIAL_DIR, winslash = "/", mustWork = TRUE)
artifact_relative_paths <- substring(
  artifact_paths,
  nchar(partial_normalized) + 2L
)
artifact_sha256 <- stats::setNames(
  vapply(artifact_paths, sha256_file, character(1)),
  artifact_relative_paths
)
manifest <- list(
  schema_version = "r4-fashr0143-manifest-v1",
  result_id = RESULT_ID,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  package_provenance = package_provenance,
  source_provenance = list(
    paths = list(
      simulation_functions = SIMULATION_FUNCTIONS_PATH,
      real_genotype_helper = REAL_GENOTYPE_HELPER_PATH,
      r4_real_genotype_helper = R4_REAL_GENOTYPE_HELPER_PATH,
      full_driver = FULL_DRIVER_PATH,
      sweep_driver = SWEEP_DRIVER_PATH,
      plotting = PLOTTING_PATH,
      genotype_cache = GENOTYPE_CACHE_PATH,
      matrix_cache = MATRIX_CACHE_PATH,
      r1_cache = R1_CACHE_DIR,
      runner = RUNNER_PATH
    ),
    sha256 = observed_source_hashes
  ),
  configuration = list(
    J = J,
    seed_list = SEED_LIST,
    rho_grid = RHO_GRID,
    true_pi0 = TRUE_PI0,
    full_output_id = FULL_OUTPUT_ID,
    sweep_output_id = SWEEP_OUTPUT_ID
  ),
  artifact_sha256 = artifact_sha256
)
atomic_save_rds(manifest, file.path(PARTIAL_DIR, "manifest.rds"))
atomic_write_lines(
  c(
    paste0("result_id=", RESULT_ID),
    paste0("completed_at=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("fashr_version=", EXPECTED_FASHR_VERSION),
    paste0("fashr_remote_sha=", EXPECTED_FASHR_REMOTE_SHA),
    "full_matrix_replicates=5",
    "lag1_sweep_replicates=5"
  ),
  file.path(PARTIAL_DIR, "complete.flag")
)
if (!file.rename(PARTIAL_DIR, FINAL_DIR)) {
  stop("Unable to promote the validated R4 cache to ", FINAL_DIR, ".")
}
log_message("R4 fashr 0.1.43 run completed: ", FINAL_DIR)
