#!/usr/bin/env Rscript

# Run and validate the center-aligned, relative-clearance, paired-posterior R3
# refresh. The underlying driver remains resumable at the replicate level;
# this wrapper pins every scientific input and atomically promotes a validated
# partial cache to its final result directory.

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[[1L]] == length(args)) return(default)
  args[[hit[[1L]] + 1L]]
}

as_flag <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

parse_numeric_pair <- function(x, name) {
  values <- suppressWarnings(as.numeric(
    trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  ))
  if (length(values) != 2L ||
      any(!is.finite(values)) ||
      values[[1L]] >= values[[2L]]) {
    stop(name, " must contain two increasing finite values.", call. = FALSE)
  }
  values
}

parse_temporal_category_probs <- function(x, name) {
  values <- suppressWarnings(as.numeric(
    trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  ))
  if (length(values) != 3L ||
      any(!is.finite(values)) ||
      any(values <= 0) ||
      abs(sum(values) - 1) > 1e-8) {
    stop(
      name,
      " must contain positive Early,Middle,Late probabilities that sum to one.",
      call. = FALSE
    )
  }
  stats::setNames(values, c("early", "middle", "late"))
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
  path <- require_file(path, "R3 summary")
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
EXPECTED_R3_DRIVER_SHA256 <-
  "8263ea489afc8b406de071020c180b8b5687de6da81993fe9bf40accca292e35"
EXPECTED_SIMULATION_FUNCTIONS_SHA256 <-
  "94a966ee2d3a71014b0ac52dc3c6aef5c4119f45fe821b889b23d34804874448"
EXPECTED_REAL_GENOTYPE_HELPER_SHA256 <-
  "c03c01a188503336a77793c96f4e2d3ac7e0cbd56f4028b552da4b2f88e6b9d7"
EXPECTED_GENOTYPE_CACHE_SHA256 <-
  "81bbef5f323a0bab2ca993c782d8a9b7c63518b83c2cdb46ef7ed1d46f65af24"
EXPECTED_GENOTYPE_CONTENT_MD5 <- c(
  `12345` = "526a7318aa2af901e09252f5a6ca3c46",
  `22345` = "517faa30d5218a956f1be84f2567369c",
  `32345` = "9b9be3205d7db54dac31763492bcb2eb",
  `42345` = "9ef73aa94a061df868b1a951fe495d9f",
  `52345` = "7dab159b8453e2f66188ae313bfbd611"
)
EXPECTED_RAISED_COSINE_CENTER_RANGES <- list(
  early = c(1.5, 2.5),
  middle = c(4.5, 10.5),
  late = c(12.5, 13.5)
)
EXPECTED_RAISED_COSINE_HALF_WIDTH <- 1.5
EXPECTED_LOCATION_TRUTH_MARGIN <- 0.10
EXPECTED_LOCATION_TRUTH_MIN_RANGE_FRACTION <- 0.10
EXPECTED_FUNCTIONAL_POSTERIOR_PAIRING <- "common_random_seed_raw_bf"

DEFAULT_RESULT_ID <- paste0(
  "r3_real_genotype_one_per_gene_J6362_",
  "matched_functional_open_middle_3_12_center_aligned_equal_cells_",
  "relative_location_clearance_",
  "paired_posterior_fashr0143_pilot5"
)
RESULT_ID <- Sys.getenv("FASH_R3_RESULT_ID", unset = DEFAULT_RESULT_ID)
MIDDLE_WINDOW <- parse_numeric_pair(
  Sys.getenv("FASH_R3_MIDDLE_WINDOW", unset = "3,12"),
  "FASH_R3_MIDDLE_WINDOW"
)
MIDDLE_BOUNDARY <- match.arg(
  Sys.getenv("FASH_R3_MIDDLE_BOUNDARY", unset = "open"),
  c("closed", "open")
)
TEMPORAL_CATEGORY_PROBS <- parse_temporal_category_probs(
  Sys.getenv(
    "FASH_R3_TEMPORAL_CATEGORY_PROBS",
    unset = "0.3333333333333333,0.3333333333333333,0.3333333333333333"
  ),
  "FASH_R3_TEMPORAL_CATEGORY_PROBS"
)
SEED_LIST <- c(12345L, 22345L, 32345L, 42345L, 52345L)
TRUTH_MECHANISMS <- c("random_bspline", "raised_cosine")
J <- 6362L
NUM_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
PREFLIGHT_ONLY <- as_flag(get_arg("--preflight-only", "false"))

if (is.na(NUM_CORES) || NUM_CORES < 1L) {
  stop("SLURM_CPUS_PER_TASK must be a positive integer.", call. = FALSE)
}
if (!nzchar(RESULT_ID) || grepl("[/\\\\]", RESULT_ID)) {
  stop("FASH_R3_RESULT_ID must be a non-empty directory name.", call. = FALSE)
}

WORKSPACE_DIR <- Sys.getenv(
  "FASH_WORKSPACE",
  unset = "/project/mstephens/ziangzhang/fash/workspace"
)
WORKSPACE_DIR <- normalizePath(WORKSPACE_DIR, winslash = "/", mustWork = TRUE)
RESULT_PARENT <- Sys.getenv(
  "FASH_R3_RESULT_PARENT",
  unset = file.path(WORKSPACE_DIR, "results", "revision_simulations")
)
FINAL_DIR <- file.path(RESULT_PARENT, RESULT_ID)
PARTIAL_DIR <- file.path(RESULT_PARENT, paste0(RESULT_ID, "_partial"))
R3_DRIVER_PATH <- require_file(
  Sys.getenv(
    "FASH_R3_DRIVER",
    unset = file.path(
      WORKSPACE_DIR,
      "code", "revision_simulations",
      "r3_center_aligned_relative_location_paired_fashr0143",
      "r3_center_aligned_relative_location_paired_driver.R"
    )
  ),
  "R3 driver"
)
SIMULATION_FUNCTIONS_PATH <- require_file(
  Sys.getenv(
    "FASH_R3_SIMULATION_FUNCTIONS",
    unset = file.path(
      WORKSPACE_DIR,
      "code", "revision_simulations",
      "r3_center_aligned_relative_location_paired_fashr0143",
      "r3_center_aligned_relative_location_paired_simulation_functions.R"
    )
  ),
  "Shared simulation functions"
)
REAL_GENOTYPE_HELPER_PATH <- require_file(
  Sys.getenv(
    "FASH_R3_REAL_GENOTYPE_HELPER",
    unset = file.path(
      WORKSPACE_DIR,
      "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
    )
  ),
  "Real-genotype helper"
)
GENOTYPE_CACHE_PATH <- require_file(
  Sys.getenv(
    "FASH_R3_GENOTYPE_CACHE",
    unset = file.path(
      WORKSPACE_DIR,
      "inputs", "r1_r2_fashr0143", "genotype_samples.rds"
    )
  ),
  "Formal genotype cache"
)
RUNNER_PATH <- normalizePath(
  sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][[1L]]),
  winslash = "/",
  mustWork = TRUE
)
WRAPPER_CORE_PATH <- require_file(
  Sys.getenv(
    "FASH_R3_WRAPPER_CORE",
    unset = file.path(
      WORKSPACE_DIR,
      "scripts", "21_r3_center_aligned_relative_location_paired_core.R"
    )
  ),
  "R3 wrapper core"
)

observed_source_hashes <- list(
  r3_driver = sha256_file(R3_DRIVER_PATH),
  simulation_functions = sha256_file(SIMULATION_FUNCTIONS_PATH),
  real_genotype_helper = sha256_file(REAL_GENOTYPE_HELPER_PATH),
  genotype_cache = sha256_file(GENOTYPE_CACHE_PATH),
  wrapper_core = sha256_file(WRAPPER_CORE_PATH),
  runner = sha256_file(RUNNER_PATH)
)
expected_source_hashes <- list(
  r3_driver = EXPECTED_R3_DRIVER_SHA256,
  simulation_functions = EXPECTED_SIMULATION_FUNCTIONS_SHA256,
  real_genotype_helper = EXPECTED_REAL_GENOTYPE_HELPER_SHA256,
  genotype_cache = EXPECTED_GENOTYPE_CACHE_SHA256
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
TRUTH_GROUP_LEVELS <- c(
  "early / switch",
  "early / non-switch",
  "middle / switch",
  "middle / non-switch",
  "late / switch",
  "late / non-switch"
)
EXPECTED_TRUTH_GROUP_COUNTS <- stats::setNames(
  as.integer(exact_temporal_truth_group_counts(
    1272L,
    temporal_category_probs = TEMPORAL_CATEGORY_PROBS
  )[TRUTH_GROUP_LEVELS]),
  TRUTH_GROUP_LEVELS
)
evaluation_grid <- seq(0, 15, by = 0.1)
expected_middle_membership <- temporal_middle_membership(
  smooth_var = evaluation_grid,
  middle_window = MIDDLE_WINDOW,
  middle_boundary = MIDDLE_BOUNDARY
)
EXPECTED_MIDDLE_GRID <- evaluation_grid[expected_middle_membership]
EXPECTED_MIDDLE_EXPRESSION <- if (identical(MIDDLE_BOUNDARY, "open")) {
  sprintf("%g < t < %g", MIDDLE_WINDOW[[1L]], MIDDLE_WINDOW[[2L]])
} else {
  sprintf("%g <= t <= %g", MIDDLE_WINDOW[[1L]], MIDDLE_WINDOW[[2L]])
}
genotype_cache <- readRDS(GENOTYPE_CACHE_PATH)
if (!is.list(genotype_cache) ||
    !all(c("configuration", "sample_ids", "samples") %in% names(genotype_cache)) ||
    !identical(genotype_cache$configuration$n_genes, J) ||
    !identical(genotype_cache$configuration$n_donors, 19L) ||
    !all(SEED_LIST %in% genotype_cache$configuration$seed_list)) {
  stop("The formal genotype cache has an unexpected schema.", call. = FALSE)
}
for (seed in SEED_LIST) {
  sample <- validate_real_genotype_sample(
    genotype_cache$samples[[as.character(seed)]],
    expected_genes = J,
    expected_donors = 19L,
    maf_min = genotype_cache$configuration$maf_min
  )
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

log_message("R3 result: ", FINAL_DIR)
log_message("fashr ", package_provenance$version, " at ", package_provenance$remote_sha)
log_message("Middle functional: ", EXPECTED_MIDDLE_EXPRESSION)
log_message(
  "Temporal truth-category probabilities: ",
  paste(
    names(TEMPORAL_CATEGORY_PROBS),
    format(TEMPORAL_CATEGORY_PROBS, trim = TRUE),
    sep = "=",
    collapse = "; "
  )
)
log_message(
  paste(
    "Raised-cosine centers: early 1.5-2.5; middle 4.5-10.5; late 12.5-13.5.",
    "Location truth clearance is scale adaptive for all three time regions."
  )
)
log_message(
  "Location clearance: max(", EXPECTED_LOCATION_TRUTH_MARGIN, ", ",
  EXPECTED_LOCATION_TRUTH_MIN_RANGE_FRACTION, " x effect range)."
)
log_message("Functional posterior pairing: ", EXPECTED_FUNCTIONAL_POSTERIOR_PAIRING)
log_message("Validated R3 source files and the five-seed genotype cache.")

old_working_directory <- setwd(WORKSPACE_DIR)
on.exit(setwd(old_working_directory), add = TRUE)
arguments <- c(
  "--vanilla",
  R3_DRIVER_PATH,
  "--J", as.character(J),
  "--n-donors", "19",
  "--n-covariates", "5",
  "--noise-sd", "1",
  "--dynamic-main-effect-sd", "1",
  "--truth-mechanisms", paste(TRUTH_MECHANISMS, collapse = ","),
  "--seed-list", paste(SEED_LIST, collapse = ","),
  "--genotype-cache", GENOTYPE_CACHE_PATH,
  "--output-id", RESULT_ID,
  "--output-dir", PARTIAL_DIR,
  "--num-cores", as.character(NUM_CORES),
  "--expected-fashr-version", EXPECTED_FASHR_VERSION,
  "--expected-fashr-remote-sha", EXPECTED_FASHR_REMOTE_SHA,
  "--expected-genotype-content-md5",
  paste(unname(EXPECTED_GENOTYPE_CONTENT_MD5), collapse = ","),
  "--middle-window", paste(MIDDLE_WINDOW, collapse = ","),
  "--middle-boundary", MIDDLE_BOUNDARY,
  "--temporal-category-probs",
  paste(sprintf("%.17g", TEMPORAL_CATEGORY_PROBS), collapse = ","),
  "--location-truth-min-range-fraction",
  as.character(EXPECTED_LOCATION_TRUTH_MIN_RANGE_FRACTION)
)

if (PREFLIGHT_ONLY) {
  log_message("Running the R3 driver input preflight.")
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(arguments, "--preflight-only", "true")
  )
  if (!identical(status, 0L)) {
    stop("The R3 driver preflight failed with status ", status, ".", call. = FALSE)
  }
  log_message("R3 preflight completed successfully.")
  quit(save = "no", status = 0L)
}
if (dir.exists(FINAL_DIR) || file.exists(FINAL_DIR)) {
  stop("Refusing to overwrite completed R3 output: ", FINAL_DIR, call. = FALSE)
}
ensure_directory(RESULT_PARENT)

log_message("Starting or resuming the ten formal R3 replicates.")
status <- system2(file.path(R.home("bin"), "Rscript"), arguments)
if (!identical(status, 0L)) {
  stop("The R3 driver failed with status ", status, ".", call. = FALSE)
}

configuration <- readRDS(require_file(
  file.path(PARTIAL_DIR, "configuration.rds"),
  "R3 configuration"
))
temporal_category_probs_match <- isTRUE(all.equal(
  configuration$temporal_category_probs,
  TEMPORAL_CATEGORY_PROBS,
  tolerance = 1e-12,
  check.attributes = TRUE
))
if (!identical(configuration$output_id, RESULT_ID) ||
    !identical(configuration$J, J) ||
    !identical(configuration$seed_list, SEED_LIST) ||
    !identical(configuration$truth_mechanisms, TRUTH_MECHANISMS) ||
    !identical(configuration$middle_window, MIDDLE_WINDOW) ||
    !identical(configuration$middle_boundary, MIDDLE_BOUNDARY) ||
    !identical(
      configuration$location_truth_margin,
      EXPECTED_LOCATION_TRUTH_MARGIN
    ) ||
    !identical(
      configuration$location_truth_min_range_fraction,
      EXPECTED_LOCATION_TRUTH_MIN_RANGE_FRACTION
    ) ||
    !identical(
      configuration$functional_posterior_pairing,
      EXPECTED_FUNCTIONAL_POSTERIOR_PAIRING
    ) ||
    !identical(
      configuration$temporal_category_design,
      "equal temporal categories"
    ) ||
    !temporal_category_probs_match ||
    !identical(
      configuration$expected_truth_group_counts,
      EXPECTED_TRUTH_GROUP_COUNTS
    ) ||
    !identical(
      configuration$raised_cosine$center_ranges,
      EXPECTED_RAISED_COSINE_CENTER_RANGES
    ) ||
    !identical(
      configuration$raised_cosine$width_half,
      EXPECTED_RAISED_COSINE_HALF_WIDTH
    ) ||
    !isTRUE(all.equal(
      configuration$middle_grid,
      EXPECTED_MIDDLE_GRID,
      tolerance = 1e-12,
      check.attributes = FALSE
    )) ||
    !identical(configuration$middle_expression, EXPECTED_MIDDLE_EXPRESSION) ||
    !identical(
      configuration$genotype_digest_method,
      "fash-genotype-content-md5-v1"
    ) ||
    !identical(
      configuration$genotype_content_digests,
      EXPECTED_GENOTYPE_CONTENT_MD5
    ) ||
    !identical(configuration$package_provenance$version, EXPECTED_FASHR_VERSION) ||
    !identical(
      configuration$package_provenance$remote_sha,
      EXPECTED_FASHR_REMOTE_SHA
    )) {
  stop("The completed R3 configuration is invalid.", call. = FALSE)
}

replicate_paths <- unlist(lapply(TRUTH_MECHANISMS, function(mechanism) {
  file.path(
    PARTIAL_DIR,
    "replicates",
    paste0(mechanism, "_seed_", SEED_LIST, ".rds")
  )
}), use.names = FALSE)
for (path in replicate_paths) require_file(path, "R3 replicate")
replicates <- lapply(replicate_paths, readRDS)
expected_keys <- data.frame(
  truth_mechanism = rep(TRUTH_MECHANISMS, each = length(SEED_LIST)),
  seed = rep(SEED_LIST, times = length(TRUTH_MECHANISMS)),
  stringsAsFactors = FALSE
)
observed_keys <- data.frame(
  truth_mechanism = vapply(replicates, `[[`, character(1), "truth_mechanism"),
  seed = vapply(replicates, `[[`, integer(1), "seed"),
  stringsAsFactors = FALSE
)
if (!identical(observed_keys, expected_keys)) {
  stop("The R3 replicate keys are incomplete or misordered.", call. = FALSE)
}
for (index in seq_along(replicates)) {
  replicate <- replicates[[index]]
  seed <- replicate$seed
  if (!isTRUE(all.equal(replicate$configuration, configuration)) ||
      !identical(
        replicate$genotype_digest,
        EXPECTED_GENOTYPE_CONTENT_MD5[[as.character(seed)]]
      ) ||
      !identical(
        replicate$selected_pair_keys,
        genotype_cache$samples[[as.character(seed)]]$selection$pair_key
      ) ||
      nrow(replicate$functional_alpha) != 320L ||
      nrow(replicate$functional_alpha_005) != 8L ||
      any(!is.finite(replicate$functional_alpha$power)) ||
      any(!is.finite(replicate$functional_alpha$empirical_fsr))) {
    stop("An R3 replicate failed the final contract: ", replicate_paths[[index]], call. = FALSE)
  }
}

summary_dir <- file.path(PARTIAL_DIR, "summary")
all_alpha <- read_required_csv(
  file.path(summary_dir, "all_replicate_functional_alpha_curves.csv"), 3200L
)
all_alpha_005 <- read_required_csv(
  file.path(summary_dir, "all_replicate_functional_alpha005.csv"), 80L
)
all_truth_groups <- read_required_csv(
  file.path(summary_dir, "all_truth_group_counts.csv"), 60L
)
all_pi0 <- read_required_csv(
  file.path(summary_dir, "all_replicate_pi0.csv"), 20L
)
read_required_csv(file.path(summary_dir, "genotype_selection_summary.csv"), 10L)
read_required_csv(file.path(summary_dir, "truth_maf_balance.csv"), 30L)
mc_alpha <- read_required_csv(
  file.path(summary_dir, "functional_testing_mc_alpha_curve.csv"), 640L
)
read_required_csv(
  file.path(summary_dir, "functional_testing_mc_alpha005_summary.csv"), 16L
)
read_required_csv(
  file.path(summary_dir, "functional_testing_mc_pi0_summary.csv"), 4L
)
require_file(file.path(PARTIAL_DIR, "example_curves.rds"), "R3 example curves")
if (any(!is.finite(all_alpha$power)) ||
    any(!is.finite(all_alpha$empirical_fsr)) ||
    any(!is.finite(all_alpha_005$power)) ||
    any(!is.finite(all_pi0$estimated_pi0)) ||
    any(!is.finite(mc_alpha$mean_power)) ||
    nrow(all_truth_groups) != 60L) {
  stop("The R3 summaries contain invalid numerical values.", call. = FALSE)
}
for (mechanism in TRUTH_MECHANISMS) {
  for (seed in SEED_LIST) {
    rows <- all_truth_groups[
      all_truth_groups$truth_mechanism == mechanism &
        all_truth_groups$seed == seed,
      ,
      drop = FALSE
    ]
    observed_counts <- stats::setNames(
      as.integer(rows$n_dynamic),
      rows$truth_group
    )[TRUTH_GROUP_LEVELS]
    if (!identical(observed_counts, EXPECTED_TRUTH_GROUP_COUNTS)) {
      stop(
        "Unexpected truth-group counts for ", mechanism,
        " seed ", seed, ".",
        call. = FALSE
      )
    }
  }
}

middle_gate_rows <- mc_alpha[
  mc_alpha$target == "middle" &
    mc_alpha$method == "FASH-IWP1-BF" &
    mc_alpha$alpha >= 0.05,
  ,
  drop = FALSE
]
middle_gate_rows$truth_mechanism <- ifelse(
  grepl("^r3a_", middle_gate_rows$scenario),
  "random_bspline",
  ifelse(
    grepl("^r3b_", middle_gate_rows$scenario),
    "raised_cosine",
    NA_character_
  )
)
if (nrow(middle_gate_rows) != 62L ||
    anyNA(middle_gate_rows$truth_mechanism)) {
  stop("Could not resolve the formal Middle calibration-gate rows.", call. = FALSE)
}
scientific_validation <- do.call(rbind, lapply(
  split(middle_gate_rows, middle_gate_rows$truth_mechanism),
  function(rows) {
    excess <- rows$mean_empirical_fsr - rows$alpha
    maximum_index <- which.max(excess)
    data.frame(
      truth_mechanism = rows$truth_mechanism[[1L]],
      method = "FASH-IWP1-BF",
      target = "middle",
      alpha_min = min(rows$alpha),
      alpha_max = max(rows$alpha),
      maximum_excess = excess[[maximum_index]],
      alpha_at_maximum_excess = rows$alpha[[maximum_index]],
      mean_empirical_fsr_at_maximum =
        rows$mean_empirical_fsr[[maximum_index]],
      prespecified_maximum_excess = 0.03,
      passed = excess[[maximum_index]] <= 0.03,
      stringsAsFactors = FALSE
    )
  }
))
utils::write.csv(
  scientific_validation,
  file.path(PARTIAL_DIR, "scientific_validation.csv"),
  row.names = FALSE
)
CALIBRATION_GATE_PASSED <- all(scientific_validation$passed)
log_message(
  "Prespecified Middle calibration gate passed: ",
  CALIBRATION_GATE_PASSED
)

artifact_paths <- list.files(
  PARTIAL_DIR,
  recursive = TRUE,
  full.names = TRUE,
  all.files = FALSE
)
artifact_paths <- artifact_paths[file.info(artifact_paths)$isdir %in% FALSE]
artifact_relative_paths <- substring(
  artifact_paths,
  nchar(normalizePath(PARTIAL_DIR, winslash = "/", mustWork = TRUE)) + 2L
)
artifact_sha256 <- stats::setNames(
  vapply(artifact_paths, sha256_file, character(1)),
  artifact_relative_paths
)
manifest <- list(
  schema_version = paste0(
    "r3-fashr0143-manifest-v7-center-aligned-",
    "relative-location-paired"
  ),
  result_id = RESULT_ID,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  package_provenance = package_provenance,
  source_provenance = list(
    paths = list(
      r3_driver = R3_DRIVER_PATH,
      simulation_functions = SIMULATION_FUNCTIONS_PATH,
      real_genotype_helper = REAL_GENOTYPE_HELPER_PATH,
      genotype_cache = GENOTYPE_CACHE_PATH,
      wrapper_core = WRAPPER_CORE_PATH,
      runner = RUNNER_PATH
    ),
    sha256 = observed_source_hashes
  ),
  configuration = configuration,
  artifact_sha256 = artifact_sha256
)
atomic_save_rds(manifest, file.path(PARTIAL_DIR, "manifest.rds"))
atomic_write_lines(
  c(
    paste0("result_id=", RESULT_ID),
    paste0("completed_at=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("fashr_version=", EXPECTED_FASHR_VERSION),
    paste0("fashr_remote_sha=", EXPECTED_FASHR_REMOTE_SHA),
    paste0("middle_definition=", EXPECTED_MIDDLE_EXPRESSION),
    paste0(
      "temporal_category_probs=",
      paste(
        names(TEMPORAL_CATEGORY_PROBS),
        TEMPORAL_CATEGORY_PROBS,
        sep = ":",
        collapse = ";"
      )
    ),
    paste0(
      "raised_cosine_center_ranges=",
      "early:1.5,2.5;middle:4.5,10.5;late:12.5,13.5"
    ),
    paste0(
      "raised_cosine_half_width=",
      EXPECTED_RAISED_COSINE_HALF_WIDTH
    ),
    paste0("location_truth_margin=", EXPECTED_LOCATION_TRUTH_MARGIN),
    paste0(
      "location_truth_min_range_fraction=",
      EXPECTED_LOCATION_TRUTH_MIN_RANGE_FRACTION
    ),
    paste0(
      "functional_posterior_pairing=",
      EXPECTED_FUNCTIONAL_POSTERIOR_PAIRING
    ),
    paste0("calibration_gate_passed=", CALIBRATION_GATE_PASSED),
    "replicates=10"
  ),
  file.path(PARTIAL_DIR, "complete.flag")
)
if (!file.rename(PARTIAL_DIR, FINAL_DIR)) {
  stop("Unable to promote the validated R3 cache to ", FINAL_DIR, ".", call. = FALSE)
}
log_message("R3 fashr 0.1.43 run completed: ", FINAL_DIR)
