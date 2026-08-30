#!/usr/bin/env Rscript

# Run the resumable R3 ideal-Gaussian measurement isolation experiment under
# the canonical IWP1 temporal mixture. The frozen R3 truth and t-adjusted
# standard errors are regenerated; only beta_hat uses the ideal observation.

find_workflowr_root <- function() {
  candidates <- c(".", "coderepo-local", "../..", "../../..", "../../../..")
  marker <- file.path(
    candidates,
    "code", "revision_simulations", "shared", "simulation_functions.R"
  )
  hits <- candidates[file.exists(marker)]
  if (length(hits) == 0L) {
    stop("Could not locate the workflowr project root.", call. = FALSE)
  }
  normalizePath(hits[[1L]], winslash = "/", mustWork = TRUE)
}

get_arg <- function(name, default = NULL) {
  arguments <- commandArgs(trailingOnly = TRUE)
  index <- which(arguments == name)
  if (length(index) == 0L || index[[1L]] == length(arguments)) return(default)
  arguments[[index[[1L]] + 1L]]
}

as_flag <- function(value) {
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

parse_integer_list <- function(value, name) {
  output <- suppressWarnings(as.integer(
    trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  ))
  if (length(output) == 0L || anyNA(output) || anyDuplicated(output)) {
    stop(name, " must contain unique comma-separated integers.", call. = FALSE)
  }
  output
}

parse_character_list <- function(value, name) {
  output <- trimws(strsplit(value, ",", fixed = TRUE)[[1L]])
  if (length(output) == 0L || any(!nzchar(output)) || anyDuplicated(output)) {
    stop(name, " must contain unique comma-separated values.", call. = FALSE)
  }
  output
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
  arguments <- if (identical(command, "shasum")) c("-a", "256", path) else path
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

atomic_write_csv <- function(object, path) {
  temporary_path <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)
  utils::write.csv(object, temporary_path, row.names = FALSE)
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

runner_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (length(runner_argument) != 1L) {
  stop("Could not resolve the runner path.", call. = FALSE)
}
RUNNER_PATH <- normalizePath(
  sub("^--file=", "", runner_argument[[1L]]),
  winslash = "/",
  mustWork = TRUE
)
PROJECT_ROOT <- find_workflowr_root()

EXPECTED_FASHR_VERSION <- "0.1.43"
EXPECTED_FASHR_REMOTE_SHA <-
  "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
EXPECTED_HELPER_SHA256 <-
  "84bc06be91531f1587b0f298bea82f0bd937444663419c0b77d589b48d3fe84e"
EXPECTED_SIMULATION_FUNCTIONS_SHA256 <-
  "45267b0884168e5ae33cc4f14e3f05b711d961b65bf9c6fbd880e748de064a6e"
EXPECTED_REAL_GENOTYPE_HELPER_SHA256 <-
  "c03c01a188503336a77793c96f4e2d3ac7e0cbd56f4028b552da4b2f88e6b9d7"
EXPECTED_GENOTYPE_CACHE_SHA256 <-
  "81bbef5f323a0bab2ca993c782d8a9b7c63518b83c2cdb46ef7ed1d46f65af24"
EXPECTED_TEMPORAL_MIXTURE_CONTRACT_SHA256 <-
  "8e2b4c527b6c17fd2645b9d46e47a7b99410d2e9be1db1a38d541380d91ad723"
EXPECTED_GENOTYPE_CONTENT_MD5 <- c(
  `12345` = "526a7318aa2af901e09252f5a6ca3c46",
  `22345` = "517faa30d5218a956f1be84f2567369c",
  `32345` = "9b9be3205d7db54dac31763492bcb2eb",
  `42345` = "9ef73aa94a061df868b1a951fe495d9f",
  `52345` = "7dab159b8453e2f66188ae313bfbd611"
)
RESULT_ID <- paste0(
  "r3_ideal_gaussian_known_t_adjusted_se_",
  paste0(
    "matched_truth_open_middle_3_12_center_aligned_",
    "iwp1_geometry_mixture_full_universe_"
  ),
  "fashr0143_pilot5"
)
BASE_R3_RESULT_ID <- paste0(
  "r3_real_genotype_one_per_gene_J6362_",
  paste0(
    "matched_functional_open_middle_3_12_center_aligned_",
    "iwp1_geometry_mixture_"
  ),
  "relative_location_clearance_full_universe_",
  "paired_posterior_fashr0143_pilot5"
)
SCHEMA_VERSION <- "r3-ideal-gaussian-measurement-iwp1-mixture-v2"
SEED_LIST <- c(12345L, 22345L, 32345L, 42345L, 52345L)
TRUTH_MECHANISMS <- c("random_bspline", "raised_cosine")
RUN_SEED_LIST <- parse_integer_list(
  get_arg("--run-seed-list", paste(SEED_LIST, collapse = ",")),
  "--run-seed-list"
)
RUN_TRUTH_MECHANISMS <- parse_character_list(
  get_arg("--run-truth-mechanisms", paste(TRUTH_MECHANISMS, collapse = ",")),
  "--run-truth-mechanisms"
)
PREFLIGHT_ONLY <- as_flag(get_arg("--preflight-only", "false"))
NUM_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
if (is.na(NUM_CORES) || NUM_CORES < 1L) {
  stop("SLURM_CPUS_PER_TASK must be a positive integer.", call. = FALSE)
}
if (!all(RUN_SEED_LIST %in% SEED_LIST)) {
  stop("--run-seed-list contains a seed outside the frozen R3 design.", call. = FALSE)
}
if (!all(RUN_TRUTH_MECHANISMS %in% TRUTH_MECHANISMS)) {
  stop(
    "--run-truth-mechanisms contains a mechanism outside the frozen R3 design.",
    call. = FALSE
  )
}

HELPER_PATH <- require_file(
  Sys.getenv(
    "FASH_R3_IDEAL_HELPER",
    unset = file.path(dirname(RUNNER_PATH), "ideal_gaussian_measurement.R")
  ),
  "Ideal-Gaussian helper"
)
SIMULATION_FUNCTIONS_PATH <- require_file(
  Sys.getenv(
    "FASH_R3_SIMULATION_FUNCTIONS",
    unset = file.path(
      PROJECT_ROOT,
      "code", "revision_simulations", "r3_r4_fashr0143", "source_snapshots",
      "r3_full_universe_functional_simulation_functions.R"
    )
  ),
  "Frozen R3 simulation functions"
)
REAL_GENOTYPE_HELPER_PATH <- require_file(
  Sys.getenv(
    "FASH_R3_REAL_GENOTYPE_HELPER",
    unset = file.path(
      PROJECT_ROOT,
      "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
    )
  ),
  "Real-genotype helper"
)
GENOTYPE_CACHE_PATH <- require_file(
  Sys.getenv(
    "FASH_R3_GENOTYPE_CACHE",
    unset = file.path(
      PROJECT_ROOT,
      "output", "revision_simulations", "shared",
      "real_genotype_one_per_gene_J6362_pilot5", "genotype_samples.rds"
    )
  ),
  "Formal R3 genotype cache"
)
TEMPORAL_MIXTURE_CONTRACT_PATH <- require_file(
  Sys.getenv(
    "FASH_R3_TEMPORAL_MIXTURE_CONTRACT",
    unset = file.path(
      PROJECT_ROOT,
      "code", "revision_simulations", "shared",
      "r3_iwp1_temporal_mixture_contract.R"
    )
  ),
  "R3 temporal-mixture contract"
)
RESULT_PARENT <- Sys.getenv(
  "FASH_R3_IDEAL_RESULT_PARENT",
  unset = file.path(
    PROJECT_ROOT, "output", "revision_simulations", "internal"
  )
)
FINAL_DIR <- file.path(RESULT_PARENT, RESULT_ID)
PARTIAL_DIR <- file.path(RESULT_PARENT, paste0(RESULT_ID, "_partial"))

observed_source_hashes <- c(
  helper = sha256_file(HELPER_PATH),
  simulation_functions = sha256_file(SIMULATION_FUNCTIONS_PATH),
  real_genotype_helper = sha256_file(REAL_GENOTYPE_HELPER_PATH),
  genotype_cache = sha256_file(GENOTYPE_CACHE_PATH),
  temporal_mixture_contract = sha256_file(TEMPORAL_MIXTURE_CONTRACT_PATH),
  runner = sha256_file(RUNNER_PATH)
)
expected_source_hashes <- c(
  helper = EXPECTED_HELPER_SHA256,
  simulation_functions = EXPECTED_SIMULATION_FUNCTIONS_SHA256,
  real_genotype_helper = EXPECTED_REAL_GENOTYPE_HELPER_SHA256,
  genotype_cache = EXPECTED_GENOTYPE_CACHE_SHA256,
  temporal_mixture_contract = EXPECTED_TEMPORAL_MIXTURE_CONTRACT_SHA256
)
if (!identical(
  observed_source_hashes[names(expected_source_hashes)],
  expected_source_hashes
)) {
  stop(
    "One or more frozen R3 inputs have an unexpected SHA-256 digest.",
    call. = FALSE
  )
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
source(HELPER_PATH)
source(TEMPORAL_MIXTURE_CONTRACT_PATH)

J <- 6362L
N_DONORS <- 19L
N_COVARIATES <- 5L
TIME_GRID <- 0:15
EVALUATION_GRID <- seq(0, 15, by = 0.1)
MIDDLE_WINDOW <- c(3, 12)
MIDDLE_BOUNDARY <- "open"
MIDDLE_EXPRESSION <- "3 < t < 12"
ALPHA_GRID <- seq(0.005, 0.20, by = 0.005)
CLASS_PROBS <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
EXPECTED_CLASS_COUNTS <- c(
  dynamic_bspline = 1272L,
  constant = 2545L,
  zero = 2545L
)
TEMPORAL_MIXTURE_CONTRACT <- r3_iwp1_temporal_mixture_contract(
  EXPECTED_CLASS_COUNTS[["dynamic_bspline"]]
)
TEMPORAL_CATEGORY_PROBS <-
  TEMPORAL_MIXTURE_CONTRACT$temporal_category_probs
TRUTH_GROUP_LEVELS <- c(
  "early / switch",
  "early / non-switch",
  "middle / switch",
  "middle / non-switch",
  "late / switch",
  "late / non-switch"
)
EXPECTED_TRUTH_GROUP_COUNTS <-
  TEMPORAL_MIXTURE_CONTRACT$truth_group_counts[TRUTH_GROUP_LEVELS]
TARGETS <- c("early", "middle", "late", "switch")
METHODS <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
RAISED_COSINE_CENTER_RANGES <- list(
  early = c(1.5, 2.5),
  middle = c(4.5, 10.5),
  late = c(12.5, 13.5)
)

genotype_cache <- readRDS(GENOTYPE_CACHE_PATH)
if (!is.list(genotype_cache) ||
    !all(c("configuration", "sample_ids", "samples") %in% names(genotype_cache)) ||
    !identical(genotype_cache$configuration$n_genes, J) ||
    !identical(genotype_cache$configuration$n_donors, N_DONORS) ||
    !all(SEED_LIST %in% genotype_cache$configuration$seed_list) ||
    !all(as.character(SEED_LIST) %in% names(genotype_cache$samples))) {
  stop("The formal genotype cache does not match the frozen R3 design.", call. = FALSE)
}
genotype_content_digests <- stats::setNames(
  character(length(SEED_LIST)), as.character(SEED_LIST)
)
for (seed in SEED_LIST) {
  sample <- validate_real_genotype_sample(
    genotype_cache$samples[[as.character(seed)]],
    expected_genes = J,
    expected_donors = N_DONORS,
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
  existing_component_seeds <- unname(revision_component_seeds(seed))
  if (ideal_gaussian_measurement_seed(seed) %in% existing_component_seeds) {
    stop("The ideal measurement seed collides with an R3 component seed.", call. = FALSE)
  }
  genotype_content_digests[[as.character(seed)]] <- observed_digest
}

configuration <- list(
  schema_version = SCHEMA_VERSION,
  output_id = RESULT_ID,
  base_r3_result_id = BASE_R3_RESULT_ID,
  J = J,
  n_donors = N_DONORS,
  n_covariates = N_COVARIATES,
  time_grid = TIME_GRID,
  evaluation_grid = EVALUATION_GRID,
  middle_window = MIDDLE_WINDOW,
  middle_boundary = MIDDLE_BOUNDARY,
  middle_expression = MIDDLE_EXPRESSION,
  alpha_grid = ALPHA_GRID,
  class_probs = CLASS_PROBS,
  expected_class_counts = EXPECTED_CLASS_COUNTS,
  temporal_category_probs = TEMPORAL_CATEGORY_PROBS,
  temporal_category_design = TEMPORAL_MIXTURE_CONTRACT$name,
  temporal_category_counts =
    TEMPORAL_MIXTURE_CONTRACT$temporal_category_counts,
  switch_status_counts = TEMPORAL_MIXTURE_CONTRACT$switch_status_counts,
  temporal_mixture_contract_sha256 =
    observed_source_hashes[["temporal_mixture_contract"]],
  expected_truth_group_counts = EXPECTED_TRUTH_GROUP_COUNTS,
  truth_mechanisms = TRUTH_MECHANISMS,
  random_bspline = list(amplitude = 2, df = 6, coefficient_sd = 1),
  raised_cosine = list(
    width_half = 1.5,
    spike_counts = 1:3,
    relative_amplitude_range = c(0.35, 0.75),
    target_centered_rms = 0.90,
    center_ranges = RAISED_COSINE_CENTER_RANGES
  ),
  dynamic_main_effect_sd = 1,
  expression_noise_sd = 1,
  covariate_effect_sd = 0.5,
  intercept_sd = 0,
  switch_threshold = 0.25,
  location_truth_margin = 0.10,
  location_truth_min_range_fraction = 0.10,
  switch_truth_margin = 0.10,
  non_switch_min_abs = 0.10,
  non_switch_min_range_fraction = 0.10,
  observation_model = "ideal_gaussian_known_t_adjusted_se",
  observation_equation = "beta_hat = true_beta + adjusted_se * iid N(0,1)",
  standard_error_source = "R3 final t-adjusted SE regenerated from frozen seeds",
  ideal_measurement_seed_offset = 23003L,
  functional_candidate_scope = "full_universe",
  functional_candidate_universe_size = J,
  functional_posterior_pairing = "common_random_seed_raw_bf",
  fash_order = 1L,
  num_basis = 20L,
  penalty = 10,
  pred_step = 1,
  genotype_source = "paper-derived YRI DS dosage",
  genotype_selection_rule = genotype_cache$configuration$selection_rule,
  genotype_dosage_field = genotype_cache$configuration$dosage_field,
  genotype_maf_min = genotype_cache$configuration$maf_min,
  genotype_sample_ids = genotype_cache$sample_ids,
  genotype_digest_method = "fash-genotype-content-md5-v1",
  genotype_content_digests = genotype_content_digests,
  seed_list = SEED_LIST,
  package_provenance = package_provenance,
  source_sha256 = observed_source_hashes
)

if (PREFLIGHT_ONLY) {
  log_message(
    "Validated the ideal-Gaussian package, frozen sources, genotype cache, ",
    "seed separation, and canonical IWP1 temporal-mixture configuration."
  )
  log_message(
    "Temporal probabilities: ",
    paste(
      names(TEMPORAL_CATEGORY_PROBS), TEMPORAL_CATEGORY_PROBS,
      sep = "=", collapse = "; "
    ),
    ". Truth-cell counts: ",
    paste(EXPECTED_TRUTH_GROUP_COUNTS, collapse = ","),
    "."
  )
  quit(save = "no", status = 0L)
}

if (dir.exists(FINAL_DIR) || file.exists(FINAL_DIR)) {
  stop("Refusing to overwrite completed output: ", FINAL_DIR, call. = FALSE)
}
ensure_directory(RESULT_PARENT)
ensure_directory(PARTIAL_DIR)
REPLICATE_DIR <- file.path(PARTIAL_DIR, "replicates")
SUMMARY_DIR <- file.path(PARTIAL_DIR, "summary")
ensure_directory(REPLICATE_DIR)
ensure_directory(SUMMARY_DIR)
configuration_path <- file.path(PARTIAL_DIR, "configuration.rds")
if (file.exists(configuration_path)) {
  cached_configuration <- readRDS(configuration_path)
  if (!isTRUE(all.equal(cached_configuration, configuration))) {
    stop("The partial output has a different configuration.", call. = FALSE)
  }
} else {
  atomic_save_rds(configuration, configuration_path)
}

formal_truth_scenario <- function(truth_mechanism) {
  prefix <- if (identical(truth_mechanism, "random_bspline")) "r3a_" else "r3b_"
  paste0(
    prefix,
    "real_genotype_one_per_gene_matched_functional_",
    truth_mechanism,
    paste0(
      "_open_middle_3_12_center_aligned_",
      "iwp1_geometry_mixture_"
    ),
    "relative_location_clearance_full_universe_",
    "paired_posterior_main_effect"
  )
}

ideal_scenario <- function(truth_mechanism) {
  paste0(
    formal_truth_scenario(truth_mechanism),
    "_ideal_gaussian_known_t_adjusted_se"
  )
}

extract_example_observations <- function(effect_sim,
                                         regression_summary,
                                         ideal_summary,
                                         seed,
                                         truth_mechanism) {
  dynamic <- effect_sim$unit_info$effect_class == "dynamic_bspline"
  rows <- lapply(TRUTH_GROUP_LEVELS, function(truth_group) {
    candidates <- which(dynamic & effect_sim$unit_info$truth_group == truth_group)
    if (length(candidates) == 0L) {
      stop("No example is available for truth group ", truth_group, ".", call. = FALSE)
    }
    index <- candidates[[1L]]
    data.frame(
      seed = seed,
      truth_mechanism = truth_mechanism,
      truth_group = truth_group,
      unit_index = index,
      variant_id = effect_sim$unit_info$variant_id[[index]],
      time = TIME_GRID,
      true_beta = effect_sim$beta_matrix[index, ],
      regression_beta_hat = regression_summary$beta_hat[index, ],
      adjusted_se = regression_summary$se[index, ],
      ideal_beta_hat = ideal_summary$beta_hat[index, ],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_replicate <- function(seed, truth_mechanism) {
  start_time <- proc.time()[["elapsed"]]
  component_seeds <- revision_component_seeds(seed)
  measurement_seed <- ideal_gaussian_measurement_seed(seed)
  genotype_sample <- validate_real_genotype_sample(
    genotype_cache$samples[[as.character(seed)]],
    expected_genes = J,
    expected_donors = N_DONORS,
    maf_min = genotype_cache$configuration$maf_min
  )
  covariates <- simulate_covariate_matrix(
    n_donors = N_DONORS,
    n_covariates = N_COVARIATES,
    seed = component_seeds[["covariates"]]
  )
  truth_scenario <- formal_truth_scenario(truth_mechanism)
  effect_sim <- simulate_matched_functional_effect_set(
    n_variants = J,
    truth_mechanism = truth_mechanism,
    time_grid = TIME_GRID,
    evaluation_grid = EVALUATION_GRID,
    class_probs = CLASS_PROBS,
    dynamic_main_effect_sd = 1,
    cosine_center_ranges = RAISED_COSINE_CENTER_RANGES,
    switch_threshold = 0.25,
    location_truth_margin = 0.10,
    location_truth_min_range_fraction = 0.10,
    switch_truth_margin = 0.10,
    non_switch_min_abs = 0.10,
    non_switch_min_range_fraction = 0.10,
    temporal_category_probs = TEMPORAL_CATEGORY_PROBS,
    seed = seed,
    class_seed = component_seeds[["classes"]],
    constant_seed = component_seeds[["constant_effects"]],
    shape_seed = component_seeds[["functional_truth"]],
    scenario = truth_scenario,
    middle_window = MIDDLE_WINDOW,
    middle_boundary = MIDDLE_BOUNDARY
  )
  if (!identical(
    effect_sim$settings$raised_cosine_center_ranges,
    RAISED_COSINE_CENTER_RANGES
  )) {
    stop("The generated truth uses unexpected center ranges.", call. = FALSE)
  }
  effect_sim <- reassign_effect_simulation_by_maf(
    effect_sim = effect_sim,
    maf = genotype_sample$variant_info$observed_maf,
    class_probs = CLASS_PROBS,
    seed = component_seeds[["classes"]],
    n_strata = 10L
  )
  for (field in c("beta_matrix", "beta_evaluation", "true_functionals")) {
    rownames(effect_sim[[field]]) <- genotype_sample$selection$pair_key
  }
  effect_sim$unit_info$variant_id <- genotype_sample$selection$pair_key

  expression_sim <- simulate_eqtl_expression_from_genotypes(
    G = genotype_sample$G,
    beta_matrix = effect_sim$beta_matrix,
    time_grid = TIME_GRID,
    covariates = covariates,
    expression_noise_sd = 1,
    covariate_effect_sd = 0.5,
    intercept_sd = 0,
    seed = component_seeds[["expression"]]
  )
  regression_summary <- estimate_eqtl_summaries_from_genotypes(
    G = genotype_sample$G,
    expression = expression_sim$expression,
    covariates = covariates,
    apply_t_se_correction = TRUE
  )
  rm(expression_sim)
  invisible(gc(verbose = FALSE))

  ideal_summary <- simulate_ideal_gaussian_eqtl_summary(
    true_beta = effect_sim$beta_matrix,
    adjusted_se = regression_summary$se,
    seed = measurement_seed,
    regression_summary = regression_summary
  )
  regression_z <-
    (regression_summary$beta_hat - effect_sim$beta_matrix) / regression_summary$se
  ideal_z <- ideal_summary$standardized_error
  observation_moments <- rbind(
    summarize_standardized_error_matrix(
      regression_z, seed, truth_mechanism, "regression_t_adjusted"
    ),
    summarize_standardized_error_matrix(
      ideal_z, seed, truth_mechanism, "ideal_gaussian"
    )
  )
  observation_histogram <- rbind(
    summarize_standardized_error_histogram(
      regression_z, seed, truth_mechanism, "regression_t_adjusted"
    ),
    summarize_standardized_error_histogram(
      ideal_z, seed, truth_mechanism, "ideal_gaussian"
    )
  )
  observation_lag_correlation <- rbind(
    summarize_standardized_error_lag_correlation(
      regression_z, seed, truth_mechanism, "regression_t_adjusted"
    ),
    summarize_standardized_error_lag_correlation(
      ideal_z, seed, truth_mechanism, "ideal_gaussian"
    )
  )

  datasets <- make_fash_datasets_from_eqtl_summary(
    beta_hat = ideal_summary$beta_hat,
    se = ideal_summary$se,
    true_beta = effect_sim$beta_matrix,
    time_grid = TIME_GRID,
    unit_info = effect_sim$unit_info,
    scenario = ideal_scenario(truth_mechanism)
  )
  unit_info <- attr(datasets, "unit_info")
  # Match the RNG state established by the formal R3 wrapper immediately
  # before its FASH fitting stage. The fit is expected to be deterministic,
  # but this prevents an undocumented package RNG dependency from drifting.
  set.seed(seed)
  fash_fits <- fit_fash_for_revision(
    datasets = datasets,
    orders = 1,
    grid = default_revision_grid(),
    num_basis = 20,
    penalty = 10,
    pred_step = 1,
    num_cores = NUM_CORES,
    apply_bf = TRUE,
    verbose = FALSE
  )
  true_functionals <- effect_sim$true_functionals[, TARGETS, drop = FALSE]
  true_dynamic <- unit_info$effect_class == "dynamic_bspline"
  functional_posterior_seed <- component_seeds[["functional_posterior"]]
  raw_results <- evaluate_fash_functional_testing(
    fit = fash_fits$fash_iwp1_raw,
    true_functionals = true_functionals,
    evaluation_grid = EVALUATION_GRID,
    alpha_grid = ALPHA_GRID,
    method = "FASH-IWP1-Raw",
    scenario = ideal_scenario(truth_mechanism),
    switch_threshold = 0.25,
    candidate_scope = "full_universe",
    true_dynamic = true_dynamic,
    num_cores = NUM_CORES,
    seed = functional_posterior_seed,
    middle_window = MIDDLE_WINDOW,
    middle_boundary = MIDDLE_BOUNDARY
  )
  bf_results <- evaluate_fash_functional_testing(
    fit = fash_fits$fash_iwp1_bf,
    true_functionals = true_functionals,
    evaluation_grid = EVALUATION_GRID,
    alpha_grid = ALPHA_GRID,
    method = "FASH-IWP1-BF",
    scenario = ideal_scenario(truth_mechanism),
    switch_threshold = 0.25,
    candidate_scope = "full_universe",
    true_dynamic = true_dynamic,
    num_cores = NUM_CORES,
    seed = functional_posterior_seed,
    middle_window = MIDDLE_WINDOW,
    middle_boundary = MIDDLE_BOUNDARY
  )
  functional_alpha <- rbind(raw_results$alpha_curve, bf_results$alpha_curve)
  functional_alpha$seed <- seed
  functional_alpha$truth_mechanism <- truth_mechanism

  truth_group_counts <- as.data.frame(table(factor(
    unit_info$truth_group[true_dynamic],
    levels = TRUTH_GROUP_LEVELS
  )))
  names(truth_group_counts) <- c("truth_group", "n_dynamic")
  truth_group_counts$seed <- seed
  truth_group_counts$truth_mechanism <- truth_mechanism
  observed_group_counts <- stats::setNames(
    as.integer(truth_group_counts$n_dynamic),
    truth_group_counts$truth_group
  )
  if (!identical(observed_group_counts, EXPECTED_TRUTH_GROUP_COUNTS)) {
    stop("The dynamic truth groups do not match the frozen counts.", call. = FALSE)
  }
  truth_maf_balance <- summarize_truth_maf_balance(
    variant_info = genotype_sample$variant_info,
    unit_info = unit_info,
    seed = seed
  )
  truth_maf_balance$truth_mechanism <- truth_mechanism
  example_observations <- extract_example_observations(
    effect_sim,
    regression_summary,
    ideal_summary,
    seed,
    truth_mechanism
  )
  digests <- data.frame(
    seed = seed,
    truth_mechanism = truth_mechanism,
    genotype_content_md5 = genotype_content_digests[[as.character(seed)]],
    true_beta_md5 = serialized_object_md5(effect_sim$beta_matrix),
    adjusted_se_md5 = serialized_object_md5(regression_summary$se),
    regression_beta_hat_md5 = serialized_object_md5(regression_summary$beta_hat),
    ideal_beta_hat_md5 = serialized_object_md5(ideal_summary$beta_hat),
    stringsAsFactors = FALSE
  )
  elapsed_seconds <- proc.time()[["elapsed"]] - start_time

  list(
    configuration = configuration,
    seed = seed,
    truth_mechanism = truth_mechanism,
    component_seeds = c(
      component_seeds,
      ideal_gaussian_measurement = measurement_seed
    ),
    genotype_digest = genotype_content_digests[[as.character(seed)]],
    selected_pair_keys = genotype_sample$selection$pair_key,
    functional_alpha = functional_alpha,
    functional_alpha_005 = functional_alpha[
      abs(functional_alpha$alpha - 0.05) < 1e-12,
      ,
      drop = FALSE
    ],
    estimated_pi0 = data.frame(
      seed = seed,
      truth_mechanism = truth_mechanism,
      method = METHODS,
      estimated_pi0 = c(
        constant_component_prior_weight(fash_fits$fash_iwp1_raw),
        constant_component_prior_weight(fash_fits$fash_iwp1_bf)
      ),
      stringsAsFactors = FALSE
    ),
    truth_group_counts = truth_group_counts,
    truth_maf_balance = truth_maf_balance,
    observation_moments = observation_moments,
    observation_histogram = observation_histogram,
    observation_lag_correlation = observation_lag_correlation,
    example_observations = example_observations,
    digests = digests,
    elapsed_seconds = elapsed_seconds
  )
}

validate_replicate <- function(object, seed, truth_mechanism) {
  required <- c(
    "configuration", "seed", "truth_mechanism", "component_seeds",
    "genotype_digest", "selected_pair_keys", "functional_alpha",
    "functional_alpha_005", "estimated_pi0", "truth_group_counts",
    "truth_maf_balance", "observation_moments", "observation_histogram",
    "observation_lag_correlation", "example_observations", "digests",
    "elapsed_seconds"
  )
  expected_alpha_rows <- length(METHODS) * length(TARGETS) * length(ALPHA_GRID)
  is.list(object) &&
    all(required %in% names(object)) &&
    identical(object$seed, seed) &&
    identical(object$truth_mechanism, truth_mechanism) &&
    isTRUE(all.equal(object$configuration, configuration)) &&
    identical(
      object$genotype_digest,
      EXPECTED_GENOTYPE_CONTENT_MD5[[as.character(seed)]]
    ) &&
    identical(object$selected_pair_keys, genotype_cache$samples[[
      as.character(seed)
    ]]$selection$pair_key) &&
    nrow(object$functional_alpha) == expected_alpha_rows &&
    nrow(object$functional_alpha_005) == length(METHODS) * length(TARGETS) &&
    all(object$functional_alpha$candidate_scope == "full_universe") &&
    all(object$functional_alpha$candidate_count == J) &&
    all(object$functional_alpha$first_stage_null_calls == 0L) &&
    all(is.finite(object$functional_alpha$power)) &&
    all(is.finite(object$functional_alpha$empirical_fsr)) &&
    nrow(object$estimated_pi0) == length(METHODS) &&
    all(is.finite(object$estimated_pi0$estimated_pi0)) &&
    nrow(object$truth_group_counts) == length(TRUTH_GROUP_LEVELS) &&
    identical(
      stats::setNames(
        as.integer(object$truth_group_counts$n_dynamic),
        as.character(object$truth_group_counts$truth_group)
      )[TRUTH_GROUP_LEVELS],
      EXPECTED_TRUTH_GROUP_COUNTS
    ) &&
    nrow(object$observation_moments) == 2L &&
    setequal(
      object$observation_moments$observation_model,
      c("regression_t_adjusted", "ideal_gaussian")
    ) &&
    nrow(object$observation_histogram) == 244L &&
    sum(object$observation_histogram$count) == 2L * J * length(TIME_GRID) &&
    nrow(object$observation_lag_correlation) ==
      2L * (length(TIME_GRID) - 1L) &&
    nrow(object$example_observations) ==
      length(TRUTH_GROUP_LEVELS) * length(TIME_GRID) &&
    nrow(object$digests) == 1L &&
    all(grepl(
      "^[[:xdigit:]]{32}$",
      unlist(object$digests[1L, grep("md5$", names(object$digests))])
    )) &&
    is.finite(object$elapsed_seconds) && object$elapsed_seconds > 0
}

log_message(
  "Starting or resuming R3 ideal-Gaussian replicates in ", PARTIAL_DIR, "."
)
for (truth_mechanism in RUN_TRUTH_MECHANISMS) {
  for (seed in RUN_SEED_LIST) {
    replicate_path <- file.path(
      REPLICATE_DIR,
      paste0(truth_mechanism, "_seed_", seed, ".rds")
    )
    if (file.exists(replicate_path)) {
      cached <- readRDS(replicate_path)
      if (!validate_replicate(cached, seed, truth_mechanism)) {
        stop("An existing replicate is invalid: ", replicate_path, call. = FALSE)
      }
      log_message("Reusing completed ", truth_mechanism, " seed ", seed, ".")
      next
    }
    log_message("Starting ", truth_mechanism, " seed ", seed, ".")
    replicate <- make_replicate(seed, truth_mechanism)
    if (!validate_replicate(replicate, seed, truth_mechanism)) {
      stop("A newly generated replicate failed validation.", call. = FALSE)
    }
    atomic_save_rds(replicate, replicate_path)
    log_message(
      "Completed ", truth_mechanism, " seed ", seed, " in ",
      sprintf("%.1f", replicate$elapsed_seconds), " seconds."
    )
    rm(replicate)
    invisible(gc(verbose = FALSE))
  }
}

replicate_paths <- unlist(lapply(TRUTH_MECHANISMS, function(truth_mechanism) {
  file.path(
    REPLICATE_DIR,
    paste0(truth_mechanism, "_seed_", SEED_LIST, ".rds")
  )
}), use.names = FALSE)
if (any(!file.exists(replicate_paths))) {
  log_message(
    "The requested subset is complete, but the ten-replicate cache is not yet complete."
  )
  quit(save = "no", status = 0L)
}
replicates <- lapply(replicate_paths, readRDS)
expected_keys <- expand.grid(
  seed = SEED_LIST,
  truth_mechanism = TRUTH_MECHANISMS,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
expected_keys <- expected_keys[order(
  match(expected_keys$truth_mechanism, TRUTH_MECHANISMS),
  match(expected_keys$seed, SEED_LIST)
), ]
rownames(expected_keys) <- NULL
observed_keys <- data.frame(
  seed = vapply(replicates, `[[`, integer(1), "seed"),
  truth_mechanism = vapply(replicates, `[[`, character(1), "truth_mechanism"),
  stringsAsFactors = FALSE
)
if (!identical(observed_keys, expected_keys)) {
  stop("The ten replicate keys are incomplete or misordered.", call. = FALSE)
}
for (index in seq_along(replicates)) {
  if (!validate_replicate(
    replicates[[index]],
    observed_keys$seed[[index]],
    observed_keys$truth_mechanism[[index]]
  )) {
    stop("A completed replicate failed final validation.", call. = FALSE)
  }
}

all_alpha <- do.call(rbind, lapply(replicates, `[[`, "functional_alpha"))
all_alpha_005 <- do.call(rbind, lapply(replicates, `[[`, "functional_alpha_005"))
all_pi0 <- do.call(rbind, lapply(replicates, `[[`, "estimated_pi0"))
all_truth_groups <- do.call(rbind, lapply(replicates, `[[`, "truth_group_counts"))
all_truth_maf <- do.call(rbind, lapply(replicates, `[[`, "truth_maf_balance"))
all_moments <- do.call(rbind, lapply(replicates, `[[`, "observation_moments"))
all_histogram <- do.call(rbind, lapply(replicates, `[[`, "observation_histogram"))
all_lag_correlation <- do.call(
  rbind,
  lapply(replicates, `[[`, "observation_lag_correlation")
)
all_digests <- do.call(rbind, lapply(replicates, `[[`, "digests"))
example_observations <- do.call(rbind, lapply(
  replicates[vapply(replicates, `[[`, integer(1), "seed") == SEED_LIST[[1L]]],
  `[[`,
  "example_observations"
))
mc_alpha <- summarize_mc_functional_alpha_curves(all_alpha)

middle_rows <- all_alpha[
  all_alpha$method == "FASH-IWP1-BF" &
    all_alpha$target == "middle" &
    is.finite(all_alpha$alpha) &
    all_alpha$alpha >= 0.05 - 1e-12 &
    all_alpha$alpha <= 0.20 + 1e-12,
  ,
  drop = FALSE
]
middle_rows$alpha <- round(middle_rows$alpha, digits = 12L)
if (nrow(middle_rows) !=
      length(SEED_LIST) * length(TRUTH_MECHANISMS) * 31L) {
  stop("The tolerance-aware Middle summary does not contain 310 rows.",
       call. = FALSE)
}
middle_curve <- do.call(rbind, lapply(
  split(middle_rows, list(middle_rows$truth_mechanism, middle_rows$alpha)),
  function(rows) {
    data.frame(
      truth_mechanism = rows$truth_mechanism[[1L]],
      alpha = rows$alpha[[1L]],
      mean_empirical_fsr = mean(rows$empirical_fsr),
      min_empirical_fsr = min(rows$empirical_fsr),
      max_empirical_fsr = max(rows$empirical_fsr),
      mean_power = mean(rows$power),
      stringsAsFactors = FALSE
    )
  }
))
middle_curve <- middle_curve[order(
  match(middle_curve$truth_mechanism, TRUTH_MECHANISMS),
  middle_curve$alpha
), ]
rownames(middle_curve) <- NULL
if (nrow(middle_curve) != 62L ||
    min(middle_curve$alpha) != 0.05 ||
    max(middle_curve$alpha) != 0.20 ||
    any(table(middle_curve$truth_mechanism) != 31L)) {
  stop("The production Middle curve failed its 62-row alpha contract.",
       call. = FALSE)
}
primary_summary <- do.call(rbind, lapply(
  split(middle_curve, middle_curve$truth_mechanism),
  function(rows) {
    excess <- rows$mean_empirical_fsr - rows$alpha
    index <- which.max(excess)
    data.frame(
      truth_mechanism = rows$truth_mechanism[[1L]],
      method = "FASH-IWP1-BF",
      target = "middle",
      alpha_min = min(rows$alpha),
      alpha_max = max(rows$alpha),
      maximum_mean_fsr_excess = excess[[index]],
      alpha_at_maximum = rows$alpha[[index]],
      mean_empirical_fsr_at_maximum = rows$mean_empirical_fsr[[index]],
      interpretation_threshold = NA_real_,
      stringsAsFactors = FALSE
    )
  }
))
primary_summary <- primary_summary[match(
  TRUTH_MECHANISMS, primary_summary$truth_mechanism
), ]
rownames(primary_summary) <- NULL

summary_objects <- list(
  all_replicate_functional_alpha_curves.csv = all_alpha,
  all_replicate_functional_alpha005.csv = all_alpha_005,
  functional_testing_mc_alpha_curve.csv = mc_alpha,
  all_replicate_pi0.csv = all_pi0,
  all_truth_group_counts.csv = all_truth_groups,
  truth_maf_balance.csv = all_truth_maf,
  standardized_error_moments.csv = all_moments,
  standardized_error_histogram.csv = all_histogram,
  standardized_error_lag_correlation.csv = all_lag_correlation,
  replicate_input_digests.csv = all_digests,
  ideal_middle_curve.csv = middle_curve,
  primary_middle_summary.csv = primary_summary
)
for (name in names(summary_objects)) {
  atomic_write_csv(summary_objects[[name]], file.path(SUMMARY_DIR, name))
}
atomic_save_rds(
  example_observations,
  file.path(PARTIAL_DIR, "example_observations.rds")
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
  schema_version = SCHEMA_VERSION,
  result_id = RESULT_ID,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  package_provenance = package_provenance,
  source_provenance = list(
    paths = c(
      helper = HELPER_PATH,
      simulation_functions = SIMULATION_FUNCTIONS_PATH,
      real_genotype_helper = REAL_GENOTYPE_HELPER_PATH,
      genotype_cache = GENOTYPE_CACHE_PATH,
      temporal_mixture_contract = TEMPORAL_MIXTURE_CONTRACT_PATH,
      runner = RUNNER_PATH
    ),
    sha256 = observed_source_hashes
  ),
  configuration = configuration,
  replicate_keys = observed_keys,
  artifact_sha256 = artifact_sha256
)
atomic_save_rds(manifest, file.path(PARTIAL_DIR, "manifest.rds"))
atomic_write_lines(
  c(
    paste0("result_id=", RESULT_ID),
    paste0("completed_at=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("schema_version=", SCHEMA_VERSION),
    paste0("fashr_version=", EXPECTED_FASHR_VERSION),
    paste0("fashr_remote_sha=", EXPECTED_FASHR_REMOTE_SHA),
    paste0("base_r3_result_id=", BASE_R3_RESULT_ID),
    "observation_model=ideal_gaussian_known_t_adjusted_se",
    "standard_error_source=R3 final t-adjusted SE regenerated from frozen seeds",
    paste0("middle_definition=", MIDDLE_EXPRESSION),
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
      "truth_group_counts=",
      paste(
        names(EXPECTED_TRUTH_GROUP_COUNTS),
        EXPECTED_TRUTH_GROUP_COUNTS,
        sep = ":",
        collapse = ";"
      )
    ),
    paste0(
      "temporal_mixture_contract_sha256=",
      observed_source_hashes[["temporal_mixture_contract"]]
    ),
    "summary_alpha_min=0.05",
    "summary_alpha_max=0.20",
    "summary_alpha_tolerance=1e-12",
    "summary_middle_rows=62",
    "functional_candidate_scope=full_universe",
    paste0("functional_candidate_universe_size=", J),
    "functional_posterior_pairing=common_random_seed_raw_bf",
    "replicates=10"
  ),
  file.path(PARTIAL_DIR, "complete.flag")
)
if (!file.rename(PARTIAL_DIR, FINAL_DIR)) {
  stop("Unable to promote the validated cache to ", FINAL_DIR, ".", call. = FALSE)
}
log_message("R3 ideal-Gaussian measurement experiment completed: ", FINAL_DIR)
