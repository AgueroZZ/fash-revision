#!/usr/bin/env Rscript

# Recompute only the FASH portions of the formal five-seed R1 and R2
# simulations with the pinned corrected fashr package. Direct-interaction
# tests are intentionally excluded and will be reused from the historical
# formal caches after this result family is downloaded.

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[[1L]] == length(args)) {
    return(default)
  }
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

require_input_file <- function(path, label) {
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
  path <- require_input_file(path, "SHA-256 input")
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

atomic_write_csv <- function(object, path) {
  temporary_path <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)
  utils::write.csv(object, temporary_path, row.names = FALSE, quote = TRUE)
  if (!file.rename(temporary_path, path)) {
    stop("Unable to atomically replace ", path, ".", call. = FALSE)
  }
  invisible(path)
}

capture_warnings <- function(expression) {
  warnings <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

EXPECTED_FASHR_VERSION <- "0.1.43"
EXPECTED_FASHR_REMOTE_SHA <-
  "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
EXPECTED_SIMULATION_FUNCTIONS_SHA256 <-
  "93f9a2c5606ae74763fb53e189af62c4d3b7c973aaf38d5b3eb3bf32f97a6487"
EXPECTED_REAL_GENOTYPE_HELPER_SHA256 <-
  "c03c01a188503336a77793c96f4e2d3ac7e0cbd56f4028b552da4b2f88e6b9d7"
EXPECTED_GENOTYPE_CACHE_SHA256 <-
  "81bbef5f323a0bab2ca993c782d8a9b7c63518b83c2cdb46ef7ed1d46f65af24"

RESULT_ID <- "r1_r2_fashr0143"
SCHEMA_VERSION <- "r1-r2-fashr0143-seed-v2"
SEED_LIST <- c(12345L, 22345L, 32345L, 42345L, 52345L)
GENOTYPE_DIGEST_METHOD <- "canonical-genotype-md5-v1"
EXPECTED_GENOTYPE_CONTENT_MD5 <- c(
  `12345` = "526a7318aa2af901e09252f5a6ca3c46",
  `22345` = "517faa30d5218a956f1be84f2567369c",
  `32345` = "9b9be3205d7db54dac31763492bcb2eb",
  `42345` = "9ef73aa94a061df868b1a951fe495d9f",
  `52345` = "7dab159b8453e2f66188ae313bfbd611"
)
J <- 6362L
N_DONORS <- 19L
N_COVARIATES <- 5L
TIME_GRID <- 0:15
CLASS_PROBS <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
EXPRESSION_NOISE_SD <- 1
DYNAMIC_MAIN_EFFECT_SD <- 1
NUM_BASIS <- 20L
PRED_STEP <- 1
PENALTY <- 10L
ALPHA_GRID <- seq(0, 0.20, by = 0.005)
NUM_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))
PREFLIGHT_ONLY <- as_flag(get_arg("--preflight-only", "false"))

if (is.na(NUM_CORES) || NUM_CORES < 1L) {
  stop("SLURM_CPUS_PER_TASK must be a positive integer.", call. = FALSE)
}

WORKSPACE_DIR <- Sys.getenv(
  "FASH_WORKSPACE",
  unset = "/project/mstephens/ziangzhang/fash/workspace"
)
WORKSPACE_DIR <- normalizePath(
  WORKSPACE_DIR,
  winslash = "/",
  mustWork = TRUE
)
RESULT_PARENT <- Sys.getenv(
  "FASH_R1_R2_RESULT_PARENT",
  unset = file.path(WORKSPACE_DIR, "results", "revision_simulations")
)
FINAL_DIR <- file.path(RESULT_PARENT, RESULT_ID)
PARTIAL_DIR <- file.path(RESULT_PARENT, paste0(RESULT_ID, "_partial"))
SIMULATION_FUNCTIONS_PATH <- Sys.getenv(
  "FASH_R1_R2_SIMULATION_FUNCTIONS",
  unset = file.path(
    WORKSPACE_DIR,
    "code", "revision_simulations", "shared", "simulation_functions.R"
  )
)
REAL_GENOTYPE_HELPER_PATH <- Sys.getenv(
  "FASH_R1_R2_REAL_GENOTYPE_HELPER",
  unset = file.path(
    WORKSPACE_DIR,
    "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
  )
)
GENOTYPE_CACHE_PATH <- Sys.getenv(
  "FASH_R1_R2_GENOTYPE_CACHE",
  unset = file.path(
    WORKSPACE_DIR,
    "inputs", "r1_r2_fashr0143", "genotype_samples.rds"
  )
)
RUNNER_PATH <- normalizePath(
  sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())][[1L]]),
  winslash = "/",
  mustWork = TRUE
)

SIMULATION_FUNCTIONS_PATH <- require_input_file(
  SIMULATION_FUNCTIONS_PATH,
  "Shared simulation helper"
)
REAL_GENOTYPE_HELPER_PATH <- require_input_file(
  REAL_GENOTYPE_HELPER_PATH,
  "Real-genotype helper"
)
GENOTYPE_CACHE_PATH <- require_input_file(
  GENOTYPE_CACHE_PATH,
  "Formal five-seed genotype cache"
)

observed_source_hashes <- list(
  simulation_functions = sha256_file(SIMULATION_FUNCTIONS_PATH),
  real_genotype_helper = sha256_file(REAL_GENOTYPE_HELPER_PATH),
  genotype_cache = sha256_file(GENOTYPE_CACHE_PATH),
  runner = sha256_file(RUNNER_PATH)
)
expected_source_hashes <- list(
  simulation_functions = EXPECTED_SIMULATION_FUNCTIONS_SHA256,
  real_genotype_helper = EXPECTED_REAL_GENOTYPE_HELPER_SHA256,
  genotype_cache = EXPECTED_GENOTYPE_CACHE_SHA256
)
for (name in names(expected_source_hashes)) {
  if (!identical(observed_source_hashes[[name]], expected_source_hashes[[name]])) {
    stop(
      "Unexpected SHA-256 for ", name, ": expected ",
      expected_source_hashes[[name]], "; found ",
      observed_source_hashes[[name]], ".",
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
    find.package("fashr"),
    winslash = "/",
    mustWork = TRUE
  ),
  r_version = R.version.string,
  platform = R.version$platform
)
if (!identical(package_provenance$version, EXPECTED_FASHR_VERSION) ||
    !identical(package_provenance$remote_sha, EXPECTED_FASHR_REMOTE_SHA)) {
  stop(
    "Expected fashr ", EXPECTED_FASHR_VERSION, " at ",
    EXPECTED_FASHR_REMOTE_SHA, "; found version ",
    package_provenance$version, " at ", package_provenance$remote_sha, ".",
    call. = FALSE
  )
}

source(SIMULATION_FUNCTIONS_PATH)
source(REAL_GENOTYPE_HELPER_PATH)

genotype_cache <- readRDS(GENOTYPE_CACHE_PATH)
if (!is.list(genotype_cache) ||
    !all(c("configuration", "sample_ids", "samples") %in%
      names(genotype_cache)) ||
    !identical(genotype_cache$configuration$n_genes, J) ||
    !identical(genotype_cache$configuration$n_donors, N_DONORS) ||
    !all(SEED_LIST %in% genotype_cache$configuration$seed_list) ||
    !all(as.character(SEED_LIST) %in% names(genotype_cache$samples))) {
  stop("The formal genotype cache has an unexpected schema.", call. = FALSE)
}
if (!identical(
  names(EXPECTED_GENOTYPE_CONTENT_MD5),
  as.character(SEED_LIST)
)) {
  stop("The expected genotype-content digest vector is misaligned.", call. = FALSE)
}
genotype_content_digests <- stats::setNames(
  character(length(SEED_LIST)),
  as.character(SEED_LIST)
)
for (seed in SEED_LIST) {
  genotype_sample <- validate_real_genotype_sample(
    genotype_cache$samples[[as.character(seed)]],
    expected_genes = J,
    expected_donors = N_DONORS,
    maf_min = genotype_cache$configuration$maf_min
  )
  observed_digest <- genotype_content_md5(
    pair_key = genotype_sample$selection$pair_key,
    sample_ids = rownames(genotype_sample$G),
    G = genotype_sample$G
  )
  expected_digest <- EXPECTED_GENOTYPE_CONTENT_MD5[[as.character(seed)]]
  if (!identical(observed_digest, expected_digest)) {
    stop(
      "Unexpected canonical genotype-content digest for seed ", seed,
      ": expected ", expected_digest, "; found ", observed_digest, ".",
      call. = FALSE
    )
  }
  genotype_content_digests[[as.character(seed)]] <- observed_digest
}

COMMON_SD_GRID <- default_revision_grid()
EXPECTED_CLASS_COUNTS <- exact_proportional_counts(J, CLASS_PROBS)
TRUE_PI0 <- unname(
  sum(EXPECTED_CLASS_COUNTS[c("constant", "zero")]) / J
)
FASH_METHODS <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "FASH-linear-Raw",
  "FASH-linear-BF"
)

source_provenance <- list(
  paths = list(
    simulation_functions = SIMULATION_FUNCTIONS_PATH,
    real_genotype_helper = REAL_GENOTYPE_HELPER_PATH,
    genotype_cache = GENOTYPE_CACHE_PATH,
    runner = RUNNER_PATH
  ),
  sha256 = observed_source_hashes
)

log_message("R1/R2 targeted result: ", FINAL_DIR)
log_message(
  "fashr ", package_provenance$version,
  " at ", package_provenance$remote_sha
)
log_message("Using ", NUM_CORES, " CPU cores.")
log_message("Validated formal genotype cache for five seeds.")

if (PREFLIGHT_ONLY) {
  log_message("Preflight completed successfully.")
  quit(save = "no", status = 0L)
}

if (dir.exists(FINAL_DIR) || file.exists(FINAL_DIR)) {
  stop("Refusing to overwrite completed output: ", FINAL_DIR, call. = FALSE)
}
ensure_directory(RESULT_PARENT)
ensure_directory(PARTIAL_DIR)
ensure_directory(file.path(PARTIAL_DIR, "replicates", "r1"))
ensure_directory(file.path(PARTIAL_DIR, "replicates", "r2"))
ensure_directory(file.path(PARTIAL_DIR, "summary"))

exact_null_weight <- function(fit) {
  if (!is.list(fit) ||
      !is.numeric(fit$psd_grid) ||
      any(!is.finite(fit$psd_grid)) ||
      sum(fit$psd_grid == 0) != 1L ||
      !is.data.frame(fit$prior_weights) ||
      !all(c("psd", "prior_weight") %in% names(fit$prior_weights))) {
    stop("A FASH fit has an invalid exact-null representation.", call. = FALSE)
  }
  null_row <- which(fit$prior_weights$psd == 0)
  if (length(null_row) == 0L) {
    return(0)
  }
  if (length(null_row) != 1L) {
    stop("A FASH fit contains multiple exact-null prior rows.", call. = FALSE)
  }
  weight <- as.numeric(fit$prior_weights$prior_weight[[null_row]])
  if (!is.finite(weight) || weight < 0 || weight > 1) {
    stop("A FASH fit contains an invalid exact-null weight.", call. = FALSE)
  }
  weight
}

fit_fash_only <- function(genotype_sample,
                          covariates,
                          effect_sim,
                          expression_sim,
                          scenario,
                          seed) {
  eqtl_summary <- estimate_eqtl_summaries_from_genotypes(
    G = genotype_sample$G,
    expression = expression_sim$expression,
    covariates = covariates,
    apply_t_se_correction = TRUE
  )
  datasets <- make_fash_datasets_from_eqtl_summary(
    beta_hat = eqtl_summary$beta_hat,
    se = eqtl_summary$se,
    true_beta = effect_sim$beta_matrix,
    time_grid = TIME_GRID,
    unit_info = effect_sim$unit_info,
    scenario = scenario
  )
  unit_info <- attr(datasets, "unit_info")

  captured <- capture_warnings({
    fash_fits <- fit_fash_for_revision(
      datasets = datasets,
      orders = 1L,
      grid = COMMON_SD_GRID,
      num_basis = NUM_BASIS,
      penalty = PENALTY,
      pred_step = PRED_STEP,
      num_cores = NUM_CORES,
      apply_bf = TRUE,
      verbose = FALSE
    )
    linear_raw <- fit_linear_mixture_fash(
      datasets = datasets,
      grid = COMMON_SD_GRID,
      pred_step = PRED_STEP,
      penalty = PENALTY
    )
    linear_bf <- BF_update_linear_mixture_fash(linear_raw)
    list(
      fash_fits = fash_fits,
      linear_raw = linear_raw,
      linear_bf = linear_bf
    )
  })
  fitted <- captured$value

  for (fit_name in c("fash_iwp1_raw", "fash_iwp1_bf")) {
    fit <- fitted$fash_fits[[fit_name]]
    if (is.null(fit) ||
        !isTRUE(all.equal(fit$psd_grid, COMMON_SD_GRID, tolerance = 0)) ||
        !isTRUE(all.equal(fit$settings$pred_step, PRED_STEP, tolerance = 0)) ||
        !identical(as.integer(fit$settings$penalty), PENALTY) ||
        length(get_fash_lfdr(fit)) != J ||
        any(!is.finite(get_fash_lfdr(fit)))) {
      stop("An IWP1 fit failed the targeted R1/R2 contract.", call. = FALSE)
    }
  }
  validate_linear_mixture_fash(
    fitted$linear_raw,
    expected_grid = COMMON_SD_GRID,
    expected_pred_step = PRED_STEP,
    expected_penalty = PENALTY
  )
  validate_linear_mixture_fash(
    fitted$linear_bf,
    expected_grid = COMMON_SD_GRID,
    expected_pred_step = PRED_STEP,
    expected_penalty = PENALTY
  )

  result_table <- rbind(
    evaluate_lfdr_method(
      lfdr = get_fash_lfdr(fitted$fash_fits$fash_iwp1_raw),
      unit_info = unit_info,
      method = "FASH-IWP1-Raw",
      target = "dynamic",
      alpha = 0.05
    ),
    evaluate_lfdr_method(
      lfdr = get_fash_lfdr(fitted$fash_fits$fash_iwp1_bf),
      unit_info = unit_info,
      method = "FASH-IWP1-BF",
      target = "dynamic",
      alpha = 0.05
    ),
    evaluate_simplified_fash_fit(
      fit = fitted$linear_raw,
      unit_info = unit_info,
      alpha = 0.05,
      method = "FASH-linear-Raw"
    ),
    evaluate_simplified_fash_fit(
      fit = fitted$linear_bf,
      unit_info = unit_info,
      alpha = 0.05,
      method = "FASH-linear-BF"
    )
  )
  alpha_curve <- compute_alpha_curve(
    result_table = result_table,
    alpha_grid = ALPHA_GRID
  )
  alpha_curve$seed <- seed

  list(
    datasets = datasets,
    unit_info = unit_info,
    genotype = genotype_sample$G,
    variant_info = genotype_sample$variant_info,
    covariates = covariates,
    true_beta = effect_sim$beta_matrix,
    expression = expression_sim$expression,
    expression_simulation = expression_sim,
    eqtl_summary = eqtl_summary,
    se_correction_summary = summarize_se_correction(eqtl_summary),
    fash_fits = fitted$fash_fits,
    simplified_fit = fitted$linear_raw,
    simplified_fit_bf = fitted$linear_bf,
    result_table = result_table,
    alpha_curve = alpha_curve,
    warnings = captured$warnings
  )
}

base_configuration <- list(
  result_id = RESULT_ID,
  J = J,
  n_donors = N_DONORS,
  n_covariates = N_COVARIATES,
  time_grid = TIME_GRID,
  class_probs = CLASS_PROBS,
  expected_class_counts = EXPECTED_CLASS_COUNTS,
  expression_noise_sd = EXPRESSION_NOISE_SD,
  covariate_effect_sd = 0.5,
  intercept_sd = 0,
  dynamic_main_effect_sd = DYNAMIC_MAIN_EFFECT_SD,
  common_sd_grid = COMMON_SD_GRID,
  pred_step = PRED_STEP,
  penalty = PENALTY,
  num_basis = NUM_BASIS,
  seed_list = SEED_LIST,
  true_pi0 = TRUE_PI0,
  genotype_source = "paper-derived YRI DS dosage",
  genotype_selection_rule = genotype_cache$configuration$selection_rule,
  genotype_repeated_variant_rule =
    genotype_cache$configuration$repeated_variant_rule,
  genotype_maf_min = genotype_cache$configuration$maf_min,
  genotype_digest_method = GENOTYPE_DIGEST_METHOD,
  genotype_sample_ids = genotype_cache$sample_ids,
  package_provenance = package_provenance,
  source_provenance = source_provenance,
  direct_interaction_tests = "excluded; reuse formal historical direct caches"
)

build_artifact <- function(scenario_id, seed) {
  component_seeds <- revision_component_seeds(seed)
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

  if (identical(scenario_id, "r1")) {
    scenario <- paste0(
      "r1_real_genotype_one_per_gene_",
      "random_bspline_main_effect_dynamic_eqtl"
    )
    effect_sim <- simulate_variant_effect_curves(
      n_variants = J,
      time_grid = TIME_GRID,
      class_probs = CLASS_PROBS,
      scenario = scenario,
      dynamic_amplitude = 2,
      bspline_df = 6,
      bspline_coefficient_sd = 1,
      constant_sd = 1,
      dynamic_main_effect_sd = DYNAMIC_MAIN_EFFECT_SD,
      exact_class_counts = TRUE,
      seed = component_seeds[["functional_truth"]]
    )
    effect_sim <- reassign_effect_simulation_by_maf(
      effect_sim = effect_sim,
      maf = genotype_sample$variant_info$observed_maf,
      class_probs = CLASS_PROBS,
      seed = component_seeds[["classes"]],
      n_strata = 10L
    )
    scenario_configuration <- list(
      scenario = scenario,
      dynamic_amplitude = 2,
      bspline_df = 6,
      bspline_coefficient_sd = 1
    )
  } else if (identical(scenario_id, "r2")) {
    scenario <- paste0(
      "r2_real_genotype_one_per_gene_",
      "timed_cosine_one_two_three_peak_main_effect_dynamic_eqtl"
    )
    evaluation_grid <- seq(0, 15, by = 0.1)
    shape_cell_probs <- c(
      k1__spiky__single = 1 / 3,
      `k2__spiky__same-sign` = 1 / 6,
      `k2__spiky__alternating-sign` = 1 / 6,
      `k3__spiky__same-sign` = 1 / 6,
      `k3__spiky__alternating-sign` = 1 / 6
    )
    primary_time_groups <- c("early", "middle", "late")
    effect_sim <- simulate_raised_cosine_multipeak_effect_set(
      n_variants = J,
      time_grid = TIME_GRID,
      evaluation_grid = evaluation_grid,
      class_probs = CLASS_PROBS,
      width_levels = c(spiky = 1.5),
      spike_counts = 1:3,
      shape_cell_probs = shape_cell_probs,
      primary_time_groups = primary_time_groups,
      center_by_observed_mean = FALSE,
      validate_functional_labels = FALSE,
      switch_threshold = 0.25,
      relative_amplitude_range = c(0.35, 0.75),
      target_centered_rms = 0.9,
      baseline_sd = 1,
      constant_sd = 1,
      dynamic_baseline_sd = DYNAMIC_MAIN_EFFECT_SD,
      exact_class_counts = TRUE,
      seed = seed,
      class_seed = component_seeds[["classes"]],
      constant_seed = component_seeds[["constant_effects"]],
      shape_seed = component_seeds[["functional_truth"]],
      scenario = scenario
    )
    effect_sim <- reassign_effect_simulation_by_maf(
      effect_sim = effect_sim,
      maf = genotype_sample$variant_info$observed_maf,
      class_probs = CLASS_PROBS,
      seed = component_seeds[["classes"]],
      n_strata = 10L
    )
    dynamic <- effect_sim$unit_info$effect_class == "dynamic_bspline"
    observed_cell_counts <- table(effect_sim$unit_info$cell_id[dynamic])
    expected_shape_counts <- exact_proportional_counts(
      EXPECTED_CLASS_COUNTS[["dynamic_bspline"]],
      shape_cell_probs
    )
    timing_counts <- table(
      effect_sim$unit_info$cell_id[dynamic],
      effect_sim$unit_info$time_group[dynamic]
    )
    centered_rms <- sqrt(rowMeans(
      (
        effect_sim$beta_matrix[dynamic, , drop = FALSE] -
          rowMeans(effect_sim$beta_matrix[dynamic, , drop = FALSE])
      )^2
    ))
    peak_lengths <- vapply(
      effect_sim$unit_info$peak_centers[dynamic],
      length,
      integer(1)
    )
    if (sum(dynamic) != EXPECTED_CLASS_COUNTS[["dynamic_bspline"]] ||
        !identical(
          as.integer(observed_cell_counts[names(shape_cell_probs)]),
          as.integer(expected_shape_counts)
        ) ||
        !setequal(colnames(timing_counts), primary_time_groups) ||
        any(apply(timing_counts, 1, function(x) diff(range(x))) > 1L) ||
        any(peak_lengths != effect_sim$unit_info$spike_count[dynamic]) ||
        any(abs(centered_rms - 0.9) > 1e-10)) {
      stop("Invalid formal R2 truth for seed ", seed, ".", call. = FALSE)
    }
    scenario_configuration <- list(
      scenario = scenario,
      evaluation_grid = evaluation_grid,
      width_half = 1.5,
      target_centered_rms = 0.9,
      shape_cell_probs = shape_cell_probs,
      shape_cell_counts = expected_shape_counts,
      primary_time_groups = primary_time_groups,
      relative_amplitude_range = c(0.35, 0.75),
      center_by_observed_mean = FALSE,
      switch_threshold = 0.25
    )
  } else {
    stop("Unknown scenario id: ", scenario_id, call. = FALSE)
  }

  expression_sim <- simulate_eqtl_expression_from_genotypes(
    G = genotype_sample$G,
    beta_matrix = effect_sim$beta_matrix,
    time_grid = TIME_GRID,
    covariates = covariates,
    expression_noise_sd = EXPRESSION_NOISE_SD,
    covariate_effect_sd = 0.5,
    intercept_sd = 0,
    seed = component_seeds[["expression"]]
  )
  fitted <- fit_fash_only(
    genotype_sample = genotype_sample,
    covariates = covariates,
    effect_sim = effect_sim,
    expression_sim = expression_sim,
    scenario = scenario,
    seed = seed
  )
  true_pi0 <- dynamic_null_proportion(fitted$unit_info, target = "dynamic")
  if (!isTRUE(all.equal(true_pi0, TRUE_PI0))) {
    stop("The simulated dynamic-null proportion is invalid.", call. = FALSE)
  }
  truth_maf_balance <- summarize_truth_maf_balance(
    variant_info = genotype_sample$variant_info,
    unit_info = fitted$unit_info,
    seed = seed
  )
  pi0 <- data.frame(
    seed = seed,
    method = c("FASH-IWP1", "FASH-IWP1", "FASH-linear", "FASH-linear"),
    fit = c("Raw", "BF-corrected", "Raw", "BF-corrected"),
    estimated_pi0 = c(
      exact_null_weight(fitted$fash_fits$fash_iwp1_raw),
      exact_null_weight(fitted$fash_fits$fash_iwp1_bf),
      exact_null_weight(fitted$simplified_fit),
      exact_null_weight(fitted$simplified_fit_bf)
    ),
    stringsAsFactors = FALSE
  )
  linear_prior_weights <- rbind(
    extract_linear_mixture_prior_table(
      fitted$simplified_fit,
      seed = seed,
      fit_label = "Raw"
    ),
    extract_linear_mixture_prior_table(
      fitted$simplified_fit_bf,
      seed = seed,
      fit_label = "BF-corrected"
    )
  )
  linear_prior_summary <- rbind(
    summarize_linear_mixture_prior_fit(
      fitted$simplified_fit,
      seed = seed,
      fit_label = "Raw"
    ),
    summarize_linear_mixture_prior_fit(
      fitted$simplified_fit_bf,
      seed = seed,
      fit_label = "BF-corrected"
    )
  )

  artifact <- list(
    schema_version = SCHEMA_VERSION,
    result_id = RESULT_ID,
    scenario_id = scenario_id,
    scenario = scenario,
    seed = seed,
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package_provenance = package_provenance,
    source_provenance = source_provenance,
    configuration = c(base_configuration, scenario_configuration),
    component_seeds = component_seeds,
    genotype_digest = genotype_content_digests[[as.character(seed)]],
    selected_pair_keys = genotype_sample$selection$pair_key,
    genotype_selection_summary = genotype_sample$selection_summary,
    truth_maf_balance = truth_maf_balance,
    unit_info = fitted$unit_info,
    fash_fits = fitted$fash_fits,
    simplified_fit = fitted$simplified_fit,
    simplified_fit_bf = fitted$simplified_fit_bf,
    result_table = fitted$result_table,
    alpha_curve = fitted$alpha_curve,
    alpha_005 = fitted$alpha_curve[
      abs(fitted$alpha_curve$alpha - 0.05) < 1e-12,
      ,
      drop = FALSE
    ],
    pi0 = pi0,
    linear_prior_weights = linear_prior_weights,
    linear_prior_summary = linear_prior_summary,
    warnings = fitted$warnings
  )

  if (identical(scenario_id, "r2")) {
    peak_alpha_curve <- compute_dynamic_subgroup_alpha_curve(
      result_table = fitted$result_table,
      unit_info = fitted$unit_info,
      subgroup_var = "spike_count",
      alpha_grid = ALPHA_GRID
    )
    peak_alpha_curve$seed <- seed
    artifact$peak_alpha_curve <- peak_alpha_curve
    artifact$peak_alpha_005 <- peak_alpha_curve[
      abs(peak_alpha_curve$alpha - 0.05) < 1e-12,
      ,
      drop = FALSE
    ]
    artifact$geometry <- cbind(
      seed = seed,
      summarize_dynamic_effect_geometry(
        effect_sim,
        dynamic_class = "dynamic_bspline"
      )
    )
  }

  if (identical(seed, SEED_LIST[[1L]])) {
    artifact$example_payload <- list(
      datasets = fitted$datasets,
      unit_info = fitted$unit_info,
      genotype = fitted$genotype,
      variant_info = fitted$variant_info,
      covariates = fitted$covariates,
      true_beta = fitted$true_beta,
      expression = fitted$expression,
      expression_simulation = fitted$expression_simulation,
      eqtl_summary = fitted$eqtl_summary,
      se_correction_summary = fitted$se_correction_summary,
      true_beta_evaluation = if (identical(scenario_id, "r2")) {
        effect_sim$beta_evaluation
      } else {
        NULL
      },
      evaluation_grid = if (identical(scenario_id, "r2")) {
        scenario_configuration$evaluation_grid
      } else {
        NULL
      },
      true_functionals = if (identical(scenario_id, "r2")) {
        effect_sim$true_functionals
      } else {
        NULL
      }
    )
  }

  artifact
}

validate_artifact <- function(artifact, scenario_id, seed) {
  required <- c(
    "schema_version", "result_id", "scenario_id", "scenario", "seed",
    "package_provenance", "source_provenance", "configuration",
    "component_seeds", "genotype_digest", "selected_pair_keys",
    "genotype_selection_summary", "truth_maf_balance", "unit_info",
    "fash_fits", "simplified_fit", "simplified_fit_bf", "result_table",
    "alpha_curve", "alpha_005", "pi0", "linear_prior_weights",
    "linear_prior_summary", "warnings"
  )
  if (!is.list(artifact) || !all(required %in% names(artifact)) ||
      !identical(artifact$schema_version, SCHEMA_VERSION) ||
      !identical(artifact$result_id, RESULT_ID) ||
      !identical(artifact$scenario_id, scenario_id) ||
      !identical(artifact$seed, seed) ||
      !identical(
        artifact$package_provenance$version,
        EXPECTED_FASHR_VERSION
      ) ||
      !identical(
        artifact$package_provenance$remote_sha,
        EXPECTED_FASHR_REMOTE_SHA
      ) ||
      !identical(
        artifact$source_provenance$sha256$simulation_functions,
        EXPECTED_SIMULATION_FUNCTIONS_SHA256
      ) ||
      !identical(
        artifact$source_provenance$sha256$real_genotype_helper,
        EXPECTED_REAL_GENOTYPE_HELPER_SHA256
      ) ||
      !identical(
        artifact$source_provenance$sha256$genotype_cache,
        EXPECTED_GENOTYPE_CACHE_SHA256
      ) ||
      !identical(
        artifact$source_provenance$sha256$runner,
        observed_source_hashes$runner
      ) ||
      !identical(
        artifact$configuration$genotype_digest_method,
        GENOTYPE_DIGEST_METHOD
      ) ||
      !identical(
        artifact$genotype_digest,
        genotype_content_digests[[as.character(seed)]]
      ) ||
      !identical(
        artifact$selected_pair_keys,
        genotype_cache$samples[[as.character(seed)]]$selection$pair_key
      ) ||
      nrow(artifact$unit_info) != J ||
      !setequal(unique(artifact$result_table$method), FASH_METHODS) ||
      !setequal(unique(artifact$alpha_curve$method), FASH_METHODS) ||
      nrow(artifact$alpha_005) != length(FASH_METHODS) ||
      nrow(artifact$pi0) != 4L ||
      any(!is.finite(artifact$pi0$estimated_pi0)) ||
      any(artifact$pi0$estimated_pi0 < 0) ||
      any(artifact$pi0$estimated_pi0 > 1) ||
      !is.character(artifact$warnings)) {
    return(FALSE)
  }
  if (identical(seed, SEED_LIST[[1L]]) && is.null(artifact$example_payload)) {
    return(FALSE)
  }
  if (identical(scenario_id, "r2") &&
      !all(c("peak_alpha_curve", "peak_alpha_005", "geometry") %in%
        names(artifact))) {
    return(FALSE)
  }
  for (fit_name in c("fash_iwp1_raw", "fash_iwp1_bf")) {
    fit <- artifact$fash_fits[[fit_name]]
    if (is.null(fit) ||
        !isTRUE(all.equal(fit$psd_grid, COMMON_SD_GRID, tolerance = 0)) ||
        length(get_fash_lfdr(fit)) != J ||
        any(!is.finite(get_fash_lfdr(fit)))) {
      return(FALSE)
    }
  }
  linear_validation <- tryCatch({
    validate_linear_mixture_fash(
      artifact$simplified_fit,
      expected_grid = COMMON_SD_GRID,
      expected_pred_step = PRED_STEP,
      expected_penalty = PENALTY
    )
    validate_linear_mixture_fash(
      artifact$simplified_fit_bf,
      expected_grid = COMMON_SD_GRID,
      expected_pred_step = PRED_STEP,
      expected_penalty = PENALTY
    )
    TRUE
  }, error = function(condition) FALSE)
  isTRUE(linear_validation)
}

artifact_paths <- list(r1 = character(), r2 = character())
for (scenario_id in c("r1", "r2")) {
  for (seed in SEED_LIST) {
    artifact_path <- file.path(
      PARTIAL_DIR,
      "replicates",
      scenario_id,
      paste0("seed_", seed, ".rds")
    )
    artifact_paths[[scenario_id]] <- c(
      artifact_paths[[scenario_id]],
      artifact_path
    )
    if (file.exists(artifact_path)) {
      artifact <- readRDS(artifact_path)
      if (!validate_artifact(artifact, scenario_id, seed)) {
        stop(
          "Existing partial artifact failed validation: ", artifact_path,
          call. = FALSE
        )
      }
      log_message("Reusing completed ", scenario_id, " seed ", seed, ".")
      rm(artifact)
      invisible(gc(verbose = FALSE))
      next
    }

    started <- proc.time()[["elapsed"]]
    log_message("Starting ", scenario_id, " seed ", seed, ".")
    artifact <- build_artifact(scenario_id, seed)
    if (!validate_artifact(artifact, scenario_id, seed)) {
      stop(
        "New artifact failed validation for ", scenario_id,
        " seed ", seed, ".",
        call. = FALSE
      )
    }
    atomic_save_rds(artifact, artifact_path)
    elapsed <- proc.time()[["elapsed"]] - started
    log_message(
      "Completed ", scenario_id, " seed ", seed,
      " in ", sprintf("%.1f", elapsed), " seconds; warnings=",
      length(artifact$warnings), "."
    )
    rm(artifact)
    invisible(gc(verbose = FALSE))
  }
}

summarize_peak_power <- function(rows, confidence_level = 0.95) {
  groups <- split(
    rows,
    list(rows$subgroup_value, rows$method, rows$alpha),
    drop = TRUE
  )
  out <- lapply(groups, function(x) {
    power <- summarize_mc_values(x$power, confidence_level)
    true_positives <- summarize_mc_values(x$true_positives, confidence_level)
    data.frame(
      spike_count = as.integer(x$subgroup_value[[1L]]),
      shape_profile = paste0(x$subgroup_value[[1L]], "-peak"),
      method = x$method[[1L]],
      alpha = x$alpha[[1L]],
      n_dynamic = x$n_dynamic[[1L]],
      n_replications = length(unique(x$seed)),
      mean_true_positives = true_positives[["mean"]],
      mean_power = power[["mean"]],
      power_sd = power[["sd"]],
      power_mc_se = power[["se"]],
      power_ci_lower = pmax(0, power[["lower"]]),
      power_ci_upper = pmin(1, power[["upper"]]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

summary_paths <- character()
for (scenario_id in c("r1", "r2")) {
  alpha_parts <- list()
  pi0_parts <- list()
  prior_weight_parts <- list()
  prior_summary_parts <- list()
  selection_parts <- list()
  maf_parts <- list()
  peak_parts <- list()
  geometry_parts <- list()

  for (index in seq_along(SEED_LIST)) {
    artifact <- readRDS(artifact_paths[[scenario_id]][[index]])
    if (!validate_artifact(artifact, scenario_id, SEED_LIST[[index]])) {
      stop("Artifact became invalid during aggregation.", call. = FALSE)
    }
    alpha_parts[[index]] <- artifact$alpha_curve
    pi0_parts[[index]] <- artifact$pi0
    prior_weight_parts[[index]] <- artifact$linear_prior_weights
    prior_summary_parts[[index]] <- artifact$linear_prior_summary
    selection <- artifact$genotype_selection_summary
    selection$seed <- artifact$seed
    selection_parts[[index]] <- selection
    maf_parts[[index]] <- artifact$truth_maf_balance
    if (identical(scenario_id, "r2")) {
      peak_parts[[index]] <- artifact$peak_alpha_curve
      geometry_parts[[index]] <- artifact$geometry
    }
    rm(artifact)
  }

  all_alpha <- do.call(rbind, alpha_parts)
  all_pi0 <- do.call(rbind, pi0_parts)
  all_prior_weights <- do.call(rbind, prior_weight_parts)
  all_prior_summary <- do.call(rbind, prior_summary_parts)
  all_selection <- do.call(rbind, selection_parts)
  all_maf <- do.call(rbind, maf_parts)
  mc_alpha <- summarize_mc_alpha_curves(all_alpha)
  mc_alpha_005 <- mc_alpha[abs(mc_alpha$alpha - 0.05) < 1e-12, ]
  mc_pi0 <- summarize_mc_pi0(all_pi0)

  outputs <- list(
    all_replicate_fash_alpha_curves = all_alpha,
    fash_mc_alpha_curve = mc_alpha,
    fash_mc_alpha005_summary = mc_alpha_005,
    all_replicate_pi0 = all_pi0,
    mc_pi0_summary = mc_pi0,
    all_replicate_linear_prior_weights = all_prior_weights,
    all_replicate_linear_prior_summary = all_prior_summary,
    genotype_selection_summary = all_selection,
    truth_maf_balance = all_maf
  )
  if (identical(scenario_id, "r2")) {
    all_peak <- do.call(rbind, peak_parts)
    peak_mc <- summarize_peak_power(all_peak)
    outputs$all_replicate_peak_alpha_curves <- all_peak
    outputs$peak_mc_alpha_curve <- peak_mc
    outputs$peak_mc_alpha005_summary <- peak_mc[
      abs(peak_mc$alpha - 0.05) < 1e-12,
    ]
    outputs$all_replicate_geometry <- do.call(rbind, geometry_parts)
  }

  for (name in names(outputs)) {
    path <- file.path(
      PARTIAL_DIR,
      "summary",
      paste0(scenario_id, "_", name, ".csv")
    )
    atomic_write_csv(outputs[[name]], path)
    summary_paths <- c(summary_paths, path)
  }
}

all_artifact_paths <- unlist(artifact_paths, use.names = FALSE)
manifest <- list(
  schema_version = "r1-r2-fashr0143-manifest-v1",
  result_id = RESULT_ID,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  package_provenance = package_provenance,
  source_provenance = source_provenance,
  configuration = base_configuration,
  artifact_sha256 = setNames(
    lapply(all_artifact_paths, sha256_file),
    sub(paste0("^", PARTIAL_DIR, "/"), "", all_artifact_paths)
  ),
  summary_sha256 = setNames(
    lapply(summary_paths, sha256_file),
    sub(paste0("^", PARTIAL_DIR, "/"), "", summary_paths)
  )
)
manifest_path <- file.path(PARTIAL_DIR, "manifest.rds")
atomic_save_rds(manifest, manifest_path)
writeLines(
  c(
    paste0("result_id=", RESULT_ID),
    paste0("completed_at=", manifest$generated_at),
    paste0("fashr_version=", package_provenance$version),
    paste0("fashr_remote_sha=", package_provenance$remote_sha),
    paste0("r1_replicates=", length(artifact_paths$r1)),
    paste0("r2_replicates=", length(artifact_paths$r2)),
    "direct_interaction_tests=reused_from_historical_formal_caches"
  ),
  con = file.path(PARTIAL_DIR, "complete.flag")
)

if (!file.rename(PARTIAL_DIR, FINAL_DIR)) {
  stop(
    "All outputs were written, but the partial directory could not be promoted to ",
    FINAL_DIR, ".",
    call. = FALSE
  )
}
log_message("Targeted R1/R2 fashr 0.1.43 run completed: ", FINAL_DIR)
