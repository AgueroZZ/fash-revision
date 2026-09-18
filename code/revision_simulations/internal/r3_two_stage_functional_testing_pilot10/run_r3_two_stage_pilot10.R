#!/usr/bin/env Rscript

# Reuse five completed records and compute five added seeds for formal R3.

get_arg <- function(name, default) {
  arguments <- commandArgs(trailingOnly = TRUE)
  index <- match(name, arguments)
  if (is.na(index)) return(default)
  if (index == length(arguments)) stop("Missing value for ", name, call. = FALSE)
  arguments[[index + 1L]]
}

root <- normalizePath(
  get_arg("--project-root", "."),
  winslash = "/",
  mustWork = TRUE
)
experiment <- Sys.getenv(
  "FASH_R3_PILOT10_CODE",
  unset = file.path(
    root,
    "code/revision_simulations/internal/r3_two_stage_functional_testing_pilot10"
  )
)
base_experiment <- Sys.getenv(
  "FASH_R3_TWO_STAGE_BASE_CODE",
  unset = file.path(
    root,
    "code/revision_simulations/internal/r3_two_stage_functional_testing"
  )
)
source(file.path(experiment, "pilot10_helpers.R"))
contract <- r3ts10_contract()
mode <- get_arg("--mode", "run")
r3ts10_assert(mode %in% c("preflight", "run"), "Mode must be preflight or run.")
cores <- as.integer(get_arg(
  "--num-cores",
  Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1")
))
r3ts10_assert(
  length(cores) == 1L && is.finite(cores) && cores >= 1L,
  "The worker count must be a positive integer."
)
if (mode == "run" && identical(Sys.info()[["sysname"]], "Darwin")) {
  stop("Full R3 production must run on Midway3.", call. = FALSE)
}

log_message <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")
  flush.console()
}

sources <- c(
  simulation = Sys.getenv(
    "FASH_R3_SIMULATION_FUNCTIONS",
    file.path(
      root,
      "code/revision_simulations/r3_r4_fashr0143/source_snapshots/r3_full_universe_functional_simulation_functions.R"
    )
  ),
  genotype_helper = Sys.getenv(
    "FASH_R3_REAL_GENOTYPE_HELPER",
    file.path(root, "code/revision_simulations/shared/real_genotype_one_per_gene.R")
  ),
  digest_helper = Sys.getenv(
    "FASH_R3_DIGEST_HELPER",
    file.path(
      root,
      "code/revision_simulations/internal/r3_ideal_gaussian_measurement/ideal_gaussian_measurement.R"
    )
  ),
  base_genotypes = Sys.getenv(
    "FASH_R3_BASE_GENOTYPE_CACHE",
    file.path(
      root,
      "output/revision_simulations/shared/real_genotype_one_per_gene_J6362_pilot5/genotype_samples.rds"
    )
  ),
  base_two_stage_helper = file.path(base_experiment, "two_stage_helpers.R"),
  base_reconstruction = file.path(base_experiment, "reconstruct_r3.R")
)
source_hashes <- vapply(sources, r3ts10_sha256, character(1))
r3ts10_assert(
  identical(source_hashes, contract$base_source_sha[names(source_hashes)]),
  "A frozen R3 source or base genotype cache hash changed."
)
for (name in c("simulation", "genotype_helper", "digest_helper")) {
  source(sources[[name]])
}
source(file.path(experiment, "genotype_cache_helpers.R"))

base_environment <- new.env(parent = globalenv())
sys.source(sources[["base_two_stage_helper"]], envir = base_environment)
r3ts10_require_base_environment(base_environment)
reconstruction_environment <- new.env(parent = globalenv())
sys.source(sources[["base_reconstruction"]], envir = reconstruction_environment)
reconstruction_environment$r3ts_contract <- r3ts10_contract
reconstruction_environment$r3ts_assert <- r3ts10_assert
reconstruction_environment$r3ts_stage1 <- base_environment$r3ts_stage1
reconstruction_environment$r3ts_stage2 <- base_environment$r3ts_stage2
reconstruction_environment$r3ts_metrics <- base_environment$r3ts_metrics

description <- utils::packageDescription("fashr")
r3ts10_assert(
  identical(description$Version, contract$package_version) &&
    identical(description$RemoteSha, contract$package_sha),
  "The installed fashr build differs from the frozen R3 build."
)
posterior_draws <- formals(
  getS3method("predict", "fash", envir = asNamespace("fashr"))
)$M
r3ts10_assert(
  isTRUE(all.equal(posterior_draws, contract$posterior_draws)),
  "The posterior draw count differs from the frozen R3 sampler."
)

formal_dir <- Sys.getenv(
  "FASH_R3_FORMAL_DIR",
  file.path(root, "output/revision_simulations/mc", contract$formal_id)
)
base_result_dir <- Sys.getenv(
  "FASH_R3_TWO_STAGE_BASE_RESULT",
  file.path(
    root,
    "output/revision_simulations/internal",
    contract$base_result_id
  )
)
pair_source_path <- Sys.getenv(
  "FASH_R3_PILOT10_PAIR_SOURCE",
  file.path(
    root,
    "output/revision_simulations/internal/r3_two_stage_pilot10_development_20260904/pair_key_source.rds"
  )
)
genotype_dir <- Sys.getenv(
  "FASH_R3_PILOT10_GENOTYPE_DIR",
  file.path(
    root,
    "output/revision_simulations/shared",
    contract$genotype_id
  )
)
genotype_path <- file.path(genotype_dir, "genotype_samples.rds")

validate_manifest_artifacts <- function(directory, manifest) {
  paths <- file.path(directory, names(manifest$artifact_sha256))
  actual <- vapply(paths, r3ts10_sha256, character(1))
  r3ts10_assert(
    identical(unname(actual), unname(manifest$artifact_sha256)),
    paste("Artifact hashes changed in", directory)
  )
  invisible(TRUE)
}

r3ts10_assert(
  file.exists(file.path(formal_dir, "complete.flag")) &&
    file.exists(file.path(formal_dir, "manifest.rds")) &&
    file.exists(file.path(formal_dir, "configuration.rds")),
  "The frozen formal R3 cache is incomplete."
)
formal_manifest_path <- file.path(formal_dir, "manifest.rds")
r3ts10_assert(
  identical(
    r3ts10_sha256(formal_manifest_path),
    contract$base_source_sha[["formal_result_manifest"]]
  ),
  "The frozen formal R3 manifest changed."
)
formal_manifest <- readRDS(formal_manifest_path)
r3ts10_assert(
  identical(formal_manifest$result_id, contract$formal_id),
  "The formal cache has the wrong result identity."
)
validate_manifest_artifacts(formal_dir, formal_manifest)
formal_configuration <- readRDS(file.path(formal_dir, "configuration.rds"))
r3ts10_validate_formal_configuration(formal_configuration)

r3ts10_assert(
  file.exists(file.path(base_result_dir, "complete.flag")) &&
    file.exists(file.path(base_result_dir, "manifest.rds")),
  "The completed pilot-5 two-stage cache is missing."
)
base_manifest_path <- file.path(base_result_dir, "manifest.rds")
r3ts10_assert(
  identical(
    r3ts10_sha256(base_manifest_path),
    contract$base_source_sha[["base_result_manifest"]]
  ),
  "The completed pilot-5 two-stage manifest changed."
)
base_manifest <- readRDS(base_manifest_path)
r3ts10_assert(
  identical(base_manifest$result_id, contract$base_result_id) &&
    setequal(base_manifest$replicate_keys, base_environment$r3ts_expected_keys()),
  "The pilot-5 two-stage cache lacks a completed replicate."
)
validate_manifest_artifacts(base_result_dir, base_manifest)

r3ts10_assert(
  file.exists(pair_source_path) &&
    file.exists(genotype_path) &&
    file.exists(file.path(genotype_dir, "complete.flag")) &&
    file.exists(file.path(genotype_dir, "manifest.rds")),
  "The completed pilot-10 genotype cache is missing."
)
base_genotypes <- readRDS(sources[["base_genotypes"]])
pair_source <- readRDS(pair_source_path)
genotypes <- readRDS(genotype_path)
genotype_validation <- r3ts10g_validate_cache(
  genotypes,
  base_genotypes,
  pair_source,
  pair_source_sha256 = r3ts10_sha256(pair_source_path),
  base_cache_sha256 = source_hashes[["base_genotypes"]]
)
genotype_manifest <- readRDS(file.path(genotype_dir, "manifest.rds"))
validate_manifest_artifacts(genotype_dir, genotype_manifest)
for (seed in contract$base_seeds) {
  key <- as.character(seed)
  r3ts10_assert(
    identical(
      genotype_content_md5(
        genotypes$samples[[key]]$selection$pair_key,
        rownames(genotypes$samples[[key]]$G),
        genotypes$samples[[key]]$G
      ),
      formal_configuration$genotype_content_digests[[key]]
    ),
    paste("A retained formal genotype changed for seed", seed, ".")
  )
}

runner_argument <- grep("^--file=", commandArgs(), value = TRUE)
r3ts10_assert(length(runner_argument) == 1L, "The runner source path is unavailable.")
runner_path <- normalizePath(
  sub("^--file=", "", runner_argument),
  winslash = "/",
  mustWork = TRUE
)
own_sources <- c(
  pilot10_helpers = file.path(experiment, "pilot10_helpers.R"),
  genotype_cache_helpers = file.path(experiment, "genotype_cache_helpers.R"),
  runner = runner_path
)
own_hashes <- vapply(own_sources, r3ts10_sha256, character(1))
package <- list(
  version = description$Version,
  remote_sha = description$RemoteSha,
  r_version = R.version.string,
  platform = R.version$platform
)
configuration <- list(
  contract = contract,
  formal_configuration = formal_configuration,
  package = package,
  frozen_source_sha256 = source_hashes,
  experiment_source_sha256 = own_hashes,
  formal_manifest_sha256 = r3ts10_sha256(formal_manifest_path),
  base_result_manifest_sha256 = r3ts10_sha256(base_manifest_path),
  genotype_manifest_sha256 = r3ts10_sha256(file.path(genotype_dir, "manifest.rds")),
  genotype_cache_sha256 = r3ts10_sha256(genotype_path),
  pair_key_source_sha256 = r3ts10_sha256(pair_source_path),
  fit_universe = "all_units",
  prior_refit_after_screening = FALSE,
  functional_candidate_scope = "dynamic_fdr_screen",
  stage1_method = "FASH-IWP1-BF",
  stage1_alpha_grid = contract$alpha_grid,
  stage2_q_dyn = contract$q_dyn,
  stage2_q_func = contract$q_func,
  posterior_draws = posterior_draws,
  num_cores = cores,
  reuse = list(
    migrated_base_replicates = 10L,
    newly_computed_replicates = 10L
  )
)
configuration_token <- r3ts10_serialized_object_md5(configuration)
log_message(
  "Preflight passed: formal configuration, pilot-5 records, ten-seed genotypes, ",
  "fashr build, and frozen sources."
)
log_message(
  "Stage 1 uses 40 alpha values; Stage 2 uses q_dyn=q_func=0.05; posterior M=",
  posterior_draws,
  "; workers=",
  cores,
  "."
)
if (mode == "preflight") quit(save = "no", status = 0L)

result_parent <- Sys.getenv(
  "FASH_R3_PILOT10_RESULT_PARENT",
  file.path(root, "output/revision_simulations/mc")
)
final_dir <- file.path(result_parent, contract$result_id)
partial_dir <- paste0(final_dir, ".partial")
checkpoint_root <- Sys.getenv("FASH_R3_PILOT10_CHECKPOINT_ROOT", unset = "")
r3ts10_assert(
  nzchar(checkpoint_root),
  "Set a server-only FASH_R3_PILOT10_CHECKPOINT_ROOT."
)
r3ts10_assert(
  !dir.exists(final_dir),
  "Refusing to overwrite an existing completed pilot-10 result directory."
)
dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)
configuration_path <- file.path(partial_dir, "configuration.rds")
if (file.exists(configuration_path)) {
  r3ts10_assert(
    identical(readRDS(configuration_path), configuration),
    "The partial result belongs to a different configuration."
  )
} else {
  r3ts10_write_rds(configuration, configuration_path)
}
checkpoint_token <- r3ts10_serialized_object_md5(list(
  formal_configuration = formal_configuration,
  frozen_source_sha256 = source_hashes,
  experiment_source_sha256 = own_hashes,
  genotype_cache_sha256 = configuration$genotype_cache_sha256,
  package = package,
  num_cores = cores
))

for (mechanism in contract$mechanisms) {
  for (seed in contract$seeds) {
    key <- paste0(mechanism, "_seed_", seed)
    result_path <- file.path(partial_dir, "replicates", paste0(key, ".rds"))
    if (file.exists(result_path)) {
      saved <- readRDS(result_path)
      r3ts10_assert(
        identical(saved$configuration_token, configuration_token),
        "A completed partial replicate has an incompatible configuration."
      )
      r3ts10_validate_record(
        saved,
        base_environment = base_environment,
        cfsr_function = r3ts10_functional_cfsr_table
      )
      log_message("Reusing completed pilot-10 record: ", key)
      next
    }

    if (seed %in% contract$base_seeds) {
      source_name <- file.path("replicates", paste0(key, ".rds"))
      source_path <- file.path(base_result_dir, source_name)
      source_record <- readRDS(source_path)
      source_sha256 <- base_manifest$artifact_sha256[[source_name]]
      r3ts10_assert(
        identical(r3ts10_sha256(source_path), unname(source_sha256)),
        paste("The reusable pilot-5 record changed:", key)
      )
      migrated <- r3ts10_migrate_base_record(
        source_record,
        configuration_token = configuration_token,
        source_sha256 = source_sha256,
        base_environment = base_environment,
        cfsr_function = r3ts10_functional_cfsr_table
      )
      r3ts10_write_rds(migrated, result_path)
      log_message("Validated and migrated completed pilot-5 record: ", key)
      next
    }

    started <- proc.time()[["elapsed"]]
    checkpoint_dir <- file.path(checkpoint_root, key)
    input_path <- file.path(checkpoint_dir, "inputs.rds")
    fit_path <- file.path(checkpoint_dir, "fits.rds")
    genotype_sample <- genotypes$samples[[as.character(seed)]]
    reused_inputs <- file.exists(input_path)
    if (reused_inputs) {
      saved <- readRDS(input_path)
      r3ts10_assert(
        identical(saved$token, checkpoint_token),
        "An input checkpoint has an incompatible configuration."
      )
      inputs <- saved$inputs
    } else {
      log_message("Generating frozen truth and genotype-level observations: ", key)
      inputs <- reconstruction_environment$r3ts_reconstruct(
        formal_configuration,
        genotype_sample,
        seed,
        mechanism
      )
      r3ts10_validate_added_inputs(inputs, genotype_sample, formal_configuration)
      r3ts10_write_rds(list(token = checkpoint_token, inputs = inputs), input_path)
    }
    r3ts10_validate_added_inputs(inputs, genotype_sample, formal_configuration)

    reused_fit <- file.exists(fit_path)
    if (reused_fit) {
      saved <- readRDS(fit_path)
      r3ts10_assert(
        identical(saved$token, checkpoint_token),
        "A fit checkpoint has an incompatible configuration."
      )
      fits <- saved$fits
    } else {
      log_message("Fitting IWP1 to all 6,362 units and applying the BF update: ", key)
      fits <- reconstruction_environment$r3ts_fit_all(
        inputs,
        formal_configuration,
        num_cores = cores
      )
      r3ts10_write_rds(list(token = checkpoint_token, fits = fits), fit_path)
    }
    fit <- fits$fash_iwp1_bf
    fitted_pi0 <- constant_component_prior_weight(fit)
    r3ts10_assert(
      length(fitted_pi0) == 1L && is.finite(fitted_pi0) &&
        fitted_pi0 >= 0 && fitted_pi0 <= 1,
      "The BF-adjusted constant-component weight is invalid."
    )
    log_message("Applying dynamic screening and conditional functional testing: ", key)
    inference <- reconstruction_environment$r3ts_infer(
      fit,
      inputs,
      formal_configuration,
      num_cores = cores
    )
    result <- list(
      schema = contract$schema,
      configuration_token = configuration_token,
      seed = seed,
      truth_mechanism = mechanism,
      selected_pair_keys = genotype_sample$selection$pair_key,
      input_digests = inputs$input_digests,
      inference = inference,
      fitted_n_units = length(fit$fash_data$data_list),
      fitted_bf_pi0 = fitted_pi0,
      reuse = c(
        genotype_cache = TRUE,
        frozen_source_configuration = TRUE,
        input_checkpoint = reused_inputs,
        fit_checkpoint = reused_fit
      ),
      validation = list(
        kind = "new_seed_frozen_generation",
        input_validation = "passed",
        fit_universe = contract$n_units,
        functional_candidate_scope = "dynamic_fdr_screen"
      ),
      checkpoint_paths = c(inputs = input_path, fits = fit_path),
      elapsed_seconds = proc.time()[["elapsed"]] - started
    )
    r3ts10_validate_record(
      result,
      base_environment = base_environment,
      cfsr_function = r3ts10_functional_cfsr_table
    )
    r3ts10_write_rds(result, result_path)
    log_message(
      "Completed ", key,
      ": Stage-1 calls=", length(inference$candidates),
      "; functional calls=", paste(inference$metrics$functional$calls, collapse = ",")
    )
    rm(list = intersect(
      c("saved", "inputs", "fits", "fit", "inference", "result"),
      ls()
    ))
    invisible(gc())
  }
}

replicate_paths <- file.path(partial_dir, "replicates", r3ts10_replicate_names())
r3ts10_assert(all(file.exists(replicate_paths)), "One or more pilot-10 records are missing.")
replicates <- lapply(replicate_paths, readRDS)
invisible(lapply(replicates, function(record) {
  r3ts10_validate_record(
    record,
    base_environment = base_environment,
    cfsr_function = r3ts10_functional_cfsr_table
  )
}))

stage1_alpha_by_seed <- do.call(rbind, lapply(
  replicates,
  r3ts10_stage1_alpha_by_seed,
  base_environment = base_environment
))
rownames(stage1_alpha_by_seed) <- NULL
functional_by_seed <- do.call(rbind, lapply(
  replicates,
  function(record) record$inference$metrics$functional
))
rownames(functional_by_seed) <- NULL
summaries <- r3ts10_summarize(stage1_alpha_by_seed, functional_by_seed)
stage1_at_005 <- r3ts10_stage1_at_alpha(summaries$stage1, contract$q_dyn)

summary_dir <- file.path(partial_dir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
tables <- list(
  stage1_alpha_by_seed = stage1_alpha_by_seed,
  stage1_mc_summary = summaries$stage1,
  stage1_alpha005_summary = stage1_at_005,
  functional_by_seed = functional_by_seed,
  functional_mc_summary = summaries$stage2,
  genotype_validation = genotype_validation
)
for (mechanism in contract$mechanisms) {
  tables[[paste0("stage2_", mechanism)]] <- r3ts10_stage2_display_values(
    summaries$stage2,
    mechanism
  )
}
tables$unconditional_power_audit <- functional_by_seed[c(
  "seed", "truth_mechanism", "target", "candidate_count", "calls",
  "true_calls", "true_target_total", "true_target_in_stage1", "empirical_power"
)]
tables$input_validation <- do.call(rbind, lapply(replicates, function(record) {
  data.frame(
    seed = record$seed,
    truth_mechanism = record$truth_mechanism,
    record_kind = record$validation$kind,
    genotype_content_md5 = unname(record$input_digests[["genotype_content_md5"]]),
    true_beta_md5 = unname(record$input_digests[["true_beta_md5"]]),
    adjusted_se_md5 = unname(record$input_digests[["adjusted_se_md5"]]),
    regression_beta_hat_md5 = unname(
      record$input_digests[["regression_beta_hat_md5"]]
    ),
    fitted_n_units = record$fitted_n_units,
    fitted_bf_pi0 = record$fitted_bf_pi0,
    stage1_calls_alpha005 = length(record$inference$candidates),
    source_replicate_sha256 = if (
      identical(record$validation$kind, "validated_pilot5_reuse")
    ) record$validation$source_replicate_sha256 else NA_character_,
    stringsAsFactors = FALSE
  )
}))
tables$all_stage1_units <- do.call(rbind, lapply(replicates, function(record) {
  fdr <- record$inference$fdr_table
  fdr <- fdr[match(seq_len(contract$n_units), fdr$index), , drop = FALSE]
  data.frame(
    seed = record$seed,
    truth_mechanism = record$truth_mechanism,
    pair_key = record$selected_pair_keys,
    unit_index = fdr$index,
    lfdr = fdr$lfdr,
    cumulative_fdr = fdr$FDR,
    selected_dynamic_alpha005 = fdr$index %in% record$inference$candidates,
    true_dynamic = record$inference$true_dynamic,
    record$inference$true_functionals,
    stringsAsFactors = FALSE
  )
}))
tables$all_stage2_candidates <- do.call(rbind, lapply(replicates, function(record) {
  table <- record$inference$stage2
  table$seed <- rep(record$seed, nrow(table))
  table$truth_mechanism <- rep(record$truth_mechanism, nrow(table))
  table$pair_key <- record$selected_pair_keys[table$index]
  table$true_functional <- record$inference$true_functionals[cbind(
    table$index,
    match(table$target, contract$targets)
  )]
  table
}))

for (name in names(tables)) {
  write.csv(
    tables[[name]],
    file.path(summary_dir, paste0(name, ".csv")),
    row.names = FALSE
  )
}
r3ts10_write_rds(tables, file.path(partial_dir, "summary_tables.rds"))
writeLines(capture.output(sessionInfo()), file.path(partial_dir, "sessionInfo.txt"))
artifacts <- list.files(partial_dir, recursive = TRUE, full.names = FALSE)
artifact_sha256 <- setNames(
  vapply(file.path(partial_dir, artifacts), r3ts10_sha256, character(1)),
  artifacts
)
manifest <- list(
  schema = contract$schema,
  result_id = contract$result_id,
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  configuration = configuration,
  configuration_token = configuration_token,
  replicate_keys = r3ts10_expected_keys(),
  artifact_sha256 = artifact_sha256
)
r3ts10_write_rds(manifest, file.path(partial_dir, "manifest.rds"))
writeLines(
  c(
    paste0("result_id=", contract$result_id),
    "replicates=20",
    "seeds=10",
    "truth_mechanisms=2",
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
  ),
  file.path(partial_dir, "complete.flag")
)
r3ts10_assert(
  file.rename(partial_dir, final_dir),
  "Failed to promote the completed pilot-10 result atomically."
)
log_message("Completed formal two-stage R3 production cache: ", final_dir)
