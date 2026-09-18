#!/usr/bin/env Rscript

# Reconstruct, fit, and evaluate the internal two-stage R3 experiment on Midway3.
get_arg <- function(name, default) {
  args <- commandArgs(trailingOnly = TRUE)
  index <- match(name, args)
  if (is.na(index)) return(default)
  if (index == length(args)) stop("Missing value for ", name, call. = FALSE)
  args[[index + 1L]]
}
root <- normalizePath(get_arg("--project-root", "."), mustWork = TRUE)
experiment <- Sys.getenv("FASH_R3_TWO_STAGE_CODE", unset = file.path(
  root, "code/revision_simulations/internal/r3_two_stage_functional_testing"
))
source(file.path(experiment, "two_stage_helpers.R"))
source(file.path(experiment, "reconstruct_r3.R"))
contract <- r3ts_contract()
mode <- get_arg("--mode", "run")
r3ts_assert(mode %in% c("preflight", "run"), "Mode must be preflight or run.")
cores <- as.integer(get_arg("--num-cores", Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
r3ts_assert(length(cores) == 1L && is.finite(cores) && cores >= 1L, "Invalid worker count.")
if (mode == "run" && identical(Sys.info()[["sysname"]], "Darwin")) {
  stop("Full R3 production must run on Midway3. Use the separate reduced smoke script locally.", call. = FALSE)
}
log_message <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")
  flush.console()
}
sources <- c(
  simulation = Sys.getenv("FASH_R3_SIMULATION_FUNCTIONS", file.path(root,
    "code/revision_simulations/r3_r4_fashr0143/source_snapshots/r3_full_universe_functional_simulation_functions.R")),
  genotype_helper = Sys.getenv("FASH_R3_REAL_GENOTYPE_HELPER", file.path(root,
    "code/revision_simulations/shared/real_genotype_one_per_gene.R")),
  digest_helper = Sys.getenv("FASH_R3_DIGEST_HELPER", file.path(root,
    "code/revision_simulations/internal/r3_ideal_gaussian_measurement/ideal_gaussian_measurement.R")),
  genotypes = Sys.getenv("FASH_R3_GENOTYPE_CACHE", file.path(root,
    "output/revision_simulations/shared/real_genotype_one_per_gene_J6362_pilot5/genotype_samples.rds"))
)
source_hashes <- vapply(sources, r3ts_sha256, character(1))
r3ts_assert(identical(source_hashes, contract$source_sha), "Frozen R3 source/input hash mismatch.")
for (name in c("simulation", "genotype_helper", "digest_helper")) source(sources[[name]])
description <- utils::packageDescription("fashr")
r3ts_assert(identical(description$Version, contract$package_version) &&
               identical(description$RemoteSha, contract$package_sha), "Unexpected fashr build.")
posterior_draws <- formals(getS3method("predict", "fash", envir = asNamespace("fashr")))$M
r3ts_assert(isTRUE(all.equal(posterior_draws, contract$posterior_draws)),
             "Posterior draw count differs from the frozen R3 sampler.")
formal_dir <- Sys.getenv("FASH_R3_FORMAL_DIR", file.path(root,
  "output/revision_simulations/mc", contract$formal_id))
ideal_dir <- Sys.getenv("FASH_R3_IDEAL_DIR", file.path(root,
  "output/revision_simulations/internal", contract$ideal_id))
replicate_names <- as.vector(outer(contract$mechanisms, contract$seeds,
                                   function(m, s) paste0(m, "_seed_", s, ".rds")))
check_reference_cache <- function(directory, id, artifacts) {
  flag <- readLines(file.path(directory, "complete.flag"), warn = FALSE)
  manifest <- readRDS(file.path(directory, "manifest.rds"))
  r3ts_assert(identical(manifest$result_id, id) && paste0("result_id=", id) %in% flag,
               "Incorrect or incomplete reference cache.")
  hashes <- manifest$artifact_sha256[artifacts]
  r3ts_assert(!anyNA(hashes) && identical(unname(hashes),
    unname(vapply(file.path(directory, artifacts), r3ts_sha256, character(1)))),
    "Reference cache artifact hash mismatch.")
  manifest
}
formal_manifest <- check_reference_cache(formal_dir, contract$formal_id,
  c("configuration.rds", file.path("replicates", replicate_names)))
ideal_manifest <- check_reference_cache(ideal_dir, contract$ideal_id,
  c("configuration.rds", "summary/replicate_input_digests.csv"))
config <- readRDS(file.path(formal_dir, "configuration.rds"))
r3ts_validate_configuration(config)
reference_digests <- read.csv(file.path(ideal_dir, "summary/replicate_input_digests.csv"))
reference_keys <- paste(reference_digests$truth_mechanism, reference_digests$seed, sep = ":")
r3ts_assert(!anyDuplicated(reference_keys) && setequal(reference_keys, r3ts_expected_keys()),
             "The reference digest table lacks one or more formal replicates.")
genotypes <- readRDS(sources[["genotypes"]])
for (seed in contract$seeds) {
  g <- validate_real_genotype_sample(genotypes$samples[[as.character(seed)]],
    expected_genes = config$J, expected_donors = config$n_donors,
    maf_min = genotypes$configuration$maf_min)
  r3ts_assert(identical(genotype_content_md5(g$selection$pair_key, rownames(g$G), g$G),
                        config$genotype_content_digests[[as.character(seed)]]),
               "Unexpected genotype selection or dosage matrix.")
}
runner_path <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])
own_sources <- c(helpers = file.path(experiment, "two_stage_helpers.R"),
                 reconstruction = file.path(experiment, "reconstruct_r3.R"), runner = runner_path)
own_hashes <- vapply(own_sources, r3ts_sha256, character(1))
package <- list(version = description$Version, remote_sha = description$RemoteSha,
                r_version = R.version.string, platform = R.version$platform)
configuration <- list(
  contract = contract, formal_configuration = config, package = package,
  frozen_source_sha256 = source_hashes, experiment_source_sha256 = own_hashes,
  reference_manifest_sha256 = c(formal = r3ts_sha256(file.path(formal_dir, "manifest.rds")),
                                ideal = r3ts_sha256(file.path(ideal_dir, "manifest.rds"))),
  fit_universe = "all_units", prior_refit_after_screening = FALSE,
  functional_candidate_scope = "dynamic_fdr_screen", num_cores = cores,
  posterior_draws = posterior_draws, q_dyn = contract$q_dyn, q_func = contract$q_func
)
token <- serialized_object_md5(configuration)
log_message("Preflight passed: frozen inputs, package, both reference caches, ten replicate keys, q_dyn=q_func=0.05.")
log_message("Posterior sampler: frozen testing_functional/predict.fash; M=", posterior_draws, "; workers=", cores)
if (mode == "preflight") quit(save = "no", status = 0L)

parent <- Sys.getenv("FASH_R3_TWO_STAGE_RESULT_PARENT", file.path(root,
  "output/revision_simulations/internal"))
final <- file.path(parent, contract$result_id)
partial <- paste0(final, ".partial")
checkpoint_root <- Sys.getenv("FASH_R3_TWO_STAGE_CHECKPOINT_ROOT", "")
r3ts_assert(nzchar(checkpoint_root), "Set a server-only FASH_R3_TWO_STAGE_CHECKPOINT_ROOT.")
r3ts_assert(!dir.exists(final), "Refusing to overwrite an existing final result directory.")
dir.create(partial, recursive = TRUE, showWarnings = FALSE)
config_path <- file.path(partial, "configuration.rds")
if (file.exists(config_path)) {
  r3ts_assert(identical(readRDS(config_path), configuration),
               "Partial output belongs to a different configuration or source revision.")
} else r3ts_write_rds(configuration, config_path)
checkpoint_token <- serialized_object_md5(list(
  formal_configuration = config, source_sha256 = source_hashes,
  reconstruction_sha256 = own_hashes[["reconstruction"]], package = package, num_cores = cores
))

for (mechanism in contract$mechanisms) {
  for (seed in contract$seeds) {
    key <- paste0(mechanism, "_seed_", seed)
    result_path <- file.path(partial, "replicates", paste0(key, ".rds"))
    if (file.exists(result_path)) {
      saved <- readRDS(result_path)
      r3ts_assert(identical(saved$configuration_token, token), "Incompatible completed replicate.")
      r3ts_validate_replicate(saved)
      log_message("Reusing completed internal replicate: ", key)
      next
    }
    started <- proc.time()[["elapsed"]]
    input_path <- file.path(checkpoint_root, key, "inputs.rds")
    fit_path <- file.path(checkpoint_root, key, "fits.rds")
    formal_replicate <- readRDS(file.path(formal_dir, "replicates", paste0(key, ".rds")))
    reused_inputs <- file.exists(input_path)
    if (reused_inputs) {
      saved <- readRDS(input_path)
      r3ts_assert(identical(saved$token, checkpoint_token), "Incompatible input checkpoint.")
      inputs <- saved$inputs
    } else {
      log_message("Reconstructing frozen truth and observations: ", key)
      inputs <- r3ts_reconstruct(config, genotypes$samples[[as.character(seed)]], seed, mechanism)
      r3ts_verify_reconstruction(inputs, reference_digests, config, formal_replicate)
      r3ts_write_rds(list(token = checkpoint_token, inputs = inputs), input_path)
    }
    r3ts_verify_reconstruction(inputs, reference_digests, config, formal_replicate)
    log_message("Four original-input digests and retained truth checks passed: ", key)
    reused_fit <- file.exists(fit_path)
    if (reused_fit) {
      saved <- readRDS(fit_path)
      r3ts_assert(identical(saved$token, checkpoint_token), "Incompatible fit checkpoint.")
      fits <- saved$fits
    } else {
      log_message("Fitting IWP1 to all 6362 units and applying BF update: ", key)
      fits <- r3ts_fit_all(inputs, config, cores)
      r3ts_write_rds(list(token = checkpoint_token, fits = fits), fit_path)
    }
    fit <- fits$fash_iwp1_bf
    old_pi0 <- formal_replicate$estimated_pi0$estimated_pi0[
      formal_replicate$estimated_pi0$fit == "BF-corrected"]
    new_pi0 <- constant_component_prior_weight(fit)
    r3ts_assert(length(old_pi0) == 1L && abs(new_pi0 - old_pi0) < 1e-9,
                 "Recomputed BF prior weight differs from formal R3.")
    log_message("Applying two-stage selection: ", key)
    inference <- r3ts_infer(fit, inputs, config, cores)
    old_calls <- unique(formal_replicate$functional_alpha_005$dynamic_discoveries[
      formal_replicate$functional_alpha_005$method == "FASH-IWP1-BF"])
    r3ts_assert(length(old_calls) == 1L && length(inference$candidates) == old_calls,
                 "Stage-1 discoveries differ from the retained full-data R3 fit.")
    result <- list(
      schema = contract$schema, configuration_token = token,
      seed = seed, truth_mechanism = mechanism, selected_pair_keys = formal_replicate$selected_pair_keys,
      input_digests = inputs$input_digests, inference = inference,
      fitted_n_units = length(fit$fash_data$data_list), fitted_bf_pi0 = new_pi0,
      original_dynamic_calls = old_calls, original_bf_pi0 = old_pi0,
      reuse = c(genotype_cache = TRUE, source_configuration = TRUE,
                input_checkpoint = reused_inputs, fit_checkpoint = reused_fit),
      checkpoint_paths = c(inputs = input_path, fits = fit_path),
      elapsed_seconds = proc.time()[["elapsed"]] - started
    )
    r3ts_validate_replicate(result)
    r3ts_write_rds(result, result_path)
    log_message("Completed ", key, ": Stage-1 calls=", length(inference$candidates),
                 "; functional calls=", paste(inference$metrics$functional$calls, collapse = ","))
    rm(list = intersect(c("inputs", "fits", "fit", "inference", "result", "saved"), ls()))
    invisible(gc())
  }
}

replicates <- lapply(file.path(partial, "replicates", replicate_names), readRDS)
invisible(lapply(replicates, r3ts_validate_replicate))
stage1 <- do.call(rbind, lapply(replicates, function(x) x$inference$metrics$stage1))
functional <- do.call(rbind, lapply(replicates, function(x) x$inference$metrics$functional))
summaries <- r3ts_summaries(stage1, functional)
summary_dir <- file.path(partial, "summary")
dir.create(summary_dir, showWarnings = FALSE)
tables <- list(stage1_by_seed = stage1, functional_by_seed = functional,
               stage1_mc_summary = summaries$stage1, functional_mc_summary = summaries$primary)
for (mechanism in contract$mechanisms) {
  tables[[paste0("primary_", mechanism)]] <- r3ts_primary_values(summaries$primary, mechanism)
}
tables$input_validation <- do.call(rbind, lapply(replicates, function(x) {
  data.frame(seed = x$seed, truth_mechanism = x$truth_mechanism,
    as.list(x$input_digests), fitted_n_units = x$fitted_n_units,
    original_bf_pi0 = x$original_bf_pi0, fitted_bf_pi0 = x$fitted_bf_pi0,
    original_dynamic_calls = x$original_dynamic_calls,
    stage1_calls = length(x$inference$candidates), stringsAsFactors = FALSE)
}))
tables$all_stage1_units <- do.call(rbind, lapply(replicates, function(x) {
  fdr <- x$inference$fdr_table
  fdr <- fdr[match(seq_len(contract$n_units), fdr$index), , drop = FALSE]
  data.frame(seed = x$seed, truth_mechanism = x$truth_mechanism,
    pair_key = x$selected_pair_keys, unit_index = fdr$index,
    lfdr = fdr$lfdr, cumulative_fdr = fdr$FDR,
    selected_dynamic = fdr$index %in% x$inference$candidates,
    true_dynamic = x$inference$true_dynamic, x$inference$true_functionals)
}))
tables$all_stage2_candidates <- do.call(rbind, lapply(replicates, function(x) {
  tab <- x$inference$stage2
  tab$seed <- rep(x$seed, nrow(tab))
  tab$truth_mechanism <- rep(x$truth_mechanism, nrow(tab))
  tab$pair_key <- x$selected_pair_keys[tab$index]
  tab$true_functional <- x$inference$true_functionals[
    cbind(tab$index, match(tab$target, contract$targets))]
  tab
}))
for (name in names(tables)) {
  write.csv(tables[[name]], file.path(summary_dir, paste0(name, ".csv")), row.names = FALSE)
}
r3ts_write_rds(tables, file.path(partial, "summary_tables.rds"))
writeLines(capture.output(sessionInfo()), file.path(partial, "sessionInfo.txt"))
artifacts <- list.files(partial, recursive = TRUE, full.names = FALSE)
hashes <- setNames(vapply(file.path(partial, artifacts), r3ts_sha256, character(1)), artifacts)
manifest <- list(schema = contract$schema, result_id = contract$result_id,
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  configuration = configuration, configuration_token = token,
  replicate_keys = r3ts_expected_keys(), artifact_sha256 = hashes)
r3ts_write_rds(manifest, file.path(partial, "manifest.rds"))
writeLines(c(paste0("result_id=", contract$result_id), "replicates=10",
             "q_dyn=0.05", "q_func=0.05", "input_digest_checks=passed",
             "fit_universe=6362", "functional_candidate_scope=dynamic_fdr_screen"),
           file.path(partial, "complete.flag"))
r3ts_assert(file.rename(partial, final), "Failed to promote the completed cache atomically.")
log_message("Completed internal R3 two-stage experiment: ", final)
