#!/usr/bin/env Rscript

# Extend the frozen R3 replicates without rerunning IWP1 or Stage 2.
argument <- function(name, default) {
  args <- commandArgs(trailingOnly = TRUE)
  index <- match(name, args)
  if (is.na(index)) return(default)
  if (index == length(args)) stop("Missing value for ", name, call. = FALSE)
  args[[index + 1L]]
}
log_message <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ..., "\n", sep = "")
  flush.console()
}
root <- normalizePath(argument("--project-root", "."), mustWork = TRUE)
experiment <- file.path(root, "code/revision_simulations/internal/r3_stage1_method_comparison")
source(file.path(experiment, "comparison_helpers.R"))
contract <- r3c_contract()
mode <- argument("--mode", "preflight")
cores <- as.integer(argument("--num-cores", Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
r3c_assert(mode %in% c("preflight", "run") && is.finite(cores) && cores >= 1L,
           "Use preflight or run and a positive worker count.")
if (mode == "run") {
  r3c_assert(identical(Sys.info()[["sysname"]], "Linux") && nzchar(Sys.getenv("SLURM_JOB_ID")),
             "Full production must run inside a Midway3 Slurm job.")
  r3c_assert(identical(as.character(getRversion()), "4.4.1"),
             "R 4.4.1 is required to reproduce the retained RDS input digests.")
}

base <- r3c_base(root)
context <- r3c_context(root, base)
genotypes <- readRDS(base$genotype_path)
if (identical(Sys.info()[["sysname"]], "Linux") && nzchar(Sys.getenv("SLURM_JOB_ID"))) {
  fit_paths <- vapply(base$records, function(x) x$checkpoint_paths[["fits"]], character(1))
  r3c_assert(all(file.exists(fit_paths)), paste(
    "Original IWP1 fit checkpoints are required for Raw-result retention. Missing:",
    paste(fit_paths[!file.exists(fit_paths)], collapse = ", ")))
}
runner <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]]))
own_paths <- c(helpers = file.path(experiment, "comparison_helpers.R"), runner = runner)
configuration <- list(
  contract = contract, baseline_manifest_sha256 = r3c_sha(file.path(base$directory, "manifest.rds")),
  baseline_configuration_sha256 = r3c_sha(file.path(base$directory, "configuration.rds")),
  genotype_sha256 = r3c_sha(base$genotype_path),
  scientific_source_sha256 = context$hashes,
  experiment_source_sha256 = vapply(own_paths, r3c_sha, character(1)),
  package_version = contract$package_version, package_sha = contract$package_sha,
  r_version = as.character(getRversion()), num_cores = cores,
  stage2 = "Reuse the original BF-IWP1 Stage-1 candidate sets and Stage-2 results unchanged."
)
token <- digest::digest(configuration, algo = "sha256")
log_message("Preflight passed: 20 R3 records, ten-seed genotypes, pinned sources and fashr build.")
log_message("Retaining six method results; displaying four. Direct tests: 100 paired donor permutations.")
if (mode == "preflight") quit(save = "no", status = 0L)

checkpoint_root <- Sys.getenv("FASH_R3_COMPARE_CHECKPOINT_ROOT", "")
result_parent <- Sys.getenv("FASH_R3_COMPARE_RESULT_PARENT", "")
r3c_assert(startsWith(checkpoint_root, "/project/mstephens/ziangzhang/fash/full_results/") &&
             startsWith(result_parent, "/project/mstephens/ziangzhang/fash/workspace/results/"),
           "Set the documented server-only checkpoint and compact-result directories.")
final_dir <- file.path(result_parent, contract$result_id)
partial_dir <- paste0(final_dir, ".partial")
r3c_assert(!dir.exists(final_dir), "Refusing to overwrite a completed comparison result.")
dir.create(partial_dir, recursive = TRUE, showWarnings = FALSE)
lock <- file.path(partial_dir, ".running")
r3c_assert(dir.create(lock, showWarnings = FALSE),
           "This result is locked. Verify the previous job has ended before removing its .running lock.")
# A normal error removes the lock; scheduler termination may leave it for manual review.
with_lock <- function(work) {
  on.exit(if (dir.exists(lock)) system2("rmdir", shQuote(lock)), add = TRUE)
  work()
}

with_lock(function() {
  configuration_path <- file.path(partial_dir, "configuration.rds")
  if (file.exists(configuration_path)) {
    r3c_assert(identical(readRDS(configuration_path), configuration),
               "The partial result uses a different configuration or source revision.")
  } else r3c_write(configuration, configuration_path)

  checkpoint <- function(path, compute) {
    if (file.exists(path)) {
      saved <- readRDS(path)
      r3c_assert(identical(saved$token, token), paste("Incompatible checkpoint:", path))
      return(saved$value)
    }
    value <- compute()
    r3c_write(list(token = token, value = value), path)
    value
  }

  for (key in names(base$records)) {
    reference <- base$records[[key]]
    stem <- paste0(reference$truth_mechanism, "_seed_", reference$seed)
    result_path <- file.path(partial_dir, "replicates", paste0(stem, ".rds"))
    if (file.exists(result_path)) {
      saved <- readRDS(result_path)
      r3c_assert(identical(saved$configuration_token, token), "Incompatible completed replicate.")
      r3c_validate_record(saved, reference)
      log_message("Reusing completed comparison: ", stem)
      next
    }
    started <- proc.time()[["elapsed"]]
    current_checkpoint <- file.path(checkpoint_root, stem)
    genotype <- genotypes$samples[[as.character(reference$seed)]]
    input_path <- reference$checkpoint_paths[["inputs"]]
    if (file.exists(input_path)) {
      inputs <- readRDS(input_path)$inputs
      input_origin <- "original R3 input checkpoint"
    } else {
      input_path <- file.path(current_checkpoint, "reconstructed_inputs.rds")
      inputs <- checkpoint(input_path, function() context$sim$r3ts_reconstruct(
        base$config$formal_configuration, genotype, reference$seed, reference$truth_mechanism))
      input_origin <- "reconstructed and matched to retained R3 input digests"
    }
    log_message("Validating exact genotype, truth, summaries, covariates and expression: ", stem)
    input_digests <- r3c_validate_inputs(inputs, genotype, reference,
                                       base$config$formal_configuration, context$sim)

    iwp <- checkpoint(file.path(current_checkpoint, "iwp_scores.rds"), function() {
      path <- reference$checkpoint_paths[["fits"]]
      r3c_assert(file.exists(path), paste("Original full IWP1 fits are required to retain Raw results:", path))
      fits <- readRDS(path)$fits
      r3c_assert(all(c("fash_iwp1_raw", "fash_iwp1_bf") %in% names(fits)),
                 "The original checkpoint lacks Raw or BF IWP1 fits.")
      raw <- r3c_scores_from_fdr(context$sim$get_fash_fdr_table(fits$fash_iwp1_raw), contract$n_units)
      bf <- r3c_scores_from_fdr(context$sim$get_fash_fdr_table(fits$fash_iwp1_bf), contract$n_units)
      expected <- r3c_scores_from_fdr(reference$inference$fdr_table, contract$n_units)
      r3c_assert(isTRUE(all.equal(bf, expected, tolerance = 0)),
                 "The IWP1 full-fit checkpoint differs from the completed R3 BF result.")
      list(scores = list("FASH-IWP1-Raw" = raw, "FASH-IWP1-BF" = bf),
           prior_weights = list("FASH-IWP1-Raw" = fits$fash_iwp1_raw$prior_weights,
                                "FASH-IWP1-BF" = fits$fash_iwp1_bf$prior_weights),
           source_fit_path = path, source_fit_sha256 = r3c_sha(path))
    })
    invisible(gc())
    log_message("Fitting or reusing matched FASH-linear Raw and BF: ", stem)
    linear_path <- file.path(current_checkpoint, "linear_fits.rds")
    linear <- checkpoint(linear_path, function() r3c_fit_linear(
      inputs, base$config$formal_configuration, context$cmp))
    linear_scores <- r3c_linear_scores(linear, inputs, context$cmp)
    prior_weights <- c(iwp$prior_weights, list(
      "FASH-linear-Raw" = linear$raw$prior_weights,
      "FASH-linear-BF" = linear$bf$prior_weights))
    linear_prior_summary <- do.call(rbind, lapply(c("raw", "bf"), function(label) {
      context$cmp$summarize_linear_mixture_prior_fit(linear[[label]], reference$seed,
        if (label == "raw") "FASH-linear-Raw" else "FASH-linear-BF")
    }))
    rm(linear)
    invisible(gc())
    payload <- r3c_direct_payload(inputs, genotype, base$config$formal_configuration)
    permutation_path <- file.path(current_checkpoint, "direct_permutations.rds")
    log_message("Computing or reusing 100 donor permutations for both direct tests: ", stem)
    permutation <- checkpoint(permutation_path, function() {
      context$cmp$compute_direct_interaction_permutation_null(
        out = payload, n_permutations = contract$n_permutations, interaction_degrees = c(1, 2),
        seed = reference$seed + contract$permutation_seed_offset,
        permute_covariates_with_expression = TRUE, num_cores = cores, verbose = TRUE)
    })
    r3c_assert(all(vapply(permutation$null_pvalues,
      function(x) identical(dim(x), c(contract$n_permutations, contract$n_units)) &&
        all(is.finite(x)) && all(x >= 0 & x <= 1), logical(1))), "Incomplete permutation null.")
    direct <- r3c_direct_scores(payload, permutation, context$cmp)
    scores <- c(iwp$scores, linear_scores, direct$scores)[contract$methods]
    record <- list(schema = contract$schema, configuration_token = token,
      seed = reference$seed, truth_mechanism = reference$truth_mechanism,
      selected_pair_keys = reference$selected_pair_keys, input_digests = input_digests,
      true_dynamic = reference$inference$true_dynamic, scores = scores,
      prior_weights = prior_weights, linear_prior_summary = linear_prior_summary,
      direct_diagnostics = direct$diagnostics, permutation_index = permutation$permutation_index,
      permutation_settings = permutation$settings, input_origin = input_origin,
      checkpoint_paths = c(inputs = input_path, iwp_fits = iwp$source_fit_path,
                           linear_fits = linear_path, direct_permutations = permutation_path),
      checkpoint_sha256 = c(inputs = r3c_sha(input_path), iwp_fits = iwp$source_fit_sha256,
        linear_fits = r3c_sha(linear_path), direct_permutations = r3c_sha(permutation_path)),
      elapsed_seconds = proc.time()[["elapsed"]] - started)
    r3c_validate_record(record, reference)
    r3c_write(record, result_path)
    log_message("Completed ", stem, "; retained Raw/BF and all direct scores; seconds=",
                 round(record$elapsed_seconds, 1))
    rm(inputs, iwp, payload, permutation, direct, scores, record, linear_scores)
    invisible(gc())
  }

  records <- lapply(base$records, function(reference) {
    path <- file.path(partial_dir, "replicates",
      paste0(reference$truth_mechanism, "_seed_", reference$seed, ".rds"))
    record <- readRDS(path)
    r3c_validate_record(record, reference)
    record
  })
  curves <- do.call(rbind, lapply(records, function(x) r3c_curves(
    x$scores, x$true_dynamic, x$seed, x$truth_mechanism)))
  rownames(curves) <- NULL
  summary <- r3c_summarize(curves)
  tables <- list(stage1_alpha_by_seed = curves, stage1_mc_summary = summary,
    stage1_alpha005_summary = summary[abs(summary$alpha - 0.05) < 1e-12, , drop = FALSE],
    functional_by_seed = base$tables$functional_by_seed,
    functional_mc_summary = base$tables$functional_mc_summary,
    stage2_random_bspline = base$tables$stage2_random_bspline,
    stage2_raised_cosine = base$tables$stage2_raised_cosine)
  r3c_write(tables, file.path(partial_dir, "summary_tables.rds"))
  dir.create(file.path(partial_dir, "summary"), showWarnings = FALSE)
  for (name in names(tables)) utils::write.csv(tables[[name]],
    file.path(partial_dir, "summary", paste0(name, ".csv")), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(partial_dir, "sessionInfo.txt"))
  files <- c("configuration.rds", "summary_tables.rds", "sessionInfo.txt",
    file.path("summary", paste0(names(tables), ".csv")),
    file.path("replicates", paste0(vapply(records, function(x)
      paste0(x$truth_mechanism, "_seed_", x$seed), character(1)), ".rds")))
  hashes <- setNames(vapply(file.path(partial_dir, files), r3c_sha, character(1)), files)
  r3c_write(list(schema = contract$schema, result_id = contract$result_id,
    configuration = configuration, configuration_token = token,
    artifact_sha256 = hashes, generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    file.path(partial_dir, "manifest.rds"))
  writeLines(c(paste0("result_id=", contract$result_id), "replicates=20", "seeds=10",
    "retained_methods=6", "display_methods=4", "stage1_alpha_rows=4800",
    "raw_and_bf_retained=IWP1,FASH-linear", "stage2=unchanged", "donor_permutations=100"),
    file.path(partial_dir, "complete.flag"))
  r3c_manifest(partial_dir)
  # Remove only this job's empty lock before atomically promoting the result.
  r3c_assert(system2("rmdir", shQuote(lock)) == 0L, "Could not remove the empty job lock.")
  r3c_assert(file.rename(partial_dir, final_dir), "Could not finalize the completed comparison cache.")
  log_message("Completed paired R3 comparison cache: ", final_dir)
})
