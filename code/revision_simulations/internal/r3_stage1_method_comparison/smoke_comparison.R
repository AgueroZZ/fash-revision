#!/usr/bin/env Rscript

# Reduced integration only: one R process, at most two workers, no scientific claims.
experiment <- "code/revision_simulations/internal/r3_stage1_method_comparison"
source(file.path(experiment, "comparison_helpers.R"))
source(file.path(experiment, "reporting.R"))
base <- r3c_base(".")
context <- r3c_context(".", base)
genotype <- readRDS(base$genotype_path)$samples[["12345"]]
keep <- seq_len(60L)
genotype$G <- genotype$G[, keep, drop = FALSE]
genotype$selection <- genotype$selection[keep, , drop = FALSE]
genotype$variant_info <- genotype$variant_info[keep, , drop = FALSE]
config <- base$config$formal_configuration
output <- "output/revision_simulations/internal/r3_stage1_method_comparison_development_20260917"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
rows <- list()
for (mechanism in r3c_contract()$mechanisms) {
  started <- proc.time()[["elapsed"]]
  inputs <- context$sim$r3ts_reconstruct(config, genotype, 12345L, mechanism)
  datasets <- context$sim$make_fash_datasets_from_eqtl_summary(inputs$regression$beta_hat,
    inputs$regression$se, inputs$effects$beta_matrix, config$time_grid,
    inputs$effects$unit_info, inputs$scenario)
  fits <- context$sim$fit_fash_for_revision(datasets, orders = 1,
    num_cores = 2L, num_basis = config$num_basis, penalty = 10, pred_step = 1,
    apply_bf = TRUE, verbose = FALSE)
  truth <- inputs$effects$unit_info$effect_class == "dynamic_bspline"
  reference <- list(seed = 12345L, truth_mechanism = mechanism,
    input_digests = inputs$input_digests, selected_pair_keys = genotype$selection$pair_key,
    inference = list(true_dynamic = truth, true_functionals = inputs$effects$true_functionals,
      fdr_table = context$sim$get_fash_fdr_table(fits$fash_iwp1_bf)))
  r3c_validate_inputs(inputs, genotype, reference, config, context$sim)
  broken <- inputs
  broken$expression$expression[1, 1, 1] <- broken$expression$expression[1, 1, 1] + 1
  stopifnot(inherits(try(r3c_validate_inputs(broken, genotype, reference, config, context$sim),
                           silent = TRUE), "try-error"))
  linear <- r3c_fit_linear(inputs, config, context$cmp)
  payload <- r3c_direct_payload(inputs, genotype, config)
  permutation <- context$cmp$compute_direct_interaction_permutation_null(payload,
    n_permutations = 3L, seed = 22345L, num_cores = 2L, verbose = FALSE)
  serial <- context$cmp$compute_direct_interaction_permutation_null(payload,
    n_permutations = 3L, seed = 22345L, num_cores = 1L, verbose = FALSE)
  stopifnot(identical(permutation$permutation_index, serial$permutation_index),
            identical(permutation$null_pvalues, serial$null_pvalues))
  direct <- r3c_direct_scores(payload, permutation, context$cmp, true_pi0 = mean(!truth))
  scores <- c(list(
    "FASH-IWP1-Raw" = r3c_scores_from_fdr(context$sim$get_fash_fdr_table(fits$fash_iwp1_raw), 60L),
    "FASH-IWP1-BF" = r3c_scores_from_fdr(reference$inference$fdr_table, 60L)),
    r3c_linear_scores(linear, inputs, context$cmp), direct$scores)
  record <- list(schema = r3c_contract()$schema, seed = 12345L, truth_mechanism = mechanism,
    input_digests = inputs$input_digests, selected_pair_keys = reference$selected_pair_keys,
    true_dynamic = truth, scores = scores, direct_diagnostics = direct$diagnostics,
    permutation_index = permutation$permutation_index, permutation_settings = permutation$settings,
    prior_weights = list("FASH-IWP1-Raw" = fits$fash_iwp1_raw$prior_weights,
      "FASH-IWP1-BF" = fits$fash_iwp1_bf$prior_weights,
      "FASH-linear-Raw" = linear$raw$prior_weights, "FASH-linear-BF" = linear$bf$prior_weights))
  smoke_contract <- r3c_contract()
  smoke_contract$n_units <- 60L
  smoke_contract$n_permutations <- 3L
  smoke_contract$true_pi0 <- mean(!truth)
  r3c_validate_record(record, reference, smoke_contract)
  curve <- r3c_curves(scores, truth, 12345L, mechanism)
  r3c_write(record, file.path(output, paste0("smoke_", mechanism, ".rds")))
  rows[[mechanism]] <- data.frame(mechanism = mechanism, units = 60L, permutations = 3L,
    retained_methods = length(scores), curve_rows = nrow(curve),
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    scope = "Reduced integration smoke only; excluded from production reporting")
  rm(inputs, datasets, fits, linear, payload, permutation, serial, direct, record, broken)
  invisible(gc())
}
utils::write.csv(do.call(rbind, rows), file.path(output, "smoke_validation.csv"), row.names = FALSE)
cat("Both 60-unit truth mechanisms passed Raw/BF, direct LRT/eFDR, exact input pairing, and serial/parallel permutation checks.\n")
