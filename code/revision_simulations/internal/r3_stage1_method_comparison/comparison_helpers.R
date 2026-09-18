# Paired R3 comparisons using the frozen R1/R2 comparator implementations.

r3c_contract <- function() {
  list(
    schema = "r3-stage1-method-comparison-v1",
    result_id = "r3_stage1_iwp1_linear_direct_paired_raw_bf_fashr0143_pilot10",
    base_id = paste0("r3_two_stage_functional_testing_iwp1_mixture_alpha_grid_",
                     "qdyn005_qfunc005_fashr0143_pilot10"),
    base_manifest_sha256 = "16f93e24ba3dec3f459962d031a242a0f30539708e53320638acbeb0bddbbf97",
    comparator_sha256 = "93f9a2c5606ae74763fb53e189af62c4d3b7c973aaf38d5b3eb3bf32f97a6487",
    seeds = seq(12345L, 102345L, by = 10000L),
    mechanisms = c("random_bspline", "raised_cosine"),
    n_units = 6362L, n_dynamic = 1272L,
    alpha_grid = seq(0.005, 0.20, by = 0.005),
    methods = c("FASH-IWP1-Raw", "FASH-IWP1-BF", "FASH-linear-Raw",
                "FASH-linear-BF", "Direct-linear-LRT-eFDR-true-pi0",
                "Direct-quadratic-LRT-eFDR-true-pi0"),
    display_methods = c("FASH-IWP1-BF", "FASH-linear-BF",
                        "Direct-linear-LRT-eFDR-true-pi0",
                        "Direct-quadratic-LRT-eFDR-true-pi0"),
    n_permutations = 100L, permutation_seed_offset = 10000L,
    true_pi0 = 5090 / 6362, pred_step = 1, penalty = 10,
    grid = sort(c(0, exp(-0.5 * seq(0, 10, by = 0.2)))),
    package_version = "0.1.43",
    package_sha = "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
  )
}

r3c_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
  invisible(TRUE)
}

r3c_sha <- function(path) {
  r3c_assert(file.exists(path) && !dir.exists(path), paste("Missing file:", path))
  digest::digest(file = path, algo = "sha256")
}

r3c_write <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3)
  r3c_assert(file.rename(temporary, path), paste("Atomic write failed:", path))
}

r3c_manifest <- function(directory, expected_sha = NULL) {
  path <- file.path(directory, "manifest.rds")
  r3c_assert(file.exists(file.path(directory, "complete.flag")),
             paste("Completed cache required:", directory))
  if (!is.null(expected_sha)) {
    r3c_assert(identical(r3c_sha(path), expected_sha), "Baseline manifest changed.")
  }
  manifest <- readRDS(path)
  hashes <- manifest$artifact_sha256
  r3c_assert(is.character(hashes) && length(hashes) > 0L &&
               !anyDuplicated(names(hashes)), "Invalid artifact manifest.")
  actual <- vapply(file.path(directory, names(hashes)), r3c_sha, character(1))
  r3c_assert(identical(unname(actual), unname(hashes)),
             paste("Cache artifact hashes changed:", directory))
  manifest
}

r3c_base <- function(root, directory = Sys.getenv("FASH_R3_COMPARE_BASE", ""),
                     genotype_path = Sys.getenv("FASH_R3_COMPARE_GENOTYPES", "")) {
  contract <- r3c_contract()
  if (!nzchar(directory)) directory <- file.path(root, "output/revision_simulations/mc", contract$base_id)
  manifest <- r3c_manifest(directory, contract$base_manifest_sha256)
  config <- readRDS(file.path(directory, "configuration.rds"))
  r3c_assert(identical(config, manifest$configuration), "Baseline configuration changed.")
  paths <- c(
    pilot10 = file.path(root, "code/revision_simulations/internal/r3_two_stage_functional_testing_pilot10/pilot10_helpers.R"),
    two_stage = file.path(root, "code/revision_simulations/internal/r3_two_stage_functional_testing/two_stage_helpers.R")
  )
  expected <- c(pilot10 = config$experiment_source_sha256[["pilot10_helpers"]],
                two_stage = config$contract$base_source_sha[["base_two_stage_helper"]])
  r3c_assert(identical(vapply(paths, r3c_sha, character(1)), expected),
             "Baseline validation helpers changed.")
  env <- new.env(parent = globalenv())
  for (path in paths) sys.source(path, envir = env)
  r3c_assert(identical(config$contract, env$r3ts10_contract()) &&
               identical(manifest$result_id, contract$base_id), "Incorrect baseline contract.")
  files <- file.path(directory, "replicates", env$r3ts10_replicate_names())
  records <- lapply(files, readRDS)
  for (record in records) {
    r3c_assert(identical(record$configuration_token, manifest$configuration_token),
               "Baseline record has a different configuration token.")
    env$r3ts10_validate_record(record, env, env$r3ts10_functional_cfsr_table)
  }
  names(records) <- vapply(records, function(x) paste(x$truth_mechanism, x$seed, sep = ":"), character(1))
  tables <- readRDS(file.path(directory, "summary_tables.rds"))
  curves <- do.call(rbind, lapply(records, env$r3ts10_stage1_alpha_by_seed, base_environment = env))
  rownames(curves) <- NULL
  r3c_assert(isTRUE(all.equal(curves, tables$stage1_alpha_by_seed, check.attributes = FALSE)),
             "Baseline Stage-1 curves do not reproduce.")
  if (!nzchar(genotype_path)) genotype_path <- file.path(
    root, "output/revision_simulations/shared", config$contract$genotype_id, "genotype_samples.rds")
  r3c_assert(identical(r3c_sha(genotype_path), config$genotype_cache_sha256),
             "Genotypes differ from the completed R3 experiment.")
  list(directory = directory, manifest = manifest, config = config, env = env,
       records = records, tables = tables, genotype_path = genotype_path)
}

r3c_context <- function(root, base) {
  local_path <- function(variable, relative) Sys.getenv(variable, file.path(root, relative))
  paths <- c(
    simulation = local_path("FASH_R3_COMPARE_SIMULATION",
      "code/revision_simulations/r3_r4_fashr0143/source_snapshots/r3_full_universe_functional_simulation_functions.R"),
    genotype_helper = local_path("FASH_R3_COMPARE_GENOTYPE_HELPER",
      "code/revision_simulations/shared/real_genotype_one_per_gene.R"),
    digest_helper = local_path("FASH_R3_COMPARE_DIGEST_HELPER",
      "code/revision_simulations/internal/r3_ideal_gaussian_measurement/ideal_gaussian_measurement.R"),
    base_reconstruction = file.path(root,
      "code/revision_simulations/internal/r3_two_stage_functional_testing/reconstruct_r3.R"),
    comparator = file.path(root,
      "code/revision_simulations/r1_r2_fashr0143/source_snapshots/r1_r2_fashr0143_simulation_functions.R")
  )
  hashes <- vapply(paths, r3c_sha, character(1))
  expected <- c(base$config$contract$base_source_sha[names(paths)[1:4]],
                comparator = r3c_contract()$comparator_sha256)
  r3c_assert(identical(hashes, expected), "A frozen scientific source changed.")
  simulation <- new.env(parent = globalenv())
  for (path in paths[1:4]) sys.source(path, envir = simulation)
  comparison <- new.env(parent = globalenv())
  sys.source(paths[["comparator"]], envir = comparison)
  description <- utils::packageDescription("fashr")
  r3c_assert(identical(description$Version, r3c_contract()$package_version) &&
               identical(description$RemoteSha, r3c_contract()$package_sha),
             "The installed fashr build differs from the pinned R3 build.")
  list(sim = simulation, cmp = comparison, paths = paths, hashes = hashes)
}

r3c_scores_from_fdr <- function(table, n) {
  r3c_assert(all(c("index", "lfdr", "FDR") %in% names(table)) && nrow(table) == n &&
               !anyDuplicated(table$index) && setequal(table$index, seq_len(n)),
             "FDR table must contain every unit index exactly once.")
  table <- table[match(seq_len(n), table$index), , drop = FALSE]
  data.frame(unit_index = seq_len(n), score = as.numeric(table$lfdr),
             adjusted_score = as.numeric(table$FDR))
}

r3c_validate_scores <- function(scores, n) {
  r3c_assert(is.data.frame(scores) && nrow(scores) == n &&
               identical(scores$unit_index, seq_len(n)) &&
               all(is.finite(scores$score)) && all(is.finite(scores$adjusted_score)) &&
               all(scores$score >= 0 & scores$score <= 1) &&
               all(scores$adjusted_score >= 0 & scores$adjusted_score <= 1),
             "Scores must be finite probabilities in the original all-unit order.")
}

r3c_curves <- function(scores, truth, seed, mechanism, alpha_grid = r3c_contract()$alpha_grid) {
  r3c_assert(is.logical(truth) && !anyNA(truth) && any(truth), "Invalid dynamic truth.")
  rows <- lapply(names(scores), function(method) {
    x <- scores[[method]]
    r3c_validate_scores(x, length(truth))
    do.call(rbind, lapply(alpha_grid, function(alpha) {
      selected <- x$adjusted_score <= alpha
      calls <- sum(selected)
      true_calls <- sum(selected & truth)
      data.frame(seed = seed, truth_mechanism = mechanism, method = method, alpha = alpha,
        calls = calls, true_calls = true_calls, false_calls = calls - true_calls,
        power = true_calls / sum(truth), empirical_fdr = (calls - true_calls) / max(calls, 1L))
    }))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

r3c_summarize <- function(curves, contract = r3c_contract()) {
  expected <- expand.grid(seed = contract$seeds, truth_mechanism = contract$mechanisms,
    method = contract$methods, alpha = contract$alpha_grid, stringsAsFactors = FALSE)
  key <- function(x) paste(x$seed, x$truth_mechanism, x$method, sprintf("%.3f", x$alpha), sep = ":")
  r3c_assert(!anyDuplicated(key(curves)) && setequal(key(curves), key(expected)),
             "Curves require every seed, truth, method, and alpha exactly once.")
  groups <- expand.grid(truth_mechanism = contract$mechanisms, method = contract$methods,
                        alpha = contract$alpha_grid, stringsAsFactors = FALSE)
  rows <- lapply(seq_len(nrow(groups)), function(i) {
    g <- groups[i, , drop = FALSE]
    x <- curves[curves$truth_mechanism == g$truth_mechanism & curves$method == g$method &
                  abs(curves$alpha - g$alpha) < 1e-12, , drop = FALSE]
    r3c_assert(all(is.finite(x$power)) && all(is.finite(x$empirical_fdr)), "Non-finite rates.")
    cbind(g, data.frame(n_seeds = nrow(x), mean_power = mean(x$power),
      power_min = min(x$power), power_max = max(x$power),
      mean_empirical_fdr = mean(x$empirical_fdr), empirical_fdr_min = min(x$empirical_fdr),
      empirical_fdr_max = max(x$empirical_fdr), mean_calls = mean(x$calls)))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

r3c_validate_inputs <- function(inputs, genotype, record, config, sim) {
  r3c_assert(identical(inputs$seed, record$seed) &&
               identical(inputs$truth_mechanism, record$truth_mechanism) &&
               identical(genotype$selection$pair_key, record$selected_pair_keys) &&
               identical(inputs$effects$unit_info$variant_id, record$selected_pair_keys) &&
               identical(inputs$effects$unit_info$effect_class == "dynamic_bspline",
                         record$inference$true_dynamic), "R3 input identity or pair order changed.")
  actual <- c(
    genotype_content_md5 = sim$genotype_content_md5(genotype$selection$pair_key,
                                                   rownames(genotype$G), genotype$G),
    true_beta_md5 = sim$serialized_object_md5(inputs$effects$beta_matrix),
    adjusted_se_md5 = sim$serialized_object_md5(inputs$regression$se),
    regression_beta_hat_md5 = sim$serialized_object_md5(inputs$regression$beta_hat)
  )
  r3c_assert(identical(actual, record$input_digests),
             "Comparator inputs do not match the four retained R3 digests.")
  covariates <- sim$simulate_covariate_matrix(config$n_donors, config$n_covariates,
                                             seed = inputs$component_seeds[["covariates"]])
  expression <- sim$simulate_eqtl_expression_from_genotypes(
    G = genotype$G, beta_matrix = inputs$effects$beta_matrix, time_grid = config$time_grid,
    covariates = covariates, expression_noise_sd = config$expression_noise_sd,
    covariate_effect_sd = config$covariate_effect_sd, intercept_sd = config$intercept_sd,
    seed = inputs$component_seeds[["expression"]])
  r3c_assert(isTRUE(all.equal(inputs$covariates, covariates, tolerance = 0)) &&
               isTRUE(all.equal(inputs$expression$expression, expression$expression, tolerance = 0)),
             "Retained covariates or individual-level expression differ from frozen R3 generation.")
  r3c_assert(isTRUE(all.equal(inputs$effects$true_functionals[, c("early", "middle", "late", "switch"), drop = FALSE],
                               record$inference$true_functionals, tolerance = 0)),
             "Functional truth changed during comparator pairing.")
  actual
}

r3c_fit_linear <- function(inputs, config, cmp) {
  contract <- r3c_contract()
  datasets <- cmp$make_fash_datasets_from_eqtl_summary(
    inputs$regression$beta_hat, inputs$regression$se, inputs$effects$beta_matrix,
    config$time_grid, inputs$effects$unit_info, inputs$scenario)
  raw <- cmp$fit_linear_mixture_fash(datasets, grid = contract$grid,
                                    pred_step = contract$pred_step, penalty = contract$penalty)
  bf <- cmp$BF_update_linear_mixture_fash(raw)
  for (fit in list(raw, bf)) cmp$validate_linear_mixture_fash(
    fit, expected_grid = contract$grid, expected_pred_step = contract$pred_step,
    expected_penalty = contract$penalty)
  list(raw = raw, bf = bf)
}

r3c_linear_scores <- function(fits, inputs, cmp) {
  methods <- c(raw = "FASH-linear-Raw", bf = "FASH-linear-BF")
  out <- lapply(names(methods), function(label) {
    result <- cmp$evaluate_lfdr_method(fits[[label]]$lfdr,
      unit_info = inputs$effects$unit_info, method = methods[[label]], target = "dynamic", alpha = 0.05)
    result[, c("unit_index", "score", "adjusted_score")]
  })
  names(out) <- unname(methods)
  out
}

r3c_direct_payload <- function(inputs, genotype, config) {
  list(genotype = genotype$G, expression = inputs$expression$expression,
       covariates = inputs$covariates, unit_info = inputs$effects$unit_info,
       settings = list(time_grid = config$time_grid, alpha = 0.05))
}

r3c_direct_scores <- function(payload, permutation, cmp, true_pi0 = r3c_contract()$true_pi0) {
  tests <- lapply(1:2, function(degree) cmp$evaluate_direct_interaction_efdr(
    payload, permutation, interaction_degree = degree, alpha = 0.05,
    pi0_method = "fixed", fixed_pi0 = true_pi0, method_adjustment = "eFDR_true_pi0"))
  names(tests) <- c("Direct-linear-LRT-eFDR-true-pi0", "Direct-quadratic-LRT-eFDR-true-pi0")
  list(scores = lapply(tests, function(x) x$result[, c("unit_index", "score", "adjusted_score")]),
       diagnostics = lapply(tests, function(x) list(lrt = x$lrt, threshold_table = x$efdr$threshold_table)))
}

r3c_validate_direct <- function(scores, diagnostics, n_null, pi0) {
  for (method in names(diagnostics)) {
    x <- scores[[method]]
    table <- diagnostics[[method]]$threshold_table
    r3c_assert(nrow(table) > 0L && all(diff(table$threshold) > 0) &&
                 all(table$null_count >= 0 & table$null_count <= n_null) &&
                 all(diff(table$null_count) >= 0) &&
                 identical(table$pi0_method, rep("fixed", nrow(table))) &&
                 all(table$pi0_multiplier == pi0), "Invalid direct-test threshold diagnostics.")
    observed_count <- findInterval(table$threshold, sort(x$score))
    efdr <- pmin(1, pmax(0, pi0 * (table$null_count / n_null) / (observed_count / nrow(x))))
    qvalue <- rev(cummin(rev(efdr)))
    r3c_assert(identical(as.integer(table$observed_count), observed_count) &&
                 isTRUE(all.equal(table$null_cdf, table$null_count / n_null)) &&
                 isTRUE(all.equal(table$observed_cdf, observed_count / nrow(x))) &&
                 isTRUE(all.equal(table$efdr, efdr)) &&
                 isTRUE(all.equal(table$qvalue, qvalue)) &&
                 isTRUE(all.equal(x$adjusted_score, qvalue[match(x$score, table$threshold)])),
               "Direct-test eFDR scores do not reproduce from saved counts.")
  }
  invisible(TRUE)
}

r3c_validate_record <- function(record, reference, contract = r3c_contract()) {
  r3c_assert(identical(record$schema, contract$schema) &&
               identical(record$seed, reference$seed) &&
               identical(record$truth_mechanism, reference$truth_mechanism) &&
               identical(record$selected_pair_keys, reference$selected_pair_keys) &&
               identical(record$input_digests, reference$input_digests) &&
               identical(record$true_dynamic, reference$inference$true_dynamic) &&
               identical(names(record$scores), contract$methods),
             "Comparison record identity, truth, input pairing, or method set is invalid.")
  for (scores in record$scores) r3c_validate_scores(scores, contract$n_units)
  expected <- r3c_scores_from_fdr(reference$inference$fdr_table, contract$n_units)
  r3c_assert(isTRUE(all.equal(record$scores[["FASH-IWP1-BF"]], expected, tolerance = 0)),
             "The preserved IWP1 BF decisions changed.")
  for (method in contract$methods[1:4]) {
    x <- record$scores[[method]]
    order <- order(x$score, x$unit_index)
    expected <- cumsum(x$score[order]) / seq_len(nrow(x))
    r3c_assert(isTRUE(all.equal(x$adjusted_score[order], expected, tolerance = 1e-12)),
               paste("Cumulative lfdr does not reproduce:", method))
  }
  r3c_assert(identical(names(record$direct_diagnostics), contract$methods[5:6]),
             "Both direct-test diagnostics are required.")
  r3c_validate_direct(record$scores, record$direct_diagnostics,
                      contract$n_units * contract$n_permutations, contract$true_pi0)
  permutation <- record$permutation_index
  r3c_assert(identical(dim(permutation), c(contract$n_permutations, 19L)) &&
               all(apply(permutation, 1, function(x) identical(sort(as.integer(x)), 1:19))) &&
               identical(record$permutation_settings$seed,
                         record$seed + contract$permutation_seed_offset) &&
               isTRUE(record$permutation_settings$permute_covariates_with_expression),
             "The donor permutation contract changed.")
  r3c_assert(identical(names(record$prior_weights), contract$methods[1:4]),
             "Raw and BF prior-weight diagnostics must be retained for both FASH methods.")
  invisible(TRUE)
}
