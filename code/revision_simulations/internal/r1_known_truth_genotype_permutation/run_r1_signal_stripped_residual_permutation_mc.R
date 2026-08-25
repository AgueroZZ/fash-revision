#!/usr/bin/env Rscript

# Run the five-seed R1 signal-stripped donor-residual-permutation experiment.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1L] == length(args)) {
    return(default)
  }
  args[hit[1L] + 1L]
}

as_flag <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

capture_warnings <- function(expression) {
  warning_messages <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(warning) {
      warning_messages <<- c(warning_messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warning_messages))
}

save_rds_atomically <- function(object, path) {
  temporary_path <- paste0(path, ".tmp_", Sys.getpid())
  saveRDS(object, temporary_path)
  if (!file.rename(temporary_path, path)) {
    unlink(temporary_path)
    stop("Could not atomically save replicate cache: ", path)
  }
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

rebuild_r1_source <- function(seed,
                              num_cores,
                              output_directory,
                              configuration) {
  run_genotype_level_bspline_eqtl_simulation(
    n_donors = configuration$n_donors,
    n_variants = configuration$J,
    time_grid = configuration$time_grid,
    n_covariates = configuration$n_covariates,
    class_probs = configuration$class_probs,
    expression_noise_sd = configuration$expression_noise_sd,
    dynamic_main_effect_sd = configuration$dynamic_main_effect_sd,
    scenario = configuration$scenario,
    alpha = 0.05,
    seed = seed,
    estimate_sigma = TRUE,
    sigma_beta_grid = configuration$linear_sigma_grid,
    num_cores = num_cores,
    num_basis = configuration$num_basis,
    output_dir = output_directory,
    save_outputs = FALSE,
    verbose = FALSE
  )
}

validate_source_structure <- function(source_object,
                                      expected_grid,
                                      n_alternatives) {
  source_fit <- source_object$fash_fits$fash_iwp1_raw
  effect_class <- as.character(source_object$unit_info$effect_class)
  checks <- c(
    identical(dim(source_object$genotype), c(19L, 1000L)),
    identical(dim(source_object$expression), c(19L, 1000L, 16L)),
    identical(dim(source_object$covariates), c(19L, 5L)),
    sum(effect_class == "dynamic_bspline") == 200L,
    sum(effect_class == "zero") == 400L,
    identical(source_fit$settings$order, 1),
    identical(source_fit$settings$num_basis, 20L),
    identical(source_fit$settings$pred_step, 1),
    identical(source_fit$settings$penalty, 10),
    isTRUE(all.equal(source_fit$psd_grid, expected_grid, tolerance = 0)),
    nrow(source_fit$L_matrix) == 1000L,
    ncol(source_fit$L_matrix) == length(expected_grid),
    n_alternatives <= sum(effect_class == "dynamic_bspline")
  )
  if (!all(checks)) {
    stop("A reconstructed R1 source object failed structural validation.")
  }
  invisible(TRUE)
}

compare_source_to_cache <- function(rebuilt, cached) {
  data.frame(
    component = c(
      "genotype", "expression", "true_beta", "beta_hat", "se",
      "IWP likelihood", "IWP lfdr"
    ),
    maximum_absolute_difference = c(
      max(abs(rebuilt$genotype - cached$genotype)),
      max(abs(rebuilt$expression - cached$expression)),
      max(abs(rebuilt$true_beta - cached$true_beta)),
      max(abs(rebuilt$eqtl_summary$beta_hat - cached$eqtl_summary$beta_hat)),
      max(abs(rebuilt$eqtl_summary$se - cached$eqtl_summary$se)),
      max(abs(
        rebuilt$fash_fits$fash_iwp1_raw$L_matrix -
          cached$fash_fits$fash_iwp1_raw$L_matrix
      )),
      max(abs(
        rebuilt$fash_fits$fash_iwp1_raw$lfdr -
          cached$fash_fits$fash_iwp1_raw$lfdr
      ))
    ),
    tolerance = rep(1e-12, 7L),
    stringsAsFactors = FALSE
  )
}

summarize_known_truth_mc <- function(alpha_rows) {
  if (anyDuplicated(alpha_rows[c("seed", "arm", "fit_stage", "alpha")])) {
    stop("The replicate alpha curves contain duplicated keys.")
  }
  groups <- split(
    alpha_rows,
    list(alpha_rows$arm, alpha_rows$fit_stage, alpha_rows$alpha),
    drop = TRUE
  )
  rows <- lapply(groups, function(x) {
    fdr <- summarize_mc_values(x$realized_fdp)
    power <- summarize_mc_values(x$power)
    discoveries <- summarize_mc_values(x$n_discoveries)
    false_discoveries <- summarize_mc_values(x$false_discoveries)
    data.frame(
      arm = x$arm[1L],
      fit_stage = x$fit_stage[1L],
      alpha = x$alpha[1L],
      n_replications = length(unique(x$seed)),
      mean_discoveries = discoveries[["mean"]],
      mean_false_discoveries = false_discoveries[["mean"]],
      mean_power = power[["mean"]],
      power_sd = power[["sd"]],
      power_mc_se = power[["se"]],
      power_ci_lower = pmax(0, power[["lower"]]),
      power_ci_upper = pmin(1, power[["upper"]]),
      mean_fdr = fdr[["mean"]],
      fdr_sd = fdr[["sd"]],
      fdr_mc_se = fdr[["se"]],
      fdr_ci_lower = pmax(0, fdr[["lower"]]),
      fdr_ci_upper = pmin(1, fdr[["upper"]]),
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, rows)
  output <- output[order(output$fit_stage, output$arm, output$alpha), ]
  rownames(output) <- NULL
  output
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation",
  "selected_signal_genotype_permutation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "r1_known_truth_genotype_permutation",
  "r1_known_truth_genotype_permutation_helpers.R"
))

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required.")
}

source_seeds <- c(12345L, 22345L, 32345L, 42345L, 52345L)
n_alternatives <- as.integer(get_arg("--n-alternatives", "200"))
num_cores <- as.integer(get_arg("--num-cores", "8"))
overwrite <- as_flag(get_arg("--overwrite", "false"))
if (length(n_alternatives) != 1L || is.na(n_alternatives) ||
    n_alternatives < 2L || n_alternatives > 200L ||
    length(num_cores) != 1L || is.na(num_cores) || num_cores < 1L) {
  stop("Invalid alternative count or core count.")
}

formal_r1_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "mc",
  "r1_random_bspline_main_effect_profile_sigma_pilot5"
)
formal_configuration_path <- file.path(
  formal_r1_directory, "configuration.rds"
)
cached_source_path <- file.path(
  formal_r1_directory, "full_fits", "seed_12345.rds"
)
control_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  paste0("r1_known_truth_genotype_permutation_mc5_J", n_alternatives, "_v1")
)
control_replicate_paths <- file.path(
  control_directory, "replicates", paste0("seed_", source_seeds, ".rds")
)
required_paths <- c(
  formal_configuration_path, cached_source_path, control_replicate_paths
)
if (any(!file.exists(required_paths))) {
  stop("A formal R1 or validated control-experiment artifact is missing.")
}
formal_configuration <- readRDS(formal_configuration_path)
expected_configuration <- list(
  J = 1000L,
  n_donors = 19L,
  n_covariates = 5L,
  expression_noise_sd = 1,
  dynamic_main_effect_sd = 1,
  num_basis = 20L,
  scenario = "genotype_random_bspline_main_effect_dynamic_eqtl",
  class_probs = c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40),
  seed_list = source_seeds
)
for (field in names(expected_configuration)) {
  if (!isTRUE(all.equal(
    formal_configuration[[field]], expected_configuration[[field]],
    tolerance = 0
  ))) {
    stop("The formal R1 configuration has changed at field: ", field)
  }
}

output_id <- paste0(
  "r1_signal_stripped_unadjusted_residual_permutation_mc5_J",
  n_alternatives,
  "_v1"
)
output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal", output_id
)
replicate_directory <- file.path(output_directory, "replicates")
summary_directory <- file.path(output_directory, "summary")
figure_directory <- file.path(output_directory, "figures")
invisible(lapply(
  c(output_directory, replicate_directory, summary_directory, figure_directory),
  dir.create, recursive = TRUE, showWarnings = FALSE
))

configuration <- list(
  experiment = paste(
    "Five-seed R1 signal-stripped synchronized donor-residual-permutation",
    "calibration diagnostic"
  ),
  output_id = output_id,
  source_seeds = source_seeds,
  permutation_seeds = vapply(
    source_seeds,
    function(seed) revision_component_seeds(seed)[["permutations"]],
    integer(1)
  ),
  n_alternatives = n_alternatives,
  n_null_per_arm = n_alternatives,
  n_units_per_arm = 2L * n_alternatives,
  true_pi0 = 0.5,
  null_construction = paste(
    "Full per-time OLS residuals after intercept, genotype, and five",
    "covariates; unadjusted full-model residuals; synchronized donor permutation",
    "across all units and time points; add back nuisance fitted values only;",
    "analyze against original genotypes."
  ),
  leverage_adjustment = "none",
  donor_maps = paste(
    "Exact maps from the validated naive genotype-permutation control",
    "experiment."
  ),
  control_output_id = basename(control_directory),
  source_configuration = formal_configuration,
  alpha_grid = c(0, seq(0.001, 0.20, by = 0.001)),
  num_cores = num_cores,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  interpretation_boundary = paste(
    "This fitted-model mechanism experiment isolates source-signal leakage;",
    "it does not establish exchangeability or exact real-data validity."
  )
)
configuration_path <- file.path(output_directory, "configuration.rds")
if (file.exists(configuration_path) && !overwrite) {
  cached_configuration <- readRDS(configuration_path)
  stable_fields <- setdiff(names(configuration), c("generated_at", "num_cores"))
  if (!isTRUE(all.equal(
    cached_configuration[stable_fields], configuration[stable_fields]
  ))) {
    stop("The existing residual-permutation output has incompatible settings.")
  }
} else {
  saveRDS(configuration, configuration_path)
}

cached_source <- readRDS(cached_source_path)
expected_grid <- default_revision_grid()
alpha_grid <- configuration$alpha_grid
replicates <- vector("list", length(source_seeds))
names(replicates) <- as.character(source_seeds)

for (seed_index in seq_along(source_seeds)) {
  seed <- source_seeds[seed_index]
  replicate_path <- file.path(
    replicate_directory, paste0("seed_", seed, ".rds")
  )
  if (file.exists(replicate_path) && !overwrite) {
    message("Reusing residual-permutation replicate cache: ", replicate_path)
    replicate <- readRDS(replicate_path)
    expected_fields <- c(
      "seed", "permutation_seed", "alpha_curve", "pi0_summary",
      "null_bf_summary", "unit_results", "donor_map",
      "construction_diagnostics", "validation"
    )
    if (!all(expected_fields %in% names(replicate)) ||
        !identical(replicate$seed, seed)) {
      stop("An existing residual-permutation replicate cache is invalid.")
    }
    replicates[[as.character(seed)]] <- replicate
    next
  }

  message("Reconstructing formal R1 source seed ", seed, ".")
  replicate_start <- proc.time()[["elapsed"]]
  source_capture <- capture_warnings(rebuild_r1_source(
    seed = seed,
    num_cores = num_cores,
    output_directory = output_directory,
    configuration = formal_configuration
  ))
  source_object <- source_capture$value
  validate_source_structure(source_object, expected_grid, n_alternatives)

  source_reconstruction <- NULL
  if (seed == source_seeds[1L]) {
    source_reconstruction <- compare_source_to_cache(
      source_object, cached_source
    )
    if (any(
      source_reconstruction$maximum_absolute_difference >
        source_reconstruction$tolerance
    ) || !identical(
      source_object$unit_info$effect_class,
      cached_source$unit_info$effect_class
    )) {
      stop("The current generator did not reproduce cached R1 seed 12345.")
    }
  }

  control_replicate <- readRDS(control_replicate_paths[seed_index])
  permutation_seed <- revision_component_seeds(seed)[["permutations"]]
  if (!identical(control_replicate$seed, seed) ||
      !identical(control_replicate$permutation_seed, permutation_seed) ||
      !all(control_replicate$validation$pass)) {
    stop("The naive control replicate failed provenance validation.")
  }

  source_fit <- source_object$fash_fits$fash_iwp1_raw
  effect_class <- as.character(source_object$unit_info$effect_class)
  alternative_indices <- which(effect_class == "dynamic_bspline")[
    seq_len(n_alternatives)
  ]
  alternative_keys <- as.character(
    source_object$unit_info$unit_id[alternative_indices]
  )
  alternative_genotype <- source_object$genotype[
    , alternative_indices, drop = FALSE
  ]
  colnames(alternative_genotype) <- alternative_keys
  alternative_expression <- source_object$expression[
    , alternative_indices, , drop = FALSE
  ]
  dimnames(alternative_expression)[[2L]] <- alternative_keys
  expected_permutation <- make_shared_genotype_permutation(
    alternative_genotype, seed = permutation_seed
  )
  if (!identical(
    expected_permutation$donor_map, control_replicate$donor_map
  )) {
    stop("The residual arm does not reproduce the naive arm donor map.")
  }

  message("Constructing signal-stripped residual nulls for seed ", seed, ".")
  residual_null <- make_signal_stripped_residual_null(
    genotype = alternative_genotype,
    expression = alternative_expression,
    covariates = source_object$covariates,
    donor_map = control_replicate$donor_map,
    leverage_adjustment = "none"
  )
  residual_eqtl <- estimate_eqtl_summaries_from_genotypes(
    G = alternative_genotype,
    expression = residual_null$null_expression,
    covariates = source_object$covariates,
    apply_t_se_correction = TRUE
  )
  residual_unit_info <- source_object$unit_info[
    alternative_indices, , drop = FALSE
  ]
  residual_unit_info$unit_id <- paste0(
    alternative_keys,
    "__signal_stripped_residual_null_seed",
    permutation_seed
  )
  residual_unit_info$effect_class <- "signal_stripped_residual_null"
  residual_unit_info$genetic_main_effect <- NA_real_
  residual_unit_info$scenario <-
    "r1_known_truth_signal_stripped_residual_permutation_mc"
  residual_true_beta <- matrix(
    0,
    nrow = n_alternatives,
    ncol = length(formal_configuration$time_grid),
    dimnames = dimnames(residual_eqtl$beta_hat)
  )
  residual_datasets <- make_fash_datasets_from_eqtl_summary(
    beta_hat = residual_eqtl$beta_hat,
    se = residual_eqtl$se,
    true_beta = residual_true_beta,
    time_grid = formal_configuration$time_grid,
    unit_info = residual_unit_info,
    scenario = "r1_known_truth_signal_stripped_residual_permutation_mc"
  )
  residual_keys <- as.character(residual_unit_info$unit_id)
  names(residual_datasets) <- residual_keys

  message("Computing residual-null likelihood rows for seed ", seed, ".")
  residual_capture <- capture_warnings(fashr::fash(
    Y = "y",
    smooth_var = "x",
    S = "sd",
    data_list = residual_datasets,
    num_basis = source_fit$settings$num_basis,
    order = source_fit$settings$order,
    betaprec = source_fit$settings$betaprec,
    pred_step = source_fit$settings$pred_step,
    penalty = source_fit$settings$penalty,
    grid = source_fit$psd_grid,
    num_cores = num_cores,
    verbose = FALSE
  ))
  residual_fit <- residual_capture$value
  names(residual_fit$fash_data$data_list) <- residual_keys
  names(residual_fit$fash_data$S) <- residual_keys
  rownames(residual_fit$L_matrix) <- residual_keys

  alternative_data <- source_fit$fash_data$data_list[alternative_indices]
  alternative_se <- source_fit$fash_data$S[alternative_indices]
  alternative_likelihood <- source_fit$L_matrix[
    alternative_indices, , drop = FALSE
  ]
  names(alternative_data) <- alternative_keys
  names(alternative_se) <- alternative_keys
  rownames(alternative_likelihood) <- alternative_keys
  unit_keys <- c(alternative_keys, residual_keys)
  true_null <- c(
    rep(FALSE, n_alternatives), rep(TRUE, n_alternatives)
  )
  merged_data <- c(
    alternative_data, residual_fit$fash_data$data_list
  )
  merged_se <- c(alternative_se, residual_fit$fash_data$S)
  merged_likelihood <- rbind(
    alternative_likelihood, residual_fit$L_matrix
  )

  raw_capture <- capture_warnings(refit_fash_from_likelihood(
    source_fit = source_fit,
    data_list = merged_data,
    se_list = merged_se,
    likelihood_matrix = merged_likelihood,
    unit_keys = unit_keys,
    penalty = source_fit$settings$penalty
  ))
  bf_capture <- capture_warnings(fashr::BF_update(
    raw_capture$value, plot = FALSE
  ))
  fits <- list(Raw = raw_capture$value, BF = bf_capture$value)
  alpha_rows <- list()
  pi0_rows <- list()
  unit_rows <- list()
  null_bf_summary <- NULL
  for (fit_stage in names(fits)) {
    fit <- fits[[fit_stage]]
    names(fit$lfdr) <- unit_keys
    alpha_rows[[fit_stage]] <- known_truth_alpha_curve(
      lfdr = fit$lfdr,
      true_null = true_null,
      alpha_grid = alpha_grid,
      arm = "signal_stripped_residual_permutation",
      fit_stage = fit_stage
    )
    pi0_rows[[fit_stage]] <- data.frame(
      seed = seed,
      permutation_seed = permutation_seed,
      arm = "signal_stripped_residual_permutation",
      fit_stage = fit_stage,
      estimated_pi0 = extract_pi0(fit),
      true_pi0 = mean(true_null),
      stringsAsFactors = FALSE
    )
    unit_rows[[fit_stage]] <- data.frame(
      seed = seed,
      permutation_seed = permutation_seed,
      arm = "signal_stripped_residual_permutation",
      fit_stage = fit_stage,
      unit_key = unit_keys,
      source_unit_id = rep(alternative_keys, 2L),
      true_null = true_null,
      lfdr = as.numeric(fit$lfdr),
      bayes_factor = if (fit_stage == "BF") {
        as.numeric(fit$BF)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
    if (fit_stage == "BF") {
      null_bf_summary <- cbind(
        seed = seed,
        permutation_seed = permutation_seed,
        summarize_null_bf(
          bf = fit$BF,
          true_null = true_null,
          arm = "signal_stripped_residual_permutation"
        )
      )
    }
  }
  alpha_curve <- do.call(rbind, alpha_rows)
  alpha_curve$seed <- seed
  alpha_curve$permutation_seed <- permutation_seed
  alpha_curve <- alpha_curve[, c(
    "seed", "permutation_seed",
    setdiff(names(alpha_curve), c("seed", "permutation_seed"))
  )]
  pi0_summary <- do.call(rbind, pi0_rows)
  unit_results <- do.call(rbind, unit_rows)
  construction_diagnostics <- residual_null$diagnostics
  construction_diagnostics$seed <- seed
  construction_diagnostics$permutation_seed <- permutation_seed
  construction_diagnostics$null_beta_hat_rms <- sqrt(mean(
    residual_eqtl$beta_hat^2
  ))
  construction_diagnostics$null_z_rms <- sqrt(mean(
    (residual_eqtl$beta_hat / residual_eqtl$se)^2
  ))
  construction_diagnostics$null_maximum_absolute_z <- max(abs(
    residual_eqtl$beta_hat / residual_eqtl$se
  ))
  construction_diagnostics$source_beta_hat_rms <- sqrt(mean(
    source_object$eqtl_summary$beta_hat[alternative_indices, ]^2
  ))
  construction_diagnostics$removed_genotype_fitted_rms <- sqrt(mean(
    residual_null$source_genotype_fitted^2
  ))

  warnings <- unique(c(
    source_capture$warnings,
    residual_capture$warnings,
    raw_capture$warnings,
    bf_capture$warnings
  ))
  validation <- data.frame(
    check = c(
      "alpha rows", "pi0 rows", "unit rows", "true pi0",
      "finite lfdr", "positive BF", "donor map",
      "full residual orthogonality", "nuisance genotype coefficient"
    ),
    pass = c(
      nrow(alpha_curve) == 2L * length(alpha_grid),
      nrow(pi0_summary) == 2L,
      nrow(unit_results) == 2L * 2L * n_alternatives,
      all(pi0_summary$true_pi0 == 0.5),
      all(is.finite(unit_results$lfdr)),
      all(
        is.na(unit_results$bayes_factor) |
          (is.finite(unit_results$bayes_factor) &
             unit_results$bayes_factor > 0)
      ),
      identical(residual_null$donor_map, control_replicate$donor_map),
      construction_diagnostics$
        maximum_full_residual_design_cross_product < 1e-10,
      construction_diagnostics$
        maximum_nuisance_partial_genotype_coefficient < 1e-10
    ),
    stringsAsFactors = FALSE
  )
  if (!all(validation$pass)) {
    stop("A residual-permutation replicate failed internal validation.")
  }

  replicate <- list(
    seed = seed,
    permutation_seed = permutation_seed,
    alpha_curve = alpha_curve,
    pi0_summary = pi0_summary,
    null_bf_summary = null_bf_summary,
    unit_results = unit_results,
    donor_map = residual_null$donor_map,
    construction_diagnostics = construction_diagnostics,
    validation = validation,
    source_reconstruction = source_reconstruction,
    warnings = warnings,
    elapsed_seconds = unname(proc.time()[["elapsed"]] - replicate_start)
  )
  save_rds_atomically(replicate, replicate_path)
  replicates[[as.character(seed)]] <- replicate
}

all_alpha <- do.call(rbind, lapply(replicates, `[[`, "alpha_curve"))
all_pi0 <- do.call(rbind, lapply(replicates, `[[`, "pi0_summary"))
all_null_bf <- do.call(rbind, lapply(replicates, `[[`, "null_bf_summary"))
all_units <- do.call(rbind, lapply(replicates, `[[`, "unit_results"))
all_diagnostics <- do.call(
  rbind,
  lapply(replicates, `[[`, "construction_diagnostics")
)
all_donor_maps <- do.call(rbind, lapply(replicates, function(x) {
  cbind(
    seed = x$seed,
    permutation_seed = x$permutation_seed,
    x$donor_map
  )
}))
mc_alpha <- summarize_known_truth_mc(all_alpha)
mc_alpha005 <- mc_alpha[abs(mc_alpha$alpha - 0.05) < 1e-12, ]
replicate_alpha005 <- all_alpha[
  abs(all_alpha$alpha - 0.05) < 1e-12, , drop = FALSE
]
if (nrow(all_alpha) != length(source_seeds) * 2L * length(alpha_grid) ||
    any(mc_alpha$n_replications != length(source_seeds)) ||
    nrow(mc_alpha005) != 2L || nrow(replicate_alpha005) != 10L ||
    !setequal(unique(all_alpha$seed), source_seeds)) {
  stop("The combined residual-permutation MC summaries are incomplete.")
}

write_csv(all_alpha, file.path(summary_directory, "replicate_alpha_curves.csv"))
write_csv(
  replicate_alpha005,
  file.path(summary_directory, "replicate_alpha005.csv")
)
write_csv(all_pi0, file.path(summary_directory, "replicate_pi0.csv"))
write_csv(
  all_null_bf,
  file.path(summary_directory, "replicate_null_bf_summary.csv")
)
write_csv(all_units, file.path(summary_directory, "replicate_unit_results.csv"))
write_csv(
  all_diagnostics,
  file.path(summary_directory, "construction_diagnostics.csv")
)
write_csv(
  all_donor_maps,
  file.path(summary_directory, "replicate_donor_permutations.csv")
)
write_csv(mc_alpha, file.path(summary_directory, "mc_alpha_curve.csv"))
write_csv(mc_alpha005, file.path(summary_directory, "mc_alpha005_summary.csv"))
write_csv(
  replicates[[as.character(source_seeds[1L])]]$source_reconstruction,
  file.path(summary_directory, "source_reconstruction_validation.csv")
)
runtime <- data.frame(
  seed = source_seeds,
  permutation_seed = vapply(replicates, `[[`, integer(1), "permutation_seed"),
  elapsed_seconds = vapply(replicates, `[[`, numeric(1), "elapsed_seconds"),
  warning_count = vapply(replicates, function(x) length(x$warnings), integer(1)),
  stringsAsFactors = FALSE
)
write_csv(runtime, file.path(summary_directory, "runtime.csv"))

control_alpha005 <- do.call(rbind, lapply(control_replicate_paths, function(path) {
  rows <- readRDS(path)$alpha_curve
  rows[abs(rows$alpha - 0.05) < 1e-12, , drop = FALSE]
}))
genuine_alpha005 <- control_alpha005[
  control_alpha005$arm == "genuine_null_baseline",
  c("seed", "fit_stage", "realized_fdp"),
  drop = FALSE
]
naive_alpha005 <- control_alpha005[
  control_alpha005$arm == "shared_genotype_permutation",
  c("seed", "fit_stage", "realized_fdp"),
  drop = FALSE
]
names(genuine_alpha005)[3L] <- "genuine_null_fdp"
names(naive_alpha005)[3L] <- "naive_genotype_permutation_fdp"
residual_alpha005 <- replicate_alpha005[
  , c("seed", "fit_stage", "realized_fdp"), drop = FALSE
]
names(residual_alpha005)[3L] <- "signal_stripped_residual_fdp"
paired_alpha005 <- Reduce(
  function(x, y) merge(x, y, by = c("seed", "fit_stage")),
  list(genuine_alpha005, naive_alpha005, residual_alpha005)
)
paired_alpha005$residual_minus_genuine <-
  paired_alpha005$signal_stripped_residual_fdp -
  paired_alpha005$genuine_null_fdp
paired_alpha005$naive_minus_genuine <-
  paired_alpha005$naive_genotype_permutation_fdp -
  paired_alpha005$genuine_null_fdp
paired_difference_rows <- list()
row_index <- 0L
for (fit_stage in c("Raw", "BF")) {
  rows <- paired_alpha005[paired_alpha005$fit_stage == fit_stage, ]
  for (contrast in c("residual_minus_genuine", "naive_minus_genuine")) {
    row_index <- row_index + 1L
    values <- summarize_mc_values(rows[[contrast]])
    paired_difference_rows[[row_index]] <- data.frame(
      fit_stage = fit_stage,
      contrast = contrast,
      n_replications = nrow(rows),
      mean_difference = values[["mean"]],
      difference_sd = values[["sd"]],
      difference_mc_se = values[["se"]],
      difference_ci_lower = values[["lower"]],
      difference_ci_upper = values[["upper"]],
      stringsAsFactors = FALSE
    )
  }
}
paired_difference_summary <- do.call(rbind, paired_difference_rows)
write_csv(
  paired_alpha005,
  file.path(summary_directory, "paired_fdr_alpha005.csv")
)
write_csv(
  paired_difference_summary,
  file.path(summary_directory, "paired_fdr_difference_alpha005.csv")
)

print(mc_alpha005)
print(all_pi0)
cat("R1 signal-stripped residual-permutation MC completed.\n")
