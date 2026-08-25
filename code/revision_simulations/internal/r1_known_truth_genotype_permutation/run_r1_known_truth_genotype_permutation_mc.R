#!/usr/bin/env Rscript

# Run the five-seed Monte Carlo extension of the R1 known-truth
# shared-genotype-permutation diagnostic.

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

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

save_rds_atomically <- function(object, path) {
  temporary_path <- paste0(path, ".tmp_", Sys.getpid())
  saveRDS(object, temporary_path)
  if (!file.rename(temporary_path, path)) {
    unlink(temporary_path)
    stop("Could not atomically save replicate cache: ", path)
  }
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
    n_alternatives <= sum(effect_class == "dynamic_bspline"),
    n_alternatives <= sum(effect_class == "zero")
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
  required_columns <- c(
    "seed", "arm", "fit_stage", "alpha", "n_discoveries",
    "false_discoveries", "true_positives", "realized_fdp", "power"
  )
  missing_columns <- setdiff(required_columns, names(alpha_rows))
  if (length(missing_columns) > 0L ||
      anyDuplicated(alpha_rows[c("seed", "arm", "fit_stage", "alpha")])) {
    stop("The replicate alpha curves are incomplete or duplicated.")
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
  output <- output[order(
    output$fit_stage, output$arm, output$alpha
  ), , drop = FALSE]
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
if (!file.exists(formal_configuration_path) || !file.exists(cached_source_path)) {
  stop("The formal R1 configuration or seed-12345 full fit is missing.")
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
  "r1_known_truth_genotype_permutation_mc5_J", n_alternatives, "_v1"
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
    "Five-seed Monte Carlo extension of the R1 known-truth",
    "shared-genotype-permutation diagnostic"
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
  source_configuration = formal_configuration,
  source_reconstruction = paste(
    "Exact formal R1 generator; validated against cached seed 12345 for",
    "genotype, expression, truth, beta/SE, IWP likelihood, and lfdr."
  ),
  alpha_grid = c(0, seq(0.001, 0.20, by = 0.001)),
  num_cores = num_cores,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  interpretation_boundary = paste(
    "Five source seeds provide a Monte Carlo diagnostic, but uncertainty",
    "intervals remain wide and do not constitute formal calibration evidence."
  )
)
configuration_path <- file.path(output_directory, "configuration.rds")
if (file.exists(configuration_path) && !overwrite) {
  cached_configuration <- readRDS(configuration_path)
  stable_fields <- setdiff(names(configuration), c("generated_at", "num_cores"))
  if (!isTRUE(all.equal(
    cached_configuration[stable_fields], configuration[stable_fields]
  ))) {
    stop("The existing MC output has incompatible settings.")
  }
} else {
  saveRDS(configuration, configuration_path)
}

cached_source <- readRDS(cached_source_path)
expected_grid <- default_revision_grid()
alpha_grid <- configuration$alpha_grid
replicates <- vector("list", length(source_seeds))
names(replicates) <- as.character(source_seeds)

for (seed in source_seeds) {
  replicate_path <- file.path(
    replicate_directory, paste0("seed_", seed, ".rds")
  )
  if (file.exists(replicate_path) && !overwrite) {
    message("Reusing MC replicate cache: ", replicate_path)
    replicate <- readRDS(replicate_path)
    expected_fields <- c(
      "seed", "permutation_seed", "alpha_curve", "pi0_summary",
      "null_bf_summary", "unit_results", "donor_map", "validation"
    )
    if (!all(expected_fields %in% names(replicate)) ||
        !identical(replicate$seed, seed)) {
      stop("An existing MC replicate cache is invalid: ", replicate_path)
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

  source_fit <- source_object$fash_fits$fash_iwp1_raw
  effect_class <- as.character(source_object$unit_info$effect_class)
  alternative_indices <- which(effect_class == "dynamic_bspline")[
    seq_len(n_alternatives)
  ]
  genuine_null_indices <- which(effect_class == "zero")[
    seq_len(n_alternatives)
  ]
  alternative_keys <- as.character(
    source_object$unit_info$unit_id[alternative_indices]
  )
  genuine_null_keys <- as.character(
    source_object$unit_info$unit_id[genuine_null_indices]
  )
  alternative_genotype <- source_object$genotype[
    , alternative_indices, drop = FALSE
  ]
  colnames(alternative_genotype) <- alternative_keys
  alternative_expression <- source_object$expression[
    , alternative_indices, , drop = FALSE
  ]

  permutation_seed <- revision_component_seeds(seed)[["permutations"]]
  shared_permutation <- make_shared_genotype_permutation(
    alternative_genotype, seed = permutation_seed
  )
  if (all(shared_permutation$donor_map$fixed_point)) {
    stop("A fixed-seed donor permutation was the identity permutation.")
  }
  permuted_eqtl <- estimate_eqtl_summaries_from_genotypes(
    G = shared_permutation$genotype,
    expression = alternative_expression,
    covariates = source_object$covariates,
    apply_t_se_correction = TRUE
  )
  permuted_unit_info <- source_object$unit_info[
    alternative_indices, , drop = FALSE
  ]
  permuted_unit_info$unit_id <- paste0(
    alternative_keys, "__permuted_null_seed", permutation_seed
  )
  permuted_unit_info$effect_class <- "permuted_null"
  permuted_unit_info$genetic_main_effect <- NA_real_
  permuted_unit_info$scenario <-
    "r1_known_truth_shared_genotype_permutation_mc"
  permuted_true_beta <- matrix(
    0,
    nrow = n_alternatives,
    ncol = length(formal_configuration$time_grid),
    dimnames = dimnames(permuted_eqtl$beta_hat)
  )
  permuted_datasets <- make_fash_datasets_from_eqtl_summary(
    beta_hat = permuted_eqtl$beta_hat,
    se = permuted_eqtl$se,
    true_beta = permuted_true_beta,
    time_grid = formal_configuration$time_grid,
    unit_info = permuted_unit_info,
    scenario = "r1_known_truth_shared_genotype_permutation_mc"
  )
  permuted_keys <- as.character(permuted_unit_info$unit_id)
  names(permuted_datasets) <- permuted_keys

  message("Computing permuted likelihood rows for source seed ", seed, ".")
  permuted_capture <- capture_warnings(fashr::fash(
    Y = "y",
    smooth_var = "x",
    S = "sd",
    data_list = permuted_datasets,
    num_basis = source_fit$settings$num_basis,
    order = source_fit$settings$order,
    betaprec = source_fit$settings$betaprec,
    pred_step = source_fit$settings$pred_step,
    penalty = source_fit$settings$penalty,
    grid = source_fit$psd_grid,
    num_cores = num_cores,
    verbose = FALSE
  ))
  permuted_fit <- permuted_capture$value
  names(permuted_fit$fash_data$data_list) <- permuted_keys
  names(permuted_fit$fash_data$S) <- permuted_keys
  rownames(permuted_fit$L_matrix) <- permuted_keys

  alternative_data <- source_fit$fash_data$data_list[alternative_indices]
  alternative_se <- source_fit$fash_data$S[alternative_indices]
  alternative_likelihood <- source_fit$L_matrix[
    alternative_indices, , drop = FALSE
  ]
  names(alternative_data) <- alternative_keys
  names(alternative_se) <- alternative_keys
  rownames(alternative_likelihood) <- alternative_keys
  genuine_null_data <- source_fit$fash_data$data_list[genuine_null_indices]
  genuine_null_se <- source_fit$fash_data$S[genuine_null_indices]
  genuine_null_likelihood <- source_fit$L_matrix[
    genuine_null_indices, , drop = FALSE
  ]
  names(genuine_null_data) <- genuine_null_keys
  names(genuine_null_se) <- genuine_null_keys
  rownames(genuine_null_likelihood) <- genuine_null_keys

  arm_inputs <- list(
    genuine_null_baseline = list(
      unit_keys = c(alternative_keys, genuine_null_keys),
      likelihood = rbind(
        alternative_likelihood, genuine_null_likelihood
      ),
      data_list = c(alternative_data, genuine_null_data),
      se_list = c(alternative_se, genuine_null_se)
    ),
    shared_genotype_permutation = list(
      unit_keys = c(alternative_keys, permuted_keys),
      likelihood = rbind(
        alternative_likelihood, permuted_fit$L_matrix
      ),
      data_list = c(
        alternative_data, permuted_fit$fash_data$data_list
      ),
      se_list = c(alternative_se, permuted_fit$fash_data$S)
    )
  )
  true_null <- c(
    rep(FALSE, n_alternatives), rep(TRUE, n_alternatives)
  )
  alpha_rows <- list()
  pi0_rows <- list()
  null_bf_rows <- list()
  unit_rows <- list()
  result_index <- 0L
  for (arm in names(arm_inputs)) {
    input <- arm_inputs[[arm]]
    raw_capture <- capture_warnings(refit_fash_from_likelihood(
      source_fit = source_fit,
      data_list = input$data_list,
      se_list = input$se_list,
      likelihood_matrix = input$likelihood,
      unit_keys = input$unit_keys,
      penalty = source_fit$settings$penalty
    ))
    bf_capture <- capture_warnings(fashr::BF_update(
      raw_capture$value, plot = FALSE
    ))
    fits <- list(Raw = raw_capture$value, BF = bf_capture$value)
    for (fit_stage in names(fits)) {
      result_index <- result_index + 1L
      fit <- fits[[fit_stage]]
      names(fit$lfdr) <- input$unit_keys
      alpha_rows[[result_index]] <- known_truth_alpha_curve(
        lfdr = fit$lfdr,
        true_null = true_null,
        alpha_grid = alpha_grid,
        arm = arm,
        fit_stage = fit_stage
      )
      pi0_rows[[result_index]] <- data.frame(
        seed = seed,
        permutation_seed = permutation_seed,
        arm = arm,
        fit_stage = fit_stage,
        estimated_pi0 = extract_pi0(fit),
        true_pi0 = mean(true_null),
        stringsAsFactors = FALSE
      )
      unit_rows[[result_index]] <- data.frame(
        seed = seed,
        permutation_seed = permutation_seed,
        arm = arm,
        fit_stage = fit_stage,
        unit_key = input$unit_keys,
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
        null_bf_rows[[arm]] <- cbind(
          seed = seed,
          permutation_seed = permutation_seed,
          summarize_null_bf(
            bf = fit$BF,
            true_null = true_null,
            arm = arm
          )
        )
      }
    }
  }

  alpha_curve <- do.call(rbind, alpha_rows)
  alpha_curve$seed <- seed
  alpha_curve$permutation_seed <- permutation_seed
  alpha_curve <- alpha_curve[, c(
    "seed", "permutation_seed", setdiff(
      names(alpha_curve), c("seed", "permutation_seed")
    )
  )]
  pi0_summary <- do.call(rbind, pi0_rows)
  null_bf_summary <- do.call(rbind, null_bf_rows)
  unit_results <- do.call(rbind, unit_rows)
  warnings <- unique(c(
    source_capture$warnings,
    permuted_capture$warnings,
    raw_capture$warnings,
    bf_capture$warnings
  ))
  validation <- data.frame(
    check = c(
      "alpha rows", "pi0 rows", "unit rows", "true pi0",
      "finite lfdr", "positive BF", "nonidentity permutation"
    ),
    pass = c(
      nrow(alpha_curve) == 4L * length(alpha_grid),
      nrow(pi0_summary) == 4L,
      nrow(unit_results) == 4L * 2L * n_alternatives,
      all(pi0_summary$true_pi0 == 0.5),
      all(is.finite(unit_results$lfdr)),
      all(
        is.na(unit_results$bayes_factor) |
          (is.finite(unit_results$bayes_factor) &
             unit_results$bayes_factor > 0)
      ),
      any(!shared_permutation$donor_map$fixed_point)
    ),
    stringsAsFactors = FALSE
  )
  if (!all(validation$pass)) {
    stop("A five-seed MC replicate failed internal validation.")
  }
  replicate <- list(
    seed = seed,
    permutation_seed = permutation_seed,
    alpha_curve = alpha_curve,
    pi0_summary = pi0_summary,
    null_bf_summary = null_bf_summary,
    unit_results = unit_results,
    donor_map = shared_permutation$donor_map,
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
all_donor_maps <- do.call(rbind, lapply(replicates, function(x) {
  cbind(
    seed = x$seed,
    permutation_seed = x$permutation_seed,
    x$donor_map
  )
}))
mc_alpha <- summarize_known_truth_mc(all_alpha)
mc_alpha005 <- mc_alpha[abs(mc_alpha$alpha - 0.05) < 1e-12, , drop = FALSE]
replicate_alpha005 <- all_alpha[
  abs(all_alpha$alpha - 0.05) < 1e-12, , drop = FALSE
]

expected_alpha_rows <- length(source_seeds) * 4L * length(alpha_grid)
if (nrow(all_alpha) != expected_alpha_rows ||
    any(mc_alpha$n_replications != length(source_seeds)) ||
    nrow(mc_alpha005) != 4L || nrow(replicate_alpha005) != 20L ||
    !setequal(unique(all_alpha$seed), source_seeds)) {
  stop("The combined five-seed Monte Carlo summaries are incomplete.")
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
  all_donor_maps,
  file.path(summary_directory, "replicate_donor_permutations.csv")
)
write_csv(mc_alpha, file.path(summary_directory, "mc_alpha_curve.csv"))
write_csv(mc_alpha005, file.path(summary_directory, "mc_alpha005_summary.csv"))
source_validation <- replicates[[as.character(source_seeds[1L])]]$source_reconstruction
write_csv(
  source_validation,
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

print(mc_alpha005)
cat("R1 known-truth five-seed permutation MC completed.\n")
