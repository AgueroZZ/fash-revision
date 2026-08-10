#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/internal/omega_precision_deep_sanity_check/omega_precision_deep_helpers.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/internal/omega_precision_deep_sanity_check/omega_precision_deep_helpers.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

get_argument <- function(name, default = NULL) {
  arguments <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_match <- which(startsWith(arguments, equals_prefix))
  if (length(equals_match) > 0L) {
    return(substring(arguments[[equals_match[[1]]]], nchar(equals_prefix) + 1L))
  }
  match_index <- which(arguments == name)
  if (length(match_index) == 0L || match_index[[1]] == length(arguments)) {
    return(default)
  }
  arguments[[match_index[[1]] + 1L]]
}

as_flag <- function(value) {
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

parse_integer_vector <- function(value, name) {
  result <- suppressWarnings(as.integer(strsplit(value, ",", fixed = TRUE)[[1]]))
  if (length(result) == 0L || anyNA(result) || any(result < 1L) || anyDuplicated(result)) {
    stop(name, " must contain unique comma-separated positive integers.")
  }
  sort(result)
}

write_csv <- function(value, path) {
  utils::write.csv(value, file = path, row.names = FALSE)
}

workflowr_root <- find_workflowr_root()
helper_path <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "omega_precision_deep_sanity_check",
  "omega_precision_deep_helpers.R"
)
runner_path <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "omega_precision_deep_sanity_check",
  "run_omega_precision_deep_sanity_check.R"
)
source(helper_path)

suppressPackageStartupMessages({
  library(fashr)
  library(Matrix)
})

default_output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "omega_precision_deep_sanity_check",
  "main"
)
output_dir <- get_argument("--output-dir", default_output_dir)
n_replications <- as.integer(get_argument("--n-replications", "200"))
seed_start <- as.integer(get_argument("--seed-start", "2026080701"))
num_cores <- as.integer(get_argument(
  "--num-cores",
  "1"
))
j_grid <- parse_integer_vector(
  get_argument("--j-grid", "100,300,1000,3000"),
  "--j-grid"
)
bootstrap_replications <- as.integer(get_argument("--bootstrap-replications", "2000"))
overwrite <- as_flag(get_argument("--overwrite", "false"))

if (anyNA(c(n_replications, seed_start, num_cores, bootstrap_replications)) ||
    n_replications < 1L || seed_start < 3L || num_cores < 1L ||
    bootstrap_replications < 100L) {
  stop("Replication, seed, core, and bootstrap arguments are invalid.")
}
if (length(j_grid) < 2L) {
  stop("--j-grid must contain at least two values for convergence slopes.")
}
if (num_cores != 1L) {
  stop(
    "This simulation must run with --num-cores 1. Forked BLAS/LAPACK calls ",
    "were not reliable on macOS during validation."
  )
}

package_path <- find.package("fashr")
package_artifact_paths <- list.files(
  package_path,
  pattern = "\\.(rdb|rdx|so|dylib)$",
  recursive = TRUE,
  full.names = TRUE,
  ignore.case = TRUE
)
package_artifact_md5 <- tools::md5sum(package_artifact_paths)
names(package_artifact_md5) <- sub(
  paste0("^", package_path, "/"),
  "",
  names(package_artifact_md5)
)

configuration <- list(
  n_time = 16L,
  x = 0:15,
  order = 1L,
  num_basis = 6L,
  betaprec = 0.8,
  pred_step = 1,
  grid = c(0, 0.45, 0.9),
  true_weights = c(0.50, 0.30, 0.20),
  error_sd = 0.65,
  rho = 0.75,
  j_grid = j_grid,
  alpha_grid = c(0.01, 0.025, 0.05, 0.10),
  n_replications = n_replications,
  seed_start = seed_start,
  num_cores = num_cores,
  bootstrap_replications = bootstrap_replications,
  rng_kind = RNGkind(),
  simulation_code_md5 = setNames(
    unname(tools::md5sum(c(runner_path, helper_path))),
    c("runner", "helper")
  ),
  package_artifact_md5 = package_artifact_md5,
  package_version = as.character(utils::packageVersion("fashr")),
  package_path = package_path
)
seed_list <- seed_start + seq_len(n_replications) - 1L
max_units <- max(j_grid)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
replicate_dir <- file.path(output_dir, "replicates")
dir.create(replicate_dir, recursive = TRUE, showWarnings = FALSE)

observation_covariance <- configuration$error_sd^2 * ar1_correlation(
  configuration$n_time,
  configuration$rho
)
diagonal_observation_covariance <- diag(diag(observation_covariance))
validate_symmetric_positive_definite(
  observation_covariance,
  "observation_covariance"
)
observation_precision <- solve(observation_covariance)
prior_covariances <- make_component_prior_covariances(
  x = configuration$x,
  grid = configuration$grid,
  num_basis = configuration$num_basis,
  betaprec = configuration$betaprec,
  order = configuration$order,
  pred_step = configuration$pred_step
)
full_component_covariances <- make_component_marginal_covariances(
  prior_covariances,
  observation_covariance
)
diagonal_component_covariances <- make_component_marginal_covariances(
  prior_covariances,
  diagonal_observation_covariance
)

validate_package_likelihood <- function() {
  simulated <- simulate_marginal_mixture(
    seed = configuration$seed_start - 1L,
    n_units = 4L,
    true_weights = configuration$true_weights,
    component_covariances = full_component_covariances
  )
  full_reference <- component_log_likelihood_matrix(
    simulated$response,
    full_component_covariances
  )
  diagonal_reference <- component_log_likelihood_matrix(
    simulated$response,
    diagonal_component_covariances
  )
  rows <- list()
  row_index <- 1L
  for (unit in seq_len(nrow(simulated$response))) {
    data_i <- data.frame(
      y = simulated$response[unit, ],
      x = configuration$x,
      offset = 0
    )
    for (component_index in seq_along(configuration$grid)) {
      psd <- configuration$grid[[component_index]]
      package_full <- fashr:::compute_L_gaussian_helper(
        data_i = data_i,
        Si = NULL,
        Omegai = observation_precision,
        psd_iwp = psd,
        num_basis = configuration$num_basis,
        betaprec = configuration$betaprec,
        order = configuration$order,
        pred_step = configuration$pred_step
      )
      package_diagonal <- fashr:::compute_L_gaussian_helper(
        data_i = data_i,
        Si = rep(configuration$error_sd, configuration$n_time),
        Omegai = NULL,
        psd_iwp = psd,
        num_basis = configuration$num_basis,
        betaprec = configuration$betaprec,
        order = configuration$order,
        pred_step = configuration$pred_step
      )
      rows[[row_index]] <- data.frame(
        unit = unit,
        psd = psd,
        observation_model = "Omega",
        package_value = package_full,
        vectorized_value = full_reference[unit, component_index],
        absolute_error = abs(package_full - full_reference[unit, component_index]),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        unit = unit,
        psd = psd,
        observation_model = "Diagonal",
        package_value = package_diagonal,
        vectorized_value = diagonal_reference[unit, component_index],
        absolute_error = abs(package_diagonal - diagonal_reference[unit, component_index]),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
  result <- do.call(rbind, rows)
  result$tolerance <- 1e-8
  result$passed <- result$absolute_error < result$tolerance
  if (!all(result$passed)) {
    stop("The vectorized likelihood failed package validation.")
  }
  result
}

validate_independent_identity <- function() {
  independent_observation_covariance <- configuration$error_sd^2 *
    diag(configuration$n_time)
  full_covariances <- make_component_marginal_covariances(
    prior_covariances,
    independent_observation_covariance
  )
  diagonal_covariances <- make_component_marginal_covariances(
    prior_covariances,
    diag(diag(independent_observation_covariance))
  )
  maximum_covariance_difference <- max(vapply(
    seq_along(full_covariances),
    function(index) max(abs(full_covariances[[index]] - diagonal_covariances[[index]])),
    numeric(1)
  ))
  simulated <- simulate_marginal_mixture(
    seed = configuration$seed_start - 2L,
    n_units = min(max_units, 300L),
    true_weights = configuration$true_weights,
    component_covariances = full_covariances
  )
  full_likelihood <- component_log_likelihood_matrix(
    simulated$response,
    full_covariances
  )
  diagonal_likelihood <- component_log_likelihood_matrix(
    simulated$response,
    diagonal_covariances
  )
  maximum_likelihood_difference <- max(abs(full_likelihood - diagonal_likelihood))
  full_raw <- fit_raw_from_likelihood(full_likelihood, configuration$grid, penalty = 1)
  diagonal_raw <- fit_raw_from_likelihood(diagonal_likelihood, configuration$grid, penalty = 1)
  maximum_weight_difference <- max(abs(
    expand_prior_weights(full_raw, configuration$grid) -
      expand_prior_weights(diagonal_raw, configuration$grid)
  ))
  maximum_lfdr_difference <- max(abs(full_raw$lfdr - diagonal_raw$lfdr))
  result <- data.frame(
    maximum_covariance_difference = maximum_covariance_difference,
    maximum_likelihood_difference = maximum_likelihood_difference,
    maximum_weight_difference = maximum_weight_difference,
    maximum_lfdr_difference = maximum_lfdr_difference,
    tolerance = 1e-10,
    passed = max(
      maximum_covariance_difference,
      maximum_likelihood_difference,
      maximum_weight_difference,
      maximum_lfdr_difference
    ) < 1e-10,
    stringsAsFactors = FALSE
  )
  if (!result$passed) {
    stop("The rho = 0 full/diagonal identity control failed.")
  }
  result
}

package_likelihood_validation <- validate_package_likelihood()
independent_identity_check <- validate_independent_identity()

make_replication <- function(seed) {
  simulated <- simulate_marginal_mixture(
    seed = seed,
    n_units = max_units,
    true_weights = configuration$true_weights,
    component_covariances = full_component_covariances
  )
  full_likelihood <- component_log_likelihood_matrix(
    simulated$response,
    full_component_covariances
  )
  diagonal_likelihood <- component_log_likelihood_matrix(
    simulated$response,
    diagonal_component_covariances
  )

  prior_parts <- list()
  fdr_parts <- list()
  bf_parts <- list()
  prior_index <- 1L
  fdr_index <- 1L
  bf_index <- 1L

  for (n_units in configuration$j_grid) {
    selected <- seq_len(n_units)
    truth <- simulated$component[selected]
    likelihoods <- list(
      Omega = full_likelihood[selected, , drop = FALSE],
      Diagonal = diagonal_likelihood[selected, , drop = FALSE]
    )

    for (observation_model in names(likelihoods)) {
      likelihood <- likelihoods[[observation_model]]
      raw_fit <- fit_raw_from_likelihood(
        likelihood,
        configuration$grid,
        penalty = 1
      )
      bf_fit <- fit_bf_from_raw(raw_fit, configuration$grid)
      oracle_fit <- fit_oracle_from_likelihood(
        likelihood,
        configuration$grid,
        configuration$true_weights
      )
      fit_list <- list(Raw = raw_fit, BF = bf_fit, Oracle = oracle_fit)

      for (correction in names(fit_list)) {
        fit <- fit_list[[correction]]
        prior_parts[[prior_index]] <- prior_weight_rows(
          fit = fit,
          truth_component = truth,
          true_weights = configuration$true_weights,
          grid = configuration$grid,
          seed = seed,
          n_units = n_units,
          observation_model = observation_model,
          correction = correction
        )
        prior_index <- prior_index + 1L
        fdr_parts[[fdr_index]] <- evaluate_fdr_fit(
          fit = fit,
          truth_component = truth,
          alpha_grid = configuration$alpha_grid,
          seed = seed,
          n_units = n_units,
          observation_model = observation_model,
          correction = correction
        )
        fdr_index <- fdr_index + 1L
      }

      bf_parts[[bf_index]] <- bf_diagnostic_row(
        fit = bf_fit,
        true_weights = configuration$true_weights,
        grid = configuration$grid,
        seed = seed,
        n_units = n_units,
        observation_model = observation_model
      )
      bf_index <- bf_index + 1L
    }
  }

  list(
    configuration = configuration,
    seed = seed,
    prior_weights = do.call(rbind, prior_parts),
    fdr = do.call(rbind, fdr_parts),
    bf = do.call(rbind, bf_parts)
  )
}

validate_replication <- function(replication, seed) {
  required_fields <- c("configuration", "seed", "prior_weights", "fdr", "bf")
  if (!all(required_fields %in% names(replication)) ||
      !identical(replication$seed, seed) ||
      !isTRUE(all.equal(replication$configuration, configuration))) {
    return(FALSE)
  }
  expected_prior_rows <- length(configuration$j_grid) * 2L * 3L *
    length(configuration$grid)
  expected_fdr_rows <- length(configuration$j_grid) * 2L * 3L *
    length(configuration$alpha_grid)
  expected_bf_rows <- length(configuration$j_grid) * 2L
  nrow(replication$prior_weights) == expected_prior_rows &&
    nrow(replication$fdr) == expected_fdr_rows &&
    nrow(replication$bf) == expected_bf_rows &&
    !anyDuplicated(replication$prior_weights[
      c("seed", "n_units", "observation_model", "correction", "psd")
    ]) &&
    !anyDuplicated(replication$fdr[
      c("seed", "n_units", "observation_model", "correction", "alpha")
    ]) &&
    !anyDuplicated(replication$bf[
      c("seed", "n_units", "observation_model")
    ]) &&
    all(is.finite(replication$prior_weights$estimated_weight)) &&
    all(replication$prior_weights$estimated_weight >= 0) &&
    all(replication$prior_weights$estimated_weight <= 1) &&
    all(is.finite(replication$fdr$fdp)) &&
    all(replication$fdr$fdp >= 0) &&
    all(replication$fdr$fdp <= 1) &&
    all(is.finite(replication$bf$estimated_pi0))
}

run_or_reuse_replication <- function(seed) {
  cache_path <- file.path(replicate_dir, paste0("seed_", seed, ".rds"))
  if (file.exists(cache_path) && !overwrite) {
    cached <- readRDS(cache_path)
    if (!validate_replication(cached, seed)) {
      stop("Cached replication does not match the current configuration: ", cache_path)
    }
    return(cached)
  }
  replication <- make_replication(seed)
  if (!validate_replication(replication, seed)) {
    stop("New replication failed validation for seed ", seed)
  }
  temporary_cache_path <- tempfile(
    pattern = paste0("seed_", seed, "_"),
    tmpdir = replicate_dir,
    fileext = ".rds.tmp"
  )
  on.exit(unlink(temporary_cache_path), add = TRUE)
  saveRDS(replication, temporary_cache_path)
  if (!file.rename(temporary_cache_path, cache_path)) {
    stop("Could not atomically install replication cache: ", cache_path)
  }
  replication
}

message(
  "Running ", n_replications, " paired replications with Jmax = ", max_units,
  " on ", num_cores, " core(s)."
)
replications <- lapply(seed_list, run_or_reuse_replication)
if (!all(vapply(
  seq_along(replications),
  function(index) validate_replication(replications[[index]], seed_list[[index]]),
  logical(1)
))) {
  stop("At least one returned replication failed validation.")
}

prior_weights_long <- do.call(rbind, lapply(replications, `[[`, "prior_weights"))
fdr_replications <- do.call(rbind, lapply(replications, `[[`, "fdr"))
bf_diagnostics <- do.call(rbind, lapply(replications, `[[`, "bf"))
rownames(prior_weights_long) <- NULL
rownames(fdr_replications) <- NULL
rownames(bf_diagnostics) <- NULL

prior_results <- summarize_prior_accuracy(prior_weights_long)
prior_accuracy_summary <- prior_results$summary
prior_vector_by_seed <- prior_results$vector_by_seed
fdr_summary <- summarize_fdr_replications(fdr_replications)
bf_summary <- summarize_bf_diagnostics(bf_diagnostics)

make_convergence_slopes <- function(vector_rows) {
  vector_rows <- vector_rows[vector_rows$correction != "Oracle", , drop = FALSE]
  key <- interaction(
    vector_rows$seed,
    vector_rows$observation_model,
    vector_rows$correction,
    drop = TRUE
  )
  slope_rows <- lapply(split(vector_rows, key), function(group) {
    fit <- stats::lm(log(pmax(group$tv, 1e-12)) ~ log(group$n_units))
    data.frame(
      seed = group$seed[[1]],
      observation_model = group$observation_model[[1]],
      correction = group$correction[[1]],
      slope = unname(stats::coef(fit)[[2]]),
      stringsAsFactors = FALSE
    )
  })
  slope_rows <- do.call(rbind, slope_rows)
  summary_key <- interaction(
    slope_rows$observation_model,
    slope_rows$correction,
    drop = TRUE
  )
  summary_rows <- lapply(seq_along(split(slope_rows, summary_key)), function(index) {
    group <- split(slope_rows, summary_key)[[index]]
    interval <- bootstrap_mean_ci(
      group$slope,
      n_bootstrap = bootstrap_replications,
      seed = configuration$seed_start + 8000L + index
    )
    data.frame(
      observation_model = group$observation_model[[1]],
      correction = group$correction[[1]],
      n_replications = nrow(group),
      mean_log_tv_slope = interval[["mean"]],
      bootstrap_ci_lower = interval[["lower"]],
      bootstrap_ci_upper = interval[["upper"]],
      stringsAsFactors = FALSE
    )
  })
  list(by_seed = slope_rows, summary = do.call(rbind, summary_rows))
}

make_paired_fdr_differences <- function(fdr_rows) {
  omega <- fdr_rows[fdr_rows$observation_model == "Omega", , drop = FALSE]
  diagonal <- fdr_rows[fdr_rows$observation_model == "Diagonal", , drop = FALSE]
  paired <- merge(
    omega,
    diagonal,
    by = c("seed", "n_units", "correction", "alpha"),
    suffixes = c("_omega", "_diagonal")
  )
  paired$delta_fdp_diagonal_minus_omega <- paired$fdp_diagonal - paired$fdp_omega
  paired$delta_power_diagonal_minus_omega <- paired$power_diagonal - paired$power_omega
  key <- interaction(paired$correction, paired$n_units, paired$alpha, drop = TRUE)
  groups <- split(paired, key)
  rows <- lapply(seq_along(groups), function(index) {
    group <- groups[[index]]
    fdr_interval <- bootstrap_mean_ci(
      group$delta_fdp_diagonal_minus_omega,
      n_bootstrap = bootstrap_replications,
      seed = configuration$seed_start + 9000L + index
    )
    power_interval <- bootstrap_mean_ci(
      group$delta_power_diagonal_minus_omega,
      n_bootstrap = bootstrap_replications,
      seed = configuration$seed_start + 10000L + index
    )
    data.frame(
      correction = group$correction[[1]],
      n_units = group$n_units[[1]],
      alpha = group$alpha[[1]],
      n_replications = nrow(group),
      mean_delta_fdr_diagonal_minus_omega = fdr_interval[["mean"]],
      fdr_bootstrap_ci_lower = fdr_interval[["lower"]],
      fdr_bootstrap_ci_upper = fdr_interval[["upper"]],
      mean_delta_power_diagonal_minus_omega = power_interval[["mean"]],
      power_bootstrap_ci_lower = power_interval[["lower"]],
      power_bootstrap_ci_upper = power_interval[["upper"]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_paired_prior_differences <- function(vector_rows) {
  vector_rows <- vector_rows[vector_rows$correction != "Oracle", , drop = FALSE]
  omega <- vector_rows[vector_rows$observation_model == "Omega", , drop = FALSE]
  diagonal <- vector_rows[vector_rows$observation_model == "Diagonal", , drop = FALSE]
  paired <- merge(
    omega,
    diagonal,
    by = c("seed", "n_units", "correction"),
    suffixes = c("_omega", "_diagonal")
  )
  paired$delta_tv_diagonal_minus_omega <- paired$tv_diagonal - paired$tv_omega
  key <- interaction(paired$correction, paired$n_units, drop = TRUE)
  groups <- split(paired, key)
  rows <- lapply(seq_along(groups), function(index) {
    group <- groups[[index]]
    interval <- bootstrap_mean_ci(
      group$delta_tv_diagonal_minus_omega,
      n_bootstrap = bootstrap_replications,
      seed = configuration$seed_start + 11000L + index
    )
    data.frame(
      correction = group$correction[[1]],
      n_units = group$n_units[[1]],
      n_replications = nrow(group),
      mean_delta_tv_diagonal_minus_omega = interval[["mean"]],
      bootstrap_ci_lower = interval[["lower"]],
      bootstrap_ci_upper = interval[["upper"]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_paired_correction_differences <- function(fdr_rows) {
  raw <- fdr_rows[fdr_rows$correction == "Raw", , drop = FALSE]
  bf <- fdr_rows[fdr_rows$correction == "BF", , drop = FALSE]
  paired <- merge(
    raw,
    bf,
    by = c("seed", "n_units", "observation_model", "alpha"),
    suffixes = c("_raw", "_bf")
  )
  paired$delta_fdp_bf_minus_raw <- paired$fdp_bf - paired$fdp_raw
  paired$delta_power_bf_minus_raw <- paired$power_bf - paired$power_raw
  key <- interaction(
    paired$observation_model,
    paired$n_units,
    paired$alpha,
    drop = TRUE
  )
  groups <- split(paired, key)
  rows <- lapply(seq_along(groups), function(index) {
    group <- groups[[index]]
    fdr_interval <- bootstrap_mean_ci(
      group$delta_fdp_bf_minus_raw,
      n_bootstrap = bootstrap_replications,
      seed = configuration$seed_start + 12000L + index
    )
    power_interval <- bootstrap_mean_ci(
      group$delta_power_bf_minus_raw,
      n_bootstrap = bootstrap_replications,
      seed = configuration$seed_start + 13000L + index
    )
    data.frame(
      observation_model = group$observation_model[[1]],
      n_units = group$n_units[[1]],
      alpha = group$alpha[[1]],
      n_replications = nrow(group),
      mean_delta_fdr_bf_minus_raw = fdr_interval[["mean"]],
      fdr_bootstrap_ci_lower = fdr_interval[["lower"]],
      fdr_bootstrap_ci_upper = fdr_interval[["upper"]],
      mean_delta_power_bf_minus_raw = power_interval[["mean"]],
      power_bootstrap_ci_lower = power_interval[["lower"]],
      power_bootstrap_ci_upper = power_interval[["upper"]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

convergence_slopes <- make_convergence_slopes(prior_vector_by_seed)
paired_fdr_differences <- make_paired_fdr_differences(fdr_replications)
paired_prior_differences <- make_paired_prior_differences(prior_vector_by_seed)
paired_correction_differences <- make_paired_correction_differences(
  fdr_replications
)

realized_prior_rows <- prior_weights_long
realized_prior_rows$population_error <- realized_prior_rows$realized_error
realized_prior_results <- summarize_prior_accuracy(realized_prior_rows)
prior_realized_accuracy_summary <- realized_prior_results$summary
prior_realized_accuracy_summary$target <- "realized_prefix_proportion"

expected_prior_rows <- n_replications * length(j_grid) * 2L * 3L *
  length(configuration$grid)
expected_fdr_rows <- n_replications * length(j_grid) * 2L * 3L *
  length(configuration$alpha_grid)
expected_bf_rows <- n_replications * length(j_grid) * 2L
if (nrow(prior_weights_long) != expected_prior_rows ||
    nrow(fdr_replications) != expected_fdr_rows ||
    nrow(bf_diagnostics) != expected_bf_rows ||
    length(unique(prior_weights_long$seed)) != n_replications) {
  stop("Combined Monte Carlo output has incomplete coverage.")
}

write_csv(package_likelihood_validation, file.path(output_dir, "package_likelihood_validation.csv"))
write_csv(independent_identity_check, file.path(output_dir, "independent_identity_check.csv"))
write_csv(prior_weights_long, file.path(output_dir, "prior_weights_long.csv"))
write_csv(prior_vector_by_seed, file.path(output_dir, "prior_vector_accuracy_by_seed.csv"))
write_csv(prior_accuracy_summary, file.path(output_dir, "prior_accuracy_summary.csv"))
write_csv(
  prior_realized_accuracy_summary,
  file.path(output_dir, "prior_realized_accuracy_summary.csv")
)
write_csv(bf_diagnostics, file.path(output_dir, "bf_diagnostics.csv"))
write_csv(bf_summary, file.path(output_dir, "bf_summary.csv"))
write_csv(fdr_replications, file.path(output_dir, "fdr_replications.csv"))
write_csv(fdr_summary, file.path(output_dir, "fdr_summary.csv"))
write_csv(convergence_slopes$by_seed, file.path(output_dir, "convergence_slopes_by_seed.csv"))
write_csv(convergence_slopes$summary, file.path(output_dir, "convergence_slope_summary.csv"))
write_csv(paired_fdr_differences, file.path(output_dir, "paired_fdr_differences.csv"))
write_csv(paired_prior_differences, file.path(output_dir, "paired_prior_differences.csv"))
write_csv(
  paired_correction_differences,
  file.path(output_dir, "paired_bf_raw_differences.csv")
)
saveRDS(configuration, file.path(output_dir, "configuration.rds"))
writeLines(capture.output(str(configuration)), file.path(output_dir, "configuration.txt"))
writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "session_info.txt"))

plot_prior_convergence <- function(summary_table, output_path) {
  plot_data <- summary_table[
    summary_table$summary_level == "vector" &
      summary_table$correction %in% c("Raw", "BF"),
    ,
    drop = FALSE
  ]
  method_levels <- c("Omega-Raw", "Omega-BF", "Diagonal-Raw", "Diagonal-BF")
  colors <- c("#0072B2", "#56B4E9", "#D55E00", "#E69F00")
  line_types <- c(1, 2, 1, 2)
  grDevices::png(output_path, width = 1500, height = 1000, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(
    range(plot_data$n_units),
    range(c(plot_data$tv_ci_lower, plot_data$tv_ci_upper), finite = TRUE),
    type = "n",
    log = "x",
    xlab = "Number of units (J, log scale)",
    ylab = "Mean total-variation error in prior weights",
    main = "Prior-weight learning under correct and misspecified error models"
  )
  for (method_index in seq_along(method_levels)) {
    method_data <- plot_data[plot_data$method == method_levels[[method_index]], ]
    method_data <- method_data[order(method_data$n_units), ]
    graphics::lines(
      method_data$n_units,
      method_data$mean_tv,
      type = "b",
      pch = 19,
      col = colors[[method_index]],
      lty = line_types[[method_index]],
      lwd = 2
    )
    graphics::segments(
      method_data$n_units,
      method_data$tv_ci_lower,
      method_data$n_units,
      method_data$tv_ci_upper,
      col = colors[[method_index]]
    )
  }
  graphics::legend(
    "topright",
    legend = method_levels,
    col = colors,
    lty = line_types,
    pch = 19,
    lwd = 2,
    bty = "n"
  )
}

plot_fdr_calibration <- function(summary_table, output_path) {
  plot_data <- summary_table[summary_table$n_units == max(summary_table$n_units), ]
  method_levels <- c(
    "Omega-Oracle", "Omega-Raw", "Omega-BF",
    "Diagonal-Oracle", "Diagonal-Raw", "Diagonal-BF"
  )
  colors <- c("#009E73", "#0072B2", "#56B4E9", "#CC79A7", "#D55E00", "#E69F00")
  line_types <- c(1, 1, 2, 3, 1, 2)
  y_max <- max(c(plot_data$fdr_ci_upper, plot_data$alpha), na.rm = TRUE)
  grDevices::png(output_path, width = 1500, height = 1000, res = 160)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(
    range(plot_data$alpha),
    c(0, y_max * 1.05),
    type = "n",
    xlab = "Nominal alpha",
    ylab = "Monte Carlo FDR (mean FDP)",
    main = paste0("FDR calibration at J = ", max(plot_data$n_units))
  )
  graphics::abline(a = 0, b = 1, col = "gray40", lty = 3, lwd = 2)
  for (method_index in seq_along(method_levels)) {
    method_data <- plot_data[plot_data$method == method_levels[[method_index]], ]
    method_data <- method_data[order(method_data$alpha), ]
    graphics::lines(
      method_data$alpha,
      method_data$mean_fdr,
      type = "b",
      pch = 19,
      col = colors[[method_index]],
      lty = line_types[[method_index]],
      lwd = 2
    )
  }
  graphics::legend(
    "topleft",
    legend = method_levels,
    col = colors,
    lty = line_types,
    pch = 19,
    lwd = 2,
    bty = "n"
  )
}

plot_prior_convergence(
  prior_accuracy_summary,
  file.path(output_dir, "prior_weight_convergence.png")
)
plot_fdr_calibration(
  fdr_summary,
  file.path(output_dir, "fdr_calibration.png")
)

saveRDS(
  list(
    configuration = configuration,
    package_likelihood_validation = package_likelihood_validation,
    independent_identity_check = independent_identity_check,
    prior_accuracy_summary = prior_accuracy_summary,
    bf_summary = bf_summary,
    fdr_summary = fdr_summary,
    convergence_slope_summary = convergence_slopes$summary,
    paired_fdr_differences = paired_fdr_differences,
    paired_prior_differences = paired_prior_differences,
    paired_correction_differences = paired_correction_differences,
    prior_realized_accuracy_summary = prior_realized_accuracy_summary
  ),
  file.path(output_dir, "omega_precision_deep_sanity_summary.rds")
)

cat("Completed Omega deep sanity check.\n")
cat("Output directory:", normalizePath(output_dir), "\n")
cat("Replications:", n_replications, "\n")
cat("Package likelihood maximum error:", max(package_likelihood_validation$absolute_error), "\n")
print(prior_accuracy_summary[
  prior_accuracy_summary$summary_level == "vector" &
    prior_accuracy_summary$n_units %in% c(min(j_grid), max(j_grid)),
  c("method", "n_units", "mean_tv", "tv_ci_lower", "tv_ci_upper")
])
print(fdr_summary[
  fdr_summary$alpha == 0.05 & fdr_summary$n_units == max(j_grid),
  c("method", "mean_fdr", "fdr_ci_lower", "fdr_ci_upper", "mean_power", "mean_discoveries")
])
