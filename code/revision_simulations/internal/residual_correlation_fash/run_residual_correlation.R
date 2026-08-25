#!/usr/bin/env Rscript

# Residual-based common-correlation FASH experiment.
#
# Stages, each runnable on its own so the expensive ones can be checkpointed:
#
#   A  extract the fixed 6,362-unit subset, refit the diagonal mixture, BF update
#   B  closed-form BMA posterior means, residuals, residual z-scores
#   C  the common correlation matrix and its comparators
#   D  independence calibration by parametric bootstrap
#   E  dependent-likelihood refit driven by the estimated correlation
#   F  exact analytic calibration, the cheap companion to stage D
#   G  linear deconvolution of the shrinkage operator, and its conditioning
#
# Usage:
#   Rscript --vanilla run_residual_correlation.R --stage A
#   Rscript --vanilla run_residual_correlation.R --stage D --replicates 1

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  arguments <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(arguments, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(arguments[equals_hit[1L]], nchar(equals_prefix) + 1L))
  }
  hit <- which(arguments == name)
  if (length(hit) == 0L || hit[1L] == length(arguments)) {
    return(default)
  }
  arguments[hit[1L] + 1L]
}

log_message <- function(...) {
  message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
}

load_exact_object <- function(path, expected_name) {
  object_environment <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = object_environment)
  if (!identical(loaded_names, expected_name)) {
    stop(path, " must contain only ", expected_name, ".")
  }
  object_environment[[expected_name]]
}

file_metadata <- function(path) {
  information <- file.info(path)
  if (nrow(information) != 1L || is.na(information$size)) {
    stop("Could not read source-file metadata: ", path)
  }
  list(
    path = normalizePath(path, mustWork = TRUE),
    size_bytes = unname(information$size),
    modification_time = format(information$mtime, tz = "UTC", usetz = TRUE)
  )
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

matrix_to_long <- function(matrix, matrix_id, matrix_label) {
  matrix <- as.matrix(matrix)
  grid <- expand.grid(
    time_a = seq_len(nrow(matrix)) - 1L,
    time_b = seq_len(ncol(matrix)) - 1L
  )
  grid$correlation <- matrix[cbind(grid$time_a + 1L, grid$time_b + 1L)]
  grid$matrix_id <- matrix_id
  grid$matrix_label <- matrix_label
  grid[, c("matrix_id", "matrix_label", "time_a", "time_b", "correlation")]
}

workflowr_root <- find_workflowr_root()
# Order matters: the shared helpers are sourced first so that this
# experiment's own definitions take precedence on any name they share.
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "one_variant_per_gene_refit", "one_variant_per_gene_refit_helpers.R"
))
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "correlated_likelihood_sensitivity", "correlated_likelihood_helpers.R"
))
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "residual_correlation_fash", "residual_correlation_helpers.R"
))

stage <- toupper(as.character(get_arg("--stage", "A")))
num_cores <- as.integer(get_arg("--cores", "8"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
n_replicates <- as.integer(get_arg("--replicates", "1"))
correlation_method <- as.character(get_arg("--correlation", "pearson"))
weight_source <- as.character(get_arg("--weights", "bf"))
# Which correlation matrix Stage E should feed into Omega. "residual" is the
# naive matrix straight out of Stage C; "deconvolved" is the Stage G inversion
# that removes the shrinkage artifact and is the scientifically meaningful one.
matrix_choice <- as.character(get_arg("--matrix", "residual"))
if (!stage %in% c("A", "B", "C", "D", "E", "F", "G") || is.na(num_cores) ||
    num_cores < 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
    is.na(n_replicates) || n_replicates < 1L ||
    !weight_source %in% c("bf", "raw") ||
    !matrix_choice %in% c("residual", "deconvolved", "deconvolved_ridge",
                          "permutation")) {
  stop("Invalid command-line arguments.")
}

cache_directory <- file.path(
  workflowr_root, "output", "revision_simulations", "internal",
  "residual_correlation_fash"
)
summary_directory <- file.path(cache_directory, "summary")
dir.create(summary_directory, recursive = TRUE, showWarnings = FALSE)

pair_metadata_path <- file.path(
  workflowr_root, "output", "revision_simulations", "internal",
  "correlated_likelihood_sensitivity", "summary", "pair_metadata.csv"
)
raw_fit_path <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real", "fash_fit1_all.RData"
)

# ---------------------------------------------------------------- Stage A ----

if (stage == "A") {
  log_message("Stage A: extracting the fixed one-variant-per-gene subset.")
  pair_metadata <- utils::read.csv(pair_metadata_path, stringsAsFactors = FALSE)
  if (!all(c("seed", "fash_index", "pair_key", "gene_id", "variant_id") %in%
             names(pair_metadata)) ||
      length(unique(pair_metadata$seed)) != 1L ||
      anyDuplicated(pair_metadata$pair_key) ||
      anyDuplicated(pair_metadata$gene_id) ||
      anyDuplicated(pair_metadata$fash_index)) {
    stop("The cached pair metadata is not a valid one-variant-per-gene subset.")
  }
  log_message("  reusing ", nrow(pair_metadata), " pairs, thinning seed ",
              unique(pair_metadata$seed), ".")

  log_message("  loading the full raw FASH(1) fit (this needs several GB).")
  full_raw <- load_exact_object(raw_fit_path, "fash_fit1")
  settings <- validate_fash_settings(full_raw$settings)
  psd_grid <- as.numeric(full_raw$psd_grid)
  full_pair_keys <- names(full_raw$fash_data$data_list)
  selected_indices <- as.integer(pair_metadata$fash_index)
  if (!identical(full_pair_keys[selected_indices], pair_metadata$pair_key)) {
    stop("The cached selection is not aligned with the source FASH fit.")
  }
  log_message("  full fit: ", length(full_pair_keys), " pairs, ",
              length(psd_grid), " grid points.")

  beta_hat <- t(vapply(
    full_raw$fash_data$data_list[selected_indices],
    function(dataset) as.numeric(dataset$y), numeric(16)
  ))
  adjusted_se <- t(vapply(
    full_raw$fash_data$S[selected_indices],
    function(se) as.numeric(se), numeric(16)
  ))
  time_grid <- as.numeric(
    full_raw$fash_data$data_list[[selected_indices[1]]]$x
  )
  offsets <- vapply(
    full_raw$fash_data$data_list[selected_indices],
    function(dataset) {
      if (is.null(dataset$offset)) 0 else max(abs(as.numeric(dataset$offset)))
    }, numeric(1)
  )
  if (max(offsets) > 0) {
    stop("This experiment assumes a zero offset in the source fit.")
  }
  rownames(beta_hat) <- pair_metadata$pair_key
  rownames(adjusted_se) <- pair_metadata$pair_key
  colnames(beta_hat) <- paste0("t", time_grid)
  colnames(adjusted_se) <- paste0("t", time_grid)

  log_message("  re-estimating the diagonal mixture on the subset.")
  raw_fit <- refit_fash_from_cached_likelihood(
    full_raw, selected_indices, penalty = settings$penalty
  )
  source_metadata <- file_metadata(raw_fit_path)
  full_prior <- full_raw$prior_weights
  rm(full_raw)
  invisible(gc())

  log_message("  applying BF_update().")
  bf_fit <- fashr::BF_update(raw_fit)
  names(bf_fit$lfdr) <- pair_metadata$pair_key
  rownames(bf_fit$posterior_weights) <- pair_metadata$pair_key

  saveRDS(
    list(
      pair_metadata = pair_metadata,
      beta_hat = beta_hat,
      adjusted_se = adjusted_se,
      time_grid = time_grid,
      settings = settings,
      psd_grid = psd_grid,
      full_prior_weights = full_prior,
      source_metadata = source_metadata
    ),
    file.path(cache_directory, "subset_inputs.rds")
  )
  saveRDS(
    list(raw = raw_fit, bf = bf_fit),
    file.path(cache_directory, "diagonal_fits.rds")
  )

  summary_table <- data.frame(
    stage = c("raw", "bf"),
    pi0 = c(
      raw_fit$prior_weights$prior_weight[raw_fit$prior_weights$psd == 0],
      bf_fit$prior_weights$prior_weight[bf_fit$prior_weights$psd == 0]
    ),
    mean_lfdr = c(mean(raw_fit$lfdr), mean(bf_fit$lfdr)),
    calls = c(
      length(cumulative_lfdr_calls(raw_fit$lfdr, alpha)),
      length(cumulative_lfdr_calls(bf_fit$lfdr, alpha))
    ),
    retained_components = c(
      nrow(raw_fit$prior_weights), nrow(bf_fit$prior_weights)
    ),
    stringsAsFactors = FALSE
  )
  write_csv(summary_table, file.path(summary_directory, "stage_a_diagonal.csv"))
  print(summary_table)
  log_message("Stage A complete.")
}

# ---------------------------------------------------------------- Stage B ----

compute_residuals <- function(inputs, fit, num_cores) {
  shared_design <- build_shared_design(
    time_grid = inputs$time_grid,
    num_basis = inputs$settings$num_basis,
    order = inputs$settings$order,
    betaprec = inputs$settings$betaprec
  )
  retained <- retained_psd_values(fit)
  scales <- psd_to_prior_scale(
    retained, inputs$settings$order, inputs$settings$pred_step
  )$prior_scale
  posterior_means <- bma_posterior_means(
    beta_hat = inputs$beta_hat,
    standard_errors = inputs$adjusted_se,
    posterior_weights = fit$posterior_weights,
    shared_design = shared_design,
    prior_scale = scales,
    num_cores = num_cores
  )
  residual_z_scores(inputs$beta_hat, posterior_means, inputs$adjusted_se)
}

if (stage == "B") {
  log_message("Stage B: closed-form posterior means and residual z-scores.")
  inputs <- readRDS(file.path(cache_directory, "subset_inputs.rds"))
  fits <- readRDS(file.path(cache_directory, "diagonal_fits.rds"))

  results <- list()
  for (name in c("raw", "bf")) {
    started <- Sys.time()
    results[[name]] <- compute_residuals(inputs, fits[[name]], num_cores)
    log_message("  ", name, " posterior means in ",
                format(difftime(Sys.time(), started, units = "secs"),
                       digits = 4), ".")
  }

  # Validate the closed form against the package's own Monte Carlo predictor
  # on a fixed sample of real units.
  set.seed(20260820)
  check_indices <- sort(sample(nrow(inputs$beta_hat), 12))
  monte_carlo <- t(vapply(check_indices, function(index) {
    stats::predict(fits$bf, index = index, M = 6000)$mean
  }, numeric(length(inputs$time_grid))))
  closed_form <- inputs$beta_hat[check_indices, , drop = FALSE] -
    results$bf$residual[check_indices, , drop = FALSE]
  validation <- data.frame(
    n_units = length(check_indices),
    monte_carlo_draws = 6000,
    maximum_absolute_difference = max(abs(closed_form - monte_carlo)),
    mean_absolute_difference = mean(abs(closed_form - monte_carlo)),
    typical_posterior_sd = mean(inputs$adjusted_se[check_indices, ]),
    stringsAsFactors = FALSE
  )
  write_csv(validation,
            file.path(summary_directory, "stage_b_predict_validation.csv"))
  print(validation)

  saveRDS(results, file.path(cache_directory, "residuals.rds"))

  scale_summary <- do.call(rbind, lapply(names(results), function(name) {
    z <- results[[name]]$z
    data.frame(
      weights = name,
      mean_z = mean(z),
      sd_z = stats::sd(z),
      mean_column_sd = mean(apply(z, 2, stats::sd)),
      minimum_column_sd = min(apply(z, 2, stats::sd)),
      maximum_column_sd = max(apply(z, 2, stats::sd)),
      excess_kurtosis = mean((z - mean(z))^4) / stats::var(as.numeric(z))^2 - 3,
      stringsAsFactors = FALSE
    )
  }))
  write_csv(scale_summary, file.path(summary_directory, "stage_b_z_scale.csv"))
  print(scale_summary)
  log_message("Stage B complete.")
}

# ---------------------------------------------------------------- Stage C ----

if (stage == "C") {
  log_message("Stage C: the common correlation matrix.")
  inputs <- readRDS(file.path(cache_directory, "subset_inputs.rds"))
  residuals <- readRDS(file.path(cache_directory, "residuals.rds"))

  matrices <- list()
  for (weights in c("bf", "raw")) {
    for (method in c("pearson", "uncentred", "spearman")) {
      matrices[[paste0(weights, "_", method)]] <-
        estimate_common_correlation(residuals[[weights]]$z, method)
    }
  }
  # Raw observed z on the same units, for contrast with the residual version.
  raw_z <- inputs$beta_hat / inputs$adjusted_se
  matrices[["observed_raw_z"]] <- estimate_common_correlation(raw_z, "pearson")

  # Comparators from the earlier experiments on the identical subset.
  comparator_path <- file.path(
    workflowr_root, "output", "revision_simulations", "internal",
    "correlated_likelihood_sensitivity", "summary",
    "correlation_matrices_long.csv"
  )
  if (file.exists(comparator_path)) {
    comparator_long <- utils::read.csv(comparator_path,
                                       stringsAsFactors = FALSE)
    for (id in unique(comparator_long$matrix_id)) {
      subset_rows <- comparator_long[comparator_long$matrix_id == id, ]
      n_time <- max(subset_rows$time_a) + 1L
      matrix_form <- matrix(0, n_time, n_time)
      matrix_form[cbind(subset_rows$time_a + 1L, subset_rows$time_b + 1L)] <-
        subset_rows$correlation
      matrices[[id]] <- matrix_form
    }
  }

  diagnostics <- do.call(rbind, lapply(names(matrices), function(name) {
    correlation_diagnostics(matrices[[name]], name)
  }))
  lag_profiles <- do.call(rbind, lapply(names(matrices), function(name) {
    correlation_lag_profile(matrices[[name]], name)
  }))
  primary <- matrices[[paste0(weight_source, "_", correlation_method)]]
  comparisons <- do.call(rbind, lapply(
    setdiff(names(matrices), paste0(weight_source, "_", correlation_method)),
    function(name) {
      compare_correlation_matrices(
        primary, matrices[[name]],
        paste0(weight_source, "_", correlation_method), name
      )
    }
  ))
  long_form <- do.call(rbind, lapply(names(matrices), function(name) {
    matrix_to_long(matrices[[name]], name, name)
  }))

  saveRDS(matrices, file.path(cache_directory, "correlation_matrices.rds"))
  write_csv(diagnostics, file.path(summary_directory,
                                   "stage_c_diagnostics.csv"))
  write_csv(lag_profiles, file.path(summary_directory,
                                    "stage_c_lag_profiles.csv"))
  write_csv(comparisons, file.path(summary_directory,
                                   "stage_c_comparisons.csv"))
  write_csv(long_form, file.path(summary_directory,
                                 "stage_c_matrices_long.csv"))
  print(diagnostics)
  cat("\nLag profile of the primary matrix:\n")
  print(lag_profiles[lag_profiles$name ==
                       paste0(weight_source, "_", correlation_method), ])
  log_message("Stage C complete.")
}

# ---------------------------------------------------------------- Stage D ----

if (stage == "D") {
  log_message("Stage D: independence calibration by parametric bootstrap.")
  inputs <- readRDS(file.path(cache_directory, "subset_inputs.rds"))
  residuals <- readRDS(file.path(cache_directory, "residuals.rds"))
  fits <- readRDS(file.path(cache_directory, "diagonal_fits.rds"))
  posterior_means <- inputs$beta_hat - residuals[[weight_source]]$residual
  seeds <- c(20260901L, 20260902L, 20260903L)[seq_len(n_replicates)]

  replicate_results <- list()
  for (seed in seeds) {
    log_message("  replicate seed ", seed, ": simulating independent errors.")
    simulated <- simulate_independent_replicate(
      posterior_means, inputs$adjusted_se, seed
    )
    data_list <- make_fash_data_list(simulated, inputs$time_grid)
    data_list <- Map(function(dataset, index) {
      dataset$SE <- as.numeric(inputs$adjusted_se[index, ])
      dataset
    }, data_list, seq_along(data_list))

    started <- Sys.time()
    log_message("    fitting FASH on the synthetic subset.")
    synthetic_raw <- fashr::fash(
      Y = "beta", smooth_var = "time", S = "SE", data_list = data_list,
      num_basis = inputs$settings$num_basis,
      order = inputs$settings$order,
      betaprec = inputs$settings$betaprec,
      pred_step = inputs$settings$pred_step,
      penalty = inputs$settings$penalty,
      grid = inputs$psd_grid,
      num_cores = num_cores, verbose = FALSE
    )
    names(synthetic_raw$lfdr) <- rownames(simulated)
    rownames(synthetic_raw$posterior_weights) <- rownames(simulated)
    log_message("    fitted in ",
                format(difftime(Sys.time(), started, units = "secs"),
                       digits = 4), ".")
    synthetic_bf <- fashr::BF_update(synthetic_raw)
    synthetic_fit <- if (weight_source == "bf") synthetic_bf else synthetic_raw

    synthetic_inputs <- inputs
    synthetic_inputs$beta_hat <- simulated
    synthetic_residuals <- compute_residuals(
      synthetic_inputs, synthetic_fit, num_cores
    )
    replicate_results[[as.character(seed)]] <- list(
      seed = seed,
      correlation = estimate_common_correlation(
        synthetic_residuals$z, correlation_method
      ),
      pi0_raw = synthetic_raw$prior_weights$prior_weight[
        synthetic_raw$prior_weights$psd == 0
      ],
      pi0_bf = synthetic_bf$prior_weights$prior_weight[
        synthetic_bf$prior_weights$psd == 0
      ],
      z_sd = stats::sd(synthetic_residuals$z)
    )
  }

  saveRDS(replicate_results,
          file.path(cache_directory, "calibration_replicates.rds"))

  calibration_summary <- do.call(rbind, lapply(replicate_results, function(x) {
    diagnostics <- correlation_diagnostics(
      x$correlation, paste0("independence_seed", x$seed)
    )
    diagnostics$pi0_raw <- if (length(x$pi0_raw) == 1L) x$pi0_raw else NA_real_
    diagnostics$pi0_bf <- if (length(x$pi0_bf) == 1L) x$pi0_bf else NA_real_
    diagnostics$residual_z_sd <- x$z_sd
    diagnostics
  }))
  calibration_lag <- do.call(rbind, lapply(replicate_results, function(x) {
    correlation_lag_profile(x$correlation,
                            paste0("independence_seed", x$seed))
  }))
  write_csv(calibration_summary,
            file.path(summary_directory, "stage_d_calibration.csv"))
  write_csv(calibration_lag,
            file.path(summary_directory, "stage_d_calibration_lag.csv"))
  print(calibration_summary)
  print(calibration_lag)
  log_message("Stage D complete.")
}

# ---------------------------------------------------------------- Stage E ----

if (stage == "E") {
  log_message("Stage E: dependent-likelihood refit.")
  inputs <- readRDS(file.path(cache_directory, "subset_inputs.rds"))
  fits <- readRDS(file.path(cache_directory, "diagonal_fits.rds"))
  if (matrix_choice == "permutation") {
    # The design-based donor-residual-block permutation null from 2026-08-10.
    # Refitting with it checks that the collapse below is a property of the
    # error correlation itself, not of the estimator used to recover it.
    permutation_path <- file.path(
      workflowr_root, "output", "revision_simulations", "internal",
      "multigene_null_beta_covariance", "summary",
      "primary_donor_residual_block_C_m_g_gt_0p900.csv"
    )
    permutation_frame <- utils::read.csv(permutation_path,
                                         stringsAsFactors = FALSE)
    correlation <- as.matrix(
      permutation_frame[, grep("^time_", names(permutation_frame))]
    )
    correlation <- (correlation + t(correlation)) / 2
    diag(correlation) <- 1
    matrix_label <- "permutation_null"
  } else if (matrix_choice == "residual") {
    matrices <- readRDS(file.path(cache_directory,
                                  "correlation_matrices.rds"))
    correlation <- matrices[[paste0(weight_source, "_", correlation_method)]]
    matrix_label <- "residual_correlation"
  } else {
    deconvolution <- readRDS(file.path(cache_directory, "deconvolution.rds"))
    solution_name <- if (matrix_choice == "deconvolved") {
      "unregularised"
    } else {
      "ridge_1e2"
    }
    correlation <- deconvolution$solutions[[solution_name]]$correlation
    matrix_label <- paste0("deconvolved_", solution_name)
  }
  if (is.null(correlation)) {
    stop("Could not resolve the requested correlation matrix.")
  }
  log_message("  using correlation matrix: ", matrix_label, ".")
  dimnames(correlation) <- list(colnames(inputs$beta_hat),
                                colnames(inputs$beta_hat))
  tag <- paste0("stage_e_", matrix_label)

  started <- Sys.time()
  dependent <- fit_fash_with_shared_correlation(
    beta_hat = inputs$beta_hat,
    adjusted_se = inputs$adjusted_se,
    time_grid = inputs$time_grid,
    correlation = correlation,
    settings = inputs$settings,
    psd_grid = inputs$psd_grid,
    num_cores = num_cores,
    verbose = FALSE
  )
  log_message("  dependent fit in ",
              format(difftime(Sys.time(), started, units = "secs"),
                     digits = 5), ".")
  dependent_bf <- tryCatch(
    fashr::BF_update(dependent$fit),
    error = function(error) {
      log_message("  BF_update failed: ", conditionMessage(error))
      NULL
    }
  )
  saveRDS(
    list(raw = dependent$fit, bf = dependent_bf,
         diagnostics = dependent$correlation_diagnostics,
         correlation = correlation),
    file.path(cache_directory, paste0("dependent_fits_", matrix_label, ".rds"))
  )

  pi0_of <- function(fit) {
    if (is.null(fit)) return(NA_real_)
    value <- fit$prior_weights$prior_weight[fit$prior_weights$psd == 0]
    if (length(value) == 1L) value else 0
  }
  rows <- list(
    list("diagonal", "raw", fits$raw),
    list("diagonal", "bf", fits$bf),
    list(matrix_label, "raw", dependent$fit),
    list(matrix_label, "bf", dependent_bf)
  )
  summary_table <- do.call(rbind, lapply(rows, function(row) {
    fit <- row[[3]]
    data.frame(
      likelihood = row[[1]], stage = row[[2]],
      pi0 = pi0_of(fit),
      mean_lfdr = if (is.null(fit)) NA_real_ else mean(fit$lfdr),
      calls = if (is.null(fit)) {
        NA_integer_
      } else {
        length(cumulative_lfdr_calls(fit$lfdr, alpha))
      },
      stringsAsFactors = FALSE
    )
  }))
  write_csv(summary_table,
            file.path(summary_directory, paste0(tag, "_summary.csv")))

  paired <- data.frame(
    pair_key = names(fits$bf$lfdr),
    diagonal_bf_lfdr = as.numeric(fits$bf$lfdr),
    dependent_bf_lfdr = if (is.null(dependent_bf)) {
      NA_real_
    } else {
      as.numeric(dependent_bf$lfdr)
    },
    diagonal_raw_lfdr = as.numeric(fits$raw$lfdr),
    dependent_raw_lfdr = as.numeric(dependent$fit$lfdr),
    stringsAsFactors = FALSE
  )
  write_csv(paired,
            file.path(summary_directory, paste0(tag, "_paired_lfdr.csv")))

  agreement <- data.frame(
    stage = c("raw", "bf"),
    spearman = c(
      safe_correlation(paired$diagonal_raw_lfdr, paired$dependent_raw_lfdr,
                       "spearman"),
      if (is.null(dependent_bf)) {
        NA_real_
      } else {
        safe_correlation(paired$diagonal_bf_lfdr, paired$dependent_bf_lfdr,
                         "spearman")
      }
    ),
    mean_absolute_difference = c(
      mean(abs(paired$diagonal_raw_lfdr - paired$dependent_raw_lfdr)),
      if (is.null(dependent_bf)) {
        NA_real_
      } else {
        mean(abs(paired$diagonal_bf_lfdr - paired$dependent_bf_lfdr))
      }
    ),
    stringsAsFactors = FALSE
  )
  write_csv(agreement,
            file.path(summary_directory, paste0(tag, "_agreement.csv")))
  print(summary_table)
  print(agreement)
  log_message("Stage E complete.")
}

# ---------------------------------------------------------------- Stage F ----

if (stage == "F") {
  log_message("Stage F: exact analytic independence calibration.")
  inputs <- readRDS(file.path(cache_directory, "subset_inputs.rds"))
  fits <- readRDS(file.path(cache_directory, "diagonal_fits.rds"))
  residuals <- readRDS(file.path(cache_directory, "residuals.rds"))
  matrices <- readRDS(file.path(cache_directory, "correlation_matrices.rds"))
  observed <- matrices[[paste0(weight_source, "_", correlation_method)]]

  shared_design <- build_shared_design(
    time_grid = inputs$time_grid,
    num_basis = inputs$settings$num_basis,
    order = inputs$settings$order,
    betaprec = inputs$settings$betaprec
  )
  fit <- fits[[weight_source]]
  scales <- psd_to_prior_scale(
    retained_psd_values(fit), inputs$settings$order, inputs$settings$pred_step
  )$prior_scale
  started <- Sys.time()
  expected <- expected_residual_z_moment(
    beta_hat = inputs$beta_hat,
    standard_errors = inputs$adjusted_se,
    posterior_weights = fit$posterior_weights,
    shared_design = shared_design,
    prior_scale = scales,
    num_cores = num_cores
  )
  log_message("  analytic calibration in ",
              format(difftime(Sys.time(), started, units = "secs"),
                     digits = 4), ".")

  saveRDS(expected, file.path(cache_directory, "analytic_calibration.rds"))

  observed_z_sd <- apply(residuals[[weight_source]]$z, 2, stats::sd)
  scale_table <- data.frame(
    time = inputs$time_grid,
    expected_z_sd_under_independence = expected$expected_z_sd,
    observed_z_sd = as.numeric(observed_z_sd),
    ratio = as.numeric(observed_z_sd) / expected$expected_z_sd,
    stringsAsFactors = FALSE
  )
  write_csv(scale_table, file.path(summary_directory, "stage_f_z_scale.csv"))

  diagnostics <- rbind(
    correlation_diagnostics(expected$correlation, "analytic_independence"),
    correlation_diagnostics(observed, "observed_residual")
  )
  lag_profiles <- rbind(
    correlation_lag_profile(expected$correlation, "analytic_independence"),
    correlation_lag_profile(observed, "observed_residual")
  )
  # The excess is what the observed matrix has that pure shrinkage cannot
  # explain. It is the only part that could be real error correlation.
  excess <- observed - expected$correlation
  write_csv(matrix_to_long(expected$correlation, "analytic_independence",
                           "analytic independence calibration"),
            file.path(summary_directory, "stage_f_analytic_long.csv"))
  write_csv(matrix_to_long(excess, "excess", "observed minus analytic"),
            file.path(summary_directory, "stage_f_excess_long.csv"))
  write_csv(diagnostics, file.path(summary_directory,
                                   "stage_f_diagnostics.csv"))
  write_csv(lag_profiles, file.path(summary_directory,
                                    "stage_f_lag_profiles.csv"))
  write_csv(
    compare_correlation_matrices(observed, expected$correlation,
                                 "observed_residual",
                                 "analytic_independence"),
    file.path(summary_directory, "stage_f_comparison.csv")
  )
  print(scale_table, digits = 4)
  print(diagnostics, digits = 4)
  wide <- reshape(lag_profiles, idvar = "lag", timevar = "name",
                  direction = "wide")
  names(wide) <- sub("mean_correlation.", "", names(wide), fixed = TRUE)
  wide$excess <- wide$observed_residual - wide$analytic_independence
  print(wide, digits = 3, row.names = FALSE)
  log_message("Stage F complete.")
}

# ---------------------------------------------------------------- Stage G ----

if (stage == "G") {
  log_message("Stage G: deconvolving the shrinkage operator.")
  inputs <- readRDS(file.path(cache_directory, "subset_inputs.rds"))
  fits <- readRDS(file.path(cache_directory, "diagonal_fits.rds"))
  residuals <- readRDS(file.path(cache_directory, "residuals.rds"))

  shared_design <- build_shared_design(
    time_grid = inputs$time_grid,
    num_basis = inputs$settings$num_basis,
    order = inputs$settings$order,
    betaprec = inputs$settings$betaprec
  )
  fit <- fits[[weight_source]]
  scales <- psd_to_prior_scale(
    retained_psd_values(fit), inputs$settings$order, inputs$settings$pred_step
  )$prior_scale

  started <- Sys.time()
  kernel <- residual_transfer_kronecker(
    standard_errors = inputs$adjusted_se,
    posterior_weights = fit$posterior_weights,
    shared_design = shared_design,
    prior_scale = scales
  )
  log_message("  accumulated the transfer kernel in ",
              format(difftime(Sys.time(), started, units = "secs"),
                     digits = 4), ".")

  z <- residuals[[weight_source]]$z
  observed_moment <- crossprod(z) / nrow(z)

  # Self-test: feeding the kernel its own prediction for P = I must return I.
  identity_target <- kernel %*% as.numeric(diag(1, length(inputs$time_grid)))
  identity_target <- matrix(identity_target, length(inputs$time_grid))
  identity_check <- deconvolve_error_correlation(kernel, identity_target)
  self_test <- data.frame(
    maximum_absolute_deviation_from_identity = max(abs(
      identity_check$solution - diag(1, length(inputs$time_grid))
    )),
    condition_number = identity_check$condition_number,
    effective_rank = identity_check$effective_rank,
    free_parameters = length(inputs$time_grid) *
      (length(inputs$time_grid) + 1) / 2,
    stringsAsFactors = FALSE
  )
  write_csv(self_test, file.path(summary_directory, "stage_g_self_test.csv"))
  print(self_test, digits = 4)

  solutions <- list(
    unregularised = deconvolve_error_correlation(kernel, observed_moment),
    ridge_1e3 = deconvolve_error_correlation(kernel, observed_moment,
                                             ridge = 1e-3),
    ridge_1e2 = deconvolve_error_correlation(kernel, observed_moment,
                                             ridge = 1e-2)
  )
  saveRDS(list(kernel = kernel, observed_moment = observed_moment,
               solutions = solutions, self_test = self_test),
          file.path(cache_directory, "deconvolution.rds"))

  solution_summary <- do.call(rbind, lapply(names(solutions), function(name) {
    result <- solutions[[name]]
    frame <- correlation_diagnostics(result$correlation,
                                     paste0("deconvolved_", name))
    frame$solution_minimum_eigenvalue <- result$solution_minimum_eigenvalue
    frame$minimum_implied_diagonal <- min(result$implied_diagonal)
    frame$maximum_implied_diagonal <- max(result$implied_diagonal)
    frame$lag1 <- mean(result$correlation[cbind(
      seq_len(nrow(result$correlation) - 1L),
      seq_len(nrow(result$correlation) - 1L) + 1L
    )])
    frame
  }))
  write_csv(solution_summary,
            file.path(summary_directory, "stage_g_solutions.csv"))
  write_csv(
    do.call(rbind, lapply(names(solutions), function(name) {
      correlation_lag_profile(solutions[[name]]$correlation,
                              paste0("deconvolved_", name))
    })),
    file.path(summary_directory, "stage_g_lag_profiles.csv")
  )
  spectrum <- data.frame(
    index = seq_along(solutions$unregularised$singular_values),
    singular_value = solutions$unregularised$singular_values,
    stringsAsFactors = FALSE
  )
  spectrum$relative <- spectrum$singular_value /
    max(spectrum$singular_value)
  write_csv(spectrum, file.path(summary_directory, "stage_g_spectrum.csv"))
  write_csv(
    do.call(rbind, lapply(
      seq_along(solutions$unregularised$least_identified_directions),
      function(index) {
        entry <- solutions$unregularised$least_identified_directions[[index]]
        frame <- matrix_to_long(
          entry$direction, paste0("weak", index),
          paste0("singular value ", signif(entry$singular_value, 4))
        )
        frame$singular_value <- entry$singular_value
        frame
      }
    )),
    file.path(summary_directory, "stage_g_weak_directions.csv")
  )
  print(solution_summary, digits = 4)
  cat("\nSingular-value spectrum of the symmetric transfer operator:\n")
  print(spectrum[c(1:5, 130:136), ], digits = 4, row.names = FALSE)
  log_message("Stage G complete.")
}
