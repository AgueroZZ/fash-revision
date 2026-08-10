#!/usr/bin/env Rscript

# Run the paired R4 sensitivity analysis for time-correlated expression errors.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the shared revision simulation functions.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1] == length(args)) {
    return(default)
  }
  args[hit[1] + 1L]
}

as_flag <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

parse_seed_list <- function(x) {
  pieces <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  if (length(pieces) == 0L || any(!nzchar(pieces))) {
    stop("--seed-list must contain one or more comma-separated integer seeds.")
  }
  seeds <- suppressWarnings(as.integer(pieces))
  if (any(is.na(seeds)) || anyDuplicated(seeds)) {
    stop("--seed-list must contain unique integer seeds.")
  }
  seeds
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

lag_average_correlations <- function(x) {
  x <- as.matrix(x)
  if (nrow(x) < 3L || ncol(x) < 2L || any(!is.finite(x))) {
    stop("The lag-correlation input must be a finite matrix with multiple rows and columns.")
  }
  correlation <- stats::cor(x)
  n_time <- ncol(correlation)
  vapply(seq_len(n_time - 1L), function(lag) {
    mean(correlation[cbind(
      seq_len(n_time - lag),
      (lag + 1L):n_time
    )])
  }, numeric(1))
}

make_lag_rows <- function(x, seed, condition, diagnostic) {
  lag_correlation <- lag_average_correlations(x)
  data.frame(
    seed = seed,
    condition = condition,
    diagnostic = diagnostic,
    lag = seq_along(lag_correlation),
    mean_correlation = lag_correlation,
    n_units = nrow(x),
    stringsAsFactors = FALSE
  )
}

compute_dependence_diagnostics <- function(out, seed, condition) {
  beta_hat <- as.matrix(out$eqtl_summary$beta_hat)
  se <- as.matrix(out$eqtl_summary$se)
  true_beta <- as.matrix(out$true_beta)
  if (!identical(dim(beta_hat), dim(se)) ||
      !identical(dim(beta_hat), dim(true_beta)) ||
      any(!is.finite(beta_hat)) || any(!is.finite(se)) ||
      any(!is.finite(true_beta)) || any(se <= 0)) {
    stop("Invalid beta-hat, SE, or truth matrix in an R4 replicate.")
  }

  truth_known_error <- (beta_hat - true_beta) / se
  dynamic_null <- out$unit_info$effect_class %in% c("constant", "zero")
  if (sum(dynamic_null) < 3L) {
    stop("The R4 replicate does not contain enough dynamic-null variants.")
  }
  null_beta_hat <- beta_hat[dynamic_null, , drop = FALSE]
  null_se <- se[dynamic_null, , drop = FALSE]
  weights <- 1 / null_se^2
  weighted_mean <- rowSums(null_beta_hat * weights) / rowSums(weights)
  centered_null_residual <- sweep(
    null_beta_hat,
    1,
    weighted_mean,
    `-`
  ) / null_se

  rbind(
    make_lag_rows(
      truth_known_error,
      seed = seed,
      condition = condition,
      diagnostic = "Truth-known standardized error"
    ),
    make_lag_rows(
      centered_null_residual,
      seed = seed,
      condition = condition,
      diagnostic = "Centered dynamic-null residual"
    )
  )
}

extract_expression_errors <- function(out) {
  G <- out$genotype
  covariates <- out$covariates
  beta <- out$true_beta
  expression <- out$expression
  intercepts <- out$expression_simulation$intercepts
  covariate_effects <- out$expression_simulation$covariate_effects
  n_donors <- nrow(G)
  n_variants <- ncol(G)
  n_time <- ncol(beta)
  errors <- array(
    NA_real_,
    dim = c(n_donors, n_variants, n_time),
    dimnames = dimnames(expression)
  )

  for (tt in seq_len(n_time)) {
    genetic_mean <- sweep(G, 2, beta[, tt], `*`)
    covariate_mean <- if (is.null(covariates)) {
      matrix(0, nrow = n_donors, ncol = n_variants)
    } else {
      covariates %*% covariate_effects[, , tt]
    }
    systematic_mean <- sweep(
      genetic_mean + covariate_mean,
      2,
      intercepts[, tt],
      `+`
    )
    errors[, , tt] <- expression[, , tt] - systematic_mean
  }
  errors
}

validate_paired_simulations <- function(independent, correlated, correlation) {
  paired_fields <- list(
    genotype = list(independent$genotype, correlated$genotype),
    covariates = list(independent$covariates, correlated$covariates),
    true_beta = list(independent$true_beta, correlated$true_beta),
    unit_info = list(independent$unit_info, correlated$unit_info),
    intercepts = list(
      independent$expression_simulation$intercepts,
      correlated$expression_simulation$intercepts
    ),
    covariate_effects = list(
      independent$expression_simulation$covariate_effects,
      correlated$expression_simulation$covariate_effects
    )
  )
  mismatched <- names(paired_fields)[!vapply(paired_fields, function(x) {
    identical(x[[1]], x[[2]])
  }, logical(1))]
  if (length(mismatched) > 0L) {
    stop(
      "The paired R4 simulations differ before the error transformation: ",
      paste(mismatched, collapse = ", ")
    )
  }

  independent_errors <- extract_expression_errors(independent)
  correlated_errors <- extract_expression_errors(correlated)
  n_series <- dim(independent_errors)[1] * dim(independent_errors)[2]
  n_time <- dim(independent_errors)[3]
  expected_correlated <- matrix(
    independent_errors,
    nrow = n_series,
    ncol = n_time
  ) %*% chol(correlation)
  observed_correlated <- matrix(
    correlated_errors,
    nrow = n_series,
    ncol = n_time
  )
  maximum_difference <- max(abs(expected_correlated - observed_correlated))
  if (!is.finite(maximum_difference) || maximum_difference > 1e-10) {
    stop(
      "The paired expression errors do not share the intended Gaussian innovations; maximum difference = ",
      format(maximum_difference, scientific = TRUE)
    )
  }
  maximum_difference
}

extract_fash_metrics <- function(out, seed, condition, methods) {
  available_methods <- unique(out$alpha_curve$method)
  missing_methods <- setdiff(methods, available_methods)
  if (length(missing_methods) > 0L) {
    stop("Missing R4 FASH methods: ", paste(missing_methods, collapse = ", "))
  }
  alpha_curve <- out$alpha_curve[out$alpha_curve$method %in% methods, ]
  alpha_curve$seed <- seed
  alpha_curve$condition <- condition
  alpha_005 <- alpha_curve[abs(alpha_curve$alpha - 0.05) < 1e-12, ]
  if (nrow(alpha_005) != length(methods)) {
    stop("The R4 alpha-0.05 result is incomplete.")
  }

  pi0 <- data.frame(
    seed = seed,
    condition = condition,
    method = c("FASH-IWP1", "FASH-IWP1", "FASH-linear", "FASH-linear"),
    fit = c("Raw", "BF-corrected", "Raw", "BF-corrected"),
    estimated_pi0 = c(
      constant_component_prior_weight(out$fash_fits$fash_iwp1_raw),
      constant_component_prior_weight(out$fash_fits$fash_iwp1_bf),
      constant_component_prior_weight(out$simplified_fit),
      constant_component_prior_weight(out$simplified_fit_bf)
    ),
    stringsAsFactors = FALSE
  )
  list(alpha_curve = alpha_curve, alpha_005 = alpha_005, pi0 = pi0)
}

summarize_condition_alpha <- function(alpha_rows) {
  conditions <- unique(alpha_rows$condition)
  out <- lapply(conditions, function(condition) {
    summary <- summarize_mc_alpha_curves(
      alpha_rows[alpha_rows$condition == condition, , drop = FALSE]
    )
    summary$condition <- condition
    summary
  })
  out <- do.call(rbind, out)
  out <- out[, c("condition", setdiff(names(out), "condition")), drop = FALSE]
  rownames(out) <- NULL
  out
}

summarize_condition_pi0 <- function(pi0_rows) {
  conditions <- unique(pi0_rows$condition)
  out <- lapply(conditions, function(condition) {
    summary <- summarize_mc_pi0(
      pi0_rows[pi0_rows$condition == condition, , drop = FALSE]
    )
    summary$condition <- condition
    summary
  })
  out <- do.call(rbind, out)
  out <- out[, c("condition", setdiff(names(out), "condition")), drop = FALSE]
  rownames(out) <- NULL
  out
}

summarize_paired_alpha <- function(alpha_rows) {
  key <- c("seed", "scenario", "target", "method", "alpha")
  metrics <- c(
    "n_discoveries", "false_discoveries", "true_positives",
    "empirical_fdr", "power"
  )
  independent <- alpha_rows[
    alpha_rows$condition == "Independent",
    c(key, metrics),
    drop = FALSE
  ]
  correlated <- alpha_rows[
    alpha_rows$condition == "Correlated",
    c(key, metrics),
    drop = FALSE
  ]
  paired <- merge(
    independent,
    correlated,
    by = key,
    suffixes = c("_independent", "_correlated"),
    all = FALSE,
    sort = FALSE
  )
  expected_rows <- nrow(independent)
  if (nrow(paired) != expected_rows || nrow(correlated) != expected_rows) {
    stop("Independent and correlated alpha curves are not completely paired.")
  }
  for (metric in metrics) {
    paired[[paste0("delta_", metric)]] <-
      paired[[paste0(metric, "_correlated")]] -
      paired[[paste0(metric, "_independent")]]
  }

  groups <- split(
    paired,
    list(paired$scenario, paired$target, paired$method, paired$alpha),
    drop = TRUE
  )
  summaries <- lapply(groups, function(x) {
    values <- lapply(metrics, function(metric) {
      summarize_mc_values(x[[paste0("delta_", metric)]])
    })
    names(values) <- metrics
    data.frame(
      scenario = x$scenario[1],
      target = x$target[1],
      method = x$method[1],
      alpha = x$alpha[1],
      n_replications = length(unique(x$seed)),
      mean_delta_discoveries = values$n_discoveries[["mean"]],
      discoveries_delta_ci_lower = values$n_discoveries[["lower"]],
      discoveries_delta_ci_upper = values$n_discoveries[["upper"]],
      mean_delta_false_discoveries = values$false_discoveries[["mean"]],
      false_discoveries_delta_ci_lower = values$false_discoveries[["lower"]],
      false_discoveries_delta_ci_upper = values$false_discoveries[["upper"]],
      mean_delta_true_positives = values$true_positives[["mean"]],
      true_positives_delta_ci_lower = values$true_positives[["lower"]],
      true_positives_delta_ci_upper = values$true_positives[["upper"]],
      mean_delta_fdr = values$empirical_fdr[["mean"]],
      fdr_delta_ci_lower = values$empirical_fdr[["lower"]],
      fdr_delta_ci_upper = values$empirical_fdr[["upper"]],
      mean_delta_power = values$power[["mean"]],
      power_delta_ci_lower = values$power[["lower"]],
      power_delta_ci_upper = values$power[["upper"]],
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, summaries)
  out$method_rank <- rank_revision_methods(out$method)
  out <- out[order(out$method_rank, out$method, out$alpha), , drop = FALSE]
  out$method_rank <- NULL
  rownames(out) <- NULL
  out
}

summarize_paired_pi0 <- function(pi0_rows) {
  key <- c("seed", "method", "fit")
  independent <- pi0_rows[
    pi0_rows$condition == "Independent",
    c(key, "estimated_pi0"),
    drop = FALSE
  ]
  correlated <- pi0_rows[
    pi0_rows$condition == "Correlated",
    c(key, "estimated_pi0"),
    drop = FALSE
  ]
  paired <- merge(
    independent,
    correlated,
    by = key,
    suffixes = c("_independent", "_correlated"),
    all = FALSE,
    sort = FALSE
  )
  if (nrow(paired) != nrow(independent) || nrow(correlated) != nrow(independent)) {
    stop("Independent and correlated pi0 estimates are not completely paired.")
  }
  paired$delta_pi0 <- paired$estimated_pi0_correlated -
    paired$estimated_pi0_independent
  groups <- split(paired, list(paired$method, paired$fit), drop = TRUE)
  out <- do.call(rbind, lapply(groups, function(x) {
    values <- summarize_mc_values(x$delta_pi0)
    data.frame(
      method = x$method[1],
      fit = x$fit[1],
      n_replications = length(unique(x$seed)),
      mean_delta_pi0 = values[["mean"]],
      pi0_delta_sd = values[["sd"]],
      pi0_delta_mc_se = values[["se"]],
      pi0_delta_ci_lower = values[["lower"]],
      pi0_delta_ci_upper = values[["upper"]],
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

summarize_lag_diagnostics <- function(lag_rows) {
  groups <- split(
    lag_rows,
    list(lag_rows$condition, lag_rows$diagnostic, lag_rows$lag),
    drop = TRUE
  )
  out <- do.call(rbind, lapply(groups, function(x) {
    values <- summarize_mc_values(x$mean_correlation)
    data.frame(
      condition = x$condition[1],
      diagnostic = x$diagnostic[1],
      lag = x$lag[1],
      n_replications = length(unique(x$seed)),
      n_units = unique(x$n_units)[1],
      mean_correlation = values[["mean"]],
      correlation_sd = values[["sd"]],
      correlation_mc_se = values[["se"]],
      ci_lower = values[["lower"]],
      ci_upper = values[["upper"]],
      stringsAsFactors = FALSE
    )
  }))
  out <- out[order(out$diagnostic, out$condition, out$lag), , drop = FALSE]
  rownames(out) <- NULL
  out
}

validate_against_r1_cache <- function(alpha_rows,
                                      pi0_rows,
                                      workflowr_root,
                                      methods,
                                      expected_seeds) {
  reference_dir <- file.path(
    workflowr_root,
    "output",
    "revision_simulations",
    "mc",
    "r1_random_bspline_main_effect_pilot5",
    "summary"
  )
  alpha_path <- file.path(reference_dir, "all_replicate_alpha_curves.csv")
  pi0_path <- file.path(reference_dir, "all_replicate_pi0.csv")
  if (!file.exists(alpha_path) || !file.exists(pi0_path)) {
    stop("The default R1 reference cache is missing.")
  }

  reference_alpha <- utils::read.csv(alpha_path, stringsAsFactors = FALSE)
  reference_alpha <- reference_alpha[
    reference_alpha$method %in% methods &
      reference_alpha$seed %in% expected_seeds,
    ,
    drop = FALSE
  ]
  observed_alpha <- alpha_rows[
    alpha_rows$condition == "Independent",
    ,
    drop = FALSE
  ]
  alpha_key <- c("seed", "method", "alpha")
  alpha_metrics <- c(
    "n_discoveries", "false_discoveries", "true_positives",
    "empirical_fdr", "power"
  )
  reference_alpha <- reference_alpha[
    do.call(order, unname(reference_alpha[alpha_key])),
    c(alpha_key, alpha_metrics),
    drop = FALSE
  ]
  observed_alpha <- observed_alpha[
    do.call(order, unname(observed_alpha[alpha_key])),
    c(alpha_key, alpha_metrics),
    drop = FALSE
  ]
  rownames(reference_alpha) <- NULL
  rownames(observed_alpha) <- NULL
  if (!isTRUE(all.equal(
    reference_alpha[alpha_key],
    observed_alpha[alpha_key],
    check.attributes = FALSE
  ))) {
    stop("The paired independent alpha-curve keys do not match the R1 cache.")
  }
  alpha_difference <- max(abs(
    as.matrix(reference_alpha[alpha_metrics]) -
      as.matrix(observed_alpha[alpha_metrics])
  ))
  if (!is.finite(alpha_difference) || alpha_difference > 1e-12) {
    stop("The paired independent alpha curves do not reproduce the R1 cache.")
  }

  reference_pi0 <- utils::read.csv(pi0_path, stringsAsFactors = FALSE)
  reference_pi0 <- reference_pi0[reference_pi0$seed %in% expected_seeds, ]
  observed_pi0 <- pi0_rows[pi0_rows$condition == "Independent", ]
  pi0_key <- c("seed", "method", "fit")
  reference_pi0 <- reference_pi0[
    do.call(order, unname(reference_pi0[pi0_key])),
    c(pi0_key, "estimated_pi0"),
    drop = FALSE
  ]
  observed_pi0 <- observed_pi0[
    do.call(order, unname(observed_pi0[pi0_key])),
    c(pi0_key, "estimated_pi0"),
    drop = FALSE
  ]
  rownames(reference_pi0) <- NULL
  rownames(observed_pi0) <- NULL
  if (!isTRUE(all.equal(
    reference_pi0[pi0_key],
    observed_pi0[pi0_key],
    check.attributes = FALSE
  ))) {
    stop("The paired independent pi0 keys do not match the R1 cache.")
  }
  pi0_difference <- max(abs(
    reference_pi0$estimated_pi0 - observed_pi0$estimated_pi0
  ))
  if (!is.finite(pi0_difference) || pi0_difference > 1e-12) {
    stop("The paired independent pi0 estimates do not reproduce the R1 cache.")
  }

  data.frame(
    component = c("alpha curves", "pi0 estimates"),
    maximum_absolute_difference = c(alpha_difference, pi0_difference),
    tolerance = 1e-12,
    passed = TRUE,
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "shared",
  "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "r4_correlated_errors",
  "plotting.R"
))

J <- as.integer(get_arg("--J", "1000"))
n_donors <- as.integer(get_arg("--n-donors", "19"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
expression_noise_sd <- as.numeric(get_arg("--noise-sd", "1"))
rho <- as.numeric(get_arg("--rho", "-0.09"))
correlation_structure <- get_arg("--correlation-structure", "lag1_only")
num_cores <- as.integer(get_arg("--num-cores", "4"))
num_basis <- as.integer(get_arg("--num-basis", "20"))
seed_list <- parse_seed_list(get_arg(
  "--seed-list",
  "12345,22345,32345,42345,52345"
))
output_id <- get_arg(
  "--output-id",
  "r4_correlated_errors_lag1_rho_m0p09_pilot5"
)
overwrite <- as_flag(get_arg("--overwrite", "false"))

if (is.na(J) || J < 10L || is.na(n_donors) ||
    n_donors < n_covariates + 3L || is.na(n_covariates) ||
    n_covariates < 0L || !is.finite(expression_noise_sd) ||
    expression_noise_sd <= 0 || !is.finite(rho) ||
    is.na(num_cores) || num_cores < 1L || is.na(num_basis) ||
    num_basis < 2L || !nzchar(output_id)) {
  stop("Invalid R4 simulation arguments.")
}
if (!correlation_structure %in% c("lag1_only", "ar1")) {
  stop("--correlation-structure must be either lag1_only or ar1.")
}

time_grid <- make_time_grid()
error_correlation <- if (correlation_structure == "lag1_only") {
  make_lag1_correlation(length(time_grid), rho)
} else {
  if (abs(rho) >= 1) {
    stop("AR(1) rho must satisfy abs(rho) < 1.")
  }
  make_ar1_correlation(length(time_grid), rho)
}
error_correlation <- validate_time_correlation(
  error_correlation,
  n_time = length(time_grid)
)

class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
scenario <- "genotype_random_bspline_main_effect_dynamic_eqtl"
fash_methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "FASH-linear-Raw",
  "FASH-linear-BF"
)
conditions <- c("Independent", "Correlated")
empirical_centered_lag1 <- -0.1520468
centering_only_lag1 <- -0.062427

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "mc",
  output_id
)
replicate_dir <- file.path(output_dir, "replicates")
summary_dir <- file.path(output_dir, "summary")
figure_dir <- file.path(output_dir, "figures")
full_fit_dir <- file.path(output_dir, "full_fits")
invisible(lapply(
  c(output_dir, replicate_dir, summary_dir, figure_dir, full_fit_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

configuration <- list(
  output_id = output_id,
  scenario = scenario,
  J = J,
  n_donors = n_donors,
  time_grid = time_grid,
  n_covariates = n_covariates,
  expression_noise_sd = expression_noise_sd,
  dynamic_main_effect_sd = 1,
  num_basis = num_basis,
  class_probs = class_probs,
  dynamic_amplitude = 2,
  bspline_df = 6,
  bspline_coefficient_sd = 1,
  true_pi0 = 0.8,
  seed_list = seed_list,
  full_fit_seed = seed_list[1],
  conditions = conditions,
  correlation_structure = correlation_structure,
  rho = rho,
  expression_error_correlation = error_correlation,
  empirical_centered_lag1_target = empirical_centered_lag1,
  centering_only_lag1_benchmark = centering_only_lag1,
  r1_reference_output_id = "r1_random_bspline_main_effect_pilot5"
)
configuration_path <- file.path(output_dir, "configuration.rds")
if (file.exists(configuration_path) && !overwrite) {
  cached_configuration <- readRDS(configuration_path)
  if (!isTRUE(all.equal(cached_configuration, configuration))) {
    stop(
      "The existing output id has different settings. Choose a new --output-id or use --overwrite true."
    )
  }
} else {
  saveRDS(configuration, configuration_path)
}

run_condition <- function(seed, condition) {
  condition_correlation <- if (condition == "Independent") {
    NULL
  } else {
    error_correlation
  }
  message("Running R4 seed ", seed, ", condition: ", condition, ".")
  out <- run_genotype_level_bspline_eqtl_simulation(
    n_donors = n_donors,
    n_variants = J,
    time_grid = time_grid,
    n_covariates = n_covariates,
    class_probs = class_probs,
    expression_noise_sd = expression_noise_sd,
    dynamic_main_effect_sd = 1,
    scenario = scenario,
    alpha = 0.05,
    seed = seed,
    num_cores = num_cores,
    num_basis = num_basis,
    save_outputs = FALSE,
    verbose = FALSE,
    expression_error_correlation = condition_correlation
  )
  true_pi0 <- dynamic_null_proportion(out$unit_info, target = "dynamic")
  if (!isTRUE(all.equal(true_pi0, configuration$true_pi0))) {
    stop("The R4 dynamic-null proportion does not match the configuration.")
  }
  metrics <- extract_fash_metrics(
    out,
    seed = seed,
    condition = condition,
    methods = fash_methods
  )
  diagnostics <- compute_dependence_diagnostics(
    out,
    seed = seed,
    condition = condition
  )
  list(out = out, metrics = metrics, diagnostics = diagnostics)
}

make_replicate <- function(seed) {
  condition_results <- lapply(conditions, function(condition) {
    run_condition(seed, condition)
  })
  names(condition_results) <- conditions
  maximum_pairing_difference <- validate_paired_simulations(
    condition_results$Independent$out,
    condition_results$Correlated$out,
    correlation = error_correlation
  )

  if (identical(seed, seed_list[1])) {
    saveRDS(
      condition_results$Independent$out,
      file.path(full_fit_dir, paste0("seed_", seed, "_independent.rds"))
    )
    saveRDS(
      condition_results$Correlated$out,
      file.path(full_fit_dir, paste0("seed_", seed, "_correlated.rds"))
    )
  }

  list(
    configuration = configuration,
    seed = seed,
    alpha_curve = do.call(rbind, lapply(
      condition_results,
      function(x) x$metrics$alpha_curve
    )),
    alpha_005 = do.call(rbind, lapply(
      condition_results,
      function(x) x$metrics$alpha_005
    )),
    pi0 = do.call(rbind, lapply(
      condition_results,
      function(x) x$metrics$pi0
    )),
    lag_diagnostics = do.call(rbind, lapply(
      condition_results,
      `[[`,
      "diagnostics"
    )),
    maximum_pairing_difference = maximum_pairing_difference
  )
}

validate_replicate <- function(replicate, seed) {
  required_fields <- c(
    "configuration", "seed", "alpha_curve", "alpha_005", "pi0",
    "lag_diagnostics", "maximum_pairing_difference"
  )
  if (!all(required_fields %in% names(replicate)) ||
      !identical(replicate$seed, seed) ||
      !isTRUE(all.equal(replicate$configuration, configuration)) ||
      !is.finite(replicate$maximum_pairing_difference) ||
      replicate$maximum_pairing_difference > 1e-10) {
    return(FALSE)
  }
  all(fash_methods %in% unique(replicate$alpha_curve$method)) &&
    identical(sort(unique(replicate$alpha_curve$condition)), sort(conditions)) &&
    nrow(replicate$alpha_005) == length(fash_methods) * length(conditions) &&
    nrow(replicate$pi0) == 4L * length(conditions) &&
    nrow(replicate$lag_diagnostics) ==
      2L * (length(time_grid) - 1L) * length(conditions)
}

replicates <- lapply(seed_list, function(seed) {
  replicate_path <- file.path(replicate_dir, paste0("seed_", seed, ".rds"))
  if (file.exists(replicate_path) && !overwrite) {
    cached <- readRDS(replicate_path)
    if (validate_replicate(cached, seed)) {
      message("Reusing paired R4 replicate cache: ", replicate_path)
      return(cached)
    }
    stop("Cached R4 replicate does not match the requested settings: ", replicate_path)
  }
  replicate <- make_replicate(seed)
  saveRDS(replicate, replicate_path)
  replicate
})

all_alpha <- do.call(rbind, lapply(replicates, `[[`, "alpha_curve"))
all_alpha_005 <- do.call(rbind, lapply(replicates, `[[`, "alpha_005"))
all_pi0 <- do.call(rbind, lapply(replicates, `[[`, "pi0"))
all_lag <- do.call(rbind, lapply(replicates, `[[`, "lag_diagnostics"))
condition_alpha <- summarize_condition_alpha(all_alpha)
condition_alpha_005 <- condition_alpha[
  abs(condition_alpha$alpha - 0.05) < 1e-12,
  ,
  drop = FALSE
]
condition_pi0 <- summarize_condition_pi0(all_pi0)
paired_alpha <- summarize_paired_alpha(all_alpha)
paired_alpha_005 <- paired_alpha[
  abs(paired_alpha$alpha - 0.05) < 1e-12,
  ,
  drop = FALSE
]
paired_pi0 <- summarize_paired_pi0(all_pi0)
lag_summary <- summarize_lag_diagnostics(all_lag)

default_reference_setting <-
  identical(J, 1000L) &&
  identical(n_donors, 19L) &&
  identical(n_covariates, 5L) &&
  isTRUE(all.equal(expression_noise_sd, 1)) &&
  identical(seed_list, c(12345L, 22345L, 32345L, 42345L, 52345L)) &&
  identical(num_basis, 20L)
reference_check <- if (default_reference_setting) {
  validate_against_r1_cache(
    all_alpha,
    all_pi0,
    workflowr_root = workflowr_root,
    methods = fash_methods,
    expected_seeds = seed_list
  )
} else {
  data.frame(
    component = "R1 cache cross-check",
    maximum_absolute_difference = NA_real_,
    tolerance = NA_real_,
    passed = NA,
    stringsAsFactors = FALSE
  )
}

pairing_check <- data.frame(
  seed = seed_list,
  maximum_absolute_error_difference = vapply(
    replicates,
    `[[`,
    numeric(1),
    "maximum_pairing_difference"
  ),
  tolerance = 1e-10,
  passed = TRUE,
  stringsAsFactors = FALSE
)
correlation_target <- data.frame(
  source_subset = "Top 500 highest-lfdr variants, one variant per gene",
  observed_centered_lag1_correlation = empirical_centered_lag1,
  centering_only_lag1_correlation = centering_only_lag1,
  inferred_excess_lag1_correlation =
    empirical_centered_lag1 - centering_only_lag1,
  simulated_expression_error_lag1_correlation = rho,
  correlation_structure = correlation_structure,
  stringsAsFactors = FALSE
)

write_csv(all_alpha, file.path(summary_dir, "all_condition_alpha_curves.csv"))
write_csv(all_alpha_005, file.path(summary_dir, "all_condition_alpha005.csv"))
write_csv(condition_alpha, file.path(summary_dir, "condition_mc_alpha_curve.csv"))
write_csv(condition_alpha_005, file.path(summary_dir, "condition_mc_alpha005_summary.csv"))
write_csv(paired_alpha, file.path(summary_dir, "paired_alpha_differences.csv"))
write_csv(paired_alpha_005, file.path(summary_dir, "paired_alpha005_differences.csv"))
write_csv(all_pi0, file.path(summary_dir, "all_condition_pi0.csv"))
write_csv(condition_pi0, file.path(summary_dir, "condition_mc_pi0_summary.csv"))
write_csv(paired_pi0, file.path(summary_dir, "paired_pi0_differences.csv"))
write_csv(all_lag, file.path(summary_dir, "all_lag_correlations.csv"))
write_csv(lag_summary, file.path(summary_dir, "lag_correlation_summary.csv"))
write_csv(pairing_check, file.path(summary_dir, "paired_innovation_check.csv"))
write_csv(reference_check, file.path(summary_dir, "r1_reference_check.csv"))
write_csv(correlation_target, file.path(summary_dir, "correlation_target.csv"))

plot_r4_condition_curves(
  condition_alpha,
  metric = "power",
  file = file.path(figure_dir, "condition_mc_power.png")
)
plot_r4_condition_curves(
  condition_alpha,
  metric = "fdr",
  file = file.path(figure_dir, "condition_mc_fdr.png")
)
plot_r4_lag_diagnostics(
  lag_summary,
  file = file.path(figure_dir, "lag_correlation_diagnostics.png"),
  empirical_centered_target = empirical_centered_lag1
)

cat("\nR4 condition summaries at alpha = 0.05:\n")
print(condition_alpha_005)
cat("\nPaired correlated-minus-independent differences at alpha = 0.05:\n")
print(paired_alpha_005)
cat("\nR4 lag-1 dependence diagnostics:\n")
print(lag_summary[lag_summary$lag == 1L, ])
cat("\nR1 reference-cache check:\n")
print(reference_check)
