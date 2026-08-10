#!/usr/bin/env Rscript

# Run paired R1 simulations under two complete real-data correlation matrices.

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
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(args, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(args[equals_hit[1]], nchar(equals_prefix) + 1L))
  }
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
  seeds <- suppressWarnings(as.integer(pieces))
  if (length(seeds) == 0L || anyNA(seeds) || anyDuplicated(seeds)) {
    stop("--seed-list must contain unique comma-separated integer seeds.")
  }
  seeds
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

extract_expression_errors <- function(out) {
  genotype <- out$genotype
  covariates <- out$covariates
  true_beta <- out$true_beta
  expression <- out$expression
  intercepts <- out$expression_simulation$intercepts
  covariate_effects <- out$expression_simulation$covariate_effects
  n_donors <- nrow(genotype)
  n_variants <- ncol(genotype)
  n_time <- ncol(true_beta)
  errors <- array(
    NA_real_,
    dim = c(n_donors, n_variants, n_time),
    dimnames = dimnames(expression)
  )

  for (time_index in seq_len(n_time)) {
    genetic_mean <- sweep(
      genotype,
      2L,
      true_beta[, time_index],
      `*`
    )
    covariate_mean <- if (is.null(covariates)) {
      matrix(0, nrow = n_donors, ncol = n_variants)
    } else {
      covariates %*% covariate_effects[, , time_index]
    }
    systematic_mean <- sweep(
      genetic_mean + covariate_mean,
      2L,
      intercepts[, time_index],
      `+`
    )
    errors[, , time_index] <- expression[, , time_index] - systematic_mean
  }
  errors
}

validate_against_independent <- function(independent,
                                         candidate,
                                         correlation,
                                         independent_errors) {
  paired_fields <- list(
    genotype = list(independent$genotype, candidate$genotype),
    covariates = list(independent$covariates, candidate$covariates),
    true_beta = list(independent$true_beta, candidate$true_beta),
    unit_info = list(independent$unit_info, candidate$unit_info),
    intercepts = list(
      independent$expression_simulation$intercepts,
      candidate$expression_simulation$intercepts
    ),
    covariate_effects = list(
      independent$expression_simulation$covariate_effects,
      candidate$expression_simulation$covariate_effects
    )
  )
  mismatched <- names(paired_fields)[!vapply(paired_fields, function(x) {
    identical(x[[1]], x[[2]])
  }, logical(1))]
  if (length(mismatched) > 0L) {
    stop(
      "The paired full-matrix simulations differ before the error transform: ",
      paste(mismatched, collapse = ", ")
    )
  }

  n_series <- dim(independent_errors)[1] * dim(independent_errors)[2]
  n_time <- dim(independent_errors)[3]
  expected_errors <- matrix(
    independent_errors,
    nrow = n_series,
    ncol = n_time
  ) %*% chol(correlation)
  observed_errors <- matrix(
    extract_expression_errors(candidate),
    nrow = n_series,
    ncol = n_time
  )
  maximum_difference <- max(abs(expected_errors - observed_errors))
  if (!is.finite(maximum_difference) || maximum_difference > 1e-10) {
    stop(
      "The paired full-matrix errors do not share the intended innovations; maximum difference = ",
      format(maximum_difference, scientific = TRUE)
    )
  }
  maximum_difference
}

extract_iwp_metrics <- function(out, seed, condition) {
  methods <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
  alpha_005 <- out$alpha_curve[
    out$alpha_curve$method %in% methods &
      abs(out$alpha_curve$alpha - 0.05) < 1e-12,
    ,
    drop = FALSE
  ]
  if (nrow(alpha_005) != length(methods) ||
      !identical(sort(alpha_005$method), sort(methods))) {
    stop("The IWP1 alpha-0.05 results are incomplete.")
  }
  alpha_005$seed <- seed
  alpha_005$condition <- condition
  alpha_005 <- alpha_005[
    ,
    c("seed", "condition", setdiff(
      names(alpha_005),
      c("seed", "condition")
    )),
    drop = FALSE
  ]

  pi0 <- data.frame(
    seed = seed,
    condition = condition,
    method = methods,
    estimated_pi0 = c(
      constant_component_prior_weight(out$fash_fits$fash_iwp1_raw),
      constant_component_prior_weight(out$fash_fits$fash_iwp1_bf)
    ),
    stringsAsFactors = FALSE
  )
  list(alpha_005 = alpha_005, pi0 = pi0)
}

matrix_to_long <- function(correlation,
                           seed,
                           condition,
                           diagnostic) {
  n_time <- ncol(correlation)
  grid <- expand.grid(
    time_a = seq_len(n_time) - 1L,
    time_b = seq_len(n_time) - 1L,
    stringsAsFactors = FALSE
  )
  grid$seed <- seed
  grid$condition <- condition
  grid$diagnostic <- diagnostic
  grid$correlation <- as.vector(correlation)
  grid[, c(
    "seed", "condition", "diagnostic",
    "time_a", "time_b", "correlation"
  )]
}

make_lag_rows <- function(correlation,
                          seed,
                          condition,
                          diagnostic) {
  n_time <- ncol(correlation)
  lag_correlation <- vapply(seq_len(n_time - 1L), function(lag) {
    mean(correlation[cbind(
      seq_len(n_time - lag),
      (lag + 1L):n_time
    )])
  }, numeric(1))
  data.frame(
    seed = seed,
    condition = condition,
    diagnostic = diagnostic,
    lag = seq_along(lag_correlation),
    mean_correlation = lag_correlation,
    semivariogram = 1 - lag_correlation,
    stringsAsFactors = FALSE
  )
}

compute_full_matrix_diagnostics <- function(out,
                                            seed,
                                            condition,
                                            target_correlation) {
  beta_hat <- as.matrix(out$eqtl_summary$beta_hat)
  se <- as.matrix(out$eqtl_summary$se)
  true_beta <- as.matrix(out$true_beta)
  expression_errors <- extract_expression_errors(out)
  expression_error_matrix <- matrix(
    expression_errors,
    nrow = dim(expression_errors)[1] * dim(expression_errors)[2],
    ncol = dim(expression_errors)[3]
  )
  if (!identical(dim(beta_hat), dim(se)) ||
      !identical(dim(beta_hat), dim(true_beta)) ||
      any(!is.finite(beta_hat)) || any(!is.finite(se)) ||
      any(!is.finite(true_beta)) || any(se <= 0)) {
    stop("Invalid beta-hat, SE, or truth matrices in a full-matrix replicate.")
  }

  truth_known_error <- (beta_hat - true_beta) / se
  dynamic_null <- out$unit_info$effect_class %in% c("constant", "zero")
  null_beta_hat <- beta_hat[dynamic_null, , drop = FALSE]
  null_se <- se[dynamic_null, , drop = FALSE]
  weights <- 1 / null_se^2
  constant_estimate <- rowSums(null_beta_hat * weights) / rowSums(weights)
  centered_null_error <- sweep(
    null_beta_hat,
    1L,
    constant_estimate,
    `-`
  ) / null_se

  diagnostic_data <- list(
    "Expression error" = expression_error_matrix,
    "Truth-known standardized beta-hat error" = truth_known_error,
    "Centered dynamic-null beta-hat residual" = centered_null_error
  )
  correlations <- lapply(diagnostic_data, stats::cor)
  matrix_rows <- do.call(rbind, lapply(names(correlations), function(name) {
    matrix_to_long(correlations[[name]], seed, condition, name)
  }))
  lag_rows <- do.call(rbind, lapply(names(correlations), function(name) {
    make_lag_rows(correlations[[name]], seed, condition, name)
  }))
  match_rows <- do.call(rbind, lapply(names(correlations), function(name) {
    difference <- correlations[[name]] - target_correlation
    data.frame(
      seed = seed,
      condition = condition,
      diagnostic = name,
      maximum_absolute_difference = max(abs(difference)),
      frobenius_difference = sqrt(sum(difference^2)),
      target_lag1_correlation = mean(target_correlation[cbind(
        seq_len(ncol(target_correlation) - 1L),
        2:ncol(target_correlation)
      )]),
      realized_lag1_correlation = mean(correlations[[name]][cbind(
        seq_len(ncol(target_correlation) - 1L),
        2:ncol(target_correlation)
      )]),
      stringsAsFactors = FALSE
    )
  }))
  list(matrix_rows = matrix_rows, lag_rows = lag_rows, match_rows = match_rows)
}

summarize_alpha <- function(alpha_rows) {
  groups <- split(
    alpha_rows,
    list(alpha_rows$condition, alpha_rows$method),
    drop = TRUE
  )
  out <- do.call(rbind, lapply(groups, function(x) {
    power <- summarize_mc_values(x$power)
    fdr <- summarize_mc_values(x$empirical_fdr)
    discoveries <- summarize_mc_values(x$n_discoveries)
    data.frame(
      condition = x$condition[1],
      scenario = x$scenario[1],
      target = x$target[1],
      method = x$method[1],
      alpha = x$alpha[1],
      n_replications = length(unique(x$seed)),
      mean_discoveries = discoveries[["mean"]],
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
  }))
  rownames(out) <- NULL
  out
}

summarize_paired_alpha <- function(alpha_rows) {
  key <- c("seed", "scenario", "target", "method", "alpha")
  metrics <- c("n_discoveries", "empirical_fdr", "power")
  independent <- alpha_rows[
    alpha_rows$condition == "Independent",
    c(key, metrics),
    drop = FALSE
  ]
  current <- alpha_rows[
    alpha_rows$condition != "Independent",
    c("condition", key, metrics),
    drop = FALSE
  ]
  paired <- merge(
    current,
    independent,
    by = key,
    suffixes = c("_correlated", "_independent"),
    all = FALSE,
    sort = FALSE
  )
  if (nrow(paired) != nrow(current)) {
    stop("The full-matrix alpha results are not completely paired.")
  }
  for (metric in metrics) {
    paired[[paste0("delta_", metric)]] <-
      paired[[paste0(metric, "_correlated")]] -
      paired[[paste0(metric, "_independent")]]
  }
  groups <- split(
    paired,
    list(paired$condition, paired$method),
    drop = TRUE
  )
  out <- do.call(rbind, lapply(groups, function(x) {
    power <- summarize_mc_values(x$delta_power)
    fdr <- summarize_mc_values(x$delta_empirical_fdr)
    discoveries <- summarize_mc_values(x$delta_n_discoveries)
    data.frame(
      condition = x$condition[1],
      scenario = x$scenario[1],
      target = x$target[1],
      method = x$method[1],
      alpha = x$alpha[1],
      n_replications = length(unique(x$seed)),
      mean_delta_discoveries = discoveries[["mean"]],
      discoveries_delta_ci_lower = discoveries[["lower"]],
      discoveries_delta_ci_upper = discoveries[["upper"]],
      mean_delta_power = power[["mean"]],
      power_delta_ci_lower = power[["lower"]],
      power_delta_ci_upper = power[["upper"]],
      mean_delta_fdr = fdr[["mean"]],
      fdr_delta_ci_lower = fdr[["lower"]],
      fdr_delta_ci_upper = fdr[["upper"]],
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

summarize_pi0 <- function(pi0_rows) {
  groups <- split(
    pi0_rows,
    list(pi0_rows$condition, pi0_rows$method),
    drop = TRUE
  )
  out <- do.call(rbind, lapply(groups, function(x) {
    values <- summarize_mc_values(x$estimated_pi0)
    data.frame(
      condition = x$condition[1],
      method = x$method[1],
      n_replications = length(unique(x$seed)),
      mean_estimated_pi0 = values[["mean"]],
      pi0_sd = values[["sd"]],
      pi0_mc_se = values[["se"]],
      pi0_ci_lower = pmax(0, values[["lower"]]),
      pi0_ci_upper = pmin(1, values[["upper"]]),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

summarize_paired_pi0 <- function(pi0_rows) {
  independent <- pi0_rows[
    pi0_rows$condition == "Independent",
    c("seed", "method", "estimated_pi0"),
    drop = FALSE
  ]
  current <- pi0_rows[
    pi0_rows$condition != "Independent",
    ,
    drop = FALSE
  ]
  paired <- merge(
    current,
    independent,
    by = c("seed", "method"),
    suffixes = c("_correlated", "_independent"),
    all = FALSE,
    sort = FALSE
  )
  if (nrow(paired) != nrow(current)) {
    stop("The full-matrix pi0 results are not completely paired.")
  }
  paired$delta_pi0 <- paired$estimated_pi0_correlated -
    paired$estimated_pi0_independent
  groups <- split(
    paired,
    list(paired$condition, paired$method),
    drop = TRUE
  )
  out <- do.call(rbind, lapply(groups, function(x) {
    values <- summarize_mc_values(x$delta_pi0)
    data.frame(
      condition = x$condition[1],
      method = x$method[1],
      n_replications = length(unique(x$seed)),
      mean_delta_pi0 = values[["mean"]],
      pi0_delta_ci_lower = values[["lower"]],
      pi0_delta_ci_upper = values[["upper"]],
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

summarize_matrix_rows <- function(matrix_rows) {
  groups <- split(
    matrix_rows,
    list(
      matrix_rows$condition,
      matrix_rows$diagnostic,
      matrix_rows$time_a,
      matrix_rows$time_b
    ),
    drop = TRUE
  )
  out <- do.call(rbind, lapply(groups, function(x) {
    values <- summarize_mc_values(x$correlation)
    data.frame(
      condition = x$condition[1],
      diagnostic = x$diagnostic[1],
      time_a = x$time_a[1],
      time_b = x$time_b[1],
      n_replications = length(unique(x$seed)),
      mean_correlation = values[["mean"]],
      correlation_ci_lower = values[["lower"]],
      correlation_ci_upper = values[["upper"]],
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

summarize_lag_rows <- function(lag_rows) {
  groups <- split(
    lag_rows,
    list(lag_rows$condition, lag_rows$diagnostic, lag_rows$lag),
    drop = TRUE
  )
  out <- do.call(rbind, lapply(groups, function(x) {
    correlation <- summarize_mc_values(x$mean_correlation)
    semivariogram <- summarize_mc_values(x$semivariogram)
    data.frame(
      condition = x$condition[1],
      diagnostic = x$diagnostic[1],
      lag = x$lag[1],
      n_replications = length(unique(x$seed)),
      mean_correlation = correlation[["mean"]],
      correlation_ci_lower = correlation[["lower"]],
      correlation_ci_upper = correlation[["upper"]],
      mean_semivariogram = semivariogram[["mean"]],
      semivariogram_ci_lower = semivariogram[["lower"]],
      semivariogram_ci_upper = semivariogram[["upper"]],
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

validate_r1_reference <- function(alpha_rows,
                                  pi0_rows,
                                  workflowr_root,
                                  expected_seeds) {
  reference_dir <- file.path(
    workflowr_root,
    "output",
    "revision_simulations",
    "mc",
    "r1_random_bspline_main_effect_pilot5",
    "summary"
  )
  alpha_reference <- utils::read.csv(
    file.path(reference_dir, "all_replicate_alpha_curves.csv"),
    stringsAsFactors = FALSE
  )
  pi0_reference <- utils::read.csv(
    file.path(reference_dir, "all_replicate_pi0.csv"),
    stringsAsFactors = FALSE
  )
  methods <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
  alpha_reference <- alpha_reference[
    alpha_reference$seed %in% expected_seeds &
      alpha_reference$method %in% methods &
      abs(alpha_reference$alpha - 0.05) < 1e-12,
    ,
    drop = FALSE
  ]
  alpha_observed <- alpha_rows[alpha_rows$condition == "Independent", ]
  alpha_key <- c("seed", "method", "alpha")
  alpha_metrics <- c(
    "n_discoveries", "false_discoveries", "true_positives",
    "empirical_fdr", "power"
  )
  alpha_reference <- alpha_reference[
    do.call(order, unname(alpha_reference[alpha_key])),
    c(alpha_key, alpha_metrics),
    drop = FALSE
  ]
  alpha_observed <- alpha_observed[
    do.call(order, unname(alpha_observed[alpha_key])),
    c(alpha_key, alpha_metrics),
    drop = FALSE
  ]
  rownames(alpha_reference) <- rownames(alpha_observed) <- NULL

  pi0_reference <- pi0_reference[
    pi0_reference$seed %in% expected_seeds &
      pi0_reference$method == "FASH-IWP1",
    ,
    drop = FALSE
  ]
  pi0_reference$method <- ifelse(
    pi0_reference$fit == "Raw",
    "FASH-IWP1-Raw",
    "FASH-IWP1-BF"
  )
  pi0_observed <- pi0_rows[pi0_rows$condition == "Independent", ]
  pi0_key <- c("seed", "method")
  pi0_reference <- pi0_reference[
    do.call(order, unname(pi0_reference[pi0_key])),
    c(pi0_key, "estimated_pi0"),
    drop = FALSE
  ]
  pi0_observed <- pi0_observed[
    do.call(order, unname(pi0_observed[pi0_key])),
    c(pi0_key, "estimated_pi0"),
    drop = FALSE
  ]
  rownames(pi0_reference) <- rownames(pi0_observed) <- NULL

  if (!isTRUE(all.equal(alpha_reference[alpha_key], alpha_observed[alpha_key])) ||
      !isTRUE(all.equal(pi0_reference[pi0_key], pi0_observed[pi0_key]))) {
    stop("The independent full-matrix keys do not match the R1 cache.")
  }
  alpha_difference <- max(abs(
    as.matrix(alpha_reference[alpha_metrics]) -
      as.matrix(alpha_observed[alpha_metrics])
  ))
  pi0_difference <- max(abs(
    pi0_reference$estimated_pi0 - pi0_observed$estimated_pi0
  ))
  data.frame(
    component = c("R1 alpha-0.05 IWP1 metrics", "R1 IWP1 pi0"),
    maximum_absolute_difference = c(alpha_difference, pi0_difference),
    tolerance = 1e-12,
    passed = c(alpha_difference <= 1e-12, pi0_difference <= 1e-12),
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

J <- as.integer(get_arg("--J", "1000"))
n_donors <- as.integer(get_arg("--n-donors", "19"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
expression_noise_sd <- as.numeric(get_arg("--noise-sd", "1"))
seed_list <- parse_seed_list(get_arg(
  "--seed-list",
  "12345,22345,32345,42345,52345"
))
num_cores <- as.integer(get_arg("--num-cores", "4"))
num_basis <- as.integer(get_arg("--num-basis", "20"))
output_id <- get_arg(
  "--output-id",
  "r4_full_empirical_correlations_top500_pilot5"
)
matrix_cache <- get_arg(
  "--matrix-cache",
  file.path(
    workflowr_root,
    "output",
    "revision_simulations",
    "real_data",
    "r4_null_like_top500_full_correlations",
    "simulation_correlation_matrices.rds"
  )
)
overwrite <- as_flag(get_arg("--overwrite", "false"))
if (is.na(J) || J < 20L || is.na(n_donors) || is.na(n_covariates) ||
    n_covariates < 0L || n_donors < n_covariates + 3L ||
    !is.finite(expression_noise_sd) || expression_noise_sd <= 0 ||
    is.na(num_cores) || num_cores < 1L || is.na(num_basis) ||
    num_basis < 2L || !nzchar(output_id) || !file.exists(matrix_cache)) {
  stop("Invalid full empirical-correlation simulation arguments.")
}

time_grid <- make_time_grid()
loaded_matrices <- readRDS(matrix_cache)
required_matrix_names <- c("direct_centered", "pairwise_difference")
if (!identical(sort(names(loaded_matrices)), sort(required_matrix_names))) {
  stop("The matrix cache does not contain the two required full matrices.")
}
direct_correlation <- validate_time_correlation(
  loaded_matrices$direct_centered,
  n_time = length(time_grid)
)
pairwise_correlation <- validate_time_correlation(
  loaded_matrices$pairwise_difference,
  n_time = length(time_grid)
)
condition_matrices <- list(
  "Independent" = diag(length(time_grid)),
  "Direct centered full matrix" = direct_correlation,
  "Pairwise-difference full matrix" = pairwise_correlation
)
matrix_checks <- do.call(rbind, lapply(names(condition_matrices), function(name) {
  correlation <- condition_matrices[[name]]
  data.frame(
    condition = name,
    minimum_eigenvalue = min(eigen(
      correlation,
      symmetric = TRUE,
      only.values = TRUE
    )$values),
    maximum_diagonal_difference = max(abs(diag(correlation) - 1)),
    mean_lag1_correlation = mean(correlation[cbind(
      seq_len(ncol(correlation) - 1L),
      2:ncol(correlation)
    )]),
    mean_off_diagonal_correlation = mean(correlation[upper.tri(correlation)]),
    positive_definite = TRUE,
    stringsAsFactors = FALSE
  )
}))

class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
scenario <- "genotype_random_bspline_main_effect_dynamic_eqtl"
methods <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
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
  nominal_alpha = 0.05,
  methods = methods,
  conditions = names(condition_matrices),
  correlation_matrix_cache = matrix_cache,
  condition_matrices = condition_matrices,
  pairing_tolerance = 1e-10,
  r1_reference_output_id = "r1_random_bspline_main_effect_pilot5"
)

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "mc",
  output_id
)
replicate_dir <- file.path(output_dir, "replicates")
summary_dir <- file.path(output_dir, "summary")
invisible(lapply(
  c(output_dir, replicate_dir, summary_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))
configuration_path <- file.path(output_dir, "configuration.rds")
if (file.exists(configuration_path) && !overwrite) {
  cached_configuration <- readRDS(configuration_path)
  if (!isTRUE(all.equal(cached_configuration, configuration))) {
    stop("The existing output id has different settings. Use a new id or --overwrite true.")
  }
} else {
  saveRDS(configuration, configuration_path)
}

run_condition <- function(seed, condition) {
  message("Running full-matrix R4 seed ", seed, ": ", condition, ".")
  correlation <- if (condition == "Independent") {
    NULL
  } else {
    condition_matrices[[condition]]
  }
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
    expression_error_correlation = correlation
  )
  if (!isTRUE(all.equal(
    dynamic_null_proportion(out$unit_info, target = "dynamic"),
    configuration$true_pi0
  ))) {
    stop("The full-matrix dynamic-null proportion is incorrect.")
  }
  out
}

make_replicate <- function(seed) {
  independent <- run_condition(seed, "Independent")
  independent_errors <- extract_expression_errors(independent)
  independent_metrics <- extract_iwp_metrics(
    independent,
    seed,
    "Independent"
  )
  independent_diagnostics <- compute_full_matrix_diagnostics(
    independent,
    seed,
    "Independent",
    condition_matrices$Independent
  )
  alpha_parts <- list(independent_metrics$alpha_005)
  pi0_parts <- list(independent_metrics$pi0)
  matrix_parts <- list(independent_diagnostics$matrix_rows)
  lag_parts <- list(independent_diagnostics$lag_rows)
  match_parts <- list(independent_diagnostics$match_rows)
  pairing_parts <- list(data.frame(
    seed = seed,
    condition = "Independent",
    maximum_absolute_error_difference = 0,
    tolerance = configuration$pairing_tolerance,
    passed = TRUE,
    stringsAsFactors = FALSE
  ))

  for (condition in names(condition_matrices)[-1L]) {
    candidate <- run_condition(seed, condition)
    maximum_difference <- validate_against_independent(
      independent,
      candidate,
      correlation = condition_matrices[[condition]],
      independent_errors = independent_errors
    )
    metrics <- extract_iwp_metrics(candidate, seed, condition)
    diagnostics <- compute_full_matrix_diagnostics(
      candidate,
      seed,
      condition,
      condition_matrices[[condition]]
    )
    alpha_parts[[length(alpha_parts) + 1L]] <- metrics$alpha_005
    pi0_parts[[length(pi0_parts) + 1L]] <- metrics$pi0
    matrix_parts[[length(matrix_parts) + 1L]] <- diagnostics$matrix_rows
    lag_parts[[length(lag_parts) + 1L]] <- diagnostics$lag_rows
    match_parts[[length(match_parts) + 1L]] <- diagnostics$match_rows
    pairing_parts[[length(pairing_parts) + 1L]] <- data.frame(
      seed = seed,
      condition = condition,
      maximum_absolute_error_difference = maximum_difference,
      tolerance = configuration$pairing_tolerance,
      passed = maximum_difference <= configuration$pairing_tolerance,
      stringsAsFactors = FALSE
    )
    rm(candidate)
    invisible(gc(verbose = FALSE))
  }
  rm(independent)
  invisible(gc(verbose = FALSE))

  list(
    configuration = configuration,
    seed = seed,
    alpha_005 = do.call(rbind, alpha_parts),
    pi0 = do.call(rbind, pi0_parts),
    matrix_diagnostics = do.call(rbind, matrix_parts),
    lag_diagnostics = do.call(rbind, lag_parts),
    matrix_match = do.call(rbind, match_parts),
    pairing_check = do.call(rbind, pairing_parts)
  )
}

validate_replicate <- function(replicate, seed) {
  required_fields <- c(
    "configuration", "seed", "alpha_005", "pi0", "matrix_diagnostics",
    "lag_diagnostics", "matrix_match", "pairing_check"
  )
  if (!all(required_fields %in% names(replicate)) ||
      !identical(replicate$seed, seed) ||
      !isTRUE(all.equal(replicate$configuration, configuration))) {
    return(FALSE)
  }
  n_conditions <- length(condition_matrices)
  n_diagnostics <- 3L
  n_time <- length(time_grid)
  nrow(replicate$alpha_005) == n_conditions * length(methods) &&
    nrow(replicate$pi0) == n_conditions * length(methods) &&
    nrow(replicate$matrix_diagnostics) ==
      n_conditions * n_diagnostics * n_time^2 &&
    nrow(replicate$lag_diagnostics) ==
      n_conditions * n_diagnostics * (n_time - 1L) &&
    nrow(replicate$matrix_match) == n_conditions * n_diagnostics &&
    nrow(replicate$pairing_check) == n_conditions &&
    identical(sort(unique(replicate$alpha_005$condition)),
              sort(names(condition_matrices))) &&
    identical(sort(unique(replicate$alpha_005$method)), sort(methods)) &&
    !anyDuplicated(replicate$alpha_005[c("seed", "condition", "method")]) &&
    !anyDuplicated(replicate$pi0[c("seed", "condition", "method")]) &&
    all(replicate$pairing_check$passed) &&
    all(replicate$pairing_check$maximum_absolute_error_difference <=
      replicate$pairing_check$tolerance) &&
    all(is.finite(replicate$alpha_005$power)) &&
    all(is.finite(replicate$alpha_005$empirical_fdr)) &&
    all(is.finite(replicate$pi0$estimated_pi0)) &&
    all(is.finite(replicate$matrix_diagnostics$correlation))
}

replicates <- lapply(seed_list, function(seed) {
  replicate_path <- file.path(replicate_dir, paste0("seed_", seed, ".rds"))
  if (file.exists(replicate_path) && !overwrite) {
    cached <- readRDS(replicate_path)
    if (validate_replicate(cached, seed)) {
      message("Reusing full-matrix R4 cache: ", replicate_path)
      return(cached)
    }
    stop("Cached full-matrix replicate does not match: ", replicate_path)
  }
  replicate <- make_replicate(seed)
  if (!validate_replicate(replicate, seed)) {
    stop("The newly generated full-matrix replicate failed validation.")
  }
  saveRDS(replicate, replicate_path)
  replicate
})

all_alpha <- do.call(rbind, lapply(replicates, `[[`, "alpha_005"))
all_pi0 <- do.call(rbind, lapply(replicates, `[[`, "pi0"))
all_matrix <- do.call(rbind, lapply(replicates, `[[`, "matrix_diagnostics"))
all_lag <- do.call(rbind, lapply(replicates, `[[`, "lag_diagnostics"))
all_match <- do.call(rbind, lapply(replicates, `[[`, "matrix_match"))
pairing_check <- do.call(rbind, lapply(replicates, `[[`, "pairing_check"))
condition_alpha <- summarize_alpha(all_alpha)
paired_alpha <- summarize_paired_alpha(all_alpha)
condition_pi0 <- summarize_pi0(all_pi0)
paired_pi0 <- summarize_paired_pi0(all_pi0)
matrix_summary <- summarize_matrix_rows(all_matrix)
lag_summary <- summarize_lag_rows(all_lag)

default_reference_setting <-
  identical(J, 1000L) && identical(n_donors, 19L) &&
  identical(n_covariates, 5L) &&
  isTRUE(all.equal(expression_noise_sd, 1)) &&
  identical(seed_list, c(12345L, 22345L, 32345L, 42345L, 52345L)) &&
  identical(num_basis, 20L)
r1_reference_check <- if (default_reference_setting) {
  validate_r1_reference(
    all_alpha,
    all_pi0,
    workflowr_root,
    expected_seeds = seed_list
  )
} else {
  data.frame(
    component = c("R1 alpha-0.05 IWP1 metrics", "R1 IWP1 pi0"),
    maximum_absolute_difference = NA_real_,
    tolerance = NA_real_,
    passed = NA,
    stringsAsFactors = FALSE
  )
}

write_csv(all_alpha, file.path(summary_dir, "all_condition_alpha005.csv"))
write_csv(condition_alpha, file.path(summary_dir, "condition_mc_alpha005_summary.csv"))
write_csv(paired_alpha, file.path(summary_dir, "paired_alpha005_differences.csv"))
write_csv(all_pi0, file.path(summary_dir, "all_condition_pi0.csv"))
write_csv(condition_pi0, file.path(summary_dir, "condition_mc_pi0_summary.csv"))
write_csv(paired_pi0, file.path(summary_dir, "paired_pi0_differences.csv"))
write_csv(all_matrix, file.path(summary_dir, "all_realized_correlation_matrices.csv"))
write_csv(matrix_summary, file.path(summary_dir, "realized_correlation_matrix_summary.csv"))
write_csv(all_lag, file.path(summary_dir, "all_realized_lag_summaries.csv"))
write_csv(lag_summary, file.path(summary_dir, "realized_lag_summary.csv"))
write_csv(all_match, file.path(summary_dir, "matrix_match_diagnostics.csv"))
write_csv(pairing_check, file.path(summary_dir, "pairing_check.csv"))
write_csv(matrix_checks, file.path(summary_dir, "target_matrix_checks.csv"))
write_csv(r1_reference_check, file.path(summary_dir, "r1_reference_check.csv"))

cat("\nFull empirical-correlation IWP1 results at alpha = 0.05:\n")
print(condition_alpha)
cat("\nEstimated pi0:\n")
print(condition_pi0)
cat("\nMaximum pairing discrepancy:\n")
print(max(pairing_check$maximum_absolute_error_difference))
cat("\nR1 reference check:\n")
print(r1_reference_check)
