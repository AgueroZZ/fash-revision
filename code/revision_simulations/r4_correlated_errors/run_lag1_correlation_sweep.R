#!/usr/bin/env Rscript

# Run the paired R4 sensitivity sweep over lag-1-only error correlations.

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

parse_numeric_list <- function(x, argument_name) {
  pieces <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  if (length(pieces) == 0L || any(!nzchar(pieces))) {
    stop(argument_name, " must contain one or more comma-separated values.")
  }
  values <- suppressWarnings(as.numeric(pieces))
  if (any(!is.finite(values)) || anyDuplicated(values)) {
    stop(argument_name, " must contain unique finite numeric values.")
  }
  values
}

parse_seed_list <- function(x) {
  values <- parse_numeric_list(x, "--seed-list")
  seeds <- as.integer(values)
  if (any(values != seeds)) {
    stop("--seed-list must contain integer seeds.")
  }
  seeds
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

rho_text <- function(rho) {
  formatC(rho, format = "f", digits = 1)
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

validate_against_zero_reference <- function(reference,
                                            candidate,
                                            correlation,
                                            reference_errors) {
  paired_fields <- list(
    genotype = list(reference$genotype, candidate$genotype),
    covariates = list(reference$covariates, candidate$covariates),
    true_beta = list(reference$true_beta, candidate$true_beta),
    unit_info = list(reference$unit_info, candidate$unit_info),
    intercepts = list(
      reference$expression_simulation$intercepts,
      candidate$expression_simulation$intercepts
    ),
    covariate_effects = list(
      reference$expression_simulation$covariate_effects,
      candidate$expression_simulation$covariate_effects
    )
  )
  mismatched <- names(paired_fields)[!vapply(paired_fields, function(x) {
    identical(x[[1]], x[[2]])
  }, logical(1))]
  if (length(mismatched) > 0L) {
    stop(
      "The paired sweep simulations differ before the error transform: ",
      paste(mismatched, collapse = ", ")
    )
  }

  n_series <- dim(reference_errors)[1] * dim(reference_errors)[2]
  n_time <- dim(reference_errors)[3]
  expected_errors <- matrix(
    reference_errors,
    nrow = n_series,
    ncol = n_time
  ) %*% chol(correlation)
  candidate_errors <- matrix(
    extract_expression_errors(candidate),
    nrow = n_series,
    ncol = n_time
  )
  maximum_difference <- max(abs(expected_errors - candidate_errors))
  if (!is.finite(maximum_difference) || maximum_difference > 1e-10) {
    stop(
      "The sweep errors do not share the intended innovations; maximum difference = ",
      format(maximum_difference, scientific = TRUE)
    )
  }
  maximum_difference
}

mean_lag1_correlation <- function(x) {
  x <- as.matrix(x)
  if (nrow(x) < 3L || ncol(x) < 2L || any(!is.finite(x))) {
    stop("The lag-1 diagnostic input must be a finite matrix.")
  }
  correlation <- stats::cor(x)
  mean(correlation[cbind(
    seq_len(ncol(correlation) - 1L),
    2:ncol(correlation)
  )])
}

compute_lag1_diagnostics <- function(out, seed, rho) {
  beta_hat <- as.matrix(out$eqtl_summary$beta_hat)
  se <- as.matrix(out$eqtl_summary$se)
  true_beta <- as.matrix(out$true_beta)
  if (!identical(dim(beta_hat), dim(se)) ||
      !identical(dim(beta_hat), dim(true_beta)) ||
      any(!is.finite(beta_hat)) || any(!is.finite(se)) ||
      any(!is.finite(true_beta)) || any(se <= 0)) {
    stop("Invalid beta-hat, SE, or truth matrix in an R4 sweep fit.")
  }

  truth_known_error <- (beta_hat - true_beta) / se
  dynamic_null <- out$unit_info$effect_class %in% c("constant", "zero")
  null_beta_hat <- beta_hat[dynamic_null, , drop = FALSE]
  null_se <- se[dynamic_null, , drop = FALSE]
  if (nrow(null_beta_hat) < 3L) {
    stop("The R4 sweep fit has too few dynamic-null variants.")
  }
  weights <- 1 / null_se^2
  weighted_mean <- rowSums(null_beta_hat * weights) / rowSums(weights)
  centered_null_residual <- sweep(
    null_beta_hat,
    1,
    weighted_mean,
    `-`
  ) / null_se

  data.frame(
    seed = seed,
    rho = rho,
    diagnostic = c(
      "Truth-known standardized error",
      "Centered dynamic-null residual"
    ),
    mean_lag1_correlation = c(
      mean_lag1_correlation(truth_known_error),
      mean_lag1_correlation(centered_null_residual)
    ),
    n_units = c(nrow(truth_known_error), nrow(centered_null_residual)),
    stringsAsFactors = FALSE
  )
}

extract_alpha005 <- function(out, seed, rho, methods) {
  alpha_rows <- out$alpha_curve[
    out$alpha_curve$method %in% methods &
      abs(out$alpha_curve$alpha - 0.05) < 1e-12,
    ,
    drop = FALSE
  ]
  if (nrow(alpha_rows) != length(methods) ||
      !identical(sort(alpha_rows$method), sort(methods))) {
    stop("The alpha-0.05 sweep result is incomplete.")
  }
  alpha_rows$seed <- seed
  alpha_rows$rho <- rho
  alpha_rows <- alpha_rows[
    ,
    c("seed", "rho", setdiff(names(alpha_rows), c("seed", "rho"))),
    drop = FALSE
  ]
  rownames(alpha_rows) <- NULL
  alpha_rows
}

extract_iwp_pi0 <- function(out, seed, rho) {
  data.frame(
    seed = seed,
    rho = rho,
    method = c("FASH-IWP1-Raw", "FASH-IWP1-BF"),
    estimated_pi0 = c(
      constant_component_prior_weight(out$fash_fits$fash_iwp1_raw),
      constant_component_prior_weight(out$fash_fits$fash_iwp1_bf)
    ),
    stringsAsFactors = FALSE
  )
}

summarize_sweep_alpha <- function(alpha_rows) {
  groups <- split(
    alpha_rows,
    list(alpha_rows$method, alpha_rows$rho),
    drop = TRUE
  )
  out <- do.call(rbind, lapply(groups, function(x) {
    discoveries <- summarize_mc_values(x$n_discoveries)
    false_discoveries <- summarize_mc_values(x$false_discoveries)
    true_positives <- summarize_mc_values(x$true_positives)
    power <- summarize_mc_values(x$power)
    fdr <- summarize_mc_values(x$empirical_fdr)
    data.frame(
      rho = x$rho[1],
      scenario = x$scenario[1],
      target = x$target[1],
      method = x$method[1],
      alpha = x$alpha[1],
      n_replications = length(unique(x$seed)),
      mean_discoveries = discoveries[["mean"]],
      mean_false_discoveries = false_discoveries[["mean"]],
      mean_true_positives = true_positives[["mean"]],
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
  out$method_rank <- rank_revision_methods(out$method)
  out <- out[order(out$method_rank, out$method, out$rho), , drop = FALSE]
  out$method_rank <- NULL
  rownames(out) <- NULL
  out
}

summarize_sweep_pi0 <- function(pi0_rows) {
  groups <- split(
    pi0_rows,
    list(pi0_rows$method, pi0_rows$rho),
    drop = TRUE
  )
  out <- do.call(rbind, lapply(groups, function(x) {
    values <- summarize_mc_values(x$estimated_pi0)
    data.frame(
      rho = x$rho[1],
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
  out$method_rank <- rank_revision_methods(out$method)
  out <- out[order(out$method_rank, out$method, out$rho), , drop = FALSE]
  out$method_rank <- NULL
  rownames(out) <- NULL
  out
}

summarize_paired_vs_zero <- function(alpha_rows) {
  key <- c("seed", "scenario", "target", "method", "alpha")
  metrics <- c(
    "n_discoveries", "false_discoveries", "true_positives",
    "empirical_fdr", "power"
  )
  baseline <- alpha_rows[
    abs(alpha_rows$rho) < 1e-12,
    c(key, metrics),
    drop = FALSE
  ]
  current <- alpha_rows[, c("rho", key, metrics), drop = FALSE]
  paired <- merge(
    current,
    baseline,
    by = key,
    suffixes = c("_rho", "_zero"),
    all = FALSE,
    sort = FALSE
  )
  if (nrow(paired) != nrow(current)) {
    stop("The correlation-sweep results are not completely paired to rho zero.")
  }
  for (metric in metrics) {
    paired[[paste0("delta_", metric)]] <-
      paired[[paste0(metric, "_rho")]] -
      paired[[paste0(metric, "_zero")]]
  }

  groups <- split(
    paired,
    list(paired$method, paired$rho),
    drop = TRUE
  )
  out <- do.call(rbind, lapply(groups, function(x) {
    values <- lapply(metrics, function(metric) {
      summarize_mc_values(x[[paste0("delta_", metric)]])
    })
    names(values) <- metrics
    data.frame(
      rho = x$rho[1],
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
  }))
  out$method_rank <- rank_revision_methods(out$method)
  out <- out[order(out$method_rank, out$method, out$rho), , drop = FALSE]
  out$method_rank <- NULL
  rownames(out) <- NULL
  out
}

summarize_paired_pi0_vs_zero <- function(pi0_rows) {
  baseline <- pi0_rows[
    abs(pi0_rows$rho) < 1e-12,
    c("seed", "method", "estimated_pi0"),
    drop = FALSE
  ]
  current <- pi0_rows[, c(
    "rho", "seed", "method", "estimated_pi0"
  ), drop = FALSE]
  paired <- merge(
    current,
    baseline,
    by = c("seed", "method"),
    suffixes = c("_rho", "_zero"),
    all = FALSE,
    sort = FALSE
  )
  if (nrow(paired) != nrow(current)) {
    stop("The correlation-sweep pi0 results are not completely paired to rho zero.")
  }
  paired$delta_pi0 <- paired$estimated_pi0_rho - paired$estimated_pi0_zero
  groups <- split(
    paired,
    list(paired$method, paired$rho),
    drop = TRUE
  )
  out <- do.call(rbind, lapply(groups, function(x) {
    values <- summarize_mc_values(x$delta_pi0)
    data.frame(
      rho = x$rho[1],
      method = x$method[1],
      n_replications = length(unique(x$seed)),
      mean_delta_pi0 = values[["mean"]],
      pi0_delta_ci_lower = values[["lower"]],
      pi0_delta_ci_upper = values[["upper"]],
      stringsAsFactors = FALSE
    )
  }))
  out$method_rank <- rank_revision_methods(out$method)
  out <- out[order(out$method_rank, out$method, out$rho), , drop = FALSE]
  out$method_rank <- NULL
  rownames(out) <- NULL
  out
}

summarize_lag1_diagnostics <- function(lag_rows) {
  groups <- split(
    lag_rows,
    list(lag_rows$diagnostic, lag_rows$rho),
    drop = TRUE
  )
  out <- do.call(rbind, lapply(groups, function(x) {
    values <- summarize_mc_values(x$mean_lag1_correlation)
    data.frame(
      rho = x$rho[1],
      diagnostic = x$diagnostic[1],
      n_replications = length(unique(x$seed)),
      n_units = unique(x$n_units)[1],
      mean_lag1_correlation = values[["mean"]],
      lag1_sd = values[["sd"]],
      lag1_mc_se = values[["se"]],
      lag1_ci_lower = values[["lower"]],
      lag1_ci_upper = values[["upper"]],
      stringsAsFactors = FALSE
    )
  }))
  out <- out[order(out$diagnostic, out$rho), , drop = FALSE]
  rownames(out) <- NULL
  out
}

validate_zero_against_r1 <- function(alpha_rows,
                                     r1_cache_dir,
                                     methods,
                                     expected_seeds) {
  reference_path <- file.path(
    r1_cache_dir,
    "summary",
    "r1_all_replicate_fash_alpha_curves.csv"
  )
  if (!file.exists(reference_path)) {
    stop("The default R1 alpha-curve cache is missing.")
  }
  reference <- utils::read.csv(reference_path, stringsAsFactors = FALSE)
  reference <- reference[
    reference$method %in% methods &
      reference$seed %in% expected_seeds &
      abs(reference$alpha - 0.05) < 1e-12,
    ,
    drop = FALSE
  ]
  observed <- alpha_rows[abs(alpha_rows$rho) < 1e-12, , drop = FALSE]
  key <- c("seed", "method", "alpha")
  metrics <- c(
    "n_discoveries", "false_discoveries", "true_positives",
    "empirical_fdr", "power"
  )
  reference <- reference[
    do.call(order, unname(reference[key])),
    c(key, metrics),
    drop = FALSE
  ]
  observed <- observed[
    do.call(order, unname(observed[key])),
    c(key, metrics),
    drop = FALSE
  ]
  rownames(reference) <- NULL
  rownames(observed) <- NULL
  if (!isTRUE(all.equal(
    reference[key],
    observed[key],
    check.attributes = FALSE
  ))) {
    stop("The rho-zero sweep keys do not match the R1 cache.")
  }
  maximum_difference <- max(abs(
    as.matrix(reference[metrics]) - as.matrix(observed[metrics])
  ))
  if (!is.finite(maximum_difference) || maximum_difference > 1e-12) {
    stop("The rho-zero sweep metrics do not reproduce the R1 cache.")
  }
  data.frame(
    component = "rho-zero alpha-0.05 metrics",
    maximum_absolute_difference = maximum_difference,
    tolerance = 1e-12,
    passed = TRUE,
    stringsAsFactors = FALSE
  )
}

validate_zero_pi0_against_r1 <- function(pi0_rows,
                                         r1_cache_dir,
                                         expected_seeds) {
  reference_path <- file.path(
    r1_cache_dir,
    "summary",
    "r1_all_replicate_pi0.csv"
  )
  if (!file.exists(reference_path)) {
    stop("The default R1 pi0 cache is missing.")
  }
  reference <- utils::read.csv(reference_path, stringsAsFactors = FALSE)
  reference <- reference[
    reference$seed %in% expected_seeds & reference$method == "FASH-IWP1",
    ,
    drop = FALSE
  ]
  reference$method <- ifelse(
    reference$fit == "Raw",
    "FASH-IWP1-Raw",
    "FASH-IWP1-BF"
  )
  observed <- pi0_rows[abs(pi0_rows$rho) < 1e-12, , drop = FALSE]
  key <- c("seed", "method")
  reference <- reference[
    do.call(order, unname(reference[key])),
    c(key, "estimated_pi0"),
    drop = FALSE
  ]
  observed <- observed[
    do.call(order, unname(observed[key])),
    c(key, "estimated_pi0"),
    drop = FALSE
  ]
  rownames(reference) <- rownames(observed) <- NULL
  if (!isTRUE(all.equal(
    reference[key],
    observed[key],
    check.attributes = FALSE
  ))) {
    stop("The rho-zero pi0 keys do not match the R1 cache.")
  }
  maximum_difference <- max(abs(
    reference$estimated_pi0 - observed$estimated_pi0
  ))
  if (!is.finite(maximum_difference) || maximum_difference > 1e-12) {
    stop("The rho-zero pi0 estimates do not reproduce the R1 cache.")
  }
  data.frame(
    component = "rho-zero IWP1 pi0",
    maximum_absolute_difference = maximum_difference,
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
  "shared",
  "real_genotype_one_per_gene.R"
))
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "r4_correlated_errors",
  "real_genotype_r1_helpers.R"
))
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "r4_correlated_errors",
  "plotting.R"
))

J <- as.integer(get_arg("--J", "6362"))
n_donors <- as.integer(get_arg("--n-donors", "19"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
expression_noise_sd <- as.numeric(get_arg("--noise-sd", "1"))
rho_grid <- parse_numeric_list(
  get_arg("--rho-list", "-0.3,-0.2,-0.1,0,0.1,0.2,0.3"),
  "--rho-list"
)
rho_grid[abs(rho_grid) < 1e-12] <- 0
rho_grid <- sort(rho_grid)
seed_list <- parse_seed_list(get_arg(
  "--seed-list",
  "12345,22345,32345,42345,52345"
))
num_cores <- as.integer(get_arg("--num-cores", "4"))
num_basis <- as.integer(get_arg("--num-basis", "20"))
output_id <- get_arg(
  "--output-id",
  "r4_lag1_correlation_sweep_J6362_m0p3_to_p0p3_fashr0143_pilot5"
)
genotype_cache_path <- get_arg(
  "--genotype-cache",
  file.path(
    workflowr_root,
    "output", "revision_simulations", "shared",
    "real_genotype_one_per_gene_J6362_pilot5",
    "genotype_samples.rds"
  )
)
r1_cache_dir <- get_arg(
  "--r1-cache-dir",
  file.path(
    workflowr_root,
    "output", "revision_simulations", "mc", "r1_r2_fashr0143"
  )
)
output_dir_argument <- get_arg("--output-dir", "")
expected_fashr_version <- get_arg("--expected-fashr-version", "")
expected_fashr_remote_sha <- get_arg("--expected-fashr-remote-sha", "")
overwrite <- as_flag(get_arg("--overwrite", "false"))

if (is.na(J) || J < 10L || is.na(n_donors) ||
    is.na(n_covariates) || n_covariates < 0L ||
    n_donors < n_covariates + 3L ||
    !is.finite(expression_noise_sd) || expression_noise_sd <= 0 ||
    is.na(num_cores) || num_cores < 1L ||
    is.na(num_basis) || num_basis < 2L ||
    length(rho_grid) < 3L || anyDuplicated(rho_grid) ||
    !any(rho_grid == 0) || !nzchar(output_id) ||
    !file.exists(genotype_cache_path) || !dir.exists(r1_cache_dir)) {
  stop("Invalid R4 correlation-sweep arguments.")
}
genotype_cache_path <- normalizePath(
  genotype_cache_path, winslash = "/", mustWork = TRUE
)
r1_cache_dir <- normalizePath(r1_cache_dir, winslash = "/", mustWork = TRUE)

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required for the R4 simulation.")
}
fashr_description <- utils::packageDescription("fashr")
package_provenance <- list(
  package = "fashr",
  version = as.character(utils::packageVersion("fashr")),
  remote_sha = if (is.null(fashr_description$RemoteSha)) {
    NA_character_
  } else {
    as.character(fashr_description$RemoteSha)
  },
  library_path = normalizePath(
    find.package("fashr"), winslash = "/", mustWork = TRUE
  ),
  r_version = R.version.string,
  platform = R.version$platform
)
if (nzchar(expected_fashr_version) &&
    !identical(package_provenance$version, expected_fashr_version)) {
  stop(
    "Expected fashr ", expected_fashr_version,
    "; found ", package_provenance$version, "."
  )
}
if (nzchar(expected_fashr_remote_sha) &&
    !identical(package_provenance$remote_sha, expected_fashr_remote_sha)) {
  stop(
    "Expected fashr RemoteSha ", expected_fashr_remote_sha,
    "; found ", package_provenance$remote_sha, "."
  )
}
genotype_cache <- validate_r4_real_genotype_cache(
  readRDS(genotype_cache_path),
  seed_list = seed_list,
  J = J,
  n_donors = n_donors
)

time_grid <- make_time_grid()
correlation_matrices <- lapply(rho_grid, function(rho) {
  validate_time_correlation(
    make_lag1_correlation(length(time_grid), rho),
    n_time = length(time_grid)
  )
})
names(correlation_matrices) <- vapply(rho_grid, rho_text, character(1))
correlation_checks <- do.call(rbind, lapply(seq_along(rho_grid), function(i) {
  correlation <- correlation_matrices[[i]]
  data.frame(
    rho = rho_grid[i],
    minimum_eigenvalue = min(eigen(
      correlation,
      symmetric = TRUE,
      only.values = TRUE
    )$values),
    maximum_diagonal_difference = max(abs(diag(correlation) - 1)),
    maximum_higher_lag_correlation = max(abs(
      correlation[abs(row(correlation) - col(correlation)) > 1L]
    )),
    positive_definite = TRUE,
    stringsAsFactors = FALSE
  )
}))

class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
scenario <- paste0(
  "r1_real_genotype_one_per_gene_",
  "random_bspline_main_effect_dynamic_eqtl"
)
expected_class_counts <- exact_proportional_counts(J, class_probs)
true_pi0 <- unname(
  sum(expected_class_counts[c("constant", "zero")]) / J
)
fash_methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "FASH-linear-Raw",
  "FASH-linear-BF"
)
pi0_methods <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
nominal_alpha <- 0.05
empirical_generating_rho <- -0.09

output_dir <- if (nzchar(output_dir_argument)) {
  output_dir_argument
} else {
  file.path(
    workflowr_root,
    "output",
    "revision_simulations",
    "mc",
    output_id
  )
}
replicate_dir <- file.path(output_dir, "replicates")
summary_dir <- file.path(output_dir, "summary")
figure_dir <- file.path(output_dir, "figures")
invisible(lapply(
  c(output_dir, replicate_dir, summary_dir, figure_dir),
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
  expected_class_counts = expected_class_counts,
  dynamic_amplitude = 2,
  bspline_df = 6,
  bspline_coefficient_sd = 1,
  true_pi0 = true_pi0,
  seed_list = seed_list,
  rho_grid = rho_grid,
  correlation_structure = "lag1_only",
  nominal_alpha = nominal_alpha,
  empirical_generating_rho = empirical_generating_rho,
  pi0_methods = pi0_methods,
  genotype_cache_path = genotype_cache_path,
  genotype_cache_fingerprint = artifact_fingerprint(genotype_cache_path),
  genotype_source = "paper-derived YRI DS dosage",
  genotype_selection_rule = genotype_cache$configuration$selection_rule,
  maf_truth_balance_method =
    "exact global class counts with within-MAF-decile permutation",
  package_provenance = package_provenance,
  pairing_tolerance = 1e-10,
  r1_reference_output_id = "r1_r2_fashr0143",
  r1_reference_cache_dir = r1_cache_dir
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

run_condition <- function(seed_inputs, rho) {
  message(
    "Running R4 sweep seed ", seed_inputs$seed,
    ", rho = ", rho_text(rho), "."
  )
  correlation <- if (rho == 0) {
    NULL
  } else {
    correlation_matrices[[rho_text(rho)]]
  }
  out <- run_r4_real_genotype_r1_condition(
    seed_inputs = seed_inputs,
    expression_error_correlation = correlation,
    num_cores = num_cores,
    num_basis = num_basis,
    nominal_alpha = nominal_alpha,
    output_dir = output_dir,
    verbose = FALSE
  )
  if (!isTRUE(all.equal(
    dynamic_null_proportion(out$unit_info, target = "dynamic"),
    configuration$true_pi0
  ))) {
    stop("The R4 sweep dynamic-null proportion is incorrect.")
  }
  out
}

make_replicate <- function(seed) {
  execution_rhos <- c(0, rho_grid[rho_grid != 0])
  seed_inputs <- prepare_r4_real_genotype_r1_seed(
    genotype_cache = genotype_cache,
    seed = seed,
    J = J,
    n_donors = n_donors,
    n_covariates = n_covariates,
    time_grid = time_grid,
    class_probs = class_probs,
    expression_noise_sd = expression_noise_sd,
    dynamic_main_effect_sd = 1
  )
  reference <- run_condition(seed_inputs, 0)
  reference_errors <- extract_expression_errors(reference)
  alpha_parts <- list(extract_alpha005(
    reference,
    seed = seed,
    rho = 0,
    methods = fash_methods
  ))
  pi0_parts <- list(extract_iwp_pi0(reference, seed = seed, rho = 0))
  lag_parts <- list(compute_lag1_diagnostics(reference, seed, 0))
  pairing_parts <- list(data.frame(
    seed = seed,
    rho = 0,
    maximum_absolute_error_difference = 0,
    tolerance = configuration$pairing_tolerance,
    passed = TRUE,
    stringsAsFactors = FALSE
  ))

  for (rho in execution_rhos[execution_rhos != 0]) {
    candidate <- run_condition(seed_inputs, rho)
    maximum_difference <- validate_against_zero_reference(
      reference,
      candidate,
      correlation = correlation_matrices[[rho_text(rho)]],
      reference_errors = reference_errors
    )
    alpha_parts[[length(alpha_parts) + 1L]] <- extract_alpha005(
      candidate,
      seed = seed,
      rho = rho,
      methods = fash_methods
    )
    pi0_parts[[length(pi0_parts) + 1L]] <- extract_iwp_pi0(
      candidate,
      seed = seed,
      rho = rho
    )
    lag_parts[[length(lag_parts) + 1L]] <- compute_lag1_diagnostics(
      candidate,
      seed,
      rho
    )
    pairing_parts[[length(pairing_parts) + 1L]] <- data.frame(
      seed = seed,
      rho = rho,
      maximum_absolute_error_difference = maximum_difference,
      tolerance = configuration$pairing_tolerance,
      passed = maximum_difference <= configuration$pairing_tolerance,
      stringsAsFactors = FALSE
    )
    rm(candidate)
    invisible(gc(verbose = FALSE))
  }

  alpha_rows <- do.call(rbind, alpha_parts)
  pi0_rows <- do.call(rbind, pi0_parts)
  lag_rows <- do.call(rbind, lag_parts)
  pairing_rows <- do.call(rbind, pairing_parts)
  alpha_rows <- alpha_rows[order(alpha_rows$rho, alpha_rows$method), ]
  pi0_rows <- pi0_rows[order(pi0_rows$rho, pi0_rows$method), ]
  lag_rows <- lag_rows[order(lag_rows$rho, lag_rows$diagnostic), ]
  pairing_rows <- pairing_rows[order(pairing_rows$rho), ]
  rownames(alpha_rows) <- NULL
  rownames(pi0_rows) <- NULL
  rownames(lag_rows) <- NULL
  rownames(pairing_rows) <- NULL

  list(
    configuration = configuration,
    seed = seed,
    genotype_digest = seed_inputs$genotype_sample$genotype_digest,
    selected_pair_keys = seed_inputs$genotype_sample$selection$pair_key,
    alpha_005 = alpha_rows,
    pi0 = pi0_rows,
    lag1_diagnostics = lag_rows,
    pairing_check = pairing_rows
  )
}

validate_replicate <- function(replicate, seed) {
  required_fields <- c(
    "configuration", "seed", "alpha_005", "pi0", "lag1_diagnostics",
    "pairing_check", "genotype_digest", "selected_pair_keys"
  )
  if (!all(required_fields %in% names(replicate)) ||
      !identical(replicate$seed, seed) ||
      !isTRUE(all.equal(replicate$configuration, configuration)) ||
      !identical(
        replicate$genotype_digest,
        genotype_cache$samples[[as.character(seed)]]$genotype_digest
      ) ||
      !identical(
        replicate$selected_pair_keys,
        genotype_cache$samples[[as.character(seed)]]$selection$pair_key
      )) {
    return(FALSE)
  }
  expected_alpha_rows <- length(rho_grid) * length(fash_methods)
  expected_pi0_rows <- length(rho_grid) * length(pi0_methods)
  expected_lag_rows <- length(rho_grid) * 2L
  nrow(replicate$alpha_005) == expected_alpha_rows &&
    nrow(replicate$pi0) == expected_pi0_rows &&
    nrow(replicate$lag1_diagnostics) == expected_lag_rows &&
    nrow(replicate$pairing_check) == length(rho_grid) &&
    identical(sort(unique(replicate$alpha_005$rho)), rho_grid) &&
    identical(sort(unique(replicate$alpha_005$method)), sort(fash_methods)) &&
    identical(sort(unique(replicate$pi0$rho)), rho_grid) &&
    identical(sort(unique(replicate$pi0$method)), sort(pi0_methods)) &&
    !anyDuplicated(replicate$alpha_005[c("seed", "rho", "method")]) &&
    !anyDuplicated(replicate$pi0[c("seed", "rho", "method")]) &&
    all(replicate$pairing_check$passed) &&
    all(replicate$pairing_check$maximum_absolute_error_difference <=
      replicate$pairing_check$tolerance) &&
    all(is.finite(replicate$alpha_005$power)) &&
    all(is.finite(replicate$alpha_005$empirical_fdr)) &&
    all(is.finite(replicate$pi0$estimated_pi0)) &&
    all(is.finite(replicate$lag1_diagnostics$mean_lag1_correlation))
}

replicates <- lapply(seed_list, function(seed) {
  replicate_path <- file.path(replicate_dir, paste0("seed_", seed, ".rds"))
  if (file.exists(replicate_path) && !overwrite) {
    cached <- readRDS(replicate_path)
    if (validate_replicate(cached, seed)) {
      message("Reusing R4 sweep replicate cache: ", replicate_path)
      return(cached)
    }
    stop("Cached R4 sweep replicate does not match: ", replicate_path)
  }
  replicate <- make_replicate(seed)
  if (!validate_replicate(replicate, seed)) {
    stop("The newly generated R4 sweep replicate failed validation.")
  }
  saveRDS(replicate, replicate_path)
  replicate
})

all_alpha_005 <- do.call(rbind, lapply(replicates, `[[`, "alpha_005"))
all_pi0 <- do.call(rbind, lapply(replicates, `[[`, "pi0"))
all_lag1 <- do.call(rbind, lapply(replicates, `[[`, "lag1_diagnostics"))
pairing_check <- do.call(rbind, lapply(replicates, `[[`, "pairing_check"))
sweep_summary <- summarize_sweep_alpha(all_alpha_005)
paired_summary <- summarize_paired_vs_zero(all_alpha_005)
pi0_summary <- summarize_sweep_pi0(all_pi0)
paired_pi0_summary <- summarize_paired_pi0_vs_zero(all_pi0)
lag1_summary <- summarize_lag1_diagnostics(all_lag1)

default_reference_setting <-
  identical(J, 6362L) &&
  identical(n_donors, 19L) &&
  identical(n_covariates, 5L) &&
  isTRUE(all.equal(expression_noise_sd, 1)) &&
  identical(seed_list, c(12345L, 22345L, 32345L, 42345L, 52345L)) &&
  identical(num_basis, 20L)
r1_reference_check <- if (default_reference_setting) {
  validate_zero_against_r1(
    all_alpha_005,
    r1_cache_dir = r1_cache_dir,
    methods = fash_methods,
    expected_seeds = seed_list
  )
} else {
  data.frame(
    component = "rho-zero R1 cache check",
    maximum_absolute_difference = NA_real_,
    tolerance = NA_real_,
    passed = NA,
    stringsAsFactors = FALSE
  )
}
r1_pi0_reference_check <- if (default_reference_setting) {
  validate_zero_pi0_against_r1(
    all_pi0,
    r1_cache_dir = r1_cache_dir,
    expected_seeds = seed_list
  )
} else {
  data.frame(
    component = "rho-zero IWP1 pi0",
    maximum_absolute_difference = NA_real_,
    tolerance = NA_real_,
    passed = NA,
    stringsAsFactors = FALSE
  )
}

write_csv(all_alpha_005, file.path(summary_dir, "all_alpha005.csv"))
write_csv(sweep_summary, file.path(summary_dir, "mc_alpha005_summary.csv"))
write_csv(
  paired_summary,
  file.path(summary_dir, "paired_vs_zero_alpha005_summary.csv")
)
write_csv(all_pi0, file.path(summary_dir, "all_pi0.csv"))
write_csv(pi0_summary, file.path(summary_dir, "mc_pi0_summary.csv"))
write_csv(
  paired_pi0_summary,
  file.path(summary_dir, "paired_vs_zero_pi0_summary.csv")
)
write_csv(all_lag1, file.path(summary_dir, "all_lag1_correlations.csv"))
write_csv(lag1_summary, file.path(summary_dir, "lag1_correlation_summary.csv"))
write_csv(pairing_check, file.path(summary_dir, "pairing_check.csv"))
write_csv(
  correlation_checks,
  file.path(summary_dir, "correlation_matrix_check.csv")
)
write_csv(
  r1_reference_check,
  file.path(summary_dir, "r1_reference_check.csv")
)
write_csv(
  r1_pi0_reference_check,
  file.path(summary_dir, "r1_pi0_reference_check.csv")
)

plot_r4_correlation_sweep(
  sweep_summary,
  file = file.path(figure_dir, "correlation_sweep_alpha005.png"),
  empirical_rho = empirical_generating_rho,
  nominal_alpha = nominal_alpha,
  methods = pi0_methods
)
plot_r4_correlation_sweep_pi0(
  pi0_summary,
  file = file.path(figure_dir, "correlation_sweep_pi0.png"),
  empirical_rho = empirical_generating_rho,
  true_pi0 = configuration$true_pi0
)

cat("\nR4 lag-1 correlation sweep at alpha = 0.05:\n")
print(sweep_summary[sweep_summary$method == "FASH-IWP1-BF", ])
cat("\nR4 lag-1 correlation sweep estimated pi0:\n")
print(pi0_summary)
cat("\nRealized truth-known lag-1 correlations:\n")
print(lag1_summary[
  lag1_summary$diagnostic == "Truth-known standardized error",
  c("rho", "mean_lag1_correlation", "lag1_ci_lower", "lag1_ci_upper")
])
cat("\nMaximum pairing error:\n")
print(max(pairing_check$maximum_absolute_error_difference))
cat("\nR1 rho-zero reference check:\n")
print(r1_reference_check)
print(r1_pi0_reference_check)
