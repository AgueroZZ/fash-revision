#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr project root.")
}

mean_ci <- function(x, level = 0.95) {
  estimate <- mean(x)
  standard_error <- stats::sd(x) / sqrt(length(x))
  multiplier <- stats::qt((1 + level) / 2, df = length(x) - 1L)
  c(
    estimate = estimate,
    lower = estimate - multiplier * standard_error,
    upper = estimate + multiplier * standard_error
  )
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "r0_linear_iwp_power_mechanism", "power_mechanism_helpers.R"
))

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

J <- 4000L
pi0 <- 0.95
time_grid <- 0:15
endpoint_scale <- 3
standard_error <- 1
equicorrelations <- c(0, 0.2, 0.4, 0.6)
confirmation_seeds <- seq.int(920001L, by = 1000L, length.out = 20L)
grid <- default_revision_grid()
alpha <- 0.05
num_basis <- 20L

screen_cache_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r0_linear_iwp_power_mechanism_error_correlation_screen",
  "analysis_cache.rds"
)
screen_cache <- readRDS(screen_cache_path)
if (screen_cache$configuration$screen_seed %in% confirmation_seeds ||
    !all(equicorrelations %in%
      screen_cache$configuration$screen_grid$equicorrelation) ||
    !endpoint_scale %in%
      screen_cache$configuration$screen_grid$endpoint_scale) {
  stop("The held-out confirmation is not separated from the screen.")
}

linear_kernel <- tcrossprod(time_grid - min(time_grid))
iwp_kernel <- compute_iwp1_kernel(time_grid, num_basis)

run_one <- function(equicorrelation, seed) {
  simulated <- simulate_linear_power_mechanism(
    J = J,
    pi0 = pi0,
    time_grid = time_grid,
    endpoint_scale = endpoint_scale,
    standard_error = standard_error,
    alternative_distribution = "fixed",
    equicorrelation = equicorrelation,
    seed = seed
  )
  linear_likelihood <- compute_common_intercept_log_likelihood(
    simulated$y_matrix,
    standard_error,
    linear_kernel,
    grid
  )
  iwp_likelihood <- compute_common_intercept_log_likelihood(
    simulated$y_matrix,
    standard_error,
    iwp_kernel,
    grid
  )
  result <- rbind(
    fit_power_mechanism_family(
      linear_likelihood,
      grid,
      simulated$true_null,
      "FASH-linear",
      alpha
    ),
    fit_power_mechanism_family(
      iwp_likelihood,
      grid,
      simulated$true_null,
      "FASH-IWP1",
      alpha
    )
  )
  result$equicorrelation <- equicorrelation
  result$seed <- seed
  result
}

start_time <- proc.time()[["elapsed"]]
seed_results <- do.call(rbind, lapply(equicorrelations, function(rho) {
  do.call(rbind, lapply(confirmation_seeds, function(seed) {
    run_one(rho, seed)
  }))
}))
elapsed_seconds <- proc.time()[["elapsed"]] - start_time
rownames(seed_results) <- NULL

bf_results <- seed_results[seed_results$adjustment == "BF", , drop = FALSE]
paired <- merge(
  bf_results[bf_results$family == "FASH-linear", ],
  bf_results[bf_results$family == "FASH-IWP1", ],
  by = c("equicorrelation", "seed", "adjustment"),
  suffixes = c("_linear", "_iwp"),
  sort = TRUE
)
paired$power_difference <- paired$power_linear - paired$power_iwp
paired$fdp_difference <-
  paired$realized_fdp_linear - paired$realized_fdp_iwp
paired$pi0_difference <-
  paired$estimated_pi0_linear - paired$estimated_pi0_iwp

confirmation_summary <- do.call(rbind, lapply(
  split(paired, paired$equicorrelation),
  function(x) {
    linear_ci <- mean_ci(x$power_linear)
    iwp_ci <- mean_ci(x$power_iwp)
    difference_ci <- mean_ci(x$power_difference)
    data.frame(
      equicorrelation = x$equicorrelation[[1L]],
      n_seeds = nrow(x),
      mean_power_linear = linear_ci[["estimate"]],
      mean_power_linear_lower = linear_ci[["lower"]],
      mean_power_linear_upper = linear_ci[["upper"]],
      mean_power_iwp = iwp_ci[["estimate"]],
      mean_power_iwp_lower = iwp_ci[["lower"]],
      mean_power_iwp_upper = iwp_ci[["upper"]],
      mean_power_difference = difference_ci[["estimate"]],
      power_difference_lower = difference_ci[["lower"]],
      power_difference_upper = difference_ci[["upper"]],
      mean_power_ratio = mean(x$power_linear) / mean(x$power_iwp),
      empirical_fdr_linear = mean(x$realized_fdp_linear),
      empirical_fdr_iwp = mean(x$realized_fdp_iwp),
      mean_pi0_linear = mean(x$estimated_pi0_linear),
      mean_pi0_iwp = mean(x$estimated_pi0_iwp),
      stringsAsFactors = FALSE
    )
  }
))
rownames(confirmation_summary) <- NULL
confirmation_summary <- confirmation_summary[
  match(equicorrelations, confirmation_summary$equicorrelation),
  ,
  drop = FALSE
]

validation <- data.frame(
  check = c(
    "confirmation seeds are held out",
    "complete seed-method results",
    "all BF fits are available",
    "finite confirmation estimates",
    "one-job thread cap"
  ),
  passed = c(
    !screen_cache$configuration$screen_seed %in% confirmation_seeds,
    nrow(seed_results) ==
      length(equicorrelations) * length(confirmation_seeds) * 4L,
    all(bf_results$bf_available),
    all(is.finite(as.matrix(confirmation_summary[
      vapply(confirmation_summary, is.numeric, logical(1))
    ]))),
    Sys.getenv("VECLIB_MAXIMUM_THREADS") == "1"
  ),
  stringsAsFactors = FALSE
)
if (any(!validation$passed)) {
  print(validation)
  stop("The error-correlation confirmation failed validation.")
}

output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r0_linear_iwp_power_mechanism_error_correlation_confirmation"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
analysis_cache <- list(
  configuration = list(
    scope = "internal held-out equicorrelated-error confirmation",
    J = J,
    pi0 = pi0,
    time_grid = time_grid,
    endpoint_scale = endpoint_scale,
    standard_error = standard_error,
    equicorrelations = equicorrelations,
    confirmation_seeds = confirmation_seeds,
    working_likelihood = "diagonal standard errors",
    truth_error_model = "equicorrelated Gaussian errors",
    grid = grid,
    penalty = 10L,
    alpha = alpha,
    screen_cache_path = normalizePath(screen_cache_path)
  ),
  seed_results = seed_results,
  paired = paired,
  confirmation_summary = confirmation_summary,
  validation = validation,
  elapsed_seconds = elapsed_seconds
)
saveRDS(analysis_cache, file.path(output_directory, "analysis_cache.rds"))
utils::write.csv(
  seed_results,
  file.path(output_directory, "seed_results.csv"),
  row.names = FALSE
)
utils::write.csv(
  paired,
  file.path(output_directory, "paired_seed_results.csv"),
  row.names = FALSE
)
utils::write.csv(
  confirmation_summary,
  file.path(output_directory, "confirmation_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  validation,
  file.path(output_directory, "validation.csv"),
  row.names = FALSE
)

print(confirmation_summary)
print(validation)
message("Elapsed seconds: ", round(elapsed_seconds, 3))
message("Saved output to: ", normalizePath(output_directory))
