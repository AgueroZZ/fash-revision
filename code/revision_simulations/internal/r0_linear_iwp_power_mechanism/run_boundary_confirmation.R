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
  x <- as.numeric(x)
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
standard_error <- 1
grid <- default_revision_grid()
alpha <- 0.05
num_basis <- 20L
confirmation_seeds <- seq.int(910001L, by = 1000L, length.out = 20L)
candidate_settings <- data.frame(
  candidate = c(
    "fixed_relative_boundary",
    "fixed_absolute_boundary",
    "narrow_relative_boundary",
    "gaussian_same_scale"
  ),
  alternative_distribution = c(
    "fixed",
    "fixed",
    "narrow_normal",
    "gaussian"
  ),
  endpoint_scale = c(2.25, 3.00, 2.25, 2.25),
  stringsAsFactors = FALSE
)

screen_cache_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r0_linear_iwp_power_mechanism_boundary_screen",
  "analysis_cache.rds"
)
screen_cache <- readRDS(screen_cache_path)
if (screen_cache$configuration$screen_seed %in% confirmation_seeds ||
    !all(c(2.25, 3.00) %in%
      screen_cache$configuration$endpoint_changes) ||
    !16L %in% screen_cache$configuration$time_counts ||
    !0.95 %in% screen_cache$configuration$pi0_values) {
  stop("The confirmation settings are not separated from the screen.")
}

linear_kernel <- tcrossprod(time_grid - min(time_grid))
iwp_kernel <- compute_iwp1_kernel(time_grid, num_basis)

run_one <- function(candidate,
                    alternative_distribution,
                    endpoint_scale,
                    seed) {
  simulated <- simulate_linear_power_mechanism(
    J = J,
    pi0 = pi0,
    time_grid = time_grid,
    endpoint_scale = endpoint_scale,
    standard_error = standard_error,
    alternative_distribution = alternative_distribution,
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
  result$candidate <- candidate
  result$alternative_distribution <- alternative_distribution
  result$endpoint_scale <- endpoint_scale
  result$seed <- seed
  result
}

start_time <- proc.time()[["elapsed"]]
seed_results <- lapply(seq_len(nrow(candidate_settings)), function(index) {
  setting <- candidate_settings[index, ]
  do.call(rbind, lapply(confirmation_seeds, function(seed) {
    run_one(
      candidate = setting$candidate,
      alternative_distribution = setting$alternative_distribution,
      endpoint_scale = setting$endpoint_scale,
      seed = seed
    )
  }))
})
elapsed_seconds <- proc.time()[["elapsed"]] - start_time
seed_results <- do.call(rbind, seed_results)
rownames(seed_results) <- NULL

bf_results <- seed_results[seed_results$adjustment == "BF", , drop = FALSE]
paired <- merge(
  bf_results[bf_results$family == "FASH-linear", ],
  bf_results[bf_results$family == "FASH-IWP1", ],
  by = c(
    "candidate", "alternative_distribution", "endpoint_scale", "seed",
    "adjustment"
  ),
  suffixes = c("_linear", "_iwp"),
  sort = TRUE
)
paired$power_difference <- paired$power_linear - paired$power_iwp
paired$fdp_difference <-
  paired$realized_fdp_linear - paired$realized_fdp_iwp
paired$pi0_difference <-
  paired$estimated_pi0_linear - paired$estimated_pi0_iwp

confirmation_summary <- do.call(rbind, lapply(
  split(paired, paired$candidate),
  function(x) {
    linear_power_ci <- mean_ci(x$power_linear)
    iwp_power_ci <- mean_ci(x$power_iwp)
    difference_ci <- mean_ci(x$power_difference)
    data.frame(
      candidate = x$candidate[[1L]],
      alternative_distribution = x$alternative_distribution[[1L]],
      endpoint_scale = x$endpoint_scale[[1L]],
      n_seeds = nrow(x),
      mean_power_linear = linear_power_ci[["estimate"]],
      mean_power_linear_lower = linear_power_ci[["lower"]],
      mean_power_linear_upper = linear_power_ci[["upper"]],
      mean_power_iwp = iwp_power_ci[["estimate"]],
      mean_power_iwp_lower = iwp_power_ci[["lower"]],
      mean_power_iwp_upper = iwp_power_ci[["upper"]],
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
  match(candidate_settings$candidate, confirmation_summary$candidate),
  ,
  drop = FALSE
]

validation <- data.frame(
  check = c(
    "confirmation seeds are held out",
    "complete seed-method results",
    "paired BF summaries are complete",
    "finite confirmation estimates",
    "one-job thread cap"
  ),
  passed = c(
    !screen_cache$configuration$screen_seed %in% confirmation_seeds,
    nrow(seed_results) ==
      nrow(candidate_settings) * length(confirmation_seeds) * 4L,
    nrow(paired) == nrow(candidate_settings) * length(confirmation_seeds) &&
      all(paired$bf_available_linear) &&
      all(paired$bf_available_iwp),
    all(is.finite(as.matrix(confirmation_summary[
      vapply(confirmation_summary, is.numeric, logical(1))
    ]))),
    Sys.getenv("VECLIB_MAXIMUM_THREADS") == "1"
  ),
  stringsAsFactors = FALSE
)
if (any(!validation$passed)) {
  print(validation)
  stop("The boundary confirmation failed validation.")
}

output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r0_linear_iwp_power_mechanism_boundary_confirmation"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
analysis_cache <- list(
  configuration = list(
    scope = "internal held-out confirmation",
    J = J,
    pi0 = pi0,
    time_grid = time_grid,
    standard_error = standard_error,
    candidate_settings = candidate_settings,
    confirmation_seeds = confirmation_seeds,
    grid = grid,
    penalty = 10L,
    alpha = alpha,
    num_basis = num_basis,
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
