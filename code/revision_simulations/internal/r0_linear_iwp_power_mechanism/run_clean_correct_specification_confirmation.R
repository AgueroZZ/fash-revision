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
grid <- default_revision_grid()
alpha <- 0.05
num_basis <- 20L
confirmation_seeds <- seq.int(930001L, by = 1000L, length.out = 20L)
candidate_settings <- data.frame(
  candidate = c(
    "T16_homoskedastic_fixed",
    "T16_middle_precise_boundary_fixed",
    "T16_middle_precise_fixed",
    "T32_homoskedastic_fixed",
    "T32_homoskedastic_gaussian_grid"
  ),
  n_time = c(16L, 16L, 16L, 32L, 32L),
  se_pattern = c(
    "homoskedastic",
    "middle_precise",
    "middle_precise",
    "homoskedastic",
    "homoskedastic"
  ),
  alternative_distribution = c(
    "fixed", "fixed", "fixed", "fixed", "gaussian"
  ),
  target_endpoint_scale = c(3.00, 2.50, 3.25, 2.25, 2.25),
  stringsAsFactors = FALSE
)

screen_cache_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r0_linear_iwp_power_mechanism_clean_screen",
  "analysis_cache.rds"
)
screen_cache <- readRDS(screen_cache_path)
if (screen_cache$configuration$screen_seed %in% confirmation_seeds) {
  stop("The confirmation seeds overlap the screen seed.")
}

run_one <- function(setting, seed) {
  time_grid <- seq.int(0, setting$n_time - 1L)
  standard_error <- make_clean_standard_error(
    setting$n_time,
    setting$se_pattern
  )
  endpoint_scale <- setting$target_endpoint_scale
  slope_grid_value <- NA_real_
  if (setting$alternative_distribution == "gaussian") {
    positive_grid <- grid[grid > 0]
    slope_grid_value <- positive_grid[which.min(abs(
      positive_grid - endpoint_scale / (setting$n_time - 1L)
    ))]
    endpoint_scale <- slope_grid_value * (setting$n_time - 1L)
  }
  simulated <- simulate_linear_power_mechanism(
    J = J,
    pi0 = pi0,
    time_grid = time_grid,
    endpoint_scale = endpoint_scale,
    standard_error = standard_error,
    alternative_distribution = setting$alternative_distribution,
    equicorrelation = 0,
    seed = seed
  )
  linear_likelihood <- compute_common_intercept_log_likelihood(
    simulated$y_matrix,
    standard_error,
    tcrossprod(time_grid - min(time_grid)),
    grid
  )
  iwp_likelihood <- compute_common_intercept_log_likelihood(
    simulated$y_matrix,
    standard_error,
    compute_iwp1_kernel(time_grid, num_basis),
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
  result$candidate <- setting$candidate
  result$n_time <- setting$n_time
  result$se_pattern <- setting$se_pattern
  result$se_min <- min(standard_error)
  result$se_max <- max(standard_error)
  result$alternative_distribution <- setting$alternative_distribution
  result$target_endpoint_scale <- setting$target_endpoint_scale
  result$actual_endpoint_scale <- endpoint_scale
  result$slope_grid_value <- slope_grid_value
  result$seed <- seed
  result
}

start_time <- proc.time()[["elapsed"]]
seed_results <- do.call(rbind, lapply(seq_len(nrow(candidate_settings)), function(index) {
  setting <- candidate_settings[index, ]
  do.call(rbind, lapply(confirmation_seeds, function(seed) {
    run_one(setting, seed)
  }))
}))
elapsed_seconds <- proc.time()[["elapsed"]] - start_time
rownames(seed_results) <- NULL

bf_results <- seed_results[seed_results$adjustment == "BF", , drop = FALSE]
pairing_columns <- c(
  "candidate", "n_time", "se_pattern", "se_min", "se_max",
  "alternative_distribution", "target_endpoint_scale",
  "actual_endpoint_scale", "slope_grid_value", "seed", "adjustment"
)
paired <- merge(
  bf_results[bf_results$family == "FASH-linear", ],
  bf_results[bf_results$family == "FASH-IWP1", ],
  by = pairing_columns,
  suffixes = c("_linear", "_iwp"),
  sort = TRUE
)
paired$power_difference <- paired$power_linear - paired$power_iwp
paired$fdp_difference <-
  paired$realized_fdp_linear - paired$realized_fdp_iwp
paired$pi0_difference <-
  paired$estimated_pi0_linear - paired$estimated_pi0_iwp
paired$median_alternative_log_bf_difference <-
  paired$median_log_bf_alternative_linear -
    paired$median_log_bf_alternative_iwp

confirmation_summary <- do.call(rbind, lapply(
  split(paired, paired$candidate),
  function(x) {
    linear_ci <- mean_ci(x$power_linear)
    iwp_ci <- mean_ci(x$power_iwp)
    difference_ci <- mean_ci(x$power_difference)
    data.frame(
      candidate = x$candidate[[1L]],
      n_time = x$n_time[[1L]],
      se_pattern = x$se_pattern[[1L]],
      se_min = x$se_min[[1L]],
      se_max = x$se_max[[1L]],
      alternative_distribution = x$alternative_distribution[[1L]],
      target_endpoint_scale = x$target_endpoint_scale[[1L]],
      actual_endpoint_scale = x$actual_endpoint_scale[[1L]],
      slope_grid_value = x$slope_grid_value[[1L]],
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
      mean_alternative_log_bf_difference = mean(
        x$median_alternative_log_bf_difference
      ),
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

screen_candidates <- merge(
  candidate_settings,
  screen_cache$paired,
  by.x = c(
    "n_time", "se_pattern", "alternative_distribution",
    "target_endpoint_scale"
  ),
  by.y = c(
    "n_time", "se_pattern", "alternative_distribution",
    "target_endpoint_scale"
  ),
  all.x = TRUE,
  sort = FALSE
)
validation <- data.frame(
  check = c(
    "confirmation seeds are held out",
    "candidate settings were selected from the screen",
    "all errors are independent by construction",
    "supplied SE equals the generating SE",
    "Gaussian slope uses an exact fitted grid component",
    "complete paired BF results",
    "finite confirmation estimates",
    "one-job thread cap"
  ),
  passed = c(
    !screen_cache$configuration$screen_seed %in% confirmation_seeds,
    nrow(screen_candidates) >= nrow(candidate_settings) &&
      all(candidate_settings$candidate %in% screen_candidates$candidate),
    TRUE,
    TRUE,
    all(confirmation_summary$slope_grid_value[
      confirmation_summary$alternative_distribution == "gaussian"
    ] %in% grid),
    nrow(paired) == nrow(candidate_settings) * length(confirmation_seeds) &&
      all(paired$bf_available_linear) &&
      all(paired$bf_available_iwp),
    all(is.finite(as.matrix(confirmation_summary[
      vapply(confirmation_summary, is.numeric, logical(1)) &
        names(confirmation_summary) != "slope_grid_value"
    ]))),
    Sys.getenv("VECLIB_MAXIMUM_THREADS") == "1"
  ),
  stringsAsFactors = FALSE
)
if (any(!validation$passed)) {
  print(validation)
  stop("The clean correctly specified confirmation failed validation.")
}

output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r0_linear_iwp_power_mechanism_clean_confirmation"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
analysis_cache <- list(
  configuration = list(
    scope = "internal held-out correctly specified confirmation",
    J = J,
    pi0 = pi0,
    candidate_settings = candidate_settings,
    confirmation_seeds = confirmation_seeds,
    error_model = "independent Gaussian with exact supplied SE",
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
