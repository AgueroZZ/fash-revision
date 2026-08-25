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

J <- 2000L
pi0 <- 0.95
time_grid <- 0:15
standard_error <- 1
grid <- default_revision_grid()
alpha <- 0.05
num_basis <- 20L
screen_seed <- 20260821L
screen_grid <- expand.grid(
  equicorrelation = c(0, 0.2, 0.4, 0.6),
  endpoint_scale = seq(1.5, 3.25, by = 0.25),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

linear_kernel <- tcrossprod(time_grid - min(time_grid))
iwp_kernel <- compute_iwp1_kernel(time_grid, num_basis)

run_one <- function(equicorrelation, endpoint_scale) {
  simulated <- simulate_linear_power_mechanism(
    J = J,
    pi0 = pi0,
    time_grid = time_grid,
    endpoint_scale = endpoint_scale,
    standard_error = standard_error,
    alternative_distribution = "fixed",
    equicorrelation = equicorrelation,
    seed = screen_seed + as.integer(round(equicorrelation * 100))
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
  result$endpoint_scale <- endpoint_scale
  result
}

start_time <- proc.time()[["elapsed"]]
results <- do.call(rbind, lapply(seq_len(nrow(screen_grid)), function(index) {
  run_one(
    equicorrelation = screen_grid$equicorrelation[[index]],
    endpoint_scale = screen_grid$endpoint_scale[[index]]
  )
}))
elapsed_seconds <- proc.time()[["elapsed"]] - start_time
rownames(results) <- NULL

bf_results <- results[results$adjustment == "BF", , drop = FALSE]
paired <- merge(
  bf_results[bf_results$family == "FASH-linear", ],
  bf_results[bf_results$family == "FASH-IWP1", ],
  by = c("equicorrelation", "endpoint_scale", "adjustment"),
  suffixes = c("_linear", "_iwp"),
  sort = TRUE
)
paired$power_difference <- paired$power_linear - paired$power_iwp
paired$power_ratio <- ifelse(
  paired$power_iwp > 0,
  paired$power_linear / paired$power_iwp,
  NA_real_
)
paired$fdp_difference <-
  paired$realized_fdp_linear - paired$realized_fdp_iwp
paired <- paired[
  order(-paired$power_difference, paired$realized_fdp_linear),
  ,
  drop = FALSE
]
rownames(paired) <- NULL

validation <- data.frame(
  check = c(
    "screen is complete",
    "each error-correlation level has an available BF comparison",
    "all reported estimates are finite",
    "one-job thread cap"
  ),
  passed = c(
    nrow(results) == nrow(screen_grid) * 4L,
    all(vapply(
      split(bf_results$bf_available, bf_results$equicorrelation),
      any,
      logical(1)
    )),
    all(is.finite(results$power)) &&
      all(is.finite(results$realized_fdp)) &&
      all(is.finite(results$estimated_pi0)),
    Sys.getenv("VECLIB_MAXIMUM_THREADS") == "1"
  ),
  stringsAsFactors = FALSE
)
if (any(!validation$passed)) {
  print(validation)
  stop("The error-correlation screen failed validation.")
}

output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r0_linear_iwp_power_mechanism_error_correlation_screen"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
analysis_cache <- list(
  configuration = list(
    scope = "internal equicorrelated-error screen",
    J = J,
    pi0 = pi0,
    time_grid = time_grid,
    standard_error = standard_error,
    screen_seed = screen_seed,
    screen_grid = screen_grid,
    working_likelihood = "diagonal standard errors",
    truth_error_model = "equicorrelated Gaussian errors",
    grid = grid,
    penalty = 10L,
    alpha = alpha
  ),
  results = results,
  paired = paired,
  validation = validation,
  elapsed_seconds = elapsed_seconds
)
saveRDS(analysis_cache, file.path(output_directory, "analysis_cache.rds"))
utils::write.csv(
  results,
  file.path(output_directory, "screen_results.csv"),
  row.names = FALSE
)
utils::write.csv(
  paired,
  file.path(output_directory, "paired_bf_results.csv"),
  row.names = FALSE
)
utils::write.csv(
  validation,
  file.path(output_directory, "validation.csv"),
  row.names = FALSE
)

print(utils::head(
  paired[, c(
    "equicorrelation", "endpoint_scale",
    "power_linear", "power_iwp", "power_difference", "power_ratio",
    "realized_fdp_linear", "realized_fdp_iwp",
    "estimated_pi0_linear", "estimated_pi0_iwp"
  )],
  15L
))
print(validation)
message("Elapsed seconds: ", round(elapsed_seconds, 3))
message("Saved output to: ", normalizePath(output_directory))
