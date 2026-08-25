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

J <- 3000L
pi0 <- 0.95
grid <- default_revision_grid()
alpha <- 0.05
num_basis <- 20L
screen_seed <- 20260822L
target_endpoint_scales <- seq(1.5, 4, by = 0.25)

configuration_grid <- rbind(
  data.frame(
    configuration = paste0("T", c(16L, 32L, 64L), "_homoskedastic"),
    n_time = c(16L, 32L, 64L),
    se_pattern = "homoskedastic",
    stringsAsFactors = FALSE
  ),
  data.frame(
    configuration = paste0(
      "T16_",
      c("endpoint_mild", "endpoint_strong", "edge_quartet", "middle_precise")
    ),
    n_time = 16L,
    se_pattern = c(
      "endpoint_mild", "endpoint_strong", "edge_quartet", "middle_precise"
    ),
    stringsAsFactors = FALSE
  )
)
screen_grid <- merge(
  configuration_grid,
  expand.grid(
    alternative_distribution = c("fixed", "gaussian"),
    target_endpoint_scale = target_endpoint_scales,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ),
  by = NULL
)

run_one <- function(configuration,
                    n_time,
                    se_pattern,
                    alternative_distribution,
                    target_endpoint_scale) {
  time_grid <- seq.int(0, n_time - 1L)
  standard_error <- make_clean_standard_error(n_time, se_pattern)
  endpoint_scale <- target_endpoint_scale
  slope_grid_value <- NA_real_
  if (alternative_distribution == "gaussian") {
    positive_grid <- grid[grid > 0]
    slope_grid_value <- positive_grid[which.min(abs(
      positive_grid - target_endpoint_scale / (n_time - 1L)
    ))]
    endpoint_scale <- slope_grid_value * (n_time - 1L)
  }
  simulated <- simulate_linear_power_mechanism(
    J = J,
    pi0 = pi0,
    time_grid = time_grid,
    endpoint_scale = endpoint_scale,
    standard_error = standard_error,
    alternative_distribution = alternative_distribution,
    equicorrelation = 0,
    seed = screen_seed +
      match(configuration, configuration_grid$configuration) * 100L +
      match(alternative_distribution, c("fixed", "gaussian"))
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
  result$configuration <- configuration
  result$n_time <- n_time
  result$se_pattern <- se_pattern
  result$se_min <- min(standard_error)
  result$se_max <- max(standard_error)
  result$alternative_distribution <- alternative_distribution
  result$target_endpoint_scale <- target_endpoint_scale
  result$actual_endpoint_scale <- endpoint_scale
  result$slope_grid_value <- slope_grid_value
  result
}

start_time <- proc.time()[["elapsed"]]
screen_results <- do.call(rbind, lapply(seq_len(nrow(screen_grid)), function(index) {
  current <- screen_grid[index, ]
  run_one(
    configuration = current$configuration,
    n_time = current$n_time,
    se_pattern = current$se_pattern,
    alternative_distribution = current$alternative_distribution,
    target_endpoint_scale = current$target_endpoint_scale
  )
}))
elapsed_seconds <- proc.time()[["elapsed"]] - start_time
rownames(screen_results) <- NULL

bf_results <- screen_results[screen_results$adjustment == "BF", , drop = FALSE]
pairing_columns <- c(
  "configuration", "n_time", "se_pattern", "se_min", "se_max",
  "alternative_distribution", "target_endpoint_scale",
  "actual_endpoint_scale", "slope_grid_value", "adjustment"
)
paired <- merge(
  bf_results[bf_results$family == "FASH-linear", ],
  bf_results[bf_results$family == "FASH-IWP1", ],
  by = pairing_columns,
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
paired$pi0_difference <-
  paired$estimated_pi0_linear - paired$estimated_pi0_iwp
paired$median_alternative_log_bf_difference <-
  paired$median_log_bf_alternative_linear -
    paired$median_log_bf_alternative_iwp
paired <- paired[
  order(-paired$power_difference, paired$realized_fdp_linear),
  ,
  drop = FALSE
]
rownames(paired) <- NULL

validation <- data.frame(
  check = c(
    "all errors are independent by construction",
    "supplied SE equals the generating SE",
    "Gaussian slopes use an exact fitted grid component",
    "screen is complete",
    "every configuration has an available BF comparison",
    "all reported estimates are finite",
    "one-job thread cap"
  ),
  passed = c(
    TRUE,
    TRUE,
    all(is.na(screen_results$slope_grid_value[
      screen_results$alternative_distribution == "fixed"
    ])) &&
      all(screen_results$slope_grid_value[
        screen_results$alternative_distribution == "gaussian"
      ] %in% grid),
    nrow(screen_results) == nrow(screen_grid) * 4L,
    all(vapply(
      split(bf_results$bf_available, bf_results$configuration),
      any,
      logical(1)
    )),
    all(is.finite(screen_results$power)) &&
      all(is.finite(screen_results$realized_fdp)) &&
      all(is.finite(screen_results$estimated_pi0)),
    Sys.getenv("VECLIB_MAXIMUM_THREADS") == "1"
  ),
  stringsAsFactors = FALSE
)
if (any(!validation$passed)) {
  print(validation)
  stop("The clean correctly specified screen failed validation.")
}

output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r0_linear_iwp_power_mechanism_clean_screen"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
analysis_cache <- list(
  configuration = list(
    scope = "internal correctly specified independent-error screen",
    J = J,
    pi0 = pi0,
    screen_seed = screen_seed,
    configuration_grid = configuration_grid,
    target_endpoint_scales = target_endpoint_scales,
    error_model = "independent Gaussian with exact supplied SE",
    grid = grid,
    penalty = 10L,
    alpha = alpha,
    num_basis = num_basis
  ),
  screen_results = screen_results,
  paired = paired,
  validation = validation,
  elapsed_seconds = elapsed_seconds
)
saveRDS(analysis_cache, file.path(output_directory, "analysis_cache.rds"))
utils::write.csv(
  screen_results,
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
    "configuration", "alternative_distribution", "actual_endpoint_scale",
    "power_linear", "power_iwp", "power_difference", "power_ratio",
    "realized_fdp_linear", "realized_fdp_iwp",
    "estimated_pi0_linear", "estimated_pi0_iwp"
  )],
  20L
))
print(validation)
message("Elapsed seconds: ", round(elapsed_seconds, 3))
message("Saved output to: ", normalizePath(output_directory))
