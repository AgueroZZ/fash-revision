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

J <- 10000L
pi0 <- 0.95
time_grid <- 0:31
standard_error <- rep(1, length(time_grid))
grid <- default_revision_grid()
slope_grid_value <- grid[which.min(abs(grid - 0.0742735782143339))]
if (!slope_grid_value %in% grid) {
  stop("The requested generative slope SD is not on the fitted grid.")
}
endpoint_sd <- slope_grid_value * diff(range(time_grid))
alpha <- 0.05
num_basis <- 20L
confirmation_seeds <- seq.int(960001L, by = 1000L, length.out = 20L)
effect_breaks <- c(0, 1, 2, 2.5, 3, 4, Inf)
effect_labels <- c("[0,1)", "[1,2)", "[2,2.5)", "[2.5,3)", "[3,4)", "[4,Inf)")

linear_kernel <- tcrossprod(time_grid - min(time_grid))
iwp_kernel <- compute_iwp1_kernel(time_grid, num_basis)

fit_bf <- function(likelihood) {
  raw <- fit_linear_mixture_fash_from_log_likelihood(
    L_matrix = likelihood,
    grid = grid,
    pred_step = 1,
    penalty = 10
  )
  BF_update_linear_mixture_fash(raw)
}

run_one <- function(seed) {
  simulated <- simulate_linear_power_mechanism(
    J = J,
    pi0 = pi0,
    time_grid = time_grid,
    endpoint_scale = endpoint_sd,
    standard_error = standard_error,
    alternative_distribution = "gaussian",
    equicorrelation = 0,
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
  linear_fit <- fit_bf(linear_likelihood)
  iwp_fit <- fit_bf(iwp_likelihood)
  linear_selected <- select_cumulative_lfdr_power_mechanism(
    linear_fit$lfdr,
    alpha
  )
  iwp_selected <- select_cumulative_lfdr_power_mechanism(
    iwp_fit$lfdr,
    alpha
  )
  true_alternative <- !simulated$true_null
  effect_bin <- cut(
    abs(simulated$endpoint),
    breaks = effect_breaks,
    labels = effect_labels,
    right = FALSE
  )

  bin_rows <- do.call(rbind, lapply(effect_labels, function(label) {
    in_bin <- true_alternative & effect_bin == label
    data.frame(
      seed = seed,
      effect_bin = label,
      n_alternatives = sum(in_bin),
      selected_linear = sum(linear_selected & in_bin),
      selected_iwp = sum(iwp_selected & in_bin),
      power_linear = sum(linear_selected & in_bin) / sum(in_bin),
      power_iwp = sum(iwp_selected & in_bin) / sum(in_bin),
      stringsAsFactors = FALSE
    )
  }))
  bin_rows$power_difference <- bin_rows$power_linear - bin_rows$power_iwp

  overall <- data.frame(
    seed = seed,
    power_linear = sum(linear_selected & true_alternative) /
      sum(true_alternative),
    power_iwp = sum(iwp_selected & true_alternative) /
      sum(true_alternative),
    power_difference = sum(linear_selected & true_alternative) /
      sum(true_alternative) -
      sum(iwp_selected & true_alternative) / sum(true_alternative),
    realized_fdp_linear = if (sum(linear_selected) == 0L) 0 else {
      sum(linear_selected & !true_alternative) / sum(linear_selected)
    },
    realized_fdp_iwp = if (sum(iwp_selected) == 0L) 0 else {
      sum(iwp_selected & !true_alternative) / sum(iwp_selected)
    },
    pi0_linear = constant_component_prior_weight(linear_fit),
    pi0_iwp = constant_component_prior_weight(iwp_fit),
    median_alternative_log_bf_difference = stats::median(
      log(linear_fit$BF[true_alternative]) -
        log(iwp_fit$BF[true_alternative])
    ),
    stringsAsFactors = FALSE
  )
  list(bin_rows = bin_rows, overall = overall)
}

start_time <- proc.time()[["elapsed"]]
seed_outputs <- lapply(confirmation_seeds, run_one)
elapsed_seconds <- proc.time()[["elapsed"]] - start_time
bin_seed_results <- do.call(rbind, lapply(seed_outputs, `[[`, "bin_rows"))
overall_seed_results <- do.call(rbind, lapply(seed_outputs, `[[`, "overall"))
rownames(bin_seed_results) <- NULL
rownames(overall_seed_results) <- NULL

bin_summary <- do.call(rbind, lapply(
  split(bin_seed_results, bin_seed_results$effect_bin),
  function(x) {
    linear_ci <- mean_ci(x$power_linear)
    iwp_ci <- mean_ci(x$power_iwp)
    difference_ci <- mean_ci(x$power_difference)
    data.frame(
      effect_bin = x$effect_bin[[1L]],
      total_alternatives = sum(x$n_alternatives),
      mean_power_linear = linear_ci[["estimate"]],
      mean_power_linear_lower = linear_ci[["lower"]],
      mean_power_linear_upper = linear_ci[["upper"]],
      mean_power_iwp = iwp_ci[["estimate"]],
      mean_power_iwp_lower = iwp_ci[["lower"]],
      mean_power_iwp_upper = iwp_ci[["upper"]],
      mean_power_difference = difference_ci[["estimate"]],
      power_difference_lower = difference_ci[["lower"]],
      power_difference_upper = difference_ci[["upper"]],
      mean_power_ratio = if (mean(x$power_iwp) == 0) NA_real_ else {
        mean(x$power_linear) / mean(x$power_iwp)
      },
      stringsAsFactors = FALSE
    )
  }
))
rownames(bin_summary) <- NULL
bin_summary <- bin_summary[
  match(effect_labels, bin_summary$effect_bin),
  ,
  drop = FALSE
]

overall_linear_ci <- mean_ci(overall_seed_results$power_linear)
overall_iwp_ci <- mean_ci(overall_seed_results$power_iwp)
overall_difference_ci <- mean_ci(overall_seed_results$power_difference)
overall_summary <- data.frame(
  n_seeds = length(confirmation_seeds),
  mean_power_linear = overall_linear_ci[["estimate"]],
  mean_power_linear_lower = overall_linear_ci[["lower"]],
  mean_power_linear_upper = overall_linear_ci[["upper"]],
  mean_power_iwp = overall_iwp_ci[["estimate"]],
  mean_power_iwp_lower = overall_iwp_ci[["lower"]],
  mean_power_iwp_upper = overall_iwp_ci[["upper"]],
  mean_power_difference = overall_difference_ci[["estimate"]],
  power_difference_lower = overall_difference_ci[["lower"]],
  power_difference_upper = overall_difference_ci[["upper"]],
  mean_power_ratio = mean(overall_seed_results$power_linear) /
    mean(overall_seed_results$power_iwp),
  empirical_fdr_linear = mean(overall_seed_results$realized_fdp_linear),
  empirical_fdr_iwp = mean(overall_seed_results$realized_fdp_iwp),
  mean_pi0_linear = mean(overall_seed_results$pi0_linear),
  mean_pi0_iwp = mean(overall_seed_results$pi0_iwp),
  mean_median_alternative_log_bf_difference = mean(
    overall_seed_results$median_alternative_log_bf_difference
  ),
  stringsAsFactors = FALSE
)

validation <- data.frame(
  check = c(
    "generative slope SD is on the fitted grid",
    "all errors are independent by construction",
    "supplied SE equals the generating SE",
    "all true alternatives are assigned to one bin",
    "every seed-bin has alternatives",
    "finite confirmation estimates",
    "one-job thread cap"
  ),
  passed = c(
    slope_grid_value %in% grid,
    TRUE,
    TRUE,
    sum(bin_seed_results$n_alternatives) ==
      length(confirmation_seeds) * as.integer(round(J * (1 - pi0))),
    all(bin_seed_results$n_alternatives > 0),
    all(is.finite(as.matrix(overall_summary))) &&
      all(is.finite(as.matrix(bin_summary[
        setdiff(names(bin_summary), c("effect_bin", "mean_power_ratio"))
      ]))),
    Sys.getenv("VECLIB_MAXIMUM_THREADS") == "1"
  ),
  stringsAsFactors = FALSE
)
if (any(!validation$passed)) {
  print(validation)
  stop("The clean Gaussian-grid effect-bin confirmation failed validation.")
}

output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r0_linear_iwp_power_mechanism_clean_gaussian_grid_bins"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
analysis_cache <- list(
  configuration = list(
    scope = "internal fully prior-matched Gaussian-slope confirmation",
    J = J,
    pi0 = pi0,
    time_grid = time_grid,
    standard_error = standard_error,
    slope_grid_value = slope_grid_value,
    endpoint_sd = endpoint_sd,
    confirmation_seeds = confirmation_seeds,
    effect_breaks = effect_breaks,
    error_model = "independent Gaussian with exact supplied SE",
    linear_truth_model = paste(
      "point null plus one zero-centered Gaussian slope component",
      "exactly on the fitted grid"
    ),
    grid = grid,
    penalty = 10L,
    alpha = alpha,
    num_basis = num_basis
  ),
  bin_seed_results = bin_seed_results,
  overall_seed_results = overall_seed_results,
  bin_summary = bin_summary,
  overall_summary = overall_summary,
  validation = validation,
  elapsed_seconds = elapsed_seconds
)
saveRDS(analysis_cache, file.path(output_directory, "analysis_cache.rds"))
utils::write.csv(
  bin_seed_results,
  file.path(output_directory, "bin_seed_results.csv"),
  row.names = FALSE
)
utils::write.csv(
  overall_seed_results,
  file.path(output_directory, "overall_seed_results.csv"),
  row.names = FALSE
)
utils::write.csv(
  bin_summary,
  file.path(output_directory, "bin_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  overall_summary,
  file.path(output_directory, "overall_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  validation,
  file.path(output_directory, "validation.csv"),
  row.names = FALSE
)

print(overall_summary)
print(bin_summary)
print(validation)
message("Elapsed seconds: ", round(elapsed_seconds, 3))
message("Saved output to: ", normalizePath(output_directory))
