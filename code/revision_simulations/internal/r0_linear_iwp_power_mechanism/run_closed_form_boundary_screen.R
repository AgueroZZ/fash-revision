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

compute_iwp1_kernel <- function(time_grid, num_basis = 20L) {
  template <- data.frame(
    x = time_grid,
    y = rep(0, length(time_grid)),
    offset = rep(0, length(time_grid))
  )
  tmb_data <- fashr:::fash_set_tmbdat(
    template,
    Si = rep(1, length(time_grid)),
    num_basis = num_basis,
    order = 1
  )
  basis <- as.matrix(tmb_data$B)
  precision <- as.matrix(tmb_data$P)
  basis %*% solve(precision, t(basis))
}

compute_common_intercept_log_likelihood <- function(y_matrix,
                                                    standard_error,
                                                    kernel,
                                                    grid) {
  y_matrix <- as.matrix(y_matrix)
  n_time <- ncol(y_matrix)
  if (!is.matrix(kernel) ||
      any(dim(kernel) != n_time) ||
      length(standard_error) != 1L ||
      !is.finite(standard_error) ||
      standard_error <= 0) {
    stop("The likelihood inputs are invalid.")
  }
  intercept <- rep(1, n_time)
  likelihood <- vapply(grid, function(scale) {
    covariance <- diag(standard_error^2, n_time) + scale^2 * kernel
    chol_covariance <- chol(covariance)
    inverse_covariance <- chol2inv(chol_covariance)
    inverse_intercept <- as.numeric(inverse_covariance %*% intercept)
    intercept_precision <- sum(inverse_intercept)
    residual_precision <- inverse_covariance -
      tcrossprod(inverse_intercept) / intercept_precision
    quadratic <- rowSums((y_matrix %*% residual_precision) * y_matrix)
    log_determinant <- 2 * sum(log(diag(chol_covariance)))
    -0.5 * (
      (n_time - 1) * log(2 * pi) +
        log_determinant +
        log(intercept_precision) +
        quadratic
    )
  }, numeric(nrow(y_matrix)))
  dimnames(likelihood) <- list(rownames(y_matrix), as.character(grid))
  likelihood
}

select_cumulative_lfdr <- function(lfdr, alpha = 0.05) {
  ordering <- order(lfdr, seq_along(lfdr))
  cumulative_fdr <- cumsum(lfdr[ordering]) / seq_along(ordering)
  last_selected <- max(c(0L, which(cumulative_fdr <= alpha)))
  selected <- rep(FALSE, length(lfdr))
  if (last_selected > 0L) {
    selected[ordering[seq_len(last_selected)]] <- TRUE
  }
  selected
}

summarize_fit <- function(fit, true_null, family, adjustment, alpha = 0.05) {
  selected <- select_cumulative_lfdr(fit$lfdr, alpha)
  discoveries <- sum(selected)
  false_discoveries <- sum(selected & true_null)
  true_positives <- sum(selected & !true_null)
  bf_available <- isTRUE(fit$bf_adjusted) &&
    !is.null(fit$BF) &&
    all(is.finite(fit$BF)) &&
    all(fit$BF > 0)
  data.frame(
    family = family,
    adjustment = adjustment,
    n_discoveries = discoveries,
    false_discoveries = false_discoveries,
    true_positives = true_positives,
    power = true_positives / sum(!true_null),
    realized_fdp = if (discoveries == 0L) 0 else {
      false_discoveries / discoveries
    },
    estimated_pi0 = constant_component_prior_weight(fit),
    bf_available = bf_available,
    median_log_bf_alternative = if (bf_available) {
      stats::median(log(fit$BF[!true_null]))
    } else {
      NA_real_
    },
    median_log_bf_null = if (bf_available) {
      stats::median(log(fit$BF[true_null]))
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

simulate_fixed_linear_matrix <- function(J,
                                         pi0,
                                         time_grid,
                                         endpoint_change,
                                         standard_error,
                                         seed) {
  set.seed(seed)
  n_null <- as.integer(round(J * pi0))
  true_null <- c(rep(TRUE, n_null), rep(FALSE, J - n_null))
  intercept <- stats::rnorm(J)
  direction <- sample(c(-1, 1), J - n_null, replace = TRUE)
  endpoint <- numeric(J)
  endpoint[!true_null] <- endpoint_change * direction
  mean_matrix <- intercept + outer(
    endpoint,
    (time_grid - min(time_grid)) / diff(range(time_grid))
  )
  noise <- matrix(
    stats::rnorm(J * length(time_grid), sd = standard_error),
    nrow = J,
    ncol = length(time_grid)
  )
  y_matrix <- mean_matrix + noise
  rownames(y_matrix) <- sprintf("unit_%05d", seq_len(J))
  list(y_matrix = y_matrix, true_null = true_null)
}

fit_family <- function(likelihood, grid, true_null, family, alpha) {
  raw <- fit_linear_mixture_fash_from_log_likelihood(
    L_matrix = likelihood,
    grid = grid,
    pred_step = 1,
    penalty = 10
  )
  bf <- tryCatch(
    suppressWarnings(BF_update_linear_mixture_fash(raw)),
    error = function(condition) NULL
  )
  if (is.null(bf) || is.null(bf$BF)) {
    bf <- list(
      lfdr = rep(1, nrow(likelihood)),
      prior_weights = data.frame(psd = 0, prior_weight = 1),
      BF = rep(NA_real_, nrow(likelihood)),
      bf_adjusted = TRUE
    )
  }
  rbind(
    summarize_fit(raw, true_null, family, "Raw", alpha),
    summarize_fit(bf, true_null, family, "BF", alpha)
  )
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

J <- 2000L
screen_seed <- 20260820L
time_counts <- c(8L, 16L, 32L)
pi0_values <- c(0.80, 0.95)
endpoint_changes <- seq(0.75, 4, by = 0.25)
standard_error <- 1
grid <- default_revision_grid()
alpha <- 0.05
num_basis <- 20L

validation_time <- 0:15
validation_data <- simulate_fixed_linear_matrix(
  J = 6L,
  pi0 = 0.5,
  time_grid = validation_time,
  endpoint_change = 2,
  standard_error = standard_error,
  seed = 901L
)
validation_datasets <- lapply(seq_len(nrow(validation_data$y_matrix)), function(i) {
  data.frame(
    x = validation_time,
    y = validation_data$y_matrix[i, ],
    sd = rep(standard_error, length(validation_time))
  )
})
names(validation_datasets) <- rownames(validation_data$y_matrix)

validation_iwp_public <- fashr::fash(
  Y = "y",
  smooth_var = "x",
  S = "sd",
  data_list = validation_datasets,
  order = 1,
  grid = grid,
  num_basis = num_basis,
  pred_step = 1,
  penalty = 10,
  num_cores = 1,
  verbose = FALSE
)
validation_iwp_closed <- compute_common_intercept_log_likelihood(
  validation_data$y_matrix,
  standard_error,
  compute_iwp1_kernel(validation_time, num_basis),
  grid
)
validation_linear_public <- compute_linear_mixture_log_likelihood(
  validation_datasets,
  grid = grid,
  pred_step = 1
)
validation_linear_closed <- compute_common_intercept_log_likelihood(
  validation_data$y_matrix,
  standard_error,
  tcrossprod(validation_time - min(validation_time)),
  grid
)
row_center <- function(x) x - x[, 1L]
max_iwp_likelihood_difference <- max(abs(
  row_center(validation_iwp_public$L_matrix) -
    row_center(validation_iwp_closed)
))
max_linear_likelihood_difference <- max(abs(
  row_center(validation_linear_public) -
    row_center(validation_linear_closed)
))

screen_grid <- expand.grid(
  n_time = time_counts,
  pi0 = pi0_values,
  endpoint_change = endpoint_changes,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)

screen_one <- function(n_time, pi0, endpoint_change) {
  time_grid <- seq.int(0, n_time - 1L)
  simulated <- simulate_fixed_linear_matrix(
    J = J,
    pi0 = pi0,
    time_grid = time_grid,
    endpoint_change = endpoint_change,
    standard_error = standard_error,
    seed = screen_seed + n_time * 100L + as.integer(round(pi0 * 100))
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
    fit_family(
      linear_likelihood,
      grid,
      simulated$true_null,
      "FASH-linear",
      alpha
    ),
    fit_family(
      iwp_likelihood,
      grid,
      simulated$true_null,
      "FASH-IWP1",
      alpha
    )
  )
  result$n_time <- n_time
  result$pi0 <- pi0
  result$endpoint_change <- endpoint_change
  result
}

start_time <- proc.time()[["elapsed"]]
screen_results <- lapply(seq_len(nrow(screen_grid)), function(index) {
  current <- screen_grid[index, ]
  screen_one(
    n_time = current$n_time,
    pi0 = current$pi0,
    endpoint_change = current$endpoint_change
  )
})
elapsed_seconds <- proc.time()[["elapsed"]] - start_time
screen_results <- do.call(rbind, screen_results)
rownames(screen_results) <- NULL

bf_results <- screen_results[screen_results$adjustment == "BF", , drop = FALSE]
paired_bf <- merge(
  bf_results[bf_results$family == "FASH-linear", ],
  bf_results[bf_results$family == "FASH-IWP1", ],
  by = c("n_time", "pi0", "endpoint_change", "adjustment"),
  suffixes = c("_linear", "_iwp"),
  sort = TRUE
)
paired_bf$power_difference <- paired_bf$power_linear - paired_bf$power_iwp
paired_bf$fdp_difference <-
  paired_bf$realized_fdp_linear - paired_bf$realized_fdp_iwp
paired_bf$pi0_difference <-
  paired_bf$estimated_pi0_linear - paired_bf$estimated_pi0_iwp
paired_bf$median_alternative_log_bf_difference <-
  paired_bf$median_log_bf_alternative_linear -
    paired_bf$median_log_bf_alternative_iwp
paired_bf <- paired_bf[
  order(-paired_bf$power_difference, paired_bf$realized_fdp_linear),
  ,
  drop = FALSE
]
rownames(paired_bf) <- NULL

validation <- data.frame(
  check = c(
    "closed-form IWP1 likelihood matches public FASH",
    "closed-form linear likelihood matches shared implementation",
    "screen is complete and finite",
    "one-job thread cap"
  ),
  passed = c(
    max_iwp_likelihood_difference <= 1e-5,
    max_linear_likelihood_difference <= 1e-8,
    nrow(screen_results) == nrow(screen_grid) * 4L &&
      all(is.finite(screen_results$power)) &&
      all(is.finite(screen_results$realized_fdp)) &&
      all(is.finite(screen_results$estimated_pi0)),
    Sys.getenv("VECLIB_MAXIMUM_THREADS") == "1"
  ),
  observed = c(
    format(max_iwp_likelihood_difference, scientific = TRUE),
    format(max_linear_likelihood_difference, scientific = TRUE),
    paste(nrow(screen_results), "method-setting rows"),
    "one R job; one BLAS thread"
  ),
  stringsAsFactors = FALSE
)
if (any(!validation$passed)) {
  print(validation)
  stop("The closed-form boundary screen failed validation.")
}

output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r0_linear_iwp_power_mechanism_boundary_screen"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
analysis_cache <- list(
  configuration = list(
    scope = "internal held-in boundary screen",
    J = J,
    screen_seed = screen_seed,
    time_counts = time_counts,
    pi0_values = pi0_values,
    endpoint_changes = endpoint_changes,
    standard_error = standard_error,
    alternative_distribution = "fixed magnitude with random sign",
    grid = grid,
    penalty = 10L,
    alpha = alpha,
    num_basis = num_basis
  ),
  screen_results = screen_results,
  paired_bf = paired_bf,
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
  paired_bf,
  file.path(output_directory, "paired_bf_results.csv"),
  row.names = FALSE
)
utils::write.csv(
  validation,
  file.path(output_directory, "validation.csv"),
  row.names = FALSE
)

print(utils::head(
  paired_bf[, c(
    "n_time", "pi0", "endpoint_change",
    "power_linear", "power_iwp", "power_difference",
    "realized_fdp_linear", "realized_fdp_iwp",
    "estimated_pi0_linear", "estimated_pi0_iwp",
    "median_alternative_log_bf_difference"
  )],
  15L
))
print(validation)
message("Elapsed seconds: ", round(elapsed_seconds, 3))
message("Saved output to: ", normalizePath(output_directory))
