#!/usr/bin/env Rscript

# Validate E0(BF) = 1 using an ordinary null BF sample mean in R2 and in a
# prespecified finite-variance control when the fitted R2 alternative remains
# too broad for practical ordinary Monte Carlo convergence.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}

log_mvn_density <- function(y, covariance) {
  cholesky <- chol(covariance)
  standardized <- forwardsolve(t(cholesky), y)
  -0.5 * (
    length(y) * log(2 * pi) +
      2 * sum(log(diag(cholesky))) +
      sum(standardized^2)
  )
}

gaussian_bf_second_moment <- function(null_covariance,
                                      alternative_covariance) {
  precision <-
    2 * solve(alternative_covariance) - solve(null_covariance)
  precision <- (precision + t(precision)) / 2
  precision_eigenvalues <- eigen(
    precision,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  if (min(precision_eigenvalues) <= 0) {
    return(Inf)
  }
  log_determinant <- function(matrix) {
    determinant(matrix, logarithm = TRUE)$modulus[1]
  }
  exp(
    -log_determinant(alternative_covariance) +
      0.5 * log_determinant(null_covariance) -
      0.5 * log_determinant(precision)
  )
}

summarize_bf <- function(bf, scenario, theoretical_second_moment) {
  quantiles <- stats::quantile(
    bf,
    probs = c(0.50, 0.90, 0.95, 0.99, 0.999, 0.9999),
    type = 8,
    names = FALSE
  )
  data.frame(
    scenario = scenario,
    n = length(bf),
    sample_mean = mean(bf),
    sample_median = quantiles[1L],
    q90 = quantiles[2L],
    q95 = quantiles[3L],
    q99 = quantiles[4L],
    q999 = quantiles[5L],
    q9999 = quantiles[6L],
    maximum = max(bf),
    proportion_greater_than_one = mean(bf > 1),
    theoretical_mean = 1,
    theoretical_second_moment = theoretical_second_moment,
    theoretical_variance = theoretical_second_moment - 1,
    stringsAsFactors = FALSE
  )
}

simulate_exact_null_bf <- function(null_covariance,
                                   alternative_covariances,
                                   alternative_weights,
                                   scenario,
                                   seed,
                                   n_simulations,
                                   batch_size,
                                   checkpoints) {
  null_cholesky <- chol(null_covariance)
  null_log_determinant <- 2 * sum(log(diag(null_cholesky)))
  alternative_cholesky <- lapply(alternative_covariances, chol)
  alternative_log_determinants <- vapply(
    alternative_cholesky,
    function(cholesky) 2 * sum(log(diag(cholesky))),
    numeric(1)
  )
  alternative_cholesky_inverse <- lapply(
    alternative_cholesky, solve
  )

  set.seed(seed)
  bayes_factor <- numeric(n_simulations)
  position <- 0L
  while (position < n_simulations) {
    current_n <- min(batch_size, n_simulations - position)
    z <- matrix(
      stats::rnorm(current_n * nrow(null_covariance)),
      nrow = current_n,
      ncol = nrow(null_covariance)
    )
    y <- z %*% null_cholesky
    null_quadratic <- rowSums(z^2)
    log_component_bf <- vapply(
      seq_along(alternative_covariances),
      function(i) {
        alternative_standardized <-
          y %*% alternative_cholesky_inverse[[i]]
        alternative_quadratic <- rowSums(
          alternative_standardized^2
        )
        -0.5 * (
          alternative_log_determinants[i] - null_log_determinant +
            alternative_quadratic - null_quadratic
        )
      },
      numeric(current_n)
    )
    if (is.null(dim(log_component_bf))) {
      log_component_bf <- matrix(log_component_bf, ncol = 1L)
    }
    weighted_log_component_bf <- sweep(
      log_component_bf,
      2L,
      log(alternative_weights),
      `+`
    )
    row_maximum <- apply(weighted_log_component_bf, 1L, max)
    batch_bf <- exp(row_maximum) * rowSums(exp(
      weighted_log_component_bf - row_maximum
    ))
    indices <- position + seq_len(current_n)
    bayes_factor[indices] <- batch_bf
    position <- position + current_n
  }

  cumulative_sum <- cumsum(bayes_factor)
  running <- data.frame(
    scenario = scenario,
    n = checkpoints,
    running_mean_bf = cumulative_sum[checkpoints] / checkpoints,
    stringsAsFactors = FALSE
  )
  list(bayes_factor = bayes_factor, running = running)
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation",
  "selected_signal_genotype_permutation_helpers.R"
))

if (!requireNamespace("fashr", quietly = TRUE) ||
    !requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("scales", quietly = TRUE)) {
  stop("The fashr, ggplot2, and scales packages are required.")
}

source_path <- file.path(
  workflowr_root, "output", "revision_simulations", "mc",
  "r2_timed_cosine_one_two_three_peak_main_effect_profile_sigma_pilot5",
  "full_fits", "seed_12345.rds"
)
output_id <- "r2_exact_null_bf_mild_setting_validation_v1"
output_directory <- file.path(
  workflowr_root, "output", "revision_simulations", "internal", output_id
)
figure_directory <- file.path(output_directory, "figures")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

source_object <- readRDS(source_path)
source_fit <- source_object$fash_fits$fash_iwp1_raw
effect_class <- as.character(source_object$unit_info$effect_class)
alternative_indices <- which(effect_class == "dynamic_bspline")[
  seq_len(200L)
]
null_indices <- which(effect_class == "zero")[seq_len(200L)]
selected_indices <- c(alternative_indices, null_indices)
unit_keys <- as.character(source_object$unit_info$unit_id[selected_indices])

raw_fit <- refit_fash_from_likelihood(
  source_fit = source_fit,
  data_list = source_fit$fash_data$data_list[selected_indices],
  se_list = source_fit$fash_data$S[selected_indices],
  likelihood_matrix = source_fit$L_matrix[selected_indices, , drop = FALSE],
  unit_keys = unit_keys,
  penalty = 10
)
full_weights <- numeric(length(source_fit$psd_grid))
full_weights[match(
  raw_fit$prior_weights$psd, source_fit$psd_grid
)] <- raw_fit$prior_weights$prior_weight
active_alternative <- which(full_weights[-1L] > 1e-10) + 1L
alternative_weights <- full_weights[active_alternative] /
  sum(full_weights[active_alternative])
alternative_psd <- source_fit$psd_grid[active_alternative]
if (length(active_alternative) < 1L) {
  stop("The R2 reference fit has no active alternative components.")
}

oracle_se <- vapply(null_indices, function(unit_index) {
  design <- cbind(
    intercept = 1,
    G = source_object$genotype[, unit_index],
    source_object$covariates
  )
  sqrt(
    source_object$settings$expression_noise_sd^2 *
      solve(crossprod(design))["G", "G"]
  )
}, numeric(1))
representative_se <- stats::median(oracle_se)
representative_data <- source_fit$fash_data$data_list[[null_indices[1L]]]
representative_s <- rep(representative_se, nrow(representative_data))
tmb_data <- fashr:::fash_set_tmbdat(
  representative_data,
  Si = representative_s,
  Omegai = NULL,
  num_basis = source_fit$settings$num_basis,
  betaprec = source_fit$settings$betaprec,
  order = source_fit$settings$order
)
fixed_design <- as.matrix(tmb_data$X)
random_design <- as.matrix(tmb_data$B)
random_precision <- as.matrix(tmb_data$P)

null_covariance <-
  diag(representative_s^2) +
  fixed_design %*% t(fixed_design) / source_fit$settings$betaprec
alternative_covariances <- lapply(alternative_psd, function(psd) {
  sigma_iwp <- psd / sqrt(
    source_fit$settings$pred_step ^
      ((2 * source_fit$settings$order) - 1) /
      (
        ((2 * source_fit$settings$order) - 1) *
          factorial(source_fit$settings$order - 1)^2
      )
  )
  null_covariance +
    random_design %*%
      (solve(random_precision) * sigma_iwp^2) %*%
      t(random_design)
})

analytic_log_likelihood <- c(
  log_mvn_density(representative_data$y, null_covariance),
  vapply(alternative_covariances, function(covariance) {
    log_mvn_density(representative_data$y, covariance)
  }, numeric(1))
)
tmb_log_likelihood <- fashr:::compute_L_gaussian_helper_seq(
  data_i = representative_data,
  Si = representative_s,
  Omegai = NULL,
  grid = c(0, alternative_psd),
  num_basis = source_fit$settings$num_basis,
  betaprec = source_fit$settings$betaprec,
  order = source_fit$settings$order,
  pred_step = source_fit$settings$pred_step
)
maximum_likelihood_error <- max(abs(
  analytic_log_likelihood - tmb_log_likelihood
))

r2_component_diagnostics <- do.call(rbind, lapply(
  seq_along(alternative_covariances),
  function(i) {
    relative_eigenvalues <- Re(eigen(
      solve(null_covariance, alternative_covariances[[i]]),
      only.values = TRUE
    )$values)
    second_moment <- gaussian_bf_second_moment(
      null_covariance,
      alternative_covariances[[i]]
    )
    data.frame(
      scenario = "Actual R2 fitted alternative",
      psd = alternative_psd[i],
      conditional_alternative_weight = alternative_weights[i],
      minimum_relative_covariance_eigenvalue =
        min(relative_eigenvalues),
      maximum_relative_covariance_eigenvalue =
        max(relative_eigenvalues),
      n_relative_eigenvalues_at_least_two =
        sum(relative_eigenvalues >= 2),
      theoretical_bf_second_moment = second_moment,
      finite_bf_second_moment = is.finite(second_moment),
      stringsAsFactors = FALSE
    )
  }
))

# A positive-weight mixture has infinite second moment if any component has an
# infinite second moment. The R2 diagnostic determines whether the planned
# finite-variance control is required.
r2_mixture_has_finite_second_moment <- all(
  r2_component_diagnostics$finite_bf_second_moment
)

target_maximum_relative_variance <- 1.5
broadest_component <- which.max(
  r2_component_diagnostics$maximum_relative_covariance_eigenvalue
)
broadest_covariance <- alternative_covariances[[broadest_component]]
broadest_maximum <- r2_component_diagnostics[
  broadest_component,
  "maximum_relative_covariance_eigenvalue"
]
covariance_increment_scale <-
  (target_maximum_relative_variance - 1) /
  (broadest_maximum - 1)
mild_covariance <- null_covariance +
  covariance_increment_scale *
    (broadest_covariance - null_covariance)
mild_relative_eigenvalues <- Re(eigen(
  solve(null_covariance, mild_covariance),
  only.values = TRUE
)$values)
mild_second_moment <- gaussian_bf_second_moment(
  null_covariance,
  mild_covariance
)
mild_component_diagnostics <- data.frame(
  scenario = "Prespecified mild one-component control",
  psd = NA_real_,
  conditional_alternative_weight = 1,
  minimum_relative_covariance_eigenvalue =
    min(mild_relative_eigenvalues),
  maximum_relative_covariance_eigenvalue =
    max(mild_relative_eigenvalues),
  n_relative_eigenvalues_at_least_two =
    sum(mild_relative_eigenvalues >= 2),
  theoretical_bf_second_moment = mild_second_moment,
  finite_bf_second_moment = is.finite(mild_second_moment),
  stringsAsFactors = FALSE
)
component_diagnostics <- rbind(
  r2_component_diagnostics,
  mild_component_diagnostics
)

simulation_seed <- 20260814L
n_simulations <- 1000000L
batch_size <- 100000L
checkpoints <- unique(as.integer(round(10 ^ seq(
  log10(100), log10(n_simulations), length.out = 100L
))))

r2_simulation <- simulate_exact_null_bf(
  null_covariance = null_covariance,
  alternative_covariances = alternative_covariances,
  alternative_weights = alternative_weights,
  scenario = "Actual R2 fitted alternative",
  seed = simulation_seed,
  n_simulations = n_simulations,
  batch_size = batch_size,
  checkpoints = checkpoints
)
mild_simulation <- simulate_exact_null_bf(
  null_covariance = null_covariance,
  alternative_covariances = list(mild_covariance),
  alternative_weights = 1,
  scenario = "Prespecified mild one-component control",
  seed = simulation_seed,
  n_simulations = n_simulations,
  batch_size = batch_size,
  checkpoints = checkpoints
)

bf_summary <- rbind(
  summarize_bf(
    r2_simulation$bayes_factor,
    "Actual R2 fitted alternative",
    if (r2_mixture_has_finite_second_moment) NA_real_ else Inf
  ),
  summarize_bf(
    mild_simulation$bayes_factor,
    "Prespecified mild one-component control",
    mild_second_moment
  )
)
running_summary <- rbind(
  r2_simulation$running,
  mild_simulation$running
)

mild_standard_error <- sqrt(
  (mild_second_moment - 1) / n_simulations
)
validation <- data.frame(
  check = c(
    "Analytic likelihood reproduces FASH/TMB",
    "Actual R2 mixture has infinite BF second moment",
    "Mild control maximum relative variance is 1.5",
    "Mild control has finite BF second moment",
    "Mild ordinary sample mean is within four theoretical SEs of one"
  ),
  observed = c(
    maximum_likelihood_error,
    as.numeric(!r2_mixture_has_finite_second_moment),
    max(mild_relative_eigenvalues),
    as.numeric(is.finite(mild_second_moment)),
    mean(mild_simulation$bayes_factor)
  ),
  reference = c(1e-7, 1, 1.5, 1, 1),
  pass = c(
    maximum_likelihood_error <= 1e-7,
    !r2_mixture_has_finite_second_moment,
    abs(max(mild_relative_eigenvalues) - 1.5) <= 1e-8,
    is.finite(mild_second_moment),
    abs(mean(mild_simulation$bayes_factor) - 1) <=
      4 * mild_standard_error
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$pass)) {
  stop("At least one mild-setting validation failed.")
}

utils::write.csv(
  component_diagnostics,
  file.path(output_directory, "component_diagnostics.csv"),
  row.names = FALSE
)
utils::write.csv(
  bf_summary,
  file.path(output_directory, "bf_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  running_summary,
  file.path(output_directory, "running_mean_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  validation,
  file.path(output_directory, "validation.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    settings = list(
      source_path = source_path,
      simulation_seed = simulation_seed,
      n_simulations = n_simulations,
      batch_size = batch_size,
      representative_se = representative_se,
      alternative_psd = alternative_psd,
      alternative_weights = alternative_weights,
      target_maximum_relative_variance =
        target_maximum_relative_variance,
      covariance_increment_scale = covariance_increment_scale
    ),
    component_diagnostics = component_diagnostics,
    bf_summary = bf_summary,
    running_summary = running_summary,
    validation = validation
  ),
  file.path(output_directory, "results.rds")
)

running_figure <- ggplot2::ggplot(
  running_summary,
  ggplot2::aes(x = n, y = running_mean_bf, color = scenario)
) +
  ggplot2::geom_hline(
    yintercept = 1,
    color = "#777777",
    linetype = "dotted"
  ) +
  ggplot2::geom_line(linewidth = 0.95) +
  ggplot2::scale_x_log10(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  ggplot2::scale_color_manual(values = c(
    "Actual R2 fitted alternative" = "#D55E00",
    "Prespecified mild one-component control" = "#0072B2"
  )) +
  ggplot2::labs(
    title = "Ordinary exact-null BF sample means",
    subtitle = paste(
      "Same exact FASH null and seed; the mild control has",
      "maximum relative variance 1.5 and a finite BF variance."
    ),
    x = "Number of exact-null draws (log scale)",
    y = "Running sample mean BF",
    color = NULL
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )
ggplot2::ggsave(
  file.path(figure_directory, "ordinary_null_bf_convergence.png"),
  running_figure,
  width = 8.4,
  height = 5.2,
  dpi = 300,
  bg = "white"
)

message("Wrote R2/mild exact-null validation to: ", output_directory)
