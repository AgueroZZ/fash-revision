#!/usr/bin/env Rscript

# Demonstrate how the exact-null BF expectation is recovered when the rare
# upper tail is sampled deliberately, and quantify why ordinary null Monte
# Carlo cannot display the same convergence at feasible sample sizes.

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

row_log_sum_exp <- function(log_values, log_weights) {
  weighted <- sweep(log_values, 2L, log_weights, `+`)
  row_maximum <- apply(weighted, 1L, max)
  row_maximum + log(rowSums(exp(weighted - row_maximum)))
}

log_mvn_density_rows <- function(y, cholesky, log_determinant) {
  standardized <- y %*% solve(cholesky)
  -0.5 * (
    ncol(y) * log(2 * pi) +
      log_determinant +
      rowSums(standardized^2)
  )
}

draw_mvn_rows <- function(n, cholesky) {
  matrix(
    stats::rnorm(n * nrow(cholesky)),
    nrow = n,
    ncol = nrow(cholesky)
  ) %*% cholesky
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
  "r1_random_bspline_main_effect_profile_sigma_pilot5",
  "full_fits", "seed_12345.rds"
)
previous_output_directory <- file.path(
  workflowr_root, "output", "revision_simulations", "internal",
  "r1_exact_null_bf_heavy_tail_v1"
)
output_id <- "r1_exact_null_bf_convergence_v1"
output_directory <- file.path(
  workflowr_root, "output", "revision_simulations", "internal", output_id
)
figure_directory <- file.path(output_directory, "figures")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

source_object <- readRDS(source_path)
source_fit <- source_object$fash_fits$fash_iwp1_raw
effect_class <- as.character(source_object$unit_info$effect_class)
alternative_indices <- which(effect_class == "dynamic_bspline")[seq_len(200L)]
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
  stop("The reference fit has no active alternative components.")
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

null_cholesky <- chol(null_covariance)
alternative_cholesky <- lapply(alternative_covariances, chol)
null_log_determinant <- 2 * sum(log(diag(null_cholesky)))
alternative_log_determinants <- vapply(
  alternative_cholesky,
  function(cholesky) 2 * sum(log(diag(cholesky))),
  numeric(1)
)

compute_log_bf <- function(y) {
  null_log_density <- log_mvn_density_rows(
    y, null_cholesky, null_log_determinant
  )
  alternative_log_density <- vapply(
    seq_along(alternative_cholesky),
    function(i) {
      log_mvn_density_rows(
        y,
        alternative_cholesky[[i]],
        alternative_log_determinants[i]
      )
    },
    numeric(nrow(y))
  )
  if (is.null(dim(alternative_log_density))) {
    alternative_log_density <- matrix(
      alternative_log_density, ncol = 1L
    )
  }
  row_log_sum_exp(
    alternative_log_density,
    log(alternative_weights)
  ) - null_log_density
}

# Use a mixture proposal q = (1 - rho) p0 + rho p1. The importance-weighted
# integrand for E0[BF] is bounded above by 1 / rho:
# BF * p0 / q = BF / ((1 - rho) + rho * BF).
importance_seed <- 20260813L
importance_rho <- 0.5
importance_n <- 1000000L
importance_batch_size <- 50000L
set.seed(importance_seed)
importance_estimate <- numeric(importance_n)
position <- 0L
while (position < importance_n) {
  current_n <- min(importance_batch_size, importance_n - position)
  proposal_source <- sample.int(
    length(alternative_weights) + 1L,
    size = current_n,
    replace = TRUE,
    prob = c(
      1 - importance_rho,
      importance_rho * alternative_weights
    )
  )
  y <- matrix(
    NA_real_,
    nrow = current_n,
    ncol = nrow(null_covariance)
  )
  for (source_index in sort(unique(proposal_source))) {
    selected <- which(proposal_source == source_index)
    selected_cholesky <- if (source_index == 1L) {
      null_cholesky
    } else {
      alternative_cholesky[[source_index - 1L]]
    }
    y[selected, ] <- draw_mvn_rows(length(selected), selected_cholesky)
  }
  log_bf <- compute_log_bf(y)
  log_denominator <- row_log_sum_exp(
    cbind(
      rep(log(1 - importance_rho), current_n),
      log(importance_rho) + log_bf
    ),
    c(0, 0)
  )
  indices <- position + seq_len(current_n)
  importance_estimate[indices] <- exp(log_bf - log_denominator)
  position <- position + current_n
}

importance_checkpoints <- unique(as.integer(round(10 ^ seq(
  log10(100), log10(importance_n), length.out = 50L
))))
importance_cumulative_sum <- cumsum(importance_estimate)
importance_cumulative_sum_squares <- cumsum(importance_estimate^2)
importance_running_mean <-
  importance_cumulative_sum[importance_checkpoints] /
  importance_checkpoints
importance_running_variance <- pmax(
  (
    importance_cumulative_sum_squares[importance_checkpoints] -
      importance_checkpoints * importance_running_mean^2
  ) / pmax(importance_checkpoints - 1L, 1L),
  0
)
importance_running <- data.frame(
  estimator = "Tail-aware mixture importance sampling",
  n = importance_checkpoints,
  estimate = importance_running_mean,
  standard_error = sqrt(
    importance_running_variance / importance_checkpoints
  ),
  stringsAsFactors = FALSE
)
importance_running$lower_95 <-
  importance_running$estimate - 1.96 * importance_running$standard_error
importance_running$upper_95 <-
  importance_running$estimate + 1.96 * importance_running$standard_error

previous_running <- utils::read.csv(file.path(
  previous_output_directory, "running_mean_summary.csv"
))
ordinary_running <- data.frame(
  estimator = "Ordinary IID sampling from p0",
  n = previous_running$n,
  estimate = previous_running$running_mean_bf,
  standard_error = NA_real_,
  lower_95 = NA_real_,
  upper_95 = NA_real_,
  stringsAsFactors = FALSE
)
running_comparison <- rbind(ordinary_running, importance_running)

# Quantify the p0 sample size required to resolve the BF expectation along the
# broadest single Gaussian alternative direction. If lambda is the alternative
# variance relative to p0 in that direction, then the fraction r of that
# component's mean contributed by |Z| <= c is
# P_{N(0, lambda)}(|Z| <= c) = r.
relative_eigenvalues <- lapply(alternative_covariances, function(covariance) {
  Re(eigen(
    solve(null_covariance, covariance),
    only.values = TRUE
  )$values)
})
maximum_relative_eigenvalues <- vapply(
  relative_eigenvalues, max, numeric(1)
)
leading_component <- which.max(maximum_relative_eigenvalues)
leading_lambda <- maximum_relative_eigenvalues[leading_component]
leading_weight <- alternative_weights[leading_component]

captured_grid <- seq(0.01, 0.999, length.out = 500L)
boundary_grid <- sqrt(leading_lambda) * stats::qnorm(
  (1 + captured_grid) / 2
)
log_p0_tail_grid <- log(2) + stats::pnorm(
  -boundary_grid, log.p = TRUE
)
tail_resolution_curve <- data.frame(
  captured_component_mean = captured_grid,
  missing_component_mean = 1 - captured_grid,
  boundary = boundary_grid,
  log10_expected_p0_draws_for_one_tail_event =
    -log_p0_tail_grid / log(10),
  leading_relative_variance = leading_lambda,
  leading_component_weight = leading_weight,
  stringsAsFactors = FALSE
)

captured_targets <- c(0.50, 0.75, 0.90, 0.95, 0.99)
target_boundaries <- sqrt(leading_lambda) * stats::qnorm(
  (1 + captured_targets) / 2
)
target_log_p0_tail <- log(2) + stats::pnorm(
  -target_boundaries, log.p = TRUE
)
tail_resolution_targets <- data.frame(
  captured_component_mean = captured_targets,
  missing_component_mean = 1 - captured_targets,
  leading_component_weight = leading_weight,
  weighted_missing_mean_contribution =
    leading_weight * (1 - captured_targets),
  boundary = target_boundaries,
  log10_expected_p0_draws_for_one_tail_event =
    -target_log_p0_tail / log(10),
  stringsAsFactors = FALSE
)

importance_summary <- data.frame(
  estimator = "Tail-aware mixture importance sampling",
  n = importance_n,
  estimate = mean(importance_estimate),
  standard_error = stats::sd(importance_estimate) / sqrt(importance_n),
  lower_95 = mean(importance_estimate) -
    1.96 * stats::sd(importance_estimate) / sqrt(importance_n),
  upper_95 = mean(importance_estimate) +
    1.96 * stats::sd(importance_estimate) / sqrt(importance_n),
  minimum = min(importance_estimate),
  maximum = max(importance_estimate),
  theoretical_upper_bound = 1 / importance_rho,
  stringsAsFactors = FALSE
)

validation <- data.frame(
  check = c(
    "Importance estimate contains theoretical mean",
    "Importance weights respect theoretical upper bound",
    "Ordinary path reproduces previous one-million-draw result"
  ),
  observed = c(
    importance_summary$estimate,
    importance_summary$maximum,
    ordinary_running$estimate[nrow(ordinary_running)]
  ),
  reference = c(1, 1 / importance_rho, 0.263764996595593),
  pass = c(
    importance_summary$lower_95 <= 1 &&
      importance_summary$upper_95 >= 1,
    importance_summary$maximum <= (1 / importance_rho) + 1e-10,
    abs(
      ordinary_running$estimate[nrow(ordinary_running)] -
        0.263764996595593
    ) < 1e-12
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$pass)) {
  stop("At least one convergence diagnostic validation failed.")
}

utils::write.csv(
  running_comparison,
  file.path(output_directory, "running_mean_comparison.csv"),
  row.names = FALSE
)
utils::write.csv(
  importance_summary,
  file.path(output_directory, "importance_sampling_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  tail_resolution_curve,
  file.path(output_directory, "tail_resolution_curve.csv"),
  row.names = FALSE
)
utils::write.csv(
  tail_resolution_targets,
  file.path(output_directory, "tail_resolution_targets.csv"),
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
      importance_seed = importance_seed,
      importance_rho = importance_rho,
      importance_n = importance_n,
      alternative_psd = alternative_psd,
      alternative_weights = alternative_weights,
      maximum_relative_eigenvalues = maximum_relative_eigenvalues,
      theoretical_mean_bf = 1
    ),
    running_comparison = running_comparison,
    importance_summary = importance_summary,
    tail_resolution_curve = tail_resolution_curve,
    tail_resolution_targets = tail_resolution_targets,
    validation = validation
  ),
  file.path(output_directory, "results.rds")
)

running_figure <- ggplot2::ggplot(
  running_comparison,
  ggplot2::aes(x = n, y = estimate, color = estimator)
) +
  ggplot2::geom_hline(
    yintercept = 1, color = "#777777", linetype = "dotted"
  ) +
  ggplot2::geom_ribbon(
    data = importance_running,
    ggplot2::aes(
      ymin = lower_95,
      ymax = upper_95,
      fill = estimator
    ),
    color = NA,
    alpha = 0.18,
    show.legend = FALSE
  ) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(
    data = ordinary_running,
    size = 1.7
  ) +
  ggplot2::scale_x_log10(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  ggplot2::scale_color_manual(values = c(
    "Ordinary IID sampling from p0" = "#D55E00",
    "Tail-aware mixture importance sampling" = "#0072B2"
  )) +
  ggplot2::scale_fill_manual(values = c(
    "Tail-aware mixture importance sampling" = "#0072B2"
  )) +
  ggplot2::labs(
    title = "The BF mean is recovered when the rare tail is represented",
    subtitle = paste(
      "Both curves target E0(BF) = 1; importance sampling is bounded",
      "and is not an ordinary p0 sample path."
    ),
    x = "Number of Monte Carlo draws (log scale)",
    y = "Running estimate of E0(BF)",
    color = NULL
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )
ggplot2::ggsave(
  file.path(figure_directory, "ordinary_vs_tail_aware_convergence.png"),
  running_figure,
  width = 8.3,
  height = 5.2,
  dpi = 300,
  bg = "white"
)

tail_figure <- ggplot2::ggplot(
  tail_resolution_curve,
  ggplot2::aes(
    x = log10_expected_p0_draws_for_one_tail_event,
    y = captured_component_mean
  )
) +
  ggplot2::geom_hline(
    yintercept = 1, color = "#777777", linetype = "dotted"
  ) +
  ggplot2::geom_line(color = "#009E73", linewidth = 1.0) +
  ggplot2::geom_point(
    data = tail_resolution_targets,
    color = "#009E73",
    size = 2.0
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1.01),
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Why an ordinary null sample path cannot visibly reach one",
    subtitle = paste0(
      "Broadest alternative direction has variance ratio ",
      format(round(leading_lambda, 1), nsmall = 1),
      ".\nThe x-axis is the null draws needed for one event beyond the boundary."
    ),
    x = expression(log[10] * " expected p0 draws for one tail event"),
    y = "Captured mean of the leading BF component"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(
  file.path(figure_directory, "ordinary_sampling_tail_resolution.png"),
  tail_figure,
  width = 8.3,
  height = 5.2,
  dpi = 300,
  bg = "white"
)

message("Wrote exact-null BF convergence diagnostic to: ", output_directory)
