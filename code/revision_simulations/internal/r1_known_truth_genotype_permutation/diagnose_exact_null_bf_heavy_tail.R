#!/usr/bin/env Rscript

# Diagnose finite-sample convergence of the BF mean under the exact Gaussian
# FASH null with a fixed conditional alternative distribution.

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
  chol_covariance <- chol(covariance)
  standardized <- forwardsolve(t(chol_covariance), y)
  -0.5 * (
    length(y) * log(2 * pi) +
      2 * sum(log(diag(chol_covariance))) +
      sum(standardized^2)
  )
}

summarize_bf <- function(values) {
  probabilities <- c(0.50, 0.90, 0.95, 0.99, 0.999, 0.9999)
  quantiles <- stats::quantile(
    values, probs = probabilities, type = 8, names = FALSE
  )
  descending <- order(values, decreasing = TRUE)
  data.frame(
    n = length(values),
    mean = mean(values),
    median = quantiles[1L],
    q90 = quantiles[2L],
    q95 = quantiles[3L],
    q99 = quantiles[4L],
    q999 = quantiles[5L],
    q9999 = quantiles[6L],
    maximum = max(values),
    proportion_greater_than_one = mean(values > 1),
    top_10_mass_share = sum(values[descending[seq_len(10L)]]) / sum(values),
    stringsAsFactors = FALSE
  )
}

summarize_block_means <- function(values,
                                  block_size,
                                  observed_reference) {
  n_blocks <- length(values) %/% block_size
  used <- values[seq_len(n_blocks * block_size)]
  block_means <- rowMeans(matrix(
    used, nrow = n_blocks, ncol = block_size, byrow = TRUE
  ))
  quantiles <- stats::quantile(
    block_means,
    probs = c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99),
    type = 8,
    names = FALSE
  )
  list(
    summary = data.frame(
      block_size = block_size,
      n_blocks = n_blocks,
      observed_reference = observed_reference,
      proportion_at_or_below_observed = mean(
        block_means <= observed_reference
      ),
      mean_of_block_means = mean(block_means),
      q01 = quantiles[1L],
      q05 = quantiles[2L],
      q25 = quantiles[3L],
      median = quantiles[4L],
      q75 = quantiles[5L],
      q95 = quantiles[6L],
      q99 = quantiles[7L],
      maximum = max(block_means),
      stringsAsFactors = FALSE
    ),
    block_means = data.frame(
      block_size = block_size,
      block = seq_along(block_means),
      mean_bf = block_means,
      stringsAsFactors = FALSE
    )
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
output_id <- "r1_exact_null_bf_heavy_tail_v1"
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
likelihood_validation <- data.frame(
  maximum_absolute_log_likelihood_error = max(abs(
    analytic_log_likelihood - tmb_log_likelihood
  )),
  tolerance = 1e-7,
  stringsAsFactors = FALSE
)
likelihood_validation$pass <- with(
  likelihood_validation,
  maximum_absolute_log_likelihood_error <= tolerance
)
if (!likelihood_validation$pass) {
  stop("The analytic Gaussian likelihood did not reproduce FASH.")
}

second_moment_rows <- lapply(seq_along(alternative_covariances), function(i) {
  relative_eigenvalues <- Re(eigen(
    solve(null_covariance, alternative_covariances[[i]]),
    only.values = TRUE
  )$values)
  second_moment_precision <-
    2 * solve(alternative_covariances[[i]]) - solve(null_covariance)
  minimum_precision_eigenvalue <- min(eigen(
    (second_moment_precision + t(second_moment_precision)) / 2,
    symmetric = TRUE,
    only.values = TRUE
  )$values)
  data.frame(
    psd = alternative_psd[i],
    conditional_alternative_weight = alternative_weights[i],
    minimum_relative_covariance_eigenvalue = min(relative_eigenvalues),
    maximum_relative_covariance_eigenvalue = max(relative_eigenvalues),
    n_relative_eigenvalues_at_least_two = sum(relative_eigenvalues >= 2),
    minimum_second_moment_precision_eigenvalue =
      minimum_precision_eigenvalue,
    finite_bf_second_moment = minimum_precision_eigenvalue > 0,
    stringsAsFactors = FALSE
  )
})
second_moment_diagnostics <- do.call(rbind, second_moment_rows)

simulation_seed <- 20260812L
n_simulations <- 1000000L
batch_size <- 100000L
set.seed(simulation_seed)
null_cholesky <- chol(null_covariance)
null_log_determinant <- 2 * sum(log(diag(null_cholesky)))
alternative_cholesky <- lapply(alternative_covariances, chol)
alternative_log_determinants <- vapply(
  alternative_cholesky,
  function(cholesky) 2 * sum(log(diag(cholesky))),
  numeric(1)
)
alternative_cholesky_inverse <- lapply(alternative_cholesky, solve)

bayes_factor <- numeric(n_simulations)
running_rows <- list()
position <- 0L
for (batch in seq_len(n_simulations / batch_size)) {
  z <- matrix(
    stats::rnorm(batch_size * nrow(null_covariance)),
    nrow = batch_size,
    ncol = nrow(null_covariance)
  )
  y <- z %*% null_cholesky
  null_quadratic <- rowSums(z^2)
  log_component_bf <- sapply(seq_along(alternative_covariances), function(i) {
    alternative_standardized <-
      y %*% alternative_cholesky_inverse[[i]]
    alternative_quadratic <- rowSums(alternative_standardized^2)
    -0.5 * (
      alternative_log_determinants[i] - null_log_determinant +
        alternative_quadratic - null_quadratic
    )
  })
  if (is.null(dim(log_component_bf))) {
    log_component_bf <- matrix(log_component_bf, ncol = 1L)
  }
  row_maximum <- apply(log_component_bf, 1L, max)
  batch_bf <- exp(row_maximum) * rowSums(
    sweep(exp(log_component_bf - row_maximum), 2L,
          alternative_weights, `*`)
  )
  indices <- position + seq_len(batch_size)
  bayes_factor[indices] <- batch_bf
  position <- position + batch_size
  running_rows[[batch]] <- data.frame(
    n = position,
    running_mean_bf = mean(bayes_factor[seq_len(position)]),
    running_maximum_bf = max(bayes_factor[seq_len(position)]),
    stringsAsFactors = FALSE
  )
}
running_summary <- do.call(rbind, running_rows)
bf_summary <- summarize_bf(bayes_factor)

block_200 <- summarize_block_means(
  bayes_factor,
  block_size = 200L,
  observed_reference = 0.229589081068532
)
block_1000 <- summarize_block_means(
  bayes_factor,
  block_size = 1000L,
  observed_reference = 0.194743204532263
)
block_summary <- rbind(block_200$summary, block_1000$summary)
block_means <- rbind(block_200$block_means, block_1000$block_means)

utils::write.csv(
  likelihood_validation,
  file.path(output_directory, "likelihood_validation.csv"),
  row.names = FALSE
)
utils::write.csv(
  second_moment_diagnostics,
  file.path(output_directory, "second_moment_diagnostics.csv"),
  row.names = FALSE
)
utils::write.csv(
  running_summary,
  file.path(output_directory, "running_mean_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  bf_summary,
  file.path(output_directory, "exact_null_bf_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  block_summary,
  file.path(output_directory, "block_mean_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  block_means,
  file.path(output_directory, "block_means.csv"),
  row.names = FALSE
)
saveRDS(
  list(
    settings = list(
      simulation_seed = simulation_seed,
      n_simulations = n_simulations,
      batch_size = batch_size,
      representative_se = representative_se,
      alternative_psd = alternative_psd,
      alternative_weights = alternative_weights,
      theoretical_mean_bf = 1
    ),
    likelihood_validation = likelihood_validation,
    second_moment_diagnostics = second_moment_diagnostics,
    running_summary = running_summary,
    bf_summary = bf_summary,
    block_summary = block_summary
  ),
  file.path(output_directory, "results.rds")
)

running_figure <- ggplot2::ggplot(
  running_summary,
  ggplot2::aes(x = n, y = running_mean_bf)
) +
  ggplot2::geom_hline(
    yintercept = 1, color = "#777777", linetype = "dotted"
  ) +
  ggplot2::geom_line(color = "#0072B2", linewidth = 1.0) +
  ggplot2::geom_point(color = "#0072B2", size = 2.0) +
  ggplot2::scale_x_continuous(
    labels = scales::label_number(scale_cut = scales::cut_short_scale())
  ) +
  ggplot2::labs(
    title = "Exact-null BF sample mean converges extremely slowly",
    subtitle = paste(
      "Fixed FASH alternative; exact Gaussian p0; theoretical mean = 1;",
      "BF second moment is infinite."
    ),
    x = "Number of exact-null draws",
    y = "Running sample mean BF"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
ggplot2::ggsave(
  file.path(figure_directory, "running_mean_bf.png"),
  running_figure, width = 8.2, height = 5.0, dpi = 300, bg = "white"
)

message("Wrote exact-null BF heavy-tail diagnostic to: ", output_directory)
