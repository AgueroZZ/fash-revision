#!/usr/bin/env Rscript

# Estimate cross-time covariance from matched synchronized and independent-time
# donor residual-block permutations.

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

matrix_key <- function(threshold_name, se_scale, estimator) {
  paste(threshold_name, se_scale, estimator, sep = "__")
}

estimate_direct_zero_matrices <- function(beta_hat, se) {
  z_score <- beta_hat / se
  covariance <- crossprod(z_score) / nrow(z_score)
  correlation <- stats::cov2cor(covariance)
  list(covariance = covariance, correlation = correlation)
}

estimate_within_unit_centered_matrices <- function(beta_hat, se) {
  centered <- weighted_center_standardize(beta_hat, se)$standardized_residual
  covariance <- crossprod(centered) / nrow(centered)
  correlation <- stats::cor(centered)
  list(covariance = covariance, correlation = correlation)
}

semivariogram_matrix <- function(covariance) {
  covariance <- as.matrix(covariance)
  diagonal <- diag(covariance)
  semivariogram <- outer(diagonal, diagonal, "+") / 2 - covariance
  diag(semivariogram) <- 0
  semivariogram
}

lag_average_matrix <- function(matrix) {
  n_time <- ncol(matrix)
  vapply(seq_len(n_time - 1L), function(lag) {
    mean(matrix[cbind(
      seq_len(n_time - lag),
      (lag + 1L):n_time
    )])
  }, numeric(1))
}

aggregate_values <- function(data, value_name, grouping_names) {
  split_key <- interaction(data[grouping_names], drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(data)), split_key)
  output <- do.call(rbind, lapply(groups, function(index) {
    values <- data[[value_name]][index]
    cbind(
      data[index[1L], grouping_names, drop = FALSE],
      data.frame(
        metric = value_name,
        mean = mean(values),
        median = stats::median(values),
        lower = unname(stats::quantile(values, 0.025, names = FALSE)),
        upper = unname(stats::quantile(values, 0.975, names = FALSE)),
        sd = stats::sd(values),
        n_replications = length(values),
        stringsAsFactors = FALSE
      )
    )
  }))
  rownames(output) <- NULL
  output
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "r4_correlated_errors",
  "real_data_correlation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "covariance_estimation",
  "donor_null_permutation_helpers.R"
))

pilot_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "donor_null_permutation_pilot"
)
comparison_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "null_permutation_estimator_comparison"
)
input_path <- file.path(pilot_dir, "input", "selected_raw_data.rds")
comparison_path <- file.path(comparison_dir, "estimator_comparison.rds")
if (!file.exists(input_path) || !file.exists(comparison_path)) {
  stop("The donor-permutation input or synchronized comparison cache is missing.")
}

input <- readRDS(input_path)
synchronized <- readRDS(comparison_path)
if (!identical(
  synchronized$configuration$analysis_id,
  "null_permutation_estimator_comparison"
) || synchronized$configuration$n_replications != 100L) {
  stop("The synchronized estimator-comparison cache is incompatible.")
}

output_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "permuted_variogram_covariance_summary"
)
summary_dir <- file.path(output_dir, "summary")
figure_dir <- file.path(output_dir, "figure")
invisible(lapply(
  c(output_dir, summary_dir, figure_dir),
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

threshold_names <- names(input$threshold_indices)
threshold_values <- input$threshold_values
threshold_counts <- lengths(input$threshold_indices)
se_scales <- c("raw_regression", "t_adjusted")
estimator_names <- c(
  "pairwise_ols",
  "within_unit_centered",
  "direct_zero_null"
)
estimator_labels <- c(
  pairwise_ols = "Pairwise OLS",
  within_unit_centered = "Within-unit centered",
  direct_zero_null = "Direct zero-null"
)
time_grid <- input$time_grid
n_time <- length(time_grid)
n_replications <- synchronized$configuration$n_replications
independent_seeds <- 20560831L + seq_len(n_replications) - 1L

key_grid <- expand.grid(
  threshold_name = threshold_names,
  se_scale = se_scales,
  estimator = estimator_names,
  stringsAsFactors = FALSE
)
keys <- mapply(
  matrix_key,
  key_grid$threshold_name,
  key_grid$se_scale,
  key_grid$estimator,
  USE.NAMES = FALSE
)
make_draw_array <- function() {
  array(
    NA_real_,
    dim = c(n_time, n_time, n_replications),
    dimnames = list(
      paste0("time_", time_grid),
      paste0("time_", time_grid),
      paste0("rep_", seq_len(n_replications))
    )
  )
}
independent_covariance_draws <- setNames(lapply(keys, function(unused) {
  make_draw_array()
}), keys)

fit_independent_time_permutation <- function(seed) {
  set.seed(seed)
  beta_hat <- matrix(
    NA_real_,
    nrow = nrow(input$unit_table),
    ncol = n_time,
    dimnames = list(
      input$unit_table$pair_key,
      paste0("time_", time_grid)
    )
  )
  raw_se <- beta_hat
  residual_df <- integer(n_time)
  for (time_index in seq_along(input$time_inputs)) {
    time_input <- input$time_inputs[[time_index]]
    source_rows <- sample.int(
      nrow(time_input$expression_residual),
      replace = FALSE
    )
    permuted_residual <- time_input$expression_residual[
      source_rows,
      ,
      drop = FALSE
    ]
    expression_residual <- time_input$projection$residualizer %*%
      permuted_residual
    genotype <- input$unit_dosage[
      time_input$donors,
      ,
      drop = FALSE
    ]
    fit <- fit_residualized_genotype_regressions(
      expression_residual,
      genotype,
      time_input$projection$residualizer,
      time_input$projection$rank
    )
    beta_hat[, time_index] <- fit$beta
    raw_se[, time_index] <- fit$standard_error
    residual_df[time_index] <- fit$residual_df
  }
  list(
    beta_hat = beta_hat,
    raw_se = raw_se,
    adjusted_se = convert_raw_to_t_adjusted_se(
      beta_hat,
      raw_se,
      residual_df
    )
  )
}

message("Running 100 independent-time donor residual permutations.")
for (replication in seq_len(n_replications)) {
  fit <- fit_independent_time_permutation(independent_seeds[replication])
  for (threshold_index in seq_along(threshold_names)) {
    selected_index <- input$threshold_indices[[threshold_index]]
    selected_beta <- fit$beta_hat[selected_index, , drop = FALSE]
    for (se_scale in se_scales) {
      selected_se <- if (se_scale == "raw_regression") {
        fit$raw_se[selected_index, , drop = FALSE]
      } else {
        fit$adjusted_se[selected_index, , drop = FALSE]
      }
      centered <- estimate_within_unit_centered_matrices(
        selected_beta,
        selected_se
      )
      direct <- estimate_direct_zero_matrices(selected_beta, selected_se)
      pairwise <- estimate_pairwise_difference_correlation(
        selected_beta,
        selected_se
      )
      estimator_results <- list(
        pairwise_ols = list(covariance = pairwise),
        within_unit_centered = centered,
        direct_zero_null = direct
      )
      for (estimator in estimator_names) {
        key <- matrix_key(
          threshold_names[threshold_index],
          se_scale,
          estimator
        )
        independent_covariance_draws[[key]][, , replication] <-
          estimator_results[[estimator]]$covariance
      }
    }
  }
  if (replication %% 10L == 0L) {
    message("Completed ", replication, " of ", n_replications, ".")
  }
}

synchronized_variogram_draws <- setNames(lapply(keys, function(unused) {
  make_draw_array()
}), keys)
independent_variogram_draws <- setNames(lapply(keys, function(unused) {
  make_draw_array()
}), keys)
covariance_signal_draws <- setNames(lapply(keys, function(unused) {
  make_draw_array()
}), keys)

for (grid_index in seq_len(nrow(key_grid))) {
  key <- keys[grid_index]
  direct_key <- matrix_key(
    key_grid$threshold_name[grid_index],
    key_grid$se_scale[grid_index],
    "direct_zero_null"
  )
  for (replication in seq_len(n_replications)) {
    synchronized_covariance <-
      synchronized$covariance_draws[[key]][, , replication]
    independent_covariance <-
      independent_covariance_draws[[key]][, , replication]
    synchronized_variogram <- semivariogram_matrix(
      synchronized_covariance
    )
    independent_variogram <- semivariogram_matrix(
      independent_covariance
    )
    covariance_signal <- independent_variogram - synchronized_variogram
    diagonal_reference <- diag(
      synchronized$covariance_draws[[direct_key]][, , replication]
    )
    diag(covariance_signal) <- diagonal_reference
    synchronized_variogram_draws[[key]][, , replication] <-
      synchronized_variogram
    independent_variogram_draws[[key]][, , replication] <-
      independent_variogram
    covariance_signal_draws[[key]][, , replication] <- covariance_signal
  }
}

lag_rows <- vector(
  "list",
  nrow(key_grid) * n_replications
)
output_index <- 1L
for (grid_index in seq_len(nrow(key_grid))) {
  key <- keys[grid_index]
  threshold_index <- match(
    key_grid$threshold_name[grid_index],
    threshold_names
  )
  for (replication in seq_len(n_replications)) {
    synchronized_lag <- lag_average_matrix(
      synchronized_variogram_draws[[key]][, , replication]
    )
    independent_lag <- lag_average_matrix(
      independent_variogram_draws[[key]][, , replication]
    )
    covariance_lag <- lag_average_matrix(
      covariance_signal_draws[[key]][, , replication]
    )
    lag_rows[[output_index]] <- data.frame(
      replication = replication,
      synchronized_seed = synchronized$configuration$seeds[replication],
      independent_seed = independent_seeds[replication],
      threshold = unname(threshold_values[threshold_index]),
      n_units = unname(threshold_counts[threshold_index]),
      se_scale = key_grid$se_scale[grid_index],
      estimator = key_grid$estimator[grid_index],
      lag = seq_along(synchronized_lag),
      synchronized_variogram = synchronized_lag,
      independent_variogram = independent_lag,
      covariance = covariance_lag,
      stringsAsFactors = FALSE
    )
    output_index <- output_index + 1L
  }
}
replicate_lag_summaries <- do.call(rbind, lag_rows)
rownames(replicate_lag_summaries) <- NULL

lag_grouping <- c(
  "threshold", "n_units", "se_scale", "estimator", "lag"
)
aggregate_lag_summaries <- do.call(rbind, lapply(
  c("synchronized_variogram", "independent_variogram", "covariance"),
  function(metric) {
    aggregate_values(replicate_lag_summaries, metric, lag_grouping)
  }
))
rownames(aggregate_lag_summaries) <- NULL

contrast_rows <- list()
contrast_index <- 1L
for (threshold in threshold_values) {
  for (se_scale in se_scales) {
    reference <- replicate_lag_summaries[
      replicate_lag_summaries$threshold == threshold &
        replicate_lag_summaries$se_scale == se_scale &
        replicate_lag_summaries$estimator == "direct_zero_null",
    ]
    for (estimator in c("within_unit_centered", "pairwise_ols")) {
      comparison <- replicate_lag_summaries[
        replicate_lag_summaries$threshold == threshold &
          replicate_lag_summaries$se_scale == se_scale &
          replicate_lag_summaries$estimator == estimator,
      ]
      comparison <- comparison[match(
        paste(reference$replication, reference$lag),
        paste(comparison$replication, comparison$lag)
      ), ]
      difference <- comparison$covariance - reference$covariance
      contrast_data <- cbind(
        comparison[c(
          "replication", "threshold", "n_units", "se_scale", "lag"
        )],
        data.frame(
          estimator = estimator,
          reference = "direct_zero_null",
          covariance_difference = difference,
          stringsAsFactors = FALSE
        )
      )
      contrast_rows[[contrast_index]] <- contrast_data
      contrast_index <- contrast_index + 1L
    }
  }
}
replicate_method_contrasts <- do.call(rbind, contrast_rows)
rownames(replicate_method_contrasts) <- NULL
aggregate_method_contrasts <- aggregate_values(
  replicate_method_contrasts,
  "covariance_difference",
  c(
    "threshold", "n_units", "se_scale", "estimator", "reference", "lag"
  )
)

mean_matrix_rows <- list()
diagnostic_rows <- list()
mean_matrix_index <- 1L
diagnostic_index <- 1L
for (grid_index in seq_len(nrow(key_grid))) {
  key <- keys[grid_index]
  threshold_index <- match(
    key_grid$threshold_name[grid_index],
    threshold_names
  )
  mean_matrix <- apply(
    covariance_signal_draws[[key]],
    c(1L, 2L),
    mean
  )
  mean_matrix_rows[[mean_matrix_index]] <- data.frame(
    threshold = unname(threshold_values[threshold_index]),
    n_units = unname(threshold_counts[threshold_index]),
    se_scale = key_grid$se_scale[grid_index],
    estimator = key_grid$estimator[grid_index],
    time_a = rep(time_grid, times = n_time),
    time_b = rep(time_grid, each = n_time),
    covariance = as.vector(mean_matrix),
    stringsAsFactors = FALSE
  )
  mean_matrix_index <- mean_matrix_index + 1L
  off_diagonal <- mean_matrix[upper.tri(mean_matrix)]
  eigenvalues <- eigen(mean_matrix, symmetric = TRUE, only.values = TRUE)$values
  diagnostic_rows[[diagnostic_index]] <- data.frame(
    threshold = unname(threshold_values[threshold_index]),
    n_units = unname(threshold_counts[threshold_index]),
    se_scale = key_grid$se_scale[grid_index],
    estimator = key_grid$estimator[grid_index],
    mean_diagonal = mean(diag(mean_matrix)),
    minimum_diagonal = min(diag(mean_matrix)),
    maximum_diagonal = max(diag(mean_matrix)),
    mean_off_diagonal = mean(off_diagonal),
    minimum_off_diagonal = min(off_diagonal),
    maximum_off_diagonal = max(off_diagonal),
    minimum_eigenvalue = min(eigenvalues),
    leading_eigenvalue_fraction = max(eigenvalues) / sum(eigenvalues),
    stringsAsFactors = FALSE
  )
  diagnostic_index <- diagnostic_index + 1L
}
mean_covariance_matrices_long <- do.call(rbind, mean_matrix_rows)
mean_matrix_diagnostics <- do.call(rbind, diagnostic_rows)
rownames(mean_covariance_matrices_long) <- NULL
rownames(mean_matrix_diagnostics) <- NULL

agreement_rows <- list()
agreement_index <- 1L
for (threshold in threshold_values) {
  for (se_scale in se_scales) {
    direct_key <- matrix_key(
      threshold_names[match(threshold, threshold_values)],
      se_scale,
      "direct_zero_null"
    )
    direct_mean <- apply(
      covariance_signal_draws[[direct_key]],
      c(1L, 2L),
      mean
    )
    direct_off <- direct_mean[upper.tri(direct_mean)]
    for (estimator in c("within_unit_centered", "pairwise_ols")) {
      comparison_key <- matrix_key(
        threshold_names[match(threshold, threshold_values)],
        se_scale,
        estimator
      )
      comparison_mean <- apply(
        covariance_signal_draws[[comparison_key]],
        c(1L, 2L),
        mean
      )
      comparison_off <- comparison_mean[upper.tri(comparison_mean)]
      agreement_rows[[agreement_index]] <- data.frame(
        threshold = threshold,
        n_units = unname(threshold_counts[match(threshold, threshold_values)]),
        se_scale = se_scale,
        estimator = estimator,
        reference = "direct_zero_null",
        off_diagonal_rmse = sqrt(mean((comparison_off - direct_off)^2)),
        off_diagonal_mean_difference = mean(comparison_off - direct_off),
        off_diagonal_maximum_absolute_difference =
          max(abs(comparison_off - direct_off)),
        off_diagonal_matrix_correlation = stats::cor(
          comparison_off,
          direct_off
        ),
        stringsAsFactors = FALSE
      )
      agreement_index <- agreement_index + 1L
    }
  }
}
matrix_agreement <- do.call(rbind, agreement_rows)
rownames(matrix_agreement) <- NULL

configuration <- list(
  analysis_id = "permuted_variogram_covariance_summary",
  created_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  synchronized_cache = comparison_path,
  input_cache = input_path,
  thresholds = threshold_values,
  selected_counts = threshold_counts,
  primary_threshold = 0.95,
  primary_se_scale = "t_adjusted",
  n_replications = n_replications,
  synchronized_seeds = synchronized$configuration$seeds,
  independent_time_seeds = independent_seeds,
  synchronized_generator = paste(
    "One donor residual-block mapping shared across all 16 time points,",
    "preserving cross-time donor pairing"
  ),
  independent_generator = paste(
    "A separate donor residual-row permutation at each time point,",
    "preserving marginal residual distributions while breaking cross-time pairing"
  ),
  covariance_estimand = paste(
    "Method-matched independent-time variogram minus synchronized variogram;",
    "shared diagonal from the direct synchronized zero-null covariance"
  ),
  estimators = estimator_labels,
  r_version = R.version.string
)

result <- list(
  configuration = configuration,
  replicate_lag_summaries = replicate_lag_summaries,
  aggregate_lag_summaries = aggregate_lag_summaries,
  replicate_method_contrasts = replicate_method_contrasts,
  aggregate_method_contrasts = aggregate_method_contrasts,
  mean_covariance_matrices_long = mean_covariance_matrices_long,
  mean_matrix_diagnostics = mean_matrix_diagnostics,
  matrix_agreement = matrix_agreement,
  independent_covariance_draws = independent_covariance_draws,
  synchronized_variogram_draws = synchronized_variogram_draws,
  independent_variogram_draws = independent_variogram_draws,
  covariance_signal_draws = covariance_signal_draws
)
saveRDS(configuration, file.path(output_dir, "configuration.rds"))
saveRDS(
  result,
  file.path(output_dir, "permuted_variogram_covariance_summary.rds"),
  compress = "xz"
)
utils::write.csv(
  replicate_lag_summaries,
  file.path(summary_dir, "replicate_lag_summaries.csv"),
  row.names = FALSE
)
utils::write.csv(
  aggregate_lag_summaries,
  file.path(summary_dir, "aggregate_lag_summaries.csv"),
  row.names = FALSE
)
utils::write.csv(
  aggregate_method_contrasts,
  file.path(summary_dir, "aggregate_method_contrasts.csv"),
  row.names = FALSE
)
utils::write.csv(
  mean_covariance_matrices_long,
  file.path(summary_dir, "mean_covariance_matrices_long.csv"),
  row.names = FALSE
)
utils::write.csv(
  mean_matrix_diagnostics,
  file.path(summary_dir, "mean_matrix_diagnostics.csv"),
  row.names = FALSE
)
utils::write.csv(
  matrix_agreement,
  file.path(summary_dir, "matrix_agreement.csv"),
  row.names = FALSE
)

primary_key <- matrix_key("lfdr_0p95", "t_adjusted", "direct_zero_null")
primary_matrix <- apply(
  covariance_signal_draws[[primary_key]],
  c(1L, 2L),
  mean
)
rownames(primary_matrix) <- colnames(primary_matrix) <- paste0("time_", time_grid)
utils::write.csv(
  primary_matrix,
  file.path(summary_dir, "primary_direct_empirical_covariance_matrix.csv")
)

colors <- c(
  pairwise_ols = "#0072B2",
  within_unit_centered = "#D55E00",
  direct_zero_null = "#009E73"
)

plot_variogram_covariance <- function(path) {
  plot_data <- aggregate_lag_summaries[
    aggregate_lag_summaries$threshold == 0.95 &
      aggregate_lag_summaries$se_scale == "t_adjusted",
  ]
  synchronized_data <- plot_data[
    plot_data$metric == "synchronized_variogram",
  ]
  independent_data <- plot_data[
    plot_data$metric == "independent_variogram",
  ]
  covariance_data <- plot_data[plot_data$metric == "covariance", ]
  grDevices::png(path, width = 1400, height = 1200, res = 160)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::layout(matrix(1:2, ncol = 1L))
  graphics::par(mar = c(3.0, 5.1, 3.4, 1.2))
  top_range <- range(c(
    synchronized_data$lower,
    synchronized_data$upper,
    independent_data$lower,
    independent_data$upper
  ))
  top_range[2L] <- max(top_range[2L], 1.22)
  graphics::plot(
    NA,
    xlim = c(1, 15),
    ylim = top_range,
    xlab = "",
    ylab = "Covariance-scale semivariogram",
    main = "Synchronized variogram and matched independence benchmark",
    xaxt = "n"
  )
  graphics::axis(1, at = 1:15, labels = FALSE)
  for (estimator in estimator_names) {
    synchronized_estimator <- synchronized_data[
      synchronized_data$estimator == estimator,
    ]
    synchronized_estimator <- synchronized_estimator[
      order(synchronized_estimator$lag),
    ]
    independent_estimator <- independent_data[
      independent_data$estimator == estimator,
    ]
    independent_estimator <- independent_estimator[
      order(independent_estimator$lag),
    ]
    graphics::polygon(
      c(
        synchronized_estimator$lag,
        rev(synchronized_estimator$lag)
      ),
      c(
        synchronized_estimator$lower,
        rev(synchronized_estimator$upper)
      ),
      col = grDevices::adjustcolor(colors[estimator], alpha.f = 0.12),
      border = NA
    )
    graphics::lines(
      synchronized_estimator$lag,
      synchronized_estimator$mean,
      col = colors[estimator],
      lwd = 2.5
    )
    graphics::lines(
      independent_estimator$lag,
      independent_estimator$mean,
      col = colors[estimator],
      lwd = 2.1,
      lty = 2
    )
  }
  graphics::legend(
    "topleft",
    legend = c(
      estimator_labels[estimator_names],
      "Solid: synchronized",
      "Dashed: independent-time benchmark"
    ),
    col = c(colors[estimator_names], "black", "black"),
    lwd = c(rep(2.5, 3L), 2.2, 2.2),
    lty = c(rep(1, 4L), 2),
    bty = "n",
    cex = 0.82
  )

  graphics::par(mar = c(4.8, 5.1, 3.4, 1.2))
  bottom_range <- range(c(covariance_data$lower, covariance_data$upper, 0))
  graphics::plot(
    NA,
    xlim = c(1, 15),
    ylim = bottom_range,
    xlab = "Lag",
    ylab = "Estimated covariance",
    main = "Independence benchmark minus synchronized variogram",
    xaxt = "n"
  )
  graphics::axis(1, at = 1:15)
  graphics::abline(h = 0, lty = 3, col = "grey45")
  for (estimator in estimator_names) {
    estimator_data <- covariance_data[covariance_data$estimator == estimator, ]
    estimator_data <- estimator_data[order(estimator_data$lag), ]
    graphics::polygon(
      c(estimator_data$lag, rev(estimator_data$lag)),
      c(estimator_data$lower, rev(estimator_data$upper)),
      col = grDevices::adjustcolor(colors[estimator], alpha.f = 0.14),
      border = NA
    )
    graphics::lines(
      estimator_data$lag,
      estimator_data$mean,
      col = colors[estimator],
      lwd = 2.5
    )
  }
  graphics::legend(
    "topright",
    legend = estimator_labels[estimator_names],
    col = colors[estimator_names],
    lwd = 2.5,
    bty = "n"
  )
}

get_mean_covariance <- function(estimator) {
  key <- matrix_key("lfdr_0p95", "t_adjusted", estimator)
  apply(covariance_signal_draws[[key]], c(1L, 2L), mean)
}

plot_covariance_matrices <- function(path) {
  matrices <- lapply(estimator_names, get_mean_covariance)
  names(matrices) <- estimator_names
  z_range <- range(c(0, unlist(matrices)))
  palette <- grDevices::hcl.colors(201L, "Reds 3", rev = TRUE)
  grDevices::png(path, width = 2250, height = 700, res = 150)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::layout(matrix(1:4, nrow = 1L), widths = c(1, 1, 1, 0.20))
  graphics::par(oma = c(0, 0, 2.7, 0))
  for (estimator in estimator_names) {
    graphics::par(mar = c(4.2, 4.2, 3.2, 1.0))
    graphics::image(
      x = time_grid,
      y = time_grid,
      z = matrices[[estimator]],
      col = palette,
      zlim = z_range,
      xlab = "Time",
      ylab = "Time",
      main = estimator_labels[estimator],
      axes = FALSE,
      useRaster = TRUE
    )
    graphics::axis(1, at = c(0, 5, 10, 15))
    graphics::axis(2, at = c(0, 5, 10, 15))
    graphics::box()
  }
  graphics::par(mar = c(4.2, 0.2, 3.2, 2.5))
  color_values <- seq(z_range[1L], z_range[2L], length.out = length(palette))
  graphics::image(
    x = c(0, 1),
    y = color_values,
    z = matrix(rep(color_values, each = 2L), nrow = 2L),
    col = palette,
    xlab = "",
    ylab = "",
    axes = FALSE,
    useRaster = TRUE
  )
  graphics::axis(4, las = 1)
  graphics::box()
  graphics::mtext(
    paste(
      "Mean empirical covariance matrices from matched permutations;",
      "lfdr > 0.95, t-adjusted SE"
    ),
    outer = TRUE,
    line = 0.7,
    cex = 1.08
  )
}

plot_annotated_primary_matrix <- function(path) {
  palette <- grDevices::colorRampPalette(
    c("white", "#FEE8C8", "#F46D43", "#A50026")
  )(201L)
  grDevices::png(path, width = 1500, height = 1350, res = 160)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mar = c(5.0, 5.0, 4.2, 1.5))
  graphics::image(
    x = time_grid,
    y = time_grid,
    z = primary_matrix,
    col = palette,
    zlim = c(0, max(primary_matrix)),
    xlab = "Time",
    ylab = "Time",
    main = paste(
      "Direct zero-null empirical covariance matrix\n",
      "matched donor permutations; lfdr > 0.95, t-adjusted SE"
    ),
    axes = FALSE,
    useRaster = TRUE
  )
  graphics::axis(1, at = time_grid)
  graphics::axis(2, at = time_grid)
  graphics::box()
  for (row_index in seq_len(n_time)) {
    for (column_index in seq_len(n_time)) {
      value <- primary_matrix[column_index, row_index]
      text_color <- if (value > 0.55 * max(primary_matrix)) "white" else "black"
      graphics::text(
        time_grid[row_index],
        time_grid[column_index],
        labels = formatC(value, format = "f", digits = 2L),
        cex = 0.66,
        col = text_color
      )
    }
  }
}

plot_variogram_covariance(file.path(
  figure_dir,
  "matched_variogram_and_covariance_t_adjusted.png"
))
plot_covariance_matrices(file.path(
  figure_dir,
  "matched_empirical_covariance_matrices_t_adjusted.png"
))
plot_annotated_primary_matrix(file.path(
  figure_dir,
  "direct_empirical_covariance_matrix_t_adjusted.png"
))

primary_output <- aggregate_lag_summaries[
  aggregate_lag_summaries$threshold == 0.95 &
    aggregate_lag_summaries$se_scale == "t_adjusted" &
    aggregate_lag_summaries$metric == "covariance" &
    aggregate_lag_summaries$lag %in% c(1L, 9L, 15L),
]
print(primary_output, row.names = FALSE)
cat("Permuted variogram/covariance summary completed.\n")
