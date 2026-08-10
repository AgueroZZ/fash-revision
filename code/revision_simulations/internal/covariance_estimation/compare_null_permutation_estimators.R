#!/usr/bin/env Rscript

# Compare residual-correlation estimators under synchronized donor permutations.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
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

lag_average_semivariance <- function(covariance) {
  n_time <- ncol(covariance)
  diagonal <- diag(covariance)
  vapply(seq_len(n_time - 1L), function(lag) {
    index_a <- seq_len(n_time - lag)
    index_b <- (lag + 1L):n_time
    mean(
      0.5 * (diagonal[index_a] + diagonal[index_b]) -
        covariance[cbind(index_a, index_b)]
    )
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
input_path <- file.path(pilot_dir, "input", "selected_raw_data.rds")
pilot_result_path <- file.path(
  pilot_dir,
  "donor_null_permutation_analysis.rds"
)
required_inputs <- c(input_path, pilot_result_path)
if (any(!file.exists(required_inputs))) {
  stop("The completed donor-permutation pilot inputs are unavailable.")
}

input <- readRDS(input_path)
pilot_result <- readRDS(pilot_result_path)
n_replications <- pilot_result$configuration$n_replications
replication_seeds <- pilot_result$configuration$residual_permutation_seeds
if (n_replications != 100L || length(replication_seeds) != n_replications) {
  stop("The expected 100-replication donor-permutation cache was not found.")
}

output_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "null_permutation_estimator_comparison"
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
correlation_draws <- setNames(lapply(keys, function(unused) {
  make_draw_array()
}), keys)
covariance_draws <- setNames(lapply(keys, function(unused) {
  make_draw_array()
}), keys)

fit_residual_permutation <- function(seed) {
  set.seed(seed)
  donor_mapping <- make_shared_donor_block_permutation(
    input$vcf_donors,
    input$donor_observation_patterns
  )
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
    source_donors <- unname(donor_mapping[time_input$donors])
    source_rows <- match(source_donors, time_input$donors)
    if (anyNA(source_rows)) {
      stop("A residual donor block is unavailable at one time point.")
    }
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

lag_rows <- vector(
  "list",
  n_replications * length(threshold_names) *
    length(se_scales) * length(estimator_names)
)
diagnostic_rows <- lag_rows
output_index <- 1L
message("Running paired estimator comparisons on 100 donor permutations.")
for (replication in seq_len(n_replications)) {
  fit <- fit_residual_permutation(replication_seeds[replication])
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
      pairwise_correlation <- estimate_pairwise_difference_correlation(
        selected_beta,
        selected_se
      )
      estimator_results <- list(
        pairwise_ols = list(
          covariance = pairwise_correlation,
          correlation = pairwise_correlation
        ),
        within_unit_centered = centered,
        direct_zero_null = direct
      )

      for (estimator in estimator_names) {
        result <- estimator_results[[estimator]]
        key <- matrix_key(
          threshold_names[threshold_index],
          se_scale,
          estimator
        )
        correlation_draws[[key]][, , replication] <- result$correlation
        covariance_draws[[key]][, , replication] <- result$covariance
        lag_correlation <- lag_average_correlation(result$correlation)
        lag_semivariance <- lag_average_semivariance(result$covariance)
        lag_rows[[output_index]] <- data.frame(
          replication = replication,
          seed = replication_seeds[replication],
          threshold = threshold_values[threshold_index],
          n_units = length(selected_index),
          se_scale = se_scale,
          estimator = estimator,
          lag = seq_along(lag_correlation),
          correlation = lag_correlation,
          standardized_variogram = 1 - lag_correlation,
          covariance_semivariogram = lag_semivariance,
          stringsAsFactors = FALSE
        )
        correlation_summary <- summarize_raw_correlation_matrix(
          result$correlation
        )
        diagnostic_rows[[output_index]] <- cbind(
          data.frame(
            replication = replication,
            seed = replication_seeds[replication],
            threshold = threshold_values[threshold_index],
            n_units = length(selected_index),
            se_scale = se_scale,
            estimator = estimator,
            mean_standardized_variance = mean(diag(result$covariance)),
            lag1_minus_long = correlation_summary$lag1 -
              correlation_summary$mean_lags_9_15,
            stringsAsFactors = FALSE
          ),
          correlation_summary
        )
        output_index <- output_index + 1L
      }
    }
  }
  if (replication %% 10L == 0L) {
    message("Completed ", replication, " of ", n_replications, ".")
  }
}

replicate_lag_summaries <- do.call(rbind, lag_rows)
replicate_diagnostics <- do.call(rbind, diagnostic_rows)
rownames(replicate_lag_summaries) <- NULL
rownames(replicate_diagnostics) <- NULL

lag_grouping <- c(
  "threshold", "n_units", "se_scale", "estimator", "lag"
)
aggregate_lag_summaries <- do.call(rbind, lapply(
  c("correlation", "standardized_variogram", "covariance_semivariogram"),
  function(metric) {
    aggregate_values(replicate_lag_summaries, metric, lag_grouping)
  }
))
rownames(aggregate_lag_summaries) <- NULL

diagnostic_grouping <- c(
  "threshold", "n_units", "se_scale", "estimator"
)
diagnostic_metrics <- c(
  "mean_standardized_variance", "lag1", "mean_lags_9_15", "lag15",
  "mean_off_diagonal", "lag1_minus_long"
)
aggregate_diagnostics <- do.call(rbind, lapply(
  diagnostic_metrics,
  function(metric) {
    aggregate_values(replicate_diagnostics, metric, diagnostic_grouping)
  }
))
rownames(aggregate_diagnostics) <- NULL

mean_matrix_rows <- list()
mean_matrix_index <- 1L
for (grid_index in seq_len(nrow(key_grid))) {
  key <- keys[grid_index]
  for (matrix_type in c("covariance", "correlation")) {
    draws <- if (matrix_type == "covariance") {
      covariance_draws[[key]]
    } else {
      correlation_draws[[key]]
    }
    mean_matrix <- apply(draws, c(1L, 2L), mean)
    mean_matrix_rows[[mean_matrix_index]] <- data.frame(
      threshold = unname(threshold_values[
        match(key_grid$threshold_name[grid_index], threshold_names)
      ]),
      n_units = unname(threshold_counts[
        match(key_grid$threshold_name[grid_index], threshold_names)
      ]),
      se_scale = key_grid$se_scale[grid_index],
      estimator = key_grid$estimator[grid_index],
      matrix_type = matrix_type,
      time_a = rep(time_grid, times = n_time),
      time_b = rep(time_grid, each = n_time),
      value = as.vector(mean_matrix),
      stringsAsFactors = FALSE
    )
    mean_matrix_index <- mean_matrix_index + 1L
  }
}
mean_matrices_long <- do.call(rbind, mean_matrix_rows)
rownames(mean_matrices_long) <- NULL

paired_rows <- list()
paired_index <- 1L
paired_metrics <- c(
  "lag1", "mean_lags_9_15", "lag15", "mean_off_diagonal",
  "lag1_minus_long"
)
for (threshold in threshold_values) {
  for (se_scale in se_scales) {
    subset_index <- replicate_diagnostics$threshold == threshold &
      replicate_diagnostics$se_scale == se_scale
    subset_data <- replicate_diagnostics[subset_index, ]
    pairwise <- subset_data[
      subset_data$estimator == "pairwise_ols",
      c("replication", paired_metrics)
    ]
    for (comparison_estimator in c(
      "within_unit_centered", "direct_zero_null"
    )) {
      comparison <- subset_data[
        subset_data$estimator == comparison_estimator,
        c("replication", paired_metrics)
      ]
      comparison <- comparison[match(pairwise$replication,
                                     comparison$replication), ]
      for (metric in paired_metrics) {
        difference <- pairwise[[metric]] - comparison[[metric]]
        paired_rows[[paired_index]] <- data.frame(
          threshold = threshold,
          se_scale = se_scale,
          contrast = paste0("pairwise_minus_", comparison_estimator),
          metric = metric,
          mean = mean(difference),
          median = stats::median(difference),
          lower = unname(stats::quantile(difference, 0.025, names = FALSE)),
          upper = unname(stats::quantile(difference, 0.975, names = FALSE)),
          sd = stats::sd(difference),
          n_replications = length(difference),
          stringsAsFactors = FALSE
        )
        paired_index <- paired_index + 1L
      }
    }
  }
}
paired_estimator_differences <- do.call(rbind, paired_rows)
rownames(paired_estimator_differences) <- NULL

pairwise_reproduction_errors <- numeric()
for (threshold_index in seq_along(threshold_names)) {
  for (se_scale in se_scales) {
    new_key <- matrix_key(
      threshold_names[threshold_index],
      se_scale,
      "pairwise_ols"
    )
    old_key <- paste(
      "donor_residual_block_permutation",
      threshold_names[threshold_index],
      se_scale,
      "ols_weighted",
      sep = "__"
    )
    pairwise_reproduction_errors <- c(
      pairwise_reproduction_errors,
      max(abs(
        correlation_draws[[new_key]] -
          pilot_result$matrix_draws[[old_key]]
      ))
    )
  }
}
maximum_pairwise_reproduction_error <- max(pairwise_reproduction_errors)
if (maximum_pairwise_reproduction_error > 1e-12) {
  stop("The paired rerun did not reproduce the existing pairwise matrices.")
}

configuration <- list(
  analysis_id = "null_permutation_estimator_comparison",
  created_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  input_cache = input_path,
  pilot_result = pilot_result_path,
  thresholds = threshold_values,
  selected_counts = threshold_counts,
  n_replications = n_replications,
  seeds = replication_seeds,
  estimators = c(
    pairwise_ols = "OLS-through-origin pairwise-difference estimator",
    within_unit_centered = paste(
      "Inverse-variance constant estimate removed within each unit,",
      "followed by empirical correlation across units"
    ),
    direct_zero_null = paste(
      "Known-zero-null standardized second moment normalized to correlation;",
      "no within-unit centering"
    )
  ),
  matrix_types = c(
    "Standardized empirical covariance",
    "Corresponding correlation matrix"
  ),
  primary_variogram = paste(
    "Covariance-scale semivariogram:",
    "0.5 * (V_tt + V_ss - 2 * V_ts), averaged within lag"
  ),
  maximum_pairwise_reproduction_error =
    maximum_pairwise_reproduction_error,
  r_version = R.version.string
)

result <- list(
  configuration = configuration,
  replicate_lag_summaries = replicate_lag_summaries,
  replicate_diagnostics = replicate_diagnostics,
  aggregate_lag_summaries = aggregate_lag_summaries,
  aggregate_diagnostics = aggregate_diagnostics,
  paired_estimator_differences = paired_estimator_differences,
  mean_matrices_long = mean_matrices_long,
  correlation_draws = correlation_draws,
  covariance_draws = covariance_draws
)
saveRDS(configuration, file.path(output_dir, "configuration.rds"))
saveRDS(
  result,
  file.path(output_dir, "estimator_comparison.rds"),
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
  aggregate_diagnostics,
  file.path(summary_dir, "aggregate_diagnostics.csv"),
  row.names = FALSE
)
utils::write.csv(
  paired_estimator_differences,
  file.path(summary_dir, "paired_estimator_differences.csv"),
  row.names = FALSE
)
utils::write.csv(
  mean_matrices_long,
  file.path(summary_dir, "mean_matrices_long.csv"),
  row.names = FALSE
)

plot_mean_matrices <- function(se_scale, matrix_type, path) {
  estimator_order <- estimator_names
  display_se <- if (se_scale == "raw_regression") {
    "raw regression"
  } else {
    "t-adjusted"
  }
  matrices <- lapply(estimator_order, function(estimator) {
    key <- matrix_key("lfdr_0p95", se_scale, estimator)
    draws <- if (matrix_type == "covariance") {
      covariance_draws[[key]]
    } else {
      correlation_draws[[key]]
    }
    apply(draws, c(1L, 2L), mean)
  })
  z_limit <- range(unlist(matrices))
  grDevices::png(path, width = 2250, height = 680, res = 150)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::layout(
    matrix(1:4, nrow = 1L),
    widths = c(1, 1, 1, 0.12)
  )
  graphics::par(oma = c(0, 0, 2.5, 0))
  palette <- grDevices::hcl.colors(101L, "Blue-Red 3", rev = TRUE)
  for (estimator_index in seq_along(estimator_order)) {
    graphics::par(mar = c(4.2, 4.2, 3.2, 1.2))
    graphics::image(
      x = time_grid,
      y = time_grid,
      z = matrices[[estimator_index]],
      col = palette,
      zlim = z_limit,
      xlab = "Time",
      ylab = "Time",
      main = estimator_labels[estimator_order[estimator_index]],
      axes = FALSE,
      useRaster = TRUE
    )
    graphics::axis(1, at = c(0, 5, 10, 15))
    graphics::axis(2, at = c(0, 5, 10, 15))
    graphics::box()
  }
  graphics::par(mar = c(4.2, 0.5, 3.2, 3.8))
  color_values <- seq(z_limit[1L], z_limit[2L], length.out = length(palette))
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
      "Mean",
      matrix_type,
      "matrices across 100 synchronized donor permutations;",
      display_se,
      "SE"
    ),
    outer = TRUE,
    line = 0.7,
    cex = 1.05
  )
}

plot_covariance_variogram <- function(se_scale, path) {
  plot_data <- aggregate_lag_summaries[
    aggregate_lag_summaries$threshold == 0.95 &
      aggregate_lag_summaries$se_scale == se_scale &
      aggregate_lag_summaries$metric == "covariance_semivariogram",
  ]
  colors <- c(
    pairwise_ols = "#0072B2",
    within_unit_centered = "#D55E00",
    direct_zero_null = "#009E73"
  )
  y_range <- range(c(plot_data$lower, plot_data$upper))
  grDevices::png(path, width = 1200, height = 800, res = 150)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mar = c(4.7, 5.0, 3.4, 1.2))
  graphics::plot(
    NA,
    xlim = c(1, 15),
    ylim = y_range,
    xlab = "Lag",
    ylab = expression(
      "Semivariogram  " * frac(1, 2) *
        plain(Var) * "{" * Z(t + h) - Z(t) * "}"
    ),
    main = paste0(
      "Covariance-scale variogram under synchronized donor permutations\n",
      if (se_scale == "raw_regression") {
        "Raw regression SE"
      } else {
        "t-adjusted SE"
      }
    ),
    xaxt = "n"
  )
  graphics::axis(1, at = 1:15)
  for (estimator in estimator_names) {
    estimator_data <- plot_data[plot_data$estimator == estimator, ]
    estimator_data <- estimator_data[order(estimator_data$lag), ]
    fill_color <- grDevices::adjustcolor(colors[estimator], alpha.f = 0.16)
    graphics::polygon(
      c(estimator_data$lag, rev(estimator_data$lag)),
      c(estimator_data$lower, rev(estimator_data$upper)),
      col = fill_color,
      border = NA
    )
    graphics::lines(
      estimator_data$lag,
      estimator_data$mean,
      col = colors[estimator],
      lwd = 2.4
    )
  }
  graphics::legend(
    "bottomright",
    legend = estimator_labels[estimator_names],
    col = colors[estimator_names],
    lwd = 2.4,
    bty = "n"
  )
}

for (se_scale in se_scales) {
  for (matrix_type in c("covariance", "correlation")) {
    plot_mean_matrices(
      se_scale,
      matrix_type,
      file.path(
        figure_dir,
        paste0("mean_", matrix_type, "_matrices_", se_scale, ".png")
      )
    )
  }
  plot_covariance_variogram(
    se_scale,
    file.path(
      figure_dir,
      paste0("covariance_variogram_", se_scale, ".png")
    )
  )
}

key_output <- aggregate_diagnostics[
  aggregate_diagnostics$threshold == 0.95 &
    aggregate_diagnostics$metric %in%
      c("lag1", "mean_lags_9_15", "lag15", "lag1_minus_long"),
]
print(key_output, row.names = FALSE)
cat(
  "Maximum pairwise reproduction error: ",
  format(maximum_pairwise_reproduction_error, scientific = TRUE),
  "\n",
  sep = ""
)
cat("Null-permutation estimator comparison completed.\n")
