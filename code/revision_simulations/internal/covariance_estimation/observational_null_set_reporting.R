# Reporting helpers for the observational null-set covariance pilot.

pilot_cache_path <- file.path(
  "output", "revision_simulations", "internal",
  "observational_null_set_pilot", "observational_null_set_pilot.rds"
)
if (!file.exists(pilot_cache_path)) {
  stop("The observational null-set pilot cache is missing: ", pilot_cache_path)
}
pilot_result <- readRDS(pilot_cache_path)
required_components <- c(
  "configuration", "selection_counts", "selected_units",
  "time_specific_screen_summary", "discovery_summary", "selected_data",
  "matrices", "centered_correlations",
  "beta_scale_matrices", "matrix_diagnostics",
  "covariance_matrices_long", "lag_summaries",
  "bootstrap_lag_intervals"
)
if (!all(required_components %in% names(pilot_result))) {
  stop("The observational null-set pilot cache is incomplete.")
}

pilot_configuration <- pilot_result$configuration
selection_counts <- pilot_result$selection_counts
selected_units <- pilot_result$selected_units
time_specific_screen_summary <- pilot_result$time_specific_screen_summary
discovery_summary <- pilot_result$discovery_summary
matrix_diagnostics <- pilot_result$matrix_diagnostics
covariance_matrices_long <- pilot_result$covariance_matrices_long
lag_summaries <- pilot_result$lag_summaries
bootstrap_lag_intervals <- pilot_result$bootstrap_lag_intervals
selected_data <- pilot_result$selected_data
centered_correlations <- pilot_result$centered_correlations

if (!identical(pilot_configuration$status, "Internal small-scale pilot") ||
    pilot_configuration$pilot_size != 200L ||
    pilot_configuration$n_bootstrap != 200L ||
    nrow(selected_units) != 400L ||
    any(selected_units$variant_gene_count != 1L) ||
    any(selected_units$gene_has_dynamic_discovery) ||
    any(selected_units$set_id == "B_zero_eqtl_null" &
          selected_units$gene_has_time_specific_discovery)) {
  stop("The cached pilot does not satisfy its recorded design invariants.")
}

set_labels <- c(
  A_dynamic_null = "A: no FASH dynamic discovery",
  B_zero_eqtl_null = "B: no dynamic or time-specific discovery"
)
estimator_labels <- c(
  `Pairwise difference` = "Pairwise difference",
  `Within-unit centered covariance` = "Within-unit centered covariance",
  `Direct zero mean` = "Direct zero mean"
)

format_integer <- function(x) {
  format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

format_decimal <- function(x, digits = 3L) {
  formatC(as.numeric(x), format = "f", digits = digits)
}

render_scrollable_table <- function(x,
                                    align = NULL,
                                    minimum_width = "720px",
                                    digits = 3L) {
  table_html <- knitr::kable(
    x,
    format = "html",
    escape = TRUE,
    align = align,
    digits = digits
  )
  cat(
    paste0(
      '<div class="null-pilot-table-scroll" style="min-width:',
      minimum_width,
      ';">',
      table_html,
      "</div>"
    )
  )
}

selection_table <- selection_counts
selection_table$set_id <- unname(set_labels[selection_table$set_id])
names(selection_table) <- c(
  "Set", "Definition", "Eligible genes",
  "Eligible genes with a gene-unique variant", "Selected pairs",
  "Candidate pairs"
)

discovery_table <- data.frame(
  Quantity = c(
    "All tested FASH pairs",
    "All tested genes",
    "FASH dynamic-discovered pairs",
    "FASH dynamic-discovered genes",
    "Largest selected dynamic lfdr",
    "Genes discovered at one or more individual days"
  ),
  Value = c(
    format_integer(pilot_configuration$n_fash_pairs),
    format_integer(pilot_configuration$n_genes),
    format_integer(discovery_summary$n_dynamic_discovered_pairs),
    format_integer(discovery_summary$n_dynamic_discovered_genes),
    format_decimal(discovery_summary$dynamic_lfdr_cutoff, 4L),
    format_integer(discovery_summary$n_time_specific_discovered_genes_union)
  ),
  stringsAsFactors = FALSE
)

selection_lfdr_summary <- do.call(rbind, lapply(
  split(selected_units$bf_adjusted_lfdr, selected_units$set_id),
  function(values) {
    data.frame(
      Minimum = min(values),
      `First quartile` = stats::quantile(values, 0.25, names = FALSE),
      Median = stats::median(values),
      Mean = mean(values),
      `Third quartile` = stats::quantile(values, 0.75, names = FALSE),
      Maximum = max(values),
      check.names = FALSE
    )
  }
))
selection_lfdr_summary$Set <- unname(set_labels[rownames(selection_lfdr_summary)])
rownames(selection_lfdr_summary) <- NULL
selection_lfdr_summary <- selection_lfdr_summary[, c(
  "Set", "Minimum", "First quartile", "Median", "Mean",
  "Third quartile", "Maximum"
)]

diagnostic_table <- matrix_diagnostics
diagnostic_table$set_id <- unname(set_labels[diagnostic_table$set_id])
diagnostic_table$estimator <- unname(
  estimator_labels[diagnostic_table$estimator]
)
names(diagnostic_table) <- c(
  "Set", "Estimator", "Minimum eigenvalue", "Maximum eigenvalue",
  "Negative eigenvalues", "Minimum diagonal", "Maximum diagonal",
  "Mean off-diagonal", "Minimum off-diagonal", "Maximum off-diagonal"
)

key_lags <- c(1L, 5L, 10L, 15L)
key_lag_table <- bootstrap_lag_intervals[
  bootstrap_lag_intervals$lag %in% key_lags,
  c(
    "set_id", "estimator", "lag", "observed_covariance",
    "covariance_ci_lower", "covariance_ci_upper",
    "observed_semivariogram", "semivariogram_ci_lower",
    "semivariogram_ci_upper"
  )
]
key_lag_table$set_id <- unname(set_labels[key_lag_table$set_id])
key_lag_table$estimator <- unname(
  estimator_labels[key_lag_table$estimator]
)
names(key_lag_table) <- c(
  "Set", "Estimator", "Lag", "Covariance",
  "Covariance 2.5%", "Covariance 97.5%", "Semivariogram",
  "Semivariogram 2.5%", "Semivariogram 97.5%"
)

pairwise_b <- pilot_result$matrices$B_zero_eqtl_null__pairwise
direct_b <- pilot_result$matrices$B_zero_eqtl_null__direct
off_diagonal_index <- upper.tri(pairwise_b)
method_agreement <- data.frame(
  Quantity = c(
    "Off-diagonal Pearson correlation",
    "Off-diagonal mean absolute difference",
    "Off-diagonal root mean squared difference",
    "Off-diagonal maximum absolute difference",
    "Maximum absolute difference in lag-averaged covariance"
  ),
  Value = c(
    stats::cor(pairwise_b[off_diagonal_index], direct_b[off_diagonal_index]),
    mean(abs(pairwise_b[off_diagonal_index] - direct_b[off_diagonal_index])),
    sqrt(mean((pairwise_b[off_diagonal_index] - direct_b[off_diagonal_index])^2)),
    max(abs(pairwise_b[off_diagonal_index] - direct_b[off_diagonal_index])),
    max(abs(
      lag_average_correlation(pairwise_b) -
        lag_average_correlation(direct_b)
    ))
  ),
  stringsAsFactors = FALSE
)

plot_standardized_covariance_heatmaps <- function() {
  matrix_ids <- c(
    "A_dynamic_null__pairwise",
    "B_zero_eqtl_null__pairwise",
    "B_zero_eqtl_null__direct"
  )
  plot_data <- covariance_matrices_long[
    covariance_matrices_long$scale == "Standardized residual" &
      covariance_matrices_long$matrix_id %in% matrix_ids,
    ,
    drop = FALSE
  ]
  plot_data$Panel <- factor(
    plot_data$matrix_id,
    levels = c(
      "A_dynamic_null__pairwise",
      "B_zero_eqtl_null__pairwise",
      "B_zero_eqtl_null__direct"
    ),
    labels = c(
      "A: pairwise difference",
      "B: pairwise difference",
      "B: direct zero mean"
    )
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = time_a, y = time_b, fill = covariance)
  ) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~Panel, nrow = 1) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "Covariance"
    ) +
    ggplot2::scale_x_continuous(breaks = c(0, 5, 10, 15)) +
    ggplot2::scale_y_continuous(breaks = c(0, 5, 10, 15)) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Day", y = "Day") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

plot_beta_scale_covariance_heatmaps <- function() {
  matrix_ids <- c(
    "A_dynamic_null__pairwise",
    "B_zero_eqtl_null__pairwise",
    "B_zero_eqtl_null__direct"
  )
  plot_data <- covariance_matrices_long[
    covariance_matrices_long$scale == "Median beta-hat scale" &
      covariance_matrices_long$matrix_id %in% matrix_ids,
    ,
    drop = FALSE
  ]
  plot_data$Panel <- factor(
    plot_data$matrix_id,
    levels = c(
      "A_dynamic_null__pairwise",
      "B_zero_eqtl_null__pairwise",
      "B_zero_eqtl_null__direct"
    ),
    labels = c(
      "A: pairwise difference",
      "B: pairwise difference",
      "B: direct zero mean"
    )
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = time_a, y = time_b, fill = covariance)
  ) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~Panel, nrow = 1) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "Covariance"
    ) +
    ggplot2::scale_x_continuous(breaks = c(0, 5, 10, 15)) +
    ggplot2::scale_y_continuous(breaks = c(0, 5, 10, 15)) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Day", y = "Day") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

prepare_interval_plot_data <- function() {
  plot_data <- bootstrap_lag_intervals
  plot_data$Method <- factor(
    paste(plot_data$set_id, plot_data$estimator, sep = "__"),
    levels = c(
      "A_dynamic_null__Pairwise difference",
      "B_zero_eqtl_null__Pairwise difference",
      "B_zero_eqtl_null__Direct zero mean"
    ),
    labels = c(
      "A: pairwise difference",
      "B: pairwise difference",
      "B: direct zero mean"
    )
  )
  plot_data[!is.na(plot_data$Method), , drop = FALSE]
}

plot_variograms <- function() {
  plot_data <- prepare_interval_plot_data()
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = lag,
      y = observed_semivariogram,
      color = Method,
      fill = Method
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = semivariogram_ci_lower,
        ymax = semivariogram_ci_upper
      ),
      alpha = 0.13,
      linewidth = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13, 15)) +
    ggplot2::labs(
      x = "Lag (days)",
      y = "Semivariogram",
      color = NULL
    ) +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_lag_covariances <- function() {
  plot_data <- prepare_interval_plot_data()
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = lag,
      y = observed_covariance,
      color = Method,
      fill = Method
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = covariance_ci_lower,
        ymax = covariance_ci_upper
      ),
      alpha = 0.13,
      linewidth = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13, 15)) +
    ggplot2::labs(
      x = "Lag (days)",
      y = "Mean standardized covariance",
      color = NULL
    ) +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_set_b_estimator_agreement <- function() {
  time_grid <- as.numeric(sub("^time_", "", colnames(pairwise_b)))
  grid <- expand.grid(
    time_a = time_grid,
    time_b = time_grid,
    stringsAsFactors = FALSE
  )
  grid$pairwise <- as.vector(pairwise_b)
  grid$direct <- as.vector(direct_b)
  grid <- grid[grid$time_a < grid$time_b, , drop = FALSE]
  grid$lag <- grid$time_b - grid$time_a

  ggplot2::ggplot(
    grid,
    ggplot2::aes(x = pairwise, y = direct, color = lag)
  ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      color = "grey45",
      linetype = "dashed"
    ) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::scale_color_viridis_c(name = "Lag") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = "Pairwise-difference covariance",
      y = "Direct zero-mean covariance"
    ) +
    ggplot2::theme_bw(base_size = 12)
}

# Set A follow-up: trajectory gallery and within-unit centering comparison.

set_a_metadata <- selected_units[
  selected_units$set_id == "A_dynamic_null",
  ,
  drop = FALSE
]
set_a_beta_hat <- selected_data$set_a$beta_hat
set_a_adjusted_se <- selected_data$set_a$adjusted_se
set_a_time_grid <- selected_data$set_a$time_grid
set_a_order <- match(rownames(set_a_beta_hat), set_a_metadata$pair_key)
if (anyNA(set_a_order) || anyDuplicated(set_a_order) ||
    !identical(dim(set_a_beta_hat), dim(set_a_adjusted_se)) ||
    nrow(set_a_beta_hat) != nrow(set_a_metadata) ||
    ncol(set_a_beta_hat) != length(set_a_time_grid)) {
  stop("The cached Set A metadata and fitted trajectories are misaligned.")
}
set_a_metadata <- set_a_metadata[set_a_order, , drop = FALSE]
if (!identical(set_a_metadata$pair_key, rownames(set_a_beta_hat))) {
  stop("The Set A pair-key alignment failed.")
}

set_a_centered_fit <- weighted_center_standardize(
  set_a_beta_hat,
  set_a_adjusted_se
)
set_a_gallery <- sample_lfdr_quartile_gallery(
  set_a_metadata,
  n_per_quartile = 25L,
  seed = 20260891L
)
set_a_gallery$constant_estimate <- unname(
  set_a_centered_fit$constant_estimate[
    match(set_a_gallery$pair_key, rownames(set_a_beta_hat))
  ]
)

set_a_pairwise_covariance <-
  pilot_result$matrices$A_dynamic_null__pairwise
set_a_centered_covariance <-
  pilot_result$matrices$A_dynamic_null__within_unit_centered_covariance
set_a_centered_correlation <- centered_correlations[[
  "A_dynamic_null__within_unit_centered_correlation"
]]
if (!identical(dim(set_a_pairwise_covariance),
               dim(set_a_centered_covariance)) ||
    !identical(dim(set_a_pairwise_covariance),
               dim(set_a_centered_correlation))) {
  stop("The Set A comparison matrices have incompatible dimensions.")
}

set_a_matrix_plot_data <- rbind(
  transform(
    covariance_matrix_to_long(
      set_a_pairwise_covariance,
      "A_dynamic_null",
      "Pairwise difference",
      "Standardized residual"
    ),
    Panel = "Pairwise-difference covariance"
  ),
  transform(
    covariance_matrix_to_long(
      set_a_centered_covariance,
      "A_dynamic_null",
      "Within-unit centered covariance",
      "Standardized residual"
    ),
    Panel = "Within-unit centered covariance"
  ),
  transform(
    covariance_matrix_to_long(
      set_a_centered_correlation,
      "A_dynamic_null",
      "Within-unit centered correlation",
      "Standardized residual"
    ),
    Panel = "Within-unit centered correlation"
  )
)
set_a_matrix_plot_data$Panel <- factor(
  set_a_matrix_plot_data$Panel,
  levels = c(
    "Pairwise-difference covariance",
    "Within-unit centered covariance",
    "Within-unit centered correlation"
  )
)

set_a_bootstrap_comparison <- bootstrap_lag_intervals[
  bootstrap_lag_intervals$set_id == "A_dynamic_null" &
    bootstrap_lag_intervals$estimator %in% c(
      "Pairwise difference",
      "Within-unit centered covariance"
    ),
  ,
  drop = FALSE
]
set_a_bootstrap_comparison$Estimator <- factor(
  set_a_bootstrap_comparison$estimator,
  levels = c(
    "Pairwise difference",
    "Within-unit centered covariance"
  ),
  labels = c(
    "Pairwise difference",
    "Within-unit centered"
  )
)

set_a_centered_correlation_lags <- lag_covariance_variogram(
  set_a_centered_correlation
)
set_a_key_lag_table <- set_a_bootstrap_comparison[
  set_a_bootstrap_comparison$lag %in% key_lags,
  c(
    "estimator", "lag", "observed_covariance",
    "covariance_ci_lower", "covariance_ci_upper",
    "observed_semivariogram", "semivariogram_ci_lower",
    "semivariogram_ci_upper"
  ),
  drop = FALSE
]
set_a_key_lag_table$estimator <- unname(
  estimator_labels[set_a_key_lag_table$estimator]
)
names(set_a_key_lag_table) <- c(
  "Estimator", "Lag", "Covariance", "Covariance 2.5%",
  "Covariance 97.5%", "Semivariogram", "Semivariogram 2.5%",
  "Semivariogram 97.5%"
)

set_a_centered_correlation_table <- data.frame(
  Lag = key_lags,
  `Centered correlation` = set_a_centered_correlation_lags[
    match(key_lags, set_a_centered_correlation_lags$lag),
    "mean_standardized_covariance"
  ],
  `Correlation-scale semivariogram` = set_a_centered_correlation_lags[
    match(key_lags, set_a_centered_correlation_lags$lag),
    "mean_semivariogram"
  ],
  check.names = FALSE
)

set_a_off_diagonal <- upper.tri(set_a_pairwise_covariance)
set_a_pairwise_entries <- set_a_pairwise_covariance[set_a_off_diagonal]
set_a_centered_entries <- set_a_centered_covariance[set_a_off_diagonal]
set_a_pairwise_lags <- lag_covariance_variogram(set_a_pairwise_covariance)
set_a_centered_lags <- lag_covariance_variogram(set_a_centered_covariance)
set_a_method_agreement <- data.frame(
  Quantity = c(
    "Off-diagonal Pearson correlation",
    "Mean centered-minus-pairwise difference",
    "Off-diagonal mean absolute difference",
    "Off-diagonal root mean squared difference",
    "Off-diagonal maximum absolute difference",
    "Maximum absolute lag-averaged covariance difference"
  ),
  Value = c(
    stats::cor(set_a_pairwise_entries, set_a_centered_entries),
    mean(set_a_centered_entries - set_a_pairwise_entries),
    mean(abs(set_a_centered_entries - set_a_pairwise_entries)),
    sqrt(mean((set_a_centered_entries - set_a_pairwise_entries)^2)),
    max(abs(set_a_centered_entries - set_a_pairwise_entries)),
    max(abs(
      set_a_pairwise_lags$mean_standardized_covariance -
        set_a_centered_lags$mean_standardized_covariance
    ))
  ),
  stringsAsFactors = FALSE
)

set_a_centered_diagnostic_table <- diagnostic_table[
  diagnostic_table$Set == set_labels[["A_dynamic_null"]] &
    diagnostic_table$Estimator %in% c(
      "Pairwise difference",
      "Within-unit centered covariance"
    ),
  ,
  drop = FALSE
]

set_a_gallery_table <- function(quartile) {
  quartile <- as.integer(quartile)
  if (length(quartile) != 1L || is.na(quartile) ||
      !quartile %in% seq_len(4L)) {
    stop("quartile must be one of 1, 2, 3, or 4.")
  }
  output <- set_a_gallery[
    set_a_gallery$lfdr_quartile == quartile,
    c(
      "gallery_page_position", "set_a_lfdr_rank", "gene_id",
      "variant_id", "bf_adjusted_lfdr", "constant_estimate"
    ),
    drop = FALSE
  ]
  names(output) <- c(
    "Panel", "Set A lfdr rank", "Gene", "Variant",
    "BF-adjusted lfdr", "Weighted constant estimate"
  )
  output
}

plot_set_a_gallery_page <- function(quartile) {
  quartile <- as.integer(quartile)
  page <- set_a_gallery[set_a_gallery$lfdr_quartile == quartile, , drop = FALSE]
  if (nrow(page) != 25L) {
    stop("Every gallery page must contain exactly 25 units.")
  }

  old_parameters <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_parameters), add = TRUE)
  graphics::par(
    mfrow = c(5, 5),
    mar = c(1.8, 2.2, 2.6, 0.6),
    oma = c(3.0, 3.2, 3.2, 0.8),
    mgp = c(1.2, 0.35, 0),
    tcl = -0.2
  )

  for (panel in seq_len(nrow(page))) {
    unit_index <- match(page$pair_key[panel], rownames(set_a_beta_hat))
    beta_value <- set_a_beta_hat[unit_index, ]
    se_value <- set_a_adjusted_se[unit_index, ]
    lower <- beta_value - se_value
    upper <- beta_value + se_value
    y_limits <- range(c(
      lower,
      upper,
      page$constant_estimate[panel]
    ))
    y_padding <- max(diff(y_limits) * 0.06, 1e-6)
    y_limits <- y_limits + c(-y_padding, y_padding)

    graphics::plot(
      set_a_time_grid,
      beta_value,
      type = "n",
      xlim = range(set_a_time_grid),
      ylim = y_limits,
      xlab = "",
      ylab = "",
      axes = FALSE,
      main = sprintf(
        "Rank %d: %s\n%s; lfdr %.3f",
        page$set_a_lfdr_rank[panel],
        page$gene_id[panel],
        page$variant_id[panel],
        page$bf_adjusted_lfdr[panel]
      ),
      cex.main = 0.50
    )
    graphics::axis(
      1,
      at = c(0, 5, 10, 15),
      labels = c(0, 5, 10, 15),
      cex.axis = 0.55
    )
    graphics::axis(2, las = 1, cex.axis = 0.48)
    graphics::box(col = "grey65")
    graphics::segments(
      set_a_time_grid,
      lower,
      set_a_time_grid,
      upper,
      col = "grey65",
      lwd = 0.7
    )
    graphics::abline(
      h = page$constant_estimate[panel],
      col = "#D55E00",
      lty = 2,
      lwd = 0.9
    )
    graphics::lines(
      set_a_time_grid,
      beta_value,
      col = "#0072B2",
      lwd = 0.8
    )
    graphics::points(
      set_a_time_grid,
      beta_value,
      pch = 16,
      col = "#0072B2",
      cex = 0.55
    )
  }

  rank_range <- range(page$set_a_lfdr_rank)
  lfdr_range <- range(page$bf_adjusted_lfdr)
  graphics::mtext("Day", side = 1, outer = TRUE, line = 1.3, cex = 0.85)
  graphics::mtext(
    expression(hat(beta)),
    side = 2,
    outer = TRUE,
    line = 1.7,
    cex = 0.85
  )
  graphics::mtext(
    sprintf(
      "Set A lfdr quartile %d: 25 fixed-random units (ranks %d-%d; lfdr %.3f-%.3f)",
      quartile,
      rank_range[1],
      rank_range[2],
      lfdr_range[1],
      lfdr_range[2]
    ),
    side = 3,
    outer = TRUE,
    line = 1.2,
    cex = 0.9,
    font = 2
  )
  invisible(page)
}

plot_set_a_matrix_comparison <- function() {
  ggplot2::ggplot(
    set_a_matrix_plot_data,
    ggplot2::aes(x = time_a, y = time_b, fill = covariance)
  ) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~Panel, nrow = 1) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-0.3, 1),
      oob = scales::squish,
      name = "Entry"
    ) +
    ggplot2::scale_x_continuous(breaks = c(0, 5, 10, 15)) +
    ggplot2::scale_y_continuous(breaks = c(0, 5, 10, 15)) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Day", y = "Day") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 9)
    )
}

plot_set_a_covariance_comparison <- function() {
  ggplot2::ggplot(
    set_a_bootstrap_comparison,
    ggplot2::aes(
      x = lag,
      y = observed_covariance,
      color = Estimator,
      fill = Estimator
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey55", linewidth = 0.4) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = covariance_ci_lower, ymax = covariance_ci_upper),
      alpha = 0.14,
      linewidth = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00")) +
    ggplot2::scale_fill_manual(values = c("#0072B2", "#D55E00")) +
    ggplot2::scale_x_continuous(breaks = seq(1, 15, by = 2)) +
    ggplot2::labs(
      x = "Lag (days)",
      y = "Mean standardized covariance",
      color = NULL
    ) +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_set_a_variogram_comparison <- function() {
  ggplot2::ggplot(
    set_a_bootstrap_comparison,
    ggplot2::aes(
      x = lag,
      y = observed_semivariogram,
      color = Estimator,
      fill = Estimator
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        ymin = semivariogram_ci_lower,
        ymax = semivariogram_ci_upper
      ),
      alpha = 0.14,
      linewidth = 0,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.8) +
    ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00")) +
    ggplot2::scale_fill_manual(values = c("#0072B2", "#D55E00")) +
    ggplot2::scale_x_continuous(breaks = seq(1, 15, by = 2)) +
    ggplot2::labs(
      x = "Lag (days)",
      y = "Semivariogram",
      color = NULL
    ) +
    ggplot2::guides(fill = "none") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}

plot_set_a_covariance_entry_agreement <- function() {
  time_grid <- as.numeric(sub(
    "^time_",
    "",
    colnames(set_a_pairwise_covariance)
  ))
  grid <- expand.grid(
    time_a = time_grid,
    time_b = time_grid,
    stringsAsFactors = FALSE
  )
  grid$pairwise <- as.vector(set_a_pairwise_covariance)
  grid$centered <- as.vector(set_a_centered_covariance)
  grid <- grid[grid$time_a < grid$time_b, , drop = FALSE]
  grid$lag <- grid$time_b - grid$time_a

  ggplot2::ggplot(
    grid,
    ggplot2::aes(x = pairwise, y = centered, color = lag)
  ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      color = "grey45",
      linetype = "dashed"
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
    ggplot2::geom_point(size = 2, alpha = 0.8) +
    ggplot2::scale_color_viridis_c(name = "Lag") +
    ggplot2::coord_equal() +
    ggplot2::labs(
      x = "Pairwise-difference covariance",
      y = "Within-unit centered covariance"
    ) +
    ggplot2::theme_bw(base_size = 12)
}
