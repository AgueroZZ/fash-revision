# Reporting helpers for the mashr-style mean-Z null-correlation comparison.

source(file.path(
  "code",
  "revision_simulations",
  "r4_correlated_errors",
  "real_data_correlation_helpers.R"
))
source(file.path(
  "code",
  "revision_simulations",
  "internal",
  "covariance_estimation",
  "zero_intercept_correlation_helpers.R"
))

mashr_mean_z_cache_dir <- file.path(
  "output",
  "revision_simulations",
  "internal",
  "mashr_mean_z_null_correlation"
)
mashr_mean_z_configuration <- readRDS(file.path(
  mashr_mean_z_cache_dir,
  "configuration.rds"
))
mashr_mean_z_analysis <- readRDS(file.path(
  mashr_mean_z_cache_dir,
  "mashr_mean_z_null_correlation.rds"
))
mashr_screen_fit <- readRDS(file.path(
  mashr_mean_z_cache_dir,
  "mashr_screen_fit.rds"
))

maximum_filter_id <- "maximum_z"
combined_filter_id <- "maximum_z_and_mean_z"
mashr_filter_id <- "maximum_z_and_mashr_lfdr"
expected_filter_ids <- c(
  maximum_filter_id,
  combined_filter_id,
  mashr_filter_id
)

validate_mashr_mean_z_cache <- function() {
  configuration <- mashr_mean_z_configuration
  analysis <- mashr_mean_z_analysis
  metadata <- analysis$candidate_metadata
  z <- analysis$candidate_z

  if (!identical(configuration$analysis_id, "mashr_mean_z_null_correlation") ||
      !identical(as.integer(configuration$thinning_seed), 20260811L) ||
      !isTRUE(all.equal(configuration$max_threshold, 2)) ||
      !isTRUE(all.equal(configuration$mean_z_threshold, 2)) ||
      !isTRUE(all.equal(configuration$mashr_pair_lfdr_threshold, 0.05)) ||
      !identical(as.integer(configuration$mashr_seed), 123L) ||
      !identical(
        configuration$mashr_null_correlation_function,
        "mashr::estimate_null_correlation_simple"
      ) ||
      !identical(
        configuration$mashr_pair_lfdr_definition,
        "posterior weight on the exact all-zero component"
      ) ||
      !identical(
        configuration$mashr_cov_methods,
        c("identity", "singletons", "equal_effects", "simple_het")
      ) ||
      !identical(configuration$mashr_prior, "nullbiased") ||
      !isTRUE(all.equal(configuration$mashr_nullweight, 10)) ||
      !identical(configuration$mashr_optmethod, "mixSQP") ||
      !identical(as.integer(configuration$bootstrap_reps), 1000L) ||
      !identical(
        as.integer(configuration$bootstrap_seeds),
        c(20260831L, 20260832L, 20260831L)
      ) ||
      !identical(as.integer(configuration$gallery_size), 25L) ||
      !identical(as.integer(configuration$gallery_seeds), 20260841:20260843) ||
      !identical(as.integer(configuration$n_time), 16L) ||
      !isTRUE(all.equal(configuration$time_grid, 0:15)) ||
      !identical(configuration$mean_z_definition, "sqrt(n_time) * rowMeans(z)") ||
      !identical(configuration$estimator, "stats::cor(z_selected)")) {
    stop("The mashr mean-Z configuration is incomplete or unexpected.")
  }

  fit_info <- file.info(configuration$fit_path)
  if (!file.exists(configuration$fit_path) ||
      unname(fit_info$size) != configuration$fit_size ||
      as.character(fit_info$mtime) != configuration$fit_mtime) {
    stop("The source FASH fit no longer matches the cached analysis.")
  }

  required_metadata_columns <- c(
    "seed", "fash_index", "pair_key", "gene_id", "variant_id", "lfdr",
    "max_absolute_z", "mean_z_score", "passes_maximum_z",
    "passes_mean_z", "passes_maximum_z_and_mean_z", "mashr_pair_lfdr",
    "mashr_min_condition_lfdr", "mashr_max_condition_lfdr",
    "passes_mashr_pair_lfdr", "passes_maximum_z_and_mashr_lfdr",
    "candidate_row"
  )
  if (!is.data.frame(metadata) ||
      !all(required_metadata_columns %in% names(metadata)) ||
      nrow(metadata) != configuration$n_genes ||
      anyDuplicated(metadata$gene_id) || anyDuplicated(metadata$pair_key) ||
      any(!is.finite(metadata$lfdr)) || any(metadata$lfdr < 0 | metadata$lfdr > 1) ||
      !is.matrix(z) || !identical(dim(z), c(configuration$n_genes, 16L)) ||
      !identical(rownames(z), metadata$pair_key) || any(!is.finite(z))) {
    stop("The thinned candidate data are incomplete or invalid.")
  }

  expected_maximum <- apply(abs(z), 1L, max)
  expected_mean_z <- sqrt(ncol(z)) * rowMeans(z)
  expected_keep_maximum <- expected_maximum < configuration$max_threshold
  expected_keep_mean <- abs(expected_mean_z) < configuration$mean_z_threshold
  expected_keep_combined <- expected_keep_maximum & expected_keep_mean
  expected_keep_mashr <- metadata$mashr_pair_lfdr >
    configuration$mashr_pair_lfdr_threshold
  expected_keep_screen3 <- expected_keep_maximum & expected_keep_mashr
  if (max(abs(expected_maximum - metadata$max_absolute_z)) > 1e-12 ||
      max(abs(expected_mean_z - metadata$mean_z_score)) > 1e-12 ||
      max(abs(metadata$mashr_pair_lfdr - mashr_screen_fit$pair_lfdr)) > 1e-12 ||
      max(abs(
        metadata$mashr_min_condition_lfdr -
          apply(mashr_screen_fit$condition_lfdr, 1L, min)
      )) > 1e-12 ||
      max(abs(
        metadata$mashr_max_condition_lfdr -
          apply(mashr_screen_fit$condition_lfdr, 1L, max)
      )) > 1e-12 ||
      !isTRUE(all.equal(
        metadata$passes_maximum_z,
        expected_keep_maximum,
        check.attributes = FALSE
      )) ||
      !isTRUE(all.equal(
        metadata$passes_mean_z,
        expected_keep_mean,
        check.attributes = FALSE
      )) ||
      !isTRUE(all.equal(
        metadata$passes_maximum_z_and_mean_z,
        expected_keep_combined,
        check.attributes = FALSE
      )) ||
      !isTRUE(all.equal(
        metadata$passes_mashr_pair_lfdr,
        expected_keep_mashr,
        check.attributes = FALSE
      )) ||
      !isTRUE(all.equal(
        metadata$passes_maximum_z_and_mashr_lfdr,
        expected_keep_screen3,
        check.attributes = FALSE
      )) || any(expected_keep_combined & !expected_keep_maximum) ||
      any(expected_keep_screen3 & !expected_keep_maximum)) {
    stop("The cached screening statistics do not match their definitions.")
  }

  expected_null_correlation <- stats::cor(
    z[expected_keep_maximum, , drop = FALSE]
  )
  if (!identical(dim(mashr_screen_fit$condition_lfdr), dim(z)) ||
      length(mashr_screen_fit$pair_lfdr) != nrow(z) ||
      max(abs(
        mashr_screen_fit$null_correlation - expected_null_correlation
      )) > 1e-12 ||
      max(abs(
        analysis$mashr_null_correlation - expected_null_correlation
      )) > 1e-12 ||
      !isTRUE(all.equal(
        configuration$mashr_fitted_pi0,
        mashr_screen_fit$fitted_pi0
      )) ||
      configuration$mashr_n_nullish_for_correlation != sum(
        expected_keep_maximum
      )) {
    stop("The cached mashr fit or null correlation is inconsistent.")
  }

  if (!identical(
    as.character(analysis$filter_definitions$filter_id),
    expected_filter_ids
  ) || nrow(analysis$selection_counts) != 3L ||
      !identical(
        as.character(analysis$selection_counts$filter_id),
        expected_filter_ids
      ) ||
      !identical(
        as.integer(analysis$selection_counts$n_selected),
        c(
          sum(expected_keep_maximum),
          sum(expected_keep_combined),
          sum(expected_keep_screen3)
        )
      )) {
    stop("The cached filter definitions or selection counts are invalid.")
  }

  expected_keeps <- list(
    expected_keep_maximum,
    expected_keep_combined,
    expected_keep_screen3
  )
  names(expected_keeps) <- expected_filter_ids
  for (filter_id in expected_filter_ids) {
    estimate <- analysis$estimates[[filter_id]]
    selected_z <- z[expected_keeps[[filter_id]], , drop = FALSE]
    observed_correlation <- estimate$estimate$sample_correlation
    if (!identical(estimate$selected_pair_keys, rownames(selected_z)) ||
        !isTRUE(all.equal(estimate$z_selected, selected_z, tolerance = 0)) ||
        !isTRUE(all.equal(
          observed_correlation,
          stats::cor(selected_z),
          tolerance = 1e-12
        ))) {
      stop("A cached selected set or empirical correlation is inconsistent.")
    }
    validate_correlation_matrix(observed_correlation)
  }

  bootstrap <- analysis$bootstrap_intervals
  if (!is.data.frame(bootstrap) || nrow(bootstrap) != 90L ||
      !setequal(bootstrap$filter_id, expected_filter_ids) ||
      !setequal(bootstrap$summary_type, c("adjacent_pair", "lag_average")) ||
      any(!is.finite(unlist(bootstrap[c(
        "point_estimate", "mean", "median", "lower", "upper",
        "semivariogram_point", "semivariogram_lower",
        "semivariogram_upper"
      )]))) || any(bootstrap$lower > bootstrap$upper) ||
      any(bootstrap$semivariogram_lower > bootstrap$semivariogram_upper)) {
    stop("The gene-bootstrap intervals are incomplete or invalid.")
  }

  gallery <- analysis$gallery_membership
  gallery_counts <- table(factor(gallery$filter_id, levels = expected_filter_ids))
  gallery_indices <- match(gallery$pair_key, metadata$pair_key)
  expected_gallery_keep <- vapply(seq_len(nrow(gallery)), function(index) {
    expected_keeps[[gallery$filter_id[index]]][gallery_indices[index]]
  }, logical(1))
  observed_overlaps <- analysis$gallery_overlaps
  recomputed_overlaps <- observed_overlaps
  for (index in seq_len(nrow(recomputed_overlaps))) {
    first_id <- recomputed_overlaps$first_filter_id[index]
    second_id <- recomputed_overlaps$second_filter_id[index]
    recomputed_overlaps$n_overlap[index] <- length(intersect(
      gallery$pair_key[gallery$filter_id == first_id],
      gallery$pair_key[gallery$filter_id == second_id]
    ))
  }
  if (!is.data.frame(gallery) || nrow(gallery) != 75L ||
      anyNA(gallery_indices) || any(gallery_counts != 25L) ||
      !all(expected_gallery_keep) ||
      anyDuplicated(interaction(gallery$filter_id, gallery$pair_key)) ||
      !identical(observed_overlaps, recomputed_overlaps) ||
      !identical(configuration$gallery_overlaps, observed_overlaps)) {
    stop("The fixed trajectory galleries are incomplete or invalid.")
  }

  if (!is.data.frame(analysis$matrix_comparison) ||
      nrow(analysis$matrix_comparison) != 3L ||
      !is.data.frame(analysis$correlation_difference_long) ||
      nrow(analysis$correlation_difference_long) != 768L ||
      any(!is.finite(unlist(analysis$matrix_comparison[, 4:8]))) ||
      any(!is.finite(analysis$correlation_difference_long$correlation_difference))) {
    stop("The cached correlation-matrix comparison is invalid.")
  }

  bad_example <- analysis$bad_example_audit
  bad_index <- match(bad_example$pair_key, metadata$pair_key)
  if (!is.data.frame(bad_example) || nrow(bad_example) != 1L ||
      is.na(bad_index) || !isTRUE(bad_example$passes_screen1) ||
      isTRUE(bad_example$passes_screen2) ||
      !isTRUE(bad_example$passes_screen3) ||
      abs(bad_example$mashr_pair_lfdr - metadata$mashr_pair_lfdr[bad_index]) >
        1e-12 ||
      analysis$mashr_screen_summary$n_screen1_pair_lfdr_at_or_below_threshold !=
        0L ||
      !identical(expected_keep_screen3, expected_keep_maximum)) {
    stop("The prespecified bad-example or Screen-3 audit is inconsistent.")
  }

  invisible(TRUE)
}

validate_mashr_mean_z_cache()

format_decimal <- function(x, digits = 3L) {
  formatC(x, format = "f", digits = digits)
}

format_integer <- function(x) {
  formatC(as.integer(x), format = "d", big.mark = ",")
}

format_percent <- function(x, digits = 1L) {
  paste0(format_decimal(100 * x, digits), "%")
}

render_scrollable_table <- function(x,
                                    align = NULL,
                                    minimum_width = "900px",
                                    digits = NULL) {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("The knitr package is required to render tables.")
  }
  arguments <- list(
    x = x,
    format = "html",
    escape = TRUE,
    align = align
  )
  if (!is.null(digits)) {
    arguments$digits <- digits
  }
  rendered <- do.call(knitr::kable, arguments)
  cat(
    '<div class="mean-z-table-scroll">\n',
    '<div style="min-width: ', minimum_width, ';">\n',
    rendered,
    '\n</div>\n',
    '</div>\n',
    sep = ""
  )
  invisible(x)
}

selection_counts <- mashr_mean_z_analysis$selection_counts
n_candidates <- selection_counts$n_candidates[1]
n_maximum <- selection_counts$n_selected[
  selection_counts$filter_id == maximum_filter_id
]
n_combined <- selection_counts$n_selected[
  selection_counts$filter_id == combined_filter_id
]
n_mashr <- selection_counts$n_selected[
  selection_counts$filter_id == mashr_filter_id
]
n_removed_by_mean_z <- n_maximum - n_combined
fraction_removed_from_maximum <- n_removed_by_mean_z / n_maximum
n_removed_by_mashr <- n_maximum - n_mashr

selection_table <- data.frame(
  Screen = selection_counts$filter_label,
  Candidates = selection_counts$n_candidates,
  Retained = selection_counts$n_selected,
  `Fraction of candidates retained` = selection_counts$selected_fraction,
  `Maximum absolute time-column mean z` =
    selection_counts$maximum_absolute_column_mean,
  `Minimum mashr pair lfdr` = selection_counts$minimum_mashr_pair_lfdr,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

bootstrap_lags <- mashr_mean_z_analysis$bootstrap_intervals[
  mashr_mean_z_analysis$bootstrap_intervals$summary_type == "lag_average",
  ,
  drop = FALSE
]
key_lag_table <- bootstrap_lags[
  bootstrap_lags$index %in% c(1L, 5L, 10L, 15L),
  ,
  drop = FALSE
]
key_lag_table <- data.frame(
  Screen = key_lag_table$filter_label,
  Lag = key_lag_table$index,
  Correlation = key_lag_table$point_estimate,
  `Correlation 95% interval` = paste0(
    "[",
    format_decimal(key_lag_table$lower),
    ", ",
    format_decimal(key_lag_table$upper),
    "]"
  ),
  Semivariogram = key_lag_table$semivariogram_point,
  `Semivariogram 95% interval` = paste0(
    "[",
    format_decimal(key_lag_table$semivariogram_lower),
    ", ",
    format_decimal(key_lag_table$semivariogram_upper),
    "]"
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

filter_label_lookup <- stats::setNames(
  mashr_mean_z_analysis$filter_definitions$filter_label,
  mashr_mean_z_analysis$filter_definitions$filter_id
)
matrix_comparison_table <- data.frame(
  Reference = unname(filter_label_lookup[
    mashr_mean_z_analysis$matrix_comparison$reference_filter_id
  ]),
  Comparison = unname(filter_label_lookup[
    mashr_mean_z_analysis$matrix_comparison$comparison_filter_id
  ]),
  `Off-diagonal matrix correlation` =
    mashr_mean_z_analysis$matrix_comparison$off_diagonal_matrix_correlation,
  `Mean comparison-minus-reference difference` =
    mashr_mean_z_analysis$matrix_comparison$mean_off_diagonal_difference,
  `Mean absolute difference` =
    mashr_mean_z_analysis$matrix_comparison$mean_absolute_off_diagonal_difference,
  `Maximum absolute difference` =
    mashr_mean_z_analysis$matrix_comparison$maximum_absolute_off_diagonal_difference,
  `Frobenius difference` =
    mashr_mean_z_analysis$matrix_comparison$full_frobenius_difference,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

mashr_screen_summary_table <- data.frame(
  Metric = c(
    "All thinned pairs",
    "Screen 1 pairs",
    "All-pair mashr lfdr <= 0.05",
    "Screen-1 mashr lfdr <= 0.05",
    "Minimum Screen-1 mashr pair lfdr",
    "Mashr fitted pi0",
    "Null-C versus Screen-1 maximum difference"
  ),
  Value = c(
    mashr_mean_z_analysis$mashr_screen_summary$n_all_thinned_pairs,
    mashr_mean_z_analysis$mashr_screen_summary$n_screen1_pairs,
    mashr_mean_z_analysis$mashr_screen_summary$n_pair_lfdr_at_or_below_threshold,
    mashr_mean_z_analysis$mashr_screen_summary$n_screen1_pair_lfdr_at_or_below_threshold,
    mashr_mean_z_analysis$mashr_screen_summary$minimum_screen1_pair_lfdr,
    mashr_mean_z_analysis$mashr_screen_summary$fitted_pi0,
    mashr_mean_z_analysis$mashr_screen_summary$null_correlation_maximum_difference
  ),
  stringsAsFactors = FALSE
)

bad_example_table <- data.frame(
  Gene = mashr_mean_z_analysis$bad_example_audit$gene_id,
  Variant = mashr_mean_z_analysis$bad_example_audit$variant_id,
  `Beta-hat range` = paste0(
    "[",
    format_decimal(mashr_mean_z_analysis$bad_example_audit$minimum_beta_hat),
    ", ",
    format_decimal(mashr_mean_z_analysis$bad_example_audit$maximum_beta_hat),
    "]"
  ),
  `IVW constant beta` =
    mashr_mean_z_analysis$bad_example_audit$inverse_variance_weighted_constant_beta,
  `max|z|` = mashr_mean_z_analysis$bad_example_audit$maximum_absolute_z,
  `Z_mean` = mashr_mean_z_analysis$bad_example_audit$mean_z_score,
  `Mashr pair lfdr` = mashr_mean_z_analysis$bad_example_audit$mashr_pair_lfdr,
  `Minimum condition lfdr` =
    mashr_mean_z_analysis$bad_example_audit$mashr_min_condition_lfdr,
  `Screen 1` = mashr_mean_z_analysis$bad_example_audit$passes_screen1,
  `Screen 2` = mashr_mean_z_analysis$bad_example_audit$passes_screen2,
  `Screen 3` = mashr_mean_z_analysis$bad_example_audit$passes_screen3,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

maximum_correlation <- mashr_mean_z_analysis$estimates[[
  maximum_filter_id
]]$estimate$sample_correlation
combined_correlation <- mashr_mean_z_analysis$estimates[[
  combined_filter_id
]]$estimate$sample_correlation
mashr_correlation <- mashr_mean_z_analysis$estimates[[
  mashr_filter_id
]]$estimate$sample_correlation
maximum_lag <- lag_average_correlation(maximum_correlation)
combined_lag <- lag_average_correlation(combined_correlation)
mashr_lag <- lag_average_correlation(mashr_correlation)

plot_correlation_panel <- function(correlation, title, z_limit) {
  display_matrix <- correlation
  diag(display_matrix) <- NA_real_
  palette <- grDevices::colorRampPalette(
    c("#2166AC", "white", "#B2182B")
  )(201)
  graphics::image(
    x = 0:15,
    y = 0:15,
    z = display_matrix,
    zlim = c(-z_limit, z_limit),
    col = palette,
    axes = FALSE,
    xlab = "Time",
    ylab = "Time",
    main = title,
    cex.main = 0.9,
    asp = 1
  )
  graphics::axis(
    1,
    at = 0:15,
    labels = 0:15,
    las = 2,
    cex.axis = 0.62
  )
  graphics::axis(2, at = 0:15, labels = 0:15, las = 1, cex.axis = 0.68)
  graphics::box()
}

plot_mashr_mean_z_heatmaps <- function() {
  off_diagonal <- upper.tri(maximum_correlation)
  z_limit <- ceiling(20 * max(abs(c(
    maximum_correlation[off_diagonal],
    combined_correlation[off_diagonal],
    mashr_correlation[off_diagonal]
  )))) / 20
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(1, 3),
    mar = c(4.8, 4.0, 5.2, 0.7),
    oma = c(2.0, 0.3, 3.0, 0.3)
  )
  plot_correlation_panel(
    maximum_correlation,
    paste0("max|z(t)| < 2\n(n = ", format_integer(n_maximum), ")"),
    z_limit
  )
  plot_correlation_panel(
    combined_correlation,
    paste0(
      "max|z(t)| < 2 and |Z_mean| < 2\n(n = ",
      format_integer(n_combined),
      ")"
    ),
    z_limit
  )
  plot_correlation_panel(
    mashr_correlation,
    paste0(
      "max|z(t)| < 2 and\nmashr pair lfdr > 0.05\n(n = ",
      format_integer(n_mashr),
      ")"
    ),
    z_limit
  )
  graphics::mtext(
    "Empirical time correlation after one-variant-per-gene thinning",
    side = 3,
    outer = TRUE,
    line = 1.0,
    cex = 1.35,
    font = 2
  )
  graphics::mtext(
    paste0(
      "Shared off-diagonal scale: blue = -",
      format_decimal(z_limit, 2),
      ", white = 0, red = +",
      format_decimal(z_limit, 2),
      "; unit diagonal omitted"
    ),
    side = 1,
    outer = TRUE,
    line = 0.3,
    cex = 0.82
  )
}

plot_lag_metric <- function(metric = c("correlation", "semivariogram")) {
  metric <- match.arg(metric)
  bootstrap <- bootstrap_lags
  colors <- c(
    maximum_z = "#0072B2",
    maximum_z_and_mean_z = "#D55E00",
    maximum_z_and_mashr_lfdr = "#009E73"
  )
  point_characters <- c(
    maximum_z = 19,
    maximum_z_and_mean_z = 1,
    maximum_z_and_mashr_lfdr = 17
  )
  line_types <- c(
    maximum_z = 1,
    maximum_z_and_mean_z = 1,
    maximum_z_and_mashr_lfdr = 2
  )
  if (metric == "correlation") {
    point_column <- "point_estimate"
    lower_column <- "lower"
    upper_column <- "upper"
    y_label <- "Lag-averaged correlation"
    title <- "Lag correlation"
    reference <- 0
  } else {
    point_column <- "semivariogram_point"
    lower_column <- "semivariogram_lower"
    upper_column <- "semivariogram_upper"
    y_label <- "Standardized semivariogram"
    title <- "Variogram"
    reference <- 1
  }
  y_limits <- range(c(
    bootstrap[[lower_column]],
    bootstrap[[upper_column]],
    reference
  ))
  padding <- 0.05 * diff(y_limits)
  y_limits <- y_limits + c(-padding, padding)
  graphics::plot(
    1:15,
    rep(NA_real_, 15L),
    type = "n",
    xlab = "Time lag",
    ylab = y_label,
    ylim = y_limits,
    main = title,
    xaxt = "n"
  )
  graphics::axis(1, at = seq(1, 15, by = 2))
  graphics::abline(h = reference, lty = 3, col = "gray50")
  for (filter_id in rev(expected_filter_ids)) {
    subset <- bootstrap[bootstrap$filter_id == filter_id, , drop = FALSE]
    graphics::polygon(
      c(subset$index, rev(subset$index)),
      c(
        subset[[lower_column]],
        rev(subset[[upper_column]])
      ),
      col = grDevices::adjustcolor(colors[filter_id], alpha.f = 0.14),
      border = NA
    )
  }
  for (filter_id in expected_filter_ids) {
    subset <- bootstrap[bootstrap$filter_id == filter_id, , drop = FALSE]
    graphics::lines(
      subset$index,
      subset[[point_column]],
      type = "b",
      pch = point_characters[filter_id],
      lty = line_types[filter_id],
      lwd = 2,
      col = colors[filter_id]
    )
  }
  graphics::legend(
    "topright",
    legend = mashr_mean_z_analysis$filter_definitions$filter_label,
    col = colors[expected_filter_ids],
    pch = point_characters[expected_filter_ids],
    lty = line_types[expected_filter_ids],
    lwd = 2,
    bty = "n",
    cex = 0.72
  )
}

plot_lag_correlation_comparison <- function() {
  plot_lag_metric("correlation")
}

plot_variogram_comparison <- function() {
  plot_lag_metric("semivariogram")
}

get_gallery <- function(filter_id) {
  if (!filter_id %in% expected_filter_ids) {
    stop("Unknown gallery filter identifier.")
  }
  metadata <- mashr_mean_z_analysis$gallery_membership[
    mashr_mean_z_analysis$gallery_membership$filter_id == filter_id,
    ,
    drop = FALSE
  ]
  metadata <- metadata[order(metadata$gallery_position), , drop = FALSE]
  z_index <- match(metadata$pair_key, rownames(mashr_mean_z_analysis$candidate_z))
  z <- mashr_mean_z_analysis$candidate_z[z_index, , drop = FALSE]
  if (nrow(metadata) != 25L || anyNA(z_index) ||
      !identical(rownames(z), metadata$pair_key)) {
    stop("The requested gallery cannot be reconstructed.")
  }
  list(metadata = metadata, z = z)
}

gallery_table <- function(filter_id) {
  gallery <- get_gallery(filter_id)
  data.frame(
    Position = gallery$metadata$gallery_position,
    Gene = gallery$metadata$gene_id,
    Variant = gallery$metadata$variant_id,
    lfdr = gallery$metadata$lfdr,
    `max|z|` = gallery$metadata$max_absolute_z,
    `Z_mean` = gallery$metadata$mean_z_score,
    `Mashr pair lfdr` = gallery$metadata$mashr_pair_lfdr,
    `Passes Screen 2` = gallery$metadata$passes_maximum_z_and_mean_z,
    `Passes Screen 3` =
      gallery$metadata$passes_maximum_z_and_mashr_lfdr,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

plot_gallery <- function(filter_id) {
  gallery <- get_gallery(filter_id)
  metadata <- gallery$metadata
  z <- gallery$z
  color <- c(
    maximum_z = "#0072B2",
    maximum_z_and_mean_z = "#D55E00",
    maximum_z_and_mashr_lfdr = "#009E73"
  )[filter_id]
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(5, 5),
    mar = c(1.35, 1.45, 2.55, 0.45),
    oma = c(2.8, 3.0, 3.5, 0.6),
    mgp = c(1.4, 0.35, 0),
    tcl = -0.2
  )
  time <- 0:15
  for (panel in seq_len(nrow(metadata))) {
    graphics::plot(
      time,
      z[panel, ],
      type = "o",
      pch = 16,
      cex = 0.48,
      lwd = 1.0,
      col = color,
      xlim = c(0, 15),
      ylim = c(-2.05, 2.05),
      axes = FALSE,
      xlab = "",
      ylab = ""
    )
    graphics::abline(h = c(-2, 2), lty = 3, col = "gray70")
    graphics::abline(h = 0, lty = 1, col = "gray82")
    graphics::abline(h = mean(z[panel, ]), lty = 2, col = "#009E73")
    row_number <- ceiling(panel / 5)
    column_number <- ((panel - 1L) %% 5L) + 1L
    if (row_number == 5L) {
      graphics::axis(1, at = c(0, 5, 10, 15), cex.axis = 0.56)
    } else {
      graphics::axis(1, at = c(0, 5, 10, 15), labels = FALSE)
    }
    if (column_number == 1L) {
      graphics::axis(2, at = c(-2, 0, 2), las = 1, cex.axis = 0.56)
    } else {
      graphics::axis(2, at = c(-2, 0, 2), labels = FALSE)
    }
    graphics::box(col = "gray55")
    graphics::title(
      main = paste0(
        metadata$gene_id[panel],
        "\nmax=",
        format_decimal(metadata$max_absolute_z[panel], 2),
        "; Zm=",
        format_decimal(metadata$mean_z_score[panel], 2),
        "; mLFDR=",
        format_decimal(metadata$mashr_pair_lfdr[panel], 2)
      ),
      cex.main = 0.57,
      line = 0.45
    )
  }
  graphics::mtext("Time", side = 1, outer = TRUE, line = 1.2, cex = 0.9)
  graphics::mtext("z score", side = 2, outer = TRUE, line = 1.4, cex = 0.9)
  graphics::mtext(
    unique(metadata$filter_label),
    side = 3,
    outer = TRUE,
    line = 1.5,
    cex = 1.2,
    font = 2
  )
  graphics::mtext(
    "Green dashed line = within-trajectory mean z",
    side = 3,
    outer = TRUE,
    line = 0.1,
    cex = 0.78
  )
}
