# Reporting helpers for the internal zero-intercept correlation exploration.

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

internal_correlation_cache <- file.path(
  "output",
  "revision_simulations",
  "internal",
  "zero_intercept_correlation"
)
internal_correlation_summary_dir <- file.path(
  internal_correlation_cache,
  "summary"
)

configuration <- readRDS(file.path(
  internal_correlation_cache,
  "configuration.rds"
))
internal_analysis <- readRDS(file.path(
  internal_correlation_cache,
  "zero_intercept_correlation_analysis.rds"
))

read_summary_csv <- function(filename) {
  path <- file.path(internal_correlation_summary_dir, filename)
  if (!file.exists(path)) {
    stop("Missing internal correlation summary: ", path)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

selection_counts <- read_summary_csv("selection_counts.csv")
estimator_diagnostics <- read_summary_csv("estimator_diagnostics.csv")
correlation_matrices_long <- read_summary_csv("correlation_matrices_long.csv")
second_moments_long <- read_summary_csv("second_moments_long.csv")
adjacent_correlations <- read_summary_csv("adjacent_correlations.csv")
lag_summaries <- read_summary_csv("lag_summaries.csv")
bootstrap_intervals <- read_summary_csv("bootstrap_intervals.csv")
random_variant_stability <- read_summary_csv("random_variant_stability.csv")
random_lag_stability <- read_summary_csv("random_lag_stability.csv")
independence_benchmark <- read_summary_csv("independence_benchmark.csv")
truncation_calibration <- read_summary_csv("truncation_calibration.csv")

top_design <- "Top-1000 highest lfdr per gene"
random_design <- "Mean across 50 pre-z random selections"
primary_random_design <- paste0(
  "Random variant per gene before z, seed ",
  configuration$random_seeds[1]
)
expected_thresholds <- c(1.5, 2, 2.5)
expected_random_seeds <- 20260811:20260860

validate_internal_correlation_cache <- function() {
  if (!identical(configuration$analysis_id, "zero_intercept_correlation") ||
      !identical(as.integer(configuration$top_n), 1000L) ||
      !isTRUE(all.equal(configuration$thresholds, expected_thresholds)) ||
      !identical(as.integer(configuration$random_seeds), expected_random_seeds) ||
      !identical(as.integer(configuration$bootstrap_reps), 1000L) ||
      !identical(as.integer(configuration$benchmark_reps), 1000L) ||
      !identical(as.integer(configuration$calibration_n), 200000L) ||
      !isTRUE(all.equal(configuration$rho_grid, seq(-0.3, 0.3, by = 0.1))) ||
      configuration$n_genes != 6362L ||
      configuration$n_fash_pairs != 1009173L) {
    stop("The internal correlation configuration is incomplete or unexpected.")
  }

  fit_info <- file.info(configuration$fit_path)
  if (!file.exists(configuration$fit_path) ||
      unname(fit_info$size) != configuration$fit_size ||
      as.character(fit_info$mtime) != configuration$fit_mtime) {
    stop("The source FASH fit no longer matches the cached exploration.")
  }

  required_rows <- c(
    selection_counts = 153L,
    estimator_diagnostics = 153L,
    correlation_matrices_long = 3072L,
    second_moments_long = 1536L,
    adjacent_correlations = 180L,
    lag_summaries = 180L,
    bootstrap_intervals = 60L,
    random_variant_stability = 150L,
    random_lag_stability = 2250L,
    independence_benchmark = 18L,
    truncation_calibration = 21L
  )
  observed_rows <- c(
    selection_counts = nrow(selection_counts),
    estimator_diagnostics = nrow(estimator_diagnostics),
    correlation_matrices_long = nrow(correlation_matrices_long),
    second_moments_long = nrow(second_moments_long),
    adjacent_correlations = nrow(adjacent_correlations),
    lag_summaries = nrow(lag_summaries),
    bootstrap_intervals = nrow(bootstrap_intervals),
    random_variant_stability = nrow(random_variant_stability),
    random_lag_stability = nrow(random_lag_stability),
    independence_benchmark = nrow(independence_benchmark),
    truncation_calibration = nrow(truncation_calibration)
  )
  if (!identical(observed_rows, required_rows)) {
    stop("One or more internal correlation summary tables are incomplete.")
  }

  if (!setequal(unique(selection_counts$threshold), expected_thresholds) ||
      !setequal(random_variant_stability$seed, expected_random_seeds) ||
      any(!is.finite(unlist(list(
        selection_counts$n_selected,
        estimator_diagnostics$sample_minimum_eigenvalue,
        correlation_matrices_long$correlation,
        second_moments_long$value,
        bootstrap_intervals$lower,
        bootstrap_intervals$upper,
        random_variant_stability$mean_lag1_correlation,
        independence_benchmark$mean,
        truncation_calibration$mean_lag1_correlation
      )))) ||
      any(bootstrap_intervals$lower > bootstrap_intervals$upper) ||
      any(estimator_diagnostics$sample_minimum_eigenvalue <= 0) ||
      any(estimator_diagnostics$sample_projection_maximum_change > 1e-8)) {
    stop("The internal correlation cache failed numerical validation.")
  }

  matrix_groups <- split(
    correlation_matrices_long,
    interaction(
      correlation_matrices_long$design,
      correlation_matrices_long$threshold,
      correlation_matrices_long$estimator,
      drop = TRUE
    )
  )
  valid_groups <- vapply(matrix_groups, function(group) {
    matrix <- matrix(
      group$correlation,
      nrow = 16L,
      ncol = 16L
    )
    isTRUE(tryCatch({
      validate_correlation_matrix(matrix)
      TRUE
    }, error = function(unused) FALSE))
  }, logical(1))
  if (!all(valid_groups)) {
    stop("At least one cached full correlation matrix is invalid.")
  }
  invisible(TRUE)
}

validate_internal_correlation_cache()

format_decimal <- function(x, digits = 3) {
  formatC(x, format = "f", digits = digits)
}

format_interval <- function(mean, lower, upper, digits = 3) {
  paste0(
    format_decimal(mean, digits),
    " [",
    format_decimal(lower, digits),
    ", ",
    format_decimal(upper, digits),
    "]"
  )
}

render_scrollable_table <- function(table,
                                    align = NULL,
                                    minimum_width = "900px",
                                    digits = NULL) {
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("The knitr package is required to render tables.")
  }
  kable_arguments <- list(
    x = table,
    format = "html",
    escape = TRUE,
    align = align
  )
  if (!is.null(digits)) {
    kable_arguments$digits <- digits
  }
  rendered_table <- do.call(knitr::kable, kable_arguments)
  cat(
    '<div class="internal-table-scroll">\n',
    '<div style="min-width: ', minimum_width, ';">\n',
    rendered_table,
    '\n</div>\n',
    '</div>\n',
    sep = ""
  )
  invisible(table)
}

prepare_top_lfdr_gallery <- function(analysis,
                                     n_per_quartile = 25L,
                                     seed = 20260891L) {
  metadata <- analysis$top_selection
  z_nullish <- analysis$top_results$threshold_2$z_nullish
  required_columns <- c(
    "rank", "pair_key", "gene_id", "variant_id", "lfdr",
    "max_absolute_z", "passes_threshold_2"
  )
  if (!is.data.frame(metadata) ||
      !all(required_columns %in% names(metadata)) ||
      !is.matrix(z_nullish) || ncol(z_nullish) != 16L ||
      is.null(rownames(z_nullish)) || anyDuplicated(rownames(z_nullish))) {
    stop("The cached top-lfdr data cannot support the trajectory gallery.")
  }

  retained <- metadata[metadata$passes_threshold_2, , drop = FALSE]
  retained <- retained[match(rownames(z_nullish), retained$pair_key), , drop = FALSE]
  observed_maximum <- apply(abs(z_nullish), 1L, max)
  if (nrow(retained) != 752L || anyNA(retained$pair_key) ||
      !identical(retained$pair_key, rownames(z_nullish)) ||
      any(observed_maximum >= 2) ||
      max(abs(observed_maximum - retained$max_absolute_z)) > 1e-10) {
    stop("The retained metadata does not match the exact threshold-2 z matrix.")
  }

  selected <- stratified_lfdr_gallery_sample(
    retained,
    n_per_quartile = n_per_quartile,
    seed = seed
  )
  z_index <- match(selected$pair_key, rownames(z_nullish))
  selected_z <- z_nullish[z_index, , drop = FALSE]
  selected$original_top1000_rank <- selected$rank
  selected$signed_mean_z <- rowMeans(selected_z)
  selected$observed_maximum_absolute_z <- apply(abs(selected_z), 1L, max)
  expected_bounds <- cbind(
    lower = c(1L, 189L, 377L, 565L),
    upper = c(188L, 376L, 564L, 752L)
  )
  valid_quartiles <- vapply(1:4, function(quartile) {
    ranks <- selected$retained_lfdr_rank[
      selected$lfdr_quartile == quartile
    ]
    length(ranks) == n_per_quartile &&
      all(ranks >= expected_bounds[quartile, "lower"]) &&
      all(ranks <= expected_bounds[quartile, "upper"])
  }, logical(1))
  if (nrow(selected) != 4L * n_per_quartile ||
      anyDuplicated(selected$pair_key) || anyNA(z_index) ||
      !all(valid_quartiles) ||
      max(abs(selected$observed_maximum_absolute_z -
        selected$max_absolute_z)) > 1e-10) {
    stop("The stratified trajectory gallery failed validation.")
  }

  list(
    seed = as.integer(seed),
    n_retained = nrow(retained),
    quartile_bounds = expected_bounds,
    metadata = selected,
    z = selected_z
  )
}

top_lfdr_gallery <- prepare_top_lfdr_gallery(internal_analysis)

top_lfdr_gallery_table <- function(quartile) {
  quartile <- as.integer(quartile)
  if (length(quartile) != 1L || is.na(quartile) ||
      !quartile %in% 1:4) {
    stop("quartile must be one integer from 1 through 4.")
  }
  selected <- top_lfdr_gallery$metadata[
    top_lfdr_gallery$metadata$lfdr_quartile == quartile,
    ,
    drop = FALSE
  ]
  data.frame(
    `Page position` = selected$gallery_page_position,
    `Retained lfdr rank` = selected$retained_lfdr_rank,
    `Original top-1,000 rank` = selected$original_top1000_rank,
    Gene = selected$gene_id,
    Variant = selected$variant_id,
    lfdr = format_decimal(selected$lfdr, 6),
    `max|z|` = format_decimal(selected$observed_maximum_absolute_z, 3),
    `Mean z` = format_decimal(selected$signed_mean_z, 3),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

plot_top_lfdr_gallery_page <- function(quartile) {
  quartile <- as.integer(quartile)
  if (length(quartile) != 1L || is.na(quartile) ||
      !quartile %in% 1:4) {
    stop("quartile must be one integer from 1 through 4.")
  }
  keep <- top_lfdr_gallery$metadata$lfdr_quartile == quartile
  metadata <- top_lfdr_gallery$metadata[keep, , drop = FALSE]
  z <- top_lfdr_gallery$z[keep, , drop = FALSE]
  if (nrow(metadata) != 25L || !identical(rownames(z), metadata$pair_key)) {
    stop("The requested trajectory gallery page is incomplete.")
  }

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(5, 5),
    mar = c(1.35, 1.45, 2.45, 0.45),
    oma = c(2.8, 3.0, 3.2, 0.6),
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
      col = "#0072B2",
      xlim = c(0, 15),
      ylim = c(-2.05, 2.05),
      axes = FALSE,
      xlab = "",
      ylab = ""
    )
    graphics::abline(h = c(-2, 2), lty = 3, col = "gray70")
    graphics::abline(h = 0, lty = 1, col = "gray82")
    graphics::lines(
      time,
      z[panel, ],
      type = "o",
      pch = 16,
      cex = 0.48,
      lwd = 1.0,
      col = "#0072B2"
    )
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
    panel_title <- paste0(
      "R", metadata$retained_lfdr_rank[panel], " | ",
      metadata$gene_id[panel], "\n",
      metadata$variant_id[panel],
      " | lfdr=", format_decimal(metadata$lfdr[panel], 4),
      " | max=", format_decimal(
        metadata$observed_maximum_absolute_z[panel],
        2
      )
    )
    graphics::title(main = panel_title, cex.main = 0.47, line = 0.45)
  }
  bounds <- top_lfdr_gallery$quartile_bounds[quartile, ]
  graphics::mtext(
    paste0(
      "Retained lfdr quartile ", quartile,
      ": ranks ", bounds["lower"], "-", bounds["upper"],
      "; 25 fixed-seed samples"
    ),
    side = 3,
    outer = TRUE,
    line = 1.2,
    cex = 1.05,
    font = 2
  )
  graphics::mtext("Time", side = 1, outer = TRUE, line = 1.4, cex = 0.9)
  graphics::mtext(
    "z = beta-hat / adjusted SE",
    side = 2,
    outer = TRUE,
    line = 1.5,
    cex = 0.9
  )
  invisible(metadata)
}

get_matrix_from_long <- function(table,
                                 design,
                                 threshold,
                                 estimator,
                                 value_column = "correlation") {
  subset <- table[
    table$design == design &
      table$threshold == threshold &
      table$estimator == estimator,
    ,
    drop = FALSE
  ]
  if (nrow(subset) != 256L) {
    stop("Could not recover one complete 16 by 16 matrix.")
  }
  matrix(
    subset[[value_column]],
    nrow = 16L,
    ncol = 16L
  )
}

get_primary_lag <- function(design, threshold, lag) {
  hit <- lag_summaries[
    lag_summaries$design == design &
      lag_summaries$threshold == threshold &
      lag_summaries$estimator == "mashr cor(z)" &
      lag_summaries$lag == lag,
    "mean_correlation"
  ]
  if (length(hit) != 1L) {
    stop("Could not find one requested lag summary.")
  }
  hit
}

get_bootstrap_interval <- function(design, summary_type, index) {
  hit <- bootstrap_intervals[
    bootstrap_intervals$design == design &
      bootstrap_intervals$summary_type == summary_type &
      bootstrap_intervals$index == index,
    ,
    drop = FALSE
  ]
  if (nrow(hit) != 1L) {
    stop("Could not find one requested bootstrap interval.")
  }
  hit
}

summarize_random_selection_lags <- function(threshold = 2) {
  subset <- random_lag_stability[random_lag_stability$threshold == threshold, ]
  do.call(rbind, lapply(split(subset, subset$lag), function(group) {
    data.frame(
      lag = group$lag[1],
      mean = mean(group$mean_correlation),
      lower = stats::quantile(group$mean_correlation, 0.025),
      upper = stats::quantile(group$mean_correlation, 0.975),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
}

top_threshold2_count <- selection_counts[
  selection_counts$design == top_design & selection_counts$threshold == 2,
]
random_threshold2_counts <- selection_counts[
  selection_counts$design == "Random variant per gene before z" &
    selection_counts$threshold == 2,
]
top_threshold2_diag <- estimator_diagnostics[
  estimator_diagnostics$design == top_design &
    estimator_diagnostics$threshold == 2,
]
random_threshold2_diag <- estimator_diagnostics[
  estimator_diagnostics$design == "Random variant per gene before z" &
    estimator_diagnostics$threshold == 2,
]

top_day01 <- adjacent_correlations$correlation[
  adjacent_correlations$design == top_design &
    adjacent_correlations$threshold == 2 &
    adjacent_correlations$estimator == "mashr cor(z)" &
    adjacent_correlations$time_a == 0
]
random_day01 <- adjacent_correlations$correlation[
  adjacent_correlations$design == random_design &
    adjacent_correlations$threshold == 2 &
    adjacent_correlations$estimator == "mashr cor(z)" &
    adjacent_correlations$time_a == 0
]
top_lag1_bootstrap <- get_bootstrap_interval(
  top_design,
  "lag_average",
  1L
)
primary_random_lag1_bootstrap <- get_bootstrap_interval(
  primary_random_design,
  "lag_average",
  1L
)
random_lag1_selection <- random_variant_stability[
  random_variant_stability$threshold == 2,
  "mean_lag1_correlation"
]

primary_results_table <- data.frame(
  Design = c(
    "Top-1,000 highest lfdr per gene",
    "Pre-z random variant per gene (50-selection mean)"
  ),
  `Retained at max|z| < 2` = c(
    as.character(top_threshold2_count$n_selected),
    paste0(
      format_decimal(mean(random_threshold2_counts$n_selected), 1),
      " [",
      min(random_threshold2_counts$n_selected),
      ", ",
      max(random_threshold2_counts$n_selected),
      "]"
    )
  ),
  `Day 0-1 correlation` = format_decimal(c(top_day01, random_day01)),
  `Mean lag-1 correlation` = format_decimal(c(
    get_primary_lag(top_design, 2, 1),
    get_primary_lag(random_design, 2, 1)
  )),
  `Lag-1 interval (source)` = c(
    paste0(
      format_interval(
        top_lag1_bootstrap$mean,
        top_lag1_bootstrap$lower,
        top_lag1_bootstrap$upper
      ),
      " (gene bootstrap)"
    ),
    paste0(
      format_interval(
        mean(random_lag1_selection),
        stats::quantile(random_lag1_selection, 0.025),
        stats::quantile(random_lag1_selection, 0.975)
      ),
      " (50 selections)"
    )
  ),
  `Lag-15 correlation` = format_decimal(c(
    get_primary_lag(top_design, 2, 15),
    get_primary_lag(random_design, 2, 15)
  )),
  `Max |column mean z|` = format_decimal(c(
    top_threshold2_diag$maximum_absolute_column_mean,
    mean(random_threshold2_diag$maximum_absolute_column_mean)
  )),
  `Minimum eigenvalue` = format_decimal(c(
    top_threshold2_diag$sample_minimum_eigenvalue,
    mean(random_threshold2_diag$sample_minimum_eigenvalue)
  )),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

threshold_table_rows <- list()
for (threshold in expected_thresholds) {
  top_count <- selection_counts[
    selection_counts$design == top_design &
      selection_counts$threshold == threshold,
  ]
  random_counts <- selection_counts[
    selection_counts$design == "Random variant per gene before z" &
      selection_counts$threshold == threshold,
  ]
  random_stability <- random_variant_stability[
    random_variant_stability$threshold == threshold,
  ]
  benchmark <- independence_benchmark[
    independence_benchmark$candidate_design ==
      "One-candidate-per-gene count" &
      independence_benchmark$threshold == threshold &
      independence_benchmark$statistic == "mean_lag1",
  ]
  threshold_table_rows[[length(threshold_table_rows) + 1L]] <- data.frame(
    Threshold = threshold,
    `Top-lfdr retained` = top_count$n_selected,
    `Top-lfdr lag 1` = get_primary_lag(top_design, threshold, 1),
    `Random retained, mean [range]` = paste0(
      format_decimal(mean(random_counts$n_selected), 1),
      " [", min(random_counts$n_selected), ", ",
      max(random_counts$n_selected), "]"
    ),
    `Random lag 1, mean [range]` = paste0(
      format_decimal(mean(random_stability$mean_lag1_correlation)),
      " [",
      format_decimal(min(random_stability$mean_lag1_correlation)),
      ", ",
      format_decimal(max(random_stability$mean_lag1_correlation)),
      "]"
    ),
    `Independent lag-1 95% interval` = paste0(
      "[",
      format_decimal(benchmark$lower),
      ", ",
      format_decimal(benchmark$upper),
      "]"
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
threshold_sensitivity_table <- do.call(rbind, threshold_table_rows)

estimator_comparison_rows <- list()
for (design in c(top_design, random_design)) {
  sample_matrix <- get_matrix_from_long(
    correlation_matrices_long,
    design,
    2,
    "mashr cor(z)"
  )
  normalized_matrix <- get_matrix_from_long(
    correlation_matrices_long,
    design,
    2,
    "normalized second moment"
  )
  second_moment <- get_matrix_from_long(
    second_moments_long,
    design,
    2,
    "literal nullish-z second moment",
    value_column = "value"
  )
  estimator_comparison_rows[[length(estimator_comparison_rows) + 1L]] <-
    data.frame(
      Design = design,
      `cor(z), day 0-1` = sample_matrix[1, 2],
      `Normalized second moment, day 0-1` = normalized_matrix[1, 2],
      `Literal second moment, day 0-1` = second_moment[1, 2],
      `Second-moment diagonal range` = paste0(
        "[",
        format_decimal(min(diag(second_moment))),
        ", ",
        format_decimal(max(diag(second_moment))),
        "]"
      ),
      `Max correlation-scale difference` = max(abs(
        sample_matrix - normalized_matrix
      )),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
}
estimator_comparison_table <- do.call(rbind, estimator_comparison_rows)

calibration_table <- truncation_calibration[
  truncation_calibration$threshold == 2,
  c("generating_rho", "selected_fraction", "mean_lag1_correlation")
]
names(calibration_table) <- c(
  "Generating lag-1 correlation",
  "Fraction retained",
  "Observed lag-1 correlation after max|z| < 2"
)
calibration_rho03_threshold2 <- truncation_calibration$mean_lag1_correlation[
  truncation_calibration$threshold == 2 &
    abs(truncation_calibration$generating_rho - 0.3) < 1e-8
]
if (length(calibration_rho03_threshold2) != 1L) {
  stop("Could not identify the rho=0.3, threshold-2 calibration result.")
}

plot_matrix_panel <- function(correlation, title) {
  palette <- grDevices::colorRampPalette(c("#2166AC", "white", "#B2182B"))(201)
  graphics::image(
    x = 0:15,
    y = 0:15,
    z = correlation,
    zlim = c(-1, 1),
    col = palette,
    axes = FALSE,
    xlab = "Time",
    ylab = "Time",
    main = title,
    cex.main = 0.95,
    asp = 1
  )
  graphics::axis(1, at = 0:15, labels = 0:15, cex.axis = 0.75)
  graphics::axis(2, at = 0:15, labels = 0:15, las = 1, cex.axis = 0.75)
  graphics::box()
}

plot_internal_correlation_heatmaps <- function() {
  top_matrix <- get_matrix_from_long(
    correlation_matrices_long,
    top_design,
    2,
    "mashr cor(z)"
  )
  random_matrix <- get_matrix_from_long(
    correlation_matrices_long,
    random_design,
    2,
    "mashr cor(z)"
  )
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(
    mfrow = c(1, 2),
    mar = c(4.2, 4.2, 4.2, 1.2),
    oma = c(1.8, 0, 3.0, 0)
  )
  plot_matrix_panel(top_matrix, "Top-1,000 lfdr, then max|z| < 2")
  plot_matrix_panel(
    random_matrix,
    "Pre-z random per gene\n(50-selection mean)"
  )
  graphics::mtext(
    "Internal zero-intercept z-score correlations",
    side = 3,
    outer = TRUE,
    line = 1.1,
    cex = 1.45,
    font = 2
  )
  graphics::mtext(
    "Common scale: blue = -1, white = 0, red = +1; both matrices are already positive definite",
    side = 1,
    outer = TRUE,
    line = 0.2,
    cex = 0.88
  )
}

plot_internal_adjacent_correlations <- function() {
  top <- adjacent_correlations[
    adjacent_correlations$design == top_design &
      adjacent_correlations$threshold == 2 &
      adjacent_correlations$estimator == "mashr cor(z)",
  ]
  random <- adjacent_correlations[
    adjacent_correlations$design == random_design &
      adjacent_correlations$threshold == 2 &
      adjacent_correlations$estimator == "mashr cor(z)",
  ]
  top_boot <- bootstrap_intervals[
    bootstrap_intervals$design == top_design &
      bootstrap_intervals$summary_type == "adjacent_pair",
  ]
  random_array <- internal_analysis$random_mean_results[[
    "threshold_2"
  ]]$sample_array
  random_adjacent <- vapply(seq_len(15L), function(index) {
    random_array[index, index + 1L, ]
  }, numeric(dim(random_array)[3]))
  random_lower <- apply(random_adjacent, 2L, stats::quantile, 0.025)
  random_upper <- apply(random_adjacent, 2L, stats::quantile, 0.975)
  y_limits <- range(c(
    top_boot$lower,
    top_boot$upper,
    random_lower,
    random_upper,
    0
  ))
  y_limits[2] <- y_limits[2] + 0.1
  graphics::plot(
    1:15,
    top$correlation,
    type = "b",
    pch = 19,
    lwd = 2,
    col = "#0072B2",
    ylim = y_limits,
    xaxt = "n",
    xlab = "Adjacent time pair",
    ylab = "Correlation",
    main = "Adjacent correlations at max|z| < 2"
  )
  graphics::axis(1, at = 1:15, labels = paste0(0:14, "-", 1:15), las = 2)
  graphics::abline(h = 0, lty = 3, col = "gray45")
  graphics::arrows(
    1:15 - 0.06,
    top_boot$lower,
    1:15 - 0.06,
    top_boot$upper,
    angle = 90,
    code = 3,
    length = 0.03,
    col = grDevices::adjustcolor("#0072B2", 0.7)
  )
  graphics::lines(
    1:15 + 0.06,
    random$correlation,
    type = "b",
    pch = 1,
    lwd = 2,
    col = "#D55E00"
  )
  graphics::arrows(
    1:15 + 0.06,
    random_lower,
    1:15 + 0.06,
    random_upper,
    angle = 90,
    code = 3,
    length = 0.03,
    col = grDevices::adjustcolor("#D55E00", 0.7)
  )
  graphics::legend(
    "topleft",
    legend = c(
      "Top-1,000 lfdr (95% gene bootstrap)",
      "Pre-z random, 50-selection interval"
    ),
    col = c("#0072B2", "#D55E00"),
    pch = c(19, 1),
    lwd = 2,
    bty = "n"
  )
}

plot_internal_lag_diagnostics <- function() {
  top <- lag_summaries[
    lag_summaries$design == top_design &
      lag_summaries$threshold == 2 &
      lag_summaries$estimator == "mashr cor(z)",
  ]
  random <- summarize_random_selection_lags(2)
  top_boot <- bootstrap_intervals[
    bootstrap_intervals$design == top_design &
      bootstrap_intervals$summary_type == "lag_average",
  ]
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.4, 3.7, 1.0))

  y_limits <- range(c(top_boot$lower, top_boot$upper, random$lower, random$upper, 0))
  graphics::plot(
    top$lag,
    top$mean_correlation,
    type = "n",
    ylim = y_limits,
    xlab = "Time lag",
    ylab = "Mean correlation",
    main = "Lag correlation"
  )
  graphics::polygon(
    c(top$lag, rev(top$lag)),
    c(top_boot$lower, rev(top_boot$upper)),
    col = grDevices::adjustcolor("#0072B2", 0.14),
    border = NA
  )
  graphics::polygon(
    c(random$lag, rev(random$lag)),
    c(random$lower, rev(random$upper)),
    col = grDevices::adjustcolor("#D55E00", 0.14),
    border = NA
  )
  graphics::lines(top$lag, top$mean_correlation, type = "b", pch = 19, lwd = 2, col = "#0072B2")
  graphics::lines(random$lag, random$mean, type = "b", pch = 1, lwd = 2, col = "#D55E00")
  graphics::abline(h = 0, lty = 3, col = "gray45")
  graphics::legend(
    "topright",
    legend = c("Top-1,000 lfdr", "Pre-z random"),
    col = c("#0072B2", "#D55E00"),
    pch = c(19, 1),
    lwd = 2,
    bty = "n"
  )

  gamma_top <- 1 - top$mean_correlation
  gamma_top_lower <- 1 - top_boot$upper
  gamma_top_upper <- 1 - top_boot$lower
  gamma_random <- 1 - random$mean
  gamma_random_lower <- 1 - random$upper
  gamma_random_upper <- 1 - random$lower
  gamma_limits <- range(c(
    gamma_top_lower,
    gamma_top_upper,
    gamma_random_lower,
    gamma_random_upper,
    1
  ))
  graphics::plot(
    top$lag,
    gamma_top,
    type = "n",
    ylim = gamma_limits,
    xlab = "Time lag",
    ylab = "Standardized semivariogram",
    main = "Variogram"
  )
  graphics::polygon(
    c(top$lag, rev(top$lag)),
    c(gamma_top_lower, rev(gamma_top_upper)),
    col = grDevices::adjustcolor("#0072B2", 0.14),
    border = NA
  )
  graphics::polygon(
    c(random$lag, rev(random$lag)),
    c(gamma_random_lower, rev(gamma_random_upper)),
    col = grDevices::adjustcolor("#D55E00", 0.14),
    border = NA
  )
  graphics::lines(top$lag, gamma_top, type = "b", pch = 19, lwd = 2, col = "#0072B2")
  graphics::lines(random$lag, gamma_random, type = "b", pch = 1, lwd = 2, col = "#D55E00")
  graphics::abline(h = 1, lty = 3, col = "gray45")
}

plot_internal_threshold_sensitivity <- function() {
  top <- lag_summaries[
    lag_summaries$design == top_design &
      lag_summaries$estimator == "mashr cor(z)" &
      lag_summaries$lag == 1,
  ]
  random <- do.call(rbind, lapply(split(
    random_variant_stability,
    random_variant_stability$threshold
  ), function(group) {
    data.frame(
      threshold = group$threshold[1],
      mean = mean(group$mean_lag1_correlation),
      lower = stats::quantile(group$mean_lag1_correlation, 0.025),
      upper = stats::quantile(group$mean_lag1_correlation, 0.975),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }))
  benchmark <- independence_benchmark[
    independence_benchmark$candidate_design == "One-candidate-per-gene count" &
      independence_benchmark$statistic == "mean_lag1",
  ]
  y_limits <- range(c(
    top$mean_correlation,
    random$lower,
    random$upper,
    benchmark$lower,
    benchmark$upper
  ))
  graphics::plot(
    top$threshold,
    top$mean_correlation,
    type = "b",
    pch = 19,
    lwd = 2,
    col = "#0072B2",
    ylim = y_limits,
    xlab = "max|z| threshold",
    ylab = "Mean lag-1 correlation",
    main = "Threshold and preselection sensitivity"
  )
  graphics::polygon(
    c(random$threshold, rev(random$threshold)),
    c(random$lower, rev(random$upper)),
    col = grDevices::adjustcolor("#D55E00", 0.14),
    border = NA
  )
  graphics::lines(random$threshold, random$mean, type = "b", pch = 1, lwd = 2, col = "#D55E00")
  graphics::polygon(
    c(benchmark$threshold, rev(benchmark$threshold)),
    c(benchmark$lower, rev(benchmark$upper)),
    col = grDevices::adjustcolor("gray40", 0.12),
    border = NA
  )
  graphics::abline(h = 0, lty = 3, col = "gray40")
  graphics::legend(
    "topleft",
    legend = c(
      "Top-1,000 lfdr",
      "Pre-z random, 50-selection interval",
      "Independent truncation benchmark"
    ),
    col = c("#0072B2", "#D55E00", "gray45"),
    pch = c(19, 1, NA),
    lty = c(1, 1, 3),
    lwd = c(2, 2, 1.5),
    bty = "n"
  )
}

plot_internal_truncation_calibration <- function() {
  colors <- c("1.5" = "#009E73", "2" = "#0072B2", "2.5" = "#D55E00")
  graphics::plot(
    range(truncation_calibration$generating_rho),
    range(truncation_calibration$mean_lag1_correlation),
    type = "n",
    xlab = "Generating lag-1 correlation",
    ylab = "Observed correlation after truncation",
    main = "Rectangular z-filter attenuation calibration"
  )
  graphics::abline(a = 0, b = 1, lty = 3, col = "gray45")
  for (threshold in expected_thresholds) {
    subset <- truncation_calibration[
      truncation_calibration$threshold == threshold,
    ]
    graphics::lines(
      subset$generating_rho,
      subset$mean_lag1_correlation,
      type = "b",
      pch = 19,
      lwd = 2,
      col = colors[as.character(threshold)]
    )
  }
  graphics::legend(
    "topleft",
    legend = paste0("max|z| < ", expected_thresholds),
    col = colors,
    pch = 19,
    lwd = 2,
    bty = "n"
  )
}
