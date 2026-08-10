# Load, validate, and format cached results for the R3 workflowr report.

mc_output_id <-
  "r3_matched_functional_relative_clearance_main_effect_pilot5"
mc_output_dir <- file.path(
  "output",
  "revision_simulations",
  "mc",
  mc_output_id
)
summary_dir <- file.path(mc_output_dir, "summary")
configuration <- readRDS(file.path(mc_output_dir, "configuration.rds"))
example_curves <- readRDS(file.path(
  mc_output_dir,
  "example_curves.rds"
))
mc_alpha <- read.csv(
  file.path(summary_dir, "functional_testing_mc_alpha_curve.csv"),
  stringsAsFactors = FALSE
)
mc_alpha_005 <- read.csv(
  file.path(
    summary_dir,
    "functional_testing_mc_alpha005_summary.csv"
  ),
  stringsAsFactors = FALSE
)
mc_pi0 <- read.csv(
  file.path(summary_dir, "functional_testing_mc_pi0_summary.csv"),
  stringsAsFactors = FALSE
)

method_order <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
target_order <- c("early", "middle", "late", "switch")
mechanism_order <- c("random_bspline", "raised_cosine")
mechanism_labels <- c(
  random_bspline = "R3A: broad random B-spline",
  raised_cosine = "R3B: compact raised cosine"
)

if (!isTRUE(all.equal(configuration$J, 1000L)) ||
    !isTRUE(all.equal(configuration$n_donors, 19L)) ||
    !isTRUE(all.equal(configuration$n_covariates, 5L)) ||
    !isTRUE(all.equal(configuration$time_grid, 0:15)) ||
    !isTRUE(all.equal(configuration$expression_noise_sd, 1)) ||
    !isTRUE(all.equal(configuration$maf_range, c(0.1, 0.5))) ||
    !isTRUE(all.equal(configuration$covariate_effect_sd, 0.5)) ||
    !isTRUE(all.equal(configuration$intercept_sd, 0)) ||
    !isTRUE(all.equal(configuration$dynamic_main_effect_sd, 1)) ||
    !identical(configuration$truth_mechanisms, mechanism_order) ||
    !isTRUE(all.equal(configuration$random_bspline$amplitude, 2)) ||
    !isTRUE(all.equal(configuration$random_bspline$df, 6)) ||
    !isTRUE(all.equal(configuration$raised_cosine$width_half, 1.5)) ||
    !isTRUE(all.equal(
      configuration$raised_cosine$spike_counts,
      1:3
    )) ||
    !isTRUE(all.equal(configuration$switch_threshold, 0.25)) ||
    !isTRUE(all.equal(configuration$location_truth_margin, 0.10)) ||
    !isTRUE(all.equal(configuration$switch_truth_margin, 0.10)) ||
    !isTRUE(all.equal(configuration$non_switch_min_abs, 0.10)) ||
    !isTRUE(all.equal(
      configuration$non_switch_min_range_fraction,
      0.10
    )) ||
    length(configuration$seed_list) != 5L ||
    !identical(names(example_curves), mechanism_order)) {
  stop("The matched functional-testing cache has unexpected settings.")
}
if (!all(method_order %in% mc_alpha$method) ||
    !all(target_order %in% mc_alpha$target) ||
    any(mc_alpha$n_replications != length(configuration$seed_list)) ||
    !all(mechanism_order %in% mc_pi0$truth_mechanism)) {
  stop("The matched functional-testing summaries are incomplete.")
}

r3a_alpha <- mc_alpha[grepl("^r3a_", mc_alpha$scenario), ]
r3b_alpha <- mc_alpha[grepl("^r3b_", mc_alpha$scenario), ]
r3a_alpha_005 <- mc_alpha_005[
  grepl("^r3a_", mc_alpha_005$scenario),
]
r3b_alpha_005 <- mc_alpha_005[
  grepl("^r3b_", mc_alpha_005$scenario),
]

format_decimal <- function(x, digits = 3) {
  formatC(x, format = "f", digits = digits)
}

format_mc_interval <- function(mean, lower, upper, digits = 3) {
  paste0(
    format_decimal(mean, digits),
    " [",
    format_decimal(lower, digits),
    ", ",
    format_decimal(upper, digits),
    "]"
  )
}

format_functional_table <- function(summary_table) {
  table <- summary_table
  table$target <- factor(table$target, levels = target_order)
  table$method <- factor(table$method, levels = method_order)
  table <- table[order(table$target, table$method), , drop = FALSE]
  data.frame(
    Target = tools::toTitleCase(as.character(table$target)),
    Method = as.character(table$method),
    `Mean calls` = format_decimal(table$mean_discoveries, 1),
    `Power (95% MC CI)` = format_mc_interval(
      table$mean_power,
      table$power_ci_lower,
      table$power_ci_upper
    ),
    `Empirical FSR (95% MC CI)` = format_mc_interval(
      table$mean_empirical_fsr,
      table$empirical_fsr_ci_lower,
      table$empirical_fsr_ci_upper
    ),
    check.names = FALSE
  )
}

format_pi0_table <- function(summary_table) {
  summary_table$truth_mechanism <- factor(
    summary_table$truth_mechanism,
    levels = mechanism_order
  )
  summary_table$fit <- factor(
    summary_table$fit,
    levels = c("Raw", "BF-corrected")
  )
  summary_table <- summary_table[
    order(summary_table$truth_mechanism, summary_table$fit),
    ,
    drop = FALSE
  ]
  data.frame(
    Mechanism = unname(mechanism_labels[
      as.character(summary_table$truth_mechanism)
    ]),
    Fit = as.character(summary_table$fit),
    `Mean estimated pi0 (95% MC CI)` = format_mc_interval(
      summary_table$mean_estimated_pi0,
      summary_table$pi0_ci_lower,
      summary_table$pi0_ci_upper
    ),
    check.names = FALSE
  )
}

metric_at <- function(summary_table, target, method, metric) {
  row <- summary_table[
    summary_table$target == target &
      summary_table$method == method,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L || !metric %in% names(row)) {
    stop("Could not extract a unique functional-testing metric.")
  }
  row[[metric]]
}

plot_truth_examples <- function(examples, mechanism_label) {
  group_order <- c(
    "early / switch",
    "early / non-switch",
    "middle / switch",
    "middle / non-switch",
    "late / switch",
    "late / non-switch"
  )
  panels <- unlist(examples[group_order], recursive = FALSE)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(
    mfrow = c(6, 2),
    mar = c(3.5, 3.8, 3.0, 0.8),
    oma = c(0, 0, 2.2, 0)
  )
  for (panel_index in seq_along(panels)) {
    panel <- panels[[panel_index]]
    observed <- panel$observed
    target_window <- switch(
      panel$time_group,
      early = c(0, 3),
      middle = c(4, 11),
      late = c(12, 15)
    )
    window_color <- switch(
      panel$time_group,
      early = "#d7eaf7",
      middle = "#dff0d8",
      late = "#fde6c9"
    )
    y_limits <- range(
      observed$estimate - 2 * observed$se,
      observed$estimate + 2 * observed$se,
      panel$true_curve$true_effect
    )
    y_padding <- max(0.15, 0.06 * diff(y_limits))
    y_limits <- y_limits + c(-y_padding, y_padding)
    shape_suffix <- if (is.na(panel$spike_count)) {
      ""
    } else {
      paste0("; ", panel$spike_count, " peak",
        if (panel$spike_count == 1L) "" else "s"
      )
    }
    plot(
      observed$time,
      observed$estimate,
      type = "n",
      xlim = range(configuration$time_grid),
      ylim = y_limits,
      xlab = if (panel_index > length(panels) - 2L) "Time" else "",
      ylab = "Genetic effect",
      main = paste0(
        panel$truth_group,
        shape_suffix,
        "\n",
        panel$variant_id
      )
    )
    rect(
      target_window[1],
      y_limits[1],
      target_window[2],
      y_limits[2],
      col = window_color,
      border = NA
    )
    arrows(
      observed$time,
      observed$estimate - 2 * observed$se,
      observed$time,
      observed$estimate + 2 * observed$se,
      angle = 90,
      code = 3,
      length = 0.035,
      col = "gray45"
    )
    points(observed$time, observed$estimate, pch = 19)
    lines(
      panel$true_curve$time,
      panel$true_curve$true_effect,
      col = "#D55E00",
      lwd = 2.2
    )
    abline(h = 0, col = "gray75", lty = 3)
  }
  mtext(
    paste0(
      mechanism_label,
      ": estimates with two-SE bars and continuous truth"
    ),
    outer = TRUE,
    cex = 1.05
  )
}

r3a_table <- format_functional_table(r3a_alpha_005)
r3b_table <- format_functional_table(r3b_alpha_005)
pi0_table <- format_pi0_table(mc_pi0)
