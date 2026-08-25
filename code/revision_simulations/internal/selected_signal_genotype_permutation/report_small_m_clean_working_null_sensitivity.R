# Reporting helpers for the small-M clean working-model null sensitivity page.

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

workflowr_root <- find_workflowr_root()
clean_output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  paste0(
    "all_gene_random_variant_small_m_clean_working_null_",
    "selection20260817_nullseed20260819"
  )
)
residual_output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  paste0(
    "all_gene_random_variant_small_m_signal_stripped_residual_sensitivity_",
    "selection20260817_seed20260811_subsets20260819"
  )
)

required_clean_files <- c(
  "configuration.rds", "alpha005.csv", "alpha_curves.csv",
  "pi0_summary.csv", "prior_weights.csv", "null_z_overall.csv",
  "null_z_by_time.csv", "null_z_correlations.csv", "runtime.csv",
  "validation.csv"
)
required_residual_files <- c("replicate_alpha005.csv")
if (any(!file.exists(file.path(clean_output_directory, required_clean_files))) ||
    any(!file.exists(
      file.path(residual_output_directory, required_residual_files)
    ))) {
  stop("At least one clean-null reporting input is missing.")
}

clean_configuration <- readRDS(file.path(
  clean_output_directory,
  "configuration.rds"
))
clean_alpha005 <- utils::read.csv(file.path(
  clean_output_directory,
  "alpha005.csv"
), stringsAsFactors = FALSE)
clean_alpha_curves <- utils::read.csv(file.path(
  clean_output_directory,
  "alpha_curves.csv"
), stringsAsFactors = FALSE)
clean_pi0 <- utils::read.csv(file.path(
  clean_output_directory,
  "pi0_summary.csv"
), stringsAsFactors = FALSE)
clean_prior_weights <- utils::read.csv(file.path(
  clean_output_directory,
  "prior_weights.csv"
), stringsAsFactors = FALSE)
clean_null_z_overall <- utils::read.csv(file.path(
  clean_output_directory,
  "null_z_overall.csv"
), stringsAsFactors = FALSE)
clean_null_z_by_time <- utils::read.csv(file.path(
  clean_output_directory,
  "null_z_by_time.csv"
), stringsAsFactors = FALSE)
clean_null_z_correlations <- utils::read.csv(file.path(
  clean_output_directory,
  "null_z_correlations.csv"
), stringsAsFactors = FALSE)
clean_runtime <- utils::read.csv(file.path(
  clean_output_directory,
  "runtime.csv"
), stringsAsFactors = FALSE)
clean_validation <- utils::read.csv(file.path(
  clean_output_directory,
  "validation.csv"
), stringsAsFactors = FALSE)
residual_alpha005 <- utils::read.csv(file.path(
  residual_output_directory,
  "replicate_alpha005.csv"
), stringsAsFactors = FALSE)

if (!all(clean_validation$passed) ||
    clean_configuration$J != 6362L ||
    clean_configuration$n_simulation_draws != 1L ||
    nrow(clean_alpha005) != 8L ||
    nrow(clean_alpha_curves) != 1600L ||
    !all(c("Raw", "BF-adjusted") %in% clean_alpha005$fit_stage) ||
    !all(c(0.05, 0.10, 0.20, 1.00) %in% clean_alpha005$m_ratio)) {
  stop("The clean-null reporting inputs failed validation.")
}

if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("scales", quietly = TRUE) ||
    !requireNamespace("knitr", quietly = TRUE)) {
  stop("ggplot2, scales, and knitr are required for reporting.")
}

format_integer <- function(value) {
  format(round(value), big.mark = ",", scientific = FALSE, trim = TRUE)
}

format_decimal <- function(value, digits = 3L) {
  formatC(value, format = "f", digits = digits)
}

format_percent <- function(value, digits = 1L) {
  ifelse(
    is.finite(value),
    paste0(formatC(100 * value, format = "f", digits = digits), "%"),
    "NA"
  )
}

format_interval <- function(values, digits = 2L) {
  quantiles <- stats::quantile(
    values,
    probs = c(0.25, 0.50, 0.75),
    names = FALSE,
    type = 8
  )
  paste0(
    format_percent(quantiles[2], digits),
    " [",
    format_percent(quantiles[1], digits),
    ", ",
    format_percent(quantiles[3], digits),
    "]"
  )
}

render_clean_table <- function(value, align = NULL) {
  print(knitr::kable(
    value,
    format = "html",
    align = align,
    escape = TRUE,
    row.names = FALSE,
    table.attr = 'class="table table-striped table-condensed"'
  ))
}

build_alpha005_table <- function(fit_stage) {
  current <- clean_alpha005[clean_alpha005$fit_stage == fit_stage, ]
  current <- current[order(current$m_ratio), ]
  data.frame(
    `M / J` = format_decimal(current$m_ratio, 2L),
    M = format_integer(current$m_size),
    `Merged pi0` = format_decimal(current$pi0_merged, 3L),
    `Recovered target pi0` = format_decimal(
      current$pi0_target_unbounded,
      3L
    ),
    `Target calls` = format_integer(current$target_calls),
    V = format_integer(current$permuted_null_calls),
    `R merged` = format_integer(current$merged_calls),
    `V / R merged` = format_percent(
      current$known_null_discovery_fraction,
      2L
    ),
    `V / M` = format_percent(current$permuted_null_call_rate, 3L),
    `Target-pi0 plug-in FDR` = format_percent(
      current$target_pi0_plugin_fdr,
      1L
    ),
    `Merged-pi0 plug-in FDR` = format_percent(
      current$merged_pi0_plugin_fdr,
      1L
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

clean_bf_alpha005_table <- build_alpha005_table("BF-adjusted")
clean_raw_alpha005_table <- build_alpha005_table("Raw")

build_residual_comparison_table <- function() {
  residual <- residual_alpha005[
    residual_alpha005$fit_stage == "BF-adjusted" &
      residual_alpha005$m_ratio %in% c(0.05, 0.10, 0.20, 1.00),
  ]
  clean <- clean_alpha005[clean_alpha005$fit_stage == "BF-adjusted", ]
  rows <- lapply(c(0.05, 0.10, 0.20, 1.00), function(ratio) {
    residual_current <- residual[abs(residual$m_ratio - ratio) < 1e-12, ]
    clean_current <- clean[abs(clean$m_ratio - ratio) < 1e-12, ]
    data.frame(
      `M / J` = format_decimal(ratio, 2L),
      `Clean-null V` = format_integer(clean_current$permuted_null_calls),
      `Clean-null V / R merged` = format_percent(
        clean_current$known_null_discovery_fraction,
        2L
      ),
      `Synchronized-residual V / R merged` = if (ratio < 1) {
        format_interval(residual_current$known_null_discovery_fraction, 2L)
      } else {
        format_percent(residual_current$known_null_discovery_fraction, 2L)
      },
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
  do.call(rbind, rows)
}

clean_residual_comparison_table <- build_residual_comparison_table()

ratio_colors <- c(
  "M/J = 0.05" = "#1B9E77",
  "M/J = 0.10" = "#D95F02",
  "M/J = 0.20" = "#7570B3",
  "M/J = 1.00" = "#333333"
)

make_ratio_label <- function(value) {
  factor(
    sprintf("M/J = %.2f", value),
    levels = names(ratio_colors)
  )
}

clean_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 10.2, color = "#4D4D4D"),
      axis.title = ggplot2::element_text(size = 11),
      axis.text = ggplot2::element_text(size = 9.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.justification = "left",
      strip.text = ggplot2::element_text(face = "bold")
    )
}

plot_clean_known_null_fdp <- function() {
  plot_data <- clean_alpha_curves[
    clean_alpha_curves$fit_stage == "BF-adjusted",
  ]
  plot_data$ratio_label <- make_ratio_label(plot_data$m_ratio)
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = nominal_alpha,
      y = known_null_discovery_fraction,
      color = ratio_label
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      color = "#888888",
      linetype = "dotted",
      linewidth = 0.7
    ) +
    ggplot2::scale_color_manual(values = ratio_colors, name = NULL) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 0.20, by = 0.025),
      labels = scales::label_percent(accuracy = 0.1)
    ) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
    ggplot2::labs(
      title = "Clean-null known-null FDP lower bound",
      subtitle = "BF-adjusted single draw; dotted gray is the nominal diagonal",
      x = "Nominal alpha",
      y = "V / R_merged"
    ) +
    clean_theme()
}

plot_clean_plugin_fdr <- function() {
  bf <- clean_alpha_curves[clean_alpha_curves$fit_stage == "BF-adjusted", ]
  plot_data <- rbind(
    data.frame(
      nominal_alpha = bf$nominal_alpha,
      value = bf$target_pi0_plugin_fdr,
      m_ratio = bf$m_ratio,
      estimand = "Target-pi0 plug-in"
    ),
    data.frame(
      nominal_alpha = bf$nominal_alpha,
      value = bf$merged_pi0_plugin_fdr,
      m_ratio = bf$m_ratio,
      estimand = "Merged-pi0 plug-in"
    )
  )
  plot_data$ratio_label <- make_ratio_label(plot_data$m_ratio)
  plot_data <- plot_data[is.finite(plot_data$value), ]
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = nominal_alpha, y = value, color = ratio_label)
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      color = "#888888",
      linetype = "dotted",
      linewidth = 0.7
    ) +
    ggplot2::facet_wrap(~ estimand, nrow = 1) +
    ggplot2::scale_color_manual(values = ratio_colors, name = NULL) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 0.20, by = 0.025),
      labels = scales::label_percent(accuracy = 0.1)
    ) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggplot2::labs(
      title = "Clean-null plug-in FDR diagnostics",
      subtitle = "BF-adjusted single draw; values are untruncated",
      x = "Nominal alpha",
      y = "Plug-in estimated FDR"
    ) +
    clean_theme()
}

plot_clean_pi0 <- function() {
  plot_data <- clean_pi0
  plot_data$ratio_label <- factor(
    sprintf("%.2f", plot_data$m_ratio),
    levels = sprintf("%.2f", c(0, 0.05, 0.10, 0.20, 1.00))
  )
  plot_data <- rbind(
    data.frame(
      ratio_label = plot_data$ratio_label,
      fit_stage = plot_data$fit_stage,
      estimand = "Merged pi0",
      value = plot_data$pi0_merged
    ),
    data.frame(
      ratio_label = plot_data$ratio_label,
      fit_stage = plot_data$fit_stage,
      estimand = "Recovered target pi0",
      value = plot_data$pi0_target_unbounded
    )
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = ratio_label,
      y = value,
      color = estimand,
      group = estimand
    )
  ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::facet_wrap(~ fit_stage, ncol = 1, scales = "free_y") +
    ggplot2::scale_color_manual(
      values = c("Merged pi0" = "#0072B2", "Recovered target pi0" = "#D55E00"),
      name = NULL
    ) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggplot2::labs(
      title = "Null weights rise as clean null controls are appended",
      subtitle = "All ratios use nested subsets of the same single clean-null draw",
      x = "M / J",
      y = "Estimated null proportion"
    ) +
    clean_theme()
}
