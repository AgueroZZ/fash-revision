# Reporting helpers for the small-M matched-null EB sensitivity page.

if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("scales", quietly = TRUE) ||
    !requireNamespace("knitr", quietly = TRUE)) {
  stop("ggplot2, scales, and knitr are required for small-M reporting.")
}

small_m_output_id <- paste0(
  "all_gene_random_variant_small_m_signal_stripped_residual_sensitivity_",
  "selection20260817_seed20260811_subsets20260819"
)
small_m_output_directory <- file.path(
  "output", "revision_simulations", "internal", small_m_output_id
)
small_m_required_files <- c(
  "configuration.rds",
  "input_provenance.csv",
  "subset_membership.csv",
  "replicate_pi0.csv",
  "replicate_alpha005.csv",
  "replicate_alpha_curves.csv",
  "replicate_prior_weights.csv",
  "reference_summaries.csv",
  "runtime.csv",
  "validation.csv"
)
small_m_required_paths <- file.path(
  small_m_output_directory,
  small_m_required_files
)
if (any(!file.exists(small_m_required_paths))) {
  stop("The completed small-M sensitivity artifacts are incomplete.")
}

small_m_configuration <- readRDS(file.path(
  small_m_output_directory,
  "configuration.rds"
))
small_m_input_provenance <- utils::read.csv(file.path(
  small_m_output_directory,
  "input_provenance.csv"
), stringsAsFactors = FALSE)
small_m_subset_membership <- utils::read.csv(file.path(
  small_m_output_directory,
  "subset_membership.csv"
), stringsAsFactors = FALSE)
small_m_pi0 <- utils::read.csv(file.path(
  small_m_output_directory,
  "replicate_pi0.csv"
), stringsAsFactors = FALSE)
small_m_alpha005 <- utils::read.csv(file.path(
  small_m_output_directory,
  "replicate_alpha005.csv"
), stringsAsFactors = FALSE)
small_m_alpha_curves <- utils::read.csv(file.path(
  small_m_output_directory,
  "replicate_alpha_curves.csv"
), stringsAsFactors = FALSE)
small_m_prior_weights <- utils::read.csv(file.path(
  small_m_output_directory,
  "replicate_prior_weights.csv"
), stringsAsFactors = FALSE)
small_m_runtime <- utils::read.csv(file.path(
  small_m_output_directory,
  "runtime.csv"
), stringsAsFactors = FALSE)
small_m_validation <- utils::read.csv(file.path(
  small_m_output_directory,
  "validation.csv"
), stringsAsFactors = FALSE)

small_m_ratios <- as.numeric(small_m_configuration$m_ratios)
small_m_ratio_levels <- c(0, small_m_ratios, 1)
small_m_ratio_labels <- c(
  "0" = "0 (target only)",
  "0.05" = "0.05",
  "0.1" = "0.10",
  "0.2" = "0.20",
  "1" = "1.00"
)
small_m_fit_stage_levels <- c("Raw", "BF-adjusted")

make_ratio_factor <- function(value) {
  factor(
    as.character(value),
    levels = as.character(small_m_ratio_levels),
    labels = unname(small_m_ratio_labels[
      as.character(small_m_ratio_levels)
    ])
  )
}

if (!identical(small_m_configuration$n_target, 6362L) ||
    !identical(small_m_configuration$n_replicates, 20L) ||
    !identical(small_m_configuration$permutation_method,
               "signal_stripped_residual_block") ||
    !identical(small_m_configuration$target_selection_method,
               "random_all_genes") ||
    !all(small_m_validation$passed) ||
    !setequal(unique(small_m_pi0$m_ratio), small_m_ratio_levels) ||
    !setequal(unique(small_m_pi0$fit_stage), small_m_fit_stage_levels) ||
    nrow(small_m_subset_membership) !=
      small_m_configuration$n_replicates *
        sum(small_m_configuration$m_sizes) ||
    any(!is.finite(small_m_pi0$pi0_merged)) ||
    any(small_m_pi0$pi0_merged < 0 | small_m_pi0$pi0_merged > 1) ||
    any(!is.finite(small_m_prior_weights$prior_weight)) ||
    any(small_m_prior_weights$prior_weight < 0)) {
  stop("The small-M reporting contract failed.")
}

small_curve <- small_m_alpha_curves[!small_m_alpha_curves$is_reference, ]
small_curve_counts <- aggregate(
  nominal_alpha ~ replicate_id + m_ratio + fit_stage,
  data = small_curve,
  FUN = length
)
if (nrow(small_curve_counts) !=
      small_m_configuration$n_replicates * length(small_m_ratios) * 2L ||
    any(small_curve_counts$nominal_alpha !=
      length(small_m_configuration$alpha_grid)) ||
    any(abs(small_m_alpha005$nominal_alpha - 0.05) > 1e-12)) {
  stop("The small-M alpha-curve grid is incomplete.")
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

format_interval <- function(values,
                            formatter = function(x) format_decimal(x, 3L)) {
  quantiles <- stats::quantile(
    values,
    probs = c(0.25, 0.50, 0.75),
    na.rm = TRUE,
    names = FALSE,
    type = 8
  )
  paste0(
    formatter(quantiles[2]),
    " [",
    formatter(quantiles[1]),
    ", ",
    formatter(quantiles[3]),
    "]"
  )
}

render_small_m_table <- function(value, align = NULL) {
  print(knitr::kable(
    value,
    format = "html",
    align = align,
    escape = TRUE,
    row.names = FALSE,
    table.attr = 'class="table table-striped table-condensed"'
  ))
}

summarize_metric_curve <- function(data, metric_name, estimand_label) {
  values <- data[[metric_name]]
  valid <- is.finite(values)
  data <- data[valid, c("m_ratio", "nominal_alpha"), drop = FALSE]
  data$value <- values[valid]
  split_key <- interaction(
    data$m_ratio,
    data$nominal_alpha,
    drop = TRUE,
    lex.order = TRUE
  )
  rows <- lapply(split(data, split_key), function(group_data) {
    quantiles <- stats::quantile(
      group_data$value,
      probs = c(0.25, 0.50, 0.75),
      names = FALSE,
      type = 8
    )
    data.frame(
      m_ratio = group_data$m_ratio[1],
      nominal_alpha = group_data$nominal_alpha[1],
      q25 = quantiles[1],
      median = quantiles[2],
      q75 = quantiles[3],
      estimand = estimand_label,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

prepare_pi0_plot_data <- function() {
  long <- rbind(
    data.frame(
      small_m_pi0[c(
        "replicate_id", "m_ratio", "m_size", "fit_stage", "is_reference"
      )],
      estimand = "Merged pi0",
      value = small_m_pi0$pi0_merged,
      stringsAsFactors = FALSE
    ),
    data.frame(
      small_m_pi0[c(
        "replicate_id", "m_ratio", "m_size", "fit_stage", "is_reference"
      )],
      estimand = "Recovered target pi0",
      value = small_m_pi0$pi0_target_unbounded,
      stringsAsFactors = FALSE
    )
  )
  split_key <- interaction(
    long$m_ratio,
    long$fit_stage,
    long$estimand,
    drop = TRUE,
    lex.order = TRUE
  )
  rows <- lapply(split(long, split_key), function(group_data) {
    quantiles <- stats::quantile(
      group_data$value,
      probs = c(0.25, 0.50, 0.75),
      names = FALSE,
      type = 8
    )
    data.frame(
      m_ratio = group_data$m_ratio[1],
      fit_stage = group_data$fit_stage[1],
      estimand = group_data$estimand[1],
      q25 = quantiles[1],
      median = quantiles[2],
      q75 = quantiles[3],
      is_reference = all(group_data$is_reference),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  result$fit_stage <- factor(
    result$fit_stage,
    levels = small_m_fit_stage_levels
  )
  result$m_ratio_label <- make_ratio_factor(result$m_ratio)
  result
}

plot_small_m_pi0 <- function() {
  plot_data <- prepare_pi0_plot_data()
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = m_ratio_label,
      y = median,
      color = estimand,
      group = estimand
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = q25, ymax = q75, fill = estimand),
      alpha = 0.13,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.95) +
    ggplot2::geom_point(
      ggplot2::aes(shape = is_reference),
      size = 2.7,
      stroke = 0.9
    ) +
    ggplot2::facet_wrap(~ fit_stage, scales = "free_y", ncol = 1) +
    ggplot2::scale_color_manual(
      values = c("Merged pi0" = "#0072B2", "Recovered target pi0" = "#D55E00"),
      name = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = c("Merged pi0" = "#0072B2", "Recovered target pi0" = "#D55E00"),
      guide = "none"
    ) +
    ggplot2::scale_shape_manual(
      values = c("FALSE" = 16, "TRUE" = 21),
      labels = c("Small-M median", "Reference fit"),
      name = NULL
    ) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggplot2::labs(
      title = "Estimated null proportions change with the appended-null ratio",
      subtitle = "Lines join medians; ribbons are interquartile ranges across 20 nested null subsets",
      x = "M / J",
      y = "Estimated null proportion"
    ) +
    small_m_theme()
}

prepare_conditional_alternative_weights <- function() {
  prior <- small_m_prior_weights
  fit_keys <- unique(prior[c(
    "replicate_id", "m_ratio", "m_size", "fit_stage", "is_reference"
  )])
  psd_values <- sort(unique(prior$psd[prior$psd > 0]))
  completed <- lapply(seq_len(nrow(fit_keys)), function(index) {
    key <- fit_keys[index, ]
    current <- prior[
      prior$replicate_id == key$replicate_id &
        prior$m_ratio == key$m_ratio &
        prior$fit_stage == key$fit_stage &
        prior$is_reference == key$is_reference &
        prior$psd > 0,
      c("psd", "prior_weight"),
      drop = FALSE
    ]
    result <- merge(
      data.frame(psd = psd_values),
      current,
      by = "psd",
      all.x = TRUE,
      sort = TRUE
    )
    result$prior_weight[is.na(result$prior_weight)] <- 0
    alternative_mass <- sum(result$prior_weight)
    result$conditional_weight <- if (alternative_mass > 0) {
      result$prior_weight / alternative_mass
    } else {
      0
    }
    result$replicate_id <- key$replicate_id
    result$m_ratio <- key$m_ratio
    result$m_size <- key$m_size
    result$fit_stage <- key$fit_stage
    result$is_reference <- key$is_reference
    result
  })
  completed <- do.call(rbind, completed)
  split_key <- interaction(
    completed$m_ratio,
    completed$fit_stage,
    completed$psd,
    drop = TRUE,
    lex.order = TRUE
  )
  rows <- lapply(split(completed, split_key), function(group_data) {
    data.frame(
      m_ratio = group_data$m_ratio[1],
      fit_stage = group_data$fit_stage[1],
      psd = group_data$psd[1],
      median_conditional_weight = stats::median(
        group_data$conditional_weight
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

small_m_alternative_weights <- prepare_conditional_alternative_weights()

build_alternative_summary_table <- function() {
  prior <- small_m_prior_weights
  fit_key <- interaction(
    prior$replicate_id,
    prior$m_ratio,
    prior$fit_stage,
    prior$is_reference,
    drop = TRUE
  )
  per_fit <- do.call(rbind, lapply(split(prior, fit_key), function(current) {
    alternative <- current[current$psd > 0, , drop = FALSE]
    alternative_mass <- sum(alternative$prior_weight)
    data.frame(
      m_ratio = current$m_ratio[1],
      fit_stage = current$fit_stage[1],
      alternative_mass = alternative_mass,
      conditional_mean_psd = if (alternative_mass > 0) {
        sum(alternative$psd * alternative$prior_weight) / alternative_mass
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }))
  rows <- lapply(small_m_fit_stage_levels, function(fit_stage) {
    lapply(small_m_ratio_levels, function(m_ratio) {
      current <- per_fit[
        per_fit$fit_stage == fit_stage &
          abs(per_fit$m_ratio - m_ratio) < 1e-12,
        ,
        drop = FALSE
      ]
      data.frame(
        Stage = fit_stage,
        `M / J` = unname(small_m_ratio_labels[as.character(m_ratio)]),
        `Alternative mass` = format_interval(
          current$alternative_mass,
          function(x) format_percent(x, 2L)
        ),
        `Conditional mean PSD` = format_interval(
          current$conditional_mean_psd,
          function(x) format_decimal(x, 4L)
        ),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })
  })
  do.call(rbind, unlist(rows, recursive = FALSE))
}

small_m_alternative_summary_table <- build_alternative_summary_table()

plot_small_m_alternative_weights <- function() {
  plot_data <- small_m_alternative_weights
  plot_data$m_ratio_label <- make_ratio_factor(plot_data$m_ratio)
  plot_data$psd_label <- factor(
    sprintf("%.4g", plot_data$psd),
    levels = sprintf("%.4g", sort(unique(plot_data$psd), decreasing = TRUE))
  )
  plot_data$fit_stage <- factor(
    plot_data$fit_stage,
    levels = small_m_fit_stage_levels
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = m_ratio_label,
      y = psd_label,
      fill = median_conditional_weight
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.2) +
    ggplot2::facet_wrap(~ fit_stage, nrow = 1) +
    ggplot2::scale_fill_viridis_c(
      option = "C",
      trans = "sqrt",
      breaks = c(0, 0.10, 0.20, 0.40),
      labels = scales::label_percent(accuracy = 1),
      name = "Median conditional\nweight"
    ) +
    ggplot2::labs(
      title = "Alternative-mixture allocation changes after conditioning out pi0",
      subtitle = "Each replicate is normalized over non-null PSD scales before taking medians; omitted scales are filled with zero",
      x = "M / J",
      y = "Alternative PSD scale"
    ) +
    small_m_theme() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 7.8),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right",
      legend.justification = "center"
    ) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(
      barheight = grid::unit(55, "mm"),
      barwidth = grid::unit(5, "mm")
    ))
}

plot_small_m_plugin_fdr <- function() {
  bf_small <- small_m_alpha_curves[
    !small_m_alpha_curves$is_reference &
      small_m_alpha_curves$fit_stage == "BF-adjusted",
    ,
    drop = FALSE
  ]
  plot_data <- rbind(
    summarize_metric_curve(
      bf_small,
      "target_pi0_plugin_fdr",
      "Target-pi0 plug-in"
    ),
    summarize_metric_curve(
      bf_small,
      "merged_pi0_plugin_fdr",
      "Merged-pi0 plug-in"
    )
  )
  plot_data$m_ratio_label <- factor(
    sprintf("M/J = %.2f", plot_data$m_ratio),
    levels = sprintf("M/J = %.2f", small_m_ratios)
  )
  reference <- small_m_alpha_curves[
    small_m_alpha_curves$is_reference &
      small_m_alpha_curves$m_ratio == 1 &
      small_m_alpha_curves$fit_stage == "BF-adjusted",
    ,
    drop = FALSE
  ]
  reference <- rbind(
    data.frame(
      nominal_alpha = reference$nominal_alpha,
      value = reference$target_pi0_plugin_fdr,
      estimand = "Target-pi0 plug-in"
    ),
    data.frame(
      nominal_alpha = reference$nominal_alpha,
      value = reference$merged_pi0_plugin_fdr,
      estimand = "Merged-pi0 plug-in"
    )
  )
  reference <- reference[is.finite(reference$value), ]
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = nominal_alpha,
      y = median,
      color = m_ratio_label,
      fill = m_ratio_label
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = q25, ymax = q75),
      alpha = 0.10,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_line(
      data = reference,
      ggplot2::aes(x = nominal_alpha, y = value),
      inherit.aes = FALSE,
      color = "#333333",
      linetype = "dashed",
      linewidth = 0.75
    ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      color = "#777777",
      linetype = "dotted",
      linewidth = 0.65
    ) +
    ggplot2::facet_wrap(~ estimand, nrow = 1) +
    ggplot2::scale_color_brewer(palette = "Dark2", name = NULL) +
    ggplot2::scale_fill_brewer(palette = "Dark2", name = NULL) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 0.20, by = 0.025),
      labels = scales::label_percent(accuracy = 0.1)
    ) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 1)) +
    ggplot2::labs(
      title = "Plug-in estimated FDR remains far above nominal at small M",
      subtitle = "Solid lines are medians and ribbons are IQRs; dashed black is the saved M/J = 1 reference",
      x = "Nominal alpha",
      y = "Plug-in estimated FDR"
    ) +
    small_m_theme()
}

plot_small_m_known_null_fdp <- function() {
  bf_small <- small_m_alpha_curves[
    !small_m_alpha_curves$is_reference &
      small_m_alpha_curves$fit_stage == "BF-adjusted",
    ,
    drop = FALSE
  ]
  plot_data <- summarize_metric_curve(
    bf_small,
    "known_null_discovery_fraction",
    "Known-null FDP lower bound"
  )
  plot_data$m_ratio_label <- factor(
    sprintf("M/J = %.2f", plot_data$m_ratio),
    levels = sprintf("M/J = %.2f", small_m_ratios)
  )
  reference <- small_m_alpha_curves[
    small_m_alpha_curves$is_reference &
      small_m_alpha_curves$m_ratio == 1 &
      small_m_alpha_curves$fit_stage == "BF-adjusted",
    ,
    drop = FALSE
  ]
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = nominal_alpha,
      y = median,
      color = m_ratio_label,
      fill = m_ratio_label
    )
  ) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = q25, ymax = q75),
      alpha = 0.10,
      color = NA
    ) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_line(
      data = reference,
      ggplot2::aes(
        x = nominal_alpha,
        y = known_null_discovery_fraction
      ),
      inherit.aes = FALSE,
      color = "#333333",
      linetype = "dashed",
      linewidth = 0.75
    ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      color = "#777777",
      linetype = "dotted",
      linewidth = 0.65
    ) +
    ggplot2::scale_color_brewer(palette = "Dark2", name = NULL) +
    ggplot2::scale_fill_brewer(palette = "Dark2", name = NULL) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, 0.20, by = 0.025),
      labels = scales::label_percent(accuracy = 0.1)
    ) +
    ggplot2::scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
    ggplot2::labs(
      title = "Known-null discoveries form a larger share as M increases",
      subtitle = "Solid lines are medians and ribbons are IQRs; dashed black is the saved M/J = 1 reference",
      x = "Nominal alpha",
      y = "V / R_merged"
    ) +
    small_m_theme()
}

small_m_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 14,
        color = "#111111",
        margin = ggplot2::margin(b = 5)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 10.2,
        color = "#4D4D4D",
        margin = ggplot2::margin(b = 10)
      ),
      axis.title = ggplot2::element_text(size = 11, color = "#222222"),
      axis.text = ggplot2::element_text(size = 9.5, color = "#333333"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        color = "#E5E5E5",
        linewidth = 0.35
      ),
      legend.position = "top",
      legend.justification = "left",
      legend.text = ggplot2::element_text(size = 9.5),
      strip.text = ggplot2::element_text(face = "bold", size = 10.5),
      plot.margin = ggplot2::margin(12, 14, 10, 12)
    )
}

build_alpha005_table <- function() {
  bf_small <- small_m_alpha005[
    !small_m_alpha005$is_reference &
      small_m_alpha005$fit_stage == "BF-adjusted",
    ,
    drop = FALSE
  ]
  rows <- lapply(small_m_ratios, function(m_ratio) {
    current <- bf_small[abs(bf_small$m_ratio - m_ratio) < 1e-12, ]
    data.frame(
      `M / J` = sprintf("%.2f", m_ratio),
      M = format_integer(unique(current$m_size)),
      `Merged pi0` = format_interval(
        current$pi0_merged,
        function(x) format_decimal(x, 3L)
      ),
      `Recovered target pi0` = format_interval(
        current$pi0_target_unbounded,
        function(x) format_decimal(x, 3L)
      ),
      `Target calls` = format_interval(current$target_calls, format_integer),
      `Control calls V` = format_interval(
        current$permuted_null_calls,
        format_integer
      ),
      `V / R_merged` = format_interval(
        current$known_null_discovery_fraction,
        function(x) format_percent(x, 2L)
      ),
      `Target-pi0 plug-in FDR` = format_interval(
        current$target_pi0_plugin_fdr,
        function(x) format_percent(x, 1L)
      ),
      `Merged-pi0 plug-in FDR` = format_interval(
        current$merged_pi0_plugin_fdr,
        function(x) format_percent(x, 1L)
      ),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })
  full_reference <- small_m_alpha005[
    small_m_alpha005$is_reference &
      small_m_alpha005$m_ratio == 1 &
      small_m_alpha005$fit_stage == "BF-adjusted",
    ,
    drop = FALSE
  ]
  rows[[length(rows) + 1L]] <- data.frame(
    `M / J` = "1.00 reference",
    M = format_integer(full_reference$m_size),
    `Merged pi0` = format_decimal(full_reference$pi0_merged, 3L),
    `Recovered target pi0` = format_decimal(
      full_reference$pi0_target_unbounded,
      3L
    ),
    `Target calls` = format_integer(full_reference$target_calls),
    `Control calls V` = format_integer(full_reference$permuted_null_calls),
    `V / R_merged` = format_percent(
      full_reference$known_null_discovery_fraction,
      2L
    ),
    `Target-pi0 plug-in FDR` = format_percent(
      full_reference$target_pi0_plugin_fdr,
      1L
    ),
    `Merged-pi0 plug-in FDR` = format_percent(
      full_reference$merged_pi0_plugin_fdr,
      1L
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  do.call(rbind, rows)
}

small_m_alpha005_table <- build_alpha005_table()

build_pi0_summary_table <- function() {
  rows <- lapply(small_m_fit_stage_levels, function(fit_stage) {
    lapply(small_m_ratio_levels, function(m_ratio) {
      current <- small_m_pi0[
        small_m_pi0$fit_stage == fit_stage &
          abs(small_m_pi0$m_ratio - m_ratio) < 1e-12,
        ,
        drop = FALSE
      ]
      data.frame(
        Stage = fit_stage,
        `M / J` = unname(small_m_ratio_labels[as.character(m_ratio)]),
        M = format_integer(unique(current$m_size)),
        `Merged pi0` = format_interval(
          current$pi0_merged,
          function(x) format_decimal(x, 4L)
        ),
        `Recovered target pi0` = format_interval(
          current$pi0_target_unbounded,
          function(x) format_decimal(x, 4L)
        ),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    })
  })
  do.call(rbind, unlist(rows, recursive = FALSE))
}

small_m_pi0_summary_table <- build_pi0_summary_table()

small_m_total_runtime <- small_m_configuration$total_seconds
small_m_median_fit_runtime <- stats::median(
  small_m_runtime$total_seconds[!small_m_runtime$is_reference]
)
