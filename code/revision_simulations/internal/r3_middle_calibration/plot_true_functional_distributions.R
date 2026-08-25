#!/usr/bin/env Rscript

# Reconstruct and plot true functional distributions for the frozen R3 design.

options(stringsAsFactors = FALSE)

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[[1L]] == length(args)) return(default)
  args[[hit[[1L]] + 1L]]
}

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }
  if (file.exists(file.path(
    "coderepo-local", "code", "revision_simulations", "shared",
    "simulation_functions.R"
  ))) {
    return(normalizePath(
      "coderepo-local", winslash = "/", mustWork = TRUE
    ))
  }
  stop("Could not find the workflowr repository root.", call. = FALSE)
}

summarize_values <- function(data) {
  groups <- split(
    data,
    list(data$truth_mechanism, data$target, data$truth_status),
    drop = TRUE
  )
  rows <- lapply(groups, function(x) {
    values <- x$functional_value
    quantiles <- stats::quantile(
      values,
      probs = c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99),
      names = FALSE
    )
    data.frame(
      truth_mechanism = x$truth_mechanism[[1L]],
      target = x$target[[1L]],
      truth_status = x$truth_status[[1L]],
      n = length(values),
      minimum = min(values),
      q01 = quantiles[[1L]],
      q05 = quantiles[[2L]],
      q25 = quantiles[[3L]],
      median = quantiles[[4L]],
      q75 = quantiles[[5L]],
      q95 = quantiles[[6L]],
      q99 = quantiles[[7L]],
      maximum = max(values),
      minimum_absolute_value = min(abs(values)),
      proportion_abs_below_0_15 = mean(abs(values) < 0.15),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$truth_mechanism, out$target, out$truth_status), ]
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The plot requires the existing ggplot2 installation.", call. = FALSE)
}

workflowr_root <- find_workflowr_root()
result_argument <- get_arg("--result-dir", "")
if (!nzchar(result_argument)) {
  stop("--result-dir is required.", call. = FALSE)
}
result_dir <- normalizePath(
  result_argument, winslash = "/", mustWork = TRUE
)
configuration_path <- file.path(result_dir, "configuration.rds")
complete_flag_path <- file.path(result_dir, "complete.flag")
if (!file.exists(configuration_path) || !file.exists(complete_flag_path)) {
  stop("The supplied R3 result is not complete.", call. = FALSE)
}

default_snapshot <- file.path(
  workflowr_root, "code", "revision_simulations", "r3_r4_fashr0143",
  "source_snapshots", "r3_prior_geometry_simulation_functions.R"
)
snapshot_path <- normalizePath(
  get_arg("--source-snapshot", default_snapshot),
  winslash = "/",
  mustWork = TRUE
)
source(snapshot_path)

configuration <- readRDS(configuration_path)
required_configuration <- c(
  "J", "time_grid", "evaluation_grid", "class_probs",
  "expected_class_counts", "expected_truth_group_counts",
  "truth_mechanisms", "dynamic_main_effect_sd", "random_bspline",
  "raised_cosine", "middle_window", "middle_boundary",
  "switch_threshold", "location_truth_margin", "switch_truth_margin",
  "non_switch_min_abs", "non_switch_min_range_fraction", "seed_list",
  "temporal_category_probs"
)
missing_configuration <- setdiff(
  required_configuration, names(configuration)
)
if (length(missing_configuration) > 0L) {
  stop(
    "The R3 configuration is missing: ",
    paste(missing_configuration, collapse = ", "),
    call. = FALSE
  )
}

expected_centers <- list(
  early = c(1.5, 1.5),
  middle = c(4.5, 10.5),
  late = c(13.5, 13.5)
)
expected_temporal_probs <- c(
  early = 0.29, middle = 0.42, late = 0.29
)
if (!identical(configuration$raised_cosine$center_ranges, expected_centers) ||
    !isTRUE(all.equal(
      configuration$temporal_category_probs,
      expected_temporal_probs,
      tolerance = 1e-12
    )) ||
    !identical(configuration$middle_window, c(3, 12)) ||
    !identical(configuration$middle_boundary, "open")) {
  stop(
    "The supplied result is not the frozen full-support mixture design.",
    call. = FALSE
  )
}

targets <- c("early", "middle", "late", "switch")
truth_rows <- list()
row_index <- 1L
for (truth_mechanism in configuration$truth_mechanisms) {
  for (seed in configuration$seed_list) {
    component_seeds <- revision_component_seeds(seed)
    effect_sim <- simulate_matched_functional_effect_set(
      n_variants = configuration$J,
      truth_mechanism = truth_mechanism,
      time_grid = configuration$time_grid,
      evaluation_grid = configuration$evaluation_grid,
      class_probs = configuration$class_probs,
      dynamic_main_effect_sd = configuration$dynamic_main_effect_sd,
      bspline_amplitude = configuration$random_bspline$amplitude,
      bspline_df = configuration$random_bspline$df,
      bspline_coefficient_sd =
        configuration$random_bspline$coefficient_sd,
      cosine_width_half = configuration$raised_cosine$width_half,
      cosine_spike_counts = configuration$raised_cosine$spike_counts,
      cosine_relative_amplitude_range =
        configuration$raised_cosine$relative_amplitude_range,
      cosine_target_centered_rms =
        configuration$raised_cosine$target_centered_rms,
      switch_threshold = configuration$switch_threshold,
      location_truth_margin = configuration$location_truth_margin,
      switch_truth_margin = configuration$switch_truth_margin,
      non_switch_min_abs = configuration$non_switch_min_abs,
      non_switch_min_range_fraction =
        configuration$non_switch_min_range_fraction,
      temporal_category_probs = configuration$temporal_category_probs,
      seed = seed,
      class_seed = component_seeds[["classes"]],
      constant_seed = component_seeds[["constant_effects"]],
      shape_seed = component_seeds[["functional_truth"]],
      scenario = paste0("r3_truth_distribution_", truth_mechanism),
      middle_window = configuration$middle_window,
      middle_boundary = configuration$middle_boundary
    )

    dynamic <- effect_sim$unit_info$effect_class == "dynamic_bspline"
    observed_counts <- table(factor(
      effect_sim$unit_info$truth_group[dynamic],
      levels = names(configuration$expected_truth_group_counts)
    ))
    if (!identical(
      stats::setNames(as.integer(observed_counts), names(observed_counts)),
      configuration$expected_truth_group_counts
    )) {
      stop("Reconstructed truth-group counts do not match the cache.")
    }

    for (target in targets) {
      values <- effect_sim$true_functionals[dynamic, target]
      truth_rows[[row_index]] <- data.frame(
        seed = seed,
        truth_mechanism = truth_mechanism,
        curve_index = which(dynamic),
        truth_group = effect_sim$unit_info$truth_group[dynamic],
        assigned_time_group =
          effect_sim$unit_info$time_group[dynamic],
        assigned_switch_status =
          effect_sim$unit_info$switch_status[dynamic],
        target = target,
        functional_value = values,
        truth_status = ifelse(
          values > 0,
          "Functional positive",
          "Functional null"
        ),
        stringsAsFactors = FALSE
      )
      row_index <- row_index + 1L
    }
  }
}

plot_data <- do.call(rbind, truth_rows)
rownames(plot_data) <- NULL
expected_rows <-
  length(configuration$truth_mechanisms) *
  length(configuration$seed_list) *
  unname(configuration$expected_class_counts[["dynamic_bspline"]]) *
  length(targets)
if (nrow(plot_data) != expected_rows || any(!is.finite(
  plot_data$functional_value
))) {
  stop("The reconstructed functional distribution is incomplete.")
}

plot_data$mechanism_label <- factor(
  plot_data$truth_mechanism,
  levels = c("random_bspline", "raised_cosine"),
  labels = c(
    "R3A: broad random B-spline",
    "R3B: compact raised cosine"
  )
)
plot_data$target_label <- factor(
  plot_data$target,
  levels = targets,
  labels = c("Early", "Middle", "Late", "Switch")
)
panel_levels <- as.vector(t(outer(
  levels(plot_data$mechanism_label),
  levels(plot_data$target_label),
  paste,
  sep = " - "
)))
plot_data$panel_label <- factor(
  paste(plot_data$mechanism_label, plot_data$target_label, sep = " - "),
  levels = panel_levels
)
plot_data$truth_status <- factor(
  plot_data$truth_status,
  levels = c("Functional null", "Functional positive")
)

output_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "diagnostics",
  "r3_middle_calibration"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
data_path <- file.path(
  output_dir, "r3_full_support_true_functional_values.csv"
)
summary_path <- file.path(
  output_dir, "r3_full_support_true_functional_summary.csv"
)
provenance_path <- file.path(
  output_dir, "r3_full_support_true_functional_provenance.txt"
)
utils::write.csv(
  plot_data[, c(
    "seed", "truth_mechanism", "curve_index", "truth_group",
    "assigned_time_group", "assigned_switch_status", "target",
    "functional_value", "truth_status"
  )],
  data_path,
  row.names = FALSE
)
utils::write.csv(
  summarize_values(plot_data),
  summary_path,
  row.names = FALSE
)
writeLines(
  c(
    paste0("result_dir=", result_dir),
    paste0("configuration_md5=", unname(tools::md5sum(
      configuration_path
    ))),
    paste0("source_snapshot=", snapshot_path),
    paste0("source_snapshot_md5=", unname(tools::md5sum(
      snapshot_path
    ))),
    paste0("seeds=", paste(configuration$seed_list, collapse = ",")),
    paste0("n_dynamic_per_seed=", unname(
      configuration$expected_class_counts[["dynamic_bspline"]]
    )),
    "truth_scope=dynamic_bspline_only",
    "grid=seq(0,15,by=0.1)",
    "middle_definition=3<t<12"
  ),
  provenance_path
)

colors <- c(
  "Functional null" = "#999999",
  "Functional positive" = "#0072B2"
)
plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = functional_value, fill = truth_status)
) +
  ggplot2::annotate(
    "rect",
    xmin = -0.10,
    xmax = 0.10,
    ymin = -Inf,
    ymax = Inf,
    fill = "#F4E3A1",
    alpha = 0.38
  ) +
  ggplot2::geom_histogram(
    binwidth = 0.05,
    boundary = 0,
    color = "white",
    linewidth = 0.12
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    color = "#111111",
    linewidth = 0.5
  ) +
  ggplot2::geom_vline(
    xintercept = c(-0.10, 0.10),
    color = "#8C6D1F",
    linetype = "dotted",
    linewidth = 0.45
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(panel_label),
    ncol = 4,
    scales = "free"
  ) +
  ggplot2::scale_fill_manual(values = colors, name = NULL) +
  ggplot2::labs(
    title = "The frozen R3 truth design excludes ambiguous functional values",
    subtitle = paste(
      "Five seeds; 1,272 dynamic truth curves per seed and mechanism;",
      "bin width = 0.05"
    ),
    x = "True functional value F(beta)",
    y = "Number of dynamic truth curves (panel-specific scale)",
    caption = paste(
      "The shaded interval (-0.10, 0.10) is excluded by the truth-acceptance",
      "margins. The Switch null point mass at -0.25 represents same-sign curves."
    )
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(
      color = "#E2E2E2", linewidth = 0.35
    ),
    strip.text = ggplot2::element_text(face = "bold", size = 10.5),
    strip.background = ggplot2::element_rect(
      fill = "#F2F2F2", color = NA
    ),
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(color = "#444444"),
    plot.caption = ggplot2::element_text(
      color = "#444444", hjust = 0, size = 9
    ),
    legend.position = "top",
    axis.title = ggplot2::element_text(size = 10.5),
    panel.spacing = grid::unit(0.9, "lines")
  )

png_path <- file.path(
  output_dir, "r3_full_support_true_functional_histograms.png"
)
pdf_path <- file.path(
  output_dir, "r3_full_support_true_functional_histograms.pdf"
)
ggplot2::ggsave(
  png_path, plot = plot, width = 14.2, height = 8.2, dpi = 220,
  bg = "white"
)
ggplot2::ggsave(
  pdf_path, plot = plot, width = 14.2, height = 8.2,
  device = grDevices::cairo_pdf, bg = "white"
)

cat(png_path, "\n")
cat(pdf_path, "\n")
cat(data_path, "\n")
cat(summary_path, "\n")
