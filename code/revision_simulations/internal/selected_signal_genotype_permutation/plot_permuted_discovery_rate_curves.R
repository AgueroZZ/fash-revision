#!/usr/bin/env Rscript

# Plot the permuted-control discovery rate over nominal cumulative-FDR levels.

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

if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("scales", quietly = TRUE)) {
  stop("The ggplot2 and scales packages are required.")
}

workflowr_root <- find_workflowr_root()
output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "selected_signal_genotype_permutation_seed20260811"
)
lfdr_path <- file.path(output_directory, "unit_lfdr.csv")
figure_directory <- file.path(output_directory, "figures")
if (!file.exists(lfdr_path)) {
  stop("The fixed-seed unit-level lfdr table is missing.")
}
dir.create(figure_directory, showWarnings = FALSE)

unit_lfdr <- utils::read.csv(lfdr_path, stringsAsFactors = FALSE)
expected_stages <- c("Raw", "BF-adjusted")
if (!identical(unique(unit_lfdr$fit_stage), expected_stages) ||
    !setequal(unique(unit_lfdr$group), c("target", "permuted_null")) ||
    any(!is.finite(unit_lfdr$lfdr)) ||
    any(unit_lfdr$lfdr < 0 | unit_lfdr$lfdr > 1)) {
  stop("The unit-level lfdr table is invalid.")
}

make_curve <- function(stage_data, stage) {
  ordering <- order(stage_data$lfdr, method = "radix")
  ordered_lfdr <- stage_data$lfdr[ordering]
  ordered_null <- stage_data$group[ordering] == "permuted_null"
  n_permuted <- sum(stage_data$group == "permuted_null")
  data.frame(
    fit_stage = stage,
    discovery_rank = c(0L, seq_along(ordering)),
    alpha_boundary = c(0, cumsum(ordered_lfdr) / seq_along(ordering)),
    permuted_discoveries = c(0L, cumsum(ordered_null)),
    permuted_discovery_rate = c(0, cumsum(ordered_null) / n_permuted),
    maximum_selected_lfdr = c(NA_real_, ordered_lfdr),
    stringsAsFactors = FALSE
  )
}

curve_data <- do.call(rbind, lapply(expected_stages, function(stage) {
  make_curve(unit_lfdr[unit_lfdr$fit_stage == stage, ], stage)
}))
rownames(curve_data) <- NULL

alpha_grid <- c(0.001, 0.005, 0.01, 0.025, 0.05, 0.075, 0.10, 0.15, 0.20)
grid_summary <- do.call(rbind, lapply(expected_stages, function(stage) {
  curve <- curve_data[curve_data$fit_stage == stage, ]
  rows <- vapply(alpha_grid, function(alpha) {
    max(which(curve$alpha_boundary <= alpha))
  }, integer(1))
  data.frame(
    fit_stage = stage,
    alpha = alpha_grid,
    discovery_rank = curve$discovery_rank[rows],
    permuted_discoveries = curve$permuted_discoveries[rows],
    permuted_discovery_rate = curve$permuted_discovery_rate[rows],
    maximum_selected_lfdr = curve$maximum_selected_lfdr[rows],
    stringsAsFactors = FALSE
  )
}))
rownames(grid_summary) <- NULL

pi0_merged <- c(
  Raw = 0.0734625472866963,
  `BF-adjusted` = 0.465590484282073
)
curve_data$pi0_merged <- unname(pi0_merged[curve_data$fit_stage])
curve_data$merged_discovery_rate <- curve_data$discovery_rank / 2354
curve_data$estimated_fdr <- with(
  curve_data,
  pi0_merged * permuted_discovery_rate / merged_discovery_rate
)
curve_data$estimated_fdr[curve_data$discovery_rank == 0L] <- NA_real_
grid_summary$pi0_merged <- unname(pi0_merged[grid_summary$fit_stage])
grid_summary$merged_discovery_rate <- grid_summary$discovery_rank / 2354
grid_summary$estimated_fdr <- with(
  grid_summary,
  pi0_merged * permuted_discovery_rate / merged_discovery_rate
)

utils::write.csv(
  curve_data,
  file.path(output_directory, "permuted_discovery_rate_curve.csv"),
  row.names = FALSE
)
utils::write.csv(
  grid_summary,
  file.path(output_directory, "permuted_discovery_rate_alpha_grid.csv"),
  row.names = FALSE
)

plot_settings <- list(
  Raw = list(
    color = "#D55E00",
    title = "Raw fit: 88.5% of permuted controls enter by alpha = 0.05",
    filename = "permuted_discovery_rate_vs_alpha_raw"
  ),
  `BF-adjusted` = list(
    color = "#0072B2",
    title = paste(
      "BF-adjusted fit: 27.5% of permuted controls enter",
      "by alpha = 0.05"
    ),
    filename = "permuted_discovery_rate_vs_alpha_bf_adjusted"
  )
)

for (stage in expected_stages) {
  curve <- curve_data[
    curve_data$fit_stage == stage & curve_data$alpha_boundary <= 0.20,
  ]
  if (max(curve$alpha_boundary) < 0.20) {
    endpoint <- curve[nrow(curve), , drop = FALSE]
    endpoint$alpha_boundary <- 0.20
    curve <- rbind(curve, endpoint)
  }
  alpha_005 <- grid_summary[
    grid_summary$fit_stage == stage & grid_summary$alpha == 0.05,
  ]
  settings <- plot_settings[[stage]]
  label <- paste0(
    alpha_005$permuted_discoveries,
    " / 1,177 (",
    scales::percent(alpha_005$permuted_discovery_rate, accuracy = 0.1),
    ")"
  )

  figure <- ggplot2::ggplot(
    curve,
    ggplot2::aes(x = alpha_boundary, y = permuted_discovery_rate)
  ) +
    ggplot2::geom_step(
      direction = "hv",
      color = settings$color,
      linewidth = 1.15
    ) +
    ggplot2::geom_vline(
      xintercept = 0.05,
      color = "#4D4D4D",
      linewidth = 0.55,
      linetype = "22"
    ) +
    ggplot2::geom_point(
      data = alpha_005,
      ggplot2::aes(x = alpha, y = permuted_discovery_rate),
      inherit.aes = FALSE,
      color = settings$color,
      fill = "white",
      shape = 21,
      stroke = 1.1,
      size = 3.2
    ) +
    ggplot2::annotate(
      "text",
      x = 0.055,
      y = alpha_005$permuted_discovery_rate,
      label = label,
      hjust = 0,
      vjust = -0.65,
      size = 3.7,
      color = "#222222"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 0.20),
      breaks = seq(0, 0.20, by = 0.025),
      labels = scales::label_number(accuracy = 0.001),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.20),
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      title = settings$title,
      subtitle = paste(
        "Merged cumulative-lfdr rule applied to 1,177 selected targets and",
        "1,177 permuted controls; seed 20260811"
      ),
      x = expression(paste("Nominal cumulative-FDR level, ", alpha)),
      y = "Permuted controls selected",
      caption = paste(
        "One shared donor permutation; no uncertainty band.",
        "The y-axis is a conditional null inclusion rate, not empirical FDR."
      )
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 14,
        color = "#111111",
        margin = ggplot2::margin(b = 5)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 10.5,
        color = "#4D4D4D",
        margin = ggplot2::margin(b = 12)
      ),
      plot.caption = ggplot2::element_text(
        size = 9,
        color = "#5A5A5A",
        hjust = 0,
        margin = ggplot2::margin(t = 10)
      ),
      axis.title = ggplot2::element_text(size = 11.5, color = "#222222"),
      axis.text = ggplot2::element_text(size = 10.5, color = "#333333"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
      panel.grid.major.y = ggplot2::element_line(
        color = "#E6E6E6",
        linewidth = 0.4
      ),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(14, 18, 12, 14)
    )

  ggplot2::ggsave(
    file.path(figure_directory, paste0(settings$filename, ".png")),
    figure,
    width = 8.2,
    height = 5.8,
    dpi = 320,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_directory, paste0(settings$filename, ".pdf")),
    figure,
    width = 8.2,
    height = 5.8,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
}

fdr_plot_settings <- list(
  Raw = list(
    color = "#D55E00",
    title = "Raw fit: estimated FDR is 6.9% at alpha = 0.05",
    filename = "scaled_permutation_fdr_vs_alpha_raw"
  ),
  `BF-adjusted` = list(
    color = "#0072B2",
    title = "BF-adjusted fit: estimated FDR is 20.1% at alpha = 0.05",
    filename = "scaled_permutation_fdr_vs_alpha_bf_adjusted"
  )
)

for (stage in expected_stages) {
  curve <- curve_data[
    curve_data$fit_stage == stage &
      curve_data$alpha_boundary <= 0.20 &
      is.finite(curve_data$estimated_fdr),
  ]
  if (max(curve$alpha_boundary) < 0.20) {
    endpoint <- curve[nrow(curve), , drop = FALSE]
    endpoint$alpha_boundary <- 0.20
    curve <- rbind(curve, endpoint)
  }
  alpha_005 <- grid_summary[
    grid_summary$fit_stage == stage & grid_summary$alpha == 0.05,
  ]
  settings <- fdr_plot_settings[[stage]]
  label <- scales::percent(alpha_005$estimated_fdr, accuracy = 0.1)

  figure <- ggplot2::ggplot(
    curve,
    ggplot2::aes(x = alpha_boundary, y = estimated_fdr)
  ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      color = "#777777",
      linewidth = 0.65,
      linetype = "33"
    ) +
    ggplot2::geom_step(
      direction = "hv",
      color = settings$color,
      linewidth = 1.15
    ) +
    ggplot2::geom_vline(
      xintercept = 0.05,
      color = "#4D4D4D",
      linewidth = 0.55,
      linetype = "22"
    ) +
    ggplot2::geom_point(
      data = alpha_005,
      ggplot2::aes(x = alpha, y = estimated_fdr),
      inherit.aes = FALSE,
      color = settings$color,
      fill = "white",
      shape = 21,
      stroke = 1.1,
      size = 3.2
    ) +
    ggplot2::annotate(
      "text",
      x = 0.055,
      y = alpha_005$estimated_fdr,
      label = paste0("Estimated FDR = ", label),
      hjust = 0,
      vjust = -0.70,
      size = 3.7,
      color = "#222222"
    ) +
    ggplot2::annotate(
      "text",
      x = 0.173,
      y = 0.173,
      label = "Nominal: y = x",
      angle = 34,
      hjust = 0.5,
      vjust = -0.55,
      size = 3.25,
      color = "#666666"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, 0.20),
      breaks = seq(0, 0.20, by = 0.025),
      labels = scales::label_number(accuracy = 0.001),
      expand = ggplot2::expansion(mult = c(0, 0.01))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0, 0.45),
      breaks = seq(0, 0.45, by = 0.05),
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      title = settings$title,
      subtitle = paste0(
        "Estimated FDR = pi0 x P(permuted selected) / P(merged selected); ",
        "merged pi0 = ",
        format(pi0_merged[[stage]], digits = 4)
      ),
      x = expression(paste("Nominal cumulative-FDR level, ", alpha)),
      y = "Scaled permutation-FDR estimate",
      caption = paste(
        "One shared donor permutation; no uncertainty band.",
        "Both fitted merged pi0 estimates are below the known design bound of 0.5."
      )
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 14,
        color = "#111111",
        margin = ggplot2::margin(b = 5)
      ),
      plot.subtitle = ggplot2::element_text(
        size = 10.5,
        color = "#4D4D4D",
        margin = ggplot2::margin(b = 12)
      ),
      plot.caption = ggplot2::element_text(
        size = 9,
        color = "#5A5A5A",
        hjust = 0,
        margin = ggplot2::margin(t = 10)
      ),
      axis.title = ggplot2::element_text(size = 11.5, color = "#222222"),
      axis.text = ggplot2::element_text(size = 10.5, color = "#333333"),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
      panel.grid.major.y = ggplot2::element_line(
        color = "#E6E6E6",
        linewidth = 0.4
      ),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(14, 18, 12, 14)
    )

  ggplot2::ggsave(
    file.path(figure_directory, paste0(settings$filename, ".png")),
    figure,
    width = 8.2,
    height = 5.8,
    dpi = 320,
    bg = "white"
  )
  ggplot2::ggsave(
    file.path(figure_directory, paste0(settings$filename, ".pdf")),
    figure,
    width = 8.2,
    height = 5.8,
    device = grDevices::cairo_pdf,
    bg = "white"
  )
}

cat("Permuted-control discovery-rate and scaled-FDR figures created.\n")
