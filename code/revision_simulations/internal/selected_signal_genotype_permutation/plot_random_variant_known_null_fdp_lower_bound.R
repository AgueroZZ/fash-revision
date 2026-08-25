#!/usr/bin/env Rscript

# Plot the known-null FDP lower bound for the random-variant merged procedure.

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

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(args, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(args[equals_hit[1L]], nchar(equals_prefix) + 1L))
  }
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1L] == length(args)) {
    return(default)
  }
  args[hit[1L] + 1L]
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required.")
}
if (!requireNamespace("scales", quietly = TRUE)) {
  stop("The scales package is required.")
}

cumulative_lfdr_calls <- function(lfdr, alpha) {
  ordering <- order(lfdr, method = "radix")
  cumulative_mean <- cumsum(lfdr[ordering]) / seq_along(ordering)
  accepted <- which(cumulative_mean <= alpha)
  if (length(accepted) == 0L) {
    return(integer())
  }
  ordering[accepted]
}

workflowr_root <- find_workflowr_root()
default_output_id <- paste0(
  "selected_signal_random_all_tested_residual_block_permutation_",
  "selection20260817_seed20260811"
)
output_id <- get_arg("--output-id", default_output_id)
if (length(output_id) != 1L || !nzchar(output_id) ||
    grepl("/", output_id, fixed = TRUE)) {
  stop("The requested output ID is invalid.")
}
output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  output_id
)
unit_lfdr_path <- file.path(output_directory, "unit_lfdr.csv")
fit_path <- file.path(output_directory, "merged_fash_fit.rds")
calibration_path <- file.path(
  output_directory,
  "calibration_diagnostics.csv"
)
if (!file.exists(unit_lfdr_path) || !file.exists(calibration_path) ||
    !file.exists(fit_path)) {
  stop("The completed random-variant pilot cache is incomplete.")
}

unit_lfdr <- utils::read.csv(unit_lfdr_path, stringsAsFactors = FALSE)
calibration <- utils::read.csv(
  calibration_path,
  stringsAsFactors = FALSE
)
configuration <- readRDS(fit_path)$configuration
expected_group_size <- as.integer(configuration$n_target_units)
method_label <- if (identical(
  configuration$permutation_method,
  "signal_stripped_residual_block"
)) {
  "signal-stripped residual-block"
} else if (identical(
  configuration$permutation_method,
  "genotype_label_independent_time"
)) {
  "independent-time genotype-label"
} else if (identical(
  configuration$permutation_method,
  "signal_stripped_independent_time_residual"
)) {
  "signal-stripped independent-time residual"
} else if (identical(
  configuration$permutation_method,
  "signal_stripped_unit_specific_residual_block"
)) {
  "signal-stripped unit-specific residual-block"
} else if (identical(configuration$permutation_method, "residual_block")) {
  "covariate-only residual-block"
} else {
  configuration$permutation_method
}
randomization_caption <- if (identical(
  configuration$permutation_method,
  "genotype_label_independent_time"
)) {
  paste(
    "one independently sampled genotype-label map per time point;\n",
    "the response retains the original genotype signal."
  )
} else if (identical(
  configuration$permutation_method,
  "signal_stripped_independent_time_residual"
)) {
  "one independently sampled donor map per time point."
} else if (identical(
  configuration$permutation_method,
  "signal_stripped_unit_specific_residual_block"
)) {
  "one independently sampled trajectory-preserving donor map per unit."
} else {
  "one shared residual-block donor map."
}
stage_data <- unit_lfdr[
  unit_lfdr$fit_stage == "BF-adjusted",
  ,
  drop = FALSE
]
calibration_row <- calibration[
  calibration$fit_stage == "BF-adjusted",
  ,
  drop = FALSE
]
if (length(expected_group_size) != 1L || is.na(expected_group_size) ||
    expected_group_size < 1L ||
    !configuration$target_selection_method %in%
      c("random_all_tested", "random_all_genes") ||
    nrow(stage_data) != 2L * expected_group_size ||
    nrow(calibration_row) != 1L ||
    !setequal(unique(stage_data$group), c("target", "permuted_null")) ||
    any(!is.finite(stage_data$lfdr)) ||
    any(stage_data$lfdr < 0 | stage_data$lfdr > 1)) {
  stop("The saved random-variant BF-adjusted data are invalid.")
}

target <- stage_data$group == "target"
permuted_null <- stage_data$group == "permuted_null"
n_target <- sum(target)
n_null <- sum(permuted_null)
if (n_target != expected_group_size || n_null != expected_group_size) {
  stop("The merged target and known-null group sizes are invalid.")
}

alpha_grid <- seq(0.001, 0.20, by = 0.001)
curve_rows <- lapply(alpha_grid, function(alpha) {
  selected_indices <- cumulative_lfdr_calls(stage_data$lfdr, alpha)
  selected <- seq_len(nrow(stage_data)) %in% selected_indices
  target_calls <- sum(selected & target)
  null_calls <- sum(selected & permuted_null)
  merged_calls <- target_calls + null_calls
  data.frame(
    nominal_alpha = alpha,
    target_calls = target_calls,
    permuted_null_calls = null_calls,
    merged_calls = merged_calls,
    known_null_fdp_lower_bound = if (merged_calls > 0L) {
      null_calls / merged_calls
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
})
curve_data <- do.call(rbind, curve_rows)
rownames(curve_data) <- NULL

alpha_005 <- curve_data[
  abs(curve_data$nominal_alpha - 0.05) < 1e-12,
  ,
  drop = FALSE
]
if (nrow(alpha_005) != 1L ||
    alpha_005$target_calls != calibration_row$target_calls ||
    alpha_005$permuted_null_calls != calibration_row$permuted_null_calls ||
    alpha_005$merged_calls != calibration_row$total_calls ||
    abs(
      alpha_005$known_null_fdp_lower_bound -
        calibration_row$known_null_discovery_fraction
    ) > 1e-14) {
  stop("The known-null FDP lower-bound curve failed validation.")
}

curve_path <- file.path(
  output_directory,
  "random_variant_known_null_fdp_lower_bound_vs_nominal_alpha.csv"
)
utils::write.csv(curve_data, curve_path, row.names = FALSE)

curve_color <- "#0072B2"
y_upper <- max(
  0.22,
  1.08 * max(curve_data$known_null_fdp_lower_bound, na.rm = TRUE)
)
curve_data$ribbon_lower <- pmin(
  curve_data$nominal_alpha,
  curve_data$known_null_fdp_lower_bound
)
curve_data$ribbon_upper <- pmax(
  curve_data$nominal_alpha,
  curve_data$known_null_fdp_lower_bound
)
figure <- ggplot2::ggplot(
  curve_data,
  ggplot2::aes(x = nominal_alpha, y = known_null_fdp_lower_bound)
) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = ribbon_lower,
      ymax = ribbon_upper
    ),
    fill = curve_color,
    alpha = 0.10
  ) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    color = "#6B6B6B",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  ggplot2::geom_step(
    color = curve_color,
    linewidth = 1.2,
    direction = "hv"
  ) +
  ggplot2::geom_vline(
    xintercept = 0.05,
    color = "#9A9A9A",
    linewidth = 0.55,
    linetype = "dotted"
  ) +
  ggplot2::geom_point(
    data = alpha_005,
    color = curve_color,
    size = 3,
    fill = "white",
    shape = 21,
    stroke = 1
  ) +
  ggplot2::annotate(
    "label",
    x = 0.054,
    y = alpha_005$known_null_fdp_lower_bound,
    label = sprintf(
      "At 5%%: %.1f%% (%d/%d)",
      100 * alpha_005$known_null_fdp_lower_bound,
      alpha_005$permuted_null_calls,
      alpha_005$merged_calls
    ),
    hjust = 0,
    vjust = -0.25,
    size = 3.5,
    color = curve_color,
    fill = "white",
    linewidth = 0.2
  ) +
  ggplot2::annotate(
    "text",
    x = 0.175,
    y = 0.175,
    label = "Nominal reference",
    angle = 31,
    hjust = 0.5,
    vjust = -0.6,
    size = 3.2,
    color = "#5F5F5F"
  ) +
  ggplot2::scale_x_continuous(
    limits = c(0, 0.20),
    breaks = seq(0, 0.20, by = 0.025),
    labels = scales::label_percent(accuracy = 0.1),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, y_upper),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    title = "Known-null FDP lower bound versus nominal alpha",
    subtitle = paste(
      paste0("Random-variant ", method_label, " experiment;"),
      "BF-adjusted merged discoveries"
    ),
    x = "Nominal alpha",
    y = "Known-null FDP lower bound: V / R[merged]",
    caption = paste(
      paste(
        "V counts selected permuted-null units; R[merged] counts all selected",
        "target and permuted-null units."
      ),
      paste(
        "This is a lower bound on realized merged FDP because false target",
        "discoveries are unobserved;", randomization_caption
      ),
      sep = "\n"
    )
  ) +
  ggplot2::theme_minimal(base_size = 12) +
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
      margin = ggplot2::margin(b = 11)
    ),
    plot.caption = ggplot2::element_text(
      size = 8.8,
      color = "#5A5A5A",
      hjust = 0,
      lineheight = 1.05,
      margin = ggplot2::margin(t = 10)
    ),
    axis.title = ggplot2::element_text(size = 11.5, color = "#222222"),
    axis.text = ggplot2::element_text(size = 10, color = "#333333"),
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(
      color = "#E4E4E4",
      linewidth = 0.4
    ),
    plot.margin = ggplot2::margin(14, 18, 12, 14)
  )

figure_directory <- file.path(output_directory, "figures")
dir.create(figure_directory, showWarnings = FALSE)
png_path <- file.path(
  figure_directory,
  "random_variant_known_null_fdp_lower_bound_vs_nominal_alpha.png"
)
pdf_path <- file.path(
  figure_directory,
  "random_variant_known_null_fdp_lower_bound_vs_nominal_alpha.pdf"
)
ggplot2::ggsave(
  png_path,
  figure,
  width = 9.5,
  height = 6.4,
  dpi = 320,
  bg = "white"
)
ggplot2::ggsave(
  pdf_path,
  figure,
  width = 9.5,
  height = 6.4,
  device = grDevices::cairo_pdf,
  bg = "white"
)

cat("Random-variant known-null FDP lower-bound plot created.\n")
cat("Curve data: ", curve_path, "\n", sep = "")
cat("PNG: ", png_path, "\n", sep = "")
cat("PDF: ", pdf_path, "\n", sep = "")
print(alpha_005, row.names = FALSE)
