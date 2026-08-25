#!/usr/bin/env Rscript

# Plot target- and merged-pi0 plug-in FDR curves for random target variants.

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
    "One independently sampled genotype-label map per time point;\n",
    "the response retains the original genotype signal;"
  )
} else if (identical(
  configuration$permutation_method,
  "signal_stripped_independent_time_residual"
)) {
  "One independently sampled donor map per time point;"
} else if (identical(
  configuration$permutation_method,
  "signal_stripped_unit_specific_residual_block"
)) {
  "One independently sampled trajectory-preserving donor map per unit;"
} else {
  "One shared residual-block donor map;"
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
n_merged <- nrow(stage_data)
pi0_target <- if (isTRUE(calibration_row$pi0_target_valid)) {
  calibration_row$pi0_target_unbounded
} else {
  NA_real_
}
pi0_merged <- calibration_row$pi0_merged
if (n_target != expected_group_size || n_null != expected_group_size ||
    (is.finite(pi0_target) && (pi0_target < 0 || pi0_target > 1)) ||
    !is.finite(pi0_merged) || pi0_merged < 0 || pi0_merged > 1) {
  stop("The target or merged null-fraction estimate is invalid.")
}

alpha_grid <- seq(0.001, 0.20, by = 0.001)
curve_rows <- lapply(alpha_grid, function(alpha) {
  selected_indices <- cumulative_lfdr_calls(stage_data$lfdr, alpha)
  selected <- seq_len(nrow(stage_data)) %in% selected_indices
  target_calls <- sum(selected & target)
  null_calls <- sum(selected & permuted_null)
  merged_calls <- target_calls + null_calls
  null_call_rate <- null_calls / n_null
  target_call_rate <- target_calls / n_target
  merged_call_rate <- merged_calls / n_merged
  data.frame(
    nominal_alpha = alpha,
    n_target = n_target,
    n_permuted_null = n_null,
    n_merged = n_merged,
    target_calls = target_calls,
    permuted_null_calls = null_calls,
    merged_calls = merged_calls,
    target_call_rate = target_call_rate,
    permuted_null_call_rate = null_call_rate,
    merged_call_rate = merged_call_rate,
    pi0_target = pi0_target,
    pi0_merged = pi0_merged,
    target_pi0_plugin_fdr = if (target_call_rate > 0 &&
        is.finite(pi0_target)) {
      pi0_target * null_call_rate / target_call_rate
    } else {
      NA_real_
    },
    merged_pi0_plugin_fdr = if (merged_call_rate > 0) {
      pi0_merged * null_call_rate / merged_call_rate
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
target_plugin_matches <- if (isTRUE(calibration_row$pi0_target_valid)) {
  isTRUE(all.equal(
    alpha_005$target_pi0_plugin_fdr,
    calibration_row$post_selection_fdr_target_from_pi0,
    tolerance = 1e-14
  ))
} else {
  is.na(alpha_005$target_pi0_plugin_fdr) &&
    is.na(calibration_row$post_selection_fdr_target_from_pi0)
}
if (nrow(alpha_005) != 1L ||
    alpha_005$target_calls != calibration_row$target_calls ||
    alpha_005$permuted_null_calls != calibration_row$permuted_null_calls ||
    alpha_005$merged_calls != calibration_row$total_calls ||
    !target_plugin_matches || abs(
      alpha_005$merged_pi0_plugin_fdr -
        calibration_row$scaled_fdr_merged_from_estimated_pi0
    ) > 1e-14) {
  stop("The alpha-0.05 plug-in curves do not reproduce saved diagnostics.")
}

curve_path <- file.path(
  output_directory,
  "random_variant_plugin_fdr_vs_nominal_alpha.csv"
)
utils::write.csv(curve_data, curve_path, row.names = FALSE)

plot_data <- rbind(
  data.frame(
    nominal_alpha = curve_data$nominal_alpha,
    estimand = "Target-set pi0",
    estimated_fdr = curve_data$target_pi0_plugin_fdr,
    stringsAsFactors = FALSE
  ),
  data.frame(
    nominal_alpha = curve_data$nominal_alpha,
    estimand = "Merged-set pi0",
    estimated_fdr = curve_data$merged_pi0_plugin_fdr,
    stringsAsFactors = FALSE
  )
)
plot_data <- plot_data[is.finite(plot_data$estimated_fdr), , drop = FALSE]
plot_data$estimand <- factor(
  plot_data$estimand,
  levels = c("Target-set pi0", "Merged-set pi0")
)
point_data <- data.frame(
  nominal_alpha = 0.05,
  estimand = factor(
    c("Target-set pi0", "Merged-set pi0"),
    levels = levels(plot_data$estimand)
  ),
  estimated_fdr = c(
    alpha_005$target_pi0_plugin_fdr,
    alpha_005$merged_pi0_plugin_fdr
  )
)
point_data <- point_data[is.finite(point_data$estimated_fdr), , drop = FALSE]
point_data$label <- ifelse(
  point_data$estimand == "Target-set pi0",
  sprintf("Target pi0: %.1f%%", 100 * point_data$estimated_fdr),
  sprintf("Merged pi0: %.1f%%", 100 * point_data$estimated_fdr)
)
method_colors <- c(
  "Target-set pi0" = "#0072B2",
  "Merged-set pi0" = "#CC79A7"
)
y_upper <- max(0.22, 1.08 * max(plot_data$estimated_fdr))

figure <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = nominal_alpha, y = estimated_fdr, color = estimand)
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    color = "#6B6B6B",
    linewidth = 0.8,
    linetype = "dashed"
  ) +
  ggplot2::geom_step(linewidth = 1.15, direction = "hv") +
  ggplot2::geom_vline(
    xintercept = 0.05,
    color = "#9A9A9A",
    linewidth = 0.55,
    linetype = "dotted"
  ) +
  ggplot2::geom_point(
    data = point_data,
    size = 2.8,
    fill = "white",
    shape = 21,
    stroke = 1
  ) +
  ggplot2::geom_label(
    data = point_data,
    ggplot2::aes(label = label),
    nudge_x = 0.004,
    hjust = 0,
    size = 3.4,
    fill = "white",
    linewidth = 0.2,
    show.legend = FALSE
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
  ggplot2::scale_color_manual(
    values = method_colors,
    labels = c(
      if (is.finite(pi0_target)) {
        sprintf("Target-set pi0 = %.3f", pi0_target)
      } else {
        "Target-set pi0 invalid"
      },
      sprintf("Merged-set pi0 = %.3f", pi0_merged)
    ),
    name = NULL
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
    title = "Plug-in estimated FDR versus nominal alpha",
    subtitle = paste(
      paste0("Random-variant ", method_label, " experiment;"),
      "BF-adjusted cumulative-lfdr calls"
    ),
    x = "Nominal alpha",
    y = "Plug-in estimated FDR",
    caption = if (identical(
      configuration$permutation_method,
      "genotype_label_independent_time"
    )) {
      paste(
        paste(
          "Target curve: pi0[target] x (V/M) / (R[target]/N[target]);",
          "merged curve: pi0[merged] x (V/M) / (R[merged]/N[merged])."
        ),
        "One independently sampled genotype-label map per time point.",
        paste(
          "The response retains the original genotype signal; curves are",
          "conditional diagnostics, not formal repeated-sampling FDR."
        ),
        sep = "\n"
      )
    } else {
      paste(
        paste(
          "Target curve: pi0[target] x (V/M) / (R[target]/N[target]);",
          "merged curve: pi0[merged] x (V/M) / (R[merged]/N[merged])."
        ),
        paste(
          randomization_caption, "curves are conditional",
          "diagnostics, not formal repeated-sampling FDR."
        ),
        sep = "\n"
      )
    }
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
    legend.position = "top",
    legend.justification = "left",
    legend.text = ggplot2::element_text(size = 10.5),
    plot.margin = ggplot2::margin(14, 18, 12, 14)
  )

figure_directory <- file.path(output_directory, "figures")
dir.create(figure_directory, showWarnings = FALSE)
png_path <- file.path(
  figure_directory,
  "random_variant_plugin_fdr_vs_nominal_alpha.png"
)
pdf_path <- file.path(
  figure_directory,
  "random_variant_plugin_fdr_vs_nominal_alpha.pdf"
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

cat("Random-variant plug-in FDR curve plot created.\n")
cat("Curve data: ", curve_path, "\n", sep = "")
cat("PNG: ", png_path, "\n", sep = "")
cat("PDF: ", pdf_path, "\n", sep = "")
print(alpha_005, row.names = FALSE)
