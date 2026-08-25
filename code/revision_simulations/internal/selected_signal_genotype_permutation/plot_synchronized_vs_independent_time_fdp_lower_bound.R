#!/usr/bin/env Rscript

# Compare synchronized and independent-time residual null calibration.

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
output_parent <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal"
)
synchronized_id <- paste0(
  "all_gene_random_variant_signal_stripped_residual_block_permutation_",
  "selection20260817_seed20260811"
)
independent_id <- paste0(
  "all_gene_random_variant_signal_stripped_independent_time_residual_",
  "permutation_selection20260817_seed20260811"
)
experiment_ids <- c(
  "Synchronized donor trajectories" = synchronized_id,
  "Independent map at each time" = independent_id
)

curve_data <- do.call(rbind, lapply(names(experiment_ids), function(method) {
  path <- file.path(
    output_parent,
    experiment_ids[[method]],
    "random_variant_known_null_fdp_lower_bound_vs_nominal_alpha.csv"
  )
  if (!file.exists(path)) {
    stop("A required calibration curve is missing: ", path)
  }
  data <- utils::read.csv(path, stringsAsFactors = FALSE)
  data$method <- method
  data
}))
rownames(curve_data) <- NULL
curve_data$method <- factor(curve_data$method, levels = names(experiment_ids))
required_columns <- c(
  "nominal_alpha", "target_calls", "permuted_null_calls", "merged_calls",
  "known_null_fdp_lower_bound", "method"
)
if (!all(required_columns %in% names(curve_data)) ||
    any(!is.finite(curve_data$nominal_alpha)) ||
    any(!is.finite(curve_data$known_null_fdp_lower_bound)) ||
    !identical(
      curve_data$nominal_alpha[curve_data$method == levels(curve_data$method)[1L]],
      curve_data$nominal_alpha[curve_data$method == levels(curve_data$method)[2L]]
    )) {
  stop("The paired calibration curves failed structural validation.")
}

alpha_five <- curve_data[curve_data$nominal_alpha == 0.05, , drop = FALSE]
if (nrow(alpha_five) != 2L) {
  stop("The paired curves do not each contain nominal alpha 0.05.")
}
alpha_five$label <- paste0(
  scales::percent(alpha_five$known_null_fdp_lower_bound, accuracy = 0.1),
  " (",
  alpha_five$permuted_null_calls,
  "/",
  alpha_five$merged_calls,
  ")"
)

colors <- c(
  "Synchronized donor trajectories" = "#D55E00",
  "Independent map at each time" = "#0072B2"
)
comparison_plot <- ggplot2::ggplot(
  curve_data,
  ggplot2::aes(
    x = nominal_alpha,
    y = known_null_fdp_lower_bound,
    color = method
  )
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linewidth = 0.7,
    linetype = "dashed",
    color = "grey45"
  ) +
  ggplot2::geom_step(linewidth = 1.0, direction = "hv") +
  ggplot2::geom_point(
    data = alpha_five,
    size = 2.8,
    shape = 21,
    fill = "white",
    stroke = 1.0
  ) +
  ggplot2::geom_text(
    data = alpha_five,
    ggplot2::aes(label = label),
    hjust = -0.08,
    vjust = c(-0.7, 1.5),
    size = 3.4,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = colors, name = NULL) +
  ggplot2::scale_x_continuous(
    limits = c(0, 0.20),
    breaks = seq(0, 0.20, by = 0.025),
    labels = scales::label_percent(accuracy = 0.1),
    expand = ggplot2::expansion(mult = c(0, 0.01))
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 0.45),
    breaks = seq(0, 0.40, by = 0.10),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    title = "Breaking donor trajectories removes most known-null inflation",
    subtitle = paste(
      "All-gene random-variant signal-stripped residual controls;",
      "BF-adjusted merged discoveries"
    ),
    x = "Nominal alpha",
    y = "Known-null FDP lower bound: V / R[merged]",
    caption = paste(
      paste(
        "Both curves use the same 6,362 targets, selection seed, and master",
        "permutation seed."
      ),
      paste(
        "This paired fixed-map contrast is a mechanism diagnostic, not a",
        "repeated-sampling FDR estimate."
      ),
      sep = "\n"
    )
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    legend.position = "top",
    panel.grid.minor = ggplot2::element_blank(),
    plot.caption = ggplot2::element_text(hjust = 0, color = "grey35")
  )

figure_directory <- file.path(
  output_parent,
  independent_id,
  "figures"
)
if (!dir.exists(figure_directory)) {
  dir.create(figure_directory, recursive = TRUE)
}
curve_path <- file.path(
  output_parent,
  independent_id,
  "synchronized_vs_independent_time_fdp_lower_bound.csv"
)
png_path <- file.path(
  figure_directory,
  "synchronized_vs_independent_time_fdp_lower_bound.png"
)
pdf_path <- file.path(
  figure_directory,
  "synchronized_vs_independent_time_fdp_lower_bound.pdf"
)
utils::write.csv(curve_data, curve_path, row.names = FALSE, quote = TRUE)
ggplot2::ggsave(png_path, comparison_plot, width = 9.2, height = 6.0, dpi = 220)
ggplot2::ggsave(pdf_path, comparison_plot, width = 9.2, height = 6.0)

cat("Synchronized-versus-independent-time comparison plot created.\n")
cat("Curve data: ", curve_path, "\n", sep = "")
cat("PNG: ", png_path, "\n", sep = "")
cat("PDF: ", pdf_path, "\n", sep = "")
print(alpha_five[, c(
  "method", "nominal_alpha", "target_calls", "permuted_null_calls",
  "merged_calls", "known_null_fdp_lower_bound"
)], row.names = FALSE)
