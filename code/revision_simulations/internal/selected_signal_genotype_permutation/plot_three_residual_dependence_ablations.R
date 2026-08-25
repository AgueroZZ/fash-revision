#!/usr/bin/env Rscript

# Compare the three signal-stripped residual dependence ablations.

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
experiment_ids <- c(
  "Shared map; trajectories preserved" = paste0(
    "all_gene_random_variant_signal_stripped_residual_block_permutation_",
    "selection20260817_seed20260811"
  ),
  "Map per unit; trajectories preserved" = paste0(
    "all_gene_random_variant_signal_stripped_unit_specific_residual_block_",
    "permutation_selection20260817_seed20260811"
  ),
  "Map per time; trajectories broken" = paste0(
    "all_gene_random_variant_signal_stripped_independent_time_residual_",
    "permutation_selection20260817_seed20260811"
  )
)

selection_tables <- lapply(experiment_ids, function(experiment_id) {
  path <- file.path(output_parent, experiment_id, "selection.csv")
  if (!file.exists(path)) {
    stop("A required selection table is missing: ", path)
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
})
if (!all(vapply(selection_tables[-1L], function(selection) {
  identical(selection, selection_tables[[1L]])
}, logical(1)))) {
  stop("The three experiments do not use the exact same target selection.")
}

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
alpha_grids <- lapply(levels(curve_data$method), function(method) {
  curve_data$nominal_alpha[curve_data$method == method]
})
if (!all(required_columns %in% names(curve_data)) ||
    any(!is.finite(curve_data$nominal_alpha)) ||
    any(!is.finite(curve_data$known_null_fdp_lower_bound)) ||
    !all(vapply(alpha_grids[-1L], identical, logical(1), alpha_grids[[1L]]))) {
  stop("The three calibration curves failed structural validation.")
}

alpha_five <- curve_data[curve_data$nominal_alpha == 0.05, , drop = FALSE]
if (nrow(alpha_five) != length(experiment_ids)) {
  stop("The three curves do not each contain nominal alpha 0.05.")
}
alpha_five$label <- paste0(
  scales::percent(alpha_five$known_null_fdp_lower_bound, accuracy = 0.1),
  " (",
  alpha_five$permuted_null_calls,
  "/",
  alpha_five$merged_calls,
  ")"
)
alpha_five$label_y <- alpha_five$known_null_fdp_lower_bound +
  c(0.014, -0.018, -0.014)

colors <- c(
  "Shared map; trajectories preserved" = "#D55E00",
  "Map per unit; trajectories preserved" = "#009E73",
  "Map per time; trajectories broken" = "#0072B2"
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
    ggplot2::aes(y = label_y, label = label),
    hjust = -0.08,
    size = 3.3,
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
    title = "Known-null inflation tracks preserved temporal trajectories",
    subtitle = paste(
      "All-gene random-variant signal-stripped residual controls;",
      "BF-adjusted merged discoveries"
    ),
    x = "Nominal alpha",
    y = "Known-null FDP lower bound: V / R[merged]",
    caption = paste(
      paste(
        "All experiments use the same 6,362 targets, target likelihoods,",
        "selection seed, and master permutation seed."
      ),
      paste(
        "Each curve is a fixed-randomization mechanism diagnostic, not a",
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

unit_specific_id <- experiment_ids[[
  "Map per unit; trajectories preserved"
]]
figure_directory <- file.path(output_parent, unit_specific_id, "figures")
if (!dir.exists(figure_directory)) {
  dir.create(figure_directory, recursive = TRUE)
}
curve_path <- file.path(
  output_parent,
  unit_specific_id,
  "three_residual_dependence_ablations.csv"
)
png_path <- file.path(
  figure_directory,
  "three_residual_dependence_ablations.png"
)
pdf_path <- file.path(
  figure_directory,
  "three_residual_dependence_ablations.pdf"
)
utils::write.csv(curve_data, curve_path, row.names = FALSE, quote = TRUE)
ggplot2::ggsave(png_path, comparison_plot, width = 10.0, height = 6.2, dpi = 220)
ggplot2::ggsave(pdf_path, comparison_plot, width = 10.0, height = 6.2)

cat("Three-design residual-dependence comparison plot created.\n")
cat("Curve data: ", curve_path, "\n", sep = "")
cat("PNG: ", png_path, "\n", sep = "")
cat("PDF: ", pdf_path, "\n", sep = "")
print(alpha_five[, c(
  "method", "nominal_alpha", "target_calls", "permuted_null_calls",
  "merged_calls", "known_null_fdp_lower_bound"
)], row.names = FALSE)
