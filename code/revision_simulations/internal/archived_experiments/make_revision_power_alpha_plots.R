#!/usr/bin/env Rscript

# Generate power-vs-alpha and empirical-FDR-vs-alpha curves from saved targeted
# simulation runs.
#
# This script uses existing RDS outputs and does not refit any models. The
# selection rule is the same lfdr/q-value rule used by the simulation summary:
# discoveries are units with adjusted_score <= alpha.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find code/revision_simulations/shared/simulation_functions.R from the current working directory.")
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

dirs <- revision_output_dirs(file.path(workflowr_root, "output", "revision_simulations"))
alpha_grid <- seq(0, 0.20, by = 0.005)

curve_specs <- list(
  list(
    output_id = "J300_bspline_medium_seed12345",
    sample_size = 300,
    raw_file = "constant_vs_bspline_simplified_sd0p3-0p5-0p8_sigmabeta1_seed12345_J300.rds",
    mechanism_label = "Cubic B-spline random curve",
    se_label = "Medium SE: 0.3, 0.5, 0.8"
  ),
  list(
    output_id = "J300_bspline_high_seed12345",
    sample_size = 300,
    raw_file = "constant_vs_bspline_simplified_sd0p5-0p8-1p2_sigmabeta1_seed12345_J300.rds",
    mechanism_label = "Cubic B-spline random curve",
    se_label = "High SE: 0.5, 0.8, 1.2"
  ),
  list(
    output_id = "J300_spiky_bspline_transient_medium_seed12345",
    sample_size = 300,
    raw_file = "constant_vs_local_bspline_transient_simplified_sd0p3-0p5-0p8_sigmabeta1_amp2_bsdf6_tdf16_width0p8_seed12345_J300.rds",
    mechanism_label = "Spiky local B-spline transient",
    se_label = "Medium SE: 0.3, 0.5, 0.8"
  ),
  list(
    output_id = "J300_spiky_bspline_transient_high_seed12345",
    sample_size = 300,
    raw_file = "constant_vs_local_bspline_transient_simplified_sd0p5-0p8-1p2_sigmabeta1_amp2_bsdf6_tdf16_width0p8_seed12345_J300.rds",
    mechanism_label = "Spiky local B-spline transient",
    se_label = "High SE: 0.5, 0.8, 1.2"
  ),
  list(
    output_id = "J1000_bspline_medium_seed12345",
    sample_size = 1000,
    raw_file = "constant_vs_bspline_simplified_sd0p3-0p5-0p8_sigmabeta1_amp2_bsdf6_tdf16_width0p8_seed12345_J1000.rds",
    mechanism_label = "Cubic B-spline random curve",
    se_label = "Medium SE: 0.3, 0.5, 0.8"
  ),
  list(
    output_id = "J1000_bspline_high_seed12345",
    sample_size = 1000,
    raw_file = "constant_vs_bspline_simplified_sd0p5-0p8-1p2_sigmabeta1_amp2_bsdf6_tdf16_width0p8_seed12345_J1000.rds",
    mechanism_label = "Cubic B-spline random curve",
    se_label = "High SE: 0.5, 0.8, 1.2"
  ),
  list(
    output_id = "J1000_spiky_bspline_transient_medium_seed12345",
    sample_size = 1000,
    raw_file = "constant_vs_local_bspline_transient_simplified_sd0p3-0p5-0p8_sigmabeta1_amp2_bsdf6_tdf16_width0p8_seed12345_J1000.rds",
    mechanism_label = "Spiky local B-spline transient",
    se_label = "Medium SE: 0.3, 0.5, 0.8"
  ),
  list(
    output_id = "J1000_spiky_bspline_transient_high_seed12345",
    sample_size = 1000,
    raw_file = "constant_vs_local_bspline_transient_simplified_sd0p5-0p8-1p2_sigmabeta1_amp2_bsdf6_tdf16_width0p8_seed12345_J1000.rds",
    mechanism_label = "Spiky local B-spline transient",
    se_label = "High SE: 0.5, 0.8, 1.2"
  )
)

summarize_calibration <- function(curve) {
  calibration_summary <- do.call(rbind, by(
    curve[curve$alpha > 0, ],
    list(curve$output_id[curve$alpha > 0], curve$method[curve$alpha > 0]),
    function(d) {
      d <- d[order(d$alpha), ]
      alpha_005_index <- which.min(abs(d$alpha - 0.05))
      excess_fdr <- d$empirical_fdr - d$alpha
      data.frame(
        output_id = d$output_id[1],
        sample_size = d$sample_size[1],
        mechanism_label = d$mechanism_label[1],
        se_label = d$se_label[1],
        mechanism_rank = d$mechanism_rank[1],
        se_rank = d$se_rank[1],
        method = d$method[1],
        fdr_at_005 = d$empirical_fdr[alpha_005_index],
        power_at_005 = d$power[alpha_005_index],
        max_fdr_minus_alpha = max(excess_fdr, na.rm = TRUE),
        alpha_at_max_excess = d$alpha[which.max(excess_fdr)],
        max_empirical_fdr = max(d$empirical_fdr, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  ))
  calibration_summary <- calibration_summary[order(
    calibration_summary$sample_size,
    calibration_summary$mechanism_rank,
    calibration_summary$se_rank,
    rank_revision_methods(calibration_summary$method)
  ), ]
  rownames(calibration_summary) <- NULL
  calibration_summary
}

write_curve_outputs <- function(curve, suffix, figure_title_suffix) {
  combined_csv <- file.path(
    dirs$summary,
    paste0("power_alpha_curve_bspline_and_spiky_transient_", suffix, ".csv")
  )
  write.csv(curve, combined_csv, row.names = FALSE)

  combined_figure <- file.path(
    dirs$figures,
    paste0("power_alpha_curve_bspline_and_spiky_transient_", suffix, ".png")
  )
  panel_var <- if ("panel_label_with_J" %in% colnames(curve) &&
    length(unique(curve$sample_size)) > 1) {
    "panel_label_with_J"
  } else {
    "panel_label"
  }

  plot_power_alpha_curve_grid(
    curve,
    panel_var = panel_var,
    file = combined_figure,
    title = paste("Power curves for B-spline and spiky transient simulations", figure_title_suffix)
  )

  combined_fdr_figure <- file.path(
    dirs$figures,
    paste0("empirical_fdr_alpha_curve_bspline_and_spiky_transient_", suffix, ".png")
  )
  plot_empirical_fdr_alpha_curve_grid(
    curve,
    panel_var = panel_var,
    file = combined_fdr_figure,
    title = paste("Empirical FDR curves for B-spline and spiky transient simulations", figure_title_suffix)
  )

  calibration_summary <- summarize_calibration(curve)
  calibration_csv <- file.path(
    dirs$summary,
    paste0("fdr_calibration_summary_bspline_and_spiky_transient_", suffix, ".csv")
  )
  write.csv(calibration_summary, calibration_csv, row.names = FALSE)

  list(
    curve_csv = combined_csv,
    power_figure = combined_figure,
    fdr_figure = combined_fdr_figure,
    calibration_csv = calibration_csv,
    calibration_summary = calibration_summary
  )
}

all_curves <- lapply(curve_specs, function(spec) {
  raw_path <- file.path(dirs$raw, spec$raw_file)
  if (!file.exists(raw_path)) {
    stop("Missing saved simulation output: ", raw_path)
  }

  simulation_output <- readRDS(raw_path)
  alpha_curve <- compute_alpha_curve(
    result_table = simulation_output$result_table,
    alpha_grid = alpha_grid
  )
  alpha_curve$sample_size <- spec$sample_size
  alpha_curve$mechanism_label <- spec$mechanism_label
  alpha_curve$se_label <- spec$se_label
  alpha_curve$mechanism_rank <- if (grepl("^Cubic", spec$mechanism_label)) 1 else 2
  alpha_curve$se_rank <- if (grepl("^Medium", spec$se_label)) 1 else 2
  alpha_curve$output_id <- spec$output_id
  alpha_curve$panel_label <- paste(
    spec$mechanism_label,
    spec$se_label,
    sep = "\n"
  )
  alpha_curve$panel_label_with_J <- paste(
    paste0("J = ", spec$sample_size),
    spec$mechanism_label,
    spec$se_label,
    sep = "\n"
  )

  csv_path <- file.path(dirs$summary, paste0("power_alpha_curve_", spec$output_id, ".csv"))
  write.csv(alpha_curve, csv_path, row.names = FALSE)

  figure_path <- file.path(dirs$figures, paste0("power_alpha_curve_", spec$output_id, ".png"))
  plot_power_alpha_curves(
    alpha_curve,
    file = figure_path,
    title = paste0("J = ", spec$sample_size, ": ", spec$mechanism_label),
    subtitle = spec$se_label
  )

  fdr_figure_path <- file.path(dirs$figures, paste0("empirical_fdr_alpha_curve_", spec$output_id, ".png"))
  plot_empirical_fdr_alpha_curves(
    alpha_curve,
    file = fdr_figure_path,
    title = paste0("J = ", spec$sample_size, ": ", spec$mechanism_label),
    subtitle = spec$se_label
  )

  alpha_curve
})

combined_curve <- do.call(rbind, all_curves)
combined_curve <- combined_curve[order(
  combined_curve$sample_size,
  combined_curve$mechanism_rank,
  combined_curve$se_rank,
  rank_revision_methods(combined_curve$method),
  combined_curve$alpha
), ]
rownames(combined_curve) <- NULL

all_outputs <- write_curve_outputs(
  combined_curve,
  suffix = "J300_J1000_seed12345",
  figure_title_suffix = "(J = 300 and J = 1000)"
)

j1000_curve <- combined_curve[combined_curve$sample_size == 1000, ]
j1000_outputs <- write_curve_outputs(
  j1000_curve,
  suffix = "J1000_seed12345",
  figure_title_suffix = "(J = 1000)"
)

alpha_005 <- combined_curve[abs(combined_curve$alpha - 0.05) < 1e-12, ]
alpha_005 <- alpha_005[order(
  alpha_005$sample_size,
  alpha_005$mechanism_label,
  alpha_005$se_label,
  rank_revision_methods(alpha_005$method)
), c(
  "output_id",
  "sample_size",
  "method",
  "alpha",
  "n_discoveries",
  "false_discoveries",
  "empirical_fdr",
  "power"
)]

print(alpha_005, row.names = FALSE)
print(j1000_outputs$calibration_summary[, c(
  "output_id",
  "method",
  "fdr_at_005",
  "max_fdr_minus_alpha",
  "alpha_at_max_excess"
)], row.names = FALSE)
message("Wrote all-J combined CSV: ", all_outputs$curve_csv)
message("Wrote J=1000 combined CSV: ", j1000_outputs$curve_csv)
message("Wrote J=1000 power figure: ", j1000_outputs$power_figure)
message("Wrote J=1000 empirical FDR figure: ", j1000_outputs$fdr_figure)
message("Wrote J=1000 calibration summary CSV: ", j1000_outputs$calibration_csv)
