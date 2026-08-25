#!/usr/bin/env Rscript

# Figures for the residual-based common-correlation experiment.
#
# Produces, into output/.../residual_correlation_fash/figures/:
#   residual_correlation_heatmaps.png  the estimated matrix beside its
#                                      comparators and its calibration null
#   residual_correlation_excess.png    observed minus calibration null
#   residual_correlation_lag.png       lag profiles on one axis
#   residual_z_scale.png               the shrinkage deflation of the z scale

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

suppressPackageStartupMessages({
  library(ggplot2)
})

workflowr_root <- find_workflowr_root()
cache_directory <- file.path(
  workflowr_root, "output", "revision_simulations", "internal",
  "residual_correlation_fash"
)
summary_directory <- file.path(cache_directory, "summary")
figure_directory <- file.path(cache_directory, "figures")
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

read_required_csv <- function(name) {
  path <- file.path(summary_directory, name)
  if (!file.exists(path)) {
    stop("Missing required summary file: ", path)
  }
  utils::read.csv(path, stringsAsFactors = FALSE)
}

matrices <- readRDS(file.path(cache_directory, "correlation_matrices.rds"))
analytic_path <- file.path(cache_directory, "analytic_calibration.rds")
analytic <- if (file.exists(analytic_path)) readRDS(analytic_path) else NULL
deconvolution_path <- file.path(cache_directory, "deconvolution.rds")
deconvolution <- if (file.exists(deconvolution_path)) {
  readRDS(deconvolution_path)
} else {
  NULL
}
calibration_path <- file.path(cache_directory, "calibration_replicates.rds")
calibration <- if (file.exists(calibration_path)) {
  readRDS(calibration_path)
} else {
  NULL
}

# The design-based donor-residual-block permutation null from 2026-08-10 is a
# completely independent estimate of the same quantity: it uses genotype and
# expression data and 400 permutations, not summary-statistic residuals. It is
# the external validation for the deconvolved matrix.
permutation_path <- file.path(
  workflowr_root, "output", "revision_simulations", "internal",
  "multigene_null_beta_covariance", "summary",
  "primary_donor_residual_block_C_m_g_gt_0p900.csv"
)
permutation <- if (file.exists(permutation_path)) {
  frame <- utils::read.csv(permutation_path, stringsAsFactors = FALSE)
  as.matrix(frame[, grep("^time_", names(frame))])
} else {
  NULL
}

to_long <- function(matrix, label) {
  matrix <- as.matrix(matrix)
  n_time <- nrow(matrix)
  frame <- expand.grid(row = seq_len(n_time), column = seq_len(n_time))
  frame$value <- matrix[cbind(frame$row, frame$column)]
  frame$time_a <- frame$row - 1L
  frame$time_b <- frame$column - 1L
  frame$panel <- label
  frame
}

# Ordered so the story reads left to right: what you can measure directly,
# then what the two independence calibrations say the estimator must return
# anyway, then the deconvolved answer beside its independent validation.
panels <- list(
  to_long(matrices$observed_raw_z, "1. Observed raw z\n(no signal removal)"),
  if (!is.null(matrices$C1)) {
    to_long(matrices$C1, "2. C1 null screen\n(2026-08-07)")
  },
  to_long(matrices$bf_pearson, "3. Residual z, naive\n(this experiment)"),
  if (!is.null(analytic)) {
    to_long(analytic$correlation,
            "4. Independence null\n(analytic shrinkage artifact)")
  },
  if (!is.null(calibration) && length(calibration) > 0L) {
    to_long(calibration[[1]]$correlation,
            "5. Independence null\n(bootstrap refit)")
  },
  if (!is.null(deconvolution)) {
    to_long(deconvolution$solutions$unregularised$correlation,
            "6. Deconvolved error correlation\n(shrinkage inverted)")
  },
  if (!is.null(permutation)) {
    to_long(permutation,
            "7. Design permutation null\n(2026-08-10, independent)")
  }
)
panels <- panels[!vapply(panels, is.null, logical(1))]
heatmap_data <- do.call(rbind, panels)
heatmap_data$panel <- factor(heatmap_data$panel,
                             levels = unique(heatmap_data$panel))

limit <- max(abs(heatmap_data$value[heatmap_data$time_a != heatmap_data$time_b]))
heatmap_plot <- ggplot(heatmap_data, aes(time_b, time_a, fill = value)) +
  geom_raster() +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B",
    midpoint = 0, limits = c(-limit, limit), name = "Correlation",
    oob = scales::squish
  ) +
  scale_y_reverse(breaks = seq(0, 15, by = 3)) +
  scale_x_continuous(breaks = seq(0, 15, by = 3)) +
  coord_equal() +
  facet_wrap(~panel, ncol = 4) +
  labs(
    title = "Common cross-time correlation of dynamic-eQTL summary statistics",
    subtitle = paste0(
      "One variant per gene, 6,362 pairs, thinning seed 20260811. ",
      "Panels 6 and 7 are independent estimates of the same quantity.\n",
      "The diagonal is 1 by construction and is clipped to the top of the ",
      "colour scale, which is set by the largest off-diagonal entry."
    ),
    x = "Day", y = "Day"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 8.5, lineheight = 1.05),
    panel.grid = element_blank()
  )
ggsave(file.path(figure_directory, "residual_correlation_heatmaps.png"),
       heatmap_plot, width = 14, height = 7.5, dpi = 200)

if (!is.null(analytic)) {
  excess <- matrices$bf_pearson - analytic$correlation
  excess_data <- to_long(excess, "Observed minus calibration null")
  excess_limit <- max(abs(
    excess_data$value[excess_data$time_a != excess_data$time_b]
  ))
  excess_plot <- ggplot(excess_data, aes(time_b, time_a, fill = value)) +
    geom_raster() +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-excess_limit, excess_limit), name = "Excess",
      oob = scales::squish
    ) +
    scale_y_reverse(breaks = seq(0, 15, by = 3)) +
    scale_x_continuous(breaks = seq(0, 15, by = 3)) +
    coord_equal() +
    labs(
      title = "Residual correlation not explained by FASH shrinkage",
      subtitle = paste(
        "Observed residual-z correlation minus the exact analytic",
        "independence calibration.\nOnly this excess could be genuine",
        "error correlation."
      ),
      x = "Day", y = "Day"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), panel.grid = element_blank())
  ggsave(file.path(figure_directory, "residual_correlation_excess.png"),
         excess_plot, width = 6, height = 5.5, dpi = 200)
}

lag_frames <- list(read_required_csv("stage_c_lag_profiles.csv"))
if (file.exists(file.path(summary_directory, "stage_g_lag_profiles.csv"))) {
  lag_frames[[length(lag_frames) + 1L]] <-
    read_required_csv("stage_g_lag_profiles.csv")
}
if (file.exists(file.path(summary_directory, "stage_f_lag_profiles.csv"))) {
  lag_frames[[length(lag_frames) + 1L]] <-
    read_required_csv("stage_f_lag_profiles.csv")
}
if (file.exists(file.path(summary_directory, "stage_d_calibration_lag.csv"))) {
  lag_frames[[length(lag_frames) + 1L]] <-
    read_required_csv("stage_d_calibration_lag.csv")
}
lag_data <- do.call(rbind, lag_frames)
if (!is.null(permutation)) {
  permutation_lag <- data.frame(
    name = "permutation_null",
    lag = seq_len(nrow(permutation) - 1L),
    mean_correlation = vapply(seq_len(nrow(permutation) - 1L), function(lag) {
      index <- seq_len(nrow(permutation) - lag)
      mean(permutation[cbind(index, index + lag)])
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  lag_data <- rbind(lag_data[, names(permutation_lag)], permutation_lag)
}
keep <- c("bf_pearson", "observed_raw_z", "C1", "C2",
          "analytic_independence", "deconvolved_unregularised",
          "permutation_null",
          grep("^independence_seed", unique(lag_data$name), value = TRUE))
lag_data <- lag_data[lag_data$name %in% keep, ]
lag_data <- lag_data[!duplicated(lag_data[, c("name", "lag")]), ]
label_map <- c(
  bf_pearson = "Residual z (this experiment)",
  observed_raw_z = "Observed raw z",
  C1 = "C1 null screen",
  C2 = "C2 null screen",
  analytic_independence = "Calibration null (analytic)",
  deconvolved_unregularised = "Deconvolved error correlation",
  permutation_null = "Design-based permutation null (independent)"
)
lag_data$label <- ifelse(
  lag_data$name %in% names(label_map),
  label_map[lag_data$name],
  "Calibration null (bootstrap)"
)
lag_plot <- ggplot(lag_data, aes(lag, mean_correlation,
                                 colour = label, linetype = label)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  scale_x_continuous(breaks = 1:15) +
  labs(
    title = "Mean correlation by day separation",
    subtitle = paste(
      "The raw residual estimator sits far below the truth because FASH",
      "shrinkage distorts it. Inverting that\noperator recovers a profile",
      "close to the independent design-based permutation estimate."
    ),
    x = "Lag in days", y = "Mean correlation",
    colour = NULL, linetype = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )
ggsave(file.path(figure_directory, "residual_correlation_lag.png"),
       lag_plot, width = 9, height = 6, dpi = 200)

if (file.exists(file.path(summary_directory, "stage_f_z_scale.csv"))) {
  scale_data <- read_required_csv("stage_f_z_scale.csv")
  long_scale <- rbind(
    data.frame(time = scale_data$time, value = scale_data$observed_z_sd,
               series = "Observed residual z"),
    data.frame(time = scale_data$time,
               value = scale_data$expected_z_sd_under_independence,
               series = "Predicted under independence")
  )
  scale_plot <- ggplot(long_scale, aes(time, value, colour = series)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    scale_x_continuous(breaks = seq(0, 15, by = 1)) +
    expand_limits(y = c(0, 1.05)) +
    labs(
      title = "The residual z scale is deflated by shrinkage, not equal to one",
      subtitle = paste(
        "The dashed line is the scale the residual would have if the",
        "posterior mean were a fixed known truth."
      ),
      x = "Day", y = "Standard deviation of residual z", colour = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(figure_directory, "residual_z_scale.png"),
         scale_plot, width = 8, height = 5, dpi = 200)
}

if (!is.null(permutation) && !is.null(deconvolution)) {
  candidates <- list(
    deconvolved = deconvolution$solutions$unregularised$correlation,
    deconvolved_ridge_1e2 = deconvolution$solutions$ridge_1e2$correlation,
    naive_residual_z = matrices$bf_pearson,
    observed_raw_z = matrices$observed_raw_z
  )
  if (!is.null(matrices$C1)) candidates$C1 <- matrices$C1
  if (!is.null(matrices$C2)) candidates$C2 <- matrices$C2
  upper <- upper.tri(permutation)
  agreement <- do.call(rbind, lapply(names(candidates), function(name) {
    candidate <- as.matrix(candidates[[name]])
    data.frame(
      estimator = name,
      mean_absolute_difference = mean(abs(candidate[upper] -
                                            permutation[upper])),
      maximum_absolute_difference = max(abs(candidate[upper] -
                                              permutation[upper])),
      offdiagonal_pearson = stats::cor(candidate[upper], permutation[upper]),
      offdiagonal_spearman = stats::cor(candidate[upper], permutation[upper],
                                        method = "spearman"),
      mean_offdiagonal = mean(candidate[upper]),
      lag1 = mean(candidate[cbind(seq_len(nrow(candidate) - 1L),
                                  seq_len(nrow(candidate) - 1L) + 1L)]),
      stringsAsFactors = FALSE
    )
  }))
  agreement <- agreement[order(agreement$mean_absolute_difference), ]
  agreement$permutation_mean_offdiagonal <- mean(permutation[upper])
  agreement$permutation_lag1 <- mean(permutation[cbind(
    seq_len(nrow(permutation) - 1L), seq_len(nrow(permutation) - 1L) + 1L
  )])
  utils::write.csv(
    agreement,
    file.path(summary_directory, "cross_experiment_agreement.csv"),
    row.names = FALSE
  )
  cat("\nAgreement with the independent 2026-08-10 permutation null:\n")
  print(agreement, digits = 3, row.names = FALSE)
}

cat("Wrote figures to", normalizePath(figure_directory), "\n")
print(list.files(figure_directory))

# ---------------------------------------------------------------------------
# Readable, number-annotated view of the deconvolved matrix beside the
# independent 2026-08-10 permutation estimate.
# ---------------------------------------------------------------------------
if (!is.null(deconvolution)) {
  annotated <- list(
    to_long(deconvolution$solutions$unregularised$correlation,
            "Deconvolved from FASH residuals (this experiment)")
  )
  if (!is.null(permutation)) {
    annotated[[2]] <- to_long(
      permutation,
      "Donor-residual-block permutation (2026-08-10, independent)"
    )
  }
  annotated_data <- do.call(rbind, annotated)
  annotated_data$panel <- factor(annotated_data$panel,
                                 levels = unique(annotated_data$panel))
  annotated_data$label <- ifelse(
    annotated_data$time_a == annotated_data$time_b, "",
    sub("^0", "", formatC(annotated_data$value, format = "f", digits = 2))
  )
  offdiag_limit <- max(abs(
    annotated_data$value[annotated_data$time_a != annotated_data$time_b]
  ))
  annotated_plot <- ggplot(annotated_data, aes(time_b, time_a, fill = value)) +
    geom_raster() +
    geom_text(aes(label = label), size = 2.1, colour = "grey15") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      limits = c(-offdiag_limit, offdiag_limit), name = "Correlation",
      oob = scales::squish
    ) +
    scale_y_reverse(breaks = 0:15) +
    scale_x_continuous(breaks = 0:15, position = "top") +
    coord_equal() +
    facet_wrap(~panel, ncol = 2) +
    labs(
      title = "Estimated cross-time error correlation of the eQTL effect estimates",
      subtitle = paste0(
        "Two estimators sharing no data path. Left: invert the FASH shrinkage ",
        "operator on residual z-scores from 6,362\nrandomly thinned ",
        "gene-variant pairs. Right: 400 donor-residual-block permutations ",
        "under beta = 0 on genotype and\nexpression data. Diagonal is 1 and ",
        "is left unlabelled."
      ),
      x = "Day", y = "Day"
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold"),
      panel.grid = element_blank()
    )
  ggsave(file.path(figure_directory, "deconvolved_matrix_annotated.png"),
         annotated_plot, width = 15, height = 8, dpi = 200)
  cat("Wrote deconvolved_matrix_annotated.png\n")
}

# ---------------------------------------------------------------------------
# Shape versus level. Two matrices can agree on the decay pattern and still
# disagree badly on how strong the correlation is. Normalising each lag profile
# by its own lag-1 value separates the two questions.
# ---------------------------------------------------------------------------
if (!is.null(deconvolution) && !is.null(permutation)) {
  shape_set <- list(
    "Naive residual z"          = matrices$bf_pearson,
    "Deconvolved"               = deconvolution$solutions$unregularised$correlation,
    "C1 null screen"            = matrices$C1,
    "C2 null screen"            = matrices$C2,
    "Design permutation (ref.)" = permutation,
    "Observed raw z"            = matrices$observed_raw_z
  )
  shape_set <- shape_set[!vapply(shape_set, is.null, logical(1))]
  lag_of <- function(x) {
    x <- as.matrix(x)
    n <- nrow(x)
    vapply(seq_len(n - 1L), function(k) {
      mean(x[cbind(seq_len(n - k), seq_len(n - k) + k)])
    }, numeric(1))
  }
  shape_frame <- do.call(rbind, lapply(names(shape_set), function(nm) {
    profile <- lag_of(shape_set[[nm]])
    data.frame(estimator = nm, lag = seq_along(profile),
               normalised = profile / profile[1], stringsAsFactors = FALSE)
  }))
  shape_frame$estimator <- factor(shape_frame$estimator,
                                  levels = names(shape_set))
  shape_plot <- ggplot(shape_frame,
                       aes(lag, normalised, colour = estimator,
                           linetype = estimator)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.6) +
    scale_x_continuous(breaks = 1:15) +
    labs(
      title = "Decay shape only: every lag profile divided by its own lag-1 value",
      subtitle = paste0(
        "Curves that overlap have the same shape and differ only in overall ",
        "level. The deconvolved, C1, permutation and\nraw-z estimators form ",
        "one family; the naive residual and C2 estimators cross zero, which ",
        "no rescaling of a\npositive matrix can produce."
      ),
      x = "Lag in days", y = "Mean correlation / lag-1 correlation",
      colour = NULL, linetype = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
  ggsave(file.path(figure_directory, "correlation_shape_comparison.png"),
         shape_plot, width = 9, height = 6, dpi = 200)
  cat("Wrote correlation_shape_comparison.png\n")
}
