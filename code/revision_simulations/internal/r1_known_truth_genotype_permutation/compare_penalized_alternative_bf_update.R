#!/usr/bin/env Rscript

# Compare the current FASH BF update with the manuscript-specified update that
# retains the conditional alternative weights from the penalized raw fit.

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

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

summarize_mc_curve <- function(curve) {
  keys <- unique(curve[, c("experiment", "arm", "method", "alpha")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    selected <- curve[
      curve$experiment == key$experiment &
        curve$arm == key$arm &
        curve$method == key$method &
        abs(curve$alpha - key$alpha) < 1e-12,
      , drop = FALSE
    ]
    n_replications <- nrow(selected)
    summarize_metric <- function(values) {
      mean_value <- mean(values)
      standard_error <- if (n_replications > 1L) {
        stats::sd(values) / sqrt(n_replications)
      } else {
        NA_real_
      }
      critical_value <- if (n_replications > 1L) {
        stats::qt(0.975, df = n_replications - 1L)
      } else {
        NA_real_
      }
      c(
        mean = mean_value,
        se = standard_error,
        lower = mean_value - critical_value * standard_error,
        upper = mean_value + critical_value * standard_error
      )
    }
    fdr <- summarize_metric(selected$realized_fdp)
    power <- summarize_metric(selected$power)
    data.frame(
      experiment = key$experiment,
      arm = key$arm,
      method = key$method,
      alpha = key$alpha,
      n_replications = n_replications,
      mean_fdr = fdr[["mean"]],
      fdr_se = fdr[["se"]],
      fdr_ci_lower = max(0, fdr[["lower"]]),
      fdr_ci_upper = min(1, fdr[["upper"]]),
      mean_power = power[["mean"]],
      power_se = power[["se"]],
      power_ci_lower = max(0, power[["lower"]]),
      power_ci_upper = min(1, power[["upper"]]),
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, rows)
  output[order(
    output$experiment, output$arm, output$method, output$alpha
  ), , drop = FALSE]
}

summarize_pooled_null_bf <- function(unit_results) {
  null_units <- unit_results[unit_results$true_null, , drop = FALSE]
  keys <- unique(null_units[, c("experiment", "arm", "method")])
  rows <- lapply(seq_len(nrow(keys)), function(i) {
    key <- keys[i, , drop = FALSE]
    selected <- null_units[
      null_units$experiment == key$experiment &
        null_units$arm == key$arm &
        null_units$method == key$method,
      , drop = FALSE
    ]
    values <- selected$bayes_factor
    quantiles <- stats::quantile(
      values, probs = c(0.50, 0.90, 0.95, 0.99), type = 8,
      names = FALSE
    )
    descending <- order(values, decreasing = TRUE)
    top_one_percent <- ceiling(0.01 * length(values))
    data.frame(
      experiment = key$experiment,
      arm = key$arm,
      method = key$method,
      n_seeds = length(unique(selected$seed)),
      n_null = length(values),
      mean_bf = mean(values),
      median_bf = quantiles[1L],
      q90_bf = quantiles[2L],
      q95_bf = quantiles[3L],
      q99_bf = quantiles[4L],
      maximum_bf = max(values),
      proportion_bf_greater_than_one = mean(values > 1),
      top_1_percent_bf_mass_share =
        sum(values[descending[seq_len(top_one_percent)]]) / sum(values),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation",
  "selected_signal_genotype_permutation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "r1_known_truth_genotype_permutation",
  "r1_known_truth_genotype_permutation_helpers.R"
))

if (!requireNamespace("fashr", quietly = TRUE) ||
    !requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("scales", quietly = TRUE)) {
  stop("The fashr, ggplot2, and scales packages are required.")
}

input_specs <- data.frame(
  experiment = c("genotype_permutation", "residual_permutation"),
  output_id = c(
    "r1_known_truth_genotype_permutation_mc5_J200_v1",
    "r1_signal_stripped_unadjusted_residual_permutation_mc5_J200_v1"
  ),
  stringsAsFactors = FALSE
)
output_id <- "r1_penalized_alternative_bf_comparison_v1"
output_directory <- file.path(
  workflowr_root, "output", "revision_simulations", "internal", output_id
)
summary_directory <- file.path(output_directory, "summary")
figure_directory <- file.path(output_directory, "figures")
dir.create(summary_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

alpha_grid <- seq(0, 0.20, by = 0.001)
unit_rows <- list()
replicate_rows <- list()
curve_rows <- list()
validation_rows <- list()
result_index <- 0L

for (spec_index in seq_len(nrow(input_specs))) {
  experiment <- input_specs$experiment[spec_index]
  input_directory <- file.path(
    workflowr_root, "output", "revision_simulations", "internal",
    input_specs$output_id[spec_index], "summary"
  )
  unit_path <- file.path(input_directory, "replicate_unit_results.csv")
  pi0_path <- file.path(input_directory, "replicate_pi0.csv")
  if (!file.exists(unit_path) || !file.exists(pi0_path)) {
    stop("A required cached five-seed summary is missing for ", experiment)
  }
  units <- utils::read.csv(unit_path, stringsAsFactors = FALSE)
  pi0 <- utils::read.csv(pi0_path, stringsAsFactors = FALSE)
  seeds <- sort(unique(units$seed))

  for (seed in seeds) {
    arms <- unique(units$arm[units$seed == seed])
    for (arm in arms) {
      raw <- units[
        units$seed == seed & units$arm == arm & units$fit_stage == "Raw",
        , drop = FALSE
      ]
      current <- units[
        units$seed == seed & units$arm == arm & units$fit_stage == "BF",
        , drop = FALSE
      ]
      current <- current[match(raw$unit_key, current$unit_key), , drop = FALSE]
      raw_pi0 <- pi0$estimated_pi0[
        pi0$seed == seed & pi0$arm == arm & pi0$fit_stage == "Raw"
      ]
      current_pi0 <- pi0$estimated_pi0[
        pi0$seed == seed & pi0$arm == arm & pi0$fit_stage == "BF"
      ]
      if (nrow(raw) != 400L || nrow(current) != 400L ||
          length(raw_pi0) != 1L || length(current_pi0) != 1L ||
          anyNA(current$unit_key) ||
          !identical(raw$unit_key, current$unit_key) ||
          !identical(raw$true_null, current$true_null) ||
          !is.finite(raw_pi0) || raw_pi0 <= 0 || raw_pi0 >= 1 ||
          !is.finite(current_pi0) || current_pi0 <= 0 ||
          current_pi0 >= 1 || any(raw$lfdr <= 0 | raw$lfdr >= 1) ||
          any(!is.finite(current$bayes_factor)) ||
          any(current$bayes_factor <= 0)) {
        stop("A cached seed-arm block is invalid: ", seed, " / ", arm)
      }

      current_lfdr_reconstructed <- stats::plogis(
        stats::qlogis(current_pi0) - log(current$bayes_factor)
      )
      current_lfdr_error <- max(abs(
        current_lfdr_reconstructed - current$lfdr
      ))

      penalized_bf <- raw_pi0 / (1 - raw_pi0) *
        (1 - raw$lfdr) / raw$lfdr
      penalized_pi0 <- fashr:::BF_control(
        penalized_bf, plot = FALSE
      )$pi0_hat_star
      penalized_lfdr <- stats::plogis(
        stats::qlogis(penalized_pi0) - log(penalized_bf)
      )

      methods <- list(
        current_unpenalized_alternative = list(
          bayes_factor = current$bayes_factor,
          lfdr = current$lfdr,
          bf_pi0 = current_pi0
        ),
        corrected_penalized_alternative = list(
          bayes_factor = penalized_bf,
          lfdr = penalized_lfdr,
          bf_pi0 = penalized_pi0
        )
      )

      for (method in names(methods)) {
        result_index <- result_index + 1L
        method_result <- methods[[method]]
        unit_result <- data.frame(
          experiment = experiment,
          seed = seed,
          arm = arm,
          method = method,
          unit_key = raw$unit_key,
          true_null = raw$true_null,
          raw_penalized_pi0 = raw_pi0,
          bf_adjusted_pi0 = method_result$bf_pi0,
          bayes_factor = method_result$bayes_factor,
          lfdr = method_result$lfdr,
          stringsAsFactors = FALSE
        )
        unit_rows[[result_index]] <- unit_result

        curve <- known_truth_alpha_curve(
          lfdr = method_result$lfdr,
          true_null = raw$true_null,
          alpha_grid = alpha_grid,
          arm = arm,
          fit_stage = method
        )
        curve$experiment <- experiment
        curve$seed <- seed
        names(curve)[names(curve) == "fit_stage"] <- "method"
        curve_rows[[result_index]] <- curve

        null_summary <- summarize_null_bf(
          bf = method_result$bayes_factor,
          true_null = raw$true_null,
          arm = arm
        )
        alpha005 <- curve[abs(curve$alpha - 0.05) < 1e-12, , drop = FALSE]
        replicate_rows[[result_index]] <- data.frame(
          experiment = experiment,
          seed = seed,
          arm = arm,
          method = method,
          raw_penalized_pi0 = raw_pi0,
          bf_adjusted_pi0 = method_result$bf_pi0,
          null_summary[, setdiff(names(null_summary), "arm")],
          alpha005_n_discoveries = alpha005$n_discoveries,
          alpha005_false_discoveries = alpha005$false_discoveries,
          alpha005_realized_fdp = alpha005$realized_fdp,
          alpha005_power = alpha005$power,
          stringsAsFactors = FALSE
        )
      }

      validation_rows[[length(validation_rows) + 1L]] <- data.frame(
        experiment = experiment,
        seed = seed,
        arm = arm,
        check = "Current BF-stage lfdr reconstructs from stored BF and pi0",
        maximum_absolute_error = current_lfdr_error,
        tolerance = 1e-12,
        pass = current_lfdr_error <= 1e-12,
        stringsAsFactors = FALSE
      )
    }
  }
}

unit_results <- do.call(rbind, unit_rows)
replicate_summary <- do.call(rbind, replicate_rows)
alpha_curves <- do.call(rbind, curve_rows)
validation <- do.call(rbind, validation_rows)
if (!all(validation$pass)) {
  stop("At least one cached current-BF reconstruction check failed.")
}

mc_curve <- summarize_mc_curve(alpha_curves)
pooled_null_bf <- summarize_pooled_null_bf(unit_results)

current_summary <- replicate_summary[
  replicate_summary$method == "current_unpenalized_alternative",
  , drop = FALSE
]
corrected_summary <- replicate_summary[
  replicate_summary$method == "corrected_penalized_alternative",
  , drop = FALSE
]
paired_keys <- c("experiment", "seed", "arm")
paired <- merge(
  current_summary, corrected_summary,
  by = paired_keys, suffixes = c("_current", "_corrected"), sort = TRUE
)
paired$delta_bf_adjusted_pi0 <-
  paired$bf_adjusted_pi0_corrected - paired$bf_adjusted_pi0_current
paired$delta_mean_null_bf <- paired$mean_bf_corrected - paired$mean_bf_current
paired$delta_alpha005_realized_fdp <-
  paired$alpha005_realized_fdp_corrected -
    paired$alpha005_realized_fdp_current
paired$delta_alpha005_power <-
  paired$alpha005_power_corrected - paired$alpha005_power_current

rank_selection_rows <- lapply(seq_len(nrow(paired)), function(i) {
  key <- paired[i, paired_keys, drop = FALSE]
  current_units <- unit_results[
    unit_results$experiment == key$experiment &
      unit_results$seed == key$seed &
      unit_results$arm == key$arm &
      unit_results$method == "current_unpenalized_alternative",
    , drop = FALSE
  ]
  corrected_units <- unit_results[
    unit_results$experiment == key$experiment &
      unit_results$seed == key$seed &
      unit_results$arm == key$arm &
      unit_results$method == "corrected_penalized_alternative",
    , drop = FALSE
  ]
  corrected_units <- corrected_units[
    match(current_units$unit_key, corrected_units$unit_key), , drop = FALSE
  ]
  current_curve <- alpha_curves[
    alpha_curves$experiment == key$experiment &
      alpha_curves$seed == key$seed &
      alpha_curves$arm == key$arm &
      alpha_curves$method == "current_unpenalized_alternative",
    , drop = FALSE
  ]
  corrected_curve <- alpha_curves[
    alpha_curves$experiment == key$experiment &
      alpha_curves$seed == key$seed &
      alpha_curves$arm == key$arm &
      alpha_curves$method == "corrected_penalized_alternative",
    , drop = FALSE
  ]
  corrected_curve <- corrected_curve[
    match(current_curve$alpha, corrected_curve$alpha), , drop = FALSE
  ]
  current_calls <- current_units$unit_key[cumulative_lfdr_calls(
    current_units$lfdr, alpha = 0.05
  )]
  corrected_calls <- corrected_units$unit_key[cumulative_lfdr_calls(
    corrected_units$lfdr, alpha = 0.05
  )]
  data.frame(
    experiment = key$experiment,
    seed = key$seed,
    arm = key$arm,
    spearman_log_bf = stats::cor(
      log(current_units$bayes_factor),
      log(corrected_units$bayes_factor),
      method = "spearman"
    ),
    alpha005_discovery_jaccard =
      length(intersect(current_calls, corrected_calls)) /
        length(union(current_calls, corrected_calls)),
    maximum_absolute_fdp_difference = max(abs(
      corrected_curve$realized_fdp - current_curve$realized_fdp
    )),
    maximum_absolute_discovery_count_difference = max(abs(
      corrected_curve$n_discoveries - current_curve$n_discoveries
    )),
    stringsAsFactors = FALSE
  )
})
rank_selection_summary <- do.call(rbind, rank_selection_rows)

comparison_metrics <- c(
  "delta_bf_adjusted_pi0", "delta_mean_null_bf",
  "delta_alpha005_realized_fdp", "delta_alpha005_power"
)
paired_groups <- unique(paired[, c("experiment", "arm")])
paired_mc_rows <- lapply(seq_len(nrow(paired_groups)), function(i) {
  key <- paired_groups[i, , drop = FALSE]
  selected <- paired[
    paired$experiment == key$experiment & paired$arm == key$arm,
    , drop = FALSE
  ]
  output <- data.frame(
    experiment = key$experiment,
    arm = key$arm,
    n_seeds = nrow(selected),
    stringsAsFactors = FALSE
  )
  for (metric in comparison_metrics) {
    output[[paste0("mean_", metric)]] <- mean(selected[[metric]])
    output[[paste0("min_", metric)]] <- min(selected[[metric]])
    output[[paste0("max_", metric)]] <- max(selected[[metric]])
  }
  output
})
paired_mc_summary <- do.call(rbind, paired_mc_rows)

write_csv(input_specs, file.path(output_directory, "input_specs.csv"))
write_csv(unit_results, file.path(summary_directory, "unit_results.csv"))
write_csv(
  replicate_summary,
  file.path(summary_directory, "replicate_method_summary.csv")
)
write_csv(
  alpha_curves,
  file.path(summary_directory, "replicate_alpha_curves.csv")
)
write_csv(mc_curve, file.path(summary_directory, "mc_alpha_curve.csv"))
write_csv(
  pooled_null_bf,
  file.path(summary_directory, "pooled_null_bf_summary.csv")
)
write_csv(paired, file.path(summary_directory, "paired_method_comparison.csv"))
write_csv(
  rank_selection_summary,
  file.path(summary_directory, "paired_rank_selection_summary.csv")
)
write_csv(
  paired_mc_summary,
  file.path(summary_directory, "paired_mc_summary.csv")
)
write_csv(validation, file.path(summary_directory, "validation.csv"))

pilot_fit_path <- file.path(
  workflowr_root, "output", "revision_simulations", "internal",
  "r1_known_truth_genotype_permutation_seed12345_perm20260811_J200",
  "arm_fits.rds"
)
if (!file.exists(pilot_fit_path)) {
  stop("The saved pilot arm fits are required for the penalty-one check.")
}
pilot_fits <- readRDS(pilot_fit_path)
penalty_one_rows <- lapply(names(pilot_fits), function(arm) {
  raw_fit <- pilot_fits[[arm]]$raw_fit
  current_bf_fit <- pilot_fits[[arm]]$bf_fit
  penalty_one_fit <- fashr::fash_eb_est(
    L_matrix = raw_fit$L_matrix,
    penalty = 1,
    grid = raw_fit$psd_grid
  )
  full_weights <- numeric(length(raw_fit$psd_grid))
  full_weights[match(
    penalty_one_fit$prior_weight$psd, raw_fit$psd_grid
  )] <- penalty_one_fit$prior_weight$prior_weight
  penalty_one_alt <- full_weights[-1L] / sum(full_weights[-1L])
  current_alt <- fashr:::collapse_L(
    exp(raw_fit$L_matrix), log = FALSE
  )$pi_hat_star
  penalty_one_bf <- as.vector(
    exp(raw_fit$L_matrix[, -1L, drop = FALSE]) %*% penalty_one_alt /
      exp(raw_fit$L_matrix[, 1L])
  )
  penalty_one_bf_pi0 <- fashr:::BF_control(
    penalty_one_bf, plot = FALSE
  )$pi0_hat_star
  current_bf_pi0 <- current_bf_fit$prior_weights$prior_weight[
    current_bf_fit$prior_weights$psd == 0
  ]
  penalty_one_lfdr <- stats::plogis(
    stats::qlogis(penalty_one_bf_pi0) - log(penalty_one_bf)
  )
  data.frame(
    arm = arm,
    n_units = nrow(raw_fit$L_matrix),
    penalty_one_raw_pi0 = full_weights[1L],
    maximum_absolute_alternative_weight_difference = max(abs(
      penalty_one_alt - current_alt
    )),
    maximum_absolute_log_bf_difference = max(abs(
      log(penalty_one_bf) - log(current_bf_fit$BF)
    )),
    maximum_relative_bf_difference = max(
      abs(penalty_one_bf - current_bf_fit$BF) /
        pmax(abs(current_bf_fit$BF), .Machine$double.xmin)
    ),
    current_bf_pi0 = current_bf_pi0,
    penalty_one_corrected_bf_pi0 = penalty_one_bf_pi0,
    maximum_absolute_lfdr_difference = max(abs(
      penalty_one_lfdr - current_bf_fit$lfdr
    )),
    stringsAsFactors = FALSE
  )
})
penalty_one_equivalence <- do.call(rbind, penalty_one_rows)
penalty_one_equivalence$pass <- with(
  penalty_one_equivalence,
  maximum_absolute_alternative_weight_difference < 1e-6 &
    maximum_absolute_log_bf_difference < 1e-6 &
    abs(current_bf_pi0 - penalty_one_corrected_bf_pi0) < 1e-12 &
    maximum_absolute_lfdr_difference < 1e-8
)
if (!all(penalty_one_equivalence$pass)) {
  stop("The penalty-one equivalence sanity check failed.")
}
write_csv(
  penalty_one_equivalence,
  file.path(summary_directory, "penalty_one_equivalence_check.csv")
)

plot_curve <- mc_curve[mc_curve$alpha <= 0.10 + 1e-12, , drop = FALSE]
plot_curve$method_label <- factor(
  plot_curve$method,
  levels = c(
    "current_unpenalized_alternative",
    "corrected_penalized_alternative"
  ),
  labels = c(
    "Current BF update: re-estimated unpenalized alternative",
    "Corrected BF update: retained penalized alternative"
  )
)
arm_labels <- c(
  genuine_null_baseline = "Genuine null baseline",
  shared_genotype_permutation = "Genotype-permuted alternatives",
  signal_stripped_residual_permutation = "Signal-stripped residual permutation"
)
plot_curve$arm_label <- factor(
  arm_labels[plot_curve$arm], levels = unname(arm_labels)
)
alpha005_plot <- plot_curve[
  abs(plot_curve$alpha - 0.05) < 1e-12, , drop = FALSE
]

figure <- ggplot2::ggplot(
  plot_curve,
  ggplot2::aes(
    x = alpha, y = mean_fdr, color = method_label,
    linetype = method_label, group = method_label
  )
) +
  ggplot2::geom_abline(
    intercept = 0, slope = 1, color = "#777777", linetype = "dotted"
  ) +
  ggplot2::geom_line(linewidth = 1.0) +
  ggplot2::geom_point(data = alpha005_plot, size = 2.2) +
  ggplot2::facet_wrap(~arm_label, nrow = 1L) +
  ggplot2::scale_color_manual(values = c("#D55E00", "#0072B2")) +
  ggplot2::scale_linetype_manual(values = c("longdash", "solid")) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 0.10, by = 0.025),
    labels = scales::label_percent(accuracy = 0.1)
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0.01, 0.05))
  ) +
  ggplot2::labs(
    title = "Retaining the penalized alternative has little effect on R1 calibration",
    subtitle = paste(
      "Paired comparison on identical likelihood rows across five source seeds;",
      "points mark nominal alpha = 5%."
    ),
    x = "Nominal cumulative-FDR level",
    y = "Mean empirical FDR across five seeds",
    color = NULL,
    linetype = NULL,
    caption = paste(
      "The corrected update normalizes alternative weights from the raw",
      "penalty-10 fit before applying the unchanged BF-control rule."
    )
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 13),
    plot.subtitle = ggplot2::element_text(size = 10, color = "#4D4D4D"),
    strip.background = ggplot2::element_rect(fill = "#F2F2F2", color = NA),
    strip.text = ggplot2::element_text(face = "bold", size = 10),
    legend.position = "top",
    legend.text = ggplot2::element_text(size = 9),
    panel.grid.major.y = ggplot2::element_line(
      color = "#E5E5E5", linewidth = 0.35
    ),
    panel.spacing.x = grid::unit(1.0, "lines"),
    plot.caption = ggplot2::element_text(
      hjust = 0, size = 8.5, color = "#5A5A5A"
    )
  )

ggplot2::ggsave(
  file.path(figure_directory, "current_vs_penalized_alternative_fdr.png"),
  figure, width = 13.0, height = 5.4, dpi = 300, bg = "white"
)
ggplot2::ggsave(
  file.path(figure_directory, "current_vs_penalized_alternative_fdr.pdf"),
  figure, width = 13.0, height = 5.4, device = grDevices::cairo_pdf
)

saveRDS(
  list(
    input_specs = input_specs,
    unit_results = unit_results,
    replicate_summary = replicate_summary,
    alpha_curves = alpha_curves,
    mc_curve = mc_curve,
    pooled_null_bf = pooled_null_bf,
    paired = paired,
    rank_selection_summary = rank_selection_summary,
    paired_mc_summary = paired_mc_summary,
    validation = validation,
    penalty_one_equivalence = penalty_one_equivalence
  ),
  file.path(output_directory, "comparison_results.rds")
)

message("Wrote penalized-alternative BF comparison to: ", output_directory)
