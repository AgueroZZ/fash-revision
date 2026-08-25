#!/usr/bin/env Rscript

# Compare the completed R3 relative-clearance run with its frozen predecessor.

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
    return(normalizePath("coderepo-local", winslash = "/", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.", call. = FALSE)
}

require_complete_result <- function(path, label) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  required <- c(
    "complete.flag",
    "configuration.rds",
    "manifest.rds",
    file.path("summary", "functional_testing_mc_alpha_curve.csv"),
    file.path("summary", "all_replicate_functional_alpha_curves.csv"),
    file.path("summary", "all_functional_calls_alpha005.csv")
  )
  missing <- required[!file.exists(file.path(path, required))]
  if (length(missing) > 0L) {
    stop(label, " is missing: ", paste(missing, collapse = ", "))
  }
  path
}

read_curves <- function(result_dir, design) {
  replicate_curves <- utils::read.csv(file.path(
    result_dir, "summary", "all_replicate_functional_alpha_curves.csv"
  ))
  mc_curves <- utils::read.csv(file.path(
    result_dir, "summary", "functional_testing_mc_alpha_curve.csv"
  ))
  required_replicate <- c(
    "scenario", "target", "method", "alpha", "dynamic_discoveries",
    "n_discoveries", "false_discoveries", "conditional_discoveries",
    "conditional_false_discoveries", "first_stage_null_calls",
    "true_positives", "estimated_fsr", "empirical_fsr", "power", "seed",
    "truth_mechanism"
  )
  required_mc <- c(
    "scenario", "target", "method", "alpha", "mean_discoveries",
    "mean_false_discoveries", "mean_first_stage_null_calls",
    "mean_true_positives", "mean_power", "mean_estimated_fsr",
    "mean_empirical_fsr", "empirical_fsr_ci_lower",
    "empirical_fsr_ci_upper"
  )
  if (!all(required_replicate %in% names(replicate_curves)) ||
      !all(required_mc %in% names(mc_curves))) {
    stop("An R3 summary file has an unexpected schema.")
  }
  replicate_curves$design <- design
  mc_curves$design <- design
  mc_curves$truth_mechanism <- ifelse(
    grepl("^r3a_", mc_curves$scenario),
    "random_bspline",
    ifelse(grepl("^r3b_", mc_curves$scenario), "raised_cosine", NA_character_)
  )
  if (anyNA(mc_curves$truth_mechanism)) {
    stop("Could not identify truth mechanisms from scenario names.")
  }
  list(replicate = replicate_curves, mc = mc_curves)
}

make_shell_rows <- function(curves, cutoffs = c(0.05, 0.10, 0.15, 0.20)) {
  group_columns <- c(
    "design", "truth_mechanism", "target", "method", "seed"
  )
  groups <- split(
    curves,
    interaction(curves[group_columns], drop = TRUE, lex.order = TRUE)
  )
  rows <- lapply(groups, function(x) {
    x <- x[order(x$alpha), , drop = FALSE]
    if (!all(vapply(
      cutoffs,
      function(cutoff) sum(abs(x$alpha - cutoff) < 1e-12) == 1L,
      logical(1)
    ))) {
      stop("A replicate curve is missing a shell cutoff.")
    }
    selected <- x[match(cutoffs, x$alpha), , drop = FALSE]
    prepend_zero <- function(values) c(0, values)
    estimated_false_mass <- selected$estimated_fsr * selected$n_discoveries
    data.frame(
      design = selected$design,
      truth_mechanism = selected$truth_mechanism,
      target = selected$target,
      method = selected$method,
      seed = selected$seed,
      alpha_lower = c(0, head(cutoffs, -1L)),
      alpha_upper = cutoffs,
      shell_calls = diff(prepend_zero(selected$n_discoveries)),
      shell_false = diff(prepend_zero(selected$false_discoveries)),
      shell_conditional_false = diff(prepend_zero(
        selected$conditional_false_discoveries
      )),
      shell_first_stage_null = diff(prepend_zero(
        selected$first_stage_null_calls
      )),
      shell_true = diff(prepend_zero(selected$true_positives)),
      shell_estimated_false_mass = diff(prepend_zero(estimated_false_mass)),
      cumulative_calls = selected$n_discoveries,
      cumulative_false = selected$false_discoveries,
      cumulative_true = selected$true_positives,
      cumulative_estimated_fsr = selected$estimated_fsr,
      cumulative_empirical_fsr = selected$empirical_fsr,
      cumulative_power = selected$power,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  if (any(out$shell_calls < 0) ||
      any(out$shell_false < 0) ||
      any(out$shell_true < 0) ||
      any(abs(
        out$shell_false -
          out$shell_conditional_false -
          out$shell_first_stage_null
      ) > 1e-8)) {
    stop("The shell decomposition is not nested or internally consistent.")
  }
  out
}

summarize_shells <- function(shell_rows) {
  group_columns <- c(
    "design", "truth_mechanism", "target", "method",
    "alpha_lower", "alpha_upper"
  )
  groups <- split(
    shell_rows,
    interaction(shell_rows[group_columns], drop = TRUE, lex.order = TRUE)
  )
  rows <- lapply(groups, function(x) {
    calls <- sum(x$shell_calls)
    false_calls <- sum(x$shell_false)
    estimated_mass <- sum(x$shell_estimated_false_mass)
    data.frame(
      design = x$design[[1L]],
      truth_mechanism = x$truth_mechanism[[1L]],
      target = x$target[[1L]],
      method = x$method[[1L]],
      alpha_lower = x$alpha_lower[[1L]],
      alpha_upper = x$alpha_upper[[1L]],
      n_replications = nrow(x),
      pooled_shell_calls = calls,
      mean_shell_calls = mean(x$shell_calls),
      pooled_shell_false = false_calls,
      mean_shell_false = mean(x$shell_false),
      pooled_shell_true = sum(x$shell_true),
      mean_shell_true = mean(x$shell_true),
      pooled_shell_conditional_false = sum(x$shell_conditional_false),
      pooled_shell_first_stage_null = sum(x$shell_first_stage_null),
      actual_shell_fsr = if (calls == 0) NA_real_ else false_calls / calls,
      estimated_shell_fsr = if (calls == 0) NA_real_ else {
        estimated_mass / calls
      },
      mean_cumulative_power = mean(x$cumulative_power),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[do.call(order, unname(out[group_columns])), ]
}

summarize_alpha005_false_calls <- function(result_dir, design) {
  calls <- utils::read.csv(file.path(
    result_dir, "summary", "all_functional_calls_alpha005.csv"
  ))
  calls <- calls[
    calls$method == "FASH-IWP1-BF" &
      calls$target == "middle" &
      calls$false_discovery,
    ,
    drop = FALSE
  ]
  groups <- split(
    calls,
    interaction(calls$truth_mechanism, calls$truth_group, drop = TRUE)
  )
  rows <- lapply(groups, function(x) {
    data.frame(
      design = design,
      truth_mechanism = x$truth_mechanism[[1L]],
      truth_group = x$truth_group[[1L]],
      false_calls = nrow(x),
      true_functional_min = min(x$true_functional),
      true_functional_median = stats::median(x$true_functional),
      true_functional_max = max(x$true_functional),
      lfsr_min = min(x$lfsr),
      lfsr_median = stats::median(x$lfsr),
      lfsr_max = max(x$lfsr),
      stringsAsFactors = FALSE
    )
  })
  if (length(rows) == 0L) return(data.frame())
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$truth_mechanism, -out$false_calls), ]
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The diagnostic requires the existing ggplot2 installation.")
}

workflowr_root <- find_workflowr_root()
previous_result <- require_complete_result(
  get_arg("--previous-result-dir"),
  "Previous R3 result"
)
new_result <- require_complete_result(
  get_arg("--new-result-dir"),
  "New R3 result"
)
output_dir <- get_arg(
  "--output-dir",
  file.path(
    workflowr_root, "output", "revision_simulations", "diagnostics",
    "r3_relative_location_clearance_formal"
  )
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

design_levels <- c(
  "Previous fixed margin",
  "Relative margin + paired seed"
)
previous <- read_curves(previous_result, design_levels[[1L]])
new <- read_curves(new_result, design_levels[[2L]])
replicate_curves <- rbind(previous$replicate, new$replicate)
mc_curves <- rbind(previous$mc, new$mc)
replicate_curves$design <- factor(
  replicate_curves$design, levels = design_levels
)
mc_curves$design <- factor(mc_curves$design, levels = design_levels)

shell_rows <- make_shell_rows(replicate_curves)
shell_summary <- summarize_shells(shell_rows)
alpha005_false_calls <- rbind(
  summarize_alpha005_false_calls(previous_result, design_levels[[1L]]),
  summarize_alpha005_false_calls(new_result, design_levels[[2L]])
)

new_configuration <- readRDS(file.path(new_result, "configuration.rds"))
new_validation <- utils::read.csv(file.path(new_result, "scientific_validation.csv"))
if (!identical(
      new_configuration$functional_posterior_pairing,
      "common_random_seed_raw_bf"
    ) ||
    !isTRUE(all.equal(
      new_configuration$location_truth_min_range_fraction,
      0.10
    )) ||
    !all(!new_validation$passed)) {
  stop("The new result does not match the expected failed-gate design.")
}

utils::write.csv(
  mc_curves,
  file.path(output_dir, "old_vs_new_cumulative_curves.csv"),
  row.names = FALSE
)
utils::write.csv(
  shell_summary,
  file.path(output_dir, "old_vs_new_shell_diagnostics.csv"),
  row.names = FALSE
)
utils::write.csv(
  alpha005_false_calls,
  file.path(output_dir, "middle_false_calls_alpha005_by_truth_group.csv"),
  row.names = FALSE
)

mechanism_labels <- c(
  random_bspline = "R3A: broad random B-spline",
  raised_cosine = "R3B: compact raised cosine"
)
target_labels <- c(
  early = "Early", middle = "Middle", late = "Late", switch = "Switch"
)
design_colors <- c(
  "Previous fixed margin" = "#777777",
  "Relative margin + paired seed" = "#0072B2"
)

cumulative_plot_data <- mc_curves[
  mc_curves$method == "FASH-IWP1-BF",
  ,
  drop = FALSE
]
cumulative_plot_data$mechanism_label <- factor(
  mechanism_labels[cumulative_plot_data$truth_mechanism],
  levels = unname(mechanism_labels)
)
cumulative_plot_data$target_label <- factor(
  target_labels[cumulative_plot_data$target],
  levels = unname(target_labels)
)

cumulative_plot <- ggplot2::ggplot(
  cumulative_plot_data,
  ggplot2::aes(
    x = alpha,
    y = mean_empirical_fsr,
    color = design,
    linetype = design
  )
) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    color = "#555555",
    linewidth = 0.35,
    linetype = 3
  ) +
  ggplot2::geom_ribbon(
    ggplot2::aes(
      ymin = empirical_fsr_ci_lower,
      ymax = empirical_fsr_ci_upper,
      fill = design
    ),
    color = NA,
    alpha = 0.10
  ) +
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(mechanism_label),
    cols = ggplot2::vars(target_label)
  ) +
  ggplot2::scale_color_manual(values = design_colors, drop = FALSE) +
  ggplot2::scale_fill_manual(values = design_colors, drop = FALSE) +
  ggplot2::scale_linetype_manual(
    values = c("Previous fixed margin" = 2, "Relative margin + paired seed" = 1),
    drop = FALSE
  ) +
  ggplot2::coord_cartesian(xlim = c(0, 0.20), ylim = c(0, 0.30)) +
  ggplot2::labs(
    x = "Nominal alpha",
    y = "Mean empirical FSR",
    color = NULL,
    fill = NULL,
    linetype = NULL,
    title = "Relative truth clearance does not remove Middle tail miscalibration",
    subtitle = "BF-updated FASH; bands are five-seed 95% intervals"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold")
  )
ggplot2::ggsave(
  file.path(output_dir, "old_vs_new_cumulative_fsr.png"),
  cumulative_plot,
  width = 13.5,
  height = 6.8,
  dpi = 180
)

middle_shells <- shell_summary[
  shell_summary$method == "FASH-IWP1-BF" &
    shell_summary$target == "middle",
  ,
  drop = FALSE
]
middle_shells$design <- factor(middle_shells$design, levels = design_levels)
middle_shells$mechanism_label <- factor(
  mechanism_labels[middle_shells$truth_mechanism],
  levels = unname(mechanism_labels)
)
middle_shells$shell <- factor(
  sprintf("(%0.2f, %0.2f]", middle_shells$alpha_lower, middle_shells$alpha_upper),
  levels = sprintf(
    "(%0.2f, %0.2f]",
    c(0, 0.05, 0.10, 0.15),
    c(0.05, 0.10, 0.15, 0.20)
  )
)
shell_long <- rbind(
  data.frame(
    middle_shells[c("design", "mechanism_label", "shell")],
    estimand = "Actual false fraction",
    value = middle_shells$actual_shell_fsr
  ),
  data.frame(
    middle_shells[c("design", "mechanism_label", "shell")],
    estimand = "Mean posterior LFSR",
    value = middle_shells$estimated_shell_fsr
  )
)
shell_long$estimand <- factor(
  shell_long$estimand,
  levels = c("Actual false fraction", "Mean posterior LFSR")
)

shell_plot <- ggplot2::ggplot(
  shell_long,
  ggplot2::aes(
    x = shell,
    y = value,
    color = estimand,
    group = estimand,
    shape = estimand
  )
) +
  ggplot2::geom_line(linewidth = 0.75) +
  ggplot2::geom_point(size = 2.0) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(mechanism_label),
    cols = ggplot2::vars(design)
  ) +
  ggplot2::scale_color_manual(values = c(
    "Actual false fraction" = "#D55E00",
    "Mean posterior LFSR" = "#0072B2"
  )) +
  ggplot2::scale_shape_manual(values = c(
    "Actual false fraction" = 16,
    "Mean posterior LFSR" = 1
  )) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    x = "Newly admitted nominal-alpha shell",
    y = "Pooled fraction across five seeds",
    color = NULL,
    shape = NULL,
    title = "R3A Middle breaks down in newly admitted high-alpha shells",
    subtitle = paste(
      "Actual shell errors approach one after Middle power saturates;",
      "posterior LFSR remains materially lower"
    )
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    strip.text = ggplot2::element_text(face = "bold")
  )
ggplot2::ggsave(
  file.path(output_dir, "middle_shell_actual_vs_estimated_fsr.png"),
  shell_plot,
  width = 11.5,
  height = 7.0,
  dpi = 180
)

new_middle <- cumulative_plot_data[
  cumulative_plot_data$design == design_levels[[2L]] &
    cumulative_plot_data$target == "middle",
  ,
  drop = FALSE
]
count_long <- rbind(
  data.frame(
    new_middle[c("alpha", "truth_mechanism")],
    call_type = "True Middle calls",
    mean_calls = new_middle$mean_true_positives
  ),
  data.frame(
    new_middle[c("alpha", "truth_mechanism")],
    call_type = "False Middle calls",
    mean_calls = new_middle$mean_false_discoveries
  )
)
count_long$mechanism_label <- factor(
  mechanism_labels[count_long$truth_mechanism],
  levels = unname(mechanism_labels)
)
count_long$call_type <- factor(
  count_long$call_type,
  levels = c("True Middle calls", "False Middle calls")
)

count_plot <- ggplot2::ggplot(
  count_long,
  ggplot2::aes(
    x = alpha,
    y = mean_calls,
    color = call_type,
    linetype = call_type
  )
) +
  ggplot2::geom_line(linewidth = 0.85) +
  ggplot2::facet_wrap(ggplot2::vars(mechanism_label), nrow = 1) +
  ggplot2::scale_color_manual(values = c(
    "True Middle calls" = "#0072B2",
    "False Middle calls" = "#D55E00"
  )) +
  ggplot2::scale_linetype_manual(values = c(
    "True Middle calls" = 1,
    "False Middle calls" = 2
  )) +
  ggplot2::labs(
    x = "Nominal alpha",
    y = "Mean cumulative calls across five seeds",
    color = NULL,
    linetype = NULL,
    title = "R3A exhausts true Middle discoveries before alpha reaches 0.20",
    subtitle = "The remaining high-alpha additions are predominantly wrong-category curves"
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold")
  )
ggplot2::ggsave(
  file.path(output_dir, "new_middle_true_and_false_call_saturation.png"),
  count_plot,
  width = 10.5,
  height = 4.8,
  dpi = 180
)

writeLines(
  c(
    paste0("generated_at=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("previous_result=", previous_result),
    paste0("new_result=", new_result),
    "comparison_scope=FASH-IWP1-BF cumulative and marginal-shell calibration",
    "shell_cutoffs=0.05,0.10,0.15,0.20",
    "aggregation=pooled counts across five fixed seeds",
    "new_formal_gate_passed=false"
  ),
  file.path(output_dir, "provenance.txt")
)

message("Saved formal R3 calibration diagnostics to: ", output_dir)
