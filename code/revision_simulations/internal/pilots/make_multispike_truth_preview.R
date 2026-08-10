#!/usr/bin/env Rscript

# Preview a versioned spiky B-spline truth mechanism before rerunning the
# genotype-level simulation and downstream methods.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) return(".")
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find code/revision_simulations/shared/simulation_functions.R.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0 || hit[1] == length(args)) return(default)
  args[hit[1] + 1]
}

target_window_limits <- function(time_group) {
  switch(
    time_group,
    early = c(0, 3),
    middle = c(4, 11),
    late = c(12, 15),
    stop("Unsupported time group.")
  )
}

plot_multispike_preview <- function(examples, output_path, device = c("png", "pdf")) {
  device <- match.arg(device)
  time_groups <- c("early", "middle", "late")
  pattern_order <- c("single", "same-sign double", "opposite-sign double")
  examples_per_pattern <- 2L
  if (device == "png") {
    grDevices::png(output_path, width = 3600, height = 2100, res = 220)
  } else {
    grDevices::pdf(output_path, width = 18, height = 10.5, onefile = TRUE)
  }
  on.exit(grDevices::dev.off(), add = TRUE)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(
    mfrow = c(length(time_groups), length(pattern_order) * examples_per_pattern),
    mar = c(3.1, 3.5, 3.6, 0.5),
    oma = c(2.2, 2.0, 3.2, 0.5)
  )
  group_colors <- c(early = "#d8edf3", middle = "#e6efd8", late = "#f7e1c5")
  curve_colors <- c(
    "single" = "#007f73",
    "same-sign double" = "#007f73",
    "opposite-sign double" = "#b23a48"
  )

  for (time_group in time_groups) {
    for (pattern in pattern_order) {
      subset_examples <- examples[
        vapply(examples, function(x) {
          x$time_group == time_group && x$spike_pattern == pattern
        }, logical(1))
      ]
      if (length(subset_examples) != examples_per_pattern) {
        stop("Each time-group/pattern cell must contain exactly two examples.")
      }
      for (example_index in seq_along(subset_examples)) {
        example <- subset_examples[[example_index]]
        y_limit <- 1.08 * max(abs(example$beta_evaluation))
        plot(
          example$evaluation_grid,
          example$beta_evaluation,
          type = "n",
          xlab = "Time",
          ylab = "True genetic effect",
          ylim = c(-y_limit, y_limit),
          main = ""
        )
        limits <- target_window_limits(time_group)
        rect(
          limits[1],
          -y_limit,
          limits[2],
          y_limit,
          col = group_colors[[time_group]],
          border = NA
        )
        abline(h = 0, col = "grey55", lty = 3)
        lines(
          example$evaluation_grid,
          example$beta_evaluation,
          col = curve_colors[[pattern]],
          lwd = 2.6
        )
        points(
          example$time_grid,
          example$beta_observed,
          pch = 16,
          cex = 0.52,
          col = "#1f2933"
        )
        mtext(
          paste(time_group, pattern, sep = " / "),
          side = 3,
          line = 1.9,
          cex = 0.78,
          font = 2
        )
        mtext(
          paste0(
            "example ", example_index,
            "; ratio ", formatC(
              example$diagnostics$target_to_outside_ratio,
              format = "f",
              digits = 2
            ),
            "; transitions ", example$diagnostics$effective_sign_transitions
          ),
          side = 3,
          line = 0.55,
          cex = 0.57
        )
      }
    }
  }
  mtext(
    "Multi-spike local cubic B-spline truth preview",
    side = 3,
    outer = TRUE,
    line = 1.4,
    cex = 1.35,
    font = 2
  )
  mtext(
    "Non-switch truths mix one and two same-sign spikes; switch truths use two opposite-sign spikes.",
    side = 3,
    outer = TRUE,
    line = 0.1,
    cex = 0.88
  )
  mtext(
    "Shading marks the assigned target window; points mark the 16 observed time values.",
    side = 1,
    outer = TRUE,
    line = 0.55,
    cex = 0.82
  )
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

seed <- as.integer(get_arg("--seed", "20260726"))
output_id <- get_arg("--output-id", "multispike_truth_v2")
non_switch_baseline_fraction <- as.numeric(get_arg(
  "--non-switch-baseline-fraction",
  "0"
))
if (!is.finite(seed) || !nzchar(output_id) ||
    !is.finite(non_switch_baseline_fraction) ||
    non_switch_baseline_fraction < 0) {
  stop("Invalid preview settings.")
}

time_grid <- make_time_grid()
evaluation_grid <- seq(min(time_grid), max(time_grid), by = 0.05)
time_groups <- c("early", "middle", "late")
preview_design <- expand.grid(
  example_index = 1:2,
  spike_pattern = c("single", "same-sign double", "opposite-sign double"),
  time_group = time_groups,
  stringsAsFactors = FALSE
)
preview_design$switch_status <- ifelse(
  preview_design$spike_pattern == "opposite-sign double",
  "switch",
  "non-switch"
)
preview_design$spike_count <- ifelse(
  preview_design$spike_pattern == "single",
  1L,
  2L
)

set.seed(seed)
examples <- vector("list", nrow(preview_design))
diagnostic_rows <- vector("list", nrow(preview_design))
for (row_index in seq_len(nrow(preview_design))) {
  design_row <- preview_design[row_index, , drop = FALSE]
  truth <- sample_multispike_local_bspline_truth(
    time_group = design_row$time_group,
    switch_status = design_row$switch_status,
    spike_count = design_row$spike_count,
    time_grid = time_grid,
    evaluation_grid = evaluation_grid,
    df = 16,
    degree = 3,
    secondary_fraction = c(0.40, 0.65),
    switch_threshold = 0.25,
    minimum_location_margin = 0.60,
    minimum_location_ratio = 1.4,
    target_centered_rms = 0.90,
    non_switch_baseline_fraction = non_switch_baseline_fraction,
    direction = if (design_row$example_index == 1L) 1 else -1
  )
  examples[[row_index]] <- c(
    truth,
    list(
      time_group = design_row$time_group,
      switch_status = design_row$switch_status,
      spike_pattern = design_row$spike_pattern,
      example_index = design_row$example_index,
      time_grid = time_grid,
      evaluation_grid = evaluation_grid
    )
  )
  diagnostic_rows[[row_index]] <- cbind(
    example_index = design_row$example_index,
    truth$diagnostics,
    early_functional = truth$contrasts[["early"]],
    middle_functional = truth$contrasts[["middle"]],
    late_functional = truth$contrasts[["late"]],
    switch_functional = truth$contrasts[["switch"]],
    generation_attempt = truth$attempt,
    stringsAsFactors = FALSE
  )
}
diagnostics <- do.call(rbind, diagnostic_rows)
row.names(diagnostics) <- NULL

if (any(abs(diagnostics$centered_rms - 0.90) > 1e-8) ||
    any(diagnostics$target_to_outside_ratio < 1.4) ||
    any(diagnostics$effective_sign_transitions[
      diagnostics$switch_status == "non-switch"
    ] != 0L) ||
    any(diagnostics$effective_sign_transitions[
      diagnostics$switch_status == "switch"
    ] < 1L) ||
    any(diagnostics$switch_functional[
      diagnostics$switch_status == "non-switch"
    ] > 0) ||
    any(diagnostics$switch_functional[
      diagnostics$switch_status == "switch"
    ] <= 0)) {
  stop("The generated preview failed its truth-shape validation.")
}

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "previews",
  output_id
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    examples = examples,
    diagnostics = diagnostics,
    settings = list(
      truth_mechanism = "mixed_single_double_v2",
      seed = seed,
      df = 16,
      degree = 3,
      target_centered_rms = 0.90,
      minimum_location_ratio = 1.4,
      secondary_fraction = c(0.40, 0.65),
      non_switch_baseline_fraction = non_switch_baseline_fraction
    )
  ),
  file.path(output_dir, "multispike_truth_preview.rds")
)
write.csv(
  diagnostics,
  file.path(output_dir, "multispike_truth_diagnostics.csv"),
  row.names = FALSE
)
plot_multispike_preview(
  examples,
  file.path(output_dir, "multispike_truth_preview.png"),
  device = "png"
)
plot_multispike_preview(
  examples,
  file.path(output_dir, "multispike_truth_preview.pdf"),
  device = "pdf"
)

cat("Saved multi-spike truth preview to:", normalizePath(output_dir), "\n")
cat(
  "Pattern counts:\n",
  paste(capture.output(print(table(diagnostics$spike_pattern))), collapse = "\n"),
  "\n"
)
cat(
  "Minimum target-to-outside ratio:",
  formatC(min(diagnostics$target_to_outside_ratio), format = "f", digits = 3),
  "\n"
)
