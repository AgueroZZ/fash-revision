#!/usr/bin/env Rscript

# Generate and audit visually unambiguous local cubic B-spline truth curves
# before connecting this alternative truth mechanism to the full FASH workflow.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) return(".")
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find code/revision_simulations/shared/simulation_functions.R from the current working directory.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0 || hit[1] == length(args)) return(default)
  args[hit[1] + 1]
}

write_csv <- function(x, path) write.csv(x, file = path, row.names = FALSE)

format_group_name <- function(group) gsub("_", " / ", group, fixed = TRUE)

make_group_summary <- function(diagnostics) {
  split_rows <- split(diagnostics, diagnostics$truth_group)
  summary_rows <- lapply(split_rows, function(dat) {
    target_name <- unique(dat$time_group)
    data.frame(
      truth_group = unique(dat$truth_group),
      n_examples = nrow(dat),
      target_peak_min = min(dat$target_peak),
      outside_peak_max = max(dat$outside_peak),
      target_to_outside_ratio_min = min(dat$target_to_outside_ratio),
      target_functional_min = min(dat[[paste0(target_name, "_functional")]]),
      switch_functional_min = min(dat$switch_functional),
      sign_transition_min = min(dat$effective_sign_transitions),
      sign_transition_max = max(dat$effective_sign_transitions),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, summary_rows)
}

plot_targeted_truth_examples <- function(curves,
                                         time_grid,
                                         output_path,
                                         device = c("png", "pdf")) {
  device <- match.arg(device)
  group_order <- c(
    "early_switch", "early_non-switch",
    "middle_switch", "middle_non-switch",
    "late_switch", "late_non-switch"
  )
  n_examples <- max(vapply(curves, length, integer(1)))
  if (device == "png") {
    grDevices::png(output_path, width = 2800, height = 3900, res = 220)
  } else {
    grDevices::pdf(output_path, width = 13, height = 18, onefile = TRUE)
  }
  on.exit(grDevices::dev.off(), add = TRUE)

  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(
    mfrow = c(length(group_order), n_examples),
    mar = c(3.6, 3.7, 4.3, 0.6),
    oma = c(2.2, 0, 2.4, 0)
  )
  window_colors <- c(early = "#d7eaf7", middle = "#dff0d8", late = "#fde6c9")

  for (group in group_order) {
    examples <- curves[[group]]
    for (example_index in seq_len(n_examples)) {
      example <- examples[[example_index]]
      group_color <- window_colors[[example$time_group]]
      y_limit <- max(2.6, 1.10 * max(abs(example$beta_evaluation)))
      plot(
        example$evaluation_grid,
        example$beta_evaluation,
        type = "n",
        xlab = "Time",
        ylab = "True genetic effect",
        ylim = c(-y_limit, y_limit),
        main = ""
      )
      if (example$time_group == "early") {
        rect(0, -y_limit, 3, y_limit, col = group_color, border = NA)
      } else if (example$time_group == "middle") {
        rect(4, -y_limit, 11, y_limit, col = group_color, border = NA)
      } else {
        rect(12, -y_limit, 15, y_limit, col = group_color, border = NA)
      }
      abline(h = 0, lty = 3, col = "grey55")
      lines(
        example$evaluation_grid,
        example$beta_evaluation,
        col = "#d1495b",
        lwd = 2.5
      )
      points(
        time_grid,
        example$beta_observed,
        pch = 16,
        cex = 0.55,
        col = "#1f2933"
      )
      mtext(
        format_group_name(group),
        side = 3,
        line = 2.25,
        cex = 0.92,
        font = 2
      )
      mtext(
        paste0(
          "example ", example_index,
          "; peak ratio = ", formatC(example$diagnostics$target_to_outside_ratio, format = "f", digits = 2),
          "; sign transitions = ", example$diagnostics$effective_sign_transitions
        ),
        side = 3,
        line = 0.85,
        cex = 0.64
      )
    }
  }
  mtext(
    "Targeted local cubic B-spline truth pilot",
    side = 3,
    outer = TRUE,
    cex = 1.30,
    font = 2
  )
  mtext(
    "Shaded region is the early, middle, or late target window. Points are the 16 observed time-grid values.",
    side = 1,
    outer = TRUE,
    line = 0.5,
    cex = 0.88
  )
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

seed <- as.integer(get_arg("--seed", "20260724"))
n_examples_per_group <- as.integer(get_arg("--n-examples-per-group", "3"))
amplitude <- as.numeric(get_arg("--amplitude", "2"))
switch_threshold <- as.numeric(get_arg("--switch-threshold", "0.25"))
minimum_location_margin <- as.numeric(get_arg("--minimum-location-margin", "0.60"))
minimum_location_ratio <- as.numeric(get_arg("--minimum-location-ratio", "2"))
targeted_profile <- get_arg("--targeted-profile", "narrow")
non_switch_baseline_fraction <- as.numeric(get_arg("--non-switch-baseline-fraction", "0.15"))
non_switch_background_fraction <- as.numeric(get_arg("--non-switch-background-fraction", "0.05"))
output_id <- get_arg("--output-id", "targeted_local_bspline_truth_v1")

if (!is.finite(seed) || n_examples_per_group < 1 || !is.finite(amplitude) || amplitude <= 0 ||
    !is.finite(switch_threshold) || switch_threshold <= 0 ||
    !is.finite(minimum_location_margin) || minimum_location_margin <= 0 ||
    !is.finite(minimum_location_ratio) || minimum_location_ratio <= 1 ||
    !targeted_profile %in% c("narrow", "broad") ||
    !is.finite(non_switch_baseline_fraction) || non_switch_baseline_fraction <= 0 ||
    !is.finite(non_switch_background_fraction) || non_switch_background_fraction < 0 ||
    !nzchar(output_id)) {
  stop("Invalid pilot settings.")
}

time_grid <- make_time_grid()
evaluation_grid <- seq(min(time_grid), max(time_grid), by = 0.1)
group_table <- expand.grid(
  time_group = c("early", "middle", "late"),
  switch_status = c("switch", "non-switch"),
  stringsAsFactors = FALSE
)
group_table$truth_group <- paste(group_table$time_group, group_table$switch_status, sep = "_")
group_table <- group_table[match(
  c("early_switch", "early_non-switch", "middle_switch", "middle_non-switch", "late_switch", "late_non-switch"),
  group_table$truth_group
), , drop = FALSE]

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "pilot",
  output_id
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(seed)
curves <- stats::setNames(vector("list", nrow(group_table)), group_table$truth_group)
diagnostic_rows <- vector("list", nrow(group_table) * n_examples_per_group)
counter <- 0L
for (group_index in seq_len(nrow(group_table))) {
  group <- group_table[group_index, , drop = FALSE]
  group_curves <- vector("list", n_examples_per_group)
  for (example_index in seq_len(n_examples_per_group)) {
    truth <- sample_targeted_local_bspline_truth(
      time_group = group$time_group,
      switch_status = group$switch_status,
      time_grid = time_grid,
      evaluation_grid = evaluation_grid,
      amplitude = amplitude,
      switch_threshold = switch_threshold,
      minimum_location_margin = minimum_location_margin,
      minimum_location_ratio = minimum_location_ratio,
      non_switch_baseline_fraction = non_switch_baseline_fraction,
      non_switch_background_fraction = non_switch_background_fraction,
      profile = targeted_profile
    )
    group_curves[[example_index]] <- c(
      truth,
      list(
        time_group = group$time_group,
        switch_status = group$switch_status,
        truth_group = group$truth_group,
        evaluation_grid = evaluation_grid
      )
    )
    counter <- counter + 1L
    diagnostic_rows[[counter]] <- cbind(
      data.frame(
        truth_group = group$truth_group,
        example_id = example_index,
        attempt = truth$attempt,
        stringsAsFactors = FALSE
      ),
      truth$diagnostics,
      as.data.frame(as.list(stats::setNames(
        truth$contrasts,
        paste0(names(truth$contrasts), "_functional")
      )), check.names = FALSE)
    )
  }
  curves[[group$truth_group]] <- group_curves
}
diagnostics <- do.call(rbind, diagnostic_rows)
row.names(diagnostics) <- NULL

target_contrast <- vapply(
  seq_len(nrow(diagnostics)),
  function(index) diagnostics[[paste0(diagnostics$time_group[index], "_functional")]][index],
  numeric(1)
)
switch_rows <- diagnostics$switch_status == "switch"
non_switch_rows <- !switch_rows
constraint_failure <- diagnostics$target_to_outside_ratio < minimum_location_ratio |
  target_contrast < minimum_location_margin |
  (switch_rows & diagnostics$switch_functional <= 0) |
  (switch_rows & diagnostics$effective_sign_transitions != 1L) |
  (non_switch_rows & diagnostics$switch_functional > 0) |
  (non_switch_rows & diagnostics$effective_sign_transitions != 0L)
if (any(constraint_failure)) {
  print(diagnostics[constraint_failure, , drop = FALSE])
  stop("Targeted local B-spline pilot failed its truth constraints.")
}

configuration <- list(
  seed = seed,
  n_examples_per_group = n_examples_per_group,
  amplitude = amplitude,
  switch_threshold = switch_threshold,
  minimum_location_margin = minimum_location_margin,
  minimum_location_ratio = minimum_location_ratio,
  targeted_profile = targeted_profile,
  non_switch_baseline_fraction = non_switch_baseline_fraction,
  non_switch_background_fraction = non_switch_background_fraction,
  time_grid = time_grid,
  evaluation_grid = evaluation_grid
)
saveRDS(configuration, file.path(output_dir, "configuration.rds"))
saveRDS(curves, file.path(output_dir, "targeted_truth_curves.rds"))
write_csv(diagnostics, file.path(output_dir, "targeted_truth_diagnostics.csv"))
write_csv(make_group_summary(diagnostics), file.path(output_dir, "targeted_truth_group_summary.csv"))

plot_targeted_truth_examples(
  curves = curves,
  time_grid = time_grid,
  output_path = file.path(output_dir, "targeted_local_bspline_truth_examples.png"),
  device = "png"
)
plot_targeted_truth_examples(
  curves = curves,
  time_grid = time_grid,
  output_path = file.path(output_dir, "targeted_local_bspline_truth_examples.pdf"),
  device = "pdf"
)

message("Saved targeted local B-spline truth pilot to: ", output_dir)
message("Minimum target-to-outside peak ratio: ", formatC(
  min(diagnostics$target_to_outside_ratio), format = "f", digits = 2
))
