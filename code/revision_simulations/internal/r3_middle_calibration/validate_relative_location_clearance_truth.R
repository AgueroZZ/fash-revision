#!/usr/bin/env Rscript

# Validate and visualize truth under scale-adaptive location clearance.

options(stringsAsFactors = FALSE)

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[[1L]] == length(args)) return(default)
  args[[hit[[1L]] + 1L]]
}

parse_integer_list <- function(x) {
  values <- suppressWarnings(as.integer(
    trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  ))
  if (length(values) < 1L || anyNA(values) || anyDuplicated(values)) {
    stop("--seed-list must contain unique comma-separated integers.")
  }
  values
}

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }
  if (file.exists(file.path(
    "coderepo-local", "code", "revision_simulations", "shared",
    "simulation_functions.R"
  ))) {
    return(normalizePath(
      "coderepo-local", winslash = "/", mustWork = TRUE
    ))
  }
  stop("Could not find the workflowr repository root.", call. = FALSE)
}

summarize_distribution <- function(data) {
  groups <- split(
    data,
    list(data$truth_mechanism, data$target, data$truth_status),
    drop = TRUE
  )
  rows <- lapply(groups, function(x) {
    values <- x$functional_value
    quantiles <- stats::quantile(
      values,
      probs = c(0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99),
      names = FALSE
    )
    data.frame(
      truth_mechanism = x$truth_mechanism[[1L]],
      target = x$target[[1L]],
      truth_status = x$truth_status[[1L]],
      n = length(values),
      minimum = min(values),
      q01 = quantiles[[1L]],
      q05 = quantiles[[2L]],
      q25 = quantiles[[3L]],
      median = quantiles[[4L]],
      q75 = quantiles[[5L]],
      q95 = quantiles[[6L]],
      q99 = quantiles[[7L]],
      maximum = max(values),
      minimum_absolute_value = min(abs(values)),
      proportion_abs_below_0_15 = mean(abs(values) < 0.15),
      proportion_abs_below_0_25 = mean(abs(values) < 0.25),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$truth_mechanism, out$target, out$truth_status), ]
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The diagnostic requires the existing ggplot2 installation.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code/revision_simulations/shared/simulation_functions.R"
))

n_variants <- as.integer(get_arg("--n-variants", "6362"))
seed_list <- parse_integer_list(get_arg(
  "--seed-list",
  "12345,22345,32345,42345,52345"
))
location_truth_margin <- as.numeric(get_arg(
  "--location-truth-margin",
  "0.10"
))
location_truth_min_range_fraction <- as.numeric(get_arg(
  "--location-truth-min-range-fraction",
  "0.10"
))
if (n_variants < 30L ||
    !is.finite(location_truth_margin) ||
    location_truth_margin <= 0 ||
    !is.finite(location_truth_min_range_fraction) ||
    location_truth_min_range_fraction < 0 ||
    location_truth_min_range_fraction >= 1) {
  stop("Invalid truth-only diagnostic settings.")
}

output_dir <- get_arg(
  "--output-dir",
  file.path(
    workflowr_root, "output", "revision_simulations", "diagnostics",
    "r3_relative_location_clearance"
  )
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

time_grid <- make_time_grid()
evaluation_grid <- seq(0, 15, by = 0.1)
class_probs <- c(
  dynamic_bspline = 0.20,
  constant = 0.40,
  zero = 0.40
)
temporal_category_probs <- c(
  early = 0.29,
  middle = 0.42,
  late = 0.29
)
truth_mechanisms <- c("random_bspline", "raised_cosine")
targets <- c("early", "middle", "late", "switch")

functional_rows <- list()
clearance_rows <- list()
support_rows <- list()
functional_index <- 1L
clearance_index <- 1L
support_index <- 1L

for (truth_mechanism in truth_mechanisms) {
  for (seed in seed_list) {
    message("Generating ", truth_mechanism, " truth for seed ", seed, ".")
    component_seeds <- revision_component_seeds(seed)
    effect_sim <- simulate_matched_functional_effect_set(
      n_variants = n_variants,
      truth_mechanism = truth_mechanism,
      time_grid = time_grid,
      evaluation_grid = evaluation_grid,
      class_probs = class_probs,
      dynamic_main_effect_sd = 1,
      switch_threshold = 0.25,
      location_truth_margin = location_truth_margin,
      location_truth_min_range_fraction =
        location_truth_min_range_fraction,
      switch_truth_margin = 0.10,
      non_switch_min_abs = 0.10,
      non_switch_min_range_fraction = 0.10,
      temporal_category_probs = temporal_category_probs,
      seed = seed,
      class_seed = component_seeds[["classes"]],
      constant_seed = component_seeds[["constant_effects"]],
      shape_seed = component_seeds[["functional_truth"]],
      scenario = paste0(
        "r3_truth_only_relative_location_clearance_",
        truth_mechanism
      ),
      middle_window = c(3, 12),
      middle_boundary = "open"
    )

    dynamic_indices <- which(
      effect_sim$unit_info$effect_class == "dynamic_bspline"
    )
    dynamic_info <- effect_sim$unit_info[dynamic_indices, , drop = FALSE]
    dynamic_functionals <- effect_sim$true_functionals[
      dynamic_indices,
      targets,
      drop = FALSE
    ]

    if (truth_mechanism == "raised_cosine") {
      primary_centers <- vapply(
        dynamic_info$peak_centers,
        function(x) {
          if (length(x) < 1L) NA_real_ else x[[1L]]
        },
        numeric(1)
      )
      half_width <- effect_sim$settings$cosine_width_half
      support_left <- primary_centers - half_width
      support_right <- primary_centers + half_width
      target_intervals <- list(
        early = c(0, 3),
        middle = c(3, 12),
        late = c(12, 15)
      )
      target_left <- vapply(
        dynamic_info$time_group,
        function(x) target_intervals[[x]][[1L]],
        numeric(1)
      )
      target_right <- vapply(
        dynamic_info$time_group,
        function(x) target_intervals[[x]][[2L]],
        numeric(1)
      )
      support_rows[[support_index]] <- data.frame(
        seed = seed,
        curve_index = dynamic_indices,
        truth_group = dynamic_info$truth_group,
        assigned_time_group = dynamic_info$time_group,
        assigned_switch_status = dynamic_info$switch_status,
        primary_center = primary_centers,
        primary_support_left = support_left,
        primary_support_right = support_right,
        target_interval_left = target_left,
        target_interval_right = target_right,
        support_within_observed_domain =
          support_left >= min(time_grid) & support_right <= max(time_grid),
        support_within_target_interval =
          support_left >= target_left & support_right <= target_right,
        observed_domain_truncation =
          pmax(0, min(time_grid) - support_left) +
          pmax(0, support_right - max(time_grid)),
        stringsAsFactors = FALSE
      )
      support_index <- support_index + 1L
    }

    target_columns <- match(dynamic_info$time_group, targets)
    target_values <- dynamic_functionals[cbind(
      seq_along(dynamic_indices),
      target_columns
    )]
    largest_competing_values <- vapply(
      seq_along(dynamic_indices),
      function(i) {
        competing <- setdiff(
          c("early", "middle", "late"),
          dynamic_info$time_group[[i]]
        )
        max(dynamic_functionals[i, competing])
      },
      numeric(1)
    )
    required <- dynamic_info$required_location_clearance
    if (any(target_values < required) ||
        any(largest_competing_values > -required)) {
      stop("Generated truth violated relative location clearance.")
    }

    clearance_rows[[clearance_index]] <- data.frame(
      seed = seed,
      truth_mechanism = truth_mechanism,
      curve_index = dynamic_indices,
      truth_group = dynamic_info$truth_group,
      assigned_time_group = dynamic_info$time_group,
      assigned_switch_status = dynamic_info$switch_status,
      effect_range = dynamic_info$effect_range,
      required_location_clearance = required,
      assigned_target_functional = target_values,
      target_clearance_excess = target_values - required,
      largest_competing_functional = largest_competing_values,
      competing_clearance_excess = -required - largest_competing_values,
      stringsAsFactors = FALSE
    )
    clearance_index <- clearance_index + 1L

    for (target in targets) {
      values <- dynamic_functionals[, target]
      functional_rows[[functional_index]] <- data.frame(
        seed = seed,
        truth_mechanism = truth_mechanism,
        curve_index = dynamic_indices,
        truth_group = dynamic_info$truth_group,
        assigned_time_group = dynamic_info$time_group,
        assigned_switch_status = dynamic_info$switch_status,
        target = target,
        functional_value = values,
        truth_status = ifelse(
          values > 0,
          "Functional positive",
          "Functional null"
        ),
        stringsAsFactors = FALSE
      )
      functional_index <- functional_index + 1L
    }
  }
}

functional_data <- do.call(rbind, functional_rows)
clearance_data <- do.call(rbind, clearance_rows)
support_data <- do.call(rbind, support_rows)
rownames(functional_data) <- NULL
rownames(clearance_data) <- NULL
rownames(support_data) <- NULL

support_groups <- split(
  support_data,
  support_data$assigned_time_group,
  drop = TRUE
)
support_summary <- do.call(rbind, lapply(support_groups, function(x) {
  data.frame(
    assigned_time_group = x$assigned_time_group[[1L]],
    dynamic_truths = nrow(x),
    minimum_primary_center = min(x$primary_center),
    maximum_primary_center = max(x$primary_center),
    proportion_support_within_observed_domain =
      mean(x$support_within_observed_domain),
    proportion_support_within_target_interval =
      mean(x$support_within_target_interval),
    mean_observed_domain_truncation = mean(x$observed_domain_truncation),
    maximum_observed_domain_truncation = max(x$observed_domain_truncation),
    stringsAsFactors = FALSE
  )
}))
support_summary <- support_summary[
  match(c("early", "middle", "late"), support_summary$assigned_time_group),
  ,
  drop = FALSE
]
rownames(support_summary) <- NULL

functional_path <- file.path(
  output_dir,
  "r3_relative_location_clearance_true_functional_values.csv"
)
summary_path <- file.path(
  output_dir,
  "r3_relative_location_clearance_true_functional_summary.csv"
)
clearance_path <- file.path(
  output_dir,
  "r3_relative_location_clearance_unit_diagnostics.csv"
)
support_path <- file.path(
  output_dir,
  "r3_raised_cosine_primary_support_audit.csv"
)
support_summary_path <- file.path(
  output_dir,
  "r3_raised_cosine_primary_support_summary.csv"
)
provenance_path <- file.path(
  output_dir,
  "r3_relative_location_clearance_provenance.txt"
)
plot_path <- file.path(
  output_dir,
  "r3_relative_location_clearance_true_functional_histograms.png"
)

utils::write.csv(functional_data, functional_path, row.names = FALSE)
utils::write.csv(
  summarize_distribution(functional_data),
  summary_path,
  row.names = FALSE
)
utils::write.csv(clearance_data, clearance_path, row.names = FALSE)
utils::write.csv(support_data, support_path, row.names = FALSE)
utils::write.csv(support_summary, support_summary_path, row.names = FALSE)
writeLines(
  c(
    paste0("generated_at=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("n_variants=", n_variants),
    paste0("seeds=", paste(seed_list, collapse = ",")),
    paste0("truth_mechanisms=", paste(truth_mechanisms, collapse = ",")),
    paste0("location_truth_margin=", location_truth_margin),
    paste0(
      "location_truth_min_range_fraction=",
      location_truth_min_range_fraction
    ),
    "middle_definition=3<t<12",
    "posterior_inference=not_run",
    "local_parallel_jobs=1",
    "local_cpu_threads=1"
  ),
  provenance_path
)

functional_data$mechanism_label <- factor(
  functional_data$truth_mechanism,
  levels = truth_mechanisms,
  labels = c(
    "R3A: broad random B-spline",
    "R3B: compact raised cosine"
  )
)
functional_data$target_label <- factor(
  functional_data$target,
  levels = targets,
  labels = c("Early", "Middle", "Late", "Switch")
)
functional_data$truth_status <- factor(
  functional_data$truth_status,
  levels = c("Functional null", "Functional positive")
)

colors <- c(
  "Functional null" = "#999999",
  "Functional positive" = "#0072B2"
)
plot <- ggplot2::ggplot(
  functional_data,
  ggplot2::aes(x = functional_value, fill = truth_status)
) +
  ggplot2::geom_histogram(
    bins = 70,
    position = "identity",
    alpha = 0.65,
    color = "white",
    linewidth = 0.10
  ) +
  ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35) +
  ggplot2::facet_grid(
    rows = ggplot2::vars(mechanism_label),
    cols = ggplot2::vars(target_label),
    scales = "free_y"
  ) +
  ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
  ggplot2::labs(
    x = "True functional value on the 0.1-day evaluation grid",
    y = "Dynamic-variant count",
    fill = NULL,
    title = "R3 truth after scale-adaptive location clearance",
    subtitle = paste0(
      "Required Early/Middle/Late clearance = max(0.10, 0.10 x effect range); ",
      "five fixed seeds"
    )
  ) +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold")
  )
ggplot2::ggsave(
  plot_path,
  plot = plot,
  width = 13.5,
  height = 6.8,
  dpi = 180
)

message("Saved truth-only relative-clearance diagnostics to: ", output_dir)
