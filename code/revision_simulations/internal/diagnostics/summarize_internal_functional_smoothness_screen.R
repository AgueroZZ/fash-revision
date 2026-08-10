#!/usr/bin/env Rscript

# Summarize the internal switch-functional smoothness and null-separation
# experiment. This script only reads cached five-seed results.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) return(".")
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository.")
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

cache_specs <- data.frame(
  setting_id = c(
    "sharp_zero_boundary",
    "moderate_zero_boundary",
    "broad_zero_boundary",
    "sharp_minabs_050",
    "broad_minabs_050",
    "broad_minabs_0925"
  ),
  setting_label = c(
    "Sharp, min |beta| = 0",
    "Moderate, min |beta| = 0",
    "Broad, min |beta| = 0",
    "Sharp, min |beta| = 0.5",
    "Broad, min |beta| = 0.5",
    "Broad, min |beta| = 0.925"
  ),
  width_half = c(1.50, 2.25, 3.00, 1.50, 3.00, 3.00),
  non_switch_min_abs_effect = c(0, 0, 0, 0.5, 0.5, 0.925),
  output_id = c(
    "sparse_timed_cosine_functional_pilot5",
    "internal_functional_smoothness_w225_pilot5",
    "internal_functional_smoothness_w300_pilot5",
    "internal_functional_smoothness_w150_minabs050_pilot5",
    "internal_functional_smoothness_w300_minabs050_pilot5",
    "internal_functional_smoothness_w300_minabs0925_pilot5"
  ),
  stringsAsFactors = FALSE
)

read_setting <- function(specification) {
  output_dir <- file.path(
    workflowr_root,
    "output",
    "revision_simulations",
    "mc",
    specification$output_id
  )
  configuration <- readRDS(file.path(output_dir, "configuration.rds"))
  curve <- read.csv(
    file.path(output_dir, "summary", "functional_testing_mc_alpha_curve.csv"),
    stringsAsFactors = FALSE
  )
  pi0 <- read.csv(
    file.path(output_dir, "summary", "functional_testing_mc_pi0_summary.csv"),
    stringsAsFactors = FALSE
  )
  expected_minimum <- if (
    "non_switch_min_abs_effect" %in% names(configuration)
  ) {
    configuration$non_switch_min_abs_effect
  } else {
    0
  }
  if (!isTRUE(all.equal(configuration$J, 1000L)) ||
      !isTRUE(all.equal(configuration$n_donors, 19L)) ||
      !isTRUE(all.equal(configuration$n_covariates, 5L)) ||
      !isTRUE(all.equal(configuration$expression_noise_sd, 1)) ||
      !isTRUE(all.equal(configuration$width_half, specification$width_half)) ||
      !isTRUE(all.equal(
        expected_minimum,
        specification$non_switch_min_abs_effect
      )) ||
      !isTRUE(all.equal(configuration$target_centered_rms, 0.9)) ||
      !isTRUE(all.equal(configuration$shape_cell_counts, c(
        k1__spiky__single = 100,
        `k2__spiky__same-sign` = 50,
        `k2__spiky__alternating-sign` = 50
      ))) ||
      length(configuration$seed_list) != 5L) {
    stop("A functional smoothness cache has unexpected settings: ",
      specification$output_id)
  }
  curve <- curve[curve$target == "switch", , drop = FALSE]
  curve$setting_id <- specification$setting_id
  curve$setting_label <- specification$setting_label
  curve$width_half <- specification$width_half
  curve$non_switch_min_abs_effect <-
    specification$non_switch_min_abs_effect
  pi0$setting_id <- specification$setting_id
  pi0$setting_label <- specification$setting_label
  pi0$width_half <- specification$width_half
  pi0$non_switch_min_abs_effect <-
    specification$non_switch_min_abs_effect
  list(curve = curve, pi0 = pi0)
}

results <- lapply(
  seq_len(nrow(cache_specs)),
  function(index) read_setting(cache_specs[index, , drop = FALSE])
)
combined_curve <- do.call(rbind, lapply(results, `[[`, "curve"))
combined_pi0 <- do.call(rbind, lapply(results, `[[`, "pi0"))
alpha_005 <- combined_curve[
  abs(combined_curve$alpha - 0.05) < 1e-12,
  ,
  drop = FALSE
]

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "functional_smoothness_lfsr_screen"
)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(cache_specs, file.path(output_dir, "setting_manifest.csv"))
write_csv(combined_curve, file.path(output_dir, "switch_alpha_curves.csv"))
write_csv(alpha_005, file.path(output_dir, "switch_alpha005_summary.csv"))
write_csv(combined_pi0, file.path(output_dir, "pi0_summary.csv"))

method_colors <- c(
  "FASH-IWP1-Raw" = "#009E73",
  "FASH-IWP1-BF" = "#0072B2"
)
setting_colors <- c(
  "Sharp, min |beta| = 0" = "#D55E00",
  "Moderate, min |beta| = 0" = "#CC79A7",
  "Broad, min |beta| = 0" = "#7A5195",
  "Sharp, min |beta| = 0.5" = "#E69F00",
  "Broad, min |beta| = 0.5" = "#56B4E9",
  "Broad, min |beta| = 0.925" = "#009E73"
)

plot_metric_panels <- function(data,
                               settings,
                               file,
                               title) {
  metrics <- c(
    mean_power = "Switch power",
    mean_conditional_empirical_fsr = "Conditional empirical FSR",
    mean_empirical_fsr = "End-to-end false-call proportion"
  )
  png(file, width = 1800, height = 620, res = 180)
  old_par <- par(
    mfrow = c(1, 3),
    mar = c(4.6, 4.8, 3.8, 1.0),
    oma = c(0, 0, 2.2, 0)
  )
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  for (metric in names(metrics)) {
    y_limit <- if (metric == "mean_power") c(0, 1) else c(0, 0.8)
    plot(
      NA,
      xlim = range(data$alpha),
      ylim = y_limit,
      xlab = "Nominal FSR level alpha",
      ylab = metrics[[metric]],
      main = metrics[[metric]],
      las = 1
    )
    grid(col = "#E6E6E6")
    if (metric != "mean_power") {
      abline(0, 1, lty = 3, lwd = 1.5, col = "#555555")
    }
    for (setting in settings) {
      subset_data <- data[
        data$setting_label == setting &
          data$method == "FASH-IWP1-BF",
        ,
        drop = FALSE
      ]
      subset_data <- subset_data[order(subset_data$alpha), ]
      lines(
        subset_data$alpha,
        subset_data[[metric]],
        col = setting_colors[[setting]],
        lwd = 2.4
      )
    }
    if (metric == "mean_power") {
      legend(
        "bottomright",
        legend = settings,
        col = setting_colors[settings],
        lwd = 2.4,
        bty = "n",
        cex = 0.73
      )
    }
  }
  mtext(title, outer = TRUE, cex = 1.15, font = 2)
}

plot_metric_panels(
  data = combined_curve,
  settings = cache_specs$setting_label[1:3],
  file = file.path(figure_dir, "width_ladder_bf.png"),
  title = "BF-updated switch testing: widening zero-boundary peaks"
)
plot_metric_panels(
  data = combined_curve,
  settings = cache_specs$setting_label[c(1, 3, 4, 5, 6)],
  file = file.path(figure_dir, "null_separation_bf.png"),
  title = "BF-updated switch testing: smoothness and null separation"
)

make_truth <- function(width_half, non_switch_min_abs_effect, seed = 12345) {
  component_seeds <- revision_component_seeds(seed)
  effect_sim <- simulate_raised_cosine_multipeak_effect_set(
    n_variants = 1000,
    time_grid = make_time_grid(),
    evaluation_grid = seq(0, 15, by = 0.1),
    class_probs = c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40),
    width_levels = c(spiky = width_half),
    spike_counts = 1:2,
    shape_cell_probs = c(
      k1__spiky__single = 0.50,
      `k2__spiky__same-sign` = 0.25,
      `k2__spiky__alternating-sign` = 0.25
    ),
    primary_time_groups = c("early", "middle", "late"),
    center_by_observed_mean = FALSE,
    switch_threshold = 0.25,
    relative_amplitude_range = c(0.35, 0.60),
    target_centered_rms = 0.9,
    constant_sd = 1,
    dynamic_baseline_sd = 0,
    exact_class_counts = TRUE,
    seed = seed,
    class_seed = component_seeds[["classes"]],
    constant_seed = component_seeds[["constant_effects"]],
    shape_seed = component_seeds[["functional_truth"]]
  )
  if (non_switch_min_abs_effect > 0) {
    non_switch <- effect_sim$unit_info$effect_class == "dynamic_bspline" &
      effect_sim$unit_info$sign_pattern != "alternating-sign"
    for (index in which(non_switch)) {
      direction <- sign(effect_sim$unit_info$peak_signs[[index]][1])
      current_minimum <- min(
        direction * effect_sim$beta_evaluation[index, ]
      )
      offset <- max(0, non_switch_min_abs_effect - current_minimum)
      effect_sim$beta_evaluation[index, ] <-
        effect_sim$beta_evaluation[index, ] + direction * offset
    }
  }
  effect_sim
}

example_specs <- data.frame(
  label = c(
    "Sharp, zero boundary",
    "Broad, zero boundary",
    "Broad, separated null"
  ),
  width_half = c(1.5, 3.0, 3.0),
  non_switch_min_abs_effect = c(0, 0, 0.925),
  stringsAsFactors = FALSE
)
png(
  file.path(figure_dir, "representative_truth_curves.png"),
  width = 1300,
  height = 1250,
  res = 170
)
old_par <- par(mfrow = c(3, 2), mar = c(3.8, 4.2, 3.0, 1.0))
for (row in seq_len(nrow(example_specs))) {
  truth <- make_truth(
    width_half = example_specs$width_half[row],
    non_switch_min_abs_effect =
      example_specs$non_switch_min_abs_effect[row]
  )
  for (pattern in c("single", "alternating-sign")) {
    available <- which(
      truth$unit_info$effect_class == "dynamic_bspline" &
        truth$unit_info$sign_pattern == pattern
    )
    index <- available[1]
    plot(
      truth$evaluation_grid,
      truth$beta_evaluation[index, ],
      type = "l",
      lwd = 2.4,
      col = if (pattern == "single") "#D55E00" else "#0072B2",
      xlab = "Time",
      ylab = "True effect",
      main = paste0(
        example_specs$label[row],
        ": ",
        if (pattern == "single") "non-switch" else "switch"
      ),
      las = 1
    )
    abline(h = c(-0.25, 0, 0.25), lty = c(3, 2, 3), col = "#777777")
  }
}
par(old_par)
dev.off()

print(alpha_005[
  ,
  c(
    "setting_label", "method", "mean_power",
    "mean_conditional_empirical_fsr", "mean_empirical_fsr"
  )
])
