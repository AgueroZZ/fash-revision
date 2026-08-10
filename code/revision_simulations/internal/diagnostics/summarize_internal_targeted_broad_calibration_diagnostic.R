#!/usr/bin/env Rscript

# Diagnose why the targeted broad B-spline functional simulation has calibrated
# switch FSR while the compact zero-background peak simulation does not.

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

targeted_id <- "internal_targeted_broad_calibration_pilot5"
old_targeted_id <- "functional_targeted_broad_baseline075_pilot5"
spiky_id <- "sparse_timed_cosine_functional_pilot5"
mc_root <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "mc"
)
targeted_dir <- file.path(mc_root, targeted_id)
old_targeted_dir <- file.path(mc_root, old_targeted_id)
spiky_dir <- file.path(mc_root, spiky_id)

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "targeted_broad_functional_calibration_diagnostic"
)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

read_summary <- function(output_dir, filename) {
  read.csv(
    file.path(output_dir, "summary", filename),
    stringsAsFactors = FALSE
  )
}

targeted_alpha <- read_summary(
  targeted_dir,
  "functional_testing_mc_alpha_curve.csv"
)
targeted_alpha_005 <- read_summary(
  targeted_dir,
  "functional_testing_mc_alpha005_summary.csv"
)
old_targeted_alpha_005 <- read_summary(
  old_targeted_dir,
  "functional_testing_mc_alpha005_summary.csv"
)
spiky_alpha <- read_summary(
  spiky_dir,
  "functional_testing_mc_alpha_curve.csv"
)
spiky_alpha_005 <- read_summary(
  spiky_dir,
  "functional_testing_mc_alpha005_summary.csv"
)
targeted_calls <- read_summary(
  targeted_dir,
  "all_replicate_functional_call_diagnostics_alpha005.csv"
)

reproduction_columns <- c(
  "run",
  "target",
  "method",
  "mean_dynamic_discoveries",
  "mean_discoveries",
  "mean_false_discoveries",
  "mean_power",
  "mean_empirical_fsr"
)
reproduction_rows <- rbind(
  transform(
    old_targeted_alpha_005,
    run = "Original cached run"
  )[, reproduction_columns],
  transform(
    targeted_alpha_005,
    run = "Current-code reproduction"
  )[, reproduction_columns]
)
write_csv(
  reproduction_rows,
  file.path(output_dir, "original_vs_reproduced_alpha005.csv")
)

false_switch_calls <- targeted_calls[
  targeted_calls$target == "switch" &
    targeted_calls$method == "FASH-IWP1-BF" &
    targeted_calls$true_null,
  ,
  drop = FALSE
]
false_switch_calls$false_call_source <- ifelse(
  false_switch_calls$effect_class == "dynamic_bspline",
  "True-dynamic non-switch",
  "First-stage dynamic null"
)
false_switch_counts <- aggregate(
  list(n_false_calls = rep(1L, nrow(false_switch_calls))),
  by = false_switch_calls[, c("seed", "false_call_source")],
  FUN = sum
)
write_csv(
  false_switch_calls,
  file.path(output_dir, "bf_false_switch_calls_alpha005.csv")
)
write_csv(
  false_switch_counts,
  file.path(output_dir, "bf_false_switch_counts_by_seed.csv")
)
bf_switch_calls <- targeted_calls[
  targeted_calls$target == "switch" &
    targeted_calls$method == "FASH-IWP1-BF",
  ,
  drop = FALSE
]
bf_switch_calls$call_truth_class <- ifelse(
  !bf_switch_calls$true_null,
  "True switch",
  ifelse(
    bf_switch_calls$effect_class == "dynamic_bspline",
    "True-dynamic non-switch",
    "First-stage dynamic null"
  )
)
call_lfsr_summary <- do.call(
  rbind,
  lapply(split(bf_switch_calls, bf_switch_calls$call_truth_class), function(x) {
    data.frame(
      call_truth_class = x$call_truth_class[1],
      n_calls = nrow(x),
      lfsr_mean = mean(x$lfsr),
      lfsr_median = stats::median(x$lfsr),
      lfsr_minimum = min(x$lfsr),
      lfsr_maximum = max(x$lfsr),
      stringsAsFactors = FALSE
    )
  })
)
rownames(call_lfsr_summary) <- NULL
write_csv(
  call_lfsr_summary,
  file.path(output_dir, "bf_switch_call_lfsr_by_truth.csv")
)

make_geometry <- function(beta_evaluation,
                          evaluation_grid,
                          unit_info,
                          true_functionals,
                          seed,
                          setting) {
  dynamic <- unit_info$effect_class == "dynamic_bspline"
  beta <- beta_evaluation[dynamic, , drop = FALSE]
  info <- unit_info[dynamic, , drop = FALSE]
  functionals <- true_functionals[dynamic, , drop = FALSE]
  curve_min <- apply(beta, 1, min)
  curve_max <- apply(beta, 1, max)
  same_sign_margin <- ifelse(
    curve_min >= 0,
    curve_min,
    ifelse(curve_max <= 0, -curve_max, 0)
  )
  data.frame(
    setting = setting,
    seed = seed,
    variant_id = info$variant_id,
    time_group = info$time_group,
    switch_status = ifelse(
      functionals[, "switch"] > 0,
      "switch",
      "non-switch"
    ),
    shape_pattern = if ("spike_pattern" %in% names(info)) {
      info$spike_pattern
    } else {
      NA_character_
    },
    maximum_absolute_effect = apply(abs(beta), 1, max),
    minimum_absolute_effect = apply(abs(beta), 1, min),
    same_sign_margin = same_sign_margin,
    fraction_grid_within_switch_threshold =
      rowMeans(abs(beta) <= 0.25),
    maximum_positive_effect = pmax(curve_max, 0),
    maximum_negative_magnitude = pmax(-curve_min, 0),
    true_switch_functional = functionals[, "switch"],
    stringsAsFactors = FALSE
  )
}

seed_list <- readRDS(
  file.path(targeted_dir, "configuration.rds")
)$seed_list
targeted_effects <- list()
targeted_geometry <- lapply(seed_list, function(seed) {
  compact <- readRDS(
    file.path(targeted_dir, "replicates", paste0("seed_", seed, ".rds"))
  )
  config <- compact$configuration
  effect_sim <- simulate_targeted_local_bspline_effect_set(
    n_variants = config$J,
    time_grid = config$time_grid,
    evaluation_grid = config$evaluation_grid,
    class_probs = config$class_probs,
    dynamic_amplitude = config$dynamic_amplitude,
    switch_threshold = config$switch_threshold,
    minimum_location_margin = config$minimum_location_margin,
    minimum_location_ratio = config$minimum_location_ratio,
    non_switch_baseline_fraction = config$non_switch_baseline_fraction,
    non_switch_background_fraction =
      config$non_switch_background_fraction,
    profile = config$targeted_profile,
    seed = compact$component_seeds[["functional_truth"]],
    scenario = config$scenario
  )
  targeted_effects[[as.character(seed)]] <<- effect_sim
  make_geometry(
    beta_evaluation = effect_sim$beta_evaluation,
    evaluation_grid = effect_sim$evaluation_grid,
    unit_info = effect_sim$unit_info,
    true_functionals = effect_sim$true_functionals,
    seed = seed,
    setting = "Targeted broad; baseline 0.75"
  )
})
targeted_geometry <- do.call(rbind, targeted_geometry)

find_spiky_raw <- function(seed) {
  candidates <- list.files(
    file.path(
      workflowr_root,
      "output",
      "revision_simulations",
      "raw"
    ),
    pattern = paste0(
      "^genotype_sparse_timed_cosine_one_two_peak_dynamic_eqtl_",
      ".*_w1p5_rms0p9_.*_seed",
      seed,
      "[.]rds$"
    ),
    full.names = TRUE
  )
  if (length(candidates) != 1) {
    stop("Expected one width-1.5 spiky raw output for seed ", seed, ".")
  }
  candidates
}

spiky_outputs <- list()
spiky_geometry <- lapply(seed_list, function(seed) {
  out <- readRDS(find_spiky_raw(seed))
  spiky_outputs[[as.character(seed)]] <<- out
  make_geometry(
    beta_evaluation = out$true_beta_evaluation,
    evaluation_grid = out$evaluation_grid,
    unit_info = out$unit_info,
    true_functionals = out$true_functionals,
    seed = seed,
    setting = "Compact peaks; zero background"
  )
})
spiky_geometry <- do.call(rbind, spiky_geometry)
geometry <- rbind(targeted_geometry, spiky_geometry)
write_csv(geometry, file.path(output_dir, "dynamic_truth_geometry.csv"))

quantiles <- function(x) {
  stats::quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE)
}
geometry_summary <- do.call(
  rbind,
  lapply(
    split(geometry, list(geometry$setting, geometry$switch_status)),
    function(rows) {
      if (nrow(rows) == 0) return(NULL)
      near_zero <- quantiles(rows$fraction_grid_within_switch_threshold)
      sign_margin <- quantiles(rows$same_sign_margin)
      max_abs <- quantiles(rows$maximum_absolute_effect)
      data.frame(
        setting = rows$setting[1],
        switch_status = rows$switch_status[1],
        n_curves = nrow(rows),
        near_zero_fraction_q25 = near_zero[1],
        near_zero_fraction_median = near_zero[2],
        near_zero_fraction_q75 = near_zero[3],
        same_sign_margin_q25 = sign_margin[1],
        same_sign_margin_median = sign_margin[2],
        same_sign_margin_q75 = sign_margin[3],
        maximum_absolute_effect_q25 = max_abs[1],
        maximum_absolute_effect_median = max_abs[2],
        maximum_absolute_effect_q75 = max_abs[3],
        stringsAsFactors = FALSE
      )
    }
  )
)
rownames(geometry_summary) <- NULL
write_csv(
  geometry_summary,
  file.path(output_dir, "dynamic_truth_geometry_summary.csv")
)

plot_rows <- rbind(
  transform(
    targeted_alpha[
      targeted_alpha$target == "switch" &
        targeted_alpha$method == "FASH-IWP1-BF",
      ,
      drop = FALSE
    ],
    setting = "Targeted broad; baseline 0.75"
  ),
  transform(
    spiky_alpha[
      spiky_alpha$target == "switch" &
        spiky_alpha$method == "FASH-IWP1-BF",
      ,
      drop = FALSE
    ],
    setting = "Compact peaks; zero background"
  )
)
comparison_colors <- c(
  "Targeted broad; baseline 0.75" = "#0072B2",
  "Compact peaks; zero background" = "#D55E00"
)
png(
  file.path(figure_dir, "bf_switch_power_fsr_comparison.png"),
  width = 1300,
  height = 620,
  res = 180
)
old_par <- par(mfrow = c(1, 2), mar = c(4.8, 4.8, 3.8, 1.0))
for (metric in c("mean_power", "mean_empirical_fsr")) {
  plot(
    NA,
    xlim = range(plot_rows$alpha),
    ylim = c(0, 1),
    xlab = "Nominal FSR level alpha",
    ylab = if (metric == "mean_power") "Switch power" else "Empirical FSR",
    main = if (metric == "mean_power") "BF switch power" else "BF switch FSR",
    las = 1
  )
  grid(col = "#E6E6E6")
  if (metric == "mean_empirical_fsr") {
    abline(0, 1, lty = 3, col = "#555555")
  }
  for (setting in names(comparison_colors)) {
    rows <- plot_rows[plot_rows$setting == setting, , drop = FALSE]
    rows <- rows[order(rows$alpha), ]
    lines(
      rows$alpha,
      rows[[metric]],
      col = comparison_colors[[setting]],
      lwd = 2.5
    )
  }
  if (metric == "mean_power") {
    legend(
      "bottomright",
      legend = names(comparison_colors),
      col = comparison_colors,
      lwd = 2.5,
      bty = "n",
      cex = 0.82
    )
  }
}
par(old_par)
dev.off()

non_switch_geometry <- geometry[
  geometry$switch_status == "non-switch",
  ,
  drop = FALSE
]
setting_order <- c(
  "Targeted broad; baseline 0.75",
  "Compact peaks; zero background"
)
short_setting_labels <- c(
  "Targeted broad\nbaseline 0.75",
  "Compact peaks\nzero background"
)
non_switch_geometry$setting_factor <- factor(
  non_switch_geometry$setting,
  levels = setting_order
)
png(
  file.path(figure_dir, "non_switch_truth_geometry_comparison.png"),
  width = 1700,
  height = 620,
  res = 180
)
old_par <- par(mfrow = c(1, 3), mar = c(8.2, 4.8, 3.8, 1.0))
geometry_metrics <- c(
  fraction_grid_within_switch_threshold =
    "Fraction of grid with |effect| <= 0.25",
  same_sign_margin = "Same-sign margin",
  maximum_absolute_effect = "Maximum absolute effect"
)
for (metric in names(geometry_metrics)) {
  boxplot(
    non_switch_geometry[[metric]] ~ non_switch_geometry$setting_factor,
    names = short_setting_labels,
    col = unname(comparison_colors[setting_order]),
    border = "#444444",
    ylab = geometry_metrics[[metric]],
    xlab = "",
    main = geometry_metrics[[metric]],
    las = 1,
    cex.axis = 0.8
  )
  grid(nx = NA, ny = NULL, col = "#E6E6E6")
}
par(old_par)
dev.off()

call_class_order <- c(
  "True switch",
  "True-dynamic non-switch",
  "First-stage dynamic null"
)
bf_switch_calls$call_truth_factor <- factor(
  bf_switch_calls$call_truth_class,
  levels = call_class_order
)
png(
  file.path(figure_dir, "bf_switch_call_lfsr_by_truth.png"),
  width = 900,
  height = 720,
  res = 180
)
boxplot(
  lfsr ~ call_truth_factor,
  data = bf_switch_calls,
  names = c("True switch", "Dynamic\nnon-switch", "First-stage\nnull"),
  col = c("#009E73", "#E69F00", "#999999"),
  border = "#444444",
  ylim = c(0, 1),
  xlab = "",
  ylab = "Individual switch lfsr",
  main = "BF switch calls at alpha = 0.05",
  las = 1
)
grid(nx = NA, ny = NULL, col = "#E6E6E6")
dev.off()

false_switch_calls <- false_switch_calls[
  order(false_switch_calls$seed, false_switch_calls$variant_id),
  ,
  drop = FALSE
]
png(
  file.path(figure_dir, "bf_false_switch_call_examples.png"),
  width = 1500,
  height = 1900,
  res = 180
)
n_panels <- nrow(false_switch_calls)
old_par <- par(
  mfrow = c(ceiling(n_panels / 2), 2),
  mar = c(3.5, 4.1, 3.5, 1.0)
)
for (index in seq_len(n_panels)) {
  call <- false_switch_calls[index, , drop = FALSE]
  seed <- as.character(call$seed)
  out <- readRDS(
    file.path(
      targeted_dir,
      "full_fits",
      paste0("seed_", seed, ".rds")
    )
  )
  effect_sim <- targeted_effects[[seed]]
  unit_index <- match(call$variant_id, out$unit_info$variant_id)
  effect_index <- match(call$variant_id, effect_sim$unit_info$variant_id)
  observed <- out$eqtl_summary$beta_hat[unit_index, ]
  se <- out$eqtl_summary$se[unit_index, ]
  true_curve <- effect_sim$beta_evaluation[effect_index, ]
  y_limits <- range(c(observed - 2 * se, observed + 2 * se, true_curve))
  plot(
    out$settings$time_grid,
    observed,
    type = "n",
    ylim = y_limits,
    xlab = "Time",
    ylab = "Genetic effect",
    main = paste0(
      "Seed ", seed, "; ", call$variant_id, "; ", call$false_call_source,
      "\nlfsr=", formatC(call$lfsr, digits = 3, format = "f"),
      "; cFSR=", formatC(call$cfsr, digits = 3, format = "f"),
      "; true functional=",
      formatC(call$true_functional, digits = 3, format = "f")
    )
  )
  arrows(
    out$settings$time_grid,
    observed - 2 * se,
    out$settings$time_grid,
    observed + 2 * se,
    angle = 90,
    code = 3,
    length = 0.025,
    col = "#777777"
  )
  points(out$settings$time_grid, observed, pch = 19, cex = 0.65)
  lines(
    effect_sim$evaluation_grid,
    true_curve,
    col = "#D55E00",
    lwd = 2.2
  )
  abline(h = c(-0.25, 0, 0.25), lty = c(3, 2, 3), col = "#666666")
}
par(old_par)
dev.off()

targeted_bf <- targeted_alpha_005[
  targeted_alpha_005$method == "FASH-IWP1-BF",
  ,
  drop = FALSE
]
spiky_bf <- spiky_alpha_005[
  spiky_alpha_005$method == "FASH-IWP1-BF",
  ,
  drop = FALSE
]
print(targeted_bf[
  ,
  c(
    "target",
    "mean_discoveries",
    "mean_false_discoveries",
    "mean_power",
    "mean_conditional_empirical_fsr",
    "mean_empirical_fsr"
  )
])
print(false_switch_counts)
print(call_lfsr_summary)
print(geometry_summary)
print(spiky_bf[
  spiky_bf$target == "switch",
  c(
    "target",
    "mean_discoveries",
    "mean_false_discoveries",
    "mean_power",
    "mean_conditional_empirical_fsr",
    "mean_empirical_fsr"
  )
])
