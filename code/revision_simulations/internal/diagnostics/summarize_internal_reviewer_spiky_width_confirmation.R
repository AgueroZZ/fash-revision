#!/usr/bin/env Rscript

# Aggregate the support-width sensitivity analysis that exactly preserves the
# reviewer-facing broad/spiky DGP. The existing df=16 reviewer cache is reused.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) return(".")
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find code/revision_simulations/shared/simulation_functions.R.")
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

truth_geometry <- function(beta_matrix) {
  centered <- t(apply(beta_matrix, 1, function(x) x - mean(x)))
  maximum <- apply(abs(centered), 1, max)
  data.frame(
    support_25 = rowSums(abs(centered) >= 0.25 * maximum),
    second_difference_roughness = apply(
      beta_matrix,
      1,
      function(x) sum(diff(x, differences = 2)^2) / sum((x - mean(x))^2)
    ),
    stringsAsFactors = FALSE
  )
}

summarize_curve <- function(rows, metric) {
  groups <- split(
    rows,
    list(rows$spiky_df, rows$method, rows$alpha),
    drop = TRUE
  )
  out <- lapply(groups, function(x) {
    value_summary <- summarize_mc_values(x[[metric]])
    data.frame(
      spiky_df = x$spiky_df[1],
      method = x$method[1],
      alpha = x$alpha[1],
      n_replications = length(unique(x$seed)),
      mean = value_summary[["mean"]],
      lower = pmax(0, value_summary[["lower"]]),
      upper = pmin(1, value_summary[["upper"]]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

plot_width_curve <- function(curve, file, y_label, title,
                             include_nominal = FALSE) {
  methods <- c(
    "FASH-IWP1-Raw",
    "FASH-IWP1-BF",
    "Direct-quadratic-LRT-eFDR-true-pi0"
  )
  method_colors <- c(
    "FASH-IWP1-Raw" = "#1B9E77",
    "FASH-IWP1-BF" = "#1B9E77",
    "Direct-quadratic-LRT-eFDR-true-pi0" = "#984EA3"
  )
  method_lty <- c(
    "FASH-IWP1-Raw" = 1,
    "FASH-IWP1-BF" = 2,
    "Direct-quadratic-LRT-eFDR-true-pi0" = 1
  )
  widths <- c(12L, 14L, 16L)
  png(file, width = 1800, height = 850, res = 180)
  par(
    mfrow = c(1, 3),
    mar = c(4.5, 4.4, 3.3, 1.0),
    oma = c(7.5, 0, 3.0, 0)
  )
  for (width_index in seq_along(widths)) {
    spiky_df <- widths[width_index]
    panel <- curve[curve$spiky_df == spiky_df, ]
    plot(
      NA,
      xlim = range(panel$alpha),
      ylim = c(0, if (include_nominal) max(panel$alpha) * 1.05 else 1),
      xlab = "Nominal FDR level alpha",
      ylab = if (width_index == 1) y_label else "",
      main = paste0("Local B-spline df = ", spiky_df)
    )
    grid(col = "grey88", lty = 3)
    abline(v = 0.05, col = "grey45", lty = 3)
    if (include_nominal) abline(a = 0, b = 1, col = "grey35", lty = 3)
    for (method in methods) {
      rows <- panel[panel$method == method, ]
      rows <- rows[order(rows$alpha), ]
      polygon(
        c(rows$alpha, rev(rows$alpha)),
        c(rows$lower, rev(rows$upper)),
        col = adjustcolor(method_colors[[method]], alpha.f = 0.10),
        border = NA
      )
      lines(
        rows$alpha,
        rows$mean,
        col = method_colors[[method]],
        lty = method_lty[[method]],
        lwd = 2.2
      )
    }
  }
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = rep(0, 4))
  plot.new()
  legend(
    x = 0.5,
    y = 0.035,
    xjust = 0.5,
    yjust = 0,
    horiz = TRUE,
    legend = c(
      "FASH-IWP1-Raw",
      "FASH-IWP1-BF",
      "Direct quadratic eFDR (true pi0)",
      if (include_nominal) "Nominal alpha" else NULL
    ),
    col = c(
      unname(method_colors[methods]),
      if (include_nominal) "grey35" else NULL
    ),
    lty = c(
      unname(method_lty[methods]),
      if (include_nominal) 3 else NULL
    ),
    lwd = c(rep(2.2, 3), if (include_nominal) 1.4 else NULL),
    bty = "n",
    cex = 0.80
  )
  mtext(title, outer = TRUE, cex = 1.15)
  dev.off()
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

seeds <- c(12345L, 22345L, 32345L, 42345L, 52345L)
new_dfs <- c(12L, 14L)
all_dfs <- c(12L, 14L, 16L)
methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
alpha_grid <- seq(0.005, 0.20, by = 0.005)

candidate_paths <- unlist(lapply(seeds, function(seed) {
  file.path(
    workflowr_root,
    "output",
    "revision_simulations",
    "internal",
    paste0(
      "confirm_reviewer_mixed_spiky_width_df12_14_seed",
      seed,
      "_B100"
    ),
    "candidates",
    paste0("spiky_df", new_dfs, ".rds")
  )
}))
if (any(!file.exists(candidate_paths))) {
  stop("One or more exact-DGP width caches are missing.")
}
candidates <- lapply(candidate_paths, readRDS)
candidate_keys <- expand.grid(
  spiky_df = new_dfs,
  seed = seeds,
  KEEP.OUT.ATTRS = FALSE
)
candidate_keys <- candidate_keys[order(candidate_keys$seed, candidate_keys$spiky_df), ]
names(candidates) <- paste(candidate_keys$seed, candidate_keys$spiky_df, sep = "_")

valid_candidate <- vapply(names(candidates), function(key) {
  x <- candidates[[key]]
  config <- x$configuration
  identical(config$base_dgp, "reviewer_mixed") &&
    identical(config$secondary_selection, "random") &&
    isTRUE(all.equal(config$secondary_fraction, c(0.40, 0.65))) &&
    identical(as.integer(config$efdr_permutations), 100L) &&
    identical(as.integer(x$spiky_df), as.integer(strsplit(key, "_")[[1]][2]))
}, logical(1))
if (!all(valid_candidate)) {
  stop("One or more exact-DGP width caches have invalid settings.")
}

old_summary_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "mc",
  "mixed_shape_multispike_margin020_pilot5",
  "summary"
)
old_shape <- read.csv(
  file.path(old_summary_dir, "all_replicate_global_shape_power.csv"),
  stringsAsFactors = FALSE
)
old_global <- read.csv(
  file.path(old_summary_dir, "all_replicate_global_alpha_curves.csv"),
  stringsAsFactors = FALSE
)
old_pi0 <- read.csv(
  file.path(old_summary_dir, "all_replicate_global_pi0.csv"),
  stringsAsFactors = FALSE
)

new_spiky_alpha <- do.call(rbind, lapply(names(candidates), function(key) {
  key_parts <- strsplit(key, "_", fixed = TRUE)[[1]]
  seed <- as.integer(key_parts[1])
  spiky_df <- as.integer(key_parts[2])
  candidate <- candidates[[key]]
  spiky_units <- candidate$effect_sim$unit_info$unit_index[
    candidate$effect_sim$unit_info$effect_class == "dynamic_bspline" &
      candidate$effect_sim$unit_info$shape_profile == "spiky"
  ]
  do.call(rbind, lapply(methods, function(method) {
    method_result <- candidate$result_table[
      candidate$result_table$method == method,
    ]
    do.call(rbind, lapply(alpha_grid, function(alpha) {
      selected <- method_result$unit_index[
        is.finite(method_result$adjusted_score) &
          method_result$adjusted_score <= alpha
      ]
      data.frame(
        seed = seed,
        spiky_df = spiky_df,
        method = method,
        alpha = alpha,
        power = sum(spiky_units %in% selected) / length(spiky_units),
        stringsAsFactors = FALSE
      )
    }))
  }))
}))
old_spiky_alpha <- old_shape[
  old_shape$shape_profile == "spiky" &
    old_shape$method %in% methods,
  c("seed", "method", "alpha", "power")
]
old_spiky_alpha$spiky_df <- 16L
all_spiky_alpha <- rbind(
  new_spiky_alpha,
  old_spiky_alpha[, names(new_spiky_alpha)]
)
spiky_power_curve <- summarize_curve(all_spiky_alpha, "power")

new_global_alpha <- do.call(rbind, lapply(names(candidates), function(key) {
  key_parts <- strsplit(key, "_", fixed = TRUE)[[1]]
  out <- candidates[[key]]$alpha_curve
  out$seed <- as.integer(key_parts[1])
  out$spiky_df <- as.integer(key_parts[2])
  out
}))
old_global <- old_global[old_global$method %in% methods, ]
old_global$spiky_df <- 16L
all_global_alpha <- rbind(
  new_global_alpha[, c(
    "seed", "spiky_df", "method", "alpha", "empirical_fdr", "power"
  )],
  old_global[, c(
    "seed", "spiky_df", "method", "alpha", "empirical_fdr", "power"
  )]
)
global_fdr_curve <- summarize_curve(all_global_alpha, "empirical_fdr")

alpha_005 <- all_spiky_alpha[abs(all_spiky_alpha$alpha - 0.05) < 1e-12, ]
wide_alpha_005 <- reshape(
  alpha_005,
  idvar = c("seed", "spiky_df"),
  timevar = "method",
  direction = "wide"
)
names(wide_alpha_005) <- sub("^power\\.", "", names(wide_alpha_005))
wide_alpha_005$bf_minus_direct_quadratic <-
  wide_alpha_005[["FASH-IWP1-BF"]] -
  wide_alpha_005[["Direct-quadratic-LRT-eFDR-true-pi0"]]

power_summary <- do.call(rbind, lapply(all_dfs, function(spiky_df) {
  rows <- wide_alpha_005[wide_alpha_005$spiky_df == spiky_df, ]
  raw_summary <- summarize_mc_values(rows[["FASH-IWP1-Raw"]])
  bf_summary <- summarize_mc_values(rows[["FASH-IWP1-BF"]])
  direct_summary <- summarize_mc_values(
    rows[["Direct-quadratic-LRT-eFDR-true-pi0"]]
  )
  advantage_summary <- summarize_mc_values(rows$bf_minus_direct_quadratic)
  data.frame(
    spiky_df = spiky_df,
    n_replications = nrow(rows),
    raw_power = raw_summary[["mean"]],
    bf_power = bf_summary[["mean"]],
    bf_power_lower = pmax(0, bf_summary[["lower"]]),
    bf_power_upper = pmin(1, bf_summary[["upper"]]),
    direct_quadratic_power = direct_summary[["mean"]],
    direct_power_lower = pmax(0, direct_summary[["lower"]]),
    direct_power_upper = pmin(1, direct_summary[["upper"]]),
    bf_minus_direct = advantage_summary[["mean"]],
    advantage_lower = advantage_summary[["lower"]],
    advantage_upper = advantage_summary[["upper"]],
    n_seeds_bf_above_direct = sum(rows$bf_minus_direct_quadratic > 0),
    stringsAsFactors = FALSE
  )
}))

geometry_rows <- do.call(rbind, lapply(names(candidates), function(key) {
  key_parts <- strsplit(key, "_", fixed = TRUE)[[1]]
  candidate <- candidates[[key]]
  spiky <- candidate$effect_sim$unit_info$shape_profile == "spiky"
  data.frame(
    seed = as.integer(key_parts[1]),
    spiky_df = as.integer(key_parts[2]),
    truth_geometry(candidate$effect_sim$beta_matrix[spiky, ]),
    stringsAsFactors = FALSE
  )
}))
for (seed in seeds) {
  component_seeds <- revision_component_seeds(seed)
  effect_sim <- simulate_targeted_local_bspline_effect_set(
    n_variants = 1000,
    class_probs = c(
      dynamic_bspline = 0.20,
      constant = 0.40,
      zero = 0.40
    ),
    dynamic_amplitude = 2,
    switch_threshold = 0.25,
    minimum_location_margin = 0.60,
    minimum_location_ratio = 1.4,
    non_switch_baseline_fraction = 0.75,
    non_switch_background_fraction = 0.05,
    profile = "mixed",
    spiky_truth_version = "mixed_single_double_v2",
    spiky_secondary_fraction = c(0.40, 0.65),
    spiky_bspline_df = 16,
    spiky_secondary_selection = "random",
    spiky_minimum_peak_separation = 3,
    spiky_non_switch_baseline_fraction = 0.20,
    target_centered_rms = 0.90,
    exact_class_counts = TRUE,
    seed = component_seeds[["functional_truth"]]
  )
  spiky <- effect_sim$unit_info$shape_profile == "spiky"
  geometry_rows <- rbind(
    geometry_rows,
    data.frame(
      seed = seed,
      spiky_df = 16L,
      truth_geometry(effect_sim$beta_matrix[spiky, ]),
      stringsAsFactors = FALSE
    )
  )
}
geometry_summary <- aggregate(
  cbind(support_25, second_difference_roughness) ~ spiky_df,
  data = geometry_rows,
  FUN = mean
)

bf_pi0_rows <- do.call(rbind, lapply(names(candidates), function(key) {
  key_parts <- strsplit(key, "_", fixed = TRUE)[[1]]
  data.frame(
    seed = as.integer(key_parts[1]),
    spiky_df = as.integer(key_parts[2]),
    bf_pi0 = candidates[[key]]$summary$bf_pi0,
    stringsAsFactors = FALSE
  )
}))
old_bf_pi0 <- old_pi0[
  old_pi0$method == "FASH-IWP1" &
    old_pi0$fit == "BF-corrected",
  c("seed", "estimated_pi0")
]
names(old_bf_pi0)[2] <- "bf_pi0"
old_bf_pi0$spiky_df <- 16L
bf_pi0_rows <- rbind(
  bf_pi0_rows,
  old_bf_pi0[, names(bf_pi0_rows)]
)
pi0_summary <- aggregate(bf_pi0 ~ spiky_df, data = bf_pi0_rows, FUN = mean)

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "confirm_reviewer_mixed_spiky_width_df12_14_16_B100_5seed"
)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(all_spiky_alpha, file.path(output_dir, "all_spiky_alpha.csv"))
write_csv(spiky_power_curve, file.path(output_dir, "mc_spiky_power_curve.csv"))
write_csv(all_global_alpha, file.path(output_dir, "all_global_alpha.csv"))
write_csv(global_fdr_curve, file.path(output_dir, "mc_global_fdr_curve.csv"))
write_csv(wide_alpha_005, file.path(output_dir, "paired_power_alpha005.csv"))
write_csv(power_summary, file.path(output_dir, "mc_power_alpha005.csv"))
write_csv(geometry_rows, file.path(output_dir, "truth_geometry.csv"))
write_csv(geometry_summary, file.path(output_dir, "truth_geometry_summary.csv"))
write_csv(bf_pi0_rows, file.path(output_dir, "bf_pi0_by_seed.csv"))
write_csv(pi0_summary, file.path(output_dir, "bf_pi0_summary.csv"))

plot_width_curve(
  curve = spiky_power_curve,
  file = file.path(figure_dir, "spiky_power_across_alpha_by_width.png"),
  y_label = "Mean spiky-effect power",
  title = "Exact reviewer DGP: spiky-effect power by support width"
)
plot_width_curve(
  curve = global_fdr_curve,
  file = file.path(figure_dir, "global_fdr_across_alpha_by_width.png"),
  y_label = "Monte Carlo estimate of FDR",
  title = "Exact reviewer DGP: global empirical FDR by support width",
  include_nominal = TRUE
)

png(
  file.path(figure_dir, "geometry_and_power_alpha005.png"),
  width = 1500,
  height = 650,
  res = 180
)
par(mfrow = c(1, 3), mar = c(4.3, 4.5, 3.0, 1.0),
    oma = c(0, 0, 2.5, 0))
plot(
  geometry_summary$spiky_df,
  geometry_summary$support_25,
  type = "b",
  pch = 19,
  lwd = 2,
  col = "#377EB8",
  xlab = "Local B-spline df",
  ylab = "Mean 25%-peak support",
  main = "Temporal support"
)
plot(
  geometry_summary$spiky_df,
  geometry_summary$second_difference_roughness,
  type = "b",
  pch = 19,
  lwd = 2,
  col = "#E41A1C",
  xlab = "Local B-spline df",
  ylab = "Mean second-difference roughness",
  main = "Curve roughness"
)
plot(
  power_summary$spiky_df,
  power_summary$bf_power,
  type = "b",
  pch = 19,
  lwd = 2,
  col = "#1B9E77",
  lty = 2,
  ylim = c(0, 1),
  xlab = "Local B-spline df",
  ylab = "Mean spiky-effect power",
  main = "Power at alpha = 0.05"
)
lines(
  power_summary$spiky_df,
  power_summary$direct_quadratic_power,
  type = "b",
  pch = 19,
  lwd = 2,
  col = "#984EA3"
)
legend(
  "bottomleft",
  legend = c("FASH-IWP1-BF", "Direct quadratic eFDR"),
  col = c("#1B9E77", "#984EA3"),
  lty = c(2, 1),
  lwd = 2,
  pch = 19,
  bty = "n",
  cex = 0.72
)
mtext(
  "Exact reviewer-facing DGP: support-width sensitivity",
  outer = TRUE,
  cex = 1.1
)
dev.off()

saveRDS(
  list(
    seeds = seeds,
    spiky_dfs = all_dfs,
    all_spiky_alpha = all_spiky_alpha,
    spiky_power_curve = spiky_power_curve,
    all_global_alpha = all_global_alpha,
    global_fdr_curve = global_fdr_curve,
    paired_power_alpha005 = wide_alpha_005,
    power_summary_alpha005 = power_summary,
    geometry_rows = geometry_rows,
    geometry_summary = geometry_summary,
    bf_pi0_rows = bf_pi0_rows,
    pi0_summary = pi0_summary
  ),
  file.path(output_dir, "confirmation_summary.rds")
)

message("Saved exact reviewer-DGP width confirmation to: ", output_dir)
