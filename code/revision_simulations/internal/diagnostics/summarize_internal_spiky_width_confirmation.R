#!/usr/bin/env Rscript

# Aggregate the locked five-seed support-width sensitivity analysis for
# localized spiky effects. This script writes internal diagnostics only.

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

summarize_grouped_values <- function(data, group_names, value_names,
                                     confidence_level = 0.95) {
  group_key <- interaction(data[group_names], drop = TRUE, lex.order = TRUE)
  groups <- split(data, group_key)
  out <- lapply(groups, function(x) {
    row <- x[1, group_names, drop = FALSE]
    for (value_name in value_names) {
      value_summary <- summarize_mc_values(x[[value_name]], confidence_level)
      row[[paste0("mean_", value_name)]] <- value_summary[["mean"]]
      row[[paste0(value_name, "_sd")]] <- value_summary[["sd"]]
      row[[paste0(value_name, "_mc_se")]] <- value_summary[["se"]]
      row[[paste0(value_name, "_ci_lower")]] <- value_summary[["lower"]]
      row[[paste0(value_name, "_ci_upper")]] <- value_summary[["upper"]]
    }
    row$n_replications <- length(unique(x$seed))
    row
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

plot_width_power <- function(mc_curve, file, methods, method_colors,
                             method_lty) {
  widths <- sort(unique(mc_curve$spiky_df))
  display_labels <- c(
    "FASH-IWP1-Raw",
    "FASH-IWP1-BF",
    "Direct quadratic eFDR (true pi0)"
  )
  png(file, width = 1800, height = 850, res = 180)
  par(mfrow = c(1, length(widths)), mar = c(4.5, 4.4, 3.3, 1.0),
      oma = c(7.5, 0, 3.0, 0))
  for (width_index in seq_along(widths)) {
    spiky_df <- widths[width_index]
    panel <- mc_curve[mc_curve$spiky_df == spiky_df, ]
    plot(
      NA,
      xlim = range(panel$alpha),
      ylim = c(0, 1),
      xlab = "Nominal FDR level alpha",
      ylab = if (width_index == 1) "Mean spiky-effect power" else "",
      main = paste0("Local B-spline df = ", spiky_df)
    )
    grid(col = "grey88", lty = 3)
    abline(v = 0.05, col = "grey45", lty = 3)
    for (method in methods) {
      rows <- panel[panel$method == method, ]
      rows <- rows[order(rows$alpha), ]
      polygon(
        c(rows$alpha, rev(rows$alpha)),
        c(pmax(0, rows$power_ci_lower), rev(pmin(1, rows$power_ci_upper))),
        col = adjustcolor(method_colors[[method]], alpha.f = 0.10),
        border = NA
      )
      lines(
        rows$alpha,
        rows$mean_power,
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
    legend = display_labels,
    col = unname(method_colors[methods]),
    lty = unname(method_lty[methods]),
    lwd = 2.2,
    bty = "n",
    cex = 0.82
  )
  mtext(
    "Paired support-width sensitivity: spiky-effect power",
    outer = TRUE,
    cex = 1.15
  )
  dev.off()
}

plot_width_fdr <- function(mc_curve, file, methods, method_colors,
                           method_lty) {
  widths <- sort(unique(mc_curve$spiky_df))
  display_labels <- c(
    "FASH-IWP1-Raw",
    "FASH-IWP1-BF",
    "Direct quadratic eFDR (true pi0)"
  )
  png(file, width = 1800, height = 850, res = 180)
  par(mfrow = c(1, length(widths)), mar = c(4.5, 4.4, 3.3, 1.0),
      oma = c(7.5, 0, 3.0, 0))
  for (width_index in seq_along(widths)) {
    spiky_df <- widths[width_index]
    panel <- mc_curve[mc_curve$spiky_df == spiky_df, ]
    plot(
      NA,
      xlim = range(panel$alpha),
      ylim = c(0, max(panel$alpha) * 1.05),
      xlab = "Nominal FDR level alpha",
      ylab = if (width_index == 1) "Monte Carlo estimate of FDR" else "",
      main = paste0("Local B-spline df = ", spiky_df)
    )
    grid(col = "grey88", lty = 3)
    abline(a = 0, b = 1, col = "grey35", lty = 3)
    abline(v = 0.05, col = "grey45", lty = 3)
    for (method in methods) {
      rows <- panel[panel$method == method, ]
      rows <- rows[order(rows$alpha), ]
      polygon(
        c(rows$alpha, rev(rows$alpha)),
        c(pmax(0, rows$fdr_ci_lower), rev(pmin(1, rows$fdr_ci_upper))),
        col = adjustcolor(method_colors[[method]], alpha.f = 0.10),
        border = NA
      )
      lines(
        rows$alpha,
        rows$mean_fdr,
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
    legend = c(display_labels, "Nominal alpha"),
    col = c(unname(method_colors[methods]), "grey35"),
    lty = c(unname(method_lty[methods]), 3),
    lwd = c(rep(2.2, length(methods)), 1.4),
    bty = "n",
    cex = 0.80
  )
  mtext(
    "Paired support-width sensitivity: global empirical FDR",
    outer = TRUE,
    cex = 1.15
  )
  dev.off()
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

seeds <- c(12345L, 22345L, 32345L, 42345L, 52345L)
spiky_dfs <- c(12L, 14L, 16L)
methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
alpha_grid <- seq(0.005, 0.20, by = 0.005)

input_dirs <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  paste0(
    "confirm_spiky_paired_width_df12_14_16_seed",
    seeds,
    "_B100"
  )
)
configuration_paths <- file.path(input_dirs, "configuration.rds")
if (any(!file.exists(configuration_paths))) {
  stop("One or more width-confirmation configurations are missing.")
}
configurations <- lapply(configuration_paths, readRDS)
valid_configuration <- vapply(configurations, function(x) {
  identical(as.integer(x$spiky_dfs), spiky_dfs) &&
    identical(x$secondary_selection, "nearest") &&
    isTRUE(all.equal(x$secondary_fraction, c(0.55, 0.70))) &&
    isTRUE(x$pair_truth_rng_across_widths) &&
    identical(as.integer(x$efdr_permutations), 100L) &&
    identical(as.integer(x$J), 1000L) &&
    identical(as.integer(x$n_donors), 19L) &&
    identical(as.integer(x$n_covariates), 5L)
}, logical(1))
if (!all(valid_configuration)) {
  stop("One or more width-confirmation configurations are invalid.")
}

replicates <- list()
for (seed_index in seq_along(seeds)) {
  for (spiky_df in spiky_dfs) {
    path <- file.path(
      input_dirs[seed_index],
      "candidates",
      paste0("spiky_df", spiky_df, ".rds")
    )
    if (!file.exists(path)) stop("Missing confirmation cache: ", path)
    key <- paste(seeds[seed_index], spiky_df, sep = "_")
    replicates[[key]] <- readRDS(path)
    if (!identical(as.integer(replicates[[key]]$spiky_df), spiky_df)) {
      stop("A confirmation cache has the wrong spiky df: ", path)
    }
  }
}

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "confirm_spiky_paired_width_df12_14_16_B100_5seed"
)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

replicate_summary <- do.call(rbind, lapply(names(replicates), function(key) {
  key_parts <- strsplit(key, "_", fixed = TRUE)[[1]]
  row <- replicates[[key]]$summary
  row$seed <- as.integer(key_parts[1])
  row
}))
rownames(replicate_summary) <- NULL

mc_width_summary <- summarize_grouped_values(
  data = replicate_summary,
  group_names = "spiky_df",
  value_names = c(
    "mean_support_25",
    "mean_support_50",
    "mean_second_difference_roughness",
    "mean_quadratic_projection",
    "median_log_bayes_factor",
    "raw_spiky_power",
    "bf_spiky_power",
    "direct_quadratic_spiky_power",
    "raw_pi0",
    "bf_pi0"
  )
)

all_global_alpha <- do.call(rbind, lapply(names(replicates), function(key) {
  key_parts <- strsplit(key, "_", fixed = TRUE)[[1]]
  out <- replicates[[key]]$alpha_curve
  out$seed <- as.integer(key_parts[1])
  out$spiky_df <- as.integer(key_parts[2])
  out
}))
rownames(all_global_alpha) <- NULL

mc_global_alpha <- do.call(rbind, lapply(
  split(
    all_global_alpha,
    list(
      all_global_alpha$spiky_df,
      all_global_alpha$method,
      all_global_alpha$alpha
    ),
    drop = TRUE
  ),
  function(x) {
    power_summary <- summarize_mc_values(x$power)
    fdr_summary <- summarize_mc_values(x$empirical_fdr)
    data.frame(
      spiky_df = x$spiky_df[1],
      method = x$method[1],
      alpha = x$alpha[1],
      n_replications = length(unique(x$seed)),
      mean_power = power_summary[["mean"]],
      power_ci_lower = pmax(0, power_summary[["lower"]]),
      power_ci_upper = pmin(1, power_summary[["upper"]]),
      mean_fdr = fdr_summary[["mean"]],
      fdr_ci_lower = pmax(0, fdr_summary[["lower"]]),
      fdr_ci_upper = pmin(1, fdr_summary[["upper"]]),
      stringsAsFactors = FALSE
    )
  }
))
rownames(mc_global_alpha) <- NULL

all_spiky_alpha <- do.call(rbind, lapply(names(replicates), function(key) {
  key_parts <- strsplit(key, "_", fixed = TRUE)[[1]]
  seed <- as.integer(key_parts[1])
  spiky_df <- as.integer(key_parts[2])
  replicate <- replicates[[key]]
  unit_info <- replicate$effect_sim$unit_info
  spiky_units <- unit_info$unit_index[
    unit_info$effect_class == "dynamic_bspline" &
      unit_info$shape_profile == "spiky"
  ]
  do.call(rbind, lapply(methods, function(method) {
    method_result <- replicate$result_table[
      replicate$result_table$method == method,
    ]
    do.call(rbind, lapply(alpha_grid, function(alpha) {
      selected_units <- method_result$unit_index[
        is.finite(method_result$adjusted_score) &
          method_result$adjusted_score <= alpha
      ]
      true_positives <- sum(spiky_units %in% selected_units)
      data.frame(
        seed = seed,
        spiky_df = spiky_df,
        method = method,
        alpha = alpha,
        n_dynamic = length(spiky_units),
        true_positives = true_positives,
        power = true_positives / length(spiky_units),
        stringsAsFactors = FALSE
      )
    }))
  }))
}))
rownames(all_spiky_alpha) <- NULL

mc_spiky_alpha <- do.call(rbind, lapply(
  split(
    all_spiky_alpha,
    list(
      all_spiky_alpha$spiky_df,
      all_spiky_alpha$method,
      all_spiky_alpha$alpha
    ),
    drop = TRUE
  ),
  function(x) {
    power_summary <- summarize_mc_values(x$power)
    data.frame(
      spiky_df = x$spiky_df[1],
      method = x$method[1],
      alpha = x$alpha[1],
      n_dynamic = x$n_dynamic[1],
      n_replications = length(unique(x$seed)),
      mean_true_positives = mean(x$true_positives),
      mean_power = power_summary[["mean"]],
      power_ci_lower = pmax(0, power_summary[["lower"]]),
      power_ci_upper = pmin(1, power_summary[["upper"]]),
      stringsAsFactors = FALSE
    )
  }
))
rownames(mc_spiky_alpha) <- NULL

spiky_alpha_005 <- all_spiky_alpha[
  abs(all_spiky_alpha$alpha - 0.05) < 1e-12,
]
wide_alpha_005 <- reshape(
  spiky_alpha_005[, c("seed", "spiky_df", "method", "power")],
  idvar = c("seed", "spiky_df"),
  timevar = "method",
  direction = "wide"
)
names(wide_alpha_005) <- sub("^power\\.", "", names(wide_alpha_005))
wide_alpha_005$bf_minus_direct_quadratic <-
  wide_alpha_005[["FASH-IWP1-BF"]] -
  wide_alpha_005[["Direct-quadratic-LRT-eFDR-true-pi0"]]

paired_advantage_summary <- summarize_grouped_values(
  data = wide_alpha_005,
  group_names = "spiky_df",
  value_names = "bf_minus_direct_quadratic"
)

write_csv(
  replicate_summary,
  file.path(output_dir, "replicate_width_summary_alpha005.csv")
)
write_csv(
  mc_width_summary,
  file.path(output_dir, "mc_width_summary_alpha005.csv")
)
write_csv(
  all_global_alpha,
  file.path(output_dir, "all_replicate_global_alpha.csv")
)
write_csv(
  mc_global_alpha,
  file.path(output_dir, "mc_global_alpha_curve.csv")
)
write_csv(
  all_spiky_alpha,
  file.path(output_dir, "all_replicate_spiky_alpha.csv")
)
write_csv(
  mc_spiky_alpha,
  file.path(output_dir, "mc_spiky_alpha_curve.csv")
)
write_csv(
  wide_alpha_005,
  file.path(output_dir, "paired_spiky_power_alpha005.csv")
)
write_csv(
  paired_advantage_summary,
  file.path(output_dir, "paired_advantage_summary_alpha005.csv")
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

plot_width_power(
  mc_curve = mc_spiky_alpha,
  file = file.path(figure_dir, "spiky_power_across_alpha_by_width.png"),
  methods = methods,
  method_colors = method_colors,
  method_lty = method_lty
)
plot_width_fdr(
  mc_curve = mc_global_alpha,
  file = file.path(figure_dir, "global_fdr_across_alpha_by_width.png"),
  methods = methods,
  method_colors = method_colors,
  method_lty = method_lty
)

png(
  file.path(figure_dir, "geometry_and_power_alpha005.png"),
  width = 1500,
  height = 1200,
  res = 180
)
par(mfrow = c(2, 2), mar = c(4.2, 4.5, 2.7, 1.0), oma = c(0, 0, 2.4, 0))
plot(
  mc_width_summary$spiky_df,
  mc_width_summary$mean_mean_support_25,
  type = "b",
  pch = 19,
  lwd = 2,
  col = "#377EB8",
  xlab = "Local B-spline df",
  ylab = "Mean 25%-peak support",
  main = "Temporal support"
)
plot(
  mc_width_summary$spiky_df,
  mc_width_summary$mean_mean_second_difference_roughness,
  type = "b",
  pch = 19,
  lwd = 2,
  col = "#E41A1C",
  xlab = "Local B-spline df",
  ylab = "Mean second-difference roughness",
  main = "Curve roughness"
)
plot(
  mc_width_summary$spiky_df,
  mc_width_summary$mean_median_log_bayes_factor,
  type = "b",
  pch = 19,
  lwd = 2,
  col = "#4DAF4A",
  xlab = "Local B-spline df",
  ylab = "Median log Bayes factor",
  main = "IWP evidence"
)
plot(
  NA,
  xlim = range(mc_width_summary$spiky_df),
  ylim = c(0, 1),
  xlab = "Local B-spline df",
  ylab = "Mean spiky-effect power",
  main = "Power at alpha = 0.05"
)
for (method in methods) {
  value_name <- switch(
    method,
    "FASH-IWP1-Raw" = "mean_raw_spiky_power",
    "FASH-IWP1-BF" = "mean_bf_spiky_power",
    "Direct-quadratic-LRT-eFDR-true-pi0" =
      "mean_direct_quadratic_spiky_power"
  )
  lines(
    mc_width_summary$spiky_df,
    mc_width_summary[[value_name]],
    type = "b",
    pch = 19,
    lwd = 2,
    col = method_colors[[method]],
    lty = method_lty[[method]]
  )
}
legend(
  "bottomleft",
  legend = methods,
  col = unname(method_colors[methods]),
  lty = unname(method_lty[methods]),
  lwd = 2,
  pch = 19,
  bty = "n",
  cex = 0.72
)
mtext(
  "Paired support-width sensitivity at alpha = 0.05",
  outer = TRUE,
  cex = 1.15
)
dev.off()

example_seed <- seeds[1]
example_replicates <- lapply(spiky_dfs, function(spiky_df) {
  replicates[[paste(example_seed, spiky_df, sep = "_")]]$effect_sim
})
example_info <- example_replicates[[1]]$unit_info
example_index <- unlist(lapply(
  c("single", "same-sign double", "opposite-sign double"),
  function(pattern) {
    candidates <- which(
      example_info$effect_class == "dynamic_bspline" &
        example_info$shape_profile == "spiky" &
        example_info$spike_pattern == pattern
    )
    candidates[1]
  }
))
png(
  file.path(figure_dir, "paired_truth_examples_by_width.png"),
  width = 1500,
  height = 1100,
  res = 180
)
par(mfrow = c(length(example_index), length(spiky_dfs)),
    mar = c(3.2, 3.5, 2.5, 0.8), oma = c(0, 0, 2.5, 0))
for (row_index in seq_along(example_index)) {
  j <- example_index[row_index]
  row_values <- unlist(lapply(
    example_replicates,
    function(x) x$beta_matrix[j, ]
  ))
  row_ylim <- range(row_values)
  for (df_index in seq_along(spiky_dfs)) {
    plot(
      0:15,
      example_replicates[[df_index]]$beta_matrix[j, ],
      type = "l",
      lwd = 2,
      col = "#D95F02",
      ylim = row_ylim,
      xlab = "Time",
      ylab = if (df_index == 1) "True effect" else "",
      main = paste0(
        example_info$spike_pattern[j],
        "; df=",
        spiky_dfs[df_index]
      ),
      cex.main = 0.78
    )
    abline(h = 0, col = "grey70", lty = 3)
  }
}
mtext(
  "Paired localized truth examples from wider to sharper support",
  outer = TRUE,
  cex = 1.1
)
dev.off()

saveRDS(
  list(
    seeds = seeds,
    spiky_dfs = spiky_dfs,
    configurations = configurations,
    replicate_summary = replicate_summary,
    mc_width_summary = mc_width_summary,
    all_global_alpha = all_global_alpha,
    mc_global_alpha = mc_global_alpha,
    all_spiky_alpha = all_spiky_alpha,
    mc_spiky_alpha = mc_spiky_alpha,
    paired_spiky_power_alpha005 = wide_alpha_005,
    paired_advantage_summary_alpha005 = paired_advantage_summary
  ),
  file.path(output_dir, "confirmation_summary.rds")
)

message("Saved paired spiky-width confirmation summary to: ", output_dir)
