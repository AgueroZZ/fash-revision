#!/usr/bin/env Rscript

# Aggregate the matched internal comparison of the old local B-spline truth
# and two raised-cosine replacements.

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

summarize_grouped_values <- function(data, group_columns, value_columns) {
  key <- interaction(data[group_columns], drop = TRUE, lex.order = TRUE)
  rows <- lapply(split(data, key), function(x) {
    out <- x[1, group_columns, drop = FALSE]
    out$n_replications <- length(unique(x$seed))
    for (value in value_columns) {
      summary <- summarize_mc_values(x[[value]])
      out[[paste0("mean_", value)]] <- summary[["mean"]]
      out[[paste0(value, "_sd")]] <- summary[["sd"]]
      out[[paste0(value, "_mc_se")]] <- summary[["se"]]
      out[[paste0(value, "_ci_lower")]] <- summary[["lower"]]
      out[[paste0(value, "_ci_upper")]] <- summary[["upper"]]
    }
    out
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

truth_labels <- c(
  old_local_bspline_single = "Old local B-spline: one peak",
  raised_cosine_single = "Raised cosine: one peak",
  raised_cosine_multipeak = "Raised cosine: 1/2/3 peaks"
)
truth_axis_labels <- c(
  old_local_bspline_single = "Old B-spline",
  raised_cosine_single = "Raised single",
  raised_cosine_multipeak = "Raised multi"
)

compute_geometry <- function(beta_matrix, time_grid, truth_mechanism, seed) {
  centered <- beta_matrix - rowMeans(beta_matrix)
  maximum <- apply(abs(centered), 1, max)
  linear_basis <- scaled_time_polynomial(time_grid, degree = 1)
  quadratic_basis <- scaled_time_polynomial(time_grid, degree = 2)
  projection_fraction <- function(x, basis) {
    fitted <- qr.fitted(qr(basis), x)
    sum(fitted^2) / sum(x^2)
  }
  data.frame(
    truth_mechanism = truth_mechanism,
    seed = seed,
    n_dynamic = nrow(beta_matrix),
    centered_rms_mean = mean(sqrt(rowMeans(centered^2))),
    support_25_mean = mean(rowSums(abs(centered) >= 0.25 * maximum)),
    support_50_mean = mean(rowSums(abs(centered) >= 0.50 * maximum)),
    second_difference_roughness_mean = mean(apply(
      centered,
      1,
      function(x) sum(diff(x, differences = 2)^2) / sum(x^2)
    )),
    linear_projection_mean = mean(apply(
      centered,
      1,
      projection_fraction,
      basis = linear_basis
    )),
    quadratic_projection_mean = mean(apply(
      centered,
      1,
      projection_fraction,
      basis = quadratic_basis
    )),
    stringsAsFactors = FALSE
  )
}

plot_mechanism_alpha <- function(mc_curve, metric, file) {
  mechanisms <- names(truth_labels)
  styles <- revision_method_styles(unique(mc_curve$method))
  value_column <- if (metric == "power") "mean_power" else "mean_fdr"
  lower_column <- if (metric == "power") "power_ci_lower" else "fdr_ci_lower"
  upper_column <- if (metric == "power") "power_ci_upper" else "fdr_ci_upper"
  png(file, width = 2100, height = 720, res = 180)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  })
  par(mfrow = c(1, 3), mar = c(4.3, 4.3, 4.1, 1.0))
  for (mechanism in mechanisms) {
    panel <- mc_curve[mc_curve$scenario == mechanism, ]
    y_limits <- if (metric == "power") {
      c(0, 1)
    } else {
      c(0, max(0.21, panel[[upper_column]], na.rm = TRUE))
    }
    plot(
      NA,
      xlim = range(panel$alpha),
      ylim = y_limits,
      xlab = "Nominal FDR level alpha",
      ylab = if (metric == "power") "Mean power" else "Monte Carlo E[FDP]",
      main = truth_labels[mechanism],
      cex.main = 0.90
    )
    grid(col = "gray90")
    abline(v = 0.05, col = "gray55", lty = 3)
    if (metric == "fdr") abline(a = 0, b = 1, col = "gray35", lty = 3)
    for (i in seq_along(styles$methods)) {
      method <- styles$methods[i]
      curve <- panel[panel$method == method, ]
      curve <- curve[order(curve$alpha), ]
      lower <- pmax(0, curve[[lower_column]])
      upper <- pmin(y_limits[2], curve[[upper_column]])
      polygon(
        c(curve$alpha, rev(curve$alpha)),
        c(lower, rev(upper)),
        col = grDevices::adjustcolor(styles$col[i], alpha.f = 0.12),
        border = NA
      )
      lines(
        curve$alpha,
        curve[[value_column]],
        col = styles$col[i],
        lty = styles$lty[i],
        lwd = 2.3
      )
    }
    if (mechanism == mechanisms[1]) {
      legend(
        "bottomright",
        legend = c("FASH-BF", "Linear LRT", "Quadratic LRT"),
        col = styles$col,
        lty = styles$lty,
        lwd = 2.3,
        bty = "n",
        cex = 0.76
      )
    }
  }
}

plot_advantages <- function(summary, file) {
  mechanisms <- names(truth_labels)
  comparators <- c(
    "Direct-linear-LRT-eFDR-true-pi0",
    "Direct-quadratic-LRT-eFDR-true-pi0"
  )
  colors <- c("#4D4D4D", "#984EA3")
  png(file, width = 1500, height = 820, res = 180)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  })
  par(mar = c(6.5, 4.8, 4.0, 1.2))
  x <- seq_along(mechanisms)
  plot(
    NA,
    xlim = c(0.6, length(mechanisms) + 0.4),
    ylim = range(
      summary$power_advantage_ci_lower,
      summary$power_advantage_ci_upper,
      0
    ),
    xaxt = "n",
    xlab = "",
    ylab = "Power advantage at alpha = 0.05",
    main = "FASH-IWP1-BF minus direct interaction"
  )
  abline(h = 0, col = "gray40", lty = 3)
  grid(nx = NA, ny = NULL, col = "gray90")
  offsets <- c(-0.08, 0.08)
  for (i in seq_along(comparators)) {
    rows <- summary[match(
      paste(mechanisms, comparators[i]),
      paste(summary$truth_mechanism, summary$comparator)
    ), ]
    segments(
      x + offsets[i],
      rows$power_advantage_ci_lower,
      x + offsets[i],
      rows$power_advantage_ci_upper,
      col = colors[i],
      lwd = 2
    )
    points(
      x + offsets[i],
      rows$mean_power_advantage,
      pch = 15 + i,
      col = colors[i],
      cex = 1.25
    )
  }
  axis(1, at = x, labels = unname(truth_axis_labels), las = 2)
  legend(
    "topleft",
    legend = c("Versus linear LRT", "Versus quadratic LRT"),
    col = colors,
    pch = 16:17,
    bty = "n"
  )
}

plot_geometry <- function(geometry_summary, file) {
  mechanisms <- names(truth_labels)
  metrics <- c(
    "centered_rms_mean",
    "support_50_mean",
    "second_difference_roughness_mean",
    "quadratic_projection_mean"
  )
  labels <- c(
    "Centered RMS",
    "Support above 50% maximum",
    "Second-difference roughness",
    "Quadratic projection fraction"
  )
  colors <- c("#D95F02", "#1B9E77", "#7570B3")
  png(file, width = 2000, height = 650, res = 180)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  })
  par(mfrow = c(1, 4), mar = c(8.0, 4.3, 3.5, 0.8))
  for (i in seq_along(metrics)) {
    rows <- geometry_summary[match(
      mechanisms,
      geometry_summary$truth_mechanism
    ), ]
    values <- rows[[paste0("mean_", metrics[i])]]
    lower <- rows[[paste0(metrics[i], "_ci_lower")]]
    upper <- rows[[paste0(metrics[i], "_ci_upper")]]
    plot(
      seq_along(mechanisms),
      values,
      type = "n",
      xaxt = "n",
      xlab = "",
      ylab = labels[i],
      main = labels[i],
      ylim = range(lower, upper)
    )
    grid(nx = NA, ny = NULL, col = "gray90")
    segments(seq_along(mechanisms), lower, seq_along(mechanisms), upper, lwd = 2)
    points(
      seq_along(mechanisms),
      values,
      pch = 16,
      cex = 1.25,
      col = colors
    )
    axis(
      1,
      at = seq_along(mechanisms),
      labels = unname(truth_axis_labels),
      las = 2
    )
  }
}

plot_truth_examples <- function(old_effect, old_latent, single, multi, file) {
  old_dynamic <- which(
    old_effect$unit_info$effect_class == "dynamic_local_bspline_transient"
  )
  old_order <- order(old_latent$transient_locations)
  old_positions <- old_order[round(seq(1, length(old_order), length.out = 5))]
  old_selected <- old_dynamic[old_positions]
  single_dynamic <- which(single$unit_info$effect_class == "dynamic_bspline")
  single_selected <- single_dynamic[seq_len(5)]
  multi_dynamic <- which(multi$unit_info$effect_class == "dynamic_bspline")
  multi_cells <- unique(multi$unit_info$cell_id[multi_dynamic])
  multi_selected <- vapply(multi_cells, function(cell) {
    multi_dynamic[which(multi$unit_info$cell_id[multi_dynamic] == cell)[1]]
  }, integer(1))

  png(file, width = 2100, height = 1250, res = 180)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  })
  par(mfrow = c(3, 5), mar = c(2.8, 3.0, 2.8, 0.6), oma = c(0, 3.0, 0, 0))
  plot_one <- function(x, y, x_observed, y_observed, title) {
    plot(
      x,
      y,
      type = "l",
      lwd = 2.2,
      col = "#D95F02",
      xlab = "Time",
      ylab = "",
      main = title,
      cex.main = 0.76
    )
    abline(h = 0, col = "gray70", lty = 3)
    points(x_observed, y_observed, pch = 16, cex = 0.5, col = "#1F78B4")
  }
  for (index in old_selected) {
    y <- old_effect$beta_matrix[index, ]
    plot_one(0:15, y, 0:15, y, "Old local B-spline")
  }
  for (index in single_selected) {
    plot_one(
      single$evaluation_grid,
      single$beta_evaluation[index, ],
      0:15,
      single$beta_matrix[index, ],
      "Raised cosine: one peak"
    )
  }
  for (index in multi_selected) {
    plot_one(
      multi$evaluation_grid,
      multi$beta_evaluation[index, ],
      0:15,
      multi$beta_matrix[index, ],
      gsub("__", "; ", multi$unit_info$cell_id[index])
    )
  }
  mtext("True effect", side = 2, outer = TRUE, line = 1.0)
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

seeds <- c(12345L, 22345L, 32345L, 42345L, 52345L)
mechanisms <- names(truth_labels)
methods <- c(
  "FASH-IWP1-BF",
  "Direct-linear-LRT-eFDR-true-pi0",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
old_paths <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "mc",
  "spiky_transient_pilot5",
  "replicates",
  paste0("seed_", seeds, ".rds")
)
new_paths <- unlist(lapply(
  c("raised_cosine_single", "raised_cosine_multipeak"),
  function(mechanism) {
    file.path(
      workflowr_root,
      "output",
      "revision_simulations",
      "internal",
      paste0(
        "spiky_truth_replacement_",
        mechanism,
        "_seed",
        seeds,
        "_B100"
      ),
      "replacement_result.rds"
    )
  }
))
if (any(!file.exists(c(old_paths, new_paths)))) {
  stop("At least one required old or replacement cache is missing.")
}

old_replicates <- lapply(old_paths, readRDS)
single_replicates <- lapply(new_paths[seq_along(seeds)], readRDS)
multi_replicates <- lapply(new_paths[length(seeds) + seq_along(seeds)], readRDS)

collect_alpha <- function(replicates, mechanism, old = FALSE) {
  do.call(rbind, lapply(seq_along(replicates), function(i) {
    x <- replicates[[i]]$alpha_curve
    x <- x[x$method %in% methods, ]
    x$seed <- seeds[i]
    x$scenario <- mechanism
    x$truth_mechanism <- mechanism
    x
  }))
}
all_alpha <- rbind(
  collect_alpha(old_replicates, "old_local_bspline_single", old = TRUE),
  collect_alpha(single_replicates, "raised_cosine_single"),
  collect_alpha(multi_replicates, "raised_cosine_multipeak")
)
mc_alpha <- summarize_mc_alpha_curves(all_alpha)

alpha_005 <- all_alpha[abs(all_alpha$alpha - 0.05) < 1e-12, ]
mc_alpha_005 <- mc_alpha[abs(mc_alpha$alpha - 0.05) < 1e-12, ]
alpha_wide <- reshape(
  alpha_005[, c("seed", "truth_mechanism", "method", "power")],
  idvar = c("seed", "truth_mechanism"),
  timevar = "method",
  direction = "wide"
)
comparators <- c(
  "Direct-linear-LRT-eFDR-true-pi0",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
advantages <- do.call(rbind, lapply(comparators, function(comparator) {
  data.frame(
    seed = alpha_wide$seed,
    truth_mechanism = alpha_wide$truth_mechanism,
    comparator = comparator,
    power_advantage =
      alpha_wide[["power.FASH-IWP1-BF"]] -
      alpha_wide[[paste0("power.", comparator)]],
    stringsAsFactors = FALSE
  )
}))
advantage_summary <- summarize_grouped_values(
  advantages,
  c("truth_mechanism", "comparator"),
  "power_advantage"
)

mechanism_wide <- reshape(
  alpha_005[, c("seed", "method", "truth_mechanism", "power")],
  idvar = c("seed", "method"),
  timevar = "truth_mechanism",
  direction = "wide"
)
replacement_changes <- do.call(rbind, lapply(
  c("raised_cosine_single", "raised_cosine_multipeak"),
  function(mechanism) {
    data.frame(
      seed = mechanism_wide$seed,
      method = mechanism_wide$method,
      replacement = mechanism,
      power_change =
        mechanism_wide[[paste0("power.", mechanism)]] -
        mechanism_wide[["power.old_local_bspline_single"]],
      stringsAsFactors = FALSE
    )
  }
))
replacement_change_summary <- summarize_grouped_values(
  replacement_changes,
  c("method", "replacement"),
  "power_change"
)

geometry_rows <- list()
for (i in seq_along(seeds)) {
  old_paired <- simulate_paired_local_bspline_effect_sets(
    candidate_settings = data.frame(
      setting_id = "old",
      transient_bspline_df = 16L,
      dynamic_amplitude = 2.75,
      normalization = "center_then_scale",
      stringsAsFactors = FALSE
    ),
    n_variants = 1000,
    class_probs = c(
      dynamic_local_bspline_transient = 0.20,
      constant = 0.40,
      zero = 0.40
    ),
    transient_bspline_degree = 3L,
    seed = seeds[i]
  )
  old_effect <- old_paired$effect_sets$old
  old_dynamic <- old_effect$unit_info$effect_class ==
    "dynamic_local_bspline_transient"
  old_geometry <- compute_geometry(
    old_effect$beta_matrix[old_dynamic, , drop = FALSE],
    make_time_grid(),
    "old_local_bspline_single",
    seeds[i]
  )
  geometry_rows[[length(geometry_rows) + 1L]] <- old_geometry
  geometry_rows[[length(geometry_rows) + 1L]] <-
    single_replicates[[i]]$geometry[, names(old_geometry)]
  geometry_rows[[length(geometry_rows) + 1L]] <-
    multi_replicates[[i]]$geometry[, names(old_geometry)]
}
all_geometry <- do.call(rbind, geometry_rows)
geometry_metrics <- c(
  "centered_rms_mean", "support_25_mean", "support_50_mean",
  "second_difference_roughness_mean", "linear_projection_mean",
  "quadratic_projection_mean"
)
geometry_summary <- summarize_grouped_values(
  all_geometry,
  "truth_mechanism",
  geometry_metrics
)

pi0_rows <- do.call(rbind, lapply(seq_along(seeds), function(i) {
  old_pi0 <- old_replicates[[i]]$pi0
  old_value <- old_pi0$estimated_pi0[
    old_pi0$method == "FASH-IWP1" & old_pi0$fit == "BF-corrected"
  ]
  data.frame(
    seed = seeds[i],
    truth_mechanism = mechanisms,
    estimated_pi0 = c(
      old_value,
      single_replicates[[i]]$bf_pi0,
      multi_replicates[[i]]$bf_pi0
    ),
    stringsAsFactors = FALSE
  )
}))
pi0_summary <- summarize_grouped_values(
  pi0_rows,
  "truth_mechanism",
  "estimated_pi0"
)

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "spiky_truth_replacement_B100_5seed"
)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(all_alpha, file.path(output_dir, "all_replicate_alpha.csv"))
write_csv(mc_alpha, file.path(output_dir, "mc_alpha.csv"))
write_csv(alpha_005, file.path(output_dir, "all_replicate_alpha005.csv"))
write_csv(mc_alpha_005, file.path(output_dir, "mc_alpha005.csv"))
write_csv(advantages, file.path(output_dir, "paired_method_advantages_alpha005.csv"))
write_csv(
  advantage_summary,
  file.path(output_dir, "mc_method_advantages_alpha005.csv")
)
write_csv(
  replacement_changes,
  file.path(output_dir, "paired_replacement_power_changes_alpha005.csv")
)
write_csv(
  replacement_change_summary,
  file.path(output_dir, "mc_replacement_power_changes_alpha005.csv")
)
write_csv(all_geometry, file.path(output_dir, "all_replicate_geometry.csv"))
write_csv(geometry_summary, file.path(output_dir, "geometry_summary.csv"))
write_csv(pi0_rows, file.path(output_dir, "all_replicate_bf_pi0.csv"))
write_csv(pi0_summary, file.path(output_dir, "bf_pi0_summary.csv"))

plot_mechanism_alpha(
  mc_alpha,
  "power",
  file.path(figure_dir, "power_across_alpha.png")
)
plot_mechanism_alpha(
  mc_alpha,
  "fdr",
  file.path(figure_dir, "fdr_across_alpha.png")
)
plot_advantages(
  advantage_summary,
  file.path(figure_dir, "method_power_advantages_alpha005.png")
)
plot_geometry(
  geometry_summary,
  file.path(figure_dir, "truth_geometry_comparison.png")
)

old_seed_one <- simulate_paired_local_bspline_effect_sets(
  candidate_settings = data.frame(
    setting_id = "old",
    transient_bspline_df = 16L,
    dynamic_amplitude = 2.75,
    normalization = "center_then_scale",
    stringsAsFactors = FALSE
  ),
  n_variants = 1000,
  class_probs = c(
    dynamic_local_bspline_transient = 0.20,
    constant = 0.40,
    zero = 0.40
  ),
  transient_bspline_degree = 3L,
  seed = seeds[1]
)
plot_truth_examples(
  old_seed_one$effect_sets$old,
  old_seed_one$latent,
  single_replicates[[1]]$effect_sim,
  multi_replicates[[1]]$effect_sim,
  file.path(figure_dir, "truth_mechanism_examples.png")
)

saveRDS(
  list(
    seeds = seeds,
    truth_labels = truth_labels,
    all_alpha = all_alpha,
    mc_alpha = mc_alpha,
    advantages = advantages,
    advantage_summary = advantage_summary,
    replacement_changes = replacement_changes,
    replacement_change_summary = replacement_change_summary,
    all_geometry = all_geometry,
    geometry_summary = geometry_summary,
    pi0_rows = pi0_rows,
    pi0_summary = pi0_summary
  ),
  file.path(output_dir, "comparison_summary.rds")
)

message("Saved matched spiky truth replacement comparison to: ", output_dir)
