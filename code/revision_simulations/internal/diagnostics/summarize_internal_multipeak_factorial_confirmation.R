#!/usr/bin/env Rscript

# Aggregate the locked five-seed raised-cosine multi-peak experiment.
# This internal analysis does not modify any workflowr page.

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

plot_power_heatmaps <- function(summary_005, file) {
  width_order <- c("not_spiky", "mildly_spiky", "spiky", "very_spiky")
  width_labels <- c("Not\nspiky", "Mildly\nspiky", "Spiky", "Very\nspiky")
  methods <- revision_method_styles(unique(summary_005$method))$methods
  png(file, width = 2100, height = 760, res = 180)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  })
  par(mfrow = c(1, length(methods)), mar = c(6.2, 4.5, 4.3, 1.2))
  palette <- colorRampPalette(c("#F7FBFF", "#6BAED6", "#08306B"))(101)
  for (method in methods) {
    x <- summary_005[summary_005$method == method, ]
    matrix_values <- matrix(
      NA_real_,
      nrow = 3,
      ncol = length(width_order),
      dimnames = list(1:3, width_order)
    )
    for (i in seq_len(nrow(x))) {
      matrix_values[
        as.character(x$spike_count[i]),
        x$width_label[i]
      ] <- x$mean_power[i]
    }
    image(
      x = seq_along(width_order),
      y = 1:3,
      z = t(matrix_values),
      zlim = c(0, 1),
      col = palette,
      axes = FALSE,
      xlab = "",
      ylab = "Number of peaks",
      main = method
    )
    axis(1, at = seq_along(width_order), labels = width_labels)
    axis(2, at = 1:3, labels = 1:3, las = 1)
    box()
    for (i in 1:3) {
      for (j in seq_along(width_order)) {
        value <- matrix_values[i, j]
        text(
          j,
          i,
          sprintf("%.2f", value),
          col = if (is.finite(value) && value >= 0.58) "white" else "black",
          font = 2
        )
      }
    }
  }
  mtext(
    "Mean power at alpha = 0.05 across five seeds",
    outer = TRUE,
    line = -1.2,
    cex = 1.15,
    font = 2
  )
}

plot_advantage_heatmaps <- function(advantage_summary, file) {
  width_order <- c("not_spiky", "mildly_spiky", "spiky", "very_spiky")
  width_labels <- c("Not\nspiky", "Mildly\nspiky", "Spiky", "Very\nspiky")
  comparators <- c(
    "Direct-linear-LRT-eFDR-true-pi0",
    "Direct-quadratic-LRT-eFDR-true-pi0"
  )
  png(file, width = 1800, height = 760, res = 180)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  })
  par(mfrow = c(1, 2), mar = c(6.2, 4.5, 4.3, 1.2))
  palette <- colorRampPalette(c("#B2182B", "#F7F7F7", "#2166AC"))(101)
  limit <- max(abs(c(
    advantage_summary$mean_power_advantage,
    advantage_summary$power_advantage_ci_lower,
    advantage_summary$power_advantage_ci_upper
  )), na.rm = TRUE)
  limit <- max(0.05, limit)
  for (comparator in comparators) {
    x <- advantage_summary[advantage_summary$comparator == comparator, ]
    matrix_values <- matrix(
      NA_real_,
      nrow = 3,
      ncol = length(width_order),
      dimnames = list(1:3, width_order)
    )
    for (i in seq_len(nrow(x))) {
      matrix_values[
        as.character(x$spike_count[i]),
        x$width_label[i]
      ] <- x$mean_power_advantage[i]
    }
    image(
      x = seq_along(width_order),
      y = 1:3,
      z = t(matrix_values),
      zlim = c(-limit, limit),
      col = palette,
      axes = FALSE,
      xlab = "",
      ylab = "Number of peaks",
      main = paste("FASH-IWP1-BF minus", sub("-LRT.*", "", comparator)),
      cex.main = 0.90
    )
    axis(1, at = seq_along(width_order), labels = width_labels)
    axis(2, at = 1:3, labels = 1:3, las = 1)
    box()
    for (i in 1:3) {
      for (j in seq_along(width_order)) {
        text(j, i, sprintf("%+.2f", matrix_values[i, j]), font = 2)
      }
    }
  }
  mtext(
    "Paired mean power advantage at alpha = 0.05",
    outer = TRUE,
    line = -1.2,
    cex = 1.15,
    font = 2
  )
}

plot_cell_alpha_grid <- function(cell_summary, file) {
  width_order <- c("not_spiky", "mildly_spiky", "spiky", "very_spiky")
  styles <- revision_method_styles(unique(cell_summary$method))
  png(file, width = 2200, height = 1700, res = 180)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  })
  par(
    mfrow = c(3, 4),
    mar = c(3.4, 3.6, 2.7, 0.8),
    oma = c(2.2, 2.2, 3.5, 0)
  )
  for (spike_count in 1:3) {
    for (width_label in width_order) {
      panel <- cell_summary[
        cell_summary$spike_count == spike_count &
          cell_summary$width_label == width_label,
      ]
      plot(
        NA,
        xlim = range(panel$alpha),
        ylim = c(0, 1),
        xlab = "",
        ylab = "",
        main = paste0(spike_count, " peak", if (spike_count > 1) "s" else "",
                      "; ", gsub("_", " ", width_label))
      )
      grid(col = "gray90")
      abline(v = 0.05, col = "gray55", lty = 3)
      for (i in seq_along(styles$methods)) {
        method <- styles$methods[i]
        curve <- panel[panel$method == method, ]
        curve <- curve[order(curve$alpha), ]
        lines(
          curve$alpha,
          curve$mean_power,
          col = styles$col[i],
          lty = styles$lty[i],
          lwd = 2.2
        )
      }
    }
  }
  legend(
    "bottom",
    inset = -0.10,
    xpd = NA,
    horiz = TRUE,
    legend = c("FASH-BF", "Linear LRT", "Quadratic LRT"),
    col = styles$col,
    lty = styles$lty,
    lwd = 2.2,
    bty = "n",
    cex = 0.85
  )
  mtext("Nominal FDR level alpha", side = 1, outer = TRUE, line = 0.4)
  mtext("Mean power", side = 2, outer = TRUE, line = 0.5)
  mtext(
    "Power across alpha by peak count and width",
    side = 3,
    outer = TRUE,
    line = 1.5,
    font = 2,
    cex = 1.25
  )
}

plot_geometry_trends <- function(geometry_summary, file) {
  width_order <- c("not_spiky", "mildly_spiky", "spiky", "very_spiky")
  width_labels <- c("Not spiky", "Mildly spiky", "Spiky", "Very spiky")
  metrics <- c(
    "support_50",
    "second_difference_roughness",
    "quadratic_projection"
  )
  metric_labels <- c(
    "Observed support above 50% of maximum",
    "Second-difference roughness",
    "Quadratic projection fraction"
  )
  colors <- c("#1B9E77", "#D95F02", "#7570B3")
  png(file, width = 1900, height = 650, res = 180)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  })
  par(mfrow = c(1, 3), mar = c(5.8, 4.6, 3.5, 1.0))
  for (m in seq_along(metrics)) {
    metric <- metrics[m]
    y_limits <- range(geometry_summary[[metric]], finite = TRUE)
    plot(
      NA,
      xlim = c(1, 4),
      ylim = y_limits,
      xaxt = "n",
      xlab = "",
      ylab = metric_labels[m],
      main = metric_labels[m]
    )
    axis(1, at = 1:4, labels = width_labels, las = 2)
    grid(col = "gray90")
    for (spike_count in 1:3) {
      x <- geometry_summary[geometry_summary$spike_count == spike_count, ]
      x$width_rank <- match(x$width_label, width_order)
      x <- x[order(x$width_rank), ]
      lines(
        x$width_rank,
        x[[metric]],
        type = "b",
        pch = 14 + spike_count,
        lwd = 2.2,
        col = colors[spike_count]
      )
    }
    if (m == 1) {
      legend(
        "topright",
        legend = paste0(1:3, " peak", ifelse(1:3 > 1, "s", "")),
        col = colors,
        pch = 15:17,
        lwd = 2.2,
        bty = "n"
      )
    }
  }
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

seeds <- c(12345L, 22345L, 32345L, 42345L, 52345L)
input_dirs <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  paste0("multipeak_factorial_seed", seeds, "_B100")
)
input_paths <- file.path(input_dirs, "factorial_result.rds")
if (any(!file.exists(input_paths))) {
  stop("Missing factorial result: ", paste(input_paths[!file.exists(input_paths)], collapse = ", "))
}
replicates <- lapply(input_paths, readRDS)

reference <- replicates[[1]]$configuration
fixed_fields <- c(
  "J", "n_donors", "n_covariates", "expression_noise_sd",
  "efdr_permutations", "class_probs", "true_pi0", "width_levels",
  "spike_counts", "relative_amplitude_range", "target_centered_rms"
)
for (field in fixed_fields) {
  if (any(vapply(
    replicates[-1],
    function(x) !isTRUE(all.equal(x$configuration[[field]], reference[[field]])),
    logical(1)
  ))) {
    stop("Replicates disagree on configuration field: ", field)
  }
}

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "multipeak_factorial_B100_5seed"
)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

all_global <- do.call(rbind, lapply(replicates, `[[`, "global_alpha"))
mc_global <- summarize_mc_alpha_curves(all_global)

all_cell <- do.call(rbind, lapply(replicates, `[[`, "cell_alpha"))
seed_cell <- aggregate(
  cbind(true_positives, power) ~
    seed + method + alpha + spike_count + width_label + width_half,
  data = all_cell,
  FUN = mean
)
mc_cell <- summarize_grouped_values(
  seed_cell,
  group_columns = c(
    "method", "alpha", "spike_count", "width_label", "width_half"
  ),
  value_columns = c("true_positives", "power")
)
mc_cell$power_ci_lower <- pmax(0, mc_cell$power_ci_lower)
mc_cell$power_ci_upper <- pmin(1, mc_cell$power_ci_upper)

mc_sign_cell <- summarize_grouped_values(
  all_cell,
  group_columns = c(
    "method", "alpha", "spike_count", "width_label", "width_half",
    "sign_pattern"
  ),
  value_columns = c("true_positives", "power")
)
mc_sign_cell$power_ci_lower <- pmax(0, mc_sign_cell$power_ci_lower)
mc_sign_cell$power_ci_upper <- pmin(1, mc_sign_cell$power_ci_upper)

cell_005 <- seed_cell[abs(seed_cell$alpha - 0.05) < 1e-12, ]
mc_cell_005 <- mc_cell[abs(mc_cell$alpha - 0.05) < 1e-12, ]
cell_wide <- reshape(
  cell_005[, c(
    "seed", "spike_count", "width_label", "width_half", "method", "power"
  )],
  idvar = c("seed", "spike_count", "width_label", "width_half"),
  timevar = "method",
  direction = "wide"
)
comparators <- c(
  "Direct-linear-LRT-eFDR-true-pi0",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
advantages <- do.call(rbind, lapply(comparators, function(comparator) {
  data.frame(
    seed = cell_wide$seed,
    spike_count = cell_wide$spike_count,
    width_label = cell_wide$width_label,
    width_half = cell_wide$width_half,
    comparator = comparator,
    power_advantage =
      cell_wide[["power.FASH-IWP1-BF"]] -
      cell_wide[[paste0("power.", comparator)]],
    stringsAsFactors = FALSE
  )
}))
advantage_summary <- summarize_grouped_values(
  advantages,
  group_columns = c(
    "spike_count", "width_label", "width_half", "comparator"
  ),
  value_columns = "power_advantage"
)

all_geometry <- do.call(rbind, lapply(seq_along(replicates), function(i) {
  x <- replicates[[i]]$geometry_summary
  x$seed <- seeds[i]
  x
}))
seed_geometry <- aggregate(
  cbind(
    support_25, support_50, second_difference_roughness,
    linear_projection, quadratic_projection, log_bayes_factor
  ) ~ seed + spike_count + width_label + width_half,
  data = all_geometry,
  FUN = mean
)
geometry_summary <- aggregate(
  cbind(
    support_25, support_50, second_difference_roughness,
    linear_projection, quadratic_projection, log_bayes_factor
  ) ~ spike_count + width_label + width_half,
  data = seed_geometry,
  FUN = mean
)

pi0_rows <- data.frame(
  seed = seeds,
  estimated_pi0 = vapply(replicates, `[[`, numeric(1), "bf_pi0")
)
pi0_summary <- summarize_mc_values(pi0_rows$estimated_pi0)
pi0_summary <- data.frame(
  n_replications = length(seeds),
  mean_estimated_pi0 = pi0_summary[["mean"]],
  estimated_pi0_sd = pi0_summary[["sd"]],
  estimated_pi0_mc_se = pi0_summary[["se"]],
  estimated_pi0_ci_lower = pi0_summary[["lower"]],
  estimated_pi0_ci_upper = pi0_summary[["upper"]]
)

write_csv(all_global, file.path(output_dir, "all_replicate_global_alpha.csv"))
write_csv(mc_global, file.path(output_dir, "mc_global_alpha.csv"))
write_csv(all_cell, file.path(output_dir, "all_replicate_sign_cell_alpha.csv"))
write_csv(seed_cell, file.path(output_dir, "all_replicate_cell_alpha.csv"))
write_csv(mc_cell, file.path(output_dir, "mc_cell_alpha.csv"))
write_csv(mc_cell_005, file.path(output_dir, "mc_cell_alpha005.csv"))
write_csv(mc_sign_cell, file.path(output_dir, "mc_sign_cell_alpha.csv"))
write_csv(advantages, file.path(output_dir, "paired_power_advantages_alpha005.csv"))
write_csv(
  advantage_summary,
  file.path(output_dir, "mc_power_advantage_alpha005.csv")
)
write_csv(all_geometry, file.path(output_dir, "all_replicate_geometry.csv"))
write_csv(geometry_summary, file.path(output_dir, "geometry_summary.csv"))
write_csv(pi0_rows, file.path(output_dir, "all_replicate_bf_pi0.csv"))
write_csv(pi0_summary, file.path(output_dir, "bf_pi0_summary.csv"))

plot_power_heatmaps(
  mc_cell_005,
  file.path(figure_dir, "power_heatmaps_alpha005.png")
)
plot_advantage_heatmaps(
  advantage_summary,
  file.path(figure_dir, "power_advantage_heatmaps_alpha005.png")
)
plot_cell_alpha_grid(
  mc_cell,
  file.path(figure_dir, "power_across_alpha_factorial_grid.png")
)
plot_mc_alpha_curves(
  mc_curve = mc_global,
  metric = "fdr",
  file = file.path(figure_dir, "global_fdr_across_alpha.png"),
  title = "Global empirical FDR across five seeds",
  subtitle = "Raised-cosine multi-peak effects; direct eFDR uses true pi0 and B=100",
  style_profile = "combined"
)
plot_geometry_trends(
  geometry_summary,
  file.path(figure_dir, "shape_geometry_trends.png")
)

file.copy(
  file.path(input_dirs[1], "figures", "truth_examples_same_sign.png"),
  file.path(figure_dir, "truth_examples_same_sign.png"),
  overwrite = TRUE
)
file.copy(
  file.path(input_dirs[1], "figures", "truth_examples_alternating_sign.png"),
  file.path(figure_dir, "truth_examples_alternating_sign.png"),
  overwrite = TRUE
)

saveRDS(
  list(
    configuration = reference,
    seeds = seeds,
    all_global = all_global,
    mc_global = mc_global,
    all_cell = all_cell,
    seed_cell = seed_cell,
    mc_cell = mc_cell,
    mc_sign_cell = mc_sign_cell,
    advantages = advantages,
    advantage_summary = advantage_summary,
    all_geometry = all_geometry,
    geometry_summary = geometry_summary,
    pi0_rows = pi0_rows,
    pi0_summary = pi0_summary
  ),
  file.path(output_dir, "confirmation_summary.rds")
)

message("Saved multi-peak factorial confirmation summary to: ", output_dir)
