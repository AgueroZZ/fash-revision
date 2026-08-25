#!/usr/bin/env Rscript

# Diagnose why refining the evaluation grid changes the middle-functional LFSR.
# The fitted FASH model is held fixed. Posterior draws are generated once on
# the 0.05-day grid and reused for the nested 0.10-day comparison.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

middle_components <- function(samples, evaluation_grid) {
  tolerance <- sqrt(.Machine$double.eps)
  inside <- evaluation_grid >= 4 - tolerance & evaluation_grid <= 11 + tolerance
  absolute_samples <- abs(samples)
  inside_max <- matrixStats::colMaxs(absolute_samples, rows = which(inside))
  outside_max <- matrixStats::colMaxs(absolute_samples, rows = which(!inside))
  statistic <- inside_max - outside_max
  list(
    inside_max = inside_max,
    outside_max = outside_max,
    statistic = statistic,
    error = statistic <= 0
  )
}

binomial_summary <- function(indicator) {
  estimate <- mean(indicator)
  n <- length(indicator)
  interval <- binom.test(sum(indicator), n)$conf.int
  c(
    estimate = estimate,
    mcse = sqrt(estimate * (1 - estimate) / n),
    ci_lower = interval[1],
    ci_upper = interval[2]
  )
}

summarize_trajectory <- function(samples, evaluation_grid) {
  data.frame(
    time = evaluation_grid,
    posterior_mean = rowMeans(samples),
    posterior_median = matrixStats::rowMedians(samples),
    lower = as.numeric(matrixStats::rowQuantiles(samples, probs = 0.025)),
    upper = as.numeric(matrixStats::rowQuantiles(samples, probs = 0.975)),
    stringsAsFactors = FALSE
  )
}

plot_pair_diagnostic <- function(
    pair_id,
    observed,
    trajectory,
    fine_grid,
    coarse_rows,
    fine_components,
    coarse_components,
    output_path,
    posterior_draws) {
  coarse_grid <- fine_grid[coarse_rows]
  coarse_trajectory <- trajectory[coarse_rows, , drop = FALSE]
  added_rows <- setdiff(seq_along(fine_grid), coarse_rows)

  absolute_samples <- abs(attr(trajectory, "posterior_samples"))
  outside_rows <- which(fine_grid < 4 | fine_grid > 11)
  outside_argmax <- outside_rows[max.col(t(absolute_samples[outside_rows, , drop = FALSE]), ties.method = "first")]
  added_outside_argmax <- outside_argmax[outside_argmax %in% added_rows]
  dominant_added_time <- if (length(added_outside_argmax) > 0L) {
    as.numeric(names(sort(table(fine_grid[added_outside_argmax]), decreasing = TRUE))[1])
  } else {
    3.95
  }
  zoom_limits <- pmax(c(dominant_added_time - 0.55, dominant_added_time + 0.55), 0)
  zoom_limits <- pmin(zoom_limits, 15)

  lfsr_fine <- binomial_summary(fine_components$error)
  lfsr_coarse <- binomial_summary(coarse_components$error)
  paired_difference <- as.numeric(fine_components$error) - as.numeric(coarse_components$error)
  paired_mcse <- sd(paired_difference) / sqrt(posterior_draws)

  png(output_path, width = 2400, height = 1900, res = 220, type = "cairo")
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(mfrow = c(2, 2), mar = c(4.4, 4.6, 3.5, 1.1), oma = c(0, 0, 2.2, 0))

  y_range <- range(
    trajectory$lower,
    trajectory$upper,
    observed$beta - 2 * observed$SE,
    observed$beta + 2 * observed$SE,
    finite = TRUE
  )
  plot(
    NA,
    xlim = range(fine_grid),
    ylim = y_range,
    xlab = "Time (days)",
    ylab = "eQTL effect",
    main = "A. Same posterior fit on nested grids"
  )
  rect(4, y_range[1], 11, y_range[2], col = "#E8F1E8", border = NA)
  abline(h = 0, col = "#777777", lty = 3)
  polygon(
    c(trajectory$time, rev(trajectory$time)),
    c(trajectory$lower, rev(trajectory$upper)),
    col = adjustcolor("#D55E00", alpha.f = 0.18),
    border = NA
  )
  arrows(
    observed$time,
    observed$beta - 2 * observed$SE,
    observed$time,
    observed$beta + 2 * observed$SE,
    angle = 90,
    code = 3,
    length = 0.035,
    col = "#4D4D4D"
  )
  points(observed$time, observed$beta, pch = 16, col = "#222222")
  lines(trajectory$time, trajectory$posterior_mean, col = "#D55E00", lwd = 2.3)
  points(
    coarse_trajectory$time,
    coarse_trajectory$posterior_mean,
    pch = 1,
    cex = 0.55,
    col = "#0072B2"
  )
  abline(v = c(4, 11), col = "#4D4D4D", lty = 2)
  legend(
    "topright",
    legend = c("0.05-grid posterior mean", "0.10-grid nodes", "95% credible interval", "Middle window"),
    col = c("#D55E00", "#0072B2", adjustcolor("#D55E00", alpha.f = 0.35), "#8AA58A"),
    lty = c(1, NA, 1, 1),
    pch = c(NA, 1, NA, NA),
    lwd = c(2.3, NA, 7, 7),
    bty = "n",
    cex = 0.72
  )

  zoom_rows <- trajectory$time >= zoom_limits[1] & trajectory$time <= zoom_limits[2]
  zoom_y <- range(trajectory$posterior_mean[zoom_rows], finite = TRUE)
  zoom_padding <- max(diff(zoom_y) * 0.25, 0.02)
  zoom_y <- zoom_y + c(-zoom_padding, zoom_padding)
  plot(
    NA,
    xlim = zoom_limits,
    ylim = zoom_y,
    xlab = "Time (days)",
    ylab = "Posterior mean",
    main = "B. Zoom at the influential added node"
  )
  rect(4, zoom_y[1], 11, zoom_y[2], col = "#E8F1E8", border = NA)
  lines(trajectory$time, trajectory$posterior_mean, col = "#D55E00", lwd = 2.5)
  lines(coarse_trajectory$time, coarse_trajectory$posterior_mean, col = "#0072B2", lwd = 1.7, lty = 2)
  points(coarse_trajectory$time, coarse_trajectory$posterior_mean, pch = 1, col = "#0072B2")
  added_zoom <- added_rows[fine_grid[added_rows] >= zoom_limits[1] & fine_grid[added_rows] <= zoom_limits[2]]
  points(
    trajectory$time[added_zoom],
    trajectory$posterior_mean[added_zoom],
    pch = 16,
    cex = 0.75,
    col = "#D55E00"
  )
  abline(v = c(4, 11), col = "#4D4D4D", lty = 2)
  legend(
    "topright",
    legend = c("0.05 grid", "linear connection of 0.10 nodes", "new 0.05 nodes"),
    col = c("#D55E00", "#0072B2", "#D55E00"),
    lty = c(1, 2, NA),
    pch = c(NA, NA, 16),
    bty = "n",
    cex = 0.72
  )

  density_coarse <- density(coarse_components$statistic)
  density_fine <- density(fine_components$statistic)
  density_ylim <- range(c(density_coarse$y, density_fine$y))
  plot(
    density_coarse,
    col = "#0072B2",
    lwd = 2,
    xlab = "Middle statistic: max inside - max outside",
    ylab = "Posterior density",
    ylim = density_ylim,
    main = "C. Functional changes despite the same fit"
  )
  lines(density_fine, col = "#D55E00", lwd = 2)
  abline(v = 0, col = "#333333", lty = 2)
  legend(
    "topright",
    legend = c(
      sprintf("0.10 grid: lfsr %.3f", lfsr_coarse["estimate"]),
      sprintf("0.05 grid: lfsr %.3f", lfsr_fine["estimate"])
    ),
    col = c("#0072B2", "#D55E00"),
    lwd = 2,
    bty = "n",
    cex = 0.78
  )

  plot_rows <- unique(round(seq(1, posterior_draws, length.out = min(10000L, posterior_draws))))
  changed_to_error <- !coarse_components$error[plot_rows] & fine_components$error[plot_rows]
  point_color <- ifelse(changed_to_error, "#D55E00", adjustcolor("#777777", alpha.f = 0.28))
  plot(
    coarse_components$statistic[plot_rows],
    fine_components$statistic[plot_rows],
    pch = 16,
    cex = 0.45,
    col = point_color,
    xlab = "Middle statistic on 0.10 grid",
    ylab = "Middle statistic on 0.05 grid",
    main = "D. Posterior draws that cross zero"
  )
  abline(a = 0, b = 1, col = "#0072B2", lty = 2)
  abline(h = 0, v = 0, col = "#333333", lty = 3)
  legend(
    "topleft",
    legend = c("Becomes non-Middle on 0.05 grid", "Other posterior draws"),
    col = c("#D55E00", "#777777"),
    pch = 16,
    bty = "n",
    cex = 0.72
  )
  mtext(pair_id, outer = TRUE, cex = 1.15, font = 2)

  invisible(c(
    lfsr_0p10 = lfsr_coarse["estimate"],
    lfsr_0p05 = lfsr_fine["estimate"],
    paired_difference_mcse = paired_mcse,
    dominant_added_outside_time = dominant_added_time
  ))
}

workflowr_root <- find_workflowr_root()
posterior_draws <- as.integer(get_arg("--posterior-draws", "30000"))
seed <- as.integer(get_arg("--seed", "20260819"))
output_id <- get_arg("--output-id", "evaluation_grid_middle_pair_visual_diagnostic")

if (posterior_draws < 3000L || is.na(seed) || !nzchar(output_id)) {
  stop("Invalid diagnostic arguments.")
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

suppressPackageStartupMessages(library(fashr))
if (!requireNamespace("matrixStats", quietly = TRUE)) {
  stop("The matrixStats package is required.")
}

fit_path <- file.path(workflowr_root, "output", "dynamic_eQTL_real", "fash_fit1_update.RData")
datasets_path <- file.path(workflowr_root, "output", "dynamic_eQTL_real", "datasets_corrected.RData")
reclassification_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "evaluation_grid_full_0p05_reclassification_pilot",
  "category_reclassification_by_pair.csv"
)
output_dir <- file.path(workflowr_root, "output", "revision_simulations", "internal", output_id)

if (!file.exists(fit_path) || !file.exists(datasets_path) || !file.exists(reclassification_path)) {
  stop("At least one required input is missing.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

reclassification <- read.csv(reclassification_path, stringsAsFactors = FALSE, check.names = FALSE)
middle <- reclassification[reclassification$category == "middle", , drop = FALSE]
middle <- middle[order(middle$lfsr_0p05, middle$pair_id), , drop = FALSE]
ripk2_id <- "ENSG00000104312_rs58780282"
target_ids <- unique(c(ripk2_id, middle$pair_id[1]))
targets <- middle[match(target_ids, middle$pair_id), , drop = FALSE]
if (anyNA(targets$index) || nrow(targets) != 2L) {
  stop("Could not identify both diagnostic pairs.")
}
targets$diagnostic_role <- c("manuscript_middle_example", "minimum_0p05_lfsr_in_full_run")

data_environment <- new.env(parent = emptyenv())
load(datasets_path, envir = data_environment)
if (!exists("datasets", envir = data_environment, inherits = FALSE)) {
  stop("The observed-data file does not contain datasets.")
}
observed_list <- lapply(seq_len(nrow(targets)), function(i) {
  data <- data_environment$datasets[[targets$index[i]]]
  if (!identical(names(data_environment$datasets)[targets$index[i]], targets$pair_id[i])) {
    stop("Observed-data indexing does not match the pair identifier.")
  }
  data.frame(time = data$time, beta = data$beta, SE = data$SE)
})
rm(data_environment)
gc()

load(fit_path)
if (!exists("fash_fit1_update") || !inherits(fash_fit1_update, "fash")) {
  stop("The fitted-model file does not contain fash_fit1_update.")
}

fine_grid <- seq(0, 15, by = 0.05)
coarse_rows <- which(abs(fine_grid / 0.10 - round(fine_grid / 0.10)) < 1e-8)
coarse_grid <- fine_grid[coarse_rows]

diagnostic_rows <- vector("list", nrow(targets))
draw_outputs <- vector("list", nrow(targets))
trajectory_outputs <- vector("list", nrow(targets))

for (i in seq_len(nrow(targets))) {
  pair_index <- targets$index[i]
  set.seed(seed + pair_index)
  samples <- predict(
    fash_fit1_update,
    index = pair_index,
    smooth_var = fine_grid,
    only.samples = TRUE,
    M = posterior_draws
  )
  if (!is.matrix(samples) || nrow(samples) != length(fine_grid) ||
      ncol(samples) != posterior_draws || any(!is.finite(samples))) {
    stop("Posterior sampling failed for ", targets$pair_id[i], ".")
  }

  fine_components <- middle_components(samples, fine_grid)
  coarse_components <- middle_components(samples[coarse_rows, , drop = FALSE], coarse_grid)
  fine_lfsr <- binomial_summary(fine_components$error)
  coarse_lfsr <- binomial_summary(coarse_components$error)
  paired_difference <- as.numeric(fine_components$error) - as.numeric(coarse_components$error)
  delta_inside <- fine_components$inside_max - coarse_components$inside_max
  delta_outside <- fine_components$outside_max - coarse_components$outside_max

  trajectory <- summarize_trajectory(samples, fine_grid)
  attr(trajectory, "posterior_samples") <- samples
  figure_path <- file.path(output_dir, paste0(gsub("[^A-Za-z0-9_]+", "_", targets$pair_id[i]), "_diagnostic.png"))
  plot_metrics <- plot_pair_diagnostic(
    pair_id = targets$pair_id[i],
    observed = observed_list[[i]],
    trajectory = trajectory,
    fine_grid = fine_grid,
    coarse_rows = coarse_rows,
    fine_components = fine_components,
    coarse_components = coarse_components,
    output_path = figure_path,
    posterior_draws = posterior_draws
  )
  attr(trajectory, "posterior_samples") <- NULL

  diagnostic_rows[[i]] <- data.frame(
    diagnostic_role = targets$diagnostic_role[i],
    pair_id = targets$pair_id[i],
    index = pair_index,
    posterior_draws = posterior_draws,
    saved_lfsr_0p10 = targets$saved_lfsr_0p10[i],
    saved_cfsr_0p10 = targets$saved_cfsr_0p10[i],
    full_run_lfsr_0p05_m3000 = targets$lfsr_0p05[i],
    full_run_cfsr_0p05_m3000 = targets$cfsr_0p05[i],
    high_draw_lfsr_0p10 = coarse_lfsr["estimate"],
    high_draw_lfsr_0p10_mcse = coarse_lfsr["mcse"],
    high_draw_lfsr_0p10_ci_lower = coarse_lfsr["ci_lower"],
    high_draw_lfsr_0p10_ci_upper = coarse_lfsr["ci_upper"],
    high_draw_lfsr_0p05 = fine_lfsr["estimate"],
    high_draw_lfsr_0p05_mcse = fine_lfsr["mcse"],
    high_draw_lfsr_0p05_ci_lower = fine_lfsr["ci_lower"],
    high_draw_lfsr_0p05_ci_upper = fine_lfsr["ci_upper"],
    paired_lfsr_difference_0p05_minus_0p10 = mean(paired_difference),
    paired_difference_mcse = sd(paired_difference) / sqrt(posterior_draws),
    discordant_draw_fraction = mean(fine_components$error != coarse_components$error),
    fraction_becoming_non_middle = mean(!coarse_components$error & fine_components$error),
    fraction_becoming_middle = mean(coarse_components$error & !fine_components$error),
    mean_added_inside_max = mean(delta_inside),
    mean_added_outside_max = mean(delta_outside),
    fraction_added_outside_exceeds_added_inside = mean(delta_outside > delta_inside),
    dominant_added_outside_time = unname(plot_metrics["dominant_added_outside_time"]),
    figure_path = normalizePath(figure_path),
    stringsAsFactors = FALSE
  )

  draw_outputs[[i]] <- data.frame(
    pair_id = targets$pair_id[i],
    statistic_0p10 = coarse_components$statistic,
    statistic_0p05 = fine_components$statistic,
    error_0p10 = coarse_components$error,
    error_0p05 = fine_components$error,
    delta_inside_max = delta_inside,
    delta_outside_max = delta_outside,
    stringsAsFactors = FALSE
  )
  trajectory$pair_id <- targets$pair_id[i]
  trajectory_outputs[[i]] <- trajectory
  rm(samples, fine_components, coarse_components)
  gc()
}

diagnostic_summary <- do.call(rbind, diagnostic_rows)
draw_diagnostics <- do.call(rbind, draw_outputs)
trajectory_summary <- do.call(rbind, trajectory_outputs)

write_csv(diagnostic_summary, file.path(output_dir, "middle_pair_grid_diagnostic_summary.csv"))
write_csv(trajectory_summary, file.path(output_dir, "middle_pair_posterior_trajectory_summary.csv"))
saveRDS(draw_diagnostics, file.path(output_dir, "middle_pair_functional_draw_diagnostics.rds"))
saveRDS(
  list(
    posterior_draws = posterior_draws,
    seed = seed,
    fine_grid_step = 0.05,
    coarse_grid_step = 0.10,
    middle_window = c(4, 11),
    fitted_model = normalizePath(fit_path),
    reclassification_input = normalizePath(reclassification_path)
  ),
  file.path(output_dir, "configuration.rds")
)

print(diagnostic_summary)
