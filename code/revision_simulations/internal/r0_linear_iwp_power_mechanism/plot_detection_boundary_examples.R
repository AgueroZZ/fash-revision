#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("output/revision_simulations")) {
    return(".")
  }
  if (file.exists("coderepo-local/output/revision_simulations")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr project root.")
}

workflowr_root <- find_workflowr_root()
input_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r0_linear_iwp_power_mechanism_clean_gaussian_grid_bins",
  "bin_summary.csv"
)
output_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r0_linear_iwp_power_mechanism_clean_gaussian_grid_bins",
  "detection_boundary_examples.png"
)

power_summary <- utils::read.csv(input_path, check.names = FALSE)
time_grid <- 0:31
endpoint_examples <- c(1.5, 2.25, 3.5)
example_bins <- c("[1,2)", "[2,2.5)", "[3,4)")
example_labels <- c("Below boundary", "At boundary", "Above boundary")

# Reuse one fixed noise realization so the visual comparison isolates effect size.
set.seed(20260820)
shared_noise <- stats::rnorm(length(time_grid))
observed_examples <- lapply(endpoint_examples, function(endpoint) {
  true_mean <- endpoint * time_grid / max(time_grid)
  list(true_mean = true_mean, observed = true_mean + shared_noise)
})

linear_color <- "#0072B2"
iwp_color <- "#D55E00"
truth_color <- "#222222"
point_color <- "#8A8A8A"

grDevices::png(output_path, width = 2100, height = 1450, res = 180)
layout(matrix(c(1, 2, 3, 4, 4, 4), nrow = 2L, byrow = TRUE), heights = c(1, 1.18))
graphics::par(
  family = "sans",
  mar = c(4.4, 4.6, 3.8, 1.2),
  mgp = c(2.6, 0.8, 0),
  tcl = -0.25
)

for (index in seq_along(endpoint_examples)) {
  endpoint <- endpoint_examples[[index]]
  example <- observed_examples[[index]]
  bin_row <- power_summary[power_summary$effect_bin == example_bins[[index]], ]
  graphics::plot(
    time_grid,
    example$observed,
    type = "n",
    ylim = c(-2.5, 5.2),
    xlab = "Time",
    ylab = if (index == 1L) "Observed effect estimate" else "",
    main = sprintf("%s: |endpoint change| = %.2f", example_labels[[index]], endpoint),
    cex.main = 0.95
  )
  graphics::abline(h = 0, col = "#DDDDDD", lty = 3)
  graphics::segments(
    x0 = time_grid,
    y0 = example$true_mean - 1,
    x1 = time_grid,
    y1 = example$true_mean + 1,
    col = grDevices::adjustcolor(point_color, alpha.f = 0.18)
  )
  graphics::points(
    time_grid,
    example$observed,
    pch = 16,
    cex = 0.55,
    col = grDevices::adjustcolor(point_color, alpha.f = 0.75)
  )
  graphics::lines(time_grid, example$true_mean, lwd = 2.6, col = truth_color)
  graphics::mtext(
    sprintf(
      "Bin power: linear %.0f%%; IWP %.0f%%",
      100 * bin_row$mean_power_linear,
      100 * bin_row$mean_power_iwp
    ),
    side = 3,
    line = 0.25,
    cex = 0.78,
    col = "#555555"
  )
}

graphics::par(mar = c(5.8, 4.8, 3.6, 1.2))
x_positions <- seq_len(nrow(power_summary))
graphics::plot(
  x_positions,
  power_summary$mean_power_linear,
  type = "n",
  ylim = c(0, 1.04),
  xaxt = "n",
  xlab = "Absolute true endpoint change, |f(31) - f(0)|",
  ylab = "Discovery probability (power)",
  main = "Empirical transition across 20 held-out seeds",
  cex.main = 1.05
)
graphics::rect(
  2.5,
  0,
  4.5,
  1.04,
  col = grDevices::adjustcolor("#F0E442", alpha.f = 0.18),
  border = NA
)
graphics::abline(h = c(0.25, 0.5, 0.75), col = "#E5E5E5", lty = 3)

draw_ci <- function(x, estimate, lower, upper, color) {
  graphics::segments(x, lower, x, upper, col = color, lwd = 1.6)
  graphics::segments(x - 0.04, lower, x + 0.04, lower, col = color, lwd = 1.3)
  graphics::segments(x - 0.04, upper, x + 0.04, upper, col = color, lwd = 1.3)
  graphics::lines(x, estimate, col = color, lwd = 2.6)
  graphics::points(x, estimate, col = color, bg = "white", pch = 21, cex = 1.05, lwd = 1.7)
}

linear_x <- x_positions - 0.07
iwp_x <- x_positions + 0.07
draw_ci(
  linear_x,
  power_summary$mean_power_linear,
  power_summary$mean_power_linear_lower,
  power_summary$mean_power_linear_upper,
  linear_color
)
draw_ci(
  iwp_x,
  power_summary$mean_power_iwp,
  power_summary$mean_power_iwp_lower,
  power_summary$mean_power_iwp_upper,
  iwp_color
)
graphics::axis(1, at = x_positions, labels = power_summary$effect_bin)
graphics::legend(
  "bottomright",
  legend = c("FASH-linear", "FASH-IWP", "Detection-boundary region"),
  col = c(linear_color, iwp_color, NA),
  pt.bg = c("white", "white", grDevices::adjustcolor("#F0E442", alpha.f = 0.35)),
  pch = c(21, 21, 22),
  lwd = c(2.6, 2.6, NA),
  pt.cex = c(1, 1, 1.4),
  bty = "n"
)
graphics::mtext(
  "Clean R0 setting: 32 time points, independent N(0,1) errors, exact SEs, pi0 = 0.95, FDR target = 0.05",
  side = 1,
  line = 4.25,
  cex = 0.82,
  col = "#555555"
)

grDevices::dev.off()
message(normalizePath(output_path))
