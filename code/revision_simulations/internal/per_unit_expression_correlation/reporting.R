#!/usr/bin/env Rscript

# Figure: is the common C more correlated than the discoveries' own units are?

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}
workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "internal",
                 "per_unit_expression_correlation",
                 "per_unit_expression_correlation_helpers.R"))

output_dir <- file.path(workflowr_root, "output", "revision_simulations",
                        "internal", "per_unit_expression_correlation")
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
result <- readRDS(file.path(output_dir, "per_unit_expression_correlation.rds"))

col_common <- "#1b1b1b"
col_null <- "#0072B2"
col_disc <- "#D55E00"
col_surv <- "#009E73"

per_unit <- result$per_unit
beta <- per_unit[per_unit$estimator == "beta_scale", ]
panel_of <- function(keys, frame) frame$lag1[frame$pair_key %in% keys]
disc_lag1 <- panel_of(result$panels$discovery1169, beta)
null_lag1 <- panel_of(result$panels$null874, beta)
surv_lag1 <- panel_of(result$panels$survivor43, beta)
common_lag <- lag_profile(result$common_c)

lag_of <- function(panel, estimator) {
  rows <- result$lag_table
  rows$correlation[rows$panel == panel & rows$estimator == estimator]
}

png(file.path(figure_dir, "common_c_vs_per_unit_correlation.png"),
    width = 2200, height = 1750, res = 190)
layout(matrix(1:4, nrow = 2, byrow = TRUE))
par(mar = c(4.4, 4.6, 3.2, 1.2), mgp = c(2.7, 0.8, 0), cex.lab = 1.05)

# --- A: lag profiles, beta scale -------------------------------------------
plot(1:15, common_lag, type = "n", ylim = c(-0.05, 0.52),
     xlab = "Lag (days apart)", ylab = "Mean null correlation",
     main = "A. Beta-scale correlation by lag")
abline(h = 0, col = "grey80", lty = 3)
lines(1:15, common_lag, col = col_common, lwd = 4)
lines(1:15, lag_of("null874", "beta_scale"), col = col_null, lwd = 2.4)
lines(1:15, lag_of("discovery1169", "beta_scale"), col = col_disc, lwd = 2.4)
lines(1:15, lag_of("survivor43", "beta_scale"), col = col_surv, lwd = 2.2,
      lty = 2)
points(1:15, common_lag, pch = 16, col = col_common, cex = 0.8)
legend("topright", bty = "n", cex = 0.86, lwd = c(4, 2.4, 2.4, 2.2),
       lty = c(1, 1, 1, 2),
       col = c(col_common, col_null, col_disc, col_surv),
       legend = c("common C imposed on every unit",
                  "null panel, own units (n = 874)",
                  "discovery panel, own units (n = 1,169)",
                  "the 43 survivors of the common-C fit"))

# --- B: per-unit lag-1 distribution ----------------------------------------
density_null <- density(null_lag1, adjust = 1.1)
density_disc <- density(disc_lag1, adjust = 1.1)
plot(0, 0, type = "n", xlim = c(0.05, 0.8),
     ylim = c(0, max(density_null$y, density_disc$y) * 1.18),
     xlab = "Own beta-scale lag-1 correlation", ylab = "Density",
     main = "B. Per-unit lag-1, against the single imposed value")
polygon(c(density_null$x, rev(density_null$x)),
        c(density_null$y, rep(0, length(density_null$y))),
        col = adjustcolor <- rgb(0, 114, 178, 45, maxColorValue = 255),
        border = col_null, lwd = 2)
polygon(c(density_disc$x, rev(density_disc$x)),
        c(density_disc$y, rep(0, length(density_disc$y))),
        col = rgb(213, 94, 0, 45, maxColorValue = 255),
        border = col_disc, lwd = 2)
abline(v = common_lag[1], col = col_common, lwd = 3)
rug(surv_lag1, col = col_surv, lwd = 2, ticksize = 0.055)
text(common_lag[1], max(density_null$y, density_disc$y) * 1.13,
     paste0("common C = ", format(round(common_lag[1], 3), nsmall = 3)),
     pos = 2, cex = 0.85, col = col_common)
legend("topleft", bty = "n", cex = 0.84,
       fill = c(rgb(0, 114, 178, 45, maxColorValue = 255),
                rgb(213, 94, 0, 45, maxColorValue = 255), NA),
       border = c(col_null, col_disc, NA),
       col = c(NA, NA, col_surv), lwd = c(NA, NA, 2),
       legend = c("null panel", "discovery panel",
                  "43 survivors (rug)"))

# --- C: the genotype confound ---------------------------------------------
scales <- c("expr_pc_only", "expr_pc_genotype", "beta_scale")
labels <- c("expression,\nPCs out,\ngenotype IN",
            "expression,\nPCs + genotype\nout",
            "propagated to\nbeta scale")
means <- vapply(scales, function(nm) {
  frame <- per_unit[per_unit$estimator == nm, ]
  c(discovery = mean(panel_of(result$panels$discovery1169, frame)),
    null = mean(panel_of(result$panels$null874, frame)))
}, numeric(2L))
plot(0, 0, type = "n", xlim = c(0.7, 3.3), ylim = c(0.34, 0.60),
     xaxt = "n", xlab = "", ylab = "Mean lag-1 correlation",
     main = "C. Leaving genotype in hides the difference")
axis(1, at = 1:3, labels = FALSE)
mtext(labels, side = 1, at = 1:3, line = 2.3, cex = 0.72)
abline(h = common_lag[1], col = col_common, lwd = 2.4, lty = 2)
lines(1:3, means["null", ], col = col_null, lwd = 2.6)
lines(1:3, means["discovery", ], col = col_disc, lwd = 2.6)
points(1:3, means["null", ], pch = 16, col = col_null, cex = 1.5)
points(1:3, means["discovery", ], pch = 16, col = col_disc, cex = 1.5)
for (index in 1:3) {
  text(index, means["null", index], format(round(means["null", index], 3),
       nsmall = 3), pos = 3, cex = 0.8, col = col_null)
  text(index, means["discovery", index], format(round(means["discovery", index], 3),
       nsmall = 3), pos = 1, cex = 0.8, col = col_disc)
}
text(1.75, common_lag[1], "common C imposed on every unit", pos = 3,
     cex = 0.8, col = col_common)
legend("bottomleft", bty = "n", cex = 0.84, lwd = 2.6,
       col = c(col_null, col_disc),
       legend = c("null panel (n = 874)", "discovery panel (n = 1,169)"))

# --- D: over-correction map ------------------------------------------------
discovery_matrix <- result$panel_matrices[["discovery1169__beta_scale"]]
excess <- result$common_c - discovery_matrix
diag(excess) <- NA
limit <- max(abs(excess), na.rm = TRUE)
breaks <- seq(-limit, limit, length.out = 51)
palette <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(50)
image(0:15, 0:15, excess, col = palette, breaks = breaks,
      xlab = "Day", ylab = "Day",
      main = "D. common C minus discovery panel's own correlation")
box()
mtext(paste0("red = common C imposes more correlation than the unit has; ",
             "mean excess ",
             format(round(mean(excess, na.rm = TRUE), 3), nsmall = 3)),
      side = 3, line = 0.1, cex = 0.7)
dev.off()

message("Figure written to ",
        file.path(figure_dir, "common_c_vs_per_unit_correlation.png"))
