#!/usr/bin/env Rscript

# Plot fixed-seed middle-functional LFSR agreement at 3,000 and 10,000
# posterior draws against the 30,000-draw reference on the 0.05-day grid.

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

summarize_agreement <- function(reference, comparison, posterior_draws) {
  absolute_difference <- abs(comparison - reference)
  data.frame(
    posterior_draws = posterior_draws,
    n_pairs = length(reference),
    mean_absolute_difference = mean(absolute_difference),
    median_absolute_difference = median(absolute_difference),
    q90_absolute_difference = unname(quantile(absolute_difference, 0.90)),
    maximum_absolute_difference = max(absolute_difference),
    spearman = suppressWarnings(cor(reference, comparison, method = "spearman")),
    pearson = suppressWarnings(cor(reference, comparison, method = "pearson")),
    fraction_within_0p005 = mean(absolute_difference <= 0.005),
    fraction_within_0p01 = mean(absolute_difference <= 0.01),
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
input_id <- get_arg("--input-id", "evaluation_grid_middle_grid_mc_convergence")
output_name <- get_arg("--output-name", "middle_lfsr_mc_sample_scatter")
output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  input_id
)
input_path <- file.path(output_dir, "pair_lfsr_by_grid_and_mc_size.csv")

if (!file.exists(input_path) || !nzchar(output_name)) {
  stop("The convergence input is missing or the output name is invalid.")
}

results <- read.csv(input_path, stringsAsFactors = FALSE, check.names = FALSE)
data <- results[
  results$method == "direct_fixed_seed" &
    results$grid_step == "0.05" &
    results$posterior_draws %in% c(3000L, 10000L, 30000L),
  ,
  drop = FALSE
]

pair_metadata <- unique(data[, c(
  "pair_id",
  "gene_symbol",
  "variant_id",
  "was_selected_0p10"
)])
if (nrow(pair_metadata) != 53L || anyDuplicated(pair_metadata$pair_id)) {
  stop("The expected 53-pair candidate universe was not recovered.")
}

lfsr_by_draw <- lapply(c(3000L, 10000L, 30000L), function(posterior_draws) {
  current <- data[data$posterior_draws == posterior_draws, c("pair_id", "lfsr")]
  setNames(current$lfsr, current$pair_id)
})
names(lfsr_by_draw) <- c("3000", "10000", "30000")

pair_order <- pair_metadata$pair_id
comparison <- data.frame(
  pair_id = pair_order,
  gene_symbol = pair_metadata$gene_symbol,
  variant_id = pair_metadata$variant_id,
  was_selected_0p10 = pair_metadata$was_selected_0p10,
  lfsr_3000 = unname(lfsr_by_draw[["3000"]][pair_order]),
  lfsr_10000 = unname(lfsr_by_draw[["10000"]][pair_order]),
  lfsr_30000 = unname(lfsr_by_draw[["30000"]][pair_order]),
  stringsAsFactors = FALSE
)

if (anyNA(comparison[, c("lfsr_3000", "lfsr_10000", "lfsr_30000")])) {
  stop("At least one pair is missing from a sample-size comparison.")
}

agreement <- rbind(
  summarize_agreement(comparison$lfsr_30000, comparison$lfsr_3000, 3000L),
  summarize_agreement(comparison$lfsr_30000, comparison$lfsr_10000, 10000L)
)

plot_path <- file.path(output_dir, paste0(output_name, ".png"))
png(plot_path, width = 2300, height = 1150, res = 220, type = "cairo")
old_par <- par(no.readonly = TRUE)
on.exit({
  par(old_par)
  dev.off()
}, add = TRUE)
par(mfrow = c(1, 2), mar = c(4.7, 4.9, 3.7, 1.2), oma = c(0, 0, 2.2, 0))

axis_limit <- c(0, max(comparison[, c("lfsr_3000", "lfsr_10000", "lfsr_30000")]) * 1.06)
point_colors <- ifelse(comparison$was_selected_0p10, "#D55E00", "#8C8C8C")
point_order <- order(comparison$was_selected_0p10)

for (panel_index in seq_along(c(3000L, 10000L))) {
  posterior_draws <- c(3000L, 10000L)[panel_index]
  y <- comparison[[if (posterior_draws == 3000L) "lfsr_3000" else "lfsr_10000"]]
  current_agreement <- agreement[agreement$posterior_draws == posterior_draws, , drop = FALSE]

  plot(
    NA,
    xlim = axis_limit,
    ylim = axis_limit,
    xlab = "LFSR with 30,000 posterior draws",
    ylab = sprintf("LFSR with %s posterior draws", format(posterior_draws, big.mark = ",")),
    main = if (posterior_draws == 3000L) {
      "A. 3,000 versus 30,000 draws"
    } else {
      "B. 10,000 versus 30,000 draws"
    }
  )
  polygon(
    x = c(axis_limit[1], axis_limit[2], axis_limit[2], axis_limit[1]),
    y = c(
      axis_limit[1] - 0.005,
      axis_limit[2] - 0.005,
      axis_limit[2] + 0.005,
      axis_limit[1] + 0.005
    ),
    col = adjustcolor("#56B4E9", alpha.f = 0.12),
    border = NA
  )
  abline(a = 0, b = 1, col = "#333333", lty = 2, lwd = 1.5)
  abline(h = 0.05, v = 0.05, col = "#999999", lty = 3)
  points(
    comparison$lfsr_30000[point_order],
    y[point_order],
    pch = 16,
    cex = 0.92,
    col = adjustcolor(point_colors[point_order], alpha.f = 0.84)
  )

  annotation <- sprintf(
    "Mean |difference| = %.4f\n90th percentile = %.4f\nSpearman rho = %.3f\nWithin 0.005 = %.0f%%",
    current_agreement$mean_absolute_difference,
    current_agreement$q90_absolute_difference,
    current_agreement$spearman,
    100 * current_agreement$fraction_within_0p005
  )
  legend(
    "topleft",
    legend = annotation,
    bty = "n",
    cex = 0.78,
    text.col = "#333333"
  )
}

legend(
  "bottomright",
  legend = c("Original 0.10-grid Middle discovery", "Other near-threshold candidate", "Absolute difference <= 0.005"),
  col = c("#D55E00", "#8C8C8C", adjustcolor("#56B4E9", alpha.f = 0.28)),
  pch = c(16, 16, 15),
  pt.cex = c(0.9, 0.9, 1.5),
  bty = "n",
  cex = 0.72
)

mtext(
  "Middle-functional LFSR agreement on the 0.05-day grid",
  outer = TRUE,
  font = 2,
  cex = 1.15
)

write.csv(comparison, file.path(output_dir, paste0(output_name, "_data.csv")), row.names = FALSE)
write.csv(agreement, file.path(output_dir, paste0(output_name, "_summary.csv")), row.names = FALSE)
cat("Wrote:", normalizePath(plot_path), "\n")
print(agreement)
