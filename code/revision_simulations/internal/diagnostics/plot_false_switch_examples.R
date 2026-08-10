#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit)) {
    return(default)
  }
  if (hit == length(args)) {
    stop("Missing value after ", flag, ".")
  }
  args[[hit + 1L]]
}

workflowr_root <- normalizePath(
  get_arg("--workflowr-root", "."),
  mustWork = TRUE
)
seed <- as.integer(get_arg("--seed", "12345"))
num_samples <- as.integer(get_arg("--num-samples", "3000"))
num_per_group <- as.integer(get_arg("--num-per-group", "3"))
output_path <- get_arg(
  "--output-path",
  file.path(
    "output",
    "revision_simulations",
    "diagnostics",
    "switch_miscalibration",
    "seed12345_bf_representative_false_switch_examples.png"
  )
)

if (!is.finite(num_samples) || num_samples < 500 ||
    !is.finite(num_per_group) || num_per_group < 1) {
  stop("Invalid plotting arguments.")
}

setwd(workflowr_root)
source("code/revision_simulations/shared/simulation_functions.R")
suppressPackageStartupMessages(library(fashr))

raw_candidates <- list.files(
  file.path("output", "revision_simulations", "raw"),
  pattern = paste0(
    "^genotype_sparse_timed_cosine_one_two_peak.*seed",
    seed,
    "[.]rds$"
  ),
  full.names = TRUE
)
if (length(raw_candidates) != 1L) {
  stop("Expected exactly one sparse timed cosine object for seed ", seed, ".")
}
out <- readRDS(raw_candidates)

diagnostic_path <- file.path(
  "output",
  "revision_simulations",
  "diagnostics",
  "switch_miscalibration",
  paste0("seed", seed, "_bf_alpha0p05_per_unit.csv")
)
if (!file.exists(diagnostic_path)) {
  stop("Run investigate_switch_miscalibration.R before plotting examples.")
}
per_unit <- read.csv(diagnostic_path, stringsAsFactors = FALSE)
false_calls <- per_unit[
  per_unit$false_functional_call &
    per_unit$truth_group %in% c("one-peak", "same-sign two-peak"),
  ,
  drop = FALSE
]

select_group <- function(group) {
  rows <- false_calls[false_calls$truth_group == group, , drop = FALSE]
  rows <- rows[order(rows$functional_lfsr, rows$index), , drop = FALSE]
  head(rows, num_per_group)
}
selected <- rbind(
  select_group("one-peak"),
  select_group("same-sign two-peak")
)
if (nrow(selected) != 2 * num_per_group) {
  stop("Not enough false switch calls in each requested truth group.")
}

fit <- out$fash_fits$fash_iwp1_bf
evaluation_grid <- out$evaluation_grid
observed_grid <- out$settings$time_grid

posterior <- lapply(selected$index, function(index) {
  set.seed(970000L + seed + index)
  samples <- predict(
    fit,
    index = index,
    smooth_var = evaluation_grid,
    only.samples = TRUE,
    M = num_samples
  )
  list(
    mean = rowMeans(samples),
    lower = apply(samples, 1, stats::quantile, probs = 0.025),
    upper = apply(samples, 1, stats::quantile, probs = 0.975)
  )
})

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
png(
  output_path,
  width = 2500,
  height = 1550,
  res = 200
)
old_par <- par(no.readonly = TRUE)
par(
  mfrow = c(2, num_per_group),
  mar = c(4.2, 4.3, 4.2, 1),
  oma = c(0, 0, 2.4, 0)
)

for (i in seq_len(nrow(selected))) {
  index <- selected$index[i]
  posterior_i <- posterior[[i]]
  estimate <- out$eqtl_summary$beta_hat[index, ]
  standard_error <- out$eqtl_summary$se[index, ]
  truth <- out$true_beta_evaluation[index, ]
  ylim <- range(
    truth,
    estimate - 1.96 * standard_error,
    estimate + 1.96 * standard_error,
    posterior_i$lower,
    posterior_i$upper
  )
  plot(
    evaluation_grid,
    posterior_i$mean,
    type = "n",
    ylim = ylim,
    xlab = "Time",
    ylab = "Genetic effect",
    main = paste0(
      selected$truth_group[i],
      "\nvariant ",
      index,
      "; lfsr = ",
      formatC(selected$functional_lfsr[i], format = "f", digits = 3),
      "; P(switch | data) = ",
      formatC(
        selected$posterior_switch_probability[i],
        format = "f",
        digits = 3
      )
    )
  )
  polygon(
    c(evaluation_grid, rev(evaluation_grid)),
    c(posterior_i$lower, rev(posterior_i$upper)),
    border = NA,
    col = grDevices::adjustcolor("#4C78A8", alpha.f = 0.18)
  )
  abline(h = c(-0.25, 0, 0.25), lty = c(3, 2, 3), col = "gray50")
  lines(evaluation_grid, truth, col = "#E45756", lwd = 2.8)
  lines(evaluation_grid, posterior_i$mean, col = "#4C78A8", lwd = 2.4)
  arrows(
    observed_grid,
    estimate - 1.96 * standard_error,
    observed_grid,
    estimate + 1.96 * standard_error,
    angle = 90,
    code = 3,
    length = 0.025,
    col = "gray45"
  )
  points(observed_grid, estimate, pch = 19, cex = 0.65)
  if (i == 1) {
    legend(
      "topright",
      legend = c(
        "True effect",
        "Posterior mean",
        "Estimate +/- 1.96 SE",
        "Switch thresholds"
      ),
      col = c("#E45756", "#4C78A8", "gray35", "gray50"),
      lty = c(1, 1, NA, 3),
      pch = c(NA, NA, 19, NA),
      lwd = c(2.8, 2.4, NA, 1),
      pt.cex = 0.7,
      bty = "n",
      cex = 0.78
    )
  }
}
mtext(
  "Representative false switch discoveries under the manuscript c = 0.25 functional",
  outer = TRUE,
  side = 3,
  line = 0.5,
  cex = 1.2,
  font = 2
)
par(old_par)
dev.off()

write.csv(
  selected,
  sub("[.]png$", ".csv", output_path),
  row.names = FALSE
)
cat(normalizePath(output_path), "\n")
