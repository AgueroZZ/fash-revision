#!/usr/bin/env Rscript

# Derive a posterior-outcome-independent temporal-category reference mixture
# from the canonical IWP1 prior used by the R3 FASH fits. The fitted PSD weights
# are deliberately not used. Final formal probabilities are rounded to two
# decimals before any new simulation is run.

options(stringsAsFactors = FALSE)

project_root <- normalizePath(getwd(), mustWork = TRUE)
output_dir <- file.path(
  project_root,
  "output/revision_simulations/diagnostics/r3_middle_calibration"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

evaluation_grid <- seq(0, 15, by = 0.1)
category_by_grid <- ifelse(
  evaluation_grid <= 3,
  "early",
  ifelse(evaluation_grid < 12, "middle", "late")
)
category_levels <- c("early", "middle", "late")

total_draws <- 1000000L
chunk_size <- 10000L
reference_seed <- 8675309L
canonical_settings <- list(
  n_basis = 20L,
  psd = 1,
  sd_poly = 1 / sqrt(1e-6),
  p = 1L,
  pred_step = 1,
  x_range = c(0, 15)
)

counts <- setNames(integer(length(category_levels)), category_levels)
set.seed(reference_seed)
remaining <- total_draws
while (remaining > 0L) {
  current_chunk <- min(chunk_size, remaining)
  draws <- fashr:::simulate_IWP(
    n_samps = current_chunk,
    n_basis = canonical_settings$n_basis,
    psd = canonical_settings$psd,
    sd_poly = canonical_settings$sd_poly,
    p = canonical_settings$p,
    pred_step = canonical_settings$pred_step,
    x_range = canonical_settings$x_range,
    x_new = evaluation_grid
  )$samples
  maximum_index <- max.col(t(abs(draws)), ties.method = "first")
  labels <- factor(
    category_by_grid[maximum_index],
    levels = category_levels
  )
  counts <- counts + as.integer(table(labels))
  remaining <- remaining - current_chunk
}

probability <- counts / total_draws
mc_se <- sqrt(probability * (1 - probability) / total_draws)
rounded_probability <- c(early = 0.29, middle = 0.42, late = 0.29)

reference <- data.frame(
  category = category_levels,
  count = as.integer(counts),
  total_draws = total_draws,
  canonical_probability = as.numeric(probability),
  mc_se = as.numeric(mc_se),
  ci_lower = pmax(0, probability - 1.96 * mc_se),
  ci_upper = pmin(1, probability + 1.96 * mc_se),
  frozen_formal_probability = as.numeric(rounded_probability),
  reference_seed = reference_seed,
  num_basis = canonical_settings$n_basis,
  psd = canonical_settings$psd,
  sd_poly = canonical_settings$sd_poly,
  order = canonical_settings$p,
  pred_step = canonical_settings$pred_step,
  stringsAsFactors = FALSE
)

write.csv(
  reference,
  file.path(output_dir, "canonical_iwp1_reference_mixture.csv"),
  row.names = FALSE
)

cat("Canonical IWP1 temporal-category reference:\n")
print(reference, row.names = FALSE)
