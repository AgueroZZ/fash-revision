#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

repository_root <- normalizePath(getwd(), mustWork = TRUE)
input_path <- file.path(
  repository_root,
  "output/revision_simulations/mc",
  "r1_random_bspline_main_effect_profile_sigma_pilot5",
  "full_fits/seed_12345.rds"
)
output_dir <- file.path(
  repository_root,
  "output/revision_simulations/internal",
  "r1_permuted_genotype_alignment_vs_n"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(20260816)
n_grid <- c(8L, 10L, 15L, 19L, 25L, 40L, 75L, 150L, 300L)
n_permutations <- 20000L

simulate_fixed_vector_alignment <- function(n_donors, n_permutations) {
  genotype <- rnorm(n_donors)
  genotype <- genotype - mean(genotype)
  genotype_ss <- sum(genotype^2)
  correlations <- replicate(
    n_permutations,
    sum(genotype * genotype[sample.int(n_donors)]) / genotype_ss
  )
  data.frame(
    n_donors = n_donors,
    n_permutations = n_permutations,
    empirical_mean = mean(correlations),
    empirical_sd = sd(correlations),
    theoretical_sd = 1 / sqrt(n_donors - 1),
    mean_absolute_correlation = mean(abs(correlations)),
    probability_abs_gt_0_2 = mean(abs(correlations) > 0.2),
    probability_abs_gt_0_3 = mean(abs(correlations) > 0.3),
    q025 = unname(quantile(correlations, 0.025)),
    q50 = unname(quantile(correlations, 0.5)),
    q975 = unname(quantile(correlations, 0.975)),
    stringsAsFactors = FALSE
  )
}

scaling_summary <- do.call(
  rbind,
  lapply(
    n_grid,
    simulate_fixed_vector_alignment,
    n_permutations = n_permutations
  )
)

if (!file.exists(input_path)) {
  stop("The saved R1 full-fit object is missing: ", input_path)
}
r1_fit <- readRDS(input_path)
dynamic_index <- which(r1_fit$unit_info$effect_class == "dynamic_bspline")
genotype <- r1_fit$genotype[, dynamic_index, drop = FALSE]
covariates <- r1_fit$covariates
design <- cbind(intercept = 1, covariates)
design_qr <- qr(design)
if (design_qr$rank != ncol(design)) {
  stop("The R1 nuisance-covariate design is rank deficient.")
}
residualizer <- diag(nrow(design)) - tcrossprod(qr.Q(design_qr))
original_residual <- residualizer %*% genotype
original_ss <- colSums(original_residual^2)

set.seed(20260817)
n_r1_permutations <- 5000L
r1_partial_correlations <- matrix(
  NA_real_,
  nrow = n_r1_permutations,
  ncol = ncol(genotype)
)
r1_leakage_coefficients <- matrix(
  NA_real_,
  nrow = n_r1_permutations,
  ncol = ncol(genotype)
)
for (permutation_index in seq_len(n_r1_permutations)) {
  donor_permutation <- sample.int(nrow(genotype))
  permuted_residual <- residualizer %*%
    genotype[donor_permutation, , drop = FALSE]
  permuted_ss <- colSums(permuted_residual^2)
  cross_product <- colSums(original_residual * permuted_residual)
  r1_partial_correlations[permutation_index, ] <-
    cross_product / sqrt(original_ss * permuted_ss)
  r1_leakage_coefficients[permutation_index, ] <-
    cross_product / permuted_ss
}

r1_correlation_vector <- as.vector(r1_partial_correlations)
r1_leakage_vector <- as.vector(r1_leakage_coefficients)
r1_summary <- data.frame(
  seed = 12345L,
  n_donors = nrow(genotype),
  n_covariates = ncol(covariates),
  residual_subspace_dimension = nrow(genotype) - ncol(design),
  n_dynamic_units = ncol(genotype),
  n_shared_donor_permutations = n_r1_permutations,
  mean_partial_correlation = mean(r1_correlation_vector),
  sd_partial_correlation = sd(r1_correlation_vector),
  mean_absolute_partial_correlation = mean(abs(r1_correlation_vector)),
  probability_abs_gt_0_2 = mean(abs(r1_correlation_vector) > 0.2),
  probability_abs_gt_0_3 = mean(abs(r1_correlation_vector) > 0.3),
  q025_partial_correlation = unname(quantile(r1_correlation_vector, 0.025)),
  q50_partial_correlation = unname(quantile(r1_correlation_vector, 0.5)),
  q975_partial_correlation = unname(quantile(r1_correlation_vector, 0.975)),
  rms_leakage_coefficient = sqrt(mean(r1_leakage_vector^2)),
  q95_absolute_leakage_coefficient = unname(
    quantile(abs(r1_leakage_vector), 0.95)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  scaling_summary,
  file.path(output_dir, "fixed_centered_genotype_scaling_summary.csv"),
  row.names = FALSE
)
write.csv(
  r1_summary,
  file.path(output_dir, "r1_covariate_adjusted_alignment_summary.csv"),
  row.names = FALSE
)

scaling_plot <- ggplot(scaling_summary, aes(x = n_donors)) +
  geom_line(
    aes(y = theoretical_sd, color = "Exact SD: 1/sqrt(N - 1)"),
    linewidth = 0.8
  ) +
  geom_point(
    aes(y = empirical_sd, color = "Permutation simulation"),
    size = 2.2
  ) +
  geom_point(
    data = r1_summary,
    aes(
      x = n_donors,
      y = sd_partial_correlation,
      color = "R1 partial correlation (5 PCs)"
    ),
    shape = 17,
    size = 3.1,
    inherit.aes = FALSE
  ) +
  scale_x_log10(breaks = n_grid) +
  scale_y_log10() +
  scale_color_manual(
    values = c(
      "Exact SD: 1/sqrt(N - 1)" = "#1B4965",
      "Permutation simulation" = "#2A9D8F",
      "R1 partial correlation (5 PCs)" = "#D1495B"
    )
  ) +
  labs(
    x = "Number of donors (N; log scale)",
    y = "SD of cor(G, PG) (log scale)",
    color = NULL,
    title = "Random permutation alignment decreases only at the N^(-1/2) rate",
    subtitle = "The exact curve assumes one fixed centered genotype vector and no additional covariates"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
ggsave(
  file.path(output_dir, "permuted_genotype_alignment_sd_vs_n.png"),
  scaling_plot,
  width = 7.5,
  height = 5.2,
  dpi = 220
)

distribution_data <- rbind(
  data.frame(
    correlation = replicate(
      n_permutations,
      {
        genotype_n19 <- scale(rnorm(19L), center = TRUE, scale = FALSE)[, 1L]
        sum(genotype_n19 * genotype_n19[sample.int(19L)]) /
          sum(genotype_n19^2)
      }
    ),
    source = "Centered genotype; no PCs"
  ),
  data.frame(
    correlation = r1_correlation_vector,
    source = "R1 dynamic genotypes; adjusted for 5 PCs"
  )
)
distribution_plot <- ggplot(
  distribution_data,
  aes(x = correlation, fill = source)
) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 70,
    alpha = 0.45,
    position = "identity"
  ) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5) +
  scale_fill_manual(values = c("#2A9D8F", "#D1495B")) +
  labs(
    x = "Signed genotype-permutation correlation",
    y = "Density",
    fill = NULL,
    title = "At N = 19, a fixed permutation is not close to orthogonal",
    subtitle = "PC adjustment changes the distribution but does not force zero alignment"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )
ggsave(
  file.path(output_dir, "permuted_genotype_alignment_distribution_n19.png"),
  distribution_plot,
  width = 7.5,
  height = 5.2,
  dpi = 220
)

print(scaling_summary)
print(r1_summary)
