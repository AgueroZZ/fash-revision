#!/usr/bin/env Rscript

# Plot matched original and permuted FASH likelihood profiles.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}

if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("scales", quietly = TRUE)) {
  stop("The ggplot2 and scales packages are required.")
}

workflowr_root <- find_workflowr_root()
output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "selected_signal_genotype_permutation_seed20260811"
)
fit_path <- file.path(output_directory, "merged_fash_fit.rds")
figure_directory <- file.path(output_directory, "figures")
if (!file.exists(fit_path)) {
  stop("The fixed-seed merged FASH fit is missing.")
}
dir.create(figure_directory, showWarnings = FALSE)

sample_seed <- 20260811L
n_sample <- 100L
fit_bundle <- readRDS(fit_path)
raw_fit <- fit_bundle$raw_fit
n_target <- fit_bundle$configuration$n_target_units
if (n_target != 1177L || nrow(raw_fit$L_matrix) != 2L * n_target ||
    ncol(raw_fit$L_matrix) != length(raw_fit$psd_grid) ||
    !identical(raw_fit$psd_grid[1L], 0)) {
  stop("The saved merged likelihood matrix is not the expected version.")
}

set.seed(sample_seed)
sampled_indices <- sort(sample.int(n_target, n_sample, replace = FALSE))
sampled_pairs <- fit_bundle$selection$pair_key[sampled_indices]
groups <- c("Selected original", "Matched permuted control")
group_colors <- c(
  `Selected original` = "#D55E00",
  `Matched permuted control` = "#0072B2"
)

make_profiles <- function(indices, group) {
  likelihood <- raw_fit$L_matrix[indices, , drop = FALSE]
  relative_likelihood <- sweep(
    likelihood,
    1L,
    likelihood[, 1L],
    `-`
  )
  data.frame(
    matched_pair_index = rep(seq_len(n_sample), each = ncol(likelihood)),
    source_pair_key = rep(sampled_pairs, each = ncol(likelihood)),
    group = group,
    sigma = rep(raw_fit$psd_grid, times = n_sample),
    relative_log_likelihood = as.vector(t(relative_likelihood)),
    stringsAsFactors = FALSE
  )
}

profile_data <- rbind(
  make_profiles(sampled_indices, groups[1L]),
  make_profiles(n_target + sampled_indices, groups[2L])
)
profile_data$group <- factor(profile_data$group, levels = groups)

probabilities <- c(0.10, 0.50, 0.90)
summary_data <- do.call(rbind, lapply(groups, function(group) {
  group_data <- profile_data[profile_data$group == group, ]
  rows <- lapply(split(group_data, group_data$sigma), function(values) {
    quantiles <- stats::quantile(
      values$relative_log_likelihood,
      probs = probabilities,
      names = FALSE,
      type = 8
    )
    data.frame(
      group = group,
      sigma = values$sigma[1L],
      lower = quantiles[1L],
      median = quantiles[2L],
      upper = quantiles[3L],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}))
rownames(summary_data) <- NULL
summary_data$group <- factor(summary_data$group, levels = groups)

sample_table <- data.frame(
  matched_pair_index = seq_len(n_sample),
  source_index = sampled_indices,
  source_pair_key = sampled_pairs,
  sample_seed = sample_seed,
  stringsAsFactors = FALSE
)
utils::write.csv(
  sample_table,
  file.path(output_directory, "matched_likelihood_profile_sample.csv"),
  row.names = FALSE
)
utils::write.csv(
  profile_data,
  file.path(output_directory, "matched_likelihood_profiles.csv"),
  row.names = FALSE
)
utils::write.csv(
  summary_data,
  file.path(output_directory, "matched_likelihood_profile_summary.csv"),
  row.names = FALSE
)

figure <- ggplot2::ggplot(
  profile_data,
  ggplot2::aes(
    x = sigma,
    y = relative_log_likelihood,
    group = source_pair_key,
    color = group
  )
) +
  ggplot2::geom_hline(
    yintercept = 0,
    color = "#666666",
    linewidth = 0.5,
    linetype = "33"
  ) +
  ggplot2::geom_line(linewidth = 0.38, alpha = 0.16) +
  ggplot2::geom_ribbon(
    data = summary_data,
    ggplot2::aes(
      x = sigma,
      ymin = lower,
      ymax = upper,
      fill = group
    ),
    inherit.aes = FALSE,
    alpha = 0.16,
    color = NA
  ) +
  ggplot2::geom_line(
    data = summary_data,
    ggplot2::aes(x = sigma, y = median, color = group),
    inherit.aes = FALSE,
    linewidth = 1.25
  ) +
  ggplot2::facet_wrap(~group, nrow = 1L) +
  ggplot2::scale_color_manual(values = group_colors, guide = "none") +
  ggplot2::scale_fill_manual(values = group_colors, guide = "none") +
  ggplot2::scale_x_continuous(
    trans = scales::pseudo_log_trans(base = 10, sigma = 0.0025),
    breaks = c(0, 0.01, 0.03, 0.10, 0.30, 1),
    labels = c("0", "0.01", "0.03", "0.10", "0.30", "1.00"),
    expand = ggplot2::expansion(mult = c(0.01, 0.02))
  ) +
  ggplot2::scale_y_continuous(
    limits = c(-30, 30),
    breaks = seq(-30, 30, by = 10),
    expand = ggplot2::expansion(mult = c(0.02, 0.02))
  ) +
  ggplot2::labs(
    title = "Likelihood profiles separate selected originals from permuted controls",
    subtitle = paste(
      "The same 100 source pairs are shown in original and permuted form;",
      "thick line = median, band = middle 80%"
    ),
    x = expression(paste("IWP process-SD grid, ", sigma)),
    y = expression(Delta * " log likelihood  " *
      (ell[j](sigma) - ell[j](0))),
    caption = paste(
      "Matched-pair sample seed 20260811.",
      "All profiles are centered at the exact-null component (sigma = 0)."
    )
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 14,
      color = "#111111",
      margin = ggplot2::margin(b = 5)
    ),
    plot.subtitle = ggplot2::element_text(
      size = 10.5,
      color = "#4D4D4D",
      margin = ggplot2::margin(b = 10)
    ),
    plot.caption = ggplot2::element_text(
      size = 9,
      color = "#5A5A5A",
      hjust = 0,
      margin = ggplot2::margin(t = 10)
    ),
    strip.background = ggplot2::element_rect(
      fill = "#F3F3F3",
      color = NA
    ),
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 11.5,
      color = "#222222",
      margin = ggplot2::margin(7, 7, 7, 7)
    ),
    axis.title = ggplot2::element_text(size = 11.5, color = "#222222"),
    axis.text = ggplot2::element_text(size = 10, color = "#333333"),
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
    panel.grid.major.y = ggplot2::element_line(
      color = "#E6E6E6",
      linewidth = 0.4
    ),
    panel.grid.minor = ggplot2::element_blank(),
    panel.spacing.x = grid::unit(1.2, "lines"),
    plot.margin = ggplot2::margin(14, 18, 12, 14)
  )

ggplot2::ggsave(
  file.path(figure_directory, "matched_original_permuted_likelihood_profiles.png"),
  figure,
  width = 11.5,
  height = 6.3,
  dpi = 320,
  bg = "white"
)
ggplot2::ggsave(
  file.path(figure_directory, "matched_original_permuted_likelihood_profiles.pdf"),
  figure,
  width = 11.5,
  height = 6.3,
  device = grDevices::cairo_pdf,
  bg = "white"
)

cat("Matched original/permuted likelihood-profile figure created.\n")
