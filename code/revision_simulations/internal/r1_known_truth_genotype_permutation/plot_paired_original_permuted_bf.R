#!/usr/bin/env Rscript

# Plot paired original and genotype-permuted Bayes factors from the five-seed
# R1 known-truth permutation experiment.

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
    !requireNamespace("scales", quietly = TRUE) ||
    !requireNamespace("patchwork", quietly = TRUE)) {
  stop("The ggplot2, scales, and patchwork packages are required.")
}

workflowr_root <- find_workflowr_root()
experiment_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r1_known_truth_genotype_permutation_mc5_J200_v1"
)
unit_path <- file.path(
  experiment_directory, "summary", "replicate_unit_results.csv"
)
summary_directory <- file.path(experiment_directory, "summary")
figure_directory <- file.path(experiment_directory, "figures")
if (!file.exists(unit_path)) {
  stop("The five-seed unit-level result table is missing.")
}
dir.create(figure_directory, showWarnings = FALSE)

units <- utils::read.csv(unit_path, stringsAsFactors = FALSE)
required_columns <- c(
  "seed", "arm", "fit_stage", "unit_key", "true_null", "bayes_factor"
)
missing_columns <- setdiff(required_columns, names(units))
if (length(missing_columns) > 0L) {
  stop("The unit-level result table is missing required columns.")
}
units <- units[
  units$arm == "shared_genotype_permutation" & units$fit_stage == "BF",
  , drop = FALSE
]
is_permuted <- grepl(
  "__permuted_null_seed", units$unit_key, fixed = TRUE
)
original <- units[!is_permuted, c("seed", "unit_key", "bayes_factor")]
names(original) <- c("seed", "source_unit_id", "original_bf")
permuted <- units[is_permuted, c("seed", "unit_key", "bayes_factor")]
permuted$source_unit_id <- sub(
  "__permuted_null_seed.*$", "", permuted$unit_key
)
names(permuted)[names(permuted) == "bayes_factor"] <- "permuted_bf"
paired <- merge(
  original,
  permuted[, c("seed", "source_unit_id", "permuted_bf")],
  by = c("seed", "source_unit_id"),
  all = FALSE,
  sort = TRUE
)
if (nrow(paired) != 1000L ||
    anyDuplicated(paired[c("seed", "source_unit_id")]) ||
    !setequal(unique(paired$seed), c(12345L, 22345L, 32345L, 42345L, 52345L)) ||
    any(!is.finite(paired$original_bf)) || any(paired$original_bf <= 0) ||
    any(!is.finite(paired$permuted_bf)) || any(paired$permuted_bf <= 0)) {
  stop("The original-permuted BF pairing failed validation.")
}
paired$log10_original_bf <- log10(paired$original_bf)
paired$log10_permuted_bf <- log10(paired$permuted_bf)
paired$permuted_bf_above_one <- paired$permuted_bf > 1
paired$permuted_tail <- factor(
  ifelse(paired$permuted_bf_above_one, "Permuted BF > 1", "Permuted BF <= 1"),
  levels = c("Permuted BF <= 1", "Permuted BF > 1")
)

pearson_correlation <- stats::cor(
  paired$log10_original_bf, paired$log10_permuted_bf
)
spearman_correlation <- stats::cor(
  paired$log10_original_bf,
  paired$log10_permuted_bf,
  method = "spearman"
)
decile_breaks <- stats::quantile(
  paired$log10_original_bf,
  probs = seq(0, 1, by = 0.1),
  names = FALSE,
  type = 7
)
if (anyDuplicated(decile_breaks)) {
  stop("Original log-BF deciles are not uniquely defined.")
}
paired$original_bf_decile <- cut(
  paired$log10_original_bf,
  breaks = decile_breaks,
  include.lowest = TRUE,
  labels = FALSE
)

decile_groups <- split(paired, paired$original_bf_decile)
decile_summary <- do.call(rbind, lapply(decile_groups, function(rows) {
  successes <- sum(rows$permuted_bf_above_one)
  interval <- stats::binom.test(
    successes, nrow(rows), conf.level = 0.95
  )$conf.int
  data.frame(
    original_bf_decile = rows$original_bf_decile[1L],
    n_pairs = nrow(rows),
    median_log10_original_bf = stats::median(rows$log10_original_bf),
    median_log10_permuted_bf = stats::median(rows$log10_permuted_bf),
    q90_log10_permuted_bf = unname(stats::quantile(
      rows$log10_permuted_bf, 0.90, type = 8
    )),
    proportion_permuted_bf_above_one = successes / nrow(rows),
    proportion_ci_lower = interval[1L],
    proportion_ci_upper = interval[2L],
    stringsAsFactors = FALSE
  )
}))
rownames(decile_summary) <- NULL
decile_summary$decile_label <- factor(
  paste0("D", decile_summary$original_bf_decile),
  levels = paste0("D", seq_len(10L))
)

utils::write.csv(
  paired,
  file.path(summary_directory, "paired_original_permuted_bf.csv"),
  row.names = FALSE
)
utils::write.csv(
  decile_summary,
  file.path(summary_directory, "paired_bf_decile_summary.csv"),
  row.names = FALSE
)

scatter <- ggplot2::ggplot(
  paired,
  ggplot2::aes(x = log10_original_bf, y = log10_permuted_bf)
) +
  ggplot2::geom_hline(
    yintercept = 0,
    color = "#777777",
    linewidth = 0.65,
    linetype = "dashed"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(color = permuted_tail),
    alpha = 0.48,
    size = 1.45
  ) +
  ggplot2::geom_smooth(
    method = "loess",
    formula = y ~ x,
    se = TRUE,
    span = 0.8,
    linewidth = 0.9,
    color = "#0072B2",
    fill = "#56B4E9",
    alpha = 0.16
  ) +
  ggplot2::annotate(
    "label",
    x = 0,
    y = 6.35,
    hjust = 0,
    vjust = 1,
    size = 3.45,
    linewidth = 0.2,
    label = sprintf(
      "Pearson r = %.3f\nSpearman rho = %.3f",
      pearson_correlation,
      spearman_correlation
    ),
    color = "#333333",
    fill = "white"
  ) +
  ggplot2::scale_color_manual(
    values = c(
      `Permuted BF <= 1` = "#A6A6A6",
      `Permuted BF > 1` = "#D55E00"
    ),
    name = NULL
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 30, by = 5),
    expand = ggplot2::expansion(mult = c(0.02, 0.03))
  ) +
  ggplot2::scale_y_continuous(
    breaks = seq(-2, 6, by = 2),
    expand = ggplot2::expansion(mult = c(0.03, 0.05))
  ) +
  ggplot2::labs(
    title = "A. Paired log-BF has a sparse upper tail, not a tight line",
    subtitle = "Each point pairs one known alternative with its genotype-permuted copy.",
    x = expression(log[10](BF[original])),
    y = expression(log[10](BF[permuted]))
  ) +
  ggplot2::theme_classic(base_size = 11.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12.5),
    plot.subtitle = ggplot2::element_text(size = 9.5, color = "#555555"),
    legend.position = "top",
    legend.text = ggplot2::element_text(size = 9.5),
    panel.grid.major.y = ggplot2::element_line(
      color = "#E7E7E7", linewidth = 0.35
    ),
    panel.grid.minor = ggplot2::element_blank()
  )

tail_plot <- ggplot2::ggplot(
  decile_summary,
  ggplot2::aes(
    x = decile_label,
    y = proportion_permuted_bf_above_one,
    group = 1
  )
) +
  ggplot2::geom_line(color = "#D55E00", linewidth = 0.85) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = proportion_ci_lower, ymax = proportion_ci_upper),
    color = "#D55E00",
    width = 0.16,
    linewidth = 0.65
  ) +
  ggplot2::geom_point(
    color = "#D55E00",
    fill = "white",
    shape = 21,
    stroke = 0.9,
    size = 2.8
  ) +
  ggplot2::annotate(
    "text",
    x = 1,
    y = decile_summary$proportion_permuted_bf_above_one[1L] + 0.025,
    label = scales::percent(
      decile_summary$proportion_permuted_bf_above_one[1L], accuracy = 1
    ),
    size = 3.2,
    color = "#7A3B00"
  ) +
  ggplot2::annotate(
    "text",
    x = 10,
    y = decile_summary$proportion_permuted_bf_above_one[10L] + 0.025,
    label = scales::percent(
      decile_summary$proportion_permuted_bf_above_one[10L], accuracy = 1
    ),
    size = 3.2,
    color = "#7A3B00"
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 0.36),
    breaks = seq(0, 0.30, by = 0.10),
    labels = scales::label_percent(accuracy = 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    title = "B. The leakage tail is enriched among stronger originals",
    subtitle = "Exact 95% binomial intervals; 100 pairs per original-BF decile.",
    x = "Decile of original log-BF (weakest to strongest)",
    y = "Permuted copies with BF > 1"
  ) +
  ggplot2::theme_classic(base_size = 11.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12.5),
    plot.subtitle = ggplot2::element_text(size = 9.5, color = "#555555"),
    panel.grid.major.y = ggplot2::element_line(
      color = "#E7E7E7", linewidth = 0.35
    ),
    panel.grid.minor = ggplot2::element_blank()
  )

figure <- scatter + tail_plot +
  patchwork::plot_layout(widths = c(1.45, 1)) +
  patchwork::plot_annotation(
    title = "Original signal strength predicts only the upper tail of permuted evidence",
    subtitle = paste(
      "Five R1 source seeds; BFs are paired within the same BF-adjusted",
      "400-unit genotype-permutation fit."
    ),
    caption = paste(
      "The weak marginal correlation is expected when unit-specific alignment",
      "between residualized G and PG dominates the leakage magnitude."
    ),
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = 15, color = "#111111"
      ),
      plot.subtitle = ggplot2::element_text(
        size = 10.5, color = "#4D4D4D"
      ),
      plot.caption = ggplot2::element_text(
        size = 9, color = "#5A5A5A", hjust = 0
      )
    )
  )

png_path <- file.path(
  figure_directory, "paired_original_vs_permuted_log_bf.png"
)
pdf_path <- file.path(
  figure_directory, "paired_original_vs_permuted_log_bf.pdf"
)
ggplot2::ggsave(
  png_path,
  figure,
  width = 13.2,
  height = 6.8,
  dpi = 320,
  bg = "white"
)
ggplot2::ggsave(
  pdf_path,
  figure,
  width = 13.2,
  height = 6.8,
  device = grDevices::cairo_pdf,
  bg = "white"
)

summary_row <- data.frame(
  n_pairs = nrow(paired),
  pearson_log10_bf = pearson_correlation,
  spearman_log10_bf = spearman_correlation,
  lowest_decile_proportion_permuted_bf_above_one =
    decile_summary$proportion_permuted_bf_above_one[1L],
  highest_decile_proportion_permuted_bf_above_one =
    decile_summary$proportion_permuted_bf_above_one[10L],
  stringsAsFactors = FALSE
)
utils::write.csv(
  summary_row,
  file.path(summary_directory, "paired_bf_summary.csv"),
  row.names = FALSE
)

print(summary_row)
cat("Paired original-versus-permuted BF figure created.\n")
