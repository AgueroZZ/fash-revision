#!/usr/bin/env Rscript

# Side-by-side boxplots of the competition statistic W for published
# discoveries and non-discoveries.
#
# W = log BF(target) - log BF(decoy), both scored by the same merged prior.
# W > 0 means the observed genotype carries more evidence for a non-constant
# effect than its permuted decoy does.
#
# IMPORTANT: the grouping conditions on the target. The published rule selects
# pairs with a high target Bayes factor, and a high target BF mechanically
# implies W > 0 unless the decoy also happens to be high. So the separation
# between these two boxes is NOT evidence that the discoveries are real; a null
# pair that reached a high BF by chance would land in the same box. The figure
# describes where the competition statistic sits, it does not validate the
# calls.

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

for (package in c("ggplot2", "patchwork")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The ", package, " package is required.")
  }
}
library(ggplot2)

workflowr_root <- find_workflowr_root()
statistics <- readRDS(file.path(
  workflowr_root, "output", "dynamic_eQTL_real",
  "fash_knockoff_full_seed20260823", "pair_statistics.rds"
))
figure_directory <- file.path(
  workflowr_root, "output", "revision_simulations", "internal",
  "fash_knockoff_full_seed20260823_from_midway3", "figures"
)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)

# Reconstruct the published cumulative-lfdr rule at alpha 0.05.
alpha <- 0.05
ordering <- order(statistics$source_bf_lfdr)
cumulative <- cumsum(statistics$source_bf_lfdr[ordering]) /
  seq_along(ordering)
n_call <- sum(cumulative <= alpha)
discovered <- rep(FALSE, nrow(statistics))
discovered[ordering[seq_len(n_call)]] <- TRUE
if (n_call != 9214L || length(unique(statistics$gene_id[discovered])) != 1176L) {
  stop("The reconstructed discovery set does not match the published counts.")
}

plot_data <- data.frame(
  group = factor(
    ifelse(discovered, "discovered", "not discovered"),
    levels = c("discovered", "not discovered")
  ),
  W = statistics$W_updated_weights
)

group_summary <- do.call(rbind, lapply(levels(plot_data$group), function(g) {
  values <- plot_data$W[plot_data$group == g]
  data.frame(
    group = g,
    n = length(values),
    min = min(values),
    q25 = quantile(values, 0.25, names = FALSE),
    median = median(values),
    q75 = quantile(values, 0.75, names = FALSE),
    max = max(values),
    frac_positive = mean(values > 0),
    stringsAsFactors = FALSE
  )
}))
print(group_summary, row.names = FALSE, digits = 4)
utils::write.csv(
  group_summary, file.path(figure_directory, "w_by_discovery_summary.csv"),
  row.names = FALSE
)

group_labels <- setNames(
  sprintf("%s\nn = %s", group_summary$group,
          formatC(group_summary$n, big.mark = ",", format = "d")),
  group_summary$group
)
group_colours <- c("discovered" = "#1f6f8b", "not discovered" = "#9aa5ab")

base_theme <- theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9),
    axis.text.x = element_text(size = 10)
  )

# Panel a draws the outlying points, because the extreme tails on BOTH sides
# are the substantive content: the decoy reaches -28 while the target reaches
# +33, and knockoff+ operates only out there.
full_panel <- ggplot(plot_data, aes(x = group, y = W, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey35",
             linewidth = 0.4) +
  geom_jitter(aes(colour = group), width = 0.28, height = 0,
              size = 0.15, alpha = 0.04) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.85,
               colour = "grey15", linewidth = 0.4) +
  scale_fill_manual(values = group_colours) +
  scale_colour_manual(values = group_colours) +
  scale_x_discrete(labels = group_labels) +
  labs(
    x = NULL, y = expression(W[j] == log ~ BF[target] - log ~ BF[decoy]),
    title = "a  Full range, outlying points shown",
    subtitle = sprintf(
      "range: discovered %.1f to %.1f;  not discovered %.1f to %.1f",
      group_summary$min[1], group_summary$max[1],
      group_summary$min[2], group_summary$max[2]
    )
  ) +
  base_theme

zoom_panel <- ggplot(plot_data, aes(x = group, y = W, fill = group)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey35",
             linewidth = 0.4) +
  geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.85,
               colour = "grey15", linewidth = 0.4) +
  scale_fill_manual(values = group_colours) +
  scale_x_discrete(labels = group_labels) +
  coord_cartesian(ylim = c(-4, 13)) +
  labs(
    x = NULL, y = NULL,
    title = "b  Quartiles, zoomed",
    subtitle = sprintf(
      "medians %.2f vs %.5f;  W > 0 in %.2f%% vs %.2f%%",
      group_summary$median[1], group_summary$median[2],
      100 * group_summary$frac_positive[1],
      100 * group_summary$frac_positive[2]
    )
  ) +
  base_theme

combined <- (full_panel | zoom_panel) +
  patchwork::plot_annotation(
    caption = paste(
      "Boxes are quartiles with 1.5 IQR whiskers. Points in panel a are jittered at alpha 0.04.",
      "The grouping conditions on the target: the published rule selects pairs with a high",
      "\ntarget Bayes factor, and a high target BF mechanically forces W > 0 unless the decoy is also high.",
      "A null pair that reached a high BF by chance would sit in the same box.",
      "\nThe figure locates W by discovery status; it does not validate the calls.",
      sep = " "
    ),
    theme = theme(plot.caption = element_text(hjust = 0, size = 8,
                                             colour = "grey25"))
  )
ggsave(file.path(figure_directory, "w_boxplot_by_discovery.png"),
       combined, width = 9.5, height = 5.8, dpi = 200)
ggsave(file.path(figure_directory, "w_boxplot_by_discovery.pdf"),
       combined, width = 9.5, height = 5.8)
message("Wrote figures to ", figure_directory)
