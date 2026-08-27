# Canonical reporting helpers for the full-data FASH target-decoy analysis.
#
# Everything on the internal workflowr page is produced here, so the page, the
# standalone plot scripts, and any later re-analysis share one implementation.

if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("patchwork", quietly = TRUE)) {
  stop("The ggplot2 and patchwork packages are required.")
}

FASH_KNOCKOFF_RUN_ID <- "fash_knockoff_full_seed20260823"
FASH_KNOCKOFF_ALPHA <- 0.05
FASH_KNOCKOFF_ARM_COLOURS <- c(
  "target (observed genotype)" = "#1f6f8b",
  "decoy (permuted genotype)" = "#c1462f"
)
FASH_KNOCKOFF_GROUP_COLOURS <- c(
  "discovered" = "#1f6f8b", "not discovered" = "#9aa5ab"
)

fash_knockoff_paths <- function(workflowr_root = ".") {
  list(
    bulk = file.path(
      workflowr_root, "output", "dynamic_eQTL_real", FASH_KNOCKOFF_RUN_ID
    ),
    summary = file.path(
      workflowr_root, "output", "revision_simulations", "internal",
      paste0(FASH_KNOCKOFF_RUN_ID, "_from_midway3")
    )
  )
}

# Loads the pair-level statistics and reconstructs the published discovery set
# from the cumulative-lfdr rule, checking it against the known counts.
load_fash_knockoff_statistics <- function(workflowr_root = ".") {
  paths <- fash_knockoff_paths(workflowr_root)
  statistics <- readRDS(file.path(paths$bulk, "pair_statistics.rds"))
  ordering <- order(statistics$source_bf_lfdr)
  cumulative <- cumsum(statistics$source_bf_lfdr[ordering]) /
    seq_along(ordering)
  n_call <- sum(cumulative <= FASH_KNOCKOFF_ALPHA)
  statistics$discovered <- FALSE
  statistics$discovered[ordering[seq_len(n_call)]] <- TRUE
  if (n_call != 9214L ||
      length(unique(statistics$gene_id[statistics$discovered])) != 1176L) {
    stop("The reconstructed discovery set does not match the published counts.")
  }
  statistics$group <- factor(
    ifelse(statistics$discovered, "discovered", "not discovered"),
    levels = c("discovered", "not discovered")
  )
  statistics
}

fash_knockoff_theme <- function() {
  ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = 11),
      plot.subtitle = ggplot2::element_text(size = 9),
      legend.title = ggplot2::element_blank()
    )
}

fash_knockoff_group_summary <- function(statistics) {
  do.call(rbind, lapply(levels(statistics$group), function(g) {
    values <- statistics$W_updated_weights[statistics$group == g]
    data.frame(
      group = g, n = length(values), min = min(values),
      q25 = stats::quantile(values, 0.25, names = FALSE),
      median = stats::median(values),
      q75 = stats::quantile(values, 0.75, names = FALSE),
      max = max(values), frac_positive = mean(values > 0),
      stringsAsFactors = FALSE
    )
  }))
}

# The headline figure. Panel a shows the tails, which is where knockoff+ acts.
plot_w_by_discovery <- function(statistics) {
  library(ggplot2)
  group_summary <- fash_knockoff_group_summary(statistics)
  labels <- setNames(
    sprintf("%s\nn = %s", group_summary$group,
            formatC(group_summary$n, big.mark = ",", format = "d")),
    group_summary$group
  )
  full <- ggplot(statistics, aes(x = group, y = W_updated_weights,
                                 fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey35",
               linewidth = 0.4) +
    geom_jitter(aes(colour = group), width = 0.28, height = 0, size = 0.15,
                alpha = 0.04) +
    geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.85,
                 colour = "grey15", linewidth = 0.4) +
    scale_fill_manual(values = FASH_KNOCKOFF_GROUP_COLOURS) +
    scale_colour_manual(values = FASH_KNOCKOFF_GROUP_COLOURS) +
    scale_x_discrete(labels = labels) +
    labs(
      x = NULL, y = expression(W[j] == log ~ BF[target] - log ~ BF[decoy]),
      title = "a  Full range, outlying points shown",
      subtitle = sprintf(
        "range: discovered %.1f to %.1f;  not discovered %.1f to %.1f",
        group_summary$min[1], group_summary$max[1],
        group_summary$min[2], group_summary$max[2]
      )
    ) +
    fash_knockoff_theme() + theme(legend.position = "none")
  zoom <- ggplot(statistics, aes(x = group, y = W_updated_weights,
                                 fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey35",
               linewidth = 0.4) +
    geom_boxplot(width = 0.4, outlier.shape = NA, alpha = 0.85,
                 colour = "grey15", linewidth = 0.4) +
    scale_fill_manual(values = FASH_KNOCKOFF_GROUP_COLOURS) +
    scale_x_discrete(labels = labels) +
    coord_cartesian(ylim = c(-4, 13)) +
    labs(
      x = NULL, y = NULL, title = "b  Quartiles, zoomed",
      subtitle = sprintf(
        "medians %.2f vs %.5f;  W > 0 in %.2f%% vs %.2f%%",
        group_summary$median[1], group_summary$median[2],
        100 * group_summary$frac_positive[1],
        100 * group_summary$frac_positive[2]
      )
    ) +
    fash_knockoff_theme() + theme(legend.position = "none")
  (full | zoom)
}

fash_knockoff_lfdr_long <- function(statistics) {
  n_pair <- nrow(statistics)
  out <- data.frame(
    arm = rep(c("target (observed genotype)", "decoy (permuted genotype)"),
              each = n_pair),
    lfdr = c(statistics$target_lfdr, statistics$decoy_lfdr)
  )
  out$arm <- factor(out$arm, levels = names(FASH_KNOCKOFF_ARM_COLOURS))
  out
}

FASH_KNOCKOFF_LFDR_BREAKS <- c(-Inf, -12, -9, -6, -4, -3, -2, log10(0.05), -1, 0)
FASH_KNOCKOFF_LFDR_LABELS <- c(
  "<1e-12", "1e-12..1e-9", "1e-9..1e-6", "1e-6..1e-4", "1e-4..1e-3",
  "1e-3..0.01", "0.01..0.05", "0.05..0.1", "0.1..1"
)

fash_knockoff_lfdr_bands <- function(statistics) {
  long <- fash_knockoff_lfdr_long(statistics)
  long$bin <- cut(log10(long$lfdr), breaks = FASH_KNOCKOFF_LFDR_BREAKS,
                  labels = FASH_KNOCKOFF_LFDR_LABELS)
  long <- subset(long, !is.na(bin) & bin != "0.1..1")
  long$bin <- droplevels(long$bin)
  counts <- as.data.frame(table(arm = long$arm, bin = long$bin))
  wide <- stats::reshape(counts, idvar = "bin", timevar = "arm",
                         direction = "wide")
  names(wide) <- c("lfdr_band", "target", "decoy")
  wide$decoy_over_target <- round(wide$decoy / wide$target, 3)
  total <- data.frame(
    lfdr_band = "all < 0.05",
    target = sum(statistics$target_lfdr < 0.05),
    decoy = sum(statistics$decoy_lfdr < 0.05)
  )
  total$decoy_over_target <- round(total$decoy / total$target, 3)
  rbind(wide[wide$lfdr_band != "0.05..0.1", ], total,
        wide[wide$lfdr_band == "0.05..0.1", ])
}

plot_lfdr_bands <- function(statistics) {
  library(ggplot2)
  bands <- fash_knockoff_lfdr_bands(statistics)
  bands <- bands[bands$lfdr_band != "all < 0.05", ]
  long <- rbind(
    data.frame(bin = bands$lfdr_band, arm = names(FASH_KNOCKOFF_ARM_COLOURS)[1],
               count = bands$target),
    data.frame(bin = bands$lfdr_band, arm = names(FASH_KNOCKOFF_ARM_COLOURS)[2],
               count = bands$decoy)
  )
  long$bin <- factor(long$bin, levels = bands$lfdr_band)
  long$arm <- factor(long$arm, levels = names(FASH_KNOCKOFF_ARM_COLOURS))
  ggplot(long, aes(x = bin, y = count, fill = arm)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7,
             alpha = 0.6) +
    geom_text(aes(label = count), position = position_dodge(width = 0.75),
              vjust = -0.35, size = 2.7) +
    scale_fill_manual(values = FASH_KNOCKOFF_ARM_COLOURS) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      x = "lfdr band (merged fit)", y = "gene-variant pairs",
      title = "Left tail of the lfdr distribution, counted",
      subtitle = sprintf(
        "lfdr < 0.05: %s target vs %s decoy (ratio %.3f)",
        format(sum(statistics$target_lfdr < 0.05), big.mark = ","),
        format(sum(statistics$decoy_lfdr < 0.05), big.mark = ","),
        sum(statistics$decoy_lfdr < 0.05) / sum(statistics$target_lfdr < 0.05)
      )
    ) +
    fash_knockoff_theme() +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 30, hjust = 1))
}

# Genes are treated as independent clusters; no within-gene correlation is
# assumed. The naive sqrt(n) sign test is reported alongside to show how much
# it overstates the evidence.
fash_knockoff_sign_evidence <- function(statistics) {
  signs <- sign(statistics$W_updated_weights)
  keep <- signs != 0
  gene <- factor(statistics$gene_id[keep])
  per_gene <- tapply(signs[keep], gene, sum)
  sizes <- tapply(signs[keep], gene, length)
  total <- sum(per_gene)
  n <- sum(sizes)
  mean_sign <- total / n
  robust_variance <- sum((per_gene - sizes * mean_sign)^2)
  design_effect <- robust_variance / n
  data.frame(
    quantity = c("pairs", "genes", "target-win excess",
                 "naive sign-test z (independent pairs)",
                 "cluster-robust z (genes as clusters)",
                 "design effect", "implied within-gene sign correlation"),
    value = c(n, length(per_gene), total, total / sqrt(n),
              total / sqrt(robust_variance), design_effect,
              (design_effect - 1) / (mean(sizes) - 1)),
    stringsAsFactors = FALSE
  )
}

fash_knockoff_selection_table <- function(statistics, q_grid) {
  target <- log(statistics$target_bf_updated_weights)
  decoy <- log(statistics$decoy_bf_updated_weights)
  gene <- factor(statistics$gene_id)
  gene_W <- as.numeric(tapply(target, gene, max) - tapply(decoy, gene, max))
  rbind(
    data.frame(unit = "gene-variant pair",
               knockoff_selection_summary(statistics$W_updated_weights, q_grid)),
    data.frame(unit = "gene (max vs max)",
               knockoff_selection_summary(gene_W, q_grid))
  )
}

fash_knockoff_gene_statistic_comparison <- function(statistics, q_grid) {
  target <- log(statistics$target_bf_updated_weights)
  decoy <- log(statistics$decoy_bf_updated_weights)
  gene <- factor(statistics$gene_id)
  top_k <- function(x, k) sum(sort(x, decreasing = TRUE)[seq_len(min(k, length(x)))])
  candidates <- list(
    "max - max" = tapply(target, gene, max) - tapply(decoy, gene, max),
    "sum - sum" = tapply(target, gene, sum) - tapply(decoy, gene, sum),
    "mean - mean" = tapply(target, gene, mean) - tapply(decoy, gene, mean),
    "top5 - top5" = tapply(target, gene, top_k, 5) -
      tapply(decoy, gene, top_k, 5),
    "top20 - top20" = tapply(target, gene, top_k, 20) -
      tapply(decoy, gene, top_k, 20)
  )
  do.call(rbind, lapply(names(candidates), function(nm) {
    x <- as.numeric(candidates[[nm]])
    path <- knockoff_estimated_fdr_path(x)
    selection <- knockoff_selection_summary(x, q_grid)
    data.frame(
      statistic = nm, frac_positive = mean(x > 0),
      t_statistic = mean(x) / (stats::sd(x) / sqrt(length(x))),
      min_estimated_fdr = min(path$estimated_fdr),
      discoveries_q20 = selection$n_discoveries[
        which.min(abs(selection$q - 0.20))
      ],
      stringsAsFactors = FALSE
    )
  }))
}

fash_knockoff_tail_concentration <- function(statistics, thresholds) {
  W <- statistics$W_updated_weights
  gene <- as.character(statistics$gene_id)
  do.call(rbind, lapply(thresholds, function(t) {
    target <- W >= t
    decoy <- W <= -t
    data.frame(
      threshold = t,
      target_pairs = sum(target), target_genes = length(unique(gene[target])),
      decoy_pairs = sum(decoy), decoy_genes = length(unique(gene[decoy])),
      pair_ratio = (1 + sum(decoy)) / max(sum(target), 1L),
      gene_ratio = (1 + length(unique(gene[decoy]))) /
        max(length(unique(gene[target])), 1L),
      stringsAsFactors = FALSE
    )
  }))
}

fash_knockoff_selected_pairs <- function(statistics, q) {
  W <- statistics$W_updated_weights
  threshold <- knockoff_plus_threshold(W, q)
  if (!is.finite(threshold)) {
    return(statistics[0, ])
  }
  selected <- statistics[W >= threshold, ]
  out <- data.frame(
    gene_id = selected$gene_id,
    variant_id = selected$variant_id,
    published_lfdr = signif(selected$source_bf_lfdr, 3),
    log10_target_BF = round(log10(selected$target_bf_updated_weights), 2),
    log10_decoy_BF = round(log10(selected$decoy_bf_updated_weights), 2),
    W = round(selected$W_updated_weights, 2),
    stringsAsFactors = FALSE
  )
  out[order(-out$W), ]
}
