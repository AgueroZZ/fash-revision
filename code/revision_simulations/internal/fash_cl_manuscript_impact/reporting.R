# Reporting helpers for the internal current-FASH versus FASH-CL page.

find_workflowr_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "_workflowr.yml"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the workflowr project root.")
    }
    current <- parent
  }
}

require_reporting_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The ", package, " package is required to render this page.")
  }
  invisible(TRUE)
}

require_reporting_package("ggplot2")
require_reporting_package("ggVennDiagram")
require_reporting_package("ggrastr")
require_reporting_package("knitr")
require_reporting_package("patchwork")
require_reporting_package("scales")

workflowr_root <- find_workflowr_root()
cache_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fash_cl_manuscript_impact",
  "analysis_cache.rds"
)
if (!file.exists(cache_path)) {
  stop(
    "The retained manuscript-impact cache is missing. Run ",
    "code/revision_simulations/internal/fash_cl_manuscript_impact/",
    "run_fash_cl_manuscript_impact.R first."
  )
}
cache <- readRDS(cache_path)
required_fields <- c(
  "configuration", "input_provenance", "discovery_summary",
  "method_overlap", "lfdr_scatter_all", "lfdr_scatter_top",
  "lfdr_comparison_summary", "prior_weight_comparison",
  "prior_weight_summary", "venn_sets", "venn_regions", "category_counts",
  "current_classification", "cl_classification", "figure3_status",
  "figure3_proposed", "figure5_status", "figure5_proposed",
  "figure6_status", "figure6_proposed", "hallmark_comparison",
  "manuscript_hallmark_reference", "posterior_plot_data",
  "observed_plot_data", "parametric_plot_data", "claim_impact",
  "runtime_summary"
)
if (!all(required_fields %in% names(cache)) ||
    !identical(
      cache$configuration$analysis_id,
      "revision_internal_fash_cl_manuscript_impact"
    ) ||
    cache$configuration$alpha != 0.05 ||
    cache$configuration$posterior_sample_size != 3000L ||
    cache$configuration$switch_threshold != 0.25) {
  stop("The retained manuscript-impact cache failed structural validation.")
}

current_md5 <- unname(tools::md5sum(cache$input_provenance$path))
if (anyNA(current_md5) ||
    !identical(current_md5, cache$input_provenance$md5)) {
  stop("At least one retained input changed after cache creation.")
}

configuration <- cache$configuration
discovery_summary <- cache$discovery_summary
method_overlap <- cache$method_overlap
lfdr_scatter_all <- cache$lfdr_scatter_all
lfdr_scatter_top <- cache$lfdr_scatter_top
lfdr_comparison_summary <- cache$lfdr_comparison_summary
prior_weight_comparison <- cache$prior_weight_comparison
prior_weight_summary <- cache$prior_weight_summary
venn_sets <- cache$venn_sets
venn_regions <- cache$venn_regions
category_counts <- cache$category_counts
current_classification <- cache$current_classification
cl_classification <- cache$cl_classification
figure3_status <- cache$figure3_status
figure3_proposed <- cache$figure3_proposed
figure5_status <- cache$figure5_status
figure5_proposed <- cache$figure5_proposed
figure6_status <- cache$figure6_status
figure6_proposed <- cache$figure6_proposed
hallmark_comparison <- cache$hallmark_comparison
manuscript_hallmark_reference <- cache$manuscript_hallmark_reference
posterior_plot_data <- cache$posterior_plot_data
observed_plot_data <- cache$observed_plot_data
parametric_plot_data <- cache$parametric_plot_data
claim_impact <- cache$claim_impact
runtime_summary <- cache$runtime_summary
input_provenance <- cache$input_provenance

expected_discovery <- data.frame(
  method = rep(c("Current FASH", "FASH-CL"), each = 2L),
  test = rep(c("IWP1 dynamic", "IWP2 nonlinear"), 2L),
  pair_count = c(9205L, 44L, 5393L, 60L),
  gene_count = c(1177L, 9L, 686L, 6L),
  stringsAsFactors = FALSE
)
if (!identical(discovery_summary, expected_discovery) ||
    nrow(cl_classification) != 5393L * 4L ||
    nrow(lfdr_scatter_all) != 1009173L ||
    nrow(lfdr_scatter_top) !=
      lfdr_comparison_summary$pair_count[
        lfdr_comparison_summary$selection ==
          "Current-FASH top pair per gene"
      ] ||
    nrow(lfdr_scatter_top) != length(unique(lfdr_scatter_top$gene_id)) ||
    any(!is.finite(lfdr_scatter_all$current_lfdr)) ||
    any(!is.finite(lfdr_scatter_all$cl_lfdr)) ||
    any(!is.finite(lfdr_scatter_top$current_lfdr)) ||
    any(!is.finite(lfdr_scatter_top$cl_lfdr)) ||
    nrow(prior_weight_summary) != 4L ||
    any(abs(prior_weight_summary$weight_sum - 1) > 1e-7) ||
    !setequal(
      unique(prior_weight_comparison$method),
      c("Current FASH", "FASH-CL")
    ) ||
    !setequal(
      unique(prior_weight_comparison$adjustment),
      c("Raw", "BF-adjusted")
    ) ||
    !identical(
      names(venn_sets),
      c("Current FASH", "FASH-CL")
    ) ||
    !all(vapply(
      venn_sets,
      function(method_sets) identical(
        names(method_sets),
        c("Genes", "Gene-variant pairs", "Variants")
      ),
      logical(1)
    )) ||
    nrow(venn_regions) != 42L ||
    !setequal(
      unique(venn_regions$unit),
      c("Genes", "Gene-variant pairs", "Variants")
    ) ||
    !identical(
      method_overlap$unit,
      c("Gene-variant pairs", "Genes", "Variants")
    ) ||
    !setequal(unique(cl_classification$category), c(
      "early", "middle", "late", "switch"
    )) ||
    nrow(figure3_proposed) != 6L ||
    nrow(figure5_proposed) != 4L ||
    nrow(figure6_proposed) != 4L ||
    !identical(
      as.logical(figure6_proposed$category_significant),
      as.logical(figure6_proposed$cfsr <= configuration$alpha)
    ) ||
    nrow(hallmark_comparison) != 8L ||
    nrow(posterior_plot_data) == 0L ||
    nrow(observed_plot_data) == 0L) {
  stop("The retained manuscript-impact cache failed numerical validation.")
}

format_integer <- function(value) {
  format(as.integer(round(value)), big.mark = ",", scientific = FALSE)
}

format_decimal <- function(value, digits = 3L) {
  formatC(as.numeric(value), format = "f", digits = digits)
}

format_scientific <- function(value, digits = 2L) {
  formatC(as.numeric(value), format = "e", digits = digits)
}

format_percent <- function(value, digits = 1L) {
  paste0(format_decimal(100 * value, digits), "%")
}

render_scrollable_table <- function(data,
                                    digits = 3L,
                                    align = NULL,
                                    minimum_width = "760px") {
  table_html <- knitr::kable(
    data,
    format = "html",
    escape = TRUE,
    align = align,
    row.names = FALSE,
    digits = digits,
    table.attr = paste0(
      'class="table table-striped table-hover" style="min-width:',
      minimum_width,
      ';"'
    )
  )
  cat(
    '<div class="impact-table-scroll">\n',
    as.character(table_html),
    "\n</div>\n",
    sep = ""
  )
  invisible(table_html)
}

merge_method_table <- function(data, key_columns, value_columns) {
  current <- data[data$method == "Current FASH", c(key_columns, value_columns)]
  cl <- data[data$method == "FASH-CL", c(key_columns, value_columns)]
  names(current)[match(value_columns, names(current))] <- paste0(
    "current_", value_columns
  )
  names(cl)[match(value_columns, names(cl))] <- paste0("cl_", value_columns)
  merge(current, cl, by = key_columns, sort = FALSE)
}

discovery_change_table <- merge_method_table(
  discovery_summary,
  key_columns = "test",
  value_columns = c("pair_count", "gene_count")
)
discovery_change_table$pair_change <- 100 * (
  discovery_change_table$cl_pair_count - discovery_change_table$current_pair_count
) / discovery_change_table$current_pair_count
discovery_change_table$gene_change <- 100 * (
  discovery_change_table$cl_gene_count - discovery_change_table$current_gene_count
) / discovery_change_table$current_gene_count
discovery_change_table <- discovery_change_table[, c(
  "test", "current_pair_count", "cl_pair_count", "pair_change",
  "current_gene_count", "cl_gene_count", "gene_change"
)]
names(discovery_change_table) <- c(
  "Test", "Current pairs", "FASH-CL pairs", "Pair change (%)",
  "Current genes", "FASH-CL genes", "Gene change (%)"
)

category_change_table <- merge_method_table(
  category_counts,
  key_columns = "category",
  value_columns = c("pair_count", "gene_count")
)
category_change_table$pair_change <- 100 * (
  category_change_table$cl_pair_count - category_change_table$current_pair_count
) / category_change_table$current_pair_count
category_change_table$gene_change <- 100 * (
  category_change_table$cl_gene_count - category_change_table$current_gene_count
) / category_change_table$current_gene_count
category_order <- c("nonlinear", "early", "middle", "late", "switch")
category_change_table <- category_change_table[
  match(category_order, category_change_table$category),
  ,
  drop = FALSE
]

status_label <- function(significant) {
  ifelse(significant, "Retained", "Lost")
}

figure3_status_table <- data.frame(
  Panel = figure3_status$panel,
  Unit = paste0(figure3_status$gene_symbol, " / ", figure3_status$variant_id),
  Role = figure3_status$role,
  `Current lfdr` = figure3_status$current_lfdr,
  `FASH-CL lfdr` = figure3_status$cl_lfdr,
  `FASH-CL status` = status_label(figure3_status$cl_significant),
  check.names = FALSE
)
figure3_replacement_table <- data.frame(
  Panel = figure3_proposed$panel,
  Unit = paste0(figure3_proposed$gene_symbol, " / ", figure3_proposed$variant_id),
  Role = figure3_proposed$replacement_role,
  lfdr = figure3_proposed$lfdr,
  `Linear p` = figure3_proposed$p_linear,
  `Quadratic p` = figure3_proposed$p_quadratic,
  check.names = FALSE
)

figure5_status_table <- data.frame(
  Panel = figure5_status$panel,
  Unit = paste0(figure5_status$gene_symbol, " / ", figure5_status$variant_id),
  `Current lfdr` = figure5_status$current_lfdr,
  `FASH-CL lfdr` = figure5_status$cl_lfdr,
  `FASH-CL status` = status_label(figure5_status$cl_significant),
  check.names = FALSE
)
figure5_replacement_table <- data.frame(
  Panel = figure5_proposed$panel,
  Unit = paste0(figure5_proposed$gene_symbol, " / ", figure5_proposed$variant_id),
  Status = figure5_proposed$replacement_status,
  lfdr = figure5_proposed$lfdr,
  check.names = FALSE
)

figure6_status_table <- data.frame(
  Panel = figure6_status$panel,
  Category = figure6_status$role,
  Unit = paste0(figure6_status$gene_symbol, " / ", figure6_status$variant_id),
  `Current lfsr` = figure6_status$current_lfsr,
  `FASH-CL lfsr` = figure6_status$cl_lfsr,
  `FASH-CL category status` = status_label(
    !is.na(figure6_status$cl_category_significant) &
      figure6_status$cl_category_significant
  ),
  check.names = FALSE
)
figure6_replacement_table <- data.frame(
  Panel = figure6_proposed$panel,
  Category = figure6_proposed$category,
  Unit = paste0(figure6_proposed$gene_symbol, " / ", figure6_proposed$variant_id),
  Status = figure6_proposed$replacement_status,
  lfsr = figure6_proposed$lfsr,
  cFSR = figure6_proposed$cfsr,
  `Significant at FSR 0.05` = figure6_proposed$category_significant,
  check.names = FALSE
)

term_labels <- c(
  HALLMARK_HYPOXIA = "Hypoxia",
  HALLMARK_KRAS_SIGNALING_UP = "KRAS signaling up"
)
hallmark_display <- hallmark_comparison
hallmark_display$Hallmark <- unname(term_labels[hallmark_display$term])
hallmark_display <- hallmark_display[, c(
  "method", "gene_set", "Hallmark", "overlap_count", "selected_size",
  "p_value", "q_value"
)]
names(hallmark_display) <- c(
  "Method", "Gene set", "Hallmark", "Overlap", "Selected genes",
  "p-value", "BH q-value"
)
manuscript_hallmark_display <- manuscript_hallmark_reference
manuscript_hallmark_display$Hallmark <- unname(
  term_labels[manuscript_hallmark_display$term]
)
manuscript_hallmark_display <- manuscript_hallmark_display[, c(
  "method", "gene_set", "Hallmark", "overlap_count", "selected_size",
  "p_value", "q_value"
)]
names(manuscript_hallmark_display) <- names(hallmark_display)

make_unit_label_table <- function() {
  status_units <- rbind(
    figure3_status[, c("key", "gene_symbol", "variant_id")],
    figure5_status[, c("key", "gene_symbol", "variant_id")],
    figure6_status[, c("key", "gene_symbol", "variant_id")]
  )
  proposed_units <- rbind(
    figure3_proposed[, c("key", "gene_symbol", "variant_id")],
    figure5_proposed[, c("key", "gene_symbol", "variant_id")],
    figure6_proposed[, c("key", "gene_symbol", "variant_id")]
  )
  units <- unique(rbind(status_units, proposed_units))
  units$unit_label <- paste0(units$gene_symbol, " / ", units$variant_id)
  units
}

unit_labels <- make_unit_label_table()
nonsignificant_figure6_keys <- figure6_proposed$key[
  !figure6_proposed$category_significant
]
unit_labels$unit_label[
  unit_labels$key %in% nonsignificant_figure6_keys
] <- paste0(
  unit_labels$unit_label[unit_labels$key %in% nonsignificant_figure6_keys],
  " [not significant]"
)

prepare_unit_plot_data <- function(keys, methods, order_label) {
  posterior <- posterior_plot_data[
    posterior_plot_data$key %in% keys &
      posterior_plot_data$method %in% methods &
      posterior_plot_data$order == order_label,
    ,
    drop = FALSE
  ]
  observed <- observed_plot_data[
    observed_plot_data$key %in% keys &
      observed_plot_data$method %in% methods &
      observed_plot_data$order == order_label,
    ,
    drop = FALSE
  ]
  parametric <- parametric_plot_data[
    parametric_plot_data$key %in% keys &
      parametric_plot_data$method %in% methods &
      parametric_plot_data$order == order_label,
    ,
    drop = FALSE
  ]
  labels <- unit_labels[match(keys, unit_labels$key), , drop = FALSE]
  ordered_facets <- unlist(lapply(methods, function(method) {
    paste(method, labels$unit_label, sep = "\n")
  }))
  add_facet <- function(data) {
    data$unit_label <- unit_labels$unit_label[match(data$key, unit_labels$key)]
    data$facet <- factor(
      paste(data$method, data$unit_label, sep = "\n"),
      levels = ordered_facets
    )
    data
  }
  list(
    posterior = add_facet(posterior),
    observed = add_facet(observed),
    parametric = add_facet(parametric)
  )
}

plot_unit_comparison <- function(keys,
                                 methods,
                                 order_label,
                                 ncol,
                                 show_parametric = TRUE,
                                 include_zero_line = FALSE) {
  data <- prepare_unit_plot_data(keys, methods, order_label)
  plot <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = data$posterior,
      ggplot2::aes(x = time, ymin = lower, ymax = upper),
      fill = "#E64B35",
      alpha = 0.16
    ) +
    ggplot2::geom_line(
      data = data$posterior,
      ggplot2::aes(x = time, y = posterior_mean),
      color = "#D73027",
      linewidth = 0.75
    ) +
    ggplot2::geom_errorbar(
      data = data$observed,
      ggplot2::aes(
        x = time,
        ymin = beta - 2 * standard_error,
        ymax = beta + 2 * standard_error
      ),
      width = 0.15,
      color = "grey25",
      linewidth = 0.34
    ) +
    ggplot2::geom_point(
      data = data$observed,
      ggplot2::aes(x = time, y = beta),
      size = 1.25,
      color = "black"
    )
  if (show_parametric) {
    plot <- plot +
      ggplot2::geom_line(
        data = data$parametric,
        ggplot2::aes(x = time, y = linear),
        color = "#1B9E77",
        linetype = "dashed",
        linewidth = 0.45
      ) +
      ggplot2::geom_line(
        data = data$parametric,
        ggplot2::aes(x = time, y = quadratic),
        color = "#7570B3",
        linetype = "dashed",
        linewidth = 0.45
      )
  }
  if (include_zero_line) {
    plot <- plot + ggplot2::geom_hline(
      yintercept = 0,
      color = "#377EB8",
      linetype = "dashed",
      linewidth = 0.4
    )
  }
  plot +
    ggplot2::facet_wrap(ggplot2::vars(facet), ncol = ncol, scales = "free_y") +
    ggplot2::labs(x = "Time", y = "eQTL effect") +
    ggplot2::theme_bw(base_size = 10.5) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold", size = 9),
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(size = 8.5)
    )
}

plot_figure3_same_units <- function() {
  plot_unit_comparison(
    figure3_status$key,
    c("Current FASH", "FASH-CL"),
    "IWP1",
    ncol = 3L,
    show_parametric = TRUE
  )
}

plot_figure3_proposed <- function() {
  plot_unit_comparison(
    figure3_proposed$key,
    "FASH-CL",
    "IWP1",
    ncol = 2L,
    show_parametric = TRUE
  )
}

plot_figure5_same_units <- function() {
  plot_unit_comparison(
    figure5_status$key,
    c("Current FASH", "FASH-CL"),
    "IWP2",
    ncol = 2L,
    show_parametric = TRUE
  )
}

plot_figure5_proposed <- function() {
  plot_unit_comparison(
    figure5_proposed$key,
    "FASH-CL",
    "IWP2",
    ncol = 2L,
    show_parametric = TRUE
  )
}

plot_figure6_same_units <- function() {
  plot_unit_comparison(
    figure6_status$key,
    c("Current FASH", "FASH-CL"),
    "IWP1",
    ncol = 2L,
    show_parametric = FALSE,
    include_zero_line = TRUE
  )
}

plot_figure6_proposed <- function() {
  plot_unit_comparison(
    figure6_proposed$key,
    "FASH-CL",
    "IWP1",
    ncol = 2L,
    show_parametric = FALSE,
    include_zero_line = TRUE
  )
}

plot_lfdr_scatter <- function() {
  status_colors <- c(
    "Neither" = "grey72",
    "Current FASH only" = "#0072B2",
    "FASH-CL only" = "#D55E00",
    "Both" = "#009E73"
  )
  status_alpha <- c(
    "Neither" = 0.055,
    "Current FASH only" = 0.55,
    "FASH-CL only" = 0.55,
    "Both" = 0.8
  )
  make_panel <- function(data, selection, point_size) {
    summary <- lfdr_comparison_summary[
      lfdr_comparison_summary$selection == selection,
      ,
      drop = FALSE
    ]
    ggplot2::ggplot(
      data,
      ggplot2::aes(
        x = current_lfdr,
        y = cl_lfdr,
        color = discovery_status,
        alpha = discovery_status
      )
    ) +
      ggrastr::geom_point_rast(
        size = point_size,
        stroke = 0,
        raster.dpi = 300
      ) +
      ggplot2::geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        color = "grey20",
        linewidth = 0.45
      ) +
      ggplot2::annotate(
        "label",
        x = 0.04,
        y = 0.96,
        hjust = 0,
        vjust = 1,
        size = 3.1,
        label = paste0(
          "n = ", format(summary$pair_count, big.mark = ","),
          "\nSpearman = ",
          format_decimal(summary$spearman_correlation, 3L)
        ),
        linewidth = 0.2,
        fill = scales::alpha("white", 0.82)
      ) +
      ggplot2::scale_color_manual(values = status_colors, drop = FALSE) +
      ggplot2::scale_alpha_manual(
        values = status_alpha,
        guide = "none",
        drop = FALSE
      ) +
      ggplot2::scale_x_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2)
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2)
      ) +
      ggplot2::coord_equal() +
      ggplot2::labs(
        title = selection,
        x = "Current-FASH BF-adjusted lfdr",
        y = "FASH-CL BF-adjusted lfdr",
        color = "FDR 0.05 call"
      ) +
      ggplot2::theme_bw(base_size = 10.5) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", size = 11),
        legend.position = "bottom"
      )
  }
  patchwork::wrap_plots(
    make_panel(
      lfdr_scatter_all,
      "All tested gene-variant pairs",
      point_size = 0.22
    ),
    make_panel(
      lfdr_scatter_top,
      "Current-FASH top pair per gene",
      point_size = 0.55
    ),
    ncol = 2L,
    guides = "collect"
  ) & ggplot2::theme(legend.position = "bottom")
}

plot_iwp1_prior_weights <- function() {
  data <- prior_weight_comparison
  data$method <- factor(data$method, levels = c("Current FASH", "FASH-CL"))
  data$adjustment <- factor(
    data$adjustment,
    levels = c("Raw", "BF-adjusted")
  )
  positive_psd <- sort(unique(data$psd[data$psd > 0]))
  positive_weight <- data$prior_weight[data$prior_weight > 0]
  x_breaks <- unique(c(
    0,
    scales::breaks_log(n = 6)(range(positive_psd))
  ))
  x_breaks <- x_breaks[x_breaks == 0 | (
    x_breaks >= min(positive_psd) & x_breaks <= max(positive_psd)
  )]
  y_breaks <- unique(c(
    0,
    scales::breaks_log(n = 6)(range(positive_weight))
  ))
  y_breaks <- y_breaks[y_breaks == 0 | (
    y_breaks >= min(positive_weight) & y_breaks <= max(positive_weight)
  )]
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = psd,
      y = prior_weight,
      color = method,
      group = method
    )
  ) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_point(size = 1.7) +
    ggplot2::facet_wrap(ggplot2::vars(adjustment), ncol = 2L) +
    ggplot2::scale_color_manual(values = c(
      "Current FASH" = "#0072B2",
      "FASH-CL" = "#D55E00"
    )) +
    ggplot2::scale_x_continuous(
      trans = scales::pseudo_log_trans(
        sigma = min(positive_psd) / 2,
        base = 10
      ),
      breaks = x_breaks,
      labels = scales::label_number(accuracy = 0.001)
    ) +
    ggplot2::scale_y_continuous(
      trans = scales::pseudo_log_trans(
        sigma = max(min(positive_weight) / 2, 1e-12),
        base = 10
      ),
      breaks = y_breaks,
      labels = scales::label_number(accuracy = 0.0001)
    ) +
    ggplot2::labs(
      x = "IWP1 PSD grid (0 is the null component)",
      y = "Prior weight",
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 10.5) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

plot_venn_diagram_comparison <- function(unit) {
  if (!unit %in% c("Genes", "Gene-variant pairs", "Variants")) {
    stop("Unknown Venn unit.")
  }
  plots <- lapply(names(venn_sets), function(method) {
    display_sets <- venn_sets[[method]][[unit]]
    names(display_sets) <- c(
      "Quadratic",
      "Linear",
      if (identical(method, "Current FASH")) "FASH" else "FASH-CL"
    )
    ggVennDiagram::ggVennDiagram(
      display_sets,
      label = "both",
      label_alpha = 0,
      set_color = "grey20",
      set_size = 3.8,
      label_size = 3.4,
      edge_size = 0.8
    ) +
      ggplot2::scale_fill_gradient(
        low = "grey95",
        high = "#D73027",
        limits = c(
          0,
          max(venn_regions$count[venn_regions$unit == unit])
        ),
        oob = scales::squish,
        labels = scales::label_comma()
      ) +
      ggplot2::labs(title = method) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold",
          hjust = 0.5,
          size = 13
        ),
        legend.position = "none",
        plot.margin = ggplot2::margin(12, 18, 22, 18)
      )
  })
  patchwork::wrap_plots(plots, ncol = 2L)
}

current_dynamic_pairs <- discovery_summary$pair_count[
  discovery_summary$method == "Current FASH" &
    discovery_summary$test == "IWP1 dynamic"
]
cl_dynamic_pairs <- discovery_summary$pair_count[
  discovery_summary$method == "FASH-CL" &
    discovery_summary$test == "IWP1 dynamic"
]
current_dynamic_genes <- discovery_summary$gene_count[
  discovery_summary$method == "Current FASH" &
    discovery_summary$test == "IWP1 dynamic"
]
cl_dynamic_genes <- discovery_summary$gene_count[
  discovery_summary$method == "FASH-CL" &
    discovery_summary$test == "IWP1 dynamic"
]
figure3_retained_count <- sum(figure3_status$cl_significant)
figure5_retained_count <- sum(figure5_status$cl_significant)
figure6_retained_count <- sum(
  !is.na(figure6_status$cl_category_significant) &
    figure6_status$cl_category_significant
)
current_fash_only_pairs <- venn_regions$count[
  venn_regions$method == "Current FASH" &
    venn_regions$unit == "Gene-variant pairs" &
    venn_regions$region_code == "001"
]
cl_fash_only_pairs <- venn_regions$count[
  venn_regions$method == "FASH-CL" &
    venn_regions$unit == "Gene-variant pairs" &
    venn_regions$region_code == "001"
]

lookup_hallmark <- function(method, gene_set, term, column) {
  value <- hallmark_comparison[
    hallmark_comparison$method == method &
      hallmark_comparison$gene_set == gene_set &
      hallmark_comparison$term == term,
    column
  ]
  if (length(value) != 1L) {
    stop("Hallmark lookup failed.")
  }
  value
}

current_all_hypoxia_q <- lookup_hallmark(
  "Current FASH", "All dynamic", "HALLMARK_HYPOXIA", "q_value"
)
cl_all_hypoxia_q <- lookup_hallmark(
  "FASH-CL", "All dynamic", "HALLMARK_HYPOXIA", "q_value"
)
current_switch_hypoxia_q <- lookup_hallmark(
  "Current FASH", "Switch", "HALLMARK_HYPOXIA", "q_value"
)
cl_switch_hypoxia_q <- lookup_hallmark(
  "FASH-CL", "Switch", "HALLMARK_HYPOXIA", "q_value"
)
current_switch_kras_q <- lookup_hallmark(
  "Current FASH", "Switch", "HALLMARK_KRAS_SIGNALING_UP", "q_value"
)
cl_switch_kras_q <- lookup_hallmark(
  "FASH-CL", "Switch", "HALLMARK_KRAS_SIGNALING_UP", "q_value"
)
