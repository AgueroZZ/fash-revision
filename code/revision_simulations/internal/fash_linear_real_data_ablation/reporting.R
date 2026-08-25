# Reporting helpers for the predictive-SD mixture FASH-linear real-data ablation.

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

invisible(lapply(
  c("ggplot2", "ggrastr", "ggVennDiagram", "knitr", "patchwork", "scales"),
  require_reporting_package
))

workflowr_root <- find_workflowr_root()
shared_path <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "shared",
  "simulation_functions.R"
)
if (!file.exists(shared_path)) {
  stop("The shared simulation functions are missing.")
}
source(shared_path)

output_id <- "fash_linear_real_data_ablation_mixture_predstep1_penalty10"
analysis_id <- paste0(
  "revision_internal_fash_linear_real_data_ablation_",
  "mixture_predstep1_penalty10"
)
gallery_analysis_id <- paste0(
  "fash_linear_discordant_example_gallery_",
  "mixture_predstep1_penalty10"
)
output_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  output_id
)
cache_path <- file.path(output_directory, "analysis_cache.rds")
run_status_path <- file.path(output_directory, "run_status.rds")
linear_raw_path <- file.path(output_directory, "linear_fit_raw.rds")
linear_bf_path <- file.path(output_directory, "linear_fit_bf.rds")
if (!file.exists(cache_path)) {
  stop(
    "The retained predictive-SD mixture FASH-linear cache is missing. ",
    "Run the page-specific analysis runner first."
  )
}
if (!file.exists(run_status_path) ||
    !file.exists(linear_raw_path) ||
    !file.exists(linear_bf_path)) {
  stop("The retained run status or compact FASH-linear fits are missing.")
}
cache <- readRDS(cache_path)
required_fields <- c(
  "configuration", "input_provenance", "fit_provenance",
  "unchanged_iwp_validation", "discovery_counts", "overlap_summary",
  "venn_sets", "lfdr_scatter_all", "lfdr_scatter_top", "lfdr_summary",
  "prior_summary", "prior_weights", "linear_version_comparison",
  "example_summary", "example_observed", "example_iwp",
  "example_linear", "runtime_summary"
)
if (!identical(names(cache), required_fields)) {
  stop("The retained FASH-linear cache failed structural validation.")
}

expected_grid <- default_revision_grid()
configuration <- cache$configuration
if (!identical(configuration$analysis_id, analysis_id) ||
    !identical(configuration$output_id, output_id) ||
    !identical(configuration$alpha, 0.05) ||
    !identical(configuration$pc_correction, "Time-specific PCs") ||
    !identical(configuration$linear_prior_mode, "mixture_grid") ||
    !identical(
      configuration$scale_definition,
      "sd_linear_departure_at_pred_step"
    ) ||
    !isTRUE(all.equal(configuration$psd_grid, expected_grid, tolerance = 0)) ||
    length(configuration$psd_grid) != 52L ||
    !identical(configuration$pred_step, 1) ||
    !identical(configuration$penalty, 10L) ||
    !identical(configuration$time_grid, 0:15) ||
    !identical(configuration$historical_output_id, "fash_linear_real_data_ablation")) {
  stop("The retained FASH-linear cache failed configuration validation.")
}

run_status <- readRDS(run_status_path)
if (!is.list(run_status) ||
    !identical(run_status$analysis_id, analysis_id) ||
    !identical(run_status$status, "complete") ||
    !identical(
      unname(tools::md5sum(cache_path)),
      run_status$analysis_cache_md5
    )) {
  stop("The retained FASH-linear run is not complete or its cache changed.")
}

input_provenance <- cache$input_provenance
current_input_md5 <- unname(tools::md5sum(input_provenance$path))
if (anyNA(current_input_md5) ||
    !identical(current_input_md5, input_provenance$md5)) {
  stop("At least one retained scientific input changed after cache creation.")
}

fit_provenance <- cache$fit_provenance
expected_fit_paths <- normalizePath(
  c(linear_raw_path, linear_bf_path),
  winslash = "/",
  mustWork = TRUE
)
if (!identical(fit_provenance$label, c("linear_raw", "linear_bf")) ||
    !identical(
      normalizePath(
        fit_provenance$path,
        winslash = "/",
        mustWork = TRUE
      ),
      expected_fit_paths
    )) {
  stop("The compact FASH-linear fit provenance is not the expected version.")
}
current_fit_md5 <- unname(tools::md5sum(expected_fit_paths))
if (anyNA(current_fit_md5) ||
    !identical(current_fit_md5, fit_provenance$md5) ||
    !identical(
      stats::setNames(current_fit_md5, fit_provenance$label),
      run_status$fit_md5
    )) {
  stop("At least one compact FASH-linear fit changed after cache creation.")
}

linear_fit_raw <- readRDS(linear_raw_path)
linear_fit_bf <- readRDS(linear_bf_path)
validate_compact_linear_mixture_fash(
  linear_fit_raw,
  expected_grid = expected_grid,
  expected_pred_step = 1,
  expected_penalty = 10L
)
validate_compact_linear_mixture_fash(
  linear_fit_bf,
  expected_grid = expected_grid,
  expected_pred_step = 1,
  expected_penalty = 10L
)
if (isTRUE(linear_fit_raw$bf_adjusted) ||
    !isTRUE(linear_fit_bf$bf_adjusted) ||
    !identical(linear_fit_raw$unit_ids, linear_fit_bf$unit_ids) ||
    length(linear_fit_raw$unit_ids) != 1009173L) {
  stop("The compact FASH-linear fits failed cross-fit validation.")
}
rm(linear_fit_raw, linear_fit_bf)

expected_validation_checks <- c(
  "input_md5_datasets",
  "input_md5_current_iwp1_raw",
  "input_md5_current_iwp1_bf",
  "input_md5_historical_profile_raw",
  "input_md5_historical_profile_bf",
  "input_md5_historical_profile_cache",
  "pair_count_and_order",
  "matched_grid_predstep_penalty",
  "raw_discovery_counts",
  "bf_discovery_counts",
  "historical_bf_lfdr"
)
unchanged_iwp_validation <- cache$unchanged_iwp_validation
if (!identical(
      unchanged_iwp_validation$check,
      expected_validation_checks
    ) ||
    nrow(unchanged_iwp_validation) != 11L ||
    !all(unchanged_iwp_validation$pass)) {
  stop("The unchanged-IWP validation contract did not pass all 11 checks.")
}

gallery_cache_path <- file.path(
  output_directory,
  "example_gallery.rds"
)
if (!file.exists(gallery_cache_path)) {
  stop(
    "The retained example-gallery cache is missing. Run ",
    "build_example_gallery.R before rendering this page."
  )
}
gallery_cache <- readRDS(gallery_cache_path)
gallery_required <- c(
  "configuration", "provenance", "gallery_summary", "gallery_observed",
  "gallery_iwp", "gallery_linear", "runtime_seconds"
)
if (!all(gallery_required %in% names(gallery_cache)) ||
    !identical(
      gallery_cache$configuration$analysis_id,
      gallery_analysis_id
    ) ||
    !identical(gallery_cache$configuration$pair_count, 6L) ||
    !identical(gallery_cache$configuration$linear_only_count, 3L) ||
    !identical(gallery_cache$configuration$iwp1_only_count, 3L) ||
    !identical(gallery_cache$configuration$pred_step, 1) ||
    !identical(gallery_cache$configuration$penalty, 10L) ||
    !identical(
      gallery_cache$configuration$linear_model,
      "Predictive-SD finite-mixture FASH-linear"
    )) {
  stop("The retained example-gallery cache failed structural validation.")
}
gallery_current_md5 <- unname(tools::md5sum(gallery_cache$provenance$path))
expected_gallery_provenance_labels <- c(
  "current_bf", "linear_bf", "linear_run_status", "explorer_index",
  "shared_functions", "helper", "explorer_helper", "builder"
)
if (!identical(
      gallery_cache$provenance$label,
      expected_gallery_provenance_labels
    ) ||
    anyNA(gallery_current_md5) ||
    !identical(gallery_current_md5, gallery_cache$provenance$md5)) {
  stop("At least one example-gallery input changed after cache creation.")
}
gallery_linear_fit_md5 <- gallery_cache$provenance$md5[
  gallery_cache$provenance$label == "linear_bf"
]
if (length(gallery_linear_fit_md5) != 1L ||
    !identical(gallery_linear_fit_md5, current_fit_md5[2])) {
  stop("The example gallery does not use the retained compact linear fit.")
}

discovery_counts <- cache$discovery_counts
overlap_summary <- cache$overlap_summary
venn_sets <- cache$venn_sets
lfdr_scatter_all <- cache$lfdr_scatter_all
lfdr_scatter_top <- cache$lfdr_scatter_top
lfdr_summary <- cache$lfdr_summary
prior_summary <- cache$prior_summary
prior_weights <- cache$prior_weights
linear_version_comparison <- cache$linear_version_comparison
example_summary <- cache$example_summary
example_observed <- cache$example_observed
example_iwp <- cache$example_iwp
example_linear <- cache$example_linear
gallery_summary <- gallery_cache$gallery_summary
gallery_observed <- gallery_cache$gallery_observed
gallery_iwp <- gallery_cache$gallery_iwp
gallery_linear <- gallery_cache$gallery_linear
runtime_summary <- cache$runtime_summary

expected_discovery_counts <- data.frame(
  method = rep(c("Current FASH", "FASH-linear"), each = 2L),
  adjustment = rep(c("Raw", "BF-adjusted"), times = 2L),
  pair_count = c(43860L, 9205L, 60188L, 14900L),
  gene_count = c(3258L, 1177L, 3863L, 1663L),
  variant_count = c(42893L, 9139L, 58570L, 14759L),
  stringsAsFactors = FALSE
)
if (!identical(discovery_counts, expected_discovery_counts) ||
    !identical(overlap_summary$current_count, c(9205L, 1177L, 9139L)) ||
    !identical(overlap_summary$linear_count, c(14900L, 1663L, 14759L)) ||
    !identical(
      overlap_summary$intersection_count,
      c(8528L, 1113L, 8466L)
    ) ||
    !identical(overlap_summary$current_only_count, c(677L, 64L, 673L)) ||
    !identical(overlap_summary$linear_only_count, c(6372L, 550L, 6293L)) ||
    nrow(lfdr_scatter_all) != 1009173L ||
    nrow(lfdr_scatter_top) != 6362L ||
    any(!is.finite(lfdr_scatter_all$current_lfdr)) ||
    any(!is.finite(lfdr_scatter_all$linear_lfdr)) ||
    nrow(prior_weights) != 4L * length(expected_grid) ||
    !isTRUE(all.equal(
      sort(unique(prior_weights$predstep_sd)),
      expected_grid,
      tolerance = 0
    )) ||
    any(abs(stats::aggregate(
      prior_weight ~ method + adjustment,
      data = prior_weights,
      FUN = sum
    )$prior_weight - 1) > 1e-6) ||
    any(prior_summary$null_weight < 0 | prior_summary$null_weight > 1) ||
    !identical(
      prior_summary$active_nonnull_components,
      rep(8L, 4L)
    ) ||
    nrow(linear_version_comparison) != 2L ||
    !identical(
      as.character(linear_version_comparison$adjustment),
      c("Raw", "BF-adjusted")
    ) ||
    !identical(
      linear_version_comparison$profile_pair_count,
      c(45060L, 15865L)
    ) ||
    !identical(
      linear_version_comparison$mixture_pair_count,
      c(60188L, 14900L)
    ) ||
    !identical(
      linear_version_comparison$intersection_pair_count,
      c(45060L, 12655L)
    ) ||
    !identical(
      linear_version_comparison$profile_only_pair_count,
      c(0L, 3210L)
    ) ||
    !identical(
      linear_version_comparison$mixture_only_pair_count,
      c(15128L, 2245L)
    ) ||
    !isTRUE(all.equal(
      linear_version_comparison$pair_jaccard,
      c(0.7486542, 0.6987852),
      tolerance = 5e-7
    )) ||
    !isTRUE(all.equal(
      linear_version_comparison$spearman_lfdr,
      c(0.9832737, 0.9832618),
      tolerance = 5e-7
    )) ||
    !isTRUE(all.equal(
      lfdr_summary$spearman_correlation,
      c(0.9073864, 0.9051038),
      tolerance = 5e-7
    )) ||
    !isTRUE(all.equal(
      prior_summary$null_weight,
      c(0.4289903, 0.9381533, 0.4572348, 0.9192507),
      tolerance = 5e-7
    )) ||
    !isTRUE(all.equal(
      prior_summary$alternative_rms_predstep_sd,
      c(0.1601786, 0.1600649, 0.0955470, 0.0955321),
      tolerance = 5e-7
    )) ||
    nrow(example_summary) != 1L ||
    !identical(example_summary$gene_symbol, "MMP15") ||
    !identical(example_summary$variant_id, "rs12926803") ||
    nrow(example_observed) != 16L ||
    nrow(example_iwp) != 151L || nrow(example_linear) != 151L ||
    nrow(gallery_summary) != 6L ||
    nrow(gallery_observed) != 96L ||
    nrow(gallery_iwp) != 906L || nrow(gallery_linear) != 906L ||
    length(unique(gallery_summary$gene_id)) != 6L ||
    !identical(
      as.character(gallery_summary$direction),
      rep(c("FASH-linear only", "IWP1 only"), each = 3L)
    ) ||
    any(!is.finite(as.matrix(gallery_summary[, c(
      "iwp1_lfdr_bf", "linear_lfdr_bf", "strober_linear_pvalue",
      "strober_nonlinear_pvalue"
    )]))) ||
    any(!is.finite(as.matrix(example_observed[, c(
      "time", "beta", "standard_error"
    )]))) ||
    any(!is.finite(as.matrix(example_iwp[, c(
      "time", "posterior_mean", "lower", "upper"
    )]))) ||
    any(!is.finite(as.matrix(example_linear[, c(
      "time", "posterior_mean", "lower", "upper"
    )]))) ||
    !identical(
      names(venn_sets),
      c("Gene-variant pairs", "Genes", "Variants")
    )) {
  stop("The retained FASH-linear cache failed numerical validation.")
}

format_integer <- function(value) {
  format(as.integer(round(value)), big.mark = ",", scientific = FALSE)
}

format_decimal <- function(value, digits = 3L) {
  formatC(as.numeric(value), format = "f", digits = digits)
}

format_percent <- function(value, digits = 1L) {
  paste0(format_decimal(100 * value, digits), "%")
}

render_scrollable_table <- function(data,
                                    digits = 3L,
                                    minimum_width = "760px") {
  table_html <- knitr::kable(
    data,
    format = "html",
    escape = TRUE,
    row.names = FALSE,
    digits = digits,
    table.attr = paste0(
      'class="table table-striped table-hover" style="min-width:',
      minimum_width,
      ';"'
    )
  )
  cat(
    '<div class="ablation-table-scroll">\n',
    as.character(table_html),
    "\n</div>\n",
    sep = ""
  )
  invisible(table_html)
}

discovery_status_colors <- c(
  "Neither" = "#bdbdbd",
  "Current FASH only" = "#3182bd",
  "FASH-linear only" = "#e6550d",
  "Both" = "#31a354"
)

make_lfdr_panel <- function(data, title, point_size, point_alpha, summary_row) {
  annotation <- paste0(
    "n = ", format_integer(nrow(data)),
    "\nSpearman = ", format_decimal(summary_row$spearman_correlation, 3L)
  )
  ggplot2::ggplot(data, ggplot2::aes(
    x = current_lfdr,
    y = linear_lfdr,
    color = discovery_status
  )) +
    ggrastr::geom_point_rast(
      size = point_size,
      alpha = point_alpha,
      raster.dpi = 300
    ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    ggplot2::annotate(
      "label",
      x = 0.05,
      y = 0.94,
      label = annotation,
      hjust = 0,
      vjust = 1,
      size = 3.7,
      linewidth = 0.2
    ) +
    ggplot2::scale_color_manual(
      values = discovery_status_colors,
      drop = FALSE,
      name = "FDR 0.05 call"
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(
      override.aes = list(size = 2, alpha = 1)
    )) +
    ggplot2::coord_cartesian(xlim = c(-0.05, 1.05), ylim = c(-0.05, 1.05)) +
    ggplot2::labs(
      title = title,
      x = "IWP1 FASH BF-adjusted lfdr",
      y = "FASH-linear BF-adjusted lfdr"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold")
    )
}

plot_lfdr_scatter <- function() {
  all_summary <- lfdr_summary[
    lfdr_summary$selection == "All tested gene-variant pairs",
    ,
    drop = FALSE
  ]
  top_summary <- lfdr_summary[
    lfdr_summary$selection == "Current-FASH top pair per gene",
    ,
    drop = FALSE
  ]
  left <- make_lfdr_panel(
    lfdr_scatter_all,
    "All tested gene-variant pairs",
    point_size = 0.35,
    point_alpha = 0.16,
    summary_row = all_summary
  )
  right <- make_lfdr_panel(
    lfdr_scatter_top,
    "Current-FASH top pair per gene",
    point_size = 0.8,
    point_alpha = 0.6,
    summary_row = top_summary
  )
  left + right + patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}

plot_overlap_venns <- function() {
  plots <- lapply(names(venn_sets), function(unit) {
    display_sets <- venn_sets[[unit]]
    names(display_sets) <- c("IWP1", "Linear")
    ggVennDiagram::ggVennDiagram(
      display_sets,
      label_alpha = 0,
      label = "count",
      set_size = 0,
      label_size = 4
    ) +
      ggplot2::annotate(
        "text",
        x = -3,
        y = 4,
        label = "Linear",
        fontface = "bold",
        size = 4
      ) +
      ggplot2::annotate(
        "text",
        x = -3,
        y = 0,
        label = "IWP1",
        fontface = "bold",
        size = 4
      ) +
      ggplot2::scale_fill_gradient(low = "#f7fbff", high = "#2171b5") +
      ggplot2::labs(title = unit) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(
          face = "bold",
          hjust = 0.5,
          size = 12
        )
      )
  })
  patchwork::wrap_plots(plots, nrow = 1L)
}

plot_prior_weights <- function() {
  adjustment_levels <- c("Raw", "BF-adjusted")
  method_levels <- c("Current FASH", "FASH-linear")
  method_colors <- c(
    "Current FASH" = "#d7301f",
    "FASH-linear" = "#2c7fb8"
  )

  null_weights <- prior_weights[prior_weights$is_null, , drop = FALSE]
  null_weights$adjustment <- factor(
    null_weights$adjustment,
    levels = adjustment_levels
  )
  null_weights$method <- factor(
    null_weights$method,
    levels = method_levels
  )
  null_weights$pi0_label <- sprintf(
    "hat(pi)[0] == %.3f",
    null_weights$prior_weight
  )

  alternative <- prior_weights[!prior_weights$is_null, , drop = FALSE]
  alternative$adjustment <- factor(
    alternative$adjustment,
    levels = adjustment_levels
  )
  alternative$method <- factor(
    alternative$method,
    levels = method_levels
  )
  alternative$conditional_weight <- ave(
    alternative$prior_weight,
    alternative$method,
    alternative$adjustment,
    FUN = function(weight) weight / sum(weight)
  )

  null_panel <- ggplot2::ggplot(
    null_weights,
    ggplot2::aes(x = method, y = prior_weight, fill = method)
  ) +
    ggplot2::geom_col(width = 0.58) +
    ggplot2::geom_text(
      ggplot2::aes(
        y = pmin(prior_weight + 0.055, 1.02),
        label = pi0_label
      ),
      parse = TRUE,
      size = 3.6
    ) +
    ggplot2::facet_wrap(~adjustment, nrow = 1L) +
    ggplot2::scale_fill_manual(values = method_colors, guide = "none") +
    ggplot2::scale_y_continuous(
      limits = c(0, 1.05),
      breaks = c(0, 0.5, 1),
      labels = scales::label_number(accuracy = 0.1),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      x = NULL,
      y = expression("Exact-null weight " * hat(pi)[0])
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )

  alternative_panel <- ggplot2::ggplot(
    alternative,
    ggplot2::aes(
      x = predstep_sd,
      y = conditional_weight,
      color = method,
      group = method
    )
  ) +
    ggplot2::geom_line(linewidth = 0.85) +
    ggplot2::geom_point(
      data = alternative[alternative$active, , drop = FALSE],
      size = 2
    ) +
    ggplot2::facet_wrap(~adjustment, nrow = 1L) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_continuous(
      labels = scales::label_number(accuracy = 0.01)
    ) +
    ggplot2::scale_color_manual(values = method_colors, name = NULL) +
    ggplot2::labs(
      x = "One-step predictive SD",
      y = "Conditional alternative weight"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold")
    )

  null_panel / alternative_panel +
    patchwork::plot_layout(heights = c(0.9, 2.1), guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}

plot_linear_only_example <- function() {
  curve_colors <- c(
    "IWP1 FASH" = "#d7301f",
    "FASH-linear" = "#2c7fb8"
  )
  ggplot2::ggplot() +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "grey65",
      linewidth = 0.45
    ) +
    ggplot2::geom_ribbon(
      data = example_iwp,
      ggplot2::aes(x = time, ymin = lower, ymax = upper, fill = method),
      alpha = 0.15
    ) +
    ggplot2::geom_ribbon(
      data = example_linear,
      ggplot2::aes(x = time, ymin = lower, ymax = upper, fill = method),
      alpha = 0.15
    ) +
    ggplot2::geom_line(
      data = example_iwp,
      ggplot2::aes(x = time, y = posterior_mean, color = method),
      linewidth = 1
    ) +
    ggplot2::geom_line(
      data = example_linear,
      ggplot2::aes(x = time, y = posterior_mean, color = method),
      linewidth = 1
    ) +
    ggplot2::geom_errorbar(
      data = example_observed,
      ggplot2::aes(
        x = time,
        ymin = beta - 2 * standard_error,
        ymax = beta + 2 * standard_error
      ),
      width = 0.12,
      color = "grey30",
      linewidth = 0.45
    ) +
    ggplot2::geom_point(
      data = example_observed,
      ggplot2::aes(x = time, y = beta),
      color = "black",
      size = 1.8
    ) +
    ggplot2::scale_color_manual(
      values = curve_colors,
      breaks = names(curve_colors),
      name = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = curve_colors,
      breaks = names(curve_colors),
      name = NULL
    ) +
    ggplot2::labs(
      title = paste0(
        example_summary$gene_symbol,
        " / ",
        example_summary$variant_id
      ),
      subtitle = paste0(
        "IWP1 lfdr = ", format_decimal(example_summary$current_lfdr, 4L),
        "; FASH-linear lfdr = ",
        format_decimal(example_summary$linear_lfdr, 4L)
      ),
      x = "Differentiation time",
      y = "Estimated eQTL effect"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(face = "bold")
    )
}

plot_example_gallery <- function() {
  curve_colors <- c(
    "IWP1 FASH" = "#d7301f",
    "FASH-linear" = "#2c7fb8"
  )
  ggplot2::ggplot() +
    ggplot2::geom_hline(
      yintercept = 0,
      color = "grey70",
      linewidth = 0.4
    ) +
    ggplot2::geom_ribbon(
      data = gallery_iwp,
      ggplot2::aes(x = time, ymin = lower, ymax = upper, fill = method),
      alpha = 0.13
    ) +
    ggplot2::geom_ribbon(
      data = gallery_linear,
      ggplot2::aes(x = time, ymin = lower, ymax = upper, fill = method),
      alpha = 0.13
    ) +
    ggplot2::geom_line(
      data = gallery_iwp,
      ggplot2::aes(x = time, y = posterior_mean, color = method),
      linewidth = 0.75
    ) +
    ggplot2::geom_line(
      data = gallery_linear,
      ggplot2::aes(x = time, y = posterior_mean, color = method),
      linewidth = 0.75
    ) +
    ggplot2::geom_errorbar(
      data = gallery_observed,
      ggplot2::aes(
        x = time,
        ymin = beta - 2 * standard_error,
        ymax = beta + 2 * standard_error
      ),
      width = 0.1,
      color = "grey35",
      linewidth = 0.32
    ) +
    ggplot2::geom_point(
      data = gallery_observed,
      ggplot2::aes(x = time, y = beta),
      color = "black",
      size = 1.1
    ) +
    ggplot2::facet_wrap(~facet_label, ncol = 3L, scales = "free_y") +
    ggplot2::scale_color_manual(
      values = curve_colors,
      breaks = names(curve_colors),
      name = NULL
    ) +
    ggplot2::scale_fill_manual(
      values = curve_colors,
      breaks = names(curve_colors),
      name = NULL
    ) +
    ggplot2::labs(
      x = "Differentiation time",
      y = "Estimated eQTL effect"
    ) +
    ggplot2::theme_bw(base_size = 10.5) +
    ggplot2::theme(
      legend.position = "bottom",
      strip.text = ggplot2::element_text(
        face = "bold",
        size = 8.4,
        lineheight = 1.05
      ),
      panel.grid.minor = ggplot2::element_blank()
    )
}

current_raw_count <- discovery_counts[
  discovery_counts$method == "Current FASH" &
    discovery_counts$adjustment == "Raw",
  ,
  drop = FALSE
]
linear_raw_count <- discovery_counts[
  discovery_counts$method == "FASH-linear" &
    discovery_counts$adjustment == "Raw",
  ,
  drop = FALSE
]
current_bf_count <- discovery_counts[
  discovery_counts$method == "Current FASH" &
    discovery_counts$adjustment == "BF-adjusted",
  ,
  drop = FALSE
]
linear_bf_count <- discovery_counts[
  discovery_counts$method == "FASH-linear" &
    discovery_counts$adjustment == "BF-adjusted",
  ,
  drop = FALSE
]
pair_overlap <- overlap_summary[
  overlap_summary$unit == "Gene-variant pairs",
  ,
  drop = FALSE
]
variant_overlap <- overlap_summary[
  overlap_summary$unit == "Variants",
  ,
  drop = FALSE
]
raw_version_comparison <- linear_version_comparison[
  linear_version_comparison$adjustment == "Raw",
  ,
  drop = FALSE
]
bf_version_comparison <- linear_version_comparison[
  linear_version_comparison$adjustment == "BF-adjusted",
  ,
  drop = FALSE
]
linear_raw_prior <- prior_summary[
  prior_summary$method == "FASH-linear" &
    prior_summary$adjustment == "Raw",
  ,
  drop = FALSE
]
linear_bf_prior <- prior_summary[
  prior_summary$method == "FASH-linear" &
    prior_summary$adjustment == "BF-adjusted",
  ,
  drop = FALSE
]
active_prior_weights <- prior_weights[
  prior_weights$active,
  c("method", "adjustment", "predstep_sd", "prior_weight", "is_null"),
  drop = FALSE
]
active_prior_weights$component <- ifelse(
  active_prior_weights$is_null,
  "Exact null",
  "Positive SD"
)
active_prior_weights$is_null <- NULL
active_prior_weights <- active_prior_weights[
  c("method", "adjustment", "component", "predstep_sd", "prior_weight")
]
