# Reporting helpers for the R6 enhancer enrichment page.
#
# The page reports enhancer annotations only, for current FASH against the two
# Strober analyses, in two sections: an unfiltered method comparison and a
# variant-level FASH-only view. It deliberately does not surface BH q-values;
# the analysis is exploratory and a within-set correction over four or five
# annotations was a narrower family than the page actually tests.

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
require_reporting_package("knitr")

workflowr_root <- find_workflowr_root()
cache_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fash_strober_enhancer_comparison",
  "analysis_cache.rds"
)
if (!file.exists(cache_path)) {
  stop(
    "The retained comparison cache is missing. Run ",
    "code/revision_simulations/internal/fash_strober_enhancer_comparison/",
    "run_fash_strober_enhancer_comparison.R first."
  )
}
cache <- readRDS(cache_path)
required_fields <- c(
  "configuration", "input_provenance", "set_metadata",
  "discovery_set_summary", "discovery_overlap", "strict_exclusivity_audit",
  "discovery_export", "coverage_by_set", "enrichment_results",
  "unmatched_enrichment", "enhancer_results", "enhancer_summary",
  "matching_balance", "matching_relaxation", "maximum_absolute_smd",
  "runtime_summary", "session_info"
)
expected_sets <- c(
  "current_all", "current_lead", "linear_all", "linear_lead",
  "quadratic_all", "quadratic_lead", "current_only_all", "current_only_lead"
)
if (!all(required_fields %in% names(cache)) ||
    !identical(
      cache$configuration$analysis_id,
      "revision_internal_fash_strober_enhancer_comparison"
    ) ||
    cache$configuration$matching_seed_count != 100L ||
    cache$configuration$controls_per_variant != 5L ||
    !identical(cache$set_metadata$discovery_set, expected_sets) ||
    !identical(cache$configuration$jackknife_blocks, as.character(1:22))) {
  stop("The retained comparison cache failed structural validation.")
}

current_md5 <- unname(tools::md5sum(cache$input_provenance$path))
if (anyNA(current_md5) ||
    !identical(current_md5, cache$input_provenance$md5)) {
  stop("At least one retained analysis input changed after cache creation.")
}

configuration <- cache$configuration
input_provenance <- cache$input_provenance
set_metadata <- cache$set_metadata
discovery_set_summary <- cache$discovery_set_summary
strict_exclusivity_audit <- cache$strict_exclusivity_audit
discovery_export <- cache$discovery_export
coverage_by_set <- cache$coverage_by_set
enrichment_results <- cache$enrichment_results
enhancer_results <- cache$enhancer_results
enhancer_summary <- cache$enhancer_summary
matching_balance <- cache$matching_balance
runtime_summary <- cache$runtime_summary

n_sets <- length(expected_sets)
if (nrow(enrichment_results) != n_sets * (29L + 83L) ||
    nrow(enhancer_results) != n_sets * (4L + 5L) ||
    nrow(matching_balance) != 2L * n_sets * 100L * 4L ||
    strict_exclusivity_audit$overlap_with_strober_union[1L] != 0L ||
    !setequal(unique(enrichment_results$annotation_system),
              configuration$annotation_systems) ||
    any(!is.finite(matching_balance$standardized_mean_difference))) {
  stop("The retained comparison cache failed result validation.")
}

n_enhancer_tests <- sum(is.finite(enhancer_results$p_value))

format_integer <- function(value) {
  format(as.integer(round(value)), big.mark = ",", scientific = FALSE)
}

format_decimal <- function(value, digits = 2L) {
  formatC(as.numeric(value), format = "f", digits = digits)
}

format_percent <- function(value, digits = 1L) {
  paste0(format_decimal(100 * value, digits), "%")
}

format_fold <- function(log2_enrichment, digits = 2L) {
  paste0(format_decimal(2^log2_enrichment, digits), "x")
}

lookup_runtime <- function(stage) {
  value <- runtime_summary$elapsed_seconds[runtime_summary$stage == stage]
  if (length(value) != 1L) {
    stop("Runtime lookup failed for: ", stage)
  }
  value
}

render_scrollable_table <- function(data,
                                    digits = 3L,
                                    align = NULL,
                                    minimum_width = "760px",
                                    escape = TRUE) {
  table_html <- knitr::kable(
    data,
    format = "html",
    escape = escape,
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
    '<div class="fash-strober-table-scroll">\n',
    as.character(table_html),
    "\n</div>\n",
    sep = ""
  )
  invisible(table_html)
}

system_levels <- c("Custom regulatory", "baselineLD v2.2")
strategy_levels <- c("All variants", "One lead variant per gene")
method_colors <- c(
  "Current FASH" = "#0072B2",
  "Strober linear" = "#D55E00",
  "Strober quadratic/nonlinear" = "#8E5BB7"
)
reference_color <- "#9a9a9a"
enhancer_labels <- c(
  "ENCODE cCRE enhancer-like" = "ENCODE enhancer-like cCRE",
  "Roadmap E020 iPS-20b: Enhancer" = "E020 enhancer (iPSC)",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer" =
    "E013 enhancer (mesoderm)",
  "Roadmap E095 left ventricle: Enhancer" = "E095 enhancer (left ventricle)",
  "Enhancer_Andersson" = "Enhancer (Andersson)",
  "Enhancer_Hoffman" = "Enhancer (Hoffman)",
  "WeakEnhancer_Hoffman" = "Weak enhancer (Hoffman)",
  "SuperEnhancer_Hnisz" = "Super-enhancer (Hnisz)",
  "Human_Enhancer_Villar" = "Human enhancer (Villar)"
)

prepare_enhancer_plot_data <- function(data) {
  data$annotation_label <- unname(enhancer_labels[data$annotation])
  data$annotation_label <- factor(
    data$annotation_label,
    levels = rev(unname(enhancer_labels))
  )
  data$annotation_system <- factor(
    data$annotation_system,
    levels = system_levels
  )
  data$selection_strategy <- factor(
    data$selection_strategy,
    levels = strategy_levels
  )
  data[
    is.finite(data$log2_enrichment) & is.finite(data$ci_lower) &
      is.finite(data$ci_upper),
    ,
    drop = FALSE
  ]
}

# Sparse categories such as baselineLD Andersson can produce intervals several
# times wider than everything else, which would compress every other row into
# illegibility. Clip the axis to a range that still contains every point
# estimate, and flag clipped interval ends with an arrow so nothing is hidden.
clip_intervals <- function(data, margin = 0.08, keep = 0.90) {
  # Clip only genuine outliers: keep the central `keep` share of interval bounds
  # and never clip so hard that a point estimate would fall off the panel.
  limits <- c(
    min(
      stats::quantile(data$ci_lower, 1 - keep, na.rm = TRUE, names = FALSE),
      min(data$log2_enrichment, na.rm = TRUE) - margin
    ),
    max(
      stats::quantile(data$ci_upper, keep, na.rm = TRUE, names = FALSE),
      max(data$log2_enrichment, na.rm = TRUE) + margin
    )
  )
  data$clipped_low <- data$ci_lower < limits[1L]
  data$clipped_high <- data$ci_upper > limits[2L]
  data$ci_lower_clipped <- pmax(data$ci_lower, limits[1L])
  data$ci_upper_clipped <- pmin(data$ci_upper, limits[2L])
  attr(data, "x_limits") <- limits
  data
}

# One shared figure grammar for all three sections: annotation system down the
# rows, all-variants versus lead-per-gene across the columns, series in colour.
enhancer_forest <- function(data, series_colors) {
  data$series <- factor(data$series, levels = names(series_colors))
  data <- clip_intervals(data)
  x_limits <- attr(data, "x_limits")
  dodge <- ggplot2::position_dodge(width = 0.72)
  # `group = series` is load-bearing. Without it, geom_point below dodges by
  # interaction(series, shape) while geom_errorbar dodges by series alone, so
  # the markers land on a different row from their own interval.
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = log2_enrichment,
      y = annotation_label,
      color = series,
      group = series
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.45
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lower_clipped, xmax = ci_upper_clipped),
      orientation = "y",
      width = 0.30,
      linewidth = 0.48,
      position = dodge
    ) +
    # Chevrons rather than segment arrowheads: geom_point dodges identically to
    # the estimate points, so the marker always lands on its own series row.
    ggplot2::geom_point(
      data = function(values) values[values$clipped_high, , drop = FALSE],
      ggplot2::aes(x = ci_upper_clipped),
      shape = 62L,
      size = 2.6,
      position = dodge,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = function(values) values[values$clipped_low, , drop = FALSE],
      ggplot2::aes(x = ci_lower_clipped),
      shape = 60L,
      size = 2.6,
      position = dodge,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = p_value < 0.05, fill = series),
      size = 2.5,
      stroke = 0.58,
      position = dodge
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(annotation_system),
      cols = ggplot2::vars(selection_strategy),
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_color_manual(values = series_colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = series_colors, drop = FALSE) +
    ggplot2::scale_shape_manual(
      values = c(`TRUE` = 21, `FALSE` = 1),
      labels = c(
        `TRUE` = "nominal p < 0.05 (uncorrected)",
        `FALSE` = "nominal p >= 0.05"
      )
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1, ncol = 1),
      fill = ggplot2::guide_legend(order = 1, ncol = 1),
      shape = ggplot2::guide_legend(
        order = 2,
        override.aes = list(fill = "grey35", color = "grey35")
      )
    ) +
    ggplot2::coord_cartesian(xlim = x_limits) +
    ggplot2::labs(
      x = "log2 enrichment relative to matched tested variants",
      y = NULL,
      color = NULL,
      fill = NULL,
      shape = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11.3) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 9.4)
    )
}

# Section 1: focal method against both Strober analyses, no exclusivity filter.
plot_method_comparison <- function(focal_method) {
  keep <- c(focal_method, "Strober linear", "Strober quadratic/nonlinear")
  data <- enhancer_results[
    enhancer_results$section == 1L & enhancer_results$method %in% keep,
    ,
    drop = FALSE
  ]
  data$series <- data$method
  colors <- method_colors[keep]
  enhancer_forest(prepare_enhancer_plot_data(data), colors)
}

# Section 2: the FASH-only set, with the unfiltered FASH set drawn behind it in
# grey so the reader can see what the exclusion actually changed.
plot_fash_only <- function() {
  exclusive_label <- "FASH only (Strober rsIDs removed)"
  reference_label <- "FASH, all discoveries (reference)"
  exclusive <- enhancer_results[
    enhancer_results$exclusivity == "Variant-level",
    ,
    drop = FALSE
  ]
  reference <- enhancer_results[
    enhancer_results$section == 1L &
      enhancer_results$method == "Current FASH",
    ,
    drop = FALSE
  ]
  if (!nrow(exclusive) || !nrow(reference)) {
    stop("FASH-only plot data is missing.")
  }
  exclusive$series <- exclusive_label
  reference$series <- reference_label
  colors <- stats::setNames(
    c(unname(method_colors[["Current FASH"]]), reference_color),
    c(exclusive_label, reference_label)
  )
  enhancer_forest(
    prepare_enhancer_plot_data(rbind(exclusive, reference)),
    colors
  )
}

# ---- Tables -----------------------------------------------------------------

# Rendered with escape = FALSE so the download column is clickable. Every value
# here is a controlled literal from the cache configuration, not user input.
link <- function(url, label) {
  ifelse(
    is.na(url),
    "",
    paste0(
      '<a href="', url, '" target="_blank" rel="noopener">', label, "</a>"
    )
  )
}
provenance_table <- configuration$annotation_provenance
provenance_table$download <- link(
  provenance_table$download_url,
  provenance_table$download_file
)
provenance_table$mirror <- link(
  provenance_table$retrieval_url,
  ifelse(is.na(provenance_table$retrieval_url), "", "Song Lab mirror")
)
provenance_table <- provenance_table[, c(
  "annotation_system", "annotation", "source", "accession_or_record",
  "genome_build", "repository", "download", "mirror"
)]
names(provenance_table) <- c(
  "Annotation system", "Annotation", "Source", "Accession or record",
  "Genome build", "Repository", "Download", "Retrieved from"
)

discovery_table <- discovery_set_summary
discovery_table <- merge(
  discovery_table,
  set_metadata[, c("discovery_set", "section", "exclusivity")],
  by = "discovery_set",
  all.x = TRUE,
  sort = FALSE
)
discovery_table <- discovery_table[
  order(discovery_table$section, match(
    discovery_table$discovery_set,
    expected_sets
  )),
  c(
    "section", "discovery_set", "method", "selection_strategy",
    "pair_count", "unique_variant_count", "unique_gene_count"
  )
]
names(discovery_table) <- c(
  "Section", "Set ID", "Method", "Selection strategy",
  "Gene-variant pairs", "Unique variants", "Unique genes"
)

coverage_table <- merge(
  coverage_by_set,
  set_metadata[, c("discovery_set", "method", "selection_strategy")],
  by = "discovery_set",
  all.x = TRUE,
  sort = FALSE
)
coverage_table$coverage_proportion <- format_percent(
  coverage_table$coverage_proportion,
  1L
)
coverage_table <- coverage_table[, c(
  "annotation_system", "discovery_set", "method", "selection_strategy",
  "original_variant_count", "covered_variant_count", "coverage_proportion"
)]
names(coverage_table) <- c(
  "Annotation system", "Set ID", "Method", "Selection strategy",
  "Original variants", "Covered variants", "Coverage"
)

# ---- Inline lookups ---------------------------------------------------------

summary_lookup <- function(system, set_name) {
  row <- enhancer_summary[
    enhancer_summary$annotation_system == system &
      enhancer_summary$discovery_set == set_name,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L) {
    stop("Enhancer summary lookup failed for ", system, " / ", set_name, ".")
  }
  row
}

median_fold_text <- function(system, set_name, digits = 2L) {
  row <- summary_lookup(system, set_name)
  paste0(
    format_fold(row$median_log2_enrichment, digits),
    " (median of ",
    row$finite_annotation_count,
    " enhancer estimates)"
  )
}

enhancer_lookup <- function(system, set_name, annotation) {
  row <- enhancer_results[
    enhancer_results$annotation_system == system &
      enhancer_results$discovery_set == set_name &
      enhancer_results$annotation == annotation,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L) {
    stop("Enhancer lookup failed for ", set_name, " / ", annotation, ".")
  }
  row
}

describe_estimate <- function(system, set_name, annotation) {
  row <- enhancer_lookup(system, set_name, annotation)
  if (!is.finite(row$enrichment) || !is.finite(row$ci_lower)) {
    return(paste0(
      "not estimable (only ",
      format_integer(row$selected_overlap),
      " of ",
      format_integer(row$selected_total),
      " selected variants overlap, below the minimum of 10)"
    ))
  }
  paste0(
    format_decimal(row$enrichment, 2L),
    "-fold (95% CI ",
    format_decimal(2^row$ci_lower, 2L),
    "-",
    format_decimal(2^row$ci_upper, 2L),
    ", nominal p = ",
    format(signif(row$p_value, 2L), scientific = FALSE),
    ", ",
    format_integer(row$selected_overlap),
    "/",
    format_integer(row$selected_total),
    " variants)"
  )
}

# Name whichever enhancers clear a nominal threshold, rather than hard-coding
# which one won; if the numbers move, the sentence moves with them.
describe_significant <- function(set_name, threshold = 0.05) {
  rows <- enhancer_results[
    enhancer_results$discovery_set == set_name &
      is.finite(enhancer_results$p_value) &
      enhancer_results$p_value < threshold,
    ,
    drop = FALSE
  ]
  if (!nrow(rows)) {
    return("no enhancer annotation reaches nominal p < 0.05")
  }
  rows <- rows[order(rows$p_value), , drop = FALSE]
  paste(
    vapply(seq_len(nrow(rows)), function(i) {
      paste0(
        unname(enhancer_labels[rows$annotation[i]]), " (",
        rows$annotation_system[i], ") at ",
        format_decimal(rows$enrichment[i], 2L), "-fold, 95% CI ",
        format_decimal(2^rows$ci_lower[i], 2L), "-",
        format_decimal(2^rows$ci_upper[i], 2L), ", p = ",
        format(signif(rows$p_value[i], 2L), scientific = FALSE), ", ",
        format_integer(rows$selected_overlap[i]), "/",
        format_integer(rows$selected_total[i]), " variants"
      )
    }, character(1)),
    collapse = "; "
  )
}

lead_gene_count <- discovery_set_summary$unique_gene_count[
  discovery_set_summary$discovery_set == "current_lead"
]

e013 <- "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer"
set_size <- function(set_name) {
  format_integer(
    set_metadata$original_variant_count[
      set_metadata$discovery_set == set_name
    ]
  )
}
