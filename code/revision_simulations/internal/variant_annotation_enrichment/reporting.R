# Reporting helpers for the internal FASH variant annotation enrichment page.

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
  "variant_annotation_enrichment",
  "analysis_cache.rds"
)
if (!file.exists(cache_path)) {
  stop(
    "The retained analysis cache is missing. Run ",
    "code/revision_simulations/internal/variant_annotation_enrichment/",
    "run_variant_annotation_enrichment.R first."
  )
}
cache <- readRDS(cache_path)
required_cache_fields <- c(
  "configuration", "input_provenance", "annotation_download_manifest",
  "discovery_summary", "set_metadata", "annotation_coverage",
  "enrichment_results", "matching_balance", "matching_seed_results",
  "matching_relaxation_audit", "ld_classification", "session_info"
)
if (!all(required_cache_fields %in% names(cache)) ||
    !identical(
      cache$configuration$analysis_id,
      "revision_internal_variant_annotation_enrichment"
    ) || cache$configuration$matching_seed_count != 100L ||
    cache$configuration$controls_per_variant != 5L) {
  stop("The retained analysis cache failed reporting validation.")
}

configuration <- cache$configuration
discovery_summary <- cache$discovery_summary
set_metadata <- cache$set_metadata
annotation_coverage <- cache$annotation_coverage
enrichment_results <- cache$enrichment_results
matching_balance <- cache$matching_balance
matching_seed_results <- cache$matching_seed_results
matching_relaxation_audit <- cache$matching_relaxation_audit
annotation_download_manifest <- cache$annotation_download_manifest
input_provenance <- cache$input_provenance

format_integer <- function(value) {
  format(as.integer(round(value)), big.mark = ",", scientific = FALSE)
}

format_decimal <- function(value, digits = 2L) {
  formatC(as.numeric(value), format = "f", digits = digits)
}

format_fold <- function(log2_enrichment, digits = 2L) {
  paste0(format_decimal(2^log2_enrichment, digits), "x")
}

lookup_summary <- function(quantity) {
  value <- discovery_summary$value[discovery_summary$quantity == quantity]
  if (length(value) != 1L) {
    stop("Discovery summary lookup failed for: ", quantity)
  }
  value
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
    '<div class="variant-enrichment-table-scroll">\n',
    as.character(table_html),
    "\n</div>\n",
    sep = ""
  )
  invisible(table_html)
}

core_annotations <- c(
  "GENCODE promoter +/-2 kb",
  "GENCODE exon",
  "GENCODE intron",
  "GENCODE intergenic",
  "ENCODE cCRE promoter-like",
  "ENCODE cCRE enhancer-like",
  "ENCODE cCRE CTCF-only",
  "Roadmap E020 iPS-20b: Active TSS",
  "Roadmap E020 iPS-20b: Enhancer",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Active TSS",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer",
  "Roadmap E095 left ventricle: Active TSS",
  "Roadmap E095 left ventricle: Enhancer"
)

core_regulatory_annotations <- c(
  "GENCODE promoter +/-2 kb",
  "ENCODE cCRE promoter-like",
  "ENCODE cCRE enhancer-like",
  "Roadmap E020 iPS-20b: Active TSS",
  "Roadmap E020 iPS-20b: Enhancer",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Active TSS",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer"
)

annotation_display <- c(
  "GENCODE promoter +/-2 kb" = "Promoter (+/-2 kb)",
  "GENCODE exon" = "Exon",
  "GENCODE intron" = "Intron",
  "GENCODE intergenic" = "Intergenic",
  "ENCODE cCRE promoter-like" = "Promoter-like cCRE",
  "ENCODE cCRE enhancer-like" = "Enhancer-like cCRE",
  "ENCODE cCRE CTCF-only" = "CTCF-only cCRE",
  "Roadmap E020 iPS-20b: Active TSS" = "E020 active TSS",
  "Roadmap E020 iPS-20b: Enhancer" = "E020 enhancer",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Active TSS" =
    "E013 active TSS",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer" = "E013 enhancer",
  "Roadmap E095 left ventricle: Active TSS" = "E095 active TSS",
  "Roadmap E095 left ventricle: Enhancer" = "E095 enhancer"
)

method_colors <- c(
  "FASH" = "#0072B2",
  "Strober linear" = "#D55E00",
  "Strober nonlinear" = "#009E73",
  "FASH exact-specific" = "#56B4E9",
  "FASH LD-specific" = "#0072B2",
  "Strober exact-specific" = "#E69F00",
  "Strober LD-specific" = "#D55E00"
)

prepare_plot_data <- function(panel, annotations) {
  data <- enrichment_results[
    enrichment_results$panel == panel &
      enrichment_results$annotation %in% annotations &
      is.finite(enrichment_results$log2_enrichment) &
      is.finite(enrichment_results$ci_lower) &
      is.finite(enrichment_results$ci_upper),
    ,
    drop = FALSE
  ]
  data$annotation_label <- unname(annotation_display[data$annotation])
  data$annotation_label <- factor(
    data$annotation_label,
    levels = rev(unname(annotation_display[annotations]))
  )
  data$method <- factor(
    data$method,
    levels = unique(set_metadata$method[set_metadata$panel == panel])
  )
  data
}

plot_enrichment_forest <- function(panel,
                                   annotations = core_annotations,
                                   title = panel,
                                   subtitle = NULL) {
  data <- prepare_plot_data(panel, annotations)
  dodge <- ggplot2::position_dodge(width = 0.58)
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = log2_enrichment,
      y = annotation_label,
      color = method
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.45,
      linetype = "dashed",
      color = "grey45"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
      orientation = "y",
      position = dodge,
      width = 0.34,
      linewidth = 0.55
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        shape = q_value_within_panel < 0.05,
        fill = method
      ),
      position = dodge,
      size = 2.65,
      stroke = 0.55
    ) +
    ggplot2::scale_color_manual(values = method_colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = method_colors, drop = FALSE) +
    ggplot2::scale_shape_manual(
      values = c(`TRUE` = 21, `FALSE` = 1),
      breaks = c("TRUE", "FALSE"),
      labels = c(`TRUE` = "Panel BH q < 0.05", `FALSE` = "q >= 0.05")
    ) +
    ggplot2::guides(
      shape = ggplot2::guide_legend(
        override.aes = list(
          color = "grey25",
          fill = c("grey25", NA_character_)
        )
      )
    ) +
    ggplot2::labs(
      x = "log2 enrichment relative to matched tested variants",
      y = NULL,
      color = NULL,
      fill = NULL,
      shape = NULL,
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      plot.title.position = "plot",
      axis.text.y = ggplot2::element_text(size = 9.5)
    )
}

plot_primary_enrichment <- function() {
  plot_enrichment_forest(
    panel = "Primary reported-threshold sets",
    annotations = core_annotations,
    title = "Functional annotation enrichment at the reported thresholds",
    subtitle = paste(
      "Points are pooled across 100 matched-control seeds; bars are",
      "delete-one-chromosome 95% intervals."
    )
  )
}

plot_specific_enrichment <- function() {
  plot_enrichment_forest(
    panel = "Method-specific sensitivity",
    annotations = core_annotations,
    title = "Does FASH-specific enrichment survive LD-aware reclassification?",
    subtitle = paste(
      "Exact-specific and within-study LD-specific definitions are shown;",
      "the LD estimate uses 19 donors."
    )
  )
}

plot_design_sensitivities <- function() {
  data <- enrichment_results[
    enrichment_results$panel %in% c(
      "Equal top-K sensitivity",
      "One lead variant per gene"
    ) & enrichment_results$annotation %in% core_regulatory_annotations &
      is.finite(enrichment_results$log2_enrichment),
    ,
    drop = FALSE
  ]
  data$annotation_label <- unname(annotation_display[data$annotation])
  data$annotation_label <- factor(
    data$annotation_label,
    levels = rev(unname(annotation_display[core_regulatory_annotations]))
  )
  dodge <- ggplot2::position_dodge(width = 0.58)
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = log2_enrichment,
      y = annotation_label,
      color = method
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.45,
      linetype = "dashed",
      color = "grey45"
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
      orientation = "y",
      position = dodge,
      width = 0.34,
      linewidth = 0.55
    ) +
    ggplot2::geom_point(position = dodge, size = 2.5) +
    ggplot2::facet_wrap(~panel, ncol = 1L, scales = "free_y") +
    ggplot2::scale_color_manual(values = method_colors, drop = FALSE) +
    ggplot2::labs(
      x = "log2 enrichment relative to matched tested variants",
      y = NULL,
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 9.3)
    )
}

plot_matching_balance <- function() {
  data <- merge(
    matching_balance,
    set_metadata[, c("discovery_set", "panel", "method")],
    by = "discovery_set",
    all.x = TRUE,
    sort = FALSE
  )
  data$covariate <- factor(
    data$covariate,
    levels = c(
      "minor_allele_frequency",
      "minimum_target_tss_distance",
      "local_tested_variant_count_1mb",
      "n_tested_genes"
    ),
    labels = c(
      "Cohort MAF",
      "Minimum target-TSS distance",
      "Local tested-variant density",
      "Number of tested genes"
    )
  )
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = covariate,
      y = standardized_mean_difference,
      color = method
    )
  ) +
    ggplot2::geom_hline(
      yintercept = c(-0.1, 0.1),
      linetype = "dotted",
      color = "grey45"
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      width = 0.62,
      linewidth = 0.45,
      position = ggplot2::position_dodge(width = 0.72)
    ) +
    ggplot2::facet_wrap(~panel, scales = "free_x") +
    ggplot2::scale_color_manual(values = method_colors, drop = FALSE) +
    ggplot2::labs(
      x = NULL,
      y = "Standardized mean difference",
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 10.5) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

signal_sentence <- function(discovery_set, label) {
  data <- enrichment_results[
    enrichment_results$discovery_set == discovery_set &
      enrichment_results$annotation %in% core_regulatory_annotations &
      is.finite(enrichment_results$log2_enrichment),
    ,
    drop = FALSE
  ]
  supported <- data[
    data$log2_enrichment > 0 & data$q_value_within_panel < 0.05,
    ,
    drop = FALSE
  ]
  if (!nrow(supported)) {
    return(paste0(
      label,
      " has no positive core-regulatory category with panel-level BH q < 0.05."
    ))
  }
  strongest <- supported[which.max(supported$log2_enrichment), , drop = FALSE]
  paste0(
    label,
    " is positively enriched in ",
    nrow(supported),
    " of the seven core-regulatory categories; the largest estimate is ",
    annotation_display[strongest$annotation],
    " (",
    format_fold(strongest$log2_enrichment),
    ", 95% interval ",
    format_fold(strongest$ci_lower),
    " to ",
    format_fold(strongest$ci_upper),
    ")."
  )
}

displayed_signal_sentence <- function(discovery_set, label) {
  data <- enrichment_results[
    enrichment_results$discovery_set == discovery_set &
      enrichment_results$annotation %in% core_annotations &
      is.finite(enrichment_results$log2_enrichment),
    ,
    drop = FALSE
  ]
  supported <- data[
    data$log2_enrichment > 0 & data$q_value_within_panel < 0.05,
    ,
    drop = FALSE
  ]
  if (!nrow(supported)) {
    return(paste0(
      label,
      " has no positive displayed annotation with panel-level BH q < 0.05."
    ))
  }
  strongest <- supported[which.max(supported$log2_enrichment), , drop = FALSE]
  paste0(
    "Across the broader displayed annotations, ",
    label,
    " has ",
    nrow(supported),
    " positive category with panel-level BH q < 0.05: ",
    annotation_display[strongest$annotation],
    " (",
    format_fold(strongest$log2_enrichment),
    ", 95% interval ",
    format_fold(strongest$ci_lower),
    " to ",
    format_fold(strongest$ci_upper),
    ", q = ",
    format_decimal(strongest$q_value_within_panel, 3L),
    ")."
  )
}

primary_signal_sentence <- signal_sentence(
  "primary_fash",
  "The all-FASH discovery set"
)
specific_signal_sentence <- signal_sentence(
  "specific_fash_ld",
  "The LD-aware FASH-specific discovery set"
)
comparison_signal_sentence <- signal_sentence(
  "primary_nonlinear",
  "For comparison, the primary Strober-nonlinear set"
)
specific_displayed_signal_sentence <- displayed_signal_sentence(
  "specific_fash_ld",
  "the LD-aware FASH-specific set"
)

key_result_table <- enrichment_results[
  enrichment_results$discovery_set %in% c(
    "primary_fash", "specific_fash_exact", "specific_fash_ld"
  ) & enrichment_results$annotation %in% c(
    "GENCODE intron",
    core_regulatory_annotations
  ),
  c(
    "method", "n_variants", "annotation", "enrichment", "ci_lower",
    "ci_upper", "q_value_within_panel"
  ),
  drop = FALSE
]
key_result_table$annotation <- unname(annotation_display[
  key_result_table$annotation
])
key_result_table$ci_lower_fold <- 2^key_result_table$ci_lower
key_result_table$ci_upper_fold <- 2^key_result_table$ci_upper
key_result_table <- key_result_table[, c(
  "method", "n_variants", "annotation", "enrichment",
  "ci_lower_fold", "ci_upper_fold", "q_value_within_panel"
)]
names(key_result_table) <- c(
  "Set", "Variants", "Annotation", "Fold enrichment",
  "95% CI lower", "95% CI upper", "Panel BH q"
)

matching_balance_summary <- aggregate(
  abs(matching_balance$standardized_mean_difference),
  by = list(
    discovery_set = matching_balance$discovery_set,
    covariate = matching_balance$covariate
  ),
  FUN = max
)
names(matching_balance_summary)[3L] <- "maximum_absolute_SMD"
matching_balance_summary <- merge(
  matching_balance_summary,
  set_metadata[, c("discovery_set", "panel", "method")],
  by = "discovery_set",
  all.x = TRUE,
  sort = FALSE
)
matching_balance_summary <- matching_balance_summary[, c(
  "panel", "method", "covariate", "maximum_absolute_SMD"
)]
names(matching_balance_summary) <- c(
  "Panel", "Set", "Covariate", "Maximum absolute SMD"
)

relaxation_summary <- merge(
  matching_relaxation_audit,
  set_metadata[, c("discovery_set", "panel", "method")],
  by = "discovery_set",
  all.x = TRUE,
  sort = FALSE
)
relaxation_summary$relaxation_level <- as.integer(
  as.character(relaxation_summary$relaxation_level)
)
relaxation_summary <- relaxation_summary[, c(
  "panel", "method", "relaxation_level", "Freq", "proportion"
)]
names(relaxation_summary) <- c(
  "Panel", "Set", "Relaxation level", "Variants", "Proportion"
)

source_table <- annotation_download_manifest[, c(
  "source_id", "source_name", "assembly", "accession", "byte_size", "md5"
)]
names(source_table) <- c(
  "Source ID", "Source", "Assembly", "Accession", "Bytes", "MD5"
)

reproduction_command <- paste(
  "Rscript --vanilla",
  paste0(
    "code/revision_simulations/internal/variant_annotation_enrichment/",
    "run_variant_annotation_enrichment.R"
  )
)
