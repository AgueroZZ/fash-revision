# Reporting helpers for the internal FASH versus FASH-CL enrichment page.

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
  "fash_cl_variant_enrichment_comparison",
  "analysis_cache.rds"
)
if (!file.exists(cache_path)) {
  stop(
    "The retained comparison cache is missing. Run ",
    "code/revision_simulations/internal/",
    "fash_cl_variant_enrichment_comparison/",
    "run_fash_cl_variant_enrichment_comparison.R first."
  )
}
cache <- readRDS(cache_path)
required_fields <- c(
  "configuration", "input_provenance", "bf_validation",
  "discovery_set_summary", "method_overlap", "set_metadata",
  "coverage_by_annotation_system", "enrichment_results",
  "enhancer_results", "enhancer_summary", "matching_balance",
  "matching_relaxation", "runtime_summary", "maximum_absolute_smd",
  "session_info"
)
expected_sets <- c(
  "current_all", "current_one_lead", "fash_cl_all", "fash_cl_one_lead"
)
if (!all(required_fields %in% names(cache)) ||
    !identical(
      cache$configuration$analysis_id,
      "revision_internal_fash_cl_variant_enrichment_comparison"
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
bf_validation <- cache$bf_validation
discovery_set_summary <- cache$discovery_set_summary
method_overlap <- cache$method_overlap
set_metadata <- cache$set_metadata
coverage_by_annotation_system <- cache$coverage_by_annotation_system
enrichment_results <- cache$enrichment_results
enhancer_results <- cache$enhancer_results
enhancer_summary <- cache$enhancer_summary
matching_balance <- cache$matching_balance
matching_relaxation <- cache$matching_relaxation
runtime_summary <- cache$runtime_summary
input_provenance <- cache$input_provenance

if (nrow(enrichment_results) != 4L * (29L + 83L) ||
    nrow(enhancer_results) != 4L * (4L + 5L) ||
    nrow(matching_balance) != 2L * 4L * 100L * 4L ||
    !setequal(unique(enrichment_results$annotation_system),
              configuration$annotation_systems) ||
    any(!is.finite(matching_balance$standardized_mean_difference)) ||
    bf_validation$stored_lfdr_reconstruction_maximum_absolute_difference >=
      1e-12) {
  stop("The retained comparison cache failed result validation.")
}

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
    '<div class="fash-cl-table-scroll">\n',
    as.character(table_html),
    "\n</div>\n",
    sep = ""
  )
  invisible(table_html)
}

method_colors <- c("Current FASH" = "#0072B2", "FASH-CL" = "#D55E00")
strategy_levels <- c("All discoveries", "One lead variant per gene")
system_levels <- c("Custom regulatory", "baselineLD v2.2")
enhancer_labels <- c(
  "ENCODE cCRE enhancer-like" = "ENCODE enhancer-like cCRE",
  "Roadmap E020 iPS-20b: Enhancer" = "E020 enhancer",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer" = "E013 enhancer",
  "Roadmap E095 left ventricle: Enhancer" = "E095 enhancer",
  "Enhancer_Andersson" = "Enhancer (Andersson)",
  "Enhancer_Hoffman" = "Enhancer (Hoffman)",
  "WeakEnhancer_Hoffman" = "Weak enhancer (Hoffman)",
  "SuperEnhancer_Hnisz" = "Super-enhancer (Hnisz)",
  "Human_Enhancer_Villar" = "Human enhancer (Villar)"
)

prepare_enhancer_plot_data <- function() {
  data <- enhancer_results
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
  data$method <- factor(data$method, levels = names(method_colors))
  data
}

plot_enhancer_enrichment <- function() {
  data <- prepare_enhancer_plot_data()
  finite <- data[
    is.finite(data$log2_enrichment) & is.finite(data$ci_lower) &
      is.finite(data$ci_upper),
    ,
    drop = FALSE
  ]
  dodge <- ggplot2::position_dodge(width = 0.58)
  ggplot2::ggplot(
    finite,
    ggplot2::aes(x = log2_enrichment, y = annotation_label, color = method)
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.45
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lower, xmax = ci_upper),
      orientation = "y",
      width = 0.36,
      linewidth = 0.52,
      position = dodge
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        shape = q_value_within_set < 0.05,
        fill = method
      ),
      size = 2.7,
      stroke = 0.6,
      position = dodge
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(annotation_system),
      cols = ggplot2::vars(selection_strategy),
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_color_manual(values = method_colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = method_colors, drop = FALSE) +
    ggplot2::scale_shape_manual(
      values = c(`TRUE` = 21, `FALSE` = 1),
      labels = c(`TRUE` = "Within-set BH q < 0.05", `FALSE` = "q >= 0.05")
    ) +
    ggplot2::labs(
      x = "log2 enrichment relative to matched tested variants",
      y = NULL,
      color = NULL,
      fill = NULL,
      shape = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 9.6)
    )
}

prepare_method_scatter_data <- function() {
  key <- c("annotation_system", "selection_strategy", "annotation")
  current <- enrichment_results[
    enrichment_results$method == "Current FASH",
    c(key, "log2_enrichment"),
    drop = FALSE
  ]
  fash_cl <- enrichment_results[
    enrichment_results$method == "FASH-CL",
    c(key, "log2_enrichment"),
    drop = FALSE
  ]
  data <- merge(
    current,
    fash_cl,
    by = key,
    suffixes = c("_current", "_fash_cl"),
    sort = FALSE
  )
  data$is_enhancer <- data$annotation %in% names(enhancer_labels)
  data$annotation_system <- factor(data$annotation_system, levels = system_levels)
  data$selection_strategy <- factor(
    data$selection_strategy,
    levels = strategy_levels
  )
  data
}

plot_method_scatter <- function() {
  data <- prepare_method_scatter_data()
  finite <- data[
    is.finite(data$log2_enrichment_current) &
      is.finite(data$log2_enrichment_fash_cl),
    ,
    drop = FALSE
  ]
  ggplot2::ggplot(
    finite,
    ggplot2::aes(
      x = log2_enrichment_current,
      y = log2_enrichment_fash_cl
    )
  ) +
    ggplot2::geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.5
    ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey82", linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0, color = "grey82", linewidth = 0.35) +
    ggplot2::geom_point(
      ggplot2::aes(color = is_enhancer),
      alpha = 0.78,
      size = 2.1
    ) +
    ggplot2::scale_color_manual(
      values = c(`TRUE` = "#D55E00", `FALSE` = "grey55"),
      labels = c(`TRUE` = "Pre-specified enhancer", `FALSE` = "Other")
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(annotation_system),
      cols = ggplot2::vars(selection_strategy),
      scales = "free"
    ) +
    ggplot2::labs(
      x = "Current FASH log2 enrichment",
      y = "FASH-CL log2 enrichment",
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold")
    )
}

plot_matching_balance <- function() {
  data <- matching_balance
  data <- merge(
    data,
    set_metadata[, c("discovery_set", "method", "selection_strategy")],
    by = "discovery_set",
    all.x = TRUE,
    sort = FALSE
  )
  covariate_labels <- c(
    "minor_allele_frequency" = "Cohort MAF",
    "minimum_target_tss_distance" = "Minimum target-TSS distance",
    "local_tested_variant_count_1mb" = "Local tested-variant density",
    "n_tested_genes" = "Number of tested genes"
  )
  data$covariate_label <- unname(covariate_labels[data$covariate])
  data$covariate_label <- factor(
    data$covariate_label,
    levels = unname(covariate_labels)
  )
  data$method <- factor(data$method, levels = names(method_colors))
  data$selection_strategy <- factor(
    data$selection_strategy,
    levels = strategy_levels
  )
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = standardized_mean_difference,
      y = covariate_label,
      color = method
    )
  ) +
    ggplot2::geom_vline(
      xintercept = c(-0.1, 0.1),
      linetype = "dotted",
      color = "grey45",
      linewidth = 0.45
    ) +
    ggplot2::geom_vline(xintercept = 0, color = "grey80", linewidth = 0.4) +
    ggplot2::geom_boxplot(
      orientation = "y",
      outlier.shape = NA,
      width = 0.58,
      position = ggplot2::position_dodge(width = 0.68)
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(annotation_system),
      cols = ggplot2::vars(selection_strategy)
    ) +
    ggplot2::scale_color_manual(values = method_colors, drop = FALSE) +
    ggplot2::labs(
      x = "Standardized mean difference",
      y = NULL,
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11.2) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold")
    )
}

discovery_table <- discovery_set_summary
discovery_table$selection_strategy <- ifelse(
  discovery_table$selection_strategy == "all",
  "All discoveries",
  "One lead variant per gene"
)
names(discovery_table) <- c(
  "Method", "Selection strategy", "Gene-variant pairs", "Unique variants",
  "Unique genes"
)

overlap_table <- method_overlap
overlap_table$selection_strategy <- ifelse(
  overlap_table$selection_strategy == "all",
  "All discoveries",
  "One lead variant per gene"
)
overlap_table$jaccard <- format_decimal(overlap_table$jaccard, 3L)
names(overlap_table) <- c(
  "Selection strategy", "Current variants", "FASH-CL variants",
  "Intersection", "Union", "Jaccard"
)

coverage_table <- merge(
  coverage_by_annotation_system,
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
  "annotation_system", "method", "selection_strategy",
  "original_variant_count", "covered_variant_count", "coverage_proportion"
)]
names(coverage_table) <- c(
  "Annotation system", "Method", "Selection strategy", "Original variants",
  "Covered variants", "Coverage"
)

enhancer_summary_table <- enhancer_summary
enhancer_summary_table$median_fold_enrichment <-
  2^enhancer_summary_table$median_log2_enrichment
enhancer_summary_table <- enhancer_summary_table[, c(
  "annotation_system", "method", "selection_strategy",
  "finite_annotation_count", "positive_annotation_count",
  "median_fold_enrichment", "maximum_enrichment", "minimum_p_value",
  "minimum_q_value"
)]
names(enhancer_summary_table) <- c(
  "Annotation system", "Method", "Selection strategy", "Finite enhancers",
  "Positive enhancers", "Median fold enrichment", "Maximum fold enrichment",
  "Minimum nominal p", "Minimum within-set BH q"
)

enhancer_numeric_table <- enhancer_results[, c(
  "annotation_system", "method", "selection_strategy", "annotation",
  "selected_overlap", "selected_total", "enrichment", "ci_lower", "ci_upper",
  "p_value", "q_value_within_set"
)]
enhancer_numeric_table$annotation <- unname(
  enhancer_labels[enhancer_numeric_table$annotation]
)
enhancer_numeric_table$fold_ci_lower <- 2^enhancer_numeric_table$ci_lower
enhancer_numeric_table$fold_ci_upper <- 2^enhancer_numeric_table$ci_upper
enhancer_numeric_table <- enhancer_numeric_table[, c(
  "annotation_system", "method", "selection_strategy", "annotation",
  "selected_overlap", "selected_total", "enrichment", "fold_ci_lower",
  "fold_ci_upper", "p_value", "q_value_within_set"
)]
names(enhancer_numeric_table) <- c(
  "Annotation system", "Method", "Selection strategy", "Enhancer annotation",
  "Selected overlap", "Selected total", "Fold enrichment", "Fold CI lower",
  "Fold CI upper", "Nominal p", "Within-set BH q"
)

full_results_table <- enrichment_results[, c(
  "annotation_system", "method", "selection_strategy", "annotation",
  "selected_overlap", "selected_total", "enrichment", "log2_enrichment",
  "ci_lower", "ci_upper", "p_value", "q_value_within_set"
)]
names(full_results_table) <- c(
  "Annotation system", "Method", "Selection strategy", "Annotation",
  "Selected overlap", "Selected total", "Fold enrichment", "log2 enrichment",
  "log2 CI lower", "log2 CI upper", "Nominal p", "Within-set BH q"
)

runtime_table <- runtime_summary
runtime_table$elapsed_seconds <- format_decimal(runtime_table$elapsed_seconds, 1L)
names(runtime_table) <- c("Stage", "Elapsed seconds")

bf_validation_table <- data.frame(
  Diagnostic = c(
    "Stored current lfdr reconstruction: maximum absolute difference",
    "Current optimizer rerun: null-weight absolute difference",
    "Current optimizer rerun: mean lfdr absolute difference",
    "Current optimizer rerun: lfdr correlation",
    "Current optimizer rerun: discovery-call Jaccard"
  ),
  Value = c(
    format(bf_validation$stored_lfdr_reconstruction_maximum_absolute_difference,
           scientific = TRUE, digits = 4L),
    format(bf_validation$null_weight_absolute_difference,
           scientific = TRUE, digits = 4L),
    format(bf_validation$mean_lfdr_absolute_difference,
           scientific = TRUE, digits = 4L),
    format_decimal(bf_validation$lfdr_correlation, 7L),
    format_decimal(bf_validation$discovery_call_jaccard, 6L)
  ),
  stringsAsFactors = FALSE
)

summary_lookup <- function(system, method, strategy) {
  row <- enhancer_summary[
    enhancer_summary$annotation_system == system &
      enhancer_summary$method == method &
      enhancer_summary$selection_strategy == strategy,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L) {
    stop("Enhancer summary lookup failed.")
  }
  row
}

cl_custom_all <- summary_lookup(
  "Custom regulatory", "FASH-CL", "All discoveries"
)
cl_custom_lead <- summary_lookup(
  "Custom regulatory", "FASH-CL", "One lead variant per gene"
)
cl_baseline_all <- summary_lookup(
  "baselineLD v2.2", "FASH-CL", "All discoveries"
)
cl_baseline_lead <- summary_lookup(
  "baselineLD v2.2", "FASH-CL", "One lead variant per gene"
)
current_custom_lead <- summary_lookup(
  "Custom regulatory", "Current FASH", "One lead variant per gene"
)
current_baseline_lead <- summary_lookup(
  "baselineLD v2.2", "Current FASH", "One lead variant per gene"
)

lead_pattern_answer <- paste0(
  "For FASH-CL, one lead variant per gene has the larger median enhancer ",
  "point estimate in both annotation systems: ",
  format_fold(cl_custom_all$median_log2_enrichment), " to ",
  format_fold(cl_custom_lead$median_log2_enrichment),
  " in the custom panel and ",
  format_fold(cl_baseline_all$median_log2_enrichment), " to ",
  format_fold(cl_baseline_lead$median_log2_enrichment),
  " in baselineLD."
)

method_comparison_answer <- paste0(
  "FASH-CL has larger median enhancer point estimates for all discoveries ",
  "in both systems. Under one-lead-per-gene selection, FASH-CL is larger in ",
  "baselineLD (", format_fold(cl_baseline_lead$median_log2_enrichment),
  " versus ", format_fold(current_baseline_lead$median_log2_enrichment),
  "), whereas current FASH is slightly larger in the custom panel (",
  format_fold(current_custom_lead$median_log2_enrichment), " versus ",
  format_fold(cl_custom_lead$median_log2_enrichment), ")."
)

smallest_enhancer_q <- min(enhancer_results$q_value_within_set, na.rm = TRUE)
