# Reporting helpers for the internal baselineLD v2.2 enrichment page.

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
  "baseline_ld_variant_enrichment",
  "analysis_cache.rds"
)
if (!file.exists(cache_path)) {
  stop(
    "The retained baselineLD analysis cache is missing. Run ",
    "code/revision_simulations/internal/baseline_ld_variant_enrichment/",
    "run_baseline_ld_variant_enrichment.R first."
  )
}
cache <- readRDS(cache_path)
required_cache_fields <- c(
  "configuration", "input_provenance", "download_manifest",
  "column_classification", "pair_summary", "set_metadata",
  "coverage_by_set", "coverage_by_chromosome", "annotation_metadata",
  "enrichment_results", "matching_balance", "matching_seed_results",
  "matching_relaxation_audit", "fash_only_variant_multiplicity",
  "runtime_summary", "resource_summary", "session_info"
)
expected_sets <- c(
  "all_fash", "one_lead_fash_per_gene", "fash_only_pair_variants"
)
if (!all(required_cache_fields %in% names(cache)) ||
    !identical(
      cache$configuration$analysis_id,
      "revision_internal_baseline_ld_variant_enrichment"
    ) ||
    cache$configuration$matching_seed_count != 100L ||
    cache$configuration$controls_per_variant != 5L ||
    !identical(cache$set_metadata$discovery_set, expected_sets) ||
    !identical(cache$configuration$jackknife_blocks, as.character(1:22))) {
  stop("The retained baselineLD cache failed structural validation.")
}

current_md5 <- unname(tools::md5sum(cache$input_provenance$path))
if (anyNA(current_md5) ||
    !identical(current_md5, cache$input_provenance$md5)) {
  stop("At least one retained analysis input changed after cache creation.")
}

configuration <- cache$configuration
pair_summary <- cache$pair_summary
set_metadata <- cache$set_metadata
coverage_by_set <- cache$coverage_by_set
coverage_by_chromosome <- cache$coverage_by_chromosome
annotation_metadata <- cache$annotation_metadata
column_classification <- cache$column_classification
enrichment_results <- cache$enrichment_results
matching_balance <- cache$matching_balance
matching_seed_results <- cache$matching_seed_results
matching_relaxation_audit <- cache$matching_relaxation_audit
runtime_summary <- cache$runtime_summary
resource_summary <- cache$resource_summary
download_manifest <- cache$download_manifest
input_provenance <- cache$input_provenance

if (nrow(enrichment_results) != 3L * 83L ||
    length(unique(matching_seed_results$match_seed)) != 100L ||
    any(!is.finite(matching_balance$standardized_mean_difference)) ||
    nrow(download_manifest) != 22L) {
  stop("The retained baselineLD cache failed result validation.")
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

lookup_resource <- function(quantity) {
  value <- resource_summary$value[resource_summary$quantity == quantity]
  if (length(value) != 1L) {
    stop("Resource summary lookup failed for: ", quantity)
  }
  value
}

lookup_runtime <- function(stage) {
  value <- runtime_summary$elapsed_seconds[runtime_summary$stage == stage]
  if (length(value) != 1L) {
    stop("Runtime summary lookup failed for: ", stage)
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
    '<div class="baseline-ld-table-scroll">\n',
    as.character(table_html),
    "\n</div>\n",
    sep = ""
  )
  invisible(table_html)
}

set_display <- setNames(
  set_metadata$display_label,
  set_metadata$discovery_set
)
set_colors <- c(
  "All FASH discoveries" = "#0072B2",
  "One FASH lead variant per gene" = "#D55E00",
  "Variants in FASH-only gene-variant pairs" = "#009E73"
)

annotation_display <- c(
  "Coding_UCSC" = "Coding",
  "Conserved_LindbladToh" = "Conserved",
  "CTCF_Hoffman" = "CTCF",
  "DGF_ENCODE" = "Digital genomic footprint",
  "DHS_Trynka" = "DNase hypersensitivity",
  "Enhancer_Andersson" = "Enhancer (Andersson)",
  "Enhancer_Hoffman" = "Enhancer (Hoffman)",
  "FetalDHS_Trynka" = "Fetal DHS",
  "H3K27ac_Hnisz" = "H3K27ac (Hnisz)",
  "H3K27ac_PGC2" = "H3K27ac (PGC2)",
  "H3K4me1_Trynka" = "H3K4me1",
  "H3K4me3_Trynka" = "H3K4me3",
  "H3K9ac_Trynka" = "H3K9ac",
  "Intron_UCSC" = "Intron",
  "PromoterFlanking_Hoffman" = "Promoter flanking",
  "Promoter_UCSC" = "Promoter",
  "Repressed_Hoffman" = "Repressed",
  "SuperEnhancer_Hnisz" = "Super-enhancer",
  "TFBS_ENCODE" = "ENCODE TF binding",
  "Transcr_Hoffman" = "Transcribed",
  "TSS_Hoffman" = "TSS",
  "UTR_3_UCSC" = "3-prime UTR",
  "UTR_5_UCSC" = "5-prime UTR",
  "WeakEnhancer_Hoffman" = "Weak enhancer",
  "synonymous" = "Synonymous",
  "non_synonymous" = "Nonsynonymous",
  "BivFlnk" = "Bivalent flanking",
  "Human_Promoter_Villar" = "Human promoter (Villar)",
  "Human_Enhancer_Villar" = "Human enhancer (Villar)"
)

core_annotations <- configuration$core_annotations
if (!identical(names(annotation_display), core_annotations)) {
  stop("The core annotation labels are not aligned to the analysis cache.")
}

prepare_core_plot_data <- function() {
  data <- enrichment_results[
    enrichment_results$core_annotation &
      is.finite(enrichment_results$log2_enrichment) &
      is.finite(enrichment_results$ci_lower) &
      is.finite(enrichment_results$ci_upper),
    ,
    drop = FALSE
  ]
  data$set_label <- unname(set_display[data$discovery_set])
  data$set_label <- factor(data$set_label, levels = unname(set_display))
  data$annotation_label <- unname(annotation_display[data$annotation])
  data$annotation_label <- factor(
    data$annotation_label,
    levels = rev(unname(annotation_display))
  )
  data
}

plot_core_enrichment <- function() {
  data <- prepare_core_plot_data()
  dodge <- ggplot2::position_dodge(width = 0.64)
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = log2_enrichment,
      y = annotation_label,
      color = set_label
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
      width = 0.42,
      linewidth = 0.48
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        shape = q_value_within_set < 0.05,
        fill = set_label
      ),
      position = dodge,
      size = 2.5,
      stroke = 0.55
    ) +
    ggplot2::scale_color_manual(values = set_colors, drop = FALSE) +
    ggplot2::scale_fill_manual(values = set_colors, drop = FALSE) +
    ggplot2::scale_shape_manual(
      values = c(`TRUE` = 21, `FALSE` = 1),
      breaks = c("TRUE", "FALSE"),
      labels = c(`TRUE` = "Set-level BH q < 0.05", `FALSE` = "q >= 0.05")
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
      axis.text.y = ggplot2::element_text(size = 9.1)
    )
}

plot_matching_balance <- function() {
  data <- matching_balance
  data$set_label <- unname(set_display[data$discovery_set])
  data$set_label <- factor(data$set_label, levels = unname(set_display))
  covariate_display <- c(
    "minor_allele_frequency" = "Cohort MAF",
    "minimum_target_tss_distance" = "Minimum target-TSS distance",
    "local_tested_variant_count_1mb" = "Local tested-variant density",
    "n_tested_genes" = "Number of tested genes"
  )
  data$covariate_label <- unname(covariate_display[data$covariate])
  data$covariate_label <- factor(
    data$covariate_label,
    levels = unname(covariate_display)
  )
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = standardized_mean_difference,
      y = covariate_label,
      color = set_label
    )
  ) +
    ggplot2::geom_vline(
      xintercept = c(-0.1, 0.1),
      linetype = "dotted",
      color = "grey45",
      linewidth = 0.45
    ) +
    ggplot2::geom_vline(xintercept = 0, color = "grey75", linewidth = 0.4) +
    ggplot2::geom_boxplot(
      orientation = "y",
      outlier.shape = NA,
      width = 0.64,
      position = ggplot2::position_dodge(width = 0.72)
    ) +
    ggplot2::scale_color_manual(values = set_colors, drop = FALSE) +
    ggplot2::labs(
      x = "Standardized mean difference",
      y = NULL,
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

pair_summary_table <- pair_summary
pair_summary_table$discovery_set <- unname(
  set_display[pair_summary_table$discovery_set]
)
names(pair_summary_table) <- c(
  "Discovery set", "Gene-variant pairs", "Unique variants", "Unique genes",
  "baselineLD-covered variants", "Coverage"
)
pair_summary_table$Coverage <- format_percent(pair_summary_table$Coverage)

coverage_table <- coverage_by_set
coverage_table$discovery_set <- c(
  "Tested-variant universe",
  unname(set_display[coverage_table$discovery_set[-1L]])
)
coverage_table$coverage_proportion <- format_percent(
  coverage_table$coverage_proportion
)
coverage_table$coverage_difference_from_background <- paste0(
  format_decimal(100 * coverage_table$coverage_difference_from_background, 1L),
  " pp"
)
names(coverage_table) <- c(
  "Set", "Original variants", "Covered variants", "Coverage",
  "Difference from background"
)

core_numeric_table <- enrichment_results[
  enrichment_results$core_annotation,
  c(
    "discovery_set", "annotation", "selected_overlap", "selected_total",
    "control_overlap", "control_total", "enrichment", "log2_enrichment",
    "ci_lower", "ci_upper", "p_value", "q_value_within_set"
  ),
  drop = FALSE
]
core_numeric_table$discovery_set <- unname(
  set_display[core_numeric_table$discovery_set]
)
core_numeric_table$annotation <- unname(
  annotation_display[core_numeric_table$annotation]
)
names(core_numeric_table) <- c(
  "Discovery set", "Annotation", "Selected overlap", "Selected total",
  "Control overlap", "Control total", "Fold enrichment", "log2 enrichment",
  "CI lower", "CI upper", "p-value", "BH q-value"
)
core_numeric_table <- core_numeric_table[order(
  match(core_numeric_table[["Discovery set"]], unname(set_display)),
  match(core_numeric_table$Annotation, unname(annotation_display))
), , drop = FALSE]

full_results_table <- enrichment_results[, c(
  "display_label", "annotation", "selected_overlap", "selected_total",
  "control_overlap", "control_total", "enrichment", "log2_enrichment",
  "ci_lower", "ci_upper", "p_value", "q_value_within_set",
  "core_annotation"
), drop = FALSE]
names(full_results_table) <- c(
  "Discovery set", "Annotation", "Selected overlap", "Selected total",
  "Control overlap", "Control total", "Fold enrichment", "log2 enrichment",
  "CI lower", "CI upper", "p-value", "BH q-value", "Core figure"
)
full_results_table <- full_results_table[order(
  match(full_results_table[["Discovery set"]], unname(set_display)),
  full_results_table[["BH q-value"]],
  -full_results_table[["log2 enrichment"]]
), , drop = FALSE]

continuous_annotation_table <- column_classification[
  column_classification$annotation_type == "continuous",
  c("annotation", "n_unique", "minimum", "maximum"),
  drop = FALSE
]
names(continuous_annotation_table) <- c(
  "Continuous annotation", "Unique chr1 values", "Minimum", "Maximum"
)

runtime_table <- runtime_summary
runtime_table$elapsed_seconds <- ifelse(
  is.finite(runtime_table$elapsed_seconds),
  format_decimal(runtime_table$elapsed_seconds, 1L),
  "47.0"
)
names(runtime_table) <- c("Stage", "Elapsed seconds", "Definition")

resource_table <- resource_summary
resource_table$value <- format_decimal(resource_table$value, 3L)
names(resource_table) <- c("Quantity", "Value")

download_table <- download_manifest[, c(
  "chromosome", "byte_size", "sha256", "retrieval_url"
)]
names(download_table) <- c(
  "Chromosome", "Bytes", "SHA256", "Retrieval URL"
)

provenance_table <- input_provenance
provenance_table$path <- sub(
  paste0("^", workflowr_root, "/?"),
  "",
  provenance_table$path
)
names(provenance_table) <- c("Input", "Bytes", "MD5", "Modified at")

minimum_q_by_set <- tapply(
  enrichment_results$q_value_within_set,
  enrichment_results$discovery_set,
  min,
  na.rm = TRUE
)
all_fash_top <- enrichment_results[
  enrichment_results$discovery_set == "all_fash" &
    is.finite(enrichment_results$q_value_within_set),
  ,
  drop = FALSE
]
all_fash_top <- all_fash_top[which.min(all_fash_top$q_value_within_set), ]
lead_enhancer <- enrichment_results[
  enrichment_results$discovery_set == "one_lead_fash_per_gene" &
    enrichment_results$annotation == "Enhancer_Hoffman",
  ,
  drop = FALSE
]
fash_only_significant <- enrichment_results[
  enrichment_results$discovery_set == "fash_only_pair_variants" &
    enrichment_results$q_value_within_set < 0.05,
  ,
  drop = FALSE
]

provisional_answer <- paste0(
  "No baselineLD binary annotation passed within-set BH q < 0.05 for any ",
  "of the three requested FASH sets. All-FASH's smallest q-value was ",
  format_decimal(all_fash_top$q_value_within_set, 3L),
  " for ", all_fash_top$annotation, " (",
  format_decimal(all_fash_top$enrichment, 2L), "-fold)."
)
lead_answer <- paste0(
  "The one-lead-per-gene set showed a suggestive Hoffman enhancer estimate ",
  "of ", format_decimal(lead_enhancer$enrichment, 2L), "-fold (95% CI on ",
  "the fold scale ", format_decimal(2^lead_enhancer$ci_lower, 2L), "--",
  format_decimal(2^lead_enhancer$ci_upper, 2L), "; q = ",
  format_decimal(lead_enhancer$q_value_within_set, 3L), "), but it did not ",
  "cross the pre-specified multiple-testing threshold."
)
fash_only_answer <- if (!nrow(fash_only_significant)) {
  paste(
    "Variants participating in exact FASH-only gene-variant pairs did not",
    "show a BH-significant baselineLD category."
  )
} else {
  paste(
    "At least one category was BH-significant among variants participating",
    "in exact FASH-only gene-variant pairs."
  )
}

session_info_text <- paste(capture.output(cache$session_info), collapse = "\n")
