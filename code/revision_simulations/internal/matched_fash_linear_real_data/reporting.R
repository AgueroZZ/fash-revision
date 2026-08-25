# Reporting helpers for the matched FASH-linear real-data internal page.

find_workflowr_root_reporting <- function(start = getwd()) {
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

require_reporting_package_matched <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The ", package, " package is required to render this page.")
  }
  invisible(TRUE)
}

invisible(lapply(
  c("ggplot2", "ggVennDiagram", "knitr", "scales"),
  require_reporting_package_matched
))

workflowr_root <- find_workflowr_root_reporting()
analysis_id <- "matched_fash_linear_real_data_fashr_0_1_43"
output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal", analysis_id
)
cache_path <- file.path(output_directory, "analysis_cache.rds")
run_status_path <- file.path(output_directory, "run_status.rds")
linear_raw_path <- file.path(output_directory, "linear_fit_raw.rds")
linear_bf_path <- file.path(output_directory, "linear_fit_bf.rds")
overlap_cache_path <- file.path(
  output_directory,
  "iwp_overlap_proportions.rds"
)
required_paths <- c(
  cache_path,
  run_status_path,
  linear_raw_path,
  linear_bf_path,
  overlap_cache_path
)
if (any(!file.exists(required_paths))) {
  stop("The versioned matched FASH-linear cache or fit artifact is missing.")
}

cache <- readRDS(cache_path)
expected_cache_fields <- c(
  "configuration", "package_provenance", "input_provenance",
  "fit_provenance", "matched_settings_table", "discovery_counts",
  "prior_weights", "venn_sets", "venn_region_counts",
  "pairwise_overlap", "validation", "runtime_summary"
)
if (!identical(names(cache), expected_cache_fields)) {
  stop("The matched FASH-linear cache failed structural validation.")
}

configuration <- cache$configuration
expected_grid <- sort(c(0, exp(-0.5 * seq(0, 10, by = 0.2))))
if (!identical(configuration$analysis_id, analysis_id) ||
    !identical(configuration$package_version, "0.1.43") ||
    !identical(
      configuration$package_sha,
      "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
    ) ||
    !identical(configuration$pc_correction, "Time-specific PCs") ||
    !identical(configuration$expected_pairs, 1009173L) ||
    !identical(configuration$time_grid, 0:15) ||
    !identical(configuration$likelihood, "gaussian") ||
    !identical(configuration$betaprec, 0) ||
    !identical(configuration$reference_iwp_order, 1L) ||
    !identical(configuration$reference_iwp_num_basis, 20L) ||
    !identical(configuration$pred_step, 1) ||
    !identical(configuration$penalty, 10L) ||
    !isTRUE(all.equal(configuration$grid, expected_grid, tolerance = 0)) ||
    !identical(configuration$alpha, 0.05) ||
    !identical(configuration$strober_rule, "eFDR <= 0.05")) {
  stop("The matched FASH-linear cache failed configuration validation.")
}

run_status <- readRDS(run_status_path)
if (!is.list(run_status) ||
    !identical(run_status$analysis_id, analysis_id) ||
    !identical(run_status$status, "complete") ||
    !identical(
      run_status$analysis_cache_md5,
      unname(tools::md5sum(cache_path))
    ) ||
    !identical(run_status$package_sha, configuration$package_sha)) {
  stop("The retained run status is incomplete or stale.")
}

fit_paths <- stats::setNames(
  c(linear_raw_path, linear_bf_path),
  c("linear_raw", "linear_bf")
)
current_fit_md5 <- unname(tools::md5sum(fit_paths))
expected_fit_md5 <- cache$fit_provenance$md5[match(
  names(fit_paths),
  cache$fit_provenance$label
)]
if (anyNA(expected_fit_md5) || !identical(current_fit_md5, expected_fit_md5)) {
  stop("A retained FASH-linear fit changed after the cache was created.")
}

current_input_md5 <- unname(tools::md5sum(cache$input_provenance$path))
if (anyNA(current_input_md5) ||
    !identical(current_input_md5, cache$input_provenance$md5)) {
  stop("A matched FASH-linear input changed after the rerun.")
}

validation <- cache$validation
if (!is.data.frame(validation) ||
    !identical(names(validation), c("check", "passed")) ||
    nrow(validation) != 13L ||
    any(!validation$passed)) {
  stop("The retained matched FASH-linear validation did not pass.")
}

package_provenance <- cache$package_provenance
matched_settings_table <- cache$matched_settings_table
discovery_counts <- cache$discovery_counts
prior_weights <- cache$prior_weights
venn_sets <- cache$venn_sets
venn_region_counts <- cache$venn_region_counts
pairwise_overlap <- cache$pairwise_overlap
runtime_summary <- cache$runtime_summary

expected_units <- c("Gene-variant pairs", "Genes", "Variants")
expected_methods <- c(
  "Strober quadratic", "Strober linear", "FASH-IWP1 BF", "FASH-linear BF"
)
if (!identical(names(venn_sets), expected_units) ||
    any(!vapply(venn_sets, function(sets) {
      identical(names(sets), expected_methods) &&
        all(lengths(sets) > 0L) &&
        all(vapply(sets, function(values) !anyDuplicated(values), logical(1)))
    }, logical(1)))) {
  stop("The retained four-method Venn sets are invalid.")
}

expected_discovery_rows <- data.frame(
  method = c("FASH-IWP1", "FASH-IWP1", "FASH-linear", "FASH-linear"),
  adjustment = c("Raw", "BF-adjusted", "Raw", "BF-adjusted"),
  pair_count = c(43860L, 9214L, 60188L, 14902L),
  gene_count = c(3258L, 1176L, 3863L, 1663L),
  variant_count = c(42893L, 9148L, 58570L, 14761L),
  stringsAsFactors = FALSE
)
observed_discovery_rows <- discovery_counts[, names(expected_discovery_rows)]
rownames(observed_discovery_rows) <- NULL
if (!identical(observed_discovery_rows, expected_discovery_rows)) {
  stop("The retained discovery-count invariants changed.")
}

overlap_cache <- readRDS(overlap_cache_path)
expected_overlap_fields <- c(
  "analysis_id", "created_at", "configuration", "input_provenance",
  "iwp1_linear", "iwp2_linear"
)
if (!identical(names(overlap_cache), expected_overlap_fields) ||
    !identical(
      overlap_cache$analysis_id,
      "fash_iwp_linear_directional_overlap_fashr_0_1_43"
    ) ||
    !identical(overlap_cache$configuration$package_version, "0.1.43") ||
    !identical(
      overlap_cache$configuration$package_sha,
      "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
    ) ||
    !identical(overlap_cache$configuration$adjustment, "BF-adjusted") ||
    !identical(
      overlap_cache$configuration$discovery_rule,
      "cumulative-lfdr FDR 0.05"
    ) ||
    !identical(overlap_cache$configuration$alpha, 0.05)) {
  stop("The directional overlap cache failed structural validation.")
}
current_overlap_input_md5 <- unname(tools::md5sum(
  overlap_cache$input_provenance$path
))
if (anyNA(current_overlap_input_md5) ||
    !identical(
      current_overlap_input_md5,
      overlap_cache$input_provenance$md5
    )) {
  stop("A directional overlap input changed after cache construction.")
}

expected_iwp1_linear <- data.frame(
  unit = expected_units,
  reference_count = c(14902L, 1663L, 14761L),
  comparison_count = c(9214L, 1176L, 9148L),
  intersection_count = c(8530L, 1112L, 8468L),
  comparison_covered_by_reference = c(
    8530 / 9214,
    1112 / 1176,
    8468 / 9148
  ),
  stringsAsFactors = FALSE
)
expected_iwp2_linear <- data.frame(
  unit = expected_units,
  reference_count = c(14902L, 1663L, 14761L),
  comparison_count = c(44L, 9L, 44L),
  intersection_count = c(21L, 6L, 21L),
  comparison_covered_by_reference = c(21 / 44, 6 / 9, 21 / 44),
  stringsAsFactors = FALSE
)
if (!identical(overlap_cache$iwp1_linear, expected_iwp1_linear) ||
    !identical(overlap_cache$iwp2_linear, expected_iwp2_linear)) {
  stop("The directional overlap-count invariants changed.")
}

format_integer_matched <- function(x) {
  format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

format_decimal_matched <- function(x, digits = 3L) {
  formatC(x, format = "f", digits = digits)
}

render_scrollable_table_matched <- function(data,
                                            caption,
                                            digits = 4L,
                                            minimum_width = "760px") {
  formatted <- data
  numeric_columns <- vapply(formatted, is.numeric, logical(1))
  formatted[numeric_columns] <- lapply(formatted[numeric_columns], function(values) {
    ifelse(
      is.na(values),
      "",
      formatC(values, format = "fg", digits = digits, big.mark = ",")
    )
  })
  html <- knitr::kable(
    formatted,
    format = "html",
    escape = TRUE,
    row.names = FALSE,
    caption = caption,
    table.attr = paste0(
      "class='table table-striped table-condensed' style='min-width:",
      minimum_width,
      ";'"
    )
  )
  knitr::asis_output(paste0(
    "<div class='matched-table-scroll'>\n",
    html,
    "\n</div>\n"
  ))
}

discovery_counts_display <- discovery_counts[, c(
  "method", "adjustment", "pair_count", "gene_count", "variant_count",
  "estimated_pi0"
)]
names(discovery_counts_display) <- c(
  "Method", "Adjustment", "Pairs", "Genes", "Unique variants",
  "Estimated pi0"
)

iwp_linear_overlap <- pairwise_overlap[
  pairwise_overlap$method_1 == "FASH-IWP1 BF" &
    pairwise_overlap$method_2 == "FASH-linear BF",
  c(
    "unit", "method_1_count", "method_2_count", "intersection_count",
    "union_count", "jaccard"
  ),
  drop = FALSE
]
rownames(iwp_linear_overlap) <- NULL
names(iwp_linear_overlap) <- c(
  "Unit", "FASH-IWP1 BF", "FASH-linear BF", "Intersection", "Union",
  "Jaccard"
)

format_directional_overlap_matched <- function(data, comparison_label) {
  display <- data[, c(
    "unit", "reference_count", "comparison_count", "intersection_count",
    "comparison_covered_by_reference"
  )]
  display$unit <- c("Pairs", "Genes", "Variants")
  display$comparison_covered_by_reference <- scales::percent(
    display$comparison_covered_by_reference,
    accuracy = 0.1
  )
  names(display) <- c(
    "Level", "FASH-linear", comparison_label, "Intersection",
    paste0(comparison_label, " covered by linear")
  )
  display
}

iwp1_linear_proportion_display <- format_directional_overlap_matched(
  overlap_cache$iwp1_linear,
  "FASH-IWP1"
)
iwp2_linear_proportion_display <- format_directional_overlap_matched(
  overlap_cache$iwp2_linear,
  "FASH-IWP2"
)

package_display <- package_provenance[, c(
  "package", "version", "remote_ref", "remote_sha", "built"
)]
names(package_display) <- c(
  "Package", "Version", "Git ref", "Git SHA", "Built"
)

iwp_bf_counts <- discovery_counts[
  discovery_counts$method == "FASH-IWP1" &
    discovery_counts$adjustment == "BF-adjusted",
  ,
  drop = FALSE
]
linear_bf_counts <- discovery_counts[
  discovery_counts$method == "FASH-linear" &
    discovery_counts$adjustment == "BF-adjusted",
  ,
  drop = FALSE
]
