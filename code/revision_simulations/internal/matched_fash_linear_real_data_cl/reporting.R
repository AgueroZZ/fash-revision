# Reporting helpers for the matched CL-PC FASH-linear internal page.

find_workflowr_root_reporting_cl <- function(start = getwd()) {
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

require_reporting_package_cl <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The ", package, " package is required to render this page.")
  }
  invisible(TRUE)
}

invisible(lapply(
  c("ggplot2", "ggVennDiagram", "knitr", "scales"),
  require_reporting_package_cl
))

workflowr_root <- find_workflowr_root_reporting_cl()
analysis_id <- "matched_fash_linear_real_data_cl_fashr_0_1_43"
output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal", analysis_id
)
cache_path <- file.path(output_directory, "analysis_cache.rds")
run_status_path <- file.path(output_directory, "run_status.rds")
fit_paths <- stats::setNames(
  file.path(output_directory, c(
    "linear_fit_raw.rds",
    "linear_fit_bf.rds",
    "cl_iwp1_bf_adjustment.rds",
    "cl_iwp2_bf_adjustment.rds"
  )),
  c("linear_raw", "linear_bf", "cl_iwp1_bf", "cl_iwp2_bf")
)
required_paths <- c(cache_path, run_status_path, fit_paths)
if (any(!file.exists(required_paths))) {
  stop("A versioned matched CL-PC cache or fit artifact is missing.")
}

cache <- readRDS(cache_path)
expected_cache_fields <- c(
  "configuration", "package_provenance", "input_provenance",
  "fit_provenance", "matched_settings_table", "discovery_counts",
  "prior_weights", "venn_sets", "venn_region_counts",
  "pairwise_overlap", "iwp1_linear", "iwp2_linear", "validation",
  "runtime_summary"
)
if (!identical(names(cache), expected_cache_fields)) {
  stop("The matched CL-PC cache failed structural validation.")
}

configuration <- cache$configuration
expected_grid <- sort(c(0, exp(-0.5 * seq(0, 10, by = 0.2))))
if (!identical(configuration$analysis_id, analysis_id) ||
    !identical(configuration$package_version, "0.1.43") ||
    !identical(
      configuration$package_sha,
      "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
    ) ||
    !identical(
      configuration$pc_correction,
      "Cell-line-collapsed PCs repeated across time"
    ) ||
    !identical(configuration$expected_pairs, 1009173L) ||
    !identical(configuration$time_grid, 0:15) ||
    !identical(configuration$likelihood, "gaussian") ||
    !identical(configuration$betaprec, 0) ||
    !identical(configuration$reference_iwp_order, 1L) ||
    !identical(configuration$reference_iwp_num_basis, 20L) ||
    !identical(configuration$pred_step, 1) ||
    !identical(configuration$penalty, 10L) ||
    !identical(configuration$grid, expected_grid) ||
    !identical(configuration$alpha, 0.05) ||
    !identical(configuration$strober_rule, "eFDR <= 0.05")) {
  stop("The matched CL-PC cache failed configuration validation.")
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
  stop("The retained matched CL-PC run status is incomplete or stale.")
}

current_fit_md5 <- unname(tools::md5sum(fit_paths))
expected_fit_md5 <- cache$fit_provenance$md5[match(
  names(fit_paths),
  cache$fit_provenance$label
)]
if (anyNA(expected_fit_md5) || !identical(current_fit_md5, expected_fit_md5)) {
  stop("A retained matched CL-PC fit changed after cache construction.")
}

current_input_md5 <- unname(tools::md5sum(cache$input_provenance$path))
if (anyNA(current_input_md5) ||
    !identical(current_input_md5, cache$input_provenance$md5)) {
  stop("A matched CL-PC input changed after cache construction.")
}

validation <- cache$validation
if (!is.data.frame(validation) ||
    !identical(names(validation), c("check", "passed")) ||
    nrow(validation) != 15L ||
    any(!validation$passed)) {
  stop("The retained matched CL-PC validation did not pass.")
}

package_provenance <- cache$package_provenance
matched_settings_table <- cache$matched_settings_table
discovery_counts <- cache$discovery_counts
prior_weights <- cache$prior_weights
venn_sets <- cache$venn_sets
venn_region_counts <- cache$venn_region_counts
pairwise_overlap <- cache$pairwise_overlap
iwp1_linear <- cache$iwp1_linear
iwp2_linear <- cache$iwp2_linear
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
  stop("The retained CL-PC four-method Venn sets are invalid.")
}

expected_discovery_rows <- data.frame(
  method = c("FASH-IWP1", "FASH-IWP1", "FASH-linear", "FASH-linear"),
  adjustment = c("Raw", "BF-adjusted", "Raw", "BF-adjusted"),
  pair_count = c(30344L, 5395L, 38634L, 9059L),
  gene_count = c(2403L, 686L, 2809L, 1012L),
  variant_count = c(29771L, 5363L, 37741L, 8989L),
  stringsAsFactors = FALSE
)
observed_discovery_rows <- discovery_counts[, names(expected_discovery_rows)]
rownames(observed_discovery_rows) <- NULL
if (!identical(observed_discovery_rows, expected_discovery_rows)) {
  stop("The retained CL-PC discovery-count invariants changed.")
}

expected_iwp1_linear <- data.frame(
  unit = expected_units,
  reference_count = c(9059L, 1012L, 8989L),
  comparison_count = c(5395L, 686L, 5363L),
  intersection_count = c(4793L, 628L, 4761L),
  comparison_covered_by_reference = c(
    4793 / 5395,
    628 / 686,
    4761 / 5363
  ),
  stringsAsFactors = FALSE
)
expected_iwp2_linear <- data.frame(
  unit = expected_units,
  reference_count = c(9059L, 1012L, 8989L),
  comparison_count = c(60L, 6L, 60L),
  intersection_count = c(6L, 3L, 6L),
  comparison_covered_by_reference = c(0.1, 0.5, 0.1),
  stringsAsFactors = FALSE
)
if (!identical(iwp1_linear, expected_iwp1_linear) ||
    !identical(iwp2_linear, expected_iwp2_linear)) {
  stop("The retained CL-PC directional-overlap invariants changed.")
}

format_integer_matched_cl <- function(x) {
  format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

format_decimal_matched_cl <- function(x, digits = 3L) {
  formatC(x, format = "f", digits = digits)
}

render_scrollable_table_matched_cl <- function(data,
                                               caption,
                                               digits = 4L,
                                               minimum_width = "760px") {
  formatted <- data
  numeric_columns <- vapply(formatted, is.numeric, logical(1))
  formatted[numeric_columns] <- lapply(
    formatted[numeric_columns],
    function(values) {
      ifelse(
        is.na(values),
        "",
        formatC(values, format = "fg", digits = digits, big.mark = ",")
      )
    }
  )
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
    "<div class='matched-cl-table-scroll'>\n",
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

format_directional_overlap_matched_cl <- function(data, comparison_label) {
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

iwp1_linear_proportion_display <- format_directional_overlap_matched_cl(
  iwp1_linear,
  "FASH-IWP1"
)
iwp2_linear_proportion_display <- format_directional_overlap_matched_cl(
  iwp2_linear,
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
