# Reporting helpers for the R6 enhancer-enrichment page.
#
# The page is a short, exploratory, orthogonal check: do the variants FASH
# discovers fall inside Roadmap enhancers of the three epigenomes that bracket
# this iPSC-to-cardiomyocyte differentiation more often than the tested-variant
# background does?
#
# Displayed scope, deliberately narrow:
#   * three Roadmap enhancer annotations (E020, E013, E095);
#   * fold enrichment against all tested variants, no matched controls;
#   * leave-one-autosome-out jackknife 95% intervals.
#
# The matched-control machinery, the ENCODE/GENCODE columns, and the baselineLD
# system all remain available through enrichment_api.R for other analyses; this
# file simply does not surface them.

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

workflowr_root <- find_workflowr_root()

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required to render this page.")
}

# load_annotation_groups(), load_background(), load_variant_sets() and
# compute_enrichment() all come from here. The page's reproducibility section
# calls exactly these four functions.
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "fash_strober_enhancer_comparison", "enrichment_api.R"
))

# ---- Displayed scope --------------------------------------------------------

# Starting, intermediate, and terminal states of the differentiation.
roadmap_enhancers <- c(
  "Roadmap E020 iPS-20b: Enhancer" = "E020 iPSC",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer" = "E013 mesoderm",
  "Roadmap E095 left ventricle: Enhancer" = "E095 left ventricle"
)

expected_sets <- c(
  "current_all", "current_lead", "linear_all", "linear_lead",
  "quadratic_all", "quadratic_lead", "current_only_all", "current_only_lead"
)
expected_variant_counts <- c(
  current_all = 9148L,
  current_lead = 1169L,
  linear_all = 5387L,
  linear_lead = 548L,
  quadratic_all = 6797L,
  quadratic_lead = 690L,
  current_only_all = 7981L,
  current_only_lead = 1029L
)

strategy_levels <- c("All variants", "One lead variant per gene")

# ---- Discovery-set definitions ---------------------------------------------

# The variant sets themselves come from discovery_sets.rds via
# load_variant_sets(); this cache supplies only the labels and the per-set
# variant and gene counts quoted in the text.
cache_path <- file.path(
  workflowr_root, "output", "revision_simulations", "internal",
  "fash_strober_enhancer_comparison_fashr0143", "analysis_cache.rds"
)
if (!file.exists(cache_path)) {
  stop(
    "The retained comparison cache is missing. Run ",
    "code/revision_simulations/internal/fash_strober_enhancer_comparison/",
    "run_fash_strober_enhancer_comparison.R first."
  )
}
cache <- readRDS(cache_path)
configuration <- cache$configuration
set_metadata <- cache$set_metadata
discovery_set_summary <- cache$discovery_set_summary
strict_exclusivity_audit <- cache$strict_exclusivity_audit

variant_sets <- load_variant_sets()

# One consistency check, on the thing the displayed numbers actually depend on:
# the eight variant sets must be the eight published sets.
discovery_order <- match(expected_sets, discovery_set_summary$discovery_set)
if (!identical(configuration$cache_id, "fash_strober_enhancer_comparison_fashr0143") ||
    !identical(set_metadata$discovery_set, expected_sets) ||
    anyNA(discovery_order) ||
    !identical(
      stats::setNames(
        as.integer(discovery_set_summary$unique_variant_count[discovery_order]),
        expected_sets
      ),
      expected_variant_counts
    ) ||
    !identical(
      stats::setNames(lengths(variant_sets[expected_sets]), expected_sets),
      expected_variant_counts
    ) ||
    strict_exclusivity_audit$overlap_with_strober_union[1L] != 0L) {
  stop("The retained discovery sets failed validation.")
}

# ---- The displayed estimates -----------------------------------------------

annotations <- load_annotation_groups("custom", groups = names(roadmap_enhancers))
background <- load_background()
background_variant_count <- nrow(background)

enrichment <- compute_enrichment(
  variant_sets,
  background = background,
  annotations = annotations,
  controls = "background"
)
enrichment <- merge(
  enrichment,
  set_metadata[, c(
    "discovery_set", "method", "section", "exclusivity", "selection_strategy"
  )],
  by.x = "variant_set",
  by.y = "discovery_set",
  all.x = TRUE,
  sort = FALSE
)
if (nrow(enrichment) != length(expected_sets) * length(roadmap_enhancers) ||
    anyNA(enrichment$method) ||
    any(!is.finite(enrichment$fold_enrichment)) ||
    any(!is.finite(enrichment$ci_lower_fold)) ||
    any(!is.finite(enrichment$ci_upper_fold))) {
  stop("The displayed enrichment estimates are incomplete.")
}

# ---- Formatting ------------------------------------------------------------

format_integer <- function(value) {
  format(as.integer(round(value)), big.mark = ",", scientific = FALSE)
}

format_decimal <- function(value, digits = 2L) {
  formatC(as.numeric(value), format = "f", digits = digits)
}

enrichment_lookup <- function(set_name, annotation) {
  row <- enrichment[
    enrichment$variant_set == set_name & enrichment$annotation == annotation,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L) {
    stop("Enrichment lookup failed for ", set_name, " / ", annotation, ".")
  }
  row
}

# "1.38-fold (95% CI 1.12-1.71)"
describe_fold <- function(set_name, annotation, digits = 2L) {
  row <- enrichment_lookup(set_name, annotation)
  paste0(
    format_decimal(row$fold_enrichment, digits), "-fold (95% CI ",
    format_decimal(row$ci_lower_fold, digits), "-",
    format_decimal(row$ci_upper_fold, digits), ")"
  )
}

fold_value <- function(set_name, annotation, digits = 2L) {
  format_decimal(enrichment_lookup(set_name, annotation)$fold_enrichment, digits)
}

set_size <- function(set_name) {
  format_integer(unname(expected_variant_counts[[set_name]]))
}

gene_count <- function(set_name) {
  format_integer(discovery_set_summary$unique_gene_count[
    discovery_set_summary$discovery_set == set_name
  ])
}

e020 <- "Roadmap E020 iPS-20b: Enhancer"
e013 <- "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer"
e095 <- "Roadmap E095 left ventricle: Enhancer"

# ---- Figures ---------------------------------------------------------------

method_colors <- c(
  "Current FASH" = "#0072B2",
  "Strober linear" = "#D55E00",
  "Strober quadratic/nonlinear" = "#8E5BB7"
)
reference_color <- "#9A9A9A"

# One shared grammar for both figures: annotations down the rows in
# differentiation order, all-variants versus lead-per-gene across the columns,
# series in colour, fold enrichment on the x-axis with a reference line at 1.
enhancer_forest <- function(data, series_colors) {
  data$annotation_label <- factor(
    unname(roadmap_enhancers[data$annotation]),
    levels = rev(unname(roadmap_enhancers))
  )
  data$selection_strategy <- factor(
    data$selection_strategy,
    levels = strategy_levels
  )
  data$series <- factor(data$series, levels = names(series_colors))
  dodge <- ggplot2::position_dodge(width = 0.55)
  # `group = series` keeps the points dodged onto the same rows as their own
  # intervals.
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = fold_enrichment,
      y = annotation_label,
      color = series,
      group = series
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      color = "grey45",
      linewidth = 0.45
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lower_fold, xmax = ci_upper_fold),
      orientation = "y",
      width = 0.22,
      linewidth = 0.5,
      position = dodge
    ) +
    ggplot2::geom_point(size = 2.1, position = dodge) +
    ggplot2::facet_wrap(ggplot2::vars(selection_strategy), nrow = 1L) +
    ggplot2::scale_color_manual(values = series_colors, drop = FALSE) +
    ggplot2::labs(
      x = "Fold enrichment relative to all tested variants",
      y = NULL,
      color = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10)
    )
}

# Section 1: FASH against both Strober analyses, no exclusivity filter.
plot_method_comparison <- function() {
  data <- enrichment[enrichment$section == 1L, , drop = FALSE]
  data$series <- data$method
  enhancer_forest(data, method_colors[names(method_colors)])
}

# Section 2: the FASH-only set, with the unfiltered FASH set behind it in grey
# so the effect of the exclusion is visible.
plot_fash_only <- function() {
  exclusive_label <- "FASH only (Strober rsIDs removed)"
  reference_label <- "FASH, all discoveries"
  exclusive <- enrichment[
    enrichment$exclusivity == "Variant-level", , drop = FALSE
  ]
  reference <- enrichment[
    enrichment$section == 1L & enrichment$method == "Current FASH", , drop = FALSE
  ]
  exclusive$series <- exclusive_label
  reference$series <- reference_label
  colors <- stats::setNames(
    c(unname(method_colors[["Current FASH"]]), reference_color),
    c(exclusive_label, reference_label)
  )
  enhancer_forest(rbind(exclusive, reference), colors)
}
