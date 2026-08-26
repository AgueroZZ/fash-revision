# Reporting helpers for the R6 enhancer-enrichment page.
#
# The page is a short, exploratory, orthogonal check: do the lead variants of
# the genes FASH discovers fall inside established Roadmap enhancer annotations
# more often than the tested-variant background does?
#
# Displayed scope, deliberately narrow:
#   * thirteen Roadmap epigenomes in four biological groups (iPSC, germ layer,
#     cardiac, and aorta as a related vascular annotation);
#   * the Strober et al. (2019) four-state enhancer definition;
#   * one lead variant per discovered gene, and nothing else;
#   * fold enrichment against all tested variants, no matched controls;
#   * leave-one-autosome-out jackknife 95% intervals.
#
# The estimates are read from the cache built by
# roadmap_enhancer_exploration/run_roadmap_enhancer_exploration.R rather than
# recomputed here: that script downloads and checksums the thirteen
# segmentations, builds the indicator table, and verifies that its own estimates
# reproduce the previously published three-epigenome R6 numbers exactly, so the
# cache is not a place where the estimator could quietly drift. The panel, the
# state definitions and the estimator itself are all reused from the validated
# modules rather than restated here.
#
# The all-variant view, both enhancer definitions side by side, and the full
# per-cell tables remain available on the internal exploration page; this file
# simply does not surface them.

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

for (package in c("ggplot2", "knitr")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The ", package, " package is required to render this page.")
  }
}

# The epigenome manifest, the two enhancer-state definitions, and the location
# of the validated cache all come from here. Nothing is redefined.
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "roadmap_enhancer_exploration", "roadmap_enhancer_exploration_helpers.R"
))

# load_variant_sets() comes from here, and it is the same API the page's
# reproducibility section calls.
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "fash_strober_enhancer_comparison", "enrichment_api.R"
))

# ---- Displayed scope --------------------------------------------------------

# Strober et al. (2019), get_valid_markers() in
# dynamic_eqtl_calling/perform_tissue_specific_chrom_hmm_enrichment_analysis.py:
# 7_Enh, 6_EnhG, 12_EnhBiv, 11_BivFlnk. Held against the constant the
# exploration module built its annotation table from, so the page cannot show
# estimates from a different state set than the one it documents.
enhancer_definition <- "Strober"
strober_enhancer_states <- c("7_Enh", "6_EnhG", "12_EnhBiv", "11_BivFlnk")

# One lead variant per discovered gene: the three Section 1 discovery sets plus
# the FASH-only set of Section 2. The all-variant sets are deliberately absent.
displayed_sets <- c(
  "current_lead", "linear_lead", "quadratic_lead", "current_only_lead"
)
expected_variant_counts <- c(
  current_lead = 1169L,
  linear_lead = 548L,
  quadratic_lead = 690L,
  current_only_lead = 1029L
)

method_levels <- c(
  "Current FASH", "Strober linear", "Strober quadratic/nonlinear"
)

# Display names for the manifest's biological groups. Roadmap files E065 aorta
# under its `Heart` GROUP, but its `ANATOMY` is `VASCULAR`; it keeps a group of
# its own so that no statement about cardiac epigenomes is silently carried by a
# vascular tissue.
group_display_labels <- c(
  "iPSC" = "iPSC",
  "Germ layer" = "Germ layer",
  "Heart" = "Cardiac",
  "Vascular" = "Vascular"
)

# ---- The verified panel -----------------------------------------------------

panel <- load_epigenome_panel()
expected_epigenomes <- c(
  "E018", "E019", "E020", "E021", "E022",
  "E011", "E012", "E013",
  "E083", "E095", "E104", "E105",
  "E065"
)

# ---- The displayed estimates ------------------------------------------------

cache_path <- file.path(roadmap_output_directory(), "analysis_cache.rds")
if (!file.exists(cache_path)) {
  stop(
    "The Roadmap enhancer cache is missing. Run ",
    "code/revision_simulations/internal/roadmap_enhancer_exploration/",
    "run_roadmap_enhancer_exploration.R first."
  )
}
cache <- readRDS(cache_path)

configuration <- cache$configuration
minimum_overlap <- cache$minimum_overlap
background_variant_count <- configuration$background_variant_count

enrichment <- cache$enrichment
enrichment <- enrichment[
  enrichment$enhancer_definition == enhancer_definition &
    enrichment$variant_set %in% displayed_sets,
  ,
  drop = FALSE
]

variant_sets <- load_variant_sets()

# ---- One consistency gate on everything the page claims ---------------------
#
# Each clause below corresponds to a statement the page makes in prose, so a
# stale cache or an edited manifest fails the render rather than publishing a
# number that no longer matches its own description.
strober_union <- union(variant_sets$linear_all, variant_sets$quadratic_all)
if (
  # The enhancer definition is exactly Strober's four states.
  !setequal(strober_enhancer_states, STROBER_ENHANCER_STATES) ||
    !setequal(strober_enhancer_states,
              ENHANCER_DEFINITIONS[[enhancer_definition]]) ||
    # The panel is the thirteen verified epigenomes, in four groups.
    nrow(panel) != 13L ||
    !identical(panel$epigenome_id, expected_epigenomes) ||
    !setequal(unique(panel$biological_group), names(group_display_labels)) ||
    !identical(
      biological_group_levels(panel),
      intersect(names(group_display_labels), panel$biological_group)
    ) ||
    # Every displayed cell is present, once: thirteen epigenomes x four sets.
    !identical(sort(unique(enrichment$epigenome_id)),
               sort(expected_epigenomes)) ||
    !identical(sort(unique(enrichment$variant_set)), sort(displayed_sets)) ||
    nrow(enrichment) != 13L * length(displayed_sets) ||
    # Lead-per-gene set sizes are the published ones.
    !identical(
      configuration$expected_variant_counts[displayed_sets],
      expected_variant_counts
    ) ||
    !identical(
      stats::setNames(lengths(variant_sets[displayed_sets]), displayed_sets),
      expected_variant_counts
    ) ||
    # The FASH-only definition: rsID-level, applied after the lead variant of
    # each gene has been chosen, so it is a subset of the FASH lead set that
    # shares no rsID with either Strober analysis.
    !all(variant_sets$current_only_lead %in% variant_sets$current_lead) ||
    length(intersect(variant_sets$current_only_lead, strober_union)) != 0L ||
    # A cell may lack an estimate only for the documented reason: fewer than
    # `minimum_overlap` selected variants inside the annotation.
    any(!enrichment$estimable &
          enrichment$selected_overlap >= minimum_overlap) ||
    any(enrichment$estimable &
          (!is.finite(enrichment$fold_enrichment) |
             !is.finite(enrichment$ci_lower_fold) |
             !is.finite(enrichment$ci_upper_fold))) ||
    # The cache reproduced the previously published R6 estimates exactly, which
    # is what licenses reading estimates from it instead of recomputing them.
    any(cache$r6_agreement$n_disagreements != 0L) ||
    max(abs(cache$r6_estimate_reproduction$fold_enrichment_published -
              cache$r6_estimate_reproduction$fold_enrichment_here),
        na.rm = TRUE) > 1e-10
) {
  stop("The displayed Roadmap enhancer estimates failed validation.")
}

enrichment$display_group <- factor(
  unname(group_display_labels[as.character(enrichment$biological_group)]),
  levels = unname(group_display_labels)
)
# Rows read top to bottom in panel order, so `rev()` on the factor levels.
enrichment$display_label <- factor(
  enrichment$display_label,
  levels = rev(panel$display_label)
)

# ---- Formatting -------------------------------------------------------------

format_integer <- function(value) {
  format(as.integer(round(value)), big.mark = ",", scientific = FALSE)
}

format_decimal <- function(value, digits = 2L) {
  formatC(as.numeric(value), format = "f", digits = digits)
}

#' One row of the estimate table, addressed the way the prose addresses it.
enrichment_lookup <- function(set_name, epigenome) {
  row <- enrichment[
    enrichment$variant_set == set_name & enrichment$epigenome_id == epigenome,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L || !row$estimable) {
    stop("Enrichment lookup failed for ", set_name, " / ", epigenome, ".")
  }
  row
}

# "1.48-fold (95% CI 1.15-1.90)"
describe_fold <- function(set_name, epigenome, digits = 2L) {
  row <- enrichment_lookup(set_name, epigenome)
  paste0(
    format_decimal(row$fold_enrichment, digits), "-fold (95% CI ",
    format_decimal(row$ci_lower_fold, digits), "-",
    format_decimal(row$ci_upper_fold, digits), ")"
  )
}

fold_value <- function(set_name, epigenome, digits = 2L) {
  format_decimal(enrichment_lookup(set_name, epigenome)$fold_enrichment, digits)
}

set_size <- function(set_name) {
  format_integer(unname(expected_variant_counts[[set_name]]))
}

#' The estimable rows of one discovery set, in panel order.
set_rows <- function(set_name) {
  rows <- enrichment[
    enrichment$variant_set == set_name & enrichment$estimable, , drop = FALSE
  ]
  rows[order(rows$group_order, rows$epigenome_order), , drop = FALSE]
}

#' How many epigenomes have a jackknife interval entirely above 1, of how many.
#'
#' The page reads the panel by counts and by which epigenomes clear 1, never by
#' picking a single best cell.
count_above_one <- function(set_name) {
  rows <- set_rows(set_name)
  paste0(sum(rows$ci_lower_fold > 1), " of ", nrow(rows))
}

#' Which epigenomes clear 1, named, as "E021, E022, E011".
labels_above_one <- function(set_name, what = c("id", "label")) {
  what <- match.arg(what)
  rows <- set_rows(set_name)
  hits <- rows[rows$ci_lower_fold > 1, , drop = FALSE]
  if (!nrow(hits)) {
    return("none")
  }
  paste(
    if (what == "id") hits$epigenome_id else as.character(hits$display_label),
    collapse = ", "
  )
}

#' Range of fold enrichment over the epigenomes that clear 1, for one set.
range_above_one <- function(set_name, digits = 2L) {
  rows <- set_rows(set_name)
  hits <- rows[rows$ci_lower_fold > 1, , drop = FALSE]
  if (!nrow(hits)) {
    stop("No interval above 1 for ", set_name, ".")
  }
  paste0(
    format_decimal(min(hits$fold_enrichment), digits), "-",
    format_decimal(max(hits$fold_enrichment), digits), "-fold"
  )
}

#' The cells with too few overlapping variants for an estimate, named.
not_estimable_note <- function() {
  rows <- enrichment[!enrichment$estimable, , drop = FALSE]
  if (!nrow(rows)) {
    return("")
  }
  paste(
    paste0(
      rows$epigenome_id, " (",
      unname(set_display_labels[as.character(rows$variant_set)]), ")"
    ),
    collapse = ", "
  )
}

set_display_labels <- c(
  current_lead = "Current FASH",
  linear_lead = "Strober linear",
  quadratic_lead = "Strober quadratic/nonlinear",
  current_only_lead = "FASH only"
)

# ---- Figures ----------------------------------------------------------------

method_colors <- c(
  "Current FASH" = "#0072B2",
  "Strober linear" = "#D55E00",
  "Strober quadratic/nonlinear" = "#8E5BB7"
)
reference_color <- "#9A9A9A"

# One shared grammar for both figures: epigenomes down the rows in panel order,
# grouped into biological blocks; series in colour; fold enrichment on the
# x-axis with a reference line at 1.
#
# `facet_grid` with `space = "free_y"` keeps every row the same height whatever
# the group sizes are, so the five-member iPSC block does not visually dominate
# the single vascular row.
enhancer_forest <- function(data, series_colors) {
  # Cells below the minimum-overlap threshold have no estimate to draw; the
  # figure caption names them rather than leaving an unexplained gap.
  data <- data[data$estimable, , drop = FALSE]
  data$series <- factor(data$series, levels = names(series_colors))
  dodge <- ggplot2::position_dodge(width = 0.6)
  # `group = series` keeps the points dodged onto the same rows as their own
  # intervals.
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = fold_enrichment,
      y = display_label,
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
    ggplot2::facet_grid(
      rows = ggplot2::vars(display_group),
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
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
      panel.spacing.y = grid::unit(0.35, "lines"),
      legend.position = "bottom",
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold"),
      strip.placement = "outside",
      axis.text.y = ggplot2::element_text(size = 10)
    )
}

# Section 1: FASH against both Strober analyses, lead variant per gene.
plot_method_comparison <- function() {
  data <- enrichment[enrichment$section == 1L, , drop = FALSE]
  data$series <- factor(data$method, levels = method_levels)
  enhancer_forest(data, method_colors)
}

# Section 2: the FASH-only lead set, with the unfiltered FASH lead set behind it
# in grey so the effect of the exclusion is visible.
plot_fash_only <- function() {
  exclusive_label <- "FASH only (Strober rsIDs removed)"
  reference_label <- "FASH, all lead variants"
  exclusive <- enrichment[
    enrichment$exclusivity == "Variant-level", , drop = FALSE
  ]
  reference <- enrichment[
    enrichment$section == 1L & enrichment$method == "Current FASH", ,
    drop = FALSE
  ]
  exclusive$series <- exclusive_label
  reference$series <- reference_label
  colors <- stats::setNames(
    c(unname(method_colors[["Current FASH"]]), reference_color),
    c(exclusive_label, reference_label)
  )
  enhancer_forest(rbind(exclusive, reference), colors)
}

# ---- Tables -----------------------------------------------------------------

#' The displayed panel: biological group, Roadmap ID, official Roadmap name.
#'
#' `roadmap_std_name` is Roadmap's own label from `EID_metadata.tab`, carried
#' through the manifest rather than retyped.
panel_table <- function() {
  knitr::kable(
    data.frame(
      Group = unname(group_display_labels[panel$biological_group]),
      `Roadmap ID` = panel$epigenome_id,
      `Roadmap tissue / cell type` = panel$roadmap_std_name,
      check.names = FALSE
    ),
    row.names = FALSE
  )
}
