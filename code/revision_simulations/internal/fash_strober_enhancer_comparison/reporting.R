# Reporting helpers for the R6 enhancer-enrichment page.
#
# Displayed scope:
#   * enhancers = Roadmap 15-state ChromHMM states 6, 7 and 12 (He & Wang 2017);
#   * thirteen epigenomes in four biological groups;
#   * one lead variant per discovered gene;
#   * fold enrichment against the remaining tested variants, no matched controls;
#   * leave-one-autosome-out jackknife, pointwise 95% intervals.
#
# Estimates are read from the cache built by run_r6_roadmap_panel.R, which also
# checks that the three-state indicators sit between the {6,7} and {6,7,11,12}
# indicator sets the internal exploration already validated. The broader
# four-state definition is carried in the same cache and surfaces on the page as
# a single sensitivity sentence, not as a second figure.

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

# ---- The cache ---------------------------------------------------------------

cache_path <- file.path(
  workflowr_root, "output", "revision_simulations", "internal",
  "r6_roadmap_enhancer_panel", "analysis_cache.rds"
)
if (!file.exists(cache_path)) {
  stop(
    "The R6 panel cache is missing. Run ",
    "code/revision_simulations/internal/fash_strober_enhancer_comparison/",
    "run_r6_roadmap_panel.R first."
  )
}
cache <- readRDS(cache_path)

configuration <- cache$configuration
panel <- cache$panel
sensitivity <- cache$sensitivity
minimum_overlap <- configuration$minimum_overlap
background_variant_count <- configuration$background_variant_count
expected_variant_counts <- configuration$expected_variant_counts

primary_states <- c("6_EnhG", "7_Enh", "12_EnhBiv")
displayed_sets <- c(
  "current_lead", "linear_lead", "quadratic_lead", "current_only_lead"
)
method_levels <- c(
  "Current FASH", "Strober linear", "Strober quadratic/nonlinear"
)

# Display names for the manifest's biological groups. Roadmap files E065 aorta
# under its `Heart` GROUP but its `ANATOMY` is `VASCULAR`, so it keeps its own
# group and is never counted as cardiac.
group_display_labels <- c(
  "iPSC" = "iPSC", "Germ layer" = "Germ layer",
  "Heart" = "Cardiac", "Vascular" = "Vascular"
)

set_metadata <- data.frame(
  variant_set = displayed_sets,
  method = c(method_levels, NA_character_),
  section = c(1L, 1L, 1L, 2L),
  stringsAsFactors = FALSE
)

enrichment <- cache$enrichment
enrichment <- enrichment[
  enrichment$enhancer_definition == "EpiCompare", , drop = FALSE
]
enrichment <- merge(
  enrichment, set_metadata, by = "variant_set", all.x = TRUE, sort = FALSE
)

# ---- One consistency gate on what the page claims ----------------------------

if (!identical(configuration$primary_states, primary_states) ||
    !identical(sort(configuration$sensitivity_states),
               sort(c(primary_states, "11_BivFlnk"))) ||
    nrow(panel) != 13L ||
    !setequal(unique(panel$biological_group), names(group_display_labels)) ||
    nrow(enrichment) != 13L * length(displayed_sets) ||
    !setequal(unique(enrichment$epigenome_id), panel$epigenome_id) ||
    !setequal(unique(enrichment$variant_set), displayed_sets) ||
    anyNA(enrichment$section) ||
    # a cell may lack an estimate only for the documented reason
    any(!enrichment$estimable &
          enrichment$selected_overlap >= minimum_overlap) ||
    # the estimator must be the unmatched, set-excluded background
    !all(enrichment$controls == "background") ||
    # every displayed set is a lead-per-gene set of unique rsIDs
    !identical(names(expected_variant_counts), displayed_sets) ||
    any(enrichment$selected_total !=
          unname(expected_variant_counts[enrichment$variant_set]))) {
  stop("The displayed R6 enrichment estimates failed validation.")
}

enrichment$display_group <- factor(
  unname(group_display_labels[as.character(enrichment$biological_group)]),
  levels = unname(group_display_labels)
)
# Rows read top to bottom in panel order, so `rev()` on the factor levels.
enrichment$display_label <- factor(
  enrichment$display_label, levels = rev(panel$display_label)
)

# ---- Formatting -------------------------------------------------------------

format_integer <- function(value) {
  format(as.integer(round(value)), big.mark = ",", scientific = FALSE)
}

format_decimal <- function(value, digits = 2L) {
  formatC(as.numeric(value), format = "f", digits = digits)
}

enrichment_lookup <- function(set_name, epigenome) {
  row <- enrichment[
    enrichment$variant_set == set_name & enrichment$epigenome_id == epigenome, ,
    drop = FALSE
  ]
  if (nrow(row) != 1L || !row$estimable) {
    stop("Enrichment lookup failed for ", set_name, " / ", epigenome, ".")
  }
  row
}

# "1.47-fold (95% CI 1.15-1.89)"
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

set_rows <- function(set_name) {
  rows <- enrichment[
    enrichment$variant_set == set_name & enrichment$estimable, , drop = FALSE
  ]
  rows[order(rows$group_order, rows$epigenome_order), , drop = FALSE]
}

#' How many annotations have a pointwise interval entirely above 1, of how many.
#'
#' The page reads the panel by counts and by which annotations clear 1, never by
#' picking a single best cell, and the surrounding prose says "pointwise" so this
#' is not read as thirteen independent tests.
count_above_one <- function(set_name) {
  rows <- set_rows(set_name)
  paste0(sum(rows$ci_lower_fold > 1), " of the ", nrow(rows))
}

labels_above_one <- function(set_name) {
  rows <- set_rows(set_name)
  hits <- rows$epigenome_id[rows$ci_lower_fold > 1]
  if (!length(hits)) "none" else paste(hits, collapse = ", ")
}

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

#' Set-annotation pairs below the minimum-overlap floor, named for the caption.
not_estimable_note <- function() {
  rows <- enrichment[!enrichment$estimable, , drop = FALSE]
  if (!nrow(rows)) {
    return("none")
  }
  paste(
    paste0(rows$epigenome_id, " (",
           unname(set_display_labels[as.character(rows$variant_set)]), ")"),
    collapse = ", "
  )
}

set_display_labels <- c(
  current_lead = "Current FASH",
  linear_lead = "Strober linear",
  quadratic_lead = "Strober quadratic/nonlinear",
  current_only_lead = "FASH only"
)

#' Agreement between the primary three-state and the broader four-state
#' definition, for the one-sentence sensitivity statement.
sensitivity_summary <- function() {
  paste0(
    "fold enrichments correlated at ",
    format_decimal(sensitivity$correlation, 2), " across the ",
    sensitivity$n_cells, " cells estimable under both"
  )
}

# ---- Figures ----------------------------------------------------------------

method_colors <- c(
  "Current FASH" = "#0072B2",
  "Strober linear" = "#D55E00",
  "Strober quadratic/nonlinear" = "#8E5BB7"
)
reference_color <- "#9A9A9A"

# One shared grammar for both figures: annotations down the rows in panel order,
# grouped into biological blocks; series in colour; fold enrichment on the
# x-axis with a reference line at 1.
enhancer_forest <- function(data, series_colors) {
  data <- data[data$estimable, , drop = FALSE]
  data$series <- factor(data$series, levels = names(series_colors))
  if (anyNA(data$series)) {
    stop("A series label did not map onto the colour scale.")
  }
  dodge <- ggplot2::position_dodge(width = 0.6)
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = fold_enrichment, y = display_label, color = series, group = series
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.45
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lower_fold, xmax = ci_upper_fold),
      orientation = "y", width = 0.22, linewidth = 0.5, position = dodge
    ) +
    ggplot2::geom_point(size = 2.1, position = dodge) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(display_group),
      scales = "free_y", space = "free_y", switch = "y"
    ) +
    ggplot2::scale_color_manual(values = series_colors, drop = FALSE) +
    ggplot2::labs(
      x = "Fold enrichment relative to remaining tested variants",
      y = NULL, color = NULL
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

# Section 1: FASH against both Strober analyses.
plot_method_comparison <- function() {
  data <- enrichment[enrichment$section == 1L, , drop = FALSE]
  data$series <- factor(data$method, levels = method_levels)
  enhancer_forest(data, method_colors)
}

# Section 2: the FASH-only lead set, with all FASH lead variants behind it.
plot_fash_only <- function() {
  labels <- c("FASH only (Strober rsIDs removed)", "FASH, all lead variants")
  exclusive <- enrichment[enrichment$variant_set == "current_only_lead", ,
                          drop = FALSE]
  reference <- enrichment[enrichment$variant_set == "current_lead", ,
                          drop = FALSE]
  exclusive$series <- labels[[1L]]
  reference$series <- labels[[2L]]
  enhancer_forest(
    rbind(exclusive, reference),
    stats::setNames(
      c(unname(method_colors[["Current FASH"]]), reference_color), labels
    )
  )
}

# ---- Table ------------------------------------------------------------------

#' The displayed panel: group, Roadmap ID, official Roadmap name.
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
