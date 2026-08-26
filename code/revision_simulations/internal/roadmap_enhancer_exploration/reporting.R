#!/usr/bin/env Rscript

# Reporting helpers for the internal Roadmap enhancer exploration page.
#
# This page is internal and exploratory. It is deliberately not optimised for
# brevity: it shows all thirteen epigenomes, all eight discovery sets, and both
# enhancer-state definitions, so that we can see the whole pattern before
# deciding whether the reviewer-facing R6 page should stay narrow.
#
# The estimates are read from the cache built by
# run_roadmap_enhancer_exploration.R rather than recomputed at render time:
# 8 sets x 26 annotation columns over 745,867 background variants is minutes of
# work, not seconds. The run script verifies that its own R6-definition
# estimates reproduce the published R6 numbers exactly, so the cache is not a
# place where the estimator could quietly drift.

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

source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "roadmap_enhancer_exploration", "roadmap_enhancer_exploration_helpers.R"
))

cache_path <- file.path(roadmap_output_directory(), "analysis_cache.rds")
if (!file.exists(cache_path)) {
  stop(
    "The exploration cache is missing. Run ",
    "code/revision_simulations/internal/roadmap_enhancer_exploration/",
    "run_roadmap_enhancer_exploration.R first."
  )
}
cache <- readRDS(cache_path)

configuration <- cache$configuration
panel <- cache$panel
set_metadata <- cache$set_metadata
annotation_provenance <- cache$annotation_provenance
state_summary <- cache$state_summary
definition_coverage <- cache$definition_coverage
r6_agreement <- cache$r6_agreement
r6_estimate_reproduction <- cache$r6_estimate_reproduction
minimum_overlap <- cache$minimum_overlap
not_estimable <- cache$not_estimable
enrichment <- cache$enrichment
summary_table <- cache$summary_table

background_variant_count <- configuration$background_variant_count

# ---- One consistency gate on what the page displays --------------------------

group_levels <- biological_group_levels(panel)
definition_levels <- c("R6", "Strober")
strategy_levels <- c("All variants", "One lead variant per gene")
method_levels <- c(
  "Current FASH", "Strober linear", "Strober quadratic/nonlinear"
)

if (nrow(panel) != 13L ||
    !identical(sort(unique(enrichment$epigenome_id)), sort(panel$epigenome_id)) ||
    !identical(sort(unique(enrichment$enhancer_definition)),
               sort(definition_levels)) ||
    !identical(sort(unique(enrichment$variant_set)),
               sort(configuration$expected_sets)) ||
    nrow(enrichment) != 13L * 2L * length(configuration$expected_sets) ||
    !all(c("estimable") %in% names(enrichment)) ||
    # A cell is allowed to be missing an estimate only for the documented
    # reason: fewer than `minimum_overlap` selected variants inside the
    # annotation. Anything else non-finite is a bug, not a data limit.
    any(!enrichment$estimable &
          enrichment$selected_overlap >= minimum_overlap) ||
    any(enrichment$estimable &
          (!is.finite(enrichment$fold_enrichment) |
             !is.finite(enrichment$ci_lower_fold) |
             !is.finite(enrichment$ci_upper_fold))) ||
    nrow(not_estimable) != sum(!enrichment$estimable) ||
    any(r6_agreement$n_disagreements != 0L) ||
    max(abs(r6_estimate_reproduction$fold_enrichment_published -
              r6_estimate_reproduction$fold_enrichment_here),
        na.rm = TRUE) > 1e-10) {
  stop("The exploration cache failed validation.")
}

enrichment$biological_group <- factor(
  enrichment$biological_group, levels = group_levels
)
enrichment$enhancer_definition <- factor(
  enrichment$enhancer_definition, levels = definition_levels
)
enrichment$selection_strategy <- factor(
  enrichment$selection_strategy, levels = strategy_levels
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
lookup <- function(set_name, epigenome, definition = "R6") {
  row <- enrichment[
    enrichment$variant_set == set_name &
      enrichment$epigenome_id == epigenome &
      enrichment$enhancer_definition == definition,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L) {
    stop(
      "Lookup failed for ", set_name, " / ", epigenome, " / ", definition, "."
    )
  }
  row
}

# "1.38-fold (95% CI 1.12-1.71)"
describe_fold <- function(set_name, epigenome, definition = "R6", digits = 2L) {
  row <- lookup(set_name, epigenome, definition)
  paste0(
    format_decimal(row$fold_enrichment, digits), "-fold (95% CI ",
    format_decimal(row$ci_lower_fold, digits), "-",
    format_decimal(row$ci_upper_fold, digits), ")"
  )
}

fold_value <- function(set_name, epigenome, definition = "R6", digits = 2L) {
  format_decimal(lookup(set_name, epigenome, definition)$fold_enrichment, digits)
}

set_size <- function(set_name) {
  format_integer(unname(configuration$expected_variant_counts[[set_name]]))
}

#' Range of fold enrichment over a group of epigenomes, for one set.
#'
#' The point of this page is whether a group behaves coherently, so the prose
#' quotes ranges over groups far more often than single estimates.
group_fold_range <- function(set_name, group, definition = "R6", digits = 2L) {
  rows <- enrichment[
    enrichment$variant_set == set_name &
      enrichment$biological_group == group &
      enrichment$enhancer_definition == definition,
    ,
    drop = FALSE
  ]
  rows <- rows[rows$estimable, , drop = FALSE]
  if (!nrow(rows)) {
    stop("No estimates for group ", group, ".")
  }
  paste0(
    format_decimal(min(rows$fold_enrichment), digits), "-",
    format_decimal(max(rows$fold_enrichment), digits), "-fold"
  )
}

#' How many epigenomes in a group have a jackknife interval entirely above 1.
group_above_one <- function(set_name, group, definition = "R6") {
  rows <- enrichment[
    enrichment$variant_set == set_name &
      enrichment$biological_group == group &
      enrichment$enhancer_definition == definition &
      enrichment$estimable,
    ,
    drop = FALSE
  ]
  paste0(sum(rows$ci_lower_fold > 1), " of ", nrow(rows))
}

#' Epigenomes in a group whose jackknife interval lies entirely above 1, named.
#'
#' The page leans on this rather than on a single "the group is enriched"
#' statement, because in several groups only some members are.
group_above_one_labels <- function(set_name, group, definition = "R6") {
  rows <- enrichment[
    enrichment$variant_set == set_name &
      enrichment$biological_group == group &
      enrichment$enhancer_definition == definition &
      enrichment$estimable,
    ,
    drop = FALSE
  ]
  rows <- rows[order(rows$epigenome_order), , drop = FALSE]
  hits <- rows$epigenome_id[rows$ci_lower_fold > 1]
  if (!length(hits)) "none" else paste(hits, collapse = ", ")
}

#' Fold enrichment for every member of a group, as "E018 1.04, E019 1.07, ...".
group_fold_list <- function(set_name, group, definition = "R6", digits = 2L) {
  rows <- enrichment[
    enrichment$variant_set == set_name &
      enrichment$biological_group == group &
      enrichment$enhancer_definition == definition,
    ,
    drop = FALSE
  ]
  rows <- rows[order(rows$epigenome_order), , drop = FALSE]
  paste(
    paste0(
      rows$epigenome_id, " ",
      ifelse(
        rows$estimable,
        format_decimal(rows$fold_enrichment, digits),
        "not estimable"
      )
    ),
    collapse = ", "
  )
}

# ---- Definition sensitivity, summarised --------------------------------------

#' How far the two definitions move any single estimate.
#'
#' Returns the Pearson correlation of fold enrichment across the whole design,
#' the largest absolute shift, and how many cells change their qualitative
#' reading (point estimate crossing 1, or interval crossing 1).
definition_agreement <- function() {
  wide <- merge(
    enrichment[
      enrichment$enhancer_definition == "R6",
      c("variant_set", "epigenome_id", "fold_enrichment", "ci_lower_fold",
        "estimable")
    ],
    enrichment[
      enrichment$enhancer_definition == "Strober",
      c("variant_set", "epigenome_id", "fold_enrichment", "ci_lower_fold",
        "estimable")
    ],
    by = c("variant_set", "epigenome_id"),
    suffixes = c("_r6", "_strober")
  )
  wide <- wide[wide$estimable_r6 & wide$estimable_strober, , drop = FALSE]
  list(
    n_cells = nrow(wide),
    correlation = stats::cor(wide$fold_enrichment_r6,
                             wide$fold_enrichment_strober),
    max_absolute_shift = max(abs(
      wide$fold_enrichment_r6 - wide$fold_enrichment_strober
    )),
    point_crossings = sum(xor(wide$fold_enrichment_r6 > 1,
                              wide$fold_enrichment_strober > 1)),
    interval_crossings = sum(xor(wide$ci_lower_fold_r6 > 1,
                                 wide$ci_lower_fold_strober > 1))
  )
}

# ---- Figures ----------------------------------------------------------------

method_colors <- c(
  "Current FASH" = "#0072B2",
  "Strober linear" = "#D55E00",
  "Strober quadratic/nonlinear" = "#8E5BB7"
)
fash_only_colors <- c(
  "FASH only (Strober rsIDs removed)" = "#0072B2",
  "FASH, all discoveries" = "#9A9A9A"
)
definition_colors <- c(
  "R6 (EnhG+Enh)" = "#0072B2",
  "Strober (EnhG+Enh+EnhBiv+BivFlnk)" = "#009E73"
)

# The three epigenomes R6 shows are drawn as filled circles, the ten new ones as
# open circles, so the reader can see at a glance whether R6's choice sits inside
# or outside its group's pattern.
r6_shapes <- c("Shown in R6" = 16L, "Added here" = 21L)

r6_membership <- function(data) {
  factor(
    ifelse(data$in_r6_page, "Shown in R6", "Added here"),
    levels = names(r6_shapes)
  )
}

#' The shared forest grammar: epigenomes down the rows, grouped into biological
#' blocks; all-variants and lead-per-gene across the columns; series in colour.
#'
#' `facet_grid` with `space = "free_y"` keeps every row the same height whatever
#' the group sizes are, so the iPSC block does not visually dominate the germ
#' layer block.
enhancer_forest <- function(data,
                            series_colors,
                            columns = "selection_strategy") {
  # Cells below the minimum-overlap threshold have no estimate to draw. They are
  # listed explicitly in the not-estimable table rather than shown as a gap the
  # reader has to interpret.
  data <- data[data$estimable, , drop = FALSE]
  data$series <- factor(data$series, levels = names(series_colors))
  data$r6_status <- r6_membership(data)
  dodge <- ggplot2::position_dodge(width = 0.62)
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
      xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.45
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lower_fold, xmax = ci_upper_fold),
      orientation = "y", width = 0.24, linewidth = 0.45, position = dodge
    ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = r6_status),
      size = 1.9, fill = "white", stroke = 0.7, position = dodge
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(biological_group),
      cols = ggplot2::vars(.data[[columns]]),
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    ggplot2::scale_color_manual(values = series_colors, drop = FALSE) +
    ggplot2::scale_shape_manual(values = r6_shapes, drop = FALSE) +
    ggplot2::labs(
      x = "Fold enrichment relative to all tested variants",
      y = NULL, color = NULL, shape = NULL
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1L, nrow = 2L),
      shape = ggplot2::guide_legend(order = 2L, nrow = 2L)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.spacing.y = grid::unit(0.35, "lines"),
      legend.position = "bottom",
      legend.box = "horizontal",
      strip.text.x = ggplot2::element_text(face = "bold"),
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold"),
      strip.placement = "outside",
      axis.text.y = ggplot2::element_text(size = 9)
    )
}

#' Figure A: FASH against both Strober analyses, across the whole panel.
plot_method_comparison <- function(definition = "R6") {
  data <- enrichment[
    enrichment$section == 1L & enrichment$enhancer_definition == definition,
    ,
    drop = FALSE
  ]
  data$series <- factor(data$method, levels = method_levels)
  enhancer_forest(data, method_colors)
}

#' Figure B: the FASH-only sets, with all FASH discoveries behind them in grey.
plot_fash_only <- function(definition = "R6") {
  exclusive <- enrichment[
    enrichment$exclusivity == "Variant-level" &
      enrichment$enhancer_definition == definition,
    ,
    drop = FALSE
  ]
  reference <- enrichment[
    enrichment$section == 1L & enrichment$method == "Current FASH" &
      enrichment$enhancer_definition == definition,
    ,
    drop = FALSE
  ]
  exclusive$series <- names(fash_only_colors)[1L]
  reference$series <- names(fash_only_colors)[2L]
  enhancer_forest(rbind(exclusive, reference), fash_only_colors)
}

#' Figure C: the two enhancer-state definitions side by side.
#'
#' Same rows, same estimator, only the state set differs. Kept as its own figure
#' so nothing about the R6 definition is silently replaced.
plot_definition_sensitivity <- function(
  sets = c("current_all", "current_lead", "current_only_lead")
) {
  data <- enrichment[enrichment$variant_set %in% sets, , drop = FALSE]
  data$series <- factor(
    ifelse(
      data$enhancer_definition == "R6",
      names(definition_colors)[1L],
      names(definition_colors)[2L]
    ),
    levels = names(definition_colors)
  )
  data$set_label <- factor(
    set_display_labels[as.character(data$variant_set)],
    levels = unname(set_display_labels[sets])
  )
  enhancer_forest(data, definition_colors, columns = "set_label")
}

# Human-readable names for the eight sets, used by the heatmap and the
# sensitivity figure.
set_display_labels <- c(
  current_all = "FASH, all",
  current_lead = "FASH, lead",
  linear_all = "Strober linear, all",
  linear_lead = "Strober linear, lead",
  quadratic_all = "Strober quadratic, all",
  quadratic_lead = "Strober quadratic, lead",
  current_only_all = "FASH only, all",
  current_only_lead = "FASH only, lead"
)

#' A fold-enrichment heatmap over all eight sets and thirteen epigenomes.
#'
#' Complements the forest plots: it drops the intervals but fits the entire
#' design on one screen, which is the only way to see the block structure.
plot_fold_heatmap <- function(definition = "R6") {
  data <- enrichment[
    enrichment$enhancer_definition == definition, , drop = FALSE
  ]
  data$set_label <- factor(
    set_display_labels[as.character(data$variant_set)],
    levels = unname(set_display_labels[configuration$expected_sets])
  )
  data <- data[data$estimable, , drop = FALSE]
  # Diverging around 1 on the log scale, so "no enrichment" is the neutral
  # colour and depletion and enrichment are visually symmetric.
  limit <- max(abs(log2(data$fold_enrichment)))
  ggplot2::ggplot(
    data,
    ggplot2::aes(x = set_label, y = display_label, fill = log2(fold_enrichment))
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = format_decimal(fold_enrichment, 2L),
        fontface = ifelse(ci_lower_fold > 1, "bold", "plain")
      ),
      size = 2.7,
      color = "grey15"
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(biological_group),
      scales = "free_y", space = "free_y", switch = "y"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#3B6FB6", mid = "grey96", high = "#C0392B",
      midpoint = 0, limits = c(-limit, limit),
      name = expression(log[2] * " fold")
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.spacing.y = grid::unit(0.3, "lines"),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 9),
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold"),
      strip.placement = "outside",
      legend.position = "right"
    )
}

# ---- Tables -----------------------------------------------------------------

#' The verified epigenome panel.
panel_table <- function() {
  knitr::kable(
    data.frame(
      Group = panel$biological_group,
      `Epigenome ID` = panel$epigenome_id,
      `Roadmap STD_NAME` = panel$roadmap_std_name,
      Mnemonic = panel$roadmap_mnemonic,
      `Roadmap GROUP` = panel$roadmap_group,
      Anatomy = panel$roadmap_anatomy,
      `In R6` = ifelse(panel$in_r6_page, "yes", ""),
      check.names = FALSE
    ),
    row.names = FALSE
  )
}

#' Genome coverage of the two definitions, per epigenome.
definition_coverage_table <- function() {
  rows <- definition_coverage[match(panel$epigenome_id,
                                    definition_coverage$epigenome_id), ,
                              drop = FALSE]
  knitr::kable(
    data.frame(
      `Epigenome` = panel$display_label,
      `R6 states (Mb)` = round(rows$r6_covered_bases / 1e6, 1),
      `Strober states (Mb)` = round(rows$strober_covered_bases / 1e6, 1),
      `Ratio` = round(rows$strober_over_r6, 2),
      check.names = FALSE
    ),
    row.names = FALSE
  )
}

#' Where the two definitions differ, state by state, on one epigenome.
state_table <- function(epigenome = "E095") {
  rows <- state_summary[state_summary$epigenome_id == epigenome, , drop = FALSE]
  rows <- rows[order(as.integer(sub("_.*$", "", rows$state))), , drop = FALSE]
  knitr::kable(
    data.frame(
      State = rows$state,
      `R6` = ifelse(rows$in_r6_definition, "yes", ""),
      `Strober` = ifelse(rows$in_strober_definition, "yes", ""),
      `Intervals` = format_integer(rows$interval_count),
      `Autosomal Mb` = round(rows$autosomal_covered_bases / 1e6, 1),
      check.names = FALSE
    ),
    row.names = FALSE
  )
}

#' The full estimate table: every epigenome x discovery set, one definition.
full_summary_table <- function(definition = "R6") {
  rows <- enrichment[
    enrichment$enhancer_definition == definition, , drop = FALSE
  ]
  rows <- rows[order(
    rows$group_order, rows$epigenome_order,
    match(rows$variant_set, configuration$expected_sets)
  ), , drop = FALSE]
  knitr::kable(
    data.frame(
      ID = rows$epigenome_id,
      `Roadmap name` = rows$roadmap_std_name,
      Group = as.character(rows$biological_group),
      `Discovery set` = unname(
        set_display_labels[as.character(rows$variant_set)]
      ),
      `N variants` = format_integer(rows$selected_total),
      `Overlap prop.` = format_decimal(
        rows$selected_overlap / rows$selected_total, 4L
      ),
      `Background prop.` = format_decimal(
        rows$control_overlap / rows$control_total, 4L
      ),
      Fold = ifelse(rows$estimable,
                    format_decimal(rows$fold_enrichment, 3L), "-"),
      `CI lower` = ifelse(rows$estimable,
                          format_decimal(rows$ci_lower_fold, 3L), "-"),
      `CI upper` = ifelse(rows$estimable,
                          format_decimal(rows$ci_upper_fold, 3L), "-"),
      check.names = FALSE
    ),
    row.names = FALSE
  )
}

#' The cells with too few overlapping variants to estimate.
not_estimable_table <- function() {
  if (!nrow(not_estimable)) {
    return("Every cell in the design cleared the minimum-overlap threshold.")
  }
  knitr::kable(
    data.frame(
      Definition = not_estimable$enhancer_definition,
      `Discovery set` = unname(
        set_display_labels[as.character(not_estimable$variant_set)]
      ),
      Epigenome = not_estimable$epigenome_id,
      `Roadmap name` = not_estimable$roadmap_std_name,
      `N variants` = format_integer(not_estimable$selected_total),
      `Overlapping` = format_integer(not_estimable$selected_overlap),
      check.names = FALSE
    ),
    row.names = FALSE
  )
}

#' The annotation files actually used.
provenance_table <- function() {
  rows <- annotation_provenance[
    match(panel$epigenome_id, annotation_provenance$epigenome_id), ,
    drop = FALSE
  ]
  knitr::kable(
    data.frame(
      ID = rows$epigenome_id,
      File = basename(rows$local_path),
      `Bytes` = format_integer(rows$byte_size),
      MD5 = rows$md5,
      check.names = FALSE
    ),
    row.names = FALSE
  )
}
