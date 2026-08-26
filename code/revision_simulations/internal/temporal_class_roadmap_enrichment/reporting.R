#!/usr/bin/env Rscript

# Reporting layer for the temporal-class Roadmap enrichment page.
#
# Internal and exploratory, and deliberately not compressed: four classes, two
# selection strategies, thirteen epigenomes, two comparators and two sensitivity
# axes all get shown.
#
# Estimates are read from the cache built by
# run_temporal_class_roadmap_enrichment.R. That script checks that pooling the
# three timing classes reproduces the thirteen-epigenome page's `current_all`
# estimates exactly, so the numbers here cannot drift away from the sister page.

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
  "temporal_class_roadmap_enrichment", "temporal_class_roadmap_helpers.R"
))

cache_path <- file.path(tcr_output_directory(), "analysis_cache.rds")
if (!file.exists(cache_path)) {
  stop(
    "The cache is missing. Run ",
    "code/revision_simulations/internal/temporal_class_roadmap_enrichment/",
    "run_temporal_class_roadmap_enrichment.R first."
  )
}
cache <- readRDS(cache_path)

configuration <- cache$configuration
panel <- cache$panel
class_summary <- cache$class_summary
switch_crosstab <- cache$switch_crosstab
set_sizes <- cache$set_sizes
lead_rule_agreement <- cache$lead_rule_agreement
estimability <- cache$estimability_forecast
pooled_crosscheck <- cache$pooled_timing_crosscheck
middle_definition <- cache$middle_definition
not_estimable <- cache$not_estimable
enrichment <- cache$enrichment

minimum_overlap <- configuration$minimum_overlap
background_variant_count <- configuration$background_variant_count
discovered_pairs <- configuration$discovered_pairs
discovered_variants <- configuration$discovered_variants

# ---- Consistency gate -------------------------------------------------------

group_levels <- biological_group_levels(panel)
class_levels <- c(TCR_CLASSES, "no_switch")
strategy_levels <- unname(SELECTION_STRATEGIES)
definition_levels <- c("R6", "Strober")

if (nrow(panel) != 13L ||
    !identical(sort(unique(enrichment$epigenome_id)), sort(panel$epigenome_id)) ||
    !identical(sort(unique(enrichment$enhancer_definition)),
               sort(definition_levels)) ||
    !all(TCR_CLASSES %in% enrichment$class) ||
    nrow(enrichment) != 676L ||
    # A missing estimate is allowed only for the documented reason: fewer than
    # `minimum_overlap` selected variants inside the annotation.
    any(!enrichment$estimable &
          enrichment$selected_overlap >= minimum_overlap) ||
    any(enrichment$estimable &
          (!is.finite(enrichment$fold_enrichment) |
             !is.finite(enrichment$ci_lower_fold) |
             !is.finite(enrichment$ci_upper_fold))) ||
    nrow(not_estimable) != sum(!enrichment$estimable) ||
    # The classification must be the open (3, 12) middle window, not the stale
    # closed [4, 11] one.
    !identical(middle_definition$window, "open (3, 12)") ||
    abs(middle_definition$mean_residual_open) > MIDDLE_WINDOW_TOLERANCE ||
    # Pooling the three timing classes must reproduce the sister page exactly.
    max(abs(pooled_crosscheck$fold_enrichment_pooled -
              pooled_crosscheck$fold_enrichment_published), na.rm = TRUE) >
      1e-10) {
  stop("The temporal-class cache failed validation.")
}

enrichment$biological_group <- factor(
  enrichment$biological_group, levels = group_levels
)
enrichment$class <- factor(enrichment$class, levels = class_levels)
CLASS_DISPLAY_LABELS <- c(TCR_CLASS_SHORT_LABELS, no_switch = "Non-switch")
if (!all(class_levels %in% names(CLASS_DISPLAY_LABELS))) {
  stop("A class has no display label.")
}
enrichment$class_label <- factor(
  unname(CLASS_DISPLAY_LABELS[as.character(enrichment$class)]),
  levels = unname(CLASS_DISPLAY_LABELS)
)
if (anyNA(enrichment$class_label)) {
  stop("A class could not be mapped to a display label.")
}
enrichment$strategy_label <- factor(
  unname(SELECTION_STRATEGIES[enrichment$strategy]), levels = strategy_levels
)
enrichment$enhancer_definition <- factor(
  enrichment$enhancer_definition, levels = definition_levels
)
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

#' Rows matching a query, with every filter spelled out.
#'
#' The page's prose addresses cells by (class, strategy, epigenome, comparator),
#' so the accessor takes exactly those and fails loudly on a miss rather than
#' silently returning the wrong row.
select_rows <- function(class_name,
                        strategy = "all",
                        epigenome = NULL,
                        group = NULL,
                        comparator = "Tested background",
                        definition = "R6",
                        lead_rule = "lfsr") {
  rows <- enrichment[
    as.character(enrichment$class) == class_name &
      enrichment$strategy == strategy &
      enrichment$comparator == comparator &
      as.character(enrichment$enhancer_definition) == definition &
      enrichment$lead_rule == lead_rule,
    ,
    drop = FALSE
  ]
  if (!is.null(epigenome)) {
    rows <- rows[rows$epigenome_id %in% epigenome, , drop = FALSE]
  }
  if (!is.null(group)) {
    rows <- rows[as.character(rows$biological_group) %in% group, , drop = FALSE]
  }
  if (!nrow(rows)) {
    stop("No rows for ", class_name, " / ", strategy, " / ",
         paste(c(epigenome, group), collapse = ","), ".")
  }
  rows[order(rows$group_order, rows$epigenome_order), , drop = FALSE]
}

# "0.43-fold (95% CI 0.25-0.74)"
describe_fold <- function(class_name, epigenome, ..., digits = 2L) {
  row <- select_rows(class_name, epigenome = epigenome, ...)
  if (nrow(row) != 1L) {
    stop("describe_fold matched ", nrow(row), " rows.")
  }
  if (!row$estimable) {
    return("not estimable")
  }
  paste0(
    format_decimal(row$fold_enrichment, digits), "-fold (95% CI ",
    format_decimal(row$ci_lower_fold, digits), "-",
    format_decimal(row$ci_upper_fold, digits), ")"
  )
}

fold_value <- function(class_name, epigenome, ..., digits = 2L) {
  row <- select_rows(class_name, epigenome = epigenome, ...)
  if (nrow(row) != 1L) {
    stop("fold_value matched ", nrow(row), " rows.")
  }
  if (!row$estimable) "not estimable" else {
    format_decimal(row$fold_enrichment, digits)
  }
}

#' Range of fold enrichment over the estimable cells of a query.
fold_range <- function(class_name, ..., digits = 2L) {
  rows <- select_rows(class_name, ...)
  rows <- rows[rows$estimable, , drop = FALSE]
  if (!nrow(rows)) {
    return("not estimable")
  }
  paste0(
    format_decimal(min(rows$fold_enrichment), digits), "-",
    format_decimal(max(rows$fold_enrichment), digits), "-fold"
  )
}

#' "4 of 4" — estimable cells whose interval lies entirely above 1.
count_above_one <- function(class_name, ...) {
  rows <- select_rows(class_name, ...)
  rows <- rows[rows$estimable, , drop = FALSE]
  paste0(sum(rows$ci_lower_fold > 1), " of ", nrow(rows))
}

#' "9 of 12" — estimable cells whose interval lies entirely below 1.
count_below_one <- function(class_name, ...) {
  rows <- select_rows(class_name, ...)
  rows <- rows[rows$estimable, , drop = FALSE]
  paste0(sum(rows$ci_upper_fold < 1), " of ", nrow(rows))
}

#' "3 of 13" — how many cells could be estimated at all.
count_estimable <- function(class_name, ...) {
  rows <- select_rows(class_name, ...)
  paste0(sum(rows$estimable), " of ", nrow(rows))
}

#' "E018 0.47, E019 0.46, ..." over the members of a group.
fold_list <- function(class_name, ..., digits = 2L) {
  rows <- select_rows(class_name, ...)
  paste(
    paste0(
      rows$epigenome_id, " ",
      ifelse(rows$estimable, format_decimal(rows$fold_enrichment, digits), "n/e")
    ),
    collapse = ", "
  )
}

class_size <- function(class_name, strategy = "all", lead_rule = "lfsr") {
  rows <- set_sizes[
    set_sizes$lead_rule == lead_rule &
      set_sizes$variant_set == paste0(class_name, "_", strategy),
  ]
  if (nrow(rows) != 1L) {
    stop("No set size for ", class_name, "_", strategy, ".")
  }
  format_integer(rows$n_variants)
}

class_pairs <- function(class_name) {
  format_integer(class_summary$pairs[class_summary$class == class_name])
}

class_genes <- function(class_name) {
  format_integer(class_summary$genes[class_summary$class == class_name])
}

# ---- Sensitivity summaries --------------------------------------------------

#' Agreement between the two enhancer-state definitions over the whole design.
definition_agreement <- function() {
  wide <- merge(
    enrichment[
      enrichment$enhancer_definition == "R6",
      c("run", "variant_set", "epigenome_id", "fold_enrichment",
        "ci_lower_fold", "estimable")
    ],
    enrichment[
      enrichment$enhancer_definition == "Strober",
      c("run", "variant_set", "epigenome_id", "fold_enrichment",
        "ci_lower_fold", "estimable")
    ],
    by = c("run", "variant_set", "epigenome_id"),
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

#' Agreement between the two lead-selection rules, on the lead sets only.
lead_rule_sensitivity <- function(definition = "R6") {
  rows <- enrichment[
    enrichment$comparator == "Tested background" &
      as.character(enrichment$enhancer_definition) == definition &
      enrichment$strategy == "lead",
    ,
    drop = FALSE
  ]
  wide <- merge(
    rows[rows$lead_rule == "lfsr",
         c("variant_set", "epigenome_id", "fold_enrichment", "ci_lower_fold",
           "estimable")],
    rows[rows$lead_rule == "lfdr",
         c("variant_set", "epigenome_id", "fold_enrichment", "ci_lower_fold",
           "estimable")],
    by = c("variant_set", "epigenome_id"),
    suffixes = c("_lfsr", "_lfdr")
  )
  both <- wide[wide$estimable_lfsr & wide$estimable_lfdr, , drop = FALSE]
  list(
    n_cells = nrow(both),
    n_total = nrow(wide),
    correlation = stats::cor(both$fold_enrichment_lfsr,
                             both$fold_enrichment_lfdr),
    max_absolute_shift = max(abs(
      both$fold_enrichment_lfsr - both$fold_enrichment_lfdr
    )),
    point_crossings = sum(xor(both$fold_enrichment_lfsr > 1,
                              both$fold_enrichment_lfdr > 1)),
    interval_crossings = sum(xor(both$ci_lower_fold_lfsr > 1,
                                 both$ci_lower_fold_lfdr > 1)),
    min_jaccard = min(lead_rule_agreement$jaccard),
    max_jaccard = max(lead_rule_agreement$jaccard)
  )
}

# ---- Figures ----------------------------------------------------------------

class_colors <- c(
  Early = "#E69F00",
  Middle = "#009E73",
  Late = "#0072B2",
  Switch = "#CC79A7",
  `Non-switch` = "#9A9A9A"
)

definition_colors <- c(
  "R6 (EnhG+Enh)" = "#0072B2",
  "Strober (EnhG+Enh+EnhBiv+BivFlnk)" = "#009E73"
)

lead_rule_colors <- c(
  "Lead by class lfsr" = "#0072B2",
  "Lead by overall lfdr" = "#D55E00"
)

r6_shapes <- c("Shown on the 3-epigenome page" = 16L, "Added here" = 21L)

#' Epigenomes down the rows in biological blocks, strategy across the columns,
#' series in colour. Same grammar as the thirteen-epigenome page.
class_forest <- function(data, series_colors, columns = "strategy_label") {
  data <- data[data$estimable, , drop = FALSE]
  data$series <- factor(data$series, levels = names(series_colors))
  data$r6_status <- factor(
    ifelse(data$in_r6_page, names(r6_shapes)[1L], names(r6_shapes)[2L]),
    levels = names(r6_shapes)
  )
  dodge <- ggplot2::position_dodge(width = 0.68)
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
      orientation = "y", width = 0.26, linewidth = 0.45, position = dodge
    ) +
    ggplot2::geom_point(
      ggplot2::aes(shape = r6_status),
      size = 1.9, fill = "white", stroke = 0.7, position = dodge
    ) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(biological_group),
      cols = ggplot2::vars(.data[[columns]]),
      scales = "free_y", space = "free_y", switch = "y"
    ) +
    ggplot2::scale_color_manual(values = series_colors, drop = FALSE) +
    ggplot2::scale_shape_manual(values = r6_shapes, drop = FALSE) +
    ggplot2::scale_x_continuous(trans = "log2") +
    ggplot2::labs(
      x = "Fold enrichment (log2 scale)", y = NULL, color = NULL, shape = NULL
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
      strip.text.x = ggplot2::element_text(face = "bold"),
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold"),
      strip.placement = "outside",
      axis.text.y = ggplot2::element_text(size = 9)
    )
}

#' Figure A: the four classes against the tested background.
plot_classes_vs_background <- function(definition = "R6", lead_rule = "lfsr") {
  data <- enrichment[
    enrichment$comparator == "Tested background" &
      as.character(enrichment$enhancer_definition) == definition &
      enrichment$lead_rule == lead_rule &
      as.character(enrichment$class) %in% TCR_CLASSES,
    ,
    drop = FALSE
  ]
  data$series <- as.character(data$class_label)
  class_forest(data, class_colors[unname(TCR_CLASS_SHORT_LABELS)])
}

#' Figure B: the four classes against the other FASH discoveries.
plot_classes_vs_discoveries <- function(definition = "R6") {
  data <- enrichment[
    enrichment$comparator == "Other FASH discoveries" &
      as.character(enrichment$enhancer_definition) == definition &
      as.character(enrichment$class) %in% TCR_CLASSES,
    ,
    drop = FALSE
  ]
  data$series <- as.character(data$class_label)
  class_forest(data, class_colors[unname(TCR_CLASS_SHORT_LABELS)])
}

#' Figure B2: switch against its own complement, within the discoveries.
plot_switch_contrast <- function(definition = "R6") {
  data <- enrichment[
    enrichment$comparator == "Other FASH discoveries" &
      as.character(enrichment$enhancer_definition) == definition &
      as.character(enrichment$class) %in% c("switch", "no_switch"),
    ,
    drop = FALSE
  ]
  data$series <- as.character(data$class_label)
  class_forest(data, class_colors[c("Switch", "Non-switch")])
}

#' Figure C: the two lead-selection rules side by side.
plot_lead_rule_sensitivity <- function(definition = "R6") {
  data <- enrichment[
    enrichment$comparator == "Tested background" &
      as.character(enrichment$enhancer_definition) == definition &
      enrichment$strategy == "lead" &
      as.character(enrichment$class) %in% TCR_CLASSES,
    ,
    drop = FALSE
  ]
  data$series <- factor(
    ifelse(data$lead_rule == "lfsr",
           names(lead_rule_colors)[1L], names(lead_rule_colors)[2L]),
    levels = names(lead_rule_colors)
  )
  class_forest(data, lead_rule_colors, columns = "class_label")
}

#' Figure D: the two enhancer-state definitions side by side.
plot_definition_sensitivity <- function(strategy = "all") {
  data <- enrichment[
    enrichment$comparator == "Tested background" &
      enrichment$lead_rule == "lfsr" &
      enrichment$strategy == strategy &
      as.character(enrichment$class) %in% TCR_CLASSES,
    ,
    drop = FALSE
  ]
  data$series <- factor(
    ifelse(data$enhancer_definition == "R6",
           names(definition_colors)[1L], names(definition_colors)[2L]),
    levels = names(definition_colors)
  )
  class_forest(data, definition_colors, columns = "class_label")
}

#' A heatmap over class x strategy and epigenome, for one comparator.
plot_fold_heatmap <- function(comparator = "Tested background",
                             definition = "R6",
                             lead_rule = "lfsr") {
  data <- enrichment[
    enrichment$comparator == comparator &
      as.character(enrichment$enhancer_definition) == definition &
      enrichment$lead_rule == lead_rule,
    ,
    drop = FALSE
  ]
  data <- data[data$estimable, , drop = FALSE]
  data$column_label <- factor(
    paste0(as.character(data$class_label), "\n", data$strategy),
    levels = unique(
      paste0(
        as.character(data$class_label[order(data$class, data$strategy)]), "\n",
        data$strategy[order(data$class, data$strategy)]
      )
    )
  )
  limit <- max(abs(log2(data$fold_enrichment)))
  ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = column_label, y = display_label, fill = log2(fold_enrichment)
    )
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = format_decimal(fold_enrichment, 2L),
        fontface = ifelse(ci_lower_fold > 1 | ci_upper_fold < 1,
                          "bold", "plain")
      ),
      size = 2.7, color = "grey15"
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
      axis.text.x = ggplot2::element_text(size = 8.5),
      axis.text.y = ggplot2::element_text(size = 9),
      strip.text.y.left = ggplot2::element_text(angle = 0, face = "bold"),
      strip.placement = "outside",
      legend.position = "right"
    )
}

# ---- Tables -----------------------------------------------------------------

class_table <- function() {
  knitr::kable(
    data.frame(
      Class = class_summary$label,
      `Mutually exclusive` = ifelse(class_summary$mutually_exclusive,
                                    "yes (with the other two windows)", "no"),
      Pairs = format_integer(class_summary$pairs),
      Variants = format_integer(class_summary$variants),
      Genes = format_integer(class_summary$genes),
      `Leads (lfsr)` = format_integer(class_summary$leads),
      `Median class lfsr` = format_decimal(class_summary$median_class_lfsr, 3L),
      check.names = FALSE
    ),
    row.names = FALSE
  )
}

switch_crosstab_table <- function() {
  knitr::kable(
    data.frame(
      `Timing window` = switch_crosstab$label,
      Pairs = format_integer(switch_crosstab$pairs),
      `Also switch` = format_integer(switch_crosstab$switch_pairs),
      `Share` = paste0(format_decimal(100 * switch_crosstab$switch_share, 1), "%"),
      check.names = FALSE
    ),
    row.names = FALSE
  )
}

estimability_table <- function() {
  knitr::kable(
    data.frame(
      Epigenome = estimability$display_label,
      Group = estimability$biological_group,
      `Background rate` = paste0(
        format_decimal(100 * estimability$background_rate, 2L), "%"
      ),
      `Variants needed` = format_integer(estimability$variants_needed),
      check.names = FALSE
    ),
    row.names = FALSE
  )
}

lead_rule_table <- function() {
  knitr::kable(
    data.frame(
      `Lead set` = lead_rule_agreement$variant_set,
      `By lfsr` = format_integer(lead_rule_agreement$n_lfsr),
      `By lfdr` = format_integer(lead_rule_agreement$n_lfdr),
      Shared = format_integer(lead_rule_agreement$n_shared),
      Jaccard = format_decimal(lead_rule_agreement$jaccard, 2L),
      check.names = FALSE
    ),
    row.names = FALSE
  )
}

not_estimable_table <- function() {
  if (!nrow(not_estimable)) {
    return("Every cell cleared the minimum-overlap threshold.")
  }
  rows <- not_estimable[
    not_estimable$enhancer_definition == "R6" &
      not_estimable$lead_rule == "lfsr",
    ,
    drop = FALSE
  ]
  summary_rows <- do.call(rbind, lapply(
    unique(rows$variant_set),
    function(set_name) {
      subset <- rows[rows$variant_set == set_name, , drop = FALSE]
      data.frame(
        `Variant set` = set_name,
        `N variants` = format_integer(subset$selected_total[1L]),
        `Comparators affected` = paste(
          unique(subset$comparator), collapse = "; "
        ),
        `Epigenomes not estimable` = paste(
          sort(unique(subset$epigenome_id)), collapse = ", "
        ),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }
  ))
  knitr::kable(summary_rows, row.names = FALSE)
}

#' The full estimate table for one comparator.
full_summary_table <- function(comparator = "Tested background",
                               definition = "R6",
                               lead_rule = "lfsr") {
  rows <- enrichment[
    enrichment$comparator == comparator &
      as.character(enrichment$enhancer_definition) == definition &
      enrichment$lead_rule == lead_rule,
    ,
    drop = FALSE
  ]
  rows <- rows[order(rows$class, rows$strategy, rows$group_order,
                     rows$epigenome_order), , drop = FALSE]
  knitr::kable(
    data.frame(
      Class = as.character(rows$class_label),
      Selection = rows$strategy,
      ID = rows$epigenome_id,
      `Roadmap name` = rows$roadmap_std_name,
      Group = as.character(rows$biological_group),
      `N variants` = format_integer(rows$selected_total),
      `Overlap prop.` = format_decimal(
        rows$selected_overlap / rows$selected_total, 4L
      ),
      `Comparator prop.` = format_decimal(
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
