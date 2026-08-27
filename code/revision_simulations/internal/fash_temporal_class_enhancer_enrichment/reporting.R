# Reporting layer for the internal temporal-class enhancer-enrichment page.
#
# R6 established that FASH's dynamic-eQTL discoveries fall inside Roadmap
# enhancers of this differentiation more often than the tested background does.
# This page splits those same discoveries into the manuscript's early, middle,
# and late timing classes and asks whether the enrichment sits in a specific
# enhancer category.
#
# Nothing is refit here. The two cached inputs are
#   * the R6 discovery pair table (BF-adjusted FASH-IWP1, cumulative-lfdr 0.05)
#   * the cached testing_functional() results for the three timing windows
# and the enrichment itself runs through enrichment_api.R, unchanged, so these
# estimates are on exactly the same footing as the published R6 page.

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
  "fash_strober_enhancer_comparison", "enrichment_api.R"
))
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "fash_temporal_class_enhancer_enrichment", "temporal_class_helpers.R"
))

# ---- Fixed scope ------------------------------------------------------------

analysis_alpha <- 0.05
confident_probability_floor <- 0.7
expected_discovered_pairs <- 9214L
expected_discovered_variants <- 9148L

# Starting, intermediate, and terminal states of the differentiation, in the
# same labelling R6 uses.
roadmap_enhancers <- c(
  "Roadmap E020 iPS-20b: Enhancer" = "E020 iPSC enhancer",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer" = "E013 mesoderm enhancer",
  "Roadmap E095 left ventricle: Enhancer" = "E095 left-ventricle enhancer"
)

# The context panel: if a class is depleted at enhancers, these say where it
# sits instead.
context_annotations <- c(
  "Roadmap E020 iPS-20b: Enhancer" = "E020 enhancer",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer" = "E013 enhancer",
  "Roadmap E095 left ventricle: Enhancer" = "E095 enhancer",
  "ENCODE cCRE enhancer-like" = "cCRE enhancer-like",
  "Roadmap E020 iPS-20b: Active TSS" = "E020 active TSS",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Active TSS" = "E013 active TSS",
  "Roadmap E095 left ventricle: Active TSS" = "E095 active TSS",
  "ENCODE cCRE promoter-like" = "cCRE promoter-like",
  "GENCODE promoter +/-2 kb" = "GENCODE promoter"
)
context_groups <- c("Enhancer", "Enhancer", "Enhancer", "Enhancer",
                    "Promoter / TSS", "Promoter / TSS", "Promoter / TSS",
                    "Promoter / TSS", "Promoter / TSS")

strategy_labels <- c(all = "All variants", lead = "One lead variant per gene")

e020 <- "Roadmap E020 iPS-20b: Enhancer"
e013 <- "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer"
e095 <- "Roadmap E095 left ventricle: Enhancer"

# ---- Cached inputs ----------------------------------------------------------

discovery_pair_path <- file.path(
  workflowr_root, "output", "revision_simulations", "internal",
  "fash_strober_enhancer_comparison_fashr0143", "discovery_pair_tables.rds"
)
if (!file.exists(discovery_pair_path)) {
  stop(
    "The R6 discovery pair tables are missing. Run ",
    "code/revision_simulations/internal/fash_strober_enhancer_comparison/",
    "run_fash_strober_enhancer_comparison.R first."
  )
}
discovered_pairs <- readRDS(discovery_pair_path)$current_all

classification_directory <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real"
)
testing_results <- lapply(TEMPORAL_CLASSES, function(class_name) {
  path <- file.path(
    classification_directory,
    paste0("classify_dyn_eQTLs_", class_name, ".RData")
  )
  if (!file.exists(path)) {
    stop(
      "Missing cached classification: ", path,
      "\nRun code/02_dyn_lfsr.R first."
    )
  }
  environment_for_load <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment_for_load)
  object_name <- paste0("testing_", class_name, "_dyn")
  if (!identical(loaded, object_name)) {
    stop("Unexpected contents in ", path, ": ", paste(loaded, collapse = ", "))
  }
  get(object_name, envir = environment_for_load)
})
names(testing_results) <- TEMPORAL_CLASSES

# ---- Classification ---------------------------------------------------------

window_probabilities <- derive_window_probabilities(testing_results)
class_assignment <- assign_temporal_class(window_probabilities)
# label_discovered_pairs() fails unless the classified pairs are exactly the R6
# discovered pairs, which is what ties the cached lfsr values to this fit.
labelled_pairs <- label_discovered_pairs(discovered_pairs, class_assignment)

class_summary <- summarise_temporal_classes(labelled_pairs)
confident_call_counts <- count_confident_calls(testing_results, analysis_alpha)
probability_row_sums <- rowSums(window_probabilities)

variant_sets <- build_temporal_variant_sets(labelled_pairs)
confident_variant_sets <- build_temporal_variant_sets(
  labelled_pairs,
  minimum_probability = confident_probability_floor
)
discovered_variants <- unique(labelled_pairs$variant_id)

if (nrow(labelled_pairs) != expected_discovered_pairs ||
    length(discovered_variants) != expected_discovered_variants ||
    sum(class_summary$pairs) != expected_discovered_pairs ||
    any(class_summary$pairs < 1L) ||
    any(probability_row_sums > 1.1) ||
    any(probability_row_sums < 0.5)) {
  stop("The temporal classification of the R6 discovery set failed validation.")
}

# ---- Enrichment -------------------------------------------------------------

background <- load_background()
background_variant_count <- nrow(background)
enhancer_annotations <- load_annotation_groups(
  "custom", groups = names(roadmap_enhancers)
)
context_annotation_table <- load_annotation_groups(
  "custom", groups = names(context_annotations)
)

# The comparator restricted to FASH's own discoveries. A class is then measured
# against the other discovered variants rather than against the whole testing
# universe, which is what "specific to this class" means.
discovery_background <- background[
  background$variant_id %in% discovered_variants, , drop = FALSE
]

decorate <- function(enrichment, labels) {
  parts <- strsplit(enrichment$variant_set, "_", fixed = TRUE)
  enrichment$class <- factor(
    vapply(parts, `[[`, character(1), 1L),
    levels = TEMPORAL_CLASSES,
    labels = unname(TEMPORAL_CLASS_LABELS[TEMPORAL_CLASSES])
  )
  enrichment$strategy <- factor(
    unname(strategy_labels[vapply(parts, `[[`, character(1), 2L)]),
    levels = unname(strategy_labels)
  )
  enrichment$annotation_label <- factor(
    unname(labels[enrichment$annotation]),
    levels = unname(labels)
  )
  enrichment
}

background_enrichment <- decorate(
  compute_enrichment(
    variant_sets, background, enhancer_annotations, controls = "background"
  ),
  roadmap_enhancers
)

within_discovery_enrichment <- decorate(
  compute_enrichment(
    variant_sets, discovery_background, enhancer_annotations,
    controls = "background"
  ),
  roadmap_enhancers
)

# The context panel uses the all-variant view only: the lead sets are too small
# to support nine annotations at the 10-overlap floor.
context_enrichment <- decorate(
  compute_enrichment(
    variant_sets[paste0(TEMPORAL_CLASSES, "_all")],
    background, context_annotation_table, controls = "background"
  ),
  context_annotations
)
context_enrichment$annotation_group <- factor(
  context_groups[match(context_enrichment$annotation, names(context_annotations))],
  levels = c("Enhancer", "Promoter / TSS")
)

confident_enrichment <- decorate(
  compute_enrichment(
    confident_variant_sets, background, enhancer_annotations,
    controls = "background"
  ),
  roadmap_enhancers
)

# Matched controls: 5 comparable variants per selected variant on the MAF,
# TSS-distance, variant-density and tested-gene strata, pooled over 100 seeds.
# Same sampler R6 retains but does not display; used here only to check that
# the early depletion is not a MAF / TSS-distance artefact.
matched_enrichment <- decorate(
  compute_enrichment(
    variant_sets, background, enhancer_annotations, controls = "matched"
  ),
  roadmap_enhancers
)

if (nrow(background_enrichment) != length(variant_sets) * length(roadmap_enhancers) ||
    nrow(within_discovery_enrichment) != nrow(background_enrichment) ||
    nrow(matched_enrichment) != nrow(background_enrichment) ||
    nrow(context_enrichment) !=
      length(TEMPORAL_CLASSES) * length(context_annotations) ||
    anyNA(background_enrichment$class) ||
    anyNA(background_enrichment$strategy) ||
    anyNA(background_enrichment$annotation_label)) {
  stop("The computed enrichment tables are incomplete.")
}

# ---- Formatting -------------------------------------------------------------

format_integer <- function(value) {
  format(as.integer(round(value)), big.mark = ",", scientific = FALSE)
}

format_decimal <- function(value, digits = 2L) {
  formatC(as.numeric(value), format = "f", digits = digits)
}

enrichment_lookup <- function(table, set_name, annotation) {
  row <- table[
    table$variant_set == set_name & table$annotation == annotation, , drop = FALSE
  ]
  if (nrow(row) != 1L) {
    stop("Enrichment lookup failed for ", set_name, " / ", annotation, ".")
  }
  row
}

# "1.38-fold (95% CI 1.12-1.68)", or an explicit note when the overlap floor
# suppressed the estimate.
describe_fold <- function(set_name,
                          annotation,
                          table = background_enrichment,
                          digits = 2L) {
  row <- enrichment_lookup(table, set_name, annotation)
  if (!is.finite(row$fold_enrichment)) {
    return(paste0(
      "not estimated (", row$selected_overlap, " overlapping variants, below ",
      "the floor of ", MINIMUM_OVERLAP, ")"
    ))
  }
  paste0(
    format_decimal(row$fold_enrichment, digits), "-fold (95% CI ",
    format_decimal(row$ci_lower_fold, digits), "-",
    format_decimal(row$ci_upper_fold, digits), ")"
  )
}

fold_value <- function(set_name,
                       annotation,
                       table = background_enrichment,
                       digits = 2L) {
  row <- enrichment_lookup(table, set_name, annotation)
  if (!is.finite(row$fold_enrichment)) {
    return("n/a")
  }
  format_decimal(row$fold_enrichment, digits)
}

set_size <- function(set_name) {
  format_integer(length(variant_sets[[set_name]]))
}

class_count <- function(class_name, column) {
  format_integer(class_summary[[column]][class_summary$class == class_name])
}

# One display table per enrichment run: fold and interval as text, classes down
# the rows, annotations across the columns.
enrichment_display_table <- function(table, labels = roadmap_enhancers) {
  entry <- ifelse(
    is.finite(table$fold_enrichment),
    paste0(
      format_decimal(table$fold_enrichment), " (",
      format_decimal(table$ci_lower_fold), "-",
      format_decimal(table$ci_upper_fold), ")"
    ),
    "not estimated"
  )
  wide <- data.frame(
    Class = as.character(table$class),
    Selection = as.character(table$strategy),
    Annotation = as.character(table$annotation_label),
    Estimate = entry,
    stringsAsFactors = FALSE
  )
  reshaped <- stats::reshape(
    wide,
    idvar = c("Class", "Selection"),
    timevar = "Annotation",
    direction = "wide"
  )
  names(reshaped) <- sub("^Estimate\\.", "", names(reshaped))
  ordering <- order(
    match(reshaped$Class, unname(TEMPORAL_CLASS_LABELS[TEMPORAL_CLASSES])),
    match(reshaped$Selection, unname(strategy_labels))
  )
  reshaped <- reshaped[ordering, c("Class", "Selection", unname(labels)),
                       drop = FALSE]
  row.names(reshaped) <- NULL
  reshaped
}

class_summary_display <- data.frame(
  Class = class_summary$label,
  `Discovered pairs` = vapply(class_summary$pairs, format_integer, character(1)),
  `Unique variants` = vapply(class_summary$variants, format_integer, character(1)),
  `Unique genes` = vapply(class_summary$genes, format_integer, character(1)),
  `Median assignment probability` = format_decimal(
    class_summary$median_probability
  ),
  `Pairs with cfsr <= 0.05` = vapply(
    confident_call_counts$pairs, format_integer, character(1)
  ),
  `Genes with cfsr <= 0.05` = vapply(
    confident_call_counts$genes, format_integer, character(1)
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

render_internal_table <- function(table, align = NULL, minimum_width = "760px") {
  rendered <- knitr::kable(
    table,
    format = "html",
    align = align,
    escape = TRUE,
    table.attr = paste0(
      'class="table table-striped table-hover" style="min-width:',
      minimum_width, ';"'
    )
  )
  cat(
    '<div class="temporal-class-table-scroll">\n',
    as.character(rendered),
    '\n</div>\n',
    sep = ""
  )
  invisible(rendered)
}

# ---- Figures ----------------------------------------------------------------

class_colors <- c("#0072B2", "#009E73", "#D55E00")
names(class_colors) <- unname(TEMPORAL_CLASS_LABELS[TEMPORAL_CLASSES])

# One shared grammar for every figure: annotations down the rows, classes in
# colour, fold enrichment on the x-axis, dashed reference line at 1, and
# leave-one-autosome-out jackknife intervals. Non-estimable cells simply have
# no point, which is why they are also stated in the text.
class_forest <- function(data, facet_by = "strategy", x_label,
                         free_y = FALSE) {
  data <- data[is.finite(data$fold_enrichment), , drop = FALSE]
  data$annotation_label <- factor(
    as.character(data$annotation_label),
    levels = rev(levels(data$annotation_label))
  )
  dodge <- ggplot2::position_dodge(width = 0.6)
  plot <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = fold_enrichment,
      y = annotation_label,
      color = class,
      group = class
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 1, linetype = "dashed", color = "grey45", linewidth = 0.45
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = ci_lower_fold, xmax = ci_upper_fold),
      orientation = "y", width = 0.25, linewidth = 0.5, position = dodge
    ) +
    ggplot2::geom_point(size = 2.1, position = dodge) +
    ggplot2::scale_color_manual(values = class_colors, drop = FALSE) +
    ggplot2::labs(x = x_label, y = NULL, color = NULL) +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10)
    )
  if (!is.null(facet_by)) {
    plot <- plot + ggplot2::facet_wrap(
      stats::as.formula(paste("~", facet_by)),
      nrow = 1L,
      scales = if (free_y) "free_y" else "fixed"
    )
  }
  plot
}

plot_background_enrichment <- function() {
  class_forest(
    background_enrichment,
    x_label = "Fold enrichment relative to all tested variants"
  )
}

plot_within_discovery_enrichment <- function() {
  class_forest(
    within_discovery_enrichment,
    x_label = "Fold enrichment relative to the other FASH discoveries"
  )
}

plot_regulatory_context <- function() {
  class_forest(
    context_enrichment,
    facet_by = "annotation_group",
    x_label = "Fold enrichment relative to all tested variants",
    free_y = TRUE
  )
}

plot_matched_enrichment <- function() {
  class_forest(
    matched_enrichment,
    x_label = "Fold enrichment relative to matched control variants"
  )
}
