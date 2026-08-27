#!/usr/bin/env Rscript

# Cache-dependent checks for the internal temporal-class enhancer-enrichment
# page. Sources reporting.R, which does all the work, and then pins the values
# the page quotes in prose.

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
original_working_directory <- getwd()
on.exit(setwd(original_working_directory), add = TRUE)
setwd(workflowr_root)

source(file.path(
  "code", "revision_simulations", "internal",
  "fash_temporal_class_enhancer_enrichment", "reporting.R"
))

stopifnot(
  # The classification covers exactly the R6 discovery set.
  nrow(labelled_pairs) == 9214L,
  length(discovered_variants) == 9148L,
  identical(class_summary$pairs, c(1157L, 1943L, 6114L)),
  identical(class_summary$variants, c(1157L, 1940L, 6062L)),
  identical(class_summary$genes, c(140L, 346L, 936L)),
  sum(class_summary$pairs) == nrow(labelled_pairs),

  # The strict cfsr classification, quoted as the power note.
  identical(confident_call_counts$pairs, c(126L, 58L, 21L)),
  identical(confident_call_counts$genes, c(8L, 14L, 11L)),

  # Variant sets and their sizes.
  identical(names(variant_sets), c(
    "early_all", "early_lead", "middle_all", "middle_lead",
    "late_all", "late_lead"
  )),
  identical(
    lengths(variant_sets),
    c(early_all = 1157L, early_lead = 140L, middle_all = 1940L,
      middle_lead = 344L, late_all = 6062L, late_lead = 932L)
  ),
  identical(
    lengths(confident_variant_sets),
    c(early_all = 523L, early_lead = 56L, middle_all = 443L,
      middle_lead = 100L, late_all = 2451L, late_lead = 484L)
  ),

  # Every class-and-strategy set is drawn from the discovery set.
  all(vapply(
    variant_sets,
    function(ids) all(ids %in% discovered_variants),
    logical(1)
  )),
  all(vapply(
    seq_along(TEMPORAL_CLASSES),
    function(i) {
      class_name <- TEMPORAL_CLASSES[i]
      all(variant_sets[[paste0(class_name, "_lead")]] %in%
            variant_sets[[paste0(class_name, "_all")]])
    },
    logical(1)
  )),

  # Comparators.
  nrow(background) == 745867L,
  background_variant_count == 745867L,
  nrow(discovery_background) == length(discovered_variants),

  # Table shapes.
  nrow(background_enrichment) == 18L,
  nrow(within_discovery_enrichment) == 18L,
  nrow(matched_enrichment) == 18L,
  nrow(confident_enrichment) == 18L,
  nrow(context_enrichment) == 27L,

  # The headline estimates the prose quotes, against the tested background.
  isTRUE(all.equal(
    enrichment_lookup(background_enrichment, "early_all", e020)$fold_enrichment,
    0.4254795, tolerance = 1e-6
  )),
  isTRUE(all.equal(
    enrichment_lookup(background_enrichment, "early_all", e013)$fold_enrichment,
    0.4356877, tolerance = 1e-6
  )),
  isTRUE(all.equal(
    enrichment_lookup(background_enrichment, "late_all", e013)$fold_enrichment,
    1.2853120, tolerance = 1e-6
  )),
  isTRUE(all.equal(
    enrichment_lookup(background_enrichment, "late_all", e095)$fold_enrichment,
    1.3759418, tolerance = 1e-6
  )),
  isTRUE(all.equal(
    enrichment_lookup(background_enrichment, "late_lead", e013)$fold_enrichment,
    1.4644825, tolerance = 1e-6
  )),
  isTRUE(all.equal(
    enrichment_lookup(background_enrichment, "middle_lead", e013)$fold_enrichment,
    1.4126117, tolerance = 1e-6
  )),

  # The early lead set falls below the overlap floor at E020 and E013, which
  # the page states rather than hides.
  !is.finite(
    enrichment_lookup(background_enrichment, "early_lead", e020)$fold_enrichment
  ),
  !is.finite(
    enrichment_lookup(background_enrichment, "early_lead", e013)$fold_enrichment
  ),
  enrichment_lookup(background_enrichment, "early_lead", e020)$selected_overlap <
    MINIMUM_OVERLAP,

  isTRUE(all.equal(
    enrichment_lookup(
      confident_enrichment, "middle_lead", e013
    )$fold_enrichment,
    2.2427396, tolerance = 1e-6
  )),
  isTRUE(all.equal(
    enrichment_lookup(confident_enrichment, "late_lead", e095)$fold_enrichment,
    1.6868627, tolerance = 1e-6
  )),

  # The within-discovery contrast: early depleted, late enriched.
  enrichment_lookup(
    within_discovery_enrichment, "early_all", e013
  )$ci_upper_fold < 1,
  enrichment_lookup(
    within_discovery_enrichment, "late_all", e095
  )$ci_lower_fold > 1,

  # The matched-control sensitivity keeps the early depletion.
  enrichment_lookup(matched_enrichment, "early_all", e013)$ci_upper_fold < 1,

  # The regulatory-context reading: early depleted at enhancers, enriched at
  # active TSS.
  enrichment_lookup(
    context_enrichment, "early_all", "Roadmap E020 iPS-20b: Active TSS"
  )$ci_lower_fold > 1,
  enrichment_lookup(
    context_enrichment, "early_all", e013
  )$ci_upper_fold < 1,

  # Display tables.
  identical(dim(enrichment_display_table(background_enrichment)), c(6L, 5L)),
  identical(
    enrichment_display_table(background_enrichment)$Class,
    rep(unname(TEMPORAL_CLASS_LABELS[TEMPORAL_CLASSES]), each = 2L)
  ),
  identical(
    dim(enrichment_display_table(context_enrichment, context_annotations)),
    c(3L, 11L)
  ),
  identical(dim(class_summary_display), c(3L, 7L))
)

for (builder in list(
  plot_background_enrichment, plot_within_discovery_enrichment,
  plot_regulatory_context, plot_matched_enrichment
)) {
  plot <- builder()
  stopifnot(inherits(plot, "ggplot"), nrow(plot$data) > 0L)
}

message("fash_temporal_class_enhancer_enrichment/reporting.R: all checks passed.")
