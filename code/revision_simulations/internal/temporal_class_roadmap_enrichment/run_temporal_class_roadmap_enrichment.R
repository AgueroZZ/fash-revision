#!/usr/bin/env Rscript

# =============================================================================
# Temporal-class enrichment across the thirteen-epigenome Roadmap panel.
#
# Steps:
#   1. Load the four cached classifications and the R6 discovery pair table.
#   2. Build class membership and the variant sets, under both lead rules.
#   3. Reuse the 13-epigenome x 2-definition annotation matrix already built by
#      roadmap_enhancer_exploration. Nothing is re-downloaded or re-intersected.
#   4. Run enrichment_api.R::compute_enrichment() three times:
#        (a) 8 class sets vs the tested background, both definitions;
#        (b) 8 class sets + 2 no-switch sets vs the other FASH discoveries;
#        (c) 8 class sets under the lfdr lead rule vs the tested background.
#   5. Cache under
#      output/revision_simulations/internal/temporal_class_roadmap_enrichment/.
#
# Everything read here is read-only: the R6 caches, the two sister experiments'
# code, and the cached FASH fit outputs. Only this experiment's own output
# directory is written.
#
#   Rscript --vanilla run_temporal_class_roadmap_enrichment.R
# =============================================================================

started_at <- Sys.time()

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

source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "fash_strober_enhancer_comparison", "enrichment_api.R"
))
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "temporal_class_roadmap_enrichment", "temporal_class_roadmap_helpers.R"
))

output_directory <- tcr_output_directory()
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

step <- function(index, total, text) {
  message("[", index, "/", total, "] ", text)
}
TOTAL_STEPS <- 6L

EXPECTED_DISCOVERED_PAIRS <- 9214L
EXPECTED_DISCOVERED_VARIANTS <- 9148L
EXPECTED_BACKGROUND_VARIANTS <- 745867L
EXPECTED_CLASS_PAIRS <- c(early = 1157L, middle = 1943L, late = 6114L,
                          switch = 981L)

# -----------------------------------------------------------------------------

step(1, TOTAL_STEPS, "Loading the cached classifications and the R6 discovery set.")

classifications <- load_all_classifications()
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

# Which middle window produced this cache? Tested, not assumed -- the repository
# carries both an open (3, 12) and a stale closed [4, 11] version of the
# functional, and the two are indistinguishable from the saved lfsr alone.
middle_definition <- verify_middle_definition(classifications)
message(
  "  middle window: ", middle_definition$window,
  " (mean residual ", format(middle_definition$mean_residual_open, digits = 2),
  " against 1 - pi0; median pi0 ",
  format(middle_definition$median_pi0, digits = 3), ")"
)

membership <- build_class_membership(discovered_pairs, classifications)
class_summary <- summarise_classes(membership, lead_rule = "lfsr")
switch_crosstab <- crosstab_switch_by_timing(membership)

observed_class_pairs <- stats::setNames(
  as.integer(class_summary$pairs[match(names(EXPECTED_CLASS_PAIRS),
                                       class_summary$class)]),
  names(EXPECTED_CLASS_PAIRS)
)
if (nrow(membership) != EXPECTED_DISCOVERED_PAIRS ||
    length(unique(membership$variant_id)) != EXPECTED_DISCOVERED_VARIANTS ||
    !identical(observed_class_pairs, EXPECTED_CLASS_PAIRS)) {
  stop("The class membership does not match the pre-specified counts.")
}
message(
  "  ", nrow(membership), " discovered pairs; classes ",
  paste(paste0(names(observed_class_pairs), " ", observed_class_pairs),
        collapse = ", "), "."
)

# -----------------------------------------------------------------------------

step(2, TOTAL_STEPS, "Building the variant sets under both lead rules.")

variant_sets <- lapply(names(LEAD_RULES), function(rule) {
  c(
    build_class_variant_sets(membership, lead_rule = rule),
    build_no_switch_sets(membership, lead_rule = rule)
  )
})
names(variant_sets) <- names(LEAD_RULES)

CLASS_SETS <- paste0(rep(TCR_CLASSES, each = 2L), c("_all", "_lead"))
CONTRAST_SETS <- c("no_switch_all", "no_switch_lead")

# The `_all` sets carry the class membership and must not depend on the lead
# rule at all. The `_lead` sets are the whole point of the sensitivity, so they
# are expected to differ.
if (!identical(
  variant_sets$lfsr[c(paste0(TCR_CLASSES, "_all"), "no_switch_all")],
  variant_sets$lfdr[c(paste0(TCR_CLASSES, "_all"), "no_switch_all")]
)) {
  stop("The all-variant sets depend on the lead rule; they must not.")
}
set_sizes <- do.call(rbind, lapply(names(LEAD_RULES), function(rule) {
  data.frame(
    lead_rule = rule,
    variant_set = names(variant_sets[[rule]]),
    n_variants = lengths(variant_sets[[rule]]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  set_sizes, file.path(output_directory, "set_sizes.csv"), row.names = FALSE
)

lead_rule_agreement <- do.call(rbind, lapply(
  c(paste0(TCR_CLASSES, "_lead"), "no_switch_lead"),
  function(set_name) {
    by_lfsr <- variant_sets$lfsr[[set_name]]
    by_lfdr <- variant_sets$lfdr[[set_name]]
    data.frame(
      variant_set = set_name,
      n_lfsr = length(by_lfsr),
      n_lfdr = length(by_lfdr),
      n_shared = length(intersect(by_lfsr, by_lfdr)),
      jaccard = length(intersect(by_lfsr, by_lfdr)) /
        length(union(by_lfsr, by_lfdr)),
      stringsAsFactors = FALSE
    )
  }
))
utils::write.csv(
  lead_rule_agreement,
  file.path(output_directory, "lead_rule_agreement.csv"),
  row.names = FALSE
)
message(
  "  lead sets share ",
  sprintf("%.0f%%-%.0f%%", 100 * min(lead_rule_agreement$jaccard),
          100 * max(lead_rule_agreement$jaccard)),
  " of their variants between the two rules."
)

# -----------------------------------------------------------------------------

step(3, TOTAL_STEPS, "Reusing the thirteen-epigenome annotation matrix.")

panel <- load_epigenome_panel()
annotation_matrix_path <- file.path(
  roadmap_output_directory(), "roadmap_enhancer_annotation_matrix.rds"
)
if (!file.exists(annotation_matrix_path)) {
  stop(
    "The thirteen-epigenome annotation matrix is missing. Run ",
    "code/revision_simulations/internal/roadmap_enhancer_exploration/",
    "run_roadmap_enhancer_exploration.R first."
  )
}
annotation_matrix <- readRDS(annotation_matrix_path)
expected_columns <- unlist(lapply(seq_len(nrow(panel)), function(index) {
  vapply(names(ENHANCER_DEFINITIONS), function(definition) {
    enhancer_annotation_name(
      panel$epigenome_id[index], panel$roadmap_std_name[index], definition
    )
  }, character(1), USE.NAMES = FALSE)
}), use.names = FALSE)
if (!identical(
  names(annotation_matrix), c("variant_id", "chromosome", expected_columns)
)) {
  stop("The annotation matrix is not the thirteen-epigenome table this expects.")
}

background <- load_background()
if (nrow(background) != EXPECTED_BACKGROUND_VARIANTS ||
    !identical(annotation_matrix$variant_id, background$variant_id)) {
  stop("The background and the annotation matrix are not aligned.")
}

# The within-discovery comparator: the same background table restricted to the
# FASH discoveries, so `enrichment_versus_background()` draws its controls from
# the other discovered variants rather than from every tested variant. The API
# is untouched; only what is handed to it changes.
discovered_variants <- unique(membership$variant_id)
discovery_background <- background[
  background$variant_id %in% discovered_variants, , drop = FALSE
]
if (nrow(discovery_background) != EXPECTED_DISCOVERED_VARIANTS) {
  stop("The within-discovery comparator is not the FASH discovery set.")
}

forecast <- estimability_forecast(annotation_matrix, panel = panel)
utils::write.csv(
  forecast, file.path(output_directory, "estimability_forecast.csv"),
  row.names = FALSE
)
message(
  "  a 10-overlap estimate needs ", min(forecast$variants_needed), "-",
  max(forecast$variants_needed), " variants depending on the epigenome."
)

# -----------------------------------------------------------------------------

step(4, TOTAL_STEPS, "Computing enrichment: three runs.")

# compute_enrichment() returns NA when fewer than MINIMUM_OVERLAP selected
# variants land inside an annotation. On this design that is a genuine power
# limit of the smaller classes against the sparser epigenomes, so those cells
# are recorded and reported rather than treated as an error.
finalise <- function(enrichment, comparator, lead_rule) {
  enrichment <- decorate_enrichment(enrichment, panel = panel)
  parsed <- parse_set_name(unique(enrichment$variant_set))
  enrichment <- merge(
    enrichment, parsed, by = "variant_set", all.x = TRUE, sort = FALSE
  )
  if (anyNA(enrichment$class)) {
    stop("A variant set could not be parsed into a class.")
  }
  enrichment$comparator <- comparator
  enrichment$lead_rule <- lead_rule
  enrichment$estimable <- is.finite(enrichment$fold_enrichment) &
    is.finite(enrichment$ci_lower_fold) & is.finite(enrichment$ci_upper_fold)
  unexplained <- !enrichment$estimable &
    enrichment$selected_overlap >= MINIMUM_OVERLAP &
    enrichment$control_overlap >= MINIMUM_OVERLAP
  if (any(unexplained)) {
    stop(
      "Non-finite estimate with adequate overlap for: ",
      paste(paste(enrichment$variant_set[unexplained],
                  enrichment$annotation[unexplained]), collapse = "; ")
    )
  }
  enrichment
}

runs <- list(
  list(
    label = "background",
    sets = variant_sets$lfsr[CLASS_SETS],
    background = background,
    comparator = "Tested background",
    lead_rule = "lfsr"
  ),
  list(
    label = "within_discovery",
    sets = variant_sets$lfsr[c(CLASS_SETS, CONTRAST_SETS)],
    background = discovery_background,
    comparator = "Other FASH discoveries",
    lead_rule = "lfsr"
  ),
  list(
    label = "background_lfdr_leads",
    sets = variant_sets$lfdr[CLASS_SETS],
    background = background,
    comparator = "Tested background",
    lead_rule = "lfdr"
  )
)

results <- lapply(runs, function(run) {
  message("  run: ", run$label, " (", length(run$sets), " sets)")
  finalise(
    compute_enrichment(
      run$sets,
      background = run$background,
      annotations = annotation_matrix,
      controls = "background"
    ),
    comparator = run$comparator,
    lead_rule = run$lead_rule
  )
})
names(results) <- vapply(runs, function(run) run$label, character(1))

enrichment <- do.call(rbind, lapply(names(results), function(label) {
  output <- results[[label]]
  output$run <- label
  output
}))
enrichment <- enrichment[order(
  match(enrichment$run, names(results)),
  enrichment$enhancer_definition,
  match(enrichment$class, c(TCR_CLASSES, "no_switch")),
  enrichment$strategy,
  enrichment$group_order, enrichment$epigenome_order
), , drop = FALSE]
rownames(enrichment) <- NULL

expected_rows <- (length(CLASS_SETS) * 2L +
                    length(CLASS_SETS) + length(CONTRAST_SETS)) *
  nrow(panel) * length(ENHANCER_DEFINITIONS)
if (nrow(enrichment) != expected_rows) {
  stop("Expected ", expected_rows, " estimates, produced ", nrow(enrichment), ".")
}

not_estimable <- enrichment[!enrichment$estimable, c(
  "run", "comparator", "lead_rule", "enhancer_definition", "variant_set",
  "epigenome_id", "roadmap_std_name", "biological_group", "selected_total",
  "selected_overlap"
), drop = FALSE]
utils::write.csv(
  enrichment, file.path(output_directory, "enrichment_results.csv"),
  row.names = FALSE
)
utils::write.csv(
  not_estimable, file.path(output_directory, "not_estimable_cells.csv"),
  row.names = FALSE
)
message(
  "  ", nrow(enrichment), " estimates; ", nrow(not_estimable),
  " below the ", MINIMUM_OVERLAP, "-variant minimum-overlap threshold."
)

# -----------------------------------------------------------------------------

step(5, TOTAL_STEPS, "Cross-checking against the sister pages.")

# The union of the three timing classes is the whole R6 discovery set, so
# pooling them must reproduce R6's own `current_all` estimate exactly. This ties
# this page's numbers to the published page through the estimator itself.
pooled <- compute_enrichment(
  list(pooled_timing = unique(unlist(
    variant_sets$lfsr[paste0(TEMPORAL_CLASSES, "_all")],
    use.names = FALSE
  ))),
  background = background,
  annotations = annotation_matrix,
  controls = "background"
)
published <- utils::read.csv(
  file.path(roadmap_output_directory(), "enrichment_results.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
published <- published[published$variant_set == "current_all", , drop = FALSE]
comparison <- merge(
  pooled[, c("annotation", "fold_enrichment")],
  published[, c("annotation", "fold_enrichment")],
  by = "annotation", suffixes = c("_pooled", "_published")
)
if (nrow(comparison) != nrow(panel) * length(ENHANCER_DEFINITIONS)) {
  stop("The pooled-timing cross-check did not match every annotation.")
}
worst_gap <- max(abs(
  comparison$fold_enrichment_pooled - comparison$fold_enrichment_published
), na.rm = TRUE)
utils::write.csv(
  comparison, file.path(output_directory, "pooled_timing_crosscheck.csv"),
  row.names = FALSE
)
if (worst_gap > 1e-10) {
  stop(
    "Pooling the three timing classes does not reproduce R6's current_all ",
    "estimate; worst gap ", worst_gap, "."
  )
}
message(
  "  pooling early+middle+late reproduces the 13-epigenome current_all ",
  "estimates exactly (worst gap ", format(worst_gap, digits = 3), ")."
)

# -----------------------------------------------------------------------------

step(6, TOTAL_STEPS, "Writing the cache.")

elapsed_minutes <- as.numeric(difftime(Sys.time(), started_at, units = "mins"))
cache <- list(
  configuration = list(
    experiment_id = EXPERIMENT_ID_TCR,
    classes = TCR_CLASSES,
    class_labels = TCR_CLASS_LABELS,
    class_short_labels = TCR_CLASS_SHORT_LABELS,
    switch_threshold = SWITCH_THRESHOLD,
    switch_cfsr_alpha = SWITCH_CFSR_ALPHA,
    lead_rules = LEAD_RULES,
    selection_strategies = SELECTION_STRATEGIES,
    class_sets = CLASS_SETS,
    contrast_sets = CONTRAST_SETS,
    minimum_overlap = MINIMUM_OVERLAP,
    discovered_pairs = nrow(membership),
    discovered_variants = length(discovered_variants),
    background_variant_count = nrow(background),
    r_version = R.version.string,
    built_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    elapsed_minutes = elapsed_minutes
  ),
  panel = panel,
  middle_definition = middle_definition,
  class_summary = class_summary,
  switch_crosstab = switch_crosstab,
  set_sizes = set_sizes,
  lead_rule_agreement = lead_rule_agreement,
  estimability_forecast = forecast,
  pooled_timing_crosscheck = comparison,
  not_estimable = not_estimable,
  enrichment = enrichment
)
saveRDS(cache, file.path(output_directory, "analysis_cache.rds"))
message(
  "Done in ", format(round(elapsed_minutes, 1)), " minutes. Cache: ",
  output_directory
)
