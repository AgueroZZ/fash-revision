#!/usr/bin/env Rscript

# =============================================================================
# The estimates the reviewer-facing R6 page displays.
#
# Enhancers are ChromHMM states 6, 7 and 12 of the Roadmap 15-state core-marks
# model -- `6_EnhG`, `7_Enh`, `12_EnhBiv` -- the definition of He & Wang (2017,
# Bioinformatics 33:3268), across thirteen epigenomes, for one lead variant per
# discovered gene.
#
# Nothing here is new machinery. The epigenome panel, the segmentation reader and
# the annotation builder come from roadmap_enhancer_exploration_helpers.R; the
# background, the discovery sets, the fold-enrichment estimator and the
# leave-one-autosome-out jackknife come from enrichment_api.R. This script only
# selects a different state set and caches what the page reads.
#
# It also cross-checks the new indicators against the already-validated
# two-definition matrix built by the internal exploration: states {6,7,12} must
# contain {6,7} and be contained in {6,7,11,12}, variant by variant, on every
# epigenome. That bracket is what rules out a silent state-selection error.
# =============================================================================

start_time <- Sys.time()

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

WORKFLOWR_ROOT <- find_workflowr_root()
setwd(WORKFLOWR_ROOT)

source(file.path(
  WORKFLOWR_ROOT, "code", "revision_simulations", "internal",
  "roadmap_enhancer_exploration", "roadmap_enhancer_exploration_helpers.R"
))
source(file.path(
  WORKFLOWR_ROOT, "code", "revision_simulations", "internal",
  "fash_strober_enhancer_comparison", "enrichment_api.R"
))

# ---- Scope -------------------------------------------------------------------

# He & Wang (2017) define enhancers in the Roadmap 15-state model as states
# 6, 7 and 12. This is the page's primary definition.
EPICOMPARE_ENHANCER_STATES <- c("6_EnhG", "7_Enh", "12_EnhBiv")

# The broader set in the companion code for Strober et al. (2019), reported on
# the page as a one-sentence sensitivity statement only.
STROBER_SENSITIVITY_STATES <- c("6_EnhG", "7_Enh", "12_EnhBiv", "11_BivFlnk")

DISPLAYED_SETS <- c(
  "current_lead", "linear_lead", "quadratic_lead", "current_only_lead"
)

OUTPUT_DIRECTORY <- file.path(
  WORKFLOWR_ROOT, "output", "revision_simulations", "internal",
  "r6_roadmap_enhancer_panel"
)
dir.create(OUTPUT_DIRECTORY, recursive = TRUE, showWarnings = FALSE)

definitions <- list(
  EpiCompare = EPICOMPARE_ENHANCER_STATES,
  Strober = STROBER_SENSITIVITY_STATES
)
stopifnot(
  all(unlist(definitions) %in% ROADMAP_15_STATE_MNEMONICS),
  # nesting the cross-check below relies on
  all(EPICOMPARE_ENHANCER_STATES %in% STROBER_SENSITIVITY_STATES)
)

message("R6 Roadmap enhancer panel")
message("  primary definition : ", paste(EPICOMPARE_ENHANCER_STATES, collapse = ", "))
message("  sensitivity        : ", paste(STROBER_SENSITIVITY_STATES, collapse = ", "))

# ---- Inputs ------------------------------------------------------------------

panel <- load_epigenome_panel()
background <- load_background()
variant_sets <- load_variant_sets()[DISPLAYED_SETS]

if (nrow(panel) != 13L ||
    anyDuplicated(background$variant_id) ||
    any(vapply(variant_sets, function(x) anyDuplicated(x) > 0L, logical(1L))) ||
    !all(unlist(variant_sets) %in% background$variant_id)) {
  stop("Inputs to the R6 panel are not in the expected shape.")
}
message("  background         : ", nrow(background), " autosomal tested variants")
message("  lead sets          : ",
        paste(sprintf("%s=%d", names(variant_sets), lengths(variant_sets)),
              collapse = ", "))

# ---- Annotation indicators ---------------------------------------------------

annotation_path <- file.path(OUTPUT_DIRECTORY, "annotation_matrix.rds")
if (file.exists(annotation_path)) {
  message("  reusing cached indicator table")
  annotation_matrix <- readRDS(annotation_path)
} else {
  message("  building indicators from the cached segmentations ...")
  annotation_matrix <- build_enhancer_annotation_matrix(
    background[, c("variant_id", "chromosome", "position")],
    panel = panel,
    definitions = definitions
  )
  saveRDS(annotation_matrix, annotation_path)
}

epicompare_columns <- enhancer_annotation_name(
  panel$epigenome_id, panel$roadmap_std_name, "EpiCompare"
)
strober_columns <- enhancer_annotation_name(
  panel$epigenome_id, panel$roadmap_std_name, "Strober"
)
if (!all(c(epicompare_columns, strober_columns) %in% names(annotation_matrix)) ||
    !identical(annotation_matrix$variant_id, background$variant_id)) {
  stop("The indicator table does not match the panel or the background.")
}

# ---- Cross-check against the already-validated internal matrix ---------------
#
# {6,7} subset {6,7,12} subset {6,7,11,12}, variant by variant, per epigenome.
# The outer two indicator sets were built and checked by the internal
# exploration, so bracketing the new one between them is an independent test of
# the state selection rather than a restatement of it.

internal_matrix_path <- file.path(
  roadmap_output_directory(), "roadmap_enhancer_annotation_matrix.rds"
)
if (!file.exists(internal_matrix_path)) {
  stop("Missing the internal exploration indicator table for the cross-check.")
}
internal_matrix <- readRDS(internal_matrix_path)
if (!identical(internal_matrix$variant_id, annotation_matrix$variant_id)) {
  stop("The internal indicator table is not aligned to this one.")
}

nesting <- do.call(rbind, lapply(seq_len(nrow(panel)), function(index) {
  id <- panel$epigenome_id[index]
  name <- panel$roadmap_std_name[index]
  two <- as.logical(internal_matrix[[
    enhancer_annotation_name(id, name, "R6")
  ]])
  three <- as.logical(annotation_matrix[[
    enhancer_annotation_name(id, name, "EpiCompare")
  ]])
  four <- as.logical(internal_matrix[[
    enhancer_annotation_name(id, name, "Strober")
  ]])
  data.frame(
    epigenome_id = id,
    n_two = sum(two), n_three = sum(three), n_four = sum(four),
    two_not_in_three = sum(two & !three),
    three_not_in_four = sum(three & !four),
    stringsAsFactors = FALSE
  )
}))
if (any(nesting$two_not_in_three != 0L) ||
    any(nesting$three_not_in_four != 0L) ||
    any(nesting$n_three < nesting$n_two) ||
    any(nesting$n_three > nesting$n_four)) {
  stop("The three-state indicators are not bracketed by {6,7} and {6,7,11,12}.")
}
message("  nesting check      : {6,7} <= {6,7,12} <= {6,7,11,12} on all 13")
rm(internal_matrix)
invisible(gc())

# ---- Estimates ---------------------------------------------------------------

# `decorate_enrichment()` recovers the epigenome and the definition from the
# annotation column name, so the two definitions are estimated, stacked, and
# labelled in one pass rather than tagged by hand.
estimate <- function(columns) {
  compute_enrichment(
    variant_sets,
    background = background,
    annotations = annotation_matrix[, c("variant_id", "chromosome", columns)],
    controls = "background"
  )
}

message("  computing enrichment ...")
enrichment <- decorate_enrichment(
  rbind(estimate(epicompare_columns), estimate(strober_columns)),
  panel = panel,
  definitions = definitions
)

# A cell has no estimate only when fewer than MINIMUM_OVERLAP selected variants
# fall inside the annotation -- a real limit of the sparsest epigenomes against
# the smallest sets. Anything else non-finite would be a bug, so it stops here.
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

expected_rows <- 13L * length(DISPLAYED_SETS) * 2L
if (nrow(enrichment) != expected_rows ||
    !setequal(unique(enrichment$epigenome_id), panel$epigenome_id) ||
    !setequal(unique(enrichment$variant_set), DISPLAYED_SETS)) {
  stop("The enrichment table is not the expected ", expected_rows, " rows.")
}

primary <- enrichment[enrichment$enhancer_definition == "EpiCompare", ]
message("  estimable cells    : ", sum(primary$estimable), " of ", nrow(primary))

# ---- How far the sensitivity definition moves the primary estimates ----------

wide <- merge(
  primary[, c("variant_set", "epigenome_id", "fold_enrichment", "ci_lower_fold",
              "estimable")],
  enrichment[enrichment$enhancer_definition == "Strober",
             c("variant_set", "epigenome_id", "fold_enrichment", "ci_lower_fold",
               "estimable")],
  by = c("variant_set", "epigenome_id"), suffixes = c("_primary", "_broader")
)
wide <- wide[wide$estimable_primary & wide$estimable_broader, , drop = FALSE]
sensitivity <- list(
  n_cells = nrow(wide),
  correlation = stats::cor(wide$fold_enrichment_primary,
                           wide$fold_enrichment_broader),
  max_absolute_shift = max(abs(wide$fold_enrichment_primary -
                                 wide$fold_enrichment_broader)),
  interval_crossings = sum(xor(wide$ci_lower_fold_primary > 1,
                               wide$ci_lower_fold_broader > 1))
)
message("  sensitivity        : r=", round(sensitivity$correlation, 3),
        ", max shift ", round(sensitivity$max_absolute_shift, 3), "-fold, ",
        sensitivity$interval_crossings, " interval crossing(s)")

# ---- Cache -------------------------------------------------------------------

configuration <- list(
  built_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  r_version = R.version.string,
  elapsed_minutes = as.numeric(
    difftime(Sys.time(), start_time, units = "mins")
  ),
  primary_states = EPICOMPARE_ENHANCER_STATES,
  sensitivity_states = STROBER_SENSITIVITY_STATES,
  displayed_sets = DISPLAYED_SETS,
  expected_variant_counts = lengths(variant_sets),
  background_variant_count = nrow(background),
  minimum_overlap = MINIMUM_OVERLAP
)

utils::write.csv(nesting, file.path(OUTPUT_DIRECTORY, "state_nesting_check.csv"),
                 row.names = FALSE)
utils::write.csv(enrichment, file.path(OUTPUT_DIRECTORY, "enrichment_results.csv"),
                 row.names = FALSE)
saveRDS(
  list(
    configuration = configuration,
    panel = panel,
    enrichment = enrichment,
    nesting = nesting,
    sensitivity = sensitivity
  ),
  file.path(OUTPUT_DIRECTORY, "analysis_cache.rds")
)

message("Done in ", round(configuration$elapsed_minutes, 1), " min -> ",
        OUTPUT_DIRECTORY)
