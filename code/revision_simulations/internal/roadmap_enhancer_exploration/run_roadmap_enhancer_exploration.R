#!/usr/bin/env Rscript

# =============================================================================
# Internal Roadmap enhancer exploration: build the cache the page reads.
#
# Steps:
#   1. Download and checksum thirteen Roadmap 15-state ChromHMM segmentations.
#   2. Record per-state interval counts and covered bases, for the record.
#   3. Build a variant x (epigenome, enhancer definition) indicator table over
#      the same tested-variant universe R6 uses.
#   4. Run enrichment_api.R::compute_enrichment() on the eight published
#      discovery sets against all 26 annotation columns, with the tested-variant
#      background and the leave-one-autosome-out jackknife. Nothing about the
#      estimator changes; only the annotation table is wider.
#   5. Cache everything under
#      output/revision_simulations/internal/roadmap_enhancer_exploration/.
#
# Nothing outside that output directory and
# data/revision_simulations/internal/roadmap_enhancer_exploration/ is written.
# The R6 page, its code, and its retained cache are read-only here.
#
#   Rscript --vanilla run_roadmap_enhancer_exploration.R
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

# The estimator, the background and the discovery sets, all unmodified.
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "fash_strober_enhancer_comparison", "enrichment_api.R"
))
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "roadmap_enhancer_exploration", "roadmap_enhancer_exploration_helpers.R"
))

output_directory <- roadmap_output_directory()
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

step <- function(index, total, text) {
  message("[", index, "/", total, "] ", text)
}

TOTAL_STEPS <- 6L

# The eight sets the R6 page reports, in its display order.
EXPECTED_SETS <- c(
  "current_all", "current_lead", "linear_all", "linear_lead",
  "quadratic_all", "quadratic_lead", "current_only_all", "current_only_lead"
)
EXPECTED_VARIANT_COUNTS <- c(
  current_all = 9148L, current_lead = 1169L,
  linear_all = 5387L, linear_lead = 548L,
  quadratic_all = 6797L, quadratic_lead = 690L,
  current_only_all = 7981L, current_only_lead = 1029L
)
EXPECTED_BACKGROUND_VARIANTS <- 745867L

# -----------------------------------------------------------------------------

step(1, TOTAL_STEPS, "Loading the epigenome panel and the R6 inputs.")

panel <- load_epigenome_panel()
background <- load_background()
variant_sets <- load_variant_sets()

if (nrow(background) != EXPECTED_BACKGROUND_VARIANTS) {
  stop(
    "The tested-variant background has ", nrow(background),
    " rows, expected ", EXPECTED_BACKGROUND_VARIANTS, "."
  )
}
if (!identical(
  stats::setNames(lengths(variant_sets[EXPECTED_SETS]), EXPECTED_SETS),
  EXPECTED_VARIANT_COUNTS
)) {
  stop("The published discovery sets are not the ones this script expects.")
}
message(
  "  ", nrow(panel), " epigenomes, ", length(ENHANCER_DEFINITIONS),
  " enhancer definitions, ", length(EXPECTED_SETS), " discovery sets, ",
  format(nrow(background), big.mark = ","), " background variants."
)

# -----------------------------------------------------------------------------

step(2, TOTAL_STEPS, "Downloading and checksumming the Roadmap segmentations.")

annotation_provenance <- download_roadmap_segmentations(panel)
utils::write.csv(
  annotation_provenance,
  file.path(output_directory, "annotation_provenance.csv"),
  row.names = FALSE
)
message(
  "  ", nrow(annotation_provenance), " segmentation files, ",
  format(round(sum(annotation_provenance$byte_size) / 1024^2, 1)), " MiB total."
)

# -----------------------------------------------------------------------------

step(3, TOTAL_STEPS, "Summarising the ChromHMM states in each segmentation.")

state_summary <- do.call(rbind, lapply(seq_len(nrow(panel)), function(index) {
  summarise_roadmap_states(
    file.path(roadmap_data_directory(), panel$local_filename[index]),
    panel$epigenome_id[index]
  )
}))
utils::write.csv(
  state_summary,
  file.path(output_directory, "chromhmm_state_summary.csv"),
  row.names = FALSE
)

# How much the Strober definition adds over the R6 one, in genome coverage.
definition_coverage <- do.call(rbind, lapply(panel$epigenome_id, function(id) {
  rows <- state_summary[state_summary$epigenome_id == id, , drop = FALSE]
  data.frame(
    epigenome_id = id,
    r6_covered_bases = sum(rows$autosomal_covered_bases[rows$in_r6_definition]),
    strober_covered_bases =
      sum(rows$autosomal_covered_bases[rows$in_strober_definition]),
    stringsAsFactors = FALSE
  )
}))
definition_coverage$strober_over_r6 <-
  definition_coverage$strober_covered_bases / definition_coverage$r6_covered_bases
utils::write.csv(
  definition_coverage,
  file.path(output_directory, "definition_coverage.csv"),
  row.names = FALSE
)
message(
  "  Strober's four states cover ",
  sprintf("%.3f", min(definition_coverage$strober_over_r6)), "-",
  sprintf("%.3f", max(definition_coverage$strober_over_r6)),
  " times the autosomal bases of R6's two, across the panel."
)

# -----------------------------------------------------------------------------

step(4, TOTAL_STEPS, "Building the variant x annotation indicator table.")

# Building this table is the expensive step (thirteen genome-wide interval
# intersections), and it depends only on the panel and the tested-variant
# universe, neither of which changes between runs. So a valid cached copy is
# reused and re-validated below rather than rebuilt.
annotation_matrix_path <- file.path(
  output_directory, "roadmap_enhancer_annotation_matrix.rds"
)
expected_annotation_columns <- unlist(lapply(
  seq_len(nrow(panel)),
  function(index) {
    vapply(names(ENHANCER_DEFINITIONS), function(definition) {
      enhancer_annotation_name(
        panel$epigenome_id[index], panel$roadmap_std_name[index], definition
      )
    }, character(1), USE.NAMES = FALSE)
  }
), use.names = FALSE)

annotation_matrix <- NULL
if (file.exists(annotation_matrix_path)) {
  candidate <- readRDS(annotation_matrix_path)
  if (identical(
    names(candidate),
    c("variant_id", "chromosome", expected_annotation_columns)
  ) && identical(candidate$variant_id, background$variant_id)) {
    message("  reusing the cached annotation matrix.")
    annotation_matrix <- candidate
  }
  rm(candidate)
  invisible(gc())
}
if (is.null(annotation_matrix)) {
  annotation_matrix <- build_enhancer_annotation_matrix(
    background[, c("variant_id", "chromosome", "position"), drop = FALSE],
    panel = panel
  )
  saveRDS(annotation_matrix, annotation_matrix_path)
}

annotation_columns <- setdiff(
  names(annotation_matrix), c("variant_id", "chromosome")
)
if (length(annotation_columns) !=
      nrow(panel) * length(ENHANCER_DEFINITIONS)) {
  stop("Unexpected number of annotation columns.")
}
# The R6 definition is nested inside the Strober one, so every R6 indicator must
# imply its Strober counterpart. A violation would mean the state filter or the
# interval assembly is wrong.
for (index in seq_len(nrow(panel))) {
  r6_column <- enhancer_annotation_name(
    panel$epigenome_id[index], panel$roadmap_std_name[index], "R6"
  )
  strober_column <- enhancer_annotation_name(
    panel$epigenome_id[index], panel$roadmap_std_name[index], "Strober"
  )
  if (any(annotation_matrix[[r6_column]] & !annotation_matrix[[strober_column]])) {
    stop(
      "R6 enhancer indicator is not nested inside the Strober one for ",
      panel$epigenome_id[index], "."
    )
  }
}
message(
  "  ", format(nrow(annotation_matrix), big.mark = ","), " variants x ",
  length(annotation_columns), " annotation columns; nesting check passed."
)

# Cross-check against the three columns the R6 page already uses. The existing
# matrix collapses states 6 and 7 into one "Enhancer" category, which is exactly
# the R6 definition here, so the two indicators must agree variant for variant.
r6_published <- load_annotation_groups("custom", groups = c(
  "Roadmap E020 iPS-20b: Enhancer",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer",
  "Roadmap E095 left ventricle: Enhancer"
))
r6_crosswalk <- c(
  "Roadmap E020 iPS-20b: Enhancer" = "E020",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer" = "E013",
  "Roadmap E095 left ventricle: Enhancer" = "E095"
)
rows <- match(r6_published$variant_id, annotation_matrix$variant_id)
if (anyNA(rows)) {
  stop("The published annotation matrix covers variants this one does not.")
}
r6_agreement <- do.call(rbind, lapply(names(r6_crosswalk), function(published) {
  id <- r6_crosswalk[[published]]
  std_name <- panel$roadmap_std_name[panel$epigenome_id == id]
  here <- annotation_matrix[[enhancer_annotation_name(id, std_name, "R6")]][rows]
  there <- as.logical(r6_published[[published]])
  data.frame(
    epigenome_id = id,
    published_column = published,
    n_variants = length(there),
    published_overlap = sum(there),
    reproduced_overlap = sum(here),
    n_disagreements = sum(xor(here, there)),
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  r6_agreement,
  file.path(output_directory, "r6_definition_agreement.csv"),
  row.names = FALSE
)
if (any(r6_agreement$n_disagreements != 0L)) {
  stop(
    "The R6-definition indicators here disagree with the published matrix: ",
    paste(r6_agreement$epigenome_id[r6_agreement$n_disagreements != 0L],
          collapse = ", ")
  )
}
message(
  "  E020/E013/E095 R6-definition indicators reproduce the published matrix ",
  "exactly (0 disagreements over ",
  format(r6_agreement$n_variants[1L], big.mark = ","), " variants)."
)
rm(r6_published)
invisible(gc())

# -----------------------------------------------------------------------------

step(5, TOTAL_STEPS, "Computing enrichment for 8 discovery sets x 26 annotations.")

enrichment <- compute_enrichment(
  variant_sets[EXPECTED_SETS],
  background = background,
  annotations = annotation_matrix,
  controls = "background"
)
enrichment <- decorate_enrichment(enrichment, panel = panel)

# The labels the page's figures need: which method produced the set, whether it
# is the all-variant or lead-per-gene view, and whether the Strober rsIDs were
# removed. Derived from the set names, not from the R6 cache, so this script has
# no write dependency on anything R6 owns.
set_metadata <- data.frame(
  variant_set = EXPECTED_SETS,
  method = c(
    "Current FASH", "Current FASH", "Strober linear", "Strober linear",
    "Strober quadratic/nonlinear", "Strober quadratic/nonlinear",
    "Current FASH", "Current FASH"
  ),
  selection_strategy = c(
    "All variants", "One lead variant per gene", "All variants",
    "One lead variant per gene", "All variants", "One lead variant per gene",
    "All variants", "One lead variant per gene"
  ),
  exclusivity = c(
    "None", "None", "None", "None", "None", "None",
    "Variant-level", "Variant-level"
  ),
  section = c(1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L),
  stringsAsFactors = FALSE
)
enrichment <- merge(
  enrichment, set_metadata, by = "variant_set", all.x = TRUE, sort = FALSE
)
expected_rows <- length(EXPECTED_SETS) * nrow(panel) *
  length(ENHANCER_DEFINITIONS)
if (nrow(enrichment) != expected_rows || anyNA(enrichment$method)) {
  stop("The enrichment estimates are incomplete.")
}

# compute_enrichment() returns NA rather than an estimate when fewer than
# MINIMUM_OVERLAP selected variants fall inside an annotation, which is a real
# limit of the sparsest epigenomes against the smallest discovery sets, not a
# failure. Those cells are recorded, reported on the page as not estimable, and
# excluded from the figures. Anything else non-finite would be a bug.
enrichment$estimable <- is.finite(enrichment$fold_enrichment) &
  is.finite(enrichment$ci_lower_fold) & is.finite(enrichment$ci_upper_fold)
unexplained <- !enrichment$estimable &
  enrichment$selected_overlap >= MINIMUM_OVERLAP &
  enrichment$control_overlap >= MINIMUM_OVERLAP
if (any(unexplained)) {
  stop(
    "Non-finite estimate with adequate overlap for: ",
    paste(
      paste(enrichment$variant_set[unexplained],
            enrichment$annotation[unexplained]),
      collapse = "; "
    )
  )
}
not_estimable <- enrichment[!enrichment$estimable, c(
  "enhancer_definition", "variant_set", "epigenome_id", "roadmap_std_name",
  "biological_group", "selected_total", "selected_overlap"
), drop = FALSE]
utils::write.csv(
  not_estimable,
  file.path(output_directory, "not_estimable_cells.csv"),
  row.names = FALSE
)
message(
  "  ", nrow(not_estimable), " of ", expected_rows,
  " cells fall below the ", MINIMUM_OVERLAP,
  "-variant minimum-overlap threshold and are not estimable."
)
enrichment <- enrichment[order(
  enrichment$enhancer_definition, enrichment$section,
  match(enrichment$variant_set, EXPECTED_SETS),
  enrichment$group_order, enrichment$epigenome_order
), , drop = FALSE]
rownames(enrichment) <- NULL

utils::write.csv(
  enrichment,
  file.path(output_directory, "enrichment_results.csv"),
  row.names = FALSE
)
summary_table <- build_summary_table(enrichment)
utils::write.csv(
  summary_table,
  file.path(output_directory, "enrichment_summary_table.csv"),
  row.names = FALSE
)
message("  ", nrow(enrichment), " estimates written.")

# One independent check on the estimator: for the three epigenomes and the R6
# state definition, re-running through the published annotation matrix must give
# the same fold enrichment the R6 page prints.
published_check <- compute_enrichment(
  variant_sets[EXPECTED_SETS],
  background = background,
  annotations = load_annotation_groups("custom", groups = names(r6_crosswalk)),
  controls = "background"
)
published_check$epigenome_id <- unname(r6_crosswalk[published_check$annotation])
comparison <- merge(
  published_check[, c("variant_set", "epigenome_id", "fold_enrichment")],
  enrichment[
    enrichment$enhancer_definition == "R6" &
      enrichment$epigenome_id %in% unname(r6_crosswalk),
    c("variant_set", "epigenome_id", "fold_enrichment")
  ],
  by = c("variant_set", "epigenome_id"),
  suffixes = c("_published", "_here")
)
if (nrow(comparison) != length(EXPECTED_SETS) * length(r6_crosswalk)) {
  stop("The R6 reproduction check did not match every estimate.")
}
worst_gap <- max(
  abs(comparison$fold_enrichment_published - comparison$fold_enrichment_here),
  na.rm = TRUE
)
if (!identical(
  is.finite(comparison$fold_enrichment_published),
  is.finite(comparison$fold_enrichment_here)
)) {
  stop("The R6 reproduction check disagrees about which cells are estimable.")
}
utils::write.csv(
  comparison,
  file.path(output_directory, "r6_estimate_reproduction.csv"),
  row.names = FALSE
)
if (worst_gap > 1e-10) {
  stop("The R6 estimates do not reproduce; worst gap ", worst_gap, ".")
}
message(
  "  R6's 24 published estimates reproduce exactly (worst gap ",
  format(worst_gap, digits = 3), ")."
)

# -----------------------------------------------------------------------------

step(6, TOTAL_STEPS, "Writing the cache.")

elapsed_minutes <- as.numeric(
  difftime(Sys.time(), started_at, units = "mins")
)
cache <- list(
  configuration = list(
    experiment_id = EXPERIMENT_ID,
    enhancer_definitions = ENHANCER_DEFINITIONS,
    enhancer_definition_tags = ENHANCER_DEFINITION_TAGS,
    expected_sets = EXPECTED_SETS,
    expected_variant_counts = EXPECTED_VARIANT_COUNTS,
    background_variant_count = nrow(background),
    r_version = R.version.string,
    built_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    elapsed_minutes = elapsed_minutes
  ),
  panel = panel,
  set_metadata = set_metadata,
  annotation_provenance = annotation_provenance,
  state_summary = state_summary,
  definition_coverage = definition_coverage,
  r6_agreement = r6_agreement,
  r6_estimate_reproduction = comparison,
  minimum_overlap = MINIMUM_OVERLAP,
  not_estimable = not_estimable,
  enrichment = enrichment,
  summary_table = summary_table
)
saveRDS(cache, file.path(output_directory, "analysis_cache.rds"))
utils::write.csv(
  data.frame(
    elapsed_minutes = elapsed_minutes,
    built_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  ),
  file.path(output_directory, "runtime_summary.csv"),
  row.names = FALSE
)
message(
  "Done in ", format(round(elapsed_minutes, 1)), " minutes. Cache: ",
  output_directory
)
