#!/usr/bin/env Rscript

# =============================================================================
# Helpers for the internal Roadmap enhancer exploration.
#
# The reviewer-facing R6 page shows three Roadmap epigenomes (E020, E013, E095).
# This experiment widens the panel to thirteen, in four pre-specified biological
# groups, and reports two enhancer-state definitions side by side:
#
#   R6       6_EnhG + 7_Enh                             (what R6 uses today)
#   Strober  6_EnhG + 7_Enh + 12_EnhBiv + 11_BivFlnk    (Strober et al. 2019)
#
# Nothing here changes the estimator. Fold enrichment, the tested-variant
# background and the leave-one-autosome-out jackknife all come from
# fash_strober_enhancer_comparison/enrichment_api.R, unmodified. This file only
# builds a wider annotation table for that API to consume, and the bookkeeping
# needed to describe it.
# =============================================================================

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

# normalize_autosome() and annotate_variant_overlaps() are shared with the other
# enrichment pages, so they are sourced rather than duplicated.
source(file.path(
  WORKFLOWR_ROOT, "code", "revision_simulations", "internal",
  "variant_annotation_enrichment", "variant_annotation_enrichment_helpers.R"
))

# The column-binding step lives in variant_annotation_enrichment's *run* script
# rather than in its helper file, and sourcing a run script would execute an
# unrelated analysis. It is four lines, so it is restated here with the same two
# guards: rows must already be aligned, and names must not collide.
append_enhancer_columns <- function(annotation_matrix, overlap_table) {
  if (!identical(annotation_matrix$variant_id, overlap_table$variant_id)) {
    stop("Annotation overlap rows are not aligned to the variant universe.")
  }
  new_columns <- setdiff(names(overlap_table), c("variant_id", "chromosome"))
  if (any(new_columns %in% names(annotation_matrix))) {
    stop("Annotation names are duplicated.")
  }
  for (column in new_columns) {
    annotation_matrix[[column]] <- as.logical(overlap_table[[column]])
  }
  annotation_matrix
}

EXPERIMENT_ID <- "roadmap_enhancer_exploration"

roadmap_code_directory <- function() {
  file.path(
    WORKFLOWR_ROOT, "code", "revision_simulations", "internal", EXPERIMENT_ID
  )
}

roadmap_data_directory <- function() {
  file.path(
    WORKFLOWR_ROOT, "data", "revision_simulations", "internal", EXPERIMENT_ID
  )
}

roadmap_output_directory <- function() {
  file.path(
    WORKFLOWR_ROOT, "output", "revision_simulations", "internal", EXPERIMENT_ID
  )
}

# -----------------------------------------------------------------------------
# The two enhancer-state definitions, spelled out
# -----------------------------------------------------------------------------

# Roadmap 15-state core-marks model, mnemonics exactly as they appear in column
# 4 of <EID>_15_coreMarks_mnemonics.bed.gz.
ROADMAP_15_STATE_MNEMONICS <- c(
  "1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh",
  "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC",
  "14_ReprPCWk", "15_Quies"
)

# Definition 1: what analysis/revision_enhancer_enrichment.rmd uses today. It
# comes from collapse_roadmap_state() in
# variant_annotation_enrichment/variant_annotation_enrichment_helpers.R, which
# maps states 6 and 7 to the single category "Enhancer".
R6_ENHANCER_STATES <- c("6_EnhG", "7_Enh")

# Definition 2: verified from Strober et al. (2019) source, not inferred from
# the paper text, which does not enumerate states.
#
#   BennyStrobes/ipsc_cardiomyocyte_differentiation@master
#   dynamic_eqtl_calling/perform_tissue_specific_chrom_hmm_enrichment_analysis.py
#   get_valid_markers(), lines 236-254:
#
#       elif marker_type == 'enhancer':
#           valid_markers['7_Enh'] = 1
#           valid_markers['6_EnhG'] = 1
#           valid_markers['12_EnhBiv'] = 1
#           valid_markers['11_BivFlnk'] = 1
#
# Note that 11_BivFlnk is also in their 'promotor' set, so their marker classes
# are not disjoint. We reproduce their enhancer set as written.
STROBER_ENHANCER_STATES <- c("6_EnhG", "7_Enh", "12_EnhBiv", "11_BivFlnk")

ENHANCER_DEFINITIONS <- list(
  R6 = R6_ENHANCER_STATES,
  Strober = STROBER_ENHANCER_STATES
)

# Short tags that go into annotation column names, so a column name alone says
# which state set produced it.
# `EpiCompare` is the three-state definition the reviewer-facing R6 page adopted
# on 2026-08-26 (He & Wang 2017: ChromHMM states 6, 7, 12). It is a naming entry
# only -- deliberately absent from ENHANCER_DEFINITIONS, so the default state set
# this experiment builds is unchanged and its cache and page are unaffected. The
# R6 run script passes its own definition list explicitly.
ENHANCER_DEFINITION_TAGS <- c(
  R6 = "R6 EnhG+Enh",
  Strober = "Strober EnhG+Enh+EnhBiv+BivFlnk",
  EpiCompare = "EpiCompare EnhG+Enh+EnhBiv"
)

#' The annotation column name for one epigenome under one definition.
#'
#' Deliberately verbose: these strings are the keys `compute_enrichment()`
#' reports back, and they end up in the page's summary table.
enhancer_annotation_name <- function(epigenome_id, std_name, definition) {
  if (!all(definition %in% names(ENHANCER_DEFINITION_TAGS))) {
    stop("Unknown enhancer definition: ", paste(definition, collapse = ", "))
  }
  paste0(
    "Roadmap ", epigenome_id, " ", std_name, ": Enhancer [",
    unname(ENHANCER_DEFINITION_TAGS[definition]), "]"
  )
}

# -----------------------------------------------------------------------------
# The epigenome panel
# -----------------------------------------------------------------------------

ROADMAP_CHROMHMM_BASE_URL <- paste0(
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/",
  "ChmmModels/coreMarks/jointModel/final/"
)

roadmap_segmentation_filename <- function(epigenome_id) {
  paste0(epigenome_id, "_15_coreMarks_mnemonics.bed.gz")
}

roadmap_segmentation_url <- function(epigenome_id) {
  paste0(ROADMAP_CHROMHMM_BASE_URL, roadmap_segmentation_filename(epigenome_id))
}

#' The pre-specified panel of thirteen epigenomes, with verified Roadmap labels.
#'
#' `roadmap_std_name` and `roadmap_mnemonic` were taken from Roadmap's own
#' metadata table, `EID_metadata.tab`, not from any secondary description.
load_epigenome_panel <- function(path = file.path(
                                   roadmap_code_directory(),
                                   "roadmap_epigenome_manifest.csv"
                                 )) {
  if (!file.exists(path)) {
    stop("Missing epigenome manifest: ", path)
  }
  panel <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(
    "epigenome_id", "biological_group", "group_order", "epigenome_order",
    "roadmap_std_name", "roadmap_mnemonic", "roadmap_group", "roadmap_anatomy",
    "display_label", "in_r6_page"
  )
  if (!all(required %in% names(panel)) || !nrow(panel) ||
      anyDuplicated(panel$epigenome_id) ||
      anyDuplicated(panel$display_label) ||
      !all(grepl("^E[0-9]{3}$", panel$epigenome_id))) {
    stop("The epigenome manifest is invalid.")
  }
  panel$in_r6_page <- as.logical(panel$in_r6_page)
  if (anyNA(panel$in_r6_page) || sum(panel$in_r6_page) != 3L) {
    stop("Exactly three panel members should be flagged as used by R6.")
  }
  panel$local_filename <- roadmap_segmentation_filename(panel$epigenome_id)
  panel$url <- roadmap_segmentation_url(panel$epigenome_id)
  panel <- panel[order(panel$group_order, panel$epigenome_order), , drop = FALSE]
  rownames(panel) <- NULL
  panel
}

#' Display order for the biological groups, and which one is the aside.
#'
#' E065 aorta sits in Roadmap's `Heart` GROUP but its ANATOMY is `VASCULAR`. It
#' is kept in its own group so that no statement about "heart" epigenomes is
#' silently carried by a vascular tissue.
biological_group_levels <- function(panel = load_epigenome_panel()) {
  unique(panel$biological_group[order(panel$group_order)])
}

# -----------------------------------------------------------------------------
# Download and verify the segmentations
# -----------------------------------------------------------------------------

validate_gzip_file <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Segmentation file is missing or empty: ", path)
  }
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  magic <- readBin(connection, what = "raw", n = 2L)
  if (!identical(as.integer(magic), c(31L, 139L))) {
    stop("Segmentation file is not gzip: ", path)
  }
  invisible(TRUE)
}

#' Fetch the thirteen segmentation files, skipping any already present and valid.
#'
#' Returns a provenance table: URL, byte size, md5, sha256, retrieval time. The
#' page records it so that "which annotation files were used" has an exact
#' answer.
download_roadmap_segmentations <- function(panel = load_epigenome_panel(),
                                           target_directory =
                                             roadmap_data_directory(),
                                           quiet = FALSE) {
  dir.create(target_directory, recursive = TRUE, showWarnings = FALSE)
  rows <- lapply(seq_len(nrow(panel)), function(index) {
    epigenome_id <- panel$epigenome_id[index]
    target_path <- file.path(target_directory, panel$local_filename[index])
    already_valid <- FALSE
    if (file.exists(target_path)) {
      already_valid <- !inherits(
        try(validate_gzip_file(target_path), silent = TRUE), "try-error"
      )
    }
    if (!already_valid) {
      if (!quiet) message("Downloading Roadmap segmentation: ", epigenome_id)
      temporary_path <- tempfile(
        pattern = paste0(epigenome_id, "-"),
        tmpdir = target_directory,
        fileext = ".download"
      )
      on.exit(unlink(temporary_path), add = TRUE)
      status <- utils::download.file(
        url = panel$url[index],
        destfile = temporary_path,
        method = "libcurl",
        mode = "wb",
        quiet = quiet
      )
      if (!identical(status, 0L)) {
        stop("Download failed for ", epigenome_id, ".")
      }
      validate_gzip_file(temporary_path)
      if (file.exists(target_path)) unlink(target_path)
      if (!file.rename(temporary_path, target_path)) {
        stop("Could not install segmentation: ", target_path)
      }
    }
    validate_gzip_file(target_path)
    data.frame(
      epigenome_id = epigenome_id,
      roadmap_std_name = panel$roadmap_std_name[index],
      biological_group = panel$biological_group[index],
      url = panel$url[index],
      local_path = target_path,
      byte_size = unname(file.info(target_path)$size),
      md5 = unname(tools::md5sum(target_path)),
      sha256 = if (requireNamespace("digest", quietly = TRUE)) {
        digest::digest(file = target_path, algo = "sha256", serialize = FALSE)
      } else {
        NA_character_
      },
      retrieved_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# -----------------------------------------------------------------------------
# Segmentation -> enhancer intervals, at state resolution
# -----------------------------------------------------------------------------

#' Read one segmentation and emit enhancer intervals for both definitions.
#'
#' The R6 annotation matrix collapses states 6 and 7 into one "Enhancer"
#' category before any indicator is built, which makes a four-state definition
#' impossible to recover from it. This reads the mnemonics directly instead, so
#' both definitions come out of the same file in one pass. States 6 and 7 belong
#' to both definitions and therefore contribute one interval row per definition.
#'
#' @return data frame with `chromosome`, `start`, `end`, `annotation`, in the
#'   1-based inclusive convention `annotate_variant_overlaps()` expects.
read_roadmap_enhancer_intervals <- function(path,
                                            epigenome_id,
                                            std_name,
                                            definitions = ENHANCER_DEFINITIONS) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The data.table package is required to read Roadmap segmentations.")
  }
  segmentation <- data.table::fread(
    path,
    sep = "\t",
    header = FALSE,
    col.names = c("chromosome", "start0", "end", "state"),
    showProgress = FALSE,
    data.table = FALSE
  )
  unexpected <- setdiff(unique(segmentation$state), ROADMAP_15_STATE_MNEMONICS)
  if (length(unexpected)) {
    stop(
      "Unexpected ChromHMM state(s) in ", basename(path), ": ",
      paste(unexpected, collapse = ", ")
    )
  }
  wanted <- unique(unlist(definitions, use.names = FALSE))
  segmentation <- segmentation[segmentation$state %in% wanted, , drop = FALSE]
  if (!nrow(segmentation)) {
    stop("No enhancer-state intervals found in ", basename(path), ".")
  }
  # BED is 0-based half-open; the shared overlap helper takes 1-based inclusive.
  start <- as.integer(segmentation$start0) + 1L
  end <- as.integer(segmentation$end)
  chromosome <- normalize_autosome(segmentation$chromosome)
  intervals <- lapply(names(definitions), function(definition) {
    keep <- segmentation$state %in% definitions[[definition]]
    if (!any(keep)) {
      stop("Definition ", definition, " selected no intervals in ", basename(path), ".")
    }
    data.frame(
      chromosome = chromosome[keep],
      start = start[keep],
      end = end[keep],
      annotation = enhancer_annotation_name(epigenome_id, std_name, definition),
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, intervals)
  output[!is.na(output$chromosome), , drop = FALSE]
}

#' Per-state interval counts and covered base pairs, for the record.
summarise_roadmap_states <- function(path, epigenome_id) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The data.table package is required to read Roadmap segmentations.")
  }
  segmentation <- data.table::fread(
    path,
    sep = "\t",
    header = FALSE,
    col.names = c("chromosome", "start0", "end", "state"),
    showProgress = FALSE,
    data.table = FALSE
  )
  autosomal <- !is.na(normalize_autosome(segmentation$chromosome))
  segmentation <- segmentation[autosomal, , drop = FALSE]
  width <- as.numeric(segmentation$end) - as.numeric(segmentation$start0)
  states <- sort(unique(segmentation$state))
  data.frame(
    epigenome_id = epigenome_id,
    state = states,
    in_r6_definition = states %in% R6_ENHANCER_STATES,
    in_strober_definition = states %in% STROBER_ENHANCER_STATES,
    interval_count = as.integer(tabulate(match(segmentation$state, states),
                                         length(states))),
    autosomal_covered_bases = vapply(states, function(s) {
      sum(width[segmentation$state == s])
    }, numeric(1)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Build the variant x (epigenome, definition) indicator table.
#'
#' Same shape as the annotation table `load_annotation_groups()` returns —
#' `variant_id`, `chromosome`, then one logical column per annotation — so
#' `compute_enrichment()` accepts it directly and the estimator is untouched.
build_enhancer_annotation_matrix <- function(variants,
                                             panel = load_epigenome_panel(),
                                             data_directory =
                                               roadmap_data_directory(),
                                             definitions = ENHANCER_DEFINITIONS,
                                             verbose = TRUE) {
  required <- c("variant_id", "chromosome", "position")
  if (!all(required %in% names(variants)) || !nrow(variants) ||
      anyDuplicated(variants$variant_id)) {
    stop("Variant coordinates have an unexpected structure.")
  }
  coordinates <- variants[, required, drop = FALSE]
  annotation_matrix <- coordinates[, c("variant_id", "chromosome"), drop = FALSE]
  for (index in seq_len(nrow(panel))) {
    epigenome_id <- panel$epigenome_id[index]
    if (verbose) message("  annotating ", epigenome_id, " ...")
    intervals <- read_roadmap_enhancer_intervals(
      file.path(data_directory, panel$local_filename[index]),
      epigenome_id,
      panel$roadmap_std_name[index],
      definitions = definitions
    )
    overlaps <- annotate_variant_overlaps(coordinates, intervals)
    annotation_matrix <- append_enhancer_columns(annotation_matrix, overlaps)
    rm(intervals, overlaps)
    invisible(gc())
  }
  expected_columns <- as.vector(t(outer(
    seq_len(nrow(panel)), names(definitions),
    Vectorize(function(i, d) {
      enhancer_annotation_name(
        panel$epigenome_id[i], panel$roadmap_std_name[i], d
      )
    })
  )))
  missing <- setdiff(expected_columns, names(annotation_matrix))
  if (length(missing)) {
    stop("Annotation matrix is missing columns: ",
         paste(missing, collapse = ", "))
  }
  annotation_matrix[, c("variant_id", "chromosome", expected_columns),
                    drop = FALSE]
}

# -----------------------------------------------------------------------------
# Labelling the estimates
# -----------------------------------------------------------------------------

#' Split an annotation column name back into epigenome and definition.
#'
#' `compute_enrichment()` returns the column name verbatim; this recovers the
#' fields the figures facet and colour by, and fails loudly on anything it does
#' not recognise rather than dropping rows.
decorate_enrichment <- function(enrichment,
                                panel = load_epigenome_panel(),
                                definitions = ENHANCER_DEFINITIONS) {
  key <- do.call(rbind, lapply(seq_len(nrow(panel)), function(index) {
    do.call(rbind, lapply(names(definitions), function(definition) {
      data.frame(
        annotation = enhancer_annotation_name(
          panel$epigenome_id[index], panel$roadmap_std_name[index], definition
        ),
        epigenome_id = panel$epigenome_id[index],
        roadmap_std_name = panel$roadmap_std_name[index],
        roadmap_mnemonic = panel$roadmap_mnemonic[index],
        biological_group = panel$biological_group[index],
        group_order = panel$group_order[index],
        epigenome_order = panel$epigenome_order[index],
        display_label = panel$display_label[index],
        in_r6_page = panel$in_r6_page[index],
        enhancer_definition = definition,
        stringsAsFactors = FALSE
      )
    }))
  }))
  unknown <- setdiff(enrichment$annotation, key$annotation)
  if (length(unknown)) {
    stop("Unrecognised annotation(s): ", paste(unknown, collapse = ", "))
  }
  merged <- merge(enrichment, key, by = "annotation", all.x = TRUE, sort = FALSE)
  if (nrow(merged) != nrow(enrichment) || anyNA(merged$epigenome_id)) {
    stop("Failed to label the enrichment estimates.")
  }
  merged[order(
    merged$enhancer_definition, merged$variant_set,
    merged$group_order, merged$epigenome_order
  ), , drop = FALSE]
}

#' The compact per-(epigenome, discovery set) table the page ends with.
build_summary_table <- function(enrichment) {
  required <- c(
    "epigenome_id", "roadmap_std_name", "biological_group",
    "enhancer_definition", "variant_set", "selected_total", "selected_overlap",
    "control_total", "control_overlap", "fold_enrichment", "ci_lower_fold",
    "ci_upper_fold"
  )
  missing <- setdiff(required, names(enrichment))
  if (length(missing)) {
    stop("Cannot build the summary table; missing: ",
         paste(missing, collapse = ", "))
  }
  data.frame(
    epigenome_id = enrichment$epigenome_id,
    roadmap_name = enrichment$roadmap_std_name,
    biological_group = enrichment$biological_group,
    enhancer_definition = enrichment$enhancer_definition,
    discovery_set = enrichment$variant_set,
    n_variants = enrichment$selected_total,
    enhancer_overlap_proportion =
      enrichment$selected_overlap / enrichment$selected_total,
    background_overlap_proportion =
      enrichment$control_overlap / enrichment$control_total,
    fold_enrichment = enrichment$fold_enrichment,
    jackknife_ci_lower = enrichment$ci_lower_fold,
    jackknife_ci_upper = enrichment$ci_upper_fold,
    stringsAsFactors = FALSE
  )
}
