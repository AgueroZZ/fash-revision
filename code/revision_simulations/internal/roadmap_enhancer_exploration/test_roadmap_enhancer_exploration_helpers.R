#!/usr/bin/env Rscript

# Unit tests for the internal Roadmap enhancer exploration helpers.
#
# The parts worth testing are the parts that could silently change a number:
# the two enhancer-state definitions, the BED-to-1-based coordinate shift, the
# state filter, and the round trip from annotation column name back to
# (epigenome, definition). Everything downstream of that is enrichment_api.R,
# which has its own check against the published R6 estimates.
#
#   Rscript --vanilla test_roadmap_enhancer_exploration_helpers.R

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
  "roadmap_enhancer_exploration", "roadmap_enhancer_exploration_helpers.R"
))

failures <- 0L
check <- function(description, condition) {
  passed <- isTRUE(condition)
  if (!passed) failures <<- failures + 1L
  message(if (passed) "PASS  " else "FAIL  ", description)
  invisible(passed)
}

# ---- The definitions themselves --------------------------------------------

check(
  "R6 enhancer definition is exactly EnhG + Enh",
  identical(sort(R6_ENHANCER_STATES), sort(c("6_EnhG", "7_Enh")))
)

check(
  "Strober enhancer definition is EnhG + Enh + EnhBiv + BivFlnk",
  identical(
    sort(STROBER_ENHANCER_STATES),
    sort(c("6_EnhG", "7_Enh", "12_EnhBiv", "11_BivFlnk"))
  )
)

check(
  "the Strober definition strictly contains the R6 definition",
  all(R6_ENHANCER_STATES %in% STROBER_ENHANCER_STATES) &&
    length(setdiff(STROBER_ENHANCER_STATES, R6_ENHANCER_STATES)) == 2L
)

check(
  "every definition state is a real 15-state mnemonic",
  all(unlist(ENHANCER_DEFINITIONS) %in% ROADMAP_15_STATE_MNEMONICS)
)

check(
  "annotation names carry the epigenome and the definition",
  identical(
    enhancer_annotation_name("E095", "Left Ventricle", "R6"),
    "Roadmap E095 Left Ventricle: Enhancer [R6 EnhG+Enh]"
  ) &&
    enhancer_annotation_name("E095", "Left Ventricle", "R6") !=
      enhancer_annotation_name("E095", "Left Ventricle", "Strober")
)

check(
  "an unknown definition is refused",
  inherits(
    try(enhancer_annotation_name("E095", "Left Ventricle", "loose"),
        silent = TRUE),
    "try-error"
  )
)

# ---- The panel --------------------------------------------------------------

panel <- load_epigenome_panel()

check("the panel has thirteen epigenomes", nrow(panel) == 13L)

check(
  "the panel contains the three epigenomes R6 shows",
  identical(sort(panel$epigenome_id[panel$in_r6_page]),
            c("E013", "E020", "E095"))
)

check(
  "the panel matches Strober's iPSC and heart cell-line lists",
  identical(
    sort(panel$epigenome_id[panel$biological_group == "iPSC"]),
    c("E018", "E019", "E020", "E021", "E022")
  ) &&
    identical(
      sort(panel$epigenome_id[panel$biological_group %in%
                                c("Heart", "Vascular")]),
      c("E065", "E083", "E095", "E104", "E105")
    )
)

check(
  "the aorta is kept out of the heart group",
  panel$biological_group[panel$epigenome_id == "E065"] == "Vascular"
)

check(
  "groups are ordered iPSC, germ layer, heart, vascular",
  identical(biological_group_levels(panel),
            c("iPSC", "Germ layer", "Heart", "Vascular"))
)

check(
  "segmentation URLs point at the 15-state core-marks mnemonics",
  all(grepl("_15_coreMarks_mnemonics\\.bed\\.gz$", panel$url)) &&
    all(startsWith(panel$url, ROADMAP_CHROMHMM_BASE_URL))
)

# ---- Segmentation reading ---------------------------------------------------

# A miniature segmentation covering one state from each relevant class: an
# R6-and-Strober state (7_Enh), two Strober-only states (12_EnhBiv, 11_BivFlnk),
# and two states in neither (1_TssA, 15_Quies).
write_segmentation <- function(lines) {
  path <- tempfile(fileext = ".bed.gz")
  connection <- gzfile(path, open = "wt")
  writeLines(lines, connection)
  close(connection)
  path
}

fixture <- write_segmentation(c(
  "chr1\t1000\t1200\t7_Enh",
  "chr1\t2000\t2100\t12_EnhBiv",
  "chr1\t3000\t3050\t11_BivFlnk",
  "chr1\t4000\t4100\t1_TssA",
  "chr2\t5000\t5100\t6_EnhG",
  "chr2\t6000\t6100\t15_Quies",
  "chrX\t7000\t7100\t7_Enh"
))

intervals <- read_roadmap_enhancer_intervals(fixture, "E999", "Test Tissue")

r6_name <- enhancer_annotation_name("E999", "Test Tissue", "R6")
strober_name <- enhancer_annotation_name("E999", "Test Tissue", "Strober")

check(
  "only the two definition columns are produced",
  identical(sort(unique(intervals$annotation)), sort(c(r6_name, strober_name)))
)

check(
  "the R6 definition keeps only the two Enh states",
  sum(intervals$annotation == r6_name) == 2L
)

check(
  "the Strober definition adds EnhBiv and BivFlnk",
  sum(intervals$annotation == strober_name) == 4L
)

check(
  "non-enhancer states are dropped under both definitions",
  !any(intervals$start == 4001L) && !any(intervals$start == 6001L)
)

check(
  "the sex chromosome interval is dropped",
  !any(intervals$chromosome %in% c("X", "23")) &&
    identical(sort(unique(intervals$chromosome)), c("1", "2"))
)

check(
  "BED half-open coordinates become 1-based inclusive",
  {
    row <- intervals[intervals$annotation == r6_name &
                       intervals$chromosome == "1", , drop = FALSE]
    nrow(row) == 1L && row$start == 1001L && row$end == 1200L
  }
)

check(
  "an unrecognised ChromHMM state is refused rather than ignored",
  {
    bad <- write_segmentation("chr1\t10\t20\t16_Invented")
    inherits(
      try(read_roadmap_enhancer_intervals(bad, "E999", "Test Tissue"),
          silent = TRUE),
      "try-error"
    )
  }
)

# ---- Variant indicators -----------------------------------------------------

# One variant just inside the 7_Enh interval, one just outside it, one inside
# the Strober-only 12_EnhBiv interval, one in 1_TssA (neither definition), and
# one inside the chr2 6_EnhG interval.
variants <- data.frame(
  variant_id = c("inside_enh", "outside_enh", "inside_enhbiv", "inside_tssa",
                 "inside_enhg"),
  chromosome = c("1", "1", "1", "1", "2"),
  position = c(1001L, 1201L, 2050L, 4050L, 5050L),
  stringsAsFactors = FALSE
)
overlaps <- annotate_variant_overlaps(variants, intervals)

check(
  "the R6 indicator fires only inside EnhG/Enh",
  identical(
    as.logical(overlaps[[r6_name]]),
    c(TRUE, FALSE, FALSE, FALSE, TRUE)
  )
)

check(
  "the Strober indicator additionally fires inside EnhBiv",
  identical(
    as.logical(overlaps[[strober_name]]),
    c(TRUE, FALSE, TRUE, FALSE, TRUE)
  )
)

check(
  "the Strober indicator is never smaller than the R6 indicator",
  all(as.logical(overlaps[[strober_name]]) >= as.logical(overlaps[[r6_name]]))
)

# ---- Labelling and summarising ----------------------------------------------

fake_enrichment <- data.frame(
  variant_set = c("current_lead", "current_lead"),
  annotation = c(
    enhancer_annotation_name("E095", "Left Ventricle", "R6"),
    enhancer_annotation_name("E013", "hESC Derived CD56+ Mesoderm Cultured Cells",
                             "Strober")
  ),
  selected_overlap = c(200L, 300L),
  selected_total = c(1000L, 1000L),
  control_overlap = c(150000L, 250000L),
  control_total = c(1000000L, 1000000L),
  fold_enrichment = c(1.33, 1.2),
  ci_lower_fold = c(1.1, 1.0),
  ci_upper_fold = c(1.6, 1.44),
  stringsAsFactors = FALSE
)
labelled <- decorate_enrichment(fake_enrichment)

check(
  "annotation names round-trip back to epigenome and definition",
  identical(labelled$epigenome_id[labelled$enhancer_definition == "R6"], "E095") &&
    identical(
      labelled$epigenome_id[labelled$enhancer_definition == "Strober"], "E013"
    ) &&
    identical(
      labelled$biological_group[labelled$enhancer_definition == "Strober"],
      "Germ layer"
    )
)

check(
  "an annotation outside the panel is refused rather than dropped",
  inherits(
    try(
      decorate_enrichment(
        transform(fake_enrichment[1L, , drop = FALSE], annotation = "made up")
      ),
      silent = TRUE
    ),
    "try-error"
  )
)

summary_table <- build_summary_table(labelled)

check(
  "the summary table reports the overlap proportions it claims to",
  all(abs(
    summary_table$enhancer_overlap_proportion -
      labelled$selected_overlap / labelled$selected_total
  ) < 1e-12) &&
    all(abs(
      summary_table$background_overlap_proportion -
        labelled$control_overlap / labelled$control_total
    ) < 1e-12)
)

check(
  "the summary table carries every requested column",
  all(c(
    "epigenome_id", "roadmap_name", "biological_group", "enhancer_definition",
    "discovery_set", "n_variants", "enhancer_overlap_proportion",
    "background_overlap_proportion", "fold_enrichment", "jackknife_ci_lower",
    "jackknife_ci_upper"
  ) %in% names(summary_table))
)

message("")
if (failures > 0L) {
  stop(failures, " test(s) failed.")
}
message("All tests passed.")
