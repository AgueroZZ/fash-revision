#!/usr/bin/env Rscript

# Unit tests for temporal_class_helpers.R, on synthetic inputs only. No cached
# result is read here; the cache-dependent checks live in test_reporting.R.

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
  "fash_temporal_class_enhancer_enrichment", "temporal_class_helpers.R"
))

expect <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop("FAILED: ", message, call. = FALSE)
  }
  invisible(TRUE)
}

expect_error <- function(expression, message) {
  outcome <- tryCatch({
    force(expression)
    FALSE
  }, error = function(e) TRUE)
  expect(outcome, message)
}

make_testing_result <- function(lfsr, cfsr, keys) {
  data.frame(
    indices = seq_along(lfsr),
    lfsr = lfsr,
    cfsr = cfsr,
    row.names = keys
  )
}

pair_keys <- c(
  "GENE1_rs1", "GENE1_rs2", "GENE2_rs3", "GENE2_rs4", "GENE3_rs5"
)

# Pair 1 peaks early, pairs 2 and 4 late, pair 3 middle, pair 5 early but with
# a weak call (winning probability 0.4).
testing_results <- list(
  early = make_testing_result(
    lfsr = c(0.10, 0.90, 0.70, 0.95, 0.60),
    cfsr = c(0.10, 0.50, 0.57, 0.66, 0.65),
    keys = pair_keys
  ),
  middle = make_testing_result(
    lfsr = c(0.80, 0.60, 0.20, 0.80, 0.75),
    cfsr = c(0.20, 0.50, 0.53, 0.60, 0.63),
    keys = pair_keys
  ),
  late = make_testing_result(
    lfsr = c(0.95, 0.05, 0.85, 0.02, 0.90),
    cfsr = c(0.02, 0.04, 0.30, 0.24, 0.35),
    keys = pair_keys
  )
)

# ---- derive_window_probabilities -------------------------------------------

probabilities <- derive_window_probabilities(testing_results)
expect(
  identical(dim(probabilities), c(5L, 3L)),
  "derive_window_probabilities returns one row per pair and three columns"
)
expect(
  identical(colnames(probabilities), TEMPORAL_CLASSES),
  "the probability columns are in time order"
)
expect(
  identical(rownames(probabilities), pair_keys),
  "the pair keys survive"
)
expect(
  isTRUE(all.equal(probabilities[, "early"], 1 - testing_results$early$lfsr,
                   check.attributes = FALSE)),
  "the probability is one minus lfsr"
)

expect_error(
  derive_window_probabilities(testing_results[c("early", "middle")]),
  "a missing window is rejected"
)
misaligned <- testing_results
rownames(misaligned$late) <- rev(pair_keys)
expect_error(
  derive_window_probabilities(misaligned),
  "misaligned row names are rejected"
)
out_of_range <- testing_results
out_of_range$middle$lfsr[1L] <- 1.5
expect_error(
  derive_window_probabilities(out_of_range),
  "an lfsr outside [0, 1] is rejected"
)

# ---- assign_temporal_class --------------------------------------------------

classes <- assign_temporal_class(probabilities)
expect(
  identical(as.character(classes$class),
            c("early", "late", "middle", "late", "early")),
  "each pair is assigned to its most probable window"
)
expect(
  isTRUE(all.equal(classes$probability, c(0.90, 0.95, 0.80, 0.98, 0.40))),
  "the reported probability is the winning one"
)
expect(
  identical(levels(classes$class), TEMPORAL_CLASSES),
  "the class factor keeps time order"
)

tied <- matrix(
  c(0.4, 0.4, 0.2), nrow = 1L, dimnames = list("GENE9_rs9", TEMPORAL_CLASSES)
)
expect(
  identical(as.character(assign_temporal_class(tied)$class), "early"),
  "ties resolve to the earliest window"
)

# ---- label_discovered_pairs -------------------------------------------------

pair_table <- data.frame(
  pair_key = rev(pair_keys),
  gene_id = rev(sub("_.*$", "", pair_keys)),
  variant_id = rev(sub("^[^_]+_", "", pair_keys)),
  score = c(0.05, 0.04, 0.01, 0.03, 0.02),
  stringsAsFactors = FALSE
)
labelled <- label_discovered_pairs(pair_table, classes)
expect(
  identical(labelled$pair_key,
            c("GENE2_rs3", "GENE1_rs1", "GENE1_rs2", "GENE2_rs4", "GENE3_rs5")),
  "the labelled pairs come back ordered by lfdr"
)
expect(
  identical(as.character(labelled$class),
            c("middle", "early", "late", "late", "early")),
  "the class travels with its own pair key, not with row position"
)
expect_error(
  label_discovered_pairs(pair_table[-1L, , drop = FALSE], classes),
  "a discovery set that is not the classified set is rejected"
)

# ---- build_temporal_variant_sets --------------------------------------------

sets <- build_temporal_variant_sets(labelled)
expect(
  identical(names(sets), c(
    "early_all", "early_lead", "middle_all", "middle_lead",
    "late_all", "late_lead"
  )),
  "both strategies are produced for all three classes"
)
expect(
  identical(sort(sets$late_all), c("rs2", "rs4")),
  "the all-variant set holds every variant of the class"
)
# GENE1_rs2 (late, lfdr 0.03) and GENE2_rs4 (late, lfdr 0.04) are the only late
# pairs of their genes, so both survive lead selection.
expect(
  identical(sort(sets$late_lead), c("rs2", "rs4")),
  "the lead set keeps the lowest-lfdr pair of each gene within the class"
)
# GENE1 contributes rs1 to early and rs2 to late: one lead per gene per class.
expect(
  identical(sort(sets$early_lead), c("rs1", "rs5")),
  "a gene splitting across classes contributes a lead to each of them"
)
expect(
  identical(sets$middle_all, "rs3") && identical(sets$middle_lead, "rs3"),
  "a single-variant class agrees across strategies"
)

confident <- build_temporal_variant_sets(labelled, minimum_probability = 0.7)
expect(
  identical(confident$early_all, "rs1"),
  "the probability floor drops the weakly assigned pair"
)
expect(
  identical(sort(confident$late_all), c("rs2", "rs4")),
  "the probability floor leaves confident pairs alone"
)
expect_error(
  build_temporal_variant_sets(labelled, minimum_probability = 2),
  "an out-of-range probability floor is rejected"
)

# ---- summaries --------------------------------------------------------------

summary_table <- summarise_temporal_classes(labelled)
expect(
  identical(summary_table$class, TEMPORAL_CLASSES),
  "the summary keeps one row per class in time order"
)
expect(
  identical(summary_table$pairs, c(2L, 1L, 2L)),
  "the summary counts pairs per class"
)
expect(
  identical(summary_table$genes, c(2L, 1L, 2L)),
  "the summary counts unique genes per class"
)
expect(
  isTRUE(all.equal(summary_table$median_probability, c(0.65, 0.80, 0.965))),
  "the summary reports the median winning probability"
)

confident_counts <- count_confident_calls(testing_results, alpha = 0.05)
expect(
  identical(confident_counts$pairs, c(0L, 0L, 2L)),
  "the cfsr counts use the cumulative false sign rate"
)
# The two confident late calls are two variants of GENE1, so the gene count is
# one: the pair key, not the row, decides.
expect(
  identical(confident_counts$variants, c(0L, 0L, 2L)) &&
    identical(confident_counts$genes, c(0L, 0L, 1L)),
  "the cfsr counts resolve variants and genes from the pair key"
)

message("temporal_class_helpers.R: all tests passed.")
