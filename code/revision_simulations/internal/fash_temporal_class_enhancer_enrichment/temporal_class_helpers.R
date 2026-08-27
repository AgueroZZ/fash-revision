# Temporal classification of the FASH dynamic-eQTL discoveries.
#
# The manuscript classifies discovered gene-variant pairs by where the peak of
# |beta(t)| falls in the 15-day differentiation:
#
#   early   sup_{t<=3}|beta|      > sup_{t>3}|beta|
#   middle  sup_{3<t<12}|beta|    > sup_{t<=3 or t>=12}|beta|
#   late    sup_{t>=12}|beta|     > sup_{t<12}|beta|
#
# `testing_functional()` returns, per pair, `lfsr` = posterior probability that
# the functional is <= 0, i.e. that the peak is *not* in that window. This file
# turns the three cached lfsr vectors into one class per pair and into the
# variant sets the enrichment API consumes. It computes nothing new about the
# fit; it only reshapes results that already exist.

TEMPORAL_CLASSES <- c("early", "middle", "late")

TEMPORAL_CLASS_LABELS <- c(
  early = "Early (peak on days 0-3)",
  middle = "Middle (peak at 3 < t < 12)",
  late = "Late (peak on days 12-15)"
)

#' Posterior probability that the peak lies in each timing window.
#'
#' @param testing_results named list with elements `early`, `middle`, `late`,
#'   each a `testing_functional()` data frame carrying an `lfsr` column and
#'   gene_variant row names.
#' @return numeric matrix, one row per pair, three columns, row names carried
#'   through from the inputs.
derive_window_probabilities <- function(testing_results) {
  if (!is.list(testing_results) ||
      !identical(sort(names(testing_results)), sort(TEMPORAL_CLASSES))) {
    stop("testing_results must be a named list of early, middle, and late.")
  }
  reference_names <- rownames(testing_results[[TEMPORAL_CLASSES[1L]]])
  if (is.null(reference_names) || !length(reference_names) ||
      anyDuplicated(reference_names)) {
    stop("The testing results must carry unique gene_variant row names.")
  }
  columns <- lapply(TEMPORAL_CLASSES, function(class_name) {
    result <- testing_results[[class_name]]
    if (!is.data.frame(result) || !"lfsr" %in% names(result) ||
        !identical(rownames(result), reference_names)) {
      stop(
        "The ", class_name,
        " testing result is not aligned with the others or has no lfsr column."
      )
    }
    lfsr <- as.numeric(result$lfsr)
    if (any(!is.finite(lfsr)) || any(lfsr < 0 | lfsr > 1)) {
      stop("The ", class_name, " lfsr column is not a probability.")
    }
    1 - lfsr
  })
  probabilities <- do.call(cbind, columns)
  dimnames(probabilities) <- list(reference_names, TEMPORAL_CLASSES)
  probabilities
}

#' Assign every pair to the window with the largest posterior probability.
#'
#' The three windows partition the evaluation grid. Exact-null posterior draws
#' belong to no timing class, so their probabilities sum to the posterior
#' non-null probability rather than one; independently sampled cached results
#' also introduce Monte Carlo error.
#'
#' @return data frame with `pair_key`, `class` (factor over the three classes
#'   in time order), `probability` (the winning probability), and the three
#'   per-window probabilities.
assign_temporal_class <- function(window_probabilities) {
  if (!is.matrix(window_probabilities) ||
      !identical(colnames(window_probabilities), TEMPORAL_CLASSES) ||
      !nrow(window_probabilities)) {
    stop("window_probabilities must be a pairs-by-three matrix.")
  }
  winner <- max.col(window_probabilities, ties.method = "first")
  data.frame(
    pair_key = rownames(window_probabilities),
    class = factor(TEMPORAL_CLASSES[winner], levels = TEMPORAL_CLASSES),
    probability = window_probabilities[
      cbind(seq_len(nrow(window_probabilities)), winner)
    ],
    early_probability = window_probabilities[, "early"],
    middle_probability = window_probabilities[, "middle"],
    late_probability = window_probabilities[, "late"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Attach the temporal class to the R6 discovered-pair table.
#'
#' @param pair_table the `current_all` table from
#'   `discovery_pair_tables.rds`: `pair_key`, `gene_id`, `variant_id`, `score`
#'   (the cumulative lfdr the pair was ranked by).
#' @param class_table output of `assign_temporal_class()`.
#' @return `pair_table` ordered by lfdr, with `class` and `probability` added.
#'   Requires an exact one-to-one match between the two pair-key sets, which is
#'   what ties the cached classification to the R6 discovery set.
label_discovered_pairs <- function(pair_table, class_table) {
  required <- c("pair_key", "gene_id", "variant_id", "score")
  if (!all(required %in% names(pair_table)) ||
      anyDuplicated(pair_table$pair_key) ||
      !setequal(pair_table$pair_key, class_table$pair_key)) {
    stop(
      "The discovered pairs and the classified pairs are not the same set; ",
      "the cached classification does not belong to this discovery set."
    )
  }
  matched <- match(pair_table$pair_key, class_table$pair_key)
  labelled <- pair_table[, required, drop = FALSE]
  labelled$class <- class_table$class[matched]
  labelled$probability <- class_table$probability[matched]
  labelled <- labelled[
    order(labelled$score, labelled$variant_id, labelled$pair_key,
          method = "radix"),
    ,
    drop = FALSE
  ]
  row.names(labelled) <- NULL
  labelled
}

#' The per-class variant sets, in both selection strategies.
#'
#' Lead selection happens *within* a class: for each class, keep the
#' lowest-lfdr pair of each gene. A gene whose variants split across classes
#' therefore contributes one lead to each class it appears in.
#'
#' @param labelled_pairs output of `label_discovered_pairs()`, already ordered
#'   by lfdr.
#' @param minimum_probability keep only pairs whose winning probability is at
#'   least this value.
#' @return named list of unique variant-ID vectors,
#'   `<class>_all` and `<class>_lead` for each class in time order.
build_temporal_variant_sets <- function(labelled_pairs,
                                        minimum_probability = 0) {
  if (!all(c("pair_key", "gene_id", "variant_id", "class", "probability") %in%
             names(labelled_pairs))) {
    stop("labelled_pairs is missing required columns.")
  }
  if (length(minimum_probability) != 1L ||
      !is.finite(minimum_probability) ||
      minimum_probability < 0 || minimum_probability > 1) {
    stop("minimum_probability must be a single value in [0, 1].")
  }
  kept <- labelled_pairs[
    labelled_pairs$probability >= minimum_probability, , drop = FALSE
  ]
  leads <- kept[
    !duplicated(paste(as.character(kept$class), kept$gene_id)), , drop = FALSE
  ]
  sets <- list()
  for (class_name in TEMPORAL_CLASSES) {
    sets[[paste0(class_name, "_all")]] <-
      unique(kept$variant_id[kept$class == class_name])
    sets[[paste0(class_name, "_lead")]] <-
      unique(leads$variant_id[leads$class == class_name])
  }
  sets
}

#' One row per class: pairs, unique variants, unique genes, and the median
#' winning probability.
summarise_temporal_classes <- function(labelled_pairs) {
  rows <- lapply(TEMPORAL_CLASSES, function(class_name) {
    subset <- labelled_pairs[labelled_pairs$class == class_name, , drop = FALSE]
    data.frame(
      class = class_name,
      label = unname(TEMPORAL_CLASS_LABELS[[class_name]]),
      pairs = nrow(subset),
      variants = length(unique(subset$variant_id)),
      genes = length(unique(subset$gene_id)),
      median_probability = if (nrow(subset)) {
        stats::median(subset$probability)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Pairs whose class call is confident at the manuscript's own cfsr threshold.
#'
#' Retained for the page's power note: these are the sets a strict reading
#' would use, and they are too small for an enrichment estimate.
count_confident_calls <- function(testing_results, alpha = 0.05) {
  rows <- lapply(TEMPORAL_CLASSES, function(class_name) {
    result <- testing_results[[class_name]]
    if (!"cfsr" %in% names(result)) {
      stop("The ", class_name, " testing result has no cfsr column.")
    }
    keys <- rownames(result)[as.numeric(result$cfsr) <= alpha]
    data.frame(
      class = class_name,
      label = unname(TEMPORAL_CLASS_LABELS[[class_name]]),
      pairs = length(keys),
      variants = length(unique(sub("^[^_]+_", "", keys))),
      genes = length(unique(sub("_.*$", "", keys))),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
