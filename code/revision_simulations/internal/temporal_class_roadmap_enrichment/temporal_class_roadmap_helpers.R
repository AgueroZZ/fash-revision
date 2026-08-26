#!/usr/bin/env Rscript

# =============================================================================
# FASH temporal classes across the thirteen-epigenome Roadmap panel.
#
# Two things this file adds to the existing temporal-class experiment:
#
#   1. the SWITCH class, from the fourth cached functional produced by
#      code/02_dyn_lfsr.R;
#   2. lead-per-gene selection by the CLASS-SPECIFIC lfsr rather than by the
#      overall cumulative lfdr.
#
# Everything else is borrowed. The early/middle/late assignment comes from
# fash_temporal_class_enhancer_enrichment/temporal_class_helpers.R unchanged, the
# annotation table from roadmap_enhancer_exploration, and the estimator from
# fash_strober_enhancer_comparison/enrichment_api.R. Nothing is refit.
#
# A structural point that the page repeats and that this file encodes: the three
# timing windows partition "where does |beta(t)| peak", while switch asks "does
# beta(t) cross zero with magnitude at least 0.25 on both sides". Those are
# different questions, so SWITCH IS NOT A FOURTH MUTUALLY EXCLUSIVE CLASS. A
# pair can be both late and switch, and many are.
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

# TEMPORAL_CLASSES, derive_window_probabilities(), assign_temporal_class(),
# label_discovered_pairs(). Sourced, not duplicated, so the timing classes on
# this page are byte-for-byte the ones on the three-epigenome page.
source(file.path(
  WORKFLOWR_ROOT, "code", "revision_simulations", "internal",
  "fash_temporal_class_enhancer_enrichment", "temporal_class_helpers.R"
))

# roadmap_output_directory(), load_epigenome_panel(),
# enhancer_annotation_name(), ENHANCER_DEFINITIONS, decorate_enrichment().
source(file.path(
  WORKFLOWR_ROOT, "code", "revision_simulations", "internal",
  "roadmap_enhancer_exploration", "roadmap_enhancer_exploration_helpers.R"
))

EXPERIMENT_ID_TCR <- "temporal_class_roadmap_enrichment"

tcr_output_directory <- function() {
  file.path(
    WORKFLOWR_ROOT, "output", "revision_simulations", "internal",
    EXPERIMENT_ID_TCR
  )
}

# -----------------------------------------------------------------------------
# The four classes
# -----------------------------------------------------------------------------

# The switch functional in code/02_dyn_lfsr.R, quoted so the definition is
# legible here without opening that file:
#
#   switch_threshold <- 0.25
#   functional_switch <- function(x) {
#     x_pos <- x[x > 0]; x_neg <- x[x < 0]
#     if (length(x_pos) == 0 || length(x_neg) == 0) return(0)
#     min(max(abs(x_pos)), max(abs(x_neg))) - switch_threshold
#   }
#
# so the functional is positive when beta(t) reaches at least 0.25 in absolute
# value on BOTH sides of zero over the 15-day course, and `lfsr` is the posterior
# probability that it is not.
SWITCH_CLASS <- "switch"
SWITCH_THRESHOLD <- 0.25

# WHICH MIDDLE WINDOW DOES THE CACHE USE?
#
# The middle functional is the open interval,
#
#   sup_{3 < t < 12}|b| - sup_{t <= 3 or t >= 12}|b|
#
# so early [0,3], middle (3,12) and late [12,15] partition the time course
# exactly. A closed [4,11] variant used to circulate here; it picked the same
# sixteen measured days but left (3,4) and (11,12) uncovered on the
# 0.1-resolution grid the functionals are evaluated on, and the two cannot be
# told apart from the saved lfsr alone. So this is checked at run time.
#
# The functional is exactly 0 under the null (b(t) == 0 gives 0 - 0), and
# `lfsr_cal = mean(x <= 0)` counts that against membership, so the three
# membership probabilities must satisfy
#
#     P(early) + P(middle) + P(late) = 1 - pi0
#
# where pi0 is the posterior null probability. A closed window would fall short
# by the mass in the two gaps. verify_middle_definition() tests this identity,
# which is what actually pins the definition -- comparing the sum against 1
# instead of 1 - pi0 mistakes pi0 for a gap and reaches the wrong conclusion.
MIDDLE_WINDOW_TOLERANCE <- 0.005

TCR_CLASSES <- c(TEMPORAL_CLASSES, SWITCH_CLASS)

TCR_CLASS_LABELS <- c(
  early = "Early (peak on days 0-3)",
  middle = "Middle (peak at 3 < t < 12)",
  late = "Late (peak on days 12-15)",
  switch = "Switch (sign change, |beta| >= 0.25 both sides)"
)

TCR_CLASS_SHORT_LABELS <- c(
  early = "Early",
  middle = "Middle",
  late = "Late",
  switch = "Switch"
)

# The manuscript's own cumulative call. `cfsr` is the cumulative false-sign rate
# testing_functional() reports alongside `lfsr`; thresholding it at 0.05 is the
# same rule the manuscript uses to call a class, and for switch it selects 981
# of the 9,214 discovered pairs -- the only one of the four classes where a
# strict call leaves enough variants to estimate anything.
SWITCH_CFSR_ALPHA <- 0.05

SELECTION_STRATEGIES <- c(all = "All variants", lead = "One lead variant per gene")

# Lead-selection rules. `lfsr` is what this page uses and what was asked for:
# within a class, each gene keeps the variant whose lfsr for THAT class's
# functional is smallest, i.e. the variant most confidently in the class. `lfdr`
# is the rule the three-epigenome sister page uses -- the overall cumulative
# lfdr the pair was discovered by -- kept so the two pages can be reconciled.
LEAD_RULES <- c(
  lfsr = "Lead by class-specific lfsr",
  lfdr = "Lead by overall cumulative lfdr"
)

# -----------------------------------------------------------------------------
# Cached inputs
# -----------------------------------------------------------------------------

#' Load one cached `testing_functional()` result.
#'
#' @param class_name one of early, middle, late, switch.
#' @return the data frame as saved by `code/02_dyn_lfsr.R`, with `lfsr`, `cfsr`, and
#'   `gene_variant` row names.
load_classification <- function(class_name,
                                directory = file.path(
                                  WORKFLOWR_ROOT, "output", "dynamic_eQTL_real"
                                )) {
  if (!class_name %in% TCR_CLASSES) {
    stop("Unknown class: ", class_name)
  }
  path <- file.path(
    directory, paste0("classify_dyn_eQTLs_", class_name, ".RData")
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
  result <- get(object_name, envir = environment_for_load)
  required <- c("lfsr", "cfsr")
  if (!is.data.frame(result) || !all(required %in% names(result)) ||
      is.null(rownames(result)) || anyDuplicated(rownames(result))) {
    stop("The cached ", class_name, " classification has an unexpected shape.")
  }
  result
}

#' All four cached classifications, checked to be row-aligned.
load_all_classifications <- function(...) {
  results <- lapply(TCR_CLASSES, load_classification, ...)
  names(results) <- TCR_CLASSES
  reference <- rownames(results[[1L]])
  aligned <- vapply(results, function(result) {
    identical(rownames(result), reference)
  }, logical(1))
  if (!all(aligned)) {
    stop(
      "The cached classifications are not row-aligned: ",
      paste(names(aligned)[!aligned], collapse = ", ")
    )
  }
  results
}

#' The per-pair posterior null probability, pi0.
#'
#' `fash_fit1_update$lfdr` is the posterior weight on the psd = 0 component,
#' i.e. the probability the pair has no effect at all. It is cached as a small
#' RDS so that verifying the middle definition does not require loading the
#' 554 MB fit every time; pass `refresh = TRUE` to re-extract it.
load_posterior_null_probability <- function(
  cache_path = file.path(tcr_output_directory(), "posterior_null_probability.rds"),
  fit_path = file.path(
    WORKFLOWR_ROOT, "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
  ),
  refresh = FALSE
) {
  if (!refresh && file.exists(cache_path)) {
    cached <- readRDS(cache_path)
    if (!all(c("indices", "posterior_null_probability") %in% names(cached))) {
      stop("The cached posterior null probabilities have an unexpected shape.")
    }
    return(cached)
  }
  if (!file.exists(fit_path)) {
    stop("Cannot extract pi0: missing ", fit_path)
  }
  environment_for_load <- new.env(parent = emptyenv())
  loaded <- load(fit_path, envir = environment_for_load)
  fit <- get(loaded[1L], envir = environment_for_load)
  if (!"lfdr" %in% names(fit) || !"prior_weights" %in% names(fit) ||
      fit$prior_weights$psd[1L] != 0) {
    stop("The fit does not expose a psd = 0 null component.")
  }
  indices <- load_classification("early")$indices
  output <- list(
    indices = indices,
    posterior_null_probability = unname(fit$lfdr[indices]),
    source = paste(
      "fash_fit1_update$lfdr (== posterior_weights[, psd = 0])",
      "at the classification indices"
    ),
    psd_grid_null_component = fit$prior_weights$psd[1L],
    extracted_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(output, cache_path)
  output
}

#' Test which middle window produced the cached classifications.
#'
#' Returns the diagnostic and, unless `strict = FALSE`, stops if the cached
#' probabilities are inconsistent with the open (3, 12) window. This is the guard
#' that would fire if a legacy closed-Middle cache were supplied.
verify_middle_definition <- function(classifications,
                                     null_probability =
                                       load_posterior_null_probability(),
                                     tolerance = MIDDLE_WINDOW_TOLERANCE,
                                     strict = TRUE) {
  indices <- classifications[[TEMPORAL_CLASSES[1L]]]$indices
  if (!identical(unname(indices), unname(null_probability$indices))) {
    stop("The cached pi0 does not belong to these classifications.")
  }
  probabilities <- vapply(
    TEMPORAL_CLASSES,
    function(class_name) 1 - as.numeric(classifications[[class_name]]$lfsr),
    numeric(length(indices))
  )
  total <- rowSums(probabilities)
  pi0 <- null_probability$posterior_null_probability
  # Against the correct reference the residual is Monte Carlo noise centred at
  # zero; against 1 it is centred at -pi0, which is the trap.
  residual_open <- total - (1 - pi0)
  residual_versus_one <- total - 1
  diagnostic <- list(
    n_pairs = length(total),
    mean_residual_open = mean(residual_open),
    median_residual_open = stats::median(residual_open),
    max_absolute_residual_open = max(abs(residual_open)),
    mean_residual_versus_one = mean(residual_versus_one),
    median_pi0 = stats::median(pi0),
    window = if (abs(mean(residual_open)) <= tolerance) "open (3, 12)" else
      "inconsistent with the open (3, 12) window"
  )
  if (strict && abs(diagnostic$mean_residual_open) > tolerance) {
    stop(
      "The cached classifications do not satisfy ",
      "P(early)+P(middle)+P(late) = 1 - pi0 (mean residual ",
      format(diagnostic$mean_residual_open, digits = 3),
      "). They were probably produced by the historical closed [4, 11] ",
      "Middle window rather than the current open (3, 12) window."
    )
  }
  diagnostic
}

# -----------------------------------------------------------------------------
# Class membership
# -----------------------------------------------------------------------------

#' Attach class membership and per-class lfsr to the discovered-pair table.
#'
#' Timing membership is the argmax of the three window probabilities, delegated
#' to `assign_temporal_class()`. Switch membership is the independent
#' `cfsr <= alpha` call. The output therefore has one `in_<class>` column per
#' class, and the timing three are mutually exclusive while switch is not.
#'
#' @param pair_table R6's `current_all` table: `pair_key`, `gene_id`,
#'   `variant_id`, `score` (the cumulative lfdr the pair was ranked by).
#' @param classifications output of `load_all_classifications()`.
#' @return `pair_table` with `timing_class`, `timing_probability`, one
#'   `in_<class>` logical per class, and one `lfsr_<class>` numeric per class.
build_class_membership <- function(pair_table,
                                  classifications,
                                  switch_alpha = SWITCH_CFSR_ALPHA) {
  required <- c("pair_key", "gene_id", "variant_id", "score")
  if (!all(required %in% names(pair_table)) ||
      anyDuplicated(pair_table$pair_key)) {
    stop("The discovered-pair table has an unexpected shape.")
  }
  if (length(switch_alpha) != 1L || !is.finite(switch_alpha) ||
      switch_alpha <= 0 || switch_alpha >= 1) {
    stop("switch_alpha must be a single value in (0, 1).")
  }

  timing <- assign_temporal_class(
    derive_window_probabilities(classifications[TEMPORAL_CLASSES])
  )
  # This is the join that ties the cached classification to the R6 discovery
  # set: it fails unless the two pair-key sets are exactly equal.
  labelled <- label_discovered_pairs(pair_table, timing)
  names(labelled)[names(labelled) == "class"] <- "timing_class"
  names(labelled)[names(labelled) == "probability"] <- "timing_probability"

  # Every classification is indexed below with a single `rows` vector, so they
  # must all carry the same row names in the same order. derive_window_
  # probabilities() only checks the three timing tables, so switch would
  # otherwise be read against the wrong pairs without complaint.
  keys <- rownames(classifications[[1L]])
  misaligned <- names(classifications)[!vapply(
    classifications, function(result) identical(rownames(result), keys),
    logical(1)
  )]
  if (length(misaligned)) {
    stop(
      "These cached classifications are not row-aligned with the others: ",
      paste(misaligned, collapse = ", ")
    )
  }
  rows <- match(labelled$pair_key, keys)
  if (anyNA(rows)) {
    stop("A discovered pair has no cached classification.")
  }

  for (class_name in TCR_CLASSES) {
    labelled[[paste0("lfsr_", class_name)]] <-
      as.numeric(classifications[[class_name]]$lfsr[rows])
    labelled[[paste0("cfsr_", class_name)]] <-
      as.numeric(classifications[[class_name]]$cfsr[rows])
  }
  for (class_name in TEMPORAL_CLASSES) {
    labelled[[paste0("in_", class_name)]] <-
      labelled$timing_class == class_name
  }
  labelled$in_switch <- labelled$cfsr_switch <= switch_alpha

  # The timing three must partition the discoveries; switch must not be forced
  # to.
  timing_membership <- rowSums(as.matrix(
    labelled[, paste0("in_", TEMPORAL_CLASSES), drop = FALSE]
  ))
  if (any(timing_membership != 1L)) {
    stop("The timing classes do not partition the discovered pairs.")
  }
  labelled
}

#' The variant sets the enrichment API consumes.
#'
#' Lead selection happens *within* a class, so a gene whose variants split
#' across classes contributes one lead to each class it appears in. Under the
#' `lfsr` rule the lead is the variant with the smallest lfsr for that class's
#' own functional; under `lfdr` it is the lowest overall cumulative lfdr, which
#' is what the three-epigenome sister page uses.
#'
#' @return named list of unique variant-ID vectors, `<class>_all` and
#'   `<class>_lead`.
build_class_variant_sets <- function(membership, lead_rule = "lfsr") {
  if (!lead_rule %in% names(LEAD_RULES)) {
    stop("lead_rule must be one of: ", paste(names(LEAD_RULES), collapse = ", "))
  }
  sets <- list()
  for (class_name in TCR_CLASSES) {
    inside <- membership[membership[[paste0("in_", class_name)]], , drop = FALSE]
    if (!nrow(inside)) {
      stop("Class ", class_name, " is empty.")
    }
    # The lead rule changes only which variant represents a gene, never the
    # membership of the class. So the `_all` set is ordered by the discovery
    # score regardless, and is identical across rules down to vector order.
    ordering <- if (identical(lead_rule, "lfsr")) {
      order(inside[[paste0("lfsr_", class_name)]], inside$score,
            inside$variant_id, method = "radix")
    } else {
      order(inside$score, inside$variant_id, method = "radix")
    }
    leads <- inside[ordering, , drop = FALSE]
    leads <- leads[!duplicated(leads$gene_id), , drop = FALSE]
    sets[[paste0(class_name, "_all")]] <- unique(
      inside$variant_id[order(inside$score, inside$variant_id, method = "radix")]
    )
    sets[[paste0(class_name, "_lead")]] <- unique(leads$variant_id)
  }
  sets
}

#' The non-switch discoveries, as the contrast for the switch class.
build_no_switch_sets <- function(membership, lead_rule = "lfsr") {
  outside <- membership[!membership$in_switch, , drop = FALSE]
  if (!nrow(outside)) {
    stop("Every discovered pair is in the switch class.")
  }
  ordering <- if (identical(lead_rule, "lfsr")) {
    # Furthest from the switch call, i.e. largest lfsr, is the most
    # representative non-switch variant for a gene.
    order(-outside$lfsr_switch, outside$score, outside$variant_id,
          method = "radix")
  } else {
    order(outside$score, outside$variant_id, method = "radix")
  }
  leads <- outside[ordering, , drop = FALSE]
  leads <- leads[!duplicated(leads$gene_id), , drop = FALSE]
  list(
    no_switch_all = unique(
      outside$variant_id[
        order(outside$score, outside$variant_id, method = "radix")
      ]
    ),
    no_switch_lead = unique(leads$variant_id)
  )
}

#' Split a set name back into its class and selection strategy.
parse_set_name <- function(set_names) {
  strategy <- sub("^.*_", "", set_names)
  class_name <- sub("_[^_]+$", "", set_names)
  if (!all(strategy %in% names(SELECTION_STRATEGIES))) {
    stop("Unparseable selection strategy in: ", paste(set_names, collapse = ", "))
  }
  data.frame(
    variant_set = set_names,
    class = class_name,
    strategy = strategy,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Descriptive tables
# -----------------------------------------------------------------------------

#' One row per class: pairs, variants, genes, leads, and the median class lfsr.
summarise_classes <- function(membership, lead_rule = "lfsr") {
  sets <- build_class_variant_sets(membership, lead_rule = lead_rule)
  rows <- lapply(TCR_CLASSES, function(class_name) {
    inside <- membership[membership[[paste0("in_", class_name)]], , drop = FALSE]
    data.frame(
      class = class_name,
      label = unname(TCR_CLASS_LABELS[[class_name]]),
      mutually_exclusive = class_name %in% TEMPORAL_CLASSES,
      pairs = nrow(inside),
      variants = length(unique(inside$variant_id)),
      genes = length(unique(inside$gene_id)),
      leads = length(sets[[paste0(class_name, "_lead")]]),
      median_class_lfsr = stats::median(
        inside[[paste0("lfsr_", class_name)]]
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Switch against the timing classes, as counts.
#'
#' The point of this table is to stop the reader treating switch as a fourth
#' bucket: it overlaps all three timing windows.
crosstab_switch_by_timing <- function(membership) {
  rows <- lapply(TEMPORAL_CLASSES, function(class_name) {
    inside <- membership[membership$timing_class == class_name, , drop = FALSE]
    data.frame(
      timing_class = class_name,
      label = unname(TCR_CLASS_LABELS[[class_name]]),
      pairs = nrow(inside),
      switch_pairs = sum(inside$in_switch),
      switch_share = if (nrow(inside)) mean(inside$in_switch) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, rows)
  rbind(
    output,
    data.frame(
      timing_class = "all",
      label = "All discovered pairs",
      pairs = nrow(membership),
      switch_pairs = sum(membership$in_switch),
      switch_share = mean(membership$in_switch),
      stringsAsFactors = FALSE
    )
  )
}

#' How many variants a set needs before a 10-overlap estimate is possible.
#'
#' Reported up front so blank cells in the figures read as a power limit rather
#' than as an absence of effect.
estimability_forecast <- function(annotation_matrix,
                                  panel = load_epigenome_panel(),
                                  definition = "R6",
                                  minimum_overlap = 10L) {
  rows <- lapply(seq_len(nrow(panel)), function(index) {
    column <- enhancer_annotation_name(
      panel$epigenome_id[index], panel$roadmap_std_name[index], definition
    )
    rate <- mean(as.logical(annotation_matrix[[column]]))
    data.frame(
      epigenome_id = panel$epigenome_id[index],
      display_label = panel$display_label[index],
      biological_group = panel$biological_group[index],
      background_rate = rate,
      variants_needed = ceiling(minimum_overlap / rate),
      stringsAsFactors = FALSE
    )
  })
  output <- do.call(rbind, rows)
  output[order(output$variants_needed), , drop = FALSE]
}
