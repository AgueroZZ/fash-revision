#!/usr/bin/env Rscript

# =============================================================================
# Variant-set enrichment, in three steps.
#
#   1. ANNOTATIONS  a table saying which variants belong to which group
#                   (one logical column per group, e.g. "is this SNP inside a
#                   mesoderm enhancer?")
#
#   2. INPUT        a variant set to test  (e.g. variants discovered only by
#                                           FASH, one per gene)
#                   plus a background      (every variant either method tested)
#
#   3. OUTPUT       one enrichment score per variant set per annotation group
#
# Everything else in this directory is bookkeeping around those three objects.
# This file is the whole analysis; read it and you know what was done.
#
#   annotations <- load_annotation_groups("custom")
#   background  <- load_background()
#   sets        <- load_variant_sets()
#
#   compute_enrichment(sets["current_only_lead"], background, annotations)
#
# `compute_enrichment` takes any named list of variant vectors, so a reader can
# substitute their own sets without touching anything else.
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

# The matched-control sampler and the chromosome jackknife live here. They are
# shared with the other enrichment pages, so they are sourced rather than
# duplicated.
source(file.path(
  WORKFLOWR_ROOT, "code", "revision_simulations", "internal",
  "variant_annotation_enrichment", "variant_annotation_enrichment_helpers.R"
))

MATCHING_SEEDS <- seq.int(20260807L, length.out = 100L)
CONTROLS_PER_VARIANT <- 5L
AUTOSOMES <- as.character(1:22)
MINIMUM_OVERLAP <- 10L

ENHANCER_GROUPS <- list(
  custom = c(
    "ENCODE cCRE enhancer-like",
    "Roadmap E020 iPS-20b: Enhancer",
    "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer",
    "Roadmap E095 left ventricle: Enhancer"
  ),
  baselineld = c(
    "Enhancer_Andersson", "Enhancer_Hoffman", "WeakEnhancer_Hoffman",
    "SuperEnhancer_Hnisz", "Human_Enhancer_Villar"
  )
)

internal_path <- function(...) {
  file.path(
    WORKFLOWR_ROOT, "output", "revision_simulations", "internal", ...
  )
}

# -----------------------------------------------------------------------------
# Step 1: annotations
# -----------------------------------------------------------------------------

#' Which variants belong to which annotation group.
#'
#' @param system "custom" (GENCODE / ENCODE cCRE / Roadmap ChromHMM, built by
#'   variant_annotation_enrichment/) or "baselineld" (baselineLD v2.2 as
#'   published).
#' @param groups annotation columns to keep; defaults to the pre-specified
#'   enhancer panel for that system. Pass `NULL` after `all_groups = TRUE` to
#'   see everything available.
#' @return data frame: `variant_id`, `chromosome`, then one logical column per
#'   group. A row exists only for variants the system can annotate, so
#'   baselineLD returns fewer rows than custom.
load_annotation_groups <- function(system = c("custom", "baselineld"),
                                   groups = NULL,
                                   all_groups = FALSE) {
  system <- match.arg(system)
  path <- switch(
    system,
    custom = internal_path(
      "variant_annotation_enrichment", "variant_annotation_matrix.rds"
    ),
    baselineld = internal_path(
      "baseline_ld_variant_enrichment", "baseline_ld_binary_annotation_matrix.rds"
    )
  )
  if (!file.exists(path)) {
    stop("Missing annotation matrix: ", path)
  }
  annotations <- readRDS(path)
  available <- setdiff(names(annotations), c("variant_id", "chromosome"))
  if (is.null(groups)) {
    groups <- if (all_groups) available else ENHANCER_GROUPS[[system]]
  }
  missing <- setdiff(groups, available)
  if (length(missing)) {
    stop("Unknown annotation group(s): ", paste(missing, collapse = ", "))
  }
  keep <- intersect(c("variant_id", "chromosome"), names(annotations))
  annotations[, c(keep, groups), drop = FALSE]
}

#' List every annotation group a system offers, for discovery.
list_annotation_groups <- function(system = c("custom", "baselineld")) {
  names(load_annotation_groups(system, all_groups = TRUE))[-(1:2)]
}

# -----------------------------------------------------------------------------
# Step 2: inputs
# -----------------------------------------------------------------------------

#' The background: every variant tested by either method, with the covariates
#' used to draw comparable controls (cohort MAF, distance to the nearest tested
#' target-gene TSS, local tested-variant density, number of tested genes).
load_background <- function() {
  path <- internal_path(
    "variant_annotation_enrichment", "variant_matching_covariates.rds"
  )
  if (!file.exists(path)) {
    stop("Missing background table: ", path)
  }
  readRDS(path)
}

#' The eight published variant sets, as a named list of variant-ID vectors.
load_variant_sets <- function() {
  path <- internal_path(
    "fash_strober_enhancer_comparison", "discovery_sets.rds"
  )
  if (!file.exists(path)) {
    stop(
      "Missing discovery sets: ", path,
      "\nRun run_fash_strober_enhancer_comparison.R first."
    )
  }
  readRDS(path)
}

# -----------------------------------------------------------------------------
# Step 3: output
# -----------------------------------------------------------------------------

#' Delete-one-chromosome jackknife for a log2 ratio of two proportions.
#'
#' Given per-chromosome hit and total counts for the selected variants and for
#' whatever plays the role of controls, return the point estimate, its standard
#' error, a 95% interval, and a nominal p-value. Chromosome blocks absorb
#' within-chromosome correlation, which is why these p-values are far more
#' conservative than a Fisher test that treats variants as independent.
jackknife_log2_ratio <- function(selected_hits,
                                 selected_total,
                                 control_hits,
                                 control_total) {
  ratio <- function(hits, total) {
    if (total <= 0) NA_real_ else hits / total
  }
  point <- log2(
    ratio(sum(selected_hits), sum(selected_total)) /
      ratio(sum(control_hits), sum(control_total))
  )
  leave_one_out <- vapply(seq_along(selected_hits), function(block) {
    log2(
      ratio(
        sum(selected_hits) - selected_hits[block],
        sum(selected_total) - selected_total[block]
      ) /
        ratio(
          sum(control_hits) - control_hits[block],
          sum(control_total) - control_total[block]
        )
    )
  }, numeric(1))
  leave_one_out <- leave_one_out[is.finite(leave_one_out)]
  blocks <- length(leave_one_out)
  standard_error <- if (blocks < 2L) {
    NA_real_
  } else {
    sqrt(
      (blocks - 1) / blocks *
        sum((leave_one_out - mean(leave_one_out))^2)
    )
  }
  list(
    log2_enrichment = point,
    jackknife_se = standard_error,
    ci_lower = point - 1.96 * standard_error,
    ci_upper = point + 1.96 * standard_error,
    p_value = 2 * stats::pnorm(-abs(point / standard_error)),
    n_blocks = blocks
  )
}

# The plain-background comparison, written out longhand because it is the
# clearest statement of what "enrichment" means here:
#
#     fold = (rate in the variant set) / (rate in the background)
#
# The matched version below answers the same question after removing the
# confounding between discovery status and MAF / TSS distance / variant density.
enrichment_versus_background <- function(selected,
                                         background,
                                         annotations,
                                         groups,
                                         minimum_overlap) {
  annotated <- annotations$variant_id
  selected <- intersect(selected, annotated)
  # The background excludes the set being tested, so a set is never part of its
  # own comparator. This mirrors the matched path, where controls are drawn only
  # from non-selected variants.
  control <- setdiff(intersect(background$variant_id, annotated), selected)
  selected_rows <- match(selected, annotated)
  control_rows <- match(control, annotated)
  chromosome <- as.character(annotations$chromosome)
  do.call(rbind, lapply(groups, function(group) {
    hit <- as.logical(annotations[[group]])
    selected_chromosome <- chromosome[selected_rows]
    control_chromosome <- chromosome[control_rows]
    per_block <- function(rows, chromosomes, mask) {
      vapply(AUTOSOMES, function(autosome) {
        sum(mask[rows][chromosomes == autosome], na.rm = TRUE)
      }, numeric(1))
    }
    selected_hits <- per_block(selected_rows, selected_chromosome, hit)
    control_hits <- per_block(control_rows, control_chromosome, hit)
    selected_total <- vapply(AUTOSOMES, function(a) sum(selected_chromosome == a), numeric(1))
    control_total <- vapply(AUTOSOMES, function(a) sum(control_chromosome == a), numeric(1))
    estimate <- if (sum(selected_hits) >= minimum_overlap &&
                    sum(control_hits) >= minimum_overlap) {
      jackknife_log2_ratio(
        selected_hits, selected_total, control_hits, control_total
      )
    } else {
      list(log2_enrichment = NA_real_, jackknife_se = NA_real_,
           ci_lower = NA_real_, ci_upper = NA_real_, p_value = NA_real_,
           n_blocks = NA_integer_)
    }
    data.frame(
      annotation = group,
      selected_overlap = sum(selected_hits),
      selected_total = sum(selected_total),
      control_overlap = sum(control_hits),
      control_total = sum(control_total),
      fold_enrichment = 2^estimate$log2_enrichment,
      log2_enrichment = estimate$log2_enrichment,
      ci_lower_fold = 2^estimate$ci_lower,
      ci_upper_fold = 2^estimate$ci_upper,
      p_value = estimate$p_value,
      stringsAsFactors = FALSE
    )
  }))
}

#' Enrichment of one or more variant sets across annotation groups.
#'
#' @param variant_sets named list of variant-ID vectors. Names become the
#'   `variant_set` column of the output. Your own sets work here.
#' @param background data frame from `load_background()`, or any table with
#'   `variant_id`, `chromosome`, and the four matching covariates.
#' @param annotations data frame from `load_annotation_groups()`.
#' @param controls `"matched"` draws `controls_per_variant` comparable variants
#'   per selected variant from the background, repeated over `seeds`, and pools
#'   the counts. `"background"` compares against the whole background with no
#'   matching - simpler to read, but confounded by the fact that discovered
#'   variants are systematically higher-MAF and closer to a TSS than average,
#'   and enhancers are themselves TSS-proximal.
#' @return one row per (variant set, annotation group), with fold enrichment,
#'   a 95% delete-one-chromosome jackknife interval, and a nominal p-value.
#'   No multiplicity correction is applied.
compute_enrichment <- function(variant_sets,
                               background = load_background(),
                               annotations = load_annotation_groups("custom"),
                               controls = c("matched", "background"),
                               seeds = MATCHING_SEEDS,
                               controls_per_variant = CONTROLS_PER_VARIANT,
                               minimum_overlap = MINIMUM_OVERLAP,
                               verbose = FALSE) {
  controls <- match.arg(controls)
  if (!is.list(variant_sets) || is.null(names(variant_sets))) {
    stop("variant_sets must be a named list of variant-ID vectors.")
  }
  groups <- setdiff(names(annotations), c("variant_id", "chromosome"))
  if (!length(groups)) {
    stop("The annotation table has no group columns.")
  }

  # Both the sets and the background are restricted to variants this annotation
  # system can speak about. For baselineLD that drops roughly a quarter of the
  # tested variants, which are absent from its EUR-anchored SNP list.
  annotated <- annotations$variant_id
  background <- background[background$variant_id %in% annotated, , drop = FALSE]
  variant_sets <- lapply(variant_sets, function(ids) {
    intersect(unique(as.character(ids)), annotated)
  })
  too_small <- names(variant_sets)[lengths(variant_sets) < minimum_overlap]
  if (length(too_small)) {
    stop(
      "Fewer than ", minimum_overlap, " annotated variants in: ",
      paste(too_small, collapse = ", ")
    )
  }

  if (identical(controls, "background")) {
    rows <- lapply(names(variant_sets), function(name) {
      output <- enrichment_versus_background(
        variant_sets[[name]], background, annotations, groups, minimum_overlap
      )
      cbind(variant_set = name, controls = "background", output)
    })
    return(do.call(rbind, rows))
  }

  if (!"full_stratum" %in% names(background)) {
    background <- derive_matching_strata(background, n_bins = 10L)
  }
  analysis <- run_streaming_matched_enrichment(
    variant_table = background,
    annotation_matrix = annotations,
    selected_sets = variant_sets,
    seeds = seeds,
    controls_per_variant = controls_per_variant,
    chromosomes = AUTOSOMES,
    minimum_overlap_count = minimum_overlap,
    verbose = verbose
  )
  results <- analysis$results
  data.frame(
    variant_set = results$discovery_set,
    controls = "matched",
    annotation = results$annotation,
    selected_overlap = results$selected_overlap,
    selected_total = results$selected_total,
    control_overlap = results$control_overlap,
    control_total = results$control_total,
    fold_enrichment = results$enrichment,
    log2_enrichment = results$log2_enrichment,
    ci_lower_fold = 2^results$ci_lower,
    ci_upper_fold = 2^results$ci_upper,
    p_value = results$p_value,
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Check that this interface reproduces the published page
# -----------------------------------------------------------------------------

#' Recompute the eight published sets through `compute_enrichment` and diff the
#' result against the estimates the page reports. Guards against this readable
#' wrapper drifting away from the analysis it is meant to describe.
verify_against_published <- function(tolerance = 1e-8) {
  published <- utils::read.csv(
    internal_path("fash_strober_enhancer_comparison", "published_estimates.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  sets <- load_variant_sets()
  background <- load_background()
  systems <- c(custom = "Custom regulatory", baselineld = "baselineLD v2.2")
  reproduced <- do.call(rbind, lapply(names(systems), function(system) {
    output <- compute_enrichment(
      sets,
      background = background,
      annotations = load_annotation_groups(system),
      controls = "matched",
      verbose = TRUE
    )
    output$annotation_system <- systems[[system]]
    output
  }))
  key <- function(system, set, annotation) paste(system, set, annotation)
  merged <- merge(
    data.frame(
      key = key(published$annotation_system, published$discovery_set,
                published$annotation),
      published_fold = published$enrichment,
      published_p = published$p_value,
      stringsAsFactors = FALSE
    ),
    data.frame(
      key = key(reproduced$annotation_system, reproduced$variant_set,
                reproduced$annotation),
      fold = reproduced$fold_enrichment,
      p = reproduced$p_value,
      stringsAsFactors = FALSE
    ),
    by = "key"
  )
  if (nrow(merged) != nrow(published)) {
    stop("Matched ", nrow(merged), " of ", nrow(published), " published estimates.")
  }
  gap <- function(a, b) {
    delta <- abs(a - b)
    delta[!is.finite(a) & !is.finite(b)] <- 0
    delta
  }
  worst <- max(
    gap(merged$published_fold, merged$fold),
    gap(merged$published_p, merged$p),
    na.rm = TRUE
  )
  failures <- sum(
    gap(merged$published_fold, merged$fold) > tolerance |
      gap(merged$published_p, merged$p) > tolerance,
    na.rm = TRUE
  )
  message(
    "compute_enrichment reproduced ", nrow(merged), " published estimates, ",
    failures, " mismatches, maximum absolute difference ",
    format(worst, digits = 3), "."
  )
  if (failures > 0L) {
    stop("compute_enrichment does not match the published page.")
  }
  invisible(merged)
}

if (!interactive() && identical(sys.nframe(), 0L)) {
  arguments <- commandArgs(trailingOnly = TRUE)
  if ("--verify" %in% arguments) {
    verify_against_published()
  } else {
    cat(
      "Three steps:\n",
      "  annotations <- load_annotation_groups(\"custom\")\n",
      "  background  <- load_background()\n",
      "  sets        <- load_variant_sets()\n",
      "  compute_enrichment(sets, background, annotations)\n\n",
      "Check this interface against the published page:\n",
      "  Rscript enrichment_api.R --verify\n",
      sep = ""
    )
  }
}
