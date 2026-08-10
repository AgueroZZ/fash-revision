#!/usr/bin/env Rscript

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
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "fash_strober_enhancer_comparison",
  "fash_strober_enhancer_helpers.R"
))

calls <- select_cumulative_lfdr_calls(c(0.01, 0.04, 0.20), alpha = 0.05)
stopifnot(identical(sort(calls), 1:2))

pair_table <- data.frame(
  pair_key = c(
    "gene1_rs_shared", "gene1_rs_unique_b", "gene1_rs_unique_a",
    "gene2_rs_cross_gene", "gene3_rs_only"
  ),
  gene_id = c("gene1", "gene1", "gene1", "gene2", "gene3"),
  variant_id = c(
    "rs_shared", "rs_unique_b", "rs_unique_a", "rs_cross_gene", "rs_only"
  ),
  stringsAsFactors = FALSE
)
ranked <- derive_ranked_discovery_sets(
  pair_table,
  score = c(0.001, 0.01, 0.01, 0.02, 0.03)
)
stopifnot(
  ranked$lead_pairs$variant_id[ranked$lead_pairs$gene_id == "gene1"] ==
    "rs_shared"
)

# Variant-level exclusivity applied to the lead set: a gene whose best variant
# is shared with Strober drops out entirely rather than falling back.
lead_only <- ranked$lead_pairs[
  !ranked$lead_pairs$variant_id %in% c("rs_shared", "rs_cross_gene"),
  ,
  drop = FALSE
]
stopifnot(
  !"gene1" %in% lead_only$gene_id,   # gene1's lead was rs_shared -> gene dropped
  !"gene2" %in% lead_only$gene_id,
  identical(lead_only$variant_id, "rs_only")
)

unmatched_matrix <- data.frame(
  variant_id = c("v1", "v2", "v3", "v4"),
  enhancer = c(TRUE, FALSE, TRUE, FALSE),
  stringsAsFactors = FALSE
)
unmatched <- compute_unmatched_enrichment(
  unmatched_matrix,
  list(both = c("v1", "v3"), mixed = c("v1", "v2"), absent = c("v2", "v9")),
  "enhancer"
)
stopifnot(
  nrow(unmatched) == 3L,
  abs(unmatched$unmatched_background_rate[1L] - 0.5) < 1e-12,
  abs(unmatched$unmatched_enrichment[
    unmatched$discovery_set == "both"
  ] - 2) < 1e-12,
  abs(unmatched$unmatched_enrichment[
    unmatched$discovery_set == "mixed"
  ] - 1) < 1e-12,
  # Variants absent from the annotation matrix are dropped, not counted as zero.
  unmatched$unmatched_selected_total[
    unmatched$discovery_set == "absent"
  ] == 1L
)

metadata <- data.frame(
  discovery_set = c("a_all", "b_all", "a_lead", "b_lead"),
  view = "Ordinary",
  selection_strategy = c("All variants", "All variants", "Lead variants",
                         "Lead variants"),
  stringsAsFactors = FALSE
)
overlap <- build_discovery_overlap(
  list(
    a_all = c("x", "y"),
    b_all = c("y", "z"),
    a_lead = "x",
    b_lead = c("x", "z")
  ),
  metadata
)
target <- overlap[
  overlap$first_set == "a_all" & overlap$second_set == "b_all",
  ,
  drop = FALSE
]
stopifnot(nrow(target) == 1L, abs(target$jaccard - 1 / 3) < 1e-12)

enhancer_toy <- data.frame(
  annotation_system = "custom",
  discovery_set = rep(c("a", "b"), each = 2L),
  method = rep(c("FASH", "Linear"), each = 2L),
  view = "Ordinary",
  selection_strategy = "All variants",
  log2_enrichment = c(0.1, -0.1, 0.2, 0.3),
  enrichment = 2^c(0.1, -0.1, 0.2, 0.3),
  p_value = c(0.8, 0.9, 0.2, 0.3),
  q_value_within_set = c(0.9, 0.9, 0.5, 0.5),
  stringsAsFactors = FALSE
)
enhancer_summary <- summarize_enhancer_panel(enhancer_toy)
stopifnot(
  nrow(enhancer_summary) == 2L,
  enhancer_summary$positive_annotation_count[
    enhancer_summary$discovery_set == "a"
  ] == 1L,
  enhancer_summary$positive_annotation_count[
    enhancer_summary$discovery_set == "b"
  ] == 2L
)

message("FASH/Strober enhancer helper tests passed.")
