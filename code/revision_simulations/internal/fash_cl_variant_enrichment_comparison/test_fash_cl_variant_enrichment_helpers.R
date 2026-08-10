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
  "fash_cl_variant_enrichment_comparison",
  "fash_cl_variant_enrichment_helpers.R"
))

toy_matrix <- matrix(c(
  0, -2, -4,
  -1000, -1001, -1003,
  4, 3, 2
), nrow = 3L, byrow = TRUE)
toy_log_sum <- row_log_sum_exp(toy_matrix)
expected_log_sum <- apply(toy_matrix, 1L, function(row) {
  maximum <- max(row)
  maximum + log(sum(exp(row - maximum)))
})
stopifnot(max(abs(toy_log_sum - expected_log_sum)) < 1e-12)

toy_likelihood <- matrix(c(
  0.0, -2.0, -4.0,
  0.0, -1.5, -3.0,
  0.0, 0.5, -1.0,
  0.0, 1.0, 0.0,
  0.0, 2.0, 0.5,
  0.0, 3.0, 1.0,
  0.0, -0.5, -2.0,
  0.0, 1.5, 0.2
), nrow = 8L, byrow = TRUE)
toy_adjustment <- compute_bf_adjusted_lfdr(
  toy_likelihood,
  chunk_size = 3L,
  verbose = FALSE
)
toy_log_alt <- compute_log_bayes_factor(
  toy_likelihood,
  toy_adjustment$alternative_weights,
  chunk_size = 3L
)
toy_direct <- if (toy_adjustment$null_weight >= 1) {
  rep(1, nrow(toy_likelihood))
} else {
  stats::plogis(
    stats::qlogis(toy_adjustment$null_weight) -
      toy_log_alt
  )
}
stopifnot(max(abs(toy_adjustment$lfdr - toy_direct)) < 1e-12)

pair_keys <- c(
  "gene1_rs2", "gene1_rs1", "gene2_rs1", "gene2_rs3", "gene3_rs4"
)
lfdr <- c(0.01, 0.01, 0.02, 0.30, 0.90)
sets <- derive_all_and_lead_sets(pair_keys, lfdr, alpha = 0.05)
stopifnot(
  identical(sets$all_variants, c("rs1", "rs2")),
  identical(sets$lead_variants, "rs1"),
  nrow(sets$discovered_pairs) == 3L,
  nrow(sets$lead_pairs) == 2L,
  sets$lead_pairs$variant_id[sets$lead_pairs$gene_id == "gene1"] == "rs1"
)

overlap <- summarize_method_overlap(
  current_sets = list(all = c("a", "b"), one_lead_per_gene = "a"),
  fash_cl_sets = list(all = c("b", "c"), one_lead_per_gene = c("a", "c"))
)
stopifnot(
  overlap$intersection_count[overlap$selection_strategy == "all"] == 1L,
  abs(overlap$jaccard[overlap$selection_strategy == "all"] - 1 / 3) < 1e-12
)

enhancer_toy <- data.frame(
  annotation_system = rep("custom", 4L),
  method = rep(c("Current FASH", "FASH-CL"), each = 2L),
  selection_strategy = "all",
  annotation = rep(c("enhancer_a", "enhancer_b"), 2L),
  log2_enrichment = c(0.1, -0.1, 0.2, 0.3),
  enrichment = 2^c(0.1, -0.1, 0.2, 0.3),
  p_value = c(0.5, 0.6, 0.2, 0.1),
  q_value_within_set = c(0.8, 0.8, 0.4, 0.4),
  stringsAsFactors = FALSE
)
enhancer_summary <- summarize_enhancer_panel(enhancer_toy)
stopifnot(
  nrow(enhancer_summary) == 2L,
  enhancer_summary$positive_annotation_count[
    enhancer_summary$method == "Current FASH"
  ] == 1L,
  enhancer_summary$positive_annotation_count[
    enhancer_summary$method == "FASH-CL"
  ] == 2L
)

message("FASH-CL enrichment helper tests passed.")
