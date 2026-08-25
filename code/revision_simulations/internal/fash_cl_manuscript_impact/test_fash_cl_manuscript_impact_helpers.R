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
helper_path <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "fash_cl_manuscript_impact",
  "fash_cl_manuscript_impact_helpers.R"
)
source(helper_path)

parsed <- parse_pair_keys(c("ENSG1_rs1", "ENSG2_rs_two"))
stopifnot(
  identical(parsed$gene_id, c("ENSG1", "ENSG2")),
  identical(parsed$variant_id, c("rs1", "rs_two"))
)

toy_log <- matrix(c(
  -1000, -1001,
  1000, 999,
  -2, -2
), nrow = 3L, byrow = TRUE)
expected_log_sum <- apply(toy_log, 1L, function(row) {
  maximum <- max(row)
  maximum + log(sum(exp(row - maximum)))
})
stopifnot(max(abs(row_log_sum_exp(toy_log) - expected_log_sum)) < 1e-12)

toy_lfdr <- c(0.01, 0.02, 0.20, 0.90)
stopifnot(identical(sort(select_cumulative_lfdr_calls(toy_lfdr, 0.05)), 1:2))

toy_fit <- list(
  L_matrix = matrix(c(
    0.0, -1.0, -2.0,
    0.0, 1.0, -1.0,
    0.0, -0.5, 0.5
  ), nrow = 3L, byrow = TRUE),
  psd_grid = c(0, 0.5, 1),
  posterior_weights = matrix(
    1 / 3,
    nrow = 3L,
    ncol = 3L,
    dimnames = list(c("G1_rs1", "G2_rs2", "G3_rs3"), NULL)
  ),
  fash_data = list(
    data_list = setNames(lapply(1:3, function(index) {
      data.frame(x = 0:2, y = index + 0:2)
    }), c("G1_rs1", "G2_rs2", "G3_rs3")),
    S = lapply(1:3, function(index) rep(1, 3L)),
    Omega = lapply(1:3, function(index) diag(3L))
  ),
  settings = list(),
  lfdr = rep(1 / 3, 3L)
)
class(toy_fit) <- "fash"
toy_subfit <- build_adjusted_subfit(
  toy_fit,
  selected_indices = c(3L, 1L),
  null_weight = 0.8,
  alternative_weights = c(0.25, 0.75)
)
stopifnot(
  identical(toy_subfit$original_indices, c(1L, 3L)),
  identical(names(toy_subfit$fash_data$data_list), c("G1_rs1", "G3_rs3")),
  max(abs(rowSums(toy_subfit$posterior_weights) - 1)) < 1e-12,
  max(abs(toy_subfit$lfdr - toy_subfit$posterior_weights[, 1L])) < 1e-12
)
toy_existing_subset <- subset_fash_fit(toy_fit, c(2L, 3L))
stopifnot(
  identical(toy_existing_subset$original_indices, c(2L, 3L)),
  identical(
    names(toy_existing_subset$fash_data$data_list),
    c("G2_rs2", "G3_rs3")
  ),
  identical(dim(toy_existing_subset$posterior_weights), c(2L, 3L))
)
toy_nested_subset <- subset_fash_fit(toy_subfit, 2L)
stopifnot(identical(toy_nested_subset$original_indices, 3L))

smooth_var <- seq(0, 15, by = 1)
early_curve <- ifelse(smooth_var <= 3, 2, 0)
middle_curve <- ifelse(smooth_var >= 4 & smooth_var <= 11, 2, 0)
late_curve <- ifelse(smooth_var >= 12, 2, 0)
switch_curve <- seq(-1, 1, length.out = length(smooth_var))
replicate_curve <- function(curve) matrix(rep(curve, 10L), ncol = 10L)
stopifnot(
  classify_functional_draws(replicate_curve(early_curve), smooth_var)["early"] == 0,
  classify_functional_draws(replicate_curve(middle_curve), smooth_var)["middle"] == 0,
  classify_functional_draws(replicate_curve(late_curve), smooth_var)["late"] == 0,
  classify_functional_draws(replicate_curve(switch_curve), smooth_var)["switch"] == 0
)

classification <- data.frame(
  original_index = c(1L, 2L, 3L, 4L),
  key = paste0("G", 1:4, "_rs", 1:4),
  category = rep("early", 4L),
  lfsr = c(0.04, 0.01, 0.20, 0.02),
  stringsAsFactors = FALSE
)
ranked <- add_cumulative_fsr(classification)
stopifnot(
  identical(ranked$key, c("G2_rs2", "G4_rs4", "G1_rs1", "G3_rs3")),
  max(abs(ranked$cfsr - c(0.01, 0.015, 0.07 / 3, 0.27 / 4))) < 1e-12
)

candidates <- data.frame(
  key = c("G1_rs2", "G1_rs1", "G2_rs3", "G3_rs4"),
  gene_id = c("G1", "G1", "G2", "G3"),
  score = c(0.02, 0.01, 0.03, 0.04),
  stringsAsFactors = FALSE
)
selected <- select_distinct_gene_examples(candidates, 2L, "score")
stopifnot(identical(selected$key, c("G1_rs1", "G2_rs3")))

regions <- three_set_venn_regions(
  first = c("a", "b", "d", "g"),
  second = c("b", "c", "d", "f"),
  third = c("d", "e", "f", "g"),
  labels = c("A", "B", "C")
)
expected_counts <- c(1L, 1L, 1L, 1L, 1L, 1L, 1L)
stopifnot(identical(regions$count, expected_counts))

toy_prior_fit <- list(
  psd_grid = c(0, 0.25, 1, 4),
  prior_weights = data.frame(
    psd = c(0, 1, 4),
    prior_weight = c(0.8, 0.15, 0.05)
  )
)
aligned_prior <- extract_full_prior_weights(toy_prior_fit)
stopifnot(
  identical(aligned_prior$psd, toy_prior_fit$psd_grid),
  isTRUE(all.equal(aligned_prior$prior_weight, c(0.8, 0, 0.15, 0.05)))
)

top_indices <- select_top_pair_indices_per_gene(
  gene_id = c("G2", "G1", "G1", "G2", "G3"),
  variant_id = c("rs2", "rs2", "rs1", "rs1", "rs3"),
  lfdr = c(0.2, 0.1, 0.1, 0.2, 0.8)
)
stopifnot(identical(top_indices, c(3L, 4L, 5L)))

term_to_gene <- data.frame(
  term = c("T1", "T1", "T2", "T2"),
  gene = c("G1", "G2", "G3", "G4"),
  stringsAsFactors = FALSE
)
enrichment <- hallmark_hypergeometric(
  selected_genes = c("G1", "G2"),
  universe_genes = paste0("G", 1:4),
  term_to_gene = term_to_gene
)
stopifnot(
  enrichment$overlap_count[enrichment$term == "T1"] == 2L,
  enrichment$overlap_count[enrichment$term == "T2"] == 0L,
  all(is.finite(enrichment$p_value)),
  all(is.finite(enrichment$q_value))
)

cat("FASH-CL manuscript-impact helper tests passed.\n")
