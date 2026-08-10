#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath("."))
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath("coderepo-local"))
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "r5_balanced_thinning",
  "balanced_thinning_helpers.R"
))

pair_keys <- c(
  paste0("geneA_v", 1:3),
  paste0("geneB_v", 1:4),
  paste0("geneC_v", 1:5)
)
gene_index <- make_gene_index(pair_keys)
selection <- select_balanced_variants_per_gene(
  pair_keys,
  seed = 12345L,
  target_per_gene = 4L,
  gene_index = gene_index
)
repeated_selection <- select_balanced_variants_per_gene(
  pair_keys,
  seed = 12345L,
  target_per_gene = 4L,
  gene_index = gene_index
)
alternative_selection <- select_balanced_variants_per_gene(
  pair_keys,
  seed = 22345L,
  target_per_gene = 4L,
  gene_index = gene_index
)

stopifnot(
  identical(selection, repeated_selection),
  !identical(selection$pair_key, alternative_selection$pair_key),
  identical(sort(unique(selection$gene_id)), c("geneB", "geneC")),
  all(table(selection$gene_id) == 4L),
  nrow(selection) == 8L,
  anyDuplicated(selection$fash_index) == 0L,
  anyDuplicated(selection$pair_key) == 0L,
  all(selection$target_per_gene == 4L)
)

prior_comparison <- compare_prior_weights(
  data.frame(psd = c(0, 1), prior_weight = c(0.8, 0.2)),
  data.frame(psd = c(0, 0.5, 1), prior_weight = c(0.75, 0.05, 0.2))
)
stopifnot(
  isTRUE(all.equal(prior_comparison$summary$full_pi0, 0.8)),
  isTRUE(all.equal(prior_comparison$summary$thinned_pi0, 0.75)),
  isTRUE(all.equal(prior_comparison$summary$prior_total_variation, 0.05))
)

lfdr_comparison <- compare_paired_lfdr(
  full_lfdr = c(0.001, 0.01, 0.9, 0.95),
  thinned_lfdr = c(0.002, 0.02, 0.85, 0.96),
  pair_keys = paste0("gene", 1:4, "_v1"),
  alpha = 0.05
)
stopifnot(
  nrow(lfdr_comparison$table) == 4L,
  lfdr_comparison$summary$full_fdr_calls == 2L,
  lfdr_comparison$summary$thinned_fdr_calls == 2L,
  lfdr_comparison$summary$fdr_call_jaccard == 1
)

gene_sets <- list(c("g1", "g2"), c("g2", "g3"), c("g1", "g4"))
curve <- cumulative_gene_union(
  gene_sets,
  seeds = c(1L, 2L, 3L),
  full_discovered_genes = c("g1", "g2", "g3")
)
envelope <- seed_subset_union_envelope(
  gene_sets,
  seeds = c(1L, 2L, 3L),
  full_discovered_genes = c("g1", "g2", "g3")
)
stopifnot(
  identical(curve$cumulative_unique_genes, c(2L, 3L, 4L)),
  identical(curve$cumulative_full_genes_recovered, c(2L, 3L, 3L)),
  all(diff(curve$cumulative_unique_genes) >= 0L),
  nrow(envelope) == 3L,
  envelope$unique_min[1] == 2L,
  envelope$unique_max[3] == 4L
)

cat("R5 balanced-thinning helper tests passed.\n")
