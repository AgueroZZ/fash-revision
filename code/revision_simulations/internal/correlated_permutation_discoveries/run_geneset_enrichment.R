#!/usr/bin/env Rscript

# Hallmark gene-set enrichment for the three discovery sets produced by
# run_three_likelihoods.R, using the same `enrich_set` construction as
# analysis/dynamic_eQTL_real_full.rmd: the universe is forced, and background
# genes with no Hallmark annotation are carried by a dummy pathway so that the
# hypergeometric universe is exactly the intended one.
#
# Two universes, answering two different questions:
#   all tested genes   is each set enriched at all?
#   the 1,169 panel    given the diagonal fit already selected these, does the
#                      correlated likelihood's subselection pick out the more
#                      biologically coherent ones? This one controls for the
#                      outcome-dependent selection of the panel.

suppressPackageStartupMessages({
  library(dplyr); library(msigdbr); library(clusterProfiler)
})
cache <- "output/revision_simulations/internal/correlated_permutation_discoveries"
alpha <- 0.05

u <- utils::read.csv(file.path(cache, "three_likelihood_units.csv"),
                     stringsAsFactors = FALSE)
call_set <- function(lfdr, keys) {
  o <- order(lfdr, method = "radix")
  keys[o[cumsum(lfdr[o]) / seq_along(o) <= alpha]]
}
sets <- list(
  diagonal = call_set(u$diagonal_bf_lfdr, u$pair_key),
  common   = call_set(u$common_bf_lfdr,   u$pair_key),
  per_unit = call_set(u$per_unit_bf_lfdr, u$pair_key)
)
genes <- lapply(sets, function(k) unique(u$gene_id[match(k, u$pair_key)]))
cat("set sizes (genes):", paste(names(genes), lengths(genes), sep = "="), "\n")
cat("nested? per_unit in common:", all(genes$per_unit %in% genes$common),
    " common in diagonal:", all(genes$common %in% genes$diagonal), "\n\n")

m_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  mutate(ensembl_use = dplyr::coalesce(ensembl_gene, db_ensembl_gene)) %>%
  filter(!is.na(ensembl_use)) %>% select(gs_name, ensembl_use) %>% distinct()
hallmark_genes <- unique(m_t2g$ensembl_use)

enrich_set <- function(selected, universe) {
  selected <- intersect(unique(selected), universe)
  missing <- setdiff(universe, hallmark_genes)
  t2g <- if (length(missing)) {
    bind_rows(m_t2g, tibble::tibble(gs_name = "__DUMMY__", ensembl_use = missing))
  } else m_t2g
  res <- enricher(gene = selected, TERM2GENE = t2g, universe = universe,
                  pAdjustMethod = "BH", qvalueCutoff = 1, pvalueCutoff = 1)
  if (is.null(res) || nrow(res@result) == 0L) return(NULL)
  r <- res@result %>% filter(ID != "__DUMMY__")
  if (nrow(r) == 0L) return(NULL)
  num <- function(x) as.numeric(sub("/.*$", "", x))
  den <- function(x) as.numeric(sub("^.*/", "", x))
  r %>% transmute(
    term = sub("^HALLMARK_", "", Description),
    hits = num(GeneRatio), set_size = den(GeneRatio),
    bg_hits = num(BgRatio), bg_size = den(BgRatio),
    fold = (num(GeneRatio)/den(GeneRatio)) / (num(BgRatio)/den(BgRatio)),
    pvalue, q = p.adjust
  ) %>% arrange(pvalue)
}

universes <- list(
  all_tested = utils::read.csv(paste0("output/revision_simulations/internal/",
    "residual_correlation_fash/all_tested_genes.csv"),
    stringsAsFactors = FALSE)$gene_id,
  panel_1169 = unique(u$gene_id)
)

all_rows <- list(); summary_rows <- list()
for (un in names(universes)) {
  U <- universes[[un]]
  for (nm in names(genes)) {
    if (un == "panel_1169" && nm == "diagonal") next   # diagonal is the panel
    r <- enrich_set(genes[[nm]], U)
    n_in_hallmark <- length(intersect(genes[[nm]], hallmark_genes))
    if (is.null(r)) {
      summary_rows[[paste(un, nm)]] <- data.frame(
        universe = un, set = nm, n_genes = length(genes[[nm]]),
        n_in_hallmark = n_in_hallmark, terms_tested = 0L,
        terms_q05 = 0L, best_p = NA_real_, best_q = NA_real_,
        median_fold = NA_real_, stringsAsFactors = FALSE)
      next
    }
    r$universe <- un; r$set <- nm
    all_rows[[paste(un, nm)]] <- r
    summary_rows[[paste(un, nm)]] <- data.frame(
      universe = un, set = nm, n_genes = length(genes[[nm]]),
      n_in_hallmark = n_in_hallmark, terms_tested = nrow(r),
      terms_q05 = sum(r$q < 0.05), best_p = min(r$pvalue), best_q = min(r$q),
      median_fold = median(r$fold), stringsAsFactors = FALSE)
  }
}
summary_rows <- do.call(rbind, summary_rows)
all_rows <- do.call(rbind, all_rows)
options(width = 200)
cat("=== summary ===\n"); print(summary_rows, digits = 3, row.names = FALSE)

cat("\n=== annotation density: share of each set's genes carrying any Hallmark annotation ===\n")
bg_rate <- mean(universes$all_tested %in% hallmark_genes)
dens <- do.call(rbind, lapply(names(genes), function(nm) {
  k <- sum(genes[[nm]] %in% hallmark_genes); n <- length(genes[[nm]])
  tt <- matrix(c(k, n - k, round(bg_rate*length(universes$all_tested)),
                 round((1-bg_rate)*length(universes$all_tested))), nrow = 2)
  data.frame(set = nm, n = n, in_hallmark = k, rate = k/n,
             background_rate = bg_rate, fold = (k/n)/bg_rate,
             fisher_p = fisher.test(tt)$p.value, stringsAsFactors = FALSE)
}))
print(dens, digits = 3, row.names = FALSE)

cat("\n=== fold enrichment on a common set of terms (size-robust comparison) ===\n")
cat("    terms reaching q < 0.05 in any set, universe = all tested genes\n")
aa <- all_rows[all_rows$universe == "all_tested", ]
keep <- unique(aa$term[aa$q < 0.05])
if (length(keep) == 0L) {
  cat("    no term reached q < 0.05 in any set\n")
} else {
  w <- reshape(aa[aa$term %in% keep, c("term","set","fold")],
               idvar = "term", timevar = "set", direction = "wide")
  names(w) <- sub("fold.", "fold_", names(w))
  wq <- reshape(aa[aa$term %in% keep, c("term","set","q")],
                idvar = "term", timevar = "set", direction = "wide")
  names(wq) <- sub("q.", "q_", names(wq))
  print(merge(w, wq, by = "term"), digits = 3, row.names = FALSE)
}

utils::write.csv(summary_rows, file.path(cache, "geneset_summary.csv"), row.names = FALSE)
utils::write.csv(all_rows, file.path(cache, "geneset_terms.csv"), row.names = FALSE)
utils::write.csv(dens, file.path(cache, "geneset_annotation_density.csv"), row.names = FALSE)
