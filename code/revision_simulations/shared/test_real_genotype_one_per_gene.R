#!/usr/bin/env Rscript

# Focused tests for deterministic one-variant-per-gene real-genotype sampling.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/real_genotype_one_per_gene.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/real_genotype_one_per_gene.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the real-genotype one-per-gene helper.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
))

pair_ids <- c(
  "geneA_rs1", "geneA_rs2", "geneB_rs2", "geneB_rs3",
  "geneC_rs4", "geneC_rs5"
)
parsed <- parse_tested_pair_ids(pair_ids)
stopifnot(
  identical(parsed$gene_id, c("geneA", "geneA", "geneB", "geneB", "geneC", "geneC")),
  identical(parsed$variant_id, c("rs1", "rs2", "rs2", "rs3", "rs4", "rs5")),
  identical(parsed$pair_key, pair_ids)
)

selection <- sample_one_tested_variant_per_gene(pair_ids, seed = 12345L)
repeated_selection <- sample_one_tested_variant_per_gene(
  pair_ids,
  seed = 12345L
)
alternative_selection <- sample_one_tested_variant_per_gene(
  pair_ids,
  seed = 22345L
)
stopifnot(
  nrow(selection) == 3L,
  identical(selection$gene_id, c("geneA", "geneB", "geneC")),
  !anyDuplicated(selection$gene_id),
  identical(selection, repeated_selection),
  !identical(selection$pair_key, alternative_selection$pair_key),
  all(selection$pair_key %in% pair_ids)
)

extracted <- list(
  sample_ids = c("donor1", "donor2", "donor3", "donor4"),
  metadata = data.frame(
    chromosome = c("1", "1", "2", "2", "3"),
    position = c(10L, 20L, 30L, 40L, 50L),
    variant_id = c("rs1", "rs2", "rs3", "rs4", "rs5"),
    ref = rep("A", 5L),
    alt = rep("G", 5L),
    stringsAsFactors = FALSE
  ),
  dosage = rbind(
    rs1 = c(0, 0, 1, 1),
    rs2 = c(0, 1, 1, 2),
    rs3 = c(0, 1, 2, 2),
    rs4 = c(0, 0, 0, 1),
    rs5 = c(0, 1, 1, 1)
  ),
  requested_variant_count = 5L,
  matched_variant_count = 5L
)
colnames(extracted$dosage) <- extracted$sample_ids

duplicate_selection <- data.frame(
  gene_id = c("geneA", "geneB", "geneC"),
  variant_id = c("rs2", "rs2", "rs5"),
  pair_key = c("geneA_rs2", "geneB_rs2", "geneC_rs5"),
  stringsAsFactors = FALSE
)
sample <- assemble_one_per_gene_genotype_sample(
  selection = duplicate_selection,
  extracted = extracted,
  seed = 12345L
)
validated <- validate_real_genotype_sample(
  sample,
  expected_genes = 3L,
  expected_donors = 4L,
  maf_min = 0
)
stopifnot(
  identical(dim(validated$G), c(4L, 3L)),
  identical(rownames(validated$G), extracted$sample_ids),
  identical(colnames(validated$G), duplicate_selection$pair_key),
  identical(validated$G[, 1L], validated$G[, 2L]),
  identical(validated$variant_info$gene_id, duplicate_selection$gene_id),
  identical(validated$variant_info$variant_id, duplicate_selection$variant_id),
  identical(validated$variant_info$observed_maf, c(0.5, 0.5, 0.375)),
  validated$selection_summary$unique_variant_ids == 2L,
  validated$selection_summary$repeated_cross_gene_assignments == 1L
)

digest_G <- matrix(
  c(0, 1, 2, 2, 1, 0),
  nrow = 2L,
  dimnames = list(
    c("donor1", "donor2"),
    c("geneA_rs1", "geneB_rs2", "geneC_rs3")
  )
)
digest_pair_key <- colnames(digest_G)
digest_sample_ids <- rownames(digest_G)
expected_content_digest <- "2134bb6e90380af33cecf7060ca43b94"
observed_content_digest <- genotype_content_md5(
  pair_key = digest_pair_key,
  sample_ids = digest_sample_ids,
  G = digest_G
)
changed_G <- digest_G
changed_G[1L, 1L] <- 1
stopifnot(
  identical(observed_content_digest, expected_content_digest),
  identical(
    observed_content_digest,
    genotype_content_md5(digest_pair_key, digest_sample_ids, digest_G)
  ),
  !identical(
    observed_content_digest,
    genotype_content_md5(
      c("geneA_rs9", digest_pair_key[-1L]),
      digest_sample_ids,
      digest_G
    )
  ),
  !identical(
    observed_content_digest,
    genotype_content_md5(
      digest_pair_key,
      c("donor9", digest_sample_ids[-1L]),
      digest_G
    )
  ),
  !identical(
    observed_content_digest,
    genotype_content_md5(digest_pair_key, digest_sample_ids, changed_G)
  )
)
nonfinite_G <- digest_G
nonfinite_G[1L, 1L] <- NA_real_
stopifnot(
  inherits(try(
    genotype_content_md5(
      digest_pair_key[-1L],
      digest_sample_ids,
      digest_G
    ),
    silent = TRUE
  ), "try-error"),
  inherits(try(
    genotype_content_md5(
      digest_pair_key,
      digest_sample_ids,
      nonfinite_G
    ),
    silent = TRUE
  ), "try-error"),
  inherits(try(
    genotype_content_md5(
      digest_pair_key,
      digest_sample_ids,
      matrix(as.character(digest_G), nrow = nrow(digest_G))
    ),
    silent = TRUE
  ), "try-error")
)

invalid <- sample
invalid$G[, 1L] <- 0
stopifnot(inherits(try(
  validate_real_genotype_sample(
    invalid,
    expected_genes = 3L,
    expected_donors = 4L,
    maf_min = 0
  ),
  silent = TRUE
), "try-error"))

maf <- c(0.11, 0.12, 0.14, 0.18, 0.21, 0.24, 0.27, 0.31, 0.36, 0.42,
         0.13, 0.16, 0.19, 0.22, 0.25, 0.29, 0.33, 0.38, 0.45, 0.49)
class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
targets <- make_maf_balanced_class_targets(
  maf = maf,
  class_probs = class_probs,
  seed = 12345L,
  n_strata = 4L
)
stopifnot(
  identical(
    as.integer(table(factor(targets, levels = names(class_probs)))),
    as.integer(exact_proportional_counts(length(maf), class_probs))
  ),
  identical(
    targets,
    make_maf_balanced_class_targets(maf, class_probs, 12345L, 4L)
  )
)

set.seed(91)
effect_sim <- simulate_variant_effect_curves(
  n_variants = length(maf),
  time_grid = 0:15,
  class_probs = class_probs,
  dynamic_amplitude = 2,
  dynamic_main_effect_sd = 1,
  exact_class_counts = TRUE,
  seed = 91L
)
original_rows <- unname(sort(apply(
  effect_sim$beta_matrix,
  1L,
  paste,
  collapse = "\r"
)))
reassigned <- reassign_effect_simulation_by_maf(
  effect_sim = effect_sim,
  maf = maf,
  class_probs = class_probs,
  seed = 12345L,
  n_strata = 4L
)
stopifnot(
  identical(reassigned$unit_info$effect_class, targets),
  identical(
    unname(sort(apply(
      reassigned$beta_matrix,
      1L,
      paste,
      collapse = "\r"
    ))),
    original_rows
  ),
  identical(reassigned$unit_info$unit_index, seq_along(maf)),
  !anyDuplicated(reassigned$unit_info$unit_id),
  identical(rownames(reassigned$beta_matrix), reassigned$unit_info$variant_id)
)

balance <- summarize_truth_maf_balance(
  variant_info = data.frame(observed_maf = maf),
  unit_info = reassigned$unit_info,
  seed = 12345L
)
stopifnot(
  nrow(balance) == length(class_probs),
  setequal(balance$effect_class, names(class_probs)),
  all(
    balance$n == exact_proportional_counts(
      length(maf),
      class_probs
    )[balance$effect_class]
  )
)

message("Real-genotype one-per-gene tests passed.")
