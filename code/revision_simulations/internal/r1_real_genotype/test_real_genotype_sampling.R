#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/internal/r1_real_genotype/real_genotype_sampling.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/internal/r1_real_genotype/real_genotype_sampling.R")) {
    return("coderepo-local")
  }
  stop("Could not find the R1 real-genotype sampling module.")
}

expect_error <- function(expression, pattern) {
  observed <- tryCatch(
    {
      force(expression)
      NULL
    },
    error = function(error) conditionMessage(error)
  )
  if (is.null(observed) || !grepl(pattern, observed, fixed = TRUE)) {
    stop(
      "Expected an error containing '", pattern,
      "' but observed: ", if (is.null(observed)) "no error" else observed
    )
  }
  invisible(observed)
}

write_synthetic_vcf <- function(path, sample_ids, rows) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  writeLines(c(
    "##fileformat=VCFv4.2",
    "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Expected genotype dosage\">",
    paste(
      c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample_ids),
      collapse = "\t"
    ),
    rows
  ), connection)
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "r1_real_genotype",
  "real_genotype_sampling.R"
))

temporary_dir <- tempfile("r1-real-genotype-test-")
dir.create(temporary_dir)
on.exit(unlink(temporary_dir, recursive = TRUE), add = TRUE)

sample_ids <- sprintf("sample_%02d", 1:6)
common_patterns <- list(
  c(0, 0, 1, 1, 2, 2),
  c(0, 1, 0, 1, 2, 2),
  c(0, 0, 0, 1, 1, 2),
  c(0, 1, 1, 1, 2, 2),
  c(0.1, 0.2, 0.9, 1.1, 1.8, 1.9),
  c(0, 0, 1, 1, 1, 1)
)
rare_pattern <- c(0, 0, 0, 0, 0, 1)
monomorphic_pattern <- rep(0, 6)

vcf_rows <- character(0)
pair_ids <- character(0)
position <- 1000L
for (gene_index in 1:4) {
  gene_id <- sprintf("ENSGTEST%d", gene_index)
  for (variant_offset in 1:8) {
    variant_id <- sprintf("rs%d%02d", gene_index, variant_offset)
    dosage <- if (variant_offset <= 6L) {
      common_patterns[[variant_offset]]
    } else if (variant_offset == 7L) {
      rare_pattern
    } else {
      monomorphic_pattern
    }
    vcf_rows <- c(vcf_rows, paste(
      c(
        gene_index,
        position,
        variant_id,
        "A",
        "G",
        "100",
        "PASS",
        ".",
        "DS",
        dosage
      ),
      collapse = "\t"
    ))
    pair_ids <- c(pair_ids, paste(gene_id, variant_id, sep = "_"))
    position <- position + 100L
  }
}
vcf_rows <- c(vcf_rows, paste(
  c(1, position, "rs99999", "C", "T", "100", "PASS", ".", "DS", common_patterns[[1L]]),
  collapse = "\t"
))

vcf_path <- file.path(temporary_dir, "synthetic.vcf.gz")
write_synthetic_vcf(vcf_path, sample_ids, vcf_rows)
stopifnot(identical(read_vcf_sample_ids(vcf_path), sample_ids))

seed_list <- 1:20
source_data <- prepare_real_genotype_source(
  pair_ids = pair_ids,
  vcf_path = vcf_path,
  seed_list = seed_list,
  n_loci = 2,
  variants_per_locus = 3,
  candidate_genes_per_seed = 4,
  minimum_raw_variants = 6,
  maf_min = 0.10,
  work_dir = temporary_dir
)

sampled <- sample_real_genotype_blocks(source_data, seed = seed_list[[1L]])
repeated <- sample_real_genotype_blocks(source_data, seed = seed_list[[1L]])
ld <- summarize_real_genotype_ld(sampled$G, sampled$variant_info)
stopifnot(
  identical(dim(sampled$G), c(6L, 6L)),
  identical(sampled$G, repeated$G),
  identical(sampled$variant_info, repeated$variant_info),
  length(unique(sampled$variant_info$gene_id)) == 2L,
  all(table(sampled$variant_info$gene_id) == 3L),
  all(sampled$variant_info$observed_maf >= 0.10),
  all(apply(sampled$G, 2L, stats::sd) > 0),
  all(is.finite(ld$locus$median_abs_r)),
  nrow(ld$locus) == 2L,
  nrow(ld$seed) == 1L,
  !"rs99999" %in% sampled$variant_info$variant_id,
  !any(sampled$variant_info$variant_id %in% c("rs107", "rs108", "rs207", "rs208"))
)

seed_signatures <- vapply(seed_list, function(seed) {
  candidate <- sample_real_genotype_blocks(source_data, seed = seed)
  paste(candidate$variant_info$unit_id, collapse = "|")
}, character(1))
stopifnot(length(unique(seed_signatures)) > 1L)

expect_error(
  parse_paper_pair_ids(c("ENSGTEST1_rs101", "malformed")),
  "Malformed paper pair identifier"
)
expect_error(
  prepare_real_genotype_source(
    pair_ids = pair_ids,
    vcf_path = vcf_path,
    seed_list = 1L,
    n_loci = 5L,
    variants_per_locus = 3L,
    candidate_genes_per_seed = 4L,
    minimum_raw_variants = 6L,
    maf_min = 0.10,
    work_dir = temporary_dir
  ),
  "5 loci are required"
)

duplicate_vcf_path <- file.path(temporary_dir, "duplicate.vcf.gz")
write_synthetic_vcf(
  duplicate_vcf_path,
  sample_ids,
  c(vcf_rows, vcf_rows[[1L]])
)
expect_error(
  extract_target_vcf_dosages(
    duplicate_vcf_path,
    target_variants = sub("^ENSGTEST1_", "", pair_ids[[1L]]),
    work_dir = temporary_dir
  ),
  "duplicate retained variant ID"
)

cat("R1 real-genotype sampling tests passed.\n")
