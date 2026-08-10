#!/usr/bin/env Rscript

# Focused tests for the internal variant annotation enrichment analysis.

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
  "variant_annotation_enrichment",
  "variant_annotation_enrichment_helpers.R"
)
source(helper_path)

stopifnot(identical(
  normalize_autosome(c("chr1", "X", "chr22", "chrM")),
  c("1", NA_character_, "22", NA_character_)
))

stopifnot(identical(
  extract_gtf_gene_id(c(
    'gene_id "ENSG000001.12"; gene_type "protein_coding";',
    'transcript_id "ENST000001.1";'
  )),
  c("ENSG000001", NA_character_)
))

keys <- c("ENSG000001_rs1", "ENSG000002_rs2")
parsed <- parse_pair_keys(keys)
stopifnot(
  identical(parsed$gene_id, c("ENSG000001", "ENSG000002")),
  identical(parsed$variant_id, c("rs1", "rs2"))
)

calls <- cumulative_lfdr_calls(c(0.01, 0.04, 0.20), alpha = 0.05)
stopifnot(identical(calls, c(1L, 2L)))

temporary_directory <- tempfile("variant-enrichment-tests-")
dir.create(temporary_directory)
on.exit(unlink(temporary_directory, recursive = TRUE), add = TRUE)

vcf_path <- file.path(temporary_directory, "toy.vcf.gz")
vcf_connection <- gzfile(vcf_path, open = "wt")
writeLines(c(
  "##fileformat=VCFv4.2",
  "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Expected dosage\">",
  paste(
    c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", "sample1", "sample2", "sample3", "sample4"),
    collapse = "\t"
  ),
  paste(c("1", "100", "rs1", "A", "G", "100", "PASS", ".", "DS",
          "0", "1", "2", "1"), collapse = "\t"),
  paste(c("1", "200", "rs2", "C", "T", "100", "PASS", ".", "DS",
          "2", "2", "1", "1"), collapse = "\t"),
  paste(c("2", "300", "rs3", "G", "A", "100", "PASS", ".", "DS",
          "0", ".", "0", "1"), collapse = "\t"),
  paste(c("2", "400", "rs4", "T", "C", "100", "PASS", ".", "DS",
          "0", "0", "0", "0"), collapse = "\t")
), vcf_connection)
close(vcf_connection)

vcf_result <- read_vcf_variant_metadata(
  vcf_path = vcf_path,
  keep_ids = c("rs1", "rs2", "rs3"),
  retain_dosage_ids = c("rs1", "rs2")
)
stopifnot(
  nrow(vcf_result$metadata) == 3L,
  identical(vcf_result$metadata$variant_id, c("rs1", "rs2", "rs3")),
  isTRUE(all.equal(vcf_result$metadata$minor_allele_frequency,
                   c(0.5, 0.25, 1 / 6))),
  identical(dim(vcf_result$dosage), c(2L, 4L)),
  identical(rownames(vcf_result$dosage), c("rs1", "rs2"))
)

variants <- data.frame(
  variant_id = c("rs1", "rs2", "rs3"),
  chromosome = c("1", "1", "2"),
  position = c(100L, 250L, 50L),
  stringsAsFactors = FALSE
)
intervals <- data.frame(
  chromosome = c("1", "2"),
  start = c(90L, 1L),
  end = c(110L, 100L),
  annotation = c("promoter", "enhancer"),
  stringsAsFactors = FALSE
)
overlaps <- annotate_variant_overlaps(variants, intervals)
stopifnot(
  overlaps$promoter[overlaps$variant_id == "rs1"],
  overlaps$enhancer[overlaps$variant_id == "rs3"],
  !overlaps$promoter[overlaps$variant_id == "rs2"]
)

matching_universe <- data.frame(
  variant_id = paste0("rs", seq_len(18L)),
  chromosome = rep(c("1", "2"), each = 9L),
  minor_allele_frequency = rep(c(0.10, 0.11, 0.12), 6L),
  minimum_target_tss_distance = rep(c(1000, 2000, 3000), 6L),
  local_tested_variant_count_1mb = rep(c(10, 12, 14), 6L),
  n_tested_genes = 1L,
  stringsAsFactors = FALSE
)
matching_universe <- derive_matching_strata(matching_universe, n_bins = 3L)
toy_matches_a <- sample_matched_controls(
  matching_universe,
  selected_ids = c("rs1", "rs10"),
  controls_per_variant = 5L,
  seed = 123L
)
toy_matches_b <- sample_matched_controls(
  matching_universe,
  selected_ids = c("rs1", "rs10"),
  controls_per_variant = 5L,
  seed = 123L
)
stopifnot(
  nrow(toy_matches_a) == 10L,
  identical(toy_matches_a, toy_matches_b),
  all(toy_matches_a$selected_id != toy_matches_a$control_id),
  all(toy_matches_a$relaxation_level %in% 0:3)
)

streaming_annotations <- data.frame(
  variant_id = matching_universe$variant_id,
  chromosome = matching_universe$chromosome,
  ubiquitous = rep(TRUE, nrow(matching_universe)),
  stringsAsFactors = FALSE
)
streaming_result <- run_streaming_matched_enrichment(
  variant_table = matching_universe,
  annotation_matrix = streaming_annotations,
  selected_sets = list(toy = c("rs1", "rs10")),
  seeds = c(11L, 12L),
  controls_per_variant = 2L,
  chromosomes = c("1", "2"),
  minimum_overlap_count = 1L
)
stopifnot(
  nrow(streaming_result$results) == 1L,
  isTRUE(all.equal(streaming_result$results$log2_enrichment, 0)),
  nrow(streaming_result$matching_balance) == 8L,
  sum(streaming_result$relaxation_audit$Freq) == 2L
)

fash_variants <- data.frame(
  variant_id = c("f1", "f2", "f3"),
  chromosome = c("1", "1", "1"),
  position = c(100L, 1000L, 3000000L),
  stringsAsFactors = FALSE
)
strober_variants <- data.frame(
  variant_id = c("s1", "s2", "s3"),
  chromosome = c("1", "1", "1"),
  position = c(150L, 1100L, 4500000L),
  stringsAsFactors = FALSE
)
dosage_matrix <- rbind(
  f1 = c(0, 0, 1, 1, 2, 2),
  f2 = c(0, 1, 0, 1, 0, 1),
  f3 = c(0, 0, 1, 1, 2, 2),
  s1 = c(0, 0, 1, 1, 2, 2),
  s2 = c(0, 0, 0, 1, 1, 1),
  s3 = c(0, 0, 1, 1, 2, 2)
)
ld_classification <- classify_cross_method_ld(
  fash_variants = fash_variants,
  strober_variants = strober_variants,
  dosage_matrix = dosage_matrix,
  window_bp = 1000000L,
  r2_threshold = 0.8
)
stopifnot(
  ld_classification$fash$ld_shared[ld_classification$fash$variant_id == "f1"],
  !ld_classification$fash$ld_shared[ld_classification$fash$variant_id == "f2"],
  !ld_classification$fash$ld_shared[ld_classification$fash$variant_id == "f3"]
)

annotation_matrix <- data.frame(
  variant_id = paste0("v", 1:12),
  chromosome = rep(c("1", "2"), each = 6L),
  active = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE,
             TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)
match_table <- data.frame(
  discovery_set = "toy",
  match_seed = 1L,
  selected_id = c("v1", "v2", "v7", "v8"),
  selected_chromosome = c("1", "1", "2", "2"),
  control_id = c("v3", "v4", "v9", "v10"),
  stringsAsFactors = FALSE
)
point_result <- compute_matched_enrichment(
  annotation_matrix,
  match_table,
  minimum_overlap_count = 1L
)
jackknife_result <- chromosome_jackknife_enrichment(
  annotation_matrix,
  match_table,
  chromosomes = c("1", "2"),
  minimum_overlap_count = 1L
)
stopifnot(
  nrow(point_result) == 1L,
  is.finite(point_result$log2_enrichment),
  nrow(jackknife_result) == 1L,
  is.finite(jackknife_result$jackknife_se),
  jackknife_result$ci_lower <= point_result$log2_enrichment,
  jackknife_result$ci_upper >= point_result$log2_enrichment
)

message("Variant annotation enrichment helper tests passed.")
