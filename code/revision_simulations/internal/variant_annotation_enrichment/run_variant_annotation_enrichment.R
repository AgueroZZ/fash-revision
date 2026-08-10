#!/usr/bin/env Rscript

# Build retained caches for the internal FASH variant annotation enrichment page.

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

require_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The ", package, " package is required.")
  }
  invisible(TRUE)
}

file_provenance <- function(path) {
  data.frame(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    byte_size = unname(file.info(path)$size),
    md5 = unname(tools::md5sum(path)),
    modified_at = format(file.info(path)$mtime, tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}

rank_unique_variants <- function(variant_id, score) {
  score_table <- data.table::data.table(
    variant_id = as.character(variant_id),
    score = as.numeric(score)
  )
  score_table <- score_table[, .(score = min(score)), by = variant_id]
  data.table::setorder(score_table, score, variant_id)
  score_table
}

lead_variant_per_gene <- function(gene_id, variant_id, score) {
  table <- data.frame(
    gene_id = as.character(gene_id),
    variant_id = as.character(variant_id),
    score = as.numeric(score),
    stringsAsFactors = FALSE
  )
  table <- table[order(table$score, table$variant_id, method = "radix"), ]
  unique(table$variant_id[!duplicated(table$gene_id)])
}

append_annotation_columns <- function(annotation_matrix, overlap_table) {
  if (!identical(annotation_matrix$variant_id, overlap_table$variant_id)) {
    stop("Annotation overlap rows are not aligned to the variant universe.")
  }
  new_columns <- setdiff(
    names(overlap_table),
    c("variant_id", "chromosome")
  )
  if (any(new_columns %in% names(annotation_matrix))) {
    stop("Annotation names are duplicated.")
  }
  for (column in new_columns) {
    annotation_matrix[[column]] <- as.logical(overlap_table[[column]])
  }
  annotation_matrix
}

message_step <- function(index, total, text) {
  message("[", index, "/", total, "] ", text)
}

workflowr_root <- find_workflowr_root()
require_package("data.table")
require_package("GenomicRanges")
require_package("IRanges")

helper_path <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "variant_annotation_enrichment",
  "variant_annotation_enrichment_helpers.R"
)
source(helper_path)

fit_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "fash_fit1_update.RData"
)
linear_path <- file.path(
  workflowr_root,
  "data",
  "dynamic_eQTL_real",
  "strober_linear",
  "linear_dynamic_eqtls_5_pc.txt"
)
nonlinear_path <- file.path(
  workflowr_root,
  "data",
  "dynamic_eQTL_real",
  "strober_nonlinear",
  "non_linear_dynamic_eqtls_5_pc.txt"
)
vcf_path <- file.path(
  dirname(workflowr_root),
  "iPSC-data",
  "genotype-data",
  "YRI_genotype.vcf.gz"
)
annotation_directory <- file.path(
  workflowr_root,
  "data",
  "revision_simulations",
  "internal",
  "variant_annotation_enrichment"
)
output_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "variant_annotation_enrichment"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

input_paths <- c(fit_path, linear_path, nonlinear_path, vcf_path)
annotation_paths <- c(
  file.path(annotation_directory, "gencode.v19.annotation.gtf.gz"),
  file.path(annotation_directory, "ENCFF788SJC.bed.gz"),
  file.path(annotation_directory, "E020_15_coreMarks_mnemonics.bed.gz"),
  file.path(annotation_directory, "E013_15_coreMarks_mnemonics.bed.gz"),
  file.path(annotation_directory, "E095_15_coreMarks_mnemonics.bed.gz"),
  file.path(annotation_directory, "download_manifest.csv")
)
if (any(!file.exists(c(input_paths, annotation_paths)))) {
  stop("One or more required input or annotation files are missing.")
}

expected <- list(
  tested_pairs = 1009173L,
  tested_variants = 745867L,
  tested_genes = 6362L,
  fash_pairs = 9205L,
  fash_variants = 9139L,
  fash_genes = 1177L,
  linear_pairs = 5404L,
  linear_variants = 5387L,
  nonlinear_pairs = 6824L,
  nonlinear_variants = 6797L
)
alpha <- 0.05
matching_seeds <- seq.int(20260807L, length.out = 100L)
controls_per_variant <- 5L
ld_window_bp <- 1000000L
ld_r2_threshold <- 0.8

message_step(1, 8, "Extracting FASH pair keys, lfdr values, and discoveries.")
fit_environment <- new.env(parent = emptyenv())
loaded_names <- load(fit_path, envir = fit_environment)
if (!identical(loaded_names, "fash_fit1_update")) {
  stop("The FASH fit file did not contain only fash_fit1_update.")
}
fit <- fit_environment$fash_fit1_update
pair_keys <- names(fit$fash_data$data_list)
lfdr <- as.numeric(fit$lfdr)
if (length(pair_keys) != expected$tested_pairs ||
    length(lfdr) != expected$tested_pairs || anyNA(pair_keys) ||
    any(!is.finite(lfdr))) {
  stop("The FASH fit dimensions do not match the fixed analysis invariants.")
}
pair_table <- parse_pair_keys(pair_keys)
if (length(unique(pair_table$variant_id)) != expected$tested_variants ||
    length(unique(pair_table$gene_id)) != expected$tested_genes) {
  stop("The tested FASH universe does not match the expected dimensions.")
}
fash_call_indices <- cumulative_lfdr_calls(lfdr, alpha = alpha)
fash_discovered_pairs <- pair_table[fash_call_indices, , drop = FALSE]
if (length(fash_call_indices) != expected$fash_pairs ||
    length(unique(fash_discovered_pairs$variant_id)) != expected$fash_variants ||
    length(unique(fash_discovered_pairs$gene_id)) != expected$fash_genes) {
  stop("The BF-adjusted FASH discovery counts changed unexpectedly.")
}
fash_variants <- unique(fash_discovered_pairs$variant_id)
fash_rank <- rank_unique_variants(pair_table$variant_id, lfdr)
fash_lead_variants <- lead_variant_per_gene(
  fash_discovered_pairs$gene_id,
  fash_discovered_pairs$variant_id,
  lfdr[fash_call_indices]
)
rm(fit, fit_environment)
invisible(gc())

message_step(2, 8, "Reading and validating Strober linear and nonlinear results.")
read_strober_result <- function(path, expected_pairs, expected_variants) {
  result <- data.table::fread(path, data.table = FALSE, showProgress = FALSE)
  required_columns <- c("rs_id", "ensamble_id", "pvalue", "eFDR")
  if (!identical(names(result), required_columns) ||
      nrow(result) != expected$tested_pairs || anyDuplicated(paste(
        result$ensamble_id, result$rs_id, sep = "_"
      ))) {
    stop("A Strober result table has an unexpected schema or pair universe.")
  }
  result_key <- paste(result$ensamble_id, result$rs_id, sep = "_")
  if (anyNA(match(result_key, pair_keys))) {
    stop("A Strober result table is not aligned to the FASH pair universe.")
  }
  discovered <- result[result$eFDR <= alpha, , drop = FALSE]
  if (nrow(discovered) != expected_pairs ||
      length(unique(discovered$rs_id)) != expected_variants) {
    stop("A Strober discovery count changed unexpectedly.")
  }
  list(
    variants = unique(discovered$rs_id),
    rank = rank_unique_variants(result$rs_id, result$pvalue),
    lead_variants = lead_variant_per_gene(
      discovered$ensamble_id,
      discovered$rs_id,
      discovered$pvalue
    ),
    discovered_pairs = data.frame(
      gene_id = discovered$ensamble_id,
      variant_id = discovered$rs_id,
      score = discovered$pvalue,
      stringsAsFactors = FALSE
    )
  )
}

linear <- read_strober_result(
  linear_path,
  expected_pairs = expected$linear_pairs,
  expected_variants = expected$linear_variants
)
nonlinear <- read_strober_result(
  nonlinear_path,
  expected_pairs = expected$nonlinear_pairs,
  expected_variants = expected$nonlinear_variants
)
linear_variants <- linear$variants
nonlinear_variants <- nonlinear$variants
strober_union_variants <- union(linear_variants, nonlinear_variants)
equal_k <- min(
  length(fash_variants),
  length(linear_variants),
  length(nonlinear_variants)
)
equal_k_sets <- list(
  equal_k_fash = head(fash_rank$variant_id, equal_k),
  equal_k_linear = head(linear$rank$variant_id, equal_k),
  equal_k_nonlinear = head(nonlinear$rank$variant_id, equal_k)
)
linear_lead_variants <- linear$lead_variants
nonlinear_lead_variants <- nonlinear$lead_variants
rm(linear, nonlinear, fash_rank)
invisible(gc())

message_step(3, 8, "Reading study VCF coordinates, cohort MAF, and discovery dosages.")
tested_variants <- unique(pair_table$variant_id)
ld_variant_ids <- Reduce(
  union,
  list(fash_variants, linear_variants, nonlinear_variants)
)
vcf <- read_vcf_variant_metadata(
  vcf_path = vcf_path,
  keep_ids = tested_variants,
  retain_dosage_ids = ld_variant_ids
)
if (nrow(vcf$metadata) != expected$tested_variants ||
    !setequal(vcf$metadata$variant_id, tested_variants) ||
    nrow(vcf$dosage) != length(ld_variant_ids)) {
  stop("The study VCF does not contain the complete tested variant universe.")
}

message_step(4, 8, "Deriving target-distance and local-density matching covariates.")
gencode <- read_gencode_v19_annotations(annotation_paths[1L])
gene_tss_index <- match(pair_table$gene_id, gencode$gene_tss$gene_id)
variant_index <- match(pair_table$variant_id, vcf$metadata$variant_id)
pair_covariates <- data.table::data.table(
  variant_id = pair_table$variant_id,
  gene_id = pair_table$gene_id,
  variant_chromosome = vcf$metadata$chromosome[variant_index],
  variant_position = vcf$metadata$position[variant_index],
  gene_chromosome = gencode$gene_tss$chromosome[gene_tss_index],
  gene_tss = gencode$gene_tss$tss[gene_tss_index]
)
pair_covariates[, valid_distance :=
  !is.na(variant_chromosome) & !is.na(gene_chromosome) &
    is.finite(variant_position) & is.finite(gene_tss) &
    variant_chromosome == gene_chromosome]
variant_aggregates <- pair_covariates[, .(
  minimum_target_tss_distance = if (any(valid_distance)) {
    min(abs(variant_position[valid_distance] - gene_tss[valid_distance]))
  } else {
    NA_real_
  },
  n_tested_genes = data.table::uniqueN(gene_id)
), by = variant_id]
aggregate_index <- match(vcf$metadata$variant_id, variant_aggregates$variant_id)
variant_table <- data.frame(
  vcf$metadata,
  minimum_target_tss_distance =
    variant_aggregates$minimum_target_tss_distance[aggregate_index],
  n_tested_genes = variant_aggregates$n_tested_genes[aggregate_index],
  stringsAsFactors = FALSE,
  check.names = FALSE
)
variant_table$local_tested_variant_count_1mb <-
  compute_local_variant_density(variant_table, window_bp = 1000000L)
eligible <- !is.na(variant_table$chromosome) &
  is.finite(variant_table$position) &
  is.finite(variant_table$minor_allele_frequency) &
  is.finite(variant_table$minimum_target_tss_distance) &
  is.finite(variant_table$n_tested_genes)
variant_table <- variant_table[eligible, , drop = FALSE]
row.names(variant_table) <- NULL
variant_table <- derive_matching_strata(variant_table, n_bins = 10L)
rm(pair_covariates, variant_aggregates, gene_tss_index, variant_index,
   aggregate_index)
invisible(gc())

retain_eligible <- function(ids) {
  intersect(as.character(ids), variant_table$variant_id)
}
fash_variants <- retain_eligible(fash_variants)
linear_variants <- retain_eligible(linear_variants)
nonlinear_variants <- retain_eligible(nonlinear_variants)
strober_union_variants <- retain_eligible(strober_union_variants)
equal_k_sets <- lapply(equal_k_sets, retain_eligible)
lead_sets <- list(
  one_lead_fash = retain_eligible(fash_lead_variants),
  one_lead_linear = retain_eligible(linear_lead_variants),
  one_lead_nonlinear = retain_eligible(nonlinear_lead_variants)
)

message_step(5, 8, "Annotating eligible variants with GENCODE, ENCODE, and Roadmap states.")
annotation_matrix <- variant_table[, c("variant_id", "chromosome"), drop = FALSE]
gencode_overlaps <- annotate_variant_overlaps(
  variant_table,
  gencode$intervals
)
annotation_matrix <- append_annotation_columns(
  annotation_matrix,
  gencode_overlaps
)
annotation_matrix[["GENCODE intergenic"]] <-
  !annotation_matrix[["GENCODE gene body"]]
rm(gencode_overlaps, gencode)
invisible(gc())

ccre <- read_encode_ccre_annotations(annotation_paths[2L])
ccre_overlaps <- annotate_variant_overlaps(variant_table, ccre)
annotation_matrix <- append_annotation_columns(annotation_matrix, ccre_overlaps)
rm(ccre, ccre_overlaps)
invisible(gc())

roadmap_sources <- data.frame(
  path = annotation_paths[3:5],
  epigenome_id = c("E020", "E013", "E095"),
  epigenome_label = c(
    "iPS-20b",
    "hESC-derived CD56+ mesoderm",
    "left ventricle"
  ),
  stringsAsFactors = FALSE
)
for (index in seq_len(nrow(roadmap_sources))) {
  roadmap <- read_roadmap_chromhmm_annotations(
    roadmap_sources$path[index],
    roadmap_sources$epigenome_id[index],
    roadmap_sources$epigenome_label[index]
  )
  roadmap_overlaps <- annotate_variant_overlaps(variant_table, roadmap)
  annotation_matrix <- append_annotation_columns(
    annotation_matrix,
    roadmap_overlaps
  )
  rm(roadmap, roadmap_overlaps)
  invisible(gc())
}

message_step(6, 8, "Classifying exact and approximate LD-aware method specificity.")
coordinate_columns <- c("variant_id", "chromosome", "position")
fash_coordinate_rows <- match(fash_variants, variant_table$variant_id)
strober_coordinate_rows <- match(strober_union_variants,
                                 variant_table$variant_id)
ld_dosage_ids <- union(fash_variants, strober_union_variants)
ld_dosage <- vcf$dosage[match(ld_dosage_ids, rownames(vcf$dosage)), , drop = FALSE]
rownames(ld_dosage) <- ld_dosage_ids
ld_classification <- classify_cross_method_ld(
  fash_variants = variant_table[
    fash_coordinate_rows,
    coordinate_columns,
    drop = FALSE
  ],
  strober_variants = variant_table[
    strober_coordinate_rows,
    coordinate_columns,
    drop = FALSE
  ],
  dosage_matrix = ld_dosage,
  window_bp = ld_window_bp,
  r2_threshold = ld_r2_threshold
)
fash_exact_specific <- ld_classification$fash$variant_id[
  !ld_classification$fash$exact_shared
]
fash_ld_specific <- ld_classification$fash$variant_id[
  !ld_classification$fash$ld_shared
]
strober_exact_specific <- ld_classification$strober$variant_id[
  !ld_classification$strober$exact_shared
]
strober_ld_specific <- ld_classification$strober$variant_id[
  !ld_classification$strober$ld_shared
]
rm(ld_dosage)
invisible(gc())

selected_sets <- c(
  list(
    primary_fash = fash_variants,
    primary_linear = linear_variants,
    primary_nonlinear = nonlinear_variants,
    specific_fash_exact = fash_exact_specific,
    specific_fash_ld = fash_ld_specific,
    specific_strober_exact = strober_exact_specific,
    specific_strober_ld = strober_ld_specific
  ),
  equal_k_sets,
  lead_sets
)
selected_sets <- lapply(selected_sets, unique)
if (any(lengths(selected_sets) == 0L)) {
  stop("At least one planned discovery set is empty after eligibility filtering.")
}

message_step(7, 8, "Running 100 repeated matched-control analyses and chromosome jackknives.")
matched_analysis <- run_streaming_matched_enrichment(
  variant_table = variant_table,
  annotation_matrix = annotation_matrix,
  selected_sets = selected_sets,
  seeds = matching_seeds,
  controls_per_variant = controls_per_variant,
  chromosomes = as.character(1:22),
  minimum_overlap_count = 10L,
  verbose = TRUE
)

set_metadata <- data.frame(
  discovery_set = names(selected_sets),
  panel = c(
    rep("Primary reported-threshold sets", 3L),
    rep("Method-specific sensitivity", 4L),
    rep("Equal top-K sensitivity", 3L),
    rep("One lead variant per gene", 3L)
  ),
  method = c(
    "FASH", "Strober linear", "Strober nonlinear",
    "FASH exact-specific", "FASH LD-specific",
    "Strober exact-specific", "Strober LD-specific",
    "FASH", "Strober linear", "Strober nonlinear",
    "FASH", "Strober linear", "Strober nonlinear"
  ),
  n_variants = lengths(selected_sets),
  stringsAsFactors = FALSE
)
matched_analysis$results <- merge(
  matched_analysis$results,
  set_metadata,
  by = "discovery_set",
  all.x = TRUE,
  sort = FALSE
)
matched_analysis$results$q_value_within_panel <- ave(
  matched_analysis$results$p_value,
  matched_analysis$results$panel,
  FUN = function(value) stats::p.adjust(value, method = "BH")
)

annotation_coverage <- data.frame(
  annotation = setdiff(names(annotation_matrix),
                       c("variant_id", "chromosome")),
  overlap_count = colSums(annotation_matrix[, setdiff(
    names(annotation_matrix), c("variant_id", "chromosome")
  ), drop = FALSE]),
  universe_size = nrow(annotation_matrix),
  stringsAsFactors = FALSE
)
annotation_coverage$overlap_proportion <-
  annotation_coverage$overlap_count / annotation_coverage$universe_size

discovery_summary <- data.frame(
  quantity = c(
    "Tested pairs", "Tested unique variants", "Eligible matched universe",
    "FASH pairs", "FASH unique variants", "FASH genes",
    "Strober linear pairs", "Strober linear unique variants",
    "Strober nonlinear pairs", "Strober nonlinear unique variants",
    "Strober union unique variants", "Exact FASH-Strober overlap",
    "Exact FASH-specific variants", "LD-aware FASH-specific variants",
    "Exact Strober-specific variants", "LD-aware Strober-specific variants",
    "Equal top-K"
  ),
  value = c(
    expected$tested_pairs,
    expected$tested_variants,
    nrow(variant_table),
    expected$fash_pairs,
    length(fash_variants),
    expected$fash_genes,
    expected$linear_pairs,
    length(linear_variants),
    expected$nonlinear_pairs,
    length(nonlinear_variants),
    length(strober_union_variants),
    sum(ld_classification$fash$exact_shared),
    length(fash_exact_specific),
    length(fash_ld_specific),
    length(strober_exact_specific),
    length(strober_ld_specific),
    equal_k
  ),
  stringsAsFactors = FALSE
)

configuration <- list(
  analysis_id = "revision_internal_variant_annotation_enrichment",
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  assembly = "GRCh37/hg19",
  unit = "unique variant rsID",
  fdr_threshold = alpha,
  controls_per_variant = controls_per_variant,
  matching_seed_count = length(matching_seeds),
  matching_seeds = matching_seeds,
  matching_covariates = c(
    "chromosome",
    "cohort minor-allele frequency",
    "minimum tested target-gene TSS distance",
    "tested-variant density within +/-1 Mb",
    "number of tested genes per variant"
  ),
  jackknife_blocks = as.character(1:22),
  ld_window_bp = ld_window_bp,
  ld_r2_threshold = ld_r2_threshold,
  donor_count = ncol(vcf$dosage),
  claims = paste(
    "Exploratory functional plausibility only; this analysis does not",
    "validate causality or discovery FDR."
  )
)

input_provenance <- do.call(rbind, lapply(c(input_paths, annotation_paths),
                                         file_provenance))
cache <- list(
  configuration = configuration,
  expected_invariants = expected,
  input_provenance = input_provenance,
  annotation_download_manifest = utils::read.csv(
    annotation_paths[6L],
    stringsAsFactors = FALSE,
    check.names = FALSE
  ),
  discovery_summary = discovery_summary,
  set_metadata = set_metadata,
  annotation_coverage = annotation_coverage,
  enrichment_results = matched_analysis$results,
  matching_balance = matched_analysis$matching_balance,
  matching_seed_results = matched_analysis$seed_results,
  matching_relaxation_audit = matched_analysis$relaxation_audit,
  ld_classification = ld_classification,
  session_info = utils::sessionInfo()
)

message_step(8, 8, "Writing retained cache and reproducibility artifacts.")
saveRDS(
  cache,
  file.path(output_directory, "analysis_cache.rds"),
  compress = "gzip"
)
saveRDS(
  selected_sets,
  file.path(output_directory, "discovery_sets.rds"),
  compress = "gzip"
)
saveRDS(
  variant_table,
  file.path(output_directory, "variant_matching_covariates.rds"),
  compress = "gzip"
)
saveRDS(
  annotation_matrix,
  file.path(output_directory, "variant_annotation_matrix.rds"),
  compress = "gzip"
)
utils::write.csv(
  matched_analysis$results,
  file.path(output_directory, "enrichment_results.csv"),
  row.names = FALSE,
  quote = TRUE
)
utils::write.csv(
  discovery_summary,
  file.path(output_directory, "discovery_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)
utils::write.csv(
  input_provenance,
  file.path(output_directory, "input_provenance.csv"),
  row.names = FALSE,
  quote = TRUE
)

message(
  "Completed internal variant annotation enrichment cache: ",
  file.path(output_directory, "analysis_cache.rds")
)
