#!/usr/bin/env Rscript

# Build retained caches for the internal baselineLD v2.2 enrichment page.

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

file_provenance <- function(path) {
  data.frame(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    byte_size = unname(file.info(path)$size),
    md5 = unname(tools::md5sum(path)),
    modified_at = format(file.info(path)$mtime, tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}

message_step <- function(index, total, text) {
  message("[", index, "/", total, "] ", text)
}

workflowr_root <- find_workflowr_root()
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.")
}

analysis_directory <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "baseline_ld_variant_enrichment"
)
source(file.path(
  analysis_directory,
  "baseline_ld_variant_enrichment_helpers.R"
))
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "variant_annotation_enrichment",
  "variant_annotation_enrichment_helpers.R"
))

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
matching_covariate_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "variant_annotation_enrichment",
  "variant_matching_covariates.rds"
)
annotation_directory <- file.path(
  workflowr_root,
  "data",
  "revision_simulations",
  "internal",
  "baseline_ld_variant_enrichment"
)
download_manifest_path <- file.path(
  annotation_directory,
  "download_manifest.csv"
)
output_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "baseline_ld_variant_enrichment"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

download_manifest <- utils::read.csv(
  download_manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (!identical(sort(as.integer(download_manifest$chromosome)), 1:22) ||
    any(!file.exists(download_manifest$local_path)) ||
    sum(download_manifest$byte_size) != 416186492) {
  stop("The validated baselineLD download manifest is incomplete.")
}
annotation_paths <- setNames(
  download_manifest$local_path,
  as.character(download_manifest$chromosome)
)

alpha <- 0.05
matching_seeds <- seq.int(20260807L, length.out = 100L)
controls_per_variant <- 5L
analysis_started <- proc.time()[["elapsed"]]

message_step(1, 7, "Deriving the three requested FASH discovery sets.")
fit_environment <- new.env(parent = emptyenv())
loaded_names <- load(fit_path, envir = fit_environment)
if (!identical(loaded_names, "fash_fit1_update")) {
  stop("The FASH fit file did not contain only fash_fit1_update.")
}
fit <- fit_environment$fash_fit1_update
pair_keys <- names(fit$fash_data$data_list)
linear_results <- data.table::fread(
  linear_path,
  data.table = FALSE,
  showProgress = FALSE
)
nonlinear_results <- data.table::fread(
  nonlinear_path,
  data.table = FALSE,
  showProgress = FALSE
)
requested <- derive_requested_fash_sets(
  pair_keys = pair_keys,
  lfdr = fit$lfdr,
  linear_results = linear_results,
  nonlinear_results = nonlinear_results,
  alpha = alpha
)
expected_pair_counts <- c(9205L, 1177L, 8062L)
expected_variant_counts <- c(9139L, 1170L, 8001L)
if (!identical(requested$pair_summary$pair_count, expected_pair_counts) ||
    !identical(
      requested$pair_summary$unique_variant_count,
      expected_variant_counts
    )) {
  stop("The requested FASH discovery-set counts changed unexpectedly.")
}
rm(fit, fit_environment, linear_results, nonlinear_results)
invisible(gc())

message_step(2, 7, "Loading the retained tested-variant matching covariates.")
variant_table <- readRDS(matching_covariate_path)
if (nrow(variant_table) != 745867L ||
    anyDuplicated(variant_table$variant_id) ||
    !all(c(
      "variant_id", "chromosome", "position", "minor_allele_frequency",
      "minimum_target_tss_distance", "local_tested_variant_count_1mb",
      "n_tested_genes", "full_stratum", "drop_density_stratum",
      "drop_distance_stratum"
    ) %in% names(variant_table))) {
  stop("The retained variant matching covariates are invalid.")
}

message_step(3, 7, "Classifying baselineLD columns and streaming annotations.")
annotation_started <- proc.time()[["elapsed"]]
headers <- lapply(annotation_paths, read_baseline_ld_header)
if (!all(vapply(headers, identical, logical(1), headers[[1L]])) ||
    length(headers[[1L]]) != 101L) {
  stop("baselineLD chromosome files do not share the expected 101-column schema.")
}
chromosome_one <- data.table::fread(
  annotation_paths[["1"]],
  showProgress = FALSE,
  data.table = FALSE
)
column_classification <- classify_baseline_ld_columns(chromosome_one)
binary_annotations <- column_classification$annotation[
  column_classification$annotation_type == "binary"
]
continuous_annotations <- column_classification$annotation[
  column_classification$annotation_type == "continuous"
]
if (length(binary_annotations) != 84L ||
    length(continuous_annotations) != 13L ||
    !"base" %in% binary_annotations) {
  stop("The baselineLD v2.2 annotation-type classification changed.")
}
rm(chromosome_one)
invisible(gc())

baseline_matrix <- build_baseline_ld_matrix(
  paths = annotation_paths,
  variant_table = variant_table,
  binary_annotation_columns = binary_annotations
)
annotation_elapsed <- proc.time()[["elapsed"]] - annotation_started
covered_ids <- baseline_matrix$variant_id
coverage_by_set <- summarize_set_coverage(
  selected_sets = requested$selected_sets,
  covered_ids = covered_ids,
  universe_ids = variant_table$variant_id
)
background_coverage <- coverage_by_set$coverage_proportion[
  coverage_by_set$discovery_set == "tested_variant_universe"
]
set_coverage <- coverage_by_set$coverage_proportion[
  coverage_by_set$discovery_set != "tested_variant_universe"
]
if (any(set_coverage < 0.60) ||
    any(abs(set_coverage - background_coverage) > 0.10)) {
  stop("baselineLD coverage is too low or too imbalanced for analysis.")
}

covered_variant_table <- variant_table[
  match(covered_ids, variant_table$variant_id),
  ,
  drop = FALSE
]
if (anyNA(covered_variant_table$variant_id) ||
    !identical(covered_variant_table$variant_id, baseline_matrix$variant_id)) {
  stop("The baselineLD and matching-covariate rows are not aligned.")
}
covered_selected_sets <- lapply(requested$selected_sets, function(ids) {
  intersect(as.character(ids), covered_ids)
})
if (any(lengths(covered_selected_sets) == 0L)) {
  stop("At least one requested set is empty after baselineLD coverage filtering.")
}

coverage_by_chromosome <- do.call(rbind, lapply(
  c(list(tested_variant_universe = variant_table$variant_id),
    requested$selected_sets),
  function(ids) {
    rows <- match(unique(as.character(ids)), variant_table$variant_id)
    chromosomes <- variant_table$chromosome[rows]
    data.frame(
      chromosome = as.character(1:22),
      original_count = as.integer(table(factor(
        chromosomes,
        levels = as.character(1:22)
      ))),
      covered_count = as.integer(table(factor(
        chromosomes[unique(as.character(ids)) %in% covered_ids],
        levels = as.character(1:22)
      ))),
      stringsAsFactors = FALSE
    )
  }
))
coverage_by_chromosome$discovery_set <- rep(
  c("tested_variant_universe", names(requested$selected_sets)),
  each = 22L
)
coverage_by_chromosome$coverage_proportion <- with(
  coverage_by_chromosome,
  covered_count / original_count
)

analysis_annotations <- setdiff(binary_annotations, "base")
annotation_matrix <- baseline_matrix[, c(
  "variant_id", "chromosome", analysis_annotations
), drop = FALSE]
if (!all(vapply(
  annotation_matrix[analysis_annotations],
  is.logical,
  logical(1)
))) {
  stop("The retained baselineLD analysis annotations must be logical.")
}

core_annotations <- c(
  "Coding_UCSC", "Conserved_LindbladToh", "CTCF_Hoffman", "DGF_ENCODE",
  "DHS_Trynka", "Enhancer_Andersson", "Enhancer_Hoffman",
  "FetalDHS_Trynka", "H3K27ac_Hnisz", "H3K27ac_PGC2",
  "H3K4me1_Trynka", "H3K4me3_Trynka", "H3K9ac_Trynka",
  "Intron_UCSC", "PromoterFlanking_Hoffman", "Promoter_UCSC",
  "Repressed_Hoffman", "SuperEnhancer_Hnisz", "TFBS_ENCODE",
  "Transcr_Hoffman", "TSS_Hoffman", "UTR_3_UCSC", "UTR_5_UCSC",
  "WeakEnhancer_Hoffman", "synonymous", "non_synonymous", "BivFlnk",
  "Human_Promoter_Villar", "Human_Enhancer_Villar"
)
if (!all(core_annotations %in% analysis_annotations)) {
  stop("At least one pre-specified core annotation is absent.")
}

message_step(4, 7, "Running 100-seed matched enrichment and jackknives.")
matching_started <- proc.time()[["elapsed"]]
matched_analysis <- run_streaming_matched_enrichment(
  variant_table = covered_variant_table,
  annotation_matrix = annotation_matrix,
  selected_sets = covered_selected_sets,
  seeds = matching_seeds,
  controls_per_variant = controls_per_variant,
  chromosomes = as.character(1:22),
  minimum_overlap_count = 10L,
  verbose = TRUE
)
matching_elapsed <- proc.time()[["elapsed"]] - matching_started

message_step(5, 7, "Adding discovery-set, annotation, and coverage metadata.")
set_metadata <- data.frame(
  discovery_set = names(covered_selected_sets),
  display_label = c(
    "All FASH discoveries",
    "One FASH lead variant per gene",
    "Variants in FASH-only gene-variant pairs"
  ),
  original_variant_count = lengths(requested$selected_sets),
  covered_variant_count = lengths(covered_selected_sets),
  stringsAsFactors = FALSE
)
matched_analysis$results <- merge(
  matched_analysis$results,
  set_metadata,
  by = "discovery_set",
  all.x = TRUE,
  sort = FALSE
)
matched_analysis$results$q_value_within_set <- ave(
  matched_analysis$results$p_value,
  matched_analysis$results$discovery_set,
  FUN = function(value) stats::p.adjust(value, method = "BH")
)
matched_analysis$results$core_annotation <-
  matched_analysis$results$annotation %in% core_annotations

annotation_metadata <- data.frame(
  annotation = analysis_annotations,
  annotation_type = "binary",
  core_annotation = analysis_annotations %in% core_annotations,
  flanking_500 = grepl("[.]flanking[.]500$", analysis_annotations),
  maf_bin = grepl("^MAFbin[0-9]+$", analysis_annotations),
  overlap_count = vapply(
    annotation_matrix[analysis_annotations],
    sum,
    numeric(1)
  ),
  universe_size = nrow(annotation_matrix),
  stringsAsFactors = FALSE
)
annotation_metadata$overlap_proportion <-
  annotation_metadata$overlap_count / annotation_metadata$universe_size

coverage_by_set$coverage_difference_from_background <-
  coverage_by_set$coverage_proportion - background_coverage
requested$pair_summary$covered_unique_variant_count <- lengths(
  covered_selected_sets
)
requested$pair_summary$coverage_proportion <-
  requested$pair_summary$covered_unique_variant_count /
  requested$pair_summary$unique_variant_count

maximum_absolute_smd <- max(
  abs(matched_analysis$matching_balance$standardized_mean_difference),
  na.rm = TRUE
)
if (!is.finite(maximum_absolute_smd) || maximum_absolute_smd > 0.10) {
  stop("Matched controls failed the pre-specified balance diagnostic.")
}

message_step(6, 7, "Building provenance and local resource diagnostics.")
input_paths <- c(
  fit_path,
  linear_path,
  nonlinear_path,
  matching_covariate_path,
  download_manifest_path,
  unname(annotation_paths)
)
input_provenance <- do.call(rbind, lapply(input_paths, file_provenance))
total_elapsed <- proc.time()[["elapsed"]] - analysis_started
runtime_summary <- data.frame(
  stage = c(
    "annotation_download_validation",
    "annotation_parse_and_subset",
    "matched_enrichment",
    "analysis_through_cache_assembly"
  ),
  elapsed_seconds = c(NA_real_, annotation_elapsed, matching_elapsed,
                      total_elapsed),
  note = c(
    "Measured separately during the validated download run",
    "Includes column classification and 22 chromosome subsets",
    "Three sets, 100 seeds, five controls, and chromosome jackknives",
    "Excludes final compressed cache serialization"
  ),
  stringsAsFactors = FALSE
)
resource_summary <- data.frame(
  quantity = c(
    "Downloaded annotation bytes",
    "Downloaded annotation MiB",
    "Covered tested variants",
    "Binary annotation columns analyzed",
    "Continuous annotation columns documented but not analyzed",
    "Binary annotation matrix bytes in R",
    "Covered matching table bytes in R",
    "Maximum absolute matching SMD"
  ),
  value = c(
    sum(download_manifest$byte_size),
    sum(download_manifest$byte_size) / 1024^2,
    nrow(annotation_matrix),
    length(analysis_annotations),
    length(continuous_annotations),
    as.numeric(object.size(annotation_matrix)),
    as.numeric(object.size(covered_variant_table)),
    maximum_absolute_smd
  ),
  stringsAsFactors = FALSE
)

configuration <- list(
  analysis_id = "revision_internal_baseline_ld_variant_enrichment",
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  assembly = "GRCh37/hg19",
  canonical_record = unique(download_manifest$canonical_record),
  canonical_archive = unique(download_manifest$canonical_archive),
  annotation_version = "baselineLD v2.2",
  annotation_input = "Chromosome-specific .annot.gz columns",
  analysis_unit = "Unique rsID after discovery-set construction",
  fash_only_definition = paste(
    "Exact gene-variant pair discovered by FASH but by neither Strober",
    "linear nor Strober nonlinear; variants are deduplicated for enrichment."
  ),
  fdr_threshold = alpha,
  controls_per_variant = controls_per_variant,
  matching_seed_count = length(matching_seeds),
  matching_seeds = matching_seeds,
  jackknife_blocks = as.character(1:22),
  core_annotations = core_annotations,
  claims = paste(
    "Exploratory functional plausibility only; annotation enrichment does",
    "not validate causality or discovery FDR."
  )
)

cache <- list(
  configuration = configuration,
  input_provenance = input_provenance,
  download_manifest = download_manifest,
  column_classification = column_classification,
  pair_summary = requested$pair_summary,
  set_metadata = set_metadata,
  coverage_by_set = coverage_by_set,
  coverage_by_chromosome = coverage_by_chromosome,
  annotation_metadata = annotation_metadata,
  enrichment_results = matched_analysis$results,
  matching_balance = matched_analysis$matching_balance,
  matching_seed_results = matched_analysis$seed_results,
  matching_relaxation_audit = matched_analysis$relaxation_audit,
  fash_only_variant_multiplicity = requested$fash_only_variant_multiplicity,
  runtime_summary = runtime_summary,
  resource_summary = resource_summary,
  session_info = utils::sessionInfo()
)

message_step(7, 7, "Writing retained caches and reproducibility artifacts.")
saveRDS(
  cache,
  file.path(output_directory, "analysis_cache.rds"),
  compress = "gzip"
)
saveRDS(
  covered_selected_sets,
  file.path(output_directory, "requested_discovery_sets.rds"),
  compress = "gzip"
)
saveRDS(
  annotation_matrix,
  file.path(output_directory, "baseline_ld_binary_annotation_matrix.rds"),
  compress = "gzip"
)
utils::write.csv(
  matched_analysis$results,
  file.path(output_directory, "enrichment_results.csv"),
  row.names = FALSE,
  quote = TRUE
)
utils::write.csv(
  coverage_by_set,
  file.path(output_directory, "coverage_by_set.csv"),
  row.names = FALSE,
  quote = TRUE
)
utils::write.csv(
  matched_analysis$matching_balance,
  file.path(output_directory, "matching_balance.csv"),
  row.names = FALSE,
  quote = TRUE
)
utils::write.csv(
  runtime_summary,
  file.path(output_directory, "runtime_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)
utils::write.csv(
  resource_summary,
  file.path(output_directory, "resource_summary.csv"),
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
  "Completed internal baselineLD enrichment cache in ",
  round(total_elapsed, 1),
  " seconds: ",
  file.path(output_directory, "analysis_cache.rds")
)
