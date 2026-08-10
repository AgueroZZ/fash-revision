#!/usr/bin/env Rscript

# Build retained caches for the FASH versus Strober enhancer page.

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
  information <- file.info(path)
  data.frame(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    byte_size = unname(information$size),
    md5 = unname(tools::md5sum(path)),
    modified_at = format(information$mtime, tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}

message_step <- function(index, total, text) {
  message("[", index, "/", total, "] ", text)
}

workflowr_root <- find_workflowr_root()
analysis_directory <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "fash_strober_enhancer_comparison"
)
helper_path <- file.path(
  analysis_directory,
  "fash_strober_enhancer_helpers.R"
)
runner_path <- file.path(
  analysis_directory,
  "run_fash_strober_enhancer_comparison.R"
)
matching_helper_path <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "variant_annotation_enrichment",
  "variant_annotation_enrichment_helpers.R"
)
source(helper_path)
source(matching_helper_path)

current_fit_path <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
)
linear_path <- file.path(
  workflowr_root,
  "data",
  "dynamic_eQTL_real",
  "strober_linear",
  "linear_dynamic_eqtls_5_pc.txt"
)
quadratic_path <- file.path(
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
custom_annotation_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "variant_annotation_enrichment",
  "variant_annotation_matrix.rds"
)
baseline_annotation_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "baseline_ld_variant_enrichment",
  "baseline_ld_binary_annotation_matrix.rds"
)
output_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fash_strober_enhancer_comparison"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

input_paths <- c(
  current_fit_path,
  linear_path,
  quadratic_path,
  matching_covariate_path,
  custom_annotation_path,
  baseline_annotation_path,
  helper_path,
  runner_path,
  matching_helper_path
)
if (any(!file.exists(input_paths))) {
  stop("At least one retained analysis input is missing.")
}

alpha <- 0.05
matching_seeds <- seq.int(20260807L, length.out = 100L)
controls_per_variant <- 5L
analysis_started <- proc.time()[["elapsed"]]

message_step(1, 9, "Loading the current FASH result.")
current_environment <- new.env(parent = emptyenv())
loaded <- load(current_fit_path, envir = current_environment)
if (!identical(loaded, "fash_fit1_update")) {
  stop("The current FASH fit has an unexpected object name.")
}
current_fit <- current_environment$fash_fit1_update
pair_keys <- names(current_fit$fash_data$data_list)
current_lfdr <- as.numeric(current_fit$lfdr)
rm(current_fit, current_environment)
invisible(gc())

current <- derive_fash_discovery_sets(pair_keys, current_lfdr, alpha)
rm(current_lfdr)
invisible(gc())

message_step(2, 9, "Loading Strober linear and quadratic/nonlinear results.")
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("The data.table package is required.")
}
linear_result <- data.table::fread(
  linear_path,
  data.table = FALSE,
  showProgress = FALSE
)
quadratic_result <- data.table::fread(
  quadratic_path,
  data.table = FALSE,
  showProgress = FALSE
)
linear <- derive_strober_discovery_sets(linear_result, pair_keys, alpha)
quadratic <- derive_strober_discovery_sets(quadratic_result, pair_keys, alpha)
rm(linear_result, quadratic_result)
invisible(gc())

expected_counts <- list(
  current = c(9205L, 9139L, 1177L, 1170L),
  linear = c(5404L, 5387L, 550L, 548L),
  quadratic = c(6824L, 6797L, 693L, 690L)
)
observed_counts <- lapply(
  list(current = current, linear = linear, quadratic = quadratic),
  function(set) {
    c(
      nrow(set$all_pairs),
      length(set$all_variants),
      nrow(set$lead_pairs),
      length(set$lead_variants)
    )
  }
)
if (!all(vapply(
  names(expected_counts),
  function(name) identical(observed_counts[[name]], expected_counts[[name]]),
  logical(1)
))) {
  stop("At least one retained method discovery count changed unexpectedly.")
}

message_step(3, 9, "Constructing variant-level FASH-only sets.")
strober_union_variants <- union(linear$all_variants, quadratic$all_variants)

# Two FASH-only sets, both defined at the variant (rsID) level: a variant is
# FASH-only if it appears in no Strober linear or quadratic discovery, whatever
# gene it was paired with.
#
#   all  - every FASH discovery at a variant Strober never found.
#   lead - FASH's per-gene best call, restricted to the ones Strober never found.
#
# Lead selection happens BEFORE the exclusion. A gene whose best variant is also
# a Strober discovery drops out entirely rather than falling back to a weaker
# exclusive variant, so every variant here is genuinely FASH's top call for its
# gene. (The reverse ordering was tried and gives near-identical enrichment while
# substituting 78 weaker variants; see the 2026-08-09 log.)
current_only_all_pairs <- current$all_pairs[
  !current$all_pairs$variant_id %in% strober_union_variants,
  ,
  drop = FALSE
]
current_only_lead_pairs <- current$lead_pairs[
  !current$lead_pairs$variant_id %in% strober_union_variants,
  ,
  drop = FALSE
]
current_only <- list(
  all_pairs = current_only_all_pairs,
  lead_pairs = current_only_lead_pairs,
  all_variants = unique(current_only_all_pairs$variant_id),
  lead_variants = unique(current_only_lead_pairs$variant_id)
)
if (any(current_only$all_variants %in% strober_union_variants) ||
    any(current_only$lead_variants %in% strober_union_variants) ||
    !all(current_only$lead_variants %in% current_only$all_variants)) {
  stop("The FASH-only sets are not variant-level exclusive and nested.")
}
if (!identical(
  c(nrow(current_only$all_pairs), length(current_only$all_variants),
    nrow(current_only$lead_pairs), length(current_only$lead_variants)),
  c(8033L, 7972L, 1036L, 1030L)
)) {
  stop("The FASH-only discovery counts changed unexpectedly.")
}

selected_sets <- list(
  current_all = current$all_variants,
  current_lead = current$lead_variants,
  linear_all = linear$all_variants,
  linear_lead = linear$lead_variants,
  quadratic_all = quadratic$all_variants,
  quadratic_lead = quadratic$lead_variants,
  current_only_all = current_only$all_variants,
  current_only_lead = current_only$lead_variants
)
set_metadata <- data.frame(
  discovery_set = names(selected_sets),
  method = c(
    "Current FASH", "Current FASH", "Strober linear", "Strober linear",
    "Strober quadratic/nonlinear", "Strober quadratic/nonlinear",
    "Current FASH", "Current FASH"
  ),
  section = c(rep(1L, 6L), rep(2L, 2L)),
  exclusivity = c(rep("None", 6L), rep("Variant-level", 2L)),
  view = c(rep("Ordinary", 6L), rep("FASH-only variants", 2L)),
  selection_strategy = rep(
    c("All variants", "One lead variant per gene"),
    4L
  ),
  original_variant_count = lengths(selected_sets),
  stringsAsFactors = FALSE
)

pair_tables <- list(
  current_all = current$all_pairs,
  current_lead = current$lead_pairs,
  linear_all = linear$all_pairs,
  linear_lead = linear$lead_pairs,
  quadratic_all = quadratic$all_pairs,
  quadratic_lead = quadratic$lead_pairs,
  current_only_all = current_only$all_pairs,
  current_only_lead = current_only$lead_pairs
)
discovery_set_summary <- do.call(rbind, lapply(names(selected_sets), function(
    set_name) {
  metadata <- set_metadata[
    set_metadata$discovery_set == set_name,
    ,
    drop = FALSE
  ]
  summarize_discovery_set(
    method = metadata$method,
    view = metadata$view,
    selection_strategy = metadata$selection_strategy,
    pairs = pair_tables[[set_name]],
    variants = selected_sets[[set_name]]
  )
}))
discovery_set_summary$discovery_set <- names(selected_sets)
discovery_set_summary <- discovery_set_summary[, c(
  "discovery_set", "method", "view", "selection_strategy", "pair_count",
  "unique_variant_count", "unique_gene_count"
)]
discovery_overlap <- build_discovery_overlap(selected_sets, set_metadata)
strict_exclusivity_audit <- data.frame(
  method = "Current FASH",
  original_variant_count = length(current$all_variants),
  excluded_variant_count = length(
    intersect(current$all_variants, strober_union_variants)
  ),
  fash_only_variant_count = length(current_only$all_variants),
  original_lead_variant_count = length(current$lead_variants),
  excluded_lead_variant_count = length(
    intersect(current$lead_variants, strober_union_variants)
  ),
  fash_only_lead_variant_count = length(current_only$lead_variants),
  fash_only_lead_gene_count = length(unique(current_only$lead_pairs$gene_id)),
  overlap_with_strober_union =
    length(intersect(current_only$all_variants, strober_union_variants)),
  stringsAsFactors = FALSE
)

message_step(4, 9, "Loading and validating retained annotation matrices.")
variant_table <- readRDS(matching_covariate_path)
custom_matrix <- readRDS(custom_annotation_path)
baseline_matrix <- readRDS(baseline_annotation_path)
if (nrow(variant_table) != 745867L || anyDuplicated(variant_table$variant_id) ||
    !identical(variant_table$variant_id, custom_matrix$variant_id) ||
    nrow(baseline_matrix) != 557405L || anyDuplicated(baseline_matrix$variant_id)) {
  stop("Retained annotation or matching inputs are not aligned.")
}
for (matrix in list(custom_matrix, baseline_matrix)) {
  annotations <- setdiff(names(matrix), c("variant_id", "chromosome"))
  if (!all(vapply(matrix[annotations], is.logical, logical(1)))) {
    stop("Retained annotation columns must be logical.")
  }
}
custom_enhancers <- c(
  "ENCODE cCRE enhancer-like",
  "Roadmap E020 iPS-20b: Enhancer",
  "Roadmap E013 hESC-derived CD56+ mesoderm: Enhancer",
  "Roadmap E095 left ventricle: Enhancer"
)
baseline_enhancers <- c(
  "Enhancer_Andersson", "Enhancer_Hoffman", "WeakEnhancer_Hoffman",
  "SuperEnhancer_Hnisz", "Human_Enhancer_Villar"
)
if (!all(custom_enhancers %in% names(custom_matrix)) ||
    !all(baseline_enhancers %in% names(baseline_matrix))) {
  stop("A pre-specified enhancer annotation is missing.")
}

build_coverage <- function(annotation_system, covered_ids) {
  do.call(rbind, lapply(names(selected_sets), function(set_name) {
    ids <- unique(selected_sets[[set_name]])
    data.frame(
      annotation_system = annotation_system,
      discovery_set = set_name,
      original_variant_count = length(ids),
      covered_variant_count = sum(ids %in% covered_ids),
      coverage_proportion = mean(ids %in% covered_ids),
      stringsAsFactors = FALSE
    )
  }))
}
coverage_by_set <- rbind(
  build_coverage("Custom regulatory", custom_matrix$variant_id),
  build_coverage("baselineLD v2.2", baseline_matrix$variant_id)
)

message_step(5, 9, "Running custom-regulatory matched enrichment.")
custom_sets <- lapply(selected_sets, function(ids) {
  intersect(ids, custom_matrix$variant_id)
})
custom_started <- proc.time()[["elapsed"]]
custom_analysis <- run_streaming_matched_enrichment(
  variant_table = variant_table,
  annotation_matrix = custom_matrix,
  selected_sets = custom_sets,
  seeds = matching_seeds,
  controls_per_variant = controls_per_variant,
  chromosomes = as.character(1:22),
  minimum_overlap_count = 10L,
  verbose = TRUE
)
custom_elapsed <- proc.time()[["elapsed"]] - custom_started

message_step(6, 9, "Running baselineLD matched enrichment.")
baseline_variant_table <- variant_table[
  match(baseline_matrix$variant_id, variant_table$variant_id),
  ,
  drop = FALSE
]
if (anyNA(baseline_variant_table$variant_id) ||
    !identical(baseline_variant_table$variant_id, baseline_matrix$variant_id)) {
  stop("baselineLD annotations and matching covariates are not aligned.")
}
baseline_sets <- lapply(selected_sets, function(ids) {
  intersect(ids, baseline_matrix$variant_id)
})
baseline_started <- proc.time()[["elapsed"]]
baseline_analysis <- run_streaming_matched_enrichment(
  variant_table = baseline_variant_table,
  annotation_matrix = baseline_matrix,
  selected_sets = baseline_sets,
  seeds = matching_seeds,
  controls_per_variant = controls_per_variant,
  chromosomes = as.character(1:22),
  minimum_overlap_count = 10L,
  verbose = TRUE
)
baseline_elapsed <- proc.time()[["elapsed"]] - baseline_started

message_step(7, 9, "Adding metadata, enhancer summaries, and diagnostics.")
add_metadata <- function(analysis, annotation_system) {
  results <- merge(
    analysis$results,
    set_metadata,
    by = "discovery_set",
    all.x = TRUE,
    sort = FALSE
  )
  results$annotation_system <- annotation_system
  results$q_value_within_set <- ave(
    results$p_value,
    results$discovery_set,
    FUN = function(values) stats::p.adjust(values, method = "BH")
  )
  results
}
custom_results <- add_metadata(custom_analysis, "Custom regulatory")
baseline_results <- add_metadata(baseline_analysis, "baselineLD v2.2")
enrichment_results <- rbind(custom_results, baseline_results)

# Unmatched over-representation against the full annotated tested-variant
# universe. This is reported alongside the matched estimate so that the effect
# of covariate matching is visible rather than implicit. No p-value is attached:
# a Fisher test here would be badly anti-conservative under LD.
unmatched_custom <- compute_unmatched_enrichment(
  custom_matrix,
  selected_sets,
  custom_enhancers
)
unmatched_custom$annotation_system <- "Custom regulatory"
unmatched_baseline <- compute_unmatched_enrichment(
  baseline_matrix,
  selected_sets,
  baseline_enhancers
)
unmatched_baseline$annotation_system <- "baselineLD v2.2"
unmatched_enrichment <- rbind(unmatched_custom, unmatched_baseline)

enhancer_results <- rbind(
  custom_results[custom_results$annotation %in% custom_enhancers, ],
  baseline_results[baseline_results$annotation %in% baseline_enhancers, ]
)
enhancer_results <- merge(
  enhancer_results,
  unmatched_enrichment,
  by = c("annotation_system", "discovery_set", "annotation"),
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(enhancer_results$unmatched_enrichment) ||
    nrow(enhancer_results) != length(selected_sets) *
      (length(custom_enhancers) + length(baseline_enhancers))) {
  stop("The unmatched enrichment merge did not align with enhancer results.")
}
enhancer_summary <- summarize_enhancer_panel(enhancer_results)

custom_balance <- custom_analysis$matching_balance
custom_balance$annotation_system <- "Custom regulatory"
baseline_balance <- baseline_analysis$matching_balance
baseline_balance$annotation_system <- "baselineLD v2.2"
matching_balance <- rbind(custom_balance, baseline_balance)
custom_relaxation <- custom_analysis$relaxation_audit
custom_relaxation$annotation_system <- "Custom regulatory"
baseline_relaxation <- baseline_analysis$relaxation_audit
baseline_relaxation$annotation_system <- "baselineLD v2.2"
matching_relaxation <- rbind(custom_relaxation, baseline_relaxation)
maximum_absolute_smd <- max(
  abs(matching_balance$standardized_mean_difference),
  na.rm = TRUE
)
if (!is.finite(maximum_absolute_smd) || maximum_absolute_smd >= 0.10) {
  stop("Matched controls failed the pre-specified balance diagnostic.")
}

message_step(8, 9, "Writing tidy variant-list and estimate exports.")
# Local only: output/ is gitignored. The repository carries code, not data.
# The full annotated universe is not re-exported here because it already exists
# as variant_annotation_matrix.rds / variant_matching_covariates.rds; these two
# small files are the tidy form of what this experiment adds on top.
score_types <- c(
  "Current FASH" = "cumulative local false discovery rate (lower is stronger)",
  "Strober linear" = "p-value (lower is stronger)",
  "Strober quadratic/nonlinear" = "p-value (lower is stronger)"
)
discovery_export <- do.call(rbind, lapply(names(pair_tables), function(name) {
  pairs <- pair_tables[[name]]
  metadata <- set_metadata[set_metadata$discovery_set == name, , drop = FALSE]
  data.frame(
    discovery_set = name,
    method = metadata$method,
    exclusivity = metadata$exclusivity,
    selection_strategy = metadata$selection_strategy,
    gene_id = pairs$gene_id,
    variant_id = pairs$variant_id,
    ranking_score = pairs$score,
    ranking_score_type = unname(score_types[[metadata$method]]),
    stringsAsFactors = FALSE
  )
}))

published_estimates <- enhancer_results[, c(
  "annotation_system", "discovery_set", "method", "exclusivity",
  "selection_strategy", "annotation", "selected_overlap", "selected_total",
  "control_overlap", "control_total", "enrichment", "log2_enrichment",
  "jackknife_se", "ci_lower", "ci_upper", "p_value", "unmatched_enrichment"
)]

total_elapsed <- proc.time()[["elapsed"]] - analysis_started
runtime_summary <- data.frame(
  stage = c(
    "custom_matched_enrichment",
    "baseline_ld_matched_enrichment",
    "analysis_through_cache_assembly"
  ),
  elapsed_seconds = c(custom_elapsed, baseline_elapsed, total_elapsed),
  stringsAsFactors = FALSE
)
input_provenance <- do.call(rbind, lapply(input_paths, file_provenance))
roadmap_base <- paste0(
  "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/",
  "ChmmModels/coreMarks/jointModel/final/"
)
configuration <- list(
  analysis_id = "revision_internal_fash_strober_enhancer_comparison",
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  fdr_threshold = alpha,
  controls_per_variant = controls_per_variant,
  matching_seed_count = length(matching_seeds),
  matching_seeds = matching_seeds,
  jackknife_blocks = as.character(1:22),
  strict_only_definition = paste(
    "A focal-method rsID is removed if it appears in any Strober linear or",
    "quadratic/nonlinear discovery, regardless of the associated gene."
  ),
  strict_lead_definition = paste(
    "Exclude all Strober-discovered rsIDs first, then select the lowest-lfdr",
    "remaining focal-method pair per gene."
  ),
  multiplicity_statement = paste(
    "Exploratory analysis. Fold enrichment, 95% delete-one-autosome intervals,",
    "and nominal p-values are reported without multiplicity correction across",
    "the 8 discovery sets and 2 annotation systems."
  ),
  annotation_provenance = data.frame(
    annotation_system = c(
      rep("Custom regulatory", 4L), rep("baselineLD v2.2", 5L)
    ),
    annotation = c(custom_enhancers, baseline_enhancers),
    source = c(
      "ENCODE candidate cis-regulatory elements (cell-type agnostic)",
      "NIH Roadmap 15-state ChromHMM coreMarks, iPS-20b",
      "NIH Roadmap 15-state ChromHMM coreMarks, hESC-derived CD56+ mesoderm",
      "NIH Roadmap 15-state ChromHMM coreMarks, left ventricle",
      rep("baselineLD v2.2 binary annotation columns", 5L)
    ),
    accession_or_record = c(
      "ENCFF788SJC", "E020", "E013", "E095",
      rep("Zenodo 10515792", 5L)
    ),
    genome_build = c(rep("hg19", 4L), rep("GRCh37/hg19", 5L)),
    repository = c(
      "encodeproject.org",
      rep("egg2.wustl.edu/roadmap", 3L),
      rep("zenodo.org (retrieved from Song Lab mirror)", 5L)
    ),
    download_file = c(
      "ENCFF788SJC.bed.gz",
      "E020_15_coreMarks_mnemonics.bed.gz",
      "E013_15_coreMarks_mnemonics.bed.gz",
      "E095_15_coreMarks_mnemonics.bed.gz",
      rep("baselineLD.{1..22}.annot.gz", 5L)
    ),
    download_url = c(
      paste0(
        "https://www.encodeproject.org/files/ENCFF788SJC/@@download/",
        "ENCFF788SJC.bed.gz"
      ),
      paste0(roadmap_base, "E020_15_coreMarks_mnemonics.bed.gz"),
      paste0(roadmap_base, "E013_15_coreMarks_mnemonics.bed.gz"),
      paste0(roadmap_base, "E095_15_coreMarks_mnemonics.bed.gz"),
      rep("https://zenodo.org/records/10515792", 5L)
    ),
    retrieval_url = c(
      rep(NA_character_, 4L),
      rep(
        paste0(
          "https://huggingface.co/datasets/songlab/ldsc/tree/main/",
          "baselineLD_v2.2"
        ),
        5L
      )
    ),
    stringsAsFactors = FALSE
  ),
  custom_enhancers = custom_enhancers,
  baseline_enhancers = baseline_enhancers,
  annotation_systems = c("Custom regulatory", "baselineLD v2.2")
)
cache <- list(
  configuration = configuration,
  input_provenance = input_provenance,
  set_metadata = set_metadata,
  discovery_set_summary = discovery_set_summary,
  discovery_overlap = discovery_overlap,
  strict_exclusivity_audit = strict_exclusivity_audit,
  discovery_export = discovery_export,
  coverage_by_set = coverage_by_set,
  enrichment_results = enrichment_results,
  unmatched_enrichment = unmatched_enrichment,
  enhancer_results = enhancer_results,
  enhancer_summary = enhancer_summary,
  matching_balance = matching_balance,
  matching_relaxation = matching_relaxation,
  maximum_absolute_smd = maximum_absolute_smd,
  runtime_summary = runtime_summary,
  session_info = utils::sessionInfo()
)

message_step(9, 9, "Writing retained cache and reproducibility artifacts.")
saveRDS(cache, file.path(output_directory, "analysis_cache.rds"),
        compress = "gzip")
saveRDS(selected_sets, file.path(output_directory, "discovery_sets.rds"),
        compress = "gzip")
saveRDS(pair_tables, file.path(output_directory, "discovery_pair_tables.rds"),
        compress = "gzip")
write_output <- function(value, filename) {
  utils::write.csv(
    value,
    file.path(output_directory, filename),
    row.names = FALSE,
    quote = TRUE
  )
}
write_output(discovery_set_summary, "discovery_set_summary.csv")
write_output(discovery_overlap, "discovery_overlap.csv")
write_output(strict_exclusivity_audit, "strict_exclusivity_audit.csv")
write_output(coverage_by_set, "coverage_by_set.csv")
write_output(enrichment_results, "enrichment_results.csv")
write_output(enhancer_results, "enhancer_results.csv")
write_output(enhancer_summary, "enhancer_summary.csv")
write_output(matching_balance, "matching_balance.csv")
write_output(runtime_summary, "runtime_summary.csv")
write_output(input_provenance, "input_provenance.csv")
write_output(configuration$annotation_provenance, "annotation_provenance.csv")
write_output(discovery_export, "discovery_sets.csv")
write_output(published_estimates, "published_estimates.csv")
message(
  "Completed FASH/FASH-CL versus Strober enhancer comparison in ",
  round(total_elapsed, 1),
  " seconds: ",
  file.path(output_directory, "analysis_cache.rds")
)
