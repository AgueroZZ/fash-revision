#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

workflowr_root <- if (file.exists("_workflowr.yml")) {
  "."
} else if (file.exists("coderepo-local/_workflowr.yml")) {
  "coderepo-local"
} else {
  stop("Run this script from the workflowr root or its parent.")
}

analysis_id <- "matched_fash_linear_real_data_fashr_0_1_43"
code_directory <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "matched_fash_linear_real_data"
)
helper_path <- file.path(code_directory, "matched_fash_linear_helpers.R")
overlap_helper_path <- file.path(
  code_directory,
  "iwp_overlap_proportion_helpers.R"
)
selection_helper_path <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "fash_strober_enhancer_comparison", "fash_strober_enhancer_helpers.R"
)
builder_path <- file.path(code_directory, "build_iwp_overlap_proportion_cache.R")
output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal", analysis_id
)
matched_cache_path <- file.path(output_directory, "analysis_cache.rds")
iwp2_path <- file.path(
  workflowr_root,
  "output", "dynamic_eQTL_real", "fash_fit2_update.RData"
)
output_path <- file.path(output_directory, "iwp_overlap_proportions.rds")

source(helper_path)
source(overlap_helper_path)
source(selection_helper_path)

input_paths <- c(
  matched_analysis_cache = matched_cache_path,
  current_iwp2_bf = iwp2_path,
  matched_helpers = helper_path,
  overlap_helpers = overlap_helper_path,
  selection_helpers = selection_helper_path,
  overlap_cache_builder = builder_path
)
if (any(!file.exists(input_paths))) {
  stop("At least one directional-overlap input is missing.")
}

package_provenance <- extract_package_provenance_matched()
expected_package_sha <- "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
if (!identical(package_provenance$version, "0.1.43") ||
    !identical(package_provenance$remote_sha, expected_package_sha)) {
  stop("The installed fashr package is not the approved 0.1.43 build.")
}

matched_cache <- readRDS(matched_cache_path)
if (!identical(matched_cache$configuration$analysis_id, analysis_id) ||
    !identical(matched_cache$configuration$alpha, 0.05) ||
    any(!matched_cache$validation$passed)) {
  stop("The matched FASH-linear cache is invalid.")
}

linear_sets <- lapply(matched_cache$venn_sets, function(sets) {
  sets[["FASH-linear BF"]]
})
iwp1_sets <- lapply(matched_cache$venn_sets, function(sets) {
  sets[["FASH-IWP1 BF"]]
})
iwp1_linear <- directional_overlap_summary_matched(
  reference_sets = linear_sets,
  comparison_sets = iwp1_sets
)

iwp2 <- load_exact_object_matched(iwp2_path, "fash_fit2_update")
pair_keys <- names(iwp2$fash_data$data_list)
if (!inherits(iwp2, "fash") || length(iwp2$lfdr) != 1009173L ||
    length(pair_keys) != length(iwp2$lfdr) || anyDuplicated(pair_keys)) {
  stop("The current BF-adjusted IWP2 fit is invalid or misaligned.")
}
iwp2_indices <- select_cumulative_lfdr_calls(iwp2$lfdr, alpha = 0.05)
iwp2_pairs <- pair_keys[iwp2_indices]
iwp2_sets <- list(
  `Gene-variant pairs` = iwp2_pairs,
  Genes = unique(sub("_.*$", "", iwp2_pairs)),
  Variants = unique(sub("^[^_]+_", "", iwp2_pairs))
)
iwp2_linear <- directional_overlap_summary_matched(
  reference_sets = linear_sets,
  comparison_sets = iwp2_sets
)

stopifnot(
  identical(iwp1_linear$reference_count, c(14902L, 1663L, 14761L)),
  identical(iwp1_linear$comparison_count, c(9214L, 1176L, 9148L)),
  identical(iwp1_linear$intersection_count, c(8530L, 1112L, 8468L)),
  identical(iwp2_linear$reference_count, c(14902L, 1663L, 14761L)),
  identical(iwp2_linear$comparison_count, c(44L, 9L, 44L)),
  identical(iwp2_linear$intersection_count, c(21L, 6L, 21L))
)

overlap_cache <- list(
  analysis_id = "fash_iwp_linear_directional_overlap_fashr_0_1_43",
  created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  configuration = list(
    package_version = package_provenance$version,
    package_sha = package_provenance$remote_sha,
    adjustment = "BF-adjusted",
    discovery_rule = "cumulative-lfdr FDR 0.05",
    alpha = 0.05
  ),
  input_provenance = make_file_provenance_matched(input_paths),
  iwp1_linear = iwp1_linear,
  iwp2_linear = iwp2_linear
)
atomic_save_rds_matched(overlap_cache, output_path)

print(iwp1_linear, row.names = FALSE)
print(iwp2_linear, row.names = FALSE)
message("Saved directional overlap cache to ", output_path)
