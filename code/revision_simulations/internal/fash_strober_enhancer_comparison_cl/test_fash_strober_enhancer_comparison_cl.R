#!/usr/bin/env Rscript

# Focused contract tests for the CL-PC R6 enhancer analysis.

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
shared_directory <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "fash_strober_enhancer_comparison"
)
source(file.path(shared_directory, "fash_strober_enhancer_helpers.R"))

matched_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "matched_fash_linear_real_data_cl_fashr_0_1_43"
)
matched_cache_path <- file.path(matched_directory, "analysis_cache.rds")
adjustment_path <- file.path(matched_directory, "cl_iwp1_bf_adjustment.rds")
statistics_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fash_linear_real_data_ablation_cl",
  "sufficient_statistics.rds"
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
required_paths <- c(
  matched_cache_path,
  adjustment_path,
  statistics_path,
  linear_path,
  quadratic_path
)
stopifnot(all(file.exists(required_paths)))

matched_cache <- readRDS(matched_cache_path)
adjustment <- readRDS(adjustment_path)
statistics_cache <- readRDS(statistics_path)
pair_keys <- as.character(statistics_cache$statistics$unit_id)
stopifnot(
  identical(matched_cache$configuration$package_version, "0.1.43"),
  identical(
    matched_cache$configuration$package_sha,
    "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
  ),
  all(matched_cache$validation$passed),
  identical(adjustment$package_version, "0.1.43"),
  identical(
    adjustment$package_sha,
    "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
  ),
  length(pair_keys) == 1009173L,
  !anyDuplicated(pair_keys),
  length(adjustment$lfdr) == length(pair_keys)
)

cl <- derive_fash_discovery_sets(pair_keys, adjustment$lfdr, alpha = 0.05)
stopifnot(identical(
  c(
    nrow(cl$all_pairs),
    length(cl$all_variants),
    nrow(cl$lead_pairs),
    length(cl$lead_variants)
  ),
  c(5395L, 5363L, 686L, 686L)
))

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
linear <- derive_strober_discovery_sets(linear_result, pair_keys, alpha = 0.05)
quadratic <- derive_strober_discovery_sets(
  quadratic_result,
  pair_keys,
  alpha = 0.05
)
strober_union <- union(linear$all_variants, quadratic$all_variants)
cl_only_all <- cl$all_pairs[
  !cl$all_pairs$variant_id %in% strober_union,
  ,
  drop = FALSE
]
cl_only_lead <- cl$lead_pairs[
  !cl$lead_pairs$variant_id %in% strober_union,
  ,
  drop = FALSE
]
stopifnot(
  identical(
    c(
      nrow(cl_only_all),
      length(unique(cl_only_all$variant_id)),
      nrow(cl_only_lead),
      length(unique(cl_only_lead$variant_id))
    ),
    c(3462L, 3456L, 406L, 406L)
  ),
  !any(cl_only_all$variant_id %in% strober_union),
  !any(cl_only_lead$variant_id %in% strober_union),
  all(cl_only_lead$pair_key %in% cl$lead_pairs$pair_key)
)

retained_cache_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fash_strober_enhancer_comparison_cl_fashr0143",
  "analysis_cache.rds"
)
if (file.exists(retained_cache_path)) {
  retained <- readRDS(retained_cache_path)
  expected_sets <- c(
    "cl_all", "cl_lead", "linear_all", "linear_lead",
    "quadratic_all", "quadratic_lead", "cl_only_all", "cl_only_lead"
  )
  observed <- retained$discovery_set_summary
  observed <- observed[match(expected_sets, observed$discovery_set), ]
  stopifnot(
    identical(
      retained$configuration$analysis_id,
      "revision_internal_fash_strober_enhancer_comparison_cl_fashr0143"
    ),
    identical(retained$set_metadata$discovery_set, expected_sets),
    identical(
      as.integer(observed$pair_count),
      c(5395L, 686L, 5404L, 550L, 6824L, 693L, 3462L, 406L)
    ),
    identical(
      as.integer(observed$unique_variant_count),
      c(5363L, 686L, 5387L, 548L, 6797L, 690L, 3456L, 406L)
    ),
    retained$strict_exclusivity_audit$overlap_with_strober_union == 0L,
    retained$maximum_absolute_smd < 0.10
  )
}

message("All CL-PC R6 enhancer comparison tests passed.")
