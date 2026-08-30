#!/usr/bin/env Rscript

# Configure the formal full-universe R3 functional calibration refresh with
# the canonical IWP1 temporal mixture and paired Raw/BF posterior summaries.

file_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (length(file_argument) != 1L) {
  stop("Could not resolve the full-universe R3 entry point.", call. = FALSE)
}
runner_path <- normalizePath(
  sub("^--file=", "", file_argument[[1L]]),
  winslash = "/",
  mustWork = TRUE
)
wrapper_core_path <- file.path(
  dirname(runner_path),
  "24_r3_full_universe_iwp1_mixture_core.R"
)
if (!file.exists(wrapper_core_path) || dir.exists(wrapper_core_path)) {
  stop(
    "The full-universe R3 wrapper core is missing beside the entry point.",
    call. = FALSE
  )
}
wrapper_core_path <- normalizePath(
  wrapper_core_path,
  winslash = "/",
  mustWork = TRUE
)

Sys.setenv(
  FASH_R3_RESULT_ID = paste0(
    "r3_real_genotype_one_per_gene_J6362_",
    paste0(
      "matched_functional_open_middle_3_12_center_aligned_",
      "iwp1_geometry_mixture_"
    ),
    "relative_location_clearance_",
    "full_universe_",
    "paired_posterior_fashr0143_pilot5"
  ),
  FASH_R3_MIDDLE_WINDOW = "3,12",
  FASH_R3_MIDDLE_BOUNDARY = "open",
  FASH_R3_TEMPORAL_CATEGORY_PROBS = "0.29,0.42,0.29",
  FASH_R3_WRAPPER_CORE = wrapper_core_path
)

source(wrapper_core_path, local = globalenv(), chdir = FALSE)
