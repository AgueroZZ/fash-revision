#!/usr/bin/env Rscript

# Configure the formal R3 calibration refresh with the exhaustive open-Middle
# estimand, full-support raised-cosine centers, and a canonical IWP1
# temporal truth-category mixture. Historical balanced-category outputs remain
# immutable and serve as a separate prior-shift stress test.

file_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (length(file_argument) != 1L) {
  stop("Could not resolve the prior-geometry R3 entry point.", call. = FALSE)
}
runner_path <- normalizePath(
  sub("^--file=", "", file_argument[[1L]]),
  winslash = "/",
  mustWork = TRUE
)
wrapper_core_path <- file.path(
  dirname(runner_path),
  "18_r3_prior_geometry_core.R"
)
if (!file.exists(wrapper_core_path) || dir.exists(wrapper_core_path)) {
  stop(
    "The prior-geometry R3 wrapper core is missing beside the entry point.",
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
    "matched_functional_open_middle_3_12_full_support_",
    "iwp1_geometry_mixture_relative_clearance_main_effect_",
    "fashr0143_pilot5"
  ),
  FASH_R3_MIDDLE_WINDOW = "3,12",
  FASH_R3_MIDDLE_BOUNDARY = "open",
  FASH_R3_TEMPORAL_CATEGORY_PROBS = "0.29,0.42,0.29",
  FASH_R3_WRAPPER_CORE = wrapper_core_path
)

source(wrapper_core_path, local = globalenv(), chdir = FALSE)
