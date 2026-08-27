#!/usr/bin/env Rscript

# Configure the center-aligned formal R3 refresh with scale-adaptive location
# clearance and common-random-number posterior summaries for Raw and BF.

file_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (length(file_argument) != 1L) {
  stop("Could not resolve the center-aligned R3 entry point.", call. = FALSE)
}
runner_path <- normalizePath(
  sub("^--file=", "", file_argument[[1L]]),
  winslash = "/",
  mustWork = TRUE
)
wrapper_core_path <- file.path(
  dirname(runner_path),
  "21_r3_center_aligned_relative_location_paired_core.R"
)
if (!file.exists(wrapper_core_path) || dir.exists(wrapper_core_path)) {
  stop(
    "The center-aligned R3 wrapper core is missing beside the entry point.",
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
    "matched_functional_open_middle_3_12_center_aligned_equal_cells_",
    "relative_location_clearance_",
    "paired_posterior_fashr0143_pilot5"
  ),
  FASH_R3_MIDDLE_WINDOW = "3,12",
  FASH_R3_MIDDLE_BOUNDARY = "open",
  FASH_R3_TEMPORAL_CATEGORY_PROBS = paste(
    rep("0.3333333333333333", 3),
    collapse = ","
  ),
  FASH_R3_WRAPPER_CORE = wrapper_core_path
)

source(wrapper_core_path, local = globalenv(), chdir = FALSE)
