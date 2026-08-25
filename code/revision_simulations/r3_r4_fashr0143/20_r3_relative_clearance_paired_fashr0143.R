#!/usr/bin/env Rscript

# Configure the formal R3 refresh with scale-adaptive location clearance and
# common-random-number posterior summaries for the Raw and BF fits.

file_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (length(file_argument) != 1L) {
  stop("Could not resolve the relative-clearance R3 entry point.", call. = FALSE)
}
runner_path <- normalizePath(
  sub("^--file=", "", file_argument[[1L]]),
  winslash = "/",
  mustWork = TRUE
)
wrapper_core_path <- file.path(
  dirname(runner_path),
  "20_r3_relative_clearance_paired_core.R"
)
if (!file.exists(wrapper_core_path) || dir.exists(wrapper_core_path)) {
  stop(
    "The relative-clearance R3 wrapper core is missing beside the entry point.",
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
    "iwp1_geometry_mixture_relative_location_clearance_",
    "paired_posterior_fashr0143_pilot5"
  ),
  FASH_R3_MIDDLE_WINDOW = "3,12",
  FASH_R3_MIDDLE_BOUNDARY = "open",
  FASH_R3_TEMPORAL_CATEGORY_PROBS = "0.29,0.42,0.29",
  FASH_R3_WRAPPER_CORE = wrapper_core_path
)

source(wrapper_core_path, local = globalenv(), chdir = FALSE)
