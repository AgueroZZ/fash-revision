#!/usr/bin/env Rscript

# Configure the formal R3 wrapper for the manuscript-consistent open-Middle
# definition. The shared wrapper core performs all package, source, input,
# replicate, summary, manifest, and atomic-promotion validation.

file_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (length(file_argument) != 1L) {
  stop("Could not resolve the open-Middle R3 entry point.", call. = FALSE)
}
runner_path <- normalizePath(
  sub("^--file=", "", file_argument[[1L]]),
  winslash = "/",
  mustWork = TRUE
)
wrapper_core_candidates <- file.path(
  dirname(runner_path),
  c("15_r3_open_middle_core.R", "13_r3_fashr0143.R")
)
existing_wrapper_cores <- wrapper_core_candidates[
  file.exists(wrapper_core_candidates) & !dir.exists(wrapper_core_candidates)
]
if (length(existing_wrapper_cores) == 0L) {
  stop(
    "The R3 wrapper core is missing beside the open-Middle entry point.",
    call. = FALSE
  )
}
wrapper_core_path <- existing_wrapper_cores[[1L]]
wrapper_core_path <- normalizePath(
  wrapper_core_path,
  winslash = "/",
  mustWork = TRUE
)

Sys.setenv(
  FASH_R3_RESULT_ID = paste0(
    "r3_real_genotype_one_per_gene_J6362_",
    "matched_functional_open_middle_3_12_",
    "relative_clearance_main_effect_fashr0143_pilot5"
  ),
  FASH_R3_MIDDLE_WINDOW = "3,12",
  FASH_R3_MIDDLE_BOUNDARY = "open",
  FASH_R3_WRAPPER_CORE = wrapper_core_path
)

source(wrapper_core_path, local = globalenv(), chdir = FALSE)
