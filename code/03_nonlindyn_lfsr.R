#!/usr/bin/env Rscript

# Canonical order-2 nonlinear dynamic-eQTL functional-classification entrypoint.
# The implementation is shared with the order-1 runner so functional
# definitions, posterior sampling, provenance, and atomic-output safeguards
# remain identical across both analyses.

command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
if (length(file_argument) != 1L) {
  stop("Could not identify the 03_nonlindyn_lfsr.R script path.")
}
script_path <- normalizePath(
  sub("^--file=", "", file_argument),
  winslash = "/",
  mustWork = TRUE
)

Sys.setenv(FASH_LFSR_ANALYSIS_ORDER = "2")
source(file.path(dirname(script_path), "02_dyn_lfsr.R"), chdir = FALSE)
