#!/usr/bin/env Rscript
# Validates the base-R DS reader added to Midway3 script
# scripts/22_per_unit_design_correlation.R against the independently written
# read_selected_vcf_dosages(), on the real YRI VCF.

source("code/revision_simulations/internal/covariance_estimation/donor_null_permutation_helpers.R")

# Extract the reader from the staged Midway3 script so the tested code is
# literally the code that will run on the cluster.
staged <- file.path(Sys.getenv("HOME"), "Desktop", "midway3", "fash",
                    "scripts", "22_per_unit_design_correlation.R")
lines <- readLines(staged)
start <- grep("^extract_colon_field <- function", lines)
end <- grep("^# ---+ inputs ---", lines)
stopifnot(length(start) == 1L, length(end) == 1L, end > start)
log_message <- function(...) invisible(NULL)
eval(parse(text = paste(lines[start:(end - 1L)], collapse = "\n")))
message("Loaded read_vcf_dosages() from the staged cluster script, lines ",
        start, "-", end - 1L, ".")

ok <- function(condition, label) {
  if (!isTRUE(condition)) stop("FAILED: ", label)
  message("ok  ", label)
}

vcf <- file.path("..", "iPSC-data", "genotype-data", "YRI_genotype.vcf.gz")
meta <- utils::read.csv(paste0("output/revision_simulations/internal/",
                              "correlated_permutation_discoveries/",
                              "panel_topgene_metadata.csv"),
                        stringsAsFactors = FALSE)
set.seed(20260826)
probe <- unique(meta$variant_id)[sort(sample(length(unique(meta$variant_id)), 250L))]

started <- Sys.time()
new_reader <- read_vcf_dosages(vcf, probe)
new_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
started <- Sys.time()
reference <- read_selected_vcf_dosages(vcf, probe)
reference_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

ok(identical(dim(new_reader), dim(t(reference))) ||
     identical(dim(new_reader), dim(reference)),
   "both readers return the same shape")
message("    new ", paste(dim(new_reader), collapse = " x "),
        ", reference ", paste(dim(reference), collapse = " x "))

# read_selected_vcf_dosages returns donors x variants; transpose to compare.
reference_aligned <- t(reference)[rownames(new_reader), colnames(new_reader),
                                 drop = FALSE]
ok(identical(dimnames(new_reader), dimnames(reference_aligned)),
   "variant and donor names agree exactly")
ok(max(abs(new_reader - reference_aligned)) == 0,
   "every dosage is bit-identical between the two readers")
ok(nrow(new_reader) == 250L && ncol(new_reader) == 19L,
   "250 requested variants and 19 donors returned")
ok(all(new_reader >= 0 & new_reader <= 2),
   "every dosage lies in [0, 2]")
ok(identical(rownames(new_reader), probe),
   "rows come back in the requested order")

message("    new reader ", round(new_seconds, 1), " s, reference ",
        round(reference_seconds, 1), " s for the same full-file scan")

# A missing variant must fail loudly rather than return a short matrix.
failed <- tryCatch({
  read_vcf_dosages(vcf, c(probe[1:3], "rs_not_a_real_variant_xyz"))
  FALSE
}, error = function(condition) grepl("absent from the VCF",
                                    conditionMessage(condition)))
ok(isTRUE(failed), "an absent variant raises a clear error")

message("All VCF dosage reader checks passed.")
