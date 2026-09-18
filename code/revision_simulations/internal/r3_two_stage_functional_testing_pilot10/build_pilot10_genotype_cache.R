#!/usr/bin/env Rscript

# Build the separate ten-seed real-genotype cache for the two-stage R3 analysis.

get_arg <- function(name, default = NULL) {
  arguments <- commandArgs(trailingOnly = TRUE)
  index <- match(name, arguments)
  if (is.na(index)) return(default)
  if (index == length(arguments)) stop("Missing value for ", name, call. = FALSE)
  arguments[[index + 1L]]
}

find_project_root <- function() {
  supplied <- get_arg("--project-root", ".")
  root <- normalizePath(supplied, winslash = "/", mustWork = TRUE)
  required <- file.path(
    root,
    "code/revision_simulations/shared/real_genotype_one_per_gene.R"
  )
  if (!file.exists(required)) {
    stop("Could not find the workflowr project root.", call. = FALSE)
  }
  root
}

root <- find_project_root()
experiment <- Sys.getenv(
  "FASH_R3_PILOT10_CODE",
  unset = file.path(
    root,
    "code/revision_simulations/internal/r3_two_stage_functional_testing_pilot10"
  )
)
source(file.path(experiment, "pilot10_helpers.R"))
source(file.path(
  root,
  "code/revision_simulations/shared/real_genotype_one_per_gene.R"
))
source(file.path(experiment, "genotype_cache_helpers.R"))

mode <- get_arg("--mode", "build")
r3ts10_assert(
  mode %in% c("pair-source", "preflight", "build"),
  "Mode must be pair-source, preflight, or build."
)
contract <- r3ts10_contract()
base_cache_path <- get_arg(
  "--base-cache",
  Sys.getenv(
    "FASH_R3_BASE_GENOTYPE_CACHE",
    unset = file.path(
      root,
      "output/revision_simulations/shared/real_genotype_one_per_gene_J6362_pilot5/genotype_samples.rds"
    )
  )
)
pair_source_path <- get_arg(
  "--pair-source",
  Sys.getenv(
    "FASH_R3_PILOT10_PAIR_SOURCE",
    unset = file.path(
      root,
      "output/revision_simulations/internal/r3_two_stage_pilot10_development_20260904/pair_key_source.rds"
    )
  )
)
output_dir <- get_arg(
  "--output-dir",
  Sys.getenv(
    "FASH_R3_PILOT10_GENOTYPE_DIR",
    unset = file.path(
      root,
      "output/revision_simulations/shared",
      contract$genotype_id
    )
  )
)
vcf_path <- get_arg(
  "--vcf-path",
  Sys.getenv("FASH_R3_VCF", unset = "")
)

r3ts10_assert(file.exists(base_cache_path), "The frozen five-seed genotype cache is missing.")
r3ts10_assert(
  identical(r3ts10_sha256(base_cache_path), contract$base_source_sha[["base_genotypes"]]),
  "The frozen five-seed genotype cache hash changed."
)
base_cache <- readRDS(base_cache_path)
r3ts10g_validate_base_cache(base_cache)

if (mode == "pair-source") {
  pair_summary_path <- get_arg("--pair-summary", "")
  r3ts10_assert(file.exists(pair_summary_path), "The frozen pair summary is missing.")
  fingerprint <- artifact_fingerprint(pair_summary_path)
  r3ts10_assert(
    isTRUE(all.equal(fingerprint, base_cache$configuration$pair_summary_fingerprint)),
    "The supplied pair summary does not match the frozen genotype-cache provenance."
  )
  message("Reading ordered pair keys from the frozen pair summary.")
  pair_summary <- readRDS(pair_summary_path)
  pair_ids <- names(pair_summary)
  rm(pair_summary)
  invisible(gc())
  pair_source <- r3ts10g_make_pair_source(pair_ids, fingerprint)
  r3ts10g_validate_pair_source(pair_source, base_cache)
  r3ts10_write_rds(pair_source, pair_source_path, version = 2L)
  message("Wrote exact ordered pair-key source: ", pair_source_path)
  message("SHA-256: ", r3ts10_sha256(pair_source_path))
  quit(save = "no", status = 0L)
}

r3ts10_assert(file.exists(pair_source_path), "The exact pair-key source is missing.")
pair_source <- readRDS(pair_source_path)
r3ts10g_validate_pair_source(pair_source, base_cache)
pair_source_sha256 <- r3ts10_sha256(pair_source_path)
base_cache_sha256 <- r3ts10_sha256(base_cache_path)
r3ts10_assert(nzchar(vcf_path) && file.exists(vcf_path), "The server YRI VCF is missing.")
vcf_fingerprint <- artifact_fingerprint(vcf_path)
expected_vcf <- base_cache$configuration$vcf_fingerprint
r3ts10_assert(
  identical(vcf_fingerprint$file_name, expected_vcf$file_name) &&
    identical(vcf_fingerprint$size_bytes, expected_vcf$size_bytes) &&
    identical(vcf_fingerprint$md5, expected_vcf$md5) &&
    identical(read_vcf_sample_ids(vcf_path), base_cache$sample_ids),
  "The server VCF or donor order does not match the frozen source."
)

cache_path <- file.path(output_dir, "genotype_samples.rds")
flag_path <- file.path(output_dir, "complete.flag")
if (file.exists(cache_path)) {
  cache <- readRDS(cache_path)
  validation <- r3ts10g_validate_cache(
    cache, base_cache, pair_source,
    pair_source_sha256 = pair_source_sha256,
    base_cache_sha256 = base_cache_sha256
  )
  message("Validated existing ten-seed genotype cache: ", cache_path)
} else if (mode == "preflight") {
  message(
    "Preflight passed: exact tested universe, five retained selections, ",
    "VCF fingerprint, and donor order."
  )
  quit(save = "no", status = 0L)
} else {
  message("Extracting VCF dosage for the five added seed selections.")
  cache <- r3ts10g_build_cache(
    base_cache = base_cache,
    pair_source = pair_source,
    vcf_path = vcf_path,
    pair_source_sha256 = pair_source_sha256,
    base_cache_sha256 = base_cache_sha256,
    work_dir = Sys.getenv("TMPDIR", unset = tempdir())
  )
  validation <- r3ts10g_validate_cache(
    cache, base_cache, pair_source,
    pair_source_sha256 = pair_source_sha256,
    base_cache_sha256 = base_cache_sha256
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  r3ts10_write_rds(cache, cache_path, version = 2L)
  message("Wrote separate ten-seed genotype cache: ", cache_path)
}

if (mode == "preflight") {
  message("Preflight passed, including validation of the existing ten-seed cache.")
  quit(save = "no", status = 0L)
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  validation,
  file.path(output_dir, "seed_validation.csv"),
  row.names = FALSE
)
r3ts10_write_rds(cache$configuration, file.path(output_dir, "configuration.rds"))
artifacts <- c("genotype_samples.rds", "configuration.rds", "seed_validation.csv")
manifest <- list(
  schema = "r3-two-stage-genotype-cache-manifest-v1",
  genotype_id = contract$genotype_id,
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  pair_key_source_sha256 = pair_source_sha256,
  base_cache_sha256 = base_cache_sha256,
  artifact_sha256 = setNames(
    vapply(file.path(output_dir, artifacts), r3ts10_sha256, character(1)),
    artifacts
  )
)
r3ts10_write_rds(manifest, file.path(output_dir, "manifest.rds"))
writeLines(
  c(
    paste0("genotype_id=", contract$genotype_id),
    "seeds=10",
    "base_samples_retained_exactly=5",
    "added_vcf_extractions=5",
    "genes_per_seed=6362",
    "donors=19",
    "maf_min=0.10"
  ),
  flag_path
)
message("Completed and validated ten-seed genotype cache: ", output_dir)
