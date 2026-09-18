#!/usr/bin/env Rscript

# Validate the exact pair-key source against all five frozen genotype samples.

project_root <- if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
  "."
} else {
  "coderepo-local"
}
experiment <- file.path(
  project_root,
  "code/revision_simulations/internal/r3_two_stage_functional_testing_pilot10"
)
source(file.path(experiment, "pilot10_helpers.R"))
source(file.path(
  project_root,
  "code/revision_simulations/shared/real_genotype_one_per_gene.R"
))
source(file.path(experiment, "genotype_cache_helpers.R"))

contract <- r3ts10_contract()
base_path <- file.path(
  project_root,
  "output/revision_simulations/shared/real_genotype_one_per_gene_J6362_pilot5/genotype_samples.rds"
)
pair_source_path <- file.path(
  project_root,
  "output/revision_simulations/internal/r3_two_stage_pilot10_development_20260904/pair_key_source.rds"
)
stopifnot(
  identical(r3ts10_sha256(base_path), contract$base_source_sha[["base_genotypes"]]),
  identical(
    r3ts10_sha256(pair_source_path),
    "a4b83e9828b0caeebffa226cc14100116599853a5c94f97ff4fcd2b8f8310145"
  )
)
base_cache <- readRDS(base_path)
pair_source <- readRDS(pair_source_path)
r3ts10g_validate_base_cache(base_cache)
pair_map <- r3ts10g_validate_pair_source(pair_source, base_cache)
stopifnot(
  nrow(pair_map) == 1009173L,
  length(unique(pair_map$gene_id)) == 6362L,
  length(unique(pair_map$variant_id)) == 745867L
)

configuration <- r3ts10g_configuration(
  base_cache = base_cache,
  pair_source_sha256 = r3ts10_sha256(pair_source_path),
  base_cache_sha256 = r3ts10_sha256(base_path),
  vcf_fingerprint = base_cache$configuration$vcf_fingerprint
)
stopifnot(
  identical(configuration$output_id, contract$genotype_id),
  identical(configuration$seed_list, contract$seeds),
  identical(configuration$base_seed_list, contract$base_seeds),
  identical(configuration$added_seed_list, contract$added_seeds)
)

cat("Pilot-10 pair-key and frozen-genotype checks passed.\n")

# A different R writer-version header must not change genotype identity.
sample <- base_cache$samples[[1L]]
serialized <- serialize(sample, NULL, version = 2L)
changed_header <- serialized
changed_header[7:10] <- writeBin(
  as.integer(4L * 65536L + 4L * 256L + 1L), raw(), size = 4L, endian = "big"
)
restored <- unserialize(changed_header)
stopifnot(!identical(serialized, changed_header), identical(sample, restored))
content_digest <- function(x) {
  genotype_content_md5(x$selection$pair_key, rownames(x$G), x$G)
}
stopifnot(identical(content_digest(sample), content_digest(restored)))
mutated <- base_cache
mutated$samples[[1L]]$G[1L, 1L] <- mutated$samples[[1L]]$G[1L, 1L] + 0.001
stopifnot(inherits(try(r3ts10g_validate_base_cache(mutated), silent = TRUE), "try-error"))
mutated_pairs <- pair_source
mutated_pairs$pair_ids[1:2] <- rev(mutated_pairs$pair_ids[1:2])
stopifnot(inherits(try(r3ts10g_validate_pair_source(mutated_pairs), silent = TRUE), "try-error"))
cat("Portable digest checks passed; altered dosage and ordered pair keys were rejected.\n")
