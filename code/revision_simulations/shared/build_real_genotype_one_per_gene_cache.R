#!/usr/bin/env Rscript

# Build the shared five-seed real-genotype cache used by formal R1 and R2.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/real_genotype_one_per_gene.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/real_genotype_one_per_gene.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[[1L]] == length(args)) return(default)
  args[hit[[1L]] + 1L]
}

as_flag <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

parse_seed_list <- function(x) {
  pieces <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  seeds <- suppressWarnings(as.integer(pieces))
  if (length(seeds) == 0L || anyNA(seeds) || anyDuplicated(seeds)) {
    stop("--seed-list must contain unique comma-separated integer seeds.")
  }
  seeds
}

resolve_input_path <- function(path, workflowr_root, argument_name) {
  candidates <- unique(c(path, file.path(workflowr_root, path)))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop(argument_name, " does not exist: ", path)
  }
  normalizePath(existing[[1L]], winslash = "/", mustWork = TRUE)
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
))

default_vcf_path <- file.path(
  workflowr_root,
  "..",
  "iPSC-data",
  "genotype-data",
  "YRI_genotype.vcf.gz"
)
default_pair_summary_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "eqtl_summary.rds"
)
vcf_path <- resolve_input_path(
  get_arg("--vcf-path", default_vcf_path),
  workflowr_root,
  "--vcf-path"
)
pair_summary_path <- resolve_input_path(
  get_arg("--pair-summary-path", default_pair_summary_path),
  workflowr_root,
  "--pair-summary-path"
)
seed_list <- parse_seed_list(get_arg(
  "--seed-list",
  "12345,22345,32345,42345,52345"
))
output_id <- get_arg("--output-id", "real_genotype_one_per_gene_J6362_pilot5")
expected_genes <- as.integer(get_arg("--expected-genes", "6362"))
expected_donors <- as.integer(get_arg("--expected-donors", "19"))
maf_min <- as.numeric(get_arg("--maf-min", "0.10"))
overwrite <- as_flag(get_arg("--overwrite", "false"))

expected_genes <- validate_integer_scalar(expected_genes, "expected_genes")
expected_donors <- validate_integer_scalar(expected_donors, "expected_donors")
if (!is.finite(maf_min) || maf_min < 0 || maf_min >= 0.5 || !nzchar(output_id)) {
  stop("Invalid cache-builder arguments.")
}

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "shared",
  output_id
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
cache_path <- file.path(output_dir, "genotype_samples.rds")

message("Reading paper gene-variant keys from ", pair_summary_path, ".")
pair_summary <- readRDS(pair_summary_path)
pair_ids <- names(pair_summary)
if (is.null(pair_ids) || length(pair_ids) == 0L) {
  stop("The paper pair summary must be a named object.")
}
rm(pair_summary)
invisible(gc())

pair_map <- parse_tested_pair_ids(pair_ids)
observed_genes <- length(unique(pair_map$gene_id))
observed_variants <- length(unique(pair_map$variant_id))
if (observed_genes != expected_genes) {
  stop(
    "The tested universe contains ", observed_genes,
    " genes; expected ", expected_genes, "."
  )
}
sample_ids <- read_vcf_sample_ids(vcf_path)
if (length(sample_ids) != expected_donors) {
  stop(
    "The VCF contains ", length(sample_ids),
    " donors; expected ", expected_donors, "."
  )
}

configuration <- list(
  format_version = 1L,
  output_id = output_id,
  n_tested_pairs = nrow(pair_map),
  n_tested_genes = observed_genes,
  n_tested_variants = observed_variants,
  n_genes = expected_genes,
  n_donors = expected_donors,
  seed_list = seed_list,
  selection_rule = "one uniformly sampled tested variant per gene",
  repeated_variant_rule = "retain repeated rsIDs assigned to different genes",
  dosage_field = "DS",
  maf_min = maf_min,
  vcf_fingerprint = artifact_fingerprint(vcf_path),
  pair_summary_fingerprint = artifact_fingerprint(pair_summary_path)
)

validate_cache <- function(cache) {
  if (!is.list(cache) ||
      !all(c("configuration", "sample_ids", "samples") %in% names(cache)) ||
      !isTRUE(all.equal(cache$configuration, configuration)) ||
      !identical(cache$sample_ids, sample_ids) ||
      !identical(names(cache$samples), as.character(seed_list))) {
    return(FALSE)
  }
  all(vapply(seed_list, function(seed) {
    sample <- cache$samples[[as.character(seed)]]
    validated <- try(
      validate_real_genotype_sample(
        sample,
        expected_genes = expected_genes,
        expected_donors = expected_donors,
        maf_min = maf_min
      ),
      silent = TRUE
    )
    if (inherits(validated, "try-error")) return(FALSE)
    identical(
      validated$genotype_digest,
      object_md5(list(
        pair_key = validated$selection$pair_key,
        sample_ids = rownames(validated$G),
        G = validated$G
      ))
    )
  }, logical(1)))
}

if (file.exists(cache_path) && !overwrite) {
  cache <- readRDS(cache_path)
  if (!validate_cache(cache)) {
    stop("The existing shared genotype cache does not match the requested design.")
  }
  message("Reusing validated shared genotype cache: ", cache_path)
  quit(save = "no", status = 0L)
}

message("Selecting one tested variant per gene for each formal seed.")
selections <- lapply(seed_list, function(seed) {
  sample_one_tested_variant_per_gene(pair_ids, seed = seed)
})
names(selections) <- as.character(seed_list)
target_variants <- unique(unlist(
  lapply(selections, `[[`, "variant_id"),
  use.names = FALSE
))
message(
  "Streaming ", format(length(target_variants), big.mark = ","),
  " selected unique variant IDs from the YRI VCF."
)
extracted <- extract_target_vcf_dosages(
  vcf_path = vcf_path,
  target_variants = target_variants,
  work_dir = tempdir()
)
if (extracted$matched_variant_count != length(target_variants)) {
  stop(
    "The VCF matched ", extracted$matched_variant_count,
    " of ", length(target_variants), " selected variants."
  )
}

samples <- lapply(seed_list, function(seed) {
  message("Assembling the real-genotype matrix for seed ", seed, ".")
  sample <- assemble_one_per_gene_genotype_sample(
    selection = selections[[as.character(seed)]],
    extracted = extracted,
    seed = seed
  )
  sample <- validate_real_genotype_sample(
    sample,
    expected_genes = expected_genes,
    expected_donors = expected_donors,
    maf_min = maf_min
  )
  sample$genotype_digest <- object_md5(list(
    pair_key = sample$selection$pair_key,
    sample_ids = rownames(sample$G),
    G = sample$G
  ))
  sample
})
names(samples) <- as.character(seed_list)

existing_seed_12345_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "one_variant_per_gene_refit_seed12345",
  "selection.csv"
)
if (12345L %in% seed_list && file.exists(existing_seed_12345_path)) {
  previous <- read.csv(existing_seed_12345_path, stringsAsFactors = FALSE)
  current <- samples[["12345"]]$selection
  if (!identical(previous$pair_key, current$pair_key)) {
    stop("Seed 12345 does not reproduce the existing one-per-gene selection.")
  }
}

cache <- list(
  configuration = configuration,
  sample_ids = sample_ids,
  samples = samples
)
if (!validate_cache(cache)) {
  stop("The newly assembled shared genotype cache failed validation.")
}
saveRDS(cache, cache_path, version = 2)

selection_table <- do.call(rbind, lapply(seed_list, function(seed) {
  selection <- samples[[as.character(seed)]]$selection
  data.frame(seed = seed, selection, stringsAsFactors = FALSE)
}))
selection_summary <- do.call(rbind, lapply(samples, `[[`, "selection_summary"))
maf_table <- do.call(rbind, lapply(seed_list, function(seed) {
  info <- samples[[as.character(seed)]]$variant_info
  data.frame(
    seed = seed,
    gene_id = info$gene_id,
    variant_id = info$variant_id,
    pair_key = info$unit_id,
    chromosome = info$chromosome,
    position = info$position,
    observed_maf = info$observed_maf,
    genotype_sd = info$genotype_sd,
    stringsAsFactors = FALSE
  )
}))
input_provenance <- data.frame(
  artifact = c("YRI genotype VCF", "Paper gene-variant summary", "Shared genotype cache"),
  file_name = c(
    configuration$vcf_fingerprint$file_name,
    configuration$pair_summary_fingerprint$file_name,
    basename(cache_path)
  ),
  size_bytes = c(
    configuration$vcf_fingerprint$size_bytes,
    configuration$pair_summary_fingerprint$size_bytes,
    file.info(cache_path)$size
  ),
  md5 = c(
    configuration$vcf_fingerprint$md5,
    configuration$pair_summary_fingerprint$md5,
    unname(tools::md5sum(cache_path))
  ),
  stringsAsFactors = FALSE
)

write_csv(selection_table, file.path(output_dir, "selection.csv"))
write_csv(selection_summary, file.path(output_dir, "seed_selection_summary.csv"))
write_csv(maf_table, file.path(output_dir, "selected_variant_metadata.csv"))
write_csv(input_provenance, file.path(output_dir, "input_provenance.csv"))

message("Wrote validated shared genotype cache: ", cache_path)
print(selection_summary)
