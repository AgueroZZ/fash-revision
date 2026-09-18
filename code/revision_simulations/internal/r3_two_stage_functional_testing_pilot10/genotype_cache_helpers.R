# Extend the frozen five-seed real-genotype cache without changing its samples.

# Canonical content digests exclude R serialization headers and compression.
r3ts10g_expected_content_md5 <- function() {
  setNames(c(
    "526a7318aa2af901e09252f5a6ca3c46", "517faa30d5218a956f1be84f2567369c",
    "9b9be3205d7db54dac31763492bcb2eb", "9ef73aa94a061df868b1a951fe495d9f",
    "7dab159b8453e2f66188ae313bfbd611", "da61aecbc01dab1e95ee2b944e347fd0",
    "29adfdc22a0438e60ad078c2fdfe02df", "092eca9c43981ebec4c6124dc2a5b41d",
    "eca57a24485e853cfbd2eed856e8cfc1", "1906f2879f103dd565248c0e30aaf3af"
  ), as.character(r3ts10_contract()$seeds))
}

r3ts10g_pair_content_md5 <- function(pair_ids) {
  path <- tempfile("r3-pair-keys-", fileext = ".txt")
  on.exit(unlink(path), add = TRUE)
  writeLines(enc2utf8(pair_ids), path, useBytes = TRUE)
  unname(tools::md5sum(path))
}

r3ts10g_require_dependencies <- function() {
  required <- c(
    "r3ts10_contract", "r3ts10_assert", "sample_one_tested_variant_per_gene",
    "parse_tested_pair_ids", "read_vcf_sample_ids", "extract_target_vcf_dosages",
    "assemble_one_per_gene_genotype_sample", "validate_real_genotype_sample",
    "artifact_fingerprint", "object_md5", "genotype_content_md5"
  )
  r3ts10_assert(
    all(vapply(required, exists, logical(1), inherits = TRUE)),
    "Pilot-10 genotype-cache dependencies are incomplete."
  )
}

r3ts10g_base_cache_id <- function() {
  "real_genotype_one_per_gene_J6362_pilot5"
}

r3ts10g_pair_source_schema <- function() {
  "r3-two-stage-pair-key-source-v1"
}

r3ts10g_make_pair_source <- function(pair_ids, source_fingerprint) {
  r3ts10g_require_dependencies()
  out <- list(
    schema = r3ts10g_pair_source_schema(),
    source_fingerprint = source_fingerprint,
    pair_ids = as.character(pair_ids),
    pair_ids_md5 = object_md5(as.character(pair_ids))
  )
  r3ts10g_validate_pair_source(out)
  out
}

r3ts10g_validate_pair_source <- function(pair_source, base_cache = NULL) {
  r3ts10g_require_dependencies()
  r3ts10_assert(
    is.list(pair_source) &&
      identical(pair_source$schema, r3ts10g_pair_source_schema()) &&
      is.character(pair_source$pair_ids) &&
      length(pair_source$pair_ids) == 1009173L &&
      !anyNA(pair_source$pair_ids) &&
      all(nzchar(pair_source$pair_ids)) &&
      !anyDuplicated(pair_source$pair_ids) &&
      identical(r3ts10g_pair_content_md5(pair_source$pair_ids),
                "0868f2c0cb004261bc300845f4e39580"),
    "The pair-key source is incomplete or has changed."
  )
  pair_map <- parse_tested_pair_ids(pair_source$pair_ids)
  r3ts10_assert(
    nrow(pair_map) == 1009173L &&
      length(unique(pair_map$gene_id)) == 6362L &&
      length(unique(pair_map$variant_id)) == 745867L,
    "The pair-key source does not reproduce the tested universe."
  )
  if (!is.null(base_cache)) {
    r3ts10_assert(
      isTRUE(all.equal(
        pair_source$source_fingerprint,
        base_cache$configuration$pair_summary_fingerprint
      )),
      "The pair-key source does not come from the frozen pair summary."
    )
    for (seed in r3ts10_contract()$base_seeds) {
      selected <- sample_one_tested_variant_per_gene(pair_source$pair_ids, seed)
      r3ts10_assert(
        identical(selected, base_cache$samples[[as.character(seed)]]$selection),
        paste("The pair-key source does not reproduce base seed", seed, "exactly.")
      )
    }
  }
  invisible(pair_map)
}

r3ts10g_validate_base_cache <- function(base_cache) {
  r3ts10g_require_dependencies()
  contract <- r3ts10_contract()
  configuration <- base_cache$configuration
  r3ts10_assert(
    is.list(base_cache) &&
      identical(names(base_cache$samples), as.character(contract$base_seeds)) &&
      identical(configuration$output_id, r3ts10g_base_cache_id()) &&
      identical(configuration$seed_list, contract$base_seeds) &&
      identical(configuration$n_tested_pairs, 1009173L) &&
      identical(configuration$n_tested_genes, 6362L) &&
      identical(configuration$n_tested_variants, 745867L) &&
      identical(configuration$n_genes, contract$n_units) &&
      identical(configuration$n_donors, 19L) &&
      identical(configuration$maf_min, 0.10) &&
      identical(configuration$dosage_field, "DS") &&
      identical(configuration$vcf_fingerprint$file_name, "YRI_genotype.vcf.gz") &&
      identical(configuration$vcf_fingerprint$size_bytes, 118088371) &&
      identical(configuration$vcf_fingerprint$md5, "4f9eb383ce3512d867b42ab806d451a8") &&
      identical(configuration$pair_summary_fingerprint$file_name, "eqtl_summary.rds") &&
      identical(configuration$pair_summary_fingerprint$size_bytes, 93591244) &&
      identical(configuration$pair_summary_fingerprint$md5, "23e7f6a0093309059424207872dfa1e0") &&
      identical(base_cache$sample_ids, rownames(base_cache$samples[[1L]]$G)),
    "The frozen five-seed genotype cache has an unexpected contract."
  )
  for (seed in contract$base_seeds) {
    sample <- validate_real_genotype_sample(
      base_cache$samples[[as.character(seed)]],
      expected_genes = contract$n_units,
      expected_donors = configuration$n_donors,
      maf_min = configuration$maf_min
    )
    expected_digest <- r3ts10g_expected_content_md5()[[as.character(seed)]]
    actual_digest <- genotype_content_md5(
      sample$selection$pair_key, rownames(sample$G), sample$G
    )
    r3ts10_assert(
      identical(actual_digest, expected_digest),
      paste("The frozen genotype digest changed for seed", seed, ".")
    )
  }
  invisible(TRUE)
}

r3ts10g_configuration <- function(base_cache, pair_source_sha256,
                                    base_cache_sha256, vcf_fingerprint) {
  contract <- r3ts10_contract()
  base <- base_cache$configuration
  list(
    format_version = 2L,
    output_id = contract$genotype_id,
    parent_output_id = base$output_id,
    parent_cache_sha256 = unname(base_cache_sha256),
    pair_key_source_sha256 = unname(pair_source_sha256),
    n_tested_pairs = base$n_tested_pairs,
    n_tested_genes = base$n_tested_genes,
    n_tested_variants = base$n_tested_variants,
    n_genes = base$n_genes,
    n_donors = base$n_donors,
    seed_list = contract$seeds,
    base_seed_list = contract$base_seeds,
    added_seed_list = contract$added_seeds,
    selection_rule = base$selection_rule,
    repeated_variant_rule = base$repeated_variant_rule,
    dosage_field = base$dosage_field,
    maf_min = base$maf_min,
    vcf_fingerprint = vcf_fingerprint,
    pair_summary_fingerprint = base$pair_summary_fingerprint,
    reuse_rule = "retain all five base samples exactly and extract dosage only for added seeds"
  )
}

r3ts10g_validate_cache <- function(cache, base_cache, pair_source,
                                    pair_source_sha256, base_cache_sha256) {
  r3ts10g_require_dependencies()
  r3ts10g_validate_base_cache(base_cache)
  r3ts10g_validate_pair_source(pair_source, base_cache)
  contract <- r3ts10_contract()
  r3ts10_assert(
    is.list(cache) && all(c("configuration", "sample_ids", "samples") %in% names(cache)) &&
      identical(cache$configuration$output_id, contract$genotype_id) &&
      identical(cache$configuration$parent_output_id, r3ts10g_base_cache_id()) &&
      identical(cache$configuration$parent_cache_sha256, unname(base_cache_sha256)) &&
      identical(cache$configuration$pair_key_source_sha256, unname(pair_source_sha256)) &&
      identical(cache$configuration$seed_list, contract$seeds) &&
      identical(cache$configuration$base_seed_list, contract$base_seeds) &&
      identical(cache$configuration$added_seed_list, contract$added_seeds) &&
      identical(cache$configuration$n_tested_pairs, 1009173L) &&
      identical(cache$configuration$n_tested_genes, 6362L) &&
      identical(cache$configuration$n_tested_variants, 745867L) &&
      identical(cache$configuration$n_genes, contract$n_units) &&
      identical(cache$configuration$n_donors, 19L) &&
      identical(cache$configuration$maf_min, 0.10) &&
      identical(cache$configuration$vcf_fingerprint$file_name, "YRI_genotype.vcf.gz") &&
      identical(cache$configuration$vcf_fingerprint$size_bytes, 118088371) &&
      identical(cache$configuration$vcf_fingerprint$md5, "4f9eb383ce3512d867b42ab806d451a8") &&
      identical(cache$sample_ids, base_cache$sample_ids) &&
      identical(names(cache$samples), as.character(contract$seeds)),
    "The pilot-10 genotype cache has an unexpected configuration."
  )
  rows <- lapply(contract$seeds, function(seed) {
    key <- as.character(seed)
    sample <- validate_real_genotype_sample(
      cache$samples[[key]],
      expected_genes = contract$n_units,
      expected_donors = cache$configuration$n_donors,
      maf_min = cache$configuration$maf_min
    )
    if (seed %in% contract$base_seeds) {
      r3ts10_assert(
        identical(sample, base_cache$samples[[key]]),
        paste("A retained base genotype sample changed for seed", seed, ".")
      )
      source <- "retained pilot-5 sample"
    } else {
      expected_selection <- sample_one_tested_variant_per_gene(pair_source$pair_ids, seed)
      r3ts10_assert(
        identical(sample$selection, expected_selection),
        paste("An added genotype selection is not reproducible for seed", seed, ".")
      )
      source <- "new VCF extraction"
    }
    expected_digest <- r3ts10g_expected_content_md5()[[as.character(seed)]]
    actual_digest <- genotype_content_md5(
      sample$selection$pair_key, rownames(sample$G), sample$G
    )
    r3ts10_assert(
      identical(actual_digest, expected_digest),
      paste("The genotype digest is invalid for seed", seed, ".")
    )
    data.frame(
      seed = seed,
      source = source,
      genes = ncol(sample$G),
      donors = nrow(sample$G),
      unique_variant_ids = length(unique(sample$selection$variant_id)),
      repeated_cross_gene_assignments = nrow(sample$selection) -
        length(unique(sample$selection$variant_id)),
      maf_min = min(sample$variant_info$observed_maf),
      maf_median = stats::median(sample$variant_info$observed_maf),
      maf_max = max(sample$variant_info$observed_maf),
      genotype_content_md5 = genotype_content_md5(
        sample$selection$pair_key, rownames(sample$G), sample$G
      ),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

r3ts10g_build_cache <- function(base_cache, pair_source, vcf_path,
                                 pair_source_sha256, base_cache_sha256,
                                 work_dir = tempdir()) {
  r3ts10g_require_dependencies()
  r3ts10g_validate_base_cache(base_cache)
  r3ts10g_validate_pair_source(pair_source, base_cache)
  contract <- r3ts10_contract()
  actual_vcf <- artifact_fingerprint(vcf_path)
  expected_vcf <- base_cache$configuration$vcf_fingerprint
  r3ts10_assert(
    identical(actual_vcf$file_name, expected_vcf$file_name) &&
      identical(actual_vcf$size_bytes, expected_vcf$size_bytes) &&
      identical(actual_vcf$md5, expected_vcf$md5),
    "The server VCF does not match the frozen real-genotype source."
  )
  r3ts10_assert(
    identical(read_vcf_sample_ids(vcf_path), base_cache$sample_ids),
    "The server VCF donor order changed."
  )
  selections <- lapply(contract$added_seeds, function(seed) {
    sample_one_tested_variant_per_gene(pair_source$pair_ids, seed)
  })
  names(selections) <- as.character(contract$added_seeds)
  target_variants <- sort(unique(unlist(
    lapply(selections, `[[`, "variant_id"),
    use.names = FALSE
  )), method = "radix")
  extracted <- extract_target_vcf_dosages(
    vcf_path = vcf_path,
    target_variants = target_variants,
    work_dir = work_dir
  )
  r3ts10_assert(
    identical(extracted$matched_variant_count, length(target_variants)),
    "The VCF did not contain every variant selected for the added seeds."
  )
  added_samples <- lapply(contract$added_seeds, function(seed) {
    sample <- assemble_one_per_gene_genotype_sample(
      selection = selections[[as.character(seed)]],
      extracted = extracted,
      seed = seed
    )
    sample <- validate_real_genotype_sample(
      sample,
      expected_genes = contract$n_units,
      expected_donors = base_cache$configuration$n_donors,
      maf_min = base_cache$configuration$maf_min
    )
    sample$genotype_digest <- object_md5(list(
      pair_key = sample$selection$pair_key,
      sample_ids = rownames(sample$G),
      G = sample$G
    ))
    sample
  })
  names(added_samples) <- as.character(contract$added_seeds)
  samples <- c(base_cache$samples, added_samples)
  samples <- samples[as.character(contract$seeds)]
  cache <- list(
    configuration = r3ts10g_configuration(
      base_cache = base_cache,
      pair_source_sha256 = pair_source_sha256,
      base_cache_sha256 = base_cache_sha256,
      vcf_fingerprint = actual_vcf
    ),
    sample_ids = base_cache$sample_ids,
    samples = samples
  )
  r3ts10g_validate_cache(
    cache, base_cache, pair_source,
    pair_source_sha256 = pair_source_sha256,
    base_cache_sha256 = base_cache_sha256
  )
  cache
}
