# Recover original generating peak labels without expression generation or fitting.
# Run once from the workflowr root; page builds only read the resulting cache.
experiment <- "code/revision_simulations/internal/r3_stage1_method_comparison"
source(file.path(experiment, "comparison_helpers.R"))
source(file.path(experiment, "reporting.R"))
base <- r3c_base(".")
context <- r3c_context(".", base)
sim <- context$sim
config <- base$config$formal_configuration
genotypes <- readRDS(base$genotype_path)
directory <- "output/revision_simulations/internal/r3_peak_stratified_power_20260918"
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
path <- file.path(directory, "peak_labels.rds")
if (file.exists(path)) stop("Completed peak labels already exist; validate them instead of overwriting.")
records <- list()
formal_directory <- file.path("output/revision_simulations/mc", base$config$contract$formal_id)
r3c_manifest(formal_directory, base$config$formal_manifest_sha256)
for (seed in base$config$contract$seeds) {
  genotype <- genotypes$samples[[as.character(seed)]]
  reference <- base$records[[paste("raised_cosine", seed, sep = ":")]]
  cs <- sim$revision_component_seeds(seed)
  # Call the same frozen generator and MAF reassignment used by r3ts_reconstruct.
  effects <- sim$simulate_matched_functional_effect_set(
    n_variants = config$J, truth_mechanism = "raised_cosine",
    time_grid = config$time_grid, evaluation_grid = config$evaluation_grid,
    class_probs = config$class_probs, dynamic_main_effect_sd = config$dynamic_main_effect_sd,
    cosine_center_ranges = config$raised_cosine$center_ranges,
    switch_threshold = config$switch_threshold,
    location_truth_margin = config$location_truth_margin,
    location_truth_min_range_fraction = config$location_truth_min_range_fraction,
    switch_truth_margin = config$switch_truth_margin,
    non_switch_min_abs = config$non_switch_min_abs,
    non_switch_min_range_fraction = config$non_switch_min_range_fraction,
    temporal_category_probs = config$temporal_category_probs,
    seed = seed, class_seed = cs[["classes"]], constant_seed = cs[["constant_effects"]],
    shape_seed = cs[["functional_truth"]],
    scenario = paste0("r3b_real_genotype_one_per_gene_matched_functional_raised_cosine_",
      "open_middle_3_12_center_aligned_equal_cells_relative_location_clearance_",
      "full_universe_paired_posterior_main_effect"),
    middle_window = config$middle_window, middle_boundary = config$middle_boundary
  )
  effects <- sim$reassign_effect_simulation_by_maf(
    effect_sim = effects, maf = genotype$variant_info$observed_maf,
    class_probs = config$class_probs, seed = cs[["classes"]], n_strata = 10L
  )
  for (field in c("beta_matrix", "beta_evaluation", "true_functionals")) {
    rownames(effects[[field]]) <- genotype$selection$pair_key
  }
  dynamic <- effects$unit_info$effect_class == "dynamic_bspline"
  digest <- sim$serialized_object_md5(effects$beta_matrix)
  functional_error <- max(abs(effects$true_functionals - reference$inference$true_functionals))
  r3c_assert(identical(genotype$selection$pair_key, reference$selected_pair_keys) &&
    identical(dynamic, reference$inference$true_dynamic) &&
    identical(dimnames(effects$true_functionals), dimnames(reference$inference$true_functionals)) &&
    is.finite(functional_error) && functional_error <= 1e-12,
    "Recovered truth differs from completed R3 beyond numerical tolerance.")
  peaks <- effects$unit_info$spike_count
  r3c_assert(all(peaks[dynamic] %in% 1:3) && all(is.na(peaks[!dynamic])),
             "Generating peak count is invalid or assigned to null units.")
  # Cross-platform floating-point tails can change serialized matrix hashes.
  # Check all unit-level functionals and every retained full curve where available.
  example_error <- NA_real_
  original_path <- file.path(formal_directory, "replicates", paste0("raised_cosine_seed_", seed, ".rds"))
  if (file.exists(original_path)) {
    original <- readRDS(original_path)
    examples <- unlist(original$example_curves, recursive = FALSE, use.names = FALSE)
    errors <- vapply(examples, function(example) {
      index <- match(example$variant_id, reference$selected_pair_keys)
      r3c_assert(!is.na(index) && identical(peaks[index], example$spike_count),
                 "Recovered peak count differs from a retained truth example.")
      max(abs(effects$beta_evaluation[index, ] - example$true_curve$true_effect),
          abs(effects$beta_matrix[index, ] - example$observed$true_effect))
    }, numeric(1))
    example_error <- max(errors)
    r3c_assert(length(examples) == 12L && example_error <= 1e-12,
               "Recovered curves differ from the retained dense or observed-time truth.")
  }
  records[[as.character(seed)]] <- list(seed = seed,
    true_beta_md5 = unname(reference$input_digests[["true_beta_md5"]]),
    reconstructed_beta_md5 = digest, max_functional_error = functional_error,
    max_example_curve_error = example_error, true_functionals = effects$true_functionals,
    units = data.frame(unit_index = seq_len(config$J), pair_key = reference$selected_pair_keys,
      true_dynamic = dynamic, peak_count = peaks))
  message("Verified original truth and peak labels for seed ", seed)
}
cache <- list(schema = "r3-generating-peak-labels-v1", records = records,
  baseline_manifest_sha256 = r3c_contract()$base_manifest_sha256,
  genotype_sha256 = base$config$genotype_cache_sha256,
  scientific_source_sha256 = context$hashes,
  producer_sha256 = r3c_sha(file.path(experiment, "prepare_peak_stratification.R")))
saveRDS(cache, path, version = 3)
writeLines(r3c_sha(path), file.path(directory, "peak_labels.sha256"))
message("Completed verified peak-label cache: ", path)
