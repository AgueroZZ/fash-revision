# Recover one missing illustration using the original seed and frozen input generator.
# No FASH fits, direct tests, permutations, or discovery decisions are recomputed.
experiment <- "code/revision_simulations/internal/r3_stage1_method_comparison"
source(file.path(experiment, "comparison_helpers.R"))
base <- r3c_base(".")
context <- r3c_context(".", base)
seed <- 12345L
mechanism <- "raised_cosine"
directory <- "output/revision_simulations/internal/r3_six_examples_20260918"
path <- file.path(directory, "two_peak_example.rds")
if (file.exists(path)) stop("The completed illustration cache already exists.")
genotype <- readRDS(base$genotype_path)$samples[[as.character(seed)]]
reference <- base$records[[paste(mechanism, seed, sep = ":")]]
inputs <- context$sim$r3ts_reconstruct(base$config$formal_configuration,
                                      genotype, seed, mechanism)
effects <- inputs$effects
r3c_assert(identical(effects$unit_info$variant_id, reference$selected_pair_keys) &&
  identical(effects$unit_info$effect_class == "dynamic_bspline",
            reference$inference$true_dynamic), "Reconstructed unit order or truth changed.")
functional_error <- max(abs(effects$true_functionals - reference$inference$true_functionals))
r3c_assert(is.finite(functional_error) && functional_error <= 1e-12,
           "Reconstructed functionals differ from retained truth.")
formal_directory <- file.path("output/revision_simulations/mc", base$config$contract$formal_id)
r3c_manifest(formal_directory, base$config$formal_manifest_sha256)
original <- readRDS(file.path(formal_directory, "replicates", "raised_cosine_seed_12345.rds"))
examples <- unlist(original$example_curves, recursive = FALSE, use.names = FALSE)
errors <- vapply(examples, function(example) {
  i <- match(example$variant_id, reference$selected_pair_keys)
  r3c_assert(identical(effects$unit_info$spike_count[i], example$spike_count),
             "Reconstructed generating peak count differs from retained examples.")
  c(truth = max(abs(effects$beta_evaluation[i, ] - example$true_curve$true_effect)),
    estimate = max(abs(inputs$regression$beta_hat[i, ] - example$observed$estimate)),
    se = max(abs(inputs$regression$se[i, ] - example$observed$se)))
}, numeric(3))
r3c_assert(length(examples) == 12L && all(is.finite(errors)) && max(errors) <= 1e-12,
           "Reconstructed truth, estimates, or standard errors differ from retained examples.")
index <- which(reference$inference$true_dynamic & effects$unit_info$spike_count == 2L)[1L]
r3c_assert(!is.na(index), "No two-peak unit was generated.")
example <- list(mechanism = mechanism, variant_id = reference$selected_pair_keys[index],
  spike_count = 2L, true_functionals = effects$true_functionals[index, ],
  true_curve = data.frame(time = base$config$formal_configuration$evaluation_grid,
                          true_effect = effects$beta_evaluation[index, ]),
  observed = data.frame(time = as.numeric(base$config$formal_configuration$time_grid),
    estimate = inputs$regression$beta_hat[index, ], se = inputs$regression$se[index, ],
    true_effect = effects$beta_matrix[index, ]))
cache <- list(schema = "r3-two-peak-illustration-v1", seed = seed, unit_index = index,
  baseline_manifest_sha256 = r3c_contract()$base_manifest_sha256,
  scientific_source_sha256 = context$hashes,
  genotype_sha256 = base$config$genotype_cache_sha256,
  original_input_digests = reference$input_digests,
  reconstructed_input_digests = inputs$input_digests,
  max_functional_error = functional_error, example_validation_errors = errors,
  example = example, producer_sha256 = r3c_sha(file.path(experiment, "prepare_two_peak_example.R")))
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
saveRDS(cache, path, version = 3)
writeLines(r3c_sha(path), file.path(directory, "two_peak_example.sha256"))
print(list(pair_key = example$variant_id, unit_index = index,
           max_functional_error = functional_error, max_example_errors = apply(errors, 1, max)))
