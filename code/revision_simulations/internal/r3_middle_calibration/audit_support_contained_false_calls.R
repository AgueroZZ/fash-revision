#!/usr/bin/env Rscript

# Reconstruct the immutable support-contained R3 truth and characterize the
# geometry of alpha-0.05 Middle calls. This is a read-only diagnostic: it does
# not refit FASH or modify the retained formal cache.

options(stringsAsFactors = FALSE)

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/real_genotype_one_per_gene.R")) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }
  if (file.exists(file.path(
    "coderepo-local", "code", "revision_simulations", "shared",
    "real_genotype_one_per_gene.R"
  ))) {
    return(normalizePath(
      "coderepo-local", winslash = "/", mustWork = TRUE
    ))
  }
  stop("Could not find the workflowr repository root.", call. = FALSE)
}

expect_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

collapse_numeric <- function(x) {
  if (length(x) == 0L) return("")
  paste(formatC(x, digits = 8L, format = "fg"), collapse = ";")
}

workflowr_root <- find_workflowr_root()
simulation_functions_path <- file.path(
  workflowr_root, "code", "revision_simulations", "r3_r4_fashr0143",
  "source_snapshots", "r3_prior_geometry_simulation_functions.R"
)
genotype_helper_path <- file.path(
  workflowr_root, "code", "revision_simulations", "shared",
  "real_genotype_one_per_gene.R"
)
source(simulation_functions_path)
source(genotype_helper_path)

cache_id <- paste0(
  "r3_real_genotype_one_per_gene_J6362_matched_functional_",
  "open_middle_3_12_support_contained_relative_clearance_",
  "main_effect_fashr0143_pilot5"
)
cache_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "mc", cache_id
)
configuration <- readRDS(file.path(cache_dir, "configuration.rds"))
genotype_cache <- readRDS(file.path(
  workflowr_root, "output", "revision_simulations", "shared",
  "real_genotype_one_per_gene_J6362_pilot5", "genotype_samples.rds"
))

expect_true(
  identical(configuration$raised_cosine$center_ranges, list(
    early = c(0.5, 1.5), middle = c(4.5, 10.5), late = c(13.5, 14.5)
  )),
  "The retained cache is not the support-contained design."
)
expect_true(
  identical(configuration$middle_window, c(3, 12)) &&
    identical(configuration$middle_boundary, "open"),
  "The retained cache does not use the open 3 < t < 12 definition."
)

replicate_paths <- sort(list.files(
  file.path(cache_dir, "replicates"),
  pattern = "^(random_bspline|raised_cosine)_seed_[0-9]+[.]rds$",
  full.names = TRUE
))
expect_true(length(replicate_paths) == 10L, "Expected ten formal replicates.")
replicates <- lapply(replicate_paths, readRDS)

call_rows <- do.call(rbind, lapply(replicates, function(x) {
  x$call_diagnostics_alpha005
}))
call_rows <- call_rows[
  call_rows$method == "FASH-IWP1-BF" & call_rows$target == "middle",
  , drop = FALSE
]

reconstruct_truth <- function(seed, truth_mechanism) {
  component_seeds <- revision_component_seeds(seed)
  effect_sim <- simulate_matched_functional_effect_set(
    n_variants = configuration$J,
    truth_mechanism = truth_mechanism,
    time_grid = configuration$time_grid,
    evaluation_grid = configuration$evaluation_grid,
    class_probs = configuration$class_probs,
    dynamic_main_effect_sd = configuration$dynamic_main_effect_sd,
    bspline_amplitude = configuration$random_bspline$amplitude,
    bspline_df = configuration$random_bspline$df,
    bspline_coefficient_sd = configuration$random_bspline$coefficient_sd,
    cosine_width_half = configuration$raised_cosine$width_half,
    cosine_spike_counts = configuration$raised_cosine$spike_counts,
    cosine_relative_amplitude_range =
      configuration$raised_cosine$relative_amplitude_range,
    cosine_target_centered_rms =
      configuration$raised_cosine$target_centered_rms,
    switch_threshold = configuration$switch_threshold,
    location_truth_margin = configuration$location_truth_margin,
    switch_truth_margin = configuration$switch_truth_margin,
    non_switch_min_abs = configuration$non_switch_min_abs,
    non_switch_min_range_fraction =
      configuration$non_switch_min_range_fraction,
    temporal_category_probs = NULL,
    seed = seed,
    class_seed = component_seeds[["classes"]],
    constant_seed = component_seeds[["constant_effects"]],
    shape_seed = component_seeds[["functional_truth"]],
    middle_window = configuration$middle_window,
    middle_boundary = configuration$middle_boundary
  )
  genotype_sample <- genotype_cache$samples[[as.character(seed)]]
  effect_sim <- reassign_effect_simulation_by_maf(
    effect_sim = effect_sim,
    maf = genotype_sample$variant_info$observed_maf,
    class_probs = configuration$class_probs,
    seed = component_seeds[["classes"]],
    n_strata = 10L
  )
  effect_sim$unit_info$observed_maf <-
    genotype_sample$variant_info$observed_maf
  effect_sim
}

requested <- unique(call_rows[c("seed", "truth_mechanism")])
geometry_rows <- vector("list", nrow(requested))
curve_rows <- vector("list", nrow(requested))
population_rows <- vector("list", nrow(requested))
for (i in seq_len(nrow(requested))) {
  seed <- requested$seed[[i]]
  mechanism <- requested$truth_mechanism[[i]]
  effect_sim <- reconstruct_truth(seed, mechanism)
  selected_calls <- call_rows[
    call_rows$seed == seed & call_rows$truth_mechanism == mechanism,
    , drop = FALSE
  ]
  indices <- selected_calls$unit_index
  unit_info <- effect_sim$unit_info[indices, , drop = FALSE]
  reconstructed_functional <- effect_sim$true_functionals[indices, "middle"]
  expect_true(
    max(abs(reconstructed_functional - selected_calls$true_functional)) < 1e-12,
    paste("Truth reconstruction failed for", mechanism, "seed", seed)
  )
  primary_center <- vapply(
    unit_info$peak_centers,
    function(x) if (length(x) == 0L) NA_real_ else x[[1L]],
    numeric(1)
  )
  all_centers <- vapply(unit_info$peak_centers, collapse_numeric, character(1))
  support_left <- primary_center - configuration$raised_cosine$width_half
  support_right <- primary_center + configuration$raised_cosine$width_half
  domain_truncation <- pmax(0, configuration$time_grid[[1L]] - support_left) +
    pmax(0, support_right - tail(configuration$time_grid, 1L))
  geometry_rows[[i]] <- data.frame(
    selected_calls,
    time_group = unit_info$time_group,
    switch_status = unit_info$switch_status,
    spike_count = unit_info$spike_count,
    primary_center = primary_center,
    all_centers = all_centers,
    primary_support_left = support_left,
    primary_support_right = support_right,
    primary_domain_truncation = domain_truncation,
    generation_attempt = unit_info$generation_attempt,
    centered_rms = unit_info$centered_rms,
    minimum_absolute_effect = unit_info$minimum_absolute_effect,
    effect_range = unit_info$effect_range,
    relative_clearance = unit_info$relative_clearance,
    observed_maf = unit_info$observed_maf,
    stringsAsFactors = FALSE
  )
  curve_rows[[i]] <- do.call(rbind, lapply(seq_along(indices), function(k) {
    data.frame(
      seed = seed,
      truth_mechanism = mechanism,
      unit_index = indices[[k]],
      false_discovery = selected_calls$false_discovery[[k]],
      truth_group = selected_calls$truth_group[[k]],
      time = configuration$evaluation_grid,
      true_effect = effect_sim$beta_evaluation[indices[[k]], ],
      stringsAsFactors = FALSE
    )
  }))
  dynamic_indices <- which(
    effect_sim$unit_info$effect_class == "dynamic_bspline"
  )
  dynamic_info <- effect_sim$unit_info[dynamic_indices, , drop = FALSE]
  dynamic_primary_center <- vapply(
    dynamic_info$peak_centers,
    function(x) if (length(x) == 0L) NA_real_ else x[[1L]],
    numeric(1)
  )
  dynamic_support_left <- dynamic_primary_center -
    configuration$raised_cosine$width_half
  dynamic_support_right <- dynamic_primary_center +
    configuration$raised_cosine$width_half
  dynamic_domain_truncation <- pmax(
    0, configuration$time_grid[[1L]] - dynamic_support_left
  ) + pmax(
    0, dynamic_support_right - tail(configuration$time_grid, 1L)
  )
  false_indices <- selected_calls$unit_index[selected_calls$false_discovery]
  population_rows[[i]] <- data.frame(
    seed = seed,
    truth_mechanism = mechanism,
    unit_index = dynamic_indices,
    time_group = dynamic_info$time_group,
    switch_status = dynamic_info$switch_status,
    spike_count = dynamic_info$spike_count,
    primary_center = dynamic_primary_center,
    primary_domain_truncation = dynamic_domain_truncation,
    true_middle_functional =
      effect_sim$true_functionals[dynamic_indices, "middle"],
    false_middle_call_alpha005 = dynamic_indices %in% false_indices,
    stringsAsFactors = FALSE
  )
}

geometry <- do.call(rbind, geometry_rows)
curves <- do.call(rbind, curve_rows)
population <- do.call(rbind, population_rows)
rownames(geometry) <- NULL
rownames(curves) <- NULL
rownames(population) <- NULL

diagnostic_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "diagnostics",
  "r3_middle_calibration"
)
dir.create(diagnostic_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  geometry,
  file.path(diagnostic_dir, "support_contained_call_geometry_alpha005.csv"),
  row.names = FALSE
)
utils::write.csv(
  curves,
  file.path(diagnostic_dir, "support_contained_call_truth_curves_alpha005.csv"),
  row.names = FALSE
)
utils::write.csv(
  population,
  file.path(
    diagnostic_dir, "support_contained_dynamic_truth_population.csv"
  ),
  row.names = FALSE
)

raised <- geometry[geometry$truth_mechanism == "raised_cosine", , drop = FALSE]
raised$call_status <- ifelse(
  raised$false_discovery, "false early/late call", "true middle call"
)
summary_groups <- split(
  raised,
  interaction(raised$call_status, raised$spike_count, drop = TRUE)
)
summary_table <- do.call(rbind, lapply(summary_groups, function(rows) {
  data.frame(
    call_status = rows$call_status[[1L]],
    spike_count = rows$spike_count[[1L]],
    calls = nrow(rows),
    mean_true_functional = mean(rows$true_functional),
    median_true_functional = stats::median(rows$true_functional),
    mean_primary_center = mean(rows$primary_center),
    mean_domain_truncation = mean(rows$primary_domain_truncation),
    proportion_domain_truncated = mean(rows$primary_domain_truncation > 0),
    mean_abs_main_effect = mean(abs(rows$genetic_main_effect)),
    mean_observed_maf = mean(rows$observed_maf),
    stringsAsFactors = FALSE
  )
}))
summary_table <- summary_table[
  order(summary_table$call_status, summary_table$spike_count),
  , drop = FALSE
]
utils::write.csv(
  summary_table,
  file.path(diagnostic_dir, "support_contained_call_geometry_summary.csv"),
  row.names = FALSE
)

raised_population <- population[
  population$truth_mechanism == "raised_cosine" &
    population$time_group %in% c("early", "late"),
  , drop = FALSE
]
raised_population$truncation_bin <- cut(
  raised_population$primary_domain_truncation,
  breaks = c(-Inf, 0.25, 0.50, 0.75, Inf),
  labels = c("[0,0.25)", "[0.25,0.50)", "[0.50,0.75)", "[0.75,1.00]"),
  right = FALSE
)
population_groups <- split(
  raised_population,
  interaction(
    raised_population$spike_count,
    raised_population$truncation_bin,
    drop = TRUE
  )
)
population_summary <- do.call(rbind, lapply(population_groups, function(rows) {
  data.frame(
    spike_count = rows$spike_count[[1L]],
    truncation_bin = as.character(rows$truncation_bin[[1L]]),
    dynamic_truths = nrow(rows),
    false_middle_calls = sum(rows$false_middle_call_alpha005),
    false_middle_call_rate = mean(rows$false_middle_call_alpha005),
    mean_truncation = mean(rows$primary_domain_truncation),
    stringsAsFactors = FALSE
  )
}))
population_summary <- population_summary[
  order(population_summary$spike_count, population_summary$mean_truncation),
  , drop = FALSE
]
utils::write.csv(
  population_summary,
  file.path(
    diagnostic_dir, "support_contained_boundary_false_call_rates.csv"
  ),
  row.names = FALSE
)

cat("Reconstructed and validated all formal alpha-0.05 BF Middle calls.\n")
print(summary_table, row.names = FALSE)
cat("\nRaised-cosine early/late false-call rates by boundary truncation:\n")
print(population_summary, row.names = FALSE)
cat("Diagnostic artifacts:", diagnostic_dir, "\n")
