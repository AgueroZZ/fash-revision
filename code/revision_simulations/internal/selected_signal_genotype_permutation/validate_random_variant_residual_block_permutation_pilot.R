#!/usr/bin/env Rscript

# Independently validate a random-variant residual-block permutation cache.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(args, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(args[equals_hit[1L]], nchar(equals_prefix) + 1L))
  }
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1L] == length(args)) {
    return(default)
  }
  args[hit[1L] + 1L]
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "covariance_estimation",
  "donor_null_permutation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation",
  "selected_signal_genotype_permutation_helpers.R"
))

expected_permutation_method <- get_arg(
  "--expected-permutation-method",
  "residual_block"
)
allowed_permutation_methods <- c(
  "genotype_label_independent_time",
  "residual_block",
  "signal_stripped_residual_block",
  "signal_stripped_independent_time_residual",
  "signal_stripped_unit_specific_residual_block"
)
if (length(expected_permutation_method) != 1L ||
    !expected_permutation_method %in% allowed_permutation_methods) {
  stop("The expected residual permutation method is invalid.")
}
expected_target_selection <- get_arg(
  "--expected-target-selection",
  "random_all_tested"
)
allowed_target_selections <- c("random_all_tested", "random_all_genes")
if (length(expected_target_selection) != 1L ||
    !expected_target_selection %in% allowed_target_selections) {
  stop("The expected target-selection method is invalid.")
}
expected_target_count <- if (identical(
  expected_target_selection,
  "random_all_genes"
)) 6362L else 1177L
expected_candidate_pool <- if (identical(
  expected_target_selection,
  "random_all_genes"
)) {
  "all_tested_variants_within_all_genes"
} else {
  "all_tested_variants_within_discovered_gene"
}
default_output_id <- paste0(
  "selected_signal_random_all_tested_residual_block_permutation_",
  "selection20260817_seed20260811"
)
output_id <- get_arg("--output-id", default_output_id)
if (length(output_id) != 1L || !nzchar(output_id) ||
    grepl("/", output_id, fixed = TRUE)) {
  stop("The requested output ID is invalid.")
}

output_parent <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal"
)
output_directory <- file.path(
  output_parent,
  output_id
)
reference_donor_map_path <- file.path(
  output_parent,
  "selected_signal_residual_block_permutation_seed20260811",
  "donor_permutation.csv"
)
full_bf_path <- file.path(
  workflowr_root,
  "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
)
fit_path <- file.path(output_directory, "merged_fash_fit.rds")
effect_path <- file.path(output_directory, "effect_estimates.rds")
selection_path <- file.path(output_directory, "selection.csv")
unit_lfdr_path <- file.path(output_directory, "unit_lfdr.csv")
diagnostic_path <- file.path(output_directory, "calibration_diagnostics.csv")
donor_map_path <- file.path(output_directory, "donor_permutation.csv")
source_path <- file.path(output_directory, "source_information.csv")
required_paths <- c(
  fit_path,
  effect_path,
  selection_path,
  unit_lfdr_path,
  diagnostic_path,
  donor_map_path,
  source_path,
  full_bf_path
)
generated_donor_map_method <- expected_permutation_method %in% c(
  "genotype_label_independent_time",
  "signal_stripped_independent_time_residual",
  "signal_stripped_unit_specific_residual_block"
)
if (!generated_donor_map_method) {
  required_paths <- c(required_paths, reference_donor_map_path)
}
if (any(!file.exists(required_paths))) {
  stop("The random-variant residual-block pilot cache is incomplete.")
}

fit_bundle <- readRDS(fit_path)
effect_estimates <- readRDS(effect_path)
selection <- utils::read.csv(selection_path, stringsAsFactors = FALSE)
unit_lfdr <- utils::read.csv(unit_lfdr_path, stringsAsFactors = FALSE)
saved_diagnostics <- utils::read.csv(
  diagnostic_path,
  stringsAsFactors = FALSE
)
donor_map <- utils::read.csv(donor_map_path, stringsAsFactors = FALSE)
donor_map$target_donor <- as.character(donor_map$target_donor)
donor_map$source_donor <- as.character(donor_map$source_donor)
if ("time_value" %in% names(donor_map)) {
  donor_map$time_value <- as.numeric(donor_map$time_value)
}
reference_donor_map <- if (!identical(
  generated_donor_map_method,
  TRUE
)) {
  utils::read.csv(reference_donor_map_path, stringsAsFactors = FALSE)
} else {
  NULL
}
if (!is.null(reference_donor_map)) {
  reference_donor_map$target_donor <- as.character(
    reference_donor_map$target_donor
  )
  reference_donor_map$source_donor <- as.character(
    reference_donor_map$source_donor
  )
}
source_information <- utils::read.csv(source_path, stringsAsFactors = FALSE)

configuration <- fit_bundle$configuration
expected_sample_counts <- c(
  19L, 19L, 16L, 19L, 16L, 19L, 19L, 19L,
  19L, 19L, 19L, 19L, 19L, 18L, 19L, 19L
)
if (!identical(
      configuration$permutation_method,
      expected_permutation_method
    ) ||
    !identical(
      configuration$target_selection_method,
      expected_target_selection
    ) ||
    configuration$seed != 20260811L ||
    configuration$selection_seed != 20260817L ||
    configuration$alpha != 0.05 ||
    configuration$n_pair_level_discoveries != 9205L ||
    configuration$n_target_units != expected_target_count ||
    configuration$n_permuted_null_units != expected_target_count ||
    !identical(
      configuration$selection_candidate_pool,
      expected_candidate_pool
    ) ||
    !identical(as.integer(configuration$sample_counts), expected_sample_counts) ||
    !identical(
      as.integer(configuration$residual_df),
      expected_sample_counts - 7L
    )) {
  stop("The random-variant residual-block configuration is invalid.")
}
if (expected_permutation_method %in% c(
  "signal_stripped_residual_block",
  "signal_stripped_independent_time_residual",
  "signal_stripped_unit_specific_residual_block"
)) {
  saved_correlation <- as.numeric(
    configuration$source_residual_genotype_correlation_by_time
  )
  temporal_alignment_valid <- if (is.null(
    configuration$temporal_residual_alignment_preserved
  )) {
    !identical(
      expected_permutation_method,
      "signal_stripped_independent_time_residual"
    )
  } else {
    identical(
      configuration$temporal_residual_alignment_preserved,
      !identical(
        expected_permutation_method,
        "signal_stripped_independent_time_residual"
      )
    )
  }
  expected_cross_unit_alignment <- !identical(
    expected_permutation_method,
    "signal_stripped_unit_specific_residual_block"
  )
  cross_unit_alignment_valid <- if (is.null(
    configuration$cross_unit_residual_alignment_preserved
  )) {
    expected_cross_unit_alignment
  } else {
    identical(
      configuration$cross_unit_residual_alignment_preserved,
      expected_cross_unit_alignment
    )
  }
  if (!identical(
        configuration$residual_source,
        "full_model_Y_on_G_and_PCs"
      ) || !isTRUE(configuration$fitted_genotype_removed) ||
      !identical(configuration$leverage_adjustment, "none") ||
      length(saved_correlation) != 16L ||
      any(!is.finite(saved_correlation)) ||
      any(saved_correlation < 0 | saved_correlation > 1e-10) ||
      !isTRUE(all.equal(
        max(saved_correlation),
        configuration$maximum_source_residual_genotype_correlation,
        tolerance = 1e-15
      )) || !temporal_alignment_valid || !cross_unit_alignment_valid) {
    stop("The signal-stripped residual provenance is invalid.")
  }
}
if (identical(
  expected_permutation_method,
  "genotype_label_independent_time"
)) {
  saved_correlation <- as.numeric(
    configuration$source_residual_genotype_correlation_by_time
  )
  if (!is.na(configuration$residual_source) ||
      isTRUE(configuration$fitted_genotype_removed) ||
      !identical(configuration$temporal_genotype_alignment_preserved, FALSE) ||
      !identical(
        configuration$cross_unit_genotype_alignment_preserved,
        TRUE
      ) || length(saved_correlation) != 16L ||
      any(!is.na(saved_correlation)) ||
      !is.na(configuration$maximum_source_residual_genotype_correlation)) {
    stop("The independent-time genotype provenance is invalid.")
  }
}
if (identical(expected_target_selection, "random_all_genes")) {
  source_roles <- source_information$role
  required_clean_roles <- c(
    "prepared_null_input_cache",
    "precomputed_null_fit_cache"
  )
  if (!identical(configuration$likelihood_source,
                 "precomputed_clean_session") ||
      !all(required_clean_roles %in% source_roles) ||
      !identical(
        configuration$prepared_null_input_md5,
        source_information$md5[
          match("prepared_null_input_cache", source_roles)
        ]
      ) ||
      !identical(
        configuration$precomputed_null_fit_md5,
        source_information$md5[
          match("precomputed_null_fit_cache", source_roles)
        ]
      )) {
    stop("The clean-session likelihood provenance is invalid.")
  }
}

required_selection_columns <- c(
  "target_index", "source_fash_index", "pair_key", "gene_id", "variant_id",
  "candidate_variant_count", "source_bf_lfdr",
  "source_pair_level_discovery"
)
if (!identical(names(selection), required_selection_columns) ||
    nrow(selection) != expected_target_count ||
    anyDuplicated(selection$pair_key) ||
    anyDuplicated(selection$gene_id) ||
    any(selection$candidate_variant_count < 1L) ||
    any(!is.finite(selection$source_bf_lfdr)) ||
    any(selection$source_bf_lfdr < 0 | selection$source_bf_lfdr > 1)) {
  stop("The saved random target selection is invalid.")
}

message("Reconstructing the fixed-seed selection from the full BF fit.")
bf_environment <- new.env(parent = emptyenv())
loaded_names <- load(full_bf_path, envir = bf_environment)
if (!identical(loaded_names, "fash_fit1_update")) {
  stop("The BF-adjusted source file has unexpected objects.")
}
full_bf <- bf_environment$fash_fit1_update
pair_keys <- names(full_bf$fash_data$data_list)
expected_selection <- if (identical(
  expected_target_selection,
  "random_all_genes"
)) {
  select_random_tested_pair_per_gene(
    pair_keys = pair_keys,
    bf_lfdr = full_bf$lfdr,
    alpha = 0.05,
    selection_seed = 20260817L
  )
} else {
  select_random_tested_pair_per_discovered_gene(
    pair_keys = pair_keys,
    bf_lfdr = full_bf$lfdr,
    alpha = 0.05,
    selection_seed = 20260817L
  )
}
selection_character_columns <- c(
  "pair_key", "gene_id", "variant_id", "source_pair_level_discovery"
)
selection_numeric_columns <- c(
  "target_index", "source_fash_index", "candidate_variant_count",
  "source_bf_lfdr"
)
for (column in selection_character_columns) {
  if (!identical(selection[[column]], expected_selection[[column]])) {
    stop("The saved selection disagrees with reconstruction: ", column)
  }
}
for (column in selection_numeric_columns) {
  if (!isTRUE(all.equal(
    selection[[column]],
    expected_selection[[column]],
    tolerance = 1e-14,
    check.attributes = FALSE
  ))) {
    stop("The saved selection disagrees with reconstruction: ", column)
  }
}
if (!isTRUE(all.equal(
  fit_bundle$selection,
  expected_selection,
  tolerance = 1e-14,
  check.attributes = FALSE
))) {
  stop("The saved fit-bundle selection disagrees with reconstruction.")
}
rm(full_bf, bf_environment, pair_keys, expected_selection)
gc(verbose = FALSE)

if (generated_donor_map_method) {
  observation_patterns <- as.character(
    configuration$donor_observation_patterns
  )
  donor_ids <- names(configuration$donor_observation_patterns)
}
if (expected_permutation_method %in% c(
  "genotype_label_independent_time",
  "signal_stripped_independent_time_residual"
)) {
  observation_matrix <- t(vapply(observation_patterns, function(pattern) {
    tokens <- regmatches(pattern, gregexpr("TRUE|FALSE", pattern))[[1L]]
    if (length(tokens) != length(configuration$time_grid)) {
      stop("A saved donor observation pattern has the wrong length.")
    }
    tokens == "TRUE"
  }, logical(length(configuration$time_grid))))
  rownames(observation_matrix) <- donor_ids
  expected_donor_map <- make_independent_time_donor_permutations(
    donor_observation_matrix = observation_matrix,
    time_grid = configuration$time_grid,
    seed = configuration$seed
  )
  reused_source_row <- source_information$role == "reused_residual_donor_map"
  valid_time_maps <- vapply(seq_along(configuration$time_grid), function(j) {
    current <- donor_map[donor_map$time_index == j, , drop = FALSE]
    observed <- donor_ids[observation_matrix[, j]]
    nrow(current) == configuration$sample_counts[j] &&
      !anyDuplicated(current$target_donor) &&
      !anyDuplicated(current$source_donor) &&
      setequal(current$target_donor, observed) &&
      setequal(current$source_donor, observed) &&
      any(!current$fixed_point)
  }, logical(1))
  if (!identical(donor_map, expected_donor_map) ||
      !all(valid_time_maps) || sum(reused_source_row) != 0L ||
      !identical(configuration$donor_map_source_path, "") ||
      !is.na(configuration$donor_map_source_md5)) {
    stop("The independent time-specific donor maps failed validation.")
  }
} else if (identical(
  expected_permutation_method,
  "signal_stripped_unit_specific_residual_block"
)) {
  expected_donor_map <- make_unit_specific_donor_block_permutations(
    donor_ids = donor_ids,
    observation_patterns = setNames(observation_patterns, donor_ids),
    unit_keys = selection$pair_key,
    seed = configuration$seed
  )
  reused_source_row <- source_information$role == "reused_residual_donor_map"
  unit_counts <- table(donor_map$unit_index)
  nonidentity_by_unit <- tapply(
    !donor_map$fixed_point,
    donor_map$unit_index,
    any
  )
  if (!identical(donor_map, expected_donor_map) ||
      length(unit_counts) != expected_target_count ||
      any(unit_counts != length(donor_ids)) ||
      length(nonidentity_by_unit) != expected_target_count ||
      !all(nonidentity_by_unit) ||
      anyDuplicated(paste(donor_map$unit_index, donor_map$target_donor)) ||
      anyDuplicated(paste(donor_map$unit_index, donor_map$source_donor)) ||
      sum(reused_source_row) != 0L ||
      !identical(configuration$donor_map_source_path, "") ||
      !is.na(configuration$donor_map_source_md5)) {
    stop("The unit-specific donor-block maps failed validation.")
  }
} else {
  reference_map_md5 <- unname(tools::md5sum(reference_donor_map_path))
  saved_map_md5 <- unname(tools::md5sum(donor_map_path))
  reused_source_row <- source_information$role == "reused_residual_donor_map"
  if (!identical(donor_map, reference_donor_map) ||
      saved_map_md5 != reference_map_md5 ||
      !identical(
        normalizePath(configuration$donor_map_source_path, mustWork = TRUE),
        normalizePath(reference_donor_map_path, mustWork = TRUE)
      ) ||
      configuration$donor_map_source_md5 != reference_map_md5 ||
      sum(reused_source_row) != 1L ||
      source_information$md5[reused_source_row] != reference_map_md5) {
    stop("The saved residual donor map is not the exact reused reference map.")
  }
  source_pattern <- donor_map$observation_pattern[
    match(donor_map$source_donor, donor_map$target_donor)
  ]
  if (nrow(donor_map) != 19L ||
      anyDuplicated(donor_map$target_donor) ||
      anyDuplicated(donor_map$source_donor) ||
      !setequal(donor_map$target_donor, donor_map$source_donor) ||
      all(donor_map$fixed_point) ||
      !identical(source_pattern, donor_map$observation_pattern)) {
    stop("The reused residual donor map failed its structural invariants.")
  }
}

expected_stages <- c("Raw", "BF-adjusted")
if (!identical(unique(unit_lfdr$fit_stage), expected_stages) ||
    nrow(unit_lfdr) != 4L * expected_target_count ||
    anyDuplicated(unit_lfdr[, c("fit_stage", "unit_key")]) ||
    !setequal(unique(unit_lfdr$group), c("target", "permuted_null")) ||
    any(!is.finite(unit_lfdr$lfdr)) ||
    any(unit_lfdr$lfdr < 0 | unit_lfdr$lfdr > 1)) {
  stop("The saved unit-level lfdr table is invalid.")
}

fits <- list(
  Raw = fit_bundle$raw_fit,
  `BF-adjusted` = fit_bundle$bf_adjusted_fit
)
diagnostics <- do.call(rbind, lapply(expected_stages, function(stage) {
  rows <- unit_lfdr$fit_stage == stage
  fit <- fits[[stage]]
  saved_lfdr <- unit_lfdr$lfdr[rows]
  if (max(abs(saved_lfdr - as.numeric(fit$lfdr))) > 1e-14) {
    stop("Saved and fitted lfdr values disagree for stage: ", stage)
  }
  summarize_matched_null_calibration(
    lfdr = saved_lfdr,
    group = unit_lfdr$group[rows],
    pi0_merged = extract_pi0(fit),
    fit_stage = stage,
    alpha = configuration$alpha
  )
}))
rownames(diagnostics) <- NULL
numeric_diagnostic_columns <- names(diagnostics)[vapply(
  diagnostics,
  is.numeric,
  logical(1)
)]
for (column in numeric_diagnostic_columns) {
  saved_diagnostics[[column]] <- as.numeric(saved_diagnostics[[column]])
}
if (!isTRUE(all.equal(
  diagnostics,
  saved_diagnostics,
  tolerance = 1e-14,
  check.attributes = FALSE
))) {
  stop("The independently recomputed calibration diagnostics disagree.")
}
if (any(abs(
  diagnostics$known_null_discovery_fraction -
    diagnostics$permuted_null_calls / diagnostics$total_calls
) > 1e-15, na.rm = TRUE) ||
    any(abs(
      diagnostics$scaled_fdr_merged_from_design_lower_bound -
        diagnostics$known_null_discovery_fraction
    ) > 1e-15, na.rm = TRUE) ||
    any(
      diagnostics$pi0_target_valid !=
        is.finite(diagnostics$post_selection_fdr_target_from_pi0)
    )) {
  stop("The recomputed calibration identities failed validation.")
}

effect_matrices <- effect_estimates[c(
  "target_beta",
  "target_adjusted_se",
  "observed_refit_beta",
  "observed_refit_raw_se",
  "observed_refit_adjusted_se",
  "permuted_beta",
  "permuted_raw_se",
  "permuted_adjusted_se"
)]
if (any(vapply(effect_matrices, function(x) {
  !identical(dim(x), c(expected_target_count, 16L)) ||
    any(!is.finite(x))
}, logical(1))) ||
    any(effect_estimates$target_adjusted_se <= 0) ||
    any(effect_estimates$observed_refit_raw_se <= 0) ||
    any(effect_estimates$observed_refit_adjusted_se <= 0) ||
    any(effect_estimates$permuted_raw_se <= 0) ||
    any(effect_estimates$permuted_adjusted_se <= 0)) {
  stop("The saved effect-estimate matrices are invalid.")
}

if (identical(
  expected_permutation_method,
  "genotype_label_independent_time"
)) {
  message("Independently reconstructing time-specific genotype permutations.")
  source_paths <- setNames(
    source_information$path,
    source_information$role
  )
  required_source_roles <- c(
    "genotype_vcf",
    "expression_matrix",
    "time_specific_pc_data"
  )
  if (!all(required_source_roles %in% names(source_paths))) {
    stop("The genotype-permutation reconstruction sources are incomplete.")
  }
  expression_data <- utils::read.csv(
    source_paths[["expression_matrix"]],
    sep = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  pc_data <- utils::read.delim(
    source_paths[["time_specific_pc_data"]],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expression_sample_ids <- names(expression_data)[-1L]
  gene_rows <- match(selection$gene_id, expression_data$Gene_id)
  if (anyNA(gene_rows) ||
      !setequal(expression_sample_ids, pc_data$Sample_id)) {
    stop("The genotype-permutation reconstruction inputs are misaligned.")
  }
  unique_variant_ids <- unique(selection$variant_id)
  unique_dosage <- read_selected_vcf_dosages(
    source_paths[["genotype_vcf"]],
    unique_variant_ids,
    chunk_size = 100000L
  )
  unit_dosage <- unique_dosage[
    ,
    match(selection$variant_id, unique_variant_ids),
    drop = FALSE
  ]
  colnames(unit_dosage) <- selection$pair_key
  reconstructed_observed_beta <- matrix(
    NA_real_,
    nrow = nrow(selection),
    ncol = 16L,
    dimnames = dimnames(effect_estimates$observed_refit_beta)
  )
  reconstructed_observed_raw_se <- reconstructed_observed_beta
  reconstructed_null_beta <- reconstructed_observed_beta
  reconstructed_null_raw_se <- reconstructed_observed_beta
  reconstructed_residual_df <- integer(16L)
  for (time_index in seq_along(configuration$time_grid)) {
    time_value <- configuration$time_grid[time_index]
    sample_ids <- grep(
      paste0("_", time_value, "$"),
      expression_sample_ids,
      value = TRUE
    )
    donors <- sub(paste0("_", time_value, "$"), "", sample_ids)
    expression_columns <- match(sample_ids, names(expression_data))
    pc_rows <- match(sample_ids, pc_data$Sample_id)
    expression_matrix <- t(as.matrix(
      expression_data[gene_rows, expression_columns, drop = FALSE]
    ))
    storage.mode(expression_matrix) <- "double"
    rownames(expression_matrix) <- donors
    colnames(expression_matrix) <- selection$pair_key
    covariates <- as.matrix(
      pc_data[pc_rows, paste0("PC", 1:5), drop = FALSE]
    )
    storage.mode(covariates) <- "double"
    rownames(covariates) <- donors
    current_donor_map <- donor_map[
      donor_map$time_index == time_index,
      ,
      drop = FALSE
    ]
    permuted_time_dosage <- apply_donor_map_to_genotype(
      genotype = unit_dosage,
      donor_map = current_donor_map,
      target_donors = donors
    )
    observed_fit <- fit_many_genotype_regressions(
      expression_matrix,
      unit_dosage[donors, , drop = FALSE],
      covariates
    )
    null_fit <- fit_many_genotype_regressions(
      expression_matrix,
      permuted_time_dosage,
      covariates
    )
    reconstructed_observed_beta[, time_index] <- observed_fit$beta
    reconstructed_observed_raw_se[, time_index] <-
      observed_fit$standard_error
    reconstructed_null_beta[, time_index] <- null_fit$beta
    reconstructed_null_raw_se[, time_index] <- null_fit$standard_error
    reconstructed_residual_df[time_index] <- observed_fit$residual_df
    if (observed_fit$residual_df != null_fit$residual_df) {
      stop("The reconstructed residual df differ at time ", time_value, ".")
    }
  }
  reconstructed_observed_adjusted_se <-
    convert_raw_to_original_t_adjusted_se(
      reconstructed_observed_beta,
      reconstructed_observed_raw_se,
      reconstructed_residual_df
    )
  reconstructed_null_adjusted_se <- convert_raw_to_original_t_adjusted_se(
    reconstructed_null_beta,
    reconstructed_null_raw_se,
    reconstructed_residual_df
  )
  reconstruction_differences <- c(
    observed_beta = max(abs(
      reconstructed_observed_beta - effect_estimates$observed_refit_beta
    )),
    observed_raw_se = max(abs(
      reconstructed_observed_raw_se - effect_estimates$observed_refit_raw_se
    )),
    observed_adjusted_se = max(abs(
      reconstructed_observed_adjusted_se -
        effect_estimates$observed_refit_adjusted_se
    )),
    null_beta = max(abs(
      reconstructed_null_beta - effect_estimates$permuted_beta
    )),
    null_raw_se = max(abs(
      reconstructed_null_raw_se - effect_estimates$permuted_raw_se
    )),
    null_adjusted_se = max(abs(
      reconstructed_null_adjusted_se - effect_estimates$permuted_adjusted_se
    ))
  )
  if (any(!is.finite(reconstruction_differences)) ||
      any(reconstruction_differences > 1e-12) ||
      !identical(
        reconstructed_residual_df,
        as.integer(effect_estimates$residual_df)
      )) {
    stop(
      "The independently reconstructed genotype effects disagree: ",
      paste(
        names(reconstruction_differences),
        format(reconstruction_differences, digits = 6),
        collapse = ", "
      )
    )
  }
  rm(
    expression_data,
    pc_data,
    unique_dosage,
    unit_dosage,
    reconstructed_observed_beta,
    reconstructed_observed_raw_se,
    reconstructed_null_beta,
    reconstructed_null_raw_se
  )
  gc(verbose = FALSE)
}

if (expected_permutation_method %in% c(
  "signal_stripped_residual_block",
  "signal_stripped_independent_time_residual",
  "signal_stripped_unit_specific_residual_block"
)) {
  message("Independently reconstructing signal-stripped null effect estimates.")
  source_paths <- setNames(
    source_information$path,
    source_information$role
  )
  required_source_roles <- c(
    "genotype_vcf",
    "expression_matrix",
    "time_specific_pc_data"
  )
  if (!all(required_source_roles %in% names(source_paths))) {
    stop("The signal-stripped reconstruction sources are incomplete.")
  }
  expression_data <- utils::read.csv(
    source_paths[["expression_matrix"]],
    sep = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  pc_data <- utils::read.delim(
    source_paths[["time_specific_pc_data"]],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  expression_sample_ids <- names(expression_data)[-1L]
  gene_rows <- match(selection$gene_id, expression_data$Gene_id)
  if (anyNA(gene_rows) ||
      !setequal(expression_sample_ids, pc_data$Sample_id)) {
    stop("The signal-stripped reconstruction inputs are misaligned.")
  }
  unique_variant_ids <- unique(selection$variant_id)
  unique_dosage <- read_selected_vcf_dosages(
    source_paths[["genotype_vcf"]],
    unique_variant_ids,
    chunk_size = 100000L
  )
  unit_dosage <- unique_dosage[
    ,
    match(selection$variant_id, unique_variant_ids),
    drop = FALSE
  ]
  colnames(unit_dosage) <- selection$pair_key
  unit_source_donor_matrix <- if (identical(
    expected_permutation_method,
    "signal_stripped_unit_specific_residual_block"
  )) {
    matrix(
      donor_map$source_donor,
      nrow = nrow(unit_dosage),
      ncol = nrow(selection),
      dimnames = list(rownames(unit_dosage), selection$pair_key)
    )
  } else {
    NULL
  }
  reconstructed_observed_beta <- matrix(
    NA_real_,
    nrow = nrow(selection),
    ncol = 16L,
    dimnames = dimnames(effect_estimates$observed_refit_beta)
  )
  reconstructed_observed_raw_se <- reconstructed_observed_beta
  reconstructed_null_beta <- reconstructed_observed_beta
  reconstructed_null_raw_se <- reconstructed_observed_beta
  reconstructed_correlation <- rep(NA_real_, 16L)
  reconstructed_residual_df <- integer(16L)
  for (time_index in seq_along(configuration$time_grid)) {
    time_value <- configuration$time_grid[time_index]
    sample_ids <- grep(
      paste0("_", time_value, "$"),
      expression_sample_ids,
      value = TRUE
    )
    donors <- sub(paste0("_", time_value, "$"), "", sample_ids)
    expression_columns <- match(sample_ids, names(expression_data))
    pc_rows <- match(sample_ids, pc_data$Sample_id)
    expression_matrix <- t(as.matrix(
      expression_data[gene_rows, expression_columns, drop = FALSE]
    ))
    storage.mode(expression_matrix) <- "double"
    rownames(expression_matrix) <- donors
    colnames(expression_matrix) <- selection$pair_key
    covariates <- as.matrix(
      pc_data[pc_rows, paste0("PC", 1:5), drop = FALSE]
    )
    storage.mode(covariates) <- "double"
    rownames(covariates) <- donors
    if (identical(
      expected_permutation_method,
      "signal_stripped_unit_specific_residual_block"
    )) {
      source_donors <- unit_source_donor_matrix[
        match(donors, rownames(unit_source_donor_matrix)),
        ,
        drop = FALSE
      ]
      source_rows <- matrix(
        match(as.vector(source_donors), donors),
        nrow = length(donors),
        ncol = nrow(selection)
      )
    } else {
      current_donor_map <- if (identical(
        expected_permutation_method,
        "signal_stripped_independent_time_residual"
      )) {
        donor_map[donor_map$time_index == time_index, , drop = FALSE]
      } else {
        donor_map
      }
      source_donors <- current_donor_map$source_donor[
        match(donors, current_donor_map$target_donor)
      ]
      source_rows <- match(source_donors, donors)
    }
    if (anyNA(source_rows)) {
      stop("The donor map cannot reconstruct time ", time_value, ".")
    }
    reconstruction <- make_signal_stripped_residual_block_null(
      expression = expression_matrix,
      genotype = unit_dosage[donors, , drop = FALSE],
      covariates = covariates,
      source_rows = source_rows
    )
    reconstructed_observed_beta[, time_index] <-
      reconstruction$observed_fit$beta
    reconstructed_observed_raw_se[, time_index] <-
      reconstruction$observed_fit$standard_error
    reconstructed_null_beta[, time_index] <- reconstruction$null_fit$beta
    reconstructed_null_raw_se[, time_index] <-
      reconstruction$null_fit$standard_error
    reconstructed_correlation[time_index] <-
      reconstruction$maximum_residual_genotype_correlation
    reconstructed_residual_df[time_index] <-
      reconstruction$observed_fit$residual_df
  }
  reconstructed_null_adjusted_se <- convert_raw_to_original_t_adjusted_se(
    reconstructed_null_beta,
    reconstructed_null_raw_se,
    reconstructed_residual_df
  )
  reconstruction_differences <- c(
    observed_beta = max(abs(
      reconstructed_observed_beta - effect_estimates$observed_refit_beta
    )),
    observed_raw_se = max(abs(
      reconstructed_observed_raw_se - effect_estimates$observed_refit_raw_se
    )),
    null_beta = max(abs(
      reconstructed_null_beta - effect_estimates$permuted_beta
    )),
    null_raw_se = max(abs(
      reconstructed_null_raw_se - effect_estimates$permuted_raw_se
    )),
    null_adjusted_se = max(abs(
      reconstructed_null_adjusted_se - effect_estimates$permuted_adjusted_se
    )),
    residual_genotype_correlation = max(abs(
      reconstructed_correlation -
        effect_estimates$source_residual_genotype_correlation_by_time
    ))
  )
  if (any(!is.finite(reconstruction_differences)) ||
      any(reconstruction_differences > 1e-12) ||
      !identical(
        reconstructed_residual_df,
        as.integer(effect_estimates$residual_df)
      )) {
    stop(
      "The independently reconstructed signal-stripped effects disagree: ",
      paste(
        names(reconstruction_differences),
        format(reconstruction_differences, digits = 6),
        collapse = ", "
      )
    )
  }
}

for (stage in expected_stages) {
  weights <- fits[[stage]]$prior_weights$prior_weight
  if (any(!is.finite(weights)) || any(weights < 0) ||
      abs(sum(weights) - 1) > 1e-8) {
    stop("The fitted mixture weights are invalid for stage: ", stage)
  }
}

if (any(!file.exists(source_information$path)) ||
    !identical(
      unname(tools::md5sum(source_information$path)),
      source_information$md5
    )) {
  stop("At least one immutable source file no longer matches its saved hash.")
}

cat(
  "Random-variant ",
  expected_target_selection,
  " ",
  expected_permutation_method,
  " permutation pilot validation passed.\n",
  sep = ""
)
print(diagnostics, row.names = FALSE)
