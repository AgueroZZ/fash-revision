#!/usr/bin/env Rscript

# Run one selected-signal matched-null FASH permutation pilot.

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

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

file_metadata <- function(paths, roles) {
  paths <- normalizePath(paths, mustWork = TRUE)
  information <- file.info(paths)
  data.frame(
    role = roles,
    path = paths,
    size_bytes = unname(information$size),
    modification_time = format(information$mtime, tz = "UTC", usetz = TRUE),
    md5 = unname(tools::md5sum(paths)),
    stringsAsFactors = FALSE
  )
}

capture_warnings <- function(expression) {
  warning_messages <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(warning) {
      warning_messages <<- c(warning_messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warning_messages))
}

extract_effect_matrix <- function(data_list, field, time_grid) {
  matrix_result <- t(vapply(data_list, function(unit_data) {
    if (!identical(as.numeric(unit_data$x), as.numeric(time_grid))) {
      stop("A source FASH unit does not use the expected time grid.")
    }
    as.numeric(unit_data[[field]])
  }, numeric(length(time_grid))))
  matrix_result
}

workflowr_root <- find_workflowr_root()
project_root <- normalizePath(file.path(workflowr_root, ".."), mustWork = TRUE)
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

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required.")
}

seed <- as.integer(get_arg("--seed", "20260811"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
num_cores <- as.integer(get_arg("--num-cores", "8"))
permutation_method <- get_arg("--permutation-method", "genotype_label")
allowed_permutation_methods <- c(
  "genotype_label",
  "genotype_label_independent_time",
  "residual_block",
  "signal_stripped_residual_block",
  "signal_stripped_independent_time_residual",
  "signal_stripped_unit_specific_residual_block"
)
target_selection_method <- get_arg("--target-selection", "top_discovered")
allowed_target_selection_methods <- c(
  "top_discovered",
  "random_all_tested",
  "random_all_genes"
)
selection_seed <- as.integer(get_arg("--selection-seed", "20260817"))
donor_map_path_argument <- get_arg("--donor-map-path", "")
donor_map_source_path <- if (nzchar(donor_map_path_argument)) {
  normalizePath(donor_map_path_argument, mustWork = TRUE)
} else {
  ""
}
write_null_input_argument <- get_arg("--write-null-input-cache", "")
read_null_input_argument <- get_arg("--read-null-input-cache", "")
precomputed_null_fit_argument <- get_arg(
  "--precomputed-null-fit-cache",
  ""
)
normalize_output_path <- function(path) {
  file.path(
    normalizePath(dirname(path), mustWork = TRUE),
    basename(path)
  )
}
write_null_input_path <- if (nzchar(write_null_input_argument)) {
  normalize_output_path(write_null_input_argument)
} else {
  ""
}
read_null_input_path <- if (nzchar(read_null_input_argument)) {
  normalizePath(read_null_input_argument, mustWork = TRUE)
} else {
  ""
}
precomputed_null_fit_path <- if (nzchar(precomputed_null_fit_argument)) {
  normalizePath(precomputed_null_fit_argument, mustWork = TRUE)
} else {
  ""
}
default_output_id <- if (identical(
  target_selection_method,
  "random_all_tested"
)) {
  paste0(
    "selected_signal_random_all_tested_",
    permutation_method,
    "_permutation_selection",
    selection_seed,
    "_seed",
    seed
  )
} else if (identical(target_selection_method, "random_all_genes")) {
  paste0(
    "all_gene_random_variant_",
    permutation_method,
    "_permutation_selection",
    selection_seed,
    "_seed",
    seed
  )
} else if (identical(permutation_method, "genotype_label")) {
  paste0("selected_signal_genotype_permutation_seed", seed)
} else if (identical(permutation_method, "residual_block")) {
  paste0("selected_signal_residual_block_permutation_seed", seed)
} else {
  paste0(
    "selected_signal_signal_stripped_residual_block_permutation_seed",
    seed
  )
}
output_id <- get_arg(
  "--output-id",
  default_output_id
)
if (length(seed) != 1L || is.na(seed) || length(alpha) != 1L ||
    !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
    length(num_cores) != 1L || is.na(num_cores) || num_cores < 1L ||
    length(permutation_method) != 1L ||
    !permutation_method %in% allowed_permutation_methods ||
    length(target_selection_method) != 1L ||
    !target_selection_method %in% allowed_target_selection_methods ||
    length(selection_seed) != 1L || is.na(selection_seed) ||
    (nzchar(write_null_input_path) && file.exists(write_null_input_path)) ||
    (nzchar(write_null_input_path) &&
      (nzchar(read_null_input_path) || nzchar(precomputed_null_fit_path))) ||
    (nzchar(precomputed_null_fit_path) && !nzchar(read_null_input_path)) ||
    (permutation_method %in% c(
      "genotype_label_independent_time",
      "signal_stripped_independent_time_residual",
      "signal_stripped_unit_specific_residual_block"
    ) && nzchar(donor_map_source_path)) ||
    !nzchar(output_id) || grepl("/", output_id, fixed = TRUE)) {
  stop(paste(
    "Invalid seed, alpha, core count, permutation method, target selection,",
    "selection seed, donor-map path, or output ID."
  ))
}

raw_fit_path <- file.path(
  workflowr_root,
  "output", "dynamic_eQTL_real", "fash_fit1_all.RData"
)
bf_fit_path <- file.path(
  workflowr_root,
  "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
)
top_pair_bundle_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "one_variant_per_gene_best_lfdr", "thinned_fash_fit.rds"
)
top_pair_summary_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "one_variant_per_gene_best_lfdr", "gene_discovery_summary.csv"
)
vcf_path <- file.path(
  project_root,
  "iPSC-data", "genotype-data", "YRI_genotype.vcf.gz"
)
expression_path <- file.path(
  project_root,
  "iPSC-data", "expression-data",
  "quantile_normalized_no_projection.txt"
)
pc_path <- file.path(
  workflowr_root,
  "data", "dynamic_eQTL_real", "principal_components_10.txt"
)
required_paths <- c(
  raw_fit_path,
  bf_fit_path,
  top_pair_bundle_path,
  top_pair_summary_path,
  vcf_path,
  expression_path,
  pc_path
)
required_roles <- c(
  "full_raw_fash",
  "full_bf_adjusted_fash",
  "validated_top_pair_bundle",
  "validated_top_pair_summary",
  "genotype_vcf",
  "expression_matrix",
  "time_specific_pc_data"
)
if (nzchar(donor_map_source_path)) {
  required_paths <- c(required_paths, donor_map_source_path)
  required_roles <- c(required_roles, "reused_residual_donor_map")
}
if (nzchar(read_null_input_path)) {
  required_paths <- c(required_paths, read_null_input_path)
  required_roles <- c(required_roles, "prepared_null_input_cache")
}
if (nzchar(precomputed_null_fit_path)) {
  required_paths <- c(required_paths, precomputed_null_fit_path)
  required_roles <- c(required_roles, "precomputed_null_fit_cache")
}
if (any(!file.exists(required_paths))) {
  stop(
    "At least one required input is missing: ",
    paste(required_paths[!file.exists(required_paths)], collapse = ", ")
  )
}

output_parent <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal"
)
final_directory <- file.path(output_parent, output_id)
staging_directory <- file.path(
  output_parent,
  paste0(".", output_id, "_staging_", Sys.getpid())
)
if (file.exists(final_directory) || file.exists(staging_directory)) {
  stop("Refusing to overwrite an existing output or staging directory.")
}
dir.create(staging_directory, recursive = FALSE)
staging_complete <- FALSE
on.exit({
  if (!staging_complete && dir.exists(staging_directory)) {
    unlink(staging_directory, recursive = TRUE)
  }
}, add = TRUE)

analysis_start <- proc.time()[["elapsed"]]
message("Recording source provenance and immutable-input hashes.")
source_information_before <- file_metadata(required_paths, required_roles)

message("Loading the full BF-adjusted fit and selecting target units.")
bf_environment <- new.env(parent = emptyenv())
loaded_names <- load(bf_fit_path, envir = bf_environment)
if (!identical(loaded_names, "fash_fit1_update")) {
  stop("The BF-adjusted source file must contain only fash_fit1_update.")
}
full_bf <- bf_environment$fash_fit1_update
pair_keys <- names(full_bf$fash_data$data_list)
if (length(pair_keys) != 1009173L || length(full_bf$lfdr) != length(pair_keys) ||
    !identical(names(full_bf$lfdr), pair_keys) || anyDuplicated(pair_keys) ||
    length(full_bf$psd_grid) != 52L ||
    !identical(full_bf$settings$order, 1) ||
    !identical(full_bf$settings$pred_step, 1) ||
    !identical(full_bf$settings$penalty, 10)) {
  stop("The full BF-adjusted source fit does not match the primary analysis.")
}
selection <- if (identical(target_selection_method, "top_discovered")) {
  select_discovered_gene_top_pairs(
    pair_keys = pair_keys,
    bf_lfdr = full_bf$lfdr,
    alpha = alpha
  )
} else if (identical(target_selection_method, "random_all_tested")) {
  select_random_tested_pair_per_discovered_gene(
    pair_keys = pair_keys,
    bf_lfdr = full_bf$lfdr,
    alpha = alpha,
    selection_seed = selection_seed
  )
} else {
  select_random_tested_pair_per_gene(
    pair_keys = pair_keys,
    bf_lfdr = full_bf$lfdr,
    alpha = alpha,
    selection_seed = selection_seed
  )
}
n_pair_discoveries <- attr(selection, "n_pair_level_discoveries")
original_psd_grid <- full_bf$psd_grid
original_settings <- full_bf$settings
if (target_selection_method %in% c("random_all_tested", "random_all_genes")) {
  if (!identical(names(full_bf$fash_data$data_list), pair_keys) ||
      !identical(names(full_bf$fash_data$S), pair_keys) ||
      nrow(full_bf$L_matrix) != length(pair_keys) ||
      (!is.null(rownames(full_bf$L_matrix)) &&
        !identical(rownames(full_bf$L_matrix), pair_keys))) {
    stop("The full FASH target inputs are not aligned by pair key.")
  }
  message("Extracting random target likelihood rows and beta/SE datasets.")
  target_rows <- selection$source_fash_index
  target_data_list <- full_bf$fash_data$data_list[target_rows]
  target_se_list <- full_bf$fash_data$S[target_rows]
  target_likelihood <- full_bf$L_matrix[target_rows, , drop = FALSE]
  source_raw <- list(
    psd_grid = original_psd_grid,
    settings = original_settings
  )
}

top_pair_summary <- utils::read.csv(
  top_pair_summary_path,
  stringsAsFactors = FALSE
)
expected_target_count <- if (identical(
  target_selection_method,
  "random_all_genes"
)) 6362L else 1177L
if (nrow(top_pair_summary) != 1L ||
    top_pair_summary$full_discovered_pairs != n_pair_discoveries ||
    top_pair_summary$full_discovered_genes != 1177L ||
    n_pair_discoveries != 9205L ||
    nrow(selection) != expected_target_count) {
  stop("The primary discovery counts do not match the validated reference.")
}

if (identical(target_selection_method, "top_discovered")) {
  message("Loading validated target likelihood rows and beta/SE datasets.")
  top_pair_bundle <- readRDS(top_pair_bundle_path)
  source_raw <- top_pair_bundle$raw_fit
  bundle_selection <- top_pair_bundle$selection
  if (!identical(
    names(source_raw$fash_data$data_list),
    bundle_selection$pair_key
  ) || !identical(rownames(source_raw$L_matrix), bundle_selection$pair_key) ||
      !isTRUE(all.equal(source_raw$psd_grid, original_psd_grid)) ||
      !identical(source_raw$settings, original_settings)) {
    stop("The validated top-pair bundle is not aligned with the primary fit.")
  }
  bundle_rows <- match(selection$pair_key, bundle_selection$pair_key)
  if (anyNA(bundle_rows) || anyDuplicated(bundle_rows)) {
    stop("At least one selected target is unavailable in the top-pair bundle.")
  }
  target_data_list <- source_raw$fash_data$data_list[bundle_rows]
  target_se_list <- source_raw$fash_data$S[bundle_rows]
  target_likelihood <- source_raw$L_matrix[bundle_rows, , drop = FALSE]
}
names(target_data_list) <- selection$pair_key
names(target_se_list) <- selection$pair_key
rownames(target_likelihood) <- selection$pair_key
if (nrow(target_likelihood) != nrow(selection) ||
    ncol(target_likelihood) != length(original_psd_grid) ||
    length(target_data_list) != nrow(selection) ||
    length(target_se_list) != nrow(selection)) {
  stop("The extracted target likelihood and beta/SE inputs are incomplete.")
}
rm(full_bf, bf_environment, pair_keys)
gc(verbose = FALSE)

time_grid <- 0:15
target_beta <- extract_effect_matrix(target_data_list, "y", time_grid)
target_adjusted_se <- t(vapply(
  target_se_list,
  as.numeric,
  numeric(length(time_grid))
))
dimnames(target_beta) <- list(selection$pair_key, paste0("time_", time_grid))
dimnames(target_adjusted_se) <- dimnames(target_beta)
if (any(!is.finite(target_beta)) || any(!is.finite(target_adjusted_se)) ||
    any(target_adjusted_se <= 0)) {
  stop("The selected target beta/SE matrices are invalid.")
}

message("Reading selected expression rows and time-specific PCs.")
expression_data <- utils::read.csv(
  expression_path,
  sep = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
pc_data <- utils::read.delim(
  pc_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!identical(names(expression_data)[1L], "Gene_id") ||
    !all(c("Sample_id", paste0("PC", 1:5)) %in% names(pc_data)) ||
    anyDuplicated(expression_data$Gene_id) || anyDuplicated(pc_data$Sample_id)) {
  stop("The expression or PC input has invalid identifiers.")
}
expression_sample_ids <- names(expression_data)[-1L]
if (length(expression_sample_ids) != 297L || nrow(pc_data) != 297L ||
    !setequal(expression_sample_ids, pc_data$Sample_id)) {
  stop("Expression and PC donor-time IDs do not match exactly.")
}
gene_rows <- match(selection$gene_id, expression_data$Gene_id)
if (anyNA(gene_rows)) {
  stop("At least one selected target gene is missing from expression data.")
}

message("Streaming selected variants from the VCF.")
unique_variant_ids <- unique(selection$variant_id)
unique_dosage <- read_selected_vcf_dosages(
  vcf_path,
  unique_variant_ids,
  chunk_size = 100000L
)
unit_dosage <- unique_dosage[
  ,
  match(selection$variant_id, unique_variant_ids),
  drop = FALSE
]
colnames(unit_dosage) <- selection$pair_key
vcf_donors <- rownames(unit_dosage)
if (nrow(unit_dosage) != 19L || anyDuplicated(vcf_donors)) {
  stop("The selected dosage matrix does not contain 19 unique donors.")
}

donor_observation_matrix <- vapply(time_grid, function(time_value) {
  paste0(vcf_donors, "_", time_value) %in% expression_sample_ids
}, logical(length(vcf_donors)))
rownames(donor_observation_matrix) <- vcf_donors
colnames(donor_observation_matrix) <- paste0("time_", time_grid)
donor_observation_patterns <- apply(
  donor_observation_matrix,
  1L,
  paste0,
  collapse = ""
)

unit_source_donor_matrix <- NULL
if (identical(permutation_method, "genotype_label")) {
  message("Applying one shared donor permutation to the genotype matrix.")
  shared_permutation <- make_shared_genotype_permutation(unit_dosage, seed)
  permuted_dosage <- shared_permutation$genotype
  donor_map <- shared_permutation$donor_map
  donor_map$observation_pattern <- unname(
    donor_observation_patterns[donor_map$target_donor]
  )
} else if (permutation_method %in% c(
  "genotype_label_independent_time",
  "signal_stripped_independent_time_residual"
)) {
  message(
    "Applying independent observed-donor permutations at each time point ",
    "for method ", permutation_method, "."
  )
  donor_map <- make_independent_time_donor_permutations(
    donor_observation_matrix = donor_observation_matrix,
    time_grid = time_grid,
    seed = seed
  )
  permuted_dosage <- unit_dosage
} else if (identical(
  permutation_method,
  "signal_stripped_unit_specific_residual_block"
)) {
  message(
    "Applying one independent donor-trajectory permutation per target unit."
  )
  donor_map <- make_unit_specific_donor_block_permutations(
    donor_ids = vcf_donors,
    observation_patterns = donor_observation_patterns,
    unit_keys = selection$pair_key,
    seed = seed
  )
  unit_source_donor_matrix <- matrix(
    donor_map$source_donor,
    nrow = length(vcf_donors),
    ncol = nrow(selection),
    dimnames = list(vcf_donors, selection$pair_key)
  )
  permuted_dosage <- unit_dosage
} else {
  if (nzchar(donor_map_source_path)) {
    message("Reusing the supplied synchronized residual donor-block map.")
    donor_map <- utils::read.csv(
      donor_map_source_path,
      stringsAsFactors = FALSE
    )
    required_donor_map_columns <- c(
      "target_donor", "source_donor", "fixed_point", "observation_pattern"
    )
    if (!all(required_donor_map_columns %in% names(donor_map))) {
      stop("The supplied residual donor map is missing required columns.")
    }
    donor_map <- donor_map[, required_donor_map_columns]
    donor_map$target_donor <- as.character(donor_map$target_donor)
    donor_map$source_donor <- as.character(donor_map$source_donor)
    donor_map$fixed_point <- as.logical(donor_map$fixed_point)
    donor_map$observation_pattern <- as.character(
      donor_map$observation_pattern
    )
    donor_rows <- match(vcf_donors, donor_map$target_donor)
    if (anyNA(donor_rows)) {
      stop("The supplied residual donor map does not cover every VCF donor.")
    }
    donor_map <- donor_map[donor_rows, , drop = FALSE]
    rownames(donor_map) <- NULL
  } else {
    message(
      "Applying one synchronized residual donor-block permutation within ",
      "missingness strata."
    )
    set.seed(seed)
    source_donor <- make_shared_donor_block_permutation(
      vcf_donors,
      donor_observation_patterns
    )
    donor_map <- data.frame(
      target_donor = vcf_donors,
      source_donor = unname(source_donor[vcf_donors]),
      fixed_point = vcf_donors == unname(source_donor[vcf_donors]),
      observation_pattern = unname(donor_observation_patterns[vcf_donors]),
      stringsAsFactors = FALSE
    )
  }
  if (nrow(donor_map) != length(vcf_donors) ||
      anyDuplicated(donor_map$target_donor) ||
      anyDuplicated(donor_map$source_donor) ||
      !setequal(donor_map$target_donor, vcf_donors) ||
      !setequal(donor_map$source_donor, vcf_donors) ||
      any(donor_map$fixed_point !=
            (donor_map$target_donor == donor_map$source_donor))) {
    stop("The residual donor map is not a valid donor bijection.")
  }
  permuted_dosage <- unit_dosage
  source_patterns <- unname(
    donor_observation_patterns[donor_map$source_donor]
  )
  target_patterns <- unname(
    donor_observation_patterns[donor_map$target_donor]
  )
  if (!identical(source_patterns, target_patterns) ||
      !identical(target_patterns, donor_map$observation_pattern)) {
    stop("The residual donor map crossed missingness-pattern strata.")
  }
}
if (all(donor_map$fixed_point)) {
  stop("The fixed-seed donor permutation is the identity permutation.")
}

expected_sample_counts <- c(
  19L, 19L, 16L, 19L, 16L, 19L, 19L, 19L,
  19L, 19L, 19L, 19L, 19L, 18L, 19L, 19L
)
permuted_beta <- matrix(
  NA_real_,
  nrow = nrow(selection),
  ncol = length(time_grid),
  dimnames = dimnames(target_beta)
)
permuted_raw_se <- permuted_beta
observed_refit_beta <- permuted_beta
observed_refit_raw_se <- permuted_beta
residual_df <- integer(length(time_grid))
sample_counts <- integer(length(time_grid))
source_residual_genotype_correlation <- rep(
  NA_real_,
  length(time_grid)
)
genotype_crossproduct_difference_by_time <- rep(
  NA_real_,
  length(time_grid)
)
genotype_sorted_dosage_difference_by_time <- rep(
  NA_real_,
  length(time_grid)
)

message("Refitting the original time-specific PC model for observed and null data.")
for (time_index in seq_along(time_grid)) {
  time_value <- time_grid[time_index]
  sample_ids <- grep(
    paste0("_", time_value, "$"),
    expression_sample_ids,
    value = TRUE
  )
  donors <- sub(paste0("_", time_value, "$"), "", sample_ids)
  if (length(sample_ids) != expected_sample_counts[time_index] ||
      anyDuplicated(donors) || any(!donors %in% vcf_donors)) {
    stop("Unexpected donor coverage at time ", time_value, ".")
  }
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

  source_rows <- NULL
  if (!identical(permutation_method, "genotype_label")) {
    if (identical(
      permutation_method,
      "signal_stripped_unit_specific_residual_block"
    )) {
      source_donors <- unit_source_donor_matrix[
        match(donors, vcf_donors),
        ,
        drop = FALSE
      ]
      source_rows <- matrix(
        match(as.vector(source_donors), donors),
        nrow = length(donors),
        ncol = nrow(selection)
      )
    } else {
      current_donor_map <- if (permutation_method %in% c(
        "genotype_label_independent_time",
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
      stop("A residual donor block is unavailable at time ", time_value, ".")
    }
  }

  if (permutation_method %in% c(
    "signal_stripped_residual_block",
    "signal_stripped_independent_time_residual",
    "signal_stripped_unit_specific_residual_block"
  )) {
    signal_stripped_null <- make_signal_stripped_residual_block_null(
      expression = expression_matrix,
      genotype = unit_dosage[donors, , drop = FALSE],
      covariates = covariates,
      source_rows = source_rows
    )
    observed_fit <- signal_stripped_null$observed_fit
    permuted_fit <- signal_stripped_null$null_fit
    source_residual_genotype_correlation[time_index] <-
      signal_stripped_null$maximum_residual_genotype_correlation
  } else {
    observed_fit <- fit_many_genotype_regressions(
      expression_matrix,
      unit_dosage[donors, , drop = FALSE],
      covariates
    )
  }
  if (identical(permutation_method, "genotype_label")) {
    permuted_fit <- fit_many_genotype_regressions(
      expression_matrix,
      permuted_dosage[donors, , drop = FALSE],
      covariates
    )
  } else if (identical(
    permutation_method,
    "genotype_label_independent_time"
  )) {
    permuted_time_dosage <- apply_donor_map_to_genotype(
      genotype = unit_dosage,
      donor_map = current_donor_map,
      target_donors = donors
    )
    permuted_fit <- fit_many_genotype_regressions(
      expression_matrix,
      permuted_time_dosage,
      covariates
    )
    genotype_crossproduct_difference_by_time[time_index] <- max(abs(
      crossprod(unit_dosage[donors, , drop = FALSE]) -
        crossprod(permuted_time_dosage)
    ))
    genotype_sorted_dosage_difference_by_time[time_index] <- max(vapply(
      seq_len(ncol(unit_dosage)),
      function(column_index) {
        max(abs(
          sort(unit_dosage[donors, column_index]) -
            sort(permuted_time_dosage[, column_index])
        ))
      },
      numeric(1)
    ))
  } else if (identical(permutation_method, "residual_block")) {
    permuted_expression_residual <- observed_fit$expression_residual[
      source_rows,
      ,
      drop = FALSE
    ]
    permuted_expression_residual <-
      observed_fit$projection$residualizer %*%
      permuted_expression_residual
    permuted_fit <- fit_residualized_genotype_regressions(
      expression_residual = permuted_expression_residual,
      genotype = unit_dosage[donors, , drop = FALSE],
      residualizer = observed_fit$projection$residualizer,
      covariate_rank = observed_fit$projection$rank
    )
  }
  observed_refit_beta[, time_index] <- observed_fit$beta
  observed_refit_raw_se[, time_index] <- observed_fit$standard_error
  permuted_beta[, time_index] <- permuted_fit$beta
  permuted_raw_se[, time_index] <- permuted_fit$standard_error
  if (observed_fit$residual_df != permuted_fit$residual_df) {
    stop("Observed and permuted residual degrees of freedom differ.")
  }
  residual_df[time_index] <- observed_fit$residual_df
  sample_counts[time_index] <- length(sample_ids)
}

observed_refit_adjusted_se <- convert_raw_to_original_t_adjusted_se(
  observed_refit_beta,
  observed_refit_raw_se,
  residual_df
)
permuted_adjusted_se <- convert_raw_to_original_t_adjusted_se(
  permuted_beta,
  permuted_raw_se,
  residual_df
)
maximum_observed_beta_difference <- max(abs(observed_refit_beta - target_beta))
maximum_observed_se_difference <- max(
  abs(observed_refit_adjusted_se - target_adjusted_se)
)
beta_maximum_index <- arrayInd(
  which.max(abs(observed_refit_beta - target_beta)),
  dim(observed_refit_beta)
)
se_maximum_index <- arrayInd(
  which.max(abs(observed_refit_adjusted_se - target_adjusted_se)),
  dim(observed_refit_adjusted_se)
)
message(
  "Observed-input reproduction maxima: beta=",
  format(maximum_observed_beta_difference, digits = 8),
  " at ", rownames(observed_refit_beta)[beta_maximum_index[1L]],
  "/", colnames(observed_refit_beta)[beta_maximum_index[2L]],
  "; adjusted SE=",
  format(maximum_observed_se_difference, digits = 8),
  " at ", rownames(observed_refit_adjusted_se)[se_maximum_index[1L]],
  "/", colnames(observed_refit_adjusted_se)[se_maximum_index[2L]],
  "."
)
if (maximum_observed_beta_difference > 1e-10 ||
    maximum_observed_se_difference > 1e-10 ||
    !identical(sample_counts, expected_sample_counts) ||
    !identical(residual_df, expected_sample_counts - 7L)) {
  stop("Raw-data alignment did not reproduce the selected source beta/SE values.")
}

null_unit_keys <- paste0(
  selection$pair_key,
  "__",
  target_selection_method,
  if (identical(target_selection_method, "random_all_tested")) {
    paste0("_selection", selection_seed)
  } else {
    ""
  },
  "_",
  permutation_method,
  "_null_seed",
  seed
)
null_datasets <- lapply(seq_len(nrow(selection)), function(unit_index) {
  data.frame(
    beta = permuted_beta[unit_index, ],
    time = time_grid,
    SE = permuted_adjusted_se[unit_index, ],
    stringsAsFactors = FALSE
  )
})
names(null_datasets) <- null_unit_keys

null_input_bundle <- list(
  cache_version = "selected_signal_null_input_v1",
  output_id = output_id,
  permutation_method = permutation_method,
  target_selection_method = target_selection_method,
  seed = seed,
  selection_seed = selection_seed,
  selected_pair_keys = selection$pair_key,
  null_unit_keys = null_unit_keys,
  null_datasets = null_datasets,
  original_settings = original_settings,
  original_psd_grid = original_psd_grid,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)

if (nzchar(write_null_input_path)) {
  temporary_null_input_path <- paste0(
    write_null_input_path,
    ".tmp_",
    Sys.getpid()
  )
  saveRDS(null_input_bundle, temporary_null_input_path)
  if (!file.rename(temporary_null_input_path, write_null_input_path)) {
    unlink(temporary_null_input_path)
    stop("Could not atomically save the prepared null-input cache.")
  }
  unlink(staging_directory, recursive = TRUE)
  staging_complete <- TRUE
  cat("Prepared null-input cache: ", write_null_input_path, "\n", sep = "")
  cat("Prepared null units: ", nrow(selection), "\n", sep = "")
  quit(save = "no", status = 0L)
}

message("Computing FASH likelihood rows for the permuted negative controls.")
null_likelihood_start <- proc.time()[["elapsed"]]
likelihood_source <- "direct_current_session"
prepared_null_input_md5 <- NA_character_
precomputed_null_fit_md5 <- NA_character_
if (nzchar(precomputed_null_fit_path)) {
  message("Loading a verified clean-session null likelihood fit.")
  prepared_bundle <- readRDS(read_null_input_path)
  precomputed_bundle <- readRDS(precomputed_null_fit_path)
  prepared_null_input_md5 <- unname(tools::md5sum(read_null_input_path))
  precomputed_null_fit_md5 <- unname(tools::md5sum(precomputed_null_fit_path))
  if (!identical(prepared_bundle$cache_version,
                 "selected_signal_null_input_v1") ||
      !identical(precomputed_bundle$cache_version,
                 "selected_signal_null_fit_v1") ||
      !identical(precomputed_bundle$source_input_md5,
                 prepared_null_input_md5) ||
      !identical(prepared_bundle$null_datasets, null_datasets) ||
      !identical(prepared_bundle$null_unit_keys, null_unit_keys) ||
      !identical(prepared_bundle$selected_pair_keys, selection$pair_key) ||
      !identical(prepared_bundle$original_settings, original_settings) ||
      !isTRUE(all.equal(
        prepared_bundle$original_psd_grid,
        original_psd_grid,
        tolerance = 0
      ))) {
    stop("The precomputed null fit does not match reconstructed null inputs.")
  }
  null_capture <- list(
    value = precomputed_bundle$fit,
    warnings = precomputed_bundle$warnings
  )
  null_likelihood_elapsed <- precomputed_bundle$elapsed_seconds
  likelihood_source <- "precomputed_clean_session"
} else {
  null_capture <- capture_warnings(fashr::fash(
    Y = "beta",
    smooth_var = "time",
    S = "SE",
    data_list = null_datasets,
    num_basis = original_settings$num_basis,
    order = original_settings$order,
    betaprec = original_settings$betaprec,
    pred_step = original_settings$pred_step,
    penalty = original_settings$penalty,
    grid = original_psd_grid,
    num_cores = num_cores,
    verbose = TRUE
  ))
  null_likelihood_elapsed <- proc.time()[["elapsed"]] -
    null_likelihood_start
}
null_fit <- null_capture$value
if (!isTRUE(all.equal(null_fit$psd_grid, original_psd_grid)) ||
    nrow(null_fit$L_matrix) != nrow(selection) ||
    ncol(null_fit$L_matrix) != length(original_psd_grid) ||
    length(null_fit$lfdr) != nrow(selection)) {
  stop("The permuted-null FASH likelihood fit is not aligned with the source fit.")
}

null_data_list <- null_fit$fash_data$data_list
null_se_list <- null_fit$fash_data$S
names(null_data_list) <- null_unit_keys
names(null_se_list) <- null_unit_keys
rownames(null_fit$L_matrix) <- null_unit_keys
unit_keys <- c(selection$pair_key, null_unit_keys)
group <- rep(c("target", "permuted_null"), each = nrow(selection))
merged_data_list <- c(target_data_list, null_data_list)
merged_se_list <- c(target_se_list, null_se_list)
merged_likelihood <- rbind(target_likelihood, null_fit$L_matrix)

message("Refitting the merged raw empirical-Bayes mixture.")
raw_start <- proc.time()[["elapsed"]]
raw_capture <- capture_warnings(refit_fash_from_likelihood(
  source_fit = source_raw,
  data_list = merged_data_list,
  se_list = merged_se_list,
  likelihood_matrix = merged_likelihood,
  unit_keys = unit_keys,
  penalty = original_settings$penalty
))
raw_elapsed <- proc.time()[["elapsed"]] - raw_start
merged_raw <- raw_capture$value

message("Applying the BF prior update to the merged fit.")
bf_start <- proc.time()[["elapsed"]]
bf_capture <- capture_warnings(fashr::BF_update(merged_raw, plot = FALSE))
bf_elapsed <- proc.time()[["elapsed"]] - bf_start
merged_bf <- bf_capture$value
names(merged_bf$lfdr) <- unit_keys

raw_pi0 <- extract_pi0(merged_raw)
bf_pi0 <- extract_pi0(merged_bf)
unit_lfdr <- rbind(
  data.frame(
    fit_stage = "Raw",
    unit_key = unit_keys,
    source_pair_key = rep(selection$pair_key, 2L),
    group = group,
    lfdr = as.numeric(merged_raw$lfdr),
    stringsAsFactors = FALSE
  ),
  data.frame(
    fit_stage = "BF-adjusted",
    unit_key = unit_keys,
    source_pair_key = rep(selection$pair_key, 2L),
    group = group,
    lfdr = as.numeric(merged_bf$lfdr),
    stringsAsFactors = FALSE
  )
)
calibration_diagnostics <- rbind(
  summarize_matched_null_calibration(
    merged_raw$lfdr,
    group,
    raw_pi0,
    fit_stage = "Raw",
    alpha = alpha
  ),
  summarize_matched_null_calibration(
    merged_bf$lfdr,
    group,
    bf_pi0,
    fit_stage = "BF-adjusted",
    alpha = alpha
  )
)
lfdr_quantiles <- rbind(
  summarize_group_lfdr(merged_raw$lfdr, group, "Raw"),
  summarize_group_lfdr(merged_bf$lfdr, group, "BF-adjusted")
)
prior_weights <- rbind(
  transform(merged_raw$prior_weights, fit_stage = "Raw"),
  transform(merged_bf$prior_weights, fit_stage = "BF-adjusted")
)
prior_weights <- prior_weights[, c("fit_stage", "psd", "prior_weight")]

if (identical(permutation_method, "genotype_label_independent_time")) {
  maximum_crossproduct_difference <- max(
    genotype_crossproduct_difference_by_time
  )
  maximum_sorted_dosage_difference <- max(
    genotype_sorted_dosage_difference_by_time
  )
} else {
  maximum_crossproduct_difference <- max(abs(
    crossprod(unit_dosage) - crossprod(permuted_dosage)
  ))
  maximum_sorted_dosage_difference <- max(vapply(
    seq_len(ncol(unit_dosage)),
    function(column_index) {
      max(abs(
        sort(unit_dosage[, column_index]) -
          sort(permuted_dosage[, column_index])
      ))
    },
    numeric(1)
  ))
}
selection_discovery_count <- sum(selection$source_pair_level_discovery)
candidate_count_summary <- if (target_selection_method %in% c(
  "random_all_tested",
  "random_all_genes"
)) {
  c(
    minimum = min(selection$candidate_variant_count),
    median = stats::median(selection$candidate_variant_count),
    mean = mean(selection$candidate_variant_count),
    maximum = max(selection$candidate_variant_count)
  )
} else {
  c(minimum = NA_real_, median = NA_real_, mean = NA_real_, maximum = NA_real_)
}
validation <- data.frame(
  check = c(
    "Full pair-level BF discoveries",
    "Discovered genes and target units",
    "Permuted-null units",
    "Unique target genes",
    "Unique target pair keys",
    "Target units that were original pair-level discoveries",
    "Minimum candidate variants per target gene",
    "Median candidate variants per target gene",
    "Mean candidate variants per target gene",
    "Maximum candidate variants per target gene",
    "VCF donors",
    "Donor permutation fixed points",
    "Maximum analysis-genotype cross-product difference",
    "Maximum sorted analysis-genotype dosage difference",
    "Missingness-pattern strata",
    "Residual donor-map pattern violations",
    "Maximum observed beta reproduction difference",
    "Maximum observed adjusted-SE reproduction difference",
    "Minimum residual degrees of freedom",
    "Maximum residual degrees of freedom",
    "Raw fit warnings",
    "BF update warnings",
    "Null likelihood warnings"
  ),
  value = c(
    n_pair_discoveries,
    nrow(selection),
    nrow(selection),
    length(unique(selection$gene_id)),
    length(unique(selection$pair_key)),
    selection_discovery_count,
    candidate_count_summary[["minimum"]],
    candidate_count_summary[["median"]],
    candidate_count_summary[["mean"]],
    candidate_count_summary[["maximum"]],
    length(vcf_donors),
    sum(donor_map$fixed_point),
    maximum_crossproduct_difference,
    maximum_sorted_dosage_difference,
    length(unique(donor_observation_patterns)),
    if (permutation_method %in% c(
      "genotype_label_independent_time",
      "signal_stripped_independent_time_residual"
    )) 0L else sum(
      donor_observation_patterns[donor_map$source_donor] !=
        donor_observation_patterns[donor_map$target_donor]
    ),
    maximum_observed_beta_difference,
    maximum_observed_se_difference,
    min(residual_df),
    max(residual_df),
    length(raw_capture$warnings),
    length(bf_capture$warnings),
    length(null_capture$warnings)
  ),
  stringsAsFactors = FALSE
)

source_information_after <- file_metadata(required_paths, required_roles)
if (!identical(source_information_before, source_information_after)) {
  stop("At least one immutable source file changed during the pilot.")
}

analysis_elapsed <- proc.time()[["elapsed"]] - analysis_start
experiment_label <- if (identical(permutation_method, "genotype_label")) {
  "Selected-signal shared-genotype-permutation negative-control pilot"
} else if (identical(
  permutation_method,
  "genotype_label_independent_time"
)) {
  paste(
    "Selected-signal independent-time genotype-label-permutation",
    "negative-control pilot"
  )
} else if (identical(permutation_method, "residual_block")) {
  "Selected-signal synchronized residual-block-permutation negative-control pilot"
} else if (identical(
  permutation_method,
  "signal_stripped_independent_time_residual"
)) {
  paste(
    "Selected-signal independent-time signal-stripped residual-permutation",
    "negative-control pilot"
  )
} else if (identical(
  permutation_method,
  "signal_stripped_unit_specific_residual_block"
)) {
  paste(
    "Selected-signal unit-specific signal-stripped residual-block-permutation",
    "negative-control pilot"
  )
} else {
  paste(
    "Selected-signal synchronized signal-stripped residual-block-permutation",
    "negative-control pilot"
  )
}
if (identical(target_selection_method, "random_all_tested")) {
  experiment_label <- paste(
    "Random-within-discovered-gene",
    experiment_label
  )
} else if (identical(target_selection_method, "random_all_genes")) {
  experiment_label <- paste("Random-within-all-tested-genes", experiment_label)
}
permutation_description <- if (identical(permutation_method, "genotype_label")) {
  paste(
    "One permutation of all 19 genotype-matrix donor rows, shared across",
    "all selected variants and all 16 time points."
  )
} else if (identical(
  permutation_method,
  "genotype_label_independent_time"
)) {
  paste(
    "At each time point, independently permute genotype labels among donors",
    "observed at that time. The time-specific donor map is shared across all",
    "selected units, but the maps differ across time points. Expression, PCs,",
    "and the original genotype signal in the response remain unchanged."
  )
} else if (identical(permutation_method, "residual_block")) {
  paste(
    "One covariate-only residual donor-block permutation within observed",
    "missingness-pattern strata, shared across all selected genes and all",
    "16 time points; residuals are reprojected against each time-specific",
    "PC design before regression on the unchanged genotype matrix."
  )
} else if (identical(
  permutation_method,
  "signal_stripped_independent_time_residual"
)) {
  paste(
    "At each time point, fit expression on genotype, an intercept, and five",
    "PCs; retain the fitted nuisance component, remove the fitted genotype",
    "component, and independently permute the full-model residuals among",
    "donors observed at that time. The time-specific donor map is shared",
    "across all selected genes but donor trajectories are not preserved."
  )
} else if (identical(
  permutation_method,
  "signal_stripped_unit_specific_residual_block"
)) {
  paste(
    "At each time point, fit expression on genotype, an intercept, and five",
    "PCs; retain the fitted nuisance component and remove the fitted genotype",
    "component. For each target unit, apply one independently sampled",
    "within-missingness-stratum donor map to its full-model residual trajectory",
    "across all 16 time points before refitting against unchanged genotype."
  )
} else {
  paste(
    "At each time point, fit expression on genotype, an intercept, and five",
    "PCs; retain the fitted nuisance component, remove the fitted genotype",
    "component, and apply one shared within-missingness-stratum donor map to",
    "the unadjusted full-model residuals before refitting against the unchanged",
    "genotype matrix."
  )
}
target_selection_description <- if (identical(
  target_selection_method,
  "top_discovered"
)) {
  paste(
    "Minimum full-data BF-adjusted lfdr pair within each gene represented",
    "among the original alpha-0.05 pair-level discoveries."
  )
} else {
  if (identical(target_selection_method, "random_all_tested")) paste(
    "One variant sampled uniformly from all tested variants within each gene",
    "represented among the original alpha-0.05 pair-level discoveries; the",
    "sampled variant need not itself be an original discovery."
  ) else paste(
    "One variant sampled uniformly from all tested variants within every",
    "tested gene, without conditioning gene eligibility or variant selection",
    "on the original discovery results."
  )
}
configuration <- list(
  experiment = experiment_label,
  seed = seed,
  alpha = alpha,
  permutation_method = permutation_method,
  permutation = permutation_description,
  residual_source = if (permutation_method %in% c(
    "signal_stripped_residual_block",
    "signal_stripped_independent_time_residual",
    "signal_stripped_unit_specific_residual_block"
  )) {
    "full_model_Y_on_G_and_PCs"
  } else if (identical(permutation_method, "residual_block")) {
    "covariate_only_Y_on_PCs"
  } else {
    NA_character_
  },
  fitted_genotype_removed = permutation_method %in% c(
    "signal_stripped_residual_block",
    "signal_stripped_independent_time_residual",
    "signal_stripped_unit_specific_residual_block"
  ),
  leverage_adjustment = if (permutation_method %in% c(
    "signal_stripped_residual_block",
    "signal_stripped_independent_time_residual",
    "signal_stripped_unit_specific_residual_block"
  )) "none" else NA_character_,
  temporal_residual_alignment_preserved = if (permutation_method %in% c(
    "residual_block",
    "signal_stripped_residual_block",
    "signal_stripped_unit_specific_residual_block"
  )) {
    TRUE
  } else if (identical(
    permutation_method,
    "signal_stripped_independent_time_residual"
  )) {
    FALSE
  } else {
    NA
  },
  temporal_genotype_alignment_preserved = if (identical(
    permutation_method,
    "genotype_label"
  )) {
    TRUE
  } else if (identical(
    permutation_method,
    "genotype_label_independent_time"
  )) {
    FALSE
  } else {
    NA
  },
  cross_unit_genotype_alignment_preserved = if (permutation_method %in% c(
    "genotype_label",
    "genotype_label_independent_time"
  )) {
    TRUE
  } else {
    NA
  },
  cross_unit_residual_alignment_preserved = if (permutation_method %in% c(
    "residual_block",
    "signal_stripped_residual_block",
    "signal_stripped_independent_time_residual"
  )) {
    TRUE
  } else if (identical(
    permutation_method,
    "signal_stripped_unit_specific_residual_block"
  )) {
    FALSE
  } else {
    NA
  },
  source_residual_genotype_correlation_by_time =
    source_residual_genotype_correlation,
  maximum_source_residual_genotype_correlation = if (permutation_method %in% c(
    "signal_stripped_residual_block",
    "signal_stripped_independent_time_residual",
    "signal_stripped_unit_specific_residual_block"
  )) {
    max(source_residual_genotype_correlation)
  } else {
    NA_real_
  },
  target_selection_method = target_selection_method,
  target_selection = target_selection_description,
  selection_seed = if (target_selection_method %in% c(
    "random_all_tested",
    "random_all_genes"
  )) selection_seed else NA_integer_,
  selection_candidate_pool = attr(selection, "candidate_pool"),
  n_target_units_that_were_original_discoveries = selection_discovery_count,
  candidate_variant_count_summary = candidate_count_summary,
  donor_map_source_path = donor_map_source_path,
  donor_map_source_md5 = if (nzchar(donor_map_source_path)) {
    unname(tools::md5sum(donor_map_source_path))
  } else {
    NA_character_
  },
  n_pair_level_discoveries = n_pair_discoveries,
  n_target_units = nrow(selection),
  n_permuted_null_units = nrow(selection),
  time_grid = time_grid,
  sample_counts = sample_counts,
  residual_df = residual_df,
  donor_observation_patterns = donor_observation_patterns,
  original_settings = original_settings,
  original_psd_grid = original_psd_grid,
  num_cores = num_cores,
  likelihood_source = likelihood_source,
  prepared_null_input_md5 = prepared_null_input_md5,
  precomputed_null_fit_md5 = precomputed_null_fit_md5,
  null_likelihood_elapsed_seconds = unname(null_likelihood_elapsed),
  raw_refit_elapsed_seconds = unname(raw_elapsed),
  bf_update_elapsed_seconds = unname(bf_elapsed),
  total_elapsed_seconds = unname(analysis_elapsed),
  null_likelihood_warnings = null_capture$warnings,
  raw_refit_warnings = raw_capture$warnings,
  bf_update_warnings = bf_capture$warnings,
  source_information = source_information_before,
  r_version = R.version.string,
  package_versions = c(fashr = as.character(utils::packageVersion("fashr"))),
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  interpretation_boundary = paste(
    "The permuted-null call rate is a conditional negative-control rejection",
    "rate. The target-set FDR quantity is post-selection diagnostic, not",
    "formal FDR."
  )
)
fit_bundle <- list(
  configuration = configuration,
  selection = selection,
  donor_map = donor_map,
  raw_fit = merged_raw,
  bf_adjusted_fit = merged_bf
)
effect_estimates <- list(
  time_grid = time_grid,
  target_beta = target_beta,
  target_adjusted_se = target_adjusted_se,
  observed_refit_beta = observed_refit_beta,
  observed_refit_raw_se = observed_refit_raw_se,
  observed_refit_adjusted_se = observed_refit_adjusted_se,
  permuted_beta = permuted_beta,
  permuted_raw_se = permuted_raw_se,
  permuted_adjusted_se = permuted_adjusted_se,
  source_residual_genotype_correlation_by_time =
    source_residual_genotype_correlation,
  residual_df = residual_df,
  sample_counts = sample_counts
)

saveRDS(configuration, file.path(staging_directory, "configuration.rds"))
saveRDS(fit_bundle, file.path(staging_directory, "merged_fash_fit.rds"))
saveRDS(effect_estimates, file.path(staging_directory, "effect_estimates.rds"))
write_csv(selection, file.path(staging_directory, "selection.csv"))
write_csv(donor_map, file.path(staging_directory, "donor_permutation.csv"))
write_csv(unit_lfdr, file.path(staging_directory, "unit_lfdr.csv"))
write_csv(
  calibration_diagnostics,
  file.path(staging_directory, "calibration_diagnostics.csv")
)
write_csv(lfdr_quantiles, file.path(staging_directory, "lfdr_quantiles.csv"))
write_csv(prior_weights, file.path(staging_directory, "prior_weights.csv"))
write_csv(validation, file.path(staging_directory, "validation.csv"))
write_csv(
  source_information_before,
  file.path(staging_directory, "source_information.csv")
)

if (!file.rename(staging_directory, final_directory)) {
  stop("Could not finalize output directory: ", final_directory)
}
staging_complete <- TRUE

cat("\nSelected-signal matched-null permutation pilot completed.\n")
cat("Output: ", final_directory, "\n", sep = "")
print(calibration_diagnostics, row.names = FALSE)
