#!/usr/bin/env Rscript

# Run cached-likelihood small-M matched-null empirical-Bayes refits.

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

write_csv <- function(value, path) {
  utils::write.csv(value, path, row.names = FALSE, na = "")
}

capture_warnings <- function(expression) {
  warnings <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

make_input_provenance <- function(paths, roles) {
  information <- file.info(paths)
  data.frame(
    role = roles,
    path = normalizePath(paths, mustWork = TRUE),
    size_bytes = as.numeric(information$size),
    modified_at = format(information$mtime, tz = "UTC", usetz = TRUE),
    md5 = unname(tools::md5sum(paths)),
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
script_directory <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation"
)
source(file.path(
  script_directory,
  "selected_signal_genotype_permutation_helpers.R"
))
source(file.path(script_directory, "small_m_eb_sensitivity_helpers.R"))

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required.")
}

full_output_id <- paste0(
  "all_gene_random_variant_signal_stripped_residual_block_permutation_",
  "selection20260817_seed20260811"
)
full_output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  full_output_id
)
full_fit_path <- file.path(full_output_directory, "merged_fash_fit.rds")
full_calibration_path <- file.path(
  full_output_directory,
  "calibration_diagnostics.csv"
)
full_selection_path <- file.path(full_output_directory, "selection.csv")
null_fit_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  paste0(
    "all_gene_random_variant_signal_stripped_residual_null_fit_",
    "selection20260817_seed20260811.rds"
  )
)
input_paths <- c(
  full_fit_path,
  full_calibration_path,
  full_selection_path,
  null_fit_path
)
input_roles <- c(
  "full_merged_fit",
  "full_calibration",
  "target_selection",
  "clean_session_null_likelihood_fit"
)
if (any(!file.exists(input_paths))) {
  stop("At least one required small-M sensitivity input is missing.")
}

output_id <- paste0(
  "all_gene_random_variant_small_m_signal_stripped_residual_sensitivity_",
  "selection20260817_seed20260811_subsets20260819"
)
output_parent <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal"
)
final_directory <- file.path(output_parent, output_id)
if (dir.exists(final_directory)) {
  stop("The completed small-M sensitivity output already exists: ",
       final_directory)
}
staging_directory <- paste0(
  final_directory,
  ".staging_",
  Sys.getpid()
)
if (dir.exists(staging_directory)) {
  stop("The small-M staging directory already exists: ", staging_directory)
}
dir.create(staging_directory, recursive = TRUE, showWarnings = FALSE)
staging_complete <- FALSE
on.exit({
  if (!staging_complete && dir.exists(staging_directory)) {
    unlink(staging_directory, recursive = TRUE)
  }
}, add = TRUE)

analysis_start <- proc.time()[["elapsed"]]
input_provenance <- make_input_provenance(input_paths, input_roles)
full_bundle <- readRDS(full_fit_path)
null_bundle <- readRDS(null_fit_path)
selection <- utils::read.csv(full_selection_path, stringsAsFactors = FALSE)
saved_calibration <- utils::read.csv(
  full_calibration_path,
  stringsAsFactors = FALSE
)

J <- as.integer(full_bundle$configuration$n_target_units)
expected_J <- 6362L
raw_full <- full_bundle$raw_fit
bf_full <- full_bundle$bf_adjusted_fit
all_keys <- rownames(raw_full$L_matrix)
full_data_list <- raw_full$fash_data$data_list
full_se_list <- raw_full$fash_data$S
if (!identical(full_bundle$configuration$target_selection_method,
               "random_all_genes") ||
    !identical(full_bundle$configuration$permutation_method,
               "signal_stripped_residual_block") ||
    J != expected_J || nrow(selection) != J ||
    !identical(null_bundle$cache_version, "selected_signal_null_fit_v1") ||
    length(null_bundle$null_unit_keys) != J ||
    nrow(raw_full$L_matrix) != 2L * J ||
    ncol(raw_full$L_matrix) != 52L ||
    length(raw_full$lfdr) != 2L * J ||
    length(bf_full$lfdr) != 2L * J ||
    length(full_data_list) != 2L * J ||
    length(full_se_list) != 2L * J ||
    length(all_keys) != 2L * J || anyDuplicated(all_keys) ||
    !identical(names(full_data_list), all_keys) ||
    !identical(names(full_se_list), all_keys) ||
    !identical(all_keys[seq_len(J)], selection$pair_key) ||
    !identical(all_keys[J + seq_len(J)], null_bundle$null_unit_keys) ||
    !isTRUE(all.equal(raw_full$psd_grid, bf_full$psd_grid, tolerance = 0))) {
  stop("The cached all-gene synchronized residual experiment is invalid.")
}

ratio_to_size <- c(
  "0.05" = as.integer(round(0.05 * J)),
  "0.10" = as.integer(round(0.10 * J)),
  "0.20" = as.integer(round(0.20 * J))
)
n_replicates <- 20L
subset_seed <- 20260819L
alpha_grid <- seq(0.001, 0.200, by = 0.001)
subset_membership <- make_nested_null_subsets(
  n_null = J,
  ratio_to_size = ratio_to_size,
  n_replicates = n_replicates,
  seed = subset_seed
)
subset_membership$null_unit_key <- all_keys[
  J + subset_membership$null_index
]

fit_one_subset <- function(indices) {
  unit_keys <- all_keys[indices]
  raw_start <- proc.time()[["elapsed"]]
  raw_capture <- capture_warnings(refit_fash_from_likelihood(
    source_fit = raw_full,
    data_list = full_data_list[indices],
    se_list = full_se_list[indices],
    likelihood_matrix = raw_full$L_matrix[indices, , drop = FALSE],
    unit_keys = unit_keys,
    penalty = raw_full$settings$penalty
  ))
  raw_seconds <- proc.time()[["elapsed"]] - raw_start
  bf_start <- proc.time()[["elapsed"]]
  bf_capture <- capture_warnings(fashr::BF_update(
    raw_capture$value,
    plot = FALSE
  ))
  bf_seconds <- proc.time()[["elapsed"]] - bf_start
  names(bf_capture$value$lfdr) <- unit_keys
  list(
    raw_fit = raw_capture$value,
    bf_fit = bf_capture$value,
    raw_seconds = raw_seconds,
    bf_seconds = bf_seconds,
    raw_warnings = raw_capture$warnings,
    bf_warnings = bf_capture$warnings
  )
}

make_pi0_rows <- function(raw_fit,
                          bf_fit,
                          m_ratio,
                          m_size,
                          replicate_id,
                          is_reference) {
  fits <- list("Raw" = raw_fit, "BF-adjusted" = bf_fit)
  rows <- lapply(names(fits), function(fit_stage) {
    pi0_merged <- extract_pi0(fits[[fit_stage]])
    pi0_target_unbounded <- if (m_size == 0L) {
      pi0_merged
    } else {
      ((J + m_size) * pi0_merged - m_size) / J
    }
    data.frame(
      replicate_id = as.integer(replicate_id),
      m_ratio = as.numeric(m_ratio),
      m_size = as.integer(m_size),
      fit_stage = fit_stage,
      pi0_merged = pi0_merged,
      pi0_target_unbounded = pi0_target_unbounded,
      pi0_target_valid = pi0_target_unbounded >= 0 &&
        pi0_target_unbounded <= 1,
      pi0_target_bounded = min(1, max(0, pi0_target_unbounded)),
      is_reference = is_reference,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_prior_rows <- function(raw_fit,
                            bf_fit,
                            m_ratio,
                            m_size,
                            replicate_id,
                            is_reference) {
  fits <- list("Raw" = raw_fit, "BF-adjusted" = bf_fit)
  rows <- lapply(names(fits), function(fit_stage) {
    prior <- fits[[fit_stage]]$prior_weights
    data.frame(
      replicate_id = as.integer(replicate_id),
      m_ratio = as.numeric(m_ratio),
      m_size = as.integer(m_size),
      fit_stage = fit_stage,
      psd = as.numeric(prior$psd),
      prior_weight = as.numeric(prior$prior_weight),
      component_type = ifelse(prior$psd == 0, "Exact null", "Alternative"),
      is_reference = is_reference,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

message("Refitting the target-only M/J = 0 reference.")
target_reference <- fit_one_subset(seq_len(J))
pi0_rows <- list(make_pi0_rows(
  target_reference$raw_fit,
  target_reference$bf_fit,
  m_ratio = 0,
  m_size = 0L,
  replicate_id = 0L,
  is_reference = TRUE
))
prior_rows <- list(make_prior_rows(
  target_reference$raw_fit,
  target_reference$bf_fit,
  m_ratio = 0,
  m_size = 0L,
  replicate_id = 0L,
  is_reference = TRUE
))
runtime_rows <- list(data.frame(
  replicate_id = 0L,
  m_ratio = 0,
  m_size = 0L,
  raw_seconds = target_reference$raw_seconds,
  bf_seconds = target_reference$bf_seconds,
  total_seconds = target_reference$raw_seconds + target_reference$bf_seconds,
  raw_warning_count = length(target_reference$raw_warnings),
  bf_warning_count = length(target_reference$bf_warnings),
  is_reference = TRUE,
  stringsAsFactors = FALSE
))

pi0_rows[[length(pi0_rows) + 1L]] <- make_pi0_rows(
  raw_full,
  bf_full,
  m_ratio = 1,
  m_size = J,
  replicate_id = 0L,
  is_reference = TRUE
)
prior_rows[[length(prior_rows) + 1L]] <- make_prior_rows(
  raw_full,
  bf_full,
  m_ratio = 1,
  m_size = J,
  replicate_id = 0L,
  is_reference = TRUE
)

curve_rows <- list()
alpha005_rows <- list()
full_group <- rep(c("target", "permuted_null"), each = J)
for (fit_stage in c("Raw", "BF-adjusted")) {
  reference_fit <- if (fit_stage == "Raw") raw_full else bf_full
  reference_curve <- summarize_small_m_curve(
    lfdr = reference_fit$lfdr,
    group = full_group,
    pi0_merged = extract_pi0(reference_fit),
    fit_stage = fit_stage,
    m_ratio = 1,
    replicate_id = 0L,
    alpha_grid = alpha_grid
  )
  reference_curve$is_reference <- TRUE
  curve_rows[[length(curve_rows) + 1L]] <- reference_curve
  alpha005_rows[[length(alpha005_rows) + 1L]] <- reference_curve[
    abs(reference_curve$nominal_alpha - 0.05) < 1e-12,
    ,
    drop = FALSE
  ]
}

small_fit_counter <- 0L
for (replicate_id in seq_len(n_replicates)) {
  message(
    "Running small-M replicate ", replicate_id, " of ", n_replicates, "."
  )
  for (m_ratio in as.numeric(names(ratio_to_size))) {
    m_size <- unname(ratio_to_size[
      match(m_ratio, as.numeric(names(ratio_to_size)))
    ])
    selected_null_rows <- subset_membership[
      subset_membership$replicate_id == replicate_id &
        abs(subset_membership$m_ratio - m_ratio) < 1e-12,
      ,
      drop = FALSE
    ]
    if (nrow(selected_null_rows) != m_size ||
        anyDuplicated(selected_null_rows$null_index)) {
      stop("The deterministic null subset is invalid.")
    }
    indices <- c(seq_len(J), J + selected_null_rows$null_index)
    subset_fit <- fit_one_subset(indices)
    small_fit_counter <- small_fit_counter + 1L
    pi0_rows[[length(pi0_rows) + 1L]] <- make_pi0_rows(
      subset_fit$raw_fit,
      subset_fit$bf_fit,
      m_ratio = m_ratio,
      m_size = m_size,
      replicate_id = replicate_id,
      is_reference = FALSE
    )
    prior_rows[[length(prior_rows) + 1L]] <- make_prior_rows(
      subset_fit$raw_fit,
      subset_fit$bf_fit,
      m_ratio = m_ratio,
      m_size = m_size,
      replicate_id = replicate_id,
      is_reference = FALSE
    )
    runtime_rows[[length(runtime_rows) + 1L]] <- data.frame(
      replicate_id = replicate_id,
      m_ratio = m_ratio,
      m_size = m_size,
      raw_seconds = subset_fit$raw_seconds,
      bf_seconds = subset_fit$bf_seconds,
      total_seconds = subset_fit$raw_seconds + subset_fit$bf_seconds,
      raw_warning_count = length(subset_fit$raw_warnings),
      bf_warning_count = length(subset_fit$bf_warnings),
      is_reference = FALSE,
      stringsAsFactors = FALSE
    )
    subset_group <- c(rep("target", J), rep("permuted_null", m_size))
    for (fit_stage in c("Raw", "BF-adjusted")) {
      current_fit <- if (fit_stage == "Raw") {
        subset_fit$raw_fit
      } else {
        subset_fit$bf_fit
      }
      current_curve <- summarize_small_m_curve(
        lfdr = current_fit$lfdr,
        group = subset_group,
        pi0_merged = extract_pi0(current_fit),
        fit_stage = fit_stage,
        m_ratio = m_ratio,
        replicate_id = replicate_id,
        alpha_grid = alpha_grid
      )
      current_curve$is_reference <- FALSE
      curve_rows[[length(curve_rows) + 1L]] <- current_curve
      alpha005_rows[[length(alpha005_rows) + 1L]] <- current_curve[
        abs(current_curve$nominal_alpha - 0.05) < 1e-12,
        ,
        drop = FALSE
      ]
    }
    rm(subset_fit)
    gc(verbose = FALSE)
  }
}

replicate_pi0 <- do.call(rbind, pi0_rows)
replicate_prior_weights <- do.call(rbind, prior_rows)
replicate_alpha_curves <- do.call(rbind, curve_rows)
replicate_alpha005 <- do.call(rbind, alpha005_rows)
runtime <- do.call(rbind, runtime_rows)
rownames(replicate_pi0) <- NULL
rownames(replicate_prior_weights) <- NULL
rownames(replicate_alpha_curves) <- NULL
rownames(replicate_alpha005) <- NULL
rownames(runtime) <- NULL

saved_bf <- saved_calibration[
  saved_calibration$fit_stage == "BF-adjusted",
  ,
  drop = FALSE
]
full_bf_alpha005 <- replicate_alpha005[
  replicate_alpha005$is_reference &
    replicate_alpha005$m_ratio == 1 &
    replicate_alpha005$fit_stage == "BF-adjusted",
  ,
  drop = FALSE
]
expected_small_combinations <- n_replicates * length(ratio_to_size) * 2L
prior_sums <- aggregate(
  prior_weight ~ replicate_id + m_ratio + fit_stage + is_reference,
  data = replicate_prior_weights,
  FUN = sum
)
prior_groups <- split(
  replicate_prior_weights,
  interaction(
    replicate_prior_weights$replicate_id,
    replicate_prior_weights$m_ratio,
    replicate_prior_weights$fit_stage,
    replicate_prior_weights$is_reference,
    drop = TRUE
  )
)
prior_component_contract_valid <- all(vapply(
  prior_groups,
  function(prior) {
    sum(prior$psd == 0) == 1L &&
      !anyDuplicated(prior$psd) &&
      all(is.finite(prior$prior_weight)) &&
      all(prior$prior_weight >= 0)
  },
  logical(1)
))
validation <- data.frame(
  check = c(
    "expected_small_fit_count",
    "expected_small_pi0_combinations",
    "expected_small_alpha005_combinations",
    "expected_small_curve_rows",
    "valid_active_mixture_components",
    "prior_weights_sum_to_one",
    "all_pi0_merged_valid",
    "full_reference_matches_saved_calls",
    "full_reference_matches_saved_target_plugin",
    "full_reference_matches_saved_merged_plugin",
    "no_refit_warnings"
  ),
  passed = c(
    small_fit_counter == n_replicates * length(ratio_to_size),
    sum(!replicate_pi0$is_reference) == expected_small_combinations,
    sum(!replicate_alpha005$is_reference) == expected_small_combinations,
    sum(!replicate_alpha_curves$is_reference) ==
      expected_small_combinations * length(alpha_grid),
    prior_component_contract_valid,
    all(abs(prior_sums$prior_weight - 1) < 1e-8),
    all(replicate_pi0$pi0_merged >= 0 & replicate_pi0$pi0_merged <= 1),
    nrow(saved_bf) == 1L && nrow(full_bf_alpha005) == 1L &&
      full_bf_alpha005$target_calls == saved_bf$target_calls &&
      full_bf_alpha005$permuted_null_calls ==
        saved_bf$permuted_null_calls,
    nrow(saved_bf) == 1L && nrow(full_bf_alpha005) == 1L &&
      abs(
        full_bf_alpha005$target_pi0_plugin_fdr -
          saved_bf$post_selection_fdr_target_from_pi0
      ) < 1e-14,
    nrow(saved_bf) == 1L && nrow(full_bf_alpha005) == 1L &&
      abs(
        full_bf_alpha005$merged_pi0_plugin_fdr -
          saved_bf$scaled_fdr_merged_from_estimated_pi0
      ) < 1e-14,
    all(runtime$raw_warning_count == 0L) &&
      all(runtime$bf_warning_count == 0L)
  ),
  stringsAsFactors = FALSE
)
if (any(!validation$passed)) {
  print(validation, row.names = FALSE)
  stop("At least one small-M sensitivity validation check failed.")
}

reference_summaries <- rbind(
  transform(
    replicate_pi0[replicate_pi0$is_reference, ],
    nominal_alpha = NA_real_,
    target_calls = NA_integer_,
    permuted_null_calls = NA_integer_,
    merged_calls = NA_integer_,
    permuted_null_call_rate = NA_real_,
    target_pi0_plugin_fdr = NA_real_,
    merged_pi0_plugin_fdr = NA_real_
  ),
  transform(
    replicate_alpha005[replicate_alpha005$is_reference, c(
      "replicate_id", "m_ratio", "m_size", "fit_stage", "pi0_merged",
      "pi0_target_unbounded", "pi0_target_valid", "pi0_target_bounded",
      "nominal_alpha", "target_calls", "permuted_null_calls",
      "merged_calls", "permuted_null_call_rate",
      "target_pi0_plugin_fdr", "merged_pi0_plugin_fdr"
    )],
    is_reference = TRUE
  )
)
reference_summaries <- reference_summaries[, c(
  "replicate_id", "m_ratio", "m_size", "fit_stage", "pi0_merged",
  "pi0_target_unbounded", "pi0_target_valid", "pi0_target_bounded",
  "is_reference", "nominal_alpha", "target_calls",
  "permuted_null_calls", "merged_calls", "permuted_null_call_rate",
  "target_pi0_plugin_fdr", "merged_pi0_plugin_fdr"
)]

total_seconds <- proc.time()[["elapsed"]] - analysis_start
configuration <- list(
  experiment = "Small-M synchronized signal-stripped residual sensitivity",
  output_id = output_id,
  target_selection_method = full_bundle$configuration$target_selection_method,
  permutation_method = full_bundle$configuration$permutation_method,
  source_selection_seed = full_bundle$configuration$selection_seed,
  source_permutation_seed = full_bundle$configuration$seed,
  subset_seed = subset_seed,
  n_target = J,
  m_ratios = as.numeric(names(ratio_to_size)),
  m_sizes = unname(ratio_to_size),
  n_replicates = n_replicates,
  nested_subsets_within_replicate = TRUE,
  alpha_grid = alpha_grid,
  primary_fit_stage = "BF-adjusted",
  full_psd_grid = raw_full$psd_grid,
  input_provenance = input_provenance,
  total_seconds = total_seconds,
  r_version = R.version.string,
  fashr_version = as.character(utils::packageVersion("fashr")),
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  interpretation_boundary = paste(
    "All small-M refits reuse subsets of one fixed synchronized residual",
    "permutation. V/M is a conditional permuted-control rejection rate,",
    "and the plug-in FDR quantities are sensitivity diagnostics rather than",
    "formal repeated-sampling calibration estimates."
  )
)

saveRDS(configuration, file.path(staging_directory, "configuration.rds"))
write_csv(input_provenance, file.path(
  staging_directory,
  "input_provenance.csv"
))
write_csv(subset_membership, file.path(
  staging_directory,
  "subset_membership.csv"
))
write_csv(replicate_pi0, file.path(staging_directory, "replicate_pi0.csv"))
write_csv(replicate_alpha005, file.path(
  staging_directory,
  "replicate_alpha005.csv"
))
write_csv(replicate_alpha_curves, file.path(
  staging_directory,
  "replicate_alpha_curves.csv"
))
write_csv(replicate_prior_weights, file.path(
  staging_directory,
  "replicate_prior_weights.csv"
))
write_csv(reference_summaries, file.path(
  staging_directory,
  "reference_summaries.csv"
))
write_csv(runtime, file.path(staging_directory, "runtime.csv"))
write_csv(validation, file.path(staging_directory, "validation.csv"))

if (!file.rename(staging_directory, final_directory)) {
  stop("Could not finalize the small-M sensitivity output directory.")
}
staging_complete <- TRUE

cat("Small-M matched-null EB sensitivity completed.\n")
cat("Output: ", final_directory, "\n", sep = "")
cat(sprintf("Total elapsed seconds: %.2f\n", total_seconds))
cat("BF-adjusted alpha-0.05 medians by M/J:\n")
bf_alpha005 <- replicate_alpha005[
  !replicate_alpha005$is_reference &
    replicate_alpha005$fit_stage == "BF-adjusted",
  ,
  drop = FALSE
]
print(aggregate(
  cbind(
    pi0_merged,
    pi0_target_unbounded,
    permuted_null_call_rate,
    target_pi0_plugin_fdr,
    merged_pi0_plugin_fdr
  ) ~ m_ratio,
  data = bf_alpha005,
  FUN = stats::median,
  na.rm = TRUE
), row.names = FALSE)
