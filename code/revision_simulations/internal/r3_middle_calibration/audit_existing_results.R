#!/usr/bin/env Rscript

# Audit the completed R3A Middle calibration curves and alpha-0.05 false calls.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.", call. = FALSE)
}

expect_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))

helper_source <- paste(deparse(body(compute_functional_lfsr)), collapse = "\n")
expect_true(
  grepl("lfsr_cal = function(x) mean(x <= 0)", helper_source, fixed = TRUE),
  "The functional helper is not using the expected one-sided posterior error."
)

cache_ids <- c(
  closed = paste0(
    "r3_real_genotype_one_per_gene_J6362_",
    "matched_functional_relative_clearance_main_effect_fashr0143_pilot5"
  ),
  open = paste0(
    "r3_real_genotype_one_per_gene_J6362_",
    "matched_functional_open_middle_3_12_",
    "relative_clearance_main_effect_fashr0143_pilot5"
  ),
  center_aligned = paste0(
    "r3_real_genotype_one_per_gene_J6362_",
    "matched_functional_open_middle_3_12_center_aligned_",
    "relative_clearance_main_effect_fashr0143_pilot5"
  )
)
seed_list <- c(12345L, 22345L, 32345L, 42345L, 52345L)
output_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "diagnostics", "r3_middle_calibration"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_replicates <- function(cache_id, mechanism = "random_bspline") {
  replicate_dir <- file.path(
    workflowr_root,
    "output", "revision_simulations", "mc", cache_id, "replicates"
  )
  paths <- file.path(
    replicate_dir,
    sprintf("%s_seed_%d.rds", mechanism, seed_list)
  )
  expect_true(all(file.exists(paths)), "At least one required replicate is missing.")
  lapply(paths, readRDS)
}

result_rows <- lapply(names(cache_ids), function(result_name) {
  cache_id <- cache_ids[[result_name]]
  cache_dir <- file.path(
    workflowr_root,
    "output", "revision_simulations", "mc", cache_id
  )
  configuration <- readRDS(file.path(cache_dir, "configuration.rds"))
  middle_window <- if (is.null(configuration$middle_window)) {
    c(4, 11)
  } else {
    configuration$middle_window
  }
  middle_boundary <- if (is.null(configuration$middle_boundary)) {
    "closed"
  } else {
    configuration$middle_boundary
  }
  middle_membership <- temporal_middle_membership(
    configuration$evaluation_grid,
    middle_window = middle_window,
    middle_boundary = middle_boundary
  )
  definition <- if (middle_boundary == "open") {
    sprintf("%g < t < %g", middle_window[[1L]], middle_window[[2L]])
  } else {
    sprintf("%g <= t <= %g", middle_window[[1L]], middle_window[[2L]])
  }
  replicates <- read_replicates(cache_id)
  alpha_rows <- do.call(rbind, lapply(replicates, `[[`, "functional_alpha"))
  alpha_rows <- alpha_rows[
    alpha_rows$target == "middle" &
      alpha_rows$method %in% c("FASH-IWP1-Raw", "FASH-IWP1-BF"),
    ,
    drop = FALSE
  ]
  split_rows <- split(
    alpha_rows,
    interaction(alpha_rows$method, alpha_rows$alpha, drop = TRUE)
  )
  do.call(rbind, lapply(split_rows, function(rows) {
    empirical_fsr <- rows$empirical_fsr
    mc_se <- if (length(empirical_fsr) > 1L) {
      stats::sd(empirical_fsr) / sqrt(length(empirical_fsr))
    } else {
      NA_real_
    }
    critical_value <- stats::qt(0.975, df = length(empirical_fsr) - 1L)
    data.frame(
      result = result_name,
      cache_id = cache_id,
      middle_definition = definition,
      middle_grid_points = sum(middle_membership),
      outside_grid_points = sum(!middle_membership),
      method = rows$method[[1L]],
      alpha = rows$alpha[[1L]],
      replications = nrow(rows),
      mean_discoveries = mean(rows$n_discoveries),
      total_discoveries = sum(rows$n_discoveries),
      mean_false_discoveries = mean(rows$false_discoveries),
      total_false_discoveries = sum(rows$false_discoveries),
      mean_first_stage_null_calls = mean(rows$first_stage_null_calls),
      mean_conditional_false_discoveries =
        mean(rows$conditional_false_discoveries),
      mean_estimated_fsr = mean(rows$estimated_fsr),
      mean_empirical_fsr = mean(empirical_fsr),
      pooled_empirical_fsr =
        sum(rows$false_discoveries) / sum(rows$n_discoveries),
      replicate_min_fsr = min(empirical_fsr),
      replicate_max_fsr = max(empirical_fsr),
      empirical_fsr_mc_se = mc_se,
      empirical_fsr_ci_lower = max(
        0,
        mean(empirical_fsr) - critical_value * mc_se
      ),
      empirical_fsr_ci_upper = min(
        1,
        mean(empirical_fsr) + critical_value * mc_se
      ),
      stringsAsFactors = FALSE
    )
  }))
})
existing_comparison <- do.call(rbind, result_rows)
rownames(existing_comparison) <- NULL
existing_comparison <- existing_comparison[
  order(
    existing_comparison$result,
    existing_comparison$method,
    existing_comparison$alpha
  ),
  ,
  drop = FALSE
]
utils::write.csv(
  existing_comparison,
  file.path(output_dir, "existing_result_comparison.csv"),
  row.names = FALSE
)

aligned_replicates <- read_replicates(cache_ids[["center_aligned"]])
call_diagnostics <- do.call(rbind, lapply(
  aligned_replicates,
  `[[`,
  "call_diagnostics_alpha005"
))
false_middle <- call_diagnostics[
  call_diagnostics$method == "FASH-IWP1-BF" &
    call_diagnostics$target == "middle" &
    call_diagnostics$false_discovery,
  ,
  drop = FALSE
]
breaks <- c(-Inf, -0.50, -0.25, -0.15, -0.10, -1e-12, Inf)
labels <- c(
  "less_than_-0.50",
  "-0.50_to_-0.25",
  "-0.25_to_-0.15",
  "-0.15_to_-0.10",
  "-0.10_to_zero",
  "zero_or_positive"
)
false_middle$truth_distance_bin <- cut(
  false_middle$true_functional,
  breaks = breaks,
  labels = labels,
  right = FALSE
)
false_call_decomposition <- as.data.frame(table(
  factor(false_middle$truth_distance_bin, levels = labels),
  factor(false_middle$truth_group),
  useNA = "ifany"
))
names(false_call_decomposition) <- c(
  "truth_distance_bin", "truth_group", "false_calls"
)
false_call_decomposition <- false_call_decomposition[
  false_call_decomposition$false_calls > 0,
  ,
  drop = FALSE
]
false_call_decomposition$total_false_calls <- nrow(false_middle)
false_call_decomposition$proportion <-
  false_call_decomposition$false_calls / nrow(false_middle)
utils::write.csv(
  false_call_decomposition,
  file.path(output_dir, "false_call_decomposition.csv"),
  row.names = FALSE
)

cat("Validated one-sided functional posterior error: P(F <= 0 | data).\n")
cat("Wrote existing R3 comparison and false-call decomposition to:\n")
cat(normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
