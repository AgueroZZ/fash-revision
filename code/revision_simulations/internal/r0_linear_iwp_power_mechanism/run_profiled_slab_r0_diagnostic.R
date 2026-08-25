#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr project root.")
}

summarize_alpha005 <- function(result_table, se_regime) {
  groups <- split(result_table, result_table$method)
  rows <- lapply(groups, function(x) {
    discoveries <- sum(x$selected)
    false_discoveries <- sum(x$selected & x$true_null)
    true_alternatives <- sum(!x$true_null)
    true_positives <- sum(x$selected & !x$true_null)
    data.frame(
      se_regime = se_regime,
      method = x$method[[1L]],
      n_true_alternatives = true_alternatives,
      n_discoveries = discoveries,
      false_discoveries = false_discoveries,
      true_positives = true_positives,
      power = true_positives / true_alternatives,
      realized_fdp = if (discoveries == 0L) 0 else {
        false_discoveries / discoveries
      },
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

seed <- 12345L
J <- 1000L
time_grid <- 0:15
class_probs <- c(constant = 0.80, linear = 0.20)
se_specs <- list(
  medium = c(0.3, 0.5, 0.8),
  high = c(0.5, 0.8, 1.2),
  very_high = c(0.8, 1.2, 1.6)
)
sigma_beta_grid <- exp(seq(log(0.05), log(5), length.out = 25))
matched_grid <- default_revision_grid()
alpha <- 0.05

reference_cache_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r0_all_linear_summary_statistic_seed12345_quick",
  "analysis_cache.rds"
)
reference_cache <- readRDS(reference_cache_path)
if (!identical(reference_cache$configuration$seed, seed) ||
    !identical(reference_cache$configuration$J, J) ||
    !identical(reference_cache$configuration$time_grid, time_grid) ||
    !isTRUE(all.equal(
      reference_cache$configuration$class_probs,
      class_probs,
      tolerance = 0
    )) ||
    !isTRUE(all.equal(
      reference_cache$configuration$se_specs,
      se_specs,
      tolerance = 0
    ))) {
  stop("The retained R0 cache does not match the requested diagnostic.")
}

fit_regime <- function(regime) {
  datasets <- simulate_revision_datasets(
    J = J,
    class_probs = class_probs,
    time_grid = time_grid,
    sd_values = se_specs[[regime]],
    scenario = paste0("r0_all_linear_", regime, "_se"),
    exact_class_counts = TRUE,
    seed = seed
  )
  unit_info <- attr(datasets, "unit_info")

  profiled_raw <- fit_simplified_fash(
    datasets = datasets,
    estimate_sigma = TRUE,
    sigma_beta_grid = sigma_beta_grid,
    scale_time = TRUE
  )
  validate_simplified_sigma_profile(profiled_raw)
  profiled_bf <- BF_update_simplified_fash(profiled_raw)

  matched_raw <- fit_linear_mixture_fash(
    datasets = datasets,
    grid = matched_grid,
    pred_step = 1,
    penalty = 10
  )
  matched_bf <- BF_update_linear_mixture_fash(matched_raw)

  result_table <- rbind(
    evaluate_lfdr_method(
      profiled_raw$lfdr,
      unit_info,
      "FASH-linear-profiled-Raw",
      "dynamic",
      alpha
    ),
    evaluate_lfdr_method(
      profiled_bf$lfdr,
      unit_info,
      "FASH-linear-profiled-BF",
      "dynamic",
      alpha
    ),
    evaluate_lfdr_method(
      matched_raw$lfdr,
      unit_info,
      "FASH-linear-mixture-Raw",
      "dynamic",
      alpha
    ),
    evaluate_lfdr_method(
      matched_bf$lfdr,
      unit_info,
      "FASH-linear-mixture-BF",
      "dynamic",
      alpha
    )
  )
  result_table$se_regime <- regime

  prior_summary <- data.frame(
    se_regime = regime,
    method = c(
      "FASH-linear-profiled-Raw",
      "FASH-linear-profiled-BF",
      "FASH-linear-mixture-Raw",
      "FASH-linear-mixture-BF"
    ),
    estimated_pi0 = c(
      profiled_raw$prior_weights$prior_weight[[1L]],
      profiled_bf$prior_weights$prior_weight[[1L]],
      constant_component_prior_weight(matched_raw),
      constant_component_prior_weight(matched_bf)
    ),
    selected_endpoint_slope_sd = c(
      profiled_raw$sigma_beta,
      profiled_bf$sigma_beta,
      NA_real_,
      NA_real_
    ),
    stringsAsFactors = FALSE
  )

  matched_prior <- rbind(
    transform(
      extract_linear_mixture_prior_table(matched_raw, seed, "Raw"),
      se_regime = regime
    ),
    transform(
      extract_linear_mixture_prior_table(matched_bf, seed, "BF"),
      se_regime = regime
    )
  )

  list(
    result_table = result_table,
    summary = summarize_alpha005(result_table, regime),
    prior_summary = prior_summary,
    matched_prior = matched_prior,
    sigma_profile = transform(profiled_raw$sigma_profile, se_regime = regime)
  )
}

start_time <- proc.time()[["elapsed"]]
results <- lapply(names(se_specs), fit_regime)
names(results) <- names(se_specs)
elapsed_seconds <- proc.time()[["elapsed"]] - start_time

summary_table <- do.call(rbind, lapply(results, `[[`, "summary"))
reference_summary <- reference_cache$alpha005_summary
reference_summary <- reference_summary[
  reference_summary$method %in% c("FASH-IWP1-Raw", "FASH-IWP1-BF"),
  c(
    "se_regime", "method", "n_true_alternatives", "n_discoveries",
    "false_discoveries", "true_positives", "power", "realized_fdp"
  ),
  drop = FALSE
]
summary_table <- rbind(reference_summary, summary_table)
summary_table <- summary_table[
  order(summary_table$se_regime, summary_table$method),
  ,
  drop = FALSE
]
rownames(summary_table) <- NULL

prior_summary <- do.call(rbind, lapply(results, `[[`, "prior_summary"))
reference_pi0 <- reference_cache$pi0_summary
reference_pi0$method <- paste0(
  reference_pi0$method,
  ifelse(reference_pi0$fit == "Raw", "-Raw", "-BF")
)
reference_pi0 <- reference_pi0[
  reference_pi0$method %in% c("FASH-IWP1-Raw", "FASH-IWP1-BF"),
  c("se_regime", "method", "estimated_pi0"),
  drop = FALSE
]
reference_pi0$selected_endpoint_slope_sd <- NA_real_
prior_summary <- rbind(reference_pi0, prior_summary)

matched_reproduction <- merge(
  summary_table[
    summary_table$method == "FASH-linear-mixture-BF",
    c("se_regime", "power")
  ],
  reference_cache$alpha005_summary[
    reference_cache$alpha005_summary$method == "FASH-linear-BF",
    c("se_regime", "power")
  ],
  by = "se_regime",
  suffixes = c("_recomputed", "_retained"),
  sort = TRUE
)

validation <- data.frame(
  check = c(
    "reference cache validated",
    "matched mixture reproduces retained R0 power",
    "all summaries finite",
    "one-job thread cap"
  ),
  passed = c(
    TRUE,
    nrow(matched_reproduction) == length(se_specs) &&
      max(abs(
        matched_reproduction$power_recomputed -
          matched_reproduction$power_retained
      )) <= 1e-12,
    all(is.finite(summary_table$power)) &&
      all(is.finite(summary_table$realized_fdp)) &&
      all(is.finite(prior_summary$estimated_pi0)),
    Sys.getenv("VECLIB_MAXIMUM_THREADS") == "1"
  ),
  stringsAsFactors = FALSE
)
if (any(!validation$passed)) {
  print(validation)
  stop("The profiled-slab R0 diagnostic failed validation.")
}

output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r0_linear_iwp_power_mechanism_profiled_slab_seed12345"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

analysis_cache <- list(
  configuration = list(
    scope = "internal mechanism diagnostic",
    seed = seed,
    J = J,
    time_grid = time_grid,
    class_probs = class_probs,
    se_specs = se_specs,
    profiled_sigma_beta_grid = sigma_beta_grid,
    matched_grid = matched_grid,
    alpha = alpha,
    reference_cache_path = normalizePath(reference_cache_path),
    retained_real_data_linear_object_class = "profiled_linear_fash",
    current_matched_linear_object_class = "linear_mixture_fash"
  ),
  summary = summary_table,
  prior_summary = prior_summary,
  matched_prior = do.call(rbind, lapply(results, `[[`, "matched_prior")),
  sigma_profile = do.call(rbind, lapply(results, `[[`, "sigma_profile")),
  validation = validation,
  elapsed_seconds = elapsed_seconds
)
saveRDS(
  analysis_cache,
  file.path(output_directory, "analysis_cache.rds")
)
utils::write.csv(
  summary_table,
  file.path(output_directory, "alpha005_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  prior_summary,
  file.path(output_directory, "prior_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  analysis_cache$sigma_profile,
  file.path(output_directory, "sigma_profile.csv"),
  row.names = FALSE
)
utils::write.csv(
  validation,
  file.path(output_directory, "validation.csv"),
  row.names = FALSE
)

print(summary_table)
print(prior_summary)
print(validation)
message("Elapsed seconds: ", round(elapsed_seconds, 3))
message("Saved output to: ", normalizePath(output_directory))
