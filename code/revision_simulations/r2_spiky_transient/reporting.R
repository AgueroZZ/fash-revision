# Load, validate, and format cached results for the R2 workflowr report.

expected_output_id <- paste0(
  "r2_real_genotype_one_per_gene_J6362_",
  "timed_cosine_one_two_three_peak_main_effect_",
  "linear_mixture_predstep1_penalty10_pilot5"
)
mc_output_dir <- file.path(
  "output",
  "revision_simulations",
  "mc",
  expected_output_id
)
summary_dir <- file.path(mc_output_dir, "summary")
raw_path <- file.path(
  mc_output_dir,
  "full_fits",
  "seed_12345.rds"
)
configuration_path <- file.path(mc_output_dir, "configuration.rds")
pairing_validation_path <- file.path(
  summary_dir,
  "real_genotype_pairing_validation.csv"
)

required_cache_paths <- c(
  raw_path,
  configuration_path,
  file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha_curve.csv"),
  file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha005_summary.csv"),
  file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_mc_alpha_curve.csv"),
  file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_mc_alpha005_summary.csv"),
  file.path(summary_dir, "mc_peak_alpha_curve.csv"),
  file.path(summary_dir, "mc_peak_alpha005_summary.csv"),
  file.path(summary_dir, "mc_pi0_summary.csv"),
  file.path(summary_dir, "all_replicate_linear_prior_weights.csv"),
  file.path(summary_dir, "all_replicate_linear_prior_summary.csv"),
  file.path(summary_dir, "genotype_selection_summary.csv"),
  file.path(summary_dir, "truth_maf_balance.csv"),
  pairing_validation_path
)
missing_cache_paths <- required_cache_paths[!file.exists(required_cache_paths)]
if (length(missing_cache_paths) > 0L) {
  stop(
    "The formal R2 linear-mixture cache is incomplete. Missing: ",
    paste(missing_cache_paths, collapse = ", ")
  )
}

read_summary_csv <- function(file_name) {
  read.csv(
    file.path(summary_dir, file_name),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

out <- readRDS(raw_path)
configuration <- readRDS(configuration_path)
fash_alpha <- read_summary_csv("iwp_vs_linear_fash_mc_alpha_curve.csv")
fash_alpha_005 <- read_summary_csv(
  "iwp_vs_linear_fash_mc_alpha005_summary.csv"
)
direct_alpha <- read_summary_csv(
  "iwp_fash_vs_direct_true_pi0_mc_alpha_curve.csv"
)
direct_alpha_005 <- read_summary_csv(
  "iwp_fash_vs_direct_true_pi0_mc_alpha005_summary.csv"
)
peak_alpha <- read_summary_csv("mc_peak_alpha_curve.csv")
peak_alpha_005 <- read_summary_csv("mc_peak_alpha005_summary.csv")
pi0_summary <- read_summary_csv("mc_pi0_summary.csv")
linear_prior_weights <- read_summary_csv(
  "all_replicate_linear_prior_weights.csv"
)
linear_prior_summary <- read_summary_csv(
  "all_replicate_linear_prior_summary.csv"
)
genotype_selection_summary <- read_summary_csv(
  "genotype_selection_summary.csv"
)
truth_maf_balance <- read_summary_csv("truth_maf_balance.csv")
pairing_validation <- read.csv(
  pairing_validation_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

historical_out <- out
historical_fash_alpha <- fash_alpha
historical_direct_alpha <- direct_alpha
historical_direct_alpha_005 <- direct_alpha_005
historical_peak_alpha <- peak_alpha
historical_peak_alpha_005 <- peak_alpha_005
source(
  "code/revision_simulations/r1_r2_fashr0143/reporting_overlay.R"
)
r2_reporting_overlay <- load_r1_r2_fashr0143_overlay(
  scenario_id = "r2",
  historical_out = historical_out,
  historical_configuration = configuration,
  historical_fash_alpha = historical_fash_alpha,
  historical_direct_alpha = historical_direct_alpha,
  historical_direct_alpha_005 = historical_direct_alpha_005,
  historical_peak_alpha = historical_peak_alpha,
  historical_peak_alpha_005 = historical_peak_alpha_005
)
out <- r2_reporting_overlay$out
fash_alpha <- r2_reporting_overlay$fash_alpha
fash_alpha_005 <- r2_reporting_overlay$fash_alpha_005
direct_alpha <- r2_reporting_overlay$direct_alpha
direct_alpha_005 <- r2_reporting_overlay$direct_alpha_005
peak_alpha <- r2_reporting_overlay$peak_alpha
peak_alpha_005 <- r2_reporting_overlay$peak_alpha_005
pi0_summary <- r2_reporting_overlay$pi0_summary
linear_prior_weights <- r2_reporting_overlay$linear_prior_weights
linear_prior_summary <- r2_reporting_overlay$linear_prior_summary
genotype_selection_summary <-
  r2_reporting_overlay$genotype_selection_summary
truth_maf_balance <- r2_reporting_overlay$truth_maf_balance
result_provenance_table <- r2_reporting_overlay$provenance_table
current_fash_manifest <- r2_reporting_overlay$manifest

class_probs <- c(
  dynamic_bspline = 0.20,
  constant = 0.40,
  zero = 0.40
)
expected_J <- 6362L
expected_class_counts <- exact_proportional_counts(expected_J, class_probs)
expected_true_pi0 <- unname(
  sum(expected_class_counts[c("constant", "zero")]) / expected_J
)
shape_cell_probs <- c(
  k1__spiky__single = 1 / 3,
  `k2__spiky__same-sign` = 1 / 6,
  `k2__spiky__alternating-sign` = 1 / 6,
  `k3__spiky__same-sign` = 1 / 6,
  `k3__spiky__alternating-sign` = 1 / 6
)
expected_shape_counts <- exact_proportional_counts(
  expected_class_counts[["dynamic_bspline"]],
  shape_cell_probs
)
fash_methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "FASH-linear-Raw",
  "FASH-linear-BF"
)
direct_methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "Direct-linear-LRT-eFDR-true-pi0",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
peak_methods <- c(
  "FASH-IWP1-BF",
  "Direct-linear-LRT-eFDR-true-pi0",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
expected_seed_list <- c(12345L, 22345L, 32345L, 42345L, 52345L)
expected_sd_grid <- default_revision_grid()
expected_pred_step <- 1
expected_penalty <- 10L

if (!identical(configuration$output_id, expected_output_id) ||
    !identical(
      configuration$scenario,
      paste0(
        "r2_real_genotype_one_per_gene_",
        "timed_cosine_one_two_three_peak_",
        "main_effect_dynamic_eqtl"
      )
    ) ||
    !isTRUE(all.equal(configuration$J, expected_J)) ||
    !isTRUE(all.equal(configuration$n_donors, 19L)) ||
    !isTRUE(all.equal(configuration$n_covariates, 5L)) ||
    !isTRUE(all.equal(configuration$time_grid, 0:15)) ||
    !isTRUE(all.equal(configuration$class_probs, class_probs)) ||
    !isTRUE(all.equal(configuration$expected_class_counts, expected_class_counts)) ||
    !isTRUE(all.equal(configuration$spike_counts, 1:3)) ||
    !isTRUE(all.equal(configuration$shape_cell_probs, shape_cell_probs)) ||
    !isTRUE(all.equal(
      configuration$shape_cell_counts,
      expected_shape_counts
    )) ||
    !isTRUE(all.equal(
      configuration$primary_time_groups,
      c("early", "middle", "late")
    )) ||
    !isTRUE(all.equal(
      configuration$relative_amplitude_range,
      c(0.35, 0.75)
    )) ||
    !identical(configuration$center_by_observed_mean, FALSE) ||
    !isTRUE(all.equal(configuration$width_half, 1.5)) ||
    !isTRUE(all.equal(configuration$target_centered_rms, 0.9)) ||
    !isTRUE(all.equal(configuration$dynamic_main_effect_sd, 1)) ||
    !isTRUE(all.equal(configuration$dynamic_baseline_sd, 1)) ||
    !isTRUE(all.equal(configuration$expression_noise_sd, 1)) ||
    !isTRUE(all.equal(configuration$num_basis, 20L)) ||
    !isTRUE(all.equal(configuration$efdr_permutations, 100L)) ||
    !isTRUE(all.equal(configuration$true_pi0, expected_true_pi0)) ||
    !identical(configuration$linear_prior_mode, "mixture_grid") ||
    !isTRUE(all.equal(
      configuration$common_sd_grid,
      expected_sd_grid,
      tolerance = 0
    )) ||
    !isTRUE(all.equal(
      configuration$common_pred_step,
      expected_pred_step,
      tolerance = 0
    )) ||
    !identical(as.integer(configuration$common_penalty), expected_penalty) ||
    !identical(configuration$genotype_source, "paper-derived YRI DS dosage") ||
    !identical(
      configuration$genotype_selection_rule,
      "one uniformly sampled tested variant per gene"
    ) ||
    !isTRUE(all.equal(configuration$genotype_maf_min, 0.10)) ||
    !identical(length(configuration$genotype_sample_ids), 19L) ||
    !isTRUE(all.equal(configuration$full_fit_seed, 12345L)) ||
    !identical(as.integer(configuration$seed_list), expected_seed_list)) {
  stop("The spiky-transient Monte Carlo cache has unexpected settings.")
}
if (!all(fash_methods %in% fash_alpha$method) ||
    !all(direct_methods %in% direct_alpha$method) ||
    !all(peak_methods %in% peak_alpha$method) ||
    any(fash_alpha_005$n_replications != length(configuration$seed_list)) ||
    any(direct_alpha_005$n_replications != length(configuration$seed_list)) ||
    any(peak_alpha_005$n_replications != length(configuration$seed_list))) {
  stop("The spiky-transient Monte Carlo summaries are incomplete.")
}
if (!isTRUE(all.equal(out$settings$n_variants, expected_J)) ||
    !isTRUE(all.equal(out$settings$n_donors, 19L)) ||
    !isTRUE(all.equal(out$settings$n_covariates, 5L)) ||
    !isTRUE(all.equal(out$settings$class_probs, class_probs)) ||
    !isTRUE(all.equal(out$settings$covariate_effect_sd, 0.5)) ||
    !isTRUE(all.equal(out$settings$intercept_sd, 0)) ||
    !identical(out$settings$linear_prior_mode, "mixture_grid") ||
    !isTRUE(all.equal(
      out$settings$pred_step,
      expected_pred_step,
      tolerance = 0
    )) ||
    !identical(as.integer(out$settings$penalty), expected_penalty) ||
    !identical(out$settings$genotype_source, configuration$genotype_source) ||
    !identical(
      out$settings$genotype_selection_rule,
      configuration$genotype_selection_rule
    ) ||
    !identical(out$genotype_digest, out$settings$genotype_digest) ||
    is.null(out$genotype_digest) ||
    !nzchar(out$genotype_digest) ||
    !isTRUE(all.equal(out$settings$seed, 12345L)) ||
    !isTRUE(out$settings$apply_t_se_correction)) {
  stop("The spiky-transient example object has unexpected settings.")
}

validate_linear_mixture_fash(
  out$simplified_fit,
  expected_grid = expected_sd_grid,
  expected_pred_step = expected_pred_step,
  expected_penalty = expected_penalty
)
validate_linear_mixture_fash(
  out$simplified_fit_bf,
  expected_grid = expected_sd_grid,
  expected_pred_step = expected_pred_step,
  expected_penalty = expected_penalty
)
if (!identical(out$simplified_fit$bf_adjusted, FALSE) ||
    !isTRUE(out$simplified_fit_bf$bf_adjusted)) {
  stop("The linear-mixture raw/BF cache stages are inconsistent.")
}

for (fit_name in c("fash_iwp1_raw", "fash_iwp1_bf")) {
  iwp_fit <- out$fash_fits[[fit_name]]
  if (is.null(iwp_fit) ||
      !isTRUE(all.equal(iwp_fit$psd_grid, expected_sd_grid, tolerance = 0)) ||
      !isTRUE(all.equal(
        iwp_fit$settings$pred_step,
        expected_pred_step,
        tolerance = 0
      )) ||
      !identical(as.integer(iwp_fit$settings$penalty), expected_penalty)) {
    stop("The cached IWP and linear fits do not share grid, pred_step, and penalty.")
  }
}

required_pairing_columns <- c(
  "check", "seed", "mismatched_items",
  "max_absolute_numeric_difference", "status", "detail"
)
required_pairing_checks <- c(
  "R1_R2_shared_configuration",
  "R1_R2_R3_shared_genotype_configuration",
  "shared_genotype_cache_fingerprint",
  "per_seed_genotype_digest",
  "per_seed_selected_pair_keys",
  "per_seed_truth_class_maf_balance"
)
if (!all(required_pairing_columns %in% names(pairing_validation)) ||
    !all(required_pairing_checks %in% pairing_validation$check) ||
    anyNA(pairing_validation$status) ||
    any(pairing_validation$status != "PASS") ||
    anyNA(pairing_validation$mismatched_items) ||
    any(pairing_validation$mismatched_items != 0) ||
    anyNA(pairing_validation$max_absolute_numeric_difference) ||
    any(pairing_validation$max_absolute_numeric_difference != 0)) {
  stop("The required exact R1/R2/R3 real-genotype pairing validation did not pass.")
}

if (nrow(genotype_selection_summary) != length(expected_seed_list) ||
    !setequal(genotype_selection_summary$seed, expected_seed_list) ||
    any(genotype_selection_summary$genes != expected_J) ||
    any(genotype_selection_summary$selected_pairs != expected_J) ||
    any(genotype_selection_summary$maf_min < 0.10) ||
    any(genotype_selection_summary$maf_max > 0.5) ||
    nrow(truth_maf_balance) !=
      length(expected_seed_list) * length(class_probs) ||
    !setequal(truth_maf_balance$seed, expected_seed_list) ||
    !setequal(truth_maf_balance$effect_class, names(class_probs)) ||
    any(abs(truth_maf_balance$standardized_mean_difference) > 0.01)) {
  stop("The real-genotype selection or truth-class MAF diagnostics are invalid.")
}

prior_weight_columns <- c(
  "seed", "fit", "predstep_sd", "prior_weight", "is_null", "active"
)
prior_summary_columns <- c(
  "seed", "fit", "estimated_pi0", "active_nonnull_components",
  "alternative_rms_predstep_sd"
)
expected_fit_labels <- c("Raw", "BF-corrected")
make_seed_fit_key <- function(seed, fit) paste(seed, fit, sep = "\r")
expected_seed_fit_keys <- as.vector(outer(
  expected_seed_list,
  expected_fit_labels,
  make_seed_fit_key
))
if (!all(prior_weight_columns %in% names(linear_prior_weights)) ||
    nrow(linear_prior_weights) !=
      length(expected_seed_fit_keys) * length(expected_sd_grid) ||
    !setequal(unique(linear_prior_weights$seed), expected_seed_list) ||
    !setequal(unique(linear_prior_weights$fit), expected_fit_labels)) {
  stop("The full-grid linear-mixture prior-weight cache is incomplete.")
}

linear_prior_weights$seed_fit_key <- make_seed_fit_key(
  linear_prior_weights$seed,
  linear_prior_weights$fit
)
prior_weight_groups <- split(
  linear_prior_weights,
  linear_prior_weights$seed_fit_key
)
weights_valid <- length(prior_weight_groups) == length(expected_seed_fit_keys) &&
  setequal(names(prior_weight_groups), expected_seed_fit_keys) &&
  all(vapply(prior_weight_groups, function(rows) {
    nrow(rows) == length(expected_sd_grid) &&
      isTRUE(all.equal(
        rows$predstep_sd,
        expected_sd_grid,
        tolerance = 1e-12
      )) &&
      isTRUE(all.equal(rows$is_null, expected_sd_grid == 0, tolerance = 0)) &&
      all(is.finite(rows$prior_weight)) &&
      all(rows$prior_weight >= 0) &&
      abs(sum(rows$prior_weight) - 1) < 1e-8 &&
      identical(rows$active, rows$prior_weight > 0)
  }, logical(1)))
if (!weights_valid) {
  stop("The full-grid linear-mixture prior weights are invalid.")
}

if (!all(prior_summary_columns %in% names(linear_prior_summary))) {
  stop("The linear-mixture prior summary is missing required columns.")
}
linear_prior_summary$seed_fit_key <- make_seed_fit_key(
  linear_prior_summary$seed,
  linear_prior_summary$fit
)
if (nrow(linear_prior_summary) != length(expected_seed_fit_keys) ||
    anyDuplicated(linear_prior_summary$seed_fit_key) ||
    !setequal(linear_prior_summary$seed_fit_key, expected_seed_fit_keys) ||
    any(!is.finite(linear_prior_summary$estimated_pi0)) ||
    any(linear_prior_summary$estimated_pi0 < 0) ||
    any(linear_prior_summary$estimated_pi0 > 1) ||
    any(!is.finite(linear_prior_summary$active_nonnull_components)) ||
    any(linear_prior_summary$active_nonnull_components < 0) ||
    any(linear_prior_summary$active_nonnull_components !=
      floor(linear_prior_summary$active_nonnull_components))) {
  stop("The linear-mixture prior summary is incomplete or invalid.")
}

prior_summary_matches_weights <- vapply(
  expected_seed_fit_keys,
  function(key) {
    weights <- prior_weight_groups[[key]]
    summary_row <- linear_prior_summary[
      linear_prior_summary$seed_fit_key == key,
      ,
      drop = FALSE
    ]
    alternative <- !weights$is_null
    alternative_mass <- sum(weights$prior_weight[alternative])
    expected_rms <- if (alternative_mass > 0) {
      sqrt(sum(
        weights$prior_weight[alternative] / alternative_mass *
          weights$predstep_sd[alternative]^2
      ))
    } else {
      NA_real_
    }
    rms_matches <- if (is.na(expected_rms)) {
      is.na(summary_row$alternative_rms_predstep_sd)
    } else {
      isTRUE(all.equal(
        summary_row$alternative_rms_predstep_sd,
        expected_rms,
        tolerance = 1e-10
      ))
    }
    isTRUE(all.equal(
      summary_row$estimated_pi0,
      weights$prior_weight[weights$is_null],
      tolerance = 1e-10
    )) &&
      summary_row$active_nonnull_components == sum(weights$active[alternative]) &&
      rms_matches
  },
  logical(1)
)
if (!all(prior_summary_matches_weights)) {
  stop("The linear-mixture prior summary does not match the full-grid weights.")
}

expected_pi0_keys <- as.vector(outer(
  c("FASH-IWP1", "FASH-linear"),
  expected_fit_labels,
  paste,
  sep = "\r"
))
observed_pi0_keys <- paste(pi0_summary$method, pi0_summary$fit, sep = "\r")
if (nrow(pi0_summary) != length(expected_pi0_keys) ||
    anyDuplicated(observed_pi0_keys) ||
    !setequal(observed_pi0_keys, expected_pi0_keys) ||
    any(pi0_summary$n_replications != length(expected_seed_list)) ||
    any(!is.finite(pi0_summary$mean_estimated_pi0))) {
  stop("The Monte Carlo pi0 summary is incomplete or invalid.")
}

dynamic <- out$unit_info$effect_class == "dynamic_bspline"
observed_shape_counts <- table(out$unit_info$cell_id[dynamic])
if (!identical(
      as.integer(observed_shape_counts[names(expected_shape_counts)]),
      as.integer(expected_shape_counts)
    ) ||
    any(!is.finite(out$unit_info$genetic_main_effect[dynamic])) ||
    !all(vapply(
      out$unit_info$peak_centers[dynamic],
      length,
      integer(1)
    ) == out$unit_info$spike_count[dynamic])) {
  stop("The cached spiky truth does not match its recorded allocation.")
}

format_decimal <- function(x, digits = 3) {
  formatC(x, format = "f", digits = digits)
}

format_mc_interval <- function(mean, lower, upper, digits = 3) {
  paste0(
    format_decimal(mean, digits),
    " [",
    format_decimal(lower, digits),
    ", ",
    format_decimal(upper, digits),
    "]"
  )
}

format_mean_range <- function(values,
                              mean_digits = 3,
                              range_digits = mean_digits) {
  values <- as.numeric(values)
  finite_values <- values[is.finite(values)]
  if (length(finite_values) == 0L) {
    return("\u2014")
  }
  label <- paste0(
    format_decimal(mean(finite_values), mean_digits),
    " [",
    format_decimal(min(finite_values), range_digits),
    ", ",
    format_decimal(max(finite_values), range_digits),
    "]"
  )
  if (length(finite_values) < length(values)) {
    label <- paste0(
      label,
      "; ",
      length(finite_values),
      "/",
      length(values),
      " defined"
    )
  }
  label
}

format_mc_table <- function(summary_table, method_order) {
  summary_table$method <- factor(
    summary_table$method,
    levels = method_order
  )
  summary_table <- summary_table[
    order(summary_table$method),
    ,
    drop = FALSE
  ]
  data.frame(
    Method = as.character(summary_table$method),
    `Mean discoveries` = format_decimal(
      summary_table$mean_discoveries,
      1
    ),
    `Power (95% MC CI)` = format_mc_interval(
      summary_table$mean_power,
      summary_table$power_ci_lower,
      summary_table$power_ci_upper
    ),
    `Empirical FDR: E[FDP] (95% MC CI)` = format_mc_interval(
      summary_table$mean_fdr,
      summary_table$fdr_ci_lower,
      summary_table$fdr_ci_upper
    ),
    check.names = FALSE
  )
}

format_pi0_table <- function(summary_table, mixture_summary) {
  summary_table$method <- factor(
    summary_table$method,
    levels = c("FASH-IWP1", "FASH-linear")
  )
  summary_table$fit <- factor(
    summary_table$fit,
    levels = c("Raw", "BF-corrected")
  )
  summary_table <- summary_table[
    order(summary_table$method, summary_table$fit),
    ,
    drop = FALSE
  ]
  diagnostic_label <- function(method,
                               fit,
                               column,
                               mean_digits,
                               range_digits) {
    if (method != "FASH-linear") {
      return("\u2014")
    }
    rows <- mixture_summary[mixture_summary$fit == fit, , drop = FALSE]
    if (nrow(rows) != length(expected_seed_list)) {
      stop("Could not match the linear-mixture diagnostic rows.")
    }
    format_mean_range(
      rows[[column]],
      mean_digits = mean_digits,
      range_digits = range_digits
    )
  }
  active_component_label <- vapply(
    seq_len(nrow(summary_table)),
    function(index) diagnostic_label(
      method = as.character(summary_table$method[index]),
      fit = as.character(summary_table$fit[index]),
      column = "active_nonnull_components",
      mean_digits = 1,
      range_digits = 0
    ),
    character(1)
  )
  alternative_rms_label <- vapply(
    seq_len(nrow(summary_table)),
    function(index) diagnostic_label(
      method = as.character(summary_table$method[index]),
      fit = as.character(summary_table$fit[index]),
      column = "alternative_rms_predstep_sd",
      mean_digits = 3,
      range_digits = 3
    ),
    character(1)
  )
  data.frame(
    Method = as.character(summary_table$method),
    Fit = as.character(summary_table$fit),
    `Mean estimated pi0 (95% MC CI)` = format_mc_interval(
      summary_table$mean_estimated_pi0,
      summary_table$pi0_ci_lower,
      summary_table$pi0_ci_upper
    ),
    `Active non-null components: mean [range]` = active_component_label,
    `Alternative RMS pred-step SD: mean [range]` = alternative_rms_label,
    check.names = FALSE
  )
}

format_peak_table <- function(summary_table) {
  table <- summary_table[
    summary_table$method %in% peak_methods,
    ,
    drop = FALSE
  ]
  table$method <- factor(table$method, levels = peak_methods)
  table <- table[order(table$spike_count, table$method), , drop = FALSE]
  data.frame(
    Peaks = table$spike_count,
    Method = as.character(table$method),
    `Dynamic eQTLs per replicate` = table$n_dynamic,
    `Power (95% MC CI)` = format_mc_interval(
      table$mean_power,
      table$power_ci_lower,
      table$power_ci_upper
    ),
    check.names = FALSE
  )
}

metric_at <- function(summary_table, method, metric) {
  row <- summary_table[
    summary_table$method == method,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L || !metric %in% names(row)) {
    stop("Could not extract a unique Monte Carlo metric.")
  }
  row[[metric]]
}

cell_labels <- c(
  k1__spiky__single = "1 peak",
  `k2__spiky__same-sign` = "2 peaks; same sign",
  `k2__spiky__alternating-sign` = "2 peaks; alternating sign",
  `k3__spiky__same-sign` = "3 peaks; same sign",
  `k3__spiky__alternating-sign` = "3 peaks; alternating sign"
)
timing_counts <- table(
  out$unit_info$cell_id[dynamic],
  out$unit_info$time_group[dynamic]
)
allocation_table <- data.frame(
  Shape = unname(cell_labels[rownames(timing_counts)]),
  Early = timing_counts[, "early"],
  Middle = timing_counts[, "middle"],
  Late = timing_counts[, "late"],
  Total = rowSums(timing_counts),
  check.names = FALSE
)

fash_table <- format_mc_table(fash_alpha_005, fash_methods)
direct_table <- format_mc_table(direct_alpha_005, direct_methods)
pi0_table <- format_pi0_table(pi0_summary, linear_prior_summary)
peak_table <- format_peak_table(peak_alpha_005)

format_integer <- function(x) {
  format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

genotype_provenance_table <- data.frame(
  Artifact = c("YRI genotype VCF", "Paper gene-variant summary", "Shared genotype cache"),
  File = c(
    configuration$genotype_source_configuration$vcf_fingerprint$file_name,
    configuration$genotype_source_configuration$pair_summary_fingerprint$file_name,
    configuration$genotype_cache_fingerprint$file_name
  ),
  `Size (bytes)` = format_integer(c(
    configuration$genotype_source_configuration$vcf_fingerprint$size_bytes,
    configuration$genotype_source_configuration$pair_summary_fingerprint$size_bytes,
    configuration$genotype_cache_fingerprint$size_bytes
  )),
  MD5 = c(
    configuration$genotype_source_configuration$vcf_fingerprint$md5,
    configuration$genotype_source_configuration$pair_summary_fingerprint$md5,
    configuration$genotype_cache_fingerprint$md5
  ),
  check.names = FALSE
)

genotype_sampling_table <- data.frame(
  Seed = genotype_selection_summary$seed,
  Genes = format_integer(genotype_selection_summary$genes),
  `Unique rsIDs` = format_integer(genotype_selection_summary$unique_variant_ids),
  `Repeated cross-gene assignments` =
    format_integer(genotype_selection_summary$repeated_cross_gene_assignments),
  `MAF range` = paste0(
    format_decimal(genotype_selection_summary$maf_min, 3),
    "-",
    format_decimal(genotype_selection_summary$maf_max, 3)
  ),
  `Median MAF` = format_decimal(genotype_selection_summary$maf_median, 3),
  check.names = FALSE
)

truth_maf_balance_table <- do.call(rbind, lapply(
  names(class_probs),
  function(effect_class) {
    rows <- truth_maf_balance[
      truth_maf_balance$effect_class == effect_class,
      ,
      drop = FALSE
    ]
    data.frame(
      `Truth class` = c(
        dynamic_bspline = "Dynamic",
        constant = "Constant",
        zero = "Zero"
      )[[effect_class]],
      `Units per seed` = format_integer(unique(rows$n)),
      `Mean MAF across seeds` = format_decimal(mean(rows$maf_mean), 3),
      `Maximum absolute MAF SMD` = format_decimal(
        max(abs(rows$standardized_mean_difference)),
        4
      ),
      check.names = FALSE
    )
  }
))
rownames(truth_maf_balance_table) <- NULL

dynamic_count <- unname(expected_class_counts[["dynamic_bspline"]])
constant_count <- unname(expected_class_counts[["constant"]])
zero_count <- unname(expected_class_counts[["zero"]])
true_pi0 <- expected_true_pi0
