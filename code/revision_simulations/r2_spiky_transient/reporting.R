# Load, validate, and format cached results for the R2 workflowr report.

mc_output_id <- paste0(
  "r2_timed_cosine_one_two_three_peak_main_effect_",
  "profile_sigma_pilot5"
)
mc_output_dir <- file.path(
  "output",
  "revision_simulations",
  "mc",
  mc_output_id
)
summary_dir <- file.path(mc_output_dir, "summary")
configuration <- readRDS(file.path(mc_output_dir, "configuration.rds"))
out <- readRDS(file.path(
  mc_output_dir,
  "full_fits",
  "seed_12345.rds"
))

fash_alpha <- read.csv(
  file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha_curve.csv"),
  stringsAsFactors = FALSE
)
fash_alpha_005 <- read.csv(
  file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha005_summary.csv"),
  stringsAsFactors = FALSE
)
direct_alpha <- read.csv(
  file.path(
    summary_dir,
    "iwp_fash_vs_direct_true_pi0_mc_alpha_curve.csv"
  ),
  stringsAsFactors = FALSE
)
direct_alpha_005 <- read.csv(
  file.path(
    summary_dir,
    "iwp_fash_vs_direct_true_pi0_mc_alpha005_summary.csv"
  ),
  stringsAsFactors = FALSE
)
peak_alpha <- read.csv(
  file.path(summary_dir, "mc_peak_alpha_curve.csv"),
  stringsAsFactors = FALSE
)
peak_alpha_005 <- read.csv(
  file.path(summary_dir, "mc_peak_alpha005_summary.csv"),
  stringsAsFactors = FALSE
)
pi0_summary <- read.csv(
  file.path(summary_dir, "mc_pi0_summary.csv"),
  stringsAsFactors = FALSE
)
linear_sigma_profiles <- read.csv(
  file.path(summary_dir, "all_replicate_linear_sigma_profiles.csv"),
  stringsAsFactors = FALSE
)
linear_sigma_summary <- read.csv(
  file.path(summary_dir, "linear_sigma_summary.csv"),
  stringsAsFactors = FALSE
)

class_probs <- c(
  dynamic_bspline = 0.20,
  constant = 0.40,
  zero = 0.40
)
shape_cell_probs <- c(
  k1__spiky__single = 1 / 3,
  `k2__spiky__same-sign` = 1 / 6,
  `k2__spiky__alternating-sign` = 1 / 6,
  `k3__spiky__same-sign` = 1 / 6,
  `k3__spiky__alternating-sign` = 1 / 6
)
expected_shape_counts <- exact_proportional_counts(
  200,
  shape_cell_probs
)
expected_linear_sigma_grid <- exp(seq(log(0.05), log(5), length.out = 25))
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

if (!identical(
      configuration$scenario,
      paste0(
        "r2_genotype_timed_cosine_one_two_three_peak_",
        "main_effect_dynamic_eqtl"
      )
    ) ||
    !isTRUE(all.equal(configuration$J, 1000L)) ||
    !isTRUE(all.equal(configuration$n_donors, 19L)) ||
    !isTRUE(all.equal(configuration$n_covariates, 5L)) ||
    !isTRUE(all.equal(configuration$time_grid, 0:15)) ||
    !isTRUE(all.equal(configuration$class_probs, class_probs)) ||
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
    !isTRUE(all.equal(configuration$efdr_permutations, 100L)) ||
    !isTRUE(all.equal(configuration$true_pi0, 0.8)) ||
    !identical(configuration$linear_sigma_estimation, "profile_grid") ||
    !isTRUE(all.equal(
      configuration$linear_sigma_grid,
      expected_linear_sigma_grid,
      tolerance = 0
    )) ||
    !isTRUE(all.equal(configuration$full_fit_seed, 12345L)) ||
    length(configuration$seed_list) != 5L) {
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
if (!isTRUE(all.equal(out$settings$maf_range, c(0.1, 0.5))) ||
    !isTRUE(all.equal(out$settings$covariate_effect_sd, 0.5)) ||
    !isTRUE(all.equal(out$settings$intercept_sd, 0)) ||
    !isTRUE(out$settings$estimate_sigma) ||
    !isTRUE(all.equal(
      out$settings$sigma_beta_grid,
      expected_linear_sigma_grid,
      tolerance = 0
    )) ||
    !isTRUE(out$settings$apply_t_se_correction)) {
  stop("The spiky-transient example object has unexpected settings.")
}

validate_simplified_sigma_profile(
  out$simplified_fit,
  require_interior = TRUE
)
validate_simplified_sigma_profile(
  out$simplified_fit_bf,
  require_interior = TRUE
)
if (!isTRUE(all.equal(
      out$simplified_fit$sigma_beta,
      out$simplified_fit_bf$sigma_beta,
      tolerance = 0
    ))) {
  stop("The BF update changed the selected linear slope scale.")
}

selected_linear_sigma <- linear_sigma_profiles[
  linear_sigma_profiles$selected,
  ,
  drop = FALSE
]
profile_counts <- table(linear_sigma_profiles$seed)
selected_counts <- table(selected_linear_sigma$seed)
profile_grids_match <- vapply(
  split(linear_sigma_profiles$sigma_beta, linear_sigma_profiles$seed),
  function(x) isTRUE(all.equal(x, expected_linear_sigma_grid, tolerance = 1e-12)),
  logical(1)
)
sigma_summary_fields <- c(
  "mean_selected_sigma", "selected_sigma_sd", "selected_sigma_mc_se",
  "selected_sigma_ci_lower", "selected_sigma_ci_upper",
  "min_selected_sigma", "max_selected_sigma"
)
if (nrow(linear_sigma_profiles) !=
      length(configuration$seed_list) * length(expected_linear_sigma_grid) ||
    !identical(sort(unique(linear_sigma_profiles$seed)),
               sort(configuration$seed_list)) ||
    any(profile_counts != length(expected_linear_sigma_grid)) ||
    nrow(selected_linear_sigma) != length(configuration$seed_list) ||
    any(selected_counts != 1L) ||
    any(selected_linear_sigma$grid_boundary) ||
    !all(profile_grids_match) ||
    nrow(linear_sigma_summary) != 1L ||
    !identical(linear_sigma_summary$estimation, "profile_grid") ||
    linear_sigma_summary$n_replications != length(configuration$seed_list) ||
    !all(vapply(
      linear_sigma_summary[, sigma_summary_fields, drop = FALSE],
      function(x) all(is.finite(x)),
      logical(1)
    )) ||
    !isTRUE(all.equal(
      linear_sigma_summary$mean_selected_sigma,
      mean(selected_linear_sigma$sigma_beta)
    )) ||
    !isTRUE(all.equal(
      linear_sigma_summary$min_selected_sigma,
      min(selected_linear_sigma$sigma_beta)
    )) ||
    !isTRUE(all.equal(
      linear_sigma_summary$max_selected_sigma,
      max(selected_linear_sigma$sigma_beta)
    ))) {
  stop("The linear-FASH slope-scale profile cache is incomplete or invalid.")
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

format_pi0_table <- function(summary_table, sigma_summary) {
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
  selected_sigma_label <- paste0(
    format_decimal(sigma_summary$mean_selected_sigma),
    " [",
    format_decimal(sigma_summary$min_selected_sigma),
    ", ",
    format_decimal(sigma_summary$max_selected_sigma),
    "]"
  )
  data.frame(
    Method = as.character(summary_table$method),
    Fit = as.character(summary_table$fit),
    `Mean estimated pi0 (95% MC CI)` = format_mc_interval(
      summary_table$mean_estimated_pi0,
      summary_table$pi0_ci_lower,
      summary_table$pi0_ci_upper
    ),
    `Selected slope SD: mean [range]` = ifelse(
      as.character(summary_table$method) == "FASH-linear",
      selected_sigma_label,
      "\u2014"
    ),
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
pi0_table <- format_pi0_table(pi0_summary, linear_sigma_summary)
peak_table <- format_peak_table(peak_alpha_005)
