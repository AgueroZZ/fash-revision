# Load, validate, and format cached results for the R1 workflowr report.

mc_output_dir <- file.path(
  "output",
  "revision_simulations",
  "mc",
  "r1_random_bspline_main_effect_profile_sigma_pilot5"
)
summary_dir <- file.path(mc_output_dir, "summary")
raw_path <- file.path(mc_output_dir, "full_fits", "seed_12345.rds")

out <- readRDS(raw_path)
configuration <- readRDS(file.path(mc_output_dir, "configuration.rds"))
fash_alpha <- read.csv(
  file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha_curve.csv"),
  stringsAsFactors = FALSE
)
fash_alpha_005 <- read.csv(
  file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha005_summary.csv"),
  stringsAsFactors = FALSE
)
direct_alpha <- read.csv(
  file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_mc_alpha_curve.csv"),
  stringsAsFactors = FALSE
)
direct_alpha_005 <- read.csv(
  file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_mc_alpha005_summary.csv"),
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

expected_class_probs <- c(
  dynamic_bspline = 0.20,
  constant = 0.40,
  zero = 0.40
)
expected_linear_sigma_grid <- exp(seq(log(0.05), log(5), length.out = 25))
if (!identical(
      configuration$scenario,
      "genotype_random_bspline_main_effect_dynamic_eqtl"
    ) ||
    !isTRUE(all.equal(configuration$J, 1000L)) ||
    !isTRUE(all.equal(configuration$n_donors, 19L)) ||
    !isTRUE(all.equal(configuration$n_covariates, 5L)) ||
    !isTRUE(all.equal(configuration$time_grid, 0:15)) ||
    !isTRUE(all.equal(configuration$class_probs, expected_class_probs)) ||
    !isTRUE(all.equal(configuration$dynamic_amplitude, 2)) ||
    !isTRUE(all.equal(configuration$bspline_df, 6)) ||
    !isTRUE(all.equal(configuration$bspline_coefficient_sd, 1)) ||
    !isTRUE(all.equal(configuration$expression_noise_sd, 1)) ||
    !isTRUE(all.equal(configuration$dynamic_main_effect_sd, 1)) ||
    !isTRUE(all.equal(configuration$num_basis, 20L)) ||
    !isTRUE(all.equal(configuration$true_pi0, 0.8)) ||
    !isTRUE(all.equal(configuration$efdr_permutations, 100L)) ||
    !identical(configuration$linear_sigma_estimation, "profile_grid") ||
    !isTRUE(all.equal(
      configuration$linear_sigma_grid,
      expected_linear_sigma_grid,
      tolerance = 0
    )) ||
    !isTRUE(all.equal(configuration$full_fit_seed, 12345L)) ||
    length(configuration$seed_list) != 5L) {
  stop("The random B-spline Monte Carlo cache does not match the intended setting.")
}
if (!isTRUE(all.equal(out$settings$n_variants, 1000L)) ||
    !isTRUE(all.equal(out$settings$n_donors, 19L)) ||
    !isTRUE(all.equal(out$settings$n_covariates, 5L)) ||
    !isTRUE(all.equal(out$settings$class_probs, expected_class_probs)) ||
    !isTRUE(all.equal(out$settings$maf_range, c(0.1, 0.5))) ||
    !isTRUE(all.equal(out$settings$covariate_effect_sd, 0.5)) ||
    !isTRUE(all.equal(out$settings$intercept_sd, 0)) ||
    !isTRUE(all.equal(out$settings$dynamic_main_effect_sd, 1)) ||
    !isTRUE(out$settings$estimate_sigma) ||
    !isTRUE(all.equal(
      out$settings$sigma_beta_grid,
      expected_linear_sigma_grid,
      tolerance = 0
    )) ||
    !isTRUE(out$settings$apply_t_se_correction)) {
  stop("The example simulation object does not match the Monte Carlo setting.")
}
if (!all(fash_methods %in% fash_alpha$method) ||
    !all(direct_methods %in% direct_alpha$method) ||
    any(fash_alpha_005$n_replications != length(configuration$seed_list)) ||
    any(direct_alpha_005$n_replications != length(configuration$seed_list))) {
  stop("The Monte Carlo summaries are incomplete.")
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
  summary_table$method <- factor(summary_table$method, levels = method_order)
  summary_table <- summary_table[
    order(summary_table$method),
    ,
    drop = FALSE
  ]
  data.frame(
    Method = as.character(summary_table$method),
    `Mean discoveries` = format_decimal(summary_table$mean_discoveries, 1),
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

get_mc_metric <- function(summary_table, method, metric) {
  row <- summary_table[summary_table$method == method, , drop = FALSE]
  if (nrow(row) != 1L || !metric %in% names(row)) {
    stop("Could not extract a unique Monte Carlo metric.")
  }
  row[[metric]]
}

fash_table <- format_mc_table(fash_alpha_005, fash_methods)
direct_table <- format_mc_table(direct_alpha_005, direct_methods)
pi0_table <- format_pi0_table(pi0_summary, linear_sigma_summary)
