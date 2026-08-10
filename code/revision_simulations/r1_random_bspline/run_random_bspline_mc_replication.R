#!/usr/bin/env Rscript

# Run compact multi-seed replications for the genotype-level random B-spline
# revision simulation. The default output is a five-seed pilot; each replicate
# stores metrics only, while the existing single-seed cache retains full data.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find code/revision_simulations/shared/simulation_functions.R from the current working directory.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0 || hit[1] == length(args)) {
    return(default)
  }
  args[hit[1] + 1]
}

as_flag <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

parse_seed_list <- function(x) {
  pieces <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  if (length(pieces) == 0 || any(!nzchar(pieces))) {
    stop("--seed-list must contain one or more comma-separated integer seeds.")
  }
  seeds <- suppressWarnings(as.integer(pieces))
  if (any(is.na(seeds)) || anyDuplicated(seeds)) {
    stop("--seed-list must contain unique integer seeds.")
  }
  seeds
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

J <- as.integer(get_arg("--J", "1000"))
n_donors <- as.integer(get_arg("--n-donors", "19"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
expression_noise_sd <- as.numeric(get_arg("--noise-sd", "1"))
dynamic_main_effect_sd <- as.numeric(get_arg("--dynamic-main-effect-sd", "1"))
num_basis <- as.integer(get_arg("--num-basis", "20"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
efdr_permutations <- as.integer(get_arg("--efdr-permutations", "100"))
seed_list <- parse_seed_list(get_arg(
  "--seed-list",
  "12345,22345,32345,42345,52345"
))
output_id <- get_arg(
  "--output-id",
  "r1_random_bspline_main_effect_profile_sigma_pilot5"
)
overwrite <- as_flag(get_arg("--overwrite", "false"))

if (J < 10 || n_donors < n_covariates + 3 || n_covariates < 0 ||
    !is.finite(expression_noise_sd) || expression_noise_sd <= 0 ||
    !is.finite(dynamic_main_effect_sd) || dynamic_main_effect_sd <= 0 ||
    num_basis < 2 || num_cores < 1 || efdr_permutations < 1 || !nzchar(output_id)) {
  stop("Invalid simulation arguments.")
}

time_grid <- make_time_grid()
class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
scenario <- "genotype_random_bspline_main_effect_dynamic_eqtl"
true_pi0_expected <- unname(class_probs["constant"] + class_probs["zero"])
linear_sigma_estimation <- "profile_grid"
linear_sigma_grid <- exp(seq(log(0.05), log(5), length.out = 25))
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
all_methods <- unique(c(fash_methods, direct_methods))

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "mc",
  output_id
)
replicate_dir <- file.path(output_dir, "replicates")
summary_dir <- file.path(output_dir, "summary")
figure_dir <- file.path(output_dir, "figures")
full_fit_dir <- file.path(output_dir, "full_fits")
invisible(lapply(c(
  output_dir,
  replicate_dir,
  summary_dir,
  figure_dir,
  full_fit_dir
), dir.create,
  recursive = TRUE, showWarnings = FALSE
))

configuration <- list(
  output_id = output_id,
  scenario = scenario,
  J = J,
  n_donors = n_donors,
  time_grid = time_grid,
  n_covariates = n_covariates,
  expression_noise_sd = expression_noise_sd,
  dynamic_main_effect_sd = dynamic_main_effect_sd,
  num_basis = num_basis,
  class_probs = class_probs,
  dynamic_amplitude = 2,
  bspline_df = 6,
  bspline_coefficient_sd = 1,
  efdr_permutations = efdr_permutations,
  true_pi0 = true_pi0_expected,
  linear_sigma_estimation = linear_sigma_estimation,
  linear_sigma_grid = linear_sigma_grid,
  full_fit_seed = seed_list[1],
  seed_list = seed_list
)
configuration_path <- file.path(output_dir, "configuration.rds")
if (file.exists(configuration_path) && !overwrite) {
  cached_configuration <- readRDS(configuration_path)
  if (!isTRUE(all.equal(cached_configuration, configuration))) {
    stop(
      "The existing output id has different settings. Choose a new --output-id or use --overwrite true."
    )
  }
} else {
  saveRDS(configuration, configuration_path)
}

make_replicate <- function(seed) {
  out <- run_genotype_level_bspline_eqtl_simulation(
    n_donors = n_donors,
    n_variants = J,
    time_grid = time_grid,
    n_covariates = n_covariates,
    class_probs = class_probs,
    expression_noise_sd = expression_noise_sd,
    dynamic_main_effect_sd = dynamic_main_effect_sd,
    scenario = scenario,
    alpha = 0.05,
    seed = seed,
    estimate_sigma = TRUE,
    sigma_beta_grid = linear_sigma_grid,
    num_cores = num_cores,
    num_basis = num_basis,
    output_dir = output_dir,
    save_outputs = FALSE,
    verbose = FALSE
  )
  true_pi0 <- dynamic_null_proportion(out$unit_info, target = "dynamic")
  if (!isTRUE(all.equal(true_pi0, true_pi0_expected))) {
    stop("The simulated dynamic-null proportion does not match the configured mixture.")
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
  linear_sigma_profile <- out$simplified_fit$sigma_profile
  linear_sigma_profile$seed <- seed
  linear_sigma_profile <- linear_sigma_profile[, c(
    "seed", "sigma_beta", "estimated_pi0", "loglik", "selected",
    "grid_boundary"
  )]
  out <- add_direct_interaction_efdr_results_to_genotype_output(
    out,
    n_permutations = efdr_permutations,
    alpha = 0.05,
    seed = seed + 10000L,
    lambda = 0.5,
    pi0_method = "conservative",
    true_pi0 = true_pi0,
    include_true_pi0 = TRUE,
    permute_covariates_with_expression = TRUE,
    num_cores = num_cores,
    overwrite = TRUE,
    verbose = FALSE
  )
  available_methods <- unique(out$result_table$method)
  missing_methods <- setdiff(all_methods, available_methods)
  if (length(missing_methods) > 0) {
    stop("Missing reviewer-facing methods: ", paste(missing_methods, collapse = ", "))
  }
  alpha_curve <- out$alpha_curve[out$alpha_curve$method %in% all_methods, ]
  alpha_curve$seed <- seed
  alpha_005 <- alpha_curve[abs(alpha_curve$alpha - 0.05) < 1e-12, ]
  pi0 <- data.frame(
    seed = seed,
    method = c("FASH-IWP1", "FASH-IWP1", "FASH-linear", "FASH-linear"),
    fit = c("Raw", "BF-corrected", "Raw", "BF-corrected"),
    estimated_pi0 = c(
      constant_component_prior_weight(out$fash_fits$fash_iwp1_raw),
      constant_component_prior_weight(out$fash_fits$fash_iwp1_bf),
      constant_component_prior_weight(out$simplified_fit),
      constant_component_prior_weight(out$simplified_fit_bf)
    ),
    stringsAsFactors = FALSE
  )
  if (identical(seed, seed_list[1])) {
    saveRDS(
      out,
      file.path(full_fit_dir, paste0("seed_", seed, ".rds"))
    )
  }
  list(
    configuration = configuration,
    seed = seed,
    permutation_seed = seed + 10000L,
    true_pi0 = true_pi0,
    alpha_curve = alpha_curve,
    alpha_005 = alpha_005,
    pi0 = pi0,
    linear_sigma_profile = linear_sigma_profile
  )
}

validate_replicate <- function(replicate, seed) {
  required_fields <- c(
    "configuration", "seed", "permutation_seed", "true_pi0",
    "alpha_curve", "alpha_005", "pi0", "linear_sigma_profile"
  )
  if (!all(required_fields %in% names(replicate)) ||
      !identical(replicate$seed, seed) ||
      !isTRUE(all.equal(replicate$configuration, configuration)) ||
      !isTRUE(all.equal(replicate$true_pi0, true_pi0_expected))) {
    return(FALSE)
  }
  observed_methods <- unique(replicate$alpha_curve$method)
  profile <- replicate$linear_sigma_profile
  all(all_methods %in% observed_methods) &&
    nrow(replicate$alpha_005) == length(all_methods) &&
    nrow(replicate$pi0) == 4 &&
    is.data.frame(profile) &&
    nrow(profile) == length(linear_sigma_grid) &&
    identical(profile$seed, rep(seed, length(linear_sigma_grid))) &&
    isTRUE(all.equal(profile$sigma_beta, linear_sigma_grid, tolerance = 0)) &&
    all(is.finite(profile$estimated_pi0)) &&
    all(is.finite(profile$loglik)) &&
    sum(profile$selected) == 1L &&
    !any(profile$selected & profile$grid_boundary)
}

replicates <- lapply(seed_list, function(seed) {
  replicate_path <- file.path(replicate_dir, paste0("seed_", seed, ".rds"))
  if (file.exists(replicate_path) && !overwrite) {
    cached <- readRDS(replicate_path)
    if (validate_replicate(cached, seed)) {
      message("Reusing compact replicate cache: ", replicate_path)
      return(cached)
    }
    stop("Cached replicate does not match the requested settings: ", replicate_path)
  }
  message("Running random B-spline replicate with seed ", seed, ".")
  replicate <- make_replicate(seed)
  saveRDS(replicate, replicate_path)
  replicate
})

all_alpha <- do.call(rbind, lapply(replicates, `[[`, "alpha_curve"))
all_alpha_005 <- do.call(rbind, lapply(replicates, `[[`, "alpha_005"))
all_pi0 <- do.call(rbind, lapply(replicates, `[[`, "pi0"))
all_linear_sigma_profiles <- do.call(
  rbind,
  lapply(replicates, `[[`, "linear_sigma_profile")
)
selected_linear_sigma <- all_linear_sigma_profiles[
  all_linear_sigma_profiles$selected,
  ,
  drop = FALSE
]
if (nrow(selected_linear_sigma) != length(seed_list) ||
    any(selected_linear_sigma$grid_boundary) ||
    !identical(sort(selected_linear_sigma$seed), sort(seed_list))) {
  stop("The selected linear slope scales are incomplete or on a grid boundary.")
}
linear_sigma_values <- summarize_mc_values(selected_linear_sigma$sigma_beta)
linear_sigma_summary <- data.frame(
  estimation = linear_sigma_estimation,
  n_replications = length(seed_list),
  mean_selected_sigma = linear_sigma_values[["mean"]],
  selected_sigma_sd = linear_sigma_values[["sd"]],
  selected_sigma_mc_se = linear_sigma_values[["se"]],
  selected_sigma_ci_lower = pmax(0, linear_sigma_values[["lower"]]),
  selected_sigma_ci_upper = linear_sigma_values[["upper"]],
  min_selected_sigma = min(selected_linear_sigma$sigma_beta),
  max_selected_sigma = max(selected_linear_sigma$sigma_beta),
  stringsAsFactors = FALSE
)
mc_alpha <- summarize_mc_alpha_curves(all_alpha)
mc_alpha_005 <- mc_alpha[abs(mc_alpha$alpha - 0.05) < 1e-12, ]
mc_pi0 <- summarize_mc_pi0(all_pi0)

fash_all_alpha <- all_alpha[all_alpha$method %in% fash_methods, ]
direct_all_alpha <- all_alpha[all_alpha$method %in% direct_methods, ]
fash_mc_alpha <- mc_alpha[mc_alpha$method %in% fash_methods, ]
direct_mc_alpha <- mc_alpha[mc_alpha$method %in% direct_methods, ]
fash_mc_alpha_005 <- mc_alpha_005[mc_alpha_005$method %in% fash_methods, ]
direct_mc_alpha_005 <- mc_alpha_005[mc_alpha_005$method %in% direct_methods, ]

write_csv(all_alpha, file.path(summary_dir, "all_replicate_alpha_curves.csv"))
write_csv(all_alpha_005, file.path(summary_dir, "all_replicate_alpha005.csv"))
write_csv(all_pi0, file.path(summary_dir, "all_replicate_pi0.csv"))
write_csv(
  all_linear_sigma_profiles,
  file.path(summary_dir, "all_replicate_linear_sigma_profiles.csv")
)
write_csv(
  linear_sigma_summary,
  file.path(summary_dir, "linear_sigma_summary.csv")
)
write_csv(mc_alpha, file.path(summary_dir, "mc_alpha_curve.csv"))
write_csv(mc_alpha_005, file.path(summary_dir, "mc_alpha005_summary.csv"))
write_csv(mc_pi0, file.path(summary_dir, "mc_pi0_summary.csv"))
write_csv(fash_all_alpha, file.path(summary_dir, "iwp_vs_linear_fash_replicate_alpha_curves.csv"))
write_csv(fash_mc_alpha, file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha_curve.csv"))
write_csv(fash_mc_alpha_005, file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha005_summary.csv"))
write_csv(direct_all_alpha, file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_replicate_alpha_curves.csv"))
write_csv(direct_mc_alpha, file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_mc_alpha_curve.csv"))
write_csv(direct_mc_alpha_005, file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_mc_alpha005_summary.csv"))

plot_subtitle <- paste0(
  length(seed_list), " independent seeds; N = ", n_donors,
  ", J = ", J, ", five covariates"
)
plot_mc_alpha_curves(
  fash_mc_alpha,
  metric = "power",
  file = file.path(figure_dir, "iwp_vs_linear_fash_mc_power.png"),
  title = "IWP versus linear FASH: Monte Carlo power",
  subtitle = plot_subtitle,
  style_profile = "combined"
)
plot_mc_alpha_curves(
  fash_mc_alpha,
  metric = "fdr",
  file = file.path(figure_dir, "iwp_vs_linear_fash_mc_fdr.png"),
  title = "IWP versus linear FASH: Monte Carlo FDR estimate",
  subtitle = plot_subtitle,
  legend_position = "topleft",
  style_profile = "combined"
)
plot_mc_alpha_curves(
  direct_mc_alpha,
  metric = "power",
  file = file.path(figure_dir, "iwp_fash_vs_direct_true_pi0_mc_power.png"),
  title = "IWP FASH versus direct interaction: Monte Carlo power",
  subtitle = paste0(plot_subtitle, "; direct eFDR uses true pi0 = ", true_pi0_expected),
  style_profile = "combined"
)
plot_mc_alpha_curves(
  direct_mc_alpha,
  metric = "fdr",
  file = file.path(figure_dir, "iwp_fash_vs_direct_true_pi0_mc_fdr.png"),
  title = "IWP FASH versus direct interaction: Monte Carlo FDR estimate",
  subtitle = paste0(plot_subtitle, "; direct eFDR uses true pi0 = ", true_pi0_expected),
  legend_position = "bottomright",
  style_profile = "combined"
)

print(fash_mc_alpha_005)
print(direct_mc_alpha_005)
print(mc_pi0)
print(linear_sigma_summary)
