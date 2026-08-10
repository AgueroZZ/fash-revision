#!/usr/bin/env Rscript

# Run the genotype-level B-spline dynamic eQTL simulation for the FASH revision.
#
# The default setting mirrors the corrected real-data pipeline: 19 donors,
# 16 time points, 1000 variants, five covariates, and t-statistic SE correction
# with df = n - 7.

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

add_discovery_rate <- function(summary_table) {
  summary_table$discovery_rate <- summary_table$n_discoveries / summary_table$n_units
  summary_table
}

save_method_comparison <- function(out, methods, label, stem, summary_dir) {
  missing_methods <- setdiff(methods, unique(out$result_table$method))
  if (length(missing_methods) > 0) {
    stop(
      "Cannot save ", label, " comparison; missing methods: ",
      paste(missing_methods, collapse = ", ")
    )
  }

  result_table <- out$result_table[out$result_table$method %in% methods, ]
  alpha_curve <- out$alpha_curve[out$alpha_curve$method %in% methods, ]
  alpha_005 <- alpha_curve[abs(alpha_curve$alpha - 0.05) < 1e-12, ]
  alpha_005 <- add_discovery_rate(alpha_005)

  write.csv(
    result_table,
    file = file.path(summary_dir, paste0(stem, "_", label, "_unit_results.csv")),
    row.names = FALSE
  )
  write.csv(
    alpha_curve,
    file = file.path(summary_dir, paste0(stem, "_", label, "_alpha_curve.csv")),
    row.names = FALSE
  )
  write.csv(
    alpha_005,
    file = file.path(summary_dir, paste0(stem, "_", label, "_alpha005_summary.csv")),
    row.names = FALSE
  )

  list(
    result_table = result_table,
    alpha_curve = alpha_curve,
    alpha_005 = alpha_005
  )
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

J <- as.integer(get_arg("--J", "1000"))
n_donors <- as.integer(get_arg("--n-donors", "19"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
seed <- as.integer(get_arg("--seed", "12345"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
num_basis <- as.integer(get_arg("--num-basis", "20"))
expression_noise_sd <- as.numeric(get_arg("--noise-sd", "1"))
dynamic_main_effect_sd <- as.numeric(get_arg(
  "--dynamic-main-effect-sd",
  "1"
))
overwrite <- as_flag(get_arg("--overwrite", "false"))
efdr_permutations <- as.integer(get_arg("--efdr-permutations", "100"))
efdr_seed <- as.integer(get_arg("--efdr-seed", as.character(seed + 10000L)))
efdr_overwrite <- as_flag(get_arg("--efdr-overwrite", "false"))
skip_efdr <- as_flag(get_arg("--skip-efdr", "false"))
include_true_pi0 <- as_flag(get_arg("--include-true-pi0", "true"))
permute_covariates_with_expression <- as_flag(
  get_arg("--permute-covariates-with-expression", "true")
)
if (!is.finite(dynamic_main_effect_sd) || dynamic_main_effect_sd <= 0) {
  stop("--dynamic-main-effect-sd must be positive and finite.")
}

time_grid <- make_time_grid()
class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
scenario <- "genotype_random_bspline_main_effect_dynamic_eqtl"
output_dir <- file.path(workflowr_root, "output", "revision_simulations")
dirs <- revision_output_dirs(output_dir)

stem <- genotype_bspline_eqtl_output_stem(
  n_donors = n_donors,
  n_variants = J,
  time_grid = time_grid,
  n_covariates = n_covariates,
  expression_noise_sd = expression_noise_sd,
  class_probs = class_probs,
  seed = seed,
  scenario = scenario
)
raw_path <- file.path(dirs$raw, paste0(stem, ".rds"))

if (file.exists(raw_path) && !overwrite) {
  message("Skipping existing output: ", raw_path)
  out <- readRDS(raw_path)
} else {
  message("Running genotype-level B-spline eQTL simulation: ", stem)
  out <- run_genotype_level_bspline_eqtl_simulation(
    n_donors = n_donors,
    n_variants = J,
    time_grid = time_grid,
    n_covariates = n_covariates,
    class_probs = class_probs,
    expression_noise_sd = expression_noise_sd,
    dynamic_main_effect_sd = dynamic_main_effect_sd,
    alpha = 0.05,
    seed = seed,
    num_cores = num_cores,
    num_basis = num_basis,
    output_dir = output_dir,
    save_outputs = TRUE,
    verbose = FALSE
  )
}

dynamic <- out$unit_info$effect_class == "dynamic_bspline"
if (!"genetic_main_effect" %in% names(out$unit_info) ||
    any(!is.finite(out$unit_info$genetic_main_effect[dynamic]))) {
  stop("The random B-spline truth is missing its genetic main effects.")
}

true_pi0 <- dynamic_null_proportion(out$unit_info, target = "dynamic")
expected_true_pi0 <- unname(class_probs["constant"] + class_probs["zero"])
if (abs(true_pi0 - expected_true_pi0) > 1e-12) {
  stop("Observed dynamic-null proportion does not match the configured class mixture.")
}

direct_methods <- c("Direct-linear-LRT", "Direct-quadratic-LRT")
if (!all(direct_methods %in% out$result_table$method)) {
  message("Adding direct interaction LRT baselines to cached output.")
  out <- add_direct_interaction_results_to_genotype_output(out, alpha = 0.05)
  saveRDS(out, raw_path)
}

efdr_methods <- c(
  "Direct-linear-LRT-eFDR-true-pi0",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
has_requested_efdr <- !is.null(out$direct_interaction_efdr) &&
  !is.null(out$direct_interaction_efdr$settings$n_permutations) &&
  out$direct_interaction_efdr$settings$n_permutations >= efdr_permutations &&
  all(efdr_methods %in% out$result_table$method) &&
  isTRUE(all.equal(out$direct_interaction_efdr$settings$true_pi0, true_pi0))

if (!skip_efdr && (!has_requested_efdr || efdr_overwrite)) {
  message(
    "Adding direct interaction eFDR baselines with ",
    efdr_permutations,
    " donor-level permutations."
  )
  out <- add_direct_interaction_efdr_results_to_genotype_output(
    out,
    n_permutations = efdr_permutations,
    alpha = 0.05,
    seed = efdr_seed,
    lambda = 0.5,
    pi0_method = "conservative",
    true_pi0 = true_pi0,
    include_true_pi0 = include_true_pi0,
    permute_covariates_with_expression = permute_covariates_with_expression,
    overwrite = efdr_overwrite,
    verbose = TRUE
  )
  saveRDS(out, raw_path)
}

if (!include_true_pi0) {
  stop("The reviewer-facing B-spline comparison requires --include-true-pi0 true.")
}

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
fash_comparison <- save_method_comparison(
  out = out,
  methods = fash_methods,
  label = "iwp_vs_linear_fash",
  stem = stem,
  summary_dir = dirs$summary
)
direct_comparison <- save_method_comparison(
  out = out,
  methods = direct_methods,
  label = "iwp_fash_vs_direct_true_pi0",
  stem = stem,
  summary_dir = dirs$summary
)

write.csv(
  out$result_table,
  file = file.path(dirs$summary, paste0(stem, "_unit_results.csv")),
  row.names = FALSE
)
write.csv(
  out$summary_table,
  file = file.path(dirs$summary, paste0(stem, "_method_summary.csv")),
  row.names = FALSE
)
write.csv(
  out$alpha_curve,
  file = file.path(dirs$summary, paste0(stem, "_alpha_curve.csv")),
  row.names = FALSE
)
write.csv(out$se_correction_summary,
  file = file.path(dirs$summary, paste0(stem, "_se_correction_summary.csv")),
  row.names = FALSE
)
if (!is.null(out$direct_interaction_efdr$threshold_tables)) {
  for (degree_label in names(out$direct_interaction_efdr$threshold_tables)) {
    write.csv(
      out$direct_interaction_efdr$threshold_tables[[degree_label]],
      file = file.path(
        dirs$summary,
        paste0(stem, "_direct_", degree_label, "_efdr_thresholds.csv")
      ),
      row.names = FALSE
    )
  }
}

plot_genotype_eqtl_examples(
  out,
  n_per_class = 3,
  file = file.path(dirs$figures, paste0(stem, "_examples.png")),
  seed = 1,
  true_curve_n = 500
)

plot_subtitle <- paste0(
  "N = ",
  out$settings$n_donors,
  ", J = ",
  out$settings$n_variants,
  ", PCs = ",
  out$settings$n_covariates,
  ", t-corrected SEs"
)

plot_power_alpha_curves(
  fash_comparison$alpha_curve,
  file = file.path(dirs$figures, paste0(stem, "_iwp_vs_linear_fash_power.png")),
  title = "IWP versus linear FASH: power",
  subtitle = plot_subtitle,
  style_profile = "combined"
)

plot_empirical_fdr_alpha_curves(
  fash_comparison$alpha_curve,
  file = file.path(dirs$figures, paste0(stem, "_iwp_vs_linear_fash_fdr.png")),
  title = "IWP versus linear FASH: realized FDP",
  subtitle = plot_subtitle,
  y_label = "Realized FDP",
  style_profile = "combined"
)

plot_power_alpha_curves(
  direct_comparison$alpha_curve,
  file = file.path(dirs$figures, paste0(stem, "_iwp_fash_vs_direct_true_pi0_power.png")),
  title = "IWP FASH versus direct interaction: power",
  subtitle = paste0(plot_subtitle, "; direct eFDR uses true pi0 = ", true_pi0),
  style_profile = "combined"
)

plot_empirical_fdr_alpha_curves(
  direct_comparison$alpha_curve,
  file = file.path(dirs$figures, paste0(stem, "_iwp_fash_vs_direct_true_pi0_fdr.png")),
  title = "IWP FASH versus direct interaction: realized FDP",
  subtitle = paste0(plot_subtitle, "; direct eFDR uses true pi0 = ", true_pi0),
  y_label = "Realized FDP",
  style_profile = "combined"
)

print(fash_comparison$alpha_005)
print(direct_comparison$alpha_005)
print(out$se_correction_summary)
