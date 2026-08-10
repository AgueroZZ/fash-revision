#!/usr/bin/env Rscript

# Run the targeted revision simulations for the B-spline and spiky transient
# dynamic-effect mechanisms.

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

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

J <- as.integer(get_arg("--J", "1000"))
seed <- as.integer(get_arg("--seed", "12345"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
overwrite <- as_flag(get_arg("--overwrite", "false"))

output_dir <- file.path(workflowr_root, "output", "revision_simulations")
dirs <- revision_output_dirs(output_dir)

targeted_specs <- list(
  list(
    label = "Cubic B-spline random curve, medium SE",
    dynamic_class = "bspline",
    sd_values = c(0.3, 0.5, 0.8),
    transient_bspline_df = 16
  ),
  list(
    label = "Cubic B-spline random curve, high SE",
    dynamic_class = "bspline",
    sd_values = c(0.5, 0.8, 1.2),
    transient_bspline_df = 16
  ),
  list(
    label = "Spiky local B-spline transient, medium SE",
    dynamic_class = "local_bspline_transient",
    sd_values = c(0.3, 0.5, 0.8),
    transient_bspline_df = 16
  ),
  list(
    label = "Spiky local B-spline transient, high SE",
    dynamic_class = "local_bspline_transient",
    sd_values = c(0.5, 0.8, 1.2),
    transient_bspline_df = 16
  )
)

expected_raw_path <- function(spec) {
  scenario <- paste0("constant_vs_", spec$dynamic_class, "_simplified")
  stem <- paste0(
    scenario,
    "_sd", numeric_vector_label(spec$sd_values),
    "_sigmabeta1",
    "_amp2",
    "_bsdf6",
    "_tdf", spec$transient_bspline_df,
    "_width0p8",
    "_seed", seed,
    "_J", J
  )
  file.path(dirs$raw, paste0(stem, ".rds"))
}

example_file_path <- function(spec) {
  se_label <- if (identical(spec$sd_values, c(0.3, 0.5, 0.8))) "medium" else "high"
  dynamic_label <- if (identical(spec$dynamic_class, "local_bspline_transient")) {
    "local_bspline_transient"
  } else {
    spec$dynamic_class
  }
  file.path(
    dirs$figures,
    paste0(
      "constant_vs_",
      dynamic_label,
      "_examples_",
      se_label,
      "_sd",
      numeric_vector_label(spec$sd_values),
      "_tdf",
      spec$transient_bspline_df,
      "_J",
      J,
      "_seed",
      seed,
      ".png"
    )
  )
}

for (spec in targeted_specs) {
  raw_path <- expected_raw_path(spec)
  if (file.exists(raw_path) && !overwrite) {
    message("Skipping existing output: ", raw_path)
    out <- readRDS(raw_path)
  } else {
    message("Running: ", spec$label, "; J = ", J)
    out <- run_constant_iwp2_simplified_comparison(
      J = J,
      pi_dynamic = 0.20,
      dynamic_class = spec$dynamic_class,
      alpha = 0.05,
      seed = seed,
      sigma_beta = 1,
      estimate_sigma = FALSE,
      num_cores = num_cores,
      sd_values = spec$sd_values,
      dynamic_amplitude = 2,
      bspline_df = 6,
      bspline_coefficient_sd = 1,
      smooth_abrupt_width = 0.8,
      transient_bspline_df = spec$transient_bspline_df,
      output_dir = output_dir,
      save_outputs = TRUE,
      verbose = FALSE
    )
  }

  plot_constant_dynamic_examples(
    out,
    n_per_class = 3,
    file = example_file_path(spec),
    seed = 1
  )
  print(out$summary_table)
}
