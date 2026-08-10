#!/usr/bin/env Rscript

# Aggregate the fixed five-seed confirmation for the internal mixed-shape
# direct-interaction comparison. This script does not build a workflowr page.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) return(".")
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find code/revision_simulations/shared/simulation_functions.R.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0 || hit[1] == length(args)) return(default)
  args[hit[1] + 1]
}

parse_seed_list <- function(x) {
  out <- suppressWarnings(as.integer(trimws(strsplit(x, ",", fixed = TRUE)[[1]])))
  if (length(out) == 0 || anyNA(out) || anyDuplicated(out)) {
    stop("--seed-list must contain unique comma-separated integers.")
  }
  out
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

summarize_shape_power <- function(shape_rows, confidence_level = 0.95) {
  groups <- split(
    shape_rows,
    list(shape_rows$shape_profile, shape_rows$method, shape_rows$alpha),
    drop = TRUE
  )
  out <- lapply(groups, function(x) {
    power_summary <- summarize_mc_values(x$power, confidence_level)
    tp_summary <- summarize_mc_values(x$true_positives, confidence_level)
    data.frame(
      shape_profile = x$shape_profile[1],
      method = x$method[1],
      alpha = x$alpha[1],
      n_dynamic = x$n_dynamic[1],
      n_replications = length(unique(x$seed)),
      mean_true_positives = tp_summary[["mean"]],
      mean_power = power_summary[["mean"]],
      power_ci_lower = pmax(0, power_summary[["lower"]]),
      power_ci_upper = pmin(1, power_summary[["upper"]]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

seeds <- parse_seed_list(get_arg(
  "--seed-list",
  "12345,22345,32345,42345,52345"
))
candidate_id <- get_arg("--candidate-id", "broad075_df10")
input_prefix <- get_arg("--input-prefix", "confirm_broad075_df10_seed")
input_suffix <- get_arg("--input-suffix", "_B100")
output_id <- get_arg("--output-id", "confirm_broad075_df10_B100_5seed")

input_paths <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  paste0(input_prefix, seeds, input_suffix),
  "candidates",
  paste0(candidate_id, ".rds")
)
missing_paths <- input_paths[!file.exists(input_paths)]
if (length(missing_paths) > 0) {
  stop("Missing confirmation caches: ", paste(missing_paths, collapse = ", "))
}
replicates <- lapply(input_paths, readRDS)
if (any(vapply(replicates, function(x) !identical(x$candidate_id, candidate_id), logical(1)))) {
  stop("A confirmation cache has the wrong candidate id.")
}
settings <- lapply(replicates, `[[`, "setting")
if (any(vapply(
  settings[-1],
  function(x) !isTRUE(all.equal(x, settings[[1]])),
  logical(1)
))) {
  stop("Confirmation caches do not share the same candidate setting.")
}

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  output_id
)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

all_alpha <- do.call(rbind, lapply(seq_along(replicates), function(i) {
  out <- replicates[[i]]$alpha_curve
  out$seed <- seeds[i]
  out
}))
mc_alpha <- summarize_mc_alpha_curves(all_alpha)
all_alpha_005 <- all_alpha[abs(all_alpha$alpha - 0.05) < 1e-12, ]
mc_alpha_005 <- mc_alpha[abs(mc_alpha$alpha - 0.05) < 1e-12, ]

all_shape <- do.call(rbind, lapply(seq_along(replicates), function(i) {
  out <- replicates[[i]]$shape_power
  out$seed <- seeds[i]
  out
}))
mc_shape <- summarize_shape_power(all_shape)
mc_shape_005 <- mc_shape[abs(mc_shape$alpha - 0.05) < 1e-12, ]

all_projection <- do.call(rbind, lapply(seq_along(replicates), function(i) {
  out <- replicates[[i]]$projection_summary
  out$seed <- seeds[i]
  out
}))
projection_summary <- aggregate(
  cbind(mean_quadratic_projection, median_quadratic_projection) ~ shape_profile,
  data = all_projection,
  FUN = mean
)

all_pi0 <- do.call(rbind, lapply(seq_along(replicates), function(i) {
  out <- replicates[[i]]$estimated_pi0
  out$seed <- seeds[i]
  out
}))
pi0_summary <- aggregate(estimated_pi0 ~ fit, data = all_pi0, FUN = mean)

alpha_005_wide <- reshape(
  all_alpha_005[, c("seed", "method", "power")],
  idvar = "seed",
  timevar = "method",
  direction = "wide"
)
per_seed_advantage <- data.frame(
  seed = alpha_005_wide$seed,
  raw_iwp_minus_direct_quadratic =
    alpha_005_wide[["power.FASH-IWP1-Raw"]] -
    alpha_005_wide[["power.Direct-quadratic-LRT-eFDR-true-pi0"]],
  bf_iwp_minus_direct_quadratic =
    alpha_005_wide[["power.FASH-IWP1-BF"]] -
    alpha_005_wide[["power.Direct-quadratic-LRT-eFDR-true-pi0"]],
  stringsAsFactors = FALSE
)

write_csv(all_alpha, file.path(output_dir, "all_replicate_alpha_curves.csv"))
write_csv(mc_alpha, file.path(output_dir, "mc_alpha_curve.csv"))
write_csv(all_alpha_005, file.path(output_dir, "all_replicate_alpha005.csv"))
write_csv(mc_alpha_005, file.path(output_dir, "mc_alpha005_summary.csv"))
write_csv(all_shape, file.path(output_dir, "all_replicate_shape_power.csv"))
write_csv(mc_shape, file.path(output_dir, "mc_shape_power_curve.csv"))
write_csv(mc_shape_005, file.path(output_dir, "mc_shape_power_alpha005.csv"))
write_csv(all_projection, file.path(output_dir, "all_replicate_projection_summary.csv"))
write_csv(projection_summary, file.path(output_dir, "projection_summary.csv"))
write_csv(all_pi0, file.path(output_dir, "all_replicate_pi0.csv"))
write_csv(pi0_summary, file.path(output_dir, "pi0_summary.csv"))
write_csv(per_seed_advantage, file.path(output_dir, "per_seed_power_advantage_alpha005.csv"))

subtitle <- paste0(
  length(seeds),
  " seeds; 75% random broad B-spline (df=10), 25% multi-spike; ",
  "direct eFDR uses true pi0 and B=100"
)
plot_mc_alpha_curves(
  mc_curve = mc_alpha,
  metric = "power",
  file = file.path(figure_dir, "power_across_alpha.png"),
  title = "Internal confirmation: IWP FASH versus direct interactions",
  subtitle = subtitle,
  style_profile = "combined"
)
plot_mc_alpha_curves(
  mc_curve = mc_alpha,
  metric = "fdr",
  file = file.path(figure_dir, "fdr_across_alpha.png"),
  title = "Internal confirmation: empirical FDR",
  subtitle = subtitle,
  style_profile = "combined"
)
plot_mc_shape_power_grid(
  mc_curve = mc_shape,
  shape_order = c("broad_random", "spiky"),
  file = file.path(figure_dir, "shape_stratified_power_across_alpha.png"),
  title = "Internal confirmation: power stratified by true-effect shape",
  subtitle = subtitle,
  style_profile = "combined"
)

saveRDS(
  list(
    seeds = seeds,
    candidate_id = candidate_id,
    setting = settings[[1]],
    all_alpha = all_alpha,
    mc_alpha = mc_alpha,
    all_shape = all_shape,
    mc_shape = mc_shape,
    all_projection = all_projection,
    projection_summary = projection_summary,
    all_pi0 = all_pi0,
    pi0_summary = pi0_summary,
    per_seed_advantage = per_seed_advantage
  ),
  file.path(output_dir, "confirmation_summary.rds")
)

message("Saved internal confirmation summary to: ", output_dir)
