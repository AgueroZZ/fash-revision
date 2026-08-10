#!/usr/bin/env Rscript

# Run one raised-cosine truth replacement under the reviewer-facing spiky
# genotype simulation settings.

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

make_replacement_effect_set <- function(old_effect,
                                        truth_mechanism,
                                        time_grid,
                                        evaluation_grid,
                                        seed) {
  spike_counts <- switch(
    truth_mechanism,
    raised_cosine_single = 1L,
    raised_cosine_multipeak = 1:3,
    stop("Unsupported truth mechanism.")
  )
  raised <- simulate_raised_cosine_multipeak_effect_set(
    n_variants = nrow(old_effect$beta_matrix),
    time_grid = time_grid,
    evaluation_grid = evaluation_grid,
    class_probs = c(
      dynamic_bspline = 0.20,
      constant = 0.40,
      zero = 0.40
    ),
    width_levels = c(spiky = 1.50),
    spike_counts = spike_counts,
    relative_amplitude_range = c(0.75, 1.00),
    target_centered_rms = 0.90,
    baseline_sd = 0,
    exact_class_counts = TRUE,
    seed = seed,
    scenario = paste0("internal_spiky_truth_replacement_", truth_mechanism)
  )

  old_dynamic <- which(
    old_effect$unit_info$effect_class == "dynamic_local_bspline_transient"
  )
  raised_dynamic <- which(
    raised$unit_info$effect_class == "dynamic_bspline"
  )
  if (length(old_dynamic) != length(raised_dynamic)) {
    stop("Old and replacement truths have different dynamic counts.")
  }

  beta_matrix <- old_effect$beta_matrix
  beta_matrix[old_dynamic, ] <- raised$beta_matrix[raised_dynamic, ]
  beta_evaluation <- matrix(
    0,
    nrow = nrow(beta_matrix),
    ncol = length(evaluation_grid),
    dimnames = list(
      rownames(beta_matrix),
      sprintf("evaluation_%03d", seq_along(evaluation_grid))
    )
  )
  constant <- which(old_effect$unit_info$effect_class == "constant")
  if (length(constant) > 0) {
    beta_evaluation[constant, ] <- old_effect$beta_matrix[constant, 1]
  }
  beta_evaluation[old_dynamic, ] <-
    raised$beta_evaluation[raised_dynamic, ]

  unit_info <- old_effect$unit_info
  unit_info$effect_class[old_dynamic] <- "dynamic_bspline"
  unit_info$scenario <- paste0(
    "internal_spiky_truth_replacement_",
    truth_mechanism
  )
  metadata_columns <- c(
    "width_label", "width_half", "spike_count", "sign_pattern",
    "cell_id", "latent_id", "baseline", "centered_rms"
  )
  for (column in metadata_columns) {
    template <- raised$unit_info[[column]]
    unit_info[[column]] <- if (is.character(template)) {
      rep(NA_character_, nrow(unit_info))
    } else if (is.integer(template)) {
      rep(NA_integer_, nrow(unit_info))
    } else {
      rep(NA_real_, nrow(unit_info))
    }
    unit_info[[column]][old_dynamic] <- template[raised_dynamic]
  }
  unit_info$peak_centers <- vector("list", nrow(unit_info))
  unit_info$peak_signs <- vector("list", nrow(unit_info))
  unit_info$peak_relative_amplitudes <- vector("list", nrow(unit_info))
  unit_info$peak_centers[old_dynamic] <-
    raised$unit_info$peak_centers[raised_dynamic]
  unit_info$peak_signs[old_dynamic] <-
    raised$unit_info$peak_signs[raised_dynamic]
  unit_info$peak_relative_amplitudes[old_dynamic] <-
    raised$unit_info$peak_relative_amplitudes[raised_dynamic]

  list(
    beta_matrix = beta_matrix,
    beta_evaluation = beta_evaluation,
    evaluation_grid = evaluation_grid,
    unit_info = unit_info,
    settings = list(
      truth_mechanism = truth_mechanism,
      width_half = 1.50,
      spike_counts = spike_counts,
      target_centered_rms = 0.90,
      old_constant_effects_preserved = TRUE,
      dynamic_baseline = 0
    )
  )
}

summarize_curve_geometry <- function(effect_sim, time_grid) {
  dynamic <- effect_sim$unit_info$effect_class == "dynamic_bspline"
  beta <- effect_sim$beta_matrix[dynamic, , drop = FALSE]
  centered <- beta - rowMeans(beta)
  maximum <- apply(abs(centered), 1, max)
  linear_basis <- scaled_time_polynomial(time_grid, degree = 1)
  quadratic_basis <- scaled_time_polynomial(time_grid, degree = 2)
  projection_fraction <- function(x, basis) {
    fitted <- qr.fitted(qr(basis), x)
    sum(fitted^2) / sum(x^2)
  }
  data.frame(
    truth_mechanism = effect_sim$settings$truth_mechanism,
    n_dynamic = nrow(beta),
    centered_rms_mean = mean(sqrt(rowMeans(centered^2))),
    centered_rms_sd = stats::sd(sqrt(rowMeans(centered^2))),
    support_25_mean = mean(rowSums(abs(centered) >= 0.25 * maximum)),
    support_50_mean = mean(rowSums(abs(centered) >= 0.50 * maximum)),
    second_difference_roughness_mean = mean(apply(
      centered,
      1,
      function(x) sum(diff(x, differences = 2)^2) / sum(x^2)
    )),
    linear_projection_mean = mean(apply(
      centered,
      1,
      projection_fraction,
      basis = linear_basis
    )),
    quadratic_projection_mean = mean(apply(
      centered,
      1,
      projection_fraction,
      basis = quadratic_basis
    )),
    stringsAsFactors = FALSE
  )
}

plot_replacement_examples <- function(effect_sim, file) {
  dynamic <- which(effect_sim$unit_info$effect_class == "dynamic_bspline")
  cells <- unique(effect_sim$unit_info$cell_id[dynamic])
  selected <- if (length(cells) == 1L) {
    dynamic[seq_len(min(6L, length(dynamic)))]
  } else {
    vapply(cells, function(cell) {
      dynamic[which(effect_sim$unit_info$cell_id[dynamic] == cell)[1]]
    }, integer(1))
  }
  png(file, width = 1500, height = 950, res = 180)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  })
  par(mfrow = c(2, 3), mar = c(3.2, 3.6, 3.0, 0.8))
  for (index in selected) {
    plot(
      effect_sim$evaluation_grid,
      effect_sim$beta_evaluation[index, ],
      type = "l",
      lwd = 2.3,
      col = "#D95F02",
      xlab = "Time",
      ylab = "True effect",
      main = gsub("__", "; ", effect_sim$unit_info$cell_id[index]),
      cex.main = 0.85
    )
    abline(h = 0, col = "gray70", lty = 3)
    points(
      0:15,
      effect_sim$beta_matrix[index, ],
      pch = 16,
      cex = 0.55,
      col = "#1F78B4"
    )
  }
  if (length(selected) < 6L) {
    for (i in seq_len(6L - length(selected))) plot.new()
  }
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

truth_mechanism <- get_arg("--truth-mechanism", "raised_cosine_single")
allowed_truths <- c("raised_cosine_single", "raised_cosine_multipeak")
if (!truth_mechanism %in% allowed_truths) {
  stop("--truth-mechanism must be one of: ", paste(allowed_truths, collapse = ", "))
}
seed <- as.integer(get_arg("--seed", "12345"))
efdr_permutations <- as.integer(get_arg("--efdr-permutations", "30"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
output_id <- get_arg(
  "--output-id",
  paste0(
    "spiky_truth_replacement_",
    truth_mechanism,
    "_seed",
    seed,
    "_B",
    efdr_permutations
  )
)
if (!is.finite(seed) || efdr_permutations < 1 || num_cores < 1 ||
    !nzchar(output_id)) {
  stop("Invalid truth replacement arguments.")
}

J <- 1000L
n_donors <- 19L
n_covariates <- 5L
time_grid <- make_time_grid()
evaluation_grid <- seq(0, 15, by = 0.1)
class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
true_pi0 <- 0.8
expression_noise_sd <- 1
num_basis <- 20L
component_seeds <- revision_component_seeds(seed)
methods <- c(
  "FASH-IWP1-BF",
  "Direct-linear-LRT-eFDR-true-pi0",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  output_id
)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

genotype_sim <- simulate_genotype_matrix(
  n_donors = n_donors,
  n_variants = J,
  seed = component_seeds[["genotype"]]
)
covariates <- simulate_covariate_matrix(
  n_donors = n_donors,
  n_covariates = n_covariates,
  seed = component_seeds[["covariates"]]
)
old_paired <- simulate_paired_local_bspline_effect_sets(
  candidate_settings = data.frame(
    setting_id = "old_local_bspline_single",
    transient_bspline_df = 16L,
    dynamic_amplitude = 2.75,
    normalization = "center_then_scale",
    stringsAsFactors = FALSE
  ),
  n_variants = J,
  time_grid = time_grid,
  class_probs = c(
    dynamic_local_bspline_transient = 0.20,
    constant = 0.40,
    zero = 0.40
  ),
  transient_bspline_degree = 3L,
  exact_class_counts = TRUE,
  seed = seed,
  scenario_prefix = "internal_spiky_truth_replacement"
)
old_effect <- old_paired$effect_sets$old_local_bspline_single
effect_sim <- make_replacement_effect_set(
  old_effect = old_effect,
  truth_mechanism = truth_mechanism,
  time_grid = time_grid,
  evaluation_grid = evaluation_grid,
  seed = component_seeds[["functional_truth"]]
)

dynamic <- effect_sim$unit_info$effect_class == "dynamic_bspline"
constant <- effect_sim$unit_info$effect_class == "constant"
if (sum(dynamic) != 200L || sum(constant) != 400L ||
    any(abs(rowMeans(effect_sim$beta_matrix[dynamic, , drop = FALSE])) > 1e-10) ||
    any(abs(
      sqrt(rowMeans(effect_sim$beta_matrix[dynamic, , drop = FALSE]^2)) -
        0.90
    ) > 1e-10) ||
    !isTRUE(all.equal(
      effect_sim$beta_matrix[constant, , drop = FALSE],
      old_effect$beta_matrix[constant, , drop = FALSE]
    ))) {
  stop("The paired truth replacement invariants failed.")
}

expression_sim <- simulate_eqtl_expression_from_genotypes(
  G = genotype_sim$G,
  beta_matrix = effect_sim$beta_matrix,
  time_grid = time_grid,
  covariates = covariates,
  expression_noise_sd = expression_noise_sd,
  seed = component_seeds[["expression"]]
)

message("Fitting ", truth_mechanism, " for seed ", seed, ".")
fit <- run_genotype_level_dynamic_eqtl_simulation(
  G = genotype_sim$G,
  time_grid = time_grid,
  covariates = covariates,
  class_probs = class_probs,
  expression_noise_sd = expression_noise_sd,
  alpha = 0.05,
  seed = seed,
  num_cores = num_cores,
  num_basis = num_basis,
  scenario = paste0("internal_spiky_truth_replacement_", truth_mechanism),
  output_dir = output_dir,
  save_outputs = FALSE,
  verbose = FALSE,
  effect_sim = effect_sim,
  expression_sim = expression_sim
)

message("Running direct eFDR with B=", efdr_permutations, ".")
fit <- add_direct_interaction_efdr_results_to_genotype_output(
  out = fit,
  n_permutations = efdr_permutations,
  alpha = 0.05,
  seed = component_seeds[["permutations"]],
  lambda = 0.5,
  pi0_method = "conservative",
  true_pi0 = true_pi0,
  include_true_pi0 = TRUE,
  permute_covariates_with_expression = TRUE,
  num_cores = num_cores,
  overwrite = TRUE,
  verbose = FALSE
)

result_table <- fit$result_table[fit$result_table$method %in% methods, ]
alpha_curve <- fit$alpha_curve[fit$alpha_curve$method %in% methods, ]
alpha_curve$seed <- seed
alpha_curve$truth_mechanism <- truth_mechanism
geometry <- summarize_curve_geometry(effect_sim, time_grid)
geometry$seed <- seed

if (!all(methods %in% unique(result_table$method)) ||
    any(table(result_table$method) != J)) {
  stop("A required method is missing or misaligned.")
}

configuration <- list(
  truth_mechanism = truth_mechanism,
  seed = seed,
  J = J,
  n_donors = n_donors,
  n_covariates = n_covariates,
  expression_noise_sd = expression_noise_sd,
  class_probs = class_probs,
  true_pi0 = true_pi0,
  width_half = 1.50,
  target_centered_rms = 0.90,
  efdr_permutations = efdr_permutations,
  old_truth = list(
    mechanism = "local_cubic_bspline_single",
    df = 16L,
    degree = 3L,
    amplitude = 2.75,
    normalization = "center_then_scale"
  ),
  paired_nonshape_components = TRUE
)
saveRDS(configuration, file.path(output_dir, "configuration.rds"))
write.csv(result_table, file.path(output_dir, "result_table.csv"), row.names = FALSE)
write.csv(alpha_curve, file.path(output_dir, "alpha_curve.csv"), row.names = FALSE)
write.csv(geometry, file.path(output_dir, "geometry.csv"), row.names = FALSE)
plot_replacement_examples(
  effect_sim,
  file.path(figure_dir, "truth_examples.png")
)

saveRDS(
  list(
    configuration = configuration,
    result_table = result_table,
    alpha_curve = alpha_curve,
    geometry = geometry,
    bf_pi0 = constant_component_prior_weight(fit$fash_fits$fash_iwp1_bf),
    effect_sim = effect_sim
  ),
  file.path(output_dir, "replacement_result.rds")
)

message("Saved internal truth replacement result to: ", output_dir)
