#!/usr/bin/env Rscript

# Run one seed of the internal raised-cosine multi-peak factorial experiment.

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

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

cumulative_lfdr_score <- function(lfdr) {
  ordering <- order(lfdr)
  cumulative <- cumsum(lfdr[ordering]) / seq_along(lfdr)
  cumulative[match(seq_along(lfdr), ordering)]
}

curve_geometry <- function(beta_matrix, time_grid) {
  centered <- t(apply(beta_matrix, 1, function(x) x - mean(x)))
  maximum <- apply(abs(centered), 1, max)
  nonconstant <- maximum > sqrt(.Machine$double.eps)
  linear_basis <- scaled_time_polynomial(time_grid, degree = 1)
  quadratic_basis <- scaled_time_polynomial(time_grid, degree = 2)
  linear_qr <- qr(linear_basis)
  quadratic_qr <- qr(quadratic_basis)
  projection_fraction <- function(x, basis_qr) {
    if (sum(x^2) <= .Machine$double.eps) return(0)
    fitted <- qr.fitted(basis_qr, x)
    sum(fitted^2) / sum(x^2)
  }
  support_25 <- integer(nrow(beta_matrix))
  support_50 <- integer(nrow(beta_matrix))
  support_25[nonconstant] <- rowSums(
    abs(centered[nonconstant, , drop = FALSE]) >=
      0.25 * maximum[nonconstant]
  )
  support_50[nonconstant] <- rowSums(
    abs(centered[nonconstant, , drop = FALSE]) >=
      0.50 * maximum[nonconstant]
  )
  roughness <- numeric(nrow(beta_matrix))
  roughness[nonconstant] <- apply(
    beta_matrix[nonconstant, , drop = FALSE],
    1,
    function(x) sum(diff(x, differences = 2)^2) / sum((x - mean(x))^2)
  )
  data.frame(
    support_25 = support_25,
    support_50 = support_50,
    second_difference_roughness = roughness,
    linear_projection = apply(
      centered,
      1,
      projection_fraction,
      basis_qr = linear_qr
    ),
    quadratic_projection = apply(
      centered,
      1,
      projection_fraction,
      basis_qr = quadratic_qr
    ),
    stringsAsFactors = FALSE
  )
}

plot_truth_grid <- function(effect_sim, sign_mode, file) {
  width_order <- names(effect_sim$settings$width_levels)
  spike_order <- 1:3
  png(file, width = 1500, height = 1500, res = 180)
  par(
    mfrow = c(length(width_order), length(spike_order)),
    mar = c(3.0, 3.5, 2.5, 0.8),
    oma = c(0, 0, 2.8, 0)
  )
  for (width_label in width_order) {
    for (spike_count in spike_order) {
      pattern <- if (spike_count == 1L) "single" else sign_mode
      candidates <- which(
        effect_sim$unit_info$effect_class == "dynamic_bspline" &
          effect_sim$unit_info$width_label == width_label &
          effect_sim$unit_info$spike_count == spike_count &
          effect_sim$unit_info$sign_pattern == pattern
      )
      if (length(candidates) == 0) {
        plot.new()
        next
      }
      index <- candidates[1]
      centered <- effect_sim$beta_evaluation[index, ] -
        mean(effect_sim$beta_matrix[index, ])
      plot(
        effect_sim$evaluation_grid,
        centered,
        type = "l",
        lwd = 2.2,
        col = "#D95F02",
        xlab = "Time",
        ylab = "Centered true effect",
        main = paste0(
          width_label,
          "; ",
          spike_count,
          if (spike_count == 1L) " peak" else " peaks"
        ),
        cex.main = 0.78
      )
      abline(h = 0, col = "grey70", lty = 3)
      points(
        0:15,
        effect_sim$beta_matrix[index, ] -
          mean(effect_sim$beta_matrix[index, ]),
        pch = 16,
        cex = 0.45,
        col = "#1F78B4"
      )
    }
  }
  mtext(
    paste0(
      "Raised-cosine multi-peak examples: ",
      gsub("-", " ", sign_mode)
    ),
    outer = TRUE,
    cex = 1.1
  )
  dev.off()
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

J <- as.integer(get_arg("--J", "2400"))
n_donors <- as.integer(get_arg("--n-donors", "19"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
expression_noise_sd <- as.numeric(get_arg("--noise-sd", "1"))
num_basis <- as.integer(get_arg("--num-basis", "20"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
efdr_permutations <- as.integer(get_arg("--efdr-permutations", "30"))
seed <- as.integer(get_arg("--seed", "12345"))
output_id <- get_arg(
  "--output-id",
  paste0("multipeak_factorial_seed", seed, "_B", efdr_permutations)
)

if (J != 2400 || n_donors != 19 || n_covariates != 5 ||
    !is.finite(expression_noise_sd) || expression_noise_sd <= 0 ||
    num_basis < 2 || num_cores < 1 ||
    efdr_permutations < 1 || !nzchar(output_id)) {
  stop("Invalid internal multi-peak factorial arguments.")
}

time_grid <- make_time_grid()
evaluation_grid <- seq(min(time_grid), max(time_grid), by = 0.1)
alpha_grid <- seq(0.005, 0.20, by = 0.005)
class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
true_pi0 <- 0.8
width_levels <- c(
  not_spiky = 3.00,
  mildly_spiky = 2.25,
  spiky = 1.50,
  very_spiky = 0.90
)
methods <- c(
  "FASH-IWP1-BF",
  "Direct-linear-LRT-eFDR-true-pi0",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
component_seeds <- revision_component_seeds(seed)

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  output_id
)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

configuration <- list(
  output_id = output_id,
  internal_only = TRUE,
  seed = seed,
  J = J,
  n_donors = n_donors,
  n_covariates = n_covariates,
  time_grid = time_grid,
  evaluation_grid = evaluation_grid,
  expression_noise_sd = expression_noise_sd,
  num_basis = num_basis,
  efdr_permutations = efdr_permutations,
  class_probs = class_probs,
  true_pi0 = true_pi0,
  width_levels = width_levels,
  spike_counts = 1:3,
  sign_patterns = c("single", "same-sign", "alternating-sign"),
  relative_amplitude_range = c(0.75, 1.00),
  target_centered_rms = 0.90,
  baseline_sd = 1,
  direct_efdr_uses_true_pi0 = TRUE
)
saveRDS(configuration, file.path(output_dir, "configuration.rds"))

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
effect_sim <- simulate_raised_cosine_multipeak_effect_set(
  n_variants = J,
  time_grid = time_grid,
  evaluation_grid = evaluation_grid,
  class_probs = class_probs,
  width_levels = width_levels,
  spike_counts = 1:3,
  relative_amplitude_range = c(0.75, 1.00),
  target_centered_rms = 0.90,
  baseline_sd = 1,
  exact_class_counts = TRUE,
  seed = component_seeds[["functional_truth"]]
)
dynamic <- effect_sim$unit_info$effect_class == "dynamic_bspline"
if (sum(dynamic) != 480L ||
    length(unique(effect_sim$unit_info$cell_id[dynamic])) != 20L ||
    any(table(effect_sim$unit_info$cell_id[dynamic]) != 24L) ||
    any(abs(effect_sim$unit_info$centered_rms[dynamic] - 0.90) > 1e-10)) {
  stop("The raised-cosine factorial allocation is invalid.")
}

expression_sim <- simulate_eqtl_expression_from_genotypes(
  G = genotype_sim$G,
  beta_matrix = effect_sim$beta_matrix,
  time_grid = time_grid,
  covariates = covariates,
  expression_noise_sd = expression_noise_sd,
  seed = component_seeds[["expression"]]
)

message("Fitting FASH and per-time eQTL summaries.")
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
  scenario = "internal_raised_cosine_multipeak_factorial",
  output_dir = output_dir,
  save_outputs = FALSE,
  verbose = FALSE,
  effect_sim = effect_sim,
  expression_sim = expression_sim
)

message("Running direct interaction eFDR with B=", efdr_permutations, ".")
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
if (!all(methods %in% unique(result_table$method))) {
  stop("A required factorial-comparison method is missing.")
}

geometry <- curve_geometry(effect_sim$beta_matrix, time_grid)
unit_diagnostics <- cbind(effect_sim$unit_info, geometry)
unit_diagnostics$bf_lfdr <- fit$fash_fits$fash_iwp1_bf$lfdr
unit_diagnostics$bf_cumulative_fdr <- cumulative_lfdr_score(
  unit_diagnostics$bf_lfdr
)
unit_diagnostics$log_bayes_factor <- log(fit$fash_fits$fash_iwp1_bf$BF)
for (method in methods) {
  method_rows <- result_table[result_table$method == method, ]
  score_name <- switch(
    method,
    "FASH-IWP1-BF" = "fash_bf_score",
    "Direct-linear-LRT-eFDR-true-pi0" = "direct_linear_score",
    "Direct-quadratic-LRT-eFDR-true-pi0" = "direct_quadratic_score"
  )
  unit_diagnostics[[score_name]] <- method_rows$adjusted_score[
    match(unit_diagnostics$unit_index, method_rows$unit_index)
  ]
}

cell_alpha <- do.call(rbind, lapply(methods, function(method) {
  method_rows <- result_table[result_table$method == method, ]
  do.call(rbind, lapply(alpha_grid, function(alpha) {
    selected <- method_rows$unit_index[
      is.finite(method_rows$adjusted_score) &
        method_rows$adjusted_score <= alpha
    ]
    do.call(rbind, lapply(
      split(
        effect_sim$unit_info[dynamic, ],
        effect_sim$unit_info$cell_id[dynamic]
      ),
      function(cell) {
        data.frame(
          seed = seed,
          method = method,
          alpha = alpha,
          cell_id = cell$cell_id[1],
          spike_count = cell$spike_count[1],
          width_label = cell$width_label[1],
          width_half = cell$width_half[1],
          sign_pattern = cell$sign_pattern[1],
          n_dynamic = nrow(cell),
          true_positives = sum(cell$unit_index %in% selected),
          power = sum(cell$unit_index %in% selected) / nrow(cell),
          stringsAsFactors = FALSE
        )
      }
    ))
  }))
}))
rownames(cell_alpha) <- NULL

global_alpha <- fit$alpha_curve[fit$alpha_curve$method %in% methods, ]
global_alpha$seed <- seed
alpha_005 <- cell_alpha[abs(cell_alpha$alpha - 0.05) < 1e-12, ]
geometry_summary <- aggregate(
  cbind(
    support_25,
    support_50,
    second_difference_roughness,
    linear_projection,
    quadratic_projection,
    log_bayes_factor
  ) ~ spike_count + width_label + width_half + sign_pattern,
  data = unit_diagnostics[dynamic, ],
  FUN = mean
)

write_csv(result_table, file.path(output_dir, "result_table.csv"))
unit_diagnostics_csv <- unit_diagnostics
for (column in c(
  "peak_centers",
  "peak_signs",
  "peak_relative_amplitudes"
)) {
  unit_diagnostics_csv[[column]] <- vapply(
    unit_diagnostics_csv[[column]],
    paste,
    collapse = ";",
    FUN.VALUE = character(1)
  )
}
write_csv(
  unit_diagnostics_csv,
  file.path(output_dir, "unit_diagnostics.csv")
)
write_csv(cell_alpha, file.path(output_dir, "cell_power_alpha.csv"))
write_csv(alpha_005, file.path(output_dir, "cell_power_alpha005.csv"))
write_csv(global_alpha, file.path(output_dir, "global_alpha.csv"))
write_csv(geometry_summary, file.path(output_dir, "cell_geometry_summary.csv"))

plot_truth_grid(
  effect_sim,
  sign_mode = "same-sign",
  file = file.path(figure_dir, "truth_examples_same_sign.png")
)
plot_truth_grid(
  effect_sim,
  sign_mode = "alternating-sign",
  file = file.path(figure_dir, "truth_examples_alternating_sign.png")
)

compact <- list(
  configuration = configuration,
  effect_sim = effect_sim,
  result_table = result_table,
  unit_diagnostics = unit_diagnostics,
  cell_alpha = cell_alpha,
  global_alpha = global_alpha,
  geometry_summary = geometry_summary,
  bf_pi0 = constant_component_prior_weight(fit$fash_fits$fash_iwp1_bf),
  bf_prior_weights = fit$fash_fits$fash_iwp1_bf$prior_weights
)
saveRDS(compact, file.path(output_dir, "factorial_result.rds"))

message("Saved internal multi-peak factorial result to: ", output_dir)
