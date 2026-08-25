#!/usr/bin/env Rscript

# Compare the genuine-null BF distribution with and without the marginal
# t-to-normal standard-error correction in the fixed R1 seed-12345 dataset.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

apply_penalized_alternative_bf_update <- function(raw_fit) {
  grid <- as.numeric(raw_fit$psd_grid)
  full_weights <- numeric(length(grid))
  full_weights[match(raw_fit$prior_weights$psd, grid)] <-
    raw_fit$prior_weights$prior_weight
  null_index <- which(grid == 0)
  if (length(null_index) != 1L || null_index != 1L ||
      full_weights[null_index] <= 0 || full_weights[null_index] >= 1 ||
      sum(full_weights[-null_index]) <= 0) {
    stop("The raw fit does not contain valid penalized mixture weights.")
  }
  pi0_raw <- full_weights[null_index]
  pi_alt <- full_weights[-null_index] / sum(full_weights[-null_index])
  likelihood <- exp(raw_fit$L_matrix)
  bayes_factor <- as.vector(
    likelihood[, -null_index, drop = FALSE] %*% pi_alt /
      likelihood[, null_index]
  )
  pi0_bf <- fashr:::BF_control(
    bayes_factor, plot = FALSE
  )$pi0_hat_star
  lfdr <- stats::plogis(
    stats::qlogis(pi0_bf) - log(bayes_factor)
  )
  list(
    pi0_raw = pi0_raw,
    pi_alt = pi_alt,
    bayes_factor = bayes_factor,
    pi0_bf = pi0_bf,
    lfdr = lfdr
  )
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation",
  "selected_signal_genotype_permutation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "r1_known_truth_genotype_permutation",
  "r1_known_truth_genotype_permutation_helpers.R"
))

if (!requireNamespace("fashr", quietly = TRUE) ||
    !requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("scales", quietly = TRUE)) {
  stop("The fashr, ggplot2, and scales packages are required.")
}

source_path <- file.path(
  workflowr_root, "output", "revision_simulations", "mc",
  "r1_random_bspline_main_effect_profile_sigma_pilot5",
  "full_fits", "seed_12345.rds"
)
output_id <- "r1_t_adjusted_se_null_bf_seed12345_v1"
output_directory <- file.path(
  workflowr_root, "output", "revision_simulations", "internal", output_id
)
figure_directory <- file.path(output_directory, "figures")
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_directory, recursive = TRUE, showWarnings = FALSE)
if (!file.exists(source_path)) {
  stop("The cached R1 seed-12345 source fit is missing.")
}

source_object <- readRDS(source_path)
source_fit <- source_object$fash_fits$fash_iwp1_raw
effect_class <- as.character(source_object$unit_info$effect_class)
alternative_indices <- which(effect_class == "dynamic_bspline")[seq_len(200L)]
null_indices <- which(effect_class == "zero")[seq_len(200L)]
selected_indices <- c(alternative_indices, null_indices)
selected_info <- source_object$unit_info[selected_indices, , drop = FALSE]
unit_keys <- as.character(selected_info$unit_id)
true_null <- c(rep(FALSE, 200L), rep(TRUE, 200L))

beta_hat <- source_object$eqtl_summary$beta_hat[
  selected_indices, , drop = FALSE
]
se_adjusted <- source_object$eqtl_summary$se[
  selected_indices, , drop = FALSE
]
se_unadjusted <- source_object$eqtl_summary$se_uncorrected[
  selected_indices, , drop = FALSE
]
oracle_noise_sd <- source_object$settings$expression_noise_sd
oracle_se_by_unit <- vapply(selected_indices, function(unit_index) {
  design <- cbind(
    intercept = 1,
    G = source_object$genotype[, unit_index],
    source_object$covariates
  )
  sqrt(
    oracle_noise_sd^2 * solve(crossprod(design))["G", "G"]
  )
}, numeric(1))
se_oracle <- matrix(
  oracle_se_by_unit,
  nrow = length(selected_indices),
  ncol = ncol(beta_hat)
)
dimnames(se_oracle) <- dimnames(beta_hat)
true_beta <- source_object$true_beta[selected_indices, , drop = FALSE]
time_grid <- source_object$settings$time_grid
if (is.null(time_grid)) {
  time_grid <- source_fit$fash_data$data_list[[1L]]$x
}

adjusted_data <- source_fit$fash_data$data_list[selected_indices]
adjusted_se_list <- source_fit$fash_data$S[selected_indices]
adjusted_likelihood <- source_fit$L_matrix[selected_indices, , drop = FALSE]
names(adjusted_data) <- unit_keys
names(adjusted_se_list) <- unit_keys
rownames(adjusted_likelihood) <- unit_keys
adjusted_raw_fit <- refit_fash_from_likelihood(
  source_fit = source_fit,
  data_list = adjusted_data,
  se_list = adjusted_se_list,
  likelihood_matrix = adjusted_likelihood,
  unit_keys = unit_keys,
  penalty = 10
)

unadjusted_data <- make_fash_datasets_from_eqtl_summary(
  beta_hat = beta_hat,
  se = se_unadjusted,
  true_beta = true_beta,
  time_grid = time_grid,
  unit_info = selected_info,
  scenario = "r1_t_adjusted_se_null_bf_diagnostic"
)
unadjusted_raw_fit <- fashr::fash(
  Y = "y",
  smooth_var = "x",
  S = "sd",
  data_list = unadjusted_data,
  num_basis = source_fit$settings$num_basis,
  order = source_fit$settings$order,
  betaprec = source_fit$settings$betaprec,
  pred_step = source_fit$settings$pred_step,
  penalty = 10,
  grid = source_fit$psd_grid,
  num_cores = 8,
  verbose = FALSE
)
names(unadjusted_raw_fit$fash_data$data_list) <- unit_keys
names(unadjusted_raw_fit$fash_data$S) <- unit_keys
rownames(unadjusted_raw_fit$L_matrix) <- unit_keys

oracle_data <- make_fash_datasets_from_eqtl_summary(
  beta_hat = beta_hat,
  se = se_oracle,
  true_beta = true_beta,
  time_grid = time_grid,
  unit_info = selected_info,
  scenario = "r1_oracle_se_null_bf_diagnostic"
)
oracle_raw_fit <- fashr::fash(
  Y = "y",
  smooth_var = "x",
  S = "sd",
  data_list = oracle_data,
  num_basis = source_fit$settings$num_basis,
  order = source_fit$settings$order,
  betaprec = source_fit$settings$betaprec,
  pred_step = source_fit$settings$pred_step,
  penalty = 10,
  grid = source_fit$psd_grid,
  num_cores = 8,
  verbose = FALSE
)
names(oracle_raw_fit$fash_data$data_list) <- unit_keys
names(oracle_raw_fit$fash_data$S) <- unit_keys
rownames(oracle_raw_fit$L_matrix) <- unit_keys

fits <- list(
  t_adjusted = adjusted_raw_fit,
  unadjusted = unadjusted_raw_fit,
  oracle_known_sigma = oracle_raw_fit
)
alpha_grid <- seq(0, 0.20, by = 0.001)
unit_rows <- list()
summary_rows <- list()
curve_rows <- list()

for (se_method in names(fits)) {
  raw_fit <- fits[[se_method]]
  update <- apply_penalized_alternative_bf_update(raw_fit)
  null_summary <- summarize_null_bf(
    bf = update$bayes_factor,
    true_null = true_null,
    arm = se_method
  )
  curve <- known_truth_alpha_curve(
    lfdr = update$lfdr,
    true_null = true_null,
    alpha_grid = alpha_grid,
    arm = "genuine_null_baseline",
    fit_stage = se_method
  )
  names(curve)[names(curve) == "fit_stage"] <- "se_method"
  alpha005 <- curve[abs(curve$alpha - 0.05) < 1e-12, , drop = FALSE]
  summary_rows[[se_method]] <- data.frame(
    se_method = se_method,
    raw_penalized_pi0 = update$pi0_raw,
    bf_adjusted_pi0 = update$pi0_bf,
    null_summary[, setdiff(names(null_summary), "arm")],
    alpha005_n_discoveries = alpha005$n_discoveries,
    alpha005_false_discoveries = alpha005$false_discoveries,
    alpha005_realized_fdp = alpha005$realized_fdp,
    alpha005_power = alpha005$power,
    stringsAsFactors = FALSE
  )
  unit_rows[[se_method]] <- data.frame(
    unit_key = unit_keys,
    true_null = true_null,
    se_method = se_method,
    bayes_factor = update$bayes_factor,
    lfdr = update$lfdr,
    stringsAsFactors = FALSE
  )
  curve_rows[[se_method]] <- curve
}

summary_table <- do.call(rbind, summary_rows)
unit_results <- do.call(rbind, unit_rows)
alpha_curves <- do.call(rbind, curve_rows)

summarize_ratio <- function(values, group, ratio_name) {
  quantiles <- stats::quantile(
    values, probs = c(0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
    type = 8, names = FALSE
  )
  data.frame(
    group = group,
    ratio = ratio_name,
    n = length(values),
    mean = mean(values),
    q25 = quantiles[1L],
    median = quantiles[2L],
    q75 = quantiles[3L],
    q90 = quantiles[4L],
    q95 = quantiles[5L],
    q99 = quantiles[6L],
    maximum = max(values),
    stringsAsFactors = FALSE
  )
}
se_matrices <- list(
  t_adjusted = se_adjusted,
  unadjusted = se_unadjusted,
  oracle_known_sigma = se_oracle
)
ratio_pairs <- list(
  t_adjusted_over_unadjusted = c("t_adjusted", "unadjusted"),
  unadjusted_over_oracle = c("unadjusted", "oracle_known_sigma"),
  t_adjusted_over_oracle = c("t_adjusted", "oracle_known_sigma")
)
se_ratio_summary <- do.call(rbind, lapply(names(ratio_pairs), function(name) {
  pair <- ratio_pairs[[name]]
  ratio_matrix <- se_matrices[[pair[1L]]] / se_matrices[[pair[2L]]]
  rbind(
    summarize_ratio(as.vector(ratio_matrix), "all_selected_units", name),
    summarize_ratio(
      as.vector(ratio_matrix[true_null, , drop = FALSE]),
      "genuine_null_units",
      name
    )
  )
}))

z_summary <- do.call(rbind, lapply(names(se_matrices), function(se_method) {
  z_values <- as.vector(
    beta_hat[true_null, , drop = FALSE] /
      se_matrices[[se_method]][true_null, , drop = FALSE]
  )
  quantiles <- stats::quantile(
    z_values, probs = c(0.01, 0.05, 0.50, 0.95, 0.99),
    type = 8, names = FALSE
  )
  data.frame(
    se_method = se_method,
    n = length(z_values),
    mean = mean(z_values),
    sd = stats::sd(z_values),
    q01 = quantiles[1L],
    q05 = quantiles[2L],
    median = quantiles[3L],
    q95 = quantiles[4L],
    q99 = quantiles[5L],
    stringsAsFactors = FALSE
  )
}))

adjusted_units <- unit_results[
  unit_results$se_method == "t_adjusted", , drop = FALSE
]
unadjusted_units <- unit_results[
  unit_results$se_method == "unadjusted", , drop = FALSE
]
unadjusted_units <- unadjusted_units[
  match(adjusted_units$unit_key, unadjusted_units$unit_key), , drop = FALSE
]
oracle_units <- unit_results[
  unit_results$se_method == "oracle_known_sigma", , drop = FALSE
]
oracle_units <- oracle_units[
  match(adjusted_units$unit_key, oracle_units$unit_key), , drop = FALSE
]
paired_units <- data.frame(
  unit_key = adjusted_units$unit_key,
  true_null = adjusted_units$true_null,
  t_adjusted_bf = adjusted_units$bayes_factor,
  unadjusted_bf = unadjusted_units$bayes_factor,
  oracle_bf = oracle_units$bayes_factor,
  log10_bf_ratio_unadjusted_over_adjusted = log10(
    unadjusted_units$bayes_factor / adjusted_units$bayes_factor
  ),
  t_adjusted_lfdr = adjusted_units$lfdr,
  unadjusted_lfdr = unadjusted_units$lfdr,
  oracle_lfdr = oracle_units$lfdr,
  stringsAsFactors = FALSE
)

validation <- data.frame(
  check = c(
    "Cached adjusted SE matrix matches selected FASH data",
    "Unadjusted FASH data use the cached uncorrected SE matrix",
    "Oracle FASH data use the exact known-sigma SE matrix",
    "Adjusted raw pi0 reproduces the prior five-seed diagnostic",
    "Adjusted BF pi0 reproduces the prior corrected diagnostic",
    "Adjusted null mean BF reproduces the prior corrected diagnostic"
  ),
  maximum_absolute_error = c(
    max(abs(
      do.call(rbind, adjusted_se_list) - se_adjusted
    )),
    max(abs(
      do.call(rbind, lapply(unadjusted_data, function(x) x$sd)) -
        se_unadjusted
    )),
    max(abs(
      do.call(rbind, lapply(oracle_data, function(x) x$sd)) - se_oracle
    )),
    abs(summary_table$raw_penalized_pi0[
      summary_table$se_method == "t_adjusted"
    ] - 0.4876576),
    abs(summary_table$bf_adjusted_pi0[
      summary_table$se_method == "t_adjusted"
    ] - 0.555),
    abs(summary_table$mean_bf[
      summary_table$se_method == "t_adjusted"
    ] - 0.2295891)
  ),
  tolerance = c(1e-12, 1e-12, 1e-12, 1e-7, 1e-12, 1e-7),
  stringsAsFactors = FALSE
)
validation$pass <- validation$maximum_absolute_error <= validation$tolerance
if (!all(validation$pass)) {
  stop("At least one t-adjusted SE diagnostic validation failed.")
}

write_csv(summary_table, file.path(output_directory, "method_summary.csv"))
write_csv(se_ratio_summary, file.path(output_directory, "se_ratio_summary.csv"))
write_csv(z_summary, file.path(output_directory, "null_z_summary.csv"))
write_csv(unit_results, file.path(output_directory, "unit_results.csv"))
write_csv(paired_units, file.path(output_directory, "paired_units.csv"))
write_csv(alpha_curves, file.path(output_directory, "alpha_curves.csv"))
write_csv(validation, file.path(output_directory, "validation.csv"))

plot_units <- unit_results[unit_results$true_null, , drop = FALSE]
plot_units$se_label <- factor(
  plot_units$se_method,
  levels = c("t_adjusted", "unadjusted", "oracle_known_sigma"),
  labels = c(
    "t-to-normal adjusted SE",
    "Unadjusted OLS SE",
    "Oracle SE with known sigma"
  )
)
bf_figure <- ggplot2::ggplot(
  plot_units,
  ggplot2::aes(x = log10(bayes_factor), color = se_label, fill = se_label)
) +
  ggplot2::geom_density(alpha = 0.16, linewidth = 1.0) +
  ggplot2::geom_vline(
    xintercept = 0, color = "#777777", linetype = "dotted"
  ) +
  ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
  ggplot2::scale_fill_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
  ggplot2::labs(
    title = "Removing the t-to-normal SE correction shifts null BF upward",
    subtitle = "R1 seed 12345; 200 genuine null units with identical beta estimates",
    x = expression(log[10](BF)),
    y = "Density",
    color = NULL,
    fill = NULL
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "top"
  )
ggplot2::ggsave(
  file.path(figure_directory, "null_bf_density.png"),
  bf_figure, width = 8.4, height = 5.2, dpi = 300, bg = "white"
)

curve_plot <- alpha_curves[alpha_curves$alpha <= 0.10, , drop = FALSE]
curve_plot$se_label <- factor(
  curve_plot$se_method,
  levels = c("t_adjusted", "unadjusted", "oracle_known_sigma"),
  labels = c(
    "t-to-normal adjusted SE",
    "Unadjusted OLS SE",
    "Oracle SE with known sigma"
  )
)
fdr_figure <- ggplot2::ggplot(
  curve_plot,
  ggplot2::aes(
    x = alpha, y = realized_fdp, color = se_label,
    linetype = se_label
  )
) +
  ggplot2::geom_abline(
    intercept = 0, slope = 1, color = "#777777", linetype = "dotted"
  ) +
  ggplot2::geom_line(linewidth = 1.0) +
  ggplot2::scale_color_manual(values = c("#0072B2", "#D55E00", "#009E73")) +
  ggplot2::scale_linetype_manual(values = c("solid", "longdash", "dotdash")) +
  ggplot2::scale_x_continuous(
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::label_percent(accuracy = 1)
  ) +
  ggplot2::labs(
    title = "Single-seed empirical FDP is sensitive to the SE scale",
    subtitle = "R1 seed 12345; corrected penalized-alternative BF update",
    x = "Nominal cumulative-FDR level",
    y = "Realized FDP",
    color = NULL,
    linetype = NULL
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "top"
  )
ggplot2::ggsave(
  file.path(figure_directory, "fdr_curve.png"),
  fdr_figure, width = 8.4, height = 5.2, dpi = 300, bg = "white"
)

saveRDS(
  list(
    summary_table = summary_table,
    se_ratio_summary = se_ratio_summary,
    z_summary = z_summary,
    unit_results = unit_results,
    paired_units = paired_units,
    alpha_curves = alpha_curves,
    validation = validation
  ),
  file.path(output_directory, "results.rds")
)

message("Wrote t-adjusted SE null-BF diagnostic to: ", output_directory)
