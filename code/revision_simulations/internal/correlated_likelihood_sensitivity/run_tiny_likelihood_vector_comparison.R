#!/usr/bin/env Rscript

# Compare FASH likelihood-vector shapes under diagonal, C1, C2, and C3 errors.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath("."))
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return(normalizePath("coderepo-local"))
  }
  stop("Could not find the workflowr repository root.")
}

get_argument <- function(name, default = NULL) {
  arguments <- commandArgs(trailingOnly = TRUE)
  prefix <- paste0(name, "=")
  matches <- startsWith(arguments, prefix)
  if (!any(matches)) {
    return(default)
  }
  if (sum(matches) > 1L) {
    stop("Argument was supplied more than once: ", name)
  }
  sub(prefix, "", arguments[matches])
}

softmax <- function(x) {
  shifted <- x - max(x)
  weights <- exp(shifted)
  weights / sum(weights)
}

safe_vector_correlation <- function(x, y) {
  if (stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }
  unname(stats::cor(x, y))
}

selection_mode <- get_argument("--selection", "random")
allowed_selection_modes <- c("random", "highest_diagonal_raw_lfdr")
if (!selection_mode %in% allowed_selection_modes) {
  stop(
    "Unknown selection mode: ", selection_mode,
    ". Expected one of: ", paste(allowed_selection_modes, collapse = ", ")
  )
}
output_id <- get_argument(
  "--output-id",
  "tiny_likelihood_vector_comparison"
)
if (!grepl("^[A-Za-z0-9_.-]+$", output_id)) {
  stop("The output ID may only contain letters, digits, dots, dashes, and underscores.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "r4_correlated_errors",
  "real_data_correlation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "correlated_likelihood_sensitivity",
  "correlated_likelihood_helpers.R"
))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required.")
}

analysis_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "correlated_likelihood_sensitivity"
)
covariance_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "mashr_mean_z_null_correlation"
)
output_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  output_id
)
if (file.exists(output_directory)) {
  stop("Refusing to overwrite the existing tiny-experiment output directory.")
}

fit_bundle <- readRDS(file.path(analysis_directory, "fit_bundle.rds"))
covariance_analysis <- readRDS(file.path(
  covariance_directory,
  "mashr_mean_z_null_correlation.rds"
))

pair_metadata <- fit_bundle$pair_metadata
candidate_metadata <- covariance_analysis$candidate_metadata
candidate_beta_hat <- as.matrix(covariance_analysis$candidate_beta_hat)
candidate_adjusted_se <- as.matrix(
  covariance_analysis$candidate_adjusted_se
)
if (!identical(pair_metadata$pair_key, candidate_metadata$pair_key) ||
    !identical(rownames(candidate_beta_hat), pair_metadata$pair_key) ||
    !identical(rownames(candidate_adjusted_se), pair_metadata$pair_key) ||
    !identical(names(fit_bundle$raw_fits), c("diagonal", "C1", "C2"))) {
  stop("The fixed thinning, covariance cache, and fitted likelihoods disagree.")
}

screen1_indices <- as.integer(
  covariance_analysis$estimates$maximum_z$selected_candidate_rows
)
centered_screen1 <- weighted_center_standardize(
  candidate_beta_hat[screen1_indices, , drop = FALSE],
  candidate_adjusted_se[screen1_indices, , drop = FALSE]
)$standardized_residual
C3 <- stats::cor(centered_screen1)
C3 <- (C3 + t(C3)) / 2
C3_diagnostics <- validate_shared_correlation(C3, name = "C3")

sampling_seed <- 20260817L
n_draws <- 3L
units_per_draw <- 8L
if (selection_mode == "random") {
  set.seed(sampling_seed)
  selected_rows <- sample.int(
    nrow(pair_metadata),
    size = n_draws * units_per_draw,
    replace = FALSE
  )
  selection_description <- paste(
    "A reproducible simple random sample of 24 units from the fixed",
    "one-variant-per-gene thinning"
  )
  display_group_prefix <- "D"
  display_group_legend <- "Random draw"
} else {
  diagonal_raw_lfdr <- as.numeric(fit_bundle$raw_fits$diagonal$lfdr)
  if (length(diagonal_raw_lfdr) != nrow(pair_metadata) ||
      any(!is.finite(diagonal_raw_lfdr))) {
    stop("The diagonal raw lfdr vector is incomplete or misaligned.")
  }
  selected_rows <- order(
    diagonal_raw_lfdr,
    decreasing = TRUE,
    method = "radix"
  )[seq_len(n_draws * units_per_draw)]
  selection_description <- paste(
    "The 24 units with the highest raw lfdr under the diagonal fit;",
    "ties are resolved by fixed thinning-row order"
  )
  display_group_prefix <- "B"
  display_group_legend <- "Display block"
}
selected_metadata <- pair_metadata[selected_rows, c(
  "fash_index", "pair_key", "gene_id", "variant_id"
)]
selected_metadata$diagonal_raw_lfdr <- as.numeric(
  fit_bundle$raw_fits$diagonal$lfdr[selected_rows]
)
selected_metadata$C1_raw_lfdr <- as.numeric(
  fit_bundle$raw_fits$C1$lfdr[selected_rows]
)
selected_metadata$C2_raw_lfdr <- as.numeric(
  fit_bundle$raw_fits$C2$lfdr[selected_rows]
)
selected_metadata$draw <- rep(seq_len(n_draws), each = units_per_draw)
selected_metadata$unit_in_draw <- rep(seq_len(units_per_draw), n_draws)
selected_metadata$display_label <- sprintf(
  paste0(display_group_prefix, "%d-U%d: %s / %s"),
  selected_metadata$draw,
  selected_metadata$unit_in_draw,
  sub("^ENSG0+", "ENSG", selected_metadata$gene_id),
  selected_metadata$variant_id
)
rownames(selected_metadata) <- NULL

selected_beta_hat <- candidate_beta_hat[selected_rows, , drop = FALSE]
selected_adjusted_se <- candidate_adjusted_se[selected_rows, , drop = FALSE]
settings <- fit_bundle$configuration$raw_settings
psd_grid <- as.numeric(fit_bundle$configuration$psd_grid)
time_grid <- as.numeric(fit_bundle$configuration$time_grid)

message("Computing the C3 likelihood for 24 selected units.")
fit_start <- proc.time()[["elapsed"]]
C3_result <- fit_fash_with_shared_correlation(
  beta_hat = selected_beta_hat,
  adjusted_se = selected_adjusted_se,
  time_grid = time_grid,
  correlation = C3,
  settings = settings,
  psd_grid = psd_grid,
  num_cores = 4L,
  verbose = FALSE
)
C3_elapsed_seconds <- proc.time()[["elapsed"]] - fit_start

likelihood_matrices <- lapply(
  fit_bundle$raw_fits,
  function(fit) as.matrix(fit$L_matrix[selected_rows, , drop = FALSE])
)
likelihood_matrices$C3 <- as.matrix(C3_result$fit$L_matrix)
method_ids <- c("diagonal", "C1", "C2", "C3")
method_labels <- c(
  diagonal = "Diagonal",
  C1 = "C1: Screen 1",
  C2 = "C2: Screen 2",
  C3 = "C3: Screen 1 + within-unit centering"
)
if (!identical(names(likelihood_matrices), method_ids) ||
    !all(vapply(likelihood_matrices, function(matrix) {
      identical(dim(matrix), c(nrow(selected_metadata), length(psd_grid))) &&
        all(is.finite(matrix))
    }, logical(1)))) {
  stop("The selected likelihood matrices are incomplete or misaligned.")
}

null_column <- which(psd_grid == 0)
if (length(null_column) != 1L) {
  stop("The PSD grid must contain exactly one null component.")
}
centered_likelihoods <- lapply(likelihood_matrices, function(matrix) {
  sweep(matrix, 1L, matrix[, null_column], `-`)
})

likelihood_long <- do.call(rbind, lapply(method_ids, function(method_id) {
  matrix <- centered_likelihoods[[method_id]]
  data.frame(
    pair_key = rep(selected_metadata$pair_key, each = length(psd_grid)),
    draw = rep(selected_metadata$draw, each = length(psd_grid)),
    unit_in_draw = rep(
      selected_metadata$unit_in_draw,
      each = length(psd_grid)
    ),
    display_label = rep(
      selected_metadata$display_label,
      each = length(psd_grid)
    ),
    method_id = method_id,
    method_label = unname(method_labels[method_id]),
    psd = rep(psd_grid, times = nrow(selected_metadata)),
    prior_sd = sqrt(rep(psd_grid, times = nrow(selected_metadata))),
    centered_log_likelihood = as.vector(t(matrix)),
    stringsAsFactors = FALSE
  )
}))

method_unit_summary <- do.call(rbind, lapply(method_ids, function(method_id) {
  matrix <- centered_likelihoods[[method_id]]
  alternative_columns <- setdiff(seq_along(psd_grid), null_column)
  best_alternative_column <- apply(
    matrix[, alternative_columns, drop = FALSE],
    1L,
    which.max
  )
  best_alternative_column <- alternative_columns[best_alternative_column]
  data.frame(
    pair_key = selected_metadata$pair_key,
    draw = selected_metadata$draw,
    unit_in_draw = selected_metadata$unit_in_draw,
    method_id = method_id,
    method_label = unname(method_labels[method_id]),
    best_alternative_log_likelihood_gain = matrix[cbind(
      seq_len(nrow(matrix)),
      best_alternative_column
    )],
    best_alternative_psd = psd_grid[best_alternative_column],
    alternative_beats_null = matrix[cbind(
      seq_len(nrow(matrix)),
      best_alternative_column
    )] > 0,
    stringsAsFactors = FALSE
  )
}))

method_profile_summary <- do.call(rbind, lapply(method_ids, function(
    method_id) {
  table <- method_unit_summary[method_unit_summary$method_id == method_id, ]
  data.frame(
    method_id = method_id,
    method_label = unname(method_labels[method_id]),
    median_best_alternative_log_likelihood_gain = stats::median(
      table$best_alternative_log_likelihood_gain
    ),
    lower_quartile_best_alternative_log_likelihood_gain = unname(
      stats::quantile(table$best_alternative_log_likelihood_gain, 0.25)
    ),
    upper_quartile_best_alternative_log_likelihood_gain = unname(
      stats::quantile(table$best_alternative_log_likelihood_gain, 0.75)
    ),
    units_with_alternative_beating_null = sum(table$alternative_beats_null),
    n_units = nrow(table),
    stringsAsFactors = FALSE
  )
}))

method_pairs <- utils::combn(method_ids, 2L, simplify = FALSE)
unit_pairwise_metrics <- do.call(rbind, lapply(method_pairs, function(pair) {
  reference <- centered_likelihoods[[pair[1L]]]
  comparison <- centered_likelihoods[[pair[2L]]]
  rows <- lapply(seq_len(nrow(reference)), function(unit_index) {
    difference <- comparison[unit_index, ] - reference[unit_index, ]
    reference_profile <- softmax(reference[unit_index, ])
    comparison_profile <- softmax(comparison[unit_index, ])
    data.frame(
      pair_key = selected_metadata$pair_key[unit_index],
      draw = selected_metadata$draw[unit_index],
      unit_in_draw = selected_metadata$unit_in_draw[unit_index],
      reference_method_id = pair[1L],
      comparison_method_id = pair[2L],
      comparison = paste(
        unname(method_labels[pair[1L]]),
        "vs",
        unname(method_labels[pair[2L]])
      ),
      pearson_correlation = safe_vector_correlation(
        reference[unit_index, ],
        comparison[unit_index, ]
      ),
      mean_absolute_difference = mean(abs(difference)),
      root_mean_squared_difference = sqrt(mean(difference^2)),
      maximum_absolute_difference = max(abs(difference)),
      normalized_profile_total_variation =
        0.5 * sum(abs(reference_profile - comparison_profile)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}))

pairwise_summary <- do.call(rbind, lapply(
  split(unit_pairwise_metrics, unit_pairwise_metrics$comparison),
  function(table) {
    data.frame(
      reference_method_id = table$reference_method_id[1L],
      comparison_method_id = table$comparison_method_id[1L],
      comparison = table$comparison[1L],
      median_pearson_correlation = stats::median(
        table$pearson_correlation,
        na.rm = TRUE
      ),
      median_root_mean_squared_difference = stats::median(
        table$root_mean_squared_difference
      ),
      lower_quartile_root_mean_squared_difference = unname(stats::quantile(
        table$root_mean_squared_difference,
        0.25
      )),
      upper_quartile_root_mean_squared_difference = unname(stats::quantile(
        table$root_mean_squared_difference,
        0.75
      )),
      maximum_root_mean_squared_difference = max(
        table$root_mean_squared_difference
      ),
      median_normalized_profile_total_variation = stats::median(
        table$normalized_profile_total_variation
      ),
      maximum_normalized_profile_total_variation = max(
        table$normalized_profile_total_variation
      ),
      stringsAsFactors = FALSE
    )
  }
))
rownames(pairwise_summary) <- NULL
pairwise_summary <- pairwise_summary[order(
  match(pairwise_summary$reference_method_id, method_ids),
  match(pairwise_summary$comparison_method_id, method_ids)
), ]

configuration <- list(
  experiment = paste(
    "Tiny fixed-subset comparison of row-null-centered FASH likelihood",
    "vectors under diagonal, C1, C2, and within-unit-centered C3 errors"
  ),
  sampling_seed = sampling_seed,
  selection_mode = selection_mode,
  selection_description = selection_description,
  n_draws = n_draws,
  units_per_draw = units_per_draw,
  n_units = nrow(selected_metadata),
  thinning_seed = fit_bundle$configuration$thinning_seed,
  C3_estimation_set = paste(
    "The 3,436 Screen-1 units used for C1; inverse-variance-weighted",
    "within-unit constant removal on the beta-hat scale, followed by",
    "SE standardization and empirical correlation"
  ),
  C3_diagnostics = C3_diagnostics[setdiff(
    names(C3_diagnostics),
    "matrix"
  )],
  psd_grid = psd_grid,
  time_grid = time_grid,
  C3_elapsed_seconds = C3_elapsed_seconds,
  generated_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE)
)

dir.create(output_directory, recursive = FALSE)
saveRDS(configuration, file.path(output_directory, "configuration.rds"))
saveRDS(C3, file.path(output_directory, "C3_correlation.rds"))
saveRDS(C3_result$fit, file.path(output_directory, "C3_tiny_fit.rds"))
utils::write.csv(
  selected_metadata,
  file.path(output_directory, "selected_units.csv"),
  row.names = FALSE
)
utils::write.csv(
  likelihood_long,
  file.path(output_directory, "likelihood_vectors_long.csv"),
  row.names = FALSE
)
utils::write.csv(
  method_unit_summary,
  file.path(output_directory, "method_unit_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  method_profile_summary,
  file.path(output_directory, "method_profile_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  unit_pairwise_metrics,
  file.path(output_directory, "unit_pairwise_metrics.csv"),
  row.names = FALSE
)
utils::write.csv(
  pairwise_summary,
  file.path(output_directory, "pairwise_summary.csv"),
  row.names = FALSE
)

method_colors <- c(
  diagonal = "#0072B2",
  C1 = "#D55E00",
  C2 = "#009E73",
  C3 = "#CC79A7"
)
likelihood_long$method_label <- factor(
  likelihood_long$method_label,
  levels = unname(method_labels)
)
likelihood_long$display_label <- factor(
  likelihood_long$display_label,
  levels = selected_metadata$display_label
)
curve_plot <- ggplot2::ggplot(
  likelihood_long,
  ggplot2::aes(
    x = prior_sd,
    y = centered_log_likelihood,
    color = method_id,
    group = method_id
  )
) +
  ggplot2::geom_hline(yintercept = 0, color = "#777777", linewidth = 0.3) +
  ggplot2::geom_line(linewidth = 0.55, alpha = 0.9) +
  ggplot2::facet_wrap(~display_label, ncol = 4, scales = "free_y") +
  ggplot2::scale_color_manual(
    values = method_colors,
    breaks = method_ids,
    labels = unname(method_labels)
  ) +
  ggplot2::labs(
    x = "Prior standard-deviation grid",
    y = "Log likelihood minus null-component log likelihood",
    color = NULL
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size = 7.5),
    legend.position = "bottom"
  )
ggplot2::ggsave(
  file.path(output_directory, "likelihood_vector_curves.png"),
  curve_plot,
  width = 16,
  height = 18,
  dpi = 180
)

metric_plot <- ggplot2::ggplot(
  unit_pairwise_metrics,
  ggplot2::aes(
    x = comparison,
    y = root_mean_squared_difference,
    color = factor(draw)
  )
) +
  ggplot2::geom_boxplot(
    ggplot2::aes(group = comparison),
    outlier.shape = NA,
    color = "#555555",
    fill = "white"
  ) +
  ggplot2::geom_point(
    position = ggplot2::position_jitter(width = 0.12, height = 0),
    alpha = 0.8,
    size = 1.8
  ) +
  ggplot2::scale_color_brewer(
    palette = "Dark2",
    name = display_group_legend
  ) +
  ggplot2::labs(
    x = NULL,
    y = "RMSE between row-null-centered log-likelihood vectors"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    legend.position = "bottom"
  )
ggplot2::ggsave(
  file.path(output_directory, "pairwise_likelihood_rmse.png"),
  metric_plot,
  width = 11,
  height = 6.8,
  dpi = 180
)

cat("\nTiny likelihood-vector comparison completed.\n")
cat("C3 elapsed seconds:", C3_elapsed_seconds, "\n")
cat("C3 minimum eigenvalue:", C3_diagnostics$minimum_eigenvalue, "\n\n")
print(method_profile_summary)
cat("\n")
print(pairwise_summary)
