#!/usr/bin/env Rscript

# Run the three-screen mashr-style null-correlation comparison.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  arguments <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(arguments, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(arguments[equals_hit[1]], nchar(equals_prefix) + 1L))
  }
  hit <- which(arguments == name)
  if (length(hit) == 0L || hit[1] == length(arguments)) {
    return(default)
  }
  arguments[hit[1] + 1L]
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

matrix_to_long <- function(x, value_name = "correlation") {
  x <- as.matrix(x)
  grid <- expand.grid(
    time_a = seq_len(nrow(x)) - 1L,
    time_b = seq_len(ncol(x)) - 1L
  )
  grid[[value_name]] <- x[cbind(grid$time_a + 1L, grid$time_b + 1L)]
  grid
}

sample_gallery <- function(candidate_metadata,
                           keep,
                           filter_id,
                           filter_label,
                           sample_size,
                           seed) {
  eligible <- which(keep)
  if (length(eligible) < sample_size) {
    stop("The selected set is too small for the requested gallery.")
  }
  set.seed(seed)
  sampled_rows <- sample(eligible, size = sample_size, replace = FALSE)
  gallery <- candidate_metadata[sampled_rows, , drop = FALSE]
  gallery$filter_id <- filter_id
  gallery$filter_label <- filter_label
  gallery$gallery_seed <- as.integer(seed)
  gallery$gallery_position <- seq_len(nrow(gallery))
  gallery$candidate_row <- sampled_rows
  rownames(gallery) <- NULL
  gallery
}

package_version_or_na <- function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    NA_character_
  }
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
  "covariance_estimation",
  "zero_intercept_correlation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "covariance_estimation",
  "mashr_screen_helpers.R"
))

thinning_seed <- as.integer(get_arg("--thinning-seed", "20260811"))
max_threshold <- as.numeric(get_arg("--max-threshold", "2"))
mean_z_threshold <- as.numeric(get_arg("--mean-z-threshold", "2"))
mashr_pair_lfdr_threshold <- as.numeric(get_arg(
  "--mashr-pair-lfdr-threshold",
  "0.05"
))
mashr_seed <- as.integer(get_arg("--mashr-seed", "123"))
bootstrap_reps <- as.integer(get_arg("--bootstrap-reps", "1000"))
bootstrap_seed_start <- as.integer(get_arg("--bootstrap-seed-start", "20260831"))
gallery_size <- as.integer(get_arg("--gallery-size", "25"))
gallery_seed_start <- as.integer(get_arg("--gallery-seed-start", "20260841"))
output_id <- get_arg("--output-id", "mashr_mean_z_null_correlation")

if (anyNA(c(
  thinning_seed,
  max_threshold,
  mean_z_threshold,
  mashr_pair_lfdr_threshold,
  mashr_seed,
  bootstrap_reps,
  bootstrap_seed_start,
  gallery_size,
  gallery_seed_start
)) || max_threshold <= 0 || mean_z_threshold <= 0 ||
    mashr_pair_lfdr_threshold <= 0 || mashr_pair_lfdr_threshold >= 1 ||
    bootstrap_reps < 20L || gallery_size < 1L || !nzchar(output_id)) {
  stop("Invalid analysis arguments.")
}

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  output_id
)
summary_dir <- file.path(output_dir, "summary")
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

fit_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "fash_fit1_update.RData"
)
message("Loading the BF-adjusted real-data FASH fit.")
fash_fit <- load_bf_adjusted_real_fash(fit_path)
pair_keys <- names(fash_fit$fash_data$data_list)
gene_index <- make_gene_index(pair_keys)

message("Drawing one variant per gene before inspecting z scores.")
selection <- select_random_variant_per_gene(
  pair_keys,
  seed = thinning_seed,
  gene_index = gene_index
)
extracted <- extract_z_matrix(fash_fit, selection$fash_index)
if (!isTRUE(all.equal(extracted$time_grid, 0:15)) ||
    !identical(rownames(extracted$z), selection$pair_key)) {
  stop("The selected data do not have the expected pair order or time grid.")
}

message("Applying the maximum-z and maximum-plus-mean-z screens.")
screen <- filter_zero_intercept_z_with_mean(
  extracted$z,
  max_threshold = max_threshold,
  mean_z_threshold = mean_z_threshold
)
if (any(screen$keep & !screen$keep_maximum)) {
  stop("The combined selected set is not nested within the maximum-z set.")
}

message("Fitting mashr to all thinned pairs with a mashr-estimated null correlation.")
mashr_screen <- fit_mashr_pair_screen(
  extracted$beta_hat,
  extracted$adjusted_se,
  z_threshold = max_threshold,
  seed = mashr_seed,
  verbose = TRUE
)
mashr_screen3_keep <- screen$keep_maximum &
  mashr_screen$pair_lfdr > mashr_pair_lfdr_threshold
if (any(mashr_screen3_keep & !screen$keep_maximum)) {
  stop("Mashr Screen 3 is not nested within the maximum-z screen.")
}

filter_definitions <- data.frame(
  filter_id = c(
    "maximum_z",
    "maximum_z_and_mean_z",
    "maximum_z_and_mashr_lfdr"
  ),
  filter_label = c(
    "max|z(t)| < 2",
    "max|z(t)| < 2 and |Z_mean| < 2",
    "max|z(t)| < 2 and mashr pair lfdr > 0.05"
  ),
  stringsAsFactors = FALSE
)
filter_keep <- list(
  maximum_z = screen$keep_maximum,
  maximum_z_and_mean_z = screen$keep,
  maximum_z_and_mashr_lfdr = mashr_screen3_keep
)

candidate_metadata <- selection
candidate_metadata$lfdr <- as.numeric(fash_fit$lfdr[selection$fash_index])
candidate_metadata$max_absolute_z <- screen$max_absolute_z
candidate_metadata$mean_z_score <- screen$mean_z_score
candidate_metadata$passes_maximum_z <- screen$keep_maximum
candidate_metadata$passes_mean_z <- screen$keep_mean_z
candidate_metadata$passes_maximum_z_and_mean_z <- screen$keep
candidate_metadata$mashr_pair_lfdr <- mashr_screen$pair_lfdr
candidate_metadata$mashr_min_condition_lfdr <- apply(
  mashr_screen$condition_lfdr,
  1L,
  min
)
candidate_metadata$mashr_max_condition_lfdr <- apply(
  mashr_screen$condition_lfdr,
  1L,
  max
)
candidate_metadata$passes_mashr_pair_lfdr <-
  mashr_screen$pair_lfdr > mashr_pair_lfdr_threshold
candidate_metadata$passes_maximum_z_and_mashr_lfdr <- mashr_screen3_keep
candidate_metadata$candidate_row <- seq_len(nrow(candidate_metadata))

estimates <- list()
selection_counts <- list()
correlation_matrices_long <- list()
adjacent_correlations <- list()
lag_summaries <- list()
bootstrap_intervals <- list()
bootstrap_seeds_used <- bootstrap_seed_start + seq_len(nrow(filter_definitions)) - 1L
if (identical(
  filter_keep$maximum_z,
  filter_keep$maximum_z_and_mashr_lfdr
)) {
  bootstrap_seeds_used[3] <- bootstrap_seeds_used[1]
}

for (filter_index in seq_len(nrow(filter_definitions))) {
  filter_id <- filter_definitions$filter_id[filter_index]
  filter_label <- filter_definitions$filter_label[filter_index]
  keep <- filter_keep[[filter_id]]
  z_selected <- extracted$z[keep, , drop = FALSE]
  message("Estimating correlation for ", filter_label, ".")
  estimate <- estimate_zero_intercept_correlation(z_selected)
  correlation <- estimate$sample_correlation
  estimates[[filter_id]] <- list(
    filter_id = filter_id,
    filter_label = filter_label,
    selected_candidate_rows = which(keep),
    selected_pair_keys = candidate_metadata$pair_key[keep],
    z_selected = z_selected,
    estimate = estimate
  )

  selection_counts[[filter_id]] <- data.frame(
    filter_id = filter_id,
    filter_label = filter_label,
    n_candidates = nrow(extracted$z),
    n_selected = sum(keep),
    selected_fraction = mean(keep),
    maximum_absolute_column_mean = max(abs(colMeans(z_selected))),
    maximum_absolute_mean_z_score = max(abs(screen$mean_z_score[keep])),
    minimum_mashr_pair_lfdr = min(mashr_screen$pair_lfdr[keep]),
    maximum_mashr_pair_lfdr = max(mashr_screen$pair_lfdr[keep]),
    stringsAsFactors = FALSE
  )

  matrix_long <- matrix_to_long(correlation)
  matrix_long$filter_id <- filter_id
  matrix_long$filter_label <- filter_label
  correlation_matrices_long[[filter_id]] <- matrix_long[, c(
    "filter_id", "filter_label", "time_a", "time_b", "correlation"
  )]

  adjacent <- correlation[cbind(1:15, 2:16)]
  adjacent_correlations[[filter_id]] <- data.frame(
    filter_id = filter_id,
    filter_label = filter_label,
    time_a = 0:14,
    time_b = 1:15,
    correlation = as.numeric(adjacent),
    stringsAsFactors = FALSE
  )

  lag_correlation <- lag_average_correlation(correlation)
  lag_summaries[[filter_id]] <- data.frame(
    filter_id = filter_id,
    filter_label = filter_label,
    lag = 1:15,
    mean_correlation = as.numeric(lag_correlation),
    semivariogram = 1 - as.numeric(lag_correlation),
    stringsAsFactors = FALSE
  )

  message("Bootstrapping selected genes for ", filter_label, ".")
  bootstrap <- bootstrap_zero_intercept_correlations(
    z_selected,
    n_bootstrap = bootstrap_reps,
    seed = bootstrap_seeds_used[filter_index]
  )
  bootstrap$filter_id <- filter_id
  bootstrap$filter_label <- filter_label
  bootstrap$point_estimate <- c(adjacent, lag_correlation)
  bootstrap$semivariogram_point <- 1 - bootstrap$point_estimate
  bootstrap$semivariogram_lower <- 1 - bootstrap$upper
  bootstrap$semivariogram_upper <- 1 - bootstrap$lower
  bootstrap_intervals[[filter_id]] <- bootstrap[, c(
    "filter_id", "filter_label", "summary_type", "index", "label",
    "point_estimate", "mean", "median", "lower", "upper",
    "semivariogram_point", "semivariogram_lower", "semivariogram_upper"
  )]
}

selection_counts <- do.call(rbind, selection_counts)
correlation_matrices_long <- do.call(rbind, correlation_matrices_long)
adjacent_correlations <- do.call(rbind, adjacent_correlations)
lag_summaries <- do.call(rbind, lag_summaries)
bootstrap_intervals <- do.call(rbind, bootstrap_intervals)
rownames(selection_counts) <- NULL
rownames(correlation_matrices_long) <- NULL
rownames(adjacent_correlations) <- NULL
rownames(lag_summaries) <- NULL
rownames(bootstrap_intervals) <- NULL

maximum_correlation <- estimates$maximum_z$estimate$sample_correlation
combined_correlation <- estimates$maximum_z_and_mean_z$estimate$sample_correlation
mashr_correlation <- estimates$maximum_z_and_mashr_lfdr$estimate$sample_correlation
correlation_list <- list(
  maximum_z = maximum_correlation,
  maximum_z_and_mean_z = combined_correlation,
  maximum_z_and_mashr_lfdr = mashr_correlation
)
comparison_specs <- data.frame(
  reference_filter_id = c(
    "maximum_z",
    "maximum_z",
    "maximum_z_and_mean_z"
  ),
  comparison_filter_id = c(
    "maximum_z_and_mean_z",
    "maximum_z_and_mashr_lfdr",
    "maximum_z_and_mashr_lfdr"
  ),
  stringsAsFactors = FALSE
)
matrix_comparison_rows <- list()
correlation_difference_rows <- list()
for (comparison_index in seq_len(nrow(comparison_specs))) {
  reference_id <- comparison_specs$reference_filter_id[comparison_index]
  comparison_id <- comparison_specs$comparison_filter_id[comparison_index]
  reference_matrix <- correlation_list[[reference_id]]
  comparison_matrix <- correlation_list[[comparison_id]]
  difference <- comparison_matrix - reference_matrix
  off_diagonal <- upper.tri(difference)
  comparison_label <- paste0(comparison_id, " minus ", reference_id)
  matrix_comparison_rows[[comparison_index]] <- data.frame(
    comparison = comparison_label,
    reference_filter_id = reference_id,
    comparison_filter_id = comparison_id,
    off_diagonal_matrix_correlation = stats::cor(
      reference_matrix[off_diagonal],
      comparison_matrix[off_diagonal]
    ),
    mean_off_diagonal_difference = mean(difference[off_diagonal]),
    mean_absolute_off_diagonal_difference = mean(abs(difference[off_diagonal])),
    maximum_absolute_off_diagonal_difference = max(abs(
      difference[off_diagonal]
    )),
    full_frobenius_difference = sqrt(sum(difference^2)),
    stringsAsFactors = FALSE
  )
  difference_long <- matrix_to_long(
    difference,
    value_name = "correlation_difference"
  )
  difference_long$comparison <- comparison_label
  difference_long$reference_filter_id <- reference_id
  difference_long$comparison_filter_id <- comparison_id
  correlation_difference_rows[[comparison_index]] <- difference_long[, c(
    "comparison", "reference_filter_id", "comparison_filter_id",
    "time_a", "time_b", "correlation_difference"
  )]
}
matrix_comparison <- do.call(rbind, matrix_comparison_rows)
correlation_difference_long <- do.call(rbind, correlation_difference_rows)
rownames(matrix_comparison) <- NULL
rownames(correlation_difference_long) <- NULL
key_lags <- lag_summaries[lag_summaries$lag %in% c(1L, 5L, 10L, 15L), ]

message("Selecting fixed trajectory galleries.")
gallery_membership <- rbind(
  sample_gallery(
    candidate_metadata,
    filter_keep$maximum_z,
    filter_id = "maximum_z",
    filter_label = filter_definitions$filter_label[1],
    sample_size = gallery_size,
    seed = gallery_seed_start
  ),
  sample_gallery(
    candidate_metadata,
    filter_keep$maximum_z_and_mean_z,
    filter_id = "maximum_z_and_mean_z",
    filter_label = filter_definitions$filter_label[2],
    sample_size = gallery_size,
    seed = gallery_seed_start + 1L
  ),
  sample_gallery(
    candidate_metadata,
    filter_keep$maximum_z_and_mashr_lfdr,
    filter_id = "maximum_z_and_mashr_lfdr",
    filter_label = filter_definitions$filter_label[3],
    sample_size = gallery_size,
    seed = gallery_seed_start + 2L
  )
)
gallery_overlap_rows <- list()
gallery_combinations <- utils::combn(filter_definitions$filter_id, 2L)
for (overlap_index in seq_len(ncol(gallery_combinations))) {
  first_filter <- gallery_combinations[1L, overlap_index]
  second_filter <- gallery_combinations[2L, overlap_index]
  gallery_overlap_rows[[overlap_index]] <- data.frame(
    first_filter_id = first_filter,
    second_filter_id = second_filter,
    n_overlap = length(intersect(
      gallery_membership$pair_key[gallery_membership$filter_id == first_filter],
      gallery_membership$pair_key[gallery_membership$filter_id == second_filter]
    )),
    stringsAsFactors = FALSE
  )
}
gallery_overlaps <- do.call(rbind, gallery_overlap_rows)

mashr_lfdr_probabilities <- c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)
mashr_lfdr_quantiles <- data.frame(
  population = c(
    rep("All thinned pairs", length(mashr_lfdr_probabilities)),
    rep("Screen 1 pairs", length(mashr_lfdr_probabilities))
  ),
  probability = rep(mashr_lfdr_probabilities, 2L),
  pair_lfdr = c(
    as.numeric(stats::quantile(
      mashr_screen$pair_lfdr,
      probs = mashr_lfdr_probabilities,
      names = FALSE
    )),
    as.numeric(stats::quantile(
      mashr_screen$pair_lfdr[screen$keep_maximum],
      probs = mashr_lfdr_probabilities,
      names = FALSE
    ))
  ),
  stringsAsFactors = FALSE
)
mashr_screen_summary <- data.frame(
  n_all_thinned_pairs = nrow(extracted$z),
  n_screen1_pairs = sum(screen$keep_maximum),
  n_pair_lfdr_at_or_below_threshold = sum(
    mashr_screen$pair_lfdr <= mashr_pair_lfdr_threshold
  ),
  n_screen1_pair_lfdr_at_or_below_threshold = sum(
    screen$keep_maximum &
      mashr_screen$pair_lfdr <= mashr_pair_lfdr_threshold
  ),
  fitted_pi0 = mashr_screen$fitted_pi0,
  minimum_screen1_pair_lfdr = min(
    mashr_screen$pair_lfdr[screen$keep_maximum]
  ),
  null_correlation_maximum_difference =
    mashr_screen$diagnostics$null_correlation_maximum_difference,
  stringsAsFactors = FALSE
)

bad_example_key <- "ENSG00000078618_rs7544633"
bad_example_index <- match(bad_example_key, candidate_metadata$pair_key)
if (is.na(bad_example_index)) {
  stop("The prespecified constant-eQTL-like gallery example is missing.")
}
bad_beta <- extracted$beta_hat[bad_example_index, ]
bad_se <- extracted$adjusted_se[bad_example_index, ]
bad_weights <- 1 / bad_se^2
bad_example_audit <- data.frame(
  pair_key = bad_example_key,
  gene_id = candidate_metadata$gene_id[bad_example_index],
  variant_id = candidate_metadata$variant_id[bad_example_index],
  fash_lfdr = candidate_metadata$lfdr[bad_example_index],
  minimum_beta_hat = min(bad_beta),
  maximum_beta_hat = max(bad_beta),
  inverse_variance_weighted_constant_beta =
    sum(bad_weights * bad_beta) / sum(bad_weights),
  maximum_absolute_z = screen$max_absolute_z[bad_example_index],
  mean_z_score = screen$mean_z_score[bad_example_index],
  mashr_pair_lfdr = mashr_screen$pair_lfdr[bad_example_index],
  mashr_min_condition_lfdr = min(
    mashr_screen$condition_lfdr[bad_example_index, ]
  ),
  passes_screen1 = screen$keep_maximum[bad_example_index],
  passes_screen2 = screen$keep[bad_example_index],
  passes_screen3 = mashr_screen3_keep[bad_example_index],
  stringsAsFactors = FALSE
)

fit_info <- file.info(fit_path)
configuration <- list(
  analysis_id = output_id,
  created_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  fit_path = normalizePath(fit_path, mustWork = TRUE),
  fit_size = unname(fit_info$size),
  fit_mtime = as.character(fit_info$mtime),
  n_fash_pairs = length(pair_keys),
  n_genes = length(gene_index),
  n_time = ncol(extracted$z),
  time_grid = extracted$time_grid,
  thinning_seed = thinning_seed,
  max_threshold = max_threshold,
  mean_z_threshold = mean_z_threshold,
  mean_z_definition = "sqrt(n_time) * rowMeans(z)",
  mashr_pair_lfdr_threshold = mashr_pair_lfdr_threshold,
  mashr_seed = mashr_seed,
  mashr_null_correlation_function =
    "mashr::estimate_null_correlation_simple",
  mashr_null_correlation_z_threshold = max_threshold,
  mashr_null_correlation_rule =
    "cor(z) among rows with max(abs(z)) < z_threshold",
  mashr_cov_methods = mashr_screen$diagnostics$cov_methods,
  mashr_prior = mashr_screen$diagnostics$prior,
  mashr_nullweight = mashr_screen$diagnostics$nullweight,
  mashr_optmethod = mashr_screen$diagnostics$optmethod,
  mashr_pair_lfdr_definition =
    "posterior weight on the exact all-zero component",
  mashr_fitted_pi0 = mashr_screen$fitted_pi0,
  mashr_n_nullish_for_correlation =
    mashr_screen$diagnostics$n_nullish_for_correlation,
  mashr_null_correlation_maximum_difference =
    mashr_screen$diagnostics$null_correlation_maximum_difference,
  mashr_fit_filename = "mashr_screen_fit.rds",
  estimator = "stats::cor(z_selected)",
  bootstrap_reps = bootstrap_reps,
  bootstrap_seeds = bootstrap_seeds_used,
  gallery_size = gallery_size,
  gallery_seeds = gallery_seed_start + 0:2,
  gallery_overlaps = gallery_overlaps,
  r_version = R.version.string,
  package_versions = c(
    fashr = package_version_or_na("fashr"),
    mashr = package_version_or_na("mashr"),
    workflowr = package_version_or_na("workflowr"),
    knitr = package_version_or_na("knitr")
  )
)

analysis <- list(
  configuration = configuration,
  filter_definitions = filter_definitions,
  candidate_metadata = candidate_metadata,
  candidate_z = extracted$z,
  candidate_beta_hat = extracted$beta_hat,
  candidate_adjusted_se = extracted$adjusted_se,
  screen = screen,
  estimates = estimates,
  selection_counts = selection_counts,
  correlation_matrices_long = correlation_matrices_long,
  adjacent_correlations = adjacent_correlations,
  lag_summaries = lag_summaries,
  bootstrap_intervals = bootstrap_intervals,
  correlation_difference_long = correlation_difference_long,
  matrix_comparison = matrix_comparison,
  key_lags = key_lags,
  gallery_membership = gallery_membership,
  gallery_overlaps = gallery_overlaps,
  mashr_null_correlation = mashr_screen$null_correlation,
  mashr_screen_summary = mashr_screen_summary,
  mashr_lfdr_quantiles = mashr_lfdr_quantiles,
  bad_example_audit = bad_example_audit
)

mashr_fit_cache <- list(
  fit = mashr_screen$fit,
  null_correlation = mashr_screen$null_correlation,
  nullish_keep = mashr_screen$nullish_keep,
  pair_lfdr = mashr_screen$pair_lfdr,
  condition_lfdr = mashr_screen$condition_lfdr,
  fitted_pi0 = mashr_screen$fitted_pi0,
  canonical_covariance_names = mashr_screen$canonical_covariance_names,
  diagnostics = mashr_screen$diagnostics
)
mashr_null_correlation_long <- matrix_to_long(
  mashr_screen$null_correlation,
  value_name = "correlation"
)
mashr_fitted_mixture <- data.frame(
  component = names(mashr_screen$fit$fitted_g$pi),
  fitted_weight = as.numeric(mashr_screen$fit$fitted_g$pi),
  stringsAsFactors = FALSE
)

saveRDS(configuration, file.path(output_dir, "configuration.rds"))
saveRDS(analysis, file.path(output_dir, "mashr_mean_z_null_correlation.rds"))
saveRDS(mashr_fit_cache, file.path(output_dir, "mashr_screen_fit.rds"))
write_csv(selection_counts, file.path(summary_dir, "selection_counts.csv"))
write_csv(candidate_metadata, file.path(summary_dir, "candidate_metadata.csv"))
write_csv(
  candidate_metadata[candidate_metadata$passes_maximum_z, ],
  file.path(summary_dir, "selected_maximum_z_units.csv")
)
write_csv(
  candidate_metadata[candidate_metadata$passes_maximum_z_and_mean_z, ],
  file.path(summary_dir, "selected_maximum_z_and_mean_z_units.csv")
)
write_csv(
  candidate_metadata[
    candidate_metadata$passes_maximum_z_and_mashr_lfdr,
  ],
  file.path(summary_dir, "selected_maximum_z_and_mashr_lfdr_units.csv")
)
write_csv(
  correlation_matrices_long,
  file.path(summary_dir, "correlation_matrices_long.csv")
)
write_csv(
  adjacent_correlations,
  file.path(summary_dir, "adjacent_correlations.csv")
)
write_csv(lag_summaries, file.path(summary_dir, "lag_summaries.csv"))
write_csv(
  bootstrap_intervals,
  file.path(summary_dir, "bootstrap_intervals.csv")
)
write_csv(
  correlation_difference_long,
  file.path(summary_dir, "correlation_difference_long.csv")
)
write_csv(matrix_comparison, file.path(summary_dir, "matrix_comparison.csv"))
write_csv(key_lags, file.path(summary_dir, "key_lags.csv"))
write_csv(
  gallery_membership,
  file.path(summary_dir, "gallery_membership.csv")
)
write_csv(gallery_overlaps, file.path(summary_dir, "gallery_overlaps.csv"))
write_csv(
  mashr_null_correlation_long,
  file.path(summary_dir, "mashr_null_correlation_long.csv")
)
write_csv(
  mashr_screen_summary,
  file.path(summary_dir, "mashr_screen_summary.csv")
)
write_csv(
  mashr_lfdr_quantiles,
  file.path(summary_dir, "mashr_lfdr_quantiles.csv")
)
write_csv(
  mashr_fitted_mixture,
  file.path(summary_dir, "mashr_fitted_mixture.csv")
)
write_csv(
  bad_example_audit,
  file.path(summary_dir, "bad_example_audit.csv")
)

message(
  "Completed. Retained ",
  selection_counts$n_selected[selection_counts$filter_id == "maximum_z"],
  " units under max|z| < 2 and ",
  selection_counts$n_selected[
    selection_counts$filter_id == "maximum_z_and_mean_z"
  ],
  " units after adding |Z_mean| < 2, and ",
  selection_counts$n_selected[
    selection_counts$filter_id == "maximum_z_and_mashr_lfdr"
  ],
  " units after adding mashr pair lfdr > ",
  mashr_pair_lfdr_threshold,
  "."
)
