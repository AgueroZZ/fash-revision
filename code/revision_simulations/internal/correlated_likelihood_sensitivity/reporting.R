# Load, validate, summarize, and plot the correlated-likelihood sensitivity cache.

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required correlated-likelihood cache file is missing: ", path)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

format_decimal <- function(x, digits = 4) {
  ifelse(
    is.finite(x),
    formatC(x, format = "f", digits = digits),
    "Unavailable"
  )
}

format_integer <- function(x) {
  ifelse(
    is.finite(x),
    format(as.integer(x), big.mark = ",", scientific = FALSE),
    "Unavailable"
  )
}

format_scientific <- function(x, digits = 2) {
  formatC(x, format = "e", digits = digits)
}

render_scrollable_table <- function(table,
                                    align = NULL,
                                    minimum_width = "760px",
                                    digits = NULL) {
  arguments <- list(
    x = table,
    format = "html",
    align = align,
    escape = TRUE,
    table.attr = paste0(
      'class="table table-striped table-hover" style="min-width:',
      minimum_width,
      ';"'
    )
  )
  if (!is.null(digits)) {
    arguments$digits <- digits
  }
  rendered <- do.call(knitr::kable, arguments)
  cat(
    '<div class="correlated-likelihood-table-scroll">\n',
    as.character(rendered),
    '\n</div>\n',
    sep = ""
  )
  invisible(rendered)
}

assert_data_frame_equal <- function(observed, expected, label) {
  rownames(observed) <- NULL
  rownames(expected) <- NULL
  structural_match <-
    identical(dim(observed), dim(expected)) &&
    identical(names(observed), names(expected))
  column_match <- structural_match && all(vapply(
    names(expected),
    function(column_name) {
      observed_values <- observed[[column_name]]
      expected_values <- expected[[column_name]]
      if ((is.numeric(observed_values) || is.logical(observed_values)) &&
          (is.numeric(expected_values) || is.logical(expected_values))) {
        return(isTRUE(all.equal(
          as.numeric(observed_values),
          as.numeric(expected_values),
          tolerance = 1e-12,
          check.attributes = FALSE
        )))
      }
      identical(observed_values, expected_values)
    },
    logical(1)
  ))
  if (!isTRUE(column_match)) {
    stop(label, " does not match the retained analysis cache.")
  }
}

expected_analysis_id <- "correlated_likelihood_sensitivity"
expected_seed <- 20260811L
expected_n_units <- 6362L
expected_n_time <- 16L
expected_alpha <- 0.05
expected_method_ids <- c("diagonal", "C1", "C2")
expected_method_labels <- c(
  diagonal = "Diagonal SE",
  C1 = "C1: Screen 1",
  C2 = "C2: Screen 2"
)
expected_stages <- c("Raw", "BF-adjusted")
expected_available_method_stage <- data.frame(
  fit_stage = c(rep("Raw", 3L), rep("BF-adjusted", 2L)),
  method_id = c("diagonal", "C1", "C2", "diagonal", "C1"),
  stringsAsFactors = FALSE
)
expected_source_sha256 <- c(
  covariance_configuration =
    "1c297e951db9a10444605bd99445d242771e5fc2e753fb8bcb04bcf8830a2a91",
  covariance_analysis =
    "85259f63fd82ac119477e75c0cb9e1c1815169b79ccfb6f3df6ecab8842b13db",
  raw_fash_fit =
    "9e9aaf8a405f7ca83439990656666035fc0a74bb1ae2858ace79944fac6ec929",
  bf_adjusted_fash_fit =
    "3e3d6735b9da734a3ab16ed53713908c02d0af4e15a659760fc5f892ef64023b"
)

cache_directory <- file.path(
  "output",
  "revision_simulations",
  "internal",
  expected_analysis_id
)
summary_directory <- file.path(cache_directory, "summary")
required_rds_files <- c("analysis.rds", "configuration.rds", "fit_bundle.rds")
required_csv_files <- c(
  "bf_rebuild_validation.csv",
  "bf_update_status.csv",
  "correlation_diagnostics.csv",
  "correlation_matrices_long.csv",
  "discovery_summary.csv",
  "identity_path_validation.csv",
  "lfdr_pairwise_metrics.csv",
  "method_stage_summary.csv",
  "pair_metadata.csv",
  "prior_pairwise_metrics.csv",
  "prior_weights.csv",
  "top_lfdr_discrepancies.csv",
  "unit_lfdr_long.csv",
  "unit_lfdr_wide.csv"
)
required_paths <- c(
  file.path(cache_directory, required_rds_files),
  file.path(summary_directory, required_csv_files)
)
if (!dir.exists(cache_directory) || any(!file.exists(required_paths))) {
  stop("The retained correlated-likelihood sensitivity cache is incomplete.")
}

correlated_likelihood_analysis <- readRDS(file.path(
  cache_directory,
  "analysis.rds"
))
correlated_likelihood_configuration <- readRDS(file.path(
  cache_directory,
  "configuration.rds"
))

if (!isTRUE(all.equal(
  correlated_likelihood_analysis$configuration,
  correlated_likelihood_configuration,
  tolerance = 0,
  check.attributes = FALSE
))) {
  stop("The analysis and configuration caches disagree.")
}

configuration <- correlated_likelihood_configuration
configuration_valid <-
  identical(configuration$analysis_id, expected_analysis_id) &&
  identical(configuration$thinning_seed, expected_seed) &&
  identical(configuration$n_units, expected_n_units) &&
  identical(configuration$n_genes, expected_n_units) &&
  identical(configuration$n_time, expected_n_time) &&
  identical(as.numeric(configuration$time_grid), as.numeric(0:15)) &&
  isTRUE(all.equal(configuration$alpha, expected_alpha)) &&
  identical(configuration$method_ids, expected_method_ids) &&
  identical(unname(configuration$method_labels),
            unname(expected_method_labels)) &&
  identical(configuration$precision_formula,
            "diag(1 / se_j) %*% solve(C) %*% diag(1 / se_j)") &&
  identical(configuration$raw_settings$penalty, 10) &&
  identical(length(configuration$psd_grid), 52L) &&
  identical(configuration$raw_penalty, 10) &&
  identical(configuration$validation_units, 24L) &&
  is.list(configuration$source_files) &&
  identical(sort(names(configuration$source_files)),
            sort(names(expected_source_sha256)))
if (!isTRUE(configuration_valid)) {
  stop("The correlated-likelihood configuration failed validation.")
}

recorded_source_sha256 <- vapply(
  configuration$source_files,
  `[[`,
  character(1),
  "sha256"
)
recorded_source_sizes <- vapply(
  configuration$source_files,
  `[[`,
  numeric(1),
  "size_bytes"
)
recorded_source_paths <- vapply(
  configuration$source_files,
  `[[`,
  character(1),
  "path"
)
source_file_info <- file.info(recorded_source_paths)
if (!identical(
  unname(recorded_source_sha256[names(expected_source_sha256)]),
  unname(expected_source_sha256)
) || any(!file.exists(recorded_source_paths)) ||
    any(source_file_info$size != recorded_source_sizes)) {
  stop("The recorded source-file provenance failed validation.")
}

pair_metadata <- correlated_likelihood_analysis$pair_metadata
correlations <- correlated_likelihood_analysis$correlations
correlation_diagnostics <-
  correlated_likelihood_analysis$correlation_diagnostics
correlation_matrices_long <-
  correlated_likelihood_analysis$correlation_matrices_long
identity_validation <-
  correlated_likelihood_analysis$identity_validation_table
method_stage_summary <-
  correlated_likelihood_analysis$method_stage_summary
prior_weights <- correlated_likelihood_analysis$prior_weights
prior_pairwise_metrics <-
  correlated_likelihood_analysis$prior_pairwise_metrics
lfdr_wide <- correlated_likelihood_analysis$lfdr_wide
lfdr_long <- correlated_likelihood_analysis$lfdr_long
lfdr_pairwise_metrics <-
  correlated_likelihood_analysis$lfdr_pairwise_metrics
discovery_summary <- correlated_likelihood_analysis$discovery_summary
top_lfdr_discrepancies <-
  correlated_likelihood_analysis$top_lfdr_discrepancies
bf_update_status <- correlated_likelihood_analysis$bf_update_status
bf_rebuild_validation <-
  correlated_likelihood_analysis$bf_rebuild_validation

csv_object_map <- list(
  bf_rebuild_validation.csv = bf_rebuild_validation,
  bf_update_status.csv = bf_update_status,
  correlation_diagnostics.csv = correlation_diagnostics,
  correlation_matrices_long.csv = correlation_matrices_long,
  discovery_summary.csv = discovery_summary,
  identity_path_validation.csv = identity_validation,
  lfdr_pairwise_metrics.csv = lfdr_pairwise_metrics,
  method_stage_summary.csv = method_stage_summary,
  pair_metadata.csv = pair_metadata,
  prior_pairwise_metrics.csv = prior_pairwise_metrics,
  prior_weights.csv = prior_weights,
  top_lfdr_discrepancies.csv = top_lfdr_discrepancies,
  unit_lfdr_long.csv = lfdr_long,
  unit_lfdr_wide.csv = lfdr_wide
)
for (csv_name in names(csv_object_map)) {
  assert_data_frame_equal(
    read_required_csv(file.path(summary_directory, csv_name)),
    csv_object_map[[csv_name]],
    csv_name
  )
}

required_pair_columns <- c(
  "seed", "fash_index", "pair_key", "gene_id", "variant_id",
  "candidate_row"
)
pair_metadata_valid <-
  all(required_pair_columns %in% names(pair_metadata)) &&
  nrow(pair_metadata) == expected_n_units &&
  all(pair_metadata$seed == expected_seed) &&
  anyDuplicated(pair_metadata$fash_index) == 0L &&
  anyDuplicated(pair_metadata$pair_key) == 0L &&
  anyDuplicated(pair_metadata$gene_id) == 0L &&
  identical(pair_metadata$candidate_row, seq_len(expected_n_units)) &&
  identical(sub("_.*$", "", pair_metadata$pair_key),
            pair_metadata$gene_id) &&
  identical(sub("^[^_]+_", "", pair_metadata$pair_key),
            pair_metadata$variant_id)
if (!isTRUE(pair_metadata_valid)) {
  stop("The retained thinning is not exactly one unique pair per gene.")
}

if (!identical(sort(names(correlations)), c("C1", "C2"))) {
  stop("The C1/C2 correlation matrices are missing.")
}
for (matrix_id in c("C1", "C2")) {
  correlation_matrix <- correlations[[matrix_id]]
  eigenvalues <- eigen(correlation_matrix, symmetric = TRUE,
                       only.values = TRUE)$values
  if (!is.matrix(correlation_matrix) ||
      !identical(dim(correlation_matrix), c(expected_n_time, expected_n_time)) ||
      any(!is.finite(correlation_matrix)) ||
      max(abs(correlation_matrix - t(correlation_matrix))) > 1e-12 ||
      max(abs(diag(correlation_matrix) - 1)) > 1e-12 ||
      min(eigenvalues) <= 0) {
    stop(matrix_id, " is not a finite positive-definite correlation matrix.")
  }
}

expected_row_counts <- c(
  correlation_diagnostics = 2L,
  correlation_matrices_long = 512L,
  identity_validation = 1L,
  method_stage_summary = 6L,
  prior_weights = 260L,
  prior_pairwise_metrics = 4L,
  lfdr_wide = 6362L,
  lfdr_long = 31810L,
  lfdr_pairwise_metrics = 4L,
  discovery_summary = 6L,
  top_lfdr_discrepancies = 400L,
  bf_update_status = 3L,
  bf_rebuild_validation = 2L
)
row_count_objects <- list(
  correlation_diagnostics = correlation_diagnostics,
  correlation_matrices_long = correlation_matrices_long,
  identity_validation = identity_validation,
  method_stage_summary = method_stage_summary,
  prior_weights = prior_weights,
  prior_pairwise_metrics = prior_pairwise_metrics,
  lfdr_wide = lfdr_wide,
  lfdr_long = lfdr_long,
  lfdr_pairwise_metrics = lfdr_pairwise_metrics,
  discovery_summary = discovery_summary,
  top_lfdr_discrepancies = top_lfdr_discrepancies,
  bf_update_status = bf_update_status,
  bf_rebuild_validation = bf_rebuild_validation
)
observed_row_counts <- vapply(row_count_objects, nrow, integer(1))
if (!identical(unname(observed_row_counts), unname(expected_row_counts))) {
  stop("One or more retained summary tables has an unexpected row count.")
}

raw_rows <- method_stage_summary$fit_stage == "Raw"
bf_rows <- method_stage_summary$fit_stage == "BF-adjusted"
c2_bf_row <- bf_rows & method_stage_summary$method_id == "C2"
available_rows <- method_stage_summary$result_status == "Available"
method_stage_valid <-
  identical(method_stage_summary$fit_stage,
            c(rep("Raw", 3L), rep("BF-adjusted", 3L))) &&
  identical(method_stage_summary$method_id,
            rep(expected_method_ids, 2L)) &&
  all(method_stage_summary$n_units == expected_n_units) &&
  all(method_stage_summary$alpha == expected_alpha) &&
  all(is.finite(method_stage_summary$pi0[available_rows])) &&
  all(method_stage_summary$pi0[available_rows] >= 0 &
        method_stage_summary$pi0[available_rows] <= 1) &&
  sum(c2_bf_row) == 1L &&
  grepl("Unavailable", method_stage_summary$result_status[c2_bf_row],
        fixed = TRUE) &&
  all(is.na(method_stage_summary[c2_bf_row, c(
    "mean_lfdr", "median_lfdr", "discovered_units", "discovered_genes",
    "pi0"
  )]))
if (!isTRUE(method_stage_valid)) {
  stop("Method/stage availability or result values are inconsistent.")
}

available_prior_keys <- paste(
  expected_available_method_stage$fit_stage,
  expected_available_method_stage$method_id,
  sep = "::"
)
observed_prior_keys <- paste(
  prior_weights$fit_stage,
  prior_weights$method_id,
  sep = "::"
)
prior_sums <- tapply(prior_weights$prior_weight, observed_prior_keys, sum)
prior_valid <-
  identical(sort(unique(observed_prior_keys)), sort(available_prior_keys)) &&
  all(table(observed_prior_keys) == 52L) &&
  all(is.finite(prior_weights$prior_weight)) &&
  all(prior_weights$prior_weight >= -1e-12) &&
  max(abs(prior_sums - 1)) <= 1e-8
if (!isTRUE(prior_valid)) {
  stop("Available prior weights failed grid, normalization, or range checks.")
}

required_lfdr_columns <- c(
  "pair_key", "gene_id", "variant_id", "diagonal_raw", "diagonal_bf",
  "C1_raw", "C1_bf", "C2_raw", "C2_bf"
)
finite_lfdr_columns <- setdiff(required_lfdr_columns, c(
  "pair_key", "gene_id", "variant_id", "C2_bf"
))
lfdr_valid <-
  identical(names(lfdr_wide), required_lfdr_columns) &&
  identical(lfdr_wide$pair_key, pair_metadata$pair_key) &&
  all(vapply(
    lfdr_wide[finite_lfdr_columns],
    function(values) all(is.finite(values) & values >= 0 & values <= 1),
    logical(1)
  )) &&
  all(is.na(lfdr_wide$C2_bf)) &&
  identical(sort(unique(paste(lfdr_long$fit_stage, lfdr_long$method_id,
                              sep = "::"))),
            sort(available_prior_keys)) &&
  all(lfdr_long$lfdr >= 0 & lfdr_long$lfdr <= 1)
if (!isTRUE(lfdr_valid)) {
  stop("Per-unit lfdr values failed ordering, availability, or range checks.")
}

identity_valid <-
  identity_validation$n_units == 24L &&
  identity_validation$row_centered_likelihood_maximum_difference <= 1e-7 &&
  identity_validation$prior_weight_maximum_difference <= 1e-6 &&
  identity_validation$lfdr_maximum_difference <= 1e-6
bf_rebuild_valid <-
  identical(bf_rebuild_validation$method_id, c("diagonal", "C1")) &&
  all(bf_rebuild_validation$lfdr_maximum_difference <= 1e-8) &&
  all(bf_rebuild_validation$prior_weight_maximum_difference <= 1e-8)
bf_status_valid <-
  identical(bf_update_status$method_id, expected_method_ids) &&
  identical(bf_update_status$bf_update_available, c(TRUE, TRUE, FALSE)) &&
  bf_update_status$raw_alternative_prior_mass[3] == 0 &&
  grepl("conditional alternative mixture is undefined",
        bf_update_status$bf_update_status[3], fixed = TRUE)
if (!isTRUE(identity_valid && bf_rebuild_valid && bf_status_valid)) {
  stop("Identity-path or BF-finalization validation failed.")
}

method_stage_summary_display <- data.frame(
  Stage = method_stage_summary$fit_stage,
  Method = method_stage_summary$method_label,
  Status = method_stage_summary$result_status,
  `pi0` = format_decimal(method_stage_summary$pi0, 4),
  `Mean lfdr` = format_decimal(method_stage_summary$mean_lfdr, 4),
  `Median lfdr` = format_decimal(method_stage_summary$median_lfdr, 4),
  `FDR calls` = format_integer(method_stage_summary$discovered_units),
  `Elapsed seconds` = format_decimal(method_stage_summary$elapsed_seconds, 1),
  check.names = FALSE
)

correlation_diagnostics_display <- data.frame(
  Matrix = correlation_diagnostics$matrix_label,
  `Estimating pairs` = format_integer(
    correlation_diagnostics$n_estimating_pairs
  ),
  `Minimum eigenvalue` = format_decimal(
    correlation_diagnostics$minimum_eigenvalue,
    4
  ),
  `Maximum eigenvalue` = format_decimal(
    correlation_diagnostics$maximum_eigenvalue,
    4
  ),
  `Condition number` = format_decimal(
    correlation_diagnostics$eigenvalue_condition_number,
    2
  ),
  check.names = FALSE
)

identity_validation_display <- data.frame(
  `Validated units` = identity_validation$n_units,
  `Row-centered likelihood max difference` = format_scientific(
    identity_validation$row_centered_likelihood_maximum_difference
  ),
  `Prior-weight max difference` = format_scientific(
    identity_validation$prior_weight_maximum_difference
  ),
  `lfdr max difference` = format_scientific(
    identity_validation$lfdr_maximum_difference
  ),
  check.names = FALSE
)

prior_pairwise_metrics_display <- data.frame(
  Stage = prior_pairwise_metrics$fit_stage,
  Comparison = paste(
    prior_pairwise_metrics$reference_method_label,
    "vs",
    prior_pairwise_metrics$comparison_method_label
  ),
  `Reference pi0` = format_decimal(prior_pairwise_metrics$reference_pi0, 4),
  `Comparison pi0` = format_decimal(
    prior_pairwise_metrics$comparison_pi0,
    4
  ),
  `pi0 difference` = format_decimal(
    prior_pairwise_metrics$pi0_difference,
    4
  ),
  `Prior total variation` = format_decimal(
    prior_pairwise_metrics$prior_total_variation,
    4
  ),
  check.names = FALSE
)

lfdr_pairwise_metrics_display <- data.frame(
  Stage = lfdr_pairwise_metrics$fit_stage,
  Comparison = paste(
    lfdr_pairwise_metrics$reference_method_label,
    "vs",
    lfdr_pairwise_metrics$comparison_method_label
  ),
  Pearson = format_decimal(lfdr_pairwise_metrics$pearson_lfdr, 3),
  Spearman = format_decimal(lfdr_pairwise_metrics$spearman_lfdr, 3),
  `Mean signed difference` = format_decimal(
    lfdr_pairwise_metrics$mean_signed_lfdr_difference,
    4
  ),
  MAE = format_decimal(
    lfdr_pairwise_metrics$mean_absolute_lfdr_difference,
    4
  ),
  RMSE = format_decimal(lfdr_pairwise_metrics$rmse_lfdr, 4),
  `Maximum difference` = format_decimal(
    lfdr_pairwise_metrics$maximum_absolute_lfdr_difference,
    4
  ),
  Jaccard = format_decimal(lfdr_pairwise_metrics$discovery_jaccard, 3),
  check.names = FALSE
)

bf_update_status_display <- data.frame(
  Method = bf_update_status$method_label,
  `Raw alternative prior mass` = format_decimal(
    bf_update_status$raw_alternative_prior_mass,
    6
  ),
  `BF result available` = ifelse(
    bf_update_status$bf_update_available,
    "Yes",
    "No"
  ),
  Status = bf_update_status$bf_update_status,
  check.names = FALSE
)

top_discrepancy_display <- do.call(
  rbind,
  lapply(
    split(
      top_lfdr_discrepancies,
      paste(
        top_lfdr_discrepancies$fit_stage,
        top_lfdr_discrepancies$reference_method_id,
        top_lfdr_discrepancies$comparison_method_id,
        sep = "::"
      )
    ),
    function(table) table[table$discrepancy_rank <= 5L, , drop = FALSE]
  )
)
top_discrepancy_display <- data.frame(
  Stage = top_discrepancy_display$fit_stage,
  Comparison = paste(
    expected_method_labels[top_discrepancy_display$reference_method_id],
    "vs",
    expected_method_labels[top_discrepancy_display$comparison_method_id]
  ),
  Rank = top_discrepancy_display$discrepancy_rank,
  Gene = top_discrepancy_display$gene_id,
  Variant = top_discrepancy_display$variant_id,
  `Reference lfdr` = format_decimal(
    top_discrepancy_display$reference_lfdr,
    4
  ),
  `Comparison lfdr` = format_decimal(
    top_discrepancy_display$comparison_lfdr,
    4
  ),
  `Absolute difference` = format_decimal(
    top_discrepancy_display$absolute_lfdr_difference,
    4
  ),
  check.names = FALSE
)

source_provenance_display <- data.frame(
  Source = names(configuration$source_files),
  Path = basename(recorded_source_paths),
  `Size (bytes)` = format(recorded_source_sizes, big.mark = ",",
                          scientific = FALSE),
  `SHA-256` = recorded_source_sha256,
  check.names = FALSE
)

method_colors <- c(
  diagonal = "#0072B2",
  C1 = "#D55E00",
  C2 = "#009E73"
)
stage_levels <- c("Raw", "BF-adjusted")

base_plot_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      plot.title.position = "plot"
    )
}

plot_correlation_heatmaps <- function() {
  plot_data <- correlation_matrices_long
  plot_data$matrix_label <- factor(
    plot_data$matrix_label,
    levels = unname(expected_method_labels[c("C1", "C2")])
  )
  off_diagonal <- plot_data$time_a != plot_data$time_b
  color_limit <- max(abs(plot_data$correlation[off_diagonal]))
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = factor(time_b), y = factor(time_a), fill = correlation)
  ) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~matrix_label, nrow = 1) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-color_limit, color_limit),
      oob = scales::squish,
      name = "Correlation"
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "Time", y = "Time") +
    base_plot_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )
}

plot_prior_weights <- function() {
  plot_data <- prior_weights
  plot_data$fit_stage <- factor(plot_data$fit_stage, levels = stage_levels)
  plot_data$method_label <- factor(
    plot_data$method_label,
    levels = unname(expected_method_labels)
  )
  plot_data$prior_standard_deviation <- sqrt(plot_data$psd)
  unavailable_annotation <- data.frame(
    fit_stage = factor("BF-adjusted", levels = stage_levels),
    x = 0.52,
    y = max(plot_data$prior_weight) * 0.84,
    label = "C2 unavailable:\nraw alternative mass = 0"
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = prior_standard_deviation,
      y = prior_weight,
      color = method_id,
      group = method_id
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.2, alpha = 0.8) +
    ggplot2::geom_text(
      data = unavailable_annotation,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      color = unname(method_colors["C2"]),
      size = 3.6
    ) +
    ggplot2::facet_wrap(~fit_stage, nrow = 1) +
    ggplot2::scale_color_manual(
      values = method_colors,
      breaks = expected_method_ids,
      labels = unname(expected_method_labels)
    ) +
    ggplot2::scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      limits = c(0, 1)
    ) +
    ggplot2::labs(
      x = "Prior standard-deviation grid",
      y = "Estimated prior weight"
    ) +
    base_plot_theme()
}

build_lfdr_comparison_plot_data <- function() {
  comparisons <- lfdr_pairwise_metrics[, c(
    "fit_stage", "reference_method_id", "comparison_method_id",
    "reference_method_label", "comparison_method_label"
  )]
  rows <- vector("list", nrow(comparisons))
  for (comparison_index in seq_len(nrow(comparisons))) {
    comparison <- comparisons[comparison_index, ]
    suffix <- if (comparison$fit_stage == "Raw") "raw" else "bf"
    reference_column <- paste0(
      comparison$reference_method_id,
      "_",
      suffix
    )
    comparison_column <- paste0(
      comparison$comparison_method_id,
      "_",
      suffix
    )
    rows[[comparison_index]] <- data.frame(
      fit_stage = comparison$fit_stage,
      comparison = paste(
        comparison$reference_method_label,
        "vs",
        comparison$comparison_method_label
      ),
      pair_key = lfdr_wide$pair_key,
      gene_id = lfdr_wide$gene_id,
      variant_id = lfdr_wide$variant_id,
      reference_lfdr = lfdr_wide[[reference_column]],
      comparison_lfdr = lfdr_wide[[comparison_column]],
      stringsAsFactors = FALSE
    )
  }
  plot_data <- do.call(rbind, rows)
  plot_data$panel <- paste0(plot_data$fit_stage, ": ", plot_data$comparison)
  panel_order <- vapply(
    rows,
    function(table) paste0(table$fit_stage[1], ": ", table$comparison[1]),
    character(1)
  )
  plot_data$panel <- factor(plot_data$panel, levels = panel_order)
  plot_data$difference <-
    plot_data$comparison_lfdr - plot_data$reference_lfdr
  plot_data
}

lfdr_comparison_plot_data <- build_lfdr_comparison_plot_data()

plot_lfdr_comparisons <- function() {
  ggplot2::ggplot(
    lfdr_comparison_plot_data,
    ggplot2::aes(x = reference_lfdr, y = comparison_lfdr)
  ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      linewidth = 0.55,
      color = "#555555"
    ) +
    ggplot2::geom_point(size = 0.55, alpha = 0.22, color = "#0072B2") +
    ggplot2::facet_wrap(~panel, ncol = 2) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      x = "Reference-method lfdr",
      y = "Comparison-method lfdr"
    ) +
    base_plot_theme() +
    ggplot2::theme(legend.position = "none")
}

plot_lfdr_difference_distributions <- function() {
  ggplot2::ggplot(
    lfdr_comparison_plot_data,
    ggplot2::aes(x = difference, fill = fit_stage)
  ) +
    ggplot2::geom_histogram(
      bins = 60,
      boundary = 0,
      color = "white",
      linewidth = 0.15
    ) +
    ggplot2::geom_vline(xintercept = 0, color = "#333333", linewidth = 0.5) +
    ggplot2::facet_wrap(~panel, ncol = 2, scales = "free_y") +
    ggplot2::scale_fill_manual(
      values = c(Raw = "#56B4E9", `BF-adjusted` = "#CC79A7")
    ) +
    ggplot2::labs(
      x = "Comparison lfdr minus reference lfdr",
      y = "Number of units"
    ) +
    base_plot_theme() +
    ggplot2::theme(legend.position = "none")
}

plot_discovery_counts <- function() {
  plot_data <- discovery_summary
  plot_data$fit_stage <- factor(plot_data$fit_stage, levels = stage_levels)
  plot_data$method_label <- factor(
    plot_data$method_label,
    levels = unname(expected_method_labels)
  )
  plot_data$available <- is.finite(plot_data$discovered_units)
  plot_data$plot_count <- ifelse(
    plot_data$available,
    plot_data$discovered_units,
    0
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = method_label, y = plot_count, fill = method_id)
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = ifelse(available, as.character(discovered_units), "Unavailable")
      ),
      vjust = -0.35,
      size = 3.6
    ) +
    ggplot2::facet_wrap(~fit_stage, nrow = 1) +
    ggplot2::scale_fill_manual(values = method_colors) +
    ggplot2::scale_y_continuous(
      limits = c(0, max(plot_data$plot_count) * 1.12),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(x = NULL, y = "Cumulative-lfdr calls at alpha = 0.05") +
    base_plot_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1),
      legend.position = "none"
    )
}

raw_diagonal_pi0 <- method_stage_summary$pi0[
  raw_rows & method_stage_summary$method_id == "diagonal"
]
raw_c1_pi0 <- method_stage_summary$pi0[
  raw_rows & method_stage_summary$method_id == "C1"
]
raw_c2_pi0 <- method_stage_summary$pi0[
  raw_rows & method_stage_summary$method_id == "C2"
]
bf_diagonal_pi0 <- method_stage_summary$pi0[
  bf_rows & method_stage_summary$method_id == "diagonal"
]
bf_c1_pi0 <- method_stage_summary$pi0[
  bf_rows & method_stage_summary$method_id == "C1"
]
raw_diagonal_calls <- method_stage_summary$discovered_units[
  raw_rows & method_stage_summary$method_id == "diagonal"
]
bf_diagonal_calls <- method_stage_summary$discovered_units[
  bf_rows & method_stage_summary$method_id == "diagonal"
]
c1_estimating_pairs <- correlation_diagnostics$n_estimating_pairs[
  correlation_diagnostics$matrix_id == "C1"
]
c2_estimating_pairs <- correlation_diagnostics$n_estimating_pairs[
  correlation_diagnostics$matrix_id == "C2"
]
