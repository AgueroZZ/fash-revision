# Load, validate, format, and plot the one-variant-per-gene refit caches.

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required one-variant-per-gene cache file is missing: ", path)
  }
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

format_decimal <- function(x, digits = 4) {
  formatC(x, format = "f", digits = digits)
}

format_integer <- function(x) {
  format(as.integer(x), big.mark = ",", scientific = FALSE)
}

format_range <- function(x, digits = 4, signed = FALSE) {
  x <- range(as.numeric(x))
  formatter <- if (signed) {
    function(value) formatC(value, format = "f", digits = digits, flag = "+")
  } else {
    function(value) formatC(value, format = "f", digits = digits)
  }
  paste0(formatter(x[1]), "\u2013", formatter(x[2]))
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
    '<div class="one-gene-table-scroll">\n',
    as.character(rendered),
    '\n</div>\n',
    sep = ""
  )
  invisible(rendered)
}

expected_seeds <- seq(12345L, 102345L, by = 10000L)
expected_fit_stages <- c("Raw", "BF-adjusted")
expected_full_units <- 1009173L
expected_genes <- 6362L
expected_alpha <- 0.05
expected_grid <- sort(c(0, exp(-0.5 * seq(0, 10, by = 0.2))))
expected_settings <- list(
  num_basis = 20,
  betaprec = 0,
  order = 1,
  pred_step = 1,
  likelihood = "gaussian",
  penalty = 10
)
expected_source_sizes <- c(
  raw_fit = 576110361,
  bf_adjusted_fit = 581122554
)
expected_source_sha256 <- c(
  raw_fit = "9e9aaf8a405f7ca83439990656666035fc0a74bb1ae2858ace79944fac6ec929",
  bf_adjusted_fit = "3e3d6735b9da734a3ab16ed53713908c02d0af4e15a659760fc5f892ef64023b"
)
expected_required_files <- c(
  "configuration.rds",
  "selection.csv",
  "prior_weight_comparison.csv",
  "lfdr_comparison.csv",
  "comparison_summary.csv",
  "thinned_fash_fit.rds"
)

cache_parent <- file.path(
  "output",
  "revision_simulations",
  "internal"
)
aggregate_path <- file.path(
  cache_parent,
  "one_variant_per_gene_refit_seed_summary.csv"
)
aggregate_summary <- read_required_csv(aggregate_path)

required_summary_columns <- c(
  "seed", "fit_stage", "n_full_units", "n_genes", "n_units",
  "full_pi0", "thinned_pi0", "pi0_difference",
  "prior_total_variation", "pearson_lfdr", "spearman_lfdr",
  "mean_absolute_lfdr_difference", "median_absolute_lfdr_difference",
  "rmse_lfdr", "full_mean_lfdr", "thinned_mean_lfdr",
  "full_fdr_calls", "thinned_fdr_calls", "fdr_call_intersection",
  "fdr_call_union", "fdr_call_jaccard", "alpha",
  "raw_eb_refit_elapsed_seconds", "bf_update_elapsed_seconds"
)
if (!all(required_summary_columns %in% names(aggregate_summary)) ||
    nrow(aggregate_summary) != length(expected_seeds) * 2L ||
    !identical(sort(unique(aggregate_summary$seed)), expected_seeds) ||
    !identical(sort(unique(aggregate_summary$fit_stage)),
               sort(expected_fit_stages)) ||
    any(table(aggregate_summary$seed) != 2L) ||
    any(table(aggregate_summary$fit_stage) != length(expected_seeds)) ||
    any(aggregate_summary$n_full_units != expected_full_units) ||
    any(aggregate_summary$n_genes != expected_genes) ||
    any(aggregate_summary$n_units != expected_genes) ||
    any(abs(aggregate_summary$alpha - expected_alpha) > 1e-12)) {
  stop("The aggregate one-variant-per-gene summary has unexpected structure.")
}
numeric_summary <- aggregate_summary[, vapply(
  aggregate_summary,
  is.numeric,
  logical(1)
), drop = FALSE]
if (any(!is.finite(as.matrix(numeric_summary))) ||
    any(aggregate_summary$full_pi0 < 0 | aggregate_summary$full_pi0 > 1) ||
    any(aggregate_summary$thinned_pi0 < 0 |
        aggregate_summary$thinned_pi0 > 1) ||
    any(aggregate_summary$prior_total_variation < 0 |
        aggregate_summary$prior_total_variation > 1) ||
    any(aggregate_summary$fdr_call_jaccard < 0 |
        aggregate_summary$fdr_call_jaccard > 1)) {
  stop("The aggregate one-variant-per-gene summary contains invalid values.")
}

full_discovery_path <- file.path(
  cache_parent,
  "one_variant_per_gene_full_bf_discovery_summary.csv"
)
full_bf_discovery_summary <- read_required_csv(full_discovery_path)
required_full_discovery_columns <- c(
  "fit", "n_tested_pairs", "n_tested_genes", "discovered_pairs",
  "discovered_genes", "alpha"
)
full_discovery_valid <-
  all(required_full_discovery_columns %in% names(full_bf_discovery_summary)) &&
  nrow(full_bf_discovery_summary) == 1L &&
  identical(
    full_bf_discovery_summary$fit,
    "Full-data BF-adjusted FASH(1)"
  ) &&
  identical(full_bf_discovery_summary$n_tested_pairs, expected_full_units) &&
  identical(full_bf_discovery_summary$n_tested_genes, expected_genes) &&
  identical(full_bf_discovery_summary$discovered_pairs, 9205L) &&
  identical(full_bf_discovery_summary$discovered_genes, 1177L) &&
  abs(full_bf_discovery_summary$alpha - expected_alpha) <= 1e-12
if (!isTRUE(full_discovery_valid)) {
  stop("The fixed full-data BF discovery summary failed validation.")
}

variant_count_path <- file.path(
  cache_parent,
  "one_variant_per_gene_variant_counts.csv"
)
variant_count_by_gene <- read_required_csv(variant_count_path)
required_variant_count_columns <- c(
  "gene_id", "n_tested_variants"
)
variant_count_values <- variant_count_by_gene$n_tested_variants
variant_count_valid <-
  all(required_variant_count_columns %in% names(variant_count_by_gene)) &&
  nrow(variant_count_by_gene) == expected_genes &&
  all(nzchar(variant_count_by_gene$gene_id)) &&
  anyDuplicated(variant_count_by_gene$gene_id) == 0L &&
  is.numeric(variant_count_values) &&
  all(is.finite(variant_count_values)) &&
  all(variant_count_values >= 1L) &&
  all(variant_count_values == as.integer(variant_count_values)) &&
  sum(variant_count_values) == expected_full_units
if (!isTRUE(variant_count_valid)) {
  stop("The full-data per-gene variant-count cache failed validation.")
}

cache_directories <- setNames(
  file.path(
    cache_parent,
    paste0("one_variant_per_gene_refit_seed", expected_seeds)
  ),
  as.character(expected_seeds)
)
configurations <- vector("list", length(expected_seeds))
selections <- vector("list", length(expected_seeds))
names(configurations) <- names(selections) <- as.character(expected_seeds)
selection_rows <- vector("list", length(expected_seeds))
validation_rows <- vector("list", length(expected_seeds))
fit_bundle_rows <- vector("list", length(expected_seeds))

for (seed_index in seq_along(expected_seeds)) {
  seed <- expected_seeds[seed_index]
  seed_id <- as.character(seed)
  cache_directory <- cache_directories[[seed_id]]
  required_paths <- file.path(cache_directory, expected_required_files)
  if (!dir.exists(cache_directory) || any(!file.exists(required_paths))) {
    stop("The seed ", seed, " one-variant-per-gene cache is incomplete.")
  }

  configuration <- readRDS(file.path(cache_directory, "configuration.rds"))
  selection <- read_required_csv(file.path(cache_directory, "selection.csv"))
  seed_summary <- read_required_csv(file.path(
    cache_directory,
    "comparison_summary.csv"
  ))
  expected_seed_summary <- aggregate_summary[
    aggregate_summary$seed == seed,
    ,
    drop = FALSE
  ]
  seed_summary <- seed_summary[
    order(match(seed_summary$fit_stage, expected_fit_stages)),
    ,
    drop = FALSE
  ]
  expected_seed_summary <- expected_seed_summary[
    order(match(expected_seed_summary$fit_stage, expected_fit_stages)),
    ,
    drop = FALSE
  ]
  rownames(seed_summary) <- rownames(expected_seed_summary) <- NULL
  if (!isTRUE(all.equal(seed_summary, expected_seed_summary, tolerance = 1e-12))) {
    stop("The seed ", seed, " summary does not match the aggregate cache.")
  }

  validation <- configuration$likelihood_validation
  source_files <- configuration$source_files
  configuration_valid <-
    identical(
      configuration$experiment,
      "One uniformly random tested variant per gene FASH(1) refit"
    ) &&
    identical(configuration$seed, seed) &&
    identical(configuration$n_full_units, expected_full_units) &&
    identical(configuration$n_genes, expected_genes) &&
    identical(configuration$n_selected_units, expected_genes) &&
    isTRUE(all.equal(configuration$alpha, expected_alpha)) &&
    isTRUE(all.equal(configuration$full_raw_settings, expected_settings)) &&
    isTRUE(all.equal(configuration$full_psd_grid, expected_grid)) &&
    length(configuration$raw_fit_warnings) == 0L &&
    length(configuration$bf_update_warnings) == 0L &&
    grepl("fash_eb_est", configuration$refit_strategy, fixed = TRUE) &&
    is.list(validation) &&
    validation$selection_seed %in% expected_seeds &&
    identical(validation$n_units, 24L) &&
    validation$likelihood_max_absolute_difference <= 1e-8 &&
    validation$raw_lfdr_max_absolute_difference <= 1e-6 &&
    validation$prior_total_variation <= 1e-6 &&
    length(validation$warnings) == 0L &&
    is.list(source_files) &&
    identical(
      unname(vapply(source_files, `[[`, numeric(1), "size_bytes")),
      unname(expected_source_sizes)
    )
  if (!isTRUE(configuration_valid)) {
    stop("The seed ", seed, " configuration failed validation.")
  }

  required_selection_columns <- c(
    "seed", "fash_index", "pair_key", "gene_id", "variant_id"
  )
  if (!all(required_selection_columns %in% names(selection)) ||
      nrow(selection) != expected_genes ||
      any(selection$seed != seed) ||
      any(!is.finite(selection$fash_index)) ||
      any(selection$fash_index < 1 | selection$fash_index > expected_full_units) ||
      anyDuplicated(selection$fash_index) ||
      any(!nzchar(selection$pair_key)) || anyDuplicated(selection$pair_key) ||
      any(!nzchar(selection$gene_id)) || anyDuplicated(selection$gene_id) ||
      any(!nzchar(selection$variant_id)) ||
      !identical(sub("_.*$", "", selection$pair_key), selection$gene_id) ||
      !identical(sub("^[^_]+_", "", selection$pair_key),
                 selection$variant_id)) {
    stop("The seed ", seed, " selection is not one unique pair per gene.")
  }

  variant_counts <- table(selection$variant_id)
  selection_rows[[seed_index]] <- data.frame(
    Seed = seed,
    `Gene-variant pairs` = nrow(selection),
    `Unique genes` = length(unique(selection$gene_id)),
    `Unique source variants` = length(variant_counts),
    `Repeated cross-gene assignments` = sum(variant_counts - 1L),
    `Maximum genes per variant` = max(variant_counts),
    check.names = FALSE
  )
  validation_rows[[seed_index]] <- data.frame(
    cache_seed = seed,
    validation_seed = validation$selection_seed,
    n_units = validation$n_units,
    likelihood_max_absolute_difference =
      validation$likelihood_max_absolute_difference,
    raw_lfdr_max_absolute_difference =
      validation$raw_lfdr_max_absolute_difference,
    prior_total_variation = validation$prior_total_variation,
    elapsed_seconds = validation$elapsed_seconds,
    stringsAsFactors = FALSE
  )
  fit_info <- file.info(file.path(cache_directory, "thinned_fash_fit.rds"))
  fit_bundle_rows[[seed_index]] <- data.frame(
    Seed = seed,
    `Fit bundle size (MB)` = unname(fit_info$size) / 1024^2,
    `Generated at (UTC)` = configuration$generated_at,
    check.names = FALSE
  )
  configurations[[seed_id]] <- configuration
  selections[[seed_id]] <- selection
}
if (!identical(
  sort(variant_count_by_gene$gene_id),
  sort(selections[[as.character(expected_seeds[1])]]$gene_id)
)) {
  stop("The variant-count cache and random-thinning caches use different genes.")
}

specific_seed <- expected_seeds[1]
specific_seed_id <- as.character(specific_seed)
specific_seed_lfdr <- read_required_csv(file.path(
  cache_directories[[specific_seed_id]],
  "lfdr_comparison.csv"
))
required_specific_lfdr_columns <- c(
  "fit_stage", "pair_key", "full_lfdr", "thinned_lfdr",
  "lfdr_difference", "absolute_difference", "full_fdr_call",
  "thinned_fdr_call"
)
specific_seed_lfdr_valid <-
  all(required_specific_lfdr_columns %in% names(specific_seed_lfdr)) &&
  nrow(specific_seed_lfdr) == expected_genes * length(expected_fit_stages) &&
  identical(
    sort(unique(specific_seed_lfdr$fit_stage)),
    sort(expected_fit_stages)
  ) &&
  all(table(specific_seed_lfdr$fit_stage) == expected_genes) &&
  all(vapply(expected_fit_stages, function(fit_stage) {
    stage_rows <- specific_seed_lfdr$fit_stage == fit_stage
    identical(
      specific_seed_lfdr$pair_key[stage_rows],
      selections[[specific_seed_id]]$pair_key
    )
  }, logical(1))) &&
  all(is.finite(specific_seed_lfdr$full_lfdr)) &&
  all(is.finite(specific_seed_lfdr$thinned_lfdr)) &&
  all(specific_seed_lfdr$full_lfdr >= 0 &
        specific_seed_lfdr$full_lfdr <= 1) &&
  all(specific_seed_lfdr$thinned_lfdr >= 0 &
        specific_seed_lfdr$thinned_lfdr <= 1)
if (!isTRUE(specific_seed_lfdr_valid)) {
  stop("The seed-specific lfdr comparison cache failed validation.")
}

selection_audit <- do.call(rbind, selection_rows)
rownames(selection_audit) <- NULL
validation_records <- do.call(rbind, validation_rows)
validation_records <- validation_records[
  order(validation_records$validation_seed, validation_records$cache_seed),
  ,
  drop = FALSE
]
for (validation_seed in unique(validation_records$validation_seed)) {
  rows <- validation_records[
    validation_records$validation_seed == validation_seed,
    c(
      "n_units", "likelihood_max_absolute_difference",
      "raw_lfdr_max_absolute_difference", "prior_total_variation",
      "elapsed_seconds"
    ),
    drop = FALSE
  ]
  if (nrow(unique(rows)) != 1L) {
    stop("Repeated validation records disagree for seed ", validation_seed, ".")
  }
}
validation_summary <- validation_records[
  !duplicated(validation_records$validation_seed),
  ,
  drop = FALSE
]
if (!identical(
  validation_summary$validation_seed,
  c(12345L, 22345L, 62345L)
)) {
  stop("The expected three direct-likelihood validations are not present.")
}
fit_bundle_audit <- do.call(rbind, fit_bundle_rows)
rownames(fit_bundle_audit) <- NULL

source_reference <- configurations[["12345"]]$source_files
for (seed_id in names(configurations)) {
  if (!isTRUE(all.equal(
    configurations[[seed_id]]$source_files,
    source_reference
  ))) {
    stop("Source-file provenance differs across seed caches.")
  }
}
source_provenance_table <- data.frame(
  Artifact = c(
    "Full raw FASH(1), penalty = 10",
    "Full BF-adjusted FASH(1)"
  ),
  `File name` = basename(c(
    source_reference$raw_fit$path,
    source_reference$bf_adjusted_fit$path
  )),
  `Size (MB)` = format_decimal(c(
    source_reference$raw_fit$size_bytes,
    source_reference$bf_adjusted_fit$size_bytes
  ) / 1024^2, 1),
  `Recorded mtime (UTC)` = c(
    source_reference$raw_fit$modification_time,
    source_reference$bf_adjusted_fit$modification_time
  ),
  `SHA-256 verified` = unname(expected_source_sha256),
  check.names = FALSE
)

design_table <- data.frame(
  Quantity = c(
    "Full gene-variant pairs",
    "Unique genes",
    "Pairs retained per thinning",
    "Random thinning seeds",
    "FASH order",
    "PSD grid points",
    "Raw EB penalty",
    "BF update",
    "Comparison alpha"
  ),
  Value = c(
    format_integer(expected_full_units),
    format_integer(expected_genes),
    format_integer(expected_genes),
    paste(expected_seeds, collapse = ", "),
    as.character(expected_settings$order),
    as.character(length(expected_grid)),
    "10 (fixed to the original analysis)",
    "Applied after every raw refit",
    format_decimal(expected_alpha, 2)
  ),
  stringsAsFactors = FALSE
)

variant_count_quartiles <- as.numeric(stats::quantile(
  variant_count_values,
  probs = c(0.25, 0.5, 0.75),
  names = FALSE
))
variant_count_min <- min(variant_count_values)
variant_count_q1 <- variant_count_quartiles[1]
variant_count_median <- variant_count_quartiles[2]
variant_count_mean <- mean(variant_count_values)
variant_count_q3 <- variant_count_quartiles[3]
variant_count_max <- max(variant_count_values)
variant_count_summary_table <- data.frame(
  Statistic = c("Minimum", "Q1", "Median", "Mean", "Q3", "Maximum"),
  `Tested variants per gene` = c(
    format_integer(variant_count_min),
    format_integer(variant_count_q1),
    format_integer(variant_count_median),
    format_decimal(variant_count_mean, 1),
    format_integer(variant_count_q3),
    format_integer(variant_count_max)
  ),
  check.names = FALSE
)

validation_table <- data.frame(
  `Validation selection seed` = validation_summary$validation_seed,
  Units = validation_summary$n_units,
  `Maximum |likelihood difference|` = format(
    validation_summary$likelihood_max_absolute_difference,
    scientific = TRUE,
    digits = 3
  ),
  `Maximum |raw lfdr difference|` = format(
    validation_summary$raw_lfdr_max_absolute_difference,
    scientific = TRUE,
    digits = 3
  ),
  `Prior total-variation distance` = format(
    validation_summary$prior_total_variation,
    scientific = TRUE,
    digits = 3
  ),
  `Elapsed seconds` = format_decimal(validation_summary$elapsed_seconds, 2),
  check.names = FALSE
)

aggregate_summary$fit_stage <- factor(
  aggregate_summary$fit_stage,
  levels = expected_fit_stages
)
aggregate_summary <- aggregate_summary[
  order(aggregate_summary$seed, aggregate_summary$fit_stage),
  ,
  drop = FALSE
]
raw_metrics <- aggregate_summary[aggregate_summary$fit_stage == "Raw", ]
bf_metrics <- aggregate_summary[
  aggregate_summary$fit_stage == "BF-adjusted",
]

make_stage_range_row <- function(rows, interpretation) {
  data.frame(
    Fit = as.character(rows$fit_stage[1]),
    `Full pi0` = format_decimal(unique(rows$full_pi0), 4),
    `Thinned pi0 range` = format_range(rows$thinned_pi0, 4),
    `Delta pi0 range` = format_range(rows$pi0_difference, 4, signed = TRUE),
    `Prior TV range` = format_range(rows$prior_total_variation, 4),
    `lfdr Spearman range` = format_range(rows$spearman_lfdr, 4),
    `lfdr MAE range` = format_range(
      rows$mean_absolute_lfdr_difference,
      4
    ),
    `Discovery Jaccard range` = format_range(rows$fdr_call_jaccard, 4),
    Interpretation = interpretation,
    check.names = FALSE
  )
}

fit_range_table <- rbind(
  make_stage_range_row(
    raw_metrics,
    "Ranks are stable, but absolute calibration changes."
  ),
  make_stage_range_row(
    bf_metrics,
    "Null weight, lfdr values, and calls are practically stable."
  )
)
rownames(fit_range_table) <- NULL

seed_result_table <- data.frame(
  Seed = aggregate_summary$seed,
  Fit = as.character(aggregate_summary$fit_stage),
  `Full pi0` = format_decimal(aggregate_summary$full_pi0, 4),
  `Thinned pi0` = format_decimal(aggregate_summary$thinned_pi0, 4),
  `Delta pi0` = formatC(
    aggregate_summary$pi0_difference,
    format = "f",
    digits = 4,
    flag = "+"
  ),
  `lfdr Spearman` = format_decimal(aggregate_summary$spearman_lfdr, 4),
  `lfdr MAE` = format_decimal(
    aggregate_summary$mean_absolute_lfdr_difference,
    4
  ),
  `Call Jaccard` = format_decimal(aggregate_summary$fdr_call_jaccard, 4),
  check.names = FALSE
)

specific_seed_metrics <- aggregate_summary[
  aggregate_summary$seed == specific_seed,
  ,
  drop = FALSE
]
specific_seed_metrics <- specific_seed_metrics[
  order(specific_seed_metrics$fit_stage),
  ,
  drop = FALSE
]
specific_seed_table <- data.frame(
  Fit = as.character(specific_seed_metrics$fit_stage),
  `Full pi0` = format_decimal(specific_seed_metrics$full_pi0, 4),
  `Thinned pi0` = format_decimal(specific_seed_metrics$thinned_pi0, 4),
  `lfdr Spearman` = format_decimal(specific_seed_metrics$spearman_lfdr, 4),
  `lfdr MAE` = format_decimal(
    specific_seed_metrics$mean_absolute_lfdr_difference,
    4
  ),
  `Reference calls on sampled pairs` = specific_seed_metrics$full_fdr_calls,
  `Thinned-refit calls` = specific_seed_metrics$thinned_fdr_calls,
  `Discovery Jaccard` = format_decimal(
    specific_seed_metrics$fdr_call_jaccard,
    4
  ),
  check.names = FALSE
)

thinned_gene_discovery_table <- data.frame(
  Seed = bf_metrics$seed,
  `Genes discovered by the thinned BF-adjusted refit` =
    bf_metrics$thinned_fdr_calls,
  check.names = FALSE
)

full_data_discovery_table <- data.frame(
  Fit = full_bf_discovery_summary$fit,
  `Tested gene-variant pairs` = format_integer(
    full_bf_discovery_summary$n_tested_pairs
  ),
  `Tested genes` = format_integer(full_bf_discovery_summary$n_tested_genes),
  `Discovered gene-variant pairs` = format_integer(
    full_bf_discovery_summary$discovered_pairs
  ),
  `Genes with at least one discovered dynamic eQTL` = format_integer(
    full_bf_discovery_summary$discovered_genes
  ),
  check.names = FALSE
)

full_raw_pi0 <- unique(raw_metrics$full_pi0)
full_bf_pi0 <- unique(bf_metrics$full_pi0)
raw_pi0_range <- range(raw_metrics$thinned_pi0)
bf_pi0_range <- range(bf_metrics$thinned_pi0)
raw_lfdr_mae_range <- range(raw_metrics$mean_absolute_lfdr_difference)
bf_lfdr_mae_range <- range(bf_metrics$mean_absolute_lfdr_difference)
raw_spearman_range <- range(raw_metrics$spearman_lfdr)
bf_spearman_range <- range(bf_metrics$spearman_lfdr)
raw_jaccard_range <- range(raw_metrics$fdr_call_jaccard)
bf_jaccard_range <- range(bf_metrics$fdr_call_jaccard)
bf_pi0_difference_range <- range(bf_metrics$pi0_difference)
unique_variant_range <- range(selection_audit$`Unique source variants`)
repeated_assignment_range <- range(
  selection_audit$`Repeated cross-gene assignments`
)
full_bf_pair_discoveries <- full_bf_discovery_summary$discovered_pairs
full_bf_gene_discoveries <- full_bf_discovery_summary$discovered_genes
bf_thinned_trained_gene_range <- range(bf_metrics$thinned_fdr_calls)
bf_thinned_trained_gene_range_width <- diff(bf_thinned_trained_gene_range)
bf_thinned_gene_median <- stats::median(bf_metrics$thinned_fdr_calls)
bf_thinned_gene_iqr <- as.numeric(stats::quantile(
  bf_metrics$thinned_fdr_calls,
  probs = c(0.25, 0.75),
  names = FALSE
))
n_seeds <- length(expected_seeds)

stability_colors <- c("Raw" = "#D55E00", "BF-adjusted" = "#0072B2")
plot_data <- aggregate_summary

stability_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.title = ggplot2::element_blank(),
      plot.title.position = "plot"
    )
}

plot_variants_per_gene_histogram <- function() {
  ggplot2::ggplot(
    variant_count_by_gene,
    ggplot2::aes(x = n_tested_variants)
  ) +
    ggplot2::geom_histogram(
      bins = 40,
      color = "white",
      fill = "#0072B2",
      linewidth = 0.25,
      alpha = 0.88
    ) +
    ggplot2::geom_vline(
      xintercept = variant_count_median,
      color = "#D55E00",
      linewidth = 0.8,
      linetype = 2
    ) +
    ggplot2::scale_x_log10(
      breaks = c(2, 5, 10, 25, 50, 100, 250, 500, 1000, 2000),
      labels = scales::label_comma()
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      title = "Tested variants per gene in the full FASH analysis",
      subtitle = paste0(
        "Median = ", format_integer(variant_count_median),
        "; IQR = ", format_integer(variant_count_q1),
        "-", format_integer(variant_count_q3),
        "; dashed line marks the median"
      ),
      x = "Tested variants per gene (log scale)",
      y = "Number of genes"
    ) +
    stability_theme() +
    ggplot2::theme(legend.position = "none")
}

plot_specific_seed_lfdr <- function() {
  specific_plot_data <- specific_seed_lfdr
  specific_plot_data$fit_stage <- factor(
    specific_plot_data$fit_stage,
    levels = expected_fit_stages
  )
  annotation_data <- data.frame(
    fit_stage = factor(
      as.character(specific_seed_metrics$fit_stage),
      levels = expected_fit_stages
    ),
    x = 0.03,
    y = 0.97,
    label = paste0(
      "Spearman = ",
      format_decimal(specific_seed_metrics$spearman_lfdr, 4),
      "\nMAE = ",
      format_decimal(
        specific_seed_metrics$mean_absolute_lfdr_difference,
        4
      )
    ),
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(
    specific_plot_data,
    ggplot2::aes(
      x = full_lfdr,
      y = thinned_lfdr,
      color = fit_stage
    )
  ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      color = "#5B5B5B",
      linewidth = 0.65,
      linetype = 2
    ) +
    ggplot2::geom_point(size = 0.7, alpha = 0.24) +
    ggplot2::geom_label(
      data = annotation_data,
      ggplot2::aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.2,
      linewidth = 0.2,
      fill = "white"
    ) +
    ggplot2::facet_wrap(~fit_stage, nrow = 1) +
    ggplot2::scale_color_manual(values = stability_colors) +
    ggplot2::coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = paste0(
        "Seed ",
        specific_seed,
        ": full-fit reference versus thinned-refit lfdr"
      ),
      subtitle = "Each point is one selected gene-variant pair; the dashed line is equality",
      x = "lfdr from the full-data fit",
      y = "lfdr from the thinned refit"
    ) +
    stability_theme() +
    ggplot2::theme(
      legend.position = "none",
      panel.spacing.x = grid::unit(1.5, "lines")
    )
}

plot_pi0_stability <- function() {
  reference <- unique(plot_data[, c("fit_stage", "full_pi0")])
  plot_data$seed_set <- "Ten thinned refits"
  reference$seed_set <- "Ten thinned refits"
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = seed_set,
      y = thinned_pi0,
      color = fit_stage,
      fill = fit_stage
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.48,
      alpha = 0.16,
      outlier.shape = NA,
      linewidth = 0.65
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(
        width = 0.08,
        height = 0,
        seed = 12345
      ),
      size = 2.2,
      alpha = 0.85
    ) +
    ggplot2::geom_point(
      data = reference,
      ggplot2::aes(
        x = seed_set,
        y = full_pi0,
        color = fit_stage
      ),
      inherit.aes = FALSE,
      shape = 23,
      fill = "white",
      stroke = 1,
      size = 3.8
    ) +
    ggplot2::facet_wrap(~fit_stage, scales = "free_y") +
    ggplot2::scale_color_manual(values = stability_colors) +
    ggplot2::scale_fill_manual(values = stability_colors) +
    ggplot2::labs(
      title = "Null-weight variability across ten thinnings",
      subtitle = "Each point is one seed; diamonds are the corresponding full-data estimates",
      x = NULL,
      y = "Estimated pi0"
    ) +
    stability_theme() +
    ggplot2::theme(legend.position = "none")
}

plot_lfdr_mae <- function() {
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = fit_stage,
      y = mean_absolute_lfdr_difference,
      color = fit_stage,
      fill = fit_stage
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.52,
      alpha = 0.16,
      outlier.shape = NA,
      linewidth = 0.65
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(
        width = 0.08,
        height = 0,
        seed = 12345
      ),
      size = 2.2,
      alpha = 0.85
    ) +
    ggplot2::scale_color_manual(values = stability_colors) +
    ggplot2::scale_fill_manual(values = stability_colors) +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      title = "Paired lfdr disagreement across ten thinnings",
      subtitle = "Each point is one seed; the vertical axis is logarithmic",
      x = NULL,
      y = "Mean absolute lfdr difference"
    ) +
    stability_theme() +
    ggplot2::theme(legend.position = "none")
}

plot_discovery_jaccard <- function() {
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = fit_stage,
      y = fdr_call_jaccard,
      color = fit_stage,
      fill = fit_stage
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.52,
      alpha = 0.16,
      outlier.shape = NA,
      linewidth = 0.65
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(
        width = 0.08,
        height = 0,
        seed = 12345
      ),
      size = 2.2,
      alpha = 0.85
    ) +
    ggplot2::scale_color_manual(values = stability_colors) +
    ggplot2::scale_fill_manual(values = stability_colors) +
    ggplot2::coord_cartesian(ylim = c(0.75, 1)) +
    ggplot2::labs(
      title = "Conditional discovery overlap across ten thinnings",
      subtitle = "Each point is one seed; both fits are evaluated on the same sampled pairs",
      x = NULL,
      y = "Discovery-set Jaccard"
    ) +
    stability_theme() +
    ggplot2::theme(legend.position = "none")
}

plot_thinned_gene_discoveries <- function() {
  gene_plot_data <- data.frame(
    seed = bf_metrics$seed,
    fit = "Thinned BF-adjusted refit",
    discovered_genes = bf_metrics$thinned_fdr_calls,
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(
    gene_plot_data,
    ggplot2::aes(
      x = fit,
      y = discovered_genes
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.42,
      color = stability_colors[["BF-adjusted"]],
      fill = stability_colors[["BF-adjusted"]],
      alpha = 0.16,
      outlier.shape = NA,
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(
        width = 0.07,
        height = 0,
        seed = 12345
      ),
      color = stability_colors[["BF-adjusted"]],
      size = 2.4,
      alpha = 0.9
    ) +
    ggplot2::scale_y_continuous(breaks = scales::breaks_pretty(n = 6)) +
    ggplot2::labs(
      title = "Genes discovered by the thinned refit across ten seeds",
      subtitle = "Each point is one independent one-variant-per-gene selection",
      x = NULL,
      y = "Discovered genes at cumulative-lfdr FDR 0.05"
    ) +
    stability_theme() +
    ggplot2::theme(legend.position = "none")
}
