# Load, validate, format, and plot the retained R5 balanced-thinning caches.

source(
  "code/revision_simulations/r5_balanced_thinning/balanced_thinning_helpers.R"
)

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required R5 cache file is missing: ", path)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

read_required_rds <- function(path) {
  if (!file.exists(path)) {
    stop("Required R5 cache file is missing: ", path)
  }
  readRDS(path)
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
                                    caption,
                                    align = NULL,
                                    minimum_width = "760px",
                                    digits = NULL) {
  arguments <- list(
    x = table,
    format = "html",
    caption = caption,
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
    '<div class="r5-table-scroll">\n',
    as.character(rendered),
    '\n</div>\n',
    sep = ""
  )
  invisible(rendered)
}

expected_experiment <-
  "Revision Real-data Sensitivity R5: balanced variant thinning"
expected_seed_experiment <- paste(
  "R5 BF-updated FASH(1) refit after uniformly sampling exactly",
  "20 tested variants per eligible gene"
)
expected_seeds <- seq(12345L, 102345L, by = 10000L)
expected_target <- 20L
expected_full_units <- 1009173L
expected_full_genes <- 6362L
expected_eligible_genes <- 6352L
expected_excluded_genes <- 10L
expected_selected_units <- 127040L
expected_alpha <- 0.05
expected_penalty <- 10
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

cache_directory <- file.path(
  "output",
  "revision_simulations",
  "r5_balanced_thinning"
)
configuration <- read_required_rds(file.path(
  cache_directory,
  "configuration.rds"
))

source_sizes <- unname(vapply(
  configuration$source_files,
  `[[`,
  numeric(1),
  "size_bytes"
))
source_hashes <- unname(vapply(
  configuration$source_files,
  `[[`,
  character(1),
  "sha256"
))
configuration_valid <-
  identical(configuration$experiment, expected_experiment) &&
  identical(configuration$seeds, expected_seeds) &&
  identical(configuration$target_per_gene, expected_target) &&
  isTRUE(all.equal(configuration$alpha, expected_alpha)) &&
  identical(configuration$fit_stage, "BF-updated") &&
  identical(configuration$n_full_units, expected_full_units) &&
  identical(configuration$n_full_genes, expected_full_genes) &&
  identical(configuration$n_eligible_genes, expected_eligible_genes) &&
  identical(configuration$n_excluded_genes, expected_excluded_genes) &&
  identical(configuration$n_selected_units_per_seed, expected_selected_units) &&
  isTRUE(all.equal(configuration$full_raw_settings, expected_settings)) &&
  isTRUE(all.equal(configuration$full_psd_grid, expected_grid)) &&
  identical(source_sizes, unname(expected_source_sizes)) &&
  identical(source_hashes, unname(expected_source_sha256))
if (!isTRUE(configuration_valid)) {
  stop("The shared R5 configuration failed validation.")
}

variant_count_by_gene <- read_required_csv(file.path(
  cache_directory,
  "full_variant_counts.csv"
))
variant_count_values <- variant_count_by_gene$n_tested_variants
variant_count_valid <-
  identical(names(variant_count_by_gene), c("gene_id", "n_tested_variants")) &&
  nrow(variant_count_by_gene) == expected_full_genes &&
  all(nzchar(variant_count_by_gene$gene_id)) &&
  anyDuplicated(variant_count_by_gene$gene_id) == 0L &&
  all(is.finite(variant_count_values)) &&
  all(variant_count_values == as.integer(variant_count_values)) &&
  all(variant_count_values >= 1L) &&
  sum(variant_count_values) == expected_full_units &&
  sum(variant_count_values >= expected_target) == expected_eligible_genes
if (!isTRUE(variant_count_valid)) {
  stop("The full-data per-gene variant counts failed validation.")
}

excluded_genes <- read_required_csv(file.path(
  cache_directory,
  "excluded_genes.csv"
))
excluded_valid <-
  identical(
    names(excluded_genes),
    c("gene_id", "n_tested_variants", "full_data_discovered_gene")
  ) &&
  nrow(excluded_genes) == expected_excluded_genes &&
  all(excluded_genes$n_tested_variants < expected_target) &&
  !any(excluded_genes$full_data_discovered_gene) &&
  identical(
    sort(excluded_genes$gene_id),
    sort(variant_count_by_gene$gene_id[
      variant_count_by_gene$n_tested_variants < expected_target
    ])
  )
if (!isTRUE(excluded_valid)) {
  stop("The R5 excluded-gene cache failed validation.")
}

full_discovery_summary <- read_required_csv(file.path(
  cache_directory,
  "full_discovery_summary.csv"
))
full_discovered_genes <- read_required_csv(file.path(
  cache_directory,
  "full_discovered_genes.csv"
))
full_summary_valid <-
  nrow(full_discovery_summary) == 1L &&
  identical(
    full_discovery_summary$fit,
    "Full-data BF-adjusted FASH(1)"
  ) &&
  identical(full_discovery_summary$n_tested_pairs, expected_full_units) &&
  identical(full_discovery_summary$n_tested_genes, expected_full_genes) &&
  identical(full_discovery_summary$discovered_pairs, 9205L) &&
  identical(full_discovery_summary$discovered_genes, 1177L) &&
  identical(full_discovery_summary$eligible_genes, expected_eligible_genes) &&
  identical(full_discovery_summary$eligible_full_discovered_genes, 1177L) &&
  identical(full_discovery_summary$excluded_genes, expected_excluded_genes) &&
  identical(full_discovery_summary$excluded_full_discovered_genes, 0L) &&
  isTRUE(all.equal(full_discovery_summary$alpha, expected_alpha)) &&
  identical(names(full_discovered_genes), "gene_id") &&
  nrow(full_discovered_genes) == 1177L &&
  all(nzchar(full_discovered_genes$gene_id)) &&
  anyDuplicated(full_discovered_genes$gene_id) == 0L &&
  all(full_discovered_genes$gene_id %in% variant_count_by_gene$gene_id) &&
  !any(full_discovered_genes$gene_id %in% excluded_genes$gene_id)
if (!isTRUE(full_summary_valid)) {
  stop("The fixed full-data BF discovery caches failed validation.")
}

seed_summary <- read_required_csv(file.path(
  cache_directory,
  "seed_summary.csv"
))
required_seed_columns <- c(
  "seed", "fit_stage", "target_per_gene", "n_full_units",
  "n_full_genes", "n_eligible_genes", "n_selected_units", "full_pi0",
  "thinned_pi0", "pi0_difference", "prior_total_variation", "n_units",
  "pearson_lfdr", "spearman_lfdr", "mean_absolute_lfdr_difference",
  "median_absolute_lfdr_difference", "rmse_lfdr", "full_fdr_calls",
  "thinned_fdr_calls", "fdr_call_intersection", "fdr_call_union",
  "fdr_call_jaccard", "alpha", "thinned_discovered_genes",
  "thinned_discovered_full_genes",
  "thinned_discovered_additional_genes", "elapsed_seconds", "warning_count"
)
numeric_seed_summary <- seed_summary[, vapply(
  seed_summary,
  is.numeric,
  logical(1)
), drop = FALSE]
seed_summary_valid <-
  all(required_seed_columns %in% names(seed_summary)) &&
  nrow(seed_summary) == length(expected_seeds) &&
  identical(seed_summary$seed, expected_seeds) &&
  all(seed_summary$fit_stage == "BF-updated") &&
  all(seed_summary$target_per_gene == expected_target) &&
  all(seed_summary$n_full_units == expected_full_units) &&
  all(seed_summary$n_full_genes == expected_full_genes) &&
  all(seed_summary$n_eligible_genes == expected_eligible_genes) &&
  all(seed_summary$n_selected_units == expected_selected_units) &&
  all(seed_summary$n_units == expected_selected_units) &&
  all(abs(seed_summary$alpha - expected_alpha) <= 1e-12) &&
  all(seed_summary$warning_count == 0L) &&
  all(is.finite(as.matrix(numeric_seed_summary))) &&
  all(seed_summary$full_pi0 >= 0 & seed_summary$full_pi0 <= 1) &&
  all(seed_summary$thinned_pi0 >= 0 & seed_summary$thinned_pi0 <= 1) &&
  all(seed_summary$spearman_lfdr >= -1 & seed_summary$spearman_lfdr <= 1) &&
  all(seed_summary$fdr_call_jaccard >= 0 & seed_summary$fdr_call_jaccard <= 1)
if (!isTRUE(seed_summary_valid)) {
  stop("The aggregate R5 seed summary failed validation.")
}

eligible_gene_ids <- sort(variant_count_by_gene$gene_id[
  variant_count_by_gene$n_tested_variants >= expected_target
])
seed_directories <- setNames(
  file.path(cache_directory, paste0("seed", expected_seeds)),
  as.character(expected_seeds)
)
discovered_gene_sets <- vector("list", length(expected_seeds))
names(discovered_gene_sets) <- as.character(expected_seeds)

for (seed_index in seq_along(expected_seeds)) {
  seed <- expected_seeds[seed_index]
  seed_id <- as.character(seed)
  seed_directory <- seed_directories[[seed_id]]
  required_files <- c(
    "configuration.rds", "selection.rds", "lfdr_comparison.rds",
    "prior_weight_comparison.csv", "comparison_summary.csv",
    "discovered_genes.csv"
  )
  required_paths <- file.path(seed_directory, required_files)
  if (!dir.exists(seed_directory) || any(!file.exists(required_paths))) {
    stop("The retained R5 cache is incomplete for seed ", seed, ".")
  }

  seed_configuration <- read_required_rds(file.path(
    seed_directory,
    "configuration.rds"
  ))
  seed_source_sizes <- unname(vapply(
    seed_configuration$source_files,
    `[[`,
    numeric(1),
    "size_bytes"
  ))
  seed_source_hashes <- unname(vapply(
    seed_configuration$source_files,
    `[[`,
    character(1),
    "sha256"
  ))
  seed_configuration_valid <-
    identical(seed_configuration$experiment, expected_seed_experiment) &&
    identical(seed_configuration$seed, seed) &&
    identical(seed_configuration$target_per_gene, expected_target) &&
    identical(seed_configuration$n_full_units, expected_full_units) &&
    identical(seed_configuration$n_full_genes, expected_full_genes) &&
    identical(seed_configuration$n_eligible_genes, expected_eligible_genes) &&
    identical(seed_configuration$n_excluded_genes, expected_excluded_genes) &&
    identical(seed_configuration$n_selected_units, expected_selected_units) &&
    isTRUE(all.equal(seed_configuration$alpha, expected_alpha)) &&
    identical(seed_configuration$fit_stage, "BF-updated") &&
    isTRUE(all.equal(seed_configuration$raw_penalty, expected_penalty)) &&
    isTRUE(all.equal(seed_configuration$full_raw_settings, expected_settings)) &&
    isTRUE(all.equal(seed_configuration$full_psd_grid, expected_grid)) &&
    length(seed_configuration$warnings) == 0L &&
    identical(seed_source_sizes, unname(expected_source_sizes)) &&
    identical(seed_source_hashes, unname(expected_source_sha256)) &&
    grepl("BF_update", seed_configuration$refit_strategy, fixed = TRUE)
  if (!isTRUE(seed_configuration_valid)) {
    stop("The R5 configuration failed validation for seed ", seed, ".")
  }

  selection <- read_required_rds(file.path(seed_directory, "selection.rds"))
  selection_valid <-
    identical(
      names(selection),
      c(
        "seed", "target_per_gene", "fash_index", "pair_key", "gene_id",
        "variant_id"
      )
    ) &&
    nrow(selection) == expected_selected_units &&
    all(selection$seed == seed) &&
    all(selection$target_per_gene == expected_target) &&
    all(is.finite(selection$fash_index)) &&
    all(selection$fash_index >= 1L & selection$fash_index <= expected_full_units) &&
    anyDuplicated(selection$fash_index) == 0L &&
    anyDuplicated(selection$pair_key) == 0L &&
    all(nzchar(selection$pair_key)) &&
    all(nzchar(selection$gene_id)) &&
    all(nzchar(selection$variant_id)) &&
    identical(sub("_.*$", "", selection$pair_key), selection$gene_id) &&
    identical(sub("^[^_]+_", "", selection$pair_key), selection$variant_id) &&
    identical(sort(unique(selection$gene_id)), eligible_gene_ids) &&
    all(table(selection$gene_id) == expected_target)
  if (!isTRUE(selection_valid)) {
    stop("The balanced selection failed validation for seed ", seed, ".")
  }

  seed_comparison <- read_required_csv(file.path(
    seed_directory,
    "comparison_summary.csv"
  ))
  expected_comparison <- seed_summary[seed_summary$seed == seed, , drop = FALSE]
  rownames(seed_comparison) <- rownames(expected_comparison) <- NULL
  if (!isTRUE(all.equal(
    seed_comparison,
    expected_comparison,
    tolerance = 1e-12,
    check.attributes = FALSE
  ))) {
    stop("The seed summary disagrees with the aggregate cache for seed ", seed, ".")
  }

  prior_comparison <- read_required_csv(file.path(
    seed_directory,
    "prior_weight_comparison.csv"
  ))
  prior_valid <-
    all(c(
      "psd", "full_weight", "thinned_weight", "difference",
      "absolute_difference"
    ) %in% names(prior_comparison)) &&
    nrow(prior_comparison) >= 2L &&
    all(is.finite(as.matrix(prior_comparison))) &&
    sum(prior_comparison$psd == 0) == 1L &&
    abs(sum(prior_comparison$full_weight) - 1) <= 1e-6 &&
    abs(sum(prior_comparison$thinned_weight) - 1) <= 1e-6 &&
    abs(prior_comparison$full_weight[prior_comparison$psd == 0] -
          expected_comparison$full_pi0) <= 1e-12 &&
    abs(prior_comparison$thinned_weight[prior_comparison$psd == 0] -
          expected_comparison$thinned_pi0) <= 1e-12
  if (!isTRUE(prior_valid)) {
    stop("The prior-weight comparison failed validation for seed ", seed, ".")
  }

  discovered_genes <- read_required_csv(file.path(
    seed_directory,
    "discovered_genes.csv"
  ))
  discovered_valid <-
    identical(
      names(discovered_genes),
      c("gene_id", "discovered_pairs", "full_data_discovered_gene")
    ) &&
    nrow(discovered_genes) == expected_comparison$thinned_discovered_genes &&
    all(nzchar(discovered_genes$gene_id)) &&
    anyDuplicated(discovered_genes$gene_id) == 0L &&
    all(discovered_genes$gene_id %in% eligible_gene_ids) &&
    all(discovered_genes$discovered_pairs >= 1L) &&
    all(discovered_genes$discovered_pairs <= expected_target) &&
    sum(discovered_genes$discovered_pairs) == expected_comparison$thinned_fdr_calls &&
    sum(discovered_genes$full_data_discovered_gene) ==
      expected_comparison$thinned_discovered_full_genes &&
    sum(!discovered_genes$full_data_discovered_gene) ==
      expected_comparison$thinned_discovered_additional_genes &&
    identical(
      discovered_genes$full_data_discovered_gene,
      discovered_genes$gene_id %in% full_discovered_genes$gene_id
    )
  if (!isTRUE(discovered_valid)) {
    stop("The discovered-gene cache failed validation for seed ", seed, ".")
  }
  discovered_gene_sets[[seed_id]] <- discovered_genes$gene_id
}

specific_seed <- expected_seeds[1]
specific_seed_id <- as.character(specific_seed)
specific_seed_lfdr <- read_required_rds(file.path(
  seed_directories[[specific_seed_id]],
  "lfdr_comparison.rds"
))
specific_selection <- read_required_rds(file.path(
  seed_directories[[specific_seed_id]],
  "selection.rds"
))
required_lfdr_columns <- c(
  "pair_key", "gene_id", "variant_id", "full_lfdr", "thinned_lfdr",
  "lfdr_difference", "absolute_difference", "full_fdr_call",
  "thinned_fdr_call"
)
specific_lfdr_valid <-
  identical(names(specific_seed_lfdr), required_lfdr_columns) &&
  nrow(specific_seed_lfdr) == expected_selected_units &&
  identical(specific_seed_lfdr$pair_key, specific_selection$pair_key) &&
  identical(specific_seed_lfdr$gene_id, specific_selection$gene_id) &&
  identical(specific_seed_lfdr$variant_id, specific_selection$variant_id) &&
  all(is.finite(specific_seed_lfdr$full_lfdr)) &&
  all(is.finite(specific_seed_lfdr$thinned_lfdr)) &&
  all(specific_seed_lfdr$full_lfdr >= 0 & specific_seed_lfdr$full_lfdr <= 1) &&
  all(specific_seed_lfdr$thinned_lfdr >= 0 &
        specific_seed_lfdr$thinned_lfdr <= 1) &&
  max(abs(
    specific_seed_lfdr$lfdr_difference -
      (specific_seed_lfdr$thinned_lfdr - specific_seed_lfdr$full_lfdr)
  )) <= 1e-12 &&
  max(abs(
    specific_seed_lfdr$absolute_difference -
      abs(specific_seed_lfdr$lfdr_difference)
  )) <= 1e-12
if (!isTRUE(specific_lfdr_valid)) {
  stop("The seed 12345 paired-lfdr cache failed validation.")
}

cumulative_discoveries <- cumulative_gene_union(
  discovered_gene_sets = discovered_gene_sets,
  seeds = expected_seeds,
  full_discovered_genes = full_discovered_genes$gene_id
)
discovery_envelope <- seed_subset_union_envelope(
  discovered_gene_sets = discovered_gene_sets,
  seeds = expected_seeds,
  full_discovered_genes = full_discovered_genes$gene_id
)
if (any(diff(cumulative_discoveries$cumulative_unique_genes) < 0L) ||
    any(diff(cumulative_discoveries$cumulative_full_genes_recovered) < 0L) ||
    cumulative_discoveries$cumulative_unique_genes[10] != 1113L ||
    cumulative_discoveries$cumulative_full_genes_recovered[10] != 1075L ||
    any(discovery_envelope$unique_min > discovery_envelope$unique_median) ||
    any(discovery_envelope$unique_median > discovery_envelope$unique_max)) {
  stop("The cumulative discovered-gene summaries failed validation.")
}

variant_count_quantiles <- as.numeric(stats::quantile(
  variant_count_values,
  probs = c(0.25, 0.5, 0.75),
  names = FALSE
))
variant_count_min <- min(variant_count_values)
variant_count_q1 <- variant_count_quantiles[1]
variant_count_median <- variant_count_quantiles[2]
variant_count_mean <- mean(variant_count_values)
variant_count_q3 <- variant_count_quantiles[3]
variant_count_max <- max(variant_count_values)

design_table <- data.frame(
  Quantity = c(
    "Full gene-variant pairs",
    "Genes in the full fit",
    "Variants retained per eligible gene",
    "Eligible genes",
    "Genes excluded because they have fewer than 20 variants",
    "Gene-variant pairs retained per seed",
    "Prespecified random seeds",
    "Reported fit stage",
    "Raw EB penalty before BF update",
    "Cumulative-lfdr FDR threshold"
  ),
  Value = c(
    format_integer(expected_full_units),
    format_integer(expected_full_genes),
    format_integer(expected_target),
    format_integer(expected_eligible_genes),
    paste0(format_integer(expected_excluded_genes),
           " (none discovered in the full analysis)"),
    format_integer(expected_selected_units),
    paste(expected_seeds, collapse = ", "),
    "BF-updated only",
    format_integer(expected_penalty),
    format_decimal(expected_alpha, 2)
  ),
  stringsAsFactors = FALSE
)

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

specific_seed_metrics <- seed_summary[seed_summary$seed == specific_seed, ]
specific_seed_table <- data.frame(
  Seed = specific_seed,
  `Retained pairs` = format_integer(specific_seed_metrics$n_selected_units),
  `Eligible genes` = format_integer(specific_seed_metrics$n_eligible_genes),
  `Full pi0` = format_decimal(specific_seed_metrics$full_pi0, 4),
  `Thinned pi0` = format_decimal(specific_seed_metrics$thinned_pi0, 4),
  `lfdr Spearman` = format_decimal(specific_seed_metrics$spearman_lfdr, 4),
  `lfdr MAE` = format_decimal(
    specific_seed_metrics$mean_absolute_lfdr_difference,
    4
  ),
  `Full-reference calls on retained pairs` =
    format_integer(specific_seed_metrics$full_fdr_calls),
  `Thinned-refit calls` = format_integer(specific_seed_metrics$thinned_fdr_calls),
  `Call Jaccard` = format_decimal(specific_seed_metrics$fdr_call_jaccard, 4),
  check.names = FALSE
)

across_seed_summary_table <- data.frame(
  Quantity = c(
    "Estimated pi0",
    "Change in pi0 from the full fit",
    "Paired lfdr Spearman correlation",
    "Paired lfdr mean absolute difference",
    "Genes discovered per thinned refit",
    "Full-data discovered genes recovered per thinned refit",
    "Additional genes discovered per thinned refit"
  ),
  `Full-data reference` = c(
    format_decimal(unique(seed_summary$full_pi0), 4),
    "0",
    "1",
    "0",
    format_integer(full_discovery_summary$discovered_genes),
    format_integer(full_discovery_summary$discovered_genes),
    "0"
  ),
  `Range across ten thinned refits` = c(
    format_range(seed_summary$thinned_pi0, 4),
    format_range(seed_summary$pi0_difference, 4, signed = TRUE),
    format_range(seed_summary$spearman_lfdr, 4),
    format_range(seed_summary$mean_absolute_lfdr_difference, 4),
    paste0(
      format_integer(min(seed_summary$thinned_discovered_genes)),
      "\u2013",
      format_integer(max(seed_summary$thinned_discovered_genes))
    ),
    paste0(
      format_integer(min(seed_summary$thinned_discovered_full_genes)),
      "\u2013",
      format_integer(max(seed_summary$thinned_discovered_full_genes))
    ),
    paste0(
      format_integer(min(seed_summary$thinned_discovered_additional_genes)),
      "\u2013",
      format_integer(max(seed_summary$thinned_discovered_additional_genes))
    )
  ),
  check.names = FALSE
)

cumulative_discovery_table <- data.frame(
  `Seeds considered` = cumulative_discoveries$n_seeds,
  `Added seed` = cumulative_discoveries$added_seed,
  `Genes discovered in the added seed` =
    cumulative_discoveries$seed_discovered_genes,
  `Cumulative unique genes discovered` =
    cumulative_discoveries$cumulative_unique_genes,
  `Cumulative full-data discovered genes recovered` =
    cumulative_discoveries$cumulative_full_genes_recovered,
  check.names = FALSE
)

full_pi0 <- unique(seed_summary$full_pi0)
thinned_pi0_range <- range(seed_summary$thinned_pi0)
pi0_difference_range <- range(seed_summary$pi0_difference)
lfdr_spearman_range <- range(seed_summary$spearman_lfdr)
lfdr_mae_range <- range(seed_summary$mean_absolute_lfdr_difference)
per_seed_discovery_range <- range(seed_summary$thinned_discovered_genes)
full_discovered_gene_count <- full_discovery_summary$discovered_genes
final_cumulative_unique <- tail(cumulative_discoveries$cumulative_unique_genes, 1)
final_cumulative_recovered <- tail(
  cumulative_discoveries$cumulative_full_genes_recovered,
  1
)
final_additional_unique <- final_cumulative_unique - final_cumulative_recovered
final_recovery_fraction <- final_cumulative_recovered / full_discovered_gene_count
n_seeds <- length(expected_seeds)

r5_colors <- c(
  "Thinned refits" = "#0072B2",
  "Full-data reference" = "#D55E00",
  "All thinned discoveries" = "#0072B2",
  "Full-data discoveries recovered" = "#E69F00"
)

r5_theme <- function() {
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
      bins = 45,
      color = "white",
      fill = r5_colors[["Thinned refits"]],
      linewidth = 0.25,
      alpha = 0.88
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(
        xintercept = expected_target,
        linetype = "Retention target: 20 variants"
      ),
      color = r5_colors[["Full-data reference"]],
      linewidth = 0.85
    ) +
    ggplot2::scale_linetype_manual(values = c("Retention target: 20 variants" = 2)) +
    ggplot2::scale_x_log10(
      breaks = c(2, 5, 10, 20, 50, 100, 250, 500, 1000, 2000),
      labels = scales::label_comma()
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      title = "Tested variants contributed by each gene",
      subtitle = paste0(
        format_integer(expected_eligible_genes),
        " of ",
        format_integer(expected_full_genes),
        " genes have at least 20 tested variants"
      ),
      x = "Tested variants per gene (log scale)",
      y = "Number of genes",
      linetype = NULL
    ) +
    r5_theme()
}

plot_seed_12345_lfdr <- function() {
  annotation <- paste0(
    "Spearman = ",
    format_decimal(specific_seed_metrics$spearman_lfdr, 4),
    "\nMAE = ",
    format_decimal(specific_seed_metrics$mean_absolute_lfdr_difference, 4),
    "\nFull pi0 = ",
    format_decimal(specific_seed_metrics$full_pi0, 4),
    "\nThinned pi0 = ",
    format_decimal(specific_seed_metrics$thinned_pi0, 4)
  )
  ggplot2::ggplot(
    specific_seed_lfdr,
    ggplot2::aes(x = full_lfdr, y = thinned_lfdr)
  ) +
    ggplot2::geom_abline(
      ggplot2::aes(
        intercept = 0,
        slope = 1,
        linetype = "Equality"
      ),
      color = "#5B5B5B",
      linewidth = 0.75
    ) +
    ggplot2::geom_point(
      color = r5_colors[["Thinned refits"]],
      size = 0.65,
      alpha = 0.23
    ) +
    ggplot2::annotate(
      "label",
      x = 0.035,
      y = 0.965,
      label = annotation,
      hjust = 0,
      vjust = 1,
      size = 3.25,
      linewidth = 0.2,
      fill = "white"
    ) +
    ggplot2::scale_linetype_manual(values = c("Equality" = 2)) +
    ggplot2::coord_equal(
      xlim = c(0, 1),
      ylim = c(0, 1),
      expand = FALSE
    ) +
    ggplot2::labs(
      title = "Seed 12345: paired lfdr values",
      subtitle = paste0(
        "The same ",
        format_integer(expected_selected_units),
        " retained gene-variant pairs appear on both axes"
      ),
      x = "lfdr from the full-data BF-updated fit",
      y = "lfdr after balanced thinning and BF update",
      linetype = NULL
    ) +
    r5_theme()
}

plot_null_weight_stability <- function() {
  plot_data <- data.frame(
    Analysis = "20 variants per eligible gene",
    pi0 = seed_summary$thinned_pi0,
    Series = "Thinned refits",
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Analysis, y = pi0, color = Series, fill = Series)
  ) +
    ggplot2::geom_boxplot(
      width = 0.42,
      alpha = 0.16,
      outlier.shape = NA,
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(
        width = 0.065,
        height = 0,
        seed = 12345
      ),
      size = 2.4,
      alpha = 0.9
    ) +
    ggplot2::geom_hline(
      ggplot2::aes(
        yintercept = full_pi0,
        linetype = "Full-data reference"
      ),
      color = r5_colors[["Full-data reference"]],
      linewidth = 0.85
    ) +
    ggplot2::scale_color_manual(values = r5_colors["Thinned refits"]) +
    ggplot2::scale_fill_manual(values = r5_colors["Thinned refits"]) +
    ggplot2::scale_linetype_manual(values = c("Full-data reference" = 2)) +
    ggplot2::scale_y_continuous(
      limits = c(0.93, 0.94),
      breaks = seq(0.93, 0.94, by = 0.002),
      labels = scales::label_number(accuracy = 0.001),
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::labs(
      title = "Estimated null weight across ten balanced thinnings",
      subtitle = "Each point is one prespecified seed; the dashed line is the fixed full-data estimate",
      x = NULL,
      y = "Estimated pi0",
      color = NULL,
      fill = NULL,
      linetype = NULL
    ) +
    r5_theme()
}

plot_paired_lfdr_agreement <- function() {
  agreement_data <- rbind(
    data.frame(
      Seed = seed_summary$seed,
      Metric = "Mean absolute difference",
      Value = seed_summary$mean_absolute_lfdr_difference,
      stringsAsFactors = FALSE
    ),
    data.frame(
      Seed = seed_summary$seed,
      Metric = "Spearman correlation",
      Value = seed_summary$spearman_lfdr,
      stringsAsFactors = FALSE
    )
  )
  agreement_data$Metric <- factor(
    agreement_data$Metric,
    levels = c("Mean absolute difference", "Spearman correlation")
  )
  ggplot2::ggplot(
    agreement_data,
    ggplot2::aes(x = Metric, y = Value, color = Metric, fill = Metric)
  ) +
    ggplot2::geom_boxplot(
      width = 0.44,
      alpha = 0.16,
      outlier.shape = NA,
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(
        width = 0.065,
        height = 0,
        seed = 12345
      ),
      size = 2.3,
      alpha = 0.9
    ) +
    ggplot2::facet_wrap(~Metric, scales = "free_y") +
    ggplot2::scale_color_manual(values = c(
      "Mean absolute difference" = "#0072B2",
      "Spearman correlation" = "#009E73"
    )) +
    ggplot2::scale_fill_manual(values = c(
      "Mean absolute difference" = "#0072B2",
      "Spearman correlation" = "#009E73"
    )) +
    ggplot2::labs(
      title = "Paired lfdr agreement across ten balanced thinnings",
      subtitle = "Each point summarizes all retained pairs from one seed",
      x = NULL,
      y = NULL
    ) +
    r5_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "none"
    )
}

plot_cumulative_discovered_genes <- function() {
  curve_data <- rbind(
    data.frame(
      n_seeds = cumulative_discoveries$n_seeds,
      Count = cumulative_discoveries$cumulative_unique_genes,
      Series = "All thinned discoveries",
      stringsAsFactors = FALSE
    ),
    data.frame(
      n_seeds = cumulative_discoveries$n_seeds,
      Count = cumulative_discoveries$cumulative_full_genes_recovered,
      Series = "Full-data discoveries recovered",
      stringsAsFactors = FALSE
    )
  )
  ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      data = discovery_envelope,
      ggplot2::aes(
        x = n_seeds,
        ymin = unique_min,
        ymax = unique_max,
        fill = "Range over all seed subsets"
      ),
      alpha = 0.15
    ) +
    ggplot2::geom_line(
      data = curve_data,
      ggplot2::aes(x = n_seeds, y = Count, color = Series),
      linewidth = 1
    ) +
    ggplot2::geom_point(
      data = curve_data,
      ggplot2::aes(x = n_seeds, y = Count, color = Series),
      size = 2.2
    ) +
    ggplot2::geom_hline(
      ggplot2::aes(
        yintercept = full_discovered_gene_count,
        linetype = "Full-data discovered genes"
      ),
      color = "#5B5B5B",
      linewidth = 0.85
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(expected_seeds),
      minor_breaks = NULL
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::label_comma(),
      breaks = scales::breaks_pretty(n = 7),
      expand = ggplot2::expansion(mult = c(0.03, 0.06))
    ) +
    ggplot2::scale_color_manual(values = r5_colors[c(
      "All thinned discoveries",
      "Full-data discoveries recovered"
    )]) +
    ggplot2::scale_fill_manual(values = c(
      "Range over all seed subsets" = "#0072B2"
    )) +
    ggplot2::scale_linetype_manual(values = c(
      "Full-data discovered genes" = 2
    )) +
    ggplot2::labs(
      title = "Cumulative unique genes discovered as seeds are added",
      subtitle = paste0(
        "Solid curves follow the prespecified seed order; the shaded band summarizes all seed subsets of each size"
      ),
      x = "Number of seeds considered",
      y = "Cumulative unique genes",
      color = NULL,
      fill = NULL,
      linetype = NULL
    ) +
    r5_theme()
}
