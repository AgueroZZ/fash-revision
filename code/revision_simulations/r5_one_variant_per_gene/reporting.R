# Load, validate, summarize, and plot the R5 one-variant-per-gene cache.

find_workflowr_root <- function() {
  if (file.exists("analysis/index.Rmd")) {
    return(normalizePath("."))
  }
  if (file.exists("coderepo-local/analysis/index.Rmd")) {
    return(normalizePath("coderepo-local"))
  }
  stop("Could not find the workflowr repository root.")
}

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

sha256_file <- function(path) {
  output <- system2(
    "shasum",
    args = c("-a", "256", normalizePath(path, mustWork = TRUE)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || length(output) != 1L) {
    stop("Could not compute SHA-256 for ", path, ".")
  }
  hash <- strsplit(output, "[[:space:]]+")[[1]][1]
  if (!grepl("^[0-9a-f]{64}$", hash)) {
    stop("Unexpected SHA-256 output for ", path, ".")
  }
  hash
}

require_reporting_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The ", package, " package is required for R5 reporting.")
  }
}

format_integer <- function(x) {
  format(as.integer(round(x)), big.mark = ",", scientific = FALSE, trim = TRUE)
}

format_decimal <- function(x, digits = 4L) {
  formatC(as.numeric(x), format = "f", digits = digits)
}

format_range <- function(x, digits = 4L) {
  paste0(
    format_decimal(min(x), digits),
    "--",
    format_decimal(max(x), digits)
  )
}

metric_five_number <- function(x) {
  values <- stats::quantile(
    as.numeric(x),
    probs = c(0, 0.25, 0.5, 0.75, 1),
    names = FALSE,
    type = 7
  )
  names(values) <- c("Minimum", "Q1", "Median", "Q3", "Maximum")
  values
}

format_metric_five_number <- function(x, digits) {
  values <- metric_five_number(x)
  stats::setNames(
    vapply(values, format_decimal, character(1), digits = digits),
    names(values)
  )
}

render_scrollable_table <- function(data,
                                    caption = NULL,
                                    align = NULL,
                                    minimum_width = "900px") {
  if (length(minimum_width) != 1L ||
      !grepl("^[0-9]+(px|em|rem|%)$", minimum_width)) {
    stop("minimum_width must be one CSS length using px, em, rem, or %.")
  }
  table_html <- knitr::kable(
    data,
    format = "html",
    escape = TRUE,
    caption = caption,
    align = align,
    row.names = FALSE
  )
  table_html <- kableExtra::kable_styling(
    table_html,
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE,
    position = "left"
  )
  knitr::asis_output(paste0(
    '<div class="r5-table-scroll">',
    '<div style="min-width:', minimum_width, ';">',
    as.character(table_html),
    "</div></div>"
  ))
}

require_reporting_package("ggplot2")
require_reporting_package("knitr")
require_reporting_package("kableExtra")
require_reporting_package("scales")

workflowr_root <- find_workflowr_root()
cache_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "r5_one_variant_per_gene_100seed_fashr0143"
)
expected_cache_id <- "r5_one_variant_per_gene_100seed_fashr0143"
expected_seeds <- seq(12345L, by = 10000L, length.out = 100L)
expected_full_units <- 1009173L
expected_full_genes <- 6362L
expected_selected_units <- expected_full_genes
expected_full_discovered_pairs <- 9214L
expected_full_discovered_genes <- 1176L
expected_full_pi0 <- 0.938159265061590
expected_alpha <- 0.05
expected_penalty <- 10
expected_grid <- sort(c(0, exp(-0.5 * seq(0, 10, by = 0.2))))
expected_source_sizes <- c(
  raw_fit = 576110361,
  bf_adjusted_fit = 581123504
)
expected_source_hashes <- c(
  raw_fit = "9e9aaf8a405f7ca83439990656666035fc0a74bb1ae2858ace79944fac6ec929",
  bf_adjusted_fit = "7f0ca9ab0fbeab89a13c83d2a0fb7c24195f7b5a5835f209399cf0e359001f50"
)
expected_fashr_provenance <- list(
  package = "fashr",
  version = "0.1.43",
  remote_sha = "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
)

shared_configuration <- read_required_rds(file.path(
  cache_directory,
  "configuration.rds"
))
variant_count_by_gene <- read_required_csv(file.path(
  cache_directory,
  "full_variant_counts.csv"
))
full_discovery_summary <- read_required_csv(file.path(
  cache_directory,
  "full_discovery_summary.csv"
))
full_discovered_genes <- read_required_csv(file.path(
  cache_directory,
  "full_discovered_genes.csv"
))
seed_summary <- read_required_csv(file.path(
  cache_directory,
  "seed_summary.csv"
))
gene_discovery_frequency <- read_required_csv(file.path(
  cache_directory,
  "gene_discovery_frequency.csv"
))

required_seed_columns <- c(
  "seed", "fit_stage", "target_per_gene", "n_full_units",
  "n_full_genes", "n_selected_units", "full_pi0", "thinned_pi0",
  "pi0_difference", "prior_total_variation", "n_units", "pearson_lfdr",
  "spearman_lfdr", "mean_absolute_lfdr_difference",
  "median_absolute_lfdr_difference", "rmse_lfdr", "full_fdr_calls",
  "thinned_fdr_calls", "fdr_call_intersection", "fdr_call_union",
  "fdr_call_jaccard", "alpha", "thinned_discovered_genes",
  "thinned_discovered_full_genes", "thinned_discovered_additional_genes",
  "elapsed_seconds", "warning_count"
)
stored_source_sizes <- vapply(
  shared_configuration$source_files,
  `[[`,
  numeric(1),
  "size_bytes"
)
stored_source_hashes <- vapply(
  shared_configuration$source_files,
  `[[`,
  character(1),
  "sha256"
)
stored_source_paths <- vapply(
  shared_configuration$source_files,
  `[[`,
  character(1),
  "path"
)
if (!identical(shared_configuration$cache_id, expected_cache_id) ||
    !identical(shared_configuration$seeds, expected_seeds) ||
    !identical(shared_configuration$target_per_gene, 1L) ||
    !identical(shared_configuration$fit_stage, "BF-updated") ||
    !isTRUE(all.equal(shared_configuration$alpha, expected_alpha)) ||
    !identical(shared_configuration$n_full_units, expected_full_units) ||
    !identical(shared_configuration$n_full_genes, expected_full_genes) ||
    !identical(
      shared_configuration$n_selected_units_per_seed,
      expected_selected_units
    ) ||
    !isTRUE(all.equal(shared_configuration$full_psd_grid, expected_grid)) ||
    !identical(stored_source_sizes, expected_source_sizes) ||
    !identical(stored_source_hashes, expected_source_hashes) ||
    !identical(
      shared_configuration$package_provenance,
      expected_fashr_provenance
    ) ||
    !identical(shared_configuration$full_raw_settings$penalty, expected_penalty)) {
  stop("The shared R5 configuration failed the retained-cache contract.")
}
if (any(!file.exists(stored_source_paths)) ||
    !identical(
      vapply(stored_source_paths, sha256_file, character(1)),
      expected_source_hashes
    )) {
  stop("At least one retained R5 source file changed after cache creation.")
}
if (!identical(names(variant_count_by_gene), c("gene_id", "n_tested_variants")) ||
    nrow(variant_count_by_gene) != expected_full_genes ||
    anyDuplicated(variant_count_by_gene$gene_id) ||
    any(variant_count_by_gene$n_tested_variants < 1L) ||
    sum(variant_count_by_gene$n_tested_variants) != expected_full_units) {
  stop("The full per-gene variant-count cache failed validation.")
}
if (nrow(full_discovery_summary) != 1L ||
    full_discovery_summary$n_tested_pairs != expected_full_units ||
    full_discovery_summary$n_tested_genes != expected_full_genes ||
    full_discovery_summary$discovered_pairs != expected_full_discovered_pairs ||
    full_discovery_summary$discovered_genes != expected_full_discovered_genes ||
    full_discovery_summary$representative_genes_per_seed != expected_full_genes ||
    !isTRUE(all.equal(full_discovery_summary$alpha, expected_alpha)) ||
    nrow(full_discovered_genes) != expected_full_discovered_genes ||
    anyDuplicated(full_discovered_genes$gene_id)) {
  stop("The full-data R5 discovery cache failed validation.")
}
if (!all(required_seed_columns %in% names(seed_summary)) ||
    nrow(seed_summary) != length(expected_seeds) ||
    !identical(seed_summary$seed, expected_seeds) ||
    any(seed_summary$fit_stage != "BF-updated") ||
    any(seed_summary$target_per_gene != 1L) ||
    any(seed_summary$n_full_units != expected_full_units) ||
    any(seed_summary$n_full_genes != expected_full_genes) ||
    any(seed_summary$n_selected_units != expected_selected_units) ||
    any(seed_summary$n_units != expected_selected_units) ||
    any(seed_summary$warning_count != 0L) ||
    any(!is.finite(as.matrix(seed_summary[, setdiff(
      required_seed_columns,
      "fit_stage"
    )]))) ||
    any(seed_summary$alpha != expected_alpha)) {
  stop("The 100-seed R5 summary failed validation.")
}
if (!identical(
  names(gene_discovery_frequency),
  c(
    "gene_id", "n_seeds_discovered", "discovery_frequency",
    "full_data_discovered_gene"
  )
) ||
    nrow(gene_discovery_frequency) != expected_full_genes ||
    anyDuplicated(gene_discovery_frequency$gene_id) ||
    any(gene_discovery_frequency$n_seeds_discovered < 0L) ||
    any(gene_discovery_frequency$n_seeds_discovered > length(expected_seeds)) ||
    any(gene_discovery_frequency$discovery_frequency < 0) ||
    any(gene_discovery_frequency$discovery_frequency > 1) ||
    sum(gene_discovery_frequency$full_data_discovered_gene) !=
      expected_full_discovered_genes) {
  stop("The gene-discovery frequency cache failed validation.")
}

seed_directories <- file.path(
  cache_directory,
  paste0("seed", expected_seeds)
)
seed_configurations <- lapply(seed_directories, function(directory) {
  read_required_rds(file.path(directory, "configuration.rds"))
})
seed_selections <- lapply(seed_directories, function(directory) {
  read_required_rds(file.path(directory, "selection.rds"))
})
for (seed_index in seq_along(expected_seeds)) {
  seed <- expected_seeds[seed_index]
  configuration <- seed_configurations[[seed_index]]
  selection <- seed_selections[[seed_index]]
  required_seed_files <- file.path(
    seed_directories[seed_index],
    c(
      "lfdr_comparison.rds", "prior_weight_comparison.csv",
      "discovered_genes.csv", "comparison_summary.csv"
    )
  )
  if (!all(file.exists(required_seed_files)) ||
      !identical(configuration$cache_id, expected_cache_id) ||
      !identical(configuration$seed, seed) ||
      !identical(configuration$target_per_gene, 1L) ||
      !identical(configuration$fit_stage, "BF-updated") ||
      !identical(
        configuration$package_provenance,
        expected_fashr_provenance
      ) ||
      length(configuration$warnings) != 0L ||
      nrow(selection) != expected_selected_units ||
      any(selection$seed != seed) ||
      any(selection$target_per_gene != 1L) ||
      anyDuplicated(selection$gene_id) ||
      anyDuplicated(selection$pair_key) ||
      length(unique(selection$gene_id)) != expected_full_genes) {
    stop("Seed ", seed, " failed the complete-cache validation.")
  }
}

repeated_variant_assignments <- vapply(seed_selections, function(selection) {
  nrow(selection) - length(unique(selection$variant_id))
}, integer(1))
maximum_variant_multiplicity <- max(vapply(seed_selections, function(selection) {
  max(table(selection$variant_id))
}, integer(1)))
repeated_variant_assignment_range <- range(repeated_variant_assignments)

seed_12345_lfdr <- read_required_rds(file.path(
  seed_directories[1],
  "lfdr_comparison.rds"
))
required_lfdr_columns <- c(
  "pair_key", "gene_id", "variant_id", "full_lfdr", "thinned_lfdr",
  "lfdr_difference", "absolute_difference", "full_fdr_call",
  "thinned_fdr_call"
)
if (!identical(names(seed_12345_lfdr), required_lfdr_columns) ||
    nrow(seed_12345_lfdr) != expected_selected_units ||
    anyDuplicated(seed_12345_lfdr$pair_key) ||
    any(!is.finite(seed_12345_lfdr$full_lfdr)) ||
    any(!is.finite(seed_12345_lfdr$thinned_lfdr))) {
  stop("The seed-12345 paired-lfdr cache failed validation.")
}

full_pi0 <- unique(seed_summary$full_pi0)
if (length(full_pi0) != 1L ||
    !isTRUE(all.equal(full_pi0, expected_full_pi0))) {
  stop("The full-data BF-updated null weight failed validation.")
}
pi0_range <- range(seed_summary$thinned_pi0)
pi0_difference_range <- range(seed_summary$pi0_difference)
lfdr_spearman_range <- range(seed_summary$spearman_lfdr)
lfdr_mae_range <- range(seed_summary$mean_absolute_lfdr_difference)
conditional_jaccard_range <- range(seed_summary$fdr_call_jaccard)
thinned_discovery_range <- range(seed_summary$thinned_discovered_genes)
ever_discovered <- sum(gene_discovery_frequency$n_seeds_discovered > 0L)
ever_discovered_full <- sum(
  gene_discovery_frequency$n_seeds_discovered > 0L &
    gene_discovery_frequency$full_data_discovered_gene
)
ever_discovered_additional <- sum(
  gene_discovery_frequency$n_seeds_discovered > 0L &
    !gene_discovery_frequency$full_data_discovered_gene
)

variant_count_values <- variant_count_by_gene$n_tested_variants
variant_count_quantiles <- stats::quantile(
  variant_count_values,
  probs = c(0, 0.25, 0.5, 0.75, 1),
  names = FALSE
)
variant_count_summary_table <- data.frame(
  Statistic = c("Minimum", "Q1", "Median", "Mean", "Q3", "Maximum"),
  `Tested variants per gene` = c(
    variant_count_quantiles[1],
    variant_count_quantiles[2],
    variant_count_quantiles[3],
    mean(variant_count_values),
    variant_count_quantiles[4],
    variant_count_quantiles[5]
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
variant_count_summary_table$`Tested variants per gene` <- vapply(
  variant_count_summary_table$`Tested variants per gene`,
  function(value) {
    if (abs(value - round(value)) < 1e-10) {
      format_integer(value)
    } else {
      format_decimal(value, 1)
    }
  },
  character(1)
)

design_table <- data.frame(
  Item = c(
    "Full-data empirical-Bayes units",
    "Genes represented in every refit",
    "Representative variants per gene",
    "Prespecified thinning seeds",
    "Selection rule",
    "Reported fit",
    "Raw empirical-Bayes penalty",
    "Discovery threshold",
    "Interpretation of seed variation"
  ),
  Value = c(
    paste0(format_integer(expected_full_units), " gene-variant pairs"),
    format_integer(expected_full_genes),
    "Exactly one",
    paste0(
      length(expected_seeds),
      " (",
      expected_seeds[1],
      " to ",
      expected_seeds[length(expected_seeds)],
      ")"
    ),
    paste(
      "Uniform within gene and independent of beta, SE, likelihood, lfdr,",
      "and calls"
    ),
    "BF-updated FASH(1) only",
    format_integer(expected_penalty),
    "Cumulative-lfdr FDR 0.05 within each specified universe",
    paste(
      "Monte Carlo sensitivity to representative-variant choice; not a",
      "cell-line bootstrap or confidence interval"
    )
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

seed_12345_summary <- seed_summary[seed_summary$seed == 12345L, , drop = FALSE]
seed_12345_table <- data.frame(
  Metric = c(
    "Full-data BF-updated null weight",
    "Thinned BF-updated null weight",
    "Paired lfdr mean absolute difference",
    "Paired lfdr Spearman correlation",
    "Reference calls on the sampled pairs",
    "Thinned-refit calls on the sampled pairs",
    "Conditional discovery Jaccard"
  ),
  Value = c(
    format_decimal(seed_12345_summary$full_pi0, 4),
    format_decimal(seed_12345_summary$thinned_pi0, 4),
    format_decimal(seed_12345_summary$mean_absolute_lfdr_difference, 4),
    format_decimal(seed_12345_summary$spearman_lfdr, 4),
    format_integer(seed_12345_summary$full_fdr_calls),
    format_integer(seed_12345_summary$thinned_fdr_calls),
    format_decimal(seed_12345_summary$fdr_call_jaccard, 4)
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

across_seed_specs <- list(
  list("Thinned BF-updated null weight", seed_summary$thinned_pi0, 4L),
  list("Absolute null-weight change", abs(seed_summary$pi0_difference), 4L),
  list(
    "Paired lfdr mean absolute difference",
    seed_summary$mean_absolute_lfdr_difference,
    4L
  ),
  list("Paired lfdr Spearman correlation", seed_summary$spearman_lfdr, 4L),
  list("Conditional discovery Jaccard", seed_summary$fdr_call_jaccard, 4L),
  list("Reference calls on sampled pairs", seed_summary$full_fdr_calls, 1L),
  list("Thinned-refit calls / genes", seed_summary$thinned_fdr_calls, 1L)
)
across_seed_summary_table <- do.call(rbind, lapply(across_seed_specs, function(spec) {
  values <- format_metric_five_number(spec[[2]], spec[[3]])
  data.frame(
    Metric = spec[[1]],
    Minimum = values[["Minimum"]],
    Q1 = values[["Q1"]],
    Median = values[["Median"]],
    Q3 = values[["Q3"]],
    Maximum = values[["Maximum"]],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))
rownames(across_seed_summary_table) <- NULL

full_status <- gene_discovery_frequency$full_data_discovered_gene
discovery_frequency_summary_table <- data.frame(
  `Gene set` = c(
    "All tested genes",
    "Full-data discovered genes",
    "Genes not discovered in the full data"
  ),
  Genes = c(
    nrow(gene_discovery_frequency),
    sum(full_status),
    sum(!full_status)
  ),
  `Discovered in at least one thinning seed` = c(
    ever_discovered,
    ever_discovered_full,
    ever_discovered_additional
  ),
  `Mean seed discovery frequency` = c(
    mean(gene_discovery_frequency$discovery_frequency),
    mean(gene_discovery_frequency$discovery_frequency[full_status]),
    mean(gene_discovery_frequency$discovery_frequency[!full_status])
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
discovery_frequency_summary_table$Genes <- format_integer(
  discovery_frequency_summary_table$Genes
)
discovery_frequency_summary_table$`Discovered in at least one thinning seed` <-
  format_integer(
    discovery_frequency_summary_table$`Discovered in at least one thinning seed`
  )
discovery_frequency_summary_table$`Mean seed discovery frequency` <-
  scales::label_percent(accuracy = 0.1)(
    discovery_frequency_summary_table$`Mean seed discovery frequency`
  )

r5_blue <- "#0072B2"
r5_orange <- "#E69F00"
r5_green <- "#009E73"
r5_theme <- ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    plot.subtitle = ggplot2::element_text(color = "grey30"),
    legend.position = "bottom",
    panel.grid.minor = ggplot2::element_blank()
  )

plot_variants_per_gene_histogram <- function() {
  ggplot2::ggplot(
    variant_count_by_gene,
    ggplot2::aes(x = n_tested_variants)
  ) +
    ggplot2::geom_histogram(
      bins = 45,
      fill = r5_blue,
      color = "white",
      linewidth = 0.25
    ) +
    ggplot2::scale_x_log10(
      breaks = c(2, 5, 10, 20, 50, 100, 250, 500, 1000, 2000),
      labels = scales::label_comma()
    ) +
    ggplot2::labs(
      title = "Unequal variant representation in the full fit",
      subtitle = paste(
        "The sensitivity refit gives each of the 6,362 genes one",
        "outcome-independent representative"
      ),
      x = "Tested variants per gene (log scale)",
      y = "Number of genes"
    ) +
    r5_theme
}

plot_seed_12345_lfdr <- function() {
  ggplot2::ggplot(
    seed_12345_lfdr,
    ggplot2::aes(x = full_lfdr, y = thinned_lfdr)
  ) +
    ggplot2::geom_abline(
      intercept = 0,
      slope = 1,
      color = "grey45",
      linetype = 2,
      linewidth = 0.7
    ) +
    ggplot2::geom_point(
      color = r5_blue,
      alpha = 0.45,
      size = 1.25
    ) +
    ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = "Seed 12345: paired lfdr values",
      subtitle = "Both axes use the identical 6,362 sampled gene-variant pairs",
      x = "Full-data BF-updated lfdr on sampled pair",
      y = "Thinned-refit BF-updated lfdr"
    ) +
    r5_theme
}

plot_null_weight_stability <- function() {
  plot_data <- data.frame(
    Analysis = "One variant per gene",
    seed = seed_summary$seed,
    pi0 = seed_summary$thinned_pi0,
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(plot_data, ggplot2::aes(x = Analysis, y = pi0)) +
    ggplot2::geom_hline(
      yintercept = full_pi0,
      color = r5_orange,
      linetype = 2,
      linewidth = 0.8
    ) +
    ggplot2::geom_boxplot(
      width = 0.32,
      outlier.shape = NA,
      fill = scales::alpha(r5_blue, 0.25),
      color = r5_blue
    ) +
    ggplot2::geom_jitter(
      width = 0.12,
      height = 0,
      color = r5_blue,
      alpha = 0.65,
      size = 1.5
    ) +
    ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 0.002)) +
    ggplot2::labs(
      title = "BF-updated null-weight stability",
      subtitle = paste0(
        "100 thinning seeds; dashed line is the full-data estimate (",
        format_decimal(full_pi0, 4),
        ")"
      ),
      x = NULL,
      y = "Estimated null weight"
    ) +
    r5_theme
}

plot_paired_lfdr_agreement <- function() {
  plot_data <- rbind(
    data.frame(
      seed = seed_summary$seed,
      Metric = "Mean absolute lfdr difference",
      Value = seed_summary$mean_absolute_lfdr_difference
    ),
    data.frame(
      seed = seed_summary$seed,
      Metric = "Spearman lfdr correlation",
      Value = seed_summary$spearman_lfdr
    )
  )
  ggplot2::ggplot(plot_data, ggplot2::aes(x = Metric, y = Value)) +
    ggplot2::geom_boxplot(
      width = 0.32,
      outlier.shape = NA,
      fill = scales::alpha(r5_blue, 0.25),
      color = r5_blue
    ) +
    ggplot2::geom_jitter(
      width = 0.12,
      height = 0,
      color = r5_blue,
      alpha = 0.65,
      size = 1.35
    ) +
    ggplot2::facet_wrap(~ Metric, scales = "free_y") +
    ggplot2::labs(
      title = "Paired lfdr agreement across 100 thinnings",
      subtitle = "Each point compares full and refitted lfdr on the same sampled pairs",
      x = NULL,
      y = NULL
    ) +
    r5_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

plot_discovery_robustness <- function() {
  plot_data <- rbind(
    data.frame(
      seed = seed_summary$seed,
      Metric = "Reference calls on sampled pairs",
      Value = seed_summary$full_fdr_calls
    ),
    data.frame(
      seed = seed_summary$seed,
      Metric = "Thinned-refit calls / genes",
      Value = seed_summary$thinned_fdr_calls
    ),
    data.frame(
      seed = seed_summary$seed,
      Metric = "Conditional discovery Jaccard",
      Value = seed_summary$fdr_call_jaccard
    )
  )
  plot_data$Metric <- factor(
    plot_data$Metric,
    levels = c(
      "Reference calls on sampled pairs",
      "Thinned-refit calls / genes",
      "Conditional discovery Jaccard"
    )
  )
  metric_colors <- c(
    "Reference calls on sampled pairs" = r5_orange,
    "Thinned-refit calls / genes" = r5_blue,
    "Conditional discovery Jaccard" = r5_green
  )
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = Metric, y = Value, color = Metric, fill = Metric)
  ) +
    ggplot2::geom_boxplot(
      width = 0.32,
      outlier.shape = NA,
      alpha = 0.18
    ) +
    ggplot2::geom_jitter(
      width = 0.12,
      height = 0,
      alpha = 0.65,
      size = 1.3
    ) +
    ggplot2::facet_wrap(~ Metric, scales = "free_y", nrow = 1) +
    ggplot2::scale_color_manual(values = metric_colors) +
    ggplot2::scale_fill_manual(values = metric_colors) +
    ggplot2::labs(
      title = "Discovery robustness within matched sampled universes",
      subtitle = paste(
        "Counts are not compared directly with the complete",
        "1,009,173-pair discovery universe"
      ),
      x = NULL,
      y = NULL
    ) +
    r5_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )
}
