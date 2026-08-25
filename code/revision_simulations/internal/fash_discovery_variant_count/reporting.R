# Load, validate, summarize, and plot the gene-level tested-variant counts by
# BF-adjusted FASH-IWP1 discovery status.

find_workflowr_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "_workflowr.yml"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the workflowr project root.")
    }
    current <- parent
  }
}

require_reporting_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The ", package, " package is required to render this page.")
  }
  invisible(TRUE)
}

load_exact_object <- function(path, object_name) {
  if (!file.exists(path)) {
    stop("Required input file is missing: ", path)
  }
  environment <- new.env(parent = emptyenv())
  loaded_objects <- load(path, envir = environment)
  if (!identical(loaded_objects, object_name)) {
    stop("Expected only object ", object_name, " in ", path, ".")
  }
  environment[[object_name]]
}

format_integer <- function(value) {
  format(as.integer(value), big.mark = ",", scientific = FALSE)
}

format_decimal <- function(value, digits = 3L) {
  formatC(as.numeric(value), format = "f", digits = digits)
}

format_p_value <- function(value, digits = 3L) {
  formatC(as.numeric(value), format = "e", digits = digits)
}

render_internal_table <- function(table,
                                  align = NULL,
                                  minimum_width = "760px") {
  rendered <- knitr::kable(
    table,
    format = "html",
    align = align,
    escape = TRUE,
    table.attr = paste0(
      'class="table table-striped table-hover" style="min-width:',
      minimum_width,
      ';"'
    )
  )
  cat(
    '<div class="variant-count-table-scroll">\n',
    as.character(rendered),
    '\n</div>\n',
    sep = ""
  )
  invisible(rendered)
}

invisible(lapply(
  c("ggplot2", "knitr", "patchwork", "scales"),
  require_reporting_package
))

workflowr_root <- find_workflowr_root()
helper_path <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "one_variant_per_gene_refit",
  "one_variant_per_gene_refit_helpers.R"
)
fit_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "fash_fit1_update.RData"
)
count_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "one_variant_per_gene_variant_counts.csv"
)

if (!file.exists(helper_path)) {
  stop("The canonical cumulative-lfdr helper file is missing.")
}
source(helper_path)

input_paths <- c(
  "BF-adjusted FASH-IWP1 fit" = fit_path,
  "Per-gene tested-variant counts" = count_path
)
input_information <- file.info(input_paths)
input_provenance <- data.frame(
  Input = names(input_paths),
  Path = normalizePath(
    unname(input_paths),
    winslash = "/",
    mustWork = TRUE
  ),
  Bytes = as.numeric(input_information$size),
  Modified = format(input_information$mtime, tz = "UTC", usetz = TRUE),
  MD5 = unname(tools::md5sum(input_paths)),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

fash_fit <- load_exact_object(fit_path, "fash_fit1_update")
if (!is.list(fash_fit$fash_data) ||
    !is.list(fash_fit$fash_data$data_list) ||
    !is.numeric(fash_fit$lfdr)) {
  stop("The BF-adjusted FASH-IWP1 fit has an unexpected structure.")
}

pair_keys <- names(fash_fit$fash_data$data_list)
expected_pair_count <- 1009173L
expected_gene_count <- 6362L
expected_discovered_pair_count <- 9205L
expected_discovered_gene_count <- 1177L
analysis_alpha <- 0.05

if (length(pair_keys) != expected_pair_count ||
    length(fash_fit$lfdr) != expected_pair_count ||
    any(!nzchar(pair_keys)) ||
    anyDuplicated(pair_keys) ||
    any(!is.finite(fash_fit$lfdr)) ||
    any(fash_fit$lfdr < 0 | fash_fit$lfdr > 1)) {
  stop("The BF-adjusted fit failed pair-key or lfdr validation.")
}

pair_gene_ids <- sub("_.*$", "", pair_keys)
if (length(unique(pair_gene_ids)) != expected_gene_count) {
  stop("The BF-adjusted fit has an unexpected number of tested genes.")
}

discovered_pair_indices <- cumulative_fdr_calls(
  fash_fit$lfdr,
  alpha = analysis_alpha
)
discovered_gene_ids <- sort(unique(pair_gene_ids[discovered_pair_indices]))
if (length(discovered_pair_indices) != expected_discovered_pair_count ||
    length(discovered_gene_ids) != expected_discovered_gene_count) {
  stop("The BF-adjusted cumulative-lfdr discovery counts changed.")
}
rm(fash_fit)
gc(verbose = FALSE)

gene_variant_counts <- utils::read.csv(
  count_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_count_columns <- c("gene_id", "n_tested_variants")
if (!identical(names(gene_variant_counts), required_count_columns) ||
    nrow(gene_variant_counts) != expected_gene_count ||
    any(!nzchar(gene_variant_counts$gene_id)) ||
    anyDuplicated(gene_variant_counts$gene_id) ||
    !is.numeric(gene_variant_counts$n_tested_variants) ||
    any(!is.finite(gene_variant_counts$n_tested_variants)) ||
    any(gene_variant_counts$n_tested_variants < 1) ||
    any(gene_variant_counts$n_tested_variants !=
        as.integer(gene_variant_counts$n_tested_variants)) ||
    sum(gene_variant_counts$n_tested_variants) != expected_pair_count ||
    !setequal(gene_variant_counts$gene_id, unique(pair_gene_ids))) {
  stop("The per-gene tested-variant count cache failed validation.")
}

discovery_levels <- c("Not discovered", "FASH-discovered")
gene_variant_counts$`Discovery status` <- ifelse(
  gene_variant_counts$gene_id %in% discovered_gene_ids,
  "FASH-discovered",
  "Not discovered"
)
gene_variant_counts$`Discovery status` <- factor(
  gene_variant_counts$`Discovery status`,
  levels = discovery_levels
)
gene_variant_counts$log1p_n_tested_variants <- log1p(
  gene_variant_counts$n_tested_variants
)

if (!identical(
      as.integer(table(gene_variant_counts$`Discovery status`)),
      c(5185L, 1177L)
    )) {
  stop("The gene-level discovery groups have unexpected sizes.")
}

group_summary <- do.call(
  rbind,
  lapply(discovery_levels, function(group_label) {
    group_values <- gene_variant_counts$n_tested_variants[
      gene_variant_counts$`Discovery status` == group_label
    ]
    quartiles <- unname(stats::quantile(
      group_values,
      probs = c(0.25, 0.75),
      names = FALSE
    ))
    data.frame(
      `Discovery status` = group_label,
      Genes = as.integer(length(group_values)),
      Mean = mean(group_values),
      `Standard deviation` = stats::sd(group_values),
      Median = as.numeric(stats::median(group_values)),
      Q1 = as.numeric(quartiles[1]),
      Q3 = as.numeric(quartiles[2]),
      IQR = as.numeric(quartiles[2] - quartiles[1]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  })
)
rownames(group_summary) <- NULL

not_discovered_counts <- gene_variant_counts$n_tested_variants[
  gene_variant_counts$`Discovery status` == "Not discovered"
]
discovered_counts <- gene_variant_counts$n_tested_variants[
  gene_variant_counts$`Discovery status` == "FASH-discovered"
]
wilcoxon_test <- stats::wilcox.test(
  discovered_counts,
  not_discovered_counts,
  alternative = "greater",
  exact = FALSE,
  correct = TRUE
)

group_summary_display <- data.frame(
  `Discovery status` = group_summary$`Discovery status`,
  Genes = vapply(group_summary$Genes, format_integer, character(1)),
  Mean = format_decimal(group_summary$Mean, digits = 2L),
  `Standard deviation` = format_decimal(
    group_summary$`Standard deviation`,
    digits = 2L
  ),
  Median = format_decimal(group_summary$Median, digits = 0L),
  `IQR (Q1-Q3)` = paste0(
    format_decimal(group_summary$Q1, digits = 0L),
    "-",
    format_decimal(group_summary$Q3, digits = 0L)
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

wilcoxon_test_table <- data.frame(
  Method = unname(wilcoxon_test$method),
  `Null hypothesis` = "The two group distributions are equal",
  Alternative = "FASH-discovered values tend to be larger",
  Calculation = "Normal approximation",
  `Continuity correction` = "Yes",
  `W statistic` = formatC(
    as.numeric(wilcoxon_test$statistic),
    format = "f",
    digits = 1L,
    big.mark = ","
  ),
  `p-value` = format_p_value(wilcoxon_test$p.value, digits = 3L),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

relative_input_paths <- sub(
  paste0("^", workflowr_root, "/"),
  "",
  input_provenance$Path
)
input_provenance_table <- data.frame(
  Input = input_provenance$Input,
  Path = relative_input_paths,
  Bytes = vapply(input_provenance$Bytes, format_integer, character(1)),
  Modified = input_provenance$Modified,
  MD5 = input_provenance$MD5,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

group_colors <- c(
  "Not discovered" = "#B7BCC5",
  "FASH-discovered" = "#D55E00"
)
group_outline_colors <- c(
  "Not discovered" = "#6F7680",
  "FASH-discovered" = "#A94700"
)

plot_variant_count_by_discovery <- function() {
  x_axis_labels <- c(
    "Not discovered" = paste0(
      "Not discovered\n(n = ",
      format_integer(group_summary$Genes[1]),
      " genes)"
    ),
    "FASH-discovered" = paste0(
      "FASH-discovered\n(n = ",
      format_integer(group_summary$Genes[2]),
      " genes)"
    )
  )
  ggplot2::ggplot(
    gene_variant_counts,
    ggplot2::aes(
      x = `Discovery status`,
      y = log1p_n_tested_variants,
      fill = `Discovery status`
    )
  ) +
    ggplot2::geom_boxplot(
      width = 0.56,
      linewidth = 0.8,
      color = "#333333",
      outlier.shape = 21,
      outlier.size = 1.8,
      outlier.stroke = 0.35,
      outlier.color = "#555555",
      outlier.fill = "white",
      outlier.alpha = 0.75
    ) +
    ggplot2::scale_fill_manual(values = group_colors, guide = "none") +
    ggplot2::scale_x_discrete(labels = x_axis_labels) +
    ggplot2::scale_y_continuous(
      breaks = 1:8,
      expand = ggplot2::expansion(mult = c(0.06, 0.12))
    ) +
    ggplot2::labs(
      title = "Tested variants per gene by FASH discovery status",
      subtitle = paste0(
        "Raw medians: ",
        format_decimal(group_summary$Median[1], digits = 0L),
        " not discovered vs ",
        format_decimal(group_summary$Median[2], digits = 0L),
        " FASH-discovered"
      ),
      x = NULL,
      y = "log(1 + number of tested variants per gene)",
      caption = paste(
        "Boxplots use the standard 1.5 x IQR rule after log(1 + x)",
        "transformation; all transformed-scale outliers are shown."
      )
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 17),
      plot.subtitle = ggplot2::element_text(size = 12.5, color = "#555555"),
      plot.caption = ggplot2::element_text(
        size = 9.5,
        color = "#666666",
        hjust = 0
      ),
      axis.title.y = ggplot2::element_text(
        face = "bold",
        margin = ggplot2::margin(r = 10)
      ),
      axis.text.x = ggplot2::element_text(
        face = "bold",
        size = 12,
        margin = ggplot2::margin(t = 8)
      ),
      axis.text.y = ggplot2::element_text(color = "#333333"),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(16, 22, 14, 16)
    )
}

plot_variant_count_histogram <- function() {
  ggplot2::ggplot(
    gene_variant_counts,
    ggplot2::aes(
      x = log1p_n_tested_variants,
      fill = `Discovery status`,
      color = `Discovery status`
    )
  ) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      binwidth = 0.2,
      boundary = 0,
      position = "identity",
      alpha = 0.45,
      linewidth = 0.4
    ) +
    ggplot2::scale_fill_manual(
      values = group_colors,
      name = "Discovery status"
    ) +
    ggplot2::scale_color_manual(
      values = group_outline_colors,
      guide = "none"
    ) +
    ggplot2::scale_x_continuous(
      breaks = 1:8,
      expand = ggplot2::expansion(mult = c(0.02, 0.04))
    ) +
    ggplot2::labs(
      title = "Overlapping histogram",
      subtitle = "Each group is normalized to density",
      x = "log(1 + number of tested variants per gene)",
      y = "Density",
      caption = paste(
        "Histogram bin width is 0.2 on the log(1 + x) scale;",
        "bars are overlaid with 45% opacity."
      )
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle = ggplot2::element_text(size = 12.5, color = "#555555"),
      plot.caption = ggplot2::element_text(
        size = 9.5,
        color = "#666666",
        hjust = 0
      ),
      axis.title = ggplot2::element_text(face = "bold"),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 10)),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 10)),
      axis.text = ggplot2::element_text(color = "#333333"),
      legend.position = "top",
      legend.justification = "left",
      legend.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(10, 22, 14, 16)
    )
}

plot_variant_count_distribution_panel <- function() {
  boxplot_panel <- plot_variant_count_by_discovery() +
    ggplot2::labs(caption = NULL)
  histogram_panel <- plot_variant_count_histogram()
  combined_plot <- patchwork::wrap_plots(
    boxplot_panel,
    histogram_panel,
    ncol = 1L,
    heights = c(1.05, 1)
  )
  attr(combined_plot, "panel_count") <- 2L
  combined_plot
}
