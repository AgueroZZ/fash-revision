#!/usr/bin/env Rscript

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

workflowr_root <- find_workflowr_root()
original_working_directory <- getwd()
on.exit(setwd(original_working_directory), add = TRUE)
setwd(workflowr_root)

source(file.path(
  "code",
  "revision_simulations",
  "internal",
  "fash_discovery_variant_count",
  "reporting.R"
))

stopifnot(
  nrow(gene_variant_counts) == 6362L,
  sum(gene_variant_counts$n_tested_variants) == 1009173L,
  length(discovered_pair_indices) == 9205L,
  length(discovered_gene_ids) == 1177L,
  identical(
    as.character(group_summary$`Discovery status`),
    c("Not discovered", "FASH-discovered")
  ),
  identical(group_summary$Genes, c(5185L, 1177L)),
  isTRUE(all.equal(
    group_summary$Mean,
    c(155.303568, 173.257434),
    tolerance = 1e-6
  )),
  isTRUE(all.equal(group_summary$`Standard deviation`, c(
    71.856459,
    67.224394
  ), tolerance = 1e-6)),
  identical(group_summary$Median, c(147, 164)),
  identical(group_summary$Q1, c(107, 125)),
  identical(group_summary$Q3, c(193, 207)),
  identical(group_summary$IQR, c(86, 82)),
  isTRUE(all.equal(
    unname(wilcoxon_test$statistic),
    3549728.5,
    tolerance = 0
  )),
  isTRUE(all.equal(
    wilcoxon_test$p.value,
    9.696674e-19,
    tolerance = 1e-6
  )),
  identical(wilcoxon_test$alternative, "greater"),
  identical(
    unname(wilcoxon_test$method),
    "Wilcoxon rank sum test with continuity correction"
  ),
  nrow(group_summary_display) == 2L,
  nrow(wilcoxon_test_table) == 1L,
  identical(
    wilcoxon_test_table$Method,
    "Wilcoxon rank sum test with continuity correction"
  ),
  identical(wilcoxon_test_table$`W statistic`, "3,549,728.5"),
  nrow(input_provenance_table) == 2L,
  all(file.exists(input_provenance$Path)),
  inherits(plot_variant_count_by_discovery(), "ggplot"),
  inherits(plot_variant_count_histogram(), "ggplot"),
  inherits(plot_variant_count_distribution_panel(), "patchwork"),
  identical(
    attr(plot_variant_count_distribution_panel(), "panel_count"),
    2L
  )
)

plot_build <- ggplot2::ggplot_build(plot_variant_count_by_discovery())
stopifnot(
  length(plot_build$data) == 1L,
  nrow(plot_build$data[[1]]) == 2L,
  isTRUE(all.equal(
    sort(plot_build$data[[1]]$middle),
    sort(log1p(c(147, 164))),
    tolerance = 1e-12
  ))
)

histogram_build <- ggplot2::ggplot_build(plot_variant_count_histogram())
histogram_data <- histogram_build$data[[1]]
histogram_density_integrals <- tapply(
  histogram_data$density * (histogram_data$xmax - histogram_data$xmin),
  histogram_data$group,
  sum
)
stopifnot(
  length(histogram_build$data) == 1L,
  identical(length(histogram_density_integrals), 2L),
  isTRUE(all.equal(
    as.numeric(histogram_density_integrals),
    c(1, 1),
    tolerance = 1e-12
  ))
)

cat("FASH discovery variant-count reporting tests passed.\n")
