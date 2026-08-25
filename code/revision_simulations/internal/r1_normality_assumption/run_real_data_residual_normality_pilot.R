#!/usr/bin/env Rscript

# Run a five-unit real-data residual-normality pilot at time zero.

find_workflowr_root <- function() {
  if (file.exists("analysis/index.Rmd")) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }
  if (file.exists("coderepo-local/analysis/index.Rmd")) {
    return(normalizePath(
      "coderepo-local",
      winslash = "/",
      mustWork = TRUE
    ))
  }
  stop("Could not find the workflowr repository root.")
}

require_pilot_package <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The ", package, " package is required for this pilot.")
  }
}

require_pilot_package("ggplot2")
require_pilot_package("patchwork")

workflowr_root <- find_workflowr_root()
project_root <- normalizePath(
  file.path(workflowr_root, ".."),
  winslash = "/",
  mustWork = TRUE
)

normality_cache <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r1_normality_assumption_seed12345_t5_v1"
)
configuration_path <- file.path(normality_cache, "configuration.rds")
expression_path <- file.path(
  project_root,
  "iPSC-data", "expression-data",
  "quantile_normalized_no_projection.txt"
)
pc_path <- file.path(
  project_root,
  "iPSC-data", "pc-data",
  "principal_components_10.txt"
)
required_inputs <- c(configuration_path, expression_path, pc_path)
if (any(!file.exists(required_inputs))) {
  stop(
    "At least one required input is missing: ",
    paste(required_inputs[!file.exists(required_inputs)], collapse = ", ")
  )
}

configuration <- readRDS(configuration_path)
formal_r1_path <- configuration$source_path
if (!file.exists(formal_r1_path)) {
  output_suffix <- sub(
    "^.*[/\\\\]output[/\\\\]",
    "",
    configuration$source_path
  )
  formal_r1_path <- file.path(workflowr_root, "output", output_suffix)
}
if (!file.exists(formal_r1_path)) {
  stop("The formal R1 source object is missing.")
}

message("Loading the formal one-variant-per-gene genotype subset.")
formal_r1 <- readRDS(formal_r1_path)
if (!is.matrix(formal_r1$genotype) ||
    nrow(formal_r1$genotype) != 19L ||
    ncol(formal_r1$genotype) != 6362L ||
    anyDuplicated(colnames(formal_r1$genotype))) {
  stop("The formal R1 genotype subset is invalid.")
}

message("Loading real expression and time-specific principal components.")
expression_data <- utils::read.csv(
  expression_path,
  sep = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
pc_data <- utils::read.delim(
  pc_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!identical(names(expression_data)[1L], "Gene_id") ||
    anyDuplicated(expression_data$Gene_id) ||
    anyDuplicated(pc_data$Sample_id) ||
    !all(c("Sample_id", paste0("PC", 1:5)) %in% names(pc_data))) {
  stop("The expression or PC data have invalid identifiers.")
}

# Time zero has 19 donors, so the full seven-parameter regression has 12 df.
time_value <- 0L
time_sample_ids <- grep(
  paste0("_", time_value, "$"),
  names(expression_data)[-1L],
  value = TRUE
)
time_donors <- sub(paste0("_", time_value, "$"), "", time_sample_ids)
pc_rows <- match(time_sample_ids, pc_data$Sample_id)
if (length(time_sample_ids) != 19L || anyNA(pc_rows) ||
    !setequal(time_donors, rownames(formal_r1$genotype))) {
  stop("The time-zero donors do not match the formal genotype subset.")
}
pc_matrix <- as.matrix(
  pc_data[pc_rows, paste0("PC", 1:5), drop = FALSE]
)
storage.mode(pc_matrix) <- "double"

unit_keys <- colnames(formal_r1$genotype)
unit_gene_ids <- sub("_.*$", "", unit_keys)
unit_variant_ids <- sub("^[^_]+_", "", unit_keys)
gene_rows <- match(unit_gene_ids, expression_data$Gene_id)

# Randomize only the candidate order; accept units using design validity alone.
set.seed(20260819)
candidate_indices <- sample.int(length(unit_keys))
selected_indices <- integer(0L)
for (candidate_index in candidate_indices) {
  if (is.na(gene_rows[candidate_index])) {
    next
  }
  genotype <- as.numeric(
    formal_r1$genotype[time_donors, candidate_index]
  )
  design <- cbind(
    intercept = 1,
    genotype = genotype,
    pc_matrix
  )
  if (qr(design)$rank == ncol(design)) {
    selected_indices <- c(selected_indices, candidate_index)
  }
  if (length(selected_indices) == 5L) {
    break
  }
}
if (length(selected_indices) != 5L) {
  stop("Could not select five full-rank real-data units.")
}

message("Fitting five time-zero regressions and extracting residuals.")
residual_rows <- vector("list", length(selected_indices))
selected_unit_rows <- vector("list", length(selected_indices))
for (selection_index in seq_along(selected_indices)) {
  unit_index <- selected_indices[selection_index]
  gene_id <- unit_gene_ids[unit_index]
  variant_id <- unit_variant_ids[unit_index]
  expression_values <- as.numeric(expression_data[
    gene_rows[unit_index],
    match(time_sample_ids, names(expression_data)),
    drop = TRUE
  ])
  regression_data <- data.frame(
    expression = expression_values,
    genotype = as.numeric(
      formal_r1$genotype[time_donors, unit_index]
    ),
    PC1 = pc_matrix[, 1L],
    PC2 = pc_matrix[, 2L],
    PC3 = pc_matrix[, 3L],
    PC4 = pc_matrix[, 4L],
    PC5 = pc_matrix[, 5L],
    donor = time_donors,
    stringsAsFactors = FALSE
  )
  fit <- stats::lm(
    expression ~ genotype + PC1 + PC2 + PC3 + PC4 + PC5,
    data = regression_data
  )
  if (fit$df.residual != 12L) {
    stop("A selected regression does not have residual df 12.")
  }
  external_residual <- stats::rstudent(fit)
  if (any(!is.finite(external_residual))) {
    stop("A selected regression produced invalid studentized residuals.")
  }
  residual_rows[[selection_index]] <- data.frame(
    unit_label = paste0("Unit ", selection_index),
    unit_index = unit_index,
    unit_key = unit_keys[unit_index],
    gene_id = gene_id,
    variant_id = variant_id,
    donor = time_donors,
    fitted_expression = stats::fitted(fit),
    raw_residual = stats::residuals(fit),
    externally_studentized_residual = external_residual,
    stringsAsFactors = FALSE
  )
  selected_unit_rows[[selection_index]] <- data.frame(
    unit_label = paste0("Unit ", selection_index),
    unit_index = unit_index,
    unit_key = unit_keys[unit_index],
    gene_id = gene_id,
    variant_id = variant_id,
    genotype_minor_allele_count = sum(regression_data$genotype),
    residual_df = fit$df.residual,
    external_reference_df = fit$df.residual - 1L,
    stringsAsFactors = FALSE
  )
}
residual_data <- do.call(rbind, residual_rows)
selected_units <- do.call(rbind, selected_unit_rows)
residual_data$unit_label <- factor(
  residual_data$unit_label,
  levels = selected_units$unit_label
)

# Build one QQ curve per unit against the exact external-residual t reference.
qq_data <- do.call(rbind, lapply(
  split(residual_data, residual_data$unit_label),
  function(unit_data) {
    n_residuals <- nrow(unit_data)
    data.frame(
      unit_label = unit_data$unit_label[1L],
      theoretical_quantile = stats::qt(
        stats::ppoints(n_residuals),
        df = 11L
      ),
      empirical_quantile = sort(
        unit_data$externally_studentized_residual
      )
    )
  }
))
qq_data$unit_label <- factor(
  qq_data$unit_label,
  levels = selected_units$unit_label
)

tail_cutoff_005 <- stats::qt(0.975, df = 11L)
tail_cutoff_001 <- stats::qt(0.995, df = 11L)
summarize_residuals <- function(data, label) {
  data.frame(
    unit_label = label,
    n = nrow(data),
    mean = mean(data$externally_studentized_residual),
    sd = stats::sd(data$externally_studentized_residual),
    maximum_absolute = max(abs(
      data$externally_studentized_residual
    )),
    tail_rate_0_05 = mean(abs(
      data$externally_studentized_residual
    ) > tail_cutoff_005),
    tail_rate_0_01 = mean(abs(
      data$externally_studentized_residual
    ) > tail_cutoff_001),
    stringsAsFactors = FALSE
  )
}
residual_summary <- do.call(rbind, c(
  lapply(
    split(residual_data, residual_data$unit_label),
    function(unit_data) {
      summarize_residuals(
        unit_data,
        as.character(unit_data$unit_label[1L])
      )
    }
  ),
  list(summarize_residuals(residual_data, "Pooled"))
))

pilot_theme <- function() {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey94"),
      plot.title.position = "plot"
    )
}

qq_plot <- ggplot2::ggplot(
  qq_data,
  ggplot2::aes(
    x = theoretical_quantile,
    y = empirical_quantile
  )
) +
  ggplot2::geom_abline(
    intercept = 0,
    slope = 1,
    color = "grey45",
    linewidth = 0.7,
    linetype = 2
  ) +
  ggplot2::geom_point(
    color = "#2166ac",
    size = 1.4,
    alpha = 0.80
  ) +
  ggplot2::facet_wrap(ggplot2::vars(unit_label), nrow = 1) +
  ggplot2::labs(
    title = "Externally studentized residual QQ plots",
    subtitle = paste(
      "Five outcome-independent real-data units at time 0;",
      "dashed line is the t(11) reference"
    ),
    x = "Theoretical t(11) quantile",
    y = "Empirical residual quantile"
  ) +
  pilot_theme()

reference_density <- data.frame(
  residual = seq(-5, 5, length.out = 1001L)
)
reference_density$density <- stats::dt(
  reference_density$residual,
  df = 11L
)
histogram_plot <- ggplot2::ggplot(
  residual_data,
  ggplot2::aes(x = externally_studentized_residual)
) +
  ggplot2::geom_histogram(
    ggplot2::aes(y = ggplot2::after_stat(density)),
    binwidth = 0.25,
    boundary = 0,
    fill = "#92c5de",
    color = "white",
    linewidth = 0.2
  ) +
  ggplot2::geom_line(
    data = reference_density,
    ggplot2::aes(x = residual, y = density),
    inherit.aes = FALSE,
    color = "#b2182b",
    linewidth = 0.9
  ) +
  ggplot2::coord_cartesian(xlim = c(-5, 5)) +
  ggplot2::labs(
    title = "Pooled residual distribution",
    subtitle = paste(
      "N = 95 residuals; red curve is the t(11) density;",
      "pooling is descriptive only"
    ),
    x = "Externally studentized residual",
    y = "Density"
  ) +
  pilot_theme()

pilot_figure <- qq_plot / histogram_plot +
  patchwork::plot_layout(heights = c(1.15, 0.85))

output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r1_real_data_residual_normality_pilot5_time0_v1"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
figure_path <- file.path(
  output_directory,
  "real_data_residual_normality_pilot.png"
)
ggplot2::ggsave(
  filename = figure_path,
  plot = pilot_figure,
  width = 12,
  height = 8,
  dpi = 180
)

utils::write.csv(
  selected_units,
  file.path(output_directory, "selected_units.csv"),
  row.names = FALSE
)
utils::write.csv(
  residual_data,
  file.path(output_directory, "residuals.csv"),
  row.names = FALSE
)
utils::write.csv(
  residual_summary,
  file.path(output_directory, "residual_summary.csv"),
  row.names = FALSE
)
result <- list(
  configuration = list(
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    selection_seed = 20260819L,
    n_units = 5L,
    time = time_value,
    n_donors = length(time_donors),
    regression = "Expression ~ genotype + PC1 + ... + PC5",
    residual = "Externally studentized residual",
    reference_distribution = "Student t with 11 degrees of freedom",
    scope = "Quick descriptive real-data pilot; not a formal normality test"
  ),
  selected_units = selected_units,
  residuals = residual_data,
  residual_summary = residual_summary,
  figure_path = figure_path
)
saveRDS(result, file.path(output_directory, "pilot_result.rds"))

print(selected_units, row.names = FALSE)
print(residual_summary, row.names = FALSE)
cat("Residual-normality pilot figure: ", figure_path, "\n", sep = "")
