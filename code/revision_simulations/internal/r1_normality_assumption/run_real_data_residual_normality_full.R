#!/usr/bin/env Rscript

# Build a day-stratified cache for traditional real-data residual diagnostics.

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

parse_runner_arguments <- function(arguments) {
  values <- list(
    max_units = Inf,
    output_id = paste0(
      "r1_real_data_standardized_residual_normality_",
      "days1_3_9_v1"
    )
  )
  for (argument in arguments) {
    if (grepl("^--max-units=", argument)) {
      values$max_units <- as.integer(sub("^--max-units=", "", argument))
    } else if (grepl("^--output-id=", argument)) {
      values$output_id <- sub("^--output-id=", "", argument)
    } else {
      stop("Unknown runner argument: ", argument)
    }
  }
  if ((!is.infinite(values$max_units) &&
       (!is.finite(values$max_units) || values$max_units < 1L)) ||
      !grepl("^[A-Za-z0-9._-]+$", values$output_id)) {
    stop("The runner arguments are invalid.")
  }
  values
}

summarize_file <- function(path, label) {
  information <- file.info(path)
  data.frame(
    input = label,
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    size_bytes = as.numeric(information$size),
    modified_at = format(information$mtime, tz = "UTC", usetz = TRUE),
    md5 = unname(tools::md5sum(path)),
    stringsAsFactors = FALSE
  )
}

runner_arguments <- parse_runner_arguments(commandArgs(trailingOnly = TRUE))
workflowr_root <- find_workflowr_root()
project_root <- normalizePath(
  file.path(workflowr_root, ".."),
  winslash = "/",
  mustWork = TRUE
)
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "r1_normality_assumption",
  "real_data_residual_normality_helpers.R"
))

normality_configuration_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "r1_normality_assumption_seed12345_t5_v1",
  "configuration.rds"
)
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
required_paths <- c(
  normality_configuration_path,
  expression_path,
  pc_path
)
if (any(!file.exists(required_paths))) {
  stop(
    "At least one required input is missing: ",
    paste(required_paths[!file.exists(required_paths)], collapse = ", ")
  )
}

normality_configuration <- readRDS(normality_configuration_path)
formal_r1_path <- normality_configuration$source_path
if (!file.exists(formal_r1_path)) {
  formal_output_suffix <- sub(
    "^.*[/\\\\]output[/\\\\]",
    "",
    normality_configuration$source_path
  )
  formal_r1_path <- file.path(
    workflowr_root,
    "output",
    formal_output_suffix
  )
}
if (!file.exists(formal_r1_path)) {
  stop("The formal R1 source object is missing.")
}

message("Loading the fixed one-variant-per-gene R1 subset.")
formal_r1 <- readRDS(formal_r1_path)
if (!is.matrix(formal_r1$genotype) ||
    nrow(formal_r1$genotype) != 19L ||
    ncol(formal_r1$genotype) != 6362L ||
    is.null(rownames(formal_r1$genotype)) ||
    is.null(colnames(formal_r1$genotype)) ||
    anyDuplicated(rownames(formal_r1$genotype)) ||
    anyDuplicated(colnames(formal_r1$genotype))) {
  stop("The fixed R1 genotype subset is invalid.")
}

all_unit_keys <- colnames(formal_r1$genotype)
all_gene_ids <- sub("_.*$", "", all_unit_keys)
all_variant_ids <- sub("^[^_]+_", "", all_unit_keys)
if (anyDuplicated(all_gene_ids) || any(!nzchar(all_variant_ids))) {
  stop("The formal R1 subset is not one unique variant per gene.")
}
selected_count <- min(length(all_unit_keys), runner_arguments$max_units)
selected_indices <- seq_len(selected_count)
unit_keys <- all_unit_keys[selected_indices]
gene_ids <- all_gene_ids[selected_indices]
variant_ids <- all_variant_ids[selected_indices]
genotype_matrix <- formal_r1$genotype[, selected_indices, drop = FALSE]

message("Loading expression and time-specific principal components.")
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
  stop("The expression or PC input identifiers are invalid.")
}
gene_rows <- match(gene_ids, expression_data$Gene_id)
if (anyNA(gene_rows)) {
  stop("At least one selected R1 gene is absent from the expression matrix.")
}

selected_units <- data.frame(
  unit_index = selected_indices,
  unit_key = unit_keys,
  gene_id = gene_ids,
  variant_id = variant_ids,
  stringsAsFactors = FALSE
)
analysis_days <- c(1L, 3L, 9L)
residual_rows <- vector("list", length(analysis_days))
fit_status_rows <- vector("list", length(analysis_days))
analysis_start <- proc.time()[["elapsed"]]

for (day_index in seq_along(analysis_days)) {
  day <- analysis_days[day_index]
  day_label <- paste("Day", day)
  message("Fitting ", format(selected_count, big.mark = ","), " units at ", day_label, ".")

  sample_ids <- grep(
    paste0("_", day, "$"),
    names(expression_data)[-1L],
    value = TRUE
  )
  donors <- sub(paste0("_", day, "$"), "", sample_ids)
  pc_rows <- match(sample_ids, pc_data$Sample_id)
  expression_columns <- match(sample_ids, names(expression_data))
  genotype_rows <- match(donors, rownames(genotype_matrix))
  if (length(sample_ids) != 19L ||
      anyNA(pc_rows) ||
      anyNA(expression_columns) ||
      anyNA(genotype_rows) ||
      anyDuplicated(donors) ||
      !setequal(donors, rownames(genotype_matrix))) {
    stop(day_label, " samples do not align across expression, PC, and genotype inputs.")
  }

  pc_matrix <- as.matrix(
    pc_data[pc_rows, paste0("PC", 1:5), drop = FALSE]
  )
  expression_matrix <- as.matrix(
    expression_data[gene_rows, expression_columns, drop = FALSE]
  )
  day_genotype <- genotype_matrix[genotype_rows, , drop = FALSE]
  storage.mode(pc_matrix) <- "double"
  storage.mode(expression_matrix) <- "double"
  storage.mode(day_genotype) <- "double"

  day_residual_rows <- vector("list", selected_count)
  day_fit_status <- vector("list", selected_count)
  for (unit_offset in seq_len(selected_count)) {
    design <- cbind(
      intercept = 1,
      genotype = day_genotype[, unit_offset],
      PC1 = pc_matrix[, 1L],
      PC2 = pc_matrix[, 2L],
      PC3 = pc_matrix[, 3L],
      PC4 = pc_matrix[, 4L],
      PC5 = pc_matrix[, 5L]
    )
    fitted_residuals <- fit_standardized_residuals(
      expression = expression_matrix[unit_offset, ],
      design = design
    )
    if (is.null(fitted_residuals)) {
      day_fit_status[[unit_offset]] <- data.frame(
        day = day,
        day_label = day_label,
        unit_index = selected_indices[unit_offset],
        unit_key = unit_keys[unit_offset],
        valid_fit = FALSE,
        exclusion_reason = "nonfinite or rank-deficient regression",
        stringsAsFactors = FALSE
      )
      next
    }

    day_residual_rows[[unit_offset]] <- data.frame(
      day = day,
      day_label = day_label,
      unit_index = selected_indices[unit_offset],
      unit_key = unit_keys[unit_offset],
      gene_id = gene_ids[unit_offset],
      variant_id = variant_ids[unit_offset],
      donor = donors,
      fitted_expression = fitted_residuals$fitted,
      raw_residual = fitted_residuals$residual,
      leverage = fitted_residuals$leverage,
      standardized_residual = fitted_residuals$standardized_residual,
      residual_df = fitted_residuals$residual_df,
      stringsAsFactors = FALSE
    )
    day_fit_status[[unit_offset]] <- data.frame(
      day = day,
      day_label = day_label,
      unit_index = selected_indices[unit_offset],
      unit_key = unit_keys[unit_offset],
      valid_fit = TRUE,
      exclusion_reason = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  residual_rows[[day_index]] <- do.call(
    rbind,
    day_residual_rows[!vapply(day_residual_rows, is.null, logical(1L))]
  )
  fit_status_rows[[day_index]] <- do.call(rbind, day_fit_status)
}

standardized_residuals <- do.call(rbind, residual_rows)
fit_status <- do.call(rbind, fit_status_rows)
rownames(standardized_residuals) <- NULL
rownames(fit_status) <- NULL
if (nrow(standardized_residuals) == 0L ||
    any(!is.finite(standardized_residuals$standardized_residual))) {
  stop("No valid finite standardized residuals were produced.")
}

qq_probability <- c(
  seq(0.001, 0.010, by = 0.001),
  seq(0.015, 0.985, by = 0.005),
  seq(0.990, 0.999, by = 0.001)
)
residual_qq <- do.call(rbind, lapply(analysis_days, function(day) {
  day_residual <- standardized_residuals$standardized_residual[
    standardized_residuals$day == day
  ]
  data.frame(
    day = day,
    day_label = paste("Day", day),
    probability = qq_probability,
    normal_reference_quantile = stats::qnorm(qq_probability),
    empirical_quantile = as.numeric(stats::quantile(
      day_residual,
      probs = qq_probability,
      names = FALSE,
      type = 8
    )),
    stringsAsFactors = FALSE
  )
}))

histogram_breaks <- seq(-5, 5, by = 0.10)
if (any(standardized_residuals$standardized_residual <= min(histogram_breaks) |
        standardized_residuals$standardized_residual >= max(histogram_breaks))) {
  stop("At least one standardized residual falls outside the fixed histogram range.")
}
residual_histogram <- do.call(rbind, lapply(analysis_days, function(day) {
  day_residual <- standardized_residuals$standardized_residual[
    standardized_residuals$day == day
  ]
  histogram <- graphics::hist(
    day_residual,
    breaks = histogram_breaks,
    plot = FALSE,
    right = FALSE
  )
  data.frame(
    day = day,
    day_label = paste("Day", day),
    bin_left = head(histogram$breaks, -1L),
    bin_right = tail(histogram$breaks, -1L),
    bin_midpoint = histogram$mids,
    count = histogram$counts,
    density = histogram$density,
    stringsAsFactors = FALSE
  )
}))

residual_day_summary <- do.call(rbind, lapply(analysis_days, function(day) {
  day_residual <- standardized_residuals$standardized_residual[
    standardized_residuals$day == day
  ]
  day_fit_status <- fit_status[fit_status$day == day, , drop = FALSE]
  data.frame(
    day = day,
    day_label = paste("Day", day),
    n_selected_units = selected_count,
    n_valid_fits = sum(day_fit_status$valid_fit),
    n_excluded_fits = sum(!day_fit_status$valid_fit),
    n_residuals = length(day_residual),
    residual_df = 12L,
    mean = mean(day_residual),
    median = stats::median(day_residual),
    sd = stats::sd(day_residual),
    iqr = stats::IQR(day_residual),
    skewness = calculate_skewness(day_residual),
    excess_kurtosis = calculate_excess_kurtosis(day_residual),
    maximum_absolute = max(abs(day_residual)),
    stringsAsFactors = FALSE
  )
}))

tail_levels <- c(0.05, 0.01, 0.001)
residual_tail_summary <- do.call(rbind, lapply(analysis_days, function(day) {
  day_residual <- standardized_residuals$standardized_residual[
    standardized_residuals$day == day
  ]
  do.call(rbind, lapply(tail_levels, function(nominal_rate) {
    cutoff <- stats::qnorm(1 - nominal_rate / 2)
    observed_rate <- mean(abs(day_residual) > cutoff)
    data.frame(
      day = day,
      day_label = paste("Day", day),
      n_residuals = length(day_residual),
      nominal_two_sided_rate = nominal_rate,
      normal_cutoff = cutoff,
      observed_two_sided_rate = observed_rate,
      observed_to_nominal_ratio = observed_rate / nominal_rate,
      stringsAsFactors = FALSE
    )
  }))
}))

spotcheck_keys <- unique(standardized_residuals$unit_key)[
  seq_len(min(10L, length(unique(standardized_residuals$unit_key))))
]
spotcheck_residuals <- standardized_residuals[
  standardized_residuals$unit_key %in% spotcheck_keys,
  ,
  drop = FALSE
]

input_provenance <- do.call(rbind, list(
  summarize_file(
    normality_configuration_path,
    "R1 normality-pilot configuration"
  ),
  summarize_file(formal_r1_path, "Formal R1 seed-12345 source"),
  summarize_file(expression_path, "Quantile-normalized expression"),
  summarize_file(pc_path, "Time-specific principal components")
))
elapsed_seconds <- proc.time()[["elapsed"]] - analysis_start
configuration <- list(
  experiment = paste(
    "Traditional standardized-residual normality diagnostic",
    "for the fixed R1 one-variant-per-gene subset"
  ),
  analysis_days = analysis_days,
  regression_formula = paste(
    "expression ~ genotype + PC1 + PC2 + PC3 + PC4 + PC5"
  ),
  residual_definition = "e_i / {s * sqrt(1 - h_ii)}",
  residual_reference = "Standard Normal, used as a conventional visual diagnostic",
  interpretation = paste(
    "Descriptive diagnostic only; repeated donors and within-fit residual",
    "dependence preclude a naive pooled normality p-value."
  ),
  n_available_units = length(all_unit_keys),
  n_selected_units = selected_count,
  n_donors_per_day = 19L,
  n_parameters = 7L,
  residual_df = 12L,
  qq_probability = qq_probability,
  histogram_breaks = histogram_breaks,
  input_provenance = input_provenance,
  elapsed_seconds = elapsed_seconds,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  r_version = R.version.string
)

expected_residual_rows <- sum(fit_status$valid_fit) * 19L
validation <- data.frame(
  check = c(
    "one unique selected variant per gene",
    "three requested days are present",
    "nineteen donors per valid regression",
    "residual df equals twelve",
    "residual-row denominator is exact",
    "all standardized residuals are finite",
    "QQ groups share one probability grid",
    "histogram counts equal residual denominators",
    "tail rates lie in the unit interval"
  ),
  pass = c(
    !anyDuplicated(selected_units$gene_id) &&
      !anyDuplicated(selected_units$unit_key),
    setequal(unique(standardized_residuals$day), analysis_days),
    all(table(
      standardized_residuals$day,
      standardized_residuals$unit_key
    )[table(
      standardized_residuals$day,
      standardized_residuals$unit_key
    ) > 0] == 19L),
    all(standardized_residuals$residual_df == 12L),
    nrow(standardized_residuals) == expected_residual_rows,
    all(is.finite(standardized_residuals$standardized_residual)),
    length(unique(table(residual_qq$day))) == 1L,
    all(vapply(analysis_days, function(day) {
      sum(residual_histogram$count[residual_histogram$day == day]) ==
        sum(standardized_residuals$day == day)
    }, logical(1L))),
    all(residual_tail_summary$observed_two_sided_rate >= 0 &
          residual_tail_summary$observed_two_sided_rate <= 1)
  ),
  stringsAsFactors = FALSE
)
if (!all(validation$pass)) {
  stop(
    "At least one output validation failed: ",
    paste(validation$check[!validation$pass], collapse = "; ")
  )
}

output_parent <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal"
)
output_directory <- file.path(output_parent, runner_arguments$output_id)
if (file.exists(output_directory)) {
  stop("The versioned output directory already exists: ", output_directory)
}
staging_directory <- tempfile(
  pattern = paste0(".", runner_arguments$output_id, "-staging-"),
  tmpdir = output_parent
)
dir.create(staging_directory, recursive = TRUE, showWarnings = FALSE)

saveRDS(configuration, file.path(staging_directory, "configuration.rds"))
saveRDS(
  standardized_residuals,
  file.path(staging_directory, "standardized_residuals.rds"),
  compress = "xz"
)
utils::write.csv(
  selected_units,
  file.path(staging_directory, "selected_units.csv"),
  row.names = FALSE
)
utils::write.csv(
  fit_status,
  file.path(staging_directory, "fit_status.csv"),
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  residual_qq,
  file.path(staging_directory, "residual_qq.csv"),
  row.names = FALSE
)
utils::write.csv(
  residual_histogram,
  file.path(staging_directory, "residual_histogram.csv"),
  row.names = FALSE
)
utils::write.csv(
  residual_day_summary,
  file.path(staging_directory, "residual_day_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  residual_tail_summary,
  file.path(staging_directory, "residual_tail_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  spotcheck_residuals,
  file.path(staging_directory, "spotcheck_residuals.csv"),
  row.names = FALSE
)
utils::write.csv(
  validation,
  file.path(staging_directory, "validation.csv"),
  row.names = FALSE
)
if (!file.rename(staging_directory, output_directory)) {
  stop("Could not atomically publish the residual-normality cache.")
}

message(
  "Residual cache complete: ",
  format(nrow(standardized_residuals), big.mark = ","),
  " residuals from ",
  format(sum(fit_status$valid_fit), big.mark = ","),
  " valid regressions in ",
  format(round(elapsed_seconds, 2), nsmall = 2),
  " seconds."
)
message("Output: ", output_directory)
