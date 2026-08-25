#!/usr/bin/env Rscript

# Compare a cached unit-specific residual permutation with one new donor map.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(args, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(args[equals_hit[1L]], nchar(equals_prefix) + 1L))
  }
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1L] == length(args)) {
    return(default)
  }
  args[hit[1L] + 1L]
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "covariance_estimation",
  "donor_null_permutation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation",
  "selected_signal_genotype_permutation_helpers.R"
))

default_experiment_id <- paste0(
  "all_gene_random_variant_signal_stripped_unit_specific_",
  "residual_block_permutation_selection20260817_seed20260811"
)
experiment_id <- get_arg("--experiment-id", default_experiment_id)
pair_key <- get_arg("--pair-key", "ENSG00000140396_rs34280145")
alternate_seed <- as.integer(get_arg("--alternate-seed", "20260812"))
experiment_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal", experiment_id
)
if (!dir.exists(experiment_dir) || is.na(alternate_seed)) {
  stop("The experiment directory or alternate seed is invalid.")
}

selection <- utils::read.csv(
  file.path(experiment_dir, "selection.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
source_information <- utils::read.csv(
  file.path(experiment_dir, "source_information.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
effect_estimates <- readRDS(file.path(experiment_dir, "effect_estimates.rds"))
current_maps <- utils::read.csv(
  file.path(experiment_dir, "donor_permutation.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
unit_index <- match(pair_key, selection$pair_key)
if (is.na(unit_index)) {
  stop("The requested pair key was not found.")
}
unit <- selection[unit_index, , drop = FALSE]

source_paths <- setNames(source_information$path, source_information$role)
required_roles <- c(
  "expression_matrix", "time_specific_pc_data", "genotype_vcf"
)
if (!all(required_roles %in% names(source_paths))) {
  stop("The saved source information is incomplete.")
}
expression_data <- utils::read.csv(
  source_paths[["expression_matrix"]],
  sep = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
pc_data <- utils::read.delim(
  source_paths[["time_specific_pc_data"]],
  check.names = FALSE,
  stringsAsFactors = FALSE
)
dosage <- read_selected_vcf_dosages(
  source_paths[["genotype_vcf"]],
  unit$variant_id,
  chunk_size = 100000L
)
colnames(dosage) <- pair_key
vcf_donors <- rownames(dosage)
time_grid <- as.numeric(effect_estimates$time_grid)
expression_sample_ids <- names(expression_data)[-1L]

donor_observation_matrix <- vapply(time_grid, function(time_value) {
  paste0(vcf_donors, "_", time_value) %in% expression_sample_ids
}, logical(length(vcf_donors)))
rownames(donor_observation_matrix) <- vcf_donors
observation_patterns <- apply(
  donor_observation_matrix,
  1L,
  paste0,
  collapse = ""
)
alternate_map <- make_unit_specific_donor_block_permutations(
  donor_ids = vcf_donors,
  observation_patterns = observation_patterns,
  unit_keys = pair_key,
  seed = alternate_seed
)
current_map <- current_maps[current_maps$unit_key == pair_key, , drop = FALSE]
current_map <- current_map[
  match(vcf_donors, current_map$target_donor),
  ,
  drop = FALSE
]
alternate_map <- alternate_map[
  match(vcf_donors, alternate_map$target_donor),
  ,
  drop = FALSE
]

fit_map <- function(donor_map) {
  beta <- raw_se <- numeric(length(time_grid))
  residual_df <- integer(length(time_grid))
  gene_row <- match(unit$gene_id, expression_data$Gene_id)
  if (is.na(gene_row)) {
    stop("The selected gene is missing from the expression matrix.")
  }
  for (time_index in seq_along(time_grid)) {
    time_value <- time_grid[time_index]
    sample_ids <- grep(
      paste0("_", time_value, "$"),
      expression_sample_ids,
      value = TRUE
    )
    donors <- sub(paste0("_", time_value, "$"), "", sample_ids)
    expression <- matrix(
      as.numeric(expression_data[
        gene_row,
        match(sample_ids, names(expression_data))
      ]),
      ncol = 1L,
      dimnames = list(donors, pair_key)
    )
    genotype <- dosage[donors, , drop = FALSE]
    pc_rows <- match(sample_ids, pc_data$Sample_id)
    covariates <- as.matrix(
      pc_data[pc_rows, paste0("PC", 1:5), drop = FALSE]
    )
    storage.mode(covariates) <- "double"
    rownames(covariates) <- donors
    source_donors <- donor_map$source_donor[
      match(donors, donor_map$target_donor)
    ]
    source_rows <- match(source_donors, donors)
    if (anyNA(source_rows)) {
      stop("The donor map is incompatible with observed donors.")
    }
    fit <- make_signal_stripped_residual_block_null(
      expression = expression,
      genotype = genotype,
      covariates = covariates,
      source_rows = source_rows
    )
    beta[time_index] <- fit$null_fit$beta
    raw_se[time_index] <- fit$null_fit$standard_error
    residual_df[time_index] <- fit$null_fit$residual_df
  }
  adjusted_se <- as.numeric(convert_raw_to_original_t_adjusted_se(
    matrix(beta, nrow = 1L),
    matrix(raw_se, nrow = 1L),
    residual_df
  ))
  list(beta = beta, adjusted_se = adjusted_se, residual_df = residual_df)
}

current_refit <- fit_map(current_map)
alternate_refit <- fit_map(alternate_map)
maximum_current_beta_difference <- max(abs(
  current_refit$beta - effect_estimates$permuted_beta[unit_index, ]
))
maximum_current_se_difference <- max(abs(
  current_refit$adjusted_se -
    effect_estimates$permuted_adjusted_se[unit_index, ]
))
if (maximum_current_beta_difference > 1e-12 ||
    maximum_current_se_difference > 1e-12) {
  stop("The current permutation reconstruction did not match the cache.")
}

make_rows <- function(label, fit, seed) {
  data.frame(
    permutation = label,
    seed = seed,
    time = time_grid,
    estimate = fit$beta,
    standard_error = fit$adjusted_se,
    lower = fit$beta - 1.96 * fit$adjusted_se,
    upper = fit$beta + 1.96 * fit$adjusted_se,
    stringsAsFactors = FALSE
  )
}
plot_data <- rbind(
  make_rows(
    "Current permutation used in the experiment",
    current_refit,
    20260811L
  ),
  make_rows(
    "New single-unit permutation",
    alternate_refit,
    alternate_seed
  )
)
plot_data$interval_excludes_zero <-
  plot_data$lower > 0 | plot_data$upper < 0
plot_data$permutation <- factor(
  plot_data$permutation,
  levels = unique(plot_data$permutation)
)

output_stem <- paste0("unit_", pair_key, "_alternate_permutation")
output_csv <- file.path(experiment_dir, paste0(output_stem, "_intervals.csv"))
utils::write.csv(plot_data, output_csv, row.names = FALSE)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("Package ggplot2 is required.")
}
plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = time, y = estimate)
) +
  ggplot2::geom_hline(
    yintercept = 0,
    color = "grey45",
    linewidth = 0.45,
    linetype = "dashed"
  ) +
  ggplot2::geom_line(color = "grey55", linewidth = 0.55) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = lower, ymax = upper),
    width = 0.16,
    color = "#0072B2",
    linewidth = 0.65
  ) +
  ggplot2::geom_point(
    ggplot2::aes(shape = interval_excludes_zero),
    color = "#0072B2",
    fill = "white",
    stroke = 0.9,
    size = 2.6
  ) +
  ggplot2::facet_wrap(~permutation, ncol = 1, scales = "fixed") +
  ggplot2::scale_shape_manual(
    values = c(`FALSE` = 21, `TRUE` = 19),
    labels = c(`FALSE` = "Includes zero", `TRUE` = "Excludes zero"),
    name = "95% interval"
  ) +
  ggplot2::scale_x_continuous(
    breaks = seq(0, 15, by = 1),
    minor_breaks = NULL
  ) +
  ggplot2::labs(
    title = "The apparent trajectory changes under a new donor permutation",
    subtitle = paste0(
      gsub("_", "-", pair_key),
      "; all inputs except the donor map are fixed"
    ),
    x = "Time",
    y = "Estimated genotype effect",
    caption = paste0(
      "Intervals are pointwise estimate +/- 1.96 x adjusted SE. ",
      "New map uses single-unit seed ", alternate_seed, "."
    )
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 14),
    plot.subtitle = ggplot2::element_text(color = "grey25"),
    strip.background = ggplot2::element_rect(fill = "grey94", color = NA),
    strip.text = ggplot2::element_text(face = "bold", hjust = 0),
    axis.text.x = ggplot2::element_text(size = 9),
    legend.position = "bottom",
    panel.spacing = grid::unit(0.8, "lines")
  )

output_png <- file.path(
  experiment_dir,
  "figures",
  paste0(output_stem, "_intervals.png")
)
output_pdf <- file.path(
  experiment_dir,
  "figures",
  paste0(output_stem, "_intervals.pdf")
)
ggplot2::ggsave(
  output_png,
  plot,
  width = 8.2,
  height = 7.2,
  dpi = 220,
  bg = "white"
)
ggplot2::ggsave(
  output_pdf,
  plot,
  width = 8.2,
  height = 7.2,
  device = grDevices::cairo_pdf
)

split_data <- split(plot_data, plot_data$permutation)
summary_data <- data.frame(
  permutation = names(split_data),
  slope = vapply(split_data, function(x) {
    unname(stats::coef(stats::lm(estimate ~ time, data = x))[2L])
  }, numeric(1)),
  positive_estimates = vapply(
    split_data,
    function(x) sum(x$estimate > 0),
    integer(1)
  ),
  intervals_excluding_zero = vapply(
    split_data,
    function(x) sum(x$interval_excludes_zero),
    integer(1)
  ),
  positive_intervals_excluding_zero = vapply(
    split_data,
    function(x) sum(x$lower > 0),
    integer(1)
  ),
  negative_intervals_excluding_zero = vapply(
    split_data,
    function(x) sum(x$upper < 0),
    integer(1)
  ),
  stringsAsFactors = FALSE
)
output_summary <- file.path(
  experiment_dir,
  paste0(output_stem, "_summary.csv")
)
utils::write.csv(summary_data, output_summary, row.names = FALSE)

cat("Current reconstruction maximum beta difference: ",
    format(maximum_current_beta_difference, scientific = TRUE), "\n", sep = "")
cat("Current reconstruction maximum SE difference: ",
    format(maximum_current_se_difference, scientific = TRUE), "\n", sep = "")
cat("Correlation between the two effect trajectories: ",
    format(stats::cor(current_refit$beta, alternate_refit$beta), digits = 6),
    "\n", sep = "")
cat("Donor maps identical: ",
    identical(current_map$source_donor, alternate_map$source_donor),
    "\n", sep = "")
print(summary_data)
cat("PNG: ", output_png, "\n", sep = "")
cat("PDF: ", output_pdf, "\n", sep = "")
cat("Intervals: ", output_csv, "\n", sep = "")
cat("Summary: ", output_summary, "\n", sep = "")
