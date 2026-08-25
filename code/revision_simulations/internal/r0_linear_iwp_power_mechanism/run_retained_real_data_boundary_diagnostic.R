#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr project root.")
}

summarize_numeric <- function(x, prefix) {
  values <- c(
    median = stats::median(x),
    p25 = stats::quantile(x, 0.25, names = FALSE),
    p75 = stats::quantile(x, 0.75, names = FALSE)
  )
  names(values) <- paste0(prefix, "_", names(values))
  values
}

select_cumulative_lfdr <- function(lfdr, alpha = 0.05) {
  ordering <- order(lfdr, seq_along(lfdr))
  cumulative_fdr <- cumsum(lfdr[ordering]) / seq_along(ordering)
  last_selected <- max(c(0L, which(cumulative_fdr <= alpha)))
  selected <- rep(FALSE, length(lfdr))
  if (last_selected > 0L) {
    selected[ordering[seq_len(last_selected)]] <- TRUE
  }
  selected
}

workflowr_root <- find_workflowr_root()
output_root <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal"
)
real_data_directory <- file.path(
  output_root,
  "fash_linear_real_data_ablation"
)
analysis_path <- file.path(real_data_directory, "analysis_cache.rds")
statistics_path <- file.path(real_data_directory, "sufficient_statistics.rds")
retained_linear_fit_path <- file.path(real_data_directory, "linear_fit_bf.rds")
matched_analysis_path <- file.path(
  output_root,
  "r0_linear_iwp_power_mechanism_current_matched_real_data",
  "analysis_cache.rds"
)
clean_confirmation_path <- file.path(
  output_root,
  "r0_linear_iwp_power_mechanism_clean_confirmation",
  "confirmation_summary.csv"
)

retained_analysis <- readRDS(analysis_path)
statistics_cache <- readRDS(statistics_path)
retained_linear_fit <- readRDS(retained_linear_fit_path)
matched_analysis <- readRDS(matched_analysis_path)
matched_fit <- matched_analysis$compact_bf
clean_confirmation <- utils::read.csv(
  clean_confirmation_path,
  stringsAsFactors = FALSE
)
statistics <- statistics_cache$statistics
retained_status <- retained_analysis$lfdr_scatter_all$discovery_status
iwp_selected <- retained_status %in% c("Current FASH only", "Both")
linear_selected <- select_cumulative_lfdr(matched_fit$lfdr, alpha = 0.05)
status <- factor(
  ifelse(
    iwp_selected & linear_selected,
    "Both",
    ifelse(
      iwp_selected,
      "Current FASH only",
      ifelse(linear_selected, "FASH-linear only", "Neither")
    )
  ),
  levels = c("Neither", "Current FASH only", "FASH-linear only", "Both")
)

if (!inherits(retained_linear_fit, "profiled_linear_fash") ||
    !inherits(matched_fit, "compact_linear_mixture_fash") ||
    !grepl(
      "one Gaussian linear slope slab",
      retained_analysis$configuration$comparator_model
    ) ||
    nrow(statistics) != length(status) ||
    nrow(statistics) != nrow(retained_analysis$lfdr_scatter_all) ||
    !identical(
      as.character(statistics$unit_id),
      as.character(matched_fit$unit_ids)
    )) {
  stop("The retained and current real-data artifacts are not aligned.")
}

determinant <- statistics$sum_w * statistics$sum_wxx -
  statistics$sum_wx^2
slope_hat <- (
  statistics$sum_w * statistics$sum_wxy -
    statistics$sum_wx * statistics$sum_wy
) / determinant
slope_variance <- statistics$sum_w / determinant
absolute_slope_z <- abs(slope_hat / sqrt(slope_variance))
linear_residual_quadratic <- statistics$sum_wyy - (
  statistics$sum_wxx * statistics$sum_wy^2 -
    2 * statistics$sum_wx * statistics$sum_wy * statistics$sum_wxy +
    statistics$sum_w * statistics$sum_wxy^2
) / determinant
linear_residual_rms_z <- sqrt(pmax(linear_residual_quadratic, 0) / 14)

category_levels <- levels(status)
category_summary <- do.call(rbind, lapply(category_levels, function(category) {
  selected <- status == category
  data.frame(
    category = category,
    pair_count = sum(selected),
    as.list(summarize_numeric(
      absolute_slope_z[selected],
      "absolute_slope_z"
    )),
    as.list(summarize_numeric(
      linear_residual_rms_z[selected],
      "linear_residual_rms_z"
    )),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}))
rownames(category_summary) <- NULL

gene_id <- sub("_.*$", "", statistics$unit_id)

summarize_multiplicity <- function(selected, method) {
  count_by_gene <- table(gene_id[selected])
  data.frame(
    method = method,
    pair_count = sum(selected),
    gene_count = length(count_by_gene),
    mean_pairs_per_discovered_gene = mean(count_by_gene),
    median_pairs_per_discovered_gene = stats::median(count_by_gene),
    p90_pairs_per_discovered_gene = stats::quantile(
      count_by_gene,
      0.90,
      names = FALSE
    ),
    stringsAsFactors = FALSE
  )
}

multiplicity_summary <- rbind(
  summarize_multiplicity(iwp_selected, "Retained FASH-IWP1"),
  summarize_multiplicity(linear_selected, "Retained FASH-linear")
)
pair_count_ratio <-
  multiplicity_summary$pair_count[[2L]] /
    multiplicity_summary$pair_count[[1L]]
gene_count_ratio <-
  multiplicity_summary$gene_count[[2L]] /
    multiplicity_summary$gene_count[[1L]]
multiplicity_ratio <-
  multiplicity_summary$mean_pairs_per_discovered_gene[[2L]] /
    multiplicity_summary$mean_pairs_per_discovered_gene[[1L]]
ratio_decomposition <- data.frame(
  quantity = c(
    "pair_count_ratio",
    "gene_count_ratio",
    "pairs_per_discovered_gene_ratio",
    "gene_ratio_times_multiplicity_ratio"
  ),
  value = c(
    pair_count_ratio,
    gene_count_ratio,
    multiplicity_ratio,
    gene_count_ratio * multiplicity_ratio
  ),
  stringsAsFactors = FALSE
)

provenance <- data.frame(
  artifact = c(
    "retained analysis cache",
    "retained profiled linear fit",
    "current matched linear diagnostic",
    "current real-data runner",
    "current internal page"
  ),
  path = normalizePath(c(
    analysis_path,
    retained_linear_fit_path,
    matched_analysis_path,
    file.path(
      workflowr_root,
      "code", "revision_simulations", "internal",
      "fash_linear_real_data_ablation",
      "run_fash_linear_real_data_ablation.R"
    ),
    file.path(
      workflowr_root,
      "analysis", "revision_internal_fash_linear_real_data_ablation.rmd"
    )
  )),
  modified = as.character(file.info(c(
    analysis_path,
    retained_linear_fit_path,
    matched_analysis_path,
    file.path(
      workflowr_root,
      "code", "revision_simulations", "internal",
      "fash_linear_real_data_ablation",
      "run_fash_linear_real_data_ablation.R"
    ),
    file.path(
      workflowr_root,
      "analysis", "revision_internal_fash_linear_real_data_ablation.rmd"
    )
  ))$mtime),
  retained_object_or_definition = c(
    retained_analysis$configuration$comparator_model,
    class(retained_linear_fit)[[1L]],
    class(matched_fit)[[1L]],
    "matched predictive-SD finite mixture",
    "matched predictive-SD finite mixture"
  ),
  stringsAsFactors = FALSE
)

validation <- data.frame(
  check = c(
    "retained comparator is a profiled single slab",
    "current comparator is the matched finite mixture",
    "current pair alignment is exact",
    "discovery counts reproduce current matched analysis",
    "ratio decomposition is exact",
    "diagnostics are finite"
  ),
  passed = c(
    inherits(retained_linear_fit, "profiled_linear_fash"),
    inherits(matched_fit, "compact_linear_mixture_fash"),
    identical(
      as.character(statistics$unit_id),
      as.character(matched_fit$unit_ids)
    ),
    sum(iwp_selected) == 9205L && sum(linear_selected) == 14900L,
    abs(pair_count_ratio - gene_count_ratio * multiplicity_ratio) <= 1e-12,
    all(is.finite(absolute_slope_z)) &&
      all(is.finite(linear_residual_rms_z))
  ),
  stringsAsFactors = FALSE
)
if (any(!validation$passed)) {
  print(validation)
  stop("The retained real-data boundary diagnostic failed validation.")
}

output_directory <- file.path(
  output_root,
  "r0_linear_iwp_power_mechanism_real_data_diagnostic"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

figure_path <- file.path(output_directory, "mechanism_summary.png")
grDevices::png(figure_path, width = 2400, height = 900, res = 180)
old_par <- graphics::par(no.readonly = TRUE)
graphics::par(mfrow = c(1, 3), mar = c(8, 4.5, 4.5, 1))

clean_plot <- clean_confirmation[
  clean_confirmation$candidate %in% c(
    "T16_homoskedastic_fixed",
    "T16_middle_precise_boundary_fixed",
    "T16_middle_precise_fixed",
    "T32_homoskedastic_gaussian_grid"
  ),
  ,
  drop = FALSE
]
clean_plot <- clean_plot[match(
  c(
    "T16_homoskedastic_fixed",
    "T16_middle_precise_boundary_fixed",
    "T16_middle_precise_fixed",
    "T32_homoskedastic_gaussian_grid"
  ),
  clean_plot$candidate
), , drop = FALSE]
power_matrix <- rbind(
  clean_plot$mean_power_linear,
  clean_plot$mean_power_iwp
)
positions <- graphics::barplot(
  power_matrix,
  beside = TRUE,
  col = c("#CC79A7", "#0072B2"),
  ylim = c(0, 0.62),
  names.arg = c(
    "T16 homo\nfixed 3.0",
    "T16 middle\nfixed 2.5",
    "T16 middle\nfixed 3.25",
    "T32 homo\nGaussian"
  ),
  ylab = "Mean BF power",
  main = "A. Held-out clean simulation",
  cex.names = 0.78
)
graphics::arrows(
  positions,
  rbind(
    clean_plot$mean_power_linear_lower,
    clean_plot$mean_power_iwp_lower
  ),
  positions,
  rbind(
    clean_plot$mean_power_linear_upper,
    clean_plot$mean_power_iwp_upper
  ),
  angle = 90,
  code = 3,
  length = 0.04
)
graphics::legend(
  "topright",
  legend = c("FASH-linear", "FASH-IWP1"),
  fill = c("#CC79A7", "#0072B2"),
  bty = "n",
  cex = 0.82
)

category_colors <- c("#999999", "#0072B2", "#CC79A7", "#009E73")
graphics::barplot(
  category_summary$absolute_slope_z_median,
  col = category_colors,
  names.arg = c("Neither", "IWP only", "Linear only", "Both"),
  las = 2,
  ylab = "Median absolute slope z",
  main = ""
)
graphics::mtext("B. Current real-data calls", side = 3, line = 1.2, font = 2)
graphics::abline(
  h = category_summary$absolute_slope_z_median[
    category_summary$category == "FASH-linear only"
  ],
  lty = 3,
  col = "gray30"
)

graphics::barplot(
  category_summary$linear_residual_rms_z_median,
  col = category_colors,
  names.arg = c("Neither", "IWP only", "Linear only", "Both"),
  las = 2,
  ylab = "Median linear-residual RMS z",
  main = ""
)
graphics::mtext("C. Departure from a line", side = 3, line = 1.2, font = 2)
graphics::par(old_par)
grDevices::dev.off()

analysis_cache <- list(
  configuration = list(
    scope = "descriptive diagnostic of retained real-data calls",
    residual_degrees_of_freedom = 14L,
    caveat = paste(
      "Real-data categories have unknown truth and dependent gene-variant",
      "pairs; no power or inferential significance is estimated."
    )
  ),
  category_summary = category_summary,
  multiplicity_summary = multiplicity_summary,
  ratio_decomposition = ratio_decomposition,
  provenance = provenance,
  validation = validation
)
saveRDS(analysis_cache, file.path(output_directory, "analysis_cache.rds"))
utils::write.csv(
  category_summary,
  file.path(output_directory, "category_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  multiplicity_summary,
  file.path(output_directory, "multiplicity_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  ratio_decomposition,
  file.path(output_directory, "ratio_decomposition.csv"),
  row.names = FALSE
)
utils::write.csv(
  provenance,
  file.path(output_directory, "provenance.csv"),
  row.names = FALSE
)
utils::write.csv(
  validation,
  file.path(output_directory, "validation.csv"),
  row.names = FALSE
)

print(category_summary)
print(multiplicity_summary)
print(ratio_decomposition)
print(provenance[, c("artifact", "modified", "retained_object_or_definition")])
print(validation)
message("Saved output to: ", normalizePath(output_directory))
