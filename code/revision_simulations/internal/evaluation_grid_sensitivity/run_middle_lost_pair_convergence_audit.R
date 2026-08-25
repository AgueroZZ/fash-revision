#!/usr/bin/env Rscript

# Audit whether middle-category losses after refining the evaluation grid are
# explained by Monte Carlo error or by a stable change in the functional.
# Posterior draws are generated on the 0.025-day grid and reused for nested
# 0.05-day and 0.10-day evaluations.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1] == length(args)) return(default)
  args[hit[1] + 1L]
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

grid_rows <- function(fine_grid, step) {
  scaled <- fine_grid / step
  which(abs(scaled - round(scaled)) < 1e-8)
}

middle_indicator <- function(samples, evaluation_grid) {
  tolerance <- sqrt(.Machine$double.eps)
  inside <- evaluation_grid >= 4 - tolerance & evaluation_grid <= 11 + tolerance
  if (!any(inside) || !any(!inside)) {
    stop("The grid does not support the middle functional.")
  }
  absolute_samples <- abs(samples)
  statistic <- matrixStats::colMaxs(absolute_samples, rows = which(inside)) -
    matrixStats::colMaxs(absolute_samples, rows = which(!inside))
  as.integer(statistic <= 0)
}

paired_mcse <- function(x, y) {
  sd(x - y) / sqrt(length(x))
}

calculate_cfsr <- function(pair_ids, lfsr) {
  ordered <- order(lfsr, pair_ids)
  data.frame(
    pair_id = pair_ids[ordered],
    lfsr = lfsr[ordered],
    cfsr = cumsum(lfsr[ordered]) / seq_along(ordered),
    rank = seq_along(ordered),
    stringsAsFactors = FALSE
  )
}

simulate_selection_counts <- function(lfsr, n_draws, n_replicates, seed) {
  set.seed(seed)
  vapply(seq_len(n_replicates), function(i) {
    estimated_lfsr <- rbinom(length(lfsr), size = n_draws, prob = lfsr) / n_draws
    ordered_lfsr <- sort(estimated_lfsr)
    sum(cumsum(ordered_lfsr) / seq_along(ordered_lfsr) <= 0.05)
  }, integer(1))
}

plot_convergence <- function(lost_results, selection_mc, output_path) {
  gene_levels <- unique(lost_results$gene_symbol)
  palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00")
  gene_colors <- setNames(rep(palette, length.out = length(gene_levels)), gene_levels)
  point_colors <- unname(gene_colors[lost_results$gene_symbol])

  png(output_path, width = 2300, height = 1850, res = 220, type = "cairo")
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(mfrow = c(2, 2), mar = c(4.6, 4.8, 3.5, 1.2), oma = c(0, 0, 2.0, 0))

  scatter_limit <- range(
    lost_results$high_draw_lfsr_0p10,
    lost_results$high_draw_lfsr_0p05,
    lost_results$high_draw_lfsr_0p025
  )
  scatter_limit <- c(0, max(scatter_limit) * 1.06)
  plot(
    lost_results$high_draw_lfsr_0p10,
    lost_results$high_draw_lfsr_0p05,
    pch = 16,
    col = point_colors,
    xlim = scatter_limit,
    ylim = scatter_limit,
    xlab = "LFSR on 0.10-day grid",
    ylab = "LFSR on 0.05-day grid",
    main = "A. Original Middle discoveries"
  )
  abline(a = 0, b = 1, col = "#555555", lty = 2)
  abline(h = 0.05, v = 0.05, col = "#999999", lty = 3)
  legend(
    "topleft",
    legend = gene_levels,
    col = unname(gene_colors[gene_levels]),
    pch = 16,
    bty = "n",
    cex = 0.72
  )

  plot(
    lost_results$high_draw_lfsr_0p05,
    lost_results$high_draw_lfsr_0p025,
    pch = 16,
    col = point_colors,
    xlim = scatter_limit,
    ylim = scatter_limit,
    xlab = "LFSR on 0.05-day grid",
    ylab = "LFSR on 0.025-day grid",
    main = "B. Refinement after 0.05"
  )
  abline(a = 0, b = 1, col = "#555555", lty = 2)
  abline(h = 0.05, v = 0.05, col = "#999999", lty = 3)

  differences <- list(
    `0.10 to 0.05` = lost_results$difference_0p05_minus_0p10,
    `0.05 to 0.025` = lost_results$difference_0p025_minus_0p05
  )
  boxplot(
    differences,
    col = c("#D55E00", "#0072B2"),
    border = "#555555",
    ylab = "Paired LFSR change",
    main = "C. Grid change compared with MC error",
    outline = FALSE
  )
  stripchart(
    differences,
    vertical = TRUE,
    method = "jitter",
    add = TRUE,
    pch = 16,
    col = adjustcolor("#222222", alpha.f = 0.55),
    jitter = 0.12
  )
  abline(h = 0, col = "#555555", lty = 2)
  legend(
    "topright",
    legend = sprintf(
      "Median paired MCSE: %.4f",
      median(lost_results$paired_mcse_0p05_minus_0p10)
    ),
    bty = "n",
    cex = 0.75
  )

  grid_labels <- c("0.10", "0.05", "0.025")
  selection_mc_3000 <- selection_mc[selection_mc$n_draws == 3000L, , drop = FALSE]
  block_by_grid <- split(
    selection_mc_3000$n_selected,
    factor(selection_mc_3000$grid_step, levels = grid_labels)
  )
  boxplot(
    block_by_grid,
    names = grid_labels,
    col = c("#56B4E9", "#E69F00", "#009E73"),
    border = "#555555",
    xlab = "Evaluation-grid step",
    ylab = "Selected Middle pairs",
    main = "D. Simulated selection at 3,000 draws",
    outline = FALSE
  )
  stripchart(
    block_by_grid,
    vertical = TRUE,
    method = "jitter",
    add = TRUE,
    pch = 16,
    col = adjustcolor("#222222", alpha.f = 0.06),
    jitter = 0.13
  )

  mtext(
    "Middle functional: high-draw convergence and Monte Carlo audit",
    outer = TRUE,
    font = 2,
    cex = 1.15
  )
}

plot_representative_pairs <- function(representatives, output_path) {
  grid_values <- c(0.10, 0.05, 0.025)
  grid_names <- c("0.10", "0.05", "0.025")
  colors <- c("#0072B2", "#D55E00", "#009E73")

  png(output_path, width = 2300, height = 1550, res = 220, type = "cairo")
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(mfrow = c(2, 3), mar = c(4.4, 4.6, 3.4, 1.0), oma = c(0, 0, 2.0, 0))

  for (i in seq_len(nrow(representatives))) {
    estimates <- c(
      representatives$high_draw_lfsr_0p10[i],
      representatives$high_draw_lfsr_0p05[i],
      representatives$high_draw_lfsr_0p025[i]
    )
    mcse <- c(
      representatives$mcse_0p10[i],
      representatives$mcse_0p05[i],
      representatives$mcse_0p025[i]
    )
    y_limit <- range(c(0.05, estimates - 2 * mcse, estimates + 2 * mcse, 0), finite = TRUE)
    y_limit[2] <- y_limit[2] * 1.12 + 0.005
    plot(
      seq_along(grid_values),
      estimates,
      type = "b",
      pch = 16,
      lwd = 2,
      col = "#444444",
      xaxt = "n",
      xlab = "Grid step",
      ylab = "Middle-functional LFSR",
      ylim = y_limit,
      main = paste0(representatives$gene_symbol[i], " / ", representatives$variant_id[i])
    )
    axis(1, at = seq_along(grid_values), labels = grid_names)
    arrows(
      seq_along(grid_values),
      estimates - 2 * mcse,
      seq_along(grid_values),
      estimates + 2 * mcse,
      angle = 90,
      code = 3,
      length = 0.045,
      col = colors,
      lwd = 1.7
    )
    points(seq_along(grid_values), estimates, pch = 16, col = colors, cex = 1.15)
    abline(h = 0.05, col = "#777777", lty = 2)
  }

  plot.new()
  legend(
    "center",
    legend = c("0.10 grid", "0.05 grid", "0.025 grid", "LFSR = 0.05"),
    col = c(colors, "#777777"),
    pch = c(16, 16, 16, NA),
    lty = c(NA, NA, NA, 2),
    bty = "n",
    cex = 1.0
  )

  mtext(
    "One original Middle discovery per gene: identical posterior fit, refined functional grid",
    outer = TRUE,
    font = 2,
    cex = 1.10
  )
}

workflowr_root <- find_workflowr_root()
posterior_draws <- as.integer(get_arg("--posterior-draws", "30000"))
block_size <- as.integer(get_arg("--block-size", "3000"))
num_cores <- as.integer(get_arg("--num-cores", "2"))
seed <- as.integer(get_arg("--seed", "20260819"))
screen_lfsr <- as.numeric(get_arg("--screen-lfsr", "0.20"))
selection_mc_replicates <- as.integer(get_arg("--selection-mc-replicates", "10000"))
output_id <- get_arg("--output-id", "evaluation_grid_middle_lost_pair_convergence_audit")

if (posterior_draws < 3000L || block_size < 100L ||
    posterior_draws %% block_size != 0L || num_cores < 1L || num_cores > 2L ||
    is.na(seed) || !is.finite(screen_lfsr) || screen_lfsr <= 0 ||
    screen_lfsr >= 1 || selection_mc_replicates < 1000L || !nzchar(output_id)) {
  stop("Invalid convergence-audit arguments.")
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

suppressPackageStartupMessages(library(fashr))
if (!requireNamespace("matrixStats", quietly = TRUE)) {
  stop("The matrixStats package is required.")
}

fit_path <- file.path(workflowr_root, "output", "dynamic_eQTL_real", "fash_fit1_update.RData")
gene_map_path <- file.path(workflowr_root, "output", "dynamic_eQTL_real", "cache_gene_map.rds")
reclassification_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "evaluation_grid_full_0p05_reclassification_pilot",
  "category_reclassification_by_pair.csv"
)
output_dir <- file.path(workflowr_root, "output", "revision_simulations", "internal", output_id)

if (!file.exists(fit_path) || !file.exists(gene_map_path) || !file.exists(reclassification_path)) {
  stop("At least one required input is missing.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

reclassification <- read.csv(reclassification_path, stringsAsFactors = FALSE, check.names = FALSE)
middle <- reclassification[reclassification$category == "middle", , drop = FALSE]
lost <- middle[middle$selected_0p10 & !middle$selected_0p05, , drop = FALSE]
candidate <- middle[middle$pair_id %in% unique(c(
  lost$pair_id,
  middle$pair_id[middle$lfsr_0p05 <= screen_lfsr]
)), , drop = FALSE]
candidate <- candidate[order(candidate$index), , drop = FALSE]
rownames(candidate) <- NULL

if (nrow(lost) != 24L || nrow(candidate) < nrow(lost) || anyDuplicated(candidate$pair_id)) {
  stop("The expected lost-pair or candidate universe was not recovered.")
}

gene_map <- readRDS(gene_map_path)
candidate$gene_id <- sub("_(rs[^_]+)$", "", candidate$pair_id)
candidate$variant_id <- sub("^.*_(rs[^_]+)$", "\\1", candidate$pair_id)
candidate$gene_symbol <- gene_map$hgnc_symbol[
  match(candidate$gene_id, gene_map$ensembl_gene_id)
]
candidate$gene_symbol[is.na(candidate$gene_symbol) | !nzchar(candidate$gene_symbol)] <-
  candidate$gene_id[is.na(candidate$gene_symbol) | !nzchar(candidate$gene_symbol)]
candidate$was_selected_0p10 <- candidate$pair_id %in% lost$pair_id

load(fit_path)
if (!exists("fash_fit1_update") || !inherits(fash_fit1_update, "fash")) {
  stop("The fitted-model file does not contain fash_fit1_update.")
}

fine_grid <- seq(0, 15, by = 0.025)
grid_steps <- c(0.10, 0.05, 0.025)
rows_by_grid <- lapply(grid_steps, function(step) grid_rows(fine_grid, step))
names(rows_by_grid) <- c("0.10", "0.05", "0.025")
expected_sizes <- as.integer(round(15 / grid_steps)) + 1L
if (!all(lengths(rows_by_grid) == expected_sizes)) {
  stop("The nested evaluation grids are invalid.")
}

sample_pair <- function(row_index) {
  pair <- candidate[row_index, , drop = FALSE]
  set.seed(seed + pair$index)
  samples <- predict(
    fash_fit1_update,
    index = pair$index,
    smooth_var = fine_grid,
    only.samples = TRUE,
    M = posterior_draws
  )
  if (!is.matrix(samples) || nrow(samples) != length(fine_grid) ||
      ncol(samples) != posterior_draws || any(!is.finite(samples))) {
    stop("Posterior sampling failed for ", pair$pair_id, ".")
  }

  indicators <- lapply(seq_along(grid_steps), function(i) {
    rows <- rows_by_grid[[i]]
    middle_indicator(samples[rows, , drop = FALSE], fine_grid[rows])
  })
  names(indicators) <- names(rows_by_grid)

  set.seed(seed + pair$index + 10000000L)
  draw_order <- sample.int(posterior_draws)
  block_id <- rep(seq_len(posterior_draws / block_size), each = block_size)
  block_results <- do.call(rbind, lapply(names(indicators), function(grid_name) {
    block_lfsr <- tapply(indicators[[grid_name]][draw_order], block_id, mean)
    data.frame(
      pair_id = pair$pair_id,
      block = as.integer(names(block_lfsr)),
      grid_step = grid_name,
      lfsr = as.numeric(block_lfsr),
      stringsAsFactors = FALSE
    )
  }))

  estimates <- vapply(indicators, mean, numeric(1))
  mcse <- sqrt(estimates * (1 - estimates) / posterior_draws)
  summary <- data.frame(
    pair_id = pair$pair_id,
    index = pair$index,
    high_draw_lfsr_0p10 = estimates["0.10"],
    mcse_0p10 = mcse["0.10"],
    high_draw_lfsr_0p05 = estimates["0.05"],
    mcse_0p05 = mcse["0.05"],
    high_draw_lfsr_0p025 = estimates["0.025"],
    mcse_0p025 = mcse["0.025"],
    difference_0p05_minus_0p10 = mean(indicators[["0.05"]] - indicators[["0.10"]]),
    paired_mcse_0p05_minus_0p10 = paired_mcse(indicators[["0.05"]], indicators[["0.10"]]),
    difference_0p025_minus_0p05 = mean(indicators[["0.025"]] - indicators[["0.05"]]),
    paired_mcse_0p025_minus_0p05 = paired_mcse(indicators[["0.025"]], indicators[["0.05"]]),
    discordance_0p10_vs_0p05 = mean(indicators[["0.10"]] != indicators[["0.05"]]),
    discordance_0p05_vs_0p025 = mean(indicators[["0.05"]] != indicators[["0.025"]]),
    stringsAsFactors = FALSE
  )
  list(summary = summary, blocks = block_results)
}

start_time <- proc.time()[["elapsed"]]
row_indices <- seq_len(nrow(candidate))
if (num_cores > 1L) {
  sampled <- parallel::mclapply(
    row_indices,
    sample_pair,
    mc.cores = num_cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
} else {
  sampled <- lapply(row_indices, sample_pair)
}
sampling_seconds <- proc.time()[["elapsed"]] - start_time

if (any(vapply(sampled, inherits, logical(1), "try-error"))) {
  stop("At least one posterior-sampling task failed.")
}

high_draw <- do.call(rbind, lapply(sampled, `[[`, "summary"))
block_lfsr <- do.call(rbind, lapply(sampled, `[[`, "blocks"))
results <- merge(candidate, high_draw, by = c("pair_id", "index"), all.x = TRUE, sort = FALSE)
results <- results[match(candidate$pair_id, results$pair_id), , drop = FALSE]

if (nrow(results) != nrow(candidate) || anyNA(results$high_draw_lfsr_0p025) ||
    any(results$high_draw_lfsr_0p025 < 0) || any(results$high_draw_lfsr_0p025 > 1)) {
  stop("The high-draw results failed validation.")
}

high_draw_by_grid <- list(
  `0.10` = setNames(results$high_draw_lfsr_0p10, results$pair_id),
  `0.05` = setNames(results$high_draw_lfsr_0p05, results$pair_id),
  `0.025` = setNames(results$high_draw_lfsr_0p025, results$pair_id)
)
baseline_by_grid <- list(
  `0.10` = middle$saved_lfsr_0p10,
  `0.05` = middle$lfsr_0p05,
  `0.025` = middle$lfsr_0p05
)

hybrid_rankings <- lapply(names(high_draw_by_grid), function(grid_name) {
  lfsr <- baseline_by_grid[[grid_name]]
  replacement <- high_draw_by_grid[[grid_name]]
  hit <- match(names(replacement), middle$pair_id)
  lfsr[hit] <- unname(replacement)
  ranking <- calculate_cfsr(middle$pair_id, lfsr)
  ranking$grid_step <- grid_name
  ranking$selected <- ranking$cfsr <= 0.05
  ranking
})
names(hybrid_rankings) <- names(high_draw_by_grid)
hybrid_ranking <- do.call(rbind, hybrid_rankings)
rownames(hybrid_ranking) <- NULL

for (grid_name in names(high_draw_by_grid)) {
  ranking <- hybrid_rankings[[grid_name]]
  mapping <- setNames(seq_len(nrow(ranking)), ranking$pair_id)
  result_rows <- mapping[results$pair_id]
  suffix <- gsub("\\.", "p", grid_name)
  results[[paste0("hybrid_cfsr_", suffix)]] <- ranking$cfsr[result_rows]
  results[[paste0("hybrid_selected_", suffix)]] <- ranking$selected[result_rows]
}

block_counts <- do.call(rbind, lapply(sort(unique(block_lfsr$block)), function(block) {
  do.call(rbind, lapply(names(high_draw_by_grid), function(grid_name) {
    current <- block_lfsr[
      block_lfsr$block == block & block_lfsr$grid_step == grid_name,
      ,
      drop = FALSE
    ]
    replacement <- setNames(current$lfsr, current$pair_id)
    lfsr <- baseline_by_grid[[grid_name]]
    lfsr[match(names(replacement), middle$pair_id)] <- unname(replacement)
    ranking <- calculate_cfsr(middle$pair_id, lfsr)
    selected <- ranking$cfsr <= 0.05
    data.frame(
      block = block,
      grid_step = grid_name,
      n_selected = sum(selected),
      n_selected_original_0p10_pairs = sum(
        ranking$pair_id[selected] %in% lost$pair_id
      ),
      stringsAsFactors = FALSE
    )
  }))
}))

lost_results <- results[results$was_selected_0p10, , drop = FALSE]
lost_results <- lost_results[order(lost_results$saved_lfsr_0p10, lost_results$pair_id), , drop = FALSE]
representatives <- lost_results[
  !duplicated(lost_results$gene_symbol),
  ,
  drop = FALSE
]

selected_summary <- do.call(rbind, lapply(names(hybrid_rankings), function(grid_name) {
  ranking <- hybrid_rankings[[grid_name]]
  selected <- ranking$cfsr <= 0.05
  data.frame(
    grid_step = grid_name,
    n_selected = sum(selected),
    n_selected_original_0p10_pairs = sum(ranking$pair_id[selected] %in% lost$pair_id),
    selected_pair_ids = paste(ranking$pair_id[selected], collapse = ";"),
    stringsAsFactors = FALSE
  )
}))

convergence_summary <- data.frame(
  group = c("all_targeted_candidates", "original_0p10_middle_discoveries"),
  n_pairs = c(nrow(results), nrow(lost_results)),
  mean_change_0p10_to_0p05 = c(
    mean(results$difference_0p05_minus_0p10),
    mean(lost_results$difference_0p05_minus_0p10)
  ),
  median_change_0p10_to_0p05 = c(
    median(results$difference_0p05_minus_0p10),
    median(lost_results$difference_0p05_minus_0p10)
  ),
  mean_change_0p05_to_0p025 = c(
    mean(results$difference_0p025_minus_0p05),
    mean(lost_results$difference_0p025_minus_0p05)
  ),
  median_change_0p05_to_0p025 = c(
    median(results$difference_0p025_minus_0p05),
    median(lost_results$difference_0p025_minus_0p05)
  ),
  median_paired_mcse_0p10_to_0p05 = c(
    median(results$paired_mcse_0p05_minus_0p10),
    median(lost_results$paired_mcse_0p05_minus_0p10)
  ),
  fraction_change_0p10_to_0p05_gt_2mcse = c(
    mean(abs(results$difference_0p05_minus_0p10) > 2 * results$paired_mcse_0p05_minus_0p10),
    mean(abs(lost_results$difference_0p05_minus_0p10) > 2 * lost_results$paired_mcse_0p05_minus_0p10)
  ),
  fraction_change_0p05_to_0p025_gt_2mcse = c(
    mean(abs(results$difference_0p025_minus_0p05) > 2 * results$paired_mcse_0p025_minus_0p05),
    mean(abs(lost_results$difference_0p025_minus_0p05) > 2 * lost_results$paired_mcse_0p025_minus_0p05)
  ),
  stringsAsFactors = FALSE
)

selection_mc <- do.call(rbind, lapply(seq_along(high_draw_by_grid), function(grid_index) {
  grid_name <- names(high_draw_by_grid)[grid_index]
  do.call(rbind, lapply(c(3000L, posterior_draws), function(n_draws) {
    data.frame(
      replicate = seq_len(selection_mc_replicates),
      grid_step = grid_name,
      n_draws = n_draws,
      n_selected = simulate_selection_counts(
        lfsr = unname(high_draw_by_grid[[grid_name]]),
        n_draws = n_draws,
        n_replicates = selection_mc_replicates,
        seed = seed + 1000L * grid_index + n_draws
      ),
      stringsAsFactors = FALSE
    )
  }))
}))

selection_mc_summary <- do.call(rbind, lapply(
  split(selection_mc, interaction(selection_mc$grid_step, selection_mc$n_draws, drop = TRUE)),
  function(data) {
    data.frame(
      grid_step = data$grid_step[1],
      n_draws = data$n_draws[1],
      n_replicates = nrow(data),
      mean_selected = mean(data$n_selected),
      median_selected = median(data$n_selected),
      q025_selected = unname(quantile(data$n_selected, 0.025)),
      q975_selected = unname(quantile(data$n_selected, 0.975)),
      probability_zero_selected = mean(data$n_selected == 0L),
      probability_at_most_two_selected = mean(data$n_selected <= 2L),
      stringsAsFactors = FALSE
    )
  }
))
rownames(selection_mc_summary) <- NULL

figure_path <- file.path(output_dir, "middle_lost_pair_convergence_audit.png")
representative_figure_path <- file.path(output_dir, "middle_representative_pair_convergence.png")
plot_convergence(lost_results, selection_mc, figure_path)
plot_representative_pairs(representatives, representative_figure_path)

write_csv(results, file.path(output_dir, "targeted_candidate_high_draw_lfsr.csv"))
write_csv(lost_results, file.path(output_dir, "original_middle_discovery_high_draw_lfsr.csv"))
write_csv(block_lfsr, file.path(output_dir, "candidate_block_lfsr.csv"))
write_csv(block_counts, file.path(output_dir, "hybrid_block_selection_counts.csv"))
write_csv(hybrid_ranking, file.path(output_dir, "hybrid_full_ranking_by_grid.csv"))
write_csv(selected_summary, file.path(output_dir, "hybrid_selection_summary.csv"))
write_csv(convergence_summary, file.path(output_dir, "convergence_summary.csv"))
write_csv(selection_mc, file.path(output_dir, "selection_mc_counts.csv"))
write_csv(selection_mc_summary, file.path(output_dir, "selection_mc_summary.csv"))
write_csv(
  data.frame(
    posterior_draws = posterior_draws,
    block_size = block_size,
    n_blocks = posterior_draws / block_size,
    num_cores = num_cores,
    screen_lfsr = screen_lfsr,
    selection_mc_replicates = selection_mc_replicates,
    n_candidates = nrow(candidate),
    n_original_middle_discoveries = nrow(lost),
    sampling_seconds = sampling_seconds,
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "runtime_and_configuration.csv")
)

saveRDS(
  list(
    fitted_model = normalizePath(fit_path),
    reclassification_input = normalizePath(reclassification_path),
    posterior_draws = posterior_draws,
    block_size = block_size,
    num_cores = num_cores,
    seed = seed,
    screen_lfsr = screen_lfsr,
    selection_mc_replicates = selection_mc_replicates,
    fine_grid_step = 0.025,
    compared_grid_steps = grid_steps,
    middle_window = c(4, 11),
    hybrid_note = paste(
      "High-draw LFSR replaces targeted candidates; non-targeted pairs retain",
      "their saved 0.10- or 0.05-grid LFSR. The 0.025 hybrid uses saved",
      "0.05-grid LFSR for non-targeted pairs."
    )
  ),
  file.path(output_dir, "configuration.rds")
)

print(convergence_summary)
print(selected_summary)
print(aggregate(n_selected ~ grid_step, block_counts, range))
print(selection_mc_summary)
