#!/usr/bin/env Rscript

# Jointly evaluate middle-functional sensitivity to the evaluation-grid step
# and the posterior Monte Carlo sample size. Each sample-size run uses the same
# pair-specific seed and evaluates all grid resolutions from draws generated on
# the 0.025-day grid. A nested-prefix analysis from the 30,000-draw run is also
# included to isolate sample-size convergence from rerun-to-rerun randomness.

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
  absolute_samples <- abs(samples)
  statistic <- matrixStats::colMaxs(absolute_samples, rows = which(inside)) -
    matrixStats::colMaxs(absolute_samples, rows = which(!inside))
  as.integer(statistic <= 0)
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

jaccard_index <- function(x, y) {
  union_size <- length(union(x, y))
  if (union_size == 0L) return(1)
  length(intersect(x, y)) / union_size
}

summarize_convergence <- function(data, metadata, group_name, pair_ids = NULL) {
  if (!is.null(pair_ids)) {
    data <- data[data$pair_id %in% pair_ids, , drop = FALSE]
  }
  reference <- data[data$posterior_draws == max(data$posterior_draws), , drop = FALSE]
  reference_key <- paste(reference$method, reference$grid_step, reference$pair_id, sep = "|")
  reference_lfsr <- setNames(reference$lfsr, reference_key)
  current_key <- paste(data$method, data$grid_step, data$pair_id, sep = "|")
  data$reference_lfsr <- unname(reference_lfsr[current_key])
  data$absolute_difference <- abs(data$lfsr - data$reference_lfsr)

  split_data <- split(
    data,
    interaction(data$method, data$posterior_draws, data$grid_step, drop = TRUE)
  )
  result <- do.call(rbind, lapply(split_data, function(current) {
    data.frame(
      group = group_name,
      method = current$method[1],
      posterior_draws = current$posterior_draws[1],
      grid_step = current$grid_step[1],
      n_pairs = nrow(current),
      mean_absolute_difference_vs_30000 = mean(current$absolute_difference),
      median_absolute_difference_vs_30000 = median(current$absolute_difference),
      q90_absolute_difference_vs_30000 = unname(quantile(current$absolute_difference, 0.90)),
      maximum_absolute_difference_vs_30000 = max(current$absolute_difference),
      spearman_vs_30000 = suppressWarnings(cor(
        current$lfsr,
        current$reference_lfsr,
        method = "spearman"
      )),
      fraction_absolute_difference_gt_0p005 = mean(current$absolute_difference > 0.005),
      fraction_absolute_difference_gt_0p01 = mean(current$absolute_difference > 0.01),
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  result
}

plot_joint_convergence <- function(convergence, selection, output_path) {
  grid_order <- c("0.10", "0.05", "0.025")
  grid_colors <- c("0.10" = "#0072B2", "0.05" = "#D55E00", "0.025" = "#009E73")
  draw_values <- sort(unique(convergence$posterior_draws))

  png(output_path, width = 2300, height = 1800, res = 220, type = "cairo")
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(mfrow = c(2, 2), mar = c(4.5, 4.7, 3.5, 1.2), oma = c(0, 0, 2.0, 0))

  for (method_index in seq_along(c("direct_fixed_seed", "nested_prefix"))) {
    method <- c("direct_fixed_seed", "nested_prefix")[method_index]
    current <- convergence[
      convergence$group == "targeted_candidates" & convergence$method == method,
      ,
      drop = FALSE
    ]
    plot(
      NA,
      xlim = range(draw_values),
      ylim = c(0, max(current$mean_absolute_difference_vs_30000) * 1.10),
      log = "x",
      xaxt = "n",
      xlab = "Posterior Monte Carlo draws",
      ylab = "Mean absolute LFSR difference vs 30,000",
      main = if (method == "direct_fixed_seed") {
        "A. Fixed-seed reruns"
      } else {
        "B. Nested prefixes from 30,000 draws"
      }
    )
    axis(1, at = draw_values, labels = format(draw_values, big.mark = ",", scientific = FALSE))
    for (grid_name in grid_order) {
      rows <- current[current$grid_step == grid_name, , drop = FALSE]
      rows <- rows[order(rows$posterior_draws), , drop = FALSE]
      lines(
        rows$posterior_draws,
        rows$mean_absolute_difference_vs_30000,
        type = "b",
        pch = 16,
        lwd = 2,
        col = grid_colors[grid_name]
      )
    }
    legend(
      "topright",
      legend = grid_order,
      col = unname(grid_colors[grid_order]),
      pch = 16,
      lwd = 2,
      title = "Grid step",
      bty = "n",
      cex = 0.78
    )
  }

  direct_selection <- selection[selection$method == "direct_fixed_seed", , drop = FALSE]
  selected_ylim <- c(0, max(direct_selection$n_selected) * 1.12)
  plot(
    NA,
    xlim = range(draw_values),
    ylim = selected_ylim,
    log = "x",
    xaxt = "n",
    xlab = "Posterior Monte Carlo draws",
    ylab = "Selected Middle pairs",
    main = "C. Cumulative-FSR discoveries"
  )
  axis(1, at = draw_values, labels = format(draw_values, big.mark = ",", scientific = FALSE))
  for (grid_name in grid_order) {
    rows <- direct_selection[direct_selection$grid_step == grid_name, , drop = FALSE]
    rows <- rows[order(rows$posterior_draws), , drop = FALSE]
    lines(
      rows$posterior_draws,
      rows$n_selected,
      type = "b",
      pch = 16,
      lwd = 2,
      col = grid_colors[grid_name]
    )
  }

  plot(
    NA,
    xlim = range(draw_values),
    ylim = c(0, 1.03),
    log = "x",
    xaxt = "n",
    xlab = "Posterior Monte Carlo draws",
    ylab = "Jaccard index vs 30,000-draw set",
    main = "D. Discovery-set agreement"
  )
  axis(1, at = draw_values, labels = format(draw_values, big.mark = ",", scientific = FALSE))
  abline(h = 1, col = "#777777", lty = 3)
  for (grid_name in grid_order) {
    rows <- direct_selection[direct_selection$grid_step == grid_name, , drop = FALSE]
    rows <- rows[order(rows$posterior_draws), , drop = FALSE]
    lines(
      rows$posterior_draws,
      rows$jaccard_vs_30000,
      type = "b",
      pch = 16,
      lwd = 2,
      col = grid_colors[grid_name]
    )
  }

  mtext(
    "Middle functional: joint grid-resolution and Monte Carlo convergence",
    outer = TRUE,
    font = 2,
    cex = 1.15
  )
}

plot_representatives <- function(results, representatives, output_path) {
  grid_order <- c("0.10", "0.05", "0.025")
  grid_colors <- c("0.10" = "#0072B2", "0.05" = "#D55E00", "0.025" = "#009E73")
  draw_values <- sort(unique(results$posterior_draws))

  png(output_path, width = 2300, height = 1550, res = 220, type = "cairo")
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  }, add = TRUE)
  par(mfrow = c(2, 3), mar = c(4.4, 4.6, 3.4, 1.0), oma = c(0, 0, 2.0, 0))

  direct <- results[results$method == "direct_fixed_seed", , drop = FALSE]
  for (i in seq_len(nrow(representatives))) {
    current <- direct[direct$pair_id == representatives$pair_id[i], , drop = FALSE]
    y_limit <- range(c(0.05, current$lfsr - 2 * current$mcse, current$lfsr + 2 * current$mcse, 0))
    y_limit[2] <- y_limit[2] * 1.12 + 0.003
    plot(
      NA,
      xlim = range(draw_values),
      ylim = y_limit,
      log = "x",
      xaxt = "n",
      xlab = "Posterior draws",
      ylab = "Middle-functional LFSR",
      main = paste0(representatives$gene_symbol[i], " / ", representatives$variant_id[i])
    )
    axis(1, at = draw_values, labels = c("3k", "10k", "30k"))
    abline(h = 0.05, col = "#777777", lty = 2)
    for (grid_name in grid_order) {
      rows <- current[current$grid_step == grid_name, , drop = FALSE]
      rows <- rows[order(rows$posterior_draws), , drop = FALSE]
      lines(
        rows$posterior_draws,
        rows$lfsr,
        type = "b",
        pch = 16,
        lwd = 2,
        col = grid_colors[grid_name]
      )
    }
  }

  plot.new()
  legend(
    "center",
    legend = c(grid_order, "LFSR = 0.05 reference"),
    col = c(unname(grid_colors[grid_order]), "#777777"),
    pch = c(16, 16, 16, NA),
    lty = c(1, 1, 1, 2),
    lwd = c(2, 2, 2, 1),
    bty = "n",
    cex = 1.0
  )

  mtext(
    "One original Middle discovery per gene: fixed-seed Monte Carlo convergence",
    outer = TRUE,
    font = 2,
    cex = 1.10
  )
}

workflowr_root <- find_workflowr_root()
num_cores <- as.integer(get_arg("--num-cores", "2"))
seed <- as.integer(get_arg("--seed", "20260819"))
screen_lfsr <- as.numeric(get_arg("--screen-lfsr", "0.20"))
output_id <- get_arg("--output-id", "evaluation_grid_middle_grid_mc_convergence")
posterior_draw_values <- c(3000L, 10000L, 30000L)
grid_steps <- c(0.10, 0.05, 0.025)
grid_names <- c("0.10", "0.05", "0.025")

if (num_cores < 1L || num_cores > 2L || is.na(seed) ||
    !is.finite(screen_lfsr) || screen_lfsr <= 0 || screen_lfsr >= 1 ||
    !nzchar(output_id)) {
  stop("Invalid convergence arguments.")
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
original_selected <- middle[middle$selected_0p10, , drop = FALSE]
candidate <- middle[
  middle$pair_id %in% unique(c(
    original_selected$pair_id,
    middle$pair_id[middle$lfsr_0p05 <= screen_lfsr]
  )),
  ,
  drop = FALSE
]
candidate <- candidate[order(candidate$index), , drop = FALSE]
rownames(candidate) <- NULL

if (nrow(original_selected) != 24L || nrow(candidate) != 53L || anyDuplicated(candidate$pair_id)) {
  stop("The expected candidate universe was not recovered.")
}

gene_map <- readRDS(gene_map_path)
candidate$gene_id <- sub("_(rs[^_]+)$", "", candidate$pair_id)
candidate$variant_id <- sub("^.*_(rs[^_]+)$", "\\1", candidate$pair_id)
candidate$gene_symbol <- gene_map$hgnc_symbol[
  match(candidate$gene_id, gene_map$ensembl_gene_id)
]
candidate$gene_symbol[is.na(candidate$gene_symbol) | !nzchar(candidate$gene_symbol)] <-
  candidate$gene_id[is.na(candidate$gene_symbol) | !nzchar(candidate$gene_symbol)]
candidate$was_selected_0p10 <- candidate$pair_id %in% original_selected$pair_id

load(fit_path)
if (!exists("fash_fit1_update") || !inherits(fash_fit1_update, "fash")) {
  stop("The fitted-model file does not contain fash_fit1_update.")
}

fine_grid <- seq(0, 15, by = 0.025)
rows_by_grid <- lapply(grid_steps, function(step) grid_rows(fine_grid, step))
names(rows_by_grid) <- grid_names

sample_pair <- function(row_index) {
  pair <- candidate[row_index, , drop = FALSE]
  direct_rows <- vector("list", length(posterior_draw_values) * length(grid_steps))
  direct_index <- 0L
  full_indicators <- NULL

  for (posterior_draws in posterior_draw_values) {
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

    indicators <- lapply(seq_along(grid_steps), function(grid_index) {
      rows <- rows_by_grid[[grid_index]]
      middle_indicator(samples[rows, , drop = FALSE], fine_grid[rows])
    })
    names(indicators) <- grid_names

    for (grid_name in grid_names) {
      direct_index <- direct_index + 1L
      estimate <- mean(indicators[[grid_name]])
      direct_rows[[direct_index]] <- data.frame(
        pair_id = pair$pair_id,
        method = "direct_fixed_seed",
        posterior_draws = posterior_draws,
        grid_step = grid_name,
        lfsr = estimate,
        mcse = sqrt(estimate * (1 - estimate) / posterior_draws),
        stringsAsFactors = FALSE
      )
    }

    if (posterior_draws == max(posterior_draw_values)) {
      full_indicators <- indicators
    }
    rm(samples, indicators)
    gc()
  }

  set.seed(seed + pair$index + 10000000L)
  draw_order <- sample.int(max(posterior_draw_values))
  nested_rows <- vector("list", length(posterior_draw_values) * length(grid_steps))
  nested_index <- 0L
  for (posterior_draws in posterior_draw_values) {
    keep <- draw_order[seq_len(posterior_draws)]
    for (grid_name in grid_names) {
      nested_index <- nested_index + 1L
      estimate <- mean(full_indicators[[grid_name]][keep])
      nested_rows[[nested_index]] <- data.frame(
        pair_id = pair$pair_id,
        method = "nested_prefix",
        posterior_draws = posterior_draws,
        grid_step = grid_name,
        lfsr = estimate,
        mcse = sqrt(estimate * (1 - estimate) / posterior_draws),
        stringsAsFactors = FALSE
      )
    }
  }

  rbind(do.call(rbind, direct_rows), do.call(rbind, nested_rows))
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

results <- do.call(rbind, sampled)
results <- merge(
  results,
  candidate[, c("pair_id", "index", "gene_id", "gene_symbol", "variant_id", "was_selected_0p10")],
  by = "pair_id",
  all.x = TRUE,
  sort = FALSE
)
results <- results[order(results$method, results$posterior_draws, results$grid_step, results$index), , drop = FALSE]
rownames(results) <- NULL

expected_rows <- nrow(candidate) * 2L * length(posterior_draw_values) * length(grid_steps)
if (nrow(results) != expected_rows || anyNA(results$lfsr) ||
    any(results$lfsr < 0) || any(results$lfsr > 1)) {
  stop("The pair-level convergence results failed validation.")
}

baseline_by_grid <- list(
  `0.10` = middle$saved_lfsr_0p10,
  `0.05` = middle$lfsr_0p05,
  `0.025` = middle$lfsr_0p05
)

ranking_rows <- list()
selection_rows <- list()
ranking_index <- 0L
selection_index <- 0L
for (method in unique(results$method)) {
  for (posterior_draws in posterior_draw_values) {
    for (grid_name in grid_names) {
      current <- results[
        results$method == method &
          results$posterior_draws == posterior_draws &
          results$grid_step == grid_name,
        ,
        drop = FALSE
      ]
      lfsr <- baseline_by_grid[[grid_name]]
      lfsr[match(current$pair_id, middle$pair_id)] <- current$lfsr
      ranking <- calculate_cfsr(middle$pair_id, lfsr)
      ranking$selected <- ranking$cfsr <= 0.05
      ranking$method <- method
      ranking$posterior_draws <- posterior_draws
      ranking$grid_step <- grid_name
      ranking_index <- ranking_index + 1L
      ranking_rows[[ranking_index]] <- ranking

      selected_ids <- ranking$pair_id[ranking$selected]
      selection_index <- selection_index + 1L
      selection_rows[[selection_index]] <- data.frame(
        method = method,
        posterior_draws = posterior_draws,
        grid_step = grid_name,
        n_selected = length(selected_ids),
        selected_pair_ids = paste(selected_ids, collapse = ";"),
        stringsAsFactors = FALSE
      )
    }
  }
}

rankings <- do.call(rbind, ranking_rows)
selection <- do.call(rbind, selection_rows)
rownames(rankings) <- NULL
rownames(selection) <- NULL

selection$jaccard_vs_30000 <- NA_real_
selection$n_overlap_vs_30000 <- NA_integer_
for (i in seq_len(nrow(selection))) {
  reference <- selection[
    selection$method == selection$method[i] &
      selection$grid_step == selection$grid_step[i] &
      selection$posterior_draws == max(posterior_draw_values),
    ,
    drop = FALSE
  ]
  current_ids <- if (nzchar(selection$selected_pair_ids[i])) {
    strsplit(selection$selected_pair_ids[i], ";", fixed = TRUE)[[1]]
  } else {
    character()
  }
  reference_ids <- if (nzchar(reference$selected_pair_ids[1])) {
    strsplit(reference$selected_pair_ids[1], ";", fixed = TRUE)[[1]]
  } else {
    character()
  }
  selection$jaccard_vs_30000[i] <- jaccard_index(current_ids, reference_ids)
  selection$n_overlap_vs_30000[i] <- length(intersect(current_ids, reference_ids))
}

convergence <- rbind(
  summarize_convergence(results, candidate, "targeted_candidates"),
  summarize_convergence(
    results,
    candidate,
    "original_0p10_middle_discoveries",
    pair_ids = original_selected$pair_id
  )
)

representatives_source <- candidate[candidate$was_selected_0p10, , drop = FALSE]
representatives_source <- representatives_source[
  order(representatives_source$saved_lfsr_0p10, representatives_source$pair_id),
  ,
  drop = FALSE
]
representatives <- representatives_source[!duplicated(representatives_source$gene_symbol), , drop = FALSE]

joint_figure_path <- file.path(output_dir, "middle_grid_mc_joint_convergence.png")
representative_figure_path <- file.path(output_dir, "middle_grid_mc_representative_pairs.png")
plot_joint_convergence(convergence, selection, joint_figure_path)
plot_representatives(results, representatives, representative_figure_path)

write_csv(results, file.path(output_dir, "pair_lfsr_by_grid_and_mc_size.csv"))
write_csv(convergence, file.path(output_dir, "lfsr_convergence_summary.csv"))
write_csv(selection, file.path(output_dir, "hybrid_selection_convergence.csv"))
write_csv(rankings, file.path(output_dir, "hybrid_rankings_by_grid_and_mc_size.csv"))
write_csv(
  data.frame(
    seed = seed,
    num_cores = num_cores,
    screen_lfsr = screen_lfsr,
    n_candidates = nrow(candidate),
    n_original_selected = nrow(original_selected),
    grid_steps = paste(grid_names, collapse = ";"),
    posterior_draw_values = paste(posterior_draw_values, collapse = ";"),
    sampling_seconds = sampling_seconds,
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "runtime_and_configuration.csv")
)
saveRDS(
  list(
    fitted_model = normalizePath(fit_path),
    reclassification_input = normalizePath(reclassification_path),
    seed = seed,
    num_cores = num_cores,
    screen_lfsr = screen_lfsr,
    grid_steps = grid_steps,
    posterior_draw_values = posterior_draw_values,
    direct_method = paste(
      "Each sample-size run resets the same pair-specific seed and calls",
      "predict() on the 0.025-day grid."
    ),
    nested_method = paste(
      "The 30,000-draw indicators are deterministically shuffled per pair;",
      "the first 3,000, 10,000, and 30,000 indicators form nested prefixes."
    ),
    hybrid_note = paste(
      "Targeted candidates use the current high-draw estimate; non-targeted",
      "pairs retain the saved 0.10- or 0.05-grid LFSR."
    )
  ),
  file.path(output_dir, "configuration.rds")
)

print(convergence[convergence$group == "targeted_candidates", ])
print(selection)
