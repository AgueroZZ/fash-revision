#!/usr/bin/env Rscript

# Diagnose switch-functional calibration by reusing fitted FASH objects and
# evaluating posterior switch evidence over threshold and support-duration
# grids.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) return(".")
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find code/revision_simulations/shared/simulation_functions.R.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0 || hit[1] == length(args)) return(default)
  args[hit[1] + 1]
}

parse_numeric_grid <- function(x) {
  values <- suppressWarnings(as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1]])))
  if (length(values) == 0 || anyNA(values) || any(!is.finite(values))) {
    stop("Grid arguments must contain finite comma-separated numbers.")
  }
  sort(unique(values))
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

alpha <- as.numeric(get_arg("--alpha", "0.05"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
num_samples <- as.integer(get_arg("--num-samples", "3000"))
threshold_grid <- parse_numeric_grid(get_arg(
  "--threshold-grid",
  "0.25,0.375,0.5,0.625,0.75,0.875,1,1.125,1.25"
))
duration_grid <- parse_numeric_grid(get_arg(
  "--duration-grid",
  "0,0.25,0.5,0.75,1"
))
output_dir <- get_arg(
  "--output-dir",
  file.path("output", "revision_simulations", "diagnostics", "switch_functional")
)

if (!is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
    num_cores < 1 || num_samples < 100 ||
    any(threshold_grid <= 0) || any(duration_grid < 0)) {
  stop("Invalid switch-functional diagnostic arguments.")
}

raw_candidates <- list.files(
  file.path(workflowr_root, "output", "revision_simulations", "raw"),
  pattern = "^genotype_sparse_timed_cosine_one_two_peak.*seed12345[.]rds$",
  full.names = TRUE
)
if (length(raw_candidates) != 1L) {
  stop("Expected exactly one seed-12345 sparse timed cosine fitted object.")
}
out <- readRDS(raw_candidates)
evaluation_grid <- out$evaluation_grid
grid_step <- stats::median(diff(evaluation_grid))
observed_grid_index <- match(out$settings$time_grid, evaluation_grid)
true_curves <- out$true_beta_evaluation[out$unit_info$variant_id, , drop = FALSE]
if (!isTRUE(all.equal(colnames(true_curves), colnames(out$true_beta_evaluation))) ||
    any(!is.finite(true_curves)) || anyNA(observed_grid_index) ||
    !is.finite(grid_step) || grid_step <= 0) {
  stop("The fitted object does not contain valid refined-grid truth curves.")
}

switch_radius <- function(curve) {
  positive_radius <- max(c(curve[curve > 0], 0))
  negative_radius <- max(c(-curve[curve < 0], 0))
  min(positive_radius, negative_radius)
}

longest_run_duration <- function(indicator) {
  runs <- rle(as.logical(indicator))
  if (!any(runs$values)) {
    return(0)
  }
  max(runs$lengths[runs$values]) * grid_step
}

switch_duration <- function(curve, threshold) {
  positive_duration <- longest_run_duration(curve > threshold)
  negative_duration <- longest_run_duration(curve < -threshold)
  min(positive_duration, negative_duration)
}

true_radius <- apply(true_curves, 1, switch_radius)
true_radius_observed <- apply(
  true_curves[, observed_grid_index, drop = FALSE],
  1,
  switch_radius
)
true_duration <- lapply(threshold_grid, function(threshold) {
  apply(true_curves, 1, switch_duration, threshold = threshold)
})
names(true_duration) <- as.character(threshold_grid)

fits <- list(
  `FASH-IWP1-Raw` = out$fash_fits$fash_iwp1_raw,
  `FASH-IWP1-BF` = out$fash_fits$fash_iwp1_bf
)

sample_switch_summaries <- function(fit, method_index) {
  fdr_table <- get_fash_fdr_table(fit)
  selected_indices <- sort(unique(fdr_table$index[fdr_table$FDR <= alpha]))
  summaries <- parallel::mclapply(
    selected_indices,
    function(index) {
      set.seed(910000L + method_index * 10000L + index)
      samples <- predict(
        fit,
        index = index,
        smooth_var = evaluation_grid,
        only.samples = TRUE,
        M = num_samples
      )
      radii <- apply(samples, 2, switch_radius)
      radii_observed <- apply(
        samples[observed_grid_index, , drop = FALSE],
        2,
        switch_radius
      )
      durations <- vapply(
        threshold_grid,
        function(threshold) {
          apply(samples, 2, switch_duration, threshold = threshold)
        },
        numeric(num_samples)
      )
      list(
        index = index,
        radii = radii,
        radii_observed = radii_observed,
        durations = durations
      )
    },
    mc.cores = num_cores,
    mc.set.seed = FALSE
  )
  names(summaries) <- as.character(selected_indices)
  summaries
}

posterior_summaries <- Map(
  sample_switch_summaries,
  fits,
  seq_along(fits)
)

evaluate_rule <- function(method,
                          summaries,
                          threshold,
                          minimum_duration,
                          grid_scope = c("refined", "observed")) {
  grid_scope <- match.arg(grid_scope)
  indices <- as.integer(names(summaries))
  threshold_index <- match(threshold, threshold_grid)
  lfsr <- vapply(
    summaries,
    function(item) {
      if (grid_scope == "observed") {
        mean(item$radii_observed <= threshold)
      } else if (minimum_duration == 0) {
        mean(item$radii <= threshold)
      } else {
        mean(item$durations[, threshold_index] <= minimum_duration)
      }
    },
    numeric(1)
  )
  cfsr_table <- functional_cfsr_table(indices, lfsr)
  discoveries <- cfsr_table$index[cfsr_table$cfsr <= alpha]
  true_positive_rule <- if (grid_scope == "observed") {
    true_radius_observed > threshold
  } else if (minimum_duration == 0) {
    true_radius > threshold
  } else {
    true_duration[[as.character(threshold)]] > minimum_duration
  }
  false_discoveries <- sum(!true_positive_rule[discoveries])
  true_positives <- sum(true_positive_rule[discoveries])
  data.frame(
    method = method,
    alpha = alpha,
    grid_scope = grid_scope,
    effect_threshold = threshold,
    minimum_duration = minimum_duration,
    dynamic_discoveries = length(indices),
    functional_discoveries = length(discoveries),
    false_discoveries = false_discoveries,
    true_positives = true_positives,
    n_true_alternatives = sum(true_positive_rule),
    power = if (sum(true_positive_rule) == 0) NA_real_ else {
      true_positives / sum(true_positive_rule)
    },
    empirical_fsr = if (length(discoveries) == 0) 0 else {
      false_discoveries / length(discoveries)
    },
    mean_selected_lfsr = if (length(discoveries) == 0) 0 else {
      mean(cfsr_table$lfsr[cfsr_table$index %in% discoveries])
    },
    stringsAsFactors = FALSE
  )
}

diagnostic_rows <- do.call(rbind, lapply(names(posterior_summaries), function(method) {
  do.call(rbind, lapply(threshold_grid, function(threshold) {
    rbind(
      do.call(rbind, lapply(duration_grid, function(minimum_duration) {
        evaluate_rule(
          method = method,
          summaries = posterior_summaries[[method]],
          threshold = threshold,
          minimum_duration = minimum_duration,
          grid_scope = "refined"
        )
      })),
      evaluate_rule(
        method = method,
        summaries = posterior_summaries[[method]],
        threshold = threshold,
        minimum_duration = 0,
        grid_scope = "observed"
      )
    )
  }))
}))
rownames(diagnostic_rows) <- NULL

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  diagnostic_rows,
  file.path(output_dir, "switch_threshold_duration_diagnostic.csv"),
  row.names = FALSE
)

png(
  file.path(output_dir, "switch_threshold_duration_diagnostic.png"),
  width = 2200,
  height = 900,
  res = 180
)
old_par <- par(no.readonly = TRUE)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.2, 1.0))
styles <- revision_method_styles(names(fits), style_profile = "combined")
for (metric in c("power", "empirical_fsr")) {
  plot(
    NA,
    xlim = range(threshold_grid),
    ylim = c(0, 1),
    xlab = "Switch effect threshold",
    ylab = if (metric == "power") "Power" else "Empirical FSR",
    main = if (metric == "power") "Threshold sensitivity: power" else {
      "Threshold sensitivity: empirical FSR"
    }
  )
  grid(col = "gray90")
  if (metric == "empirical_fsr") abline(h = alpha, col = "gray35", lty = 3)
  rows <- diagnostic_rows[
    diagnostic_rows$minimum_duration == 0 &
      diagnostic_rows$grid_scope == "refined",
    ,
    drop = FALSE
  ]
  for (i in seq_along(styles$methods)) {
    method_rows <- rows[rows$method == styles$methods[i], , drop = FALSE]
    method_rows <- method_rows[order(method_rows$effect_threshold), , drop = FALSE]
    lines(
      method_rows$effect_threshold,
      method_rows[[metric]],
      col = styles$col[i],
      lty = styles$lty[i],
      lwd = 2.5
    )
    points(
      method_rows$effect_threshold,
      method_rows[[metric]],
      col = styles$col[i],
      pch = 19
    )
  }
  legend(
    if (metric == "power") "bottomleft" else "topright",
    legend = styles$methods,
    col = styles$col,
    lty = styles$lty,
    lwd = 2.5,
    bty = "n"
  )
}
par(old_par)
dev.off()

print(diagnostic_rows)
