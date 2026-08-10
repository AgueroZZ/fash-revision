#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  hit <- match(flag, args)
  if (is.na(hit)) {
    return(default)
  }
  if (hit == length(args)) {
    stop("Missing value after ", flag, ".")
  }
  args[[hit + 1L]]
}

workflowr_root <- normalizePath(
  get_arg("--workflowr-root", "."),
  mustWork = TRUE
)
seed <- as.integer(get_arg("--seed", "12345"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
num_samples <- as.integer(get_arg("--num-samples", "3000"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
method <- get_arg("--method", "FASH-IWP1-BF")
candidate_scope <- get_arg("--candidate-scope", "fdr-selected")
output_dir <- get_arg(
  "--output-dir",
  file.path(
    "output",
    "revision_simulations",
    "diagnostics",
    "switch_miscalibration"
  )
)

if (!method %in% c("FASH-IWP1-Raw", "FASH-IWP1-BF") ||
    !candidate_scope %in% c("fdr-selected", "oracle-dynamic", "all") ||
    !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
    !is.finite(num_samples) || num_samples < 500 ||
    !is.finite(num_cores) || num_cores < 1) {
  stop("Invalid diagnostic arguments.")
}

setwd(workflowr_root)
source("code/revision_simulations/shared/simulation_functions.R")

raw_candidates <- list.files(
  file.path("output", "revision_simulations", "raw"),
  pattern = paste0(
    "^genotype_sparse_timed_cosine_one_two_peak.*seed",
    seed,
    "[.]rds$"
  ),
  full.names = TRUE
)
if (length(raw_candidates) != 1L) {
  stop("Expected exactly one fitted sparse timed cosine object for seed ", seed, ".")
}
out <- readRDS(raw_candidates)
fit <- switch(
  method,
  "FASH-IWP1-Raw" = out$fash_fits$fash_iwp1_raw,
  "FASH-IWP1-BF" = out$fash_fits$fash_iwp1_bf
)

evaluation_grid <- out$evaluation_grid
observed_grid <- out$settings$time_grid
observed_grid_index <- match(observed_grid, evaluation_grid)
if (anyNA(observed_grid_index)) {
  stop("Observed time points are missing from the evaluation grid.")
}

switch_radius <- function(curve) {
  positive_radius <- max(c(curve[curve > 0], 0))
  negative_radius <- max(c(-curve[curve < 0], 0))
  min(positive_radius, negative_radius)
}

fdr_table <- get_fash_fdr_table(fit)
dynamic_indices <- switch(
  candidate_scope,
  "fdr-selected" = sort(unique(fdr_table$index[fdr_table$FDR <= alpha])),
  "oracle-dynamic" = which(out$unit_info$effect_class == "dynamic_bspline"),
  "all" = seq_len(nrow(out$unit_info))
)
if (length(dynamic_indices) == 0) {
  stop("No dynamic discoveries at the requested alpha.")
}

posterior_summaries <- parallel::mclapply(
  dynamic_indices,
  function(index) {
    set.seed(880000L + seed + index)
    samples <- predict(
      fit,
      index = index,
      smooth_var = evaluation_grid,
      only.samples = TRUE,
      M = num_samples
    )
    sample_positive <- apply(samples, 2, max)
    sample_negative <- apply(-samples, 2, max)
    sample_radius <- pmin(sample_positive, sample_negative)
    posterior_mean <- rowMeans(samples)
    list(
      index = index,
      lfsr = mean(sample_radius <= 0.25),
      posterior_switch_probability = mean(sample_radius > 0.25),
      posterior_mean_radius = switch_radius(posterior_mean),
      sample_radius_median = stats::median(sample_radius),
      sample_radius_q025 = unname(stats::quantile(sample_radius, 0.025)),
      sample_radius_q975 = unname(stats::quantile(sample_radius, 0.975)),
      sample_positive_median = stats::median(sample_positive),
      sample_negative_median = stats::median(sample_negative),
      posterior_mean = posterior_mean,
      posterior_lower = apply(samples, 1, stats::quantile, probs = 0.025),
      posterior_upper = apply(samples, 1, stats::quantile, probs = 0.975)
    )
  },
  mc.cores = num_cores,
  mc.set.seed = FALSE
)
names(posterior_summaries) <- as.character(dynamic_indices)

per_unit <- do.call(rbind, lapply(posterior_summaries, function(item) {
  index <- item$index
  truth <- out$true_beta_evaluation[index, ]
  estimate <- out$eqtl_summary$beta_hat[index, ]
  posterior_null_weight <- fit$posterior_weights[index, "0"]
  data.frame(
    index = index,
    unit_id = out$unit_info$unit_id[index],
    effect_class = out$unit_info$effect_class[index],
    cell_id = out$unit_info$cell_id[index],
    spike_count = out$unit_info$spike_count[index],
    sign_pattern = out$unit_info$sign_pattern[index],
    true_switch = out$true_functionals[index, "switch"] > 0,
    true_radius = switch_radius(truth),
    observed_radius = switch_radius(estimate),
    posterior_mean_radius = item$posterior_mean_radius,
    posterior_switch_probability = item$posterior_switch_probability,
    functional_lfsr = item$lfsr,
    global_lfdr = fit$lfdr[index],
    posterior_null_weight = posterior_null_weight,
    sample_radius_median = item$sample_radius_median,
    sample_radius_q025 = item$sample_radius_q025,
    sample_radius_q975 = item$sample_radius_q975,
    sample_positive_median = item$sample_positive_median,
    sample_negative_median = item$sample_negative_median,
    stringsAsFactors = FALSE
  )
}))

cfsr_table <- functional_cfsr_table(
  per_unit$index,
  per_unit$functional_lfsr
)
per_unit$cfsr <- cfsr_table$cfsr[match(per_unit$index, cfsr_table$index)]
per_unit$functional_call <- per_unit$cfsr <= alpha
per_unit$false_functional_call <-
  per_unit$functional_call & !per_unit$true_switch

per_unit$truth_group <- ifelse(
  per_unit$effect_class != "dynamic_bspline",
  "first-stage dynamic null",
  ifelse(
    per_unit$true_switch,
    "opposite-sign two-peak",
    ifelse(per_unit$spike_count == 1, "one-peak", "same-sign two-peak")
  )
)

summarize_group <- function(rows) {
  data.frame(
    selected_dynamic = nrow(rows),
    switch_calls = sum(rows$functional_call),
    false_switch_calls = sum(rows$false_functional_call),
    call_rate = mean(rows$functional_call),
    median_observed_radius = stats::median(rows$observed_radius),
    median_posterior_mean_radius = stats::median(rows$posterior_mean_radius),
    median_posterior_switch_probability =
      stats::median(rows$posterior_switch_probability),
    median_functional_lfsr = stats::median(rows$functional_lfsr),
    median_global_lfdr = stats::median(rows$global_lfdr)
  )
}

group_summary <- do.call(
  rbind,
  lapply(split(per_unit, per_unit$truth_group), summarize_group)
)
group_summary$truth_group <- rownames(group_summary)
rownames(group_summary) <- NULL
group_summary <- group_summary[, c("truth_group", setdiff(
  names(group_summary),
  "truth_group"
))]

overall_summary <- data.frame(
  seed = seed,
  method = method,
  candidate_scope = candidate_scope,
  alpha = alpha,
  dynamic_discoveries = nrow(per_unit),
  switch_calls = sum(per_unit$functional_call),
  false_switch_calls = sum(per_unit$false_functional_call),
  true_switch_calls = sum(per_unit$functional_call & per_unit$true_switch),
  empirical_fsr = mean(
    per_unit$false_functional_call[per_unit$functional_call]
  ),
  estimated_fsr = mean(
    per_unit$functional_lfsr[per_unit$functional_call]
  ),
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
method_stem <- if (method == "FASH-IWP1-BF") "bf" else "raw"
scope_stem <- gsub("-", "_", candidate_scope)
stem <- paste0(
  "seed",
  seed,
  "_",
  method_stem,
  "_",
  scope_stem,
  "_alpha",
  gsub("[.]", "p", alpha)
)
write.csv(
  per_unit,
  file.path(output_dir, paste0(stem, "_per_unit.csv")),
  row.names = FALSE
)
write.csv(
  group_summary,
  file.path(output_dir, paste0(stem, "_group_summary.csv")),
  row.names = FALSE
)
write.csv(
  overall_summary,
  file.path(output_dir, paste0(stem, "_overall_summary.csv")),
  row.names = FALSE
)

false_rows <- per_unit[per_unit$false_functional_call, , drop = FALSE]
false_rows <- false_rows[order(false_rows$functional_lfsr), , drop = FALSE]
example_indices <- head(false_rows$index, 6)
if (length(example_indices) > 0) {
  png(
    file.path(output_dir, paste0(stem, "_false_switch_examples.png")),
    width = 2400,
    height = 1500,
    res = 200
  )
  old_par <- par(no.readonly = TRUE)
  par(mfrow = c(2, 3), mar = c(4, 4.2, 3.2, 1))
  for (index in example_indices) {
    item <- posterior_summaries[[as.character(index)]]
    ylim <- range(
      out$true_beta_evaluation[index, ],
      out$eqtl_summary$beta_hat[index, ],
      item$posterior_lower,
      item$posterior_upper
    )
    plot(
      evaluation_grid,
      item$posterior_mean,
      type = "n",
      ylim = ylim,
      xlab = "Time",
      ylab = "Genetic effect",
      main = paste0(
        out$unit_info$cell_id[index],
        "\nindex ",
        index,
        "; lfsr = ",
        formatC(
          per_unit$functional_lfsr[per_unit$index == index],
          digits = 3,
          format = "f"
        )
      )
    )
    polygon(
      c(evaluation_grid, rev(evaluation_grid)),
      c(item$posterior_lower, rev(item$posterior_upper)),
      border = NA,
      col = grDevices::adjustcolor("#4C78A8", alpha.f = 0.18)
    )
    lines(evaluation_grid, item$posterior_mean, col = "#4C78A8", lwd = 2)
    lines(
      evaluation_grid,
      out$true_beta_evaluation[index, ],
      col = "#E45756",
      lwd = 2
    )
    points(
      observed_grid,
      out$eqtl_summary$beta_hat[index, ],
      pch = 19,
      cex = 0.65
    )
    arrows(
      observed_grid,
      out$eqtl_summary$beta_hat[index, ] -
        1.96 * out$eqtl_summary$se[index, ],
      observed_grid,
      out$eqtl_summary$beta_hat[index, ] +
        1.96 * out$eqtl_summary$se[index, ],
      angle = 90,
      code = 3,
      length = 0.03,
      col = "gray45"
    )
    abline(h = c(-0.25, 0, 0.25), lty = c(3, 2, 3), col = "gray55")
  }
  par(old_par)
  dev.off()
}

print(overall_summary)
print(group_summary)
