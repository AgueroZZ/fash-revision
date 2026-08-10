#!/usr/bin/env Rscript

# Summarize a matched five-seed sensitivity analysis for the dynamic-eQTL
# genotype main-effect standard deviation.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) return(".")
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository.")
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

workflowr_root <- find_workflowr_root()

cache_specs <- data.frame(
  dynamic_main_effect_sd = c(0, 0.5, 1, 1.5, 2),
  setting = c("SD = 0", "SD = 0.5", "SD = 1", "SD = 1.5", "SD = 2"),
  output_id = c(
    "sparse_timed_cosine_functional_pilot5",
    "internal_dynamic_main_sd050_pilot5",
    "internal_dynamic_main_n01_pilot5",
    "internal_dynamic_main_sd150_pilot5",
    "internal_dynamic_main_sd200_pilot5"
  ),
  stringsAsFactors = FALSE
)

cache_dir <- function(index) {
  file.path(
    workflowr_root,
    "output",
    "revision_simulations",
    "mc",
    cache_specs$output_id[index]
  )
}

read_switch_curve <- function(index) {
  curve <- read.csv(
    file.path(
      cache_dir(index),
      "summary",
      "functional_testing_mc_alpha_curve.csv"
    ),
    stringsAsFactors = FALSE
  )
  curve <- curve[curve$target == "switch", , drop = FALSE]
  curve$dynamic_main_effect_sd <- cache_specs$dynamic_main_effect_sd[index]
  curve$setting <- cache_specs$setting[index]
  curve
}

read_switch_truth <- function(index) {
  truth <- read.csv(
    file.path(
      cache_dir(index),
      "summary",
      "all_replicate_truth_counts.csv"
    ),
    stringsAsFactors = FALSE
  )
  truth <- truth[truth$target == "switch", , drop = FALSE]
  configuration <- readRDS(file.path(cache_dir(index), "configuration.rds"))
  n_dynamic <- sum(configuration$shape_cell_counts)
  truth$dynamic_main_effect_sd <- cache_specs$dynamic_main_effect_sd[index]
  truth$setting <- cache_specs$setting[index]
  truth$n_dynamic <- n_dynamic
  truth$true_switch_proportion <- truth$n_true_alternatives / n_dynamic
  truth
}

read_pi0 <- function(index) {
  pi0 <- read.csv(
    file.path(
      cache_dir(index),
      "summary",
      "functional_testing_mc_pi0_summary.csv"
    ),
    stringsAsFactors = FALSE
  )
  pi0$dynamic_main_effect_sd <- cache_specs$dynamic_main_effect_sd[index]
  pi0$setting <- cache_specs$setting[index]
  pi0
}

read_truth_transitions <- function(index) {
  full_fit_paths <- list.files(
    file.path(cache_dir(index), "full_fits"),
    pattern = "[.]rds$",
    full.names = TRUE
  )
  if (length(full_fit_paths) == 0) return(NULL)
  do.call(rbind, lapply(full_fit_paths, function(path) {
    out <- readRDS(path)
    dynamic <- out$unit_info$effect_class == "dynamic_bspline"
    data.frame(
      seed = as.integer(sub("seed_([0-9]+)[.]rds", "\\1", basename(path))),
      dynamic_main_effect_sd = cache_specs$dynamic_main_effect_sd[index],
      setting = cache_specs$setting[index],
      deviation_pattern = out$unit_info$sign_pattern[dynamic],
      realized_switch =
        out$true_functionals[out$unit_info$variant_id[dynamic], "switch"] > 0,
      dynamic_main_effect = out$unit_info$dynamic_main_effect[dynamic],
      stringsAsFactors = FALSE
    )
  }))
}

combined_curve <- do.call(
  rbind,
  lapply(seq_len(nrow(cache_specs)), read_switch_curve)
)
truth_counts <- do.call(
  rbind,
  lapply(seq_len(nrow(cache_specs)), read_switch_truth)
)
pi0_summary <- do.call(
  rbind,
  lapply(seq_len(nrow(cache_specs)), read_pi0)
)
transition_rows <- do.call(
  rbind,
  Filter(
    Negate(is.null),
    lapply(seq_len(nrow(cache_specs)), read_truth_transitions)
  )
)
transition_counts <- aggregate(
  list(n_variants = rep(1L, nrow(transition_rows))),
  by = transition_rows[
    ,
    c(
      "dynamic_main_effect_sd",
      "setting",
      "deviation_pattern",
      "realized_switch"
    )
  ],
  FUN = sum
)
transition_counts$proportion_within_pattern <- ave(
  transition_counts$n_variants,
  interaction(
    transition_counts$dynamic_main_effect_sd,
    transition_counts$deviation_pattern
  ),
  FUN = function(x) x / sum(x)
)

alpha_005 <- combined_curve[
  abs(combined_curve$alpha - 0.05) < 1e-12,
  ,
  drop = FALSE
]

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "dynamic_main_effect_sd_ladder"
)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

write_csv(combined_curve, file.path(output_dir, "switch_alpha_curves.csv"))
write_csv(alpha_005, file.path(output_dir, "switch_alpha005_summary.csv"))
write_csv(truth_counts, file.path(output_dir, "true_switch_counts_by_seed.csv"))
write_csv(pi0_summary, file.path(output_dir, "pi0_summary.csv"))
write_csv(transition_rows, file.path(output_dir, "truth_transition_rows.csv"))
write_csv(transition_counts, file.path(output_dir, "truth_transition_counts.csv"))

colors <- setNames(
  c("#000000", "#E69F00", "#009E73", "#0072B2", "#CC79A7"),
  cache_specs$setting
)
metrics <- c(
  mean_power = "Switch power",
  mean_conditional_empirical_fsr = "Conditional empirical FSR",
  mean_empirical_fsr = "End-to-end false-call proportion"
)

png(
  file.path(figure_dir, "bf_sd_ladder_alpha_curves.png"),
  width = 1800,
  height = 620,
  res = 180
)
old_par <- par(
  mfrow = c(1, 3),
  mar = c(4.6, 4.8, 3.8, 1.0),
  oma = c(0, 0, 2.2, 0)
)
for (metric in names(metrics)) {
  plot(
    NA,
    xlim = range(combined_curve$alpha),
    ylim = if (metric == "mean_power") c(0, 1) else c(0, 0.8),
    xlab = "Nominal FSR level alpha",
    ylab = metrics[[metric]],
    main = metrics[[metric]],
    las = 1
  )
  grid(col = "#E6E6E6")
  if (metric != "mean_power") {
    abline(0, 1, lty = 3, lwd = 1.5, col = "#555555")
  }
  for (setting in cache_specs$setting) {
    rows <- combined_curve[
      combined_curve$setting == setting &
        combined_curve$method == "FASH-IWP1-BF",
      ,
      drop = FALSE
    ]
    rows <- rows[order(rows$alpha), ]
    lines(
      rows$alpha,
      rows[[metric]],
      col = colors[[setting]],
      lwd = 2.5
    )
  }
  if (metric == "mean_power") {
    legend(
      "bottomright",
      legend = cache_specs$setting,
      col = colors[cache_specs$setting],
      lwd = 2.5,
      bty = "n",
      cex = 0.82
    )
  }
}
mtext(
  "BF-updated switch testing across dynamic main-effect scales",
  outer = TRUE,
  cex = 1.15,
  font = 2
)
par(old_par)
dev.off()

method_styles <- data.frame(
  method = c("FASH-IWP1-Raw", "FASH-IWP1-BF"),
  label = c("Raw", "BF-updated"),
  color = c("#D55E00", "#0072B2"),
  pch = c(16, 17),
  stringsAsFactors = FALSE
)
png(
  file.path(figure_dir, "alpha005_sd_ladder.png"),
  width = 2100,
  height = 760,
  res = 180
)
old_par <- par(
  mfrow = c(1, 3),
  mar = c(4.8, 5.2, 4.8, 1.2),
  oma = c(0, 0, 3.2, 0)
)
for (metric in names(metrics)) {
  plot(
    NA,
    xlim = range(cache_specs$dynamic_main_effect_sd),
    ylim = if (metric == "mean_power") c(0, 1) else c(0, 0.8),
    xlab = "Dynamic main-effect SD",
    ylab = metrics[[metric]],
    main = metrics[[metric]],
    las = 1
  )
  grid(col = "#E6E6E6")
  if (metric != "mean_power") {
    abline(h = 0.05, lty = 3, lwd = 1.5, col = "#555555")
  }
  for (index in seq_len(nrow(method_styles))) {
    method <- method_styles$method[index]
    rows <- alpha_005[alpha_005$method == method, , drop = FALSE]
    rows <- rows[order(rows$dynamic_main_effect_sd), ]
    lines(
      rows$dynamic_main_effect_sd,
      rows[[metric]],
      col = method_styles$color[index],
      lwd = 2.5
    )
    points(
      rows$dynamic_main_effect_sd,
      rows[[metric]],
      col = method_styles$color[index],
      pch = method_styles$pch[index],
      cex = 1.0
    )
  }
  if (metric == "mean_power") {
    legend(
      "topright",
      legend = method_styles$label,
      col = method_styles$color,
      pch = method_styles$pch,
      lwd = 2.5,
      bty = "n"
    )
  }
}
mtext(
  "Switch testing at alpha = 0.05",
  outer = TRUE,
  line = 1.0,
  cex = 1.15,
  font = 2
)
par(old_par)
dev.off()

truth_mean <- aggregate(
  true_switch_proportion ~ dynamic_main_effect_sd,
  data = truth_counts,
  FUN = mean
)
truth_se <- aggregate(
  true_switch_proportion ~ dynamic_main_effect_sd,
  data = truth_counts,
  FUN = function(x) stats::sd(x) / sqrt(length(x))
)
names(truth_se)[2] <- "mc_se"
truth_mean <- merge(truth_mean, truth_se, by = "dynamic_main_effect_sd")
write_csv(
  truth_mean,
  file.path(output_dir, "true_switch_prevalence_summary.csv")
)

png(
  file.path(figure_dir, "true_switch_prevalence_by_sd.png"),
  width = 900,
  height = 680,
  res = 180
)
plot(
  truth_counts$dynamic_main_effect_sd,
  truth_counts$true_switch_proportion,
  pch = 1,
  col = "#777777",
  xlab = "Dynamic main-effect SD",
  ylab = "True switch proportion among dynamic eQTLs",
  ylim = c(0, max(truth_counts$true_switch_proportion) * 1.08),
  las = 1
)
grid(col = "#E6E6E6")
lines(
  truth_mean$dynamic_main_effect_sd,
  truth_mean$true_switch_proportion,
  lwd = 2.5,
  col = "#0072B2"
)
points(
  truth_mean$dynamic_main_effect_sd,
  truth_mean$true_switch_proportion,
  pch = 16,
  cex = 1.1,
  col = "#0072B2"
)
nonzero_se <- truth_mean$mc_se > 0
arrows(
  truth_mean$dynamic_main_effect_sd[nonzero_se],
  truth_mean$true_switch_proportion[nonzero_se] -
    1.96 * truth_mean$mc_se[nonzero_se],
  truth_mean$dynamic_main_effect_sd[nonzero_se],
  truth_mean$true_switch_proportion[nonzero_se] +
    1.96 * truth_mean$mc_se[nonzero_se],
  angle = 90,
  code = 3,
  length = 0.04,
  col = "#0072B2"
)
dev.off()

print(alpha_005[
  ,
  c(
    "dynamic_main_effect_sd",
    "method",
    "mean_power",
    "mean_conditional_empirical_fsr",
    "mean_empirical_fsr"
  )
])
print(truth_mean)
print(pi0_summary)
