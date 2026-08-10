#!/usr/bin/env Rscript

# Aggregate the locked five-seed confirmation of the localized-pair spiky truth
# mechanism. This script writes internal diagnostics only.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) return(".")
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find code/revision_simulations/shared/simulation_functions.R.")
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

summarize_spiky_power <- function(rows, confidence_level = 0.95) {
  groups <- split(rows, list(rows$method, rows$alpha), drop = TRUE)
  out <- lapply(groups, function(x) {
    power_summary <- summarize_mc_values(x$power, confidence_level)
    data.frame(
      shape_profile = "spiky",
      method = x$method[1],
      alpha = x$alpha[1],
      n_dynamic = x$n_dynamic[1],
      n_replications = length(unique(x$seed)),
      mean_true_positives = mean(x$true_positives),
      mean_power = power_summary[["mean"]],
      power_ci_lower = pmax(0, power_summary[["lower"]]),
      power_ci_upper = pmin(1, power_summary[["upper"]]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

seeds <- c(12345L, 22345L, 32345L, 42345L, 52345L)
input_paths <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  paste0(
    "confirm_spiky_localpair_df12_frac055_070_seed",
    seeds,
    "_B100"
  ),
  "candidates",
  "spiky_df12.rds"
)
if (any(!file.exists(input_paths))) {
  stop("One or more locked spiky confirmation caches are missing.")
}
replicates <- lapply(input_paths, readRDS)
if (any(vapply(replicates, function(x) x$spiky_df != 12L, logical(1)))) {
  stop("A confirmation cache has the wrong spiky df.")
}

output_dir <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "confirm_spiky_localpair_df12_frac055_070_B100_5seed"
)
figure_dir <- file.path(output_dir, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
alpha_grid <- seq(0.005, 0.20, by = 0.005)

all_alpha <- do.call(rbind, lapply(seq_along(replicates), function(i) {
  out <- replicates[[i]]$alpha_curve
  out$seed <- seeds[i]
  out
}))
mc_alpha <- summarize_mc_alpha_curves(all_alpha)
mc_alpha_005 <- mc_alpha[abs(mc_alpha$alpha - 0.05) < 1e-12, ]

all_spiky_alpha <- do.call(rbind, lapply(seq_along(replicates), function(i) {
  replicate <- replicates[[i]]
  unit_info <- replicate$effect_sim$unit_info
  spiky_units <- unit_info$unit_index[
    unit_info$effect_class == "dynamic_bspline" &
      unit_info$shape_profile == "spiky"
  ]
  result_rows <- lapply(methods, function(method) {
    method_result <- replicate$result_table[
      replicate$result_table$method == method,
    ]
    do.call(rbind, lapply(alpha_grid, function(alpha) {
      selected_units <- method_result$unit_index[
        is.finite(method_result$adjusted_score) &
          method_result$adjusted_score <= alpha
      ]
      true_positives <- sum(spiky_units %in% selected_units)
      data.frame(
        seed = seeds[i],
        method = method,
        alpha = alpha,
        n_dynamic = length(spiky_units),
        true_positives = true_positives,
        power = true_positives / length(spiky_units),
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, result_rows)
}))
mc_spiky_alpha <- summarize_spiky_power(all_spiky_alpha)
mc_spiky_alpha_005 <- mc_spiky_alpha[
  abs(mc_spiky_alpha$alpha - 0.05) < 1e-12,
]

all_units <- do.call(rbind, lapply(seq_along(replicates), function(i) {
  out <- replicates[[i]]$unit_diagnostics
  out$seed <- seeds[i]
  out
}))
spiky_units <- all_units[
  all_units$effect_class == "dynamic_bspline" &
    all_units$shape_profile == "spiky",
]
pattern_summary <- aggregate(
  cbind(
    raw_selected,
    bf_selected,
    direct_quadratic_selected,
    log_bayes_factor,
    support_25,
    support_50,
    second_difference_roughness,
    quadratic_projection
  ) ~ spike_pattern,
  data = spiky_units,
  FUN = mean
)
names(pattern_summary)[names(pattern_summary) == "raw_selected"] <- "raw_power"
names(pattern_summary)[names(pattern_summary) == "bf_selected"] <- "bf_power"
names(pattern_summary)[
  names(pattern_summary) == "direct_quadratic_selected"
] <- "direct_quadratic_power"

per_seed <- do.call(rbind, lapply(seq_along(replicates), function(i) {
  out <- replicates[[i]]$unit_diagnostics
  spiky <- which(
    out$effect_class == "dynamic_bspline" &
      out$shape_profile == "spiky"
  )
  data.frame(
    seed = seeds[i],
    raw_spiky_power = mean(out$raw_selected[spiky]),
    bf_spiky_power = mean(out$bf_selected[spiky]),
    direct_quadratic_spiky_power =
      mean(out$direct_quadratic_selected[spiky]),
    bf_minus_direct_quadratic =
      mean(out$bf_selected[spiky]) -
      mean(out$direct_quadratic_selected[spiky]),
    median_log_bayes_factor = median(out$log_bayes_factor[spiky]),
    mean_support_25 = mean(out$support_25[spiky]),
    mean_second_difference_roughness =
      mean(out$second_difference_roughness[spiky]),
    mean_quadratic_projection = mean(out$quadratic_projection[spiky]),
    raw_pi0 = replicates[[i]]$summary$raw_pi0,
    bf_pi0 = replicates[[i]]$summary$bf_pi0,
    stringsAsFactors = FALSE
  )
}))

old_shape_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "confirm_broad075_df10_B100_5seed",
  "mc_shape_power_alpha005.csv"
)
old_spiky <- read.csv(old_shape_path)
old_spiky <- old_spiky[
  old_spiky$shape_profile == "spiky" &
    old_spiky$method %in% methods,
  c("method", "mean_power")
]
old_spiky$mechanism <- "previous multi-spike"
new_spiky <- mc_spiky_alpha_005[, c("method", "mean_power")]
new_spiky$mechanism <- "localized-pair multi-spike"
mechanism_comparison <- rbind(old_spiky, new_spiky)
mechanism_comparison <- mechanism_comparison[
  order(mechanism_comparison$method, mechanism_comparison$mechanism),
]

write_csv(all_alpha, file.path(output_dir, "all_replicate_global_alpha.csv"))
write_csv(mc_alpha, file.path(output_dir, "mc_global_alpha_curve.csv"))
write_csv(mc_alpha_005, file.path(output_dir, "mc_global_alpha005.csv"))
write_csv(all_spiky_alpha, file.path(output_dir, "all_replicate_spiky_alpha.csv"))
write_csv(mc_spiky_alpha, file.path(output_dir, "mc_spiky_alpha_curve.csv"))
write_csv(mc_spiky_alpha_005, file.path(output_dir, "mc_spiky_alpha005.csv"))
write_csv(all_units, file.path(output_dir, "all_replicate_unit_diagnostics.csv"))
write_csv(pattern_summary, file.path(output_dir, "spike_pattern_summary_alpha005.csv"))
write_csv(per_seed, file.path(output_dir, "per_seed_spiky_summary_alpha005.csv"))
write_csv(mechanism_comparison, file.path(output_dir, "old_vs_new_spiky_power_alpha005.csv"))

subtitle <- paste0(
  "5 seeds; local cubic B-spline df=12; nearest secondary fraction 0.55-0.70; ",
  "direct eFDR uses true pi0 and B=100"
)
plot_mc_shape_power_grid(
  mc_curve = mc_spiky_alpha,
  shape_order = "spiky",
  file = file.path(figure_dir, "spiky_power_across_alpha.png"),
  title = "Internal confirmation: localized-pair spiky effects",
  subtitle = subtitle,
  style_profile = "combined",
  width = 1500,
  height = 1050
)
plot_mc_alpha_curves(
  mc_curve = mc_alpha,
  metric = "fdr",
  file = file.path(figure_dir, "global_fdr_across_alpha.png"),
  title = "Internal confirmation: global empirical FDR",
  subtitle = subtitle,
  style_profile = "combined"
)

example <- replicates[[1]]$effect_sim
example_info <- example$unit_info
patterns <- c("single", "same-sign double", "opposite-sign double")
example_index <- unlist(lapply(patterns, function(pattern) {
  candidates <- which(
    example_info$effect_class == "dynamic_bspline" &
      example_info$shape_profile == "spiky" &
      example_info$spike_pattern == pattern
  )
  head(candidates, 4)
}))
png(
  file.path(figure_dir, "truth_examples.png"),
  width = 1800,
  height = 1350,
  res = 180
)
par(mfrow = c(3, 4), mar = c(3.1, 3.4, 2.4, 0.8), oma = c(0, 0, 2, 0))
for (j in example_index) {
  plot(
    0:15,
    example$beta_matrix[j, ],
    type = "l",
    lwd = 2,
    col = "#D95F02",
    xlab = "Time",
    ylab = "True effect",
    main = paste(example_info$time_group[j], example_info$spike_pattern[j]),
    cex.main = 0.78
  )
  abline(h = 0, col = "grey70", lty = 3)
}
mtext(
  "Localized-pair multi-spike truth examples",
  outer = TRUE,
  cex = 1.1
)
dev.off()

saveRDS(
  list(
    seeds = seeds,
    all_alpha = all_alpha,
    mc_alpha = mc_alpha,
    all_spiky_alpha = all_spiky_alpha,
    mc_spiky_alpha = mc_spiky_alpha,
    all_units = all_units,
    pattern_summary = pattern_summary,
    per_seed = per_seed,
    mechanism_comparison = mechanism_comparison
  ),
  file.path(output_dir, "confirmation_summary.rds")
)

message("Saved localized-pair spiky confirmation summary to: ", output_dir)
