#!/usr/bin/env Rscript

# Quantify how the Middle definition changes R3A proposal and accepted truths.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.", call. = FALSE)
}

row_max <- function(x) {
  apply(x, 1, max)
}

summarize_numeric <- function(x, prefix) {
  values <- c(
    mean = mean(x),
    sd = stats::sd(x),
    q05 = unname(stats::quantile(x, 0.05)),
    q10 = unname(stats::quantile(x, 0.10)),
    q25 = unname(stats::quantile(x, 0.25)),
    median = stats::median(x),
    q75 = unname(stats::quantile(x, 0.75)),
    q90 = unname(stats::quantile(x, 0.90)),
    q95 = unname(stats::quantile(x, 0.95)),
    max = max(x)
  )
  stats::setNames(as.list(values), paste(prefix, names(values), sep = "_"))
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
output_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "diagnostics", "r3_middle_calibration"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

time_grid <- 0:15
evaluation_grid <- seq(0, 15, by = 0.1)
time_scaled <- time_grid / 15
evaluation_scaled <- evaluation_grid / 15
basis_observed <- splines::bs(
  time_scaled,
  df = 6,
  degree = 3,
  intercept = TRUE
)
basis_evaluation <- predict(basis_observed, newx = evaluation_scaled)

closed_regions <- list(
  early = evaluation_grid <= 3,
  middle = evaluation_grid >= 4 & evaluation_grid <= 11,
  late = evaluation_grid >= 12,
  unassigned = (evaluation_grid > 3 & evaluation_grid < 4) |
    (evaluation_grid > 11 & evaluation_grid < 12)
)
open_regions <- list(
  early = evaluation_grid <= 3,
  middle = evaluation_grid > 3 & evaluation_grid < 12,
  late = evaluation_grid >= 12
)

classify_region <- function(abs_curves, regions) {
  maxima <- vapply(
    regions,
    function(region) row_max(abs_curves[, region, drop = FALSE]),
    numeric(nrow(abs_curves))
  )
  if (is.null(dim(maxima))) {
    maxima <- matrix(maxima, nrow = nrow(abs_curves))
  }
  colnames(maxima) <- names(regions)
  category <- colnames(maxima)[max.col(maxima, ties.method = "first")]
  ordered <- t(apply(maxima, 1, sort, decreasing = TRUE))
  margin <- ordered[, 1] - ordered[, 2]
  list(category = category, margin = margin, maxima = maxima)
}

set.seed(20260823L)
n_proposals <- 100000L
chunk_size <- 2500L
proposal_rows <- vector("list", ceiling(n_proposals / chunk_size))
proposal_index <- 1L
completed <- 0L
while (completed < n_proposals) {
  n_chunk <- min(chunk_size, n_proposals - completed)
  main_effect <- stats::rnorm(n_chunk, mean = 0, sd = 1)
  coefficients <- matrix(
    stats::rnorm(n_chunk * ncol(basis_observed), mean = 0, sd = 1),
    nrow = n_chunk,
    byrow = TRUE
  )
  raw_observed <- coefficients %*% t(basis_observed)
  raw_evaluation <- coefficients %*% t(basis_evaluation)
  observed_mean <- rowMeans(raw_observed)
  deviation_observed <- raw_observed - observed_mean
  deviation_evaluation <- raw_evaluation - observed_mean
  scale_factor <- 2 / row_max(abs(deviation_observed))
  beta_evaluation <- deviation_evaluation * scale_factor + main_effect
  abs_curves <- abs(beta_evaluation)
  closed <- classify_region(abs_curves, closed_regions)
  open <- classify_region(abs_curves, open_regions)
  proposal_rows[[proposal_index]] <- data.frame(
    proposal_index = completed + seq_len(n_chunk),
    closed_category = closed$category,
    open_category = open$category,
    closed_margin = closed$margin,
    open_margin = open$margin,
    stringsAsFactors = FALSE
  )
  completed <- completed + n_chunk
  proposal_index <- proposal_index + 1L
}
proposals <- do.call(rbind, proposal_rows)

transition <- as.data.frame(table(
  factor(
    proposals$closed_category,
    levels = c("early", "middle", "late", "unassigned")
  ),
  factor(proposals$open_category, levels = c("early", "middle", "late"))
))
names(transition) <- c("closed_category", "open_category", "count")
transition$proportion <- transition$count / n_proposals
utils::write.csv(
  transition,
  file.path(output_dir, "proposal_transition_matrix.csv"),
  row.names = FALSE
)

category_probabilities <- do.call(rbind, lapply(
  c("closed", "open"),
  function(definition) {
    category <- proposals[[paste0(definition, "_category")]]
    margin <- proposals[[paste0(definition, "_margin")]]
    levels <- if (definition == "closed") {
      c("early", "middle", "late", "unassigned")
    } else {
      c("early", "middle", "late")
    }
    do.call(rbind, lapply(levels, function(level) {
      indicator <- category == level
      probability <- mean(indicator)
      data.frame(
        definition = definition,
        category = level,
        proposals = n_proposals,
        count = sum(indicator),
        probability = probability,
        mc_se = sqrt(probability * (1 - probability) / n_proposals),
        probability_with_margin_0p10 = mean(indicator & margin >= 0.10),
        probability_with_margin_0p25 = mean(indicator & margin >= 0.25),
        stringsAsFactors = FALSE
      )
    }))
  }
))
utils::write.csv(
  category_probabilities,
  file.path(output_dir, "proposal_category_probabilities.csv"),
  row.names = FALSE
)

seed_list <- c(12345L, 22345L, 32345L, 42345L, 52345L)
definitions <- list(
  closed = list(window = c(4, 11), boundary = "closed"),
  open = list(window = c(3, 12), boundary = "open")
)
accepted_units <- list()
accepted_index <- 1L
for (definition_name in names(definitions)) {
  definition <- definitions[[definition_name]]
  for (seed in seed_list) {
    effect_sim <- simulate_matched_functional_effect_set(
      n_variants = 6362L,
      truth_mechanism = "random_bspline",
      evaluation_grid = evaluation_grid,
      class_probs = c(
        dynamic_bspline = 0.20,
        constant = 0.40,
        zero = 0.40
      ),
      dynamic_main_effect_sd = 1,
      location_truth_margin = 0.10,
      switch_truth_margin = 0.10,
      non_switch_min_abs = 0.10,
      non_switch_min_range_fraction = 0.10,
      seed = seed,
      middle_window = definition$window,
      middle_boundary = definition$boundary
    )
    dynamic_index <- which(
      effect_sim$unit_info$effect_class == "dynamic_bspline"
    )
    info <- effect_sim$unit_info[dynamic_index, , drop = FALSE]
    functionals <- effect_sim$true_functionals[dynamic_index, , drop = FALSE]
    target_columns <- match(info$time_group, colnames(functionals))
    target_functional <- functionals[cbind(seq_len(nrow(info)), target_columns)]
    nearest_competing_functional <- vapply(
      seq_len(nrow(info)),
      function(position) {
        max(functionals[
          position,
          setdiff(c("early", "middle", "late"), info$time_group[[position]])
        ])
      },
      numeric(1)
    )
    accepted_units[[accepted_index]] <- data.frame(
      definition = definition_name,
      seed = seed,
      truth_group = info$truth_group,
      time_group = info$time_group,
      switch_status = info$switch_status,
      generation_attempt = info$generation_attempt,
      target_functional = target_functional,
      nearest_competing_functional = nearest_competing_functional,
      centered_rms = info$centered_rms,
      effect_range = info$effect_range,
      absolute_main_effect = abs(info$genetic_main_effect),
      stringsAsFactors = FALSE
    )
    accepted_index <- accepted_index + 1L
  }
}
accepted_units <- do.call(rbind, accepted_units)
saveRDS(
  accepted_units,
  file.path(output_dir, "accepted_truth_unit_metrics.rds"),
  version = 3
)

group_keys <- unique(accepted_units[, c(
  "definition", "seed", "truth_group", "time_group", "switch_status"
)])
accepted_summary <- do.call(rbind, lapply(
  seq_len(nrow(group_keys)),
  function(key_index) {
    key <- group_keys[key_index, , drop = FALSE]
    rows <- accepted_units[
      accepted_units$definition == key$definition &
        accepted_units$seed == key$seed &
        accepted_units$truth_group == key$truth_group,
      ,
      drop = FALSE
    ]
    values <- c(
      as.list(key),
      list(n = nrow(rows)),
      summarize_numeric(rows$generation_attempt, "attempt"),
      summarize_numeric(rows$target_functional, "target"),
      summarize_numeric(rows$nearest_competing_functional, "competitor"),
      summarize_numeric(rows$centered_rms, "centered_rms"),
      summarize_numeric(rows$effect_range, "effect_range"),
      summarize_numeric(rows$absolute_main_effect, "absolute_main_effect")
    )
    as.data.frame(values, stringsAsFactors = FALSE)
  }
))
rownames(accepted_summary) <- NULL
utils::write.csv(
  accepted_summary,
  file.path(output_dir, "accepted_truth_summary.csv"),
  row.names = FALSE
)

cat("Wrote paired proposal and accepted-truth geometry diagnostics to:\n")
cat(normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
