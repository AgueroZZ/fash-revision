#!/usr/bin/env Rscript

# Diagnose how the fitted dynamic FASH prior allocates maxima across the old
# and revised Early/Middle/Late regions. The null component is excluded because
# constant curves create structural ties and the R3 functional test is applied
# after the dynamic screen.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(fashr)
})

project_root <- normalizePath(getwd(), mustWork = TRUE)
output_dir <- file.path(
  project_root,
  "output/revision_simulations/diagnostics/r3_middle_calibration"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

time_grid <- seq(0, 15, by = 0.1)
draws_per_fit <- 50000L

fit_sets <- c(
  w150_minabs050 = file.path(
    project_root,
    "output/revision_simulations/mc",
    "internal_functional_smoothness_w150_minabs050_pilot5/full_fits"
  ),
  w300_minabs0925 = file.path(
    project_root,
    "output/revision_simulations/mc",
    "internal_functional_smoothness_w300_minabs0925_pilot5/full_fits"
  )
)

definitions <- list(
  closed = list(
    early = time_grid <= 3,
    middle = time_grid >= 4 & time_grid <= 11,
    late = time_grid >= 12
  ),
  open = list(
    early = time_grid <= 3,
    middle = time_grid > 3 & time_grid < 12,
    late = time_grid >= 12
  )
)

classify_maxima <- function(draws, regions) {
  max_index <- max.col(t(abs(draws)), ties.method = "first")
  labels <- rep("unassigned", ncol(draws))
  for (region_name in names(regions)) {
    labels[regions[[region_name]][max_index]] <- region_name
  }
  labels
}

conditional_dynamic_fit <- function(fit) {
  keep <- fit$prior_weights$psd > 0 & fit$prior_weights$prior_weight > 0
  if (!any(keep)) {
    stop("The fitted prior has no positive-PSD component.")
  }
  fit$prior_weights <- fit$prior_weights[keep, , drop = FALSE]
  fit$prior_weights$prior_weight <-
    fit$prior_weights$prior_weight / sum(fit$prior_weights$prior_weight)
  fit
}

summarize_labels <- function(labels, definition, fit_set, seed, method,
                             component = "all") {
  categories <- if (definition == "closed") {
    c("early", "middle", "late", "unassigned")
  } else {
    c("early", "middle", "late")
  }
  counts <- table(factor(labels, levels = categories))
  n <- length(labels)
  data.frame(
    fit_set = fit_set,
    seed = seed,
    method = method,
    definition = definition,
    component = as.character(component),
    category = categories,
    draws = n,
    count = as.integer(counts),
    probability = as.integer(counts) / n,
    mc_se = sqrt((as.integer(counts) / n) *
                   (1 - as.integer(counts) / n) / n),
    stringsAsFactors = FALSE
  )
}

result_rows <- list()
weight_rows <- list()
result_index <- 0L
weight_index <- 0L

for (fit_set_name in names(fit_sets)) {
  files <- sort(list.files(
    fit_sets[[fit_set_name]],
    pattern = "^seed_[0-9]+[.]rds$",
    full.names = TRUE
  ))
  if (length(files) == 0L) {
    stop("No retained full fits found in: ", fit_sets[[fit_set_name]])
  }

  for (file in files) {
    retained <- readRDS(file)
    seed <- as.integer(sub("^seed_([0-9]+)[.]rds$", "\\1", basename(file)))

    for (method in c("raw", "bf")) {
      fit_name <- paste0("fash_iwp1_", method)
      fit <- conditional_dynamic_fit(retained$fash_fits[[fit_name]])

      weight_index <- weight_index + 1L
      weight_rows[[weight_index]] <- data.frame(
        fit_set = fit_set_name,
        seed = seed,
        method = method,
        psd = fit$prior_weights$psd,
        conditional_dynamic_weight = fit$prior_weights$prior_weight,
        stringsAsFactors = FALSE
      )

      method_offset <- if (method == "raw") 0L else 1000000L
      fit_set_offset <- match(fit_set_name, names(fit_sets)) * 10000000L
      set.seed(seed + method_offset + fit_set_offset)
      simulated <- fashr::simulate_fash_prior(
        fit,
        M = draws_per_fit,
        x_range = c(0, 15),
        x_new = time_grid
      )

      for (definition in names(definitions)) {
        labels <- classify_maxima(simulated$samples, definitions[[definition]])

        result_index <- result_index + 1L
        result_rows[[result_index]] <- summarize_labels(
          labels,
          definition,
          fit_set_name,
          seed,
          method
        )

        for (component_index in sort(unique(simulated$component))) {
          component_keep <- simulated$component == component_index
          result_index <- result_index + 1L
          result_rows[[result_index]] <- summarize_labels(
            labels[component_keep],
            definition,
            fit_set_name,
            seed,
            method,
            component = sprintf(
              "psd=%.8g",
              fit$prior_weights$psd[component_index]
            )
          )
        }
      }
    }
  }
}

results <- do.call(rbind, result_rows)
weights <- do.call(rbind, weight_rows)

write.csv(
  results,
  file.path(output_dir, "fitted_dynamic_prior_category_probabilities.csv"),
  row.names = FALSE
)
write.csv(
  weights,
  file.path(output_dir, "fitted_dynamic_prior_component_weights.csv"),
  row.names = FALSE
)

overall <- results[results$component == "all", , drop = FALSE]
aggregate_summary <- aggregate(
  overall[c("probability")],
  overall[c("fit_set", "method", "definition", "category")],
  function(x) c(mean = mean(x), sd = sd(x), min = min(x), max = max(x))
)
probability_summary <- as.data.frame(aggregate_summary$probability)
names(probability_summary) <- c("mean", "sd", "min", "max")
aggregate_summary <- data.frame(
  aggregate_summary[c("fit_set", "method", "definition", "category")],
  probability_summary,
  row.names = NULL
)
write.csv(
  aggregate_summary,
  file.path(output_dir, "fitted_dynamic_prior_category_summary.csv"),
  row.names = FALSE
)

cat("Wrote fitted dynamic-prior geometry diagnostics to:\n")
cat(output_dir, "\n")
