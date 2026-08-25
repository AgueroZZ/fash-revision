#!/usr/bin/env Rscript

# Aggregate seed-level exploratory R3 truth-specification pilot outputs without
# treating individual calls as independent Monte Carlo replicates.

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  position <- match(name, args)
  if (is.na(position) || position == length(args)) return(default)
  args[[position + 1L]]
}

parse_seed_list <- function(x) {
  values <- suppressWarnings(as.integer(strsplit(x, ",", fixed = TRUE)[[1L]]))
  if (length(values) == 0L || anyNA(values) || anyDuplicated(values)) {
    stop("--seed-list must contain unique comma-separated integer seeds.")
  }
  values
}

output_root <- get_arg(
  "--output-root",
  "output/revision_simulations/internal"
)
output_prefix <- get_arg("--output-prefix", "")
seed_list <- parse_seed_list(get_arg(
  "--seed-list", "12345,22345,32345,42345,52345"
))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
summary_file <- get_arg("--summary-file", "")

if (!nzchar(output_prefix) || !is.finite(alpha) || alpha <= 0 || alpha > 1) {
  stop("Invalid aggregation arguments.")
}

result_paths <- file.path(
  output_root,
  paste0(output_prefix, "_seed", seed_list),
  "pilot_result.rds"
)
if (any(!file.exists(result_paths))) {
  stop("At least one requested pilot_result.rds is missing.")
}
results <- lapply(result_paths, readRDS)
configurations <- lapply(results, `[[`, "configuration")
reference <- configurations[[1L]]
comparable_truth_generation_settings <- function(settings) {
  if (!is.list(settings)) return(settings)
  if (!is.null(names(settings))) {
    settings <- settings[!grepl("(^seed$|_seed$)", names(settings))]
  }
  lapply(settings, comparable_truth_generation_settings)
}
if (!all(vapply(configurations, function(x) {
  identical(x$J, reference$J) &&
    identical(x$truth_settings, reference$truth_settings) &&
    identical(
      comparable_truth_generation_settings(x$truth_generation_settings),
      comparable_truth_generation_settings(reference$truth_generation_settings)
    ) &&
    identical(x$middle_window, reference$middle_window) &&
    identical(x$middle_boundary, reference$middle_boundary) &&
    identical(x$middle_definition, reference$middle_definition) &&
    identical(x$package_provenance, reference$package_provenance)
}, logical(1)))) {
  stop("Pilot configurations are not comparable across the requested seeds.")
}
if (!identical(reference$middle_window, c(3, 12)) ||
    !identical(reference$middle_boundary, "open") ||
    !identical(reference$truth_generation_settings$middle_window, c(3, 12)) ||
    !identical(reference$truth_generation_settings$middle_boundary, "open")) {
  stop("Pilot truth generation does not use the required open Middle estimand.")
}

seed_rows <- do.call(rbind, Map(function(result, seed) {
  rows <- result$functional_alpha
  rows <- rows[
    rows$method == "FASH-IWP1-BF" & abs(rows$alpha - alpha) < 1e-12,
    ,
    drop = FALSE
  ]
  rows$seed <- seed
  rows
}, results, seed_list))
if (nrow(seed_rows) == 0L) {
  stop("No BF rows match the requested alpha.")
}

summary <- do.call(rbind, lapply(split(seed_rows, seed_rows$target), function(rows) {
  discoveries <- sum(rows$n_discoveries)
  false_discoveries <- sum(rows$false_discoveries)
  data.frame(
    target = rows$target[[1L]],
    alpha = alpha,
    n_replications = nrow(rows),
    mean_discoveries = mean(rows$n_discoveries),
    mean_power = mean(rows$power),
    mean_empirical_fsr = mean(rows$empirical_fsr),
    empirical_fsr_mc_se = if (nrow(rows) > 1L) {
      stats::sd(rows$empirical_fsr) / sqrt(nrow(rows))
    } else {
      NA_real_
    },
    maximum_seed_empirical_fsr = max(rows$empirical_fsr),
    pooled_discoveries = discoveries,
    pooled_false_discoveries = false_discoveries,
    pooled_empirical_fsr = if (discoveries == 0L) 0 else {
      false_discoveries / discoveries
    },
    stringsAsFactors = FALSE
  )
}))
summary <- summary[match(c("early", "middle", "late", "switch"), summary$target), ]
rownames(summary) <- NULL

if (nzchar(summary_file)) {
  dir.create(dirname(summary_file), recursive = TRUE, showWarnings = FALSE)
  write.csv(summary, summary_file, row.names = FALSE)
}
print(summary, row.names = FALSE)
