#!/usr/bin/env Rscript

# Validate the two-part reporting contract for the internal R1 page.

source(file.path(
  "code", "revision_simulations", "internal",
  "r1_normality_assumption", "reporting.R"
))

# Recompute every cached QQ quantile from the auditable residual-level cache.
recomputed_qq <- do.call(rbind, lapply(c(1L, 3L, 9L), function(day) {
  cached_day_qq <- residual_qq[residual_qq$day == day, , drop = FALSE]
  day_residuals <- standardized_residuals$standardized_residual[
    standardized_residuals$day == day
  ]
  data.frame(
    day = day,
    probability = cached_day_qq$probability,
    empirical_quantile = as.numeric(stats::quantile(
      day_residuals,
      probs = cached_day_qq$probability,
      names = FALSE,
      type = 8
    ))
  )
}))
residual_qq_order <- order(residual_qq$day, residual_qq$probability)
recomputed_qq_order <- order(
  recomputed_qq$day,
  recomputed_qq$probability
)

stopifnot(
  nrow(selected_units) == 6362L,
  !anyDuplicated(selected_units$gene_id),
  !anyDuplicated(selected_units$unit_key),
  nrow(standardized_residuals) == 362634L,
  setequal(unique(standardized_residuals$day), c(1L, 3L, 9L)),
  all(table(standardized_residuals$day) == 120878L),
  all(standardized_residuals$residual_df == 12L),
  max(abs(
    residual_qq$empirical_quantile[residual_qq_order] -
      recomputed_qq$empirical_quantile[recomputed_qq_order]
  )) < 1e-12,
  max(abs(
    residual_qq$normal_reference_quantile -
      stats::qnorm(residual_qq$probability)
  )) < 1e-12,
  all(residual_day_summary$n_valid_fits == 6362L),
  all(residual_day_summary$n_excluded_fits == 0L),
  max(abs(residual_day_summary$sd - c(
    1.00178837050339,
    1.01599796043499,
    1.00119090038572
  ))) < 1e-12,
  all(residual_validation$pass),
  all(simulation_validation$pass)
)

# Enforce self-explanatory figure/table code and traditional terminology.
rmd_text <- paste(readLines(
  file.path("analysis", "revision_internal_r1_normality_assumption.rmd"),
  warn = FALSE
), collapse = "\n")
stopifnot(
  grepl("residual_plot_data", rmd_text, fixed = TRUE),
  grepl("standardized_residual", rmd_text, fixed = TRUE),
  grepl("stats::qnorm", rmd_text, fixed = TRUE),
  grepl("stats::quantile", rmd_text, fixed = TRUE),
  grepl("geom_histogram", rmd_text, fixed = TRUE),
  grepl("stats::dnorm", rmd_text, fixed = TRUE),
  grepl("normal_reference_density", rmd_text, fixed = TRUE),
  grepl("Day 1", rmd_text, fixed = TRUE),
  grepl("Day 3", rmd_text, fixed = TRUE),
  grepl("Day 9", rmd_text, fixed = TRUE),
  grepl("Figure 1.", rmd_text, fixed = TRUE),
  grepl("Figure 2.", rmd_text, fixed = TRUE),
  !grepl("externally studentized", rmd_text, ignore.case = TRUE),
  !grepl("t\\(11\\)", rmd_text),
  !grepl("plot_zero_null_qq\\(\\)", rmd_text),
  !grepl("plot_realized_fdp_curves\\(\\)", rmd_text)
)

# Protect the accepted Part 2 result from accidental changes.
stopifnot(
  nrow(gaussian_performance_alpha005) == 4L,
  setequal(
    gaussian_performance_alpha005$se_scale,
    c("raw regression SE", "t-adjusted SE")
  ),
  setequal(gaussian_performance_alpha005$fit_stage, c("Raw", "BF")),
  all(abs(gaussian_performance_alpha005$alpha - 0.05) < 1e-12)
)
performance_key <- paste(
  gaussian_performance_alpha005$se_scale,
  gaussian_performance_alpha005$fit_stage,
  sep = " | "
)
observed_fdp <- stats::setNames(
  gaussian_performance_alpha005$realized_fdp,
  performance_key
)
expected_fdp <- c(
  "raw regression SE | Raw" = 0.340063761955366,
  "raw regression SE | BF" = 0.146666666666667,
  "t-adjusted SE | Raw" = 0.0576923076923077,
  "t-adjusted SE | BF" = 0.0332778702163062
)
stopifnot(
  setequal(names(observed_fdp), names(expected_fdp)),
  max(abs(observed_fdp[names(expected_fdp)] - expected_fdp)) < 1e-12
)

cat("R1 residual-normality and SE-consequence reporting tests passed.\n")
