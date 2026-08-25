source("code/revision_simulations/appendix_b/appendix_b_helpers.R")

provenance <- require_appendix_b_fashr()
stopifnot(
  identical(provenance$version, APPENDIX_B_FASHR_VERSION),
  identical(provenance$remote_sha, APPENDIX_B_FASHR_REMOTE_SHA)
)

rho_values <- seq(0.05, 0.50, by = 0.01)
for (rho_dynamic in rho_values) {
  counts <- appendix_b_class_counts(
    J = 1000L,
    rho_dynamic = rho_dynamic,
    rho_nonlinear = rho_dynamic / 2
  )
  stopifnot(
    sum(counts) == 1000L,
    counts[["nondynamic"]] == round(1000 * (1 - rho_dynamic)),
    counts[["linear"]] == round(1000 * rho_dynamic / 2),
    counts[["nonlinear"]] == round(1000 * rho_dynamic / 2)
  )
}

set.seed(12345)
bundle_a <- build_appendix_b_datasets(
  J = 20L,
  rho_dynamic = 0.4,
  rho_nonlinear = 0.2,
  sigma_vec = c(0.1, 0.3, 0.5)
)
set.seed(12345)
bundle_b <- build_appendix_b_datasets(
  J = 20L,
  rho_dynamic = 0.4,
  rho_nonlinear = 0.2,
  sigma_vec = c(0.1, 0.3, 0.5)
)
class_table <- table(bundle_a$truth$class)
stopifnot(
  identical(bundle_a, bundle_b),
  identical(unname(bundle_a$class_counts), c(12L, 4L, 4L)),
  identical(names(bundle_a$data_list), bundle_a$truth$unit_id),
  identical(as.integer(class_table), c(4L, 12L, 4L)),
  identical(names(class_table), c(
    "linear",
    "nondynamic",
    "nonlinear"
  )),
  identical(stats::setNames(as.integer(class_table), names(class_table)), c(
    linear = 4L,
    nondynamic = 12L,
    nonlinear = 4L
  ))
)

lfdr <- c(0.01, 0.02, 0.20, 0.90)
calls <- cumulative_lfdr_calls(lfdr, alpha = 0.05)
stopifnot(
  identical(calls$indices, c(1L, 2L)),
  calls$n_calls == 2L
)

metrics <- evaluate_appendix_b_discoveries(
  lfdr = lfdr,
  truth_positive = c(TRUE, FALSE, TRUE, FALSE),
  alpha = 0.05
)
stopifnot(
  metrics$discoveries == 2L,
  metrics$true_discoveries == 1L,
  metrics$false_discoveries == 1L,
  metrics$realized_fdp == 0.5,
  metrics$power == 0.5
)

zero_metrics <- evaluate_appendix_b_discoveries(
  lfdr = c(0.8, 0.9),
  truth_positive = c(FALSE, TRUE),
  alpha = 0.05
)
stopifnot(
  zero_metrics$discoveries == 0L,
  zero_metrics$realized_fdp == 0,
  zero_metrics$power == 0
)

fit_with_omitted_zero_weight <- list(
  psd_grid = c(0, 0.1, 0.2),
  prior_weights = data.frame(
    psd = c(0.1, 0.2),
    prior_weight = c(0.4, 0.6)
  )
)
stopifnot(extract_fash_null_weight(fit_with_omitted_zero_weight) == 0)

fit_with_explicit_zero_weight <- fit_with_omitted_zero_weight
fit_with_explicit_zero_weight$prior_weights <- data.frame(
  psd = c(0, 0.1, 0.2),
  prior_weight = c(0, 0.4, 0.6)
)
stopifnot(extract_fash_null_weight(fit_with_explicit_zero_weight) == 0)

invalid_counts <- tryCatch(
  {
    appendix_b_class_counts(101L, 0.33, 0.165)
    NULL
  },
  error = function(condition) condition
)
stopifnot(inherits(invalid_counts, "error"))

cat("Appendix B helper tests passed.\n")
