#!/usr/bin/env Rscript

workflowr_root <- if (file.exists("_workflowr.yml")) {
  "."
} else if (file.exists("coderepo-local/_workflowr.yml")) {
  "coderepo-local"
} else {
  stop("Run this test from the workflowr root or its parent.")
}

source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "matched_fash_linear_real_data", "scoped_shared_provenance.R"
))

cache <- readRDS(file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "matched_fash_linear_real_data_fashr_0_1_43",
  "analysis_cache.rds"
))
shared_index <- match("shared_functions", cache$input_provenance$label)
stopifnot(!is.na(shared_index))

historical_path <- file.path(
  workflowr_root,
  "code", "revision_simulations", "r3_r4_fashr0143", "source_snapshots",
  "r4_simulation_functions.R"
)
current_path <- cache$input_provenance$path[[shared_index]]
recorded_md5 <- cache$input_provenance$md5[[shared_index]]
result <- validate_matched_r7_scoped_shared_provenance(
  historical_path = historical_path,
  current_path = current_path,
  recorded_md5 = recorded_md5
)

expected_dependency_functions <- c(
  "BF_update_linear_mixture_fash",
  "compute_linear_mixture_log_likelihood",
  "default_revision_grid",
  "expand_grid_prior_weights",
  "fit_linear_mixture_fash_from_log_likelihood",
  "log_marginal_common_intercept",
  "solve_from_chol",
  "validate_linear_mixture_fash",
  "validate_linear_mixture_grid"
)
stopifnot(
  isTRUE(result$passed),
  identical(
    result$historical_file_md5,
    "d9748fb9cffa5cfa0fa0fa1f917a6954"
  ),
  setequal(result$dependency_functions, expected_dependency_functions),
  nrow(result$comparison) == length(expected_dependency_functions),
  all(result$comparison$formals_identical),
  all(result$comparison$body_identical)
)

wrong_md5 <- try(validate_matched_r7_scoped_shared_provenance(
  historical_path = historical_path,
  current_path = current_path,
  recorded_md5 = paste(rep("0", 32L), collapse = "")
), silent = TRUE)
stopifnot(inherits(wrong_md5, "try-error"))

cat("R7 scoped shared-function provenance tests passed.\n")
