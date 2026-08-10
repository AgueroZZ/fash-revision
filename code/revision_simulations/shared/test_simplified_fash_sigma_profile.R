#!/usr/bin/env Rscript

# Focused tests for the single-Gaussian slope-scale profile used by the
# reviewer-facing R1 and R2 simulations.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the shared revision simulation functions.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "shared",
  "simulation_functions.R"
))

set.seed(20260807)
n_units <- 120L
n_time <- 8L
time_grid <- seq_len(n_time) - 1
scaled_time <- time_grid / max(time_grid)
is_dynamic <- seq_len(n_units) > 80L
true_slope <- numeric(n_units)
true_slope[is_dynamic] <- stats::rnorm(sum(is_dynamic), sd = 1.5)
intercept <- stats::rnorm(n_units, sd = 0.75)
standard_error <- 0.15

datasets <- lapply(seq_len(n_units), function(j) {
  list(
    x = time_grid,
    y = intercept[j] + true_slope[j] * scaled_time +
      stats::rnorm(n_time, sd = standard_error),
    sd = rep(standard_error, n_time)
  )
})
names(datasets) <- sprintf("unit_%03d", seq_len(n_units))

sigma_grid <- c(0.25, 0.5, 1, 2, 4)
fit <- fit_simplified_fash(
  datasets = datasets,
  estimate_sigma = TRUE,
  sigma_beta_grid = sigma_grid,
  scale_time = TRUE
)

stopifnot(
  nrow(fit$sigma_profile) == length(sigma_grid),
  identical(fit$sigma_profile$sigma_beta, sigma_grid),
  all(is.finite(fit$sigma_profile$estimated_pi0)),
  all(fit$sigma_profile$estimated_pi0 >= 0),
  all(fit$sigma_profile$estimated_pi0 <= 1),
  all(is.finite(fit$sigma_profile$loglik)),
  sum(fit$sigma_profile$selected) == 1L,
  !fit$selected_sigma_on_boundary,
  isTRUE(all.equal(
    fit$sigma_beta,
    fit$sigma_profile$sigma_beta[fit$sigma_profile$selected]
  )),
  isTRUE(all.equal(fit$loglik, max(fit$sigma_profile$loglik)))
)

stopifnot(isTRUE(validate_simplified_sigma_profile(
  fit,
  require_interior = TRUE
)))

fit_bf <- BF_update_simplified_fash(fit)
stopifnot(
  identical(fit_bf$sigma_beta, fit$sigma_beta),
  identical(fit_bf$sigma_profile, fit$sigma_profile),
  identical(
    fit_bf$selected_sigma_on_boundary,
    fit$selected_sigma_on_boundary
  ),
  isTRUE(validate_simplified_sigma_profile(
    fit_bf,
    require_interior = TRUE
  ))
)

fixed_fit <- fit_simplified_fash(
  datasets = datasets,
  sigma_beta = 1,
  estimate_sigma = FALSE,
  scale_time = TRUE
)
stopifnot(
  nrow(fixed_fit$sigma_profile) == 1L,
  fixed_fit$sigma_profile$selected,
  fixed_fit$selected_sigma_on_boundary,
  isTRUE(validate_simplified_sigma_profile(
    fixed_fit,
    require_interior = FALSE
  ))
)

cat("Simplified-FASH sigma-profile tests passed.\n")
