#!/usr/bin/env Rscript

# Deterministic tests for the R4 time-correlated expression-error generator.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the shared revision simulation functions.")
}

expect_error <- function(code) {
  observed <- FALSE
  tryCatch(
    force(code),
    error = function(e) {
      observed <<- TRUE
    }
  )
  if (!observed) {
    stop("Expected an error, but the expression succeeded.")
  }
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "shared",
  "simulation_functions.R"
))

lag1 <- make_lag1_correlation(16, -0.09)
stopifnot(
  identical(dim(lag1), c(16L, 16L)),
  max(abs(diag(lag1) - 1)) < 1e-12,
  max(abs(lag1 - t(lag1))) < 1e-12,
  max(abs(lag1[row(lag1) == col(lag1) + 1] + 0.09)) < 1e-12,
  max(abs(lag1[col(lag1) == row(lag1) + 1] + 0.09)) < 1e-12,
  max(abs(lag1[abs(row(lag1) - col(lag1)) > 1])) < 1e-12,
  min(eigen(lag1, symmetric = TRUE, only.values = TRUE)$values) > 0
)

expect_error(make_lag1_correlation(1, -0.09))
expect_error(make_lag1_correlation(16, 0.5))
expect_error(validate_time_correlation(diag(3), n_time = 4))
asymmetric <- diag(4)
asymmetric[1, 2] <- 0.1
expect_error(validate_time_correlation(asymmetric, n_time = 4))
non_unit_diagonal <- diag(c(1, 1, 1, 0.9))
expect_error(validate_time_correlation(non_unit_diagonal, n_time = 4))
indefinite <- matrix(0.9, nrow = 4, ncol = 4)
diag(indefinite) <- 1
indefinite[1, 2] <- indefinite[2, 1] <- -0.9
expect_error(validate_time_correlation(indefinite, n_time = 4))

small_G <- matrix(
  c(
    0, 1, 2, 0,
    1, 2, 0, 1,
    2, 0, 1, 2
  ),
  nrow = 4,
  ncol = 3,
  dimnames = list(
    sprintf("donor_%02d", 1:4),
    sprintf("variant_%02d", 1:3)
  )
)
small_beta <- matrix(
  seq(-0.4, 0.4, length.out = 12),
  nrow = 3,
  ncol = 4,
  dimnames = list(colnames(small_G), sprintf("time_%02d", 1:4))
)
small_covariates <- matrix(
  c(-1.2, -0.3, 0.4, 1.1),
  ncol = 1,
  dimnames = list(rownames(small_G), "PC1")
)

simulate_reference_iid_expression <- function(G,
                                              beta_matrix,
                                              covariates,
                                              expression_noise_sd,
                                              covariate_effect_sd,
                                              intercept_sd,
                                              seed) {
  set.seed(seed)
  n_donors <- nrow(G)
  n_variants <- ncol(G)
  n_time <- ncol(beta_matrix)
  n_covariates <- ncol(covariates)
  expression <- array(NA_real_, dim = c(n_donors, n_variants, n_time))
  intercepts <- matrix(
    rnorm(n_variants * n_time, mean = 0, sd = intercept_sd),
    nrow = n_variants,
    ncol = n_time
  )
  covariate_effects <- array(
    NA_real_,
    dim = c(n_covariates, n_variants, n_time)
  )

  for (tt in seq_len(n_time)) {
    genetic_mean <- sweep(G, 2, beta_matrix[, tt], `*`)
    gamma <- matrix(
      rnorm(n_covariates * n_variants, mean = 0, sd = covariate_effect_sd),
      nrow = n_covariates,
      ncol = n_variants
    )
    covariate_effects[, , tt] <- gamma
    covariate_mean <- covariates %*% gamma
    noise <- matrix(
      rnorm(n_donors * n_variants, mean = 0, sd = expression_noise_sd),
      nrow = n_donors,
      ncol = n_variants
    )
    expression[, , tt] <- sweep(
      genetic_mean + covariate_mean + noise,
      2,
      intercepts[, tt],
      `+`
    )
  }

  list(
    expression = expression,
    intercepts = intercepts,
    covariate_effects = covariate_effects
  )
}

iid_default <- simulate_eqtl_expression_from_genotypes(
  G = small_G,
  beta_matrix = small_beta,
  time_grid = 0:3,
  covariates = small_covariates,
  expression_noise_sd = 1,
  covariate_effect_sd = 0.5,
  intercept_sd = 0.25,
  expression_error_correlation = NULL,
  seed = 8917
)
iid_identity <- simulate_eqtl_expression_from_genotypes(
  G = small_G,
  beta_matrix = small_beta,
  time_grid = 0:3,
  covariates = small_covariates,
  expression_noise_sd = 1,
  covariate_effect_sd = 0.5,
  intercept_sd = 0.25,
  expression_error_correlation = diag(4),
  seed = 8917
)
iid_reference <- simulate_reference_iid_expression(
  G = small_G,
  beta_matrix = small_beta,
  covariates = small_covariates,
  expression_noise_sd = 1,
  covariate_effect_sd = 0.5,
  intercept_sd = 0.25,
  seed = 8917
)
stopifnot(
  identical(iid_default$expression, iid_identity$expression),
  identical(iid_default$intercepts, iid_identity$intercepts),
  identical(iid_default$covariate_effects, iid_identity$covariate_effects),
  identical(unname(iid_default$expression), iid_reference$expression),
  identical(unname(iid_default$intercepts), iid_reference$intercepts),
  identical(
    unname(iid_default$covariate_effects),
    iid_reference$covariate_effects
  )
)

set.seed(1921)
n_donors <- 200L
n_variants <- 200L
n_time <- 16L
large_G <- matrix(
  sample(0:2, n_donors * n_variants, replace = TRUE),
  nrow = n_donors,
  ncol = n_variants,
  dimnames = list(
    sprintf("donor_%03d", seq_len(n_donors)),
    sprintf("variant_%03d", seq_len(n_variants))
  )
)
large_beta <- matrix(
  0,
  nrow = n_variants,
  ncol = n_time,
  dimnames = list(colnames(large_G), sprintf("time_%02d", seq_len(n_time)))
)
correlated <- simulate_eqtl_expression_from_genotypes(
  G = large_G,
  beta_matrix = large_beta,
  time_grid = 0:15,
  expression_noise_sd = 1,
  covariate_effect_sd = 0,
  intercept_sd = 0,
  expression_error_correlation = lag1,
  seed = 9517
)
noise_matrix <- matrix(
  correlated$expression,
  nrow = n_donors * n_variants,
  ncol = n_time
)
empirical_correlation <- stats::cor(noise_matrix)
empirical_lag1 <- mean(empirical_correlation[cbind(1:15, 2:16)])
empirical_higher_lags <- unlist(lapply(2:15, function(lag) {
  empirical_correlation[cbind(seq_len(n_time - lag), (lag + 1L):n_time)]
}))
stopifnot(
  abs(empirical_lag1 - (-0.09)) < 0.02,
  max(abs(empirical_higher_lags)) < 0.02
)

cat("R4 correlated-error generation tests passed.\n")
