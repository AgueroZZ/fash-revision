# Verify exact-integration equivalence, selection edge cases, and estimands.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
source("code/revision_simulations/internal/appendix_c1_mc10/mc10_helpers.R")
source("code/revision_simulations/internal/appendix_c1_mc10/gaussian_likelihood.R")
require_appendix_b_fashr()
reference <- readRDS("data/appendixB/fashr0143/focused_example.rds")
checks <- list()
for (order in 1:2) {
  fit <- reference$fit_bundle[[paste0("order", order)]]$raw
  setup <- c1_gaussian_setup(fit$fash_data$data_list[[1]]$x, order, 20L)
  timing <- system.time(L <- c1_gaussian_likelihood(fit$fash_data, fit$psd_grid, setup))
  error <- max(abs(L - fit$L_matrix))
  stopifnot(error < 1e-7)
  for (penalty in c(1, 10)) {
    direct <- c1_fit_from_likelihood(fit$L_matrix, fit$fash_data, fit$psd_grid, penalty, order, 20L)
    fast <- c1_fit_from_likelihood(L, fit$fash_data, fit$psd_grid, penalty, order, 20L)
    score_error <- max(abs(direct$raw$lfdr - fast$raw$lfdr), abs(direct$bf$lfdr - fast$bf$lfdr))
    stopifnot(score_error < 1e-5,
              identical(cumulative_lfdr_calls(direct$bf$lfdr)$indices,
                        cumulative_lfdr_calls(fast$bf$lfdr)$indices))
    checks[[paste(order, penalty, sep = ":")]] <- c(likelihood_max_error = error,
      lfdr_max_error = score_error, likelihood_seconds = timing[["elapsed"]])
  }
  # Independently exercise the installed TMB routine across all truth classes
  # on the denser component grid, including the exact null component.
  indices <- c(1L, 300L, 961L, 1000L, 1081L, 1200L)
  data <- fit$fash_data
  data$data_list <- data$data_list[indices]
  data$S <- data$S[indices]
  dense_grid <- sort(c(0, exp(-0.5 * seq(0, 10, by = 0.1))))
  exact <- c1_gaussian_likelihood(data, dense_grid, setup)
  native <- fashr:::fash_L_compute(data, num_cores = 1L, grid = dense_grid,
                                 num_basis = 20L, betaprec = 0, order = order, pred_step = 1)
  stopifnot(max(abs(exact - native)) < 1e-7)
  checks[[paste0("dense_order", order)]] <- max(abs(exact - native))
}
zero <- evaluate_appendix_b_discoveries(c(0.9, 0.8), c(FALSE, TRUE), 0.05)
stopifnot(zero$discoveries == 0L, zero$realized_fdp == 0, zero$power == 0)
toy <- data.frame(seed = 1:10, setting = "test", rate = c(0, rep(0.5, 9)),
                  discoveries = c(0, rep(2, 9)))
summary <- c1_summary(toy, "setting", c("rate", "discoveries"))
stopifnot(summary$mean_rate == 0.45, summary$rate_min == 0, summary$rate_max == 0.5)
dir.create("output/revision_simulations/internal/appendix_c1_mc10_validation", recursive = TRUE, showWarnings = FALSE)
attr(checks, "integration_sha256") <- c1_sha("code/revision_simulations/internal/appendix_c1_mc10/gaussian_likelihood.R")
saveRDS(checks, "output/revision_simulations/internal/appendix_c1_mc10_validation/equivalence.rds")
print(checks)
cat("Exact Gaussian integration and Monte Carlo summary checks passed.\n")
