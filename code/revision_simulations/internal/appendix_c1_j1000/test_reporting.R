# Independently recompute selection and estimands from every retained fitted score.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
source("code/revision_simulations/internal/appendix_c1_j1000/reporting.R")
report <- c1j_load_report()
results <- report$tables$alpha_by_seed
timing <- numeric()
for (seed in report$contract$seeds) {
  record <- readRDS(file.path(report$directory, "replicates", paste0("seed_", seed, ".rds")))
  stopifnot(length(record$warnings) == 0L,
    identical(record$input_sha256, digest::digest(record$data_bundle, algo = "sha256")),
    length(record$data_bundle$data_list) == 1000L)
  timing <- c(timing, record$elapsed_seconds)
  truth <- record$data_bundle$truth
  for (model in 1:2) for (stage in c("raw", "bf")) {
    positive <- if (model == 1L) truth$dynamic_true else truth$nonlinear_true
    stopifnot(sum(positive) == if (model == 1L) 200L else 100L)
    scores <- record$fit_bundle[[paste0("order", model)]][[stage]]$lfdr
    ranked <- order(scores, seq_along(scores))
    cumulative <- cumsum(scores[ranked]) / seq_along(scores)
    for (threshold in report$contract$alpha) {
      n <- max(c(0L, which(cumulative <= threshold)))
      selected <- head(ranked, n)
      actual <- c(discoveries = n, true_discoveries = sum(positive[selected]),
        false_discoveries = sum(!positive[selected]),
        realized_fdp = sum(!positive[selected]) / max(1, n),
        power = sum(positive[selected]) / sum(positive))
      row <- results[results$seed == seed & results$order == paste0("IWP", model) &
        results$stage == if (stage == "raw") "Raw EB" else "BF updated", ]
      row <- row[abs(row$alpha - threshold) < 1e-12, ]
      stopifnot(nrow(row) == 1L, max(abs(unlist(row[names(actual)]) - actual)) < 1e-14)
    }
  }
}
summary <- report$tables$alpha_summary
raw <- summary[summary$stage == "Raw EB", ]
bf <- summary[summary$stage == "BF updated", ]
key <- function(x) paste(x$order, x$alpha)
bf <- bf[match(key(raw), key(bf)), ]
loss <- raw$mean_power - bf$mean_power
excess <- bf$mean_empirical_fdr - bf$alpha
cat("Total fitting seconds:", sum(timing), "\n")
cat("Maximum mean-power loss (percentage points):", max(loss) * 100, "\n")
cat("BF-adjusted mean-FDR exceedances:", sum(excess > 1e-12), "of", length(excess), "\n")
print(bf[which.max(excess), c("order", "alpha", "mean_empirical_fdr")])
print(summary[abs(summary$alpha - 0.05) < 1e-12, ])
write.csv(data.frame(order = raw$order, alpha = raw$alpha,
  power_loss = loss, adjusted_fdr_excess = excess),
  file.path(report$directory, "claim_checks.csv"), row.names = FALSE)
cat("Independent selection, denominators, and completed-cache validation passed.\n")
