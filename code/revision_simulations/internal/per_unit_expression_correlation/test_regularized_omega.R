#!/usr/bin/env Rscript
# Checks on the regularization and Omega algebra used by Midway3 script
# scripts/23_per_unit_regularized_fash.R, run against the cached local C_j.

n_time <- 16L
r <- readRDS(paste0("output/revision_simulations/internal/",
                    "per_unit_expression_correlation/",
                    "per_unit_expression_correlation.rds"))
A <- r$matrices$beta_scale
lambda <- 0.1

ok <- function(condition, label) {
  if (!isTRUE(condition)) stop("FAILED: ", label)
  message("ok  ", label)
}
lag_one <- function(m) mean(m[cbind(seq_len(n_time - 1L), 2:n_time)])

# Target as script 23 builds it: mean over all units, symmetrized, unit diagonal.
target <- apply(A, c(2L, 3L), mean)
target <- (target + t(target)) / 2
diag(target) <- 1
ok(min(eigen(target, symmetric = TRUE, only.values = TRUE)$values) > 1e-10,
   "the all-unit mean target is strictly positive definite")

# The chunked accumulation in script 23 must equal the plain mean.
accumulated <- colSums(A, dims = 1L) / dim(A)[1L]
accumulated <- (accumulated + t(accumulated)) / 2
diag(accumulated) <- 1
ok(max(abs(accumulated - target)) < 1e-12,
   "colSums(dims = 1) accumulation equals apply(mean)")

regularize <- function(m, value) {
  shrunk <- (1 - value) * m + value * target
  shrunk <- (shrunk + t(shrunk)) / 2
  diag(shrunk) <- 1
  shrunk
}

# The convex combination already has a unit diagonal; forcing it must be a no-op.
raw_diag_error <- max(vapply(seq_len(dim(A)[1L]), function(j) {
  max(abs(diag((1 - lambda) * A[j, , ] + lambda * target) - 1))
}, numeric(1L)))
ok(raw_diag_error < 1e-12,
   "forcing the unit diagonal is a no-op for a convex combination")

# Every regularized matrix is positive definite at lambda = 0.1.
extremes <- vapply(seq_len(dim(A)[1L]), function(j) {
  values <- eigen(regularize(A[j, , ], lambda), symmetric = TRUE,
                  only.values = TRUE)$values
  c(min(values), max(values))
}, numeric(2L))
ok(min(extremes[1L, ]) > 1e-10,
   "every regularized correlation is positive definite at lambda = 0.1")
message("    min eigenvalue ", signif(min(extremes[1L, ]), 4),
        ", max condition ", signif(max(extremes[2L, ] / extremes[1L, ]), 5))

# The linear-attenuation claim the log reports must hold exactly.
own <- vapply(seq_len(dim(A)[1L]), function(j) lag_one(A[j, , ]), numeric(1L))
shrunk <- vapply(seq_len(dim(A)[1L]), function(j) {
  lag_one(regularize(A[j, , ], lambda))
}, numeric(1L))
predicted <- (1 - lambda) * own + lambda * lag_one(target)
ok(max(abs(shrunk - predicted)) < 1e-12,
   "lag-1 shrinks exactly as (1 - lambda) * own + lambda * target")

disc <- r$unit_table$pair_key %in% r$panels$discovery1169
nul <- r$unit_table$pair_key %in% r$panels$null874
gap_raw <- mean(own[nul]) - mean(own[disc])
gap_shrunk <- mean(shrunk[nul]) - mean(shrunk[disc])
ok(abs(gap_shrunk - (1 - lambda) * gap_raw) < 1e-12,
   "the discovery-versus-null gap is attenuated by exactly (1 - lambda)")
message("    gap ", signif(gap_raw, 4), " -> ", signif(gap_shrunk, 4),
        " (", round(100 * gap_shrunk / gap_raw), "% retained)")

# Omega as script 23 builds it must equal the explicit sandwich form.
set.seed(20260826)
probe <- sample(dim(A)[1L], 25L)
se <- abs(rnorm(n_time, 0.4, 0.1)) + 0.05
worst <- 0
for (j in probe) {
  shrunk_matrix <- regularize(A[j, , ], lambda)
  fast <- solve(shrunk_matrix) * tcrossprod(1 / se)
  explicit <- diag(1 / se) %*% solve(shrunk_matrix) %*% diag(1 / se)
  worst <- max(worst, max(abs(fast - explicit)))
}
ok(worst < 1e-9,
   "solve(C) * tcrossprod(1/se) equals diag(1/se) C^-1 diag(1/se)")
message("    worst absolute difference ", signif(worst, 3))

# Omega must be positive definite and must invert back to the covariance.
worst_roundtrip <- 0
for (j in probe) {
  shrunk_matrix <- regularize(A[j, , ], lambda)
  omega <- solve(shrunk_matrix) * tcrossprod(1 / se)
  omega <- (omega + t(omega)) / 2
  chol(omega)
  covariance <- diag(se) %*% shrunk_matrix %*% diag(se)
  worst_roundtrip <- max(worst_roundtrip, max(abs(solve(omega) - covariance)))
}
ok(worst_roundtrip < 1e-9,
   "Omega is Cholesky-factorable and inverts to diag(se) C diag(se)")
message("    worst round-trip difference ", signif(worst_roundtrip, 3))

message("All regularized-Omega checks passed.")
