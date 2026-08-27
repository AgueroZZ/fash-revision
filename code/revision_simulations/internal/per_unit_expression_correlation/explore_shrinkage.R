#!/usr/bin/env Rscript

# Choose the shrinkage target and the shrinkage amount for a regularized
# per-unit Omega, from conditioning alone.
#
# Candidate targets are compared on their own conditioning and on whether they
# are outcome-selected. Candidate rules are compared on three things: the
# condition-number distribution they achieve, whether the amount of shrinkage
# is confounded with discovery status, and how much of the discovery-versus-null
# lag-1 gap survives.
#
# No discovery count enters this script. The amount of regularization must be
# fixed before any fit, or the likelihood is being tuned on the outcome.

r <- readRDS("output/revision_simulations/internal/per_unit_expression_correlation/per_unit_expression_correlation.rds")
A <- r$matrices$beta_scale
keys <- r$unit_table$pair_key
n <- dim(A)[1]; p <- 16
disc <- keys %in% r$panels$discovery1169
nul  <- keys %in% r$panels$null874
surv <- keys %in% r$panels$survivor43

cond <- function(m) { e <- eigen((m+t(m))/2, symmetric=TRUE, only.values=TRUE)$values; max(e)/min(e) }
lag1 <- function(m) mean(m[cbind(1:(p-1), 2:p)])

targets <- list(
  perm_null874  = r$common_c,
  mean_all_units = apply(A, c(2,3), mean),
  mean_null874   = apply(A[nul,,,drop=FALSE], c(2,3), mean)
)
cat("=== candidate shrinkage targets ===\n")
for (nm in names(targets)) {
  t0 <- targets[[nm]]; diag(t0) <- 1
  cat(sprintf("  %-15s lag-1 %.4f  mean-offdiag %.4f  cond %6.2f\n",
              nm, lag1(t0), mean(t0[row(t0)!=col(t0)]), cond(t0)))
}

target <- targets$mean_all_units; diag(target) <- 1
grid <- seq(0, 1, by = 0.005)

min_lambda <- function(m, cap) {
  for (l in grid) if (cond((1-l)*m + l*target) <= cap) return(l)
  1
}
cat("\n=== minimal per-unit lambda to reach a condition-number cap (target = mean of all units) ===\n")
cat(sprintf("%6s %8s %8s %8s %8s %8s %8s | %s\n", "cap","mean","median","q90","q99","max","frac>0","global lambda needed"))
res <- list()
for (cap in c(20, 30, 50)) {
  lam <- vapply(seq_len(n), function(j) min_lambda(A[j,,], cap), numeric(1))
  res[[as.character(cap)]] <- lam
  cat(sprintf("%6d %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f | %.3f\n",
      cap, mean(lam), median(lam), quantile(lam,.9), quantile(lam,.99), max(lam),
      mean(lam>0), max(lam)))
}

cat("\n=== is the shrinkage amount confounded with discovery status? ===\n")
own <- vapply(seq_len(n), function(j) lag1(A[j,,]), numeric(1))
for (cap in c(20,30,50)) {
  lam <- res[[as.character(cap)]]
  cat(sprintf("  cap %2d: Spearman(lambda, own lag-1) = %+0.3f | mean lambda: discovery %.3f, null %.3f, survivors %.3f\n",
      cap, cor(lam, own, method="spearman"), mean(lam[disc]), mean(lam[nul]), mean(lam[surv])))
}

cat("\n=== does the discovery-vs-null lag-1 gap survive shrinkage? ===\n")
cat(sprintf("  %-28s discovery %.4f  null %.4f  gap %.4f\n", "raw C_j",
            mean(own[disc]), mean(own[nul]), mean(own[nul])-mean(own[disc])))
for (cap in c(20,30,50)) {
  lam <- res[[as.character(cap)]]
  sh <- vapply(seq_len(n), function(j) lag1((1-lam[j])*A[j,,] + lam[j]*target), numeric(1))
  cat(sprintf("  %-28s discovery %.4f  null %.4f  gap %.4f  (%.0f%% of raw gap)\n",
      paste0("per-unit lambda, cap ", cap),
      mean(sh[disc]), mean(sh[nul]), mean(sh[nul])-mean(sh[disc]),
      100*(mean(sh[nul])-mean(sh[disc]))/(mean(own[nul])-mean(own[disc]))))
}
for (gl in c(0.2, 0.3, 0.5)) {
  sh <- (1-gl)*own + gl*lag1(target)
  cat(sprintf("  %-28s discovery %.4f  null %.4f  gap %.4f  (%.0f%% of raw gap)\n",
      paste0("global lambda = ", gl), mean(sh[disc]), mean(sh[nul]),
      mean(sh[nul])-mean(sh[disc]),
      100*(mean(sh[nul])-mean(sh[disc]))/(mean(own[nul])-mean(own[disc]))))
}
cat("\n=== condition number under a single global lambda ===\n")
for (gl in c(0.1, 0.2, 0.3, 0.5)) {
  cn <- vapply(seq_len(n), function(j) cond((1-gl)*A[j,,] + gl*target), numeric(1))
  cat(sprintf("  lambda %.2f: cond median %6.1f  q99 %8.1f  max %9.1f  frac>50 %.3f  non-PD %d\n",
      gl, median(cn), quantile(cn,.99), max(cn), mean(cn>50), sum(cn<0)))
}
