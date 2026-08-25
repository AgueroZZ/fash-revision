#!/usr/bin/env Rscript

# Quick comparison of three observation models on one fixed, outcome-selected
# panel: the most significant variant of every gene the diagonal BF fit
# discovered.
#
#   A  per-unit Omega_j, from that unit's own 400 permutation draws
#   B  one common C for the panel, averaged over the same per-unit estimates
#   C  the diagonal S = "SE" likelihood
#
# All three use identical FASH settings, so the only difference is the
# observation model.
#
# IMPORTANT: the panel is selected on the outcome. Absolute pi0 and discovery
# counts are therefore not FDR statements. Only the contrast across the three
# likelihoods is interpretable.

suppressPackageStartupMessages(library(fashr))
log_message <- function(...) {
  message(format(Sys.time(), "%H:%M:%S"), "  ", ...); utils::flush.console()
}
cache <- "output/revision_simulations/internal/correlated_permutation_discoveries"
num_cores <- as.integer(Sys.getenv("NCORES", "8"))
alpha <- 0.05

# ------------------------------------------------------- panel and null draws ---

meta <- utils::read.csv(file.path(cache, "panel_topgene_metadata.csv"),
                        stringsAsFactors = FALSE)
perm <- readRDS(file.path(
  "output/revision_simulations/internal/panel_null_topgene1169",
  "multigene_null_beta_covariance.rds"
))
draws <- perm$beta_draws$donor_residual_block_permutation
pair_keys <- dimnames(draws)[[1]]
n_time <- dim(draws)[2]
stopifnot(setequal(pair_keys, meta$pair_key))
meta <- meta[match(pair_keys, meta$pair_key), ]
log_message(length(pair_keys), " units, ", dim(draws)[3], " null draws each.")

# Per-unit null correlation, and the common average over this same panel.
per_unit_C <- lapply(seq_along(pair_keys), function(j) {
  C <- stats::cor(t(draws[j, , ]))
  C <- (C + t(C)) / 2; diag(C) <- 1; C
})
common_C <- Reduce(`+`, per_unit_C) / length(per_unit_C)
common_C <- (common_C + t(common_C)) / 2; diag(common_C) <- 1

min_eigen <- vapply(per_unit_C, function(C) {
  min(eigen(C, symmetric = TRUE, only.values = TRUE)$values)
}, numeric(1))
log_message("per-unit C: minimum eigenvalue ranges ",
            signif(min(min_eigen), 3), " to ", signif(max(min_eigen), 3))
if (min(min_eigen) <= 1e-8) stop("A per-unit correlation is not positive definite.")

lag1 <- vapply(per_unit_C, function(C) mean(C[cbind(1:(n_time-1), 2:n_time)]),
               numeric(1))
log_message("per-unit lag-1: mean ", round(mean(lag1), 3),
            ", sd ", round(sd(lag1), 3),
            ", range [", round(min(lag1), 3), ", ", round(max(lag1), 3), "]")
log_message("common C lag-1: ", round(mean(common_C[cbind(1:(n_time-1), 2:n_time)]), 3))

# ---------------------------------------------------------- observed data ---

full <- new.env(parent = emptyenv())
load("output/dynamic_eQTL_real/fash_fit1_all.RData", envir = full)
fit_full <- full$fash_fit1
idx <- match(pair_keys, names(fit_full$fash_data$data_list))
stopifnot(!anyNA(idx))
beta_hat <- t(vapply(fit_full$fash_data$data_list[idx],
                     function(u) as.numeric(u$y), numeric(n_time)))
se <- t(vapply(fit_full$fash_data$S[idx], as.numeric, numeric(n_time)))
time_grid <- as.numeric(fit_full$fash_data$data_list[[idx[1]]]$x)
settings <- fit_full$settings
psd_grid <- as.numeric(fit_full$psd_grid)
rownames(beta_hat) <- rownames(se) <- pair_keys
rm(full, fit_full); invisible(gc())
stopifnot(identical(time_grid, as.numeric(0:15)))

data_list <- lapply(seq_along(pair_keys), function(j) {
  data.frame(beta = beta_hat[j, ], time = time_grid, SE = se[j, ])
})
names(data_list) <- pair_keys

omega_from <- function(C, j) {
  P <- solve(C) * tcrossprod(1 / se[j, ])
  (P + t(P)) / 2
}

# ------------------------------------------------------------------- fits ---

fit_one <- function(label, omega = NULL, use_se = FALSE) {
  log_message("fitting ", label)
  started <- Sys.time()
  raw <- if (use_se) {
    fashr::fash(Y = "beta", smooth_var = "time", S = "SE", data_list = data_list,
                num_basis = settings$num_basis, order = settings$order,
                betaprec = settings$betaprec, pred_step = settings$pred_step,
                penalty = settings$penalty, grid = psd_grid,
                num_cores = num_cores, verbose = FALSE)
  } else {
    fashr::fash(Y = "beta", smooth_var = "time", Omega = omega,
                data_list = data_list,
                num_basis = settings$num_basis, order = settings$order,
                betaprec = settings$betaprec, pred_step = settings$pred_step,
                penalty = settings$penalty, grid = psd_grid,
                num_cores = num_cores, verbose = FALSE)
  }
  names(raw$lfdr) <- pair_keys
  bf <- fashr::BF_update(raw)
  bf_ok <- !is.null(bf$BF) && !identical(bf$prior_weights, raw$prior_weights)
  if (bf_ok) names(bf$lfdr) <- pair_keys
  log_message("  done in ", round(as.numeric(difftime(Sys.time(), started, units = "secs"))), "s",
              if (bf_ok) "" else "  (BF_update did not apply)")
  list(raw = raw, bf = if (bf_ok) bf else NULL, bf_available = bf_ok)
}

fits <- list(
  per_unit = fit_one("A: per-unit Omega",
                     omega = lapply(seq_along(pair_keys),
                                    function(j) omega_from(per_unit_C[[j]], j))),
  common   = fit_one("B: common C",
                     omega = lapply(seq_along(pair_keys),
                                    function(j) omega_from(common_C, j))),
  diagonal = fit_one("C: diagonal SE", use_se = TRUE)
)

# ---------------------------------------------------------------- summary ---

pi0_of <- function(f) {
  v <- f$prior_weights$prior_weight[f$prior_weights$psd == 0]
  if (length(v) == 1L) v else 0
}
calls_of <- function(l) {
  o <- order(l, method = "radix"); sum(cumsum(l[o]) / seq_along(o) <= alpha)
}
summary_rows <- do.call(rbind, lapply(names(fits), function(nm) {
  f <- fits[[nm]]
  rbind(
    data.frame(model = nm, stage = "raw", bf_available = TRUE,
               pi0 = pi0_of(f$raw), mean_lfdr = mean(f$raw$lfdr),
               calls = calls_of(f$raw$lfdr),
               components = nrow(f$raw$prior_weights)),
    data.frame(model = nm, stage = "bf", bf_available = f$bf_available,
               pi0 = if (f$bf_available) pi0_of(f$bf) else NA_real_,
               mean_lfdr = if (f$bf_available) mean(f$bf$lfdr) else NA_real_,
               calls = if (f$bf_available) calls_of(f$bf$lfdr) else NA_integer_,
               components = if (f$bf_available) nrow(f$bf$prior_weights) else NA_integer_)
  )
}))
utils::write.csv(summary_rows, file.path(cache, "three_likelihood_summary.csv"),
                 row.names = FALSE)

unit_table <- data.frame(
  pair_key = pair_keys, gene_id = meta$gene_id, variant_id = meta$variant_id,
  per_unit_lag1 = lag1,
  diagonal_full_lfdr = meta$diagonal_bf_lfdr,
  stringsAsFactors = FALSE
)
for (nm in names(fits)) {
  unit_table[[paste0(nm, "_raw_lfdr")]] <- as.numeric(fits[[nm]]$raw$lfdr)
  unit_table[[paste0(nm, "_bf_lfdr")]] <-
    if (fits[[nm]]$bf_available) as.numeric(fits[[nm]]$bf$lfdr) else NA_real_
}
utils::write.csv(unit_table, file.path(cache, "three_likelihood_units.csv"),
                 row.names = FALSE)
saveRDS(list(per_unit_C = per_unit_C, common_C = common_C, lag1 = lag1,
             summary = summary_rows, units = unit_table),
        file.path(cache, "three_likelihood_fits.rds"))

print(summary_rows, digits = 5, row.names = FALSE)
log_message("wrote three_likelihood_{summary,units}.csv")
