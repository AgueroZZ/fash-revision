#!/usr/bin/env Rscript

# Extract a compact, page-ready cache from the full-data correlated FASH(1)
# fit downloaded from Midway3.
#
# The fit is 1.41 GB because `fash_data$Omega` holds one 16x16 precision matrix
# per unit. Everything the internal page needs is tiny, so it is extracted once
# here and the large object is never loaded again.
#
# Verifies against the numbers the Midway3 run reported before extracting.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

log_message <- function(...) {
  message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
  utils::flush.console()
}

workflowr_root <- find_workflowr_root()
alpha <- 0.05

# Reported by Slurm job 54276646; the extract must reproduce these exactly.
EXPECTED_UNITS <- 1009173L
EXPECTED_PI0 <- 0.9941774
EXPECTED_CALLS <- 43L

cache_directory <- file.path(
  workflowr_root, "output", "revision_simulations", "internal",
  "correlated_permutation_discoveries"
)
dir.create(cache_directory, recursive = TRUE, showWarnings = FALSE)

load_single <- function(path, label) {
  environment_holder <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment_holder)
  if (length(loaded) != 1L) {
    stop(path, " must contain exactly one object.")
  }
  log_message("Loaded ", label, ": `", loaded, "` from ", basename(path))
  list(object = environment_holder[[loaded]], name = loaded)
}

cumulative_calls <- function(lfdr, alpha) {
  ordering <- order(lfdr, method = "radix")
  ordering[cumsum(lfdr[ordering]) / seq_along(ordering) <= alpha]
}

# ------------------------------------------------------- correlated BF fit ---

correlated <- load_single(
  file.path(workflowr_root, "output", "dynamic_eQTL_real",
            "fash_fit1_corr_update.RData"),
  "correlated BF fit"
)
fit <- correlated$object

pair_keys <- names(fit$fash_data$data_list)
lfdr <- as.numeric(fit$lfdr)
names(lfdr) <- pair_keys
pi0 <- fit$prior_weights$prior_weight[fit$prior_weights$psd == 0]
selected <- cumulative_calls(lfdr, alpha)

log_message("units ", length(pair_keys), ", pi0 ", signif(pi0, 7),
            ", calls ", length(selected))
stopifnot(
  length(pair_keys) == EXPECTED_UNITS,
  abs(pi0 - EXPECTED_PI0) < 1e-6,
  length(selected) == EXPECTED_CALLS
)
log_message("Verified against the Midway3 run.")

n_time <- length(fit$fash_data$data_list[[selected[1]]]$y)
time_grid <- as.numeric(fit$fash_data$data_list[[selected[1]]]$x)

discoveries <- data.frame(
  correlated_rank = seq_along(selected),
  pair_key = pair_keys[selected],
  gene_id = sub("_.*$", "", pair_keys[selected]),
  variant_id = sub("^[^_]*_", "", pair_keys[selected]),
  correlated_lfdr = unname(lfdr[selected]),
  stringsAsFactors = FALSE
)

beta_hat <- t(vapply(fit$fash_data$data_list[selected],
                     function(unit) as.numeric(unit$y), numeric(n_time)))
# The correlated fit stores Omega, not S. Because the supplied correlation has
# a unit diagonal, Omega = D^{-1} C^{-1} D^{-1} implies diag(solve(Omega)) has
# entries se^2, so the marginal standard errors recover exactly.
standard_errors <- t(vapply(fit$fash_data$Omega[selected], function(omega) {
  sqrt(diag(solve(as.matrix(omega))))
}, numeric(n_time)))
rownames(beta_hat) <- rownames(standard_errors) <- discoveries$pair_key
colnames(beta_hat) <- colnames(standard_errors) <- paste0("t", time_grid)

posterior_weights <- as.matrix(fit$posterior_weights)[selected, , drop = FALSE]
rownames(posterior_weights) <- discoveries$pair_key

# One representative precision matrix, to confirm the correlation that was
# actually used matches the matrix we exported.
representative_omega <- as.matrix(fit$fash_data$Omega[[selected[1]]])
representative_se <- standard_errors[1, ]
# Omega = D^{-1} C^{-1} D^{-1}, so solve(Omega) = D C D and recovering C needs
# division by the standard errors, not multiplication.
implied_correlation <- diag(1 / representative_se) %*%
  solve(representative_omega) %*% diag(1 / representative_se)
implied_correlation <- (implied_correlation + t(implied_correlation)) / 2
if (max(abs(diag(implied_correlation) - 1)) > 1e-8) {
  stop("The correlation recovered from Omega does not have a unit diagonal.")
}

correlated_cache <- list(
  discoveries = discoveries,
  time_grid = time_grid,
  beta_hat = beta_hat,
  standard_errors = standard_errors,
  posterior_weights = posterior_weights,
  retained_psd = as.numeric(fit$prior_weights$psd),
  prior_weights = fit$prior_weights,
  psd_grid = as.numeric(fit$psd_grid),
  settings = fit$settings,
  pi0 = pi0,
  n_units = length(pair_keys),
  alpha = alpha,
  implied_correlation = implied_correlation
)

# Correlated lfdr and rank for every unit that the diagonal fit called, so the
# page can show both directions of the overlap.
diagonal_calls <- utils::read.csv(
  file.path(workflowr_root, "output", "revision_simulations", "internal",
            "residual_correlation_fash", "diagonal_bf_calls_alpha005.csv"),
  stringsAsFactors = FALSE
)
matched <- match(diagonal_calls$pair_key, pair_keys)
stopifnot(!anyNA(matched))
ordering <- order(lfdr, method = "radix")
rank_of <- integer(length(pair_keys))
rank_of[ordering] <- seq_along(ordering)
diagonal_calls$correlated_lfdr <- unname(lfdr[matched])
diagonal_calls$correlated_rank <- rank_of[matched]
diagonal_calls$correlated_call <- diagonal_calls$pair_key %in%
  discoveries$pair_key

top_context <- data.frame(
  correlated_rank = seq_len(5000),
  pair_key = pair_keys[ordering[seq_len(5000)]],
  correlated_lfdr = unname(lfdr[ordering[seq_len(5000)]]),
  stringsAsFactors = FALSE
)

rm(fit, correlated)
invisible(gc())

# ---------------------------------------------------- diagonal BF fit lfdr ---

diagonal <- load_single(
  file.path(workflowr_root, "output", "dynamic_eQTL_real",
            "fash_fit1_update.RData"),
  "diagonal BF fit"
)
diagonal_lfdr <- as.numeric(diagonal$object$lfdr)
names(diagonal_lfdr) <- names(diagonal$object$fash_data$data_list)
diagonal_se <- t(vapply(
  diagonal$object$fash_data$S[match(discoveries$pair_key,
                                    names(diagonal_lfdr))],
  as.numeric, numeric(n_time)
))
rownames(diagonal_se) <- discoveries$pair_key

discoveries$diagonal_lfdr <- unname(
  diagonal_lfdr[match(discoveries$pair_key, names(diagonal_lfdr))]
)
diagonal_ordering <- order(diagonal_lfdr, method = "radix")
diagonal_rank_of <- integer(length(diagonal_lfdr))
diagonal_rank_of[diagonal_ordering] <- seq_along(diagonal_ordering)
discoveries$diagonal_rank <- diagonal_rank_of[
  match(discoveries$pair_key, names(diagonal_lfdr))
]
discoveries$diagonal_call <- discoveries$pair_key %in% diagonal_calls$pair_key

# The standard errors must be identical: the correlated fit was built from the
# same t-adjusted values, only rearranged into Omega.
se_difference <- max(abs(standard_errors - diagonal_se))
log_message("max |SE recovered from Omega - diagonal SE| = ",
            signif(se_difference, 3))
stopifnot(se_difference < 1e-8)

rm(diagonal)
invisible(gc())

correlated_cache$discoveries <- discoveries
saveRDS(correlated_cache,
        file.path(cache_directory, "correlated_discovery_cache.rds"))
utils::write.csv(discoveries,
                 file.path(cache_directory, "correlated_discoveries.csv"),
                 row.names = FALSE)
utils::write.csv(diagonal_calls,
                 file.path(cache_directory,
                           "diagonal_calls_under_correlation.csv"),
                 row.names = FALSE)
utils::write.csv(top_context,
                 file.path(cache_directory, "correlated_top_units.csv"),
                 row.names = FALSE)

log_message("Discovered pairs: ", nrow(discoveries),
            "; unique genes: ", length(unique(discoveries$gene_id)),
            "; also called by the diagonal fit: ", sum(discoveries$diagonal_call))
log_message("Wrote cache to ", cache_directory)
