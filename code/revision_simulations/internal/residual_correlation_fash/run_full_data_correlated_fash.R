#!/usr/bin/env Rscript

# Full-data FASH refit with a shared cross-time error correlation.
#
# This mirrors `code/01_fash.R` exactly, with one change: the diagonal
# standard-error likelihood `S = "SE"` is replaced by a per-unit precision
#
#   Omega_j = diag(1 / se_j) %*% solve(C) %*% diag(1 / se_j),
#
# where `se_j` are the same t-adjusted standard errors already stored in
# `datasets_corrected.RData` and `C` is a shared 16 x 16 correlation matrix.
# Every other setting - grid, num_basis, order, betaprec, pred_step, penalty -
# is identical to the retained manuscript fit, so the output is directly
# comparable to `fash_fit1_all.RData` / `fash_fit1_update.RData`.
#
# Usage:
#   Rscript --vanilla run_full_data_correlated_fash.R \
#     --correlation permutation_null_correlation_mg090.csv \
#     --datasets    results/fash_fit1_all.RData \
#     --output-dir  results/correlated_permutation \
#     --order 1 --num-cores 16
#
# Validate first with a short run before committing the full job:
#   ... --smoke 2000 --output-dir results/correlated_permutation_smoke

get_arg <- function(name, default = NULL) {
  arguments <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(arguments, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(arguments[equals_hit[1L]], nchar(equals_prefix) + 1L))
  }
  hit <- which(arguments == name)
  if (length(hit) == 0L || hit[1L] == length(arguments)) {
    return(default)
  }
  arguments[hit[1L] + 1L]
}

has_flag <- function(name) {
  any(commandArgs(trailingOnly = TRUE) == name)
}

log_message <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", ...)
  utils::flush.console()
}

correlation_path <- get_arg("--correlation")
datasets_path <- get_arg("--datasets", "./results/datasets_corrected.RData")
output_directory <- get_arg("--output-dir", "./results/correlated")
order <- as.integer(get_arg("--order", "1"))
penalty <- as.numeric(get_arg("--penalty", "10"))
num_cores <- as.integer(get_arg("--num-cores", "16"))
smoke <- as.integer(get_arg("--smoke", "0"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
force <- has_flag("--force")

if (is.null(correlation_path) || !file.exists(correlation_path)) {
  stop("--correlation must point to an existing 16 x 16 correlation CSV or RDS.")
}
if (!file.exists(datasets_path)) {
  stop("--datasets must point to an existing datasets_corrected.RData ",
       "or saved FASH fit.")
}
if (is.na(order) || !order %in% c(1L, 2L) ||
    is.na(penalty) || !is.finite(penalty) || penalty < 1 ||
    is.na(num_cores) || num_cores < 1L ||
    is.na(smoke) || smoke < 0L ||
    is.na(alpha) || alpha <= 0 || alpha >= 1) {
  stop("Invalid command-line arguments.")
}

suppressPackageStartupMessages(library(fashr))
log_message("fashr ", as.character(utils::packageVersion("fashr")),
            ", ", R.version.string)

# ------------------------------------------------------------- correlation ---

read_correlation <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    matrix_form <- readRDS(path)
  } else {
    frame <- utils::read.csv(path, stringsAsFactors = FALSE,
                             check.names = FALSE)
    # Tolerate a leading label column written by write.csv(row.names = TRUE).
    numeric_columns <- vapply(frame, function(column) {
      !any(is.na(suppressWarnings(as.numeric(as.character(column)))))
    }, logical(1))
    frame <- frame[, numeric_columns, drop = FALSE]
    matrix_form <- as.matrix(frame)
  }
  storage.mode(matrix_form) <- "double"
  dimnames(matrix_form) <- NULL
  matrix_form
}

correlation <- read_correlation(correlation_path)
if (nrow(correlation) != ncol(correlation) || nrow(correlation) != 16L ||
    any(!is.finite(correlation))) {
  stop("The correlation matrix must be a finite 16 x 16 matrix.")
}
if (max(abs(correlation - t(correlation))) > 1e-8) {
  stop("The correlation matrix is not symmetric.")
}
if (max(abs(diag(correlation) - 1)) > 1e-8) {
  stop("The correlation matrix does not have a unit diagonal.")
}
correlation <- (correlation + t(correlation)) / 2
diag(correlation) <- 1
eigenvalues <- eigen(correlation, symmetric = TRUE, only.values = TRUE)$values
if (min(eigenvalues) <= 1e-10) {
  stop("The correlation matrix is not strictly positive definite.")
}
log_message("Correlation: minimum eigenvalue ", signif(min(eigenvalues), 6),
            ", condition number ",
            signif(max(eigenvalues) / min(eigenvalues), 6),
            ", mean off-diagonal ",
            signif(mean(correlation[upper.tri(correlation)]), 6))
correlation_precision <- solve(correlation)
correlation_precision <- (correlation_precision +
                            t(correlation_precision)) / 2

# ----------------------------------------------------------------- datasets ---

# `--datasets` accepts either form, detected automatically:
#   (a) datasets_corrected.RData, holding a list named `datasets`; or
#   (b) any saved FASH fit, e.g. fash_fit1_all.RData, whose
#       fash_data$data_list and fash_data$S already carry the identical betas
#       and t-adjusted standard errors.
# Form (b) is preferred: it is the same object the manuscript fit used, so the
# likelihood inputs cannot drift, and it does not depend on an intermediate
# file that may not have been retained.
log_message("Loading ", datasets_path, " (this needs several GB).")
source_environment <- new.env(parent = emptyenv())
loaded <- load(datasets_path, envir = source_environment)
if (length(loaded) != 1L) {
  stop(datasets_path, " must contain exactly one object.")
}
source_object <- source_environment[[loaded]]

if (identical(loaded, "datasets") ||
      (is.list(source_object) && !inherits(source_object, "fash") &&
         is.data.frame(source_object[[1]]) &&
         all(c("beta", "SE") %in% names(source_object[[1]])))) {
  datasets <- source_object
  log_message("Input recognised as a `datasets` list: ", length(datasets),
              " gene-variant units.")
} else if (is.list(source_object) &&
             !is.null(source_object$fash_data$data_list) &&
             !is.null(source_object$fash_data$S)) {
  unit_data <- source_object$fash_data$data_list
  unit_se <- source_object$fash_data$S
  if (length(unit_data) != length(unit_se)) {
    stop("The FASH fit has mismatched data_list and S lengths.")
  }
  log_message("Input recognised as a FASH fit (`", loaded, "`): ",
              length(unit_data), " units; rebuilding beta/time/SE.")
  datasets <- Map(function(unit, se) {
    data.frame(beta = as.numeric(unit$y), time = as.numeric(unit$x),
               SE = as.numeric(se))
  }, unit_data, unit_se)
  names(datasets) <- names(unit_data)
} else {
  stop(datasets_path, " is neither a `datasets` list nor a FASH fit ",
       "carrying fash_data$data_list and fash_data$S.")
}
rm(source_object, source_environment)
invisible(gc())

if (smoke > 0L) {
  if (smoke > length(datasets)) {
    stop("--smoke exceeds the number of available units.")
  }
  datasets <- datasets[seq_len(smoke)]
  log_message("SMOKE RUN: restricted to the first ", length(datasets),
              " units. Do not use this output as a result.")
}

expected_time <- 0:15
bad <- vapply(datasets, function(unit) {
  !all(c("beta", "time", "SE") %in% names(unit)) ||
    nrow(unit) != 16L ||
    !identical(as.numeric(unit$time), as.numeric(expected_time)) ||
    any(!is.finite(unit$beta)) ||
    any(!is.finite(unit$SE)) || any(unit$SE <= 0)
}, logical(1))
if (any(bad)) {
  stop(sum(bad), " units are invalid or are not on the day 0-15 grid ",
       "in ascending order; the correlation matrix would be misaligned. ",
       "First offender: ", names(datasets)[which(bad)[1]])
}
log_message("All units validated on the day 0-15 grid.")

# -------------------------------------------------------------------- Omega ---

log_message("Building per-unit precision matrices (about ",
            round(2264 * length(datasets) / 1024^2), " MB).")
Omega <- lapply(datasets, function(unit) {
  inverse_se <- 1 / unit$SE
  precision <- correlation_precision * tcrossprod(inverse_se)
  (precision + t(precision)) / 2
})
names(Omega) <- NULL

set.seed(20260821)
check_index <- sort(sample(length(Omega), min(1000L, length(Omega))))
not_pd <- vapply(check_index, function(index) {
  !isTRUE(tryCatch({
    chol(Omega[[index]])
    TRUE
  }, error = function(error) FALSE))
}, logical(1))
if (any(not_pd)) {
  stop(sum(not_pd), " sampled precision matrices are not positive definite.")
}
log_message("Sampled ", length(check_index),
            " precision matrices; all positive definite.")
invisible(gc())

# ---------------------------------------------------------------------- fit ---

log_precision <- seq(0, 10, by = 0.2)
fine_grid <- sort(c(0, exp(-0.5 * log_precision)))

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
raw_path <- file.path(output_directory,
                      paste0("fash_fit", order, "_corr_all.RData"))
update_path <- file.path(output_directory,
                         paste0("fash_fit", order, "_corr_update.RData"))
if (!force && (file.exists(raw_path) || file.exists(update_path))) {
  stop("Output already exists in ", output_directory,
       ". Choose a new --output-dir or pass --force.")
}

log_message("Fitting FASH(", order, "), penalty ", penalty, ", ",
            length(fine_grid), " grid points, ", num_cores, " cores.")
started <- Sys.time()
fash_fit_corr <- fash(
  Y = "beta",
  smooth_var = "time",
  Omega = Omega,
  data_list = datasets,
  num_basis = 20,
  order = order,
  betaprec = 0,
  pred_step = 1,
  penalty = penalty,
  grid = fine_grid,
  num_cores = num_cores,
  verbose = TRUE
)
elapsed <- difftime(Sys.time(), started, units = "mins")
log_message("Fit complete in ", format(elapsed, digits = 5), ".")
save(fash_fit_corr, file = raw_path)
log_message("Wrote ", raw_path)

# `BF_update()` warns and returns the ORIGINAL object when every Bayes factor
# is NA/NaN, which happens when the raw fit puts zero total prior mass on the
# non-null components so the conditional alternative mixture is 0/0. Detect
# that and report the BF stage as unavailable rather than silently copying the
# raw numbers into a column labelled "bf".
log_message("Applying BF_update().")
fash_fit_corr_update <- fashr::BF_update(fash_fit_corr)
bf_available <- !is.null(fash_fit_corr_update$BF) &&
  !identical(fash_fit_corr_update$prior_weights, fash_fit_corr$prior_weights)
if (!bf_available) {
  log_message("BF_update() did not apply: the raw fit assigns no prior mass ",
              "to any non-null component, so the BF normalisation is 0/0. ",
              "The BF stage is reported as unavailable.")
} else {
  save(fash_fit_corr_update, file = update_path)
  log_message("Wrote ", update_path)
}

# ------------------------------------------------------------------ summary ---

null_weight <- function(fit) {
  value <- fit$prior_weights$prior_weight[fit$prior_weights$psd == 0]
  if (length(value) == 1L) value else 0
}
cumulative_calls <- function(lfdr, alpha) {
  ordering <- order(lfdr, method = "radix")
  sum(cumsum(lfdr[ordering]) / seq_along(ordering) <= alpha)
}
summary_table <- data.frame(
  stage = c("raw", "bf"),
  bf_available = c(TRUE, bf_available),
  order = order,
  penalty = penalty,
  n_units = length(datasets),
  pi0 = c(null_weight(fash_fit_corr),
          if (bf_available) null_weight(fash_fit_corr_update) else NA_real_),
  mean_lfdr = c(mean(fash_fit_corr$lfdr),
                if (bf_available) {
                  mean(fash_fit_corr_update$lfdr)
                } else {
                  NA_real_
                }),
  calls_at_alpha = c(cumulative_calls(fash_fit_corr$lfdr, alpha),
                     if (bf_available) {
                       cumulative_calls(fash_fit_corr_update$lfdr, alpha)
                     } else {
                       NA_integer_
                     }),
  alpha = alpha,
  retained_components = c(nrow(fash_fit_corr$prior_weights),
                          if (bf_available) {
                            nrow(fash_fit_corr_update$prior_weights)
                          } else {
                            NA_integer_
                          }),
  correlation_source = basename(correlation_path),
  smoke_run = smoke > 0L,
  fit_minutes = as.numeric(elapsed),
  stringsAsFactors = FALSE
)
utils::write.csv(summary_table,
                 file.path(output_directory, "correlated_fit_summary.csv"),
                 row.names = FALSE)
print(summary_table, digits = 6)
log_message("Done.")
