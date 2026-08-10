#!/usr/bin/env Rscript

# Deterministic tests for donor-level null-permutation pilot helpers.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "r4_correlated_errors",
  "real_data_correlation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "covariance_estimation",
  "donor_null_permutation_helpers.R"
))

temporary_vcf <- tempfile(fileext = ".vcf.gz")
connection <- gzfile(temporary_vcf, open = "wt")
writeLines(c(
  "##fileformat=VCFv4.2",
  paste(
    c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT", "donor_b", "donor_a", "donor_c"),
    collapse = "\t"
  ),
  paste(c("1", "10", "rs_keep_2", "A", "G", ".", "PASS", ".", "DS",
          "0.2", "1.2", "2.0"), collapse = "\t"),
  paste(c("1", "20", "rs_drop", "A", "C", ".", "PASS", ".", "GT:DS",
          "0/0:0", "0/1:1", "1/1:2"), collapse = "\t"),
  paste(c("1", "30", "rs_keep_1", "T", "C", ".", "PASS", ".", "GT:DS",
          "0/1:0.7", "0/0:0.1", "1/1:1.8"), collapse = "\t")
), connection)
close(connection)

dosage <- read_selected_vcf_dosages(
  temporary_vcf,
  c("rs_keep_1", "rs_keep_2"),
  chunk_size = 100L
)
expected_dosage <- matrix(
  c(0.7, 0.1, 1.8, 0.2, 1.2, 2.0),
  nrow = 3L,
  dimnames = list(
    c("donor_b", "donor_a", "donor_c"),
    c("rs_keep_1", "rs_keep_2")
  )
)
if (!isTRUE(all.equal(dosage, expected_dosage, tolerance = 0))) {
  stop("The selected-VCF dosage extractor did not preserve IDs or values.")
}

missing_error <- try(
  read_selected_vcf_dosages(
    temporary_vcf,
    c("rs_keep_1", "rs_missing"),
    chunk_size = 100L
  ),
  silent = TRUE
)
if (!inherits(missing_error, "try-error")) {
  stop("The selected-VCF dosage extractor did not reject a missing ID.")
}

set.seed(20260806)
n <- 24L
n_units <- 5L
covariates <- cbind(
  PC1 = stats::rnorm(n),
  PC2 = stats::rnorm(n)
)
genotype <- matrix(stats::rnorm(n * n_units), nrow = n)
alpha <- matrix(stats::rnorm(3L * n_units), nrow = 3L)
beta_truth <- seq(-0.4, 0.4, length.out = n_units)
design <- cbind(1, covariates)
expression <- design %*% alpha +
  sweep(genotype, 2L, beta_truth, `*`) +
  matrix(stats::rnorm(n * n_units, sd = 0.7), nrow = n)

vectorized_fit <- fit_many_genotype_regressions(
  expression,
  genotype,
  covariates
)
for (unit_index in seq_len(n_units)) {
  data <- data.frame(
    y = expression[, unit_index],
    g = genotype[, unit_index],
    PC1 = covariates[, 1L],
    PC2 = covariates[, 2L]
  )
  reference <- stats::lm(y ~ g + PC1 + PC2, data = data)
  reference_summary <- summary(reference)$coefficients
  if (!isTRUE(all.equal(
    vectorized_fit$beta[unit_index],
    unname(reference_summary["g", "Estimate"]),
    tolerance = 1e-10
  )) || !isTRUE(all.equal(
    vectorized_fit$standard_error[unit_index],
    unname(reference_summary["g", "Std. Error"]),
    tolerance = 1e-10
  ))) {
    stop("The vectorized genotype regression disagrees with lm().")
  }
}
if (vectorized_fit$residual_df != n - 4L) {
  stop("The vectorized genotype regression returned the wrong residual df.")
}

constant_genotype <- genotype
constant_genotype[, 1L] <- 1
constant_error <- try(
  fit_many_genotype_regressions(expression, constant_genotype, covariates),
  silent = TRUE
)
rank_error <- try(
  fit_many_genotype_regressions(
    expression,
    genotype,
    cbind(covariates, duplicate_PC1 = covariates[, 1L])
  ),
  silent = TRUE
)
if (!inherits(constant_error, "try-error") ||
    !inherits(rank_error, "try-error")) {
  stop("The vectorized regression did not reject an invalid design.")
}

set.seed(202608061)
donor_ids <- paste0("donor_", 1:8)
observation_patterns <- setNames(
  c("1111", "1111", "1111", "1101", "1101", "1011", "1011", "1001"),
  donor_ids
)
source_donor <- make_shared_donor_block_permutation(
  donor_ids,
  observation_patterns
)
if (!setequal(names(source_donor), donor_ids) ||
    anyDuplicated(source_donor) ||
    any(observation_patterns[source_donor] != observation_patterns)) {
  stop("The shared donor-block permutation did not preserve missingness strata.")
}

set.seed(20260807)
n_mc <- 100000L
true_correlation <- 0.35
standard_normal_1 <- stats::rnorm(n_mc)
standard_normal_2 <- true_correlation * standard_normal_1 +
  sqrt(1 - true_correlation^2) * stats::rnorm(n_mc)
se_1 <- exp(stats::runif(n_mc, log(0.4), log(1.5)))
se_2 <- exp(stats::runif(n_mc, log(0.4), log(1.5)))
constant_effect <- stats::rnorm(n_mc, sd = 0.8)
beta_hat <- cbind(
  constant_effect + se_1 * standard_normal_1,
  constant_effect + se_2 * standard_normal_2
)
se <- cbind(se_1, se_2)
ordinary <- estimate_ordinary_pairwise_correlation(beta_hat, se)[1L, 2L]
ols <- estimate_pairwise_difference_correlation(beta_hat, se)[1L, 2L]
if (abs(ordinary - true_correlation) > 0.02 ||
    abs(ols - true_correlation) > 0.02) {
  stop("The pairwise estimators did not recover the known correlation.")
}

unequal_beta <- cbind(c(rep(sqrt(1.4), 100L), 0), 0)
unequal_se <- cbind(c(rep(1, 100L), 100), 1)
unequal_ordinary <- estimate_ordinary_pairwise_correlation(
  unequal_beta,
  unequal_se
)[1L, 2L]
unequal_ols <- estimate_pairwise_difference_correlation(
  unequal_beta,
  unequal_se
)[1L, 2L]
if (abs(unequal_ols - 0.3) > 0.01 ||
    unequal_ordinary - unequal_ols < 0.4) {
  stop("The OLS estimator did not downweight an extreme unequal-SE unit.")
}

raw_se <- matrix(stats::runif(200L * 4L, 0.1, 0.6), nrow = 200L)
beta_for_adjustment <- raw_se * matrix(stats::rnorm(length(raw_se)), nrow = 200L)
adjusted_se <- convert_raw_to_t_adjusted_se(
  beta_for_adjustment,
  raw_se,
  residual_df = c(12, 12, 9, 11)
)
if (any(adjusted_se < raw_se) || any(!is.finite(adjusted_se))) {
  stop("The t-based SE adjustment is not finite and conservative.")
}

cat("Donor null-permutation helper tests passed.\n")
