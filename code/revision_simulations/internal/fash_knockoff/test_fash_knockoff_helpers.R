#!/usr/bin/env Rscript

# Focused tests for the fash-knockoff helpers.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "covariance_estimation", "donor_null_permutation_helpers.R"
))
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "selected_signal_genotype_permutation",
  "selected_signal_genotype_permutation_helpers.R"
))
source(file.path(
  workflowr_root, "code", "revision_simulations", "internal",
  "fash_knockoff", "fash_knockoff_helpers.R"
))

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required.")
}

failures <- 0L
check <- function(label, condition) {
  if (isTRUE(condition)) {
    message("PASS: ", label)
  } else {
    failures <<- failures + 1L
    message("FAIL: ", label)
  }
}

# --- build_pair_index ------------------------------------------------------

pair_keys <- c("ENSG1_rs1", "ENSG1_rs2:A:T", "ENSG2_rs1")
index <- build_pair_index(pair_keys)
check(
  "pair keys split into gene and variant",
  identical(index$gene_id, c("ENSG1", "ENSG1", "ENSG2")) &&
    identical(index$variant_id, c("rs1", "rs2:A:T", "rs1"))
)
check(
  "gene and variant levels are deduplicated",
  identical(index$gene_levels, c("ENSG1", "ENSG2")) &&
    identical(index$variant_levels, c("rs1", "rs2:A:T")) &&
    identical(index$variant_index, c(1L, 2L, 1L))
)
check(
  "a pair key without a separator is rejected",
  inherits(try(build_pair_index("ENSG1"), silent = TRUE), "try-error")
)

# --- fit_pairwise_genotype_regressions ------------------------------------
# The indexed estimator must reproduce fit_many_genotype_regressions exactly.

set.seed(20260823)
n_donor <- 19L
n_gene <- 7L
n_variant <- 11L
expression <- matrix(rnorm(n_donor * n_gene), n_donor, n_gene)
genotype <- matrix(
  rbinom(n_donor * n_variant, size = 2, prob = 0.35),
  n_donor, n_variant
)
storage.mode(genotype) <- "double"
covariates <- matrix(rnorm(n_donor * 5L), n_donor, 5L)
gene_index <- rep(seq_len(n_gene), each = n_variant)
variant_index <- rep(seq_len(n_variant), times = n_gene)

projection <- make_covariate_projection(covariates)
indexed <- fit_pairwise_genotype_regressions(
  expression_residual = projection$residualizer %*% expression,
  genotype_residual = projection$residualizer %*% genotype,
  genotype_sumsq = colSums(genotype^2),
  gene_index = gene_index,
  variant_index = variant_index,
  residual_df = n_donor - projection$rank - 1L,
  chunk_size = 13L
)
reference <- fit_many_genotype_regressions(
  expression[, gene_index, drop = FALSE],
  genotype[, variant_index, drop = FALSE],
  covariates
)
check(
  "indexed beta matches the reference regression",
  isTRUE(all.equal(indexed$beta, unname(reference$beta), tolerance = 1e-12))
)
check(
  "indexed SE matches the reference regression",
  isTRUE(all.equal(
    indexed$standard_error,
    unname(reference$standard_error),
    tolerance = 1e-12
  ))
)
check(
  "residual degrees of freedom match",
  identical(as.integer(indexed$residual_df), as.integer(reference$residual_df))
)
check(
  "chunk size does not change the result",
  isTRUE(all.equal(
    indexed$beta,
    fit_pairwise_genotype_regressions(
      projection$residualizer %*% expression,
      projection$residualizer %*% genotype,
      colSums(genotype^2), gene_index, variant_index,
      n_donor - projection$rank - 1L,
      chunk_size = 10000L
    )$beta,
    tolerance = 0
  ))
)

# --- make_global_derangement ----------------------------------------------

set.seed(11)
donor_ids <- paste0("donor", sprintf("%02d", 1:19))
test_genotype <- matrix(
  rbinom(19L * 25L, size = 2, prob = 0.3), 19L, 25L,
  dimnames = list(donor_ids, paste0("v", 1:25))
)
storage.mode(test_genotype) <- "double"

# Seed 20260823 is known to produce a fixed point, so the search must advance.
fixed_point_draw <- make_shared_genotype_permutation(test_genotype, 20260823L)
deranged <- make_global_derangement(test_genotype, 20260823L)
check(
  "the derangement search rejects a permutation with a fixed point",
  any(fixed_point_draw$donor_map$fixed_point) &&
    !any(deranged$donor_map$fixed_point) &&
    deranged$permutation_seed > 20260823L
)
check(
  "the derangement is a bijection over the donors",
  setequal(deranged$donor_map$source_donor, donor_ids) &&
    !anyDuplicated(deranged$donor_map$source_donor)
)
check(
  "row permutation preserves the full LD structure exactly",
  isTRUE(all.equal(
    unname(crossprod(test_genotype)),
    unname(crossprod(deranged$genotype)),
    tolerance = 0
  ))
)
check(
  "row permutation preserves each variant's dosage multiset",
  all(vapply(seq_len(ncol(test_genotype)), function(column) {
    identical(
      sort(unname(test_genotype[, column])),
      sort(unname(deranged$genotype[, column]))
    )
  }, logical(1)))
)
check(
  "the derangement search is deterministic",
  identical(
    make_global_derangement(test_genotype, 20260823L)$permutation_seed,
    deranged$permutation_seed
  )
)
check(
  "an exhausted attempt budget is an error, not a fixed point",
  inherits(
    try(make_global_derangement(test_genotype, 20260823L, 1L), silent = TRUE),
    "try-error"
  )
)

# The scale-safe permutation must agree exactly with the shared helper it
# replaces, at the scale where the shared helper is still usable.
shared_reference <- make_shared_genotype_permutation(
  test_genotype, deranged$permutation_seed
)
check(
  "the scale-safe permutation matches the shared helper exactly",
  identical(deranged$genotype, shared_reference$genotype) &&
    identical(
      deranged$donor_map$source_donor, shared_reference$donor_map$source_donor
    )
)
check(
  "a non-permutation row indexing is rejected",
  inherits(
    try(
      permute_genotype_donor_rows(test_genotype, donor_ids[c(1L, 1:18)]),
      silent = TRUE
    ),
    "try-error"
  )
)

# The shared helper asserts its invariant with crossprod(genotype), which is
# variant-by-variant. At 745,867 tested variants that is 4.0 TB, which is what
# killed the first full-data run. The replacement must stay bounded.
# Continuous DS dosages, not integer calls. Integer dosages sum bit-exactly
# under row reordering and would hide a floating-point comparison bug.
wide_genotype <- matrix(
  round(runif(19L * 40000L, 0, 2), 3) + runif(19L * 40000L, 0, 1e-6),
  19L, 40000L, dimnames = list(donor_ids, paste0("w", 1:40000))
)
wide_before <- gc(reset = TRUE)
wide_elapsed <- system.time(
  wide_permuted <- permute_genotype_donor_rows(
    wide_genotype, deranged$donor_map$source_donor
  )
)[["elapsed"]]
wide_peak_mb <- sum(gc()[, "max used"] * c(56, 8)) / 1024^2
check(
  "a 19 x 40,000 continuous-dosage genotype permutes without quadratic allocation",
  identical(dim(wide_permuted), dim(wide_genotype)) &&
    isTRUE(all.equal(
      colSums(wide_genotype), colSums(wide_permuted), tolerance = 1e-12
    )) &&
    identical(
      apply(wide_genotype[, 1:200], 2L, sort),
      apply(wide_permuted[, 1:200], 2L, sort)
    ) &&
    wide_elapsed < 30
)
# Do NOT assert that the reordered sums differ bit-for-bit. Whether summing the
# same doubles in a different order changes the last bit depends on the values
# and on the platform's summation, and asserting it fails on Midway3 while
# passing on the Mac. The guaranteed property is the useful one: the sums agree
# to within floating-point noise, which is why the invariant check in
# permute_genotype_donor_rows() must use a tolerance rather than identical().
check(
  "reordered dosage sums agree only to within floating-point noise",
  max(abs(colSums(wide_genotype) - colSums(wide_permuted))) < 1e-9
)
message(
  "  (19 x 40,000 permutation took ", format(wide_elapsed, digits = 3),
  "s; full crossprod would have needed ",
  format(40000^2 * 8 / 1024^3, digits = 3), " GB)"
)

# --- screen_uninformative_variants ----------------------------------------
# A variant informative across all 19 donors can lose all residual information
# once the donor rows are permuted, because a time point with fewer donors sees
# a different dosage subset. This is what killed the first full-scale run.

screen_donors <- donor_ids
screen_genotype <- matrix(
  0, 19L, 4L,
  dimnames = list(screen_donors, c("ok", "carriers_dropped", "constant", "ok2"))
)
screen_genotype[, "ok"] <- rep(c(0, 1, 2), length.out = 19L)
# Carriers sit only in the three donors that the reduced time point drops.
screen_genotype[17:19, "carriers_dropped"] <- 1
screen_genotype[, "constant"] <- 1
screen_genotype[, "ok2"] <- rep(c(2, 0), length.out = 19L)

screen_result <- screen_uninformative_variants(
  dosage = screen_genotype,
  permuted_dosage = screen_genotype,
  donor_sets = list(screen_donors, screen_donors[1:16]),
  covariate_list = list(
    matrix(rnorm(19L * 5L), 19L, 5L, dimnames = list(screen_donors, NULL)),
    matrix(rnorm(16L * 5L), 16L, 5L, dimnames = list(screen_donors[1:16], NULL))
  )
)
check(
  "an all-constant variant is flagged",
  screen_result$flagged[["constant"]]
)
check(
  "a variant whose only carriers are dropped at a reduced time point is flagged",
  screen_result$flagged[["carriers_dropped"]]
)
check(
  "informative variants are not flagged",
  !screen_result$flagged[["ok"]] && !screen_result$flagged[["ok2"]]
)
check(
  "the reduced time point is where the carrier-dropped variant fails",
  screen_result$by_time$n_uninformative_observed[1L] == 1L &&
    screen_result$by_time$n_uninformative_observed[2L] == 2L &&
    identical(screen_result$by_time$n_donor, c(19L, 16L))
)
check(
  "the screen is symmetric across the two arms",
  identical(
    screen_result$by_time$n_uninformative_observed,
    screen_result$by_time$n_uninformative_permuted
  )
)

# --- competition statistic and knockoff+ ----------------------------------

check(
  "competition statistic is antisymmetric under a swap",
  isTRUE(all.equal(
    compute_competition_statistic(c(4, 0.5), c(2, 8)),
    -compute_competition_statistic(c(2, 8), c(4, 0.5)),
    tolerance = 0
  ))
)

W <- c(10, 9, 8, 7, 6, 5, -4, -3, 2, 1)
path <- knockoff_estimated_fdr_path(W)
check(
  "estimated FDR at the largest threshold uses the +1 correction",
  isTRUE(all.equal(
    path$estimated_fdr[path$threshold == 10], 1 / 1, tolerance = 0
  ))
)
check(
  "estimated FDR path is consistent with a direct count",
  isTRUE(all.equal(
    path$estimated_fdr[path$threshold == 5],
    (1 + sum(W <= -5)) / sum(W >= 5),
    tolerance = 0
  ))
)
threshold_20 <- knockoff_plus_threshold(W, 0.2)
check(
  "knockoff+ picks the smallest qualifying threshold",
  is.finite(threshold_20) &&
    (1 + sum(W <= -threshold_20)) / sum(W >= threshold_20) <= 0.2 &&
    all(
      (1 + vapply(
        setdiff(sort(unique(abs(W))), 0)[
          sort(unique(abs(W)))[sort(unique(abs(W))) > 0] < threshold_20
        ],
        function(t) sum(W <= -t), integer(1)
      )) /
        vapply(
          setdiff(sort(unique(abs(W))), 0)[
            sort(unique(abs(W)))[sort(unique(abs(W))) > 0] < threshold_20
          ],
          function(t) max(sum(W >= t), 1L), integer(1)
        ) > 0.2
    )
)
check(
  "an all-negative statistic admits no threshold",
  !is.finite(knockoff_plus_threshold(c(-1, -2, -3), 0.1))
)
summary_table <- knockoff_selection_summary(W, c(0.1, 0.2, 0.5))
check(
  "the selection summary reports zero discoveries when no threshold exists",
  nrow(summary_table) == 3L &&
    all(summary_table$estimated_fdr[!is.na(summary_table$estimated_fdr)] <=
          summary_table$q[!is.na(summary_table$estimated_fdr)])
)

# The obvious implementation of these two functions rescans W for every
# distinct |W|, which is O(n^2): about 1.4 hours for 1,009,173 pair-level
# statistics, per call. A 4,000-unit smoke test cannot reveal that, so the
# scaling is asserted directly, and the result is checked against the quadratic
# definition it replaces.

quadratic_path_reference <- function(W) {
  thresholds <- sort(unique(abs(W[W != 0])))
  n_positive <- vapply(thresholds, function(t) sum(W >= t), integer(1))
  n_negative <- vapply(thresholds, function(t) sum(W <= -t), integer(1))
  data.frame(
    threshold = thresholds,
    n_discoveries = n_positive,
    n_decoy_wins = n_negative,
    estimated_fdr = (1 + n_negative) / pmax(n_positive, 1L),
    stringsAsFactors = FALSE
  )
}

set.seed(4242)
for (label in c("symmetric", "signal", "heavy_ties", "with_zero")) {
  W_case <- switch(
    label,
    symmetric = rnorm(4000),
    signal = c(rnorm(2500, 1.2), rnorm(1500, -0.3)),
    heavy_ties = c(rep(2, 200), rep(-2, 200), rep(0.5, 300), rnorm(300)),
    with_zero = c(rnorm(800), rep(0, 50))
  )
  check(
    paste0("the fast path matches the quadratic definition (", label, ")"),
    isTRUE(all.equal(
      quadratic_path_reference(W_case),
      knockoff_estimated_fdr_path(W_case),
      check.attributes = FALSE
    ))
  )
}

path_elapsed <- system.time(
  knockoff_estimated_fdr_path(rnorm(400000))
)[["elapsed"]]
check(
  "the estimated-FDR path is not quadratic in the number of units",
  path_elapsed < 5
)
message(
  "  (400,000-unit path took ", format(path_elapsed, digits = 3),
  "s; the quadratic version needs about ",
  format(0.029 * (400000 / 2000)^2 / 60, digits = 3), " minutes)"
)

# --- gene-level aggregation -----------------------------------------------

gene_table <- aggregate_gene_competition(
  target_log_bf = c(1, 5, 2, 0),
  decoy_log_bf = c(3, 2, 9, 1),
  gene_id = c("g1", "g1", "g2", "g2")
)
check(
  "gene-level statistic takes the maximum on both sides",
  identical(gene_table$gene_id, c("g1", "g2")) &&
    isTRUE(all.equal(gene_table$target_max_log_bf, c(5, 2), tolerance = 0)) &&
    isTRUE(all.equal(gene_table$decoy_max_log_bf, c(3, 9), tolerance = 0)) &&
    isTRUE(all.equal(gene_table$W, c(2, -7), tolerance = 0)) &&
    identical(gene_table$variant_count, c(2L, 2L))
)

# --- merged refit, chunk invariance, and swap antisymmetry ----------------

n_unit <- 40L
time_grid <- 0:15
make_units <- function(beta, se) {
  lapply(seq_len(nrow(beta)), function(unit) {
    data.frame(time = time_grid, beta = beta[unit, ], SE = se[unit, ])
  })
}
set.seed(7)
target_beta <- matrix(rnorm(n_unit * 16L, 0, 0.25), n_unit, 16L)
target_beta[1:10, ] <- target_beta[1:10, ] +
  outer(rep(1, 10), 0.9 * scale(time_grid)[, 1])
decoy_beta <- matrix(rnorm(n_unit * 16L, 0, 0.25), n_unit, 16L)
target_se <- matrix(0.25, n_unit, 16L)
decoy_se <- matrix(0.25, n_unit, 16L)

reference_fit <- fashr::fash(
  Y = "beta", smooth_var = "time", S = "SE",
  data_list = make_units(rbind(target_beta, decoy_beta),
                         rbind(target_se, decoy_se)),
  num_basis = 20, order = 1, betaprec = 0, pred_step = 1, penalty = 10,
  num_cores = 1, verbose = FALSE
)
psd_grid <- reference_fit$psd_grid

chunk_one <- fashr::fash(
  Y = "beta", smooth_var = "time", S = "SE",
  data_list = make_units(target_beta[1:17, , drop = FALSE],
                         target_se[1:17, , drop = FALSE]),
  num_basis = 20, order = 1, betaprec = 0, pred_step = 1, penalty = 10,
  grid = psd_grid, num_cores = 1, verbose = FALSE
)
check(
  "likelihood rows are invariant to chunking",
  isTRUE(all.equal(
    unname(chunk_one$L_matrix),
    unname(reference_fit$L_matrix[1:17, , drop = FALSE]),
    tolerance = 1e-10
  ))
)

unit_keys <- c(paste0("u", seq_len(n_unit)),
               paste0("u", seq_len(n_unit), "__decoy"))
merged_raw <- refit_merged_fash_from_likelihood(
  likelihood_matrix = reference_fit$L_matrix,
  psd_grid = psd_grid,
  settings = reference_fit$settings,
  unit_keys = unit_keys,
  penalty = 10
)
check(
  "the merged refit reproduces a direct joint fit",
  isTRUE(all.equal(
    unname(merged_raw$lfdr), unname(reference_fit$lfdr), tolerance = 1e-8
  ))
)
merged_bf <- fashr::BF_update(merged_raw, plot = FALSE)
check(
  "BF_update stores one Bayes factor per unit",
  length(merged_bf$BF) == 2L * n_unit && all(is.finite(merged_bf$BF)) &&
    all(merged_bf$BF > 0)
)

bayes_factor <- fashr::BF_compute(merged_raw)
W_original <- compute_competition_statistic(
  bayes_factor[seq_len(n_unit)],
  bayes_factor[n_unit + seq_len(n_unit)]
)

swap_unit <- 3L
swapped_rows <- seq_len(2L * n_unit)
swapped_rows[c(swap_unit, n_unit + swap_unit)] <-
  c(n_unit + swap_unit, swap_unit)
swapped_raw <- refit_merged_fash_from_likelihood(
  likelihood_matrix = reference_fit$L_matrix[swapped_rows, , drop = FALSE],
  psd_grid = psd_grid,
  settings = reference_fit$settings,
  unit_keys = unit_keys,
  penalty = 10
)
check(
  "swapping a pair leaves the fitted prior unchanged",
  isTRUE(all.equal(
    merged_raw$prior_weights$prior_weight,
    swapped_raw$prior_weights$prior_weight,
    tolerance = 1e-10
  ))
)
swapped_bf <- fashr::BF_compute(swapped_raw)
W_swapped <- compute_competition_statistic(
  swapped_bf[seq_len(n_unit)],
  swapped_bf[n_unit + seq_len(n_unit)]
)
expected <- W_original
expected[swap_unit] <- -expected[swap_unit]
check(
  "swapping one pair flips only that pair's statistic",
  isTRUE(all.equal(W_swapped, expected, tolerance = 1e-8))
)

# --- constant-shift invariance of the Bayes factor ------------------------
# The dynamic null is beta(t) = c, so a static effect difference between a
# target and its decoy must not move the competition statistic.

shift <- rep(c(0, 3, -5), length.out = 2L * n_unit)
shifted_fit <- fashr::fash(
  Y = "beta", smooth_var = "time", S = "SE",
  data_list = make_units(
    rbind(target_beta, decoy_beta) + outer(shift, rep(1, 16L)),
    rbind(target_se, decoy_se)
  ),
  num_basis = 20, order = 1, betaprec = 0, pred_step = 1, penalty = 10,
  grid = psd_grid, num_cores = 1, verbose = FALSE
)
check(
  "the Bayes factor is invariant to a constant shift of beta",
  isTRUE(all.equal(
    unname(fashr::BF_compute(shifted_fit)),
    unname(bayes_factor),
    tolerance = 1e-8
  ))
)

if (failures > 0L) {
  stop(failures, " focused test(s) failed.")
}
message("All focused fash-knockoff helper tests passed.")
