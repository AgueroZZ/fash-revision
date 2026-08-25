# Helpers for the full-data FASH target-decoy (fash-knockoff) analysis.
#
# The analysis reuses the immutable full-data FASH likelihood rows for every
# target unit and computes new likelihood rows only for the permuted decoy
# units. Target and decoy units are then fitted jointly so that a single
# empirical-Bayes prior scores both sides of every pair.

# ---------------------------------------------------------------------------
# Pair bookkeeping
# ---------------------------------------------------------------------------

build_pair_index <- function(pair_keys) {
  pair_keys <- as.character(pair_keys)
  if (length(pair_keys) < 1L || any(!nzchar(pair_keys)) ||
      anyDuplicated(pair_keys)) {
    stop("The FASH pair keys are missing, empty, or duplicated.")
  }
  gene_id <- sub("_.*$", "", pair_keys)
  variant_id <- sub("^[^_]+_", "", pair_keys)
  if (any(!nzchar(gene_id)) || any(!nzchar(variant_id)) ||
      !identical(paste(gene_id, variant_id, sep = "_"), pair_keys)) {
    stop("At least one pair key does not split into a gene and a variant.")
  }
  gene_levels <- unique(gene_id)
  variant_levels <- unique(variant_id)
  list(
    pair_key = pair_keys,
    gene_id = gene_id,
    variant_id = variant_id,
    gene_levels = gene_levels,
    variant_levels = variant_levels,
    gene_index = match(gene_id, gene_levels),
    variant_index = match(variant_id, variant_levels)
  )
}

# ---------------------------------------------------------------------------
# Vectorized pairwise genotype regression
# ---------------------------------------------------------------------------

# Fits Y ~ 1 + G + PC1..PC5 for every gene-variant pair at one time point.
#
# The estimator is algebraically identical to
# `fit_many_genotype_regressions()` but never materializes a
# donor-by-pair expression matrix. With one covariate projection R,
#
#   beta = <R g, R y> / ||R g||^2
#   RSS  = ||R y||^2 - <R g, R y>^2 / ||R g||^2
#
# so only the per-pair cross product depends on the pairing. Each gene's
# residual sum of squares and each variant's residual norm are computed once.
fit_pairwise_genotype_regressions <- function(expression_residual,
                                             genotype_residual,
                                             genotype_sumsq,
                                             gene_index,
                                             variant_index,
                                             residual_df,
                                             chunk_size = 200000L) {
  expression_residual <- as.matrix(expression_residual)
  genotype_residual <- as.matrix(genotype_residual)
  genotype_sumsq <- as.numeric(genotype_sumsq)
  gene_index <- as.integer(gene_index)
  variant_index <- as.integer(variant_index)
  residual_df <- as.integer(residual_df)
  chunk_size <- as.integer(chunk_size)
  n_pair <- length(gene_index)

  # Reported individually. A single combined condition with one generic message
  # cannot say which input was wrong, and on a cluster that costs a whole
  # submission to find out.
  problems <- character()
  note <- function(condition, message) {
    if (isTRUE(condition)) problems <<- c(problems, message)
  }
  note(
    nrow(expression_residual) != nrow(genotype_residual),
    sprintf(
      "donor counts differ: expression has %d rows, genotype has %d",
      nrow(expression_residual), nrow(genotype_residual)
    )
  )
  note(
    length(genotype_sumsq) != ncol(genotype_residual),
    sprintf(
      "genotype_sumsq has length %d but the genotype has %d columns",
      length(genotype_sumsq), ncol(genotype_residual)
    )
  )
  note(
    length(variant_index) != n_pair,
    sprintf(
      "gene_index has length %d but variant_index has length %d",
      n_pair, length(variant_index)
    )
  )
  note(n_pair < 1L, "there are no pairs to fit")
  note(
    any(!is.finite(expression_residual)),
    sprintf(
      "the residualized expression has %d non-finite value(s)",
      sum(!is.finite(expression_residual))
    )
  )
  note(
    any(!is.finite(genotype_residual)),
    sprintf(
      "the residualized genotype has %d non-finite value(s)",
      sum(!is.finite(genotype_residual))
    )
  )
  note(
    any(!is.finite(genotype_sumsq)) || any(genotype_sumsq <= 0),
    sprintf(
      "%d genotype column(s) have a non-positive or non-finite sum of squares",
      sum(!is.finite(genotype_sumsq) | genotype_sumsq <= 0)
    )
  )
  note(anyNA(gene_index) || anyNA(variant_index), "an index is missing")
  if (!anyNA(gene_index)) {
    note(
      min(gene_index) < 1L || max(gene_index) > ncol(expression_residual),
      sprintf(
        "gene_index spans %d..%d but the expression has %d columns",
        min(gene_index), max(gene_index), ncol(expression_residual)
      )
    )
  }
  if (!anyNA(variant_index)) {
    note(
      min(variant_index) < 1L || max(variant_index) > ncol(genotype_residual),
      sprintf(
        "variant_index spans %d..%d but the genotype has %d columns",
        min(variant_index), max(variant_index), ncol(genotype_residual)
      )
    )
  }
  note(
    length(residual_df) != 1L || is.na(residual_df) || residual_df < 2L,
    "the residual degrees of freedom are missing or below two"
  )
  note(
    length(chunk_size) != 1L || is.na(chunk_size) || chunk_size < 1L,
    "the chunk size is missing or below one"
  )
  if (length(problems) > 0L) {
    stop(
      "Invalid pairwise genotype regression inputs:\n",
      paste0("  - ", problems, collapse = "\n")
    )
  }

  denominator <- colSums(genotype_residual^2)
  tolerance <- 1e-12 * pmax(1, genotype_sumsq)
  uninformative <- !is.finite(denominator) | denominator <= tolerance
  if (any(uninformative)) {
    stop(
      sum(uninformative), " of ", length(denominator),
      " genotype column(s) carry no residual information at this time point ",
      "(first: column ", which(uninformative)[1L],
      ", residual sum of squares ",
      format(denominator[which(uninformative)[1L]], digits = 3),
      "). Screen these variants out before fitting."
    )
  }
  gene_sumsq <- colSums(expression_residual^2)
  if (any(!is.finite(gene_sumsq)) || any(gene_sumsq <= 0)) {
    bad <- which(!is.finite(gene_sumsq) | gene_sumsq <= 0)
    stop(
      length(bad), " gene(s) have no residual expression variation at this ",
      "time point (first: column ", bad[1L], ")."
    )
  }

  cross_product <- numeric(n_pair)
  starts <- seq.int(1L, n_pair, by = chunk_size)
  for (start in starts) {
    end <- min(start + chunk_size - 1L, n_pair)
    rows <- start:end
    cross_product[rows] <- colSums(
      genotype_residual[, variant_index[rows], drop = FALSE] *
        expression_residual[, gene_index[rows], drop = FALSE]
    )
  }

  pair_denominator <- denominator[variant_index]
  beta <- cross_product / pair_denominator
  residual_sumsq <- gene_sumsq[gene_index] - cross_product^2 / pair_denominator
  negative_tolerance <- -1e-8 * pmax(1, gene_sumsq[gene_index])
  if (any(!is.finite(residual_sumsq)) ||
      any(residual_sumsq < negative_tolerance)) {
    stop("The pairwise regression produced a negative residual sum of squares.")
  }
  residual_sumsq <- pmax(residual_sumsq, 0)
  raw_se <- sqrt((residual_sumsq / residual_df) / pair_denominator)
  if (any(!is.finite(beta)) || any(!is.finite(raw_se)) || any(raw_se <= 0)) {
    stop("The pairwise genotype regression produced invalid estimates.")
  }
  list(beta = beta, standard_error = raw_se, residual_df = residual_df)
}

# ---------------------------------------------------------------------------
# The single global genotype permutation
# ---------------------------------------------------------------------------

# Permutes the donor rows of a genotype matrix, at whole-genome scale.
#
# This deliberately does NOT reuse `make_shared_genotype_permutation()`. That
# helper asserts its invariant with `crossprod(genotype)`, which is a
# variant-by-variant matrix: fine for the 1,177 variants it was written for,
# but 745,867 tested variants would need 4.0 TB.
#
# The assertion is also mathematically vacuous at any scale. For a permutation
# matrix P, (PG)^T (PG) = G^T P^T P G = G^T G identically, so the only thing it
# can detect is that the row indexing was not a permutation at all. That is
# checked directly here in O(19n): the donor map must be a bijection, and both
# the column sums and the column sums of squares must be unchanged, which no
# non-permutation row indexing can preserve by accident. A bounded random block
# of columns is additionally checked against the full crossprod identity as a
# structural spot check.
permute_genotype_donor_rows <- function(genotype,
                                       source_donor,
                                       crossprod_check_columns = 1000L) {
  genotype <- as.matrix(genotype)
  storage.mode(genotype) <- "double"
  donor_ids <- rownames(genotype)
  source_donor <- as.character(source_donor)
  if (nrow(genotype) < 2L || ncol(genotype) < 1L || is.null(donor_ids) ||
      any(!nzchar(donor_ids)) || anyDuplicated(donor_ids) ||
      any(!is.finite(genotype)) || length(source_donor) != nrow(genotype) ||
      anyDuplicated(source_donor) || !setequal(source_donor, donor_ids)) {
    stop("Invalid genotype matrix or donor permutation.")
  }
  permuted <- genotype[source_donor, , drop = FALSE]
  rownames(permuted) <- donor_ids

  # Summing the same 19 doubles in a different row order is not bit-identical,
  # so these are tolerance comparisons, not `identical()`. Dosages are
  # continuous DS values, not integer genotype calls.
  if (!isTRUE(all.equal(
        colSums(genotype), colSums(permuted), tolerance = 1e-12
      )) ||
      !isTRUE(all.equal(
        colSums(genotype^2), colSums(permuted^2), tolerance = 1e-12
      ))) {
    stop("The donor row permutation did not preserve per-variant dosages.")
  }

  # An exact multiset check and the full crossprod identity, both restricted to
  # a bounded block of columns. Run over all 745,867 variants the crossprod
  # would need 4.0 TB, and a per-column sort would be a 745,867-iteration loop.
  check_columns <- min(as.integer(crossprod_check_columns), ncol(genotype))
  if (check_columns >= 2L) {
    columns <- if (check_columns == ncol(genotype)) {
      seq_len(ncol(genotype))
    } else {
      sort(sample.int(ncol(genotype), check_columns))
    }
    original_block <- genotype[, columns, drop = FALSE]
    permuted_block <- permuted[, columns, drop = FALSE]
    if (!identical(
          apply(original_block, 2L, sort), apply(permuted_block, 2L, sort)
        )) {
      stop("The donor row permutation did not preserve per-variant dosages.")
    }
    if (!isTRUE(all.equal(
      unname(crossprod(original_block)),
      unname(crossprod(permuted_block)),
      tolerance = 1e-12
    ))) {
      stop("The donor row permutation did not preserve the LD structure.")
    }
  }
  permuted
}

# Draws one donor row permutation for the entire genotype matrix and requires
# it to be a derangement. A fixed point would leave one donor's decoy genotype
# equal to its own genotype, so that donor keeps contributing real signal to
# the decoy. With 19 donors a uniform draw has a fixed point roughly 63% of the
# time, so the seed is advanced deterministically until a derangement appears.
make_global_derangement <- function(genotype, seed, max_attempts = 100L) {
  seed <- as.integer(seed)
  max_attempts <- as.integer(max_attempts)
  donor_ids <- rownames(as.matrix(genotype))
  if (length(seed) != 1L || is.na(seed) || length(max_attempts) != 1L ||
      is.na(max_attempts) || max_attempts < 1L || is.null(donor_ids)) {
    stop("Invalid permutation seed, attempt budget, or donor labels.")
  }
  for (attempt in seq_len(max_attempts)) {
    candidate_seed <- seed + attempt - 1L
    set.seed(candidate_seed)
    source_donor <- sample(donor_ids, length(donor_ids), replace = FALSE)
    if (any(source_donor == donor_ids)) {
      next
    }
    return(list(
      genotype = permute_genotype_donor_rows(genotype, source_donor),
      donor_map = data.frame(
        target_donor = donor_ids,
        source_donor = source_donor,
        fixed_point = donor_ids == source_donor,
        stringsAsFactors = FALSE
      ),
      permutation_seed = candidate_seed,
      attempts = attempt
    ))
  }
  stop("Could not draw a derangement within the attempt budget.")
}

# Flags variants that carry no residual genotype information at some time
# point, under EITHER the observed or the permuted donor assignment.
#
# A variant can be informative in the donors present at a time point yet
# uninformative once the donor rows are permuted, because a time point with
# fewer than 19 donors sees a different subset of dosages after permutation.
# The original analysis used `lm()`, which silently aliases such a column
# rather than failing, so the cached target beta/SE may already be meaningless
# for those pairs.
#
# The screen is deliberately symmetric: a variant is dropped if it fails on
# either arm. Dropping only on the decoy arm would be a decoy-dependent filter
# and would break the exchangeability the competition relies on. The criterion
# is a property of the design, not of the expression, so it does not condition
# on any outcome.
screen_uninformative_variants <- function(dosage,
                                         permuted_dosage,
                                         donor_sets,
                                         covariate_list) {
  dosage <- as.matrix(dosage)
  permuted_dosage <- as.matrix(permuted_dosage)
  if (!identical(dim(dosage), dim(permuted_dosage)) ||
      length(donor_sets) != length(covariate_list) ||
      length(donor_sets) < 1L) {
    stop("Invalid inputs to the uninformative-variant screen.")
  }
  flagged <- logical(ncol(dosage))
  per_time <- vector("list", length(donor_sets))
  for (position in seq_along(donor_sets)) {
    donors <- donor_sets[[position]]
    projection <- make_covariate_projection(covariate_list[[position]])
    counts <- integer(2L)
    for (arm in seq_len(2L)) {
      block <- if (arm == 1L) {
        dosage[donors, , drop = FALSE]
      } else {
        permuted_dosage[donors, , drop = FALSE]
      }
      denominator <- colSums((projection$residualizer %*% block)^2)
      tolerance <- 1e-12 * pmax(1, colSums(block^2))
      bad <- !is.finite(denominator) | denominator <= tolerance
      counts[arm] <- sum(bad)
      flagged <- flagged | bad
    }
    per_time[[position]] <- data.frame(
      time_position = position,
      n_donor = length(donors),
      n_uninformative_observed = counts[1L],
      n_uninformative_permuted = counts[2L],
      stringsAsFactors = FALSE
    )
  }
  list(flagged = flagged, by_time = do.call(rbind, per_time))
}

# ---------------------------------------------------------------------------
# Target-decoy competition statistics
# ---------------------------------------------------------------------------

# W_j = log BF_j - log BF_jpi. Under the merged fit both Bayes factors are
# scored by the same prior, so swapping a pair flips the sign of W_j and
# leaves every other unit's statistic unchanged.
compute_competition_statistic <- function(target_bayes_factor,
                                          decoy_bayes_factor) {
  target_bayes_factor <- as.numeric(target_bayes_factor)
  decoy_bayes_factor <- as.numeric(decoy_bayes_factor)
  if (length(target_bayes_factor) != length(decoy_bayes_factor) ||
      length(target_bayes_factor) < 1L ||
      any(!is.finite(target_bayes_factor)) ||
      any(!is.finite(decoy_bayes_factor)) ||
      any(target_bayes_factor <= 0) || any(decoy_bayes_factor <= 0)) {
    stop("Target and decoy Bayes factors must be positive and finite.")
  }
  log(target_bayes_factor) - log(decoy_bayes_factor)
}

# Aggregates pair-level evidence to the gene level by taking the maximum on
# both sides. Both maxima run over the same variant set, so the
# maximum-over-m versus single-draw asymmetry cancels.
aggregate_gene_competition <- function(target_log_bf,
                                       decoy_log_bf,
                                       gene_id) {
  target_log_bf <- as.numeric(target_log_bf)
  decoy_log_bf <- as.numeric(decoy_log_bf)
  gene_id <- as.character(gene_id)
  if (length(target_log_bf) != length(decoy_log_bf) ||
      length(gene_id) != length(target_log_bf) ||
      length(gene_id) < 1L || any(!nzchar(gene_id)) ||
      any(!is.finite(target_log_bf)) || any(!is.finite(decoy_log_bf))) {
    stop("Invalid gene-level aggregation inputs.")
  }
  gene_levels <- unique(gene_id)
  gene_factor <- factor(gene_id, levels = gene_levels)
  target_max <- as.numeric(tapply(target_log_bf, gene_factor, max))
  decoy_max <- as.numeric(tapply(decoy_log_bf, gene_factor, max))
  variant_count <- as.integer(table(gene_factor))
  if (anyNA(target_max) || anyNA(decoy_max) || any(variant_count < 1L)) {
    stop("At least one gene has no tested variant on one side of the pairing.")
  }
  data.frame(
    gene_id = gene_levels,
    variant_count = variant_count,
    target_max_log_bf = target_max,
    decoy_max_log_bf = decoy_max,
    W = target_max - decoy_max,
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# knockoff+ selection
# ---------------------------------------------------------------------------

# The estimated FDR at threshold t is (1 + #{W <= -t}) / (#{W >= t} v 1).
#
# Both tail counts are evaluated by binary search against the sorted positive
# and negative magnitudes, so the whole path costs O(n log n). The obvious
# implementation -- loop over each distinct |W| and rescan the vector -- is
# O(n^2): about 1.4 hours for the 1,009,173 pair-level statistics, per call.
knockoff_estimated_fdr_path <- function(W) {
  W <- as.numeric(W)
  if (length(W) < 1L || any(!is.finite(W))) {
    stop("Invalid competition statistics.")
  }
  thresholds <- sort(unique(abs(W[W != 0])))
  if (length(thresholds) == 0L) {
    stop("Every competition statistic is exactly zero.")
  }
  positive_magnitudes <- sort(W[W > 0])
  negative_magnitudes <- sort(-W[W < 0])
  n_positive <- length(positive_magnitudes) -
    findInterval(thresholds, positive_magnitudes, left.open = TRUE)
  n_negative <- length(negative_magnitudes) -
    findInterval(thresholds, negative_magnitudes, left.open = TRUE)
  data.frame(
    threshold = thresholds,
    n_discoveries = as.integer(n_positive),
    n_decoy_wins = as.integer(n_negative),
    estimated_fdr = (1 + n_negative) / pmax(n_positive, 1L),
    stringsAsFactors = FALSE
  )
}

# The knockoff+ threshold is the smallest t at which that ratio is at most q.
knockoff_plus_threshold <- function(W, q, path = NULL) {
  q <- as.numeric(q)
  if (length(q) != 1L || !is.finite(q) || q <= 0 || q >= 1) {
    stop("Invalid target FDR level.")
  }
  if (is.null(path)) {
    path <- knockoff_estimated_fdr_path(W)
  }
  qualifying <- which(path$estimated_fdr <= q)
  if (length(qualifying) == 0L) {
    return(Inf)
  }
  path$threshold[qualifying[1L]]
}

knockoff_selection_summary <- function(W, q_grid) {
  W <- as.numeric(W)
  q_grid <- as.numeric(q_grid)
  if (length(q_grid) < 1L || any(!is.finite(q_grid)) ||
      any(q_grid <= 0) || any(q_grid >= 1)) {
    stop("Invalid target FDR grid.")
  }
  # One path, reused across the whole q grid.
  path <- knockoff_estimated_fdr_path(W)
  do.call(rbind, lapply(q_grid, function(q) {
    qualifying <- which(path$estimated_fdr <= q)
    if (length(qualifying) == 0L) {
      return(data.frame(
        q = q,
        threshold = NA_real_,
        n_discoveries = 0L,
        n_decoy_wins = NA_integer_,
        estimated_fdr = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    row <- path[qualifying[1L], , drop = FALSE]
    data.frame(
      q = q,
      threshold = row$threshold,
      n_discoveries = row$n_discoveries,
      n_decoy_wins = row$n_decoy_wins,
      estimated_fdr = row$estimated_fdr,
      stringsAsFactors = FALSE
    )
  }))
}

# ---------------------------------------------------------------------------
# Merged empirical-Bayes refit from a stacked likelihood matrix
# ---------------------------------------------------------------------------

# Refits the mixture on a stacked target-plus-decoy likelihood matrix without
# carrying any beta/SE data. The diagonal-likelihood FASH object needs only
# `L_matrix`, `psd_grid`, `prior_weights`, and `posterior_weights` for
# `BF_compute()` and `BF_update()`, so the one-million-unit `data_list` is
# deliberately omitted.
refit_merged_fash_from_likelihood <- function(likelihood_matrix,
                                             psd_grid,
                                             settings,
                                             unit_keys,
                                             penalty) {
  if (!requireNamespace("fashr", quietly = TRUE)) {
    stop("The fashr package is required.")
  }
  likelihood_matrix <- as.matrix(likelihood_matrix)
  psd_grid <- as.numeric(psd_grid)
  unit_keys <- as.character(unit_keys)
  if (nrow(likelihood_matrix) != length(unit_keys) ||
      ncol(likelihood_matrix) != length(psd_grid) ||
      length(unit_keys) < 2L || any(!nzchar(unit_keys)) ||
      anyDuplicated(unit_keys) || anyNA(likelihood_matrix) ||
      any(is.nan(likelihood_matrix)) || any(likelihood_matrix == Inf) ||
      sum(psd_grid == 0) != 1L) {
    stop("Invalid merged FASH likelihood inputs.")
  }
  rownames(likelihood_matrix) <- unit_keys
  empirical_bayes <- fashr::fash_eb_est(
    L_matrix = likelihood_matrix,
    grid = psd_grid,
    penalty = penalty
  )
  rownames(empirical_bayes$posterior_weight) <- unit_keys
  null_column <- which(empirical_bayes$prior_weight$psd == 0)
  if (length(null_column) != 1L) {
    stop("The merged FASH fit does not have exactly one exact-null component.")
  }
  lfdr <- empirical_bayes$posterior_weight[, null_column]
  names(lfdr) <- unit_keys
  structure(
    list(
      prior_weights = empirical_bayes$prior_weight,
      posterior_weights = empirical_bayes$posterior_weight,
      psd_grid = psd_grid,
      lfdr = lfdr,
      settings = settings,
      fash_data = list(data_list = NULL, S = NULL, Omega = NULL),
      L_matrix = likelihood_matrix,
      eb_result = empirical_bayes
    ),
    class = "fash"
  )
}

extract_prior_pi0 <- function(fash_fit) {
  null_row <- which(fash_fit$prior_weights$psd == 0)
  if (length(null_row) != 1L) {
    stop("The FASH fit does not contain exactly one null prior weight.")
  }
  as.numeric(fash_fit$prior_weights$prior_weight[null_row])
}

# ---------------------------------------------------------------------------
# Diagnostics
# ---------------------------------------------------------------------------

# Pre-flight gate for the "permute genotype only" decision: if the original
# genotype is systematically better explained by the expression PCs than a
# permuted genotype is, the decoy standard errors are systematically deflated
# and the competition becomes anti-conservative.
summarize_genotype_covariate_alignment <- function(genotype,
                                                   permuted_genotype,
                                                   covariates,
                                                   chunk_size = 20000L) {
  genotype <- as.matrix(genotype)
  permuted_genotype <- as.matrix(permuted_genotype)
  if (!identical(dim(genotype), dim(permuted_genotype)) ||
      nrow(genotype) != nrow(as.matrix(covariates))) {
    stop("Genotype, permuted genotype, and covariates are not conformable.")
  }
  projection <- make_covariate_projection(covariates)
  r_squared <- function(matrix_in) {
    total <- colSums(scale(matrix_in, center = TRUE, scale = FALSE)^2)
    residual <- colSums((projection$residualizer %*% matrix_in)^2)
    keep <- total > 0
    1 - residual[keep] / total[keep]
  }
  original <- r_squared(genotype)
  permuted <- r_squared(permuted_genotype)
  probabilities <- c(0.5, 0.75, 0.9, 0.95, 0.99)
  data.frame(
    arm = c("original", "permuted"),
    n_variant = c(length(original), length(permuted)),
    mean_r_squared = c(mean(original), mean(permuted)),
    rbind(
      stats::quantile(original, probabilities),
      stats::quantile(permuted, probabilities)
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}
