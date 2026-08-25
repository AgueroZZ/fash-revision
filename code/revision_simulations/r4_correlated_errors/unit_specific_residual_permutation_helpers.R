# Helpers for unit-specific signal-stripped residual-permutation diagnostics.

with_preserved_random_seed <- function(seed, expression) {
  seed <- as.integer(seed)
  if (length(seed) != 1L || is.na(seed)) {
    stop("The random seed must be one finite integer.")
  }
  had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_random_seed) {
    previous_random_seed <- get(".Random.seed", envir = .GlobalEnv)
  }
  on.exit({
    if (had_random_seed) {
      assign(".Random.seed", previous_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expression)
}

select_random_r4_units <- function(pair_keys, n_units, seed) {
  pair_keys <- as.character(pair_keys)
  n_units <- as.integer(n_units)
  seed <- as.integer(seed)
  if (length(pair_keys) < 1L || any(!nzchar(pair_keys)) ||
      anyDuplicated(pair_keys) || length(n_units) != 1L || is.na(n_units) ||
      n_units < 1L || n_units > length(pair_keys) || length(seed) != 1L ||
      is.na(seed)) {
    stop("Invalid pair keys, unit count, or selection seed.")
  }
  selected_indices <- with_preserved_random_seed(seed, {
    sort(sample.int(length(pair_keys), size = n_units, replace = FALSE))
  })
  data.frame(
    selected_order = seq_along(selected_indices),
    source_fash_index = selected_indices,
    pair_key = pair_keys[selected_indices],
    gene_id = sub("_.*$", "", pair_keys[selected_indices]),
    variant_id = sub("^[^_]+_", "", pair_keys[selected_indices]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

make_r4_donor_map <- function(donor_ids, observation_patterns, seed) {
  donor_ids <- as.character(donor_ids)
  observation_patterns <- setNames(
    as.character(observation_patterns), names(observation_patterns)
  )
  seed <- as.integer(seed)
  if (length(donor_ids) < 2L || any(!nzchar(donor_ids)) ||
      anyDuplicated(donor_ids) || !setequal(names(observation_patterns), donor_ids) ||
      any(!nzchar(observation_patterns)) || length(seed) != 1L || is.na(seed)) {
    stop("Invalid donor IDs, observation patterns, or permutation seed.")
  }
  source_donor <- with_preserved_random_seed(seed, {
    make_shared_donor_block_permutation(
      donor_ids = donor_ids,
      observation_patterns = observation_patterns[donor_ids]
    )
  })
  donor_map <- data.frame(
    target_donor = donor_ids,
    source_donor = unname(source_donor[donor_ids]),
    fixed_point = donor_ids == unname(source_donor[donor_ids]),
    observation_pattern = unname(observation_patterns[donor_ids]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  source_pattern <- unname(observation_patterns[donor_map$source_donor])
  if (anyDuplicated(donor_map$source_donor) ||
      !setequal(donor_map$source_donor, donor_ids) ||
      !identical(source_pattern, donor_map$observation_pattern)) {
    stop("The residual donor map is not a valid within-pattern bijection.")
  }
  donor_map
}

summarize_r4_null_draws <- function(beta_draws, time_grid = 0:15) {
  beta_draws <- as.matrix(beta_draws)
  time_grid <- as.numeric(time_grid)
  if (nrow(beta_draws) < 3L || ncol(beta_draws) != length(time_grid) ||
      any(!is.finite(beta_draws)) || anyDuplicated(time_grid)) {
    stop("Invalid null beta draws or time grid.")
  }
  correlation <- stats::cor(beta_draws)
  if (any(!is.finite(correlation)) || max(abs(diag(correlation) - 1)) > 1e-12) {
    stop("The null-draw correlation matrix is invalid.")
  }
  lags <- seq_len(ncol(beta_draws) - 1L)
  variogram <- data.frame(
    lag = lags,
    mean_correlation = vapply(lags, function(lag) {
      mean(correlation[cbind(seq_len(ncol(beta_draws) - lag),
                             seq.int(1L + lag, ncol(beta_draws)))])
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  variogram$semivariogram <- 1 - variogram$mean_correlation
  list(correlation = correlation, variogram = variogram)
}

bootstrap_r4_null_variogram_intervals <- function(beta_draws,
                                                   n_bootstrap = 1000L,
                                                   seed,
                                                   time_grid = 0:15) {
  beta_draws <- as.matrix(beta_draws)
  n_bootstrap <- as.integer(n_bootstrap)
  seed <- as.integer(seed)
  observed <- summarize_r4_null_draws(beta_draws, time_grid)$variogram
  if (length(n_bootstrap) != 1L || is.na(n_bootstrap) || n_bootstrap < 20L ||
      length(seed) != 1L || is.na(seed)) {
    stop("Invalid bootstrap replication count or seed.")
  }
  bootstrap_values <- with_preserved_random_seed(seed, {
    replicate(n_bootstrap, {
      draw_indices <- sample.int(nrow(beta_draws), nrow(beta_draws), replace = TRUE)
      summarize_r4_null_draws(beta_draws[draw_indices, , drop = FALSE], time_grid)$
        variogram$semivariogram
    })
  })
  if (!identical(dim(bootstrap_values), c(15L, n_bootstrap)) ||
      any(!is.finite(bootstrap_values))) {
    stop("The bootstrap variogram draws are invalid.")
  }
  observed$lower <- apply(bootstrap_values, 1L, stats::quantile,
    probs = 0.025, names = FALSE, type = 8
  )
  observed$upper <- apply(bootstrap_values, 1L, stats::quantile,
    probs = 0.975, names = FALSE, type = 8
  )
  observed$n_bootstrap <- n_bootstrap
  observed
}

correlation_matrix_to_long_r4 <- function(correlation, unit_label) {
  correlation <- as.matrix(correlation)
  if (!identical(dim(correlation), c(16L, 16L)) ||
      any(!is.finite(correlation))) {
    stop("Expected a finite 16 by 16 correlation matrix.")
  }
  grid <- expand.grid(
    time_row = seq_len(nrow(correlation)) - 1L,
    time_column = seq_len(ncol(correlation)) - 1L
  )
  grid$correlation <- as.vector(correlation)
  grid$unit_label <- as.character(unit_label)
  grid
}
