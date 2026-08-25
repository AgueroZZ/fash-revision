#!/usr/bin/env Rscript

# Estimate null beta-hat correlation matrices using gene-level minimum-lfdr screens.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath(".", mustWork = TRUE))
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return(normalizePath("coderepo-local", mustWork = TRUE))
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  equals_prefix <- paste0(name, "=")
  equals_hit <- which(startsWith(args, equals_prefix))
  if (length(equals_hit) > 0L) {
    return(substring(args[equals_hit[1L]], nchar(equals_prefix) + 1L))
  }
  hit <- which(args == name)
  if (length(hit) == 0L || hit[1L] == length(args)) {
    return(default)
  }
  args[hit[1L] + 1L]
}

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE)
}

rbind_fill <- function(items) {
  items <- Filter(Negate(is.null), items)
  if (length(items) == 0L) {
    return(data.frame())
  }
  columns <- unique(unlist(lapply(items, names), use.names = FALSE))
  items <- lapply(items, function(item) {
    missing_columns <- setdiff(columns, names(item))
    for (column in missing_columns) {
      item[[column]] <- NA
    }
    item[, columns, drop = FALSE]
  })
  output <- do.call(rbind, items)
  rownames(output) <- NULL
  output
}

threshold_id <- function(threshold) {
  sub("\\.", "p", format(threshold, nsmall = 3L, trim = TRUE))
}

matrix_to_long <- function(correlation, selection, method, estimator,
                           matrix_role) {
  time_labels <- sub("^time_", "", colnames(correlation))
  grid <- expand.grid(
    time_a = time_labels,
    time_b = time_labels,
    stringsAsFactors = FALSE
  )
  grid$correlation <- as.vector(correlation)
  grid$selection <- selection
  grid$method <- method
  grid$estimator <- estimator
  grid$matrix_role <- matrix_role
  grid
}

summarize_draw_array <- function(draws, selection, method, estimator) {
  n_replications <- dim(draws)[3L]
  time_count <- dim(draws)[1L]
  mean_matrix <- apply(draws, c(1L, 2L), mean)
  median_matrix <- apply(draws, c(1L, 2L), stats::median)
  lower_matrix <- apply(draws, c(1L, 2L), stats::quantile,
                        probs = 0.025, names = FALSE)
  upper_matrix <- apply(draws, c(1L, 2L), stats::quantile,
                        probs = 0.975, names = FALSE)
  dimnames(mean_matrix) <- dimnames(median_matrix) <-
    dimnames(lower_matrix) <- dimnames(upper_matrix) <-
    dimnames(draws)[1:2]
  lag_draws <- vapply(seq_len(n_replications), function(index) {
    lag_average_correlation(draws[, , index])
  }, numeric(time_count - 1L))
  lag_draws <- t(lag_draws)
  diagnostics <- do.call(rbind, lapply(seq_len(n_replications), function(index) {
    summarize_raw_correlation_matrix(draws[, , index])
  }))
  list(
    mean_matrix = mean_matrix,
    median_matrix = median_matrix,
    lower_matrix = lower_matrix,
    upper_matrix = upper_matrix,
    lag_summary = data.frame(
      selection = selection,
      method = method,
      estimator = estimator,
      lag = seq_len(ncol(lag_draws)),
      mean = colMeans(lag_draws),
      median = apply(lag_draws, 2L, stats::median),
      lower = apply(lag_draws, 2L, stats::quantile,
                    probs = 0.025, names = FALSE),
      upper = apply(lag_draws, 2L, stats::quantile,
                    probs = 0.975, names = FALSE),
      sd = apply(lag_draws, 2L, stats::sd),
      n_replications = n_replications,
      stringsAsFactors = FALSE
    ),
    diagnostics = do.call(rbind, lapply(names(diagnostics), function(metric) {
      values <- diagnostics[, metric]
      data.frame(
        selection = selection,
        method = method,
        estimator = estimator,
        metric = metric,
        mean = mean(values),
        median = stats::median(values),
        lower = stats::quantile(values, 0.025, names = FALSE),
        upper = stats::quantile(values, 0.975, names = FALSE),
        sd = stats::sd(values),
        n_replications = n_replications,
        stringsAsFactors = FALSE
      )
    }))
  )
}

summarize_mc_unit_correlations <- function(beta_draws,
                                           selection,
                                           method,
                                           n_bootstrap = 500L,
                                           seed = 20260810L) {
  beta_draws <- as.array(beta_draws)
  n_units <- dim(beta_draws)[1L]
  n_time <- dim(beta_draws)[2L]
  n_replications <- dim(beta_draws)[3L]
  n_bootstrap <- as.integer(n_bootstrap)
  if (length(dim(beta_draws)) != 3L || n_units < 2L || n_time < 2L ||
      n_replications < 20L || any(!is.finite(beta_draws)) ||
      is.na(n_bootstrap) || n_bootstrap < 100L || is.na(seed)) {
    stop("Invalid Monte Carlo beta-hat draws or bootstrap settings.")
  }

  unit_correlations <- array(
    NA_real_,
    dim = c(n_units, n_time, n_time),
    dimnames = list(
      dimnames(beta_draws)[[1L]],
      dimnames(beta_draws)[[2L]],
      dimnames(beta_draws)[[2L]]
    )
  )
  for (unit_index in seq_len(n_units)) {
    correlation <- stats::cor(t(beta_draws[unit_index, , ]))
    if (any(!is.finite(correlation))) {
      stop("A unit-level Monte Carlo correlation is non-finite.")
    }
    unit_correlations[unit_index, , ] <- correlation
  }
  mean_matrix <- apply(unit_correlations, c(2L, 3L), mean)
  diag(mean_matrix) <- 1

  set.seed(seed)
  bootstrap_matrices <- array(
    NA_real_,
    dim = c(n_time, n_time, n_bootstrap),
    dimnames = list(
      dimnames(beta_draws)[[2L]],
      dimnames(beta_draws)[[2L]],
      paste0("bootstrap_", seq_len(n_bootstrap))
    )
  )
  lag_draws <- matrix(NA_real_, nrow = n_bootstrap, ncol = n_time - 1L)
  diagnostic_draws <- vector("list", n_bootstrap)
  for (bootstrap_index in seq_len(n_bootstrap)) {
    sampled_units <- sample.int(n_units, n_units, replace = TRUE)
    matrix_draw <- apply(
      unit_correlations[sampled_units, , , drop = FALSE],
      c(2L, 3L),
      mean
    )
    diag(matrix_draw) <- 1
    bootstrap_matrices[, , bootstrap_index] <- matrix_draw
    lag_draws[bootstrap_index, ] <- lag_average_correlation(matrix_draw)
    diagnostic_draws[[bootstrap_index]] <-
      summarize_raw_correlation_matrix(matrix_draw)
  }
  diagnostic_draws <- do.call(rbind, diagnostic_draws)
  lower_matrix <- apply(bootstrap_matrices, c(1L, 2L), stats::quantile,
                        probs = 0.025, names = FALSE)
  upper_matrix <- apply(bootstrap_matrices, c(1L, 2L), stats::quantile,
                        probs = 0.975, names = FALSE)
  dimnames(lower_matrix) <- dimnames(upper_matrix) <- dimnames(mean_matrix)

  list(
    mean_matrix = mean_matrix,
    lower_matrix = lower_matrix,
    upper_matrix = upper_matrix,
    unit_correlations = unit_correlations,
    bootstrap_matrices = bootstrap_matrices,
    lag_summary = data.frame(
      selection = selection,
      method = method,
      estimator = "direct_mc_unit_mean",
      lag = seq_len(ncol(lag_draws)),
      mean = lag_average_correlation(mean_matrix),
      median = apply(lag_draws, 2L, stats::median),
      lower = apply(lag_draws, 2L, stats::quantile,
                    probs = 0.025, names = FALSE),
      upper = apply(lag_draws, 2L, stats::quantile,
                    probs = 0.975, names = FALSE),
      sd = apply(lag_draws, 2L, stats::sd),
      n_replications = n_replications,
      n_bootstrap = n_bootstrap,
      stringsAsFactors = FALSE
    ),
    diagnostics = do.call(rbind, lapply(names(diagnostic_draws), function(metric) {
      values <- diagnostic_draws[, metric]
      data.frame(
        selection = selection,
        method = method,
        estimator = "direct_mc_unit_mean",
        metric = metric,
        mean = summarize_raw_correlation_matrix(mean_matrix)[[metric]],
        median = stats::median(values),
        lower = stats::quantile(values, 0.025, names = FALSE),
        upper = stats::quantile(values, 0.975, names = FALSE),
        sd = stats::sd(values),
        n_replications = n_replications,
        n_bootstrap = n_bootstrap,
        stringsAsFactors = FALSE
      )
    }))
  )
}

workflowr_root <- find_workflowr_root()
project_root <- normalizePath(file.path(workflowr_root, ".."), mustWork = TRUE)
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "r4_correlated_errors",
  "real_data_correlation_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "covariance_estimation",
  "observational_null_set_helpers.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "covariance_estimation",
  "donor_null_permutation_helpers.R"
))

n_replications <- as.integer(get_arg("--n-reps", "400"))
seed_start <- as.integer(get_arg("--seed-start", "20260810"))
output_id <- get_arg("--output-id", "multigene_null_beta_covariance")
if (is.na(n_replications) || n_replications < 20L || is.na(seed_start) ||
    !nzchar(output_id)) {
  stop("Invalid replication count, seed, or output ID.")
}

fit_path <- file.path(
  workflowr_root, "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
)
expression_path <- file.path(
  project_root, "iPSC-data", "expression-data",
  "quantile_normalized_no_projection.txt"
)
pc_path <- file.path(
  workflowr_root, "data", "dynamic_eQTL_real", "principal_components_10.txt"
)
vcf_path <- file.path(
  project_root, "iPSC-data", "genotype-data", "YRI_genotype.vcf.gz"
)
required_paths <- c(fit_path, expression_path, pc_path, vcf_path)
if (any(!file.exists(required_paths))) {
  stop("Required input is missing: ",
       paste(required_paths[!file.exists(required_paths)], collapse = ", "))
}

output_dir <- file.path(
  workflowr_root, "output", "revision_simulations", "internal", output_id
)
input_dir <- file.path(output_dir, "input")
summary_dir <- file.path(output_dir, "summary")
figure_dir <- file.path(output_dir, "figures")
invisible(lapply(c(output_dir, input_dir, summary_dir, figure_dir), dir.create,
                 recursive = TRUE, showWarnings = FALSE))

message("Loading the BF-adjusted FASH fit and constructing gene-level screens.")
fit_environment <- new.env(parent = emptyenv())
loaded_names <- load(fit_path, envir = fit_environment)
if (!identical(loaded_names, "fash_fit1_update")) {
  stop("The fitted-object file must contain only fash_fit1_update.")
}
fash_fit <- fit_environment$fash_fit1_update
pair_keys <- names(fash_fit$fash_data$data_list)
lfdr <- as.numeric(fash_fit$lfdr)
if (length(pair_keys) != length(lfdr) || anyDuplicated(pair_keys) ||
    any(!is.finite(lfdr)) || any(lfdr < 0 | lfdr > 1)) {
  stop("The FASH fit does not contain valid pair keys and lfd r values.")
}

gene_id <- parse_gene_ids(pair_keys)
gene_min_lfdr <- tapply(lfdr, gene_id, min)
thresholds <- c(0.90, 0.925)
catalog <- make_randomized_pair_catalog(pair_keys, seed = 20260809L)
selection_list <- list()
selection_counts <- list()
for (threshold in thresholds) {
  selection_name <- paste0("m_g_gt_", threshold_id(threshold))
  eligible_genes <- names(gene_min_lfdr)[gene_min_lfdr > threshold]
  selection <- select_random_unique_variant_per_gene(
    catalog,
    eligible_genes = eligible_genes,
    n_select = length(eligible_genes),
    set_id = selection_name
  )
  selection$selected$gene_min_lfdr <- unname(
    gene_min_lfdr[selection$selected$gene_id]
  )
  if (any(selection$selected$gene_min_lfdr <= threshold) ||
      anyDuplicated(selection$selected$gene_id) ||
      anyDuplicated(selection$selected$variant_id)) {
    stop("The gene-level selection failed its invariants.")
  }
  selection_list[[selection_name]] <- selection
  selection_counts[[selection_name]] <- data.frame(
    selection = selection_name,
    threshold = threshold,
    n_genes_passing_min_lfdr = length(eligible_genes),
    n_genes_with_unique_variant = selection$n_eligible_genes_with_unique_variant,
    n_selected = selection$n_selected,
    stringsAsFactors = FALSE
  )
}
selection_counts <- do.call(rbind, selection_counts)
selection_090 <- selection_list[["m_g_gt_0p900"]]$selected
selection_0925 <- selection_list[["m_g_gt_0p925"]]$selected
if (!all(selection_0925$gene_id %in% selection_090$gene_id) ||
    !all(selection_0925$pair_key %in% selection_090$pair_key)) {
  stop("The strict m_g > 0.925 selection is not nested in m_g > 0.90.")
}
unit_table <- selection_090
unit_table$unit_index <- seq_len(nrow(unit_table))
selection_indices <- lapply(selection_list, function(selection) {
  index <- match(selection$selected$pair_key, unit_table$pair_key)
  if (anyNA(index)) {
    stop("A selected unit is missing from the m_g > 0.90 base panel.")
  }
  index
})

message("Reading expression and PC inputs for the selected genes.")
expression_data <- utils::read.csv(
  expression_path,
  sep = "",
  check.names = FALSE,
  stringsAsFactors = FALSE
)
pc_data <- utils::read.delim(
  pc_path,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
if (!identical(names(expression_data)[1L], "Gene_id") ||
    anyDuplicated(expression_data$Gene_id) ||
    !all(c("Sample_id", paste0("PC", 1:5)) %in% names(pc_data)) ||
    anyDuplicated(pc_data$Sample_id)) {
  stop("The expression or PC input has invalid identifiers.")
}
expression_sample_ids <- names(expression_data)[-1L]
if (length(expression_sample_ids) != 297L || nrow(pc_data) != 297L ||
    !setequal(expression_sample_ids, pc_data$Sample_id)) {
  stop("Expression and PC sample IDs are not the expected 297 matched samples.")
}
gene_row <- match(unit_table$gene_id, expression_data$Gene_id)
if (anyNA(gene_row)) {
  stop("At least one selected gene is missing from the expression matrix.")
}

message("Streaming selected genotype dosages from the VCF.")
unique_variant_ids <- unit_table$variant_id
if (anyDuplicated(unique_variant_ids)) {
  stop("The selected panel must have one globally unique variant per gene.")
}
unit_dosage <- read_selected_vcf_dosages(
  vcf_path,
  unique_variant_ids,
  chunk_size = 100000L
)
colnames(unit_dosage) <- unit_table$pair_key
vcf_donors <- rownames(unit_dosage)
if (nrow(unit_dosage) != 19L || anyDuplicated(vcf_donors)) {
  stop("The selected genotype matrix does not contain 19 unique donors.")
}

time_grid <- 0:15
expected_sample_counts <- c(
  19L, 19L, 16L, 19L, 16L, 19L, 19L, 19L,
  19L, 19L, 19L, 19L, 19L, 18L, 19L, 19L
)
time_inputs <- vector("list", length(time_grid))
names(time_inputs) <- paste0("time_", time_grid)
for (time_index in seq_along(time_grid)) {
  time_value <- time_grid[time_index]
  sample_ids <- grep(
    paste0("_", time_value, "$"),
    expression_sample_ids,
    value = TRUE
  )
  donors <- sub(paste0("_", time_value, "$"), "", sample_ids)
  if (length(sample_ids) != expected_sample_counts[time_index] ||
      anyDuplicated(donors) || any(!donors %in% vcf_donors)) {
    stop("Unexpected donor coverage at time ", time_value, ".")
  }
  pc_rows <- match(sample_ids, pc_data$Sample_id)
  expression_columns <- match(sample_ids, names(expression_data))
  expression_matrix <- t(as.matrix(
    expression_data[gene_row, expression_columns, drop = FALSE]
  ))
  storage.mode(expression_matrix) <- "double"
  rownames(expression_matrix) <- donors
  colnames(expression_matrix) <- unit_table$pair_key
  covariates <- as.matrix(pc_data[pc_rows, paste0("PC", 1:5), drop = FALSE])
  storage.mode(covariates) <- "double"
  rownames(covariates) <- donors
  projection <- make_covariate_projection(covariates)
  expression_residual <- projection$residualizer %*% expression_matrix
  rownames(expression_residual) <- donors
  colnames(expression_residual) <- unit_table$pair_key
  observed_residual_sd <- sqrt(colSums(expression_residual^2) /
    (nrow(expression_matrix) - projection$rank))
  if (any(!is.finite(observed_residual_sd)) || any(observed_residual_sd <= 0)) {
    stop("Invalid PC-residual expression scale at time ", time_value, ".")
  }
  time_inputs[[time_index]] <- list(
    time = time_value,
    donors = donors,
    expression_residual = expression_residual,
    projection = projection,
    observed_residual_sd = observed_residual_sd
  )
}

donor_observation_matrix <- vapply(time_inputs, function(input) {
  vcf_donors %in% input$donors
}, logical(length(vcf_donors)))
rownames(donor_observation_matrix) <- vcf_donors
colnames(donor_observation_matrix) <- paste0("time_", time_grid)
donor_observation_patterns <- apply(donor_observation_matrix, 1L, paste0,
                                     collapse = "")

selected_raw_data <- list(
  unit_table = unit_table,
  selection_indices = selection_indices,
  selection_counts = selection_counts,
  time_grid = time_grid,
  expected_sample_counts = expected_sample_counts,
  vcf_donors = vcf_donors,
  donor_observation_matrix = donor_observation_matrix,
  donor_observation_patterns = donor_observation_patterns,
  unit_dosage = unit_dosage,
  time_inputs = time_inputs
)
saveRDS(selected_raw_data, file.path(input_dir, "selected_raw_data.rds"),
        compress = "xz")

make_valid_donor_label_permutation <- function(max_attempts = 1000L) {
  max_attempts <- as.integer(max_attempts)
  if (is.na(max_attempts) || max_attempts < 1L) {
    stop("max_attempts must be a positive integer.")
  }
  for (attempt in seq_len(max_attempts)) {
    donor_permutation <- sample(vcf_donors, length(vcf_donors), replace = FALSE)
    names(donor_permutation) <- vcf_donors
    is_identifiable <- vapply(time_inputs, function(input) {
      permuted_donors <- unname(donor_permutation[input$donors])
      genotype <- unit_dosage[permuted_donors, , drop = FALSE]
      genotype_residual <- input$projection$residualizer %*% genotype
      denominator <- colSums(genotype_residual^2)
      tolerance <- 1e-12 * pmax(1, colSums(genotype^2))
      all(is.finite(denominator) & denominator > tolerance)
    }, logical(1))
    if (all(is_identifiable)) {
      return(donor_permutation)
    }
  }
  stop(
    "Could not draw an identifiable donor-label permutation after ",
    max_attempts,
    " attempts."
  )
}

fit_regressions <- function(generator = "observed", replication_seed = NA_integer_) {
  if (!identical(generator, "observed")) {
    set.seed(replication_seed)
  }
  beta_hat <- matrix(
    NA_real_,
    nrow = nrow(unit_table),
    ncol = length(time_grid),
    dimnames = list(unit_table$pair_key, paste0("time_", time_grid))
  )
  raw_se <- beta_hat
  if (identical(generator, "donor_residual_block_permutation")) {
    donor_permutation <- make_shared_donor_block_permutation(
      vcf_donors,
      donor_observation_patterns
    )
  } else if (identical(generator, "donor_label_permutation")) {
    donor_permutation <- make_valid_donor_label_permutation()
  } else if (!generator %in% c("observed", "time_independent_residual_permutation")) {
    stop("Unknown generator: ", generator)
  }

  for (time_index in seq_along(time_grid)) {
    input <- time_inputs[[time_index]]
    genotype <- unit_dosage[input$donors, , drop = FALSE]
    if (identical(generator, "observed") ||
        identical(generator, "donor_label_permutation")) {
      expression_residual <- input$expression_residual
    } else if (identical(generator, "donor_residual_block_permutation")) {
      source_donors <- unname(donor_permutation[input$donors])
      source_rows <- match(source_donors, input$donors)
      if (anyNA(source_rows)) {
        stop("A block-permuted donor is missing at time ", input$time, ".")
      }
      expression_residual <- input$projection$residualizer %*%
        input$expression_residual[source_rows, , drop = FALSE]
    } else {
      source_rows <- sample.int(length(input$donors), length(input$donors),
                                replace = FALSE)
      expression_residual <- input$projection$residualizer %*%
        input$expression_residual[source_rows, , drop = FALSE]
    }
    if (identical(generator, "donor_label_permutation")) {
      permuted_donors <- unname(donor_permutation[input$donors])
      genotype <- unit_dosage[permuted_donors, , drop = FALSE]
    }
    fit <- fit_residualized_genotype_regressions(
      expression_residual,
      genotype,
      input$projection$residualizer,
      input$projection$rank
    )
    beta_hat[, time_index] <- fit$beta
    raw_se[, time_index] <- fit$standard_error
  }
  list(beta_hat = beta_hat, raw_se = raw_se)
}

estimator_functions <- list(
  ordinary_mean = estimate_ordinary_pairwise_correlation,
  ols_weighted = estimate_pairwise_difference_correlation
)
null_generators <- c(
  "donor_residual_block_permutation",
  "donor_label_permutation",
  "time_independent_residual_permutation"
)

message("Fitting the observed raw eQTL regressions for reference.")
observed_fit <- fit_regressions("observed")
observed_matrices <- list()
for (selection_name in names(selection_indices)) {
  selected_index <- selection_indices[[selection_name]]
  for (estimator_name in names(estimator_functions)) {
    key <- paste(selection_name, "observed_raw", estimator_name, sep = "__")
    observed_matrices[[key]] <- estimator_functions[[estimator_name]](
      observed_fit$beta_hat[selected_index, , drop = FALSE],
      observed_fit$raw_se[selected_index, , drop = FALSE]
    )
  }
}

draws <- list()
beta_draws <- setNames(lapply(null_generators, function(generator) {
  array(
    NA_real_,
    dim = c(nrow(unit_table), length(time_grid), n_replications),
    dimnames = list(
      unit_table$pair_key,
      paste0("time_", time_grid),
      paste0("rep_", seq_len(n_replications))
    )
  )
}), null_generators)
for (selection_name in names(selection_indices)) {
  for (generator in null_generators) {
    for (estimator_name in names(estimator_functions)) {
      key <- paste(selection_name, generator, estimator_name, sep = "__")
      draws[[key]] <- array(
        NA_real_,
        dim = c(length(time_grid), length(time_grid), n_replications),
        dimnames = list(
          paste0("time_", time_grid),
          paste0("time_", time_grid),
          paste0("rep_", seq_len(n_replications))
        )
      )
    }
  }
}

message("Running ", n_replications, " data-faithful null replications per generator.")
for (generator_index in seq_along(null_generators)) {
  generator <- null_generators[generator_index]
  generator_seed_offset <- 100000L * generator_index
  for (replication in seq_len(n_replications)) {
    fit <- fit_regressions(
      generator = generator,
      replication_seed = seed_start + generator_seed_offset + replication - 1L
    )
    beta_draws[[generator]][, , replication] <- fit$beta_hat
    for (selection_name in names(selection_indices)) {
      selected_index <- selection_indices[[selection_name]]
      selected_beta <- fit$beta_hat[selected_index, , drop = FALSE]
      selected_se <- fit$raw_se[selected_index, , drop = FALSE]
      for (estimator_name in names(estimator_functions)) {
        key <- paste(selection_name, generator, estimator_name, sep = "__")
        draws[[key]][, , replication] <- estimator_functions[[estimator_name]](
          selected_beta,
          selected_se
        )
      }
    }
    if (replication %% max(1L, floor(n_replications / 10L)) == 0L) {
      message(generator, ": completed ", replication, " of ", n_replications,
              " replications.")
    }
  }
}

null_summaries <- list()
mc_summaries <- list()
matrix_long <- list()
matrix_interval_long <- list()
lag_summary <- list()
diagnostics <- list()
output_index <- 1L
for (key in names(draws)) {
  pieces <- strsplit(key, "__", fixed = TRUE)[[1L]]
  summary <- summarize_draw_array(
    draws[[key]],
    selection = pieces[1L],
    method = pieces[2L],
    estimator = pieces[3L]
  )
  null_summaries[[key]] <- summary
  matrix_long[[output_index]] <- matrix_to_long(
    summary$mean_matrix,
    selection = pieces[1L],
    method = pieces[2L],
    estimator = pieces[3L],
    matrix_role = "null_mean"
  )
  lag_summary[[output_index]] <- summary$lag_summary
  diagnostics[[output_index]] <- summary$diagnostics
  output_index <- output_index + 1L
}

message("Summarizing direct Monte Carlo covariance across null replications.")
for (generator_index in seq_along(null_generators)) {
  generator <- null_generators[generator_index]
  for (selection_index in seq_along(selection_indices)) {
    selection_name <- names(selection_indices)[selection_index]
    selected_index <- selection_indices[[selection_name]]
    summary <- summarize_mc_unit_correlations(
      beta_draws[[generator]][selected_index, , , drop = FALSE],
      selection = selection_name,
      method = generator,
      n_bootstrap = 500L,
      seed = seed_start + 10000L * generator_index + selection_index
    )
    key <- paste(selection_name, generator, "direct_mc_unit_mean", sep = "__")
    mc_summaries[[key]] <- summary
    matrix_long[[output_index]] <- matrix_to_long(
      summary$mean_matrix,
      selection = selection_name,
      method = generator,
      estimator = "direct_mc_unit_mean",
      matrix_role = "mc_null_mean"
    )
    lower_long <- matrix_to_long(
      summary$lower_matrix,
      selection = selection_name,
      method = generator,
      estimator = "direct_mc_unit_mean",
      matrix_role = "gene_bootstrap_lower"
    )
    upper_long <- matrix_to_long(
      summary$upper_matrix,
      selection = selection_name,
      method = generator,
      estimator = "direct_mc_unit_mean",
      matrix_role = "gene_bootstrap_upper"
    )
    matrix_interval_long[[length(matrix_interval_long) + 1L]] <- lower_long
    matrix_interval_long[[length(matrix_interval_long) + 1L]] <- upper_long
    lag_summary[[output_index]] <- summary$lag_summary
    diagnostics[[output_index]] <- summary$diagnostics
    output_index <- output_index + 1L
  }
}

observed_long <- do.call(rbind, lapply(names(observed_matrices), function(key) {
  pieces <- strsplit(key, "__", fixed = TRUE)[[1L]]
  matrix_to_long(
    observed_matrices[[key]],
    selection = pieces[1L],
    method = pieces[2L],
    estimator = pieces[3L],
    matrix_role = "observed"
  )
}))
observed_lags <- do.call(rbind, lapply(names(observed_matrices), function(key) {
  pieces <- strsplit(key, "__", fixed = TRUE)[[1L]]
  values <- lag_average_correlation(observed_matrices[[key]])
  data.frame(
    selection = pieces[1L],
    method = pieces[2L],
    estimator = pieces[3L],
    lag = seq_along(values),
    mean = values,
    median = values,
    lower = NA_real_,
    upper = NA_real_,
    sd = NA_real_,
    n_replications = 1L,
    stringsAsFactors = FALSE
  )
}))
observed_diagnostics <- do.call(rbind, lapply(names(observed_matrices), function(key) {
  pieces <- strsplit(key, "__", fixed = TRUE)[[1L]]
  values <- summarize_raw_correlation_matrix(observed_matrices[[key]])
  do.call(rbind, lapply(names(values), function(metric) {
    data.frame(
      selection = pieces[1L],
      method = pieces[2L],
      estimator = pieces[3L],
      metric = metric,
      mean = values[[metric]],
      median = values[[metric]],
      lower = NA_real_,
      upper = NA_real_,
      sd = NA_real_,
      n_replications = 1L,
      stringsAsFactors = FALSE
    )
  }))
}))

matrix_long <- rbind_fill(c(matrix_long, list(observed_long)))
matrix_interval_long <- rbind_fill(matrix_interval_long)
lag_summary <- rbind_fill(c(lag_summary, list(observed_lags)))
diagnostics <- rbind_fill(c(diagnostics, list(observed_diagnostics)))
rownames(matrix_long) <- rownames(matrix_interval_long) <-
  rownames(lag_summary) <- rownames(diagnostics) <- NULL

selection_labels <- c(
  m_g_gt_0p900 = expression(m[g] > 0.90),
  m_g_gt_0p925 = expression(m[g] > 0.925)
)
method_labels <- c(
  donor_residual_block_permutation = "Null: donor residual block",
  donor_label_permutation = "Null: donor label",
  time_independent_residual_permutation = "Null: time-independent residual"
)
primary_matrix_data <- matrix_long[
  matrix_long$estimator == "direct_mc_unit_mean" &
    matrix_long$matrix_role == "mc_null_mean" &
    matrix_long$method %in% names(method_labels),
  , drop = FALSE
]
primary_matrix_data$selection_label <- factor(
  primary_matrix_data$selection,
  levels = names(selection_labels),
  labels = selection_labels
)
primary_matrix_data$method_label <- factor(
  primary_matrix_data$method,
  levels = names(method_labels),
  labels = method_labels
)
primary_matrix_data$time_a <- factor(
  primary_matrix_data$time_a,
  levels = as.character(time_grid)
)
primary_matrix_data$time_b <- factor(
  primary_matrix_data$time_b,
  levels = rev(as.character(time_grid))
)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("The ggplot2 package is required for figures.")
}
heatmap_plot <- ggplot2::ggplot(
  primary_matrix_data,
  ggplot2::aes(x = time_a, y = time_b, fill = correlation)
) +
  ggplot2::geom_tile() +
  ggplot2::scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-0.1, 0.8), oob = scales::squish,
    name = "Correlation"
  ) +
  ggplot2::facet_grid(selection_label ~ method_label) +
  ggplot2::coord_equal() +
  ggplot2::labs(
    x = "Time point", y = "Time point",
    title = "Null beta-hat correlation matrices from donor permutations",
    subtitle = "Direct Monte Carlo covariance, averaged across gene-variant units"
  ) +
  ggplot2::theme_minimal(base_size = 10) +
  ggplot2::theme(
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 0, vjust = 0.5),
    strip.text = ggplot2::element_text(size = 9),
    plot.title = ggplot2::element_text(face = "bold")
  )
ggplot2::ggsave(
  file.path(figure_dir, "observed_and_null_beta_hat_correlation_heatmaps.png"),
  heatmap_plot,
  width = 13,
  height = 6.5,
  dpi = 220
)

primary_lag_data <- lag_summary[
  lag_summary$estimator == "direct_mc_unit_mean" &
    lag_summary$method %in% names(method_labels),
  , drop = FALSE
]
primary_lag_data$selection_label <- factor(
  primary_lag_data$selection,
  levels = names(selection_labels),
  labels = selection_labels
)
primary_lag_data$method_label <- factor(
  primary_lag_data$method,
  levels = names(method_labels),
  labels = method_labels
)
lag_plot <- ggplot2::ggplot(
  primary_lag_data,
  ggplot2::aes(x = lag, y = mean, color = method_label, fill = method_label)
) +
  ggplot2::geom_ribbon(
    data = primary_lag_data[primary_lag_data$n_replications > 1L, , drop = FALSE],
    ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.16,
    color = NA,
    inherit.aes = TRUE
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 1.4) +
  ggplot2::facet_wrap(~selection_label, nrow = 1L) +
  ggplot2::scale_x_continuous(breaks = c(1, 3, 5, 7, 9, 11, 13, 15)) +
  ggplot2::labs(
    x = "Time lag", y = "Mean correlation",
    color = NULL, fill = NULL,
    title = "Lag profile of null beta-hat correlation",
    subtitle = "Shaded bands are 95% intervals from resampling gene-variant units"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    legend.position = "bottom",
    plot.title = ggplot2::element_text(face = "bold")
  )
ggplot2::ggsave(
  file.path(figure_dir, "observed_and_null_beta_hat_lag_profiles.png"),
  lag_plot,
  width = 11,
  height = 5.2,
  dpi = 220
)

source_information <- data.frame(
  role = c("fash_fit", "expression", "time_specific_pcs", "genotype_vcf"),
  path = required_paths,
  size_bytes = unname(file.info(required_paths)$size),
  mtime_utc = format(file.info(required_paths)$mtime, tz = "UTC", usetz = TRUE),
  md5 = unname(tools::md5sum(required_paths)),
  stringsAsFactors = FALSE
)
configuration <- list(
  analysis_id = output_id,
  created_at = format(Sys.time(), tz = "America/Chicago", usetz = TRUE),
  n_replications = n_replications,
  seed_start = seed_start,
  time_grid = time_grid,
  selection = paste(
    "Genes satisfy a strict minimum BF-adjusted FASH lfdr threshold; one",
    "globally gene-unique variant is picked per eligible gene using a",
    "fixed random score independent of lfdr magnitude."
  ),
  regression = "Time-specific Y ~ 1 + G + PC1 + ... + PC5",
  primary_null = paste(
    "A single donor-level permutation of each PC-residual expression",
    "trajectory, constrained by observed missingness pattern and followed",
    "by re-projection on the time-specific PC residualizer; G and X stay fixed."
  ),
  sensitivities = c(
    donor_label = "One genotype-donor label permutation shared across all times.",
    independent_time_residual = "Independent donor residual permutations at each time point."
  ),
  primary_estimator = paste(
    "For each unit, correlate beta-hat values across null permutations and",
    "then average the unit-level correlation matrices. Pairwise estimators",
    "using raw regression standard errors are retained as a sensitivity."
  ),
  source_information = source_information,
  r_version = R.version.string,
  fashr_version = as.character(utils::packageVersion("fashr"))
)
analysis_result <- list(
  configuration = configuration,
  gene_min_lfdr = data.frame(
    gene_id = names(gene_min_lfdr),
    min_lfdr = as.numeric(gene_min_lfdr),
    stringsAsFactors = FALSE
  ),
  selection_counts = selection_counts,
  selected_units = lapply(selection_list, `[[`, "selected"),
  observed_matrices = observed_matrices,
  null_draws = draws,
  beta_draws = beta_draws,
  null_summaries = null_summaries,
  mc_summaries = mc_summaries,
  matrix_long = matrix_long,
  matrix_interval_long = matrix_interval_long,
  lag_summary = lag_summary,
  diagnostics = diagnostics,
  source_information = source_information
)
saveRDS(configuration, file.path(output_dir, "configuration.rds"))
saveRDS(analysis_result,
        file.path(output_dir, "multigene_null_beta_covariance.rds"),
        compress = "xz")
write_csv(selection_counts, file.path(summary_dir, "selection_counts.csv"))
write_csv(do.call(rbind, lapply(selection_list, `[[`, "selected")),
          file.path(summary_dir, "selected_units.csv"))
write_csv(data.frame(
  gene_id = names(gene_min_lfdr),
  min_lfdr = as.numeric(gene_min_lfdr),
  stringsAsFactors = FALSE
), file.path(summary_dir, "gene_min_lfdr.csv"))
write_csv(source_information, file.path(summary_dir, "source_information.csv"))
write_csv(matrix_long, file.path(summary_dir, "correlation_matrices_long.csv"))
write_csv(matrix_interval_long,
          file.path(summary_dir, "direct_mc_matrix_intervals_long.csv"))
write_csv(lag_summary, file.path(summary_dir, "lag_summaries.csv"))
write_csv(diagnostics, file.path(summary_dir, "matrix_diagnostics.csv"))

message("Saved multigene null beta-hat covariance exploration to ", output_dir)
