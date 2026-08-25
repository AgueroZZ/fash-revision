#!/usr/bin/env Rscript

# Build retained caches for the internal current-FASH versus FASH-CL page.

find_workflowr_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "_workflowr.yml"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the workflowr project root.")
    }
    current <- parent
  }
}

file_provenance <- function(path) {
  information <- file.info(path)
  data.frame(
    path = normalizePath(path, winslash = "/", mustWork = TRUE),
    byte_size = unname(information$size),
    md5 = unname(tools::md5sum(path)),
    modified_at = format(information$mtime, tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}

load_exact_object <- function(path, expected_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  if (!identical(loaded, expected_name)) {
    stop("Unexpected object in ", path, ".")
  }
  environment[[expected_name]]
}

atomic_save_rds <- function(object, path) {
  temporary <- paste0(path, ".tmp")
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically replace ", path, ".")
  }
  invisible(path)
}

message_step <- function(index, total, text) {
  message("[", index, "/", total, "] ", text)
}

summarize_discoveries <- function(pair_keys, lfdr, alpha = 0.05) {
  indices <- select_cumulative_lfdr_calls(lfdr, alpha)
  parsed <- parse_pair_keys(pair_keys[indices])
  parsed$original_index <- indices
  parsed$lfdr <- as.numeric(lfdr[indices])
  parsed <- parsed[order(parsed$lfdr, parsed$key), , drop = FALSE]
  rownames(parsed) <- NULL
  list(
    indices = indices,
    table = parsed,
    pair_count = length(indices),
    gene_count = length(unique(parsed$gene_id)),
    variant_count = length(unique(parsed$variant_id))
  )
}

resolve_manuscript_units <- function(specification, pair_keys, gene_map) {
  pair_table <- parse_pair_keys(pair_keys)
  pair_table$original_index <- seq_len(nrow(pair_table))
  output <- vector("list", nrow(specification))
  for (row in seq_len(nrow(specification))) {
    candidate_gene_ids <- unique(gene_map$ensembl_gene_id[
      gene_map$hgnc_symbol == specification$gene_symbol[row]
    ])
    candidates <- pair_table[
      pair_table$gene_id %in% candidate_gene_ids &
        pair_table$variant_id == specification$variant_id[row],
      ,
      drop = FALSE
    ]
    if (nrow(candidates) != 1L) {
      stop(
        "Could not uniquely resolve manuscript unit ",
        specification$gene_symbol[row], " / ", specification$variant_id[row],
        "."
      )
    }
    output[[row]] <- cbind(
      specification[row, , drop = FALSE],
      candidates[, c("key", "gene_id", "original_index"), drop = FALSE]
    )
  }
  do.call(rbind, output)
}

load_current_classifications <- function(paths, pair_keys) {
  categories <- names(paths)
  output <- lapply(categories, function(category) {
    expected_name <- paste0("testing_", category, "_dyn")
    result <- load_exact_object(paths[[category]], expected_name)
    keys <- rownames(result)
    if (is.null(keys) || anyNA(match(keys, pair_keys))) {
      stop("Current classification keys are not aligned for ", category, ".")
    }
    data.frame(
      original_index = as.integer(result$indices),
      key = keys,
      category = category,
      lfsr = as.numeric(result$lfsr),
      cfsr = as.numeric(result$cfsr),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, output)
  parsed <- parse_pair_keys(result$key)
  result$gene_id <- parsed$gene_id
  result$variant_id <- parsed$variant_id
  rownames(result) <- NULL
  result
}

make_plot_bundle <- function(fit,
                             original_indices,
                             method,
                             order_label,
                             smooth_var,
                             sample_size,
                             seed) {
  original_indices <- sort(unique(as.integer(original_indices)))
  subset_fit <- subset_fash_fit(fit, original_indices)
  local_indices <- seq_along(original_indices)
  posterior <- extract_posterior_plot_data(
    subset_fit,
    local_indices,
    method = method,
    order_label = order_label,
    smooth_var = smooth_var,
    sample_size = sample_size,
    seed = seed
  )
  observed <- do.call(rbind, lapply(local_indices, function(local_index) {
    extract_observed_data(subset_fit, local_index, method, order_label)
  }))
  parametric <- fit_parametric_curves(observed, grid = smooth_var)
  list(posterior = posterior, observed = observed, parametric = parametric)
}

compute_quadratic_misfit <- function(raw_fit, original_indices) {
  vapply(original_indices, function(index) {
    data <- raw_fit$fash_data$data_list[[index]]
    standard_error <- as.numeric(raw_fit$fash_data$S[[index]])
    weights <- 1 / standard_error^2
    fit <- stats::lm(
      y ~ x + I(x^2),
      data = data,
      weights = weights
    )
    weighted_rmse <- sqrt(stats::weighted.mean(stats::residuals(fit)^2, weights))
    weighted_rmse / max(stats::median(standard_error), 1e-8)
  }, numeric(1))
}

run_checkpointed_classification <- function(fit,
                                             progress_path,
                                             smooth_var,
                                             sample_size,
                                             switch_threshold,
                                             seed,
                                             cores,
                                             batch_size,
                                             input_md5) {
  pair_keys <- names(fit$fash_data$data_list)
  configuration <- list(
    pair_keys = pair_keys,
    original_indices = fit$original_indices,
    smooth_var = smooth_var,
    sample_size = as.integer(sample_size),
    switch_threshold = switch_threshold,
    seed = as.integer(seed),
    input_md5 = input_md5
  )
  progress <- if (file.exists(progress_path)) readRDS(progress_path) else NULL
  if (is.null(progress)) {
    progress <- list(configuration = configuration, results = data.frame())
  } else if (!identical(progress$configuration, configuration)) {
    stop("The classification checkpoint configuration has changed.")
  }

  completed <- if (nrow(progress$results) == 0L) {
    integer()
  } else {
    as.integer(progress$results$original_index)
  }
  missing_local <- which(!fit$original_indices %in% completed)
  total_batches <- ceiling(length(missing_local) / batch_size)
  if (length(missing_local) == 0L) {
    message("  Classification checkpoint is already complete.")
  }

  for (batch_number in seq_len(total_batches)) {
    first <- (batch_number - 1L) * batch_size + 1L
    last <- min(batch_number * batch_size, length(missing_local))
    batch_local <- missing_local[first:last]
    batch_started <- proc.time()[["elapsed"]]
    batch_results <- parallel::mclapply(
      batch_local,
      function(local_index) {
        original_index <- fit$original_indices[local_index]
        set.seed(seed + original_index)
        samples <- predict(
          fit,
          index = local_index,
          smooth_var = smooth_var,
          only.samples = TRUE,
          M = sample_size
        )
        lfsr <- classify_functional_draws(
          samples,
          smooth_var = smooth_var,
          switch_threshold = switch_threshold
        )
        data.frame(
          original_index = original_index,
          key = pair_keys[local_index],
          early = unname(lfsr["early"]),
          middle = unname(lfsr["middle"]),
          late = unname(lfsr["late"]),
          switch = unname(lfsr["switch"]),
          stringsAsFactors = FALSE
        )
      },
      mc.cores = cores,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
    if (any(vapply(batch_results, inherits, logical(1), "try-error"))) {
      stop("At least one posterior classification worker failed.")
    }
    batch_table <- do.call(rbind, batch_results)
    progress$results <- rbind(progress$results, batch_table)
    progress$results <- progress$results[
      !duplicated(progress$results$original_index),
      ,
      drop = FALSE
    ]
    progress$results <- progress$results[
      order(progress$results$original_index),
      ,
      drop = FALSE
    ]
    atomic_save_rds(progress, progress_path)
    message(
      "  Classification batch ", batch_number, "/", total_batches,
      ": retained ", nrow(progress$results), "/", length(pair_keys),
      " units in ",
      round(proc.time()[["elapsed"]] - batch_started, 1), " seconds."
    )
  }

  if (nrow(progress$results) != length(pair_keys) ||
      any(!is.finite(as.matrix(progress$results[, c(
        "early", "middle", "late", "switch"
      )])))) {
    stop("The completed classification checkpoint is incomplete.")
  }
  progress$results
}

classification_wide_to_long <- function(wide_table) {
  categories <- c("early", "middle", "late", "switch")
  long <- do.call(rbind, lapply(categories, function(category) {
    data.frame(
      original_index = wide_table$original_index,
      key = wide_table$key,
      category = category,
      lfsr = as.numeric(wide_table[[category]]),
      stringsAsFactors = FALSE
    )
  }))
  parsed <- parse_pair_keys(long$key)
  long$gene_id <- parsed$gene_id
  long$variant_id <- parsed$variant_id
  add_cumulative_fsr(long)
}

workflowr_root <- find_workflowr_root()
analysis_directory <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "fash_cl_manuscript_impact"
)
bf_helper_path <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "fash_cl_variant_enrichment_comparison",
  "fash_cl_variant_enrichment_helpers.R"
)
helper_path <- file.path(
  analysis_directory,
  "fash_cl_manuscript_impact_helpers.R"
)
source(bf_helper_path)
source(helper_path)
suppressPackageStartupMessages(library(fashr))

fit_directory <- file.path(workflowr_root, "output", "dynamic_eQTL_real")
current_iwp1_raw_path <- file.path(fit_directory, "fash_fit1_all.RData")
current_iwp1_path <- file.path(fit_directory, "fash_fit1_update.RData")
current_iwp2_path <- file.path(fit_directory, "fash_fit2_update.RData")
cl_iwp1_raw_path <- file.path(fit_directory, "fash_fit1_all_CL.RData")
cl_iwp2_raw_path <- file.path(fit_directory, "fash_fit2_all_CL.RData")
cl_iwp1_bf_path <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fash_cl_variant_enrichment_comparison",
  "fash_cl_bf_adjustment.rds"
)
gene_map_path <- file.path(fit_directory, "cache_gene_map.rds")
strober_linear_path <- file.path(
  workflowr_root,
  "data",
  "dynamic_eQTL_real",
  "strober_linear",
  "linear_dynamic_eqtls_5_pc.txt"
)
strober_quadratic_path <- file.path(
  workflowr_root,
  "data",
  "dynamic_eQTL_real",
  "strober_nonlinear",
  "non_linear_dynamic_eqtls_5_pc.txt"
)
classification_paths <- c(
  early = file.path(fit_directory, "classify_dyn_eQTLs_early.RData"),
  middle = file.path(fit_directory, "classify_dyn_eQTLs_middle.RData"),
  late = file.path(fit_directory, "classify_dyn_eQTLs_late.RData"),
  switch = file.path(fit_directory, "classify_dyn_eQTLs_switch.RData")
)
output_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "fash_cl_manuscript_impact"
)
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
progress_path <- file.path(output_directory, "classification_progress.rds")

input_paths <- c(
  current_iwp1_raw_path,
  current_iwp1_path,
  current_iwp2_path,
  cl_iwp1_raw_path,
  cl_iwp2_raw_path,
  cl_iwp1_bf_path,
  gene_map_path,
  strober_linear_path,
  strober_quadratic_path,
  classification_paths,
  bf_helper_path,
  helper_path,
  file.path(analysis_directory, "run_fash_cl_manuscript_impact.R")
)
if (any(!file.exists(input_paths))) {
  stop("At least one required retained input is missing.")
}

analysis_started <- proc.time()[["elapsed"]]
alpha <- 0.05
smooth_var <- seq(0, 15, by = 0.1)
posterior_sample_size <- 3000L
switch_threshold <- 0.25
classification_seed <- 20260810L
classification_cores <- as.integer(Sys.getenv(
  "FASH_CL_CLASSIFICATION_CORES",
  unset = "12"
))
classification_batch_size <- 64L
if (!is.finite(classification_cores) || classification_cores < 1L) {
  stop("FASH_CL_CLASSIFICATION_CORES must be a positive integer.")
}

gene_map <- readRDS(gene_map_path)
required_gene_columns <- c("ensembl_gene_id", "hgnc_symbol")
if (!all(required_gene_columns %in% names(gene_map))) {
  stop("The retained gene map is missing required columns.")
}
strober_linear <- utils::read.delim(strober_linear_path, stringsAsFactors = FALSE)
strober_quadratic <- utils::read.delim(
  strober_quadratic_path,
  stringsAsFactors = FALSE
)
for (name in c("strober_linear", "strober_quadratic")) {
  table <- get(name)
  table$key <- paste0(table$ensamble_id, "_", table$rs_id)
  assign(name, table)
}

figure3_spec <- data.frame(
  artifact = "Figure 3",
  panel = LETTERS[1:6],
  gene_symbol = c("VAMP3", "GINM1", "ADAMTS12", "FKBP9", "TCF12", "SRXN1"),
  variant_id = c(
    "rs12130857", "rs11963075", "rs1593084", "rs3750076",
    "rs11071295", "rs911446"
  ),
  role = c(
    "Strober overlap", "Strober overlap", "Strong FASH-only",
    "Strong FASH-only", "Threshold-near FASH-only",
    "Threshold-near FASH-only"
  ),
  stringsAsFactors = FALSE
)
figure5_spec <- data.frame(
  artifact = "Figure 5",
  panel = LETTERS[1:4],
  gene_symbol = c("AMN", "RIPK2", "LPAR2", "ST6GALNAC2"),
  variant_id = c("rs34145961", "rs58780282", "rs73004962", "rs11657534"),
  role = "Nonlinear dynamic",
  stringsAsFactors = FALSE
)
figure6_spec <- data.frame(
  artifact = "Figure 6",
  panel = LETTERS[1:4],
  gene_symbol = c("PPIP5K2", "RIPK2", "SELENOP", "GRAMD4"),
  variant_id = c("rs154356", "rs58780282", "rs315268", "rs12170707"),
  role = c("early", "middle", "late", "switch"),
  stringsAsFactors = FALSE
)

message_step(1, 8, "Loading authoritative current FASH results.")
current_iwp1 <- load_exact_object(current_iwp1_path, "fash_fit1_update")
current_iwp1_adjusted_prior <- extract_full_prior_weights(current_iwp1)
current_pair_keys <- names(current_iwp1$fash_data$data_list)
current_pair_table <- parse_pair_keys(current_pair_keys)
current_pair_table$original_index <- seq_len(nrow(current_pair_table))
figure3_units <- resolve_manuscript_units(
  figure3_spec,
  current_pair_keys,
  gene_map
)
figure5_units <- resolve_manuscript_units(
  figure5_spec,
  current_pair_keys,
  gene_map
)
figure6_units <- resolve_manuscript_units(
  figure6_spec,
  current_pair_keys,
  gene_map
)
current_iwp1_discovery <- summarize_discoveries(
  current_pair_keys,
  current_iwp1$lfdr,
  alpha
)
if (current_iwp1_discovery$pair_count != 9205L ||
    current_iwp1_discovery$gene_count != 1177L) {
  stop("Current IWP1 discovery invariants changed.")
}
current_classification <- load_current_classifications(
  classification_paths,
  current_pair_keys
)
current_classification_significant <- current_classification[
  current_classification$cfsr <= alpha,
  ,
  drop = FALSE
]
current_category_counts <- do.call(rbind, lapply(
  c("early", "middle", "late", "switch"),
  function(category) {
    subset <- current_classification_significant[
      current_classification_significant$category == category,
      ,
      drop = FALSE
    ]
    data.frame(
      method = "Current FASH",
      category = category,
      pair_count = nrow(subset),
      gene_count = length(unique(subset$gene_id)),
      stringsAsFactors = FALSE
    )
  }
))
expected_current_categories <- data.frame(
  category = c("early", "middle", "late", "switch"),
  pair_count = c(124L, 24L, 20L, 984L),
  gene_count = c(8L, 5L, 12L, 250L),
  stringsAsFactors = FALSE
)
if (!identical(
  current_category_counts[, c("category", "pair_count", "gene_count")],
  expected_current_categories
)) {
  stop("Current functional-classification invariants changed.")
}
current_iwp1_plot_indices <- unique(c(
  figure3_units$original_index,
  figure6_units$original_index
))
current_iwp1_plot_bundle <- make_plot_bundle(
  current_iwp1,
  current_iwp1_plot_indices,
  method = "Current FASH",
  order_label = "IWP1",
  smooth_var = smooth_var,
  sample_size = posterior_sample_size,
  seed = classification_seed
)
current_iwp1_lfdr <- as.numeric(current_iwp1$lfdr)
rm(current_iwp1)
invisible(gc())

current_iwp1_raw <- load_exact_object(current_iwp1_raw_path, "fash_fit1")
if (!identical(names(current_iwp1_raw$fash_data$data_list), current_pair_keys)) {
  stop("Current raw and BF-adjusted IWP1 pair keys are not aligned.")
}
current_iwp1_raw_prior <- extract_full_prior_weights(current_iwp1_raw)
rm(current_iwp1_raw)
invisible(gc())

current_iwp2 <- load_exact_object(current_iwp2_path, "fash_fit2_update")
if (!identical(names(current_iwp2$fash_data$data_list), current_pair_keys)) {
  stop("Current IWP1 and IWP2 pair keys are not aligned.")
}
current_iwp2_discovery <- summarize_discoveries(
  current_pair_keys,
  current_iwp2$lfdr,
  alpha
)
if (current_iwp2_discovery$pair_count != 44L ||
    current_iwp2_discovery$gene_count != 9L) {
  stop("Current IWP2 discovery invariants changed.")
}
current_iwp2_plot_bundle <- make_plot_bundle(
  current_iwp2,
  figure5_units$original_index,
  method = "Current FASH",
  order_label = "IWP2",
  smooth_var = smooth_var,
  sample_size = posterior_sample_size,
  seed = classification_seed
)
current_iwp2_lfdr <- as.numeric(current_iwp2$lfdr)
rm(current_iwp2)
invisible(gc())

message_step(2, 8, "Loading validated FASH-CL IWP1 BF adjustment.")
cl_iwp1_adjustment <- readRDS(cl_iwp1_bf_path)
required_adjustment <- c("pair_count", "null_weight", "alternative_weights", "lfdr")
if (!all(required_adjustment %in% names(cl_iwp1_adjustment)) ||
    cl_iwp1_adjustment$pair_count != length(current_pair_keys) ||
    length(cl_iwp1_adjustment$lfdr) != length(current_pair_keys)) {
  stop("The retained FASH-CL IWP1 BF adjustment failed validation.")
}
cl_iwp1_discovery <- summarize_discoveries(
  current_pair_keys,
  cl_iwp1_adjustment$lfdr,
  alpha
)
if (cl_iwp1_discovery$pair_count != 5393L ||
    cl_iwp1_discovery$gene_count != 686L) {
  stop("FASH-CL IWP1 discovery invariants changed.")
}

current_called <- logical(length(current_pair_keys))
current_called[current_iwp1_discovery$indices] <- TRUE
cl_called <- logical(length(current_pair_keys))
cl_called[cl_iwp1_discovery$indices] <- TRUE
discovery_status <- ifelse(
  current_called & cl_called,
  "Both",
  ifelse(
    current_called,
    "Current FASH only",
    ifelse(cl_called, "FASH-CL only", "Neither")
  )
)
discovery_status <- factor(
  discovery_status,
  levels = c("Neither", "Current FASH only", "FASH-CL only", "Both")
)
lfdr_scatter_all <- data.frame(
  current_lfdr = current_iwp1_lfdr,
  cl_lfdr = as.numeric(cl_iwp1_adjustment$lfdr),
  discovery_status = discovery_status
)
current_top_indices <- select_top_pair_indices_per_gene(
  gene_id = current_pair_table$gene_id,
  variant_id = current_pair_table$variant_id,
  lfdr = current_iwp1_lfdr,
  pair_key = current_pair_table$key
)
if (length(current_top_indices) != length(unique(current_pair_table$gene_id))) {
  stop("Current-FASH top-pair selection did not return one pair per gene.")
}
lfdr_scatter_top <- data.frame(
  key = current_pair_table$key[current_top_indices],
  gene_id = current_pair_table$gene_id[current_top_indices],
  variant_id = current_pair_table$variant_id[current_top_indices],
  current_lfdr = current_iwp1_lfdr[current_top_indices],
  cl_lfdr = as.numeric(cl_iwp1_adjustment$lfdr[current_top_indices]),
  discovery_status = discovery_status[current_top_indices],
  stringsAsFactors = FALSE
)
lfdr_summary_row <- function(label, current_lfdr, cl_lfdr, gene_count) {
  data.frame(
    selection = label,
    pair_count = length(current_lfdr),
    gene_count = gene_count,
    pearson_correlation = stats::cor(current_lfdr, cl_lfdr),
    spearman_correlation = stats::cor(
      current_lfdr,
      cl_lfdr,
      method = "spearman"
    ),
    stringsAsFactors = FALSE
  )
}
lfdr_comparison_summary <- rbind(
  lfdr_summary_row(
    "All tested gene-variant pairs",
    lfdr_scatter_all$current_lfdr,
    lfdr_scatter_all$cl_lfdr,
    length(unique(current_pair_table$gene_id))
  ),
  lfdr_summary_row(
    "Current-FASH top pair per gene",
    lfdr_scatter_top$current_lfdr,
    lfdr_scatter_top$cl_lfdr,
    nrow(lfdr_scatter_top)
  )
)

message_step(3, 8, "Running or resuming FASH-CL functional classification.")
cl_iwp1_raw <- load_exact_object(cl_iwp1_raw_path, "fash_fit1")
if (!identical(names(cl_iwp1_raw$fash_data$data_list), current_pair_keys)) {
  stop("Current FASH and FASH-CL IWP1 pair keys are not aligned.")
}
cl_iwp1_raw_prior <- extract_full_prior_weights(cl_iwp1_raw)
if (length(cl_iwp1_adjustment$alternative_weights) !=
    length(cl_iwp1_raw$psd_grid) - 1L) {
  stop("FASH-CL BF-adjusted prior weights do not align with the PSD grid.")
}
cl_iwp1_adjusted_prior <- data.frame(
  psd = as.numeric(cl_iwp1_raw$psd_grid),
  prior_weight = c(
    cl_iwp1_adjustment$null_weight,
    (1 - cl_iwp1_adjustment$null_weight) *
      cl_iwp1_adjustment$alternative_weights
  ),
  stringsAsFactors = FALSE
)
add_prior_labels <- function(data, method, adjustment) {
  data$method <- method
  data$adjustment <- adjustment
  data
}
prior_weight_comparison <- rbind(
  add_prior_labels(current_iwp1_raw_prior, "Current FASH", "Raw"),
  add_prior_labels(current_iwp1_adjusted_prior, "Current FASH", "BF-adjusted"),
  add_prior_labels(cl_iwp1_raw_prior, "FASH-CL", "Raw"),
  add_prior_labels(cl_iwp1_adjusted_prior, "FASH-CL", "BF-adjusted")
)
prior_weight_summary <- do.call(rbind, lapply(
  split(
    prior_weight_comparison,
    list(
      prior_weight_comparison$method,
      prior_weight_comparison$adjustment
    ),
    drop = TRUE
  ),
  function(data) {
    data.frame(
      method = data$method[1L],
      adjustment = data$adjustment[1L],
      null_weight = data$prior_weight[data$psd == 0],
      nonzero_component_count = sum(data$prior_weight > 0),
      weight_sum = sum(data$prior_weight),
      stringsAsFactors = FALSE
    )
  }
))
rownames(prior_weight_summary) <- NULL
if (any(abs(prior_weight_summary$weight_sum - 1) > 1e-7) ||
    any(lengths(split(
      prior_weight_comparison$prior_weight,
      list(
        prior_weight_comparison$method,
        prior_weight_comparison$adjustment
      )
    )) != length(current_iwp1_raw_prior$psd))) {
  stop("IWP1 prior-weight comparison failed structural validation.")
}
cl_iwp1_classification_fit <- build_adjusted_subfit(
  cl_iwp1_raw,
  cl_iwp1_discovery$indices,
  null_weight = cl_iwp1_adjustment$null_weight,
  alternative_weights = cl_iwp1_adjustment$alternative_weights
)
classification_wide <- run_checkpointed_classification(
  cl_iwp1_classification_fit,
  progress_path = progress_path,
  smooth_var = smooth_var,
  sample_size = posterior_sample_size,
  switch_threshold = switch_threshold,
  seed = classification_seed,
  cores = classification_cores,
  batch_size = classification_batch_size,
  input_md5 = unname(tools::md5sum(cl_iwp1_raw_path))
)
cl_classification <- classification_wide_to_long(classification_wide)
cl_classification$lfdr <- cl_iwp1_adjustment$lfdr[
  cl_classification$original_index
]
cl_classification_significant <- cl_classification[
  cl_classification$cfsr <= alpha,
  ,
  drop = FALSE
]
cl_category_counts <- do.call(rbind, lapply(
  c("early", "middle", "late", "switch"),
  function(category) {
    subset <- cl_classification_significant[
      cl_classification_significant$category == category,
      ,
      drop = FALSE
    ]
    data.frame(
      method = "FASH-CL",
      category = category,
      pair_count = nrow(subset),
      gene_count = length(unique(subset$gene_id)),
      stringsAsFactors = FALSE
    )
  }
))

message_step(4, 8, "Selecting Figure 3 and Figure 6 FASH-CL examples.")
linear_lookup <- strober_linear[, c("key", "pvalue", "eFDR")]
names(linear_lookup)[2:3] <- c("p_linear", "efdr_linear")
quadratic_lookup <- strober_quadratic[, c("key", "pvalue", "eFDR")]
names(quadratic_lookup)[2:3] <- c("p_quadratic", "efdr_quadratic")
cl_dynamic_candidates <- merge(
  cl_iwp1_discovery$table,
  linear_lookup,
  by = "key",
  all.x = TRUE,
  sort = FALSE
)
cl_dynamic_candidates <- merge(
  cl_dynamic_candidates,
  quadratic_lookup,
  by = "key",
  all.x = TRUE,
  sort = FALSE
)
for (column in c("p_linear", "efdr_linear", "p_quadratic", "efdr_quadratic")) {
  cl_dynamic_candidates[[column]][is.na(cl_dynamic_candidates[[column]])] <- 1
}
cl_dynamic_candidates$strober_significant <-
  cl_dynamic_candidates$efdr_linear <= alpha |
  cl_dynamic_candidates$efdr_quadratic <= alpha

figure3_exact_keys <- figure3_units$key
ginm1_key <- figure3_units$key[figure3_units$gene_symbol == "GINM1"]
overlap_candidates <- cl_dynamic_candidates[
  cl_dynamic_candidates$strober_significant,
  ,
  drop = FALSE
]
figure3_overlap <- overlap_candidates[
  overlap_candidates$key == ginm1_key,
  ,
  drop = FALSE
]
if (nrow(figure3_overlap) != 1L) {
  stop("GINM1 must remain a unique FASH-CL/Strober overlap example.")
}
additional_overlap <- select_distinct_gene_examples(
  overlap_candidates[
    overlap_candidates$gene_id != figure3_overlap$gene_id,
    ,
    drop = FALSE
  ],
  count = 1L,
  rank_columns = "lfdr"
)
figure3_overlap <- rbind(additional_overlap, figure3_overlap)
figure3_overlap$replacement_role <- c(
  "Strober overlap replacement",
  "Strober overlap retained from paper"
)

fash_cl_only <- cl_dynamic_candidates[
  !cl_dynamic_candidates$strober_significant &
    cl_dynamic_candidates$p_linear >= 0.20 &
    cl_dynamic_candidates$p_quadratic >= 0.20,
  ,
  drop = FALSE
]
fash_cl_only$quadratic_misfit <- compute_quadratic_misfit(
  cl_iwp1_raw,
  fash_cl_only$original_index
)
strong_cutoff <- stats::quantile(fash_cl_only$lfdr, 0.20, names = FALSE)
borderline_cutoff <- stats::quantile(fash_cl_only$lfdr, 0.80, names = FALSE)
used_genes <- unique(figure3_overlap$gene_id)
strong_candidates <- fash_cl_only[
  fash_cl_only$lfdr <= strong_cutoff &
    !fash_cl_only$gene_id %in% used_genes,
  ,
  drop = FALSE
]
figure3_strong <- select_distinct_gene_examples(
  strong_candidates,
  count = 2L,
  rank_columns = c("quadratic_misfit", "lfdr"),
  decreasing = c(TRUE, FALSE)
)
figure3_strong$replacement_role <- "Strong FASH-CL-only replacement"
used_genes <- c(used_genes, figure3_strong$gene_id)
borderline_candidates <- fash_cl_only[
  fash_cl_only$lfdr >= borderline_cutoff &
    !fash_cl_only$gene_id %in% used_genes,
  ,
  drop = FALSE
]
figure3_borderline <- select_distinct_gene_examples(
  borderline_candidates,
  count = 2L,
  rank_columns = c("quadratic_misfit", "lfdr"),
  decreasing = c(TRUE, TRUE)
)
figure3_borderline$replacement_role <- "Threshold-near FASH-CL-only replacement"
figure3_proposed <- dplyr::bind_rows(
  figure3_overlap,
  figure3_strong,
  figure3_borderline
)
figure3_proposed$artifact <- "Figure 3 proposed FASH-CL"
figure3_proposed$panel <- LETTERS[seq_len(nrow(figure3_proposed))]

figure6_proposed <- list()
used_category_genes <- character()
for (category in c("early", "middle", "late", "switch")) {
  original <- figure6_units[figure6_units$role == category, , drop = FALSE]
  significant <- cl_classification_significant[
    cl_classification_significant$category == category,
    ,
    drop = FALSE
  ]
  retained <- significant[significant$key == original$key, , drop = FALSE]
  if (nrow(retained) == 1L && !retained$gene_id %in% used_category_genes) {
    selected <- retained
    selected$replacement_status <- "Retained paper unit"
  } else if (nrow(significant) > 0L) {
    eligible <- significant[
      !significant$gene_id %in% used_category_genes,
      ,
      drop = FALSE
    ]
    selected <- select_distinct_gene_examples(
      eligible,
      count = 1L,
      rank_columns = c("lfsr", "lfdr")
    )
    selected$replacement_status <- "Same-category FASH-CL replacement"
  } else {
    eligible <- cl_classification[
      cl_classification$category == category &
        !cl_classification$gene_id %in% used_category_genes,
      ,
      drop = FALSE
    ]
    selected <- select_distinct_gene_examples(
      eligible,
      count = 1L,
      rank_columns = c("lfsr", "cfsr", "lfdr")
    )
    selected$replacement_status <-
      "No significant FASH-CL unit; closest nonsignificant candidate"
  }
  selected$category_significant <- selected$cfsr <= alpha
  selected$panel <- original$panel
  figure6_proposed[[category]] <- selected
  used_category_genes <- c(used_category_genes, selected$gene_id)
}
figure6_proposed <- do.call(rbind, figure6_proposed)
rownames(figure6_proposed) <- NULL

cl_iwp1_plot_indices <- unique(c(
  figure3_units$original_index,
  figure6_units$original_index,
  figure3_proposed$original_index,
  figure6_proposed$original_index
))
cl_iwp1_plot_fit <- build_adjusted_subfit(
  cl_iwp1_raw,
  cl_iwp1_plot_indices,
  null_weight = cl_iwp1_adjustment$null_weight,
  alternative_weights = cl_iwp1_adjustment$alternative_weights
)
cl_iwp1_plot_bundle <- make_plot_bundle(
  cl_iwp1_plot_fit,
  seq_along(cl_iwp1_plot_fit$original_indices),
  method = "FASH-CL",
  order_label = "IWP1",
  smooth_var = smooth_var,
  sample_size = posterior_sample_size,
  seed = classification_seed
)
rm(cl_iwp1_classification_fit, cl_iwp1_plot_fit, cl_iwp1_raw)
invisible(gc())

message_step(5, 8, "Applying BF adjustment to FASH-CL IWP2 and selecting Figure 5.")
cl_iwp2_raw <- load_exact_object(cl_iwp2_raw_path, "fash_fit2")
if (!identical(names(cl_iwp2_raw$fash_data$data_list), current_pair_keys)) {
  stop("Current FASH and FASH-CL IWP2 pair keys are not aligned.")
}
cl_iwp2_adjustment <- compute_bf_adjusted_lfdr(
  cl_iwp2_raw$L_matrix,
  chunk_size = 50000L,
  verbose = FALSE
)
cl_iwp2_discovery <- summarize_discoveries(
  current_pair_keys,
  cl_iwp2_adjustment$lfdr,
  alpha
)
if (cl_iwp2_discovery$pair_count != 60L ||
    cl_iwp2_discovery$gene_count != 6L) {
  stop("FASH-CL IWP2 discovery invariants changed.")
}
figure5_retained <- cl_iwp2_discovery$table[
  cl_iwp2_discovery$table$key %in% figure5_units$key,
  ,
  drop = FALSE
]
figure5_new <- select_distinct_gene_examples(
  cl_iwp2_discovery$table[
    !cl_iwp2_discovery$table$gene_id %in% figure5_retained$gene_id,
    ,
    drop = FALSE
  ],
  count = 4L - nrow(figure5_retained),
  rank_columns = "lfdr"
)
figure5_new_position <- 1L
figure5_proposed <- lapply(seq_len(nrow(figure5_units)), function(row) {
  retained <- figure5_retained[
    figure5_retained$key == figure5_units$key[row],
    ,
    drop = FALSE
  ]
  if (nrow(retained) == 1L) {
    selected <- retained
    selected$replacement_status <- "Retained paper unit"
  } else {
    selected <- figure5_new[figure5_new_position, , drop = FALSE]
    figure5_new_position <<- figure5_new_position + 1L
    selected$replacement_status <- "FASH-CL nonlinear replacement"
  }
  selected$panel <- figure5_units$panel[row]
  selected
})
figure5_proposed <- dplyr::bind_rows(figure5_proposed)
cl_iwp2_plot_indices <- unique(c(
  figure5_units$original_index,
  figure5_proposed$original_index
))
cl_iwp2_plot_fit <- build_adjusted_subfit(
  cl_iwp2_raw,
  cl_iwp2_plot_indices,
  null_weight = cl_iwp2_adjustment$null_weight,
  alternative_weights = cl_iwp2_adjustment$alternative_weights
)
cl_iwp2_plot_bundle <- make_plot_bundle(
  cl_iwp2_plot_fit,
  seq_along(cl_iwp2_plot_fit$original_indices),
  method = "FASH-CL",
  order_label = "IWP2",
  smooth_var = smooth_var,
  sample_size = posterior_sample_size,
  seed = classification_seed
)
cl_iwp2_lfdr <- cl_iwp2_adjustment$lfdr
rm(cl_iwp2_plot_fit, cl_iwp2_raw)
invisible(gc())

message_step(6, 8, "Building Figure 4 overlap and Table 1 comparisons.")
strober_linear_significant <- strober_linear[
  strober_linear$eFDR <= alpha,
  ,
  drop = FALSE
]
strober_quadratic_significant <- strober_quadratic[
  strober_quadratic$eFDR <= alpha,
  ,
  drop = FALSE
]
strober_quadratic_variants <- parse_pair_keys(
  strober_quadratic_significant$key
)$variant_id
strober_linear_variants <- parse_pair_keys(
  strober_linear_significant$key
)$variant_id

build_venn_output <- function(fash_table, method) {
  display_sets <- list(
    Genes = list(
      `Strober (Quadratic)` = unique(
        strober_quadratic_significant$ensamble_id
      ),
      `Strober (Linear)` = unique(strober_linear_significant$ensamble_id),
      `FASH (IWP1)` = unique(fash_table$gene_id)
    ),
    `Gene-variant pairs` = list(
      `Strober (Quadratic)` = strober_quadratic_significant$key,
      `Strober (Linear)` = strober_linear_significant$key,
      `FASH (IWP1)` = fash_table$key
    ),
    Variants = list(
      `Strober (Quadratic)` = unique(strober_quadratic_variants),
      `Strober (Linear)` = unique(strober_linear_variants),
      `FASH (IWP1)` = unique(fash_table$variant_id)
    )
  )
  if (identical(method, "FASH-CL")) {
    names(display_sets$Genes)[3L] <- "FASH-CL (IWP1)"
    names(display_sets$`Gene-variant pairs`)[3L] <- "FASH-CL (IWP1)"
    names(display_sets$Variants)[3L] <- "FASH-CL (IWP1)"
  }
  region_sets <- list(
    Genes = list(
      strober_quadratic_significant$ensamble_id,
      strober_linear_significant$ensamble_id,
      fash_table$gene_id
    ),
    `Gene-variant pairs` = list(
      strober_quadratic_significant$key,
      strober_linear_significant$key,
      fash_table$key
    ),
    Variants = list(
      strober_quadratic_variants,
      strober_linear_variants,
      fash_table$variant_id
    )
  )
  regions <- do.call(rbind, lapply(names(region_sets), function(unit) {
    sets <- region_sets[[unit]]
    result <- three_set_venn_regions(
      sets[[1L]],
      sets[[2L]],
      sets[[3L]],
      labels = c("Strober quadratic", "Strober linear", "FASH")
    )
    result$unit <- unit
    result
  }))
  regions$method <- method
  list(sets = display_sets, regions = regions)
}

venn_tables <- list(
  `Current FASH` = current_iwp1_discovery$table,
  `FASH-CL` = cl_iwp1_discovery$table
)
venn_output <- Map(build_venn_output, venn_tables, names(venn_tables))
venn_sets <- lapply(venn_output, `[[`, "sets")
venn_regions <- do.call(rbind, lapply(venn_output, `[[`, "regions"))
rownames(venn_regions) <- NULL

category_counts <- rbind(
  data.frame(
    method = "Current FASH",
    category = "nonlinear",
    pair_count = current_iwp2_discovery$pair_count,
    gene_count = current_iwp2_discovery$gene_count,
    stringsAsFactors = FALSE
  ),
  current_category_counts,
  data.frame(
    method = "FASH-CL",
    category = "nonlinear",
    pair_count = cl_iwp2_discovery$pair_count,
    gene_count = cl_iwp2_discovery$gene_count,
    stringsAsFactors = FALSE
  ),
  cl_category_counts
)

message_step(7, 8, "Recomputing matched-version Hallmark enrichment.")
msigdb <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
msigdb$gene <- ifelse(
  is.na(msigdb$ensembl_gene) | msigdb$ensembl_gene == "",
  msigdb$db_ensembl_gene,
  msigdb$ensembl_gene
)
term_to_gene <- unique(data.frame(
  term = msigdb$gs_name,
  gene = msigdb$gene,
  stringsAsFactors = FALSE
))
all_tested_genes <- unique(current_pair_table$gene_id)
current_switch_genes <- unique(current_classification_significant$gene_id[
  current_classification_significant$category == "switch"
])
cl_switch_genes <- unique(cl_classification_significant$gene_id[
  cl_classification_significant$category == "switch"
])
enrichment_inputs <- list(
  `Current FASH|All dynamic` = unique(current_iwp1_discovery$table$gene_id),
  `Current FASH|Switch` = current_switch_genes,
  `FASH-CL|All dynamic` = unique(cl_iwp1_discovery$table$gene_id),
  `FASH-CL|Switch` = cl_switch_genes
)
hallmark_all <- do.call(rbind, lapply(names(enrichment_inputs), function(name) {
  split_name <- strsplit(name, "|", fixed = TRUE)[[1L]]
  result <- hallmark_hypergeometric(
    enrichment_inputs[[name]],
    all_tested_genes,
    term_to_gene
  )
  result$method <- split_name[1L]
  result$gene_set <- split_name[2L]
  result
}))
hallmark_targets <- c(
  "HALLMARK_HYPOXIA",
  "HALLMARK_KRAS_SIGNALING_UP"
)
hallmark_comparison <- hallmark_all[
  hallmark_all$term %in% hallmark_targets,
  ,
  drop = FALSE
]
manuscript_hallmark_reference <- data.frame(
  method = "Current FASH manuscript",
  gene_set = rep(c("All dynamic", "Switch"), each = 2L),
  term = rep(hallmark_targets, 2L),
  overlap_count = c(25L, 11L, 11L, 7L),
  selected_size = c(1177L, 1177L, 250L, 250L),
  p_value = c(0.017, 0.24, 0.00067, 0.0022),
  q_value = c(0.436, 0.840, 0.025, 0.035),
  stringsAsFactors = FALSE
)

message_step(8, 8, "Assembling manuscript-facing cache and exports.")
symbol_lookup <- setNames(gene_map$hgnc_symbol, gene_map$ensembl_gene_id)
add_symbol <- function(table) {
  table$gene_symbol <- unname(symbol_lookup[table$gene_id])
  table$gene_symbol[is.na(table$gene_symbol) | table$gene_symbol == ""] <-
    table$gene_id[is.na(table$gene_symbol) | table$gene_symbol == ""]
  table
}
figure3_proposed <- add_symbol(figure3_proposed)
figure5_proposed <- add_symbol(figure5_proposed)
figure6_proposed <- add_symbol(figure6_proposed)

figure3_status <- figure3_units
figure3_status$current_lfdr <- current_iwp1_lfdr[figure3_status$original_index]
figure3_status$cl_lfdr <- cl_iwp1_adjustment$lfdr[figure3_status$original_index]
figure3_status$current_significant <-
  figure3_status$original_index %in% current_iwp1_discovery$indices
figure3_status$cl_significant <-
  figure3_status$original_index %in% cl_iwp1_discovery$indices

figure5_status <- figure5_units
figure5_status$current_lfdr <- current_iwp2_lfdr[figure5_status$original_index]
figure5_status$cl_lfdr <- cl_iwp2_lfdr[figure5_status$original_index]
figure5_status$current_significant <-
  figure5_status$original_index %in% current_iwp2_discovery$indices
figure5_status$cl_significant <-
  figure5_status$original_index %in% cl_iwp2_discovery$indices

figure6_status <- figure6_units
figure6_status$current_dynamic_significant <-
  figure6_status$original_index %in% current_iwp1_discovery$indices
figure6_status$cl_dynamic_significant <-
  figure6_status$original_index %in% cl_iwp1_discovery$indices
figure6_status$current_lfsr <- vapply(seq_len(nrow(figure6_status)), function(row) {
  value <- current_classification$lfsr[
    current_classification$key == figure6_status$key[row] &
      current_classification$category == figure6_status$role[row]
  ]
  if (length(value) == 1L) value else NA_real_
}, numeric(1))
figure6_status$current_cfsr <- vapply(seq_len(nrow(figure6_status)), function(row) {
  value <- current_classification$cfsr[
    current_classification$key == figure6_status$key[row] &
      current_classification$category == figure6_status$role[row]
  ]
  if (length(value) == 1L) value else NA_real_
}, numeric(1))
figure6_status$cl_lfsr <- vapply(seq_len(nrow(figure6_status)), function(row) {
  value <- cl_classification$lfsr[
    cl_classification$key == figure6_status$key[row] &
      cl_classification$category == figure6_status$role[row]
  ]
  if (length(value) == 1L) value else NA_real_
}, numeric(1))
figure6_status$cl_cfsr <- vapply(seq_len(nrow(figure6_status)), function(row) {
  value <- cl_classification$cfsr[
    cl_classification$key == figure6_status$key[row] &
      cl_classification$category == figure6_status$role[row]
  ]
  if (length(value) == 1L) value else NA_real_
}, numeric(1))
figure6_status$current_category_significant <- figure6_status$current_cfsr <= alpha
figure6_status$cl_category_significant <- figure6_status$cl_cfsr <= alpha

discovery_summary <- data.frame(
  method = rep(c("Current FASH", "FASH-CL"), each = 2L),
  test = rep(c("IWP1 dynamic", "IWP2 nonlinear"), 2L),
  pair_count = c(
    current_iwp1_discovery$pair_count,
    current_iwp2_discovery$pair_count,
    cl_iwp1_discovery$pair_count,
    cl_iwp2_discovery$pair_count
  ),
  gene_count = c(
    current_iwp1_discovery$gene_count,
    current_iwp2_discovery$gene_count,
    cl_iwp1_discovery$gene_count,
    cl_iwp2_discovery$gene_count
  ),
  stringsAsFactors = FALSE
)

current_cl_pair_intersection <- length(intersect(
  current_iwp1_discovery$table$key,
  cl_iwp1_discovery$table$key
))
current_cl_gene_intersection <- length(intersect(
  current_iwp1_discovery$table$gene_id,
  cl_iwp1_discovery$table$gene_id
))
current_cl_variant_intersection <- length(intersect(
  current_iwp1_discovery$table$variant_id,
  cl_iwp1_discovery$table$variant_id
))
method_overlap <- data.frame(
  unit = c("Gene-variant pairs", "Genes", "Variants"),
  current_count = c(
    current_iwp1_discovery$pair_count,
    current_iwp1_discovery$gene_count,
    length(unique(current_iwp1_discovery$table$variant_id))
  ),
  cl_count = c(
    cl_iwp1_discovery$pair_count,
    cl_iwp1_discovery$gene_count,
    length(unique(cl_iwp1_discovery$table$variant_id))
  ),
  intersection = c(
    current_cl_pair_intersection,
    current_cl_gene_intersection,
    current_cl_variant_intersection
  ),
  union = c(
    length(union(
      current_iwp1_discovery$table$key,
      cl_iwp1_discovery$table$key
    )),
    length(union(
      current_iwp1_discovery$table$gene_id,
      cl_iwp1_discovery$table$gene_id
    )),
    length(union(
      current_iwp1_discovery$table$variant_id,
      cl_iwp1_discovery$table$variant_id
    ))
  ),
  stringsAsFactors = FALSE
)
method_overlap$jaccard <- method_overlap$intersection / method_overlap$union

posterior_plot_data <- rbind(
  current_iwp1_plot_bundle$posterior,
  current_iwp2_plot_bundle$posterior,
  cl_iwp1_plot_bundle$posterior,
  cl_iwp2_plot_bundle$posterior
)
observed_plot_data <- rbind(
  current_iwp1_plot_bundle$observed,
  current_iwp2_plot_bundle$observed,
  cl_iwp1_plot_bundle$observed,
  cl_iwp2_plot_bundle$observed
)
parametric_plot_data <- rbind(
  current_iwp1_plot_bundle$parametric,
  current_iwp2_plot_bundle$parametric,
  cl_iwp1_plot_bundle$parametric,
  cl_iwp2_plot_bundle$parametric
)

claim_impact <- data.frame(
  claim = c(
    "FASH identifies many more dynamic eQTL pairs than Strober",
    "FASH identifies many variants missed by parametric models",
    "The flexible EB-LGP model captures diverse nonlinear trajectories",
    "Dynamic eQTLs can be classified by posterior functionals",
    "Switch genes emphasize hypoxia and KRAS-related biology"
  ),
  status = c(
    "Not supported under FASH-CL at pair level",
    "Weakened but still present",
    "Methodologically preserved; examples change",
    "Preserved; counts and examples change",
    "Overall hypoxia persists; switch-specific signal is not preserved"
  ),
  manuscript_action = c(
    "Rewrite the central discovery-count paragraph and Figure 4 interpretation.",
    "Replace absolute superiority language with set-specific novel-discovery counts.",
    "Replace lost Figure 3 and Figure 5 examples and avoid attributing all differences to model flexibility.",
    "Update Table 1, Figure 6 examples, and category-specific text.",
    "Rewrite Table 2: retain all-dynamic hypoxia, but remove the switch-specific hypoxia/KRAS emphasis."
  ),
  stringsAsFactors = FALSE
)

runtime_summary <- data.frame(
  stage = c("Total retained analysis", "FASH-CL IWP1 BF adjustment reused", "FASH-CL IWP2 BF adjustment"),
  elapsed_seconds = c(
    proc.time()[["elapsed"]] - analysis_started,
    cl_iwp1_adjustment$elapsed_seconds,
    cl_iwp2_adjustment$elapsed_seconds
  ),
  stringsAsFactors = FALSE
)
input_provenance <- do.call(rbind, lapply(input_paths, file_provenance))

cache <- list(
  configuration = list(
    analysis_id = "revision_internal_fash_cl_manuscript_impact",
    alpha = alpha,
    smooth_var = smooth_var,
    posterior_sample_size = posterior_sample_size,
    switch_threshold = switch_threshold,
    classification_seed = classification_seed,
    classification_cores = classification_cores,
    classification_batch_size = classification_batch_size,
    fashr_version = as.character(utils::packageVersion("fashr")),
    msigdbr_version = as.character(utils::packageVersion("msigdbr")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  ),
  input_provenance = input_provenance,
  discovery_summary = discovery_summary,
  method_overlap = method_overlap,
  lfdr_scatter_all = lfdr_scatter_all,
  lfdr_scatter_top = lfdr_scatter_top,
  lfdr_comparison_summary = lfdr_comparison_summary,
  prior_weight_comparison = prior_weight_comparison,
  prior_weight_summary = prior_weight_summary,
  venn_sets = venn_sets,
  venn_regions = venn_regions,
  category_counts = category_counts,
  current_classification = current_classification,
  cl_classification = cl_classification,
  figure3_status = figure3_status,
  figure3_proposed = figure3_proposed,
  figure5_status = figure5_status,
  figure5_proposed = figure5_proposed,
  figure6_status = figure6_status,
  figure6_proposed = figure6_proposed,
  hallmark_comparison = hallmark_comparison,
  hallmark_all = hallmark_all,
  manuscript_hallmark_reference = manuscript_hallmark_reference,
  posterior_plot_data = posterior_plot_data,
  observed_plot_data = observed_plot_data,
  parametric_plot_data = parametric_plot_data,
  claim_impact = claim_impact,
  runtime_summary = runtime_summary,
  cl_iwp1_adjustment_summary = cl_iwp1_adjustment[c(
    "null_weight", "alternative_weights", "bf_mixture_status",
    "posterior_mixture_status", "elapsed_seconds"
  )],
  cl_iwp2_adjustment_summary = cl_iwp2_adjustment[c(
    "null_weight", "alternative_weights", "bf_mixture_status",
    "posterior_mixture_status", "elapsed_seconds"
  )],
  session_info = utils::sessionInfo()
)

atomic_save_rds(cache, file.path(output_directory, "analysis_cache.rds"))
atomic_save_rds(
  list(wide = classification_wide, long = cl_classification),
  file.path(output_directory, "classification_results.rds")
)
exports <- list(
  discovery_summary = discovery_summary,
  method_overlap = method_overlap,
  lfdr_comparison_summary = lfdr_comparison_summary,
  lfdr_scatter_top = lfdr_scatter_top,
  iwp1_prior_weights = prior_weight_comparison,
  iwp1_prior_weight_summary = prior_weight_summary,
  venn_regions = venn_regions,
  category_counts = category_counts,
  manuscript_unit_status = rbind(
    transform(figure3_status, unit_type = "discovery"),
    transform(figure5_status, unit_type = "discovery")
  ),
  replacement_examples = rbind(
    data.frame(
      artifact = "Figure 3",
      panel = figure3_proposed$panel,
      key = figure3_proposed$key,
      gene_id = figure3_proposed$gene_id,
      gene_symbol = figure3_proposed$gene_symbol,
      variant_id = figure3_proposed$variant_id,
      status = figure3_proposed$replacement_role,
      statistic = figure3_proposed$lfdr,
      stringsAsFactors = FALSE
    ),
    data.frame(
      artifact = "Figure 5",
      panel = figure5_proposed$panel,
      key = figure5_proposed$key,
      gene_id = figure5_proposed$gene_id,
      gene_symbol = figure5_proposed$gene_symbol,
      variant_id = figure5_proposed$variant_id,
      status = figure5_proposed$replacement_status,
      statistic = figure5_proposed$lfdr,
      stringsAsFactors = FALSE
    ),
    data.frame(
      artifact = "Figure 6",
      panel = figure6_proposed$panel,
      key = figure6_proposed$key,
      gene_id = figure6_proposed$gene_id,
      gene_symbol = figure6_proposed$gene_symbol,
      variant_id = figure6_proposed$variant_id,
      status = figure6_proposed$replacement_status,
      statistic = figure6_proposed$lfsr,
      stringsAsFactors = FALSE
    )
  ),
  hallmark_comparison = hallmark_comparison,
  claim_impact = claim_impact,
  runtime_summary = runtime_summary,
  input_provenance = input_provenance
)
for (name in names(exports)) {
  utils::write.csv(
    exports[[name]],
    file.path(output_directory, paste0(name, ".csv")),
    row.names = FALSE
  )
}

message(
  "FASH-CL manuscript-impact cache completed in ",
  round(proc.time()[["elapsed"]] - analysis_started, 1),
  " seconds."
)
