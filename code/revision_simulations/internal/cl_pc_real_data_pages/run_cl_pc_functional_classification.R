#!/usr/bin/env Rscript

# Recompute all four functional classifications for a corrected CL-PC fit.

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)

find_workflowr_root <- function(start = getwd()) {
  nested_root <- file.path(start, "coderepo-local")
  if (file.exists(file.path(nested_root, "_workflowr.yml"))) {
    return(normalizePath(nested_root, winslash = "/", mustWork = TRUE))
  }
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "_workflowr.yml"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) stop("Could not locate workflowr root.")
    current <- parent
  }
}

parse_arguments <- function(arguments) {
  output <- list(
    order = NA_integer_,
    num_cores = 2L,
    posterior_draws = 3000L,
    seed = 20260820L,
    alpha = 0.05,
    switch_threshold = 0.25,
    batch_size = 64L,
    pilot_pairs = 0L
  )
  names_by_flag <- c(
    "--order" = "order",
    "--num-cores" = "num_cores",
    "--posterior-draws" = "posterior_draws",
    "--seed" = "seed",
    "--alpha" = "alpha",
    "--switch-threshold" = "switch_threshold",
    "--batch-size" = "batch_size",
    "--pilot-pairs" = "pilot_pairs"
  )
  position <- 1L
  while (position <= length(arguments)) {
    flag <- arguments[[position]]
    if (!flag %in% names(names_by_flag) || position == length(arguments)) {
      stop("Unknown or incomplete argument: ", flag)
    }
    output[[names_by_flag[[flag]]]] <- arguments[[position + 1L]]
    position <- position + 2L
  }
  integer_names <- c(
    "order", "num_cores", "posterior_draws", "seed", "batch_size",
    "pilot_pairs"
  )
  output[integer_names] <- lapply(output[integer_names], as.integer)
  output$alpha <- as.numeric(output$alpha)
  output$switch_threshold <- as.numeric(output$switch_threshold)
  output
}

load_exact_object <- function(path, expected_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  if (!identical(loaded, expected_name)) {
    stop("Unexpected object in ", path, ": ", paste(loaded, collapse = ", "))
  }
  environment[[expected_name]]
}

atomic_save_rds <- function(object, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(object, temporary, compress = "gzip")
  if (!file.rename(temporary, path)) stop("Could not atomically save ", path)
  invisible(path)
}

select_cumulative_lfdr <- function(lfdr, alpha) {
  ordering <- order(lfdr, seq_along(lfdr), method = "radix")
  accepted <- which(cumsum(lfdr[ordering]) / seq_along(ordering) <= alpha)
  if (!length(accepted)) integer() else ordering[seq_len(max(accepted))]
}

pair_seed <- function(base_seed, original_index) {
  seed <- (as.double(base_seed) + as.double(original_index)) %%
    as.double(.Machine$integer.max)
  if (seed < 1) seed <- seed + 1
  as.integer(seed)
}

classify_draws <- function(samples, time_grid, switch_threshold) {
  early <- matrixStats::colMaxs(abs(samples[
    time_grid <= 3, , drop = FALSE
  ])) - matrixStats::colMaxs(abs(samples[
    time_grid > 3, , drop = FALSE
  ]))
  middle <- matrixStats::colMaxs(abs(samples[
    time_grid > 3 & time_grid < 12, , drop = FALSE
  ])) - matrixStats::colMaxs(abs(samples[
    time_grid <= 3 | time_grid >= 12, , drop = FALSE
  ]))
  late <- matrixStats::colMaxs(abs(samples[
    time_grid >= 12, , drop = FALSE
  ])) - matrixStats::colMaxs(abs(samples[
    time_grid < 12, , drop = FALSE
  ]))
  positive <- matrixStats::colMaxs(pmax(samples, 0))
  negative <- matrixStats::colMaxs(pmax(-samples, 0))
  switch <- pmin(positive, negative) - switch_threshold
  c(
    early = mean(early <= 0),
    middle = mean(middle <= 0),
    late = mean(late <= 0),
    switch = mean(switch <= 0)
  )
}

make_testing_table <- function(results, category, discovery_order) {
  table <- data.frame(
    indices = as.integer(results$original_index),
    lfsr = as.numeric(results[[category]]),
    stringsAsFactors = FALSE,
    row.names = results$key
  )
  discovery_rank <- match(table$indices, discovery_order)
  if (anyNA(discovery_rank)) stop("A classified pair is outside the discovery set.")
  ranked_order <- order(table$lfsr, discovery_rank, method = "radix")
  table$cfsr <- NA_real_
  table$cfsr[ranked_order] <-
    cumsum(table$lfsr[ranked_order]) / seq_along(ranked_order)
  table[order(table$indices), , drop = FALSE]
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
if (is.na(arguments$order) || !arguments$order %in% c(1L, 2L) ||
    is.na(arguments$num_cores) || arguments$num_cores < 1L ||
    arguments$num_cores > 2L ||
    is.na(arguments$posterior_draws) || arguments$posterior_draws < 100L ||
    is.na(arguments$seed) || arguments$seed < 1L ||
    is.na(arguments$batch_size) || arguments$batch_size < 1L ||
    is.na(arguments$pilot_pairs) || arguments$pilot_pairs < 0L ||
    !is.finite(arguments$alpha) || arguments$alpha <= 0 ||
    arguments$alpha >= 1 ||
    !is.finite(arguments$switch_threshold) ||
    arguments$switch_threshold <= 0) {
  stop("Invalid numerical arguments; --num-cores is capped at two.")
}
if (.Platform$OS.type == "windows" && arguments$num_cores > 1L) {
  stop("Forked classification requires a Unix-like platform.")
}

suppressPackageStartupMessages(library(fashr))
workflowr_root <- find_workflowr_root()
analysis_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "cl_pc_real_data_pages_fashr0143"
)
helper_path <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "fash_cl_manuscript_impact",
  "fash_cl_manuscript_impact_helpers.R"
)
source(helper_path)

if (arguments$order == 1L) {
  fit_path <- file.path(analysis_directory, "fash_fit1_update_CL.RData")
  fit_object <- "fash_fit1_update"
  expected_pairs <- 5395L
  output_directory <- file.path(analysis_directory, "dynamic_classification")
  object_suffix <- "dyn"
  file_stem <- "classify_dyn_eQTLs_"
} else {
  fit_path <- file.path(analysis_directory, "fash_fit2_update_CL.RData")
  fit_object <- "fash_fit2_update"
  expected_pairs <- 60L
  output_directory <- file.path(analysis_directory, "nonlinear_classification")
  object_suffix <- "nonlin_dyn"
  file_stem <- "classify_nonlin_dyn_eQTLs_"
}
if (!file.exists(fit_path)) stop("The corrected CL-PC fit is missing.")
if (file.exists(file.path(output_directory, "classification_summary.csv"))) {
  stop("The finalized classification output already exists.")
}

message("Loading corrected CL-PC order-", arguments$order, " fit.")
fit <- load_exact_object(fit_path, fit_object)
selected_indices <- select_cumulative_lfdr(fit$lfdr, arguments$alpha)
if (length(selected_indices) != expected_pairs || anyDuplicated(selected_indices)) {
  stop("The reconstructed discovery universe does not match its contract.")
}
selected_keys <- names(fit$fash_data$data_list)[selected_indices]
selected_fit <- subset_fash_fit(fit, selected_indices)
rm(fit)
invisible(gc())

time_grid <- seq(0, 15, by = 0.1)
fit_md5 <- unname(tools::md5sum(fit_path))
package_description <- utils::packageDescription("fashr")
configuration <- list(
  order = arguments$order,
  fit_file = normalizePath(fit_path, winslash = "/", mustWork = TRUE),
  fit_md5 = fit_md5,
  package_version = as.character(package_description[["Version"]]),
  package_sha = as.character(package_description[["RemoteSha"]]),
  categories = c("early", "middle", "late", "switch"),
  definitions = c(
    early = "t <= 3 versus t > 3",
    middle = "3 < t < 12 versus t <= 3 or t >= 12",
    late = "t >= 12 versus t < 12",
    switch = "signed excursion on both sides of zero exceeds 0.25"
  ),
  grid = time_grid,
  grid_step = 0.1,
  posterior_draws = arguments$posterior_draws,
  num_cores = arguments$num_cores,
  seed = arguments$seed,
  pair_seed_rule = "base seed plus original fitted pair index",
  alpha = arguments$alpha,
  switch_threshold = arguments$switch_threshold,
  candidate_pairs = expected_pairs,
  original_indices = selected_indices,
  pair_keys = selected_keys
)

pilot_count <- min(arguments$pilot_pairs, expected_pairs)
if (pilot_count > 0L) {
  selected_fit <- subset_fash_fit(selected_fit, seq_len(pilot_count))
  message("Running a non-retained pilot on ", pilot_count, " pairs.")
}

compute_one <- function(local_index) {
  original_index <- selected_fit$original_indices[[local_index]]
  set.seed(pair_seed(arguments$seed, original_index))
  samples <- predict(
    selected_fit,
    index = local_index,
    smooth_var = time_grid,
    only.samples = TRUE,
    M = arguments$posterior_draws
  )
  if (!is.matrix(samples) || nrow(samples) != length(time_grid) ||
      ncol(samples) != arguments$posterior_draws || any(!is.finite(samples))) {
    stop("Unexpected posterior samples for fitted index ", original_index, ".")
  }
  lfsr <- classify_draws(samples, time_grid, arguments$switch_threshold)
  data.frame(
    original_index = original_index,
    key = names(selected_fit$fash_data$data_list)[[local_index]],
    early = unname(lfsr[["early"]]),
    middle = unname(lfsr[["middle"]]),
    late = unname(lfsr[["late"]]),
    switch = unname(lfsr[["switch"]]),
    stringsAsFactors = FALSE
  )
}

if (pilot_count > 0L) {
  started <- proc.time()[["elapsed"]]
  pilot <- parallel::mclapply(
    seq_len(pilot_count),
    compute_one,
    mc.cores = arguments$num_cores,
    mc.preschedule = TRUE,
    mc.set.seed = FALSE
  )
  if (any(vapply(pilot, inherits, logical(1L), "try-error"))) {
    stop("At least one pilot worker failed.")
  }
  elapsed <- proc.time()[["elapsed"]] - started
  message(sprintf(
    "Pilot completed in %.2f seconds (%.3f pair-seconds; projected %.1f minutes).",
    elapsed,
    elapsed / pilot_count,
    elapsed / pilot_count * expected_pairs / 60
  ))
  quit(save = "no", status = 0L)
}

dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
progress_path <- file.path(output_directory, "classification_progress.rds")
progress <- if (file.exists(progress_path)) readRDS(progress_path) else NULL
if (is.null(progress)) {
  progress <- list(
    configuration = configuration,
    results = data.frame(),
    posterior_sampling_seconds = 0
  )
} else if (!identical(progress$configuration, configuration)) {
  stop("The retained checkpoint configuration no longer matches this run.")
}
if (is.null(progress$posterior_sampling_seconds)) {
  stop("The checkpoint is missing cumulative runtime provenance.")
}
completed_indices <- if (nrow(progress$results)) {
  as.integer(progress$results$original_index)
} else {
  integer()
}
remaining <- which(!selected_fit$original_indices %in% completed_indices)
total_batches <- ceiling(length(remaining) / arguments$batch_size)
run_started <- proc.time()[["elapsed"]]

if (length(remaining)) {
  for (batch_number in seq_len(total_batches)) {
    first <- (batch_number - 1L) * arguments$batch_size + 1L
    last <- min(batch_number * arguments$batch_size, length(remaining))
    batch_indices <- remaining[first:last]
    batch_started <- proc.time()[["elapsed"]]
    batch <- parallel::mclapply(
      batch_indices,
      compute_one,
      mc.cores = arguments$num_cores,
      mc.preschedule = TRUE,
      mc.set.seed = FALSE
    )
    if (any(vapply(batch, inherits, logical(1L), "try-error"))) {
      stop("At least one posterior-classification worker failed.")
    }
    progress$results <- rbind(progress$results, do.call(rbind, batch))
    progress$results <- progress$results[
      !duplicated(progress$results$original_index), , drop = FALSE
    ]
    progress$results <- progress$results[
      match(selected_fit$original_indices[
        selected_fit$original_indices %in% progress$results$original_index
      ], progress$results$original_index),
      ,
      drop = FALSE
    ]
    batch_elapsed <- proc.time()[["elapsed"]] - batch_started
    progress$posterior_sampling_seconds <-
      progress$posterior_sampling_seconds + batch_elapsed
    atomic_save_rds(progress, progress_path)
    message(sprintf(
      "Batch %d/%d: %d/%d pairs retained in %.1f seconds.",
      batch_number,
      total_batches,
      nrow(progress$results),
      expected_pairs,
      batch_elapsed
    ))
  }
}

results <- progress$results
expected_keys_in_result_order <- selected_keys[
  match(results$original_index, selected_indices)
]
if (nrow(results) != expected_pairs ||
    !setequal(as.integer(results$original_index), selected_indices) ||
    anyNA(expected_keys_in_result_order) ||
    !identical(as.character(results$key), expected_keys_in_result_order) ||
    any(!is.finite(as.matrix(results[, c(
      "early", "middle", "late", "switch"
    )])))) {
  stop("The completed classification is incomplete or misaligned.")
}

categories <- configuration$categories
tables <- setNames(lapply(categories, function(category) {
  make_testing_table(results, category, selected_indices)
}), categories)
summary <- do.call(rbind, lapply(categories, function(category) {
  table <- tables[[category]]
  selected <- table$cfsr <= arguments$alpha
  selected_keys_category <- rownames(table)[selected]
  data.frame(
    category = category,
    candidate_pairs = nrow(table),
    selected_pairs = sum(selected),
    selected_genes = length(unique(sub("_.*$", "", selected_keys_category))),
    minimum_lfsr = min(table$lfsr),
    maximum_selected_lfsr = if (any(selected)) max(table$lfsr[selected]) else NA,
    maximum_selected_cfsr = if (any(selected)) max(table$cfsr[selected]) else NA,
    stringsAsFactors = FALSE
  )
}))
rownames(summary) <- NULL
configuration$completed_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
configuration$posterior_sampling_seconds <-
  progress$posterior_sampling_seconds

for (category in categories) {
  object_name <- paste0("testing_", category, "_", object_suffix)
  environment <- new.env(parent = emptyenv())
  assign(object_name, tables[[category]], envir = environment)
  save(
    list = object_name,
    file = file.path(output_directory, paste0(file_stem, category, ".RData")),
    envir = environment,
    compress = "gzip"
  )
}
saveRDS(configuration, file.path(output_directory, "run_configuration.rds"))
utils::write.csv(
  summary,
  file.path(output_directory, "classification_summary.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(utils::sessionInfo()),
  file.path(output_directory, "session_info.txt")
)
message("Completed retained CL-PC order-", arguments$order, " classification.")
print(summary, row.names = FALSE)
