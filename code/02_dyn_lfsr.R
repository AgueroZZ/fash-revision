#!/usr/bin/env Rscript

# Canonical functional local-false-sign-rate runner for real-data FASH results.
#
# Direct execution handles order-1 dynamic eQTLs. The order-2 entrypoint
# `03_nonlindyn_lfsr.R` selects the nonlinear specification and reuses this
# implementation so the two production workflows cannot drift apart.

analysis_order <- Sys.getenv("FASH_LFSR_ANALYSIS_ORDER", unset = "1")
if (!analysis_order %in% c("1", "2")) {
  stop("FASH_LFSR_ANALYSIS_ORDER must be 1 or 2.")
}

EXPECTED_FASHR_VERSION <- "0.1.43"
EXPECTED_FASHR_REMOTE_SHA <- "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"

analysis_spec <- if (analysis_order == "1") {
  list(
    order = 1L,
    label = "dynamic eQTL",
    runner = "code/02_dyn_lfsr.R",
    retained_object = "fash_fit1_update",
    raw_object = "fash_fit1",
    expected_pairs = 9214L,
    fit_candidates = c(
      file.path("output", "dynamic_eQTL_real", "fash_fit1_update.RData"),
      file.path("results", "fash_fit1_update.RData"),
      file.path("output", "dynamic_eQTL_real", "fash_fit1_all.RData"),
      file.path("results", "fash_fit1_all.RData")
    ),
    object_names = c(
      early = "testing_early_dyn",
      middle = "testing_middle_dyn",
      late = "testing_late_dyn",
      switch = "testing_switch_dyn"
    ),
    file_names = c(
      early = "classify_dyn_eQTLs_early.RData",
      middle = "classify_dyn_eQTLs_middle.RData",
      late = "classify_dyn_eQTLs_late.RData",
      switch = "classify_dyn_eQTLs_switch.RData"
    )
  )
} else {
  list(
    order = 2L,
    label = "nonlinear dynamic eQTL",
    runner = "code/03_nonlindyn_lfsr.R",
    retained_object = "fash_fit2_update",
    raw_object = "fash_fit2",
    expected_pairs = 44L,
    fit_candidates = c(
      file.path("output", "dynamic_eQTL_real", "fash_fit2_update.RData"),
      file.path("results", "fash_fit2_update.RData"),
      file.path("output", "dynamic_eQTL_real", "fash_fit2_all.RData"),
      file.path("results", "fash_fit2_all.RData")
    ),
    object_names = c(
      early = "testing_early_nonlin_dyn",
      middle = "testing_middle_nonlin_dyn",
      late = "testing_late_nonlin_dyn",
      switch = "testing_switch_nonlin_dyn"
    ),
    file_names = c(
      early = "classify_nonlin_dyn_eQTLs_early.RData",
      middle = "classify_nonlin_dyn_eQTLs_middle.RData",
      late = "classify_nonlin_dyn_eQTLs_late.RData",
      switch = "classify_nonlin_dyn_eQTLs_switch.RData"
    )
  )
}

parse_arguments <- function(args) {
  defaults <- list(
    fit_file = NA_character_,
    output_dir = NA_character_,
    categories = "all",
    grid_step = 0.10,
    posterior_draws = 3000L,
    num_cores = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4")),
    seed = 20260820L,
    alpha = 0.05,
    switch_threshold = 0.25,
    expected_pairs = analysis_spec$expected_pairs,
    allow_bf_update = FALSE,
    dry_run = FALSE,
    help = FALSE
  )

  value_arguments <- c(
    "--fit-file" = "fit_file",
    "--output-dir" = "output_dir",
    "--categories" = "categories",
    "--grid-step" = "grid_step",
    "--posterior-draws" = "posterior_draws",
    "--num-cores" = "num_cores",
    "--seed" = "seed",
    "--alpha" = "alpha",
    "--switch-threshold" = "switch_threshold",
    "--expected-pairs" = "expected_pairs"
  )

  i <- 1L
  while (i <= length(args)) {
    argument <- args[[i]]
    if (argument %in% c("--allow-bf-update", "--dry-run", "--help", "-h")) {
      if (argument == "--allow-bf-update") defaults$allow_bf_update <- TRUE
      if (argument == "--dry-run") defaults$dry_run <- TRUE
      if (argument %in% c("--help", "-h")) defaults$help <- TRUE
      i <- i + 1L
      next
    }
    if (!argument %in% names(value_arguments) || i == length(args)) {
      stop("Unknown or incomplete command-line argument: ", argument)
    }
    defaults[[value_arguments[[argument]]]] <- args[[i + 1L]]
    i <- i + 2L
  }

  defaults$grid_step <- as.numeric(defaults$grid_step)
  defaults$posterior_draws <- as.integer(defaults$posterior_draws)
  defaults$num_cores <- as.integer(defaults$num_cores)
  defaults$seed <- as.integer(defaults$seed)
  defaults$alpha <- as.numeric(defaults$alpha)
  defaults$switch_threshold <- as.numeric(defaults$switch_threshold)
  defaults$expected_pairs <- as.integer(defaults$expected_pairs)
  defaults
}

print_usage <- function() {
  cat(paste(
    paste0("Canonical ", analysis_spec$label, " functional classification"),
    "Middle: 3 < t < 12 versus t <= 3 or t >= 12",
    "",
    "Usage:",
    paste0("  Rscript ", analysis_spec$runner, " [options]"),
    "",
    "Core options:",
    "  --fit-file PATH          BF-adjusted or raw FASH .RData file.",
    "                           If omitted, standard output/ and results/ paths",
    paste0("                           are searched, preferring ",
           analysis_spec$retained_object, ".RData."),
    "  --output-dir PATH        New directory for classification files.",
    "                           Required unless --dry-run is used.",
    "  --categories VALUE       all (default), middle, or a comma-separated",
    "                           subset of early,middle,late,switch.",
    "  --grid-step VALUE        Evaluation-grid step; default 0.10.",
    "  --posterior-draws N      Posterior draws per pair; default 3000.",
    "  --num-cores N            Forked workers; defaults to SLURM_CPUS_PER_TASK",
    "                           or 4 outside Slurm.",
    "  --seed N                 Base seed; default 20260820.",
    "  --alpha VALUE            Discovery and cumulative-FSR cutoff; default 0.05.",
    "  --switch-threshold VALUE Switch amplitude threshold; default 0.25.",
    paste0("  --expected-pairs N       Expected discovery count; default ",
           analysis_spec$expected_pairs, "."),
    "                           Use 0 to disable this provenance check.",
    "  --allow-bf-update        Explicitly permit BF_update() when only raw",
    paste0("                           ", analysis_spec$raw_object,
           " is available. This can change the"),
    "                           retained discovery universe and is off by default.",
    "  --dry-run                Validate inputs and discoveries without sampling.",
    "  --help                    Print this message.",
    sep = "\n"
  ), "\n", sep = "")
}

normalize_existing_file <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) return(NA_character_)
  normalizePath(path, mustWork = TRUE)
}

find_default_fit <- function() {
  candidates <- analysis_spec$fit_candidates
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop(
      "Could not find a standard FASH fit. Supply --fit-file explicitly."
    )
  }
  normalizePath(existing[[1L]], mustWork = TRUE)
}

load_analysis_fit <- function(fit_file, allow_bf_update) {
  fit_environment <- new.env(parent = emptyenv())
  loaded_names <- load(fit_file, envir = fit_environment)

  if (analysis_spec$retained_object %in% loaded_names) {
    fit <- get(analysis_spec$retained_object, envir = fit_environment)
    fit_treatment <- paste(
      "retained BF-adjusted", analysis_spec$retained_object
    )
  } else if (analysis_spec$raw_object %in% loaded_names) {
    if (!allow_bf_update) {
      stop(
        "The input contains raw ", analysis_spec$raw_object,
        " but not the retained ", analysis_spec$retained_object,
        ". An exact classification-only rerun requires the ",
        "retained BF-adjusted fit. Supply that file, or use ",
        "--allow-bf-update only if changing the discovery universe is acceptable."
      )
    }
    message(
      "The input contains raw ", analysis_spec$raw_object, " but not ",
      analysis_spec$retained_object, "; ",
      "applying fashr::BF_update() once."
    )
    fit <- fashr::BF_update(
      get(analysis_spec$raw_object, envir = fit_environment),
      plot = FALSE
    )
    fit_treatment <- paste(analysis_spec$raw_object, "followed by BF_update")
  } else {
    stop(
      "The fit file must contain ", analysis_spec$retained_object, " or ",
      analysis_spec$raw_object, ". Loaded: ",
      paste(loaded_names, collapse = ", ")
    )
  }

  if (!inherits(fit, "fash")) {
    stop("The selected fitted object does not inherit from class 'fash'.")
  }
  list(fit = fit, treatment = fit_treatment)
}

parse_categories <- function(value) {
  category_order <- c("early", "middle", "late", "switch")
  requested <- trimws(tolower(strsplit(value, ",", fixed = TRUE)[[1L]]))
  if (identical(requested, "all")) requested <- category_order
  requested <- unique(requested)
  if (length(requested) == 0L || any(!requested %in% category_order)) {
    stop(
      "--categories must be middle, all, or a comma-separated subset of: ",
      paste(category_order, collapse = ", ")
    )
  }
  category_order[category_order %in% requested]
}

validate_grid <- function(grid_step) {
  if (length(grid_step) != 1L || !is.finite(grid_step) || grid_step <= 0) {
    stop("--grid-step must be one finite positive number.")
  }
  n_intervals <- 15 / grid_step
  if (abs(n_intervals - round(n_intervals)) > 1e-8) {
    stop("--grid-step must divide the 0-to-15 time range exactly.")
  }
  seq(0, 15, by = grid_step)
}

build_functionals <- function(time_grid, switch_threshold) {
  list(
    early = function(x) {
      max(abs(x[time_grid <= 3])) - max(abs(x[time_grid > 3]))
    },
    middle = function(x) {
      max(abs(x[time_grid > 3 & time_grid < 12])) -
        max(abs(x[time_grid <= 3 | time_grid >= 12]))
    },
    late = function(x) {
      max(abs(x[time_grid >= 12])) - max(abs(x[time_grid < 12]))
    },
    switch = function(x) {
      positive <- x[x > 0]
      negative <- x[x < 0]
      if (length(positive) == 0L || length(negative) == 0L) return(0)
      min(max(abs(positive)), max(abs(negative))) - switch_threshold
    }
  )
}

pair_seed <- function(base_seed, pair_index) {
  seed <- (as.double(base_seed) + as.double(pair_index)) %%
    as.double(.Machine$integer.max)
  if (seed < 1) seed <- seed + 1
  as.integer(seed)
}

compute_pair_lfsr <- function(pair_index,
                              fit,
                              time_grid,
                              posterior_draws,
                              category_functionals,
                              base_seed) {
  set.seed(pair_seed(base_seed, pair_index))
  posterior_samples <- predict(
    fit,
    index = pair_index,
    smooth_var = time_grid,
    only.samples = TRUE,
    M = posterior_draws
  )
  if (!is.matrix(posterior_samples) ||
      nrow(posterior_samples) != length(time_grid) ||
      ncol(posterior_samples) != posterior_draws) {
    stop("Unexpected posterior-sample dimensions for pair index ", pair_index)
  }

  vapply(category_functionals, function(functional) {
    statistic_draws <- apply(posterior_samples, 2L, functional)
    mean(statistic_draws <= 0)
  }, numeric(1L))
}

make_testing_table <- function(indices, pair_names, lfsr) {
  result <- data.frame(
    indices = as.integer(indices),
    lfsr = as.numeric(lfsr),
    stringsAsFactors = FALSE
  )
  rownames(result) <- pair_names

  # Match fashr::testing_functional() exactly: rank by LFSR using the original
  # discovery order for ties, compute cumulative FSR, and finally return rows
  # in fitted-index order.
  ranked_order <- order(result$lfsr)
  ranked_cfsr <- cumsum(result$lfsr[ranked_order]) / seq_along(ranked_order)
  result$cfsr <- NA_real_
  result$cfsr[ranked_order] <- ranked_cfsr
  result[order(result$indices), , drop = FALSE]
}

write_atomic_outputs <- function(output_dir,
                                 category_tables,
                                 configuration,
                                 summary_table,
                                 object_names,
                                 file_names) {
  if (file.exists(output_dir) || dir.exists(output_dir)) {
    stop("Refusing to overwrite an existing output path: ", output_dir)
  }

  parent_directory <- dirname(output_dir)
  dir.create(parent_directory, recursive = TRUE, showWarnings = FALSE)
  staging_directory <- paste0(output_dir, ".staging-", Sys.getpid())
  if (file.exists(staging_directory) || dir.exists(staging_directory)) {
    stop("Staging path already exists: ", staging_directory)
  }
  dir.create(staging_directory, recursive = TRUE, showWarnings = FALSE)
  completed <- FALSE
  on.exit({
    if (!completed && dir.exists(staging_directory)) {
      unlink(staging_directory, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  for (category in names(category_tables)) {
    output_environment <- new.env(parent = emptyenv())
    assign(
      object_names[[category]],
      category_tables[[category]],
      envir = output_environment
    )
    save(
      list = object_names[[category]],
      file = file.path(staging_directory, file_names[[category]]),
      envir = output_environment
    )
  }

  saveRDS(configuration, file.path(staging_directory, "run_configuration.rds"))
  write.csv(
    summary_table,
    file.path(staging_directory, "classification_summary.csv"),
    row.names = FALSE
  )
  writeLines(
    capture.output(sessionInfo()),
    file.path(staging_directory, "session_info.txt")
  )

  if (!file.rename(staging_directory, output_dir)) {
    stop("Could not atomically move the staging directory to: ", output_dir)
  }
  completed <- TRUE
  invisible(output_dir)
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
if (arguments$help) {
  print_usage()
  quit(save = "no", status = 0L)
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)
suppressPackageStartupMessages(library(fashr))

installed_fashr_version <- as.character(packageVersion("fashr"))
installed_fashr_sha <- packageDescription("fashr")$RemoteSha
if (!identical(installed_fashr_version, EXPECTED_FASHR_VERSION) ||
    !identical(installed_fashr_sha, EXPECTED_FASHR_REMOTE_SHA)) {
  stop(
    "This classification requires fashr ", EXPECTED_FASHR_VERSION,
    " at ", EXPECTED_FASHR_REMOTE_SHA, ". Installed version/SHA: ",
    installed_fashr_version, "/",
    if (is.null(installed_fashr_sha)) "<missing>" else installed_fashr_sha,
    "."
  )
}

if (is.na(arguments$num_cores) || arguments$num_cores < 1L ||
    is.na(arguments$posterior_draws) || arguments$posterior_draws < 100L ||
    is.na(arguments$seed) || arguments$seed < 1L ||
    !is.finite(arguments$alpha) || arguments$alpha <= 0 ||
    arguments$alpha >= 1 ||
    !is.finite(arguments$switch_threshold) ||
    arguments$switch_threshold <= 0 ||
    is.na(arguments$expected_pairs) || arguments$expected_pairs < 0L) {
  stop("One or more numerical command-line arguments are invalid.")
}

slurm_cores <- suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")))
if (!is.na(slurm_cores) && arguments$num_cores > slurm_cores) {
  stop(
    "--num-cores exceeds SLURM_CPUS_PER_TASK (", slurm_cores, ")."
  )
}
if (!arguments$dry_run &&
    (is.na(arguments$output_dir) || !nzchar(arguments$output_dir))) {
  stop("--output-dir is required unless --dry-run is used.")
}

categories <- parse_categories(arguments$categories)
time_grid <- validate_grid(arguments$grid_step)
fit_file <- if (is.na(arguments$fit_file)) {
  find_default_fit()
} else {
  normalized <- normalize_existing_file(arguments$fit_file)
  if (is.na(normalized)) {
    stop("--fit-file does not exist: ", arguments$fit_file)
  }
  normalized
}

message("Loading fitted model: ", fit_file)
load_start <- proc.time()[["elapsed"]]
fit_record <- load_analysis_fit(
  fit_file = fit_file,
  allow_bf_update = arguments$allow_bf_update
)
fit <- fit_record$fit
load_seconds <- proc.time()[["elapsed"]] - load_start

discovery_result <- fashr::fdr_control(
  fit,
  alpha = arguments$alpha,
  plot = FALSE
)
selected_indices <- as.integer(
  discovery_result$fdr_results$index[
    discovery_result$fdr_results$FDR <= arguments$alpha
  ]
)
if (length(selected_indices) == 0L || anyDuplicated(selected_indices)) {
  stop("The ", analysis_spec$label, " discovery universe is empty or duplicated.")
}
if (arguments$expected_pairs > 0L &&
    length(selected_indices) != arguments$expected_pairs) {
  stop(
    "Expected ", arguments$expected_pairs,
    " ", analysis_spec$label, " pairs but reconstructed ",
    length(selected_indices), "."
  )
}

all_pair_names <- names(fit$fash_data$data_list)
if (is.null(all_pair_names) || anyNA(all_pair_names[selected_indices]) ||
    any(!nzchar(all_pair_names[selected_indices]))) {
  stop("The selected fitted datasets do not have complete pair names.")
}
selected_pair_names <- all_pair_names[selected_indices]

message("Fit treatment: ", fit_record$treatment)
message("Analysis: order ", analysis_spec$order, " ", analysis_spec$label)
message("Discovery pairs: ", length(selected_indices))
message("Requested categories: ", paste(categories, collapse = ", "))
message(
  "Numerical setting: grid step ", arguments$grid_step,
  ", posterior draws ", arguments$posterior_draws,
  ", base seed ", arguments$seed
)
message("Middle definition: 3 < t < 12")
message(sprintf("Fit load and discovery reconstruction: %.1f seconds", load_seconds))

if (arguments$dry_run) {
  message("Dry-run complete; posterior sampling and output writing were skipped.")
  quit(save = "no", status = 0L)
}

functionals <- build_functionals(
  time_grid = time_grid,
  switch_threshold = arguments$switch_threshold
)
requested_functionals <- functionals[categories]

compute_one <- function(pair_index) {
  tryCatch(
    list(
      index = pair_index,
      lfsr = compute_pair_lfsr(
        pair_index = pair_index,
        fit = fit,
        time_grid = time_grid,
        posterior_draws = arguments$posterior_draws,
        category_functionals = requested_functionals,
        base_seed = arguments$seed
      ),
      error = NA_character_
    ),
    error = function(error) {
      list(
        index = pair_index,
        lfsr = rep(NA_real_, length(categories)),
        error = conditionMessage(error)
      )
    }
  )
}

sampling_start <- proc.time()[["elapsed"]]
if (arguments$num_cores > 1L && .Platform$OS.type != "windows") {
  pair_results <- parallel::mclapply(
    selected_indices,
    compute_one,
    mc.cores = arguments$num_cores,
    mc.set.seed = FALSE,
    mc.preschedule = TRUE
  )
} else {
  pair_results <- lapply(selected_indices, compute_one)
}
sampling_seconds <- proc.time()[["elapsed"]] - sampling_start

failed <- vapply(pair_results, function(result) {
  !is.na(result$error) || any(!is.finite(result$lfsr))
}, logical(1L))
if (any(failed)) {
  first_failure <- pair_results[[which(failed)[1L]]]
  stop(
    "Posterior sampling failed for ", sum(failed), " pair(s). First failure: ",
    first_failure$index, " — ", first_failure$error
  )
}

observed_indices <- vapply(pair_results, `[[`, integer(1L), "index")
if (!identical(observed_indices, selected_indices)) {
  stop("Parallel output order does not match the selected pair order.")
}
lfsr_matrix <- do.call(rbind, lapply(pair_results, `[[`, "lfsr"))
colnames(lfsr_matrix) <- categories

category_tables <- setNames(lapply(categories, function(category) {
  make_testing_table(
    indices = selected_indices,
    pair_names = selected_pair_names,
    lfsr = lfsr_matrix[, category]
  )
}), categories)

summary_table <- do.call(rbind, lapply(categories, function(category) {
  table <- category_tables[[category]]
  selected <- table$cfsr <= arguments$alpha
  data.frame(
    category = category,
    candidate_pairs = nrow(table),
    selected_pairs = sum(selected),
    minimum_lfsr = min(table$lfsr),
    maximum_selected_lfsr = if (any(selected)) {
      max(table$lfsr[selected])
    } else {
      NA_real_
    },
    maximum_selected_cfsr = if (any(selected)) {
      max(table$cfsr[selected])
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}))
rownames(summary_table) <- NULL

output_dir <- normalizePath(
  arguments$output_dir,
  mustWork = FALSE
)
configuration <- list(
  generated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
  analysis_order = analysis_spec$order,
  analysis_label = analysis_spec$label,
  runner = analysis_spec$runner,
  fashr_version = installed_fashr_version,
  fashr_remote_sha = installed_fashr_sha,
  fit_file = fit_file,
  fit_treatment = fit_record$treatment,
  categories = categories,
  middle_definition = "3 < t < 12",
  middle_complement = "t <= 3 or t >= 12",
  grid = time_grid,
  grid_step = arguments$grid_step,
  posterior_draws = arguments$posterior_draws,
  num_cores = arguments$num_cores,
  seed = arguments$seed,
  pair_seed_rule = "base seed plus fitted pair index",
  alpha = arguments$alpha,
  switch_threshold = arguments$switch_threshold,
  discovery_pair_count = length(selected_indices),
  fit_load_and_discovery_seconds = load_seconds,
  posterior_sampling_seconds = sampling_seconds
)

write_atomic_outputs(
  output_dir = output_dir,
  category_tables = category_tables,
  configuration = configuration,
  summary_table = summary_table,
  object_names = analysis_spec$object_names,
  file_names = analysis_spec$file_names
)

message(sprintf("Posterior sampling: %.1f seconds", sampling_seconds))
message("Completed output: ", output_dir)
print(summary_table, row.names = FALSE)
