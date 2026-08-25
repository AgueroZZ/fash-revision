#!/usr/bin/env Rscript

# Focused retained-artifact tests for the two CL-PC real-data pages.

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

load_exact_object <- function(path, expected_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  if (!identical(loaded, expected_name)) {
    stop("Unexpected object in ", path, ".")
  }
  environment[[expected_name]]
}

select_cumulative_lfdr <- function(lfdr, alpha = 0.05) {
  ordering <- order(lfdr, seq_along(lfdr), method = "radix")
  accepted <- which(cumsum(lfdr[ordering]) / seq_along(ordering) <= alpha)
  if (!length(accepted)) integer() else ordering[seq_len(max(accepted))]
}

workflowr_root <- find_workflowr_root()
raw_directory <- file.path(workflowr_root, "output", "dynamic_eQTL_real")
output_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "cl_pc_real_data_pages_fashr0143"
)
matched_directory <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "internal",
  "matched_fash_linear_real_data_cl_fashr_0_1_43"
)
raw_paths <- c(
  iwp1 = file.path(raw_directory, "fash_fit1_all_CL.RData"),
  iwp2 = file.path(raw_directory, "fash_fit2_all_CL.RData")
)
compact_paths <- c(
  iwp1 = file.path(matched_directory, "cl_iwp1_bf_adjustment.rds"),
  iwp2 = file.path(matched_directory, "cl_iwp2_bf_adjustment.rds")
)
expected_raw_md5 <- c(
  iwp1 = "57383f6068304370930f2410021808cb",
  iwp2 = "550fd3d6bb5cf7a5a6d80288da7974a2"
)
stopifnot(
  all(file.exists(c(raw_paths, compact_paths))),
  identical(unname(tools::md5sum(raw_paths)), unname(expected_raw_md5))
)

package_description <- utils::packageDescription("fashr")
stopifnot(
  identical(as.character(package_description[["Version"]]), "0.1.43"),
  identical(
    as.character(package_description[["RemoteSha"]]),
    "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
  )
)

if (dir.exists(output_directory)) {
  cache_path <- file.path(output_directory, "analysis_cache.rds")
  stopifnot(file.exists(cache_path))
  cache <- readRDS(cache_path)
  stopifnot(
    identical(cache$configuration$analysis_id, "cl_pc_real_data_pages_fashr0143"),
    identical(cache$configuration$package_version, "0.1.43"),
    identical(
      cache$configuration$package_sha,
      "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
    ),
    all(cache$validation$passed),
    identical(
      as.integer(cache$discovery_counts$pair_count),
      c(5395L, 60L)
    ),
    identical(
      as.integer(cache$discovery_counts$gene_count),
      c(686L, 6L)
    )
  )

  fit_paths <- c(
    iwp1 = file.path(output_directory, "fash_fit1_update_CL.RData"),
    iwp2 = file.path(output_directory, "fash_fit2_update_CL.RData")
  )
  expected_objects <- c(iwp1 = "fash_fit1_update", iwp2 = "fash_fit2_update")
  expected_pairs <- c(iwp1 = 5395L, iwp2 = 60L)
  expected_genes <- c(iwp1 = 686L, iwp2 = 6L)
  expected_orders <- c(iwp1 = 1L, iwp2 = 2L)
  for (label in names(fit_paths)) {
    stopifnot(file.exists(fit_paths[[label]]))
    fit <- load_exact_object(fit_paths[[label]], expected_objects[[label]])
    compact <- readRDS(compact_paths[[label]])
    selected <- select_cumulative_lfdr(fit$lfdr)
    pair_keys <- names(fit$fash_data$data_list)
    stopifnot(
      inherits(fit, "fash"),
      as.integer(fit$settings$order) == expected_orders[[label]],
      length(pair_keys) == 1009173L,
      !anyDuplicated(pair_keys),
      identical(as.numeric(fit$lfdr), as.numeric(compact$lfdr)),
      length(selected) == expected_pairs[[label]],
      length(unique(sub("_.*$", "", pair_keys[selected]))) ==
        expected_genes[[label]]
    )
    rm(fit, compact)
    invisible(gc())
  }

  classification_specs <- list(
    dynamic = list(
      directory = "dynamic_classification",
      candidate_pairs = 5395L,
      object_suffix = "dyn",
      selected_pairs = c(1L, 98L, 17L, 593L),
      selected_genes = c(1L, 22L, 7L, 150L)
    ),
    nonlinear = list(
      directory = "nonlinear_classification",
      candidate_pairs = 60L,
      object_suffix = "nonlin_dyn",
      selected_pairs = c(0L, 42L, 8L, 5L),
      selected_genes = c(0L, 3L, 3L, 2L)
    )
  )
  for (spec in classification_specs) {
    directory <- file.path(output_directory, spec$directory)
    if (!dir.exists(directory)) next
    summary_path <- file.path(directory, "classification_summary.csv")
    configuration_path <- file.path(directory, "run_configuration.rds")
    progress_path <- file.path(directory, "classification_progress.rds")
    stopifnot(
      file.exists(summary_path),
      file.exists(configuration_path),
      file.exists(progress_path)
    )
    summary <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
    configuration <- readRDS(configuration_path)
    progress <- readRDS(progress_path)
    stopifnot(
      identical(summary$category, c("early", "middle", "late", "switch")),
      all(summary$candidate_pairs == spec$candidate_pairs),
      identical(as.integer(summary$selected_pairs), spec$selected_pairs),
      identical(as.integer(summary$selected_genes), spec$selected_genes),
      configuration$posterior_draws == 3000L,
      configuration$num_cores == 2L,
      configuration$seed == 20260820L,
      configuration$alpha == 0.05,
      configuration$switch_threshold == 0.25,
      identical(configuration$categories, c("early", "middle", "late", "switch")),
      is.finite(configuration$posterior_sampling_seconds),
      configuration$posterior_sampling_seconds > 0,
      identical(
        progress$posterior_sampling_seconds,
        configuration$posterior_sampling_seconds
      )
    )
    for (category in summary$category) {
      path <- file.path(directory, paste0(
        "classify_",
        if (identical(spec$object_suffix, "dyn")) "dyn" else "nonlin_dyn",
        "_eQTLs_", category, ".RData"
      ))
      expected_object <- paste0("testing_", category, "_", spec$object_suffix)
      table <- load_exact_object(path, expected_object)
      stopifnot(
        nrow(table) == spec$candidate_pairs,
        all(c("indices", "lfsr", "cfsr") %in% names(table)),
        all(is.finite(table$lfsr)),
        all(is.finite(table$cfsr))
      )
    }
  }
}

message("All CL-PC real-data page artifact tests passed.")
