#!/usr/bin/env Rscript

# Full-data FASH target-decoy (fash-knockoff) analysis.
#
# One global donor permutation P is applied to the entire genotype matrix. For
# every one of the tested gene-variant pairs a permuted decoy unit is built,
# and target plus decoy units are fitted jointly so a single empirical-Bayes
# prior scores both sides of every pair.
#
# The likelihood rows for the target units are NOT recomputed: they are taken
# verbatim from the immutable full-data FASH object. Only the decoy units
# require new likelihood evaluations.
#
# The run is split into three stages that must execute in separate R
# processes. The full FASH object is several gigabytes, and forking likelihood
# workers from a session that has loaded it is what produced the anomalous
# runtimes recorded in earlier logs.
#
#   --stage prepare     load the source fit, permute the genotype matrix,
#                       refit the time-specific PC model, cache decoy beta/SE
#   --stage likelihood  chunked, resumable decoy likelihood evaluation
#   --stage merge       merged empirical-Bayes refit, BF update, competition
#                       statistics, knockoff+ selection
#
# Example:
#   Rscript code/revision_simulations/internal/fash_knockoff/run_fash_knockoff.R \
#     --stage prepare --seed 20260823
#   Rscript ... --stage likelihood --seed 20260823 --num-cores 8
#   Rscript ... --stage merge --seed 20260823

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
  utils::write.csv(x, file = path, row.names = FALSE)
}

file_metadata <- function(paths, roles) {
  paths <- normalizePath(paths, mustWork = TRUE)
  information <- file.info(paths)
  data.frame(
    role = roles,
    path = paths,
    size_bytes = unname(information$size),
    modification_time = format(information$mtime, tz = "UTC", usetz = TRUE),
    md5 = unname(tools::md5sum(paths)),
    stringsAsFactors = FALSE
  )
}

capture_warnings <- function(expression) {
  warning_messages <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(warning) {
      warning_messages <<- c(warning_messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warning_messages))
}

workflowr_root <- find_workflowr_root()
project_root <- normalizePath(file.path(workflowr_root, ".."), mustWork = FALSE)
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

# ---------------------------------------------------------------------------
# Arguments and fixed design constants
# ---------------------------------------------------------------------------

stage <- get_arg("--stage", "")
expected_fashr_version <- get_arg("--expected-fashr-version", "0.1.43")
expected_fashr_remote_sha <- get_arg(
  "--expected-fashr-remote-sha",
  "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
)
preflight_only <- identical(get_arg("--preflight-only", "false"), "true")
seed <- as.integer(get_arg("--seed", "20260823"))
alpha <- as.numeric(get_arg("--alpha", "0.05"))
num_cores <- as.integer(get_arg("--num-cores", "8"))
chunk_units <- as.integer(get_arg("--chunk-units", "25000"))
validation_sample <- as.integer(get_arg("--validation-sample", "2000"))
# Reduced-scale smoke runs keep whole genes so the gene-level maximum on both
# sides still runs over a gene's complete tested-variant set.
gene_limit <- as.integer(get_arg("--gene-limit", "0"))
default_output_id <- if (gene_limit > 0L) {
  paste0("fash_knockoff_genes", gene_limit, "_seed", seed)
} else {
  paste0("fash_knockoff_full_seed", seed)
}
output_id <- get_arg("--output-id", default_output_id)
allowed_stages <- c("prepare", "likelihood", "merge")
if (!stage %in% allowed_stages ||
    !nzchar(expected_fashr_version) || !nzchar(expected_fashr_remote_sha) || length(seed) != 1L || is.na(seed) ||
    !is.finite(alpha) || alpha <= 0 || alpha >= 1 ||
    is.na(num_cores) || num_cores < 1L ||
    is.na(chunk_units) || chunk_units < 100L ||
    is.na(validation_sample) || validation_sample < 0L ||
    is.na(gene_limit) || gene_limit < 0L ||
    !nzchar(output_id) || grepl("/", output_id, fixed = TRUE)) {
  stop("Invalid stage, seed, alpha, core count, chunk size, or output ID.")
}

# The primary analysis these units come from.
source_pair_count <- 1009173L
expected_grid_size <- 52L
time_grid <- 0:15
expected_sample_counts <- c(
  19L, 19L, 16L, 19L, 16L, 19L, 19L, 19L,
  19L, 19L, 19L, 19L, 19L, 18L, 19L, 19L
)
q_grid <- c(0.01, 0.05, 0.10, 0.20, 0.30, 0.50)

# Input resolution follows the pattern in code/02_dyn_lfsr.R: an explicit
# environment override wins, otherwise a documented candidate list is searched
# relative to the project root, and the resolved path is always reported. This
# lets the same script run from the local workflowr tree, from the Midway3
# workspace, and from a server-side checkout, without any absolute path being
# baked in.
path_from_environment <- function(variable, default) {
  value <- Sys.getenv(variable, unset = "")
  if (nzchar(value)) value else default
}

# Collects every unresolved input before failing. Dying on the first one costs
# a full round trip to the cluster to learn about the second.
unresolved_inputs <- list()

resolve_input <- function(role, override_variable, candidates) {
  override <- Sys.getenv(override_variable, unset = "")
  if (nzchar(override)) {
    if (!file.exists(override)) {
      stop(
        role, " was set through ", override_variable,
        " but does not exist: ", override
      )
    }
    return(normalizePath(override, mustWork = TRUE))
  }
  candidates <- unique(candidates)
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  unresolved_inputs[[length(unresolved_inputs) + 1L]] <<- list(
    role = role, variable = override_variable, candidates = candidates
  )
  NA_character_
}

report_unresolved_inputs <- function() {
  if (length(unresolved_inputs) == 0L) {
    return(invisible(NULL))
  }
  lines <- unlist(lapply(unresolved_inputs, function(entry) {
    c(
      paste0("Could not resolve ", entry$role, ". Set ", entry$variable,
             ", or place the file at one of:"),
      paste0("  ", entry$candidates)
    )
  }))
  stop(
    length(unresolved_inputs), " input(s) could not be resolved:\n",
    paste(lines, collapse = "\n"),
    call. = FALSE
  )
}

# Search roots, most specific first.
search_roots <- unique(c(
  Sys.getenv("FASH_KNOCKOFF_PROJECT_ROOT", unset = ""),
  workflowr_root,
  project_root,
  file.path(project_root, "coderepo-local"),
  "/project/mstephens/ziangzhang/fash",
  "/project/mstephens/ziangzhang/fash/workspace"
))
search_roots <- search_roots[nzchar(search_roots) & dir.exists(search_roots)]
candidates_for <- function(relative_paths) {
  as.vector(t(outer(search_roots, relative_paths, file.path)))
}

# Either the raw fit (`fash_fit1`) or the BF-updated fit (`fash_fit1_update`)
# is accepted. `BF_update()` does not touch `L_matrix`, so the target
# likelihood rows are bit-identical either way; only the reported source lfdr
# differs, and it is derived here when absent.
source_fit_path <- resolve_input(
  "the full FASH fit", "FASH_KNOCKOFF_SOURCE_FIT",
  candidates_for(c(
    file.path("output", "dynamic_eQTL_real", "fash_fit1_update.RData"),
    file.path("output", "dynamic_eQTL_real", "fash_fit1_all.RData"),
    file.path("full_results", "application", "result_fash",
              "fash_fit1_update.RData"),
    file.path("full_results", "application", "result_fash",
              "fash_fit1_all.RData"),
    file.path("results", "fash_fit1_update.RData"),
    file.path("results", "fash_fit1_all.RData")
  ))
)
vcf_path <- resolve_input(
  "the genotype VCF", "FASH_KNOCKOFF_VCF",
  candidates_for(c(
    file.path("iPSC-data", "genotype-data", "YRI_genotype.vcf.gz"),
    file.path("data", "raw", "YRI_genotype.vcf.gz"),
    file.path("data", "genotype-data", "YRI_genotype.vcf.gz"),
    file.path("data", "YRI_genotype.vcf.gz")
  ))
)
expression_path <- resolve_input(
  "the expression matrix", "FASH_KNOCKOFF_EXPRESSION",
  candidates_for(c(
    file.path("iPSC-data", "expression-data",
              "quantile_normalized_no_projection.txt"),
    file.path("data", "raw", "quantile_normalized_no_projection.txt"),
    file.path("data", "expression-data",
              "quantile_normalized_no_projection.txt"),
    file.path("data", "quantile_normalized_no_projection.txt")
  ))
)
pc_path <- resolve_input(
  "the time-specific PCs", "FASH_KNOCKOFF_PC",
  candidates_for(c(
    file.path("data", "dynamic_eQTL_real", "principal_components_10.txt"),
    file.path("data", "raw", "principal_components_10.txt"),
    file.path("data", "principal_components_10.txt")
  ))
)
output_parent <- path_from_environment(
  "FASH_KNOCKOFF_OUTPUT_PARENT",
  file.path(workflowr_root, "output", "revision_simulations", "internal")
)

# Small diagnostic tables are copied here as well, so they stay inspectable
# when the analysis writes its bulk output to a large-results directory.
summary_directory <- path_from_environment("FASH_KNOCKOFF_SUMMARY_DIR", "")

report_unresolved_inputs()

required_paths <- c(source_fit_path, vcf_path, expression_path, pc_path)
required_roles <- c(
  "full_fash_fit", "genotype_vcf", "expression_matrix",
  "time_specific_pc_data"
)
if (any(!file.exists(required_paths))) {
  stop(
    "At least one required input is missing: ",
    paste(required_paths[!file.exists(required_paths)], collapse = ", ")
  )
}

fashr_description <- utils::packageDescription("fashr")
observed_fashr_version <- as.character(fashr_description$Version)
observed_fashr_remote_sha <- if (is.null(fashr_description$RemoteSha)) {
  NA_character_
} else {
  as.character(fashr_description$RemoteSha)
}
message(
  "R ", getRversion(), " | fashr ", observed_fashr_version,
  " | RemoteSha ", observed_fashr_remote_sha
)
if (!identical(observed_fashr_version, expected_fashr_version)) {
  stop(
    "Expected fashr ", expected_fashr_version, "; found ",
    observed_fashr_version, "."
  )
}
if (!identical(observed_fashr_remote_sha, expected_fashr_remote_sha)) {
  stop(
    "Expected fashr RemoteSha ", expected_fashr_remote_sha, "; found ",
    observed_fashr_remote_sha, "."
  )
}

if (preflight_only) {
  message("Resolved input paths:")
  for (position in seq_along(required_paths)) {
    message(
      "  ", required_roles[position], ": ", required_paths[position],
      " (", format(file.info(required_paths[position])$size, big.mark = ","),
      " bytes)"
    )
  }
  message("Output parent: ", output_parent)
  message("Preflight passed; no analysis was run.")
  quit(save = "no", status = 0L)
}

work_directory <- file.path(output_parent, output_id)
chunk_directory <- file.path(work_directory, "decoy_likelihood_chunks")
prepared_path <- file.path(work_directory, "prepared_decoy_input.rds")
target_likelihood_path <- file.path(work_directory, "target_likelihood.rds")
dir.create(work_directory, recursive = TRUE, showWarnings = FALSE)

message("Stage: ", stage, " | seed: ", seed, " | output: ", output_id)

# ===========================================================================
# Stage 1: prepare the decoy beta/SE inputs
# ===========================================================================

if (identical(stage, "prepare")) {
  stage_start <- proc.time()[["elapsed"]]
  source_information_before <- file_metadata(required_paths, required_roles)

  message("Loading the full FASH fit: ", source_fit_path)
  bf_environment <- new.env(parent = emptyenv())
  loaded_names <- load(source_fit_path, envir = bf_environment)
  if (length(loaded_names) != 1L ||
      !loaded_names %in% c("fash_fit1", "fash_fit1_update")) {
    stop(
      "The source fit file must contain exactly one object named fash_fit1 ",
      "or fash_fit1_update; found: ", paste(loaded_names, collapse = ", ")
    )
  }
  full_bf <- bf_environment[[loaded_names]]
  source_fit_is_bf_updated <- identical(loaded_names, "fash_fit1_update")
  message(
    "  object: ", loaded_names,
    " (", if (source_fit_is_bf_updated) "BF-updated" else "raw", ")"
  )
  pair_keys <- names(full_bf$fash_data$data_list)
  if (length(pair_keys) != source_pair_count ||
      !identical(names(full_bf$lfdr), pair_keys) ||
      length(full_bf$psd_grid) != expected_grid_size ||
      nrow(full_bf$L_matrix) != source_pair_count ||
      ncol(full_bf$L_matrix) != expected_grid_size ||
      !identical(full_bf$settings$order, 1) ||
      !identical(full_bf$settings$pred_step, 1) ||
      !identical(full_bf$settings$penalty, 10)) {
    stop("The full BF-adjusted source fit does not match the primary analysis.")
  }
  if (!is.null(rownames(full_bf$L_matrix)) &&
      !identical(rownames(full_bf$L_matrix), pair_keys)) {
    stop("The full likelihood matrix is not aligned by pair key.")
  }

  original_psd_grid <- full_bf$psd_grid
  original_settings <- full_bf$settings

  # The published discovery rule is applied to the BF-updated lfdr, so derive
  # it when the raw fit was supplied. This is only a reporting column: the
  # competition statistic never uses it.
  if (!source_fit_is_bf_updated) {
    message("Deriving the BF-updated lfdr from the raw fit for reporting.")
    bf_lfdr_capture <- capture_warnings(
      fashr::BF_update(full_bf, plot = FALSE)$lfdr
    )
    source_bf_lfdr <- as.numeric(bf_lfdr_capture$value)
    if (length(source_bf_lfdr) != length(pair_keys)) {
      stop("The derived BF-updated lfdr is not aligned with the pair keys.")
    }
  } else {
    source_bf_lfdr <- as.numeric(full_bf$lfdr)
  }

  # Restrict to the analysis scope. A reduced run keeps every tested variant
  # of the retained genes, never a variant-level subsample.
  source_rows <- seq_len(source_pair_count)
  if (gene_limit > 0L) {
    source_gene_id <- sub("_.*$", "", pair_keys)
    retained_genes <- unique(source_gene_id)[seq_len(gene_limit)]
    if (anyNA(retained_genes)) {
      stop("The requested gene limit exceeds the number of tested genes.")
    }
    source_rows <- which(source_gene_id %in% retained_genes)
    message(
      "REDUCED SCOPE: ", gene_limit, " gene(s), ",
      format(length(source_rows), big.mark = ","),
      " pair(s). This is not the full-data analysis."
    )
  }
  pair_keys <- pair_keys[source_rows]
  n_pair <- length(pair_keys)
  target_lfdr <- source_bf_lfdr[source_rows]

  # Reuse, never recompute, the target likelihood rows. These are cached to
  # disk only after the informativeness screen below, so that the cache and the
  # decoy inputs describe exactly the same pair set.
  message("Extracting the immutable target likelihood rows.")
  target_likelihood <- full_bf$L_matrix[source_rows, , drop = FALSE]
  rownames(target_likelihood) <- pair_keys

  # Keep a small slice of the source beta/SE so the regression pipeline can be
  # validated against the published analysis before the decoy is trusted.
  set.seed(seed)
  validation_sample <- min(validation_sample, n_pair)
  validation_rows <- if (validation_sample > 0L) {
    sort(sample.int(n_pair, validation_sample))
  } else {
    integer(0)
  }
  # Both validation slices are keyed by pair key, not by row position, because
  # the informativeness screen below reindexes the pair set. A candidate pool
  # larger than the 50 rows actually needed leaves room for screened-out pairs.
  validation_keys <- pair_keys[validation_rows]
  check_pool_rows <- if (validation_sample > 0L) {
    seq_len(min(200L, n_pair))
  } else {
    integer(0)
  }
  check_pool_keys <- pair_keys[check_pool_rows]
  check_pool_likelihood <- target_likelihood[check_pool_rows, , drop = FALSE]
  check_pool_beta <- NULL
  check_pool_adjusted_se <- NULL
  if (length(check_pool_rows) > 0L) {
    check_pool_beta <- t(vapply(
      full_bf$fash_data$data_list[source_rows[check_pool_rows]],
      function(unit) as.numeric(unit$y), numeric(length(time_grid))
    ))
    check_pool_adjusted_se <- t(vapply(
      full_bf$fash_data$S[source_rows[check_pool_rows]], as.numeric,
      numeric(length(time_grid))
    ))
  }
  source_beta <- NULL
  source_adjusted_se <- NULL
  if (length(validation_rows) > 0L) {
    source_beta <- t(vapply(
      full_bf$fash_data$data_list[source_rows[validation_rows]],
      function(unit) {
        if (!identical(as.numeric(unit$x), as.numeric(time_grid))) {
          stop("A source FASH unit does not use the expected time grid.")
        }
        as.numeric(unit$y)
      },
      numeric(length(time_grid))
    ))
    source_adjusted_se <- t(vapply(
      full_bf$fash_data$S[source_rows[validation_rows]], as.numeric,
      numeric(length(time_grid))
    ))
  }
  rm(full_bf, bf_environment)
  invisible(gc(verbose = FALSE))

  pair_index <- build_pair_index(pair_keys)
  message(
    "Pairs: ", format(n_pair, big.mark = ","),
    " | genes: ", format(length(pair_index$gene_levels), big.mark = ","),
    " | unique variants: ",
    format(length(pair_index$variant_levels), big.mark = ",")
  )

  message("Reading expression and time-specific PCs.")
  expression_data <- utils::read.csv(
    expression_path, sep = "", check.names = FALSE, stringsAsFactors = FALSE
  )
  pc_data <- utils::read.delim(
    pc_path, check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!identical(names(expression_data)[1L], "Gene_id") ||
      !all(c("Sample_id", paste0("PC", 1:5)) %in% names(pc_data)) ||
      anyDuplicated(expression_data$Gene_id) ||
      anyDuplicated(pc_data$Sample_id)) {
    stop("The expression or PC input has invalid identifiers.")
  }
  expression_sample_ids <- names(expression_data)[-1L]
  if (length(expression_sample_ids) != 297L || nrow(pc_data) != 297L ||
      !setequal(expression_sample_ids, pc_data$Sample_id)) {
    stop("Expression and PC donor-time IDs do not match exactly.")
  }
  gene_rows <- match(pair_index$gene_levels, expression_data$Gene_id)
  if (anyNA(gene_rows)) {
    stop("At least one tested gene is missing from the expression matrix.")
  }

  # Streaming 745,867 variants out of the VCF dominates this stage, so the
  # dosage matrix is cached. A retry after a later failure must not pay for it
  # again. The cache is keyed on the VCF checksum and the exact variant list.
  dosage_cache_path <- file.path(work_directory, "tested_variant_dosage.rds")
  vcf_md5 <- source_information_before$md5[
    source_information_before$role == "genotype_vcf"
  ]
  dosage <- NULL
  if (file.exists(dosage_cache_path)) {
    cached_dosage <- readRDS(dosage_cache_path)
    if (identical(cached_dosage$cache_version, "fash_knockoff_dosage_v1") &&
        identical(cached_dosage$vcf_md5, vcf_md5) &&
        identical(cached_dosage$variant_levels, pair_index$variant_levels)) {
      message(
        "Reusing the cached tested-variant dosage matrix (",
        ncol(cached_dosage$dosage), " variants)."
      )
      dosage <- cached_dosage$dosage
    } else {
      message("The dosage cache does not match this run; rebuilding it.")
    }
    rm(cached_dosage)
  }
  if (is.null(dosage)) {
    message(
      "Streaming ", format(length(pair_index$variant_levels), big.mark = ","),
      " tested variants from the VCF."
    )
    vcf_start <- proc.time()[["elapsed"]]
    dosage <- read_selected_vcf_dosages(
      vcf_path, pair_index$variant_levels, chunk_size = 100000L
    )
    dosage <- dosage[, match(pair_index$variant_levels, colnames(dosage)),
                     drop = FALSE]
    message(
      "  VCF streaming took ",
      format((proc.time()[["elapsed"]] - vcf_start) / 60, digits = 4),
      " minutes; caching the dosage matrix."
    )
    saveRDS(
      list(
        cache_version = "fash_knockoff_dosage_v1",
        vcf_md5 = vcf_md5,
        variant_levels = pair_index$variant_levels,
        dosage = dosage
      ),
      dosage_cache_path,
      compress = FALSE
    )
  }
  vcf_donors <- rownames(dosage)
  if (nrow(dosage) != 19L || anyDuplicated(vcf_donors) ||
      !identical(colnames(dosage), pair_index$variant_levels)) {
    stop("The tested-variant dosage matrix is not aligned or complete.")
  }

  # One permutation for the entire genotype matrix: a single donor row map
  # shared by every variant and every time point. Row permutation leaves each
  # variant's dosage multiset and the full LD structure exactly unchanged.
  message("Building the single global genotype permutation.")
  permutation <- make_global_derangement(dosage, seed)
  permuted_dosage <- permutation$genotype
  donor_map <- permutation$donor_map
  # Invariants are asserted inside permute_genotype_donor_rows() in O(19n);
  # a full variant-by-variant crossprod would need 4 TB at this scale.
  if (any(donor_map$fixed_point) ||
      !identical(dim(dosage), dim(permuted_dosage)) ||
      !identical(colnames(dosage), colnames(permuted_dosage)) ||
      !identical(rownames(dosage), rownames(permuted_dosage))) {
    stop("The global donor permutation failed its invariants.")
  }
  message(
    "  derangement found at seed ", permutation$permutation_seed,
    " after ", permutation$attempts, " attempt(s); LD structure preserved."
  )
  genotype_sumsq <- colSums(dosage^2)

  # Assemble the per-time-point donor sets and covariates once, then screen for
  # variants that carry no residual genotype information on either arm.
  donor_sets <- vector("list", length(time_grid))
  covariate_list <- vector("list", length(time_grid))
  for (time_position in seq_along(time_grid)) {
    time_value <- time_grid[time_position]
    sample_ids <- grep(
      paste0("_", time_value, "$"), expression_sample_ids, value = TRUE
    )
    donors <- sub(paste0("_", time_value, "$"), "", sample_ids)
    if (length(sample_ids) != expected_sample_counts[time_position] ||
        anyDuplicated(donors) || any(!donors %in% vcf_donors)) {
      stop("Unexpected donor coverage at time ", time_value, ".")
    }
    donor_sets[[time_position]] <- donors
    covariates <- as.matrix(
      pc_data[match(sample_ids, pc_data$Sample_id), paste0("PC", 1:5),
              drop = FALSE]
    )
    storage.mode(covariates) <- "double"
    rownames(covariates) <- donors
    covariate_list[[time_position]] <- covariates
  }

  message("Screening variants for residual genotype information.")
  screen <- screen_uninformative_variants(
    dosage = dosage,
    permuted_dosage = permuted_dosage,
    donor_sets = donor_sets,
    covariate_list = covariate_list
  )
  write_csv(
    screen$by_time,
    file.path(work_directory, "uninformative_variant_screen.csv")
  )
  n_flagged_variant <- sum(screen$flagged)
  retained_pairs <- !screen$flagged[pair_index$variant_index]
  n_dropped_pair <- sum(!retained_pairs)
  message(
    "  variants dropped: ", format(n_flagged_variant, big.mark = ","),
    " of ", format(length(screen$flagged), big.mark = ","),
    " (", format(100 * n_flagged_variant / length(screen$flagged), digits = 3),
    "%) | pairs dropped: ", format(n_dropped_pair, big.mark = ","),
    " of ", format(length(retained_pairs), big.mark = ","),
    " (", format(100 * n_dropped_pair / length(retained_pairs), digits = 3), "%)"
  )
  message(
    "  observed-arm failures by time point: ",
    paste(screen$by_time$n_uninformative_observed, collapse = " "),
    " | permuted-arm: ",
    paste(screen$by_time$n_uninformative_permuted, collapse = " ")
  )
  if (!any(retained_pairs)) {
    stop("The informativeness screen removed every pair.")
  }
  if (n_dropped_pair / length(retained_pairs) > 0.01) {
    stop(
      "The informativeness screen would drop ",
      format(100 * n_dropped_pair / length(retained_pairs), digits = 3),
      "% of pairs, which is too many to treat as a design technicality. ",
      "Inspect uninformative_variant_screen.csv before continuing."
    )
  }

  if (n_dropped_pair > 0L) {
    pair_keys <- pair_keys[retained_pairs]
    n_pair <- length(pair_keys)
    source_rows <- source_rows[retained_pairs]
    target_lfdr <- target_lfdr[retained_pairs]
    target_likelihood <- target_likelihood[retained_pairs, , drop = FALSE]
    pair_index <- build_pair_index(pair_keys)
    # The screen removes whole variants, so surviving variants keep every
    # informative time point and the variant-level dosage columns are reindexed
    # rather than modified.
    keep_variant <- match(pair_index$variant_levels, colnames(dosage))
    if (anyNA(keep_variant)) {
      stop("A retained variant is missing from the dosage matrix.")
    }
    dosage <- dosage[, keep_variant, drop = FALSE]
    permuted_dosage <- permuted_dosage[, keep_variant, drop = FALSE]
    genotype_sumsq <- genotype_sumsq[keep_variant]
    gene_rows <- match(pair_index$gene_levels, expression_data$Gene_id)
    if (anyNA(gene_rows)) {
      stop("A retained gene is missing from the expression matrix.")
    }
    message(
      "  retained: ", format(n_pair, big.mark = ","), " pairs, ",
      format(length(pair_index$gene_levels), big.mark = ","), " genes, ",
      format(length(pair_index$variant_levels), big.mark = ","), " variants"
    )
  }

  message("Caching the immutable target likelihood rows.")
  saveRDS(
    list(
      cache_version = "fash_knockoff_target_likelihood_v1",
      seed = seed,
      gene_limit = gene_limit,
      source_rows = source_rows,
      pair_key = pair_keys,
      L_matrix = target_likelihood,
      psd_grid = original_psd_grid,
      settings = original_settings,
      bf_lfdr = target_lfdr,
      source_fit_object = loaded_names,
      n_dropped_uninformative_pair = n_dropped_pair
    ),
    target_likelihood_path,
    compress = FALSE
  )

  # Remap both validation slices onto the screened pair set.
  surviving_pool <- which(check_pool_keys %in% pair_keys)
  likelihood_check_rows <- utils::head(surviving_pool, 50L)
  cached_likelihood_check <- check_pool_likelihood[likelihood_check_rows, ,
                                                  drop = FALSE]
  check_beta <- check_pool_beta[likelihood_check_rows, , drop = FALSE]
  check_adjusted_se <- check_pool_adjusted_se[likelihood_check_rows, ,
                                              drop = FALSE]
  validation_rows <- match(validation_keys, pair_keys)
  kept_validation <- !is.na(validation_rows)
  if (any(!kept_validation)) {
    message(
      "  ", sum(!kept_validation), " of ", length(kept_validation),
      " validation pairs were screened out; the check uses the remainder."
    )
    validation_rows <- validation_rows[kept_validation]
    source_beta <- source_beta[kept_validation, , drop = FALSE]
    source_adjusted_se <- source_adjusted_se[kept_validation, , drop = FALSE]
  }
  rm(target_likelihood, check_pool_likelihood, check_pool_beta,
     check_pool_adjusted_se)
  invisible(gc(verbose = FALSE))

  cached_likelihood_deviation <- NA_real_
  if (length(likelihood_check_rows) > 0L) {
    message(
      "Reproducing ", length(likelihood_check_rows),
      " cached target likelihood rows to confirm numerical agreement."
    )
    check_capture <- capture_warnings(fashr::fash(
      Y = "beta", smooth_var = "time", S = "SE",
      data_list = lapply(seq_along(likelihood_check_rows), function(position) {
        data.frame(
          time = time_grid,
          beta = check_beta[position, ],
          SE = check_adjusted_se[position, ]
        )
      }),
      num_basis = original_settings$num_basis,
      order = original_settings$order,
      betaprec = original_settings$betaprec,
      pred_step = original_settings$pred_step,
      penalty = original_settings$penalty,
      grid = original_psd_grid,
      num_cores = 1L,
      verbose = FALSE
    ))
    cached_likelihood_deviation <- max(abs(
      unname(check_capture$value$L_matrix) - unname(cached_likelihood_check)
    ))
    message(
      "  max |recomputed - cached| log-likelihood: ",
      format(cached_likelihood_deviation, digits = 3),
      " (scale ", format(max(abs(cached_likelihood_check)), digits = 3), ")"
    )
    if (!is.finite(cached_likelihood_deviation) ||
        cached_likelihood_deviation > 1e-8) {
      stop(
        "This environment does not reproduce the cached target likelihood ",
        "rows. Target and decoy likelihoods would not be comparable; refusing ",
        "to continue."
      )
    }
    rm(check_beta, check_adjusted_se, cached_likelihood_check)
  }


  message("Refitting the time-specific PC model for the permuted genotype.")
  decoy_beta <- matrix(
    NA_real_, n_pair, length(time_grid),
    dimnames = list(pair_keys, paste0("time_", time_grid))
  )
  decoy_raw_se <- decoy_beta
  residual_df_by_time <- integer(length(time_grid))
  validation_beta <- if (length(validation_rows) > 0L) {
    matrix(NA_real_, length(validation_rows), length(time_grid))
  } else {
    NULL
  }
  validation_raw_se <- validation_beta
  alignment_tables <- vector("list", length(time_grid))

  for (time_position in seq_along(time_grid)) {
    time_value <- time_grid[time_position]
    sample_ids <- grep(
      paste0("_", time_value, "$"), expression_sample_ids, value = TRUE
    )
    donors <- donor_sets[[time_position]]
    covariates <- covariate_list[[time_position]]
    expression_matrix <- t(as.matrix(
      expression_data[gene_rows, match(sample_ids, names(expression_data)),
                      drop = FALSE]
    ))
    storage.mode(expression_matrix) <- "double"

    # The PCs are functions of the expression matrix and are therefore NOT
    # permuted; only the genotype rows move.
    projection <- make_covariate_projection(covariates)
    residual_df <- nrow(expression_matrix) - projection$rank - 1L
    residual_df_by_time[time_position] <- residual_df
    expression_residual <- projection$residualizer %*% expression_matrix
    observed_genotype <- dosage[donors, , drop = FALSE]
    decoy_genotype <- permuted_dosage[donors, , drop = FALSE]

    alignment_tables[[time_position]] <- data.frame(
      time = time_value,
      summarize_genotype_covariate_alignment(
        observed_genotype, decoy_genotype, covariates
      ),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    decoy_fit <- fit_pairwise_genotype_regressions(
      expression_residual = expression_residual,
      genotype_residual = projection$residualizer %*% decoy_genotype,
      genotype_sumsq = colSums(decoy_genotype^2),
      gene_index = pair_index$gene_index,
      variant_index = pair_index$variant_index,
      residual_df = residual_df
    )
    decoy_beta[, time_position] <- decoy_fit$beta
    decoy_raw_se[, time_position] <- decoy_fit$standard_error

    if (length(validation_rows) > 0L) {
      observed_fit <- fit_pairwise_genotype_regressions(
        expression_residual = expression_residual,
        genotype_residual = projection$residualizer %*% observed_genotype,
        genotype_sumsq = genotype_sumsq,
        gene_index = pair_index$gene_index[validation_rows],
        variant_index = pair_index$variant_index[validation_rows],
        residual_df = residual_df
      )
      validation_beta[, time_position] <- observed_fit$beta
      validation_raw_se[, time_position] <- observed_fit$standard_error
    }
    message(
      "  time ", time_value, ": n = ", nrow(expression_matrix),
      ", residual df = ", residual_df
    )
  }

  if (!identical(residual_df_by_time, expected_sample_counts - 7L)) {
    stop("The residual degrees of freedom do not match the primary analysis.")
  }

  message("Converting raw standard errors to the original t-adjusted scale.")
  decoy_adjusted_se <- convert_raw_to_original_t_adjusted_se(
    decoy_beta, decoy_raw_se, residual_df_by_time
  )

  # Reproduce the published beta/SE from raw inputs on the validation slice.
  beta_deviation <- NA_real_
  se_deviation <- NA_real_
  if (length(validation_rows) > 0L) {
    validation_adjusted_se <- convert_raw_to_original_t_adjusted_se(
      validation_beta, validation_raw_se, residual_df_by_time
    )
    beta_deviation <- max(abs(validation_beta - source_beta))
    se_deviation <- max(abs(validation_adjusted_se - source_adjusted_se))
    message(
      "Validation slice (", length(validation_rows), " pairs): max |beta| ",
      "deviation ", format(beta_deviation, digits = 3),
      ", max |SE| deviation ", format(se_deviation, digits = 3)
    )
    if (!is.finite(beta_deviation) || beta_deviation > 1e-8 ||
        !is.finite(se_deviation) || se_deviation > 1e-6) {
      stop(
        "The regression pipeline does not reproduce the published beta/SE ",
        "values; refusing to build decoy units."
      )
    }
  }

  alignment_summary <- do.call(rbind, alignment_tables)
  write_csv(
    alignment_summary,
    file.path(work_directory, "genotype_covariate_alignment.csv")
  )
  write_csv(donor_map, file.path(work_directory, "donor_permutation.csv"))
  original_alignment <- alignment_summary$mean_r_squared[
    alignment_summary$arm == "original"
  ]
  permuted_alignment <- alignment_summary$mean_r_squared[
    alignment_summary$arm == "permuted"
  ]
  message(
    "Genotype-PC alignment gate: mean R^2 original ",
    format(mean(original_alignment), digits = 4), " vs permuted ",
    format(mean(permuted_alignment), digits = 4),
    " (a systematically larger original value deflates decoy SEs and makes ",
    "the competition anti-conservative)."
  )

  saveRDS(
    list(
      cache_version = "fash_knockoff_decoy_input_v1",
      seed = seed,
      gene_limit = gene_limit,
      pair_key = pair_keys,
      gene_id = pair_index$gene_id,
      variant_id = pair_index$variant_id,
      decoy_beta = decoy_beta,
      decoy_adjusted_se = decoy_adjusted_se,
      residual_df_by_time = residual_df_by_time,
      donor_map = donor_map,
      permutation_seed = permutation$permutation_seed,
      permutation_attempts = permutation$attempts,
      psd_grid = original_psd_grid,
      settings = original_settings,
      time_grid = time_grid,
      validation_rows = validation_rows,
      validation_beta_deviation = beta_deviation,
      validation_se_deviation = se_deviation,
      cached_likelihood_deviation = cached_likelihood_deviation,
      n_dropped_uninformative_pair = n_dropped_pair,
      uninformative_screen_by_time = screen$by_time,
      genotype_covariate_alignment = alignment_summary,
      source_information = source_information_before
    ),
    prepared_path,
    compress = FALSE
  )
  write_csv(
    file_metadata(required_paths, required_roles),
    file.path(work_directory, "source_information_prepare.csv")
  )
  message(
    "Prepare stage complete in ",
    format(proc.time()[["elapsed"]] - stage_start, digits = 5), " seconds."
  )
  message("Next: rerun with --stage likelihood in a fresh R process.")
  quit(save = "no", status = 0L)
}

# ===========================================================================
# Stage 2: chunked, resumable decoy likelihood evaluation
# ===========================================================================

if (identical(stage, "likelihood")) {
  stage_start <- proc.time()[["elapsed"]]
  if (!file.exists(prepared_path)) {
    stop("The prepared decoy input cache is missing; run --stage prepare first.")
  }
  prepared <- readRDS(prepared_path)
  if (!identical(prepared$cache_version, "fash_knockoff_decoy_input_v1") ||
      !identical(prepared$seed, seed) ||
      !identical(prepared$gene_limit, gene_limit)) {
    stop("The prepared decoy input cache does not match this run.")
  }
  decoy_beta <- prepared$decoy_beta
  decoy_adjusted_se <- prepared$decoy_adjusted_se
  n_unit <- nrow(decoy_beta)
  if (n_unit != length(prepared$pair_key) || n_unit < 1L ||
      !identical(dim(decoy_beta), dim(decoy_adjusted_se)) ||
      any(!is.finite(decoy_beta)) || any(!is.finite(decoy_adjusted_se)) ||
      any(decoy_adjusted_se <= 0)) {
    stop("The prepared decoy beta/SE matrices are invalid.")
  }
  original_psd_grid <- prepared$psd_grid
  original_settings <- prepared$settings
  dir.create(chunk_directory, recursive = TRUE, showWarnings = FALSE)

  chunk_starts <- seq.int(1L, n_unit, by = chunk_units)
  message(
    "Evaluating ", format(n_unit, big.mark = ","), " decoy units in ",
    length(chunk_starts), " chunks of at most ", chunk_units,
    " on ", num_cores, " cores."
  )
  all_warnings <- character()
  for (chunk_position in seq_along(chunk_starts)) {
    start <- chunk_starts[chunk_position]
    end <- min(start + chunk_units - 1L, n_unit)
    chunk_path <- file.path(
      chunk_directory, sprintf("chunk_%06d.rds", chunk_position)
    )
    if (file.exists(chunk_path)) {
      existing <- readRDS(chunk_path)
      if (identical(existing$cache_version, "fash_knockoff_chunk_v1") &&
          identical(existing$row_start, start) &&
          identical(existing$row_end, end) &&
          identical(existing$seed, seed) &&
          nrow(existing$L_matrix) == end - start + 1L) {
        next
      }
      stop("Chunk ", chunk_position, " exists but does not match this run.")
    }
    rows <- start:end
    chunk_data <- lapply(rows, function(row) {
      data.frame(
        time = time_grid,
        beta = decoy_beta[row, ],
        SE = decoy_adjusted_se[row, ]
      )
    })
    chunk_start_time <- proc.time()[["elapsed"]]
    chunk_capture <- capture_warnings(fashr::fash(
      Y = "beta", smooth_var = "time", S = "SE",
      data_list = chunk_data,
      num_basis = original_settings$num_basis,
      order = original_settings$order,
      betaprec = original_settings$betaprec,
      pred_step = original_settings$pred_step,
      penalty = original_settings$penalty,
      grid = original_psd_grid,
      num_cores = num_cores,
      verbose = FALSE
    ))
    chunk_fit <- chunk_capture$value
    chunk_likelihood <- chunk_fit$L_matrix
    if (!isTRUE(all.equal(chunk_fit$psd_grid, original_psd_grid,
                          tolerance = 0)) ||
        nrow(chunk_likelihood) != length(rows) ||
        ncol(chunk_likelihood) != expected_grid_size ||
        anyNA(chunk_likelihood) || any(chunk_likelihood == Inf)) {
      stop("Chunk ", chunk_position, " produced an invalid likelihood block.")
    }
    rownames(chunk_likelihood) <- NULL
    all_warnings <- unique(c(all_warnings, chunk_capture$warnings))
    saveRDS(
      list(
        cache_version = "fash_knockoff_chunk_v1",
        seed = seed,
        row_start = start,
        row_end = end,
        L_matrix = chunk_likelihood,
        warnings = chunk_capture$warnings,
        elapsed_seconds = proc.time()[["elapsed"]] - chunk_start_time
      ),
      chunk_path,
      compress = FALSE
    )
    elapsed <- proc.time()[["elapsed"]] - stage_start
    message(
      "  chunk ", chunk_position, "/", length(chunk_starts),
      " rows ", start, "-", end,
      " in ", format(proc.time()[["elapsed"]] - chunk_start_time, digits = 4),
      "s | elapsed ", format(elapsed / 60, digits = 4), " min",
      " | projected total ",
      format(elapsed / chunk_position * length(chunk_starts) / 3600,
             digits = 4),
      " h"
    )
  }
  if (length(all_warnings) > 0L) {
    writeLines(
      all_warnings, file.path(work_directory, "decoy_likelihood_warnings.txt")
    )
  }
  message(
    "Likelihood stage complete in ",
    format((proc.time()[["elapsed"]] - stage_start) / 60, digits = 5),
    " minutes."
  )
  message("Next: rerun with --stage merge in a fresh R process.")
  quit(save = "no", status = 0L)
}

# ===========================================================================
# Stage 3: merged empirical-Bayes refit and knockoff+ selection
# ===========================================================================

stage_start <- proc.time()[["elapsed"]]
if (!file.exists(prepared_path) || !file.exists(target_likelihood_path)) {
  stop("The prepared caches are missing; run the earlier stages first.")
}
message("Loading cached target likelihood rows.")
target_cache <- readRDS(target_likelihood_path)
if (!identical(target_cache$cache_version,
               "fash_knockoff_target_likelihood_v1") ||
    !identical(target_cache$seed, seed) ||
    !identical(target_cache$gene_limit, gene_limit) ||
    nrow(target_cache$L_matrix) != length(target_cache$pair_key)) {
  stop("The target likelihood cache is invalid.")
}
prepared <- readRDS(prepared_path)
if (!identical(prepared$cache_version, "fash_knockoff_decoy_input_v1") ||
    !identical(prepared$seed, seed) ||
    !identical(prepared$gene_limit, gene_limit) ||
    !identical(prepared$pair_key, target_cache$pair_key)) {
  stop("The prepared decoy cache is not aligned with the target cache.")
}
pair_keys <- target_cache$pair_key
n_pair <- length(pair_keys)
original_psd_grid <- target_cache$psd_grid
original_settings <- target_cache$settings
target_bf_lfdr <- target_cache$bf_lfdr

message("Assembling the decoy likelihood matrix from cached chunks.")
chunk_paths <- sort(list.files(
  chunk_directory, pattern = "^chunk_[0-9]{6}\\.rds$", full.names = TRUE
))
if (length(chunk_paths) == 0L) {
  stop("No decoy likelihood chunks are available.")
}
decoy_likelihood <- matrix(NA_real_, n_pair, expected_grid_size)
covered <- logical(n_pair)
for (chunk_path in chunk_paths) {
  chunk <- readRDS(chunk_path)
  if (!identical(chunk$cache_version, "fash_knockoff_chunk_v1") ||
      !identical(chunk$seed, seed) ||
      nrow(chunk$L_matrix) != chunk$row_end - chunk$row_start + 1L ||
      ncol(chunk$L_matrix) != expected_grid_size) {
    stop("Chunk ", basename(chunk_path), " is invalid.")
  }
  rows <- chunk$row_start:chunk$row_end
  if (any(covered[rows])) {
    stop("Chunk ", basename(chunk_path), " overlaps an earlier chunk.")
  }
  decoy_likelihood[rows, ] <- chunk$L_matrix
  covered[rows] <- TRUE
}
if (!all(covered) || anyNA(decoy_likelihood) ||
    any(decoy_likelihood == Inf)) {
  stop(
    "The decoy likelihood matrix is incomplete; ",
    sum(!covered), " unit(s) are missing."
  )
}
rm(covered)
invisible(gc(verbose = FALSE))

decoy_keys <- paste0(pair_keys, "__decoy_seed", seed)
unit_keys <- c(pair_keys, decoy_keys)
merged_likelihood <- rbind(target_cache$L_matrix, decoy_likelihood)
rm(decoy_likelihood, target_cache)
invisible(gc(verbose = FALSE))

message(
  "Refitting the merged empirical-Bayes mixture on ",
  format(nrow(merged_likelihood), big.mark = ","), " units."
)
raw_start <- proc.time()[["elapsed"]]
raw_capture <- capture_warnings(refit_merged_fash_from_likelihood(
  likelihood_matrix = merged_likelihood,
  psd_grid = original_psd_grid,
  settings = original_settings,
  unit_keys = unit_keys,
  penalty = original_settings$penalty
))
merged_raw <- raw_capture$value
rm(merged_likelihood)
invisible(gc(verbose = FALSE))
message(
  "  merged raw refit in ",
  format(proc.time()[["elapsed"]] - raw_start, digits = 5), " seconds."
)

message("Applying the BF prior update to the merged fit.")
bf_capture <- capture_warnings(fashr::BF_update(merged_raw, plot = FALSE))
merged_bf <- bf_capture$value
names(merged_bf$lfdr) <- unit_keys

# Two admissible scoring functions. Both are computed from the merged fit, so
# both flip sign when a target-decoy pair is swapped; they differ only in
# which alternative weights normalize the collapsed Bayes factor.
raw_bayes_factor <- as.numeric(fashr::BF_compute(merged_raw))
bf_bayes_factor <- as.numeric(fashr::BF_compute(merged_bf))
target_rows <- seq_len(n_pair)
decoy_rows <- n_pair + target_rows

statistics <- data.frame(
  pair_key = pair_keys,
  gene_id = prepared$gene_id,
  variant_id = prepared$variant_id,
  source_bf_lfdr = target_bf_lfdr,
  target_bf_raw_weights = raw_bayes_factor[target_rows],
  decoy_bf_raw_weights = raw_bayes_factor[decoy_rows],
  target_bf_updated_weights = bf_bayes_factor[target_rows],
  decoy_bf_updated_weights = bf_bayes_factor[decoy_rows],
  target_lfdr = as.numeric(merged_bf$lfdr)[target_rows],
  decoy_lfdr = as.numeric(merged_bf$lfdr)[decoy_rows],
  stringsAsFactors = FALSE
)
statistics$W_raw_weights <- compute_competition_statistic(
  statistics$target_bf_raw_weights, statistics$decoy_bf_raw_weights
)
statistics$W_updated_weights <- compute_competition_statistic(
  statistics$target_bf_updated_weights, statistics$decoy_bf_updated_weights
)

# BF_update rescales the alternative weights by a common factor, so the
# collapsed two-component Bayes factor -- and therefore the competition
# statistic -- is unchanged by it. lfdr, by contrast, moves with pi0. This is
# checked rather than assumed, because it is the reason the BF statistic
# sidesteps the BF_update calibration question entirely.
bf_weight_deviation <- max(abs(
  statistics$W_raw_weights - statistics$W_updated_weights
))
message(
  "Competition statistic under raw versus BF-updated alternative weights: ",
  "max |difference| ", format(bf_weight_deviation, digits = 3), "."
)
if (!is.finite(bf_weight_deviation) || bf_weight_deviation > 1e-8) {
  warning(
    "The raw and BF-updated competition statistics disagree by more than ",
    "1e-8; report both scoring functions rather than treating them as one."
  )
}

# Gene-level competition: the maximum over the same variant set on both
# sides, so a gene's target maximum is not compared against a single decoy.
gene_raw <- aggregate_gene_competition(
  log(statistics$target_bf_raw_weights),
  log(statistics$decoy_bf_raw_weights),
  statistics$gene_id
)
gene_updated <- aggregate_gene_competition(
  log(statistics$target_bf_updated_weights),
  log(statistics$decoy_bf_updated_weights),
  statistics$gene_id
)

selection_tables <- rbind(
  data.frame(
    unit = "pair", weights = "raw",
    knockoff_selection_summary(statistics$W_raw_weights, q_grid),
    stringsAsFactors = FALSE
  ),
  data.frame(
    unit = "pair", weights = "bf_updated",
    knockoff_selection_summary(statistics$W_updated_weights, q_grid),
    stringsAsFactors = FALSE
  ),
  data.frame(
    unit = "gene", weights = "raw",
    knockoff_selection_summary(gene_raw$W, q_grid),
    stringsAsFactors = FALSE
  ),
  data.frame(
    unit = "gene", weights = "bf_updated",
    knockoff_selection_summary(gene_updated$W, q_grid),
    stringsAsFactors = FALSE
  )
)

sign_summary <- data.frame(
  unit = c("pair", "pair", "gene", "gene"),
  weights = c("raw", "bf_updated", "raw", "bf_updated"),
  n_unit = c(n_pair, n_pair, nrow(gene_raw), nrow(gene_updated)),
  n_target_wins = c(
    sum(statistics$W_raw_weights > 0), sum(statistics$W_updated_weights > 0),
    sum(gene_raw$W > 0), sum(gene_updated$W > 0)
  ),
  n_decoy_wins = c(
    sum(statistics$W_raw_weights < 0), sum(statistics$W_updated_weights < 0),
    sum(gene_raw$W < 0), sum(gene_updated$W < 0)
  ),
  minimum_estimated_fdr = c(
    min(knockoff_estimated_fdr_path(statistics$W_raw_weights)$estimated_fdr),
    min(knockoff_estimated_fdr_path(
      statistics$W_updated_weights
    )$estimated_fdr),
    min(knockoff_estimated_fdr_path(gene_raw$W)$estimated_fdr),
    min(knockoff_estimated_fdr_path(gene_updated$W)$estimated_fdr)
  ),
  stringsAsFactors = FALSE
)

message("Writing outputs.")
saveRDS(
  list(
    cache_version = "fash_knockoff_merged_fit_v1",
    seed = seed,
    gene_limit = gene_limit,
    prior_weights_raw = merged_raw$prior_weights,
    prior_weights_bf = merged_bf$prior_weights,
    merged_pi0_raw = extract_prior_pi0(merged_raw),
    merged_pi0_bf = extract_prior_pi0(merged_bf),
    psd_grid = original_psd_grid,
    settings = original_settings,
    unit_keys = unit_keys,
    warnings = unique(c(raw_capture$warnings, bf_capture$warnings))
  ),
  file.path(work_directory, "merged_fash_fit_summary.rds"),
  compress = TRUE
)
saveRDS(
  statistics, file.path(work_directory, "pair_statistics.rds"),
  compress = TRUE
)
write_csv(gene_raw, file.path(work_directory, "gene_statistics_raw.csv"))
write_csv(
  gene_updated, file.path(work_directory, "gene_statistics_bf_updated.csv")
)
write_csv(
  selection_tables, file.path(work_directory, "knockoff_selection.csv")
)
write_csv(sign_summary, file.path(work_directory, "competition_summary.csv"))
write_csv(
  knockoff_estimated_fdr_path(gene_updated$W),
  file.path(work_directory, "gene_estimated_fdr_path_bf_updated.csv")
)
write_csv(
  data.frame(
    quantity = c("merged_pi0_raw", "merged_pi0_bf", "n_pair", "n_gene"),
    value = c(
      extract_prior_pi0(merged_raw), extract_prior_pi0(merged_bf),
      n_pair, nrow(gene_updated)
    ),
    stringsAsFactors = FALSE
  ),
  file.path(work_directory, "merged_fit_diagnostics.csv")
)
write_csv(
  file_metadata(required_paths, required_roles),
  file.path(work_directory, "source_information_merge.csv")
)

message("")
message("Merged prior pi0: raw ", format(extract_prior_pi0(merged_raw),
        digits = 5), " | BF-adjusted ",
        format(extract_prior_pi0(merged_bf), digits = 5))
print(sign_summary)
message("")
print(selection_tables[selection_tables$unit == "gene" &
                         selection_tables$weights == "bf_updated", ])
if (nzchar(summary_directory)) {
  dir.create(summary_directory, recursive = TRUE, showWarnings = FALSE)
  small_tables <- list.files(
    work_directory, pattern = "\\.csv$", full.names = TRUE
  )
  copied <- file.copy(small_tables, summary_directory, overwrite = TRUE)
  message(
    "Copied ", sum(copied), " of ", length(small_tables),
    " diagnostic table(s) to ", summary_directory
  )
}

writeLines(
  c(
    paste0("seed=", seed),
    paste0("gene_limit=", gene_limit),
    paste0("n_pair=", n_pair),
    paste0("n_gene=", nrow(gene_updated)),
    paste0("fashr=", observed_fashr_version),
    paste0("fashr_remote_sha=", observed_fashr_remote_sha),
    paste0("r_version=", as.character(getRversion()))
  ),
  file.path(work_directory, "complete.flag")
)
message(
  "Merge stage complete in ",
  format((proc.time()[["elapsed"]] - stage_start) / 60, digits = 5),
  " minutes."
)
message("Wrote ", file.path(work_directory, "complete.flag"))
