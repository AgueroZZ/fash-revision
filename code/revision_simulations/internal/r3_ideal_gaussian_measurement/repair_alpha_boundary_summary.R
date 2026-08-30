#!/usr/bin/env Rscript

# Build an immutable derived cache that restores the nominal alpha = 0.05 row
# omitted by a strict floating-point comparison in the production summary.

select_alpha_interval <- function(alpha,
                                  lower,
                                  upper,
                                  tolerance = 1e-12) {
  if (length(lower) != 1L || length(upper) != 1L ||
      !is.finite(lower) || !is.finite(upper) || lower > upper) {
    stop("lower and upper must be ordered finite scalars.", call. = FALSE)
  }
  if (length(tolerance) != 1L || !is.finite(tolerance) || tolerance < 0) {
    stop("tolerance must be one non-negative finite value.", call. = FALSE)
  }
  is.finite(alpha) &
    alpha >= lower - tolerance &
    alpha <= upper + tolerance
}

canonicalize_alpha <- function(alpha, digits = 12L) {
  if (length(digits) != 1L || !is.finite(digits) ||
      digits != round(digits) || digits < 1L) {
    stop("digits must be one positive integer.", call. = FALSE)
  }
  round(as.numeric(alpha), digits = as.integer(digits))
}

rebuild_middle_summaries <- function(all_alpha,
                                     mechanisms,
                                     expected_seeds,
                                     lower = 0.05,
                                     upper = 0.20,
                                     tolerance = 1e-12) {
  required <- c(
    "method", "target", "alpha", "empirical_fsr", "power", "seed",
    "truth_mechanism"
  )
  if (!is.data.frame(all_alpha) || !all(required %in% names(all_alpha))) {
    stop("all_alpha does not contain the required columns.", call. = FALSE)
  }
  if (length(mechanisms) == 0L || any(!nzchar(mechanisms)) ||
      anyDuplicated(mechanisms)) {
    stop("mechanisms must contain unique non-empty values.", call. = FALSE)
  }
  if (length(expected_seeds) == 0L || anyNA(expected_seeds) ||
      anyDuplicated(expected_seeds)) {
    stop("expected_seeds must contain unique values.", call. = FALSE)
  }
  rows <- all_alpha[
    all_alpha$method == "FASH-IWP1-BF" &
      all_alpha$target == "middle" &
      all_alpha$truth_mechanism %in% mechanisms &
      select_alpha_interval(all_alpha$alpha, lower, upper, tolerance),
    ,
    drop = FALSE
  ]
  if (nrow(rows) == 0L ||
      any(!is.finite(rows$empirical_fsr)) ||
      any(!is.finite(rows$power))) {
    stop("No valid BF-adjusted Middle rows remain after filtering.", call. = FALSE)
  }
  rows$alpha <- canonicalize_alpha(rows$alpha)
  rows$alpha_key <- sprintf("%.12f", rows$alpha)
  group_key <- interaction(
    rows$truth_mechanism,
    rows$alpha_key,
    drop = TRUE,
    lex.order = TRUE
  )
  curve <- do.call(rbind, lapply(split(rows, group_key), function(group) {
    if (!setequal(group$seed, expected_seeds) ||
        nrow(group) != length(expected_seeds) ||
        length(unique(group$alpha)) != 1L ||
        length(unique(group$truth_mechanism)) != 1L) {
      stop("A Middle curve point does not contain one row per expected seed.",
           call. = FALSE)
    }
    data.frame(
      truth_mechanism = group$truth_mechanism[[1L]],
      alpha = group$alpha[[1L]],
      mean_empirical_fsr = mean(group$empirical_fsr),
      min_empirical_fsr = min(group$empirical_fsr),
      max_empirical_fsr = max(group$empirical_fsr),
      mean_power = mean(group$power),
      stringsAsFactors = FALSE
    )
  }))
  curve <- curve[order(
    match(curve$truth_mechanism, mechanisms),
    curve$alpha
  ), ]
  rownames(curve) <- NULL
  if (!setequal(curve$truth_mechanism, mechanisms)) {
    stop("The rebuilt curve is missing a truth mechanism.", call. = FALSE)
  }

  primary <- do.call(rbind, lapply(
    split(curve, curve$truth_mechanism),
    function(group) {
      excess <- group$mean_empirical_fsr - group$alpha
      index <- which.max(excess)
      data.frame(
        truth_mechanism = group$truth_mechanism[[1L]],
        method = "FASH-IWP1-BF",
        target = "middle",
        alpha_min = min(group$alpha),
        alpha_max = max(group$alpha),
        maximum_mean_fsr_excess = excess[[index]],
        alpha_at_maximum = group$alpha[[index]],
        mean_empirical_fsr_at_maximum =
          group$mean_empirical_fsr[[index]],
        interpretation_threshold = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  ))
  primary <- primary[match(mechanisms, primary$truth_mechanism), ]
  rownames(primary) <- NULL
  list(curve = curve, primary = primary)
}

get_arg <- function(name, default = NULL) {
  arguments <- commandArgs(trailingOnly = TRUE)
  index <- which(arguments == name)
  if (length(index) == 0L || index[[1L]] == length(arguments)) return(default)
  arguments[[index[[1L]] + 1L]]
}

as_flag <- function(value) {
  tolower(value) %in% c("1", "true", "t", "yes", "y")
}

require_file <- function(path, label) {
  if (!file.exists(path) || dir.exists(path)) {
    stop(label, " is missing: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

sha256_file <- function(path) {
  path <- require_file(path, "SHA-256 input")
  command <- if (identical(Sys.info()[["sysname"]], "Darwin")) {
    "shasum"
  } else {
    "sha256sum"
  }
  arguments <- if (identical(command, "shasum")) c("-a", "256", path) else path
  output <- system2(command, arguments, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Unable to compute SHA-256 for ", path, ".", call. = FALSE)
  }
  sub("[[:space:]].*$", "", output[[1L]])
}

relative_file_inventory <- function(directory) {
  directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  paths <- list.files(
    directory,
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE,
    include.dirs = FALSE
  )
  paths <- sort(paths[file.info(paths)$isdir %in% FALSE])
  relative <- substring(paths, nchar(directory) + 2L)
  stats::setNames(vapply(paths, sha256_file, character(1)), relative)
}

copy_directory_tree <- function(source, destination) {
  source <- normalizePath(source, winslash = "/", mustWork = TRUE)
  if (file.exists(destination) || dir.exists(destination)) {
    stop("Destination already exists: ", destination, call. = FALSE)
  }
  if (!dir.create(destination, recursive = TRUE)) {
    stop("Unable to create destination: ", destination, call. = FALSE)
  }
  source_directories <- list.dirs(
    source, recursive = TRUE, full.names = TRUE
  )
  relative_directories <- substring(
    source_directories, nchar(source) + 2L
  )
  relative_directories <- relative_directories[nzchar(relative_directories)]
  for (relative in relative_directories) {
    path <- file.path(destination, relative)
    if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) {
      stop("Unable to create copied directory: ", path, call. = FALSE)
    }
  }
  source_files <- list.files(
    source,
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE,
    include.dirs = FALSE
  )
  relative_files <- substring(source_files, nchar(source) + 2L)
  destination_files <- file.path(destination, relative_files)
  copied <- file.copy(
    from = source_files,
    to = destination_files,
    overwrite = FALSE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
  if (length(copied) != length(source_files) || any(!copied)) {
    stop("One or more parent artifacts could not be copied.", call. = FALSE)
  }
  invisible(destination)
}

atomic_save_rds <- function(object, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3)
  if (!file.rename(temporary, path)) {
    stop("Unable to atomically replace ", path, ".", call. = FALSE)
  }
  invisible(path)
}

atomic_write_csv <- function(object, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(object, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) {
    stop("Unable to atomically replace ", path, ".", call. = FALSE)
  }
  invisible(path)
}

atomic_write_lines <- function(text, path) {
  temporary <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary), add = TRUE)
  writeLines(text, temporary, useBytes = TRUE)
  if (!file.rename(temporary, path)) {
    stop("Unable to atomically replace ", path, ".", call. = FALSE)
  }
  invisible(path)
}

repair_alpha_boundary_summary <- function(input_dir,
                                          output_dir,
                                          script_path,
                                          preflight_only = FALSE) {
  parent_result_id <- paste0(
    "r3_ideal_gaussian_known_t_adjusted_se_",
    "matched_truth_open_middle_3_12_full_universe_",
    "fashr0143_pilot5"
  )
  derived_result_id <- paste0(parent_result_id, "_summaryfix1")
  parent_schema <- "r3-ideal-gaussian-measurement-v1"
  derived_schema <- "r3-ideal-gaussian-measurement-summaryfix-v1"
  expected_seeds <- c(12345L, 22345L, 32345L, 42345L, 52345L)
  mechanisms <- c("random_bspline", "raised_cosine")
  methods <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
  targets <- c("early", "middle", "late", "switch")
  expected_alpha_grid <- canonicalize_alpha(seq(0.005, 0.20, by = 0.005))

  input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
  output_dir <- normalizePath(
    dirname(output_dir), winslash = "/", mustWork = TRUE
  ) |>
    file.path(basename(output_dir))
  script_path <- require_file(script_path, "Repair script")
  if (!identical(basename(input_dir), parent_result_id)) {
    stop("The input directory does not have the expected parent result ID.",
         call. = FALSE)
  }
  if (!identical(basename(output_dir), derived_result_id)) {
    stop("The output directory does not have the expected derived result ID.",
         call. = FALSE)
  }
  if (file.exists(output_dir) || dir.exists(output_dir)) {
    stop("Refusing to overwrite derived output: ", output_dir, call. = FALSE)
  }
  temporary_dir <- paste0(output_dir, "_partial")
  if (file.exists(temporary_dir) || dir.exists(temporary_dir)) {
    stop("A derived partial directory already exists: ", temporary_dir,
         call. = FALSE)
  }

  parent_manifest_path <- require_file(
    file.path(input_dir, "manifest.rds"), "Parent manifest"
  )
  parent_configuration_path <- require_file(
    file.path(input_dir, "configuration.rds"), "Parent configuration"
  )
  parent_complete_path <- require_file(
    file.path(input_dir, "complete.flag"), "Parent completion flag"
  )
  parent_manifest <- readRDS(parent_manifest_path)
  parent_configuration <- readRDS(parent_configuration_path)
  parent_completion <- readLines(parent_complete_path, warn = FALSE)
  if (!identical(parent_manifest$schema_version, parent_schema) ||
      !identical(parent_manifest$result_id, parent_result_id) ||
      !identical(parent_configuration$output_id, parent_result_id) ||
      !identical(parent_configuration$seed_list, expected_seeds) ||
      !identical(parent_configuration$truth_mechanisms, mechanisms) ||
      !all(c(
        paste0("result_id=", parent_result_id),
        paste0("schema_version=", parent_schema),
        "replicates=10"
      ) %in% parent_completion)) {
    stop("The parent cache manifest, configuration, or flag is invalid.",
         call. = FALSE)
  }

  parent_inventory_before <- relative_file_inventory(input_dir)
  manifest_artifacts <- parent_manifest$artifact_sha256
  observed_manifest_artifacts <- parent_inventory_before[
    names(manifest_artifacts)
  ]
  if (anyNA(observed_manifest_artifacts) ||
      !identical(unname(observed_manifest_artifacts),
                 unname(manifest_artifacts))) {
    stop("The parent cache failed its recorded artifact hashes.", call. = FALSE)
  }

  replicate_paths <- unlist(lapply(mechanisms, function(mechanism) {
    file.path(
      input_dir,
      "replicates",
      paste0(mechanism, "_seed_", expected_seeds, ".rds")
    )
  }), use.names = FALSE)
  if (any(!file.exists(replicate_paths))) {
    stop("The parent cache does not contain ten replicate files.", call. = FALSE)
  }
  replicates <- lapply(replicate_paths, readRDS)
  observed_keys <- data.frame(
    seed = vapply(replicates, `[[`, integer(1), "seed"),
    truth_mechanism = vapply(
      replicates, `[[`, character(1), "truth_mechanism"
    ),
    stringsAsFactors = FALSE
  )
  if (!identical(observed_keys, parent_manifest$replicate_keys)) {
    stop("The parent replicate keys do not match its manifest.", call. = FALSE)
  }
  for (replicate in replicates) {
    alpha <- replicate$functional_alpha
    if (nrow(alpha) !=
          length(methods) * length(targets) * length(expected_alpha_grid) ||
        nrow(replicate$functional_alpha_005) !=
          length(methods) * length(targets) ||
        any(alpha$candidate_scope != "full_universe") ||
        any(alpha$candidate_count != 6362L) ||
        any(alpha$first_stage_null_calls != 0L) ||
        any(!is.finite(alpha$empirical_fsr)) ||
        any(!is.finite(alpha$power)) ||
        !identical(
          sort(unique(canonicalize_alpha(alpha$alpha))),
          expected_alpha_grid
        )) {
      stop("A parent replicate failed the full alpha-grid contract.",
           call. = FALSE)
    }
  }
  all_alpha <- do.call(rbind, lapply(replicates, `[[`, "functional_alpha"))
  if (nrow(all_alpha) != 3200L) {
    stop("The reconstructed parent alpha table does not have 3,200 rows.",
         call. = FALSE)
  }
  saved_all_alpha <- utils::read.csv(
    file.path(
      input_dir, "summary", "all_replicate_functional_alpha_curves.csv"
    ),
    stringsAsFactors = FALSE
  )
  if (!isTRUE(all.equal(
    all_alpha,
    saved_all_alpha,
    tolerance = 1e-14,
    check.attributes = FALSE
  ))) {
    stop("The saved parent alpha table does not reproduce from replicates.",
         call. = FALSE)
  }
  old_middle <- utils::read.csv(
    file.path(input_dir, "summary", "ideal_middle_curve.csv"),
    stringsAsFactors = FALSE
  )
  if (nrow(old_middle) != 60L || min(old_middle$alpha) <= 0.05) {
    stop("The parent cache does not have the expected 60-row boundary defect.",
         call. = FALSE)
  }

  rebuilt <- rebuild_middle_summaries(
    all_alpha = all_alpha,
    mechanisms = mechanisms,
    expected_seeds = expected_seeds,
    lower = 0.05,
    upper = 0.20,
    tolerance = 1e-12
  )
  if (nrow(rebuilt$curve) != 62L ||
      nrow(rebuilt$primary) != 2L ||
      min(rebuilt$curve$alpha) != 0.05 ||
      any(rebuilt$primary$alpha_min != 0.05)) {
    stop("The rebuilt summaries failed the corrected boundary contract.",
         call. = FALSE)
  }
  if (isTRUE(preflight_only)) {
    message(
      "Validated the immutable parent and rebuilt a 62-row Middle summary in memory."
    )
    return(invisible(list(
      curve = rebuilt$curve,
      primary = rebuilt$primary,
      parent_inventory = parent_inventory_before
    )))
  }

  copy_directory_tree(input_dir, temporary_dir)
  atomic_save_rds(
    parent_manifest,
    file.path(temporary_dir, "parent_manifest.rds")
  )
  atomic_save_rds(
    parent_configuration,
    file.path(temporary_dir, "parent_configuration.rds")
  )
  derived_configuration <- parent_configuration
  derived_configuration$schema_version <- derived_schema
  derived_configuration$production_output_id <- parent_result_id
  derived_configuration$output_id <- derived_result_id
  derived_configuration$summary_repair <- list(
    reason = paste(
      "Restore alpha 0.05 omitted by a strict floating-point lower-bound",
      "comparison in the production Middle summary."
    ),
    lower = 0.05,
    upper = 0.20,
    tolerance = 1e-12,
    canonical_digits = 12L,
    parent_middle_rows = 60L,
    corrected_middle_rows = 62L,
    no_refit = TRUE
  )
  atomic_save_rds(
    derived_configuration,
    file.path(temporary_dir, "configuration.rds")
  )
  atomic_write_csv(
    rebuilt$curve,
    file.path(temporary_dir, "summary", "ideal_middle_curve.csv")
  )
  atomic_write_csv(
    rebuilt$primary,
    file.path(temporary_dir, "summary", "primary_middle_summary.csv")
  )

  parent_inventory_after <- relative_file_inventory(input_dir)
  if (!identical(parent_inventory_after, parent_inventory_before)) {
    stop("The immutable parent cache changed during the repair.", call. = FALSE)
  }
  parent_manifest_sha256 <- sha256_file(parent_manifest_path)
  repair_script_sha256 <- sha256_file(script_path)
  copied_parent_replicate_hashes <- relative_file_inventory(temporary_dir)[
    substring(
      replicate_paths,
      nchar(input_dir) + 2L
    )
  ]
  original_parent_replicate_hashes <- parent_inventory_before[
    substring(replicate_paths, nchar(input_dir) + 2L)
  ]
  if (!identical(
    unname(copied_parent_replicate_hashes),
    unname(original_parent_replicate_hashes)
  )) {
    stop("A copied replicate differs from the immutable parent.", call. = FALSE)
  }

  artifact_inventory <- relative_file_inventory(temporary_dir)
  artifact_inventory <- artifact_inventory[
    !names(artifact_inventory) %in% c("manifest.rds", "complete.flag")
  ]
  derived_manifest <- list(
    schema_version = derived_schema,
    result_id = derived_result_id,
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    parent = list(
      result_id = parent_result_id,
      schema_version = parent_schema,
      manifest_sha256 = parent_manifest_sha256,
      tree_sha256 = parent_inventory_before
    ),
    production_source_provenance = parent_manifest$source_provenance,
    package_provenance = parent_manifest$package_provenance,
    repair_provenance = list(
      script = script_path,
      script_sha256 = repair_script_sha256,
      lower = 0.05,
      upper = 0.20,
      tolerance = 1e-12,
      canonical_digits = 12L,
      parent_middle_rows = 60L,
      corrected_middle_rows = 62L,
      no_refit = TRUE
    ),
    configuration = derived_configuration,
    replicate_keys = observed_keys,
    artifact_sha256 = artifact_inventory
  )
  atomic_save_rds(
    derived_manifest,
    file.path(temporary_dir, "manifest.rds")
  )
  atomic_write_lines(
    c(
      paste0("result_id=", derived_result_id),
      paste0("completed_at=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      paste0("schema_version=", derived_schema),
      paste0("parent_result_id=", parent_result_id),
      paste0("parent_manifest_sha256=", parent_manifest_sha256),
      paste0("repair_script_sha256=", repair_script_sha256),
      "summary_alpha_min=0.05",
      "summary_alpha_max=0.20",
      "summary_alpha_tolerance=1e-12",
      "summary_middle_rows=62",
      "refit=FALSE",
      "replicates=10"
    ),
    file.path(temporary_dir, "complete.flag")
  )
  if (!file.rename(temporary_dir, output_dir)) {
    stop("Unable to promote the corrected derived cache: ", output_dir,
         call. = FALSE)
  }
  message("Corrected derived cache completed: ", output_dir)
  invisible(list(
    output_dir = output_dir,
    curve = rebuilt$curve,
    primary = rebuilt$primary,
    manifest = derived_manifest
  ))
}

repair_main <- function() {
  file_argument <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(file_argument) != 1L) {
    stop("Could not resolve the repair script path.", call. = FALSE)
  }
  script_path <- normalizePath(
    sub("^--file=", "", file_argument[[1L]]),
    winslash = "/",
    mustWork = TRUE
  )
  input_dir <- get_arg("--input-dir", "")
  output_dir <- get_arg("--output-dir", "")
  if (!nzchar(input_dir) || !nzchar(output_dir)) {
    stop("--input-dir and --output-dir are required.", call. = FALSE)
  }
  repair_alpha_boundary_summary(
    input_dir = input_dir,
    output_dir = output_dir,
    script_path = script_path,
    preflight_only = as_flag(get_arg("--preflight-only", "false"))
  )
}

if (sys.nframe() == 0L) {
  repair_main()
}
