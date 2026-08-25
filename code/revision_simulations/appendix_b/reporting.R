APPENDIX_B_CACHE_DIR <- file.path("data", "appendixB", "fashr0143")
APPENDIX_B_LOCAL_HELPER <- file.path(
  "code", "revision_simulations", "appendix_b", "appendix_b_helpers.R"
)
APPENDIX_B_LOCAL_RUNNER <- file.path(
  "code", "revision_simulations", "appendix_b", "run_appendix_b_fashr0143.R"
)
APPENDIX_B_BOUNDARY_WARNING <- paste(
  "The estimated prior weight of the null component (PSD = 0) is zero;",
  "lfdr is set to 0 for all datasets."
)

source(APPENDIX_B_LOCAL_HELPER)

appendix_b_sha256_file <- function(path) {
  if (!file.exists(path) || dir.exists(path)) {
    stop("Required file does not exist: ", path, call. = FALSE)
  }
  command <- if (identical(Sys.info()[["sysname"]], "Darwin")) {
    "shasum"
  } else {
    "sha256sum"
  }
  arguments <- if (identical(command, "shasum")) {
    c("-a", "256", path)
  } else {
    path
  }
  output <- system2(command, arguments, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Unable to compute SHA-256 for ", path, ".", call. = FALSE)
  }
  sub("[[:space:]].*$", "", output[[1L]])
}

require_reporting_equal <- function(
    observed,
    expected,
    label,
    tolerance = 0) {
  comparison <- all.equal(
    observed,
    expected,
    tolerance = tolerance,
    check.attributes = TRUE
  )
  if (!isTRUE(comparison)) {
    stop(label, " failed: ", paste(comparison, collapse = "; "), call. = FALSE)
  }
  invisible(observed)
}

parse_appendix_b_complete_flag <- function(path) {
  lines <- readLines(path, warn = FALSE)
  pieces <- strsplit(lines, "=", fixed = TRUE)
  if (any(lengths(pieces) != 2L)) {
    stop("The Appendix B complete.flag has an invalid format.", call. = FALSE)
  }
  values <- vapply(pieces, `[[`, character(1), 2L)
  names(values) <- vapply(pieces, `[[`, character(1), 1L)
  values
}

load_appendix_b_reporting_cache <- function(
    cache_dir = APPENDIX_B_CACHE_DIR) {
  required_files <- c(
    "complete.flag",
    "grid_summary.rds",
    "grid_summary.csv",
    "focused_example.rds",
    "focused_alpha_curve.csv",
    "manifest.rds"
  )
  required_paths <- file.path(cache_dir, required_files)
  missing_paths <- required_paths[!file.exists(required_paths)]
  if (length(missing_paths) > 0L) {
    stop(
      "The Appendix B fashr 0.1.43 cache is incomplete. Missing: ",
      paste(missing_paths, collapse = ", "),
      call. = FALSE
    )
  }

  local_package <- require_appendix_b_fashr()
  completion <- parse_appendix_b_complete_flag(
    file.path(cache_dir, "complete.flag")
  )
  manifest <- readRDS(file.path(cache_dir, "manifest.rds"))
  grid_summary <- readRDS(file.path(cache_dir, "grid_summary.rds"))
  grid_csv <- utils::read.csv(
    file.path(cache_dir, "grid_summary.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  focused_example <- readRDS(file.path(cache_dir, "focused_example.rds"))
  focused_alpha <- utils::read.csv(
    file.path(cache_dir, "focused_alpha_curve.csv"),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  stopifnot(
    identical(manifest$schema_version, "appendix-b-fashr0143-v1"),
    identical(manifest$result_id, "appendix_b_fashr0143"),
    identical(completion[["result_id"]], manifest$result_id),
    identical(completion[["fashr_version"]], APPENDIX_B_FASHR_VERSION),
    identical(completion[["fashr_remote_sha"]], APPENDIX_B_FASHR_REMOTE_SHA),
    identical(manifest$package_provenance$version, APPENDIX_B_FASHR_VERSION),
    identical(
      manifest$package_provenance$remote_sha,
      APPENDIX_B_FASHR_REMOTE_SHA
    ),
    identical(local_package$version, manifest$package_provenance$version),
    identical(local_package$remote_sha, manifest$package_provenance$remote_sha),
    identical(
      manifest$source$helper_sha256,
      appendix_b_sha256_file(APPENDIX_B_LOCAL_HELPER)
    ),
    identical(
      manifest$source$runner_sha256,
      appendix_b_sha256_file(APPENDIX_B_LOCAL_RUNNER)
    )
  )

  output_paths <- c(
    grid_summary_rds = file.path(cache_dir, "grid_summary.rds"),
    grid_summary_csv = file.path(cache_dir, "grid_summary.csv"),
    focused_example_rds = file.path(cache_dir, "focused_example.rds"),
    focused_alpha_curve_csv = file.path(
      cache_dir,
      "focused_alpha_curve.csv"
    )
  )
  observed_hashes <- vapply(
    output_paths,
    appendix_b_sha256_file,
    character(1)
  )
  expected_hashes <- unlist(manifest$output_sha256, use.names = TRUE)
  require_reporting_equal(
    observed_hashes[names(expected_hashes)],
    expected_hashes,
    "Appendix B retained output hashes"
  )

  stopifnot(
    nrow(grid_summary) == 92L,
    identical(sort(unique(grid_summary$setting)), c("denser", "original")),
    !anyDuplicated(grid_summary[c("setting", "rho_dynamic")]),
    all(grid_summary$J == 1000L),
    all(grid_summary$stream_seed == 12345L),
    nrow(focused_alpha) == 80L,
    !anyDuplicated(focused_alpha[c("order", "stage", "alpha")]),
    identical(focused_example$schema_version, "appendix-b-focused-v1"),
    identical(focused_example$parameters$J, 1200L),
    identical(focused_example$parameters$rho_dynamic, 0.2),
    identical(focused_example$parameters$rho_nonlinear, 0.1),
    nrow(focused_example$fit_bundle$truth) == 1200L,
    length(focused_example$data_bundle$data_list) == 1200L,
    identical(
      focused_example$data_bundle$truth,
      focused_example$fit_bundle$truth
    )
  )

  warning_rows <- grid_summary$warning_count > 0L
  raw_iwp1_boundary_rows <- grid_summary$raw_pi0_iwp1 == 0
  stopifnot(
    all(grid_summary$warning_count %in% c(0L, 1L)),
    identical(warning_rows, raw_iwp1_boundary_rows),
    all(grid_summary$warning_text[warning_rows] == APPENDIX_B_BOUNDARY_WARNING),
    all(!nzchar(grid_summary$warning_text[!warning_rows])),
    length(focused_example$fit_bundle$order1$warnings) == 0L,
    length(focused_example$fit_bundle$order2$warnings) == 0L
  )

  require_reporting_equal(
    grid_csv,
    grid_summary,
    "Appendix B grid CSV/RDS agreement",
    tolerance = 1e-12
  )
  focused_alpha_recomputed <- focused_alpha_curve(
    focused_example$fit_bundle
  )
  require_reporting_equal(
    focused_alpha,
    focused_alpha_recomputed,
    "Appendix B focused alpha-curve recomputation",
    tolerance = 1e-12
  )

  list(
    completion = completion,
    manifest = manifest,
    grid_summary = grid_summary,
    focused_example = focused_example,
    focused_alpha = focused_alpha,
    boundary_warning_rows = grid_summary[warning_rows, , drop = FALSE]
  )
}

prepare_appendix_b_grid_long <- function(grid_summary) {
  required <- c(
    "setting", "grid_spacing", "rho_dynamic", "rho_nonlinear",
    "true_pi0_iwp1", "true_pi0_iwp2", "raw_pi0_iwp1",
    "bf_pi0_iwp1", "raw_pi0_iwp2", "bf_pi0_iwp2", "warning_count"
  )
  if (!is.data.frame(grid_summary) || !all(required %in% names(grid_summary))) {
    stop("grid_summary does not satisfy the Appendix B reporting contract.",
         call. = FALSE)
  }

  specifications <- list(
    list(order = "IWP1", stage = "Raw EB", truth = "true_pi0_iwp1",
         estimate = "raw_pi0_iwp1"),
    list(order = "IWP1", stage = "BF updated", truth = "true_pi0_iwp1",
         estimate = "bf_pi0_iwp1"),
    list(order = "IWP2", stage = "Raw EB", truth = "true_pi0_iwp2",
         estimate = "raw_pi0_iwp2"),
    list(order = "IWP2", stage = "BF updated", truth = "true_pi0_iwp2",
         estimate = "bf_pi0_iwp2")
  )
  rows <- lapply(specifications, function(specification) {
    data.frame(
      setting = grid_summary$setting,
      grid_spacing = grid_summary$grid_spacing,
      rho_dynamic = grid_summary$rho_dynamic,
      rho_nonlinear = grid_summary$rho_nonlinear,
      order = specification$order,
      stage = specification$stage,
      true_pi0 = grid_summary[[specification$truth]],
      estimated_pi0 = grid_summary[[specification$estimate]],
      boundary_warning = if (
        specification$order == "IWP1" && specification$stage == "Raw EB"
      ) {
        grid_summary$warning_count > 0L
      } else {
        rep(FALSE, nrow(grid_summary))
      },
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result$setting <- factor(
    result$setting,
    levels = c("original", "denser"),
    labels = c("Original grid (spacing 0.2)", "Denser grid (spacing 0.1)")
  )
  result$order <- factor(result$order, levels = c("IWP1", "IWP2"))
  result$stage <- factor(
    result$stage,
    levels = c("Raw EB", "BF updated")
  )
  result
}

prepare_appendix_b_alpha005 <- function(focused_alpha) {
  result <- focused_alpha[abs(focused_alpha$alpha - 0.05) < 1e-12, ]
  if (nrow(result) != 4L ||
      !setequal(result$order, c("IWP1", "IWP2")) ||
      !setequal(result$stage, c("Raw EB", "BF updated"))) {
    stop("The focused alpha-0.05 summary is incomplete.", call. = FALSE)
  }
  result$target <- ifelse(
    result$order == "IWP1",
    "Dynamic (B or C)",
    "Nonlinear dynamic (C)"
  )
  result$order <- factor(result$order, levels = c("IWP1", "IWP2"))
  result$stage <- factor(
    result$stage,
    levels = c("Raw EB", "BF updated")
  )
  result[order(result$order, result$stage), ]
}

prepare_appendix_b_alpha_long <- function(focused_alpha) {
  specifications <- list(
    list(metric = "Realized FDP", column = "realized_fdp"),
    list(metric = "Power", column = "power")
  )
  rows <- lapply(specifications, function(specification) {
    data.frame(
      order = focused_alpha$order,
      stage = focused_alpha$stage,
      alpha = focused_alpha$alpha,
      metric = specification$metric,
      value = focused_alpha[[specification$column]],
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result$order <- factor(result$order, levels = c("IWP1", "IWP2"))
  result$stage <- factor(
    result$stage,
    levels = c("Raw EB", "BF updated")
  )
  result$metric <- factor(
    result$metric,
    levels = c("Realized FDP", "Power")
  )
  result
}

prepare_appendix_b_prior_table <- function(focused_example) {
  specifications <- list(
    list(order = "IWP1", target = "Dynamic (B or C)", stage = "Raw EB",
         fit = focused_example$fit_bundle$order1$raw,
         true_pi0 = 1 - focused_example$parameters$rho_dynamic),
    list(order = "IWP1", target = "Dynamic (B or C)", stage = "BF updated",
         fit = focused_example$fit_bundle$order1$bf,
         true_pi0 = 1 - focused_example$parameters$rho_dynamic),
    list(order = "IWP2", target = "Nonlinear dynamic (C)", stage = "Raw EB",
         fit = focused_example$fit_bundle$order2$raw,
         true_pi0 = 1 - focused_example$parameters$rho_nonlinear),
    list(order = "IWP2", target = "Nonlinear dynamic (C)", stage = "BF updated",
         fit = focused_example$fit_bundle$order2$bf,
         true_pi0 = 1 - focused_example$parameters$rho_nonlinear)
  )
  result <- do.call(rbind, lapply(specifications, function(specification) {
    estimated_pi0 <- extract_fash_null_weight(specification$fit)
    data.frame(
      order = specification$order,
      target = specification$target,
      stage = specification$stage,
      true_pi0 = specification$true_pi0,
      estimated_pi0 = estimated_pi0,
      estimation_error = estimated_pi0 - specification$true_pi0,
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  result$order <- factor(result$order, levels = c("IWP1", "IWP2"))
  result$stage <- factor(
    result$stage,
    levels = c("Raw EB", "BF updated")
  )
  result
}

prepare_appendix_b_cumulative_lfdr <- function(
    focused_example,
    order = c("IWP1", "IWP2"),
    stage = c("Raw EB", "BF updated"),
    alpha = 0.05) {
  order <- match.arg(order)
  stage <- match.arg(stage)
  fit <- if (order == "IWP1") {
    if (stage == "Raw EB") {
      focused_example$fit_bundle$order1$raw
    } else {
      focused_example$fit_bundle$order1$bf
    }
  } else {
    if (stage == "Raw EB") {
      focused_example$fit_bundle$order2$raw
    } else {
      focused_example$fit_bundle$order2$bf
    }
  }
  truth <- focused_example$fit_bundle$truth
  lfdr <- extract_fash_lfdr(fit)
  ranked_index <- order(lfdr, seq_along(lfdr))
  ranked_lfdr <- lfdr[ranked_index]
  result <- data.frame(
    order = order,
    stage = stage,
    rank = seq_along(ranked_index),
    unit_id = truth$unit_id[ranked_index],
    class = truth$class[ranked_index],
    lfdr = ranked_lfdr,
    cumulative_lfdr = cumsum(ranked_lfdr) / seq_along(ranked_lfdr),
    stringsAsFactors = FALSE
  )
  result$called_at_alpha <- result$cumulative_lfdr <= alpha
  result$class <- factor(
    result$class,
    levels = c("nondynamic", "linear", "nonlinear"),
    labels = c("A: Constant", "B: Linear dynamic", "C: Nonlinear dynamic")
  )
  result
}

prepare_appendix_b_trajectory_examples <- function(
    focused_example,
    examples_per_class = 3L) {
  examples_per_class <- as.integer(examples_per_class)
  if (length(examples_per_class) != 1L || is.na(examples_per_class) ||
      examples_per_class < 1L) {
    stop("examples_per_class must be one positive integer.", call. = FALSE)
  }
  truth <- focused_example$data_bundle$truth
  data_list <- focused_example$data_bundle$data_list
  classes <- c("nondynamic", "linear", "nonlinear")
  selected <- unlist(lapply(classes, function(class_name) {
    head(which(truth$class == class_name), examples_per_class)
  }), use.names = FALSE)
  rows <- lapply(selected, function(index) {
    dataset <- data_list[[index]]
    data.frame(
      unit_id = truth$unit_id[[index]],
      class = truth$class[[index]],
      time = dataset$x,
      estimate = dataset$y,
      standard_error = dataset$sd,
      truth = dataset$truef,
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result$class <- factor(
    result$class,
    levels = classes,
    labels = c("A: Constant", "B: Linear dynamic", "C: Nonlinear dynamic")
  )
  result$panel <- factor(
    paste(result$class, result$unit_id, sep = " — "),
    levels = unique(paste(result$class, result$unit_id, sep = " — "))
  )
  result
}
