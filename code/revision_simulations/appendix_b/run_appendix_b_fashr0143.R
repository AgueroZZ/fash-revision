workspace_dir <- Sys.getenv(
  "FASH_WORKSPACE",
  unset = "/project/mstephens/ziangzhang/fash/workspace"
)
source(file.path(workspace_dir, "config", "paths.R"))

helper_path <- file.path(WORKSPACE_DIR, "code", "appendix_b_helpers.R")
require_input_file(helper_path)
source(helper_path)

RESULT_ID <- "appendix_b_fashr0143"
RESULT_PARENT <- file.path(
  WORKSPACE_RESULTS_DIR,
  "revision_simulations"
)
FINAL_DIR <- file.path(RESULT_PARENT, RESULT_ID)
PARTIAL_DIR <- file.path(RESULT_PARENT, paste0(RESULT_ID, "_partial"))

GRID_J <- 1000L
FOCUSED_J <- 1200L
STREAM_SEED <- 12345L
SIGMA_VEC <- c(0.1, 0.3, 0.5)
RHO_VALUES <- seq(0.05, 0.50, by = 0.01)
GRID_PENALTY <- 1
FOCUSED_PENALTY <- 10
NUM_BASIS <- 20L
NUM_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "1"))

GRID_SETTINGS <- list(
  original = list(spacing = 0.2),
  denser = list(spacing = 0.1)
)

log_message <- function(...) {
  cat(
    sprintf("[%s] ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    ...,
    "\n",
    sep = ""
  )
  flush.console()
}

sha256_file <- function(path) {
  require_input_file(path)
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

atomic_save_rds <- function(object, path) {
  temporary_path <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)
  saveRDS(object, temporary_path, version = 3)
  if (!file.rename(temporary_path, path)) {
    stop("Unable to atomically replace ", path, ".", call. = FALSE)
  }
  invisible(path)
}

atomic_write_csv <- function(object, path) {
  temporary_path <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)
  utils::write.csv(object, temporary_path, row.names = FALSE, quote = TRUE)
  if (!file.rename(temporary_path, path)) {
    stop("Unable to atomically replace ", path, ".", call. = FALSE)
  }
  invisible(path)
}

validate_checkpoint <- function(checkpoint, setting_name, spacing) {
  required_names <- c(
    "schema_version", "setting", "spacing", "stream_seed", "rows",
    "rng_state", "package_provenance"
  )
  if (!is.list(checkpoint) ||
      !all(required_names %in% names(checkpoint))) {
    stop("Invalid checkpoint schema for ", setting_name, ".", call. = FALSE)
  }
  stopifnot(
    identical(checkpoint$schema_version, "appendix-b-grid-checkpoint-v1"),
    identical(checkpoint$setting, setting_name),
    identical(checkpoint$spacing, spacing),
    identical(checkpoint$stream_seed, STREAM_SEED),
    identical(
      checkpoint$package_provenance$version,
      APPENDIX_B_FASHR_VERSION
    ),
    identical(
      checkpoint$package_provenance$remote_sha,
      APPENDIX_B_FASHR_REMOTE_SHA
    )
  )
  if (length(checkpoint$rows) > length(RHO_VALUES)) {
    stop("Checkpoint has too many completed rows for ", setting_name, ".",
         call. = FALSE)
  }
  invisible(checkpoint)
}

run_grid_setting <- function(setting_name, spacing, package_provenance) {
  checkpoint_path <- file.path(
    PARTIAL_DIR,
    paste0("checkpoint_grid_", setting_name, ".rds")
  )
  grid <- sort(c(0, exp(-0.5 * seq(0, 10, by = spacing))))

  if (file.exists(checkpoint_path)) {
    checkpoint <- readRDS(checkpoint_path)
    validate_checkpoint(checkpoint, setting_name, spacing)
    assign(".Random.seed", checkpoint$rng_state, envir = .GlobalEnv)
    rows <- checkpoint$rows
    log_message(
      "Resuming ", setting_name, " grid after ", length(rows),
      " completed rho values."
    )
  } else {
    set.seed(STREAM_SEED)
    rows <- list()
    checkpoint <- list(
      schema_version = "appendix-b-grid-checkpoint-v1",
      setting = setting_name,
      spacing = spacing,
      stream_seed = STREAM_SEED,
      rows = rows,
      rng_state = .Random.seed,
      package_provenance = package_provenance
    )
    atomic_save_rds(checkpoint, checkpoint_path)
  }

  next_index <- length(rows) + 1L
  if (next_index <= length(RHO_VALUES)) {
    for (index in seq.int(next_index, length(RHO_VALUES))) {
      rho_dynamic <- RHO_VALUES[[index]]
      rho_nonlinear <- rho_dynamic / 2
      started <- proc.time()[["elapsed"]]
      log_message(
        "Starting setting=", setting_name,
        " spacing=", spacing,
        " rho=", sprintf("%.2f", rho_dynamic),
        " (", index, "/", length(RHO_VALUES), ")."
      )

      data_bundle <- build_appendix_b_datasets(
        J = GRID_J,
        rho_dynamic = rho_dynamic,
        rho_nonlinear = rho_nonlinear,
        sigma_vec = SIGMA_VEC
      )
      fit_bundle <- fit_appendix_b_models(
        data_bundle = data_bundle,
        grid = grid,
        penalty = GRID_PENALTY,
        num_basis = NUM_BASIS,
        num_cores = NUM_CORES,
        verbose = FALSE
      )
      elapsed_seconds <- proc.time()[["elapsed"]] - started
      row <- summarize_appendix_b_grid_fit(
        fit_bundle = fit_bundle,
        setting = setting_name,
        grid_spacing = spacing,
        rho_dynamic = rho_dynamic,
        rho_nonlinear = rho_nonlinear,
        J = GRID_J,
        stream_seed = STREAM_SEED,
        elapsed_seconds = elapsed_seconds
      )
      row$warning_text <- paste(
        unique(c(
          fit_bundle$order1$warnings,
          fit_bundle$order2$warnings
        )),
        collapse = " | "
      )
      rows[[index]] <- row

      checkpoint$rows <- rows
      checkpoint$rng_state <- .Random.seed
      checkpoint$last_completed_at <- format(
        Sys.time(),
        tz = "UTC",
        usetz = TRUE
      )
      atomic_save_rds(checkpoint, checkpoint_path)
      log_message(
        "Completed setting=", setting_name,
        " rho=", sprintf("%.2f", rho_dynamic),
        " in ", sprintf("%.1f", elapsed_seconds), " seconds."
      )

      rm(data_bundle, fit_bundle)
      invisible(gc())
    }
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  stopifnot(
    nrow(result) == length(RHO_VALUES),
    identical(result$rho_dynamic, RHO_VALUES),
    all(result$J == GRID_J),
    all(result$n_nondynamic + result$n_linear + result$n_nonlinear == GRID_J)
  )
  result
}

run_focused_example <- function(package_provenance) {
  checkpoint_path <- file.path(PARTIAL_DIR, "checkpoint_focused_example.rds")
  if (file.exists(checkpoint_path)) {
    focused <- readRDS(checkpoint_path)
    required_names <- c(
      "schema_version", "data_bundle", "fit_bundle", "parameters",
      "package_provenance"
    )
    if (!is.list(focused) || !all(required_names %in% names(focused))) {
      stop("Invalid focused-example checkpoint.", call. = FALSE)
    }
    stopifnot(
      identical(focused$schema_version, "appendix-b-focused-v1"),
      identical(
        focused$package_provenance$version,
        APPENDIX_B_FASHR_VERSION
      ),
      identical(
        focused$package_provenance$remote_sha,
        APPENDIX_B_FASHR_REMOTE_SHA
      )
    )
    log_message("Loaded the completed focused-example checkpoint.")
    return(focused)
  }

  set.seed(STREAM_SEED)
  grid <- sort(c(0, exp(-0.5 * seq(0, 10, by = 0.2))))
  data_bundle <- build_appendix_b_datasets(
    J = FOCUSED_J,
    rho_dynamic = 0.2,
    rho_nonlinear = 0.1,
    sigma_vec = SIGMA_VEC
  )
  started <- proc.time()[["elapsed"]]
  fit_bundle <- fit_appendix_b_models(
    data_bundle = data_bundle,
    grid = grid,
    penalty = FOCUSED_PENALTY,
    num_basis = NUM_BASIS,
    num_cores = NUM_CORES,
    verbose = TRUE
  )
  elapsed_seconds <- proc.time()[["elapsed"]] - started

  focused <- list(
    schema_version = "appendix-b-focused-v1",
    data_bundle = data_bundle,
    fit_bundle = fit_bundle,
    parameters = list(
      J = FOCUSED_J,
      rho_dynamic = 0.2,
      rho_nonlinear = 0.1,
      sigma_vec = SIGMA_VEC,
      grid = grid,
      grid_spacing = 0.2,
      penalty = FOCUSED_PENALTY,
      num_basis = NUM_BASIS,
      stream_seed = STREAM_SEED,
      elapsed_seconds = elapsed_seconds
    ),
    package_provenance = package_provenance,
    generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  atomic_save_rds(focused, checkpoint_path)
  log_message(
    "Completed focused example in ",
    sprintf("%.1f", elapsed_seconds), " seconds."
  )
  focused
}

if (dir.exists(FINAL_DIR) || file.exists(FINAL_DIR)) {
  stop("Refusing to overwrite completed output: ", FINAL_DIR, call. = FALSE)
}
ensure_output_directory(RESULT_PARENT)
ensure_output_directory(PARTIAL_DIR)

package_provenance <- require_appendix_b_fashr()
runner_path <- normalizePath(
  file.path(WORKSPACE_DIR, "scripts", "10_appendix_b_simulation.R"),
  winslash = "/",
  mustWork = TRUE
)
helper_path <- normalizePath(helper_path, winslash = "/", mustWork = TRUE)

log_message("Appendix B output: ", FINAL_DIR)
log_message(
  "fashr ", package_provenance$version,
  " at ", package_provenance$remote_sha
)
log_message("Using ", NUM_CORES, " Slurm CPUs.")

grid_rows <- lapply(names(GRID_SETTINGS), function(setting_name) {
  run_grid_setting(
    setting_name = setting_name,
    spacing = GRID_SETTINGS[[setting_name]]$spacing,
    package_provenance = package_provenance
  )
})
grid_summary <- do.call(rbind, grid_rows)
rownames(grid_summary) <- NULL
stopifnot(
  nrow(grid_summary) == 2L * length(RHO_VALUES),
  !anyDuplicated(grid_summary[c("setting", "rho_dynamic")]),
  all(is.finite(unlist(grid_summary[c(
    "raw_pi0_iwp1", "bf_pi0_iwp1", "raw_pi0_iwp2", "bf_pi0_iwp2"
  )]))),
  all(unlist(grid_summary[c(
    "raw_pi0_iwp1", "bf_pi0_iwp1", "raw_pi0_iwp2", "bf_pi0_iwp2"
  )]) >= 0),
  all(unlist(grid_summary[c(
    "raw_pi0_iwp1", "bf_pi0_iwp1", "raw_pi0_iwp2", "bf_pi0_iwp2"
  )]) <= 1)
)

focused_example <- run_focused_example(package_provenance)
focused_alpha <- focused_alpha_curve(focused_example$fit_bundle)

grid_rds_path <- file.path(PARTIAL_DIR, "grid_summary.rds")
grid_csv_path <- file.path(PARTIAL_DIR, "grid_summary.csv")
focused_rds_path <- file.path(PARTIAL_DIR, "focused_example.rds")
focused_csv_path <- file.path(PARTIAL_DIR, "focused_alpha_curve.csv")
manifest_path <- file.path(PARTIAL_DIR, "manifest.rds")

atomic_save_rds(grid_summary, grid_rds_path)
atomic_write_csv(grid_summary, grid_csv_path)
atomic_save_rds(focused_example, focused_rds_path)
atomic_write_csv(focused_alpha, focused_csv_path)

manifest <- list(
  schema_version = "appendix-b-fashr0143-v1",
  result_id = RESULT_ID,
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  package_provenance = package_provenance,
  source = list(
    helper_path = helper_path,
    helper_sha256 = sha256_file(helper_path),
    runner_path = runner_path,
    runner_sha256 = sha256_file(runner_path)
  ),
  grid_design = list(
    J = GRID_J,
    stream_seed = STREAM_SEED,
    sigma_vec = SIGMA_VEC,
    rho_values = RHO_VALUES,
    settings = GRID_SETTINGS,
    penalty = GRID_PENALTY,
    num_basis = NUM_BASIS,
    num_cores = NUM_CORES
  ),
  focused_design = focused_example$parameters,
  output_sha256 = list(
    grid_summary_rds = sha256_file(grid_rds_path),
    grid_summary_csv = sha256_file(grid_csv_path),
    focused_example_rds = sha256_file(focused_rds_path),
    focused_alpha_curve_csv = sha256_file(focused_csv_path)
  )
)
atomic_save_rds(manifest, manifest_path)

writeLines(
  c(
    paste0("result_id=", RESULT_ID),
    paste0("completed_at=", manifest$generated_at),
    paste0("fashr_version=", package_provenance$version),
    paste0("fashr_remote_sha=", package_provenance$remote_sha),
    paste0("grid_rows=", nrow(grid_summary)),
    paste0("focused_rows=", nrow(focused_alpha))
  ),
  con = file.path(PARTIAL_DIR, "complete.flag")
)

if (!file.rename(PARTIAL_DIR, FINAL_DIR)) {
  stop(
    "All outputs were written, but the partial directory could not be promoted to ",
    FINAL_DIR, ".",
    call. = FALSE
  )
}

log_message("Appendix B fashr 0.1.43 run completed: ", FINAL_DIR)
