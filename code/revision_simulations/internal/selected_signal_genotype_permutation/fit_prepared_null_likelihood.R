#!/usr/bin/env Rscript

# Fit prepared matched-null likelihood rows in a clean R session.

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

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required.")
}

input_argument <- get_arg("--input-cache", "")
output_argument <- get_arg("--output-cache", "")
num_cores <- as.integer(get_arg("--num-cores", "8"))
if (!nzchar(input_argument) || !file.exists(input_argument) ||
    !nzchar(output_argument) || file.exists(output_argument) ||
    !dir.exists(dirname(output_argument)) ||
    length(num_cores) != 1L || is.na(num_cores) || num_cores < 1L) {
  stop("Invalid input cache, output cache, or core count.")
}
input_path <- normalizePath(input_argument, mustWork = TRUE)
output_path <- file.path(
  normalizePath(dirname(output_argument), mustWork = TRUE),
  basename(output_argument)
)
source_input_md5 <- unname(tools::md5sum(input_path))
input_bundle <- readRDS(input_path)
required_names <- c(
  "cache_version", "output_id", "permutation_method",
  "target_selection_method", "seed", "selection_seed",
  "selected_pair_keys", "null_unit_keys", "null_datasets",
  "original_settings", "original_psd_grid", "generated_at"
)
if (!identical(input_bundle$cache_version, "selected_signal_null_input_v1") ||
    !all(required_names %in% names(input_bundle)) ||
    length(input_bundle$null_datasets) < 1L ||
    length(input_bundle$null_datasets) != length(input_bundle$null_unit_keys) ||
    !identical(names(input_bundle$null_datasets),
               input_bundle$null_unit_keys) ||
    anyDuplicated(input_bundle$null_unit_keys) ||
    length(input_bundle$original_psd_grid) != 52L) {
  stop("The prepared null-input cache failed structural validation.")
}

settings <- input_bundle$original_settings
required_settings <- c(
  "num_basis", "order", "betaprec", "pred_step", "penalty"
)
if (!all(required_settings %in% names(settings))) {
  stop("The prepared null-input cache has incomplete FASH settings.")
}

message(
  "Computing ",
  length(input_bundle$null_datasets),
  " prepared null likelihood rows in a clean R session."
)
fit_start <- proc.time()[["elapsed"]]
fit_capture <- capture_warnings(fashr::fash(
  Y = "beta",
  smooth_var = "time",
  S = "SE",
  data_list = input_bundle$null_datasets,
  num_basis = settings$num_basis,
  order = settings$order,
  betaprec = settings$betaprec,
  pred_step = settings$pred_step,
  penalty = settings$penalty,
  grid = input_bundle$original_psd_grid,
  num_cores = num_cores,
  verbose = TRUE
))
fit_elapsed <- proc.time()[["elapsed"]] - fit_start
fit <- fit_capture$value
if (nrow(fit$L_matrix) != length(input_bundle$null_unit_keys) ||
    ncol(fit$L_matrix) != length(input_bundle$original_psd_grid) ||
    !isTRUE(all.equal(
      fit$psd_grid,
      input_bundle$original_psd_grid,
      tolerance = 0
    )) || anyNA(fit$L_matrix) || any(is.nan(fit$L_matrix)) ||
    any(fit$L_matrix == Inf)) {
  stop("The clean-session null likelihood fit failed validation.")
}
names(fit$fash_data$data_list) <- input_bundle$null_unit_keys
names(fit$fash_data$S) <- input_bundle$null_unit_keys
rownames(fit$L_matrix) <- input_bundle$null_unit_keys
names(fit$lfdr) <- input_bundle$null_unit_keys

output_bundle <- list(
  cache_version = "selected_signal_null_fit_v1",
  source_input_path = input_path,
  source_input_md5 = source_input_md5,
  output_id = input_bundle$output_id,
  null_unit_keys = input_bundle$null_unit_keys,
  fit = fit,
  elapsed_seconds = unname(fit_elapsed),
  warnings = fit_capture$warnings,
  num_cores = num_cores,
  r_version = R.version.string,
  package_versions = c(fashr = as.character(utils::packageVersion("fashr"))),
  generated_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
temporary_path <- paste0(output_path, ".tmp_", Sys.getpid())
saveRDS(output_bundle, temporary_path)
if (!file.rename(temporary_path, output_path)) {
  unlink(temporary_path)
  stop("Could not atomically save the clean-session null fit cache.")
}

cat("Clean-session null likelihood fit completed.\n")
cat("Output: ", output_path, "\n", sep = "")
cat("Rows: ", nrow(fit$L_matrix), "\n", sep = "")
cat("Elapsed seconds: ", format(fit_elapsed, digits = 8), "\n", sep = "")
cat("Warnings: ", length(fit_capture$warnings), "\n", sep = "")
