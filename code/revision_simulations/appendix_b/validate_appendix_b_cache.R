source("code/revision_simulations/appendix_b/appendix_b_helpers.R")

CACHE_DIR <- "data/appendixB/fashr0143"
LOCAL_HELPER <-
  "code/revision_simulations/appendix_b/appendix_b_helpers.R"
LOCAL_RUNNER <-
  "code/revision_simulations/appendix_b/run_appendix_b_fashr0143.R"

HISTORICAL_SHA256 <- c(
  "data/appendixB/fash_fit_example.RData" =
    "7de6a0569bf20ac015b1ff9230490ffab330a3bffdc7d22ad655ae9775325718",
  "data/appendixB/simulation_result_dense_grid.RData" =
    "325526dbe525ed39fa25aa92a08443549b99c2fb1b52fb5238873ab506057755",
  "data/appendixB/simulation_result_dense_grid_J1000.RData" =
    "7336bf2aeaeab6e945ff737cf46c5b35a0a31b1abb715ca8385ce93ca3bbbb2f",
  "data/appendixB/simulation_result_denser_grid.RData" =
    "8b5ece22f87f88206d044b040f711007e90e9536bccbb968403cb1fea9d3caff",
  "data/appendixB/simulation_result_denser_grid_J1000.RData" =
    "0ec2b0eb7dd339cbdaa6aba6a6980d72d22e00ef8cc5adcdbfef1140eb35b1ef",
  "data/appendixB/simulation_result_loose_grid.RData" =
    "eb298b8823f1e18d6c934808d298f035a86efe8bbd8ae4f139acb8beea66cba8"
)

sha256_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required file does not exist: ", path, call. = FALSE)
  }
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

require_equal <- function(observed, expected, label, tolerance = 0) {
  comparison <- all.equal(
    observed,
    expected,
    tolerance = tolerance,
    check.attributes = TRUE
  )
  if (!isTRUE(comparison)) {
    stop(
      label, " failed: ", paste(comparison, collapse = "; "),
      call. = FALSE
    )
  }
  invisible(observed)
}

required_files <- c(
  "complete.flag",
  "grid_summary.rds",
  "grid_summary.csv",
  "focused_example.rds",
  "focused_alpha_curve.csv",
  "manifest.rds"
)
required_paths <- file.path(CACHE_DIR, required_files)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "The fashr 0.1.43 Appendix B cache is incomplete. Missing: ",
    paste(missing_paths, collapse = ", "),
    call. = FALSE
  )
}

package_provenance <- require_appendix_b_fashr()
manifest <- readRDS(file.path(CACHE_DIR, "manifest.rds"))
grid_summary <- readRDS(file.path(CACHE_DIR, "grid_summary.rds"))
grid_csv <- utils::read.csv(
  file.path(CACHE_DIR, "grid_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
focused_example <- readRDS(file.path(CACHE_DIR, "focused_example.rds"))
focused_alpha_csv <- utils::read.csv(
  file.path(CACHE_DIR, "focused_alpha_curve.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  identical(manifest$schema_version, "appendix-b-fashr0143-v1"),
  identical(manifest$result_id, "appendix_b_fashr0143"),
  identical(
    manifest$package_provenance$version,
    APPENDIX_B_FASHR_VERSION
  ),
  identical(
    manifest$package_provenance$remote_sha,
    APPENDIX_B_FASHR_REMOTE_SHA
  ),
  identical(
    package_provenance$version,
    manifest$package_provenance$version
  ),
  identical(
    package_provenance$remote_sha,
    manifest$package_provenance$remote_sha
  )
)

source_sha <- c(
  helper = sha256_file(LOCAL_HELPER),
  runner = sha256_file(LOCAL_RUNNER)
)
stopifnot(
  identical(unname(source_sha[["helper"]]), manifest$source$helper_sha256),
  identical(unname(source_sha[["runner"]]), manifest$source$runner_sha256)
)

output_paths <- c(
  grid_summary_rds = file.path(CACHE_DIR, "grid_summary.rds"),
  grid_summary_csv = file.path(CACHE_DIR, "grid_summary.csv"),
  focused_example_rds = file.path(CACHE_DIR, "focused_example.rds"),
  focused_alpha_curve_csv = file.path(CACHE_DIR, "focused_alpha_curve.csv")
)
observed_output_sha <- vapply(output_paths, sha256_file, character(1))
expected_output_sha <- unlist(manifest$output_sha256, use.names = TRUE)
require_equal(
  observed_output_sha[names(expected_output_sha)],
  expected_output_sha,
  "Retained output SHA-256 validation"
)

required_grid_columns <- c(
  "setting", "grid_spacing", "rho_dynamic", "rho_nonlinear",
  "true_pi0_iwp1", "true_pi0_iwp2", "J", "n_nondynamic",
  "n_linear", "n_nonlinear", "stream_seed", "raw_pi0_iwp1",
  "bf_pi0_iwp1", "raw_pi0_iwp2", "bf_pi0_iwp2",
  "warning_count", "elapsed_seconds", "warning_text"
)
stopifnot(
  identical(names(grid_summary), required_grid_columns),
  nrow(grid_summary) == 92L,
  !anyDuplicated(grid_summary[c("setting", "rho_dynamic")]),
  identical(sort(unique(grid_summary$setting)), c("denser", "original")),
  all(grid_summary$J == 1000L),
  all(grid_summary$stream_seed == 12345L)
)

expected_boundary_warning <- paste(
  "The estimated prior weight of the null component (PSD = 0) is zero;",
  "lfdr is set to 0 for all datasets."
)
warning_rows <- grid_summary$warning_count > 0L
raw_iwp1_boundary_rows <- grid_summary$raw_pi0_iwp1 == 0
stopifnot(
  all(grid_summary$warning_count %in% c(0L, 1L)),
  identical(warning_rows, raw_iwp1_boundary_rows),
  all(grid_summary$warning_count[warning_rows] == 1L),
  all(grid_summary$warning_text[warning_rows] == expected_boundary_warning),
  all(!nzchar(grid_summary$warning_text[!warning_rows]))
)

expected_rho <- seq(0.05, 0.50, by = 0.01)
for (setting_name in c("original", "denser")) {
  setting_rows <- grid_summary[grid_summary$setting == setting_name, ]
  setting_rows <- setting_rows[order(setting_rows$rho_dynamic), ]
  expected_spacing <- if (setting_name == "original") 0.2 else 0.1
  stopifnot(
    identical(setting_rows$rho_dynamic, expected_rho),
    all(setting_rows$grid_spacing == expected_spacing),
    all(setting_rows$rho_nonlinear == setting_rows$rho_dynamic / 2),
    all(setting_rows$true_pi0_iwp1 == 1 - setting_rows$rho_dynamic),
    all(setting_rows$true_pi0_iwp2 == 1 - setting_rows$rho_nonlinear)
  )
}

for (row_index in seq_len(nrow(grid_summary))) {
  row <- grid_summary[row_index, ]
  expected_counts <- appendix_b_class_counts(
    J = row$J,
    rho_dynamic = row$rho_dynamic,
    rho_nonlinear = row$rho_nonlinear
  )
  observed_counts <- c(
    nondynamic = row$n_nondynamic,
    linear = row$n_linear,
    nonlinear = row$n_nonlinear
  )
  require_equal(
    as.integer(observed_counts),
    as.integer(expected_counts),
    paste0("Class counts at grid row ", row_index)
  )
}

weight_columns <- c(
  "raw_pi0_iwp1", "bf_pi0_iwp1", "raw_pi0_iwp2", "bf_pi0_iwp2"
)
weight_values <- unlist(grid_summary[weight_columns], use.names = FALSE)
stopifnot(
  all(is.finite(weight_values)),
  all(weight_values >= 0),
  all(weight_values <= 1),
  all(is.finite(grid_summary$elapsed_seconds)),
  all(grid_summary$elapsed_seconds > 0)
)

grid_csv_numeric <- names(grid_summary)[vapply(grid_summary, is.numeric, logical(1))]
for (column in grid_csv_numeric) {
  require_equal(
    grid_csv[[column]],
    grid_summary[[column]],
    paste0("grid_summary.csv column ", column),
    tolerance = 1e-12
  )
}
for (column in setdiff(names(grid_summary), grid_csv_numeric)) {
  require_equal(
    grid_csv[[column]],
    grid_summary[[column]],
    paste0("grid_summary.csv column ", column)
  )
}

stopifnot(
  identical(focused_example$schema_version, "appendix-b-focused-v1"),
  identical(focused_example$parameters$J, 1200L),
  identical(focused_example$parameters$rho_dynamic, 0.2),
  identical(focused_example$parameters$rho_nonlinear, 0.1),
  identical(focused_example$parameters$penalty, 10),
  identical(focused_example$parameters$num_basis, 20L),
  identical(focused_example$parameters$stream_seed, 12345L),
  identical(
    focused_example$package_provenance$version,
    APPENDIX_B_FASHR_VERSION
  ),
  identical(
    focused_example$package_provenance$remote_sha,
    APPENDIX_B_FASHR_REMOTE_SHA
  ),
  length(focused_example$fit_bundle$order1$warnings) == 0L,
  length(focused_example$fit_bundle$order2$warnings) == 0L,
  nrow(focused_example$fit_bundle$truth) == 1200L,
  length(focused_example$data_bundle$data_list) == 1200L,
  identical(
    focused_example$data_bundle$truth,
    focused_example$fit_bundle$truth
  )
)

focused_alpha_recomputed <- focused_alpha_curve(focused_example$fit_bundle)
require_equal(
  focused_alpha_csv,
  focused_alpha_recomputed,
  "Focused alpha-curve recomputation",
  tolerance = 1e-12
)

historical_observed_sha <- vapply(
  names(HISTORICAL_SHA256),
  sha256_file,
  character(1)
)
require_equal(
  unname(historical_observed_sha),
  unname(HISTORICAL_SHA256),
  "Historical Appendix B cache preservation"
)

load_result_df <- function(path) {
  environment <- new.env(parent = emptyenv())
  load(path, envir = environment)
  if (!exists("result_df", envir = environment, inherits = FALSE)) {
    stop(path, " does not contain result_df.", call. = FALSE)
  }
  environment$result_df
}

old_grid <- rbind(
  transform(
    load_result_df("data/appendixB/simulation_result_dense_grid.RData"),
    setting = "original"
  ),
  transform(
    load_result_df("data/appendixB/simulation_result_denser_grid.RData"),
    setting = "denser"
  )
)
old_grid$rho_dynamic <- 1 - old_grid$pi_00
old_grid$rho_key <- sprintf("%.2f", old_grid$rho_dynamic)
grid_summary$rho_key <- sprintf("%.2f", grid_summary$rho_dynamic)
grid_comparison <- merge(
  old_grid,
  grid_summary,
  by = c("setting", "rho_key"),
  all = FALSE,
  sort = TRUE
)
stopifnot(nrow(grid_comparison) == 92L)

grid_delta_summary <- data.frame(
  quantity = c(
    "Raw pi0, IWP1", "BF pi0, IWP1", "Raw pi0, IWP2", "BF pi0, IWP2"
  ),
  mean_new_minus_old = c(
    mean(grid_comparison$raw_pi0_iwp1 - grid_comparison$hat_pi_00),
    mean(grid_comparison$bf_pi0_iwp1 - grid_comparison$tilde_pi_00),
    mean(grid_comparison$raw_pi0_iwp2 - grid_comparison$hat_pi_01),
    mean(grid_comparison$bf_pi0_iwp2 - grid_comparison$tilde_pi_01)
  ),
  max_absolute_delta = c(
    max(abs(grid_comparison$raw_pi0_iwp1 - grid_comparison$hat_pi_00)),
    max(abs(grid_comparison$bf_pi0_iwp1 - grid_comparison$tilde_pi_00)),
    max(abs(grid_comparison$raw_pi0_iwp2 - grid_comparison$hat_pi_01)),
    max(abs(grid_comparison$bf_pi0_iwp2 - grid_comparison$tilde_pi_01))
  ),
  stringsAsFactors = FALSE
)

old_environment <- new.env(parent = emptyenv())
load("data/appendixB/fash_fit_example.RData", envir = old_environment)
old_focused <- list(
  order1 = list(
    raw = old_environment$fash_fit1,
    bf = old_environment$fash_fit1_update
  ),
  order2 = list(
    raw = old_environment$fash_fit2,
    bf = old_environment$fash_fit2_update
  ),
  truth = focused_example$fit_bundle$truth
)
old_alpha <- focused_alpha_curve(old_focused)
alpha_005_old <- old_alpha[abs(old_alpha$alpha - 0.05) < 1e-12, ]
alpha_005_new <- focused_alpha_recomputed[
  abs(focused_alpha_recomputed$alpha - 0.05) < 1e-12,
]
focused_comparison <- merge(
  alpha_005_old,
  alpha_005_new,
  by = c("order", "stage", "alpha"),
  suffixes = c("_old", "_new"),
  sort = TRUE
)

cat("Appendix B fashr 0.1.43 cache validation passed.\n\n")
cat("Grid comparison with the historical retained summaries:\n")
print(grid_delta_summary, row.names = FALSE)
cat("\nFocused example at cumulative-lfdr alpha = 0.05:\n")
print(
  focused_comparison[c(
    "order", "stage", "discoveries_old", "discoveries_new",
    "realized_fdp_old", "realized_fdp_new", "power_old", "power_new"
  )],
  row.names = FALSE
)
