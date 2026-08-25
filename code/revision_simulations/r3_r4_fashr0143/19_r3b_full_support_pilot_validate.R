#!/usr/bin/env Rscript

# Validate the single-seed R3B full-support pilot. This pilot is a design check,
# not the formal Monte Carlo calibration result.

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[[1L]] == length(args)) return(default)
  args[[hit[[1L]] + 1L]]
}

expect_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

atomic_write_lines <- function(text, path) {
  temporary_path <- paste0(path, ".tmp-", Sys.getpid())
  on.exit(unlink(temporary_path), add = TRUE)
  writeLines(text, temporary_path, useBytes = TRUE)
  if (!file.rename(temporary_path, path)) {
    stop("Unable to atomically write ", path, ".", call. = FALSE)
  }
}

result_dir_argument <- get_arg("--result-dir", "")
expect_true(nzchar(result_dir_argument), "--result-dir is required.")
result_dir <- normalizePath(
  result_dir_argument, winslash = "/", mustWork = TRUE
)
result_id <- basename(result_dir)
expected_result_id <- paste0(
  "r3b_real_genotype_one_per_gene_J6362_open_middle_3_12_",
  "full_support_iwp1_geometry_mixture_fashr0143_seed12345_pilot"
)
expect_true(identical(result_id, expected_result_id), "Unexpected pilot ID.")

configuration <- readRDS(file.path(result_dir, "configuration.rds"))
expect_true(
  identical(configuration$output_id, expected_result_id) &&
    identical(configuration$J, 6362L) &&
    identical(configuration$seed_list, 12345L) &&
    identical(configuration$truth_mechanisms, "raised_cosine") &&
    identical(configuration$middle_window, c(3, 12)) &&
    identical(configuration$middle_boundary, "open") &&
    identical(configuration$middle_expression, "3 < t < 12") &&
    identical(
      configuration$temporal_category_probs,
      c(early = 0.29, middle = 0.42, late = 0.29)
    ) &&
    identical(
      configuration$raised_cosine$center_ranges,
      list(
        early = c(1.5, 1.5),
        middle = c(4.5, 10.5),
        late = c(13.5, 13.5)
      )
    ) &&
    identical(configuration$package_provenance$version, "0.1.43") &&
    identical(
      configuration$package_provenance$remote_sha,
      "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
    ),
  "The pilot configuration failed its frozen contract."
)

replicate_path <- file.path(
  result_dir, "replicates", "raised_cosine_seed_12345.rds"
)
expect_true(file.exists(replicate_path), "The pilot replicate is missing.")
replicate <- readRDS(replicate_path)
expect_true(
  identical(replicate$seed, 12345L) &&
    identical(replicate$truth_mechanism, "raised_cosine") &&
    nrow(replicate$functional_alpha) == 320L &&
    nrow(replicate$functional_alpha_005) == 8L &&
    all(is.finite(replicate$functional_alpha$empirical_fsr)),
  "The pilot replicate failed its result contract."
)
expected_counts <- c(
  "early / switch" = 185L,
  "early / non-switch" = 184L,
  "middle / switch" = 267L,
  "middle / non-switch" = 267L,
  "late / switch" = 185L,
  "late / non-switch" = 184L
)
observed_counts <- stats::setNames(
  as.integer(replicate$truth_group_counts$n_dynamic),
  as.character(replicate$truth_group_counts$truth_group)
)[names(expected_counts)]
expect_true(
  identical(observed_counts, expected_counts),
  "The pilot truth-group counts are invalid."
)

middle <- replicate$functional_alpha[
  replicate$functional_alpha$method == "FASH-IWP1-BF" &
    replicate$functional_alpha$target == "middle",
  , drop = FALSE
]
excess <- middle$empirical_fsr - middle$alpha
maximum_index <- which.max(excess)
validation <- data.frame(
  seed = 12345L,
  truth_mechanism = "raised_cosine",
  method = "FASH-IWP1-BF",
  target = "middle",
  maximum_excess = excess[[maximum_index]],
  alpha_at_maximum_excess = middle$alpha[[maximum_index]],
  empirical_fsr_at_maximum = middle$empirical_fsr[[maximum_index]],
  note = paste(
    "single-seed design pilot only; formal acceptance requires five seeds"
  ),
  stringsAsFactors = FALSE
)
utils::write.csv(
  validation,
  file.path(result_dir, "pilot_scientific_validation.csv"),
  row.names = FALSE
)

selected <- middle[vapply(
  middle$alpha,
  function(alpha) any(abs(alpha - c(0.05, 0.10, 0.15, 0.20)) < 1e-12),
  logical(1)
), c(
  "alpha", "n_discoveries", "false_discoveries", "estimated_fsr",
  "empirical_fsr"
), drop = FALSE]
utils::write.csv(
  selected,
  file.path(result_dir, "pilot_selected_middle_rows.csv"),
  row.names = FALSE
)

atomic_write_lines(
  c(
    "status=complete",
    paste0("result_id=", expected_result_id),
    "scope=single-seed design pilot",
    "formal_calibration_claim=false",
    paste0("completed_at=", format(Sys.time(), tz = "UTC", usetz = TRUE))
  ),
  file.path(result_dir, "pilot_complete.flag")
)

cat("Validated single-seed R3B full-support pilot:\n")
print(selected, row.names = FALSE)
cat("Maximum empirical FSR excess in the pilot:\n")
print(validation, row.names = FALSE)
