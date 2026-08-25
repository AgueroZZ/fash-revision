#!/usr/bin/env Rscript

# Read-only validation for the balanced support-contained R3 cache and the new
# canonical-IWP1 prior-geometry R3 cache.

options(stringsAsFactors = FALSE)

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[[1L]] == length(args)) return(default)
  args[[hit[[1L]] + 1L]]
}

as_flag <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

expect_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

require_file <- function(path, label) {
  expect_true(file.exists(path) && !dir.exists(path), paste(label, "is missing:", path))
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
  expect_true(is.null(status) || status == 0L, paste("SHA-256 failed for", path))
  sub("[[:space:]].*$", "", output[[1L]])
}

result_dir_argument <- get_arg("--result-dir", "")
design <- match.arg(
  get_arg("--design", "balanced"),
  c("balanced", "prior_geometry")
)
require_gate <- as_flag(get_arg(
  "--require-calibration-gate",
  if (design == "prior_geometry") "true" else "false"
))
expect_true(nzchar(result_dir_argument), "--result-dir is required.")
result_dir <- normalizePath(
  result_dir_argument,
  winslash = "/",
  mustWork = TRUE
)
result_id <- basename(result_dir)

expected_result_id <- if (design == "balanced") {
  paste0(
    "r3_real_genotype_one_per_gene_J6362_",
    "matched_functional_open_middle_3_12_support_contained_",
    "relative_clearance_main_effect_fashr0143_pilot5"
  )
} else {
  paste0(
    "r3_real_genotype_one_per_gene_J6362_",
    "matched_functional_open_middle_3_12_full_support_",
    "iwp1_geometry_mixture_relative_clearance_main_effect_",
    "fashr0143_pilot5"
  )
}
expected_schema <- if (design == "balanced") {
  "r3-fashr0143-manifest-v4"
} else {
  "r3-fashr0143-manifest-v5-prior-geometry"
}
expected_source_hashes <- if (design == "balanced") {
  c(
    r3_driver = "7f5461d888aac8faaf3425e9ed9861aa5a1e1eaac84c3c48bee42af42c159246",
    simulation_functions = "f928897aef52c7a6bafa521dca00a732d5cddb391f0f20688646663f39c535b9",
    real_genotype_helper = "c03c01a188503336a77793c96f4e2d3ac7e0cbd56f4028b552da4b2f88e6b9d7",
    genotype_cache = "81bbef5f323a0bab2ca993c782d8a9b7c63518b83c2cdb46ef7ed1d46f65af24",
    wrapper_core = "8db000832342ae606709002e99a24bfe96437c9e131275550ede740609f83062",
    runner = "998994daabfd701de696fa200a09b296ac60f775cb60b070f743b2052bb51c89"
  )
} else {
  c(
    r3_driver = "32a0c6e3fb99e355d6c0ab0d68acbf7818b01d12f4516236addbfda2b33d48f9",
    simulation_functions = "3c5a97612591817a7edbbd3c331b380fae1d3a2fb468b4a78c0aed13ee3b8db9",
    real_genotype_helper = "c03c01a188503336a77793c96f4e2d3ac7e0cbd56f4028b552da4b2f88e6b9d7",
    genotype_cache = "81bbef5f323a0bab2ca993c782d8a9b7c63518b83c2cdb46ef7ed1d46f65af24",
    wrapper_core = "c4410bdae7ea3a84ebb8417bda182879ea3705932ee1e54e3703d770b8e904e1",
    runner = "b4dc6b42a7b30f888717ea4b708412f6393919daa54ef8524faf9dde7f05312d"
  )
}
expected_genotype_digests <- c(
  `12345` = "526a7318aa2af901e09252f5a6ca3c46",
  `22345` = "517faa30d5218a956f1be84f2567369c",
  `32345` = "9b9be3205d7db54dac31763492bcb2eb",
  `42345` = "9ef73aa94a061df868b1a951fe495d9f",
  `52345` = "7dab159b8453e2f66188ae313bfbd611"
)
seed_list <- as.integer(names(expected_genotype_digests))
truth_mechanisms <- c("random_bspline", "raised_cosine")
truth_group_levels <- c(
  "early / switch",
  "early / non-switch",
  "middle / switch",
  "middle / non-switch",
  "late / switch",
  "late / non-switch"
)
expected_truth_group_counts <- if (design == "balanced") {
  stats::setNames(rep(212L, 6L), truth_group_levels)
} else {
  stats::setNames(c(185L, 184L, 267L, 267L, 185L, 184L), truth_group_levels)
}

expect_true(identical(result_id, expected_result_id), "Unexpected result directory name.")
complete_path <- require_file(file.path(result_dir, "complete.flag"), "Complete flag")
manifest_path <- require_file(file.path(result_dir, "manifest.rds"), "Manifest")
manifest <- readRDS(manifest_path)
expect_true(identical(manifest$schema_version, expected_schema), "Unexpected manifest schema.")
expect_true(identical(manifest$result_id, expected_result_id), "Unexpected manifest result ID.")
expect_true(
  identical(manifest$package_provenance$version, "0.1.43") &&
    identical(
      manifest$package_provenance$remote_sha,
      "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
    ),
  "Unexpected fashr package provenance."
)

observed_source_hashes <- unlist(
  manifest$source_provenance$sha256,
  use.names = TRUE
)
expect_true(
  identical(observed_source_hashes[names(expected_source_hashes)], expected_source_hashes),
  "Formal source hashes do not match the frozen contract."
)

artifact_hashes <- unlist(manifest$artifact_sha256, use.names = TRUE)
actual_files <- list.files(
  result_dir,
  recursive = TRUE,
  full.names = TRUE,
  all.files = FALSE
)
actual_files <- actual_files[file.info(actual_files)$isdir %in% FALSE]
actual_relative <- substring(actual_files, nchar(result_dir) + 2L)
actual_artifacts <- setdiff(actual_relative, c("manifest.rds", "complete.flag"))
expect_true(
  identical(sort(names(artifact_hashes)), sort(actual_artifacts)),
  "The artifact manifest is incomplete or contains extra entries."
)
observed_artifact_hashes <- vapply(
  file.path(result_dir, names(artifact_hashes)),
  sha256_file,
  character(1)
)
expect_true(
  identical(unname(observed_artifact_hashes), unname(artifact_hashes)),
  "At least one formal artifact SHA-256 does not match the manifest."
)

configuration <- readRDS(require_file(
  file.path(result_dir, "configuration.rds"),
  "Configuration"
))
expect_true(
  identical(configuration$output_id, expected_result_id) &&
    identical(configuration$J, 6362L) &&
    identical(configuration$n_donors, 19L) &&
    identical(configuration$seed_list, seed_list) &&
    identical(configuration$truth_mechanisms, truth_mechanisms) &&
    identical(configuration$middle_window, c(3, 12)) &&
    identical(configuration$middle_boundary, "open") &&
    identical(configuration$middle_expression, "3 < t < 12") &&
    identical(configuration$raised_cosine$center_ranges, if (design == "balanced") {
      list(
        early = c(0.5, 1.5),
        middle = c(4.5, 10.5),
        late = c(13.5, 14.5)
      )
    } else {
      list(
        early = c(1.5, 1.5),
        middle = c(4.5, 10.5),
        late = c(13.5, 13.5)
      )
    }) &&
    identical(configuration$expected_truth_group_counts, expected_truth_group_counts) &&
    identical(configuration$genotype_content_digests, expected_genotype_digests),
  "The formal R3 configuration failed its scientific contract."
)
if (design == "prior_geometry") {
  expect_true(
    identical(
      configuration$temporal_category_design,
      "canonical IWP1 prior geometry"
    ) &&
      identical(
        configuration$temporal_category_probs,
        c(early = 0.29, middle = 0.42, late = 0.29)
      ),
    "The prior-geometry configuration is missing its frozen mixture."
  )
} else {
  expect_true(
    is.null(configuration$temporal_category_probs),
    "The historical balanced cache unexpectedly has a temporal mixture override."
  )
}

replicate_paths <- unlist(lapply(truth_mechanisms, function(mechanism) {
  file.path(
    result_dir,
    "replicates",
    paste0(mechanism, "_seed_", seed_list, ".rds")
  )
}), use.names = FALSE)
expect_true(all(file.exists(replicate_paths)), "At least one formal replicate is missing.")
replicates <- lapply(replicate_paths, readRDS)
for (replicate in replicates) {
  expect_true(
    isTRUE(all.equal(replicate$configuration, configuration)) &&
      identical(
        replicate$genotype_digest,
        expected_genotype_digests[[as.character(replicate$seed)]]
      ) &&
      nrow(replicate$functional_alpha) == 320L &&
      nrow(replicate$functional_alpha_005) == 8L &&
      all(is.finite(replicate$functional_alpha$power)) &&
      all(is.finite(replicate$functional_alpha$empirical_fsr)),
    paste("Replicate contract failed for", replicate$truth_mechanism, replicate$seed)
  )
  observed_counts <- stats::setNames(
    as.integer(replicate$truth_group_counts$n_dynamic),
    as.character(replicate$truth_group_counts$truth_group)
  )[truth_group_levels]
  expect_true(
    identical(observed_counts, expected_truth_group_counts),
    paste("Truth-group counts failed for", replicate$truth_mechanism, replicate$seed)
  )
}

summary_dir <- file.path(result_dir, "summary")
summary_contracts <- c(
  all_replicate_functional_alpha_curves.csv = 3200L,
  all_replicate_functional_alpha005.csv = 80L,
  all_truth_group_counts.csv = 60L,
  all_replicate_pi0.csv = 20L,
  genotype_selection_summary.csv = 10L,
  truth_maf_balance.csv = 30L,
  functional_testing_mc_alpha_curve.csv = 640L,
  functional_testing_mc_alpha005_summary.csv = 16L,
  functional_testing_mc_pi0_summary.csv = 4L
)
for (file_name in names(summary_contracts)) {
  object <- read.csv(require_file(
    file.path(summary_dir, file_name),
    paste("Summary", file_name)
  ))
  expect_true(
    nrow(object) == summary_contracts[[file_name]],
    paste("Unexpected row count for", file_name)
  )
}

mc_alpha <- read.csv(file.path(
  summary_dir,
  "functional_testing_mc_alpha_curve.csv"
))
middle <- mc_alpha[
  mc_alpha$target == "middle" &
    mc_alpha$method == "FASH-IWP1-BF" &
    mc_alpha$alpha >= 0.05,
  ,
  drop = FALSE
]
middle$truth_mechanism <- ifelse(
  grepl("^r3a_", middle$scenario),
  "random_bspline",
  ifelse(grepl("^r3b_", middle$scenario), "raised_cosine", NA_character_)
)
expect_true(nrow(middle) == 62L && !anyNA(middle$truth_mechanism), "Invalid gate rows.")
gate <- do.call(rbind, lapply(split(middle, middle$truth_mechanism), function(rows) {
  excess <- rows$mean_empirical_fsr - rows$alpha
  index <- which.max(excess)
  data.frame(
    truth_mechanism = rows$truth_mechanism[[1L]],
    maximum_excess = excess[[index]],
    alpha_at_maximum = rows$alpha[[index]],
    empirical_fsr_at_maximum = rows$mean_empirical_fsr[[index]],
    passed = excess[[index]] <= 0.03,
    stringsAsFactors = FALSE
  )
}))

if (design == "prior_geometry") {
  stored_gate <- read.csv(require_file(
    file.path(result_dir, "scientific_validation.csv"),
    "Scientific validation"
  ))
  expect_true(
    identical(stored_gate$truth_mechanism, gate$truth_mechanism) &&
      max(abs(stored_gate$maximum_excess - gate$maximum_excess)) < 1e-12,
    "Stored and independently recomputed calibration gates disagree."
  )
}
if (require_gate) {
  expect_true(all(gate$passed), "The prespecified Middle calibration gate failed.")
}

complete_lines <- readLines(complete_path, warn = FALSE)
expected_center_line <- if (design == "balanced") {
  paste0(
    "raised_cosine_center_ranges=",
    "early:0.5,1.5;middle:4.5,10.5;late:13.5,14.5"
  )
} else {
  paste0(
    "raised_cosine_center_ranges=",
    "early:1.5,1.5;middle:4.5,10.5;late:13.5,13.5"
  )
}
expect_true(
  paste0("result_id=", expected_result_id) %in% complete_lines &&
    "fashr_version=0.1.43" %in% complete_lines &&
    "middle_definition=3 < t < 12" %in% complete_lines &&
    expected_center_line %in% complete_lines &&
    "replicates=10" %in% complete_lines,
  "The complete flag is missing required provenance."
)

selected_alpha <- c(0.05, 0.10, 0.15, 0.20)
selected <- middle[
  vapply(
    middle$alpha,
    function(alpha) any(abs(alpha - selected_alpha) < 1e-12),
    logical(1)
  ),
  c(
    "truth_mechanism", "alpha", "mean_discoveries",
    "mean_estimated_fsr", "mean_empirical_fsr",
    "empirical_fsr_ci_lower", "empirical_fsr_ci_upper"
  ),
  drop = FALSE
]

cat("Validated formal R3 cache:", result_dir, "\n")
cat("Design:", design, "\n")
cat("Middle calibration gate:\n")
print(gate, row.names = FALSE)
cat("\nSelected BF Middle calibration rows:\n")
print(selected, row.names = FALSE)
