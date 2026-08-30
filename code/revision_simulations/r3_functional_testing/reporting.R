# Load, validate, and format cached results for the R3 workflowr report.

mc_output_id <-
  paste0(
    "r3_real_genotype_one_per_gene_J6362_",
    paste0(
      paste0(
        "matched_functional_open_middle_3_12_center_aligned_",
        "iwp1_geometry_mixture_"
      ),
      paste0(
        "relative_location_clearance_full_universe_",
        "paired_posterior_fashr0143_pilot5"
      )
    )
  )
mc_output_dir <- file.path(
  "output",
  "revision_simulations",
  "mc",
  mc_output_id
)
summary_dir <- file.path(mc_output_dir, "summary")
r3_replicate_alpha_path <- file.path(
  summary_dir,
  "all_replicate_functional_alpha_curves.csv"
)
r3_manifest_path <- file.path(mc_output_dir, "manifest.rds")
r3_complete_flag_path <- file.path(mc_output_dir, "complete.flag")
r3_scientific_validation_path <- file.path(
  mc_output_dir,
  "scientific_validation.csv"
)
r1_r2_output_dir <- file.path(
  "output",
  "revision_simulations",
  "mc",
  "r1_r2_fashr0143"
)
r1_r2_manifest_path <- file.path(r1_r2_output_dir, "manifest.rds")
shared_genotype_cache_path <- file.path(
  "output",
  "revision_simulations",
  "shared",
  "real_genotype_one_per_gene_J6362_pilot5",
  "genotype_samples.rds"
)
required_cache_paths <- c(
  r3_manifest_path,
  r3_complete_flag_path,
  r3_scientific_validation_path,
  file.path(mc_output_dir, "configuration.rds"),
  file.path(mc_output_dir, "example_curves.rds"),
  r3_replicate_alpha_path,
  file.path(summary_dir, "functional_testing_mc_alpha_curve.csv"),
  file.path(summary_dir, "functional_testing_mc_alpha005_summary.csv"),
  file.path(summary_dir, "functional_testing_mc_pi0_summary.csv"),
  file.path(summary_dir, "genotype_selection_summary.csv"),
  file.path(summary_dir, "truth_maf_balance.csv"),
  file.path(summary_dir, "all_truth_group_counts.csv"),
  r1_r2_manifest_path,
  shared_genotype_cache_path
)
if (any(!file.exists(required_cache_paths))) {
  stop("The formal real-genotype R3 cache is incomplete.")
}
r3_manifest <- readRDS(r3_manifest_path)
r3_completion <- readLines(r3_complete_flag_path, warn = FALSE)
scientific_validation <- read.csv(
  r3_scientific_validation_path,
  stringsAsFactors = FALSE
)
r1_r2_manifest <- readRDS(r1_r2_manifest_path)
shared_genotype_cache <- readRDS(shared_genotype_cache_path)
configuration <- readRDS(file.path(mc_output_dir, "configuration.rds"))
example_curves <- readRDS(file.path(
  mc_output_dir,
  "example_curves.rds"
))
mc_alpha <- read.csv(
  file.path(summary_dir, "functional_testing_mc_alpha_curve.csv"),
  stringsAsFactors = FALSE
)
mc_alpha_replicates <- read.csv(
  r3_replicate_alpha_path,
  stringsAsFactors = FALSE
)
mc_alpha_005 <- read.csv(
  file.path(
    summary_dir,
    "functional_testing_mc_alpha005_summary.csv"
  ),
  stringsAsFactors = FALSE
)
mc_pi0 <- read.csv(
  file.path(summary_dir, "functional_testing_mc_pi0_summary.csv"),
  stringsAsFactors = FALSE
)
genotype_selection_summary <- read.csv(
  file.path(summary_dir, "genotype_selection_summary.csv"),
  stringsAsFactors = FALSE
)
truth_maf_balance <- read.csv(
  file.path(summary_dir, "truth_maf_balance.csv"),
  stringsAsFactors = FALSE
)
truth_group_counts <- read.csv(
  file.path(summary_dir, "all_truth_group_counts.csv"),
  stringsAsFactors = FALSE
)

method_order <- c("FASH-IWP1-Raw", "FASH-IWP1-BF")
target_order <- c("early", "middle", "late", "switch")
mechanism_order <- c("random_bspline", "raised_cosine")
mechanism_labels <- c(
  random_bspline = "R3A: broad random B-spline",
  raised_cosine = "R3B: compact raised cosine"
)
expected_seed_list <- c(12345L, 22345L, 32345L, 42345L, 52345L)
expected_class_probs <- c(
  dynamic_bspline = 0.20,
  constant = 0.40,
  zero = 0.40
)
expected_class_counts <- c(
  dynamic_bspline = 1272L,
  constant = 2545L,
  zero = 2545L
)
expected_genotype_content_md5 <- c(
  `12345` = "526a7318aa2af901e09252f5a6ca3c46",
  `22345` = "517faa30d5218a956f1be84f2567369c",
  `32345` = "9b9be3205d7db54dac31763492bcb2eb",
  `42345` = "9ef73aa94a061df868b1a951fe495d9f",
  `52345` = "7dab159b8453e2f66188ae313bfbd611"
)
expected_fashr_remote_sha <-
  "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
expected_r3_source_sha256 <- c(
  r3_driver =
    "ca3f786ab11749b12f17e9799006a49c4264c623c3131dd18153f33daeb9da18",
  simulation_functions =
    "45267b0884168e5ae33cc4f14e3f05b711d961b65bf9c6fbd880e748de064a6e",
  real_genotype_helper =
    "c03c01a188503336a77793c96f4e2d3ac7e0cbd56f4028b552da4b2f88e6b9d7",
  genotype_cache =
    "81bbef5f323a0bab2ca993c782d8a9b7c63518b83c2cdb46ef7ed1d46f65af24",
  temporal_mixture_contract =
    "8e2b4c527b6c17fd2645b9d46e47a7b99410d2e9be1db1a38d541380d91ad723",
  wrapper_core =
    "2442bb89842c9078972210b0764f77ac9fdfb0d1e5155ce958ce503ae9f73b7a",
  runner =
    "2b6a64770d4039e92ff7cb9ce84c73148c81e6b00ace525e70234ea9e42a3d18"
)
expected_temporal_category_probs <- stats::setNames(
  c(0.29, 0.42, 0.29),
  c("early", "middle", "late")
)
expected_truth_group_counts <- setNames(
  c(185L, 184L, 267L, 267L, 185L, 184L),
  c(
    "early / switch",
    "early / non-switch",
    "middle / switch",
    "middle / non-switch",
    "late / switch",
    "late / non-switch"
  )
)

if (!identical(
      r3_manifest$schema_version,
      paste0(
        "r3-fashr0143-manifest-v9-full-universe-functional-",
        "iwp1-temporal-mixture"
      )
    ) ||
    !identical(r3_manifest$result_id, mc_output_id) ||
    !identical(r3_manifest$package_provenance$package, "fashr") ||
    !identical(r3_manifest$package_provenance$version, "0.1.43") ||
    !identical(
      r3_manifest$package_provenance$remote_sha,
      expected_fashr_remote_sha
    ) ||
    !all(c(
      paste0("result_id=", mc_output_id),
      "fashr_version=0.1.43",
      paste0("fashr_remote_sha=", expected_fashr_remote_sha),
      "middle_definition=3 < t < 12",
      "temporal_category_probs=early:0.29;middle:0.42;late:0.29",
      paste0(
        "truth_group_counts=early / switch:185;",
        "early / non-switch:184;middle / switch:267;",
        "middle / non-switch:267;late / switch:185;",
        "late / non-switch:184"
      ),
      paste0(
        "temporal_mixture_contract_sha256=",
        expected_r3_source_sha256[["temporal_mixture_contract"]]
      ),
      paste0(
        "raised_cosine_center_ranges=",
        "early:1.5,2.5;middle:4.5,10.5;late:12.5,13.5"
      ),
      "raised_cosine_half_width=1.5",
      "location_truth_margin=0.1",
      "location_truth_min_range_fraction=0.1",
      "functional_posterior_pairing=common_random_seed_raw_bf",
      "functional_candidate_scope=full_universe",
      "functional_candidate_universe_size=6362",
      "calibration_gate_passed=FALSE",
      "replicates=10"
    ) %in% r3_completion) ||
    !identical(
      unlist(
        r3_manifest$source_provenance$sha256[
          names(expected_r3_source_sha256)
        ],
        use.names = TRUE
      ),
      expected_r3_source_sha256
    ) ||
    !identical(configuration$output_id, mc_output_id) ||
    !identical(configuration$package_provenance$version, "0.1.43") ||
    !identical(
      configuration$package_provenance$remote_sha,
      expected_fashr_remote_sha
    ) ||
    !identical(
      configuration$genotype_digest_method,
      "fash-genotype-content-md5-v1"
    ) ||
    !identical(
      configuration$genotype_content_digests,
      expected_genotype_content_md5
    ) ||
    !isTRUE(all.equal(configuration$J, 6362L)) ||
    !isTRUE(all.equal(configuration$n_donors, 19L)) ||
    !isTRUE(all.equal(configuration$n_covariates, 5L)) ||
    !isTRUE(all.equal(configuration$time_grid, 0:15)) ||
    !isTRUE(all.equal(configuration$evaluation_grid, seq(0, 15, by = 0.1))) ||
    !isTRUE(all.equal(configuration$middle_window, c(3, 12))) ||
    !identical(configuration$middle_boundary, "open") ||
    !isTRUE(all.equal(configuration$middle_grid, seq(3.1, 11.9, by = 0.1))) ||
    !identical(configuration$middle_expression, "3 < t < 12") ||
    !isTRUE(all.equal(configuration$expression_noise_sd, 1)) ||
    !isTRUE(all.equal(configuration$covariate_effect_sd, 0.5)) ||
    !isTRUE(all.equal(configuration$intercept_sd, 0)) ||
    !isTRUE(all.equal(configuration$dynamic_main_effect_sd, 1)) ||
    !isTRUE(all.equal(configuration$class_probs, expected_class_probs)) ||
    !isTRUE(all.equal(
      configuration$expected_class_counts,
      expected_class_counts
    )) ||
    !identical(
      configuration$expected_truth_group_counts,
      expected_truth_group_counts
    ) ||
    !identical(
      configuration$temporal_category_design,
      "user-specified temporal-category probabilities"
    ) ||
    !isTRUE(all.equal(
      configuration$temporal_category_probs,
      expected_temporal_category_probs,
      tolerance = 1e-12,
      check.attributes = TRUE
    )) ||
    !identical(configuration$truth_mechanisms, mechanism_order) ||
    !isTRUE(all.equal(configuration$random_bspline$amplitude, 2)) ||
    !isTRUE(all.equal(configuration$random_bspline$df, 6)) ||
    !isTRUE(all.equal(configuration$raised_cosine$width_half, 1.5)) ||
    !identical(
      configuration$raised_cosine$center_ranges,
      list(
        early = c(1.5, 2.5),
        middle = c(4.5, 10.5),
        late = c(12.5, 13.5)
      )
    ) ||
    !isTRUE(all.equal(
      configuration$raised_cosine$spike_counts,
      1:3
    )) ||
    !isTRUE(all.equal(configuration$switch_threshold, 0.25)) ||
    !isTRUE(all.equal(configuration$location_truth_margin, 0.10)) ||
    !isTRUE(all.equal(
      configuration$location_truth_min_range_fraction,
      0.10
    )) ||
    !isTRUE(all.equal(configuration$switch_truth_margin, 0.10)) ||
    !isTRUE(all.equal(configuration$non_switch_min_abs, 0.10)) ||
    !isTRUE(all.equal(
      configuration$non_switch_min_range_fraction,
      0.10
    )) ||
    !identical(
      configuration$functional_posterior_pairing,
      "common_random_seed_raw_bf"
    ) ||
    !identical(configuration$functional_candidate_scope, "full_universe") ||
    !identical(configuration$functional_candidate_universe_size, 6362L) ||
    !identical(configuration$genotype_source, "paper-derived YRI DS dosage") ||
    !identical(
      configuration$genotype_selection_rule,
      "one uniformly sampled tested variant per gene"
    ) ||
    !identical(configuration$genotype_dosage_field, "DS") ||
    !isTRUE(all.equal(configuration$genotype_maf_min, 0.10)) ||
    !identical(length(configuration$genotype_sample_ids), 19L) ||
    !identical(as.integer(configuration$seed_list), expected_seed_list) ||
    !identical(names(example_curves), mechanism_order)) {
  stop("The matched functional-testing cache has unexpected settings.")
}
required_validation_columns <- c(
  "truth_mechanism",
  "method",
  "target",
  "alpha_min",
  "alpha_max",
  "maximum_excess",
  "prespecified_maximum_excess",
  "passed"
)
validation_passed <- stats::setNames(
  scientific_validation$passed,
  scientific_validation$truth_mechanism
)
expected_validation_excess <- c(
  raised_cosine = -0.0100538301311181,
  random_bspline = 0.058188878457505
)
expected_validation_alpha <- c(
  raised_cosine = 0.07,
  random_bspline = 0.20
)
expected_validation_passed <- c(
  raised_cosine = TRUE,
  random_bspline = FALSE
)
if (!all(required_validation_columns %in% names(scientific_validation)) ||
    nrow(scientific_validation) != 2L ||
    !setequal(
      scientific_validation$truth_mechanism,
      mechanism_order
    ) ||
    any(scientific_validation$method != "FASH-IWP1-BF") ||
    any(scientific_validation$target != "middle") ||
    any(abs(scientific_validation$alpha_min - 0.05) > 1e-12) ||
    any(abs(scientific_validation$alpha_max - 0.20) > 1e-12) ||
    any(abs(
      scientific_validation$prespecified_maximum_excess - 0.03
    ) > 1e-12) ||
    !identical(
      validation_passed[names(expected_validation_passed)],
      expected_validation_passed
    ) ||
    !isTRUE(all.equal(
      stats::setNames(
        scientific_validation$maximum_excess,
        scientific_validation$truth_mechanism
      )[names(expected_validation_excess)],
      expected_validation_excess,
      tolerance = 1e-12,
      check.attributes = TRUE
    )) ||
    !isTRUE(all.equal(
      stats::setNames(
        scientific_validation$alpha_at_maximum_excess,
        scientific_validation$truth_mechanism
      )[names(expected_validation_alpha)],
      expected_validation_alpha,
      tolerance = 1e-12,
      check.attributes = TRUE
    ))) {
  stop("The R3 scientific-validation record is invalid.")
}
if (!all(method_order %in% mc_alpha$method) ||
    !all(target_order %in% mc_alpha$target) ||
    any(mc_alpha$n_replications != length(configuration$seed_list)) ||
    any(mc_alpha$candidate_scope != "full_universe") ||
    any(mc_alpha$mean_candidate_count != 6362) ||
    any(mc_alpha$mean_first_stage_null_calls != 0) ||
    !all(mechanism_order %in% mc_pi0$truth_mechanism)) {
  stop("The matched functional-testing summaries are incomplete.")
}

replicate_curve_columns <- c(
  "scenario",
  "target",
  "method",
  "alpha",
  "candidate_scope",
  "candidate_count",
  "first_stage_null_calls",
  "power",
  "empirical_fsr",
  "seed",
  "truth_mechanism"
)
replicate_curve_keys <- c("scenario", "target", "method", "alpha")
expected_alpha_grid <- seq(0.005, 0.20, by = 0.005)
expected_replicate_rows <-
  length(mechanism_order) *
  length(target_order) *
  length(method_order) *
  length(expected_alpha_grid) *
  length(expected_seed_list)
if (!all(replicate_curve_columns %in% names(mc_alpha_replicates)) ||
    nrow(mc_alpha_replicates) != expected_replicate_rows ||
    anyDuplicated(mc_alpha_replicates[c(replicate_curve_keys, "seed")]) ||
    !setequal(mc_alpha_replicates$seed, expected_seed_list) ||
    !setequal(mc_alpha_replicates$target, target_order) ||
    !setequal(mc_alpha_replicates$method, method_order) ||
    !setequal(mc_alpha_replicates$truth_mechanism, mechanism_order) ||
    any(mc_alpha_replicates$candidate_scope != "full_universe") ||
    any(mc_alpha_replicates$candidate_count != 6362L) ||
    any(mc_alpha_replicates$first_stage_null_calls != 0L) ||
    !isTRUE(all.equal(
      sort(unique(mc_alpha_replicates$alpha)),
      expected_alpha_grid
    ))) {
  stop("The replicate-level functional alpha curves are incomplete.")
}

functional_curve_key <- function(x) {
  paste(
    x$scenario,
    x$target,
    x$method,
    sprintf("%.17g", x$alpha),
    sep = "\r"
  )
}
replicate_groups <- split(
  seq_len(nrow(mc_alpha_replicates)),
  functional_curve_key(mc_alpha_replicates)
)
mc_pointwise_ranges <- do.call(
  rbind,
  lapply(replicate_groups, function(row_index) {
    group <- mc_alpha_replicates[row_index, , drop = FALSE]
    data.frame(
      scenario = group$scenario[1],
      target = group$target[1],
      method = group$method[1],
      alpha = group$alpha[1],
      n_replications = nrow(group),
      mean_power_from_replicates = mean(group$power),
      power_replication_min = min(group$power),
      power_replication_max = max(group$power),
      mean_empirical_fsr_from_replicates = mean(group$empirical_fsr),
      empirical_fsr_replication_min = min(group$empirical_fsr),
      empirical_fsr_replication_max = max(group$empirical_fsr),
      stringsAsFactors = FALSE
    )
  })
)
rownames(mc_pointwise_ranges) <- NULL
summary_range_index <- match(
  functional_curve_key(mc_alpha),
  functional_curve_key(mc_pointwise_ranges)
)
if (anyNA(summary_range_index) ||
    any(mc_pointwise_ranges$n_replications != length(expected_seed_list)) ||
    !isTRUE(all.equal(
      mc_alpha$mean_power,
      mc_pointwise_ranges$mean_power_from_replicates[summary_range_index],
      tolerance = 1e-12
    )) ||
    !isTRUE(all.equal(
      mc_alpha$mean_empirical_fsr,
      mc_pointwise_ranges$mean_empirical_fsr_from_replicates[
        summary_range_index
      ],
      tolerance = 1e-12
    ))) {
  stop("The replicate-level curves do not reproduce the saved summaries.")
}
pointwise_range_columns <- c(
  "power_replication_min",
  "power_replication_max",
  "empirical_fsr_replication_min",
  "empirical_fsr_replication_max"
)
for (column in pointwise_range_columns) {
  mc_alpha[[column]] <- mc_pointwise_ranges[[column]][summary_range_index]
}

if (!identical(r1_r2_manifest$schema_version, "r1-r2-fashr0143-manifest-v1") ||
    !identical(r1_r2_manifest$result_id, "r1_r2_fashr0143") ||
    !identical(r1_r2_manifest$configuration$J, 6362L) ||
    !identical(r1_r2_manifest$package_provenance$version, "0.1.43") ||
    !identical(
      r1_r2_manifest$package_provenance$remote_sha,
      expected_fashr_remote_sha
    ) ||
    !identical(
      r1_r2_manifest$source_provenance$sha256$genotype_cache,
      r3_manifest$source_provenance$sha256$genotype_cache
    ) ||
    !identical(shared_genotype_cache$configuration$n_genes, 6362L) ||
    !identical(shared_genotype_cache$configuration$n_donors, 19L) ||
    !identical(shared_genotype_cache$configuration$seed_list, expected_seed_list)) {
  stop("The formal R1/R2/R3 real-genotype provenance is inconsistent.")
}

if (nrow(genotype_selection_summary) !=
      length(expected_seed_list) * length(mechanism_order) ||
    !setequal(genotype_selection_summary$seed, expected_seed_list) ||
    !setequal(genotype_selection_summary$truth_mechanism, mechanism_order) ||
    any(genotype_selection_summary$genes != 6362L) ||
    any(genotype_selection_summary$selected_pairs != 6362L) ||
    any(genotype_selection_summary$maf_min < 0.10) ||
    any(genotype_selection_summary$maf_max > 0.50) ||
    nrow(truth_maf_balance) !=
      length(expected_seed_list) * length(mechanism_order) *
        length(expected_class_probs) ||
    any(
      truth_maf_balance$n !=
        expected_class_counts[truth_maf_balance$effect_class]
    ) ||
    nrow(truth_group_counts) !=
      length(expected_seed_list) * length(mechanism_order) * 6L ||
    anyNA(expected_truth_group_counts[truth_group_counts$truth_group]) ||
    any(
      truth_group_counts$n_dynamic !=
        expected_truth_group_counts[truth_group_counts$truth_group]
    )) {
  stop("The R3 real-genotype sampling, MAF, or truth-count diagnostics are invalid.")
}

for (seed in expected_seed_list) {
  seed_name <- as.character(seed)
  r1_reference <- readRDS(file.path(
    r1_r2_output_dir,
    "replicates",
    "r1",
    paste0("seed_", seed, ".rds")
  ))
  r3_replicates <- lapply(mechanism_order, function(mechanism) {
    readRDS(file.path(
      mc_output_dir,
      "replicates",
      paste0(mechanism, "_seed_", seed, ".rds")
    ))
  })
  if (!all(vapply(
    r3_replicates,
    function(replicate) {
      identical(replicate$genotype_digest, expected_genotype_content_md5[[seed_name]]) &&
        identical(replicate$genotype_digest, r1_reference$genotype_digest) &&
        identical(replicate$selected_pair_keys, r1_reference$selected_pair_keys) &&
        identical(
          replicate$selected_pair_keys,
          shared_genotype_cache$samples[[seed_name]]$selection$pair_key
        )
    },
    logical(1)
  ))) {
    stop("The formal R1/R3 genotype pairing differs for seed ", seed, ".")
  }

  seed_sampling <- genotype_selection_summary[
    genotype_selection_summary$seed == seed,
    setdiff(names(genotype_selection_summary), "truth_mechanism"),
    drop = FALSE
  ]
  seed_maf <- truth_maf_balance[
    truth_maf_balance$seed == seed,
    setdiff(names(truth_maf_balance), "truth_mechanism"),
    drop = FALSE
  ]
  if (nrow(unique(seed_sampling)) != 1L ||
      nrow(unique(seed_maf)) != length(expected_class_probs)) {
    stop("R3A and R3B genotype diagnostics differ for seed ", seed, ".")
  }
}

r3a_alpha <- mc_alpha[grepl("^r3a_", mc_alpha$scenario), ]
r3b_alpha <- mc_alpha[grepl("^r3b_", mc_alpha$scenario), ]
r3a_alpha_005 <- mc_alpha_005[
  grepl("^r3a_", mc_alpha_005$scenario),
]
r3b_alpha_005 <- mc_alpha_005[
  grepl("^r3b_", mc_alpha_005$scenario),
]
plot_alpha_limits <- c(0, max(mc_alpha$alpha))
if (!isTRUE(all.equal(min(mc_alpha$alpha), 0.005)) ||
    !isTRUE(all.equal(plot_alpha_limits, c(0, 0.20))) ||
    !isTRUE(all.equal(range(r3a_alpha$alpha), c(0.005, 0.20))) ||
    !isTRUE(all.equal(range(r3b_alpha$alpha), c(0.005, 0.20)))) {
  stop("The full R3 alpha ranges are invalid.")
}

format_decimal <- function(x, digits = 3) {
  formatC(x, format = "f", digits = digits)
}

format_mc_interval <- function(mean, lower, upper, digits = 3) {
  paste0(
    format_decimal(mean, digits),
    " [",
    format_decimal(lower, digits),
    ", ",
    format_decimal(upper, digits),
    "]"
  )
}

format_functional_table <- function(summary_table) {
  table <- summary_table
  table$target <- factor(table$target, levels = target_order)
  table$method <- factor(table$method, levels = method_order)
  table <- table[order(table$target, table$method), , drop = FALSE]
  data.frame(
    Target = tools::toTitleCase(as.character(table$target)),
    Method = as.character(table$method),
    `Mean calls` = format_decimal(table$mean_discoveries, 1),
    `Power (95% MC CI)` = format_mc_interval(
      table$mean_power,
      table$power_ci_lower,
      table$power_ci_upper
    ),
    `Empirical FSR (95% MC CI)` = format_mc_interval(
      table$mean_empirical_fsr,
      table$empirical_fsr_ci_lower,
      table$empirical_fsr_ci_upper
    ),
    check.names = FALSE
  )
}

format_pi0_table <- function(summary_table) {
  summary_table$truth_mechanism <- factor(
    summary_table$truth_mechanism,
    levels = mechanism_order
  )
  summary_table$fit <- factor(
    summary_table$fit,
    levels = c("Raw", "BF-corrected")
  )
  summary_table <- summary_table[
    order(summary_table$truth_mechanism, summary_table$fit),
    ,
    drop = FALSE
  ]
  data.frame(
    Mechanism = unname(mechanism_labels[
      as.character(summary_table$truth_mechanism)
    ]),
    Fit = as.character(summary_table$fit),
    `Mean estimated pi0 (95% MC CI)` = format_mc_interval(
      summary_table$mean_estimated_pi0,
      summary_table$pi0_ci_lower,
      summary_table$pi0_ci_upper
    ),
    check.names = FALSE
  )
}

metric_at <- function(summary_table, target, method, metric) {
  row <- summary_table[
    summary_table$target == target &
      summary_table$method == method,
    ,
    drop = FALSE
  ]
  if (nrow(row) != 1L || !metric %in% names(row)) {
    stop("Could not extract a unique functional-testing metric.")
  }
  row[[metric]]
}

plot_truth_examples <- function(examples, mechanism_label) {
  group_order <- c(
    "early / switch",
    "early / non-switch",
    "middle / switch",
    "middle / non-switch",
    "late / switch",
    "late / non-switch"
  )
  panels <- unlist(examples[group_order], recursive = FALSE)
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(
    mfrow = c(6, 2),
    mar = c(3.5, 3.8, 3.0, 0.8),
    oma = c(0, 0, 2.2, 0)
  )
  for (panel_index in seq_along(panels)) {
    panel <- panels[[panel_index]]
    observed <- panel$observed
    target_window <- switch(
      panel$time_group,
      early = c(0, 3),
      middle = configuration$middle_window,
      late = c(12, 15)
    )
    window_color <- switch(
      panel$time_group,
      early = "#d7eaf7",
      middle = "#dff0d8",
      late = "#fde6c9"
    )
    y_limits <- range(
      observed$estimate - 2 * observed$se,
      observed$estimate + 2 * observed$se,
      panel$true_curve$true_effect
    )
    y_padding <- max(0.15, 0.06 * diff(y_limits))
    y_limits <- y_limits + c(-y_padding, y_padding)
    shape_suffix <- if (is.na(panel$spike_count)) {
      ""
    } else {
      paste0("; ", panel$spike_count, " peak",
        if (panel$spike_count == 1L) "" else "s"
      )
    }
    plot(
      observed$time,
      observed$estimate,
      type = "n",
      xlim = range(configuration$time_grid),
      ylim = y_limits,
      xlab = if (panel_index > length(panels) - 2L) "Time" else "",
      ylab = "Genetic effect",
      main = paste0(
        panel$truth_group,
        shape_suffix,
        "\n",
        panel$variant_id
      )
    )
    rect(
      target_window[1],
      y_limits[1],
      target_window[2],
      y_limits[2],
      col = window_color,
      border = NA
    )
    arrows(
      observed$time,
      observed$estimate - 2 * observed$se,
      observed$time,
      observed$estimate + 2 * observed$se,
      angle = 90,
      code = 3,
      length = 0.035,
      col = "gray45"
    )
    points(observed$time, observed$estimate, pch = 19)
    lines(
      panel$true_curve$time,
      panel$true_curve$true_effect,
      col = "#D55E00",
      lwd = 2.2
    )
    abline(h = 0, col = "gray75", lty = 3)
  }
  mtext(
    paste0(
      mechanism_label,
      ": estimates with two-SE bars and continuous truth"
    ),
    outer = TRUE,
    cex = 1.05
  )
}

format_integer <- function(x) {
  format(as.integer(x), big.mark = ",", scientific = FALSE, trim = TRUE)
}

genotype_sampling_display <- genotype_selection_summary[
  genotype_selection_summary$truth_mechanism == "random_bspline",
  ,
  drop = FALSE
]
genotype_sampling_display <- genotype_sampling_display[
  order(genotype_sampling_display$seed),
  ,
  drop = FALSE
]
genotype_sampling_table <- data.frame(
  Seed = genotype_sampling_display$seed,
  Genes = format_integer(genotype_sampling_display$genes),
  `Unique rsIDs` = format_integer(
    genotype_sampling_display$unique_variant_ids
  ),
  `Repeated cross-gene assignments` = format_integer(
    genotype_sampling_display$repeated_cross_gene_assignments
  ),
  `MAF range` = paste0(
    format_decimal(genotype_sampling_display$maf_min, 3),
    "-",
    format_decimal(genotype_sampling_display$maf_max, 3)
  ),
  `Median MAF` = format_decimal(
    genotype_sampling_display$maf_median,
    3
  ),
  check.names = FALSE
)

truth_maf_display <- truth_maf_balance[
  truth_maf_balance$truth_mechanism == "random_bspline",
  ,
  drop = FALSE
]
truth_maf_balance_table <- do.call(rbind, lapply(
  names(expected_class_probs),
  function(effect_class) {
    rows <- truth_maf_display[
      truth_maf_display$effect_class == effect_class,
      ,
      drop = FALSE
    ]
    data.frame(
      `Truth class` = c(
        dynamic_bspline = "Dynamic",
        constant = "Constant",
        zero = "Zero"
      )[[effect_class]],
      `Units per seed` = format_integer(expected_class_counts[[effect_class]]),
      `Mean MAF across seeds` = format_decimal(mean(rows$maf_mean), 3),
      `Maximum absolute MAF SMD` = format_decimal(
        max(abs(rows$standardized_mean_difference)),
        3
      ),
      check.names = FALSE
    )
  }
))
rownames(truth_maf_balance_table) <- NULL

r3a_table <- format_functional_table(r3a_alpha_005)
r3b_table <- format_functional_table(r3b_alpha_005)
pi0_table <- format_pi0_table(mc_pi0)
