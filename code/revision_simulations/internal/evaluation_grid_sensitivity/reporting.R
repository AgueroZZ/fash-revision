# Load and validate the current all-category numerical-sensitivity cache.

r4_sensitivity_cache_dir <- file.path(
  "output",
  "revision_simulations",
  "internal",
  "evaluation_grid_mc_all_category_current_fashr0143_20260826"
)

configuration <- readRDS(file.path(r4_sensitivity_cache_dir, "configuration.rds"))
population <- read.csv(
  file.path(r4_sensitivity_cache_dir, "baseline_discovery_population.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
population_summary <- read.csv(
  file.path(r4_sensitivity_cache_dir, "population_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
pair_lfsr <- read.csv(
  file.path(r4_sensitivity_cache_dir, "pair_lfsr_by_setting.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
comparison_summary <- read.csv(
  file.path(r4_sensitivity_cache_dir, "comparison_summary.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
runtime <- read.csv(
  file.path(r4_sensitivity_cache_dir, "runtime.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
cache_validation <- read.csv(
  file.path(r4_sensitivity_cache_dir, "validation.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

category_order <- c("early", "middle", "late", "switch")
category_labels <- c(
  early = "Early",
  middle = "Middle",
  late = "Late",
  switch = "Switch"
)
expected_pair_counts <- c(
  early = 126L,
  middle = 58L,
  late = 21L,
  switch = 981L
)
expected_gene_counts <- c(
  early = 8L,
  middle = 14L,
  late = 11L,
  switch = 250L
)
expected_setting_ids <- c(
  "grid_0p15_M3000",
  "grid_0p10_M3000",
  "grid_0p05_M3000",
  "grid_0p10_M2000",
  "grid_0p10_M5000"
)
expected_comparison_ids <- c(
  "grid_0p10_vs_0p15_M3000",
  "grid_0p10_vs_0p05_M3000",
  "M2000_vs_M3000_grid_0p10",
  "M3000_vs_M5000_grid_0p10"
)

current_input_md5 <- unname(tools::md5sum(configuration$input_paths))
if (
  configuration$seed != 20260826L ||
  configuration$num_cores != 2L ||
  configuration$alpha != 0.05 ||
  configuration$switch_threshold != 0.25 ||
  configuration$middle_definition != "open_3_12" ||
  configuration$category_pair_count != sum(expected_pair_counts) ||
  configuration$unique_pair_count != 1178L ||
  !identical(configuration$category_order, category_order) ||
  !identical(
    unname(configuration$expected_pair_counts),
    unname(expected_pair_counts)
  ) ||
  !identical(
    unname(configuration$expected_gene_counts),
    unname(expected_gene_counts)
  ) ||
  length(current_input_md5) != length(configuration$input_md5) ||
  anyNA(current_input_md5) ||
  !identical(current_input_md5, configuration$input_md5) ||
  !setequal(pair_lfsr$setting_id, expected_setting_ids) ||
  !identical(
    unique(comparison_summary$comparison_id),
    expected_comparison_ids
  ) ||
  nrow(population) != sum(expected_pair_counts) ||
  nrow(pair_lfsr) != sum(expected_pair_counts) *
    length(expected_setting_ids) ||
  nrow(comparison_summary) != length(category_order) *
    length(expected_comparison_ids) ||
  anyDuplicated(population[c("category", "pair_id")]) ||
  anyDuplicated(pair_lfsr[c("category", "pair_id", "setting_id")]) ||
  any(!is.finite(pair_lfsr$lfsr)) ||
  any(pair_lfsr$lfsr < 0 | pair_lfsr$lfsr > 1) ||
  !all(cache_validation$passed)
) {
  stop("The current all-category numerical-sensitivity cache failed validation.")
}

population_summary <- population_summary[
  match(category_order, population_summary$category),
  ,
  drop = FALSE
]
comparison_summary <- comparison_summary[
  order(
    match(comparison_summary$comparison_id, expected_comparison_ids),
    match(comparison_summary$category, category_order)
  ),
  ,
  drop = FALSE
]

grid_comparisons <- comparison_summary[
  comparison_summary$comparison_id %in% expected_comparison_ids[1:2],
  ,
  drop = FALSE
]
mc_2000_vs_3000 <- comparison_summary[
  comparison_summary$comparison_id == "M2000_vs_M3000_grid_0p10",
  ,
  drop = FALSE
]
mc_3000_vs_5000 <- comparison_summary[
  comparison_summary$comparison_id == "M3000_vs_M5000_grid_0p10",
  ,
  drop = FALSE
]
sampling_runtime_seconds <- runtime$elapsed_seconds[
  runtime$stage == "posterior_sampling_and_evaluation"
]

# Fail the render if the cache no longer supports the formal page's conclusion.
if (
  nrow(grid_comparisons) != 8L ||
  nrow(mc_2000_vs_3000) != 4L ||
  nrow(mc_3000_vs_5000) != 4L ||
  max(grid_comparisons$mean_absolute_change) >= 0.002 ||
  max(grid_comparisons$maximum_absolute_change) > 0.01 ||
  max(mc_3000_vs_5000$mean_absolute_change) >= 0.0024 ||
  min(mc_3000_vs_5000$fraction_within_0p01) < 0.998 ||
  max(abs(mc_3000_vs_5000$mean_signed_change)) >= 0.0007 ||
  length(sampling_runtime_seconds) != 1L ||
  !is.finite(sampling_runtime_seconds)
) {
  stop("The current cache no longer supports the formal R4 conclusion.")
}
