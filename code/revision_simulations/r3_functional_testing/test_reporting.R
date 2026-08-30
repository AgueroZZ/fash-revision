# Regression checks for the promoted center-aligned R3 reporting cache.

source("code/revision_simulations/shared/simulation_functions.R")
source("code/revision_simulations/r3_functional_testing/reporting.R")

expected_result_id <- paste0(
  "r3_real_genotype_one_per_gene_J6362_",
  paste0(
    "matched_functional_open_middle_3_12_center_aligned_",
    "iwp1_geometry_mixture_"
  ),
  paste0(
    "relative_location_clearance_full_universe_",
    "paired_posterior_fashr0143_pilot5"
  )
)
expected_temporal_probs <- stats::setNames(
  c(0.29, 0.42, 0.29),
  c("early", "middle", "late")
)
expected_truth_counts <- stats::setNames(
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

stopifnot(
  identical(mc_output_id, expected_result_id),
  identical(
    r3_manifest$schema_version,
    paste0(
      "r3-fashr0143-manifest-v9-full-universe-functional-",
      "iwp1-temporal-mixture"
    )
  ),
  identical(
    configuration$temporal_category_design,
    "user-specified temporal-category probabilities"
  ),
  isTRUE(all.equal(
    configuration$temporal_category_probs,
    expected_temporal_probs,
    tolerance = 1e-12,
    check.attributes = TRUE
  )),
  identical(configuration$expected_truth_group_counts, expected_truth_counts),
  identical(
    configuration$functional_posterior_pairing,
    "common_random_seed_raw_bf"
  ),
  identical(configuration$functional_candidate_scope, "full_universe"),
  identical(configuration$functional_candidate_universe_size, 6362L),
  all(mc_alpha$candidate_scope == "full_universe"),
  all(mc_alpha$mean_candidate_count == 6362),
  all(mc_alpha$mean_first_stage_null_calls == 0),
  all(mc_alpha_replicates$candidate_scope == "full_universe"),
  all(mc_alpha_replicates$candidate_count == 6362L),
  all(mc_alpha_replicates$first_stage_null_calls == 0L),
  isTRUE(all.equal(configuration$location_truth_margin, 0.10)),
  isTRUE(all.equal(
    configuration$location_truth_min_range_fraction,
    0.10
  )),
  isTRUE(all.equal(plot_alpha_limits, c(0, 0.20))),
  isTRUE(all.equal(min(mc_alpha$alpha), 0.005)),
  isTRUE(all.equal(max(mc_alpha$alpha), 0.20)),
  isTRUE(all.equal(
    sort(unique(mc_alpha$alpha)),
    seq(0.005, 0.20, by = 0.005)
  )),
  !exists("display_alpha_max", inherits = FALSE),
  !exists("r3a_alpha_display", inherits = FALSE),
  !exists("r3b_alpha_display", inherits = FALSE),
  length(list.files(
    file.path(mc_output_dir, "replicates"),
    pattern = "[.]rds$"
  )) == 10L
)

pointwise_range_columns <- c(
  "power_replication_min",
  "power_replication_max",
  "empirical_fsr_replication_min",
  "empirical_fsr_replication_max"
)
stopifnot(
  nrow(mc_alpha_replicates) == 3200L,
  length(unique(mc_alpha_replicates$seed)) == 5L,
  all(pointwise_range_columns %in% names(mc_alpha)),
  all(mc_alpha$power_replication_min <= mc_alpha$mean_power),
  all(mc_alpha$mean_power <= mc_alpha$power_replication_max),
  all(
    mc_alpha$empirical_fsr_replication_min <=
      mc_alpha$mean_empirical_fsr
  ),
  all(
    mc_alpha$mean_empirical_fsr <=
      mc_alpha$empirical_fsr_replication_max
  ),
  any(
    mc_alpha$power_replication_min <
      mc_alpha$power_replication_max
  ),
  any(
    mc_alpha$empirical_fsr_replication_min <
      mc_alpha$empirical_fsr_replication_max
  )
)

plot_test_curve <- data.frame(
  target = rep("early", 4L),
  method = rep(c("FASH-IWP1-Raw", "FASH-IWP1-BF"), each = 2L),
  alpha = rep(c(0.05, 0.10), 2L),
  mean_power = c(0.10, 0.20, 0.08, 0.16),
  power_replication_min = c(0.08, 0.18, 0.06, 0.13),
  power_replication_max = c(0.12, 0.22, 0.10, 0.19),
  stringsAsFactors = FALSE
)
polygon_call_count <- 0L
polygon <- function(...) {
  polygon_call_count <<- polygon_call_count + 1L
}
plot_device <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_device)
plot_mc_functional_curve_grid(
  mc_curve = plot_test_curve,
  metric = "power",
  target_order = "early",
  interval_columns = c(
    "power_replication_min",
    "power_replication_max"
  ),
  interval_methods = "FASH-IWP1-BF"
)
invisible(grDevices::dev.off())
rm(polygon)
unlink(plot_device)
stopifnot(polygon_call_count == 1L)

validation_passed <- stats::setNames(
  scientific_validation$passed,
  scientific_validation$truth_mechanism
)
stopifnot(
  identical(validation_passed[["raised_cosine"]], TRUE),
  identical(validation_passed[["random_bspline"]], FALSE)
)

r3a_middle_alpha005 <- r3a_alpha_005[
  r3a_alpha_005$method == "FASH-IWP1-BF" &
    r3a_alpha_005$target == "middle",
  ,
  drop = FALSE
]
stopifnot(
  nrow(r3a_middle_alpha005) == 1L,
  isTRUE(all.equal(
    r3a_middle_alpha005$mean_empirical_fsr,
    0.0353398479772499,
    tolerance = 1e-12
  ))
)

message("Promoted R3 reporting tests passed.")
