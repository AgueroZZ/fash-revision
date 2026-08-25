#!/usr/bin/env Rscript

# Verify the proposed single-lobe R3B geometry.  The test is intentionally
# independent of FASH fitting: it establishes that the generated deviation,
# rather than only its sampled center, lies strictly within its target window.

options(stringsAsFactors = FALSE)

if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
  workflowr_root <- "."
} else if (file.exists(
  "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
)) {
  workflowr_root <- "coderepo-local"
} else {
  stop("Could not locate simulation_functions.R.")
}
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "shared",
  "simulation_functions.R"
))

center_ranges <- list(
  early = c(1.5, 1.5),
  middle = c(4.5, 10.5),
  late = c(13.5, 13.5)
)
half_width <- 1.2

simulation <- simulate_matched_functional_effect_set(
  n_variants = 300L,
  truth_mechanism = "raised_cosine",
  time_grid = 0:15,
  evaluation_grid = seq(0, 15, by = 0.1),
  cosine_width_half = half_width,
  cosine_spike_counts = 1L,
  cosine_center_ranges = center_ranges,
  middle_window = c(3, 12),
  middle_boundary = "open",
  location_truth_margin = 0.10,
  location_truth_min_range_fraction = 0.10,
  non_switch_min_abs = 0.10,
  non_switch_min_range_fraction = 0.10,
  seed = 98765L
)

dynamic <- simulation$unit_info[
  simulation$unit_info$effect_class == "dynamic_bspline",
  ,
  drop = FALSE
]
primary_center <- vapply(
  dynamic$peak_centers,
  function(x) x[[1L]],
  numeric(1)
)
support_left <- primary_center - half_width
support_right <- primary_center + half_width
strictly_contained <-
  (dynamic$time_group == "early" & support_left > 0 & support_right < 3) |
  (dynamic$time_group == "middle" & support_left > 3 & support_right < 12) |
  (dynamic$time_group == "late" & support_left > 12 & support_right < 15)

if (!all(strictly_contained)) {
  print(cbind(
    dynamic[!strictly_contained, c("time_group", "truth_group")],
    support_left = support_left[!strictly_contained],
    support_right = support_right[!strictly_contained]
  ))
  stop("The proposed R3B primary-lobe supports are not target-contained.")
}

expected_support <- data.frame(
  time_group = c("early", "middle-left", "middle-right", "late"),
  support_left = c(0.3, 3.3, 9.3, 12.3),
  support_right = c(2.7, 5.7, 11.7, 14.7)
)
print(expected_support, row.names = FALSE)
cat(
  "Verified strict target containment for ",
  nrow(dynamic),
  " dynamic single-lobe R3B truths.\n",
  sep = ""
)
