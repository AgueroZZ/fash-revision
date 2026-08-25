source("code/revision_simulations/appendix_b/reporting.R")

cache <- load_appendix_b_reporting_cache()
stopifnot(
  identical(cache$manifest$result_id, "appendix_b_fashr0143"),
  identical(cache$manifest$package_provenance$version, "0.1.43"),
  identical(
    cache$manifest$package_provenance$remote_sha,
    "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"
  ),
  nrow(cache$grid_summary) == 92L,
  nrow(cache$focused_alpha) == 80L,
  nrow(cache$boundary_warning_rows) == 7L,
  all(cache$boundary_warning_rows$raw_pi0_iwp1 == 0)
)

grid_long <- prepare_appendix_b_grid_long(cache$grid_summary)
stopifnot(
  nrow(grid_long) == 368L,
  !anyDuplicated(grid_long[c(
    "setting", "rho_dynamic", "order", "stage"
  )]),
  all(is.finite(grid_long$true_pi0)),
  all(is.finite(grid_long$estimated_pi0)),
  all(grid_long$true_pi0 >= 0 & grid_long$true_pi0 <= 1),
  all(grid_long$estimated_pi0 >= 0 & grid_long$estimated_pi0 <= 1),
  sum(grid_long$boundary_warning) == 7L,
  all(
    grid_long$true_pi0[grid_long$order == "IWP1"] ==
      1 - grid_long$rho_dynamic[grid_long$order == "IWP1"]
  ),
  all(
    grid_long$true_pi0[grid_long$order == "IWP2"] ==
      1 - grid_long$rho_nonlinear[grid_long$order == "IWP2"]
  )
)

alpha005 <- prepare_appendix_b_alpha005(cache$focused_alpha)
stopifnot(
  nrow(alpha005) == 4L,
  identical(alpha005$discoveries, c(210L, 205L, 106L, 100L)),
  identical(alpha005$false_discoveries, c(12L, 9L, 7L, 2L)),
  isTRUE(all.equal(
    alpha005$power,
    c(0.825, 0.816666666666667, 0.825, 0.816666666666667),
    tolerance = 1e-12
  ))
)

alpha_long <- prepare_appendix_b_alpha_long(cache$focused_alpha)
stopifnot(
  nrow(alpha_long) == 160L,
  !anyDuplicated(alpha_long[c("order", "stage", "alpha", "metric")]),
  all(alpha_long$value >= 0 & alpha_long$value <= 1)
)

prior_table <- prepare_appendix_b_prior_table(cache$focused_example)
stopifnot(
  nrow(prior_table) == 4L,
  !anyDuplicated(prior_table[c("order", "stage")]),
  isTRUE(all.equal(prior_table$true_pi0, c(0.8, 0.8, 0.9, 0.9))),
  isTRUE(all.equal(
    prior_table$estimated_pi0,
    c(0.781644420785063, 0.853333333333333, 0.849507291710244, 0.9275),
    tolerance = 1e-12
  ))
)

for (order_name in c("IWP1", "IWP2")) {
  for (stage_name in c("Raw EB", "BF updated")) {
    ranked <- prepare_appendix_b_cumulative_lfdr(
      cache$focused_example,
      order = order_name,
      stage = stage_name,
      alpha = 0.05
    )
    expected_calls <- alpha005$discoveries[
      alpha005$order == order_name & alpha005$stage == stage_name
    ]
    stopifnot(
      nrow(ranked) == 1200L,
      identical(ranked$rank, seq_len(1200L)),
      !anyDuplicated(ranked$unit_id),
      all(diff(ranked$lfdr) >= 0),
      all(diff(ranked$cumulative_lfdr) >= -1e-14),
      sum(ranked$called_at_alpha) == expected_calls
    )
  }
}

trajectories <- prepare_appendix_b_trajectory_examples(
  cache$focused_example,
  examples_per_class = 3L
)
stopifnot(
  length(unique(trajectories$unit_id)) == 9L,
  all(table(trajectories$class) > 0L),
  all(is.finite(trajectories$estimate)),
  all(is.finite(trajectories$standard_error)),
  all(trajectories$standard_error > 0),
  all(is.finite(trajectories$truth))
)

cat("Appendix B reporting contract test passed.\n")
