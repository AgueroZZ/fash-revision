# Regression checks for the formal full-universe R3 functional calibration run.

shared_functions_path <- file.path(
  "code", "revision_simulations", "shared", "simulation_functions.R"
)
snapshot_functions_path <- file.path(
  "code", "revision_simulations", "r3_r4_fashr0143", "source_snapshots",
  "r3_full_universe_functional_simulation_functions.R"
)
driver_path <- file.path(
  "code", "revision_simulations", "r3_r4_fashr0143", "source_snapshots",
  "r3_full_universe_functional_driver.R"
)
core_path <- file.path(
  "code", "revision_simulations", "r3_r4_fashr0143",
  "22_r3_full_universe_functional_core.R"
)
wrapper_path <- file.path(
  "code", "revision_simulations", "r3_r4_fashr0143",
  "22_r3_full_universe_functional_fashr0143.R"
)
sbatch_path <- file.path(
  "code", "revision_simulations", "r3_r4_fashr0143",
  "22_r3_full_universe_functional_fashr0143.sbatch"
)

source(shared_functions_path)

observed_functional_indices <- integer()
get_fash_fdr_table <- function(fit) {
  data.frame(index = c(2L, 4L), FDR = c(0.01, 0.02))
}
compute_functional_lfsr <- function(fit,
                                    functionals,
                                    indices,
                                    smooth_var,
                                    num_cores = 1,
                                    seed = NULL) {
  observed_functional_indices <<- indices
  lfsr <- stats::setNames(
    seq(0.01, 0.04, length.out = length(indices)),
    as.character(indices)
  )
  stats::setNames(
    rep(list(lfsr), length(functionals)),
    names(functionals)
  )
}

true_functionals <- matrix(
  rep(c(1, -1, 1, -1, 1, -1), 4L),
  nrow = 6L,
  ncol = 4L,
  dimnames = list(NULL, c("early", "middle", "late", "switch"))
)
full_result <- evaluate_fash_functional_testing(
  fit = list(),
  true_functionals = true_functionals,
  evaluation_grid = seq(0, 15, by = 0.1),
  alpha_grid = c(0.05, 0.10),
  true_dynamic = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE),
  candidate_scope = "full_universe",
  middle_window = c(3, 12),
  middle_boundary = "open"
)
stopifnot(
  identical(observed_functional_indices, 1:6),
  all(full_result$alpha_curve$candidate_scope == "full_universe"),
  all(full_result$alpha_curve$candidate_count == 6L),
  all(full_result$alpha_curve$dynamic_discoveries == 2L),
  all(full_result$alpha_curve$non_dynamic_calls == 3L),
  all(full_result$alpha_curve$first_stage_null_calls == 0L),
  identical(full_result$candidate_indices, 1:6)
)

screened_result <- evaluate_fash_functional_testing(
  fit = list(),
  true_functionals = true_functionals,
  evaluation_grid = seq(0, 15, by = 0.1),
  alpha_grid = c(0.05, 0.10),
  true_dynamic = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE),
  candidate_scope = "dynamic_fdr_screen",
  middle_window = c(3, 12),
  middle_boundary = "open"
)
stopifnot(
  identical(observed_functional_indices, c(2L, 4L)),
  all(screened_result$alpha_curve$candidate_scope == "dynamic_fdr_screen"),
  all(screened_result$alpha_curve$candidate_count == 2L)
)

full_rows <- full_result$alpha_curve
full_rows$seed <- 12345L
full_summary <- summarize_mc_functional_alpha_curves(full_rows)
stopifnot(
  all(full_summary$candidate_scope == "full_universe"),
  all(full_summary$mean_candidate_count == 6)
)

legacy_screened_rows <- screened_result$alpha_curve
legacy_screened_rows$seed <- 12345L
legacy_screened_rows$candidate_scope <- NULL
legacy_screened_rows$candidate_count <- NULL
legacy_screened_rows$non_dynamic_calls <- NULL
legacy_summary <- summarize_mc_functional_alpha_curves(legacy_screened_rows)
stopifnot(
  all(legacy_summary$candidate_scope == "dynamic_fdr_screen"),
  all(
    legacy_summary$mean_candidate_count ==
      legacy_summary$mean_dynamic_discoveries
  )
)

required_production_files <- c(
  snapshot_functions_path,
  driver_path,
  core_path,
  wrapper_path,
  sbatch_path
)
stopifnot(all(file.exists(required_production_files)))
stopifnot(identical(
  readLines(shared_functions_path, warn = FALSE),
  readLines(snapshot_functions_path, warn = FALSE)
))

driver_text <- paste(readLines(driver_path, warn = FALSE), collapse = "\n")
core_text <- paste(readLines(core_path, warn = FALSE), collapse = "\n")
wrapper_text <- paste(readLines(wrapper_path, warn = FALSE), collapse = "\n")
sbatch_text <- paste(readLines(sbatch_path, warn = FALSE), collapse = "\n")
expected_result_id <- paste0(
  "r3_real_genotype_one_per_gene_J6362_",
  "matched_functional_open_middle_3_12_center_aligned_equal_cells_",
  "relative_location_clearance_full_universe_",
  "paired_posterior_fashr0143_pilot5"
)

stopifnot(
  grepl('candidate_scope = "full_universe"', driver_text, fixed = TRUE),
  grepl("functional_candidate_scope", driver_text, fixed = TRUE),
  grepl("candidate_count == J", driver_text, fixed = TRUE),
  grepl("common_random_seed_raw_bf", core_text, fixed = TRUE),
  grepl("full_universe", core_text, fixed = TRUE),
  grepl(
    "ca3f786ab11749b12f17e9799006a49c4264c623c3131dd18153f33daeb9da18",
    core_text,
    fixed = TRUE
  ),
  grepl(
    "45267b0884168e5ae33cc4f14e3f05b711d961b65bf9c6fbd880e748de064a6e",
    core_text,
    fixed = TRUE
  ),
  grepl("full_universe_", wrapper_text, fixed = TRUE),
  grepl("paired_posterior_fashr0143_pilot5", wrapper_text, fixed = TRUE),
  grepl(expected_result_id, sbatch_text, fixed = TRUE),
  grepl("22_r3_full_universe_functional_core.R", wrapper_text, fixed = TRUE),
  grepl("22_r3_full_universe_functional_core.R", sbatch_text, fixed = TRUE),
  grepl("r3_full_universe_functional_fashr0143", sbatch_text, fixed = TRUE),
  grepl(
    'mechanism_curve$method == "FASH-IWP1-BF"',
    driver_text,
    fixed = TRUE
  ),
  grepl("interval_methods = \"FASH-IWP1-BF\"", driver_text, fixed = TRUE)
)

message("Full-universe R3 functional calibration tests passed.")
