#!/usr/bin/env Rscript

# Screen paired sharp transient effect settings using power and discovery counts.
# Realized FDP is saved for later diagnostics but is not used to select a setting.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find code/revision_simulations/shared/simulation_functions.R.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0 || hit[1] == length(args)) {
    return(default)
  }
  args[hit[1] + 1]
}

as_flag <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

safe_numeric_label <- function(x) {
  gsub("\\.", "p", format(x, trim = TRUE, scientific = FALSE))
}

candidate_table <- function() {
  data.frame(
    setting_id = c(
      "legacy_df16_amp2",
      "exact_df16_amp2",
      "sharp_df16_amp2p5",
      "sharp_df16_amp2p75",
      "sharp_df16_amp3",
      "local_df12_amp2p5",
      "local_df12_amp3"
    ),
    transient_bspline_df = c(16L, 16L, 16L, 16L, 16L, 12L, 12L),
    dynamic_amplitude = c(2, 2, 2.5, 2.75, 3, 2.5, 3),
    normalization = c(
      "legacy",
      rep("center_then_scale", 6)
    ),
    selection_eligible = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
}

select_candidate_without_fdp <- function(screen_table,
                                         target_discoveries = c(80, 150),
                                         target_power = c(0.40, 0.70),
                                         reference_rms = 0.94) {
  candidates <- screen_table[
    screen_table$method == "FASH-IWP1-BF" &
      screen_table$selection_eligible,
    ,
    drop = FALSE
  ]
  eligible <- candidates$n_discoveries >= target_discoveries[1] &
    candidates$n_discoveries <= target_discoveries[2] &
    candidates$power >= target_power[1] &
    candidates$power <= target_power[2]

  if (any(eligible)) {
    candidates <- candidates[eligible, , drop = FALSE]
    selection_basis <- "within prespecified discovery and power ranges"
  } else {
    discovery_distance <- pmax(
      target_discoveries[1] - candidates$n_discoveries,
      0,
      candidates$n_discoveries - target_discoveries[2]
    ) / diff(target_discoveries)
    power_distance <- pmax(
      target_power[1] - candidates$power,
      0,
      candidates$power - target_power[2]
    ) / diff(target_power)
    candidates$range_distance <- discovery_distance + power_distance
    candidates <- candidates[
      candidates$range_distance == min(candidates$range_distance),
      ,
      drop = FALSE
    ]
    selection_basis <- "closest to prespecified discovery and power ranges"
  }

  candidates$rms_distance <- abs(candidates$rms_median - reference_rms)
  candidates$discovery_target_distance <- abs(candidates$n_discoveries - 100)
  candidates <- candidates[order(
    -candidates$transient_bspline_df,
    candidates$rms_distance,
    candidates$discovery_target_distance,
    candidates$setting_id
  ), , drop = FALSE]

  selected <- candidates[1, , drop = FALSE]
  selected$selection_basis <- selection_basis
  selected$target_discoveries <- paste(target_discoveries, collapse = "-")
  selected$target_power <- paste(target_power, collapse = "-")
  selected$reference_rms <- reference_rms
  selected
}

validate_pairing <- function(outputs, candidate_settings, tolerance = 1e-10) {
  reference_id <- candidate_settings$setting_id[1]
  reference <- outputs[[reference_id]]
  reference_null <- reference$unit_info$effect_class %in% c("constant", "zero")
  checks <- lapply(candidate_settings$setting_id, function(setting_id) {
    current <- outputs[[setting_id]]
    max_null_truth_difference <- max(abs(
      current$true_beta[reference_null, , drop = FALSE] -
        reference$true_beta[reference_null, , drop = FALSE]
    ))
    max_expression_delta_error <- 0
    for (tt in seq_along(reference$settings$time_grid)) {
      expected_delta <- sweep(
        reference$genotype,
        2,
        current$true_beta[, tt] - reference$true_beta[, tt],
        `*`
      )
      observed_delta <- current$expression[, , tt] -
        reference$expression[, , tt]
      max_expression_delta_error <- max(
        max_expression_delta_error,
        max(abs(observed_delta - expected_delta))
      )
    }
    data.frame(
      setting_id = setting_id,
      identical_genotype = identical(current$genotype, reference$genotype),
      identical_covariates = identical(current$covariates, reference$covariates),
      identical_effect_classes = identical(
        current$unit_info$effect_class,
        reference$unit_info$effect_class
      ),
      max_null_truth_difference = max_null_truth_difference,
      max_expression_delta_error = max_expression_delta_error,
      paired = max_null_truth_difference <= tolerance &&
        max_expression_delta_error <= tolerance,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, checks)
}

workflowr_root <- find_workflowr_root()
source(file.path(workflowr_root, "code", "revision_simulations", "shared", "simulation_functions.R"))

n_donors <- as.integer(get_arg("--n-donors", "19"))
n_variants <- as.integer(get_arg("--J", "1000"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
seed <- as.integer(get_arg("--seed", "12345"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
num_basis <- as.integer(get_arg("--num-basis", "20"))
expression_noise_sd <- as.numeric(get_arg("--noise-sd", "1"))
efdr_permutations <- as.integer(get_arg("--efdr-permutations", "0"))
overwrite <- as_flag(get_arg("--overwrite", "false"))
candidate_ids_arg <- get_arg("--candidate-ids", "all")

time_grid <- make_time_grid()
class_probs <- c(
  dynamic_local_bspline_transient = 0.20,
  constant = 0.40,
  zero = 0.40
)
candidates <- candidate_table()
candidate_scope <- "all"
if (!identical(candidate_ids_arg, "all")) {
  requested_candidate_ids <- strsplit(candidate_ids_arg, ",", fixed = TRUE)[[1]]
  unknown_candidate_ids <- setdiff(requested_candidate_ids, candidates$setting_id)
  if (length(unknown_candidate_ids) > 0) {
    stop(
      "Unknown candidate IDs: ",
      paste(unknown_candidate_ids, collapse = ", ")
    )
  }
  candidates <- candidates[
    match(requested_candidate_ids, candidates$setting_id),
    ,
    drop = FALSE
  ]
  candidate_scope <- paste(requested_candidate_ids, collapse = "-")
}
component_seeds <- revision_component_seeds(seed)

output_root <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "pilot",
  "spiky_signal_strength"
)
raw_dir <- file.path(output_root, "raw")
summary_dir <- file.path(output_root, "summary")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

pilot_stem <- paste0(
  "spiky_signal_pilot_v1",
  "_N", n_donors,
  "_J", n_variants,
  "_T", length(time_grid),
  "_pc", n_covariates,
  "_noise", safe_numeric_label(expression_noise_sd),
  "_seed", seed,
  if (candidate_scope == "all") "" else paste0("_subset_", candidate_scope)
)

genotype_sim <- simulate_genotype_matrix(
  n_donors = n_donors,
  n_variants = n_variants,
  seed = component_seeds[["genotype"]]
)
covariates <- simulate_covariate_matrix(
  n_donors = n_donors,
  n_covariates = n_covariates,
  seed = component_seeds[["covariates"]]
)
paired_effects <- simulate_paired_local_bspline_effect_sets(
  candidate_settings = candidates,
  n_variants = n_variants,
  time_grid = time_grid,
  class_probs = class_probs,
  seed = seed
)

outputs <- vector("list", nrow(candidates))
names(outputs) <- candidates$setting_id
geometry_rows <- vector("list", nrow(candidates))
method_rows <- vector("list", nrow(candidates))

for (candidate_row in seq_len(nrow(candidates))) {
  setting <- candidates[candidate_row, , drop = FALSE]
  setting_id <- setting$setting_id
  raw_path <- file.path(raw_dir, paste0(pilot_stem, "_", setting_id, ".rds"))

  if (file.exists(raw_path) && !overwrite) {
    message("Reusing pilot candidate: ", setting_id)
    out <- readRDS(raw_path)
  } else {
    message("Running pilot candidate: ", setting_id)
    effect_sim <- paired_effects$effect_sets[[setting_id]]
    expression_sim <- simulate_eqtl_expression_from_genotypes(
      G = genotype_sim$G,
      beta_matrix = effect_sim$beta_matrix,
      time_grid = time_grid,
      covariates = covariates,
      expression_noise_sd = expression_noise_sd,
      seed = component_seeds[["expression"]]
    )
    out <- run_genotype_level_dynamic_eqtl_simulation(
      G = genotype_sim$G,
      time_grid = time_grid,
      covariates = covariates,
      class_probs = class_probs,
      expression_noise_sd = expression_noise_sd,
      dynamic_amplitude = setting$dynamic_amplitude,
      transient_bspline_df = setting$transient_bspline_df,
      transient_bspline_degree = 3,
      alpha = 0.05,
      seed = NULL,
      sigma_beta = 1,
      estimate_sigma = FALSE,
      num_cores = num_cores,
      num_basis = num_basis,
      apply_t_se_correction = TRUE,
      scenario = effect_sim$unit_info$scenario[1],
      save_outputs = FALSE,
      verbose = FALSE,
      effect_sim = effect_sim,
      expression_sim = expression_sim
    )
    out$pilot <- list(
      setting = as.list(setting),
      component_seeds = component_seeds,
      selection_uses_realized_fdp = FALSE
    )
    saveRDS(out, raw_path)
  }

  outputs[[setting_id]] <- out
  geometry <- summarize_dynamic_effect_geometry(list(
    beta_matrix = out$true_beta,
    unit_info = out$unit_info
  ))
  geometry_rows[[candidate_row]] <- cbind(setting, geometry)
  method_summary <- out$summary_table[
    out$summary_table$method %in% c(
      "FASH-IWP1-Raw",
      "FASH-IWP1-BF",
      "FASH-linear-Raw",
      "FASH-linear-BF"
    ),
    ,
    drop = FALSE
  ]
  method_rows[[candidate_row]] <- merge(
    cbind(setting, geometry),
    method_summary,
    by = NULL
  )
}

geometry_table <- do.call(rbind, geometry_rows)
method_table <- do.call(rbind, method_rows)
pairing_table <- validate_pairing(outputs, candidates)
if (!all(pairing_table$paired)) {
  stop("At least one candidate failed the paired-data validation.")
}

selection <- select_candidate_without_fdp(method_table)
selected_id <- selection$setting_id

write.csv(
  candidates,
  file.path(summary_dir, paste0(pilot_stem, "_candidate_settings.csv")),
  row.names = FALSE
)
write.csv(
  geometry_table,
  file.path(summary_dir, paste0(pilot_stem, "_effect_geometry.csv")),
  row.names = FALSE
)
write.csv(
  method_table,
  file.path(summary_dir, paste0(pilot_stem, "_method_summary.csv")),
  row.names = FALSE
)
write.csv(
  pairing_table,
  file.path(summary_dir, paste0(pilot_stem, "_pairing_validation.csv")),
  row.names = FALSE
)
write.csv(
  selection,
  file.path(summary_dir, paste0(pilot_stem, "_selected_candidate.csv")),
  row.names = FALSE
)

display_columns <- c(
  "setting_id",
  "transient_bspline_df",
  "dynamic_amplitude",
  "normalization",
  "rms_median",
  "half_support_median",
  "method",
  "n_discoveries",
  "power"
)
print(method_table[, display_columns], row.names = FALSE)
message("Selected without using realized FDP: ", selected_id)

if (efdr_permutations > 0) {
  selected <- outputs[[selected_id]]
  true_pi0 <- dynamic_null_proportion(selected$unit_info, target = "dynamic")
  efdr_path <- file.path(
    raw_dir,
    paste0(
      pilot_stem,
      "_",
      selected_id,
      "_efdr_B",
      efdr_permutations,
      ".rds"
    )
  )
  if (file.exists(efdr_path) && !overwrite) {
    message("Reusing selected-candidate direct eFDR cache: ", efdr_path)
    selected <- readRDS(efdr_path)
  } else {
    message(
      "Running selected-candidate direct eFDR with B = ",
      efdr_permutations,
      " and true pi0 = ",
      format(true_pi0, digits = 3),
      "."
    )
    selected <- add_direct_interaction_efdr_results_to_genotype_output(
      selected,
      n_permutations = efdr_permutations,
      alpha = 0.05,
      seed = component_seeds[["permutations"]],
      lambda = 0.5,
      pi0_method = "conservative",
      true_pi0 = true_pi0,
      include_true_pi0 = TRUE,
      permute_covariates_with_expression = TRUE,
      overwrite = TRUE,
      verbose = TRUE
    )
    saveRDS(selected, efdr_path)
  }
  selected_methods <- c(
    "FASH-IWP1-Raw",
    "FASH-IWP1-BF",
    "Direct-linear-LRT-eFDR-true-pi0",
    "Direct-quadratic-LRT-eFDR-true-pi0"
  )
  selected_summary <- selected$summary_table[
    selected$summary_table$method %in% selected_methods,
    ,
    drop = FALSE
  ]
  write.csv(
    selected_summary,
    file.path(
      summary_dir,
      paste0(
        pilot_stem,
        "_",
        selected_id,
        "_efdr_B",
        efdr_permutations,
        "_summary.csv"
      )
    ),
    row.names = FALSE
  )
  print(selected_summary, row.names = FALSE)
}
