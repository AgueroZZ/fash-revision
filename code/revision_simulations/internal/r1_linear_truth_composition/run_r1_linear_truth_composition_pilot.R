#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr project root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[[1L]] == length(args)) {
    return(default)
  }
  args[[hit[[1L]] + 1L]]
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

capture_warnings <- function(stage, expression) {
  warning_messages <- character()
  value <- withCallingHandlers(
    force(expression),
    warning = function(condition) {
      warning_messages <<- c(warning_messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(
    value = value,
    warnings = unique(data.frame(
      stage = rep(stage, length(warning_messages)),
      message = warning_messages,
      stringsAsFactors = FALSE
    ))
  )
}

summarize_alpha005 <- function(result_table, source) {
  groups <- split(
    result_table,
    list(result_table$scenario, result_table$method),
    drop = TRUE
  )
  rows <- lapply(groups, function(x) {
    discoveries <- sum(x$selected)
    false_discoveries <- sum(x$selected & x$true_null)
    true_alternatives <- sum(!x$true_null)
    true_positives <- sum(x$selected & !x$true_null)
    data.frame(
      scenario = x$scenario[[1L]],
      method = x$method[[1L]],
      alpha = x$alpha[[1L]],
      n_units = nrow(x),
      n_true_alternatives = true_alternatives,
      n_discoveries = discoveries,
      false_discoveries = false_discoveries,
      true_positives = true_positives,
      power = true_positives / true_alternatives,
      realized_fdp = if (discoveries == 0L) 0 else
        false_discoveries / discoveries,
      source = source,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$method_rank <- rank_revision_methods(out$method)
  out <- out[order(out$scenario, out$method_rank, out$method), , drop = FALSE]
  out$method_rank <- NULL
  rownames(out) <- NULL
  out
}

summarize_class_specific_alpha005 <- function(result_table, source) {
  groups <- split(
    result_table,
    list(
      result_table$scenario,
      result_table$method,
      result_table$effect_class
    ),
    drop = TRUE
  )
  rows <- lapply(groups, function(x) {
    is_dynamic_alternative <- !x$true_null[[1L]]
    data.frame(
      scenario = x$scenario[[1L]],
      method = x$method[[1L]],
      effect_class = x$effect_class[[1L]],
      alpha = x$alpha[[1L]],
      n_units = nrow(x),
      n_selected = sum(x$selected),
      selection_rate = mean(x$selected),
      is_dynamic_alternative = is_dynamic_alternative,
      metric = if (is_dynamic_alternative) "class-specific power" else
        "null selection rate",
      source = source,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$method_rank <- rank_revision_methods(out$method)
  out <- out[order(
    out$scenario,
    out$method_rank,
    out$method,
    out$effect_class
  ), , drop = FALSE]
  out$method_rank <- NULL
  rownames(out) <- NULL
  out
}

make_power_plot <- function(alpha005_summary,
                            class_specific_summary,
                            path) {
  method_order <- c("FASH-IWP1-BF", "FASH-linear-BF")
  scenario_order <- c(
    "r1_original_all_bspline_dynamic_truth",
    "r1_linear90_bspline10_dynamic_truth",
    "r1_all_linear_dynamic_truth"
  )
  scenario_labels <- c(
    "Original R1\n100% B-spline",
    "90% linear\n10% B-spline",
    "100% linear"
  )
  method_colors <- c("#0072B2", "#CC79A7")

  overall <- alpha005_summary[
    alpha005_summary$scenario %in% scenario_order &
      alpha005_summary$method %in% method_order,
    ,
    drop = FALSE
  ]
  overall$scenario <- factor(overall$scenario, levels = scenario_order)
  overall$method <- factor(overall$method, levels = method_order)
  overall <- overall[order(overall$method, overall$scenario), , drop = FALSE]
  overall_matrix <- matrix(
    overall$power,
    nrow = length(method_order),
    byrow = TRUE,
    dimnames = list(method_order, scenario_labels)
  )

  mixed <- class_specific_summary[
    class_specific_summary$scenario ==
      "r1_linear90_bspline10_dynamic_truth" &
      class_specific_summary$method %in% method_order &
      class_specific_summary$effect_class %in%
        c("dynamic_linear", "dynamic_bspline"),
    ,
    drop = FALSE
  ]
  mixed$method <- factor(mixed$method, levels = method_order)
  mixed$effect_class <- factor(
    mixed$effect_class,
    levels = c("dynamic_linear", "dynamic_bspline")
  )
  mixed <- mixed[order(mixed$method, mixed$effect_class), , drop = FALSE]
  mixed_matrix <- matrix(
    mixed$selection_rate,
    nrow = length(method_order),
    byrow = TRUE,
    dimnames = list(method_order, c("Linear (n = 1,145)", "B-spline (n = 127)"))
  )

  grDevices::png(path, width = 2200, height = 1050, res = 180)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(7, 4.5, 3.5, 1), xpd = NA)

  bar_positions <- graphics::barplot(
    overall_matrix,
    beside = TRUE,
    col = method_colors,
    ylim = c(0, 1.05),
    ylab = "Overall dynamic power",
    main = "A. Overall power at cumulative-lfdr alpha = 0.05",
    las = 1,
    cex.names = 0.86
  )
  graphics::text(
    x = bar_positions,
    y = overall_matrix + 0.025,
    labels = sprintf("%.3f", overall_matrix),
    cex = 0.78
  )
  graphics::legend(
    "bottom",
    inset = c(0, -0.30),
    legend = method_order,
    fill = method_colors,
    horiz = TRUE,
    bty = "n",
    cex = 0.88
  )

  mixed_positions <- graphics::barplot(
    mixed_matrix,
    beside = TRUE,
    col = method_colors,
    ylim = c(0, 1.05),
    ylab = "Class-specific power",
    main = "B. Power within the 90/10 mixed alternatives",
    las = 1,
    cex.names = 0.86
  )
  graphics::text(
    x = mixed_positions,
    y = mixed_matrix + 0.025,
    labels = sprintf("%.3f", mixed_matrix),
    cex = 0.78
  )
  graphics::legend(
    "bottom",
    inset = c(0, -0.30),
    legend = method_order,
    fill = method_colors,
    horiz = TRUE,
    bty = "n",
    cex = 0.88
  )
  invisible(path)
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal",
  "r1_linear_truth_composition", "r1_linear_truth_composition_helpers.R"
))

seed <- 12345L
J <- 6362L
n_donors <- 19L
n_covariates <- 5L
time_grid <- 0:15
num_cores <- as.integer(get_arg("--num-cores", "2"))
output_id <- get_arg(
  "--output-id",
  "r1_linear_truth_composition_seed12345_pilot"
)
if (length(num_cores) != 1L || is.na(num_cores) ||
    num_cores < 1L || num_cores > 2L || !nzchar(output_id)) {
  stop("num_cores must be one or two and output_id must be nonempty.")
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)
options(mc.cores = num_cores)

output_parent <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal"
)
final_output_dir <- file.path(output_parent, output_id)
staging_output_dir <- file.path(
  output_parent,
  paste0(output_id, ".staging-", Sys.getpid())
)
if (dir.exists(final_output_dir)) {
  stop("The final output directory already exists: ", final_output_dir)
}
if (dir.exists(staging_output_dir)) {
  stop("The staging output directory already exists: ", staging_output_dir)
}
dir.create(staging_output_dir, recursive = TRUE, showWarnings = FALSE)

genotype_cache_path <- file.path(
  workflowr_root,
  "output", "revision_simulations", "shared",
  "real_genotype_one_per_gene_J6362_pilot5", "genotype_samples.rds"
)
formal_output_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "mc",
  paste0(
    "r1_real_genotype_one_per_gene_J6362_random_bspline_main_effect_",
    "linear_mixture_predstep1_penalty10_pilot5"
  )
)
formal_full_fit_path <- file.path(
  formal_output_dir,
  "full_fits", "seed_12345.rds"
)
formal_replicate_curve_path <- file.path(
  formal_output_dir,
  "summary", "iwp_vs_linear_fash_replicate_alpha_curves.csv"
)
required_inputs <- c(
  genotype_cache_path,
  formal_full_fit_path,
  formal_replicate_curve_path
)
if (any(!file.exists(required_inputs))) {
  stop("One or more required formal R1 inputs are missing.")
}

message("Loading and validating the formal R1 seed-12345 reference.")
genotype_cache <- readRDS(genotype_cache_path)
genotype_sample <- validate_real_genotype_sample(
  genotype_cache$samples[[as.character(seed)]],
  expected_genes = J,
  expected_donors = n_donors,
  maf_min = genotype_cache$configuration$maf_min
)
formal_out <- readRDS(formal_full_fit_path)
component_seeds <- revision_component_seeds(seed)
covariates <- simulate_covariate_matrix(
  n_donors = n_donors,
  n_covariates = n_covariates,
  seed = component_seeds[["covariates"]]
)
class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
base_effect_sim <- simulate_variant_effect_curves(
  n_variants = J,
  time_grid = time_grid,
  class_probs = class_probs,
  scenario = "r1_original_all_bspline_dynamic_truth",
  dynamic_amplitude = 2,
  bspline_df = 6,
  bspline_coefficient_sd = 1,
  constant_sd = 1,
  dynamic_main_effect_sd = 1,
  exact_class_counts = TRUE,
  seed = component_seeds[["functional_truth"]]
)
base_effect_sim <- reassign_effect_simulation_by_maf(
  effect_sim = base_effect_sim,
  maf = genotype_sample$variant_info$observed_maf,
  class_probs = class_probs,
  seed = component_seeds[["classes"]],
  n_strata = 10L
)

formal_genotype_exact <- identical(genotype_sample$G, formal_out$genotype)
formal_covariates_exact <- identical(covariates, formal_out$covariates)
formal_classes_exact <- identical(
  base_effect_sim$unit_info$effect_class,
  formal_out$unit_info$effect_class
)
formal_truth_max_difference <- max(abs(
  unname(base_effect_sim$beta_matrix) - unname(formal_out$true_beta)
))
if (!formal_genotype_exact || !formal_covariates_exact ||
    !formal_classes_exact || formal_truth_max_difference > 1e-12) {
  stop("The reconstructed seed-12345 inputs do not match formal R1.")
}

rownames(base_effect_sim$beta_matrix) <- colnames(genotype_sample$G)
base_effect_sim$unit_info$variant_id <- colnames(genotype_sample$G)
truths <- make_r1_linear_truth_scenarios(
  base_effect_sim = base_effect_sim,
  time_grid = time_grid,
  linear_amplitude = 2,
  linear_sign_seed = seed + 30001L,
  mixture_seed = seed + 30002L
)

message("Generating the two paired expression datasets.")
expression_simulations <- lapply(truths$scenarios, function(effect_sim) {
  simulate_eqtl_expression_from_genotypes(
    G = genotype_sample$G,
    beta_matrix = effect_sim$beta_matrix,
    time_grid = time_grid,
    covariates = covariates,
    expression_noise_sd = 1,
    covariate_effect_sd = 0.5,
    intercept_sd = 0,
    seed = component_seeds[["expression"]]
  )
})
paired_intercepts_exact <- identical(
  expression_simulations$all_linear$intercepts,
  expression_simulations$linear90_bspline10$intercepts
)
paired_covariate_effects_exact <- identical(
  expression_simulations$all_linear$covariate_effects,
  expression_simulations$linear90_bspline10$covariate_effects
)
paired_expression_difference_max_error <- 0
for (time_index in seq_along(time_grid)) {
  expected_difference <- sweep(
    genotype_sample$G,
    2,
    truths$scenarios$all_linear$beta_matrix[, time_index] -
      truths$scenarios$linear90_bspline10$beta_matrix[, time_index],
    `*`
  )
  observed_difference <-
    expression_simulations$all_linear$expression[, , time_index] -
    expression_simulations$linear90_bspline10$expression[, , time_index]
  paired_expression_difference_max_error <- max(
    paired_expression_difference_max_error,
    max(abs(observed_difference - expected_difference))
  )
}
if (!paired_intercepts_exact || !paired_covariate_effects_exact ||
    paired_expression_difference_max_error > 1e-12) {
  stop("The two new scenarios do not share the same nuisance realization.")
}

common_sd_grid <- default_revision_grid()
common_pred_step <- 1
common_penalty <- 10L
method_names <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "FASH-linear-Raw",
  "FASH-linear-BF"
)

fit_scenario <- function(scenario_key) {
  effect_sim <- truths$scenarios[[scenario_key]]
  expression_sim <- expression_simulations[[scenario_key]]
  runtime_rows <- list()
  warning_rows <- list()

  message("Estimating per-time effects for scenario ", scenario_key, ".")
  sufficient_time <- system.time({
    eqtl_summary <- estimate_eqtl_summaries_from_genotypes(
      G = genotype_sample$G,
      expression = expression_sim$expression,
      covariates = covariates,
      apply_t_se_correction = TRUE
    )
    datasets <- make_fash_datasets_from_eqtl_summary(
      beta_hat = eqtl_summary$beta_hat,
      se = eqtl_summary$se,
      true_beta = effect_sim$beta_matrix,
      time_grid = time_grid,
      unit_info = effect_sim$unit_info,
      scenario = effect_sim$unit_info$scenario[[1L]]
    )
  })
  unit_info <- attr(datasets, "unit_info")
  runtime_rows[["sufficient_statistics"]] <- data.frame(
    scenario = unit_info$scenario[[1L]],
    stage = "sufficient statistics",
    elapsed_seconds = unname(sufficient_time[["elapsed"]]),
    stringsAsFactors = FALSE
  )

  message("Fitting FASH-IWP1 raw for scenario ", scenario_key, ".")
  iwp_raw_time <- system.time({
    captured <- capture_warnings(
      paste(scenario_key, "FASH-IWP1-Raw"),
      fashr::fash(
        Y = "y",
        smooth_var = "x",
        S = "sd",
        data_list = datasets,
        order = 1,
        grid = common_sd_grid,
        num_basis = 20,
        pred_step = common_pred_step,
        penalty = common_penalty,
        num_cores = num_cores,
        verbose = FALSE
      )
    )
  })
  iwp_raw <- captured$value
  warning_rows[["iwp_raw"]] <- captured$warnings
  runtime_rows[["iwp_raw"]] <- data.frame(
    scenario = unit_info$scenario[[1L]],
    stage = "FASH-IWP1 raw",
    elapsed_seconds = unname(iwp_raw_time[["elapsed"]]),
    stringsAsFactors = FALSE
  )

  message("Applying the FASH-IWP1 BF update for scenario ", scenario_key, ".")
  iwp_bf_time <- system.time({
    captured <- capture_warnings(
      paste(scenario_key, "FASH-IWP1-BF"),
      fashr::BF_update(iwp_raw, plot = FALSE)
    )
  })
  iwp_bf <- captured$value
  warning_rows[["iwp_bf"]] <- captured$warnings
  runtime_rows[["iwp_bf"]] <- data.frame(
    scenario = unit_info$scenario[[1L]],
    stage = "FASH-IWP1 BF update",
    elapsed_seconds = unname(iwp_bf_time[["elapsed"]]),
    stringsAsFactors = FALSE
  )

  message("Fitting FASH-linear raw for scenario ", scenario_key, ".")
  linear_raw_time <- system.time({
    captured <- capture_warnings(
      paste(scenario_key, "FASH-linear-Raw"),
      fit_linear_mixture_fash(
        datasets = datasets,
        grid = common_sd_grid,
        pred_step = common_pred_step,
        penalty = common_penalty
      )
    )
  })
  linear_raw <- captured$value
  warning_rows[["linear_raw"]] <- captured$warnings
  runtime_rows[["linear_raw"]] <- data.frame(
    scenario = unit_info$scenario[[1L]],
    stage = "FASH-linear raw",
    elapsed_seconds = unname(linear_raw_time[["elapsed"]]),
    stringsAsFactors = FALSE
  )

  message("Applying the FASH-linear BF update for scenario ", scenario_key, ".")
  linear_bf_time <- system.time({
    captured <- capture_warnings(
      paste(scenario_key, "FASH-linear-BF"),
      BF_update_linear_mixture_fash(linear_raw)
    )
  })
  linear_bf <- captured$value
  warning_rows[["linear_bf"]] <- captured$warnings
  runtime_rows[["linear_bf"]] <- data.frame(
    scenario = unit_info$scenario[[1L]],
    stage = "FASH-linear BF update",
    elapsed_seconds = unname(linear_bf_time[["elapsed"]]),
    stringsAsFactors = FALSE
  )

  result_table <- rbind(
    evaluate_lfdr_method(
      get_fash_lfdr(iwp_raw), unit_info,
      method = "FASH-IWP1-Raw", target = "dynamic", alpha = 0.05
    ),
    evaluate_lfdr_method(
      get_fash_lfdr(iwp_bf), unit_info,
      method = "FASH-IWP1-BF", target = "dynamic", alpha = 0.05
    ),
    evaluate_lfdr_method(
      linear_raw$lfdr, unit_info,
      method = "FASH-linear-Raw", target = "dynamic", alpha = 0.05
    ),
    evaluate_lfdr_method(
      linear_bf$lfdr, unit_info,
      method = "FASH-linear-BF", target = "dynamic", alpha = 0.05
    )
  )
  alpha_curve <- compute_alpha_curve(
    result_table,
    alpha_grid = seq(0, 0.20, by = 0.005)
  )
  names(alpha_curve)[names(alpha_curve) == "empirical_fdr"] <- "realized_fdp"
  pi0 <- data.frame(
    scenario = unit_info$scenario[[1L]],
    method = c("FASH-IWP1", "FASH-IWP1", "FASH-linear", "FASH-linear"),
    fit = c("Raw", "BF-corrected", "Raw", "BF-corrected"),
    estimated_pi0 = c(
      constant_component_prior_weight(iwp_raw),
      constant_component_prior_weight(iwp_bf),
      constant_component_prior_weight(linear_raw),
      constant_component_prior_weight(linear_bf)
    ),
    stringsAsFactors = FALSE
  )
  prior_weights <- rbind(
    transform(
      extract_linear_mixture_prior_table(
        linear_raw, seed = seed, fit_label = "Raw"
      ),
      scenario = unit_info$scenario[[1L]]
    ),
    transform(
      extract_linear_mixture_prior_table(
        linear_bf, seed = seed, fit_label = "BF-corrected"
      ),
      scenario = unit_info$scenario[[1L]]
    )
  )

  list(
    unit_info = unit_info,
    beta_hat = eqtl_summary$beta_hat,
    se = eqtl_summary$se,
    result_table = result_table,
    alpha005 = summarize_alpha005(result_table, source = "new paired pilot"),
    class_specific_alpha005 = summarize_class_specific_alpha005(
      result_table,
      source = "new paired pilot"
    ),
    alpha_curve = alpha_curve,
    pi0 = pi0,
    prior_weights = prior_weights,
    runtime = do.call(rbind, runtime_rows),
    warnings = do.call(rbind, warning_rows),
    lfdr = list(
      FASH_IWP1_Raw = get_fash_lfdr(iwp_raw),
      FASH_IWP1_BF = get_fash_lfdr(iwp_bf),
      FASH_linear_Raw = linear_raw$lfdr,
      FASH_linear_BF = linear_bf$lfdr
    )
  )
}

pilot_start <- proc.time()[["elapsed"]]
scenario_results <- lapply(names(truths$scenarios), fit_scenario)
names(scenario_results) <- names(truths$scenarios)
pilot_elapsed_seconds <- proc.time()[["elapsed"]] - pilot_start

formal_result_table <- formal_out$result_table[
  formal_out$result_table$method %in% method_names,
  ,
  drop = FALSE
]
formal_result_table$scenario <- "r1_original_all_bspline_dynamic_truth"
formal_alpha005 <- summarize_alpha005(
  formal_result_table,
  source = "existing formal R1 seed-12345 cache"
)
formal_class_specific <- summarize_class_specific_alpha005(
  formal_result_table,
  source = "existing formal R1 seed-12345 cache"
)
formal_alpha_curve <- compute_alpha_curve(
  formal_result_table,
  alpha_grid = seq(0, 0.20, by = 0.005)
)
names(formal_alpha_curve)[names(formal_alpha_curve) == "empirical_fdr"] <-
  "realized_fdp"
formal_pi0 <- data.frame(
  scenario = "r1_original_all_bspline_dynamic_truth",
  method = c("FASH-IWP1", "FASH-IWP1", "FASH-linear", "FASH-linear"),
  fit = c("Raw", "BF-corrected", "Raw", "BF-corrected"),
  estimated_pi0 = c(
    constant_component_prior_weight(formal_out$fash_fits$fash_iwp1_raw),
    constant_component_prior_weight(formal_out$fash_fits$fash_iwp1_bf),
    constant_component_prior_weight(formal_out$simplified_fit),
    constant_component_prior_weight(formal_out$simplified_fit_bf)
  ),
  stringsAsFactors = FALSE
)

new_alpha005 <- do.call(rbind, lapply(scenario_results, `[[`, "alpha005"))
new_class_specific <- do.call(rbind, lapply(
  scenario_results,
  `[[`,
  "class_specific_alpha005"
))
new_alpha_curves <- do.call(rbind, lapply(
  scenario_results,
  `[[`,
  "alpha_curve"
))
new_pi0 <- do.call(rbind, lapply(scenario_results, `[[`, "pi0"))
new_prior_weights <- do.call(rbind, lapply(
  scenario_results,
  `[[`,
  "prior_weights"
))
runtime <- do.call(rbind, lapply(scenario_results, `[[`, "runtime"))
warnings <- do.call(rbind, lapply(scenario_results, `[[`, "warnings"))

alpha005_summary <- rbind(formal_alpha005, new_alpha005)
class_specific_summary <- rbind(formal_class_specific, new_class_specific)
alpha_curves <- rbind(formal_alpha_curve, new_alpha_curves)
pi0_summary <- rbind(formal_pi0, new_pi0)
runtime <- rbind(
  runtime,
  data.frame(
    scenario = "both new scenarios",
    stage = "total production pilot",
    elapsed_seconds = pilot_elapsed_seconds,
    stringsAsFactors = FALSE
  )
)

formal_reference_curve <- read.csv(
  formal_replicate_curve_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
formal_reference_alpha005 <- formal_reference_curve[
  formal_reference_curve$seed == seed &
    abs(formal_reference_curve$alpha - 0.05) < 1e-12 &
    formal_reference_curve$method %in% method_names,
  ,
  drop = FALSE
]
formal_reference_alpha005 <- formal_reference_alpha005[
  match(method_names, formal_reference_alpha005$method),
  ,
  drop = FALSE
]
formal_alpha005_ordered <- formal_alpha005[
  match(method_names, formal_alpha005$method),
  ,
  drop = FALSE
]
formal_summary_max_difference <- max(abs(cbind(
  formal_reference_alpha005$n_discoveries -
    formal_alpha005_ordered$n_discoveries,
  formal_reference_alpha005$false_discoveries -
    formal_alpha005_ordered$false_discoveries,
  formal_reference_alpha005$true_positives -
    formal_alpha005_ordered$true_positives,
  formal_reference_alpha005$power - formal_alpha005_ordered$power,
  formal_reference_alpha005$empirical_fdr -
    formal_alpha005_ordered$realized_fdp
)))

scenario_counts <- do.call(rbind, lapply(
  truths$scenarios,
  function(effect_sim) {
    counts <- table(effect_sim$unit_info$effect_class)
    data.frame(
      scenario = effect_sim$unit_info$scenario[[1L]],
      effect_class = names(counts),
      n = as.integer(counts),
      stringsAsFactors = FALSE
    )
  }
))
expected_scenario_counts <- data.frame(
  scenario = c(
    "r1_all_linear_dynamic_truth",
    "r1_all_linear_dynamic_truth",
    "r1_all_linear_dynamic_truth",
    "r1_linear90_bspline10_dynamic_truth",
    "r1_linear90_bspline10_dynamic_truth",
    "r1_linear90_bspline10_dynamic_truth",
    "r1_linear90_bspline10_dynamic_truth"
  ),
  effect_class = c(
    "constant", "dynamic_linear", "zero",
    "constant", "dynamic_bspline", "dynamic_linear", "zero"
  ),
  n = c(2545L, 1272L, 2545L, 2545L, 127L, 1145L, 2545L),
  stringsAsFactors = FALSE
)
scenario_counts <- scenario_counts[order(
  scenario_counts$scenario,
  scenario_counts$effect_class
), , drop = FALSE]
expected_scenario_counts <- expected_scenario_counts[order(
  expected_scenario_counts$scenario,
  expected_scenario_counts$effect_class
), , drop = FALSE]
rownames(scenario_counts) <- NULL
rownames(expected_scenario_counts) <- NULL

expected_scenarios <- c(
  "r1_original_all_bspline_dynamic_truth",
  "r1_all_linear_dynamic_truth",
  "r1_linear90_bspline10_dynamic_truth"
)
observed_alpha005_keys <- paste(
  alpha005_summary$scenario,
  alpha005_summary$method,
  sep = "\r"
)
expected_alpha005_keys <- as.vector(outer(
  expected_scenarios,
  method_names,
  paste,
  sep = "\r"
))

validation <- data.frame(
  check = c(
    "formal genotype exact",
    "formal covariates exact",
    "formal effect classes exact",
    "formal true beta exact",
    "formal alpha-0.05 cache reproduction",
    "scenario class counts",
    "paired intercepts exact",
    "paired covariate effects exact",
    "paired expression difference explained by truth",
    "complete alpha-0.05 scenario-method keys",
    "finite alpha-0.05 metrics",
    "two-thread cap enforced"
  ),
  passed = c(
    formal_genotype_exact,
    formal_covariates_exact,
    formal_classes_exact,
    formal_truth_max_difference <= 1e-12,
    formal_summary_max_difference <= 1e-12,
    identical(scenario_counts, expected_scenario_counts),
    paired_intercepts_exact,
    paired_covariate_effects_exact,
    paired_expression_difference_max_error <= 1e-12,
    length(observed_alpha005_keys) == length(expected_alpha005_keys) &&
      !anyDuplicated(observed_alpha005_keys) &&
      setequal(observed_alpha005_keys, expected_alpha005_keys),
    all(is.finite(alpha005_summary$power)) &&
      all(is.finite(alpha005_summary$realized_fdp)) &&
      all(alpha005_summary$power >= 0 & alpha005_summary$power <= 1) &&
      all(alpha005_summary$realized_fdp >= 0 &
        alpha005_summary$realized_fdp <= 1),
    num_cores <= 2L &&
      Sys.getenv("OMP_NUM_THREADS") == "1" &&
      Sys.getenv("VECLIB_MAXIMUM_THREADS") == "1"
  ),
  observed = c(
    as.character(formal_genotype_exact),
    as.character(formal_covariates_exact),
    as.character(formal_classes_exact),
    format(formal_truth_max_difference, scientific = TRUE),
    format(formal_summary_max_difference, scientific = TRUE),
    paste(paste(scenario_counts$effect_class, scenario_counts$n, sep = "="),
      collapse = "; "),
    as.character(paired_intercepts_exact),
    as.character(paired_covariate_effects_exact),
    format(paired_expression_difference_max_error, scientific = TRUE),
    paste(length(observed_alpha005_keys), "unique keys"),
    paste(range(alpha005_summary$power), collapse = " to "),
    paste(num_cores, "FASH workers; one BLAS thread per process")
  ),
  stringsAsFactors = FALSE
)
if (any(!validation$passed)) {
  print(validation)
  stop("One or more retained-pilot validation checks failed.")
}

configuration <- list(
  output_id = output_id,
  scope = "unlinked internal single-seed paired pilot",
  seed = seed,
  component_seeds = component_seeds,
  linear_sign_seed = seed + 30001L,
  mixture_seed = seed + 30002L,
  J = J,
  n_donors = n_donors,
  n_covariates = n_covariates,
  time_grid = time_grid,
  class_probs = class_probs,
  scenario_counts = scenario_counts,
  linear_definition = list(
    centered_time_varying_deviation = TRUE,
    maximum_absolute_deviation = 2,
    direction = "reproducible random sign",
    genetic_main_effect = "exact formal R1 N(0,1) realization"
  ),
  observation_model = list(
    expression_noise_sd = 1,
    covariate_effect_sd = 0.5,
    intercept_sd = 0,
    apply_t_se_correction = TRUE
  ),
  method_settings = list(
    grid = common_sd_grid,
    pred_step = common_pred_step,
    penalty = common_penalty,
    num_basis = 20L,
    num_cores = num_cores,
    alpha = 0.05,
    discovery_rule = "cumulative lfdr"
  ),
  input_paths = list(
    genotype_cache = normalizePath(genotype_cache_path, winslash = "/"),
    formal_full_fit = normalizePath(formal_full_fit_path, winslash = "/"),
    formal_replicate_curve = normalizePath(
      formal_replicate_curve_path,
      winslash = "/"
    )
  ),
  input_fingerprints = list(
    genotype_cache = artifact_fingerprint(genotype_cache_path),
    formal_full_fit = artifact_fingerprint(formal_full_fit_path),
    formal_replicate_curve = artifact_fingerprint(formal_replicate_curve_path)
  ),
  runtime_environment = list(
    R_version = R.version.string,
    fashr_version = as.character(utils::packageVersion("fashr")),
    OMP_NUM_THREADS = Sys.getenv("OMP_NUM_THREADS"),
    OPENBLAS_NUM_THREADS = Sys.getenv("OPENBLAS_NUM_THREADS"),
    VECLIB_MAXIMUM_THREADS = Sys.getenv("VECLIB_MAXIMUM_THREADS")
  )
)

unit_results <- rbind(
  formal_result_table,
  do.call(rbind, lapply(scenario_results, `[[`, "result_table"))
)
analysis_cache <- list(
  configuration = configuration,
  truth_membership = truths$membership,
  scenario_counts = scenario_counts,
  alpha005_summary = alpha005_summary,
  class_specific_alpha005 = class_specific_summary,
  alpha_curves = alpha_curves,
  pi0_summary = pi0_summary,
  linear_prior_weights = new_prior_weights,
  unit_results = unit_results,
  new_scenario_sufficient_statistics = lapply(
    scenario_results,
    function(result) list(
      unit_info = result$unit_info,
      beta_hat = result$beta_hat,
      se = result$se,
      lfdr = result$lfdr
    )
  ),
  runtime = runtime,
  warnings = warnings,
  validation = validation
)

saveRDS(configuration, file.path(staging_output_dir, "configuration.rds"))
saveRDS(analysis_cache, file.path(staging_output_dir, "analysis_cache.rds"))
write_csv(truths$membership, file.path(staging_output_dir, "truth_membership.csv"))
write_csv(scenario_counts, file.path(staging_output_dir, "scenario_counts.csv"))
write_csv(alpha005_summary, file.path(staging_output_dir, "alpha005_summary.csv"))
write_csv(
  class_specific_summary,
  file.path(staging_output_dir, "class_specific_alpha005.csv")
)
write_csv(alpha_curves, file.path(staging_output_dir, "alpha_curves.csv"))
write_csv(pi0_summary, file.path(staging_output_dir, "pi0_summary.csv"))
write_csv(
  new_prior_weights,
  file.path(staging_output_dir, "linear_prior_weights.csv")
)
write_csv(runtime, file.path(staging_output_dir, "runtime.csv"))
write_csv(warnings, file.path(staging_output_dir, "warnings.csv"))
write_csv(validation, file.path(staging_output_dir, "validation.csv"))
make_power_plot(
  alpha005_summary = alpha005_summary,
  class_specific_summary = class_specific_summary,
  path = file.path(staging_output_dir, "power_comparison.png")
)

required_outputs <- c(
  "configuration.rds",
  "analysis_cache.rds",
  "truth_membership.csv",
  "scenario_counts.csv",
  "alpha005_summary.csv",
  "class_specific_alpha005.csv",
  "alpha_curves.csv",
  "pi0_summary.csv",
  "linear_prior_weights.csv",
  "runtime.csv",
  "warnings.csv",
  "validation.csv",
  "power_comparison.png"
)
required_output_paths <- file.path(staging_output_dir, required_outputs)
if (any(!file.exists(required_output_paths)) ||
    any(file.info(required_output_paths)$size <= 0)) {
  stop("One or more required pilot outputs are missing or empty.")
}
if (!file.rename(staging_output_dir, final_output_dir)) {
  stop("Could not atomically finalize the internal pilot output directory.")
}

cat("\nR1 linear-truth composition pilot completed.\n")
cat("Output: ", normalizePath(final_output_dir, winslash = "/"), "\n", sep = "")
cat("Elapsed seconds: ", sprintf("%.1f", pilot_elapsed_seconds), "\n", sep = "")
print(alpha005_summary[
  alpha005_summary$method %in% c("FASH-IWP1-BF", "FASH-linear-BF"),
  c("scenario", "method", "power", "realized_fdp", "n_discoveries")
])

