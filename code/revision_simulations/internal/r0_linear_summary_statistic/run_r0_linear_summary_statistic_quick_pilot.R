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
  if (length(hit) == 0L || hit[[1L]] == length(args)) return(default)
  args[[hit[[1L]] + 1L]]
}

capture_warnings <- function(stage, expression) {
  messages <- character()
  value <- withCallingHandlers(
    force(expression),
    warning = function(condition) {
      messages <<- c(messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(
    value = value,
    warnings = unique(data.frame(
      stage = rep(stage, length(messages)),
      message = messages,
      stringsAsFactors = FALSE
    ))
  )
}

summarize_alpha005 <- function(result_table, se_regime) {
  groups <- split(result_table, result_table$method)
  rows <- lapply(groups, function(x) {
    discoveries <- sum(x$selected)
    false_discoveries <- sum(x$selected & x$true_null)
    true_alternatives <- sum(!x$true_null)
    true_positives <- sum(x$selected & !x$true_null)
    data.frame(
      se_regime = se_regime,
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
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out$method_rank <- rank_revision_methods(out$method)
  out <- out[order(out$method_rank, out$method), , drop = FALSE]
  out$method_rank <- NULL
  rownames(out) <- NULL
  out
}

make_power_fdp_plot <- function(summary_table, path) {
  methods <- c("FASH-IWP1-BF", "FASH-linear-BF")
  regimes <- c("medium", "high", "very_high")
  labels <- c("Medium\n0.3 / 0.5 / 0.8", "High\n0.5 / 0.8 / 1.2", "Very high\n0.8 / 1.2 / 1.6")
  colors <- c("#0072B2", "#CC79A7")
  selected <- summary_table[
    summary_table$method %in% methods & summary_table$se_regime %in% regimes,
    ,
    drop = FALSE
  ]
  selected$method <- factor(selected$method, levels = methods)
  selected$se_regime <- factor(selected$se_regime, levels = regimes)
  selected <- selected[order(selected$method, selected$se_regime), , drop = FALSE]
  power <- matrix(selected$power, nrow = 2L, byrow = TRUE,
    dimnames = list(methods, labels))
  fdp <- matrix(selected$realized_fdp, nrow = 2L, byrow = TRUE,
    dimnames = list(methods, labels))

  grDevices::png(path, width = 2100, height = 1000, res = 180)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(1, 2), mar = c(7, 4.5, 3.5, 1), xpd = FALSE)

  positions <- graphics::barplot(
    power, beside = TRUE, col = colors, ylim = c(0, 1), las = 1,
    ylab = "Power", main = "A. All-linear power at alpha = 0.05",
    cex.names = 0.84
  )
  graphics::text(positions, power + 0.025, sprintf("%.3f", power), cex = 0.78)
  graphics::legend("bottom", inset = c(0, -0.30), legend = methods,
    fill = colors, horiz = TRUE, bty = "n", cex = 0.86, xpd = NA)

  positions <- graphics::barplot(
    fdp, beside = TRUE, col = colors, ylim = c(0, max(0.10, fdp * 1.35)),
    las = 1, ylab = "Realized FDP", main = "B. Single-seed realized FDP",
    cex.names = 0.84
  )
  graphics::abline(h = 0.05, lty = 3, col = "gray35")
  graphics::text(positions, fdp + 0.004, sprintf("%.3f", fdp), cex = 0.78)
  graphics::legend("bottom", inset = c(0, -0.30), legend = methods,
    fill = colors, horiz = TRUE, bty = "n", cex = 0.86, xpd = NA)
  invisible(path)
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))

seed <- 12345L
J <- 1000L
time_grid <- 0:15
class_probs <- c(constant = 0.80, linear = 0.20)
num_cores <- as.integer(get_arg("--num-cores", "2"))
output_id <- "r0_all_linear_summary_statistic_seed12345_quick"
if (is.na(num_cores) || num_cores < 1L || num_cores > 2L) {
  stop("num_cores must be one or two.")
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  BLIS_NUM_THREADS = "1"
)
options(mc.cores = num_cores)

se_specs <- list(
  medium = c(0.3, 0.5, 0.8),
  high = c(0.5, 0.8, 1.2),
  very_high = c(0.8, 1.2, 1.6)
)
datasets_by_regime <- lapply(names(se_specs), function(regime) {
  simulate_revision_datasets(
    J = J,
    class_probs = class_probs,
    time_grid = time_grid,
    sd_values = se_specs[[regime]],
    scenario = paste0("r0_all_linear_", regime, "_se"),
    exact_class_counts = TRUE,
    seed = seed
  )
})
names(datasets_by_regime) <- names(se_specs)

reference <- datasets_by_regime$medium
reference_classes <- attr(reference, "unit_info")$effect_class
paired_truth_max_difference <- 0
paired_standardized_noise_max_difference <- 0
paired_se_category_exact <- TRUE
for (regime in setdiff(names(se_specs), "medium")) {
  current <- datasets_by_regime[[regime]]
  paired_truth_max_difference <- max(
    paired_truth_max_difference,
    max(vapply(seq_len(J), function(index) {
      max(abs(reference[[index]]$truef - current[[index]]$truef))
    }, numeric(1)))
  )
  paired_standardized_noise_max_difference <- max(
    paired_standardized_noise_max_difference,
    max(vapply(seq_len(J), function(index) {
      reference_z <- (reference[[index]]$y - reference[[index]]$truef) /
        reference[[index]]$sd
      current_z <- (current[[index]]$y - current[[index]]$truef) /
        current[[index]]$sd
      max(abs(reference_z - current_z))
    }, numeric(1)))
  )
  reference_categories <- unlist(lapply(reference, function(dat) {
    match(dat$sd, se_specs$medium)
  }), use.names = FALSE)
  current_categories <- unlist(lapply(current, function(dat) {
    match(dat$sd, se_specs[[regime]])
  }), use.names = FALSE)
  paired_se_category_exact <- paired_se_category_exact &&
    identical(reference_categories, current_categories) &&
    identical(reference_classes, attr(current, "unit_info")$effect_class)
}

common_grid <- default_revision_grid()
common_pred_step <- 1
common_penalty <- 10L

fit_one_regime <- function(regime) {
  datasets <- datasets_by_regime[[regime]]
  unit_info <- attr(datasets, "unit_info")
  runtime <- list()
  warnings <- list()

  message("Fitting FASH-IWP1 raw for ", regime, " SE.")
  elapsed <- system.time({
    captured <- capture_warnings(
      paste(regime, "FASH-IWP1-Raw"),
      fashr::fash(
        Y = "y", smooth_var = "x", S = "sd", data_list = datasets,
        order = 1, grid = common_grid, num_basis = 20,
        pred_step = common_pred_step, penalty = common_penalty,
        num_cores = num_cores, verbose = FALSE
      )
    )
  })
  iwp_raw <- captured$value
  warnings$iwp_raw <- captured$warnings
  runtime$iwp_raw <- data.frame(
    se_regime = regime, stage = "FASH-IWP1 raw",
    elapsed_seconds = unname(elapsed[["elapsed"]]), stringsAsFactors = FALSE
  )

  elapsed <- system.time({
    captured <- capture_warnings(
      paste(regime, "FASH-IWP1-BF"),
      fashr::BF_update(iwp_raw, plot = FALSE)
    )
  })
  iwp_bf <- captured$value
  warnings$iwp_bf <- captured$warnings
  runtime$iwp_bf <- data.frame(
    se_regime = regime, stage = "FASH-IWP1 BF update",
    elapsed_seconds = unname(elapsed[["elapsed"]]), stringsAsFactors = FALSE
  )

  message("Fitting FASH-linear raw/BF for ", regime, " SE.")
  elapsed <- system.time({
    captured <- capture_warnings(
      paste(regime, "FASH-linear-Raw"),
      fit_linear_mixture_fash(
        datasets = datasets, grid = common_grid,
        pred_step = common_pred_step, penalty = common_penalty
      )
    )
  })
  linear_raw <- captured$value
  warnings$linear_raw <- captured$warnings
  runtime$linear_raw <- data.frame(
    se_regime = regime, stage = "FASH-linear raw",
    elapsed_seconds = unname(elapsed[["elapsed"]]), stringsAsFactors = FALSE
  )

  elapsed <- system.time({
    captured <- capture_warnings(
      paste(regime, "FASH-linear-BF"),
      BF_update_linear_mixture_fash(linear_raw)
    )
  })
  linear_bf <- captured$value
  warnings$linear_bf <- captured$warnings
  runtime$linear_bf <- data.frame(
    se_regime = regime, stage = "FASH-linear BF update",
    elapsed_seconds = unname(elapsed[["elapsed"]]), stringsAsFactors = FALSE
  )

  result_table <- rbind(
    evaluate_lfdr_method(get_fash_lfdr(iwp_raw), unit_info,
      "FASH-IWP1-Raw", "dynamic", 0.05),
    evaluate_lfdr_method(get_fash_lfdr(iwp_bf), unit_info,
      "FASH-IWP1-BF", "dynamic", 0.05),
    evaluate_lfdr_method(linear_raw$lfdr, unit_info,
      "FASH-linear-Raw", "dynamic", 0.05),
    evaluate_lfdr_method(linear_bf$lfdr, unit_info,
      "FASH-linear-BF", "dynamic", 0.05)
  )
  result_table$se_regime <- regime
  pi0 <- data.frame(
    se_regime = regime,
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

  list(
    result_table = result_table,
    alpha005 = summarize_alpha005(result_table, regime),
    pi0 = pi0,
    runtime = do.call(rbind, runtime),
    warnings = do.call(rbind, warnings)
  )
}

pilot_start <- proc.time()[["elapsed"]]
results <- lapply(names(se_specs), fit_one_regime)
names(results) <- names(se_specs)
pilot_elapsed <- proc.time()[["elapsed"]] - pilot_start

unit_results <- do.call(rbind, lapply(results, `[[`, "result_table"))
alpha005_summary <- do.call(rbind, lapply(results, `[[`, "alpha005"))
pi0_summary <- do.call(rbind, lapply(results, `[[`, "pi0"))
runtime <- do.call(rbind, lapply(results, `[[`, "runtime"))
runtime <- rbind(runtime, data.frame(
  se_regime = "all", stage = "total production pilot",
  elapsed_seconds = pilot_elapsed, stringsAsFactors = FALSE
))
warnings <- do.call(rbind, lapply(results, `[[`, "warnings"))

expected_methods <- c(
  "FASH-IWP1-Raw", "FASH-IWP1-BF",
  "FASH-linear-Raw", "FASH-linear-BF"
)
observed_keys <- paste(alpha005_summary$se_regime, alpha005_summary$method)
expected_keys <- as.vector(outer(names(se_specs), expected_methods, paste))
class_counts <- table(reference_classes)
validation <- data.frame(
  check = c(
    "exact class counts",
    "paired true beta trajectories",
    "paired SE-category assignments",
    "paired standardized noise",
    "complete regime-method summaries",
    "finite power and realized FDP",
    "resource cap"
  ),
  passed = c(
    identical(as.integer(class_counts[c("constant", "linear")]), c(800L, 200L)),
    paired_truth_max_difference <= 1e-12,
    paired_se_category_exact,
    paired_standardized_noise_max_difference <= 1e-12,
    length(observed_keys) == length(expected_keys) &&
      !anyDuplicated(observed_keys) && setequal(observed_keys, expected_keys),
    all(is.finite(alpha005_summary$power)) &&
      all(is.finite(alpha005_summary$realized_fdp)) &&
      all(alpha005_summary$power >= 0 & alpha005_summary$power <= 1) &&
      all(alpha005_summary$realized_fdp >= 0 &
        alpha005_summary$realized_fdp <= 1),
    num_cores <= 2L && Sys.getenv("VECLIB_MAXIMUM_THREADS") == "1"
  ),
  observed = c(
    paste(names(class_counts), as.integer(class_counts), collapse = "; "),
    format(paired_truth_max_difference, scientific = TRUE),
    as.character(paired_se_category_exact),
    format(paired_standardized_noise_max_difference, scientific = TRUE),
    paste(length(observed_keys), "unique keys"),
    paste(range(alpha005_summary$power), collapse = " to "),
    paste(num_cores, "workers; one BLAS thread")
  ),
  stringsAsFactors = FALSE
)
if (any(!validation$passed)) {
  print(validation)
  stop("R0 validation failed.")
}

configuration <- list(
  output_id = output_id,
  scope = "unlinked internal quick single-seed pilot",
  J = J,
  seed = seed,
  time_grid = time_grid,
  class_probs = class_probs,
  linear_generator = list(
    source = "fashr::simulate_process(type = linear, sd_poly = 1)",
    intercept_sd = 1,
    slope_sd_per_time_step = 1 / 16,
    normalize = FALSE
  ),
  se_specs = se_specs,
  paired_across_se_regimes = TRUE,
  method_settings = list(
    grid = common_grid, pred_step = common_pred_step,
    penalty = common_penalty, num_basis = 20L,
    num_cores = num_cores, alpha = 0.05,
    discovery_rule = "cumulative lfdr"
  ),
  R_version = R.version.string,
  fashr_version = as.character(utils::packageVersion("fashr"))
)
analysis_cache <- list(
  configuration = configuration,
  unit_results = unit_results,
  alpha005_summary = alpha005_summary,
  pi0_summary = pi0_summary,
  runtime = runtime,
  warnings = warnings,
  validation = validation
)

output_parent <- file.path(
  workflowr_root, "output", "revision_simulations", "internal"
)
final_dir <- file.path(output_parent, output_id)
staging_dir <- file.path(output_parent, paste0(output_id, ".staging-", Sys.getpid()))
if (dir.exists(final_dir) || dir.exists(staging_dir)) {
  stop("R0 output or staging directory already exists.")
}
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(configuration, file.path(staging_dir, "configuration.rds"))
saveRDS(analysis_cache, file.path(staging_dir, "analysis_cache.rds"))
write.csv(alpha005_summary, file.path(staging_dir, "alpha005_summary.csv"),
  row.names = FALSE)
write.csv(pi0_summary, file.path(staging_dir, "pi0_summary.csv"), row.names = FALSE)
write.csv(runtime, file.path(staging_dir, "runtime.csv"), row.names = FALSE)
write.csv(warnings, file.path(staging_dir, "warnings.csv"), row.names = FALSE)
write.csv(validation, file.path(staging_dir, "validation.csv"), row.names = FALSE)
make_power_fdp_plot(
  alpha005_summary,
  file.path(staging_dir, "power_fdp_comparison.png")
)

required <- file.path(staging_dir, c(
  "configuration.rds", "analysis_cache.rds", "alpha005_summary.csv",
  "pi0_summary.csv", "runtime.csv", "warnings.csv", "validation.csv",
  "power_fdp_comparison.png"
))
if (any(!file.exists(required)) || any(file.info(required)$size <= 0)) {
  stop("R0 outputs are incomplete.")
}
if (!file.rename(staging_dir, final_dir)) {
  stop("Could not finalize the R0 output directory.")
}

cat("\nR0 quick pilot completed in ", sprintf("%.1f", pilot_elapsed),
  " seconds.\n", sep = "")
print(alpha005_summary[
  alpha005_summary$method %in% c("FASH-IWP1-BF", "FASH-linear-BF"),
  c("se_regime", "method", "power", "realized_fdp", "n_discoveries")
])
cat("Output: ", normalizePath(final_dir, winslash = "/"), "\n", sep = "")
