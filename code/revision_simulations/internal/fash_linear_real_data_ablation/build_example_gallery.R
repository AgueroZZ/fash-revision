#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

find_workflowr_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "_workflowr.yml"))) return(current)
    parent <- dirname(current)
    if (identical(parent, current)) stop("Could not locate the workflowr project root.")
    current <- parent
  }
}

load_exact_object <- function(path, object_name) {
  environment <- new.env(parent = emptyenv())
  loaded <- load(path, envir = environment)
  if (!identical(loaded, object_name)) {
    stop("Expected only object ", object_name, " in ", path)
  }
  environment[[object_name]]
}

make_provenance <- function(paths) {
  information <- file.info(paths)
  data.frame(
    label = names(paths),
    path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    bytes = as.numeric(information$size),
    modified = format(information$mtime, tz = "UTC", usetz = TRUE),
    md5 = unname(tools::md5sum(paths)),
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
analysis_directory <- file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "fash_linear_real_data_ablation"
)
app_directory <- file.path(workflowr_root, "apps", "fash_trajectory_explorer")
helper_path <- file.path(analysis_directory, "fash_linear_real_data_helpers.R")
explorer_helper_path <- file.path(app_directory, "explorer_helpers.R")
shared_path <- file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
)
source(shared_path)
source(helper_path)
source(explorer_helper_path)
if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is required for IWP1 posterior prediction.")
}

output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "fash_linear_real_data_ablation_mixture_predstep1_penalty10"
)
explorer_output_directory <- file.path(
  workflowr_root,
  "output", "revision_simulations", "internal",
  "fash_trajectory_explorer_mixture_predstep1_penalty10"
)
paths <- c(
  current_bf = file.path(
    workflowr_root, "output", "dynamic_eQTL_real", "fash_fit1_update.RData"
  ),
  linear_bf = file.path(output_directory, "linear_fit_bf.rds"),
  linear_run_status = file.path(output_directory, "run_status.rds"),
  explorer_index = file.path(explorer_output_directory, "explorer_index.rds"),
  shared_functions = shared_path,
  helper = helper_path,
  explorer_helper = explorer_helper_path,
  builder = file.path(analysis_directory, "build_example_gallery.R")
)
if (any(!file.exists(paths))) stop("At least one required gallery input is missing.")

linear_run_status <- readRDS(paths[["linear_run_status"]])
if (!identical(
      linear_run_status$analysis_id,
      paste0(
        "revision_internal_fash_linear_real_data_ablation_",
        "mixture_predstep1_penalty10"
      )
    ) ||
    !identical(linear_run_status$status, "complete") ||
    !identical(
      unname(tools::md5sum(paths[["linear_bf"]])),
      unname(linear_run_status$fit_md5[["linear_bf"]])
    )) {
  stop("The current-PC mixture run is incomplete or its BF fit changed.")
}

started <- proc.time()[["elapsed"]]
message("[1/3] Loading the compact pair index and BF-adjusted fits.")
explorer_cache <- readRDS(paths[["explorer_index"]])
explorer_index <- explorer_cache$index
explorer_presets <- explorer_cache$presets
current_bf <- load_exact_object(paths[["current_bf"]], "fash_fit1_update")
linear_bf <- readRDS(paths[["linear_bf"]])
validate_compact_linear_mixture_fash(
  linear_bf,
  expected_grid = default_revision_grid(),
  expected_pred_step = 1,
  expected_penalty = 10L
)
if (!identical(
      explorer_cache$configuration$analysis_id,
      "fash_trajectory_explorer_mixture_predstep1_penalty10"
    ) ||
    !identical(
      explorer_cache$configuration$linear_output_id,
      basename(output_directory)
    ) ||
    !isTRUE(all.equal(
      explorer_cache$configuration$pred_step,
      1,
      tolerance = 0
    )) ||
    !identical(as.integer(explorer_cache$configuration$penalty), 10L) ||
    !isTRUE(linear_bf$bf_adjusted) ||
    nrow(explorer_index) != 1009173L ||
    !identical(names(current_bf$fash_data$data_list), explorer_index$key) ||
    !identical(linear_bf$unit_ids, explorer_index$key) ||
    !isTRUE(all.equal(
      as.numeric(current_bf$lfdr),
      explorer_index$iwp1_lfdr_bf,
      tolerance = 0
    )) ||
    !isTRUE(all.equal(
      as.numeric(linear_bf$lfdr),
      explorer_index$linear_lfdr_bf,
      tolerance = 0
    ))) {
  stop("The gallery inputs are not exactly pair-key and lfdr aligned.")
}

message("[2/3] Computing six deterministic discordant trajectories.")
gallery_presets <- explorer_presets[
  explorer_presets$group %in% c(
    "Strong FASH-linear-only",
    "Strong IWP1-only"
  ) & explorer_presets$rank <= 3L,
  ,
  drop = FALSE
]
gallery_presets$direction <- ifelse(
  gallery_presets$group == "Strong FASH-linear-only",
  "FASH-linear only",
  "IWP1 only"
)
gallery_presets$direction <- factor(
  gallery_presets$direction,
  levels = c("FASH-linear only", "IWP1 only")
)
gallery_presets <- gallery_presets[order(
  gallery_presets$direction,
  gallery_presets$rank
), , drop = FALSE]
gallery_indices <- gallery_presets$original_index
gallery_summary <- cbind(
  gallery_presets[, c("direction", "rank"), drop = FALSE],
  explorer_index[gallery_indices, c(
    "original_index", "key", "gene_id", "gene_symbol", "variant_id",
    "iwp1_lfdr_bf", "linear_lfdr_bf", "strober_linear_pvalue",
    "strober_nonlinear_pvalue"
  ), drop = FALSE]
)
gallery_summary$facet_label <- paste0(
  gallery_summary$gene_symbol, " / ", gallery_summary$variant_id,
  "\n", gallery_summary$direction,
  "\nlfdr IWP1/linear: ",
  format_explorer_number(gallery_summary$iwp1_lfdr_bf, 3L), " / ",
  format_explorer_number(gallery_summary$linear_lfdr_bf, 3L),
  "\nStrober p linear/nonlinear: ",
  format_explorer_number(gallery_summary$strober_linear_pvalue, 2L), " / ",
  format_explorer_number(gallery_summary$strober_nonlinear_pvalue, 2L)
)

grid <- seq(0, 15, by = 0.1)
gallery_results <- lapply(seq_len(nrow(gallery_summary)), function(row_index) {
  original_index <- gallery_summary$original_index[row_index]
  dataset <- current_bf$fash_data$data_list[[original_index]]
  standard_error <- as.numeric(current_bf$fash_data$S[[original_index]])
  set.seed(20260810L + original_index %% 1000000L)
  iwp_prediction <- predict(
    current_bf,
    index = original_index,
    smooth_var = grid
  )
  iwp <- data.frame(
    original_index = original_index,
    facet_label = gallery_summary$facet_label[row_index],
    method = "IWP1 FASH",
    time = iwp_prediction$x,
    posterior_mean = iwp_prediction$mean,
    lower = iwp_prediction$lower,
    upper = iwp_prediction$upper,
    stringsAsFactors = FALSE
  )
  linear <- extract_linear_mixture_posterior_plot_data(
    dataset = dataset,
    standard_error = standard_error,
    fit = linear_bf,
    unit_index = original_index,
    grid = grid,
    sample_size = 10000L,
    seed = 20260810L + original_index %% 1000000L
  )
  linear$original_index <- original_index
  linear$facet_label <- gallery_summary$facet_label[row_index]
  linear$method <- "FASH-linear"
  linear <- linear[, names(iwp), drop = FALSE]
  observed <- data.frame(
    original_index = original_index,
    facet_label = gallery_summary$facet_label[row_index],
    time = as.numeric(dataset$x),
    beta = as.numeric(dataset$y),
    standard_error = standard_error,
    stringsAsFactors = FALSE
  )
  list(observed = observed, iwp = iwp, linear = linear)
})

gallery_observed <- do.call(rbind, lapply(gallery_results, `[[`, "observed"))
gallery_iwp <- do.call(rbind, lapply(gallery_results, `[[`, "iwp"))
gallery_linear <- do.call(rbind, lapply(gallery_results, `[[`, "linear"))
gallery_levels <- gallery_summary$facet_label
gallery_observed$facet_label <- factor(
  gallery_observed$facet_label,
  levels = gallery_levels
)
gallery_iwp$facet_label <- factor(gallery_iwp$facet_label, levels = gallery_levels)
gallery_linear$facet_label <- factor(
  gallery_linear$facet_label,
  levels = gallery_levels
)

if (nrow(gallery_summary) != 6L ||
    length(unique(gallery_summary$gene_id)) != 6L ||
    !identical(as.integer(gallery_summary$rank), rep(1:3, 2L)) ||
    nrow(gallery_observed) != 96L || nrow(gallery_iwp) != 906L ||
    nrow(gallery_linear) != 906L ||
    any(!is.finite(as.matrix(gallery_iwp[, c(
      "time", "posterior_mean", "lower", "upper"
    )]))) ||
    any(!is.finite(as.matrix(gallery_linear[, c(
      "time", "posterior_mean", "lower", "upper"
    )])))) {
  stop("The computed trajectory gallery failed invariant checks.")
}

message("[3/3] Saving the gallery-only cache.")
cache <- list(
  configuration = list(
    analysis_id = paste0(
      "fash_linear_discordant_example_gallery_",
      "mixture_predstep1_penalty10"
    ),
    pair_count = 6L,
    linear_only_count = 3L,
    iwp1_only_count = 3L,
    selection = "Strongest BF-adjusted calls with distinct genes",
    linear_model = "Predictive-SD finite-mixture FASH-linear",
    pred_step = 1,
    penalty = 10L,
    time_grid = grid
  ),
  provenance = make_provenance(paths),
  gallery_summary = gallery_summary,
  gallery_observed = gallery_observed,
  gallery_iwp = gallery_iwp,
  gallery_linear = gallery_linear,
  runtime_seconds = proc.time()[["elapsed"]] - started
)
temporary_path <- file.path(output_directory, "example_gallery.rds.tmp")
saveRDS(cache, temporary_path)
output_path <- file.path(output_directory, "example_gallery.rds")
if (!file.rename(temporary_path, output_path)) {
  stop("Could not atomically replace the example-gallery cache.")
}
message(
  "Gallery cache saved: ", output_path, " (",
  format(round(cache$runtime_seconds, 1), nsmall = 1), " seconds)."
)
print(gallery_summary[, c(
  "direction", "rank", "gene_symbol", "variant_id", "iwp1_lfdr_bf",
  "linear_lfdr_bf", "strober_linear_pvalue", "strober_nonlinear_pvalue"
)])
