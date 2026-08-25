script_argument <- commandArgs()[grep("^--file=", commandArgs())][[1L]]
script_path <- normalizePath(
  sub("^--file=", "", script_argument),
  winslash = "/",
  mustWork = TRUE
)
repository_root <- normalizePath(
  file.path(dirname(script_path), "..", "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

expected_fashr_version <- "0.1.43"
expected_fashr_remote_sha <- "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"

if (!requireNamespace("fashr", quietly = TRUE)) {
  stop("The fashr package is not installed.", call. = FALSE)
}

description <- utils::packageDescription("fashr")
actual_version <- as.character(utils::packageVersion("fashr"))
actual_remote_sha <- if (is.null(description$RemoteSha)) {
  NA_character_
} else {
  as.character(description$RemoteSha)
}

stopifnot(
  identical(actual_version, expected_fashr_version),
  identical(actual_remote_sha, expected_fashr_remote_sha)
)

source(file.path(
  repository_root,
  "code",
  "revision_simulations",
  "shared",
  "simulation_functions.R"
))

full_fit_paths <- c(
  r1 = file.path(
    repository_root,
    "output",
    "revision_simulations",
    "mc",
    paste0(
      "r1_real_genotype_one_per_gene_J6362_random_bspline_main_effect_",
      "linear_mixture_predstep1_penalty10_pilot5"
    ),
    "full_fits",
    "seed_12345.rds"
  ),
  r2 = file.path(
    repository_root,
    "output",
    "revision_simulations",
    "mc",
    paste0(
      "r2_real_genotype_one_per_gene_J6362_timed_cosine_one_two_three_",
      "peak_main_effect_linear_mixture_predstep1_penalty10_pilot5"
    ),
    "full_fits",
    "seed_12345.rds"
  )
)

missing_paths <- full_fit_paths[!file.exists(full_fit_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "Missing formal seed-12345 full-fit input: ",
    paste(missing_paths, collapse = ", "),
    call. = FALSE
  )
}

expected_bf_discoveries <- list(
  r1 = c(`FASH-IWP1-BF` = 1213L, `FASH-linear-BF` = 625L),
  r2 = c(`FASH-IWP1-BF` = 688L, `FASH-linear-BF` = 32L)
)

extract_null_weight <- function(fit) {
  stopifnot(
    is.numeric(fit$psd_grid),
    sum(fit$psd_grid == 0) == 1L,
    is.data.frame(fit$prior_weights),
    all(c("psd", "prior_weight") %in% names(fit$prior_weights))
  )
  null_rows <- which(fit$prior_weights$psd == 0)
  if (length(null_rows) == 0L) {
    return(0)
  }
  stopifnot(length(null_rows) == 1L)
  weight <- as.numeric(fit$prior_weights$prior_weight[[null_rows]])
  stopifnot(is.finite(weight), weight >= 0, weight <= 1)
  weight
}

for (scenario_id in names(full_fit_paths)) {
  retained <- readRDS(full_fit_paths[[scenario_id]])
  stopifnot(
    is.list(retained),
    nrow(retained$unit_info) == 6362L,
    all(c("fash_iwp1_raw", "fash_iwp1_bf") %in% names(retained$fash_fits))
  )

  iwp_raw <- retained$fash_fits$fash_iwp1_raw
  iwp_bf <- fashr::BF_update(iwp_raw, plot = FALSE)
  linear_raw <- retained$simplified_fit
  linear_bf <- BF_update_linear_mixture_fash(linear_raw)

  fits <- list(
    `FASH-IWP1-Raw` = iwp_raw,
    `FASH-IWP1-BF` = iwp_bf,
    `FASH-linear-Raw` = linear_raw,
    `FASH-linear-BF` = linear_bf
  )
  null_weights <- vapply(fits, extract_null_weight, numeric(1))
  stopifnot(
    all(is.finite(null_weights)),
    all(vapply(fits, function(fit) {
      lfdr <- get_fash_lfdr(fit)
      length(lfdr) == 6362L && all(is.finite(lfdr))
    }, logical(1)))
  )

  result_table <- rbind(
    evaluate_lfdr_method(
      lfdr = get_fash_lfdr(iwp_raw),
      unit_info = retained$unit_info,
      method = "FASH-IWP1-Raw",
      target = "dynamic",
      alpha = 0.05
    ),
    evaluate_lfdr_method(
      lfdr = get_fash_lfdr(iwp_bf),
      unit_info = retained$unit_info,
      method = "FASH-IWP1-BF",
      target = "dynamic",
      alpha = 0.05
    ),
    evaluate_simplified_fash_fit(
      fit = linear_raw,
      unit_info = retained$unit_info,
      alpha = 0.05,
      method = "FASH-linear-Raw"
    ),
    evaluate_simplified_fash_fit(
      fit = linear_bf,
      unit_info = retained$unit_info,
      alpha = 0.05,
      method = "FASH-linear-BF"
    )
  )
  alpha_curve <- compute_alpha_curve(result_table, alpha_grid = 0.05)

  stopifnot(
    nrow(result_table) == 4L * 6362L,
    nrow(alpha_curve) == 4L,
    setequal(names(expected_bf_discoveries[[scenario_id]]), alpha_curve$method[
      grepl("-BF$", alpha_curve$method)
    ])
  )

  observed_bf <- alpha_curve$n_discoveries[
    match(names(expected_bf_discoveries[[scenario_id]]), alpha_curve$method)
  ]
  stopifnot(identical(
    as.integer(observed_bf),
    unname(expected_bf_discoveries[[scenario_id]])
  ))

  cat(
    scenario_id,
    ": corrected BF discoveries = ",
    paste(names(expected_bf_discoveries[[scenario_id]]), observed_bf, sep = ":", collapse = ", "),
    "; exact-null weights = ",
    paste(names(null_weights), sprintf("%.6f", null_weights), sep = ":", collapse = ", "),
    "\n",
    sep = ""
  )
}

cat("R1/R2 retained-fit contract test passed.\n")
