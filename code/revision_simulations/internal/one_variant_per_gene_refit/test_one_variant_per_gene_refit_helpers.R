#!/usr/bin/env Rscript

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "one_variant_per_gene_refit",
  "one_variant_per_gene_refit_helpers.R"
))

synthetic_fit <- list(
  prior_weights = data.frame(psd = c(0, 1), prior_weight = c(0.5, 0.5)),
  posterior_weights = matrix(0.5, nrow = 3, ncol = 2),
  psd_grid = c(0, 1),
  lfdr = rep(0.5, 3),
  settings = list(
    num_basis = 4,
    order = 1,
    betaprec = 0,
    pred_step = 1,
    likelihood = "gaussian",
    penalty = 2
  ),
  fash_data = list(
    data_list = list(
      geneA_variant1 = data.frame(y = 1:4, x = 0:3, offset = 0),
      geneB_variant2 = data.frame(y = 5:8, x = 0:3, offset = 0),
      geneC_variant3 = data.frame(y = 9:12, x = 0:3, offset = 0)
    ),
    S = list(
      geneA_variant1 = rep(0.1, 4),
      geneB_variant2 = rep(0.2, 4),
      geneC_variant3 = rep(0.3, 4)
    ),
    Omega = NULL
  ),
  L_matrix = rbind(
    c(-0.2, -1.0),
    c(-1.0, -0.2),
    c(-0.4, -0.5)
  )
)
refit_datasets <- make_fash_refit_datasets(synthetic_fit, c(3L, 1L))
stopifnot(
  identical(names(refit_datasets), c("geneC_variant3", "geneA_variant1")),
  identical(names(refit_datasets[[1]]), c("beta", "time", "SE")),
  identical(refit_datasets[[1]]$beta, as.numeric(9:12)),
  identical(refit_datasets[[1]]$time, as.numeric(0:3)),
  identical(refit_datasets[[1]]$SE, rep(0.3, 4))
)

invalid_fit <- synthetic_fit
invalid_fit$fash_data$S[[1]][2] <- -1
invalid_error <- try(
  make_fash_refit_datasets(invalid_fit, c(1L, 2L)),
  silent = TRUE
)
stopifnot(inherits(invalid_error, "try-error"))

cached_refit <- refit_fash_from_cached_likelihood(synthetic_fit, c(3L, 1L))
manual_eb <- fashr::fash_eb_est(
  synthetic_fit$L_matrix[c(3L, 1L), , drop = FALSE],
  grid = synthetic_fit$psd_grid,
  penalty = synthetic_fit$settings$penalty
)
stopifnot(
  inherits(cached_refit, "fash"),
  identical(names(cached_refit$fash_data$data_list),
            c("geneC_variant3", "geneA_variant1")),
  identical(rownames(cached_refit$L_matrix),
            c("geneC_variant3", "geneA_variant1")),
  isTRUE(all.equal(
    unname(cached_refit$prior_weights$prior_weight),
    unname(manual_eb$prior_weight$prior_weight)
  )),
  isTRUE(all.equal(
    unname(cached_refit$posterior_weights),
    unname(manual_eb$posterior_weight)
  )),
  identical(names(cached_refit$lfdr),
            c("geneC_variant3", "geneA_variant1"))
)

unpenalized_refit <- refit_fash_from_cached_likelihood(
  synthetic_fit,
  c(3L, 1L),
  penalty = 1
)
manual_unpenalized_eb <- fashr::fash_eb_est(
  synthetic_fit$L_matrix[c(3L, 1L), , drop = FALSE],
  grid = synthetic_fit$psd_grid,
  penalty = 1
)
invalid_penalty_error <- try(
  refit_fash_from_cached_likelihood(
    synthetic_fit,
    c(3L, 1L),
    penalty = 1.5
  ),
  silent = TRUE
)
stopifnot(
  identical(unpenalized_refit$settings$penalty, 1),
  isTRUE(all.equal(
    unname(unpenalized_refit$prior_weights$prior_weight),
    unname(manual_unpenalized_eb$prior_weight$prior_weight)
  )),
  isTRUE(all.equal(
    unname(unpenalized_refit$posterior_weights),
    unname(manual_unpenalized_eb$posterior_weight)
  )),
  inherits(invalid_penalty_error, "try-error")
)

prior_comparison <- compare_prior_weights(
  data.frame(psd = c(0, 0.1), prior_weight = c(0.8, 0.2)),
  data.frame(psd = c(0, 0.2), prior_weight = c(0.7, 0.3)),
  fit_stage = "Raw"
)
stopifnot(
  identical(prior_comparison$table$psd, c(0, 0.1, 0.2)),
  isTRUE(all.equal(prior_comparison$summary$full_pi0, 0.8)),
  isTRUE(all.equal(prior_comparison$summary$thinned_pi0, 0.7)),
  isTRUE(all.equal(prior_comparison$summary$pi0_difference, -0.1)),
  isTRUE(all.equal(prior_comparison$summary$prior_total_variation, 0.3))
)

pair_keys <- paste0("gene", 1:4, "_variant", 1:4)
lfdr_comparison <- compare_paired_lfdr(
  full_lfdr = c(0.01, 0.02, 0.8, 0.9),
  thinned_lfdr = c(0.01, 0.2, 0.7, 0.95),
  pair_keys = pair_keys,
  fit_stage = "BF-adjusted",
  alpha = 0.05
)
stopifnot(
  lfdr_comparison$summary$n_units == 4L,
  lfdr_comparison$summary$full_fdr_calls == 2L,
  lfdr_comparison$summary$thinned_fdr_calls == 1L,
  lfdr_comparison$summary$fdr_call_intersection == 1L,
  lfdr_comparison$summary$fdr_call_union == 2L,
  isTRUE(all.equal(lfdr_comparison$summary$fdr_call_jaccard, 0.5)),
  identical(lfdr_comparison$table$full_fdr_call, c(TRUE, TRUE, FALSE, FALSE)),
  identical(lfdr_comparison$table$thinned_fdr_call, c(TRUE, FALSE, FALSE, FALSE))
)

empty_calls <- cumulative_fdr_calls(c(0.5, 0.6), alpha = 0.05)
stopifnot(identical(empty_calls, integer(0)))

minimum_lfdr_selection <- select_minimum_lfdr_variant_per_gene(
  pair_keys = c(
    "geneA_variant2",
    "geneA_variant1",
    "geneB_variant1",
    "geneB_variant2"
  ),
  lfdr = c(0.20, 0.20, 0.90, 0.10)
)
stopifnot(
  identical(
    minimum_lfdr_selection$pair_key,
    c("geneA_variant1", "geneB_variant2")
  ),
  identical(minimum_lfdr_selection$selection_lfdr, c(0.20, 0.10)),
  identical(
    minimum_lfdr_selection$selection_rule,
    rep("minimum_full_bf_lfdr", 2L)
  )
)

cat("One-variant-per-gene FASH refit helper tests passed.\n")
