# Focused C.1 rerun with J = 1000; reuse the validated generator and likelihood.
source("code/revision_simulations/internal/appendix_c1_mc10/reporting.R")

c1j_contract <- function() {
  contract <- c1_contract()
  contract$schema <- "appendix-c1-j1000-v1"
  contract$result_id <- "appendix_c1_j1000_fashr0143_20260918"
  contract$focused_J <- 1000L
  contract$grid_source <- c1_contract()$result_id
  contract
}

c1j_source_paths <- function() c(
  generator = "code/revision_simulations/appendix_b/appendix_b_helpers.R",
  integration = "code/revision_simulations/internal/appendix_c1_mc10/gaussian_likelihood.R",
  summaries = "code/revision_simulations/internal/appendix_c1_mc10/mc10_helpers.R",
  helpers = "code/revision_simulations/internal/appendix_c1_j1000/helpers.R",
  runner = "code/revision_simulations/internal/appendix_c1_j1000/run_appendix_c1_j1000.R"
)

c1j_validate <- function(alpha, contract = c1j_contract()) {
  expected <- length(contract$seeds) * 2L * 2L * length(contract$alpha)
  stopifnot(nrow(alpha) == expected,
    setequal(unique(alpha$seed), contract$seeds),
    !anyDuplicated(alpha[c("seed", "order", "stage", "alpha")]),
    all(alpha$J == 1000L), all(is.finite(alpha$power)),
    all(alpha$discoveries == alpha$true_discoveries + alpha$false_discoveries),
    all(abs(alpha$realized_fdp - alpha$false_discoveries / pmax(alpha$discoveries, 1)) < 1e-14))
  alternatives <- ifelse(alpha$order == "IWP1", 200L, 100L)
  stopifnot(all(abs(alpha$power - alpha$true_discoveries / alternatives) < 1e-14))
  for (seed in contract$seeds) for (order in c("IWP1", "IWP2"))
    for (stage in c("Raw EB", "BF updated")) {
      rows <- alpha[alpha$seed == seed & alpha$order == order & alpha$stage == stage, ]
      stopifnot(identical(sort(rows$alpha), contract$alpha))
    }
  invisible(TRUE)
}

c1j_summarize <- function(alpha) {
  result <- c1_summary(alpha, c("order", "stage", "alpha"),
    c("discoveries", "realized_fdp", "power"))
  names(result) <- sub("realized_fdp", "empirical_fdr", names(result), fixed = TRUE)
  result
}

c1j_select_examples <- function(focused) {
  # Keep four distinct examples, taking the first eligible numeric ID per panel.
  orders <- c(1L, 1L, 2L, 2L)
  classes <- c("nonlinear", "linear", "nonlinear", "nonlinear")
  selected <- character()
  truth <- focused$data_bundle$truth
  for (panel in seq_along(orders)) {
    fit <- focused$fit_bundle[[paste0("order", orders[[panel]])]]$bf
    calls <- cumulative_lfdr_calls(fit$lfdr, 0.05)$indices
    eligible <- which(truth$class == classes[[panel]] & truth$index %in% calls &
                      !truth$unit_id %in% selected)
    stopifnot(length(eligible) > 0L)
    numeric_id <- as.integer(sub("^[A-Z]", "", truth$unit_id[eligible]))
    selected <- c(selected, truth$unit_id[eligible[order(numeric_id)[1L]]])
  }
  selected
}
