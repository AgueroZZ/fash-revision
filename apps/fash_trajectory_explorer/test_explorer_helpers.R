#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_argument <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(script_argument) == 1L) {
  sub("^--file=", "", script_argument)
} else {
  file.path(getwd(), "test_explorer_helpers.R")
}
source(file.path(dirname(normalizePath(script_path)), "explorer_helpers.R"))

toy_index <- data.frame(
  original_index = 1:5,
  gene_symbol = c("GPR78", "GPR78", "GPR7", "ABCA1", "ENSGX"),
  gene_id = c("ENSG001", "ENSG001", "ENSG002", "ENSG003", "ENSG004"),
  variant_id = c("rs4583742", "rs100", "rs101", "rs102", "rs103"),
  discovery_status_bf = c(
    "FASH-linear only", "Neither", "IWP1 only", "Both", "Neither"
  ),
  iwp1_lfdr_bf = c(0.5, 0.4, 0.01, 0.02, 0.9),
  linear_lfdr_bf = c(0.001, 0.3, 0.5, 0.02, 0.8),
  stringsAsFactors = FALSE
)

stopifnot(
  normalize_explorer_query("  GPR78 ") == "gpr78",
  nrow(rank_explorer_matches(toy_index, "")) == 0L,
  identical(rank_explorer_matches(toy_index, "RS4583742")$original_index, 1L),
  identical(rank_explorer_matches(toy_index, "ensg003")$original_index, 4L),
  identical(rank_explorer_matches(toy_index, "gpr78")$original_index, c(1L, 2L)),
  identical(rank_explorer_matches(toy_index, "gpr", limit = 2L)$original_index, c(1L, 3L)),
  nrow(rank_explorer_matches(toy_index, "missing")) == 0L
)

pair_table <- data.frame(
  key = paste0("g", c(1, 1, 2, 3, 4), "_rs", 1:5),
  gene_id = paste0("g", c(1, 1, 2, 3, 4)),
  variant_id = paste0("rs", 1:5),
  stringsAsFactors = FALSE
)
selected <- select_distinct_gene_examples(
  indices = 1:5,
  pair_table = pair_table,
  primary_score = c(0.01, 0.001, 0.02, 0.03, 0.04),
  secondary_score = c(0.5, 0.9, 0.8, 0.7, 0.6),
  n = 3L
)
stopifnot(identical(selected, c(2L, 3L, 4L)))

row <- data.frame(
  iwp1_lfdr_raw = 0.0012,
  iwp1_lfdr_bf = 0.02,
  linear_lfdr_raw = 0.4,
  linear_lfdr_bf = 0.6,
  strober_linear_pvalue = 1e-8,
  strober_nonlinear_pvalue = 0.2,
  strober_linear_efdr = 0.01,
  strober_nonlinear_efdr = 0.4,
  iwp1_called_raw = TRUE,
  iwp1_called_bf = TRUE,
  linear_called_raw = FALSE,
  linear_called_bf = FALSE
)
stopifnot(
  grepl("IWP1 lfdr raw/BF", make_explorer_plot_subtitle(row), fixed = TRUE),
  nrow(make_explorer_metrics(row)) == 6L,
  format_explorer_number(1e-8, 2L) == "1.00e-08"
)

# Precomputed search fields must give byte-identical results to the on-the-fly
# path, because the app caches them once at startup for the live search.
toy_fields <- make_explorer_search_fields(toy_index)
stopifnot(
  identical(
    rank_explorer_matches(toy_index, "gpr78"),
    rank_explorer_matches(toy_index, "gpr78", fields = toy_fields)
  ),
  identical(match_explorer_text(toy_fields, ""), rep(TRUE, 5L)),
  identical(match_explorer_text(toy_fields, "GPR"), c(TRUE, TRUE, TRUE, FALSE, FALSE)),
  inherits(try(rank_explorer_matches(toy_index, "gpr", fields = lapply(toy_fields, `[`, 1:2)), silent = TRUE), "try-error")
)

# Multi-token search and discovery status.
stopifnot(
  identical(split_explorer_tokens("  GPR78   rs4583742 "), c("gpr78", "rs4583742")),
  identical(split_explorer_tokens("   "), character(0)),
  nrow(search_explorer_index(toy_index, "")) == 0L,
  identical(
    search_explorer_index(toy_index, "gpr78", fields = toy_fields)$original_index,
    c(1L, 2L)
  ),
  identical(
    search_explorer_index(toy_index, "gpr78 rs4583742", fields = toy_fields)$original_index,
    1L
  ),
  nrow(search_explorer_index(toy_index, "gpr78 rs102", fields = toy_fields)) == 0L,
  identical(
    make_discovery_status(
      c(TRUE, TRUE, FALSE, FALSE),
      c(TRUE, FALSE, TRUE, FALSE)
    ),
    c("Both", "IWP1 only", "FASH-linear only", "Neither")
  ),
  inherits(try(make_discovery_status(TRUE, c(TRUE, FALSE)), silent = TRUE), "try-error")
)

# Metric registry.
catalog <- explorer_metric_catalog()
stopifnot(
  nrow(catalog) == 4L,
  identical(catalog$id, c("iwp1", "linear", "strober_linear", "strober_nonlinear")),
  explorer_metric_field("iwp1", "bf", "column") == "iwp1_lfdr_bf",
  explorer_metric_field("iwp1", "raw", "call") == "iwp1_called_raw",
  explorer_metric_field("strober_linear", "raw", "column") ==
    explorer_metric_field("strober_linear", "bf", "column"),
  explorer_metric_label("linear", "raw") == "FASH-linear raw lfdr",
  explorer_metric_label("strober_nonlinear", "bf") == "Strober nonlinear p-value",
  inherits(try(explorer_metric_field("nope"), silent = TRUE), "try-error")
)

# Lead-variant selection: minimum metric per gene, ties broken on ascending key.
lead_gene <- c("g1", "g1", "g1", "g2", "g2", "g3")
lead_key <- c("g1_rsB", "g1_rsA", "g1_rsC", "g2_rsB", "g2_rsA", "g3_rsA")
lead_metric <- c(0.10, 0.10, 0.90, 0.20, 0.05, 0.50)
lead <- select_lead_variant_per_gene(lead_gene, lead_metric, lead_key)
stopifnot(
  identical(lead, c(2L, 5L, 6L)),
  identical(lead_key[lead], c("g1_rsA", "g2_rsA", "g3_rsA")),
  identical(
    select_lead_variant_per_gene(lead_gene, lead_metric, seq_along(lead_gene)),
    c(1L, 5L, 6L)
  ),
  inherits(try(select_lead_variant_per_gene(lead_gene, c(NA, lead_metric[-1]), lead_key), silent = TRUE), "try-error"),
  inherits(try(select_lead_variant_per_gene(lead_gene, lead_metric[-1], lead_key), silent = TRUE), "try-error")
)

# A blank or nonsensical threshold box must widen to 1, never collapse to 0.
stopifnot(
  sanitize_threshold(NULL) == 1,
  sanitize_threshold(NA) == 1,
  sanitize_threshold("") == 1,
  sanitize_threshold(-0.5) == 1,
  sanitize_threshold(2) == 1,
  sanitize_threshold(0.05) == 0.05,
  sanitize_threshold(0) == 0
)

# Three-state significance filter.
called_vector <- c(TRUE, FALSE, TRUE)
stopifnot(
  identical(apply_significance_filter(called_vector, "any"), rep(TRUE, 3L)),
  identical(apply_significance_filter(called_vector, "yes"), called_vector),
  identical(apply_significance_filter(called_vector, "no"), !called_vector),
  inherits(try(apply_significance_filter(c(TRUE, NA), "yes"), silent = TRUE), "try-error")
)

# Discordance filter, call bitmask, and view summary all read the same calls.
toy_calls <- list(
  iwp1 = c(TRUE, FALSE, TRUE, FALSE),
  linear = c(TRUE, TRUE, FALSE, FALSE),
  strober_linear = c(FALSE, TRUE, FALSE, FALSE),
  strober_nonlinear = c(FALSE, FALSE, FALSE, TRUE)
)
stopifnot(
  identical(apply_discordance_filter(toy_calls, "any"), rep(TRUE, 4L)),
  identical(
    apply_discordance_filter(toy_calls, "fash_models"),
    c(FALSE, TRUE, TRUE, FALSE)
  ),
  identical(
    apply_discordance_filter(toy_calls, "strober_tests"),
    c(FALSE, TRUE, FALSE, TRUE)
  ),
  identical(
    apply_discordance_filter(toy_calls, "fash_vs_strober"),
    c(TRUE, FALSE, TRUE, TRUE)
  ),
  identical(make_call_code(toy_calls), c(3L, 6L, 1L, 8L)),
  identical(decode_call_code(make_call_code(toy_calls)), toy_calls),
  inherits(try(decode_call_code(16L), silent = TRUE), "try-error"),
  inherits(try(apply_discordance_filter(toy_calls[1:2], "fash_models"), silent = TRUE), "try-error")
)

summary_all <- summarize_pair_view(c("g1", "g1", "g2", "g3"), toy_calls)
summary_one <- summarize_pair_view("g1", lapply(toy_calls, `[`, 1L))
summary_none <- summarize_pair_view(character(0), lapply(toy_calls, `[`, integer(0)))
stopifnot(
  summary_all$pairs == 4L,
  summary_all$genes == 3L,
  identical(
    summary_all$significant,
    c(iwp1 = 2L, linear = 2L, strober_linear = 1L, strober_nonlinear = 1L)
  ),
  summary_all$joint["iwp1", "linear"] == 1L,
  summary_all$joint["iwp1", "strober_linear"] == 0L,
  identical(diag(summary_all$joint), summary_all$significant),
  summary_one$pairs == 1L && summary_one$genes == 1L,
  identical(summary_one$significant["iwp1"], c(iwp1 = 1L)),
  summary_none$pairs == 0L && sum(summary_none$joint) == 0L
)

# The Strober-style overlay must recover an exactly linear signal and must
# weight precise time points more heavily than imprecise ones.
trend_time <- c(0, 1, 2, 3)
trend <- fit_weighted_linear_trend(
  time = trend_time,
  beta = 0.5 + 0.25 * trend_time,
  standard_error = rep(0.1, 4L),
  grid = c(0, 4)
)
noisy_trend <- fit_weighted_linear_trend(
  time = trend_time,
  beta = c(0.5, 0.75, 1.0, 5.0),
  standard_error = c(0.01, 0.01, 0.01, 10),
  grid = trend_time
)
stopifnot(
  isTRUE(all.equal(attr(trend, "slope"), 0.25)),
  isTRUE(all.equal(attr(trend, "intercept"), 0.5)),
  isTRUE(all.equal(trend$fitted, c(0.5, 1.5))),
  isTRUE(all.equal(attr(noisy_trend, "slope"), 0.25, tolerance = 1e-3)),
  inherits(try(fit_weighted_linear_trend(0, 1, 0.1, 0:1), silent = TRUE), "try-error"),
  inherits(try(fit_weighted_linear_trend(c(1, 1), c(1, 2), c(0.1, 0.1), 0:1), silent = TRUE), "try-error"),
  inherits(try(fit_weighted_linear_trend(0:1, c(1, 2), c(0, 0.1), 0:1), silent = TRUE), "try-error")
)

# Evidence cards: primary/secondary statistics swap with the lfdr version, and
# the chip never disagrees with the statistic it sits next to.
cards_bf <- make_explorer_metric_cards(row, "bf")
cards_raw <- make_explorer_metric_cards(row, "raw")
stopifnot(
  nrow(cards_bf) == 4L,
  identical(cards_bf$primary_value, c(0.02, 0.6, 1e-8, 0.2)),
  identical(cards_bf$secondary_value, c(0.0012, 0.4, 0.01, 0.4)),
  identical(cards_bf$called, c(TRUE, FALSE, TRUE, FALSE)),
  identical(cards_bf$primary_label, c("BF lfdr", "BF lfdr", "p-value", "p-value")),
  identical(cards_raw$primary_value, c(0.0012, 0.4, 1e-8, 0.2)),
  identical(cards_raw$secondary_value, c(0.02, 0.6, 0.01, 0.4)),
  identical(cards_raw$called, c(TRUE, FALSE, TRUE, FALSE)),
  identical(cards_bf$family, c("fash", "fash", "strober", "strober"))
)

cat("FASH trajectory explorer helper tests passed.\n")
