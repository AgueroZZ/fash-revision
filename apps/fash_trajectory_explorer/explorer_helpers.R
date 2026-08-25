# Pure helpers for the local FASH trajectory explorer.

normalize_explorer_query <- function(query) {
  query <- trimws(as.character(query)[1L])
  if (is.na(query)) "" else tolower(query)
}

make_explorer_search_fields <- function(index) {
  required <- c("gene_symbol", "gene_id", "variant_id")
  if (!is.data.frame(index) || !all(required %in% names(index))) {
    stop("index is missing the columns needed for search fields.")
  }
  list(
    symbol = tolower(index$gene_symbol),
    gene = tolower(index$gene_id),
    variant = tolower(index$variant_id)
  )
}

match_explorer_text <- function(fields, query) {
  query <- normalize_explorer_query(query)
  lengths_ok <- length(unique(lengths(fields))) == 1L
  if (!is.list(fields) || length(fields) == 0L || !lengths_ok) {
    stop("fields must be a non-empty list of equal-length character vectors.")
  }
  if (query == "") {
    return(rep(TRUE, length(fields[[1L]])))
  }
  Reduce(`|`, lapply(fields, grepl, pattern = query, fixed = TRUE))
}

rank_explorer_matches <- function(index, query, limit = 200L, fields = NULL) {
  required <- c(
    "original_index", "gene_symbol", "gene_id", "variant_id",
    "discovery_status_bf"
  )
  if (!is.data.frame(index) || !all(required %in% names(index))) {
    stop("index is missing required explorer columns.")
  }
  query <- normalize_explorer_query(query)
  limit <- as.integer(limit)
  if (query == "" || is.na(limit) || limit < 1L) {
    return(index[0L, , drop = FALSE])
  }
  if (is.null(fields)) fields <- make_explorer_search_fields(index)
  if (length(fields$symbol) != nrow(index)) {
    stop("fields are not aligned with index.")
  }

  # Exact and prefix matches are subsets of the substring match, so the two
  # cheap comparisons only run on the substring candidates.
  candidate <- which(match_explorer_text(fields, query))
  if (length(candidate) == 0L) {
    return(index[0L, , drop = FALSE])
  }
  candidate_fields <- lapply(fields, `[`, candidate)
  exact <- Reduce(`|`, lapply(candidate_fields, `==`, query))
  prefix <- Reduce(`|`, lapply(candidate_fields, startsWith, prefix = query))

  rank <- ifelse(exact, 0L, ifelse(prefix, 1L, 2L))
  best_lfdr <- if (all(c("iwp1_lfdr_bf", "linear_lfdr_bf") %in% names(index))) {
    pmin(index$iwp1_lfdr_bf[candidate], index$linear_lfdr_bf[candidate])
  } else {
    rep(Inf, length(candidate))
  }
  ordering <- order(
    rank,
    best_lfdr,
    candidate_fields$symbol,
    candidate_fields$gene,
    candidate_fields$variant,
    index$original_index[candidate]
  )
  index[candidate[ordering[seq_len(min(limit, length(ordering)))]], , drop = FALSE]
}

select_distinct_gene_examples <- function(indices,
                                          pair_table,
                                          primary_score,
                                          secondary_score,
                                          n = 3L,
                                          primary_decreasing = FALSE,
                                          secondary_decreasing = TRUE) {
  indices <- as.integer(indices)
  if (length(indices) == 0L ||
      nrow(pair_table) != length(primary_score) ||
      length(primary_score) != length(secondary_score)) {
    stop("Example-selection inputs are empty or misaligned.")
  }
  first <- primary_score[indices]
  second <- secondary_score[indices]
  if (primary_decreasing) first <- -first
  if (secondary_decreasing) second <- -second
  ordering <- order(first, second, pair_table$key[indices])
  ordered_indices <- indices[ordering]
  distinct <- !duplicated(pair_table$gene_id[ordered_indices])
  selected <- ordered_indices[distinct]
  selected[seq_len(min(as.integer(n), length(selected)))]
}

format_explorer_number <- function(value, digits = 3L) {
  value <- as.numeric(value)
  ifelse(
    is.na(value),
    "NA",
    ifelse(
      value == 0,
      "0",
      ifelse(
        abs(value) < 0.001,
        formatC(value, format = "e", digits = digits),
        formatC(value, format = "f", digits = digits)
      )
    )
  )
}

make_explorer_plot_subtitle <- function(row) {
  paste0(
    "IWP1 lfdr raw/BF: ",
    format_explorer_number(row$iwp1_lfdr_raw, 3L), " / ",
    format_explorer_number(row$iwp1_lfdr_bf, 3L),
    "\nFASH-linear lfdr raw/BF: ",
    format_explorer_number(row$linear_lfdr_raw, 3L), " / ",
    format_explorer_number(row$linear_lfdr_bf, 3L),
    "\nStrober p-values: linear ",
    format_explorer_number(row$strober_linear_pvalue, 2L),
    "; nonlinear ",
    format_explorer_number(row$strober_nonlinear_pvalue, 2L)
  )
}

make_explorer_metrics <- function(row, alpha = 0.05) {
  yes_no <- function(value) ifelse(isTRUE(value), "Yes", "No")
  data.frame(
    Analysis = c(
      "IWP1 FASH", "IWP1 FASH", "FASH-linear", "FASH-linear",
      "Strober linear", "Strober nonlinear"
    ),
    Metric = c("Raw lfdr", "BF-adjusted lfdr", "Raw lfdr", "BF-adjusted lfdr", "p-value", "p-value"),
    Value = c(
      format_explorer_number(row$iwp1_lfdr_raw, 4L),
      format_explorer_number(row$iwp1_lfdr_bf, 4L),
      format_explorer_number(row$linear_lfdr_raw, 4L),
      format_explorer_number(row$linear_lfdr_bf, 4L),
      format_explorer_number(row$strober_linear_pvalue, 3L),
      format_explorer_number(row$strober_nonlinear_pvalue, 3L)
    ),
    `Adjusted metric` = c(
      "--", "--", "--", "--",
      format_explorer_number(row$strober_linear_efdr, 4L),
      format_explorer_number(row$strober_nonlinear_efdr, 4L)
    ),
    `Called at 0.05` = c(
      yes_no(row$iwp1_called_raw),
      yes_no(row$iwp1_called_bf),
      yes_no(row$linear_called_raw),
      yes_no(row$linear_called_bf),
      yes_no(row$strober_linear_efdr <= alpha),
      yes_no(row$strober_nonlinear_efdr <= alpha)
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

split_explorer_tokens <- function(query) {
  tokens <- strsplit(normalize_explorer_query(query), "[[:space:]]+")[[1L]]
  tokens[nzchar(tokens)]
}

# Multi-token queries such as "GPR78 rs4583742" require every token to match
# somewhere in the row; ranking then follows the first token.
search_explorer_index <- function(index, query, limit = 200L, fields = NULL) {
  tokens <- split_explorer_tokens(query)
  if (length(tokens) == 0L) {
    return(index[0L, , drop = FALSE])
  }
  if (is.null(fields)) fields <- make_explorer_search_fields(index)
  if (length(tokens) == 1L) {
    return(rank_explorer_matches(index, tokens[1L], limit, fields))
  }
  mask <- Reduce(`&`, lapply(tokens, match_explorer_text, fields = fields))
  keep <- which(mask)
  if (length(keep) == 0L) {
    return(index[0L, , drop = FALSE])
  }
  rank_explorer_matches(
    index[keep, , drop = FALSE],
    tokens[1L],
    limit,
    lapply(fields, `[`, keep)
  )
}

make_discovery_status <- function(iwp1_called, linear_called) {
  iwp1_called <- as.logical(iwp1_called)
  linear_called <- as.logical(linear_called)
  if (length(iwp1_called) != length(linear_called) || length(iwp1_called) == 0L ||
      anyNA(iwp1_called) || anyNA(linear_called)) {
    stop("Discovery status needs aligned complete logical call vectors.")
  }
  ifelse(
    iwp1_called & linear_called,
    "Both",
    ifelse(
      iwp1_called,
      "IWP1 only",
      ifelse(linear_called, "FASH-linear only", "Neither")
    )
  )
}

# ---------------------------------------------------------------------------
# Method registry shared by the trajectory tab and the all-pair browser tab.
# Every statistic below is smaller-is-more-significant, so one ascending order
# serves ranking and lead-variant selection alike. Strober has no raw/BF
# distinction, so both versions point at the same reported columns.
# ---------------------------------------------------------------------------

explorer_metric_catalog <- function() {
  data.frame(
    id = c("iwp1", "linear", "strober_linear", "strober_nonlinear"),
    label = c("IWP1 FASH", "FASH-linear", "Strober linear", "Strober nonlinear"),
    family = c("fash", "fash", "strober", "strober"),
    statistic = c("lfdr", "lfdr", "p-value", "p-value"),
    column_bf = c(
      "iwp1_lfdr_bf", "linear_lfdr_bf",
      "strober_linear_pvalue", "strober_nonlinear_pvalue"
    ),
    column_raw = c(
      "iwp1_lfdr_raw", "linear_lfdr_raw",
      "strober_linear_pvalue", "strober_nonlinear_pvalue"
    ),
    call_bf = c(
      "iwp1_called_bf", "linear_called_bf",
      "strober_linear_called", "strober_nonlinear_called"
    ),
    call_raw = c(
      "iwp1_called_raw", "linear_called_raw",
      "strober_linear_called", "strober_nonlinear_called"
    ),
    stringsAsFactors = FALSE
  )
}

explorer_metric_field <- function(id, version = c("bf", "raw"), what = c("column", "call")) {
  version <- match.arg(version)
  what <- match.arg(what)
  catalog <- explorer_metric_catalog()
  position <- match(as.character(id)[1L], catalog$id)
  if (is.na(position)) stop("Unknown explorer metric id: ", id)
  catalog[[paste0(what, "_", version)]][position]
}

explorer_metric_label <- function(id, version = c("bf", "raw")) {
  version <- match.arg(version)
  catalog <- explorer_metric_catalog()
  position <- match(as.character(id)[1L], catalog$id)
  if (is.na(position)) stop("Unknown explorer metric id: ", id)
  suffix <- if (catalog$statistic[position] != "lfdr") {
    " p-value"
  } else if (version == "bf") {
    " BF lfdr"
  } else {
    " raw lfdr"
  }
  paste0(catalog$label[position], suffix)
}

# ---------------------------------------------------------------------------
# All-pair browser primitives.
# ---------------------------------------------------------------------------

select_lead_variant_per_gene <- function(gene_id, metric, tie_break) {
  gene_id <- as.character(gene_id)
  metric <- as.numeric(metric)
  count <- length(gene_id)
  if (count == 0L || length(metric) != count || length(tie_break) != count) {
    stop("gene_id, metric, and tie_break must be aligned non-empty vectors.")
  }
  if (any(!is.finite(metric))) {
    stop("metric must be finite for lead-variant selection.")
  }
  ordering <- order(metric, tie_break, method = "radix")
  sort(ordering[!duplicated(gene_id[ordering])])
}

# An empty or malformed maximum-statistic box must mean "no restriction", never
# "keep nothing", so the fallback is the inclusive upper bound of 1.
sanitize_threshold <- function(value) {
  value <- suppressWarnings(as.numeric(value)[1L])
  if (length(value) == 0L || is.na(value) || value < 0) {
    return(1)
  }
  min(value, 1)
}

apply_significance_filter <- function(called, mode = c("any", "yes", "no")) {
  mode <- match.arg(mode)
  called <- as.logical(called)
  if (length(called) == 0L || anyNA(called)) {
    stop("called must be a non-empty logical vector without missing values.")
  }
  switch(
    mode,
    any = rep(TRUE, length(called)),
    yes = called,
    no = !called
  )
}

apply_discordance_filter <- function(calls,
                                     mode = c(
                                       "any", "fash_models", "strober_tests",
                                       "fash_vs_strober"
                                     )) {
  mode <- match.arg(mode)
  required <- c("iwp1", "linear", "strober_linear", "strober_nonlinear")
  if (!is.list(calls) || !all(required %in% names(calls)) ||
      length(unique(lengths(calls[required]))) != 1L ||
      any(vapply(calls[required], function(x) anyNA(x) || !is.logical(x), logical(1)))) {
    stop("calls must be equal-length complete logical vectors for all four methods.")
  }
  count <- length(calls$iwp1)
  switch(
    mode,
    any = rep(TRUE, count),
    fash_models = xor(calls$iwp1, calls$linear),
    strober_tests = xor(calls$strober_linear, calls$strober_nonlinear),
    fash_vs_strober = xor(
      calls$iwp1 | calls$linear,
      calls$strober_linear | calls$strober_nonlinear
    )
  )
}

make_call_code <- function(calls) {
  required <- c("iwp1", "linear", "strober_linear", "strober_nonlinear")
  if (!is.list(calls) || !all(required %in% names(calls)) ||
      length(unique(lengths(calls[required]))) != 1L) {
    stop("calls must contain equal-length logical vectors for all four methods.")
  }
  as.integer(calls$iwp1) +
    2L * as.integer(calls$linear) +
    4L * as.integer(calls$strober_linear) +
    8L * as.integer(calls$strober_nonlinear)
}

decode_call_code <- function(code) {
  code <- as.integer(code)
  if (any(is.na(code)) || any(code < 0L | code > 15L)) {
    stop("code must contain integers in [0, 15].")
  }
  list(
    iwp1 = bitwAnd(code, 1L) > 0L,
    linear = bitwAnd(code, 2L) > 0L,
    strober_linear = bitwAnd(code, 4L) > 0L,
    strober_nonlinear = bitwAnd(code, 8L) > 0L
  )
}

summarize_pair_view <- function(gene_id, calls) {
  required <- c("iwp1", "linear", "strober_linear", "strober_nonlinear")
  if (!is.list(calls) || !all(required %in% names(calls))) {
    stop("calls must contain logical vectors for all four methods.")
  }
  calls <- calls[required]
  if (length(gene_id) == 0L) {
    empty <- matrix(
      0L, nrow = 4L, ncol = 4L,
      dimnames = list(required, required)
    )
    return(list(
      pairs = 0L,
      genes = 0L,
      significant = stats::setNames(integer(4L), required),
      joint = empty
    ))
  }
  if (any(lengths(calls) != length(gene_id))) {
    stop("calls are not aligned with gene_id.")
  }
  membership <- vapply(calls, as.numeric, numeric(length(gene_id)))
  if (is.null(dim(membership))) {
    membership <- matrix(membership, nrow = 1L, dimnames = list(NULL, required))
  }
  joint <- crossprod(membership)
  storage.mode(joint) <- "integer"
  list(
    pairs = length(gene_id),
    genes = length(unique(gene_id)),
    significant = stats::setNames(as.integer(diag(joint)), required),
    joint = joint
  )
}

# ---------------------------------------------------------------------------
# Strober-style reference trend. Strober et al. test genotype-level data; this
# is only an inverse-variance-weighted least-squares line through the
# time-specific-PC beta-hat summary, and every caller must label it as such.
# ---------------------------------------------------------------------------

fit_weighted_linear_trend <- function(time, beta, standard_error, grid) {
  time <- as.numeric(time)
  beta <- as.numeric(beta)
  standard_error <- as.numeric(standard_error)
  grid <- as.numeric(grid)
  if (length(time) < 2L || length(beta) != length(time) ||
      length(standard_error) != length(time) || length(grid) == 0L ||
      any(!is.finite(c(time, beta, standard_error, grid))) ||
      any(standard_error <= 0) || length(unique(time)) < 2L) {
    stop("fit_weighted_linear_trend needs at least two finite points with positive standard errors.")
  }
  weight <- 1 / standard_error^2
  sum_w <- sum(weight)
  sum_wx <- sum(weight * time)
  sum_wy <- sum(weight * beta)
  sum_wxx <- sum(weight * time^2)
  sum_wxy <- sum(weight * time * beta)
  denominator <- sum_w * sum_wxx - sum_wx^2
  if (!is.finite(denominator) || abs(denominator) < .Machine$double.eps) {
    stop("The weighted design is degenerate; no linear trend is identified.")
  }
  slope <- (sum_w * sum_wxy - sum_wx * sum_wy) / denominator
  intercept <- (sum_wy - slope * sum_wx) / sum_w
  result <- data.frame(
    time = grid,
    fitted = intercept + slope * grid,
    stringsAsFactors = FALSE
  )
  attr(result, "slope") <- slope
  attr(result, "intercept") <- intercept
  result
}

# ---------------------------------------------------------------------------
# Single-pair evidence cards for the trajectory tab.
# ---------------------------------------------------------------------------

make_explorer_metric_cards <- function(row, version = c("bf", "raw"), alpha = 0.05) {
  version <- match.arg(version)
  catalog <- explorer_metric_catalog()
  primary_column <- if (version == "bf") catalog$column_bf else catalog$column_raw
  primary_value <- vapply(primary_column, function(column) {
    as.numeric(row[[column]])[1L]
  }, numeric(1), USE.NAMES = FALSE)
  secondary_label <- ifelse(
    catalog$statistic == "lfdr",
    if (version == "bf") "raw lfdr" else "BF lfdr",
    "eFDR"
  )
  secondary_column <- ifelse(
    catalog$statistic == "lfdr",
    if (version == "bf") catalog$column_raw else catalog$column_bf,
    c("", "", "strober_linear_efdr", "strober_nonlinear_efdr")
  )
  secondary_value <- vapply(secondary_column, function(column) {
    as.numeric(row[[column]])[1L]
  }, numeric(1), USE.NAMES = FALSE)
  called <- c(
    isTRUE(row[[if (version == "bf") "iwp1_called_bf" else "iwp1_called_raw"]]),
    isTRUE(row[[if (version == "bf") "linear_called_bf" else "linear_called_raw"]]),
    isTRUE(as.numeric(row$strober_linear_efdr)[1L] <= alpha),
    isTRUE(as.numeric(row$strober_nonlinear_efdr)[1L] <= alpha)
  )
  data.frame(
    id = catalog$id,
    label = catalog$label,
    family = catalog$family,
    primary_label = ifelse(
      catalog$statistic == "lfdr",
      if (version == "bf") "BF lfdr" else "raw lfdr",
      "p-value"
    ),
    primary_value = primary_value,
    secondary_label = secondary_label,
    secondary_value = secondary_value,
    called = called,
    stringsAsFactors = FALSE
  )
}
