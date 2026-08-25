# Scientific helpers for the internal current-FASH versus FASH-CL manuscript page.

parse_pair_keys <- function(keys) {
  keys <- as.character(keys)
  separator <- regexpr("_", keys, fixed = TRUE)
  if (any(separator < 2L)) {
    stop("Every pair key must contain a gene-variant separator.")
  }
  data.frame(
    key = keys,
    gene_id = substr(keys, 1L, separator - 1L),
    variant_id = substr(keys, separator + 1L, nchar(keys)),
    stringsAsFactors = FALSE
  )
}

extract_full_prior_weights <- function(fit) {
  if (is.null(fit$psd_grid) || is.null(fit$prior_weights) ||
      !all(c("psd", "prior_weight") %in% names(fit$prior_weights))) {
    stop("The FASH fit does not contain a valid prior-weight table.")
  }
  psd_grid <- as.numeric(fit$psd_grid)
  retained <- fit$prior_weights[, c("psd", "prior_weight"), drop = FALSE]
  indices <- match(as.numeric(retained$psd), psd_grid)
  if (anyNA(indices) || anyDuplicated(indices) ||
      any(!is.finite(retained$prior_weight)) ||
      any(retained$prior_weight < 0)) {
    stop("The retained prior weights do not align with the PSD grid.")
  }
  weights <- numeric(length(psd_grid))
  weights[indices] <- as.numeric(retained$prior_weight)
  if (abs(sum(weights) - 1) > 1e-7) {
    stop("The aligned prior weights do not sum to one.")
  }
  data.frame(
    psd = psd_grid,
    prior_weight = weights,
    stringsAsFactors = FALSE
  )
}

select_top_pair_indices_per_gene <- function(gene_id,
                                             variant_id,
                                             lfdr,
                                             pair_key = NULL) {
  gene_id <- as.character(gene_id)
  variant_id <- as.character(variant_id)
  lfdr <- as.numeric(lfdr)
  if (is.null(pair_key)) {
    pair_key <- paste0(gene_id, "_", variant_id)
  }
  pair_key <- as.character(pair_key)
  if (length(gene_id) == 0L ||
      !identical(length(gene_id), length(variant_id)) ||
      !identical(length(gene_id), length(lfdr)) ||
      !identical(length(gene_id), length(pair_key)) ||
      anyNA(gene_id) || anyNA(variant_id) || anyNA(lfdr) ||
      any(!is.finite(lfdr))) {
    stop("Top-pair inputs must be aligned, finite, and nonempty.")
  }
  ordering <- order(gene_id, lfdr, variant_id, pair_key)
  selected <- ordering[!duplicated(gene_id[ordering])]
  selected[order(gene_id[selected])]
}

row_log_sum_exp <- function(x) {
  x <- as.matrix(x)
  maximum <- matrixStats::rowMaxs(x)
  maximum + log(rowSums(exp(x - maximum)))
}

select_cumulative_lfdr_calls <- function(lfdr, alpha = 0.05) {
  lfdr <- as.numeric(lfdr)
  if (length(lfdr) == 0L || any(!is.finite(lfdr)) ||
      any(lfdr < 0 | lfdr > 1)) {
    stop("lfdr must be a non-empty finite vector in [0, 1].")
  }
  ordering <- order(lfdr, seq_along(lfdr))
  cumulative_fdr <- cumsum(lfdr[ordering]) / seq_along(ordering)
  selected_rank <- which(cumulative_fdr <= alpha)
  if (length(selected_rank) == 0L) {
    return(integer())
  }
  ordering[seq_len(max(selected_rank))]
}

build_adjusted_subfit <- function(raw_fit,
                                  selected_indices,
                                  null_weight,
                                  alternative_weights) {
  selected_indices <- sort(unique(as.integer(selected_indices)))
  if (!inherits(raw_fit, "fash")) {
    stop("raw_fit must inherit from class fash.")
  }
  if (length(selected_indices) == 0L ||
      min(selected_indices) < 1L ||
      max(selected_indices) > nrow(raw_fit$L_matrix)) {
    stop("selected_indices are empty or out of range.")
  }
  if (!is.finite(null_weight) || null_weight <= 0 || null_weight >= 1) {
    stop("null_weight must be strictly between zero and one.")
  }
  alternative_weights <- as.numeric(alternative_weights)
  if (length(alternative_weights) != ncol(raw_fit$L_matrix) - 1L ||
      any(!is.finite(alternative_weights)) ||
      any(alternative_weights < 0) ||
      abs(sum(alternative_weights) - 1) > 1e-8) {
    stop("alternative_weights must be a normalized non-null mixture.")
  }

  prior <- c(null_weight, (1 - null_weight) * alternative_weights)
  selected_likelihood <- raw_fit$L_matrix[selected_indices, , drop = FALSE]
  active <- prior > 0
  log_posterior <- sweep(
    selected_likelihood[, active, drop = FALSE],
    2L,
    log(prior[active]),
    `+`
  )
  log_normalizer <- row_log_sum_exp(log_posterior)
  posterior <- exp(log_posterior - log_normalizer)
  pair_keys <- names(raw_fit$fash_data$data_list)[selected_indices]
  rownames(posterior) <- pair_keys
  colnames(posterior) <- as.character(raw_fit$psd_grid[active])

  fit <- raw_fit
  fit$L_matrix <- selected_likelihood
  rownames(fit$L_matrix) <- pair_keys
  fit$posterior_weights <- posterior
  fit$prior_weights <- data.frame(
    psd = raw_fit$psd_grid[active],
    prior_weight = prior[active],
    row.names = NULL
  )
  fit$lfdr <- posterior[, 1L]
  fit$fash_data$data_list <- raw_fit$fash_data$data_list[selected_indices]
  fit$fash_data$S <- raw_fit$fash_data$S[selected_indices]
  if (!is.null(raw_fit$fash_data$Omega)) {
    fit$fash_data$Omega <- raw_fit$fash_data$Omega[selected_indices]
  }
  fit$original_indices <- selected_indices
  class(fit) <- class(raw_fit)
  fit
}

subset_fash_fit <- function(fit, selected_indices) {
  selected_indices <- sort(unique(as.integer(selected_indices)))
  if (!inherits(fit, "fash") || length(selected_indices) == 0L ||
      min(selected_indices) < 1L ||
      max(selected_indices) > length(fit$fash_data$data_list)) {
    stop("fit or selected_indices are invalid.")
  }
  output <- fit
  output$L_matrix <- fit$L_matrix[selected_indices, , drop = FALSE]
  output$posterior_weights <- fit$posterior_weights[
    selected_indices, , drop = FALSE
  ]
  output$lfdr <- as.numeric(fit$lfdr[selected_indices])
  output$fash_data$data_list <- fit$fash_data$data_list[selected_indices]
  output$fash_data$S <- fit$fash_data$S[selected_indices]
  if (!is.null(fit$fash_data$Omega)) {
    output$fash_data$Omega <- fit$fash_data$Omega[selected_indices]
  }
  output$original_indices <- if (!is.null(fit$original_indices)) {
    as.integer(fit$original_indices[selected_indices])
  } else {
    selected_indices
  }
  class(output) <- class(fit)
  output
}

classify_functional_draws <- function(samples,
                                      smooth_var,
                                      switch_threshold = 0.25) {
  samples <- as.matrix(samples)
  smooth_var <- as.numeric(smooth_var)
  if (nrow(samples) != length(smooth_var) || ncol(samples) == 0L) {
    stop("Posterior samples and smooth_var are not aligned.")
  }

  early <- matrixStats::colMaxs(abs(samples[
    smooth_var <= 3, , drop = FALSE
  ])) - matrixStats::colMaxs(abs(samples[
    smooth_var > 3, , drop = FALSE
  ]))
  middle <- matrixStats::colMaxs(abs(samples[
    smooth_var >= 4 & smooth_var <= 11, , drop = FALSE
  ])) - pmax(
    matrixStats::colMaxs(abs(samples[
      smooth_var < 4, , drop = FALSE
    ])),
    matrixStats::colMaxs(abs(samples[
      smooth_var > 11, , drop = FALSE
    ]))
  )
  late <- matrixStats::colMaxs(abs(samples[
    smooth_var >= 12, , drop = FALSE
  ])) - matrixStats::colMaxs(abs(samples[
    smooth_var < 12, , drop = FALSE
  ]))
  positive <- matrixStats::colMaxs(pmax(samples, 0))
  negative <- matrixStats::colMaxs(pmax(-samples, 0))
  switch <- pmin(positive, negative) - switch_threshold

  c(
    early = mean(early <= 0),
    middle = mean(middle <= 0),
    late = mean(late <= 0),
    switch = mean(switch <= 0)
  )
}

add_cumulative_fsr <- function(classification_table) {
  required <- c("original_index", "key", "category", "lfsr")
  if (!all(required %in% names(classification_table))) {
    stop("classification_table is missing required columns.")
  }
  split_table <- split(classification_table, classification_table$category)
  output <- lapply(split_table, function(category_table) {
    ordering <- order(
      category_table$lfsr,
      category_table$original_index,
      category_table$key
    )
    ranked <- category_table[ordering, , drop = FALSE]
    ranked$rank <- seq_len(nrow(ranked))
    ranked$cfsr <- cumsum(ranked$lfsr) / ranked$rank
    ranked
  })
  result <- do.call(rbind, output)
  rownames(result) <- NULL
  result
}

select_distinct_gene_examples <- function(candidates,
                                          count,
                                          rank_columns,
                                          decreasing = FALSE,
                                          excluded_keys = character()) {
  required <- c("key", "gene_id", rank_columns)
  if (!all(required %in% names(candidates))) {
    stop("Replacement candidates are missing required columns.")
  }
  candidates <- candidates[
    !candidates$key %in% excluded_keys,
    ,
    drop = FALSE
  ]
  if (length(decreasing) == 1L) {
    decreasing <- rep(decreasing, length(rank_columns))
  }
  if (length(decreasing) != length(rank_columns)) {
    stop("decreasing must align with rank_columns.")
  }
  rank_values <- lapply(seq_along(rank_columns), function(index) {
    value <- candidates[[rank_columns[index]]]
    if (decreasing[index]) -value else value
  })
  ordering <- do.call(order, c(rank_values, list(candidates$key)))
  candidates <- candidates[ordering, , drop = FALSE]
  candidates <- candidates[!duplicated(candidates$gene_id), , drop = FALSE]
  utils::head(candidates, count)
}

three_set_venn_regions <- function(first, second, third, labels) {
  if (length(labels) != 3L) {
    stop("Exactly three labels are required.")
  }
  first <- unique(as.character(first))
  second <- unique(as.character(second))
  third <- unique(as.character(third))
  universe <- sort(unique(c(first, second, third)))
  membership <- data.frame(
    item = universe,
    first = universe %in% first,
    second = universe %in% second,
    third = universe %in% third,
    stringsAsFactors = FALSE
  )
  region_code <- paste0(
    as.integer(membership$first),
    as.integer(membership$second),
    as.integer(membership$third)
  )
  codes <- c("100", "010", "001", "110", "101", "011", "111")
  region_labels <- c(
    labels[1L], labels[2L], labels[3L],
    paste(labels[1:2], collapse = " + "),
    paste(labels[c(1, 3)], collapse = " + "),
    paste(labels[2:3], collapse = " + "),
    paste(labels, collapse = " + ")
  )
  data.frame(
    region_code = codes,
    region = region_labels,
    count = vapply(codes, function(code) sum(region_code == code), integer(1)),
    stringsAsFactors = FALSE
  )
}

extract_observed_data <- function(fit, local_index, method, order_label) {
  data <- fit$fash_data$data_list[[local_index]]
  standard_error <- fit$fash_data$S[[local_index]]
  key <- names(fit$fash_data$data_list)[local_index]
  parsed <- parse_pair_keys(key)
  data.frame(
    method = method,
    order = order_label,
    key = key,
    gene_id = parsed$gene_id,
    variant_id = parsed$variant_id,
    time = as.numeric(data$x),
    beta = as.numeric(data$y),
    standard_error = as.numeric(standard_error),
    stringsAsFactors = FALSE
  )
}

extract_posterior_plot_data <- function(fit,
                                        local_indices,
                                        method,
                                        order_label,
                                        smooth_var = seq(0, 15, by = 0.1),
                                        sample_size = 3000L,
                                        seed = 20260810L) {
  local_indices <- as.integer(local_indices)
  output <- vector("list", length(local_indices))
  for (position in seq_along(local_indices)) {
    local_index <- local_indices[position]
    set.seed(seed + fit$original_indices[local_index])
    samples <- predict(
      fit,
      index = local_index,
      smooth_var = smooth_var,
      only.samples = TRUE,
      M = sample_size
    )
    quantiles <- matrixStats::rowQuantiles(
      samples,
      probs = c(0.025, 0.975)
    )
    key <- names(fit$fash_data$data_list)[local_index]
    parsed <- parse_pair_keys(key)
    output[[position]] <- data.frame(
      method = method,
      order = order_label,
      key = key,
      gene_id = parsed$gene_id,
      variant_id = parsed$variant_id,
      original_index = fit$original_indices[local_index],
      time = smooth_var,
      posterior_mean = rowMeans(samples),
      lower = quantiles[, 1L],
      upper = quantiles[, 2L],
      lfdr = fit$lfdr[local_index],
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, output)
}

fit_parametric_curves <- function(observed_table,
                                  grid = seq(0, 15, length.out = 151L)) {
  required <- c("method", "order", "key", "time", "beta", "standard_error")
  if (!all(required %in% names(observed_table))) {
    stop("observed_table is missing required columns.")
  }
  split_data <- split(
    observed_table,
    interaction(observed_table$method, observed_table$order,
                observed_table$key, drop = TRUE)
  )
  result <- lapply(split_data, function(data) {
    weights <- 1 / data$standard_error^2
    linear_fit <- stats::lm(beta ~ time, data = data, weights = weights)
    quadratic_fit <- stats::lm(
      beta ~ time + I(time^2),
      data = data,
      weights = weights
    )
    base <- data[1L, c("method", "order", "key"), drop = FALSE]
    data.frame(
      base[rep(1L, length(grid)), , drop = FALSE],
      time = grid,
      linear = stats::predict(linear_fit, newdata = data.frame(time = grid)),
      quadratic = stats::predict(
        quadratic_fit,
        newdata = data.frame(time = grid)
      ),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}

hallmark_hypergeometric <- function(selected_genes,
                                    universe_genes,
                                    term_to_gene) {
  required <- c("term", "gene")
  if (!all(required %in% names(term_to_gene))) {
    stop("term_to_gene must contain term and gene columns.")
  }
  universe_genes <- unique(as.character(universe_genes))
  selected_genes <- intersect(unique(as.character(selected_genes)), universe_genes)
  term_to_gene <- unique(term_to_gene[
    term_to_gene$gene %in% universe_genes,
    required,
    drop = FALSE
  ])
  term_sets <- split(term_to_gene$gene, term_to_gene$term)
  result <- lapply(names(term_sets), function(term) {
    term_genes <- unique(term_sets[[term]])
    overlap <- intersect(selected_genes, term_genes)
    count <- length(overlap)
    term_size <- length(term_genes)
    selected_size <- length(selected_genes)
    universe_size <- length(universe_genes)
    p_value <- stats::phyper(
      count - 1L,
      term_size,
      universe_size - term_size,
      selected_size,
      lower.tail = FALSE
    )
    data.frame(
      term = term,
      overlap_count = count,
      selected_size = selected_size,
      term_size = term_size,
      universe_size = universe_size,
      p_value = p_value,
      overlap_genes = paste(sort(overlap), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  result <- do.call(rbind, result)
  result$q_value <- stats::p.adjust(result$p_value, method = "BH")
  result[order(result$p_value, result$term), , drop = FALSE]
}

percent_change <- function(current, comparison) {
  ifelse(current == 0, NA_real_, 100 * (comparison - current) / current)
}
