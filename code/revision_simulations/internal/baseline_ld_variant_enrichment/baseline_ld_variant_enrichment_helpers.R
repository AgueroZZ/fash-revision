# Helper functions for the internal baselineLD v2.2 variant enrichment analysis.

require_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The ", package, " package is required.")
  }
  invisible(TRUE)
}

parse_fash_pair_keys <- function(pair_keys) {
  pair_keys <- as.character(pair_keys)
  if (!length(pair_keys) || anyNA(pair_keys) ||
      any(!grepl("^[^_]+_.+$", pair_keys))) {
    stop("Pair keys must use the gene_variant format.")
  }
  data.frame(
    pair_key = pair_keys,
    gene_id = sub("_.*$", "", pair_keys),
    variant_id = sub("^[^_]+_", "", pair_keys),
    stringsAsFactors = FALSE
  )
}

select_cumulative_lfdr <- function(lfdr, alpha = 0.05) {
  lfdr <- as.numeric(lfdr)
  if (!length(lfdr) || any(!is.finite(lfdr)) ||
      any(lfdr < 0 | lfdr > 1) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("lfdr and alpha must be valid probabilities.")
  }
  ordering <- order(lfdr, method = "radix")
  cumulative_fdr <- cumsum(lfdr[ordering]) / seq_along(ordering)
  accepted <- which(cumulative_fdr <= alpha)
  if (!length(accepted)) {
    return(integer())
  }
  ordering[seq_len(max(accepted))]
}

validate_strober_table <- function(result, pair_keys) {
  required_columns <- c("rs_id", "ensamble_id", "pvalue", "eFDR")
  if (!identical(names(result), required_columns) || !nrow(result) ||
      anyDuplicated(paste(result$ensamble_id, result$rs_id, sep = "_"))) {
    stop("A Strober result table has an unexpected schema.")
  }
  result_keys <- paste(result$ensamble_id, result$rs_id, sep = "_")
  if (!setequal(result_keys, pair_keys)) {
    stop("A Strober result table does not match the tested pair universe.")
  }
  invisible(TRUE)
}

derive_requested_fash_sets <- function(pair_keys,
                                       lfdr,
                                       linear_results,
                                       nonlinear_results,
                                       alpha = 0.05) {
  pair_table <- parse_fash_pair_keys(pair_keys)
  if (length(lfdr) != nrow(pair_table)) {
    stop("FASH pair keys and lfdr values are not aligned.")
  }
  validate_strober_table(linear_results, pair_keys)
  validate_strober_table(nonlinear_results, pair_keys)

  fash_indices <- select_cumulative_lfdr(lfdr, alpha = alpha)
  fash_pairs <- pair_table[fash_indices, , drop = FALSE]
  fash_pairs$score <- as.numeric(lfdr[fash_indices])
  fash_pairs <- fash_pairs[
    order(fash_pairs$score, fash_pairs$variant_id, method = "radix"),
    ,
    drop = FALSE
  ]
  row.names(fash_pairs) <- NULL

  linear_pairs <- linear_results[linear_results$eFDR <= alpha, , drop = FALSE]
  nonlinear_pairs <- nonlinear_results[
    nonlinear_results$eFDR <= alpha,
    ,
    drop = FALSE
  ]
  strober_pair_keys <- union(
    paste(linear_pairs$ensamble_id, linear_pairs$rs_id, sep = "_"),
    paste(nonlinear_pairs$ensamble_id, nonlinear_pairs$rs_id, sep = "_")
  )
  fash_only_pairs <- fash_pairs[
    !fash_pairs$pair_key %in% strober_pair_keys,
    ,
    drop = FALSE
  ]
  lead_pairs <- fash_pairs[!duplicated(fash_pairs$gene_id), , drop = FALSE]

  selected_sets <- list(
    all_fash = unique(fash_pairs$variant_id),
    one_lead_fash_per_gene = unique(lead_pairs$variant_id),
    fash_only_pair_variants = unique(fash_only_pairs$variant_id)
  )
  pair_summary <- data.frame(
    discovery_set = names(selected_sets),
    pair_count = c(nrow(fash_pairs), nrow(lead_pairs), nrow(fash_only_pairs)),
    unique_variant_count = lengths(selected_sets),
    unique_gene_count = c(
      length(unique(fash_pairs$gene_id)),
      length(unique(lead_pairs$gene_id)),
      length(unique(fash_only_pairs$gene_id))
    ),
    stringsAsFactors = FALSE
  )
  variant_multiplicity <- as.data.frame(table(
    fash_only_pairs$variant_id
  ), stringsAsFactors = FALSE)
  names(variant_multiplicity) <- c("variant_id", "fash_only_pair_count")
  variant_multiplicity$fash_only_pair_count <- as.integer(
    variant_multiplicity$fash_only_pair_count
  )

  list(
    selected_sets = selected_sets,
    pair_summary = pair_summary,
    all_fash_pairs = fash_pairs,
    lead_fash_pairs = lead_pairs,
    fash_only_pairs = fash_only_pairs,
    fash_only_variant_multiplicity = variant_multiplicity
  )
}

read_baseline_ld_header <- function(path) {
  if (!file.exists(path)) {
    stop("baselineLD annotation file does not exist: ", path)
  }
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  line <- readLines(connection, n = 1L, warn = FALSE)
  if (!length(line)) {
    stop("baselineLD annotation file is empty: ", path)
  }
  strsplit(line, "\t", fixed = TRUE)[[1L]]
}

classify_baseline_ld_columns <- function(annotation_table) {
  metadata_columns <- c("CHR", "BP", "SNP", "CM")
  if (!all(metadata_columns %in% names(annotation_table))) {
    stop("The baselineLD annotation table lacks required metadata columns.")
  }
  annotation_columns <- setdiff(names(annotation_table), metadata_columns)
  if (!length(annotation_columns)) {
    stop("No baselineLD annotation columns were found.")
  }
  rows <- lapply(annotation_columns, function(annotation) {
    values <- as.numeric(annotation_table[[annotation]])
    observed <- unique(values[is.finite(values)])
    data.frame(
      annotation = annotation,
      annotation_type = if (length(observed) && all(observed %in% c(0, 1))) {
        "binary"
      } else {
        "continuous"
      },
      n_unique = length(observed),
      minimum = if (length(observed)) min(observed) else NA_real_,
      maximum = if (length(observed)) max(observed) else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

read_baseline_ld_chromosome <- function(path,
                                        chromosome,
                                        keep_variants,
                                        annotation_columns) {
  require_namespace("data.table")
  required_keep_columns <- c("variant_id", "chromosome", "position")
  if (!all(required_keep_columns %in% names(keep_variants)) ||
      anyDuplicated(keep_variants$variant_id) || !length(annotation_columns)) {
    stop("Invalid tested variants or baselineLD annotation columns.")
  }
  chromosome <- as.character(chromosome)
  keep <- keep_variants[as.character(keep_variants$chromosome) == chromosome,
                        , drop = FALSE]
  if (!nrow(keep)) {
    return(NULL)
  }
  expected_columns <- c("CHR", "BP", "SNP", "CM", annotation_columns)
  header <- read_baseline_ld_header(path)
  if (!all(expected_columns %in% header)) {
    stop("The baselineLD chromosome file has an unexpected schema: ", path)
  }
  annotation <- data.table::fread(
    path,
    select = c("CHR", "BP", "SNP", annotation_columns),
    showProgress = FALSE,
    data.table = TRUE
  )
  annotation <- annotation[SNP %chin% keep$variant_id]
  if (anyDuplicated(annotation$SNP)) {
    stop("Duplicate retained rsIDs were found in ", basename(path), ".")
  }
  retained_index <- match(annotation$SNP, keep$variant_id)
  if (anyNA(retained_index) || any(as.character(annotation$CHR) != chromosome) ||
      any(as.integer(annotation$BP) != as.integer(keep$position[retained_index]))) {
    stop("baselineLD rsID coordinates do not agree with the study VCF.")
  }
  for (annotation_name in annotation_columns) {
    values <- annotation[[annotation_name]]
    if (anyNA(values) || any(!values %in% c(0, 1))) {
      stop("Non-binary values were found in annotation ", annotation_name, ".")
    }
    annotation[[annotation_name]] <- as.logical(values)
  }
  data.frame(
    variant_id = as.character(annotation$SNP),
    chromosome = chromosome,
    position = as.integer(annotation$BP),
    annotation[, ..annotation_columns],
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

build_baseline_ld_matrix <- function(paths,
                                     variant_table,
                                     binary_annotation_columns) {
  if (is.null(names(paths)) || any(!nzchar(names(paths))) ||
      anyDuplicated(names(paths)) || !length(paths) ||
      any(!file.exists(paths)) || !length(binary_annotation_columns)) {
    stop("Named baselineLD chromosome paths are required.")
  }
  rows <- lapply(names(paths), function(chromosome) {
    read_baseline_ld_chromosome(
      path = paths[[chromosome]],
      chromosome = chromosome,
      keep_variants = variant_table,
      annotation_columns = binary_annotation_columns
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  combined <- as.data.frame(
    data.table::rbindlist(rows, use.names = TRUE),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!nrow(combined) || anyDuplicated(combined$variant_id)) {
    stop("The retained baselineLD annotation matrix is empty or duplicated.")
  }
  for (annotation in binary_annotation_columns) {
    values <- combined[[annotation]]
    if (anyNA(values) || any(!values %in% c(0, 1))) {
      stop("Non-binary values were found in annotation ", annotation, ".")
    }
    combined[[annotation]] <- as.logical(values)
  }
  variant_order <- match(combined$variant_id, variant_table$variant_id)
  combined <- combined[order(variant_order, method = "radix"), , drop = FALSE]
  row.names(combined) <- NULL
  combined
}

summarize_set_coverage <- function(selected_sets,
                                   covered_ids,
                                   universe_ids) {
  if (is.null(names(selected_sets)) || any(!nzchar(names(selected_sets))) ||
      !length(covered_ids) || !length(universe_ids)) {
    stop("Invalid discovery sets or coverage universe.")
  }
  rows <- lapply(names(selected_sets), function(set_name) {
    selected <- unique(as.character(selected_sets[[set_name]]))
    data.frame(
      discovery_set = set_name,
      original_count = length(selected),
      covered_count = sum(selected %in% covered_ids),
      coverage_proportion = mean(selected %in% covered_ids),
      stringsAsFactors = FALSE
    )
  })
  universe_row <- data.frame(
    discovery_set = "tested_variant_universe",
    original_count = length(unique(universe_ids)),
    covered_count = sum(unique(universe_ids) %in% covered_ids),
    coverage_proportion = mean(unique(universe_ids) %in% covered_ids),
    stringsAsFactors = FALSE
  )
  rbind(universe_row, do.call(rbind, rows))
}
