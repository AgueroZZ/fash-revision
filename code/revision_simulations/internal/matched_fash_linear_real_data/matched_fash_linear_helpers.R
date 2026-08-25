# Helpers for the versioned matched FASH-linear real-data rerun.

find_workflowr_root_matched <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "_workflowr.yml"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the workflowr project root.")
    }
    current <- parent
  }
}

atomic_save_rds_matched <- function(object, path) {
  temporary_path <- paste0(path, ".tmp")
  saveRDS(object, temporary_path)
  if (!file.rename(temporary_path, path)) {
    stop("Could not atomically write ", path, ".")
  }
  invisible(path)
}

load_exact_object_matched <- function(path, expected_name) {
  environment <- new.env(parent = emptyenv())
  loaded_names <- load(path, envir = environment)
  if (!identical(loaded_names, expected_name)) {
    stop("Unexpected object name in ", path, ".")
  }
  environment[[expected_name]]
}

make_file_provenance_matched <- function(paths) {
  if (is.null(names(paths)) || any(names(paths) == "") || anyDuplicated(names(paths))) {
    stop("Every provenance path must have a unique nonempty label.")
  }
  information <- file.info(paths)
  if (anyNA(information$size)) {
    stop("At least one provenance path is missing.")
  }
  data.frame(
    label = names(paths),
    path = normalizePath(paths, winslash = "/", mustWork = TRUE),
    byte_size = unname(information$size),
    md5 = unname(tools::md5sum(paths)),
    modified_at = format(information$mtime, tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}

extract_package_provenance_matched <- function(package = "fashr") {
  description <- utils::packageDescription(package)
  data.frame(
    package = package,
    version = as.character(description$Version),
    remote_ref = as.character(description$RemoteRef),
    remote_sha = as.character(description$RemoteSha),
    built = as.character(description$Built),
    library_path = normalizePath(
      find.package(package),
      winslash = "/",
      mustWork = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

extract_null_weight_matched <- function(fit) {
  expanded <- expand_grid_prior_weights(fit$prior_weights, fit$psd_grid)
  unname(expanded[fit$psd_grid == 0])
}

summarize_discoveries_matched <- function(method,
                                          adjustment,
                                          lfdr,
                                          pair_table,
                                          prior_null_weight,
                                          alpha = 0.05) {
  if (length(lfdr) != nrow(pair_table)) {
    stop("lfdr and pair_table are not aligned.")
  }
  indices <- select_cumulative_lfdr_calls_linear(lfdr, alpha)
  table <- pair_table[indices, , drop = FALSE]
  table$lfdr <- as.numeric(lfdr[indices])
  table <- table[order(table$lfdr, table$key), , drop = FALSE]
  rownames(table) <- NULL
  list(
    indices = indices,
    table = table,
    summary = data.frame(
      method = method,
      adjustment = adjustment,
      pair_count = nrow(table),
      gene_count = length(unique(table$gene_id)),
      variant_count = length(unique(table$variant_id)),
      estimated_pi0 = prior_null_weight,
      alpha = alpha,
      stringsAsFactors = FALSE
    )
  )
}

normalize_named_sets_matched <- function(sets) {
  if (!is.list(sets) || length(sets) < 2L || is.null(names(sets)) ||
      any(names(sets) == "") || anyDuplicated(names(sets))) {
    stop("sets must be a named list with at least two entries.")
  }
  lapply(sets, function(values) {
    values <- unique(as.character(values))
    values[!is.na(values) & nzchar(values)]
  })
}

exclusive_set_regions_matched <- function(sets, unit) {
  sets <- normalize_named_sets_matched(sets)
  set_names <- names(sets)
  patterns <- expand.grid(
    rep(list(c(FALSE, TRUE)), length(sets)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  patterns <- patterns[rowSums(patterns) > 0, , drop = FALSE]
  patterns <- patterns[order(-rowSums(patterns), do.call(order, patterns)), , drop = FALSE]
  universe <- unique(unlist(sets, use.names = FALSE))
  membership <- vapply(sets, function(values) universe %in% values, logical(length(universe)))
  if (is.null(dim(membership))) {
    membership <- matrix(membership, ncol = length(sets))
  }

  rows <- lapply(seq_len(nrow(patterns)), function(index) {
    pattern <- as.logical(patterns[index, ])
    matches <- apply(membership, 1L, function(row) identical(as.logical(row), pattern))
    included_names <- set_names[pattern]
    data.frame(
      unit = unit,
      region = paste(included_names, collapse = " & "),
      included_set_count = sum(pattern),
      count = sum(matches),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

pairwise_overlap_matched <- function(sets, unit) {
  sets <- normalize_named_sets_matched(sets)
  combinations <- utils::combn(names(sets), 2L, simplify = FALSE)
  rows <- lapply(combinations, function(methods) {
    first <- sets[[methods[[1L]]]]
    second <- sets[[methods[[2L]]]]
    intersection_count <- length(intersect(first, second))
    union_count <- length(union(first, second))
    data.frame(
      unit = unit,
      method_1 = methods[[1L]],
      method_2 = methods[[2L]],
      method_1_count = length(first),
      method_2_count = length(second),
      intersection_count = intersection_count,
      union_count = union_count,
      jaccard = if (union_count == 0L) NA_real_ else {
        intersection_count / union_count
      },
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

build_four_method_venn_sets_matched <- function(strober_quadratic,
                                                strober_linear,
                                                fash_iwp1,
                                                fash_linear) {
  required_strober <- c("key", "ensamble_id", "rs_id")
  required_fash <- c("key", "gene_id", "variant_id")
  if (!all(required_strober %in% names(strober_quadratic)) ||
      !all(required_strober %in% names(strober_linear)) ||
      !all(required_fash %in% names(fash_iwp1)) ||
      !all(required_fash %in% names(fash_linear))) {
    stop("A discovery table is missing required identifiers.")
  }

  list(
    `Gene-variant pairs` = list(
      `Strober quadratic` = unique(strober_quadratic$key),
      `Strober linear` = unique(strober_linear$key),
      `FASH-IWP1 BF` = unique(fash_iwp1$key),
      `FASH-linear BF` = unique(fash_linear$key)
    ),
    Genes = list(
      `Strober quadratic` = unique(strober_quadratic$ensamble_id),
      `Strober linear` = unique(strober_linear$ensamble_id),
      `FASH-IWP1 BF` = unique(fash_iwp1$gene_id),
      `FASH-linear BF` = unique(fash_linear$gene_id)
    ),
    Variants = list(
      `Strober quadratic` = unique(strober_quadratic$rs_id),
      `Strober linear` = unique(strober_linear$rs_id),
      `FASH-IWP1 BF` = unique(fash_iwp1$variant_id),
      `FASH-linear BF` = unique(fash_linear$variant_id)
    )
  )
}

validate_four_method_venn_sets_matched <- function(venn_sets) {
  expected_units <- c("Gene-variant pairs", "Genes", "Variants")
  expected_methods <- c(
    "Strober quadratic", "Strober linear", "FASH-IWP1 BF", "FASH-linear BF"
  )
  if (!identical(names(venn_sets), expected_units)) {
    stop("The Venn units are not in the expected order.")
  }
  for (unit in expected_units) {
    sets <- normalize_named_sets_matched(venn_sets[[unit]])
    if (!identical(names(sets), expected_methods) || any(lengths(sets) == 0L)) {
      stop("The Venn sets are incomplete for ", unit, ".")
    }
  }
  invisible(TRUE)
}
