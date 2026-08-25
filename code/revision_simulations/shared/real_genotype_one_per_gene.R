# Shared real-genotype utilities for the formal R1 and R2 simulations.

validate_integer_scalar <- function(x, name, minimum = 1L) {
  if (length(x) != 1L || !is.finite(x) || x != as.integer(x) || x < minimum) {
    stop(name, " must be an integer greater than or equal to ", minimum, ".")
  }
  as.integer(x)
}

parse_tested_pair_ids <- function(pair_ids) {
  if (!is.character(pair_ids) || length(pair_ids) == 0L || anyNA(pair_ids) ||
      any(!nzchar(pair_ids))) {
    stop("pair_ids must be a nonempty character vector without missing values.")
  }
  matched <- regexec("^([^_]+)_(rs[^_]+)$", pair_ids)
  pieces <- regmatches(pair_ids, matched)
  valid <- lengths(pieces) == 3L
  if (!all(valid)) {
    stop("Malformed tested pair identifier: ", pair_ids[which(!valid)[1L]], ".")
  }
  out <- data.frame(
    gene_id = vapply(pieces, `[[`, character(1), 2L),
    variant_id = vapply(pieces, `[[`, character(1), 3L),
    pair_key = pair_ids,
    stringsAsFactors = FALSE
  )
  out <- unique(out)
  rownames(out) <- NULL
  out
}

sample_one_tested_variant_per_gene <- function(pair_ids, seed) {
  seed <- validate_integer_scalar(seed, "seed", minimum = 0L)
  pair_map <- parse_tested_pair_ids(pair_ids)
  gene_levels <- sort(unique(pair_map$gene_id), method = "radix")
  gene_indices <- split(
    seq_len(nrow(pair_map)),
    factor(pair_map$gene_id, levels = gene_levels)
  )
  set.seed(seed)
  selected_indices <- vapply(
    gene_indices,
    function(indices) {
      if (length(indices) == 1L) indices else sample(indices, size = 1L)
    },
    integer(1)
  )
  selection <- pair_map[selected_indices, , drop = FALSE]
  selection <- selection[match(gene_levels, selection$gene_id), , drop = FALSE]
  rownames(selection) <- NULL
  if (nrow(selection) != length(gene_levels) ||
      anyDuplicated(selection$gene_id) ||
      anyDuplicated(selection$pair_key)) {
    stop("The one-per-gene selection is not one unique pair per tested gene.")
  }
  selection
}

read_vcf_sample_ids <- function(vcf_path) {
  if (length(vcf_path) != 1L || !file.exists(vcf_path)) {
    stop("vcf_path does not exist: ", vcf_path)
  }
  connection <- gzfile(vcf_path, open = "rt")
  on.exit(close(connection), add = TRUE)
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (length(line) == 0L) break
    if (startsWith(line, "#CHROM\t")) {
      fields <- strsplit(line, "\t", fixed = TRUE)[[1L]]
      if (length(fields) < 10L) {
        stop("The VCF #CHROM header does not contain dosage samples.")
      }
      sample_ids <- fields[10:length(fields)]
      if (any(!nzchar(sample_ids)) || anyDuplicated(sample_ids)) {
        stop("The VCF sample identifiers must be nonempty and unique.")
      }
      return(sample_ids)
    }
  }
  stop("The VCF does not contain a #CHROM header.")
}

extract_target_vcf_dosages <- function(vcf_path,
                                       target_variants,
                                       work_dir = tempdir()) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("The data.table package is required to read filtered VCF rows.")
  }
  if (length(vcf_path) != 1L || !file.exists(vcf_path)) {
    stop("vcf_path does not exist: ", vcf_path)
  }
  target_variants <- sort(unique(as.character(target_variants)), method = "radix")
  if (length(target_variants) == 0L || anyNA(target_variants) ||
      any(!nzchar(target_variants))) {
    stop("target_variants must contain one or more nonempty identifiers.")
  }
  gzip_path <- Sys.which("gzip")
  awk_path <- Sys.which("awk")
  if (!nzchar(gzip_path) || !nzchar(awk_path)) {
    stop("Both gzip and awk must be available to stream the VCF.")
  }
  dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
  target_path <- tempfile("target-variants-", tmpdir = work_dir, fileext = ".txt")
  filtered_path <- tempfile("target-vcf-", tmpdir = work_dir, fileext = ".tsv")
  on.exit(unlink(c(target_path, filtered_path), force = TRUE), add = TRUE)
  writeLines(target_variants, target_path, useBytes = TRUE)

  awk_program <- paste0(
    "BEGIN { FS=\"\\t\" } ",
    "NR==FNR { keep[$1]=1; next } ",
    "!/^#/ && ($3 in keep) { print }"
  )
  command <- paste(
    shQuote(gzip_path), "-dc", shQuote(vcf_path),
    "|", shQuote(awk_path), shQuote(awk_program), shQuote(target_path), "-",
    ">", shQuote(filtered_path)
  )
  status <- system(command, ignore.stdout = FALSE, ignore.stderr = FALSE)
  if (!identical(status, 0L)) {
    stop("Streaming VCF extraction failed with status ", status, ".")
  }
  if (!file.exists(filtered_path) || file.info(filtered_path)$size == 0) {
    stop("No target variants were recovered from the VCF.")
  }

  sample_ids <- read_vcf_sample_ids(vcf_path)
  raw <- data.table::fread(
    filtered_path,
    sep = "\t",
    header = FALSE,
    colClasses = "character",
    data.table = FALSE,
    showProgress = FALSE
  )
  expected_columns <- 9L + length(sample_ids)
  if (ncol(raw) != expected_columns) {
    stop(
      "Filtered VCF rows have ", ncol(raw),
      " columns; expected ", expected_columns, "."
    )
  }
  variant_ids <- raw[[3L]]
  if (anyDuplicated(variant_ids)) {
    stop(
      "The VCF contains duplicate retained variant ID: ",
      variant_ids[duplicated(variant_ids)][1L], "."
    )
  }
  if (any(grepl(",", raw[[5L]], fixed = TRUE))) {
    stop("The retained VCF rows must be biallelic.")
  }
  if (any(raw[[9L]] != "DS")) {
    stop("Every retained VCF row must use FORMAT=DS.")
  }

  dosage <- vapply(
    raw[10:expected_columns],
    function(column) suppressWarnings(as.numeric(column)),
    numeric(nrow(raw))
  )
  if (is.null(dim(dosage))) {
    dosage <- matrix(dosage, nrow = nrow(raw))
  }
  rownames(dosage) <- variant_ids
  colnames(dosage) <- sample_ids
  if (any(!is.finite(dosage))) {
    stop("The retained VCF dosage matrix contains missing or non-finite values.")
  }

  metadata <- data.frame(
    chromosome = raw[[1L]],
    position = suppressWarnings(as.integer(raw[[2L]])),
    variant_id = variant_ids,
    ref = raw[[4L]],
    alt = raw[[5L]],
    stringsAsFactors = FALSE
  )
  if (anyNA(metadata$position) || any(metadata$position < 1L)) {
    stop("The retained VCF positions must be positive integers.")
  }

  list(
    sample_ids = sample_ids,
    metadata = metadata,
    dosage = dosage,
    requested_variant_count = length(target_variants),
    matched_variant_count = nrow(metadata)
  )
}

assemble_one_per_gene_genotype_sample <- function(selection,
                                                  extracted,
                                                  seed) {
  seed <- validate_integer_scalar(seed, "seed", minimum = 0L)
  required_selection <- c("gene_id", "variant_id", "pair_key")
  if (!is.data.frame(selection) ||
      !all(required_selection %in% names(selection)) ||
      nrow(selection) == 0L ||
      anyNA(selection[required_selection]) ||
      anyDuplicated(selection$gene_id) ||
      anyDuplicated(selection$pair_key)) {
    stop("selection must contain one unique pair_key per unique gene_id.")
  }
  required_extracted <- c("sample_ids", "metadata", "dosage")
  if (!is.list(extracted) || !all(required_extracted %in% names(extracted))) {
    stop("extracted must contain sample_ids, metadata, and dosage.")
  }
  dosage <- as.matrix(extracted$dosage)
  storage.mode(dosage) <- "numeric"
  variant_rows <- match(selection$variant_id, rownames(dosage))
  metadata_rows <- match(selection$variant_id, extracted$metadata$variant_id)
  if (anyNA(variant_rows) || anyNA(metadata_rows)) {
    missing_ids <- unique(selection$variant_id[is.na(variant_rows) | is.na(metadata_rows)])
    stop("Selected variants are missing from extracted VCF rows: ", missing_ids[1L], ".")
  }
  G <- t(dosage[variant_rows, , drop = FALSE])
  rownames(G) <- extracted$sample_ids
  colnames(G) <- selection$pair_key
  allele_frequency <- colMeans(G) / 2
  observed_maf <- pmin(allele_frequency, 1 - allele_frequency)
  genotype_sd <- apply(G, 2L, stats::sd)
  metadata <- extracted$metadata[metadata_rows, , drop = FALSE]
  variant_info <- data.frame(
    unit_index = seq_len(nrow(selection)),
    unit_id = selection$pair_key,
    gene_id = selection$gene_id,
    variant_id = selection$variant_id,
    chromosome = metadata$chromosome,
    position = metadata$position,
    ref = metadata$ref,
    alt = metadata$alt,
    observed_maf = observed_maf,
    genotype_sd = genotype_sd,
    stringsAsFactors = FALSE
  )
  repeated_assignments <- nrow(selection) - length(unique(selection$variant_id))
  selection_summary <- data.frame(
    seed = seed,
    genes = nrow(selection),
    selected_pairs = nrow(selection),
    unique_variant_ids = length(unique(selection$variant_id)),
    repeated_cross_gene_assignments = repeated_assignments,
    maf_min = min(observed_maf),
    maf_q25 = unname(stats::quantile(observed_maf, 0.25)),
    maf_median = stats::median(observed_maf),
    maf_mean = mean(observed_maf),
    maf_q75 = unname(stats::quantile(observed_maf, 0.75)),
    maf_max = max(observed_maf),
    stringsAsFactors = FALSE
  )
  list(
    seed = seed,
    G = G,
    variant_info = variant_info,
    selection = selection,
    selection_summary = selection_summary
  )
}

validate_real_genotype_sample <- function(sample,
                                          expected_genes = 6362L,
                                          expected_donors = 19L,
                                          maf_min = 0.10) {
  expected_genes <- validate_integer_scalar(expected_genes, "expected_genes")
  expected_donors <- validate_integer_scalar(expected_donors, "expected_donors")
  if (length(maf_min) != 1L || !is.finite(maf_min) ||
      maf_min < 0 || maf_min >= 0.5) {
    stop("maf_min must be finite and in [0, 0.5).")
  }
  required <- c("seed", "G", "variant_info", "selection", "selection_summary")
  if (!is.list(sample) || !all(required %in% names(sample))) {
    stop("sample is missing required real-genotype fields.")
  }
  G <- as.matrix(sample$G)
  info <- sample$variant_info
  selection <- sample$selection
  if (!identical(dim(G), c(expected_donors, expected_genes)) ||
      nrow(info) != expected_genes ||
      nrow(selection) != expected_genes) {
    stop("The real-genotype sample dimensions do not match the expected design.")
  }
  if (any(!is.finite(G)) ||
      anyDuplicated(rownames(G)) ||
      anyDuplicated(colnames(G)) ||
      anyDuplicated(info$gene_id) ||
      anyDuplicated(selection$gene_id) ||
      !identical(colnames(G), info$unit_id) ||
      !identical(info$unit_id, selection$pair_key) ||
      !identical(info$gene_id, selection$gene_id) ||
      !identical(info$variant_id, selection$variant_id)) {
    stop("The real-genotype sample identifiers are incomplete or misaligned.")
  }
  observed_maf <- pmin(colMeans(G) / 2, 1 - colMeans(G) / 2)
  genotype_sd <- apply(G, 2L, stats::sd)
  if (any(genotype_sd <= 0) ||
      any(observed_maf < maf_min - 1e-12) ||
      max(abs(observed_maf - info$observed_maf)) > 1e-12 ||
      max(abs(genotype_sd - info$genotype_sd)) > 1e-12) {
    stop("The real-genotype sample fails MAF, polymorphism, or metadata validation.")
  }
  sample$G <- G
  sample
}

artifact_fingerprint <- function(path) {
  if (length(path) != 1L || !file.exists(path)) {
    stop("Cannot fingerprint missing file: ", path)
  }
  information <- file.info(path)
  checksum <- unname(tools::md5sum(path))
  if (is.na(checksum)) stop("Could not compute MD5 for ", path, ".")
  list(
    file_name = basename(path),
    size_bytes = unname(information$size),
    modification_time_utc = format(information$mtime, tz = "UTC", usetz = TRUE),
    md5 = checksum
  )
}

object_md5 <- function(x) {
  path <- tempfile("object-md5-", fileext = ".rds")
  on.exit(unlink(path, force = TRUE), add = TRUE)
  saveRDS(x, path, version = 2)
  checksum <- unname(tools::md5sum(path))
  if (is.na(checksum)) stop("Could not compute an object MD5 checksum.")
  checksum
}

.write_length_prefixed_utf8 <- function(connection, values) {
  values <- enc2utf8(values)
  writeBin(
    as.integer(length(values)),
    connection,
    size = 4L,
    endian = "big"
  )
  for (value in values) {
    bytes <- charToRaw(value)
    writeBin(
      as.integer(length(bytes)),
      connection,
      size = 4L,
      endian = "big"
    )
    writeBin(bytes, connection)
  }
  invisible(NULL)
}

genotype_content_md5 <- function(pair_key, sample_ids, G) {
  if (!is.character(pair_key) || length(pair_key) == 0L ||
      anyNA(pair_key) || any(!nzchar(pair_key)) || anyDuplicated(pair_key)) {
    stop("pair_key must contain unique, nonempty character identifiers.")
  }
  if (!is.character(sample_ids) || length(sample_ids) == 0L ||
      anyNA(sample_ids) || any(!nzchar(sample_ids)) ||
      anyDuplicated(sample_ids)) {
    stop("sample_ids must contain unique, nonempty character identifiers.")
  }
  if (!is.matrix(G) || !is.numeric(G) || length(dim(G)) != 2L ||
      nrow(G) == 0L || ncol(G) == 0L || any(!is.finite(G))) {
    stop("G must be a nonempty finite numeric matrix.")
  }
  if (length(pair_key) != ncol(G) || length(sample_ids) != nrow(G)) {
    stop("pair_key and sample_ids must align with the columns and rows of G.")
  }

  path <- tempfile("genotype-content-md5-", fileext = ".bin")
  connection <- file(path, open = "wb")
  connection_open <- TRUE
  on.exit({
    if (connection_open) close(connection)
    unlink(path, force = TRUE)
  }, add = TRUE)

  header <- charToRaw("fash-genotype-content-md5-v1")
  writeBin(
    as.integer(length(header)),
    connection,
    size = 4L,
    endian = "big"
  )
  writeBin(header, connection)
  .write_length_prefixed_utf8(connection, pair_key)
  .write_length_prefixed_utf8(connection, sample_ids)
  writeBin(as.integer(nrow(G)), connection, size = 4L, endian = "big")
  writeBin(as.integer(ncol(G)), connection, size = 4L, endian = "big")
  writeBin(as.double(G), connection, size = 8L, endian = "big")
  close(connection)
  connection_open <- FALSE

  checksum <- unname(tools::md5sum(path))
  if (is.na(checksum)) {
    stop("Could not compute a canonical genotype-content MD5 checksum.")
  }
  checksum
}

make_maf_balanced_class_targets <- function(maf,
                                            class_probs,
                                            seed,
                                            n_strata = 10L) {
  maf <- as.numeric(maf)
  seed <- validate_integer_scalar(seed, "seed", minimum = 0L)
  n_strata <- validate_integer_scalar(n_strata, "n_strata")
  if (length(maf) == 0L || any(!is.finite(maf)) ||
      any(maf < 0) || any(maf > 0.5) ||
      is.null(names(class_probs)) || any(!nzchar(names(class_probs))) ||
      any(!is.finite(class_probs)) || any(class_probs < 0) ||
      abs(sum(class_probs) - 1) > 1e-8) {
    stop("Invalid MAF values or class probabilities.")
  }
  counts <- exact_proportional_counts(length(maf), class_probs)
  slots <- do.call(rbind, lapply(names(counts), function(effect_class) {
    count <- counts[[effect_class]]
    if (count == 0L) return(NULL)
    data.frame(
      effect_class = rep(effect_class, count),
      position = (seq_len(count) - 0.5) / count,
      stringsAsFactors = FALSE
    )
  }))
  slots <- slots[order(slots$position, slots$effect_class), , drop = FALSE]
  maf_order <- order(maf, seq_along(maf))
  targets <- character(length(maf))
  targets[maf_order] <- slots$effect_class

  rank_index <- integer(length(maf))
  rank_index[maf_order] <- seq_along(maf)
  stratum <- pmin(
    n_strata,
    pmax(1L, ceiling(rank_index * n_strata / length(maf)))
  )
  set.seed(seed)
  for (current_stratum in sort(unique(stratum))) {
    indices <- which(stratum == current_stratum)
    targets[indices] <- sample(targets[indices], length(indices), replace = FALSE)
  }
  observed_counts <- table(factor(targets, levels = names(counts)))
  if (!identical(as.integer(observed_counts), as.integer(counts))) {
    stop("MAF-balanced class assignment changed the requested global counts.")
  }
  targets
}

reassign_effect_simulation_by_maf <- function(effect_sim,
                                               maf,
                                               class_probs,
                                               seed,
                                               n_strata = 10L) {
  if (!is.list(effect_sim) ||
      is.null(effect_sim$beta_matrix) ||
      is.null(effect_sim$unit_info) ||
      !"effect_class" %in% names(effect_sim$unit_info)) {
    stop("effect_sim must contain beta_matrix and unit_info$effect_class.")
  }
  J <- nrow(effect_sim$unit_info)
  if (length(maf) != J || nrow(effect_sim$beta_matrix) != J) {
    stop("maf and effect_sim must contain the same number of units.")
  }
  targets <- make_maf_balanced_class_targets(
    maf = maf,
    class_probs = class_probs,
    seed = seed,
    n_strata = n_strata
  )
  source_classes <- effect_sim$unit_info$effect_class
  if (!identical(
    as.integer(table(factor(source_classes, levels = names(class_probs)))),
    as.integer(table(factor(targets, levels = names(class_probs))))
  )) {
    stop("The generated truth counts do not match the MAF-balanced targets.")
  }

  set.seed(seed)
  source_order <- integer(J)
  for (effect_class in names(class_probs)) {
    source_indices <- which(source_classes == effect_class)
    target_indices <- which(targets == effect_class)
    source_order[target_indices] <- sample(
      source_indices,
      length(source_indices),
      replace = FALSE
    )
  }
  if (any(source_order < 1L) || anyDuplicated(source_order)) {
    stop("Truth reassignment did not construct a complete row permutation.")
  }

  row_aligned_names <- c("beta_matrix", "beta_evaluation", "true_functionals")
  for (field in intersect(row_aligned_names, names(effect_sim))) {
    object <- effect_sim[[field]]
    if (!is.null(dim(object)) && nrow(object) == J) {
      effect_sim[[field]] <- object[source_order, , drop = FALSE]
    }
  }
  effect_sim$unit_info <- effect_sim$unit_info[source_order, , drop = FALSE]
  effect_sim$unit_info$unit_index <- seq_len(J)
  effect_sim$unit_info$unit_id <- sprintf(
    "%s_%04d",
    effect_sim$unit_info$effect_class,
    seq_len(J)
  )
  effect_sim$unit_info$variant_id <- sprintf("variant_%04d", seq_len(J))
  rownames(effect_sim$unit_info) <- NULL
  for (field in intersect(row_aligned_names, names(effect_sim))) {
    object <- effect_sim[[field]]
    if (!is.null(dim(object)) && nrow(object) == J) {
      rownames(object) <- effect_sim$unit_info$variant_id
      effect_sim[[field]] <- object
    }
  }
  effect_sim$settings$maf_balancing <- list(
    method = "exact global class counts with within-MAF-stratum permutation",
    n_strata = as.integer(n_strata),
    seed = as.integer(seed)
  )
  effect_sim
}

summarize_truth_maf_balance <- function(variant_info, unit_info, seed) {
  seed <- validate_integer_scalar(seed, "seed", minimum = 0L)
  if (!is.data.frame(variant_info) ||
      !"observed_maf" %in% names(variant_info) ||
      !is.data.frame(unit_info) ||
      !"effect_class" %in% names(unit_info) ||
      nrow(variant_info) != nrow(unit_info)) {
    stop("variant_info and unit_info must contain aligned MAF and class rows.")
  }
  maf <- as.numeric(variant_info$observed_maf)
  if (any(!is.finite(maf))) stop("Observed MAF values must be finite.")
  overall_mean <- mean(maf)
  overall_sd <- stats::sd(maf)
  classes <- unique(unit_info$effect_class)
  rows <- lapply(classes, function(effect_class) {
    values <- maf[unit_info$effect_class == effect_class]
    data.frame(
      seed = seed,
      effect_class = effect_class,
      n = length(values),
      maf_min = min(values),
      maf_q25 = unname(stats::quantile(values, 0.25)),
      maf_median = stats::median(values),
      maf_mean = mean(values),
      maf_q75 = unname(stats::quantile(values, 0.75)),
      maf_max = max(values),
      standardized_mean_difference = if (overall_sd > 0) {
        (mean(values) - overall_mean) / overall_sd
      } else {
        0
      },
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
