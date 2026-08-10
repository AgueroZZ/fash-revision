# Deterministic sampling of paper-derived real-genotype locus blocks.

revision_real_genotype_seed <- function(seed,
                                        component = c(
                                          "gene_priority",
                                          "block_start",
                                          "across_locus"
                                        )) {
  component <- match.arg(component)
  if (length(seed) != 1L || !is.finite(seed) || seed != as.integer(seed)) {
    stop("seed must be a finite integer scalar.")
  }
  offsets <- c(
    gene_priority = 7101L,
    block_start = 7201L,
    across_locus = 7401L
  )
  max_seed <- .Machine$integer.max - 1
  as.integer(((as.double(seed) + offsets[[component]] - 1) %% max_seed) + 1)
}

validate_positive_integer <- function(x, name) {
  if (length(x) != 1L || !is.finite(x) || x < 1 || x != as.integer(x)) {
    stop(name, " must be a positive integer.")
  }
  as.integer(x)
}

parse_paper_pair_ids <- function(pair_ids) {
  if (!is.character(pair_ids) || length(pair_ids) == 0L || anyNA(pair_ids) ||
      any(!nzchar(pair_ids))) {
    stop("pair_ids must be a nonempty character vector without missing or blank values.")
  }
  matched <- regexec("^([^_]+)_(rs[^_]+)$", pair_ids)
  pieces <- regmatches(pair_ids, matched)
  valid <- lengths(pieces) == 3L
  if (!all(valid)) {
    stop("Malformed paper pair identifier: ", pair_ids[which(!valid)[1L]], ".")
  }
  out <- data.frame(
    gene_id = vapply(pieces, `[[`, character(1), 2L),
    variant_id = vapply(pieces, `[[`, character(1), 3L),
    stringsAsFactors = FALSE
  )
  unique(out)
}

read_vcf_sample_ids <- function(vcf_path) {
  if (length(vcf_path) != 1L || !file.exists(vcf_path)) {
    stop("vcf_path does not exist: ", vcf_path)
  }
  connection <- gzfile(vcf_path, open = "rt")
  on.exit(close(connection), add = TRUE)
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (length(line) == 0L) {
      break
    }
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
  target_variants <- unique(as.character(target_variants))
  if (length(target_variants) == 0L || anyNA(target_variants) ||
      any(!nzchar(target_variants))) {
    stop("target_variants must contain one or more nonempty variant identifiers.")
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
    duplicate_id <- variant_ids[duplicated(variant_ids)][1L]
    stop("The VCF contains duplicate retained variant ID: ", duplicate_id, ".")
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
  colnames(dosage) <- sample_ids
  rownames(dosage) <- variant_ids
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

prepare_real_genotype_source <- function(pair_ids,
                                         vcf_path,
                                         seed_list,
                                         n_loci = 20L,
                                         variants_per_locus = 50L,
                                         candidate_genes_per_seed = 120L,
                                         minimum_raw_variants = 80L,
                                         maf_min = 0.10,
                                         work_dir = tempdir()) {
  n_loci <- validate_positive_integer(n_loci, "n_loci")
  variants_per_locus <- validate_positive_integer(
    variants_per_locus,
    "variants_per_locus"
  )
  candidate_genes_per_seed <- validate_positive_integer(
    candidate_genes_per_seed,
    "candidate_genes_per_seed"
  )
  minimum_raw_variants <- validate_positive_integer(
    minimum_raw_variants,
    "minimum_raw_variants"
  )
  seed_list <- as.integer(seed_list)
  if (length(seed_list) == 0L || anyNA(seed_list) || anyDuplicated(seed_list)) {
    stop("seed_list must contain one or more unique integer seeds.")
  }
  if (length(maf_min) != 1L || !is.finite(maf_min) ||
      maf_min <= 0 || maf_min >= 0.5) {
    stop("maf_min must be a finite scalar strictly between 0 and 0.5.")
  }

  pair_map <- parse_paper_pair_ids(pair_ids)
  raw_gene_counts <- table(pair_map$gene_id)
  eligible_gene_names <- names(raw_gene_counts)[
    raw_gene_counts >= minimum_raw_variants
  ]
  if (length(eligible_gene_names) < n_loci) {
    stop(
      "Only ", length(eligible_gene_names),
      " genes meet the raw variant-count floor; ", n_loci, " loci are required."
    )
  }
  candidate_genes_per_seed <- min(
    candidate_genes_per_seed,
    length(eligible_gene_names)
  )

  candidate_orders <- lapply(seed_list, function(seed) {
    set.seed(revision_real_genotype_seed(seed, "gene_priority"))
    head(
      sample(eligible_gene_names, length(eligible_gene_names), replace = FALSE),
      candidate_genes_per_seed
    )
  })
  names(candidate_orders) <- as.character(seed_list)
  candidate_genes <- unique(unlist(candidate_orders, use.names = FALSE))
  candidate_pair_map <- pair_map[
    pair_map$gene_id %in% candidate_genes,
    ,
    drop = FALSE
  ]
  target_variants <- unique(candidate_pair_map$variant_id)
  extracted <- extract_target_vcf_dosages(
    vcf_path = vcf_path,
    target_variants = target_variants,
    work_dir = work_dir
  )

  allele_frequency <- rowMeans(extracted$dosage) / 2
  observed_maf <- pmin(allele_frequency, 1 - allele_frequency)
  polymorphic <- apply(extracted$dosage, 1L, stats::sd) > 0
  keep <- is.finite(observed_maf) &
    observed_maf >= maf_min &
    observed_maf <= 0.5 &
    polymorphic
  variant_table <- extracted$metadata[keep, , drop = FALSE]
  variant_table$observed_maf <- observed_maf[keep]
  dosage <- extracted$dosage[keep, , drop = FALSE]
  if (nrow(variant_table) == 0L) {
    stop("No candidate VCF variants pass the MAF and polymorphism filters.")
  }
  variant_lookup <- setNames(seq_len(nrow(variant_table)), variant_table$variant_id)
  candidate_pair_map$genotype_row <- unname(
    variant_lookup[candidate_pair_map$variant_id]
  )
  candidate_pair_map <- candidate_pair_map[
    !is.na(candidate_pair_map$genotype_row),
    ,
    drop = FALSE
  ]

  eligible_counts <- table(candidate_pair_map$gene_id)
  for (seed in seed_list) {
    priority <- candidate_orders[[as.character(seed)]]
    available <- priority[
      unname(eligible_counts[priority]) >= variants_per_locus &
        !is.na(unname(eligible_counts[priority]))
    ]
    if (length(available) < n_loci) {
      stop(
        "Seed ", seed, " has only ", length(available),
        " eligible candidate loci after VCF and MAF filtering; ",
        n_loci, " are required."
      )
    }
  }

  out <- list(
    settings = list(
      seed_list = seed_list,
      n_loci = n_loci,
      variants_per_locus = variants_per_locus,
      candidate_genes_per_seed = candidate_genes_per_seed,
      minimum_raw_variants = minimum_raw_variants,
      maf_min = maf_min
    ),
    sample_ids = extracted$sample_ids,
    candidate_orders = candidate_orders,
    pair_map = candidate_pair_map,
    variant_table = variant_table,
    dosage = dosage,
    extraction_summary = data.frame(
      requested_variants = extracted$requested_variant_count,
      matched_vcf_variants = extracted$matched_variant_count,
      retained_common_polymorphic_variants = nrow(variant_table),
      candidate_genes = length(candidate_genes),
      stringsAsFactors = FALSE
    )
  )
  class(out) <- c("r1_real_genotype_source", "list")
  out
}

sample_real_genotype_blocks <- function(source, seed) {
  if (!inherits(source, "r1_real_genotype_source")) {
    stop("source must be created by prepare_real_genotype_source().")
  }
  if (length(seed) != 1L || !is.finite(seed) || seed != as.integer(seed)) {
    stop("seed must be a finite integer scalar.")
  }
  seed <- as.integer(seed)
  priority <- source$candidate_orders[[as.character(seed)]]
  if (is.null(priority)) {
    stop("The requested seed was not prepared in the real-genotype source.")
  }
  n_loci <- source$settings$n_loci
  variants_per_locus <- source$settings$variants_per_locus
  eligible_counts <- table(source$pair_map$gene_id)
  available <- priority[
    unname(eligible_counts[priority]) >= variants_per_locus &
      !is.na(unname(eligible_counts[priority]))
  ]
  if (length(available) < n_loci) {
    stop(
      "The prepared source does not contain enough eligible loci for seed ",
      seed, "."
    )
  }
  selected_genes <- head(available, n_loci)
  variant_lookup <- setNames(
    seq_len(nrow(source$variant_table)),
    source$variant_table$variant_id
  )

  set.seed(revision_real_genotype_seed(seed, "block_start"))
  genotype_blocks <- vector("list", n_loci)
  metadata_blocks <- vector("list", n_loci)
  for (block_index in seq_len(n_loci)) {
    gene_id <- selected_genes[[block_index]]
    gene_rows <- source$pair_map[
      source$pair_map$gene_id == gene_id,
      "genotype_row",
      drop = TRUE
    ]
    gene_rows <- unique(gene_rows)
    metadata <- source$variant_table[gene_rows, , drop = FALSE]
    chromosome_order <- suppressWarnings(as.numeric(metadata$chromosome))
    if (anyNA(chromosome_order)) {
      chromosome_order <- match(metadata$chromosome, sort(unique(metadata$chromosome)))
    }
    metadata <- metadata[
      order(chromosome_order, metadata$position, metadata$variant_id),
      ,
      drop = FALSE
    ]
    if (length(unique(metadata$chromosome)) != 1L) {
      stop("Paper variants for gene ", gene_id, " span multiple chromosomes.")
    }
    maximum_start <- nrow(metadata) - variants_per_locus + 1L
    if (maximum_start < 1L) {
      stop("Gene ", gene_id, " does not contain a complete variant block.")
    }
    block_start <- sample.int(maximum_start, 1L)
    block_rows <- block_start:(block_start + variants_per_locus - 1L)
    block_metadata <- metadata[block_rows, , drop = FALSE]
    source_rows <- unname(variant_lookup[block_metadata$variant_id])
    block_genotype <- t(source$dosage[source_rows, , drop = FALSE])
    unit_ids <- paste(gene_id, block_metadata$variant_id, sep = "_")
    rownames(block_genotype) <- source$sample_ids
    colnames(block_genotype) <- unit_ids
    genotype_blocks[[block_index]] <- block_genotype

    block_metadata$unit_index <- NA_integer_
    block_metadata$unit_id <- unit_ids
    block_metadata$block_index <- block_index
    block_metadata$gene_id <- gene_id
    block_metadata$block_start_position <- min(block_metadata$position)
    block_metadata$block_end_position <- max(block_metadata$position)
    block_metadata$sampling_seed <- seed
    metadata_blocks[[block_index]] <- block_metadata[, c(
      "unit_index", "unit_id", "block_index", "gene_id", "variant_id",
      "chromosome", "position", "ref", "alt", "observed_maf",
      "block_start_position", "block_end_position", "sampling_seed"
    )]
  }

  G <- do.call(cbind, genotype_blocks)
  variant_info <- do.call(rbind, metadata_blocks)
  rownames(variant_info) <- NULL
  variant_info$unit_index <- seq_len(nrow(variant_info))
  expected_dimension <- c(
    length(source$sample_ids),
    n_loci * variants_per_locus
  )
  if (!identical(dim(G), expected_dimension) ||
      !identical(colnames(G), variant_info$unit_id) ||
      any(!is.finite(G)) || any(apply(G, 2L, stats::sd) <= 0)) {
    stop("The sampled real-genotype matrix failed its structural validation.")
  }
  if (!all(table(variant_info$block_index) == variants_per_locus)) {
    stop("The sampled real-genotype blocks do not have the requested size.")
  }

  list(
    G = G,
    variant_info = variant_info,
    selected_genes = selected_genes,
    settings = source$settings
  )
}

summarize_real_genotype_ld <- function(genotype, variant_info) {
  genotype <- as.matrix(genotype)
  storage.mode(genotype) <- "numeric"
  required_columns <- c(
    "block_index", "gene_id", "position", "observed_maf", "sampling_seed"
  )
  if (ncol(genotype) != nrow(variant_info) ||
      !all(required_columns %in% colnames(variant_info))) {
    stop("genotype and variant_info are not aligned for LD summarization.")
  }
  if (any(!is.finite(genotype)) || any(apply(genotype, 2L, stats::sd) <= 0)) {
    stop("genotype must be finite and polymorphic for LD summarization.")
  }
  correlation <- stats::cor(genotype)
  if (any(!is.finite(correlation))) {
    stop("The genotype correlation matrix contains non-finite values.")
  }

  block_ids <- unique(variant_info$block_index)
  locus_rows <- vector("list", length(block_ids))
  within_values <- numeric(0)
  for (block_offset in seq_along(block_ids)) {
    block_id <- block_ids[[block_offset]]
    indices <- which(variant_info$block_index == block_id)
    block_correlation <- correlation[indices, indices, drop = FALSE]
    off_diagonal <- block_correlation[upper.tri(block_correlation)]
    adjacent <- if (length(indices) > 1L) {
      block_correlation[cbind(
        seq_len(length(indices) - 1L),
        2:length(indices)
      )]
    } else {
      numeric(0)
    }
    within_values <- c(within_values, off_diagonal)
    locus_rows[[block_offset]] <- data.frame(
      sampling_seed = unique(variant_info$sampling_seed[indices]),
      block_index = block_id,
      gene_id = unique(variant_info$gene_id[indices]),
      chromosome = unique(variant_info$chromosome[indices]),
      start_position = min(variant_info$position[indices]),
      end_position = max(variant_info$position[indices]),
      span_bp = max(variant_info$position[indices]) -
        min(variant_info$position[indices]),
      n_variants = length(indices),
      median_maf = stats::median(variant_info$observed_maf[indices]),
      median_abs_r = stats::median(abs(off_diagonal)),
      fraction_abs_r_ge_0_8 = mean(abs(off_diagonal) >= 0.8),
      mean_adjacent_abs_r = mean(abs(adjacent)),
      stringsAsFactors = FALSE
    )
  }
  locus <- do.call(rbind, locus_rows)
  rownames(locus) <- NULL

  block_vector <- variant_info$block_index
  across_index <- which(
    upper.tri(correlation) & outer(block_vector, block_vector, `!=`),
    arr.ind = TRUE
  )
  if (nrow(across_index) == 0L || length(within_values) == 0L) {
    stop("Both within-locus and across-locus pairs are required for LD summaries.")
  }
  set.seed(revision_real_genotype_seed(
    unique(variant_info$sampling_seed),
    "across_locus"
  ))
  matched_count <- min(nrow(across_index), length(within_values))
  selected_across <- across_index[
    sample.int(nrow(across_index), matched_count, replace = FALSE),
    ,
    drop = FALSE
  ]
  across_values <- correlation[selected_across]
  seed <- data.frame(
    sampling_seed = unique(variant_info$sampling_seed),
    n_loci = length(block_ids),
    n_variants = ncol(genotype),
    median_locus_span_bp = stats::median(locus$span_bp),
    pooled_within_locus_median_abs_r = stats::median(abs(within_values)),
    pooled_within_locus_fraction_abs_r_ge_0_8 = mean(abs(within_values) >= 0.8),
    mean_locus_median_abs_r = mean(locus$median_abs_r),
    mean_locus_adjacent_abs_r = mean(locus$mean_adjacent_abs_r),
    matched_across_locus_median_abs_r = stats::median(abs(across_values)),
    within_minus_across_median_abs_r = stats::median(abs(within_values)) -
      stats::median(abs(across_values)),
    stringsAsFactors = FALSE
  )

  list(
    locus = locus,
    seed = seed,
    correlation = correlation,
    within_abs_r = abs(within_values),
    across_abs_r = abs(across_values)
  )
}
