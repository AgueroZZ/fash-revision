# Helper functions for the internal FASH variant annotation enrichment analysis.

require_namespace <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("The ", package, " package is required.")
  }
  invisible(TRUE)
}

safe_overlap_proportion <- function(overlap, total) {
  overlap <- as.numeric(overlap)
  total <- as.numeric(total)
  output_length <- max(length(overlap), length(total))
  overlap <- rep(overlap, length.out = output_length)
  total <- rep(total, length.out = output_length)
  output <- overlap / total
  zero <- is.finite(total) & total > 0 & overlap == 0
  output[zero] <- 0.5 / (total[zero] + 1)
  output
}

normalize_autosome <- function(chromosome) {
  chromosome <- sub("^chr", "", as.character(chromosome), ignore.case = TRUE)
  chromosome[!chromosome %in% as.character(1:22)] <- NA_character_
  chromosome
}

parse_pair_keys <- function(keys) {
  keys <- as.character(keys)
  if (!length(keys) || anyNA(keys) || any(!nzchar(keys)) ||
      any(!grepl("^[^_]+_.+$", keys))) {
    stop("Pair keys must be non-empty strings in gene_variant format.")
  }
  data.frame(
    pair_key = keys,
    gene_id = sub("_.*$", "", keys),
    variant_id = sub("^[^_]+_", "", keys),
    stringsAsFactors = FALSE
  )
}

cumulative_lfdr_calls <- function(lfdr, alpha = 0.05) {
  lfdr <- as.numeric(lfdr)
  alpha <- as.numeric(alpha)
  if (!length(lfdr) || any(!is.finite(lfdr)) || any(lfdr < 0 | lfdr > 1) ||
      length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
    stop("lfdr and alpha must be valid probabilities.")
  }
  ordering <- order(lfdr, method = "radix")
  cumulative_fdr <- cumsum(lfdr[ordering]) / seq_along(ordering)
  accepted <- which(cumulative_fdr <= alpha)
  if (!length(accepted)) {
    return(integer())
  }
  as.integer(ordering[seq_len(max(accepted))])
}

read_vcf_header <- function(vcf_path) {
  connection <- gzfile(vcf_path, open = "rt")
  on.exit(close(connection), add = TRUE)
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (!length(line)) {
      stop("The VCF header line was not found.")
    }
    if (startsWith(line, "#CHROM\t")) {
      return(strsplit(line, "\t", fixed = TRUE)[[1L]])
    }
  }
}

read_vcf_variant_metadata <- function(vcf_path,
                                      keep_ids,
                                      retain_dosage_ids = character()) {
  require_namespace("data.table")
  if (!file.exists(vcf_path)) {
    stop("VCF file does not exist: ", vcf_path)
  }
  keep_ids <- unique(as.character(keep_ids))
  retain_dosage_ids <- unique(as.character(retain_dosage_ids))
  if (!length(keep_ids) || anyNA(keep_ids) || any(!nzchar(keep_ids)) ||
      any(!retain_dosage_ids %in% keep_ids)) {
    stop("keep_ids and retain_dosage_ids are invalid.")
  }

  header <- read_vcf_header(vcf_path)
  if (length(header) < 10L || !identical(header[1:9], c(
      "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
      "FORMAT"
    ))) {
    stop("The VCF header has an unexpected structure.")
  }
  sample_columns <- header[10:length(header)]
  selected_columns <- c("#CHROM", "POS", "ID", "REF", "ALT", "FORMAT",
                        sample_columns)
  vcf <- data.table::fread(
    vcf_path,
    skip = "#CHROM",
    header = TRUE,
    select = selected_columns,
    na.strings = ".",
    showProgress = FALSE,
    data.table = TRUE
  )
  data.table::setnames(
    vcf,
    c("#CHROM", "POS", "ID", "REF", "ALT", "FORMAT"),
    c("chromosome", "position", "variant_id", "reference", "alternate",
      "format")
  )
  vcf <- vcf[variant_id %chin% keep_ids]
  if (!nrow(vcf)) {
    stop("No requested variants were found in the VCF.")
  }
  if (anyDuplicated(vcf$variant_id)) {
    duplicated_ids <- unique(vcf$variant_id[duplicated(vcf$variant_id)])
    stop("Requested VCF IDs are duplicated: ", paste(head(duplicated_ids),
                                                       collapse = ", "))
  }
  if (any(vcf$format != "DS")) {
    stop("The retained VCF records must use dosage-only FORMAT=DS.")
  }

  dosage <- as.matrix(vcf[, ..sample_columns])
  storage.mode(dosage) <- "double"
  observed_counts <- rowSums(!is.na(dosage))
  alternate_allele_frequency <- rowMeans(dosage, na.rm = TRUE) / 2
  alternate_allele_frequency[observed_counts == 0L] <- NA_real_
  minor_allele_frequency <- pmin(
    alternate_allele_frequency,
    1 - alternate_allele_frequency
  )

  metadata <- data.frame(
    variant_id = vcf$variant_id,
    chromosome = normalize_autosome(vcf$chromosome),
    position = as.integer(vcf$position),
    reference = as.character(vcf$reference),
    alternate = as.character(vcf$alternate),
    alternate_allele_frequency = alternate_allele_frequency,
    minor_allele_frequency = minor_allele_frequency,
    n_nonmissing_dosages = as.integer(observed_counts),
    stringsAsFactors = FALSE
  )
  ordering <- match(keep_ids, metadata$variant_id)
  ordering <- ordering[!is.na(ordering)]
  metadata <- metadata[ordering, , drop = FALSE]

  dosage_output <- matrix(
    numeric(),
    nrow = 0L,
    ncol = length(sample_columns),
    dimnames = list(character(), sample_columns)
  )
  if (length(retain_dosage_ids)) {
    dosage_order <- match(retain_dosage_ids, vcf$variant_id)
    dosage_order <- dosage_order[!is.na(dosage_order)]
    dosage_output <- dosage[dosage_order, , drop = FALSE]
    rownames(dosage_output) <- vcf$variant_id[dosage_order]
    colnames(dosage_output) <- sample_columns
  }

  list(metadata = metadata, dosage = dosage_output)
}

intervals_to_granges <- function(intervals) {
  require_namespace("GenomicRanges")
  require_namespace("IRanges")
  required_columns <- c("chromosome", "start", "end", "annotation")
  if (!all(required_columns %in% names(intervals)) || !nrow(intervals)) {
    stop("Annotation intervals have an unexpected structure.")
  }
  chromosome <- normalize_autosome(intervals$chromosome)
  keep <- !is.na(chromosome) & is.finite(intervals$start) &
    is.finite(intervals$end) & intervals$start >= 1L &
    intervals$end >= intervals$start & nzchar(intervals$annotation)
  intervals <- intervals[keep, , drop = FALSE]
  chromosome <- chromosome[keep]
  GenomicRanges::GRanges(
    seqnames = chromosome,
    ranges = IRanges::IRanges(
      start = as.integer(intervals$start),
      end = as.integer(intervals$end)
    ),
    annotation = as.character(intervals$annotation)
  )
}

annotate_variant_overlaps <- function(variants, intervals) {
  require_namespace("GenomicRanges")
  require_namespace("IRanges")
  required_variant_columns <- c("variant_id", "chromosome", "position")
  if (!all(required_variant_columns %in% names(variants)) || !nrow(variants) ||
      anyDuplicated(variants$variant_id)) {
    stop("Variant coordinates have an unexpected structure.")
  }
  annotation_names <- unique(as.character(intervals$annotation))
  if (!length(annotation_names) || anyNA(annotation_names) ||
      any(!nzchar(annotation_names))) {
    stop("No valid annotation names were supplied.")
  }

  variant_ranges <- GenomicRanges::GRanges(
    seqnames = normalize_autosome(variants$chromosome),
    ranges = IRanges::IRanges(
      start = as.integer(variants$position),
      end = as.integer(variants$position)
    )
  )
  annotation_ranges <- intervals_to_granges(intervals)
  overlap_matrix <- matrix(
    FALSE,
    nrow = nrow(variants),
    ncol = length(annotation_names),
    dimnames = list(NULL, annotation_names)
  )
  hits <- GenomicRanges::findOverlaps(
    variant_ranges,
    annotation_ranges,
    ignore.strand = TRUE
  )
  if (length(hits)) {
    annotation_index <- match(
      as.character(S4Vectors::mcols(annotation_ranges)$annotation)[
        S4Vectors::subjectHits(hits)
      ],
      annotation_names
    )
    overlap_matrix[cbind(S4Vectors::queryHits(hits), annotation_index)] <- TRUE
  }
  data.frame(
    variants[, c("variant_id", "chromosome"), drop = FALSE],
    overlap_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

extract_gtf_gene_id <- function(attributes) {
  attributes <- as.character(attributes)
  matched <- grepl('gene_id "[^"]+"', attributes)
  gene_id <- rep(NA_character_, length(attributes))
  gene_id[matched] <- sub(
    '.*gene_id "([^"]+)".*',
    "\\1",
    attributes[matched]
  )
  sub("\\.[0-9]+$", "", gene_id)
}

count_leading_comment_lines <- function(path, comment_prefix = "#") {
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  count <- 0L
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (!length(line) || !startsWith(line, comment_prefix)) {
      break
    }
    count <- count + 1L
  }
  count
}

read_gencode_v19_annotations <- function(path, promoter_width = 2000L) {
  require_namespace("data.table")
  require_namespace("GenomicRanges")
  require_namespace("IRanges")
  gtf <- data.table::fread(
    path,
    sep = "\t",
    header = FALSE,
    skip = count_leading_comment_lines(path),
    select = c(1L, 3L, 4L, 5L, 7L, 9L),
    col.names = c("chromosome", "feature", "start", "end", "strand",
                  "attributes"),
    showProgress = FALSE,
    data.table = FALSE
  )
  gtf$chromosome <- normalize_autosome(gtf$chromosome)
  gtf <- gtf[!is.na(gtf$chromosome) & gtf$feature %in% c("gene", "exon"),
             , drop = FALSE]
  gtf$gene_id <- extract_gtf_gene_id(gtf$attributes)
  if (!nrow(gtf) || anyNA(gtf$gene_id)) {
    stop("GENCODE v19 parsing failed.")
  }

  genes <- gtf[gtf$feature == "gene", , drop = FALSE]
  genes$tss <- ifelse(genes$strand == "+", genes$start, genes$end)
  gene_tss <- genes[, c("gene_id", "chromosome", "tss", "strand")]
  gene_tss <- gene_tss[!duplicated(gene_tss$gene_id), , drop = FALSE]

  gene_ranges <- GenomicRanges::reduce(GenomicRanges::GRanges(
    seqnames = genes$chromosome,
    ranges = IRanges::IRanges(genes$start, genes$end)
  ))
  exons <- gtf[gtf$feature == "exon", , drop = FALSE]
  exon_ranges <- GenomicRanges::reduce(GenomicRanges::GRanges(
    seqnames = exons$chromosome,
    ranges = IRanges::IRanges(exons$start, exons$end)
  ))
  intron_ranges <- GenomicRanges::setdiff(gene_ranges, exon_ranges,
                                           ignore.strand = TRUE)
  promoter_ranges <- GenomicRanges::reduce(GenomicRanges::GRanges(
    seqnames = genes$chromosome,
    ranges = IRanges::IRanges(
      start = pmax(1L, as.integer(genes$tss) - promoter_width),
      end = as.integer(genes$tss) + promoter_width
    )
  ))

  ranges_to_table <- function(ranges, annotation) {
    data.frame(
      chromosome = as.character(GenomicRanges::seqnames(ranges)),
      start = BiocGenerics::start(ranges),
      end = BiocGenerics::end(ranges),
      annotation = annotation,
      stringsAsFactors = FALSE
    )
  }
  intervals <- rbind(
    ranges_to_table(promoter_ranges, "GENCODE promoter +/-2 kb"),
    ranges_to_table(exon_ranges, "GENCODE exon"),
    ranges_to_table(intron_ranges, "GENCODE intron"),
    ranges_to_table(gene_ranges, "GENCODE gene body")
  )
  list(intervals = intervals, gene_tss = gene_tss)
}

read_encode_ccre_annotations <- function(path) {
  require_namespace("data.table")
  ccre <- data.table::fread(
    path,
    sep = "\t",
    header = FALSE,
    select = c(1L, 2L, 3L, 4L, 9L),
    col.names = c("chromosome", "start0", "end", "accession", "rgb"),
    showProgress = FALSE,
    data.table = FALSE
  )
  category <- unname(c(
    "255,0,0" = "ENCODE cCRE promoter-like",
    "255,205,0" = "ENCODE cCRE enhancer-like",
    "0,176,240" = "ENCODE cCRE CTCF-only"
  )[as.character(ccre$rgb)])
  if (anyNA(category)) {
    stop("ENCODE cCRE file contains an unexpected RGB classification.")
  }
  data.frame(
    chromosome = normalize_autosome(ccre$chromosome),
    start = as.integer(ccre$start0) + 1L,
    end = as.integer(ccre$end),
    annotation = category,
    accession = as.character(ccre$accession),
    stringsAsFactors = FALSE
  )
}

collapse_roadmap_state <- function(state) {
  state_number <- suppressWarnings(as.integer(sub("_.*$", "", state)))
  unname(c(
    "1" = "Active TSS",
    "2" = "Active TSS",
    "3" = "Transcribed",
    "4" = "Transcribed",
    "5" = "Transcribed",
    "6" = "Enhancer",
    "7" = "Enhancer",
    "8" = "Heterochromatin/repeats",
    "9" = "Heterochromatin/repeats",
    "10" = "Bivalent/poised",
    "11" = "Bivalent/poised",
    "12" = "Bivalent/poised",
    "13" = "Polycomb-repressed",
    "14" = "Polycomb-repressed",
    "15" = "Quiescent"
  )[as.character(state_number)])
}

read_roadmap_chromhmm_annotations <- function(path,
                                               epigenome_id,
                                               epigenome_label) {
  require_namespace("data.table")
  segmentation <- data.table::fread(
    path,
    sep = "\t",
    header = FALSE,
    col.names = c("chromosome", "start0", "end", "state"),
    showProgress = FALSE,
    data.table = FALSE
  )
  category <- collapse_roadmap_state(segmentation$state)
  if (anyNA(category)) {
    stop("Roadmap file contains an unexpected ChromHMM state.")
  }
  data.frame(
    chromosome = normalize_autosome(segmentation$chromosome),
    start = as.integer(segmentation$start0) + 1L,
    end = as.integer(segmentation$end),
    annotation = paste0(
      "Roadmap ", epigenome_id, " ", epigenome_label, ": ", category
    ),
    stringsAsFactors = FALSE
  )
}

compute_local_variant_density <- function(variants, window_bp = 1000000L) {
  required_columns <- c("variant_id", "chromosome", "position")
  if (!all(required_columns %in% names(variants)) ||
      anyDuplicated(variants$variant_id) || window_bp < 1L) {
    stop("Invalid variant coordinates or density window.")
  }
  output <- integer(nrow(variants))
  chromosome_groups <- split(seq_len(nrow(variants)), variants$chromosome)
  for (indices in chromosome_groups) {
    ordering <- indices[order(variants$position[indices], method = "radix")]
    positions <- variants$position[ordering]
    left <- 1L
    right <- 0L
    counts <- integer(length(ordering))
    for (index in seq_along(ordering)) {
      while (left < index && positions[index] - positions[left] > window_bp) {
        left <- left + 1L
      }
      if (right < index) {
        right <- index
      }
      while (right < length(ordering) &&
             positions[right + 1L] - positions[index] <= window_bp) {
        right <- right + 1L
      }
      counts[index] <- right - left + 1L
    }
    output[ordering] <- counts
  }
  output
}

quantile_bin <- function(values, n_bins) {
  values <- as.numeric(values)
  output <- rep(NA_integer_, length(values))
  observed <- which(is.finite(values))
  if (!length(observed)) {
    return(output)
  }
  unique_values <- unique(values[observed])
  if (length(unique_values) == 1L) {
    output[observed] <- 1L
    return(output)
  }
  breaks <- unique(as.numeric(stats::quantile(
    values[observed],
    probs = seq(0, 1, length.out = n_bins + 1L),
    na.rm = TRUE,
    names = FALSE,
    type = 8
  )))
  if (length(breaks) < 2L) {
    output[observed] <- 1L
    return(output)
  }
  output[observed] <- as.integer(cut(
    values[observed],
    breaks = breaks,
    include.lowest = TRUE,
    labels = FALSE
  ))
  output
}

derive_matching_strata <- function(variant_table, n_bins = 10L) {
  required_columns <- c(
    "variant_id", "chromosome", "minor_allele_frequency",
    "minimum_target_tss_distance", "local_tested_variant_count_1mb",
    "n_tested_genes"
  )
  if (!all(required_columns %in% names(variant_table)) ||
      anyDuplicated(variant_table$variant_id) || n_bins < 2L) {
    stop("Variant matching covariates have an unexpected structure.")
  }
  variant_table$maf_bin <- quantile_bin(
    variant_table$minor_allele_frequency,
    n_bins
  )
  variant_table$distance_bin <- quantile_bin(
    log10(variant_table$minimum_target_tss_distance + 1),
    n_bins
  )
  variant_table$density_bin <- quantile_bin(
    log1p(variant_table$local_tested_variant_count_1mb),
    n_bins
  )
  variant_table$tested_gene_bin <- cut(
    variant_table$n_tested_genes,
    breaks = c(-Inf, 1, 2, 4, Inf),
    labels = c("1", "2", "3-4", "5+"),
    right = TRUE
  )
  strata_from_bins(variant_table)
}

# Split out from derive_matching_strata so that a reproduction run can rebuild
# the exact strata from exported bin columns instead of re-deriving quantiles
# from floating-point covariates. Writing MAF to text loses about 5e-16, which
# is enough to move roughly 17,000 variants across a decile boundary; keying off
# the bins removes that dependence entirely.
strata_from_bins <- function(variant_table) {
  required_columns <- c(
    "chromosome", "maf_bin", "distance_bin", "density_bin", "tested_gene_bin"
  )
  if (!all(required_columns %in% names(variant_table))) {
    stop("Matching bins are missing.")
  }
  if (anyNA(variant_table[, required_columns])) {
    stop("Matching strata contain missing values.")
  }
  variant_table$full_stratum <- interaction(
    variant_table$chromosome,
    variant_table$maf_bin,
    variant_table$distance_bin,
    variant_table$density_bin,
    variant_table$tested_gene_bin,
    drop = TRUE,
    lex.order = TRUE
  )
  variant_table$drop_density_stratum <- interaction(
    variant_table$chromosome,
    variant_table$maf_bin,
    variant_table$distance_bin,
    variant_table$tested_gene_bin,
    drop = TRUE,
    lex.order = TRUE
  )
  variant_table$drop_distance_stratum <- interaction(
    variant_table$chromosome,
    variant_table$maf_bin,
    variant_table$tested_gene_bin,
    drop = TRUE,
    lex.order = TRUE
  )
  variant_table
}

prepare_matching_index <- function(variant_table, selected_ids) {
  selected_ids <- unique(as.character(selected_ids))
  if (!length(selected_ids) || any(!selected_ids %in% variant_table$variant_id)) {
    stop("Selected IDs are invalid.")
  }
  selected_rows <- match(selected_ids, variant_table$variant_id)
  control_universe <- which(!variant_table$variant_id %in% selected_ids)
  if (!length(control_universe)) {
    stop("No eligible matched-control variants remain.")
  }

  numeric_covariates <- c(
    "minor_allele_frequency", "minimum_target_tss_distance",
    "local_tested_variant_count_1mb", "n_tested_genes"
  )
  transformed <- cbind(
    variant_table$minor_allele_frequency,
    log10(variant_table$minimum_target_tss_distance + 1),
    log1p(variant_table$local_tested_variant_count_1mb),
    log1p(variant_table$n_tested_genes)
  )
  colnames(transformed) <- numeric_covariates
  scales <- apply(transformed[control_universe, , drop = FALSE], 2L, stats::sd)
  scales[!is.finite(scales) | scales == 0] <- 1

  split_pool <- function(stratum) {
    split(control_universe, as.character(stratum[control_universe]))
  }
  full_pools <- split_pool(variant_table$full_stratum)
  drop_density_pools <- split_pool(variant_table$drop_density_stratum)
  drop_distance_pools <- split_pool(variant_table$drop_distance_stratum)
  chromosome_maf_pools <- split(
    control_universe,
    paste(
      variant_table$chromosome[control_universe],
      variant_table$maf_bin[control_universe],
      sep = "::"
    )
  )

  pools <- vector("list", length(selected_rows))
  relaxation_levels <- integer(length(selected_rows))
  for (selected_index in seq_along(selected_rows)) {
    row_index <- selected_rows[selected_index]
    candidate_pools <- list(
      full_pools[[as.character(variant_table$full_stratum[row_index])]],
      drop_density_pools[[as.character(
        variant_table$drop_density_stratum[row_index]
      )]],
      drop_distance_pools[[as.character(
        variant_table$drop_distance_stratum[row_index]
      )]]
    )
    relaxation_level <- which(lengths(candidate_pools) > 0L)[1L] - 1L
    if (is.na(relaxation_level)) {
      chromosome_maf_key <- paste(
        variant_table$chromosome[row_index],
        variant_table$maf_bin[row_index],
        sep = "::"
      )
      candidate_pool <- chromosome_maf_pools[[chromosome_maf_key]]
      if (!length(candidate_pool)) {
        stop("No same-chromosome, same-MAF-bin controls for ",
             variant_table$variant_id[row_index], ".")
      }
      differences <- sweep(
        transformed[candidate_pool, , drop = FALSE],
        2L,
        transformed[row_index, ],
        "-"
      )
      distances <- rowSums(sweep(differences, 2L, scales, "/")^2)
      candidate_pool <- candidate_pool[
        order(distances, variant_table$variant_id[candidate_pool],
              method = "radix")
      ]
      pool <- head(candidate_pool, 100L)
      relaxation_level <- 3L
    } else {
      pool <- candidate_pools[[relaxation_level + 1L]]
    }
    pools[[selected_index]] <- pool
    relaxation_levels[selected_index] <- relaxation_level
  }
  list(
    selected_rows = selected_rows,
    pools = pools,
    relaxation_levels = relaxation_levels,
    variant_table = variant_table
  )
}

sample_from_matching_index <- function(matching_index,
                                       controls_per_variant,
                                       seed) {
  controls_per_variant <- as.integer(controls_per_variant)
  if (controls_per_variant < 1L) {
    stop("controls_per_variant must be positive.")
  }
  variant_table <- matching_index$variant_table
  set.seed(seed)
  sampled_rows <- matrix(
    NA_integer_,
    nrow = length(matching_index$selected_rows),
    ncol = controls_per_variant
  )
  for (selected_index in seq_along(matching_index$selected_rows)) {
    pool <- matching_index$pools[[selected_index]]
    sampled_rows[selected_index, ] <- pool[sample.int(
      length(pool),
      size = controls_per_variant,
      replace = length(pool) < controls_per_variant
    )]
  }
  selected_rows <- rep(
    matching_index$selected_rows,
    each = controls_per_variant
  )
  sampled_rows <- as.vector(t(sampled_rows))
  data.frame(
    selected_id = variant_table$variant_id[selected_rows],
    selected_chromosome = variant_table$chromosome[selected_rows],
    control_id = variant_table$variant_id[sampled_rows],
    control_chromosome = variant_table$chromosome[sampled_rows],
    relaxation_level = rep(
      as.integer(matching_index$relaxation_levels),
      each = controls_per_variant
    ),
    control_rank = rep(
      seq_len(controls_per_variant),
      times = length(matching_index$selected_rows)
    ),
    stringsAsFactors = FALSE
  )
}

sample_matched_controls <- function(variant_table,
                                    selected_ids,
                                    controls_per_variant = 5L,
                                    seed = 20260807L) {
  matching_index <- prepare_matching_index(variant_table, selected_ids)
  sample_from_matching_index(
    matching_index,
    controls_per_variant = controls_per_variant,
    seed = seed
  )
}

run_repeated_matching <- function(variant_table,
                                  selected_sets,
                                  seeds,
                                  controls_per_variant = 5L) {
  if (is.null(names(selected_sets)) || any(!nzchar(names(selected_sets))) ||
      anyDuplicated(names(selected_sets))) {
    stop("selected_sets must be a uniquely named list.")
  }
  rows <- vector("list", length(selected_sets) * length(seeds))
  row_index <- 0L
  for (set_name in names(selected_sets)) {
    matching_index <- prepare_matching_index(
      variant_table,
      selected_sets[[set_name]]
    )
    for (seed in seeds) {
      row_index <- row_index + 1L
      matched <- sample_from_matching_index(
        matching_index,
        controls_per_variant = controls_per_variant,
        seed = seed
      )
      matched$discovery_set <- set_name
      matched$match_seed <- as.integer(seed)
      rows[[row_index]] <- matched
    }
  }
  output <- do.call(rbind, rows)
  output[, c(
    "discovery_set", "match_seed", "selected_id", "selected_chromosome",
    "control_id", "control_chromosome", "relaxation_level", "control_rank"
  )]
}

summarize_covariate_balance <- function(variant_table, match_table) {
  numeric_covariates <- c(
    "minor_allele_frequency", "minimum_target_tss_distance",
    "local_tested_variant_count_1mb", "n_tested_genes"
  )
  groups <- split(
    seq_len(nrow(match_table)),
    interaction(match_table$discovery_set, match_table$match_seed, drop = TRUE)
  )
  rows <- lapply(groups, function(indices) {
    group <- match_table[indices, , drop = FALSE]
    selected <- variant_table[
      match(unique(group$selected_id), variant_table$variant_id),
      , drop = FALSE
    ]
    controls <- variant_table[
      match(group$control_id, variant_table$variant_id),
      , drop = FALSE
    ]
    do.call(rbind, lapply(numeric_covariates, function(covariate) {
      selected_values <- as.numeric(selected[[covariate]])
      control_values <- as.numeric(controls[[covariate]])
      pooled_sd <- sqrt((stats::var(selected_values) +
                           stats::var(control_values)) / 2)
      standardized_difference <- if (is.finite(pooled_sd) && pooled_sd > 0) {
        (mean(selected_values) - mean(control_values)) / pooled_sd
      } else {
        0
      }
      data.frame(
        discovery_set = group$discovery_set[1L],
        match_seed = group$match_seed[1L],
        covariate = covariate,
        selected_mean = mean(selected_values),
        control_mean = mean(control_values),
        standardized_mean_difference = standardized_difference,
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, rows)
}

run_streaming_matched_enrichment <- function(
    variant_table,
    annotation_matrix,
    selected_sets,
    seeds,
    controls_per_variant = 5L,
    chromosomes = as.character(1:22),
    minimum_overlap_count = 10L,
    verbose = FALSE) {
  annotation_columns <- setdiff(
    names(annotation_matrix),
    c("variant_id", "chromosome")
  )
  if (!length(annotation_columns) ||
      !all(vapply(annotation_matrix[annotation_columns], is.logical,
                  logical(1))) ||
      !identical(annotation_matrix$variant_id, variant_table$variant_id) ||
      is.null(names(selected_sets)) || any(!nzchar(names(selected_sets))) ||
      !length(seeds)) {
    stop("Invalid inputs for streaming matched enrichment.")
  }

  point_rows <- list()
  jackknife_rows <- list()
  balance_rows <- list()
  seed_rows <- list()
  relaxation_rows <- list()
  output_index <- 0L

  for (set_name in names(selected_sets)) {
    if (isTRUE(verbose)) {
      message(
        "  Matching set ",
        match(set_name, names(selected_sets)),
        "/",
        length(selected_sets),
        ": ",
        set_name,
        " (n = ",
        length(unique(selected_sets[[set_name]])),
        ")"
      )
    }
    selected_ids <- unique(as.character(selected_sets[[set_name]]))
    matching_index <- prepare_matching_index(variant_table, selected_ids)
    selected_rows <- match(selected_ids, annotation_matrix$variant_id)
    selected_chromosome <- variant_table$chromosome[
      match(selected_ids, variant_table$variant_id)
    ]
    selected_values <- as.matrix(
      annotation_matrix[selected_rows, annotation_columns, drop = FALSE]
    )
    storage.mode(selected_values) <- "logical"
    selected_overlap <- colSums(selected_values)
    selected_total <- length(selected_ids)

    selected_overlap_by_chromosome <- matrix(
      0,
      nrow = length(chromosomes),
      ncol = length(annotation_columns),
      dimnames = list(chromosomes, annotation_columns)
    )
    selected_total_by_chromosome <- setNames(
      integer(length(chromosomes)),
      chromosomes
    )
    for (chromosome in chromosomes) {
      keep <- selected_chromosome == chromosome
      selected_total_by_chromosome[chromosome] <- sum(keep)
      if (any(keep)) {
        selected_overlap_by_chromosome[chromosome, ] <- colSums(
          selected_values[keep, , drop = FALSE]
        )
      }
    }

    control_overlap_by_chromosome <- matrix(
      0,
      nrow = length(chromosomes),
      ncol = length(annotation_columns),
      dimnames = list(chromosomes, annotation_columns)
    )
    control_total_by_chromosome <- setNames(
      integer(length(chromosomes)),
      chromosomes
    )

    relaxation_table <- as.data.frame(table(
      relaxation_level = matching_index$relaxation_levels
    ), stringsAsFactors = FALSE)
    relaxation_table$discovery_set <- set_name
    relaxation_table$proportion <- relaxation_table$Freq / selected_total
    relaxation_rows[[length(relaxation_rows) + 1L]] <- relaxation_table[, c(
      "discovery_set", "relaxation_level", "Freq", "proportion"
    )]

    for (seed in seeds) {
      matched <- sample_from_matching_index(
        matching_index,
        controls_per_variant = controls_per_variant,
        seed = seed
      )
      matched$discovery_set <- set_name
      matched$match_seed <- as.integer(seed)
      control_rows <- match(matched$control_id, annotation_matrix$variant_id)
      control_values <- as.matrix(
        annotation_matrix[control_rows, annotation_columns, drop = FALSE]
      )
      storage.mode(control_values) <- "logical"

      for (chromosome in chromosomes) {
        keep <- matched$selected_chromosome == chromosome
        control_total_by_chromosome[chromosome] <-
          control_total_by_chromosome[chromosome] + sum(keep)
        if (any(keep)) {
          control_overlap_by_chromosome[chromosome, ] <-
            control_overlap_by_chromosome[chromosome, ] +
            colSums(control_values[keep, , drop = FALSE])
        }
      }

      balance_rows[[length(balance_rows) + 1L]] <-
        summarize_covariate_balance(variant_table, matched)
      seed_control_overlap <- colSums(control_values)
      seed_control_total <- nrow(control_values)
      seed_rows[[length(seed_rows) + 1L]] <- data.frame(
        discovery_set = set_name,
        match_seed = as.integer(seed),
        annotation = annotation_columns,
        log2_enrichment = log2(
          safe_overlap_proportion(selected_overlap, selected_total) /
            safe_overlap_proportion(seed_control_overlap, seed_control_total)
        ),
        stringsAsFactors = FALSE
      )
    }

    control_overlap <- colSums(control_overlap_by_chromosome)
    control_total <- sum(control_total_by_chromosome)
    valid_overlap <- selected_overlap >= minimum_overlap_count &
      control_overlap >= minimum_overlap_count
    enrichment <- rep(NA_real_, length(annotation_columns))
    enrichment[valid_overlap] <-
      safe_overlap_proportion(
        selected_overlap[valid_overlap],
        selected_total
      ) / safe_overlap_proportion(
        control_overlap[valid_overlap],
        control_total
      )
    point_log2 <- log2(enrichment)

    output_index <- output_index + 1L
    point_rows[[output_index]] <- data.frame(
      discovery_set = set_name,
      annotation = annotation_columns,
      selected_overlap = as.integer(selected_overlap),
      selected_total = selected_total,
      selected_proportion = selected_overlap / selected_total,
      control_overlap = as.integer(control_overlap),
      control_total = control_total,
      control_proportion = control_overlap / control_total,
      enrichment = enrichment,
      log2_enrichment = point_log2,
      stringsAsFactors = FALSE
    )

    leave_one_out <- matrix(
      NA_real_,
      nrow = length(chromosomes),
      ncol = length(annotation_columns),
      dimnames = list(chromosomes, annotation_columns)
    )
    for (chromosome in chromosomes) {
      selected_overlap_loo <- selected_overlap -
        selected_overlap_by_chromosome[chromosome, ]
      selected_total_loo <- selected_total -
        selected_total_by_chromosome[chromosome]
      control_overlap_loo <- control_overlap -
        control_overlap_by_chromosome[chromosome, ]
      control_total_loo <- control_total -
        control_total_by_chromosome[chromosome]
      leave_one_out[chromosome, ] <- log2(
        safe_overlap_proportion(
          selected_overlap_loo,
          selected_total_loo
        ) / safe_overlap_proportion(
          control_overlap_loo,
          control_total_loo
        )
      )
    }
    jackknife_se <- apply(leave_one_out, 2L, function(estimates) {
      estimates <- estimates[is.finite(estimates)]
      if (length(estimates) < 2L) {
        return(NA_real_)
      }
      sqrt(
        (length(estimates) - 1) / length(estimates) *
          sum((estimates - mean(estimates))^2)
      )
    })
    jackknife_rows[[output_index]] <- data.frame(
      discovery_set = set_name,
      annotation = annotation_columns,
      log2_enrichment = point_log2,
      jackknife_se = jackknife_se,
      ci_lower = point_log2 - 1.96 * jackknife_se,
      ci_upper = point_log2 + 1.96 * jackknife_se,
      n_jackknife_blocks = colSums(is.finite(leave_one_out)),
      stringsAsFactors = FALSE
    )
  }

  point_results <- do.call(rbind, point_rows)
  jackknife_results <- do.call(rbind, jackknife_rows)
  results <- summarize_enrichment_results(point_results, jackknife_results)
  list(
    results = results,
    point_results = point_results,
    jackknife_results = jackknife_results,
    matching_balance = do.call(rbind, balance_rows),
    seed_results = do.call(rbind, seed_rows),
    relaxation_audit = do.call(rbind, relaxation_rows)
  )
}

safe_pairwise_r2 <- function(first, second) {
  keep <- is.finite(first) & is.finite(second)
  if (sum(keep) < 3L || stats::sd(first[keep]) == 0 ||
      stats::sd(second[keep]) == 0) {
    return(NA_real_)
  }
  stats::cor(first[keep], second[keep])^2
}

classify_one_ld_direction <- function(source_variants,
                                      target_variants,
                                      dosage_matrix,
                                      window_bp,
                                      r2_threshold) {
  target_groups <- split(seq_len(nrow(target_variants)),
                         target_variants$chromosome)
  target_groups <- lapply(target_groups, function(indices) {
    indices[order(target_variants$position[indices], method = "radix")]
  })
  rows <- lapply(seq_len(nrow(source_variants)), function(index) {
    source <- source_variants[index, , drop = FALSE]
    chromosome_targets <- target_groups[[source$chromosome]]
    if (is.null(chromosome_targets)) {
      candidates <- integer()
    } else {
      target_positions <- target_variants$position[chromosome_targets]
      left <- findInterval(source$position - window_bp, target_positions) + 1L
      right <- findInterval(source$position + window_bp, target_positions)
      candidates <- if (left <= right) {
        chromosome_targets[seq.int(left, right)]
      } else {
        integer()
      }
    }
    exact_match <- source$variant_id %in% target_variants$variant_id
    if (!length(candidates)) {
      return(data.frame(
        variant_id = source$variant_id,
        exact_shared = exact_match,
        best_partner_id = NA_character_,
        best_r2 = NA_real_,
        ld_shared = exact_match,
        stringsAsFactors = FALSE
      ))
    }
    source_dosage <- dosage_matrix[source$variant_id, ]
    r2_values <- vapply(candidates, function(candidate) {
      safe_pairwise_r2(
        source_dosage,
        dosage_matrix[target_variants$variant_id[candidate], ]
      )
    }, numeric(1))
    if (all(is.na(r2_values))) {
      best_index <- NA_integer_
      best_r2 <- NA_real_
      best_partner <- NA_character_
    } else {
      best_index <- which.max(replace(r2_values, is.na(r2_values), -Inf))
      best_r2 <- r2_values[best_index]
      best_partner <- target_variants$variant_id[candidates[best_index]]
    }
    data.frame(
      variant_id = source$variant_id,
      exact_shared = exact_match,
      best_partner_id = best_partner,
      best_r2 = best_r2,
      ld_shared = exact_match || (!is.na(best_r2) && best_r2 >= r2_threshold),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

classify_cross_method_ld <- function(fash_variants,
                                     strober_variants,
                                     dosage_matrix,
                                     window_bp = 1000000L,
                                     r2_threshold = 0.8) {
  required_columns <- c("variant_id", "chromosome", "position")
  if (!all(required_columns %in% names(fash_variants)) ||
      !all(required_columns %in% names(strober_variants)) ||
      is.null(rownames(dosage_matrix)) ||
      !all(c(fash_variants$variant_id, strober_variants$variant_id) %in%
           rownames(dosage_matrix)) || window_bp < 1L ||
      r2_threshold <= 0 || r2_threshold > 1) {
    stop("Invalid inputs for cross-method LD classification.")
  }
  list(
    fash = classify_one_ld_direction(
      fash_variants, strober_variants, dosage_matrix, window_bp, r2_threshold
    ),
    strober = classify_one_ld_direction(
      strober_variants, fash_variants, dosage_matrix, window_bp, r2_threshold
    )
  )
}

estimate_annotation_enrichment <- function(annotation_matrix,
                                           selected_ids,
                                           control_ids,
                                           annotation,
                                           minimum_overlap_count) {
  selected_rows <- match(unique(selected_ids), annotation_matrix$variant_id)
  control_rows <- match(control_ids, annotation_matrix$variant_id)
  if (anyNA(selected_rows) || anyNA(control_rows)) {
    stop("Matched variant IDs are absent from the annotation matrix.")
  }
  selected_values <- as.logical(annotation_matrix[[annotation]][selected_rows])
  control_values <- as.logical(annotation_matrix[[annotation]][control_rows])
  selected_overlap <- sum(selected_values)
  control_overlap <- sum(control_values)
  selected_total <- length(selected_values)
  control_total <- length(control_values)
  if (selected_overlap < minimum_overlap_count ||
      control_overlap < minimum_overlap_count) {
    log2_enrichment <- NA_real_
    enrichment <- NA_real_
  } else {
    selected_adjusted <- safe_overlap_proportion(
      selected_overlap,
      selected_total
    )
    control_adjusted <- safe_overlap_proportion(
      control_overlap,
      control_total
    )
    enrichment <- selected_adjusted / control_adjusted
    log2_enrichment <- log2(enrichment)
  }
  data.frame(
    annotation = annotation,
    selected_overlap = selected_overlap,
    selected_total = selected_total,
    selected_proportion = selected_overlap / selected_total,
    control_overlap = control_overlap,
    control_total = control_total,
    control_proportion = control_overlap / control_total,
    enrichment = enrichment,
    log2_enrichment = log2_enrichment,
    stringsAsFactors = FALSE
  )
}

compute_matched_enrichment <- function(annotation_matrix,
                                       match_table,
                                       minimum_overlap_count = 10L) {
  annotation_columns <- setdiff(
    names(annotation_matrix),
    c("variant_id", "chromosome")
  )
  if (!length(annotation_columns) ||
      !all(vapply(annotation_matrix[annotation_columns], is.logical,
                  logical(1)))) {
    stop("The annotation matrix must contain logical annotation columns.")
  }
  if (!"discovery_set" %in% names(match_table)) {
    match_table$discovery_set <- "discovery"
  }
  groups <- split(seq_len(nrow(match_table)), match_table$discovery_set)
  rows <- lapply(groups, function(indices) {
    group <- match_table[indices, , drop = FALSE]
    estimates <- do.call(rbind, lapply(annotation_columns, function(annotation) {
      estimate_annotation_enrichment(
        annotation_matrix,
        selected_ids = group$selected_id,
        control_ids = group$control_id,
        annotation = annotation,
        minimum_overlap_count = minimum_overlap_count
      )
    }))
    estimates$discovery_set <- group$discovery_set[1L]
    estimates
  })
  output <- do.call(rbind, rows)
  output[, c("discovery_set", setdiff(names(output), "discovery_set"))]
}

chromosome_jackknife_enrichment <- function(annotation_matrix,
                                            match_table,
                                            chromosomes = as.character(1:22),
                                            minimum_overlap_count = 10L) {
  point <- compute_matched_enrichment(
    annotation_matrix,
    match_table,
    minimum_overlap_count = minimum_overlap_count
  )
  key <- interaction(point$discovery_set, point$annotation, drop = TRUE)
  leave_one_out <- vector("list", length(chromosomes))
  for (index in seq_along(chromosomes)) {
    chromosome <- chromosomes[index]
    retained <- match_table$selected_chromosome != chromosome
    estimate <- compute_matched_enrichment(
      annotation_matrix,
      match_table[retained, , drop = FALSE],
      minimum_overlap_count = 0L
    )
    estimate$omitted_chromosome <- chromosome
    estimate$key <- interaction(
      estimate$discovery_set,
      estimate$annotation,
      drop = TRUE
    )
    leave_one_out[[index]] <- estimate
  }
  leave_one_out <- do.call(rbind, leave_one_out)
  rows <- lapply(seq_len(nrow(point)), function(index) {
    estimates <- leave_one_out$log2_enrichment[leave_one_out$key == key[index]]
    estimates <- estimates[is.finite(estimates)]
    if (length(estimates) < 2L || !is.finite(point$log2_enrichment[index])) {
      jackknife_se <- NA_real_
    } else {
      jackknife_se <- sqrt(
        (length(estimates) - 1) / length(estimates) *
          sum((estimates - mean(estimates))^2)
      )
    }
    data.frame(
      discovery_set = point$discovery_set[index],
      annotation = point$annotation[index],
      log2_enrichment = point$log2_enrichment[index],
      jackknife_se = jackknife_se,
      ci_lower = point$log2_enrichment[index] - 1.96 * jackknife_se,
      ci_upper = point$log2_enrichment[index] + 1.96 * jackknife_se,
      n_jackknife_blocks = length(estimates),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_enrichment_results <- function(point_results, jackknife_results) {
  key_columns <- c("discovery_set", "annotation")
  output <- merge(point_results, jackknife_results, by = key_columns,
                  suffixes = c("", ".jackknife"), sort = FALSE)
  output$z_score <- output$log2_enrichment / output$jackknife_se
  output$p_value <- 2 * stats::pnorm(-abs(output$z_score))
  output$q_value <- stats::p.adjust(output$p_value, method = "BH")
  output
}
