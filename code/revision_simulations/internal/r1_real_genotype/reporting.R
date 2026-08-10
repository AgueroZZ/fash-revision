# Load, validate, and format cached results for the internal real-genotype R1 report.

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required real-genotype R1 cache file is missing: ", path)
  }
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

format_decimal <- function(x, digits = 3) {
  formatC(x, format = "f", digits = digits)
}

format_mc_interval <- function(mean, lower, upper, digits = 3) {
  paste0(
    format_decimal(mean, digits),
    " [",
    format_decimal(lower, digits),
    ", ",
    format_decimal(upper, digits),
    "]"
  )
}

format_mc_table <- function(summary_table, method_order) {
  summary_table$method <- factor(summary_table$method, levels = method_order)
  summary_table <- summary_table[
    order(summary_table$method),
    ,
    drop = FALSE
  ]
  data.frame(
    Method = as.character(summary_table$method),
    `Mean discoveries` = format_decimal(summary_table$mean_discoveries, 1),
    `Power (95% MC CI)` = format_mc_interval(
      summary_table$mean_power,
      summary_table$power_ci_lower,
      summary_table$power_ci_upper
    ),
    `Empirical FDR: E[FDP] (95% MC CI)` = format_mc_interval(
      summary_table$mean_fdr,
      summary_table$fdr_ci_lower,
      summary_table$fdr_ci_upper
    ),
    check.names = FALSE
  )
}

format_pi0_table <- function(summary_table) {
  summary_table$method <- factor(
    summary_table$method,
    levels = c("FASH-IWP1", "FASH-linear")
  )
  summary_table$fit <- factor(
    summary_table$fit,
    levels = c("Raw", "BF-corrected")
  )
  summary_table <- summary_table[
    order(summary_table$method, summary_table$fit),
    ,
    drop = FALSE
  ]
  data.frame(
    Method = as.character(summary_table$method),
    Fit = as.character(summary_table$fit),
    `Mean estimated pi0 (95% MC CI)` = format_mc_interval(
      summary_table$mean_estimated_pi0,
      summary_table$pi0_ci_lower,
      summary_table$pi0_ci_upper
    ),
    check.names = FALSE
  )
}

get_mc_metric <- function(summary_table, method, metric) {
  row <- summary_table[summary_table$method == method, , drop = FALSE]
  if (nrow(row) != 1L || !metric %in% names(row)) {
    stop("Could not extract a unique Monte Carlo metric.")
  }
  row[[metric]]
}

render_scrollable_table <- function(
  table,
  align = NULL,
  minimum_width = "760px",
  digits = NULL
) {
  arguments <- list(
    x = table,
    format = "html",
    align = align,
    escape = TRUE,
    table.attr = paste0(
      'class="table table-striped table-hover" style="min-width:',
      minimum_width,
      ';"'
    )
  )
  if (!is.null(digits)) {
    arguments$digits <- digits
  }
  rendered <- do.call(knitr::kable, arguments)
  cat(
    '<div class="r1-real-table-scroll">\n',
    as.character(rendered),
    '\n</div>\n',
    sep = ""
  )
  invisible(rendered)
}

expected_seeds <- c(12345L, 22345L, 32345L, 42345L, 52345L)
expected_sample_ids <- c(
  "18489", "18499", "18505", "18508", "18511", "18517", "18520",
  "18858", "18870", "18907", "18912", "19093", "19108", "19159",
  "19190", "19209", "18855", "19127", "19193"
)
expected_class_probs <- c(
  dynamic_bspline = 0.20,
  constant = 0.40,
  zero = 0.40
)
expected_class_counts <- c(
  dynamic_bspline = 200L,
  constant = 400L,
  zero = 400L
)
expected_fash_methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "FASH-linear-Raw",
  "FASH-linear-BF"
)
expected_direct_methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "Direct-linear-LRT-eFDR-true-pi0",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
expected_all_methods <- union(expected_fash_methods, expected_direct_methods)
expected_vcf_md5 <- "4f9eb383ce3512d867b42ab806d451a8"
expected_pair_summary_md5 <- "23e7f6a0093309059424207872dfa1e0"
expected_repeated_loci <- c("ENSG00000100027", "ENSG00000115170")

mc_output_dir <- file.path(
  "output",
  "revision_simulations",
  "internal",
  "r1_real_genotype_locus_blocks_pilot5"
)
summary_dir <- file.path(mc_output_dir, "summary")
configuration_path <- file.path(mc_output_dir, "configuration.rds")
genotype_cache_path <- file.path(
  mc_output_dir,
  "genotype_samples",
  "real_genotype_samples.rds"
)
full_fit_path <- file.path(mc_output_dir, "full_fits", "seed_12345.rds")

required_rds <- c(configuration_path, genotype_cache_path, full_fit_path)
if (any(!file.exists(required_rds))) {
  stop(
    "The internal real-genotype R1 cache is incomplete: ",
    paste(required_rds[!file.exists(required_rds)], collapse = ", ")
  )
}

configuration <- readRDS(configuration_path)
genotype_cache <- readRDS(genotype_cache_path)
out <- readRDS(full_fit_path)

fash_alpha <- read_required_csv(file.path(
  summary_dir,
  "iwp_vs_linear_fash_mc_alpha_curve.csv"
))
fash_alpha_005 <- read_required_csv(file.path(
  summary_dir,
  "iwp_vs_linear_fash_mc_alpha005_summary.csv"
))
direct_alpha <- read_required_csv(file.path(
  summary_dir,
  "iwp_fash_vs_direct_true_pi0_mc_alpha_curve.csv"
))
direct_alpha_005 <- read_required_csv(file.path(
  summary_dir,
  "iwp_fash_vs_direct_true_pi0_mc_alpha005_summary.csv"
))
pi0_summary <- read_required_csv(file.path(summary_dir, "mc_pi0_summary.csv"))
variant_metadata <- read_required_csv(file.path(
  summary_dir,
  "real_genotype_variant_metadata.csv"
))
locus_ld <- read_required_csv(file.path(
  summary_dir,
  "real_genotype_locus_ld.csv"
))
seed_ld <- read_required_csv(file.path(
  summary_dir,
  "real_genotype_seed_ld.csv"
))
input_provenance <- read_required_csv(file.path(
  summary_dir,
  "real_genotype_input_provenance.csv"
))
extraction_summary <- read_required_csv(file.path(
  summary_dir,
  "real_genotype_extraction_summary.csv"
))

configuration_valid <-
  identical(configuration$format_version, 1L) &&
  identical(
    configuration$scenario,
    "r1_real_genotype_locus_blocks_random_bspline"
  ) &&
  identical(configuration$J, 1000L) &&
  identical(configuration$n_donors, 19L) &&
  identical(configuration$sample_ids, expected_sample_ids) &&
  isTRUE(all.equal(configuration$time_grid, 0:15)) &&
  identical(configuration$n_covariates, 5L) &&
  isTRUE(all.equal(configuration$expression_noise_sd, 1)) &&
  isTRUE(all.equal(configuration$dynamic_main_effect_sd, 1)) &&
  identical(configuration$num_basis, 20L) &&
  isTRUE(all.equal(configuration$class_probs, expected_class_probs)) &&
  identical(configuration$expected_class_counts, expected_class_counts) &&
  isTRUE(all.equal(configuration$dynamic_amplitude, 2)) &&
  isTRUE(all.equal(configuration$bspline_df, 6)) &&
  isTRUE(all.equal(configuration$bspline_coefficient_sd, 1)) &&
  isTRUE(all.equal(configuration$constant_sd, 1)) &&
  isTRUE(all.equal(configuration$covariate_effect_sd, 0.5)) &&
  isTRUE(all.equal(configuration$intercept_sd, 0)) &&
  isTRUE(configuration$apply_t_se_correction) &&
  identical(configuration$n_loci, 20L) &&
  identical(configuration$variants_per_locus, 50L) &&
  isTRUE(all.equal(configuration$maf_min, 0.10)) &&
  identical(configuration$candidate_genes_per_seed, 120L) &&
  identical(configuration$minimum_raw_variants, 80L) &&
  identical(configuration$gene_priority_seed_offset, 7101L) &&
  identical(configuration$block_start_seed_offset, 7201L) &&
  identical(configuration$across_locus_seed_offset, 7401L) &&
  identical(configuration$efdr_permutations, 100L) &&
  identical(configuration$permutation_seed_rule, "seed + 10000") &&
  isTRUE(all.equal(configuration$true_pi0, 0.8)) &&
  identical(configuration$full_fit_seed, 12345L) &&
  identical(configuration$seed_list, expected_seeds) &&
  identical(configuration$vcf_fingerprint$md5, expected_vcf_md5) &&
  identical(
    configuration$pair_summary_fingerprint$md5,
    expected_pair_summary_md5
  )
if (!configuration_valid) {
  stop("The real-genotype R1 configuration does not match the intended pilot.")
}

if (!is.list(genotype_cache) ||
    !isTRUE(all.equal(genotype_cache$configuration, configuration)) ||
    !identical(genotype_cache$sample_ids, expected_sample_ids) ||
    !identical(names(genotype_cache$samples), as.character(expected_seeds)) ||
    length(genotype_cache$samples) != length(expected_seeds)) {
  stop("The real-genotype sample cache does not match the configuration.")
}

expected_pair_count <- 20L * choose(50L, 2L)
for (seed in expected_seeds) {
  sample <- genotype_cache$samples[[as.character(seed)]]
  if (!is.list(sample) ||
      !identical(dim(sample$G), c(19L, 1000L)) ||
      !identical(rownames(sample$G), expected_sample_ids) ||
      nrow(sample$variant_info) != 1000L ||
      !identical(colnames(sample$G), sample$variant_info$unit_id) ||
      any(!is.finite(sample$G)) ||
      any(apply(sample$G, 2L, stats::sd) <= 0) ||
      any(sample$variant_info$observed_maf < 0.10) ||
      any(sample$variant_info$observed_maf > 0.50) ||
      any(sample$variant_info$sampling_seed != seed) ||
      length(sample$selected_genes) != 20L ||
      anyDuplicated(sample$variant_info$unit_id)) {
    stop("A sampled real-genotype matrix is incomplete for seed ", seed, ".")
  }

  block_counts <- table(factor(
    sample$variant_info$block_index,
    levels = seq_len(20L)
  ))
  block_metadata <- split(
    sample$variant_info,
    sample$variant_info$block_index
  )
  block_valid <- vapply(
    block_metadata,
    function(block) {
      length(unique(block$gene_id)) == 1L &&
        length(unique(block$chromosome)) == 1L &&
        !is.unsorted(block$position)
    },
    logical(1)
  )
  observed_maf <- pmin(
    colMeans(sample$G) / 2,
    1 - colMeans(sample$G) / 2
  )
  if (length(block_counts) != 20L ||
      any(block_counts != 50L) ||
      !all(block_valid) ||
      max(abs(observed_maf - sample$variant_info$observed_maf)) > 1e-8 ||
      nrow(sample$ld$locus) != 20L ||
      nrow(sample$ld$seed) != 1L ||
      length(sample$ld$within_abs_r) != expected_pair_count ||
      length(sample$ld$across_abs_r) != expected_pair_count ||
      any(!is.finite(sample$ld$within_abs_r)) ||
      any(!is.finite(sample$ld$across_abs_r))) {
    stop("Real-genotype block or LD validation failed for seed ", seed, ".")
  }
}

if (nrow(variant_metadata) != 5000L ||
    nrow(locus_ld) != 100L ||
    nrow(seed_ld) != 5L ||
    !identical(sort(unique(variant_metadata$sampling_seed)), expected_seeds) ||
    !identical(sort(unique(locus_ld$sampling_seed)), expected_seeds) ||
    !identical(sort(seed_ld$sampling_seed), expected_seeds) ||
    any(variant_metadata$observed_maf < 0.10) ||
    any(locus_ld$n_variants != 50L) ||
    any(seed_ld$n_loci != 20L) ||
    any(seed_ld$n_variants != 1000L)) {
  stop("The saved real-genotype metadata or LD summaries are incomplete.")
}

selected_locus_counts <- table(locus_ld$gene_id)
repeated_loci <- sort(names(selected_locus_counts[selected_locus_counts > 1L]))
if (length(selected_locus_counts) != 98L ||
    !identical(repeated_loci, expected_repeated_loci) ||
    any(selected_locus_counts[repeated_loci] != 2L) ||
    any(seed_ld$pooled_within_locus_median_abs_r <=
        seed_ld$matched_across_locus_median_abs_r)) {
  stop("The five-seed locus-selection or LD diagnostic is unexpected.")
}

if (nrow(input_provenance) != 2L ||
    !all(c("YRI genotype VCF", "Paper gene-variant summary") %in%
         input_provenance$artifact)) {
  stop("The real-genotype input provenance table is incomplete.")
}
vcf_row <- input_provenance[
  input_provenance$artifact == "YRI genotype VCF",
  ,
  drop = FALSE
]
pair_row <- input_provenance[
  input_provenance$artifact == "Paper gene-variant summary",
  ,
  drop = FALSE
]
if (nrow(vcf_row) != 1L ||
    nrow(pair_row) != 1L ||
    !identical(vcf_row$md5[[1L]], expected_vcf_md5) ||
    !identical(pair_row$md5[[1L]], expected_pair_summary_md5) ||
    !identical(vcf_row$file_name[[1L]], configuration$vcf_fingerprint$file_name) ||
    !identical(
      pair_row$file_name[[1L]],
      configuration$pair_summary_fingerprint$file_name
    ) ||
    !isTRUE(all.equal(
      vcf_row$size_bytes[[1L]],
      configuration$vcf_fingerprint$size_bytes
    )) ||
    !isTRUE(all.equal(
      pair_row$size_bytes[[1L]],
      configuration$pair_summary_fingerprint$size_bytes
    ))) {
  stop("The cached input fingerprints do not match the saved provenance.")
}

expected_extraction <- c(
  requested_variants = 94496L,
  matched_vcf_variants = 94496L,
  retained_common_polymorphic_variants = 94496L,
  candidate_genes = 576L
)
if (nrow(extraction_summary) != 1L ||
    !identical(
      as.integer(extraction_summary[1L, names(expected_extraction)]),
      unname(expected_extraction)
    ) ||
    !isTRUE(all.equal(
      extraction_summary,
      genotype_cache$extraction_summary,
      check.attributes = FALSE
    ))) {
  stop("The real-genotype extraction summary is incomplete or mismatched.")
}

validate_mc_summary <- function(table, methods, n_alpha) {
  nrow(table) == length(methods) * n_alpha &&
    identical(sort(unique(table$method)), sort(methods)) &&
    all(table$n_replications == 5L) &&
    all(is.finite(table$mean_power)) &&
    all(is.finite(table$mean_fdr))
}

expected_alpha_grid <- seq(0, 0.2, by = 0.005)
if (!validate_mc_summary(fash_alpha, expected_fash_methods, 41L) ||
    !validate_mc_summary(direct_alpha, expected_direct_methods, 41L) ||
    !isTRUE(all.equal(sort(unique(fash_alpha$alpha)), expected_alpha_grid)) ||
    !isTRUE(all.equal(sort(unique(direct_alpha$alpha)), expected_alpha_grid)) ||
    !validate_mc_summary(fash_alpha_005, expected_fash_methods, 1L) ||
    !validate_mc_summary(direct_alpha_005, expected_direct_methods, 1L) ||
    any(abs(fash_alpha_005$alpha - 0.05) > 1e-12) ||
    any(abs(direct_alpha_005$alpha - 0.05) > 1e-12) ||
    nrow(pi0_summary) != 4L ||
    !identical(sort(unique(pi0_summary$method)),
               sort(c("FASH-IWP1", "FASH-linear"))) ||
    !identical(sort(unique(pi0_summary$fit)),
               sort(c("Raw", "BF-corrected"))) ||
    any(pi0_summary$n_replications != 5L)) {
  stop("The real-genotype R1 Monte Carlo summaries are incomplete.")
}

observed_class_counts <- table(out$unit_info$effect_class)
observed_class_counts <- observed_class_counts[names(expected_class_counts)]
full_fit_valid <-
  identical(out$settings$scenario,
            "r1_real_genotype_locus_blocks_random_bspline") &&
  identical(out$settings$seed, 12345L) &&
  identical(out$settings$n_donors, 19L) &&
  identical(out$settings$n_variants, 1000L) &&
  identical(out$settings$n_covariates, 5L) &&
  isTRUE(all.equal(out$settings$class_probs, expected_class_probs)) &&
  isTRUE(all.equal(out$settings$expression_noise_sd, 1)) &&
  isTRUE(all.equal(out$settings$covariate_effect_sd, 0.5)) &&
  isTRUE(all.equal(out$settings$intercept_sd, 0)) &&
  isTRUE(all.equal(out$settings$dynamic_amplitude, 2)) &&
  isTRUE(all.equal(out$settings$bspline_df, 6)) &&
  isTRUE(all.equal(out$settings$dynamic_main_effect_sd, 1)) &&
  isTRUE(out$settings$apply_t_se_correction) &&
  identical(out$settings$genotype_source,
            "paper-derived YRI dosage locus blocks") &&
  identical(out$settings$n_loci, 20L) &&
  identical(out$settings$variants_per_locus, 50L) &&
  isTRUE(all.equal(out$settings$maf_min, 0.10)) &&
  identical(out$settings$sample_ids, expected_sample_ids) &&
  identical(out$settings$vcf_md5, expected_vcf_md5) &&
  identical(out$settings$pair_summary_md5, expected_pair_summary_md5) &&
  identical(dim(out$genotype), c(19L, 1000L)) &&
  identical(rownames(out$genotype), expected_sample_ids) &&
  identical(colnames(out$genotype), out$variant_info$unit_id) &&
  identical(unname(as.integer(observed_class_counts)),
            unname(expected_class_counts)) &&
  all(out$eqtl_summary$df == 12L) &&
  all(out$eqtl_summary$rank == 7L) &&
  all(expected_all_methods %in% unique(out$alpha_curve$method))
if (!full_fit_valid) {
  stop("The full seed-12345 fit does not match the real-genotype R1 pilot.")
}

first_real_sample <- genotype_cache$samples[["12345"]]

input_provenance_table <- data.frame(
  Artifact = input_provenance$artifact,
  File = input_provenance$file_name,
  `Size (MiB)` = format_decimal(input_provenance$size_bytes / 1024^2, 1),
  `Modified (UTC)` = input_provenance$modification_time_utc,
  MD5 = input_provenance$md5,
  check.names = FALSE
)

extraction_table <- data.frame(
  Metric = c(
    "Requested candidate variants",
    "Matched VCF variants",
    "Retained common polymorphic variants",
    "Candidate genes represented"
  ),
  Value = format(
    as.integer(extraction_summary[1L, names(expected_extraction)]),
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  ),
  check.names = FALSE
)

sampling_design_table <- data.frame(
  Setting = c(
    "Donors",
    "Time points",
    "Locus blocks per seed",
    "Position-contiguous variants per block",
    "Real genotype columns per seed",
    "Observed MAF threshold",
    "Monte Carlo seeds",
    "Unique loci across 100 seed-by-block selections"
  ),
  Value = c(
    "19",
    "16 (0 through 15)",
    "20",
    "50",
    "1,000",
    "at least 0.10",
    paste(expected_seeds, collapse = ", "),
    "98 (two loci were each selected in two independent seeds)"
  ),
  check.names = FALSE
)

ld_seed_table <- data.frame(
  Seed = seed_ld$sampling_seed,
  `Median locus span (kb)` = format_decimal(
    seed_ld$median_locus_span_bp / 1000,
    1
  ),
  `Within-locus median |r|` = format_decimal(
    seed_ld$pooled_within_locus_median_abs_r
  ),
  `Pairs with |r| >= 0.8` = paste0(
    format_decimal(
      100 * seed_ld$pooled_within_locus_fraction_abs_r_ge_0_8,
      1
    ),
    "%"
  ),
  `Mean adjacent |r|` = format_decimal(seed_ld$mean_locus_adjacent_abs_r),
  `Matched across-locus median |r|` = format_decimal(
    seed_ld$matched_across_locus_median_abs_r
  ),
  `Within-minus-across` = format_decimal(
    seed_ld$within_minus_across_median_abs_r
  ),
  check.names = FALSE
)

mean_within_ld <- mean(seed_ld$pooled_within_locus_median_abs_r)
mean_across_ld <- mean(seed_ld$matched_across_locus_median_abs_r)
mean_ld_gap <- mean(seed_ld$within_minus_across_median_abs_r)
mean_high_ld_fraction <- mean(
  seed_ld$pooled_within_locus_fraction_abs_r_ge_0_8
)
mean_median_locus_span_kb <- mean(seed_ld$median_locus_span_bp) / 1000
maf_range <- range(variant_metadata$observed_maf)

ld_summary_table <- data.frame(
  Diagnostic = c(
    "Mean pooled within-locus median |r|",
    "Mean matched across-locus median |r|",
    "Mean within-minus-across median |r|",
    "Mean fraction of within-locus pairs with |r| >= 0.8",
    "Mean seed-specific median locus span"
  ),
  Value = c(
    format_decimal(mean_within_ld),
    format_decimal(mean_across_ld),
    format_decimal(mean_ld_gap),
    paste0(format_decimal(100 * mean_high_ld_fraction, 1), "%"),
    paste0(format_decimal(mean_median_locus_span_kb, 1), " kb")
  ),
  check.names = FALSE
)

repeated_loci_table <- do.call(
  rbind,
  lapply(expected_repeated_loci, function(gene_id) {
    rows <- locus_ld[locus_ld$gene_id == gene_id, , drop = FALSE]
    data.frame(
      Locus = gene_id,
      Seeds = paste(rows$sampling_seed, collapse = ", "),
      `Selections across five seeds` = nrow(rows),
      check.names = FALSE
    )
  })
)
rownames(repeated_loci_table) <- NULL

fash_table <- format_mc_table(fash_alpha_005, expected_fash_methods)
direct_table <- format_mc_table(direct_alpha_005, expected_direct_methods)
pi0_table <- format_pi0_table(pi0_summary)
fash_iwp_raw_fdr_interval <- format_mc_interval(
  get_mc_metric(fash_alpha_005, "FASH-IWP1-Raw", "mean_fdr"),
  get_mc_metric(fash_alpha_005, "FASH-IWP1-Raw", "fdr_ci_lower"),
  get_mc_metric(fash_alpha_005, "FASH-IWP1-Raw", "fdr_ci_upper")
)
fash_iwp_bf_fdr_interval <- format_mc_interval(
  get_mc_metric(fash_alpha_005, "FASH-IWP1-BF", "mean_fdr"),
  get_mc_metric(fash_alpha_005, "FASH-IWP1-BF", "fdr_ci_lower"),
  get_mc_metric(fash_alpha_005, "FASH-IWP1-BF", "fdr_ci_upper")
)
direct_iwp_bf_power <- get_mc_metric(
  direct_alpha_005,
  "FASH-IWP1-BF",
  "mean_power"
)
direct_linear_power <- get_mc_metric(
  direct_alpha_005,
  "Direct-linear-LRT-eFDR-true-pi0",
  "mean_power"
)
direct_quadratic_power <- get_mc_metric(
  direct_alpha_005,
  "Direct-quadratic-LRT-eFDR-true-pi0",
  "mean_power"
)
direct_iwp_bf_fdr <- get_mc_metric(
  direct_alpha_005,
  "FASH-IWP1-BF",
  "mean_fdr"
)
direct_linear_fdr <- get_mc_metric(
  direct_alpha_005,
  "Direct-linear-LRT-eFDR-true-pi0",
  "mean_fdr"
)
direct_quadratic_fdr <- get_mc_metric(
  direct_alpha_005,
  "Direct-quadratic-LRT-eFDR-true-pi0",
  "mean_fdr"
)

plot_real_genotype_ld_heatmap <- function(
  sample = first_real_sample,
  main = "Seed 12345 real-genotype LD"
) {
  correlation <- stats::cor(sample$G)
  block_size <- unique(table(sample$variant_info$block_index))
  if (length(block_size) != 1L || block_size != 50L) {
    stop("The LD heatmap requires equal 50-variant locus blocks.")
  }
  palette <- grDevices::colorRampPalette(
    c("#2166AC", "#F7F7F7", "#B2182B")
  )(201L)
  centers <- seq(25.5, ncol(sample$G), by = 50)
  boundaries <- seq(50.5, ncol(sample$G) - 0.5, by = 50)
  graphics::image(
    x = seq_len(ncol(sample$G)),
    y = seq_len(ncol(sample$G)),
    z = correlation,
    col = palette,
    zlim = c(-1, 1),
    useRaster = TRUE,
    axes = FALSE,
    xlab = "Locus block",
    ylab = "Locus block",
    main = main,
    sub = "Signed dosage correlation: blue = -1, white = 0, red = +1"
  )
  graphics::axis(1, at = centers, labels = seq_len(20L), cex.axis = 0.65)
  graphics::axis(2, at = centers, labels = seq_len(20L), cex.axis = 0.65)
  graphics::abline(v = boundaries, h = boundaries, col = "#FFFFFF", lwd = 0.5)
  graphics::box()
  invisible(correlation)
}

plot_real_genotype_ld_distribution <- function(
  samples = genotype_cache$samples
) {
  within <- unlist(
    lapply(samples, function(sample) sample$ld$within_abs_r),
    use.names = FALSE
  )
  across <- unlist(
    lapply(samples, function(sample) sample$ld$across_abs_r),
    use.names = FALSE
  )
  within_density <- stats::density(within, from = 0, to = 1, cut = 0)
  across_density <- stats::density(across, from = 0, to = 1, cut = 0)
  y_max <- max(within_density$y, across_density$y)
  graphics::plot(
    within_density,
    xlim = c(0, 1),
    ylim = c(0, y_max * 1.05),
    lwd = 2.5,
    col = "#B2182B",
    xlab = "Absolute dosage correlation |r|",
    ylab = "Density",
    main = "Within-locus versus matched across-locus LD",
    sub = "Five independent genotype samples; matched pair counts per seed"
  )
  graphics::lines(across_density, lwd = 2.5, col = "#2166AC")
  graphics::legend(
    "topright",
    legend = c("Within locus", "Matched across loci"),
    col = c("#B2182B", "#2166AC"),
    lwd = 2.5,
    bty = "n"
  )
  graphics::abline(
    v = c(mean_within_ld, mean_across_ld),
    col = c("#B2182B", "#2166AC"),
    lty = 3
  )
  invisible(list(within = within, across = across))
}
