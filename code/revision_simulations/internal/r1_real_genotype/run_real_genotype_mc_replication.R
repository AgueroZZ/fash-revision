#!/usr/bin/env Rscript

# Run an internal R1-style Monte Carlo experiment using real YRI genotype
# dosage blocks sampled from loci that contributed gene-variant pairs to the
# paper analysis.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists("coderepo-local/code/revision_simulations/shared/simulation_functions.R")) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[[1L]] == length(args)) {
    return(default)
  }
  args[hit[[1L]] + 1L]
}

as_flag <- function(x) {
  tolower(x) %in% c("1", "true", "t", "yes", "y")
}

parse_seed_list <- function(x) {
  pieces <- trimws(strsplit(x, ",", fixed = TRUE)[[1L]])
  if (length(pieces) == 0L || any(!nzchar(pieces))) {
    stop("--seed-list must contain one or more comma-separated integer seeds.")
  }
  seeds <- suppressWarnings(as.integer(pieces))
  if (anyNA(seeds) || anyDuplicated(seeds)) {
    stop("--seed-list must contain unique integer seeds.")
  }
  seeds
}

resolve_input_path <- function(path, workflowr_root, argument_name) {
  candidates <- unique(c(path, file.path(workflowr_root, path)))
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    stop(argument_name, " does not exist: ", path)
  }
  normalizePath(existing[[1L]], winslash = "/", mustWork = TRUE)
}

resolve_output_root <- function(path, workflowr_root) {
  if (!grepl("^/", path)) {
    path <- file.path(workflowr_root, path)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

artifact_fingerprint <- function(path) {
  information <- file.info(path)
  checksum <- unname(tools::md5sum(path))
  if (is.na(checksum)) {
    stop("Could not compute an MD5 checksum for ", path, ".")
  }
  list(
    file_name = basename(path),
    size_bytes = unname(information$size),
    modification_time_utc = format(
      information$mtime,
      tz = "UTC",
      usetz = TRUE
    ),
    md5 = checksum
  )
}

write_csv <- function(x, path) {
  write.csv(x, file = path, row.names = FALSE)
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "internal", "r1_real_genotype",
  "real_genotype_sampling.R"
))

default_vcf_path <- file.path(
  workflowr_root,
  "..",
  "genotype-data",
  "YRI_genotype.vcf.gz"
)
default_pair_summary_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "eqtl_summary.rds"
)
vcf_path <- resolve_input_path(
  get_arg("--vcf-path", default_vcf_path),
  workflowr_root,
  "--vcf-path"
)
pair_summary_path <- resolve_input_path(
  get_arg("--pair-summary-path", default_pair_summary_path),
  workflowr_root,
  "--pair-summary-path"
)
n_loci <- as.integer(get_arg("--n-loci", "20"))
variants_per_locus <- as.integer(get_arg("--variants-per-locus", "50"))
maf_min <- as.numeric(get_arg("--maf-min", "0.10"))
candidate_genes_per_seed <- as.integer(get_arg(
  "--candidate-genes-per-seed",
  "120"
))
minimum_raw_variants <- as.integer(get_arg("--minimum-raw-variants", "80"))
n_covariates <- as.integer(get_arg("--n-covariates", "5"))
expression_noise_sd <- as.numeric(get_arg("--noise-sd", "1"))
dynamic_main_effect_sd <- as.numeric(get_arg(
  "--dynamic-main-effect-sd",
  "1"
))
num_basis <- as.integer(get_arg("--num-basis", "20"))
num_cores <- as.integer(get_arg("--num-cores", "4"))
efdr_permutations <- as.integer(get_arg("--efdr-permutations", "100"))
seed_list <- parse_seed_list(get_arg(
  "--seed-list",
  "12345,22345,32345,42345,52345"
))
output_id <- get_arg(
  "--output-id",
  "r1_real_genotype_locus_blocks_pilot5"
)
output_root <- resolve_output_root(
  get_arg(
    "--output-root",
    file.path("output", "revision_simulations", "internal")
  ),
  workflowr_root
)
overwrite <- as_flag(get_arg("--overwrite", "false"))

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
J <- as.integer(n_loci * variants_per_locus)
if (J < 10L || J %% 5L != 0L ||
    n_covariates < 0L || 19L < n_covariates + 3L ||
    !is.finite(maf_min) || maf_min <= 0 || maf_min >= 0.5 ||
    !is.finite(expression_noise_sd) || expression_noise_sd <= 0 ||
    !is.finite(dynamic_main_effect_sd) || dynamic_main_effect_sd <= 0 ||
    num_basis < 2L || num_cores < 1L || efdr_permutations < 1L ||
    !nzchar(output_id)) {
  stop("Invalid real-genotype simulation arguments.")
}

sample_ids <- read_vcf_sample_ids(vcf_path)
if (length(sample_ids) != 19L) {
  stop("The production R1 real-genotype experiment requires exactly 19 VCF samples.")
}

time_grid <- make_time_grid()
class_probs <- c(dynamic_bspline = 0.20, constant = 0.40, zero = 0.40)
expected_class_counts <- exact_proportional_counts(J, class_probs)
scenario <- "r1_real_genotype_locus_blocks_random_bspline"
true_pi0_expected <- unname(class_probs[["constant"]] + class_probs[["zero"]])
fash_methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "FASH-linear-Raw",
  "FASH-linear-BF"
)
direct_methods <- c(
  "FASH-IWP1-Raw",
  "FASH-IWP1-BF",
  "Direct-linear-LRT-eFDR-true-pi0",
  "Direct-quadratic-LRT-eFDR-true-pi0"
)
all_methods <- unique(c(fash_methods, direct_methods))

output_dir <- file.path(output_root, output_id)
replicate_dir <- file.path(output_dir, "replicates")
summary_dir <- file.path(output_dir, "summary")
figure_dir <- file.path(output_dir, "figures")
full_fit_dir <- file.path(output_dir, "full_fits")
genotype_sample_dir <- file.path(output_dir, "genotype_samples")
invisible(lapply(c(
  output_dir,
  replicate_dir,
  summary_dir,
  figure_dir,
  full_fit_dir,
  genotype_sample_dir
), dir.create, recursive = TRUE, showWarnings = FALSE))

vcf_fingerprint <- artifact_fingerprint(vcf_path)
pair_fingerprint <- artifact_fingerprint(pair_summary_path)
configuration <- list(
  format_version = 1L,
  output_id = output_id,
  scenario = scenario,
  J = J,
  n_donors = length(sample_ids),
  sample_ids = sample_ids,
  time_grid = time_grid,
  n_covariates = n_covariates,
  expression_noise_sd = expression_noise_sd,
  dynamic_main_effect_sd = dynamic_main_effect_sd,
  num_basis = num_basis,
  class_probs = class_probs,
  expected_class_counts = expected_class_counts,
  dynamic_amplitude = 2,
  bspline_df = 6,
  bspline_coefficient_sd = 1,
  constant_sd = 1,
  covariate_effect_sd = 0.5,
  intercept_sd = 0,
  apply_t_se_correction = TRUE,
  n_loci = n_loci,
  variants_per_locus = variants_per_locus,
  maf_min = maf_min,
  candidate_genes_per_seed = candidate_genes_per_seed,
  minimum_raw_variants = minimum_raw_variants,
  gene_priority_seed_offset = 7101L,
  block_start_seed_offset = 7201L,
  across_locus_seed_offset = 7401L,
  efdr_permutations = efdr_permutations,
  permutation_seed_rule = "seed + 10000",
  true_pi0 = true_pi0_expected,
  full_fit_seed = seed_list[[1L]],
  seed_list = seed_list,
  vcf_fingerprint = vcf_fingerprint,
  pair_summary_fingerprint = pair_fingerprint
)
configuration_path <- file.path(output_dir, "configuration.rds")
if (file.exists(configuration_path) && !overwrite) {
  cached_configuration <- readRDS(configuration_path)
  if (!isTRUE(all.equal(cached_configuration, configuration))) {
    stop(
      "The existing output ID has different settings or input fingerprints. ",
      "Choose a new --output-id or use --overwrite true."
    )
  }
} else {
  saveRDS(configuration, configuration_path)
}

validate_genotype_samples <- function(x) {
  if (!is.list(x) ||
      !all(c("configuration", "sample_ids", "samples", "extraction_summary") %in% names(x)) ||
      !isTRUE(all.equal(x$configuration, configuration)) ||
      !identical(x$sample_ids, sample_ids) ||
      !identical(names(x$samples), as.character(seed_list))) {
    return(FALSE)
  }
  all(vapply(seq_along(seed_list), function(index) {
    seed <- seed_list[[index]]
    sample <- x$samples[[as.character(seed)]]
    if (!is.list(sample) ||
        !all(c("G", "variant_info", "selected_genes", "ld") %in% names(sample))) {
      return(FALSE)
    }
    block_counts <- table(sample$variant_info$block_index)
    identical(dim(sample$G), c(19L, J)) &&
      identical(rownames(sample$G), sample_ids) &&
      identical(colnames(sample$G), sample$variant_info$unit_id) &&
      nrow(sample$variant_info) == J &&
      length(sample$selected_genes) == n_loci &&
      length(block_counts) == n_loci &&
      all(block_counts == variants_per_locus) &&
      all(sample$variant_info$observed_maf >= maf_min) &&
      all(is.finite(sample$G)) &&
      all(apply(sample$G, 2L, stats::sd) > 0) &&
      nrow(sample$ld$locus) == n_loci &&
      nrow(sample$ld$seed) == 1L &&
      identical(sample$ld$seed$sampling_seed, seed)
  }, logical(1)))
}

genotype_samples_path <- file.path(
  genotype_sample_dir,
  "real_genotype_samples.rds"
)
if (file.exists(genotype_samples_path) && !overwrite) {
  genotype_samples <- readRDS(genotype_samples_path)
  if (!validate_genotype_samples(genotype_samples)) {
    stop("The cached real-genotype samples do not match the requested configuration.")
  }
  message("Reusing cached real-genotype samples: ", genotype_samples_path)
} else {
  message("Reading paper pair identifiers from ", pair_summary_path, ".")
  pair_summary <- readRDS(pair_summary_path)
  pair_ids <- names(pair_summary)
  if (is.null(pair_ids) || length(pair_ids) == 0L) {
    stop("The pair-summary object must be a named list of gene-variant units.")
  }
  rm(pair_summary)
  invisible(gc())
  message("Preparing candidate real-genotype loci with one VCF pass.")
  real_source <- prepare_real_genotype_source(
    pair_ids = pair_ids,
    vcf_path = vcf_path,
    seed_list = seed_list,
    n_loci = n_loci,
    variants_per_locus = variants_per_locus,
    candidate_genes_per_seed = candidate_genes_per_seed,
    minimum_raw_variants = minimum_raw_variants,
    maf_min = maf_min,
    work_dir = tempdir()
  )
  rm(pair_ids)
  invisible(gc())

  samples <- lapply(seed_list, function(seed) {
    message("Sampling real-genotype locus blocks for seed ", seed, ".")
    sample <- sample_real_genotype_blocks(real_source, seed)
    ld <- summarize_real_genotype_ld(sample$G, sample$variant_info)
    ld$correlation <- NULL
    sample$ld <- ld
    sample
  })
  names(samples) <- as.character(seed_list)
  genotype_samples <- list(
    configuration = configuration,
    sample_ids = real_source$sample_ids,
    samples = samples,
    extraction_summary = real_source$extraction_summary
  )
  if (!validate_genotype_samples(genotype_samples)) {
    stop("The newly sampled real-genotype matrices failed validation.")
  }
  saveRDS(genotype_samples, genotype_samples_path)
  rm(real_source)
  invisible(gc())
}

make_replicate <- function(seed) {
  real_sample <- genotype_samples$samples[[as.character(seed)]]
  message("Running R1 real-genotype replicate with seed ", seed, ".")
  out <- run_genotype_level_bspline_eqtl_simulation(
    G = real_sample$G,
    time_grid = time_grid,
    n_covariates = n_covariates,
    class_probs = class_probs,
    expression_noise_sd = expression_noise_sd,
    covariate_effect_sd = 0.5,
    intercept_sd = 0,
    dynamic_amplitude = 2,
    bspline_df = 6,
    bspline_coefficient_sd = 1,
    constant_sd = 1,
    dynamic_main_effect_sd = dynamic_main_effect_sd,
    alpha = 0.05,
    seed = seed,
    num_cores = num_cores,
    num_basis = num_basis,
    exact_class_counts = TRUE,
    apply_t_se_correction = TRUE,
    scenario = scenario,
    output_dir = output_dir,
    save_outputs = FALSE,
    verbose = FALSE
  )

  observed_class_counts <- table(factor(
    out$unit_info$effect_class,
    levels = names(class_probs)
  ))
  if (!identical(as.integer(observed_class_counts), unname(expected_class_counts)) ||
      any(out$eqtl_summary$df != 12) ||
      any(out$eqtl_summary$rank != n_covariates + 2L)) {
    stop("The real-genotype replicate does not match the retained R1 class or regression setting.")
  }
  true_pi0 <- dynamic_null_proportion(out$unit_info, target = "dynamic")
  if (!isTRUE(all.equal(true_pi0, true_pi0_expected))) {
    stop("The simulated dynamic-null proportion does not match the configured mixture.")
  }

  metadata <- real_sample$variant_info
  metadata$genotype_sd <- apply(real_sample$G, 2L, stats::sd)
  out$variant_info <- metadata
  out$unit_info$gene_id <- metadata$gene_id
  out$unit_info$source_variant_id <- metadata$variant_id
  out$unit_info$block_index <- metadata$block_index
  out$unit_info$chromosome <- metadata$chromosome
  out$unit_info$position <- metadata$position
  out$unit_info$observed_maf <- metadata$observed_maf
  attr(out$datasets, "unit_info") <- out$unit_info
  out$settings$genotype_source <- "paper-derived YRI dosage locus blocks"
  out$settings$n_loci <- n_loci
  out$settings$variants_per_locus <- variants_per_locus
  out$settings$maf_min <- maf_min
  out$settings$sample_ids <- sample_ids
  out$settings$vcf_md5 <- vcf_fingerprint$md5
  out$settings$pair_summary_md5 <- pair_fingerprint$md5

  out <- add_direct_interaction_efdr_results_to_genotype_output(
    out,
    n_permutations = efdr_permutations,
    alpha = 0.05,
    seed = seed + 10000L,
    lambda = 0.5,
    pi0_method = "conservative",
    true_pi0 = true_pi0,
    include_true_pi0 = TRUE,
    permute_covariates_with_expression = TRUE,
    num_cores = num_cores,
    overwrite = TRUE,
    verbose = FALSE
  )
  available_methods <- unique(out$result_table$method)
  missing_methods <- setdiff(all_methods, available_methods)
  if (length(missing_methods) > 0L) {
    stop("Missing R1 comparison methods: ", paste(missing_methods, collapse = ", "))
  }

  alpha_curve <- out$alpha_curve[out$alpha_curve$method %in% all_methods, ]
  alpha_curve$seed <- seed
  alpha_005 <- alpha_curve[abs(alpha_curve$alpha - 0.05) < 1e-12, ]
  pi0 <- data.frame(
    seed = seed,
    method = c("FASH-IWP1", "FASH-IWP1", "FASH-linear", "FASH-linear"),
    fit = c("Raw", "BF-corrected", "Raw", "BF-corrected"),
    estimated_pi0 = c(
      constant_component_prior_weight(out$fash_fits$fash_iwp1_raw),
      constant_component_prior_weight(out$fash_fits$fash_iwp1_bf),
      constant_component_prior_weight(out$simplified_fit),
      constant_component_prior_weight(out$simplified_fit_bf)
    ),
    stringsAsFactors = FALSE
  )
  if (identical(seed, seed_list[[1L]])) {
    saveRDS(out, file.path(full_fit_dir, paste0("seed_", seed, ".rds")))
  }
  list(
    configuration = configuration,
    seed = seed,
    permutation_seed = seed + 10000L,
    true_pi0 = true_pi0,
    alpha_curve = alpha_curve,
    alpha_005 = alpha_005,
    pi0 = pi0,
    selected_genes = real_sample$selected_genes,
    variant_info = real_sample$variant_info,
    locus_ld = real_sample$ld$locus,
    seed_ld = real_sample$ld$seed
  )
}

validate_replicate <- function(replicate, seed) {
  required_fields <- c(
    "configuration", "seed", "permutation_seed", "true_pi0",
    "alpha_curve", "alpha_005", "pi0", "selected_genes",
    "variant_info", "locus_ld", "seed_ld"
  )
  if (!is.list(replicate) ||
      !all(required_fields %in% names(replicate)) ||
      !identical(replicate$seed, seed) ||
      !isTRUE(all.equal(replicate$configuration, configuration)) ||
      !isTRUE(all.equal(replicate$true_pi0, true_pi0_expected)) ||
      nrow(replicate$variant_info) != J ||
      length(replicate$selected_genes) != n_loci ||
      nrow(replicate$locus_ld) != n_loci ||
      nrow(replicate$seed_ld) != 1L) {
    return(FALSE)
  }
  observed_methods <- unique(replicate$alpha_curve$method)
  block_counts <- table(replicate$variant_info$block_index)
  all(all_methods %in% observed_methods) &&
    nrow(replicate$alpha_005) == length(all_methods) &&
    nrow(replicate$pi0) == 4L &&
    length(block_counts) == n_loci &&
    all(block_counts == variants_per_locus) &&
    all(replicate$variant_info$observed_maf >= maf_min)
}

replicates <- lapply(seed_list, function(seed) {
  replicate_path <- file.path(replicate_dir, paste0("seed_", seed, ".rds"))
  full_fit_path <- file.path(full_fit_dir, paste0("seed_", seed, ".rds"))
  full_fit_required <- identical(seed, seed_list[[1L]])
  if (file.exists(replicate_path) && !overwrite) {
    cached <- readRDS(replicate_path)
    if (validate_replicate(cached, seed) &&
        (!full_fit_required || file.exists(full_fit_path))) {
      message("Reusing compact replicate cache: ", replicate_path)
      return(cached)
    }
    stop("Cached replicate does not match the requested settings: ", replicate_path)
  }
  replicate <- make_replicate(seed)
  saveRDS(replicate, replicate_path)
  replicate
})

all_alpha <- do.call(rbind, lapply(replicates, `[[`, "alpha_curve"))
all_alpha_005 <- do.call(rbind, lapply(replicates, `[[`, "alpha_005"))
all_pi0 <- do.call(rbind, lapply(replicates, `[[`, "pi0"))
mc_alpha <- summarize_mc_alpha_curves(all_alpha)
mc_alpha_005 <- mc_alpha[abs(mc_alpha$alpha - 0.05) < 1e-12, ]
mc_pi0 <- summarize_mc_pi0(all_pi0)

fash_all_alpha <- all_alpha[all_alpha$method %in% fash_methods, ]
direct_all_alpha <- all_alpha[all_alpha$method %in% direct_methods, ]
fash_mc_alpha <- mc_alpha[mc_alpha$method %in% fash_methods, ]
direct_mc_alpha <- mc_alpha[mc_alpha$method %in% direct_methods, ]
fash_mc_alpha_005 <- mc_alpha_005[mc_alpha_005$method %in% fash_methods, ]
direct_mc_alpha_005 <- mc_alpha_005[mc_alpha_005$method %in% direct_methods, ]
all_variant_metadata <- do.call(rbind, lapply(replicates, `[[`, "variant_info"))
all_locus_ld <- do.call(rbind, lapply(replicates, `[[`, "locus_ld"))
all_seed_ld <- do.call(rbind, lapply(replicates, `[[`, "seed_ld"))
rownames(all_variant_metadata) <- NULL
rownames(all_locus_ld) <- NULL
rownames(all_seed_ld) <- NULL

input_provenance <- data.frame(
  artifact = c("YRI genotype VCF", "Paper gene-variant summary"),
  file_name = c(vcf_fingerprint$file_name, pair_fingerprint$file_name),
  size_bytes = c(vcf_fingerprint$size_bytes, pair_fingerprint$size_bytes),
  modification_time_utc = c(
    vcf_fingerprint$modification_time_utc,
    pair_fingerprint$modification_time_utc
  ),
  md5 = c(vcf_fingerprint$md5, pair_fingerprint$md5),
  stringsAsFactors = FALSE
)

write_csv(all_alpha, file.path(summary_dir, "all_replicate_alpha_curves.csv"))
write_csv(all_alpha_005, file.path(summary_dir, "all_replicate_alpha005.csv"))
write_csv(all_pi0, file.path(summary_dir, "all_replicate_pi0.csv"))
write_csv(mc_alpha, file.path(summary_dir, "mc_alpha_curve.csv"))
write_csv(mc_alpha_005, file.path(summary_dir, "mc_alpha005_summary.csv"))
write_csv(mc_pi0, file.path(summary_dir, "mc_pi0_summary.csv"))
write_csv(fash_all_alpha, file.path(summary_dir, "iwp_vs_linear_fash_replicate_alpha_curves.csv"))
write_csv(fash_mc_alpha, file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha_curve.csv"))
write_csv(fash_mc_alpha_005, file.path(summary_dir, "iwp_vs_linear_fash_mc_alpha005_summary.csv"))
write_csv(direct_all_alpha, file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_replicate_alpha_curves.csv"))
write_csv(direct_mc_alpha, file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_mc_alpha_curve.csv"))
write_csv(direct_mc_alpha_005, file.path(summary_dir, "iwp_fash_vs_direct_true_pi0_mc_alpha005_summary.csv"))
write_csv(all_variant_metadata, file.path(summary_dir, "real_genotype_variant_metadata.csv"))
write_csv(all_locus_ld, file.path(summary_dir, "real_genotype_locus_ld.csv"))
write_csv(all_seed_ld, file.path(summary_dir, "real_genotype_seed_ld.csv"))
write_csv(input_provenance, file.path(summary_dir, "real_genotype_input_provenance.csv"))
write_csv(
  genotype_samples$extraction_summary,
  file.path(summary_dir, "real_genotype_extraction_summary.csv")
)

plot_subtitle <- paste0(
  length(seed_list), " seeds; 20 paper-derived loci x 50 variants; N = 19, J = ", J
)
if (n_loci != 20L || variants_per_locus != 50L) {
  plot_subtitle <- paste0(
    length(seed_list), " seeds; ", n_loci, " paper-derived loci x ",
    variants_per_locus, " variants; N = 19, J = ", J
  )
}
plot_mc_alpha_curves(
  fash_mc_alpha,
  metric = "power",
  file = file.path(figure_dir, "iwp_vs_linear_fash_mc_power.png"),
  title = "Real-genotype R1: IWP versus linear FASH power",
  subtitle = plot_subtitle,
  style_profile = "combined"
)
plot_mc_alpha_curves(
  fash_mc_alpha,
  metric = "fdr",
  file = file.path(figure_dir, "iwp_vs_linear_fash_mc_fdr.png"),
  title = "Real-genotype R1: IWP versus linear FASH FDR estimate",
  subtitle = plot_subtitle,
  legend_position = "topleft",
  style_profile = "combined"
)
plot_mc_alpha_curves(
  direct_mc_alpha,
  metric = "power",
  file = file.path(figure_dir, "iwp_fash_vs_direct_true_pi0_mc_power.png"),
  title = "Real-genotype R1: IWP FASH versus direct interaction power",
  subtitle = paste0(plot_subtitle, "; direct eFDR uses true pi0 = 0.8"),
  style_profile = "combined"
)
plot_mc_alpha_curves(
  direct_mc_alpha,
  metric = "fdr",
  file = file.path(figure_dir, "iwp_fash_vs_direct_true_pi0_mc_fdr.png"),
  title = "Real-genotype R1: IWP FASH versus direct interaction FDR estimate",
  subtitle = paste0(plot_subtitle, "; direct eFDR uses true pi0 = 0.8"),
  legend_position = "bottomright",
  style_profile = "combined"
)

print(fash_mc_alpha_005)
print(direct_mc_alpha_005)
print(mc_pi0)
print(all_seed_ld)
