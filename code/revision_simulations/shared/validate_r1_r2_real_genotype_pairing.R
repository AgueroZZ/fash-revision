#!/usr/bin/env Rscript

# Validate that formal R1, R2, R3A, and R3B used the same real-genotype
# selection. Method-specific settings remain page-specific.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/real_genotype_one_per_gene.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/real_genotype_one_per_gene.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

get_arg <- function(name, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  hit <- which(args == name)
  if (length(hit) == 0L || hit[[1L]] == length(args)) return(default)
  args[hit[[1L]] + 1L]
}

workflowr_root <- find_workflowr_root()
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "simulation_functions.R"
))
source(file.path(
  workflowr_root,
  "code", "revision_simulations", "shared", "real_genotype_one_per_gene.R"
))

r1_output_id <- get_arg(
  "--r1-output-id",
  paste0(
    "r1_real_genotype_one_per_gene_J6362_",
    "random_bspline_main_effect_",
    "linear_mixture_predstep1_penalty10_pilot5"
  )
)
r2_output_id <- get_arg(
  "--r2-output-id",
  paste0(
    "r2_real_genotype_one_per_gene_J6362_",
    "timed_cosine_one_two_three_peak_main_effect_",
    "linear_mixture_predstep1_penalty10_pilot5"
  )
)
r3_output_id <- get_arg(
  "--r3-output-id",
  paste0(
    "r3_real_genotype_one_per_gene_J6362_",
    "matched_functional_relative_clearance_main_effect_pilot5"
  )
)
genotype_cache_path <- get_arg(
  "--genotype-cache",
  file.path(
    workflowr_root,
    "output", "revision_simulations", "shared",
    "real_genotype_one_per_gene_J6362_pilot5",
    "genotype_samples.rds"
  )
)
if (!file.exists(genotype_cache_path)) {
  stop("The shared genotype cache is missing: ", genotype_cache_path)
}

r1_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "mc", r1_output_id
)
r2_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "mc", r2_output_id
)
r3_dir <- file.path(
  workflowr_root,
  "output", "revision_simulations", "mc", r3_output_id
)
r1_configuration_path <- file.path(r1_dir, "configuration.rds")
r2_configuration_path <- file.path(r2_dir, "configuration.rds")
r3_configuration_path <- file.path(r3_dir, "configuration.rds")
required_paths <- c(
  r1_configuration_path,
  r2_configuration_path,
  r3_configuration_path
)
if (any(!file.exists(required_paths))) {
  stop("Formal R1/R2/R3 configuration caches are incomplete.")
}

genotype_cache <- readRDS(genotype_cache_path)
r1_configuration <- readRDS(r1_configuration_path)
r2_configuration <- readRDS(r2_configuration_path)
r3_configuration <- readRDS(r3_configuration_path)
seed_list <- genotype_cache$configuration$seed_list

make_check <- function(check,
                       seed = NA_integer_,
                       mismatched_items = 0L,
                       max_absolute_numeric_difference = 0,
                       detail = "") {
  data.frame(
    check = check,
    seed = seed,
    mismatched_items = as.integer(mismatched_items),
    max_absolute_numeric_difference = as.numeric(max_absolute_numeric_difference),
    status = if (
      isTRUE(mismatched_items == 0L) &&
        isTRUE(max_absolute_numeric_difference == 0)
    ) "PASS" else "FAIL",
    detail = detail,
    stringsAsFactors = FALSE
  )
}

checks <- list()
add_check <- function(x) {
  checks[[length(checks) + 1L]] <<- x
}

configuration_fields <- c(
  "J", "n_donors", "n_covariates", "time_grid", "class_probs",
  "expected_class_counts", "linear_prior_mode", "common_sd_grid",
  "common_pred_step", "common_penalty", "genotype_source",
  "genotype_selection_rule", "genotype_repeated_variant_rule",
  "genotype_maf_min", "genotype_sample_ids", "genotype_cache_fingerprint",
  "genotype_source_configuration", "seed_list"
)
configuration_mismatches <- sum(vapply(configuration_fields, function(field) {
  !isTRUE(all.equal(r1_configuration[[field]], r2_configuration[[field]]))
}, logical(1)))
add_check(make_check(
  "R1_R2_shared_configuration",
  mismatched_items = configuration_mismatches,
  detail = "Shared genotype, grid, penalty, class, donor, and seed settings"
))

cross_page_configuration_fields <- c(
  "J", "n_donors", "n_covariates", "time_grid", "class_probs",
  "expected_class_counts", "genotype_source", "genotype_selection_rule",
  "genotype_repeated_variant_rule", "genotype_maf_min",
  "genotype_sample_ids", "genotype_cache_fingerprint",
  "genotype_source_configuration", "seed_list"
)
cross_page_configuration_mismatches <- sum(vapply(
  cross_page_configuration_fields,
  function(field) {
    !isTRUE(all.equal(r1_configuration[[field]], r3_configuration[[field]])) ||
      !isTRUE(all.equal(r2_configuration[[field]], r3_configuration[[field]]))
  },
  logical(1)
))
add_check(make_check(
  "R1_R2_R3_shared_genotype_configuration",
  mismatched_items = cross_page_configuration_mismatches,
  detail = "Shared genotype, class, donor, and seed settings"
))

r3_mechanism_mismatches <- sum(c(
  !identical(
    r3_configuration$truth_mechanisms,
    c("random_bspline", "raised_cosine")
  ),
  !identical(
    as.integer(r3_configuration$expected_truth_group_counts),
    rep(212L, 6L)
  )
))
add_check(make_check(
  "R3_functional_truth_configuration",
  mismatched_items = r3_mechanism_mismatches,
  detail = "R3A/R3B mechanisms and six exact dynamic truth cells"
))

actual_cache_fingerprint <- artifact_fingerprint(genotype_cache_path)
fingerprint_mismatches <- sum(c(
  !identical(r1_configuration$genotype_cache_fingerprint, actual_cache_fingerprint),
  !identical(r2_configuration$genotype_cache_fingerprint, actual_cache_fingerprint),
  !identical(r3_configuration$genotype_cache_fingerprint, actual_cache_fingerprint)
))
add_check(make_check(
  "shared_genotype_cache_fingerprint",
  mismatched_items = fingerprint_mismatches,
  detail = actual_cache_fingerprint$md5
))

for (seed in seed_list) {
  r1_path <- file.path(r1_dir, "replicates", paste0("seed_", seed, ".rds"))
  r2_path <- file.path(r2_dir, "replicates", paste0("seed_", seed, ".rds"))
  r3a_path <- file.path(
    r3_dir,
    "replicates",
    paste0("random_bspline_seed_", seed, ".rds")
  )
  r3b_path <- file.path(
    r3_dir,
    "replicates",
    paste0("raised_cosine_seed_", seed, ".rds")
  )
  if (any(!file.exists(c(r1_path, r2_path, r3a_path, r3b_path)))) {
    stop("Formal compact replicate caches are incomplete for seed ", seed, ".")
  }
  r1 <- readRDS(r1_path)
  r2 <- readRDS(r2_path)
  r3a <- readRDS(r3a_path)
  r3b <- readRDS(r3b_path)
  shared <- genotype_cache$samples[[as.character(seed)]]

  digest_mismatches <- sum(c(
    !identical(r1$genotype_digest, r2$genotype_digest),
    !identical(r1$genotype_digest, r3a$genotype_digest),
    !identical(r1$genotype_digest, r3b$genotype_digest),
    !identical(r1$genotype_digest, shared$genotype_digest),
    !identical(r2$genotype_digest, shared$genotype_digest),
    !identical(r3a$genotype_digest, shared$genotype_digest),
    !identical(r3b$genotype_digest, shared$genotype_digest)
  ))
  add_check(make_check(
    "per_seed_genotype_digest",
    seed = seed,
    mismatched_items = digest_mismatches,
    detail = shared$genotype_digest
  ))

  r1_keys <- r1$selected_pair_keys
  r2_keys <- r2$selected_pair_keys
  r3a_keys <- r3a$selected_pair_keys
  r3b_keys <- r3b$selected_pair_keys
  shared_keys <- shared$selection$pair_key
  key_mismatches <- sum(
    r1_keys != shared_keys |
      r2_keys != shared_keys |
      r3a_keys != shared_keys |
      r3b_keys != shared_keys |
      r1_keys != r2_keys |
      r1_keys != r3a_keys |
      r1_keys != r3b_keys
  )
  add_check(make_check(
    "per_seed_selected_pair_keys",
    seed = seed,
    mismatched_items = key_mismatches,
    detail = paste(length(shared_keys), "ordered gene-variant keys")
  ))

  balance_columns <- c("effect_class", "n", "maf_mean")
  r1_balance <- r1$truth_maf_balance[
    order(r1$truth_maf_balance$effect_class),
    balance_columns,
    drop = FALSE
  ]
  r2_balance <- r2$truth_maf_balance[
    order(r2$truth_maf_balance$effect_class),
    balance_columns,
    drop = FALSE
  ]
  r3a_balance <- r3a$truth_maf_balance[
    order(r3a$truth_maf_balance$effect_class),
    balance_columns,
    drop = FALSE
  ]
  r3b_balance <- r3b$truth_maf_balance[
    order(r3b$truth_maf_balance$effect_class),
    balance_columns,
    drop = FALSE
  ]
  balance_difference <- max(abs(c(
    r1_balance$maf_mean - r2_balance$maf_mean,
    r1_balance$maf_mean - r3a_balance$maf_mean,
    r1_balance$maf_mean - r3b_balance$maf_mean
  )))
  balance_mismatches <- sum(
    r1_balance$effect_class != r2_balance$effect_class |
      r1_balance$effect_class != r3a_balance$effect_class |
      r1_balance$effect_class != r3b_balance$effect_class |
      r1_balance$n != r2_balance$n |
      r1_balance$n != r3a_balance$n |
      r1_balance$n != r3b_balance$n
  )
  add_check(make_check(
    "per_seed_truth_class_maf_balance",
    seed = seed,
    mismatched_items = balance_mismatches,
    max_absolute_numeric_difference = balance_difference,
    detail = "R1, R2, R3A, and R3B share MAF-balanced class assignments"
  ))

  truth_group_mismatches <- sum(c(
    nrow(r3a$truth_group_counts) != 6L,
    nrow(r3b$truth_group_counts) != 6L,
    !identical(as.integer(r3a$truth_group_counts$n_dynamic), rep(212L, 6L)),
    !identical(as.integer(r3b$truth_group_counts$n_dynamic), rep(212L, 6L))
  ))
  add_check(make_check(
    "per_seed_R3_functional_truth_counts",
    seed = seed,
    mismatched_items = truth_group_mismatches,
    detail = "Six R3A/R3B dynamic truth cells with 212 units each"
  ))
}

full_seed <- seed_list[[1L]]
r1_full_path <- file.path(r1_dir, "full_fits", paste0("seed_", full_seed, ".rds"))
r2_full_path <- file.path(r2_dir, "full_fits", paste0("seed_", full_seed, ".rds"))
if (!file.exists(r1_full_path) || !file.exists(r2_full_path)) {
  stop("The seed-", full_seed, " full-fit caches are incomplete.")
}
r1_full <- readRDS(r1_full_path)
r2_full <- readRDS(r2_full_path)
shared_full <- genotype_cache$samples[[as.character(full_seed)]]
full_difference <- max(c(
  abs(r1_full$genotype - shared_full$G),
  abs(r2_full$genotype - shared_full$G),
  abs(r1_full$genotype - r2_full$genotype)
))
full_id_mismatches <- sum(c(
  !identical(rownames(r1_full$genotype), genotype_cache$sample_ids),
  !identical(rownames(r2_full$genotype), genotype_cache$sample_ids),
  !identical(colnames(r1_full$genotype), shared_full$selection$pair_key),
  !identical(colnames(r2_full$genotype), shared_full$selection$pair_key)
))
add_check(make_check(
  "seed_12345_full_fit_genotype_matrix",
  seed = full_seed,
  mismatched_items = full_id_mismatches,
  max_absolute_numeric_difference = full_difference,
  detail = "Full-fit genotype values and ordered identifiers"
))

validation <- do.call(rbind, checks)
rownames(validation) <- NULL
if (any(validation$status != "PASS")) {
  print(validation)
  stop("The formal R1/R2/R3 real-genotype pairing audit failed.")
}

r1_summary_dir <- file.path(r1_dir, "summary")
r2_summary_dir <- file.path(r2_dir, "summary")
r3_summary_dir <- file.path(r3_dir, "summary")
dir.create(r1_summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(r2_summary_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(r3_summary_dir, recursive = TRUE, showWarnings = FALSE)
r1_output_path <- file.path(
  r1_summary_dir,
  "real_genotype_pairing_validation.csv"
)
r2_output_path <- file.path(
  r2_summary_dir,
  "real_genotype_pairing_validation.csv"
)
r3_output_path <- file.path(
  r3_summary_dir,
  "real_genotype_pairing_validation.csv"
)
write.csv(validation, r1_output_path, row.names = FALSE)
write.csv(validation, r2_output_path, row.names = FALSE)
write.csv(validation, r3_output_path, row.names = FALSE)
validation_md5 <- unname(tools::md5sum(c(
  r1_output_path,
  r2_output_path,
  r3_output_path
)))
if (length(unique(validation_md5)) != 1L) {
  stop("The R1, R2, and R3 pairing-validation copies differ.")
}

message("Formal R1/R2/R3 real-genotype pairing validation passed.")
print(validation)
