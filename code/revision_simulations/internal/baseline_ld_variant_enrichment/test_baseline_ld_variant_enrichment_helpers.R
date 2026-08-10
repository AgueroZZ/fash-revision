#!/usr/bin/env Rscript

# Focused tests for the internal baselineLD variant enrichment helpers.

find_workflowr_root <- function(start = getwd()) {
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

workflowr_root <- find_workflowr_root()
helper_path <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "baseline_ld_variant_enrichment",
  "baseline_ld_variant_enrichment_helpers.R"
)
source(helper_path)

pair_keys <- c("gene1_rs1", "gene1_rs2", "gene2_rs3", "gene3_rs4")
lfdr <- c(0.01, 0.02, 0.03, 0.90)
linear <- data.frame(
  rs_id = c("rs1", "rs2", "rs3", "rs4"),
  ensamble_id = c("gene1", "gene1", "gene2", "gene3"),
  pvalue = c(0.01, 0.20, 0.30, 0.40),
  eFDR = c(0.01, 0.20, 0.30, 0.40),
  stringsAsFactors = FALSE
)
nonlinear <- linear
nonlinear$eFDR <- c(0.20, 0.20, 0.01, 0.40)
sets <- derive_requested_fash_sets(
  pair_keys,
  lfdr,
  linear,
  nonlinear,
  alpha = 0.05
)
stopifnot(
  identical(sets$selected_sets$all_fash, c("rs1", "rs2", "rs3")),
  identical(sets$selected_sets$one_lead_fash_per_gene, c("rs1", "rs3")),
  identical(sets$selected_sets$fash_only_pair_variants, "rs2"),
  nrow(sets$fash_only_pairs) == 1L
)

toy_path <- tempfile(fileext = ".annot.gz")
connection <- gzfile(toy_path, open = "wt")
writeLines(c(
  paste(
    c("CHR", "BP", "SNP", "CM", "Binary_A", "Binary_B", "Continuous_C"),
    collapse = "\t"
  ),
  "1\t100\trs1\t0\t1\t0\t0.25",
  "1\t200\trs2\t0\t0\t1\t0.75",
  "1\t300\trs9\t0\t1\t1\t1.25"
), connection)
close(connection)

toy_table <- data.table::fread(toy_path, showProgress = FALSE)
classification <- classify_baseline_ld_columns(toy_table)
stopifnot(
  identical(
    classification$annotation_type,
    c("binary", "binary", "continuous")
  )
)
toy_variants <- data.frame(
  variant_id = c("rs1", "rs2", "rs3"),
  chromosome = c("1", "1", "2"),
  position = c(100L, 200L, 500L),
  stringsAsFactors = FALSE
)
toy_matrix <- build_baseline_ld_matrix(
  paths = c("1" = toy_path),
  variant_table = toy_variants,
  binary_annotation_columns = c("Binary_A", "Binary_B")
)
stopifnot(
  identical(toy_matrix$variant_id, c("rs1", "rs2")),
  is.logical(toy_matrix$Binary_A),
  identical(toy_matrix$Binary_A, c(TRUE, FALSE)),
  identical(toy_matrix$Binary_B, c(FALSE, TRUE))
)
unlink(toy_path)

coverage <- summarize_set_coverage(
  selected_sets = list(first = c("rs1", "rs3"), second = "rs2"),
  covered_ids = c("rs1", "rs2"),
  universe_ids = c("rs1", "rs2", "rs3")
)
stopifnot(
  coverage$covered_count[coverage$discovery_set == "first"] == 1L,
  coverage$coverage_proportion[coverage$discovery_set == "second"] == 1
)

fit_path <- file.path(
  workflowr_root,
  "output",
  "dynamic_eQTL_real",
  "fash_fit1_update.RData"
)
linear_path <- file.path(
  workflowr_root,
  "data",
  "dynamic_eQTL_real",
  "strober_linear",
  "linear_dynamic_eqtls_5_pc.txt"
)
nonlinear_path <- file.path(
  workflowr_root,
  "data",
  "dynamic_eQTL_real",
  "strober_nonlinear",
  "non_linear_dynamic_eqtls_5_pc.txt"
)
fit_environment <- new.env(parent = emptyenv())
loaded_names <- load(fit_path, envir = fit_environment)
stopifnot(identical(loaded_names, "fash_fit1_update"))
fit <- fit_environment$fash_fit1_update
real_sets <- derive_requested_fash_sets(
  pair_keys = names(fit$fash_data$data_list),
  lfdr = fit$lfdr,
  linear_results = data.table::fread(linear_path, data.table = FALSE),
  nonlinear_results = data.table::fread(nonlinear_path, data.table = FALSE),
  alpha = 0.05
)
stopifnot(
  identical(real_sets$pair_summary$pair_count, c(9205L, 1177L, 8062L)),
  identical(
    real_sets$pair_summary$unique_variant_count,
    c(9139L, 1170L, 8001L)
  ),
  real_sets$pair_summary$unique_gene_count[1L] == 1177L,
  real_sets$pair_summary$unique_gene_count[3L] == 1119L
)

message("baselineLD variant enrichment helper tests passed.")
