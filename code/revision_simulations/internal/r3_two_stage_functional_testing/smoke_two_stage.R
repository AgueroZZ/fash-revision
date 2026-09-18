#!/usr/bin/env Rscript

# A reduced integration check; these outputs are never production results.
experiment <- "code/revision_simulations/internal/r3_two_stage_functional_testing"
source(file.path(experiment, "two_stage_helpers.R"))
source(file.path(experiment, "reconstruct_r3.R"))
source("code/revision_simulations/r3_r4_fashr0143/source_snapshots/r3_full_universe_functional_simulation_functions.R")
source("code/revision_simulations/shared/real_genotype_one_per_gene.R")
source("code/revision_simulations/internal/r3_ideal_gaussian_measurement/ideal_gaussian_measurement.R")
contract <- r3ts_contract()
config <- readRDS(file.path("output/revision_simulations/mc", contract$formal_id, "configuration.rds"))
genotypes <- readRDS("output/revision_simulations/shared/real_genotype_one_per_gene_J6362_pilot5/genotype_samples.rds")
sample <- genotypes$samples[["12345"]]
keep <- seq_len(60L)
sample$G <- sample$G[, keep, drop = FALSE]
sample$selection <- sample$selection[keep, , drop = FALSE]
sample$variant_info <- sample$variant_info[keep, , drop = FALSE]
rows <- list()
for (mechanism in contract$mechanisms) {
  started <- proc.time()[["elapsed"]]
  inputs <- r3ts_reconstruct(config, sample, 12345L, mechanism)
  stopifnot(nrow(inputs$effects$true_functionals) == 60L)
  # Preserve all fitting and posterior settings; only the smoke sample is smaller.
  fits <- r3ts_fit_all(inputs, config, num_cores = 2L)
  result <- r3ts_infer(fits$fash_iwp1_bf, inputs, config, num_cores = 1L)
  stopifnot(nrow(result$fdr_table) == 60L,
            all(result$stage2$index %in% result$candidates),
            nrow(result$stage2) == 4L * length(result$candidates),
            all(is.finite(result$stage2$lfsr)),
            all(result$metrics$functional$true_target_total > 0))
  rows[[mechanism]] <- data.frame(
    mechanism = mechanism, n_units = 60L,
    dynamic_candidates = length(result$candidates),
    functional_candidates = nrow(result$stage2),
    elapsed_seconds = proc.time()[["elapsed"]] - started,
    scope = "reduced integration smoke; not scientific results"
  )
  rm(inputs, fits, result)
  invisible(gc())
}
output <- "output/revision_simulations/internal/r3_two_stage_development_20260904"
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, rows), file.path(output, "smoke_validation.csv"), row.names = FALSE)
cat("Both reduced two-stage integration checks passed.\n")
