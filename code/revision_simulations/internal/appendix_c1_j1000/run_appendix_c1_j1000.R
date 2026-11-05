# Run sequentially from the workflowr repository root; no historical cache edits.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("code/revision_simulations/internal/appendix_c1_j1000/helpers.R")
package <- require_appendix_b_fashr()
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1L)
contract <- c1j_contract()
directory <- file.path("output/revision_simulations/internal", contract$result_id)
dir.create(file.path(directory, "replicates"), recursive = TRUE, showWarnings = FALSE)
equivalence <- readRDS("output/revision_simulations/internal/appendix_c1_mc10_validation/equivalence.rds")
stopifnot(identical(attr(equivalence, "integration_sha256"), c1_sha(c1j_source_paths()[["integration"]])))
grid_report <- c1_load_report()
configuration <- list(contract = contract, package = package,
  source_sha256 = vapply(c1j_source_paths(), c1_sha, character(1)),
  grid_manifest_sha256 = c1_sha(file.path(grid_report$directory, "manifest.rds")))
config_path <- file.path(directory, "configuration.rds")
if (file.exists(config_path)) stopifnot(identical(readRDS(config_path), configuration)) else
  c1_write(configuration, config_path)
if (file.exists(file.path(directory, "complete.flag"))) stop("The cache is already complete.")
token <- digest::digest(configuration, algo = "sha256")
records <- vector("list", length(contract$seeds))
for (i in seq_along(contract$seeds)) {
  seed <- contract$seeds[[i]]
  path <- file.path(directory, "replicates", paste0("seed_", seed, ".rds"))
  if (file.exists(path)) {
    record <- readRDS(path)
    stopifnot(identical(record$configuration_token, token))
  } else {
    started <- proc.time()[["elapsed"]]
    set.seed(seed)
    data <- build_appendix_b_datasets(contract$focused_J, 0.2, 0.1, contract$sigma)
    fit <- c1_fit_models(data, sort(c(0, exp(-0.5 * seq(0, 10, by = 0.2)))),
                        contract$focused_penalty, contract$num_basis, 1L)
    alpha <- focused_alpha_curve(fit, contract$alpha)
    alpha$seed <- seed
    alpha$J <- contract$focused_J
    record <- list(configuration_token = token, seed = seed, alpha = alpha,
      data_bundle = data, fit_bundle = fit,
      input_sha256 = digest::digest(data, algo = "sha256"),
      warnings = unique(c(fit$order1$warnings, fit$order2$warnings)),
      elapsed_seconds = proc.time()[["elapsed"]] - started)
    c1_write(record, path)
  }
  if (i == 1L) c1_write(record[c("data_bundle", "fit_bundle")],
                        file.path(directory, "focused_seed12345.rds"))
  records[[i]] <- record$alpha
  cat("Completed seed", seed, "in", record$elapsed_seconds, "seconds; warnings:",
      paste(record$warnings, collapse = " | "), "\n")
  flush.console()
}
alpha <- do.call(rbind, records)
rownames(alpha) <- NULL
c1j_validate(alpha)
summary <- c1j_summarize(alpha)
c1_write(list(alpha_by_seed = alpha, alpha_summary = summary), file.path(directory, "summary_tables.rds"))
write.csv(alpha, file.path(directory, "alpha_by_seed.csv"), row.names = FALSE)
write.csv(summary, file.path(directory, "alpha_summary.csv"), row.names = FALSE)
artifacts <- c("configuration.rds", "focused_seed12345.rds", "summary_tables.rds",
  "alpha_by_seed.csv", "alpha_summary.csv", paste0("replicates/seed_", contract$seeds, ".rds"))
manifest <- list(schema = contract$schema, configuration = configuration,
  completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
  artifact_sha256 = setNames(vapply(file.path(directory, artifacts), c1_sha, character(1)), artifacts))
c1_write(manifest, file.path(directory, "manifest.rds"))
writeLines(c("seeds=10", "J=1000", "focused_alpha_rows=1600", "workers=1"),
           file.path(directory, "complete.flag"))
print(summary[abs(summary$alpha - 0.05) < 1e-12, ])
cat("Completed and validated", directory, "\n")
