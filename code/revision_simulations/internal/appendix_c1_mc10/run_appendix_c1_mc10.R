# Run from the workflowr repository root. The author authorized four threads.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
source("code/revision_simulations/internal/appendix_c1_mc10/mc10_helpers.R")
contract <- c1_contract()
workers <- as.integer(Sys.getenv("FASH_C1_WORKERS", "4"))
stopifnot(workers >= 1L, workers <= 4L)
directory <- file.path("output/revision_simulations/internal", contract$result_id)
dir.create(directory, recursive = TRUE, showWarnings = FALSE)
source_paths <- c(helper = "code/revision_simulations/appendix_b/appendix_b_helpers.R",
                 runner = "code/revision_simulations/internal/appendix_c1_mc10/run_appendix_c1_mc10.R",
                 extension = "code/revision_simulations/internal/appendix_c1_mc10/mc10_helpers.R",
                 integration = "code/revision_simulations/internal/appendix_c1_mc10/gaussian_likelihood.R")
equivalence <- readRDS("output/revision_simulations/internal/appendix_c1_mc10_validation/equivalence.rds")
stopifnot(identical(attr(equivalence, "integration_sha256"), c1_sha(source_paths[["integration"]])))
configuration <- list(contract = contract, package = require_appendix_b_fashr(),
                      source_sha256 = vapply(source_paths, c1_sha, character(1)))
config_path <- file.path(directory, "configuration.rds")
if (file.exists(config_path)) stopifnot(identical(readRDS(config_path), configuration)) else
  c1_write(configuration, config_path)
if (file.exists(file.path(directory, "complete.flag"))) stop("The cache is already complete.")

cluster <- parallel::makePSOCKcluster(workers, outfile = file.path(directory, "workers.log"))
root <- normalizePath(".")
parallel::clusterCall(cluster, function(root) {
  setwd(root)
  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
             VECLIB_MAXIMUM_THREADS = "1", MKL_NUM_THREADS = "1")
  source("code/revision_simulations/internal/appendix_c1_mc10/mc10_helpers.R", local = .GlobalEnv)
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1L)
  require_appendix_b_fashr()
}, root)
cat("Starting ten seeds with", workers, "single-thread workers.\n"); flush.console()
outputs <- tryCatch(parallel::parLapplyLB(cluster, contract$seeds, c1_run_seed,
                                         directory = directory, configuration = configuration),
                    finally = parallel::stopCluster(cluster))
records <- lapply(outputs, readRDS)
grid <- do.call(rbind, lapply(records, `[[`, "grid"))
alpha <- do.call(rbind, lapply(records, `[[`, "alpha"))
rownames(grid) <- rownames(alpha) <- NULL
c1_validate(grid, alpha)
grid_long <- c1_grid_long(grid)
grid_summary <- c1_summary(grid_long, c("setting", "rho_dynamic", "order", "stage", "true_pi0"), "pi0")
alpha_summary <- c1_summary(alpha, c("order", "stage", "alpha"),
                             c("discoveries", "realized_fdp", "power"))
names(alpha_summary) <- sub("realized_fdp", "empirical_fdr", names(alpha_summary), fixed = TRUE)
tables <- list(grid_by_seed = grid, grid_long = grid_long, grid_summary = grid_summary,
               alpha_by_seed = alpha, alpha_summary = alpha_summary)
c1_write(tables, file.path(directory, "summary_tables.rds"))
for (name in names(tables)) write.csv(tables[[name]], file.path(directory, paste0(name, ".csv")), row.names = FALSE)
artifacts <- c("configuration.rds", "summary_tables.rds", "focused_seed12345.rds",
               paste0(names(tables), ".csv"), file.path("replicates", basename(unlist(outputs))))
manifest <- list(schema = contract$schema, configuration = configuration,
                 completed_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
                 artifact_sha256 = setNames(vapply(file.path(directory, artifacts), c1_sha, character(1)), artifacts))
c1_write(manifest, file.path(directory, "manifest.rds"))
writeLines(c("seeds=10", "grid_rows=920", "focused_alpha_rows=1600", "workers=4"),
           file.path(directory, "complete.flag"))
cat("Completed and validated", directory, "\n")
