# Read completed J = 1000 focused results and the unchanged J = 1000 grid cache.
source("code/revision_simulations/internal/appendix_c1_j1000/helpers.R")

c1j_load_report <- function(root = ".") {
  contract <- c1j_contract()
  directory <- file.path(root, "output/revision_simulations/internal", contract$result_id)
  stopifnot(file.exists(file.path(directory, "complete.flag")))
  manifest <- readRDS(file.path(directory, "manifest.rds"))
  stopifnot(identical(manifest$schema, contract$schema),
    identical(manifest$configuration$contract, contract),
    identical(unname(vapply(file.path(directory, names(manifest$artifact_sha256)), c1_sha, character(1))),
              unname(manifest$artifact_sha256)),
    identical(vapply(setNames(file.path(root, c1j_source_paths()), names(c1j_source_paths())), c1_sha, character(1)),
              manifest$configuration$source_sha256))
  original <- c1_load_report(root)
  stopifnot(identical(c1_sha(file.path(original$directory, "manifest.rds")),
                      manifest$configuration$grid_manifest_sha256))
  focused <- readRDS(file.path(directory, "focused_seed12345.rds"))
  stopifnot(length(focused$data_bundle$data_list) == 1000L,
    identical(unname(focused$data_bundle$class_counts), c(800L, 100L, 100L)))
  tables <- original$tables
  replacement <- readRDS(file.path(directory, "summary_tables.rds"))
  c1j_validate(replacement$alpha_by_seed, contract)
  stopifnot(isTRUE(all.equal(c1j_summarize(replacement$alpha_by_seed), replacement$alpha_summary)))
  tables$alpha_by_seed <- replacement$alpha_by_seed
  tables$alpha_summary <- replacement$alpha_summary
  list(directory = directory, manifest = manifest, contract = contract,
       tables = tables, focused = focused)
}

c1j_load_display <- function(report) {
  path <- file.path(report$directory, "posterior_display.rds")
  stopifnot(identical(c1_sha(path), readLines(file.path(report$directory, "posterior_display.sha256"))))
  display <- readRDS(path)
  stopifnot(identical(display$focused_sha256, report$manifest$artifact_sha256[["focused_seed12345.rds"]]),
    identical(display$ids, c1j_select_examples(report$focused)),
    nrow(display$predictions) == 604L, nrow(display$observations) == 64L)
  display
}
