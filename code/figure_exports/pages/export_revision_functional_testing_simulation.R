script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
source(file.path(repo_root, "code", "figure_exports", "export_page_figures.R"))

export_page_figures(
  page_rmd = "analysis/revision_functional_testing_simulation.rmd",
  page_html = "docs/revision_functional_testing_simulation.html",
  final_chunk = "r3b-fsr",
  repo_root = repo_root
)
