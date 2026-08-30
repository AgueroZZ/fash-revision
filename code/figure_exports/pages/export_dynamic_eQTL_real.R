script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
source(file.path(repo_root, "code", "figure_exports", "export_page_figures.R"))

export_page_figures(
  page_rmd = "analysis/dynamic_eQTL_real.rmd",
  page_html = "docs/dynamic_eQTL_real.html",
  final_chunk = "switch_examples",
  repo_root = repo_root
)
