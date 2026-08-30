script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
source(file.path(repo_root, "code", "figure_exports", "export_page_figures.R"))

export_page_figures(
  page_rmd = "analysis/revision_internal_matched_fash_linear_real_data.rmd",
  page_html = "docs/revision_internal_matched_fash_linear_real_data.html",
  final_chunk = "venn-variants",
  repo_root = repo_root
)
