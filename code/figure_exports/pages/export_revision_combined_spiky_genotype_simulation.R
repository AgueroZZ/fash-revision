script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
source(file.path(repo_root, "code", "figure_exports", "export_page_figures.R"))

export_page_figures(
  page_rmd = "analysis/revision_combined_spiky_genotype_simulation.rmd",
  page_html = "docs/revision_combined_spiky_genotype_simulation.html",
  final_chunk = "peak-power",
  repo_root = repo_root
)
