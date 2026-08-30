script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
source(file.path(repo_root, "code", "figure_exports", "export_page_figures.R"))

export_page_figures(
  page_rmd = "analysis/revision_internal_evaluation_grid_mc_sensitivity.rmd",
  page_html = "docs/revision_internal_evaluation_grid_mc_sensitivity.html",
  final_chunk = "M3000-vs-M5000",
  repo_root = repo_root,
  rasterize_point_chunks = c(
    "grid-0p10-vs-0p15",
    "grid-0p10-vs-0p05",
    "M2000-vs-M3000",
    "M3000-vs-M5000"
  )
)
