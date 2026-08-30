script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."))
source(file.path(repo_root, "code", "figure_exports", "export_page_figures.R"))

export_page_figures(
  page_rmd = "analysis/nonlinear_dynamic_eQTL_real.rmd",
  page_html = "docs/nonlinear_dynamic_eQTL_real.html",
  final_chunk = "switch_examples",
  repo_root = repo_root,
  offline_biomart_cache = "output/dynamic_eQTL_real/cache_gene_map.rds"
)
