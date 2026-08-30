script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_argument[[1L]]))
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."))
source(file.path(repo_root, "code", "figure_exports", "export_page_figures.R"))

temporary_html <- tempfile(fileext = ".html")
writeLines(
  c(
    '<img src="figure/example.rmd/first-1.png">',
    '<img src="figure/example.rmd/with%20space-1.png">',
    '<img src="figure/example.rmd/first-1.png">'
  ),
  temporary_html
)
stopifnot(identical(
  extract_displayed_figure_basenames(temporary_html),
  c("first-1.png", "with space-1.png")
))
unlink(temporary_html)

example_rmd <- c(
  "---",
  'title: "Example"',
  "---",
  "",
  "```{r setup, include=FALSE}",
  "x <- 1",
  "```",
  "",
  "```{r final-figure , echo=FALSE}",
  "plot(x, x)",
  "```",
  "",
  "Text that must be truncated."
)
truncated <- truncate_rmd_after_chunk(example_rmd, "final-figure")
stopifnot(
  identical(tail(truncated, 1L), "```"),
  !any(grepl("must be truncated", truncated, fixed = TRUE)),
  identical(chunk_label_from_header("```{r}"), "")
)
missing_chunk_error <- tryCatch(
  {
    truncate_rmd_after_chunk(example_rmd, "absent")
    FALSE
  },
  error = function(error) TRUE
)
stopifnot(missing_chunk_error)

expected_pngs <- c("a-1.png", "b-1.png")
expected_pdfs <- file.path(tempdir(), c("a-1.pdf", "b-1.pdf"))
stopifnot(identical(
  assert_expected_figure_set(expected_pngs, expected_pdfs),
  c("a-1.pdf", "b-1.pdf")
))

promotion_root <- tempfile("figure-promotion-test-")
staging_dir <- file.path(promotion_root, "staging")
final_dir <- file.path(promotion_root, "Figure", "page")
dir.create(staging_dir, recursive = TRUE)
dir.create(final_dir, recursive = TRUE)
writeBin(charToRaw("%PDF-new-a"), file.path(staging_dir, "a-1.pdf"))
writeBin(charToRaw("%PDF-new-b"), file.path(staging_dir, "b-1.pdf"))
writeBin(charToRaw("%PDF-old"), file.path(final_dir, "old.pdf"))
writeLines("legacy", file.path(final_dir, "legacy.png"))
promoted <- promote_staged_pdfs(
  file.path(staging_dir, c("a-1.pdf", "b-1.pdf")),
  final_dir
)
stopifnot(
  identical(sort(basename(promoted)), c("a-1.pdf", "b-1.pdf")),
  file.exists(file.path(final_dir, "legacy.png")),
  !file.exists(file.path(final_dir, "old.pdf"))
)
unlink(promotion_root, recursive = TRUE, force = TRUE)

index_lines <- readLines(file.path(repo_root, "analysis", "index.rmd"), warn = FALSE)
index_text <- paste(index_lines, collapse = "\n")
link_positions <- gregexpr("\\(([^)]+[.]html)\\)", index_text, perl = TRUE)[[1L]]
links <- regmatches(index_text, list(link_positions))[[1L]]
formal_html <- sub("^\\((.*)\\)$", "\\1", links)
formal_html <- formal_html[formal_html != "index.html"]
stopifnot(
  length(formal_html) == 10L,
  !any(grepl("_cl[.]html$", formal_html))
)

page_scripts <- sort(list.files(
  file.path(repo_root, "code", "figure_exports", "pages"),
  pattern = "^export_.*[.]R$"
))
expected_scripts <- sort(paste0(
  "export_",
  sub("[.]html$", "", formal_html),
  ".R"
))
stopifnot(identical(page_scripts, expected_scripts))

html_figure_counts <- vapply(formal_html, function(html_name) {
  length(extract_displayed_figure_basenames(file.path(repo_root, "docs", html_name)))
}, integer(1L))
stopifnot(sum(html_figure_counts) == 61L)

for (html_name in formal_html) {
  page_stem <- sub("[.]html$", "", html_name)
  expected_pngs <- extract_displayed_figure_basenames(
    file.path(repo_root, "docs", html_name)
  )
  generated_pdfs <- sort(list.files(
    file.path(repo_root, "Figure", page_stem),
    pattern = "[.]pdf$",
    full.names = TRUE
  ))
  assert_expected_figure_set(expected_pngs, generated_pdfs)
  validate_pdf_files(generated_pdfs)
}

new_r_files <- c(
  file.path(repo_root, "code", "figure_exports", "export_page_figures.R"),
  file.path(repo_root, "code", "figure_exports", "test_export_page_figures.R"),
  file.path(repo_root, "code", "figure_exports", "pages", page_scripts)
)
parse_results <- lapply(new_r_files, parse)
stopifnot(length(parse_results) == 12L)

cat("All figure export helper tests passed.\n")
