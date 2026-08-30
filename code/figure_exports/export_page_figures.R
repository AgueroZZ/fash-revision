trim_scalar <- function(x) {
  trimws(as.character(x)[[1L]])
}

find_export_repo_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (
      file.exists(file.path(current, "_workflowr.yml")) &&
        file.exists(file.path(current, "analysis", "index.rmd"))
    ) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the workflowr repository root from: ", start)
    }
    current <- parent
  }
}

extract_displayed_figure_basenames <- function(html_path) {
  if (!file.exists(html_path)) {
    stop("Rendered page does not exist: ", html_path)
  }
  html <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
  pattern <- "figure/[^\"'<>[:space:]]+[.](?:png|PNG)"
  positions <- gregexpr(pattern, html, perl = TRUE)[[1L]]
  if (identical(positions[[1L]], -1L)) {
    stop("No workflowr PNG references were found in: ", html_path)
  }
  references <- regmatches(html, list(positions))[[1L]]
  decoded <- utils::URLdecode(references)
  sort(unique(basename(decoded)))
}

chunk_label_from_header <- function(header) {
  body <- sub("^```\\{r", "", header)
  body <- sub("\\}[[:space:]]*$", "", body)
  fields <- strsplit(body, ",", fixed = TRUE)[[1L]]
  first_field <- if (length(fields) == 0L) "" else fields[[1L]]
  trim_scalar(first_field)
}

truncate_rmd_after_chunk <- function(lines, final_chunk) {
  chunk_starts <- grep("^```\\{r(?:[[:space:],}]|$)", lines, perl = TRUE)
  labels <- vapply(lines[chunk_starts], chunk_label_from_header, character(1L))
  matches <- chunk_starts[labels == final_chunk]
  if (length(matches) != 1L) {
    stop(
      "Expected exactly one final chunk named '", final_chunk,
      "'; found ", length(matches), "."
    )
  }
  start <- matches[[1L]]
  later_lines <- seq.int(start + 1L, length(lines))
  closing_candidates <- later_lines[grepl(
    "^```[[:space:]]*$",
    lines[later_lines],
    perl = TRUE
  )]
  if (length(closing_candidates) == 0L) {
    stop("Final chunk has no closing fence: ", final_chunk)
  }
  lines[seq_len(closing_candidates[[1L]])]
}

install_offline_biomart_adapter <- function(page_env, cache_path) {
  if (!file.exists(cache_path)) {
    stop("Offline BioMart cache does not exist: ", cache_path)
  }
  gene_map <- as.data.frame(readRDS(cache_path), stringsAsFactors = FALSE)
  required <- c("ensembl_gene_id", "hgnc_symbol")
  if (!all(required %in% names(gene_map))) {
    stop(
      "Offline BioMart cache lacks required columns: ",
      paste(setdiff(required, names(gene_map)), collapse = ", ")
    )
  }
  gene_map <- unique(gene_map[required])

  offline_mart <- structure(list(cache_path = cache_path), class = "offline_mart")
  page_env$useEnsembl <- function(...) offline_mart
  page_env$useMart <- function(...) offline_mart
  page_env$getBM <- function(attributes, filters, values, mart, ...) {
    attributes <- as.character(attributes)
    filters <- trim_scalar(filters)
    if (!filters %in% names(gene_map)) {
      stop("Offline BioMart adapter cannot filter by: ", filters)
    }
    if (!all(attributes %in% names(gene_map))) {
      stop(
        "Offline BioMart adapter cannot return: ",
        paste(setdiff(attributes, names(gene_map)), collapse = ", ")
      )
    }
    keep <- gene_map[[filters]] %in% as.character(values)
    unique(gene_map[keep, attributes, drop = FALSE])
  }
  invisible(gene_map)
}

install_plot_export_overrides <- function(
    page_env,
    rasterize_point_chunks = character(),
    raster_dpi = 300L) {
  rasterize_point_chunks <- unique(as.character(rasterize_point_chunks))
  if (length(rasterize_point_chunks) > 0L) {
    if (!requireNamespace("ggrastr", quietly = TRUE)) {
      stop("Package 'ggrastr' is required for selective point rasterization.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' is required for selective point rasterization.")
    }
    ggplot_method_class <- if (
      any(methods("print") == "print.ggplot2::ggplot")
    ) {
      "ggplot2::ggplot"
    } else {
      "ggplot"
    }
    original_print_ggplot <- getS3method(
      "print",
      ggplot_method_class,
      envir = asNamespace("ggplot2")
    )
    export_print_ggplot <- local({
      selected_chunks <- rasterize_point_chunks
      selected_dpi <- as.integer(raster_dpi)
      original_print <- original_print_ggplot
      function(x, ...) {
        label <- knitr::opts_current$get("label")
        if (!is.null(label) && label %in% selected_chunks) {
          x <- ggrastr::rasterise(x, layers = "Point", dpi = selected_dpi)
        }
        original_print(x, ...)
      }
    })
    registerS3method(
      "print",
      ggplot_method_class,
      export_print_ggplot,
      envir = asNamespace("base")
    )
  }

  page_env$ggsave <- function(filename, ...) {
    invisible(filename)
  }
  invisible(page_env)
}

validate_pdf_files <- function(paths) {
  if (length(paths) == 0L) {
    stop("No PDF files were generated.")
  }
  sizes <- file.info(paths)$size
  bad_size <- is.na(sizes) | sizes < 1000
  if (any(bad_size)) {
    stop(
      "Generated PDF files are empty or implausibly small: ",
      paste(basename(paths[bad_size]), collapse = ", ")
    )
  }
  bad_header <- vapply(paths, function(path) {
    connection <- file(path, open = "rb")
    header <- readBin(connection, what = "raw", n = 5L)
    close(connection)
    !identical(rawToChar(header), "%PDF-")
  }, logical(1L))
  if (any(bad_header)) {
    stop(
      "Generated files do not have a PDF header: ",
      paste(basename(paths[bad_header]), collapse = ", ")
    )
  }
  invisible(TRUE)
}

assert_expected_figure_set <- function(expected_pngs, generated_pdfs) {
  expected_pdfs <- sub("[.](?:png|PNG)$", ".pdf", expected_pngs, perl = TRUE)
  actual_pdfs <- basename(generated_pdfs)
  missing <- setdiff(expected_pdfs, actual_pdfs)
  extra <- setdiff(actual_pdfs, expected_pdfs)
  if (length(missing) > 0L || length(extra) > 0L) {
    details <- c(
      if (length(missing) > 0L) {
        paste0("missing: ", paste(missing, collapse = ", "))
      },
      if (length(extra) > 0L) {
        paste0("extra: ", paste(extra, collapse = ", "))
      }
    )
    stop("Generated PDF set does not match rendered page (", paste(details, collapse = "; "), ").")
  }
  invisible(sort(expected_pdfs))
}

promote_staged_pdfs <- function(staged_pdfs, final_dir) {
  parent_dir <- dirname(final_dir)
  dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
  promotion_dir <- tempfile(
    pattern = paste0(".promote-", basename(final_dir), "-"),
    tmpdir = parent_dir
  )
  backup_dir <- tempfile(
    pattern = paste0(".backup-", basename(final_dir), "-"),
    tmpdir = parent_dir
  )
  dir.create(promotion_dir, recursive = FALSE, showWarnings = FALSE)
  on.exit(unlink(promotion_dir, recursive = TRUE, force = TRUE), add = TRUE)
  on.exit(unlink(backup_dir, recursive = TRUE, force = TRUE), add = TRUE)

  if (dir.exists(final_dir)) {
    existing <- list.files(
      final_dir,
      all.files = TRUE,
      no.. = TRUE,
      full.names = TRUE
    )
    retained_non_pdf <- existing[!grepl("[.]pdf$", basename(existing))]
    if (length(retained_non_pdf) > 0L) {
      retained_copy <- file.copy(
        retained_non_pdf,
        promotion_dir,
        recursive = TRUE,
        copy.mode = TRUE,
        copy.date = TRUE
      )
      if (!all(retained_copy)) {
        stop("Failed to stage retained non-PDF page artifacts.")
      }
    }
  }

  copied <- file.copy(
    staged_pdfs,
    file.path(promotion_dir, basename(staged_pdfs)),
    overwrite = FALSE,
    copy.mode = TRUE,
    copy.date = FALSE
  )
  if (!all(copied)) {
    stop(
      "Failed to promote staged PDFs: ",
      paste(basename(staged_pdfs[!copied]), collapse = ", ")
    )
  }

  had_existing_dir <- dir.exists(final_dir)
  if (had_existing_dir && !file.rename(final_dir, backup_dir)) {
    stop("Failed to stage the previous complete page directory for replacement.")
  }
  promoted <- file.rename(promotion_dir, final_dir)
  if (!promoted) {
    if (had_existing_dir) {
      file.rename(backup_dir, final_dir)
    }
    stop("Failed to atomically promote the complete page directory.")
  }
  if (had_existing_dir) {
    unlink(backup_dir, recursive = TRUE, force = TRUE)
  }
  invisible(file.path(final_dir, basename(staged_pdfs)))
}

export_page_figures <- function(
    page_rmd,
    page_html,
    final_chunk,
    repo_root = find_export_repo_root(),
    rasterize_point_chunks = character(),
    raster_dpi = 300L,
    offline_biomart_cache = NULL) {
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  rmd_path <- normalizePath(
    file.path(repo_root, page_rmd),
    winslash = "/",
    mustWork = TRUE
  )
  html_path <- normalizePath(
    file.path(repo_root, page_html),
    winslash = "/",
    mustWork = TRUE
  )
  page_stem <- sub("[.][^.]+$", "", basename(rmd_path))
  figure_root <- file.path(repo_root, "Figure")
  dir.create(figure_root, recursive = TRUE, showWarnings = FALSE)
  staging_dir <- tempfile(
    pattern = paste0(".staging-", page_stem, "-"),
    tmpdir = figure_root
  )
  dir.create(staging_dir, recursive = FALSE, showWarnings = FALSE)
  on.exit(unlink(staging_dir, recursive = TRUE, force = TRUE), add = TRUE)

  expected_pngs <- extract_displayed_figure_basenames(html_path)
  rmd_lines <- readLines(rmd_path, warn = FALSE)
  export_lines <- truncate_rmd_after_chunk(rmd_lines, final_chunk)
  staged_rmd <- file.path(staging_dir, basename(rmd_path))
  staged_md <- file.path(staging_dir, paste0(page_stem, ".md"))
  writeLines(export_lines, staged_rmd, useBytes = TRUE)

  page_env <- new.env(parent = globalenv())
  install_plot_export_overrides(
    page_env,
    rasterize_point_chunks = rasterize_point_chunks,
    raster_dpi = raster_dpi
  )
  if (!is.null(offline_biomart_cache)) {
    install_offline_biomart_adapter(
      page_env,
      file.path(repo_root, offline_biomart_cache)
    )
  }

  old_working_directory <- getwd()
  old_root_directory <- knitr::opts_knit$get("root.dir")
  on.exit(setwd(old_working_directory), add = TRUE)
  on.exit(knitr::opts_knit$set(root.dir = old_root_directory), add = TRUE)
  setwd(repo_root)
  knitr::opts_knit$set(root.dir = repo_root)
  knitr::opts_chunk$set(
    dev = "cairo_pdf",
    fig.ext = "pdf",
    fig.path = paste0(staging_dir, "/"),
    error = FALSE,
    export_force_warning = TRUE
  )
  knitr::opts_hooks$set(export_force_warning = function(options) {
    options$warning <- TRUE
    options
  })

  message("Exporting editable PDFs for ", page_stem, " ...")
  knitr::knit(
    input = staged_rmd,
    output = staged_md,
    quiet = FALSE,
    envir = page_env,
    encoding = "UTF-8"
  )
  markdown <- paste(readLines(staged_md, warn = FALSE), collapse = "\n")
  if (grepl("Removed [0-9]+ rows containing missing values", markdown, ignore.case = TRUE)) {
    stop("A plot dropped rows during PDF export; inspect the staged knit warnings.")
  }

  staged_pdfs <- sort(list.files(
    staging_dir,
    pattern = "[.]pdf$",
    full.names = TRUE,
    recursive = FALSE
  ))
  assert_expected_figure_set(expected_pngs, staged_pdfs)
  validate_pdf_files(staged_pdfs)

  final_dir <- file.path(figure_root, page_stem)
  final_paths <- promote_staged_pdfs(staged_pdfs, final_dir)
  message(
    "Exported ", length(final_paths), " PDF figure(s) to ",
    normalizePath(final_dir, winslash = "/", mustWork = TRUE)
  )
  invisible(data.frame(
    page = page_stem,
    pdf = final_paths,
    rasterized_points = vapply(
      sub("[.]pdf$", "", basename(final_paths)),
      function(figure_name) {
        any(startsWith(figure_name, paste0(rasterize_point_chunks, "-")))
      },
      logical(1L)
    ),
    stringsAsFactors = FALSE
  ))
}
