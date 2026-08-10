#!/usr/bin/env Rscript

# Download and validate versioned hg19-compatible functional annotations.

find_workflowr_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "_workflowr.yml"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the workflowr project root.")
    }
    current <- parent
  }
}

validate_gzip_file <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Downloaded annotation is missing or empty: ", path)
  }
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  magic <- readBin(connection, what = "raw", n = 2L)
  if (!identical(as.integer(magic), c(31L, 139L))) {
    stop("Downloaded annotation is not a gzip file: ", path)
  }
  invisible(TRUE)
}

checksum_file <- function(path, algorithm) {
  if (algorithm == "md5") {
    return(unname(tools::md5sum(path)))
  }
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required for SHA256 checksums.")
  }
  digest::digest(file = path, algo = algorithm, serialize = FALSE)
}

download_one_source <- function(source, target_directory) {
  target_path <- file.path(target_directory, source$local_filename)
  expected_md5 <- source$expected_md5
  if (is.na(expected_md5)) {
    expected_md5 <- ""
  }

  existing_is_valid <- FALSE
  if (file.exists(target_path)) {
    validation <- try(validate_gzip_file(target_path), silent = TRUE)
    existing_is_valid <- !inherits(validation, "try-error")
    if (existing_is_valid && nzchar(expected_md5)) {
      existing_is_valid <- identical(
        checksum_file(target_path, "md5"),
        expected_md5
      )
    }
  }

  if (!existing_is_valid) {
    temporary_path <- tempfile(
      pattern = paste0(source$source_id, "-"),
      tmpdir = target_directory,
      fileext = ".download"
    )
    on.exit(unlink(temporary_path), add = TRUE)
    status <- utils::download.file(
      url = source$url,
      destfile = temporary_path,
      method = "libcurl",
      mode = "wb",
      quiet = FALSE
    )
    if (!identical(status, 0L)) {
      stop("Download failed for source ", source$source_id, ".")
    }
    validate_gzip_file(temporary_path)
    observed_md5 <- checksum_file(temporary_path, "md5")
    if (nzchar(expected_md5) && !identical(observed_md5, expected_md5)) {
      stop(
        "MD5 mismatch for ", source$source_id, ": expected ",
        expected_md5, ", observed ", observed_md5, "."
      )
    }
    if (file.exists(target_path)) {
      unlink(target_path)
    }
    if (!file.rename(temporary_path, target_path)) {
      stop("Could not atomically install annotation: ", target_path)
    }
  }

  validate_gzip_file(target_path)
  data.frame(
    source_id = source$source_id,
    source_name = source$source_name,
    assembly = source$assembly,
    accession = source$expected_accession,
    url = source$url,
    local_path = target_path,
    byte_size = unname(file.info(target_path)$size),
    md5 = checksum_file(target_path, "md5"),
    sha256 = checksum_file(target_path, "sha256"),
    retrieved_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}

workflowr_root <- find_workflowr_root()
script_directory <- file.path(
  workflowr_root,
  "code",
  "revision_simulations",
  "internal",
  "variant_annotation_enrichment"
)
manifest_path <- file.path(script_directory, "annotation_sources.csv")
target_directory <- file.path(
  workflowr_root,
  "data",
  "revision_simulations",
  "internal",
  "variant_annotation_enrichment"
)
dir.create(target_directory, recursive = TRUE, showWarnings = FALSE)

manifest <- utils::read.csv(
  manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("", "NA")
)
required_columns <- c(
  "source_id", "source_name", "assembly", "url", "local_filename",
  "expected_accession", "expected_md5"
)
if (!all(required_columns %in% names(manifest)) || !nrow(manifest) ||
    anyDuplicated(manifest$source_id) ||
    any(!manifest$assembly %in% c("hg19", "GRCh37")) ||
    any(!nzchar(manifest$url)) || any(!nzchar(manifest$local_filename))) {
  stop("The annotation source manifest is invalid.")
}

download_rows <- lapply(seq_len(nrow(manifest)), function(index) {
  message("Validating annotation source: ", manifest$source_id[index])
  download_one_source(manifest[index, , drop = FALSE], target_directory)
})
download_manifest <- do.call(rbind, download_rows)
utils::write.csv(
  download_manifest,
  file.path(target_directory, "download_manifest.csv"),
  row.names = FALSE,
  quote = TRUE
)
message("Validated ", nrow(download_manifest), " annotation sources.")
