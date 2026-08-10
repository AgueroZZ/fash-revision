#!/usr/bin/env Rscript

# Download and validate baselineLD v2.2 chromosome annotation files.

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

validate_gzip_magic <- function(path) {
  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("Annotation file is missing or empty: ", path)
  }
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  magic <- readBin(connection, what = "raw", n = 2L)
  if (!identical(as.integer(magic), c(31L, 139L))) {
    stop("Annotation file is not gzip encoded: ", path)
  }
  invisible(TRUE)
}

sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("The digest package is required for SHA256 validation.")
  }
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

validate_source_file <- function(path, expected_size, expected_sha256) {
  validate_gzip_magic(path)
  observed_size <- unname(file.info(path)$size)
  if (!identical(as.numeric(observed_size), as.numeric(expected_size))) {
    stop(
      "Unexpected byte size for ", basename(path), ": expected ",
      expected_size, ", observed ", observed_size, "."
    )
  }
  observed_sha256 <- sha256_file(path)
  if (!identical(observed_sha256, expected_sha256)) {
    stop(
      "SHA256 mismatch for ", basename(path), ": expected ",
      expected_sha256, ", observed ", observed_sha256, "."
    )
  }
  invisible(TRUE)
}

download_one_source <- function(source, target_directory) {
  target_path <- file.path(target_directory, source$local_filename)
  existing_validation <- try(
    validate_source_file(
      target_path,
      expected_size = source$byte_size,
      expected_sha256 = source$sha256
    ),
    silent = TRUE
  )
  if (inherits(existing_validation, "try-error")) {
    temporary_path <- tempfile(
      pattern = paste0("baselineLD-chr", source$chromosome, "-"),
      tmpdir = target_directory,
      fileext = ".download"
    )
    on.exit(unlink(temporary_path), add = TRUE)
    status <- utils::download.file(
      url = source$retrieval_url,
      destfile = temporary_path,
      method = "libcurl",
      mode = "wb",
      quiet = FALSE
    )
    if (!identical(status, 0L)) {
      stop("Download failed for chromosome ", source$chromosome, ".")
    }
    validate_source_file(
      temporary_path,
      expected_size = source$byte_size,
      expected_sha256 = source$sha256
    )
    if (file.exists(target_path)) {
      unlink(target_path)
    }
    if (!file.rename(temporary_path, target_path)) {
      stop("Could not atomically install annotation: ", target_path)
    }
  }

  validate_source_file(
    target_path,
    expected_size = source$byte_size,
    expected_sha256 = source$sha256
  )
  data.frame(
    chromosome = source$chromosome,
    canonical_record = source$canonical_record,
    canonical_archive = source$canonical_archive,
    retrieval_url = source$retrieval_url,
    local_path = normalizePath(target_path, winslash = "/", mustWork = TRUE),
    byte_size = unname(file.info(target_path)$size),
    sha256 = sha256_file(target_path),
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
  "baseline_ld_variant_enrichment"
)
source_manifest_path <- file.path(script_directory, "baseline_ld_sources.csv")
target_directory <- file.path(
  workflowr_root,
  "data",
  "revision_simulations",
  "internal",
  "baseline_ld_variant_enrichment"
)
dir.create(target_directory, recursive = TRUE, showWarnings = FALSE)

manifest <- utils::read.csv(
  source_manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
required_columns <- c(
  "chromosome", "canonical_record", "canonical_archive", "retrieval_url",
  "local_filename", "byte_size", "sha256"
)
if (!identical(names(manifest), required_columns) ||
    !identical(sort(as.integer(manifest$chromosome)), 1:22) ||
    anyDuplicated(manifest$chromosome) ||
    anyDuplicated(manifest$local_filename) ||
    any(manifest$byte_size <= 0) ||
    any(!grepl("^[0-9a-f]{64}$", manifest$sha256)) ||
    sum(manifest$byte_size) != 416186492) {
  stop("The baselineLD v2.2 source manifest is invalid.")
}

download_rows <- lapply(seq_len(nrow(manifest)), function(index) {
  message(
    "Validating baselineLD v2.2 chromosome ",
    manifest$chromosome[index],
    "/22."
  )
  download_one_source(manifest[index, , drop = FALSE], target_directory)
})
download_manifest <- do.call(rbind, download_rows)
utils::write.csv(
  download_manifest,
  file.path(target_directory, "download_manifest.csv"),
  row.names = FALSE,
  quote = TRUE
)
message(
  "Validated 22 baselineLD annotation files (",
  sum(download_manifest$byte_size),
  " bytes)."
)
