# Validate the shared-function dependency closure used by formal R7 reporting.

matched_r7_shared_root_functions <- function() {
  c(
    "default_revision_grid",
    "compute_linear_mixture_log_likelihood",
    "fit_linear_mixture_fash_from_log_likelihood",
    "validate_linear_mixture_fash",
    "BF_update_linear_mixture_fash",
    "expand_grid_prior_weights"
  )
}

matched_r7_function_names <- function(environment) {
  objects <- mget(
    ls(environment, all.names = TRUE),
    envir = environment,
    inherits = FALSE
  )
  names(objects)[vapply(objects, is.function, logical(1))]
}

matched_r7_dependency_closure <- function(environment, root_functions) {
  if (!is.environment(environment)) {
    stop("environment must be an environment.")
  }
  root_functions <- unique(as.character(root_functions))
  if (length(root_functions) == 0L || any(!nzchar(root_functions))) {
    stop("root_functions must contain at least one function name.")
  }
  available_functions <- matched_r7_function_names(environment)
  missing_roots <- setdiff(root_functions, available_functions)
  if (length(missing_roots) > 0L) {
    stop(
      "Missing shared root function(s): ",
      paste(missing_roots, collapse = ", "),
      "."
    )
  }

  closure <- character()
  queue <- root_functions
  while (length(queue) > 0L) {
    function_name <- queue[[1L]]
    queue <- queue[-1L]
    if (function_name %in% closure) {
      next
    }
    closure <- c(closure, function_name)
    globals <- codetools::findGlobals(
      environment[[function_name]],
      merge = FALSE
    )$functions
    dependencies <- intersect(globals, available_functions)
    queue <- c(queue, setdiff(dependencies, closure))
  }
  sort(unique(closure))
}

validate_matched_r7_scoped_shared_provenance <- function(
    historical_path,
    current_path,
    recorded_md5,
    root_functions = matched_r7_shared_root_functions()) {
  if (!requireNamespace("codetools", quietly = TRUE)) {
    stop("The codetools package is required for scoped provenance validation.")
  }
  historical_path <- normalizePath(
    historical_path,
    winslash = "/",
    mustWork = TRUE
  )
  current_path <- normalizePath(
    current_path,
    winslash = "/",
    mustWork = TRUE
  )
  if (!is.character(recorded_md5) || length(recorded_md5) != 1L ||
      is.na(recorded_md5) || !grepl("^[0-9a-f]{32}$", recorded_md5)) {
    stop("recorded_md5 must be one lowercase MD5 string.")
  }

  historical_file_md5 <- unname(tools::md5sum(historical_path))
  current_file_md5 <- unname(tools::md5sum(current_path))
  if (!identical(historical_file_md5, recorded_md5)) {
    stop(
      "The frozen shared-function snapshot does not match the cache-recorded MD5."
    )
  }

  historical_environment <- new.env(parent = baseenv())
  current_environment <- new.env(parent = baseenv())
  sys.source(historical_path, envir = historical_environment)
  sys.source(current_path, envir = current_environment)

  historical_closure <- matched_r7_dependency_closure(
    historical_environment,
    root_functions
  )
  current_closure <- matched_r7_dependency_closure(
    current_environment,
    root_functions
  )
  if (!identical(historical_closure, current_closure)) {
    stop("The current R7 shared-function dependency closure changed.")
  }

  comparison <- data.frame(
    function_name = historical_closure,
    formals_identical = vapply(historical_closure, function(function_name) {
      identical(
        formals(historical_environment[[function_name]]),
        formals(current_environment[[function_name]])
      )
    }, logical(1)),
    body_identical = vapply(historical_closure, function(function_name) {
      identical(
        body(historical_environment[[function_name]]),
        body(current_environment[[function_name]])
      )
    }, logical(1)),
    stringsAsFactors = FALSE
  )
  if (any(!comparison$formals_identical) || any(!comparison$body_identical)) {
    changed <- comparison$function_name[
      !comparison$formals_identical | !comparison$body_identical
    ]
    stop(
      "The current R7 shared-function contract changed for: ",
      paste(changed, collapse = ", "),
      "."
    )
  }

  list(
    passed = TRUE,
    historical_path = historical_path,
    current_path = current_path,
    historical_file_md5 = historical_file_md5,
    current_file_md5 = current_file_md5,
    root_functions = root_functions,
    dependency_functions = historical_closure,
    comparison = comparison
  )
}
