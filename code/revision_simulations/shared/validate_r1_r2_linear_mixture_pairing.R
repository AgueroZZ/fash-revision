#!/usr/bin/env Rscript

# Validate that changing only the FASH-linear prior leaves the paired R1/R2
# simulation data and every IWP/direct-comparator result exactly unchanged.
# Run this script only after the new finite-mixture caches are complete.

find_workflowr_root <- function() {
  if (file.exists("code/revision_simulations/shared/simulation_functions.R")) {
    return(".")
  }
  if (file.exists(
    "coderepo-local/code/revision_simulations/shared/simulation_functions.R"
  )) {
    return("coderepo-local")
  }
  stop("Could not find the workflowr repository root.")
}

parse_arguments <- function(args) {
  defaults <- list(
    new_r1_output_id =
      paste0(
        "r1_random_bspline_main_effect_",
        "linear_mixture_predstep1_penalty10_pilot5"
      ),
    new_r2_output_id = paste0(
      "r2_timed_cosine_one_two_three_peak_",
      "main_effect_linear_mixture_predstep1_penalty10_pilot5"
    ),
    historical_r1_output_id =
      "r1_random_bspline_main_effect_profile_sigma_pilot5",
    historical_r2_output_id = paste0(
      "r2_timed_cosine_one_two_three_peak_",
      "main_effect_profile_sigma_pilot5"
    )
  )

  option_map <- c(
    "--new-r1-output-id" = "new_r1_output_id",
    "--new-r2-output-id" = "new_r2_output_id",
    "--historical-r1-output-id" = "historical_r1_output_id",
    "--historical-r2-output-id" = "historical_r2_output_id"
  )
  if (length(args) == 1L && args[[1L]] %in% c("--help", "-h")) {
    cat(
      paste(
        "Usage: validate_r1_r2_linear_mixture_pairing.R",
        "[--new-r1-output-id ID] [--new-r2-output-id ID]",
        "[--historical-r1-output-id ID]",
        "[--historical-r2-output-id ID]"
      ),
      "\n"
    )
    quit(save = "no", status = 0L)
  }
  if (length(args) %% 2L != 0L) {
    stop("Every command-line option must be followed by a value.")
  }

  parsed <- defaults
  if (length(args) > 0L) {
    for (index in seq.int(1L, length(args), by = 2L)) {
      option <- args[[index]]
      if (!option %in% names(option_map)) {
        stop("Unknown command-line option: ", option)
      }
      value <- args[[index + 1L]]
      if (!nzchar(value)) {
        stop("Output IDs must not be empty.")
      }
      parsed[[unname(option_map[[option]])]] <- value
    }
  }
  parsed
}

read_required_csv <- function(path) {
  if (!file.exists(path)) {
    stop("Required pairing input is missing: ", path)
  }
  utils::read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

encode_key <- function(table, key_columns) {
  pieces <- lapply(table[, key_columns, drop = FALSE], function(column) {
    if (is.numeric(column)) {
      return(format(
        column,
        digits = 17L,
        scientific = TRUE,
        trim = TRUE,
        na.encode = TRUE
      ))
    }
    ifelse(is.na(column), "<NA>", as.character(column))
  })
  do.call(paste, c(pieces, sep = "\u001f"))
}

exact_cell_equal <- function(historical, candidate) {
  if (!identical(typeof(historical), typeof(candidate)) ||
      length(historical) != length(candidate)) {
    return(rep(FALSE, max(length(historical), length(candidate))))
  }
  both_missing <- is.na(historical) & is.na(candidate)
  both_present <- !is.na(historical) & !is.na(candidate)
  equal <- both_missing
  equal[both_present] <- historical[both_present] == candidate[both_present]
  equal[is.na(equal)] <- FALSE
  equal
}

empty_mismatch_table <- function() {
  data.frame(
    check = character(),
    key = character(),
    column = character(),
    historical_value = character(),
    new_value = character(),
    stringsAsFactors = FALSE
  )
}

make_failure_result <- function(check, historical_path, new_path, message) {
  list(
    summary = data.frame(
      check = check,
      historical_file = historical_path,
      new_file = new_path,
      key_columns = NA_character_,
      methods = NA_character_,
      historical_rows = NA_integer_,
      new_rows = NA_integer_,
      columns_compared = NA_integer_,
      mismatched_cells = NA_integer_,
      max_absolute_numeric_difference = NA_real_,
      status = "FAIL",
      message = message,
      stringsAsFactors = FALSE
    ),
    mismatches = data.frame(
      check = check,
      key = NA_character_,
      column = "<check error>",
      historical_value = NA_character_,
      new_value = message,
      stringsAsFactors = FALSE
    )
  )
}

compare_exact_table <- function(check,
                                historical_path,
                                new_path,
                                key_columns,
                                methods = NULL,
                                method_column = "method",
                                mismatch_limit = 100L) {
  historical <- read_required_csv(historical_path)
  candidate <- read_required_csv(new_path)

  if (!is.null(methods)) {
    if (!method_column %in% names(historical) ||
        !method_column %in% names(candidate)) {
      stop("The method column is missing from ", check, ".")
    }
    historical <- historical[historical[[method_column]] %in% methods, , drop = FALSE]
    candidate <- candidate[candidate[[method_column]] %in% methods, , drop = FALSE]
    if (!setequal(unique(historical[[method_column]]), methods) ||
        !setequal(unique(candidate[[method_column]]), methods)) {
      stop("One or more unchanged methods are missing from ", check, ".")
    }
  }
  if (!all(key_columns %in% names(historical)) ||
      !all(key_columns %in% names(candidate))) {
    stop("One or more key columns are missing from ", check, ".")
  }

  historical_key <- encode_key(historical, key_columns)
  candidate_key <- encode_key(candidate, key_columns)
  if (anyDuplicated(historical_key) || anyDuplicated(candidate_key)) {
    stop("Pairing keys are not unique in ", check, ".")
  }

  mismatches <- empty_mismatch_table()
  columns_match <- identical(names(historical), names(candidate))
  keys_match <- setequal(historical_key, candidate_key)
  mismatched_cells <- 0L
  maximum_numeric_difference <- 0

  if (!columns_match) {
    missing_from_new <- setdiff(names(historical), names(candidate))
    added_in_new <- setdiff(names(candidate), names(historical))
    mismatches <- rbind(
      mismatches,
      data.frame(
        check = check,
        key = NA_character_,
        column = "<column set>",
        historical_value = paste(missing_from_new, collapse = ";"),
        new_value = paste(added_in_new, collapse = ";"),
        stringsAsFactors = FALSE
      )
    )
    mismatched_cells <- mismatched_cells +
      length(missing_from_new) + length(added_in_new)
  }

  if (!keys_match) {
    missing_keys <- setdiff(historical_key, candidate_key)
    added_keys <- setdiff(candidate_key, historical_key)
    key_rows <- rbind(
      data.frame(
        check = check,
        key = head(missing_keys, mismatch_limit),
        column = "<missing key>",
        historical_value = "present",
        new_value = "absent",
        stringsAsFactors = FALSE
      ),
      data.frame(
        check = check,
        key = head(added_keys, mismatch_limit),
        column = "<added key>",
        historical_value = "absent",
        new_value = "present",
        stringsAsFactors = FALSE
      )
    )
    mismatches <- rbind(mismatches, head(key_rows, mismatch_limit))
    mismatched_cells <- mismatched_cells + length(missing_keys) + length(added_keys)
  }

  if (columns_match && keys_match) {
    candidate <- candidate[match(historical_key, candidate_key), , drop = FALSE]
    rownames(historical) <- NULL
    rownames(candidate) <- NULL

    for (column_name in names(historical)) {
      historical_values <- historical[[column_name]]
      candidate_values <- candidate[[column_name]]
      equal <- exact_cell_equal(historical_values, candidate_values)
      bad <- which(!equal)
      mismatched_cells <- mismatched_cells + length(bad)

      if (is.numeric(historical_values) && is.numeric(candidate_values)) {
        finite <- is.finite(historical_values) & is.finite(candidate_values)
        if (any(finite)) {
          maximum_numeric_difference <- max(
            maximum_numeric_difference,
            abs(historical_values[finite] - candidate_values[finite])
          )
        }
      }
      if (length(bad) > 0L && nrow(mismatches) < mismatch_limit) {
        keep <- head(bad, mismatch_limit - nrow(mismatches))
        mismatches <- rbind(
          mismatches,
          data.frame(
            check = check,
            key = historical_key[keep],
            column = column_name,
            historical_value = as.character(historical_values[keep]),
            new_value = as.character(candidate_values[keep]),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  passed <- columns_match && keys_match && mismatched_cells == 0L
  list(
    summary = data.frame(
      check = check,
      historical_file = historical_path,
      new_file = new_path,
      key_columns = paste(key_columns, collapse = ";"),
      methods = if (is.null(methods)) "all" else paste(methods, collapse = ";"),
      historical_rows = nrow(historical),
      new_rows = nrow(candidate),
      columns_compared = if (columns_match) ncol(historical) else NA_integer_,
      mismatched_cells = mismatched_cells,
      max_absolute_numeric_difference = if (columns_match && keys_match) {
        maximum_numeric_difference
      } else {
        NA_real_
      },
      status = if (passed) "PASS" else "FAIL",
      message = if (passed) {
        "Exact keyed equality confirmed."
      } else {
        "The historical and new paired results differ."
      },
      stringsAsFactors = FALSE
    ),
    mismatches = mismatches
  )
}

run_check <- function(..., .function = compare_exact_table) {
  arguments <- list(...)
  tryCatch(
    do.call(.function, arguments),
    error = function(error) {
      make_failure_result(
        check = arguments$check,
        historical_path = if (!is.null(arguments$historical_path)) {
          arguments$historical_path
        } else {
          arguments$historical_directory
        },
        new_path = if (!is.null(arguments$new_path)) {
          arguments$new_path
        } else {
          arguments$new_directory
        },
        message = conditionMessage(error)
      )
    }
  )
}

object_label <- function(value) {
  if (is.atomic(value) && length(value) == 1L) {
    return(as.character(value))
  }
  paste0(typeof(value), "[", length(value), "]")
}

normalize_iwp_fit_storage <- function(fits) {
  for (fit_name in names(fits)) {
    fit <- fits[[fit_name]]
    if (is.list(fit) &&
        is.list(fit$settings) &&
        !is.null(fit$settings$penalty)) {
      fit$settings$penalty <- as.numeric(fit$settings$penalty)
    }
    fits[[fit_name]] <- fit
  }
  fits
}

compare_replicate_metadata <- function(check,
                                       historical_directory,
                                       new_directory,
                                       scenario) {
  historical_files <- sort(list.files(
    historical_directory,
    pattern = "^seed_[0-9]+[.]rds$",
    full.names = FALSE
  ))
  new_files <- sort(list.files(
    new_directory,
    pattern = "^seed_[0-9]+[.]rds$",
    full.names = FALSE
  ))
  if (!identical(historical_files, new_files) || length(historical_files) == 0L) {
    stop("The paired replicate RDS file sets do not match for ", scenario, ".")
  }

  required_fields <- if (identical(scenario, "R1")) {
    c("seed", "permutation_seed", "true_pi0")
  } else {
    c("seed", "component_seeds")
  }
  mismatches <- empty_mismatch_table()
  compared <- 0L
  for (file_name in historical_files) {
    historical <- readRDS(file.path(historical_directory, file_name))
    candidate <- readRDS(file.path(new_directory, file_name))
    if (!all(required_fields %in% names(historical)) ||
        !all(required_fields %in% names(candidate)) ||
        is.null(historical$configuration) ||
        is.null(candidate$configuration)) {
      stop("Required replicate metadata is missing from ", file_name, ".")
    }

    for (field in required_fields) {
      compared <- compared + 1L
      if (!identical(historical[[field]], candidate[[field]])) {
        mismatches <- rbind(
          mismatches,
          data.frame(
            check = check,
            key = file_name,
            column = field,
            historical_value = object_label(historical[[field]]),
            new_value = object_label(candidate[[field]]),
            stringsAsFactors = FALSE
          )
        )
      }
    }

    allowed_configuration_differences <- c(
      "output_id",
      "linear_sigma_estimation",
      "linear_sigma_grid",
      "linear_prior_mode",
      "common_sd_grid",
      "common_pred_step",
      "common_penalty"
    )
    historical_configuration_fields <- setdiff(
      names(historical$configuration),
      allowed_configuration_differences
    )
    candidate_configuration_fields <- setdiff(
      names(candidate$configuration),
      allowed_configuration_differences
    )
    if (length(historical_configuration_fields) == 0L ||
        !setequal(
          historical_configuration_fields,
          candidate_configuration_fields
        )) {
      stop(
        "The non-prior scientific configuration fields differ in ",
        file_name, "."
      )
    }
    for (field in historical_configuration_fields) {
      compared <- compared + 1L
      if (!identical(
            historical$configuration[[field]],
            candidate$configuration[[field]]
          )) {
        mismatches <- rbind(
          mismatches,
          data.frame(
            check = check,
            key = file_name,
            column = paste0("configuration.", field),
            historical_value = object_label(historical$configuration[[field]]),
            new_value = object_label(candidate$configuration[[field]]),
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }

  passed <- nrow(mismatches) == 0L
  list(
    summary = data.frame(
      check = check,
      historical_file = historical_directory,
      new_file = new_directory,
      key_columns = "replicate_file;metadata_field",
      methods = "not_applicable",
      historical_rows = length(historical_files),
      new_rows = length(new_files),
      columns_compared = compared,
      mismatched_cells = nrow(mismatches),
      max_absolute_numeric_difference = if (passed) 0 else NA_real_,
      status = if (passed) "PASS" else "FAIL",
      message = if (passed) {
        "Seeds and shared scientific configuration are exactly paired."
      } else {
        "Replicate metadata differs between historical and new caches."
      },
      stringsAsFactors = FALSE
    ),
    mismatches = mismatches
  )
}

compare_full_fit_inputs <- function(check,
                                    historical_path,
                                    new_path,
                                    scenario) {
  if (!file.exists(historical_path) || !file.exists(new_path)) {
    stop("A paired full-fit RDS file is missing for ", scenario, ".")
  }
  historical <- readRDS(historical_path)
  candidate <- readRDS(new_path)
  fields <- c(
    "datasets", "unit_info", "genotype", "variant_info", "covariates",
    "true_beta", "expression", "expression_simulation", "eqtl_summary",
    "se_correction_summary", "fash_fits", "direct_interaction_lrt",
    "direct_interaction_efdr"
  )
  if (identical(scenario, "R2")) {
    fields <- c(
      fields,
      "true_beta_evaluation", "evaluation_grid", "true_functionals"
    )
  }
  if (!all(fields %in% names(historical)) ||
      !all(fields %in% names(candidate))) {
    stop("Required full-fit pairing fields are missing for ", scenario, ".")
  }

  mismatches <- empty_mismatch_table()
  full_fit_key <- tools::file_path_sans_ext(basename(new_path))
  for (field in fields) {
    historical_value <- historical[[field]]
    candidate_value <- candidate[[field]]
    if (identical(field, "fash_fits")) {
      historical_value <- normalize_iwp_fit_storage(historical_value)
      candidate_value <- normalize_iwp_fit_storage(candidate_value)
    }
    if (!identical(historical_value, candidate_value)) {
      mismatches <- rbind(
        mismatches,
        data.frame(
          check = check,
          key = full_fit_key,
          column = field,
          historical_value = object_label(historical_value),
          new_value = object_label(candidate_value),
          stringsAsFactors = FALSE
        )
      )
    }
  }
  passed <- nrow(mismatches) == 0L
  list(
    summary = data.frame(
      check = check,
      historical_file = historical_path,
      new_file = new_path,
      key_columns = "seed;object_field",
      methods = "IWP_and_direct_inputs",
      historical_rows = 1L,
      new_rows = 1L,
      columns_compared = length(fields),
      mismatched_cells = nrow(mismatches),
      max_absolute_numeric_difference = if (passed) 0 else NA_real_,
      status = if (passed) "PASS" else "FAIL",
      message = if (passed) {
        "Seed-12345 inputs, truth, IWP fits, and direct-test objects are exact."
      } else {
        "The paired seed-12345 full-fit objects differ."
      },
      stringsAsFactors = FALSE
    ),
    mismatches = mismatches
  )
}

write_csv_atomic <- function(table, path) {
  temporary_path <- tempfile(
    pattern = paste0(basename(path), "."),
    tmpdir = dirname(path)
  )
  on.exit(unlink(temporary_path), add = TRUE)
  utils::write.csv(table, temporary_path, row.names = FALSE)
  if (!file.rename(temporary_path, path)) {
    stop("Could not finalize validation output: ", path)
  }
  invisible(path)
}

validate_scenario <- function(scenario,
                              historical_summary_dir,
                              new_summary_dir,
                              include_r2_checks = FALSE) {
  if (!dir.exists(historical_summary_dir)) {
    stop("Historical summary directory is missing: ", historical_summary_dir)
  }
  if (!dir.exists(new_summary_dir)) {
    stop(
      "New summary directory is missing. Complete the new cache before ",
      "running this validator: ", new_summary_dir
    )
  }

  unchanged_methods <- c(
    "FASH-IWP1-Raw",
    "FASH-IWP1-BF",
    "Direct-linear-LRT-eFDR-true-pi0",
    "Direct-quadratic-LRT-eFDR-true-pi0"
  )
  paired_file <- function(directory, file_name) file.path(directory, file_name)

  checks <- list(
    run_check(
      check = paste0(scenario, "_replicate_metadata"),
      historical_directory = file.path(
        dirname(historical_summary_dir),
        "replicates"
      ),
      new_directory = file.path(dirname(new_summary_dir), "replicates"),
      scenario = scenario,
      .function = compare_replicate_metadata
    ),
    run_check(
      check = paste0(scenario, "_seed_12345_full_fit_inputs"),
      historical_path = file.path(
        dirname(historical_summary_dir),
        "full_fits",
        "seed_12345.rds"
      ),
      new_path = file.path(
        dirname(new_summary_dir),
        "full_fits",
        "seed_12345.rds"
      ),
      scenario = scenario,
      .function = compare_full_fit_inputs
    ),
    run_check(
      check = paste0(scenario, "_all_replicate_alpha_curves"),
      historical_path = paired_file(
        historical_summary_dir,
        "all_replicate_alpha_curves.csv"
      ),
      new_path = paired_file(new_summary_dir, "all_replicate_alpha_curves.csv"),
      key_columns = c("seed", "method", "alpha"),
      methods = unchanged_methods
    ),
    run_check(
      check = paste0(scenario, "_all_replicate_alpha005"),
      historical_path = paired_file(
        historical_summary_dir,
        "all_replicate_alpha005.csv"
      ),
      new_path = paired_file(new_summary_dir, "all_replicate_alpha005.csv"),
      key_columns = c("seed", "method", "alpha"),
      methods = unchanged_methods
    )
  )

  if (include_r2_checks) {
    checks <- c(
      checks,
      list(
        run_check(
          check = "R2_all_replicate_geometry",
          historical_path = paired_file(
            historical_summary_dir,
            "all_replicate_geometry.csv"
          ),
          new_path = paired_file(
            new_summary_dir,
            "all_replicate_geometry.csv"
          ),
          key_columns = "seed"
        ),
        run_check(
          check = "R2_all_replicate_peak_alpha_curves",
          historical_path = paired_file(
            historical_summary_dir,
            "all_replicate_peak_alpha_curves.csv"
          ),
          new_path = paired_file(
            new_summary_dir,
            "all_replicate_peak_alpha_curves.csv"
          ),
          key_columns = c(
            "seed", "method", "alpha", "subgroup_var", "subgroup_value"
          ),
          methods = unchanged_methods
        ),
        run_check(
          check = "R2_all_replicate_peak_alpha005",
          historical_path = paired_file(
            historical_summary_dir,
            "all_replicate_peak_alpha005.csv"
          ),
          new_path = paired_file(
            new_summary_dir,
            "all_replicate_peak_alpha005.csv"
          ),
          key_columns = c(
            "seed", "method", "alpha", "subgroup_var", "subgroup_value"
          ),
          methods = unchanged_methods
        )
      )
    )
  }

  validation_summary <- do.call(rbind, lapply(checks, `[[`, "summary"))
  mismatch_parts <- lapply(checks, `[[`, "mismatches")
  mismatch_parts <- mismatch_parts[vapply(mismatch_parts, nrow, integer(1)) > 0L]
  validation_mismatches <- if (length(mismatch_parts) == 0L) {
    empty_mismatch_table()
  } else {
    do.call(rbind, mismatch_parts)
  }

  write_csv_atomic(
    validation_summary,
    file.path(new_summary_dir, "unchanged_method_pairing_validation.csv")
  )
  mismatch_path <- file.path(
    new_summary_dir,
    "unchanged_method_pairing_mismatches.csv"
  )
  if (nrow(validation_mismatches) > 0L) {
    write_csv_atomic(validation_mismatches, mismatch_path)
  } else if (file.exists(mismatch_path)) {
    unlink(mismatch_path)
  }

  validation_summary
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
workflowr_root <- find_workflowr_root()
mc_root <- file.path(
  workflowr_root,
  "output",
  "revision_simulations",
  "mc"
)
summary_directory <- function(output_id) {
  file.path(mc_root, output_id, "summary")
}

r1_validation <- validate_scenario(
  scenario = "R1",
  historical_summary_dir = summary_directory(
    arguments$historical_r1_output_id
  ),
  new_summary_dir = summary_directory(arguments$new_r1_output_id),
  include_r2_checks = FALSE
)
r2_validation <- validate_scenario(
  scenario = "R2",
  historical_summary_dir = summary_directory(
    arguments$historical_r2_output_id
  ),
  new_summary_dir = summary_directory(arguments$new_r2_output_id),
  include_r2_checks = TRUE
)

all_validation <- rbind(r1_validation, r2_validation)
print(all_validation[, c(
  "check", "historical_rows", "new_rows", "mismatched_cells", "status"
)])
if (any(all_validation$status != "PASS")) {
  stop(
    "R1/R2 linear-mixture pairing validation failed. See the validation ",
    "CSV files in the new cache summary directories."
  )
}
cat("R1/R2 linear-mixture pairing validation passed exactly.\n")
