#!/usr/bin/env Rscript

# Benchmark the Gaussian Omega likelihood path against the diagonal-S path.
#
# The primary benchmark measures the implementation that fashr currently uses:
# fash_L_compute() rebuilds the TMB data and objective for every unit-by-grid
# component. Secondary benchmarks separate zero and positive PSD components,
# describe the supplied precision matrices, and examine scaling in the number
# of observations.

usage <- function() {
  paste(
    "Usage:",
    "  Rscript --vanilla run_omega_likelihood_timing.R [options]",
    "",
    "Options:",
    "  --output-dir PATH       Output directory",
    "  --paired-blocks INTEGER Number of randomized paired timing blocks",
    "  --warmups INTEGER       Untimed warm-up calls per condition",
    "  --bootstrap-reps INTEGER Paired bootstrap repetitions",
    "  --help                  Show this message",
    sep = "\n"
  )
}

parse_integer <- function(value, name, minimum) {
  parsed <- suppressWarnings(as.integer(value))
  if (length(parsed) != 1L || is.na(parsed) || as.character(parsed) != value ||
      parsed < minimum) {
    stop(name, " must be an integer greater than or equal to ", minimum, ".")
  }
  parsed
}

parse_arguments <- function(args) {
  settings <- list(
    output_dir = file.path(
      "output",
      "revision_simulations",
      "internal",
      "omega_precision_deep_sanity_check",
      "timing"
    ),
    paired_blocks = 20L,
    warmups = 2L,
    bootstrap_reps = 5000L
  )

  index <- 1L
  while (index <= length(args)) {
    argument <- args[[index]]
    if (argument == "--help") {
      cat(usage(), "\n")
      quit(save = "no", status = 0L)
    }

    if (startsWith(argument, "--output-dir=")) {
      settings$output_dir <- substring(
        argument,
        nchar("--output-dir=") + 1L
      )
    } else if (startsWith(argument, "--paired-blocks=")) {
      settings$paired_blocks <- parse_integer(
        substring(argument, nchar("--paired-blocks=") + 1L),
        "--paired-blocks",
        2L
      )
    } else if (startsWith(argument, "--warmups=")) {
      settings$warmups <- parse_integer(
        substring(argument, nchar("--warmups=") + 1L),
        "--warmups",
        0L
      )
    } else if (startsWith(argument, "--bootstrap-reps=")) {
      settings$bootstrap_reps <- parse_integer(
        substring(argument, nchar("--bootstrap-reps=") + 1L),
        "--bootstrap-reps",
        10L
      )
    } else if (argument %in% c(
      "--output-dir",
      "--paired-blocks",
      "--warmups",
      "--bootstrap-reps"
    )) {
      if (index == length(args)) {
        stop("Missing value after ", argument, ".")
      }
      value <- args[[index + 1L]]
      if (argument == "--output-dir") {
        settings$output_dir <- value
      } else if (argument == "--paired-blocks") {
        settings$paired_blocks <- parse_integer(
          value,
          "--paired-blocks",
          2L
        )
      } else if (argument == "--warmups") {
        settings$warmups <- parse_integer(value, "--warmups", 0L)
      } else {
        settings$bootstrap_reps <- parse_integer(
          value,
          "--bootstrap-reps",
          10L
        )
      }
      index <- index + 1L
    } else {
      stop("Unknown argument: ", argument, "\n\n", usage())
    }
    index <- index + 1L
  }

  if (!is.character(settings$output_dir) ||
      length(settings$output_dir) != 1L ||
      !nzchar(settings$output_dir)) {
    stop("--output-dir must be a non-empty path.")
  }
  settings
}

write_csv <- function(x, path) {
  utils::write.csv(x, file = path, row.names = FALSE)
}

elapsed_seconds <- function() {
  unname(proc.time()[["elapsed"]])
}

ar1_correlation <- function(size, rho) {
  indices <- seq_len(size)
  rho^abs(outer(indices, indices, "-"))
}

make_banded_ar1_precision <- function(size, rho, standard_error) {
  if (size < 2L || abs(rho) >= 1 || standard_error <= 0) {
    stop("Invalid inputs for an AR(1) precision matrix.")
  }
  scale <- 1 / (standard_error^2 * (1 - rho^2))
  diagonal <- rep(1 + rho^2, size)
  diagonal[c(1L, size)] <- 1
  Matrix::sparseMatrix(
    i = c(seq_len(size), 2:size, seq_len(size - 1L)),
    j = c(seq_len(size), seq_len(size - 1L), 2:size),
    x = scale * c(diagonal, rep(-rho, 2L * (size - 1L))),
    dims = c(size, size)
  )
}

make_benchmark_scenario <- function(size,
                                    n_units = 1L,
                                    rho = 0.75,
                                    standard_error = 0.65,
                                    seed = 20260807L) {
  set.seed(seed + 1009L * size + 9176L * n_units)
  x <- seq(0, 15, length.out = size)
  covariance <- standard_error^2 * ar1_correlation(size, rho)
  error_factor <- chol(covariance)
  innovations <- matrix(
    stats::rnorm(n_units * size),
    nrow = n_units,
    ncol = size
  )
  errors <- innovations %*% error_factor
  latent_mean <- vapply(
    seq_len(n_units),
    function(unit) {
      0.30 * sin(x / 3 + unit / 5) +
        0.08 * cos(x / 2 - unit / 7) +
        0.015 * (unit - 1L) * x
    },
    numeric(size)
  )
  latent_mean <- t(latent_mean)
  response <- latent_mean + errors
  standard_errors <- rep(standard_error, size)
  dense_precision <- solve(covariance)
  diagonal_precision <- diag(1 / standard_errors^2, size)
  # Keep the input as a base matrix because fash_set_tmbdat() computes its
  # determinant before converting it to a sparse TMB matrix. Exact off-band
  # zeros are retained by that conversion and are recorded in the metadata.
  banded_precision <- as.matrix(make_banded_ar1_precision(
    size = size,
    rho = rho,
    standard_error = standard_error
  ))

  list(
    n = size,
    n_units = n_units,
    x = x,
    response = response,
    covariance = covariance,
    diagonal_covariance = diag(standard_errors^2, size),
    standard_errors = standard_errors,
    dense_precision = dense_precision,
    diagonal_precision = diagonal_precision,
    banded_precision = banded_precision,
    rho = rho,
    standard_error = standard_error
  )
}

component_inputs <- function(scenario, method) {
  if (method == "standard_error") {
    return(list(Si = scenario$standard_errors, Omegai = NULL))
  }
  if (method == "diagonal_omega") {
    return(list(Si = NULL, Omegai = scenario$diagonal_precision))
  }
  if (method == "dense_omega") {
    return(list(Si = NULL, Omegai = scenario$dense_precision))
  }
  if (method == "banded_omega") {
    return(list(Si = NULL, Omegai = scenario$banded_precision))
  }
  stop("Unknown likelihood method: ", method)
}

component_data <- function(scenario, unit = 1L) {
  data.frame(
    y = scenario$response[unit, ],
    x = scenario$x,
    offset = 0
  )
}

evaluate_component <- function(scenario,
                               method,
                               psd,
                               num_basis,
                               betaprec,
                               order,
                               pred_step) {
  inputs <- component_inputs(scenario, method)
  fashr:::compute_L_gaussian_helper(
    data_i = component_data(scenario),
    Si = inputs$Si,
    Omegai = inputs$Omegai,
    psd_iwp = psd,
    num_basis = num_basis,
    betaprec = betaprec,
    order = order,
    pred_step = pred_step
  )
}

make_tmb_object <- function(tmb_data,
                            use_standard_errors,
                            psd,
                            order,
                            pred_step) {
  if (psd != 0) {
    tmb_data$sigmaIWP <- psd / sqrt(
      pred_step^(2 * order - 1) /
        ((2 * order - 1) * factorial(order - 1)^2)
    )
    parameters <- list(
      W = rep(0, ncol(tmb_data$X) + ncol(tmb_data$B))
    )
    dll <- if (use_standard_errors) "Gaussian_ind" else "Gaussian_dep"
  } else {
    parameters <- list(W = rep(0, ncol(tmb_data$X)))
    dll <- if (use_standard_errors) {
      "Gaussian_ind_fixed"
    } else {
      "Gaussian_dep_fixed"
    }
  }

  TMB::MakeADFun(
    data = tmb_data,
    parameters = parameters,
    random = "W",
    DLL = dll,
    silent = TRUE
  )
}

time_staged_component_batch <- function(scenario,
                                        method,
                                        psd,
                                        batch_size,
                                        num_basis,
                                        betaprec,
                                        order,
                                        pred_step) {
  inputs <- component_inputs(scenario, method)
  data_i <- component_data(scenario)
  stage_seconds <- c(
    tmb_data = 0,
    make_adfun = 0,
    first_marginal_evaluation = 0
  )
  checksum <- 0

  for (iteration in seq_len(batch_size)) {
    start <- elapsed_seconds()
    tmb_data <- fashr::fash_set_tmbdat(
      data_i = data_i,
      Si = inputs$Si,
      Omegai = inputs$Omegai,
      num_basis = num_basis,
      betaprec = betaprec,
      order = order
    )
    after_data <- elapsed_seconds()
    objective <- make_tmb_object(
      tmb_data = tmb_data,
      use_standard_errors = !is.null(inputs$Si),
      psd = psd,
      order = order,
      pred_step = pred_step
    )
    after_object <- elapsed_seconds()
    value <- -objective$fn()
    after_evaluation <- elapsed_seconds()

    stage_seconds[["tmb_data"]] <- stage_seconds[["tmb_data"]] +
      after_data - start
    stage_seconds[["make_adfun"]] <- stage_seconds[["make_adfun"]] +
      after_object - after_data
    stage_seconds[["first_marginal_evaluation"]] <-
      stage_seconds[["first_marginal_evaluation"]] +
      after_evaluation - after_object
    checksum <- checksum + value
  }

  list(
    stage_seconds = stage_seconds,
    total_seconds = sum(stage_seconds),
    last_value = value,
    checksum = checksum
  )
}

time_component_batch <- function(scenario,
                                 method,
                                 psd,
                                 batch_size,
                                 num_basis,
                                 betaprec,
                                 order,
                                 pred_step) {
  checksum <- 0
  start <- elapsed_seconds()
  for (iteration in seq_len(batch_size)) {
    checksum <- checksum + evaluate_component(
      scenario = scenario,
      method = method,
      psd = psd,
      num_basis = num_basis,
      betaprec = betaprec,
      order = order,
      pred_step = pred_step
    )
  }
  elapsed <- elapsed_seconds() - start
  list(elapsed_seconds = elapsed, checksum = checksum)
}

summarize_timing <- function(data, group_columns) {
  if (nrow(data) == 0L || !"ms_per_evaluation" %in% names(data)) {
    stop("Timing data are empty or malformed.")
  }
  group_key <- do.call(
    paste,
    c(lapply(data[group_columns], as.character), sep = "\r")
  )
  grouped_indices <- split(seq_len(nrow(data)), group_key, drop = TRUE)
  rows <- lapply(grouped_indices, function(indices) {
    values <- data$ms_per_evaluation[indices]
    group_values <- data[indices[[1]], group_columns, drop = FALSE]
    cbind(
      group_values,
      data.frame(
        n_measurements = length(values),
        mean_ms = mean(values),
        sd_ms = if (length(values) > 1L) stats::sd(values) else NA_real_,
        median_ms = stats::median(values),
        q25_ms = unname(stats::quantile(values, 0.25)),
        q75_ms = unname(stats::quantile(values, 0.75)),
        min_ms = min(values),
        max_ms = max(values),
        stringsAsFactors = FALSE
      )
    )
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

paired_comparisons <- function(data,
                               strata_columns,
                               baseline_method,
                               candidate_methods,
                               bootstrap_reps,
                               bootstrap_seed) {
  if (length(strata_columns) == 0L) {
    stratum_key <- rep("all", nrow(data))
  } else {
    stratum_key <- do.call(
      paste,
      c(lapply(data[strata_columns], as.character), sep = "\r")
    )
  }
  strata <- split(seq_len(nrow(data)), stratum_key, drop = TRUE)
  rows <- list()
  row_index <- 1L

  for (stratum_indices in strata) {
    stratum_data <- data[stratum_indices, , drop = FALSE]
    stratum_values <- if (length(strata_columns) == 0L) {
      NULL
    } else {
      stratum_data[1, strata_columns, drop = FALSE]
    }

    for (candidate in candidate_methods) {
      baseline <- stratum_data[
        stratum_data$method == baseline_method,
        c("block", "ms_per_evaluation"),
        drop = FALSE
      ]
      names(baseline)[2] <- "baseline_ms"
      candidate_data <- stratum_data[
        stratum_data$method == candidate,
        c("block", "ms_per_evaluation"),
        drop = FALSE
      ]
      names(candidate_data)[2] <- "candidate_ms"
      paired <- merge(baseline, candidate_data, by = "block", sort = TRUE)
      if (nrow(paired) < 2L || any(!is.finite(as.matrix(paired[, -1])))) {
        stop("Incomplete paired timing data for ", candidate, ".")
      }

      differences <- paired$candidate_ms - paired$baseline_ms
      ratios <- paired$candidate_ms / paired$baseline_ms
      set.seed(bootstrap_seed + row_index)
      bootstrap_difference <- numeric(bootstrap_reps)
      bootstrap_ratio <- numeric(bootstrap_reps)
      for (replicate_index in seq_len(bootstrap_reps)) {
        sampled <- sample.int(nrow(paired), nrow(paired), replace = TRUE)
        bootstrap_difference[[replicate_index]] <- mean(differences[sampled])
        bootstrap_ratio[[replicate_index]] <-
          mean(paired$candidate_ms[sampled]) /
          mean(paired$baseline_ms[sampled])
      }

      comparison <- data.frame(
        baseline_method = baseline_method,
        candidate_method = candidate,
        n_paired_blocks = nrow(paired),
        baseline_mean_ms = mean(paired$baseline_ms),
        candidate_mean_ms = mean(paired$candidate_ms),
        mean_difference_ms = mean(differences),
        difference_ci_lower_ms = unname(stats::quantile(
          bootstrap_difference,
          0.025
        )),
        difference_ci_upper_ms = unname(stats::quantile(
          bootstrap_difference,
          0.975
        )),
        ratio_of_means = mean(paired$candidate_ms) /
          mean(paired$baseline_ms),
        ratio_ci_lower = unname(stats::quantile(bootstrap_ratio, 0.025)),
        ratio_ci_upper = unname(stats::quantile(bootstrap_ratio, 0.975)),
        paired_geometric_mean_ratio = exp(mean(log(ratios))),
        relative_overhead_percent = 100 * (
          mean(paired$candidate_ms) / mean(paired$baseline_ms) - 1
        ),
        stringsAsFactors = FALSE
      )
      if (!is.null(stratum_values)) {
        comparison <- cbind(stratum_values, comparison)
      }
      rows[[row_index]] <- comparison
      row_index <- row_index + 1L
    }
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

omega_metadata <- function(scenario, num_basis, betaprec, order) {
  data_i <- component_data(scenario)
  omega_objects <- list(
    diagonal_omega = scenario$diagonal_precision,
    dense_omega = scenario$dense_precision,
    banded_omega = scenario$banded_precision
  )
  covariance_references <- list(
    diagonal_omega = scenario$diagonal_covariance,
    dense_omega = scenario$covariance,
    banded_omega = scenario$covariance
  )

  rows <- lapply(names(omega_objects), function(method) {
    omega <- omega_objects[[method]]
    omega_dense <- as.matrix(omega)
    tmb_data <- fashr::fash_set_tmbdat(
      data_i = data_i,
      Si = NULL,
      Omegai = omega,
      num_basis = num_basis,
      betaprec = betaprec,
      order = order
    )
    tmb_omega <- tmb_data$Omega
    size <- nrow(omega_dense)
    off_band <- abs(row(omega_dense) - col(omega_dense)) > 1L
    eigenvalues <- eigen(
      (omega_dense + t(omega_dense)) / 2,
      symmetric = TRUE,
      only.values = TRUE
    )$values
    reference_covariance <- covariance_references[[method]]

    data.frame(
      n_observations = size,
      method = method,
      input_class = paste(class(omega), collapse = "/"),
      input_bytes = as.numeric(utils::object.size(omega)),
      input_nonzeros = sum(omega_dense != 0),
      input_density = sum(omega_dense != 0) / size^2,
      input_off_band_nonzeros = sum(omega_dense[off_band] != 0),
      input_max_abs_off_band = if (any(off_band)) {
        max(abs(omega_dense[off_band]))
      } else {
        0
      },
      tmb_class = paste(class(tmb_omega), collapse = "/"),
      tmb_bytes = as.numeric(utils::object.size(tmb_omega)),
      tmb_nonzeros = Matrix::nnzero(tmb_omega),
      tmb_density = Matrix::nnzero(tmb_omega) / size^2,
      symmetry_error = max(abs(omega_dense - t(omega_dense))),
      minimum_eigenvalue = min(eigenvalues),
      condition_number = kappa(omega_dense),
      log_determinant = as.numeric(
        determinant(omega_dense, logarithm = TRUE)$modulus
      ),
      inverse_identity_error = max(abs(
        omega_dense %*% reference_covariance - diag(size)
      )),
      max_abs_difference_from_dense_omega = max(abs(
        omega_dense - scenario$dense_precision
      )),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

draw_timing_plot <- function(primary_summary,
                             component_summary,
                             scaling_summary,
                             output_path) {
  method_order <- c(
    "standard_error",
    "diagonal_omega",
    "dense_omega",
    "banded_omega"
  )
  method_labels <- c(
    standard_error = "Diagonal S",
    diagonal_omega = "Diagonal Omega",
    dense_omega = "Dense Omega",
    banded_omega = "Banded Omega"
  )
  method_colors <- c(
    standard_error = "#4C78A8",
    diagonal_omega = "#72B7B2",
    dense_omega = "#E45756",
    banded_omega = "#F2CF5B"
  )

  grDevices::pdf(
    output_path,
    width = 12,
    height = 4.2,
    family = "sans",
    useDingbats = FALSE
  )
  on.exit(grDevices::dev.off(), add = TRUE)
  old_parameters <- graphics::par(
    mfrow = c(1, 3),
    mar = c(5.5, 4.2, 2.5, 0.8),
    oma = c(0, 0, 1.2, 0),
    las = 1
  )
  on.exit(graphics::par(old_parameters), add = TRUE)

  primary_methods <- c("standard_error", "dense_omega")
  primary_values <- primary_summary$mean_ms[
    match(primary_methods, primary_summary$method)
  ]
  graphics::barplot(
    primary_values,
    names.arg = unname(method_labels[primary_methods]),
    col = unname(method_colors[primary_methods]),
    border = NA,
    ylim = c(0, max(primary_values) * 1.08),
    ylab = "Mean ms per component",
    main = "52-grid full path"
  )
  graphics::box()

  component_total <- component_summary[
    component_summary$stage == "total",
    ,
    drop = FALSE
  ]
  psd_order <- c("zero", "positive")
  component_values <- matrix(
    NA_real_,
    nrow = length(method_order),
    ncol = length(psd_order),
    dimnames = list(method_order, psd_order)
  )
  for (method in method_order) {
    for (psd_case in psd_order) {
      selected <- component_total$method == method &
        component_total$psd_case == psd_case
      component_values[method, psd_case] <- component_total$mean_ms[selected]
    }
  }
  graphics::barplot(
    component_values,
    beside = TRUE,
    names.arg = c("PSD = 0", "PSD > 0"),
    col = unname(method_colors[method_order]),
    border = NA,
    ylim = c(0, max(component_values) * 1.08),
    ylab = "Mean ms per component",
    main = "Zero/positive breakdown"
  )
  graphics::legend(
    "topleft",
    legend = unname(method_labels[method_order]),
    fill = unname(method_colors[method_order]),
    border = NA,
    bty = "n",
    cex = 0.72
  )
  graphics::box()

  scaling_range <- range(scaling_summary$mean_ms, finite = TRUE)
  scaling_n <- sort(unique(scaling_summary$n_observations))
  graphics::plot(
    scaling_n,
    rep(NA_real_, length(scaling_n)),
    type = "n",
    ylim = c(0, scaling_range[[2]] * 1.08),
    xlab = "Number of observations",
    ylab = "Mean ms per positive-PSD component",
    main = "Observation-count scaling"
  )
  for (method in method_order) {
    method_data <- scaling_summary[
      scaling_summary$method == method,
      ,
      drop = FALSE
    ]
    method_data <- method_data[order(method_data$n_observations), ]
    graphics::lines(
      method_data$n_observations,
      method_data$mean_ms,
      type = "b",
      pch = 16,
      lwd = 1.5,
      col = method_colors[[method]]
    )
  }
  graphics::legend(
    "topleft",
    legend = unname(method_labels[method_order]),
    col = unname(method_colors[method_order]),
    lty = 1,
    pch = 16,
    bty = "n",
    cex = 0.72
  )
  graphics::box()
  graphics::mtext(
    "fashr Gaussian likelihood timing",
    outer = TRUE,
    font = 2,
    cex = 1.05
  )
}

settings <- parse_arguments(commandArgs(trailingOnly = TRUE))
dir.create(settings$output_dir, recursive = TRUE, showWarnings = FALSE)

# Keep the comparison single-threaded at the process-library level as well as
# inside fashr. These variables are set before package namespaces are loaded.
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

required_packages <- c("fashr", "Matrix", "TMB")
missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]
if (length(missing_packages) > 0L) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

# fash_set_tmbdat() calls determinant() on Omega. Attaching Matrix is required
# for determinant dispatch when Omega is supplied as a sparse Matrix object.
suppressPackageStartupMessages(library(Matrix))

benchmark_grid <- sort(c(0, exp(-0.5 * seq(0, 10, by = 0.2))))
positive_psd <- 0.5
num_basis <- 20L
betaprec <- 0
order <- 1L
pred_step <- 1
primary_n <- 16L
primary_units <- 4L
component_batch_size <- 25L
scaling_batch_size <- 10L
scaling_n <- c(8L, 16L, 32L, 64L)
method_order <- c(
  "standard_error",
  "diagonal_omega",
  "dense_omega",
  "banded_omega"
)

validation_rows <- list()
add_validation <- function(check,
                           value,
                           tolerance,
                           passed,
                           details = "") {
  validation_rows[[length(validation_rows) + 1L]] <<- data.frame(
    check = check,
    value = as.numeric(value),
    tolerance = as.numeric(tolerance),
    passed = isTRUE(passed),
    details = details,
    stringsAsFactors = FALSE
  )
}

add_validation(
  check = "default_grid_has_52_components",
  value = abs(length(benchmark_grid) - 52L),
  tolerance = 0,
  passed = length(benchmark_grid) == 52L,
  details = "Grid is c(0, exp(-0.5 * seq(0, 10, by = 0.2)))."
)
add_validation(
  check = "default_grid_has_one_zero_component",
  value = abs(sum(benchmark_grid == 0) - 1L),
  tolerance = 0,
  passed = sum(benchmark_grid == 0) == 1L
)

primary_scenario <- make_benchmark_scenario(
  size = primary_n,
  n_units = primary_units
)
scaling_scenarios <- setNames(
  lapply(
    scaling_n,
    function(size) make_benchmark_scenario(size = size, n_units = 1L)
  ),
  as.character(scaling_n)
)

metadata <- do.call(
  rbind,
  lapply(
    scaling_scenarios,
    omega_metadata,
    num_basis = num_basis,
    betaprec = betaprec,
    order = order
  )
)
rownames(metadata) <- NULL
add_validation(
  check = "all_precision_matrices_are_symmetric",
  value = max(metadata$symmetry_error),
  tolerance = 1e-12,
  passed = max(metadata$symmetry_error) <= 1e-12
)
add_validation(
  check = "all_precision_matrices_are_positive_definite",
  value = min(metadata$minimum_eigenvalue),
  tolerance = 0,
  passed = min(metadata$minimum_eigenvalue) > 0,
  details = "Value is the smallest eigenvalue; it must be positive."
)
add_validation(
  check = "precision_inverse_identity_error",
  value = max(metadata$inverse_identity_error),
  tolerance = 1e-10,
  passed = max(metadata$inverse_identity_error) <= 1e-10
)
banded_difference <- metadata$max_abs_difference_from_dense_omega[
  metadata$method == "banded_omega"
]
add_validation(
  check = "analytic_banded_ar1_matches_dense_inverse",
  value = max(banded_difference),
  tolerance = 1e-10,
  passed = max(banded_difference) <= 1e-10
)

primary_data_standard_error <- fashr::fash_set_data(
  Y = primary_scenario$response,
  smooth_var = primary_scenario$x,
  offset = 0,
  S = primary_scenario$standard_errors,
  Omega = NULL
)
primary_data_dense_omega <- fashr::fash_set_data(
  Y = primary_scenario$response,
  smooth_var = primary_scenario$x,
  offset = 0,
  S = NULL,
  Omega = primary_scenario$dense_precision
)

run_primary <- function(method) {
  fash_data <- if (method == "standard_error") {
    primary_data_standard_error
  } else if (method == "dense_omega") {
    primary_data_dense_omega
  } else {
    stop("Primary benchmark only supports standard_error and dense_omega.")
  }
  fashr:::fash_L_compute(
    fash_data = fash_data,
    likelihood = "gaussian",
    num_cores = 1L,
    grid = benchmark_grid,
    pred_step = pred_step,
    num_basis = num_basis,
    betaprec = betaprec,
    order = order,
    verbose = FALSE
  )
}

primary_validation_standard_error <- run_primary("standard_error")
primary_validation_dense_omega <- run_primary("dense_omega")
expected_primary_dimension <- c(primary_units, length(benchmark_grid))
primary_dimension_error <- max(abs(
  dim(primary_validation_standard_error) - expected_primary_dimension
))
add_validation(
  check = "primary_likelihood_matrix_dimensions",
  value = primary_dimension_error,
  tolerance = 0,
  passed = primary_dimension_error == 0 &&
    identical(
      dim(primary_validation_standard_error),
      dim(primary_validation_dense_omega)
    )
)
add_validation(
  check = "primary_likelihood_matrices_are_finite",
  value = sum(!is.finite(primary_validation_standard_error)) +
    sum(!is.finite(primary_validation_dense_omega)),
  tolerance = 0,
  passed = all(is.finite(primary_validation_standard_error)) &&
    all(is.finite(primary_validation_dense_omega))
)

for (psd in c(0, positive_psd)) {
  standard_error_value <- evaluate_component(
    primary_scenario,
    "standard_error",
    psd,
    num_basis,
    betaprec,
    order,
    pred_step
  )
  diagonal_omega_value <- evaluate_component(
    primary_scenario,
    "diagonal_omega",
    psd,
    num_basis,
    betaprec,
    order,
    pred_step
  )
  dense_omega_value <- evaluate_component(
    primary_scenario,
    "dense_omega",
    psd,
    num_basis,
    betaprec,
    order,
    pred_step
  )
  banded_omega_value <- evaluate_component(
    primary_scenario,
    "banded_omega",
    psd,
    num_basis,
    betaprec,
    order,
    pred_step
  )
  diagonal_error <- abs(standard_error_value - diagonal_omega_value)
  correlated_representation_error <- abs(
    dense_omega_value - banded_omega_value
  )
  suffix <- if (psd == 0) "zero" else "positive"
  add_validation(
    check = paste0("diagonal_omega_equals_standard_error_", suffix),
    value = diagonal_error,
    tolerance = 1e-8,
    passed = diagonal_error <= 1e-8
  )
  add_validation(
    check = paste0("banded_omega_equals_dense_omega_", suffix),
    value = correlated_representation_error,
    tolerance = 1e-8,
    passed = correlated_representation_error <= 1e-8
  )

  for (method in method_order) {
    staged <- time_staged_component_batch(
      scenario = primary_scenario,
      method = method,
      psd = psd,
      batch_size = 1L,
      num_basis = num_basis,
      betaprec = betaprec,
      order = order,
      pred_step = pred_step
    )
    direct <- evaluate_component(
      primary_scenario,
      method,
      psd,
      num_basis,
      betaprec,
      order,
      pred_step
    )
    staged_error <- abs(staged$last_value - direct)
    add_validation(
      check = paste0("staged_path_matches_helper_", method, "_", suffix),
      value = staged_error,
      tolerance = 1e-10,
      passed = staged_error <= 1e-10
    )
  }
}

cat("Running primary 52-grid paired benchmark...\n")
for (method in c("standard_error", "dense_omega")) {
  for (warmup in seq_len(settings$warmups)) {
    invisible(run_primary(method))
  }
}

set.seed(202608071L)
primary_rows <- list()
primary_index <- 1L
for (block in seq_len(settings$paired_blocks)) {
  invisible(gc())
  randomized_methods <- sample(c("standard_error", "dense_omega"))
  for (method in randomized_methods) {
    start <- elapsed_seconds()
    result <- run_primary(method)
    elapsed <- elapsed_seconds() - start
    primary_rows[[primary_index]] <- data.frame(
      block = block,
      method = method,
      order_within_block = match(method, randomized_methods),
      elapsed_seconds = elapsed,
      n_units = primary_units,
      n_observations = primary_n,
      num_basis = num_basis,
      grid_size = length(benchmark_grid),
      n_component_evaluations = primary_units * length(benchmark_grid),
      ms_per_evaluation = 1000 * elapsed /
        (primary_units * length(benchmark_grid)),
      checksum = sum(result),
      stringsAsFactors = FALSE
    )
    primary_index <- primary_index + 1L
  }
}
primary_raw <- do.call(rbind, primary_rows)
primary_summary <- summarize_timing(primary_raw, "method")
primary_comparison <- paired_comparisons(
  data = primary_raw,
  strata_columns = character(),
  baseline_method = "standard_error",
  candidate_methods = "dense_omega",
  bootstrap_reps = settings$bootstrap_reps,
  bootstrap_seed = 202608072L
)

cat("Running zero/positive component breakdown...\n")
for (psd in c(0, positive_psd)) {
  for (method in method_order) {
    for (warmup in seq_len(settings$warmups)) {
      invisible(evaluate_component(
        primary_scenario,
        method,
        psd,
        num_basis,
        betaprec,
        order,
        pred_step
      ))
    }
  }
}

set.seed(202608073L)
component_rows <- list()
component_index <- 1L
for (psd in c(0, positive_psd)) {
  psd_case <- if (psd == 0) "zero" else "positive"
  for (block in seq_len(settings$paired_blocks)) {
    invisible(gc())
    randomized_methods <- sample(method_order)
    for (method in randomized_methods) {
      timed <- time_staged_component_batch(
        scenario = primary_scenario,
        method = method,
        psd = psd,
        batch_size = component_batch_size,
        num_basis = num_basis,
        betaprec = betaprec,
        order = order,
        pred_step = pred_step
      )
      stage_values <- c(timed$stage_seconds, total = timed$total_seconds)
      for (stage in names(stage_values)) {
        component_rows[[component_index]] <- data.frame(
          block = block,
          method = method,
          order_within_block = match(method, randomized_methods),
          psd_case = psd_case,
          psd = psd,
          stage = stage,
          batch_size = component_batch_size,
          elapsed_seconds = unname(stage_values[[stage]]),
          ms_per_evaluation = 1000 * unname(stage_values[[stage]]) /
            component_batch_size,
          checksum = timed$checksum,
          stringsAsFactors = FALSE
        )
        component_index <- component_index + 1L
      }
    }
  }
}
component_raw <- do.call(rbind, component_rows)
component_summary <- summarize_timing(
  component_raw,
  c("method", "psd_case", "psd", "stage")
)
component_comparison <- paired_comparisons(
  data = component_raw[component_raw$stage == "total", , drop = FALSE],
  strata_columns = c("psd_case", "psd"),
  baseline_method = "standard_error",
  candidate_methods = c(
    "diagonal_omega",
    "dense_omega",
    "banded_omega"
  ),
  bootstrap_reps = settings$bootstrap_reps,
  bootstrap_seed = 202608074L
)

cat("Running observation-count scaling benchmark...\n")
for (scenario in scaling_scenarios) {
  for (method in method_order) {
    for (warmup in seq_len(settings$warmups)) {
      invisible(evaluate_component(
        scenario,
        method,
        positive_psd,
        num_basis,
        betaprec,
        order,
        pred_step
      ))
    }
  }
}

set.seed(202608075L)
scaling_rows <- list()
scaling_index <- 1L
for (size in scaling_n) {
  scenario <- scaling_scenarios[[as.character(size)]]
  for (block in seq_len(settings$paired_blocks)) {
    invisible(gc())
    randomized_methods <- sample(method_order)
    for (method in randomized_methods) {
      timed <- time_component_batch(
        scenario = scenario,
        method = method,
        psd = positive_psd,
        batch_size = scaling_batch_size,
        num_basis = num_basis,
        betaprec = betaprec,
        order = order,
        pred_step = pred_step
      )
      scaling_rows[[scaling_index]] <- data.frame(
        block = block,
        method = method,
        order_within_block = match(method, randomized_methods),
        n_observations = size,
        num_basis = num_basis,
        psd = positive_psd,
        batch_size = scaling_batch_size,
        elapsed_seconds = timed$elapsed_seconds,
        ms_per_evaluation = 1000 * timed$elapsed_seconds /
          scaling_batch_size,
        checksum = timed$checksum,
        stringsAsFactors = FALSE
      )
      scaling_index <- scaling_index + 1L
    }
  }
}
scaling_raw <- do.call(rbind, scaling_rows)
scaling_summary <- summarize_timing(
  scaling_raw,
  c("n_observations", "method")
)
scaling_comparison <- paired_comparisons(
  data = scaling_raw,
  strata_columns = "n_observations",
  baseline_method = "standard_error",
  candidate_methods = c(
    "diagonal_omega",
    "dense_omega",
    "banded_omega"
  ),
  bootstrap_reps = settings$bootstrap_reps,
  bootstrap_seed = 202608076L
)

expected_primary_rows <- settings$paired_blocks * 2L
expected_component_rows <- settings$paired_blocks * 2L *
  length(method_order) * 4L
expected_scaling_rows <- settings$paired_blocks * length(scaling_n) *
  length(method_order)
add_validation(
  check = "primary_timing_row_count",
  value = abs(nrow(primary_raw) - expected_primary_rows),
  tolerance = 0,
  passed = nrow(primary_raw) == expected_primary_rows
)
add_validation(
  check = "component_timing_row_count",
  value = abs(nrow(component_raw) - expected_component_rows),
  tolerance = 0,
  passed = nrow(component_raw) == expected_component_rows
)
add_validation(
  check = "scaling_timing_row_count",
  value = abs(nrow(scaling_raw) - expected_scaling_rows),
  tolerance = 0,
  passed = nrow(scaling_raw) == expected_scaling_rows
)
add_validation(
  check = "primary_timings_are_positive_and_finite",
  value = sum(!is.finite(primary_raw$ms_per_evaluation) |
    primary_raw$ms_per_evaluation <= 0),
  tolerance = 0,
  passed = all(is.finite(primary_raw$ms_per_evaluation)) &&
    all(primary_raw$ms_per_evaluation > 0)
)
component_total <- component_raw[component_raw$stage == "total", , drop = FALSE]
add_validation(
  check = "component_total_timings_are_positive_and_finite",
  value = sum(!is.finite(component_total$ms_per_evaluation) |
    component_total$ms_per_evaluation <= 0),
  tolerance = 0,
  passed = all(is.finite(component_total$ms_per_evaluation)) &&
    all(component_total$ms_per_evaluation > 0)
)
add_validation(
  check = "component_stage_timings_are_nonnegative_and_finite",
  value = sum(!is.finite(component_raw$ms_per_evaluation) |
    component_raw$ms_per_evaluation < 0),
  tolerance = 0,
  passed = all(is.finite(component_raw$ms_per_evaluation)) &&
    all(component_raw$ms_per_evaluation >= 0)
)
add_validation(
  check = "scaling_timings_are_positive_and_finite",
  value = sum(!is.finite(scaling_raw$ms_per_evaluation) |
    scaling_raw$ms_per_evaluation <= 0),
  tolerance = 0,
  passed = all(is.finite(scaling_raw$ms_per_evaluation)) &&
    all(scaling_raw$ms_per_evaluation > 0)
)
add_validation(
  check = "all_paired_comparisons_are_finite",
  value = sum(!is.finite(c(
    primary_comparison$ratio_of_means,
    component_comparison$ratio_of_means,
    scaling_comparison$ratio_of_means
  ))),
  tolerance = 0,
  passed = all(is.finite(c(
    primary_comparison$ratio_of_means,
    component_comparison$ratio_of_means,
    scaling_comparison$ratio_of_means
  )))
)

package_description <- utils::packageDescription("fashr")
run_settings <- data.frame(
  script_version = "1.0.0",
  run_started_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  fashr_version = as.character(utils::packageVersion("fashr")),
  fashr_path = find.package("fashr"),
  fashr_built = if (is.null(package_description$Built)) {
    NA_character_
  } else {
    package_description$Built
  },
  R_version = R.version.string,
  platform = R.version$platform,
  blas = unname(extSoftVersion()[["BLAS"]]),
  paired_blocks = settings$paired_blocks,
  warmups = settings$warmups,
  bootstrap_reps = settings$bootstrap_reps,
  primary_n = primary_n,
  primary_units = primary_units,
  num_basis = num_basis,
  order = order,
  betaprec = betaprec,
  pred_step = pred_step,
  grid_size = length(benchmark_grid),
  positive_psd = positive_psd,
  component_batch_size = component_batch_size,
  scaling_batch_size = scaling_batch_size,
  scaling_n = paste(scaling_n, collapse = ","),
  omp_threads = Sys.getenv("OMP_NUM_THREADS"),
  openblas_threads = Sys.getenv("OPENBLAS_NUM_THREADS"),
  veclib_threads = Sys.getenv("VECLIB_MAXIMUM_THREADS"),
  stringsAsFactors = FALSE
)

output_paths <- list(
  settings = file.path(settings$output_dir, "benchmark_settings.csv"),
  metadata = file.path(settings$output_dir, "omega_metadata.csv"),
  primary_raw = file.path(settings$output_dir, "primary_full_grid_raw.csv"),
  primary_summary = file.path(
    settings$output_dir,
    "primary_full_grid_summary.csv"
  ),
  primary_comparison = file.path(
    settings$output_dir,
    "primary_full_grid_paired_comparison.csv"
  ),
  component_raw = file.path(
    settings$output_dir,
    "component_breakdown_raw.csv"
  ),
  component_summary = file.path(
    settings$output_dir,
    "component_breakdown_summary.csv"
  ),
  component_comparison = file.path(
    settings$output_dir,
    "component_breakdown_paired_comparison.csv"
  ),
  scaling_raw = file.path(settings$output_dir, "n_scaling_raw.csv"),
  scaling_summary = file.path(settings$output_dir, "n_scaling_summary.csv"),
  scaling_comparison = file.path(
    settings$output_dir,
    "n_scaling_paired_comparison.csv"
  ),
  validation = file.path(settings$output_dir, "validation_checks.csv"),
  plot = file.path(settings$output_dir, "omega_likelihood_timing.pdf"),
  session_info = file.path(settings$output_dir, "session_info.txt")
)

write_csv(run_settings, output_paths$settings)
write_csv(metadata, output_paths$metadata)
write_csv(primary_raw, output_paths$primary_raw)
write_csv(primary_summary, output_paths$primary_summary)
write_csv(primary_comparison, output_paths$primary_comparison)
write_csv(component_raw, output_paths$component_raw)
write_csv(component_summary, output_paths$component_summary)
write_csv(component_comparison, output_paths$component_comparison)
write_csv(scaling_raw, output_paths$scaling_raw)
write_csv(scaling_summary, output_paths$scaling_summary)
write_csv(scaling_comparison, output_paths$scaling_comparison)
draw_timing_plot(
  primary_summary = primary_summary,
  component_summary = component_summary,
  scaling_summary = scaling_summary,
  output_path = output_paths$plot
)
writeLines(capture.output(utils::sessionInfo()), output_paths$session_info)

required_output_paths <- unlist(output_paths[names(output_paths) != "validation"])
missing_output_count <- sum(!file.exists(required_output_paths))
add_validation(
  check = "required_outputs_exist",
  value = missing_output_count,
  tolerance = 0,
  passed = missing_output_count == 0
)
validation_checks <- do.call(rbind, validation_rows)
write_csv(validation_checks, output_paths$validation)

cat("\nPrimary full-grid comparison (milliseconds per component):\n")
print(primary_comparison, row.names = FALSE)
cat("\nOmega representation metadata at n =", primary_n, ":\n")
print(
  metadata[
    metadata$n_observations == primary_n,
    c(
      "method",
      "input_class",
      "input_nonzeros",
      "tmb_nonzeros",
      "tmb_density"
    )
  ],
  row.names = FALSE
)
cat(
  "\nValidation checks passed:",
  sum(validation_checks$passed),
  "/",
  nrow(validation_checks),
  "\n"
)
cat("Outputs written to:", normalizePath(settings$output_dir), "\n")

if (!all(validation_checks$passed)) {
  print(validation_checks[!validation_checks$passed, , drop = FALSE])
  stop("One or more timing benchmark validation checks failed.")
}
