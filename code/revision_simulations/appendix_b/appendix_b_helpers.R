APPENDIX_B_FASHR_VERSION <- "0.1.43"
APPENDIX_B_FASHR_REMOTE_SHA <-
  "bf223df75da6e41ae48607a56b4cd12d7c3b24e7"

assert_scalar_number <- function(x, name, lower = -Inf, upper = Inf) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < lower || x > upper) {
    stop(
      name, " must be one finite number in [", lower, ", ", upper, "].",
      call. = FALSE
    )
  }
  invisible(x)
}

fashr_provenance <- function() {
  if (!requireNamespace("fashr", quietly = TRUE)) {
    stop("The fashr package is not installed.", call. = FALSE)
  }

  description <- utils::packageDescription("fashr")
  list(
    package = "fashr",
    version = as.character(utils::packageVersion("fashr")),
    remote_sha = if (is.null(description$RemoteSha)) {
      NA_character_
    } else {
      as.character(description$RemoteSha)
    },
    library_path = normalizePath(
      find.package("fashr"),
      winslash = "/",
      mustWork = TRUE
    ),
    r_version = R.version.string,
    platform = R.version$platform
  )
}

require_appendix_b_fashr <- function() {
  provenance <- fashr_provenance()

  if (!identical(provenance$version, APPENDIX_B_FASHR_VERSION)) {
    stop(
      "Appendix B requires fashr ", APPENDIX_B_FASHR_VERSION,
      "; found ", provenance$version, ".",
      call. = FALSE
    )
  }
  if (!identical(provenance$remote_sha, APPENDIX_B_FASHR_REMOTE_SHA)) {
    stop(
      "Appendix B requires fashr RemoteSha ",
      APPENDIX_B_FASHR_REMOTE_SHA, "; found ",
      if (is.na(provenance$remote_sha)) "<missing>" else provenance$remote_sha,
      ".",
      call. = FALSE
    )
  }

  provenance
}

appendix_b_class_counts <- function(J, rho_dynamic, rho_nonlinear) {
  assert_scalar_number(J, "J", lower = 1)
  if (J != as.integer(J)) {
    stop("J must be an integer.", call. = FALSE)
  }
  assert_scalar_number(rho_dynamic, "rho_dynamic", lower = 0, upper = 1)
  assert_scalar_number(
    rho_nonlinear,
    "rho_nonlinear",
    lower = 0,
    upper = 1
  )
  if (rho_dynamic <= rho_nonlinear) {
    stop("rho_dynamic must be greater than rho_nonlinear.", call. = FALSE)
  }

  counts <- c(
    nondynamic = round(J * (1 - rho_dynamic)),
    linear = round(J * (rho_dynamic - rho_nonlinear)),
    nonlinear = round(J * rho_nonlinear)
  )
  if (sum(counts) != J) {
    stop(
      "Rounded Appendix B class counts do not sum to J: ",
      paste(counts, collapse = "/"), " versus J = ", J, ".",
      call. = FALSE
    )
  }

  stats::setNames(as.integer(counts), names(counts))
}

build_appendix_b_datasets <- function(
    J,
    rho_dynamic,
    rho_nonlinear,
    sigma_vec = c(0.1, 0.3, 0.5)) {
  if (!is.numeric(sigma_vec) || length(sigma_vec) < 1L ||
      any(!is.finite(sigma_vec)) || any(sigma_vec <= 0)) {
    stop("sigma_vec must contain finite positive values.", call. = FALSE)
  }

  counts <- appendix_b_class_counts(J, rho_dynamic, rho_nonlinear)
  make_group <- function(n, generator) {
    if (n == 0L) {
      return(list())
    }
    lapply(seq_len(n), function(index) generator())
  }

  nondynamic <- make_group(counts[["nondynamic"]], function() {
    fashr::simulate_process(
      sd_poly = 1,
      type = "nondynamic",
      sd = sigma_vec,
      normalize = FALSE
    )
  })
  linear <- make_group(counts[["linear"]], function() {
    fashr::simulate_process(
      sd_poly = 1,
      type = "linear",
      sd = sigma_vec,
      normalize = FALSE
    )
  })
  nonlinear <- make_group(counts[["nonlinear"]], function() {
    fashr::simulate_process(
      sd_poly = 0,
      type = "nonlinear",
      sd = sigma_vec,
      sd_fun = 5,
      p = 2,
      normalize = FALSE
    )
  })

  data_list <- c(nondynamic, linear, nonlinear)
  unit_id <- c(
    if (counts[["nondynamic"]] > 0L) {
      paste0("A", seq_len(counts[["nondynamic"]]))
    },
    if (counts[["linear"]] > 0L) {
      paste0("B", seq_len(counts[["linear"]]))
    },
    if (counts[["nonlinear"]] > 0L) {
      paste0("C", seq_len(counts[["nonlinear"]]))
    }
  )
  names(data_list) <- unit_id

  class <- rep(
    c("nondynamic", "linear", "nonlinear"),
    times = unname(counts)
  )
  truth <- data.frame(
    index = seq_len(J),
    unit_id = unit_id,
    class = class,
    dynamic_true = class != "nondynamic",
    nonlinear_true = class == "nonlinear",
    stringsAsFactors = FALSE
  )

  stopifnot(
    length(data_list) == J,
    nrow(truth) == J,
    identical(names(data_list), truth$unit_id)
  )

  list(
    data_list = data_list,
    truth = truth,
    class_counts = counts,
    parameters = list(
      J = as.integer(J),
      rho_dynamic = rho_dynamic,
      rho_nonlinear = rho_nonlinear,
      sigma_vec = sigma_vec
    )
  )
}

capture_warnings <- function(expression) {
  warnings <- character()
  value <- withCallingHandlers(
    expression,
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, warnings = unique(warnings))
}

fit_appendix_b_models <- function(
    data_bundle,
    grid,
    penalty,
    num_basis = 20L,
    num_cores = 1L,
    verbose = FALSE) {
  require_appendix_b_fashr()
  if (!is.list(data_bundle) || is.null(data_bundle$data_list) ||
      is.null(data_bundle$truth)) {
    stop("data_bundle must come from build_appendix_b_datasets().", call. = FALSE)
  }
  if (!is.numeric(grid) || length(grid) < 2L || any(!is.finite(grid)) ||
      any(grid < 0) || sum(grid == 0) != 1L) {
    stop("grid must contain one zero and at least one positive finite value.", call. = FALSE)
  }
  grid <- sort(unique(grid))
  assert_scalar_number(penalty, "penalty", lower = 1)
  assert_scalar_number(num_basis, "num_basis", lower = 5)
  assert_scalar_number(num_cores, "num_cores", lower = 1)

  fit_one_order <- function(order) {
    captured <- capture_warnings(fashr::fash(
      Y = "y",
      smooth_var = "x",
      S = "sd",
      data_list = data_bundle$data_list,
      num_basis = as.integer(num_basis),
      order = as.integer(order),
      betaprec = 0,
      pred_step = 1,
      penalty = penalty,
      grid = grid,
      num_cores = as.integer(num_cores),
      verbose = verbose
    ))
    raw_fit <- captured$value
    bf_captured <- capture_warnings(
      fashr::BF_update(raw_fit, plot = FALSE)
    )
    list(
      raw = raw_fit,
      bf = bf_captured$value,
      warnings = unique(c(captured$warnings, bf_captured$warnings))
    )
  }

  list(
    order1 = fit_one_order(1L),
    order2 = fit_one_order(2L),
    truth = data_bundle$truth,
    class_counts = data_bundle$class_counts,
    settings = list(
      grid = grid,
      penalty = penalty,
      num_basis = as.integer(num_basis),
      num_cores = as.integer(num_cores)
    )
  )
}

extract_fash_null_weight <- function(fit) {
  if (!is.list(fit) ||
      !is.numeric(fit$psd_grid) ||
      sum(fit$psd_grid == 0) != 1L ||
      !is.data.frame(fit$prior_weights) ||
      !all(c("psd", "prior_weight") %in% names(fit$prior_weights))) {
    stop(
      "The fit must define one zero-PSD grid component and a valid prior table.",
      call. = FALSE
    )
  }

  null_row <- which(fit$prior_weights$psd == 0)
  if (length(null_row) == 0L) {
    # fashr stores only positive mixture weights. An omitted zero-PSD row
    # therefore represents an estimated exact-null weight of zero.
    return(0)
  }
  if (length(null_row) > 1L) {
    stop("The fit contains more than one null prior weight.", call. = FALSE)
  }

  null_weight <- as.numeric(fit$prior_weights$prior_weight[null_row])
  if (length(null_weight) != 1L ||
      !is.finite(null_weight) ||
      null_weight < 0 ||
      null_weight > 1) {
    stop("The fit contains an invalid null prior weight.", call. = FALSE)
  }
  null_weight
}

extract_fash_lfdr <- function(fit) {
  if (!is.null(fit$lfdr)) {
    lfdr <- as.numeric(fit$lfdr)
  } else {
    null_column <- which(fit$prior_weights$psd == 0)
    if (length(null_column) != 1L) {
      stop("The fit does not contain exactly one null component.", call. = FALSE)
    }
    lfdr <- as.numeric(fit$posterior_weights[, null_column])
  }
  if (length(lfdr) == 0L || any(!is.finite(lfdr)) ||
      any(lfdr < 0) || any(lfdr > 1)) {
    stop("The fit contains invalid lfdr values.", call. = FALSE)
  }
  lfdr
}

cumulative_lfdr_calls <- function(lfdr, alpha = 0.05) {
  assert_scalar_number(alpha, "alpha", lower = 0, upper = 1)
  if (!is.numeric(lfdr) || length(lfdr) == 0L || any(!is.finite(lfdr)) ||
      any(lfdr < 0) || any(lfdr > 1)) {
    stop("lfdr must contain finite values in [0, 1].", call. = FALSE)
  }

  ordered_indices <- order(lfdr, seq_along(lfdr))
  cumulative_lfdr <- cumsum(lfdr[ordered_indices]) / seq_along(lfdr)
  accepted <- which(cumulative_lfdr <= alpha)
  n_calls <- if (length(accepted) == 0L) 0L else max(accepted)

  list(
    indices = if (n_calls == 0L) integer() else ordered_indices[seq_len(n_calls)],
    ordered_indices = ordered_indices,
    cumulative_lfdr = cumulative_lfdr,
    n_calls = n_calls
  )
}

evaluate_appendix_b_discoveries <- function(
    lfdr,
    truth_positive,
    alpha = 0.05) {
  if (!is.logical(truth_positive) || length(truth_positive) != length(lfdr) ||
      anyNA(truth_positive)) {
    stop("truth_positive must be a complete logical vector aligned to lfdr.", call. = FALSE)
  }
  call_result <- cumulative_lfdr_calls(lfdr, alpha)
  called <- call_result$indices
  true_discoveries <- sum(truth_positive[called])
  false_discoveries <- length(called) - true_discoveries
  n_positive <- sum(truth_positive)

  data.frame(
    alpha = alpha,
    discoveries = length(called),
    true_discoveries = true_discoveries,
    false_discoveries = false_discoveries,
    realized_fdp = if (length(called) == 0L) 0 else {
      false_discoveries / length(called)
    },
    power = if (n_positive == 0L) NA_real_ else true_discoveries / n_positive,
    stringsAsFactors = FALSE
  )
}

focused_alpha_curve <- function(
    fit_bundle,
    alpha_values = seq(0.01, 0.20, by = 0.01)) {
  if (!is.list(fit_bundle) || is.null(fit_bundle$order1) ||
      is.null(fit_bundle$order2) || is.null(fit_bundle$truth)) {
    stop("fit_bundle must come from fit_appendix_b_models().", call. = FALSE)
  }

  specifications <- list(
    list(order = "IWP1", stage = "Raw EB", fit = fit_bundle$order1$raw,
         truth = fit_bundle$truth$dynamic_true),
    list(order = "IWP1", stage = "BF updated", fit = fit_bundle$order1$bf,
         truth = fit_bundle$truth$dynamic_true),
    list(order = "IWP2", stage = "Raw EB", fit = fit_bundle$order2$raw,
         truth = fit_bundle$truth$nonlinear_true),
    list(order = "IWP2", stage = "BF updated", fit = fit_bundle$order2$bf,
         truth = fit_bundle$truth$nonlinear_true)
  )

  rows <- lapply(specifications, function(specification) {
    lfdr <- extract_fash_lfdr(specification$fit)
    result <- do.call(
      rbind,
      lapply(alpha_values, function(alpha) {
        evaluate_appendix_b_discoveries(
          lfdr = lfdr,
          truth_positive = specification$truth,
          alpha = alpha
        )
      })
    )
    result$order <- specification$order
    result$stage <- specification$stage
    result
  })

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result[, c(
    "order", "stage", "alpha", "discoveries", "true_discoveries",
    "false_discoveries", "realized_fdp", "power"
  )]
}

summarize_appendix_b_grid_fit <- function(
    fit_bundle,
    setting,
    grid_spacing,
    rho_dynamic,
    rho_nonlinear,
    J,
    stream_seed,
    elapsed_seconds) {
  counts <- fit_bundle$class_counts
  data.frame(
    setting = setting,
    grid_spacing = grid_spacing,
    rho_dynamic = rho_dynamic,
    rho_nonlinear = rho_nonlinear,
    true_pi0_iwp1 = 1 - rho_dynamic,
    true_pi0_iwp2 = 1 - rho_nonlinear,
    J = as.integer(J),
    n_nondynamic = counts[["nondynamic"]],
    n_linear = counts[["linear"]],
    n_nonlinear = counts[["nonlinear"]],
    stream_seed = as.integer(stream_seed),
    raw_pi0_iwp1 = extract_fash_null_weight(fit_bundle$order1$raw),
    bf_pi0_iwp1 = extract_fash_null_weight(fit_bundle$order1$bf),
    raw_pi0_iwp2 = extract_fash_null_weight(fit_bundle$order2$raw),
    bf_pi0_iwp2 = extract_fash_null_weight(fit_bundle$order2$bf),
    warning_count = length(unique(c(
      fit_bundle$order1$warnings,
      fit_bundle$order2$warnings
    ))),
    elapsed_seconds = elapsed_seconds,
    stringsAsFactors = FALSE
  )
}
